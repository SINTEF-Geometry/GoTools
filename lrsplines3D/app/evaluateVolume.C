/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/DefineRefs3D.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines3D/LRSplinePlotUtils3D.h" // For visualization purposes
#include "GoTools/geometry/PointCloud.h"     
#include "GoTools/geometry/LineCloud.h"
#include <fstream>

using namespace Go;
using std::vector;


int main(int argc, char *argv[])
{
  if (argc != 7)
    {
      std::cout << "Usage: input spline volume, number of derivatives, number of refinements, parameter (u, v, w)" << std::endl;
      exit(-1);
    }

  std::ifstream input(argv[1]);
  int ider = atoi(argv[2]);
  int numref = atoi(argv[3]);
  
  // Read spline volume from file
   // Read header specifying the type of geometry entity
  GoTools::init();
  ObjectHeader header;
  try {
    header.read(input);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(-1);
    }
  
  shared_ptr<SplineVolume> vol1(new SplineVolume());
  vol1->read(input);
  shared_ptr<LRSplineVolume> vol2(new LRSplineVolume(vol1.get(),
							1.0e-6));

  std::ofstream of("lrvol.g2");
  vol2->writeStandardHeader(of);
  vol2->write(of);

  double upar = atof(argv[4]); //0.5*(vol2->startparam_u() + vol2->endparam_u());
  double vpar = atof(argv[5]); //.5*(vol2->startparam_v() + vol2->endparam_v());
  double wpar = atof(argv[6]); //0.5*(vol2->startparam_w() + vol2->endparam_w());

  int numpt = (ider+1)*(ider+2)*(ider+3)/6;
  vector<Point> pts1(numpt), pts2(numpt);
  vol1->point(pts1, upar, vpar, wpar, ider);
  vol2->point(pts2, upar, vpar, wpar, ider);

  std::cout << "Spline volume: " << std::endl;
  for (size_t ki=0; ki<pts1.size(); ++ki)
    std::cout << pts1[ki] << " ";
  std::cout << std::endl;

  std::cout << "LR volume: " << std::endl;
  for (size_t ki=0; ki<pts2.size(); ++ki)
    std::cout << pts2[ki] << " ";
  std::cout << std::endl;

  int kder1 = ider/2;
  int kder2 = ider - kder1/2;
  int kder3 = std::max(0, ider - kder1 - kder2);

  Point p1, p2, p3;
  std::cout << "kder1: " << kder1 << ", kder2: " << kder2 << ", kder3: " << kder3 << std::endl;
  p1 = vol2->operator()(upar, vpar, wpar, kder1, kder2, kder3);
  std::cout << "p1: " << p1 << std::endl;
  p2 = vol2->operator()(upar, vpar, wpar, kder2, kder3, kder1);
  std::cout << "p2: " << p2 << std::endl;
  p3 = vol2->operator()(upar, vpar, wpar, kder3, kder1, kder2);
  std::cout << "p3: " << p3 << std::endl;

  std::ofstream ofel1("vol_elements.g2");
  vector<double> mid0, bd0;
  extractElementMidAndBoundary(vol2, mid0, bd0);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud0(mid0.begin(), mid0.size()/3);

  elem_cloud0.writeStandardHeader(ofel1);
  elem_cloud0.write(ofel1);

  // Make line cloud from element boundaries
  LineCloud elem_lines0(bd0.begin(), bd0.size()/6);
  elem_lines0.writeStandardHeader(ofel1);
  elem_lines0.write(ofel1);

  const Array<double,6> span = vol2->parameterSpan();
  std::cout << "Span: " << span[0] << " " << span[1]  << " " << span[2]  << " " << span[3]  << " " << span[4]  << " " << span[5]  << std::endl;

  vector<LRSplineVolume::Refinement3D> refs_x, refs_y, refs_z;
  std::cout << "Specify " << numref << " refinements: u, v, w, dir" << std::endl;
  for (int ki=0; ki<numref; ++ki)
    {
      double u, v, w;
      int d;
      std::cin >> u >> v >> w >> d;

      int mult = 1;
      if (d == 1 || d == 4)
	{
	  Direction3D dir = XDIR;
	  LRSplineVolume::Refinement3D ref;
	  DefineRefs3D::refineFullSpan(*vol2, u, v, w, dir, mult, ref);
	  refs_x.push_back(ref);
	}
      if (d == 2 || d == 4)
	{
	  Direction3D dir = YDIR;
	  LRSplineVolume::Refinement3D ref;
	  DefineRefs3D::refineFullSpan(*vol2, u, v, w, dir, mult, ref);
	  refs_y.push_back(ref);
	}
      if (d == 3 || d == 4)
	{
	  Direction3D dir = ZDIR;
	  LRSplineVolume::Refinement3D ref;
	  DefineRefs3D::refineFullSpan(*vol2, u, v, w, dir, mult, ref);
	  refs_z.push_back(ref);
	}
    }

  for (size_t kr=0; kr<refs_x.size(); ++kr)
    {
      vol2->refine(refs_x[kr], true);
    }

  for (size_t kr=0; kr<refs_y.size(); ++kr)
    {
      vol2->refine(refs_y[kr], true);
    }

  for (size_t kr=0; kr<refs_z.size(); ++kr)
    {
      vol2->refine(refs_z[kr], true);
    }
  
   std::ofstream ofm("ref_vol_mesh.g2");
  vector<double> mid1, bd1;
  extractElementMidAndBoundary(vol2, mid1, bd1);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud1(mid1.begin(), mid1.size()/3);

  elem_cloud1.writeStandardHeader(ofm);
  elem_cloud1.write(ofm);

  // Make line cloud from element boundaries
  LineCloud elem_lines1(bd1.begin(), bd1.size()/6);
  elem_lines1.writeStandardHeader(ofm);
  elem_lines1.write(ofm);

   std::ofstream of2("ref_vol.g2");
  vol2->writeStandardHeader(of2);
  vol2->write(of2);

  
  vector<Point> pts3(numpt);
  vol2->point(pts3, upar, vpar, wpar, ider);
  
  std::cout << "LR volume refined: " << std::endl;
  for (size_t ki=0; ki<pts3.size(); ++ki)
    std::cout << pts3[ki] << " ";
  std::cout << std::endl;

  Point p4, p5, p6;
  p4 = vol2->operator()(upar, vpar, wpar, kder1, kder2, kder3);
  std::cout << "p4: " << p4 << std::endl;
  p5 = vol2->operator()(upar, vpar, wpar, kder2, kder3, kder1);
  std::cout << "p5: " << p5 << std::endl;
  p6 = vol2->operator()(upar, vpar, wpar, kder3, kder1, kder2);
  std::cout << "p6: " << p6 << std::endl;
}

  
 
