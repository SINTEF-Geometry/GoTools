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

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <memory> // For std::unique_ptr
#include <iomanip> // For std::setprecision

// GoTools headers
#include "GoTools/lrsplines3D/LRSplineVolume.h" 
#include "GoTools/trivariate/SplineVolume.h"  // For defining the initial B-spline volume
#include "GoTools/lrsplines3D/Mesh3D.h"        // The LR mesh
#include "GoTools/lrsplines3D/Direction3D.h" // For XDIR, YDIR, ZDIR enum 
#include "GoTools/lrsplines3D/DefineRefs3D.h"
#include "GoTools/geometry/ObjectHeader.h"   // For reading the input spline volume
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/PointCloud.h"     // For visualization purposes
#include "GoTools/geometry/LineCloud.h"      // For visualization purposes
#include "GoTools/lrsplines3D/LRSplinePlotUtils3D.h" // For extracting element boundary information

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: lrspline_in (.g2), lrspline_out.g2, num refinements" << std::endl;
    return -1;
  }
  
  std::ifstream invol(argv[1]);
  std::ofstream outvol(argv[2]);
  int numpar = atoi(argv[3]);

  GoTools::init();
  ObjectHeader oh;
  oh.read(invol);
      
  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(invol);

   std::ofstream ofm0("ref_vol_mesh0.g2");
  vector<double> mid0, bd0;
  extractElementMidAndBoundary(vol, mid0, bd0);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud0(mid0.begin(), mid0.size()/3);

  elem_cloud0.writeStandardHeader(ofm0);
  elem_cloud0.write(ofm0);

  // Make line cloud from element boundaries
  LineCloud elem_lines0(bd0.begin(), bd0.size()/6);
  elem_lines0.writeStandardHeader(ofm0);
  elem_lines0.write(ofm0);

  const Array<double,6> span = vol->parameterSpan();
  std::cout << "Span: " << span[0] << " " << span[1]  << " " << span[2]  << " " << span[3]  << " " << span[4]  << " " << span[5]  << std::endl;

  vector<LRSplineVolume::Refinement3D> refs_x, refs_y, refs_z;
  std::cout << "Specify " << numpar << " refinements: u, v, w, dir" << std::endl;
  for (int ki=0; ki<numpar; ++ki)
    {
      double u, v, w;
      int d;
      std::cin >> u >> v >> w >> d;

      int mult = 1;
      if (d == 1 || d == 4)
	{
	  Direction3D dir = XDIR;
	  LRSplineVolume::Refinement3D ref;
	  DefineRefs3D::refineFullSpan(*vol, u, v, w, dir, mult, ref);
	  refs_x.push_back(ref);
	}
      if (d == 2 || d == 4)
	{
	  Direction3D dir = YDIR;
	  LRSplineVolume::Refinement3D ref;
	  DefineRefs3D::refineFullSpan(*vol, u, v, w, dir, mult, ref);
	  refs_y.push_back(ref);
	}
      if (d == 3 || d == 4)
	{
	  Direction3D dir = ZDIR;
	  LRSplineVolume::Refinement3D ref;
	  DefineRefs3D::refineFullSpan(*vol, u, v, w, dir, mult, ref);
	  refs_z.push_back(ref);
	}
    }

  for (size_t kr=0; kr<refs_x.size(); ++kr)
    {
      vol->refine(refs_x[kr], true);
    }

  for (size_t kr=0; kr<refs_y.size(); ++kr)
    {
      vol->refine(refs_y[kr], true);
    }

  for (size_t kr=0; kr<refs_z.size(); ++kr)
    {
      vol->refine(refs_z[kr], true);
    }
  
   std::ofstream ofm("ref_vol_mesh.g2");
  vector<double> mid1, bd1;
  extractElementMidAndBoundary(vol, mid1, bd1);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud1(mid1.begin(), mid1.size()/3);

  elem_cloud1.writeStandardHeader(ofm);
  elem_cloud1.write(ofm);

  // Make line cloud from element boundaries
  LineCloud elem_lines1(bd1.begin(), bd1.size()/6);
  elem_lines1.writeStandardHeader(ofm);
  elem_lines1.write(ofm);

  vol->writeStandardHeader(outvol);
  vol->write(outvol);
}
