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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;



int main(int argc, char *argv[])
{
  if (argc != 7) {
    std::cout << "Usage: surface in (.g2), point cloud(.raw), lrspline_out.g2, tol, maxiter, smoothing factor" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  double aepsge = atof(argv[4]);
  int max_iter = atoi(argv[5]);
  double smoothwg = atof(argv[6]);
  
  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  // Read parameterized points (u, v, x, y, z) or (u, v, z)
  int ki, kj;
  int dim = sf1->dimension();
  int del = dim + 2;
  int nmb_pts;
  ptsin >> nmb_pts;
  vector<double> data(del*nmb_pts);
  for (ki=0; ki<nmb_pts; ++ki)
    for (kj=0; kj<del; ++kj)
      ptsin >> data[del*ki+kj];

  BoundingBox box = sf1->boundingBox();
  Point low = box.low();
  Point high = box.high();
  Point mid = 0.5*(low + high);
  bool translate = true;
  if (translate)
    {
      sf1->translate(-mid);
      for (ki=0; ki<nmb_pts; ++ki)
	for (kj=2; kj<5; ++kj)
	  data[5*ki+kj] -= mid[kj-2];
    }

  // Write translated surface and points to g2 format
  vector<double> data2;
  data2.reserve(nmb_pts*dim);
  for (ki=0, kj=0; ki<nmb_pts; ++ki, kj+=del)
    data2.insert(data2.end(), data.begin()+kj+2, data.begin()+kj+del);
  PointCloud3D cloud(data2.begin(), nmb_pts);
  std::ofstream of1("translated_sf.g2");
  std::ofstream of2("translated_points.g2");
  sf1->writeStandardHeader(of1);
  sf1->write(of1);
  cloud.writeStandardHeader(of2);
  cloud.write(of2);
  
  bool repar = true; //false; //true;
  LRSurfApprox approx(sf1, data, aepsge, true, repar, true);
  approx.setSmoothingWeight(smoothwg);
  approx.setSmoothBoundary(/*false*/ true);
  approx.setFixCorner(true);

  double maxdist, avdist, avdist_total; // will be set below
  int nmb_out_eps;        // will be set below
  shared_ptr<LRSplineSurface> surf = 
    approx.getApproxSurf(maxdist, avdist_total, avdist, nmb_out_eps, max_iter);

  std::cout << "No. elements: " << surf->numElements();
  std::cout << ", maxdist= " << maxdist << "avdist= " << avdist_total;
  std::cout << ", avdist(out)= " << avdist;
  std::cout << ", nmb out= " << nmb_out_eps << std::endl;
  
  if (surf.get())
    {
      if (translate)
	{
	  // Translate back
	  surf->translate(mid);
	}

      surf->writeStandardHeader(fileout);
      surf->write(fileout);

      if (dim == 1)
	{
	  std::ofstream of2("surf_3D.g2");
	  surf->to3D();
	  surf->writeStandardHeader(of2);
	  surf->write(of2);
	}
    }
}

