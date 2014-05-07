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
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include <iostream>
#include <fstream>
#include <string.h>

//#define DEBUG

using namespace Go;
using std::vector;



int main(int argc, char *argv[])
{
  if (argc != 7) {
    std::cout << "Usage: surface in (.g2), point cloud(.txt), surface out.g2, info_out(.txt), tolerance, max number of iterations" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  std::ofstream infoout(argv[4]); 
  double aepsge = atof(argv[5]);
  int max_iter = atoi(argv[6]);
  
  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  if (sf1->dimension() != 1 && sf1->dimension() != 3)
    {
      infoout << "1D or 3D surface expected. Dimension = " << sf1->dimension() << std::endl;
      std::cout << "1D or 3D surface expected. Dimension = " << sf1->dimension() << std::endl;
      return -1; 
    }

  // Read parameterized points (u, v, x, y, z) or (u, v, z)
  int ki, kj;
  int dim = sf1->dimension();
  int del = dim + 2;
  int nmb_pts = 0;
  vector<double> data;
  char xx;
  while (!ptsin.eof())
    {
      double tmp;
      ptsin >> tmp;
      data.push_back(tmp);
      for (ki=1; ki<del; ++ki)
	{
	  ptsin >> xx;
	  if (xx != ',')
	    ptsin.putback(xx);
	  ptsin >> tmp;
	  data.push_back(tmp);
	}
      nmb_pts++;
      Utils::eatwhite(ptsin);
    }

  BoundingBox box = sf1->boundingBox();
  Point low = box.low();
  Point high = box.high();
  Point mid;
  double umin = sf1->paramMin(XFIXED);
  double umax = sf1->paramMax(XFIXED);
  double vmin = sf1->paramMin(YFIXED);
  double vmax = sf1->paramMax(YFIXED);
  if (sf1->dimension() == 1)
    mid.setValue(0.5*(umin+umax), 0.5*(vmin+vmax), 0.0);
  else
    {
      mid = 0.5*(low + high);
      mid[2] = 0.0;
    }

  // Translate surface and point cloud to origo
  if (sf1->dimension() == 1)
    {
      sf1->setParameterDomain(umin - mid[0], umax - mid[0],
			      vmin - mid[1], vmax - mid[1]);
    }
  else
    sf1->translate(-mid);
  for (ki=0; ki<nmb_pts; ++ki)
    for (kj=del-3; kj<del-1; ++kj)
      data[del*ki+kj] -= mid[kj-del+3];

#ifdef DEBUG
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
#endif
  
  std::cout << std::endl;
  std::cout << "Input points and surface read and pre processed. Ready for surface creation.";
  std::cout << std::endl << std::endl;
  
  double smoothwg = 0.00000001;
  bool repar = true; 
  LRSurfApprox approx(sf1, data, aepsge, true, repar, true);
  approx.setSmoothingWeight(smoothwg);
  approx.setSmoothBoundary(true);
  approx.setFixCorner(true);
  approx.setVerbose(true);

  double maxdist, avdist, avdist_total; // will be set below
  int nmb_out_eps;        // will be set below
  shared_ptr<LRSplineSurface> surf = 
    approx.getApproxSurf(maxdist, avdist_total, avdist, nmb_out_eps, max_iter);

  std::cout << std::endl;
  std::cout << "Approximation completed. Writing to output files." << std::endl;
 
  infoout << "Total number of points: " << nmb_pts << std::endl;
  infoout << "Number of elements: " << surf->numElements() << std::endl;
  infoout << "Maximum distance: " << maxdist << std::endl;
  infoout << "Average distance: " << avdist_total << std::endl;
  infoout << "Average distance for points outside of the tolerance: " << avdist << std::endl;
  infoout << "Number of points outside the tolerance: " << nmb_out_eps << std::endl;
  
  if (surf.get())
    {
      // Translate back
      if (surf->dimension() == 3)
	{
	  surf->translate(mid);
	}
      else
	{
	  // Update parameter domain
	  surf->setParameterDomain(umin + mid[0], umax + mid[0],
				   vmin + mid[1], vmax + mid[1]);
	}
      
      surf->writeStandardHeader(fileout);
      surf->write(fileout);

#ifdef DEBUG
      if (dim == 1)
	{
	  std::ofstream of2("surf_3D.g2");
	  surf->to3D();
	  surf->writeStandardHeader(of2);
	  surf->write(of2);
	}
#endif
    }
}

