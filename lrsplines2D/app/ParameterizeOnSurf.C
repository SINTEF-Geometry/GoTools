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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include <iostream>
#include <fstream>
#include <string.h>

//#define DEBUG

using namespace Go;
using std::vector;


int main(int argc, char *argv[])
{
  if (argc != 5) {
    std::cout << "Usage: surface in (.g2), point cloud (.txt/.xyz), param_points_out.txt, info_out.txt" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream ptsout(argv[3]); 
  std::ofstream infoout(argv[4]); 
  
  // Read surface
  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  if (sf1->dimension() != 3)
    {
      infoout << "3D surface expected. Dimension = " << sf1->dimension() << std::endl;
      return -1; 
    }

  // Read points
  int nmb_pts = 0;
  vector<double> data;
  char xx;
  int del = 3;
  while (!ptsin.eof())
    {
      double tmp;
      ptsin >> tmp;
      data.push_back(tmp);
      for (int ki=1; ki<del; ++ki)
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
  Point mid = 0.5*(low + high);
   bool translate = true;
  if (translate)
    {
      for (int ki=0; ki<nmb_pts; ++ki)
	for (int kj=0; kj<del; ++kj)
	  data[del*ki+kj] -= mid[kj];
      sf1->translate(-mid);
    }
  
#ifdef DEBUG
  // Write translated surface and points to file in g2 format
  std::ofstream of("translated.g2");
  sf1->writeStandardHeader(of);
  sf1->write(of);
  PointCloud3D points(data.begin(), nmb_pts);
  points.writeStandardHeader(of);
  points.write(of);
#endif

  int dim = sf1->dimension();
  int maxiter = 4;
  double aeps = 0.001;

  // double umin = sf1->paramMin(XFIXED);
  // double umax = sf1->paramMax(XFIXED);
  // double vmin = sf1->paramMin(YFIXED);
  // double vmax = sf1->paramMax(YFIXED);

  double *curr;
  double dist;

  double maxdist = 0.0;
  double avdist = 0.0;

  (void)ptsout.precision(15);

  // For each point, project onto surface
  int ki;
  for (ki=0, curr=&data[0]; ki<nmb_pts; ++ki, curr+=3)
    {
      // Get seed
      Point curr_pt(curr, curr+dim);
      LRBSpline2D *bspline = LRSplineUtils::mostComparableBspline(sf1.get(), curr_pt);
      Point seed = bspline->getGrevilleParameter();

      // Perform closest point
      double upar, vpar;
      Point close_pt;
      sf1->closestPoint(curr_pt, upar, vpar, close_pt,
			dist, aeps, maxiter, NULL, seed.begin());

      maxdist = std::max(maxdist, dist);
      avdist += fabs(dist);

      ptsout << upar << " " << vpar << " ";
      if (translate)
	ptsout << curr[0]+mid[0] << " " << curr[1]+mid[1] << " " << curr[2]+mid[2];
      else
	ptsout << curr[0] << " " << curr[1] << " " << curr[2];
      ptsout << std::endl;
    }

  avdist /= (double)nmb_pts;
  infoout << "Number of points: " << nmb_pts << std::endl;
  infoout << "Maximum distance between points and surface: " << maxdist << std::endl;
  infoout << "Average distance between points and surface: " << avdist << std::endl;
}

      




 
