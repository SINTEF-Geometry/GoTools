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
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;



int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};


int main(int argc, char *argv[])
{
  if (argc < 7) {
    std::cout << "Usage: surface in (.g2), point cloud (.g2), points_out.g2, nmb _levels, (max_level or positive levels)" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  
  int nmb_level = atoi(argv[4]);
  vector<double> limits(2*nmb_level+1);
  if (argc == 7)
    {
      double max_level = atof(argv[5]);
      double ddel = max_level/(double)nmb_level;
      limits[nmb_level] = 0;
      for (int ki=1; ki<=nmb_level; ++ki)
	{
	  limits[nmb_level-ki] = -ki*ddel;
	  limits[nmb_level+ki] = ki*ddel;
	}
    }
  else
    {
      limits[nmb_level] = 0;
      for (int ki=1; ki<=nmb_level; ++ki)
	{
	  double curr_level = atof(argv[4+ki]);
	  limits[ki-1] = -curr_level;
	  limits[nmb_level+ki] = curr_level;
	}
      std::sort(limits.begin(), limits.end());
    }

  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  // Represent the surface as tensor product
  shared_ptr<ParamSurface> tpsf(sf1->asSplineSurface());

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  int dim = tpsf->dimension();
  RectDomain rd = tpsf->containingDomain();
  int maxiter = 4;
  double aeps = 0.001;

  int ki, kj;
  double *curr;
  double dist;
  vector<vector<double> > level_points(2*nmb_level+2);

  double maxdist = 0.0;
  double mindist = 0.0;
  double avdist = 0.0;

  // For each point, classify according to distance
  Point seed;
  double upar, vpar;
  Point close_pt;
  for (ki=0, curr=&data[0]; ki<nmb_pts; ++ki, curr+=3)
    {
      if (dim == 3)
	{
	  // Get seed
	  Point curr_pt(curr, curr+dim);

	  // TEST
	  seed = Point(curr[0], curr[1]);

	  if (seed.dimension() > 0)
	    tpsf->closestPoint(curr_pt, upar, vpar, close_pt,
	  		       dist, aeps, maxiter, NULL, seed.begin());
	  else
	    tpsf->closestPoint(curr_pt, upar, vpar, close_pt,
			       dist, aeps, maxiter);
	  seed.setValue(upar, vpar);

	  // Perform closest point
	  // double seed[2];
	  // seed[0] = std::max(umin, std::min(umax, curr[0]));
	  // seed[1] = std::max(vmin, std::min(vmax, curr[1]));
	  // LRBSpline2D *bspline = LRSplineUtils::mostComparableBspline(sf1.get(), curr_pt);
	  // seed = bspline->getGrevilleParameter();

	  // // Perform closest point
	  // // double seed[2];
	  // // seed[0] = std::max(umin, std::min(umax, curr[0]));
	  // // seed[1] = std::max(vmin, std::min(vmax, curr[1]));
	  // sf1->closestPoint(curr_pt, upar, vpar, close_pt,
	  // 		    dist, aeps, maxiter, NULL, seed.begin());
	  Point vec = curr_pt - close_pt;
	  Point norm;
	  //sf1->normal(norm, upar, vpar);
	  tpsf->normal(norm, upar, vpar);
	  if (vec*norm < 0.0)
	    dist *= -1;
	}
      else
	{
	  // Evaluate
	  Point pos;
	  tpsf->point(pos, curr[0], curr[1]);
	  dist = curr[2]-pos[0];
	}

 
      maxdist = std::max(maxdist, dist);
      mindist = std::min(mindist, dist);
      avdist += fabs(dist);

      // Find classification
      for (kj=0; kj<(int)limits.size(); ++kj)
	if (dist < limits[kj])
	  {
	    level_points[kj].push_back(curr[0]);
	    level_points[kj].push_back(curr[1]);
	    level_points[kj].push_back(curr[2]);
	    break;
	  }
      if (kj == (int)limits.size())
	{
	  level_points[kj].push_back(curr[0]);
	  level_points[kj].push_back(curr[1]);
	  level_points[kj].push_back(curr[2]);
	}
    }

  // Write to file
  for (ki=0; ki<(int)level_points.size(); ++ki)
    {
      if (level_points[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_points[ki].begin(), level_points[ki].size()/3);

      double cc[3];
      if (ki <= nmb_level)
	{
	  cc[0] = ((nmb_level-ki)*colors[0][0] + ki*colors[1][0])/nmb_level;
	  cc[1] = ((nmb_level-ki)*colors[0][1] + ki*colors[1][1])/nmb_level;
	  cc[2] = ((nmb_level-ki)*colors[0][2] + ki*colors[1][2])/nmb_level;
	}
      else
	{
	  cc[0] = ((ki-nmb_level-1)*colors[2][0] + 
		   (2*nmb_level-ki+1)*colors[1][0])/nmb_level;
	  cc[1] = ((ki-nmb_level-1)*colors[2][1] + 
		   (2*nmb_level-ki+1)*colors[1][1])/nmb_level;
	  cc[2] = ((ki-nmb_level-1)*colors[2][2] + 
		   (2*nmb_level-ki+1)*colors[1][2])/nmb_level;
	}

      fileout << "400 1 0 4 " << cc[0] << " " << cc[1];
      fileout << " " << cc[2] << " 255" << std::endl;
      level_cloud.write(fileout);
    }

  avdist /= (double)nmb_pts;

  std::cout << "Max dist: " << maxdist << "Max dist below: " << mindist;
  std::cout << ", average dist: " << avdist << std::endl;
}

      




 
