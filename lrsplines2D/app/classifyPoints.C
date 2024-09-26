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
  if (argc != 7) {
    std::cout << "Usage: surface in (.g2), point cloud (.g2), points_out.g2, grid (0/1), max level, nmb _levels" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  
  int grid = atoi(argv[4]);
  double max_level = atof(argv[5]);
  int nmb_level = atoi(argv[6]);
  double min_level = -max_level;

  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  int dim = sf1->dimension();
  RectDomain rd = sf1->containingDomain();
  int maxiter = 4;
  double aeps = 0.001;

  double umin = sf1->paramMin(XFIXED);
  double umax = sf1->paramMax(XFIXED);
  double vmin = sf1->paramMin(YFIXED);
  double vmax = sf1->paramMax(YFIXED);

  int ki, kj;
  double *curr;
  double dist;
  vector<double> limits(2*nmb_level+1);
  vector<vector<double> > level_points(2*nmb_level+2);

  double limit[2];
  double cell_del[2];
  if (grid)
    {
      std::cout << "Give domain start (umin, umax): " << std::endl;
      for (ki=0; ki<2; ++ki)
	std::cin >> limit[ki];
      std::cout << "Cell size (u, v): " << std::endl;
      for (ki=0; ki<2; ++ki)
	std::cin >> cell_del[ki];
    }

  // Set distance levels 
  double del = max_level/(double)nmb_level;
  limits[nmb_level] = 0;
  for (ki=1; ki<=nmb_level; ++ki)
    {
      limits[nmb_level-ki] = -ki*del;
      limits[nmb_level+ki] = ki*del;
    }

  double maxdist = 0.0;
  double mindist = 0.0;
  double avdist = 0.0;

  // For each point, classify according to distance
  for (ki=0, curr=&data[0]; ki<nmb_pts; ++ki, curr+=3)
    {
      if (dim == 3)
	{
	  // Get seed
	  Point curr_pt(curr, curr+dim);
	  LRBSpline2D *bspline = LRSplineUtils::mostComparableBspline(sf1.get(), curr_pt);
	  Point seed = bspline->getGrevilleParameter();

	  // Perform closest point
	  double upar, vpar;
	  Point close_pt;
	  // double seed[2];
	  // seed[0] = std::max(umin, std::min(umax, curr[0]));
	  // seed[1] = std::max(vmin, std::min(vmax, curr[1]));
	  sf1->closestPoint(curr_pt, upar, vpar, close_pt,
			    dist, aeps, maxiter, NULL, NULL, seed.begin());
	  Point vec = curr_pt - close_pt;
	  Point norm;
	  sf1->normal(norm, upar, vpar);
	  if (vec*norm < 0.0)
	    dist *= -1;
	}
      else
	{
	  // Evaluate
	  Point pos;
	  sf1->point(pos, curr[0], curr[1]);
	  dist = curr[2]-pos[0];
	}

      if (grid && dim == 1)
	{
	  // Compute cell distance
	  // First identify cell
	  int idx1 = (curr[0] - limit[0])/cell_del[0];
	  int idx2 = (curr[1] - limit[1])/cell_del[1];

	  // Evaluate corners
	  Point pos1, pos2, pos3, pos4;
	  double u1 = std::max(umin, limit[0]+idx1*cell_del[0]);
	  double u2 = std::min(umax, limit[0]+(idx1+1)*cell_del[0]);
	  double v1 = std::max(vmin, limit[1]+idx2*cell_del[1]);
	  double v2 = std::min(vmax, limit[1]+(idx2+1)*cell_del[1]);
	  sf1->point(pos1, u1, v1);
	  sf1->point(pos2, u2, v1);
	  sf1->point(pos3, u1, v2);
	  sf1->point(pos4, u2, v2);
	  double dist1 = curr[2]-pos1[0];
	  double dist2 = curr[2]-pos2[0];
	  double dist3 = curr[2]-pos3[0];
	  double dist4 = curr[2]-pos4[0];

	  if (dist1*dist2<0.0 || dist1*dist3<0.0 || dist1*dist4<0.0 || 
	      dist2*dist3<0.0 || dist2*dist4<0.0 || dist3*dist4<0.0)
	    dist = 0.0;
	  else
	    {
	      int sgn = (dist1 >= 0.0) ? 1 : -1;
	      dist = std::min(std::min(fabs(dist1),fabs(dist2)),
			      std::min(fabs(dist3),fabs(dist4)));
	      dist *= sgn;
	    }
	}

      maxdist = std::max(maxdist, dist);
      mindist = std::min(mindist, dist);
      avdist += fabs(dist);

      // Find classification
      for (kj=0; kj<limits.size(); ++kj)
	if (dist < limits[kj])
	  {
	    level_points[kj].push_back(curr[0]);
	    level_points[kj].push_back(curr[1]);
	    level_points[kj].push_back(curr[2]);
	    break;
	  }
      if (kj == limits.size())
	{
	  level_points[kj].push_back(curr[0]);
	  level_points[kj].push_back(curr[1]);
	  level_points[kj].push_back(curr[2]);
	}
    }

  // Write to file
  for (ki=0; ki<level_points.size(); ++ki)
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

      




 
