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
#include "GoTools/lrsplines2D/LRApproxApp.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;



int colours[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};

int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}



int main(int argc, char *argv[])
{
  if (argc != 8) {
    std::cout << "Usage: surface in (.g2), point cloud (.g2), points_out.g2, tolerance, factor1 (positive), factor2 (negative), minimum tolerance" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  
  double tol = atof(argv[4]);
  double fac1 = atof(argv[5]);
  double fac2 = atof(argv[6]);
  double mintol = atof(argv[7]);

  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> surf(new LRSplineSurface());
  surf->read(sfin);

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  int ki, kj, kr, ka;
  vector<vector<double> > level_points(3);

  // Get all knot values in the u-direction
  const double* const uknots = surf->mesh().knotsBegin(XFIXED);
  const double* const uknots_end = surf->mesh().knotsEnd(XFIXED);
  int nmb_knots_u = surf->mesh().numDistinctKnots(XFIXED);
  const double* knotu;

  // Get all knot values in the v-direction
  const double* const vknots = surf->mesh().knotsBegin(YFIXED);
  const double* const vknots_end = surf->mesh().knotsEnd(YFIXED);
  const double* knotv;

  // Construct mesh of element pointers
  vector<Element2D*> elements;
  surf->constructElementMesh(elements);

  double max_above = 0.0, max_below = 0.0, avdist = 0.0;
  int nmb_points = 0;

  // For each point, classify according to distance
  // Sort points in v-direction
  qsort(&data[0], nmb_pts, 3*sizeof(double), compare_v_par);

  double *curr;
  double dist;

  int pp0, pp1;
  Element2D* elem = NULL;
  for (pp0=0, knotv=vknots; pp0<(int)data.size() && data[pp0+1] < (*knotv); 
       pp0+=3);
  for (kj=0, ++knotv; knotv!= vknots_end; ++knotv, ++kj)
    {
      
      for (pp1=pp0; pp1<(int)data.size() && data[pp1+1] < (*knotv); pp1+=3);
      if (knotv+1 == vknots_end)
	for (; pp1<(int)data.size() && data[pp1+1] <= (*knotv); pp1+=3);
      // 	pp1 = (int)data.size();

      // Sort the current sub set of data according to the u-parameter
      qsort(&data[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_u_par);

      // Traverse the relevant points and identify the associated element
      int pp2, pp3;
      for (pp2=pp0, knotu=uknots; pp2<pp1 && data[pp2] < (*knotu); pp2+=3);
      for (ki=0, ++knotu; knotu!=uknots_end; ++knotu, ++ki)
	{
	  for (pp3=pp2; pp3<pp1 && data[pp3] < (*knotu); pp3 += 3);
	  if (knotu+1 == uknots_end)
	    for (; pp3<pp1 && data[pp3] <= (*knotu); pp3+=3);
	  //   pp3 = pp1;
	  
	  // Fetch associated element
	  elem = elements[kj*(nmb_knots_u-1)+ki];

	  int nump = (pp3 - pp2)/3;
	  for (kr=0, curr=&data[pp2]; kr<nump; ++kr, curr+=3)
	    {
	      // Evaluate
	      Point pos;
	      surf->point(pos, curr[0], curr[1], elem);
	      dist = curr[2]-pos[0];

	      max_above = std::max(max_above, dist);
	      max_below = std::min(max_below, dist);
	      avdist += fabs(dist);
	      nmb_points++;
	  
	      // Find classification
	      double tol2 = (curr[2] < 0.0) ?
		tol - fac2*curr[2] : tol + fac1*curr[2];
	      tol2 = std::max(tol2, mintol);
	      int ix = 1;
	      if (fabs(dist) <= tol2)
		ix = 1;
	      else if (dist < 0)
		ix = 0;
	      else 
		ix = 2;
	      level_points[ix].push_back(curr[0]);
	      level_points[ix].push_back(curr[1]);
	      level_points[ix].push_back(curr[2]);
	    }
	  pp2 = pp3;
	}
      pp0 = pp1;
    }
  avdist /= nmb_points;

  // Write to file
  for (ki=0; ki<level_points.size(); ++ki)
    {
      if (level_points[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_points[ki].begin(), level_points[ki].size()/3);

      fileout << "400 1 0 4 " << colours[ki][0] << " " << colours[ki][1];
      fileout << " " << colours[ki][2] << " 255" << std::endl;
      level_cloud.write(fileout);
    }

  std::cout << "Max dist: " << max_above << ", max dist below: " << max_below;
  std::cout << ", average dist: " << avdist << std::endl;
}

      




 
