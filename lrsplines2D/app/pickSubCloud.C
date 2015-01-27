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
  if (argc != 4) {
    std::cout << "Usage: point cloud in (.g2), point cloud out(g2), rotate (0/1) " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  int rotate = atoi(argv[3]);

  ObjectHeader header;
  header.read(filein);
  PointCloud3D points;
  points.read(filein);

  BoundingBox box0 = points.boundingBox();
  printf("Domain: [ %13.3f , %13.3f ] x [ %13.3f , %13.3f ] \n",
	 box0.low()[0], box0.high()[0], box0.low()[1] ,box0.high()[1]);

  Vector3D vec1, vec2;
  if (rotate)
    {
      int ki;
      double tmp[3];
      std::cout << "from vec:" << std::endl;
      for (ki=0; ki<3; ++ki)
	std::cin >> tmp[ki];
      vec1 = Vector3D(tmp[0], tmp[1], tmp[2]);
      std::cout << "to vec: " << std::endl;
      for (ki=0; ki<3; ++ki)
	std::cin >> tmp[ki];
      vec2 = Vector3D(tmp[0], tmp[1], tmp[2]);
      vec1.normalize();
      vec2.normalize();
      points.rotate(vec1, vec2);
    }

  std::ofstream of("rotated_points.g2");
  points.writeStandardHeader(of);
  points.write(of);

  BoundingBox box = points.boundingBox();
  printf("Domain: [ %13.3f , %13.3f ] x [ %13.3f , %13.3f ] \n",
	 box.low()[0], box.high()[0], box.low()[1] ,box.high()[1]);
  printf("New domain: \n");
double u1, u2, v1, v2;
  std::cin >> u1;  
  std::cin >> u2;  
  std::cin >> v1;  
  std::cin >> v2;

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  // Sort the points according to the u-parameter
  qsort(&data[0], nmb_pts, 3*sizeof(double), compare_u_par);

  // Traverse points
  int pp0, pp1;
  for (pp0=0; pp0<(int)data.size() && data[pp0]<u1; pp0+=3);
  for (pp1=pp0; pp1<(int)data.size() && data[pp1]<u2; pp1+=3);

  // Sort the current sub set of points according to the v-parameter
  qsort(&data[0]+pp0, (pp1-pp0)/3, 3*sizeof(double), compare_v_par);

  // Traverse sub set of points
  int pp2, pp3;
  for (pp2=pp0; pp2<pp1 && data[pp2+1]<v1; pp2+=3);
  for (pp3=pp2; pp3<pp1 && data[pp3+1]<v2; pp3+=3);

  // Collect output
  PointCloud3D points2(data.begin()+pp2, (pp3-pp2)/3);

  if (rotate)
    {
      // Rotate back
      points2.rotate(vec2, vec1);
    }

  BoundingBox box2 = points2.boundingBox();
  printf("Domain: [ %13.3f , %13.3f ] x [ %13.3f , %13.3f ] \n",
	 box2.low()[0], box2.high()[0], box2.low()[1] ,box2.high()[1]);

  // Write to output
  points2.writeStandardHeader(fileout);
  points2.write(fileout);
}

