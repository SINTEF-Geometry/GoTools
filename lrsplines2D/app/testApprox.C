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
  if (argc != 6) {
    std::cout << "Usage: point cloud (.g2) lrspline_out.g2 tol maxiter to3D(-1/n)" << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  double AEPSGE = atof(argv[3]);
  int max_iter = atoi(argv[4]);
  int to3D = atoi(argv[5]);
  
  ObjectHeader header;
  header.read(filein);
  PointCloud3D points;
  points.read(filein);

  BoundingBox box = points.boundingBox();
  Point low = box.low();
  Vector3D vec(-low[0], -low[1], 0.0);
  points.translate(vec);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  int dim = 1;
  //LRSurfApprox approx(6, 4, 6, 4, data, 1, AEPSGE, true, true);
  LRSurfApprox approx(4, 4, 4, 4, data, 1, AEPSGE, true, true /*false*/);
   approx.setTurn3D(to3D);

  double maxdist, avdist; // will be set below
  int nmb_out_eps;        // will be set below
  shared_ptr<LRSplineSurface> surf = approx.getApproxSurf(maxdist, avdist, nmb_out_eps, max_iter);

  std::cout << "Maxdist= " << maxdist << ", avdist= " << avdist;
  std::cout << ", nmb out= " << nmb_out_eps << std::endl;

  if (surf.get())
    {
      surf->to3D();
      surf->writeStandardHeader(fileout);
      surf->write(fileout);
    }
}

