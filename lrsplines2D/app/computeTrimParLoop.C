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

#include "GoTools/lrsplines2D/LRSurfApproxUtils.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 6) {
    std::cout << "Usage: point cloud (.g2), trim_out.g2, max recursion, div u, div v" << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  int max_rec = atoi(argv[3]);
  int nmb_u = atoi(argv[4]);
  int nmb_v = atoi(argv[5]);

  ObjectHeader header;
  header.read(filein);
  PointCloud3D points;
  points.read(filein);

  BoundingBox box = points.boundingBox();
  Point low = box.low();
  Point high = box.high();
  Point mid = 0.5*(low + high);
  Vector3D vec(-mid[0], -mid[1], 0.0);
  points.translate(vec);

  std::ofstream of("translated_cloud.g2");
  points.writeStandardHeader(of);
  points.write(of);


  vector<double> points2(points.rawData(), 
			 points.rawData()+points.dimension()*points.numPoints());
  vector<double> seq;
  LRSurfApproxUtils::computeTrimInfo(points2, 1, max_rec, nmb_u, nmb_v, seq);

  fileout << "410 1 0 0" << std::endl;
  fileout << seq.size()/2-1 << std::endl;
  size_t ki, kj;
  for (ki=0; ki<seq.size()-4; ki+=2)
    {
      for (kj=0; kj<2; ++kj)
	fileout << seq[ki+kj] << "  ";
      fileout << 0 << " ";
      for (; kj<4; ++kj)
	fileout << seq[ki+kj] << "  ";
      fileout << 0 << std::endl;
    }
  fileout << seq[ki] << " " << seq[ki+1] << " " << 0 << " ";
  fileout << seq[0] << " " << seq[1] << " " << 0 << std::endl;
}
