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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include <iostream>
#include <fstream>
#include <algorithm>

#include <vector>

using std::vector;

int main (int argc, char *argv[]) {

  if (argc != 4) {
    std::cout << "usage: ./scale_intensity <input 4d pt cloud> <output point cloud> <multiplication factor> " << std::endl;
    return -1;
  }

  std::ifstream ifs(argv[1]);
  std::ofstream ofs(argv[2]);
  double fac = atof(argv[3]);

  int num_pts;
  ifs >> num_pts;

  vector<double> pc4d;

  double extent[6];
  extent[0] = extent[2] = extent[4] = 1.0e8; //std::numeric_limits<double>::max();
  extent[1] = extent[3] = extent[5] = -1.0e8; //std::numeric_limits<double>::lowest();
  double minv=1.0e8, maxv=-1.0e8;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifs >> p0 >> p1 >> p2 >> q0;
      q0 *= fac;
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      pc4d.push_back(q0);
    }

  // Write to output
  ofs << num_pts << std::endl;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      for (int ka=0; ka<4; ++ka)
	ofs << pc4d[4*ix+ka] << " ";
      ofs << std::endl;
    }
}
