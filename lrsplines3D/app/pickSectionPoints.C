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

  if (argc != 5 && argc != 6) {
    std::cout << "usage: ./pickSectionPoints <input 4d pt cloud> <output section (g2)> <direction> <value> (<output 4d points>)" << std::endl;
    return -1;
  }

  std::ifstream ifs(argv[1]);
  std::ofstream ofs(argv[2]);
  int dir = atoi(argv[3]);
  double val = atof(argv[4]);
  char* out2;
  bool write4d = false;
  if (argc == 6)
    {
      write4d = true;
      out2 = argv[5];
    }

  if (dir < 0 || dir > 2)
    {
      std::cout << "Direction must be 0, 1 or 2" << std::endl;
      return -1;
    }
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
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      pc4d.push_back(q0);
      extent[0] = std::min(extent[0], p0);
      extent[1] = std::max(extent[1], p0);
      extent[2] = std::min(extent[2], p1);
      extent[3] = std::max(extent[3], p1);
      extent[4] = std::min(extent[4], p2);
      extent[5] = std::max(extent[5], p2);
      minv = std::min(minv, q0);
      maxv = std::max(maxv, q0);
    }

  std::cout << "Range: " << minv << " " << maxv << std::endl;
  std::cout << "Domain: " << extent[0] << " " << extent[1] << " " << extent[2];
  std::cout << " " << extent[3] << " " << extent[4] << " " << extent[5] << std::endl;

  double eps = 1.0e-6;
  if (val < extent[2*dir]-eps || val > extent[2*dir+1]+eps)
    {
      std::cout << "Value out of range" << std::endl;
      return -1;
    }

  // Traverse points and collect points in section
  int pp0;
  vector<double> section_pt;
  for (pp0=0; pp0<(int)pc4d.size(); pp0+=4)
    {
      if (pc4d[pp0+dir] >= val-eps && pc4d[pp0+dir] <= val+eps)
	{
	  for (int ki=0; ki<4; ++ki)
	    {
	      if (ki == dir)
		continue;
	      section_pt.push_back(pc4d[pp0+ki]);
	    }
	}
    }


  // Write to output
  Go::PointCloud3D points(section_pt.begin(), section_pt.size()/3);
  points.writeStandardHeader(ofs);
  points.write(ofs);

  if (write4d)
    {
      std::ofstream of2(out2);
      int nmb_pt = (int)section_pt.size()/3;
      of2 << nmb_pt << std::endl;
      for (int ka=0; ka<nmb_pt; ++ka)
	{
	  int kb, kc;
	  for (kb=0, kc=0; kb<3; ++kb)
	    {
	      if (kb == dir)
		of2 << val;
	      else
		{
		  of2 << section_pt[3*ka+kc];
		  ++kc;
		}
	      of2 << " ";
	    }
	  of2 << section_pt[3*ka+2] << std::endl;
	}
    }
}
