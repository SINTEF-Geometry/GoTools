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

#include <iostream>
#include <fstream>
#include <algorithm>

#include <vector>

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

int compare_w_par(const void* el1, const void* el2)
{
  if (((double*)el1)[2] < ((double*)el2)[2])
    return -1;
  else if (((double*)el1)[2] > ((double*)el2)[2])
    return 1;
  else
    return 0;
}

using std::vector;

int main (int argc, char *argv[]) {

  if (argc != 3) {
    std::cout << "usage: ./subCloud3 <input 4d pt cloud> <output sub cloud>" << std::endl;
    return -1;
  }

  std::ifstream ifs(argv[1]);
  std::ofstream ofs(argv[2]);

  int num_pts;
  ifs >> num_pts;

  vector<double> pc4d;

  double extent[6];
  double u1, u2, v1, v2, w1, w2;
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
  std::cout << "New domain: " << std::endl;
  std::cin >> u1 >> u2 >> v1 >> v2 >> w1 >> w2;

  
  // Sort the points according to the u-parameter
  qsort(&pc4d[0], num_pts, 4*sizeof(double), compare_u_par);

  // Traverse points
  int pp0, pp1;
  for (pp0=0; pp0<(int)pc4d.size() && pc4d[pp0]<u1; pp0+=4);
  for (pp1=pp0; pp1<(int)pc4d.size() && pc4d[pp1]<u2; pp1+=4);

  // Sort the current sub set of points according to the v-parameter
  qsort(&pc4d[0]+pp0, (pp1-pp0)/4, 4*sizeof(double), compare_v_par);

  // Traverse sub set of points
  int pp2, pp3;
  for (pp2=pp0; pp2<pp1 && pc4d[pp2+1]<v1; pp2+=4);
  for (pp3=pp2; pp3<pp1 && pc4d[pp3+1]<v2; pp3+=4);


  // Sort the current sub set of points according to the w-parameter
  qsort(&pc4d[0]+pp2, (pp3-pp2)/4, 4*sizeof(double), compare_w_par);

  // Traverse sub set of points
  int pp4, pp5;
  for (pp4=pp2; pp4<pp3 && pc4d[pp4+2]<w1; pp4+=4);
  for (pp5=pp4; pp5<pp3 && pc4d[pp5+2]<w2; pp5+=4);

  // Write to output
  int nmb_sub = (pp5 - pp4)/4;
  ofs << nmb_sub << std::endl;
  for (int ki=0; ki<nmb_sub; ++ki)
    {
      for (int kj=0; kj<4; ++kj)
	ofs << pc4d[pp4+4*ki+kj] << " ";
      ofs << std::endl;
    }
}
