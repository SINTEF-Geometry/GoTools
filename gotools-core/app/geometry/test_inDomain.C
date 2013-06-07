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
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include <fstream>

using namespace Go;
using std::vector;

int main(int argc, char** argv)
{
  if (argc != 2 && argc != 4)
    {
      std::cout << "Usage; infile (nmb u pts) (nmb v pts)" << std::endl;
      return -1;
    }

  std::ifstream infile(argv[1]);
  ObjectHeader head;
  head.read(infile);
  BoundedSurface gosf;
  gosf.read(infile);

  const RectDomain& dom = gosf.containingDomain();
  double start_u = dom.umin();
  double end_u = dom.umax();
  double start_v = dom.vmin();
  double end_v = dom.vmax();

  int nmb_u = 20, nmb_v = 20;
  if (argc == 4)
    {
      nmb_u = atoi(argv[2]);
      nmb_v = atoi(argv[3]);
    }
  double del_u = (end_u - start_u)/(double)(nmb_u-1);
  double del_v = (end_v - start_v)/(double)(nmb_v-1);

  vector<Point> inside;
  int ki, kj;
  double u, v;
  for (v=start_v, ki=0; ki<nmb_v; v+=del_v, ki++)
    for (u=start_u, kj=0; kj<nmb_u; u+=del_u, kj++)
      {
	Point curr;

	if (!gosf.inDomain(u,v))
	  continue;

	gosf.point(curr, u, v);
	inside.push_back(curr);
      }

  std::ofstream out_file("inside_pts.g2");
  out_file << "400 1 0 4 255 0 0 255" << std::endl;
  out_file << inside.size() << std::endl;
  for (ki=0; ki < (int)inside.size(); ki++)
    {
      out_file << inside[ki][0] << " " << inside[ki][1] << " ";
      out_file << inside[ki][2] << std::endl;
    }
}


