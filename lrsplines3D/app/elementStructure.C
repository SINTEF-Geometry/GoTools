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

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/GoTools.h"

using namespace std;
using namespace Go;

int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};


int main (int argc, char *argv[]) {

  if (argc != 4) {
    cout << "usage: ./elementStructure <input lrspline(.g2)> <output element mid points> <output element boundaries>" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  ofstream ofs2(argv[3]);

  GoTools::init();

  ObjectHeader oh;
  oh.read(ifs);

  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifs);

  int num_el = vol->numElements();
  vector<double> mid_el(3*num_el);
  vector<double> bd_el(72*num_el);
  int ki=0, kj=0;
  Point pos, pos2, dir1, dir2, dir3;
  for (auto it=vol->elementsBegin(); it != vol->elementsEnd(); ++it)
    {
      const Element3D* elem = it->second.get();
      double umin = elem->umin();
      double umax = elem->umax();
      double vmin = elem->vmin();
      double vmax = elem->vmax();
      double wmin = elem->wmin();
      double wmax = elem->wmax();
      mid_el[ki++] = 0.5*(umin+umax);
      mid_el[ki++] = 0.5*(vmin+vmax);
      mid_el[ki++] = 0.5*(wmin+wmax);

      pos = Point(umin, vmin, wmin);
      dir1 = Point(umax-umin, 0.0, 0.0);
      dir2 = Point(0.0, vmax-vmin, 0.0);
      dir3 = Point(0.0, 0.0, wmax-wmin);

      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir1;
      pos2 = pos+dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir2;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos+dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos += dir3;
      pos2 = pos-dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos -= dir1;
      pos2 = pos-dir2;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos2 = pos-dir3;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];

      pos -= dir2;
      pos2 = pos+dir1;
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos[ka];
      for (int ka=0; ka<3; ++ka)
	bd_el[kj++] = pos2[ka];
     }

  // Make point cloud
  PointCloud3D elem_cloud(mid_el.begin(), mid_el.size()/3);

  elem_cloud.writeStandardHeader(ofs);
  elem_cloud.write(ofs);

  LineCloud elem_lines(bd_el.begin(), bd_el.size()/6);
  elem_lines.writeStandardHeader(ofs2);
  elem_lines.write(ofs2);
}
