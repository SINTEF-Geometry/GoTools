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
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"

using namespace std;
using namespace Go;


int main (int argc, char *argv[]) {

  if (argc != 5 && argc != 6) {
    cout << "usage: ./constParamSurf <input lrvol(.g2)> <output lrsurf (.g2)> <par. dir.> <par. val.> (<nmb sfs>)" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  int dir = atoi(argv[3]);
  double val = atof(argv[4]);
  int nmb = 1;
  if (argc == 6)
    nmb = atoi(argv[5]);

  GoTools::init();

  ObjectHeader oh;
  for (int ki=0; ki<nmb; ++ki)
    {
      oh.read(ifs);
      
      shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
      vol->read(ifs);

      Direction3D dir2 = (dir == 0) ? XDIR : ((dir == 1) ? YDIR : ZDIR);
      double t1 = vol->paramMin(dir2);
      double t2 = vol->paramMax(dir2);
      if (val >= t1 && val <= t2)
	{
	  shared_ptr<LRSplineSurface> lrsurf(vol->constParamSurface(val, dir));
      
	  lrsurf->writeStandardHeader(ofs);
	  lrsurf->write(ofs);
	}
    }
}
