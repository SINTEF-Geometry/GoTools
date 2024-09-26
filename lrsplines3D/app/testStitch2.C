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
#include "GoTools/lrsplines3D/LRVolStitch.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;

int main (int argc, char *argv[]) {

  if (argc != 7) {
    cout << "usage: ./testMBA <input lrsplines(.g2)> <output lrsplines(.g2)> <nmb u> <nmb v> <nmb w> <cont>" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);

  int nmb_u = atoi(argv[3]);
  int nmb_v = atoi(argv[4]);
  int nmb_w = atoi(argv[5]);
  int cont = atoi(argv[6]);
  double eps = 0.001;
  
  int nmb = nmb_u*nmb_v*nmb_w;

  // Read spline volumes
  vector<shared_ptr<LRSplineVolume> > vols(nmb);
  ObjectHeader header;
  for (int ki=0; ki<nmb; ++ki)
    {
      header.read(ifs);
      vols[ki] = shared_ptr<LRSplineVolume>(new LRSplineVolume());
      vols[ki]->read(ifs);
    }
  
  LRVolStitch stitch;
  stitch.stitchRegVols(vols, nmb_u, nmb_v, nmb_w, eps, cont);
  
  for (size_t ki=0; ki<vols.size(); ++ ki)
    {
      vols[ki]->writeStandardHeader(ofs);
      vols[ki]->write(ofs);
    }
}
 
