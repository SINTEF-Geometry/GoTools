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

#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
  if (argc != 3 && argc != 5)
      cout << "Usage: " << "infile eps, (swap pardir1, swap pardir2)" << endl;

  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  ObjectHeader head;
  is >> head;

  // Read volume from file
  SplineVolume vol;
  is >> vol;

  double eps = atof(argv[2]);

  int dir1=-1, dir2=-1;
  if (argc > 3)
  {
      dir1 = atoi(argv[3]);
      dir2 = atoi(argv[4]);
  }

  if (dir1 >= 0 && dir2 >=0 && dir1 != dir2)
      vol.swapParameterDirection(dir1, dir2);

  // For all parameter directions
  for (int ki=0; ki<3; ki++)
  {
      int period = vol.volumePeriodicity(ki, eps);
      std::cout << "Periodicity: " << period << std::endl;
  }
}


