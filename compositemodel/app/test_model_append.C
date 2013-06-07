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

#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
using namespace std;
using namespace Go;




int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  // Test number of input arguments
  if (argc != 12)
    {
      std::cout << "Input arguments : Input file on IGES format, Input file on IGES format, ";
      std::cout << "x-value of first point on plane, y-value f first point, ... , z-value of third point" << std::endl;
      exit(-1);
    }


  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.1;
  double approx = 0.001;


  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");
  std::ifstream file2(argv[2]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  CompositeModelFactory factory(approx, gap, neighbour, kink, bend);

  shared_ptr<SurfaceModel> sm1, sm2;
  sm1.reset((SurfaceModel*) factory.createFromIges(file1));
  sm2.reset((SurfaceModel*) factory.createFromIges(file2));

  sm1 -> append(sm2);

  Point
    p_x (atof(argv[3]), atof(argv[4]), atof(argv[5])),
    p_y (atof(argv[6]), atof(argv[7]), atof(argv[8])),
    p_z (atof(argv[9]), atof(argv[10]), atof(argv[11]));

  ftPlane pl(p_x, p_y, p_z);

  ftCurve cv = sm1 -> intersect(pl);

  ofstream osCrv ("dumpCurve.g2");
  ofstream osPl ("dumpPlane.g2");

  cv.writeSpaceCurve(osCrv);

  osPl << "200 1 0 0" << std::endl;
  osPl << "3 0" << std:: endl;
  osPl << "2 2" << std::endl;
  osPl << "0 0 1 1" << std::endl;
  osPl << "2 2" << std::endl;
  osPl << "0 0 1 1" << std::endl;
  osPl << p_x << std::endl;
  osPl << p_x << std::endl;
  osPl << p_y << std::endl;
  osPl << p_z << std::endl;

}
