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
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace Go;
using namespace std;

int main(int argc, char** argv)
{

  if (argc != 11)
    {
      cout << "Usage: " << argv[0] << " top_x top_y top_z bot_x bot_y bot_z axis_x axis_y axis_z outfile" << endl;
      exit(-1);
    }

  Point top (atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
  Point bot (atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
  Point axis (atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));

  SplineCurve* line = new SplineCurve(bot, top);
  SplineSurface* srf = SweepSurfaceCreator::rotationalSweptSurface(*line,
								   2 * M_PI,
								   top,
								   axis);

  ofstream outfile(argv[10]);
  if (outfile.bad())
    {
      cout << "Error when creating output file " << argv[1] << endl;
      exit(-1);
    }

  srf->writeStandardHeader(outfile);
  srf->write(outfile);

}
