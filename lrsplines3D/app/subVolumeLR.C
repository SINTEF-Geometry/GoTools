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
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;



int main(int argc, char* argv[] )
{

    ALWAYS_ERROR_IF(argc != 9, "Usage: " << argv[0]
		    << " volumeinfile volumeoutfile from_u from_v from_w to_u to_v to_w" << endl);

    // Open input volume file
    ifstream input(argv[1]);
    ALWAYS_ERROR_IF(input.bad(), "Bad or no input filename");

    // Open output volume swap file 1
    ofstream os(argv[2]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    // Read volume from file
    GoTools::init();
  ObjectHeader header;
  try {
    header.read(input);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(-1);
    }
  
  shared_ptr<LRSplineVolume> vol1(new LRSplineVolume());
  vol1->read(input);

  double fuzzy = 1.0e-9;
  double umin = atof(argv[3]);
  double umax = atof(argv[4]);
  double vmin = atof(argv[5]);
  double vmax = atof(argv[6]);
  double wmin = atof(argv[7]);
  double wmax = atof(argv[8]);
 
  shared_ptr<LRSplineVolume> vol2(vol1->subVolume(umin, vmin,
						  wmin, umax,
						  vmax, wmax,
						  fuzzy));

    vol2->writeStandardHeader(os);
    vol2->write(os);
}
