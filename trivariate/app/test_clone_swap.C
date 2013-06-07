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
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;


int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc != 6, "Usage: " << argv[0]
		    << " volumeinfile1 volumeinfile2 outfileswap1 outfileswap2 outfileclone" << endl);

    // Open input volume file 1
    ifstream is1(argv[1]);
    ALWAYS_ERROR_IF(is1.bad(), "Bad or no input 1 filename");

    // Open input volume file 2
    ifstream is2(argv[2]);
    ALWAYS_ERROR_IF(is2.bad(), "Bad or no input 2 filename");

    // Open output volume swap file 1
    ofstream os_swap1(argv[3]);
    ALWAYS_ERROR_IF(os_swap1.bad(), "Bad output swap 1 filename");

    // Open output volume swap file 2
    ofstream os_swap2(argv[4]);
    ALWAYS_ERROR_IF(os_swap2.bad(), "Bad output swap 2 filename");

    // Open output volume clone file
    ofstream os_clone(argv[5]);
    ALWAYS_ERROR_IF(os_clone.bad(), "Bad output clone filename");

    // Read volume 1 from file
    SplineVolume vol1;
    ObjectHeader head;
    is1 >> head >> vol1;

    // Read volume 2 from file
    SplineVolume vol2;
    is2 >> head >> vol2;

    // Volume for clone
    SplineVolume* vol_clone = vol1.clone();
    vol1.swap(vol2);

    vol1.writeStandardHeader(os_swap1);
    vol1.write(os_swap1);

    vol2.writeStandardHeader(os_swap2);
    vol2.write(os_swap2);

    vol_clone->writeStandardHeader(os_clone);
    vol_clone->write(os_clone);
}
