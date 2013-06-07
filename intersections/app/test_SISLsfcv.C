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

#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <fstream>
#include <iomanip>


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
    if (argc != 4) {
	cout << "Usage: " << argv[0] << " FileSf FileCv aepsge" << endl;
	return 0;
    }

    ObjectHeader header;

    // Read the first curve from file
    ifstream input1(argv[1]);
    if (input1.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input1);
    SplineSurface surf;
    surf.read(input1);
    input1.close();
    
    // Read the second curve from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input2);
    SplineCurve curve;
    curve.read(input2);
    input2.close();

    double aepsge;
    aepsge = atof(argv[3]);

    // Set up and run the intersection algorithm
    SISLCurve* pcurve = Curve2SISL(curve);
    SISLSurf* psurf = GoSurf2SISL(surf);

    double astart1 = curve.startparam();
    double estart2[] = { surf.startparam_u(), surf.startparam_v() };
    double aend1 = curve.endparam();
    double eend2[] = { surf.endparam_u(), surf.endparam_v() };
    cout << "astart1 = " << astart1 << endl
	 << "estart2[] = { " << estart2[0] << ", "
	 << estart2[1] << " }" << endl
	 << "aend1 = " << aend1 << endl
	 << "eend2[] = { " << eend2[0] << ", " 
	 << eend2[1] << " }" << endl;

    double anext1 = astart1;
    double enext2[] = { estart2[0], estart2[1] };
//     double anext1 = 0.0;
//     double enext2[] = { 2.0, 0.0 };
    cout << "anext1 = " << anext1 << endl
	 << "enext2[] = { " << enext2[0] << ", "
	 << enext2[1] << " }" << endl;

    double cpos1;
    double gpos2[2];
    int jstat = 0;
    cout << "jstat = " << jstat << endl;

    s1772(pcurve, psurf, aepsge, astart1, estart2, aend1, eend2,
	  anext1, enext2, &cpos1, gpos2, &jstat);

    // Write the results
    cout << "Results s1772:" << endl
	 << "jstat = " << jstat << endl
	 << "cpos1 = " << cpos1 << endl
	 << "gpos2[] = { " << gpos2[0] << ", "
	 << gpos2[1] << " }" << endl;

    return 0;

}
