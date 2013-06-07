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
#include "GoTools/geometry/SplineCurve.h"
#include <memory>
#include <fstream>
#include <iomanip>


using std::ifstream;
using std::cout;
using std::endl;
using std::cerr;
using namespace Go;


int main(int argc, char** argv)
{
    if (argc != 7) {
	cout << "Usage: test_curveSplineCurveInt FILE1 FILE2 "
	     << "pstart1 pstart2 pstop1 pstop2"
	     << endl;
	return 0;
    }
    double pstart1, pstart2, pstop1, pstop2;
    pstart1 = atof(argv[3]);
    pstart2 = atof(argv[4]);
    pstop1 = atof(argv[5]);
    pstop2 = atof(argv[6]);
    ObjectHeader header;

    // Read first curve from file
    ifstream input(argv[1]);
    if (input.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input);
    shared_ptr<SplineCurve> curve1(new SplineCurve());
    curve1->read(input);
    input.close();

   
    
    // Read second curve from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input2);
    shared_ptr<SplineCurve> curve2(new SplineCurve());
    curve2->read(input2);
    input.close();

//     int prec = std::cout.precision(16);
    cout << "\nFile : " << argv[1] << " Parameter range: "
	 << curve1->startparam() <<" " << curve1->endparam();
    cout << "  Start: " << pstart1 << "  Stop: " << pstop1;
    cout << "\nFile : " << argv[2] << " Parameter range: "
	 << curve2->startparam() <<" " << curve2->endparam();
    cout << "  Start: " << pstart2 << "  Stop: " << pstop2 << endl;

    const double eps = 1.e-4;
    cout << "  Eps= " << eps << endl;

//     SplineCurveInt sci1(curve1);

//     SplineCurveInt sci2(curve2);

//     int istat;
//     istat = sci1.checkCoincidence(pstart1, pstop1, eps, &sci2,
// 				  pstart2, pstop2);

//     cout << "\nistat = " << istat;

//     if (istat==0)
// 	cout << "  Curves are not coinciding." << endl;
//     else
// 	cout << "  Curves are coinciding." << endl;

//     std::cout.precision(prec);  

    return 0;
}
