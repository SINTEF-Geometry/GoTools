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
#include "GoTools/geometry/ObjectHeader.h"
//#include "GoTools/geometry/GoMakeGnuplot.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;

using namespace Go;

int main(int argc, char* argv[])
{
    ALWAYS_ERROR_IF(argc != 4,"Expecting 3 arguments.");


    SplineCurve the_curve, smooth_curve;
    ObjectHeader header;
    SplineCurve *the_curvep=new SplineCurve;

    ifstream infile2(argv[1]);
    infile2 >> header >> (*the_curvep);

    ifstream infile(argv[1]);
    infile >> header >> the_curve;

    smooth_curve = the_curve;
    smooth_curve.removeKnot(atof(argv[2]));

    ofstream outfile(argv[3]);
    outfile << header <<smooth_curve;

    double tpar = 0;
    Point the_tpoint, smooth_tpoint;
    while (tpar != 1001) {
	cout << "Specify tpar to test difference between curves: ";
	cin >> tpar;
	the_curve.point(the_tpoint, tpar);
	smooth_curve.point(smooth_tpoint, tpar);
	cout << " Distance evaluated in tpar: " << 
	    smooth_tpoint.dist(the_tpoint) << endl;
    }

    return 0;

}
