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
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include <memory>

using namespace Go;
using std::ifstream;
using std::ofstream;
using std::endl;

int main(int argc, char* argv[] )
{
    if (argc != 5)
    {
	std::cout << "Usage: " << argv[0]
		  << " input_ellipse tmin tmax output_subcurve" << endl;
	return -1;
    }

    // Open input surface file
    ifstream is(argv[1]);
    double tmin(atof(argv[2]));
    double tmax(atof(argv[3]));
    ofstream os(argv[4]);
    if (is.bad())
    {
	std::cout << "Bad or no input filename" << std::endl;
	return -1;
    }

    // Read surface from file
    ObjectHeader head;
    Ellipse ellipse; // Typically: centre, dir, normal, r1, r2.
    is >> head;
    ASSERT(head.classType() == Ellipse::classType());
    is >> ellipse;

    if (tmin < ellipse.startparam() || tmax > ellipse.endparam())
    {
	std::cout << "tmin or tmax outside domain of ellipse." << std::endl;
	return -1;
    }

    std::cout << "Writing to file." << std::endl;

    // Extract subcurve, write to file.
    ellipse.setParamBounds(tmin, tmax);

    shared_ptr<SplineCurve> sub_ellipse(ellipse.geometryCurve());
 
    sub_ellipse->writeStandardHeader(os);
    sub_ellipse->write(os);

    return 0;
}
