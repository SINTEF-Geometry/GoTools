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

#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/Point.h"
#include <fstream>


using namespace Go;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
    if (argc != 6) {
	MESSAGE("Usage: surface, closed_dir_u, closed_dir_v, "
		"space_curve, cv_par");
	return 0;
    }

    // Read input arguments

    ifstream sffile(argv[1]);
    ALWAYS_ERROR_IF(sffile.bad(), "Input file not found or file corrupt");
    ObjectHeader header;
    header.read(sffile);
    SplineSurface sf;
    sf.read(sffile);

    bool closed_dir_u = (atoi(argv[2]) == 0 ? false : true);
    bool closed_dir_v = (atoi(argv[3]) == 0 ? false : true);

    ifstream cvfile(argv[4]);
    ALWAYS_ERROR_IF(cvfile.bad(), "Input file not found or file corrupt");
    header.read(cvfile);
    SplineCurve space_cv;
    space_cv.read(cvfile);

    double cv_par = atof(argv[5]);

    // Run the test
    shared_ptr<Point> pt 
	= CreatorsUtils::projectCurvePoint(sf, closed_dir_u, closed_dir_v,
					   &space_cv, cv_par);

    cout << "pt:  " << *pt << endl;

    return 0;
}

