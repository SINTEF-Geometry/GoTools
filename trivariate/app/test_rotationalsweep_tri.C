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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;


int main(int argc, char* argv[] )
{

    ALWAYS_ERROR_IF(argc != 10, "Usage: " << argv[0]
		    << " surfaceinfile volumeoutfile angle point_x point_y point_z axis_x axis_y axis_z" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no surface input filename");

    // Open output volume file
    ofstream os(argv[2]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    double angle = atof(argv[3]);
    Point pt(atof(argv[4]), atof(argv[5]), atof(argv[6]));
    Point axis(atof(argv[7]), atof(argv[8]), atof(argv[9]));

    // Read surface from file
    SplineSurface surface;
    ObjectHeader head;
    is >> head >> surface;

    SplineVolume* vol = SweepVolumeCreator::rotationalSweptVolume(surface, angle, pt, axis);

    vol->writeStandardHeader(os);
    vol->write(os);




    ofstream oslines("helplines.g2");
    oslines << "100 0 9 0" << endl;
    oslines << "3 0" << endl;
    oslines << "2 2" << endl;
    oslines << "0 0 1 1" << endl;
    oslines << "1 0 0" << endl;
    oslines << "1 0 3" << endl;

    oslines << "100 0 9 0" << endl;
    oslines << "3 0" << endl;
    oslines << "2 2" << endl;
    oslines << "0 0 1 1" << endl;
    oslines << cos(angle) << " " << sin(angle) << " 0" << endl;
    oslines << cos(angle) << " " << sin(angle) << " 6" << endl;

    oslines << "100 0 9 0" << endl;
    oslines << "3 0" << endl;
    oslines << "2 2" << endl;
    oslines << "0 0 1 1" << endl;
    oslines << "0.625 0 0" << endl;
    oslines << "0.625 0 3" << endl;

    oslines << "100 0 9 0" << endl;
    oslines << "3 0" << endl;
    oslines << "2 2" << endl;
    oslines << "0 0 1 1" << endl;
    oslines << cos(angle)/1.6 << " " << sin(angle)/1.6 << " 0" << endl;
    oslines << cos(angle)/1.6 << " " << sin(angle)/1.6 << " 6" << endl;
}
