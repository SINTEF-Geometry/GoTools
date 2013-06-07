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

    ALWAYS_ERROR_IF(argc != 7, "Usage: " << argv[0]
		    << " surfaceinfile curveinfile volumeoutfile point_x point_y point_z" << endl);

    // Open input surface file
    ifstream is_surf(argv[1]);
    ALWAYS_ERROR_IF(is_surf.bad(), "Bad or no surface input filename");

    // Open input curve file
    ifstream is_crv(argv[2]);
    ALWAYS_ERROR_IF(is_crv.bad(), "Bad or no curve input filename");

    // Open output volume file
    ofstream os(argv[3]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    Point pt(atof(argv[4]), atof(argv[5]), atof(argv[6]));

    // Read surface from file
    SplineSurface surf;
    ObjectHeader head;
    is_surf >> head >> surf;

    // Read curve from file
    SplineCurve curve;
    is_crv >> head >> curve;

    SplineVolume* vol = SweepVolumeCreator::linearSweptVolume(surf, curve, pt);

    vol->writeStandardHeader(os);
    vol->write(os);
}
