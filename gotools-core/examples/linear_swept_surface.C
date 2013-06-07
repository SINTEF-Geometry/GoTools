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

#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/utils/Point.h"
#include <fstream>

using namespace std;
using namespace Go;

//===========================================================================
//                                                                           
// File: linear_swept_surface.C
//                                                                           
/// Description:
///
/// This program demonstrates the use of the static function 'linearSweptSurface'
/// in the class 'SweepSurfaceCreator'.
/// The function can generate a B-spline surface by sweeping one curve along
/// another. A given point on the sweeping curve will be swept along the other
/// curve.
/// The curves must be such that it doesn't lead to a self-intersecting surface.
///
/// This program creates two curves, one semi-circle in the xz-plane and a circle
/// segment in the yz-plane.
/// The first surface is created by sweeping the semi-circle's midpoint along
/// the circle segment.
/// The second surface is created by sweeping the circle segment's startpoint
/// along the semi-circle.
///
/// Output is a file in Go-format for plotting the curves and the surfaces.
/// The file name is hard-coded to "linear_swept_surface.g2"

//===========================================================================

int main(int argc, char** argv)
{
    cout << "\nRunning program " << argv[0]
	 << ".\nNo input from user." << endl;

    // Define the first curve. Semi circle in the xz-plane.
    double radius = 2.0;
    Point centre(0.0, 0.0, 0.0);
    Point x_axis(1.0, 0.0, 0.0);  // Direction vector
    Point y_axis(0.0, 1.0, 0.0);  // Direction vector
    Point z_axis(0.0, 0.0, 1.0);  // Direction vector
    Point normal = y_axis;
    Circle circle1(radius, centre, normal, x_axis);  // Full circle
    circle1.setParamBounds(0.0, M_PI);   // Semi circle
    SplineCurve* curve1 = circle1.geometryCurve(); // Spline representation
    cout << "\nBounding box curve 1 = " << curve1->boundingBox() << endl;

    // Define the second curve. Semi circle in the yz-plane.
    normal = x_axis;
    radius = 7.0;
    centre[1] = 8.0;
    Circle circle2(radius, centre, normal, y_axis); // Full circle
    circle2.setParamBounds(0.2*M_PI, 0.9*M_PI);     // Circle segment
    SplineCurve* curve2 = circle2.geometryCurve(); // Spline representation
    cout << "Bounding box curve 2 = " << curve2->boundingBox() << endl;

    // Create a surface by sweeping curve1 along curve2
    Point point_on_curve1; // Sweeping point
    curve1->point(point_on_curve1,
		  0.5*(curve1->startparam()+curve1->endparam()));
    SplineSurface* surf1 = 
	SweepSurfaceCreator::linearSweptSurface(*curve1, *curve2,
						point_on_curve1);
    cout << "\nCreate surface1 by sweeping curve1 along curve2\n";
    cout << "Point on curve1 " << point_on_curve1 << endl;
    cout << "Bounding box surface1 = "   << surf1->boundingBox()   << endl;

    // Create a surface by sweeping curve2 along curve1
    Point point_on_curve2; // Sweeping point
    curve2->point(point_on_curve2, curve2->startparam());
    SplineSurface* surf2 =
	SweepSurfaceCreator::linearSweptSurface(*curve1, *curve2,
						point_on_curve2);
    cout << "\nCreate surface2 by sweeping curve2 along curve1\n";
    cout << "Point on curve2 " << point_on_curve2 << endl;
    cout << "Bounding box surface2 = "   << surf2->boundingBox()   << endl;

    // Open output  file
    ofstream fout("linear_swept_surface.g2");
    // Write curve1 to file. Colour=red.
    fout << "100 1 0 4 255 0 0  255" << endl;
    curve1->write(fout);
    // Write point on curve1 to file. Colour=red.
    fout << "400 1 0 4 255 0 0  255" << endl;
    fout << "1\n";
    fout << point_on_curve1 << endl;

    // Write curve2 to file. Colour=green.
    fout << "100 1 0 4 0 255 0  255" << endl;
    curve2->write(fout);
    // Write point on curve2 to file. Colour=green.
    fout << "400 1 0 4 0 255 0  255" << endl;
    fout << "1\n";
    fout << point_on_curve2 << endl;

    surf1->writeStandardHeader(fout);
    surf1->write(fout);

    surf2->writeStandardHeader(fout);
    surf2->write(fout);
    // cout << "Open the file 'linear_swept_surface.g2' in 'goview' to look"
    //      << " at the results" << endl;

    delete curve1;
    delete curve2;
    delete surf1;
    delete surf2;

    return 0;
}
