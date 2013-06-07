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
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/utils/Point.h"

using namespace std;
using namespace Go;

//===========================================================================
//                                                                           
// File: linear_swept_volume.C
//                                                                           
/// Description:
///
/// This program demonstrates the use of the static function 'linearSweptVolume'
/// in the class 'SweepVolumeCreator'.
/// The function can generate a SplineVolume by sweeping a surface along a curve
/// or sweeping a curve along a surface. 
/// A sweeping point on the curve or the surface must be specified.
/// If the point lies on the the surface, the surface will be swept along the
/// curve. If the point lies on the the curve, the curve will be swept along the
/// surface. The curve and the surface must be such that it doesn't lead to
/// self-intersection.
///
/// This program creates a curve and a surface. The curve is a line segment
/// parallell with the z-axis and the surface is part of an xy-plane.
/// Then it uses this curve and surface to create two volumes.
/// The first volume is created by sweeping the surface along the curve.
/// The second volume is created by sweeping the curve along the surface.
///
/// Output is a file in Go-format. The file name is hard-coded to
/// "linear_swept_volume.g2". The program 'goview' can't display volumes, but you
/// can use the programs 'makeShield' or 'getBoundarySfs' to extract the boundary
/// faces. They write a new file which can be used by 'goview'.
/// Both programs have inputfilename and outputfilename as arguments, but
/// 'makeShield' has a third optional parameter. If this third parameter is set
/// to 0 (zero) and the file has more then one volume, the volumes will be
/// displayed in different colours.
/// 
//===========================================================================

int main(int argc, char** argv)
{
    cout << "\nRunning program " << argv[0]
	 << ".\nNo input from user." << endl;

    Point z_axis(0.0, 0.0, 1.0);
    Point x_axis(1.0, 0.0, 0.0);

    // Define the curve. Line parallell with the z-axis.
    Point loc(3.0, 1.0, 0.0); // Location
    Point dir(z_axis);        // Direction
    Line line(loc, z_axis);
    line.setParamBounds(1.0, 3.0);  // Line segment
    SplineCurve* curve = line.geometryCurve(); // Spline representation
    cout << "\nCurve bounding box = " << curve->boundingBox() << endl;

    // Define the surface. Plane perpendicular to the z-axis.
    loc = Point(0.0, 0.0, 0.0);
    Point normal = z_axis;
    Plane plane(loc, normal);
    plane.setParameterBounds(-4.0, -1.0, -2.0, 1.0); // (u1, v1, u2, v2)
    SplineSurface* surf = plane.geometrySurface(); // Spline representation
    cout << "Surface bounding box = " << surf->boundingBox() << endl;

   // Open output file for curve and surface.
    ofstream fout("vol_lin_sweep_curve_and_surf.g2");
    // Write curve to file. Colour=green.
    fout << "100 1 0 4 0 255 0  255\n" << *curve << endl;
    // Write surface to file. Colour=red.
    fout << "200 1 0 4 255 0 0  255\n" << *surf << endl;
    fout.close();

    // Create a linearly swept volume by sweeping the surface along the curve.
    Point point_on_surf;
    surf->point(point_on_surf, surf->endparam_u(),
	       0.5*(surf->startparam_v() + surf->endparam_v()));
    cout << "\nCreate volume by sweeping the surface along the curve\n";
    cout << "Point on surface= " << point_on_surf << endl;

    SplineVolume* vol1 = SweepVolumeCreator::linearSweptVolume(*surf, *curve, 
							       point_on_surf);
    cout << "Bounding box volume 1  = " << vol1->boundingBox() << endl;
    cout << "Volume 1 is rational?    " << boolalpha << vol1->rational() << endl;
    cout << "Volume 1 is left handed? " << boolalpha << vol1->isLeftHanded() << endl;

    // Create a linearly swept volume by sweeping the curve along the surface.
    Point point_on_curve;
    //curve->point(point_on_curve, 0.5*(curve->startparam() + curve->endparam()));
    curve->point(point_on_curve, curve->startparam());
    cout << "\nCreate volume by sweeping the curve along the surface.\n";
    cout << "Point on curve= " << point_on_curve << endl;

    SplineVolume* vol2 = SweepVolumeCreator::linearSweptVolume(*surf, *curve, 
							       point_on_curve);
    cout << "Bounding box volume 2 =  " << vol2->boundingBox() << endl;
    cout << "Volume 2 is rational?    " << boolalpha << vol2->rational() << endl;
    cout << "Volume 2 is left handed? " << boolalpha << vol2->isLeftHanded() << endl;

    // Open output file for volumes.
    ofstream vout("linear_swept_volume.g2");

    // Write volumes to file.
    vol1->writeStandardHeader(vout);
    vol1->write(vout);

    vol2->writeStandardHeader(vout);
    vol2->write(vout);
    vout.close();

    // cout << "\nRun: makeShield linear_swept_volume.g2 linear_swept_volume_sf.g2"
    //      << " 0\nand open the files 'vol_lin_sweep_curve_and_surf.g2' and "
    //      << "linear_swept_volume_sf.g2\nin this order in 'goview' to look at "
    //      << "the results.\n" << endl;

    delete curve;
    delete surf;
    delete vol1;
    delete vol2;

    return 0;
}
