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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"
#include <fstream>

using std::string;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using namespace Go;

//===========================================================================
//                                                                           
// File: rotational_swept_surface.C
//                                                                           
/// Description:
///
/// This program demonstrates the use of the static function
/// 'rotationalSweptSurface' in the class 'SweepSurfaceCreator'.
/// The function can generate a B-spline surface by rotating a curve around 
/// an axis.
/// This program reads a spline curve from file, prints it's bounding box and
/// start and end points. The curve must be such that it doesn't lead to a
/// self-intersecting surface.
/// Then it creates a SplineSurface by calling 'rotationalSweptSurface' and
/// prints it's bounding box, axis and the point on the axis.
/// The spline curve's file name is hardcoded to 'approj_curve.g2', and the the
/// axis direction and a point on the axis are also hardcoded in the program.
///
/// Output is a file in Go-format for plotting the curve and the surface.
/// The file name is hard-coded to "rotational_swept_surface.g2"
///
/// See also example program 'surface_of_revolution.C'
///
//===========================================================================

int main(int argc, char** argv)
{
    const string inp_curve_filename("data/approj_curve.g2");

    cout << "\nRunning program '" << argv[0]
	 << "'\nSpline curve filename= '"
	 << inp_curve_filename.c_str() << "'." << endl;

    // Read spline curve file
    ifstream cfile(inp_curve_filename.c_str());
    if (!cfile) {
	cerr << "\nFile error. Could not open file: "
	     << inp_curve_filename.c_str() << endl;
	return 1;
    }
    shared_ptr<SplineCurve> curve(new SplineCurve);
    ObjectHeader header;
    cfile >> header;
    if (!header.classType() == SplineCurve::classType()) {
	THROW("Object type is NOT SplineCurve.");
    }
    cfile >> (*curve);
    cfile.close();

    // Print some curve information
    Point pnt3d(3);
    curve->point(pnt3d, curve->startparam());     
    cout << "\nSplineCurve:  Dim= " << curve->dimension()
	 << "\nStart.  Param= " << curve->startparam() << "  Point= "
	 << pnt3d << endl;
    curve->point(pnt3d, curve->endparam());    
    cout << "End.  Param= " << curve->endparam() << "  Point= "
	 << pnt3d << endl;
    cout << "Bounding box =   " << curve->boundingBox() << endl;

    // Create a surface by rotating the curve around the axis an angle of 2PI.
    double angle = 2.0*M_PI;
    Point point_on_axis(0.0, 5.0, 200.0); 
    Point axis_dir(1.0, 0.0, 0.0);
    SplineSurface* surf =
	SweepSurfaceCreator::rotationalSweptSurface(*curve, angle,
						    point_on_axis, axis_dir);
    cout << "\nSurface:  Dim= " << surf->dimension() << endl;
    cout << "Bounding box =  "  << surf->boundingBox() << endl;
    cout << "Point on axis =  " << point_on_axis << endl;
    cout << "Axis direction = " << axis_dir << endl;

    // Open output  file
    ofstream fout("rotational_swept_surface.g2");
    // Write curve to file. Colour=red.
    fout << "100 1 0 4 255 0 0  255" << endl;
    curve->write(fout);

    // Write surface to file. Default colour=blue.    
    surf->writeStandardHeader(fout);
    surf->write(fout);

    // Write axis to file. Colour=green.
    double dlength = 1.2*(surf->boundingBox().high()[0] -
			  surf->boundingBox().low()[0]);
    Point endp = point_on_axis + dlength*axis_dir;
    SplineCurve* axis =  new SplineCurve(point_on_axis, endp);
    fout << "100 1 0 4 0 255 0  255" << endl;
    axis->write(fout);

    // cout << "Open the file 'rotational_swept_surface.g2' in 'goview' to look"
    //      << " at the results" << endl;

    delete surf;
    delete axis;

    return 0;
}
