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
#include "GoTools/creators/ProjectCurve.h"
#include "GoTools/geometry/SurfaceOfRevolution.h"
#include <fstream>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using namespace Go;

//===========================================================================
//                                                                           
// File: surface_of_revolution.C
//                                                                           
/// Description:
///  
/// This program demonstrates the use of the class SurfaceOfRevolution.
///
/// SurfaceOfRevolution is swept out by a SplineCurve that is rotated
/// around an axis with a complete revolution, and is thereby a
/// parametric surface. The space dimension is 3.
/// The curve must be such that it doesn't lead to a self-intersecting surface.
///
/// This program reads a spline curve from file, prints it's bounding box and
/// start and end points. Then it creates a SurfaceOfRevolution object and
/// prints it's bounding box and axis and location. 
///
/// Input/Output:
///
/// The spline curve's file name is hardcoded to 'approj_curve.g2', and the
/// location and the axis direction are also hardcoded in the program.
///
/// Output is a file in Go-format for plot of the surface.
/// The file name is hard-coded to 'surface_of_revolution.g2'
///   
//===========================================================================

int main(int argc, char** argv)
{
    string inp_curve_filename("data/approj_curve.g2"); 
    Point location(0.0, 5.0, 200.0);
    Point axis_dir(1.0, 0.0, 0.0);
    cout << "\nRunning program '" << argv[0] << "' with spline curve filename= '"
	 << inp_curve_filename.c_str() << "'." << endl;


    // Read spline curve file
    ifstream cfile(inp_curve_filename.c_str());
    if (!cfile) {
	cerr << "\nFile error. Could not open file: " << inp_curve_filename.c_str()
	     << endl;
	return 1;
    }
    shared_ptr<SplineCurve> spline_curve(new SplineCurve);
    ObjectHeader header;
    cfile >> header;
    if (!header.classType() == SplineCurve::classType()) {
	THROW("Object type is NOT SplineCurve.");
    }
    cfile >> (*spline_curve);
    cfile.close();

    // Print some curve information
    Point pnt3d(3);
    spline_curve->point(pnt3d, spline_curve->startparam());     
    cout << "\nSplineCurve:  Dim= " << spline_curve->dimension()
	 << "\nStart.  Param= " << spline_curve->startparam() << "  Point= "
	 << pnt3d << endl;
    spline_curve->point(pnt3d, spline_curve->endparam());    
    cout << "End.  Param= " << spline_curve->endparam() << "  Point= "
	 << pnt3d << endl;
    cout << "Bounding box =   " << spline_curve->boundingBox() << endl;

    // Create the SurfaceOfRevolution object.
    SurfaceOfRevolution surf_of_revolution(location, axis_dir, spline_curve);
    cout << "\nSurface:  Dim= " << surf_of_revolution.dimension() << endl;
    cout << "Bounding box =  " << surf_of_revolution.boundingBox() << endl;
    cout << "Location =  " << surf_of_revolution.getLocation() << endl;
    cout << "Axis direction =  " << surf_of_revolution.getAxisDir() << endl;


    // Make a SplineSurface representation and write to file.
    SplineSurface* spline_surf = surf_of_revolution.geometrySurface();
    ofstream fout("surface_of_revolution.g2");
    spline_surf->writeStandardHeader(fout);
    spline_surf->write(fout);
    fout.close();
    // cout << "\nOpen the files 'surface_of_revolution.g2' and 'approj_curve.g2'"
    //      << " in 'goview' to look at the results.\n" << endl;
    delete spline_surf;

    return 0;
}
