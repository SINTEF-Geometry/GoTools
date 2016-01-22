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

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace std;
using namespace Go;


//===========================================================================
//                                                                           
// File: closestpoint_surface.C                                                
///                                                                           
/// Description:
///
/// This program demonstrates the use of the function
/// 'SplineSurface::closestPoint(const Point& pt, double& clo_u, double& clo_v, 
/// Point& clo_pt, double& clo_dist, double epsilon,
/// const RectDomain* domain_of_interest = NULL, double *seed = 0)'.
/// The declaration of the function is in 'ParamSurface.h'.
/// The function compute the closest point on a surface from a specified point. 
/// It reads a spline surface file in Go-format and a file in plain ASCII format
/// with the xyz-coordinates of the points we want to find the closest point on
/// the surface to.
///
/// Input/Output:

/// The program will read the the spline surface object from 'surface.g2' and
/// the xyz-coordinates of the input points from 'inp_surf_close_points.dat'.
/// The program will write an output file 'surf_close_points.g2' and display
/// some informatation about the closest points on the screen.
/// The file 'surf_close_points.g2' contains line segments from the input points
/// to the closest points on the surface in Go-format, and can together with
/// 'surface.g2' be displayed with the program 'goview'.
///
//===========================================================================

int main(int argc, char** argv)
{
    // Read the surface from a file in Go-format.
    string surf_filename("surface.g2");
    ifstream sfile(surf_filename.c_str());
    if (!sfile) {
	cerr << "\nFile error. Could not open file: " << surf_filename.c_str() << endl;
	return 1;
    }
    ObjectHeader head;
    SplineSurface surf;
    sfile >> head;
    if (!(head.classType() == SplineSurface::classType())) {
	THROW("Object type is NOT SplineSurface.");
    }
    sfile >> surf;
    sfile.close();

    // Read the points from a file. xyz-coordinates.
    string point_filename("inp_surf_close_points.dat");
    ifstream pfile(point_filename.c_str());
    if (!pfile) {
	cerr << "\nFile error. Could not open file: " << point_filename.c_str() << endl;
	return 1;
    }
    vector<Point> points;
    while (1) {
	Point p(3);
	pfile >> p;
	if (!pfile) break;
	points.push_back(p);
    }
    pfile.close();
    int N = (int)points.size();
    
    cout << "\nProgram '" << argv[0] << "' using input files '" << surf_filename.c_str()
	 << "' and '" << point_filename.c_str()
	 << ", and output file 'surf_close_points.g2'." << endl;
    cout << "Parameter u from " << surf.startparam_u() << " to "
	 << surf.endparam_u() << endl;
    cout << "Parameter v from " << surf.startparam_v() << " to "
	 << surf.endparam_v() << endl;

 
    // Find the points on the surface closest to these points.
    double close_u;      // Closest point's u parameter.
    double close_v;      // Closest point's v parameter.
    Point  close_pt(3);  // Closest point's coordinates.
    double close_dist;   // Distance between the two points.
    double epsilon = 1e-8;  // Parameter tolerance

    // Write to file vectors from a point to the closest point on the surface.
    ofstream fout2("surf_close_points.g2");
    // Class_LineCloud=410 MAJOR_VERSION=1 MINOR_VERSION=1 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout2 << "410 1 0 4 255 0 0 255" << endl; // Header.
    fout2 << N << endl;

    // Find closest point using the whole surface. (The two last arguments
    // 'RectDomain* domain_of_interest' and 'double *seed' are by default
    // equal to 0).
    cout << "\nClosest point using the whole surface.";
    for (int i=0; i<N; ++i) {
	surf.closestPoint(points[i], close_u, close_v, close_pt, close_dist,
			  epsilon);
	fout2 << points[i] << ' ' <<  close_pt << endl;  // write vector
	cout << "\nPoint: " << points[i] << "  Closest point: " << close_pt
	     << "\nParameter values= " <<  close_u << " , " <<  close_v
	     << "  Closest distance= " << close_dist << endl;
    }
    fout2.close();

    // Find closest point. Restrict the search to the lowest half of the
    // u and v parameter range.

    double u1 = surf.startparam_u();
    double u2 = 0.5*(surf.startparam_u() + surf.endparam_u());
    double v1 = surf.startparam_v();
    double v2 = 0.5*(surf.startparam_v() + surf.endparam_v());
    Array<double, 2> corner1(u1, v1);
    Array<double, 2> corner2(u2, v2);
    RectDomain param_domain(corner1, corner2);
    cout << "\nSurface parameters U-domain: " << u1 << " - " << surf.endparam_u()
	 << "\tV-domain: " << v1 << " - " << surf.endparam_v() << endl;    
    cout << "Restricted surface parameters u-domain: " << corner1[0] << " - "
	 << corner2[0] << " \tv-domain: " << corner1[1] << " - " << corner2[1]
	 << endl;
    cout << "Restricted rectangle corners. Lower left "
	 << param_domain.lowerLeft() << " \tUpper right "
	 << param_domain.upperRight() << endl;

    cout << "\nClosest points in the lowest half of the u and v parameter range.";
    N=2;
    points.resize(N);
    points[0] = Point(70,10,300);    // Outside restricted area.
    points[1] = Point(30,-30,170);   // Inside restricted area. 
    for (int i=0; i<N; ++i) {
	surf.closestPoint(points[i], close_u, close_v, close_pt, close_dist,
			  epsilon, &param_domain);
	cout << "\nPoint: " << points[i] << "  Closest point: " << close_pt
	     << "\nParameter values= " <<  close_u << " , " <<  close_v
	     << "  Closest distance= " << close_dist << endl;
	if (close_u < u1 || close_u > u2 || close_v < v1 || close_v > v2) {
	    cout << "Error??? u,v-parameter values are outside restricted area!"
		 << endl;
	}
    }
    cout << endl;

}
