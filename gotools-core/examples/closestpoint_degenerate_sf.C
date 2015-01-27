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
// File: closestpoint_degenerate_sf.C
//                                                                           
/// Description:
///
/// This program demonstrates the use of the function
/// SplineSurface::closestPoint(const Point& pt, double& clo_u, double& clo_v, 
/// Point& clo_pt, double& clo_dist, double epsilon,
/// const RectDomain* domain_of_interest = NULL, double *seed = 0)
/// The declaration of the function is in 'ParamSurface.h'.
/// The function compute the closest point on a surface from a specified point. 
/// It reads a spline surface file in Go-format and a file in plain ASCII format
/// with the xyz-coordinates of the points we want to find the closest point on
/// the surface to.
///
///
/// Input/Output:
///
/// The program will read the the spline surface object from 'degenerate_sf.g2'
/// and the xyz-coordinates of the input points from
/// 'inp_degen_surf_close_points.dat'.
/// The program will write an output file 'surf_close_points.g2' and display
/// some informatation about the closest points on the screen.
/// The file 'degenerate_sf_close_points.g2' contains line segments from the input points
/// to the closest points on the surface in Go-format, and can together with
/// 'degenerate_sf.g2' be displayed with the program 'goview'.
///
/// The program will use the file 'degenerate_sf.g2'.
//
//===========================================================================

int main(int argc, char** argv)
{
    // Read the surface from a file in Go-format.
    string filename("data/degenerate_sf.g2");
    cout << "\nProgram " << argv[0] << " using file " << filename.c_str() << endl;
    ifstream file(filename.c_str());
    if (!file) {
	cerr << "\nFile error. Could not open file: " << filename.c_str() << endl;
	return 1;
    }
    ObjectHeader head;
    SplineSurface surf;
    file >> head;
    if (!head.classType() == SplineSurface::classType()) {
	THROW("Object type is NOT SplineSurface.");
    }
    file >> surf;
    file.close();

    // Read the points from a file. xyz-coordinates.
    string point_filename("data/inp_degen_surf_close_points.dat");
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
    
    cout << "\nProgram '" << argv[0] << "' using input files '" << filename.c_str()
	 << "' and '" << point_filename.c_str()
	 << ", and output file 'degen_surf_close_points.g2'." << endl;

    // Find the points on the surface closest to these points.
    double close_u;      // Closest point's u parameter.
    double close_v;      // Closest point's v parameter.
    Point  close_pt(3);  // Closest point's coordinates.
    double close_dist;   // Distance between the two points.
    double epsilon = 1e-8;  // Parameter tolerance

    // Write to file vectors from a point to the closest point on the surface.
    ofstream fout2("degenerate_sf_close_points.g2");
    // Class_LineCloud=410 MAJOR_VERSION=1 MINOR_VERSION=1 auxillary data=4
    // The four auxillary data values defines the colour (r g b alpha)
    fout2 << "410 1 0 4 255 0 0 255" << endl; // Header.
    fout2 << N << endl;

    // Find closest point using the whole surface. (The two last arguments
    // 'RectDomain* domain_of_interest' and 'double *seed' are by default
    // equal to 0).
    cout << "\nClosest points from inputfile points to points on the surface ";
    for (int i=0; i<N; ++i) {
	surf.closestPoint(points[i], close_u, close_v, close_pt, close_dist,
			  epsilon);
	fout2 << points[i] << ' ' <<  close_pt << endl;  // write vector
	cout << "Point: " << points[i] << "  Closest point: " << close_pt
	     << "\nParameter values= " <<  close_u << " , " <<  close_v
	     << "  Closest distance= " << close_dist << endl;
    }
    fout2.close();


    // Find closest point from points on the surface. Should be 0 + some tolerance.
    cout << "\nClosest points from points on the surface." << endl;
    const int nsp = 9;
    double du = (surf.endparam_u() - surf.startparam_u()) / (nsp-1);
    double dv = (surf.endparam_v() - surf.startparam_v()) / (nsp-1);
    cout << "Parameter u from " << surf.startparam_u() << " to " << surf.endparam_u()
	 << "  step "  << du << endl;
    cout << "Parameter v from " << surf.startparam_v() << " to " << surf.endparam_v()
	 << "  step "  << dv << endl;
    double max_dist = 0.0;
    Point point;
    for (double v=surf.startparam_v(); v<=surf.endparam_v(); v += dv) {
	for (double u=surf.startparam_u(); u<=surf.endparam_u(); u += du) {
	    surf.point(point, u, v);  // interpolate at u,v
	    surf.closestPoint(point, close_u, close_v, close_pt, close_dist,
			      epsilon);
#ifdef DEBUG
	    cout << "\n        Point: " << point << "\nClosest point: " << close_pt
		 << "\nParameter values= " <<  close_u << " , " <<  close_v
		 << "  Closest distance= " << close_dist << endl;
#endif
	}
	max_dist = std::max(close_dist, max_dist);
    }
    cout << "\nMaximum distance between an interpolated point and the "
	 << "corresponding input point is " << max_dist << '\n' << endl;

}
