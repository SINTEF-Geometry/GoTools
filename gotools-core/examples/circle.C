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

#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>

using namespace std;
using namespace Go;


//===========================================================================
//                                                                           
// File: circle.C                                                   
//                                                                           
/// Description:
///
/// This program demonstrates the use of the class Circle.
/// It is a subclass of ElementaryCurve.
/// The space dimension of a circle is either 2 or 3.
/// The default parametrization is an angle from 0 to 2*PI.
///
/// This program constructs a full circle in the xy-plane and computes and
/// prints the circle's length, bounding box and the closest point on the
/// circle to a given point, and the closest point when the parameter search
/// area is bounded.
/// Then it creates a circle segment (from 120 to 210 degrees) and computes
/// and prints the circle segments's length, bounding box and closest point.
/// Then it makes a SplineCurve representation of the circle segment. It's
/// endpoints will have the same parameter values as the circle's segment.
/// Finally the SplineCurve circle segment is written to the file
/// "spline_circle_segm.g2".
///
//===========================================================================

int main(int argc, char** argv)
{
    cout << "\nRunning program " << argv[0] << "\nNo input arguments." << endl;
    // Construct a full circle in the xy-plane.
    // The constructor takes the radius, the centre, the normal to the
    // plane of the circle, and the (approximate) direction of the
    // (local) x-axis as input.
    double twopi = 2.0 * M_PI;
    double epsilon = 1.0e-10;
    double radius = 1.0;
    Point centre(0.0, 0.0, 0.0);
    Point x_axis(1.0, 0.0, 0.0);  // Direction vector
    Point z_axis(0.0, 0.0, 1.0);  // Direction vector
    Point normal = z_axis;

    cout << "\n*** Full Circle ***" << endl;
    Circle circle1(radius, centre, normal, x_axis);
    cout << "Start parameter = " << circle1.startparam() << ".  "
	 << "End   parameter = " << circle1.endparam() << endl;
    cout << "Circle length   = " << circle1.length(epsilon) << endl;
    cout << "Bounding box    = " << circle1.boundingBox() << endl;

    // Closest point(). Full circle. Parameter search from 0 to 2*PI 
    // double tmin, tmax, clo_t, clo_u, clo_v, clo_dist;
    double tmin, tmax, clo_t, clo_dist;
    Point pnt, clo_pt;
    pnt = Point(0.5, -0.5, 0.0);
    tmin = 0.0;
    tmax = twopi;
    circle1.closestPoint(pnt, tmin, tmax, clo_t, clo_pt, clo_dist);
    cout << "\nclosestPoint() between parameters " << tmin << " and " << tmax
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nt = "    << clo_t << "   distance = "     << clo_dist << endl;

    // Closest point(). Parameter search from 120 to 210 deg.
    tmin = 4*twopi/12;    tmax = 7*twopi/12;
    circle1.closestPoint(pnt, tmin, tmax, clo_t, clo_pt, clo_dist);
    cout << "\nclosestPoint() between parameters " << tmin << " and " << tmax
	 << "\nPoint= " << pnt   << "   Closest point= " << clo_pt
	 << "\nt = "    << clo_t << "   distance = "     << clo_dist
	 << "\nThis closest point is at the end of the search area" << endl;


    // Circle segment
    Circle circle_segm(circle1);  // Default copy constructor
    circle_segm.setParamBounds(tmin, tmax);
    cout << "\n*** Circle segment ***" << endl;
    cout << "Start parameter = " << circle_segm.startparam() << ".  "
	 << "End parameter = "   << circle_segm.endparam() << endl;
    cout << "Circle length   = " << circle_segm.length(epsilon) << endl;
    cout << "Bounding box    = " << circle_segm.boundingBox() << "  ??? @bsp"
	 << endl;

    // Spline representation of a circle.
    SplineCurve* spline_circle_segm = circle_segm.geometryCurve();
    Point start_point, end_point;
    spline_circle_segm->point(start_point, tmin);    
    spline_circle_segm->point(end_point, tmax);    
    cout << "\n*** Spline curve circle segment ***" << endl;
    cout << "Start parameter = " << spline_circle_segm->startparam() << ".  "
	 << "Start point = "     << start_point << ".  " << endl
	 << "End parameter =   " << spline_circle_segm->endparam() << ".  "
	 << "End point = "       << end_point << endl;
    cout << "Curve length    = " << spline_circle_segm->length(epsilon) << endl;
    cout << "Bounding box    = " << spline_circle_segm->boundingBox() << " ??? @bsp"
	 << endl;

    // Write subcurve to file
    ofstream fout("spline_circle_segm.g2");
     // Class_SplineCurve=100 MAJOR_VERSION=1 MINOR_VERSION=1 auxillary_data=0
    spline_circle_segm->writeStandardHeader(fout); // write header.
    fout << *spline_circle_segm;    // write spline curve data.
    fout.close();
    // cout << "Open the file 'spline_circle_segm.g2' in 'goview' to look at the results"
    //      << endl;
    delete spline_circle_segm;
}










