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

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace std;
using namespace Go;


//===========================================================================
//                                                                           
// File: closestpoint_curve.C                                                   
//                                                                           
/// Description:
///
/// This program demonstrates the use of the two functions
/// ParamCurve::closestPoint(const Point& pt, double& clo_t, Point& clo_pt,
/// double& clo_dist)  and
/// SplineCurve::closestPoint(const Point& pt, double tmin, double tmax,
/// double& clo_t, Point& clo_pt, double& clo_dist, double const *seed = 0).
/// They compute the closest point on a curve from a specified point.
///
/// ParamCurve::closestPoint() takes the whole curve into account (not just an
/// interval of it). It starts searching at the start of the curve's parameter
/// values, and may return a point at a local distance minimum point instead of
/// the real closest point. Points at the end of a closed curve will get it's
/// closest point from the start of the curve.
///  
/// SplineCurve::closestPoint() may also find a local minimum, but has some 
/// additional arguments to restrict and govern the search. The parameters 
/// 'tmin' and 'tmax' restricts the search to parameter values between 
/// 'tmin' and 'tmax'. The parameter 'seed' may be provided by the
/// user to give an initial guess of the parameter value of the curves's closest
/// point.
///
/// The program will use the file 'interpol_curve2_hermite.curve.g2'.
/// If you don't have the file, run :
/// interpol_curve_hermite data/curve2.dat interpol_curve2_hermite.curve.g2 interpol_curve2_hermite.points.g2
///
//
//===========================================================================

int main(int argc, char** argv)
{
    // Read the curve from file
    string filename("interpol_curve2_hermite.curve.g2");
    cout << "\nProgram " << argv[0] << " using file " << filename.c_str() << endl;
    ifstream file(filename.c_str());
    if (!file) {
	cerr << "\nFile error. Could not open file: " << filename.c_str() << endl;
	return 1;
    }
    ObjectHeader head;
    SplineCurve cv;
    file >> head;
    if (!head.classType() == SplineCurve::classType()) {
	THROW("Object type is NOT SplineCurve.");
    }
    file >> cv;
    file.close();
 
    // Find the point on the curve closest to this point.
    Point point(0.0, 1.5, 0.0);
    double close_param;   // Closest point's parameter value.
    Point  close_pt(3);   // Closest point's coordinates.
    double close_dist;    // Distance between the two points.

    // Find closest point using the whole curve.
    ParamCurve* par_cv = static_cast<ParamCurve*>(&cv);
    par_cv->closestPoint(point, close_param, close_pt, close_dist);
    cout << "Closest point using the whole curve." << endl;
    cout << "Point: " << point << "  Closest point: " << close_pt
	 << "  Parameter value= " <<  close_param << "  Closest distance= "
	 << close_dist << endl;

    // Find closest point using the whole curve and a seed.
    double seed_param = cv.startparam() + 0.9*(cv.endparam() - cv.startparam());
    cv.closestPoint(point, cv.startparam(), cv.endparam(), close_param,
		    close_pt, close_dist,
		    &seed_param);
    cout << "\nClosest point using the whole curve and a seed.\n"
	 << "Point: " << point << "  Closest point: " << close_pt
	 <<  "  Parameter value= " <<  close_param <<  "  Closest distance= "
	 << close_dist << endl;

    // Find closest point searching from 'seed_param' to  cv.endparam().
    cv.closestPoint(point, seed_param, cv.endparam(), close_param, close_pt,
		    close_dist);
    cout << "\nClosest point  searching from 'seed_param' to  'cv.endparam()'."
	 << " (" << seed_param << " to " << cv.endparam() << ")\n"
	 << "Point: " << point << "  Closest point: " << close_pt
	 <<  "  Parameter value= " <<  close_param <<  "  Closest distance= "
	 << close_dist << endl;
}










