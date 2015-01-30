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
#include "GoTools/geometry/SplineUtils.h"

using namespace std;

namespace Go {

//==========================================================================
void SplineSurface::normal(Point& pt, double upar, double vpar) const
//==========================================================================
{
    bool succeeded = false;
    if (dim_ == 2)
      {
	succeeded = true;
	pt = Point(0.0, 0.0, 1.0);
      }
    else
      {
	try {
	  succeeded = normal_not_failsafe(pt, upar, vpar);
	} catch ( ... ) {
	  //MESSAGE("Failed finding normal, trying a new method.");
	}
      }
    if (!succeeded) {
	/*MESSAGE("Could not compute normal in standard way - surface too degenerate."
	  "  Using closest found normal instead.");*/
	// doing a search for closest defined normal
	vector<Point> neighbourhood_normals;
	neighbourhood_normals.reserve(4);
	vector<bool> interval_for_u;
	vector<pair<double, double> > param_intervals;
	param_intervals.reserve(4);
	// Determine in which directions we can search for a valid normal
	if (upar - DEFAULT_PARAMETER_EPSILON > startparam_u()) {
	    param_intervals.push_back(pair<double, double>(upar, startparam_u()));
	    interval_for_u.push_back(true);
	}
	if (upar + DEFAULT_PARAMETER_EPSILON < endparam_u()) {
	    param_intervals.push_back(pair<double, double>(upar, endparam_u()));
	    interval_for_u.push_back(true);
	}
	if (vpar - DEFAULT_PARAMETER_EPSILON > startparam_v()) {
	    param_intervals.push_back(pair<double, double>(vpar, startparam_v()));
	    interval_for_u.push_back(false);
	} 
	if (vpar + DEFAULT_PARAMETER_EPSILON < endparam_v()) {
	    param_intervals.push_back(pair<double, double>(vpar, endparam_v()));
	    interval_for_u.push_back(false);
	}

	// searching for neighbourhood normals
	for (int interv = 0; interv < int(param_intervals.size()); ++interv) {
	    // binary search in the current interval
	    Point p;
	    double fixed_val = interval_for_u[interv] ? vpar : upar;
	    bool found = search_for_normal(interval_for_u[interv],
					   fixed_val,
					   param_intervals[interv].first,
					   param_intervals[interv].second,
					   p);
	    if (found) {
		neighbourhood_normals.push_back(p);
	    }
	}
	if (neighbourhood_normals.size() == 0) {
	    pt.setValue(0.0);
	    /*MESSAGE("Could not find any defined normal on surface.  "
	      "Setting normal to 0 instead.");*/
	} else {
	    // taking average of detected normals
	    pt = neighbourhood_normals[0];
	    for (int i = 1; i < int(neighbourhood_normals.size()); ++i) {
		pt += neighbourhood_normals[i];
	    }
	    double l = pt.length();
	    if (l < DEFAULT_SPACE_EPSILON) {
		/*MESSAGE("Neighbouring normals cancel each other out.  "
		  "Setting normal to 0 instead.")*/;
	    }
	    pt /= l;
	}
    }
}

//==========================================================================
bool SplineSurface::search_for_normal(bool interval_in_u,
				      double fixed_parameter,
				      double interval_start, // normal is not defined here
				      double interval_end,
				      Point& normal) const
//==========================================================================
{
    // We want to find the closest defined normal to the point on the surface with
    // the parameter values (interval_start, fixed_parameter) if interval_in_u = true,
    // or to the parameter values (fixed_parameter, interval_start) if interval_in_u = false.
    // We suppose that the normal for that point has been detected as undefined.
    
    const int INTERVAL_PARTITION = 10;
    double running_parameter = interval_start;
    double& u = interval_in_u ? running_parameter : fixed_parameter;
    double& v = interval_in_u ? fixed_parameter : running_parameter;

    // first a quick check - if we can find a normal in immediate vincinity to the
    // point (within the DEFAULT_PARAMETER_EPSILON value), we will return it.
    if (interval_end - interval_start > 0) {
	running_parameter += DEFAULT_PARAMETER_EPSILON;
    } else {
	running_parameter -= DEFAULT_PARAMETER_EPSILON;
    }
    bool normal_defined = false;
    try { 
      normal_defined = normal_not_failsafe(normal, u, v);
    } catch ( ... ) {
      MESSAGE("Failed finding normal, trying a new method.");
    }
    if (normal_defined) {
	return true;
    }
    // if we got here, our little 'wiggling' of the moving parameter did not produce
    // a valid normal.  We will do a more thorough search.

    // trying to locate a defined normal at some moderate distance from the 
    // "difficult" point
    double stepsize = (interval_end - interval_start) / double(INTERVAL_PARTITION);
    int i;
    for (i = 1; i <= INTERVAL_PARTITION; ++i) {
	running_parameter = interval_start + i * stepsize;
	try { 
	  normal_defined = normal_not_failsafe(normal, u, v);
	} catch ( ... ) {
	  MESSAGE("Failed finding normal, trying a new method.");
	}
	if (normal_defined) { 
	    break;
	}
    }
    if (i > INTERVAL_PARTITION) {
	// we didn't find any defined normal in this interval at all
	return false;
    }

    // restricting interval
    interval_end = running_parameter; // as detected in the above loop
    // trying to find the value of running_parameter that is as close as possible
    // to interval_start while still with a defined normal.
    Point tentative_normal;
    while (fabs(running_parameter - interval_start) > DEFAULT_PARAMETER_EPSILON) {
	running_parameter = (interval_start + interval_end) * 0.5;
	bool tentative_normal_defined = false;
	try { 
	    tentative_normal_defined =
		normal_not_failsafe(tentative_normal, u, v);
	} catch ( ... ) {
	    MESSAGE("Failed finding normal, trying a new method.");
	}
	if (!tentative_normal_defined) {
	    // we have gotten so close to the original point that the normal
	    // has become undefined.  Keep previous normal
	    break;
	} else {
	    interval_end = running_parameter; // shrinking interval
	    normal = tentative_normal; // this normal was defined, and closer to orig. point
	}
    }
    
    // we have now found a defined normal that should be as close to the original point as
    // possible
    return true;
}

//==========================================================================
bool SplineSurface::normal_not_failsafe(Point& pt, double upar, double vpar) const
//==========================================================================
{
    double tol = DEFAULT_SPACE_EPSILON;

#ifdef _OPENMP
    vector<Point> derivs(3, Point(1.0, 1.0, 1.0));
#else
    static vector<Point> derivs(3, Point(1.0, 1.0, 1.0));
#endif
    point(derivs, upar, vpar, 1);
    //    vector<Point> derivs = ParamSurface::point(upar, vpar, 1);

    //    Point &der1=derivs[1];
    //    Point &der2=derivs[2];
    pt.resize(3);
    pt.setToCrossProd(derivs[1], derivs[2]);
    double l = pt.length();
    int iterations = 0;
    // @@sbr But this depends on the parametrization!!! Why not use
    // angle between cross tangents as a guiding line?
    double cross_tan_ang = derivs[1].angle_smallest(derivs[2]);
    cross_tan_ang = min(cross_tan_ang, fabs(M_PI - cross_tan_ang));
    double cross_tan_ang_tol = 1e-03;
    while (l < tol) {
	// Surface is degenerate at that point. One of the derivatives
	// must be small.
	int okderindex = 1;
	if (derivs[1].length2() < derivs[2].length2())
	    okderindex = 2;
	Point okder = derivs[okderindex];
// #ifdef _MSC_VER
// 	const RectDomain& rd 
// 	  = dynamic_cast<const RectDomain&>(parameterDomain());
// #else
	const RectDomain& rd = parameterDomain();
// #endif
	double lowx = okderindex==2 ? rd.vmin() : rd.umin();
	double x = okderindex==2 ? vpar : upar;
	bool xislow = false;
	if (abs(x-lowx) < tol) {
	    xislow = true;
	}
	double lowt = okderindex==2 ? rd.umin() : rd.vmin();
	double hight = okderindex==2 ? rd.umax() : rd.vmax();
	double t = okderindex==2 ? upar : vpar;
	double newt = lowt;
	bool lowchoice = true;
	if (hight - t > t - lowt) {
	    newt = hight;
	    lowchoice = false;
	}
	if (okderindex==2) {
	    upar = newt;
	} else {
	    vpar = newt;
	}
	point (derivs, upar, vpar, 1);
	bool okderfirst = true;
	if (okderindex==2) {
	    okderfirst = (lowchoice == xislow);
	} else {
	    okderfirst = (lowchoice != xislow);
	}
	pt = okderfirst ?
	    okder % derivs[okderindex] : derivs[okderindex] % okder;
	l = pt.length();
	++iterations;
	if (iterations == 2) break;
    }
    if (l < tol || cross_tan_ang < cross_tan_ang_tol) {
//         if (cross_tan_ang < cross_tan_ang_tol)
// 	    MESSAGE("Too small angle between cross tangents, "
// 		    "degenerate point!");
	pt.setValue(0.0);
	return false;
    }
    pt /= l;
    return true;
}

//==========================================================================
SplineSurface* SplineSurface::normal() const
//==========================================================================
{
    THROW("Not yet implemented!");

}

// //==========================================================================
// void SplineSurface::normal(Point& pt, double upar, double vpar) const
// //==========================================================================
// {
//     double tol = DEFAULT_SPACE_EPSILON;

//     static vector<Point> derivs(3, Point(1.0, 1.0, 1.0));
//     point(derivs, upar, vpar, 1);
//     //    vector<Point> derivs = ParamSurface::point(upar, vpar, 1);

//     Point &der1=derivs[1];
//     Point &der2=derivs[2];
//     pt.resize(3);
//     pt.setToCrossProd(derivs[1], derivs[2]);
//     double l = pt.length();
//     int iterations = 0;
//     while (l < tol) {
// 	// Surface is degenerate at that point. One of the derivatives
// 	// must be small.
// 	int okderindex = 1;
// 	if (derivs[1].length2() < derivs[2].length2())
// 	    okderindex = 2;
// 	Point okder = derivs[okderindex];
// // #ifdef _MSC_VER
// // 	const RectDomain& rd 
// // 	  = dynamic_cast<const RectDomain&>(parameterDomain());
// // #else
// 	const RectDomain& rd = parameterDomain();
// // #endif
// 	double lowt = okderindex==1 ? rd.umin() : rd.vmin();
// 	double hight = okderindex==1 ? rd.umax() : rd.vmax();
// 	double t = okderindex==1 ? upar : vpar;
// 	double newt = lowt;
// 	if (hight - t > t - lowt)
// 	    newt = hight;
// 	if (okderindex==1)
// 	    upar = newt;
// 	else
// 	    vpar = newt;
// 	point (derivs, upar, vpar, 1);
// 	pt = derivs[okderindex] % okder;
// 	l = pt.length();
// 	++iterations;
// 	if (iterations == 2) break;
//     }
//     MESSAGE_IF (l < tol,
// 		   "Could not compute normal in standard way - surface too degenerate."
// 		   "  Using pseudo-normal instead.");
//     if (l < tol) {
// 	// trying to calculate a "pseudo-normal"
// 	Point location;
// 	point(location, upar, vpar);
// 	// since we are in a degenerate point, we suppose that there exists a control
// 	// point whose location is exactly here
// 	int corner_1;
// 	BoundingBox bbox = boundingBox();
// 	double bbox_span = bbox.high().dist(bbox.low());

// 	if (rational_) {
// 	    MESSAGE("Calculation of pseudo-normal not implemented for rational surfaces."
// 		       " Returning zero normal.");
// 	    pt.setValue(0.0);
// 	    return;
	    
// 	}

// 	int num_points = coefs_.size() / 3;
// 	for (corner_1 = 0; corner_1 != num_points; ++corner_1) {
// 	    double dx = coefs_[3 * corner_1 + 0] - location[0];
// 	    double dy = coefs_[3 * corner_1 + 1] - location[1];
// 	    double dz = coefs_[3 * corner_1 + 2] - location[2];
// 	    double dist = sqrt(dx * dx + dy * dy + dz * dz);
// 	    if (dist < tol) {
// 		// we suppose that we have found a control point lying at the exact
// 		// point on the surface that we are examining
// 		break;
// 	    }
// 	}
// 	if (corner_1 == num_points) {
// 	    MESSAGE("Unable to calculate pseudo-normal.  Returning zero normal.");
// 	    pt.setValue(0.0);
// 	    return;
// 	}
// 	// Now we'll have to try to define a triangle of control points that could
// 	// be used to define a normal
// 	int corner_2 = -1;
// 	double min_dist = 2 * bbox_span;
// 	for (int i = 0; i < num_points; ++i) {
// 	    double dx = coefs_[3 * i + 0] - coefs_[3 * corner_1 + 0];
// 	    double dy = coefs_[3 * i + 1] - coefs_[3 * corner_1 + 1];
// 	    double dz = coefs_[3 * i + 2] - coefs_[3 * corner_1 + 2];
// 	    double dist = sqrt(dx * dx + dy * dy + dz * dz);
// 	    if (dist < min_dist && dist > tol) {
// 		min_dist = dist;
// 		corner_2 = i;
// 	    }
// 	}
// 	if (corner_2 == -1) {
// 	    MESSAGE("Unable to calculate pseudo-normal.  Returning zero normal.");
// 	    pt.setValue(0.0);
// 	    return;
// 	}
// 	// looking for third corner
// 	int corner_3 = -1;
// 	min_dist = 2 * bbox_span;
// 	for (int i = 0; i < num_points; ++i) {
// 	    double dx = coefs_[3 * i + 0] - coefs_[3 * corner_1 + 0];
// 	    double dy = coefs_[3 * i + 1] - coefs_[3 * corner_1 + 1];
// 	    double dz = coefs_[3 * i + 2] - coefs_[3 * corner_1 + 2];
// 	    double dist_1 = sqrt(dx * dx + dy * dy + dz * dz);
// 	    if (dist_1 < min_dist && dist_1 > tol) {
// 		dx = coefs_[3 * i + 0] - coefs_[3 * corner_2 + 0];
// 		dy = coefs_[3 * i + 1] - coefs_[3 * corner_2 + 1];
// 		dz = coefs_[3 * i + 2] - coefs_[3 * corner_2 + 2];
// 		double dist_2 = sqrt(dx * dx + dy * dy + dz * dz);
// 		if (dist_2 > tol) {
// 		    min_dist = dist_1;
// 		    corner_3 = i;
// 		}
// 	    }
// 	}
// 	if (corner_3 == -1) {
// 	    MESSAGE("Unable to calculate pseudo-normal.  Returning zero normal.");
// 	    pt.setValue(0.0);
// 	    return;
// 	}
// 	Point T1(&coefs_[3 * corner_2], &coefs_[3 * (corner_2 + 1)]);
// 	Point T2(&coefs_[3 * corner_3], &coefs_[3 * (corner_3 + 1)]);
// 	pt = T2.cross(T1);
// 	l = pt.length();
//     } 
//     pt /= l;

// }  



} // namespace Go



