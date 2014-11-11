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

#ifndef _GOINTERSECTIONS_H
#define _GOINTERSECTIONS_H

/// \file GoIntersections.h
/// Declaration file for a set of free intersection functions operating on
/// object belonging to GoTools but using the functionality of SISL.

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <utility>

namespace Go
{

/// Enumeration of various pretopologies that can be associated with
/// intersections.
enum {  pretop_UNDEF,
	pretop_IN,  
	pretop_OUT, 
	pretop_ON,  
	pretop_AT };


// Intersect a curve with a point
void intersectCurvePoint(const ParamCurve* crv, Point pnt, double epsge,
			 std::vector<double>& intersections,
			 std::vector<std::pair<double, double> >& int_crvs);

/// Intersect two 2D spline curves. Collect intersection parameters
/// and pretopology information.
/// \param cv1 pointer to the first 2D curve 
/// \param cv2 pointer to the second 2D curve
/// \param epsge geometrical tolerance
/// \retval intersections this vector will contain the parameter pairs
///         of the found intersections (one vector entry per intersection.
///         The two parameters in the pair<> correspond to the parameter value
///         in 'cv1' and 'cv2' for a particular intersection.
/// \retval pretopology vector containing a pretopology indicator for each
///                     detected intersection point.  There is one entry per
///                     intersection point.
void intersect2Dcurves(const ParamCurve* cv1, 
		       const ParamCurve* cv2, 
		       double epsge,
		       std::vector<std::pair<double,double> >& intersections,
		       std::vector<int>& pretopology,
		       std::vector<std::pair<std::pair<double,double>, 
		       std::pair<double,double> > >& int_crvs);

///Intersect two spline curves. Collect intersection parameters.
/// \param cv1 pointer to the first spline curve
/// \param cv2 pointer to the second spline curve
/// \param epsge geometrical tolerance
/// \retval intersections this vector will contain the parameter pairs
///         of the found intersections (one vector entry per intersection.
///         The two parameters in the pair<> correspond to the parameter value
///         in 'cv1' and 'cv2' for a particular intersection.
void intersectcurves(SplineCurve* cv1, SplineCurve* cv2, double epsge,
		     std::vector<std::pair<double,double> >& intersections);

/// Intersect a spline curve and a spline surface
 void intersectCurveSurf(const SplineCurve *cv, const SplineSurface *sf,
			 double epsge, 
			 std::vector<std::pair<double, Point> >& int_pts,
			 std::vector<int>& pretopology,
			 std::vector<std::pair<std::pair<double,Point>, 
			 std::pair<double,Point> > >& int_crvs);


/// Compute the closest point between two curves
/// \param cv1 pointer to the first curve
/// \param cv2 pointer to the second curve
/// \param epsge geometrical tolerance
/// \retval par1 parameter of the closest point in the first curve
/// \retval par2 parameter of the closest point in the second curve
/// \retval dist distance between the curves at the closest point
void closestPtCurves(SplineCurve* cv1, SplineCurve* cv2, double epsge,
		     double& par1, double& par2, double& dist);

} // namespace Go


#endif // _GOINTERSECTIONS_H
