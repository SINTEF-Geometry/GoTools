//===========================================================================
//                                                                           
// File: Intersections.
//                                                                           
// Created: Thu April 5 2001                                         
//                                                                           
// Author: Vibeke Skytt
//          
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GOINTERSECTIONS_H
#define _GOINTERSECTIONS_H

/// \file GoIntersections.h
/// Declaration file for a set of free intersection functions operating on
/// object belonging to GoTools but using the functionality of SISL.

#include "GoTools/geometry/SplineCurve.h"
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
