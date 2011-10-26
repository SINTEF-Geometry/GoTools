//===========================================================================
//                                                                           
// File: closestPtSurfSurfPlane.h
//                                                                           
// Created:
//                                                                           
// Author: B. Spjelkavik <bsp@sintef.no>
//          
// Revision:  $Id: closestPtSurfSurfPlane.h,v 1.4 2005-07-05 12:29:39 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _CLOSESTPTSURFSURFPLANE_H
#define _CLOSESTPTSURFSURFPLANE_H

/** @file closestPtSurfSurfPlane.h
 * Declaration file for an intersection function operating on 
 * object of class ParamCurve belonging to geometry.
 * Ported from the sisl function s9iterate.
 */


#include "GoTools/geometry/ParamSurface.h"

namespace Go {
    /** \enum AlgorithmChoice This enumeration denotes whether the 'closestPtSurfSurfPlane' 
     *  function should base itself upon a geometrical algorithm (marching on the surfaces) or
     *  a pure functional one (constrained minimization of a cost function).
     *  The default has been set to FUNCTIONAL, as there are flaws in the GEOMETRICAL
     *  one when it comes to close-to-degenerate cases.  
     */
    enum AlgorithmChoice {GEOMETRICAL, FUNCTIONAL};

  /** Newton iteration on the distance function between
   * two surfaces and a plane to find a closest point or an intersection point.
   * Ported from the sisl-function \c s9iterate.
   * \param epoint - Vector of Points containing parts of plane description.
   * epoint[0] contains a point in the plane.
   * epoint[1] contains the normal vector to the plane.
   * \param epnt1  0-2 Derivatives + normal of start point for iteration in
   * first surface.
   * \param epnt2  0-2 Derivatives + normal of start point for iteration in
   * second surface.
   * \param epar1 Parameter pair of start point in first surface.
   * \param epar2 Parameter pair of start point in second surface.
   * \param psurf1 Pointer to the first surface.
   * \param psurf2 Pointer to the second surface.
   * \param aepsge Absolute tolerance.
   *
   * \param  gpnt1  0-2 Derivatives + normal of result of iteration in 
   * first surface.
   * \param  gpnt2  0-2 Derivatives + normal of result of iteration in 
   * second surface.
   * \param gpar1 Parameter pair of result of iteration in first surface.
   * \param gpar2 Parameter pair of result of iteration in second surface.
   * \param jstat =1 Intersection. =2 Minimum value.  =3 Nothing found.
   * \param algo choose whether the implementation of the algorithm should be
   *        based on general function minimization, or a geometric approach
   */

void closestPtSurfSurfPlane(const std::vector<Point>& epoint,
			    const std::vector<Point>& epnt1,
			    const std::vector<Point>& epnt2,
			    const Point& epar1,
			    const Point& epar2,
			    const ParamSurface* psurf1,
			    const ParamSurface* psurf2,
			    double aepsge,
			    std::vector<Point>& gpnt1,
			    std::vector<Point>& gpnt2,
			    Point& gpar1, 
			    Point& gpar2, 
			    int& jstat,
			    AlgorithmChoice algo = FUNCTIONAL);


/// This is the geometrically-based implementation of closestPtSurfSurfPlane.
/// It is called from the function above (\ref closestPtSurfSurfPlane()) if
/// it is called with algo = GEOMETRICAL.  The user can also choose to call
/// this function directly.   The parameter description is the same as in
/// \ref closestPtSurfSurfPlane().
void closestPtSurfSurfPlaneGeometrical(const std::vector<Point>& epoint,
				       const std::vector<Point>& epnt1,
				       const std::vector<Point>& epnt2,
				       const Point& epar1,
				       const Point& epar2,
				       const ParamSurface* psurf1,
				       const ParamSurface* psurf2,
				       double aepsge,
				       std::vector<Point>& gpnt1,
				       std::vector<Point>& gpnt2,
				       Point& gpar1, Point& gpar2, int& jstat);

/// This is the functional-based implementation of closestPtSurfSurfPlane.
/// It is called from the function above (\ref closestPtSurfSurfPlane()) if
/// it is called with algo = FUNCTIONAL.  The user can also choose to call
/// this function directly.   The parameter description is the same as in
/// \ref closestPtSurfSurfPlane().
void 
closestPtSurfSurfPlaneFunctional(const std::vector<Point>& epoint, //plane description
				 const std::vector<Point>& epnt1, // start pt. in surf. 1
				 const std::vector<Point>& epnt2, // start pt. in surf. 2
				 const Point& epar1, // parameter start pt. in surf. 1
				 const Point& epar2, // parameter start pt. in surf. 2
				 const ParamSurface* psurf1, // ptr. to surf. 1
				 const ParamSurface* psurf2, // ptr. to surf. 2
				 double aepsge, // absolute tolerance
				 std::vector<Point>& gpnt1, // result of iter. in surf. 1
				 std::vector<Point>& gpnt2, // result of iter. in surf. 2
				 Point& gpar1, Point& gpar2, int& jstat); // results of param.


// Anonymous namespace
namespace {
  /** Computes the next step values along the parameters of a surface in an
   * iteration to an intersection point between two surfaces and a plane.
   * \param fpnt 0-2 Derivatives + normal of start point for iteration in
   * first surface.
   * \param gpnt 0-2 Derivatives + normal of start point for iteration in
   * second surface.
   * \param snorm Normal vector of the intersection plane.
   * \param spoint Point in the intersection plane.
   *
   * \param delta New step values for first surface.
   * \param kstat =1 The equation system is singular.  =0 The equation system
   * is not singular.
   */
void nextStep(const std::vector<Point>& fpnt,const std::vector<Point>& gpnt,
	      const Point& snorm, const Point& spoint,
	      Point& delta, int& kstat);

} // Anonymous namespace

} // namespace Go  

#endif //  _CLOSESTPTSURFSURFPLANE_H
