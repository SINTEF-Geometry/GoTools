//===========================================================================
//                                                                           
// File: ClosestPoint.h                                                      
//                                                                           
// Created: Mon Feb 13 13:36:38 2012                                         
//                                                                           
// Author: Bjørn Spjelkavik <bsp@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
//===========================================================================

#ifndef _CLOSESTPOINT_
#define _CLOSESTPOINT_


#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"


namespace Go
{
    /** \enum AlgorithmChoice This enumeration denotes whether the 'closestPtSurfSurfPlane' 
     *  function should base itself upon a geometrical algorithm (marching on the surfaces) or
     *  a pure functional one (constrained minimization of a cost function).
     *  The default has been set to FUNCTIONAL, as there are flaws in the GEOMETRICAL
     *  one when it comes to close-to-degenerate cases.  
     */
    enum AlgorithmChoice {GEOMETRICAL, FUNCTIONAL};

/// Namespace for computing closest points.
namespace ClosestPoint
{

  /** Newton iteration on the distance function between
   * two curves in 2D to find a closest point or an intersection point.
   * The solution is restricted to the intervals [tmin1,tmax1] and 
   * [tmin2,tmax2] in the first and second curves's parameter domain.
   * Ported from the sisl-function \c s1770_2D.
   * \param cv1   Curve number one.
   * \param cv2   Curve number two.
   * \param aepsge Geometry resolution.
   * \param tmin1 Min. parameter value. First curve  
   * \param tmax1 Max. parameter value. First curve
   * \param tmin2 Min. parameter value. Second curve
   * \param tmax2 Max. parameter value. Second curve
   * \param seed1 Start point for iteration along curve number one.
   * \param seed2 Start point for iteration along curve number two.
   * \param method 1 or 2. Using the the method in s1770_2D, starting
   * with order one or two. method=3 uses the method from s1770 (order one).
   * \param quick Reduce requirement on exactness.
   * \retval par1 Parameter value of the first curve's point.
   * \retval par2 Parameter value of the second curve's point.
   * \retval dist Distance in space between the points.
   * \retval ptc1 Point on curve number one.
   * \retval ptc2 Point on curve number two.
   * \retval istat =1 Intersection. =2 Minimum value.  =3 Nothing found.
   */
void closestPtCurves2D(ParamCurve* cv1, ParamCurve* cv2, double aepsge,
		       double tmin1, double tmax1, double tmin2, double tmax2,
		       double seed1, double seed2, int method, bool quick, 
		       double& par1, double& par2, double& dist, 
		       Point& ptc1, Point& ptc2, int& istat);



  /** Compute a closest point or an intersection point between two curves.
   * Ported from the sisl-functions s1770, s1770_s9corr, s1770_s9dir
   * \param cv1   Curve number one.
   * \param cv2   Curve number two.
   * \retval par1 Parameter value of the first curve's point.
   * \retval par2 Parameter value of the second curve's point.
   * \retval dist Distance between the points.
   * \retval ptc1 Point on curve number one.
   * \retval ptc2 Point on curve number two.
   */
  void closestPtCurves(const ParamCurve* cv1, const ParamCurve* cv2,
		       double& par1, double& par2, double& dist,
		       Point& ptc1, Point& ptc2);


  /** Newton iteration on the distance function between two curves to
   * find a closest point or an intersection point.
   * The solution is restricted to the intervals [tmin1,tmax1] and 
   * [tmin2,tmax2] in the first and second curves's parameter domain.
   * Ported from the sisl function s1770.
   * \param cv1 Curve number one.
   * \param cv2 Curve number two.
   * \param tmin1 Min. parameter value. First curve  
   * \param tmax1 Max. parameter value. First curve
   * \param tmin2 Min. parameter value. Second curve
   * \param tmax2 Max. parameter value. Second curve
   * \param seed1 Start point for iteration along curve number one.
   * \param seed2 Start point for iteration along curve number two.
   * \retval par1 Parameter value of the first curve's point.
   * \retval par2 Parameter value of the second curve's point.
   * \retval dist Distance between the points.
   * \retval ptc1 Point on curve number one.
   * \retval ptc2 Point on curve number two.
   */
  void closestPtCurves(const ParamCurve* cv1, const ParamCurve* cv2, double tmin1,
		       double tmax1, double tmin2, double tmax2,
		       double seed1, double seed2, double& par1, double& par2,
		       double& dist, Point& ptc1, Point& ptc2);



  /** Newton iteration on the distance function between
   * a curve and a surface to find a closest point or an intersection point.
   * The solution is restricted to the interval [astart1,aend1] in the
   * curves's parameter domain and [estart2[],eend2[]] in the surface's 
   * parameter domain.
   * Ported from the sisl-function \c s1772.
   * \param pcurve Pointer to the curve in the intersection.
   * \param psurf  Pointer to the surface in the intersection.
   * \param aepsge Geometry resolution.
   * \param astart1 Start parameter value of the curve.
   * \param aend1 End parameter value of the curve.
   * \param domain the parametric domain of interest of the surface
   * \param anext1 Start point for iteration along curve number one.
   * \param enext2[] Start parameter values of the iteration on the surface.
   * \param second_order if 'true' the function will attempt a second-order method 
   *                     directly.
   * \retval cpos1 Parameter value of the curve in the intersection point.
   * \retval gpos2[] Parameter values of the surface in the intersection point.
   * \retval dist Distance between the points.
   * \retval pt_cv The point on the curve.
   * \retval pt_su The point on the surface.
   */
    void closestPtCurveSurf(ParamCurve* pcurve, ParamSurface* psurf, double aepsge,
			    double astart1, double aend1, RectDomain *domain,
			    double anext1, double enext2[],
			    double& cpos1, double gpos2[],
			    double& dist, Point& pt_cv, Point& pt_su,
			    bool second_order=false);

  /** Newton iteration on the distance function between
   * a curve and a surface to find a closest point or an intersection point.
   * The solution is restricted to the interval [astart1,aend1] in the
   * curves's parameter domain and [estart2[],eend2[]] in the surface's 
   * parameter domain.
   * Ported from the sisl-function \c s1772.
   * \param pcurve Pointer to the curve in the intersection.
   * \param psurf  Pointer to the surface in the intersection.
   * \param aepsge Geometry resolution.
   * \param astart1 Start parameter value of the curve.
   * \param estart2[] Start parameter values of surface.
   * \param aend1 End parameter value of the curve.
   * \param eend2[] End parameter values of the surface.
   * \param anext1 Start point for iteration along curve number one.
   * \param enext2[] Start parameter values of the iteration on the surface.
   * \param second_order if 'true' the function will attempt a second-order method 
   *                     directly.
   * \retval cpos1 Parameter value of the curve in the intersection point.
   * \retval gpos2[] Parameter values of the surface in the intersection point.
   * \retval dist Distance between the points.
   * \retval pt_cv The point on the curve.
   * \retval pt_su The point on the surface.
   * \retval istat =1 Intersection. =2 Minimum value.  =3 Nothing found.
   */
void closestPtCurveSurf(ParamCurve* pcurve, ParamSurface* psurf, double aepsge,
		       double astart1, double estart2[], double aend1,
		       double eend2[], double anext1, double enext2[],
		       double& cpos1, double gpos2[],
		       double& dist, Point& pt_cv, Point& pt_su, int& istat,
			bool second_order=false);



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

} // namespace ClosestPoint

} // namespace Go

#endif // _CLOSESTPOINT_
