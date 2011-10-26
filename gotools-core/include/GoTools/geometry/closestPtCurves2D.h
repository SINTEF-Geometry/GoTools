//===========================================================================
//                                                                           
// File: closestPtCurves2D.h
//                                                                           
// Created:
//                                                                           
// Author: B. Spjelkavik <bsp@sintef.no>
//          
// Revision:  $Id: closestPtCurves2D.h,v 1.1 2004-03-11 11:15:11 bsp Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _CLOSESTPTCURVES2D_H
#define _CLOSESTPTCURVES2D_H

/** @file closestPtCurves2D.h
 * Declaration file for an intersection function operating on 
 * objects of class ParamCurve belonging to geometry.
 * For 2D curves only.
 * Ported from the sisl-functions s1770_2D, s1770_2D_s9corr, s1770_2D_s9dir,
 * s1770_2D_s6local_pretop and s1770_2D_s6sekant1.
 * The equation solver Gaussian Elimintion used in s1770_2D has been
 * replaced by Cramer's rule to solve the 2x2 equation system.
 */


#include "GoTools/geometry/ParamCurve.h"

namespace Go {

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


// Anonymous namespace
namespace {
  /** To be sure we are inside the intervals [astart1,aend1] and
   * [astart2,aend2]. If we are outside clipping is used to adjust
   * the step value.
   * Ported from the sisl function \c s1770_2D_s9corr.
   * \param gd Old and new step value.
   * \param acoef Parameter values
   * \param corr If the step value has been changed, corr is increased
   * by one. If not, corr is set to zero.
   */
void insideParamDomain2D(Point& gd, const Point& acoef, double astart1, 
			 double aend1, double astart2,double aend2, int& corr);


  /** Make and solve an equation system to compute the distance vector and
   * value beetween a point on the first curve and a point on the second
   * curve. And to compute a next step on both curves.
   * This is equivalent to the nearest way to the
   * parameter plane in the tangent plane from a point in the
   * distance surface between two curves.
   * \param dist Length of the distance vector between the points.
   * \param diff Distance vector between the points.
   * \param delta Relative step parameter values towards intersection on
   * the curve 1 delta[0] and the curve 2 delta[1].
   * \param det Value of determinant of the equation system
   * \param kstat kstat=0:Solution found.  kstat=1:Singular system.
   * \param eval1 Value and derivatives of the first curve.
   * \param eval2 Value and derivatives of the second curve.
   * \param method method=1 or 2: Use the the method in s1770_2D, starting
   * with order one or two. method=3: Use the method from s1770 (order one).
   */
void nextStep2D(double& dist, Point& diff, Point& delta, double& det,
		int& kstat, std::vector<Point>& eval1, 
		std::vector<Point>& eval2, int method);


  /** Fills the 2x2 matrix A and right hand side b of the equation system Ax=b.
   * This is an equation for computing the next step of the parameters s and t
   * for the 2D curves f and g in an iteration for computing their
   * intersection point or closest point. The equation is the one in the
   * sisl function s1770_2D_s9dir with order equal to one.
   * @verbatim
                      | ds |   |        |
    | -g'(s)  f'(t) | |    | = | d(s,t) |
                      | dt |   |        |  
   @endverbatim
   * \param A Matrix
   * \param b Right hand side
   * \param eval1 Curve g evaluated at s
   * \param eval2 Curve f evaluated at t
   * \param diff Distance vector between the space points g(s) and f(t)
   */
void makeEqSys1(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff);


  /** Fills the 2x2 matrix A and right hand side b of the equation system Ax=b.
   * This is an equation for computing the next step of the parameters s and t
   * for the 2D curves f and g in an iteration for computing their
   * intersection or closest point. The equation is the one in
   * sisl function s1770_2D_s9dir with order equal to two.
   * @verbatim

    | -<g'(s),g'(s)>-<d(s,t),g"(s)  <g'(s),f'(t)> | | ds |   | <d(s,t),g'(s)> |
    |                                             | |    | = |                |
    | -<g'(s),f'(t)>  <f'(s),f'(s)>-<d(s,t),f"(s)>| | dt |   | <d(s,t),f'(t)> |

   @endverbatim
   * \param A Matrix
   * \param b Right hand side
   * \param eval1 Curve g evaluated at s
   * \param eval2 Curve f evaluated at t
   * \param diff Distance vector between the space points g(s) and f(t)
   */
void makeEqSys2(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff);


  /** Fills the 2x2 matrix A and right hand side b of the equation system Ax=b.
   * This is an equation for computing the next step of the parameters s and t
   * for the 2D curves f and g in an iteration for computing their
   * intersection or closest point. The equation is the one in
   * sisl function s1770_s9dir.
   * @verbatim

    | -<g'(s),f'(t)>  <f'(t),f'(t)> | | ds |   | <d(s,t),f'(t)> |
    |                               | |    | = |                |
    | -<g'(s),g'(s)>  <g'(s),f'(t)> | | dt |   | <d(s,t),g'(s)> |

   @endverbatim
   * \param A Matrix
   * \param b Right hand side
   * \param eval1 Curve g evaluated at s
   * \param eval2 Curve f evaluated at t
   * \param diff Distance vector d(s,t) between the space points g(s) and f(t)
   */
void makeEqSys3(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff);


  /** Finds if we have a minimum or a maximum or a point of inflection
   * situation. This function assume that it is a singular situation.
   * Ported from the sisl function \c s1770_2D_s6local_pretop.
   * \param dist The length of the distance vector between the points.
   * \param diff Distance vector between the points.
   * \param eval1 Value and derivatives on the first curve.
   * \param eval2 Value and derivatives on the second curve.
   * \return 0:Maximum position or inflection point. 1:Minimum position. 
   */
int localPretop(double dist,const Point& diff,
		const std::vector<Point>& eval1,
		const std::vector<Point>& eval2);


  /** Singularity (closest point) found. If we have a point of inflection
   * or a maximum position situation, use the secant method. Parameters
   * should stay in the intervals [astart1,aend1] and [astart2,aend2] for
   * curve 1 and curve 2.
   * \param cv1   Curve number one.
   * \param cv2   Curve number two.
   * \param quick Reduce requirement on exactness.
   * \param aepsge Geometry resolution.
   * \param delta Parameter distance on first curve beetveen start values.
   * \param diff Distance vector between the points.
   * \param eval1 Value and derivatives of the first curve.
   * \param eval2 Value and derivatives of the second curve.
   * \retval par_val Parameter values for the two curves.
   * \retval dist Distance between the points.

   */
void singular(ParamCurve* cv1, ParamCurve* cv2,Point& par_val,double& dist,
	      int quick,double aepsge,double delta,const Point& diff,
	      const std::vector<Point>&eval1,const std::vector<Point>&eval2,
	      double astart1,double astart2,double aend1,double aend2);


  /** Secant method iteration on the distance function between
   * two curves to find a closest point or an intersection point.
   * Parameters should stay in the intervals [astart1,aend1] and
   * [astart2,aend2] for curve 1 and curve 2.
   * Ported from the sisl function \c s1770_2D_s6sekant1.
   * \param pcurve1 Curve number one.
   * \param pcurve2 Curve number two.
   * \param delta Parameter distance on first curve betveen start values.
   * \param aepsge Geometry resolution.
   * \retval par_val Parameter values of the curves points.
   * \retval dist Distance in space between the points.
   * \retval jstat=0:OK  =2:Can't find a reasonable start point  
   */
void secant2D(ParamCurve *pcurve1,ParamCurve *pcurve2,Point& par_val,
	      double& dist,int& jstat, double delta,double aepsge,
	      double astart1,double astart2,double aend1,double aend2);


  /** Set the values to be returned from the function \c closestPtCurves2D.
   * \param par_val Parameter values for the two curves.
   * \param pcurve1   Curve number one.
   * \param pcurve2   Curve number two.
   * \param aepsge Geometry resolution.
   * \retval jstat =1 Intersection. =2 Minimum value.  =3 ???.
   * \retval par1 Parameter value of the first curve's point.
   * \retval par2 Parameter value of the second curve's point.
   * \retval dist Distance in space between the points.
   * \retval ptc1 Point on curve number one.
   * \retval ptc2 Point on curve number two.
   */
void setReturnValues(const Point& par_val,ParamCurve* cv1,ParamCurve* cv2,
		     double aepsge,  int& jstat,double& par1,double& par2,
		     double& dist,Point& ptc1,Point& ptc2);


  /** Evaluate curve 1 in a start point, find the closest point on curve 2,
   * and compute the distance vector between the points, the length
   * of the distance vector and the scalar product of the distance vector
   * and the normal vector on curve 2 in the closest point.
   * \param pcurve1 Curve number one.
   * \param par1 Parameter value of start point on curve 1.
   * \param astart2 Start of curve 2's parameter interval.
   * \param aend2 End of curve 2's parameter interval. 
   * \param pcurve2   Curve number two.
   * \param aepsge Geometry resolution
   * \retval par2 Parameter value of closest point on curve 2.
   * \retval y Scalar product of the distance vector and the normal vector
   * on curve 2 in the closest point.
   * \retval dist Length of the distance vector.
   * \retval jstat 0:OK.  -1:The curve has a singularity.(Should not happen.)
   *  -2:We have found an intersection point.
   */
void newPointEval(ParamCurve *pcurve1,double par1,ParamCurve *pcurve2,
		  double astart2,double aend2,double aepsge,
		  double& par2,double& y,double& dist,int& jstat);

} // Anonymous namespace
 
} // namespace Go

#endif // _CLOSESTPTCURVES2D_H
