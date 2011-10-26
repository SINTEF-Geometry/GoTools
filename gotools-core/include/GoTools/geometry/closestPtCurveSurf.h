//===========================================================================
//                                                                           
// File: closestPtCurveSurf.h
//                                                                           
// Created:
//                                                                           
// Author: B. Spjelkavik <bsp@sintef.no>
//          
// Revision:  $Id: closestPtCurveSurf.h,v 1.5 2005-11-17 08:55:21 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _CLOSESTPTCURVESURF_H
#define _CLOSESTPTCURVESURF_H

/** @file closestPtCurveSurf.h
 * Declaration file for an intersection function operating on 
 * object of class ParamCurve and ParamSurface belonging to geometry.
 * For 3D curves only.
 * Ported from the sisl-functions s1772, s1772_s9corr, s1772_s9dir,
 * s1772_s6local_pretop and s1772_s6sekant1.
 * The equation solver Gaussian Elimintion used in s1772 has been
 * replaced by Cramer's rule to solve the 3x3 equation system.
 */


#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"

namespace Go {

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


// Anonymous namespace
namespace {
  /** To be sure we are inside the intervals [astart1,aend1] and
   * [astart2,aend2]. If we are outside clipping is used to adjust
   * the step value.
   * Ported from the sisl function \c s1772_s9corr.
   * \param gd Old and new step value.
   * \param acoef Parameter values
   * \param corr If the step value has been changed, corr is increased
   * by one. If not, corr is set to zero.
   */
void insideParamDomain(Point& gd, const Point& acoef, double astart1, 
		       double aend1, double astart2[],double aend2[],
		       int& corr);

  /** Make and solve an equation system to compute the distance vector and
   * and the length of the distance vector, and to compute a next step on
   *  all three parameter directions towards an intersection or a closest
   *  point. Ported from the sisl-function s1772_s9dir.
    * \param dist The length of the distance vector between the points.
   * \param diff Distance vector between the points.
   * \param delta Relative step parameter values towards intersection on
   * the curve delta[2] and the surface (delta[0],delta[1]).
   * \param kstat kstat=0:Solution found.  kstat=1:Singular system.
   * \param eval_cv Value and derivatives on the curve.
   * \param eval_su Value and derivatives on the surface.
   * \param order Order of the series expansion. One or two.
   */
void nextStep(double& dist,Point& diff,Point& delta,int& kstat,
	      std::vector<Point>& eval_cv,std::vector<Point>& eval_su, 
	      int order);


  /** Finds if we have a minimum or a maximum or a point of inflection
   * situation. This function assume that it is a singular situation.
   * Ported from the sisl function \c s1772_s6local_pretop.
   * \param dist Length of the distance vector between the points.
   * \param diff Distance vector between the points.
   * \param normal The normal vector on the surface.
   * \param eval_cv Value and derivatives in the point on the curve.
   * \param eval_su Value and derivatives in the point on the surface.
   * \return 0:Maximum position or inflection point. 1:Minimum position. 
   */
int localPretop(double dist,const Point& diff,const Point& normal,
		const std::vector<Point>& eval_cv,
		const std::vector<Point>& eval_su);


  /** Singularity (closest point) found. If we have a point of inflection
   * or a maximum position situation, use the secant method. Parameters
   * should stay in the intervals [astart1,aend1] and [estart2,eend2] for
   * curve 1 and curve 2.
   * \param cv1   Curve number one.
   * \param cv2   Curve number two.
   * \param quick Reduce requirement on exactness.
   * \param aepsge Geometry resolution.
   * \param delta Parameter distance on the curve beetveen start values.
   * \param diff Distance vector between the points.
   * \param norm_vec The normal vector on the surface.
   * \param eval_cv Value and derivatives on the first curve.
   * \param eval_su Value and derivatives on the second curve.
   * \retval par_val Parameter values for the curve and the surface.
   * \retval dist Distance between the points.

   */
void singular(ParamCurve* pcurve,ParamSurface* psurf,Point& par_val,
	      double& dist, double aepsge,double delta,const Point& diff,
	      const Point& norm_vec,
	      const std::vector<Point>&eval_cv,
	      const std::vector<Point>&eval_su,
	      double astart1,double estart2[],double aend1,double eend2[]);


  /** Secant method iteration on the distance function between
   * a curve and a surface to find a closest point or an intersection point.
   * Parameters should stay in the intervals [astart1,aend1] and
   * [estart2,eend2] for the curve and the surface.
   * Ported from the sisl function \c s1772_s6sekant1.
   * \param pcurve Pointer to the curve in the intersection.
   * \param psurf Pointer to the surface in the intersection.
   * \param delta Parameter distance on the curve between start values.
   * \param aepsge Geometry resolution.
   * \retval par_val  Parameter values for the surface and the curve.
   * \retval dist Distance in space between the points.
   * \retval jstat=0:OK  =2:Can't find a reasonable start point  ????
   */
void secant(ParamCurve *pcurve,ParamSurface *psurf,Point& par_val,
	    double& dist,int& jstat, double delta,double aepsge,
	    double astart1,double estart2[],double aend1,double eend2[]);


  /** Set the values to be returned from the function \c closestPtCurveSurf
   * \param par_val Parameter values the surface and the curve.
   * \param pcurve  Pointer to the curve.
   * \param psurf   Pointer to the surface.
   * \param aepsge Geometry resolution.
   * \retval jstat =1 Intersection. =2 Minimum value.  =3 ???.
   * \retval par_cv Parameter value of the first curve's point.
   * \retval par_su Parameter value of the second curve's point.
   * \retval dist Distance in space between the points.
   * \retval ptc1 The point on the curve.
   * \retval ptc2 The point on the surface.
   */
void setReturnValues(const Point& par_val,ParamCurve* pcurve,
		     ParamSurface* psurf,double aepsge,
		     int& jstat,double& par_cv,double par_su[],
		     double& dist,Point& pt_cv,Point& pt_su);


  /** Evaluate the curve in a start point, find the closest point on the
   * surface, compute the distance vector between the points, the length
   * of the distance vector, and compute the scalar product of the distance
   * vector and the normal vector on the surface in the closest point.
   * \param pcurve Pointer to the curve.
   * \param par_cv Parameter value of start point on the curve.
   * \param psurf Pointer to the surface.
   * \param estart2 Start of the surface's parameter interval.
   * \param eend2 End of the surface's parameter interval.
   * \param aepsge Geometry resolution.
   * \retval par_su Parameter value of closest point on the surface.
   * \retval y Scalar product of the distance vector and the normal vector
   * on the surface in the closest point.
   * \retval dist Length of the distance vector.
   * \retval jstat 0:OK.  -1:The curve has a singularity.(Should not happen.)
   *  -2:We have found an intersection point.
   */
void newPointEval(ParamCurve *pcurve,double par_cv,ParamSurface *psurf,
		  double estart2[],double eend2[],double aepsge,
		  Point& par_su,double& y,double& dist,int& jstat);


} // Anonymous namespace

} // namespace Go  

#endif // _CLOSESTPTCURVESURF_H
