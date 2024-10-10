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

#include <limits>
#include <algorithm>
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/utils/CurvatureUtils.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/intersections/Param1FunctionInt.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/geometry/BoundedSurface.h"

using std::numeric_limits;
using std::vector;
using std::swap;
using std::copy;

using namespace Go;

namespace {
const double LARGEST_POSSIBLE_VALUE = numeric_limits<double>::max();

shared_ptr<ParamCurve>
make_surf_curve_preserve_critical(const ParamObjectInt* pobj,
				  int running_dir,
				  double fixed_value);

// shared_ptr<ParamCurve>
// make_surf_curve_preserve_critical(shared_ptr<const SplineSurface> surf,
// 				  bool critical_values,
// 				  int running_dir,
// 				  double fixed_value);

bool measure_coincidence_region(const ParamCurve* marching_curve,
			      const ParamGeomInt* second_object,
			      double isect_param_marching,
			      const double* const isect_param_other,
			      bool step_forward_marching,
			      const shared_ptr<const GeoTol> gtol,
			      double& last_pval_inside_marching,
			      double& first_pval_outside_marching);

bool measure_coincidence_region(const ParamCurve* marching_curve_1d,
				double C,
				double isect_param_marching,
				const double* const isect_param_other,
				bool step_forward_marching,
				const shared_ptr<const GeoTol> gtol,
				double& last_pval_inside_marching,
				double& first_pval_outside_marching);


//===========================================================================
// Abstract base class for two helper classes used by the
// determineCoincidenceRegion() function.
// It is used to step along and analyse the coincidence region between two objects.  
// One of the objects is always a curve (which we are marching along), and the other 
// is either a curve or a surface (as specified in the derived classes BindToCurve and
// BindToSurface).  It is supposed that the two objects intersect.
class BindToObject
//===========================================================================
{
public:
    BindToObject(bool step_forward, 
		 const shared_ptr<const GeoTol> tol,
		 const ParamCurve* marching_cv)
	: step_forward_(step_forward), tol_(tol), march_cv_(marching_cv),
	  fuzzy_march_(fuzzy_eps_*(march_cv_->endparam()-march_cv_->startparam())) {};

    virtual ~BindToObject() {}

    // Expand the coincidence region of the two bound objects until we have stepped 
    // outside the region (if possible) or until we have reached the end of the parameter
    // domain.
    // first_pval_outside_curve - On input, this argument should be set to the parameter
    //                            value for the intersection point in the curve.  It 
    //                            will return a parameter value for the curve that is
    //                            outside the coincidence region, or on the end of the
    //                            curve's parameter domain.
    // first_pval_outside_other - On input, this argument should be set to the parameter 
    //                            value(s) for the intersection point in the second object.
    //                            It will return the parameter value(s) for the other object
    //                            that is outside the coincidence region, or at the end of
    //                            the object's parameter domain.  (Implemented as a pointer
    //                            rather than a reference, since the number of parameters
    //                            of the second object is not known by this base class.  The
    //                            user is responsable for preparing the memory area that will
    //                            contain the final values).
    // current_distance         - Returns the distance of the two points found on the two
    //                            objects outside the coincidence region (parameterized by
    //                            'first_pval_outside_curve' and 'first_pval_outside_other').
    // exited_Region            - If we were able to step outside the coincidence region 
    //                            before hitting the bounds of the parameter domains of the
    //                            objects, this argument will be satt to 'true'.  Otherwise
    //                            it will be 'false'.
    virtual void expandCoincidenceRegion(double& first_pval_outside_curve,
					 double* first_pval_outside_other,
					 double& current_distance,
					 bool& exited_region) = 0;

    // After the bound of the confidence interval has been detected and bracketed, this 
    // function can be called to narrow the bracket down to a certain tolerance, inside
    // which we know that the bound will be located.  The maximum allowed bracket size is
    // specified inside the inherited GeoTol object.  This function is typically called
    // after expandCoincidenceRegion().
    // last_pval_inside_curve   - lower bracket for the parameter interval containing the
    //                            bound of the influence area in the curve object.  If 
    //                            necessary, it will be changed in order to tighten the
    //                            bracket interval.
    // first_pval_outside_curve - upper bracket for the paramter interval containing the
    //                            bound of the influence area in the curve object.  If
    //                            necessary, it will be changed in order to tighten the
    //                            bracket interval.
    // last_pval_inside_other   - points to the value(s) containing the lower bracket(s)
    //                            for the parameter interval(s) containing the bound(s) of 
    //                            the influence area in the other object.  If necessary,
    //                            the value(s) will be changed in order to tighten the
    //                            bracket interval.
    // first_pval_outside_other - points to the value(s) containing the upper bracket(s)
    //                            for the parameter interval(s) containing the bound(s) of
    //                            the influence area in the other object.  If necessary,
    //                            the value(s) will be changed in order to tighten the
    //                            bracket interval.
    // inside_dist              - returns the distance between the two objects at the 
    //                            parameter values 'last_pval_inside_*****'. (should be less
    //                            than the geometric tolerance.
    // outside_dist             - returns the distance between the two objects at the 
    //                            parameter values 'first_pval_outside_*****'. (should be
    //                            larger than the geometric tolerance.
    void determineConfidenceInterval(double& last_pval_inside_curve,
				     double* last_pval_inside_other,
				     double& first_pval_outside_curve,
				     double* first_pval_outside_other,
				     double& inside_dist,
				     double& outside_dist);
    // should return the number of parameters in the second object (1 or 2)
    virtual int numParametersOtherObject() const = 0;
    
protected:
    // protected member data
    static const double fuzzy_eps_;
    bool step_forward_;
    const shared_ptr<const GeoTol> tol_;

    const ParamCurve* march_cv_;
    const double fuzzy_march_;

    // protected member functions
    virtual void other_closest_point(const Point& pt,
				     double* cl_par,
				     double& cl_dist,
				     Point& cl_pt,
				     bool mind_segment_values) = 0;
}; 
    
const double BindToObject::fuzzy_eps_ = 1.0e-4; //@ move elsewhere??

//===========================================================================
// helper class for the determineCoincidenceRegion function.
// Only applicable for 1-dimensional first_obj (i.e. march_cv).
class BindToPoint : public BindToObject
//===========================================================================
{
public:
    BindToPoint(const ParamCurve* marching_cv, 
		Point other_point,
		bool step_forward,
		const shared_ptr<const GeoTol> tol) 
	: BindToObject(step_forward, tol, marching_cv), other_point_(other_point)
	 {}

    virtual ~BindToPoint() {}

    virtual void expandCoincidenceRegion(double& first_pval_outside_curve,
					 double* first_pval_outside_other,
					 double& current_distance,
					 bool& exited_region);

//     virtual void determineConfidenceInterval(double& last_pval_inside_curve,
// 					     double* last_pval_inside_other,
// 					     double& first_pval_outside_curve,
// 					     double* first_pval_outside_other,
// 					     double& inside_dist,
// 					     double& outside_dist);
    virtual int numParametersOtherObject() const {return 0;}

protected:
    virtual void other_closest_point(const Point& pt,
				     double* cl_par,
				     double& cl_dist,
				     Point& cl_pt,
				     bool mind_segment_values);
private:
    const Point other_point_; // 1-dim
};

//===========================================================================
// helper class for the determineCoincidenceRegion function.
class BindToCurve : public BindToObject
//===========================================================================
{
public:
    BindToCurve(const ParamCurve* marching_cv, 
		const ParamCurve* other_cv,
		bool step_forward,
		const shared_ptr<const GeoTol> tol) 
	: BindToObject(step_forward, tol, marching_cv), other_cv_(other_cv),
	  fuzzy_other_(fuzzy_eps_ * (other_cv_->endparam() - other_cv_->startparam())) {}

    virtual ~BindToCurve() {}

    virtual void expandCoincidenceRegion(double& first_pval_outside_curve,
					 double* first_pval_outside_other,
					 double& current_distance,
					 bool& exited_region);

//     virtual void determineConfidenceInterval(double& last_pval_inside_curve,
// 					     double* last_pval_inside_other,
// 					     double& first_pval_outside_curve,
// 					     double* first_pval_outside_other,
// 					     double& inside_dist,
// 					     double& outside_dist);
    virtual int numParametersOtherObject() const {return 1;}

protected:
    virtual void other_closest_point(const Point& pt,
				     double* cl_par,
				     double& cl_dist,
				     Point& cl_pt,
				     bool mind_segment_values);
private:
    const ParamCurve* other_cv_;
    const double fuzzy_other_;
};


// Helper class for the determineCoincidenceRegion function.

//===========================================================================
class BindToSurface : public BindToObject {
//===========================================================================
public:
    BindToSurface(const ParamCurve* marching_cv, 
		  const ParamSurfaceInt* other_sf,
		  bool step_forward,
		  const shared_ptr<const GeoTol> tol) 
	: BindToObject(step_forward, tol, marching_cv), other_sf_(other_sf)
    {
	const SplineSurface* spline_sf = NULL;
	const ParamSurface* par_sf = other_sf_->getParamSurface().get();
	if (par_sf->instanceType() == Class_SplineSurface) {
	    spline_sf = dynamic_cast<const SplineSurface*>(par_sf);
	} else if (par_sf->instanceType() == Class_BoundedSurface) {
	    const BoundedSurface* bd_sf
		= dynamic_cast<const BoundedSurface*>(par_sf);
	    spline_sf =
		dynamic_cast<const SplineSurface*>
		(bd_sf->underlyingSurface().get());
	}
	fuzzy_u_ = fuzzy_eps_ * (spline_sf->endparam_u()
				 - spline_sf->startparam_u());
	fuzzy_v_ = fuzzy_eps_ * (spline_sf->endparam_v()
				 - spline_sf->startparam_v());

    }

    virtual ~BindToSurface() {}

    virtual void expandCoincidenceRegion(double& first_pval_outside_curve,
					 double* first_pval_outside_other,
					 double& current_distance,
					 bool& exited_region);

//     virtual void determineConfidenceInterval(double& last_pval_inside_curve,
// 					     double* last_pval_inside_other,
// 					     double& first_pval_outside_curve,
// 					     double* first_pval_outside_other,
// 					     double& inside_dist,
// 					     double& outside_dist);

    virtual int numParametersOtherObject() const
    { return 2; }

protected:
    virtual void other_closest_point(const Point& pt,
				     double* cl_par,
				     double& cl_dist,
				     Point& cl_pt,
				     bool mind_segment_values);

private:
    const ParamSurfaceInt* other_sf_;
    double fuzzy_u_;
    double fuzzy_v_;

};


// //===========================================================================
// // helper class for the determineCoincidenceRegion function.
// class BindToPoint : public BindToObject
// //===========================================================================
// {
// public:
//     BindToPoint(const ParamCurve* marching_cv, 
// 		Point other_pt,
// 		bool step_forward,
// 		const shared_ptr<const GeoTol> tol) 
// 	: BindToObject(step_forward, tol, marching_cv), other_cv_(other_cv),
// 	  fuzzy_other_(fuzzy_eps_ * (other_cv_->endparam() - other_cv_->startparam())) {}

//     virtual ~BindToPoint() {}

//     virtual void expandCoincidenceRegion(double& first_pval_outside_curve,
// 					 double* first_pval_outside_other,
// 					 double& current_distance,
// 					 bool& exited_region);

//     virtual void determineConfidenceInterval(double& last_pval_inside_curve,
// 					     double* last_pval_inside_other,
// 					     double& first_pval_outside_curve,
// 					     double* first_pval_outside_other,
// 					     double& inside_dist,
// 					     double& outside_dist);
//     virtual int numParametersOtherObject() const {return 0;}

// protected:
//     // Returning other_pt_, cl_par not set.
//     virtual void other_closest_point(const Point& pt,
// 				     double* cl_par,
// 				     double& cl_dist,
// 				     Point& cl_pt,
// 				     bool mind_segment_values);
// private:
//     const Point other_pt_;
//     const double fuzzy_other_;
// };


}; // end anonymous namespace 

namespace Go
{

//===========================================================================
   // Returns true if the parameter par has passed next_val.
static bool hasPassedVal(double par, double next_val, bool step_forward)
//===========================================================================
{
  bool passed;

  if (step_forward)
    passed = (par > next_val)? true : false;
  else
    passed = (par < next_val)? true : false;

  return passed;
}

//===========================================================================
    // Returns true if the parameter par has reached par_end
static bool hasReachedEnd(double par, double par_end, bool step_forward)
//===========================================================================
{
  bool ended;

  if (step_forward)
    ended = (par >= par_end)? true : false;
  else
    ended = (par <= par_end)? true : false;

  return ended;
}

//===========================================================================
    // Returns the the value of the first of the two parameters par1 and par2. 
static double nextLimit(double par1, double par2, bool step_forward)
//===========================================================================
{
  if (step_forward)
    return std::min(par1, par2);
  else
    return std::max(par1, par2);
}


//===========================================================================
    // Sets stepping direction and min and max values of the search interval.
static void setRangeAndDirection(double tstart, double tend, 
				 bool& step_forward, 
				 double& tmin, double& tmax)
//===========================================================================
{
  if (tstart < tend) {
    step_forward = true;
    tmin = tstart;
    tmax = tend;
  }
  else {
    step_forward = false;
    tmin = tend;
    tmax = tstart;
  }
}

//===========================================================================
int checkCoincide(ParamCurveInt *curve,
		  double start,double end, shared_ptr<GeoTol> tol,
		  ParamCurveInt *other,
		  double other_start, double other_end)
//===========================================================================
{
    DEBUG_ERROR_IF(curve->dimension() != other->dimension(),
	     "Dimension mismatch.");

    const double& toler = tol->getEpsge();  // @@ ?????
    const double tminstep = 1.e-10;  // @@ ?????
    const double fuzzy_eps = 1.0e4*tol->getRelParRes(); //1.e-4;  // Knot distance tolerance.  @bsp
    const static int nder=2;   // Order of derivatives to be calulated.
    vector<Point> ft;     // Partial derivatives up to order nder
    vector<Point> gs;     // Partial derivatives up to order nder
    ft.resize(nder+1);
    gs.resize(nder+1);


    bool step_forward1;   // True if we are going forward along this curve
    bool step_forward2;   // True if we are going forward along the other curve.
    bool switch_curves;   // True if we are going to switch curves.
    int march_curve;      // Indicates which curve that is marched along. 1 or 2.
    int ntrials;          // Number of trials on one step.  Equals 0, 1 or 2.
    int istat;            // Status variable.
    int niter;             // Number of iterations;
    double pmin1, pmax1;  // Searching range, curve1.
    double pmin2, pmax2;  // Searching range, curve2.
    double tstep;         // Step length.
    double tx, tx1, tx2;  // Parameter values of first curve.
    double ty, ty1, ty2;  // Parameter values of second curve.
    double next_knot1;    // Next knot value, first curve.   
    double next_knot2;    // Next knot value, second curve.
    double next_limit;    // Next knot value or end parameter.
    double pnt_dist;      // Distance between closest points.
    double fuzzy1,fuzzy2; // Knot distance tolerance.
  
    Point pnt1(curve->dimension());  // Point on curve1
    Point pnt2(other->dimension());  // Point on curve2
    Point pnt3(curve->dimension());  // Point on curve1
    Point pnt4(other->dimension());  // Point on curve2
    Point pnt5(curve->dimension());  // Point on curve1
    Point pnt6(other->dimension());  // Point on curve2

    // Tolerances for testing if a parameter is equal to a knot value
    fuzzy1 = fuzzy_eps*(curve->endparam() - curve->startparam());
    fuzzy2 = fuzzy_eps*(other->endparam() - other->startparam());

    if (start > end) {
	std::swap(start, end);
	std::swap(other_start, other_end);
    }

    // Make sure that we do not start or end outside the curves.
  
    curve->assureInRange(start); 
    curve->assureInRange(end); 

    other->assureInRange(other_start); 
    other->assureInRange(other_end); 

    // Check end points
    double mid = 0.5*(start + end);
    double other_mid = 0.5*(other_start + other_end);
    curve->point(pnt1, &start);
    other->point(pnt2, &other_start);
    curve->point(pnt3, &end);
    other->point(pnt4, &other_end);
    curve->point(pnt5, &mid);
    other->point(pnt6, &other_mid);
    if (pnt1.dist(pnt2) > toler || pnt3.dist(pnt4) > toler)
	return 0;   // No coincidence
    if (pnt1.dist(pnt3) <= toler && 
	!(pnt1.dist(pnt5) > toler || pnt3.dist(pnt5) > toler ||
	  pnt2.dist(pnt6) > toler || pnt4.dist(pnt6) > toler ))
	return 1;   // The interval is too short
    if (end - start < tol->getRelParRes())
	return 0;  // Probably a closed curve
    if (fabs(other_end - other_start) < tol->getRelParRes())
	return 0;  // Probably a closed curve
  
    // With both end points equal we then check if curves are both linear. Needed
    // to handle cases with practically unbounded domains (cylinders, planes, cones).
    bool curve_linear = (curve->getParamCurve()->instanceType() == Class_Line);
    bool other_linear = (other->getParamCurve()->instanceType() == Class_Line);
    if (curve_linear && other_linear) {
        return 1;
    }

//    double delta = 0.9*eps*eps;
    // double init_dist = std::max(pnt1.dist(pnt2), pnt3.dist(pnt4));
    //double delta = (init_dist < 0.9*toler) ? 0.9*toler : toler;
    double delta = (1 + tol->getEpsBracket()) * tol->getEpsge();

    // Set searching range and directions.
    setRangeAndDirection(start, end, step_forward1, pmin1, pmax1); 
    setRangeAndDirection(other_start, other_end, step_forward2, pmin2, pmax2); 


    // Inital step along first curve
    march_curve = 1;
    tx1 = start;
    ty1 = other_start;
    ntrials = 0;
    niter = 0;

    while (!hasReachedEnd(tx1, end, step_forward1) &&
	   !hasReachedEnd(ty1, other_end, step_forward2)) {

	ntrials++;
	niter++;
	switch_curves = false;
	next_knot1 = curve->nextSegmentVal(tx1,step_forward1, fuzzy1);
	next_knot2 = other->nextSegmentVal(ty1,step_forward2, fuzzy2);

// 	cout << "\nI " << niter <<" "<< ntrials <<" "<< march_curve <<" "
// 	     <<" Cv1 "<< tx1 <<" "<< next_knot1
// 	     <<" Cv2 "<< ty1 <<" "<< next_knot2 << endl;

	if (march_curve == 1) {

	    // Evaluate 0-2nd derivatives of both curves.
	    curve->point(ft, &tx1, nder);
	    other->point(gs, &ty1, nder);

	    istat = stepLength(ft, gs, delta, step_forward1, tstep);
	    ALWAYS_ERROR_IF(istat<1,"Step along curve1 not found");
	    if (fabs(tstep)<tminstep) {
//		cout << "\n Small step " << tstep << endl;
		tstep = step_forward1 ? tminstep : -tminstep;
	    }

	    // Make sure that we don't pass any knots or end-parameter
	    next_limit = nextLimit(next_knot1, end, step_forward1);
	    if (hasPassedVal(tx1+tstep, next_limit, step_forward1))
		tstep = next_limit-tx1;

	    tx2 = tx1 + tstep;

	    tx=0.5*(tx1+tx2);
	    for (int i=1; i<=2; ++i,tx=tx2) {
		curve->point(pnt1, &tx);

		// Avoid going back
		double tmin = (step_forward2) ? ty1 : pmin2;
		double tmax = (step_forward2) ? pmax2 : ty1;
		other->getParamCurve()->closestPoint(pnt1,tmin,tmax,ty,
						     pnt2,pnt_dist,&ty1);
		//curve->knotIntervalFuzzy(ty, fuzzy2);
		//cout << "M1 "<< i <<" " << tx <<" "<< ty <<" "<< pnt_dist << endl;

		// Make sure that we don't pass any knots or end-parameter
		next_limit = nextLimit(next_knot2, other_end, step_forward2);
		if (ntrials < 2 && hasPassedVal(ty, next_limit, step_forward2)) {
		    march_curve = 2;  // switch curves
		    switch_curves = true;
		    break;
		}
	
		if (pnt_dist > toler)
		    return 0;  // Curves not within tolerance;

	    }
	    if (switch_curves)
		continue;

	    tx1 = tx2;
	    ty1 = ty;

	    ntrials = 0;

	}
 
	else {  // march_curve == 2
	    // Evaluate 0-2nd derivatives of both curves. 
	    curve->point(ft, &tx1, nder);
	    other->point(gs, &ty1, nder);

	    istat = stepLength(gs, ft, delta, step_forward2, tstep);
	    ALWAYS_ERROR_IF(istat < 1,"Step along curve2 not found");
	    if (fabs(tstep)<tminstep) {
//		cout << "\n Small step " << tstep << endl;
		tstep = step_forward2 ? tminstep : -tminstep;
	    }
	
	    // Make sure that we don't pass any knots or end-parameter
	    next_limit = nextLimit(next_knot2, other_end, step_forward2);
	    if (hasPassedVal(ty1+tstep, next_limit, step_forward2))
		tstep = next_limit-ty1;

	    ty2 = ty1 + tstep;

	    ty=0.5*(ty1+ty2);
	    for (int i=1; i<=2; ++i,ty=ty2) {
		other->point(pnt2, &ty);

		// Avoid going back
		//double tmin = (step_forward1) ? tx1 : pmin1;
		//double tmax = (step_forward1) ? pmax1 : tx1;
		curve->getParamCurve()->closestPoint(pnt2,pmin1,pmax1,tx,
						     pnt1,pnt_dist,&tx1);
		curve->knotIntervalFuzzy(tx, fuzzy1);
//		cout << "M2 "<< i <<" " << tx <<" "<< ty <<" "<< pnt_dist << endl;
	
		// Make sure that we don't pass any knots or end-parameter
		next_limit = nextLimit(next_knot1, end, step_forward1);
		if (ntrials < 2 && hasPassedVal(tx, next_limit, step_forward1)) {
		    march_curve = 1;  // switch curves
		    switch_curves = true;
		    break;
		}
	
		if (pnt_dist > toler)
		    return 0;  // Curves not within tolerance;

	    }
	    if (switch_curves)
		continue;
      
	    tx1 = tx;
	    ty1 = ty2;

	    ntrials = 0;

	} // march_curve == 2 

    } // while (!hasReachedEnd

    // Curves are within tolerance
    return 1;

}

//===========================================================================
int checkCoincide(ParamCurveInt *curve,
		  double start, double end,
		  ParamSurfaceInt* intsurf,
		  Point su_start,
		  Point su_end,
		  shared_ptr<GeoTol> tol)
//===========================================================================
{
//  const double REL_PAR_RES=1.e-10;       // @bsp
    const double REL_PAR_RES = tol->getRelParRes();
  const static int nder=2;   // Order of derivatives to be calulated.
  const double& aepsge = tol->getEpsge();
  std::vector<Point> ft;     // Partial derivatives up to order nder
  std::vector<Point> gs;     // Partial derivatives up to order nder
  ft.resize(nder+1);
  gs.resize((nder+1)*(nder+2)/2);

//  cout << " REL_PAR_RES=" << REL_PAR_RES << endl;

  shared_ptr<ParamSurface> par_surf = intsurf->getParamSurface();

  DEBUG_ERROR_IF(curve->dimension() != par_surf->dimension(),
	   "Dimension mismatch.");

  bool step_forward;    // True if we are going forward along this curve
  int niter;            // Number of iterations;
  //  int kstat;            // Local status variable;
  double pmin, pmax;    // Searching range, this curve.
  double tstep;         // Step length.
  double tx, tx1, tx2;  // Parameter values of this curve.
  double paru, parv;    // Parameter values of the surface.
  double tdist;         // Distance between two points.
  double guess[2];      // Start parameters for closest point iteration.
  double min_step;      // Minimum step length along the curve.
  double max_step;      // Maximum step length along the curve.
  const RectDomain& domain = intsurf->getDegDomain(tol->getEpsge());

  Point clo_pt(curve->dimension());   // Closest point.
  std::vector<Point> cvder(nder+1); 

  if (start > end)
  {
      std::swap(start,end);
      std::swap(su_start, su_end);
  }

  // Make sure that we do not start or end outside the curve
  curve->assureInRange(start); 
  curve->assureInRange(end); 
    // Check end points
  double mid = 0.5*(start + end);
  Point su_mid = 0.5*(su_start + su_end);
  Point pnt1(curve->dimension()); 
  Point pnt2(curve->dimension()); 
  Point pnt3(curve->dimension()); 
  Point pnt4(curve->dimension()); 
  curve->point(pnt1, &start);
  curve->point(pnt2, &end);
  curve->point(pnt3, &mid);
  intsurf->point(pnt4, su_mid.begin());
  if (pnt1.dist(pnt2) <= aepsge && !(pnt1.dist(pnt3) > aepsge ||
				     pnt2.dist(pnt3) > aepsge ||
				     pnt1.dist(pnt4) > aepsge || 
				     pnt2.dist(pnt4) > aepsge))
      return 1;   // The interval is too short

  if (end-start <= REL_PAR_RES && su_start.dist(su_end) > REL_PAR_RES)
      return 0;  // Probably closed curve or surface
 
  if (end-start > REL_PAR_RES && su_start.dist(su_end) <= REL_PAR_RES)
      return 0;  // Probably closed curve or surface
 
  // Set searching range and directions.
  setRangeAndDirection(start, end, step_forward, pmin, pmax);

  // Set min and max step length along the curve.
  //@bsp Hvordan skal de settes? Hvis de bare er avhengige av kurven
  // kan de vaere medlem av klassen og settes i constructor.
  min_step = 0.01*(pmax-pmin);  //@bsp ??
  max_step = 0.1*(pmax-pmin);

  // Inital step along the curve
  tx1 = tx2 = start;
  paru = su_start[0];
  parv = su_start[1];

  // Evaluate 0-2nd derivatives of this curve.
  curve->point(ft, &tx1, nder);

  // Evaluate 0-2nd derivatives of the surface.
  par_surf->point(gs, paru, parv,nder);
 

  niter = 0;
  while (!hasReachedEnd(tx1, end, step_forward)) {
    niter++;
    tx1 = tx2;
    
    tstep = stepLength(ft, gs, step_forward, cvder, 
		       min_step, max_step, aepsge);
    //kstat = stepLength(forward, tstep, cvder);  @ ???
    //ALWAYS_ERROR_IF(kstat<0,"Zero determinant in evalDistCurve",ComputationError());

    // Make sure that we don't pass end-parameter
    if (hasPassedVal(tx1+tstep, end, step_forward))
      tstep = end-tx1;

    tx2 = tx1 + tstep;

    tx=0.5*(tx1+tx2);
    for (int i=1; i<=2; ++i,tx=tx2) {
      curve->point(ft, &tx, nder);
       
      //  Get closest point in the surface.
      guess[0] = paru;
      guess[1] = parv;
      par_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			     REL_PAR_RES, &domain, guess);
      //	         REL_COMP_RES, NULL, guess);    @bsp
      
      
      // Evaluate 0-2nd derivatives of the surface.
      par_surf->point(gs, paru, parv,nder);
 
     
      // Check if point on curve and surface are within positional and
      // angular tolerances.
	    
      tdist = ft[0].dist(gs[0]);

      if (tdist>aepsge) {
	// Points not within tolerances, curve and surface do not coincide.
	return 0;
      }
      
    }
    
  }
  return 1;

  // Curve and surface are within tolerance
  // Check if the marching reached the expected end point
  double ptol = 1.0e-4; //100.0*tol->getRelParRes();
  ptol = std::min(ptol, 
		  0.01*std::max(fabs(su_end[0]-su_start[0]),fabs(su_end[1]-su_start[1])));
  if (fabs(su_end[0]-paru) > ptol || fabs(su_end[1]-parv) > ptol)
      return 0;
  else
      return 1;

}


//===========================================================================
// Denne versjonen av funksjonen har en loop over de to parameter retningene
// for aa sjekke om vi har passert skjoeter.
//  for (int direction=1; direction<=2; ++direction) {
//   C.log
int checkCoincide(ParamCurveInt *curve,
		  double start, double end,
		  SplineSurfaceInt *intsurf,
		  const Point& su_start,
		  const Point& su_end,
		  shared_ptr<GeoTol> tol)
//===========================================================================
{
  const double REL_PAR_RES=1.e-10;      // @bsp
  const double& aepsge = tol->getEpsge();
  const double fuzzy_eps = 1.0e4*tol->getRelParRes(); //1.e-4;  // Knot distance tolerance.   @bsp
  const static int nder=2;   // Order of derivatives to be calulated.
  std::vector<Point> ft;     // Partial derivatives up to order nder
  std::vector<Point> gs;     // Partial derivatives up to order nder
  ft.resize(nder+1);
  gs.resize((nder+1)*(nder+2)/2);

//  std::cout << " REL_PAR_RES=" << REL_PAR_RES << std::endl;

  shared_ptr<const SplineSurface> spline_surf = intsurf->splineSurface();
  DEBUG_ERROR_IF(curve->dimension() != spline_surf->dimension(),
	   "Dimension mismatch.");

  bool step_forward;    // True if we are going forward along this curve.
  int niter;            // Number of iterations.
  //  int kstat;            // Local status variable.
  int uinterval;        // Knot interval in 1. parameter direction of surface.
  int vinterval;        // Knot interval in 2. parameter direction of surface.
  int prev_uinterval;   // Previous knot interval in 1. parameter direction.
  int prev_vinterval;   // Previous knot interval in 2. parameter direction.
  double pmin, pmax;    // Searching range, this curve.
  double tstep;         // Step length.
  double tx, tx1, tx2;  // Parameter values of this curve.
  double next_cvknot;   // Next knot value, this curve.
  double next_limit;    // Next knot value or end parameter, this curve.
  double paru, parv;    // Parameter values of the surface.
  double tdist;         // Distance between two points.
  double guess[2];      // Start parameters for closest point iteration.
  double min_step;      // Minimum step length along the curve.
  double max_step;      // Maximum step length along the curve.

  double par1, par2;    // Parameter values to closest point between two curves
  double seed1, seed2;  // Start parameters for closest point iteration.
  double fuzzy1,fuzzy2; // Knot distance tolerance.
  Point ptc1, ptc2;     //  Closest points on two curves.
  const RectDomain& domain = intsurf->getDegDomain(tol->getEpsge());

  Point clo_pt(curve->dimension());   // Closest point.
  std::vector<Point> cvder(nder+1); //  Scratch vector for derivatives.

  // Tolerances for testing if a parameter is equal to a knot value
  fuzzy1 = fuzzy_eps*( spline_surf->basis_u().endparam() - 
		       spline_surf->basis_u().startparam());
  fuzzy2 = fuzzy_eps*( spline_surf->basis_v().endparam() - 
		       spline_surf->basis_v().startparam());

 // Make sure that we do not start or end outside the curve
  curve->assureInRange(start); 
  curve->assureInRange(end); 
 
    // Check end points
  Point pnt1(curve->dimension()); 
  Point pnt2(curve->dimension()); 
  curve->point(pnt1, &start);
  curve->point(pnt2, &end);
  if (pnt1.dist(pnt2) <= aepsge)
	return 1;   // The interval is too short

  // Set searching range and directions.
  setRangeAndDirection(start, end, step_forward, pmin, pmax);

  // Set min and max step length along the curve.
  //@bsp Hvordan skal de settes? Hvis de bare er avhengige av kurven
  // kan de vaere medlem av klassen og settes i constructor.
  min_step = 0.01*(pmax-pmin);  //@bsp ??
  max_step = 0.1*(pmax-pmin);

 
  // Inital step along the curve
  tx1 = tx2 = start;
  paru = su_start[0];
  parv = su_start[1];

  // Evaluate 0-2nd derivatives of this curve.
  curve->point(ft, &tx1, nder);

  // Evaluate 0-2nd derivatives of the surface.
  spline_surf->point(gs, paru, parv,nder);
 
  // Indices to knot interval of the surface.
  uinterval = spline_surf->basis_u().knotInterval(paru);
  vinterval = spline_surf->basis_v().knotInterval(parv);

  niter = 0;
  while (!hasReachedEnd(tx2, end, step_forward)) {
    niter++;
    tx1 = tx2;
    next_cvknot = curve->nextSegmentVal(tx1, step_forward, fuzzy1);
    prev_uinterval = uinterval;
    prev_vinterval = vinterval;
     
    tstep = stepLength(ft, gs, step_forward, cvder, 
		       min_step, max_step, aepsge);
    //kstat = stepLength(forward, tstep, cvder);  @ ???
 //ALWAYS_ERROR_IF(kstat<0,"Zero determinant in evalDistCurve",ComputationError());

    // Make sure that we don't pass any knots or end-parameter
    next_limit = nextLimit(next_cvknot, end, step_forward);

//     cout << "\nI " << niter <<"  Cv "<< tx1 <<" "<< next_cvknot <<" "<< tstep <<"  Su "
// 	 << uinterval <<" "<< paru <<" "<< vinterval <<" "<< parv<< endl;

    if (hasPassedVal(tx1+tstep, next_limit, step_forward))
      tstep = next_limit-tx1;

    tx2 = tx1 + tstep;

    tx=0.5*(tx1+tx2);
    for (int ki=1; ki<=2; ++ki,tx=tx2) {
      curve->point(ft, &tx, nder);      

      //  Get closest point in the surface.
      guess[0] = paru;
      guess[1] = parv;
      spline_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			 REL_PAR_RES, &domain, guess);
      //		 REL_COMP_RES, NULL, guess);    @bsp
      spline_surf->basis_u().knotIntervalFuzzy(paru, fuzzy1);      
      spline_surf->basis_v().knotIntervalFuzzy(parv, fuzzy2);

      // Evaluate 0-2nd derivatives of the surface.
      spline_surf->point(gs, paru, parv,nder);
 
     
      // Check if point on curve and surface are within positional and
      // angular tolerances.
	    
      tdist = ft[0].dist(gs[0]);

//       cout << "M1 "<< ki <<" " << tx <<" "<< paru <<" "<< parv <<" "<< tdist
// 	   << endl;     

      if (tdist>aepsge) {
	// Points not within tolerances, curve and surface do not coincide.
	return 0;
      }
 
      // Indices to knot interval of the surface.
      uinterval = spline_surf->basis_u().knotInterval(paru);
      vinterval = spline_surf->basis_v().knotInterval(parv);


      // Check if any parameter lines of the surface is crossed in one of 
      // the parameter directions.
      //char supardir[]={'U','V'}; // @bsp for test print only
      for (int direction=1; direction<=2; ++direction) {

	SplineCurve* const_par_curve;
	if (direction == 1) {
	  if (uinterval == prev_uinterval)
	    continue;
	  // At least one parameter line is crossed. Fetch the constant
	  // parameter curve at the closest parameter line in the direction 
	  // of the marching.
	       
	  // Pick constant parameter curve.
	  if (uinterval > prev_uinterval)	  
	    paru = spline_surf->basis_u().begin()[prev_uinterval+1];
	  else if (uinterval < prev_uinterval)	  
	    paru = spline_surf->basis_u().begin()[prev_uinterval-1];

	  const_par_curve = spline_surf-> constParamCurve(paru, false);
	}
	else {
	  if (vinterval == prev_vinterval)
	    continue;
	  // At least one parameter line is crossed. Fetch the constant 
	  // parameter curve at the closest parameter line in the direction 
	  // of the marching.
	       
	  // Pick constant parameter curve.
	  if (vinterval > prev_vinterval)	  
	    parv = spline_surf->basis_v().begin()[prev_vinterval+1];
	  else if (vinterval < prev_vinterval)	  
	    parv = spline_surf->basis_v().begin()[prev_vinterval-1];
	  
	  const_par_curve = spline_surf-> constParamCurve(parv, true);
	}

	// Find the closest point between the input curve and the constant
	// parameter curve.
	ParamCurve* pcurve = curve->getParamCurve().get();
	seed1 = 0.5*(tx1+tx);      // @bsp
	//seed1 = 0.5*(tx1+tx2);
	seed2 = (direction==1)? parv : paru;
	  //seed2 = 0.5*(const_par_curve->startparam()+const_par_curve->endparam());
	ClosestPoint::closestPtCurves(pcurve, const_par_curve, pmin, pmax,
			const_par_curve->startparam(),
			const_par_curve->endparam(),
			seed1, seed2, par1, par2, tdist, ptc1, ptc2);
	delete const_par_curve;
 
// 	cout << seed1 <<" "<< seed2 <<" SP "<< par1<<" "<<  par2 << endl;
// 	cout << supardir[direction-1] << "1 "<< ki <<" " << par1 <<" "<< tdist
// 	     <<" "<< uinterval <<" "<< paru <<" "<< vinterval <<" "<< parv
// 	     << endl;    	
	
	//if (tdist>aepsge) {
	// Points not within tolerances, curve and surface do not coincide.
	//return 0;
	//}

	// Set new parameter values to the iteration.
	tx2 = par1;
	
       
	// Test midpoint of reduced step. First evaluate curve in midpoint.
	       
	tx = 0.5*(tx1 + tx2);
	curve->point(ft, &tx, nder);
	
	// Find closest point on surface.
	guess[0] = paru;   //@bsp
	guess[1] = parv;
	spline_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			       REL_PAR_RES, &domain, guess);
	//		   REL_COMP_RES, NULL, guess);    @bsp
	spline_surf->basis_u().knotIntervalFuzzy(paru, fuzzy1);      
	spline_surf->basis_v().knotIntervalFuzzy(parv, fuzzy2);     	

	// Calculate point and derivatives in surface.
	spline_surf->point(gs, paru, parv,nder);
	
	// Check if point on curve and surface are within positional and
	// angular tolerances.
	    
	tdist = ft[0].dist(gs[0]);
 
// 	cout << supardir[direction-1] << "2M "<< ki <<" " << tx <<" "
// 	     << paru <<" "<< parv <<" "<< tdist<< endl;

	if (tdist>aepsge) {

	  // Points not within tolerances, curve and surface do not coincide.
	  return 0;
	}

	// Calculate point and derivatives in the curve in the endpoint of the
	// step.
	       
	tx = tx2;
	curve->point(ft, &tx, nder);
	
	// Find closest point on surface.
	guess[0] = paru;   //@bsp
	guess[1] = parv;
	spline_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			       REL_PAR_RES, &domain, guess);
	//		   REL_COMP_RES, NULL, guess);    @bsp
	spline_surf->basis_u().knotIntervalFuzzy(paru, fuzzy1);      
	spline_surf->basis_v().knotIntervalFuzzy(parv, fuzzy2);     
	
	// Calculate point and derivatives in surface.
	spline_surf->point(gs, paru, parv,nder);
	
	// Check if point on curve and surface are within positional and
	// angular tolerances.
	    
	tdist = ft[0].dist(gs[0]);
 
// 	cout << supardir[direction-1] << "2E "<< ki <<" " << tx <<" "
// 	     << paru <<" "<< parv <<" "<< tdist<< endl;

	if (tdist>aepsge) {

	  // Points not within tolerances, curve and surface do not coincide.
	  return 0;
	}

	// Indices to knot interval of the surface.
	uinterval = spline_surf->basis_u().knotInterval(paru);
	vinterval = spline_surf->basis_v().knotInterval(parv);

	// Mark that a new step is to be initiated.
	ki = 2;
 
      } // direction=
 
    } // ki=
    
  } // while

  // Curve and surface are within tolerance
  return 1;
 
}


//===========================================================================
// Check if the two given surfaces coincide within the loop given
// by the parameter values of the intersection points representing it
int checkCoincide(ParamSurfaceInt *surf1, ParamSurfaceInt *surf2,
		  std::vector<double>& par_loop, 
		  const shared_ptr<GeoTol>  tol)
//===========================================================================
{
    // The number of parameter values in the loop must be divided by 4 
    // (the number of parameter directions in the two surfaces).
    ASSERT(par_loop.size() % 4 == 0);

    int nmb_sample = 10;  // Lower bound on number of curve-surface marches

    // Compute the length of the loop in the 4D parameter domain
    double len2 = 0.0;
    int nmb_par = (int)par_loop.size()/4;
    double *p1, *p2;
    double pmin, pcurr;
    int ki, kj, kmin;
    double* par_ptr = &par_loop[0];
    for (ki=1, p1=par_ptr+((ki-1)*4), p2=par_ptr+(ki*4); ki<nmb_par; 
	 ki++, p1=p2, p2=par_ptr+(ki*4))
	len2 += Utils::distance_squared(p1, p2, p2);

    // For each pair of intersection points (in sequence) check
    // coincidence orthogonal to the parameter direction defined by
    // the pair of intersection points by performing a number of
    // coincidence test marchings between a curve and a surface
    for (ki=1, p1=par_ptr+((ki-1)*4), p2=par_ptr+(ki*4); ki<nmb_par; 
	 ki++, p1=p2, p2=par_ptr+(ki*4))
    {
	// Decide the number of marching attempts for this pair
	double curr_len2 =  Utils::distance_squared(p1, p2, p2);
	int curr_sample = (int)(nmb_sample*curr_len2/len2)+1;

	// Find parameter direction
	double par[4];
	for (kmin=0, pmin=fabs(p2[0]-p1[0]), par[0]=0.5*(p1[0]+p2[0]), kj=1; 
	     kj<4; kj++)
	{
	    par[kj] = 0.5*(p1[kj]+p2[kj]);
	    pcurr = fabs(p2[kj]-p1[kj]);
	    if (pcurr < pmin)
	    {
		kmin = kj;
		pmin = pcurr;
	    }
	}

	// Set test surface and parameter values
	ParamSurfaceInt *curr_sf, *other_sf;
	//	RectDomain *domain;
	double tint, ta, tb, tc, td;
	double *par1, *par2, par_end[2];
	int idx;
	if (kmin < 2)
	{
	    idx = kmin;
	    par1 = par;
	    par2 = par+2;
	    curr_sf = surf1;
	    other_sf = surf2;
	}
	else
	{
	    idx = kmin-2;
	    par1 = par+2;
	    par2 = par;
	    curr_sf = surf2;
	    other_sf = surf1;
	}
	const RectDomain& other_domain = other_sf->getDegDomain(tol->getEpsge());

	ta = p1[kmin];
	tb = p2[kmin];
	tint = (tb - ta)/(double)(curr_sample+1);
	vector<double> mima = curr_sf->getMima();
	int idx2 = 1 - idx;
	for (par1[idx]=ta+0.5*tint, kj=0; kj<curr_sample; 
	     par1[idx]+=tint, kj++)
	{
	    // Define a surface curve across the current surface
	    // First set endpoint of curve
	    par_end[idx] = par1[idx];
	    par_end[idx2] = 
		(fabs(par1[idx2]-mima[2*idx2]) < fabs(par1[idx2]-mima[2*idx2+1])) 
		? mima[2*idx2+1] : mima[2*idx2];
	    tc = par1[idx2];
	    td = par_end[idx2];

	    // Make curve-on-surface curve
	    Point pp1(par1[0],par1[1]), pp2(par_end[0],par_end[1]);
	    shared_ptr<SplineCurve> pcrv;
	    if (tc < td)
		pcrv = (shared_ptr<SplineCurve>)
		(new SplineCurve(pp1, tc, pp2, td));
	    else
		pcrv = (shared_ptr<SplineCurve>)
		(new SplineCurve(pp2, td, pp1, tc));
	    shared_ptr<ParamCurve> constcrv(new CurveOnSurface(
						curr_sf->getParamSurface(), 
						pcrv, true));

	    // Project the startpoint of the curve onto the other surface
	    Point pt, clo_pt;
	    double clo_d, clo_uv[2];
	    curr_sf->point(pt, par1);
	    other_sf->getParamSurface()->closestPoint(pt, clo_uv[0], clo_uv[1],
						      clo_pt, clo_d, 
						      tol->getEpsge(), 
						      &other_domain, par2);
	    if (clo_d > tol->getEpsge())
		return 0; // Not coincidence

	    // March along the constant parameter curve until there
	    // is no coincidence or the a surface boundary is reached.
	    double last_pval, next_pval;
// 	    shared_ptr<ParamCurveInt> curve_int(new ParamCurveInt(constcrv));
	    bool coincident = false;
	    try {
		coincident = measure_coincidence_region(constcrv.get(), other_sf,
							tc, clo_uv, (tc < td), tol,
							last_pval, next_pval);
	    }
	    catch (...) {
		// May throw here. If so, keep coincident = 0. @@@vsk
		coincident = false;
	    }
	    if (!coincident)
		return 0;
	}

    }
    return 1;  // All instances of coincidence marching found coincidence
}


/*
//===========================================================================
// Denne versjonen av funksjonen har egen kode for de to parameter retningene
// naar vi sjekker om vi har passert skjoeter. (uinterval og vinterval)
// z.log
int checkCoincide(ParamCurveInt *curve,
				    double start, double end,
				    SplineSurfaceInt *intsurf,
				    const Point& su_start,
				    const Point& su_end)
//===========================================================================
{
  const double REL_PAR_RES=1.e-10;      // @bsp
  const double aepsge = 1.e-6;          // @bsp
  const double fuzzy_eps = 1.e-4;  // Knot distance tolerance.   @bsp
  const static int nder=2;   // Order of derivatives to be calulated.
  std::vector<Point> ft;     // Partial derivatives up to order nder
  std::vector<Point> gs;     // Partial derivatives up to order nder

  cout << " REL_PAR_RES=" << REL_PAR_RES << endl;

  shared_ptr<SplineSurface> spline_surf = intsurf->splineSurface();
  DEBUG_ERROR_IF(dim_ != spline_surf->dimension(),
	   "Dimension mismatch.");

  bool step_forward;    // True if we are going forward along this curve.
  int niter;            // Number of iterations.
  int kstat;            // Local status variable.
  int uinterval;        // Knot interval in 1. parameter direction of surface.
  int vinterval;        // Knot interval in 2. parameter direction of surface.
  int prev_uinterval;   // Previous knot interval in 1. parameter direction.
  int prev_vinterval;   // Previous knot interval in 2. parameter direction.
  double pmin, pmax;    // Searching range, this curve.
  double tstep;         // Step length.
  double tx, tx1, tx2;  // Parameter values of this curve.
  double next_cvknot;   // Next knot value, this curve.
  double next_limit;    // Next knot value or end parameter, this curve.
  double paru, parv;    // Parameter values of the surface.
  double tdist;         // Distance between two points.
  double guess[2];      // Start parameters for closest point iteration.
  double min_step;      // Minimum step length along the curve.
  double max_step;      // Maximum step length along the curve.

  double par1, par2;    // Parameter values to closest point between two curves
  double seed1, seed2;  // Start parameters for closest point iteration.
  double fuzzy1,fuzzy2; // Knot distance tolerance.
  Point ptc1, ptc2;     //  Closest points on two curves.

  Point clo_pt(dim_);   // Closest point.
  std::vector<Point> cvder(nder+1);  // Scratch vector for derivatives.

  // Tolerances for testing if a parameter is equal to a knot value
  fuzzy1 = fuzzy_eps*( spline_surf->basis_u().endparam() - 
		       spline_surf->basis_u().startparam());
  fuzzy2 = fuzzy_eps*( spline_surf->basis_v().endparam() - 
		       spline_surf->basis_v().startparam());

 // Make sure that we do not start or end outside the curve
  assureInRange(start); 
  assureInRange(end); 
 
  // Set searching range and directions.
  setRangeAndDirection(start, end, step_forward, pmin, pmax);

  // Set min and max step length along the curve.
  //@bsp Hvordan skal de settes? Hvis de bare er avhengige av kurven
  // kan de vaere medlem av klassen og settes i constructor.
  min_step = 0.01*(pmax-pmin);  //@bsp ??
  max_step = 0.1*(pmax-pmin);

 
  // Inital step along the curve
  tx1 = tx2 = start;
  paru = su_start[0];
  parv = su_start[1];

  // Evaluate 0-2nd derivatives of this curve.
  curve->point(ft, tx1, nder);

  // Evaluate 0-2nd derivatives of the surface.
  spline_surf->point(gs, paru, parv,nder);
 
  // Indices to knot interval of the surface.
  uinterval = spline_surf->basis_u().knotInterval(paru);
  vinterval = spline_surf->basis_v().knotInterval(parv);

  niter = 0;
  while (!hasReachedEnd(tx2, end, step_forward)) {
    niter++;
    tx1 = tx2;
    next_cvknot = nextSegmentVal(tx1, step_forward);
    prev_uinterval = uinterval;
    prev_vinterval = vinterval;
     
    tstep = stepLength(step_forward, cvder, min_step, max_step, aepsge);
    //kstat = stepLength(forward, tstep, cvder);  @ ???
 //ALWAYS_ERROR_IF(kstat<0,"Zero determinant in evalDistCurve",ComputationError());

    // Make sure that we don't pass any knots or end-parameter
    next_limit = nextLimit(next_cvknot, end, step_forward);

    cout << "\nI " << niter <<"  Cv "<< tx1 <<" "<< next_cvknot <<" "<< tstep <<"  Su "
	 << uinterval <<" "<< paru <<" "<< vinterval <<" "<< parv<< endl;

    if (hasPassedVal(tx1+tstep, next_limit, step_forward))
      tstep = next_limit-tx1;

    tx2 = tx1 + tstep;

    tx=0.5*(tx1+tx2);
    for (int ki=1; ki<=2; ++ki,tx=tx2) {
      curve->point(ft, tx, nder);      

      //  Get closest point in the surface.
      guess[0] = paru;
      guess[1] = parv;
      spline_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			 REL_PAR_RES, NULL, guess);
      //		 REL_COMP_RES, NULL, guess);    @bsp
      spline_surf->basis_u().knotIntervalFuzzy(paru, fuzzy1);      
      spline_surf->basis_v().knotIntervalFuzzy(parv, fuzzy2);

      // Evaluate 0-2nd derivatives of the surface.
      spline_surf->point(gs, paru, parv,nder);
 
     
      // Check if point on curve and surface are within positional and
      // angular tolerances.
	    
      tdist = ft[0].dist(gs[0]);

      cout << "M1 "<< ki <<" " << tx <<" "<< paru <<" "<< parv <<" "<< tdist
	   << endl;     

      if (tdist>aepsge) {
	// Points not within tolerances, curve and surface do not coincide.
	return 0;
      }
 
      // Indices to knot interval of the surface.
      uinterval = spline_surf->basis_u().knotInterval(paru);
      vinterval = spline_surf->basis_v().knotInterval(parv);


      // Check if any parameter lines of the surface is crossed in the 1. 
      // parameter direction.
 
      if (uinterval != prev_uinterval) {
	// At least one parameter line is crossed. Fetch the constant parameter
	// curve at the closest parameter line in the direction of the marching.
	       
	// Pick constant parameter curve.
	if (uinterval > prev_uinterval)	  
	  paru = spline_surf->basis_u().begin()[prev_uinterval+1];
	else if (uinterval < prev_uinterval)	  
	  paru = spline_surf->basis_u().begin()[prev_uinterval-1];

	SplineCurve* const_par_curve = spline_surf->constParamCurve(paru, false);

	// Find the closest point between the input curve and the constant
	// parameter curve.
	ParamCurve* pcurve = curve->getParamCurve().get();
	seed1 = 0.5*(tx1+tx);      // @bsp
	//seed1 = 0.5*(tx1+tx2);
	seed2 = parv;
	//seed2 = 0.5*(const_par_curve->startparam()+const_par_curve->endparam());
	ClosestPoint::closestPtCurves(pcurve, const_par_curve, pmin, pmax,
			const_par_curve->startparam(),
			const_par_curve->endparam(),
			seed1, seed2, par1, par2, tdist, ptc1, ptc2);
	delete const_par_curve;

 	cout << seed1 <<" "<< seed2 <<" SP "<< par1<<" "<<  par2 << endl;
	cout << "U1 "<< ki <<" " << par1 <<" "<< tdist <<" "<< uinterval 
	     <<" "<< paru <<" "<< vinterval <<" "<< parv << endl;

	//if (tdist>aepsge) {
	  // Points not within tolerances, curve and surface do not coincide.
	  //return 0;
	//}

	// Set new parameter values to the iteration.
	tx2 = par1;

       
	// Test midpoint of reduced step. First evaluate curve in midpoint.
	       
	tx = 0.5*(tx1 + tx2);
	curve->point(ft, tx, nder);
	
	// Find closest point on surface.
	guess[0] = paru;   //@bsp
	guess[1] = parv;
	spline_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			       REL_PAR_RES, NULL, guess);
	//		   REL_COMP_RES, NULL, guess);    @bsp
	spline_surf->basis_u().knotIntervalFuzzy(paru, fuzzy1);      
	spline_surf->basis_v().knotIntervalFuzzy(parv, fuzzy2);     
	
	// Calculate point and derivatives in surface.
	spline_surf->point(gs, paru, parv,nder);
	
	// Check if point on curve and surface are within positional and
	// angular tolerances.
	    
	tdist = ft[0].dist(gs[0]);
 
	cout << "U2M "<< ki <<" " << tx <<" "<< paru <<" "<< parv <<" "<< tdist<< endl;

	if (tdist>aepsge) {

	  // Points not within tolerances, curve and surface do not coincide.
	  return 0;
	}


	// Calculate point and derivatives in the curve in the endpoint of the
	// step.
	       
	tx = tx2;
	curve->point(ft, tx, nder);
	
	// Find closest point on surface.
	guess[0] = paru;   //@bsp
	guess[1] = parv;
	spline_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			       REL_PAR_RES, NULL, guess);
	//		   REL_COMP_RES, NULL, guess);    @bsp
	spline_surf->basis_u().knotIntervalFuzzy(paru, fuzzy1);      
	spline_surf->basis_v().knotIntervalFuzzy(parv, fuzzy2);     
	
	// Calculate point and derivatives in surface.
	spline_surf->point(gs, paru, parv,nder);
	
	// Check if point on curve and surface are within positional and
	// angular tolerances.
	    
	tdist = ft[0].dist(gs[0]);
 
	cout << "U2E "<< ki <<" " << tx <<" "<< paru <<" "<< parv <<" "<< tdist<< endl;

	if (tdist>aepsge) {

	  // Points not within tolerances, curve and surface do not coincide.
	  return 0;
	}


      // Indices to knot interval of the surface.
      uinterval = spline_surf->basis_u().knotInterval(paru);
      vinterval = spline_surf->basis_v().knotInterval(parv);

	// Mark that a new step is to be initiated.
	ki = 2;
      }


      // Check if any parameter lines of the surface is crossed in the 2. 
      // parameter direction.
      
      if (vinterval != prev_vinterval) {
	// At least one parameter line is crossed. Fetch the constant parameter
	// curve at the closest parameter line in the direction of the marching.
	
	// Pick constant parameter curve.
	if (vinterval > prev_vinterval)	  
	  parv = spline_surf->basis_v().begin()[prev_vinterval+1];
	else if (vinterval < prev_vinterval)	  
	  parv = spline_surf->basis_v().begin()[prev_vinterval-1];

	SplineCurve* const_par_curve = spline_surf->constParamCurve(parv, true);
	// Find the closest point between the input curve and the constant
	// parameter curve.
	ParamCurve* pcurve = curve->getParamCurve().get();
	seed1 = 0.5*(tx1+tx);      // @bsp
	//seed1 = 0.5*(tx1+tx2);
	seed2 = paru;
	//seed2 = 0.5*(const_par_curve->startparam()+const_par_curve->endparam());
	ClosestPoint::closestPtCurves(pcurve, const_par_curve, pmin, pmax,
			const_par_curve->startparam(),
			const_par_curve->endparam(),
			seed1, seed2, par1, par2, tdist, ptc1, ptc2);
	delete const_par_curve;

	cout << seed1 <<" "<< seed2 <<" SP "<< par1<<" "<<  par2 << endl;
	cout << "V1 "<< ki <<" " << par1 <<" "<< tdist <<" "<< uinterval
	     <<" "<< paru <<" "<< vinterval <<" "<< parv << endl;
	
	//if (tdist>aepsge) {
	  // Points not within tolerances, curve and surface do not coincide.
	  //return 0;
	//}
	
	// Set new parameter values to the iteration.
	tx2 = par1;
	
	
	// Test midpoint of reduced step. First evaluate curve in midpoint.
	
	tx = 0.5*(tx1 + tx2);
	curve->point(ft, tx, nder);
	
	// Find closest point on surface.
	guess[0] = paru;   //@bsp
	guess[1] = parv;
	spline_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			       REL_PAR_RES, NULL, guess);
	//		   REL_COMP_RES, NULL, guess);    @bsp
	spline_surf->basis_u().knotIntervalFuzzy(paru, fuzzy1);      
	spline_surf->basis_v().knotIntervalFuzzy(parv, fuzzy2);     	
	
	// Calculate point and derivatives in surface.
	spline_surf->point(gs, paru, parv,nder);
	
	// Check if point on curve and surface are within positional and
	// angular tolerances.
	
	tdist = ft[0].dist(gs[0]);

	cout << "V2M "<< ki <<" " << tx <<" "<< paru <<" "<< parv <<" "<< tdist<< endl;

	if (tdist>aepsge) {
	  // Points not within tolerances, curve and surface do not coincide.
	  return 0;
	}


	// Calculate point and derivatives in the curve in the endpoint of the
	// step

	tx = tx2;
	curve->point(ft, tx, nder);
	
	// Find closest point on surface.
	guess[0] = paru;   //@bsp
	guess[1] = parv;
	spline_surf->closestPoint(ft[0],paru,parv,clo_pt,tdist,
			       REL_PAR_RES, NULL, guess);
	//		   REL_COMP_RES, NULL, guess);    @bsp
	spline_surf->basis_u().knotIntervalFuzzy(paru, fuzzy1);      
	spline_surf->basis_v().knotIntervalFuzzy(parv, fuzzy2);     	
	
	// Calculate point and derivatives in surface.
	spline_surf->point(gs, paru, parv,nder);
	
	// Check if point on curve and surface are within positional and
	// angular tolerances.
	
	tdist = ft[0].dist(gs[0]);

	cout << "V2E "<< ki <<" " << tx <<" "<< paru <<" "<< parv <<" "<< tdist<< endl;

	if (tdist>aepsge) {
	  // Points not within tolerances, curve and surface do not coincide.
	  return 0;
	}

	// Indices to knot interval of the surface.
	uinterval = spline_surf->basis_u().knotInterval(paru);
	vinterval = spline_surf->basis_v().knotInterval(parv);	

	// Mark that a new step is to be initiated.
	ki = 2;
	
      } // crossed2
    } // ki=
    
  } // while

  //  Curve and surface are within tolerance
  return 1;
 
}
*/

//===========================================================================
void evalDistCurve(const vector<Point>& marching_curve, 
		   const vector<Point>& other_curve, 
		   vector<Point>& dist_curve)
//===========================================================================
{
    ASSERT(marching_curve.size() == 3 && other_curve.size() == 3);
    dist_curve.resize(3);

    dist_curve.resize(3);
    dist_curve[0] = marching_curve[0] - other_curve[0];

    const Point& tan_march = marching_curve[1];
    const Point& tan_other = other_curve[1];
    
    double other_l2 = tan_other.length2();
    double scal_prod = tan_march * tan_other;
    
    // this is the factor we need to multiply the other curve with in order
    // to get a projection of its tangent onto the marching curve with the same
    // length as the tangent of the marching curve
    double fac = scal_prod / other_l2; 
    
    dist_curve[1] = marching_curve[1] - fac * other_curve[1];
    dist_curve[2] = marching_curve[2] - fac * fac * other_curve[2];
}

//===========================================================================
int stepLength(const vector<Point>& ft, const vector<Point>& gs, 
	       double delta, bool forward, double& delta_t)
//===========================================================================
{
    vector<Point> D(3); // holder for the distance curve differential information
    evalDistCurve(ft, gs, D);

    double D1_norm2 = D[1].length2();
    double D2_norm2 = D[2].length2();

    const double VANISHING_LIMIT = 1.0e-8; // @ what is a reasonable number for this??
    const double RAD_EPS = 1.0e-4; // @ as taken from the older stepLength routine.
    delta = fabs(delta) - D[0].length();
    DEBUG_ERROR_IF(delta < 0, "Logical error in stepLength()");

    if (4 * D1_norm2 * D1_norm2 > D2_norm2 * delta * delta) {
	// the second derivative term of the taylor series is the most important
	//delta_t = delta / sqrt(D1_norm2);
	delta_t = sqrt(delta / D1_norm2);
	delta_t *= forward ? 1 : -1;  // changing orientation if not forward
    } else if (D2_norm2 > VANISHING_LIMIT) {
	// the fourth derivative term of the taylor series is the most important
	delta_t = 2 * sqrt(delta / D2_norm2);
	delta_t = sqrt(delta_t);
	delta_t *= forward ? 1: -1; // changing orientation if not forward
    } else {
	// both the second and the fourth derivative term have all but vanished.
	// Seek solace in alternative solutions.
	if (ft[0].dimension() == 1) {
	    delta_t = 100.0*RAD_EPS; // @@sbr Well, curv not defined, but defined in (u, v, f(u,v)).
	    delta_t *= forward ? 1.0 : -1.0;
	} else {
	    vector<Point> dummy(3);
	    double rad = std::min(curvatureRadius(ft, dummy), curvatureRadius(gs, dummy));
	    delta_t = stepLenFromRadius(rad, RAD_EPS);
	    delta_t *= forward ? 1 : -1;
	}
    }
    return 1;
}



//===========================================================================
double stepLength(const vector<Point>& ft, 
		  const vector<Point>& gs, 
		  bool forward, 
		  std::vector<Point>& cvder,
		  double min_step, 
		  double max_step, 
		  double aepsge)
//===========================================================================
{ 
    // Doesn't make sense to call this function if ft[0] and gs[2] are
    // not approx. coincident
    //ASSERT(ft[0].dist(gs[0]) < aepsge); 

  // Evaluate shortest distance curve.
  double result; 
  int kstat = evalProjectedCurve(ft, gs, cvder);
  //ALWAYS_ERROR_IF(kstat<0,"Zero determinant in evalProjectedCurve");
   if (kstat < 0)
   {
       double tmp_step;
       double len = ft[1].length();
       if (len > 0.001*aepsge)
	   tmp_step = 2.0*aepsge/len;
       else
	   tmp_step = 0.002*aepsge;
       result = forward ? tmp_step : -tmp_step;
   }
   else
   {
  
  // 'cvder' is now representing the projection of the curve onto the surface.

       kstat = stepLength(ft, cvder, aepsge, forward, result);
  //ALWAYS_ERROR_IF(kstat<1,"Error occured when calculating step length");
       if (kstat < 1)
       {
       double tmp_step;
       double len = ft[1].length();
       if (len > 0.001*aepsge)
	   tmp_step = 2.0*aepsge/len;
       else
	   tmp_step = 0.002*aepsge;
       result = forward ? tmp_step : -tmp_step;
//	   result = forward ? aepsge : -aepsge;
   }
   }
   if (result < min_step) {
      result = min_step;
  } else if (result > max_step) {
      result = max_step;
      }
  return result;
}


//===========================================================================
int evalProjectedCurve(const vector<Point>& ft, const vector<Point>& gs, 
		       vector<Point>& res)
//===========================================================================
{ 
  // Evaluate the curve representing the shortest distance between a curve
  // and a surface in the current point.
    res.resize(3);

    res[0] = gs[0]; // a point on the surface

    // Make references to derivatives
    const Point& dfdt = ft[1]; const Point& d2fdt2 = ft[2];
    const Point& dSdu   = gs[1]; const Point& dSdv = gs[2];
    const Point& d2Sdu2  = gs[3]; const Point& d2Sdudv   = gs[4]; 
    const Point& d2Sdv2 = gs[5];  
    
    // intermediate values
    double dSdu_length2 = dSdu.length2();
    double dSdv_length2 = dSdv.length2();
    double dSdu_dSdv = dSdu * dSdv;
    double dfdt_dSdu    = dfdt * dSdu;
    double dfdt_dSdv    = dfdt * dSdv;
    
    // The determinant
    //double det = (dSdu*dSdu)*(dSdv*dSdv) - (dSdu*dSdv)*(dSdu*dSdv);
    double det = dSdu_length2 * dSdv_length2 - dSdu_dSdv * dSdu_dSdv;
    if (det == 0.0)
	return -1;
    double det_inv = 1 / det;
    
    //Use Cramer's rule to compute the 2x2 linear system, where the solution
    //is the derivative of the surface in the parameter plane.
    
    // First derivatives     
    double u1 = ((dSdv_length2 * dfdt_dSdu) - (dSdu_dSdv * dfdt_dSdv)) * det_inv;
    double v1 = ((dSdu_length2 * dfdt_dSdv) - (dSdu_dSdv * dfdt_dSdu)) * det_inv;
    res[1] = u1 * dSdu - v1 * dSdv;
    
    // Second derivatives
    Point temp = u1 * u1 * d2Sdu2 + 2 * u1 * v1 * d2Sdudv + v1 * v1 * d2Sdv2;
    Point vec = d2fdt2 - temp;
    double vec_dSdu = vec * dSdu;
    double vec_dSdv = vec * dSdv;
    double u2 = ( dSdv_length2 * (vec_dSdu) - dSdu_dSdv * (vec_dSdv) ) * det_inv;
    double v2 = ( dSdu_length2 * (vec_dSdv) - dSdu_dSdv * (vec_dSdu) ) * det_inv;
    res[2] = temp + u2 * dSdu + v2 * dSdv;
    
    return 0;
}

//===========================================================================
bool determineCoincidenceRegion(const ParamObjectInt* obj1,
				const ParamObjectInt* obj2,
				shared_ptr<const GeoTol> tol,
				const double* current_params,
				int dir,
				bool forward,
				double& last_param_val_inside,
				double& first_param_val_outside)
//===========================================================================
{
    // assure that we are going to march along the first object
    int total_num_param = obj1->numParams() + obj2->numParams();
    const double* param_pointer = 0;
    vector<double> temp;
    if (dir >= obj1->numParams()) { // 'dir' is a parameter in the second object
	// swapping parameters
	temp = vector<double>(current_params, current_params + total_num_param);
	rotate(temp.begin(), temp.begin() + obj1->numParams(), temp.end());
	param_pointer = &temp[0];
	// swapping 'direction'
	dir -= obj1->numParams();
	// swapping objects
	swap(obj1, obj2);
    } else {
	param_pointer = current_params;
    }
    ASSERT(dir < obj1->numParams()); // now, 'dir' should be in the first object

    const ParamCurve* first_cv = NULL;
    shared_ptr<const ParamCurve> temp_curve;

    const GeomObject* obj2_geom = NULL;
    
    if (obj1->numParams() > 1) {
	ASSERT(obj1->numParams() == 2);
	// this should be a surface.  We must make the proper curve
	int fixed_dir = (dir + 1) % 2;
	temp_curve = make_surf_curve_preserve_critical(obj1, dir, param_pointer[fixed_dir]);
	first_cv = temp_curve.get();
    } else {
	// this is a curve
	const ParamCurveInt* tempcv = dynamic_cast<const ParamCurveInt*>(obj1);
	if (tempcv != 0) {
	    first_cv = tempcv->getParamCurve().get();
	} else { // It could also be a function.
	    const Param1FunctionInt* tempfunc = dynamic_cast<const Param1FunctionInt*>(obj1);
	    if (tempfunc != 0)
		first_cv = tempfunc->getParamCurve().get();
	}
	ASSERT(first_cv != 0);
    }

    if (obj2->numParams() > 1) {
	ASSERT(obj2->numParams() == 2);
	// this should be a surface.  We must make the proper curve
	//	int fixed_dir = (dir + 1) % 2;
	const ParamSurfaceInt* tempsurf = dynamic_cast<const ParamSurfaceInt*>(obj2);
	DEBUG_ERROR_IF(!tempsurf, "Unrecognized surface format.");

	obj2_geom = tempsurf->getParamSurface().get();
    } else if (obj2->numParams() > 0) { // IntersectorFuncConst has 0-dim obj2.
	// this is a curve
	const ParamCurveInt* tempcv = dynamic_cast<const ParamCurveInt*>(obj2);
	if (tempcv != 0) {
	    obj2_geom = tempcv->getParamCurve().get();
	} else { // It could also be a function.
	    const Param1FunctionInt* tempfunc = dynamic_cast<const Param1FunctionInt*>(obj2);
	    if (tempfunc != 0)
		obj2_geom = tempfunc->getParamCurve().get();
	}
	ASSERT(obj2_geom != 0);
    }

    return measure_coincidence_region(first_cv, 
				      dynamic_cast<const ParamGeomInt*>(obj2),
				      param_pointer[dir],
				      &param_pointer[obj1->numParams()],
				      forward,
				      tol,
				      last_param_val_inside,
				      first_param_val_outside);
}


//===========================================================================
bool determineCoincidenceRegion(const ParamFunctionInt* obj_1d,
				double C,
				shared_ptr<const GeoTol> tol,
				const double* current_params,
				int dir,
				bool forward,
				double& last_param_val_inside,
				double& first_param_val_outside)
//===========================================================================
{
    ASSERT(dir < obj_1d->numParams()); // now, 'dir' should be in the first object
    // assure that we are going to march along the first object
    const double* param_pointer = 0;
    vector<double> temp;
//     if (dir >= obj1->numParams()) { // 'dir' is a parameter in the second object
// 	// swapping parameters
// 	temp = vector<double>(current_params, current_params + total_num_param);
// 	rotate(temp.begin(), temp.begin() + obj1->numParams(), temp.end());
// 	param_pointer = &temp[0];
// 	// swapping 'direction'
// 	dir -= obj1->numParams();
// 	// swapping objects
// 	swap(obj1, obj2);
//     } else {
	param_pointer = current_params;
//     }

    const ParamCurve* first_cv = NULL;
    shared_ptr<const ParamCurve> temp_curve;

    if (obj_1d->numParams() > 1) {
	ASSERT(obj_1d->numParams() == 2);
	// this should be a surface.  We must make the proper curve
	int fixed_dir = (dir + 1) % 2;
	temp_curve = make_surf_curve_preserve_critical(obj_1d, dir, param_pointer[fixed_dir]);
	first_cv = temp_curve.get();
    } else {
	// this is a curve
	const ParamCurveInt* tempcv = dynamic_cast<const ParamCurveInt*>(obj_1d);
	if (tempcv != 0) {
	    first_cv = tempcv->getParamCurve().get();
	} else { // It could also be a function.
	    const Param1FunctionInt* tempfunc = dynamic_cast<const Param1FunctionInt*>(obj_1d);
	    if (tempfunc != 0)
		first_cv = tempfunc->getParamCurve().get();
	}
	ASSERT(first_cv != 0);
    }

//     Point other_pt(1);
//     other_pt[0] = C;
    return measure_coincidence_region(first_cv,
				      C,
				      param_pointer[dir],
				      &param_pointer[obj_1d->numParams()],
				      forward,
				      tol,
				      last_param_val_inside,
				      first_param_val_outside);
}


} // end namespace

namespace {

// 	shared_ptr<const SplineSurface> spline_sf;
// 	const ParamSurfaceInt* tempsurf = dynamic_cast<const ParamSurfaceInt*>(obj1);
// 	bool critical_values;
// 	if (tempsurf == 0) {
// 	    const Param2FunctionInt* tempfunc = dynamic_cast<const Param2FunctionInt*>(obj1);
// 	    spline_sf = dynamic_pointer_cast<const SplineSurface, const ParamSurface>
// 		(tempfunc->getSurface());
// 	    critical_values = tempfunc->hasCriticalVals(dir);
// 	} else {
// 	    spline_sf = dynamic_pointer_cast<const SplineSurface, const ParamSurface>
// 		(tempsurf->getSurface());
// 	    critical_values = tempsurf->hasCriticalVals(dir);
// 	}
// 	ASSERT(spline_sf.get() != NULL);

//===========================================================================
shared_ptr<ParamCurve>
make_surf_curve_preserve_critical(const ParamObjectInt* pobj,
				  int running_dir,
				  double fixed_value)
//===========================================================================
{
    ASSERT(pobj->numParams() == 2); // should only be run on 2-manifold (surf. or function)
    ASSERT(running_dir == 0 || running_dir == 1);
    vector<double> crit_vals = pobj->getCriticalValsAndKnots(running_dir);
    
    // extracting underlying spline surface
    shared_ptr<const ParamSurface> param_surf;
    const ParamSurfaceInt* tempsurf = dynamic_cast<const ParamSurfaceInt*>(pobj);
    const Param2FunctionInt* tempfunc = dynamic_cast<const Param2FunctionInt*>(pobj);
    if (tempsurf) {
	// this was a surface interface
	param_surf = tempsurf->getParamSurface();
    } else if (tempfunc) {
	// this was a function interface
	param_surf = tempfunc->getParamSurface();
    } else {
	THROW("Unrecognized ParamObjectInt type encountered.");
    }
    DEBUG_ERROR_IF(param_surf.get() == 0, "Could not get hold of ParamSurface from the supplied "
 	                            "ParamObjectInt to make_surf_curve_preserve_critical.");
    vector<double> knots, coefs;
    int fixed_dir = (running_dir + 1) % 2;
    if (int(crit_vals.size()) > 0) {
	ASSERT(crit_vals.size() >= 2); // we must be able to deduce an interval from this vector
	// making knotvector with 2-multiple knots at ends (for an order 2 curve)
	knots.push_back(crit_vals[0]);
	knots.insert(knots.end(), crit_vals.begin(), crit_vals.end());
	knots.push_back(knots.back());
    } else {
	knots.resize(4);
	RectDomain dom = param_surf->containingDomain();
	knots[0] = knots[1] = (running_dir == 0) ? dom.umin() : dom.vmin();
	knots[2] = knots[3] = (running_dir == 0) ? dom.umax() : dom.vmax();
    }
    const int num_coefs = (int)knots.size() - 2;
    coefs.resize(num_coefs * 2);
    for (int i = 0; i < num_coefs; ++i) {
	coefs[2 * i + fixed_dir] = fixed_value;
	coefs[2 * i + running_dir] = knots[i + 1];
    }
    shared_ptr<SplineCurve> curve(new SplineCurve(num_coefs, 2, &knots[0], &coefs[0], 2));
    shared_ptr<const ParamSurface> surface = param_surf;
    shared_ptr<ParamSurface> unconsted_surf = const_pointer_cast<ParamSurface>(surface);
    shared_ptr<ParamCurve> surfcurve(new CurveOnSurface(unconsted_surf, curve, true));
    return surfcurve;
}

//===========================================================================
bool measure_coincidence_region(const ParamCurve* marching_curve,
				const ParamGeomInt* other_obj,
				double isect_param_marching,
				const double* const isect_param_other,
				bool step_forward,
				const shared_ptr<const GeoTol> gtol,
				double& last_pval_inside_marching,
				double& first_pval_outside_marching)
//===========================================================================
{
//     double eps = gtol->getEpsge();
    double eps_bracket;
    
    // argument checks
    DEBUG_ERROR_IF(!marching_curve || !other_obj,
		   "Null pointers detected.");
    eps_bracket = gtol->getEpsBracket();
    DEBUG_ERROR_IF(eps_bracket <= 0 || eps_bracket >= 1,
		   "wrong value for eps_bracket.");

    // initializing
    bool distance_exceeds_tolerance = false;
    double prev_distance, current_distance = 0;
    double first_pval_outside_other[2];
    double last_pval_inside_other[2];
    int num_param; // = other_obj->numParams();
    first_pval_outside_marching = isect_param_marching;

    double end_march = step_forward
	? marching_curve->endparam() : marching_curve->startparam();
    
    // choosing function according to type of second object
    ParamGeomInt* other_tmp = (ParamGeomInt*)(other_obj);
    shared_ptr<BindToObject> bso;
    ParamSurfaceInt* surf_2nd = other_tmp->getParamSurfaceInt();
    const ParamCurve* curve_2nd
	= surf_2nd ? 0 : other_tmp->getParamCurveInt()->getParamCurve().get();
    if (curve_2nd) {
	num_param = 1;
	bso = shared_ptr<BindToCurve>
	    (new BindToCurve(marching_curve, curve_2nd,
			     step_forward, gtol));
	// we suppose that the distance is within geometric tolerance
	// before expanding coincidence region
	DEBUG_ERROR_IF(marching_curve->point(first_pval_outside_marching).
		       dist(curve_2nd->point(isect_param_other[0]))
		       > gtol->getEpsge(),
		       "IntersectionPoint's objects did not lie within "
		       "expected tolerance");
    } else if (surf_2nd) {
	num_param = 2;
	bso = shared_ptr<BindToSurface>
	    (new BindToSurface(marching_curve, surf_2nd,
			       step_forward, gtol));
	// we suppose that the distance is within geometric tolerance
	// before expanding coincidence region
	Point other_pt;
	surf_2nd->point(other_pt, isect_param_other);
	DEBUG_ERROR_IF(marching_curve->point(first_pval_outside_marching).
		       dist(other_pt)
		       > gtol->getEpsge(),
		       "IntersectionPoint's objects did not lie within "
		       "expected tolerance");
    } else { 
	THROW("Unrecognized/unsupported ParamObjectInt type.");
    }
    copy(isect_param_other, isect_param_other + num_param,
	 first_pval_outside_other);
    copy(isect_param_other, isect_param_other + num_param,
	 last_pval_inside_other);

    // iterating until lower bracket is passed, or end of marching
    // curve reached
    while(!distance_exceeds_tolerance) {
	// saving old parameter values
	last_pval_inside_marching = first_pval_outside_marching;
	prev_distance = current_distance;
	bso->expandCoincidenceRegion(first_pval_outside_marching,
				     first_pval_outside_other,
				     current_distance,
				     distance_exceeds_tolerance);
	if (fabs(first_pval_outside_marching-end_march)
	    < gtol->getRelParRes()) {
	    break; // arrived at boundary
	}
    }
    if (!distance_exceeds_tolerance) {
	// we arrived at end of marching curve without passing the
	// lower limit
	last_pval_inside_marching = first_pval_outside_marching;
	return false;
    }
    
    // if we got here, we have found a point on the curve that is
    // outside the influence area.  We now want a tighter estimate of
    // the locality of this area.
     bso->determineConfidenceInterval(last_pval_inside_marching,
				      last_pval_inside_other,
				      first_pval_outside_marching,
				      first_pval_outside_other,
				      prev_distance,
				      current_distance);
     return true;
}

//===========================================================================
bool measure_coincidence_region(const ParamCurve* marching_curve_1d,
				double C,
				double isect_param_marching,
				const double* const isect_param_other,
				bool step_forward,
				const shared_ptr<const GeoTol> gtol,
				double& last_pval_inside_marching,
				double& first_pval_outside_marching)
//===========================================================================
{
    ASSERT(marching_curve_1d->dimension() == 1);

//     double eps = gtol->getEpsge();
    double eps_bracket;
    
    // argument checks
    DEBUG_ERROR_IF(!marching_curve_1d,
		   "Null pointer detected.");
    eps_bracket = gtol->getEpsBracket();
    DEBUG_ERROR_IF(eps_bracket <= 0 || eps_bracket >= 1,
		   "wrong value for eps_bracket.");

    // initializing
    bool distance_exceeds_tolerance = false;
    double prev_distance, current_distance = 0;
    double first_pval_outside_other[2];
    double last_pval_inside_other[2];
//     int num_param;
    first_pval_outside_marching = isect_param_marching;

    double end_march = step_forward
	? marching_curve_1d->endparam() : marching_curve_1d->startparam();
    
    // choosing function according to type of second object
    Point C_pt(1);
    C_pt[0] = C;
    shared_ptr<BindToObject> bso =
	shared_ptr<BindToPoint>(new BindToPoint(marching_curve_1d,
						C_pt, step_forward, gtol));
//     const ParamSurface* surf_2nd
// 	= dynamic_cast<const ParamSurface*>(other_obj);
//     const ParamCurve* curve_2nd
// 	= surf_2nd ? 0 : dynamic_cast<const ParamCurve*>(other_obj);
//     if (curve_2nd) {
// 	num_param = 1;
// 	bso = shared_ptr<BindToCurve>
// 	    (new BindToCurve(marching_curve_1d, curve_2nd,
// 			     step_forward, gtol));
//     } else if (surf_2nd) {
// 	num_param = 2;
// 	bso = shared_ptr<BindToSurface>
// 	    (new BindToSurface(marching_curve_1d, surf_2nd,
// 			       step_forward, gtol));
//     } else { 
// 	THROW("Unrecognized/unsupported ParamObjectInt type.");
//     }
//     copy(isect_param_other, isect_param_other + num_param,
// 	 first_pval_outside_other);
//     copy(isect_param_other, isect_param_other + num_param,
// 	 last_pval_inside_other);

    // iterating until lower bracket is passed, or end of marching
    // curve reached
    while(!distance_exceeds_tolerance) {
	// saving old parameter values
	last_pval_inside_marching = first_pval_outside_marching;
	prev_distance = current_distance;
	bso->expandCoincidenceRegion(first_pval_outside_marching,
				     first_pval_outside_other,
				     current_distance,
				     distance_exceeds_tolerance);
	if (fabs(first_pval_outside_marching-end_march)
	    < gtol->getRelParRes()) {
	    break; // arrived at boundary
	}
    }
    if (!distance_exceeds_tolerance) {
	// we arrived at end of marching curve without passing the
	// lower limit
	last_pval_inside_marching = first_pval_outside_marching;
	return false;
    }
    
    // if we got here, we have found a point on the curve that is
    // outside the influence area.  We now want a tighter estimate of
    // the locality of this area.
     bso->determineConfidenceInterval(last_pval_inside_marching,
				      last_pval_inside_other,
				      first_pval_outside_marching,
				      first_pval_outside_other,
				      prev_distance,
				      current_distance);
     return true;
}

//===========================================================================
void BindToObject::determineConfidenceInterval(double& last_pval_inside_curve,
					       double* last_pval_inside_other,
					       double& first_pval_outside_curve,
					       double* first_pval_outside_other,
					       double& inside_dist,
					       double& outside_dist)
//===========================================================================
{
    const double& center_value = tol_->getEpsge();
    const double& bracket = tol_->getEpsBracket();
    const double lower_limit = (1 - bracket) * center_value;
    const double upper_limit = (1 + bracket) * center_value;
    
    ASSERT(inside_dist < center_value && outside_dist > center_value);

    bool outside_val_fixed = (outside_dist < upper_limit);
    bool inside_val_fixed = (inside_dist > lower_limit);
    
    Point p1, p2;
    double new_dist;

    // we knot that 'center_value' is bracketed by 'inside_dist' and 'outside_dist'.
    // Now we need to make the bracket narrow enough to fit inside [lower_limit, upper_limit].

    int num_param_other = numParametersOtherObject();
    ASSERT(num_param_other <= 2);
    int nmb_iter = 0;
    while (!inside_val_fixed || !outside_val_fixed) {
	nmb_iter++;
	double mid_val = (inside_val_fixed) ? 0.9*center_value + 0.1*upper_limit :
	    0.9*center_value + 0.1*lower_limit;
	//((outside_val_fixed) ? 0.3*center_value + 0.7*lower_limit : center_value);
	double frac1 = (outside_dist - mid_val)/(outside_dist - inside_dist);
	double frac2 = (mid_val - inside_dist)/(outside_dist - inside_dist);
	if (nmb_iter > 5)
	    frac1 = frac2 = 0.5;
	double mid_p_marching = frac1*last_pval_inside_curve + frac2*first_pval_outside_curve;
	double mid_p_other[2];
	for (int i = 0; i < num_param_other; ++i) {
	    mid_p_other[i] = frac1*last_pval_inside_other[i] + frac2*first_pval_outside_other[i];
	}
	march_cv_->point(p1, mid_p_marching);
	other_closest_point(p1, mid_p_other, new_dist, p2, true);
	if (new_dist >= outside_dist || new_dist <= inside_dist) {
	    // either the distance function is ill-behaved close to this point ("wild" surface)
	    // or we are having numerical problems with the closest point routine.  In
	    // this case, we will satify ourselves with the values already found, even though 
	    // they may define a larger interval than ideally sought.
	    return;
	}
	if (new_dist > center_value) {
	    first_pval_outside_curve = mid_p_marching;
	    copy(mid_p_other, mid_p_other + num_param_other, first_pval_outside_other);
	    outside_dist = new_dist;
	    outside_val_fixed = (outside_dist <= upper_limit);
	} else { // new_dist <= center_value
	    last_pval_inside_curve = mid_p_marching;
	    copy(mid_p_other, mid_p_other + num_param_other, last_pval_inside_other);
	    inside_dist = new_dist;
	    inside_val_fixed = (new_dist >= lower_limit);
	}
    }
}

//===========================================================================
void BindToPoint::other_closest_point(const Point& pt,
				      double* cl_par,
				      double& cl_dist,
				      Point& cl_pt,
				      bool mind_segment_values)
//===========================================================================
{
    cl_dist = pt.dist(other_point_);
    cl_pt = other_point_;
}

//===========================================================================
void BindToCurve::other_closest_point(const Point& pt,
				      double* cl_par,
				      double& cl_dist,
				      Point& cl_pt,
				      bool mind_segment_values)
//===========================================================================
{
    const double lower_bound = other_cv_->startparam();
    const double upper_bound = other_cv_->endparam();
    double& par = cl_par[0];
    const double& pareps = std::max(tol_->getRelParRes(),
				    tol_->getRelParRes() * (other_cv_->endparam() - other_cv_->startparam()));

    if (!mind_segment_values) { // @@sbr Not the opposite ... ?
	other_cv_->closestPoint(pt,          // reference point
				lower_bound, // start of search interval
				upper_bound, // end of search interval
				par,         // parameter value
				cl_pt,       // pos. of closest pt.
				cl_dist,     // distance to reference pt.
				&par);       // seed
    } else {
	// we must consider possible problems near critical values
	double pmin = std::max(other_cv_->nextSegmentVal(cl_par[0], false, pareps), lower_bound);
	double pmax = std::min(other_cv_->nextSegmentVal(cl_par[0], true, pareps), upper_bound);

	const ParamCurve* pc = other_cv_; //->getParamCurve();
	bool repeat;
	double prev_dist = cl_dist = LARGEST_POSSIBLE_VALUE;
	do {
	    repeat = false;
	    pc->closestPoint(pt, pmin, pmax, *cl_par, cl_pt, cl_dist, cl_par);
	    //other_cv_->knotIntervalFuzzy(par, fuzzy_other_);
	    
	    if (fabs(par - pmin) < fuzzy_other_) {
		if (pmin != lower_bound && cl_dist < prev_dist) {
		    // change search interval to other side of critical value
		    pmax = pmin;
		    pmin = std::max(other_cv_->nextSegmentVal(pmin, false, pareps), lower_bound);
		    repeat = true;
		}
	    } else if (fabs(par - pmax) < fuzzy_other_) {
		if (pmax != upper_bound && cl_dist < prev_dist) {
		    // change serach interval to other side of critical value
		    pmin = pmax;
		    pmax = std::min(other_cv_->nextSegmentVal(pmax, true, pareps), upper_bound);
		    repeat = true;
		}
	    }
	    prev_dist = cl_dist;
	} while (repeat);
    }
}

//===========================================================================
void BindToSurface::other_closest_point(const Point& pt,
					double* cl_par,
					double& cl_dist,
					Point& cl_pt,
					bool mind_segment_values)
//===========================================================================
{
    double& par_u = cl_par[0];
    double& par_v = cl_par[1];
    
    double epsilon = tol_->getEpsge() * 1e-2; // @@ since we are typically using this routine
                                              // to determine coincidence areas where the 
                                              // closeness of two objects is less than
                                              // geometrical tolerance, we need the tolerance
                                              // for the closest point algorithm to be a few
                                              // orders of magnitude smaller than that.

    // @@sbr We're currently assuming that the sf has a startparameter. What about trimmed sfs?
    const ParamSurface* par_sf = other_sf_->getParamSurface().get();
    RectDomain dom = other_sf_->getDomain();
    const RectDomain deg_dom = ((ParamSurfaceInt*)other_sf_)->getDegDomain(tol_->getEpsge());

    const Vector2D& lower_left = deg_dom.lowerLeft();
    const Vector2D& upper_right = deg_dom.upperRight();

    const double pareps_u = std::max(tol_->getRelParRes() * (upper_right[0] - lower_left[0]), 
				     tol_->getRelParRes());
    const double pareps_v = std::max(tol_->getRelParRes() * (upper_right[1] - lower_left[1]),
				     tol_->getRelParRes());

    // VSK, 0609 Make sure that the seed is inside the legal domain
    cl_par[0] = std::max(cl_par[0], lower_left[0]);
    cl_par[0] = std::min(cl_par[0], upper_right[0]);
    cl_par[1] = std::max(cl_par[1], lower_left[1]);
    cl_par[1] = std::min(cl_par[1], upper_right[1]);

    if (!mind_segment_values) {
	par_sf->closestPoint(pt,      // reference point
				par_u,   // u parameter
				par_v,   // v parameter
				cl_pt,   // pos. of closest pt.
				cl_dist, // distance to ref. pt.
				epsilon, // internal tolerance
				&deg_dom,     // domain to search
				cl_par); // seed
    } else {
	// we must consider possible problems near segment transitions
	double umin = std::max(lower_left[0], other_sf_->nextSegmentVal(0, par_u, false, pareps_u));
	double vmin = std::max(lower_left[1], other_sf_->nextSegmentVal(1, par_v, false, pareps_v));
	double umax = std::min(upper_right[0], other_sf_->nextSegmentVal(0, par_u, true, pareps_u));
	double vmax = std::min(upper_right[1], other_sf_->nextSegmentVal(1, par_v, true, pareps_v));
	const ParamSurface* ps = other_sf_->getParamSurface().get();

    const Vector2D& lower_left = dom.lowerLeft();
    const Vector2D& upper_right = dom.upperRight();

	bool repeat;
	double prev_dist = cl_dist = LARGEST_POSSIBLE_VALUE;
	do {
	    repeat = false;
	    RectDomain dom1(lower_left, upper_right);
	    ps->closestPoint(pt, par_u, par_v, cl_pt, cl_dist, epsilon, &dom1, cl_par);
	    if (fabs(par_u - umin) < fuzzy_u_) {
		if (par_u != lower_left[0] && cl_dist < prev_dist) {
		    // change search interval to other side of critical value
		    umax = umin;
		    umin = std::max(other_sf_->nextSegmentVal(0, par_u, false, pareps_u), 
			       lower_left[0]);
		    repeat = true;
		}
	    } else if (fabs(par_u - umax) < fuzzy_u_) {
		if (par_u != upper_right[0] && cl_dist < prev_dist) {
		    // change search interval to other side of critical value
		    umin = umax;
		    umax = std::min(other_sf_->nextSegmentVal(0, par_u, true, pareps_u), 
			       upper_right[0]);
		    repeat = true;
		}
	    }
	    if (fabs(par_v - vmin) < fuzzy_v_) {
		if (par_v != lower_left[1] && cl_dist < prev_dist) {
		    // change search interval to other side of critical value
		    vmax = vmin;
		    vmin = std::max(other_sf_->nextSegmentVal(1, par_v, false,pareps_v), 
			       lower_left[1]);
		    repeat = true;
		}
	    } else if (fabs(par_v - vmax) < fuzzy_v_) {
		if (par_v != upper_right[1] && cl_dist < prev_dist) {
		    // change search interval to other side of critical value
		    vmin = vmax;
		    vmax = std::min(other_sf_->nextSegmentVal(1, par_v, true, pareps_v), 
			       upper_right[1]);
		    repeat = true;
		}
	    }
	    prev_dist = cl_dist;
	} while (repeat);
    }
//     Point curr_pt;
//     double param[2];
//     param[0] = par_u;
//     param[1] = par_v;
//     other_sf_->point(curr_pt, param[0], param[1]);
}



// //===========================================================================
// void BindToCurve::determineConfidenceInterval(double& inside_marching_param,
// 					      double* inside_other_param,
// 					      double& outside_marching_param,
// 					      double* outside_other_param,
// 					      double& inside_dist,
// 					      double& outside_dist)
// //===========================================================================
// {
//     const double& center_value = tol_->getEpsge();
//     const double& bracket = tol_->getEpsBracket();
//     const double lower_limit = (1 - bracket) * center_value;
//     const double upper_limit = (1 + bracket) * center_value;

//     ASSERT(inside_dist < center_value && outside_dist > center_value);

//     bool outside_val_fixed = (outside_dist < upper_limit);
//     bool inside_val_fixed = (inside_dist > lower_limit);

//     //    bool forward_marching = outside_marching_param > inside_marching_param;
//     //    bool forward_other = outside_other_param[0] > inside_other_param[0];

//     Point p1, p2;
//     double new_dist;

//     // we know that 'center_value' is bracketed by 'inside_dist' and 'outside_dist'.
//     // Now we need to make the bracket narrow enough to fit inside [lower_limit, upper_limit]
//     const ParamCurve* march_curve = march_cv_; //->getParamCurve();
//     const ParamCurve* other_curve = other_cv_; //->getParamCurve();
//     while (!inside_val_fixed || !outside_val_fixed) {
// 	double mid_p_marching = (inside_marching_param + outside_marching_param) * 0.5;
// 	double mid_p_other = (inside_other_param[0] + outside_other_param[0]) * 0.5;
// 	double pmin = std::min(inside_other_param[0], outside_other_param[0]);
// 	double pmax = std::max(inside_other_param[0], outside_other_param[0]);
// 	march_curve->point(p1, mid_p_marching);
// 	other_curve->closestPoint(p1, pmin, pmax, mid_p_other, p2, new_dist, &mid_p_other);
// 	if (new_dist > center_value) {
// 	    outside_marching_param = mid_p_marching;
// 	    outside_other_param[0] = mid_p_other;
// 	    outside_dist = new_dist;
// 	    outside_val_fixed = (outside_dist <= upper_limit);
// 	} else { // new_dist <= center_value
// 	    inside_marching_param = mid_p_marching;
// 	    inside_other_param[0] = mid_p_other;
// 	    inside_dist = new_dist;
// 	    inside_val_fixed = (new_dist >= lower_limit);
// 	}
//     }
// }

//===========================================================================
void BindToPoint::expandCoincidenceRegion(double& march_par,
					  double* other_par,
					  double& dist,
					  bool& passed_tolerance)
//===========================================================================
{
    // calculating points
    vector<Point> march_pt(3), other_pt(3); // march point and other point
    const double pareps = std::max(tol_->getRelParRes(),
				   tol_->getRelParRes() * (march_cv_->endparam() - march_cv_->startparam()));
    march_cv_->point(march_pt, march_par, 2, step_forward_);
    // @@sbr The other object has a rather simple iso-curve.
    other_pt[0] = other_point_;
    other_pt[1] = Point(1);
    other_pt[1][0] = 1.0;
    other_pt[2] = Point(1);
    other_pt[2][0] = 1.0;
//     other_cv_->point(other_pt, *other_par, 2);

    // determining end of segment value (passing them in one step can be dangerous)
    double march_eos = march_cv_->nextSegmentVal(march_par, step_forward_, pareps);
    
    // calculate new step length
    double step;
    double aim = (1 + tol_->getEpsBracket()) * tol_->getEpsge();
    stepLength(march_pt, other_pt, aim, step_forward_, step);

    // VSK, 0405. Avoid too small step
    double len = march_pt[1].length();
    double delta = tol_->getEpsge()/len;
    double fac = 0.1;
    if (fabs(step) < fac*delta)
	step = step_forward_ ? fac*delta : -fac*delta;

    if (hasPassedVal(march_par + step, march_eos, step_forward_)) {
	step = march_eos - march_par;
    }
    
    // march along curves
    double tolerance = tol_->getEpsge();
    dist = 0;

    // checking midpt. and endpt. of interval
    for (int i = 0; i < 2 && dist < tolerance; ++i){
	march_par += 0.5 * step;
	march_cv_->point(march_pt[0], march_par);
	other_closest_point(march_pt[0], other_par, dist, other_pt[0], true);
    }
    passed_tolerance = (dist > tolerance);
    //march_cv_->knotIntervalFuzzy(march_par, fuzzy_march_);
}

//===========================================================================
void BindToCurve::expandCoincidenceRegion(double& march_par,
					  double* other_par,
					  double& dist,
					  bool& passed_tolerance)
//===========================================================================
{
    // calculating points
    vector<Point> march_pt(3), other_pt(3); // march point and other point
    const double pareps = std::max(tol_->getRelParRes(),
				   tol_->getRelParRes() * (march_cv_->endparam() - march_cv_->startparam()));
    march_cv_->point(march_pt, march_par, 2, step_forward_);
    other_cv_->point(other_pt, *other_par, 2);

    // determining end of segment value (passing them in one step can be dangerous)
    double march_eos = march_cv_->nextSegmentVal(march_par, step_forward_, pareps);
    
    // calculate new step length
    double step;
    double aim = (1 + tol_->getEpsBracket()) * tol_->getEpsge();
    ASSERT(march_pt[0].dist(other_pt[0]) < aim);
    stepLength(march_pt, other_pt, aim, step_forward_, step);

    // VSK, 0405. Avoid too small step
    double len = march_pt[1].length();
    double delta = tol_->getEpsge()/len;
    double fac = 0.1;
    if (fabs(step) < fac*delta)
	step = step_forward_ ? fac*delta : -fac*delta;

    if (hasPassedVal(march_par + step, march_eos, step_forward_)) {
	step = march_eos - march_par;
    }
    
    // march along curves
    double tolerance = tol_->getEpsge();
    dist = 0;

    // checking midpt. and endpt. of interval
    for (int i = 0; i < 2 && dist < tolerance; ++i){
	march_par += 0.5 * step;
	march_cv_->point(march_pt[0], march_par);
	other_closest_point(march_pt[0], other_par, dist, other_pt[0], true);
    }
    passed_tolerance = (dist > tolerance);
    //march_cv_->knotIntervalFuzzy(march_par, fuzzy_march_);
}

//===========================================================================
void BindToSurface::expandCoincidenceRegion(double& march_par,
					    double* other_par,
					    double& dist,
					    bool& passed_tolerance)
//===========================================================================
{
    // calculating points
    double geotol = tol_->getEpsge();
    vector<Point> march_pt(3), surf_pt(6);
    march_cv_->point(march_pt, march_par, 2, step_forward_);
    other_sf_->point(surf_pt, other_par, 2);

	dist = march_pt[0].dist(surf_pt[0]);
    ASSERT(dist < geotol);

    const double partol = 
	tol_->getRelParRes() * (march_cv_->endparam() - march_cv_->startparam());
 
    // determining end of segment value
    double march_eos = march_cv_->nextSegmentVal(march_par, step_forward_, partol);
    if (fabs(march_eos - march_par) <= tol_->getRelParRes())
	march_eos = march_cv_->nextSegmentVal(march_eos, step_forward_, partol);
	if (fabs(march_par - march_eos) < tol_->getRelParRes())
	{
		passed_tolerance = false;
		return;
	}


    // calculate new step length
    double aim = (1 + tol_->getEpsBracket()) * geotol;
    vector<Point> tmp; // dummy
    double min_val = (step_forward_) ? 0 : march_eos - march_par; 
    double max_val = (step_forward_) ? march_eos - march_par : 0;

    double step = stepLength(march_pt, surf_pt, step_forward_, tmp, min_val, max_val, aim);

    // VSK, 0405. Avoid too small step
    double len = march_pt[1].length();
    double delta = tol_->getEpsge()/len;
    double fac = 0.1;
    if (fabs(step) < fac*delta)
	step = step_forward_ ? fac*delta : -fac*delta;

    if (hasPassedVal(march_par + step, march_eos, step_forward_)) {
	step = march_eos - march_par;
    }

    // checking midpoint and endpoint of interval
    dist = 0;
    for (int i = 0; i < 2 && dist < geotol; ++i) {
	march_par += 0.5 * step;
	march_cv_->point(march_pt[0], march_par);
	other_closest_point(march_pt[0], other_par, dist, surf_pt[0], true);
    }
    passed_tolerance = (dist > geotol);
    //march_cv_->knotIntervalFuzzy(march_par, fuzzy_march_);
}


// //===========================================================================
// void BindToSurface::determineConfidenceInterval(double& last_pval_inside_curve,
// 						double* last_pval_inside_other,
// 						double& first_pval_outside_curve,
// 						double* first_pval_outside_other,
// 						double& inside_dist,
// 						double& outside_dist)
// //===========================================================================
// {
//     const double& center_value = tol_->getEpsge();
//     const double& bracket = tol_->getEpsBracket();
//     const double lower_limit = (1 - bracket) * center_value;
//     const double upper_limit = (1 + bracket) * center_value;
//     ASSERT(inside_dist < center_value && outside_dist > center_value);

//     bool outside_val_fixed;
// }


}; // end anonymous namespace 

// //===========================================================================
// void determine_confidence_interval(const ParamCurveInt* const marching_curve,
// 				   const ParamCurveInt* const other_curve,
// 				   double& inside_param_marching,
// 				   double& inside_param_other,
// 				   double& outside_param_marching,
// 				   double& outside_param_other,
// 				   double& inside_val,
// 				   double& outside_val, 
// 				   double center_value,
// 				   double bracket,
// 				   bool& succeeded)
// //===========================================================================
// {
//     double lower_limit = (1 - bracket) * center_value;
//     double upper_limit = (1 + bracket) * center_value;

//     bool outside_val_fixed = false;
//     bool inside_val_fixed = false;

//     bool forward_marching = outside_param_marching > inside_param_marching;
//     bool forward_other = outside_param_other > inside_param_other;
    
//     Point temp_point_1, temp_point_2;
//     double temp_val;

//     // make sure center_value is bracketed
//     if (outside_val < center_value) {
// 	double step_marching = outside_param_marching - inside_param_marching;
// 	double step_other = outside_param_other - inside_param_other;
// 	bool bracketed_center_value = false;
// 	double pmin, pmax;
// 	double endval_marching = forward_marching ? marching_curve->endparam() :
// 	                                            marching_curve->startparam();
// 	double endval_other = forward_other ? other_curve->endparam() :
// 	                                      other_curve->startparam();
// 	bool endval_reached = false;
// 	while (!bracketed_center_value) { 
// 	    // @ nb: marching in this loop does not consider critical values.
// 	    // Maybe this might cause problems in some obscure cases???
// 	    inside_val = outside_val;
// 	    inside_param_marching = outside_param_marching;
// 	    inside_param_other = outside_param_other;
// 	    outside_param_marching += step_marching;
// 	    outside_param_other += step_other;
// 	    if (hasPassedVal(outside_param_marching, endval_marching, forward_marching)) {
// 		outside_param_marching = endval_marching;
// 		endval_reached = true;
// 	    }
// 	    if (hasPassedVal(outside_param_other, endval_other, forward_other)) {
// 		outside_param_other = endval_other;
// 	    }

// 	    pmin = std::min(inside_param_other, endval_other);
// 	    pmax = std::max(inside_param_other, endval_other);
// 	    marching_curve->getParamCurve()->point(temp_point_1, outside_param_marching);
// 	    other_curve->getParamCurve()->closestPoint(temp_point_1,
// 						       pmin, pmax,
// 						       outside_param_other,
// 						       temp_point_2,
// 						       outside_val,
// 						       &outside_param_other);
// 	    bracketed_center_value = (outside_val > center_value);
// 	    if (endval_reached && !bracketed_center_value) {
// 		// unsuccessful in bracketing the center value, because we have reached
// 		// the end of the curve
// 		succeeded = false;
// 		return;
// 	    }
// 	}
// 	inside_val_fixed = true;
//     } 
//     // when we got here, we know that 'center_value' is bracketed by 
//     // 'inside_val' and 'outside_val'.  Now we need to make the bracket narrow
//     // fit inside [lower_limit, upper_limit]
//     outside_val_fixed = (outside_val <= upper_limit);

//     // making sure 'inside_val' is bigger than 'lower_limit' and 'outside_val'
//     // is smaller than 'upper_limit', using bisection
//     // @@ It can be recommended to change bisection to the 'False Position Method'
//     // (detailed in Numerical Recipes 9.2) to speed up the algorithm.
//     while (!inside_val_fixed || !outside_val_fixed) {
// 	double mid_p_marching = (inside_param_marching + outside_param_marching) * 0.5;
// 	double mid_p_other    = (inside_param_other + outside_param_other) * 0.5;
// 	double pmin = std::min(inside_param_other, outside_param_other);
// 	double pmax = std::max(inside_param_other, outside_param_other);
// 	marching_curve->getParamCurve()->point(temp_point_1, mid_p_marching);
// 	other_curve->getParamCurve()->closestPoint(temp_point_1,
// 						   pmin, pmax,
// 						   mid_p_other,
// 						   temp_point_2,
// 						   temp_val, 
// 						   &mid_p_other);
// 	if (temp_val > center_value) {
// 	    outside_param_marching = mid_p_marching;
// 	    outside_param_other = mid_p_other;
// 	    outside_val = temp_val;
// 	    outside_val_fixed = (outside_val <= upper_limit);
// 	} else { // temp_val <= center_value
// 	    inside_param_marching = mid_p_marching;
// 	    inside_param_other = mid_p_other;
// 	    inside_val = temp_val;
// 	    inside_val_fixed = (temp_val >= lower_limit);
// 	}
//     }
// }



// //===========================================================================
// void BindToSurface::expandCoincidenceRegion(double& marching_par,
// 					    double* other_par,
// 					    double& dist,
// 					    bool& passed_lower_limit)
// //===========================================================================
// {
//     // calculating points
//     double& upar = other_par[0];
//     double& vpar = other_par[1];
//     vector<Point> marching_point(3), other_point(6);
//     marching_cv_->point(marching_point, marching_par, 2, step_forward_);
//     other_sf_->point(other_point, upar, vpar, 2);

//     // determining next critical values (passing them in one step can be dangerous)
//     double march_critical = marching_cv_->nextSegmentVal(marching_par, step_forward_);
    
//     // calculate new step length
//     static vector<Point> scratch(3);
//     double min_step = 0;
//     double max_step = march_critical - marching_par;
//     double step = stepLength(marching_point, 
// 			     other_point, 
// 			     step_forward_,
// 			     scratch,
// 			     min_step,
// 			     max_step,
// 			     tol_->aepsge());
//     // the following assertion should always be true, due to the max_step set above
//     ASSERT(!hasPassedVal(marching_par + step, march_critical, step_forward_));

//     RectDomain dom; // domain of interest
//     dom.umin() = ;
//     dom.umax() = ;
//     dom.vmin() = ;
//     dom.vmax() = ;

// //     // determine search interval for other point
// //     double pmin = other_forward ? *other_par : other_critical;
// //     double pmax = other_forward ? other_critical : *other_par;
    
//     // march along curves
//     const RectDomain& whole_domain = other_sf_->getDomain();
//     double fuzzy_u = fuzzy_eps_ * (whole_domain.umax() - whole_domain.umin());
//     double fuzzy_v = fuzzy_eps_ * (whole_Domain.vmax() - whole_domain.vmin());
//     bool exited = false;
//     passed_lower_limit = false;
//     for (int i = 0; i < 2; ++i) { // checking midpoint and endpoint of interval
// 	marching_par += 0.5 * step;
// 	marching_cv_->getParamCurve()->point(marching_point[0], marching_par);
// 	other_sf_->getParamSurface()->closestPoint(marching_point[0],
// 						   upar, vpar, // par. of closest pt.
// 						   other_point[0], // closest point
// 						   dist, // distance to cl. pt.
// 						   tol_->aepsge(), // geom. tolerance
// 						   &dom, // domain of interest
// 						   other_par); // seed;
// 	other_sf_->knotIntervalFuzzy(upar, vpar, fuzzy_u, fuzzy_v);
// 	if (dist > lower_limit) {
// 	    passed_lower_limit = true;
// 	    break;
// 	}
//     }
//     double fuzzy = fuzzy_eps_ * (marching_cv_->endparam() - marching_cv_->startparam());
//     marching_cv_->knotIntervalFuzzy(marching_par, fuzzy);
// }

// //===========================================================================
// bool measureCoincidenceRegion(const ParamCurveInt* marching_curve,
// 			      const ParamCurveInt* other_curve,
// 			      double isect_param_marching,
// 			      double isect_param_other,
// 			      bool step_forward,
// 			      const shared_ptr<const GeoTol> gtol,
// 			      double& last_param_inside_marching,
// 			      double& first_param_outside_marching,
// 			      double& last_param_inside_other,
// 			      double& first_param_outside_other)
// //===========================================================================
// {
//     double eps = gtol->getEpsge();
//     double eps_bracket = gtol->getEpsBracket();

//     // argument checks
//     DEBUG_ERROR_IF(!marching_curve || !other_curve, "Null pointers found.");
//     DEBUG_ERROR_IF(eps_bracket <= 0 || eps_bracket >= 1, "wrong value for eps_bracket.");
    
//     // initializing
//     double lower_limit = (1 - eps_bracket) * eps;
//     bool passed_lower_limit = false;
//     double prev_distance, current_distance = 0;
//     first_param_outside_marching = isect_param_marching;
//     first_param_outside_other = isect_param_other;
//     double end_val = step_forward ? marching_curve->endparam() : 
// 	                            marching_curve->startparam();

    
//     // iterating until lower bracket is passed, or end of marching curve reached
//     while (!passed_lower_limit) {
// 	// saving old parameter values
// 	last_param_inside_marching = first_param_outside_marching;
// 	last_param_inside_other = first_param_outside_other;
// 	prev_distance = current_distance;
// 	expand_coincidence_region(marching_curve, 
// 				  other_curve,
// 				  step_forward, 
// 				  first_param_outside_marching,
// 				  first_param_outside_other,
// 				  current_distance,
// 				  lower_limit,
// 				  eps,
// 				  passed_lower_limit);
// 	if (first_param_outside_marching == end_val) {
// 	    // arrived at boundary
// 	    break;
// 	}
//     }
//     if (!passed_lower_limit) {
// 	// we arrived at end of marching curve without passing the lower limit
// 	return false;
//     }
//     bool succeeded;
//     determine_confidence_interval(marching_curve, other_curve,
// 				  last_param_inside_marching, last_param_inside_other,
// 				  first_param_outside_marching, first_param_outside_other,
// 				  prev_distance, current_distance,
// 				  eps, eps_bracket, succeeded);
//     return succeeded;
// }



// //===========================================================================
// void expand_coincidence_region(const ParamCurveInt* const marching_curve,
// 			       const ParamCurveInt* const secondry_curve,
// 			       bool marching_fwd,
// 			       double& marching_par,
// 			       double& secondry_par,
// 			       double& dist,
// 			       double lower_limit, // smaller than eps
// 			       double eps,
// 			       bool& passed_lower_limit)
// //===========================================================================
// {
//     // calculating points
//     vector<Point> marching_point(3), secondry_point(3);
//     marching_curve->point(marching_point, marching_par, 2, marching_fwd);
//     secondry_curve->point(secondry_point, secondry_par, 2);
//     bool secondry_fwd = 
// 	marching_point[1] * secondry_point[1] >= 0 ? marching_fwd : !marching_fwd;
//     if (secondry_fwd == false) {    
// 	// re-calculating secondary point
// 	secondry_curve->point(secondry_point, secondry_par, 2, secondry_fwd);
//     }

//     // determining next critical values (passing them in one step can be dangerous)
//     double march_critical = marching_curve->nextSegmentVal(marching_par, marching_fwd);
//     double secnd_critical = secondry_curve->nextSegmentVal(secondry_par, secondry_fwd);

//     // calculate new step length
//     double step;
//     stepLength(marching_point, secondry_point, eps, marching_fwd, step);
//     if (hasPassedVal(marching_par + step, march_critical, marching_fwd)){
// 	step = march_critical - marching_par;
//     }
//     // determine search interval for secondary point
//     double pmin = secondry_fwd ? secondry_par : secnd_critical;
//     double pmax = secondry_fwd ? secnd_critical : secondry_par;

//     // march along curves
//     const double fuzzy_eps = 1.0e-4; // @@ move this elsewhere?  globally?
//     double fuzzy = fuzzy_eps*(secondry_curve->endparam() - secondry_curve->startparam());
//     bool exited = false;
//     passed_lower_limit = false;
//     for (int i = 0; i < 2; ++i) { // checking midpoint and endpoint of interval
// 	marching_par += 0.5 * step;
// 	marching_curve->getParamCurve()->point(marching_point[0], marching_par);
// 	secondry_curve->getParamCurve()->closestPoint(marching_point[0], 
// 						      pmin, pmax, // search interval
// 						      secondry_par, // new param. value
// 						      secondry_point[0],  // closest point
// 						      dist,         // closest distance
// 						      &secondry_par); // seed
// 	secondry_curve->knotIntervalFuzzy(secondry_par, fuzzy);
// 	if (dist > lower_limit) {
// 	    passed_lower_limit = true;
// 	    break;
// 	} 
//     }
//     marching_curve->knotIntervalFuzzy(marching_par, fuzzy);
// }

// //===========================================================================
// bool expand_coincidence_region(const ParamCurveInt* const marching_curve,
// 			       const ParamCurveInt* const secondry_curve,
// 			       bool marching_fwd,
// 			       bool secondry_fwd,
// 			       double& marching_par,
// 			       double& secondry_par,
// 			       double& dist,
// 			       bool& change_marching_curve,
// 			       double lower_limit, // smaller than eps
// 			       double eps) 
// //===========================================================================
// {
//     // calculating points
//     vector<Point> marching_point(3), secondry_point(3);
//     marching_curve->point(marching_point, marching_par, 2);
//     secondry_curve->point(secondry_point, secondry_par, 2);

//     // determining next critical values
//     double marching_crtcl = marching_curve->nextSegmentVal(marching_par, marching_fwd);
//     double secondry_crtcl = secondry_curve->nextSegmentVal(secondry_par, secondry_fwd);

//     // calculate new step length
//     double step;
//     stepLength(marching_point, secondry_point, eps, marching_fwd, step);
//     if (hasPassedVal(marching_par + step, marching_crtcl, marching_fwd)){
// 	step = marching_crtcl - marching_par;
//     }
//     // determine search interval for secondary point
//     double pmin = secondry_fwd ? secondry_par : secondry_crtcl;
//     double pmax = secondry_fwd ? secondry_crtcl : secondry_par;

//     // march along curves
//     const double fuzzy_eps = 1.0e-4; // @@ move this elsewhere?  globally?
//     double fuzzy = fuzzy_eps*(secondry_curve->endparam()-secondry_curve->startparam());
//     change_marching_curve = false;
//     bool exited = false;
//     for (int i = 0; i < 2; ++i) { // checking midpoint and endpoint of interval
// 	marching_par += 0.5 * step;
// 	marching_curve->getParamCurve()->point(marching_point[0], marching_par);
// 	secondry_curve->getParamCurve()->closestPoint(marching_point[0], 
// 						      pmin, pmax, // search interval
// 						      secondry_par, // new param. value
// 						      secondry_point[0],  // closest point
// 						      dist,         // closest distance
// 						      &secondry_par); // seed
// 	if (hasPassedVal(secondry_par, secondry_crtcl, secondry_fwd)) {
// 	    secondry_par = secondry_crtcl;
// 	    change_marching_curve = true;
// 	}
// 	secondry_curve->knotIntervalFuzzy(secondry_par, fuzzy);
// 	if (dist > lower_limit) {
// 	    return true;
// 	} 
//     }
//     return false;
// }

// //===========================================================================
// void measureCoincidenceRegion(ParamCurveInt* curve1,
// 			      ParamCurveInt* curve2,
// 			      double isect_param_1,
// 			      double isect_param_2,
// 			      bool step_forward_1,
// 			      bool step_forward_2,
// 			      double eps,
// 			      double eps_bracket,
// 			      double& last_par_inside_1,
// 			      double& first_par_outside_1,
// 			      double& last_par_inside_2,
// 			      double& first_par_outside_2)
// //===========================================================================
// {
//     // argument checks
//     DEBUG_ERROR_IF(!curve1 || !curve2, "Null pointers detected in measureCoincidenceRegion()");
//     DEBUG_ERROR_IF(eps_bracket <= 0 || eps_bracket >= 1, 
// 	     "eps_bracket must be between 0 and 1 in measureCoincidenceRegion()");

//     // initializing
//     bool march_along_first_curve = true;
//     bool change_marching_curve = false;
//     bool passed_lower_limit = false;
//     double lower_limit = (1 - eps_bracket) * eps; 
//     double current_distance, prev_distance;
//     double& current_param_1 = first_par_outside_1;
//     double& current_param_2 = first_par_outside_2;
//     double& prev_param_1 = last_par_inside_1;
//     double& prev_param_2 = last_par_inside_2;
//     current_param_1 = isect_param_1;
//     current_param_2 = isect_param_2;

//     // iterating until lower bracket is passed
//     while(!passed_lower_limit) {
// 	// saving old parameter values
// 	prev_param_1 = current_param_1; 
// 	prev_param_2 = current_param_2;
// 	prev_distance = current_distance;
// 	if (march_along_first_curve) {
// 	    passed_lower_limit = 
// 		expand_coincidence_region(curve1, 
// 					  curve2,
// 					  step_forward_1, 
// 					  step_forward_2,
// 					  current_param_1, 
// 					  current_param_2,
// 					  current_distance,
// 					  change_marching_curve,
// 					  lower_limit, 
// 					  eps);
// 	} else {
// 	    passed_lower_limit = 
// 		expand_coincidence_region(curve2, curve1,
// 					  step_forward_2, step_forward_1,
// 					  current_param_2, current_param_1,
// 					  current_distance,
// 					  change_marching_curve,
// 					  lower_limit,
// 					  eps);
// 	}
//     }
//     // we have now found a point with a distance bigger than (1 - eps_bracket) * eps
//     determine_confidence_interval(curve1, curve2,
// 				  prev_param_1, prev_param_2,
// 				  current_param_1, current_param_2,
// 				  prev_distance,
// 				  current_distance,
// 				  eps, eps_bracket);
// }

// //===========================================================================
// int stepLength(const vector<Point>& ft, const vector<Point>& gs, 
// 	       double delta, bool forward, double& delta_t)
// //===========================================================================
// {

//   /* Given two curves f(t), g(s) and parameters t1 and s1 to an intersection
//    * point, this function computes the length of the next parameter step along 
//    * the curve f(t).
//    * The curves must have been evaluated in the start point (ft and gs).
//    * @verbatim
//    Assume f(t1) = g(s1).

//    Let h(t) = [f(t) - g(s(t))]^2.
//    g(s(t1)) is the point on g closest to f(t1).
   
//    We want to find delta_t so that 
//       h(t1 + delta_t) = delta = 0.9*eps^2

//    Series expansion gives :
//      delta = h(t1+delta_t) = h(t1) + delta_t*h'(t1) + delta_t^2/2*h"(t1)
//                               =0

//       h"/2*delta_t^2 + h'*delta_t - delta = 0

//    Where 
//       h'(t) = 2[f(t)-g(s(t))][f'(t)-g'(s)*s'(t)]
//       h"(t) = 2[f'(t)-g'(s)*s'(t)]^2 + 
//               2[f(t)-g(s(t))][f"(t)-g"(s)*s'(t)^2 - g'(s)*s"(t)]

//       s'(t) = (f'(t).g'(s))/(g'(s).g'(s))
//       s"(t) = (f"(t).g'(s) + f'(t).g(s)*s'(t) - 2g'(s).'g"(s)*s'(t)^2) /
//               ((g'(s).g'(s))
//    @endverbatim
//    * \param forward True if we are stepping forward along the curve.
//    * \param delta_t The calculated step length
//    * \return 0:Error.  >0 OK.
//    */


//   const double rad_eps = 1.e-4; //@@ ??  // Absolute tolerance describing the 
//   // deviation between the circle of curvature and an Hermite approximation
//   // to the circle.

//   double fdgd = ft[1]*gs[1];   // f'(t) . g'(s)
//   double gdgd = gs[1]*gs[1];   // g'(s) . g'(s)

//   if (gdgd == 0.0) { // Should not happen!
//     std::cout << "\nWARNING! g'(s) = 0" << std::endl;
//     return 0;
//   }


//   double sd1 = fdgd/gdgd; // s'(t)
//   double sd2 =
//     (ft[2]*gs[1]+(fdgd-2.0*(gs[1]*gs[2])*sd1)*sd1)/(gdgd*gdgd);  // s"(t)

//   Point p0 = ft[0]-gs[0];
//   Point p1 = ft[1]-gs[1]*sd1;

//   double hd1 = 2.0*(p0*p1); // h'(t)
//   double hd2 = 2.0*(p1*p1 + p0*(ft[2]-gs[2]*(sd1*sd1)-gs[1]*sd2)); // h"(t)

//   std::cout << "fdgd, gdgd : " << fdgd <<" "<< gdgd << endl; 
//   std::cout << "sd1, sd2   : " << sd1  <<" "<< sd2 << endl; 
//   std::cout << "hd1, hd2   : " << hd1  <<" "<< hd2 << endl; 


//   //  int nsolu = secDegreeEq(0.5*hd2, hd1,-delta, delta_t1, delta_t2);

//   if (hd2==0.0) {
 
//    if (hd1==0.0) {    // Equal values, first and second derivatives.
//                        // Use radius of curvature.
//        vector<Point> derivs;
//        //double crad = curvatureRadius(ft,gs);
//        double crad = curvatureRadius(ft,derivs);
//       delta_t = stepLenFromRadius(crad, rad_eps);
//       if (!forward)
// 	delta_t = -delta_t;
//       cout << "Step  1 curvatureRadius :" <<  delta_t << endl;
//       return 1;
//     }
//     else {
//       delta_t = delta/hd1;
//       if ((forward && delta_t > 0.0) || (!forward && delta_t < 0.0)) {
// 	cout << "Step  1 :" <<  delta_t << endl;
// 	return 1;
//       }
//       else {
// 	cout << "\nStep  0 :" <<  delta_t << endl;
// 	return 0;
//       }
//     }
//   }
//   else {
//     double d = hd1*hd1+2.0*delta*hd2;

//     if (d < 0.0) {  // No real roots. Use radius of curvature.
//       cout << "\nNo real roots. discr.:" <<  d << endl;
//       vector<Point> derivs;
//       //double crad = curvatureRadius(ft,gs);
//       double crad = curvatureRadius(ft,derivs);
//       delta_t = stepLenFromRadius(crad, rad_eps);
//       if (!forward)
// 	delta_t = -delta_t;
//       cout << "Step  1 curvatureRadius :" <<  delta_t << endl;
//       return 1;
//     }

//     d = sqrt(d);
//     double delta_t1 = -hd1/hd2 + d/hd2;
//     double delta_t2 = -hd1/hd2 - d/hd2;

//     if (forward)
//       delta_t = std::max(delta_t1,delta_t2); //@@@ test om pos. og neg. delta_t?
//     else
//       delta_t = std::min(delta_t1,delta_t2);
//     cout << "Step  2 :" << delta_t1 <<" "<<  delta_t2 << endl;
//     return 2;
//   }

// }
// //===========================================================================
// void measureCoincidenceRegion(ParamCurveInt* curve1,
// 			      ParamCurveInt* curve2,
// 			      double isect_param_1,
// 			      double isect_param_2,
// 			      bool step_forward_1,
// 			      bool step_forward_2,
// 			      double eps,
// 			      double& last_par_inside_1,
// 			      double& first_par_outside_1,
// 			      double& last_par_inside_2,
// 			      double& first_par_outside_2)
// //===========================================================================
// {
//     // argument checks
//     DEBUG_ERROR_IF(!curve1 || !curve2, "Null pointers detected in measureCoincidenceRegion()");

//     bool march_along_first_curve = true;
//     bool change_marching_curve = false;
//     bool exited_region = false;
//     last_par_inside_1 = first_par_outside_1 = isect_param_1;
//     last_par_inside_2 = first_par_outside_2 = isect_param_2;

//     // loop while the points are within proximity of each other
//     while (!exited_region) {
// 	if (march_along_first_curve) {
// 	    exited_region = 
// 		expand_coincidence_region(curve1, curve2,
// 					  step_forward_1, step_forward_2,
// 					  last_par_inside_1, last_par_inside_2,
// 					  first_par_outside_1, first_par_outside_2,
// 					  change_marching_curve,
// 					  eps);
// 	    march_along_first_curve = !change_marching_curve;

// 	} else { // march along second curve
// 	    exited_region = 
// 		expand_coincidence_region(curve2, curve1,
// 					  step_forward_2, step_forward_1,
// 					  last_par_inside_2, last_par_inside_1,
// 					  first_par_outside_2, first_par_outside_1, 
// 					  change_marching_curve,
// 					  eps);
// 	    march_along_first_curve = change_marching_curve;
// 	}
//     }
// }

// //===========================================================================
// void measureCoincidenceRegion(ParamCurveInt* curve1,
// 			      ParamCurveInt* curve2,
// 			      double isect_param_1,
// 			      double isect_param_2,
// 			      bool step_forward_1,
// 			      bool step_forward_2,
// 			      double eps,
// 			      double eps_bracket,
// 			      double& last_par_inside_1,
// 			      double& first_par_outside_1,
// 			      double& last_par_inside_2,
// 			      double& first_par_outside_2)
// //===========================================================================
// {
//     // argument checks
//     DEBUG_ERROR_IF(!curve1 || !curve2, "Null pointers detected in measureCoincidenceRegion()");
//     DEBUG_ERROR_IF(eps_bracket <= 0 || eps_bracket >= 1, 
// 	     "eps_bracket must be between 0 and 1 in measureCoincidenceRegion()");

//     // initializing
//     double lower_limit = (1 - eps_bracket) * eps; 
//     bool passed_lower_limit = false;
//     bool arrived_at_boundary = false;
//     double current_distance, prev_distance;
//     double& current_param_1 = first_par_outside_1;
//     double& current_param_2 = first_par_outside_2;
//     double& prev_param_1 = last_par_inside_1;
//     double& prev_param_2 = last_par_inside_2;
//     current_param_1 = isect_param_1;
//     current_param_2 = isect_param_2;

//     // iterating until lower bracket is passed
//     while(!passed_lower_limit) {
// 	// saving old parameter values
// 	prev_param_1 = current_param_1; 
// 	prev_param_2 = current_param_2;
// 	prev_distance = current_distance;
// 	passed_lower_limit = 
// 	    expand_coincidence_region(curve1, 
// 				      curve2,
// 				      step_forward_1, 
// 				      step_forward_2,
// 				      current_param_1, 
// 				      current_param_2,
// 				      current_distance,
// 				      change_marching_curve,
// 				      lower_limit, 
// 				      eps,
// 				      arrived_at_boundary);
// 	if (arrived_at_boundary) {
// 	    // we cannot progress anymore in this direction
	    
// 	}
//     }
//     // we have now found a point with a distance bigger than (1 - eps_bracket) * eps
//     determine_confidence_interval(curve1, curve2,
// 				  prev_param_1, prev_param_2,
// 				  current_param_1, current_param_2,
// 				  prev_distance,
// 				  current_distance,
// 				  eps, eps_bracket);
// }


// //===========================================================================
// double stepLength(vector<Point>& ft, vector<Point>& gs, 
// 		  bool forward, std::vector<Point>& cvder,
// 		  double min_step, double max_step, 
// 		  double aepsge)
// //===========================================================================
// { 

//   // Evaluate shortest distance curve.
//   int kstat = evalDistCurve(ft, gs, cvder);
//   ALWAYS_ERROR_IF(kstat<0,"Zero determinant in evalDistCurve",ComputationError());

//   // Calculate the step length.
//   double tg0 = cvder[0]*cvder[0];
//   double tg1 = 2.0*cvder[0]*cvder[1];
//   double tg2 = 2.0*(cvder[1]*cvder[1] + cvder[0]*cvder[2]);

//   double tdel;
//   if (tg2 != 0.0) {
//     double td1, td2, td3;
//     td1 = 2.0*tg2*(tg0-aepsge*aepsge);
//     if (td1 > tg1*tg1)
// 	td1 *= -1.0;
//     td1 = sqrt(tg1*tg1 - td1);
//     td2 = (-tg1 + td1)/tg2;
//     td3 = (-tg1 - td1)/tg2;
//     tdel = 0.75*std::min(fabs(td2), fabs(td3));
//   }
//   else if (tg1 != 0.0)
//     tdel = 0.75*fabs((aepsge*aepsge - tg0)/tg1);
//   else
//     tdel = 0.1*max_step;
  
//   tdel = std::max(min_step, std::min(max_step, tdel));
  
//   if (forward)
//     return tdel;
//   else
//     return -tdel;
// }

// //===========================================================================
// int evalDistCurve(vector<Point>& ft, vector<Point>& gs, vector<Point>& res)
// //===========================================================================
// { 
//   // Evaluate the curve representing the shortest distance between a curve
//   // and a surface in the current point.
//   // Ported from sh1767_evdistcrv

//   // Evaluate distance vector
//   res[0] = ft[0] - gs[0];

//   // Make references to derivatives
//   Point& gder = ft[1];  Point& gder2 = ft[2];
//   Point& gu   = gs[1];  Point& gv    = gs[2];
//   Point& guu  = gs[3];  Point& guv   = gs[4];  Point& gvv = gs[5];  
     
//   // The determinant
//   double det = (gu*gu)*(gv*gv) - (gu*gv)*(gu*gv);
//   if (det == 0.0)
//     return -1;
     
//   //Use Cramer's rule to compute the 2x2 linear system, where the solution
//   //is the derivative of the surface in the parameter plane.
 
//   // First derivatives     
//   double u1 = ( (gv*gv)*(gder*gu) - (gu*gv)*(gder*gv) )/det;
//   double v1 = (-(gv*gu)*(gder*gu) + (gu*gu)*(gder*gv) )/det;
//   res[1] = gder - u1*gu - v1*gv;

//   // Second derivatives
//   Point vec = gder2 - u1*u1*guu - 2*u1*v1*guv - v1*v1*gvv;
//   double u2 = ( (gv*gv)*(vec*gu) - (gu*gv)*(vec*gv) )/det;
//   double v2 = (-(gv*gu)*(vec*gu) + (gu*gu)*(vec*gv) )/det;
//   res[2] = vec - u2*gu - v2*gv;
  
//   return 0;
// }

// //===========================================================================
// bool measure_coincidence_region(const ParamCurveInt* marching_curve,
// 				const ParamObjectInt* second_object,
// 				double isect_param_curve,
// 				const double* isect_param_other,
// 				bool step_forward_marching,
// 				const shared_ptr<const GeoTol> gtol,
// 				double& last_pval_inside_curve,
// 				double& first_pval_outside_curve)
// //===========================================================================
// {
//     double eps = gtol->getEpsge();
//     double eps_bracket = gtol->getEpsBracket();
    
//     // argument checks
//     DEBUG_ERROR_IF(!marching_curve || !second_object, "Null pointers detected.");
//     DEBUG_ERROR_IF(eps_bracket <= 0 || eps_bracket >= 1, "wrong value for eps_bracket.");

//     // initializing
//     double lower_limit = (1 - eps_bracket) * eps;
//     bool passed_lower_limit = false;
//     double prev_distance, current_distance = 0;
//     first_pval_outside_curve = isect_param_curve;
//     copy_n(isect_param_other, first_param_outside_other, second_object->numParams());
//     double end_val = step_forward ? marching_curve->endparam() : 
// 	                            marching_curve->startparam();
    
//     // choosing function according to type of second object
//     shared_ptr<BindToObject> bso;
//     ParamSurfaceInt* surf_2nd = dynamic_cast<ParamSurfaceInt*>(second_object);
//     ParamCurveInt* curve_2nd = surf_2nd ? 0 : dynamic_cast<ParamCurveInt*>(second_object);
//     if (curve_2nd) {
// 	bso = shared_ptr<BindToCurve>
// 	    (new BindToCurve(marching_curve, curve_2nd, step_forward, lower_limit, gtol));
//     } else if (surf_2nd) {
// 	bso = shared_ptr<BindToSurface>
// 	    (new BindToSurface(marching_curve, surf_2nd, step_forward, lower_limit, gtol));
//     } else { 
// 	THROW("Unrecognized/unsupported ParamObjectInt type.");
//     }

//     // iterating until lower bracket is passed, or end of marching curve reached
//     while(!passed_lower_limit) {
// 	// saving old parameter values
// 	last_pval_inside_curve = first_pval_outside_other;
// 	prev_distance = current_distance;
// 	bso->expandCoincidenceRegion(first_pval_outside_curve,
// 				     first_pval_outside_other,
// 				     current_distance,
// 				     passed_lower_limit);
// 	if (first_pval_outside_curve == end_val) {
// 	    break; // arrived at boundary
// 	}
//     }
//     if (!passed_lower_limit) {
// 	// we arrived at end of marching curve without passing the lower limit
// 	return false;
//     }
//     bool succeeded = bso->determineConfidenceInterval(last_pval_inside_curve,
// 						      last_pval_inside_other,
// 						      first_pval_outside_curve,
// 						      first_pval_outside_other,
// 						      prev_distance,
// 						      current_distance);
//     return succeeded;
// }

// //===========================================================================
// void BindToCurve::expandCoincidenceRegion(double& march_par,
// 					  double* other_par,
// 					  double& dist,
// 					  bool& passed_tolerance)
// //===========================================================================
// {
//     // calculating points
//     vector<Point> march_pt(3), other_pt(3); // march point and other point
//     marching_cv_->point(march_pt, march_par, 2, step_forward_);
//     other_cv_->point(other_pt, *other_par, 0);
//     bool other_forward = 
// 	march_pt[1] * other_pt[1] >= 0 ? step_forward_ : !step_forward_;
//     if (other_forward == false) {
// 	// re-calculating other point
// 	other_cv_->point(other_pt, *other_par, 0, other_forward);
//     }
//     other_pt[1] = other_pt[2] = Point(0.0, 0.0, 0.0);

//     // determining next critical values (passing them in one step can be dangerous)
//     double march_critical = marching_cv_->nextSegmentVal(march_par, step_forward_);
//     double other_critical = other_cv_->nextSegmentVal(*other_par, other_forward);
    
//     // calculate new step length
//     double step;
//     double aim = step_forward_ ? (1 + tol_->getEpsBracket()) * tol_->getEpsge() :
// 	                         (1 - tol_->getEpsBracket()) * tol_->getEpsge() ;
//     stepLength(march_pt, other_pt, aim, step_forward_, step);
//     if (hasPassedVal(march_par + step, march_critical, step_forward_)) {
// 	step = march_critical - march_par;
//     }
//     // determine search interval for other point
//     double pmin = other_forward ? *other_par : other_critical;
//     double pmax = other_forward ? other_critical : *other_par;
    
//     // march along curves
//     double fuzzy = fuzzy_eps_ * (other_cv_->endparam() - other_cv_->startparam());
//     bool exited = false;
//     double tolerance = tol_->getEpsge();
//     passed_tolerance = false;
//     shared_ptr<const ParamCurve> pc_other = other_cv_->getParamCurve();

//     // checking midpt. and endpt. of interval
//     dist = 0;
//     for (int i = 0; i < 2 && dist < tolerance; ++i){
// 	march_par += 0.5 * step;
// 	marching_cv_->getParamCurve()->point(march_pt[0], march_par);
// 	double prev_dist = dist = LARGEST_POSSIBLE_VALUE;
// 	bool repeat;
// 	do {
// 	    repeat = false;
// 	    pc_other->closestPoint(march_pt[0], 
// 				   pmin, pmax, *other_par, other_pt[0], dist, other_par);
// 	    other_cv_->knotIntervalFuzzy(*other_par, fuzzy);

// 	    if (dist > tolerance) {
// 		if (other_forward && fabs(*other_par - pmax) < fuzzy) {
// 		    if (pmax != other_cv_->endparam() && dist < prev_dist) {
// 			pmin = pmax;
// 			pmax = other_cv_->nextSegmentVal(pmax, true);
// 			repeat = true; // change interval in second curve and retry
// 		    }
// 		} else if (!other_forward && fabs(*other_par - pmin) < fuzzy) {
// 		    if (pmin != other_cv_->startparam() && dist < prev_dist) {
// 			pmax = pmin;
// 			pmin = other_cv_->nextSegmentVal(pmin, false);
// 			repeat = true; // change interval in second curve and retry
// 		    }
// 		}
// 	    }
// 	    prev_dist = dist;
// 	} while (repeat);
//     }
//     passed_tolerance = (dist > tolerance);
//     fuzzy = fuzzy_eps_ * (marching_cv_->endparam() - marching_cv_->startparam());
//     marching_cv_->knotIntervalFuzzy(march_par, fuzzy);
// }

// //===========================================================================
// void BindToCurve::determineConfidenceInterval(double& inside_marching_param,
// 					      double* inside_other_param,
// 					      double& outside_marching_param,
// 					      double* outside_other_param,
// 					      double& inside_dist,
// 					      double& outside_dist)
// //===========================================================================
// {
//     const double& center_value = tol_->getEpsge();
//     const double& bracket = tol_->getEpsBracket();
//     const double lower_limit = (1 - bracket) * center_value;
//     const double upper_limit = (1 + bracket) * center_value;

//     ASSERT(inside_dist < center_value && outside_dist > center_value);

//     bool outside_val_fixed = (outside_dist < upper_limit);
//     bool inside_val_fixed = (inside_dist > lower_limit);

//     Point p1, p2;
//     double new_dist;

//     // we know that 'center_value' is bracketed by 'inside_dist' and 'outside_dist'.
//     // Now we need to make the bracket narrow enough to fit inside [lower_limit, upper_limit]
//     shared_ptr<const ParamCurve> march_curve = marching_cv_->getParamCurve();
//     shared_ptr<const ParamCurve> other_curve = other_cv_->getParamCurve();
//     while (!inside_val_fixed || !outside_val_fixed) {
// 	double mid_p_marching = (inside_marching_param + outside_marching_param) * 0.5;
// 	double mid_p_other = (inside_other_param[0] + outside_other_param[0]) * 0.5;
// 	double pmin = std::min(inside_other_param[0], outside_other_param[0]);
// 	double pmax = std::max(inside_other_param[0], outside_other_param[0]);
// 	march_curve->point(p1, mid_p_marching);
// 	other_curve->closestPoint(p1, pmin, pmax, mid_p_other, p2, new_dist, &mid_p_other);
// 	if (new_dist > center_value) {
// 	    outside_marching_param = mid_p_marching;
// 	    outside_other_param[0] = mid_p_other;
// 	    outside_dist = new_dist;
// 	    outside_val_fixed = (outside_dist <= upper_limit);
// 	} else { // new_dist <= center_value
// 	    inside_marching_param = mid_p_marching;
// 	    inside_other_param[0] = mid_p_other;
// 	    inside_dist = new_dist;
// 	    inside_val_fixed = (new_dist >= lower_limit);
// 	}
//     }
// }
