//==========================================================================
//                                                                          
// File: CurvatureFunctions.C                                                
//                                                                          
// Created: Fri Mar  3 14:18:46 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: CurvatureFunctions.C,v 1.2 2007-01-15 10:12:39 vsk Exp $
//                                                                          
// Description: This file contains curvature related functions that
// has been removed from ParamSurfaceInt.C and SfSelfIntersector.C
//                                                                          
//==========================================================================


//#if 0

#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/SfSelfIntersector.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/utils/Point.h"
//#include "GoTools/utils/GeneralFunctionMinimizer.h"
//#include "FunctionMinimizer.h"
//#include "GoTools/utils/brent_minimize.h"
#include "GoTools/geometry/GeometryTools.h"
#include <vector>

// From ParamSurfaceInt.C

using namespace Go;
using std::vector;
using std::pair;
using std::make_pair;


class curve_curvature_radius_functor {
public:
    curve_curvature_radius_functor(const ParamCurve* cv)
	: cv_(cv), x_(3,Point(cv_->dimension())) {}

    ~curve_curvature_radius_functor() {}

    double operator() (const double *t) const
    {
	// To compute the curvature, we must assume that curve is
	// curve length parametrized. Thus we're reparametrizing with
	// g such that g' = f'/|f'|.
	// g'' = (f''*|f'|^2 - (f'*f'')*f')/(|f'|^3).
	x_.resize(3);
	cv_->point(x_, *t, 2);

	Point double_der 
	    = (x_[1].length()*x_[2] 
	       - ((x_[1]*x_[2])*x_[1]) / x_[1].length())
	    / (x_[1].length()*x_[1].length());

	return 1.0/(double_der.length());
    }

    double grad(const double* arg, double* grad) const
    {
	ptvec_.resize(2);
	cv_->point(ptvec_, *arg, 1);
	*grad = 2 * (ptvec_[0] - x_[0]) * ptvec_[1];
	return ptvec_[0].dist2(x_[0]);
    }


    double minPar(int n) const {return cv_->startparam();}
    double maxPar(int n) const {return cv_->endparam();}

private:
    const ParamCurve* cv_;
    mutable std::vector<Point> x_;
    mutable std::vector<Point> ptvec_;
};


// From ParamSurfaceInt.C

// Warning from gcc: "‘double curvature(Go::ParamCurve*, double)’
//defined but not used". We out-comment this function.
// //===========================================================================
// static double curvature(ParamCurve* crv, double tpar)
// //===========================================================================
// {
//     ALWAYS_ERROR_IF(tpar < crv->startparam() || tpar > crv->endparam(),
// 		    "Evaluating outside curve.");

//     curve_curvature_radius_functor f(crv);
//     return 1/f(&tpar);
// }


// From ParamSurfaceInt.C

//===========================================================================
static double maxCurvature(ParamCurve* crv, double tmin, double tmax,
			   double guess_param, double& max_curv_param)
//===========================================================================
{
    // @@sbr Method altered, not yet fully tested.
    MESSAGE("Experimental, verify that the method works!");

    // Note: gcc complains about '"/*" within comment'. These comments
    // need cleaning.

// /*    curve_curvature_radius_functor f(crv);
//     if (false)
//     {
// 	/*const double TOL = 1.0e-08;
//     FunctionMinimizer<curve_curvature_radius_functor> funmin(1, f,
// 							     &guess_param,
// 							     TOL);
//     Point dir(1, 1.0);
// //     distfun.grad(&guess_param, &(dir[0])); // determining direction

//     bool hit_domain_edge = false;
//     funmin.minimize(dir, hit_domain_edge);

//     max_curv_param = funmin.getPar(0);
//     }

// //     FunctionMinimizer<curve_curvature_radius_functor> funmin(1, f,
// // 							     &guess_param,
// // 							     TOL);
// //     distfun.grad(&guess_param, &(dir[0])); // determining direction

// //     max_curv_param = guess_param;
// //     double max_curv;
// //     double c = 0.5*(tmin+tmax);
// //     max_curv = brent_minimize(f, tmin, c, tmax, max_curv_param, TOL);


//     else
//     {
// // Not used anymore, but kept for reference.
//     double ax = guess_param;
//     double step = (tmax - tmin)*1e-05;
//     double bx = (guess_param - step < tmin) ?
// 	guess_param + step : guess_param - step;
//     double cx, fa, fb, fc;
//     bracket_minimum(f, ax, bx, cx, fa, fb, fc); // We try to bracket solution.

//     ax = std::max(tmin, std::min(tmax, ax)); // Function may have
// 					     // iterated outside
// 					     // boundaries.
//     bx = std::max(tmin, std::min(tmax, bx));
//     cx = std::max(tmin, std::min(tmax, cx));
//     double tolerance = 1e-14;
//     double bottom_t;
//     if ((fabs(tmin-bx) < tolerance) || (fabs(bx-tmax) < tolerance)) {
// 	// Minimum is an end point, we return that value as the closest point.
// 	bottom_t = bx;
//     } else { // Solution is bracketed, we locate the minimum.
// 	try {
// 	    brent_minimize(&f, ax, bx, cx, bottom_t, tolerance);
// 	} catch (const IterationCountExceeded& ob) {
// 	    MESSAGE("Error in solving equation. Iterationcount exceeded. "
// 		    "Returning best guess.");
// 	    bottom_t = ob.guess_par;
// 	} catch (...) {
// 	    THROW("Unexpected exception occured!");
// 	}
//     }

//     // debugging
//     bool debug = false;
//     if (debug) {
// 	int nmb_samples = 10001;
// 	double from = crv->startparam();
// 	double to = crv->endparam();
// 	double step = (to - from)/(nmb_samples-1);
// 	std::ofstream outfile("data/curv_radius_plot.dta");
// 	for (int i = 0; i < nmb_samples; ++i) {
// 	    double tpar = from + i*step;
// 	    double tpar_curv = f(tpar);
// 	    outfile << tpar << " " << tpar_curv << std::endl;
// 	}
//     } // end debugging
//     }


// //     max_curv_param = bottom_t;
//     return curvature(crv, max_curv_param);*/
    return 0;
    
    // /*curve_curvature_radius_functor f(crv);

    // double ax = guess_param;
    // double step = (tmax - tmin)*1e-05;
    // double bx = (guess_param - step < tmin)
    // 	? guess_param + step : guess_param - step;
    // double cx, fa, fb, fc;
    // bracket_minimum(f, ax, bx, cx, fa, fb, fc); // We try to bracket
    // 						// solution.

    // ax = std::max(tmin, std::min(tmax, ax)); // Function may have
    // 					     // iterated outside
    // 					     // boundaries.
    // bx = std::max(tmin, std::min(tmax, bx));
    // cx = std::max(tmin, std::min(tmax, cx));
    // double tolerance = 1e-14;
    // double bottom_t;
    // if ((fabs(tmin-bx) < tolerance) || (fabs(bx-tmax) < tolerance)) {
    // 	// Minimum is an end point, we return that value as the
    // 	// closest point.
    // 	bottom_t = bx;
    // } else { // Solution is bracketed, we locate the minimum.
    // 	try {
    // 	    brent_minimize(&f, ax, bx, cx, bottom_t, tolerance);
    // 	} catch (const IterationCountExceeded& ob) {
    // 	    MESSAGE("Error in solving equation. Iterationcount exceeded. "
    // 		    "Returning best guess.");
    // 	    bottom_t = ob.guess_par;
    // 	} catch (...) {
    // 	    THROW("Unexpected exception occured!");
    // 	}
    // }

    // // debugging
    // bool debug = false;
    // if (debug) {
    // 	int nmb_samples = 10001;
    // 	double from = crv->startparam();
    // 	double to = crv->endparam();
    // 	double step = (to - from)/(nmb_samples-1);
    // 	std::ofstream outfile("data/curv_radius_plot.dta");
    // 	for (int i = 0; i < nmb_samples; ++i) {
    // 	    double tpar = from + i*step;
    // 	    double tpar_curv = f(tpar);
    // 	    outfile << tpar << " " << tpar_curv << std::endl;
    // 	}
    // } // end debugging
    
    // max_curv_param = bottom_t;
    // return curvature(crv, max_curv_param);*/
}


// Public member of ParamSurfaceInt

//===========================================================================
void ParamSurfaceInt::maxCurvatures(bool dir_is_u, int nmb_params,
				    double iso_from, double iso_to,
				    double guess_param,
				    vector<double>& max_curv_params,
				    vector<double>& max_curvatures)
//===========================================================================
{
    // Iterate to maximum curvature in a number of iso-curves

    vector<double> params(nmb_params);
    double ta = (dir_is_u) ? startParam(1) : startParam(0);
    double tb = (dir_is_u) ? endParam(1) : endParam(0);
    double tdel = (tb - ta)/(double)(nmb_params-1);
    double par = ta;

    max_curv_params.resize(nmb_params);
    max_curvatures.resize(nmb_params);
    for (int ki=0; ki<nmb_params; ki++, par+=tdel)
    {
      vector<shared_ptr<ParamCurve> > constcrvs
	  = surf_->constParamCurves(par, dir_is_u);
      double max_curv=0, max_curv_par;
      for (size_t kj=0; kj<constcrvs.size(); kj++)
      {
          max_curv = maxCurvature(constcrvs[kj].get(), iso_from, iso_to,
                                  guess_param, max_curv_par);
          if (kj==0 || max_curv > max_curvatures[ki])
          {
              max_curvatures[ki] = max_curv;
              max_curv_params[ki] = max_curv_par;
          }
      }
    }
}


// Private member of SfSelfIntersector

//===========================================================================
void
SfSelfIntersector::getMaxCurvatures(shared_ptr<ParamSurfaceInt> surf, int nsample,
				    vector<pair<double,double> >& max_curv1,
				    vector<pair<double,double> >& max_curv2)
//===========================================================================
{
    double ta1 = surf->startParam(0);
    double ta2 = surf->startParam(1);
    double tb1 = surf->endParam(0);
    double tb2 = surf->endParam(1);

    //int nsample = 3;

    vector<double> curv, par;
    surf->maxCurvatures(true, nsample, ta1, tb1, 0.5*(ta1+tb1),
			 par, curv);
    for (size_t ki=0; ki<curv.size(); ki++)
	max_curv1.push_back(make_pair(curv[ki], par[ki]));

    curv.clear();
    par.clear();
    surf->maxCurvatures(false, nsample, ta2, tb2, 0.5*(ta2+tb2),
			 par, curv);
    for (size_t ki=0; ki<curv.size(); ki++)
	max_curv2.push_back(make_pair(curv[ki], par[ki]));
}
 

//===========================================================================


//#endif // if 0
