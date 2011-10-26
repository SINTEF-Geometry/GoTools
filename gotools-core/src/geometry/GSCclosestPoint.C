//===========================================================================
//                                                                           
// File: GSCclosestPoint.C                                                   
//                                                                           
// Created: Thu Oct 19 15:49:45 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: GSCclosestPoint.C,v 1.39 2008-11-28 08:12:41 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineUtils.h"
#include <vector>

using namespace std;
using namespace Go;

namespace {

//===========================================================================
double choose_seed(const Point& pt, const SplineCurve& cv,
		   double tmin, double tmax) 
//===========================================================================
{
    const BsplineBasis& basis = cv.basis();
    int first_ind, last_ind, dummy_ind;
    basis.coefsAffectingParam(tmin, first_ind, dummy_ind);
    basis.coefsAffectingParam(tmax, dummy_ind, last_ind);
    int nmb_coefs = last_ind - first_ind + 1;
    int g1 = first_ind + closest_in_array(pt.begin(), 
					  &(*cv.coefs_begin()), 
					  nmb_coefs, 
					  cv.dimension());
    double seed =  basis.grevilleParameter(g1);
    seed = std::max(seed, tmin);
    seed = std::min(seed, tmax);
    return seed;
}

}; // end anonymous namespace 

namespace Go
{

//===========================================================================
void SplineCurve::closestPoint(const Point& pt,
			       double tmin,
			       double tmax,
			       double& clo_t,
			       Point& clo_pt,
			       double& clo_dist,
			       double const *seed) const
//===========================================================================
{
    double guess_param = seed ? *seed : choose_seed(pt, *this, tmin, tmax);
    ParamCurve::closestPointGeneric(pt, tmin, tmax, guess_param, clo_t, clo_pt, clo_dist);
}

};



// This is the old implementation.  It is now deprecated, but kept here for 
// future reference.

// namespace Go
// {


//     /** Missing doxygen documentation
//      *
//      */

// class closest_point_functor
// {
// public:
//     closest_point_functor(const SplineCurve* cv, const Point& pt)
// 	: cv_(cv), pt_(pt), x_(2,Point(cv_->dimension())) {}
//     double operator() (double t)
// 	{
// 	    cv_->point(x_, t, 1);
// 	    x_[0] -= pt_;
// 	    return x_[0]*x_[1];
// 	}
// private:
//     const SplineCurve* cv_;
//     const Point& pt_;
//     std::vector<Point> x_;
// };



// //-----------------------------------------------------------------------------
// void SplineCurve::closestPoint(const Point&   pt,
// 				 double    tmin,
// 				 double    tmax,
// 				 double&   clo_t,
// 				 Point&  clo_pt,
// 				 double&   clo_dist,
// 				 double const *seed) const
// //-----------------------------------------------------------------------------
// {
//     double guess_param;
//     if (seed != 0) {
// 	guess_param = *seed;
// 	if (guess_param < tmin || guess_param > tmax) {
// 	    MESSAGE("Suggested parameter for closest point "
// 		       "must lie inside domain!");
// 	    guess_param = (tmax < guess_param ? tmax : guess_param);
// 	    guess_param = (tmin > guess_param ? tmin : guess_param);
// //  	    guess_param = max(tmin, min(tmax, guess_param));
// 	}
//     } else {
// 	// We must make sure that we start searching inside [tmin, tmax].
// 	int first_ind, last_ind, dummy_ind;
// 	// We choose the average index.
// 	basis_.coefsAffectingParam(tmin, first_ind, dummy_ind);
// 	//int avg_first_ind = (int) 0.5*(first_ind + dummy_ind);
// 	basis_.coefsAffectingParam(tmax, dummy_ind, last_ind);
// 	//int avg_last_ind = (int) 0.5*(dummy_ind + last_ind);
// 	int nmb_coefs = last_ind - first_ind + 1;
// 	int g1 = first_ind + closest_in_array(pt.begin(), &coefs_[first_ind*dim_],
// 					      nmb_coefs, dim_);
// 	guess_param = basis_.grevilleParameter(g1);
//     }

//     if (guess_param < tmin || guess_param > tmax) {
// // 	MESSAGE("Failed finding guess_param inside legal interval, should not happen, fix!!!");
// 	// May happen for instance with cv without inner knots.
// 	guess_param = (guess_param < tmin) ? tmin : tmax;
// // 	guess_param = 0.5*(tmin + tmax); // At least it is better than outside domain.
//     }

//     // The new version of closestPointGeneric does not like guess_param at a bd parameter.
//     if (guess_param == tmin)
// 	guess_param = tmin + 0.01*(tmax - tmin);
//     else if (guess_param == tmax)
// 	guess_param = tmax - 0.01*(tmax - tmin);

//     ParamCurve::closestPointGeneric(pt, tmin, tmax, guess_param, 
// 				    clo_t, clo_pt, clo_dist);

//     if (clo_t < tmin || clo_t > tmax) {
//       MESSAGE("Returned closest point outside domain, should be moved inside?");
//       if (clo_t < tmin) {
// 	clo_t = tmin;
//       } else if (clo_t > tmax) {
// 	clo_t = tmax;
//       }
//     }

// }


// } // namespace Go








