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
    int g1 = first_ind + SplineUtils::closest_in_array(pt.begin(), 
					  &(*cv.coefs_begin()), 
					  nmb_coefs, 
					  cv.dimension());
    double seed = (cv.order() > 1) ? basis.grevilleParameter(g1) :
        0.5*(basis.begin()[g1] + basis.begin()[g1+1]);
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
// 	int g1 = first_ind + SplineUtils::closest_in_array(pt.begin(), &coefs_[first_ind*dim_],
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








