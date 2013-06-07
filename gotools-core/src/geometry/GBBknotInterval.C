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

#include "GoTools/geometry/BsplineBasis.h"
#include <algorithm>
#include <math.h>

using namespace Go;


//-----------------------------------------------------------------------------
int BsplineBasis:: knotInterval( double t) const
//-----------------------------------------------------------------------------
{
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To localize the point t in the array et_.
*              The output ileft should satisfy the relations
*                          
*                    et_[ileft] <= t < et_[ileft+1].
* 
*              There are two exceptions to this. If (ax >= et_[in_])
*              then ileft should be in_-1 (this corresponds to extending
*              the polynomial piece between et_[in_-1] and et_[in_] to the
*              right of the natural parameter interval.
*              Similarly, if (ax < et_[ik-1]) then ileft should still be
*              ik-1.
*
*
*
* INPUT      : et_     - Doublevector of dimension [in_+ik] containing
*                       the knot vector.
*              order_     - The polynomial order of the B-splines associated
*                       with et_.
*              in_     - The dimension of the spline space associated with
*                       the knot vector et_.
*              t     - The point at which the B-spline values and derivatives
*                       are to be computed.
*
*
*
* OUTPUT : ileft - Pointer to the interval in the knot vector
*                       where ax is located, check the relations above.
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The aim is to do as little work as possible in the cases
*              where ileft has the right or almost the right value.
*              First of all we make sure that ileft has a legal value
*              (a value in the range ik-1 to in_-1). Then we check
*              if the current value is OK.
*              If it is not we check that ax is in the interior of et_
*              or if the right value is obtained by either increasing
*              or decreasing ileft by 1. If the right value still has
*              not been found we do a binary search.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Knut Moerken, University of Oslo, August 1988.
* 
* MODIFIED   : Atgeirr F Rasmussen, Sintef, October 1998, August 2000
*
*********************************************************************
*/


    // Check the validity of the current BsplineBasis object.
    // Throws a CorruptData exception if something is wrong.
    // Not called if GO_NO_CHECKS was defined in
    // errormacros.h.
    //CHECK(this);

    // Make sure that last_knot_interval_ is in the legal range.
    int& ileft = last_knot_interval_;
    if (ileft < 0 || ileft > order_+num_coefs_-2)
	ileft = order_-1;

    // Check if the current value of ileft is acceptable.
  
    if (knots_[ileft] <= t && t < knots_[ileft+1])
	return ileft;  
    // Check if t is outside (knots_[order_-1],knots_[num_coefs_]).
    else if (t >= knots_[num_coefs_-1])
	ileft = num_coefs_ - 1;
    else if (t <= knots_[order_-1])
	ileft = order_ - 1;
  
    // Check if it is sufficient to increase or decrease ileft by one.

    else if (knots_[ileft+1] <= t && t < knots_[ileft+2])
	ileft += 1;
    else if (knots_[ileft-1] <= t && t < knots_[ileft])
	ileft -= 1;
  
    // Last resort - a binary search.
    else {
      
	// kmin and kmax gives the upper and lower limits on the possible 
	// values of ileft.
      
	int kmin,kmax;
      
	kmin = order_ - 1; kmax = num_coefs_ - 1;
	ileft = (kmin+kmax)/2;
	
	while (t < knots_[ileft] || knots_[ileft+1] <= t) {
	    if (t < knots_[ileft])
		kmax = ileft;
	    else
		kmin = ileft;
	    
	    ileft = (kmin+kmax)/2;
	}
    }

    return ileft;
}

//-----------------------------------------------------------------------------
int BsplineBasis:: knotIntervalFuzzy( double& t, double tol) const
//-----------------------------------------------------------------------------
{
    // Check the validity of the current BsplineBasis object.
    // Throws a CorruptData exception if something is wrong.
    // Not called if GO_NO_CHECKS was defined in
    // errormacros.h.
    
    knotInterval(t);
    if (t - knots_[last_knot_interval_] < tol) {
	t = knots_[last_knot_interval_];
    } else if (knots_[last_knot_interval_ + 1] - t < tol) {
	t = knots_[++last_knot_interval_];
	while (last_knot_interval_ < num_coefs_ &&
	       knots_[last_knot_interval_] == (knots_[last_knot_interval_+1])) {
	    ++last_knot_interval_;
	}
	if (last_knot_interval_ == num_coefs_) {
	    --last_knot_interval_;
	}
    }
    return last_knot_interval_;
}


// //-----------------------------------------------------------------------------
// int BsplineBasis:: knotIntervalFuzzy( double& t, double tol) const
// //-----------------------------------------------------------------------------
// {
//     // Check the validity of the current BsplineBasis object.
//     // Throws a CorruptData exception if something is wrong.
//     // Not called if GO_NO_CHECKS was defined in
//     // errormacros.h.
//     CHECK(this);

//     // A binary search through the knot vector
//     std::vector<double>::const_iterator lb
// 	= std::lower_bound(knots_.begin(), knots_.end(), t);
//     if (lb < knots_.end() && *lb - t < tol) {
// 	t = *lb;
//     } else if (lb > knots_.begin() && t - (*(lb-1))  < tol) {
// 	t = *(lb-1);
//     }
//     return knotInterval(t);
// }
