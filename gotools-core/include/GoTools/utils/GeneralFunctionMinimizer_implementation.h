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

#ifndef _GENERAL_FUNCTIONMINIMIZER_IMPLEMENTATION_H
#define _GENERAL_FUNCTIONMINIMIZER_IMPLEMENTATION_H

#include <limits>
#include <assert.h>
#include "GoTools/utils/errormacros.h"
#include <algorithm>
#include <iostream> // @@debug

namespace Go {

//===========================================================================
template<class Functor>
const double FunctionMinimizer<Functor>::root_machine_precision_ = 
sqrt(numerical_tolerance(1));

template<class Functor>
const double FunctionMinimizer<Functor>::rel_tol_ = 1.0e-6;

template<class Functor>
const double FunctionMinimizer<Functor>::perturbation_ = 1.0e-10; 

template<class Functor>
const double FunctionMinimizer<Functor>::golden_ratio_ = 0.61803399;

template<class Functor>
const double FunctionMinimizer<Functor>::default_partition_ = 
1.0 / 20.0; // default partition ratio of an interval

template<class Functor>
const double FunctionMinimizer<Functor>::tiny_ = 1.0e-20;  // a really tiny number

template<class Functor>
const int FunctionMinimizer<Functor>::max_iter_ = 100;
//===========================================================================

//===========================================================================
template<class Functor>
void minimise_conjugated_gradient(FunctionMinimizer<Functor>& dfmin)
//===========================================================================
{
    const double TOL = std::numeric_limits<double>::epsilon(); //1.0e-8;
    const double EPS = 1.0e-10;
    // minimising the 'dfmin' function using conjugated gradients.
    const int N = dfmin.numPars();
    Point gradient(N), old_gradient(N), dir(N);
    dfmin.grad(old_gradient);
    dir = -old_gradient;
    double old_val = dfmin.fval();
    while(true) {

	// make sure direction is not uphill (is this already guaranteed??)
	// and truncating it if we are at the border of the domain
	if (dir * old_gradient > 0) {
	    dir *= -1;
	}
	for (int i = 0; i < N; ++i) {
	    if ((dfmin.atMin(i) && dir[i] < 0) || (dfmin.atMax(i) && dir[i] > 0)) {
		dir[i] = 0;
		// conjugate gradient value is broken since we modified dir.
		// Restarting cycle.
	    }
	}

	// minimize along this direction
	bool hit_domain_edge = false;
	double new_val = dfmin.minimize(dir, hit_domain_edge); 
	if (2.0 * fabs(new_val - old_val) <= TOL  * (fabs(new_val) + fabs(old_val)+ EPS)) {
	    // we have reached a minimum
	    break;
	} else {
	    old_val = new_val;
	}

	// choose new direction 
	dfmin.grad(gradient);
	if (!hit_domain_edge) {
	    // we are still in a conjugated gradient cycle.  Choose new direction 
	    // using conjugated gradients (Polak-Ribiere variant)

	    Point diff = gradient - old_gradient;
	    double factor = gradient * diff / old_gradient.length2();
	    Point dir_saved = dir;
	    dir *= factor;
	    dir -= gradient;

	    if (dir * old_gradient > 0) {
		dir *= -1;
	    }

	    bool on_boundary = false;
	    for (int i = 0; i != N; ++i) {
		if ((dfmin.atMin(i) && dir[i] < 0) || (dfmin.atMax(i) && dir[i] > 0)) {
		    on_boundary = true;
		    gradient[i] = old_gradient[i] = 0;	
		    dir_saved[i] = 0; // we will have to recalculated 'dir' based on this vector
		}
	    }

	    if (on_boundary) {
		// the direction we choose will take us out of the domain.  Reduce problem
		// to conjugated gradient in a lower dimension.
		dir = dir_saved;
		diff = gradient - old_gradient;
		factor = gradient * diff / old_gradient.length2();
		dir *= factor;
		dir -= gradient;	
	    } 
	    old_gradient = gradient;
	} else {
	    // We ran into an edge of the domain.  Re-initialising conjugate
	    // gradient method using steepest descent.
	    dir = - gradient;
	    for (int i = 0; i != N; ++i) {
		if ((dfmin.atMin(i) && dir[i] < 0) || (dfmin.atMax(i) && dir[i] > 0)) {
		    dir[i] = 0;
		    gradient[i] = 0;
		}
	    }
	    old_gradient = gradient;
 	}
    }
}


//===========================================================================
template<class Functor>
inline FunctionMinimizer<Functor>::
FunctionMinimizer(int num_param, const Functor& fun, const double* const seed, double tol)
    : fun_(fun), par_(seed, seed + num_param), 
      at_min_(num_param), at_max_(num_param), cached_value_updated_(false), 
      cached_grad_updated_(false),  cached_grad_(num_param), 
      minimization_tol_(tol > root_machine_precision_ ? tol : root_machine_precision_), 
      param_tol_(num_param)
//===========================================================================
{
    for (int i = 0; i < num_param; ++i) {
        // We should scale only if domain size does not reflect geometry size
        // (i.e. not curve length parametrized basis).
	param_tol_[i] = rel_tol_;// * (fun_.maxPar(i) - fun_.minPar(i));
    }
    
    checkBorder(); // determine at_min_ and at_max_
    resetCache();
}

//===========================================================================
template<class Functor>
inline void FunctionMinimizer<Functor>::resetCache() const
//===========================================================================
{
    cached_value_updated_ = false;
    cached_grad_updated_ = false;
}

//===========================================================================
template<class Functor>
inline void FunctionMinimizer<Functor>::checkBorder()
//===========================================================================
{
    int n = numPars();
    for (int i = 0; i < n; ++i) {
	at_min_[i] = par_[i] < (fun_.minPar(i) + param_tol_[i]);
	at_max_[i] = par_[i] > (fun_.maxPar(i) - param_tol_[i]);
    }
}

//===========================================================================
template<class Functor>
inline double FunctionMinimizer<Functor>::minimize(const Point& dir, 
                                                   bool& hit_domain_edge, bool rerun)
//===========================================================================
{
    // flipping dir if it conflicts with gradient of function (will this ever happen??)
    Point oriented_dir(numPars());
    if (!orientDirection(dir, oriented_dir)) {
	// gradient of function was exactly perpendicular to 'dir'.  We are already at
	// a minimum along this direction
	return fval();
    }
    // deciding maximum possible steplength in the given direction
    double max_step = determineMaxSteplength(oriented_dir);
    if (max_step < tiny_) {
	// already on boundary in this direction.
	return fval();
    }
    // bracketing minimum in the given direction
    double bracket[3]; // startpoint, midpoint and endpoint of bracketing interval
    double fval_brak[3]; // function value at the three bracketing points
    
    if (bracketInterval(oriented_dir, max_step, bracket, fval_brak)) {
	// bracketing of interval successful.  Applying Brent
	hit_domain_edge = false;
	return linminBrent(oriented_dir, bracket, fval_brak);
    } 
    // We were not able to bracket the interval.  This may either be because the search
    // interval grew all the way to 'max_step' without being able to bracket a minimum,
    // or that the search interval shrank and collapsed in a singularity, indicating 
    // that although we know this to be a direction of descent, numerically the function
    // is not able to produce a smaller value in the neighbourhood of the departing point.
    if (fval_brak[1] >= fval_brak[0]) {
	ASSERT(bracket[2] < max_step);
	// the interval has collapsed onto the starting point.  We consider this to be 
	// a minimum (its directional derivative, although not zero, is so small that 
	// it does not modify the function within machine precision).
	return fval();
    }
    ASSERT(bracket[2] == max_step);
    // If we got here, we have maximally expanded the interval without being able to
    // verify that we have bracketed a minimum. The conclusion is that the minimum
    // is located very close to max_step.  Our strategy now is to find outwhat happens
    // at the end of the interval.  If the gradient is positive, we know that we have 
    // just passed the minimum, and we can search backwards from the endpoint.  If not,
    // we conclude that this is a minimum.  We check if this point gives a smaller 
    // function value than where we stand now, and if that is the case, we move there
    // before returning.
    Point end_grad(numPars());
    Point temp= par_ + max_step * oriented_dir;
    grad(temp, end_grad); 
    double max_step_val = fval(temp);

    if (scalarProductSign(end_grad, oriented_dir) == 1) {
	// derivative is POSITIVE at 'max_step'.  There MUST be a minimum that we 
	// missed.  Let us search backwards for it.
	if (!rerun) {
	    moveUV(oriented_dir, max_step); // moving to end of interval
	    double tmp = minimize(dir, hit_domain_edge, true); // search backwds
	    return tmp;
	}
	MESSAGE("Warning: unable to pinpoint exact minimum in FunctionMinimizer::minimize()");
    }
    // the derivative is NEGATIVE at 'max_step', and we conclude that this is a minimum
    // as good as any.  If it gives a smaller function value than where we stand now,
    // we will move there.  Otherwise, we stay put and consider ourselves already on
    // a minimum.
    if (max_step_val < fval()) {
	hit_domain_edge = true;
	moveUV(oriented_dir, max_step);
	return max_step_val;
    }
    // unable to find any point better than where we are currently standing
    hit_domain_edge = false;
    return fval();
}

//===========================================================================
template<class Functor> inline 
int FunctionMinimizer<Functor>::scalarProductSign(const Point& p1, const Point& p2)
//===========================================================================
{
    double prod = p1 * p2;
    double l1 = *std::max_element(p1.begin(), p1.end());
    double l2 = *std::max_element(p2.begin(), p2.end());
    // find the numerical tolerance related to the point coordinate of the largest
    // magnitude (its error dwarfs the errors of the others).
    double tol = 10 * numerical_tolerance(l1 > l2 ? l1 : l2);

    if (prod < -tol) return -1;
    if (prod > tol ) return 1;

    // indiscernable from 0
    return 0;
}

//===========================================================================
template<class Functor> inline bool 
FunctionMinimizer<Functor>::orientDirection(const Point& dir, Point& result)
//===========================================================================
{
    Point start_grad(numPars());
    grad(start_grad);

    result = dir / dir.lengthInf();

    int sign = scalarProductSign(start_grad, result);

    if (sign == 1) {
	result *= -1;
    }
    return (sign != 0);
}

//===========================================================================
template<class Functor>
inline double FunctionMinimizer<Functor>::linminBrent(const Point& dir, 
						      const double* bracket,
						      const double* fval_brak)
//===========================================================================
{
    const int MAXITER = 100;
    double a = bracket[0];  double fa = fval_brak[0];
    double b = bracket[1];  double fb = fval_brak[1];
    double c = bracket[2];  double fc = fval_brak[2];

    if (!(a <= b && b <= c)) {
	THROW("Error: Assertion 'a <= b && b <= c' failed.");
    }
    if (!(fa >= fb && fb <= fc)) {
	THROW("Error: Assertion 'fa >= fb && fb <= fc' failed.");
    }

    double b2, fb2, b3, fb3;
    if (fa < fc) {
	b2 = a; fb2 = fa;
	b3 = c; fb3 = fc;
    } else {
	b2 = c; fb2 = fc;
	b3 = a; fb3 = fa;
    }
    
    double u = b;
    double last_step = std::numeric_limits<double>::max();
    double before_last_step = last_step;
    double fu;
    int iter;
    // a step of 'delta' moves u and v by respectively delta * dir[0] and delta * dir[1].
    // therefore, if we want a precision of (minimization_tol_ * central bracket value)
    // in each parameter, we must use the tolerance below:
    const double central_value = fabs(b) > perturbation_ ? fabs(b) : perturbation_;
    double tol_1 = minimization_tol_ * central_value;
    // Since a, b, c and u are "fictive" function arguments (the real argument in R^n being 
    // 'par_ + u * dir'), we must check if the minimum tolerance of this last expression
    // exceeds the tolerance we have chosen.  In that case, increase the tolerance to this
    // new value.
    double tol_2 = numerical_tolerance(par_, dir); 
    const double tol = tol_1 > tol_2 ? tol_1 : tol_2;

    // three below variables to speed up convergence when suspecting to be really close
    // to the function minimum
    bool last_step_was_tol = false;
    bool last_step_was_positive = false;
    double btemp = b2;
    for (iter = 0; iter < MAXITER; ++iter) {
	// -------------------------------------
	// STEP 1: coming up with a point to try
	// -------------------------------------
	
	// trying parabolic interpolation
	bool use_parabolic = parabolicFit(b, fb, b2, fb2, b3, fb3, 0.5 * before_last_step, u);
	before_last_step = last_step;

	if (use_parabolic) {
	    // so far, it seems like we can use a parabolic fit.  But is it inside the interval?
	    use_parabolic = (u > a + tol && u < c - tol); // useful only if inside the interval
	}
	// here, we should know whether the parabolic fit is useful or not
	if (!use_parabolic) {
	    // we must come up with another point to try.  Resorting to golden search
	    const double rdiff = c - b;
	    const double ldiff = b - a;
	    if (rdiff < ldiff) {
		u = a + golden_ratio_ * ldiff;
	    } else {
		u = c - golden_ratio_ * rdiff;
	    }
	}
	// assuring a minimum steplength
	if (fabs(u - b) < tol) {
	    if (last_step_was_tol && btemp == b) { // exact equality is justified here!
		u = last_step_was_positive ? b - tol : b + tol;
	    } else {
		u = (u > b) ? b + tol : b - tol;
	    }
	    btemp = b;
	    last_step_was_positive = (u > b);
	    last_step_was_tol = true;

	} else {
	    last_step_was_tol = false;
	}
	last_step = fabs(b - u);
	
	// --------------------------------------------------------
	// STEP 2: trying the point and deciding what to do with it
	// --------------------------------------------------------
	
	// evaluating function in u (and evaluating gradient if premature end is allowed)
	
	fu = fval(par_ + u * dir);
	adjustBrackets(a, fa, b, fb, c, fc, b2, fb2, b3, fb3, u, fu);

	// end criteria
	if (fabs(b - (0.5 * (c + a))) < 2 * tol - 0.5 * (c - a)) {
	    break;
	}
    }
    if (iter == MAXITER) {
	MESSAGE("Failed to converge properly in linminBrent.");
	//throw runtime_error("Failed to converge in linminBrent.");
    }

    // setting current point at the found minimum
    moveUV(dir, b);    

    return fb;
}

//===========================================================================
template<class Functor>
inline void FunctionMinimizer<Functor>::adjustBrackets(double& a, double& fa,
						       double& b, double& fb,
						       double& c, double& fc,
						       double& b2, double& fb2,
						       double& b3, double& fb3,
						       double u, double fu)
//===========================================================================
{
    if (fu < fb) {
	// narrowing down brackets
	if (u < b) { 
	    c = b; fc = fb;
	} else {
	    a = b; fa = fb;
	}
	// updating minima
	b3 = b2; fb3 = fb2;
	b2 = b; fb2 = fb;
	b = u; fb = fu;
    } else {
	// we did not find a better minima, but we can still narrow down brackets
	if (u < b) {
	    a = u; fa = fu;
	} else {
	    c = u; fc = fu;
	}

	// checking if previous minima should be updated
	if (fu < fb2) {
	    b3 = b2; fb3 = fb2;
	    b2 = u; fb2 = fu;
	} else if (fu < fb3) {
	    b3 = u; fb3 = fu;
	}
    }
}

//===========================================================================
template<class Functor> inline double 
FunctionMinimizer<Functor>::determineMaxSteplength(const Point& dir) const
//===========================================================================
{
    const int n = numPars();
    std::vector<double> max_uv(n);

    for (int i = 0; i < n; ++i) {
	if (dir[i] > param_tol_[i]) {
	    double delta = fun_.maxPar(i) - par_[i];
	    max_uv[i] = delta / dir[i]; // a positive value
	} else if (dir[i] < -param_tol_[i]) {
	    double delta = fun_.minPar(i) - par_[i];
	    max_uv[i] = delta / dir[i]; // a positive value
	} else {
	    max_uv[i] = std::numeric_limits<double>::max(); // largest machine-repr. number
	}
    }
    return *(min_element(max_uv.begin(), max_uv.end()));
}

//===========================================================================
template<class Functor> inline double 
FunctionMinimizer<Functor>::parabolicEstimate(double p2, double fp2, const Point& dir)
//===========================================================================
{
    Point gradient(numPars());
    grad(gradient);
    double directional_deriv = gradient * dir;
    return -directional_deriv * p2 * p2 / (2 * (fp2 - fval() - p2 * directional_deriv));
}

//===========================================================================
template<class Functor> inline bool
FunctionMinimizer<Functor>::bracketInterval(const Point& dir,  // points to 2-array
					    const double max_step,
					    //const double oper_tol, // operating tolerance
					    double* bracket,   // points to 3-array
					    double* fval_brak) // points to 3-array
//===========================================================================
{
    // We suppose that 'dir' is a direciton of descent - anything else is an error!
    // we want to find three collinear points (a, b, c) on the line with direction 'dir'
    // and passing through the point (curU(), curV()) such that f(a) > f(b) < f(c). 

    if (!(max_step > 0)) {
	THROW("Error: Assertion 'max_step > 0' failed.");
    }
    //ASSERT(oper_tol > 0);
    //ASSERT(max_step > oper_tol);

    double& a = bracket[0]; // the three points we want to determine
    double& b = bracket[1];
    double& c = bracket[2];
    double& fa = fval_brak[0];
    double& fb = fval_brak[1];
    double& fc = fval_brak[2]; // function value in these points

    a = 0;
    fa = fval(); // current point (curU(), curV())

    b = default_partition_ * max_step;
    fb = fval(par_ + b * dir);

    // determining tolerance based on numerical precision of machine
    Point g;  grad(g);
    double num_tolerance = numerical_tolerance(par_, dir, fa, dir * g);

    if (fb >= fa) {
	// we know there is a minimum between a and b, since 'dir' is assumed to be a 
	// direction of descent.  We now only have to bisect the interval sufficiently
	// many times towards 'a' in order to bracket the minimum
	//std::cout << "Fallen into slow bracketing." << std::endl;
	do {
	    c = b; fc = fb;
	    b = parabolicEstimate(c, fc, dir);
	    fb = fval(par_ + b * dir);
	    if (!(b > a && b < c)) {
		THROW("Error: Assertion 'b > a && b < c' failed.");
	    }
	} while (fb >= fa && fabs(b-a) > num_tolerance);

	if (fb >= fa) {
	    // we did not succeed in bracketing a minimum between a and c, due to
	    // numerical issues (although 'dir' is direction of descent, we were
	    // unable to find a smaller function value within a distance from the
	    // startpoint bigger than the numerical tolerance).
	    return false;
	}
	if (!(fb < fa && fb <= fc)) {
	    THROW("Error: Assertion 'fb < fa && fb <= fc' failed.");
	}
	return true;
    }

    // if we got here, (fb < fa) --  we have (probably) not yet arrived at the minimum value
    // the below code is based on ideas from the book "Numerical Recipes in C" for bracketing
    // a minimum.
    c = b + (1 + golden_ratio_) * (b - a);
    c = (c < max_step) ? c : max_step;
    fc = fval(par_ + c * dir);

    //std::cout << "Entering normal bracketing." << std::endl;
    // we here know that a < b < c <= max_step, and that fb < fa
    while (fb > fc) {
	// coming up with suggestion for new point, based on parabolic interp. of a, b and c.
	double fu, u;
	parabolicFit(a, fa, b, fb, c, fc, max_step - b, u);

	double interval_length = c - a;
	if (fabs(u-b) > interval_length) {
	    u = (u > b) ? b + interval_length : b - interval_length;
	}
	u = (u < max_step) ? u : max_step; // not really necessary, but safeguard against roundoff
	
	// deciding what to do with new, estimated point
	if (b < u && u < c) {
	    fu = fval(par_ + u * dir);
	    if (fu < fc) { // minimum between b and c
		a = b; fa = fb;
		b = u; fb = fu;
		return true;
	    } else if (fu > fb) { // minimum between a and u
		c = u; fc = fu;
		return true;
	    }
	    // no help in parabolic fit
	    u = c + (1 + golden_ratio_) * (c - b);
	    u = (u < max_step) ? u : max_step;
	} else if (c < u) {
	    fu = fval(par_ + u * dir);
	    if (fu < fc) {
		shift3(b, c, u, u + (1 + golden_ratio_) * (u - c));
		u = (u < max_step) ? u : max_step;
		shift3(fb, fc, fu, fval(par_ + u * dir));
	    }
	} else { // reject parabolic u
	    u = c + (1 + golden_ratio_) * (c - b);
	    u = (u < max_step) ? u : max_step;
	    fu = fval(par_ + u * dir);
	}
	shift3(a, b, c, u);
	shift3(fa, fb, fc, fu);
	if (c == max_step && fc <= fb) { // yes, equality check of doubles is justified here
	    return false;
	}
// 	if (fabs(c - max_step) < minimization_tol_ && fc <= fb) {
// 		return false;
// 	}
    }
    // if we got here, then fa > fb < fc
    if (!(fb < fa && fb <= fc)) {
	THROW("Error: Assertion 'fb < fa && fb <= fc' failed.");
    }
    return true;
}

//===========================================================================
template<class Functor>
inline bool FunctionMinimizer<Functor>::parabolicFit(double a, double fa, 
						     double b, double fb, 
						     double c, double fc,
						     double max_from_b, double&u)
//===========================================================================
{
    // returns 'true' if it could determine a minimum on the parabola containing 
    // the points (a, fa), (b, fb) and (c, fc).  If no such minimum could be found
    // (ie. the function is linear), then it returns false.
    if (!(max_from_b > 0)) {
	THROW("Error: Assertion 'max_from_b > 0' failed.");
    }

    double tmp1 = (b - a) * (fb - fc);
    double tmp2 = (b - c) * (fb - fa);
    double denom = 2 * (tmp1 - tmp2); // denominator
    u = (b - a) * tmp1 - (b - c) * tmp2; // nominator
    
    if (fabs(u) < fabs(denom) * max_from_b) {
	u /= denom;
	u = b - u;
	return true;
    } else {
	// we would have to step farther than allowed.  Set to max and return false.
	u = (u * denom > 0) ? b - max_from_b : b + max_from_b;
	return false;
    }
    // never get here, but keep compiler happy
    return true;
}

//===========================================================================    
template<class Functor>
inline double FunctionMinimizer<Functor>::numerical_tolerance(double c)
//===========================================================================    
{
    if (c == double(0)) {
	return tiny_;
    }
    static const double eps = std::numeric_limits<double>::epsilon();
    return fabs(c) * eps;
}

//===========================================================================    
template<class Functor>
inline double FunctionMinimizer<Functor>::numerical_tolerance(const Point& x, 
							      const Point& dir)
//===========================================================================    
{
    // we want to determine the smallest positive double value 'a' such that 
    // (x + a * dir) != x.  If we define tol(x_i) to be the numerical precision
    // around the i-component of x, and dir_i to be the i-component of y, 
    // then we search for the value:
    // minimize_over_i (tol(x_i) / |dir_i|)

    double res = std::numeric_limits<double>::max();
    for (int i = 0; i < x.size(); ++i) {
	if (dir[i] != double(0.0)) {
	    double temp = numerical_tolerance(x[i]) / fabs(dir[i]);
	    res = res < temp ? res : temp;
	    //double temp = fabs(x[i]) * std::numeric_limits<double>::epsilon(); //tol(x_i)
	    //temp /= fabs(dir[i]); // |dir_i|
	    //res = res < temp ? res : temp;
	}
    }
    return res;
}

//===========================================================================    
template<class Functor>
inline double FunctionMinimizer<Functor>::numerical_tolerance(const Point& x, 
							      const Point& dir,
							      const double fx,
							      const double dfx)
//===========================================================================    
{
    //    we want to determine the smallest positive double value 'a' such that 
    // -> (x + a * dir) != x
    //    and such that 
    // -> f(x + a * dir) != f(x).  
    //    To satisfy the last criteria, we use the function value in x, given by
    //    'fx', as well as the directional derivative, which is given by 'dfx'. 


    double a1 = numerical_tolerance(x, dir);

    double a2 = (fabs(fx) + perturbation_) * std::numeric_limits<double>::epsilon();
    a2 /= fabs(dfx);

    return a1 > a2 ? a1 : a2;
}

//===========================================================================    
template<class Functor>
inline double FunctionMinimizer<Functor>::fval() const
//===========================================================================
{
    if (!cached_value_updated_) {
	cached_value_ = fun_(par_.begin());
	cached_value_updated_ = true;
    }
    return cached_value_;
}

//===========================================================================
template<class Functor>
inline double FunctionMinimizer<Functor>::fval(const Point& par) const
//===========================================================================
{
    return fun_(par.begin());
}

//===========================================================================
template<class Functor>
inline void FunctionMinimizer<Functor>::grad(Point& result) const
//===========================================================================
{
    if (!cached_grad_updated_) {
	fun_.grad(par_.begin(), cached_grad_.begin());
	cached_grad_updated_ = true;
    }
    result = cached_grad_;
}

//===========================================================================
template<class Functor>
inline void FunctionMinimizer<Functor>::grad(const Point& param, Point& result) const
//===========================================================================
{
    fun_.grad(param.begin(), result.begin());
}

//===========================================================================
template<class Functor>
inline void FunctionMinimizer<Functor>::moveUV(const Point& dir, double multiplier)
//===========================================================================
{
    par_ += dir * multiplier;
    checkBorder();
    resetCache();
}


}; // end namespace Go

#endif // _GENERAL_FUNCTIONMINIMIZER_IMPLEMENTATION_H

