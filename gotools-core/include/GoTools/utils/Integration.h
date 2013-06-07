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

#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_


#include <math.h>
#include <functional>
#include <numeric>
#include <limits>
#include <stdexcept>
#include "GoTools/utils/errormacros.h"

namespace Go {


/// This routine computes the n'th stage of refinement of an extended
/// trapezoidal rule. Used as a help routine for 'simpsons_rule', found
/// further down in this file.  NB: It only carries out the computation 
/// for ONE stage of the procedure, so if you want to use it to integrate
/// a function, you must call it successively for n = 1, n = 2, ....
template <typename Functor>
void trapezoidal(Functor& f, double a, double b, double& s, int n)
{
    DEBUG_ERROR_IF(n < 1, "Bad function argument to trapezoidal().");
    if (n == 1) {
	s = 0.5 * (b - a) * (f(a) + f(b));
    } else {
	int num_int_pts = 1 << (n-2);
	//int num_int_pts = std::power(2, n-2);
	double spacing = (b - a) / num_int_pts;
	double x = a + 0.5 * spacing;
	double sum = 0.0;
	for (int i = 0; i < num_int_pts; ++i) {
	    sum += f(x);
	    x += spacing;
	}
	s += (b - a) * sum / num_int_pts;
	s *= 0.5;
    }
}


/// Routine to calculate the integral of the functor f from a to b
/// using Simpson's rule. 
template <typename Functor>
double simpsons_rule(Functor& f, double a, double b,
		     const double eps = 1.0e-6, const int max_iter = 20)
{
    const double one_third = double(1) / double(3);
    double result = 0;
    double tz = 0;
    double last_tz = std::numeric_limits<double>::min();
    double last_result =  std::numeric_limits<double>::min();

    for (int j = 1; j <= max_iter; ++j) {
	trapezoidal(f, a, b, tz, j);
	result = (4.0 * tz - last_tz) * one_third;
	if ((fabs(result - last_result) < eps * fabs(last_result)) ||
	    (fabs(result) < eps && fabs(last_result) < eps && j > 6)) {
	    return result;
	}
	last_result = result;
	last_tz = tz;
    }
    MESSAGE("Too many steps in simpsons_rule.");
    throw std::runtime_error("Too many steps in simpsons_rule.");
} // 

/// Routine to calculate the integral of the functor f from a to b
/// using Gaussian quadrature with W=1 and N=10. 
template <typename Functor>
double gaussian_quadrature(Functor& f, double a, double b)
{
    static double x[] = { 0.1488743389, 
			  0.4333953941,
			  0.6794095682, 
			  0.8650633666, 
			  0.9739065285 };
    static double weight[] = { 0.2955242247, 
			       0.2692667193,
			       0.2190863625, 
			       0.1494513491, 
			       0.0666713443 };

    double midpt = 0.5 * (b + a);
    double half_length = 0.5 * (b - a);
    double scaled_result = double(0);
    double step = double(0);
    for (int i = 0; i < 5; ++i) {
	step = half_length * x[i];
	scaled_result += weight[i] * (f(midpt + step) + f(midpt - step));
    }
    // rescale result to correspond to actual interval size
    return scaled_result * half_length; 
}



/// Help functor to simpsons_rule2D and gaussian_quadrature2D
///\cond
template <typename Functor>
class Integrate2ndFunctor {
public:
    Integrate2ndFunctor(Functor& f, double a, double b)
	: f_(f), a_(a), b_(b) { }

    double operator() (double t) const
    {
	std::binder1st<Functor> fu(f_, t);
	return gaussian_quadrature(fu, a_, b_);
    }

private:
    Functor f_;
    double a_, b_;

};
///\endcond

/// Routine to integrate the two-dimensional functor f over the
/// rectangle defined by ax, bx, ay and by. Uses Simpson's rule.
template <typename Functor2D>
double simpsons_rule2D(Functor2D& f,
		       double ax, double bx, double ay, double by,
		       const double eps = 1.0e-6, const int max_iter = 20)
{
    Integrate2ndFunctor<Functor2D> fu(f, ay, by);
    return simpsons_rule(fu, ax, bx, eps, max_iter);
}


/// Routine to integrate the two-dimensional functor f over the
/// rectangle defined by ax, bx, ay and by. Uses Gaussian quadrature.
template <typename Functor2D>
double gaussian_quadrature2D(Functor2D& f,
			     double ax, double bx, double ay, double by)
{
    Integrate2ndFunctor<Functor2D> fu(f, ay, by);
    return gaussian_quadrature(fu, ax, bx);
}


} // namespace GO


#endif // _INTEGRATION_H

