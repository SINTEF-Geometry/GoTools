//===========================================================================
//                                                                           
// File: brent_minimize.h                                                    
//                                                                           
// Created: Fri Jan 20 09:06:21 2006                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: brent_minimize.h,v 1.2 2006-01-23 16:22:46 afr Exp $
//                                                                           
//===========================================================================

#ifndef _BRENT_MINIMIZE_H
#define _BRENT_MINIMIZE_H

#include "GoTools/utils/GeneralFunctionMinimizer.h"

namespace Go
{
    /** Brief description. 
     *  Detailed description.
     */

template<class Functor>
class Fun2Fun
{
public:
    Fun2Fun(const Functor& f, double a, double b) : f_(f), a_(a), b_(b) {}
    double operator()(const double* arg) const
    {
	return f_(*arg);
    }
/// void grad(const double* arg, double* grad) const;
    double minPar(int n) const { return a_; }
    double maxPar(int n) const { return b_; }
private:
    Functor f_;
    double a_;
    double b_;
};

//===========================================================================
template<class Functor>
inline double brent_minimize(const Functor& f,
			     double a, double b, double c,
			     double& parmin,
			     const double rel_tolerance = std::sqrt(std::numeric_limits<double>::epsilon()))
//===========================================================================
{
    Fun2Fun<Functor> f2(f, a, c);
    Go::FunctionMinimizer<Fun2Fun<Functor> > fmin(1, f2, &a, rel_tolerance);
    Go::Point dir(1);
    dir[0] = 1.0;
    double bracket[3];
    double fval_brak[3];
    bracket[0] = 0.0;
    bracket[1] = b-a;
    bracket[2] = c-a;
    fval_brak[0] = f(a);
    fval_brak[1] = f(b);
    fval_brak[2] = f(c);
    double minimum = fmin.linminBrent(dir, bracket, fval_brak);
    parmin = fmin.getPar(0);
    return minimum;
}


} // namespace Go

#endif // _BRENT_MINIMIZE_H

