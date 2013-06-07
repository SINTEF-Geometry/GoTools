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

#ifndef _BRENT_MINIMIZE_H
#define _BRENT_MINIMIZE_H

#include "GoTools/utils/GeneralFunctionMinimizer.h"

namespace Go
{
    /** Functor evaluation
     */

template<class Functor>
class Fun2Fun
{
public:
  /// Constructor
  /// \param f the functor
  /// \param a start parameter
  /// \param b end parameter
 Fun2Fun(const Functor& f, double a, double b) : f_(f), a_(a), b_(b) {}

  /// Evaluate functor
  double operator()(const double* arg) const
  {
    return f_(*arg);
  }
// void grad(const double* arg, double* grad) const;
  /// Return start parameter of functor
  double minPar(int n) const { return a_; }
  /// Return end parameter
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

