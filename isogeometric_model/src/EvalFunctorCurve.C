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

#include "GoTools/isogeometric_model/EvalFunctorCurve.h"

namespace Go
{

  //===========================================================================
  EvalFunctorCurve::EvalFunctorCurve(BdCondFunctor* fbd, shared_ptr<SplineCurve> geo_curve, int dimension):
    functor_(fbd),
    geo_curve_(geo_curve),
    dim_(dimension)
  //===========================================================================
  {
  }


  //===========================================================================
  EvalFunctorCurve::~EvalFunctorCurve()
  //===========================================================================
  {
  }


  //===========================================================================
  Point EvalFunctorCurve::eval( double t) const
  //===========================================================================
  {
    Point geo_pt;
    geo_curve_->point(geo_pt, t);
    return functor_->evaluate(geo_pt);
  }


  //===========================================================================
  void EvalFunctorCurve::eval( double t, int n, Point der[]) const
  //===========================================================================
  {
    if (n >= 0)
      der[0] = eval(t);
    if (n >= 1)
      {
	double delta_default = 0.01;
	double delta;
	double start_par = start();
	double end_par = end();
	if (t + delta_default <= end_par)
	  delta = delta_default;
	else if (t - delta_default >= start_par)
	  delta = -delta_default;
	else if (end_par - t >= t - start_par)
	  delta = end_par - t;
	else
	  delta = start_par - t;
	der[1] = (1.0 / delta) * (eval(t + delta) - der[0]);
      }
    if (n >= 2)
      {
	Point zero(dim());
	for (int i = 0; i < dim(); ++i)
	  zero[i] = 0.0;
	for (int i = 2; i <= n; ++i)
	  der[i] = zero;
      }
  }


  //===========================================================================
  double EvalFunctorCurve::start() const
  //===========================================================================
  {
    return geo_curve_->startparam();
  }


  //===========================================================================
  double EvalFunctorCurve::end() const
  //===========================================================================
  {
    return geo_curve_->endparam();
  }


  //===========================================================================
  int EvalFunctorCurve::dim() const
  //===========================================================================
  {
    return dim_;
  }


  //===========================================================================
  bool EvalFunctorCurve::approximationOK(double par, Point approxpos,
					 double tol1, double tol2) const
  //===========================================================================
  {
    return true;
  }


} // namespace Go
