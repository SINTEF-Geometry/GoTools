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

#include "GoTools/isogeometric_model/EvalFunctorSurface.h"

namespace Go
{

  //===========================================================================
  EvalFunctorSurface::EvalFunctorSurface(BdCondFunctor* fbd, shared_ptr<SplineSurface> geo_surface, int dimension):
    functor_(fbd),
    geo_surface_(geo_surface),
    dim_(dimension)
  //===========================================================================
  {
  }


  //===========================================================================
  EvalFunctorSurface::~EvalFunctorSurface()
  //===========================================================================
  {
  }


  //===========================================================================
    Point EvalFunctorSurface::eval( double u, double v) const
  //===========================================================================
  {
    Point geo_pt;
    geo_surface_->point(geo_pt, u, v);
    return functor_->evaluate(geo_pt);
  }


  //===========================================================================
    void EvalFunctorSurface::eval( double u, double v, int n, Point der[]) const
  //===========================================================================
  {
    if (n >= 0)
	der[0] = eval(u, v);

    if (n >= 2)
	MESSAGE("eval(): Higher than first order derivs not supported yet.");

    if (n >= 1)
      {
	double delta_default = 0.01;
	double delta_u, delta_v;
	double start_par_u = start_u();
	double end_par_u = end_u();
	if (u + delta_default <= end_par_u)
	  delta_u = delta_default;
	else if (u - delta_default >= start_par_u)
	  delta_u = -delta_default;
	else if (end_par_u - u >= u - start_par_u)
	  delta_u = end_par_u - u;
	else
	  delta_u = start_par_u - u;
	der[1] = (1.0 / delta_u) * (eval(u + delta_u, v) - der[0]);

	double start_par_v = start_v();
	double end_par_v = end_v();
	if (v + delta_default <= end_par_v)
	  delta_v = delta_default;
	else if (v - delta_default >= start_par_v)
	  delta_v = -delta_default;
	else if (end_par_v - v >= v - start_par_v)
	  delta_v = end_par_v - v;
	else
	  delta_v = start_par_v - v;
	der[2] = (1.0 / delta_v) * (eval(u, v + delta_v) - der[0]);
      }
    if (n >= 2)
      {
	Point zero(dim());
	int num_pts = (n+1)*(n+2)/2;
	for (int i = 0; i < dim(); ++i)
	  zero[i] = 0.0;
	for (int i = 3; i < num_pts; ++i)
	  der[i] = zero;
      }
  }


  //===========================================================================
  double EvalFunctorSurface::start_u() const
  //===========================================================================
  {
      return geo_surface_->basis_u().startparam();
  }



  //===========================================================================
  double EvalFunctorSurface::start_v() const
  //===========================================================================
  {
      return geo_surface_->basis_v().startparam();
  }


  //===========================================================================
  double EvalFunctorSurface::end_u() const
  //===========================================================================
  {
      return geo_surface_->basis_u().endparam();
  }


  //===========================================================================
  double EvalFunctorSurface::end_v() const
  //===========================================================================
  {
      return geo_surface_->basis_v().endparam();
  }


  //===========================================================================
  int EvalFunctorSurface::dim() const
  //===========================================================================
  {
    return dim_;
  }


  //===========================================================================
    bool EvalFunctorSurface::approximationOK(double par_u, double par_v, Point approxpos,
					     double tol1, double tol2) const
  //===========================================================================
  {
    return true;
  }


} // namespace Go
