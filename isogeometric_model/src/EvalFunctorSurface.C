//===========================================================================
//                                                                           
// File: EvalFunctorSurface.C                                                
//                                                                           
// Created: Wed Sep 19 15:45:14 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================



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
      MESSAGE("eval(): Not implemented yet.");

#if 0
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
#endif
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
