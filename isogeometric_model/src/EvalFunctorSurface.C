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
