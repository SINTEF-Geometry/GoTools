//===========================================================================
//
// File : EvalFunctorCurve.C
//
// Created: Tue Oct 12 11:33:05 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================


#include "GoTools/isogeometric_model/EvalFunctorCurve.h"

using std::shared_ptr;

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
