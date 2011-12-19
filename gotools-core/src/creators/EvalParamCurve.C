//===========================================================================
//                                                                           
// File: EvalParamCurve.C                                                    
//                                                                           
// Created: September 2010
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================



#include "GoTools/creators/EvalParamCurve.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/ParamSurface.h"

#include <fstream>

using namespace std;
using namespace Go;

//===========================================================================
EvalParamCurve::EvalParamCurve(shared_ptr<Go::ParamCurve>& crv)
    : crv_(crv)
//===========================================================================
{
}


//===========================================================================
EvalParamCurve::~EvalParamCurve()
//===========================================================================
{
}


//===========================================================================
Point EvalParamCurve::eval(double t) const
//===========================================================================
{
    return crv_->ParamCurve::point(t);
}


//===========================================================================
void EvalParamCurve::eval(double t, int n, Point der[]) const
//===========================================================================
{
  vector<Point> pts(n);
  crv_->point(pts, t, n);
  for (int ki=0; ki<n; ++ki)
    der[ki] = pts[ki];
}


//===========================================================================
double EvalParamCurve::start() const
//===========================================================================
{
  return crv_->startparam();
}


//===========================================================================
double EvalParamCurve::end() const
//===========================================================================
{
  return crv_->endparam();
}


//===========================================================================
int EvalParamCurve::dim() const
//===========================================================================
{
  return crv_->dimension(); 
}


//===========================================================================
bool EvalParamCurve::approximationOK(double par, Point approxpos,
				  double tol1, double tol2) const
//===========================================================================
{
    // Only first tolerance is used.
  Point pos = eval(par);

  double dist = pos.dist(approxpos); // Both points are in space.

  return (dist < tol1);
}

void EvalParamCurve::write(std::ostream& out) const
{
  crv_->writeStandardHeader(out);
  crv_->write(out);
  shared_ptr<CurveOnSurface> sf_cv = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv_);
  if (sf_cv.get())
    {
      sf_cv->underlyingSurface()->writeStandardHeader(out);
      sf_cv->underlyingSurface()->write(out);
    }
  }
