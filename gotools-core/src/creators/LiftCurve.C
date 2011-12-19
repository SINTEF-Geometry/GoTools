//===========================================================================
//                                                                           
// File: LiftCurve.C                                                    
//                                                                           
// Created: Sat Mar 29 09:00:22 2003                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: LiftCurve.C,v 1.3 2005-11-17 08:55:30 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================



#include "GoTools/creators/LiftCurve.h"

#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/utils/config.h"
#include "GoTools/geometry/ParamCurve.h"

#include <fstream>

using namespace std;
using namespace Go;


//===========================================================================
LiftCurve::LiftCurve(shared_ptr<Go::ParamCurve>& parameter_crv,
		     shared_ptr<Go::ParamSurface>& surf,
		     double epsgeo)
    : parameter_crv_(parameter_crv), surf_(surf), epsgeo_(epsgeo)
//===========================================================================
{
    // Test input
    ALWAYS_ERROR_IF(parameter_crv_.get() == 0 || surf_.get() == 0,
		"Missing input data.");
    ALWAYS_ERROR_IF(parameter_crv_->dimension() != 2 || surf_->dimension() != 3,
		"Dimension mismatch.");
}


//===========================================================================
LiftCurve::~LiftCurve()
//===========================================================================
{
}


//===========================================================================
Point LiftCurve::eval(double t) const
//===========================================================================
{
    Point par_pt = parameter_crv_->ParamCurve::point(t);
    Point space_pt = surf_->ParamSurface::point(par_pt[0], par_pt[1]);

    return space_pt;
}


//===========================================================================
void LiftCurve::eval(double t, int n, Point der[]) const
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  if (n == 0)
      der[0] = eval(t);
  else {
      vector<Point> par_pt(2);
      parameter_crv_->point(par_pt, t, 1); // We compute position and derivative of parameter curve.
      vector<Point> space_pt = surf_->ParamSurface::point(par_pt[0][0], par_pt[0][1], 1);
      der[0] = space_pt[0];
      der[1] = par_pt[1][0]*space_pt[1] + par_pt[1][1]*space_pt[2];
//       der[1].normalize();
  }
}


//===========================================================================
double LiftCurve::start() const
//===========================================================================
{
  return parameter_crv_->startparam();
}


//===========================================================================
double LiftCurve::end() const
//===========================================================================
{
  return parameter_crv_->endparam();
}


//===========================================================================
int LiftCurve::dim() const
//===========================================================================
{
    return 3; // Dimension of lifted curve, not that of the parameter curve.
}


//===========================================================================
bool LiftCurve::approximationOK(double par, Point approxpos,
				  double tol1, double tol2) const
//===========================================================================
{
    // Only first tolerance is used.
  Point pos = eval(par);

  double dist = pos.dist(approxpos); // Both points are in space.

  // @@sbr Currently no input tolerance is used, only the epsgeo_.
  return (dist < epsgeo_);
}
