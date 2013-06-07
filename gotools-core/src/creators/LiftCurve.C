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
