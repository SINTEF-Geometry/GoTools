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
