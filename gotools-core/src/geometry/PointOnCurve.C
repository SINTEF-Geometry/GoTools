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

#include "GoTools/geometry/PointOnCurve.h"

using namespace Go;
using std::vector;

//===========================================================================
PointOnCurve::PointOnCurve()
//===========================================================================
{
  t1_ = t2_ = par_ = 0.0;
}

//===========================================================================
PointOnCurve::PointOnCurve(shared_ptr<ParamCurve> curve, double par)
    : par_(par), crv_(curve)
//===========================================================================
{
    t1_ = crv_->startparam();
    t2_ = crv_->endparam();

    point_ = crv_->point(par_);
}

//===========================================================================
PointOnCurve::PointOnCurve(shared_ptr<ParamCurve> curve, Point pnt)
    : point_(pnt), crv_(curve)
//===========================================================================
{

    t1_ = crv_->startparam();
    t2_ = crv_->endparam();

    double clo_dist;
    Point clo_pt;
    crv_->closestPoint(point_, t1_, t2_, par_, clo_pt, clo_dist);
}

//===========================================================================
PointOnCurve::~PointOnCurve()
//===========================================================================
{
}

//===========================================================================
Point PointOnCurve::getPos() const 
//===========================================================================
{
    return point_;
}

//===========================================================================
void PointOnCurve::evaluate(int der, std::vector<Point>& deriv) const
//===========================================================================
{
  if (crv_.get())
    crv_->point(deriv, par_, der);
}

//===========================================================================
void PointOnCurve::setParInterval(double start, double end)
//===========================================================================
{
    ASSERT(start <= par_ && par_ <= end);
    t1_ = start;
    t2_ = end;
}


