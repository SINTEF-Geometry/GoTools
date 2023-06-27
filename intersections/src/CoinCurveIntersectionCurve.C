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

#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/geometry/SplineCurve.h"

using namespace Go;
using std::vector;

//===========================================================================
shared_ptr<ParamCurve> CoincCurveIntersectionCurve::getCurve() const
//===========================================================================
{
    if (geom_cached_) {
	return cached_geom_curve_;
    }

    ParamObjectInt* obj1 = (ParamObjectInt*)(ipoints_.front()->getObj1());
    ParamCurveInt* cv1 = dynamic_cast<ParamCurveInt*>(obj1);
    if (!cv1)
      return cached_geom_curve_;
    shared_ptr<ParamCurve> cv1_2 = cv1->getParamCurve();
    cached_geom_curve_ =  shared_ptr<ParamCurve>(cv1_2->subCurve(startpar_,
								 endpar_));
    return cached_geom_curve_;
}

//===========================================================================
shared_ptr<ParamCurve>
CoincCurveIntersectionCurve::getParamCurve(int obj_nmb) const
//===========================================================================
{
  if (obj_nmb == 1 && par1_cached_) {
    return cached_param_curve_1_;
  } else if (obj_nmb == 2 && par2_cached_) {
    return cached_param_curve_2_;
  }

  if (obj_nmb != 1 && obj_nmb != 2) {
    throw("Argument to getParamCurve() should be either 1 or 2.");
  }

  if (obj_nmb == 1)
    {
      Point par1(1), par2(1);
      par1.setValue(startpar_);
      par2.setValue(endpar_);
      shared_ptr<ParamCurve> cached_param_curve_1_(new SplineCurve(par1, startpar_,
								   par2, endpar_));
         return cached_param_curve_1_;
    }
  else
    {
      // Note that this parameter curve may be inaccurately parameterized
      // as it assumes corresponding parameterization for both curves
      Point par1(1), par2(1);
      par1.setValue(ipoints_.front()->getPar2()[0]);
      par2.setValue(ipoints_.back()->getPar2()[0]);
      if (ipoints_.front()->getPar1()[0] > ipoints_.back()->getPar1()[0])
	std::swap(par1, par2);
      shared_ptr<ParamCurve> cached_param_curve_2_(new SplineCurve(par1, startpar_,
								   par2, endpar_));
         return cached_param_curve_2_;
    }
}

//===========================================================================
void
CoincCurveIntersectionCurve::evaluateAt(double pval, Point& pos, Point& tan)
//===========================================================================
{
  ParamObjectInt* obj1 = (ParamObjectInt*)(ipoints_.front()->getObj1());
  ParamCurveInt* cv1 = dynamic_cast<ParamCurveInt*>(obj1);
  if (!cv1)
    return;
  shared_ptr<ParamCurve> cv1_2 = cv1->getParamCurve();
  ParamObjectInt* obj2 = (ParamObjectInt*)(ipoints_.front()->getObj2());
  ParamCurveInt* cv2 = dynamic_cast<ParamCurveInt*>(obj2);
  if (!cv2)
    return;
  shared_ptr<ParamCurve> cv2_2 = cv2->getParamCurve();

  vector<Point> der1(2);
    cv1_2->point(der1, pval, 1);

    double par1 = ipoints_.front()->getPar2()[0];
    double par2 = ipoints_.back()->getPar2()[0];
    double seed = par1 + (pval - startpar_)*(par2-par1)/(endpar_ - startpar_);
    if (par1 > par2)
      std::swap(par1, par2);

    double tpar, dist;
    Point close;
    cv2_2->closestPoint(der1[0], par1, par2, tpar, close, dist, &seed);
    pos = 0.5*(der1[0] + close);
    tan = der1[1];
}

//===========================================================================
bool CoincCurveIntersectionCurve::
getGuidePointTangent(shared_ptr<IntersectionPoint> pt, 
		     Point& tan, int type) const
//===========================================================================
{
  return false;  // Doesn't really make sense
}
