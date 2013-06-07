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

#include "GoTools/creators/CrossTangentOffset.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;

//===========================================================================

CrossTangentOffset::CrossTangentOffset(shared_ptr<SplineCurve>& poscurve,
				       shared_ptr<SplineCurve>& tangcv1,
				       shared_ptr<SplineCurve>& tangcv2,
				       shared_ptr<SplineCurve>& blend1,
				       shared_ptr<SplineCurve>& blend2)
    : poscurve_(poscurve)
//===========================================================================
{
  // Test input
  ALWAYS_ERROR_IF(poscurve.get() == 0 || tangcv1.get() == 0 || 
	      tangcv2.get() == 0 || blend1.get() == 0 || blend2.get() == 0, 
	      "Missing curve");
  ALWAYS_ERROR_IF(poscurve->dimension() != 3 ||
	      poscurve->dimension() != tangcv1->dimension() ||
	      poscurve->dimension() != tangcv2->dimension(),
	      "Dimension mismatch");
  ALWAYS_ERROR_IF(blend1->dimension() != 1 || blend2->dimension() != 1,
		  "Blending function of dimension different from 1");


  //poscurve_ = poscurve;
  tangcurves_.push_back(tangcv1);
  tangcurves_.push_back(tangcv2);
  blends_.push_back(blend1);
  blends_.push_back(blend2);

  // Create linear length function.
  double par1 = poscurve_->startparam();
  double par2 = poscurve_->endparam();

  // Compute length of cross tangent curve in the endpoints.
  Point cross1 = evalcrtan(par1);
  Point cross2 = evalcrtan(par2);
  Point l1(1), l2(1);
  l1.setValue(cross1.length());
  l2.setValue(cross2.length());
  length_ = (shared_ptr<SplineCurve>)(new SplineCurve(l1, l2));
  length_->setParameterInterval(par1, par2);
}

//===========================================================================

CrossTangentOffset::~CrossTangentOffset()
//===========================================================================
{
}

//===========================================================================

Point CrossTangentOffset::eval( double t) const 
//===========================================================================
{
  int dim = poscurve_->dimension();
  Point pos(dim), cross(dim), len(1);

  poscurve_->point(pos, t);
  cross = evalcrtan(t);
  cross.normalize();
  length_->point(len,t);

  return pos + cross*len[0];
}

//===========================================================================

void CrossTangentOffset::eval(double t, int n, Point der[]) const 
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  if (n == 0)
    der[0] = eval(t);
  else
    {
      int dim = poscurve_->dimension();
      der[0].resize(dim);
      der[1].resize(dim);
      der[0].setValue(0.0);
      der[1].setValue(0.0);
      vector<Point> pos(2), len(2);
      Point cross[2];
      Point crn(dim), crn1(dim);
      poscurve_->point(pos, t, 1);
      evalcrtan(t, 1, cross);
      double length = cross[0].length();
      crn = cross[0]/length;
      crn1 = cross[1]/length;      
      length_->point(len, t, 1);

      der[0] = pos[0] + crn*len[0][0];
      der[1] = pos[1] + crn*len[1][0] + crn1*len[0][0] - 
	crn*(crn*crn1)*len[0][0];
    }
}


//===========================================================================

double CrossTangentOffset::start() const 
//===========================================================================
{
  return poscurve_->startparam();
}


//===========================================================================

double CrossTangentOffset::end() const 
//===========================================================================
{
  return poscurve_->endparam();
}

//===========================================================================

int CrossTangentOffset::dim() const
//===========================================================================
{
  return poscurve_->dimension();
}

//===========================================================================

bool CrossTangentOffset::approximationOK(double par, Point approxpos,
					  double tol1, double tol2) const
//===========================================================================
{
  double tol3 = 0.000001*tol1;
  Point pos = eval(par);
  double dist = pos.dist(approxpos);  // Distance between original
                                      // curve and approximation
  if (dist > tol1)
    return false;   // Approximation not good enough

  if (dist < tol3)
    return true;

  Point diff = evalcrtan(par);
  Point pt1, pt2;
  tangcurves_[0]->point(pt1, par);
  tangcurves_[1]->point(pt2, par);
  Point normal = pt1 % pt2;
  double ang = normal.angle(diff);
  double pihalf = 3.141592653589793/2.0;
  if (fabs(pihalf-ang) > tol2)
    return false;   // Approximated cross tangent not in tangent plane

  return true;
}


//===========================================================================

Point CrossTangentOffset::evalcrtan( double t) const
//===========================================================================
{
  int dim = poscurve_->dimension();
  Point point1(dim), point2(dim);
  Point fac;
  point1.setValue(0.0);

  for (size_t ki=0; ki<tangcurves_.size(); ki++)
    {
      tangcurves_[ki]->point(point2, t);
      blends_[ki]->point(fac, t);
      point1 += fac[0]*point2;
    }

  return point1;
}

//===========================================================================
void CrossTangentOffset::evalcrtan(double t, int n, Point der[]) const
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  if (n == 0)
    der[0] = evalcrtan(t);
  else
    {
      int dim = poscurve_->dimension();
      der[0].resize(dim);
      der[1].resize(dim);
      der[0].setValue(0.0);
      der[1].setValue(0.0);
      vector<Point> pt(2), fac(2);
      for (size_t ki=0; ki<tangcurves_.size(); ki++)
	{
	  tangcurves_[ki]->point(pt, t, 1);
	  blends_[ki]->point(fac, t, 1);
	  der[0] += fac[0][0]*pt[0];
	  der[1] += (fac[1][0]*pt[0] + fac[0][0]*pt[1]);
	}
    }
}


