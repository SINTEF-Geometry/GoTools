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

#include "GoTools/creators/HermiteAppC.h"
#include "GoTools/creators/HermiteGrid1D.h"
#include "GoTools/creators/EvalCurve.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/HermiteInterpolator.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;


HermiteAppC::HermiteAppC(EvalCurve* crv, double tolerance1,
			     double tolerance2)
  : curve_(crv), tol1_(tolerance1),tol2_(tolerance2), min_interval_(0.0001),
    grid_(*crv, crv->start(), crv->end()), method_failed_(false)
//-------------------------------------------------------------------------
// PURPOSE: Constructor
//
// INPUT: crv	- original curve. The curve is NOT copied, only pointer set.
//        tolerance- Required accuracy of approximation.
//
//-------------------------------------------------------------------------
{
  if (tol1_ < min_interval_)
    min_interval_ = tol1_;
}

HermiteAppC::HermiteAppC(EvalCurve* crv,
			     double initpars[],
			     int n,
			     double tolerance1,
			     double tolerance2)
  : curve_(crv), tol1_(tolerance1), tol2_(tolerance2), min_interval_(0.0001),
    grid_(*crv, initpars, n), method_failed_(false)
//-------------------------------------------------------------------------
// PURPOSE: Constructor
//
// INPUT: crv	- original curve. The curve is NOT copied, only pointer set.
//        initpars - original grid, assumed to be strictly increasing.
//                   and inside domain of "crv"
//	  n     - Number of elements in original grid.
//        tolerance- Required accuracy of approximation.
//-------------------------------------------------------------------------
{
  double start = crv->start();
  double end  = crv->end();

  int i;
  for (i=1; i<n; i++)
    if (initpars[i] < start || initpars[i] > end)
      THROW("Input grid illegal");

  if (tol1_ < min_interval_)
    min_interval_ = tol1_;
}

void HermiteAppC::refineApproximation()
//-------------------------------------------------------------------------
// PURPOSE: Refine initial Hermite grid until interpolant on Hermite grid
//          approximates paramatrically original curve within tolerance
//
//-------------------------------------------------------------------------
{
  int j=0;

  while (j < grid_.size()-1) {
      j = bisectSegment(j);
      if (j == -1) {
	  method_failed_ = true;
	  MESSAGE("Method failed, possibly due to small knot interval. "
		  "Tol too strict I guess.");
	  return;
      }
  }
}

int HermiteAppC::bisectSegment(int j)
//------------------------------------------------------------------------
//   Check the Hermite curve interpolant over the segment
//   (g->_knots[jj],g->_knots[jj+1]).
//   If the interpolant violates tolerance, the Hermite data grid is refined.
//   This is done by adding more grid parameters, i.e
//   g->_knots is extended.
//
// INPUT: j  - Index of grid parameter determining start point of segment
//             before refinement.
//
// OUTPUT:bisectSegment()  - Index of grid parameter determining end point of
//                           segment after refinement.
//
//
//------------------------------------------------------------------------
{

  double new_knot;

  int isOK = testSegment(j,new_knot);

  if (isOK == 1)
    return j+1;
  else if (isOK == -1)
    return -1;

  grid_.addKnot(*curve_,new_knot); // @@sbr072009 Tolerance check?

//   // We test tolerance.
//   double spar, epar;
//   Point bezcoef[4];
//   grid_.getSegment(j,j+1,spar,epar,bezcoef);
//   // The new grid point should be within tol1_.

  // Refine new interval to the left of new knot

  j = bisectSegment(j);

  if (j == -1)
    return -1;

  // Refine new interval to the right of new knot

  j = bisectSegment(j);

  return j;
}

int HermiteAppC::testSegment(int j, double& new_knot)
//-------------------------------------------------------------------
// PURPOSE: Calculate distance from segment to evaluator.
//          Find the parameter of greatest deviation (among a sample)
//	    If segment is too small, the warning flag is set and
//          distance 0 is returned.
//
//
//
//--------------------------------------------------------------------
{
  double spar,epar;
  Point bezcoef[4];

  grid_.getSegment(j,j+1,spar,epar,bezcoef);

  double t,t0,t1,t2,t3;

  int numtest = 9;	// Should be an odd number
  int n;
  for (n = 1; n <= numtest; n++)
  {
    t = (double)n/(double)(numtest+1);

    // Calculate position on Bezier segment

    t0  = (1-t)*(1-t);  t3  = t*t;
    t1  = 3*t0*t;     	t2  = 3*t3*(1-t);
    t0 *= (1-t);	t3 *= t;

    Point bezval = bezcoef[0]*t0 + bezcoef[1]*t1 +bezcoef[2]*t2 +
      bezcoef[3]*t3;

    // Check quality of approximation point
    bool isOK =  curve_->approximationOK(spar + t*(epar-spar), bezval,
					 tol1_, tol2_);

    if (!isOK)
      break;
  }

  new_knot = 0.5*(spar + epar);

  if (n <= numtest && epar-spar < min_interval_)
    {
      MESSAGE("Knot interval too small");
      return -1;  // Do not subdivide any more
    }

  return (n > numtest);
}

shared_ptr<SplineCurve> HermiteAppC::getCurve()
//----------------------------------------------------------------------
// PURPOSE: Return the cubic spline curve Hermite interpolating the grid.
//
//
//----------------------------------------------------------------------
{
  shared_ptr<SplineCurve> cv;
  if (method_failed_)
      return cv;

  vector<double> coefs;
  HermiteInterpolator interpolator;
  vector<Point> data = grid_.getData();
  vector<double> param = grid_.getKnots();
  interpolator.interpolate(data, param, coefs);

  BsplineBasis basis = interpolator.basis();

  cv = (shared_ptr<SplineCurve>)(new SplineCurve(basis.numCoefs(), 
						     basis.order(),
						     basis.begin(), 
						     &coefs[0], 
						     grid_.dim(), false));

  return cv;
}

