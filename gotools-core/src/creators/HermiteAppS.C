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

#include "GoTools/creators/HermiteAppS.h"
#include "GoTools/creators/HermiteGrid1D.h"
#include "GoTools/creators/EvalCurveSet.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/HermiteInterpolator.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;

HermiteAppS::HermiteAppS(EvalCurveSet* surface, double tolerance1,
			     double tolerance2, std::vector<int> dims)
  : surface_(surface), grid_(*surface_, surface_->start(), surface_->end(), dims),
    tol1_(tolerance1),tol2_(tolerance2), min_interval_(0.0001)
//-------------------------------------------------------------------------
// PURPOSE: Constructor
//
// INPUT: surface	- original surface. The surface are NOT copied, only pointer set.
//        tolerance- Required accuracy of approximation.
//
//-------------------------------------------------------------------------
{
    // First we check that the parameter spaces correspond.
//     double from = surface_->start();
//     double to = surface->end();
//     double equal_tol = 1e-10;
}

HermiteAppS::HermiteAppS(EvalCurveSet* surface,
			 double initpars[],
			 int n,
			 double tolerance1,
			 double tolerance2,
			 std::vector<int> dims)
  : surface_(surface), grid_(*surface, initpars, n, dims),
    tol1_(tolerance1), tol2_(tolerance2), min_interval_(0.01)
//-------------------------------------------------------------------------
// PURPOSE: Constructor
//
// INPUT: surf	- original surface. The surface is NOT copied, only pointer set.
//        initpars - original grid, assumed to be strictly increasing.
//                   and inside domain of "surf"
//	  n     - Number of elements in original grid.
//        tolerance- Required accuracy of approximation.
//-------------------------------------------------------------------------
{
    // We check that initpars is ascending.
    for (int ki = 0; ki < n - 1; ++ki)
	if (initpars[ki] >= initpars[ki+1])
	    THROW("Input grid illegal");

    if ((initpars[0] < surface_->start()) ||
	(initpars[n-1] > surface_->end()))
	THROW("Input grid illegal");
}

void HermiteAppS::refineApproximation()
//-------------------------------------------------------------------------
// PURPOSE: Refine initial Hermite grid until interpolant on Hermite grid
//          approximates paramatrically original surface within tolerance
//
//-------------------------------------------------------------------------
{
  int j=0;

  // Reset approximation errors if stored
  surface_->resetErr();

  while (j < grid_.size()-1)
      j = bisectSegment(j);
}

int HermiteAppS::bisectSegment(int j)
//------------------------------------------------------------------------
//   Check the Hermite surface interpolant over the segment
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

  bool isOK = testSegment(j,new_knot);

  if (isOK)
    return j+1;

//   for (int ki = 0; ki < grids_.size(); ++ki)
  grid_.addKnot(*surface_, new_knot);

  // Reset approximation errors if stored
  surface_->resetErr();

  // Refine new interval to the left of new knot

  j = bisectSegment(j);

  // Refine new interval to the right of new knot

  j = bisectSegment(j);

  return j;
}

bool HermiteAppS::testSegment(int j, double& new_knot)
//-------------------------------------------------------------------
// PURPOSE: Calculate distance from segment to original surface.
//          Find the parameter of greatest deviation (among a sample)
//	    If segment is to small, the warning flag is set and
//          distance 0 is returned.
//
//
//
//--------------------------------------------------------------------
{
    // We must return parameter of greatest deviation. But such a deviation may be
    // more essential to for instance a cross tangent curve than for a pos curve...
    // We let new knot be average value of new_knot in each grid.

//     double avg_value = 0.0;
//     int nmb_cvs = grids_.size();

//     for (ki = 0; ki < grids_.size(); ++ki) {

    double spar,epar;
    vector<vector<Point> > bezcoef; //[4];

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

	    vector<Point> bezvals;
	    for (size_t ki = 0; ki < bezcoef.size(); ++ki)
		bezvals.push_back(bezcoef[ki][0]*t0 + bezcoef[ki][1]*t1 +
				  bezcoef[ki][2]*t2 + bezcoef[ki][3]*t3);

	    // Check quality of approximation point
	    bool isOK =  surface_->approximationOK(spar + t*(epar-spar), bezvals,
						   tol1_, tol2_);

	    if (!isOK)
		break;
	}

    new_knot = 0.5*(spar + epar);
//     avg_value += 0.5*(spar + epar)/nmb_cvs;

    if (n <= numtest && epar-spar < min_interval_)
	{
	    MESSAGE("Knot interval too small");
	    return 1;  // Do not subdivide any more
	}

    return (n > numtest);
}

vector<shared_ptr<SplineCurve> > HermiteAppS::getCurves()
//----------------------------------------------------------------------
// PURPOSE: Return the cubic spline curves Hermite interpolating the grid.
//
//
//----------------------------------------------------------------------
{
    vector<shared_ptr<SplineCurve> > return_cvs;

//     for (int ki = 0; ki < grids_.size(); ++ki) {
    vector<vector<Point> > data = grid_.getData();
    vector<double> param = grid_.getKnots();

    for (size_t ki = 0; ki < data.size(); ++ki) {
	HermiteInterpolator interpolator;
	vector<double> coefs;
	interpolator.interpolate(data[ki], param, coefs);

	BsplineBasis basis = interpolator.basis();

	shared_ptr<SplineCurve> cv(new SplineCurve(basis.numCoefs(), 
						       basis.order(),
						       basis.begin(), 
						       &coefs[0],
						   grid_.dims((int)ki), false));
	return_cvs.push_back(cv);
    }

    return return_cvs;
}

