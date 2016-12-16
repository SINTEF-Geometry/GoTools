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

#include "GoTools/creators/HermiteApprEvalSurf.h"


#include "GoTools/creators/HermiteApprEvalSurf.h"
#include "GoTools/creators/HermiteGrid2D.h"
#include "GoTools/creators/EvalSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/HermiteInterpolator.h"

using std::vector;
using std::max;
using std::min;

namespace Go
{

HermiteApprEvalSurf::HermiteApprEvalSurf(EvalSurface* sf, double tolerance1,
                                         double tolerance2)
    : surface_(sf), tol1_(tolerance1),tol2_(tolerance2), min_interval_(1.0e-06),
      grid_(*sf, sf->start_u(), sf->end_u(), sf->start_v(), sf->end_v()), method_failed_(false)
//-------------------------------------------------------------------------
// PURPOSE: Constructor
//
// INPUT: sf	- original surface. The surface is NOT copied, only pointer set.
//        tolerance- Required accuracy of approximation.
//
//-------------------------------------------------------------------------
{
  if (tol1_ < min_interval_)
    min_interval_ = tol1_;
}

HermiteApprEvalSurf::HermiteApprEvalSurf(EvalSurface* sf,
                                         double initpars_u[],
                                         int mm,
                                         double initpars_v[],
                                         int nn,
                                         double tolerance1,
                                         double tolerance2)
  : surface_(sf), tol1_(tolerance1), tol2_(tolerance2), min_interval_(0.0001),
    grid_(*sf, initpars_u, initpars_v, mm, nn), method_failed_(false)
//-------------------------------------------------------------------------
// PURPOSE: Constructor
//
// INPUT: sf	- original surface. The surface is NOT copied, only pointer set.
//        initpars - original grid, assumed to be strictly increasing.
//                   and inside domain of "sf"
//	  n     - Number of elements in original grid.
//        tolerance- Required accuracy of approximation.
//-------------------------------------------------------------------------
{
  double start_u = sf->start_u();
  double end_u  = sf->end_u();
  double start_v = sf->start_v();
  double end_v  = sf->end_v();

  int ki;
  for (ki = 1; ki < mm; ++ki)
  {
    if (initpars_u[ki] < start_u || initpars_u[ki] > end_u)
      THROW("Input grid illegal");
  }
  for (ki = 1; ki < nn; ++ki)
  {
    if (initpars_v[ki] < start_v || initpars_v[ki] > end_v)
      THROW("Input grid illegal");
  }

  if (tol1_ < min_interval_)
    min_interval_ = tol1_;
}

void HermiteApprEvalSurf::refineApproximation()
//-------------------------------------------------------------------------
// PURPOSE: Refine initial Hermite grid until interpolant on Hermite grid
//          approximates parametrically original surface within tolerance
//
//-------------------------------------------------------------------------
{
    int ki = 0, kj = 0;

    while (kj < grid_.size2()-1) {
        while (ki < grid_.size1()-1) {
            bool dir_is_u;
            int segment = bisectSegment(ki, kj, dir_is_u);
            if (segment == -1) {
                method_failed_ = true;
                MESSAGE("Method failed, possibly due to small knot interval. "
                        "Tol too strict I guess.");
                return;
            }
        }
    }
}

int HermiteApprEvalSurf::bisectSegment(int left1, int left2, bool& dir_is_u)
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

  int isOK = testSegment(left1, left2, new_knot, dir_is_u);

  if (isOK == 1)
      return (dir_is_u) ? left1 + 1 : left2 + 1;
  else if (isOK == -1)
    return -1;

  grid_.addKnot(*surface_, new_knot, dir_is_u); // @@sbr072009 Tolerance check?

//   // We test tolerance.
//   double spar, epar;
//   Point bezcoef[4];
//   grid_.getSegment(j,j+1,spar,epar,bezcoef);
//   // The new grid point should be within tol1_.

  // Refine new interval to the left of new knot

  isOK = bisectSegment(left1, left2, dir_is_u);

  if (isOK == -1)
    return -1;

  // Refine new interval to the right of new knot

  isOK = bisectSegment(left1, left2, dir_is_u);

  return isOK;
}

int HermiteApprEvalSurf::testSegment(int left1, int left2, double& new_knot, bool& dir_is_u)
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
    double spar1, epar1, spar2, epar2;
    Point bezcoef[4];

    grid_.getSegment(left1, left1 + 1, left2, left2 + 1,
                     spar1, epar1, spar2, epar2, bezcoef);

    double t,t0,t1,t2,t3;

    int numtest = 9;	// Should be an odd number
    double p1, p2;
    double tau1[4], tau2[4];
    int ki, kj, kr, km, kn, ix;
    double upar, vpar;
    const int dim = surface_->dim();

    for (kn=0; kn<numtest; ++kn)
    {
        p2 = (double)(kn+1)/(double)(numtest+1);
        tau2[0]  = (1-p2)*(1-p2);
        tau2[3]  = p2*p2;
        tau2[1]  = 3*tau2[0]*p2;
        tau2[2]  = 3*tau2[3]*(1-p2);
        tau2[0] *= (1-p2);
        tau2[3] *= p2;
        vpar = spar2 + p2*(epar2 - spar2);
        for (km=0; km<numtest; ++km)
        {
            // Calculate position on Bezier segment
            p1 = (double)(km+1)/(double)(numtest+1);
            tau1[0]  = (1-p1)*(1-p1);
            tau1[3]  = p1*p1;
            tau1[1]  = 3*tau1[0]*p1;
            tau1[2]  = 3*tau1[3]*(1-p1);
            tau1[0] *= (1-p1);
            tau1[3] *= p1;

            Point bezval(dim, 0.0);
            for (kj=0, ix=0; kj<4; ++kj)
            {
                vector<double> tmp(dim, 0.0);
                for (ki=0; ki<4; ++ki, ix+=dim)
                    for (kr=0; kr<dim; ++kr)
                        tmp[kr] += bezcoef[ki][kr]*tau1[ki];
                for (kr=0; kr<dim; ++kr)
                    bezval[kr] += tmp[kr]*tau2[kj];
            }

            // Calculate the position on the original surface
            upar = spar1 + p1*(epar1 - spar1);

            // Check quality of approximation point
            bool isOK =  surface_->approximationOK(upar, vpar, bezval,
                                                   tol1_, tol2_);

            if (!isOK)
                break;
        }
        if (km < numtest)
            break;   // Refinement required. Further testing not necessary
    }

    double dom1 = epar1 - spar1;
    double dom2 = epar2 - spar2;
    // We split in the direction with the largest interval.
    dir_is_u = (dom1 > dom2);
    new_knot = (dir_is_u) ? 0.5*(spar1 + epar1) : 0.5*(spar2 + epar2);

    if ((km <= numtest && dom1 < min_interval_) || (kn <= numtest && dom2 < min_interval_))
    {
        MESSAGE("Knot interval too small");
        return -1;  // Do not subdivide any more
    }

    int isOK = ((km > numtest) && (kn > numtest));
    
    return isOK;
}

shared_ptr<SplineSurface> HermiteApprEvalSurf::getSurface()
//----------------------------------------------------------------------
// PURPOSE: Return the cubic spline surface Hermite interpolating the grid.
//
//
//----------------------------------------------------------------------
{
    shared_ptr<SplineSurface> sf;
    if (method_failed_)
        return sf;

#if 0
    MESSAGE("Missing implementation of the interpolator!");
#else

    MESSAGE("Under construction!");  

    // We first interpolate in the second parameter direction, pos & der_v.
    vector<double> coefs_v;
    HermiteInterpolator interpolator_v;
    vector<Point> data = grid_.getData();
    vector<double> param_u = grid_.getKnots(true);
    vector<double> param_v = grid_.getKnots(false);
    interpolator_v.interpolate(data, param_v, coefs_v);


    vector<double> coefs;
    HermiteInterpolator interpolator_u;
    interpolator_u.interpolate(data, param_u, coefs);

  
    // Then der_u & der_uv.

    // Transpose the coefficients of pos and der_v.
  
    // Transpose the coefficients of der_u and der_uv

    // We then interpolate in the first parameter direction.

    // Transpose the coefficients.
    BsplineBasis basis_v = interpolator_v.basis();
    BsplineBasis basis_u = interpolator_u.basis();
    const int dim = surface_->dim();
    sf = (shared_ptr<SplineSurface>)(new SplineSurface(basis_u, basis_v, coefs.begin(), dim));
  
#endif
  
    return sf;
}

}
