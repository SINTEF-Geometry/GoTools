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

#include "GoTools/compositemodel/HermiteApprEvalSurf.h"


#include "GoTools/compositemodel/HermiteApprEvalSurf.h"
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
    int kj = 0;

    bool debug_mode = false;
    if (debug_mode)
    {
        MESSAGE("In debug mode!");
    }
    
    while (kj < grid_.size2()-1) {
        int ki = 0;
        while (ki < grid_.size1()-1) {
            bool dir_is_u;
            int segment = bisectSegment(ki, kj, dir_is_u);
            if (debug_mode && ((grid_.size1() > 2) || (grid_.size2() > 2)))
            {
                MESSAGE("Exiting early!");
                return;
            }
            if (segment == -1) {
                method_failed_ = true;
                MESSAGE("Method failed, possibly due to small knot interval. "
                        "Tol too strict I guess.");
                return;
            }

            ++ki;
        }

        ++kj;
    }

    return;
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

  // If isOK == 0 we should refine, in the direction with the largest knot span.
  int isOK = testSegment(left1, left2, new_knot, dir_is_u);
  if (isOK == 1)
  {
      return (dir_is_u) ? left1 + 1 : left2 + 1;
  }
  else if (isOK == -1)
  {
    return -1;
  }

  grid_.addKnot(*surface_, new_knot, dir_is_u); // @@sbr072009 Tolerance check?

  bool debug_mode = false;
  if (debug_mode && ((grid_.size1() > 3) || (grid_.size2() > 3)))
  {
      MESSAGE("Exiting early!");
      return (dir_is_u) ? left1 + 1 : left2 + 1;
  }

//   // We test tolerance.
//   double spar, epar;
//   Point bezcoef[4];
//   grid_.getSegment(j,j+1,spar,epar,bezcoef);
//   // The new grid point should be within tol1_.

  // Refine new interval to the left of new knot

  int grid_ind = bisectSegment(left1, left2, dir_is_u);

  if (grid_ind == -1)
    return -1;

  // Refine new interval to the right of new knot

  grid_ind = bisectSegment(left1, left2, dir_is_u);

  return grid_ind;
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
    Point bezcoef[16];

    grid_.getSegment(left1, left1 + 1, left2, left2 + 1,
                     spar1, epar1, spar2, epar2, bezcoef);

//    double t,t0,t1,t2,t3;

    int numtest = 9;	// Should be an odd number
    double p1, p2;
    double tau1[4], tau2[4];
    int ki, kj, km, kn;
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

            Point bezval(dim);
            bezval.setValue(0.0);
            for (kj=0; kj<4; ++kj)
            {
                Point tmp(dim);
                tmp.setValue(0.0);
                for (ki=0; ki<4; ++ki)
                {
                    //for (kr=0; kr<dim; ++kr)
                    tmp += bezcoef[kj*4+ki]*tau1[ki];
                }
                //for (kr=0; kr<dim; ++kr)
                bezval += tmp*tau2[kj];
            }

            // Calculate the position on the original surface
            upar = spar1 + p1*(epar1 - spar1);

            // Check quality of approximation point
            // tol1_ is used as tolerance in geometry space.
            // tol2_ is currently not used (as of 2017/01/05).
            bool apprOK =  surface_->approximationOK(upar, vpar, bezval,
                                                   tol1_, tol2_);

            if (!apprOK)
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
        std::cout << "dom1: " << dom1 << ", dom2: " << dom2 << ", spar1: " << spar1 << ", spar2: " << spar2 << std::endl;
        MESSAGE("Knot interval too small");
        method_failed_ = true;
        return -1;  // Do not subdivide any more
    }

    int isOK = ((km == numtest) && (kn == numtest)) ? 1 : 0;
    // if (isOK == 0)
    // {
    //     std::cout << "Not ok! km: " << km << ", kn: " << kn << ", dom1: " << dom1 << ", dom2: " << dom2 <<
    //         ", upar: " << upar << ", vpar: " << vpar << std::endl;
    // }

    return isOK;
}

shared_ptr<SplineSurface> HermiteApprEvalSurf::getSurface(bool& method_failed)
//----------------------------------------------------------------------
// PURPOSE: Return the cubic spline surface Hermite interpolating the grid.
//
//
//----------------------------------------------------------------------
{
    method_failed = method_failed_;

    // We extract the data used by the interpolator.
    // We use the version with array of double's (as opposed to Point's).
    const int mm = grid_.size1();
    const int nn = grid_.size2();
    const int dim = surface_->dim();
    vector<double> pos_der_v, der_u_der_uv;
    pos_der_v.reserve(mm*nn*dim*2);
    der_u_der_uv.reserve(mm*nn*dim*2);
    vector<Point> data = grid_.getData();
    for (int kj = 0; kj < nn; ++kj)
    {
        for (int ki = 0; ki < mm; ++ki)
        {
            pos_der_v.insert(pos_der_v.end(),
                             data[4*(kj*mm+ki)].begin(), data[4*(kj*mm+ki)].end());
            der_u_der_uv.insert(der_u_der_uv.end(),
                                data[4*(kj*mm+ki)+1].begin(), data[4*(kj*mm+ki)+1].end());
        }
        for (int ki = 0; ki < mm; ++ki)
        {
            pos_der_v.insert(pos_der_v.end(),
                             data[4*(kj*mm+ki)+2].begin(), data[4*(kj*mm+ki)+2].end());
            der_u_der_uv.insert(der_u_der_uv.end(),
                                data[4*(kj*mm+ki)+3].begin(), data[4*(kj*mm+ki)+3].end());
        }
    }

    vector<double> coefs_pos_der_v;
    vector<double> param_v = grid_.getKnots(false);
    vector<double> param_v2(param_v.size()*2);
    for (size_t ki = 0; ki < param_v.size(); ++ki)
    {
        param_v2[2*ki] = param_v[ki];
        param_v2[2*ki+1] = param_v[ki];
    }

    // We first interpolate in the second parameter direction, pos & der_v.
    HermiteInterpolator interpolator;
    // The number of points includes the tangents.
    interpolator.interpolate(2*nn, mm*dim, &param_v2[0], &pos_der_v[0], coefs_pos_der_v);

#if 0
    std::cout << "Result from interpolating pos_der_v:" << std::endl;
    for (size_t ki = 0; ki < coefs_pos_der_v.size(); ki+= dim)
    {
        std::cout << coefs_pos_der_v[ki] << " " << coefs_pos_der_v[ki+1] << " " << coefs_pos_der_v[ki+2] << std::endl;
    }
#endif
    
    vector<double> coefs_der_u_der_uv;
    // Then der_u & der_uv.
    interpolator.interpolate(2*nn, mm*dim, &param_v2[0], &der_u_der_uv[0], coefs_der_u_der_uv);
    BsplineBasis basis_v = interpolator.basis();

#if 0
    std::cout << "Result from interpolating der_u_der_uv:" << std::endl;
    for (size_t ki = 0; ki < coefs_der_u_der_uv.size(); ki+= dim)
    {
        std::cout << coefs_der_u_der_uv[ki] << " " << coefs_der_u_der_uv[ki+1] << " " << coefs_der_u_der_uv[ki+2] << std::endl;
    }
#endif
    
    // Transpose coefs_pos_der_v & coefs_der_u_der_uv.
    // vector<double> coefs_pos_der_v_tr;
    // coefs_pos_der_v_tr.reserve(coefs_pos_der_v.size());
    // vector<double> coefs_der_u_der_uv_tr;
    // coefs_der_u_der_uv_tr.reserve(coefs_der_u_der_uv.size());
    vector<double> coefs_pre_u_int;
    for (int ki = 0; ki < mm; ++ki)
    {
        for (int kj = 0; kj < 2*nn; ++kj)
        {
//            size_t ind = 2*(ki*mm + kj);
            coefs_pre_u_int.insert(coefs_pre_u_int.end(),
                                   coefs_pos_der_v.begin() + (kj*mm + ki)*dim,
                                   coefs_pos_der_v.begin() + (kj*mm + ki + 1)*dim);
        }
        for (int kj = 0; kj < 2*nn; ++kj)
        {
//            size_t ind = 2*(ki*mm + kj);
            coefs_pre_u_int.insert(coefs_pre_u_int.end(),
                                   coefs_der_u_der_uv.begin() + (kj*mm + ki)*dim,
                                   coefs_der_u_der_uv.begin() + (kj*mm + ki + 1)*dim);
        }
    }

#if 0
    std::cout << "coefs_pre_u_int.size(): " << coefs_pre_u_int.size() << std::endl;
#endif
    
    vector<double> sf_coefs;
    vector<double> param_u = grid_.getKnots(true);
    vector<double> param_u2(param_u.size()*2);
    for (size_t ki = 0; ki < param_u.size(); ++ki)
    {
        param_u2[2*ki] = param_u[ki];
        param_u2[2*ki+1] = param_u[ki];
    }
    // We then interpolate in the first parameter direction.
    interpolator.interpolate(2*mm, 2*nn*dim, &param_u2[0], &coefs_pre_u_int[0], sf_coefs);
    BsplineBasis basis_u = interpolator.basis();

    // We then transpose the sf coefs.
    vector<double> sf_coefs_tr(sf_coefs.size());
    for (int ki = 0; ki < 2*mm; ++ki)
    {
        for (int kj = 0; kj < 2*nn; ++kj)
        {
            for (int kk = 0; kk < dim; ++kk)
            {
                sf_coefs_tr[(kj*2*mm+ki)*dim+kk] = sf_coefs[(ki*2*nn+kj)*dim+kk];
            }
        }
    }

#if 0
    std::cout << "sf_coefs:\n" << std::endl;
    for (size_t ki = 0; ki < sf_coefs.size()/dim; ++ki)
    {
        std::cout << "ki = " << ki << ": " << sf_coefs[ki*3] << " " << sf_coefs[ki*3+1] << " " << sf_coefs[ki*3+2] << std::endl;  
        if ((ki+1)%(2*nn) == 0)
        {
            std::cout << std::endl;
        }
    }

    std::cout << "\nsf_coefs_tr:\n" << std::endl;
    for (size_t ki = 0; ki < sf_coefs_tr.size()/dim; ++ki)
    {
        std::cout << "ki = " << ki << ": " << sf_coefs_tr[ki*3] << " " << sf_coefs_tr[ki*3+1] << " " <<
            sf_coefs_tr[ki*3+2] << std::endl;
        if ((ki+1)%(2*mm) == 0)
        {
            std::cout << std::endl;
        }
    }
#endif
    
    // Create the Hermite interpolating surface.
    shared_ptr<SplineSurface> sf(new SplineSurface(basis_u, basis_v, sf_coefs_tr.begin(), dim));

    std::cout << "num_coefs_u: " << sf->numCoefs_u() << ", num_coefs_v: " << sf->numCoefs_v() << std::endl;

    return sf;
}

}
