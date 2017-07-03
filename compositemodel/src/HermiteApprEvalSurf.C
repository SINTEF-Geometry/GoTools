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
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"

#include <fstream>

using std::vector;
using std::max;
using std::min;

namespace Go
{

double avgAngle(Point bezcoef[16], bool dir_is_u);

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
  : surface_(sf), tol1_(tolerance1), tol2_(tolerance2), min_interval_(1.0e-06),
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

    const bool allow_failure = (no_split_cvs_2d_.size() > 0);    
    const bool debug_mode = false;
    if (debug_mode)
    {
        MESSAGE("In debug mode!");
    }
    
    while (kj < grid_.size2()-1) {
        int ki = 0;
        while (ki < grid_.size1()-1) {
            bool dir_is_u;
            int segment = bisectSegment(ki, kj, dir_is_u);
            //std::cout << "DEBUG: refineApproxiMation(): dir_is_u: " << dir_is_u << std::endl;
            if (debug_mode && ((grid_.size1() > 2) || (grid_.size2() > 2)))
            {
                MESSAGE("Exiting early!");
                return;
            }
            if ((segment == -1) && (!allow_failure)) {
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

  //std::cout << "DEBUG: bisectSegment(): dir_is_u: " << dir_is_u << std::endl;
    
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

    int elem_status = grid_.getNoSplitStatus(left1, left2);
    if (elem_status == 3)
    {
        // vector<double> knots_u = grid_.getKnots(true);
        // vector<double> knots_v = grid_.getKnots(false);
        MESSAGE("We should not split any further. spar1 = " << spar1 << ", epar1: " << epar1 <<
                ", spar2: " << spar2 << ", epar2: " << epar2);
        return -1;
    }

#if 1
    // We write to file the bezier coefs.
    std::ofstream fileout("tmp/bez_coefs.g2");
    vector<double> pts_data;
    pts_data.reserve(16*3);
    for (int ki = 0; ki < 16; ++ki) {
        pts_data.insert(pts_data.end(), bezcoef[ki].begin(), bezcoef[ki].end());
    }
    PointCloud3D pt_cloud(pts_data.begin(), pts_data.size()/surface_->dim());
    pt_cloud.writeStandardHeader(fileout);
    pt_cloud.write(fileout);
#endif
    
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
    
    int isOK = ((km == numtest) && (kn == numtest)) ? 1 : 0;
    // if (isOK == 0)
    // {
    //     std::cout << "Not ok! km: " << km << ", kn: " << kn << ", dom1: " << dom1 << ", dom2: " << dom2 <<
    //         ", upar: " << upar << ", vpar: " << vpar << std::endl;
    // }
    if (isOK)
    {
        return isOK;
    }
    else
    {
        const double dom1 = epar1 - spar1;
        const double dom2 = epar2 - spar2;

        const double max_dom = std::max(dom1, dom2);
        
        int split_status = splitDomain(spar1, epar1, spar2, epar2, bezcoef,
                                       dir_is_u, new_knot);

        if (split_status != 0)
        {
//            MESSAGE("We did not split! We should update the grid element with a no split-status!");
            if (split_status == 3)
            {
//                MESSAGE("Split failed even though we have not split in both directions! Something wrong going on!");
                grid_.setNoSplitStatus(left1, left2, split_status);
                return -1;
            }
            else
            {
                int elem_status = grid_.getNoSplitStatus(left1, left2);
                if ((elem_status == 3) || (elem_status == split_status))
                { // This means that the element is already marked as "not split", hence we should not have gotten this far.
                    MESSAGE("The element status implies that we should not have split in this direction! Fix!");
                    return -1;
                }
                else
                {
                    int new_elem_status = elem_status + split_status;
                    grid_.setNoSplitStatus(left1, left2, new_elem_status);
                    
                    // @@sbr201706 We should consider adding yet another flag to the elem to denote approximation failure.
                    return 0;//-1;
                }
            }
        }
    
        if ((km <= numtest && max_dom < min_interval_) || (kn <= numtest && max_dom < min_interval_))
        {
            MESSAGE("INFO: dom1: " << dom1 << ", dom2: " << dom2 << ", spar1: " << spar1 << ", spar2: " << spar2);
            MESSAGE("Knot interval too small: max_dom = " << max_dom << ", min_interval_ = " << min_interval_);
            method_failed_ = true;
            return -1;  // Do not subdivide any more
        }
        else
        {
            return 0;
        }
    }
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
    vector<double> knots_u = grid_.getKnots(true);
    vector<double> knots_v = grid_.getKnots(false);
    const double knot_tol = 1.0e-14; // We require exact match, i.e. close to double precision.
    const int mm = grid_.size1();
    const int nn = grid_.size2();
    const int mm_red = mm - removed_grid_u_.size();
    const int nn_red = nn - removed_grid_v_.size();
    MESSAGE("INFO: mm: " << mm << ", mm_red: " << mm_red << ", nn: " << nn << ", nn_red: " << nn_red);
    const int dim = surface_->dim();
    vector<double> pos_der_v, der_u_der_uv;
    pos_der_v.reserve(mm_red*nn_red*dim*2);
    der_u_der_uv.reserve(mm_red*nn_red*dim*2);
    vector<Point> data = grid_.getData();

#ifndef NDEBUG
    {
        std::ofstream fileout_debug("tmp/grid_data_u.g2");
        std::ofstream fileout_debug2("tmp/grid_data_v.g2");
        std::ofstream fileout_debug3("tmp/grid_pts.g2");
        vector<double> line_pts_u, line_pts_v, grid_pts;
        int num_nodes = data.size()/4;
        line_pts_u.reserve(2*num_nodes*3);
        line_pts_v.reserve(2*num_nodes*3);
        const double scaling_factor = 0.1;
        for (size_t ki = 0; ki < data.size(); ki += 4)
        {
            Point pos = data[ki];
            Point end_u = pos + data[ki+1]*scaling_factor;
            Point end_v = pos + data[ki+2]*scaling_factor;
            Point end_uv = pos + data[ki+3];
            line_pts_u.insert(line_pts_u.end(), pos.begin(), pos.end());
            line_pts_u.insert(line_pts_u.end(), end_u.begin(), end_u.end());
            line_pts_v.insert(line_pts_v.end(), pos.begin(), pos.end());
            line_pts_v.insert(line_pts_v.end(), end_v.begin(), end_v.end());
            // line_pts.insert(line_pts.end(), pos.begin(), pos.end());
            // line_pts.insert(line_pts.end(), end_uv.begin(), end_uv.end());
            int ind = ki/4;
            int ind_u = ind%mm;
            int ind_v = ind/mm;
            if (ind_u < mm - 1)
            {
                grid_pts.insert(grid_pts.end(), pos.begin(), pos.end());
                Point next_pos = data[ki+4];
                grid_pts.insert(grid_pts.end(), next_pos.begin(), next_pos.end());
            }
            if (ind_v < nn - 1)
            {
                grid_pts.insert(grid_pts.end(), pos.begin(), pos.end());
                Point next_pos = data[ki+4*mm];
                grid_pts.insert(grid_pts.end(), next_pos.begin(), next_pos.end());
            }            
        }
	LineCloud line_cloud(&line_pts_u[0], num_nodes);
        line_cloud.writeStandardHeader(fileout_debug);
        line_cloud.write(fileout_debug);
	LineCloud line_cloud2(&line_pts_v[0], num_nodes);
        line_cloud2.writeStandardHeader(fileout_debug2);
        line_cloud2.write(fileout_debug2);
	LineCloud line_cloud3(&grid_pts[0], grid_pts.size()/6);
        line_cloud3.writeStandardHeader(fileout_debug3);
        line_cloud3.write(fileout_debug3);
    }
#endif
    
    auto iter_v = removed_grid_v_.begin();
    for (int kj = 0; kj < nn; ++kj)
    {
        while ((iter_v != removed_grid_v_.end()) && (*iter_v + knot_tol < knots_v[kj]))
        {
            ++iter_v;
        }
        if ((iter_v != removed_grid_v_.end()) && (*iter_v - knots_v[kj] < knot_tol))
        {
            continue;
        }
        
        auto iter_u = removed_grid_u_.begin();
        for (int ki = 0; ki < mm; ++ki)
        {
            while ((iter_u != removed_grid_u_.end()) && (*iter_u + knot_tol < knots_u[ki]))
            {
                ++iter_u;
            }
            if ((iter_u != removed_grid_u_.end()) && (*iter_u - knots_u[ki] < knot_tol))
            {
                continue;
            }

            pos_der_v.insert(pos_der_v.end(),
                             data[4*(kj*mm+ki)].begin(), data[4*(kj*mm+ki)].end());
            der_u_der_uv.insert(der_u_der_uv.end(),
                                data[4*(kj*mm+ki)+1].begin(), data[4*(kj*mm+ki)+1].end());
        }
        iter_u = removed_grid_u_.begin();
        for (int ki = 0; ki < mm; ++ki)
        {
            while ((iter_u != removed_grid_u_.end()) && (*iter_u + knot_tol < knots_u[ki]))
            {
                ++iter_u;
            }
            if ((iter_u != removed_grid_u_.end()) && (*iter_u - knots_u[ki] < knot_tol))
            {
                continue;
            }

            pos_der_v.insert(pos_der_v.end(),
                             data[4*(kj*mm+ki)+2].begin(), data[4*(kj*mm+ki)+2].end());
            der_u_der_uv.insert(der_u_der_uv.end(),
                                data[4*(kj*mm+ki)+3].begin(), data[4*(kj*mm+ki)+3].end());
        }
    }

    vector<double> coefs_pos_der_v;
    coefs_pos_der_v.reserve(2*nn_red*mm_red*dim);
    vector<double> param_v = grid_.getKnots(false);
    vector<double> param_v2(nn_red*2);
    int cntr = 0;
    auto iter = removed_grid_v_.begin();
    for (size_t ki = 0; ki < param_v.size(); ++ki)
    {
        while ((iter != removed_grid_v_.end()) && (*iter + knot_tol < param_v[ki]))
        {
            ++iter;
        }        
        if ((iter != removed_grid_v_.end()) && (*iter - param_v[ki] < knot_tol))
        {
            continue;
        }

        param_v2[2*cntr] = param_v[ki];
        param_v2[2*cntr+1] = param_v[ki];
        ++cntr;
    }

    // We first interpolate in the second parameter direction, pos & der_v.
    HermiteInterpolator interpolator;
    // The number of points includes the tangents.
    interpolator.interpolate(2*nn_red, mm_red*dim, &param_v2[0], &pos_der_v[0], coefs_pos_der_v);

#if 0
    std::cout << "Result from interpolating pos_der_v:" << std::endl;
    for (size_t ki = 0; ki < coefs_pos_der_v.size(); ki+= dim)
    {
        std::cout << coefs_pos_der_v[ki] << " " << coefs_pos_der_v[ki+1] << " " << coefs_pos_der_v[ki+2] << std::endl;
    }
#endif
    
    vector<double> coefs_der_u_der_uv;
    // Then der_u & der_uv.
    interpolator.interpolate(2*nn_red, mm_red*dim, &param_v2[0], &der_u_der_uv[0], coefs_der_u_der_uv);
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
    coefs_pre_u_int.reserve(4*mm_red*nn_red*dim);
    for (int ki = 0; ki < mm_red; ++ki)
    {
        for (int kj = 0; kj < 2*nn_red; ++kj)
        {
//            size_t ind = 2*(ki*mm + kj);
            coefs_pre_u_int.insert(coefs_pre_u_int.end(),
                                   coefs_pos_der_v.begin() + (kj*mm_red + ki)*dim,
                                   coefs_pos_der_v.begin() + (kj*mm_red + ki + 1)*dim);
        }
        for (int kj = 0; kj < 2*nn_red; ++kj)
        {
//            size_t ind = 2*(ki*mm + kj);
            coefs_pre_u_int.insert(coefs_pre_u_int.end(),
                                   coefs_der_u_der_uv.begin() + (kj*mm_red + ki)*dim,
                                   coefs_der_u_der_uv.begin() + (kj*mm_red + ki + 1)*dim);
        }
    }

#if 0
    std::cout << "coefs_pre_u_int.size(): " << coefs_pre_u_int.size() << std::endl;
#endif
    
    vector<double> sf_coefs;
    vector<double> param_u = grid_.getKnots(true);
    vector<double> param_u2(mm_red*2);
    iter = removed_grid_u_.begin();
    cntr = 0;
    for (size_t ki = 0; ki < param_u.size(); ++ki)
    {
        while ((iter != removed_grid_u_.end()) && (*iter + knot_tol < param_u[ki]))
        {
            ++iter;
        }        
        if ((iter != removed_grid_u_.end()) && (*iter - param_u[ki] < knot_tol))
        {
            continue;
        }

        param_u2[2*cntr] = param_u[ki];
        param_u2[2*cntr+1] = param_u[ki];
        ++cntr;
    }
    // We then interpolate in the first parameter direction.
    interpolator.interpolate(2*mm_red, 2*nn_red*dim, &param_u2[0], &coefs_pre_u_int[0], sf_coefs);
    BsplineBasis basis_u = interpolator.basis();

    // We then transpose the sf coefs.
    vector<double> sf_coefs_tr(sf_coefs.size());
    for (int ki = 0; ki < 2*mm_red; ++ki)
    {
        for (int kj = 0; kj < 2*nn_red; ++kj)
        {
            for (int kk = 0; kk < dim; ++kk)
            {
                sf_coefs_tr[(kj*2*mm_red+ki)*dim+kk] = sf_coefs[(ki*2*nn_red+kj)*dim+kk];
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

    MESSAGE("INFO: num_coefs_u: " << sf->numCoefs_u() << ", num_coefs_v: " << sf->numCoefs_v());

    return sf;
}

const HermiteGrid2D& HermiteApprEvalSurf::getGrid() const
{
    return grid_;
}

void HermiteApprEvalSurf::removeGridLines(const std::vector<double>& grid_lines_u,
                                          const std::vector<double>& grid_lines_v)
{
    MESSAGE("We should verify that the input knots are unique and included in the grid!");
//    grid_.removeGridLines(grid_lines_u, grid_lines_v);
    removed_grid_u_ = grid_lines_u;
    removed_grid_v_ = grid_lines_v;
}


void HermiteApprEvalSurf::setNoSplit(const std::vector<shared_ptr<SplineCurve> >& no_split_cvs_2d,
                                     double no_split_2d_tol)
{
//    std::cout << "DEBUG: no_split_cvs_2d.size(): " << no_split_cvs_2d.size() << std::endl;
    no_split_cvs_2d_ = no_split_cvs_2d;
    no_split_2d_tol_ = no_split_2d_tol;

    // We compute the bounding box for the cvs.
    for (size_t ki = 0; ki < no_split_cvs_2d_.size(); ++ki)
    {
        BoundingBox bd_box = no_split_cvs_2d_[ki]->boundingBox();
        no_split_bd_box_.push_back(bd_box);

    }
}


// A simplified curvature function averaging the legs and calculating the radius of the corresponding
// circle. We do not care about the direction of the curvature (i.e. convex vs concave).
double avgAngle(Point bezcoef[16], bool dir_is_u)
{
    double avg_ang = 0.0; //avg_curv = 0.0;
    const int step = (dir_is_u) ? 1 : 4;
    for (int ki = 0; ki < 4; ++ki)
    {
        int first = (dir_is_u) ? 4*ki : ki;
        for (int kj = 0; kj < 2; ++kj)
        {
            Point leg1 = bezcoef[first+step]- bezcoef[first];
            Point leg2 = bezcoef[first+2*step]- bezcoef[first+step];
            //Point avg_leg = 0.5*(leg1 + leg2);
            double ang = leg1.angle(leg2); // Range [0, MPI).
            avg_ang += ang;
            //double fraction = ang/(0.5*M_PI);
            // We then compute how many legs are needed to complete the piecewise circle, using the
            // inscribed circle.
            // std::cout << "ang: " << ang << std::endl;
            // double radius;
        }
    }

    avg_ang /= 8.0;
    return avg_ang;
}


int HermiteApprEvalSurf::splitDomain(double spar1, double epar1, double spar2, double epar2, Point bezcoef[16],
                                     bool& dir_is_u, double& new_knot)
{
    int status = 0;

    double dom1 = epar1 - spar1;
    double dom2 = epar2 - spar2;

    // We multiply the dom size by a curvature weight.
    const double min_curv_weight = 1.0;
    const double max_curv_weight = 10.0;
    // Replaced curvature with average angle for grid vectors.
    double avg_ang_dir_u = avgAngle(bezcoef, true);
    double avg_ang_dir_v = avgAngle(bezcoef, false);
    // We make sure that the weight is in the range [1.0, 10.0].
    double wgt1 = std::max(min_curv_weight, std::min(max_curv_weight, avg_ang_dir_u/avg_ang_dir_v));
    double wgt2 = std::max(min_curv_weight, std::min(max_curv_weight, avg_ang_dir_v/avg_ang_dir_u));

    // We split (bisect) in the direction with the largest interval (weighted, wrt curvature).
    dir_is_u = (wgt1*dom1 > wgt2*dom2);
    new_knot = (dir_is_u) ? 0.5*(spar1 + epar1) : 0.5*(spar2 + epar2);

    if (no_split_bd_box_.size() == 0)
    {
        return status;
    }
    else
    {
        // We check for intersection with boxes in no_split_bd_box_.
        Point low(spar1, spar2);
        Point high(epar1, epar2);
        BoundingBox domain_box;
        domain_box.setFromPoints(low, high);
        // If there are kink curves inside the domain we must avoid splitting close to those kinks.
        vector<double> kinks_u, kinks_v;
        //vector<pair<double> > ranges_u, ranges_v;
        for (size_t ki = 0; ki < no_split_bd_box_.size(); ++ki)
        {
        
            bool overlap = domain_box.overlaps(no_split_bd_box_[ki]);//, no_split_2d_tol_);
            if (overlap)
            {
                // If there is an overlap we must carefully select the split direction and parameter.
                //MESSAGE("DEBUG: We have an overlap between domain and kink curve! ki = " << ki);
                Point no_low = no_split_bd_box_[ki].low();
                Point no_high = no_split_bd_box_[ki].high();
                double width_u = no_high[0] - no_low[0];
                double width_v = no_high[1] - no_low[1];
                //std::cout << "width_u: " << width_u << ", width_v: " << width_v << std::endl;
                // The direction with the smallest width is assumed to correspond to a kink along an iso curve.
                if (width_u < width_v)
                {
                    kinks_u.push_back(0.5*(no_low[0] + no_high[0]));
                }
                else
                {
                    kinks_v.push_back(0.5*(no_low[1] + no_high[1]));
                }
            }
        }
    
        if (kinks_u.size() > 0)
        {
            MESSAGE("DEBUG: Handling kinks in the v direction (constant u parameter)!");
            std::sort(kinks_u.begin(), kinks_u.end());
            vector<double> ranges(kinks_u.begin(), kinks_u.end());
            if (low[0] < ranges[0])
            {
                ranges.insert(ranges.begin(), low[0]);
            }
            if (high[0] > ranges.back())
            {
                ranges.push_back(high[0]);
            }
            int max_range_ind = -1;
            double max_range = -1.0;
            for (size_t ki = 0; ki < ranges.size() - 1; ++ki)
            {
                double range = ranges[ki+1] - ranges[ki];
                if (range > max_range)
                {
                    max_range = range;
                    max_range_ind = ki;
                }
            }
#if 0
#ifndef NDEBUG
            std::cout << "DEBUG: max_range: [" << ranges[max_range_ind] << ", " << ranges[max_range_ind+1] << "]" <<
                ", spar2: " << spar2 << ", epar2: " << epar2 << std::endl;
#endif
#endif

            double new_dom1 = ranges[max_range_ind+1] - ranges[max_range_ind];
            if (new_dom1 < no_split_2d_tol_)
            {
                MESSAGE("We should not split in the u direction!");
                dir_is_u = false;
                if (dom2 < no_split_2d_tol_)
                {
                    MESSAGE("We should not split in the v direction either!");
                    status = 3;
                }
                else
                {
                    status = 1;
                    new_knot = 0.5*(spar2 + epar2);
                }
            }
            else
            {
                // We split (bisect) in the direction with the largest interval (weighted, wrt curvature).
                if (wgt1*new_dom1 > wgt2*dom2)
                {
                    double new_knot_u = 0.5*(ranges[max_range_ind+1] + ranges[max_range_ind]);
#ifndef NDEBUG
                    MESSAGE("INFO: Replacing the u knot " << new_knot << " with " << new_knot_u <<
                            ", spar2: " << spar2 << ", epar2: " << epar2);
#endif
                    dir_is_u = true;
                    new_knot = new_knot_u;
                }
                else
                {
                    if (dom2 < no_split_2d_tol_)
                    {
                        MESSAGE("We should not split in the v direction!");
                        // status = 2;
                    }
                    dir_is_u = false;
                    new_knot = 0.5*(spar2 + epar2);
                }
            }
        }
        if (kinks_v.size() > 0)
        {
            MESSAGE("DEBUG: Missing! We must handle kinks in the u direction (constant v parameter)!");

        }
    
        return status;
    }
}
    
}
