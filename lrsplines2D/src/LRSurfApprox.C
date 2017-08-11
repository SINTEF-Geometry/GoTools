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

#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include "GoTools/lrsplines2D/LRSurfSmoothLS.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"
#include "GoTools/lrsplines2D/LinDepUtils.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRSplineMBA.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

//#define DEBUG

using std::vector;
using std::cout;
using std::endl;
using namespace Go;

//==============================================================================
LRSurfApprox::LRSurfApprox(vector<double>& points, 
			   int dim, double epsge,  bool init_mba, 
			   double mba_level,
			   bool closest_dist, bool repar)
  : nmb_pts_((int)points.size()/(2+dim)), points_(points), useMBA_(false), 
    toMBA_(4), initMBA_(init_mba), initMBA_coef_(mba_level), 
    maxdist_(-10000.0), maxdist_prev_(-10000.0), avdist_(0.0), 
    avdist_all_(0), avdist_all_prev_(0), outsideeps_(0), aepsge_(epsge), 
    smoothweight_(1.0e-3), 
    smoothbd_(false), repar_(repar), check_close_(closest_dist), 
    fix_corner_(false), to3D_(-1), grid_(false), check_init_accuracy_(false),
    initial_surface_(false), has_min_constraint_(false), has_max_constraint_(false),
  has_local_constraint_(false), verbose_(false)
//==============================================================================
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  grid_start_[0] = grid_start_[1] = 0.0;
  cell_size_[0] = cell_size_[1] = 1.0;
  usize_min_ = vsize_min_ = -1;

  fix_boundary_ = false; //true;
  make_ghost_points_ = false;

  // if (dim > 1)
  //   {
  //     initMBA_ = false;
  //     useMBA_ = false;
  //     toMBA_ = 10e4;  // A large number
  //   }

  // Create an LR B-spline surface with the domain given by the 
  // parameter domain of the points. Only one element will be
  // created
  makeInitSurf(dim);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(shared_ptr<SplineSurface>& srf,
			   vector<double>& points, 
			   double epsge, bool closest_dist,
			   bool repar)
  : points_(points), useMBA_(false), toMBA_(4), initMBA_(false), 
    initMBA_coef_(0.0), 
    maxdist_(-10000.0), maxdist_prev_(-10000.0), avdist_(0.0), 
    avdist_all_(0.0), avdist_all_prev_(0), 
    outsideeps_(0), aepsge_(epsge), smoothweight_(1.0e-3), smoothbd_(false), 
    repar_(repar), check_close_(closest_dist), 
    fix_corner_(false), to3D_(-1), grid_(false), check_init_accuracy_(false),
    initial_surface_(true), has_min_constraint_(false), has_max_constraint_(false), 
    has_local_constraint_(false), verbose_(false)
//==============================================================================
{
  nmb_pts_ = (int)points.size()/(2+srf->dimension());
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  grid_start_[0] = grid_start_[1] = 0.0;
  cell_size_[0] = cell_size_[1] = 1.0;
  fix_boundary_ = false; //true;
  make_ghost_points_ = false;
  usize_min_ = vsize_min_ = -1;

  // if (srf->dimension() > 1)
  //   {
  //     initMBA_ = false;
  //     useMBA_ = false;
  //     toMBA_ = 10e4;  // A large number
  //   }

  // Create an LR B-spline surface based on the given spline surface
  makeInitSurf(srf);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(shared_ptr<LRSplineSurface>& srf,
			   vector<double>& points, 
			   double epsge, bool closest_dist,
			   bool repar, bool check_init_accuracy)
//==============================================================================
  : points_(points), useMBA_(false), toMBA_(4), initMBA_(false), 
    initMBA_coef_(0.0),
    maxdist_(-10000.0), maxdist_prev_(-10000.0), avdist_(0.0), 
    avdist_all_(0.0), avdist_all_prev_(0), 
    outsideeps_(0), aepsge_(epsge), smoothweight_(1.0e-3), 
    smoothbd_(false), repar_(repar), check_close_(closest_dist), 
    fix_corner_(false), to3D_(-1), check_init_accuracy_(check_init_accuracy), 
    grid_(false), initial_surface_(true), has_min_constraint_(false), 
    has_max_constraint_(false), has_local_constraint_(false), verbose_(false)
{
  nmb_pts_ = (int)points.size()/(2+srf->dimension());
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  grid_start_[0] = grid_start_[1] = 0.0;
  cell_size_[0] = cell_size_[1] = 1.0;
  fix_boundary_ = false; //true;
  make_ghost_points_ = false;
  srf_ = srf;
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
  usize_min_ = vsize_min_ = -1;

  // if (srf->dimension() > 1)
  //   {
  //     initMBA_ = false;
  //     useMBA_ = false;
  //     toMBA_ = 10e4;  // A large number
  //   }
}

//==============================================================================
LRSurfApprox::LRSurfApprox(int ncoef_u, int order_u, int ncoef_v, int order_v,
			   vector<double>& points, 
			   int dim, double epsge, bool init_mba, 
			   double mba_level,
			   bool closest_dist, bool repar)
//==============================================================================
  : nmb_pts_((int)points.size()/(2+dim)), points_(points), useMBA_(false),
    toMBA_(4), initMBA_(init_mba), initMBA_coef_(mba_level), 
    maxdist_(-10000.0), maxdist_prev_(-10000.0), avdist_(0.0), 
    avdist_all_(0.0), avdist_all_prev_(0), outsideeps_(0), aepsge_(epsge), 
    smoothweight_(1.0e-3), 
    smoothbd_(false), repar_(repar), check_close_(closest_dist), 
    fix_corner_(false), to3D_(-1), grid_(false), check_init_accuracy_(false),
    initial_surface_(false), has_min_constraint_(false), has_max_constraint_(false),
    has_local_constraint_(false), verbose_(false)
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  grid_start_[0] = grid_start_[1] = 0.0;
  cell_size_[0] = cell_size_[1] = 1.0;
  usize_min_ = vsize_min_ = -1;

  fix_boundary_ = false; //true;
  make_ghost_points_ = false;

  // if (dim > 1)
  //   {
  //     initMBA_ = false;
  //     useMBA_ = false;
  //     toMBA_ = 10e4;  // A large number
  //   }

  // Create an LR B-spline surface with unset coefficients and the domain
  // given by the parameter domain of the points. The size of the spline
  // space is given. The knots will be equally spaced
  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(int order_u, vector<double>& knots_u, 
			   int order_v, vector<double>& knots_v,
			   vector<double>& points, int dim, 
			   double epsge, bool init_mba, 
			   double mba_level,
			   bool closest_dist, bool repar)
//==============================================================================
  : nmb_pts_((int)points.size()/(2+dim)), points_(points), useMBA_(false),
    toMBA_(4), initMBA_(init_mba), initMBA_coef_(mba_level), 
    maxdist_(-10000.0), maxdist_prev_(-10000.0), avdist_(0.0), 
    avdist_all_(0.0), avdist_all_prev_(0), outsideeps_(0), aepsge_(epsge), 
  smoothweight_(1.0e-3), 
    smoothbd_(false), repar_(repar), check_close_(closest_dist), 
    fix_corner_(false), to3D_(-1), grid_(false), check_init_accuracy_(false),
    initial_surface_(false), has_min_constraint_(false), has_max_constraint_(false),
    has_local_constraint_(false), verbose_(false)
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  grid_start_[0] = grid_start_[1] = 0.0;
  cell_size_[0] = cell_size_[1] = 1.0;
  usize_min_ = vsize_min_ = -1;

  fix_boundary_ = false; //true;
  make_ghost_points_ = false;

  // if (dim > 1)
  //   {
  //     initMBA_ = false;
  //     useMBA_ = false;
  //     toMBA_ = 10e4;  // A large number
  //   }

  // Compute domain
  double domain[4]; //umin, umax, vmin, vmax;
  computeParDomain(dim, domain[0], domain[1], domain[2], domain[3]);

  // Check if the knot vectors should be extended
  int ki;
  for (ki=order_u; ki<(int)knots_u.size(); ++ki)
    {
      double knot = knots_u[ki];
      if (knot > domain[0] && 0.5*(knots_u[ki-1]+knot) < domain[0])
	knots_u.insert(knots_u.begin()+ki, 0.5*(knots_u[ki-1]+knot));
      if (knot > domain[0])
	break;
    }
  for (ki=(int)knots_u.size()-order_u-1; ki>0; --ki)
    {
      double knot = knots_u[ki];
      if (knot < domain[1] && 0.5*(knots_u[ki+1]+knot) > domain[1])
	knots_u.insert(knots_u.begin()+ki+1, 0.5*(knots_u[ki+1]+knot));
      if (knot < domain[1])
	break;
    }

  for (ki=order_v; ki<(int)knots_v.size(); ++ki)
    {
      double knot = knots_v[ki];
      if (knot > domain[2] && 0.5*(knots_v[ki-1]+knot) < domain[2])
	knots_v.insert(knots_v.begin()+ki, 0.5*(knots_v[ki-1]+knot));
      if (knot > domain[2])
	break;
    }
  for (ki=(int)knots_v.size()-order_v-1; ki>0; --ki)
    {
      double knot = knots_v[ki];
      if (knot < domain[3] && 0.5*(knots_v[ki+1]+knot) > domain[3])
	knots_v.insert(knots_v.begin()+ki+1, 0.5*(knots_v[ki+1]+knot));
      if (knot < domain[3])
	break;
    }
 
  // Create an LR B-spline representation of a tensor-product spline surface 
  // with unset coefficients and the given knots
  int ncoef_u = (int)knots_u.size() - order_u;
  int ncoef_v = (int)knots_v.size() - order_v;
  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v, &knots_u[0], &knots_v[0]);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(int ncoef_u, int order_u, int ncoef_v, int order_v,
			   vector<double>& points, int dim, 
			   double domain[4], double epsge, bool init_mba, 
			   double mba_level,
			   bool closest_dist, bool repar)
//==============================================================================
  : nmb_pts_((int)points.size()/(2+dim)), points_(points), useMBA_(false),
    toMBA_(4), initMBA_(init_mba), initMBA_coef_(mba_level), 
    maxdist_(-10000.0), maxdist_prev_(-10000.0), avdist_(0.0), 
    avdist_all_(0.0), avdist_all_prev_(0), outsideeps_(0), aepsge_(epsge), 
    smoothweight_(1.0e-3), 
    smoothbd_(false), repar_(repar), check_close_(closest_dist), 
    fix_corner_(false), to3D_(-1), grid_(false), check_init_accuracy_(false),
    initial_surface_(false), has_min_constraint_(false), has_max_constraint_(false),
    has_local_constraint_(false), verbose_(false)
{
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  grid_start_[0] = grid_start_[1] = 0.0;
  cell_size_[0] = cell_size_[1] = 1.0;
  usize_min_ = vsize_min_ = -1;

  fix_boundary_ = false; //true;
  make_ghost_points_ = false;

  // if (dim > 1)
  //   {
  //     initMBA_ = false;
  //     useMBA_ = false;
  //     toMBA_ = 10e4;  // A large number
  //   }

  // Create an LR B-spline surface with unset coefficients and the domain
  // given by the parameter domain of the points. The size of the spline
  // space is given. The knots will be equally spaced
  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v, domain);
}

//==============================================================================
LRSurfApprox::~LRSurfApprox()
//==============================================================================
{
}

//==============================================================================
 shared_ptr<LRSplineSurface> LRSurfApprox::getApproxSurf(double& maxdist, 
							 double& avdist_all,
							 double& avdist,
							 int& nmb_out_eps, 
							 int max_iter)
//==============================================================================
{
//   // We start the timer.
// #ifdef _OPENMP
//   double time0 = omp_get_wtime();
// #endif

#ifdef _OPENMP
    // When using OpenMP we choose between splitting the threads on
    // the surface elements or on the points for each element.  As the
    // iteration progresses we should turn towards splitting on the
    // elements.  Initial switch threshold set to num_elem ==
    // avg_num_pnts_per_elem.
    const int num_elem = srf_->numElements();
    const int num_pts = points_.size()/(srf_->dimension());
    // We let the number of elem vs average numer of points per elem be the threshold
    // for switching the OpenMP level.
    const double pts_per_elem = num_pts/num_elem;
    const bool omp_for_elements = (num_elem > pts_per_elem); // As opposed to element points.
    const bool omp_for_mba_update = true;
#ifndef NDEBUG
    std::cout << "num_elem: " << num_elem << ", pts_per_elem: " << pts_per_elem << ", openmp_for_elements: " <<
	omp_for_elements << std::endl;
#endif
#else
    const bool omp_for_elements = true;//false; // 201503 The omp version seems to be faster even when run sequentially.
    const bool omp_for_mba_update = true;//false; // 201503 The omp version seems to be faster even when run sequentially.
#endif

#ifdef DEBUG
  std::ofstream of0("init0_sf.g2");
  shared_ptr<LRSplineSurface> tmp0(srf_->clone());
  if (tmp0->dimension() == 1)
    tmp0->to3D();
  tmp0->writeStandardHeader(of0);
  tmp0->write(of0);
  of0 << std::endl;
  // LineCloud lines0 = tmp0->getElementBds();
  // lines0.writeStandardHeader(of0);
  // lines0.write(of0);
  std::ofstream of02("init0_tpsf.g2");
  shared_ptr<SplineSurface> ssf0(tmp0->asSplineSurface());
  ssf0->writeStandardHeader(of02);
  ssf0->write(of02);
#endif

  if (srf_->dimension() == 3)
    {
      // Reparameterize to reflect the surface size
      double len1, len2;
      srf_->estimateSfSize(len1, len2);

      double umin = srf_->paramMin(XFIXED);
      double umax = srf_->paramMax(XFIXED);
      double vmin = srf_->paramMin(YFIXED);
      double vmax = srf_->paramMax(YFIXED);
      srf_->setParameterDomain(umin, umin+len1, vmin, vmin+len2);

      // Reparameterize also data points
      int del = 5;  // Parameter pair + geometric dimension
      int nmb = (int)points_.size()/del;
      for (int kj=0; kj<nmb; ++kj)
	{
	  points_[kj*del] = umin + (points_[kj*del] - umin)*len1/(umax - umin);
	  points_[kj*del+1] = vmin + (points_[kj*del+1] - vmin)*len2/(vmax - vmin);
	}
    }

  LRSurfSmoothLS LSapprox;

  if (make_ghost_points_ && !initial_surface_ && srf_->dimension() == 1 && 
      !useMBA_)
    {
      // This is experimental code and should, if kept, be integrated
      // with LRSurfSmoothLS::addDataPoints
      vector<double> ghost_points;
      constructGhostPoints(ghost_points);
      LRSplineUtils::distributeDataPoints(srf_.get(), ghost_points, true, false);
    }

  // Initiate with data points
  LRSplineUtils::distributeDataPoints(srf_.get(), points_, true, true);

  if (make_ghost_points_ && initial_surface_)
    {
      // No need to construct ghost points from extrapolation. Use
      // the input surface
      constructInnerGhostPoints();
    }
      
  vector<Element2D*> ghost_elems;
  if (check_init_accuracy_ /*|| useMBA_*/)
  {
    // Compute accuracy in data points
      if (omp_for_elements)
	  computeAccuracy_omp(ghost_elems);
      else
	  computeAccuracy(ghost_elems);
  }

  // Initiate approximation engine
  if (fix_corner_)
    setCoefKnown();

  if (fix_boundary_)
    setFixBoundary(true);

  // Initial approximation of LR B-spline surface
  if (/*useMBA_ || */initMBA_)
  {
      if (omp_for_mba_update && srf_->dimension() == 1)
      {
	  LRSplineMBA::MBADistAndUpdate_omp(srf_.get());
      }
      else
      {
	  LRSplineMBA::MBADistAndUpdate(srf_.get());
      }
      //LRSplineMBA::MBAUpdate(srf_.get());
      if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
	adaptSurfaceToConstraints();
      // computeAccuracy();
      // LRSplineMBA::MBAUpdate(srf_.get());
      if (omp_for_mba_update && srf_->dimension() == 1)
      {
	  LRSplineMBA::MBADistAndUpdate_omp(srf_.get());
      }
      else
      {
	  LRSplineMBA::MBADistAndUpdate(srf_.get());
      }
     if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
     	adaptSurfaceToConstraints();
     LSapprox.setInitSf(srf_, coef_known_);
     updateCoefKnown();
    }
  else
    {
     LSapprox.setInitSf(srf_, coef_known_);
      LSapprox.updateLocals();
      //performSmooth(&LSapprox);
     if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
     	adaptSurfaceToConstraints();
    }

#ifdef DEBUG
  std::ofstream of1("init_sf.g2");
  shared_ptr<LRSplineSurface> tmp;
  if (srf_->dimension() == 1)
    {
      tmp = shared_ptr<LRSplineSurface>(srf_->clone());
      tmp->to3D();
    }
  else
    tmp = srf_;
  tmp->writeStandardHeader(of1);
  tmp->write(of1);
  std::ofstream of12("init_tpsf.g2");
  shared_ptr<SplineSurface> ssf1(tmp->asSplineSurface());
  ssf1->writeStandardHeader(of12);
  ssf1->write(of12);
  of12 << std::endl;
  LineCloud lines = tmp->getElementBds();
  lines.writeStandardHeader(of12);
  lines.write(of12);
#endif

  // Compute accuracy in data points
  if (omp_for_elements)
      computeAccuracy_omp(ghost_elems);
  else
      computeAccuracy(ghost_elems);

  if (verbose_)
    {
      std::cout << "Number of data points: " << nmb_pts_ << std::endl;
      std::cout << "Number of coefficients: " << srf_->numBasisFunctions() << std::endl;
      std::cout << "Initial surface. Maximum distance: " << maxdist_;
      std::cout << ", average distance: " << avdist_all_ << std::endl;
      std::cout << "Number of points outside tolerance: " << outsideeps_;
      std::cout << ", average distance in outside points: " << avdist_ << std::endl;
    }

  ghost_elems.clear();
  points_.clear();  // Not used anymore TESTING
  for (int ki=0; ki<max_iter; ++ki)
    {
      // Check if the requested accuracy is reached
      if (maxdist_ <= aepsge_ || outsideeps_ == 0)
	break;

      // Refine surface
      prev_ =  shared_ptr<LRSplineSurface>(srf_->clone());

      // Check if any ghost points need to be updated
      if (!useMBA_ && ki<toMBA_ && ghost_elems.size() > 0)
	{
	  updateGhostElems(ghost_elems);
	}

      if (ki > 0 || (!initial_surface_))
	{
	  int nmb_refs = refineSurf();
	  if (nmb_refs == 0)
	    break;  // No refinements performed
	}
      //refineSurf2();
#ifdef DEBUG
      std::ofstream of2("refined_sf.g2");
      shared_ptr<LRSplineSurface> tmp2;
      if (srf_->dimension() == 1)
	{
	  tmp2 = shared_ptr<LRSplineSurface>(srf_->clone());
	  tmp2->to3D();
	}
      else
	tmp2 = srf_;
      tmp2->writeStandardHeader(of2);
      tmp2->write(of2);
      of2 << std::endl;
      LineCloud lines2 = tmp2->getElementBds();
      lines2.writeStandardHeader(of2);
      lines2.write(of2);

      std::ofstream of3("point_clouds.g2");
      int del = srf_->dimension() + 3; // Parameter pair, position and 
      // distance between surface and point
      for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
	   it != srf_->elementsEnd(); ++it)
	{
	  vector<double>& elem_data = it->second->getDataPoints();
	  int nmb = (int)elem_data.size()/del;
	  if (elem_data.size() > 0)
	    {
	      vector<double> tmppt;
	      if (srf_->dimension() == 1)
		tmppt = elem_data;
	      else
		{
		  tmppt.reserve(3*elem_data.size()/del);
		  for (int kr=0; kr<nmb; ++kr)
		    tmppt.insert(tmppt.end(), elem_data.begin()+kr*del+2, 
				 elem_data.begin()+(kr+1)*del-1);
		  
		}
	      PointCloud3D cloud(tmppt.begin(), nmb);
	      cloud.writeStandardHeader(of3);
	      cloud.write(of3);
	    }
	}
#endif

      // Update coef_known from information in LR B-splines
      updateCoefKnown();
      //unsetCoefKnown(); TEST
      if (fix_corner_)
	setCoefKnown();
      if (fix_boundary_)
	setFixBoundary(true);
  
      // Check for linear independence (overloading)
      vector<LRBSpline2D*> funs = LinDepUtils::unpeelableBasisFunctions(*srf_);
#ifdef DEBUG
      std::cout << "Number of unpeelable functions: " << funs.size() << std::endl;
#endif
     
      // Construct additional ghost points in elements with a low
      // distribution of points
      bool ghost_points_inner = true;
      //if (false /*make_ghost_points_ && ki>0*/)
      if (make_ghost_points_ && ki>0 && ghost_points_inner &&
	  !useMBA_ && ki<toMBA_)
	{
	  constructInnerGhostPoints();
	}

      // if ((has_min_constraint_ || has_max_constraint_) && 
      // 	  srf_->dimension() == 1)
      // 	{
      // 	  addConstraintGhostPoints();
      // 	}

       // Update surface
      if (useMBA_ || ki >= toMBA_)
      {
	if (srf_->dimension() == 3)
	  {
	    LRSplineMBA::MBADistAndUpdate(srf_.get());
	  }
	else if (omp_for_mba_update)
	  {
	      LRSplineMBA::MBAUpdate_omp(srf_.get());
	  }
	  else
	  {
	      LRSplineMBA::MBAUpdate(srf_.get());
	  }
	  if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
	    adaptSurfaceToConstraints();
	  if (omp_for_mba_update && srf_->dimension() == 1)
	  {
	      LRSplineMBA::MBADistAndUpdate_omp(srf_.get());
	  }
	  else
	  {
	      LRSplineMBA::MBADistAndUpdate(srf_.get());
	  }
	  // computeAccuracy();
	  // LRSplineMBA::MBAUpdate(srf_.get());
	  if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
	    adaptSurfaceToConstraints();
	}
      else
	{
	  try {
	    LSapprox.updateLocals();
	    performSmooth(&LSapprox);
	    if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
	      adaptSurfaceToConstraints();
	  }
	  catch (...)
	    {
	      // Surface update failed.
	      if (srf_->dimension() == 3)
		{
		  // Surface update failed. Return previous surface
		  srf_ = prev_;
		  break;
		}
	      else
		{
		  // Switch to MBA method
		  useMBA_ = true;
		  if (srf_->dimension() == 3)
		    {
		      LRSplineMBA::MBADistAndUpdate(srf_.get());
		    }
		  else if (omp_for_mba_update)
		    {
		      LRSplineMBA::MBAUpdate_omp(srf_.get());
		    }
		  else
		    {
		      LRSplineMBA::MBAUpdate(srf_.get());
		    }
		  if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
		    adaptSurfaceToConstraints();
		  if (0)//omp_for_mba_update)
		    {
		      LRSplineMBA::MBADistAndUpdate_omp(srf_.get());
		    }
		  else
		    {
		      LRSplineMBA::MBADistAndUpdate(srf_.get());
		    }
		  // computeAccuracy();
		  // LRSplineMBA::MBAUpdate(srf_.get());
		  if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
		    adaptSurfaceToConstraints();
		  //break;
		}
	    }
	}
  
#ifdef DEBUG
      std::ofstream of4("updated_sf.g2");
#endif
      shared_ptr<LRSplineSurface> tmp3;
      if (srf_->dimension() == 1)
	{
	  tmp3 = shared_ptr<LRSplineSurface>(srf_->clone());
	  tmp3->to3D();
	}
      else
	tmp3 = srf_;
#ifdef DEBUG
      tmp3->writeStandardHeader(of4);
      tmp3->write(of4);
      LineCloud lines3 = tmp3->getElementBds();
      lines3.writeStandardHeader(of4);
      lines3.write(of4);
      std::ofstream of42("updated_tpsf.g2");
      shared_ptr<SplineSurface> ssf4(tmp3->asSplineSurface());
      ssf4->writeStandardHeader(of42);
      ssf4->write(of42);
      of42 << std::endl;
      LineCloud lines32 = tmp3->getElementBds();
      lines32.writeStandardHeader(of42);
      lines32.write(of42);
#endif
  
      if (ki == to3D_ && srf_->dimension() == 1)
	{
	  // Turn the current function into a 3D surface
	  // before continuing the iteration
	  turnTo3D();
	}

      maxdist_prev_ = maxdist_;
      avdist_all_prev_ = avdist_all_;

      ghost_elems.clear();
  if (omp_for_elements)
      computeAccuracy_omp(ghost_elems);
  else
      computeAccuracy(ghost_elems);
      if (srf_->dimension() == 1 && (maxdist_ > 1.1*maxdist_prev_ ||
				     avdist_all_ > 1.1*avdist_all_prev_))
      	useMBA_ = true;

      if (verbose_)
	{
	  std::cout << std::endl;
	  std::cout << "Iteration number " << ki+1 <<". Maximum distance: " << maxdist_;
	  std::cout << ", average distance: " << avdist_all_ << std::endl;
	  std::cout << "Number of points outside tolerance: " << outsideeps_;
	  std::cout << ", average distance in outside points: " << avdist_ << std::endl;
	  std::cout << "Number of coefficients: " << srf_->numBasisFunctions() << std::endl;
	}
    }

  // Set accuracy information
  maxdist = maxdist_;
  avdist_all = avdist_all_;
  avdist = avdist_;
  nmb_out_eps = outsideeps_;

// #ifdef _OPENMP
//   double time1 = omp_get_wtime();
//   double time_spent = time1 - time0;
//   std::cout << "time_spent in getApproxSurf: " << time_spent << std::endl;
// #endif
  
  return srf_;
}

//==============================================================================
void LRSurfApprox::performSmooth(LRSurfSmoothLS *LSapprox)
//==============================================================================
{
// #ifdef _OPENMP
//   double time0 = omp_get_wtime();
// #endif

  //std::cout << "Smoothing weight: " << smoothweight_ << std::endl;
  double wgt1 = 0.0;//0.8*smoothweight_;
  double wgt3 = 0.8*smoothweight_;//0.0; //0.1*smoothweight_; //0.9*smoothweight_; // 0.5*smoothweight_;
  double wgt2 = (1.0 - wgt3 -wgt1)*smoothweight_;
  double fac = 100.0;

  if (smoothweight_ > 0.0)
    LSapprox->setOptimize(wgt1, wgt2, wgt3);
  
  if (smoothbd_)
    LSapprox->smoothBoundary(fac*wgt1, fac*wgt2, fac*wgt3);

  double approx_weight = 1.0-wgt1-wgt2-wgt3;  
  const bool use_omp = true;
  if (use_omp)
  {
      LSapprox->setLeastSquares_omp(approx_weight);
  }
  else
  {
      LSapprox->setLeastSquares(approx_weight);
  }

  shared_ptr<LRSplineSurface> lrsf_out;
  int isOK = LSapprox->equationSolve(lrsf_out);
#ifdef DEBUG
  std::cout << "isOK: " << isOK << std::endl;
#endif

// #ifdef _OPENMP
//   double time1 = omp_get_wtime();
//   double time_spent = time1 - time0;
//   std::cout << "time_spent in performSmooth(): " << time_spent << std::endl;
// #endif

  srf_ = lrsf_out;
}

//==============================================================================
void LRSurfApprox::computeAccuracy(vector<Element2D*>& ghost_elems)
//==============================================================================
{
//   // We start the timer.
// #ifdef _OPENMP
//   double time0 = omp_get_wtime();
//   double time_computeAccuracyElement = 0.0;
// #endif

  // Check the accuracy of all data points, element by element
  // Note that only points more distant from the surface than the tolerance
  // are considered in avdist_ 

  // Initiate accuracy information
  maxdist_ = 0.0;
  avdist_ = 0.0;
  avdist_all_ = 0.0;
  outsideeps_ = 0;

#ifdef _OPENMP
  const bool omp_for_element_pts = true;
#else
  const bool omp_for_element_pts = false;
#endif

#ifdef DEBUG
  std::ofstream of1("error_pnts1.g2");
  std::ofstream of2("error_pnts2.g2");
  std::ofstream of3("ok_pnts1.g2");
  std::ofstream of4("ok_pnts2.g2");
  vector<Point> err1, err2, ok1, ok2;

  std::ofstream of5("accuracy_info.txt");
  vector<Element2D*> elem;
#endif

  RectDomain rd = srf_->containingDomain();
  int dim = srf_->dimension();
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  LRSplineSurface::ElementMap::const_iterator it;
  int num = srf_->numElements();
  int kj;

  double ghost_fac = 0.8;
  ghost_elems.clear();

  //for (it=srf_->elementsBegin(), kj=0; it != srf_->elementsEnd(); ++it, ++kj)
  for (it=srf_->elementsBegin(), kj=0; kj<num; ++it, ++kj)
    {
      if (!it->second->hasDataPoints())
	{
	  // Reset accuracy information in element
	  it->second->resetAccuracyInfo();
	  continue;   // No points in which to check accuracy
	}

      double umin = it->second->umin();
      double umax = it->second->umax();
      double vmin = it->second->vmin();
      double vmax = it->second->vmax();
      vector<double>& points = it->second->getDataPoints();
      vector<double>& ghost_points = it->second->getGhostPoints();
      int nmb_pts = it->second->nmbDataPoints();
      int nmb_ghost = it->second->nmbGhostPoints();

       // Local error information
      double max_err = 0.0;
      double av_err = 0.0;
      double acc_err = 0.0;
      int outside = 0;
      double acc_err_sgn = 0.0;
      double av_err_sgn = 0.0;

      // Check if the accuracy can have been changed
      const vector<LRBSpline2D*>& bsplines = it->second->getSupport();
      size_t nb;
      for (nb=0; nb<bsplines.size(); ++nb)
	if (!bsplines[nb]->coefFixed())
	  break;

      if (/*useMBA_ ||*/ nb < bsplines.size())
	{
	  // Compute distances in data points and update parameter pairs
	  // if requested
// #ifdef _OPENMP
// 	    double time0_part = omp_get_wtime();
// #endif
	  if (nmb_pts > 0)
	  {
	      if (omp_for_element_pts)
		  computeAccuracyElement_omp(points, nmb_pts, del, rd, it->second.get());
	      else
		  computeAccuracyElement(points, nmb_pts, del, rd, it->second.get());
	  }
	  
	  // Compute distances in ghost points
	  if (nmb_ghost > 0 && !useMBA_)
	  {
	      if (omp_for_element_pts)
		  computeAccuracyElement_omp(ghost_points, nmb_ghost, del, rd, it->second.get());
	      else
		  computeAccuracyElement(ghost_points, nmb_ghost, del, rd, it->second.get());
	  }
// #ifdef _OPENMP
// 	    double time1_part = omp_get_wtime();
// 	    time_computeAccuracyElement += time1_part - time0_part;
// #endif
	}

      // Accumulate error information related to data points
      int ki;
      double *curr;
#ifdef DEBUG
      int n_above = 0, n_below = 0;
#endif
      for (ki=0, curr=&points[0]; ki<nmb_pts;)
	{
	  Point curr_pt(curr+(dim==3)*2, curr+del-1);

	  // Accumulate approximation error
	  double dist2 = fabs(curr[del-1]);
	  maxdist_ = std::max(maxdist_, dist2);
	  max_err = std::max(max_err, dist2);
	  acc_err += dist2;
	  acc_err_sgn += curr[del-1];
	  avdist_all_ += dist2;
	  if (dist2 > aepsge_)
	    {
	      av_err_sgn += curr[del-1];
	      avdist_ += dist2;
	      outsideeps_++;
	      av_err += dist2;
	      outside++;
		  
#ifdef DEBUG
	      // Accumulate error points
	      if (curr[del-1] > 0)
		{
		  err1.push_back(curr_pt);
		  n_above++;
		}
	      else
		{
		  err2.push_back(curr_pt);
		  n_below++;
		}
#endif		    
	    }
	  else
	    {
#ifdef DEBUG
	      if (curr[del-1] > 0)
		ok1.push_back(curr_pt);
	      else
		ok2.push_back(curr_pt);	
#endif
	    }	     	  

	  if (dim == 3 && repar_)
	    {
	      // Check if the point has moved
	      if (curr[0] < umin || curr[0] > umax || curr[1] < vmin || curr[1] > vmax)
		{
		  // Find element
		  Element2D *elem = srf_->coveringElement(curr[0], curr[1]);
		  elem->addDataPoints(points.begin()+ki*del, 
				      points.begin()+(ki+1)*del, false);
		  it->second->eraseDataPoints(points.begin()+ki*del, 
					      points.begin()+(ki+1)*del);
		  nmb_pts--;
		}
	      else
		{
		  curr += del;
		  ki++;
		}
	    }
	  else
	    {
	      curr += del;
	      ki++;
	    }
	}
      if (outside > 0)
	{
	  av_err /= (double)outside;
	  av_err_sgn /= (double)outside;
	}

      // Previous accuracy information
      double av_prev, max_prev;
      int nmb_out_prev;
      //double acc_prev = it->second->getAccumulatedError();
      it->second->getAccuracyInfo(av_prev, max_prev, nmb_out_prev);

// #ifdef DEBUG
//       of5 << "El nmb: " << elem.size() << ", div: " << (nmb_out_prev < 0);
//       of5 << ", prev max: " << max_prev;
//       of5 << ", prev average: " << acc_prev/nmb_pts << std::endl;
//       of5 << "nmb pts: " << nmb_pts << ", curr max: " << max_err;
//       of5 << ", curr average: " << acc_err/nmb_pts << ", av sgn: ";
//       of5 << acc_err_sgn/nmb_pts << std::endl << "av outside: " << av_err;
//       of5 << ", av out sgn: " << av_err_sgn;
//       of5 << ", nmb above: " << n_above;
//       of5 << ", nmb below: " << n_below << std::endl << std::endl;
//       elem.push_back(it->second.get());
// #endif

      if (max_err > aepsge_ && max_prev > 0.0 && max_err > ghost_fac*max_prev &&
	  nmb_ghost > 0.25*nmb_pts)
	{
	  // Collect element for update of ghost points
	  ghost_elems.push_back(it->second.get());
	}

      // Store updated accuracy information in the element
      it->second->setAccuracyInfo(acc_err, av_err, max_err, outside);
#ifdef DEBUG
      int write = 0;
      if (write)
	{
	  std::cout << "Element nr " << kj << ", max err: " << max_err << ", mean err: ";
	  std::cout << av_err << ", average: " << acc_err/(double)nmb_pts << ", outside: ";
	  std::cout << outside << std::endl;
	  std::ofstream of("elem_pts.g2");
	  if (nmb_pts > 0)
	    {
	      of << "400 1 0 4 50 50 155 255" << std::endl;
	      of << nmb_pts << std::endl;
	      for (int kh1=0; kh1<nmb_pts; ++kh1)
		{
		  Point tmppt(points.begin()+kh1*del+(dim==3)*2, points.begin()+(kh1+1)*del-1);
		  of << tmppt << std::endl;
		}
	    }
	  if (nmb_ghost > 0)
	    {
	      of << "400 1 0 4 50 155 50 255" << std::endl;
	      of << nmb_ghost << std::endl;
	      for (int kh1=0; kh1<nmb_ghost; ++kh1)
		{
		  Point tmppt(ghost_points.begin()+kh1*del+(dim==3)*2, 
			      ghost_points.begin()+(kh1+1)*del-1);
		  of << tmppt << std::endl;
		}
	    }
	}
#endif

    }

  avdist_all_ /= (double)nmb_pts_;
#ifdef DEBUG
      if (err1.size() > 0)
	{
	  of1 << "400 1 0 4 255 0 0 255" << std::endl;
	  of1 << err1.size() << std::endl; 
	  for (size_t kh=0; kh<err1.size(); ++kh)
	    of1 << err1[kh] << std::endl;
	}
      if (err2.size() > 0)
	{
	  of2 << "400 1 0 4 0 255 0 255" << std::endl;
	  of2 << err2.size() << std::endl; 
	  for (size_t kh=0; kh<err2.size(); ++kh)
	    of2 << err2[kh] << std::endl;
	}
      if (ok1.size() > 0)
	{
	  of3 << "400 1 0 4 0 0 255 255" << std::endl;
	  of3 << ok1.size() << std::endl; 
	  for (size_t kh=0; kh<ok1.size(); ++kh)
	    of3 << ok1[kh] << std::endl;
	}
      if (ok2.size() > 0)
	{
	  of4 << "400 1 0 4 90 90 75 255" << std::endl;
	  of4 << ok2.size() << std::endl; 
	  for (size_t kh=0; kh<ok2.size(); ++kh)
	    of4 << ok2[kh] << std::endl;
	}
#endif
  if (outsideeps_ > 0)
    avdist_ /= (double)outsideeps_;

// #ifdef _OPENMP
//   double time1 = omp_get_wtime();
//   double time_spent = time1 - time0;
//   std::cout << "time_spent in computeAccuracy: " << time_spent << std::endl;
//   std::cout << "time_spent in computeAccuracyElement: " << time_computeAccuracyElement << std::endl;
// #endif
}


//==============================================================================
void LRSurfApprox::computeAccuracy_omp(vector<Element2D*>& ghost_elems)
//==============================================================================
{
//   // We start the timer.
// #ifdef _OPENMP
//   double time0 = omp_get_wtime();
//   double time_computeAccuracyElement = 0.0;
// #endif

  // Check the accuracy of all data points, element by element
  // Note that only points more distant from the surface than the tolerance
  // are considered in avdist_ 

  // Initiate accuracy information
  maxdist_ = 0.0;
  avdist_ = 0.0;
  avdist_all_ = 0.0;
  outsideeps_ = 0;

  RectDomain rd = srf_->containingDomain();
  int dim = srf_->dimension();
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  LRSplineSurface::ElementMap::const_iterator it;
  int num = srf_->numElements();
  int kj;

  double ghost_fac = 0.8;
  ghost_elems.clear();

  //for (it=srf_->elementsBegin(), kj=0; it != srf_->elementsEnd(); ++it, ++kj)
  vector<LRSplineSurface::ElementMap::const_iterator> elem_iters;
  const int num_elem = srf_->numElements();
  elem_iters.reserve(num_elem);
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
  {
      elem_iters.push_back(it);
  }

#pragma omp parallel default(none) private(kj, it) shared(dim, elem_iters, rd, del, ghost_fac, ghost_elems)
  {
      double av_prev, max_prev;
      int nmb_out_prev;
      double umin, umax, vmin, vmax;
      double max_err;
      double av_err;
      double acc_err;
      int outside;
      double acc_err_sgn;
      double av_err_sgn;
      int nmb_pts;
      int nmb_ghost;
      size_t nb;
      int ki;
      double *curr;
      double dist2;
      Element2D *elem;
      double acc_prev;

#pragma omp for schedule(auto)//guided)//static,8)//runtime)//dynamic,4)
      for (kj = 0; kj < num_elem ; ++kj)
      {
	  it = elem_iters[kj];

	  if (!it->second->hasDataPoints())
	  {
	      // Reset accuracy information in element
	      it->second->resetAccuracyInfo();
	      continue;   // No points in which to check accuracy
	  }

	  umin = it->second->umin();
	  umax = it->second->umax();
	  vmin = it->second->vmin();
	  vmax = it->second->vmax();
	  vector<double>& points = it->second->getDataPoints();
	  vector<double>& ghost_points = it->second->getGhostPoints();
	  nmb_pts = it->second->nmbDataPoints();
	  nmb_ghost = it->second->nmbGhostPoints();

	  // Local error information
	  max_err = 0.0;
	  av_err = 0.0;
	  acc_err = 0.0;
	  outside = 0;
	  acc_err_sgn = 0.0;
	  av_err_sgn = 0.0;

	  // Check if the accuracy can have been changed
	  const vector<LRBSpline2D*>& bsplines = it->second->getSupport();
	  for (nb=0; nb<bsplines.size(); ++nb)
	      if (!bsplines[nb]->coefFixed())
		  break;

	  if (/*useMBA_ ||*/ nb < bsplines.size())
	  {
	      // Compute distances in data points and update parameter pairs
	      // if requested
// #ifdef _OPENMP
// 	    double time0_part = omp_get_wtime();
// #endif
	      if (nmb_pts > 0)
	      {
		  computeAccuracyElement(points, nmb_pts, del, rd, it->second.get());
	      }
	  
	      // Compute distances in ghost points
	      if (nmb_ghost > 0 && !useMBA_)
	      {
		  computeAccuracyElement(ghost_points, nmb_ghost, del, rd, it->second.get());
	      }
// #ifdef _OPENMP
// 	    double time1_part = omp_get_wtime();
// 	    time_computeAccuracyElement += time1_part - time0_part;
// #endif
	  }

	  // Accumulate error information related to data points
	  for (ki=0, curr=&points[0]; ki<nmb_pts;)
	  {
	      Point curr_pt(curr+(dim==3)*2, curr+del-1);

	      // Accumulate approximation error
	      dist2 = fabs(curr[del-1]);
	      maxdist_ = std::max(maxdist_, dist2);
	      max_err = std::max(max_err, dist2);
	      acc_err += dist2;
	      acc_err_sgn += curr[del-1];
	      avdist_all_ += dist2;
	      if (dist2 > aepsge_)
	      {
		  av_err_sgn += curr[del-1];
		  avdist_ += dist2;
		  outsideeps_++;
		  av_err += dist2;
		  outside++;
		  
	      }
	      else
	      {
	      }	     	  

	      if (dim == 3 && repar_)
	      {
		  // Check if the point has moved
		  if (curr[0] < umin || curr[0] > umax || curr[1] < vmin || curr[1] > vmax)
		  {
		      // Find element
		      elem = srf_->coveringElement(curr[0], curr[1]);
		      elem->addDataPoints(points.begin()+ki*del, 
					  points.begin()+(ki+1)*del, false);
		      it->second->eraseDataPoints(points.begin()+ki*del, 
						  points.begin()+(ki+1)*del);
		      nmb_pts--;
		  }
		  else
		  {
		      curr += del;
		      ki++;
		  }
	      }
	      else
	      {
		  curr += del;
		  ki++;
	      }
	  }
	  if (outside > 0)
	  {
	      av_err /= (double)outside;
	      av_err_sgn /= (double)outside;
	  }

	  // Previous accuracy information
	  acc_prev = it->second->getAccumulatedError();
	  it->second->getAccuracyInfo(av_prev, max_prev, nmb_out_prev);

	  if (max_err > aepsge_ && max_prev > 0.0 && max_err > ghost_fac*max_prev &&
	      nmb_ghost > 0.25*nmb_pts)
	  {
	      // Collect element for update of ghost points
#pragma omp critical
	      ghost_elems.push_back(it->second.get());
	  }

	  // Store updated accuracy information in the element
	  it->second->setAccuracyInfo(acc_err, av_err, max_err, outside);

      }
  }

  avdist_all_ /= (double)nmb_pts_;
  if (outsideeps_ > 0)
    avdist_ /= (double)outsideeps_;

// #ifdef _OPENMP
//   double time1 = omp_get_wtime();
//   double time_spent = time1 - time0;
//   std::cout << "time_spent in computeAccuracy: " << time_spent << std::endl;
//   std::cout << "time_spent in computeAccuracyElement: " << time_computeAccuracyElement << std::endl;
// #endif
}

//==============================================================================
  void LRSurfApprox::computeAccuracyElement(vector<double>& points, int nmb, int del,
					    RectDomain& rd, const Element2D* elem)
//==============================================================================
{
  int ki, kj, kr;
  double *curr;
  int dim = srf_->dimension();
  // double umin = srf_->paramMin(XFIXED);
  // double umax = srf_->paramMax(XFIXED);
  // double vmin = srf_->paramMin(YFIXED);
  // double vmax = srf_->paramMax(YFIXED);
  int maxiter = 3; //4;
  Element2D* elem2 = (Element2D*)elem;

  // Fetch basis functions
  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  const int nmb_bsplines = (int)bsplines.size();
  double bval, sfval;

  vector<double> grid_height;
  double elem_grid_start[2];
  int grid1, grid2, grid3, grid4;
  if (grid_)
    {
      // Compute height in grid cells corresponding to this element.
      // First identify relevenat grid cells
      double elmin_u = elem->umin();
      double elmax_u = elem->umax();
      double elmin_v = elem->vmin();
      double elmax_v = elem->vmax();
      grid1 = std::max(0, (int)((elmin_u - grid_start_[0])/cell_size_[0]));
      grid2 = (int)((elmax_u - grid_start_[0])/cell_size_[0]);
      grid3 = std::max(0, (int)((elmin_v - grid_start_[1])/cell_size_[1]));
      grid4 = (int)((elmax_v - grid_start_[1])/cell_size_[1]);
      if (grid_start_[0] + grid2*cell_size_[0] < elmax_u)
	grid2++;
      if (grid_start_[1] + grid4*cell_size_[1] < elmax_v)
	grid4++;
      grid_height.resize((grid2 - grid1 + 1)*(grid4 - grid3 + 1));
      double upar, vpar, vpar1;
      Point pos;
      elem_grid_start[0] = grid_start_[0] + grid1*cell_size_[0];
      elem_grid_start[1] = grid_start_[1] + grid3*cell_size_[1];

      // Restrict the grid to the current element. This may give a larger
      // computed error, but avoids eccessive execution times
      for (kr=0, kj=grid3, vpar=elem_grid_start[1]; kj<=grid4; 
	   ++kj, vpar+=cell_size_[1])
	{
	  vpar1 = std::max(elmin_v, std::min(vpar, elmax_v));
	  for (ki=grid1, upar=elem_grid_start[0]; ki<=grid2; 
	       ++ki, ++kr, upar+=cell_size_[0])
	    {
	      srf_->point(pos, std::max(elmin_u, std::min(upar, elmax_u)), 
			  vpar1, elem2);
	      grid_height[kr] = pos[0];
	    }
	}
    }

  int idx1, idx2, sgn;
  double dist, dist1, dist2, dist3, dist4, upar, vpar;
  Point close_pt, vec, norm, pos, curr_pt;
  const int num_threads = 8;
  const int dyn_div = nmb/num_threads;

    for (ki=0, curr=&points[0]; ki<nmb; ++ki, curr+=del)
    {
      curr_pt = Point(curr+(dim==3)*2, curr+del-1);
      if (check_close_ && dim == 3)
	{
	  // Compute closest point
	  // VSK. 052013. Should be changed to make use of the fact that we know
	  // the element (at least initially)
	  // double upar, vpar;
	  // Point close_pt;
	  srf_->setCurrentElement((Element2D*)elem);
	  srf_->closestPoint(curr_pt, upar, vpar, close_pt,
			     dist, aepsge_, maxiter, elem2, &rd, curr);
	  vec = curr_pt - close_pt;
	  // Point norm;
	  srf_->normal(norm, upar, vpar);
	  if (vec*norm < 0.0)
	    dist *= -1;
	  if (to3D_ >= 0)
	    dist = curr[del-2]-close_pt[2];
	}
      else
	{
	  if (grid_)
	    {
	      // Identify grid cell
	      idx1 = (int)((curr[0] - elem_grid_start[0])/cell_size_[0]);
	      idx2 = (int)((curr[1] - elem_grid_start[1])/cell_size_[1]);
	      
	      // Check distance in grid corners
	      dist1 = curr[2]-grid_height[idx2*(grid2-grid1+1)+idx1];
	      dist2 = curr[2]-grid_height[idx2*(grid2-grid1+1)+idx1+1];
	      dist3 = 
		curr[2]-grid_height[(idx2+1)*(grid4-grid3+1)+idx1];
	      dist4 = 
		curr[2]-grid_height[(idx2+1)*(grid4-grid3+1)+idx1+1];
	      
	      // Select minimum distance
	      if (dist1*dist2<0.0 || dist1*dist3<0.0 || dist1*dist4<0.0 || 
		  dist2*dist3<0.0 || dist2*dist4<0.0 || dist3*dist4<0.0)
		dist = 0.0;
	      else
		{
		  sgn = (dist1 >= 0.0) ? 1 : -1;
		  dist = std::min(std::min(fabs(dist1),fabs(dist2)),
				  std::min(fabs(dist3),fabs(dist4)));
		  dist *= sgn;
		}
	    }
	  else
	    {
	      // Evaluate
	      if (dim == 1)
		{
		  // Point pos;
		  // srf_->point(pos, curr[0], curr[1], elem);
		  sfval = 0.0;
		  for (kr=0; kr<nmb_bsplines; ++kr)
		    {
		      bsplines[kr]->evalpos(curr[0], curr[1], &bval);
		      sfval += bval;
		    }
	      
		  dist = curr[2] - sfval;
		  //dist = curr[2] - pos[0];
		}
	      else
		{
		  Point pos;
		  srf_->point(pos, curr[0], curr[1], elem2);
		  dist = pos.dist(Point(curr+2, curr+del));
		  vec = curr_pt - pos;
		  // Point norm;
		  srf_->normal(norm, curr[0], curr[1], elem2);
		  if (vec*norm < 0.0)
		    dist *= -1;
		}
	    }
	}
      curr[del-1] = dist;
  }

}


//==============================================================================
void LRSurfApprox::computeAccuracyElement_omp(vector<double>& points, int nmb, int del,
					      RectDomain& rd, const Element2D* elem)
//==============================================================================
{
  int ki, kj, kr;
  double *curr;
  int dim = srf_->dimension();
  // double umin = srf_->paramMin(XFIXED);
  // double umax = srf_->paramMax(XFIXED);
  // double vmin = srf_->paramMin(YFIXED);
  // double vmax = srf_->paramMax(YFIXED);
  int maxiter = 3; //4;
  Element2D* elem2 = (Element2D*)elem;

  // Fetch basis functions
  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  const int nmb_bsplines = (int)bsplines.size();
  double bval, sfval;

  vector<double> grid_height;
  double elem_grid_start[2];
  int grid1, grid2, grid3, grid4;
  if (grid_)
    {
      // Compute height in grid cells corresponding to this element.
      // First identify relevenat grid cells
      double elmin_u = elem->umin();
      double elmax_u = elem->umax();
      double elmin_v = elem->vmin();
      double elmax_v = elem->vmax();
      grid1 = std::max(0, (int)((elmin_u - grid_start_[0])/cell_size_[0]));
      grid2 = (int)((elmax_u - grid_start_[0])/cell_size_[0]);
      grid3 = std::max(0, (int)((elmin_v - grid_start_[1])/cell_size_[1]));
      grid4 = (int)((elmax_v - grid_start_[1])/cell_size_[1]);
      if (grid_start_[0] + grid2*cell_size_[0] < elmax_u)
	grid2++;
      if (grid_start_[1] + grid4*cell_size_[1] < elmax_v)
	grid4++;
      grid_height.resize((grid2 - grid1 + 1)*(grid4 - grid3 + 1));
      double upar, vpar, vpar1;
      Point pos;
      elem_grid_start[0] = grid_start_[0] + grid1*cell_size_[0];
      elem_grid_start[1] = grid_start_[1] + grid3*cell_size_[1];

      // Restrict the grid to the current element. This may give a larger
      // computed error, but avoids eccessive execution times
      for (kr=0, kj=grid3, vpar=elem_grid_start[1]; kj<=grid4; 
	   ++kj, vpar+=cell_size_[1])
	{
	  vpar1 = std::max(elmin_v, std::min(vpar, elmax_v));
	  for (ki=grid1, upar=elem_grid_start[0]; ki<=grid2; 
	       ++ki, ++kr, upar+=cell_size_[0])
	    {
	      srf_->point(pos, std::max(elmin_u, std::min(upar, elmax_u)), 
			  vpar1, elem2);
	      grid_height[kr] = pos[0];
	    }
	}
    }

  int idx1, idx2, sgn;
  double dist, dist1, dist2, dist3, dist4, upar, vpar;
  Point close_pt, vec, norm, pos, curr_pt;
  const int num_threads = 8;
  const int dyn_div = nmb/num_threads;

#ifdef _OPENMP
  pthread_attr_t attr;
  size_t stacksize;
  pthread_attr_getstacksize(&attr, &stacksize);
  //	std::cout << "stacksize (in MB): " << (double)stacksize/(1024.0*1024.0) << std::endl;
#endif
  //	omp_set_num_threads(4);
#pragma omp parallel default(none) private(ki, curr, idx1, idx2, dist, upar, vpar, close_pt, curr_pt, vec, norm, dist1, dist2, dist3, dist4, sgn, pos, sfval, kr, kj, bval) \
  shared(points, nmb, del, dim, rd, maxiter, elem_grid_start, grid2, grid1, grid_height, grid3, grid4, elem2, bsplines)
#pragma omp for schedule(dynamic, 4)//static, 4)//runtime)//guided)//auto)
  for (ki=0; ki<nmb; ++ki)
    {
      // #ifdef _OPENMP
      //       const int num_omp_threads = omp_get_num_threads();
      // 	if (num_omp_threads > 4)
      // 	{
      // 	    ;//printf("omp_get_num_threads(): %d\n",omp_get_num_threads());
      // 	    // std::cout << "num_omp_threads: " << num_omp_threads << std::endl;
      // 	}
      // #endif
      curr = &points[ki*del];
      curr_pt = Point(curr+(dim==3)*2, curr+del-1);
      if (check_close_ && dim == 3)
	{
	  // Compute closest point
	  // VSK. 052013. Should be changed to make use of the fact that we know
	  // the element (at least initially)
	  // double upar, vpar;
	  // Point close_pt;
	  srf_->closestPoint(curr_pt, upar, vpar, close_pt,
			     dist, aepsge_, maxiter, elem2, &rd, curr);
	  vec = curr_pt - close_pt;
	  // Point norm;
	  srf_->normal(norm, upar, vpar, elem2);
	  if (vec*norm < 0.0)
	    dist *= -1;
	  if (to3D_ >= 0)
	    dist = curr[del-2]-close_pt[2];
	}
      else
	{
	  if (grid_)
	    {
	      // Identify grid cell
	      idx1 = (int)((curr[0] - elem_grid_start[0])/cell_size_[0]);
	      idx2 = (int)((curr[1] - elem_grid_start[1])/cell_size_[1]);
	      
	      // Check distance in grid corners
	      dist1 = curr[2]-grid_height[idx2*(grid2-grid1+1)+idx1];
	      dist2 = curr[2]-grid_height[idx2*(grid2-grid1+1)+idx1+1];
	      dist3 = 
		curr[2]-grid_height[(idx2+1)*(grid4-grid3+1)+idx1];
	      dist4 = 
		curr[2]-grid_height[(idx2+1)*(grid4-grid3+1)+idx1+1];
	      
	      // Select minimum distance
	      if (dist1*dist2<0.0 || dist1*dist3<0.0 || dist1*dist4<0.0 || 
		  dist2*dist3<0.0 || dist2*dist4<0.0 || dist3*dist4<0.0)
		dist = 0.0;
	      else
		{
		  sgn = (dist1 >= 0.0) ? 1 : -1;
		  dist = std::min(std::min(fabs(dist1),fabs(dist2)),
				  std::min(fabs(dist3),fabs(dist4)));
		  dist *= sgn;
		}
	    }
	  else
	    {
	      // Evaluate
	      if (dim == 1)
		{
		  // Point pos;
		  // srf_->point(pos, curr[0], curr[1], elem);
		  sfval = 0.0;
		  for (kr=0; kr<nmb_bsplines; ++kr)
		    {
		      bsplines[kr]->evalpos(curr[0], curr[1], &bval);
		      sfval += bval;
		    }
	      
		  dist = curr[2] - sfval;
		  //dist = curr[2] - pos[0];
		}
	      else
		{
		  Point pos;
		  srf_->point(pos, curr[0], curr[1], elem2);
		  dist = pos.dist(Point(curr+2, curr+del));
		  vec = curr_pt - pos;
		  // Point norm;
		  srf_->normal(norm, curr[0], curr[1]);
		  if (vec*norm < 0.0)
		    dist *= -1;
		}
	    }
	}
      curr[del-1] = dist;
    }

}


//==============================================================================
int LRSurfApprox::refineSurf()
//==============================================================================
{
#ifdef DEBUG
  std::ofstream of0("element_info.dat");
  int idx=0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (it->second->hasAccuracyInfo())
  	{
  	  double av_err, max_err;
  	  int nmb_out;
  	  it->second->getAccuracyInfo(av_err, max_err, nmb_out);
  	  of0 << "Element " << idx << ", domain: [(" << it->second->umin() << ",";
  	  of0 << it->second->umax() << ")x(" << it->second->vmin() << ",";
  	  of0 << it->second->vmax() << ")]" << std::endl;
  	  of0 << "Average error: " << av_err << std::endl;
  	  of0 << "Maximum error: " << max_err << std::endl;
  	  of0 << "Number of points: " << it->second->nmbDataPoints();
  	  of0 << ", number of error points: " << nmb_out << std::endl;
  	}
      ++idx;
    }
  std::cout << "Number of elements: " << srf_->numElements() << std::endl;
  std::cout << "Maxdist= " << maxdist_ << ", avdist= " << avdist_;
  std::cout << ", nmb out= " << outsideeps_ << std::endl;
#endif

  int choice = 1;  // Strategy for knot insertion in one single B-spline

  // Construct indexed bspline array and collect related accuracy information
  int num_bspl = srf_->numBasisFunctions();
  vector<LRBSpline2D*> bsplines(num_bspl);
  vector<double> error(num_bspl, 0.0);
  vector<double> max_error(num_bspl, 0.0);
  vector<double> av_error(num_bspl, 0.0);
  vector<double> domain_size(num_bspl);
  vector<int> num_pts(num_bspl, 0);
  vector<int> num_out_pts(num_bspl, 0); 
  double mean_err = 0.0;
  size_t kr = 0;
  for (LRSplineSurface::BSplineMap::const_iterator it=srf_->basisFunctionsBegin();
       it != srf_->basisFunctionsEnd(); ++it)
    {
      LRBSpline2D* curr = it->second.get();

      for (auto it2=curr->supportedElementBegin(); 
	   it2 != curr->supportedElementEnd(); ++it2)
	{
	  num_pts[kr] += (*it2)->nmbDataPoints();
	  num_out_pts[kr] += (*it2)->getNmbOutsideTol();
	  error[kr] += (*it2)->getAccumulatedError();
	  max_error[kr] = std::max(max_error[kr], (*it2)->getMaxError());
	  av_error[kr] += (*it2)->getAverageError();  // Only counting those 
	  // points being outside of the tolerance
	}
      av_error[kr] /= (double)(curr->nmbSupportedElements());

      // Use sqrt to reduce the significance of this property compared to the
      // error
      domain_size[kr] = sqrt((curr->umax()-curr->umin())*(curr->vmax()-curr->vmin()));
      bsplines[kr] = curr;
      mean_err += av_error[kr++];
    }
  mean_err /= (double)num_bspl;

  // Sort bsplines according to average error weighted with the domain size
  vector<int> bspl_perm(num_bspl);
  int ki, kj;
  for (ki=0; ki<num_bspl; ++ki)
    bspl_perm[ki] = ki;

  // Do the sorting
  int group_fac = 3;
  double error_fac = 0.1;
  double error_fac2 = 10.0;
  for (ki=0; ki<num_bspl; ++ki)
    for (kj=ki+1; kj<num_bspl; ++kj)
      {
	double curr_err1 = error[bspl_perm[ki]]*domain_size[bspl_perm[ki]];
	double curr_err2 = error[bspl_perm[kj]]*domain_size[bspl_perm[kj]];

	if (num_out_pts[bspl_perm[ki]] > group_fac ||
	    (double)num_out_pts[bspl_perm[ki]] > 
	    error_fac*((double)num_pts[bspl_perm[ki]]))
	  curr_err1 *= error_fac2;
	if (num_out_pts[bspl_perm[kj]] > group_fac ||
	    (double)num_out_pts[bspl_perm[kj]] > 
	    error_fac*((double)num_pts[bspl_perm[kj]]))
	  curr_err2 *= error_fac2;

	// Modify if there is a significant number of large error points 
	if (curr_err1 < curr_err2)
	  std::swap(bspl_perm[ki], bspl_perm[kj]);
      }
  
  // Split the most important B-splines, but only if the maximum
  // error is larger than the tolerance
  //double fac = 0.5;
  int nmb_perm = (int)bspl_perm.size();
  int nmb_split = (int)(0.5*nmb_perm);
  //nmb_split = std::min(nmb_split, 600);  // Limit the number of refinements
  int min_nmb_pts = 1; //4;
  //double pnt_fac = 0.2;
  //int min_nmb_out = 4;

  vector<LRSplineSurface::Refinement2D> refs;
  int nmb_refs = 0;

  int nmb_fixed = 0;
  for (kr=0; kr<bspl_perm.size(); ++kr)
    {
      if (max_error[bspl_perm[kr]] < aepsge_)
	{
	  // TEST. Keep coefficient fixed
	  bsplines[bspl_perm[kr]]->setFixCoef(1);
	  nmb_fixed++;
	  continue;
	}

      // Do not split B-splines with too few points in its domain
      if (num_pts[bspl_perm[kr]] < min_nmb_pts)
	continue;

      if (nmb_refs >= nmb_split)
	break;

      //if (av_error[bspl_perm[kr]] < fac*mean_err)
      if (false /*av_error[bspl_perm[kr]] < fac*mean_err &&
	  (num_out_pts[bspl_perm[kr]] < (int)(pnt_fac*num_pts[bspl_perm[kr]]) ||
	  num_pts[bspl_perm[kr]] < min_nmb_out)*/)
	continue;  // Do not split this B-spline at this stage

      nmb_refs++;  // Split this B-spline
      
      // How to split					
      defineRefs(bsplines[bspl_perm[kr]], refs, choice);
    }
  
#ifdef DEBUG
  std::ofstream of("refine0.dat");
  //std::streamsize prev = of.precision(15);
  (void)of.precision(15);
  for (kr=0; kr<refs.size(); ++kr)
    {
      of << refs[kr].kval << "  " << refs[kr].start << "  " << refs[kr].end;
      of << "  " << refs[kr].d << "  " << refs[kr].multiplicity << std::endl;
    }
  std::cout << "Number of refinements: " << refs.size() << std::endl;
  std::cout << "Number of coef fixed: " << nmb_fixed << std::endl;
#endif

  for (kr=0; kr<refs.size(); ++kr)
    {
#ifdef DEBUG
      // std::cout << "Refine nr " << kr << ": " << refs[kr].kval << "  " << refs[kr].start << "  ";
      // std::cout << refs[kr].end << "  " << refs[kr].d << "  " << refs[kr].multiplicity << std::endl;
#endif
      // Perform refinements, one at the time to keep information stored in the elements
      srf_->refine(refs[kr], true /*false*/);
#ifdef DEBUG
      // std::ofstream of2("refined2_sf.g2");
      // srf_->writeStandardHeader(of2);
      // srf_->write(of2);
      // of2 << std::endl;

      // std::ofstream of3("refined3_sf.g2");
      // shared_ptr<LRSplineSurface> tmp2;
      // if (srf_->dimension() == 1)
      // 	{
      // 	  tmp2 = shared_ptr<LRSplineSurface>(srf_->clone());
      // 	  tmp2->to3D();
      // 	}
      // else
      // 	tmp2 = srf_;
      // tmp2->writeStandardHeader(of3);
      // tmp2->write(of3);
      // of3 << std::endl;

      //   // For all elements, check that the sum of scaled B-splines in the 
      // // midpoint is 1
      // double tol = 1.0e-10;
      // for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
      // 	   it != srf_->elementsEnd(); ++it)
      // 	{
      // 	  double upar = 0.5*(it->second->umin() + it->second->umax());
      // 	  double vpar = 0.5*(it->second->vmin() + it->second->vmax());
      // 	  double val = it->second->sumOfScaledBsplines(upar, vpar);
      // 	  if (fabs(1.0-val) > tol)
      // 	    {
      // 	      std::cout << "B-splines in element (" << upar <<  "," << vpar;
      // 	      std::cout << ") do not sum up to one" << std::endl;
      // 	    }
      // 	}
      int stop_break = 1;
#endif
    }
#ifdef DEBUG
  std::ofstream ofmesh("mesh1.eps");
  writePostscriptMesh(*srf_, ofmesh);
#endif

  return (int)refs.size();
}

//==============================================================================
void LRSurfApprox::refineSurf2()
//==============================================================================
{
#ifdef DEBUG
  int idx=0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (it->second->hasAccuracyInfo())
  	{
  	  double av_err, max_err;
  	  int nmb_out;
  	  it->second->getAccuracyInfo(av_err, max_err, nmb_out);
  	  std::cout << "Element " << idx << ", domain: [(" << it->second->umin() << ",";
  	  std::cout << it->second->umax() << ")x(" << it->second->vmin() << ",";
  	  std::cout << it->second->vmax() << ")]" << std::endl;
  	  std::cout << "Average error: " << av_err << std::endl;
  	  std::cout << "Maximum error: " << max_err << std::endl;
  	  std::cout << "Number of points: " << it->second->nmbDataPoints();
  	  std::cout << ", number of error points: " << nmb_out << std::endl;
  	}
      ++idx;
    }
#endif

  // Construct indexed element array
  int num_el = srf_->numElements();
  vector<Element2D*> elem(num_el);
  size_t kr = 0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    elem[kr++] = it->second.get();
  
  // Sort elements according to average error
  vector<int> el_perm(num_el);
  int ki, kj;
  for (ki=0; ki<num_el; ++ki)
    el_perm[ki] = ki;

  // Do the sorting
  for (ki=0; ki<num_el; ++ki)
    for (kj=ki+1; kj<num_el; ++kj)
      {
	double av_err1 = elem[el_perm[ki]]->getAverageError();
 	double av_err2 = elem[el_perm[kj]]->getAverageError();

	if (av_err1 < av_err2)
	  std::swap(el_perm[ki], el_perm[kj]);
      }
  
  // Define threshhold for refinement. The average error includes only those points that
  // are outside of the resolution
  double max_av = 0;
  double mean_av = 0.0;
  int num_err = 0;
  for (ki=0; ki<num_el; ++ki)
    {
      double av_err = elem[el_perm[ki]]->getAverageError();
      max_av = std::max(max_av, av_err);
      mean_av += av_err;
      if (elem[el_perm[ki]]->getMaxError() > aepsge_)
	num_err++;
    }
  mean_av /= (double)num_el;

  double threshhold = 0.25*mean_av;  // Need to experiment with this

  vector<LRSplineSurface::Refinement2D> refs;
  int nmb_ref = std::max(std::min((int)el_perm.size(), 4), (int)(0.75*num_err));
  
  for (kr=0; kr<el_perm.size(); )
    {
      size_t nmb_perm = el_perm.size();
      if (elem[el_perm[kr]]->getAverageError() < threshhold && (int)refs.size() > nmb_ref)
	break;  // No more refinements at the current stage
      if (elem[el_perm[kr]]->getMaxError() < aepsge_)
	break;

      // Check feasability of split
      //size_t nmb_refs = refs.size();
      vector<Element2D*> elements;  // Elements affected by the refinement(s)
      checkFeasibleRef(elem[el_perm[kr]], refs, elements);
      if (elements.size() > 0)
	{
	  // Remove affected elements from pool
	  for (size_t kh=0; kh<elements.size(); ++kh)
	    {
	      for (ki=0; ki<(int)el_perm.size();)
		{
		  if (elements[kh] == elem[el_perm[ki]])
		    {
		      el_perm.erase(el_perm.begin()+ki);
		      break;
		    }
		  else
		    ki++;
		}
	    }
	  if (el_perm.size() == nmb_perm)
	    kr++;
	}
      else
	kr++;
    }
  
#ifdef DEBUG
  std::ofstream of("refine0.dat");
  (void)of.precision(15);
  for (kr=0; kr<refs.size(); ++kr)
    {
      of << refs[kr].kval << "  " << refs[kr].start << "  " << refs[kr].end;
      of << "  " << refs[kr].d << "  " << refs[kr].multiplicity << std::endl;
    }
#endif

  for (kr=0; kr<refs.size(); ++kr)
    {
#ifdef DEBUG
      std::cout << "Refine nr " << kr << ": " << refs[kr].kval << "  " << refs[kr].start << "  ";
      std::cout << refs[kr].end << "  " << refs[kr].d << "  " << refs[kr].multiplicity << std::endl;
#endif
      // Perform refinements, one at the time to keep information stored in the elements
      srf_->refine(refs[kr], true /*false*/);
      std::ofstream of2("refined2_sf.g2");
      srf_->writeStandardHeader(of2);
      srf_->write(of2);
      of2 << std::endl;
      int stop_break = 1;
     }

  // // Update coef_known from information in LR B-splines
  // //updateCoefKnown();
  // unsetCoefKnown();
  //  if (fix_boundary_)
  //   setFixBoundary(true);

  
  // // Prelimenary implementation. Write accuracy information to standard output
  // // and get refinement informatiation
  // // For each element, check if any data points are stored
  // int idx=0;
  // for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
  //      it != srf_->elementsEnd(); ++it)
  //   {
  //     if (it->second->hasAccuracyInfo())
  // 	{
  // 	  double av_err, max_err;
  // 	  int nmb_out;
  // 	  it->second->getAccuracyInfo(av_err, max_err, nmb_out);
  // 	  std::cout << "Element " << idx << ", domain: [(" << it->second->umin() << ",";
  // 	  std::cout << it->second->umax() << ")x(" << it->second->vmin() << ",";
  // 	  std::cout << it->second->vmax() << ")]" << std::endl;
  // 	  std::cout << "Average error: " << av_err << std::endl;
  // 	  std::cout << "Maximume error: " << max_err << std::endl;
  // 	  std::cout << "Number of points: " << it->second->nmbDataPoints();
  // 	  std::cout << ", number of error points: " << nmb_out << std::endl;
  // 	}
  //     ++idx;
  //   }

  // // Specify refinement
  // while (true)
  //   {
  //     int more;
  //     std::cout << "Refine more? " << std::endl;
  //     std::cin >> more;
  //     if (!more)
  // 	break;
  //     double parval, start, end;
  //     int dir;
  //     int mult = 1;
  //     std::cout << "New knot, give pardir (0,1), fixed value, start and end: ";
  //     std::cout << std::endl;
  //     std::cin >> dir;
  //     std::cin >> parval;
  //     std::cin >> start;
  //     std::cin >> end;

  //     // The refinement must be performed inserting one knot segment at the
  //     // time in order to maintain information stored in the element
  //     srf_->refine((dir==0) ? YFIXED : XFIXED, parval, start, end, mult);

  //     // Update coef_known from information in LR B-splines
  //     updateCoefKnown();
  //   }
  
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim)
//==============================================================================
{
  // Create a cubic Bezier surface represented as an LR B-spline surface
  // Compute domain
  double umin, umax, vmin, vmax;
  computeParDomain(dim, umin, umax, vmin, vmax);

  vector<double> knots_u(8);
  vector<double> knots_v(8);

  int kj;
  for (kj=0; kj<4; ++kj)
    {
      knots_u[kj] = umin;
      knots_v[kj] = vmin;
      knots_u[kj+4] = umax;
      knots_v[kj+4] = vmax;
    }

  makeInitSurf(dim, 4, 4, 4, 4, &knots_u[0], &knots_v[0]);
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
				int order_v)
//==============================================================================
{
  // Compute domain
  double domain[4]; //umin, umax, vmin, vmax;
  computeParDomain(dim, domain[0], domain[1], domain[2], domain[3]);

  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v, domain);
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
				int order_v, double domain[4])
//==============================================================================
{
  vector<double> knots_u(ncoef_u+order_u);
  vector<double> knots_v(ncoef_v+order_v);

  double umin, umax, vmin, vmax;
  computeParDomain(dim, umin, umax, vmin, vmax);
  umin = std::min(umin, domain[0]);
  umax = std::max(umax, domain[1]);
  vmin = std::min(vmin, domain[2]);
  vmax = std::max(vmax, domain[3]);

  int kj;
  double del_u = (umax - umin)/(double)(ncoef_u-order_u+1);
  double del_v = (vmax - vmin)/(double)(ncoef_v-order_v+1);
  for (kj=0; kj<order_u; ++kj)
    knots_u[kj] = umin;
   for (; kj<ncoef_u; ++kj)
    knots_u[kj] = umin + (kj-order_u+1)*del_u;
  for (; kj<ncoef_u+order_u; ++kj)
    knots_u[kj] = umax;
  for (kj=0; kj<order_v; ++kj)
    knots_v[kj] = vmin;
  for (; kj<ncoef_v; ++kj)
    knots_v[kj] = vmin + (kj-order_v+1)*del_v;
  for (; kj<ncoef_v+order_v; ++kj)
    knots_v[kj] = vmax;

  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v, &knots_u[0], &knots_v[0]);
}

//==============================================================================
void LRSurfApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
				int order_v, double *knots_u, double *knots_v)
//==============================================================================
{
  int nmb = (int)points_.size()/(dim+2);
  shared_ptr<SplineSurface> result_surf;
  if (initMBA_)
    {
      vector<double> coefs(ncoef_u*ncoef_v*dim, initMBA_coef_);
      result_surf = shared_ptr<SplineSurface>(new SplineSurface(ncoef_u, ncoef_v, 
								order_u, order_v,  
								knots_u,     // ptr to knots in u
								knots_v,     // ptr to knots in v
								&coefs[0], // ptr to coefs
								dim));
    }
  else
    result_surf = createSurf(&points_[0], nmb, 
			     dim, ncoef_u, order_u,
			     ncoef_v, order_v, knots_u,
			     knots_v, smoothweight_,
			     maxdist_, avdist_, outsideeps_);


  // Make LR spline surface
  double knot_tol = 1.0e-6;

  srf_ = shared_ptr<LRSplineSurface>(new LRSplineSurface(result_surf.get(), knot_tol));
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}

//==============================================================================
shared_ptr<SplineSurface> LRSurfApprox::createSurf(double* points, int nmb_pts,
						   int dim, int ncoef_u, int order_u, 
						   int ncoef_v, int order_v, 
						   double *knots_u, double *knots_v,
						   double smoothweight,
						   double& maxdist, double& avdist,
						   int& nmb_outside)
//==============================================================================
{
  shared_ptr<SplineSurface> result_surf;
  vector<double> coefs(ncoef_u*ncoef_v*dim, 0.0);
  shared_ptr<SplineSurface> sf1(new SplineSurface(ncoef_u, ncoef_v, 
						  order_u, order_v,  
						  knots_u,     // ptr to knots in u
						  knots_v,     // ptr to knots in v
						  &coefs[0], // ptr to coefs
						  dim));

  // Reformatting data to the format that ApproxSurf wants (separating parameter values and 
  // data points in two distinct vectors).
  vector<double> pts, param;
  pts.reserve(dim*nmb_pts);
  param.reserve(2 * nmb_pts);
  int ki, kj;
  double* it;
  for (it=points, kj=0; kj<nmb_pts; ++kj) 
    {
      for (ki=0; ki<2; ++ki)
	param.push_back(*it++);
      for (ki=0; ki<dim; ++ki)
	pts.push_back(*it++);
    }
    
  // Approximate data points to create initial surface
  SmoothSurf asurf;  // Engine for least squares approximation with smoothing
  int stat = 0;
  int seem[2];       
  seem[0] = seem[1] = 0;  // Not a closed surface

  // Define weights
  int min_der = std::max(3, std::min(order_u, order_v)-1);
  double wgt1 = 0.0;
  double wgt3 = (min_der >= 3) ? 0.5*smoothweight : 0.0;
  double wgt2 = (1.0 - wgt3)*smoothweight;
  wgt3 *= smoothweight;
  double approxweight = 1.0 - wgt1 - wgt2 - wgt3;
  std::vector<double> pt_weight(nmb_pts, 1.0);
 
  // Prepare for approximation
  vector<int> coef_known(ncoef_u*ncoef_v, 0);
  asurf.attach(sf1, seem, &coef_known[0], 0, false);

  if (smoothweight_ > 0.0)
    asurf.setOptimize(wgt1, wgt2, wgt3);

  asurf.setLeastSquares(pts, param,
			 pt_weight, approxweight);
  
  // Approximate
  stat = asurf.equationSolve(result_surf);
  return result_surf;
}

//==============================================================================
void LRSurfApprox::makeInitSurf(shared_ptr<SplineSurface> surf)
//==============================================================================
{
  // Make LR spline surface
  double knot_tol = 1.0e-6;

  srf_ = shared_ptr<LRSplineSurface>(new LRSplineSurface(surf.get(), knot_tol));
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}

//==============================================================================
void LRSurfApprox::computeParDomain(int dim, double& umin, double& umax, 
				    double& vmin, double& vmax)
//==============================================================================
{
  // Compute domain
  umin = umax = points_[0];
  vmin = vmax = points_[1];
  int del = 2 + dim;
  for (size_t ki=del; ki<points_.size(); ki+=del)
    {
      umin = std::min(umin, points_[ki]);
      umax = std::max(umax, points_[ki]);
      vmin = std::min(vmin, points_[ki+1]);
      vmax = std::max(vmax, points_[ki+1]);
    }
}

//==============================================================================
void LRSurfApprox::setCoefKnown()
//==============================================================================
{
  if (fix_corner_)
    {
      // Fix coefficients in corners
      for (LRSplineSurface::BSplineMap::iterator it=srf_->basisFunctionsBeginNonconst(); 
	   it != srf_->basisFunctionsEndNonconst(); ++it)
	{
	  int deg1 = it->second->degree(XFIXED);
	  int deg2 = it->second->degree(YFIXED);
	  int mult1 = it->second->endmult_u(true);
	  int mult2 = it->second->endmult_u(false);
	  int mult3 = it->second->endmult_v(true);
	  int mult4 = it->second->endmult_v(false);
	  if ((mult1 == deg1+1 || mult2 == deg1+1) &&
	      (mult3 == deg2+1 || mult4 == deg2+1))
	    it->second->setFixCoef(1);
	}
    }

  // Set coef_known from boundary information
  // Note that k-multiple knots at boundaries are expected, otherwise no
  // coefficients will be fixed
  
  for (size_t ki=0; ki<4; ++ki)
    {
      if (edge_derivs_[ki] == 0)
	continue;  // No coefficients to fix
      
      // Set boundary characteristica
      Direction2D d = (ki == 1 || ki == 3) ? XFIXED : YFIXED;  // Orthogonal to the curve
      Direction2D d2 = (ki == 1 || ki == 3) ? YFIXED : XFIXED; // Along the curve
      bool atstart = (ki == 0 || ki == 3);  // Whether the constant
      // parameter curve is in the start related to the opposite parameter direction
      // of the surface

      // Define coefficients as fixed in the associated b-splines
      int fixcoef = 1;
      setCoefKnown(d, d2, atstart, fixcoef);
    }

  // Transfer information to the global vector
  updateCoefKnown();
}

//==============================================================================
void LRSurfApprox::setCoefKnown(Direction2D d, Direction2D d2, 
				bool atstart, int fixcoef)
//==============================================================================
{
#ifdef DEBUG
  std::ofstream of("mesh.eps");
  writePostscriptMesh(*srf_, of);
#endif

  // // Check if the basis has k-tupple knots in the current direction
  // const Mesh2D& mesh = srf_->mesh();
  // int ix = (atstart) ? mesh.firstMeshVecIx(d) : 
  //   mesh.lastMeshVecIx(d);  
  // int mult = mesh.minMultInLine(d, ix);
  // int deg = srf_->degree(d);  // Degree orthogonal to the constant parameter curve
  // if (mult < deg)
  //   return;  // No fixed coefficients are set

  // // Fetch LRBSplines. Traverse the knot vector to make keys 
  // // for the bspline map.
  // int dir = d;
  // int startmult[2], endmult[2];
  // double startval[2], endval[2];
  // // int num = mesh.numDistinctKnots(d);  // Number of knots orthogonal to the
  // // // constant parameter curve
  // int num = mesh.numDistinctKnots(d); //mesh.numDistinctKnots(d2);  // Number of knots orthogonal to the
  // // constant parameter curve

  // // Set parameter value at the constant parameter curve
  // if (atstart)
  //   {
  //     startval[dir] = mesh.kval(d, 0);
  //   }
  // else
  //   {
  //     endval[dir] = mesh.kval(d, mesh.numDistinctKnots(d)-1);
  //   }
  // int deg2 = srf_->degree(d2);  // Degree along the constant parameter curve

  // vector<int> knot_idx =  LRBSpline2DUtils::derive_knots(mesh, d2, 
  // 							 mesh.firstMeshVecIx(d2),
  // 							 mesh.lastMeshVecIx(d2),
  // 							 atstart ? ix : ix-1,
  // 							 atstart ? ix+1 : ix);
  // // vector<int> knot_idx =  LRBSpline2DUtils::derive_knots(mesh, d, 
  // // 							 mesh.firstMeshVecIx(d2),
  // // 							 mesh.lastMeshVecIx(d2),
  // // 							 atstart ? ix : ix-1,
  // // 							 atstart ? ix+1 : ix);

  // // Since we have multiple knots in the surface boundary and we are only interested
  // // in the basis functions being non-zero along the boundary, we know that these
  // // basis functions must have a multiplicity equal to the order at the boundary
  // // and one in the other end.
  // startmult[dir] = atstart ? deg+1 : 1;
  // endmult[dir] = atstart ? 1 : deg+1;
  // size_t k1, k2;
  // vector<double> coefs;
  // for (k1=0, k2=deg2+1; k2<knot_idx.size(); ++k1, ++k2)
  //   {
  //     // Fetch domain of first minimal LR B-spline along the constant
  //     // parameter curve
  //     // First orthogonal to the curve
  //     if (atstart)
  // 	{
  // 	  int k_idx = 
  // 	    Mesh2DUtils::search_upwards_for_nonzero_multiplicity(mesh, d, 1, 
  // 								 knot_idx[k1], knot_idx[k2]);
  // 	  endval[dir] = mesh.kval(d, k_idx);
  // 	}
  //     else
  // 	{
  // 	  int k_idx = 
  // 	    Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh, d, num-2, 
  // 								   knot_idx[k1], knot_idx[k2]);
  // 	  startval[dir] = mesh.kval(d, k_idx);
  // 	}

  //     // Along the curve
  //     startval[1-dir] = mesh.kval(d2, knot_idx[k1]);
  //     endval[1-dir] = mesh.kval(d2, knot_idx[k2]);

  //     // Count multiplicitity along the curve
  //     int km = 1;
  //     for (; km<=deg2; ++km)
  // 	if (knot_idx[k1+km] > knot_idx[k1])
  // 	  break;
  //     startmult[1-dir] = km;

  //     km = 1;
  //     for (; km<=deg2; ++km)
  // 	if (knot_idx[k2-km] < knot_idx[k2])
  // 	  break;
  //     endmult[1-dir] = km;

  //     // Fetch the associated LR B-spline
  //     LRSplineSurface::BSplineMap::iterator bm = srf_->bsplineFromDomain(startval[0], startval[1], 
  // 									 endval[0], endval[1], 
  // 									 startmult[0], startmult[1], 
  // 									 endmult[0], endmult[1]);

  //     // Set fixed coefficient information in the b-spline
  //     bm->second->setFixCoef(fixcoef);
  //   }

  // Traverse all B-splines and check whether they have maximum multiplicity along
  // the given edge
  for (LRSplineSurface::BSplineMap::iterator it=srf_->basisFunctionsBeginNonconst(); 
       it != srf_->basisFunctionsEndNonconst(); ++it)
    {
      int deg = it->second->degree(d);
      int mult = (d == XFIXED) ? it->second->endmult_u(atstart) :
	it->second->endmult_v(atstart);
      if (mult == deg+1)
	it->second->setFixCoef(fixcoef);
    }
  
}

//==============================================================================
void LRSurfApprox::updateCoefKnown()
//==============================================================================
{
  // Fetch all boundary b-splines and check if the coefficient should be fixed
  int fixcoef = 1;
  for (LRSplineSurface::BSplineMap::iterator it=srf_->basisFunctionsBeginNonconst(); 
       it != srf_->basisFunctionsEndNonconst(); ++it)
    {
      int deg1 = it->second->degree(XFIXED);
      int deg2 = it->second->degree(YFIXED);
      int mult1_1 = it->second->endmult_u(true);
      int mult1_2 = it->second->endmult_u(false);
      int mult2_1 = it->second->endmult_v(true);
      int mult2_2 = it->second->endmult_v(false);
      if (mult1_1 == deg1+1 || mult1_2 == deg1+1 || 
	  mult2_1 == deg2+1 || mult2_2 == deg2+1)
	{
	  const vector<Element2D*>& curr_el = it->second->supportedElements();
	  int nmb_pts = 0;
	  for (size_t ki=0; ki<curr_el.size(); ++ki)
	    nmb_pts += curr_el[ki]->nmbDataPoints();
	  if (nmb_pts == 0)
	    it->second->setFixCoef(fixcoef);
	}
    }
  
  coef_known_.resize(srf_->numBasisFunctions());

  LRSplineSurface::BSplineMap::const_iterator it_bs;			
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    coef_known_[ki] = it_bs->second->coefFixed();
}

//==============================================================================
void LRSurfApprox::unsetCoefKnown()
//==============================================================================
{
  coef_known_.resize(srf_->numBasisFunctions());

  LRSplineSurface::BSplineMap::const_iterator it_bs;			
  size_t ki;
  for (it_bs=srf_->basisFunctionsBegin(), ki=0; it_bs!=srf_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    {
      coef_known_[ki] = 0;
      it_bs->second->setFixCoef(0);
    }
}

//==============================================================================
void LRSurfApprox::defineRefs(LRBSpline2D* bspline,
			      vector<LRSplineSurface::Refinement2D>& refs,
			      int choice)
//==============================================================================
{
  // For each alternative (knot span) in each parameter direction, collect
  // accuracy statistic
  // Compute also average element size
  double tol = srf_->getKnotTol();
  int size1 = bspline->degree(XFIXED)+1;
  int size2 = bspline->degree(YFIXED)+1;
  vector<double> u_info(size1, 0.0);
  vector<double> v_info(size2, 0.0);
  vector<double> v_elsize(size1, 0.0);
  vector<double> u_elsize(size2, 0.0);
  
  const vector<int>& kvec_u = bspline->kvec(XFIXED);
  const vector<int>& kvec_v = bspline->kvec(YFIXED);
  const Mesh2D* mesh = bspline->getMesh();
  
  const vector<Element2D*>& elem = bspline->supportedElements();
  size_t ki;
  for (size_t ki=0; ki<elem.size(); ++ki)
    {
      // Localize element with regard to the information containers
      double umin = elem[ki]->umin();
      double umax = elem[ki]->umax();
      double vmin = elem[ki]->vmin();
      double vmax = elem[ki]->vmax();

      size_t kj1, kj2;
      for (kj1=1; kj1<kvec_u.size(); ++kj1)
	if (mesh->kval(XFIXED, kvec_u[kj1-1]) <= umin && 
	    mesh->kval(XFIXED, kvec_u[kj1]) >= umax)
	  break;
      for (kj2=1; kj2<kvec_v.size(); ++kj2)
	if (mesh->kval(YFIXED, kvec_v[kj2-1]) <= vmin && 
	    mesh->kval(YFIXED, kvec_v[kj2]) >= vmax)
	  break;

      if (elem[ki]->getNmbOutsideTol() > 0)
	{
	  u_info[kj1-1] += (umax-umin)*(vmax-vmin)*elem[ki]->getAccumulatedError();
	  v_info[kj2-1] += (umax-umin)*(vmax-vmin)*elem[ki]->getAccumulatedError();
	}

      // Element size
      u_elsize[kj2-1] += (umax-umin);
      v_elsize[kj1-1] += (vmax-vmin);
    } 

  // Modify priority information of strips to reduce the weight towards
  // the ends of the b-spline
  double fac1 = 0.25;
  double fac2 = 0.5;
  if (size1 >= 3)
    {
      u_info[0] *= fac1;
      u_info[size1-1] *= fac1;
    }
  if (size1 >= 4)
    {
      u_info[1] *= fac2;
      u_info[size1-2] *= fac2;
    }
  if (size2 >= 3)
    {
      v_info[0] *= fac1;
      v_info[size2-1] *= fac1;
    }
  if (size2 >= 4)
    {
      v_info[1] *= fac2;
      v_info[size2-2] *= fac2;
    }
    
  // Set treshhold for which strips to split
  double max_info = 0.0;
  double av_info = 0.0;
  int kj;
  for (kj=0; kj<size1; ++kj)
    {
      max_info = std::max(max_info, u_info[kj]);
      av_info += u_info[kj];
      v_elsize[kj] /= (double)size2;
    }
  for (kj=0; kj<size2; ++kj)
    {
      max_info = std::max(max_info, v_info[kj]);
      av_info += v_info[kj];
      u_elsize[kj] /= (double)size1;
    }
  av_info /= (double)(size1+size2);

  double threshhold = std::min(av_info, 0.5*max_info);
  double sizefac = 3.0;
  for (kj=0; kj<size1; ++kj)
    {
      double u1 = mesh->kval(XFIXED, kvec_u[kj]);
      double u2 = mesh->kval(XFIXED, kvec_u[kj+1]);
      if ((u_info[kj] >= threshhold || u2-u1 > sizefac*v_elsize[kj]) &&
	  (usize_min_ < 0.0 || (u2 - u1) >= 2.0*usize_min_))
	{
	  LRSplineSurface::Refinement2D curr_ref;
	  curr_ref.setVal(0.5*(u1+u2), bspline->vmin(), bspline->vmax(), XFIXED, 1);

	  // Check if the current refinement can be combined with an existing one
	  for (ki=0; ki<refs.size(); ++ki)
	    {
	      // Check direction and knot value
	      if (refs[ki].d == curr_ref.d && fabs(refs[ki].kval-curr_ref.kval) < tol)
		{
		  // Check extent of refinement
		  if (!(refs[ki].start > curr_ref.end+tol ||
			curr_ref.start > refs[ki].end+tol))
		    {
		      // Merge new knots
		      refs[ki].start = std::min(refs[ki].start, curr_ref.start);
		      refs[ki].end = std::max(refs[ki].end, curr_ref.end);
		      break;
		    }
		}
	    }

	  if (ki == refs.size())
	    refs.push_back(curr_ref);
	}
    }

  for (kj=0; kj<size2; ++kj)
    {
      double v1 = mesh->kval(YFIXED, kvec_v[kj]);
      double v2 = mesh->kval(YFIXED, kvec_v[kj+1]);
      if ((v_info[kj] >= threshhold  || v2-v1 > sizefac*u_elsize[kj])&&
	  (vsize_min_ < 0.0 || (v2 - v1) >= 2.0*vsize_min_))
	{
	  LRSplineSurface::Refinement2D curr_ref;
	  curr_ref.setVal(0.5*(v1+v2), bspline->umin(), bspline->umax(), YFIXED, 1);

	  // Check if the current refinement can be combined with an existing one
	  for (ki=0; ki<refs.size(); ++ki)
	    {
	      // Check direction and knot value
	      if (refs[ki].d == curr_ref.d && fabs(refs[ki].kval-curr_ref.kval) < tol)
		{
		  // Check extent of refinement
		  if (!(refs[ki].start > curr_ref.end+tol ||
			curr_ref.start > refs[ki].end+tol))
		    {
		      // Merge new knots
		      refs[ki].start = std::min(refs[ki].start, curr_ref.start);
		      refs[ki].end = std::max(refs[ki].end, curr_ref.end);
		      break;
		    }
		}
	    }

	  if (ki == refs.size())
	    refs.push_back(curr_ref);
	}
    }

}

//==============================================================================
void LRSurfApprox::checkFeasibleRef(Element2D* elem, 
				    vector<LRSplineSurface::Refinement2D>& refs,
				    vector<Element2D*>& affected)
//==============================================================================
{
  // Fetch B-splines
  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();

  int degree1 = srf_->degree(XFIXED);
  int degree2 = srf_->degree(YFIXED);
  int xmult = (degree1 <= 3) ? 1 : 2;
  int ymult = (degree2 <= 3) ? 1 : 2;
  
 // Refine one B-spline with support in the parent element in one or two
  // parameter directions depending on how many elements with a large error
  // that lies in its support
  // First check refinement in the u-direction
  double u_par = 0.5*(elem->umin() + elem->umax());
  int ixu = -1;
  double max_wgt = 0.0;
  size_t ki, kj;
  for (ki=0; ki<nmb; ++ki)
    {
      // Count the number of elements with large error affected
      double curr_wgt = 0.0;
      const vector<Element2D*>& curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->umax() < u_par || curr_el[kj]->umin() > u_par)
	    continue;  // Element not affected

	  // Compute weight for importance of refinement
	  double max_err, av_err;
	  int nmb_outside;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el[kj]->nmbDataPoints();
	  double wgt = av_err*(int)nmb_outside/(int)nmb_pts;
	  curr_wgt += wgt;
	}
      if (curr_wgt > max_wgt)
	{
	  max_wgt = curr_wgt;
	  ixu = (int)ki;
	}
    }
	  
  // The v-direction
  double v_par = 0.5*(elem->vmin() + elem->vmax());
  int ixv = -1;
  max_wgt = 0.0;
  for (ki=0; ki<nmb; ++ki)
    {
      // Count the number of elements with large error affected
      double curr_wgt = 0.0;
      const vector<Element2D*>& curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->vmax() < v_par || curr_el[kj]->vmin() > v_par)
	    continue;  // Element not affected

	  // Compute weight for importance of refinement
	  double max_err, av_err;
	  int nmb_outside;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el[kj]->nmbDataPoints();
	  double wgt = av_err*(int)nmb_outside/(int)nmb_pts;
	  curr_wgt += wgt;
	}
      if (curr_wgt > max_wgt)
	{
	  max_wgt = curr_wgt;
	  ixv = (int)ki;
	}
    }
  
  // Fetch the affected elements in both parameter directions
  double fac = 0.1;
  double fac3 = 0.95;
  vector<Element2D*> aff_u;
  int nmb_u = 0;
  if (ixu >= 0)
    {
      const vector<Element2D*>& curr_el_u = bsplines[ixu]->supportedElements();
      for (kj=0; kj<curr_el_u.size(); ++kj)
	{
	  if (true /*curr_el_u[kj]->umax() > u_par && curr_el_u[kj]->umin() < u_par*/)
	    {
	      double max_err, av_err;
	      int nmb_outside;
	      curr_el_u[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	      int nmb_pts = curr_el_u[kj]->nmbDataPoints();
	      aff_u.push_back(curr_el_u[kj]);
	      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
		nmb_u++;
	    }
	}
    }

  vector<Element2D*> aff_v;
  int nmb_v = 0;
  if (ixv >= 0)
    {
      const vector<Element2D*>& curr_el_v = bsplines[ixv]->supportedElements();
      for (kj=0; kj<curr_el_v.size(); ++kj)
	{
	  if (true /*curr_el_v[kj]->vmax() > v_par && curr_el_v[kj]->vmin() < v_par*/)
	    {
	      double max_err, av_err;
	      int nmb_outside;
	      curr_el_v[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	      int nmb_pts = curr_el_v[kj]->nmbDataPoints();
	      aff_v.push_back(curr_el_v[kj]);
	      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
		nmb_v++;
	    }
	}
    }

  // Assemble information
  double fac2 = 0.5;
  std::set<Element2D*> affected_combined;
  if (ixu >= 0 && nmb_u > (int)fac2*bsplines[ixu]->degree(XFIXED))
    {
      affected_combined.insert(aff_u.begin(), aff_u.end());
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(u_par, bsplines[ixu]->vmin(), bsplines[ixu]->vmax(),
		      XFIXED, xmult);
      refs.push_back(curr_ref);
    }
			       
    if (ixv >= 0 && nmb_v > (int)fac2*bsplines[ixv]->degree(YFIXED))
    {
      affected_combined.insert(aff_v.begin(), aff_v.end());
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(v_par, bsplines[ixv]->umin(), bsplines[ixv]->umax(),
		      YFIXED, ymult);
      refs.push_back(curr_ref);
    }

    affected.insert(affected.end(), affected_combined.begin(), affected_combined.end());
}

int comp_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int comp_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

//==============================================================================
void LRSurfApprox::constructGhostPoints(vector<double>& ghost_points)
//==============================================================================
{
  // Sort the points in a grid
  //int nmb_u=10, nmb_v=10;
  int nmb_u=30, nmb_v=30;
  vector<vector<int> > end_ix;
  // First sort the points in the 1. parameter directions
  int dim = srf_->dimension();
  int del = dim+2;  // Parameter pair and position
  int nmb = (int)points_.size()/del;  // Number of data points
  double facsquare = 0.2;

  nmb_u = std::max(5, std::min(nmb_u, (int)(facsquare*sqrt((double)nmb))));
  nmb_v = std::max(5, std::min(nmb_v, (int)(facsquare*sqrt((double)nmb))));

  int ki, kj, kr, kh;
  vector<double> minmax(2*dim);
  for (ki=0; ki<dim; ++ki)
    minmax[2*ki] = minmax[2*ki+1] = points_[del-dim+ki];
  for (kj=0; kj<nmb; ++kj)
    for (ki=0; ki<dim; ++ki)
      {
	minmax[2*ki] = std::min(minmax[2*ki], points_[(kj+1)*del-dim+ki]);
	minmax[2*ki+1] = std::max(minmax[2*ki+1], points_[(kj+1)*del-dim+ki]);
      }
  for (ki=0; ki<dim; ++ki)
    {
      double delta = minmax[2*ki+1] - minmax[2*ki];
      minmax[2*ki] -= 0.1*delta;
      minmax[2*ki+1] += 0.1*delta;
    }

  // // TESTING
  // for (ki=0; ki<dim; ++ki)
  //   minmax[2*ki] = 0.0;
  

  qsort(&points_[0], nmb, del*sizeof(double), comp_u_par);
  
  double u[2], v[2];
  u[0] = srf_->paramMin(XFIXED);
  u[1] = srf_->paramMax(XFIXED);
  v[0] = srf_->paramMin(YFIXED);
  v[1] = srf_->paramMax(YFIXED);
  double u_del = (u[1] - u[0])/(double)(nmb_u-1);
  double v_del = (v[1] - v[0])/(double)(nmb_v-1);
  //double fac2 = 0.01;
  //int min_nmb1 = std::max(std::min(100, (int)(0.25*nmb)), (int)(fac2*nmb));

  // Sort each u-strip in the v-direction and create grid
  int pp0, pp1, pp2, pp3;
  int ix0, ix1, ix2, ix3;
  double upar, vpar;
  for (ki=0, pp0=0, upar=u[0]+u_del, ix0=0; ki<nmb_u-1; 
       ++ki, upar+=u_del, ix0=ix1, pp0=pp1)
    {
      vector<int> ixs;
      for (pp1=pp0, ix1=ix0; pp1<(int)points_.size() && points_[pp1]<upar; 
	   pp1+=del, ++ix1);
      if (ki == nmb_u-2)
	{
	  pp1 = (int)points_.size();
	  ix1 = pp1/del;
	}

      // Sort according to the v-parameter
      qsort(&points_[pp0], ix1-ix0, del*sizeof(double), comp_v_par);
	   
      ixs.push_back(ix0);
      for (kj=0, pp2=pp0, vpar=v[0]+v_del, ix2=ix0; kj<nmb_v-1;
	   ++kj, vpar+=v_del, ix2=ix3, pp2=pp3)
	{
	  for (pp3=pp2, ix3=ix2; pp3<pp1 && points_[pp3+1]<vpar; pp3+=del, ++ix3);
	  if (kj == nmb_v-2)
	    {
	      pp3 = pp1;
	      ix3 = pp3/del;
	    }
	  ixs.push_back(ix3);
	}
      end_ix.push_back(ixs);
    }

  // Compute the average number of points in each cell (with points)
  int nmb_cell=0;
  int av_pts = 0;
  int prev = 0;
  for (ki=0; ki<(int)end_ix.size(); ++ki)
    for (kj=1, prev=end_ix[ki][0]; kj<(int)end_ix[ki].size(); ++kj)
      {
	if (end_ix[ki][kj] > prev)
	  {
	    av_pts += (end_ix[ki][kj]-prev);
	    nmb_cell++;
	  }
	prev = end_ix[ki][kj];
      }
  av_pts /= nmb_cell;
    
#ifdef DEBUG
  std::ofstream of1("ghost_pnt1.g2");
#endif

  // For all u-strips, traverse cells in the v-direction and construct ghost points 
  // in cells with non or few points
  int min_nmb = std::max(3, std::min(av_pts/200, 100));
  int prev2;
  int nmb_ev = std::min(std::max(av_pts/300, 5),10);
  vector<double> ghost1;
  for (ki=0, prev=0; ki<nmb_u-1; ++ki)
    {
      for (kj=1; kj<nmb_v; kj=kr+1)
	{
	  // Find the last cell with an acceptable number of points
	  for (kh=kj-1, prev2=prev; kh<nmb_v-1; prev2=end_ix[ki][kh+1], ++kh)
	    if (end_ix[ki][kh+1] - prev2 < min_nmb)
	      break;
	  prev = prev2;
	  kh = std::max(0, kh-1);

	  // Find the first cell with an acceptable number of points
	  for (kr=kh+1; kr<nmb_v-1; prev2=end_ix[ki][kr+1], ++kr)
	    if (end_ix[ki][kr+1] - prev2 >= min_nmb)
	      break;
	  kr = std::min(nmb_v-1, kr+1);

	  if (kr > kh+1 && end_ix[ki][kr]-end_ix[ki][kh] > 3)
	    {
	      // Construct ghost points
	      constructLocalGhostPts(&points_[del*end_ix[ki][kh]], 
				     end_ix[ki][kr]-end_ix[ki][kh],
				     del, u[0]+ki*u_del, u[0]+(ki+1)*u_del,
				     v[0]+kh*v_del, v[0]+kr*v_del,
				     2*nmb_ev, 2*(kr-kh-1)*nmb_ev, ghost1);
#ifdef DEBUG
	      int nmb_ghost1 = (int)ghost1.size()/del;
	      of1 << "400 1 0 4 0 255 0 255" << std::endl;
	      of1 << nmb_ghost1 << std::endl;
	      
	      for (int ki2=0; ki2<nmb_ghost1; ++ki2)
		{
		  for (int kj2=0; kj2<del; ++kj2)
		    of1 << ghost1[ki2*del+kj2] << " ";
		  of1 << std::endl;
		}
	      int stop_break1 = 1;
#endif
	    }
	}
    }

  // For all v-strips, traverse cells in the v-direction and construct ghost points 
  // in cells with non or few points
  vector<double> ghost2;
  for (ki=1; ki<nmb_v; ++ki)
    {
      for (kj=0; kj<nmb_u-1; kj=kr)
	{
	  // Find the last cell with an acceptable number of points
	  for (kh=kj; kh<nmb_u-1; ++kh)
	    if (end_ix[kh][ki] - end_ix[kh][ki-1] < min_nmb)
	      break;
	  kh = std::max(0, kh-1);

	  // Find the first cell with an acceptable number of points
	  for (kr=kh+1; kr<nmb_u-1; ++kr)
	    if (end_ix[kr][ki] - end_ix[kr][ki-1] >= min_nmb)
	      break;
	  kr = std::min(nmb_u-1, kr);

	  if (kr > kh+1)
	    {
	      // Collect input points
	      vector<double> points2;
	      for (int kk=kh; kk<=std::min(kr,nmb_u-2); ++kk)
		points2.insert(points2.end(), points_.begin()+del*end_ix[kk][ki-1], 
			       points_.begin()+del*end_ix[kk][ki]);
	      
	      // Construct ghost points
	      int npt = (int)points2.size()/del;
	      if (npt > 3)
		{
		  constructLocalGhostPts(&points2[0], npt, del,
					 u[0]+kh*u_del, u[0]+kr*u_del,
					 v[0]+(ki-1)*v_del, v[0]+ki*v_del,
					 2*(kr-kh-1)*nmb_ev, 2*nmb_ev, ghost2);
#ifdef DEBUG
		  int nmb_ghost2 = (int)ghost2.size()/del;
		  of1 << "400 1 0 4 255 0 0 255" << std::endl;
		  of1 << nmb_ghost2 << std::endl;
	      
		  for (int ki2=0; ki2<nmb_ghost2; ++ki2)
		    {
		      for (int kj2=0; kj2<del; ++kj2)
			of1 << ghost2[ki2*del+kj2] << " ";
		      of1 << std::endl;
		    }
		  int stop_break2 = 1;
#endif
		}
	    }
	}
    }
  
  // Distribute the intermediate ghost points to the cells
  double u1, u2, v1, v2;
  vector<vector<vector<double> > > tmp_pts;
  tmp_pts.resize(nmb_u-1);
  qsort(&ghost1[0], (int)ghost1.size()/del, del*sizeof(double), comp_u_par);
  qsort(&ghost2[0], (int)ghost2.size()/del, del*sizeof(double), comp_u_par);
  int pg1_0, pg1_1, pg1_2, pg1_3, pg1_4, pg2_0, pg2_1, pg2_2, pg2_3, pg2_4;
  pg1_1 = pg2_1 = 0;
  for (ki=0, pg1_0=0, pg2_0=0, upar=u[0]+u_del; ki<nmb_u-1; 
       ++ki, upar+=u_del, pg1_0=pg1_1, pg2_0=pg2_1)
    {
      tmp_pts[ki].resize(nmb_v-1);

      // Make extended point set of ghost points
      u1 = std::max(u[0], upar - 1.5*u_del);
      u2 = std::min(u[1], upar + 0.5*u_del);

      for (pg1_1=pg1_0; pg1_1<(int)ghost1.size() && ghost1[pg1_1]<u1; pg1_1+=del);
      for (pg1_2=pg1_1; pg1_2<(int)ghost1.size() && ghost1[pg1_2]<u2; pg1_2+=del);
      for (pg2_1=pg2_0; pg2_1<(int)ghost2.size() && ghost2[pg2_1]<u1; pg2_1+=del);
      for (pg2_2=pg2_1; pg2_2<(int)ghost2.size() && ghost2[pg2_2]<u2; pg2_2+=del);

      // Copy relevant selection of intermediate ghost points
      vector<double> ghost1_v(ghost1.begin()+pg1_1, ghost1.begin()+pg1_2);
      vector<double> ghost2_v(ghost2.begin()+pg2_1, ghost2.begin()+pg2_2);

      // Sort according to the v-parameter
      qsort(&ghost1_v[0], (int)ghost1_v.size()/del, del*sizeof(double), comp_v_par);
      qsort(&ghost2_v[0], (int)ghost2_v.size()/del, del*sizeof(double), comp_v_par);
	   
      pg1_3 = pg2_3 = 0;
      for (kj=0, pg1_2=0, pg2_2=0, vpar=v[0]+v_del; 
	   kj<nmb_v-1; ++kj, vpar+=v_del, pg1_2=pg1_3, pg2_2=pg2_3)
	{
	  if (end_ix[ki][kj+1]-end_ix[ki][kj] > min_nmb)
	    continue;

	  // Make extended point set of ghost points
	  v1 = std::max(v[0], vpar - 1.5*v_del);
	  v2 = std::min(v[1], vpar + 0.5*v_del);

	  for (pg1_3=pg1_2; pg1_3<(int)ghost1_v.size() && ghost1_v[pg1_3+1]<v1; pg1_3+=del);
	  for (pg1_4=pg1_3; pg1_4<(int)ghost1_v.size() && ghost1_v[pg1_4+1]<v2; pg1_4+=del);
	  for (pg2_3=pg2_2; pg2_3<(int)ghost2_v.size() && ghost2_v[pg2_3+1]<v1; pg2_3+=del);
	  for (pg2_4=pg2_3; pg2_4<(int)ghost2_v.size() && ghost2_v[pg2_4+1]<v2; pg2_4+=del);

	  if (pg1_4 > pg1_3)
	    tmp_pts[ki][kj].insert(tmp_pts[ki][kj].end(), ghost1_v.begin()+pg1_3,
				   ghost1_v.begin()+pg1_4);
	  if (pg2_4 > pg2_3)
	    tmp_pts[ki][kj].insert(tmp_pts[ki][kj].end(), ghost2_v.begin()+pg2_3,
				   ghost2_v.begin()+pg2_4);
	  if (tmp_pts[ki][kj].size() > 0)
	    tmp_pts[ki][kj].insert(tmp_pts[ki][kj].end(), 
				   points_.begin()+del*end_ix[ki][kj],
				   points_.begin()+del*end_ix[ki][kj+1]);
	}
    }

#ifdef DEBUG
  std::ofstream of3("approx_sf.g2");
#endif

  // For each cell with ghost points, create a surface approximating ghost points
  // and initial data points. Construct new ghost points by evaluating this surface
      
  int nmb_evu = nmb_ev;
  int nmb_evv = nmb_ev;
  if (grid_)
    {
      nmb_evu = std::max(1, std::min(nmb_evu, (int)(u_del/cell_size_[0])));
      nmb_evv = std::max(1, std::min(nmb_evv, (int)(v_del/cell_size_[1])));
    }

  double delu = u_del/(double)nmb_evu;
  double delv = v_del/(double)nmb_evu;
  for (ki=0; ki<nmb_u-1; ++ki)
    for (kj=0; kj<nmb_v-1; ++kj)
      {
	if (tmp_pts[ki][kj].size() == 0)
	  continue;

	// Make approximative surface
	double smoothweight = 0.2;  
	int order = 3;
	vector<double> knots_u(2*order), knots_v(2*order);
	u1 = std::max(u[0], u[0] + ki*u_del - 0.5*u_del);
	u2 = std::min(u[1], u[0] + (ki+1)*u_del + 0.5*u_del);
	v1 = std::max(v[0], v[0] + kj*v_del - 0.5*u_del);
	v2 = std::min(v[1], v[0] + (kj+1)*v_del + 0.5*v_del);
	for (int ix=0; ix<order; ++ix)
	  {
	    knots_u[ix] = u1;
	    knots_u[order+ix] = u2;
	    knots_v[ix] = v1;
	    knots_v[order+ix] = v2;
	  }

	  double maxdist, avdist;
	  int outsideeps;
	  int kn2 = (int)tmp_pts[ki][kj].size()/del;
	  shared_ptr<SplineSurface> surf = 
	    createSurf(&tmp_pts[ki][kj][0], kn2, dim, order, order,
		       order, order, &knots_u[0], &knots_v[0], smoothweight,
		       maxdist, avdist, outsideeps);

#ifdef DEBUG
	  std::ofstream of4("curr_approx_sf.g2");
	  of4 << "400 1 0 4 100 155 0 255 " << std::endl;
	  of4 << kn2 << std::endl;
	  for (size_t ki2=0; ki2<tmp_pts[ki][kj].size(); ki2+=del)
	    of4 << tmp_pts[ki][kj][ki2] << " " << tmp_pts[ki][kj][ki2+1] << " " << tmp_pts[ki][kj][ki2+2] << std::endl;
	  
    shared_ptr<LRSplineSurface> tmp_srf(new LRSplineSurface(surf.get(), 1.0e-6));
    if (tmp_srf->dimension() == 1)
      tmp_srf->to3D();
    tmp_srf->writeStandardHeader(of3);
    tmp_srf->write(of3);

    tmp_srf->writeStandardHeader(of4);
    tmp_srf->write(of4);
#endif

	  // Evaluate
	  int nu = nmb_evu + (ki == nmb_u-1);
	  int nv = nmb_evv + (kj == nmb_v-1);
	  for (kr=0, u1=u[0]+ki*u_del; kr<nu; ++kr, u1+=delu)
	    for (kh=0, v1=v[0]+kj*v_del; kh<nv; ++kh, v1+=delv)
	      {
		Point pos = surf->ParamSurface::point(u1, v1);

		// Limit ghost point to stay within allowed interval
		for (int ix=0; ix<dim; ++ix)
		  pos[ix] = std::max(minmax[2*ix], std::min(minmax[2*ix+1],pos[ix]));

		ghost_points.push_back(u1);
		ghost_points.push_back(v1);
		ghost_points.insert(ghost_points.end(), pos.begin(), pos.end());
		//ghost_points.push_back(0.0);  // TEST for Liguria data
	      }
      }


#ifdef DEBUG
  int nmb_ghost = (int)ghost_points.size()/del;
  std::ofstream of("ghost_pnt.g2");
  of << "400 1 0 4 55 100 100 255" << std::endl;
  of << nmb_ghost << std::endl;

  for (ki=0; ki<nmb_ghost; ++ki)
    {
      for (kj=0; kj<del; ++kj)
	of << ghost_points[ki*del+kj] << " ";
      of << std::endl;
    }
#endif

  
}


//==============================================================================
void LRSurfApprox::constructLocalGhostPts(double *startpt, int kn2,
					  int dim, double u1, double u2,
					  double v1, double v2, 
					  int nmb_u, int nmb_v,
					  vector<double>& ghostpts)
//==============================================================================
  {
#ifdef DEBUG
  std::ofstream of2("corner_cloud.g2");
  std::ofstream of3("corner_sf.g2");
  of2 << "400 1 0 4 100 100 55 255 " << std::endl;
  of2 << kn2 << std::endl;
#endif
    // Bound pointset in all directions
    vector<double> ptbound(2*dim);
    int ix, ix2;
    for (ix2=0; ix2<dim; ++ix2)
      ptbound[2*ix2] = ptbound[2*ix2+1] = startpt[ix2];
   double *curr = startpt;
   for (ix=0; ix<kn2; ++ix, curr+=dim)
     {
      for (ix2=0; ix2<dim; ++ix2)
	{
	  ptbound[2*ix2] = std::min(ptbound[2*ix2], curr[ix2]); 
	  ptbound[2*ix2+1] = std::max(ptbound[2*ix2+1], curr[ix2]); 
#ifdef DEBUG
	  of2 << curr[ix2] << " ";
#endif
	}
#ifdef DEBUG
      of2 << std::endl;
#endif
     }

   double min_size = 0.001;
   if (ptbound[1]-ptbound[0] < min_size || 
       ptbound[3]-ptbound[2] < min_size)
     return; // No method to construct ghost points

   vector<double> points;
   double *currpt;
   if (dim == 3)
     currpt = startpt;
   else
     {
       for (ix=0; ix<kn2; ++ix)
	 points.insert(points.end(), startpt+ix*dim+2, startpt+(ix+1)*dim);
       currpt = &points[0];
     }

    // Make approximative constant surface
    double smoothweight = 0.2;  // Prioritize smoothing
    int order = 1;
    vector<double> knots_u(2*order), knots_v(2*order);
    for (ix=0; ix<order; ++ix)
      {
	knots_u[ix] = ptbound[0];
	knots_u[order+ix] = ptbound[1];
	knots_v[ix] = ptbound[2];
	knots_v[order+ix] = ptbound[3];
      }

    double maxdist, avdist;
    int outsideeps;
    shared_ptr<SplineSurface> surf0;
    if (dim == 3)
      {
	surf0 = createSurf(currpt, kn2, dim-2, order, order,
			   order, order, &knots_u[0], &knots_v[0], smoothweight,
			   maxdist, avdist, outsideeps);
      }
// #ifdef DEBUG
//     shared_ptr<LRSplineSurface> tmp_srf0(new LRSplineSurface(surf0.get(), 1.0e-6));
//     tmp_srf0->to3D();
//     tmp_srf0->writeStandardHeader(of3);
//     tmp_srf0->write(of3);
// #endif

    // Make approximative linear surface
    order = 2;
    knots_u.insert(knots_u.begin(), ptbound[0]);
    knots_u.push_back(ptbound[1]);
    knots_v.insert(knots_v.begin(), ptbound[2]);
    knots_v.push_back(ptbound[3]);

    shared_ptr<SplineSurface> surf1 = 
      createSurf(currpt, kn2, dim-2, order, order,
		 order, order, &knots_u[0], &knots_v[0], smoothweight,
		 maxdist, avdist, outsideeps);
#ifdef DEBUG
    shared_ptr<LRSplineSurface> tmp_srf1(new LRSplineSurface(surf1.get(), 1.0e-6));
    if (tmp_srf1->dimension() == 1)
      tmp_srf1->to3D();
    tmp_srf1->writeStandardHeader(of3);
    tmp_srf1->write(of3);
#endif

    // Make approximative quadratic surface
    order = 3;
    knots_u.insert(knots_u.begin(), ptbound[0]);
    knots_u.push_back(ptbound[1]);
    knots_v.insert(knots_v.begin(), ptbound[2]);
    knots_v.push_back(ptbound[3]);
    shared_ptr<SplineSurface> surf2 = 
      createSurf(currpt, kn2, dim-2, order, order,
		 order, order, &knots_u[0], &knots_v[0], smoothweight,
		 maxdist, avdist, outsideeps);

#ifdef DEBUG
    shared_ptr<LRSplineSurface> tmp_srf2(new LRSplineSurface(surf2.get(), 1.0e-6));
    if (tmp_srf2->dimension() == 1)
      tmp_srf2->to3D();
    tmp_srf2->writeStandardHeader(of3);
    tmp_srf2->write(of3);
#endif

    // Evaluate ghost points
    int ki, kj;
    double upar, vpar;
    double udel = (u2 - u1)/(double)(nmb_u-1);
    double vdel = (v2 - v1)/(double)(nmb_v-1);
    for (ki=0, upar=u1; ki<nmb_u; upar+=udel, ++ki)
      for (kj=0, vpar=v1; kj<nmb_v; vpar+=vdel, ++kj)
	{
	  Point pos;
	  if (dim == 3)
	    {
	      Point pos0 = surf0->ParamSurface::point(upar, vpar);
	      Point pos1 = surf1->ParamSurface::point(upar, vpar);
	      Point pos2 = surf2->ParamSurface::point(upar, vpar);
	      pos = (pos0+pos1+pos2)/3.0;
	    }
	  else if (initMBA_)
	    {
	      Point pos0 = srf_->ParamSurface::point(upar, vpar);
	      Point pos1 = surf1->ParamSurface::point(upar, vpar);
	      Point pos2 = surf2->ParamSurface::point(upar, vpar);
	      pos = (pos0+pos1+pos2)/3.0;
	    }
	  else
	    {
	      Point pos1 = surf1->ParamSurface::point(upar, vpar);
	      Point pos2 = surf2->ParamSurface::point(upar, vpar);
	      pos = 0.5*(pos1+pos2);
	    }
	  for (ix=2; ix<dim; ++ix)
	    {
	      double tmp = ptbound[2*ix+1] - ptbound[2*ix];
	      pos[ix-2] = std::max(ptbound[2*ix]-0.5*tmp,
				   std::min(pos[ix-2], ptbound[2*ix+1]+0.5*tmp));
	      ghostpts.push_back(upar);
	      ghostpts.push_back(vpar);
	      ghostpts.insert(ghostpts.end(), pos.begin(), pos.end());
	    }
	}
  }


// //==============================================================================
// void LRSurfApprox::constructGhostPoints(vector<double>& ghost_points)
// //==============================================================================
// {
//   // 1. version. Only corner points are constructed.
//   // For each corner, fetch points in the neighbourhood of the corner
 
//   // First sort the points in the 1. parameter directions
//   int dim = srf_->dimension();
//   int del = dim+2;                   // Number of entries for each point
//   int nmb = (int)points_.size()/del;  // Number of data points

//   qsort(&points_[0], nmb, del*sizeof(double), comp_u_par);
  
//   double u[2], v[2];
//   u[0] = srf_->paramMin(XFIXED);
//   u[1] = srf_->paramMax(XFIXED);
//   v[0] = srf_->paramMin(YFIXED);
//   v[1] = srf_->paramMax(YFIXED);
//   double fac1 = 0.1;
//   double u_del = fac1*(u[1] - u[0]);
//   double v_del = fac1*(v[1] - v[0]);
//   double fac2 = 0.01;
//   int min_nmb1 = std::max(std::min(100, (int)(0.25*nmb)), (int)(fac2*nmb));

//   int u_nmb = (int)nmb/(v[1]-v[0]);
//   int v_nmb = (int)nmb/(u[1]-u[0]);

// #ifdef DEBUG
//   std::ofstream of2("corner_cloud.g2");
//   std::ofstream of3("corner_sf.g2");
// #endif

//   // Fetch point strips close to the boundaries in the u direction
//   int sgn1 = 1;
//   int pp1, pp2, ki, kj, kn1, kn2;
//   int start_ix;
//   double upar;
//   for (ki=0, pp1=0, upar=u[0]+u_del, start_ix=0; ki<2; 
//        ++ki, pp1=(int)points_.size()-del, sgn1=-1, upar=u[1]-u_del)
//     {
//       for (kn1=0; sgn1*points_[pp1]<sgn1*upar || kn1<min_nmb1; 
// 	   ++kn1, pp1+=sgn1*del);
//       if (ki == 1)
// 	{
// 	  start_ix = pp1 + sgn1*del;
// 	  upar = std::min(upar, points_[start_ix]);
// 	}
//       else
// 	upar = std::max(upar, points_[pp1 - sgn1*del]);

//       // Sort the current sub set of points according to the v-parameter
//       qsort(&points_[start_ix], kn1, del*sizeof(double), comp_v_par);

//       // Fetch point strips close to the boundaries in the v direction
//       int sgn2 = 1;
//       double vpar;
//       int min_nmb2 = std::max(10, (int)(fac2*kn1));
//       for (kj=0, pp2=start_ix, vpar=v[0]+v_del; kj<2; 
// 	   ++kj, pp2=start_ix+kn1*del, sgn2=-1, vpar=v[1]-v_del)
// 	{
// 	  for (kn2=0; sgn2*points_[pp2+1]<sgn2*vpar || kn2<min_nmb2; ++kn2, 
// 		 pp2+=sgn2*del);
// 	  if (kj == 1)
// 	    {
// 	      start_ix = pp2 + sgn2*del;
// 	      vpar = std::min(vpar, points_[start_ix+1]);
// 	    }
// 	  else
// 	    vpar = std::max(vpar, points_[pp2 - sgn2*del]);

// #ifdef DEBUG
// 	  of2 << "400 1 0 4 100 100 55 255" << std::endl;
// 	  of2 << kn2 << std::endl;
// 	  for (int kr=0; kr<kn2; ++kr)
// 	    {
// 	      for (int kh=0; kh<3; ++kh)
// 		of2 << points_[start_ix+kr*3+kh] << " ";
// 	      of2 << std::endl;
// 	    }
// #endif

// 	  // Bound pointset in all directions
// 	  vector<double> ptbound(2*del);
// 	  int ix, ix2;
// 	  for (ix2=0; ix2<del; ++ix2)
// 	    ptbound[2*ix2] = ptbound[2*ix2+1] = points_[start_ix+ix2];
// 	  for (ix=0; ix<kn2; ++ix)
// 	    for (ix2=0; ix2<del; ++ix2)
// 	      {
// 		ptbound[2*ix2] = 
// 		  std::min(ptbound[2*ix2], points_[start_ix+ix*del+ix2]); 
// 		ptbound[2*ix2+1] = 
// 		  std::max(ptbound[2*ix2+1], points_[start_ix+ix*del+ix2]); 
// 	      }

// 	  // Make approximative surface
// 	  double smoothweight = 0.2;  // Prioritize smoothing
// 	  int order = 3;
// 	  vector<double> knots_u(2*order), knots_v(2*order);
// 	  for (ix=0; ix<order; ++ix)
// 	    {
// 	      // knots_u[ix] = (ki==0) ? u[0] : upar;
// 	      // knots_u[order+ix] = (ki==0) ? upar : u[1];
// 	      // knots_v[ix] = (kj==0) ? v[0] : vpar;
// 	      // knots_v[order+ix] = (kj==0) ? vpar : v[1];
// 	      knots_u[ix] = ptbound[0];
// 	      knots_u[order+ix] = ptbound[1];
// 	      knots_v[ix] = ptbound[2];
// 	      knots_v[order+ix] = ptbound[3];
// 	    }

// 	  double maxdist, avdist;
// 	  int outsideeps;
// 	  shared_ptr<SplineSurface> surf = 
// 	    createSurf(&points_[start_ix], kn2, dim, order, order,
// 		       order, order, &knots_u[0], &knots_v[0], smoothweight,
// 		       maxdist, avdist, outsideeps);

// 	  // Fetch corner point
// 	  Point corner = surf->ParamSurface::point(u[ki], v[kj]);
// 	  for (ix2=2; ix2<del; ++ix2)
// 	    {
// 	      double tmp = ptbound[2*ix2+1]-ptbound[2*ix2];
// 	      corner[ix2-2] = std::max(ptbound[2*ix2]-0.5*tmp, 
// 				       std::min(corner[ix2-2], ptbound[2*ix2+1]+tmp));
// 	    }

// 	  ghost_points.push_back(u[ki]);
// 	  ghost_points.push_back(v[kj]);
// 	  ghost_points.insert(ghost_points.end(), corner.begin(), corner.end());
// #ifdef DEBUG
// 	  if (surf.get())
// 	    {
// 	      shared_ptr<LRSplineSurface> tmp_srf(new LRSplineSurface(surf.get(), 1.0e-6));
// 	      tmp_srf->to3D();
// 	      tmp_srf->writeStandardHeader(of3);
// 	      tmp_srf->write(of3);
// 	    }
// #endif

// 	  bool add_bd_pts = true;
// 	  if (add_bd_pts)
// 	    {
// 	  // Additional points
// 	  double fac3 = 0.25;
// 	  double dom1[4];
// 	  double ulen = (ki==0) ? upar-u[0] : u[1]-upar;
// 	  if (ki == 0 && ptbound[0] > u[0] + 0.1*ulen)
// 	    {
// 	      dom1[0] = u[0];
// 	      dom1[1] = ptbound[0];
// 	    }
// 	  else if (ki == 1 && ptbound[1] < u[1] - 0.1*ulen)
// 	    {
// 	      dom1[0] = ptbound[1];
// 	      dom1[1] = upar;
// 	    }
// 	  else
// 	    dom1[0] = dom1[1] = upar;
// 	  dom1[2] = (kj==0) ? v[0] : vpar;
// 	  dom1[3] = (kj==1) ? v[1] : vpar;
	  
// 	  if (dom1[1] - dom1[0] > 0.1*ulen)
// 	    {
// 	      // Make additional points
// 	      int unmb = std::max(1, (int)(fac3*u_nmb*(dom1[1]-dom1[0])/(u[1]-u[0])));
// 	      //(int)(0.1*kn2*(dom1[1]-dom1[0])/ulen);
// 	      int vnmb = std::max(3, (int)(fac3*v_nmb*(dom1[3]-dom1[2])/(v[1]-v[0])));
// 				  //(int)(0.1*kn2));
// 	      double ustep = (dom1[1] - dom1[0])/(double)unmb;
// 	      double vstep = (dom1[3] - dom1[2])/(double)vnmb;
// 	      double u1, v1;
// 	      int kr1, kr2;
// 	      for (kr2=0, v1=dom1[2]+(kj==1)*vstep; kr2<vnmb; ++kr2, v1+=vstep)
// 		for (kr1=0, u1=dom1[0]+(ki==1)*ustep; kr1<unmb; ++kr1, u1+=ustep)
// 		  {
// 		    Point pos = surf->ParamSurface::point(u1, v1);
// 		    for (ix2=2; ix2<del; ++ix2)
// 		      {
// 			double tmp = ptbound[2*ix2+1]-ptbound[2*ix2];
// 			pos[ix2-2] = std::max(ptbound[2*ix2]-0.5*tmp, 
// 					      std::min(pos[ix2-2], ptbound[2*ix2+1]+tmp));
// 		      }
		    
// 		    ghost_points.push_back(u1);
// 		    ghost_points.push_back(v1);
// 		    ghost_points.insert(ghost_points.end(), pos.begin(), pos.end());
// 		  }
// 	    }

// 	  double dom2[4];
// 	  double vlen = (kj==0) ? vpar-v[0] : v[1]-vpar;
// 	  if (kj == 0 && ptbound[2] > v[0] + 0.1*vlen)
// 	    {
// 	      dom2[2] = v[0];
// 	      dom2[3] = ptbound[2];
// 	    }
// 	  else if (kj == 1 && ptbound[3] < v[1] - 0.1*vlen)
// 	    {
// 	      dom2[2] = ptbound[3];
// 	      dom2[3] = vpar;
// 	    }
// 	  else
// 	    dom2[2] = dom2[3] = vpar;
// 	  dom2[0] = (ki==0) ? u[0] : upar;
// 	  dom2[1] = (ki==1) ? u[1] : upar;
	  
// 	  if (dom2[3] - dom2[2] > 0.1*vlen)
// 	    {
// 	      // Make additional points
// 	      int unmb = std::max(3, (int)(fac3*u_nmb*(dom1[1]-dom1[0])/(u[1]-u[0])));
// 	      //(int)(0.1*kn2*(dom1[1]-dom1[0])/ulen);
// 	      int vnmb = std::max(1, (int)(fac3*v_nmb*(dom1[3]-dom1[2])/(v[1]-v[0])));
// 	      double ustep = (dom2[1] - dom2[0])/(double)unmb;
// 	      double vstep = (dom2[3] - dom2[2])/(double)vnmb;
// 	      double u1, v1;
// 	      int kr1, kr2;
// 	      for (kr2=0, v1=dom2[2]+(kj==1)*vstep; kr2<vnmb; ++kr2, v1+=vstep)
// 		for (kr1=0, u1=dom2[0]+(ki==1)*ustep; kr1<unmb; ++kr1, u1+=ustep)
// 		  {
// 		    Point pos = surf->ParamSurface::point(u1, v1);
// 		    for (ix2=2; ix2<del; ++ix2)
// 		      {
// 			double tmp = ptbound[2*ix2+1]-ptbound[2*ix2];
// 			pos[ix2-2] = std::max(ptbound[2*ix2]-0.5*tmp, 
// 					      std::min(pos[ix2-2], ptbound[2*ix2+1]+tmp));
// 		      }
		    
// 		    ghost_points.push_back(u1);
// 		    ghost_points.push_back(v1);
// 		    ghost_points.insert(ghost_points.end(), pos.begin(), pos.end());
// 		  }
// 	    }
// 	    }
// 	}
//     }
// #ifdef DEBUG
//   int nmb_ghost = (int)ghost_points.size()/del;
//   std::ofstream of("ghost_pnt.g2");
//   of << "400 1 0 4 0 255 0 255" << std::endl;
//   of << nmb_ghost << std::endl;

//   for (ki=0; ki<nmb_ghost; ++ki)
//     {
//       for (kj=0; kj<del; ++kj)
// 	of << ghost_points[ki*del+kj] << " ";
//       of << std::endl;
//     }
// #endif

  
// }

//==============================================================================
void LRSurfApprox::constructInnerGhostPoints()
//==============================================================================
{
#ifdef DEBUG
  std::ofstream of("ghost_pnt2.g2");
#endif

  // Global information
  int dim = srf_->dimension();
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  int tot_nmb_pts = (int)points_.size()/(del-1);  // Distance not included
  double u1 = srf_->paramMin(XFIXED);
  double u2 = srf_->paramMax(XFIXED);
  double v1 = srf_->paramMin(YFIXED);
  double v2 = srf_->paramMax(YFIXED);
  double pdsize = (u2-u1)*(v2-v1);
  double nmb_fac = 0.2;
   
  LRSplineSurface::ElementMap::const_iterator it;
  for (it=srf_->elementsBegin(); it != srf_->elementsEnd(); ++it)
    {
      double umin = it->second->umin();
      double umax = it->second->umax();
      double vmin = it->second->vmin();
      double vmax = it->second->vmax();
      int nmb = it->second->nmbDataPoints();
      int nmb2 = it->second->nmbGhostPoints();

      // Check if the element has too few data points. In that
      // case, construct ghost points to stabilize the construction
      int nmb_init;  // Expected number of data points provided a
      // balanced distribution
      double pdsize2 = (umax-umin)*(vmax-vmin);
      nmb_init = (int)(tot_nmb_pts*pdsize2/pdsize);
      if (nmb+nmb2 < std::min(25, (int)(nmb_fac*nmb_init)))
	{
	  vector<double> ghost_points;
	  nmb_init -= (nmb + nmb2);
	  nmb_init = (int)sqrt((double)nmb_init);
	  int nmb_ghost_u = (int)(nmb_fac*nmb_init*(umax-umin)/(vmax-vmin));
	  int nmb_ghost_v = (int)(nmb_fac*nmb_init*(vmax-vmin)/(umax-umin));
	  nmb_ghost_u = std::max(3,std::min(nmb_ghost_u, 10));
	  nmb_ghost_v = std::max(3,std::min(nmb_ghost_v, 10));
	  if (grid_)
	    {
	      nmb_ghost_u = std::max(1, std::min(nmb_ghost_u, (int)((umax-umin)/cell_size_[0])));
	      nmb_ghost_v = std::max(1, std::min(nmb_ghost_v, (int)((vmax-vmin)/cell_size_[1])));
	    }

	  ghost_points.reserve(nmb_ghost_u*nmb_ghost_v*del);
	  double u_del = (umax-umin)/(double)nmb_ghost_u;
	  double v_del = (vmax-vmin)/(double)nmb_ghost_v;
	  double upar, vpar;
	  int ki, kj;
	  for (kj=0, vpar=vmin+0.5*v_del; kj<nmb_ghost_v; ++kj, vpar+=v_del)
	    for (ki=0, upar=umin+0.5*u_del; ki<nmb_ghost_u; ++ki, upar+=u_del)
	      {
		Point pos;
		srf_->point(pos, upar, vpar);
		ghost_points.push_back(upar);
		ghost_points.push_back(vpar);
		ghost_points.insert(ghost_points.end(), pos.begin(),
				    pos.end());
		ghost_points.push_back(0.0);
	      }
	  it->second->addGhostPoints(ghost_points.begin(), ghost_points.end(),
				     false);

#ifdef DEBUG
	  int nmb_ghost = (int)ghost_points.size()/del;
	  of << "400 1 0 4 100 55 100 255" << std::endl;
	  of << nmb_ghost << std::endl;
	  
	  for (ki=0; ki<nmb_ghost; ++ki)
	    {
	      for (kj=0; kj<del-1; ++kj)
		of << ghost_points[ki*del+kj] << " ";
	      of << std::endl;
	    }
#endif
	}
    }
}





//==============================================================================
void LRSurfApprox::updateGhostElems(vector<Element2D*>& elems, bool enable_omp)
//==============================================================================
{
  // Update corresponding coefficients using LR-MBA
  vector<Element2D*> elems2;  // Elements influence by update
  LRSplineMBA::MBAUpdate(srf_.get(), elems, elems2);

  // Recompute ghost points
  updateGhostPoints(elems);

  // Update distances in influenced elements
  RectDomain rd = srf_->containingDomain();
  int dim = srf_->dimension();
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point

  for (size_t ki=0; ki<elems2.size(); ++ki)
    {
      vector<double>& points = elems2[ki]->getDataPoints();
      int nmb_pts = elems2[ki]->nmbDataPoints();

      // Compute distances in data points and update parameter pairs
      // if requested
      if (nmb_pts > 0)
	{
	    if (enable_omp)
	    {
		computeAccuracyElement_omp(points, nmb_pts, del, rd, elems2[ki]);
	    }
	    else
	    {
		computeAccuracyElement(points, nmb_pts, del, rd, elems2[ki]);
	    }

	  // Local error information
	  double max_err = 0.0;
	  double av_err = 0.0;
	  double acc_err = 0.0;
	  int outside = 0;

	  double *curr;
	  int kj;
	  for (kj=0, curr=&points[0]; kj<nmb_pts; ++kj, curr+=del)
	    {
	      // Accumulate approximation error
	      double dist2 = fabs(curr[del-1]);
	      max_err = std::max(max_err, dist2);
	      acc_err += dist2;
	      if (dist2 > aepsge_)
		{
		  outside++;
		  av_err += dist2;
		}
	    }
	  if (outside > 0)
	    av_err /= (double)outside;

	  // Store updated accuracy information in the element
	  elems2[ki]->setAccuracyInfo(acc_err, av_err, max_err, outside);
	}
    }
}

//==============================================================================
void LRSurfApprox::updateGhostPoints(vector<Element2D*>& elems)
//==============================================================================
{
  // Global information
  int dim = srf_->dimension();
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  int tot_nmb_pts = (int)points_.size()/(del-1);  // Distance not included
  double u1 = srf_->paramMin(XFIXED);
  double u2 = srf_->paramMax(XFIXED);
  double v1 = srf_->paramMin(YFIXED);
  double v2 = srf_->paramMax(YFIXED);
  double pdsize = (u2-u1)*(v2-v1);
  double nmb_fac = 0.2;
   
  for (size_t ki=0; ki<elems.size(); ++ki)
    {
      double umin = elems[ki]->umin();
      double umax = elems[ki]->umax();
      double vmin = elems[ki]->vmin();
      double vmax = elems[ki]->vmax();
      int nmb = elems[ki]->nmbDataPoints();

      // Remove old ghost points
      elems[ki]->eraseGhostPoints();

      // Construct new ones
      int nmb_init;  // Expected number of data points provided a
      // balanced distribution
      double pdsize2 = (umax-umin)*(vmax-vmin);
      nmb_init = (int)(tot_nmb_pts*pdsize2/pdsize);
      nmb_init -= nmb;

      if (nmb_init <= 0)
	continue;  // No new ghost points
      
      nmb_init = (int)sqrt((double)nmb_init);
      int nmb_ghost_u = (int)(nmb_fac*nmb_init*(umax-umin)/(vmax-vmin));
      int nmb_ghost_v = (int)(nmb_fac*nmb_init*(vmax-vmin)/(umax-umin));
      nmb_ghost_u = std::max(3,std::min(nmb_ghost_u, 10));
      nmb_ghost_v = std::max(3,std::min(nmb_ghost_v, 10));

      vector<double> ghost_points;
      ghost_points.reserve(nmb_ghost_u*nmb_ghost_v*del);
      double u_del = (umax-umin)/(double)nmb_ghost_u;
      double v_del = (vmax-vmin)/(double)nmb_ghost_v;
      double upar, vpar;
      int kr, kj;
      for (kj=0, vpar=vmin+0.5*v_del; kj<nmb_ghost_v; ++kj, vpar+=v_del)
	for (kr=0, upar=umin+0.5*u_del; kr<nmb_ghost_u; ++kr, upar+=u_del)
	  {
	    Point pos;
	    srf_->point(pos, upar, vpar);
	    ghost_points.push_back(upar);
	    ghost_points.push_back(vpar);
	    ghost_points.insert(ghost_points.end(), pos.begin(),
				pos.end());
	    ghost_points.push_back(0.0);
	  }
      elems[ki]->addGhostPoints(ghost_points.begin(), ghost_points.end(),
				false);
    }

}

//==============================================================================
void LRSurfApprox::addConstraintGhostPoints()
//==============================================================================
{
  int dim = srf_->dimension();
  if (dim != 1)
    return;  // Only for 1D case

  if (!(has_min_constraint_ || has_max_constraint_))
    return;  // Nothing to do

#ifdef DEBUG
  std::ofstream of("ghost_pnt3.g2");
#endif

  // Global information
  int del = 3 + dim;  // Parameter pair, position and distance between surface and point
  int nmb_ghost_u = 2;
  int nmb_ghost_v = 2;
  LRSplineSurface::ElementMap::const_iterator it;
  for (it=srf_->elementsBegin(); it != srf_->elementsEnd(); ++it)
    {
      // Fetch element domain
      double umin = it->second->umin();
      double umax = it->second->umax();
      double vmin = it->second->vmin();
      double vmax = it->second->vmax();

      int nmbg = it->second->nmbGhostPoints();
      if (nmbg > 100)
	continue;  // Ghost points already created

       // Fetch associated B-splines
      const vector<LRBSpline2D*>& bsplines = it->second->getSupport();

      // Check if the bspline coefficient exceeds the limits in either direction
      bool exceeds_min = false, exceeds_max = false;
      for (size_t kr=0; kr<bsplines.size(); ++kr)
	{
	  Point coef = bsplines[kr]->Coef();
	  if (has_min_constraint_ && coef[0] < minval_)
	    exceeds_min = true;
	  if (has_max_constraint_ && coef[0] > maxval_)
	    exceeds_max = true;
	}
	   
      vector<double> ghost_points;
      if (has_min_constraint_ && exceeds_min && has_max_constraint_ && exceeds_max)
	{
	  // Add ghost points at medium value
	  double val = 0.5*(minval_ + maxval_);
	  ghost_points.reserve(nmb_ghost_u*nmb_ghost_v*del);
	  double u_del = (umax-umin)/(double)nmb_ghost_u;
	  double v_del = (vmax-vmin)/(double)nmb_ghost_v;
	  double upar, vpar;
	  int ki, kj;
	  for (kj=0, vpar=vmin+0.5*v_del; kj<nmb_ghost_v; ++kj, vpar+=v_del)
	    for (ki=0, upar=umin+0.5*u_del; ki<nmb_ghost_u; ++ki, upar+=u_del)
	      {
		Point pos;
		srf_->point(pos, upar, vpar);
		ghost_points.push_back(upar);
		ghost_points.push_back(vpar);
		ghost_points.push_back(val);
		ghost_points.push_back(val-pos[0]);
	      }
	}

      else if (has_min_constraint_ && exceeds_min)
	{
	  // Add ghost points at minimum value
	  ghost_points.reserve(nmb_ghost_u*nmb_ghost_v*del);
	  double u_del = (umax-umin)/(double)nmb_ghost_u;
	  double v_del = (vmax-vmin)/(double)nmb_ghost_v;
	  double upar, vpar;
	  int ki, kj;
	  for (kj=0, vpar=vmin+0.5*v_del; kj<nmb_ghost_v; ++kj, vpar+=v_del)
	    for (ki=0, upar=umin+0.5*u_del; ki<nmb_ghost_u; ++ki, upar+=u_del)
	      {
		Point pos;
		srf_->point(pos, upar, vpar);
		ghost_points.push_back(upar);
		ghost_points.push_back(vpar);
		ghost_points.push_back(minval_);
		ghost_points.push_back(minval_-pos[0]);
	      }
	}

      else if (has_max_constraint_ && exceeds_max)
	{
	  // Add ghost points at maximum value
	  ghost_points.reserve(nmb_ghost_u*nmb_ghost_v*del);
	  double u_del = (umax-umin)/(double)nmb_ghost_u;
	  double v_del = (vmax-vmin)/(double)nmb_ghost_v;
	  double upar, vpar;
	  int ki, kj;
	  for (kj=0, vpar=vmin+0.5*v_del; kj<nmb_ghost_v; ++kj, vpar+=v_del)
	    for (ki=0, upar=umin+0.5*u_del; ki<nmb_ghost_u; ++ki, upar+=u_del)
	      {
		Point pos;
		srf_->point(pos, upar, vpar);
		ghost_points.push_back(upar);
		ghost_points.push_back(vpar);
		ghost_points.push_back(maxval_);
		ghost_points.push_back(maxval_-pos[0]);
	      }
	}

      if (ghost_points.size() > 0)
	{
	  it->second->addGhostPoints(ghost_points.begin(), ghost_points.end(),
				     false);

#ifdef DEBUG
	int nmb_ghost = (int)ghost_points.size()/del;
	of << "400 1 0 4 100 55 100 255" << std::endl;
	of << nmb_ghost << std::endl;
	  
	for (int ki=0; ki<nmb_ghost; ++ki)
	  {
	    for (int kj=0; kj<del-1; ++kj)
	      of << ghost_points[ki*del+kj] << " ";
	    of << std::endl;
	  }
#endif
	int stop_break = 1;
	}
    }
  int stop_break2 = 1;
}

//==============================================================================
void LRSurfApprox::adaptSurfaceToConstraints()
//==============================================================================
{
  if (srf_->dimension() != 1)
    return;  // Only applicable for functions

  if (!(has_min_constraint_ || has_max_constraint_ || has_local_constraint_))
    return;

  LRSplineSurface::BSplineMap::const_iterator it1 = srf_->basisFunctionsBegin();
  for (; it1 != srf_->basisFunctionsEnd(); ++it1)
    {
      Point coef = it1->second->Coef();
      if (has_min_constraint_)
	coef[0] = std::max(coef[0], minval_);
      if (has_max_constraint_)
	coef[0] = std::min(coef[0], maxval_);
      srf_->setCoef(coef, it1->second.get());
      if (has_local_constraint_)
	{
	  double bb[2], curr_bb[2];  // Bounding box
	  bb[0] = HUGE;
	  bb[1] = -HUGE;
	  
	  // For all elements
	  const vector<Element2D*>& elem = it1->second->supportedElements();
	  for (size_t ki=0; ki<elem.size(); ++ki)
	    {
	      bool found = elem[ki]->getDataBoundingBox(curr_bb);
	      if (found)
		{
		  bb[0] = std::min(bb[0], curr_bb[0]);
		  bb[1] = std::max(bb[1], curr_bb[1]);
		}
	    }
	  if (bb[1] >= bb[0])
	    {
	      coef[0] = std::max(coef[0], 
				 bb[0] - constraint_fac_*(bb[1] - bb[0]));
	      coef[0] = std::min(coef[0], 
				 bb[1] + constraint_fac_*(bb[1] - bb[0]));
	    }
	}
    }
}

//==============================================================================
void LRSurfApprox::turnTo3D()
//==============================================================================
{
  if (!srf_->dimension() == 1)
    return;  // Not possible to make 3D surface

  // Make data points 3D
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      it->second->makeDataPoints3D();
    }
  
  // Turn surface into 3D
  srf_->to3D();
}
