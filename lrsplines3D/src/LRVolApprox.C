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

#include "GoTools/lrsplines3D/LRVolApprox.h"
#include "GoTools/lrsplines3D/LRBSpline3DUtils.h"
//#include "GoTools/lrsplines3D/LinDepUtils.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/LRSpline3DMBA.h"
#include "GoTools/lrsplines3D/LRSpline3DUtils.h"
#include "GoTools/lrsplines3D/LRFeature3DUtils.h"
//#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/PointCloud.h"
//#include "GoTools/lrsplines3D/LRSplinePlotUtils.h"
#include "GoTools/trivariate/SmoothVolume.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

//#define DEBUG
//#define DEBUG0

using std::vector;
using std::cout;
using std::endl;
using std::pair;
using namespace Go;

//==============================================================================
LRVolApprox::LRVolApprox(vector<double>& points, 
                         int dim, double epsge,  
                         double mba_level,
                         bool closest_dist, bool repar)
  : nmb_pts_((int)points.size()/(3+dim)), points_(points), useMBA_(true), 
    initMBA_(true), initMBA_coef_(mba_level), 
    maxdist_(-10000.0), maxdist_prev_(-10000.0), 
    avdist_(0.0), avdist_all_(0), avdist_all_prev_(0), 
    outsideeps_(0), outsideeps_prev_(0), maxout_(-10000.0), avout_(0.0),
    aepsge_(epsge), smoothweight_(0.0),
    repar_(repar), check_close_(closest_dist), 
    fix_corner_(false), to3D_(false), grid_(false), check_init_accuracy_(false),
    initial_volume_(false), has_min_constraint_(false), has_max_constraint_(false),
    has_local_constraint_(false), outfrac_(0.0), write_feature_(false), verbose_(true)
//==============================================================================
{
  face_derivs_[0] = face_derivs_[1] 
    = face_derivs_[2] = face_derivs_[3] 
    = face_derivs_[4] = face_derivs_[5] = 0;

  grid_start_[0] = grid_start_[1] = grid_start_[2] = 0.0;
  cell_size_[0] = cell_size_[1] = cell_size_[2] = 1.0;
  usize_min_ = vsize_min_ = wsize_min_ = -1;

  fix_boundary_ = false; //true;

  makeInitVol(dim);
}
/*
//==============================================================================
LRVolApprox::LRVolApprox(shared_ptr<SplineSurface>& srf,
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
  face_derivs_[0] = face_derivs_[1] = face_derivs_[2] = face_derivs_[3] = 0;
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
LRVolApprox::LRVolApprox(shared_ptr<LRSplineSurface>& srf,
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
  face_derivs_[0] = face_derivs_[1] = face_derivs_[2] = face_derivs_[3] = 0;
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
LRVolApprox::LRVolApprox(int ncoef_u, int order_u, int ncoef_v, int order_v,
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
  face_derivs_[0] = face_derivs_[1] = face_derivs_[2] = face_derivs_[3] = 0;
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
*/
/*
//==============================================================================
LRVolApprox::LRVolApprox(int order_u, vector<double>& knots_u, 
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
  face_derivs_[0] = face_derivs_[1] = face_derivs_[2] = face_derivs_[3] = 0;
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
*/

//==============================================================================
LRVolApprox::LRVolApprox(int ncoef_u, int order_u, int ncoef_v, int order_v,
			 int ncoef_w, int order_w,
			 vector<double>& points, int dim, 
			 double domain[], double epsge, //bool init_mba, 
			 double mba_level,
			 bool closest_dist, bool repar)
  //==============================================================================
  : nmb_pts_((int)points.size()/(3+dim)), points_(points), useMBA_(false),
    initMBA_(true), initMBA_coef_(mba_level), 
    maxdist_(-10000.0), maxdist_prev_(-10000.0), avdist_(0.0), 
    avdist_all_(0.0), avdist_all_prev_(0), outsideeps_(0), outsideeps_prev_(0),
    maxout_(-10000.0), avout_(0.0), aepsge_(epsge), 
    smoothweight_(0.0), repar_(repar), check_close_(closest_dist), 
    fix_corner_(false), to3D_(false), grid_(false), check_init_accuracy_(false),
    initial_volume_(false), has_min_constraint_(false), has_max_constraint_(false),
    has_local_constraint_(false), outfrac_(0.0), write_feature_(false), verbose_(true)
{
  face_derivs_[0] = face_derivs_[1] = face_derivs_[2] = face_derivs_[3] = 0;
  face_derivs_[4] = face_derivs_[5] = 0;
  grid_start_[0] = grid_start_[1] = grid_start_[2] =0.0;
  cell_size_[0] = cell_size_[1] = cell_size_[2] = 1.0;
  usize_min_ = vsize_min_ = wsize_min_ = -1;

  fix_boundary_ = false; //true;

  // Create an LR B-spline volume with unset coefficients. 
  // The size of the spline space is given. The knots will be equally spaced
  vector<double> knots_u(ncoef_u+order_u);
  vector<double> knots_v(ncoef_v+order_v);
  vector<double> knots_w(ncoef_w+order_w);
  double udel = (domain[1] - domain[0])/(double)(ncoef_u - order_u + 1);
  double vdel = (domain[3] - domain[2])/(double)(ncoef_v - order_v + 1);
  double wdel = (domain[5] - domain[4])/(double)(ncoef_w - order_w + 1);
  int kj;
  for (kj=0; kj<order_u; ++kj)
    {
      knots_u[kj] = domain[0];
      knots_u[ncoef_u+kj] = domain[1];
    }
    for (double upar=domain[0]+udel; kj<ncoef_u; ++kj, upar+=udel)
      knots_u[kj] = upar;

    for (kj=0; kj<order_v; ++kj)
    {
      knots_v[kj] = domain[2];
      knots_v[ncoef_v+kj] = domain[3];
    }
    for (double vpar=domain[2]+vdel; kj<ncoef_v; ++kj, vpar+=vdel)
      knots_v[kj] = vpar;

  for (kj=0; kj<order_w; ++kj)
    {
      knots_w[kj] = domain[4];
      knots_w[ncoef_w+kj] = domain[5];
    }
    for (double wpar=domain[4]+wdel; kj<ncoef_w; ++kj, wpar+=wdel)
      knots_w[kj] = wpar;

  makeInitVol(dim, ncoef_u, order_u, ncoef_v, order_v, ncoef_w, order_w,
	      &knots_u[0], &knots_v[0], &knots_w[0]);
}


//==============================================================================
LRVolApprox::~LRVolApprox()
//==============================================================================
{
}

//==============================================================================
shared_ptr<LRSplineVolume> LRVolApprox::getApproxVol(double& maxdist, 
                                                     double& avdist_all,
                                                     double& avdist,
                                                     int& nmb_out_eps, 
                                                     int& iterations)
//==============================================================================
{
//   // We start the timer.
// #ifdef _OPENMP
//   double time0 = omp_get_wtime();
// #endif

#ifdef _OPENMP
  std::cout << "Open MP" << std::endl;
  
    // When using OpenMP we choose between splitting the threads on
    // the surface elements or on the points for each element.  As the
    // iteration progresses we should turn towards splitting on the
    // elements.  Initial switch threshold set to num_elem ==
    // avg_num_pnts_per_elem.
    const int num_elem = vol_->numElements();
    const int num_pts = points_.size()/(vol_->dimension());
    // We let the number of elem vs average numer of points per elem be the threshold
    // for switching the OpenMP level.
    const double pts_per_elem = num_pts/num_elem;
    const bool omp_for_elements = true; //false; //(num_elem > pts_per_elem); // As opposed to element points.
    const bool omp_for_mba_update = true;
#ifndef NDEBUG
    std::cout << "num_elem: " << num_elem << ", pts_per_elem: " << pts_per_elem << ", openmp_for_elements: " <<
	omp_for_elements << std::endl;
#endif
#else
    const bool omp_for_elements = false; // 201503 The omp version seems to be faster even when run sequentially.
    const bool omp_for_mba_update = false;
#endif

#ifdef DEBUG
  std::ofstream of0("init0_vol.g2");
  shared_ptr<LRSplineVolume> tmp0(vol_->clone());
  if (tmp0->dimension() == 1)
    tmp0->to3D();
  tmp0->writeStandardHeader(of0);
  tmp0->write(of0);
  of0 << std::endl;
  // LineCloud lines0 = tmp0->getElementBds();
  // lines0.writeStandardHeader(of0);
  // lines0.write(of0);
  std::ofstream of02("init0_tpvol.g2");
  shared_ptr<SplineVolume> svol0(tmp0->asSplineVolume());
  svol0->writeStandardHeader(of02);
  svol0->write(of02);
#endif
  /* @obar // IS THIS NEEDED???
  if (vol_->dimension() == 3)
    {
      // Reparameterize to reflect the volume size
      double len1, len2;
      vol_->estimateSfSize(len1, len2);

      double umin = vol_->paramMin(XFIXED);
      double umax = vol_->paramMax(XFIXED);
      double vmin = vol_->paramMin(YFIXED);
      double vmax = vol_->paramMax(YFIXED);
      vol_->setParameterDomain(umin, umin+len1, vmin, vmin+len2);

      // Reparameterize also data points
      int del = 5;  // Parameter pair + geometric dimension
      int nmb = (int)points_.size()/del;
      for (int kj=0; kj<nmb; ++kj)
	{
	  points_[kj*del] = umin + (points_[kj*del] - umin)*len1/(umax - umin);
	  points_[kj*del+1] = vmin + (points_[kj*del+1] - vmin)*len2/(vmax - vmin);
	}
    }
  */
  
  // Distribute given data points to elements
  LRSpline3DUtils::distributeDataPoints(vol_.get(), points_, true, true);

  double mineps = aepsge_;
  for (size_t ki=0; ki<tolerances_.size(); ++ki)
    {
#ifdef DEBG
      std::cout << tolerances_[ki].box.umin() << " " << tolerances_[ki].box.umax() << " ";
      std::cout << tolerances_[ki].box.vmin() << " " << tolerances_[ki].box.vmax() << " ";
      std::cout << tolerances_[ki].box.wmin() << " " << tolerances_[ki].box.wmax() << " ";
      std::cout << tolerances_[ki].tol << std::endl;
#endif
      mineps = std::min(mineps, tolerances_[ki].tol);
    }
  double delta = 0.0; //0.25/mineps;
#ifdef DEBUG
  std::cout << "Aepsge: " << aepsge_ << ", mineps: " << mineps << std::endl;
#endif
  
  vector<Element3D*> ghost_elems;
  if (check_init_accuracy_ /*|| useMBA_*/)
  {
    // Compute accuracy in data points
      if (omp_for_elements)
          computeAccuracy_omp(ghost_elems);
      else
          computeAccuracy(ghost_elems);
  }


  // Initiate approximation engine
  //if (fix_corner_)
  //  setCoefKnown();

  //if (fix_boundary_)
  //  setFixBoundary(true);

  // Initial approximation of LR B-spline surface
  if (omp_for_mba_update && vol_->dimension() == 1)
    {
      // double delta = (avdist_ < mineps) ? 0.0 : 0.001;
      // delta = 0.005;
      LRSpline3DMBA::MBADistAndUpdate_omp(vol_.get(), mineps, delta);
    }
  else
    {
      LRSpline3DMBA::MBADistAndUpdate(vol_.get());
    }
  //if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
  //  {
  //    adaptSurfaceToConstraints();
  //  }
  int m_itermax = 2; //3; //4;
  for (int m_iter=1; m_iter<m_itermax; ++m_iter)
    {
#ifdef DEBUG
  cout << "Running MBA a second time... " << endl;
#endif
  if (omp_for_mba_update && vol_->dimension() == 1)
    {
      // double delta = (avdist_ < mineps) ? 0.0 : 0.001;
      // delta = 0.005;
      LRSpline3DMBA::MBADistAndUpdate_omp(vol_.get(), mineps, delta);
    }
  else
    {
      LRSpline3DMBA::MBADistAndUpdate(vol_.get());
    }
    }
  // if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
  //  {
  //    adaptSurfaceToConstraints();
  //  }
  updateCoefKnown();
  
#ifdef DEBUG
    std::ofstream of1("init_vol.g2");
    shared_ptr<LRSplineVolume> tmp;
    if (vol_->dimension() == 1)
      {
        tmp = shared_ptr<LRSplineVolume>(vol_->clone());
        tmp->to3D();
      }
    else
      tmp = vol_;
    tmp->writeStandardHeader(of1);
    tmp->write(of1);
    std::ofstream of12("init_tpvol.g2");
    shared_ptr<SplineVolume> svol1(tmp->asSplineVolume());
    svol1->writeStandardHeader(of12);
    svol1->write(of12);
    of12 << std::endl;
    cout << "Need to implement PlaneCloud if you want this..." << endl;
    //LineCloud lines = tmp->getElementBds();
    //lines.writeStandardHeader(of12);
    //lines.write(of12);
#endif

#ifdef DEBUG0
    std::cout << "Compute accuracy" << std::endl;
#endif

    // if (verbose_)
    //   {
	// Compute accuracy in data points
	if (omp_for_elements)
	  computeAccuracy_omp(ghost_elems);
	else
	  computeAccuracy(ghost_elems);
      // }

    // Compute accuracy in data points
    if (verbose_)
      {
        std::cout << "Number of data points: " << nmb_pts_ << std::endl;
        std::cout << "Number of coefficients: " << vol_->numBasisFunctions() << std::endl;
        std::cout << "Initial volume. Maximum distance: " << maxdist_;
        std::cout << ", average distance: " << avdist_all_ << std::endl;
        std::cout << "Number of points outside tolerance: " << outsideeps_;
        std::cout << ", average distance in outside points: " << avdist_ << std::endl;
	std::cout << "Maximum distance exceeding tolerance (dist-tol): " << maxout_ << std::endl;
	std::cout << "Average distance exceeding tolerance (dist-tol): " << avout_ << std::endl;
      }

  if (write_feature_)
    {
      bool feature = true;
      if (feature_levels_.size() > 0)
	{
	  size_t kl;
	  for (kl=0; kl<feature_levels_.size(); ++kl)
	    if (feature_levels_[kl] == 0)
	      break;
	  if (kl == feature_levels_.size())
	    feature = false;
	}

      if (feature)
	{
	  std::ofstream f_out("cellinfo0.txt");
	  LRFeature3DUtils::writeCellInfo(*vol_, aepsge_, ncell1_, ncell2_, ncell3_, f_out);
	}
    }

    ghost_elems.clear();
    //points_.clear();  // Not used anymore TESTING

    int level;
    // Start the iteration
    double threshold_prev = -1.0;
    double eps = 1.0-12;
    double epsfrac = 0.001*mineps;
    bool recompute_accuracy = false;
    for (level=0; level<iterations; ++level)
      {
#ifdef DEBUG0
	std::cout << "Current level = " << level << std::endl;
#endif
	// Check if the requested accuracy is reached
	if (maxout_ <= eps || outsideeps_ == 0)
	  break;


	// Refine volume
	if (maxout_ < epsfrac)
	  prev_ =  shared_ptr<LRSplineVolume>(vol_->clone());

	// OBAR
	//      int nn = 4;
	//      double hu = (vol_->paramMax(XDIR)-vol_->paramMin(XDIR))/nn;
	//      double hv = (vol_->paramMax(YDIR)-vol_->paramMin(YDIR))/nn;
	//      double hw = (vol_->paramMax(ZDIR)-vol_->paramMin(ZDIR))/nn;
	//      cout << "Elements: "  << vol_->numElements() << endl;
	//      for (int ix=0; ix!=nn; ++ix)
	//        {
	//          double ux = vol_->paramMin(XDIR)+hu*ix;
	//          for (int jx=0; jx!=nn; ++jx)
	//            {
	//              double vx = vol_->paramMin(YDIR)+hv*jx;
	//              for (int kx=0; kx!=nn; ++kx)
	//                {
	//                  double wx = vol_->paramMin(ZDIR)+hw*kx;
	//                  Point pt;
	//                  vol_->point(pt,ux,vx,wx);
	//                  cout << pt << " ";
	//                }
	//            }
	//        } cout << endl;

	// Check if any ghost points need to be updated
	//if (!useMBA_ && ki<toMBA_ && ghost_elems.size() > 0)
        //{
        //  updateGhostElems(ghost_elems);
        //}

#ifdef DEBUG0
	std::cout << "Refine" << std::endl;
#endif
	if (level > 0 || (!initial_volume_))
	  {
	    // double threshold = ((iterations-level-1)*maxdist_ +
	    // 		      (level+1)*aepsge_)/(double)iterations;
	    //double threshold = std::max(aepsge_, 0.5*maxdist_);
	    //double threshold = (mineps + avdist_ + maxdist_)/3.0;
	    double threshold = (mineps + 0.5*maxout_)/2.0 + avout_;
	    if (level > 0 && maxdist_/maxdist_prev_ > 0.9)
	      {
		// Slow convergence. Reduce threshold
		threshold = (mineps + avdist_ + 0.5*maxdist_)/3.0;
	      }
	    if (threshold_prev > 0.0 && threshold/threshold_prev > 0.9)
	      threshold = 0.9*threshold_prev;
	    threshold = std::max(mineps, threshold);
	    //threshold = aepsge_;
#ifdef DEBUG0
	    std::cout << "Level " << level+1 << ", threshold = " << threshold << std::endl;
#endif
	    int nmb_refs = refineVol(threshold);
	    if (nmb_refs == 0)
	      break;  // No refinements performed
#ifdef DEBUG0
	    std::cout << "Distribute data points" << std::endl;
#endif
	    LRSpline3DUtils::distributeDataPoints(vol_.get(), points_, 
						  true, true);
	    threshold_prev = threshold;
	  }

	// OBAR
	//      cout << "Elements: "  << vol_->numElements() << endl;
	//      for (int ix=0; ix!=nn; ++ix)
	//        {
	//          double ux = vol_->paramMin(XDIR)+hu*ix;
	//          for (int jx=0; jx!=nn; ++jx)
	//            {
	//              double vx = vol_->paramMin(YDIR)+hv*jx;
	//              for (int kx=0; kx!=nn; ++kx)
	//                {
	//                  double wx = vol_->paramMin(ZDIR)+hw*kx;
	//                  Point pt;
	//                  vol_->point(pt,ux,vx,wx);
	//                  cout << pt << " ";
	//                }
	//            }
	//        } cout << endl;
    
	//#ifdef DEBUG
	//      std::ofstream of2("refined_vol.g2");
	//      shared_ptr<LRSplineVolume> tmp2;
	//      if (vol_->dimension() == 1)
	//	{
	//	  tmp2 = shared_ptr<LRSplineVolume>(vol_->clone());
	//	  tmp2->to3D();
	//	}
	//      else
	//	tmp2 = vol_;
	//      tmp2->writeStandardHeader(of2);
	//      tmp2->write(of2);
	//      of2 << std::endl;
	//      cout << "Need to implement PlaneCloud if you want this..." << endl;
	//      //LineCloud lines2 = tmp2->getElementBds();
	//      //lines2.writeStandardHeader(of2);
	//      //lines2.write(of2);

	//      std::ofstream of3("point_clouds.g2");
	//      int del = vol_->dimension() + 4; // Parameter triple, position and
	//      // distance between volume and point
	//      for (LRSplineVolume::ElementMap::const_iterator it=vol_->elementsBegin();
	//	   it != vol_->elementsEnd(); ++it)
	//	{
	//	  vector<double>& elem_data = it->second->getDataPoints();
	//	  int nmb = (int)elem_data.size()/del;
	//	  if (elem_data.size() > 0)
	//	    {
	//	      vector<double> tmppt;
	//	      if (vol_->dimension() == 1)
	//		tmppt = elem_data;
	//	      else
	//		{
	//		  tmppt.reserve(3*elem_data.size()/del);
	//		  for (int kr=0; kr<nmb; ++kr)
	//		    tmppt.insert(tmppt.end(), elem_data.begin()+kr*del+2,
	//				 elem_data.begin()+(kr+1)*del-1);
		  
	//		}
	//	      PointCloud4D cloud(tmppt.begin(), nmb);
	//	      cloud.writeStandardHeader(of3);
	//	      cloud.write(of3);
	//	    }
	//	}
	//#endif

	// Update coef_known from information in LR B-splines
	updateCoefKnown();
	//unsetCoefKnown(); TEST
	//if (fix_corner_)
	//        setCoefKnown();
	//if (fix_boundary_)
	//	setFixBoundary(true);
  
	// Check for linear independence (overloading)
	// NOT AVAILABLE YET ALTHOUGH NOT NEEDED FOR MBA
	//vector<LRBSpline3D*> funs = LinDepUtils::unpeelableBasisFunctions(*vol_);
#ifdef DEBUG
	//std::cout << "Number of unpeelable functions: " << funs.size() << std::endl;
#endif
      
#ifdef DEBUG0
	std::cout << "MBA 1" << std::endl;
#endif

	// Update surface
	if (vol_->dimension() == 3)
	  {
	    LRSpline3DMBA::MBADistAndUpdate(vol_.get());
	  }
	else if (omp_for_mba_update)
	  {
	    // double delta = (avdist_ < mineps) ? 0.0 : 0.001;
	    // delta = 0.005;
	    LRSpline3DMBA::MBADistAndUpdate_omp(vol_.get(), mineps,
						maxout_ < mineps ? 0.0 : delta);
	  }
	else
	  {
	    LRSpline3DMBA::MBADistAndUpdate(vol_.get());
	    //LRSpline3DMBA::MBAUpdate(vol_.get());
	  }
	//if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
	//  adaptSurfaceToConstraints();
	for (int m_iter=1; m_iter<m_itermax; ++m_iter)
	  {
#ifdef DEBUG0
	    std::cout << "MBA 2" << std::endl;
#endif

	    if (omp_for_mba_update && vol_->dimension() == 1)
	      {
		// double delta = (avdist_ < mineps) ? 0.0 : 0.001;
		// delta = 0.005;
		LRSpline3DMBA::MBADistAndUpdate_omp(vol_.get(), mineps, 
						maxout_ < mineps ? 0.0 : delta);
	      }
	    else
	      {
		LRSpline3DMBA::MBADistAndUpdate(vol_.get());
	      }
	  }
	// computeAccuracy();
	// LRSpline3DMBA::MBAUpdate(vol_.get());
	//if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
	//  adaptSurfaceToConstraints();
      
      
	/*#ifdef DEBUG
	  std::ofstream of4("updated_vol.g2");
	  #endif*/
	/*      shared_ptr<LRSplineVolume> tmp3;
		if (vol_->dimension() == 1)
		{
		tmp3 = shared_ptr<LRSplineVolume>(vol_->clone());
		tmp3->to3D();
		}
		else
		tmp3 = vol_;*/
	/*#ifdef DEBUG
	  tmp3->writeStandardHeader(of4);
	  tmp3->write(of4);
	  cout << "Need to implement PlaneCloud if you want this..." << endl;
	  //LineCloud lines3 = tmp3->getElementBds();
	  //lines3.writeStandardHeader(of4);
	  //lines3.write(of4);
	  std::ofstream of42("updated_tpvol.g2");
	  shared_ptr<SplineVolume> ssf4(tmp3->asSplineVolume());
	  ssf4->writeStandardHeader(of42);
	  ssf4->write(of42);
	  of42 << std::endl;
	  //LineCloud lines32 = tmp3->getElementBds();
	  //lines32.writeStandardHeader(of42);
	  //lines32.write(of42);
	  #endif*/
        
	maxdist_prev_ = maxdist_;
	avdist_all_prev_ = avdist_all_;
	outsideeps_prev_ = outsideeps_;
           
#ifdef DEBUG0
	std::cout << "Compute accuracy" << std::endl;
#endif

	  ghost_elems.clear();
	  if (omp_for_elements)
	    computeAccuracy_omp(ghost_elems);
	  else
	    computeAccuracy(ghost_elems);
	  if (verbose_)
	    {
	      std::cout << std::endl;
	      std::cout << "Iteration number " << level+1 <<". Maximum distance: " << maxdist_;
	      std::cout << ", average distance: " << avdist_all_ << std::endl;
	      std::cout << "Number of points outside tolerance: " << outsideeps_;
	      std::cout << ", average distance in outside points: " << avdist_ << std::endl;
	      std::cout << "Number of coefficients: " << vol_->numBasisFunctions() << std::endl;
	      std::cout << "Maximum distance exceeding tolerance (dist-tol): " << maxout_ << std::endl;
	      std::cout << "Average distance exceeding tolerance (dist-tol): " << avout_ << std::endl;
	    }
#ifdef DEBUG0
	  std::cout << "epsfrac: " << epsfrac << std::endl;
#endif
	  if (prev_.get() && maxdist_ > maxdist_prev_ &&
	      outsideeps_ >= outsideeps_prev_)
	    {
	      vol_ = prev_;
	      recompute_accuracy = true;
	      break;
	    }
	  if (write_feature_)
	    {
	      bool feature = true;
	      if (feature_levels_.size() > 0)
		{
		  size_t kl;
		  for (kl=0; kl<feature_levels_.size(); ++kl)
		    if (feature_levels_[kl] == level+1)
		      break;
		  if (kl == feature_levels_.size())
		    feature = false;
		}
	      
	      if (feature)
		{
		  std::string body = "cellinfo";
		  std::string extension = ".txt";
		  std::string ver = std::to_string(level+1);
		  std::string outfile = body + ver + extension;
		  std::ofstream f_out2(outfile.c_str());
		  LRFeature3DUtils::writeCellInfo(*vol_, aepsge_, ncell1_,
						  ncell2_, ncell3_, f_out2);
		}
	    }
	  
      }

    if (recompute_accuracy)
      {
	LRSpline3DUtils::distributeDataPoints(vol_.get(), points_, 
					      true, true);
	if (omp_for_elements)
	  computeAccuracy_omp(ghost_elems);
	else
	  computeAccuracy(ghost_elems);
      }
      
    // Set accuracy information
    maxdist = maxdist_;
    avdist_all = avdist_all_;
    avdist = avdist_;
    nmb_out_eps = outsideeps_;
    iterations = level;
    // #ifdef _OPENMP
    //   double time1 = omp_get_wtime();
    //   double time_spent = time1 - time0;
    //   std::cout << "time_spent in getApproxVol: " << time_spent << std::endl;
    // #endif
  
    return vol_;
}


/*
//==============================================================================
void LRVolApprox::performSmooth(LRSurfSmoothLS *LSapprox)
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
*/

//==============================================================================
void LRVolApprox::computeAccuracy(vector<Element3D*>& ghost_elems)
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
  maxout_ = 0.0;
  avout_ = 0.0;

#ifdef _OPENMP
  const bool omp_for_element_pts = false; //true;
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
  vector<Element3D*> elem;
#endif


  int dim = vol_->dimension();
  int del = 4 + dim;  // Parameter triple, position and distance between surface and point
  LRSplineVolume::ElementMap::const_iterator it;
  int num = vol_->numElements();
  int kj;

  double ghost_fac = 0.8; // @obar WHAT IS THIS?
  ghost_elems.clear();

  for (it=vol_->elementsBegin(), kj=0; kj<num; ++it, ++kj)
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
      double wmin = it->second->wmin();
      double wmax = it->second->wmax();

      vector<double>& points = it->second->getDataPoints();
      assert(points.size() != 0 );//@@@obar remove
      //vector<double>& ghost_points = it->second->getGhostPoints();
      int nmb_pts = it->second->nmbDataPoints();
      //int nmb_ghost = it->second->nmbGhostPoints();

       // Local error information
      double max_err = 0.0;
      double av_err = 0.0;
      double acc_err = 0.0;
      double acc_outside = 0.0;
       int outside = 0;
      double acc_err_sgn = 0.0;
      double av_err_sgn = 0.0;

      double acc_err_pos = 0.0;
      double acc_err_neg = 0.0;
      int nmb_err_pos = 0;
      int nmb_err_neg = 0;

      // Check if the accuracy can have been changed
      const vector<LRBSpline3D*>& bsplines = it->second->getSupport();
      size_t nb;
      for (nb=0; nb<bsplines.size(); ++nb)
	if (!bsplines[nb]->coefFixed())
	  break;

      if (true /*nb < bsplines.size()*/) // || useMBA_ 
	{
	  // Compute distances in data points and update parameter pairs
	  // if requested
// #ifdef _OPENMP
// 	    double time0_part = omp_get_wtime();
// #endif
	  if (nmb_pts > 0)
	  {
	      if (omp_for_element_pts)
		  computeAccuracyElement_omp(points, nmb_pts, del, it->second.get());
	      else
		  computeAccuracyElement(points, nmb_pts, del, it->second.get());
	  }
	  
	  // Compute distances in ghost points
	  //if (nmb_ghost > 0 && !useMBA_)
	  //{
	  //   if (omp_for_element_pts)
	  //	  computeAccuracyElement_omp(ghost_points, nmb_ghost, del, rd, it->second.get());
	  //    else
	  //	  computeAccuracyElement(ghost_points, nmb_ghost, del, rd, it->second.get());
	  //}
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
      double tol;  // Local tolerance taking the possibility for a 
      // varying tolerance threshold into account
      for (ki=0, curr=&points[0]; ki<nmb_pts;)
	{
	  //Point curr_pt(curr+(dim==3)*3, curr+del-1);

	  // Accumulate approximation error
	  double dist = curr[del-1];
	  double dist2 = fabs(dist);
	  maxdist_ = std::max(maxdist_, dist2);
	  max_err = std::max(max_err, dist2);
	  acc_err += dist2;
	  acc_err_sgn += dist;
	  avdist_all_ += dist2;
	  tol = aepsge_;
	  for (size_t kr=0; kr<tolerances_.size(); ++kr)
	    {
	      if (tolerances_[kr].contains(curr[0], curr[1], curr[2]))
		{
		  tol = tolerances_[kr].tol;
		  break;
		}
	    }
	  maxout_ = std::max(maxout_, dist2-tol);
	  if (dist > tol)
	    {
	      acc_err_pos += (dist-tol);
	      nmb_err_pos++;
	    }
	  else if (dist < -tol)
	    {
	      acc_err_neg += ((-dist) - tol);
	      nmb_err_neg++;
	    }
	  if (dist2 > tol)
	    {
	      av_err_sgn += curr[del-1];
	      avdist_ += dist2;
	      outsideeps_++;
	      av_err += dist2;
	      outside++;
	      avout_ += (dist2-tol);
	      acc_outside += (dist2-tol);
		  
#ifdef DEBUG
	      // Accumulate error points
	      if (dist > 0)
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
	      if (curr[0] < umin || curr[0] > umax ||
		  curr[1] < vmin || curr[1] > vmax ||
		  curr[2] < wmin || curr[2] > wmax)
		{
		  // Find element
		  Element3D *elem = vol_->coveringElement(curr[0], curr[1], curr[2]);
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
      // it->second->getAccuracyInfo(av_prev, max_prev, nmb_out_prev);

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

      // if (max_err > aepsge_ && max_prev > 0.0 && max_err > ghost_fac*max_prev) //&& nmb_ghost > 0.25*nmb_pts
      // 	{
      // 	  // Collect element for update of ghost points
      // 	  ghost_elems.push_back(it->second.get());
      // 	}

      // Store updated accuracy information in the element
      it->second->setAccuracyInfo(acc_err, av_err, max_err, outside, acc_outside);
      it->second->setSignedAccInfo(acc_err_pos, nmb_err_pos, acc_err_neg, nmb_err_neg);
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
		  Point tmppt(points.begin()+kh1*del+(dim==3)*3, points.begin()+(kh1+1)*del-1);
		  of << tmppt << std::endl;
		}
	    }
	  //if (nmb_ghost > 0)
	  //  {
	  //    of << "400 1 0 4 50 155 50 255" << std::endl;
	  //    of << nmb_ghost << std::endl;
	  //    for (int kh1=0; kh1<nmb_ghost; ++kh1)
	  //	{
	  //	  Point tmppt(ghost_points.begin()+kh1*del+(dim==3)*2,
	  //		      ghost_points.begin()+(kh1+1)*del-1);
	  //  	  of << tmppt << std::endl;
	  //	}
	  //  }
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
    {
      avdist_ /= (double)outsideeps_;
      avout_ /= (double)outsideeps_;
    }

// #ifdef _OPENMP
//   double time1 = omp_get_wtime();
//   double time_spent = time1 - time0;
//   std::cout << "time_spent in computeAccuracy: " << time_spent << std::endl;
//   std::cout << "time_spent in computeAccuracyElement: " << time_computeAccuracyElement << std::endl;
// #endif
}


//==============================================================================
void LRVolApprox::computeAccuracy_omp(vector<Element3D*>& ghost_elems)
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
#ifdef DEBUG
  std::cout << "Compute accuracy OMP" << std::endl;
#endif
  // Initiate accuracy information
  maxdist_ = 0.0;
  avdist_ = 0.0;
  avdist_all_ = 0.0;
  outsideeps_ = 0;
  maxout_ = 0.0;
  avout_ = 0.0;

  //const Array<double,6> dom = vol_->parameterSpan();
  int dim = vol_->dimension();
  int del = 4 + dim;  // Parameter tripple, position and distance between volume and point
  LRSplineVolume::ElementMap::const_iterator it;
  int kj;

  double ghost_fac = 0.8;
  ghost_elems.clear();

  //for (it=vol_->elementsBegin(), kj=0; it != vol_->elementsEnd(); ++it, ++kj)
  vector<LRSplineVolume::ElementMap::const_iterator> elem_iters;
  const int num_elem = vol_->numElements();
  elem_iters.reserve(num_elem);
  for (LRSplineVolume::ElementMap::const_iterator it=vol_->elementsBegin();
       it != vol_->elementsEnd(); ++it)
  {
      elem_iters.push_back(it);
  }

  vector<double> elemmax(num_elem, 0.0);
  vector<double> elemmax_out(num_elem, 0.0);
  vector<double> elemacc_out(num_elem, 0.0);
  vector<double> elemacc_all(num_elem, 0.0);
  vector<double> elem_avout(num_elem, 0.0);
  vector<int> elemout(num_elem, 0);
#pragma omp parallel default(none) private(kj, it) shared(dim, elem_iters, del, elemmax, elemmax_out, elemacc_out, elemacc_all, elem_avout, elemout, num_elem)
  {
      // double av_prev, max_prev;
      // int nmb_out_prev;
      double umin, umax, vmin, vmax;
      double max_err;
      double av_err;
      double acc_err;
      int outside;
      double acc_outside;
      double acc_err_sgn;
      double av_err_sgn;
      double acc_err_pos;
      double acc_err_neg;
      int nmb_err_pos;
      int nmb_err_neg;
      int nmb_pts;
      size_t kr, nb;
      int ki;
      double *curr;
      double dist, dist2;
      Element3D *elem;
      double acc_prev;
      double tol;

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
	  nmb_pts = it->second->nmbDataPoints();

	  // Local error information
	  max_err = 0.0;
	  av_err = 0.0;
	  acc_err = 0.0;
	  outside = 0;
	  acc_outside = 0.0;
	  acc_err_sgn = 0.0;
	  av_err_sgn = 0.0;
 
	  acc_err_pos = 0.0;
	  acc_err_neg = 0.0;
	  nmb_err_pos = 0;
	  nmb_err_neg = 0;
	  
	  // Check if the accuracy can have been changed
	  const vector<LRBSpline3D*>& bsplines = it->second->getSupport();
	  for (nb=0; nb<bsplines.size(); ++nb)
	      if (!bsplines[nb]->coefFixed())
		  break;

	  if (true /* nb < bsplines.size()*/) // || useMBA_ 
	  {
	      // Compute distances in data points and update parameter pairs
	      // if requested
// #ifdef _OPENMP
// 	    double time0_part = omp_get_wtime();
// #endif
	      if (nmb_pts > 0)
	      {
		  computeAccuracyElement(points, nmb_pts, del, it->second.get());
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
	      dist = curr[del-1];
	      dist2 = fabs(curr[del-1]);
	      elemmax[kj] = std::max(elemmax[kj], dist2);
	      max_err = std::max(max_err, dist2);
	      acc_err += dist2;
	      acc_err_sgn += curr[del-1];
	      elemacc_all[kj] += dist2;
	      tol = aepsge_;
	      for (kr=0; kr<tolerances_.size(); ++kr)
		{
		  if (tolerances_[kr].contains(curr[0], curr[1], curr[2]))
		    {
		      tol = tolerances_[kr].tol;
		      break;
		    }
		}
	      elemmax_out[kj] = std::max(elemmax_out[kj], dist2-tol);
	      if (dist > tol)
		{
		  acc_err_pos += (dist-tol);
		  nmb_err_pos++;
		}
	      else if (dist < -tol)
		{
		  acc_err_neg += ((-dist) - tol);
		  nmb_err_neg++;
		}
	      if (dist2 > tol)
		{
		  av_err_sgn += curr[del-1];
		  elemacc_out[kj] += dist2;
		  elemout[kj]++;
		  av_err += dist2;
		  outside++;
		  elem_avout[kj] += (dist2-tol);
		  acc_outside += (dist2-tol);
		}

	      // if (dim == 3 && repar_)
	      // {
	      // 	  // Check if the point has moved
	      // 	  if (curr[0] < umin || curr[0] > umax || curr[1] < vmin || curr[1] > vmax)
	      // 	  {
	      // 	      // Find element
	      // 	      elem = vol_->coveringElement(curr[0], curr[1]);
	      // 	      elem->addDataPoints(points.begin()+ki*del, 
	      // 				  points.begin()+(ki+1)*del, false);
	      // 	      it->second->eraseDataPoints(points.begin()+ki*del, 
	      // 					  points.begin()+(ki+1)*del);
	      // 	      nmb_pts--;
	      // 	  }
	      // 	  else
	      // 	  {
	      // 	      curr += del;
	      // 	      ki++;
	      // 	  }
	      // }
	      curr += del;
	      ki++;
	  }
	  if (outside > 0)
	  {
	      av_err /= (double)outside;
	      av_err_sgn /= (double)outside;
	  }

	  // Previous accuracy information
	  // acc_prev = it->second->getAccumulatedError();
	  // it->second->getAccuracyInfo(av_prev, max_prev, nmb_out_prev);

// 	  if (max_err > aepsge_ && max_prev > 0.0 && max_err > ghost_fac*max_prev &&
// 	      nmb_ghost > 0.25*nmb_pts)
// 	  {
// 	      // Collect element for update of ghost points
// #pragma omp critical
// 	      ghost_elems.push_back(it->second.get());
// 	  }

	  // Store updated accuracy information in the element
	  it->second->setAccuracyInfo(acc_err, av_err, max_err, outside, acc_outside);
	  it->second->setSignedAccInfo(acc_err_pos, nmb_err_pos, acc_err_neg, nmb_err_neg);
	  //it->second->setAccuracyInfo(acc_err, av_err, max_err, outside);

      }
  }

 for (kj=0; kj<num_elem; ++kj)
   {
     maxdist_ = std::max(maxdist_, elemmax[kj]);
     maxout_ = std::max(maxout_, elemmax_out[kj]);
     avdist_ += elemacc_out[kj];
     avdist_all_ += elemacc_all[kj];
     outsideeps_ += elemout[kj];
     avout_ += elem_avout[kj];
   }
  avdist_all_ /= (double)nmb_pts_;
  if (outsideeps_ > 0)
    {
      avdist_ /= (double)outsideeps_;
      avout_ /= (double)outsideeps_;
    }

// #ifdef _OPENMP
//   double time1 = omp_get_wtime();
//   double time_spent = time1 - time0;
//   std::cout << "time_spent in computeAccuracy: " << time_spent << std::endl;
//   std::cout << "time_spent in computeAccuracyElement: " << time_computeAccuracyElement << std::endl;
// #endif
#ifdef DEBUG
  std::cout << "Finished compute accuracy OMP" << std::endl;
#endif
}


//==============================================================================
  void LRVolApprox::computeAccuracyElement(vector<double>& points, int nmb, int del,
                                           const Element3D* elem)
//==============================================================================
{
  int ki, kj, kk, kr;
  double *curr;
  int dim = vol_->dimension();
  //int maxiter = 3; //4; // @obar WHAT IS THIS??
  Element3D* elem3 = (Element3D*)elem;

  // Fetch basis functions
  const vector<LRBSpline3D*>& bsplines = elem->getSupport();
  const int nmb_bsplines = (int)bsplines.size();
  double volval;
  vector<double> bval(nmb_bsplines);
  vector<double> tmpval(3*nmb_bsplines);

  //int idx1, idx2, idx3, sgn;
  double dist;
  Point curr_pt;
  //const int num_threads = 8;
  //const int dyn_div = nmb/num_threads;

  double tol = 1.0e-12;  // Numeric tolerance

  double umax = vol_->endparam_u();
  double vmax = vol_->endparam_v();
  double wmax = vol_->endparam_w();

  for (ki=0, curr=&points[0]; ki<nmb; ++ki, curr+=del)
    {
      // Evaluate
      if (dim == 1)
	{
	  bool u_at_end = (curr[0] > umax-tol) ? true : false;
	  bool v_at_end = (curr[1] > vmax-tol) ? true : false;
          bool w_at_end = (curr[2] > wmax-tol) ? true : false;
	  LRSpline3DUtils::evalAllBSplines2(bsplines, curr[0], curr[1], curr[2],
					    u_at_end, v_at_end, w_at_end,
					    &bval[0], &tmpval[0]);
	  // vector<Point> bpos;
	  // LRSpline3DUtils::evalAllBSplinePos(bsplines, curr[0], curr[1],
	  // 				     curr[2], u_at_end, v_at_end,
	  // 				     w_at_end, bpos);
	  // vol_->point(pos, curr[0], curr[1], curr[2]);
	  volval = 0.0;
	  //double volval2 = 0.0;
	  // double bval;
	  for (kr=0; kr<nmb_bsplines; ++kr)
	    {
	      // bsplines[kr]->evalpos(curr[0], curr[1], curr[2], &bval);
	      // volval2 += bval;
	      const Point& tmp_pt = bsplines[kr]->coefTimesGamma();
	      volval += bval[kr]*tmp_pt[0];
	      // volval += bpos[kr][0];
	    }
	  // if (fabs(volval-volval2) > 1.0e-6)
	  //   std::cout << "Value mismatch " << volval << " " << volval2 << std::endl;
	  // if (fabs(pos[0]-volval) > 1.0e-4)
	  //   std::cout << "Evaluation mismatch" << std::endl;
	  dist = curr[3] - volval;
	  //dist = curr[2] - pos[0];
	}
      else
	{
	  MESSAGE("3D NOT SUPPORTED YET...");
	}
      curr[del-1] = dist;
      // curr_pt = Point(curr+(dim==3)*3, curr+del-1);
      // if (check_close_ && dim == 3)
      // 	{
      // 	  MESSAGE("3D NOT SUPPORTED YET...");
      // 	  // Compute closest point
      // 	  // VSK. 052013. Should be changed to make use of the fact that we know
      // 	  // the element (at least initially)
      // 	  // double upar, vpar;
      // 	  // Point close_pt;
      // 	  //vol_->setCurrentElement((Element3D*)elem);
      // 	  /*vol_->closestPoint(curr_pt, upar, vpar, wpar, close_pt,
      // 			     dist, aepsge_);//, maxiter, elem3, &, curr);
      // 	  vec = curr_pt - close_pt;
      // 	  // Point norm;
      // 	  vol_->normal(norm, upar, vpar, wpar);
      // 	  if (vec*norm < 0.0)
      // 	    dist *= -1;
      // 	  if (to3D_ >= 0)
      // 	    dist = curr[del-2]-close_pt[2];*/
      // 	}

    }
 }

//==============================================================================
void LRVolApprox::computeAccuracyElement_omp(vector<double>& points, int nmb, int del,
                                             const Element3D* elem)
//==============================================================================
{
  int ki, kj, kr;
  double *curr;
  int dim = vol_->dimension();
  int maxiter = 3; //4;
  Element3D* elem2 = (Element3D*)elem;

  // Fetch basis functions
  const vector<LRBSpline3D*>& bsplines = elem->getSupport();
  const int nmb_bsplines = (int)bsplines.size();
  double volval;

  double umax = vol_->endparam_u();
  double vmax = vol_->endparam_v();
  double wmax = vol_->endparam_w();

  double tol = 1.0e-12;  // Numeric tolerance
  int sgn;
  double dist, upar, vpar, wpar;
  Point close_pt, vec, norm, pos, curr_pt;
  bool u_at_end, v_at_end, w_at_end;
  const int num_threads = 8;
  const int dyn_div = nmb/num_threads;
  vector<double> bval;//(nmb_bsplines);
  vector<double> tmpval;//(3*nmb_bsplines);

#ifdef _OPENMP
  pthread_attr_t attr;
  size_t stacksize;
  pthread_attr_getstacksize(&attr, &stacksize);
#endif
  //	std::cout << "stacksize (in MB): " << (double)stacksize/(1024.0*1024.0) << std::endl;
  //	omp_set_num_threads(4);
#pragma omp parallel default(none) private(ki, curr, dist, u_at_end, v_at_end, w_at_end, volval, kr, kj, bval, tmpval) \
  shared(points, nmb, del, dim, umax, vmax, wmax, tol, maxiter, elem2, bsplines, nmb_bsplines)
  {
    bval.resize(bsplines.size());
    tmpval.resize(3*bsplines.size());
#pragma omp for schedule(auto)//static, 4)//runtime)//guided)//auto)
  //#pragma omp for schedule(dynamic, 4)//static, 4)//runtime)//guided)//auto)
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
      //curr_pt = Point(curr+(dim==3)*2, curr+del-1);
      if (check_close_ && dim == 3)
	{
	  // Compute closest point
	  // VSK. 052013. Should be changed to make use of the fact that we know
	  // the element (at least initially)
	  // double upar, vpar;
	  // Point close_pt;
	  // vol_->closestPoint(curr_pt, upar, vpar, close_pt,
	  // 		     dist, aepsge_, maxiter, elem2, &rd, curr);
	  // vec = curr_pt - close_pt;
	  // // Point norm;
	  // vol_->normal(norm, upar, vpar, elem2);
	  // if (vec*norm < 0.0)
	  //   dist *= -1;
	  // if (to3D_ >= 0)
	  //   dist = curr[del-2]-close_pt[2];
	  MESSAGE("3D NOT SUPPORTED YET...");
	}
      else
	{
	  // Evaluate
	  if (dim == 1)
	    {
	      u_at_end = (curr[0] > umax-tol) ? true : false;
	      v_at_end = (curr[1] > vmax-tol) ? true : false;
	      w_at_end = (curr[2] > wmax-tol) ? true : false;

	      LRSpline3DUtils::evalAllBSplines2(bsplines, curr[0], curr[1], curr[2],
						u_at_end, v_at_end, w_at_end,
						&bval[0], &tmpval[0]);
	      volval = 0.0;
	      for (kr=0; kr<nmb_bsplines; ++kr)
		{
		  const Point& tmp_pt = bsplines[kr]->coefTimesGamma();
		  volval += bval[kr]*tmp_pt[0];
		}
	      
	      dist = curr[3] - volval;
	      //dist = curr[2] - pos[0];
	    }
	  else
	    {
	      MESSAGE("3D NOT SUPPORTED YET...");
	      // Point pos;
	      // vol_->point(pos, curr[0], curr[1], elem2);
	      // dist = pos.dist(Point(curr+2, curr+del));
	      // vec = curr_pt - pos;
	      // // Point norm;
	      // vol_->normal(norm, curr[0], curr[1]);
	      // if (vec*norm < 0.0)
	      //   dist *= -1;
	    }
	}
    }
  curr[del-1] = dist;
  }
}

 
 
//==============================================================================
bool compare_elems(pair<Element3D*,double> el1, pair<Element3D*,double> el2)
{
  return (el1.second > el2.second);
}

int divide(double *err, int *perm, int low, int high)
{
  int p1 = perm[high];
  int ki = low - 1;
  for (int kj=low; kj<=high-1; ++kj)
    {
      if (err[perm[kj]] > err[p1])
	{
	  ++ki;
	  std::swap(perm[ki],perm[kj]);
	}
    }
  std::swap(perm[ki+1], perm[high]);
  return (ki+1);
}

void quicksort(double *err, int *perm, int low, int high)
{
  if (low < high)
    {
      int pos = divide(err, perm, low, high);

      quicksort(err, perm, low, pos-1);
      quicksort(err, perm, pos+1, high);
    }
}

//==============================================================================
int LRVolApprox::refineVol(double threshold)
//==============================================================================
{
  ref_x_ = ref_y_ = ref_z_ = 0;
  nmb1_ = nmb2_ = nmb3_ = 0;
  
  int choice = 0;  // Strategy for knot insertion in one single B-spline
  int totdeg = vol_->degree(XDIR)*vol_->degree(YDIR)*vol_->degree(ZDIR);

#ifdef DEBUG
  std::ofstream of("error_elems.txt");
  double maxerrfac = 0.8;
#endif

  double mineps = std::min(aepsge_, threshold);
  double error_fac = 0.1;
  double error_fac2 = 10.0;
  int el_out = 0;
  vector<pair<Element3D*, double> > elem_out;
  double av_wgt = 0.0;
  double nmb_frac = 0.1;
  int num_elem = vol_->numElements();
  int nmb_no_elem = 0;
  double av_nmb = (double)nmb_pts_/(double)num_elem;
  for (LRSplineVolume::ElementMap::const_iterator it=vol_->elementsBegin();
       it != vol_->elementsEnd(); ++it)
    {
      double av_err, max_err;
      int nmb_out;
      int nmb_pts = it->second->nmbDataPoints();
      if (nmb_pts == 0)
	nmb_no_elem++;
      it->second->getAccuracyInfo(av_err, max_err, nmb_out);
      if (nmb_out > 0)
	{
#ifdef DEBUG
	  if (max_err > maxerrfac*maxdist_ && max_err > 0.5*(aepsge_ + maxdist_))
	    {
	      of << it->second.get() << ", " << max_err << ", " << nmb_pts << ", " << nmb_out << ", " << av_err << std::endl;
	      of << it->second->umin() << ", " << it->second->umax() << ", ";
	      of << it->second->vmin() << ", " << it->second->vmax() << ", ";
	      of << it->second->wmin() << ", " << it->second->wmax() << std::endl;
	      of << std::endl;
	    }
#endif
	  double wgt = nmb_out + max_err + av_err;
	  // if ((double)nmb_pts < nmb_frac*av_nmb)//20)
	  //   wgt /= 2.0;

	  if ((double)nmb_out/(double)nmb_pts < outfrac_)
	    wgt *= error_fac;
	  av_wgt += wgt;
	  elem_out.push_back(std::make_pair(it->second.get(), wgt));
	  el_out++;
	}
    }
  av_wgt /= (double)el_out;
  double wgt_fac = 0.5;

  vector<LRSplineVolume::Refinement3D> refs_x, refs_y, refs_z;
  double elem_frac = 0.1; //0.005; //0.01;
   size_t kr = 0;
#ifdef DEBUG0
  std::cout << "Number of elements: " << num_elem << std::endl;
  std::cout << "Number of elements without points: " << nmb_no_elem << std::endl;
  std::cout << "Fraction of outside elements: " << (double)elem_out.size()/(double)num_elem << std::endl;
#endif

  std::set<Element3D*> elems;
  int num_bspl = vol_->numBasisFunctions();
  int nmb_refs = 0;
  double split_fac = 0.1;
  int level_bsplines = (int)(split_fac*num_bspl);
  if ((double)elem_out.size() > elem_frac*(double)num_elem)
    {
  int group_fac = 3;
  wgt_fac = 0.75;
#ifdef DEBUG0
  std::ofstream ofout("outfrac.txt");
  ofout << "Number of B-splines: " << num_bspl << std::endl;
#endif
  
  // Construct indexed bspline array and collect related accuracy information
  vector<LRBSpline3D*> bsplines(num_bspl);
  vector<double> error(num_bspl, 0.0);
  vector<double> max_error(num_bspl, 0.0);
  vector<double> av_error(num_bspl, 0.0);
  vector<double> acc_err_pos(num_bspl, 0.0);
  vector<double> acc_err_neg(num_bspl, 0.0);
  vector<int> num_pts(num_bspl, 0);
  vector<int> num_out_pts(num_bspl, 0); 
  vector<int> nmb_out_pos(num_bspl, 0); 
  vector<int> nmb_out_neg(num_bspl, 0); 
  vector<double> error2(num_bspl);
  double mean_err = 0.0;
  double domain_size;
  double average_nmb_out = 0.0;
  double average_nmb = 0.0;
  double basis_average_out = 0.0;
  int nmb_basis_out = 0;
  vector<int> bspl_perm(num_bspl, 0);
  int nmb_perm = 0;
  double frac_av_out = 0.0;
  double frac_nmb_out = 0.0;
  kr = 0;
  for (LRSplineVolume::BSplineMap::const_iterator it=vol_->basisFunctionsBegin();
       it != vol_->basisFunctionsEnd(); ++it, ++kr)
    {
      LRBSpline3D* curr = it->second.get();

      int el_out2 = 0.0;
      for (auto it2=curr->supportedElementBegin(); 
	   it2 != curr->supportedElementEnd(); ++it2)
	{
	  num_pts[kr] += (*it2)->nmbDataPoints();
	  num_out_pts[kr] += (*it2)->getNmbOutsideTol();
	  //error[kr] += (*it2)->getAccumulatedError();
	  error[kr] += (*it2)->getAccumulatedOutside();
	  max_error[kr] = std::max(max_error[kr], (*it2)->getMaxError());
	  av_error[kr] += (*it2)->getAverageError();  // Only counting those 
	  // points being outside of the tolerance
	  if ((*it2)->getNmbOutsideTol() > 0)
	    el_out2++;

	  double curr_acc_pos, curr_acc_neg;
	  int curr_nmb_pos, curr_nmb_neg;
	  (*it2)->getSignedAccInfo(curr_acc_pos, curr_nmb_pos,
				 curr_acc_neg, curr_nmb_neg);
	  acc_err_pos[kr] += curr_acc_pos;
	  acc_err_neg[kr] += curr_acc_neg;
	  nmb_out_pos[kr] += curr_nmb_pos;
	  nmb_out_neg[kr] += curr_nmb_neg;
	}
     if (el_out2 > 0)
	++nmb_basis_out;
      
      av_error[kr] /= (double)(curr->nmbSupportedElements());

      int bnmb_el = curr->nmbSupportedElements();
      basis_average_out += (double)el_out2/(double)bnmb_el;
      if (nmb_out_pos[kr] > 0 || nmb_out_neg[kr] > 0)
	{
	  frac_av_out += std::min(acc_err_pos[kr], acc_err_neg[kr])/
	    std::max(acc_err_pos[kr], acc_err_neg[kr]);
	  frac_nmb_out += std::min((int)nmb_out_pos[kr], (int)nmb_out_neg[kr])/
	    std::max((double)nmb_out_pos[kr], (double)nmb_out_neg[kr]);
	}
      
     // Use cubic root to reduce the significance of this property compared to the
      // error
      domain_size = std::pow((curr->umax()-curr->umin())
                                 *(curr->vmax()-curr->vmin())
                                 *(curr->wmax()-curr->wmin()),1/3.0);
      error2[kr] = error[kr]*domain_size;
      if (num_out_pts[kr] > group_fac /* || (double)num_out_pts[kr] > 
					 error_fac*((double)num_pts[kr])*/)
      	error2[kr] *= error_fac2;

      if ((double)num_out_pts[kr]/(double)num_pts[kr] < outfrac_)
	{
	  // Neglectable fraction of outside points
	  error2[kr] *= error_fac;
	}
      
      bsplines[kr] = curr;
      mean_err += av_error[kr];
      average_nmb_out += (double)(num_out_pts[kr]);
      average_nmb += (double)(num_pts[kr]);

      if (num_out_pts[kr] > 0 && (double)num_out_pts[kr]/(double)num_pts[kr] >= outfrac_)
	bspl_perm[nmb_perm++] = kr;
    }
  mean_err /= (double)num_bspl;
  average_nmb_out /= (double)num_bspl;
  average_nmb /= (double)num_bspl;
  basis_average_out /= (double)num_bspl;
  frac_av_out /= (double)num_bspl;
  frac_nmb_out /= (double)num_bspl;

#ifdef DEBUG0
  std::cout << "Number of Bsplines: " << num_bspl << std::endl;
  std::cout << "Number of selected Bsplines: " << nmb_perm << std::endl;
  std::cout << "Average number of points: " << average_nmb << std::endl;
  std::cout << "Average number of outside points: " << average_nmb_out << std::endl;
  std::cout << "Average fraction of outside elements in bspline: " << basis_average_out << std::endl;
  std::cout << "Elements with outside points: " << elem_out.size() << std::endl;
  std::cout << "Average distribution of outside distance (above - below): " << frac_av_out << std::endl;
  std::cout << "Average distribution of number of outside points (above - below): " << frac_nmb_out << std::endl;
#endif
 // Sort bsplines according to average error weighted with the domain size
  // int ki, kj;
  // for (ki=0; ki<num_bspl; ++ki)
  //   bspl_perm[ki] = ki;
#ifdef DEBUG0
  std::cout << "Before sorting B-splines " << std::endl;
#endif
  // Do the sorting
  quicksort(&error2[0], &bspl_perm[0], 0, nmb_perm-1);
  // for (ki=0; ki<num_bspl; ++ki)
  //   for (kj=ki+1; kj<num_bspl; ++kj)
  //     {
  // 	// double curr_err1 = error[bspl_perm[ki]]*domain_size[bspl_perm[ki]];
  // 	// double curr_err2 = error[bspl_perm[kj]]*domain_size[bspl_perm[kj]];

  // 	// if (num_out_pts[bspl_perm[ki]] > group_fac ||
  // 	//     (double)num_out_pts[bspl_perm[ki]] > error_fac*((double)num_pts[bspl_perm[ki]]) )
  // 	//   curr_err1 *= error_fac2;
  // 	// if (num_out_pts[bspl_perm[kj]] > group_fac ||
  // 	//     (double)num_out_pts[bspl_perm[kj]] > error_fac*((double)num_pts[bspl_perm[kj]]) )
  // 	//   curr_err2 *= error_fac2;

  // 	// Modify if there is a significant number of large error points 
  // 	//if (curr_err1 < curr_err2)
  // 	if (error2[bspl_perm[ki]] < error2[bspl_perm[kj]])
  // 	  std::swap(bspl_perm[ki], bspl_perm[kj]);
  //     }
#ifdef DEBUG0
  std::ofstream of("sorted_bsplines.txt");
  for (int kaa=0; kaa<nmb_perm; ++kaa)
    {
      of << bsplines[bspl_perm[kaa]] << " (" << bsplines[bspl_perm[kaa]]->umin() << ",";
      of << bsplines[bspl_perm[kaa]]->umax() << ") x (" << bsplines[bspl_perm[kaa]]->vmin();
      of << "," << bsplines[bspl_perm[kaa]]->vmax() << ") x (" << bsplines[bspl_perm[kaa]]->wmin();
      of << "," << bsplines[bspl_perm[kaa]]->wmax() << "): " << error2[bspl_perm[kaa]] << std::endl;
      of << num_pts[bspl_perm[kaa]] << " " << num_out_pts[bspl_perm[kaa]] << " ";
      of << error[bspl_perm[kaa]] << " " << max_error[bspl_perm[kaa]] << " ";
      of << av_error[bspl_perm[kaa]] << std::endl;
    }
  
  std::cout << "After sorting B-splines " << std::endl;
#endif
  
  // Split the most important B-splines, but only if the maximum
  // error is larger than the tolerance
  // We only split half the BSplines? Reduction factor...
  double red_fac = 0.5; //0.75; //0.6;
  int nmb_split = (int)(red_fac*nmb_perm);
  nmb_split = std::min(nmb_perm, std::max(nmb_split, 5*totdeg));
  int min_nmb_pts = 4;//2;//1; //4;

  // if (choice == 0) { // direct generalisation of 2D version

  int nmb_fixed = 0;
  double average_threshold = std::max(0.01*average_nmb, average_nmb_out);
  //for (int ki=0; ki<nmb_perm; ++ki)
  for (int ki=0; ki<nmb_split; ++ki)
    {
      //if (max_error[bspl_perm[kr]] < aepsge_)
      if (num_out_pts[bspl_perm[ki]] == 0)
	{
	  // TEST. Keep coefficient fixed
	  bsplines[bspl_perm[kr]]->setFixCoef(1);
	  nmb_fixed++;
	  continue;
	}
      else
	bsplines[bspl_perm[ki]]->setFixCoef(0);  // Adjacent B-splines may 
      // have changed

      // Do not split any B-spline with too few points in its domain
      if (num_pts[bspl_perm[ki]] < min_nmb_pts)
	continue;

      // if (max_error[bspl_perm[ki]] < ((nmb_refs < level_bsplines) ?
      // 				      mineps : threshold))
      if (threshold > aepsge_ && max_error[bspl_perm[ki]] < threshold)
      	continue;
	  
      // if (nmb_refs >= nmb_split)
      // 	break;
      double av_pos = (nmb_out_pos[bspl_perm[ki]] > 0) ?
	acc_err_pos[bspl_perm[ki]]/(double)nmb_out_pos[bspl_perm[ki]] : 0.0;
      double av_neg = (nmb_out_neg[bspl_perm[ki]] > 0) ?
	acc_err_neg[bspl_perm[ki]]/(double)nmb_out_neg[bspl_perm[ki]] : 0.0;
#ifdef DEBUG0
      ofout << acc_err_pos[bspl_perm[ki]] << " " << nmb_out_pos[bspl_perm[ki]] << " " << av_pos << " ";
       ofout << acc_err_neg[bspl_perm[ki]] << " " << nmb_out_neg[bspl_perm[ki]] << " " << av_neg  << std::endl;
#endif

      // How to split
      vector<Element3D*> elem_div;
      defineRefs(bsplines[bspl_perm[ki]], average_threshold,
		 refs_x, refs_y, refs_z, elem_div);
      if (elem_div.size() > 0)
	nmb_refs++;  //  B-spline is split
	  
      elems.insert(elem_div.begin(), elem_div.end());
      //defineRefs(bsplines[bspl_perm[kr]], refs, choice);
    }
#ifdef DEBUG0
  std::cout << "nmb_refs: " << nmb_refs << ", nmb_split: " << nmb_split;
  std::cout << ", nmb basis out: " << nmb_basis_out << std::endl;
  std::cout << "nmb1 = " << nmb1_ << ", nmb2 = " << nmb2_;
  std::cout << ", nmb3 = " << nmb3_ << std::endl;
  nmb1_ = nmb2_ = nmb3_ = 0;
#endif
    }
#ifdef DEBUG0
 std::cout << "Removing affected elements. elems.size(): " << elems.size();
 std::cout << ", elem_out.size(): " << elem_out.size() << std::endl;
#endif
 std::vector<Element3D*> elems2(elems.begin(), elems.end());
 for (kr=0; kr<elems2.size(); ++kr)
   {
     size_t kh;
     for (kh=0; kh<elem_out.size(); ++kh)
       if (elems2[kr] == elem_out[kh].first)
	 break;
     if (kh < elem_out.size())
       elem_out.erase(elem_out.begin() + kh);
   }
#ifdef DEBUG0
 std::cout << "End removing affected elements. elem_out.size(): " << elem_out.size() << std::endl;
  
 std::cout << "Number of refinements: " << refs_x.size()+refs_y.size()+refs_z.size() << std::endl;
 std::cout << "Remaining elements with outside points: " << elem_out.size() << std::endl;
#endif
      // Sort remaining elements
      double frac = wgt_fac*av_wgt*(double)elem_out.size()/(double)num_elem; //0.6*av_wgt;
#ifdef DEBUG0
      std::cout << "frac = " << frac << std::endl;
#endif
      std::sort(elem_out.begin(), elem_out.end(), compare_elems);
      if (elem_out.size() > 0)
	{
	  double maxeldist = 0.0;
	  double maxdistwgt;
	  for (kr=0; kr<elem_out.size(); ++kr)
	    {
	      double elmax = elem_out[kr].first->getMaxError();
	      if (elmax > maxeldist)
		{
		  maxeldist = elmax;
		  maxdistwgt = elem_out[kr].second;
		}
	    }
#ifdef DEBUG0
	  std::cout << "minwgt: " << elem_out[elem_out.size()-1].second << ", maxwgt: ";
	  std::cout << elem_out[0].second << ", average: " << av_wgt << std::endl;
	  std::cout << "median wgt: " << elem_out[elem_out.size()/2].second;
	  std::cout << ", maximum dist: " << maxeldist << ", wgt " << maxdistwgt << std::endl;
#endif
	}
      double maxeldist = 0.0;
      int nmb_dismissed = 0;
      int nmb_refinements = nmb1_ + nmb2_ + nmb3_;
      for (kr=0; kr<elem_out.size(); ++kr)
	{
	  if (elem_out[kr].second < wgt_fac*av_wgt) //frac) 
	    break;   // Not a significant element
	  // if (elem_out[kr].first->getMaxError() < ((nmb_refs < level_bsplines) ? aepsge_ : threshold))//threshold)
	  if (threshold > aepsge_ && nmb_refs >= level_bsplines &&
	      elem_out[kr].first->getMaxError() < threshold)
	    {
	      nmb_dismissed++;
	      maxeldist = std::max(elem_out[kr].first->getMaxError(), maxeldist);
	      continue;
	    }
	  
	  vector<Element3D*> elements;  // Elements affected by the refinement(s)
	  checkFeasibleRef(elem_out[kr].first, refs_x, refs_y, refs_z, elements);
	  if (nmb1_ + nmb2_ + nmb3_ > nmb_refinements) //elements.size() > 0)
	    {
	      nmb_refs++;
	      nmb_refinements = nmb1_ + nmb2_ + nmb3_;
	    }
	  else
	    {
	      nmb_dismissed++;
	      //maxeldist = std::max(elem_out[kr].first->getMaxError(), maxeldist);
	    }
	  for (size_t ki=0; ki<elements.size(); ++ki)
	    {
	      size_t kj;
	      for (kj=kr+1; kj<elem_out.size(); ++kj)
		if (elements[ki] == elem_out[kj].first)
		  break;
	      if (kj < elem_out.size())
		elem_out.erase(elem_out.begin()+kj);
	    }
	}
      // nmb_dismissed += ((int)elem_out.size() - kr);
      // for (kr=(kr>0) ? kr-1 : 0; kr<elem_out.size(); ++kr)
      // 	{
      // 	  maxeldist = std::max(elem_out[kr].first->getMaxError(), maxeldist);
      // 	}
      // std::cout << "Remaining maximum dist: " << maxeldist;
#ifdef DEBUG0
      std::cout << "Number of elements not split: " << nmb_dismissed << std::endl;

      std::cout << "nmb_refs: " << nmb_refs << std::endl;
 std::cout << "Number of refinements: " << refs_x.size()+refs_y.size()+refs_z.size() << std::endl;
#endif
      // // Flag coefficient if not necessary to update
     //  for (; kr<bspl_perm.size(); ++kr)
     // 	{
     // 	  //if (max_error[bspl_perm[kr]] < aepsge_)
     // 	  if (num_out_pts[bspl_perm[kr]] == 0)
     // 	    {
     // 	      // TEST. Keep coefficient fixed
     // 	      bsplines[bspl_perm[kr]]->setFixCoef(1);
     // 	      nmb_fixed++;
     // 	    }
	// }
   // } else { // Tensor product refinement
   //    const Mesh3D& mesh = vol_->mesh();
   //    double umin = vol_->paramMin(XDIR);
   //    double vmin = vol_->paramMin(YDIR);
   //    double wmin = vol_->paramMin(ZDIR);
   //    double umax = vol_->paramMax(XDIR);
   //    double vmax = vol_->paramMax(YDIR);
   //    double wmax = vol_->paramMax(ZDIR);
   //    LRSplineVolume::Refinement3D ref;
   //    for (int ix=0; ix!=mesh.numDistinctKnots(XDIR)-1; ++ix)
   //      {
   //        ref.setVal((mesh.kval(XDIR,ix)+mesh.kval(XDIR,ix+1))*0.5, vmin,vmax,wmin,wmax,XDIR,1);
   //        refs.push_back(ref);
   //      }
   //    for (int ix=0; ix!=mesh.numDistinctKnots(YDIR)-1; ++ix)
   //      {
   //        ref.setVal((mesh.kval(YDIR,ix)+mesh.kval(YDIR,ix+1))*0.5, wmin,wmax,umin,umax,YDIR,1);
   //        refs.push_back(ref);
   //      }
   //    for (int ix=0; ix!=mesh.numDistinctKnots(ZDIR)-1; ++ix)
   //      {
   //        ref.setVal((mesh.kval(ZDIR,ix)+mesh.kval(ZDIR,ix+1))*0.5, umin,umax,vmin,vmax,ZDIR,1);
   //        refs.push_back(ref);
   //      }
   //  }
  // for (kr=0; kr<refs.size(); ++kr)
  // {
  //   vol_->refine(refs[kr], true );
  // }
#ifdef DEBUG0
 std::cout << "ref_x = " << ref_x_ << ", ref_y = " << ref_y_;
 std::cout << ", ref_z = " << ref_z_ << std::endl;
 std::cout << "nmb1 = " << nmb1_ << ", nmb2 = " << nmb2_;
 std::cout << ", nmb3 = " << nmb3_ << std::endl;
 std::cout << "Refine x " << std::endl;
#endif
 vol_->refine(refs_x, true);
#ifdef DEBUG0
 std::cout << "Refine y " << std::endl;
#endif
 vol_->refine(refs_y, true);
#ifdef DEBUG0
 std::cout << "Refine z " << std::endl;
#endif
 vol_->refine(refs_z, true);

#ifdef DEBUG0
  std::cout << "Refs x: " << refs_x.size() << ", y: " << refs_y.size();
  std::cout << ", z: " << refs_z.size() << std::endl;
#endif
  return (int)refs_x.size() + refs_y.size() + refs_z.size();
}
#if 1 
//==============================================================================
void LRVolApprox::checkFeasibleRef(Element3D* elem, 
				   vector<LRSplineVolume::Refinement3D>& refs_x,
				   vector<LRSplineVolume::Refinement3D>& refs_y,
				   vector<LRSplineVolume::Refinement3D>& refs_z,
				   vector<Element3D*>& affected)
//==============================================================================
{
  double tol = vol_->getKnotTol();
  double eps = 0.01;

  // Fetch B-splines
  const vector<LRBSpline3D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();

  int degree1 = vol_->degree(XDIR);
  int degree2 = vol_->degree(YDIR);
  int degree3 = vol_->degree(ZDIR);
  int xmult = (degree1 <= 3) ? 1 : 2;
  int ymult = (degree2 <= 3) ? 1 : 2;
  int zmult = (degree3 <= 3) ? 1 : 2;
  
 // Refine one B-spline with support in the parent element in one or two
  // parameter directions depending on how many elements with a large error
  // that lies in its support
  // First check refinement in the u-direction
  double u_par = 0.5*(elem->umin() + elem->umax());
  double udel = elem->umax() - elem->umin();
  double minsize_u = (usize_min_ > 0.0) ? 2.0*usize_min_ : 1.0e-8;
  size_t ki, kj;
  double udelmax = 0.0;
  bool udir = true;
  std::set<Element3D*> uelems;
  for (ki=0; ki<nmb; ++ki)
    {
      udelmax = std::max(udelmax, bsplines[ki]->umax() - bsplines[ki]->umin());
      
      // Collect elements
      const vector<Element3D*>& curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->umax() < u_par || curr_el[kj]->umin() > u_par)
	    continue;  // Element not affected

	  if (curr_el[kj]->umax() - curr_el[kj]->umin() < minsize_u)
	    {
	      // Element too small to be divided
	      udir = false;
	    }

	  uelems.insert(curr_el[kj]);
	}
    }
	  
  // The v-direction
  double v_par = 0.5*(elem->vmin() + elem->vmax());
  double vdel = elem->vmax() - elem->vmin();
  double minsize_v = (vsize_min_ > 0.0) ? 2.0*vsize_min_ : 1.0e-8;
  double vdelmax = 0.0;
  bool vdir = true;
  std::set<Element3D*> velems;
  for (ki=0; ki<nmb; ++ki)
    {
       vdelmax = std::max(vdelmax, bsplines[ki]->vmax() - bsplines[ki]->vmin());
      
      // Collect elements
      const vector<Element3D*>& curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->vmax() < v_par || curr_el[kj]->vmin() > v_par)
	    continue;  // Element not affected

	  if (curr_el[kj]->vmax() - curr_el[kj]->vmin() < minsize_v)
	    {
	      // Element too small to be divided
	      vdir = false;
	    }
	  
	  velems.insert(curr_el[kj]);
	}
    }
  
   // The w-direction
  double w_par = 0.5*(elem->wmin() + elem->wmax());
  double wdel = elem->wmax() - elem->wmin();
  double minsize_w = (wsize_min_ > 0.0) ? 2.0*wsize_min_ : 1.0e-8;
  double wdelmax = 0.0;
  bool wdir = true;
  std::set<Element3D*> welems;
  for (ki=0; ki<nmb; ++ki)
    {
      wdelmax = std::max(wdelmax, bsplines[ki]->wmax() - bsplines[ki]->wmin());
      
      // Collect elements
      const vector<Element3D*>& curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->wmax() < w_par || curr_el[kj]->wmin() > w_par)
	    continue;  // Element not affected

	  if (curr_el[kj]->wmax() - curr_el[kj]->wmin() < minsize_w)
	    {
	      // Element too small to be divided
	      wdir = false;
	    }
	  
	  welems.insert(curr_el[kj]);
	}
    }
  
  // Estimate significant of split in all parameter directions
  vector<Element3D*> element_u(uelems.begin(), uelems.end());
  vector<Element3D*> element_v(velems.begin(), velems.end());
  vector<Element3D*> element_w(welems.begin(), welems.end());
  double fac = 0.1;
  double fac3 = 0.75; //0.95;
  int nmb_u = 0;
  for (ki=0; ki<element_u.size(); ++ki)
    {
      double max_err, av_err;
      int nmb_outside;
      element_u[ki]->getAccuracyInfo(av_err, max_err, nmb_outside);
      int nmb_pts = element_u[ki]->nmbDataPoints();
      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	nmb_u++;
    }

  int nmb_v = 0;
  for (ki=0; ki<element_v.size(); ++ki)
    {
      double max_err, av_err;
      int nmb_outside;
      element_v[ki]->getAccuracyInfo(av_err, max_err, nmb_outside);
      int nmb_pts = element_v[ki]->nmbDataPoints();
      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	nmb_v++;
    }

  int nmb_w = 0;
  for (ki=0; ki<element_w.size(); ++ki)
    {
      double max_err, av_err;
      int nmb_outside;
      element_w[ki]->getAccuracyInfo(av_err, max_err, nmb_outside);
      int nmb_pts = element_w[ki]->nmbDataPoints();
      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	nmb_w++;
    }

  if (udir == false && vdir == false && wdir == false)
    udir = vdir = wdir = true;   // No direction without small elements,
  // use other criteria
  
  // Assemble information
  int div_x = 0, div_y = 0, div_z = 0;
  double fac2 = 0.5;
  double sizefac = 3.0;
  if (udir &&
      //((nmb_u < 1 && nmb_u >= std::max(nmb_v, nmb_w)) ||
       ((nmb_u <= 1 && nmb_u >= std::max(nmb_v, nmb_w)) ||
      //((nmb_u > std::max(nmb_v, nmb_w)) ||
       udel > sizefac*std::min(vdel, wdel) ||
       udelmax >= std::max(vdelmax, wdelmax)))
    {
      for (ki=0; ki<nmb; ++ki)
	{
	  LRSplineVolume::Refinement3D curr_ref;
	  curr_ref.setVal(u_par, bsplines[ki]->vmin(), bsplines[ki]->vmax(),
			  bsplines[ki]->wmin(), bsplines[ki]->wmax(), XDIR, xmult);
	  appendRef(refs_x, curr_ref, tol);
	  ref_x_++;
	  div_x = 1;
	}
    }
			       
  if (vdir &&
      //((nmb_v < 1 && nmb_v >= std::max(nmb_u, nmb_w)) ||
       ((nmb_v <= 1 && nmb_v >= std::max(nmb_u, nmb_w)) ||
      //((nmb_v > std::max(nmb_u, nmb_w)) ||
       vdel > sizefac*std::min(udel, wdel) ||
       vdelmax >= std::max(udelmax, wdelmax)))
    {
      for (ki=0; ki<nmb; ++ki)
	{
	  LRSplineVolume::Refinement3D curr_ref;
	  curr_ref.setVal(v_par, bsplines[ki]->wmin(), bsplines[ki]->wmax(),
			  bsplines[ki]->umin(), bsplines[ki]->umax(), YDIR, ymult);
	  appendRef(refs_y, curr_ref, tol);
	  ref_y_++;
	  div_y = 1;
	}
    }

  if (wdir &&
      //((nmb_w < 1 && nmb_w >= std::max(nmb_u, nmb_v)) ||
       ((nmb_w <= 1 && nmb_w >= std::max(nmb_u, nmb_v)) ||
      //((nmb_w > std::max(nmb_u, nmb_v)) ||
       wdel > sizefac*std::min(udel, vdel) ||
       wdelmax >= std::max(udelmax, vdelmax)))
    {
     for (ki=0; ki<nmb; ++ki)
	{
	  LRSplineVolume::Refinement3D curr_ref;
	  curr_ref.setVal(w_par, bsplines[ki]->umin(), bsplines[ki]->umax(),
			  bsplines[ki]->vmin(), bsplines[ki]->vmax(), ZDIR, zmult);
	  appendRef(refs_z, curr_ref, tol);
	  ref_z_++;
	  div_z = 1;
	}
    }

  int nmb_div = div_x + div_y + div_z;
  if (nmb_div == 1)
    nmb1_++;
  else if (nmb_div == 2)
    nmb2_++;
  else if (nmb_div == 3)
    nmb3_++;
}
#endif
#if 0
//==============================================================================
void LRVolApprox::checkFeasibleRef(Element3D* elem, 
				   vector<LRSplineVolume::Refinement3D>& refs_x,
				   vector<LRSplineVolume::Refinement3D>& refs_y,
				   vector<LRSplineVolume::Refinement3D>& refs_z,
				   vector<Element3D*>& affected)
//==============================================================================
{
  double tol = vol_->getKnotTol();
  double eps = 0.01;

  // Fetch B-splines
  const vector<LRBSpline3D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();

  int degree1 = vol_->degree(XDIR);
  int degree2 = vol_->degree(YDIR);
  int degree3 = vol_->degree(ZDIR);
  int xmult = (degree1 <= 3) ? 1 : 2;
  int ymult = (degree2 <= 3) ? 1 : 2;
  int zmult = (degree3 <= 3) ? 1 : 2;
  
 // Refine one B-spline with support in the parent element in one or two
  // parameter directions depending on how many elements with a large error
  // that lies in its support
  // First check refinement in the u-direction
  double u_par = 0.5*(elem->umin() + elem->umax());
  double udel = elem->umax() - elem->umin();
  int ixu = -1;
  double max_wgt = 0.0;
  double minsize_u = (usize_min_ > 0.0) ? 2.0*usize_min_ : 1.0e-8;
  size_t ki, kj;
  for (ki=0; ki<nmb; ++ki)
    {
      // Count the number of elements with large error affected
      double curr_wgt = 0.0;
      const vector<Element3D*>& curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->umax() < u_par || curr_el[kj]->umin() > u_par)
	    continue;  // Element not affected

	  if (curr_el[kj]->umax() - curr_el[kj]->umin() < minsize_u)
	    {
	      // Element too small to be divided
	      break;
	    }
	  
	  // Compute weight for importance of refinement
	  double max_err, av_err;
	  int nmb_outside, nmb_out_sign;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el[kj]->nmbDataPoints();
	  double wgt = av_err*(int)nmb_outside/(int)nmb_pts;
	  curr_wgt += wgt;
	}
      curr_wgt = (bsplines[ki]->umax() - bsplines[ki]->umin());
      if (kj < curr_el.size())
	continue;  // Current B-spline is not a candidate for a split
      
        if (curr_wgt > max_wgt)
	{
	  max_wgt = curr_wgt;
	  ixu = (int)ki;
	}
    }
	  
  // The v-direction
  double v_par = 0.5*(elem->vmin() + elem->vmax());
  double vdel = elem->vmax() - elem->vmin();
  int ixv = -1;
  max_wgt = 0.0;
  double minsize_v = (vsize_min_ > 0.0) ? 2.0*vsize_min_ : 1.0e-8;
  for (ki=0; ki<nmb; ++ki)
    {
      // Count the number of elements with large error affected
      double curr_wgt = 0.0;
      const vector<Element3D*>& curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->vmax() < v_par || curr_el[kj]->vmin() > v_par)
	    continue;  // Element not affected

	  if (curr_el[kj]->vmax() - curr_el[kj]->vmin() < minsize_v)
	    {
	      // Element too small to be divided
	      break;
	    }
	  
	  // Compute weight for importance of refinement
	  double max_err, av_err;
	  int nmb_outside, nmb_out_sign;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el[kj]->nmbDataPoints();
	  double wgt = av_err*(int)nmb_outside/(int)nmb_pts;
	  curr_wgt += wgt;
	}
      curr_wgt = (bsplines[ki]->vmax() - bsplines[ki]->vmin());
      if (kj < curr_el.size())
	continue;  // Current B-spline is not a candidate for a split
      
       if (curr_wgt > max_wgt)
	{
	  max_wgt = curr_wgt;
	  ixv = (int)ki;
	}
    }
  
   // The w-direction
  double w_par = 0.5*(elem->wmin() + elem->wmax());
  double wdel = elem->wmax() - elem->wmin();
  int ixw = -1;
  max_wgt = 0.0;
  double minsize_w = (wsize_min_ > 0.0) ? 2.0*wsize_min_ : 1.0e-8;
  for (ki=0; ki<nmb; ++ki)
    {
      // Count the number of elements with large error affected
      double curr_wgt = 0.0;
      const vector<Element3D*>& curr_el = bsplines[ki]->supportedElements();
      for (kj=0; kj<curr_el.size(); ++kj)
	{
	  if (curr_el[kj]->wmax() < w_par || curr_el[kj]->wmin() > w_par)
	    continue;  // Element not affected

	  if (curr_el[kj]->wmax() - curr_el[kj]->wmin() < minsize_w)
	    {
	      // Element too small to be divided
	      break;
	    }
	  
	  // Compute weight for importance of refinement
	  double max_err, av_err;
	  int nmb_outside, nmb_out_sign;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el[kj]->nmbDataPoints();
	  double wgt = av_err*(int)nmb_outside/(int)nmb_pts;
	  curr_wgt += wgt;
	}
      curr_wgt = (bsplines[ki]->wmax() - bsplines[ki]->wmin());
     if (kj < curr_el.size())
	continue;  // Current B-spline is not a candidate for a split
      
      if (curr_wgt > max_wgt)
	{
	  max_wgt = curr_wgt;
	  ixw = (int)ki;
	}
    }
  
 // Fetch the affected elements in both parameter directions
  double fac = 0.1;
  double fac3 = 0.75; //0.95;
  vector<Element3D*> aff_u;
  int nmb_u = 0;
  if (ixu >= 0)
    {
      const vector<Element3D*>& curr_el_u = bsplines[ixu]->supportedElements();
      for (kj=0; kj<curr_el_u.size(); ++kj)
	{
	  if (curr_el_u[kj]->umin() > u_par || curr_el_u[kj]->umax() < u_par)
	    continue;
	  double max_err, av_err;
	  int nmb_outside;
	  curr_el_u[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el_u[kj]->nmbDataPoints();
	  aff_u.push_back(curr_el_u[kj]);
	  if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	    nmb_u++;
	}
    }

  vector<Element3D*> aff_v;
  int nmb_v = 0;
  if (ixv >= 0)
    {
      const vector<Element3D*>& curr_el_v = bsplines[ixv]->supportedElements();
      for (kj=0; kj<curr_el_v.size(); ++kj)
	{
	  if (curr_el_v[kj]->vmin() > v_par || curr_el_v[kj]->vmax() < v_par)
	    continue;
	  double max_err, av_err;
	  int nmb_outside;
	  curr_el_v[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el_v[kj]->nmbDataPoints();
	  aff_v.push_back(curr_el_v[kj]);
	  if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	    nmb_v++;
	}
    }

  vector<Element3D*> aff_w;
  int nmb_w = 0;
  if (ixw >= 0)
    {
      const vector<Element3D*>& curr_el_w = bsplines[ixw]->supportedElements();
      for (kj=0; kj<curr_el_w.size(); ++kj)
	{
	  if (curr_el_w[kj]->wmin() > w_par || curr_el_w[kj]->wmax() < w_par)
	    continue;
	  double max_err, av_err;
	  int nmb_outside;
	  curr_el_w[kj]->getAccuracyInfo(av_err, max_err, nmb_outside);
	  int nmb_pts = curr_el_w[kj]->nmbDataPoints();
	  aff_w.push_back(curr_el_w[kj]);
	  if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	    nmb_w++;
	}
    }

  // if (elem->getNmbOutsideTol() == 1)
  //   {
  //     std::cout << "Element: " << elem->umin() << " " << elem->umax() << " ";
  //     std::cout << elem->vmin() << " " << elem->vmax() << " ";
  //     std::cout << elem->wmin() << " " << elem->wmax() << " " << elem->getMaxError() << std::endl;;
  //     std::cout << "ixu: " << ixu << ", ixv: " << ixv << ", ixw: " << ixw << std::endl;
  //     std::cout << aff_u.size() << " " << aff_v.size() << " " << aff_w.size() << std::endl;
  //   }
  
  // Assemble information
  int nmb_div = 0;
  double fac2 = 0.5;
  double sizefac = 3.0;
  std::set<Element3D*> affected_combined;
  if (ixu >= 0 && ((nmb_u >= nmb_v && nmb_u >= nmb_w) ||
		   udel > sizefac*std::min(vdel, wdel)))
    {
      affected_combined.insert(aff_u.begin(), aff_u.end());
      LRSplineVolume::Refinement3D curr_ref;
      curr_ref.setVal(u_par, bsplines[ixu]->vmin(), bsplines[ixu]->vmax(),
		      bsplines[ixu]->wmin(), bsplines[ixu]->wmax(), XDIR, xmult);
      //refs.push_back(curr_ref);
      appendRef(refs_x, curr_ref, tol);
      ref_x_++;
      nmb_div++;
    }
  // else if (ixu >= 0)
  //   {
  //     std::cout << "u, nmbu = " << nmb_u << ", nmbv = " << nmb_v << ", nmbw = ";
  //     std::cout << nmb_w << ", udel = " << udel << ", vdel = " << vdel;
  //     std::cout << ", wdel = " << wdel << std::endl;
  //   }
			       
    if (ixv >= 0 && ((nmb_v >= nmb_u && nmb_v >= nmb_w) ||
		     vdel > sizefac*std::min(udel, wdel)))
    {
      affected_combined.insert(aff_v.begin(), aff_v.end());
      LRSplineVolume::Refinement3D curr_ref;
      // curr_ref.setVal(v_par, bsplines[ixv]->umin(), bsplines[ixv]->umax(),
      // 		      bsplines[ixv]->wmin(), bsplines[ixv]->wmax(), YDIR, ymult);
      curr_ref.setVal(v_par, bsplines[ixv]->wmin(), bsplines[ixv]->wmax(),
		      bsplines[ixv]->umin(), bsplines[ixv]->umax(), YDIR, ymult);
      //refs.push_back(curr_ref);
      appendRef(refs_y, curr_ref, tol);
      ref_y_++;
      nmb_div++;
    }
  // else if (ixv >= 0)
  //   {
  //     std::cout << "v, nmbu = " << nmb_u << ", nmbv = " << nmb_v << ", nmbw = ";
  //     std::cout << nmb_w << ", udel = " << udel << ", vdel = " << vdel;
  //     std::cout << ", wdel = " << wdel << std::endl;
  //   }

    if (ixw >= 0 && ((nmb_w >= nmb_u && nmb_w >= nmb_v) ||
		     wdel > sizefac*std::min(udel, vdel)))
    {
      affected_combined.insert(aff_w.begin(), aff_w.end());
      LRSplineVolume::Refinement3D curr_ref;
      curr_ref.setVal(w_par, bsplines[ixw]->umin(), bsplines[ixw]->umax(),
		      bsplines[ixw]->vmin(), bsplines[ixw]->vmax(), ZDIR, zmult);
      //refs.push_back(curr_ref);
      appendRef(refs_z, curr_ref, tol);
      ref_z_++;
      nmb_div++;
    }
  // else if (ixw >= 0)
  //   {
  //     std::cout << "w, nmbu = " << nmb_u << ", nmbv = " << nmb_v << ", nmbw = ";
  //     std::cout << nmb_w << ", udel = " << udel << ", vdel = " << vdel;
  //     std::cout << ", wdel = " << wdel << std::endl;
  //   }

  //affected.insert(affected.end(), affected_combined.begin(), affected_combined.end());
  if (nmb_div == 1)
    nmb1_++;
  else if (nmb_div == 2)
    nmb2_++;
  else if (nmb_div == 3)
    nmb3_++;
}
#endif
/*
//==============================================================================
void LRVolApprox::refineSurf2()
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
      srf_->refine(refs[kr], true );
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
*/
 
//==============================================================================
void LRVolApprox::makeInitVol(int dim)
//==============================================================================
{
  // Create a cubic Bezier surface represented as an LR B-spline surface
  // Compute domain
  double umin, umax, vmin, vmax, wmin, wmax;
  computeParDomain(dim, umin, umax, vmin, vmax, wmin, wmax);


  int order = 4;
  vector<double> knots_u(2*order);
  vector<double> knots_v(2*order);
  vector<double> knots_w(2*order);

  int kj;
  for (kj=0; kj<order; ++kj)
    {
      knots_u[kj] = umin;
      knots_v[kj] = vmin;
      knots_w[kj] = wmin;
      knots_u[kj+order] = umax;
      knots_v[kj+order] = vmax;
      knots_w[kj+order] = wmax;
    }

  makeInitVol(dim, order, order, order, order, order, order, &knots_u[0], &knots_v[0], &knots_w[0]);
}
 
 /*
//==============================================================================
void LRVolApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
				int order_v)
//==============================================================================
{
  // Compute domain
  double domain[4]; //umin, umax, vmin, vmax;
  computeParDomain(dim, domain[0], domain[1], domain[2], domain[3]);

  makeInitSurf(dim, ncoef_u, order_u, ncoef_v, order_v, domain);
}
*/
 /*
//==============================================================================
void LRVolApprox::makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
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
 */
 
//==============================================================================
void LRVolApprox::makeInitVol(int dim, 
                              int ncoef_u, int order_u, 
                              int ncoef_v, int order_v, 
                              int ncoef_w, int order_w, 
                              double *knots_u, double *knots_v, double *knots_w)
//==============================================================================
{
  int nmb = (int)points_.size()/(dim+3);
  shared_ptr<SplineVolume> result_vol;

  // std::cout << initMBA_ << std::endl;
  // if (!initMBA_)
  //   {
  //     try {
  // 	result_vol = createVol(&points_[0], nmb, 
  // 			       dim, ncoef_u, order_u,
  // 			       ncoef_v, order_v, ncoef_w, order_w,
  // 			       knots_u, knots_v, knots_w); 
  //     }
  //     catch (...)
  // 	{
  // 	  initMBA_ = true;
  // 	  std::cout << initMBA_ << std::endl;
  // 	}
  //   }

  // if (initMBA_)
  //   {
      vector<double> coefs(ncoef_u*ncoef_v*ncoef_w*dim, initMBA_coef_);
      result_vol = shared_ptr<SplineVolume>(new SplineVolume(ncoef_u, ncoef_v, ncoef_w,
							     order_u, order_v, order_w,
							     knots_u,     // ptr to knots in u
							     knots_v,     // ptr to knots in v
							     knots_w,     // ptr to knots in w
							     &coefs[0], // ptr to coefs
							     dim));
      //    }
      
  // Make LR spline volume
  double knot_tol = 1.0e-6;

  vol_ = shared_ptr<LRSplineVolume>(new LRSplineVolume(result_vol.get(), knot_tol));
  coef_known_.assign(vol_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}

 
//==============================================================================
shared_ptr<SplineVolume> LRVolApprox::createVol(double* points, int nmb_pts,
						int dim, int ncoef_u, int order_u, 
						int ncoef_v, int order_v, 
						int ncoef_w, int order_w, 
						double *knots_u, double *knots_v,
						double *knots_w)
//==============================================================================
{
  shared_ptr<SplineVolume> result_vol;
  vector<double> coefs(ncoef_u*ncoef_v*ncoef_w*dim, 0.0);
  shared_ptr<SplineVolume> vol1(new SplineVolume(ncoef_u, ncoef_v, ncoef_w,
						 order_u, order_v, order_w,
						 knots_u,     // ptr to knots in u
						 knots_v,     // ptr to knots in v
						 knots_w,     // ptr to knots in w
						 &coefs[0], // ptr to coefs
						 dim));

  // Reformatting data to the format that ApproxSurf wants (separating parameter values and 
  // data points in two distinct vectors).
  vector<double> pts, param;
  pts.reserve(dim*nmb_pts);
  param.reserve(3* nmb_pts);
  int ki, kj;
  double* it;
  for (it=points, kj=0; kj<nmb_pts; ++kj) 
    {
      for (ki=0; ki<3; ++ki)
	param.push_back(*it++);
      for (ki=0; ki<dim; ++ki)
	pts.push_back(*it++);
    }
    
  // Approximate data points to create initial surface
  SmoothVolume avol;  // Engine for least squares approximation with smoothing
  int stat = 0;

  // Define weights
  int min_der = std::max(3, std::min(order_u, std::min(order_v, order_w))-1);
  double wgt1 = 0.0;
  double wgt3 = (min_der >= 3) ? 0.5*smoothweight_ : 0.0;
  double wgt2 = (1.0 - wgt3)*smoothweight_;
  wgt3 *= smoothweight_;
  double approxweight = 1.0 - wgt1 - wgt2 - wgt3;
  std::vector<double> pt_weight(nmb_pts, 1.0);
 
  // Prepare for approximation
  vector<CoefStatus> coef_stat(ncoef_u*ncoef_v*ncoef_w, CoefFree);
  avol.attach(vol1, coef_stat);

  if (smoothweight_ > 0.0)
    avol.setOptimize(wgt1, wgt2, wgt3);

  avol.setLeastSquares(pts, param, pt_weight, approxweight);
  
  // Approximate
  stat = avol.equationSolve(result_vol);
  return result_vol;
  std::cout << "Initial volume: " << stat << std::endl;
}
 
 /*
//==============================================================================
void LRVolApprox::makeInitSurf(shared_ptr<SplineSurface> surf)
//==============================================================================
{
  // Make LR spline surface
  double knot_tol = 1.0e-6;

  srf_ = shared_ptr<LRSplineSurface>(new LRSplineSurface(surf.get(), knot_tol));
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}
 */
 
//==============================================================================
void LRVolApprox::computeParDomain(int dim, 
                                   double& umin, double& umax, 
                                   double& vmin, double& vmax,
                                   double& wmin, double& wmax )
//==============================================================================
{
  // Compute domain
  umin = umax = points_[0];
  vmin = vmax = points_[1];
  wmin = wmax = points_[2];
  int del = 3 + dim;
 for (size_t ki=del; ki<points_.size(); ki+=del)
    {
      umin = std::min(umin, points_[ki]);
      umax = std::max(umax, points_[ki]);
      vmin = std::min(vmin, points_[ki+1]);
      vmax = std::max(vmax, points_[ki+1]);
      wmin = std::min(wmin, points_[ki+2]);
      wmax = std::max(wmax, points_[ki+2]);
    }

}
 
/*
//==============================================================================
void LRVolApprox::setCoefKnown()
//==============================================================================
{
  MESSAGE("CHECK THIS FUNCTION: setCoefKnown()");
  if (fix_corner_)
    {
      // Fix coefficients in corners
      for (LRSplineVolume::BSplineMap::iterator it=srf_->basisFunctionsBeginNonconst();
           it != vol_->basisFunctionsEndNonconst(); ++it)
	{
	  int deg1 = it->second->degree(XDIR);
	  int deg2 = it->second->degree(YDIR);
	  int deg3 = it->second->degree(ZDIR);
	  int mult1 = it->second->endmult_u(true);
	  int mult2 = it->second->endmult_u(false);
	  int mult3 = it->second->endmult_v(true);
	  int mult4 = it->second->endmult_v(false);
	  int mult5 = it->second->endmult_w(true);
	  int mult6 = it->second->endmult_w(false);
	  if ((mult1 == deg1+1 || mult2 == deg1+1) &&
	      (mult3 == deg2+1 || mult4 == deg2+1) &&
	      (mult5 == deg3+1 || mult6 == deg3+1))
	    it->second->setFixCoef(1);
	}
    }

  // Set coef_known from boundary information
  // Note that k-multiple knots at boundaries are expected, otherwise no
  // coefficients will be fixed
  for (size_t ki=0; ki<6; ++ki)
    {
      if (face_derivs_[ki] == 0)
	continue;  // No coefficients to fix
      
      // Set boundary characteristics
      Direction3D d = (ki == 1 || ki == 3 || ki == 5) ? XFIXED : YFIXED;  // Orthogonal to the curve
      bool atstart = (ki == 0 || ki == 3);  // Whether the constant
      // parameter curve is in the start related to the opposite parameter direction
      // of the surface

      // Define coefficients as fixed in the associated b-splines
      int fixcoef = 1;
      setCoefKnown(d, atstart, fixcoef);
    }

  // Transfer information to the global vector
  updateCoefKnown();
}
*/
/*
//==============================================================================
void LRVolApprox::setCoefKnown(Direction2D d,
				bool atstart, int fixcoef)
//==============================================================================
{
//#ifdef DEBUG
//  std::ofstream of("mesh.eps");
//  writePostscriptMesh(*vol_, of);
//#endif

  for (LRSplineVolume::BSplineMap::iterator it=vol_->basisFunctionsBeginNonconst();
       it != vol_->basisFunctionsEndNonconst(); ++it)
    {
      int deg = it->second->degree(d);
      int mult = (d == XDIR) ? it->second->endmult_u(atstart) :
                               ((d == YDIR) ? it->second->endmult_v(atstart) : it->second->endmult_v(atstart));
      if (mult == deg+1)
	it->second->setFixCoef(fixcoef);
    }
  
}
*/
 
//==============================================================================
void LRVolApprox::updateCoefKnown()
//==============================================================================
{
  // Fetch all boundary b-splines and check if the coefficient should be fixed
  int fixcoef = 1;
  for (LRSplineVolume::BSplineMap::iterator it=vol_->basisFunctionsBeginNonconst(); 
       it != vol_->basisFunctionsEndNonconst(); ++it)
    {
      int deg1 = it->second->degree(XDIR);
      int deg2 = it->second->degree(YDIR);
      int deg3 = it->second->degree(ZDIR);
      int mult1_1 = it->second->endmult_u(true);
      int mult1_2 = it->second->endmult_u(false);
      int mult2_1 = it->second->endmult_v(true);
      int mult2_2 = it->second->endmult_v(false);
      int mult3_1 = it->second->endmult_w(true);
      int mult3_2 = it->second->endmult_w(false);
      if (mult1_1 == deg1+1 || mult1_2 == deg1+1 || 
	  mult2_1 == deg2+1 || mult2_2 == deg2+1 || 
          mult3_1 == deg3+1 || mult3_2 == deg3+1 )
	{
	  const vector<Element3D*>& curr_el = it->second->supportedElements();
	  int nmb_pts = 0;
	  for (size_t ki=0; ki<curr_el.size(); ++ki) nmb_pts += curr_el[ki]->nmbDataPoints();
	  if (nmb_pts == 0)
	    it->second->setFixCoef(fixcoef);
	}
    }
  
  coef_known_.resize(vol_->numBasisFunctions());

  LRSplineVolume::BSplineMap::const_iterator it_bs;			
  size_t ki;
  for (it_bs=vol_->basisFunctionsBegin(), ki=0; it_bs!=vol_->basisFunctionsEnd(); 
       ++it_bs, ++ki)
    coef_known_[ki] = it_bs->second->coefFixed();
}
 
 /*
//==============================================================================
void LRVolApprox::unsetCoefKnown()
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
 */
 
//==============================================================================
void LRVolApprox::defineRefs(LRBSpline3D* bspline, double average_out,
                             vector<LRSplineVolume::Refinement3D>& refs_x,
                             vector<LRSplineVolume::Refinement3D>& refs_y,
                             vector<LRSplineVolume::Refinement3D>& refs_z,
			     vector<Element3D*>& elem_div)
//==============================================================================
{
  // For each alternative (knot span) in each parameter direction, collect
  // accuracy statistic
  // Compute also average element size
  double tol = vol_->getKnotTol();
  int size1 = bspline->degree(XDIR)+1;
  int size2 = bspline->degree(YDIR)+1;
  int size3 = bspline->degree(ZDIR)+1;
  int sizeall = size1 + size2 + size3;
  double scratch1[24];
  int scratch2[24];
  vector<double> scratch1_2;
  vector<int> scratch2_2;
  if (sizeall > 12)
    {
      scratch1_2.resize(2*sizeall, 0.0);
      scratch2_2.resize(2*sizeall, 0);
    }
  else
    {
      memset(scratch1, 0.0, 24*sizeof(double));
      memset(scratch2, 0, 24*sizeof(int));
    }
  double *u_info = (sizeall <= 12) ? scratch1 : &scratch1_2[0];
  double *v_info = u_info + size1;
  double *w_info = v_info + size2;
  double *vw_elsize = w_info + size3;
  double *wu_elsize = vw_elsize + size1;
  double *uv_elsize = wu_elsize + size2;
  int *u_inside = (sizeall <= 12) ? scratch2 : &scratch2_2[0];
  int *v_inside = u_inside + size1;
  int *w_inside = v_inside + size2;
  int *u_outside = w_inside + size3;
  int *v_outside = u_outside + size1;
  int *w_outside = v_outside + size2;
  
  // vector<double> u_info(size1, 0.0);
  // vector<double> v_info(size2, 0.0);
  // vector<double> w_info(size3, 0.0);
  // vector<int> u_inside(size1, 0);
  // vector<int> v_inside(size2, 0);
  // vector<int> w_inside(size3, 0);
  // vector<int> u_outside(size1, 0);
  // vector<int> v_outside(size2, 0);
  // vector<int> w_outside(size3, 0);

  // vector<double> vw_elsize(size1, 0.0); //
  // vector<double> wu_elsize(size2, 0.0); //
  // vector<double> uv_elsize(size3, 0.0); //
  
  const vector<int>& kvec_u = bspline->kvec(XDIR);
  const vector<int>& kvec_v = bspline->kvec(YDIR);
  const vector<int>& kvec_w = bspline->kvec(ZDIR);
  const Mesh3D* mesh = bspline->getMesh();

  double av_kdiff_u = (mesh->kval(XDIR, kvec_u[kvec_u.size()-1]) -
		       mesh->kval(XDIR, kvec_u[0]))/(double)(kvec_u.size()-1);
  double av_kdiff_v = (mesh->kval(YDIR, kvec_v[kvec_v.size()-1]) -
		       mesh->kval(YDIR, kvec_v[0]))/(double)(kvec_v.size()-1);
  double av_kdiff_w = (mesh->kval(ZDIR, kvec_w[kvec_w.size()-1]) -
		       mesh->kval(ZDIR, kvec_w[0]))/(double)(kvec_w.size()-1);

  const vector<Element3D*>& elem = bspline->supportedElements();
  int nmb_outside_pts = 0;
  int curr_nmb_out, curr_nmb;
  double dom;
  bool refined = false;
  size_t ki;
  for (ki=0; ki<elem.size(); ++ki)
    {
      // Localize element with regard to the information containers
      double umin = elem[ki]->umin();
      double umax = elem[ki]->umax();
      double vmin = elem[ki]->vmin();
      double vmax = elem[ki]->vmax();
      double wmin = elem[ki]->wmin();
      double wmax = elem[ki]->wmax();

      size_t kj1, kj2, kj3;
      for (kj1=1; kj1<kvec_u.size(); ++kj1)
	if (mesh->kval(XDIR, kvec_u[kj1-1]) <= umin && mesh->kval(XDIR, kvec_u[kj1]) >= umax)
	  break;
      for (kj2=1; kj2<kvec_v.size(); ++kj2)
	if (mesh->kval(YDIR, kvec_v[kj2-1]) <= vmin && mesh->kval(YDIR, kvec_v[kj2]) >= vmax)
	  break;
      for (kj3=1; kj3<kvec_w.size(); ++kj3)
	if (mesh->kval(ZDIR, kvec_w[kj3-1]) <= wmin && mesh->kval(ZDIR, kvec_w[kj3]) >= wmax)
	  break;

      curr_nmb = elem[ki]->nmbDataPoints();
      curr_nmb_out = elem[ki]->getNmbOutsideTol();
      if (curr_nmb_out > 0 && (double)curr_nmb_out/(double)curr_nmb >= outfrac_)
	{
	  nmb_outside_pts += curr_nmb_out;
	  dom = (umax-umin)*(vmax-vmin)*(wmax-wmin);
	  double accout = elem[ki]->getAccumulatedOutside();
	  double maxout = elem[ki]->getMaxError();
	  double avout = accout/(double)curr_nmb_out;
	  if (maxout > avout + 0.9*(maxout-avout))
	    accout *= 2.0;
	  u_info[kj1-1] += dom*accout;
	  v_info[kj2-1] += dom*accout;
	  w_info[kj3-1] += dom*accout;
	  if (umax-umin > 0.9*(kvec_u[kj1]-kvec_u[kj1-1]))
	    u_outside[kj1-1] += curr_nmb_out;
	  if (vmax-vmin > 0.9*(kvec_v[kj2]-kvec_v[kj2-1]))
	    v_outside[kj2-1] += curr_nmb_out;
	  if (wmax-wmin > 0.9*(kvec_w[kj3]-kvec_w[kj3-1]))
	    w_outside[kj3-1] += curr_nmb_out;
	}
      else if (curr_nmb_out == 0)
	{
	  if (umax-umin > 0.9*(kvec_u[kj1]-kvec_u[kj1-1]))
	    u_inside[kj1-1]++;
	  if (vmax-vmin > 0.9*(kvec_v[kj2]-kvec_v[kj2-1]))
	    v_inside[kj2-1]++;
	  if (wmax-wmin > 0.9*(kvec_w[kj3]-kvec_w[kj3-1]))
	    w_inside[kj3-1]++;
	}
      // if (elem[ki]->getNmbOutsideTol() > 0)
      // 	{
      // 	  u_info[kj1-1] += (umax-umin)*(vmax-vmin)*(wmax-wmin)*elem[ki]->getAccumulatedError();
      // 	  v_info[kj2-1] += (umax-umin)*(vmax-vmin)*(wmax-wmin)*elem[ki]->getAccumulatedError();
      //     w_info[kj3-1] += (umax-umin)*(vmax-vmin)*(wmax-wmin)*elem[ki]->getAccumulatedError();
      // 	}

      // Element size
      vw_elsize[kj1-1] += (vmax-vmin)*(wmax-wmin);
      wu_elsize[kj2-1] += (wmax-wmin)*(umax-umin);
      uv_elsize[kj3-1] += (umax-umin)*(vmax-vmin);      
    } 

  // Set threshold for which strips to split
  double max_info = 0.0;
  double av_info = 0.0;
  int kj;
  for (kj=0; kj<size1; ++kj)
    {
      max_info = std::max(max_info, u_info[kj]);
      av_info += u_info[kj];
      vw_elsize[kj] /= (double)(size2*size3);
      vw_elsize[kj] = sqrt(vw_elsize[kj]);
    }
  for (kj=0; kj<size2; ++kj)
    {
      max_info = std::max(max_info, v_info[kj]);
      av_info += v_info[kj];
      wu_elsize[kj] /= (double)(size3*size1);
      wu_elsize[kj] = sqrt(wu_elsize[kj]);
    }
  for (kj=0; kj<size3; ++kj)
    {
      max_info = std::max(max_info, w_info[kj]);
      av_info += w_info[kj];
      uv_elsize[kj] /= (double)(size1*size2);
      uv_elsize[kj] = sqrt(uv_elsize[kj]);
    }

  // Modify priority information of strips to reduce the weight towards
  // the ends of the b-spline
  // double fac1 = 0.25; //1.0; //0.25;
  // double fac2 = 0.5;  //1.0; //0.5;
  // if (size1 >= 3)
  //   {
  //     u_info[0] *= fac1;
  //     u_info[size1-1] *= fac1;
  //   }
  // if (size1 >= 4)
  //   {
  //     u_info[1] *= fac2;
  //     u_info[size1-2] *= fac2;
  //   }
  // if (size2 >= 3)
  //   {
  //     v_info[0] *= fac1;
  //     v_info[size2-1] *= fac1;
  //   }
  // if (size2 >= 4)
  //   {
  //     v_info[1] *= fac2;
  //     v_info[size2-2] *= fac2;
  //   }
  // if (size3 >= 3)
  //   {
  //     w_info[0] *= fac1;
  //     w_info[size3-1] *= fac1;
  //   }
  // if (size3 >= 4)
  //   {
  //     w_info[1] *= fac2;
  //     w_info[size3-2] *= fac2;
  //   }

  av_info /= (double)(size1+size2+size3);
  int div_x = 0, div_y = 0, div_z = 0;
  std::set<Element3D*> curr_el;
  //double threshold = std::min(av_info, 0.5*max_info);
  double threshold = std::min(1.1*av_info, 0.8*max_info);
  double sizefac = 1.5; //3.0;
  double sizefac2 = 2.0;
  double minsize_u = std::max(2.0*usize_min_, 1.0e-8);
  for (kj=0; kj<size1; ++kj)
    {
      if (kvec_u[kj] == kvec_u[kj+1])
	continue;
      double u1 = mesh->kval(XDIR, kvec_u[kj]);
      double u2 = mesh->kval(XDIR, kvec_u[kj+1]);
      if (((u_info[kj] >= threshold || u2-u1 > sizefac*vw_elsize[kj]) &&
	   (u2 - u1) >= minsize_u && 
	   (u_inside[kj] == 0 || (double)u_outside[kj] > average_out)) ||
	  u2 - u1 > sizefac2*(av_kdiff_u + av_kdiff_v + av_kdiff_w))
	{
	  LRSplineVolume::Refinement3D curr_ref;
	  double par = 0.5*(u1+u2);
	  curr_ref.setVal(par, 
                          bspline->vmin(), bspline->vmax(),
                          bspline->wmin(), bspline->wmax(), 
                          XDIR, 1);

	  appendRef(refs_x, curr_ref, tol);
	  for (ki=0; ki<elem.size(); ++ki)
	    {
	      if (elem[ki]->umin() < par && elem[ki]->umax() > par)
		curr_el.insert(elem[ki]);
	    }
	  refined = true;
	  ref_x_++;
	  div_x = 1;
	}
    }

  double minsize_v = std::max(2.0*vsize_min_, 1.0e-8);
  for (kj=0; kj<size2; ++kj)
    {
      if (kvec_v[kj] == kvec_v[kj+1])
	continue;
      double v1 = mesh->kval(YDIR, kvec_v[kj]);
      double v2 = mesh->kval(YDIR, kvec_v[kj+1]);
      if (((v_info[kj] >= threshold  || v2-v1 > sizefac*wu_elsize[kj]) &&
	   (v2 - v1) >= minsize_v && 
	   (v_inside[kj] == 0 || (double)v_outside[kj] > average_out)) ||
	  v2 - v1 > sizefac2*(av_kdiff_u + av_kdiff_v + av_kdiff_w))
	{
	  LRSplineVolume::Refinement3D curr_ref;
	  double par = 0.5*(v1+v2);
	  curr_ref.setVal(par,
			  bspline->wmin(), bspline->wmax(),
			  bspline->umin(), bspline->umax(),
			  YDIR, 1);

	  appendRef(refs_y, curr_ref, tol);
	  for (ki=0; ki<elem.size(); ++ki)
	    {
	      if (elem[ki]->vmin() < par && elem[ki]->vmax() > par)
		curr_el.insert(elem[ki]);
	    }
	  refined = true;
	  ref_y_++;
	  div_y = 1;
        }
    }
  
  double minsize_w = std::max(2.0*wsize_min_, 1.0e-8);
  for (kj=0; kj<size3; ++kj)
    {
      if (kvec_w[kj] == kvec_w[kj+1])
	continue;
      double w1 = mesh->kval(ZDIR, kvec_w[kj]);
      double w2 = mesh->kval(ZDIR, kvec_w[kj+1]);
      if (((w_info[kj] >= threshold  || w2-w1 > sizefac*uv_elsize[kj]) &&
	   (w2 - w1) >= minsize_w && 
	   (w_inside[kj] == 0 || (double)w_outside[kj] > average_out)) ||
	  w2 - w1 > sizefac2*(av_kdiff_u + av_kdiff_v + av_kdiff_w))
	{
	  LRSplineVolume::Refinement3D curr_ref;
	  double par = 0.5*(w1+w2);
	  curr_ref.setVal(par,
                          bspline->umin(), bspline->umax(), 
                          bspline->vmin(), bspline->vmax(), 
                          ZDIR, 1);

	  appendRef(refs_z, curr_ref, tol);
	  for (ki=0; ki<elem.size(); ++ki)
	    {
	      if (elem[ki]->wmin() < par && elem[ki]->wmax() > par)
		curr_el.insert(elem[ki]);
	    }
	  refined = true;
	  ref_z_++;
 	  div_z = 1;
       }
    }

  // std::cout << "Start remove elements, refine = " << refined << ", elem size: ";
  // std::cout << elem.size() << ", elem out size: " << elem_out.size() << std::endl;
  if (refined)
    {
      elem_div.insert(elem_div.end(), curr_el.begin(), curr_el.end());
      // for (size_t ki=0; ki<elem.size(); ++ki)
      // 	{
      // 	  size_t kj;
      // 	  for (kj=0; kj<elem_out.size(); ++kj)
      // 	    if (elem[ki] == elem_out[kj].first)
      // 	      break;
      // 	  if (kj < elem_out.size())
      // 	    elem_out.erase(elem_out.begin() + kj);
      // 	}
    }
  // std::cout << "End remove elements" << std::endl;
  int nmb_div = div_x + div_y + div_z;
  if (nmb_div == 1)
    nmb1_++;
  else if (nmb_div == 2)
    nmb2_++;
  else if (nmb_div == 3)
    nmb3_++;
}
 
//==============================================================================
void LRVolApprox::appendRef(vector<LRSplineVolume::Refinement3D>& refs,
			     LRSplineVolume::Refinement3D& curr_ref,
			     double tol)
//==============================================================================
{
  // Check if the current refinement can be combined with an existing one
  size_t ki;
  for (ki=0; ki<refs.size(); ++ki)
    {
      // Check direction and knot value
      if (/*refs[ki].d == curr_ref.d &&*/ 
	  fabs(refs[ki].kval-curr_ref.kval) < tol)
	{
	  // Ensure equality
	  curr_ref.kval = refs[ki].kval;
	  
	  // Check extent std::cout refinement
	  if (!(refs[ki].start1 > curr_ref.end1+tol ||
		curr_ref.start1 > refs[ki].end1+tol || 
		refs[ki].start2 > curr_ref.end2+tol ||
		curr_ref.start2 > refs[ki].end2+tol ))
	    {
	      // Merge new knots
	      refs[ki].start1 = std::min(refs[ki].start1, curr_ref.start1);
	      refs[ki].start2 = std::min(refs[ki].start2, curr_ref.start2);
	      refs[ki].end1 = std::max(refs[ki].end1, curr_ref.end1);
	      refs[ki].end2 = std::max(refs[ki].end2, curr_ref.end2);
	      break;
	    }
	}
    }

  if (ki == refs.size())
    refs.push_back(curr_ref);
}

 /*
//==============================================================================
void LRVolApprox::checkFeasibleRef(Element2D* elem, 
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
	  if (true ) // curr_el_u[kj]->umax() > u_par && curr_el_u[kj]->umin() < u_par
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
	  if (true ) // curr_el_v[kj]->vmax() > v_par && curr_el_v[kj]->vmin() < v_par
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
*/

 /*
//==============================================================================
void LRVolApprox::adaptSurfaceToConstraints()
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
 */
/*
//==============================================================================
void LRVolApprox::turnTo3D()
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

 */
