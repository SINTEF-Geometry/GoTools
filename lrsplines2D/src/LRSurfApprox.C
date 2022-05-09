/* (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
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
#include "GoTools/lrsplines2D/LRFeatureUtils.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/lrsplines2D/LRFeatureUtils.h"
#include "GoTools/lrsplines2D/LogLikelyhood.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "stdio.h"
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

//#define DEBUG
//#define DEBUG1
//#define DEBUG2
//#define DEBUG_SURF
//#define DEBUG_DIST
//#define DEBUG_REFINE

using std::vector;
using std::set;
using std::cout;
using std::endl;
using std::pair;
using std::make_pair;
using namespace Go;

//==============================================================================
LRSurfApprox::LRSurfApprox(vector<double>& points, 
			   int dim, double epsge,  bool init_mba, 
			   double mba_level,
			   bool closest_dist, bool repar)
  : nmb_pts_((int)points.size()/(2+dim)), points_(points), useMBA_(false), 
    toMBA_(4), initMBA_(init_mba), initMBA_coef_(mba_level), 
    aepsge_(epsge), repar_(repar), check_close_(closest_dist), 
    check_init_accuracy_(false), initial_surface_(false)
//==============================================================================
{
  initDefaultParams();

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
    initMBA_coef_(0.0), aepsge_(epsge), 
    repar_(repar), check_close_(closest_dist), check_init_accuracy_(false),
    initial_surface_(true)
//==============================================================================
{
  initDefaultParams();

  nmb_pts_ = (int)points.size()/(2+srf->dimension());

  // Create an LR B-spline surface based on the given spline surface
  makeInitSurf(srf);
}

//==============================================================================
LRSurfApprox::LRSurfApprox(shared_ptr<LRSplineSurface>& srf,
			   vector<double>& points, 
			   double epsge, bool init_mba, double mba_level,
			   bool closest_dist, bool repar)
//==============================================================================
  : srf_(srf), points_(points), useMBA_(false), 
    toMBA_(4), initMBA_(init_mba), 
    initMBA_coef_(mba_level), aepsge_(epsge),  repar_(repar), 
    check_close_(closest_dist), check_init_accuracy_(false), 
    initial_surface_(true)
{
  initDefaultParams();

  nmb_pts_ = 0;  // Points already distributed

  srf_ = srf;
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed
}

//==============================================================================
LRSurfApprox::LRSurfApprox(shared_ptr<LRSplineSurface>& srf,
			   vector<double>& points, 
			   double epsge, bool closest_dist,
			   bool repar, bool check_init_accuracy)
  : points_(points), useMBA_(false), toMBA_(4), 
    initMBA_(false), initMBA_coef_(0.0), aepsge_(epsge), 
    repar_(repar), check_close_(closest_dist), 
    check_init_accuracy_(check_init_accuracy), 
    initial_surface_(true)
//==============================================================================
{
  initDefaultParams();

  nmb_pts_ = (int)points.size()/(2+srf->dimension());

  srf_ = srf;
  coef_known_.assign(srf_->numBasisFunctions(), 0.0);  // Initially nothing is fixed

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
    aepsge_(epsge), repar_(repar), check_close_(closest_dist), 
    check_init_accuracy_(false), initial_surface_(false)
{
  initDefaultParams();

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
    aepsge_(epsge), repar_(repar), check_close_(closest_dist), 
    check_init_accuracy_(false), initial_surface_(false)
{
  initDefaultParams();

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
    aepsge_(epsge), repar_(repar), check_close_(closest_dist), 
    check_init_accuracy_(false), initial_surface_(false)
{
  initDefaultParams();

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
void LRSurfApprox::getOutlierPts(vector<double>& outliers, int& nmb_outliers)
//==============================================================================
{
  LRSplineSurface::ElementMap::const_iterator it;
  int num = srf_->numElements();
  int kj;

  for (it=srf_->elementsBegin(), kj=0; kj<num; ++it, ++kj)
    {
      if (!it->second->hasDataPoints())
	continue;

      vector<double> outliers_local;
      it->second->getOutlierPts(outliers_local);
      outliers.insert(outliers.end(), outliers_local.begin(), outliers_local.end());
    }
  nmb_outliers = (int)outliers.size()/3;
}


//==============================================================================
void LRSurfApprox::getRegularPts(vector<double>& regular, int& nmb_regular)
//==============================================================================
{
  LRSplineSurface::ElementMap::const_iterator it;
  int num = srf_->numElements();
  int kj;

  for (it=srf_->elementsBegin(), kj=0; kj<num; ++it, ++kj)
    {
      if (!it->second->hasDataPoints())
	continue;

      vector<double> regular_local;
      it->second->getRegularPts(regular_local);
      regular.insert(regular.end(), regular_local.begin(), regular_local.end());
    }
  nmb_regular = (int)regular.size()/3;
}


//==============================================================================
void LRSurfApprox::getClassifiedPts(vector<double>& outliers, int& nmb_outliers,
				    vector<double>& regular, int& nmb_regular)
//==============================================================================
{
  LRSplineSurface::ElementMap::const_iterator it;
  int num = srf_->numElements();
  int kj;

  for (it=srf_->elementsBegin(), kj=0; kj<num; ++it, ++kj)
    {
      if (!it->second->hasDataPoints())
	continue;

      vector<double> outliers_local, regular_local;
      it->second->getClassifiedPts(outliers_local, regular_local);
      outliers.insert(outliers.end(), outliers_local.begin(), outliers_local.end());
      regular.insert(regular.end(), regular_local.begin(), regular_local.end());
    }
  nmb_outliers = (int)outliers.size()/3;
  nmb_regular = (int)regular.size()/3;
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
    const bool omp_for_elements = true; // 201503 The omp version seems to be faster even when run sequentially.
    const bool omp_for_mba_update = true;//false; // 201503 The omp version seems to be faster even when run sequentially.
#endif

#ifdef DEBUG_SURF
  std::ofstream of0("init0_sf.g2");
  std::ofstream of0el("init0_el.g2");
  srf_->writeStandardHeader(of0);
  srf_->write(of0);
  of0 << std::endl;
  LineCloud lines0 = srf_->getElementBds();
  lines0.writeStandardHeader(of0el);
  lines0.write(of0el);
  // std::ofstream of02("init0_tpsf.g2");
  // shared_ptr<SplineSurface> ssf0(tmp0->asSplineSurface());
  // ssf0->writeStandardHeader(of02);
  // ssf0->write(of02);
#endif

#ifdef DEBUG_REFINE
  FILE *fp = fopen("acc_stat.txt","w");
  fprintf(fp, "Max iterations = %d, tolerance = %4.2f, no pts: %d \n",max_iter, aepsge_,nmb_pts_);
  fprintf(fp,"iter, maxdist, average dist, no. pts. out, no. coefs, rel. improvement, no. pts.in, approx efficiency, rel element without-element div, rel element under-element div, max inner knots, average inner knots, average out, no. el.  \n");
#endif
  int div = 1; //(alter) ? 2 : 1;
  int currdiv = (alter_) ? 1 : 3;

  if (srf_->dimension() == 3 && initial_surface_)
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

  // Initiate with data points
  if (points_.size() > 0)
    LRSplineUtils::distributeDataPoints(srf_.get(), points_, true, 
					LRSplineUtils::REGULAR_POINTS, 
					outlier_detection_);
  if (sign_points_.size() > 0)
    LRSplineUtils::distributeDataPoints(srf_.get(), sign_points_, true, 
					LRSplineUtils::SIGNIFICANT_POINTS, 
					outlier_detection_);

  if (make_ghost_points_ && !initial_surface_ && srf_->dimension() == 1 && 
      !useMBA_)
    {
      // This is experimental code and should, if kept, be integrated
      // with LRSurfSmoothLS::addDataPoints
      vector<double> ghost_points;
      constructGhostPoints(ghost_points);
      LRSplineUtils::distributeDataPoints(srf_.get(), ghost_points, true, 
					  LRSplineUtils::GHOST_POINTS,
					  outlier_detection_);
    }

  if (make_ghost_points_ && initial_surface_)
    {
      // No need to construct ghost points from extrapolation. Use
      // the input surface
      constructInnerGhostPoints();
    }

  // TEST
  if (srf_->dimension() == 1)
    evalsrf_ = shared_ptr<Eval1D3DSurf>(new Eval1D3DSurf(srf_));
  
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
  if (initMBA_)
    {
      runMBAUpdate(false);
      LSapprox.setInitSf(srf_, coef_known_);
      updateCoefKnown();
    }
  else if (!(initial_surface_ && useMBA_))
    {
     LSapprox.setInitSf(srf_, coef_known_);
      LSapprox.updateLocals();
      //performSmooth(&LSapprox);
     if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
     	adaptSurfaceToConstraints();
    }


#ifdef DEBUG_SURF
  std::ofstream of1("init_sf.g2");
  std::ofstream of1el("init_el.g2");
  srf_->writeStandardHeader(of1);
  srf_->write(of1);
  of1 << std::endl;
  // std::ofstream of12("init_tpsf.g2");
  // shared_ptr<SplineSurface> ssf1(tmp->asSplineSurface());
  // ssf1->writeStandardHeader(of12);
  // ssf1->write(of12);
  // of12 << std::endl;
  LineCloud lines = srf_->getElementBds();
  lines.writeStandardHeader(of1el);
  lines.write(of1el);
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
      std::cout << "Maximum distance exceeding tolerance (dist-tol): " << maxout_ << std::endl;
      std::cout << "Average distance exceeding tolerance (dist-tol): " << avout_ << std::endl;
      std::cout << "Number of signicant points: " << nmb_sign_ << std::endl;
      std::cout << "Maximum distance, significant points: " << maxdist_sign_ << std::endl;
      std::cout << "Average distance, significant points: " << avdist_sign_ << std::endl;
      std::cout << "Number of significant points outside tolerance(" << sign_aepsge_ << "): " << outsideeps_sign_ << std::endl;
    }

#ifdef DEBUG_DIST
  std::ofstream r_out("residual0.txt");
  int dim = srf_->dimension();
  int del2 = dim + 3; // Parameter pair, point and distance
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (!it->second->hasDataPoints())
	continue;
      int del = it->second->getNmbValPrPoint();
      if (del == 0)
	del = del2;
      vector<double>& points = it->second->getDataPoints();
      vector<double>& sign_points = it->second->getSignificantPoints();
      int nmb_pts = (int)points.size()/del;
      int nmb_sign = (int)sign_points.size()/del;
      int nmb_all = nmb_pts + nmb_sign;
      double *curr;
      int ka, kb;
      for (ka=0, curr=&points[0]; ka<nmb_all; ++ka)
	{
	  for (kb=0; kb<del2; ++kb)
	    r_out << curr[kb] << " ";
	  r_out << std::endl;
	  if (ka == nmb_pts-1 && nmb_all > nmb_pts)
	    curr = &sign_points[0];
	  else
	    curr += del;
	}
    }
#endif

  double Tny = 5.0;
  if (compute_AIC_)
    {
      vector<double> residual;
      residual.reserve(points_.size()+sign_points_.size());
      int dim = srf_->dimension();
      int del2 = dim + 3; // Parameter pair, point and distance
      for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
	   it != srf_->elementsEnd(); ++it)
	{
	  if (!it->second->hasDataPoints())
	    continue;
	  int del = it->second->getNmbValPrPoint();
	  if (del == 0)
	    del = del2;
	  vector<double>& points = it->second->getDataPoints();
	  vector<double>& sign_points = it->second->getSignificantPoints();
	  int nmb_pts = (int)points.size()/del;
	  int nmb_sign = (int)sign_points.size()/del;
	  int nmb_all = nmb_pts + nmb_sign;
	  double *curr;
	  int ka, kb;
	  for (ka=0, curr=&points[0]; ka<nmb_all; ++ka)
	    {
	      residual.push_back(curr[del2-1]);
	      if (ka == nmb_pts-1 && nmb_all > nmb_pts)
		curr = &sign_points[0];
	      else
		curr += del;
	    }
	}

      // loglh is log likelyhood with the simplified and loglh2 with the standard
      // Mahalanobis distance using student t-distribution. logn and logn2 are the log
      // likelyhoods with close to normal distribution
      // Only log likelyhood and AIC with standard Mahalanobis distance are stored.
      double loglh2 = 0, logn2 = 0;
      double loglh = LogLikelyhood::compute(residual, Tny, true, loglh2);
      double logn = LogLikelyhood::compute(residual, 50.0, false, logn2);
      int ncoef = srf_->numBasisFunctions();
      double AIC = 2.0*(ncoef - loglh);
      double AIC2 = 2.0*(ncoef - loglh2);
      double AICn = 2.0*(ncoef - logn);
      double AICn2 = 2.0*(ncoef - logn2);
      AIC_.push_back(loglh2);
      AIC_.push_back(AIC2);
      AIC_.push_back(logn2);
      AIC_.push_back(AICn2);
      ncoef_.push_back(ncoef);
    }
  
  if (write_feature_)
    {
      std::ofstream f_out("cellinfo0.txt");
      LRFeatureUtils::writeCellInfo(*srf_, aepsge_, ncell_, f_out);
    }

#ifdef DEBUG_REFINE
  fprintf(fp,"  0 \t %9.3f \t %9.3f \t %9d \t %9d \t %9d \t %9.3f\n", maxdist_, avdist_all_, outsideeps_, srf_->numBasisFunctions(), nmb_pts_-outsideeps_, avdist_);
#endif

  ghost_elems.clear();
  points_.clear();  // Not used anymore TESTING
  double threshold_prev = -1.0;
  int prevcoef = srf_->numBasisFunctions();
  int prevelem = srf_->numElements();
  double av_prev = avdist_all_;
  double max_prev = maxdist_;
  int outsideeps_prev = outsideeps_;
  int ki;
  for (ki=0; ki<max_iter; ++ki)
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

      if (true)//ki > 0 || (!initial_surface_))
	{
	  double threshold;
	  if (threshold1_ == 1 || threshold1_ == 4)
	    {
	      threshold = (aepsge_ + 0.5*maxout_)/2.0 + avout_;
	      if (ki > 0 && maxdist_/maxdist_prev_ > 0.9)
		{
		  // Slow convergence. Reduce threshold
		  threshold = (aepsge_ + avdist_ + 0.5*maxdist_)/3.0;
		}
	      if (threshold_prev > 0.0 && threshold/threshold_prev > 0.9)
		threshold = 0.9*threshold_prev;
	      threshold = std::max(aepsge_, threshold);
	    }
	  else
	    threshold = aepsge_;
#ifdef DEBUG_REFINE
	  std::cout << "Threshold: " << threshold << std::endl;
#endif
	  int nmb_refs;
	  if (category1_ <= 4)
	    nmb_refs = refineSurf3(ki+1, currdiv, threshold);
	  else if (category1_ == 5)
	    nmb_refs = refineSurf4(currdiv, threshold);
	  else
	    nmb_refs = refineSurf(ki+1, currdiv, threshold);
	  if (nmb_refs == 0)
	    {
#ifdef DEBUG_REFINE	      
	      std::cout << "No refinements performed" << std::endl;
#endif
	      break;  // No refinements performed
	    }
	  
	  // if (points_.size() > 0)
	  //   LRSplineUtils::distributeDataPoints(srf_.get(), points_, true, 
	  // 					LRSplineUtils::REGULAR_POINTS, 
	  // 					outlier_detection_);
	  // if (sign_points_.size() > 0)
	  //   LRSplineUtils::distributeDataPoints(srf_.get(), sign_points_, true, 
	  // 					LRSplineUtils::SIGNIFICANT_POINTS, 
	  // 					outlier_detection_);

	  threshold_prev = threshold;
	}
#ifdef DEBUG
      std::ofstream of2("refined_sf.g2");
      std::ofstream of2el("refined_el.g2");
      srf_->writeStandardHeader(of2);
      srf_->write(of2);
      of2 << std::endl;
      LineCloud lines2 = srf_->getElementBds();
      lines2.writeStandardHeader(of2el);
      lines2.write(of2el);

      std::ofstream of3("point_clouds.g2");
      // distance between surface and point
      int dim = srf_->dimension();
      for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
	   it != srf_->elementsEnd(); ++it)
	{
	  vector<double>& elem_data = it->second->getDataPoints();
	  // Number of doubles for each point
	  int del = it->second->getNmbValPrPoint();
	  if (del == 0)
	    del = dim+3;  // Parameter pair, point and distance
	  int back = (del > dim+3) ? 2 : 1;
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
				 elem_data.begin()+(kr+1)*del-back);
		  
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

      // Check if continued least squared method if feasible
      if (!useMBA_)
	{
	  if (srf_->numBasisFunctions() >= maxLScoef_)
	    useMBA_ = true;  // Left side matrix too large
	}
  
//       // Check for linear independence (overloading)
//       vector<LRBSpline2D*> funs = LinDepUtils::unpeelableBasisFunctions(*srf_);
// #ifdef DEBUG
//       std::cout << "Number of unpeelable functions: " << funs.size() << std::endl;
// #endif
     
      // Construct additional ghost points in elements with a low
      // distribution of points
      bool ghost_points_inner = true;
      //if (false /*make_ghost_points_ && ki>0*/)
      if (make_ghost_points_ && ki>0 && ghost_points_inner &&
	  !useMBA_ && ki<toMBA_)
	{
#ifdef DEBUG
	  std::cout << "Ghost points " << std::endl;
#endif
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
	//#ifdef DEBUG
	std::cout << "Using MBA" << std::endl;
	//#endif
	runMBAUpdate(true);
      }
      else
	{
	  //#ifdef DEBUG
	  std::cout << "Using LS" << std::endl;
	  //#endif
	  try {
	    LSapprox.updateLocals();
	    performSmooth(&LSapprox);
	    if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
	      adaptSurfaceToConstraints();
	  }
	  catch (...)
	    {
#ifdef DEBUG
	      std::cout << "Switch to MBA" << std::endl;
#endif
	      useMBA_ = true;
	      runMBAUpdate(true);
	    }
	}

#ifdef DEBUG_SURF
      char tmp[4];
      int ki2 = ki + 1;
      sprintf(tmp, "_%d", ki2);
      char filename1[40];
      strcpy(filename1, "updated_sf");
      char filename2[40];
      strcpy(filename2, "updated_el");
      char filename3[40];
      strcpy(filename3, "updated_coef");
      strncat(filename1, tmp, 3);
      strncat(filename1, ".g2", 3);
      strncat(filename2, tmp, 3);
      strncat(filename2, ".g2", 3);
      strncat(filename3, tmp, 3);
      strncat(filename3, ".g2", 3);
      std::ofstream of4(filename1);
      std::ofstream of4el(filename2);
      std::ofstream of4coef(filename3);
      srf_->writeStandardHeader(of4);
      srf_->write(of4);
      of4 << std::endl;
      LineCloud lines3 = srf_->getElementBds();
      lines3.writeStandardHeader(of4el);
      lines3.write(of4el);
      of4coef << "400 1 0 4 100 0 155 255" << std::endl;
      of4coef << srf_->numBasisFunctions() << std::endl;
      LRSplineSurface::BSplineMap::const_iterator it1 = 
	srf_->basisFunctionsBegin();
      for (; it1 != srf_->basisFunctionsEnd(); ++it1)
	{
	  Point coef = it1->second->Coef();
	  if (srf_->dimension() == 1)
	    {
	      Point greville = it1->second->getGrevilleParameter();
	      of4coef << greville << " " << coef << std::endl;
	    }
	  else
	    of4coef << coef << std::endl;
	}
      // std::ofstream of42("updated_tpsf.g2");
      // shared_ptr<SplineSurface> ssf4(tmp3->asSplineSurface());
      // ssf4->writeStandardHeader(of42);
      // ssf4->write(of42);
      // of42 << std::endl;
      // LineCloud lines32 = tmp3->getElementBds();
      // lines32.writeStandardHeader(of42);
      // lines32.write(of42);
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
      if (ki == max_iter-1)
	outlier_detection_ = false;  // Last accuracy check

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
	  std::cout << "Maximum distance exceeding tolerance (dist-tol): " << maxout_ << std::endl;
	  std::cout << "Average distance exceeding tolerance (dist-tol): " << avout_ << std::endl;
	  std::cout << "Number of coefficients: " << srf_->numBasisFunctions() << std::endl;
	  std::cout << "Number of registered outliers: " << nmb_outliers_ << std::endl;
	  std::cout << "Number of signicant points: " << nmb_sign_ << std::endl;
	  std::cout << "Maximum distance, significant points: " << maxdist_sign_ << std::endl;
	  std::cout << "Average distance, significant points: " << avdist_sign_ << std::endl;
	  std::cout << "Number of significant points outside tolerance(" << sign_aepsge_ << "): " << outsideeps_sign_ << std::endl;
	}

      int nmb_div_el = 0, nmb_none = 0, nmb_under = 0;
      int nmbcoef = srf_->numBasisFunctions();
      int numelem = srf_->numElements();
      if (nmbcoef == prevcoef)
	{
#ifdef DEBUG
	  std::cout << "No added degrees of freedom" << std::endl;
#endif
	  break;
	}

      int max_inner = 0;
      double av_inner = 0.0;
      for (LRSplineSurface::BSplineMap::const_iterator it1=srf_->basisFunctionsBegin();
	   it1 != srf_->basisFunctionsEnd(); ++it1)
	{
	  LRBSpline2D* curr = it1->second.get();

	  int n1 = curr->suppMax(XFIXED) - curr->suppMin(XFIXED) - 1;
	  int n2 = curr->suppMax(YFIXED) - curr->suppMin(YFIXED) - 1;
	  max_inner = std::max(max_inner, std::max(n1, n2));
	  av_inner += (double)(n1 + n2);
	}
      av_inner /= (double)(nmbcoef+nmbcoef);

      double outfac = (double)(outsideeps_prev-outsideeps_)/(double)(nmbcoef-prevcoef);
#ifdef TEST_REFINE
      fprintf(fp,"%3d \t %9.3f \t %9.3f \t %9d \t %9d \t %9.3f \t %9d \t %9.3f \t %9.3f \t %9.3f \t %9d \t %9.3f \t %9.3f \t %9d \n",
	      (ki+1)/div, maxdist_, avdist_all_, outsideeps_, nmbcoef,
	      outfac, nmb_pts_-outsideeps_,
	      (double)(nmb_pts_-outsideeps_)/(double)nmbcoef,
	      (double)nmb_none/(double)nmb_div_el, (double)nmb_under/(double)nmb_div_el,
	      max_inner, av_inner,
	      // outsideeps_prev-outsideeps_,
	      // nmbcoef-prevcoef, numelem-prevelem,
	      // max_prev-maxdist_, av_prev-avdist_all_,
	      avdist_, numelem);
#endif
      prevcoef = nmbcoef;
      prevelem = numelem;
      outsideeps_prev = outsideeps_;
      max_prev = maxdist_;
      av_prev = avdist_all_;
      if (outfac < swap_)
	{
	  category1_ = category2_;
	  threshold1_ = threshold2_;
	}
      
#ifdef DEBUG_DIST
      std::string body = "residual";
      std::string extension = ".txt";
      std::string ver = std::to_string(ki+1);
      std::string outfile = body + ver + extension;
      std::ofstream r_out2(outfile.c_str());
      int dim = srf_->dimension();
      int del2 = dim + 3; // Parameter pair, point and distance
      for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
	   it != srf_->elementsEnd(); ++it)
	{
	  if (!it->second->hasDataPoints())
	    continue;
	  int del = it->second->getNmbValPrPoint();
	  if (del == 0)
	    del = del2;
	  vector<double>& points = it->second->getDataPoints();
	  vector<double>& sign_points = it->second->getSignificantPoints();
	  int nmb_pts = (int)points.size()/del;
	  int nmb_sign = (int)sign_points.size()/del;
	  int nmb_all = nmb_pts + nmb_sign;
	  double *curr;
	  int ka, kb;
	  for (ka=0, curr=&points[0]; ka<nmb_all; ++ka)
	    {
	      for (kb=0; kb<del2; ++kb)
		r_out2 << curr[kb] << " ";
	      r_out2 << std::endl;
	      if (ka == nmb_pts-1 && nmb_all > nmb_pts)
		curr = &sign_points[0];
	      else
		curr += del;
	    }
	}
#endif

  if (compute_AIC_)
    {
      vector<double> residual;
      residual.reserve(points_.size()+sign_points_.size());
      int dim = srf_->dimension();
      int del2 = dim + 3; // Parameter pair, point and distance
      for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
	   it != srf_->elementsEnd(); ++it)
	{
	  if (!it->second->hasDataPoints())
	    continue;
	  int del = it->second->getNmbValPrPoint();
	  if (del == 0)
	    del = del2;
	  vector<double>& points = it->second->getDataPoints();
	  vector<double>& sign_points = it->second->getSignificantPoints();
	  int nmb_pts = (int)points.size()/del;
	  int nmb_sign = (int)sign_points.size()/del;
	  int nmb_all = nmb_pts + nmb_sign;
	  double *curr;
	  int ka, kb;
	  for (ka=0, curr=&points[0]; ka<nmb_all; ++ka)
	    {
	      residual.push_back(curr[del2-1]);
	      if (ka == nmb_pts-1 && nmb_all > nmb_pts)
		curr = &sign_points[0];
	      else
		curr += del;
	    }
	}

      // loglh is log likelyhood with the simplified and loglh2 with the standard
      // Mahalanobis distance using student t-distribution. logn and logn2 are the log
      // likelyhoods with close to normal distribution
      // Only log likelyhood and AIC with standard Mahalanobis distance are stored.
      double loglh2 = 0, logn2 = 0;
      double loglh = LogLikelyhood::compute(residual, Tny, true, loglh2);
      double logn = LogLikelyhood::compute(residual, 50.0, false, logn2);
      int ncoef = srf_->numBasisFunctions();
      double AIC = 2.0*(ncoef - loglh);
      double AIC2 = 2.0*(ncoef - loglh2);
      double AICn = 2.0*(ncoef - logn);
      double AICn2 = 2.0*(ncoef - logn2);
      AIC_.push_back(loglh2);
      AIC_.push_back(AIC2);
      AIC_.push_back(logn2);
      AIC_.push_back(AICn2);
      ncoef_.push_back(ncoef);
    }
  
  if (write_feature_)
        {
         std::string body = "cellinfo";
         std::string extension = ".txt";
         std::string ver = std::to_string(ki+1);
         std::string outfile = body + ver + extension;
         std::ofstream f_out2(outfile.c_str());
	 LRFeatureUtils::writeCellInfo(*srf_, aepsge_, ncell_, f_out2);
	}

    }

  if (nmb_sign_ > 0 && maxdist_sign_ > sign_aepsge_)
    {
      // Perform a last MBA iteration with increased weight on
      // significant points
      double sign_fac2 = 10.0;
      if (srf_->dimension() == 3)
	{
	  LRSplineMBA::MBADistAndUpdate(srf_.get(), 
					sign_fac2*significant_fac_,
					mba_sgn_);
	}
      else if (omp_for_mba_update)
	{
	  LRSplineMBA::MBAUpdate_omp(srf_.get(), 
				     sign_fac2*significant_fac_, mba_sgn_);
	}
      else
	{
	  LRSplineMBA::MBAUpdate(srf_.get(), 
				 sign_fac2*significant_fac_, mba_sgn_);
	}

      // Update distance information
       if (omp_for_elements)
	computeAccuracy_omp(ghost_elems);
      else
	computeAccuracy(ghost_elems);

      if (verbose_)
	{
	  std::cout << std::endl;
	  std::cout << "Significant points update. Maximum distance: " << maxdist_;
	  std::cout << ", average distance: " << avdist_all_ << std::endl;
	  std::cout << "Number of points outside tolerance: " << outsideeps_;
	  std::cout << ", average distance in outside points: " << avdist_ << std::endl;
	  std::cout << "Maximum distance exceeding tolerance (dist-tol): " << maxout_ << std::endl;
	  std::cout << "Average distance exceeding tolerance (dist-tol): " << avout_ << std::endl;
	  std::cout << "Number of coefficients: " << srf_->numBasisFunctions() << std::endl;
	  std::cout << "Number of registered outliers: " << nmb_outliers_ << std::endl;
	  std::cout << "Number of signicant points: " << nmb_sign_ << std::endl;
	  std::cout << "Maximum distance, significant points: " << maxdist_sign_ << std::endl;
	  std::cout << "Average distance, significant points: " << avdist_sign_ << std::endl;
	  std::cout << "Number of significant points outside tolerance(" << sign_aepsge_ << "): " << outsideeps_sign_ << std::endl;
	}
   }
#ifdef TEST_REFINE
  std::cout << "Level: " << ki+1 << std::endl;
  if (ki%div == 1)
    {
      int nmbcoef = srf_->numBasisFunctions();
      int numelem = srf_->numElements();
      fprintf(fp,"%3.1f \t %9.3f \t %9.3f \t %9d \t %9d \t %9.3f \t %9d \t %9d \t %9d \t %9d \t %9.3f \t %9.3f \t %9.3f \t %9d \t %9d \n",
	      (double)ki/(double)div, maxdist_, avdist_all_, outsideeps_, nmbcoef,
	      (double)(outsideeps_prev-outsideeps_)/(double)(nmbcoef-prevcoef),
	      nmb_pts_-outsideeps_, outsideeps_prev-outsideeps_,
	      nmbcoef-prevcoef, numelem-prevelem,
	      max_prev-maxdist_, av_prev-avdist_all_, avdist_, numelem);
    }
#endif
  

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
#ifdef TEST_REFINE
  fclose(fp);
#endif
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
    LSapprox->setLeastSquares_omp(approx_weight, significant_fac_);
  }
  else
  {
      LSapprox->setLeastSquares(approx_weight, significant_fac_);
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
  maxout_ = 0.0;
  avout_ = 0.0;
  outsideeps_sign_ = 0;
  maxdist_sign_ = 0.0;
  avdist_sign_ = 0.0;

  double distfac = (maxdist_prev_ > 0) ? avdist_all_prev_/maxdist_prev_ : 1.0;
  double threshfac = (distfac < 0.1) ? 0.75 : 0.5;
  double outlier_threshold = (1.0-threshfac)*maxdist_prev_ + 
    threshfac*avdist_all_prev_;
  outlier_threshold = std::max(outlier_threshold, 5.0*aepsge_);
  double outlier_fac = 0.5;
  bool update_global = false;
  int nmb_outliers = 0;

  int outlierK = 100;  
  double dom_size = (srf_->endparam_u()-srf_->startparam_u())*
    (srf_->endparam_v()-srf_->startparam_v());
  double density = dom_size/nmb_pts_;
  double outlier_rad = 0.5*sqrt(outlierK*density);

#ifdef _OPENMP
  const bool omp_for_element_pts = false;//true;
#else
  const bool omp_for_element_pts = false;
#endif
  //std::cout << "Omp for elements: " << omp_for_element_pts << std::endl;
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

  LRSplineSurface::ElementMap::const_iterator it;
  int num = srf_->numElements();
  int kj;

  double ghost_fac = 0.8;
  ghost_elems.clear();


  //for (it=srf_->elementsBegin(), kj=0; it != srf_->elementsEnd(); ++it, ++kj)
  for (it=srf_->elementsBegin(), kj=0; kj<num; ++it, ++kj)
    {
      if (!(it->second->hasDataPoints() || it->second->hasSignificantPoints()))
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
      vector<double>& sign_points = it->second->getSignificantPoints();
      vector<double>& ghost_points = it->second->getGhostPoints();
      int nmb_pts = it->second->nmbDataPoints();
      int nmb_ghost = it->second->nmbGhostPoints();
      int del = it->second->getNmbValPrPoint();
      if (del == 0)
	del = dim+3;  // Parameter pair, point and distance
      int nmb_sign = (int)sign_points.size()/del;
      int nmb_all = nmb_pts + nmb_sign;


       // Local error information
      double max_err = 0.0;
      double av_err = 0.0;
      double acc_err = 0.0;
      double acc_outside = 0.0;
      int outside = 0;
      int outside_sign = 0;
      double acc_err_sgn = 0.0;
      double av_err_sgn = 0.0;

      // Check if the accuracy can have been changed
      const vector<LRBSpline2D*>& bsplines = it->second->getSupport();
      size_t nb;
      for (nb=0; nb<bsplines.size(); ++nb)
	if (!bsplines[nb]->coefFixed())
	  break;

      vector<double> prev_point_dist(nmb_pts, 0.0);
      vector<double> prev_sign_dist(nmb_sign, 0.0);
      vector<double> prev_ghost_dist(nmb_ghost, 0.0);
      if (true/*useMBA_ || nb < bsplines.size()*/)
	{
	  // Compute distances in data points and update parameter pairs
	  // if requested
// #ifdef _OPENMP
// 	    double time0_part = omp_get_wtime();
// #endif
	  if (nmb_pts > 0)
	  {
	      if (omp_for_element_pts)
		  computeAccuracyElement_omp(points, nmb_pts, del, rd, 
					     it->second.get(), prev_point_dist);
	      else
		  computeAccuracyElement(points, nmb_pts, del, rd, 
					 it->second.get(), prev_point_dist);
	  }
	  
	  // Compute distances in significant points
	  if (nmb_sign > 0)
	  {
	      if (omp_for_element_pts)
		  computeAccuracyElement_omp(sign_points, nmb_sign, del, rd, 
					     it->second.get(), prev_sign_dist);
	      else
		  computeAccuracyElement(sign_points, nmb_sign, del, rd, 
					 it->second.get(), prev_sign_dist);
	  }
	  
	  // Compute distances in ghost points
	  if (nmb_ghost > 0 && !useMBA_)
	  {
	      if (omp_for_element_pts)
		  computeAccuracyElement_omp(ghost_points, nmb_ghost, del, 
					     rd, it->second.get(), 
					     prev_ghost_dist);
	      else
		  computeAccuracyElement(ghost_points, nmb_ghost, del, rd, 
					 it->second.get(), prev_ghost_dist);
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
      std::ofstream ofp("curr_el_points.g2");
      std::ofstream ofg("curr_ghost_points.g2");
      ofp << "400 1 0 4 0 100 155 255" << std::endl;
      ofp << nmb_pts << std::endl;
      for (int ka=0; ka<nmb_pts; ++ka)
	{
	  for (int kb=0; kb<3; ++kb)
	    ofp << points[ka*del+kb] << "  ";
	  ofp << std::endl;
	}

      ofg << "400 1 0 4 155 100 0 255" << std::endl;
      ofg << nmb_ghost << std::endl;
      for (int ka=0; ka<nmb_ghost; ++ka)
	{
	  for (int kb=0; kb<3; ++kb)
	    ofg << ghost_points[ka*del+kb] << "  ";
	  ofg << std::endl;
	}      
#endif
      int ix = (del > dim+3) ? del-2 : del-1;
      double max_err_prev = 0.0;
      double acc_err_prev = 0.0;
      int nmb_pts2 = 0;  // Number of points excluding outliers
      double minheight = std::numeric_limits<double>::max();
      double maxheight = std::numeric_limits<double>::lowest();
      double tol;  // Local tolerance taking the possibility for a tolerance 
      // threshold varying with a high distant from zero into account
      for (ki=0, curr=(nmb_pts>0) ? &points[0] : &sign_points[0]; ki<nmb_all;)
	{
	  Point curr_pt(curr+(dim==3)*2, curr+ix);
	  bool outlier  = (del > dim+3 && curr[del-1] < 0.0);
	  if (!outlier)
	    {
	      // Height limits
	      double height = curr[ix-1];
	      minheight = std::min(minheight, height);
	      maxheight = std::max(maxheight, height);

	      // Accumulate approximation error
	      double dist2 = fabs(curr[ix]);
	      maxdist_ = std::max(maxdist_, dist2);
	      max_err = std::max(max_err, dist2);
	      acc_err += dist2;
	      acc_err_sgn += curr[ix];
	      avdist_all_ += dist2;
	      double dist3 = 0.0; 
	      if (ki < nmb_pts)
		{
		  ++nmb_pts2;
		  dist3 = fabs(prev_point_dist[ki]);
		}
	      max_err_prev = std::max(max_err_prev, dist3);
	      acc_err_prev += dist3;
	      if (ki < nmb_pts)
		{
		  tol = aepsge_;
		  for (size_t kr=0; kr<tolerances_.size(); ++kr)
		    {
		      if (tolerances_[kr].contains(curr[0], curr[1]))
			{
			  tol = tolerances_[kr].tol;
			  break;
			}
		    }
		}
	      else
		tol = sign_aepsge_;
	      if (has_var_tol_ && (ki<nmb_pts || has_var_tol_sign_))
		{
		  tol = (height < 0.0) ? tol - var_fac_neg_*height :
		    tol + var_fac_pos_*height;
		  tol = std::max(tol, (ki<nmb_pts) ? mintol_ :
				 mintol_*(sign_aepsge_/aepsge_));
		  maxout_ = std::max(maxout_, dist2-tol);
		}
	      else
		maxout_ = std::max(maxdist_-aepsge_, 
				   maxdist_sign_-sign_aepsge_);
	      //if (dist2 > aepsge_)
	      if (ki >= nmb_pts)
		{
		  maxdist_sign_ = std::max(maxdist_sign_, dist2);
		  avdist_sign_ += dist2;
		}
	      if (dist2 > tol)
		{
		  av_err_sgn += curr[ix];
		  avdist_ += dist2;
		  outsideeps_++;
		  av_err += dist2;
		  outside++;
		  avout_ += (dist2-tol);
		  acc_outside += (dist2-tol);
		  if (ki >= nmb_pts)
		    {
		      outside_sign++;
		      outsideeps_sign_++;
		    }

#ifdef DEBUG
		  // Accumulate error points
		  if (curr[ix] > 0)
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
		  if (curr[ix] > 0)
		    ok1.push_back(curr_pt);
		  else
		    ok2.push_back(curr_pt);	
#endif
		}	     	  
	      
	      if (dim == 3 && repar_ && ki<nmb_pts)
		{
		  // Check if the point has moved
		  // Do not parameter iterate significant points
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
		      if (ki == nmb_pts-1 && nmb_all > nmb_pts)
			curr = &sign_points[0];
		      else
			curr += del;
		      ki++;
		    }
		}
	      else
		{
		  if (ki == nmb_pts-1 && nmb_all > nmb_pts)
		    curr = &sign_points[0];
		  else
		    curr += del;
		  ki++;
		}
	    }
	  else
	    {
	      if (ki == nmb_pts-1 && nmb_all > nmb_pts)
		curr = &sign_points[0];
	      else
		curr += del;
	      ki++;
	    }
	}
      if (outside > 0)
	{
	  av_err /= (double)outside;
	  av_err_sgn /= (double)outside;
	}

      nmb_outliers += (nmb_pts - nmb_pts2);

      // Previous accuracy information
      double av_prev, max_prev;
      int nmb_out_prev;
      int nmb_out_sign_prev;
      //double acc_prev = it->second->getAccumulatedError();
      it->second->getAccuracyInfo(av_prev, max_prev, nmb_out_prev, 
				  nmb_out_sign_prev);

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

      if (max_err > aepsge_ && max_prev > 0.0 && 
	  (max_err > ghost_fac*max_prev || acc_err > acc_err_prev) &&
	  nmb_ghost > 0 /*0.25*nmb_pts*/)
	{
	  // Collect element for update of ghost points
	  ghost_elems.push_back(it->second.get());
	}

      if (outside > 0 && outlier_detection_ && max_prev > 0.0 &&
	  max_err > outlier_threshold && max_err > ghost_fac*max_err_prev &&
	  max_err > outlier_fac*acc_err/(double)nmb_pts2)
	{
	  int found = defineOutlierPts(it->second.get(), prev_point_dist, 
				       outlier_threshold, outlier_rad);
	  if (found > 0)
	    {
#ifdef DEBUG1
	      std::ofstream ofp2("curr_el_points2.g2");
	      ofp2.precision(15);
	      ofp2 << "400 1 0 4 0 100 155 255" << std::endl;
	      ofp2 << nmb_pts << std::endl;
	      for (int ka=0; ka<nmb_pts; ++ka)
		{
		  for (int kb=0; kb<3; ++kb)
		    ofp2 << points[ka*del+kb] << "  ";
		  ofp2 << std::endl;
		}
#endif

	      // Recompute local statistics
	      nmb_outliers += found;
	      max_err = 0.0;
	      acc_err = 0.0;
	      av_err = 0.0;
	      outside = 0;
	      outside_sign = 0;
	      minheight = std::numeric_limits<double>::max();
	      maxheight = std::numeric_limits<double>::lowest();
	      for (ki=0, curr=(nmb_pts>0) ? &points[0] : &sign_points[0]; 
		   ki<nmb_all; ++ki)
		{
		  Point curr_pt(curr+(dim==3)*2, curr+ix);
		  bool outlier  = (del > dim+3 && curr[del-1] < 0.0);
		  if (!outlier)
		    {
		      // Height limits
		      double height = curr[ix-1];
		      minheight = std::min(minheight, height);
		      maxheight = std::max(maxheight, height);

		      // Accumulate approximation error
		      double dist2 = fabs(curr[ix]);
		      max_err = std::max(max_err, dist2);
		      acc_err += dist2;
		      if (ki < nmb_pts)
			{
			  tol = aepsge_;
			  for (size_t kr=0; kr<tolerances_.size(); ++kr)
			    {
			      if (tolerances_[kr].contains(curr[0], curr[1]))
				{
				  tol = tolerances_[kr].tol;
				  break;
				}
			    }
			}
		      else
			tol = sign_aepsge_;
		      if (has_var_tol_ && (ki<nmb_pts || has_var_tol_sign_))
			{
			  tol = (height < 0.0) ? tol - var_fac_neg_*height :
			    tol + var_fac_pos_*height;
			  tol = std::max(tol, (ki<nmb_pts) ? mintol_ :
					 mintol_*(sign_aepsge_/aepsge_));
			}
		      if (dist2 > tol)
			{
			  av_err += dist2;
			  outside++;	
			}
		    }
		  if (ki == nmb_pts-1 && nmb_all > nmb_pts)
		    curr = &sign_points[0];
		  else
		    curr += del;
		}
	      if (outside > 0)
		av_err /= (double)outside;
	      update_global = true;
	    }
	}

      // Store updated accuracy information in the element
      it->second->setAccuracyInfo(acc_err, av_err, max_err, outside, 
				  outside_sign, acc_outside);
      it->second->setHeightInfo(minheight, maxheight);
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

  avdist_all_ /= (double)(nmb_pts_ - nmb_outliers);
  if (nmb_sign_ > 0)
    avdist_sign_ /= (double)nmb_sign_;
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

  if (update_global)
    {
      maxdist_ = 0.0;
      avdist_ = 0.0;
      avdist_all_ = 0.0;
      outsideeps_ = 0;
      for (it=srf_->elementsBegin(), kj=0; kj<num; ++it, ++kj)
	{
	  if (!(it->second->hasDataPoints() || 
		it->second->hasSignificantPoints()))
	    continue;
	  
	  double av_dist, max_dist;
	  int nmb_out, nmb_out_sign;
	  it->second->getAccuracyInfo(av_dist, max_dist, nmb_out, nmb_out_sign);

	  maxdist_ = std::max(maxdist_, max_dist);
	  avdist_ += (av_dist*nmb_out);
	  avdist_all_ += it->second->getAccumulatedError();
	  outsideeps_ += nmb_out;
	}
      avdist_ /= (double)outsideeps_;
      avdist_all_ /= (double)(nmb_pts_-nmb_outliers);
    }

  nmb_outliers_ = nmb_outliers;
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
  maxout_ = 0.0;
  avout_ = 0.0;
  outsideeps_sign_ = 0;
  maxdist_sign_ = 0.0;
  avdist_sign_ = 0.0;

  double distfac = (maxdist_prev_ > 0) ? avdist_all_prev_/maxdist_prev_ : 1.0;
  double threshfac = (distfac < 0.1) ? 0.75 : 0.5;
  double outlier_threshold = (1.0-threshfac)*maxdist_prev_ + 
    threshfac*avdist_all_prev_;
  outlier_threshold = std::max(outlier_threshold, 5.0*aepsge_);
  double outlier_fac = 0.5;
  bool update_global = false;
  int nmb_outliers = 0;

  int outlierK = 100;  
  double dom_size = (srf_->endparam_u()-srf_->startparam_u())*
    (srf_->endparam_v()-srf_->startparam_v());
  double density = dom_size/nmb_pts_;
  double outlier_rad = 0.5*sqrt(outlierK*density);


  RectDomain rd = srf_->containingDomain();
  int dim = srf_->dimension();
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

#pragma omp parallel default(none) private(kj, it) shared(dim, elem_iters, rd, ghost_fac, ghost_elems, outlier_threshold, outlier_fac, outlier_rad, update_global, nmb_outliers, num_elem)
  {
      double av_prev, max_prev;
      int nmb_out_prev;
      int nmb_out_sign_prev;
      double umin, umax, vmin, vmax;
      double max_err;
      double av_err;
      double acc_err;
      int outside;
      int outside_sign;
      double acc_err_sgn;
      double av_err_sgn;
      double acc_outside;
      int nmb_pts;
      int nmb_sign;
      int nmb_all;
      int nmb_ghost;
      size_t nb;
      int ki;
      double *curr;
      double dist2, dist3;
      Element2D *elem;
      double acc_prev;
      int del;
      double minheight, maxheight, height;

#pragma omp for schedule(auto)//guided)//static,8)//runtime)//dynamic,4)
      for (kj = 0; kj < num_elem ; ++kj)
      {
	  it = elem_iters[kj];

	  if (!(it->second->hasDataPoints() || 
		it->second->hasSignificantPoints()))
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
	  vector<double>& sign_points = it->second->getSignificantPoints();
	  vector<double>& ghost_points = it->second->getGhostPoints();
	  nmb_pts = it->second->nmbDataPoints();
	  nmb_ghost = it->second->nmbGhostPoints();
	  del = it->second->getNmbValPrPoint();
	  if (del == 0)
	    del = dim+3;  // Parameter pair, point and distance
	  nmb_sign = (int)sign_points.size()/del;
	  nmb_all = nmb_pts + nmb_sign;

	  // Local error information
	  max_err = 0.0;
	  av_err = 0.0;
	  acc_err = 0.0;
	  outside = 0;
	  outside_sign = 0;
	  acc_err_sgn = 0.0;
	  av_err_sgn = 0.0;
	  acc_outside = 0.0;
	  minheight = std::numeric_limits<double>::max();
	  maxheight = std::numeric_limits<double>::lowest();

	  // Check if the accuracy can have been changed
	  const vector<LRBSpline2D*>& bsplines = it->second->getSupport();
	  for (nb=0; nb<bsplines.size(); ++nb)
	      if (!bsplines[nb]->coefFixed())
		  break;

	  vector<double> prev_point_dist(nmb_pts, 0.0);
	  vector<double> prev_sign_dist(nmb_sign, 0.0);
	  vector<double> prev_ghost_dist(nmb_ghost, 0.0);
	  if (true/*useMBA_ || nb < bsplines.size()*/)
	  {
	      // Compute distances in data points and update parameter pairs
	      // if requested
// #ifdef _OPENMP
// 	    double time0_part = omp_get_wtime();
// #endif
	      if (nmb_pts > 0)
	      {
		  computeAccuracyElement(points, nmb_pts, del, rd, 
					 it->second.get(), prev_point_dist);
	      }
	  
	      // Compute distances in significant points
	      if (nmb_sign > 0)
		{
		  computeAccuracyElement(sign_points, nmb_sign, del, rd, 
					 it->second.get(), prev_sign_dist);
		}

	      // Compute distances in ghost points
	      if (nmb_ghost > 0 && !useMBA_)
	      {
		  computeAccuracyElement(ghost_points, nmb_ghost, del, rd, 
					 it->second.get(), prev_ghost_dist);
	      }
// #ifdef _OPENMP
// 	    double time1_part = omp_get_wtime();
// 	    time_computeAccuracyElement += time1_part - time0_part;
// #endif
	  }

	  // Accumulate error information related to data points
	  int ix = (del > dim+3) ? del-2 : del-1;
	  double max_err_prev = 0.0;
	  double acc_err_prev = 0.0;
	  int nmb_pts2 = 0;  // Number of points excluding outliers
	  double tol;  // Local tolerance taking the possibility for a tolerance 
	  // threshold varying with a high distant from zero into account
	  for (ki=0, curr=(nmb_pts>0) ? &points[0] : &sign_points[0]; 
	       ki<nmb_all;)
	  {
	      Point curr_pt(curr+(dim==3)*2, curr+ix);

	      bool outlier  = (del > dim+3 && curr[del-1] < 0.0);
	      if (!outlier)
		{

		  // Height limits
		  height = curr[ix-1];
		  minheight = std::min(minheight, height);
		  maxheight = std::max(maxheight, height);

		  // Accumulate approximation error
		  dist2 = fabs(curr[ix]);
		  maxdist_ = std::max(maxdist_, dist2);
		  max_err = std::max(max_err, dist2);
		  acc_err += dist2;
		  acc_err_sgn += curr[ix];
		  avdist_all_ += dist2;
		  dist3 = 0.0;
		  if (ki < nmb_pts)
		    {
		      ++nmb_pts2;
		      dist3 = fabs(prev_point_dist[ki]);
		    }
		  max_err_prev = std::max(max_err_prev, dist3);
		  acc_err_prev += dist3;
		  if (ki < nmb_pts)
		    {
		      tol = aepsge_;
		      for (size_t kr=0; kr<tolerances_.size(); ++kr)
			{
			  if (tolerances_[kr].contains(curr[0], curr[1]))
			    {
			      tol = tolerances_[kr].tol;
			      break;
			    }
			}
		    }
		  else
		    tol = sign_aepsge_;
		  if (has_var_tol_ && (ki<nmb_pts || has_var_tol_sign_))
		    {
		      tol = (height < 0.0) ? tol - var_fac_neg_*height :
			tol + var_fac_pos_*height;
		      tol = std::max(tol, (ki<nmb_pts) ? mintol_ :
				     mintol_*(sign_aepsge_/aepsge_));
		      maxout_ = std::max(maxout_, dist2-tol);
		    }
		  else
		    maxout_ = maxdist_ - aepsge_;
		  //if (dist2 > aepsge_)
	      //if (dist2 > aepsge_)
		  if (ki >= nmb_pts)
		    {
		      maxdist_sign_ = std::max(maxdist_sign_, dist2);
		      avdist_sign_ += dist2;
		    }

		  if (dist2 > tol)
		    {
		      av_err_sgn += curr[ix];
		      avdist_ += dist2;
		      outsideeps_++;
		      av_err += dist2;
		      outside++;
		      avout_ += (dist2-tol);
		      acc_outside += (dist2-tol);
		      if (ki >= nmb_pts)
			{
			  outside_sign++;
			  outsideeps_sign_++;
			}
		    }
		  else
		    {
		    }	     	  

		  if (dim == 3 && repar_ && ki<nmb_pts)
		    {
		      // Check if the point has moved
		      // Do not parameter iterate significant points
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
			  if (ki == nmb_pts-1 && nmb_all > nmb_pts)
			    curr = &sign_points[0];
			  else
			    curr += del;
			  ki++;
			}
		    }
		  else
		    {
		      if (ki == nmb_pts-1 && nmb_all > nmb_pts)
			curr = &sign_points[0];
		      else
			curr += del;
		      ki++;
		    }
		}
	      else
		{
		  if (ki == nmb_pts-1 && nmb_all > nmb_pts)
		    curr = &sign_points[0];
		  else
		    curr += del;
		  ki++;
		}
	  }
	  if (outside > 0)
	    {
	      av_err /= (double)outside;
	      av_err_sgn /= (double)outside;
	    }

	  nmb_outliers += (nmb_pts - nmb_pts2);

	  // Previous accuracy information
	  acc_prev = it->second->getAccumulatedError();
	  it->second->getAccuracyInfo(av_prev, max_prev, nmb_out_prev, 
				      nmb_out_sign_prev);

	  if (max_err > aepsge_ && max_prev > 0.0 && max_err > ghost_fac*max_prev &&
	      nmb_ghost > 0.25*nmb_pts)
	  {
	      // Collect element for update of ghost points
#pragma omp critical
	      ghost_elems.push_back(it->second.get());
	  }

	  if (outside > 0 && outlier_detection_ && max_prev > 0.0 &&
	      max_err > outlier_threshold && max_err > ghost_fac*max_err_prev &&
	      max_err > outlier_fac*acc_err/(double)nmb_pts2)
	    {
	      int found = defineOutlierPts(it->second.get(), prev_point_dist, 
					   outlier_threshold, outlier_rad);
	      if (found > 0)
		{
#ifdef DEBUG1
		  std::ofstream ofp2("curr_el_points2.g2");
		  ofp2.precision(15);
		  ofp2 << "400 1 0 4 0 100 155 255" << std::endl;
		  ofp2 << nmb_pts << std::endl;
		  for (int ka=0; ka<nmb_pts; ++ka)
		    {
		      for (int kb=0; kb<3; ++kb)
			ofp2 << points[ka*del+kb] << "  ";
		      ofp2 << std::endl;
		    }
#endif
		  // Recompute local statistics
		  nmb_outliers += found;
		  max_err = 0.0;
		  acc_err = 0.0;
		  av_err = 0.0;
		  outside = 0;
		  outside_sign = 0;
		  minheight = std::numeric_limits<double>::max();
		  maxheight = std::numeric_limits<double>::lowest();
		  for (ki=0, curr=(nmb_pts>0) ? &points[0] : &sign_points[0]; 
		       ki<nmb_all; ++ki)
		    {
		      Point curr_pt(curr+(dim==3)*2, curr+ix);
		      bool outlier  = (del > dim+3 && curr[del-1] < 0.0);
		      if (!outlier)
			{
			  // Height limits
			  double height = curr[ix-1];
			  minheight = std::min(minheight, height);
			  maxheight = std::max(maxheight, height);

			  // Accumulate approximation error
			  double dist2 = fabs(curr[ix]);
			  max_err = std::max(max_err, dist2);
			  acc_err += dist2;
			  if (ki < nmb_pts)
			    {
			      tol = aepsge_;
			      for (size_t kr=0; kr<tolerances_.size(); ++kr)
				{
				  if (tolerances_[kr].contains(curr[0], curr[1]))
				    {
				      tol = tolerances_[kr].tol;
				      break;
				    }
				}
			    }
			  else
			    tol = sign_aepsge_;
			  if (has_var_tol_ && (ki<nmb_pts || has_var_tol_sign_))
			    {
			      tol = (height < 0.0) ? tol - var_fac_neg_*height :
				tol + var_fac_pos_*height;
			      tol = std::max(tol, (ki<nmb_pts) ? mintol_ :
					     mintol_*(sign_aepsge_/aepsge_));
			    }
			  if (dist2 > tol)
			    {
			      av_err += dist2;
			      outside++;	
			    }
			}
		    }
		  if (outside > 0)
		    av_err /= (double)outside;
		  update_global = true;
		}
	    }

	  // Store updated accuracy information in the element
	  it->second->setAccuracyInfo(acc_err, av_err, max_err, outside, 
				      outside_sign, acc_outside);
	  it->second->setHeightInfo(minheight, maxheight);

      }
  }

 avdist_all_ /= (double)(nmb_pts_ - nmb_outliers);
 if (nmb_sign_ > 0)
   avdist_sign_ /= (double)nmb_sign_;
 if (outsideeps_ > 0)
   {
     avdist_ /= (double)outsideeps_;
     avout_ /= (double)outsideeps_;
   }

  if (update_global)
    {
      maxdist_ = 0.0;
      avdist_ = 0.0;
      avdist_all_ = 0.0;
      outsideeps_ = 0;
      for (it=srf_->elementsBegin(), kj=0; kj<num; ++it, ++kj)
	{
	  if (!(it->second->hasDataPoints() || 
		it->second->hasSignificantPoints()))
	    continue;
	  
	  double av_dist, max_dist;
	  int nmb_out, nmb_out_sign;
	  it->second->getAccuracyInfo(av_dist, max_dist, nmb_out, 
				      nmb_out_sign);

	  maxdist_ = std::max(maxdist_, max_dist);
	  avdist_ += (av_dist*nmb_out);
	  avdist_all_ += it->second->getAccumulatedError();
	  outsideeps_ += nmb_out;
	}
      avdist_ /= (double)outsideeps_;
      avdist_all_ /= (double)(nmb_pts_-nmb_outliers);
    }

  nmb_outliers_ = nmb_outliers;
// #ifdef _OPENMP
//   double time1 = omp_get_wtime();
//   double time_spent = time1 - time0;
//   std::cout << "time_spent in computeAccuracy: " << time_spent << std::endl;
//   std::cout << "time_spent in computeAccuracyElement: " << time_computeAccuracyElement << std::endl;
// #endif
}

//==============================================================================
  void LRSurfApprox::computeAccuracyElement(vector<double>& points, int nmb, int del,
					    RectDomain& rd, const Element2D* elem,
					    vector<double>& prev_point_dist)
//==============================================================================
{
  int ki, kj, kr;
  double tol = 1.0e-12;  // Numeric tolerance
  double *curr;
  int dim = srf_->dimension();
  // double umin = srf_->paramMin(XFIXED);
  double umax = srf_->paramMax(XFIXED);
  // double vmin = srf_->paramMin(YFIXED);
  double vmax = srf_->paramMax(YFIXED);
  int maxiter = 3; //4;
  Element2D* elem2 = (Element2D*)elem;

  // Fetch basis functions
  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  const int nmb_bsplines = (int)bsplines.size();
  double bval;
  long double sfval;

  vector<double> grid_height;
  double elem_grid_start[2];
  int grid1=0, grid2=0, grid3=0, grid4=0;
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
      if (grid_start_[0] + grid2*cell_size_[0] <= elmax_u)
	grid2++;
      if (grid_start_[1] + grid4*cell_size_[1] <= elmax_v)
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

  int del2 = (del > dim+3) ? del-1 : del;
    for (ki=0, curr=&points[0]; ki<nmb; ++ki, curr+=del)
    {
      curr_pt = Point(curr+(dim==3)*2, curr+del2-1);
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
	    dist = curr[del2-2]-close_pt[2];
	}
      else
	{
	  if (grid_)
	    {
	      // Identify grid cell
	      idx1 = (int)((curr[0] - elem_grid_start[0])/cell_size_[0]);
	      idx2 = (int)((curr[1] - elem_grid_start[1])/cell_size_[1]);
	      idx1 = std::min(grid2-grid1, idx1);
	      idx2 = std::min(grid4-grid3, idx2);
	      
	      // Check distance in grid corners
	      dist1 = curr[2]-grid_height[idx2*(grid2-grid1+1)+idx1];
	      dist2 = curr[2]-grid_height[idx2*(grid2-grid1+1)+idx1+1];
	      dist3 = 
		curr[2]-grid_height[(idx2+1)*(grid2-grid1+1)+idx1];
	      dist4 = 
		curr[2]-grid_height[(idx2+1)*(grid2-grid1+1)+idx1+1];
	      
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
		  // srf_->point(pos, curr[0], curr[1]/*, elem*/);
		  bool u_at_end = (curr[0] > umax-tol) ? true : false;
		  bool v_at_end = (curr[1] > vmax-tol) ? true : false;
		  vector<double> bval;
		  LRSplineUtils::evalAllBSplines(bsplines, curr[0], curr[1],
		  				   u_at_end, v_at_end, bval);
		  // // vector<Point> bpos;
		  // // LRSplineUtils::evalAllBSplinePos(bsplines, curr[0], curr[1],
		  // // 				   u_at_end, v_at_end, bpos);
		  sfval = 0.0;
		  for (kr=0; kr<nmb_bsplines; ++kr)
		    {
		  //     // bsplines[kr]->evalpos(curr[0], curr[1], &bval);
		  //     // sfval += bval;
		      const Point& tmp_pt = bsplines[kr]->coefTimesGamma();
		      long double tmpval = bval[kr]*tmp_pt[0];
		      sfval += tmpval;
		      // sfval += bpos[kr][0];
		    }
	      
		  dist = curr[2] - sfval;
		  // dist = curr[2] - pos[0];
#if 0
		  // TEST
		  if (fabs(dist) > aepsge_)
		    {
		      Point clo_pt;
		      double clo_u, clo_v, clo_dist;
		      double seed[2];
		      seed[0] = curr[0];
		      seed[1] = curr[1];
		      Point pt(curr[0], curr[1], curr[2]);
		      srf_->setCurrentElement(elem2);
		      evalsrf_->closestPoint(pt, clo_u, clo_v, clo_pt,
					     clo_dist, aepsge_, 1, seed);
		      if (clo_dist < fabs(dist))
			dist = (dist < 0.0) ? -clo_dist : clo_dist;
		      int stop_break = 1;
		    }
#endif
		}
	      else
		{
		  Point pos;
		  srf_->point(pos, curr[0], curr[1], elem2);
		  dist = pos.dist(Point(curr+2, curr+del2));
		  vec = curr_pt - pos;
		  // Point norm;
		  srf_->normal(norm, curr[0], curr[1], elem2);
		  if (vec*norm < 0.0)
		    dist *= -1;
		}
	    }
	}
      prev_point_dist[ki] = curr[del2-1];
      curr[del2-1] = dist;
  }

}


//==============================================================================
void LRSurfApprox::computeAccuracyElement_omp(vector<double>& points, int nmb, int del,
					      RectDomain& rd, const Element2D* elem,
					      vector<double>& prev_point_dist)
//==============================================================================
{
  int ki, kj, kr;
  double *curr;
  int dim = srf_->dimension();
  // double umin = srf_->paramMin(XFIXED);
  double umax = srf_->paramMax(XFIXED);
  // double vmin = srf_->paramMin(YFIXED);
  double vmax = srf_->paramMax(YFIXED);
  int maxiter = 3; //4;
  Element2D* elem2 = (Element2D*)elem;
  int del2 = (del > dim+3) ? del-1 : del;

  // Fetch basis functions
  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  const int nmb_bsplines = (int)bsplines.size();
  //double bval, sfval;

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
      if (grid_start_[0] + grid2*cell_size_[0] <= elmax_u)
	grid2++;
      if (grid_start_[1] + grid4*cell_size_[1] <= elmax_v)
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
#pragma omp parallel default(none) private(ki, curr, idx1, idx2, dist, upar, vpar, close_pt, curr_pt, vec, norm, dist1, dist2, dist3, dist4, sgn, pos, kr, kj/*, sfval, bval*/) \
  shared(points, nmb, umax, vmax, del, dim, rd, maxiter, elem_grid_start, grid2, grid1, grid_height, grid3, grid4, elem2, bsplines, del2, prev_point_dist, nmb_bsplines)
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
      curr_pt = Point(curr+(dim==3)*2, curr+del2-1);
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
	    dist = curr[del2-2]-close_pt[2];
	}
      else
	{
	  if (grid_)
	    {
	      // Identify grid cell
	      idx1 = (int)((curr[0] - elem_grid_start[0])/cell_size_[0]);
	      idx2 = (int)((curr[1] - elem_grid_start[1])/cell_size_[1]);
	      idx1 = std::min(grid2-grid1, idx1);
	      idx2 = std::min(grid4-grid3, idx2);
	      
	      // Check distance in grid corners
	      dist1 = curr[2]-grid_height[idx2*(grid2-grid1+1)+idx1];
	      dist2 = curr[2]-grid_height[idx2*(grid2-grid1+1)+idx1+1];
	      dist3 = 
		curr[2]-grid_height[(idx2+1)*(grid2-grid1+1)+idx1];
	      dist4 = 
		curr[2]-grid_height[(idx2+1)*(grid2-grid1+1)+idx1+1];
	      
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
		  double sfval = 0.0;
		  bool u_at_end = (curr[0] >= umax) ? true : false;
		  bool v_at_end = (curr[1] >= vmax) ? true : false;
		  vector<double> bval;
		  LRSplineUtils::evalAllBSplines(bsplines, curr[0], curr[1],
		  				   u_at_end, v_at_end, bval);
		  sfval = 0.0;
		  for (kr=0; kr<nmb_bsplines; ++kr)
		    {
		      // bsplines[kr]->evalpos(curr[0], curr[1], &bval);
		      // sfval += bval;
		      const Point& tmp_pt = bsplines[kr]->coefTimesGamma();
		      sfval += bval[kr]*tmp_pt[0];

		    }
	      
		  dist = curr[2] - sfval;
		  //dist = curr[2] - pos[0];
#if 0
		  // TEST
		  if (fabs(dist) > aepsge_)
		    {
		      Point clo_pt;
		      double clo_u, clo_v, clo_dist;
		      double seed[2];
		      seed[0] = curr[0];
		      seed[1] = curr[1];
		      Point pt(curr[0], curr[1], curr[2]);
		      srf_->setCurrentElement(elem2);
		      evalsrf_->closestPoint(pt, clo_u, clo_v, clo_pt,
					     clo_dist, aepsge_, 1, seed);
		      if (clo_dist < fabs(dist))
			  dist = (dist < 0.0) ? -clo_dist : clo_dist;
		      int stop_break = 1;
		    }
#endif
		}
	      else
		{
		  Point pos;
		  srf_->point(pos, curr[0], curr[1], elem2);
		  dist = pos.dist(Point(curr+2, curr+del2));
		  vec = curr_pt - pos;
		  // Point norm;
		  srf_->normal(norm, curr[0], curr[1]);
		  if (vec*norm < 0.0)
		    dist *= -1;
		}
	    }
	}
      prev_point_dist[ki] = curr[del2-1];
      curr[del2-1] = dist;
    }

}

//==============================================================================
void  LRSurfApprox::runMBAUpdate(bool computed_accuracy)
//==============================================================================
{
#ifdef _OPENMP
    const bool omp_for_mba_update = true;
#else
    const bool omp_for_mba_update = true;//false; // 201503 The omp version seems to be faster even when run sequentially.
#endif

  double sign_scale = 1.0/std::max(1.0, (double)nmb_mba_iter_);
  if (omp_for_mba_update && (!computed_accuracy) && srf_->dimension() == 1)
    {
      LRSplineMBA::MBADistAndUpdate_omp(srf_.get(), 
					sign_scale*significant_fac_, 
					mba_sgn_);
    }
  else if ((!computed_accuracy) || srf_->dimension() == 3)
    {
      LRSplineMBA::MBADistAndUpdate(srf_.get(), 
				    sign_scale*significant_fac_, mba_sgn_);
    }
  else if (omp_for_mba_update)
    {
      LRSplineMBA::MBAUpdate_omp(srf_.get(), 
				 sign_scale*significant_fac_, mba_sgn_);
    }
  else
    {
      LRSplineMBA::MBAUpdate(srf_.get(), 
			     sign_scale*significant_fac_, mba_sgn_);
    }
  
  if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
    adaptSurfaceToConstraints();

  for (int mba_iter=1; mba_iter<nmb_mba_iter_; ++mba_iter)
    {
      if (omp_for_mba_update && srf_->dimension() == 1)
	{
	  LRSplineMBA::MBADistAndUpdate_omp(srf_.get(), 
					    sign_scale*significant_fac_, 
					    mba_sgn_);
	}
      else
	{
	  LRSplineMBA::MBADistAndUpdate(srf_.get(), 
					sign_scale*significant_fac_,
					mba_sgn_);
	}
      if (has_min_constraint_ || has_max_constraint_ || has_local_constraint_)
	adaptSurfaceToConstraints();
    }
}

//==============================================================================
int LRSurfApprox::defineOutlierPts(Element2D* element, 
				   vector<double>& prev_dist, double lim,
				   double rad)
//==============================================================================
{
  vector<double>& points = element->getDataPoints();
  int nmb = element->nmbDataPoints();
  int del = element->getNmbValPrPoint();

  if (del != 5)
    return 0;  // 3D case or no outlier information possible
  int ix = 3;   // Only elevation (1D)
  int nmb_found = 0;

#ifdef DEBUG2
  std::ofstream of1("el_pts.g2");
  std::ofstream of2("cand_pts.g2");
  std::ofstream of3("outlier_pts.g2");
  std::ofstream of4("domain_pts.g2");
  std::ofstream of5("domain2_pts.g2");
  of1.precision(15);
  of2.precision(15);
  of3.precision(15);
  of4.precision(15);
  of5.precision(15);
  of1 << "400 1 0 4 0 55 200 255" << std::endl;
  of1 << nmb << std::endl;
#endif

  // Traverse point cloud and look for outlier candidates
  int ki;
  vector<int> icand;
  int icand_max = 100;  // If more outlier candidates, postpone search 
  // to next iteration level
  double *curr;
  double fac1 = 0.8;
  double mind0 = std::numeric_limits<double>::max();
  double maxd0 = std::numeric_limits<double>::lowest();
  double minh0 = std::numeric_limits<double>::max();
  double maxh0 = std::numeric_limits<double>::lowest();
  Element2D *elem2 = element;
  double phi_lim = M_PI/6.0;
  for (ki=0, curr=&points[0]; ki<nmb; ++ki, curr+=del)
    {
      if (curr[ix+1] < 0.0)
	continue;
      mind0 = std::min(mind0, curr[ix]);
      maxd0 = std::max(maxd0, curr[ix]);
      minh0 = std::min(minh0, curr[ix-1]);
      maxh0 = std::max(maxh0, curr[ix-1]);
      double val = fabs(curr[ix]);
      if (val >= lim && val > fac1*fabs(prev_dist[ki]) )
	{
	  // Check if projected distance is compatibel with vertical
	  Point clo_pt;
	  double clo_u, clo_v, clo_dist;
	  double seed[2];
	  seed[0] = curr[0];
	  seed[1] = curr[1];
	  Point pt(curr[0], curr[1], curr[2]);
	  srf_->setCurrentElement(elem2);
	  evalsrf_->closestPoint(pt, clo_u, clo_v, clo_pt,
				 clo_dist, aepsge_, 1, seed);
	  double xydist = sqrt((clo_u-curr[0])*(clo_u-curr[0]) +
			       (clo_v-curr[1])*(clo_v-curr[1]));
	  double phi = atan((fabs(val)-clo_dist)/xydist);
	  double phi2 = atan(clo_dist/xydist);
	  if (clo_dist >= lim || phi < phi_lim)
	    icand.push_back(ki);
	  if ((int)icand.size() > icand_max)
	    break;
	}
#ifdef DEBUG2
      of1 << curr[0] << " " << curr[1] << " " << curr[2] << std::endl;
#endif
    }

  if (icand.size() == 0 && icand.size() > icand_max)
    return 0;
#ifdef DEBUG2
  of2 << "400 1 0 4 55 200 0 255" << std::endl;
  of2 << icand.size() << std::endl;
#endif

  // Fetch elements possibly overlapping the outlier neighbourhoods
  vector<Element2D*> elem_cand;
  elem_cand.push_back(element);
  for (size_t ka=0; ka<icand.size(); ++ka)
    {
      int kb = icand[ka]*del;
      getCandElements(points[kb], points[kb+1], rad, element, elem_cand);
#ifdef DEBUG2
      of2 << points[kb] << " " << points[kb+1] << " " << points[kb+2] << std::endl;
#endif
    }

  double rad2 = rad*rad;
  vector<double> avdist(icand.size(), 0.0);
  vector<double> avdist2(icand.size(), 0.0);
  vector<double> avdist3(icand.size(), 0.0);
  vector<double> avsgn(icand.size(), 0.0);
  vector<double> avsgn2(icand.size(), 0.0);
  vector<double> avsgn3(icand.size(), 0.0);
  vector<int> nmb_pt(icand.size(), 0);
  vector<int> nmb_pt3(icand.size(), 0);
  vector<int> nmb_above(icand.size(), 0);
  vector<int> nmb_below(icand.size(), 0);
  vector<double> avg_above(icand.size(), 0.0);
  vector<double> avg_below(icand.size(), 0.0);
  vector<double> mind(icand.size(), std::numeric_limits<double>::max());
  vector<double> maxd(icand.size(), std::numeric_limits<double>::lowest());
  vector<double> mind2(icand.size(), std::numeric_limits<double>::max());
  vector<double> maxd2(icand.size(), std::numeric_limits<double>::lowest());
  vector<double> mind3(icand.size(), std::numeric_limits<double>::max());
  vector<double> maxd3(icand.size(), std::numeric_limits<double>::lowest());
  vector<double> minh(icand.size(), std::numeric_limits<double>::max());
  vector<double> maxh(icand.size(), std::numeric_limits<double>::lowest());
  vector<vector<double> > pt_dist(icand.size());
  vector<double> std1(icand.size(), 0.0);
  vector<double> std2(icand.size(), 0.0);
  vector<double> std3(icand.size(), 0.0);
  vector<double> close_pt1(del*icand.size(), 0.0);
  vector<double> close_d1(icand.size(), std::numeric_limits<double>::max());
  vector<double> close_pt2(del*icand.size(), 0.0);
  vector<double> close_d2(icand.size(), std::numeric_limits<double>::max());
  for (size_t ka=0; ka<elem_cand.size(); ++ka)
    {
      // Fetch candidate neighbourhood points. For ka == 0, points2 = points
      vector<double>& points2 = elem_cand[ka]->getDataPoints();
      int nmb2 = elem_cand[ka]->nmbDataPoints();

      // For all points within the given radius of each candidate outlier,
      // compute statistics
      for (ki=0, curr=&points2[0]; ki<nmb2; ++ki, curr+=del)
	{
	  if (curr[ix+1] < 0.0)
	    continue;

	  double x1 = curr[0];
	  double y1 = curr[1];
	  for (size_t kb=0; kb<icand.size(); ++kb)
	    {
	      double x0 = points[icand[kb]*del];
	      double y0 = points[icand[kb]*del+1];
	      double d2 = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0);
	      if (d2 > rad2)
		continue;

	      double d3 = sqrt(d2);
	      double val = curr[ix];
	      avdist[kb] += fabs(val);
	      avsgn[kb] += val;
	      nmb_pt[kb]++;
	      pt_dist[kb].push_back(val);
	      mind[kb] = std::min(mind[kb], val);
	      maxd[kb] = std::max(maxd[kb], val);
	      minh[kb] = std::min(minh[kb], curr[ix-1]);
	      maxh[kb] = std::max(maxh[kb], curr[ix-1]);
	      if (val > 0.0)
		{
		  nmb_above[kb]++;
		  avg_above[kb] += val;
		}
	      else if (val < 0.0)
		{
		  nmb_below[kb]++;
		  avg_below[kb] += fabs(val);
		}

	      if (!(ka == 0 && icand[kb] == ki))
		{
		  mind2[kb] = std::min(mind2[kb], val);
		  maxd2[kb] = std::max(maxd2[kb], val);
		  avdist2[kb] += fabs(val);
		  avsgn2[kb] += val;
		}

	      if (fabs(val) < lim || 
		  val*points[icand[kb]*del+ix] < 0.0)
		{
		  if (d3 < close_d1[kb])
		    {
		      std::copy(curr, curr+del, close_pt1.begin()+kb*del);
		      close_d1[kb] = d3;
		    }
		  mind3[kb] = std::min(mind3[kb], val);
		  maxd3[kb] = std::max(maxd3[kb], val);
		  avdist3[kb] += fabs(val);
		  avsgn3[kb] += val;
		  nmb_pt3[kb]++;
#ifdef DEBUG2
		  of4 << x1 << " " << y1 << " " << curr[2] << std::endl;
#endif
		}
	      else
		{
		  if ((!(ka == 0 && icand[kb] == ki)) && d3 < close_d2[kb])
		    {
		      std::copy(curr, curr+del, close_pt2.begin()+kb*del);
		      close_d2[kb] = d3;
		    }
#ifdef DEBUG2
		  of5 << x1 << " " << y1 << " " << curr[2] << std::endl;
#endif
		}
	    }
	}
    }
  for (size_t kb=0; kb<icand.size(); ++kb)
    {
      avdist[kb] /= (double)(nmb_pt[kb]);
      avdist2[kb] /= (double)(nmb_pt[kb]-1);
      avsgn[kb] /= (double)(nmb_pt[kb]);
      avsgn2[kb] /= (double)(nmb_pt[kb]-1);
      avdist3[kb] /= (double)(nmb_pt3[kb]);
      avsgn3[kb] /= (double)(nmb_pt3[kb]);
      avg_above[kb] /= (double)(nmb_above[kb]);
      avg_below[kb] /= (double)(nmb_below[kb]);
    }
	
  // Compute standard deviation
  double tmp;
  for (size_t kb=0; kb<icand.size(); ++kb)
    {
      for (size_t kc=0; kc<pt_dist[kb].size(); ++kc)
	{
	  tmp = pt_dist[kb][kc] - avsgn[kb];
	  std1[kb] += tmp*tmp;
	  tmp = pt_dist[kb][kc] - avsgn2[kb];
	  std2[kb] += tmp*tmp;
	  if (fabs(pt_dist[kb][kc]) < lim || 
	      pt_dist[kb][kc]*points[icand[kb]*del+ix] < 0.0)
	    {
	      tmp = pt_dist[kb][kc] - avsgn3[kb];
	      std3[kb] += tmp*tmp;
	    }
	}
      std1[kb] /= (double)nmb_pt[kb];
      tmp = points[icand[kb]*del+ix] - avsgn2[kb];
      std2[kb] -= tmp*tmp;
      std2[kb] /= (double)(nmb_pt[kb]-1);
      std3[kb] /= (double)nmb_pt3[kb];
    }

  double stdfac = 0.9;
  double rfac = 0.75;
  double close_ang = 0.25*M_PI;
  double clfac = 0.25;
  double nmb_fac = 0.8;

  double eps = 2.0*aepsge_;
  for (size_t kb=0; kb<icand.size(); ++kb)
    {
      double rdiv = (maxd3[kb] - mind3[kb])/(maxd[kb]-mind[kb]);
      double mult_fac = (double)nmb_pt3[kb]/(double)nmb_pt[kb];
      double dh1 = close_pt1[kb*del+ix-1]-points[icand[kb]*del+ix-1];
      double phi = atan(fabs(dh1)/close_d1[kb]);
      //double cldiv = (nmb_pt[kb] > nmb_pt3[kb]+1) ? close_d1[kb]/close_d2[kb] : 1.0;
      if (stdfac*mult_fac*std1[kb] > std3[kb] &&
	  mult_fac*avdist[kb] > avdist3[kb] &&
	  rdiv < mult_fac*rfac &&
	  phi > close_ang /* &&  // Makes sense if xy and h are given in 
	  // approximately the same unit
	  cldiv > clfac*/ &&
	  (double)nmb_pt3[kb] > nmb_fac*(double)nmb_pt[kb] &&
	  fabs(points[icand[kb]*del+ix-1] - close_pt1[kb*del+ix-1]) > eps &&
	  fabs(points[icand[kb]*del+ix] - close_pt1[kb*del+ix]) > eps /*&&
									maxh[kb]-minh[kb] > maxd[kb]-mind[kb]*/)
	{
	  points[icand[kb]*del+ix+1] = -1;
	  nmb_found++;
	}
      else
	  points[icand[kb]*del+ix+1] += 1;
	
      }

#ifdef DEBUG2
  if (nmb_found > 0)
    {
      of3 << "400 1 0 4 255 0 0 255" << std::endl;
      of3 << nmb_found << std::endl;
      for (size_t ka=0; ka<icand.size(); ++ka)
	if (points[icand[ka]*del+ix+1] < 0)
	  {
	    double *tmp = &points[icand[ka]*del];
	    of3 << tmp[0] << " " << tmp[1] << " " << tmp[2] << std::endl;
	  }
    }
#endif
  return nmb_found;
}

#if 0
//==============================================================================
bool LRSurfApprox::defineOutlierPts2(vector<double>& points, int nmb, int del,
				    vector<double>& prev_dist, double lim,
				    double max_err, double av_dist, 
				    double density)
//==============================================================================
{
  if (del != 5)
    return false;  // 3D case or no outlier information possible

#ifdef DEBUG2
  std::ofstream of("curr_el_pts.g2");
  of.precision(15);
  of << "400 1 0 4 0 100 155 255" << std::endl;
  of << nmb << std::endl;
  for (int ka=0; ka<nmb; ++ka)
    {
      for (int kb=0; kb<3; ++kb)
	of << points[ka*del+kb] << "  ";
      of << std::endl;
    }
#endif
  // Compute range, signed average and standard deviation for all points
  int ix = 3;   // Only elevation (1D)
  double av_sgn = 0.0;
  double av_sgn_prev = 0.0;
  double *curr, *curr2;
  int ki, kj;
  int nmb2 = 0;  // Number of points when already found outliers are removed
  double maxr1 = std::numeric_limits<double>::lowest();
  double minr1 = std::numeric_limits<double>::max();
  double maxh1 = std::numeric_limits<double>::lowest();
  double minh1 = std::numeric_limits<double>::max();
  for (ki=0, curr=&points[0]; ki<nmb; ++ki, curr+=del)
    {
      if (curr[ix+1] < 0.0)
	continue;
      maxr1 = std::max(maxr1, curr[ix]);
      minr1 = std::min(minr1, curr[ix]);
      maxh1 = std::max(maxh1, curr[ix-1]);
      minh1 = std::min(minh1, curr[ix-1]);
      av_sgn += curr[ix];
      av_sgn_prev += prev_dist[ki];
      nmb2++;
    }
  av_sgn /= (double)nmb2;
  av_sgn_prev /= (double)nmb2;

  if (maxh1-minh1 < maxr1-minr1)
    return false;   // Variation in distances to surface larger than 
  // variation in height

  double std1 = 0.0;
  double tmp;
  for (ki=0, curr=&points[0]; ki<nmb; ++ki, curr+=del)
    {
      if (curr[ix+1] < 0.0)
	continue;
      tmp = curr[ix] - av_sgn;
      std1 += tmp*tmp;
    }
  std1 /= (double)nmb2;

  // Traverse point cloud and look for outlier candidates
  int nmb_outliers = 0;
  double fac1 = 0.8;
  double fac2 = 0.6;
  double fac3 = 0.95;
  if (nmb2 < 20)
    {
      fac1 = 0.9;
      fac2 = 0.4;
    }
  else if (nmb2 > 500)
    fac2 = 0.9;

  double rfac = 0.75;
  double abovefac = 0.05;
  while (true) 
    {
      for (ki=0, curr=&points[0]; ki<nmb; ++ki, curr+=del)
	{
	  if (curr[ix+1] < 0.0)
	    continue;
	  
	  double val = curr[ix];
	  if (fabs(val) > lim)
	    {
	      // Candidate outlier. Test configuration
	      int nmb_above = 0;
	      int nmb_above2 = 0;
	      double max_below = 0.0;
	      double lim2 = 0.5*(fabs(val)+lim);
	      double av_dist2 = (av_dist*nmb2 - fabs(val))/(double)(nmb2-1);
	      double av_sgn2 = (av_sgn*nmb2 - val)/(double)(nmb2-1);
	      double std2 = 0.0;
	      double maxr2 = std::numeric_limits<double>::lowest();
	      double minr2 = std::numeric_limits<double>::max();
	      double maxr3 = std::numeric_limits<double>::lowest();
	      double minr3 = std::numeric_limits<double>::max();
	      double maxh2 = std::numeric_limits<double>::lowest();
	      double minh2 = std::numeric_limits<double>::max();
	      for (kj=0, curr2=&points[0]; kj<nmb; ++kj, curr2+=del)
		{
		  if (kj == ki || curr2[ix+1] < 0.0)
		    continue;
		  if (fabs(curr2[ix]) > lim)
		    nmb_above++;
		  else
		    {
		      max_below = std::max(max_below, fabs(curr2[ix]));
		      maxr3 = std::max(maxr3, curr2[ix]);
		      minr3 = std::min(minr3, curr2[ix]);
		    }
		  if (fabs(curr2[ix]) > lim2)
		    nmb_above2++;
		  maxr2 = std::max(maxr2, curr2[ix]);
		  minr2 = std::min(minr2, curr2[ix]);
		  maxh2 = std::max(maxh2, curr2[ix-1]);
		  minh2 = std::min(minh2, curr2[ix-1]);
		  tmp = curr2[ix] - av_sgn2;
		  std2 += tmp*tmp;
		}
	      std2 /= (double)(nmb2-1);

	      double rdiv = (maxr2 - minr2)/(maxr1 - minr1);
	      double rdiv3 = (maxr3 - minr3)/(maxr1 - minr1);
	      double hdiv = (maxh2 - minh2)/(maxh1 - minh1);
	      double stdfac = 0.98; //(nmb_above == 0) ? fac2 : fac3;
	      if (fabs(val) > fac1*fabs(prev_dist[ki]) &&
		  av_dist > av_dist2 && stdfac*std1 > std2 &&
		  rdiv < ((nmb_above == 0) ? rfac : 1.0) && 
		  rdiv3 < rfac && nmb_above < abovefac*nmb &&
		  rdiv3*density < fabs(val)/max_below) 
		{
#ifdef DEBUG2
		  std::ofstream of2("curr_out.g2");
		  of2.precision(15);
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << 1 << std::endl;
		  for (int kb=0; kb<3; ++kb)
		    of2 << points[ki*del+kb] << "  ";
		  of2 << std::endl;
#endif
		  curr[ix+1] = -1.0;
		  nmb_outliers++;
		  nmb2--;
		  av_dist = av_dist2;
		  av_sgn = av_sgn2;
		  std1 = std2;
		  minr1 = minr2;
		  maxr1 = maxr2;
		  minh1 = minh2;
		  maxh1 = maxh2;
		  break;
		}
	    }
	}
      if (ki == nmb)
	break;
    }

  return (nmb_outliers > 0);
}
#endif

//==============================================================================
bool compare_refs(LRSplineSurface::Refinement2D r1, 
		  LRSplineSurface::Refinement2D r2)
{
  return (r1.kval < r2.kval);
}


//==============================================================================
bool compare_elems(pair<Element2D*,double> el1, pair<Element2D*,double> el2)
{
  return (el1.second > el2.second);
}

//==============================================================================
void LRSurfApprox::getRefineExtension(Element2D *elem, Direction2D fixdir,
				      int strategy, double& ppar, double& pmin, double& pmax,
				      set<pair<Element2D*,pair<Direction2D,double> > >& unbalanced_elem)
//==============================================================================
{
  // Fetch B-splines
  const vector<LRBSpline2D*>& bsplines = elem->getSupport();
  size_t nmb = bsplines.size();

  ppar = (fixdir == XFIXED) ? 0.5*(elem->umin() + elem->umax()) :
    0.5*(elem->vmin() + elem->vmax());
  pmin = (fixdir == XFIXED) ? elem->vmin() : elem->umin();
  pmax = (fixdir == XFIXED) ? elem->vmax() : elem->umax();
  double par2 = 0.5*(pmin + pmax); 
  vector<pair<Element2D*,pair<Direction2D,double> > > unbalanced2;
  if (strategy <= 1 || strategy >=5)
    {
      // All overlapping B-splines
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  double bmin = (fixdir == XFIXED) ? bsplines[ki]->vmin() :
	    bsplines[ki]->umin();
	  double bmax = (fixdir == XFIXED) ? bsplines[ki]->vmax() :
	    bsplines[ki]->umax();
	  pmin = std::min(pmin, bmin);
	  pmax = std::max(pmax, bmax);
	}
    }
  else if (strategy == 2)
    {
      // Largest overlapping B-spline
      double tol = srf_->getKnotTol();
      double max_size = 0.0;
      double min_frac = 0.0;
      int ix = -1;
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  // Compute size of B-spline
	  double bmin = (fixdir == XFIXED) ? bsplines[ki]->vmin() :
	    bsplines[ki]->umin();
	  double bmax = (fixdir == XFIXED) ? bsplines[ki]->vmax() :
	    bsplines[ki]->umax();
	  double bsize = bmax - bmin;
	  double bdel1 = par2 - bmin;
	  double bdel2 = bmax - par2;
	  double frac = std::min(bdel1, bdel2)/std::max(bdel1,bdel2);
	  if ((fabs(max_size-bsize) < tol && frac < min_frac) ||
	      bsize > max_size)
	    {
	      max_size = bsize;
	      min_frac = frac;
	      ix = (int)ki;
	    }
	}
      pmin = (fixdir == XFIXED) ? bsplines[ix]->vmin() :
	bsplines[ix]->umin();
      pmax = (fixdir == XFIXED) ? bsplines[ix]->vmax() :
	bsplines[ix]->umax();
    }
  else if (strategy == 3)
    {
      // "Best" overlapping B-spline
      double tol = 0.1;
      double max_wgt = 0.0;
      double min_frac = 0.0;
      int ix = -1;
     for (size_t ki=0; ki<nmb; ++ki)
	{
	  // Count the number of elements with large error affected
	  double curr_wgt = 0.0;
	  const vector<Element2D*>& curr_el = bsplines[ki]->supportedElements();
	  for (size_t kj=0; kj<curr_el.size(); ++kj)
	    {
	      double emin = (fixdir == XFIXED) ?
		curr_el[kj]->umin() : curr_el[kj]->vmin();
	      double emax = (fixdir == XFIXED) ?
		curr_el[kj]->umax() : curr_el[kj]->vmax();
	      if (emax < ppar || emin > ppar)
		continue;  // Element not affected

	      // Compute weight for importance of refinement
	      double max_err, av_err;
	      int nmb_outside, nmb_out_sign;
	      curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside, 
					   nmb_out_sign);
	      int nmb_pts = curr_el[kj]->nmbDataPoints();
	      if (nmb_pts > 0)
		{
		  double wgt = av_err*(double)nmb_outside/(double)nmb_pts;
		  curr_wgt += wgt;
		}
	    }

	  double bmin = (fixdir == XFIXED) ? bsplines[ki]->vmin() :
	    bsplines[ki]->umin();
	  double bmax = (fixdir == XFIXED) ? bsplines[ki]->vmax() :
	    bsplines[ki]->umax();
	  double bdel1 = par2 - bmin;
	  double bdel2 = bmax - par2;
	  double frac = std::min(bdel1, bdel2)/std::max(bdel1,bdel2);
	  if ((fabs(max_wgt-curr_wgt) < tol && frac < min_frac) ||
	      curr_wgt > max_wgt)
	    {
	      max_wgt = curr_wgt;
	      min_frac = frac;
	      ix = (int)ki;
	    }
	}
      pmin = (fixdir == XFIXED) ? bsplines[ix]->vmin() :
	bsplines[ix]->umin();
      pmax = (fixdir == XFIXED) ? bsplines[ix]->vmax() :
	bsplines[ix]->umax();
    }
  else // if (strategy == 4)
	{
	  // Combination of 2 and 3
      double tol = 0.1;
      double max_wgt = 0.0;
      double max_size = 0.0;
      double min_frac = 0.0;
      double maxfrac_combined = 0.0;
      int ix = -1;
      for (size_t ki=0; ki<nmb; ++ki)
	{
	  // Count the number of elements with large error affected
	  double curr_wgt = 0.0;
	  const vector<Element2D*>& curr_el = bsplines[ki]->supportedElements();
	  vector<pair<Element2D*,pair<Direction2D,double> > > unbalanced;
	  for (size_t kj=0; kj<curr_el.size(); ++kj)
	    {
	      double emin = (fixdir == XFIXED) ?
		curr_el[kj]->umin() : curr_el[kj]->vmin();
	      double emax = (fixdir == XFIXED) ?
		curr_el[kj]->umax() : curr_el[kj]->vmax();
	      if (emax < ppar || emin > ppar)
		continue;  // Element not affected

	      if ((emax - ppar > 0.55*(emax-emin) || ppar - emin > 0.55*(emax-emin)) &&
		  emax-ppar > 0.0001 && ppar-emin > 0.0001)
		{
		  unbalanced.push_back(make_pair(curr_el[kj],make_pair(fixdir,ppar)));
		}


	      // Compute weight for importance of refinement
	      double max_err, av_err;
	      int nmb_outside, nmb_out_sign;
	      curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside, 
					   nmb_out_sign);
	      int nmb_pts = curr_el[kj]->nmbDataPoints();
	      if (nmb_pts > 0)
		{
		  double wgt = av_err*(double)nmb_outside/(double)nmb_pts;
		  curr_wgt += wgt;
		}
	    }

	  double bmin = (fixdir == XFIXED) ? bsplines[ki]->vmin() :
	    bsplines[ki]->umin();
	  double bmax = (fixdir == XFIXED) ? bsplines[ki]->vmax() :
	    bsplines[ki]->umax();
	  double bsize = bmax - bmin;
	  double bdel1 = par2 - bmin;
	  double bdel2 = bmax - par2;
	  double frac = std::min(bdel1, bdel2)/std::max(bdel1,bdel2);
	  double frac_combined = (ix < 0) ? 1.0 : curr_wgt/max_wgt + bsize/max_size;
	  if ((fabs(maxfrac_combined-frac_combined) < tol && frac < min_frac) ||
	      frac_combined > maxfrac_combined)
	    {
	      maxfrac_combined = frac_combined;
	      max_wgt = curr_wgt;
	      max_size = bsize;
	      min_frac = frac;
	      ix = (int)ki;
	      unbalanced2 = unbalanced;
	    }
	}
      unbalanced_elem.insert(unbalanced2.begin(), unbalanced2.end());
#if 0
      if (unbalanced2.size() > 0)
	{
	  if (fixdir == XFIXED)
	    ppar = 0.5*(unbalanced2[0].first->umin() + unbalanced2[0].first->umax());
	  else
	    ppar = 0.5*(unbalanced2[0].first->vmin() + unbalanced2[0].first->vmax());
	}
#endif
      pmin = (fixdir == XFIXED) ? bsplines[ix]->vmin() :
	bsplines[ix]->umin();
      pmax = (fixdir == XFIXED) ? bsplines[ix]->vmax() :
	bsplines[ix]->umax();
   }
}

//==============================================================================
int LRSurfApprox::refineSurf3(int iter, int& dir, double threshold)
//==============================================================================
{
  // Test
  double tol = srf_->getKnotTol();

  // Traverse all B-splines and check for elements with outside points
    for (LRSplineSurface::BSplineMap::const_iterator it1=srf_->basisFunctionsBegin();
       it1 != srf_->basisFunctionsEnd(); ++it1)
    {
      LRBSpline2D* curr = it1->second.get();

      int el_out2 = 0.0;
      bool out_pts = false;
      for (auto it2=curr->supportedElementBegin(); 
	   it2 != curr->supportedElementEnd(); ++it2)
	{
	  int num_out = (*it2)->getNmbOutsideTol();
	  if (num_out > 0)
	    {
	      out_pts = true;
	      break;
	    }
	}
      if (out_pts)
	curr->setFixCoef(0);
      else
	curr->setFixCoef(1);
    }

  int num_elem = srf_->numElements();
  double av_wgt = 0.0;
  int el_out = 0;
  double min_wgt = std::numeric_limits<double>::max();
  double max_wgt = 0.0;
  vector<double> all_wgt(num_elem, 0.0);
  size_t kr=0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it, ++kr)
    {
      double av_err, max_err;
      int nmb_out, nmb_out_sign;
      int nmb_pts = it->second->nmbDataPoints();
      it->second->getAccuracyInfo(av_err, max_err, nmb_out, nmb_out_sign);
      if (nmb_out > 0 || nmb_out_sign > 0)
	{
	  double wgt = nmb_out + 2.0*nmb_out_sign + av_err; //+ max_err 
	  all_wgt[kr] = wgt;
	  av_wgt += wgt;
	  min_wgt = std::min(min_wgt, wgt);
	  max_wgt = std::max(max_wgt, wgt);
	  el_out++;
	}
    }
  av_wgt /= (double)el_out;

  std::sort(all_wgt.begin(), all_wgt.end());
  for (kr=0; kr<all_wgt.size(); ++kr)
    if (all_wgt[kr] > 0.0)
      break;

  double med_wgt2 = all_wgt[((int)kr+num_elem)/2];
  double fac = (max_wgt > 2.0*min_wgt) ? 0.5 : 1.0;;
  double thresh2;
  if (threshold1_ == 2)
    {
      thresh2 = fac*min_wgt + (1.0-fac)*av_wgt;
      // if (thresh2 - floor(thresh2) < 0.5) //0.25)
      //   thresh2 = floor(thresh2);
      double highlim = std::max(min_wgt, 0.9*prev_thresh_);
      thresh2 = std::min(thresh2, highlim); //prev_thresh_);
      prev_thresh_ = thresh2;
    }
  else
    thresh2 = min_wgt;
#ifdef DEBUG_REFINE
  std::cout << "Num elements: " << num_elem << ", elements out: " << el_out << std::endl;
  std::cout << "thresh2 = " << thresh2 << std::endl;
  double thresh3 = (kr < 0.9*num_elem) ? med_wgt2 : min_wgt; //fac*min_wgt + (1.0-fac)*av_wgt; //min_wgt; 
  std::cout << "min_wgt = " << min_wgt << ", av_wgt = " << av_wgt << ", max_wgt = " << max_wgt << std::endl;
  std::cout << "num_elem = " << num_elem << ", first = " << kr << ", med_wgt = " << all_wgt[num_elem/2] << ", med_wgt2 = " << all_wgt[((int)kr+num_elem)/2] << std::endl;
  std::cout << "thresh3 = " << thresh3 << std::endl;
#endif
  double choose_fac1 = 0.75;
  double choose_fac2 = 0.05;
  double outel_fac = (double)el_out/(double)num_elem;
#ifdef DEBUG_REFINE
  std::cout << "outel_fac: " << outel_fac << std::endl;
#endif
  int dir2 = dir;//(outel_fac > choose_fac2 && outel_fac < choose_fac1) ? dir : 3;
  // if (iter == 1)
  //   dir2 = 3;
  if (dir2 == dir && dir2 != 3)
    dir = 3 - dir;

  prev_el_out_= el_out;

#ifdef DEBUG_REFINE
  std::cout << "Dir: " << dir2 << ", category: " << category1_ << std::endl;
#endif
  set<pair<Element2D*,pair<Direction2D,double> > > unbalanced;
  vector<LRSplineSurface::Refinement2D> refs_x, refs_y;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      double av_err, max_err;
      int nmb_out, nmb_out_sign;
      it->second->getAccuracyInfo(av_err, max_err, nmb_out, nmb_out_sign);
      if (max_err < threshold)
	continue;
      double wgt = nmb_out + 2.0*nmb_out_sign + av_err; //+ max_err 
      if (wgt < thresh2)
	continue;
      if (nmb_out > 0 || nmb_out_sign > 0)
	{
	  double umin = it->second->umin();
	  double umax = it->second->umax();
	  double vmin = it->second->vmin();
	  double vmax = it->second->vmax();

	  const vector<LRBSpline2D*>& bsplines = it->second->getSupport();
	  // std::cout << "Bsplines: " << std::endl;
	  // for (size_t kr=0; kr<bsplines.size(); ++kr)
	  //   {
	  //     std::cout << bsplines[kr]->umin() << " "  << bsplines[kr]->umax() << " ";
	  //     std::cout << bsplines[kr]->vmin() << " "  << bsplines[kr]->vmax() << std::endl;
	  //   }
	  
	  if ((dir2 == 1 || dir2 == 3) && umax-umin > 2.0*usize_min_)
	    {
	      double v1, v2, ppar;
	      getRefineExtension(it->second.get(), XFIXED, category1_,  ppar, v1, v2, unbalanced);
	      
	      LRSplineSurface::Refinement2D curr_ref1;
	      curr_ref1.setVal(ppar, v1, v2, XFIXED, 1);
	      // std::cout << "El x: " << umin << " " << umax << " " << vmin << " " << vmax;
	      // std::cout << ". Bsize: " << bsplines.size();
	      // std::cout << ". Ref: " << 0.5*(umin+umax) << " " << v1 << " " << v2 << std::endl;
	      appendRef(refs_x, curr_ref1, tol);
	    }

	  if ((dir2 == 2 || dir2 == 3) && vmax-vmin > 2.0*vsize_min_)
	    {
	      double u1, u2, ppar;
	      getRefineExtension(it->second.get(), YFIXED, category1_,  ppar, u1, u2, unbalanced);
	      
	      LRSplineSurface::Refinement2D curr_ref2;
	      curr_ref2.setVal(ppar, u1, u2, YFIXED, 1);
	      // std::cout << "El y: " << umin << " " << umax << " " << vmin << " " << vmax;
	      // std::cout << ". Bsize: " << bsplines.size();
	      // std::cout << ". Ref: " << 0.5*(vmin+vmax) << " " << u1 << " " << u2 << std::endl;
	      appendRef(refs_y, curr_ref2, tol);
	    }
	}
    }

  // srf_->refine2(refs_x, true);
  // srf_->refine2(refs_y, true);
#ifdef DEBUG
  std::cout << "Refs x: " << refs_x.size() << ", refs_y: " << refs_y.size() << std::endl;
#endif
#if 0
  std::cout << "Possible unbalanced: " << unbalanced.size() << std::endl;
  double fuzzy = 1.0e-10;
  for (auto it=unbalanced.begin(); it!=unbalanced.end(); )
    {
      Direction2D d = (*it).second.first;
      double par = (*it).second.second;
      Element2D* elem = (*it).first;
      double umin = elem->umin();
      double umax = elem->umax();
      double vmin = elem->vmin();
      double vmax = elem->vmax();
      if (d == XFIXED)
	{
	  size_t kr;
	  for (kr=0; kr<refs_x.size(); ++kr)
	    {
	      if (vmin > refs_x[kr].end+fuzzy || vmax < refs_x[kr].start-fuzzy)
		continue;
	      if (fabs(par - 0.5*(umin+umax)) < fuzzy)
		break;
	    }
	  if (kr < refs_x.size())
	    it = unbalanced.erase(it);
	  else
	    ++it;
	}
      else
	{
	  size_t kr;
	  for (kr=0; kr<refs_y.size(); ++kr)
	    {
	      if (vmin > refs_y[kr].end+fuzzy || vmax < refs_y[kr].start-fuzzy)
		continue;
	      if (fabs(par - 0.5*(umin+umax)) < fuzzy)
		break;
	    }
	  if (kr < refs_x.size())
	    it = unbalanced.erase(it);
	  else
	    ++it;
	}
    }
#endif
#ifdef DEBUG
  std::cout << "Unbalanced: " << unbalanced.size() << std::endl;
#endif
 
  for (kr=0; kr<refs_x.size(); ++kr)
    {
      // std::cout << "Ref: " << refs_x[kr].kval << " " << refs_x[kr].start << " ";
      // std::cout << refs_x[kr].end << " " << refs_x[kr].d << std::endl;
      srf_->refine(refs_x[kr], true /*false*/);
    }

  for (size_t kr=0; kr<refs_y.size(); ++kr)
    {
      // std::cout << "Ref: " << refs_y[kr].kval << " " << refs_y[kr].start << " ";
      // std::cout << refs_y[kr].end << " " << refs_y[kr].d << std::endl;
      srf_->refine(refs_y[kr], true /*false*/);
    }
  return (int)refs_x.size() + (int)refs_y.size();
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
int LRSurfApprox::refineSurf4(int& dir, double threshold)
//==============================================================================
{
  double tol = srf_->getKnotTol();
  int dir2 = dir;
  if (dir2 == dir && dir2 != 3)
    dir = 3 - dir;

  // Traverse all B-splines and check for elements with outside points
    for (LRSplineSurface::BSplineMap::const_iterator it1=srf_->basisFunctionsBegin();
       it1 != srf_->basisFunctionsEnd(); ++it1)
    {
      LRBSpline2D* curr = it1->second.get();

      int el_out2 = 0.0;
      bool out_pts = false;
      for (auto it2=curr->supportedElementBegin(); 
	   it2 != curr->supportedElementEnd(); ++it2)
	{
	  int num_out = (*it2)->getNmbOutsideTol();
	  if (num_out > 0)
	    {
	      out_pts = true;
	      break;
	    }
	}
      if (out_pts)
	curr->setFixCoef(0);
      else
	curr->setFixCoef(1);
    }

  vector<LRSplineSurface::Refinement2D> refs_x, refs_y;
  for (LRSplineSurface::BSplineMap::const_iterator it=srf_->basisFunctionsBegin();
       it != srf_->basisFunctionsEnd(); ++it)
    {
      LRBSpline2D* bspline = it->second.get();
      const Mesh2D* mesh = bspline->getMesh();
      int size1 = bspline->degree(XFIXED)+1;
      int size2 = bspline->degree(YFIXED)+1;

      bool refine = false;
      for (auto it2=bspline->supportedElementBegin(); 
	   it2 != bspline->supportedElementEnd(); ++it2)
	{
	  int num_out = (*it2)->getNmbOutsideTol();
	  double av_err, max_err;
	  int nmb_out, nmb_out_sign;
	  (*it2)->getAccuracyInfo(av_err, max_err, nmb_out, nmb_out_sign);
 	  if (num_out > 0 && max_err >= threshold)
	    {
	      refine = true;
	      break;
	    }
	}

      if (refine)
	{
	  // Refine all knot spans
	  if (dir2 == 1 || dir2 == 3)
	    {
	      const vector<int>& kvec = bspline->kvec(XFIXED);
	      for (int kj=0; kj<size1; ++kj)
		{
		  double u1 = mesh->kval(XFIXED, kvec[kj]);
		  double u2 = mesh->kval(XFIXED, kvec[kj+1]);
		  if (u2-u1 <= tol)
		    continue;
		  LRSplineSurface::Refinement2D curr_ref;
		  curr_ref.setVal(0.5*(u1+u2), bspline->vmin(),
				  bspline->vmax(), XFIXED, 1);
		  appendRef(refs_x, curr_ref, tol);
		}
	    }
	  if (dir2 == 2 || dir2 == 3)
	    {
	      const vector<int>& kvec = bspline->kvec(YFIXED);
	      for (int kj=0; kj<size2; ++kj)
		{
		  double v1 = mesh->kval(YFIXED, kvec[kj]);
		  double v2 = mesh->kval(YFIXED, kvec[kj+1]);
		  if (v2-v1 <= tol)
		    continue;
		  LRSplineSurface::Refinement2D curr_ref;
		  curr_ref.setVal(0.5*(v1+v2), bspline->umin(),
				  bspline->umax(), YFIXED, 1);
		  appendRef(refs_y, curr_ref, tol);
		}
	    }
	}
    }
  for (size_t kr=0; kr<refs_x.size(); ++kr)
    {
      srf_->refine(refs_x[kr], true /*false*/);
    }

  for (size_t kr=0; kr<refs_y.size(); ++kr)
    {
      srf_->refine(refs_y[kr], true /*false*/);
    }
  // srf_->refine(refs_x, true);
  // srf_->refine(refs_y, true);
  return (int)refs_x.size() + (int)refs_y.size();
}

//==============================================================================
//int LRSurfApprox::refineSurf(int iter)
int LRSurfApprox::refineSurf(int iter, int& dir, double threshold)
//==============================================================================
{
  //int iter=0;
  
#ifdef DEBUG
  std::ofstream of0("element_info.dat");
  int idx=0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      if (it->second->hasAccuracyInfo())
  	{
  	  double av_err, max_err;
  	  int nmb_out, nmb_out_sign;
  	  it->second->getAccuracyInfo(av_err, max_err, nmb_out, nmb_out_sign);
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
  // Traverse all B-splines and check for elements with outside points
    for (LRSplineSurface::BSplineMap::const_iterator it1=srf_->basisFunctionsBegin();
       it1 != srf_->basisFunctionsEnd(); ++it1)
    {
      LRBSpline2D* curr = it1->second.get();

      int el_out2 = 0.0;
      bool out_pts = false;
      for (auto it2=curr->supportedElementBegin(); 
	   it2 != curr->supportedElementEnd(); ++it2)
	{
	  int num_out = (*it2)->getNmbOutsideTol();
	  if (num_out > 0)
	    {
	      out_pts = true;
	      break;
	    }
	}
      if (out_pts)
	curr->setFixCoef(0);
      else
	curr->setFixCoef(1);
    }

  int num_elem = srf_->numElements();
  int el_out = 0;
  vector<pair<Element2D*, double> > elem_out;
  double av_wgt = 0.0;
  for (LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
       it != srf_->elementsEnd(); ++it)
    {
      double av_err, max_err;
      int nmb_out, nmb_out_sign;
      int nmb_pts = it->second->nmbDataPoints();
      int nmb_sign = it->second->nmbSignificantPoints();
      it->second->getAccuracyInfo(av_err, max_err, nmb_out, nmb_out_sign);
      if (nmb_out > 0 || nmb_out_sign > 0 || max_err > threshold)
	{
	  double wgt = nmb_out + 2.0*nmb_out_sign + max_err + av_err;
	  if (nmb_pts + nmb_sign < 20)
	    wgt /= 2.0;
	  av_wgt += wgt;
	  elem_out.push_back(std::make_pair(it->second.get(), wgt));
	  el_out++;
	}
    }
  av_wgt /= (double)el_out;


  int choice = 1;  // Strategy for knot insertion in one single B-spline
  double choose_fac1 = 0.75;
  double choose_fac2 = 0.05;
  double outel_fac = (double)el_out/(double)num_elem;
#ifdef DEBUG_REFINE
  std::cout << "outel_fac: " << outel_fac << std::endl;
#endif
  int dir2 = dir; //(outel_fac > choose_fac2 && outel_fac < choose_fac1) ? dir : 3;
  // if (iter == 1)
  //   dir2 = 3;
  if (dir2 == dir && dir2 != 3)
    dir = 3 - dir;
  //dir = (dir < 3) ? dir + 1 : 1;
#ifdef DEBUG_REFINE
  std::cout << "Dir: " << dir2 << std::endl;
#endif

  // Construct indexed bspline array and collect related accuracy information
  int group_fac = 3;
  double error_fac = 0.1;
  double error_fac2 = 10.0;
  std::set<Element2D*> elems;
  int num_bspl = srf_->numBasisFunctions();
  vector<LRBSpline2D*> bsplines(num_bspl);
  vector<double> error(num_bspl, 0.0);
  vector<double> max_error(num_bspl, 0.0);
  vector<double> av_error(num_bspl, 0.0);
  vector<int> num_pts(num_bspl, 0);
  vector<int> num_out_pts(num_bspl, 0); 
  vector<int> num_out_sign(num_bspl, 0); 
  vector<double> error2(num_bspl);
  double mean_err = 0.0;
  size_t kr = 0;
  double domsize;
  double average_nmb_out = 0.0;
  double average_nmb = 0.0;
  double basis_average_out = 0.0;
  vector<int> bspl_perm(num_bspl, 0);
  size_t nmb_perm= 0;
  for (LRSplineSurface::BSplineMap::const_iterator it=srf_->basisFunctionsBegin();
       it != srf_->basisFunctionsEnd(); ++it, ++kr)
    {
      LRBSpline2D* curr = it->second.get();

      int el_out2 = 0.0;
      for (auto it2=curr->supportedElementBegin(); 
	   it2 != curr->supportedElementEnd(); ++it2)
	{
	  num_pts[kr] += (*it2)->nmbDataPoints();
	  num_out_pts[kr] += std::max(0, (*it2)->getNmbOutsideTol());
	  num_out_sign[kr] += std::max(0, (*it2)->getNmbSignOutsideTol());
	  //error[kr] += (*it2)->getAccumulatedError();
	  error[kr] += (*it2)->getAccumulatedOutside();
	  max_error[kr] = std::max(max_error[kr], (*it2)->getMaxError());
	  av_error[kr] += (*it2)->getAverageError();  // Only counting those 
	  // points being outside of the tolerance
	  if ((*it2)->getNmbOutsideTol() > 0)
	    el_out2++;
	}
      av_error[kr] /= (double)(curr->nmbSupportedElements());
      bsplines[kr] = curr;
      mean_err += av_error[kr];

      int bnmb_el = curr->nmbSupportedElements();
      basis_average_out += (double)el_out2/(double)bnmb_el;

      // Use sqrt to reduce the significance of this property compared to the
      // error. What would be the effect of instead squaring the error?
      domsize = sqrt((curr->umax()-curr->umin())*(curr->vmax()-curr->vmin()));
      error2[kr] = error[kr]*domsize;
      if (num_out_pts[kr] > group_fac || (double)num_out_pts[kr] > 
	  error_fac*((double)num_pts[kr]))
	error2[kr] *= error_fac2;
      average_nmb_out += (double)(num_out_pts[kr]);
      average_nmb += (double)(num_pts[kr]);
      if (num_out_pts[kr] > 0)
	bspl_perm[nmb_perm++] = kr;
    }
  mean_err /= (double)num_bspl;
  average_nmb_out /= (double)num_bspl;
  average_nmb /= (double)num_bspl;
  basis_average_out /= (double)num_bspl;

  // Sort bsplines according to average error weighted with the domain size
  int ki, kj;
  //vector<int> bspl_perm(num_bspl);
  // for (ki=0; ki<num_bspl; ++ki)
  //   bspl_perm[ki] = ki;

  // Do the sorting
  quicksort(&error2[0], &bspl_perm[0], 0, nmb_perm-1);
  // for (ki=0; ki<num_bspl; ++ki)
  //   {
  //      for (kj=ki+1; kj<num_bspl; ++kj)
  // 	{
  // 	  // Modify if there is a significant number of large error points 
  // 	  if (error2[bspl_perm[ki]] < error2[bspl_perm[kj]])
  // 	    {
  // 	      std::swap(bspl_perm[ki], bspl_perm[kj]);
  // 	    }
  // 	}
  //   }
  
  // Split the most important B-splines, but only if the maximum
  // error is larger than the tolerance
  //double fac = 0.5;
  //int nmb_perm = (int)bspl_perm.size();
  int nmb_split = (int)(0.75*nmb_perm);  //(int)(0.5*nmb_perm);
  //nmb_split = std::min(nmb_split, 600);  // Limit the number of refinements
  int min_nmb_pts = 1; //4;
  //double pnt_fac = 0.2;
  //int min_nmb_out = 4;1

  vector<LRSplineSurface::Refinement2D> refs_x, refs_y;
  int nmb_refs = 0;

  int nmb_fixed = 0;
  //nmb_split = nmb_perm;
  min_nmb_pts = 0;
  double average_threshold = 0.0;
  if (threshold1_ >= 3)
    average_threshold = std::max(0.01*average_nmb, average_nmb_out);
  average_threshold = std::min(average_threshold, 0.9*prev_thresh_);
  prev_thresh_ = average_threshold;
#ifdef DEBUG_REFINE
  std::cout << "Average threshold: " << average_threshold << std::endl;
#endif
  for (kr=0; kr<nmb_perm; ++kr)
    {
      //if (max_error[bspl_perm[kr]] < aepsge_)
      if (num_out_pts[bspl_perm[kr]] == 0 && num_out_sign[bspl_perm[kr]] == 0)
	{
	  // Keep coefficient fixed. Will reduce computation time, but
	  // may lead to less good approximations as adjacent B-splines
	  // overlapping the same points may change
	  bsplines[bspl_perm[kr]]->setFixCoef(1);
	  nmb_fixed++;
	  continue;
	}
      else
	bsplines[bspl_perm[kr]]->setFixCoef(0);  // Adjacent B-splines may 
      // have changed

      // Do not split B-splines with too few points in its domain
      if (num_pts[bspl_perm[kr]] < min_nmb_pts && 
	  num_out_sign[bspl_perm[kr]] == 0)
	continue;

      if (max_error[bspl_perm[kr]] < threshold)
	continue;
      // if (nmb_refs >= nmb_split)
      // 	break;

      if (nmb_refs >= nmb_split && num_out_sign[bspl_perm[kr]] == 0)
	continue;

      // //if (av_error[bspl_perm[kr]] < fac*mean_err)
      // if (false /*av_error[bspl_perm[kr]] < fac*mean_err &&
      // 	  (num_out_pts[bspl_perm[kr]] < (int)(pnt_fac*num_pts[bspl_perm[kr]]) ||
      // 	  num_pts[bspl_perm[kr]] < min_nmb_out)*/)
      // 	continue;  // Do not split this B-spline at this stage

      nmb_refs++;  // Split this B-spline
      
      // How to split					
      vector<Element2D*> elem_div;
      defineRefs(bsplines[bspl_perm[kr]], average_threshold, dir2,
		 refs_x, refs_y, iter, elem_div);
       elems.insert(elem_div.begin(), elem_div.end());
    }
  
  bool elemref = (category1_ == 7);
  if (elemref)
    {
      // Removing affected elements
      std::vector<Element2D*> elems2(elems.begin(), elems.end());
      for (kr=0; kr<elems2.size(); ++kr)
	{
	  size_t kh;
	  for (kh=0; kh<elem_out.size(); ++kh)
	    if (elems2[kr] == elem_out[kh].first)
	      break;
	  if (kh < elem_out.size())
	    elem_out.erase(elem_out.begin() + kh);
	}
 
  // Sort remaining elements
  double frac = 0.6*av_wgt;
  std::sort(elem_out.begin(), elem_out.end(), compare_elems);
  for (kr=0; kr<elem_out.size(); ++kr)
    {
      if (elem_out[kr].second < frac)
	break;   // Not a significant element

      vector<Element2D*> elements;  // Elements affected by the refinement(s)
      checkFeasibleRef(elem_out[kr].first, dir2, iter, refs_x, refs_y, elements);
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

  // // Flag coefficient if not necessary to update
  // for (; kr<bspl_perm.size(); ++kr)
  //   {
  //     //if (max_error[bspl_perm[kr]] < aepsge_)
  //     if (num_out_pts[bspl_perm[kr]] == 0)
  // 	{
  // 	  // TEST. Keep coefficient fixed
  // 	  bsplines[bspl_perm[kr]]->setFixCoef(1);
  // 	  nmb_fixed++;
  // 	}
  //   }
    }
  
#ifdef DEBUG
  std::ofstream of("refine0.dat");
  //std::streamsize prev = of.precision(15);
  (void)of.precision(15);
  for (kr=0; kr<refs_x.size(); ++kr)
    {
      of << refs_x[kr].kval << "  " << refs_x[kr].start << "  " << refs_x[kr].end;
      of << "  " << refs_x[kr].d << "  " << refs_x[kr].multiplicity << std::endl;
    }
  std::cout << "Number of refinements: " << refs_x.size() << std::endl;
  std::cout << "Number of coef fixed: " << nmb_fixed << std::endl;
 #endif

  // srf_->refine(refs_x, true);
  // srf_->refine(refs_y, true);
  // Sort refinements to start from the ends of the surface to minimize
  // number of knot vector indices updates
  //std::sort(refs.begin(), refs.end(), compare_refs);
  //refs_x.clear();
  for (kr=0; kr<refs_x.size(); ++kr)
    {
#ifdef DEBUG
      //      std::cout << "Refine nr " << kr << ": " << refs[kr].kval << "  " << refs[kr].start << "  ";
      //      std::cout << refs[kr].end << "  " << refs[kr].d << "  " << refs[kr].multiplicity << std::endl;
#endif
      // Perform refinements, one at the time to keep information stored in the elements
      srf_->refine(refs_x[kr], true /*false*/);
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
  for (kr=0; kr<refs_y.size(); ++kr)
    {
      srf_->refine(refs_y[kr], true /*false*/);
      int stop_break = 1;
    }

//   #ifdef DEBUG
//   std::ofstream ofmesh("mesh1.eps");
//   writePostscriptMesh(*srf_, ofmesh);
//   #endif

  return (int)refs_x.size() + (int)refs_y.size();
}

#if 0
//==============================================================================
int LRSurfApprox::refineSurf2()
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
  	  int nmb_out, nmb_out_sign;
  	  it->second->getAccuracyInfo(av_err, max_err, nmb_out, nmb_out_sign);
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

  vector<LRSplineSurface::Refinement2D> refs_x, refs_y;
  int nmb_ref = std::max(std::min((int)el_perm.size(), 4), (int)(0.75*num_err));
  
  for (kr=0; kr<el_perm.size(); )
    {
      size_t nmb_perm = el_perm.size();
      if (elem[el_perm[kr]]->getAverageError() < threshhold && 
	  (int)refs_x.size()+(int)refs_y.size() > nmb_ref)
	break;  // No more refinements at the current stage
      //if (elem[el_perm[kr]]->getMaxError() < aepsge_)
      if (elem[el_perm[kr]]->getNmbOutsideTol() == 0)
	{
	  ++kr;
	  continue; //break;
	}

      // Check feasability of split
      //size_t nmb_refs = refs.size();
      vector<Element2D*> elements;  // Elements affected by the refinement(s)
      checkFeasibleRef(elem[el_perm[kr]], 3, 0, refs_x, refs_y, elements);
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
  for (kr=0; kr<refs_x.size(); ++kr)
    {
      of << refs_x[kr].kval << "  " << refs_x[kr].start << "  " << refs_x[kr].end;
      of << "  " << refs_x[kr].d << "  " << refs_x[kr].multiplicity << std::endl;
    }
  for (kr=0; kr<refs_y.size(); ++kr)
    {
      of << refs_y[kr].kval << "  " << refs_y[kr].start << "  " << refs_y[kr].end;
      of << "  " << refs_y[kr].d << "  " << refs_y[kr].multiplicity << std::endl;
    }
#endif

  for (kr=0; kr<refs_x.size(); ++kr)
    {
#ifdef DEBUG
      std::cout << "Refine nr " << kr << ": " << refs_x[kr].kval << "  " << refs_x[kr].start << "  ";
      std::cout << refs_x[kr].end << "  " << refs_x[kr].d << "  " << refs_x[kr].multiplicity << std::endl;
#endif
      // Perform refinements, one at the time to keep information stored in the elements
      srf_->refine(refs_x[kr], true /*false*/);
      // std::ofstream of2("refined2_sf.g2");
      // srf_->writeStandardHeader(of2);
      // srf_->write(of2);
      // of2 << std::endl;
      // int stop_break = 1;
     }

  for (kr=0; kr<refs_y.size(); ++kr)
    {
      srf_->refine(refs_y[kr], true /*false*/);
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
  
  return (int)refs_x.size() + (int)refs_y.size();
}
#endif

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
  int nmb2 = (int)sign_points_.size()/(dim+2);
  shared_ptr<SplineSurface> result_surf;
  if (!initMBA_)
    {
      try {
	result_surf = createSurf(&points_[0], nmb, 
				 dim, ncoef_u, order_u,
				 ncoef_v, order_v, knots_u,
				 knots_v, smoothweight_,
				 maxdist_, avdist_, outsideeps_,
				 (nmb2>0) ? &sign_points_[0] : 0, nmb2);
      }
      catch (...)
	{
	  initMBA_ = true;
	}
    }

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
						   int& nmb_outside,
						   double* points2, int nmb_pts2)
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
  pts.reserve(dim*(nmb_pts+nmb_pts2));
  param.reserve(2*(nmb_pts+nmb_pts2));
  int ki, kj;
  double* it;
  for (it=points, kj=0; kj<nmb_pts; ++kj) 
    {
      for (ki=0; ki<2; ++ki)
	param.push_back(*it++);
      for (ki=0; ki<dim; ++ki)
	pts.push_back(*it++);
    }
  if (points2 != 0)
    {
      for (it=points2, kj=0; kj<nmb_pts2; ++kj) 
	{
	  for (ki=0; ki<2; ++ki)
	    param.push_back(*it++);
	  for (ki=0; ki<dim; ++ki)
	    pts.push_back(*it++);
	}
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
  for (size_t ki=0; ki<sign_points_.size(); ki+=del)
    {
      umin = std::min(umin, sign_points_[ki]);
      umax = std::max(umax, sign_points_[ki]);
      vmin = std::min(vmin, sign_points_[ki+1]);
      vmax = std::max(vmax, sign_points_[ki+1]);
    }
}

//==============================================================================
void LRSurfApprox::initDefaultParams()
//==============================================================================
{
  nmb_sign_ = 0;
  sign_aepsge_ = aepsge_;
  nmb_outliers_ = 0;
  maxdist_ = std::numeric_limits<double>::lowest();
  maxdist_prev_ = std::numeric_limits<double>::lowest();
  maxdist_sign_ = std::numeric_limits<double>::lowest();
  avdist_ = 0.0;
  avdist_all_ = 0.0;
  avdist_sign_ = 0.0;
  avdist_all_prev_ = 0.0;
  outsideeps_ = 0;
  outsideeps_sign_ = 0;
  maxout_ = std::numeric_limits<double>::lowest();
  avout_ = 0.0;
  smoothweight_ = 1.0e-3;
  significant_fac_ = 5.0;
  maxLScoef_ = 25000; // 46340; Memory management to costly for large surfaces
  smoothbd_ = false;
  fix_corner_ = false;
  to3D_ = -1;
  grid_ = false;
  has_min_constraint_ = false;
  has_max_constraint_ = false;
  has_local_constraint_ = false;
  nmb_mba_iter_ = 2;
  mba_sgn_ = 0;
  outlier_detection_ = false;
  has_var_tol_ = false;
  var_fac_pos_ = 0.0;
  var_fac_neg_ = 0.0;
  mintol_ = 0.01;
  verbose_ = false;

  prev_el_out_ = -1.0;
  prev_thresh_ = std::numeric_limits<double>::max();
  
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 0;
  grid_start_[0] = grid_start_[1] = 0.0;
  cell_size_[0] = cell_size_[1] = 1.0;
  usize_min_ = vsize_min_ = -1;

  fix_boundary_ = false; //true;
  make_ghost_points_ = false;

  category1_ = 1;   // Full span
  category2_ = 0;   // Not applied
  alter_ = 1;       // Refine in both parameter directions
  threshold1_ = 2;   // No threshold
  threshold2_ = 0;   // No threshold
  swap_ = -100.0;   // Do not swap refinement strategy
  
  write_feature_ = false;
  ncell_ = 1;
  compute_AIC_ = false;
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
	  int nmb_pts = 0, nmb_sign = 0;
	  for (size_t ki=0; ki<curr_el.size(); ++ki)
	    {
	      nmb_pts += curr_el[ki]->nmbDataPoints();
	      nmb_sign += curr_el[ki]->nmbSignificantPoints();
	    }
	  if (nmb_pts == 0 && nmb_sign == 0)
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
void LRSurfApprox::defineRefs(LRBSpline2D* bspline, double average_out,
			      int dir,
			      vector<LRSplineSurface::Refinement2D>& refs_x,
			      vector<LRSplineSurface::Refinement2D>& refs_y,
			      int iter,
			      vector<Element2D*>& elem_div)
//==============================================================================
{
  // For each alternative (knot span) in each parameter direction, collect
  // accuracy statistic
  // Compute also average element size
  double tol = srf_->getKnotTol();
  int size1 = bspline->degree(XFIXED)+1;
  int size2 = bspline->degree(YFIXED)+1;
  double scratch[30];
  double *alloc = NULL;
  double *u_info, *v_info, *u_elsize, *v_elsize;
  vector<int> u_inside(size1, 0);
  vector<int> v_inside(size2, 0);
  vector<int> u_outside(size1, 0);
  vector<int> v_outside(size2, 0);
  if (size1+size2 < 15)
    {
      std::fill(scratch, scratch+30, 0.0);
      u_info = scratch;
      v_info = u_info+size1;
      v_elsize = v_info + size2;
      u_elsize = v_elsize + size1;
    }
  else
    {
      alloc = new double[2*(size1+size2)];
      std::fill(alloc, alloc+2*(size1+size2), 0.0);
      u_info = alloc;
      v_info = u_info+size1;
      v_elsize = v_info + size2;
      u_elsize = v_elsize + size1;
    }
  // vector<double> u_info(size1, 0.0);
  // vector<double> v_info(size2, 0.0);
  // vector<double> v_elsize(size1, 0.0);
  // vector<double> u_elsize(size2, 0.0);
  
  const vector<int>& kvec_u = bspline->kvec(XFIXED);
  const vector<int>& kvec_v = bspline->kvec(YFIXED);
  const Mesh2D* mesh = bspline->getMesh();
  
  double av_kdiff_u = 0.0, av_kdiff_v = 0.0;
  for (size_t kj=1; kj<kvec_u.size(); ++kj)
    av_kdiff_u += (mesh->kval(XFIXED, kvec_u[kj]) - 
		   mesh->kval(XFIXED, kvec_u[kj-1]));
  av_kdiff_u /= (double)(kvec_u.size()-1);
  for (size_t kj=1; kj<kvec_v.size(); ++kj)
    av_kdiff_v += (mesh->kval(YFIXED, kvec_v[kj]) - 
		   mesh->kval(YFIXED, kvec_v[kj-1]));
  av_kdiff_v /= (double)(kvec_v.size()-1);

  const vector<Element2D*>& elem = bspline->supportedElements();
  int nmb_outside_pts = 0;
  int curr_nmb_out;
  double dom;
  bool refined = false;
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

      curr_nmb_out = elem[ki]->getNmbOutsideTol();
      if (curr_nmb_out > 0)
	{
	  nmb_outside_pts += curr_nmb_out;
	  dom = (umax-umin)*(vmax-vmin);
	  // u_info[kj1-1] += dom*elem[ki]->getAccumulatedError();
	  // v_info[kj2-1] += dom*elem[ki]->getAccumulatedError();
	  u_info[kj1-1] += dom*elem[ki]->getAccumulatedOutside();
	  v_info[kj2-1] += dom*elem[ki]->getAccumulatedOutside();
	  if (umax-umin > 0.9*(kvec_u[kj1]-kvec_u[kj1-1]))
	    u_outside[kj1-1] += curr_nmb_out;
	  if (vmax-vmin > 0.9*(kvec_v[kj2]-kvec_v[kj2-1]))
	    v_outside[kj2-1] += curr_nmb_out;
	}
      else
	{
	  if (umax-umin > 0.9*(kvec_u[kj1]-kvec_u[kj1-1]))
	    u_inside[kj1-1]++;
	  if (vmax-vmin > 0.9*(kvec_v[kj2]-kvec_v[kj2-1]))
	    v_inside[kj2-1]++;
	}

      // Element size
      u_elsize[kj2-1] += (umax-umin);
      v_elsize[kj1-1] += (vmax-vmin);
    } 

  if (nmb_outside_pts == 0)
    return;  // Security. Should not happen

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
    
  // Set threshold for which strips to split
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

  std::set<Element2D*> curr_el;
  double threshhold = std::min(av_info, 0.5*max_info);
  double sizefac = 1.5; //3.0;
  if (dir == 1 || dir == 3)
    {
  double minsize_u = std::max(2.0*usize_min_, 1.0e-8);
  for (kj=0; kj<size1; ++kj)
    {
      double u1 = mesh->kval(XFIXED, kvec_u[kj]);
      double u2 = mesh->kval(XFIXED, kvec_u[kj+1]);
      if (((u_info[kj] >= threshhold || u2-u1 > sizefac*v_elsize[kj]) &&
	   (u2 - u1) >= minsize_u && 
	   (u_inside[kj] == 0 || (double)u_outside[kj] > average_out)) ||
	  u2 - u1 > 1.5*(av_kdiff_u + av_kdiff_v))
	{
	  // Check if a candidate knot value exists
	  double knotval = 0.5*(u1+u2);
	  if (init_knots_u_.size() > 0)
	    {
	      int kk1=0, kk2=0;
	      for (kk1=0; kk1<(int)init_knots_u_.size(); ++kk1)
		if (init_knots_u_[kk1] > u1)
		  break;
	      for (kk2=0; kk2<(int)init_knots_u_.size(); ++kk2)
		if (init_knots_u_[kk2] > u2)
		  break;
	      kk2--;
	      if (kk2 >= kk1)
		knotval = init_knots_u_[(kk1+kk2)/2];
	    }

	  LRSplineSurface::Refinement2D curr_ref;
	  curr_ref.setVal(knotval, bspline->vmin(), bspline->vmax(), XFIXED, 1);

	  // Check if the current refinement can be combined with an existing one
	  appendRef(refs_x, curr_ref, tol);
	  for (size_t ki=0; ki<elem.size(); ++ki)
	    {
	      if (elem[ki]->umin() < knotval && elem[ki]->umax() > knotval)
		curr_el.insert(elem[ki]);
	    }
	  refined = true;
	}
    }
    }

  if (dir == 2 || dir == 3)
    {
  double minsize_v = std::max(2.0*vsize_min_, 1.0e-8);
  for (kj=0; kj<size2; ++kj)
    {
      double v1 = mesh->kval(YFIXED, kvec_v[kj]);
      double v2 = mesh->kval(YFIXED, kvec_v[kj+1]);
      if (((v_info[kj] >= threshhold  || v2-v1 > sizefac*u_elsize[kj]) &&
	  (v2 - v1) >= minsize_v && 
	  (v_inside[kj] == 0 || (double)v_outside[kj] > average_out)) ||
	  v2 - v1 > 1.5*(av_kdiff_u + av_kdiff_v))
	{
	  // Check if a candidate knot value exists
	  double knotval = 0.5*(v1+v2);
	  if (init_knots_v_.size() > 0)
	    {
	      int kk1=0, kk2=0;
	      for (kk1=0; kk1<(int)init_knots_v_.size(); ++kk1)
		if (init_knots_v_[kk1] > v1)
		  break;
	      for (kk2=0; kk2<(int)init_knots_v_.size(); ++kk2)
		if (init_knots_v_[kk2] > v2)
		  break;
	      kk2--;
	      if (kk2 >= kk1)
		knotval = init_knots_v_[(kk1+kk2)/2];
	    }

	  LRSplineSurface::Refinement2D curr_ref;
	  curr_ref.setVal(knotval, bspline->umin(), bspline->umax(), YFIXED, 1);

	  // Check if the current refinement can be combined with an existing one
	  appendRef(refs_y, curr_ref, tol);
	  for (size_t ki=0; ki<elem.size(); ++ki)
	    {
	      if (elem[ki]->vmin() < knotval && elem[ki]->vmax() > knotval)
		curr_el.insert(elem[ki]);
	    }
	  refined = true;
	}
    }
    }
  if (refined)
    {
      elem_div.insert(elem_div.end(), curr_el.begin(), curr_el.end());
    }

  if (alloc)
    delete [] alloc;
}

#if 0
//==============================================================================
void LRSurfApprox::checkFeasibleRef(Element2D* elem, int dir, int iter, 
				    vector<LRSplineSurface::Refinement2D>& refs_x,
				    vector<LRSplineSurface::Refinement2D>& refs_y,
				    vector<Element2D*>& affected)
//==============================================================================
{
  double tol = srf_->getKnotTol();
  double eps = 0.01;

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
	  int nmb_outside, nmb_out_sign;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside, 
				       nmb_out_sign);
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
	  int nmb_outside, nmb_out_sign;
	  curr_el[kj]->getAccuracyInfo(av_err, max_err, nmb_outside, 
				       nmb_out_sign);
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
	  double max_err, av_err;
	  int nmb_outside, nmb_out_sign;
	  curr_el_u[kj]->getAccuracyInfo(av_err, max_err, nmb_outside,
					 nmb_out_sign);
	  int nmb_pts = curr_el_u[kj]->nmbDataPoints();
	  aff_u.push_back(curr_el_u[kj]);
	  if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	    nmb_u++;
	}
    }

  vector<Element2D*> aff_v;
  int nmb_v = 0;
  if (ixv >= 0)
    {
      const vector<Element2D*>& curr_el_v = bsplines[ixv]->supportedElements();
      for (kj=0; kj<curr_el_v.size(); ++kj)
	{
	  double max_err, av_err;
	  int nmb_outside, nmb_out_sign;
	  curr_el_v[kj]->getAccuracyInfo(av_err, max_err, nmb_outside,
					 nmb_out_sign);
	  int nmb_pts = curr_el_v[kj]->nmbDataPoints();
	  aff_v.push_back(curr_el_v[kj]);
	  if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	    nmb_v++;
	}
    }

  // Assemble information
  double fac2 = 0.5;
  std::set<Element2D*> affected_combined;
  if (ixu >= 0 && nmb_u >= nmb_v)
    {
      affected_combined.insert(aff_u.begin(), aff_u.end());
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(u_par, bsplines[ixu]->vmin(), bsplines[ixu]->vmax(),
		      XFIXED, xmult);
      //refs.push_back(curr_ref);
      appendRef(refs_x, curr_ref, tol);
    }
			       
  if (ixv >= 0 && nmb_u >= nmb_v)
    {
      affected_combined.insert(aff_v.begin(), aff_v.end());
      LRSplineSurface::Refinement2D curr_ref;
      curr_ref.setVal(v_par, bsplines[ixv]->umin(), bsplines[ixv]->umax(),
		      YFIXED, ymult);
      //refs.push_back(curr_ref);
      appendRef(refs_y, curr_ref, tol);
    }

    affected.insert(affected.end(), affected_combined.begin(), affected_combined.end());
}
#endif

//==============================================================================
void LRSurfApprox::checkFeasibleRef(Element2D* elem, int dir, int iter,
				    vector<LRSplineSurface::Refinement2D>& refs_x,
				    vector<LRSplineSurface::Refinement2D>& refs_y,
				    vector<Element2D*>& affected)
//==============================================================================
{
  double tol = srf_->getKnotTol();
  double eps = 0.01;

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
  double udel = elem->umax() - elem->umin();
  double minsize_u = (usize_min_ > 0.0) ? 2.0*usize_min_ : 1.0e-8;
  size_t ki, kj;
  double udelmax = 0.0;
  bool udir = true;
  std::set<Element2D*> uelems;
  for (ki=0; ki<nmb; ++ki)
    {
      udelmax = std::max(udelmax, bsplines[ki]->umax() - bsplines[ki]->umin());
      
      // Collect elements
      const vector<Element2D*>& curr_el = bsplines[ki]->supportedElements();
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
  std::set<Element2D*> velems;
  for (ki=0; ki<nmb; ++ki)
    {
       vdelmax = std::max(vdelmax, bsplines[ki]->vmax() - bsplines[ki]->vmin());
      
      // Collect elements
      const vector<Element2D*>& curr_el = bsplines[ki]->supportedElements();
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
  
  // Estimate significant of split in all parameter directions
  vector<Element2D*> element_u(uelems.begin(), uelems.end());
  vector<Element2D*> element_v(velems.begin(), velems.end());
  double fac = 0.1;
  double fac3 = 0.75; //0.95;
  int nmb_u = 0;
  for (ki=0; ki<element_u.size(); ++ki)
    {
      double max_err, av_err;
      int nmb_outside, nmb_out_sign;
      element_u[ki]->getAccuracyInfo(av_err, max_err, nmb_outside, nmb_out_sign);
      int nmb_pts = element_u[ki]->nmbDataPoints();
      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	nmb_u++;
    }

  int nmb_v = 0;
  for (ki=0; ki<element_v.size(); ++ki)
    {
      double max_err, av_err;
      int nmb_outside, nmb_out_sign;
      element_v[ki]->getAccuracyInfo(av_err, max_err, nmb_outside, nmb_out_sign);
      int nmb_pts = element_v[ki]->nmbDataPoints();
      if (nmb_outside > fac*nmb_pts || av_err > fac3*avdist_)
	nmb_v++;
    }


  // if (udir == false && vdir == false)
  //   udir = vdir = true;   // No direction without small elements,
  // // use other criteria
  
  // Assemble information
int div_x = 0, div_y = 0;
  double fac2 = 0.5;
  double sizefac = 3.0;
  if ((dir == 1 || dir == 3) && udir &&
       (nmb_u >= nmb_v || udel > sizefac*vdel || udelmax >= vdelmax))
    {
      for (ki=0; ki<nmb; ++ki)
	{
	  LRSplineSurface::Refinement2D curr_ref;
	  curr_ref.setVal(u_par, bsplines[ki]->vmin(), bsplines[ki]->vmax(), XFIXED, xmult);
	  appendRef(refs_x, curr_ref, tol);
	}
    }
			       
  if ((dir == 2 || dir == 3) && vdir &&
      (nmb_v >= nmb_u || vdel > sizefac*udel || vdelmax >= udelmax))
    {
      for (ki=0; ki<nmb; ++ki)
	{
	  LRSplineSurface::Refinement2D curr_ref;
	  curr_ref.setVal(v_par, bsplines[ki]->umin(), bsplines[ki]->umax(), YFIXED, ymult);
	  appendRef(refs_y, curr_ref, tol);
	}
    }
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

    shared_ptr<SplineSurface> surf1;
    try {
      surf1 = 
      createSurf(currpt, kn2, dim-2, order, order,
		 order, order, &knots_u[0], &knots_v[0], smoothweight,
		 maxdist, avdist, outsideeps);
    }
    catch (...)
      {
      }

#ifdef DEBUG
    if (surf1.get())
      {
	shared_ptr<LRSplineSurface> tmp_srf1(new LRSplineSurface(surf1.get(), 1.0e-6));
	if (tmp_srf1->dimension() == 1)
	  tmp_srf1->to3D();
	tmp_srf1->writeStandardHeader(of3);
	tmp_srf1->write(of3);
      }
#endif

    // Make approximative quadratic surface
    order = 3;
    knots_u.insert(knots_u.begin(), ptbound[0]);
    knots_u.push_back(ptbound[1]);
    knots_v.insert(knots_v.begin(), ptbound[2]);
    knots_v.push_back(ptbound[3]);
    shared_ptr<SplineSurface> surf2;
    try {
      surf2 = 
      createSurf(currpt, kn2, dim-2, order, order,
		 order, order, &knots_u[0], &knots_v[0], smoothweight,
		 maxdist, avdist, outsideeps);
    }
    catch (...)
      {
      }

#ifdef DEBUG
    if (surf2.get())
      {
	shared_ptr<LRSplineSurface> tmp_srf2(new LRSplineSurface(surf2.get(), 1.0e-6));
	if (tmp_srf2->dimension() == 1)
	  tmp_srf2->to3D();
	tmp_srf2->writeStandardHeader(of3);
	tmp_srf2->write(of3);
      }
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
	      Point pos0(dim-2), pos1(dim-2), pos2(dim-2);
	      pos0.setValue(0.0);
	      pos1.setValue(0.0);
	      pos2.setValue(0.0);
	      if (surf0.get())
		pos0 = surf0->ParamSurface::point(upar, vpar);
	      if (surf1.get())
		pos1 = surf1->ParamSurface::point(upar, vpar);
	      if (surf2.get())
		pos2 = surf2->ParamSurface::point(upar, vpar);
	      pos = (pos0+pos1+pos2)/3.0;
	    }
	  else if (initMBA_)
	    {
	      Point pos1(dim-2), pos2(dim-2);
	      pos1.setValue(0.0);
	      pos2.setValue(0.0);
	      Point pos0 = srf_->ParamSurface::point(upar, vpar);
	      if (surf1.get())
		pos1 = surf1->ParamSurface::point(upar, vpar);
	      if (surf2.get())
		pos2 = surf2->ParamSurface::point(upar, vpar);
	      pos = (pos0+pos1+pos2)/3.0;
	    }
	  else
	    {
	      Point pos1(dim-2), pos2(dim-2);
	      pos1.setValue(0.0);
	      pos2.setValue(0.0);
	      if (surf1.get())
		pos1 = surf1->ParamSurface::point(upar, vpar);
	      if (surf2.get())
		pos2 = surf2->ParamSurface::point(upar, vpar);
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

  LRSplineSurface::ElementMap::const_iterator it=srf_->elementsBegin();
  if (it == srf_->elementsEnd())
    return; // No elements


  // Global information
  int dim = srf_->dimension();
  int del = it->second->getNmbValPrPoint();
  if (del == 0)
    del = dim + 3;
  int tot_nmb_pts = (int)points_.size()/(del-1);  // Distance not included
  double u1 = srf_->paramMin(XFIXED);
  double u2 = srf_->paramMax(XFIXED);
  double v1 = srf_->paramMin(YFIXED);
  double v2 = srf_->paramMax(YFIXED);
  double pdsize = (u2-u1)*(v2-v1);
  double nmb_fac = 0.2;
   
  for (; it != srf_->elementsEnd(); ++it)
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
		if (del > dim+3)
		  ghost_points.push_back(1.0);  // Related to outlier detection
	      }
	  it->second->addGhostPoints(ghost_points.begin(), ghost_points.end(),
				     false, del);

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

  for (size_t ki=0; ki<elems2.size(); ++ki)
    {
      vector<double>& points = elems2[ki]->getDataPoints();
      int nmb_pts = elems2[ki]->nmbDataPoints();
      int del = elems2[ki]->getNmbValPrPoint();
      if (del == 0)
	del = dim+3;  // Parameter pair, point and distance

      // Compute distances in data points and update parameter pairs
      // if requested
      vector<double> prev_point_dist(nmb_pts, 0.0);
      if (nmb_pts > 0)
	{
	    if (enable_omp)
	    {
	      computeAccuracyElement_omp(points, nmb_pts, del, rd, elems2[ki],
					 prev_point_dist);
	    }
	    else
	    {
	      computeAccuracyElement(points, nmb_pts, del, rd, elems2[ki],
				     prev_point_dist);
	    }

	  // Local error information
	  double max_err = 0.0;
	  double av_err = 0.0;
	  double acc_err = 0.0;
	  double acc_outside = 0.0;
	  int outside = 0;
	  int outside_sign = 0;
	  double tol;  // Local tolerance taking the possibility for a tolerance 
	  // threshold varying with a high distant from zero into account

	  double *curr;
	  int kj;
	  for (kj=0, curr=&points[0]; kj<nmb_pts; ++kj, curr+=del)
	    {
	      // Accumulate approximation error
	      double dist2 = fabs(curr[del-1]);
	      max_err = std::max(max_err, dist2);
	      acc_err += dist2;
	      tol = aepsge_;
	      for (size_t kr=0; kr<tolerances_.size(); ++kr)
		{
		  if (tolerances_[kr].contains(curr[0], curr[1]))
		    {
		      tol = tolerances_[kr].tol;
		      break;
		    }
		}
	      if (has_var_tol_)
		{
		  double height = curr[del-2];
		  tol = (height < 0.0) ? tol + var_fac_neg_*height :
		    tol + var_fac_pos_*height;
		}
	      tol = std::max(tol, mintol_);
	      //if (dist2 > aepsge_)
	      if (dist2 > tol)
		{
		  outside++;
		  av_err += dist2;
		  acc_outside += (dist2-tol);
		}
	    }
	  if (outside > 0)
	    av_err /= (double)outside;

	  // Store updated accuracy information in the element
	  elems2[ki]->setAccuracyInfo(acc_err, av_err, max_err, outside, 
				      acc_outside, outside_sign);
	}
    }
}

//==============================================================================
void LRSurfApprox::updateGhostPoints(vector<Element2D*>& elems)
//==============================================================================
{
  if (elems.size() == 0)
    return;  // Nothing to update

  // Global information
  int dim = srf_->dimension();
  int del = elems[0]->getNmbValPrPoint();
  if (del == 0)
    del = dim+3;  // Parameter pair, point and distance
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
	    if (del > dim+3)
	      ghost_points.push_back(1.0);  // Related to outlier detection
	  }
      elems[ki]->addGhostPoints(ghost_points.begin(), ghost_points.end(),
				false, del);
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
      int del = it->second->getNmbValPrPoint();
      if (del == 0)
	del = dim+3;  // Parameter pair, point and distance

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
		if (del > dim+3)
		  ghost_points.push_back(1.0);  // Related to outlier detection
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
		if (del > dim+3)
		  ghost_points.push_back(1.0);  // Related to outlier detection
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
		if (del > dim+3)
		  ghost_points.push_back(1.0);  // Related to outlier detection
	      }
	}

      if (ghost_points.size() > 0)
	{
	  it->second->addGhostPoints(ghost_points.begin(), ghost_points.end(),
				     false, del);

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

#ifdef DEBUG
  std::ofstream of("constrained_coef.g2");
#endif

  vector<LRBSpline2D*> free_bsplines;
  LRSplineSurface::BSplineMap::const_iterator it1 = srf_->basisFunctionsBegin();
  for (; it1 != srf_->basisFunctionsEnd(); ++it1)
    {
      bool mod = false;
      Point coef = it1->second->Coef();
      double val = coef[0];
      if (has_min_constraint_)
	{
	  if (minval_ > coef[0])
	    mod = true;
	  coef[0] = std::max(coef[0], minval_);
	}
      if (has_max_constraint_)
	{
	  if (maxval_ < coef[0])
	    mod = true;
	  coef[0] = std::min(coef[0], maxval_);
	}
      if (has_local_constraint_)
	{
	  double bb[2], curr_bb[2];  // Bounding box
	  bb[0] = std::numeric_limits<double>::max();
	  bb[1] = std::numeric_limits<double>::lowest();
	  
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
	      double low = bb[0] - constraint_fac_*(bb[1] - bb[0]);
	      double high = bb[1] + constraint_fac_*(bb[1] - bb[0]);
	      if (coef[0] < low || coef[0] > high)
		mod = true;
	      coef[0] = std::max(coef[0], low);	
	      coef[0] = std::min(coef[0], high);
	    }
	  else
	    free_bsplines.push_back(it1->second.get());
#ifdef DEBUG
	  if (mod)
	    {
	      Point greville = it1->second->getGrevilleParameter();
	      of << "400 1 0 4 0 100 155 255" << std::endl;
	      of << "1" << std::endl;
	      of << greville << " " << val << std::endl;
	      of << "400 1 0 4 155 100 0 255" << std::endl;
	      of << "1" << std::endl;
	      of << greville << " " << coef[0] << std::endl;
	    }
#endif
	}
      if (mod)
	srf_->setCoef(coef, it1->second.get());
    }
  for (size_t ki=0; ki<free_bsplines.size(); ++ki)
    {
      // The B-spline does not overlap any data points. Modify coefficient
      // to improve the smoothness of the non-trimmed surface.
      // First fetch coefficient information from nearby B-splines
      double av_val = 0.0;
      int nmb_val = 0;
      const vector<Element2D*> elements = 
	free_bsplines[ki]->supportedElements();
      for (size_t kj=0; kj<elements.size(); ++kj)
	{
	  const vector<LRBSpline2D*> supp_bsplines = 
	    elements[kj]->getSupport();
	  for (size_t kr=0; kr<supp_bsplines.size(); ++kr)
	    {
	      // Check if the B-spline overlaps any data points
	      const vector<Element2D*> elements2 = 
		supp_bsplines[kr]->supportedElements();
	      size_t kh;
	      for (kh=0; kh<elements2.size(); ++kh)
		if (elements2[kh]->hasDataPoints() ||
		    elements2[kh]->hasSignificantPoints())
		  break;
	      if (kh < elements2.size())
		{
		  Point next_coef = supp_bsplines[kr]->Coef();
		  av_val += next_coef[0];
		  nmb_val++;
		}
	    }
	}
      if (nmb_val > 0)
	{
	  av_val /= (double)nmb_val;
	  Point new_coef(1);
	  new_coef[0] = av_val;
	  Point tmp = free_bsplines[ki]->Coef();
	  srf_->setCoef(new_coef, free_bsplines[ki]);
#ifdef DEBUG
	  Point greville = free_bsplines[ki]->getGrevilleParameter();
	  of << "400 1 0 4 0 100 155 255" << std::endl;
	  of << "1" << std::endl;
	  of << greville << " " << tmp << std::endl;
	  of << "400 1 0 4 155 100 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << greville << " " << av_val << std::endl;
#endif
	}
    }
}

//==============================================================================
void LRSurfApprox::getCandElements(double x, double y, double rad, 
				   Element2D* start_elem,
				   vector<Element2D*>& elems)
//==============================================================================
{
  double umin = start_elem->umin();
  double umax = start_elem->umax();
  double vmin = start_elem->vmin();
  double vmax = start_elem->vmax();
  if (umin <= x-rad && umax >= x+rad && 
      vmin <= y-rad && vmax >= y+rad)
    {
      // The neighbourhood domain is complete inside the start element
      return;
    }

  // Fetch the B-splines covering the start element
  vector<LRBSpline2D*> bsupp = start_elem->getSupport();
  for (size_t ka=0; ka<bsupp.size(); ++ka)
    {
      vector<Element2D*> esupp = bsupp[ka]->supportedElements();
      for (size_t kb=0; kb<esupp.size(); ++kb)
	{
	  // Check if the element overlaps the neighbourhood domain
	  double umin2 = esupp[kb]->umin();
	  double umax2 = esupp[kb]->umax();
	  double vmin2 = esupp[kb]->vmin();
	  double vmax2 = esupp[kb]->vmax();
	  if (umin2 >= x+rad || umax2 <= x-rad || 
	      vmin2 >= y+rad || vmax2 <= y-rad)
	    continue;

	  // Check if the element is found already
	  size_t kc;
	  for (kc=0; kc<elems.size(); ++kc)
	    if (elems[kc] == esupp[kb])		
	      break;
	  if (kc < elems.size())
	    continue;

	  // A new element inside the neighbourhood domain.
	  // Collect and continue searching
	  elems.push_back(esupp[kb]);
	  getCandElements(x, y, rad, esupp[kb], elems);
	}
    }
}

//==============================================================================
void LRSurfApprox::appendRef(vector<LRSplineSurface::Refinement2D>& refs,
			     LRSplineSurface::Refinement2D& curr_ref, double tol)
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
