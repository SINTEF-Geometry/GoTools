
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

#include "GoTools/creators/ApproxSurf.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include <algorithm>
#include <fstream>
#include <iterator>

//#define DEBUG

using namespace Go;
using std::vector;
using std::max;
using std::min;
using std::cout;
using std::endl;
using std::back_inserter;

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
#include "GoTools/geometry/Utils.h"     // make std::min and std::max work (redefined in boost/smart_ptr.hpp)
#endif


ApproxSurf::ApproxSurf()
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxSurf.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  prevdist_ = maxdist_ = -10000.0;
  prevav_ = avdist_ = 0;
  dim_ = 3;
  aepsge_ = 0.01;
  smoothweight_ = 1.0e-3; //1e-9;
  smoothfac_ = 1.0;
  outsideeps_ = 0;
  constdir_ = 0;
  use_normals_ = false;
  close_belt_ = false;
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 1;
  pts_stabil_ = 0;
  norm_stabil_ = 0;
  orig_ = false;
  repar_ = true;
  refine_ = true;
  c1fac1_ = 0.0;
  c1fac2_ = 0.0;
}

//***************************************************************************

ApproxSurf::ApproxSurf(std::vector<shared_ptr<SplineCurve> >& crvs,
			   const std::vector<double>& points, 
			   const std::vector<double>& parvals,
			   double domain[],
			   int dim, double aepsge,
		       int constdir, bool repar)
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxSurf.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  prevdist_ = maxdist_ = -10000.0;
  prevav_ = avdist_ = 0;
  outsideeps_ = 0;
  dim_ = dim;
  aepsge_ = aepsge;
  smoothweight_ = 1.0e-3; //1.0e-9; 
  constdir_ = constdir;
  use_normals_ = false;
  close_belt_ = false;
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 1;

  pts_stabil_ = 0;
  norm_stabil_ = 0;
  points_ = points;
  parvals_ = parvals;
  orig_ = false;
  repar_ = repar;
  refine_ = true;
  c1fac1_ = 0.0;
  c1fac2_ = 0.0;


  makeInitSurf(crvs, domain);
  init_srf_ = shared_ptr<SplineSurface>(curr_srf_->clone());
  smoothfac_ = 1.0/((curr_srf_->endparam_u() - curr_srf_->startparam_u()) +
		    (curr_srf_->endparam_v() - curr_srf_->startparam_v()));

}

//***************************************************************************

ApproxSurf::ApproxSurf(shared_ptr<SplineSurface>& srf,
		       const std::vector<double>& points, 
		       const std::vector<double>& parvals,
		       int dim, double aepsge, int constdir,
		       bool approx_orig,
		       bool close_belt, int nmb_stabil,
		       bool repar)
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxSurf.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  prevdist_ = maxdist_ = -10000.0;
  prevav_ = avdist_ = 0;
  outsideeps_ = 0;
  dim_ = dim;
  aepsge_ = aepsge;
  smoothweight_ = 1.0e-3; // 1.0e-9;
  constdir_ = constdir;
  use_normals_ = false;
  close_belt_ = close_belt;
  edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = 1;
  pts_stabil_ = nmb_stabil;
  norm_stabil_ = 0;
  orig_ = approx_orig;
  repar_ = repar;
  refine_ = true;
  c1fac1_ = 0.0;
  c1fac2_ = 0.0;

  points_ = points;
  parvals_ = parvals;

  smoothfac_ = 1.0/((srf->endparam_u() - srf->startparam_u()) +
		    (srf->endparam_v() - srf->startparam_v()));

  curr_srf_ = srf;
  init_srf_ = shared_ptr<SplineSurface>(srf->clone());
}

//***************************************************************************


ApproxSurf::~ApproxSurf()
   //--------------------------------------------------------------------------
   //     Destructor for class ApproxSurf.
   //
   //     Purpose : 
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
//   
}

//***************************************************************************

int ApproxSurf::makeInitSurf(std::vector<shared_ptr<SplineCurve> >& crvs, 
			       double domain[])
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial surface as a Coons patch
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  // Make surface
  vector<shared_ptr<ParamCurve> > bdcrvs(crvs.size());
  std::copy(crvs.begin(), crvs.end(), bdcrvs.begin());
  CurveLoop boundary(bdcrvs, 10e-6);
  curr_srf_ =
      shared_ptr<SplineSurface>(CoonsPatchGen::createCoonsPatch(boundary));

  // Scale surface to live on the given domain
  curr_srf_->setParameterDomain(domain[0], domain[1],
				domain[2], domain[3]);

  return 0;
}


//***************************************************************************

void
ApproxSurf::spline_space_cont(shared_ptr<SplineSurface> sf, int& nmbc1, 
			      int& nmbc2)
//--------------------------------------------------------------------------
//  Purpose : Count the continuity of the spline space corresponding to
//            a given surface.
//
//  Input   : sf    - SurfacCurve represented in one spline space in each
//                    parameter direction.
//            
//
//  Output  : nmbc1 - Continuity of the spline space in 1. par. dir.
//            nmbc2 - Continuity of the spline space in 2. par. dir.
//
//  Calls   : 
//
//  Written by : Vibeke Skytt,  SINTEF. 12.99.
//--------------------------------------------------------------------------
{
  int kcurr, kmax, ki, kj;

  int in1 = sf->numCoefs_u();
  vector<double>::const_iterator st1 = sf->basis_u().begin();
  for (ki=sf->order_u(), kmax=0; ki<in1; ki=kj)
    {
      for (kj=ki+1, kcurr=1; kj<in1 && st1[ki]==st1[kj]; kj++, kcurr++);
      kmax = std::max(kmax, kcurr);
    }
  nmbc1 = sf->order_u() - kmax - 1;

  int in2 = sf->numCoefs_v();
  vector<double>::const_iterator st2 = sf->basis_v().begin();
  for (ki=sf->order_v(), kmax=0; ki<in2; ki=kj)
    {
      for (kj=ki+1, kcurr=1; kj<in2 && st2[ki]==st2[kj]; kj++, kcurr++);
      kmax = std::max(kmax, kcurr);
    }
  nmbc2 = sf->order_v() - kmax - 1;
}


//===========================================================================

int ApproxSurf::get_min_deriv(shared_ptr<SplineSurface> sf, double support_mult)
//--------------------------------------------------------------------------
//  Purpose : Find the minimum derivative in sf (counting mult of inner
//            knots & support of neighbour intervals).
//
//  Input   : cv           : The given curve.
//            
//
//  Output  : return value : min_deriv
//
//  Written by : Vibeke Skytt,  SINTEF. 04.99.
//  Modified by : Jon Mikkelsen, SINTEF. 03.01. The parametrization is
//               changed so that the knot intervals are closer to 1.
//  Modified by : Jon Mikkelsen, SINTEF. 08.01. Removed the previous change
//  Modified by : Sverre Briseid, SINTEF. 08.04. Improved analysis of spline
//                spaces.
//--------------------------------------------------------------------------
{
  // find the maximum number of derivatives to be calculated
  int min_deriv=0, nmbc1, nmbc2;
  spline_space_cont(sf, nmbc1, nmbc2);
  min_deriv = std::min(std::min(nmbc1+1, nmbc2+1), 3);

  double dist1, dist2;
  int in1=sf->numCoefs_u();
  int in2=sf->numCoefs_v();
  
  double min_support_frac = HUGE;
  double tol = -1.0;
  int ki;
  int min_deriv_u = min_deriv;
  while (true)
    {
      min_support_frac = HUGE; // We reset the value for each iteration.
      int k1=sf->order_u() - min_deriv_u;
      vector<double>::const_iterator st1 = sf->basis_u().begin();
      for(ki=min_deriv; ki<in1-1; ki++)
	{
	  dist1=(st1[ki+k1]-st1[ki]);
	  dist2=(st1[ki+k1+1]-st1[ki+1]);
	  double local_frac = dist2/dist1;
	  if (local_frac > 1.0)
	    local_frac = 1/local_frac;
	  if (local_frac < min_support_frac)
	    {
	      min_support_frac = local_frac;
	    }
	}
      tol=pow(1e-08,1.0/min_deriv_u); //(1e-10,1.0/min_deriv);
				    ////(1e-8,1.0/min_deriv);
      tol *= support_mult;
      if (min_support_frac<tol && min_deriv_u>1)
	{
	  min_deriv_u--;
// 	  tol=pow(1e-08,1.0/min_deriv);//pow(1e-10,1.0/min_deriv);
// 	  tol *= support_mult;
	}
      else
	{
	  break;
	}
    }

  int min_deriv_v = min_deriv;
  while (true)
    {
      min_support_frac = HUGE; // We reset the value for each iteration.
      int k1=sf->order_v() - min_deriv_v;
      vector<double>::const_iterator st2 = sf->basis_v().begin();
      for(ki=min_deriv; ki<in2-1; ki++)
	{
	  dist1=(st2[ki+k1]-st2[ki]);
	  dist2=(st2[ki+k1+1]-st2[ki+1]);
	  double local_frac = dist2/dist1;
	  if (local_frac > 1.0)
	    local_frac = 1/local_frac;
	  if (local_frac < min_support_frac)
	    {
	      min_support_frac = local_frac;
	    }
	}
      tol=pow(1e-08,1.0/min_deriv_v); //(1e-10,1.0/min_deriv);
				    ////(1e-8,1.0/min_deriv);
      tol *= support_mult;
      if (min_support_frac<tol && min_deriv_v>1)
	{
	  min_deriv_v--;
// 	  tol=pow(1e-08,1.0/min_deriv);//pow(1e-10,1.0/min_deriv);
// 	  tol *= support_mult;
	}
      else
	{
	  break;
	}
    }

  if(min_support_frac<tol)
    return 0;

  min_deriv = std::min(min_deriv_u, min_deriv_v);
  return min_deriv;
}

int ApproxSurf::makeSmoothSurf()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Adapt the initial surface to approximate the given
   //               data points.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
    SmoothSurf srfgen;
    int stat = 0;

    // int ki, kj;
    int seem[2];
    seem[0] = seem[1] = 0;
    double support = 
      std::max(curr_srf_->endparam_u() - curr_srf_->startparam_u(),
	       curr_srf_->endparam_v() - curr_srf_->startparam_v());
    int min_der = get_min_deriv(curr_srf_, support);
    double wgt1 = 0.0;
    double wgt3 = (min_der >= 3) ? 0.5*smoothweight_ : 0.0;
    double wgt2 = (1.0 - wgt3)*smoothweight_;
    wgt3 *= smoothweight_;
//     double wgt2 = smoothweight_/(smoothfac_*smoothfac_);
//     double wgt3 = smoothweight_/(smoothfac_*smoothfac_*smoothfac_*smoothfac_);
    double normweight = use_normals_ ? 1e-5 : 0.0;
    double weight_sum = wgt1 + wgt2 + wgt3 + normweight;
    if (weight_sum > 1.0) {
	wgt1 /= weight_sum;
	wgt2 /= weight_sum;
	wgt3 /= weight_sum;
	normweight /= weight_sum;
    }
    double approxweight = 1.0 - wgt1 - wgt2 - wgt3 - normweight;
    std::vector<double> pt_weight(parvals_.size()/2, 1.0);
    double wgt_orig = 0.0;
    if (orig_)
      wgt_orig = 0.1*approxweight;
    approxweight -= wgt_orig;

    srfgen.attach(curr_srf_, seem, &coef_known_[0], 0, use_normals_);

    if (smoothweight_ > 0.0)
      srfgen.setOptimize(wgt1, wgt2, wgt3);

    srfgen.setLeastSquares(points_, parvals_,
				  pt_weight, approxweight);

    if (use_normals_) {
	stat = srfgen.setNormalCond(norm_points_, norm_parvals_,
				    pt_weight, normweight);
	if (stat < 0)
	    return stat;
    }

    srfgen.approxOrig(wgt_orig);

    if (c1fac1_ > 0.0)
      srfgen.setPeriodicity(1, 1, c1fac1_, 0.0);
    if (c1fac2_ > 0.0)
      srfgen.setPeriodicity(2, 1, c1fac2_, 0.0);
      
    try {
    stat = srfgen.equationSolve(curr_srf_);
    }
    catch (...)
      {
	return -1;
      }
    if (stat < 0)
	return stat;

    return 0;
}

//***************************************************************************

int ApproxSurf::checkAccuracy(vector<double>& acc_outside_u,
			      vector<int>& nmb_outside_u,
			      vector<double>& acc_outside_v,
			      vector<int>& nmb_outside_v)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Check the accuracy of the current surface compared to
   //               the given data points.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
//   double par_tol = 0.000000000001;
  maxdist_ = -10000.0;
  avdist_ = 0.0;
  outsideeps_ = 0;

  acc_outside_u.resize(curr_srf_->order_u() + curr_srf_->numCoefs_u()-1);
  nmb_outside_u.resize(curr_srf_->order_u() + curr_srf_->numCoefs_u()-1);
  acc_outside_v.resize(curr_srf_->order_v() + curr_srf_->numCoefs_v()-1);
  nmb_outside_v.resize(curr_srf_->order_v() + curr_srf_->numCoefs_v()-1);
  fill(acc_outside_u.begin(), acc_outside_u.end(), 0.0);
  fill(nmb_outside_u.begin(), nmb_outside_u.end(), 0);
  fill(acc_outside_v.begin(), acc_outside_v.end(), 0.0);
  fill(nmb_outside_v.begin(), nmb_outside_v.end(), 0);

  // Traverse the data points and check the quality of the surface
  // approximation
  int nmbpnt = (int)parvals_.size()/2;
  double dist;
  Point pos(3);

  for (int ki=0; ki<nmbpnt; ki++)
    {
      // Evaluate the surface in the current parameter_value
      curr_srf_->point(pos, parvals_[2*ki], parvals_[2*ki+1]);

      // Compute the distance to the corresponding data point
      vector<double>::const_iterator it = points_.begin() + ki * dim_;
      dist = pos.dist(Point(it, it + dim_));

      if (dist > aepsge_)
	{
	  int left1 = curr_srf_->basis_u().knotInterval(parvals_[2*ki]);
	  int left2 = curr_srf_->basis_v().knotInterval(parvals_[2*ki+1]);
	  acc_outside_u[left1] += dist;
	  nmb_outside_u[left1]++;
	  acc_outside_v[left2] += dist;
	  nmb_outside_v[left2]++;
	}

      // Statistics
      maxdist_ = std::max(maxdist_, dist);
      avdist_ += dist;
      if (dist > aepsge_)
	outsideeps_++;
    }

  avdist_ /= (double)nmbpnt;
	
  return 0;
}

//***************************************************************************

int ApproxSurf::reParam()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Reparameterize the data points
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  double cpar, cparprev=-10000.0;
  shared_ptr<SplineCurve> qc;

  // Traverse all data points
  int nbpt = (int)parvals_.size()/2;
  double *pt, *par;
  int ki, kc;
  double guess[2], clpar[2];
  double cldist;
  Point clpoint(dim_);

  for (ki=0, pt=&points_[0], par=&parvals_[0]; ki<nbpt; ki++, pt+=dim_, par+=2) {
      guess[0] = par[0];
      guess[1] = par[1];
      if (constdir_ == 0) {
	  curr_srf_->closestPoint(Point(pt,pt+dim_), clpar[0],
				  clpar[1], clpoint, cldist,
				  aepsge_, NULL, &guess[0]);
	  par[0] = clpar[0];
	  par[1] = clpar[1];
      } else if (constdir_ == 1 || constdir_ == 2) {
	  kc = 2 - constdir_;
	  cpar = par[constdir_ - 1];
	  if (qc.get() == 0 || cpar != cparprev) {
	      qc.reset(curr_srf_->constParamCurve(cpar, (constdir_ != 1)));
	  }
	  
	  qc->closestPoint(Point(pt,pt+dim_), qc->startparam(),
			   qc->endparam(), clpar[kc], clpoint,
			   cldist, guess+kc);
	  par[kc] = clpar[kc];
	  
	  cparprev = cpar;
      }
  }
  return 0;
}

//***************************************************************************

int ApproxSurf::doApprox(int max_iter, int keep_init)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Generate the approximating surface
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  int stat = 0;
#ifdef DEBUG
  std::ofstream out("sf_approx.g2");
  curr_srf_->writeStandardHeader(out);
  curr_srf_->write(out);
#endif

  vector<double> acc_outside_u;
  vector<int> nmb_outside_u;
  vector<double> acc_outside_v;
  vector<int> nmb_outside_v;

  stat = checkAccuracy(acc_outside_u, nmb_outside_u,
		       acc_outside_v, nmb_outside_v);
  if (stat < 0)
    return stat;
#ifdef DEBUG
  cout << "iter 0,  max " << maxdist_ << " average " << avdist_;
  cout << "# out " << outsideeps_ << endl;
#endif

  for (int ki=0; ki<max_iter; ki++)
    {
      // Reparameterize the data points
      if (repar_)
	{
	  stat = reParam();
	  if (stat < 0)
	    return stat; 
	}

      // Set free and fixed coefficients
      if (close_belt_)
	coefKnownFromPoints();
      else
	setCoefKnown();

      // Approximate
      prev_srf_ = shared_ptr<SplineSurface>(curr_srf_->clone());  // Store last version
      prevdist_ = maxdist_;
      prevav_ = avdist_;
      stat = makeSmoothSurf();
      if (stat < 0)
	return stat;

#ifdef DEBUG
      curr_srf_->writeStandardHeader(out);
      curr_srf_->write(out);
#endif

      stat = checkAccuracy(acc_outside_u, nmb_outside_u,
			   acc_outside_v, nmb_outside_v);
      if (stat < 0)
	return stat;


#ifdef DEBUG
      cout << "iter " << ki+1 << ", max " << maxdist_ << " average ";
      cout << avdist_ << " # out " << outsideeps_ << endl;
#endif

      if (maxdist_ < aepsge_)
	break;

      if (maxdist_ > prevdist_)
	{
	  curr_srf_ = prev_srf_;
	  maxdist_ = prevdist_;
	  avdist_ = prevav_;
	  break;
	}

      if (ki < keep_init)
	{
	  // Return to initial surface, but use the spline space of
	  // the current one
	  BsplineBasis basis_u = curr_srf_->basis_u();
	  BsplineBasis basis_v = curr_srf_->basis_v();
	  vector<double> diff1;
	  std::set_difference(basis_u.begin(), basis_u.end(),
			      init_srf_->basis_u().begin(), init_srf_->basis_u().end(),
			      std::back_inserter(diff1));
	  vector<double> diff2;
	  std::set_difference(basis_v.begin(), basis_v.end(),
			      init_srf_->basis_v().begin(), init_srf_->basis_v().end(),
			      std::back_inserter(diff2));
 	  curr_srf_ = shared_ptr<SplineSurface>(init_srf_->clone());
	  curr_srf_->insertKnot_u(diff1);
	  curr_srf_->insertKnot_v(diff2);
	}

      if (refine_ && ki < max_iter-1)
	stat = refineSplineSpace(acc_outside_u, nmb_outside_u,
				 acc_outside_v, nmb_outside_v);
      if (stat < 0)
	return stat;

    }

  return 0;
}

//***************************************************************************

shared_ptr<SplineSurface> 
ApproxSurf::getApproxSurf(double& maxdist, double& avdist,
			  int& nmb_out_eps, int max_iter,
			  int keep_init)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Return the approximating surface and parameters
   //               describing the accuracy.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
    // Generate the approximating surface
    int stat = 0;
    if (max_iter == 0) {
	stat = makeSmoothSurf();
	if (stat != 0) {
	    MESSAGE("Error in making smooth surface: " << stat);
	}
	vector<double> acc_outside_u;
	vector<int> nmb_outside_u;
	vector<double> acc_outside_v;
	vector<int> nmb_outside_v;
	stat = checkAccuracy(acc_outside_u, nmb_outside_u,
			     acc_outside_v, nmb_outside_v);
	if (stat != 0) {
	    MESSAGE("Error in checking accuracy: " << stat);
	}
    } else {
      stat = doApprox(max_iter, keep_init);
	if (stat != 0) {
	    MESSAGE("Error in approximating: " << stat);
	}
	
    }

  // Set accuracy parameters
  maxdist = maxdist_;
  avdist = avdist_;
  nmb_out_eps = outsideeps_;

  // Return surface

  return curr_srf_;
}

void
ApproxSurf::coefKnownFromPoints()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Set coef_known in order to free only the part
   //               of the surface where the approximation points are
   //               situated
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
  if (parvals_.size() == 0)
    return;  // No information

  
  // Initially all coeficients are kept fixed
  int kn1 = curr_srf_->numCoefs_u();
  int kn2 = curr_srf_->numCoefs_v();
  coef_known_.resize(kn1*kn2);
  fill(coef_known_.begin(), coef_known_.end(), 1);

  // Free the coefficients corresponding to the basis functions
  // influencing the area of the approximation conditions

  int ki;
  int c1_1, c1_2, c2_1, c2_2;
  int k1, k2;
  int nmb = (int)parvals_.size() - 2*pts_stabil_;
  for (ki=0; ki<nmb; ki+=2)
    {
      // Find knot interval in both parameter directions
      curr_srf_->basis_u().coefsAffectingParam(parvals_[ki], c1_1, c1_2);
      curr_srf_->basis_v().coefsAffectingParam(parvals_[ki+1], c2_1, c2_2);

      // Set corresponding coefficients free
      for (k1=c1_1; k1<=c1_2; ++k1)
	for (k2=c2_1; k2<=c2_2; ++k2)
	  coef_known_[k2*kn1+k1] = 0;
    }
  
  int nmb2 = (int)norm_parvals_.size() - 2*norm_stabil_;
  for (ki=0; ki<nmb2; ki+=2)
    {
      // Find knot interval in both parameter directions
      curr_srf_->basis_u().coefsAffectingParam(norm_parvals_[ki], c1_1, c1_2);
      curr_srf_->basis_v().coefsAffectingParam(norm_parvals_[ki+1], c2_1, c2_2);

      // Set corresponding coefficients free
      for (k1=c1_1; k1<=c1_2; ++k1)
	for (k2=c2_1; k2<=c2_2; ++k2)
	  coef_known_[k2*kn1+k1] = 0;
    }

  // Fix specified edges
  int kr, kj;
  for (kr = 0; kr < edge_derivs_[0]; ++kr)
    for (kj = 0; kj < kn1; ++kj)
      coef_known_[kj+kn1*kr] = 1;

  for (kr = 0; kr < edge_derivs_[1]; ++kr)
    for (kj = 0; kj < kn2; ++kj)
      coef_known_[(kj+1)*kn1-1-kr] = 1;

  for (kr = 0; kr < edge_derivs_[2]; ++kr)
    for (kj = 0; kj < kn1; ++kj)
      coef_known_[(kn2-1-kr)*kn1+kj] = 1;

  for (kr = 0; kr < edge_derivs_[3]; ++kr)
    for (kj = 0; kj < kn2; ++kj)
      coef_known_[kj*kn1+kr] = 1;
    
 }
 

void ApproxSurf::setCoefKnown()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Set coef_known according to fixed edges
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
  int kn1 = curr_srf_->numCoefs_u();
  int kn2 = curr_srf_->numCoefs_v();
  coef_known_.resize(kn1*kn2);
  fill(coef_known_.begin(), coef_known_.end(), 0);

  int ki, kj;
  for (ki = 0; ki < edge_derivs_[0]; ++ki)
    for (kj = 0; kj < kn1; ++kj)
      coef_known_[kj+kn1*ki] = 1;

  for (ki = 0; ki < edge_derivs_[1]; ++ki)
    for (kj = 0; kj < kn2; ++kj)
      coef_known_[(kj+1)*kn1-1-ki] = 1;

  for (ki = 0; ki < edge_derivs_[2]; ++ki)
    for (kj = 0; kj < kn1; ++kj)
      coef_known_[(kn2-1-ki)*kn1+kj] = 1;

  for (ki = 0; ki < edge_derivs_[3]; ++ki)
    for (kj = 0; kj < kn2; ++kj)
      coef_known_[kj*kn1+ki] = 1;


}

int ApproxSurf::refineSplineSpace(vector<double>& acc_outside_u,
				  vector<int>& nmb_outside_u,
				  vector<double>& acc_outside_v,
				  vector<int>& nmb_outside_v)
{
  // Count number of knot intervals with too high approximation error in
  // both parameter directions
  int nmb_err_u = 0, nmb_err_v = 0;
  double max_acc_u = 0.0, max_acc_v = 0.0;
  double av_acc_u = 0.0, av_acc_v = 0.0;
  size_t ki;
  for (ki=0; ki<nmb_outside_u.size(); ++ki)
    if (nmb_outside_u[ki] > 0)
      {
	nmb_err_u++;
	max_acc_u = std::max(max_acc_u, acc_outside_u[ki]);
	av_acc_u += acc_outside_u[ki];
      }
  if (nmb_err_u > 0)
    av_acc_u /= (double)(nmb_err_u);

  for (ki=0; ki<nmb_outside_v.size(); ++ki)
    if (nmb_outside_v[ki] > 0)
      {
	nmb_err_v++;
	max_acc_v = std::max(max_acc_v, acc_outside_v[ki]);
	av_acc_v += acc_outside_v[ki];
      }
  if (nmb_err_v > 0)
    av_acc_v /= (double)(nmb_err_u);

  if (nmb_err_u == 0 && nmb_err_v == 0)
    return 0;

  // Decide on criterion for knot insertion
  double limit_acc_err = 0.0;
  if (nmb_err_u + nmb_err_v > 12)
    limit_acc_err = std::max(0.5*(max_acc_u+av_acc_u), 0.5*(max_acc_v+av_acc_v));
  else if (nmb_err_u + nmb_err_v > 8)
    limit_acc_err = std::max(av_acc_u, av_acc_v);
  else if (nmb_err_u + nmb_err_v > 4)
    limit_acc_err = std::min(av_acc_u, av_acc_v);
  
  int kn1 = curr_srf_->numCoefs_u();
  int kn2 = curr_srf_->numCoefs_v();
  double min_int_u, min_int_v;
  min_int_u = (curr_srf_->endparam_u() - curr_srf_->startparam_u())/(double)(100*kn1);
  min_int_v = (curr_srf_->endparam_v() - curr_srf_->startparam_v())/(double)(100*kn2);

  // Set new knots
  vector<double> knot_u, knot_v;
  for (ki=0; ki<nmb_outside_u.size(); ++ki)
    {
      if (acc_outside_u[ki] > limit_acc_err)
	{
	  double t1 = curr_srf_->basis_u().begin()[ki];
	  double t2 = curr_srf_->basis_u().begin()[ki+1];
	  if (t2 - t1 > min_int_u)
	    knot_u.push_back(0.5*(t1+t2));
	}
    }

  for (ki=0; ki<nmb_outside_v.size(); ++ki)
    {
      if (acc_outside_v[ki] > limit_acc_err)
	{
	  double t1 = curr_srf_->basis_v().begin()[ki];
	  double t2 = curr_srf_->basis_v().begin()[ki+1];
	  if (t2 - t1 > min_int_v)
	    knot_v.push_back(0.5*(t1+t2));
	}
    }

  // Insert knots
  if (knot_u.size() > 0)
    curr_srf_->insertKnot_u(knot_u);

  if (knot_v.size() > 0)
    curr_srf_->insertKnot_v(knot_v);

  return 1;
}

void  ApproxSurf::setC1Approx(double fac1, double fac2)
{
  c1fac1_ = fac1;
  c1fac2_ = fac2;
}
