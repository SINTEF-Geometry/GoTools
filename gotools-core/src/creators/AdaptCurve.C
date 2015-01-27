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

//#define DEBUG
#include "GoTools/creators/AdaptCurve.h"
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/Utils.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineUtils.h" // debugging

using namespace Go;
using std::vector;
using std::make_pair;

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
#include "GoTools/geometry/Utils.h"     // make std::min and std::max work (redefined in boost/smart_ptr.hpp)
#endif

//***************************************************************************

AdaptCurve::AdaptCurve(const EvalCurve *evalcrv, double aepsge)
   //--------------------------------------------------------------------------
   //     Constructor for class AdaptCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
  : prev_maxdist_(10000.0), prev_avdist_(0.0), maxdist_(10000.0),
    maxdist_inpoints_(10000.0), avdist_(0.0), aepsge_(aepsge),
    smoothweight_(0.000000001), smoothfac_(1.0), evalcrv_(evalcrv),
    dim_(evalcrv->dim()), cont_(2), 
    order_(4), init_sample_(0)
{
  fix_[0] = fix_[1] = 0;
  min_sample_par_ = evalcrv->start();
  max_sample_par_ = evalcrv->end();
  makeInitCurve();
}

//***************************************************************************

AdaptCurve::AdaptCurve(const EvalCurve *evalcrv, double aepsge, int in, int ik)
   //--------------------------------------------------------------------------
   //     Constructor for class AdaptCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
  : prev_maxdist_(10000.0), prev_avdist_(0.0),
    maxdist_(10000.0), maxdist_inpoints_(10000.0), avdist_(0.0), aepsge_(aepsge),
    smoothweight_(0.000000001), evalcrv_(evalcrv),
    dim_(evalcrv->dim()), cont_(2), 
    order_(ik), init_sample_(0)
{
  fix_[0] = fix_[1] = 0;
  min_sample_par_ = evalcrv->start();
  max_sample_par_ = evalcrv->end();
  makeInitCurve(in, ik);
}

//***************************************************************************

AdaptCurve::AdaptCurve(const EvalCurve *evalcrv, double aepsge, int in, 
		       int ik, std::vector<double>& knots)
   //--------------------------------------------------------------------------
   //     Constructor for class AdaptCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
  : prev_maxdist_(10000.0), prev_avdist_(0.0),
    maxdist_(10000.0), maxdist_inpoints_(10000.0), avdist_(0.0), aepsge_(aepsge),
    smoothweight_(0.000000001), smoothfac_(1.0), evalcrv_(evalcrv),
    dim_(evalcrv->dim()), cont_(2), 
    order_(ik), init_sample_(0)
{
  fix_[0] = fix_[1] = 0;
  min_sample_par_ = evalcrv->start();
  max_sample_par_ = evalcrv->end();
  makeInitCurve(in, ik, knots);
}

//***************************************************************************

AdaptCurve::AdaptCurve(const EvalCurve *evalcrv, double aepsge, 
		       shared_ptr<SplineCurve> curve)
   //--------------------------------------------------------------------------
   //     Constructor for class AdaptCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
  : prev_maxdist_(10000.0), prev_avdist_(0.0),
    maxdist_(10000.0), maxdist_inpoints_(10000.0), avdist_(0.0), aepsge_(aepsge),
    smoothweight_(0.000000001), smoothfac_(1.0), evalcrv_(evalcrv),
    dim_(evalcrv->dim()), cont_(2), 
    order_(curve->order()), init_sample_(0)
{
  fix_[0] = fix_[1] = 1;
  min_sample_par_ = evalcrv->start();
  max_sample_par_ = evalcrv->end();
  curr_crv_ = curve;
  prev_crv_ = curve;
}

//***************************************************************************

AdaptCurve::~AdaptCurve()
   //--------------------------------------------------------------------------
   //     Destructor for class AdaptCurve.
   //
   //     Purpose : 
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
//   if (curr_crv_)
//     freeCurve(curr_crv_);
}

//***************************************************************************

void AdaptCurve::makeInitCurve()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
  // Initially a Bezier curve
  makeInitCurve(order_, order_);
}


//***************************************************************************

void AdaptCurve::makeInitCurve(int in, int ik, std::vector<double> knots)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
  // Define dummy coefficient vector
  std::vector<double> coefs(in*dim_, 0.0);

  // Make initial curve
  curr_crv_ = shared_ptr<SplineCurve>
    ( new SplineCurve(in, ik, &knots[0], &coefs[0], dim_));
 
}


//***************************************************************************

void AdaptCurve::makeInitCurve(int in, int ik)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
  std::vector<double> knots;
  knots.reserve(in + ik);

//   double ta = parvals_[0];
//   double tb = parvals_[parvals_.size() - 1];
  double ta = evalcrv_->start();
  double tb = evalcrv_->end();
  double tint = (tb - ta)/(double)(in-ik+1);
  double tdiff;

  // Define knot vector
  int ki;
  for (ki=0; ki<ik; ki++)
    knots.push_back(ta);
  for (ki=ik, tdiff=tint; ki<in; ki++, tdiff+=tint)
    knots.push_back(ta+tdiff);
  for (ki=0; ki<ik; ki++)
    knots.push_back(tb);

  // Define dummy coefficient vector
  std::vector<double> coefs(in*dim_, 0.0);

  // Make initial curve
  curr_crv_ = shared_ptr<SplineCurve>
    ( new SplineCurve(in, ik, &knots[0], &coefs[0], dim_));
 
}


//***************************************************************************

void AdaptCurve::adjustSmoothWeight()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Check distribution of knots and parameter values. If
   //               necessary increase the smoothing weight
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
  // Compute the average number of parameter values for each knot interval
  int ki, kj;
  int kn = curr_crv_->numCoefs();
  int kk = curr_crv_->order();
  std::vector<double>::const_iterator st = curr_crv_->basis().begin();
  int knpar = (int)parvals_.size();
  int minnmb = knpar;
  int currnmb=0;
  double avnmb = (double)knpar/(double)(kn-kk+1);
  for (kj=0, ki=kk; ki<kn; ki++)
    {
      currnmb = 0;
      while (kj<knpar && parvals_[kj] >= st[ki-1] && parvals_[kj] < st[ki])
	{
	  currnmb++;
	  kj++;
	}
      minnmb = std::min(minnmb, currnmb);
    }

  if (minnmb < 1 || (double)minnmb < avnmb/2.5)
    {
    smoothweight_ = 0.000001;
    if (getenv("DEBUG"))
      {
	std::cout << " Weight: " << smoothweight_ << " min = " << minnmb <<
	    " average " << avnmb << std::endl;
      }
    }
}


  //***************************************************************************

void AdaptCurve::makeSmoothCurve()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Adapt the initial curve to approximate the given
   //               data points.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
  SmoothCurve crvgen(dim_);

  int bd_fix[2];
  bd_fix[0] = bd_fix[1] = 1;
  double wgt1 = 0.0; //0.000001*smoothweight_;
  double wgt2 = smoothweight_/smoothfac_;
  double wgt3 = smoothweight_/(smoothfac_*smoothfac_);
  double approxweight = 1.0 - wgt1 - wgt2 - wgt3;
  vector<int> coef_known(curr_crv_->numCoefs(), 0);
  vector<sideConstraint> all_constraints, pt_constraints(2), tangent_constraints(2);
  getConstraints(pt_constraints, tangent_constraints);
  // We set start pt of approximated curve.
  if (start_pt_.size() == 0) // || (start_pt_[0].get() == 0))
    coef_known[0] =  1; // We lock start pt when approximating.
  else
    all_constraints.push_back(pt_constraints[0]);
  if (end_pt_.size() == 0) // || (end_pt_[0].get() == 0))
    coef_known[coef_known.size()-1] =  1; // We lock end pt when approximating.
  else
    all_constraints.push_back(pt_constraints[1]);

  // We next add possible end tangents.
  if (start_pt_.size() > 1) // && start_pt_[1].get() != 0)
      all_constraints.push_back(tangent_constraints[0]);    
  if (end_pt_.size() > 1) // && end_pt_[1].get() != 0)
      all_constraints.push_back(tangent_constraints[1]);

//   crvgen.attach(curr_crv_, bd_fix);
  //crvgen.attach(curr_crv_, 0, &coef_known[0], all_constraints.size());
  crvgen.attach(curr_crv_, &coef_known[0], (int)all_constraints.size());

  crvgen.setOptim(wgt1, wgt2, wgt3);

  crvgen.setLeastSquares(points_, parvals_, pt_weight_, approxweight);

  if (all_constraints.size() != 0)
    crvgen.setSideConstraints(all_constraints);

  // DEBUG output
#ifdef DEBUG
  std::ofstream out_file("crv_out.g2");
  curr_crv_->writeStandardHeader(out_file);
  curr_crv_->write(out_file);
  out_file << "400 1 0 4 255 0 0 255" << std::endl;
  out_file << parvals_.size() << std::endl;
  for (size_t ki=0; ki<points_.size(); ki+=3)
    out_file << points_[ki] << "  " << points_[ki+1] << "  " << points_[ki+2] << std::endl;
#endif

  // We then are ready to solve system.
  crvgen.equationSolve(curr_crv_);

#ifdef DEBUG
  // DEBUG output
  curr_crv_->writeStandardHeader(out_file);
  curr_crv_->write(out_file);
#endif

  return;
}

//***************************************************************************

void AdaptCurve::checkAccuracy(std::vector<double>& newknots, int uniform)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Check the accuracy of the current curve compared to
   //               the given data points.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
//     // debug
//     std::ofstream debug("data/debug.g2");
//     SplineDebugUtils::writeSpaceParamCurve(*curr_crv_, debug);
//     // end of debug

  bool reparam = false;

//     double par_tol = 0.000000000001;
    maxdist_ = -10000.0;
    maxdist_inpoints_= -10000.0;
    avdist_ = 0.0;
    int nmbdist = 0;

    // Traverse the data points and check the quality of the curve
    // approximation
    int nmbpnt = (int)parvals_.size();
    bool do_insert = (nmbpnt < 4*init_sample_);  // If the data set should
    // be increased
    double dist1, dist2;
    int left = 0, prevleft = 0;
    vector<double>::const_iterator st = curr_crv_->basis().begin();
    double t1 = curr_crv_->startparam(); 
    double t2 = curr_crv_->endparam();
    double t3 = evalcrv_->start();
    double t4 = evalcrv_->end();
    double parfac = (t4 - t3)/(t2 - t1);
    double ta, tb;
    double par;
    int distOK = 1;
    vector<double> newpar;
    vector<Point> newpoints;
    vector<double> newweight;
    double newknot;
    double frac = 0, pardist = 0;
    Point pos(dim_);
    Point clpos(dim_);
    double *pt, *param;
    int ki, kj;
    for (ki=0, pt=&points_[0], param=&parvals_[0]; ki<nmbpnt; 
	 ki++, pt+=dim_, param++) {
	// Compute closest point
	pos.setValue(pt);
	try {
	    curr_crv_->closestPoint(pos, t1, t2, 
				    par, clpos, dist1, param);
	} catch (...) {
	    MESSAGE("Failed finding closest point.");
	    par = *param; // We're using our input value.
	    // This should at least be an upper boundary of closest dist.
	    dist1 = (pos - curr_crv_->ParamCurve::point(par)).length();
	}
	left = curr_crv_->basis().knotInterval(par);
	if (reparam)
	    parvals_[ki] = par;
	    

	// Statistics
	if (parvals_[ki] >= min_sample_par_ && parvals_[ki] <= max_sample_par_)
	  {
	    maxdist_inpoints_ = std::max(maxdist_inpoints_, dist1);
	    maxdist_ = std::max(maxdist_, dist1);
	    avdist_ += dist1;
	    nmbdist++;
	  }

	// Check the middle point between this and the succeeding
	// sample point
	if (ki < nmbpnt-1)
	  {
	    ta = 0.5*(param[0] + param[1]);
	    
	    // Since the parameter interval of the approximating curve may
	    // differ from that of the evaluator based curve, compute
	    // the parameter corresponding to the evaluator based curve
	    tb = t3 + (ta - t1)*parfac;
	    
	    //pos = evalcrv_->eval(tb);
	    pos = evalcrv_->eval(ta);

	    // Compute closest point
	    try {
	      curr_crv_->closestPoint(pos, param[0], param[1],
				      par, clpos, dist2, &ta);
	    } catch (...) {
	      MESSAGE("Failed finding closest point.");
	      par = ta; // We're using our input value.
	      // This should at least be an upper boundary of closest dist.
	      dist2 = (pos - curr_crv_->ParamCurve::point(par)).length();
	    }
	    left = curr_crv_->basis().knotInterval(par);

	    // Statistics
	    double tc = (reparam) ? par : ta;
	    if (tc >= min_sample_par_ && tc <= max_sample_par_)
	      {
		maxdist_ = std::max(maxdist_, dist2);
		avdist_ += dist2;
		nmbdist++;
	      }

	    if (do_insert)
	      {
		// Add the new point to the sample points
		// Remember values
		newpar.push_back(tc);
		newpoints.push_back(pos);
		newweight.push_back((tb<min_sample_par_ || tb>max_sample_par_)
				    ? 0.1 : 1.0);
	      }
	  }
	else 
	  dist2 = 0.0;

	// Check if the accuracy is good enough in this point
	if (dist1 > aepsge_ || dist2 > aepsge_) {
	    // Necessary to refine the spline space. Find the knot
	    // interval to refine.

	    distOK = 0;
	    if (left > prevleft) {
		if (frac > 0.0) {
		    // Set previous new knot
		    ta = st[prevleft];
		    tb = st[prevleft+1];

		    newknot = (uniform) ? 0.5*(ta+tb) :
			st[prevleft] + pardist/frac;

		    if (newknot-ta < 0.1*(tb-ta))
			newknot = ta + 0.1*(tb-ta);

		    if (tb-newknot < 0.1*(tb-ta))
			newknot = tb - 0.1*(tb-ta);

		    newknots.push_back(newknot);
		}

		// Reset parameters
		prevleft = left;
		pardist = 0.0;
		frac = 0.0;
	    }

	    // Start computing the next new knot
	    pardist += (par-st[left])*(dist1+dist2);
	    frac += (dist1 + dist2);
	}
    }

    if (!distOK && frac > 0) {
	// Add the last new knot

	ta = st[prevleft];
	tb = st[prevleft+1];

	newknot = (uniform) ? 0.5*(ta+tb) :
	    st[prevleft] + pardist/frac;

	if (newknot-ta < 0.1*(tb-ta))
	    newknot = ta + 0.1*(tb-ta);

	if (tb-newknot < 0.1*(tb-ta))
	    newknot = tb - 0.1*(tb-ta);

	newknots.push_back(newknot);
	    
    }
    avdist_ /= (double)(nmbdist);

    if (newpar.size() > 0)
      {
	// Insert the new data points
	vector<double> tmp_par;
	vector<double> tmp_pnt;
	vector<double> tmp_weight;
	tmp_par.reserve(parvals_.size() + newpar.size());
	tmp_pnt.reserve(points_.size() + dim_*newpoints.size());
	tmp_weight.reserve(pt_weight_.size() + newweight.size());
	for (ki=0, kj=0; ki<(int)parvals_.size() && kj<(int)newpar.size(); )
	  {
	    if (parvals_[ki] < newpar[kj])
	      {
		tmp_par.push_back(parvals_[ki]);
		tmp_pnt.insert(tmp_pnt.end(), points_.begin()+ki*dim_,
			       points_.begin() + (ki+1)*dim_);
		tmp_weight.push_back(pt_weight_[ki]);
		ki++;
	      }
	    else
	      {
		// It should be no coincidence between values
		tmp_par.push_back(newpar[kj]);
		tmp_pnt.insert(tmp_pnt.end(), newpoints[kj].begin(),
			       newpoints[kj].end());
		tmp_weight.push_back(newweight[kj]);
		kj++;
	      }
	  }
	if (ki<(int)parvals_.size())
	  {
	    tmp_par.insert(tmp_par.end(), parvals_.begin()+ki, parvals_.end());
	    tmp_pnt.insert(tmp_pnt.end(), points_.begin()+ki*dim_,
			   points_.end());
	    tmp_weight.insert(tmp_weight.end(), pt_weight_.begin()+ki,
			      pt_weight_.end());
	  }
	else if (kj<(int)newpar.size())
	  {
	    tmp_par.insert(tmp_par.end(), newpar.begin()+kj, newpar.end());
	    for (; kj<(int)newpar.size(); kj++)
	      tmp_pnt.insert(tmp_pnt.end(), newpoints[kj].begin(),
			     newpoints[kj].end());
	    tmp_weight.insert(tmp_weight.end(), newweight.begin()+kj,
			      newweight.end());
	  }

	parvals_ = tmp_par;
	points_ = tmp_pnt;
	pt_weight_ = tmp_weight;
      }
}

//***************************************************************************

int AdaptCurve::approximate(int max_iter)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Generate the approximating curve
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
  if (getenv("DEBUG"))
    {
      std::cout << "crv" << std::endl;
    }

#ifdef DEBUG
  std::ofstream of("adapt_crv.g2");
  curr_crv_->writeStandardHeader(of);
  curr_crv_->write(of);
  evalcrv_->write(of);
#endif
  
  // Make initial sample points
  initSamples();

  // Check distribution of knots and parameter values. If
  // necessary increase the smoothing weight  
  adjustSmoothWeight();

  // Make the initial approximation. First fix endpoints.
  int nmbpoints = (int)points_.size()/dim_;
  std::vector<double>::iterator coef = curr_crv_->coefs_begin();
  int in = curr_crv_->numCoefs();
  int ki;
  for (ki=0; ki<dim_; ki++)
    {
      if (fix_[0] == 0)
	coef[ki] = points_[ki];
      if (fix_[1] == 0)
	coef[(in-1)*dim_+ki] = 
	  points_[(nmbpoints-1)*dim_+ki];
    }

  // Approximate
  makeSmoothCurve();

  if (max_iter == 0) { // In order to set maxerr_ & meanerr_.
      std::vector<double> newknots_dummy;
      checkAccuracy(newknots_dummy, true); // @@sbr Using uniform knot insertion.
  }

  // Check the accuracy and iterate the approximation
  prev_maxdist_ =   100000.0;
  prev_avdist_ = 100000.0;
  for (ki=0; ki<max_iter; ki++)
    {
      std::vector<double> newknots;

//       stat = checkAccuracy(newknots, (ki <= 4));
      checkAccuracy(newknots, true); // @@sbr Using uniform knot insertion.

      if (getenv("DEBUG"))
	{
	  std::cout << " # pnts " << nmbpoints << " # coef " << in << " max "
		    << maxdist_ << " average " << avdist_ <<  std::endl;
	}

      if (maxdist_ <= aepsge_ || newknots.size() == 0)
	break;   // The required accuracy is reached.

      if (maxdist_inpoints_ > 0.95*prev_maxdist_) {
	  MESSAGE("Convergence not too promising...");
// 	break;   // Not enough gain in refining
      }

      if (maxdist_inpoints_ > prev_maxdist_ && ki > 1)
	{
	  // Stop iteration. Return to previous result
	  curr_crv_ = prev_crv_;
	  maxdist_ = prev_maxdist_;
	  avdist_ = prev_avdist_;
	  break;
	}

      if (maxdist_ < prev_maxdist_)
	{
	  prev_maxdist_ = maxdist_;
	  prev_avdist_ = avdist_;
	  prev_crv_ = shared_ptr<SplineCurve>(curr_crv_->clone());
	}

      // Refine the spline space
      curr_crv_->insertKnot(newknots);
 
      // Approximate the data by a curve in the new spline space.
      makeSmoothCurve();
    }

  return (ki == max_iter);
}

//***************************************************************************

void AdaptCurve::getConstraints(std::vector<sideConstraint>& pt_constraints,
				 std::vector<sideConstraint>& tangent_constraints)
{
  ALWAYS_ERROR_IF(pt_constraints.size() != 2 || tangent_constraints.size() != 2,
	      "Input vectors are not of size 2.");

  // We start by setting pt_constraints.
  if (start_pt_.size() != 0) { // && (start_pt_[0].get() != 0)) {
    sideConstraint constraint;
    constraint.dim_ = curr_crv_->dimension();
    constraint.factor_.push_back(make_pair(0, 1.0));
    for (int ki = 0; ki < 3; ++ki)
      constraint.constant_term_[ki] = (start_pt_[0])[ki];
    pt_constraints[0] = constraint;
  }
  if (end_pt_.size() != 0) { // && (end_pt_[0].get() != 0)) {
    sideConstraint constraint;
    constraint.dim_ = curr_crv_->dimension();
    constraint.factor_.push_back(make_pair(curr_crv_->numCoefs() - 1, 1.0));
    for (int ki = 0; ki < 3; ++ki)
      constraint.constant_term_[ki] = (end_pt_[0])[ki];
    pt_constraints[1] = constraint;
  }

  // We then set tangent_constraints.
  if (start_pt_.size() > 1) { //&& (start_pt_[1].get() != 0)) {
    sideConstraint constraint;
    constraint.dim_ = curr_crv_->dimension();
    constraint.factor_.push_back(make_pair(1, 1.0));
    constraint.factor_.push_back(make_pair(0, -1.0));
    double t_inc = curr_crv_->basis().grevilleParameter(1) -
      curr_crv_->basis().grevilleParameter(0);
    for (int ki = 0; ki < 3; ++ki)
      constraint.constant_term_[ki] = ((start_pt_[1])[ki])*t_inc;
    tangent_constraints[0] = constraint;
  }
  if (end_pt_.size() > 1) { // && (end_pt_[1].get() != 0)) {
    sideConstraint constraint;
    constraint.dim_ = curr_crv_->dimension();
    constraint.factor_.push_back(make_pair(curr_crv_->numCoefs() - 1, 1.0));
    constraint.factor_.push_back(make_pair(curr_crv_->numCoefs() - 2, -1.0));
    double t_inc = curr_crv_->basis().grevilleParameter(curr_crv_->numCoefs() - 1) -
      curr_crv_->basis().grevilleParameter(curr_crv_->numCoefs() - 2);
    for (int ki = 0; ki < 3; ++ki)
      constraint.constant_term_[ki] = ((end_pt_[1])[ki])*t_inc;
    tangent_constraints[1] = constraint;
  }
}

//***************************************************************************
void AdaptCurve::setEndPoints(const std::vector<Point>& start_pt,
			       const std::vector<Point>& end_pt)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Set end tangents to be used when approximating points.
   //
   //     Calls   :
   //
   //     Written by : Sverre Briseid,  SINTEF,  02-06
   //--------------------------------------------------------------------------
{
  start_pt_ = start_pt;
  end_pt_ = end_pt;
  double ta = evalcrv_->start();
  double tb = evalcrv_->end();
  double fac = 0.0001;
  min_sample_par_ += fac*(tb-ta);
  max_sample_par_ -= fac*(tb-ta);
}

//***************************************************************************

shared_ptr<SplineCurve> AdaptCurve::getAdaptCurve(double& maxdist, 
							double& avdist,
							int max_iter)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Return the approximating curve and parameters
   //               describing the accuracy.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  09-11
   //--------------------------------------------------------------------------
{
    // Set accuracy parameters
    maxdist = maxdist_;
    avdist = avdist_;

    // Return curve

    return curr_crv_;
}

//***************************************************************************

void  AdaptCurve::unsetSmooth()
{
  smoothweight_ = 0.0;
}

//***************************************************************************
void AdaptCurve::setSmooth(double w)
{
   if (w>0.0 && w<1.0)
      smoothweight_=w;
}

void AdaptCurve::initSamples()
//-------------------------------------------------------------------------
//
// PURPOSE: Create parameter values and sample points for the first time
//
//
//-------------------------------------------------------------------------
{
  // Count number of spline segments
  int ki, kj;
  int kn = curr_crv_->numCoefs();
  vector<double>::iterator st = curr_crv_->knotsBegin();
  int nseg=0;
  for (ki=order_; ki<=kn; ki++)
    {
      if (st[ki] > st[ki-1])
	nseg++;
    }

  int nprsample = 20;
  //int nsample = std::max(std::min(nseg*nprsample, 1000),30);
  int nsample = std::max(std::min(nseg*nprsample, 2000),500);

  // Make sure that the data set arrays are empty
  init_sample_ = nsample;
  parvals_.clear();
  points_.clear();
  parvals_.reserve(nsample);
  points_.reserve(nsample*dim_);

  // Compute sample points
  int left, right;
  double ta, tb, tint;
  Point pnt;
  double weight;
  double t1 = curr_crv_->startparam();
  double t2 = curr_crv_->endparam();
  double t3 = evalcrv_->start();
  double t4 = evalcrv_->end();
  double parfac = (t4 - t3)/(t2 - t1);
  int dim = curr_crv_->dimension();

  // Point at start
  ta = t1;
  tb = t3 + (ta - t1)*parfac;

  weight = 1.0;
  if (tb < min_sample_par_)
    weight = 0.1;
  if (tb > max_sample_par_)
    weight = 0.1;
  //pnt = evalcrv_->eval(tb);
  pnt = evalcrv_->eval(ta);
  Point prevpt = pnt;
  // TESTING
  Point cv_pt = curr_crv_->ParamCurve::point(ta);

  parvals_.push_back(ta);
  points_.insert(points_.end(), pnt.begin(), pnt.end());
  
  pt_weight_.push_back(weight);

  double len = 0.0;
  for (kj=0, left=order_-1; kj<nseg; kj++)
    {
      for (right=left+1; right<kn && st[right]==st[left]; right++);
      nprsample = (int)((double)(nsample)*(st[right]-st[left])/(t2-t1));
      nprsample = std::max(nprsample,2);
      tint = (st[right] - st[left])/(double)(nprsample);
      for (ki=0, ta=st[left]+0.5*tint; ki<nprsample; ki++, ta+=tint)
	{
	  // Since the parameter interval of the approximating curve may
	  // differ from that of the evaluator based curve, compute
	  // the parameter corresponding to the evaluator based curve
	  tb = t3 + (ta - t1)*parfac;

	  weight = 1.0;
	  if (tb < min_sample_par_)
	    weight = 0.1;
	  if (tb > max_sample_par_)
	    weight = 0.1;
	  //pnt = evalcrv_->eval(tb);
	  pnt = evalcrv_->eval(ta);
	  len += prevpt.dist(pnt);
	  prevpt = pnt;

	  // TESTING
	  cv_pt = curr_crv_->ParamCurve::point(ta);

	  parvals_.push_back(ta);
	  points_.insert(points_.end(), pnt.begin(), pnt.end());

	  pt_weight_.push_back(weight);
	}
      left = right;
    }

  // Point at end
  ta = t2;
  tb = t3 + (ta - t1)*parfac;

  weight = 1.0;
  if (tb < min_sample_par_)
    weight = 0.1;
  if (tb > max_sample_par_)
    weight = 0.1;
  //pnt = evalcrv_->eval(tb);
  pnt = evalcrv_->eval(ta);
  // TESTING
  cv_pt = curr_crv_->ParamCurve::point(ta);

  parvals_.push_back(ta);
  points_.insert(points_.end(), pnt.begin(), pnt.end());
  
  pt_weight_.push_back(weight);
  len += prevpt.dist(pnt);
  smoothfac_ = 1.0/len;
 }

