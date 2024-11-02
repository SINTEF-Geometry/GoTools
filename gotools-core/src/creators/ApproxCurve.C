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

#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/utils/Point.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineUtils.h" // debugging

using namespace Go;
using std::vector;
using std::max;
using std::min;
using std::make_pair;

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
#include "GoTools/geometry/Utils.h"     // make std::min and std::max work (redefined in boost/smart_ptr.hpp)
#endif

ApproxCurve::ApproxCurve()
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  maxdist_ = -10000.0;
  avdist_ = 0;
  dim_ = 3;
  aepsge_ = 0.01;
  smoothweight_ = 0.000000001;
  smoothfac_ = 1.0;
  c1fac_ = 0.0;
}

//***************************************************************************

ApproxCurve::ApproxCurve(const std::vector<double>& points, 
			     const std::vector<double>& parvals,
			     int dim, double aepsge)
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  maxdist_ = -10000.0;
  avdist_ = 0;
  dim_ = dim;
  aepsge_ = aepsge;
  smoothweight_ = 0.000000001;
  c1fac_ = 0.0;

  points_.reserve(points.size());
  parvals_.reserve(parvals.size());

  for (size_t ki=0; ki < points.size(); ki++)
    points_.push_back(points[ki]);

  for (size_t ki=0; ki<parvals.size(); ki++)
    parvals_.push_back(parvals[ki]);

  smoothfac_ = 1.0/(parvals_[parvals_.size()-1] - parvals_[0]);

  makeInitCurve();
}

//***************************************************************************

ApproxCurve::ApproxCurve(const std::vector<double>& points, 
			     const std::vector<double>& parvals, int dim,
			     double aepsge, int in, int ik,
			     const std::vector<double>& knots)
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  maxdist_ = -10000.0;
  avdist_ = 0;
  dim_ = dim;
  aepsge_ = aepsge;
  smoothweight_ = 0.000000001;
  c1fac_ = 0.0;

  points_.reserve(points.size());
  parvals_.reserve(parvals.size());
  for (size_t ki=0; ki<points.size(); ki++)
    points_.push_back(points[ki]);

  for (size_t ki=0; ki<parvals.size(); ki++)
    parvals_.push_back(parvals[ki]);

  smoothfac_ = 1.0/(parvals_[parvals_.size()-1] - parvals_[0]);

  makeInitCurve(in, ik, knots);
}

//***************************************************************************

ApproxCurve::ApproxCurve(const std::vector<double>& points, 
			     const std::vector<double>& parvals, int dim,
			     double aepsge, int in, int ik)
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  maxdist_ = -10000.0;
  avdist_ = 0;
  dim_ = dim;
  aepsge_ = aepsge;
  smoothweight_ = 0.000000001;
  c1fac_ = 0.0;

  points_.reserve(points.size());
  parvals_.reserve(parvals.size());
  for (size_t ki=0; ki<points.size(); ki++)
    points_.push_back(points[ki]);

  for (size_t ki=0; ki<parvals.size(); ki++)
    parvals_.push_back(parvals[ki]);

  smoothfac_ = 1.0/(parvals_[parvals_.size()-1] - parvals_[0]);

  makeInitCurve(in, ik);
}

//***************************************************************************

ApproxCurve::~ApproxCurve()
   //--------------------------------------------------------------------------
   //     Destructor for class ApproxCurve.
   //
   //     Purpose : 
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
//   if (curr_crv_)
//     freeCurve(curr_crv_);
}

//***************************************************************************

void ApproxCurve::makeInitCurve()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
    int nmbpar = (int)parvals_.size();

  int ik = 4;  // Cubic curve by default
  int in = nmbpar/6; 
  in = std::min(std::max(in, 4), 10);  // Minimum of 4 and maximum of 10 coeffients.
                               // Otherwise depending on the number of points.

  makeInitCurve(in, ik);

  return;
}


//***************************************************************************

void ApproxCurve::makeInitCurve(int in, int ik, std::vector<double> knots)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  // Define dummy coefficient vector
  std::vector<double> coefs(in*dim_, 0.0);

  // Make initial curve
  curr_crv_ = shared_ptr<SplineCurve>
    ( new SplineCurve(in, ik, &knots[0], &coefs[0], dim_));
 
  // Check distribution of knots and parameter values. If
  // necessary increase the smoothing weight
  adjustSmoothWeight();

  return;
}


//***************************************************************************

void ApproxCurve::makeInitCurve(int in, int ik)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  std::vector<double> knots;
  knots.reserve(in + ik);

  double ta = parvals_[0];
  double tb = parvals_[parvals_.size() - 1];
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
 
  // Check distribution of knots and parameter values. If
  // necessary increase the smoothing weight
  adjustSmoothWeight();

  return;
}


//***************************************************************************

void ApproxCurve::adjustSmoothWeight()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Check distribution of knots and parameter values. If
   //               necessary increase the smoothing weight
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
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
        MESSAGE("INFO: Weight: " << smoothweight_ << " min = " << minnmb << " average " << avnmb);
    }
}


  //***************************************************************************

void ApproxCurve::makeSmoothCurve()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Adapt the initial curve to approximate the given
   //               data points.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  SmoothCurve crvgen(dim_);

  int bd_fix[2];
  bd_fix[0] = bd_fix[1] = 1;
  double wgt1 = 0.000001*smoothweight_;
  double wgt2 = smoothweight_/smoothfac_;
  double wgt3 = smoothweight_/(smoothfac_*smoothfac_);
  double approxweight = 1.0 - wgt1 - wgt2 - wgt3;
  std::vector<double> pt_weight(parvals_.size(), 1.0);
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

  crvgen.setLeastSquares(points_, parvals_, pt_weight, approxweight);

  if (all_constraints.size() != 0)
    crvgen.setSideConstraints(all_constraints);

  if (c1fac_ > 0.0)
    crvgen.setPeriodicity(1, c1fac_);

  // We then are ready to solve system.
  crvgen.equationSolve(curr_crv_);

  return;
}

//***************************************************************************

void ApproxCurve::checkAccuracy(std::vector<double>& newknots, int uniform)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Check the accuracy of the current curve compared to
   //               the given data points.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
//     // debug
//     std::ofstream debug("data/debug.g2");
//     SplineDebugUtils::writeSpaceParamCurve(*curr_crv_, debug);
//     // end of debug

#ifdef DEBUG
  std::ofstream of("curr_cv.g2");
  curr_crv_->writeStandardHeader(of);
  curr_crv_->write(of);
#endif
  bool reparam = (dim_ == 1) ? false : true;

//     double par_tol = 0.000000000001;
    maxdist_ = -10000.0;
    avdist_ = 0.0;

    // Traverse the data points and check the quality of the curve
    // approximation
    int nmbpnt = (int)parvals_.size();
    double dist;
    int left = 0, prevleft = 0;
    std::vector<double>::const_iterator st = curr_crv_->basis().begin();
    double ta, tb;
    double par;
    int distOK = 1;
    double newknot;
    double frac = 0, pardist = 0;
    Point pos(dim_);
    Point clpos(dim_);
    double *pt, *param;
    int ki;
    for (ki=0, pt=&points_[0], param=&parvals_[0]; ki<nmbpnt; 
	 ki++, pt+=dim_, param++)
      {
	if (reparam || dim_ > 1)
	  {
	    // Compute closest point
	    ta = curr_crv_->startparam(); 
	    tb = curr_crv_->endparam();
	    pos.setValue(pt);
	    try {
	      curr_crv_->closestPoint(pos, ta, tb, 
				      par, clpos, dist, param);
	    } catch (...) {
	      MESSAGE("Failed finding closest point.");
	      par = *param; // We're using our input value.
	      // This should at least be an upper boundary of closest dist.
	      dist = (pos - curr_crv_->ParamCurve::point(par)).length();
	    }
	  }
	else
	  {
	    par = *param; 
	    dist = (pos - curr_crv_->ParamCurve::point(par)).length();
	  }
	left = curr_crv_->basis().knotInterval(par);
	if (reparam)
	    parvals_[ki] = par;
	    

	// Statistics
	maxdist_ = std::max(maxdist_, dist);
	avdist_ += dist;

	// Check if the accuracy is good enough in this point
	if (dist > aepsge_) {
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
	    pardist += (par-st[left])*dist;
	    frac += dist;
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
    avdist_ /= (double)nmbpnt;
	
    return;
}

//***************************************************************************

int ApproxCurve::doApprox(int max_iter)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Generate the approximating curve
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  // Make the initial approximation. First fix endpoints.
  int nmbpoints = (int)points_.size()/dim_;
  std::vector<double>::iterator coef = curr_crv_->coefs_begin();
  int in = curr_crv_->numCoefs();
  int ki;
  for (ki=0; ki<dim_; ki++)
    {
      coef[ki] = points_[ki];
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
  double prevmax = 100000.0, prevav = 100000.0;
  for (ki=0; ki<max_iter; ki++)
    {
      std::vector<double> newknots;

//       stat = checkAccuracy(newknots, (ki <= 4));
      checkAccuracy(newknots, true); // @@sbr Using uniform knot insertion.

      //MESSAGE("crv # pnts " << nmbpoints << " # coef " << in << " max " << maxdist_ << " average " << avdist_);

      if (maxdist_ <= aepsge_ || newknots.size() == 0)
	break;   // The required accuracy is reached.

      if (maxdist_ > 0.95*prevmax) {
	//MESSAGE("Convergence not too promising...");
// 	break;   // Not enough gain in refining
      }

      prevmax = maxdist_;
      prevav = avdist_;

      // Refine the spline space
      curr_crv_->insertKnot(newknots);
 
      // Approximate the data by a curve in the new spline space.
      makeSmoothCurve();
    }

  return (ki == max_iter);
}

//***************************************************************************

void ApproxCurve::getConstraints(std::vector<sideConstraint>& pt_constraints,
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
void ApproxCurve::setEndPoints(const std::vector<Point>& start_pt,
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
}

//***************************************************************************

shared_ptr<SplineCurve> ApproxCurve::getApproxCurve(double& maxdist, 
							double& avdist,
							int max_iter)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Return the approximating curve and parameters
   //               describing the accuracy.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  00-04
   //--------------------------------------------------------------------------
{
  // Generate the approximating curve

    /*int stat =*/ doApprox(max_iter);

    // Set accuracy parameters
    maxdist = maxdist_;
    avdist = avdist_;

    // Return curve

    return curr_crv_;
}

//***************************************************************************

void  ApproxCurve::unsetSmooth()
{
  smoothweight_ = 0.0;
}

//***************************************************************************
void ApproxCurve::setSmooth(double w)
{
   if (w>0.0 && w<1.0)
      smoothweight_=w;
}

//***************************************************************************

void  ApproxCurve::setC1Approx(double fac)
{
  c1fac_ = fac;
}
