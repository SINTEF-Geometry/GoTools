// ----------------------------------------------------------------
//       Implementation file for class ApproxCrvToSeqs.
// ----------------------------------------------------------------
//
// This file contains the following member functions:
//             
// ----------------------------------------------------------------

#include "GoTools/creators/ApproxCrvToSeqs.h"
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/utils/Point.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineUtils.h" // debugging

using namespace Go;
using std::vector;

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
#include "GoTools/geometry/Utils.h"     // make std::min and std::max work (redefined in boost/smart_ptr.hpp)
#endif

ApproxCrvToSeqs::ApproxCrvToSeqs()
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxCrvToSeqs.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  maxdist_ = -10000.0;
  avdist_ = 0;
  dim_ = 3;
  aepsge_ = 0.01;
  smoothweight_ = 0.000000001;
  smoothfac_ = 1.0;
  c1fac_ = 0.0;
  a1_ = a2_ = 0.0;
}

//***************************************************************************

ApproxCrvToSeqs::ApproxCrvToSeqs(const vector<vector<double> >& points, 
				 const vector<vector<double> >& parvals,
				 int dim, double aepsge)
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxCrvToSeqs.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  maxdist_ = -10000.0;
  avdist_ = 0;
  dim_ = dim;
  aepsge_ = aepsge;
  smoothweight_ = 0.000000001;
  c1fac_ = 0.0;

  start_pt_.resize(parvals.size());
  end_pt_.resize(parvals.size());

  points_ = points;
  parvals_.resize(parvals.size());
  for (size_t ki=0; ki<parvals.size(); ++ki)
    parvals_[ki].insert(parvals_[ki].end(), parvals[ki].begin(), 
			parvals[ki].end());

  a1_ = parvals[0][0];
  a2_ = parvals[0][parvals_[0].size()-1];
  for (size_t ki=1; ki<parvals_.size(); ++ki)
    {
      a1_ = std::min(a1_, parvals[ki][0]);
      a2_ = std::max(a2_, parvals[ki][parvals_[ki].size()-1]);
    }
  smoothfac_ = 1.0/(a2_ - a1_);

  makeInitCurves();
}

//***************************************************************************

ApproxCrvToSeqs::ApproxCrvToSeqs(const vector<vector<double> >& points, 
				 const vector<vector<double> >& parvals, 
				 int dim, double aepsge, int in, int ik,
				 vector<double>& knots)
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxCrvToSeqs.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  maxdist_ = -10000.0;
  avdist_ = 0;
  dim_ = dim;
  aepsge_ = aepsge;
  smoothweight_ = 0.000000001;
  c1fac_ = 0.0;

  start_pt_.resize(parvals.size());
  end_pt_.resize(parvals.size());

  points_ = points;
  parvals_.resize(parvals.size());
  for (size_t ki=0; ki<parvals.size(); ++ki)
    parvals_[ki].insert(parvals_[ki].end(), parvals[ki].begin(), 
			parvals[ki].end());

  a1_ = parvals[0][0];
  a2_ = parvals[0][parvals_[0].size()-1];
  for (size_t ki=1; ki<parvals_.size(); ++ki)
    {
      a1_ = std::min(a1_, parvals[ki][0]);
      a2_ = std::max(a2_, parvals[ki][parvals_[ki].size()-1]);
    }
   smoothfac_ = 1.0/(a2_ - a1_);

  makeInitCurves(in, ik, knots);
}

//***************************************************************************

ApproxCrvToSeqs::ApproxCrvToSeqs(const vector<vector<double> >& points, 
				 const vector<vector<double> >& parvals, 
				 int dim, double aepsge, int in, int ik)
   //--------------------------------------------------------------------------
   //     Constructor for class ApproxCrvToSeqs.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  maxdist_ = -10000.0;
  avdist_ = 0;
  dim_ = dim;
  aepsge_ = aepsge;
  smoothweight_ = 0.000000001;
  c1fac_ = 0.0;

  start_pt_.resize(parvals.size());
  end_pt_.resize(parvals.size());

  points_ = points;
  parvals_.resize(parvals.size());
  for (size_t ki=0; ki<parvals.size(); ++ki)
    parvals_[ki].insert(parvals_[ki].end(), parvals[ki].begin(), 
			parvals[ki].end());

  a1_ = parvals[0][0];
  a2_ = parvals[0][parvals_[0].size()-1];
  for (size_t ki=1; ki<parvals_.size(); ++ki)
    {
      a1_ = std::min(a1_, parvals[ki][0]);
      a2_ = std::max(a2_, parvals[ki][parvals_[ki].size()-1]);
    }
   smoothfac_ = 1.0/(a2_ - a1_);

  makeInitCurves(in, ik);
}

//***************************************************************************

ApproxCrvToSeqs::~ApproxCrvToSeqs()
   //--------------------------------------------------------------------------
   //     Destructor for class ApproxCrvToSeqs.
   //
   //     Purpose : 
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
//   if (curr_crv_)
//     freeCurve(curr_crv_);
}

//***************************************************************************

void ApproxCrvToSeqs::makeInitCurves()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
    int nmbpar = (int)parvals_[0].size();
  for (size_t ki=0; ki<parvals_.size(); ++ki)
    nmbpar = std::max(nmbpar, (int)parvals_[ki].size());

  int ik = 4;  // Cubic curve by default
  int in = nmbpar/6; 
  in = std::min(std::max(in, 4), 10);  // Minimum of 4 and maximum of 10 coeffients.
                               // Otherwise depending on the number of points.

  makeInitCurves(in, ik);

  return;
}


//***************************************************************************

void ApproxCrvToSeqs::makeInitCurves(int in, int ik, vector<double> knots)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  // Define dummy coefficient vector
  vector<double> coefs(in*dim_, 0.0);

  // Make initial curve
  curr_crv_.resize(parvals_.size());
  for (size_t ki=0; ki<parvals_.size(); ++ki)
    curr_crv_[ki] = shared_ptr<SplineCurve>
      ( new SplineCurve(in, ik, &knots[0], &coefs[0], dim_));
 
  // Check distribution of knots and parameter values. If
  // necessary increase the smoothing weight
  adjustSmoothWeight();

  return;
}


//***************************************************************************

void ApproxCrvToSeqs::makeInitCurves(int in, int ik)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Make an initial spline space generating a uniform
   //               knot-vector on the given parameter interval.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  vector<double> knots;
  knots.reserve(in + ik);

  double tint = (a2_ - a1_)/(double)(in-ik+1);
  double tdiff;

  // Define knot vector
  int ki;
  for (ki=0; ki<ik; ki++)
    knots.push_back(a1_);
  for (ki=ik, tdiff=tint; ki<in; ki++, tdiff+=tint)
    knots.push_back(a1_+tdiff);
  for (ki=0; ki<ik; ki++)
    knots.push_back(a2_);

  // Define dummy coefficient vector
  vector<double> coefs(in*dim_, 0.0);

  // Make initial curves
  curr_crv_.resize(parvals_.size());
  for (size_t ki=0; ki<parvals_.size(); ++ki)
    curr_crv_[ki] = shared_ptr<SplineCurve>
      ( new SplineCurve(in, ik, &knots[0], &coefs[0], dim_));
 
  // Check distribution of knots and parameter values. If
  // necessary increase the smoothing weight
  adjustSmoothWeight();

  return;
}


//***************************************************************************

void ApproxCrvToSeqs::adjustSmoothWeight()
   //--------------------------------------------------------------------------
   //
   //     Purpose : Check distribution of knots and parameter values. If
   //               necessary increase the smoothing weight
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  // Compute the average number of parameter values for each knot interval
  int ki, kj;
  int minnmb = 1000*(int)parvals_[0].size();  // A large number
  double avnmb = 0.0;
  for (size_t kr=0; kr<parvals_.size(); ++kr)
    {
      int kn = curr_crv_[kr]->numCoefs();
      int kk = curr_crv_[kr]->order();
      vector<double>::const_iterator st = curr_crv_[kr]->basis().begin();
      int knpar = (int)parvals_[kr].size();
      int currnmb=0;
      avnmb += (double)knpar/(double)(kn-kk+1);
      for (kj=0, ki=kk; ki<kn; ki++)
	{
	  currnmb = 0;
	  while (kj<knpar && parvals_[kr][kj] >= st[ki-1] && 
		 parvals_[kr][kj] < st[ki])
	    {
	      currnmb++;
	      kj++;
	    }
	  minnmb = std::min(minnmb, currnmb);
	}
    }

  avnmb /= (double)(parvals_.size());
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

void ApproxCrvToSeqs::makeSmoothCurve(int idx)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Adapt the initial curve to approximate the given
   //               data points.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  SmoothCurve crvgen(dim_);

  int bd_fix[2];
  bd_fix[0] = bd_fix[1] = 1;
  double wgt1 = 0.000001*smoothweight_;
  double wgt2 = smoothweight_/smoothfac_;
  double wgt3 = smoothweight_/(smoothfac_*smoothfac_);
  double approxweight = 1.0 - wgt1 - wgt2 - wgt3;
  vector<double> pt_weight(parvals_[idx].size(), 1.0);
  vector<int> coef_known(curr_crv_[idx]->numCoefs(), 0);
  vector<sideConstraint> all_constraints, pt_constraints(2), tangent_constraints(2);
  getConstraints(idx, pt_constraints, tangent_constraints);
  // We set start pt of approximated curve.
  if (start_pt_[idx].size() == 0) // || (start_pt_[0].get() == 0))
    coef_known[0] =  1; // We lock start pt when approximating.
  else
    all_constraints.push_back(pt_constraints[0]);
  if (end_pt_[idx].size() == 0) // || (end_pt_[0].get() == 0))
    coef_known[coef_known.size()-1] =  1; // We lock end pt when approximating.
  else
    all_constraints.push_back(pt_constraints[1]);

  // We next add possible end tangents.
  if (start_pt_[idx].size() > 1) // && start_pt_[1].get() != 0)
      all_constraints.push_back(tangent_constraints[0]);    
  if (end_pt_[idx].size() > 1) // && end_pt_[1].get() != 0)
      all_constraints.push_back(tangent_constraints[1]);

//   crvgen.attach(curr_crv_, bd_fix);
  //crvgen.attach(curr_crv_, 0, &coef_known[0], all_constraints.size());
  crvgen.attach(curr_crv_[idx], &coef_known[0], (int)all_constraints.size());

  crvgen.setOptim(wgt1, wgt2, wgt3);

  crvgen.setLeastSquares(points_[idx], parvals_[idx], pt_weight, approxweight);

  if (all_constraints.size() != 0)
    crvgen.setSideConstraints(all_constraints);

  if (c1fac_ > 0.0)
    crvgen.setPeriodicity(1, c1fac_);

  // We then are ready to solve system.
  crvgen.equationSolve(curr_crv_[idx]);

  return;
}

//***************************************************************************

void ApproxCrvToSeqs::checkAccuracy(vector<double>& newknots)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Check the accuracy of the current curve compared to
   //               the given data points.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
//     // debug
//     std::ofstream debug("data/debug.g2");
//     SplineDebugUtils::writeSpaceParamCurve(*curr_crv_, debug);
//     // end of debug

    bool reparam = true;

//     double par_tol = 0.000000000001;
    maxdist_ = -10000.0;
    avdist_ = 0.0;

    // Traverse the data points and check the quality of the curve
    // approximations
    int kn = curr_crv_[0]->numCoefs();
    vector<int> refine(kn, 0);
    int nmb_allpnts = 0;
    for (size_t kr=0; kr<parvals_.size(); ++kr)
      {
	  int nmbpnt = (int)parvals_[kr].size();
	nmb_allpnts += nmbpnt;
	double dist;
	int left = 0, prevleft = 0;
	vector<double>::const_iterator st = curr_crv_[kr]->basis().begin();
	double ta = curr_crv_[kr]->startparam();
	double tb = curr_crv_[kr]->endparam();
	double par;
	int distOK = 1;
	double frac = 0.0;
	Point pos(dim_);
	Point clpos(dim_);
	double *pt, *param;
	int ki;
	for (ki=0, pt=&points_[kr][0], param=&parvals_[kr][0]; ki<nmbpnt; 
	     ki++, pt+=dim_, param++) {
	  // Compute closest point
	  pos.setValue(pt);
	  try {
	    curr_crv_[kr]->closestPoint(pos, ta, tb, 
					par, clpos, dist, param);
	  } catch (...) {
	    MESSAGE("Failed finding closest point.");
	    par = *param; // We're using our input value.
	    // This should at least be an upper boundary of closest dist.
	    dist = (pos - curr_crv_[kr]->ParamCurve::point(par)).length();
	  }
	  left = curr_crv_[kr]->basis().knotInterval(par);
	  if (reparam)
	    parvals_[kr][ki] = par;
	    

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
		refine[prevleft] = 1;
	      }

	      // Reset parameters
	      prevleft = left;
	      frac = 0.0;
	    }

	    // Start computing the next new knot
	    frac += dist;
	  }
	}
      

	if (!distOK && frac > 0) {
	  // Add the last new knot

	  refine[prevleft] = 1;
	}
      }
    avdist_ /= (double)nmb_allpnts;
	
    // Set the actual knots
    vector<double>::const_iterator st0 = curr_crv_[0]->basis().begin();
    for (size_t kr=0; kr<refine.size(); ++kr)
      {
	if (refine[kr])
	  {
	    double ta = st0[kr];
	    double tb = st0[kr+1];
	    newknots.push_back(0.5*(ta+tb));
	  }
      }
    return;
}

//***************************************************************************

int ApproxCrvToSeqs::doApprox(int max_iter)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Generate the approximating curves
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
   //--------------------------------------------------------------------------
{
  if (getenv("DEBUG"))
    {
      std::cout << "crv" << std::endl;
    }

  // Make the initial approximation. First fix endpoints.
  for (size_t kr=0; kr<parvals_.size(); ++kr)
    {
	int nmbpoints = (int)points_[kr].size()/dim_;
      vector<double>::iterator coef = curr_crv_[kr]->coefs_begin();
      int in = curr_crv_[kr]->numCoefs();
      int ki;
      for (ki=0; ki<dim_; ki++)
	{
	  coef[ki] = points_[kr][ki];
	  coef[(in-1)*dim_+ki] = 
	    points_[kr][(nmbpoints-1)*dim_+ki];
	}

      // Approximate
      makeSmoothCurve((int)kr);
    }

  if (max_iter == 0.0) { // In order to set maxerr_ & meanerr_.
      vector<double> newknots_dummy;
      checkAccuracy(newknots_dummy); 
  }

  // Check the accuracy and iterate the approximation
  double prevmax = 100000.0, prevav = 100000.0;
  int ki;
  for (ki=0; ki<max_iter; ki++)
    {
      vector<double> newknots;

      checkAccuracy(newknots); 

      if (getenv("DEBUG"))
	{
	  std::cout << "Iteration" << ki << " max "
		    << maxdist_ << " average " << avdist_ <<  std::endl;
	}

      if (maxdist_ <= aepsge_ || newknots.size() == 0)
	break;   // The required accuracy is reached.

      if (maxdist_ > 0.95*prevmax) {
	  MESSAGE("Convergence not too promising...");
// 	break;   // Not enough gain in refining
      }

      prevmax = maxdist_;
      prevav = avdist_;

      for (size_t kr=0; kr<curr_crv_.size(); ++kr)
	{
	  // Refine the spline spaces
	  curr_crv_[kr]->insertKnot(newknots);
 
	  // Approximate the data by a curves in the new spline space.
	  makeSmoothCurve((int)kr);
	}
    }

  return (ki == max_iter);
}

//***************************************************************************

void ApproxCrvToSeqs::getConstraints(int idx,
				     vector<sideConstraint>& pt_constraints,
				     vector<sideConstraint>& tangent_constraints)
{
  ALWAYS_ERROR_IF(pt_constraints.size() != 2 || tangent_constraints.size() != 2,
	      "Input vectors are not of size 2.");

  // We start by setting pt_constraints.
  if (start_pt_[idx].size() != 0) { // && (start_pt_[0].get() != 0)) {
    sideConstraint constraint;
    constraint.dim_ = curr_crv_[idx]->dimension();
    constraint.factor_.push_back(std::make_pair(0, 1.0));
    for (int ki = 0; ki < 3; ++ki)
      constraint.constant_term_[ki] = (start_pt_[idx][0])[ki];
    pt_constraints[0] = constraint;
  }
  if (end_pt_[idx].size() != 0) { // && (end_pt_[0].get() != 0)) {
    sideConstraint constraint;
    constraint.dim_ = curr_crv_[idx]->dimension();
    constraint.factor_.push_back(std::make_pair(curr_crv_[idx]->numCoefs() - 1, 1.0));
    for (int ki = 0; ki < 3; ++ki)
      constraint.constant_term_[ki] = (end_pt_[idx][0])[ki];
    pt_constraints[1] = constraint;
  }

  // We then set tangent_constraints.
  if (start_pt_[idx].size() > 1) { //&& (start_pt_[1].get() != 0)) {
    sideConstraint constraint;
    constraint.dim_ = curr_crv_[idx]->dimension();
    constraint.factor_.push_back(std::make_pair(1, 1.0));
    constraint.factor_.push_back(std::make_pair(0, -1.0));
    double t_inc = curr_crv_[idx]->basis().grevilleParameter(1) -
      curr_crv_[idx]->basis().grevilleParameter(0);
    for (int ki = 0; ki < 3; ++ki)
      constraint.constant_term_[ki] = ((start_pt_[idx][1])[ki])*t_inc;
    tangent_constraints[0] = constraint;
  }
  if (end_pt_[idx].size() > 1) { // && (end_pt_[1].get() != 0)) {
    sideConstraint constraint;
    constraint.dim_ = curr_crv_[idx]->dimension();
    constraint.factor_.push_back(std::make_pair(curr_crv_[idx]->numCoefs() - 1, 1.0));
    constraint.factor_.push_back(std::make_pair(curr_crv_[idx]->numCoefs() - 2, -1.0));
    double t_inc = 
      curr_crv_[idx]->basis().grevilleParameter(curr_crv_[idx]->numCoefs() - 1) -
      curr_crv_[idx]->basis().grevilleParameter(curr_crv_[idx]->numCoefs() - 2);
    for (int ki = 0; ki < 3; ++ki)
      constraint.constant_term_[ki] = ((end_pt_[idx][1])[ki])*t_inc;
    tangent_constraints[1] = constraint;
  }
}

//***************************************************************************
void ApproxCrvToSeqs::setEndPoints(const vector<vector<Point> >& start_pt,
				   const vector<vector<Point> >& end_pt)
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

vector<shared_ptr<SplineCurve> > 
ApproxCrvToSeqs::getApproxCurves(double& maxdist, 
				 double& avdist,
				 int max_iter)
   //--------------------------------------------------------------------------
   //
   //     Purpose : Return the approximating curve and parameters
   //               describing the accuracy.
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF,  10-12
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

void  ApproxCrvToSeqs::unsetSmooth()
{
  smoothweight_ = 0.0;
}

//***************************************************************************
void ApproxCrvToSeqs::setSmooth(double w)
{
   if (w>0.0 && w<1.0)
      smoothweight_=w;
}

//***************************************************************************

void  ApproxCrvToSeqs::setC1Approx(double fac)
{
  c1fac_ = fac;
}
