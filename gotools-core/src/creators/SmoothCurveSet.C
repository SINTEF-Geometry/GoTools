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

#include "GoTools/creators/SmoothCurveSet.h"
#include "GoTools/creators/SolveCGCO.h"
#include "GoTools/creators/SolveBCG.h"
#include "GoTools/creators/Integrate.h"

using std::vector;
using namespace Go;

void
SmoothCurveSet::preparePeriodicity(int cvidx, int seem)
//--------------------------------------------------------------------------
//     Purpose : Set pointers between identical coefficients at a periodic 
//               seem. If possible, update fixed coefficients at the seem.
//
//     Calls   : 
//
//     Written by : Sverre,  SINTEF Oslo,  06.2006.
//--------------------------------------------------------------------------
{
  int ki;
  int kn = cvs_[cvidx]->numCoefs();
  int idim1 = idim_; // Dimension of projective space.

  // Treat the periodicity.
  if (seem > 0)// && (coefknown_[cvidx][ki] > 0 && coefknown_[cvidx][kn-1] > 0))
    {
      if (coefknown_[cvidx][0] == 1)
	{
	  coefknown_[cvidx][kn-1] = 1;
	  for (ki=0; ki<idim1; ki++)
	    coef_array_[cvidx][(kn-1)*idim1+ki] =
	      coef_array_[cvidx][ki];
	}
      else if (coefknown_[cvidx][kn-1] == 1)
	{
	  coefknown_[cvidx][0] = 1;
	  for (ki=0; ki<idim1; ki++)
	    coef_array_[cvidx][ki] =
	      coef_array_[cvidx][(kn-1)*idim1+ki];
	}
      else
	coefknown_[cvidx][kn-1] = kpointer_;
    }
}


int SmoothCurveSet::setPeriodicity()
//--------------------------------------------------------------------------
//     Purpose : Set periodicity constraints
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF Oslo,  06.2006.
//--------------------------------------------------------------------------
{
  vector<shared_ptr<cvSetConstraint> > cv_set_constraints;
  int cvidx;
  for (cvidx = 0; cvidx < (int)cvs_.size(); ++cvidx)
    {
      if (cont_seem_[cvidx] > 1)
	{
	  for (int der = 1; der < cont_seem_[cvidx]; ++der)
	    {
	      shared_ptr<cvSetConstraint> cv_set_con =
		shared_ptr<cvSetConstraint>
		(new cvSetConstraint(cvidx,
				     cvs_[cvidx]->startparam(),
				     der,
				     cvidx,
				     cvs_[cvidx]->endparam(),
				     der,
				     false));
	      cv_set_constraints.push_back(cv_set_con);
	    }
	}
    }
  if (cv_set_constraints.size() > 0)
    {
      double dummy_wgt = 0.0;
      int kstat =
	setCvSetConstraints(cv_set_constraints,
			    false, dummy_wgt);
      if (kstat != 0)
	return kstat;
    }

  return 0;
}


int
SmoothCurveSet::set_weights(shared_ptr<SplineCurve> cv, double support_mult,
			      double weights[], double new_weights[])
//--------------------------------------------------------------------------
//  Purpose : Translate user weights into internal weights for smoothing.
//            Making sure that knot multiplicity is not in conflict with weights.
//            A 0.0 weight will remain 0.0, but other weights may alter.
//
//  Input   : cv       : The surface to be smoothed.
//            weights  : Weights representing the emphasis on various
//                       smoothing terms. The weights have the following
//                       understanding:
//                       weights[0] : Weight wrt 1st derivative, inside [0,0, 1.0]
//                       weights[1] : Weight wrt 2nd derivative, inside [0,0, 1.0]
//                       weights[2] : Weight wrt 3rd derivative, inside [0,0, 1.0]
//                       weights[3] : Approximation of orig cv.
//                       weights[0]+weights[1]+weights[2]+weights[3]=1.0
//
//  Output  : new_weights
//            return value : status 
//                                > 0 : Warning
//                                = 0 : OK
//                                < 0 : Error
//  Calls   : 
//
//  Written by : Sverre Briseid, SINTEF, May 2006 (based on version in SmoothSurf.C)
//
//--------------------------------------------------------------------------
{
  double eq_tol = 1.0e-12;
  int ki;
  int derivs_to_min;
  // Check continuity of curve.
  int inc;
  // In the rest of the routine we assume that weights sum to 1.0.
  // Possibly add as routine requirement?
  double weight_sum = weights[0] + weights[1] + weights[2] + weights[3];
  if (weight_sum != 1.0)
    {
      for (ki = 0; ki < 4; ++ki)
	weights[ki] /= weight_sum;
    }
  spline_space_cont(cv, inc);
  derivs_to_min = inc+1;
  if (inc < 1) // Only point continuous.
    {
      new_weights[0] = 1.0;
      new_weights[1] = new_weights[2] = 0.0;
      new_weights[3] = weights[3];
      return 2;       // Spline space not sufficiently continuous.
    }

  // Smoothing weights.
  double tlength;
  double fac;
  double omega_s = 1.0 - weights[3];
  double omega_a1 = 1.0 - omega_s;
  new_weights[0] = 0.1*weights[0];
  fac = (fabs(weights[1]+weights[2])<eq_tol) ? 0.5 : 
      weights[1]/(weights[1]+weights[2]);
  new_weights[1] = weights[1] + fac*0.9*weights[0];
  new_weights[2] = weights[2] + (1.0-fac)*0.9*weights[0];
  new_weights[3] = weights[3];

  // Adjust weight on 3. derivate according to the continuity of the
  // spline space. 
  if (derivs_to_min < 3)
    new_weights[2] = 0.0;
  // We reduse the weight on the second derivative
  if (derivs_to_min < 2)
    new_weights[1] /= 2.0;

  new_weights[2] = pow(new_weights[2], 0.1);
  new_weights[0] = (1.0 - new_weights[2])*(0.32*weights[0])*(0.32*weights[0]);
  new_weights[1] = new_weights[1]*new_weights[1]*new_weights[1]*(1.0-new_weights[2]);

  for (tlength=0.0, ki=0; ki<3; ki++)
    tlength += new_weights[ki];

  new_weights[0] *= omega_s*((1.0 - 0.9*weights[3])/tlength);
  new_weights[1] *= omega_s*((1.0 - 0.9*pow(weights[3], 0.08))/tlength);
  new_weights[2] *= omega_s*((1.0 - 0.9*weights[3]*sqrt(weights[3]))/tlength);

  new_weights[3] *= omega_a1;
   
  // Adjust for standard size of cv.
  double par_domain[2];
  par_domain[0] = cv->startparam();
  par_domain[1] = cv->endparam();
  double std_cv_size = 100.0;
  fac = (par_domain[1]-par_domain[0])/std_cv_size;

  if (derivs_to_min < 2)
    new_weights[1] = 0.0;
/*   new_weights[1] *= (fac*fac); */
/*   new_weights[2] *= (fac*fac*fac*fac); */
  new_weights[1] *= fac;
  new_weights[2] *= (fac*fac);

  // We also check the support of all basis functions.
  int min_deriv = get_min_deriv(cv, support_mult);
  double sum = 0.0;
  for (ki = 0; ki < 3; ++ki)
    if (ki < min_deriv)
      sum += new_weights[ki];
    else
      new_weights[ki] = 0.0;
  if (sum == 0.0) // We want at least smoothing wrt first derivs.
    new_weights[0] = 1.0;

  // Normalize weights.
  for (tlength=0.0, ki=0; ki<=3; ki++)
    tlength += new_weights[ki];
  for (ki=0; ki<=3; ki++) 
    new_weights[ki] /= tlength;

  return 0;
}


void
SmoothCurveSet::spline_space_cont(shared_ptr<SplineCurve> cv, int& nmbc)
//--------------------------------------------------------------------------
//  Purpose : Count the continuity of the spline space corresponding to
//            a given surface.
//
//  Input   : cv    - SurfacCurve represented in one spline space in each
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

  int in = cv->numCoefs();
  vector<double>::const_iterator st = cv->knotsBegin();
  for (ki=cv->order(), kmax=0; ki<in; ki=kj)
    {
      for (kj=ki+1, kcurr=1; kj<in && st[ki]==st[kj]; kj++, kcurr++);
      kmax = std::max(kmax, kcurr);
    }
  nmbc = cv->order() - kmax - 1;
}


//===========================================================================

int SmoothCurveSet::get_min_deriv(shared_ptr<SplineCurve> cv, double support_mult)
//--------------------------------------------------------------------------
//  Purpose : Find the minimum derivative in cv (counting mult of inner
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
  int min_deriv=0, nmbc;
  spline_space_cont(cv, nmbc);
  min_deriv = std::min(nmbc +1, 3);

  double dist1, dist2;
  int in=cv->numCoefs();
  
  double min_support_frac = std::numeric_limits<double>::max();
  double tol = -1.0;
  int ki;
  while (true)
    {
      min_support_frac = std::numeric_limits<double>::max(); // We reset the value for each iteration.
      int k1=cv->order() - min_deriv;
      for(ki=min_deriv; ki<in-1; ki++)
	{
	  dist1=(cv->knotsBegin()[ki+k1]-cv->knotsBegin()[ki]);
	  dist2=(cv->knotsBegin()[ki+k1+1]-cv->knotsBegin()[ki+1]);
	  double local_frac = dist2/dist1;
	  if (local_frac > 1.0)
	    local_frac = 1/local_frac;
	  if (local_frac < min_support_frac)
	    {
	      min_support_frac = local_frac;
	    }
	}
      tol=pow(1e-08,1.0/min_deriv); //(1e-10,1.0/min_deriv);
				    ////(1e-8,1.0/min_deriv);
      tol *= support_mult;
      if (min_support_frac<tol && min_deriv>1)
	{
	  min_deriv--;
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

  return min_deriv;
}


SmoothCurveSet::SmoothCurveSet()
  : idim_(3), kdim_(1), ider_(3), copyCoef_(1), kncond_(0), knconstraint_(0),
    kpointer_(3)
   //--------------------------------------------------------------------------
   //     Constructor for class SISmoothCurve.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Sverre Briseid,  SINTEF,  Nov 2005
   //--------------------------------------------------------------------------
{

}


SmoothCurveSet::~SmoothCurveSet()
{

}

int SmoothCurveSet::attach(vector<shared_ptr<SplineCurve> >& incvs,
			     vector<int>& seem,
			   vector<vector<int> >& coef_known,
			     int numSideConstraints)
{
   //--------------------------------------------------------------------------
   //
   //     Purpose : Initialize class variables related to the spline space.
   //
   //     Calls   :
   //
   //     Written by : Sverre Briseid,  SINTEF,  Nov 2005
   //     Revised by :
   //--------------------------------------------------------------------------
  if (incvs.size() < 1)
    {
      return -1;
    }
  int kstat = 0;
  int kh, ki, kj;
  idim_ = incvs[0]->dimension();
  int nmbcvs = (int)incvs.size();

  // Prepare for data storage for each curve.
  if ((int)cvs_.size() != nmbcvs)
    {
      // We must update data structure to fit the new size of the curves.
      cvs_.resize(nmbcvs);
      coef_array_.resize(nmbcvs);
      pivot_.resize(nmbcvs);
      cv_integral_.resize(nmbcvs);
      for (kh=0; kh<nmbcvs; kh++)
	cv_integral_[kh] = shared_ptr<integralInfo>(new integralInfo());
    }
  coefknown_ = coef_known;
  cont_seem_ = seem;

  kncond_ = 0;  // No free coefficients found yet.
  for (kh=0; kh<nmbcvs; kh++)
    {
      // Fetch data from B-spline curves.
      bool israt = incvs[kh]->rational();
      int rdim = (israt) ? idim_+1 : idim_;
      int kn = incvs[kh]->numCoefs();

      vector<double>::iterator scoef_begin = (israt) ? incvs[kh]->rcoefs_begin() :
	incvs[kh]->coefs_begin();
      vector<double>::iterator scoef_end = (israt) ? incvs[kh]->rcoefs_end() :
	incvs[kh]->coefs_end();

      if (copyCoef_)
	{
	  coef_array_[kh].resize(kn*rdim);
	  std::copy(scoef_begin, scoef_end, coef_array_[kh].begin());
	}

      if (israt)
	{
	  // Divide by the weight.
	  vector<double>::iterator sc = scoef_begin;
	  for (ki=0; ki<kn; ki++, sc+=rdim)
	    for (kj=0; kj<idim_; kj++)
	      sc[kj] /= sc[idim_];
	}

      if (cvs_[kh].get() != NULL &&  
	  incvs[kh]->numCoefs() == cvs_[kh]->numCoefs() && 
	  incvs[kh]->order() == cvs_[kh]->order() && 
	  numSideConstraints == knconstraint_)
	{
	  // It is assumed that the new surface (incvs[kh]) lies in the same
	  // spline space that the old surface and that there
	  // is no change concerning the existance of normal conditions.
	  // If these conditions are broken, memory faults and/or wrong
	  // results will occur. Including higher order derivatives in the
	  // smoothing functional than in previous runs, will have no effect.

	}
      else
	{
	  cv_integral_[kh]->resize(ider_, kn);
	  cvs_[kh] = incvs[kh];

	  // reset pivot_ array.
	  pivot_[kh].resize(kn);
	  std::fill(pivot_[kh].begin(), pivot_[kh].end(), 0);
	}
      cvs_[kh] = incvs[kh];

//        // Store information about free, fixed and irrelevant coefficients
//        // and create pivot_ array. Count number of free coefficients.
//        coefknown_[kh] = coef_known[kh].begin(); //.begin();

      // We must not forget to set seem for periodic cvs_.
      if (seem[kh] > 0)
	{
	  preparePeriodicity(kh, seem[kh]);
	}

      for (ki=0; ki<kn; ki++)
	{
	  if (coefknown_[kh][ki] == 0)
	    {
	      pivot_[kh][ki] = kncond_;
	      kncond_++;
	    }
	  else
	    {
	      pivot_[kh][ki] = -MAXINT;
	    }
	}
    }

  // Add side constraints
  knconstraint_ = numSideConstraints;
  kncond_ += knconstraint_;   // @@@ VSK. What about constraints with no 
  // free coefficients?

  // Allocate scratch for arrays in the equation system. 

  gmat_.resize(kdim_*kdim_*kncond_*kncond_);
  std::fill(gmat_.begin(), gmat_.end(), 0.0);
  gright_.resize(idim_*kncond_);
  std::fill(gright_.begin(), gright_.end(), 0.0);

  // Now that the equation system is set up we may add the seem constraints.
  kstat = setPeriodicity();
  if (kstat < 0)
    return kstat;

  return 0;
}

int SmoothCurveSet::setOptimize(double weight1,
				double weight2,
				double weight3)
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		    from the smoothness functional.
//
//
//     Written by : Sverre Briseid,  SINTEF,  Nov 2005
//--------------------------------------------------------------------------
{
  int grstat = 0;   // Initialize status variable.
  int ki, kp, kr;
  // Number of coefficients to compute,
  // i.e. size of matrix of equation
  // system. 
  int kl1, kl2;
  double tval;      // Contribution to equation system from smoothness term.
  vector<double>::iterator  scoef; // Pointer to curve coefficients. 


  // We start by updating the weights according to spline spaces.
  vector<double> orig_wgts(4, 0.0);
  orig_wgts[0] = weight1;
  orig_wgts[1] = weight2;
  orig_wgts[2] = weight3;
  vector<double> new_wgts(orig_wgts);
  grstat = setWeights(&orig_wgts[0], &new_wgts[0]);
  if (grstat != 0)
    return grstat;

  // Traverse all B-splines and set up matrices of equation system.
  int idxcv;
  for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
    {
      int inc=1;
      spline_space_cont(cvs_[idxcv], inc);
      if (inc < 1) // Only point continuous.
	{
	  return 2;       // Spline space not sufficiently continuous.
	}
      // Set parameter area.
      int kn = cvs_[idxcv]->numCoefs();
      double ta = cvs_[idxcv]->startparam();   // Start of par. interval.
      double tb = cvs_[idxcv]->endparam();     // End of par. interval.
      vector<double>::iterator sc;  // Pointer into coefficient array of the original curve.
      bool israt = cvs_[idxcv]->rational();
      if (copyCoef_)
	scoef = coef_array_[idxcv].begin();
      else
	scoef = (israt) ? cvs_[idxcv]->rcoefs_begin()
	  : cvs_[idxcv]->coefs_begin();
      shared_ptr<integralInfo> cig = cv_integral_[idxcv];
       if (cig->integralset_ == 0)
	 {
	   GaussQuadInner(cvs_[idxcv]->basis(), ider_, ta, tb, cig->integral_);
	   cig->integralset_ = 1;
	 }

      for (ki=0; ki<kn; ki++)
	for (kp=0; kp<kn; kp++)
	  {
	    // We must locate the index of the coef in the equation system.
	    if ((kp > ki) && (coefknown_[idxcv][ki] == 0))
	      continue;
 	    if (coefknown_[idxcv][kp] == 1 || coefknown_[idxcv][kp] == 2)
 	      continue;
	    kl1 = (coefknown_[idxcv][ki] > 2) ?
	      pivot_[idxcv][coefknown_[idxcv][ki]-kpointer_] : pivot_[idxcv][ki];
	    kl2 = (coefknown_[idxcv][kp] > 2) ?
	      pivot_[idxcv][coefknown_[idxcv][kp]-kpointer_] : pivot_[idxcv][kp];

// 	    if (kl2 > kl1)
// //  &&
// // 		!(ki < cont_bound[idxcv][0] || kn-ki-1 < cont_bound[idxcv][1]))
// 	      continue;

// 	    if (kp < cont_bound[idxcv][0] || kn-kp-1 < cont_bound[idxcv][1])
// 	      continue;

	    // Compute an element in the left side matrix.
	    // Complete term of the smoothness function.
	    tval = new_wgts[0]*cig->integral_[1][ki][kp] + 
	      new_wgts[1]*cig->integral_[2][ki][kp] +
	      new_wgts[2]*cig->integral_[3][ki][kp];

	    if (coefknown_[idxcv][ki] == 1)
// 	    if (ki < cont_bound[idxcv][0] || kn-ki-1 < cont_bound[idxcv][1])
	      {
		// The contribution of this term is added to the right
		// side of the equation system. First fetch the known
		// coefficient.

		sc = scoef + ki*idim_;
		for (kr=0; kr<idim_; kr++)
		  gright_[kr*kncond_+kl2] -= sc[kr]*tval;
	      }
	    else
	      {
		// The contribution of this term is added to the left
		//  side of the equation system.
		for(kr=0; kr<kdim_; kr++)
		  {
		    gmat_[(kr*kncond_+kl1)*kdim_*kncond_+kr*kncond_+kl2] +=tval;
		    if (kl2 < kl1)
		      gmat_[(kr*kncond_+kl2)*kdim_*kncond_+kr*kncond_+kl1] +=tval;
		  }
	      }
	  }
    }
   
  return 0;
}

int
SmoothCurveSet::setLeastSquares(const vector<vector<double> >& pnts,
				const vector<vector<double> >& param_pnts,
				const vector<vector<double> >& pnt_weights,
				double weight)
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of data points.
//
//
//     Written by : Sverre Briseid,  SINTEF,  Nov 2005
//--------------------------------------------------------------------------
{
   // int kstat = 0;
   int k1, k2, k3, k4, kr2;
   int kl1, kl2;
   int kj;
   int kleft;
   double tz;     // Help variable.  
   vector<double>::iterator sc;    // Pointer into the coefficient array of the original curve.
//    vector<double>::iterator  scoef; // Pointer to surface coefficients. 

   // Traverse all points in the pointset. 
   int idxcv;
   for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
     {
       int kk = cvs_[idxcv]->order();
       bool israt = cvs_[idxcv]->rational();
        vector<double>::iterator scoef;
       if (copyCoef_)
	 scoef = coef_array_[idxcv].begin();
       else
	 scoef = (israt) ? cvs_[idxcv]->rcoefs_begin() : cvs_[idxcv]->coefs_begin();
       // Allocate scratch for B-spline basis functions. 
       vector<double> sbasis;
       const double *pnt = &pnts[idxcv][0];
       int num_points = (int)pnts[idxcv].size()/idim_;
       for (int kr=0; kr<num_points; kr++, pnt+=idim_)
	 {
	   // Fetch B-spline basis functions different from zero. 
	   kleft = cvs_[idxcv]->basis().knotInterval(param_pnts[idxcv][kr]);
	   sbasis = cvs_[idxcv]->basis().computeBasisValues(param_pnts[idxcv][kr], 0);

	   for (k1=kleft-kk+1, k2=0; k1<=kleft; k1++, k2++)
	     {
	       // Test if the B-spline basis is too close to the boundary, i.e.
	       // the coefficient is set by the boundary conditions.           
	       if (coefknown_[idxcv][k1] == 1 || coefknown_[idxcv][k1] == 2)
		 continue;

	       kl1 = (coefknown_[idxcv][k1] > 2) ?
		 pivot_[idxcv][coefknown_[idxcv][k1]-kpointer_] :
		 pivot_[idxcv][k1];

	       tz = weight*pnt_weights[idxcv][kr]*sbasis[k2];

	       // Add contribution to right hand side.
	       for (kj=0; kj<idim_; kj++)
		 gright_[kj*kncond_+kl1] += tz*pnt[kj];

	       for (k3=kleft-kk+1, k4=0; k3<=kleft; k3++, k4++)
		 {
		   // Test if the B-spline basis is too close to the boundary.
		   if (coefknown_[idxcv][k3] == 1 || coefknown_[idxcv][k3] == 2)
		     {
		       // Adjust for boundary conditions. 
		       sc = scoef + k3*idim_;
		       for (kj=0; kj<idim_; kj++)
			 gright_[kj*kncond_+kl1] -= sc[kj]*tz*sbasis[k4];
		     }
		   else
		     {
		       // The term gives a contribution to the left hand side. 
		       kl2 = (coefknown_[idxcv][k3] > 2) ?
			 pivot_[idxcv][coefknown_[idxcv][k3]-kpointer_] :
			 pivot_[idxcv][k3];

 		       if (kl2 > kl1) continue;
		 
		       for(kr2=0; kr2<kdim_; kr2++)
			 {
			   gmat_[(kr2*kncond_+kl1)*kdim_*kncond_+kr2*kncond_+kl2]
			     += tz*sbasis[k4];
			   if (kl2 < kl1)
			     gmat_[(kr2*kncond_+kl2)*kdim_*kncond_+kr2*kncond_+kl1]
			       += tz*sbasis[k4];
			 }
}
		 }
	     }
	 }

     }
   return 0;
}

void SmoothCurveSet::setApproxOrig(double weight)
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of an original curve
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF Oslo,  Nov 2005.
//--------------------------------------------------------------------------
{
  int ki, kp, kr;
  int kl1, kl2;
  int kstart, kend;
  double tval;
  std::vector<double>::iterator sc;    // Pointer into the coefficient array of the original curve.

  // Travers all B-splines and set up matrices of equation system.
  int idxcv;
  for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
    {
      int kn = cvs_[idxcv]->numCoefs();
      int kk = cvs_[idxcv]->order();
      bool israt = cvs_[idxcv]->rational();
      std::vector<double>::iterator scoef;
       if (copyCoef_)
	 scoef = coef_array_[idxcv].begin();
       else
	 scoef = (israt) ? cvs_[idxcv]->rcoefs_begin() : cvs_[idxcv]->coefs_begin();
      for (kp=0; kp<kn; kp++)
	{
	  if (coefknown_[idxcv][kp] == 1 || coefknown_[idxcv][kp] == 2)
	    continue;
	  kstart = std::max(0,kp-kk+1);
	  kend = std::min(kp+kk, kn);
	  for (ki=kstart; ki<kend; ki++)
	    {
	      kl1 = (coefknown_[idxcv][ki] > 2) ? 
		pivot_[idxcv][coefknown_[idxcv][ki]-kpointer_] :
		pivot_[idxcv][ki];
	      kl2 = (coefknown_[idxcv][kp] > 2) ? 
		pivot_[idxcv][coefknown_[idxcv][kp]-kpointer_] :
		pivot_[idxcv][kp];

	      if (coefknown_[idxcv][ki] == 1 || coefknown_[idxcv][ki] == 2)
		continue;

	      // Compute an element in the left side matrix.
	      tval = weight*cv_integral_[idxcv]->integral_[0][ki][kp];

	      // The contribution of this term is added to the right
	      // side of the equation system. First fetch the known
	      // coefficient.

	      sc = scoef + ki*idim_;
	      for (kr=0; kr<idim_; kr++)
		gright_[kr*kncond_+kl2] += sc[kr]*tval;

	      // The contribution of this term is added to the left
	      //  side of the equation system.
	      for(kr=0; kr<kdim_; kr++)
		{
		  gmat_[(kr*kncond_+kl1)*kdim_*kncond_+kr*kncond_+kl2] +=tval;
		}

	    }
	}
    }
   
  return;
}



int
SmoothCurveSet::setOrthCond(const vector<vector<double> >& pnts,
			    const vector<vector<double> >& param_pnts,
			    double weight)
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of normal directions.
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF,  Nov 2005
//--------------------------------------------------------------------------
{
  if (kdim_ != 3)
    return 1;   // Not prepared for approximation of normals.

  //int nmbpoint = num_pnts;   // Number of normal directions
  // int kstat = 0;
  int kr1,kr2;
  int k1, k2, k3, k4;
  int kl1, kl2;
  int kleft=0;  // Parameter used to position in knot vector
  double tz1, tz2, tz3;         // Help variables.  
  vector<double>::iterator sc;  // Pointer into the coefficient array of the original surf.

  // Traverse all points in the pointset. 
  int idxcv;
  for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
    {
      int kk = cvs_[idxcv]->order();
      // int kn = cvs_[idxcv]->numCoefs();
      vector<double>::iterator st = cvs_[idxcv]->knotsBegin();
      bool israt = cvs_[idxcv]->rational();
       vector<double>::iterator scoef;
       if (copyCoef_)
	 scoef = coef_array_[idxcv].begin();
       else
	 scoef = (israt) ? cvs_[idxcv]->rcoefs_begin() : cvs_[idxcv]->coefs_begin();
      vector<double> sbasis;
      const double *pnt = &pnts[idxcv][0];
      int num_pnts = (int)pnts[idxcv].size()/idim_;
      for (int kp=0; kp<num_pnts; kp++, pnt+=idim_)
	{
	  // Fetch B-spline basis functions different from zero. 
	  kleft = cvs_[idxcv]->basis().knotInterval(param_pnts[idxcv][kp]);
	  sbasis = cvs_[idxcv]->basis().computeBasisValues(param_pnts[idxcv][kp], 0);

	  for (k1=kleft-kk+1, k2=0; k1<=kleft; k1++, k2++)
	    {
	      // Test if the B-spline basis is too close to the boundary, i.e.
	      // the coefficient is set by the boundary conditions.           
	      if (coefknown_[idxcv][k1] == 1 || coefknown_[idxcv][k1] == 2)
		continue;

	      kl1 = (coefknown_[idxcv][k1] > 2) ?
		pivot_[idxcv][coefknown_[idxcv][k1]-kpointer_] :
		pivot_[idxcv][k1];

	      tz1 = weight*sbasis[k2];

	      for (k3=kleft-kk+1, k4=0; k3<=kleft; k3++, k4++)
		{
		  // Test if the B-spline basis is too close to the boundary.
		  tz2=tz1*sbasis[k4];
		  if (coefknown_[idxcv][k3] == 1 || coefknown_[idxcv][k3] == 2)
		    {
		      // Adjust for boundary conditions. 
		      sc = scoef + k3*idim_;
		      for(kr1=0; kr1<idim_; kr1++)
			{
			  tz3=tz2*pnt[kr1];
			  for(kr2=0; kr2<idim_; kr2++)
			    gright_[kr1*kncond_+kl1]-=tz3*pnt[kr2]*sc[kr2];
			}
		    }
		  else
		    {
		      // The term gives a contribution on the left hand side. 
		 
		      kl2 = (coefknown_[idxcv][k3] > 2) ?
			pivot_[idxcv][coefknown_[idxcv][k3]-kpointer_] :
			pivot_[idxcv][k3];
		      for(kr1=0; kr1<kdim_; kr1++)
			{
			  tz3=tz2*pnt[kr1];
			  for(kr2=0; kr2<idim_; kr2++)
			    {
			      gmat_[(kr1*kncond_+kl1)*kncond_*kdim_+
				   kr2*kncond_+kl2]+=
				tz3*pnt[kr2];
			    }
			}
		    }
		}
	    }
	}
    }

  return 0;
}


int
SmoothCurveSet::updateSideConstraints(vector<shared_ptr<sideConstraintSet> >&
					constraints,
				      const vector<vector<int> >& coef_known)
//--------------------------------------------------------------------------
//     Purpose : Insert values for known coefs, remove constraints when all
//               are known
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF Oslo,  06.02.
{
  int i, j, k;
  double equality_eps = 1e-06;
  for (i = 0; i < (int)constraints.size(); ++i)
    {
      int dim = constraints[i]->dim_; // @@ Handling rational sfs?
      for (j = 0; j < (int)constraints[i]->factor_.size(); ++j)
	{
 	  int cv_ind = constraints[i]->factor_[j].first.first;
	  if (dim != idim_) {
	    return -106;
	  }

	  int coef_ind = constraints[i]->factor_[j].first.second;
	  bool israt = cvs_[cv_ind]->rational();
	  vector<double>::iterator scoef;
	  if (copyCoef_)
	    scoef = coef_array_[cv_ind].begin();
	  else
	    scoef = (israt) ? cvs_[cv_ind]->rcoefs_begin() : cvs_[cv_ind]->coefs_begin();
	  if (coef_known[cv_ind][coef_ind] != 0) {
	    // As coef was already known, we add to const term on right side.
// 	    int kl1 = (coefknown_[cv_ind][coef_ind] > 2) ?
// 	      pivot_[cv_ind][coefknown_[cv_ind][coef_ind]-kpointer_] :
// 	      pivot_[cv_ind][coef_ind];
	    for (k = 0; k < dim; ++k)
	      {
		constraints[i]->constant_term_[k] -=
		  scoef[dim*coef_ind+k]*constraints[i]->
		  factor_[j].second;
	      }
	    // We must then update other pointers after removing object.
	    constraints[i]->factor_.erase(constraints[i]->factor_.begin()+j);
	    --j;
	  }
	}
      // If all coefs in constraint were known, we remove constraint.
      if (constraints[i]->factor_.size() == 0)
	{
	  double sum = 0.0;
	  for (j = 0; j < dim; ++j)
	    sum += fabs(constraints[i]->constant_term_[j]);
	  if (sum > equality_eps*10000.0)
	    { // Check whether constraint was fulfilled.
// 	      puts("distance too large!"); // @@sbr
// 	      puts("0.0 != " && sum && "!");
	    }
	  constraints.erase(constraints.begin() + i); // We remove constraint.
	  --i;
	}
    }

  return 0;
}


void
SmoothCurveSet::setSideConstraints(vector<shared_ptr<sideConstraintSet> >&
				     constraints,
				     bool replace_constraints)
//--------------------------------------------------------------------------
//     Purpose : Set linear side constraints to the minimization problem,
//               solve using the method of Lagrange multipliers.
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF,  05.02
//--------------------------------------------------------------------------
{
  // Given an equation system based on least squares and a smoothing term:
  //
  //      ((A^T)A + wE)x = (A^T)b,
  //
  // we add linear side constraints
  //
  //      Cx = d,
  //
  // resulting in the matrix eqation
  //
  //      { (A^T)A + wE       C^T } {   x    }     { (A^T)b }
  //      {                       } {        }  =  {        }
  //      {      C             0  } { lambda }     {    d   },
  //
  // where lambda denotes the Lagrange Multipliers.
  // We thus have to add new matrices C, C^T & d to the system.
  updateSideConstraints(constraints, coefknown_);

  // Normal conditions (kdim_ = 3) should now work.
  // We check whether constraints are achievable. If coef in constraint is
  // already known, add to const term on right side. If all coefs are
  // removed, delete constraint. Issue warning if constraint is impossible.
  int nmb_constraints = (int)constraints.size();
  // If any side constraints were removed, we must update size of matrix.
  int ki;
  int keep_size = (replace_constraints) ?
    (kncond_ - knconstraint_) : kncond_; // The number of rows/columns to keep.
  if ((!replace_constraints) || (nmb_constraints != knconstraint_))
    {
      int new_kncond = (replace_constraints) ?
	kncond_ - knconstraint_ + nmb_constraints :
	kncond_ + nmb_constraints;
      // For ease of algorithm, we copy matrices to new matrices.

      vector<double> new_gmat(kdim_*kdim_*new_kncond*new_kncond, 0.0);
      vector<double> new_gright(idim_*new_kncond, 0.0);

      for (int kj = 0; kj < kncond_; ++kj)
	{ // @@sbr Verify that code is correct when not using sparse stuff.
	  std::copy(gmat_.begin() + kj*kdim_*kncond_,
		    gmat_.begin() + kj*kdim_*kncond_ + keep_size,
		    new_gmat.begin() + kj*new_kncond);
	}
      for (ki = 0; ki < idim_; ++ki)
	{
	  std::copy(gright_.begin() + ki*kncond_,
		    gright_.begin() + ki*kncond_ + keep_size,
		    new_gright.begin() + ki*new_kncond);
	}
      knconstraint_ = (replace_constraints) ? nmb_constraints :
	knconstraint_ + nmb_constraints;
      kncond_ = new_kncond;
      // We must release old values.

      gmat_ = new_gmat;
      gright_ = new_gright;
    }

  // @@ We're still assuming kdim_ == 1...
  // We update values in matries as described in the above equations.
  // Dimension of vector x (# coefs).
//   int nmb_free_coefs = kncond_ - knconstraint_;
  // We start by updating gmat_ by adding new elements given by side
  // constraints.
  int kpointer = 0;
  for (ki = 0; ki < nmb_constraints; ++ki)
    {
      for (int kj = 0; kj < (int)constraints[ki]->factor_.size(); ++kj)
	{ // We start with gmat_->
	  int cv_id = constraints[ki]->factor_[kj].first.first;
	  int coef_id = constraints[ki]->factor_[kj].first.second;
// 	    int kl1 = coef_id; // - cont_bound[0];
	  int piv_id = (coefknown_[cv_id][coef_id] > 2) ?
	    pivot_[cv_id][coefknown_[cv_id][coef_id]-kpointer] :
	    pivot_[cv_id][coef_id];
	  // We have made  sure that all elements in constraints[i] are free.

	  gmat_[kncond_*(keep_size+ki)+piv_id] =
	    constraints[ki]->factor_[kj].second;
	  gmat_[kncond_*piv_id+keep_size+ki] =
	    constraints[ki]->factor_[kj].second;
	}
    }

  // We next update gright_ by adding const values given by side constraints.
  for (ki = 0; ki < idim_; ++ki)
    {
      for (int kj = 0; kj < nmb_constraints; ++kj)
	{
	  gright_[ki*kncond_+keep_size+kj] =
	    constraints[kj]->constant_term_[ki];
	}
    }

  return;
}


int SmoothCurveSet::setApproxSideConstraints(vector<shared_ptr<sideConstraintSet> >&
					     constraints,
					     double weight)
//--------------------------------------------------------------------------
//     Purpose : Use least sqyares to approximate the constraints.
//
//     Calls   : 
//
//     Written by : Sverre Briseid, SINTEF Oslo, Dec 2005.
//--------------------------------------------------------------------------
{
//   puts("Under construction ... ");
//   return -1;

  updateSideConstraints(constraints, coefknown_);

  int kstat = 0;
  int ki, kj, kk;
  int nmb_constraints = (int)constraints.size();
  if (nmb_constraints == 0)
    {
      kstat = 1; // We issue a warning as there was no need to call routine.
      return kstat;
    }



  // As e would like to minimize ||Cx-d||, where x is the unknown coefs and
  // C & d are given by the constraints, we seek a solution to
  // Ax = (C^T)Cx = (C^T)d = e.
  // We add values to both sides of the system. We multiply both terms by
  // weight, allowing easy adjustment of factor.
  int dim = constraints[0]->dim_;
  for (ki = 0; ki < nmb_constraints; ++ki)
    {
      for (kj = 0; kj < (int)constraints[ki]->factor_.size(); ++kj)
	{
	  // We add elements given by constraints to left side of equation.
	  int kj_cv_id = constraints[ki]->factor_[kj].first.first;
	  int kj_coef_id = constraints[ki]->factor_[kj].first.second;
	  int kj_id = (coefknown_[kj_cv_id][kj_coef_id] > 2) ?
	    pivot_[kj_cv_id][coefknown_[kj_cv_id][kj_coef_id]-kpointer_] :
	    pivot_[kj_cv_id][kj_coef_id];
	  for (kk = 0; kk < (int)constraints[ki]->factor_.size(); ++kk)
	    {
	      int kk_cv_id = constraints[ki]->factor_[kk].first.first;
	      int kk_coef_id =
		constraints[ki]->factor_[kk].first.second;
	      int kk_id = (coefknown_[kk_cv_id][kk_coef_id] > 2) ?
		pivot_[kk_cv_id][coefknown_[kk_cv_id][kk_coef_id]-kpointer_] :
		pivot_[kk_cv_id][kk_coef_id];
	      double term = constraints[ki]->factor_[kj].second *
		constraints[ki]->factor_[kk].second;
	      gmat_[kj_id*kncond_ + kk_id] += term*weight;
	    }

	  // We then compute right side of equation.
	  for (kk = 0; kk < dim; ++kk)
	    {
	      gright_[kncond_*kk+kj_id] +=
		constraints[ki]->factor_[kj].second *
		constraints[ki]->constant_term_[kk]*weight;
	    }
	}
    }


  return kstat;
}


vector<std::pair<std::pair<int,int>, double > >
SmoothCurveSet::getSideConstraint(int cv_id,
				  double tpar,
				  int der,
				  int sign,
				  int *jstat)
//--------------------------------------------------------------------------
//     Purpose : Extract side constraint expression in der in input cv.
//
//     Calls   : 
//
//     Written by : Sverre Briseid, SINTEF Oslo, Nov 2005.
//--------------------------------------------------------------------------
{
  vector<std::pair<std::pair<int,int>, double > >  side_constraints;

  // We then extract the basis functions which are non-zero in par.
  // We first locate the interval of param.
  int kk = cvs_[cv_id]->order();
  int kleft = cvs_[cv_id]->basis().knotInterval(tpar);

  // We then compute the derivs.
  int kder;
  if (cvs_[cv_id]->rational())
    kder = std::min(kk-1, der);
  else
    kder = der;

  vector<double> ebder = 
    cvs_[cv_id]->basis().computeBasisValues(tpar, kder);

  // We then run through the bezier functions different from 0.0.
  int ctr = 0;
  int kj;
  for (kj=kleft-kk+1; kj<=kleft; kj++)
    {
      side_constraints.push_back(std::make_pair(std::make_pair(cv_id, kj),
					   sign*ebder[ctr*(kder+1)+kder]));
      ++ctr;
    }

  return side_constraints;
}


int
SmoothCurveSet::setWeights(double wgts[], double new_wgts[])
{

  double wgt_sum = wgts[0] + wgts[1] + wgts[2] + wgts[3];
  int idxcv, ki;
  int nmb_unknown_coefs = 0;
  for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
      for (ki = 0; ki < (int)coefknown_[idxcv].size(); ++ki)
      if (coefknown_[idxcv][ki] == 0)
	++nmb_unknown_coefs;
  for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
    {
//       int nmb_unknown_coefs = 0;
//       for (ki = 0; ki < coefknown_[idxcv].size(); ++ki)
// 	if (coefknown_[idxcv][ki] == 0)
// 	  ++nmb_unknown_coefs;
      double support_mult = 10.0*(double)nmb_unknown_coefs/1000.0;
//   if (support_mult > 1.0)
//     support_mult = 1.0;
      if (support_mult < 0.01)
	support_mult = 0.01;

      vector<double> cv_wgts(4, 0.0);
      int grstat = set_weights(cvs_[idxcv], support_mult,
			       &wgts[0], &cv_wgts[0]);
      if (grstat != 0)
	  return grstat;
      for (ki = 0; ki < 3; ++ki)
	if (cv_wgts[ki] < new_wgts[ki])
	  new_wgts[ki] = cv_wgts[ki];
    }
  // We rescale to original wgt_sum.
  double new_wgt_sum = new_wgts[0] + new_wgts[1] + new_wgts[2];
  for (ki = 0; ki < 3; ++ki)
    {
      new_wgts[ki] *= wgt_sum/new_wgt_sum;
    }

  return 0;
}


void
SmoothCurveSet::setInterpolationConditions(const vector<vector<double> >&
					     pnts,
					   const vector<vector<double> >&
					     param_pnts,
					   const vector<vector<int> >& der,
					     bool appr_constraints,
					     double appr_wgt,
					     int* jstat)
//--------------------------------------------------------------------------
//     Purpose : Set interpolation conditions using linear side constraints.
//
//     Calls   : 
//
//     Written by : Sverre Briseid, SINTEF Oslo, Nov 2005.
//--------------------------------------------------------------------------
{
  // We run through all pnts creating the linear side constraints.
  int kstat = 0;
  vector<shared_ptr<sideConstraintSet> > side_constraint_sets;

  *jstat = 0;

  int ki, kj;
  int idxcv;
  for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
    {
	int num_points = (int)pnts[idxcv].size()/idim_;
      for (ki = 0; ki < num_points; ++ki)
	{
	  shared_ptr<sideConstraintSet> curr_constraint = 
	    shared_ptr<sideConstraintSet>(new sideConstraintSet(idim_));
	  for (kj = 0; kj < idim_; ++kj)
	    {
	      curr_constraint->constant_term_[kj] = 0.0;
	    }

	  double tpar = param_pnts[idxcv][ki];
	  int sign = 1;
	  vector<std::pair<std::pair<int,int>, double > > side_constraints =
	    getSideConstraint(idxcv, tpar, der[idxcv][ki], sign, &kstat);
	  if (kstat < 0)
	    {
	      *jstat = kstat;
	      return;
	    }
	  curr_constraint->factor_.insert
	    (curr_constraint->factor_.end(),
	     side_constraints.begin(), side_constraints.end());

	  for (kj = 0; kj < idim_; ++kj)
	    {
	      curr_constraint->constant_term_[kj] += pnts[idxcv][ki*idim_+kj];
	    }
	  side_constraint_sets.push_back(curr_constraint);
	}
    }

  kstat = updateSideConstraints(side_constraint_sets, coefknown_); //@@sbr
  if (kstat < 0)
    {
      *jstat = kstat;
      return;
    }

  if (appr_constraints)
    {  
      kstat = setApproxSideConstraints(side_constraint_sets, appr_wgt);
      if (kstat < 0)
	{
	  *jstat = kstat;
	  return;
	}
    }
  else    
    {  
      setSideConstraints(side_constraint_sets, false);
    }

}


int SmoothCurveSet::setCvSetConstraints(const vector<shared_ptr<cvSetConstraint> >&
					  cv_set_constraints,
					  bool appr_constraints,
					  double appr_wgt)
//--------------------------------------------------------------------------
//     Purpose : Set constraints between derivatives in pairs of input cvs.
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF,  Nov 2005
//--------------------------------------------------------------------------
{
  // It would seem nicer to have each interpolation condition contained in an
  // object ...
  int kstat = 0;
  int ki, kj;

  vector<shared_ptr<sideConstraintSet> > side_constraint_sets;
  for (ki = 0; ki < (int)cv_set_constraints.size(); ++ki)
    {
      shared_ptr<sideConstraintSet> curr_constraint_set = 
	shared_ptr<sideConstraintSet>(new sideConstraintSet(idim_));
      for (kj = 0; kj < idim_; ++kj)
	{
	  curr_constraint_set->constant_term_[kj] = 0.0;
	}

      int sign = 1;
      vector<std::pair<std::pair<int,int>, double > > side_constraints1 =
	getSideConstraint(cv_set_constraints[ki]->cv1_id_,
			  cv_set_constraints[ki]->cv1_par_,
			  cv_set_constraints[ki]->cv1_der_,
			  sign, &kstat);
      if (kstat < 0)
	{
	  return kstat;
	}

      sign = (cv_set_constraints[ki]->opp_) ? 1 : -1;
      vector<std::pair<std::pair<int,int>, double > > side_constraints2 =
	getSideConstraint(cv_set_constraints[ki]->cv2_id_,
			  cv_set_constraints[ki]->cv2_par_,
			  cv_set_constraints[ki]->cv2_der_,
			  sign, &kstat);
      if (kstat < 0)
	{
	  return kstat;
	}

      // The two objects are assumed to be on the same side of an equation
      // sign. If der1[ki] or der2[ki] is negative it means that all coefs
      // should be flipped. It is taken care of inside getSideConstraint().

      curr_constraint_set->factor_.insert
	(curr_constraint_set->factor_.end(),
	 side_constraints1.begin(), side_constraints1.end());
      curr_constraint_set->factor_.insert
	(curr_constraint_set->factor_.end(),
	 side_constraints2.begin(), side_constraints2.end());

      side_constraint_sets.push_back(curr_constraint_set);
    }


  kstat = updateSideConstraints(side_constraint_sets, coefknown_);
  if (kstat < 0)
    {
      return kstat;
    }

  if (appr_constraints)
    {  
      kstat = setApproxSideConstraints(side_constraint_sets, appr_wgt);
      if (kstat < 0)
	{
	  return kstat;
	}
    }
  else    
    {  
      setSideConstraints(side_constraint_sets, false);
    }


  return 0;
}


int SmoothCurveSet::equationSolve(vector<shared_ptr<SplineCurve> >& curves)
//--------------------------------------------------------------------------
//     Purpose : Solve linear equation system.
//
//     Calls   : s6lufacp, s6lusolp
//
//     Written by : Sverre Briseid,  SINTEF,  Nov 2005
//--------------------------------------------------------------------------
{
  int kstat = 0;
  int ki, kr;
  vector<int> nlvec(kdim_*kncond_,0); // Pivot_ array used in equation solver.
//   int *nl = nlvec.begin();

  // Solve the equation system by Conjugate Gradient Method.
       
  vector<double> eb(idim_*kncond_);
  std::copy(gright_.begin(), gright_.begin()+idim_*kncond_, eb.begin());
       
  // Copy coefficients to array of unknowns
       
  // Set up CG-object
       
//   SIsolveCG solveCg;

  shared_ptr<SolveCG> solveSS;

  if (knconstraint_ > 0)
    {
      solveSS = 
	shared_ptr<SolveCG>(new SolveCGCO(kncond_ - knconstraint_, knconstraint_));
    }
  else
    {
      solveSS = shared_ptr<SolveCG>(new SolveCG());
    }

       
  int precond = 1;//(knconstraint_ > 0) ? 0 : 1; // true


  // Create sparse matrix.
       
  solveSS->attachMatrix(&gmat_[0], kdim_*kncond_);
       
  // Attach parameters.
       
  solveSS->setTolerance(0.00000001);
  int nmb_iter = precond ? kdim_*kncond_ : 100*kdim_*kncond_;
  // As the preconditioning matrix when using constraints is not
  // optimal (constraint part is the scaled identity) we then raise
  // number of iterations.
  if (precond && knconstraint_ > 0)
    {
      int mult_factor = 5 + (knconstraint_/50);
      nmb_iter *= mult_factor;
    }
  // The max number of iterations should be a function of the number of
  // unknowns.
  solveSS->setMaxIterations(std::min(nmb_iter, 10000));
       
  // Preconditioning
       
  //int precond = 1;
  //        printf("Precondintioning (0/1) ? ");
  //        scanf("%d",&precond);
  if (precond)
    {
      double omega = 0.1;
      // 	   printf("Omega = ");
      // 	   scanf("%lf",&omega);
      solveSS->precondRILU(omega);
    }
       
  // Solve equation systems.
  vector<double> result(idim_*kncond_,0.0);// kk*kncond_,0.0);
  // We fill result with existing coefs as a good starting point makes
  // the job easier for the solver (i.e. if initial curves are
  // somewhat close to the solution).
  int idxcv;
  for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
    {
      bool israt = cvs_[idxcv]->rational();
      vector<double>::iterator scoef;
      if (copyCoef_)
	scoef = coef_array_[idxcv].begin();
      else
	scoef = (israt) ? cvs_[idxcv]->rcoefs_begin() : cvs_[idxcv]->coefs_begin();

      int kn = cvs_[idxcv]->numCoefs();
//        kn12 = kn1*kn2;
      int kj, kk;
      for (kj=0; kj<kn; kj++)
	{
	  if (coefknown_[idxcv][kj] == 1 || coefknown_[idxcv][kj] == 2)
	    continue;
	   
	  int kl1 = (coefknown_[idxcv][kj] > 2) ?
	    pivot_[idxcv][coefknown_[idxcv][kj]-kpointer_] :
	    pivot_[idxcv][kj];

	  for (kk=0; kk<idim_; kk++)
	    result[kk*kncond_+kl1] = scoef[kj*kdim_+kk]; 
	}
    }

  if(kdim_==3)
    {
      kstat=solveSS->solve(&result[0], &gright_[0], kncond_*kdim_);
      if (kstat < 0)
 	return kstat;
    }
  else
    {
      for (kr=0; kr<idim_; kr++)
	{
	  kstat = solveSS->solve(&result[kr*kncond_],
				 &gright_[kr*kncond_], kncond_);
	  //	       printf("solveSS->solve status %d \n", kstat);
	  if (kstat != 0)
	    return kstat;
	}
    }
       
  for (idxcv = 0; idxcv < (int)cvs_.size(); ++idxcv)
    {
      int kn = cvs_[idxcv]->numCoefs();
      int kk = cvs_[idxcv]->order();
      vector<double>::iterator st = cvs_[idxcv]->knotsBegin();
      bool israt = cvs_[idxcv]->rational();
      vector<double>::iterator scoef;
      if (copyCoef_)
	scoef = coef_array_[idxcv].begin();
      else
	scoef = (israt) ? cvs_[idxcv]->rcoefs_begin() : cvs_[idxcv]->coefs_begin();
      // Copy result to output array.	   
      for (ki=0; ki < kn; ++ki)
	{
	  if (coefknown_[idxcv][ki] == 1)
	    continue;
	  int kl1 = (coefknown_[idxcv][ki] > 2) ?
	    pivot_[idxcv][coefknown_[idxcv][ki]-kpointer_] :
	    pivot_[idxcv][ki];
	  for (kr=0; kr<idim_; kr++)
	    {
	      scoef[ki*idim_+kr] = result[kr*kncond_+kl1]; 
	    }
	}
   
   
      // Create curve.
      shared_ptr<SplineCurve> curve = 
	shared_ptr<SplineCurve>(new SplineCurve(kn, kk, st, scoef, idim_, israt));
      curves.push_back(curve);
    }

  return 0;
}

