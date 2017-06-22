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

#include <algorithm>
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/creators/Integrate.h"
#include "GoTools/creators/SolveCG.h"
#include "GoTools/utils/Values.h"

#include <math.h>
#include <fstream>

using namespace Go;
using std::vector;
using std::max;
using std::min;

SmoothSurf::SmoothSurf()
    : kpointer_(3), copy_coefs_(true), omega_(0.1)
   //--------------------------------------------------------------------------
   //     Constructor for class SmoothSurf.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF SI,  09.93.
   //     Revised by : Vibeke Skytt,  SINTEF Oslo, 09.95
   //--------------------------------------------------------------------------
{
   /* Set values of variables global to the class.  */

   ider_scratch_ = ider_ = 3;
   kdim_ = idim_ = 3;
   norm_dim_ = 1;
   cont_seem_[0] = cont_seem_[1] = 0;
   kncond_ = 0;
   knconstraint_ = 0;
   integral1_ = 0;
   integral2_ = 0;
   integralset_ = false;
   rational_ = false;
}

/****************************************************************************/


SmoothSurf::SmoothSurf(bool copy_coefs)
    : kpointer_(3), copy_coefs_(copy_coefs), omega_(0.1)
   //--------------------------------------------------------------------------
   //     Constructor for class SmoothSurf.
   //
   //     Purpose : Initialize class variables
   //
   //     Calls   :
   //
   //     Written by : Vibeke Skytt,  SINTEF SI,  09.93.
   //     Revised by : Vibeke Skytt,  SINTEF Oslo, 09.95
   //--------------------------------------------------------------------------
{
   /* Set values of variables global to the class.  */

   ider_scratch_ = ider_ = 3;
   kdim_ = idim_ = 3;
   norm_dim_ = 1;
   cont_seem_[0] = cont_seem_[1] = 0;
   kncond_ = 0;
   knconstraint_ = 0;
   integral1_ = 0;
   integral2_ = 0;
   integralset_ = false;
   rational_ = false;
}

/****************************************************************************/


SmoothSurf::~SmoothSurf()
   //--------------------------------------------------------------------------
   //     Destructor for class SmoothSurf.
   //
   //     Written by : Vibeke Skytt,  SINTEF SI,  09.93.
   //--------------------------------------------------------------------------
{
  releaseScratch();
}


/****************************************************************************/


void SmoothSurf::releaseScratch()
   //--------------------------------------------------------------------------
   //     Free all memory allocated for the class members of SmoothSurf.
   //
   //     Written by : Vibeke Skytt,  SINTEF SI,  11.99
   //--------------------------------------------------------------------------
{
   int ki;

   if (integral1_ && integral2_)
     {
       for (ki=0; ki<=ider_scratch_; ki++)
	 {
	   if (integral1_[ki]) delete [] integral1_[ki];
	   integral1_[ki] = 0;
	   if (integral2_[ki]) delete [] integral2_[ki];
	   integral2_[ki] = 0;
	 }
       if (integral1_) delete [] integral1_;
       integral1_ = 0;
       if (integral2_) delete [] integral2_;
       integral2_ = 0;
     }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  SmoothSurf::prepareIntegral
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

void SmoothSurf::prepareIntegral()

//--------------------------------------------------------------------------
//
//     Purpose : Prepare storage for integrals of inner products of 
//               basis functions.
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF SI,  08.99.
//     Revised by :
//--------------------------------------------------------------------------
{
   // Allocate scratch for arrays of integrals of inner product of B-splines.

   int ider1 = ider_scratch_ + 1;
   vec1_.resize(ider1*kn1_*kn1_);
   vec2_.resize(ider1*kn2_*kn2_);
   std::fill(vec1_.begin(), vec1_.end(), 0.0);
   std::fill(vec2_.begin(), vec2_.end(), 0.0);

   // Set pointers into arrays.

   int ki, kj;
   integral1_ = new double**[ider1];
   integral2_ = new double**[ider1];
   if (integral1_ == 0 || integral2_ == 0)
     THROW("Allocation error");

   for (ki=0; ki<=ider_scratch_; ki++)
   {
      integral1_[ki] = new double*[kn1_];
      integral2_[ki] = new double*[kn2_];
      if (integral1_[ki] == 0 || integral2_[ki] == 0)
	THROW("Allocation error");

      for (kj=0; kj<kn1_; kj++)
	integral1_[ki][kj] = &vec1_[(ki*kn1_+kj)*kn1_];

      for (kj=0; kj<kn2_; kj++)
	integral2_[ki][kj] = &vec2_[(ki*kn2_+kj)*kn2_];
   }

   return;
}
   
///////////////////////////////////////////////////////////////////////////////
//
//  SmoothSurf::attach
//
///////////////////////////////////////////////////////////////////////////////

void
SmoothSurf::attach(shared_ptr<SplineSurface>& insf,  // Input surface 
		                              // representing the spline space
		   int seem[],       // Expected continuity across seem.
		   int coef_known[], // Indicates whether or not each
		                     // coefficient is known (1),
		                     // unknown (0) or not of interest(2)
		   int num_side_constraints, // Number of side constraints
		                             // to the minimization problem
		   int has_normal_cond)    // Indicates if normal conditions
                                           // will be given.

//--------------------------------------------------------------------------
//
//     Purpose : Initialize class variables related to the spline space.
//
//     Calls   :
//
// NB!
// If attach is used more than once for an instance of SmoothSurf, 
// then it is assumed that the new surface (insf) lies in the same
// spline space that the old surface (srf_) and that there
// is no change concerning the existance of normal conditions.
// If these conditions are broken, memory faults and/or wrong
// results will occur. Including higher order derivatives in the
// smoothing functional than in previous runs, will have no effect.
//
//     Written by : Vibeke Skytt,  SINTEF SI,  09.95.
//     Revised by : Vibeke Skytt,  SINTEF,  11.99.
//--------------------------------------------------------------------------
{
   // int kstat = 0;
   int ki, kj;
   // Fetch data of B-spline surface. 

   idim_ = insf->dimension();
   st1_ = insf->basis_u().begin();
   st2_ = insf->basis_v().begin();
   rational_ = insf->rational();

   if (rational_)
     {
       kdim_ = idim_ + 1;
       if (copy_coefs_)
	 {
	   coef_array_.resize((insf->numCoefs_u())*(insf->numCoefs_v())*kdim_);
	   copy(insf->rcoefs_begin(),
		insf->rcoefs_end(),
		coef_array_.begin());
	   scoef_ = coef_array_.begin();
	 }
       else
	 scoef_ = insf->rcoefs_begin();
       int in_coefs_u = insf->numCoefs_u();
       int in_coefs_v = insf->numCoefs_v();
       int in_ord_u = insf->order_u();
       int in_ord_v = insf->order_v();
       vector<double> bspl_coefs(in_coefs_u * in_coefs_v * 2, 0.0);
       for (int pos_in = idim_, pos_out = 1;
	    pos_in < kdim_ * in_coefs_u * in_coefs_v;
	    pos_in += kdim_, pos_out +=2)
	 bspl_coefs[pos_out] = scoef_[pos_in];
       bspline_surface_ = shared_ptr<SplineSurface>
	 (new SplineSurface(in_coefs_u, in_coefs_v, in_ord_u, in_ord_v,
			    st1_, st2_, bspl_coefs.begin(), 1, true));
     }
   else
     {
       kdim_ = idim_;
       if (copy_coefs_)
	 {
	   coef_array_.resize((insf->numCoefs_u())*(insf->numCoefs_v())*kdim_);
	   copy(insf->coefs_begin(),
		insf->coefs_end(),
		coef_array_.begin());
	   scoef_ = coef_array_.begin();
	 }
       else
	 scoef_ = insf->coefs_begin();
     }

   if (has_normal_cond && idim_ == 3)
     norm_dim_ = 3;          // Prepare for approximation of normal conditions.

   if (srf_.get() != 0 &&  
       insf->numCoefs_u() == kn1_ && insf->numCoefs_v() == kn2_ && 
       insf->order_u() == kk1_ && insf->order_v() == kk2_ &&
       num_side_constraints == knconstraint_)
     {
       // It is assumed that the new surface (insf) lies in the same
       // spline space that the old surface and that there
       // is no change concerning the existance of normal conditions.
       // If these conditions are broken, memory faults and/or wrong
       // results will occur. Including higher order derivatives in the
       // smoothing functional than in previous runs, will have no effect.

       // Zero out the arrays of the equation system.

       std::fill(gmat_.begin(), gmat_.end(), 0.0);
       std::fill(gright_.begin(), gright_.end(), 0.0);

       srf_ = insf;
     }
   else
     {
       kk1_ = insf->order_u();
       kk2_ = insf->order_v();
       kn1_ = insf->numCoefs_u();
       kn2_ = insf->numCoefs_v();

       if (srf_.get() != 0)
	 {
	   // Memory already allocated for the previous iteration of
	   // surface editing and smoothing. Free this memory.

	   releaseScratch(); 
	 }

       integralset_ = false;
       srf_ = insf;

       // Allocate scratch for pivot_ array.
       pivot_.resize(kn1_*kn2_);

       // Update coefknown_ according to information about continuity over
       // a seem. Also update coefficients at the seem if the required 
       // information exists.
   
       coefknown_ = coef_known;
       if (seem[0] > 0 || seem[1] > 0)
	 preparePeriodicity(seem);

       // Store information about free, fixed and irrelevant coefficients
       // and create pivot_ array. Count number of free coefficients.

       kncond_ = 0;
       int kn12 = kn1_*kn2_;
       for (kj=0; kj<kn12; kj+=kn1_)
	 for (ki=0; ki<kn1_; ki++)
	   {
	     if (coefknown_[kj+ki] == 0)
	       {
		 pivot_[kj+ki] = kncond_;
		 kncond_++;
	       }
	     else
	       pivot_[kj+ki] = -MAXINT;
	   }

       // Add side constraints
       knconstraint_ = num_side_constraints;
       kncond_ += knconstraint_; // @@@ VSK. What about constraints with no 
                                 // free coefficients? It seems appropriate to
                                 // adjust values when including constraints
                                 // in system.

       // Allocate scratch for arrays of integrals of inner product of 
       // B-splines.
       prepareIntegral();

       // Allocate scratch for arrays in the equation system. 
       MESSAGE("DEBUG: kncond_: " << kncond_);

       int ksize = norm_dim_*norm_dim_*kncond_*kncond_;
       gmat_.resize(ksize);
       gright_.resize(idim_*kncond_);
       std::fill(gmat_.begin(), gmat_.end(), 0.0);
       std::fill(gright_.begin(), gright_.end(), 0.0);
     }

}



//===========================================================================
void
SmoothSurf::setOptimize(const double wgt1,  /* Weight of 1. order term */
			const double wgt2,  /* Weight of 2. order term */
			const double wgt3)  /* Weight of 3. order term */
//===========================================================================
{
  if (rational_)
    setOptimizeRational(wgt1, wgt2, wgt3);
  else
    setOptimizeNonrational(wgt1, wgt2, wgt3);
}


/***************************************************************************/

void
SmoothSurf::getBasis(const double *sb1, const double *sb2,
		     int kleft1, int kleft2,
		     int ider, double *sbasis)
//--------------------------------------------------------------------------
//     Purpose : Given all non-zero B-spline basis functions, compute the
//               corresponding rational basis functions.
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF,  09.99.
//--------------------------------------------------------------------------
{
   int kk12 = kk1_*kk2_;
   int ider1 = ider + 1;

   int ki, kj, kd1, kd2, kl1, kl2, k1, k2, ks1;
   double *s1 = sbasis;

   ks1 = ider;
   for (kd2=0; kd2<=ider; kd2++, ks1--)
     for (kd1=ks1; kd1<=ider-kd2; kd1++, s1+=kk12)
       for (kj=kleft2-kk2_+1, kl2=kd2, k2=0; k2<kk12;
	    kj++, kl2+=ider1, k2+=kk1_)
	   for (ki=kleft1-kk1_+1, kl1=kd1, k1=0; k1<kk1_;
		ki++, kl1+=ider1, k1++)
	       s1[k2+k1] = sb1[kl1]*sb2[kl2];

   return;
}



/****************************************************************************/

void
SmoothSurf::setLeastSquares(const std::vector<double>&  pnts, // Data points.   
			    const std::vector<double>&  param_pnts, // Parametrization.
			    const std::vector<double>&   pnt_weights,
			    const double wgt)  // Weight of current term.     
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		    from the approximation of data points.
//
//     Calls   : s1220  (SISL) - Compute all non-zero B-splines in a point.
//
//     Written by : Vibeke Skytt,  SINTEF SI,  09.93.
//--------------------------------------------------------------------------
{
    int nmbpoint = (int)pnts.size()/idim_;   // Number of data points. 
  // int kstat = 0;
  int kk;
  int k1, k2, k3, k4, k5, k6, k7, k8;
  int kl1, kl2;
  int kleft1=0, kleft2=0;  // Parameter used in s1220 to be positioned
                           // in the knot vector.                           
  double tz;     // Help variable.  
  double tval;   // Contribution to the matrices of the minimization problem.
  double *sc;    // Pointer into the coefficient array of the original surf.
  double const1 = (double)2.0*wgt;

  // Allocate scratch for B-spline basis functions. 
  
  vector<double> scratch(kk1_+kk2_+kk1_*kk2_, 0.0);
  double *sb1 = &scratch[0];  // Storage of B-spline basis functions. 
  double *sb2 = sb1+kk1_;          // Storage of B-spline basis functions. 
  double *sbasis = sb2+kk2_;       // Surface basis functions.

  // Traverse all points in the pointset.  

  const double *pnt = &pnts[0];
  const double *par = &param_pnts[0];
  for (int kr=0; kr<nmbpoint; kr++, pnt+=idim_, par+=2)
   {
     // Fetch B-spline basis functions different from zero. 

      if (rational_)
	{
	  Point p(1);
	  vector<double>::iterator bspl_it = bspline_surface_->rcoefs_begin();
	  kleft1 = bspline_surface_->basis_u().knotInterval(par[0]);
	  kleft2 = bspline_surface_->basis_v().knotInterval(par[1]);
	  bspl_it += 2 * ((kleft2 -kk2_ + 1) * kn1_ + kleft1 - kk1_ + 1);
	  for (int j = 0; j < kk2_; ++j)
	    {
	      for (int i = 0; i < kk1_; ++i)
		{
		  bspl_it[2*(j*kn1_ + i)] = 1.0;
		  bspline_surface_->point(p, par[0], par[1]);
		  sbasis[j*kk1_+i] = p[0];
		  bspl_it[2*(j*kn1_ + i)] = 0.0;
		}
	    }
	}
      else
	{
	  srf_->basis_u().computeBasisValues(par[0], sb1, 0);
	  srf_->basis_v().computeBasisValues(par[1], sb2, 0);
	  kleft1 = srf_->basis_u().lastKnotInterval();
	  kleft2 = srf_->basis_v().lastKnotInterval();

	  // Compute the surface basis functions.
	  getBasis(sb1, sb2, kleft1, kleft2, 0, sbasis);
	}

     for (k1=kleft1-kk1_+1, k3=0; k1<=kleft1; k1++, k3++)
       for (k2=kleft2-kk2_+1, k4=0; k2<=kleft2; k2++, k4+=kk1_)
	 {
	   // Test if the the coefficient is free to change.

	   if (coefknown_[k2*kn1_+k1] == 1 || coefknown_[k2*kn1_+k1] == 2)
	     continue;

	   kl1 = (coefknown_[k2*kn1_+k1] > 2) ?
	       pivot_[coefknown_[k2*kn1_+k1]-kpointer_] : pivot_[k2*kn1_+k1];

	   tz = pnt_weights[kr]*sbasis[k4+k3];

	   // Add contribution to right hand side. 

	   for (kk=0; kk<idim_; kk++)
	     {
	       tval = const1*pnt[kk]*tz;
 	       gright_[kk*kncond_+kl1] += tval;
	     }

	   for (k5=kleft1-kk1_+1, k7=0; k5<=kleft1; k5++, k7++)
 	     for (k6=kleft2-kk2_+1, k8=0; k6<=kleft2; k6++, k8+=kk1_)
	       {
	         if (coefknown_[k6*kn1_+k5] == 2)
		   continue;

 		 // Compute contribution to left hand side. 

		 tval = const1*tz*sbasis[k8+k7];

		 // Test if the current coefficient is at the boundary.

		 if (coefknown_[k6*kn1_+k5] == 1)
 		   {
		     // Adjust for boundary conditions. 

		     sc = &*scoef_ + (k6*kn1_ + k5)*kdim_;
 		     for (kk=0; kk<idim_; kk++)
 			gright_[kk*kncond_+kl1] -= sc[kk]*tval;
		   }
		 else
		   {
		     // The term gives a contribution on the left hand side. 

		     kl2 = (coefknown_[k6*kn1_+k5] > 2) ?
		       pivot_[coefknown_[k6*kn1_+k5]-kpointer_] : 
		       pivot_[k6*kn1_+k5];

 		     if (kl2 > kl1) continue;

 		     for (kk=0; kk<norm_dim_; kk++)
		       {
			 gmat_[(kk*kncond_+kl1)*norm_dim_*kncond_+kk*kncond_+kl2]
			     += tval;
			 if (kl2 < kl1)
			   gmat_[(kk*kncond_+kl2)*norm_dim_*kncond_+kk*kncond_+kl1]
			       += tval;
		       }
		   }
 	       }
 	 }
   }

  return;
}



/****************************************************************************/

int
SmoothSurf::setNormalCond(const std::vector<double>& pnts,
			  const std::vector<double>& param_pnts,
			  const std::vector<double>&   pnt_weights,
			  const double weight)
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of normal directions.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF,  9910.
//--------------------------------------------------------------------------
{
  if (norm_dim_ != 3)
    return 1;   // Not prepared for approximation of normals.

  int nmbpoint = (int)pnts.size()/idim_;   // Number of normal directions
  // int kstat = 0;
  int kk, kb;
  int k1, k2, k3, k4, k5, k6, k7, k8;
  int kl1, kl2;
  int kleft1=0, kleft2=0;  // Parameter used in s1220 to be positioned
                           // in the knot vector.                           
  int kder = 1;            // Number of derivatives of B-splines to compute.
  int kder1 = kder+1; 
  int kk12 = kk1_*kk2_;
  double tz1, tz2;         // Help variables.  
  double tval;   // Contribution to the matrices of the minimization problem.
  double *sc;    // Pointer into the coefficient array of the original surf.
  double const1 = (double)2.0*weight;
  double tdum;

  // Allocate scratch for B-spline basis functions. 
  
  vector<double> scratch((kk1_+kk2_+kk1_*kk2_*kder1)*kder1, 0.0);
  double *sb1 = &scratch[0];  // Storage of B-spline basis functions. 
  double *sb2 = sb1+kk1_*kder1;    // Storage of B-spline basis functions. 
  double *sbasis = sb2+kk2_*kder1;       // Surface basis functions.

  // Traverse all normal directions.

  const double *pnt = &pnts[0];
  const double *par = &param_pnts[0];
  for (int kr=0; kr<nmbpoint; kr++, pnt+=idim_, par+=2)
   {

     // Fetch B-spline basis functions different from zero. 

     if (rational_)
       {
	 vector<Point> pts(3);
	 vector<double>::iterator bspl_it = bspline_surface_->rcoefs_begin();
	 kleft1 = bspline_surface_->basis_u().knotInterval(par[0]);
	 kleft2 = bspline_surface_->basis_v().knotInterval(par[1]);
	 bspl_it += 2 * ((kleft2 -kk2_ + 1) * kn1_ + kleft1 - kk1_ + 1);
	 int pos_du = 0;
	 int pos_dv = kk1_*kk2_;
	 for (int j = 0; j < kk2_; ++j)
	   {
	     for (int i = 0; i < kk1_; ++i, ++pos_du, ++pos_dv)
	       {
		 bspl_it[2*(j*kn1_ + i)] = 1.0;
		 bspline_surface_->point(pts, par[0], par[1], 1);
		 sbasis[pos_du] = pts[1][0];
		 sbasis[pos_dv] = pts[2][0];
		 bspl_it[2*(j*kn1_ + i)] = 0.0;
	       }
	   }
       }
     else
       {
	 srf_->basis_u().computeBasisValues(par[0], sb1, kder);
	 srf_->basis_v().computeBasisValues(par[1], sb2, kder);
	 kleft1 = srf_->basis_u().lastKnotInterval();
	 kleft2 = srf_->basis_v().lastKnotInterval();

	 // Compute the surface basis functions.
	 getBasis(sb1, sb2, kleft1, kleft2, kder, sbasis);
       }

     for (k1=kleft1-kk1_+1, k3=0; k1<=kleft1; k1++, k3++)
       for (k2=kleft2-kk2_+1, k4=0; k2<=kleft2; k2++, k4+=kk1_)
	 {
	   // Test if the the coefficient is free to be changed.

	   if (coefknown_[k2*kn1_+k1] == 1 || coefknown_[k2*kn1_+k1] == 2)
	     continue;

	   kl1 = (coefknown_[k2*kn1_+k1] > 2) ?
	       pivot_[coefknown_[k2*kn1_+k1]-kpointer_] : pivot_[k2*kn1_+k1];

	   tz1 = pnt_weights[kr]*sbasis[k4+k3];
	   tz2 = pnt_weights[kr]*sbasis[kk12+k4+k3];

	   for (k5=kleft1-kk1_+1, k7=0; k5<=kleft1; k5++, k7++)
 	     for (k6=kleft2-kk2_+1, k8=0; k6<=kleft2; k6++, k8+=kk1_)
 	       {
		  if (coefknown_[k6*kn1_+k5] == 2)
		    continue;

 		  // Compute contribution to left hand side. 

 		  tval = const1*(tz1*sbasis[k7+k8] +
				 tz2*sbasis[kk12+k7+k8]);

 		  // Test if the current coefficient is at the boundary.

		  if (coefknown_[k6*kn1_+k5] == 1)
 		  {
		     // Adjust for boundary conditions. 

 		     sc = &*scoef_ + (k6*kn1_ + k5)*kdim_;
 		     for (tdum=(double)0, kk=0; kk<idim_; kk++)
			tdum += sc[kk]*pnt[kk];

		     for (kk=0; kk<idim_; kk++)
 			gright_[kk*kncond_+kl1] -= tval*tdum*pnt[kk];
 		  }
 		  else
 		  {
		     // The term gives a contribution on the left hand side. 

		     kl2 = (coefknown_[k6*kn1_+k5] > 2) ?
			 pivot_[coefknown_[k6*kn1_+k5]-kpointer_] : 
			 pivot_[k6*kn1_+k5];

 		     if (kl2 > kl1) continue;

 		     for (kk=0; kk<norm_dim_; kk++)
		       {
			 for (kb=0; kb<norm_dim_; kb++)
			   {
			     gmat_[(kk*kncond_+kl1)*norm_dim_*kncond_+
				   kk*kncond_+kl2] +=
				 tval*pnt[kk]*pnt[kb];
			     if (kl2 < kl1)
			       gmat_[(kk*kncond_+kl2)*norm_dim_*kncond_+
				     kk*kncond_+kl1] +=
				   tval*pnt[kk]*pnt[kb];
			   }
 		     }
 		  }
 	       }
 	 }
    }

  return 0;
}

/****************************************************************************/

void
SmoothSurf::approxOrig(double weight)  // Weight of current term.     
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of an original surface, non-rational case
{
  if (rational_)
    approxOrigRational(weight);
  else
    approxOrigNonrational(weight);
}


/****************************************************************************/

void
SmoothSurf::approxOrigNonrational(double weight)  // Weight of current term.     
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of an original surface, non-rational case
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF,  9808 (originally approxOrig())
//--------------------------------------------------------------------------
{
    int kq, kp, ki, kj, kr, kjstart, kjend, kistart, kiend;

    int kl1, kl2;
    double wgt2 = (double)2*weight;
    double innerprod;

    for(kq=0; kq<kn2_; kq++)
	for(kp=0; kp<kn1_; kp++) {
	    if (coefknown_[kq*kn1_+kp] == 1 || coefknown_[kq*kn1_+kp] == 2)
		continue;

	    kl2 = (coefknown_[kq*kn1_+kp] > 2) ?
		pivot_[coefknown_[kq*kn1_+kp]-kpointer_] : pivot_[kq*kn1_+kp];

	    kjstart = max(0, kq-kk2_+1);
	    kjend = min(kq+kk2_,kn2_);

	    for(kj=kjstart; kj<kjend; kj++) {
		kistart = max(0, kp-kk1_+1);
		kiend = min(kp+kk1_,kn1_);
		for(ki=kistart; ki<kiend; ki++) {
		    // coefficient is allready known
		    if(coefknown_[kj*kn1_+ki]==1)
			continue;
		    kl1 = (coefknown_[kj*kn1_+ki] > 2) ?
			pivot_[coefknown_[kj*kn1_+ki]-kpointer_] :
			pivot_[kj*kn1_+ki];

		    innerprod=integral1_[0][ki][kp]*integral2_[0][kj][kq];
		    innerprod*=wgt2;

		    for (kr=0; kr<idim_; kr++)
			gright_[kr*kncond_+kl2] +=
			    innerprod*scoef_[(kj*kn1_+ki)*kdim_+kr];

		    for (kr=0; kr<norm_dim_; kr++) {
			gmat_[(kr*kncond_+kl2)*norm_dim_*kncond_+kr*kncond_+kl1]
			    += innerprod;
		    }
		}
	    }
	}
}


void
SmoothSurf::approxOrigRational(double weight)  // Weight of current term.     
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of an original surface, rational case
//
//     Calls   : 
//
//     Written by : Kjell Fredrik Pettersen,  SINTEF,  2009-12-11
//--------------------------------------------------------------------------
{
  // For the rational case, we can not split the double integral into a product
  // of two separate integrals for each parameter value (like the non-rational case),
  // because of the denominator function depending on both variables.
  // Instead we use the grid evaluator to evaluate each spline product in
  // points for numerical integration by Gaussian quadrature.
  //
  // The algorithm is very similar to the one used in setOptimizeRational()

  // Get points and weights for integration
  vector<double> gq_points_u, gq_points_v;
  vector<double> gq_weights_u, gq_weights_v;
  GaussQuadValues(bspline_surface_->basis_u(), gq_points_u, gq_weights_u);
  GaussQuadValues(bspline_surface_->basis_v(), gq_points_v, gq_weights_v);
  int seg_samples_u = (int)gq_weights_u.size();   // Number of samples inside each Bezier segment in first direction
  int seg_samples_v = (int)gq_weights_v.size();   // Number of samples inside each Bezier segment in second direction
  int bez_segs_u = kn1_ - kk1_ + 1;   // Number of Bezier segments in first direction
  int bez_segs_v = kn2_ - kk2_ + 1;   // Number of Bezier segments in second direction
  int samples_u = bez_segs_u * seg_samples_u;  // Number of samples in total ( == gq_points_u.size() ) in first direction
  int samples_v = bez_segs_v * seg_samples_v;  // Number of samples in total ( == gq_points_v.size() ) in second direction

  // Get B-spline evaluations in each sample point
  vector<double> basisvals_u(samples_u * kk1_);   // The B-spline values for integration points in first direction
  vector<int> left_u(samples_u);
  bspline_surface_->basis_u().computeBasisValues(&gq_points_u[0], &gq_points_u[samples_u], &basisvals_u[0], &left_u[0]);
  vector<double> basisvals_v(samples_v * kk2_);   // The B-spline values for integration points in second direction
  vector<int> left_v(samples_v);
  bspline_surface_->basis_v().computeBasisValues(&gq_points_v[0], &gq_points_v[samples_v], &basisvals_v[0], &left_v[0]);

  // For each pair (Bi,Cj) of B-splines in first and second direction, get all evaluations
  // of the function B_i*C_j/G in all integration sample points where B_i and C_j
  // have support, where G is the denominator function of the surface

  // The results are organized in the vector sample_evaluations. It has one element for each B-spline pair(Bi,Cj)
  // starting with i=0, j=0, then i=1, j=0 etc.
  // For each B-spline pair, the entry is a vector of evaluations, organized in two levels. The first
  // (outermost) level is for each sample value in second paramter directions, The second (innermost)
  // level is for each sample value in first direction.
  vector<vector<double> > sample_evaluations;
  sample_evaluations.resize(kn1_ * kn2_);

  // Some variables to describe the Bezier segments where each B-spline has support
  vector<int> left_bas_u(kn1_);  // The first Bezier segment where each B-spline has support
  vector<int> left_bas_v(kn2_);
  vector<int> segs_bas_u(kn1_);  // The number of Bezier segments where each B-spline has support
  vector<int> segs_bas_v(kn2_);
  for (int i = 0; i < kn1_; ++i)
    {
      left_bas_u[i] = std::max(0, i-kk1_+1);
      int right_bas = std::min(i, kn1_-kk1_);
      segs_bas_u[i] = right_bas - left_bas_u[i] + 1;
    }
  for (int i = 0; i < kn2_; ++i)
    {
      left_bas_v[i] = std::max(0, i-kk2_+1);
      int right_bas = std::min(i, kn2_-kk2_);
      segs_bas_v[i] = right_bas - left_bas_v[i] + 1;
    }

  // Fill sample_evaluations
  vector<double>::iterator bspl_it = bspline_surface_->rcoefs_begin();
  for (int j = 0, samp_ev_pos = 0; j < kn2_; ++j)
    {
      int curr_samp_v = seg_samples_v*segs_bas_v[j];  // Number of sample points where the spline has support
      double* start_basisvals_v = &basisvals_v[left_bas_v[j]*kk2_*seg_samples_v];
      int* start_left_v = &left_v[left_bas_v[j]*seg_samples_v];
      for (int i = 0; i < kn1_; ++i, ++samp_ev_pos, bspl_it+=2)
	{
	  int curr_samp_u = seg_samples_u*segs_bas_u[i];  // Number of sample points where the spline has support
	  double* start_basisvals_u = &basisvals_u[left_bas_u[i]*kk1_*seg_samples_u];
	  int* start_left_u = &left_u[left_bas_u[i]*seg_samples_u];
	  sample_evaluations[samp_ev_pos].resize(curr_samp_u*curr_samp_v);
	  *bspl_it = 1.0;
	  bspline_surface_->pointsGrid(curr_samp_u, curr_samp_v, 0,
				       start_basisvals_u, start_basisvals_v,
				       start_left_u, start_left_v,
				       &sample_evaluations[samp_ev_pos][0]);
	  *bspl_it = 0.0;
	}
    }

  // Travers all B-splines and set up matrices of equation system.
  for (int j1 = 0, pos_1 = 0; j1 < kn2_; ++j1)  // For each B-spline in second direction, first B-spline pair
    for (int i1 = 0; i1 < kn1_; ++i1, ++pos_1)  // For each B-spline in first direction, first B-spline pair
      {
	if (coefknown_[pos_1] == 1 || coefknown_[pos_1] == 2)
	  continue;
	int piv_1 = (coefknown_[pos_1] > 2) ?
	  pivot_[coefknown_[pos_1]-kpointer_] : pivot_[pos_1];

	for (int j2 = 0, pos_2 = 0; j2 < kn2_; ++j2)  // For each B-spline in second direction, second B-spline pair
	  {
	    if (left_bas_v[j1] + segs_bas_v[j1] <= left_bas_v[j2] ||
		left_bas_v[j2] + segs_bas_v[j2] <= left_bas_v[j1])
	      {
		pos_2 += kn1_;
		continue;     // B-splines do not overlap
	      }

	    int first_seg_v = std::max(left_bas_v[j1], left_bas_v[j2]);
	    int end_seg_v = std::min(left_bas_v[j1]+segs_bas_v[j1], left_bas_v[j2]+segs_bas_v[j2]);
	    for (int i2 = 0; i2 < kn1_; ++i2, ++pos_2)  // For each B-spline in first direction, second B-spline pair
	      {
		if (left_bas_u[i1] + segs_bas_u[i1] <= left_bas_u[i2] ||
		    left_bas_u[i2] + segs_bas_u[i2] <= left_bas_u[i1])
		  continue;     // B-splines do not overlap
		if (coefknown_[pos_2] == 1)
		  continue;     // No contribution from fixed coefficient

		int piv_2 = (coefknown_[pos_2] > 2) ?
		  pivot_[coefknown_[pos_2]-kpointer_] : pivot_[pos_2];

		if (pos_2 < pos_1)
		  continue;

		// Now we know we have a contribution, and need to evaluate the integral term

		int first_seg_u = std::max(left_bas_u[i1], left_bas_u[i2]);
		int end_seg_u = std::min(left_bas_u[i1]+segs_bas_u[i1], left_bas_u[i2]+segs_bas_u[i2]);

		double term = 0.0;    // The final contribution for these two B-spline pairs to the equation system

		for (int seg_j = first_seg_v; seg_j < end_seg_v; ++seg_j)   // For each common Bezier segment in second direction
		  for (int seg_i = first_seg_u; seg_i < end_seg_u; ++seg_i)   // For each common Bezier segment in first direction
		    {
		      double area = (st1_[seg_i+kk1_]-st1_[seg_i+kk1_-1]) * (st2_[seg_j+kk2_]-st2_[seg_j+kk2_-1]);
		      vector<double>::const_iterator samp_eval_1
			= sample_evaluations[pos_1].begin() + 
			+ seg_samples_u
			  * (seg_i-left_bas_u[i1]
			     + seg_samples_v * (seg_j-left_bas_v[j1]) * segs_bas_u[i1]);
		      vector<double>::const_iterator samp_eval_2
			= sample_evaluations[pos_2].begin() + 
			+ seg_samples_u
			  * (seg_i-left_bas_u[i2]
			     + seg_samples_v * (seg_j-left_bas_v[j2]) * segs_bas_u[i2]);
		      for (int samp_j = 0; samp_j < seg_samples_v; ++samp_j)   // For each sample in second direction
			{
			  for (int samp_i = 0;          // For each sample in first direction
			       samp_i < seg_samples_u;
			       ++samp_i, ++samp_eval_1, ++samp_eval_2)
			      term +=
				gq_weights_u[samp_i] * gq_weights_v[samp_j]
				* samp_eval_1[0] * samp_eval_2[0] * area;

			  if (samp_j < seg_samples_v-1)
			    {
			      samp_eval_1 += seg_samples_u * (segs_bas_u[i1]-1);
			      samp_eval_2 += seg_samples_u * (segs_bas_u[i2]-1);
			    }
			}  // End -- For each sample in second direction
		    }  // End -- For each common Bezier segment in first and second direction

		term *= 2.0 * weight;

		// Contribution on left side of equation system
		for (int k=0; k<norm_dim_; k++)
		  {
		    gmat_[(k*kncond_+piv_2)*norm_dim_*kncond_+k*kncond_+piv_1]
		      += term;
		    if (pos_1 != pos_2)
		      gmat_[(k*kncond_+piv_1)*norm_dim_*kncond_+k*kncond_+piv_2]
			+= term;
		  }

		// Contribution on right side of equation system
		double* sc = &*scoef_ + pos_2*kdim_;
		for (int k=0; k<idim_; ++k)
		  gright_[k*kncond_+piv_1] += sc[k]*term;
		if (pos_1 != pos_2)
		  {
		    sc = &*scoef_ + pos_1*kdim_;
		    for (int k=0; k<idim_; ++k)
		      gright_[k*kncond_+piv_2] += sc[k]*term;
		  }

	      }  // End -- For each B-spline in first direction, second B-spline pair
	  }  // End -- For each B-spline in first direction, first B-spline pair
      }  // End -- For each B-spline in first and second direction, first B-spline pair

}

/****************************************************************************/
void
SmoothSurf::setPeriodicity(int pardir,  // Parameter direction, 1 or 2.
			   int cont,
			   double weight1,
			   double weight2)
//--------------------------------------------------------------------------
//     Purpose : Set periodicity constraints in one par. dir.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF Oslo,  04.98.
//--------------------------------------------------------------------------
{
    int nmb_constraint = cont + 1;
    // Store continuity requirement
    int kk = (pardir == 1) ? kk1_ : kk2_;
    kk = std::min(kk-1, 3);              // At most C2-continuity currently
    cont_seem_[pardir-1] = std::min(nmb_constraint, kk);

    if (rational_)
      {
	if (cont_seem_[pardir-1] > 2)
	  setRationalCnAtSeem(pardir, 2, weight1, weight2);
	else if (cont_seem_[pardir-1] > 1)
	  setRationalCnAtSeem(pardir, 1, weight1, weight2);
      }
    else
      {
	// Set constraints on C1-continuity
	if (cont_seem_[pardir-1] > 1)
	  setC1AtSeem(pardir, weight1);
  
	// Set constraints on C2-continuity
	if (cont_seem_[pardir-1] > 2)
	  setC2AtSeem(pardir, weight2);
      }
  
  return;
}


/****************************************************************************/

void
SmoothSurf::setSideConstraints(std::vector<sideConstraint>& constraints)
//--------------------------------------------------------------------------
//     Purpose : Set linear side constraints to the minimization problem,
//               using Lagrange multiplier.
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF,  04.02
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
    //      { (A^T)A + wE      C^T } {   x    }     { (A^T)b }
    //      {                      } {        }  =  {        }
    //      {      C            0  } { lambda }     {   d    },
    //
    // where lambda denotes the Lagrange Multipliers.
    // We thus have to add new matrices C, C^T & d to the system.
    ALWAYS_ERROR_IF(norm_dim_ == 3,
	// @@sbr Assuming normal conditions are not included in the system!
	// If this is not the case, all calculations concerning gmat_ must
	// be rewritten. Hence norm_dim_ is assumed to be 1.
		"Not prepared for normal conditions!");

    // We check whether constraints are achievable. If coef in constraint is
    // already known, add to const term on right side.
    // If all coefs are removed,
    // delete constraint. Issue warning if constraint is impossible.
    for (size_t i = 0; i < constraints.size(); ++i) {
	int dim = constraints[i].dim_;
	ALWAYS_ERROR_IF(dim != idim_,
			"Constraint not member of the same geometric "
			"space as surface!");

	for (size_t j = 0; j < constraints[i].factor_.size(); ++j) {
	    int coef_ind = constraints[i].factor_[j].first;
	    if (coefknown_[coef_ind] != 0) {
		// As coef was already known, we add to const term on
		// right side.
		for (int k = 0; k < dim; ++k)
		    constraints[i].constant_term_[k] -=
			scoef_[dim*coef_ind+k]*
			constraints[i].factor_[j].second;
		constraints[i].factor_.erase
		    (constraints[i].factor_.begin() + j);
		--j;
	    }
	}
	// If all coefs in constraint were known, we remove constraint.
	if (constraints[i].factor_.size() == 0) {
	    double sum = 0.0;
	    for (int j = 0; j < dim; ++j)
		sum += constraints[i].constant_term_[j];
	    if (sum != 0.0) // Perform check whether constraint was fulfilled.
		MESSAGE("All coefs already known, "
			"side constraint impossible.");
	    constraints.erase(constraints.begin() + i);
	    --i;
	}
    }

    // If any side constraints were removed, we must update size of matrix.
    if (int (constraints.size()) != knconstraint_) {
	int new_knconstraint = (int)constraints.size();
	int new_kncond = kncond_ - (knconstraint_ - new_knconstraint);
	// For ease of algorithm, we copy matrices to new matrices.
	vector<double> new_gmat(norm_dim_*norm_dim_*new_kncond*new_kncond);
	vector<double> new_gright(idim_*new_kncond);
	//     for (i = 0; i < norm_dim_; ++i) // Treat one dimension at the time.
	for (int i = 0; i < new_kncond; ++i)
	    copy(gmat_.begin() + i*norm_dim_*kncond_,
		 gmat_.begin() + i*norm_dim_*kncond_ + new_kncond,
		 new_gmat.begin() + i*new_kncond);
	for (int i = 0; i < idim_; ++i)
	    copy(gright_.begin() + i*kncond_,
		 gright_.begin() + i*kncond_ + new_kncond,
		 new_gright.begin() + i*new_kncond);
	gmat_ = new_gmat;
	gright_ = new_gright;
	knconstraint_ = new_knconstraint;
	kncond_ = new_kncond;
    }

    // We update values in matries as described in the above equations.
    int nmb_free_coefs = kncond_ - knconstraint_; // Dimension of vector x
                                                  // (# coefs).
    // We start by updating gmat by adding new elements given by
    // side constraints.
    for (size_t i = 0; i < constraints.size(); ++i)
	for (size_t j = 0; j < constraints[i].factor_.size(); ++j) {
	    // We start with gmat_.
	    // We have made  sure that all elements in constraints[i] are free.
	    gmat_[kncond_*(nmb_free_coefs+i)+
		  pivot_[constraints[i].factor_[j].first]] =
		constraints[i].factor_[j].second;
	    gmat_[kncond_*pivot_[constraints[i].factor_[j].first]+
		  nmb_free_coefs+i] =
		constraints[i].factor_[j].second;
	}

    // We next update gright_ by adding const values given by side constraints.
    for (int i = 0; i < idim_; ++i)
	for (int j = 0; j < knconstraint_; ++j)
	    gright_[i*kncond_+nmb_free_coefs+j] =
		constraints[j].constant_term_[i];
}


/****************************************************************************/

int
SmoothSurf::equationSolve(shared_ptr<SplineSurface>& surf)
                                                  // Resulting surface.
//--------------------------------------------------------------------------
//     Purpose : Solve linear equation system.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF Oslo,  10.95.
//--------------------------------------------------------------------------
{
   int kstat = 0;
   int ki, kj, kk;
   int kl1;
   int kn12 = kn1_*kn2_;

#ifdef CREATORS_DEBUG
   if (0) {
       FILE *fp = 0;
       fp = fopen("fA.m", "w");
       fprintf(fp,"A=[ ");
       for (kj=0; kj<kncond_; kj++) {
	   for (ki=0; ki<kncond_; ki++)
	       fprintf(fp, "%18.7f", gmat_[kj*kncond_+ki]);
	   if (kj<kncond_-1) fprintf(fp,"\n");
       }
       fprintf(fp," ]; \n");
       fclose(fp);
   }
#endif // CREATORS_DEBUG

   // Solve the equation system by Conjugate Gradient Method.

   vector<double> eb(idim_*kncond_, 0.0);
   for (ki=0; ki<idim_*kncond_; ki++)
     eb[ki] = gright_[ki];

       // Copy coefficients to array of unknowns

   for (ki=0; ki<kn12; ki+=kn1_)
     for (kj=0; kj<kn1_; kj++)
       {
	 if (coefknown_[ki+kj] == 1 || coefknown_[ki+kj] == 2)
	   continue;
	   
	 kl1 = (coefknown_[ki+kj] > 2) ?
	   pivot_[coefknown_[ki+kj]-kpointer_] : pivot_[ki+kj];

	 for (kk=0; kk<idim_; kk++)
	   gright_[kk*kncond_+kl1] = scoef_[(ki+kj)*kdim_+kk]; 
       }
       
   // Set up CG-object

   SolveCG solveCg;

   // Create sparse matrix.

   ASSERT(gmat_.size() > 0);
   solveCg.attachMatrix(&gmat_[0], norm_dim_*kncond_);

   // Attach parameters.

   solveCg.setTolerance(0.00000001);

   // Preconditioning
   // @@sbr When using side-constraints our system is no longer guaranteed to be
   //       symmetric and positive definite! We should then use a different solver.
   int precond = (knconstraint_ > 0 ) ? 0 : 1;
   int nmb_iter = precond ? norm_dim_*kncond_ : 2*norm_dim_*kncond_;
   const int max_iter = 10000; // @@sbr201706 This limit depends on hardware & max allowed runtime ... As
                               // of now 10000 seems to make sense (was 1000).
   if (nmb_iter > max_iter)
   {
       MESSAGE("nmb_iter is reduced from " << nmb_iter << " to " << max_iter);
       nmb_iter = max_iter;
   }
   solveCg.setMaxIterations(nmb_iter);
   if (precond) {
       solveCg.precondRILU(omega_);
   }

   // Solve equation systems.
       
   if (norm_dim_ == 3)
     {
       kstat = solveCg.solve(&gright_[0], &eb[0], norm_dim_*kncond_);
       //	   printf("solveCg.solve status %d \n", kstat);
       if (kstat < 0) 
	 return kstat;
       if (kstat == 1)
	 // @@ Handle this from the outside.
	 THROW("Failed solving system (within tolerance)!");
     }
   else
     {
       for (kk=0; kk<idim_; kk++)
	 {
           kstat = solveCg.solve(&gright_[kk*kncond_], &eb[kk*kncond_],
				 kncond_);
	   //	       printf("solveCg.solve status %d \n", kstat);
	   if (kstat < 0)
	     return kstat;
	   if (kstat == 1)
           {
             // MESSAGE("Tolerance failure, continuing nonetheless! dim = " << kk);
             // continue; // For debugging, to visualize the failed result.
             THROW("Failed solving system (within tolerance)!");
           }
	 }
     }

   // Copy result to output array. 
   // @@@ VSK. Is it necessary to take the Lagrange multipliers
   // into acount?
   // The Lagrange multipliers are given by the last knconstraint_ elements in
   // solution vector.
   // By iterating through the first kncond_ - knconstraint_ elements,
   // there should not be any danger of mixing things up.

   for (ki=0; ki<kn12; ki+=kn1_)
     for (kj=0; kj<kn1_; kj++)
       {
	 if (coefknown_[ki+kj] == 1 || coefknown_[ki+kj] == 2)
	   continue;
	   
	 kl1 = (coefknown_[ki+kj] > 2) ?
	     pivot_[coefknown_[ki+kj]-kpointer_] : pivot_[ki+kj];

	 for (kk=0; kk<idim_; kk++)
  	   scoef_[(ki+kj)*kdim_+kk] = gright_[kk*kncond_+kl1]; 
       }

   // Make sure that continuity requirements at a seem are satisfied.
   if ((cont_seem_[0] > 0 || cont_seem_[1] > 0) && !rational_) {
       kstat = adjustAtSeem();
       if (kstat < 0)
	   return kstat;
   }

   // Create SplineSurface. 
   surf = shared_ptr<SplineSurface>(new SplineSurface(kn1_,kn2_,kk1_,kk2_,
						      st1_,st2_,scoef_,idim_,rational_));

   return 0;
}


/***************************************************************************/


void SmoothSurf::setRelaxParam(double omega)
{
    // We make sure that the omega is inside valid range.
    if (omega < 0.0)
    {
        omega = 0.0;
    }
    if (omega > 1.0)
    {
        omega = 1.0;
    }
    
    omega_ = omega;
}


/***************************************************************************/


void
SmoothSurf::setOptimizeNonrational(const double wgt1,  /* Weight of 1. order term */
				   const double wgt2,  /* Weight of 2. order term */
				   const double wgt3)  /* Weight of 3. order term */
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		    from the smoothness functional for the non-rational case.
//
//     Calls   : SmoothSurf::GaussQuad
//               s1220   (SISL) - Compute derivatives of B-splines
//
//     Written by : Vibeke Skytt,  SINTEF SI,  09-10.93.
//                  Renamed from setOptimize to setOptimizeNonrational 2009-12-02
//--------------------------------------------------------------------------
{
    //int grstat = 0;   // Initialize status variable.
   int ki, kj, kk, kp, kq, kr;
   int k1, k2, k3, k4;
   int kjstart, kjend;
   int kistart, kiend;

   int kleft = 0;    /* Parameter used in 1220 to be positioned in
			the knot vector.                            */
   int kstart;
   int kl1, kl2;
   int order = std::max(kk1_,kk2_);
   double const1 = wgt1*M_PI;
   double const2 = wgt2*M_PI/(double)8.0;
   double const3 = wgt3*M_PI/(double)32.0;

   // Modify number of derivatives to compute.
   int der = ider_;
   if (der == 3 && wgt3 == 0.0)
     der--;
   if (der == 2 && wgt2 == 0.0)
     der--;
   if (der == 1 && wgt1 == 0.0)
     der--;
   ider_ = der;

   // Set parameter area.

   double ta1 = st1_[kk1_-1];   /* Start of par. interval in 1. par. dir. */
   double ta2 = st1_[kn1_];     /* End of par. interval in 1. par. dir.   */
   double tb1 = st2_[kk2_-1];   /* Start of par. interval in 2. par. dir. */
   double tb2 = st2_[kn2_];     /* End of par. interval in 2. par. dir.   */

   double tval;    /* Value of the complete smoothness functional. */
   double tval1 = 0.0, tval2 = 0.0, tval3 = 0.0;  // Contributions to the 
                                                  // smoothness functional.
   double *sc;  // Pointer into coefficient array of the original surface. 

   int ksz1 = std::min(kk1_, kn1_-kk1_) + kk1_;
   int ksz2 = std::min(kk2_, kn2_-kk2_) + kk2_;
   int kk11 = ksz1*ksz1+1;
   int kk22 = ksz2*ksz2+1;
   vector<double> scratch(4*kk11 + 4*kk22 + 3*order, 0.0);

   double *boundary1 = &scratch[0];  // Contribution to the smoothness 
                                         // functional from the boundary 
                                         // of the surface in 1. par. dir. 
   double *boundary2 = boundary1 + 4*kk11;  // Contribution to the 
                                   // smoothness functional from the boundary 
                                   // of the surface in 2. par. dir.     


   // Compute all integrals of inner product of B-splines.

   //Integrate GaussQuad;

   if (integralset_ == false)
     {
	 GaussQuadInner(srf_->basis_u(), ider_, ta1, ta2, integral1_);
	 GaussQuadInner(srf_->basis_v(), ider_, tb1, tb2, integral2_);
	 integralset_ = true;
     }

   double *sbder = boundary2 + 4*kk22;  // Storage of derivatives of B-splines
				        // at the boundary of the surface.  

   // Compute boundary contibutions to the matrices
   // Compute all derivatives of B-splines up to order ider_-1
   // at the start of the parameter interval in 1. par. dir.

   srf_->basis_u().computeBasisValues(ta1, sbder, ider_-1);

   for (k1=0; k1<kk1_; k1++)
     for (k2=0; k2<kk1_; k2++)
     {
       // Compute contribution from the 1. boundary in 1. par. dir.

       boundary1[k1*ksz1+k2] -= sbder[k1*ider_]*sbder[1+k2*ider_];
       boundary1[kk11+k1*ksz1+k2] -= sbder[1+k1*ider_]*sbder[k2*ider_];
       boundary1[2*kk11+k1*ksz1+k2] -= sbder[1+k1*ider_]*sbder[2+k2*ider_];
       boundary1[3*kk11+k1*ksz1+k2] -= sbder[2+k1*ider_]*sbder[1+k2*ider_];
     }

   // Compute all derivatives of B-splines up to order ider_-1
   // at the end of the parameter interval in 1. par. dir.

   for (kr=0; kr<order*ider_; kr++)
     sbder[kr] = (double)0;

   srf_->basis_u().computeBasisValues(ta2, sbder, ider_-1);
   kleft = srf_->basis_u().lastKnotInterval();

   kstart = std::min(kleft-kk1_+1,kk1_);
   for (ki=kstart, k1=0; k1<kk1_; ki++, k1++)
     for (kp=kstart, k2=0; k2<kk1_; kp++, k2++)
     {
       // Compute contribution from the 2. boundary in 1. par. dir.

      boundary1[ki*ksz1+kp] += sbder[k1*ider_]*sbder[1+k2*ider_];
      boundary1[kk11+ki*ksz1+kp] += sbder[1+k1*ider_]*sbder[k2*ider_];
      boundary1[2*kk11+ki*ksz1+kp] += sbder[1+k1*ider_]*sbder[2+k2*ider_];
      boundary1[3*kk11+ki*ksz1+kp] += sbder[2+k1*ider_]*sbder[1+k2*ider_];
     }

   // Compute all derivatives of B-splines up to order ider_-1
   // at the start of the parameter interval in 2. par. dir.

   for (kr=0; kr<order*ider_; kr++)
     sbder[kr] = (double)0;


   srf_->basis_v().computeBasisValues(tb1, sbder, ider_-1);

   for (k3=0; k3<kk2_; k3++)
     for (k4=0; k4<kk2_; k4++)
     {
       // Compute contribution from the 1. boundary in 2. par. dir.

       boundary2[k3*ksz2+k4] -= sbder[k3*ider_]*sbder[1+k4*ider_];
       boundary2[kk22+k3*ksz2+k4] -= sbder[1+k3*ider_]*sbder[k4*ider_];
       boundary2[2*kk22+k3*ksz2+k4] -= sbder[1+k3*ider_]*sbder[2+k4*ider_];
       boundary2[3*kk22+k3*ksz2+k4] -= sbder[2+k3*ider_]*sbder[1+k4*ider_];
     }

   // Compute all derivatives of B-splines up to order ider_-1
   // at the end of the parameter interval in 2. par. dir.

   for (kr=0; kr<order*ider_; kr++)
     sbder[kr] = (double)0;

   srf_->basis_v().computeBasisValues(tb2, sbder, ider_-1);
   kleft = srf_->basis_v().lastKnotInterval();

   kstart = std::min(kleft-kk2_+1,kk2_);
   for (kj=kstart, k3=0; k3<kk2_; kj++, k3++)
     for (kq=kstart, k4=0; k4<kk2_; kq++, k4++)
     {
       // Compute contribution from the 2. boundary in 2. par. dir.

       boundary2[kj*ksz2+kq] += sbder[k3*ider_]*sbder[1+k4*ider_];
       boundary2[kk22+kj*ksz2+kq] += sbder[1+k3*ider_]*sbder[k4*ider_];
       boundary2[2*kk22+kj*ksz2+kq] += sbder[1+k3*ider_]*sbder[2+k4*ider_];
       boundary2[3*kk22+kj*ksz2+kq] += sbder[2+k3*ider_]*sbder[1+k4*ider_];
     }

   // Travers all B-splines and set up matrices of equation system.

   for (kl2=0, kq=0; kq<kn2_; kq++)
     for (kp=0; kp<kn1_; kp++)
       {
	 if (coefknown_[kq*kn1_+kp] == 1 || coefknown_[kq*kn1_+kp] == 2)
	   continue;

	 kl2 = (coefknown_[kq*kn1_+kp] > 2) ? 
	     pivot_[coefknown_[kq*kn1_+kp]-kpointer_] : pivot_[kq*kn1_+kp];

	 kjstart = std::max(0, kq-kk2_+1);
	 kjend = std::min(kq+kk2_,kn2_);

	 for (kl1=0, kj=kjstart; kj<kjend; kj++)
	   {
	     if (kj<kk2_ && kq<kk2_)
	       {
		 k3 = kj; k4 = kq;
	       }
	     else if (kj >= kn2_-kk2_ && kq >= kn2_-kk2_)
	       {
		 k3 = kj - kn2_ + kk2_ + std::min(kk2_, kn2_-kk2_);
		 k4 = kq - kn2_ + kk2_ + std::min(kk2_, kn2_-kk2_);
	       }
	     else
	       {
		 k3 = 0; k4 = kk22-1;
	       }

	     kistart = std::max(0, kp-kk1_+1);
	     kiend = std::min(kp+kk1_,kn1_);

	     for (ki=kistart; ki<kiend; ki++)
	      {
	       if (coefknown_[kj*kn1_+ki] == 2)
		 continue;

	       if (ki<kk1_ && kp<kk1_)
		 {
		   k1 = ki; k2 = kp;
		 }
	       else if (ki >= kn1_-kk1_ && kp >= kn1_-kk1_)
		 {
		   k1 = ki - kn1_ + kk1_ + std::min(kk1_, kn1_-kk1_);
		   k2 = kp - kn1_ + kk1_ + std::min(kk1_, kn1_-kk1_);
		 }
	       else
		 {
		   k1 = 0; k2 = kk11-1;
		 }

	       kl1 = (coefknown_[kj*kn1_+ki] > 2) ? 
		   pivot_[coefknown_[kj*kn1_+ki]-kpointer_] :
		   pivot_[kj*kn1_+ki];

	       if (kl2 > kl1 && coefknown_[kj*kn1_+ki] != 1) 
		 continue;

	       // Compute an element in the left side matrix.

	       // Compute contribution from 3. order term.

	       if (ider_ > 2)
		 tval3 = 5.0*(integral1_[3][ki][kp]*integral2_[0][kj][kq]
			      + integral1_[0][ki][kp]*integral2_[3][kj][kq])
		     +9.0*(integral1_[2][ki][kp]*integral2_[1][kj][kq]
			   + integral1_[1][ki][kp]*integral2_[2][kj][kq])
		     +3.0*((boundary1[3*kk11+k1*ksz1+k2]-integral1_[2][ki][kp])
			   *(boundary2[k3*ksz2+k4]-integral2_[1][kj][kq])
			   +(boundary1[2*kk11+k1*ksz1+k2]-
			     integral1_[2][ki][kp])
			   *(boundary2[kk22+k3*ksz2+k4]-integral2_[1][kj][kq])
			   + (boundary1[kk11+k1*ksz1+k2]-integral1_[1][ki][kp])
			   *(boundary2[2*kk22+k3*ksz2+k4]-
			     integral2_[2][kj][kq])
			   + (boundary1[k1*ksz1+k2]-integral1_[1][ki][kp])
			   *(boundary2[3*kk22+k3*ksz2+k4]-
			     integral2_[2][kj][kq]));

		  // Compute contribution from 2. order term.

	       if (ider_ > 1)
		 tval2 = 4.0*integral1_[1][ki][kp]*integral2_[1][kj][kq]
		     + 3.0*(integral1_[2][ki][kp]*integral2_[0][kj][kq]
			    + integral1_[0][ki][kp]*integral2_[2][kj][kq])
		     + (boundary1[kk11+k1*ksz1+k2]-integral1_[1][ki][kp])
		     *(boundary2[k3*ksz2+k4]-integral2_[1][kj][kq])
		     + (boundary1[k1*ksz1+k2]-integral1_[1][ki][kp])
		     *(boundary2[kk22+k3*ksz2+k4]-integral2_[1][kj][kq]);

	       // Compute contribution from 1. order term.

	       if (ider_ > 0)
		 tval1 = integral1_[1][ki][kp]*integral2_[0][kj][kq]
		     + integral1_[0][ki][kp]*integral2_[1][kj][kq];

	       // Complete term of the smoothness function.

	       tval = const1*tval1 + const2*tval2 + 2.0*const3*tval3;

	       if (coefknown_[kj*kn1_+ki] == 1)
	       {
		  // The contribution of this term is added to the right
		  // side of the equation system. First fetch the known
		  // coefficient.

		  sc = &*scoef_ + (kj*kn1_ + ki)*kdim_;
		  for (kr=0; kr<idim_; kr++)
		    gright_[kr*kncond_+kl2] -= sc[kr]*tval;
	       }
	       else
	       {
		  // The contribution of this term is added to the left
		  //  side of the equation system.

		  if (kl2 > kl1) continue;

		  for (kk=0; kk<norm_dim_; kk++)
		  {
		     gmat_[(kk*kncond_+kl1)*norm_dim_*kncond_+kk*kncond_+kl2]
			 += tval;
		     if (kl2 < kl1)
		       gmat_[(kk*kncond_+kl2)*norm_dim_*kncond_+kk*kncond_+kl1]
			   += tval;
		  }
	       }
	      }
	    }
       }

   return;
}



/***************************************************************************/

void
SmoothSurf::setOptimizeRational(const double wgt1,  /* Weight of 1. order term */
				const double wgt2,  /* Weight of 2. order term */
				const double wgt3)  /* Weight of 3. order term */
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		    from the smoothness functional for the rational case.
//
//     Written by : Kjell Fredrik Pettersen, SINTEF IKT, 2009-12-02
//--------------------------------------------------------------------------
{
  // For the rational case, we can not split the double integral into a product
  // of two separate integrals for each parameter value (like the non-rational case),
  // because of the denominator function depending on both variables.
  // Instead we use the grid evaluator to evaluate each spline product in
  // points for numerical integration by Gaussian quadrature.

  int derivs = 3;
  if (wgt3 == 0.0)
    {
      derivs = 2;
      if (wgt2 == 0.0)
	{
	  derivs = 1;
	  if (wgt1 == 0.0)
	    return;
	}
    }

  double const1 = wgt1*M_PI;
  double const2 = wgt2*M_PI/(double)8.0;
  double const3 = wgt3*M_PI/(double)32.0;

  // Get points and weights for integration
  vector<double> gq_points_u, gq_points_v;
  vector<double> gq_weights_u, gq_weights_v;
  GaussQuadValues(bspline_surface_->basis_u(), gq_points_u, gq_weights_u);
  GaussQuadValues(bspline_surface_->basis_v(), gq_points_v, gq_weights_v);
  int seg_samples_u = (int)gq_weights_u.size();   // Number of samples inside each Bezier segment in first direction
  int seg_samples_v = (int)gq_weights_v.size();   // Number of samples inside each Bezier segment in second direction
  int bez_segs_u = kn1_ - kk1_ + 1;   // Number of Bezier segments in first direction
  int bez_segs_v = kn2_ - kk2_ + 1;   // Number of Bezier segments in second direction
  int samples_u = bez_segs_u * seg_samples_u;  // Number of samples in total ( == gq_points_u.size() ) in first direction
  int samples_v = bez_segs_v * seg_samples_v;  // Number of samples in total ( == gq_points_v.size() ) in second direction

  // Get B-spline evaluations in each sample point
  vector<double> basisvals_u(samples_u * kk1_ * (derivs + 1));   // The B-spline values and derivatives for integration points in first direction
  vector<int> left_u(samples_u);
  bspline_surface_->basis_u().computeBasisValues(&gq_points_u[0], &gq_points_u[samples_u], &basisvals_u[0], &left_u[0], derivs);
  vector<double> basisvals_v(samples_v * kk2_ * (derivs + 1));   // The B-spline values and derivatives for integration points in second direction
  vector<int> left_v(samples_v);
  bspline_surface_->basis_v().computeBasisValues(&gq_points_v[0], &gq_points_v[samples_v], &basisvals_v[0], &left_v[0], derivs);

  // For each pair (Bi,Cj) of B-splines in first and second direction, get all needed combinations of
  // partial derivatives of the function B_i*C_j/G in all integration sample points where B_i and C_j
  // have support, where G is the denominator function of the surface

  // The results are organized in the vector sample_evaluations. It has one element for each B-spline pair(Bi,Cj)
  // starting with i=0, j=0, then i=1, j=0 etc.
  // For each B-spline pair, the entry is a vector of evaluations, organized in three levels. The first
  // (outermost) level is for each sample value in second paramter directions. The second level is for each
  // sample value in first direction. The third (innermost) level is for all derivatives in the given sample value
  // pair, starting with the function value, then d_u, d_v, d_uu, d_uv, d_vv, d_uuu etc.
  vector<vector<double> > sample_evaluations;
  sample_evaluations.resize(kn1_ * kn2_);

  // Some variables to describe the Bezier segments where each B-spline has support
  vector<int> left_bas_u(kn1_);  // The first Bezier segment where each B-spline has support
  vector<int> left_bas_v(kn2_);
  vector<int> segs_bas_u(kn1_);  // The number of Bezier segments where each B-spline has support
  vector<int> segs_bas_v(kn2_);
  for (int i = 0; i < kn1_; ++i)
    {
      left_bas_u[i] = std::max(0, i-kk1_+1);
      int right_bas = std::min(i, kn1_-kk1_);
      segs_bas_u[i] = right_bas - left_bas_u[i] + 1;
    }
  for (int i = 0; i < kn2_; ++i)
    {
      left_bas_v[i] = std::max(0, i-kk2_+1);
      int right_bas = std::min(i, kn2_-kk2_);
      segs_bas_v[i] = right_bas - left_bas_v[i] + 1;
    }

  // Fill sample_evaluations
  int part_derivs = (derivs+1)*(derivs+2)/2;
  vector<double>::iterator bspl_it = bspline_surface_->rcoefs_begin();
  for (int j = 0, samp_ev_pos = 0; j < kn2_; ++j)
    {
      int curr_samp_v = seg_samples_v*segs_bas_v[j];  // Number of sample points where the spline has support
      double* start_basisvals_v = &basisvals_v[left_bas_v[j]*kk2_*(derivs+1)*seg_samples_v];
      int* start_left_v = &left_v[left_bas_v[j]*seg_samples_v];
      for (int i = 0; i < kn1_; ++i, ++samp_ev_pos, bspl_it+=2)
	{
	  int curr_samp_u = seg_samples_u*segs_bas_u[i];  // Number of sample points where the spline has support
	  double* start_basisvals_u = &basisvals_u[left_bas_u[i]*kk1_*(derivs+1)*seg_samples_u];
	  int* start_left_u = &left_u[left_bas_u[i]*seg_samples_u];
	  sample_evaluations[samp_ev_pos].resize(curr_samp_u*curr_samp_v*part_derivs);
	  *bspl_it = 1.0;
	  bspline_surface_->pointsGrid(curr_samp_u, curr_samp_v, derivs,
				       start_basisvals_u, start_basisvals_v,
				       start_left_u, start_left_v,
				       &sample_evaluations[samp_ev_pos][0]);
	  *bspl_it = 0.0;
	}
    }

  // Travers all B-splines and set up matrices of equation system.
  for (int j1 = 0, pos_1 = 0; j1 < kn2_; ++j1)  // For each B-spline in second direction, first B-spline pair
    for (int i1 = 0; i1 < kn1_; ++i1, ++pos_1)  // For each B-spline in first direction, first B-spline pair
      {
	if (coefknown_[pos_1] == 1 || coefknown_[pos_1] == 2)
	  continue;
	int piv_1 = (coefknown_[pos_1] > 2) ?
	  pivot_[coefknown_[pos_1]-kpointer_] : pivot_[pos_1];

	for (int j2 = 0, pos_2 = 0; j2 < kn2_; ++j2)  // For each B-spline in second direction, second B-spline pair
	  {
	    if (left_bas_v[j1] + segs_bas_v[j1] <= left_bas_v[j2] ||
		left_bas_v[j2] + segs_bas_v[j2] <= left_bas_v[j1])
	      {
		pos_2 += kn1_;
		continue;     // B-splines do not overlap
	      }

	    int first_seg_v = std::max(left_bas_v[j1], left_bas_v[j2]);
	    int end_seg_v = std::min(left_bas_v[j1]+segs_bas_v[j1], left_bas_v[j2]+segs_bas_v[j2]);
	    for (int i2 = 0; i2 < kn1_; ++i2, ++pos_2)  // For each B-spline in first direction, second B-spline pair
	      {
		if (left_bas_u[i1] + segs_bas_u[i1] <= left_bas_u[i2] ||
		    left_bas_u[i2] + segs_bas_u[i2] <= left_bas_u[i1])
		  continue;     // B-splines do not overlap
		int piv_2 = (coefknown_[pos_2] > 2) ?
		  pivot_[coefknown_[pos_2]-kpointer_] : pivot_[pos_2];

		if (piv_2 < piv_1 && coefknown_[pos_2] != 1)
		  continue;

		// Now we know we have a contribution, and need to evaluate all integral terms

		int first_seg_u = std::max(left_bas_u[i1], left_bas_u[i2]);
		int end_seg_u = std::min(left_bas_u[i1]+segs_bas_u[i1], left_bas_u[i2]+segs_bas_u[i2]);

		double term = 0.0;    // The final contribution for these two B-spline pairs to the equation system

		for (int seg_j = first_seg_v; seg_j < end_seg_v; ++seg_j)   // For each common Bezier segment in second direction
		  for (int seg_i = first_seg_u; seg_i < end_seg_u; ++seg_i)   // For each common Bezier segment in first direction
		    {
		      double area = (st1_[seg_i+kk1_]-st1_[seg_i+kk1_-1]) * (st2_[seg_j+kk2_]-st2_[seg_j+kk2_-1]);
		      vector<double>::const_iterator samp_eval_1
			= sample_evaluations[pos_1].begin() + 
			+ part_derivs * seg_samples_u
			  * (seg_i-left_bas_u[i1]
			     + seg_samples_v * (seg_j-left_bas_v[j1]) * segs_bas_u[i1]);
		      vector<double>::const_iterator samp_eval_2
			= sample_evaluations[pos_2].begin() + 
			+ part_derivs * seg_samples_u
			  * (seg_i-left_bas_u[i2]
			     + seg_samples_v * (seg_j-left_bas_v[j2]) * segs_bas_u[i2]);
		      for (int samp_j = 0; samp_j < seg_samples_v; ++samp_j)   // For each sample in second direction
			{
			  for (int samp_i = 0;          // For each sample in first direction
			       samp_i < seg_samples_u;
			       ++samp_i, samp_eval_1 += part_derivs, samp_eval_2 += part_derivs)
			    {
			      // Compute contribution in given sample pair
			      double eval_1 = 0.0, eval_2 = 0.0, eval_3 = 0.0;

			      // Third order derivatives
			      if (derivs > 2)
				eval_3 =
				  // d_uuu^2 and d_vvv^2
				  5.0 * (samp_eval_1[6]*samp_eval_2[6] + samp_eval_1[9]*samp_eval_2[9])
				  // d_uuv^2 and d_uvv^2
				  + 9.0 * (samp_eval_1[7]*samp_eval_2[7] + samp_eval_1[8]*samp_eval_2[8])
				  // d_uuv*dvvv and d_uuu*d_uvv
				  + 3.0 * (samp_eval_1[7]*samp_eval_2[9] + samp_eval_1[9]*samp_eval_2[7]
					   + samp_eval_1[6]*samp_eval_2[8] + samp_eval_1[8]*samp_eval_2[6]);

			      // Second order derivatives
			      if (derivs > 1)
				eval_2 =
				  // d_uv^2
				  4.0 * samp_eval_1[4]*samp_eval_2[4]
				  // d_uu^2 and d_vv^2
				  + 3.0 * (samp_eval_1[3]*samp_eval_2[3] + samp_eval_1[5]*samp_eval_2[5])
				  // d_uu*d_vv
				  + samp_eval_1[3]*samp_eval_2[5] + samp_eval_1[5]*samp_eval_2[3];

			      // First order derivatives
			      if (derivs > 0)
				eval_1 =
				  // d_u^2 and d_v^2
				  samp_eval_1[1]*samp_eval_2[1] + samp_eval_1[2]*samp_eval_2[2];

			      // Sum together, and add to term
			      term +=
				(const1*eval_1 + const2*eval_2 + 2.0*const3*eval_3)
				* gq_weights_u[samp_i] * gq_weights_v[samp_j] * area;
			    }  // End -- For each sample in first direction

			  if (samp_j < seg_samples_v-1)
			    {
			      samp_eval_1 += part_derivs * seg_samples_u * (segs_bas_u[i1]-1);
			      samp_eval_2 += part_derivs * seg_samples_u * (segs_bas_u[i2]-1);
			    }
			}  // End -- For each sample in second direction
		    }  // End -- For each common Bezier segment in first and second direction

		// Add contribution to equation system
		if (coefknown_[pos_2] == 1)
		  {
		    // The contribution of this term is added to the right
		    // side of the equation system.

		    double* sc = &*scoef_ + pos_2*kdim_;
		    for (int k=0; k<idim_; ++k)
		      gright_[k*kncond_+piv_1] -= sc[k]*term;
		  }
		else
		  {
		    // The contribution of this term is added to the left
		    //  side of the equation system.
		    for (int k=0; k<norm_dim_; k++)
		      {
			gmat_[(k*kncond_+piv_2)*norm_dim_*kncond_+k*kncond_+piv_1]
			  += term;
			if (piv_1 != piv_2)
			  gmat_[(k*kncond_+piv_1)*norm_dim_*kncond_+k*kncond_+piv_2]
			    += term;
		      }
		  }

	      }  // End -- For each B-spline in first direction, second B-spline pair
	  }  // End -- For each B-spline in second direction, second B-spline pair
      }  // End -- For each B-spline in first and second direction, first B-spline pair

}


/****************************************************************************/

int
SmoothSurf::adjustAtSeem()
//--------------------------------------------------------------------------
//     Purpose : Ensure that the expected continuity at the seem
//               is satisfied.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF Oslo,  04.99.
//--------------------------------------------------------------------------
{
  int ki, kj, kk;
  double t1, t2;
  double tmean, frac, td;

  if (cont_seem_[0] > 1)
    {
      double del1 = 1.0/(st1_[kk1_] - st1_[1]);
      double del4 = 1.0/(st1_[kn1_+kk1_-2] - st1_[kn1_-1]);
      frac = (st1_[kk1_] - st1_[1])/(st1_[kn1_+kk1_-2] - st1_[kn1_-1] + 
				  st1_[kk1_] - st1_[1]);

      // 1. derivative.
      for (kj=0; kj<kn2_; kj++)
	{
	  if ((coefknown_[kj*kn1_+1] == 1 || coefknown_[kj*kn1_+1] == 2) && 
	      (coefknown_[kj*kn1_+kn1_-2] == 1 ||
	       coefknown_[kj*kn1_+kn1_-2] == 2))
	    continue;

	  t1 = (coefknown_[kj*kn1_+1] == 1) ? 1.0 :
	     ((coefknown_[kj*kn1_+kn1_-2] == 1) ? 0.0 : frac);
	  t2 = 1.0 - t1;

	  for (kk=0; kk<idim_; kk++)
	    {
	      tmean = t1*(scoef_[(kj*kn1_+1)*kdim_+kk] - 
			  scoef_[kj*kn1_*kdim_+kk])*del1 +
		  t2*(scoef_[(kj*kn1_+kn1_-1)*kdim_+kk] - 
		      scoef_[(kj*kn1_+kn1_-2)*kdim_+kk])*del4;
	      td = (tmean+del1*scoef_[kj*kn1_*kdim_+kk])/del1;
	      scoef_[(kj*kn1_+1)*kdim_+kk] = td;
	      td = (del4*scoef_[(kj*kn1_+kn1_-1)*kdim_+kk]-tmean)/del4;
	      scoef_[(kj*kn1_+kn1_-2)*kdim_+kk] = td;
	    }
	}

      if (cont_seem_[0] > 2)
	{
	  double del2 = 1.0/(st1_[kk1_+1] - st1_[2]);
	  double del3 = 1.0/(st1_[kn1_+kk1_-3] - st1_[kn1_-2]);
	  double del5 = 1.0/(st1_[kk1_] - st1_[2]);
	  double del6 = 1.0/(st1_[kn1_+kk1_-3] - st1_[kn1_-1]);
	  frac = (st1_[kk1_+1] - st1_[1])/(st1_[kn1_+kk1_-2] - st1_[kn1_-2] + 
					st1_[kk1_+1] - st1_[1]);
	  // 2. derivative.
	  for (kj=0; kj<kn2_; kj++)
	    {
	      if ((coefknown_[kj*kn1_+2] == 1 ||
		   coefknown_[kj*kn1_+2] == 2) && 
		  (coefknown_[kj*kn1_+kn1_-3] == 1 || 
		   coefknown_[kj*kn1_+kn1_-3] == 2))
		continue;

	      t1 = (coefknown_[kj*kn1_+2] == 1) ? 1.0 :
		  ((coefknown_[kj*kn1_+kn1_-3] == 1) ? 0.0 : frac);
	      t2 = 1.0 - t1;

	      for (kk=0; kk<idim_; kk++)
		{
		  tmean = t1*((scoef_[(kj*kn1_+2)*kdim_+kk] - 
			       scoef_[(kj*kn1_+1)*kdim_+kk])*del2 -
			      (scoef_[(kj*kn1_+1)*kdim_+kk] -
			       scoef_[kj*kn1_*kdim_+kk])*del1)*del5 +
		      t2*((scoef_[(kj*kn1_+kn1_-1)*kdim_+kk] - 
			   scoef_[(kj*kn1_+kn1_-2)*kdim_+kk])*del4 -
			  (scoef_[(kj*kn1_+kn1_-2)*kdim_+kk] - 
			   scoef_[(kj*kn1_+kn1_-3)*kdim_+kk])*del3)*del6;
		  td = (tmean/del5 + del1*(scoef_[(kj*kn1_+1)*kdim_+kk] -
					  scoef_[kj*kn1_*kdim_+kk]) +
		       del2*scoef_[(kj*kn1_+1)*kdim_+kk])/del2;
		  scoef_[(kj*kn1_+2)*kdim_+kk] = td;
		  td = (tmean/del6 - del4*(scoef_[(kj*kn1_+kn1_-1)*kdim_+kk] -
					  scoef_[(kj*kn1_+kn1_-2)*kdim_+kk]) +
		       del3*scoef_[(kj*kn1_+kn1_-2)*kdim_+kk])/del3;
		  scoef_[(kj*kn1_+kn1_-3)*kdim_+kk] = td;
		}
	    }
	}
    }

  if (cont_seem_[1] > 1)
    {
      double del1 = 1.0/(st2_[kk2_] - st2_[1]);
      double del4 = 1.0/(st2_[kn2_+kk2_-2] - st2_[kn2_-1]);
      frac = (st2_[kk2_] - st2_[1])/(st2_[kn2_+kk2_-2] - st2_[kn2_-1] +
				  st2_[kk2_] - st2_[1]);

      for (ki=0; ki<kn1_; ki++)
       {
	  if ((coefknown_[kn1_+ki] == 1 || coefknown_[kn1_+ki] == 2) && 
	      (coefknown_[(kn2_-2)*kn1_+ki] == 1 || 
	       coefknown_[(kn2_-2)*kn1_+ki] == 2))
	    continue;

	 t1 = (coefknown_[kn1_+ki] == 1) ? 1.0 :
	      ((coefknown_[(kn2_-2)*kn1_+ki] == 1) ? 0.0 : frac);
	 t2 = 1.0 - t1;

	 for (kk=0; kk<idim_; kk++)
	   {
	     tmean =
		 t1*(scoef_[(kn1_+ki)*kdim_+kk]-scoef_[ki*kdim_+kk])*del1 +
		 t2*(scoef_[((kn2_-1)*kn1_+ki)*kdim_+kk] - 
		     scoef_[((kn2_-2)*kn1_+ki)*kdim_+kk])*del4;
	     scoef_[(kn1_+ki)*kdim_+kk] =
		 (tmean+del1*scoef_[ki*kdim_+kk])/del1;
	     scoef_[((kn2_-2)*kn1_+ki)*kdim_+kk] = 
		 (del4*scoef_[((kn2_-1)*kn1_+ki)*kdim_+kk]-tmean)/del4;
	   }
       }

      if (cont_seem_[1] > 2)
	{
	  double del2 = 1.0/(st2_[kk2_+1] - st2_[2]);
	  double del3 = 1.0/(st2_[kn2_+kk2_-3] - st2_[kn2_-2]);
	  double del5 = 1.0/(st2_[kk2_] - st2_[2]);
	  double del6 = 1.0/(st2_[kn2_+kk2_-3] - st2_[kn2_-1]);
	  frac = (st2_[kk2_+1] - st2_[1])/(st2_[kn2_+kk2_-2] - st2_[kn2_-2] +
					st2_[kk2_+1] - st2_[1]);

	  for (ki=0; ki<kn1_; ki++)
	    {
	      if ((coefknown_[2*kn1_+ki] == 1 ||
		   coefknown_[2*kn1_+ki] == 2) && 
		  (coefknown_[(kn2_-3)*kn1_+ki] == 1 || 
		   coefknown_[(kn2_-3)*kn1_+ki] == 2))
		continue;

	      t1 = (coefknown_[2*kn1_+ki] == 1) ? 1.0 :
		  ((coefknown_[(kn2_-3)*kn1_+ki] == 1) ? 0.0 : frac);
	      t2 = 1.0 - t1;

	      for (kk=0; kk<idim_; kk++)
		{
		  tmean = t1*((scoef_[(2*kn1_+ki)*kdim_+kk] -
			       scoef_[(kn1_+ki)*kdim_+kk])*del2 -
			      (scoef_[(kn1_+ki)*kdim_+kk] -
			       scoef_[ki*kdim_+kk])*del1)*del5 +
		      t2*((scoef_[((kn2_-1)*kn1_+ki)*kdim_+kk] - 
			   scoef_[((kn2_-2)*kn1_+ki)*kdim_+kk])*del4 -
			  (scoef_[((kn2_-2)*kn1_+ki)*kdim_+kk] -
			   scoef_[((kn2_-3)*kn1_+ki)*kdim_+kk])*del3)*del6;
		  scoef_[(2*kn1_+ki)*kdim_+kk] = 
		      (tmean/del5 + del1*(scoef_[(kn1_+ki)*kdim_+kk] -
					  scoef_[ki*kdim_+kk]) +
		       del2*scoef_[(kn1_+ki)*kdim_+kk])/del2;
		  scoef_[((kn2_-3)*kn1_+ki)*kdim_+kk] = 
		      (tmean/del6 - del4*
		       (scoef_[((kn2_-1)*kn1_+ki)*kdim_+kk] -
			scoef_[((kn2_-2)*kn1_+ki)*kdim_+kk]) +
		       del3*scoef_[((kn2_-2)*kn1_+ki)*kdim_+kk])/del3;
		}
	    }

	}
    }
  return 0;
}



/****************************************************************************/

void
SmoothSurf::setC1AtSeem(int pardir,  // Parameter direction, 1 or 2.
			double weight)
//--------------------------------------------------------------------------
//     Purpose : Set C1-continuity constraints in one par. dir.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF Oslo,  04.98.
//--------------------------------------------------------------------------
{
  int ki, kj, kk, kp, kq;
  int sign = 1;
  
  int kl1, kl2;

  double del[2];
  if (pardir == 1)
    {
      del[0] = 1.0/(st1_[kk1_] - st1_[1]);
      del[1] = 1.0/(st1_[kn1_+kk1_-2] - st1_[kn1_-1]);
    }
  else
    {
      del[0] = 1.0/(st2_[kk2_] - st2_[1]);
      del[1] = 1.0/(st2_[kn2_+kk2_-2] - st2_[kn2_-1]);
    }

  double tdel1, tdel2;
  double tintgr;
  for (ki=0; ki<kn1_; ki++)
   {
      if (pardir == 2 && coefknown_[kn1_+ki] > 0 && 
 	 coefknown_[(kn2_-2)*kn1_+ki] > 0) 
        continue;

     if (pardir == 1 && (ki > 1 && ki < kn1_-2))
       continue;

     for (kj=0; kj<kn2_; kj++)
      {
 	if (pardir == 1 && coefknown_[kj*kn1_+1] > 0 && 
 	    coefknown_[kj*kn1_+kn1_-2] > 0) 
 	  continue;

	if (coefknown_[kj*kn1_+ki] == 2)
	  continue;

	if (pardir == 2 && (kj > 1 && kj < kn2_-2))
	  continue;

	kl1 = (coefknown_[kj*kn1_+ki] > 2) ?
	    pivot_[coefknown_[kj*kn1_+ki]-kpointer_] : pivot_[kj*kn1_+ki];

	for (kp=0; kp<kn1_; kp++)
	 {
 	   if (pardir == 2 && coefknown_[kn1_+kp] > 0 && 
 	       coefknown_[(kn2_-2)*kn1_+kp] > 0) 
 	     continue;

	   if (pardir == 1 && (kp > 1 && kp < kn1_-2))
	     continue;

	   for (kq=0; kq<kn2_; kq++)
	    {
 	      if (pardir == 1 && coefknown_[kq*kn1_+1] > 0 && 
 		  coefknown_[kq*kn1_+kn1_-2] > 0) 
 		continue;

	      if (coefknown_[kq*kn1_+kp] == 1 || coefknown_[kq*kn1_+kp] == 2)
		continue;

	      if (pardir == 2 && (kq > 1 && kq < kn2_-2))
		continue;

	      kl2 = (coefknown_[kq*kn1_+kp] > 2) ?
		  pivot_[coefknown_[kq*kn1_+kp]-kpointer_] :
		  pivot_[kq*kn1_+kp];

	      if (pardir == 1)
		{
		  sign = (ki == kp || ki+kp == kn1_-1) ? 1 : -1;
		  tdel1 = (ki <= 1) ? del[0] : del[1];
		  tdel2 = (kp <= 1) ? del[0] : del[1];
		  tintgr = integral2_[0][kj][kq];
		}
	      else
		{
		  sign = (kj == kq || kj+kq == kn2_-1) ? 1 : -1;
		  tdel1 = (kj <= 1) ? del[0] : del[1];
		  tdel2 = (kq <= 1) ? del[0] : del[1];
// 		  tdel2 = (kp <= 1) ? del[0] : del[1];
		  tintgr = integral1_[0][ki][kp];
		}

	      if (coefknown_[kj*kn1_+ki] == 1)
		{
		  for (kk=0; kk<idim_; kk++)
		    gright_[kk*kncond_+kl2] -= 
		      sign*weight*tdel1*tdel2*tintgr*
			scoef_[(kj*kn1_+ki)*kdim_+kk];
		}
	      else
		{
		  // if (kl2 > kl1) continue;

		  for (kk=0; kk<norm_dim_; kk++)
		    {
		      gmat_[(kk*kncond_+kl2)*norm_dim_*kncond_+kk*kncond_+kl1] += 
			sign*weight*tdel1*tdel2*tintgr;
		      // if (kl2 < kl1)
		    // gmat_[(kk*kncond_+kl2)*norm_dim_*kncond_+kk*kncond_+kl1] += 
			  // sign*weight;
		    }
		}
	    }
	 }
      }
   }
	
  return;
}

/****************************************************************************/

void
SmoothSurf::setC2AtSeem(int pardir,  // Parameter direction, 1 or 2.
			  double weight)
//--------------------------------------------------------------------------
//     Purpose : Set C2-continuity constraints in one par. dir.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF Oslo,  08.99.
//--------------------------------------------------------------------------
{
  int ki, kj, kk, kp, kq;
  int k1, k2;
  int sign1, sign2;
  
  int kl1, kl2;

  double del1[4], del2[2], dx[6];
  double tintgr;
  if (pardir == 1)
    {
      del1[0] = 1.0/(st1_[kk1_] - st1_[1]);
      del1[1] = 1.0/(st1_[kk1_+1] - st1_[2]);
      del1[2] = 1.0/(st1_[kn1_+kk1_-3] - st1_[kn1_-2]);
      del1[3] = 1.0/(st1_[kn1_+kk1_-2] - st1_[kn1_-1]);
      del2[0] = 1.0/(st1_[kk1_] - st1_[2]);
      del2[1] = 1.0/(st1_[kn1_+kk1_-3] - st1_[kn1_-1]);
    }
  else
    {
      del1[0] = 1.0/(st2_[kk2_] - st2_[1]);
      del1[1] = 1.0/(st2_[kk2_+1] - st2_[2]);
      del1[2] = 1.0/(st2_[kn2_+kk2_-3] - st2_[kn2_-2]);
      del1[3] = 1.0/(st2_[kn2_+kk2_-2] - st2_[kn2_-1]);
      del2[0] = 1.0/(st2_[kk2_] - st2_[2]);
      del2[1] = 1.0/(st2_[kn2_+kk2_-3] - st2_[kn2_-1]);
    }

  dx[0] = del2[0]*del1[0];
  dx[1] = del2[0]*(del1[0] + del1[1]);
  dx[2] = del2[0]*del1[1];
  dx[3] = del2[1]*del1[2];
  dx[4] = del2[1]*(del1[2]+del1[3]);
  dx[5] = del2[1]*del1[3];

  for (sign2=1, k2=0, kp=0; kp<kn1_; kp++)
    {
      if (pardir == 2 && coefknown_[kn1_+kp] > 0 && 
	  coefknown_[(kn2_-2)*kn1_+kp] > 0) 
	continue;

      if (pardir == 1 && kp > 2 && kp < kn1_-3)
	continue;

      for (kq=0; kq<kn2_; kq++)
	{
	  if (pardir == 1 && coefknown_[kq*kn1_+1] > 0 && 
	      coefknown_[kq*kn1_+kn1_-2] > 0) 
	    continue;

	  if (pardir == 2 && kq > 2 && kq < kn2_-3)
	    continue;

	  kl2 = (coefknown_[kq*kn1_+kp] > 2) ?
	      pivot_[coefknown_[kq*kn1_+kp]-kpointer_] : pivot_[kq*kn1_+kp];

	  for (sign1=1, k1=0, ki=0; ki<kn1_; ki++)
	    {
	      if (pardir == 2 && coefknown_[kn1_+ki] > 0 && 
		  coefknown_[(kn2_-2)*kn1_+ki] > 0) 
		continue;

	      if (pardir == 1 && ki > 2 && ki < kn1_-3)
		continue;

	      for (kj=0; kj<kn2_; kj++)
		{
		  if (pardir == 1 && coefknown_[kj*kn1_+1] > 0 && 
		      coefknown_[kj*kn1_+kn1_-2] > 0) 
		    continue;

		  if (pardir == 2 && kj > 2 && kj < kn2_-3)
		    continue;

		  kl1 = (coefknown_[kj*kn1_+ki] > 2) ?
		      pivot_[coefknown_[kj*kn1_+ki]-kpointer_] :
		      pivot_[kj*kn1_+ki];

		  tintgr = (pardir == 1) ? integral2_[0][kj][kq] :
		      integral1_[0][ki][kp];

		  if (coefknown_[kq*kn1_+kp] != 1 && 
		      coefknown_[kq*kn1_+kp] != 2 &&
		      coefknown_[kj*kn1_+ki] != 2)
		    {
		      if (coefknown_[kj*kn1_+ki] == 1)
			{
			  for (kk=0; kk<idim_; kk++)
			    gright_[kk*kncond_+kl2] -= 
				weight*sign1*sign2*dx[k1]*dx[k2]*
				tintgr*scoef_[(kj*kn1_+ki)*kdim_+kk];
			}
		      else
			{
			  for (kk=0; kk<norm_dim_; kk++)
			    gmat_[(kk*kncond_+kl2)*norm_dim_*kncond_+
				  kk*kncond_+kl1] += 
				weight*sign1*sign2*dx[k1]*dx[k2]*tintgr;
			}
		    }
		  if (pardir == 2)
		    {
		      sign1 *= -1;
		      k1++;
		    }
		}
	      if (pardir == 1)
		{
		  sign1 *= -1;
		  k1++;
		}
	      else
		{
		  sign1 = 1;
		  k1 = 0;
		}
	    }
	  if (pardir == 2)
	    {
	      sign2 *= -1;
	      k2++;
	    }
	}
      if (pardir == 1)
	{
	  sign2 *= -1;
	  k2++;
	}
      else
	{
	  sign2 = 1;
	  k2 = 0;
	}
    }
}



/****************************************************************************/
void SmoothSurf::setRationalCnAtSeem(int pardir, int cn, double weight1, double weight2)

//--------------------------------------------------------------------------
//     Purpose : Set Cn-continuity constraints in one par. dir (1 or 2)
//               for the rational case
//--------------------------------------------------------------------------
{
  // For pardir=1, the constraint is given by minimizing the integral over all v
  // of (S_u(u_min,v)-S_u(u_max,v))^2, and (only in case cn=2) the integral
  // of (S_uu(u_min,v)-S_uu(u_max,v))^2 where S_u and S_uu is
  // first and second order differentiation of the surface w.r.t. u.
  // For pardir=2, the same is done by integration along u and differentiating
  // w.r.t. v
  //
  // Throughout the computations, we assume that the weights along the two edge
  // curves are identical.
  //
  // The integration is done by estimation using Gaussian quadrature

  // Get some direction dependent constants
  // int ncoefs_1 = (pardir==1) ? kn1_ : kn2_;   // Along direction of Cn-continuity
  int ncoefs_2 = (pardir==1) ? kn2_ : kn1_;   // Along integration direction
  int order_1 = (pardir==1) ? kk1_ : kk2_;
  int order_2 = (pardir==1) ? kk2_ : kk1_;
  BsplineBasis basis_1 = (pardir==1) ? bspline_surface_->basis_u() : bspline_surface_->basis_v();
  BsplineBasis basis_2 = (pardir==1) ? bspline_surface_->basis_v() : bspline_surface_->basis_u();

  // Create the denominator surface function, and storage of evaluations of denominator function
  vector<double> denom(cn*2+1);  // For pardir=1 and a given value v, holds the values
                                 // D(u_min,v), D_u(u_min,v), D_u(u_max,v), D_uu(u_min,v), D_uu(u_max,v),
                                 // (the last two only if cn=2) where D is the denominator surface function.
                                 // Remember that D(u_min,v) and D(u_max,v) are assumed to be identical,
                                 //
                                 // For pardir=2 and a given value u, holds the values
                                 // D(u,v_min), D_v(u,v_min), D_v(u,v_max), D_vv(u,v_min), D_vv(u,v_max)
  int denom_eval_size = ((cn+1)*(cn+2))/2;
  vector<Point> denom_eval(denom_eval_size);   // Used for evaluation to get the values for denom
  for (int i = 0; i < denom_eval_size; ++i)
    denom_eval[i] = Point(1);
  vector<double> denomcoefs(kn1_*kn2_);
  for (int i = 0, scoef_pos = idim_; i < kn1_*kn2_; ++i, scoef_pos+=kdim_)
    denomcoefs[i] = scoef_[scoef_pos];
  shared_ptr<SplineSurface> denom_surf(new SplineSurface(kn1_, kn2_, kk1_, kk2_, st1_, st2_, denomcoefs.begin(), 1, false));

  // Get parameter sample values to evaluate at
  vector<double> gq_points;    // All sample parameter values along the integration path
  vector<double> gq_weights;   // The weights used for numerical integration
  GaussQuadValues(basis_2, gq_points, gq_weights);
  int bez_segs = ncoefs_2 - order_2 + 1;   // Number of Bezier segments in integration direction
  int seg_samples = (int)gq_weights.size();   // Number of samples inside each Bezier segment in integration direction
  // int samples = bez_segs * seg_samples;  // Number of samples in total ( == gq_points.size() ) in integration direction

  double ta1 = st1_[kk1_-1];   /* Start of par. interval in 1. par. dir. */
  double ta2 = st1_[kn1_];     /* End of par. interval in 1. par. dir.   */
  double tb1 = st2_[kk2_-1];   /* Start of par. interval in 2. par. dir. */
  double tb2 = st2_[kn2_];     /* End of par. interval in 2. par. dir.   */

  // Coefficient terms in expression to square and minimize
  vector<double> term_coefs_c1(order_2*4);  // For each sample value, the terms for the coefficients in the
                             // expression to minimize for C1-continuity. For pardir=1, the terms
                             // come as follows: First all c_(0,j) where j goes through all
                             // Bsplines in integration direction with support in the sample value.
                             // Then all c_(n-1,j) (the coefficients at the right edge of the surface)
                             // Then all c_(1,j) and then all c_(n-2,j).
                             // For pardir=2, first all c_(j,0) then all c_(j,n-1) etc.
  vector<double> term_coefs_c2;   // Same as term_coefs_c1, but for C2-continuity. Additionally holds terms for the
                   // c_(2,j) and c_(n-3,j) for pardir=1, or c_(j,2) and c_(j,n-3) for pardir=2
  if (cn==2)
    term_coefs_c2.resize(order_2*6);

  double
    d_left_0 = 0.0, d_right_0 = 0.0,
    d_left_1 = 0.0, d_right_1 = 0.0;
  if (pardir==1)
    {
      d_left_0 = 1.0 / (st1_[kk1_] - ta1);
      d_right_0 = 1.0 / (ta2 - st1_[kn1_-1]);
      if (cn==2)
	{
	  d_left_1 = 1.0 / (st1_[kk1_+1] - ta1);
	  d_right_1 = 1.0 / (ta2 - st1_[kn1_-2]);
	}
    }
  else
    {
      d_left_0 = 1.0 / (st2_[kk2_] - tb1);
      d_right_0 = 1.0 / (tb2 - st2_[kn2_-1]);
      if (cn==2)
	{
	  d_left_1 = 1.0 / (st2_[kk2_+1] - tb1);
	  d_right_1 = 1.0 / (tb2 - st2_[kn2_-2]);
	}
    }

  // Run through all sample points and add terms to equation system matrices
  vector<double>::const_iterator par_it = gq_points.begin();
  for (int seg = 0; seg < bez_segs; ++seg)      // For every Bezier segment
    {
      double length = (pardir==1)
	? st2_[seg+kk2_]-st2_[seg+kk2_-1]
	: st1_[seg+kk1_]-st1_[seg+kk1_-1];
      for (int sample = 0; sample < seg_samples; ++sample, ++par_it)    // For every sample point in given segment
	{
	  double par = *par_it;

	  // Get denominator evalutions, including appropriate differentiations
	  if (pardir==1)
	    {
	      denom_surf->point(denom_eval, ta1, par, cn, true, true);
	      denom[0] = denom_eval[0][0];
	      denom[1] = denom_eval[1][0];
	      if (cn==2)
		denom[3] = denom_eval[3][0];

	      denom_surf->point(denom_eval, ta2, par, cn, false, true);
	      denom[2] = denom_eval[1][0];
	      if (cn==2)
		denom[4] = denom_eval[3][0];
	    }
	  else
	    {
	      denom_surf->point(denom_eval, par, tb1, cn, true, true);
	      denom[0] = denom_eval[0][0];
	      denom[1] = denom_eval[2][0];
	      if (cn==2)
		denom[3] = denom_eval[5][0];

	      denom_surf->point(denom_eval, par, tb2, cn, true, false);
	      denom[2] = denom_eval[2][0];
	      if (cn==2)
		denom[4] = denom_eval[5][0];
	    }

	  vector<double> basis_vals = basis_2.computeBasisValues(par);  // The B-spline evaluations along integration direction
	  int basis_vals_left = basis_2.lastKnotInterval();

	  // Calculate C1-continuity expression terms
	  double left_0 = (-denom[1]-(order_1-1)*d_left_0*denom[0]) / (denom[0]*denom[0]);
	  double right_0 = (denom[2]-(order_1-1)*d_right_0*denom[0]) / (denom[0]*denom[0]);
	  double left_1 = (order_1-1)*d_left_0/denom[0];
	  double right_1 = (order_1-1)*d_right_0/denom[0];

	  for (int i = 0; i < order_2; ++i)
	    {
	      term_coefs_c1[i] = left_0 * basis_vals[i];
	      term_coefs_c1[i+order_2] = right_0 * basis_vals[i];
	      term_coefs_c1[i+2*order_2] = left_1 * basis_vals[i];
	      term_coefs_c1[i+3*order_2] = right_1 * basis_vals[i];
	    }

	  // Calculate C2-continuity expression terms
	  if (cn==2)
	    {
	      left_0 = ((order_1-1) * (order_1-2) * d_left_0 * d_left_0 * denom[0] * denom[0]
			+ 2.0 * denom[1] * denom[1]
			+ 2.0 * (order_1-1) * d_left_0 * denom[1] * denom[0]
			- denom[3] * denom[0]) / (denom[0]*denom[0]*denom[0]);
	      right_0 = (- (order_1-1) * (order_1-2) * d_right_0 * d_right_0 * denom[0] * denom[0]
			 - 2.0 * denom[2] * denom[2]
			 + 2.0 * (order_1-1) * d_right_0 * denom[2] * denom[0]
			 - denom[4] * denom[0]) / (denom[0]*denom[0]*denom[0]);
	      left_1 = (order_1-1) * d_left_0 * (- (order_1-2) * (d_left_0 + d_left_1) * denom[0]
						 - 2.0 * denom[1]) / (denom[0] * denom[0]);
	      right_1 = (order_1-1) * d_right_0 * ((order_1-2) * (d_right_0 + d_right_1) * denom[0]
						   - 2.0 * denom[2]) / (denom[0] * denom[0]);
	      double left_2 = (order_1-1) * (order_1-2) * d_left_0 * d_left_1 / denom[0];
	      double right_2 = - (order_1-1) * (order_1-2) * d_right_0 * d_right_1 / denom[0];

	      for (int i = 0; i < order_2; ++i)
		{
		  term_coefs_c2[i] = left_0 * basis_vals[i];
		  term_coefs_c2[i+order_2] = right_0 * basis_vals[i];
		  term_coefs_c2[i+2*order_2] = left_1 * basis_vals[i];
		  term_coefs_c2[i+3*order_2] = right_1 * basis_vals[i];
		  term_coefs_c2[i+4*order_2] = left_2 * basis_vals[i];
		  term_coefs_c2[i+5*order_2] = right_2 * basis_vals[i];
		}
	    }

	  // Combine to get terms to equation system

	  // Loop for each position in continuity direction for the first point.
	  // 0 = first position, 1 = last, 2 = next-to-first, 3 = next-to-last, etc.
	  for (int cont_dir_pos_1 = 0; cont_dir_pos_1 < 2*(cn+1); ++cont_dir_pos_1)
	    for (int int_dir_pos_1 = 0; int_dir_pos_1 < order_2; ++int_dir_pos_1)   // For each position in integration direction for first point
	      {
		int i_pos_1, j_pos_1;
		if (pardir==1)
		  {
		    i_pos_1 = ((cont_dir_pos_1 & 1) == 0)
		      ? (cont_dir_pos_1 >> 1)
		      : kn1_ - 1 - (cont_dir_pos_1 >> 1);
		    j_pos_1 = basis_vals_left - kk2_ + 1 + int_dir_pos_1;
		  }
		else
		  {
		    i_pos_1 = basis_vals_left - kk1_ + 1 + int_dir_pos_1;
		    j_pos_1 = ((cont_dir_pos_1 & 1) == 0)
		      ? (cont_dir_pos_1 >> 1)
		      : kn2_ - 1 - (cont_dir_pos_1 >> 1);
		  }

		int pos_1 = j_pos_1*kn1_ + i_pos_1;
		if (coefknown_[pos_1] == 1 || coefknown_[pos_1] == 2)
		  continue;
		int piv_1 = (coefknown_[pos_1] > 2) ?
		  pivot_[coefknown_[pos_1]-kpointer_] : pivot_[pos_1];

		// Loop for each position in continuity direction for the second point.
		for (int cont_dir_pos_2 = 0; cont_dir_pos_2 < 2*(cn+1); ++cont_dir_pos_2)
		  for (int int_dir_pos_2 = 0; int_dir_pos_2 < order_2; ++int_dir_pos_2)   // For each position in integration direction for second point
		    {
		      int i_pos_2, j_pos_2;
		      if (pardir==1)
			{
			  i_pos_2 = ((cont_dir_pos_2 & 1) == 0)
			    ? (cont_dir_pos_2 >> 1)
			    : kn1_ - 1 - (cont_dir_pos_2 >> 1);
			  j_pos_2 = basis_vals_left - kk2_ + 1 + int_dir_pos_2;
			}
		      else
			{
			  i_pos_2 = basis_vals_left - kk1_ + 1 + int_dir_pos_2;
			  j_pos_2 = ((cont_dir_pos_2 & 1) == 0)
			    ? (cont_dir_pos_2 >> 1)
			    : kn2_ - 1 - (cont_dir_pos_2 >> 1);
			}

		      int pos_2 = j_pos_2*kn1_ + i_pos_2;
		      int piv_2 = (coefknown_[pos_2] > 2) ?
			pivot_[coefknown_[pos_2]-kpointer_] : pivot_[pos_2];

		      if (piv_2 < piv_1 && coefknown_[pos_2] != 1)
			continue;

		      double term = 0.0;
		      if (cont_dir_pos_1 < 4 && cont_dir_pos_2 < 4)
			term = weight1
			  * term_coefs_c1[cont_dir_pos_1*order_2+int_dir_pos_1]
			  * term_coefs_c1[cont_dir_pos_2*order_2+int_dir_pos_2];
		      if (cn==2)
			term += weight2
			  * term_coefs_c2[cont_dir_pos_1*order_2+int_dir_pos_1]
			  * term_coefs_c2[cont_dir_pos_2*order_2+int_dir_pos_2];

		      term *= gq_weights[sample] * length;

		      // Add contribution to equation system
		      if (coefknown_[pos_2] == 1)
			{
			  // The contribution of this term is added to the right
			  // side of the equation system.

			  double* sc = &*scoef_ + pos_2*kdim_;
			  for (int k=0; k<idim_; ++k)
			    gright_[k*kncond_+piv_1] -= sc[k]*term;
			}
		      else
			{
			  // The contribution of this term is added to the left
			  //  side of the equation system.
			  for (int k=0; k<norm_dim_; k++)
			    {
			      gmat_[(k*kncond_+piv_2)*norm_dim_*kncond_+k*kncond_+piv_1]
				+= term;
			      if (piv_1 != piv_2)
				gmat_[(k*kncond_+piv_1)*norm_dim_*kncond_+k*kncond_+piv_2]
				  += term;
			    }
			}
		    }    // End -- For each second sample point
	      }    // End -- For each first sample point
	}    // End -- For every sample point in given segment
    }   // End -- For every Bezier segment

}



/****************************************************************************/

void
SmoothSurf::preparePeriodicity(int seem[]) // Information about conditions
				       // along the seem of closed surfaces.
//--------------------------------------------------------------------------
//     Purpose : Set pointers between identical coefficients at a periodic 
//               seem. If possible, update fixed coefficients at the seem.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF Oslo,  07.99.
//--------------------------------------------------------------------------
{
  int ki, kj;
  int idx[4];
  int kcfix = 0;
  double corner[4];
  corner[0] = corner[1] = corner[2] = corner[3] = 0.0;
  idx[0] = 0;
  idx[1] = kn1_-1;
  idx[2] = (kn2_-1)*kn1_;
  idx[3] = kn2_*kn1_-1;

  if (seem[0] > 0 && seem[1] > 0)
    {
      // Treat the corners of a double closed surface specifically.
      // First count the number of fixed corners.

      for (ki=0; ki<4; ki++)
	{
	  if (coefknown_[idx[ki]])
	    {
	      kcfix++;

	      for (kj=0; kj<kdim_; kj++)
	        corner[kj] += scoef_[idx[ki]*kdim_+kj];
	    }
	}

      if (kcfix == 0)
	{
	  // No fixed coefficients. Set pointers to representative
	  // coefficient.

	  for (ki=1; ki<4; ki++)
	    coefknown_[idx[ki]] = kpointer_ + idx[0];
	}
      else if (kcfix < 4)
	{
	  // Fix the corner coefficients corresponding to the fixed ones.

	  for (kj=0; kj<kdim_; kj++)
	    corner[kj] /= (double)kcfix;

	  for (ki=0; ki<4; ki++)
	    {
	      if (coefknown_[idx[ki]] == 0)
		{
		  coefknown_[idx[ki]] = 1;
		  for (kj=0; kj<kdim_; kj++)
		    scoef_[idx[ki]*kdim_+kj] = corner[kj];
		}
	    }
	}
    }

  // Treat the periodic edges.
  // 1. parameter direction.

  int kstart = (seem[0] > 0 && seem[1] > 0) ? 1 : 0;
  int kend = (seem[0] > 0 && seem[1] > 0) ? kn1_-1 : kn1_;

  if (seem[1] > 0)
    {
      for (ki=kstart; ki<kend; ki++)
	{
	  if (coefknown_[ki] > 0 && coefknown_[(kn2_-1)*kn1_+ki] > 0)
	    continue;

	  if (coefknown_[ki] == 1)
	    {
	      coefknown_[(kn2_-1)*kn1_+ki] = 1;
	      for (kj=0; kj<kdim_; kj++)
		scoef_[((kn2_-1)*kn1_+ki)*kdim_+kj] = scoef_[ki*kdim_+kj];
	    }
	  else if (coefknown_[(kn2_-1)*kn1_+ki] == 1)
	    {
	      coefknown_[ki] = 1;
	      for (kj=0; kj<kdim_; kj++)
		scoef_[ki*kdim_+kj] = scoef_[((kn2_-1)*kn1_+ki)*kdim_+kj];
	    }
	  else
	    coefknown_[(kn2_-1)*kn1_+ki] = kpointer_ + ki;
	}
    }

  // 2. parameter direction.

  kstart = (seem[0] > 0 && seem[1] > 0) ? 1 : 0;
  kend = (seem[0] > 0 && seem[1] > 0) ? kn2_-1 : kn2_;

  if (seem[0] > 0)
    {
      for (ki=kstart; ki<kend; ki++)
	{
	  if (coefknown_[ki*kn1_] > 0 && coefknown_[ki*kn1_+kn1_-1] > 0)
	    continue;

	  if (coefknown_[ki*kn1_] == 1)
	    {
	      coefknown_[ki*kn1_+kn1_-1] = 1;
	      for (kj=0; kj<kdim_; kj++)
		scoef_[(ki*kn1_+kn1_-1)*kdim_+kj] =
		    scoef_[(ki*kn1_)*kdim_+kj];
	    }
	  else if (coefknown_[ki*kn1_+kn1_-1] == 1)
	    {
	      coefknown_[ki*kn1_] = 1;
	      for (kj=0; kj<kdim_; kj++)
		scoef_[(ki*kn1_)*kdim_+kj] =
		    scoef_[(ki*kn1_+kn1_-1)*kdim_+kj];
	    }
	  else
	    coefknown_[ki*kn1_+kn1_-1] = kpointer_ + ki*kn1_;
	}
    }
}

