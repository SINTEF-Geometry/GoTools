// ----------------------------------------------------------------
//       Implementation file for class SmoothSurfSet.
// ----------------------------------------------------------------
//
// This file contains the following member functions:
//             a. SmoothSurfSet()
//             b. ~SmoothSurfSet()
//             d. attach()
//             e. setOptimize()
//             f. getBasis()
//             g. setLeastSquares()
//             h. setNormalCond()
//             i. approxOrig()
//             j. 
//             k. equationSolve()
// ----------------------------------------------------------------


#include "GoTools/creators/SmoothSurfSet.h"
#include "GoTools/creators/Integrate.h"
#include "GoTools/creators/SolveCG.h"
//#include "newmat.h"
#include "GoTools/utils/LUDecomp.h"
#include "GoTools/utils/Values.h"

#include <cmath>
#include <algorithm>


using namespace Go;
using std::vector;
using std::max;
using std::min;

SmoothSurfSet::SmoothSurfSet()
    : copy_coefs_(true)
   //--------------------------------------------------------------------------
   //     Constructor for class SmoothSurfSet.
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

   ider_ = 3;
   idim1_ = idim_ = 3;
   kdim_ = 1;
   kncond_ = 0;
   knconstraint_ = 0;
}

/****************************************************************************/


SmoothSurfSet::SmoothSurfSet(bool copy_coefs)
    : copy_coefs_(copy_coefs)
   //--------------------------------------------------------------------------
   //     Constructor for class SmoothSurfSet.
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

   ider_ = 3;
   idim1_ = idim_ = 3;
   kdim_ = 1;
   kncond_ = 0;
   knconstraint_ = 0;
}

/****************************************************************************/


SmoothSurfSet::~SmoothSurfSet()
   //--------------------------------------------------------------------------
   //     Destructor for class SmoothSurfSet.
   //
   //     Written by : Vibeke Skytt,  SINTEF SI,  09.93.
   //--------------------------------------------------------------------------
{
   if (copy_coefs_ && coef_array_.size() > 0) 
     coef_array_.erase(coef_array_.begin(), coef_array_.end());

}


   
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  SmoothSurfSet::attach
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SmoothSurfSet::attach(std::vector<shared_ptr<SplineSurface> >& insf,  
		      // Input surface representing the spline space
		      std::vector< std::vector<int> >& coef_known, 
			// Indicates whether or not each coefficient is 
			// known (1), unknown (0) or not of interest(2)
			int num_side_constraints, 
			// Number of side constraints to the minimization 
			// problem
			int has_normal_cond)  
                        // Indicates if normal conditions will be given.

//--------------------------------------------------------------------------
//
//     Purpose : Initialize class variables related to the spline space.
//
//     Calls   :
//
// NB!
// If attach is used more than once for an instance of SmoothSurfSet, 
// then it is assumed that the new surface (insf) lies in the same
// spline space that the old surface (srf) and that there
// is no change concerning the existance of normal conditions.
// If these conditions are broken, memory faults and/or wrong
// results will occur. Including higher order derivatives in the
// smoothing functional than in previous runs, will have no effect.
//
//     Written by : Vibeke Skytt,  SINTEF SI,  09.95.
//     Revised by : Vibeke Skytt,  SINTEF,  11.99.
//--------------------------------------------------------------------------
{
   int kh, ki, kj;
   int nmbsfs = (int)insf.size();
   vector<double>::iterator  scoef; // Pointer to surface coefficients. 

   // Dimension of geometry space
   idim_ = insf[0]->dimension();
   
   if (has_normal_cond && idim_ == 3)
       kdim_ = 3;          // Prepare for approximation of normal conditions.

   // Prepare for data storage for each surface.
   if (int(srfs_.size()) != nmbsfs)
     {
       srfs_.resize(nmbsfs);
       coef_array_.resize(nmbsfs);
       pivot_.resize(nmbsfs);
       surf_integral_.resize(nmbsfs);
     }

   kncond_ = 0;  // No free coefficients found yet.
   coef_known_.resize(coef_known.size());
   for (kh=0; kh<nmbsfs; kh++)
     {
       // Fetch data from B-spline surface. 
       bool israt = insf[kh]->rational();
       int dim = (israt) ? idim_+1 : idim_;
       idim1_ = dim;
       int kn1 = insf[kh]->numCoefs_u();
       int kn2 = insf[kh]->numCoefs_v();

       if (copy_coefs_)
	 {
	   coef_array_[kh].resize(kn1*kn2*dim);
	   vector<double>::const_iterator sc = israt ? insf[kh]->rcoefs_begin()
	     : insf[kh]->coefs_begin();
	   std::copy(sc, sc+kn1*kn2*dim, coef_array_[kh].begin());
	   scoef = coef_array_[kh].begin();
	 }
       else
	 scoef = (israt) ? insf[kh]->rcoefs_begin() : insf[kh]->coefs_begin();

       if (israt)
	 {
	   // Divide by the weight.

	   double *sc = &scoef[0];
	   for (ki=0; ki<kn1*kn2; ki++, sc+=idim1_)
	     for (kj=0; kj<idim_; kj++)
	       sc[kj] /= sc[idim_];
	 }

       if (srfs_[kh].get() != 0 &&  
	   insf[kh]->numCoefs_u() == srfs_[kh]->numCoefs_u() && 
	   insf[kh]->numCoefs_v() == srfs_[kh]->numCoefs_v() && 
	   insf[kh]->order_u() == srfs_[kh]->order_u() && 
	   insf[kh]->order_v() == srfs_[kh]->order_v() &&
	   num_side_constraints == knconstraint_)
	 {
	   // It is assumed that the new surface (insf[kh]) lies in the same
	   // spline space that the old surface and that there
	   // is no change concerning the existance of normal conditions.
	   // If these conditions are broken, memory faults and/or wrong
	   // results will occur. Including higher order derivatives in the
	   // smoothing functional than in previous runs, will have no effect.

	   surf_integral_[kh].integralset = true;
	   srfs_[kh] = insf[kh];
	 }
       else
	 {
	   if (srfs_[kh].get() != 0)
	     {
	       // Memory already allocated for the previous iteration of
	       // surface editing and smoothing. Free this memory.

	       surf_integral_[kh].erase();
	     }

	   surf_integral_[kh].resize(ider_, kn1, kn2);
	   srfs_[kh] = insf[kh];

	   // Allocate scratch for pivot_ array.
	   pivot_[kh].resize(kn1*kn2);

	 }

       // Store information about free, fixed and irrelevant coefficients
       // and create pivot_ array. Count number of free coefficients.
       coef_known_[kh] = coef_known[kh].begin();
       std::fill(pivot_[kh].begin(), pivot_[kh].end(), 0);
       int kn12 = kn1*kn2;
       for (kj=0; kj<kn12; kj+=kn1)
	 for (ki=0; ki<kn1; ki++)
	   {
	     if (coef_known_[kh][kj+ki] == 0)
	       {
		 pivot_[kh][kj+ki] = kncond_;
		 kncond_++;
	       }
	     else
	       pivot_[kh][kj+ki] = -MAXINT;
	   }
     }

   // Add side constraints
   knconstraint_ = num_side_constraints;
   kncond_ += knconstraint_;

   // Allocate scratch for arrays in the equation system. 

   gmat_.resize(kdim_*kdim_*kncond_*kncond_);
   gright_.resize(idim_*kncond_);
   std::fill(gmat_.begin(), gmat_.end(), 0.0);
   std::fill(gright_.begin(), gright_.end(), 0.0);

   return;
}



/***************************************************************************/

void
SmoothSurfSet::setOptimize(double wgt1,  /* Weight of 1. order term */
			     double wgt2,  /* Weight of 2. order term */
			     double wgt3)  /* Weight of 3. order term */
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		    from the smoothness functional.
//
//
//     Written by : Vibeke Skytt,  SINTEF SI,  09-10.93.
//--------------------------------------------------------------------------
{
    //int grstat = 0;   // Initialize status variable.
   int idxsf, ki, kj, kk, kp, kq, kr;
   int k1, k2, k3, k4;
   int kjstart, kjend;
   int kistart, kiend;

   int kleft = 0;    /* Parameter used in 1220 to be positioned in
			the knot vector.                            */
   int kstart;
   int kl1, kl2;
   double const1 = wgt1*M_PI;
   double const2 = wgt2*M_PI/(double)8.0;
   double const3 = wgt3*M_PI/(double)32.0;
   int numsfs = (int)srfs_.size();  // Number of surfaces in the surface set.
   vector<double>::iterator  scoef; // Pointer to surface coefficients. 
   //Integrate GaussQuad;     // Class performing integration of inner
                              // products of B-splines.

   // Modify number of derivatives to compute.
   int der = ider_;
   if (der == 3 && wgt3 == 0.0)
     der--;
   if (der == 2 && wgt2 == 0.0)
     der--;
   if (der == 1 && wgt1 == 0.0)
     der--;
   ider_ = der;

   // For each surface add the contribution of smoothing to the equation 
   // system.
   
   for (idxsf=0; idxsf<numsfs; idxsf++)
     {
       // Set parameter area.

       double ta1 = srfs_[idxsf]->startparam_u(); 
       double ta2 = srfs_[idxsf]->endparam_u(); 
       double tb1 = srfs_[idxsf]->startparam_v();
       double tb2 = srfs_[idxsf]->endparam_v(); 
       
       int kn1 = srfs_[idxsf]->numCoefs_u();
       int kn2 = srfs_[idxsf]->numCoefs_v();
       int kk1 = srfs_[idxsf]->order_u();
       int kk2 = srfs_[idxsf]->order_v();
       bool israt = srfs_[idxsf]->rational();
       if (copy_coefs_)
	 scoef = coef_array_[idxsf].begin();
       else
	 scoef = (israt) ? srfs_[idxsf]->rcoefs_begin() 
	   : srfs_[idxsf]->coefs_begin();

       integralInfo *cig = &surf_integral_[idxsf];

       double tval;    /* Value of the complete smoothness functional. */
       double tval1 = 0.0, tval2 = 0.0, tval3 = 0.0;  // Contributions to the 
                                                      // smoothness functional.
       double *sc;  // Pointer into coefficient array of the original surface. 

       int order = std::max(kk1,kk2);
       int ksz1 = std::min(kk1, kn1-kk1) + kk1;
       int ksz2 = std::min(kk2, kn2-kk2) + kk2;
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

       if (cig->integralset == false)
	 {
	     GaussQuadInner(srfs_[idxsf]->basis_u(), 
			    ider_, ta1, ta2, 
			    cig->integral1);

	     GaussQuadInner(srfs_[idxsf]->basis_v(), 
			    ider_, tb1, tb2, 
			    cig->integral2);
	 }

       double *sbder = boundary2 + 4*kk22;  // Storage of derivatives of B-splines
       // at the boundary of the surface.  

       // Compute boundary contibutions to the matrices
       // Compute all derivatives of B-splines up to order ider_-1
       // at the start of the parameter interval in 1. par. dir.

       srfs_[idxsf]->basis_u().computeBasisValues(ta1, sbder, ider_-1);

       for (k1=0; k1<kk1; k1++)
	 for (k2=0; k2<kk1; k2++)
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

       srfs_[idxsf]->basis_u().computeBasisValues(ta2, sbder, ider_-1);
       kleft = srfs_[idxsf]->basis_u().lastKnotInterval();

       kstart = std::min(kleft-kk1+1,kk1);
       for (ki=kstart, k1=0; k1<kk1; ki++, k1++)
	 for (kp=kstart, k2=0; k2<kk1; kp++, k2++)
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


       srfs_[idxsf]->basis_v().computeBasisValues(tb1, sbder, ider_-1);

       for (k3=0; k3<kk2; k3++)
	 for (k4=0; k4<kk2; k4++)
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

       srfs_[idxsf]->basis_v().computeBasisValues(tb2, sbder, ider_-1);
       kleft = srfs_[idxsf]->basis_v().lastKnotInterval();

       kstart = std::min(kleft-kk2+1,kk2);
       for (kj=kstart, k3=0; k3<kk2; kj++, k3++)
	 for (kq=kstart, k4=0; k4<kk2; kq++, k4++)
	   {
	     // Compute contribution from the 2. boundary in 2. par. dir.

	     boundary2[kj*ksz2+kq] += sbder[k3*ider_]*sbder[1+k4*ider_];
	     boundary2[kk22+kj*ksz2+kq] += sbder[1+k3*ider_]*sbder[k4*ider_];
	     boundary2[2*kk22+kj*ksz2+kq] += sbder[1+k3*ider_]*sbder[2+k4*ider_];
	     boundary2[3*kk22+kj*ksz2+kq] += sbder[2+k3*ider_]*sbder[1+k4*ider_];
	   }

       // Travers all B-splines and set up matrices of equation system.

       for (kl2=0, kq=0; kq<kn2; kq++)
	 for (kp=0; kp<kn1; kp++)
	   {
	     if (coef_known_[idxsf][kq*kn1+kp] == 1 || 
		 coef_known_[idxsf][kq*kn1+kp] == 2)
	       continue;

	     kl2 = (coef_known_[idxsf][kq*kn1+kp] > 2) ? 
	       pivot_[idxsf][coef_known_[idxsf][kq*kn1+kp]]
	       : pivot_[idxsf][kq*kn1+kp];

	     kjstart = std::max(0, kq-kk2+1);
	     kjend = std::min(kq+kk2,kn2);

	     for (kl1=0, kj=kjstart; kj<kjend; kj++)
	       {
		 if (kj<kk2 && kq<kk2)
		   {
		     k3 = kj; k4 = kq;
		   }
		 else if (kj >= kn2-kk2 && kq >= kn2-kk2)
		   {
		     k3 = kj - kn2 + kk2 + std::min(kk2, kn2-kk2);
		     k4 = kq - kn2 + kk2 + std::min(kk2, kn2-kk2);
		   }
		 else
		   {
		     k3 = 0; k4 = kk22-1;
		   }

		 kistart = std::max(0, kp-kk1+1);
		 kiend = std::min(kp+kk1,kn1);

		 for (ki=kistart; ki<kiend; ki++)
		   {
		     if (coef_known_[idxsf][kj*kn1+ki] == 2)
		       continue;

		     if (ki<kk1 && kp<kk1)
		       {
			 k1 = ki; k2 = kp;
		       }
		     else if (ki >= kn1-kk1 && kp >= kn1-kk1)
		       {
			 k1 = ki - kn1 + kk1 + std::min(kk1, kn1-kk1);
			 k2 = kp - kn1 + kk1 + std::min(kk1, kn1-kk1);
		       }
		     else
		       {
			 k1 = 0; k2 = kk11-1;
		       }

		     kl1 = (coef_known_[idxsf][kj*kn1+ki] > 2) ? 
		       pivot_[idxsf][coef_known_[idxsf][kj*kn1+ki]] 
		       : pivot_[idxsf][kj*kn1+ki];

		     if (kl2 > kl1 && coef_known_[idxsf][kj*kn1+ki] != 1) 
		       continue;

		     // Compute an element in the left side matrix.

		     // Compute contribution from 3. order term.

		     if (ider_ > 2)
		       tval3 = 5.0*
			 (cig->integral1[3][ki][kp]*cig->integral2[0][kj][kq]
			  + cig->integral1[0][ki][kp]*cig->integral2[3][kj][kq])
			 + 9.0*
			 (cig->integral1[2][ki][kp]*cig->integral2[1][kj][kq]
			  + cig->integral1[1][ki][kp]*cig->integral2[2][kj][kq])
			 + 3.0*
			 ((boundary1[3*kk11+k1*ksz1+k2]-cig->integral1[2][ki][kp])
			  *(boundary2[k3*ksz2+k4]-cig->integral2[1][kj][kq])
			  +(boundary1[2*kk11+k1*ksz1+k2]-cig->integral1[2][ki][kp])
			  *(boundary2[kk22+k3*ksz2+k4]-cig->integral2[1][kj][kq])
			  + (boundary1[kk11+k1*ksz1+k2]-cig->integral1[1][ki][kp])
			  *(boundary2[2*kk22+k3*ksz2+k4]-cig->integral2[2][kj][kq])
			  + (boundary1[k1*ksz1+k2]-cig->integral1[1][ki][kp])
			  *(boundary2[3*kk22+k3*ksz2+k4]-cig->integral2[2][kj][kq]));

		     // Compute contribution from 2. order term.

		     if (ider_ > 1)
		       tval2 = 4.0*cig->integral1[1][ki][kp]*
			 cig->integral2[1][kj][kq]
			 + 3.0*
			 (cig->integral1[2][ki][kp]*cig->integral2[0][kj][kq]
			  + cig->integral1[0][ki][kp]*cig->integral2[2][kj][kq])
			 + (boundary1[kk11+k1*ksz1+k2]-cig->integral1[1][ki][kp])
			 *(boundary2[k3*ksz2+k4]-cig->integral2[1][kj][kq])
			 + (boundary1[k1*ksz1+k2]-cig->integral1[1][ki][kp])
			 *(boundary2[kk22+k3*ksz2+k4]-cig->integral2[1][kj][kq]);

		     // Compute contribution from 1. order term.

		     if (ider_ > 0)
		       tval1 = cig->integral1[1][ki][kp]*
			 cig->integral2[0][kj][kq]
			 + cig->integral1[0][ki][kp]*cig->integral2[1][kj][kq];

		     // Complete term of the smoothness function.

		     tval = const1*tval1 + const2*tval2 + 2.0*const3*tval3;

		     if (coef_known_[idxsf][kj*kn1+ki] == 1)
		       {
			 // The contribution of this term is added to the right
			 // side of the equation system. First fetch the known
			 // coefficient.

			 sc = &*scoef + (kj*kn1 + ki)*idim1_;
			 for (kr=0; kr<idim_; kr++)
			   gright_[kr*kncond_+kl2] -= sc[kr]*tval;
		       }
		     else
		       {
			 // The contribution of this term is added to the left
			 //  side of the equation system.

			 if (kl2 > kl1) continue;

			 for (kk=0; kk<kdim_; kk++)
			   {
			     gmat_[(kk*kncond_+kl1)*kdim_*kncond_+kk*kncond_+kl2] += tval;
			     if (kl2 < kl1)
			       gmat_[(kk*kncond_+kl2)*kdim_*kncond_+kk*kncond_+kl1] += tval;
			   }
		       }
		   }
	       }
	   }
     }

   return;
}

/***************************************************************************/

void
SmoothSurfSet::getBasis(double *sb1, double *sb2, int kk1, int kk2,
			  int kleft1, int kleft2, int ider, double *sbasis)
//--------------------------------------------------------------------------
//     Purpose : Given all non-zero B-spline basis functions, compute the
//               corresponding basis functions.
//
//
//     Written by : Vibeke Skytt,  SINTEF,  09.99.
//--------------------------------------------------------------------------
{
   int kk12 = kk1*kk2;
   int ider1 = ider + 1;

   int ki, kj, kd1, kd2, kl1, kl2, k1, k2, ks1;
   double *s1 = sbasis;

   ks1 = ider;
   for (kd2=0; kd2<=ider; kd2++, ks1--)
     for (kd1=ks1; kd1<=ider-kd2; kd1++, s1+=kk12)
       for (kj=kleft2-kk2+1, kl2=kd2, k2=0; k2<kk12; kj++, kl2+=ider1, k2+=kk1)
	 for (ki=kleft1-kk1+1, kl1=kd1, k1=0; k1<kk1; ki++, kl1+=ider1, k1++)
	     s1[k2+k1] = sb1[kl1]*sb2[kl2];

   return;
}



/****************************************************************************/

void SmoothSurfSet::
setLeastSquares(std::vector<std::vector<double> >&  pnts, // Data points.   
		std::vector<std::vector<double> >&  param_pnts, // Parametrization.
		std::vector<std::vector<double> >&  pnt_weights,
		double wgt)  // Weight of current term.     
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		    from the approximation of data points.
//
//
//     Written by : Vibeke Skytt,  SINTEF SI,  09.93.
//--------------------------------------------------------------------------
{
  int kk;
  int k1, k2, k3, k4, k5, k6, k7, k8;
  int kl1, kl2;
  int kleft1=0, kleft2=0;  // Parameter used in s1220 to be positioned
                           // in the knot vector.                           
  double tz;     // Help variable.  
  double tval;   // Contribution to the matrices of the minimization problem.
  double *sc;    // Pointer into the coefficient array of the original surf.
  double const1 = (double)2.0*wgt;
  int idxsf;
  int numsfs = (int)srfs_.size();  // Number of surfaces in the surface set.
  vector<double>::iterator  scoef; // Pointer to surface coefficients. 

   // For each surface add the contribution of the least squares term to 
   // the equation system.
   
   for (idxsf=0; idxsf<numsfs; idxsf++)
     {
	 int nmbpoint = (int)pnts[idxsf].size()/idim_;   // Number of data points. 
       int kn1 = srfs_[idxsf]->numCoefs_u();
       //int kn2 = srfs_[idxsf]->numCoefs_v();
       int kk1 = srfs_[idxsf]->order_u();
       int kk2 = srfs_[idxsf]->order_v();
       bool israt = srfs_[idxsf]->rational();
       if (copy_coefs_)
	 scoef = coef_array_[idxsf].begin();
       else
	 scoef = (israt) ? srfs_[idxsf]->rcoefs_begin() 
	   : srfs_[idxsf]->coefs_begin();

       // Allocate scratch for B-spline basis functions. 
       vector<double> scratch(kk1+kk2+kk1*kk2, 0.0);
       double *sb1 = &scratch[0];  // Storage of B-spline basis functions. 
       double *sb2 = sb1+kk1;          // Storage of B-spline basis functions. 
       double *sbasis = sb2+kk2;       // Surface basis functions.

       // Traverse all points in the pointset.  
       double *pnt = &pnts[idxsf][0];
       double *par = &param_pnts[idxsf][0];
       for (int kr=0; kr<nmbpoint; kr++, pnt+=idim_, par+=2)
	 {
	   // Fetch B-spline basis functions different from zero. 
	   srfs_[idxsf]->basis_u().computeBasisValues(par[0], sb1, 0);
	   srfs_[idxsf]->basis_v().computeBasisValues(par[1], sb2, 0);
	   kleft1 = srfs_[idxsf]->basis_u().lastKnotInterval();
	   kleft2 = srfs_[idxsf]->basis_v().lastKnotInterval();

	   // Compute the surface basis functions.
	   getBasis(sb1, sb2, kk1, kk2, kleft1, kleft2, 0, sbasis);

	   for (k1=kleft1-kk1+1, k3=0; k1<=kleft1; k1++, k3++)
	     for (k2=kleft2-kk2+1, k4=0; k2<=kleft2; k2++, k4+=kk1)
	       {
		 // Test if the the coefficient is free to change.
		 if (coef_known_[idxsf][k2*kn1+k1] == 1 || 
		     coef_known_[idxsf][k2*kn1+k1] == 2)
		   continue;

		 kl1 = (coef_known_[idxsf][k2*kn1+k1] > 2) ?
		   pivot_[idxsf][coef_known_[idxsf][k2*kn1+k1]] 
		   : pivot_[idxsf][k2*kn1+k1];

		 tz = pnt_weights[idxsf][kr]*sbasis[k4+k3];

		 // Add contribution to right hand side. 
		 for (kk=0; kk<idim_; kk++)
		   {
		     tval = const1*pnt[kk]*tz;
		     gright_[kk*kncond_+kl1] += tval;
		   }

		 for (k5=kleft1-kk1+1, k7=0; k5<=kleft1; k5++, k7++)
		   for (k6=kleft2-kk2+1, k8=0; k6<=kleft2; k6++, k8+=kk1)
		     {
		       if (coef_known_[idxsf][k6*kn1+k5] == 2)
			 continue;

		       // Compute contribution to left hand side. 
		       tval = const1*tz*sbasis[k8+k7];

		       // Test if the current coefficient is at the boundary.
		       if (coef_known_[idxsf][k6*kn1+k5] == 1)
			 {
			   // Adjust for boundary conditions. 
			   sc = &*scoef + (k6*kn1 + k5)*idim1_;
			   for (kk=0; kk<idim_; kk++)
			     gright_[kk*kncond_+kl1] -= sc[kk]*tval;
			 }
		       else
			 {
			   // The term gives a contribution on the left hand side. 
			   kl2 = (coef_known_[idxsf][k6*kn1+k5] > 2) ?
			     pivot_[idxsf][coef_known_[idxsf][k6*kn1+k5]] : 
			     pivot_[idxsf][k6*kn1+k5];

			   if (kl2 > kl1) continue;

			   for (kk=0; kk<kdim_; kk++)
			     {
			       gmat_[(kk*kncond_+kl1)*kdim_*kncond_+kk*kncond_+kl2] += tval;
			       if (kl2 < kl1)
				 gmat_[(kk*kncond_+kl2)*kdim_*kncond_+kk*kncond_+kl1] += tval;
			     }
			 }
		     }
	       }
	 }
     }

    return;
 }


/****************************************************************************/

int
SmoothSurfSet::setNormalCond(std::vector<std::vector<double> >& pnts,
			       std::vector<std::vector<double> >& param_pnts,
			       std::vector<std::vector<double> >&   pnt_weights,
			       double weight)
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of normal directions.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF,  9910.
//--------------------------------------------------------------------------
{
  if (kdim_ != 3)
    return 1;   // Not prepared for approximation of normals.

  int kk, kb;
  int k1, k2, k3, k4, k5, k6, k7, k8;
  int kl1, kl2;
  int kleft1=0, kleft2=0;  // Parameter used in s1220 to be positioned
                           // in the knot vector.                           
  int kder = 1;            // Number of derivatives of B-splines to compute.
  int kder1 = kder+1; 
  double tz1, tz2;         // Help variables.  
  double tval;   // Contribution to the matrices of the minimization problem.
  double *sc;    // Pointer into the coefficient array of the original surf.
  double const1 = (double)2.0*weight;
  double tdum;

  int idxsf;
  int numsfs = (int)srfs_.size();  // Number of surfaces in the surface set.
  vector<double>::iterator  scoef; // Pointer to surface coefficients. 

   // For each surface add the contribution of the least squares term to 
   // the equation system.
   
   for (idxsf=0; idxsf<numsfs; idxsf++)
     {
       int nmbpoint = (int)pnts[idxsf].size()/idim_;   // Number of data points. 
       int kn1 = srfs_[idxsf]->numCoefs_u();
       //int kn2 = srfs_[idxsf]->numCoefs_v();
       int kk1 = srfs_[idxsf]->order_u();
       int kk2 = srfs_[idxsf]->order_v();
       int kk12 = kk1*kk2;
       bool israt = srfs_[idxsf]->rational();
       if (copy_coefs_)
	 scoef = coef_array_[idxsf].begin();
       else
	 scoef = (israt) ? srfs_[idxsf]->rcoefs_begin() 
	   : srfs_[idxsf]->coefs_begin();

       // Allocate scratch for B-spline basis functions. 
  
       vector<double> scratch((kk1+kk2+kk1*kk2*kder1)*kder1, 0.0);
       double *sb1 = &scratch[0];  // Storage of B-spline basis functions. 
       double *sb2 = sb1+kk1*kder1;    // Storage of B-spline basis functions. 
       double *sbasis = sb2+kk2*kder1;       // Surface basis functions.

       // Traverse all normal directions.

       double *pnt = &pnts[idxsf][0];
       double *par = &param_pnts[idxsf][0];
       for (int kr=0; kr<nmbpoint; kr++, pnt+=idim_, par+=2)
	 {
	   // Fetch B-spline basis functions different from zero. 

	   srfs_[idxsf]->basis_u().computeBasisValues(par[0], sb1, kder);
	   srfs_[idxsf]->basis_v().computeBasisValues(par[1], sb2, kder);
	   kleft1 = srfs_[idxsf]->basis_u().lastKnotInterval();
	   kleft2 = srfs_[idxsf]->basis_v().lastKnotInterval();

	   // Compute the surface basis functions.
	   getBasis(sb1, sb2, kk1, kk2, kleft1, kleft2, kder, sbasis);

	   for (k1=kleft1-kk1+1, k3=0; k1<=kleft1; k1++, k3++)
	     for (k2=kleft2-kk2+1, k4=0; k2<=kleft2; k2++, k4+=kk1)
	       {
		 // Test if the the coefficient is free to be changed.

		 if (coef_known_[idxsf][k2*kn1+k1] == 1 || 
		     coef_known_[idxsf][k2*kn1+k1] == 2)
		   continue;

		 kl1 = (coef_known_[idxsf][k2*kn1+k1] > 2) ?
		   pivot_[idxsf][coef_known_[idxsf][k2*kn1+k1]] 
		   : pivot_[idxsf][k2*kn1+k1];

		 tz1 = pnt_weights[idxsf][kr]*sbasis[k4+k3];
		 tz2 = pnt_weights[idxsf][kr]*sbasis[kk12+k4+k3];

		 for (k5=kleft1-kk1+1, k7=0; k5<=kleft1; k5++, k7++)
		   for (k6=kleft2-kk2+1, k8=0; k6<=kleft2; k6++, k8+=kk1)
		     {
		       if (coef_known_[idxsf][k6*kn1+k5] == 2)
			 continue;

		       // Compute contribution to left hand side. 

		       tval = const1*(tz1*sbasis[k7+k8] +
				      tz2*sbasis[kk12+k7+k8]);

		       // Test if the current coefficient is at the boundary.

		       if (coef_known_[idxsf][k6*kn1+k5] == 1)
			 {
			   // Adjust for boundary conditions. 

			   sc = &*scoef + (k6*kn1 + k5)*idim1_;
			   for (tdum=(double)0, kk=0; kk<idim_; kk++)
			     tdum += sc[kk]*pnt[kk];

			   for (kk=0; kk<idim_; kk++)
			     gright_[kk*kncond_+kl1] -= tval*tdum*pnt[kk];
			 }
		       else
			 {
			   // The term gives a contribution on the left hand side. 

			   kl2 = (coef_known_[idxsf][k6*kn1+k5] > 2) ?
			     pivot_[idxsf][coef_known_[idxsf][k6*kn1+k5]] : 
			     pivot_[idxsf][k6*kn1+k5];

			   if (kl2 > kl1) continue;

			   for (kk=0; kk<kdim_; kk++)
			     {
			       for (kb=0; kb<kdim_; kb++)
				 {
				   gmat_[(kk*kncond_+kl1)*kdim_*kncond_+kk*kncond_+kl2] 
				     += tval*pnt[kk]*pnt[kb];
				   if (kl2 < kl1)
				     gmat_[(kk*kncond_+kl2)*kdim_*kncond_+kk*kncond_+kl1] += 
				       tval*pnt[kk]*pnt[kb];
				 }
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
SmoothSurfSet::approxOrig(double weight)  // Weight of current term.     
//--------------------------------------------------------------------------
//     Purpose : Compute the contribution to the equation system
//		 from the approximation of an original surface
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF,  9808.
//--------------------------------------------------------------------------
{
  int ki, kj, kr;

  int kl1;
  double wgt2 = (double)2*weight;
  int idxsf;
  int numsfs = (int)srfs_.size();  // Number of surfaces in the surface set.
  vector<double>::iterator  scoef; // Pointer to surface coefficients. 

  // For each surface add the contribution of the approximation of
  // the original surface to the equation system.
   
  for (idxsf=0; idxsf<numsfs; idxsf++)
    {
      int kn1 = srfs_[idxsf]->numCoefs_u();
      int kn2 = srfs_[idxsf]->numCoefs_v();
      bool israt = srfs_[idxsf]->rational();
      if (copy_coefs_)
	scoef = coef_array_[idxsf].begin();
      else
	scoef = (israt) ? srfs_[idxsf]->rcoefs_begin() 
	  : srfs_[idxsf]->coefs_begin();

   
      // Traverse all coefficients and add the contribution of the
      //  approximation of the original coefficients to the matrix.   
   
      for (kj=0; kj<kn2; kj++)
	for (ki=0; ki<kn1; ki++)
	  {
	    if (coef_known_[idxsf][kj*kn1+ki] == 1 || 
		coef_known_[idxsf][kj*kn1+ki] == 2)
	      continue;

	    kl1 = (coef_known_[idxsf][kj*kn1+ki] > 2) ?
	      pivot_[idxsf][coef_known_[idxsf][kj*kn1+ki]] : 
	      pivot_[idxsf][kj*kn1+ki];
	 
	    // Add the contribution to the right hand side. 
	 
	    for (kr = 0; kr<idim_; kr++)
	      gright_[kr*kncond_+kl1] += wgt2*scoef[(kj*kn1+ki)*idim1_+kr];
	 
	    // Add the contribution to the left hand side. 
	    
	    for (kr=0; kr<kdim_; kr++)
	      gmat_[(kr*kncond_+kl1)*kdim_*kncond_+kr*kncond_+kl1] += wgt2;

	  }
    }
   
}

/****************************************************************************/

void
SmoothSurfSet::setSideConstraints(std::vector<sideConstraintSet>& constraints)
//--------------------------------------------------------------------------
//     Purpose : Set linear side constraints to the minimization problem,
//               using Lagrange multiplier.
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
    ALWAYS_ERROR_IF(kdim_ == 3,
		// @@sbr Assuming normal conditions are not included in the system!
		// If this is not the case, all calculations concerning gmat_ must
		// be rewritten. Hence kdim_ is assumed to be 1, I suppose...?
		"Not prepared for normal conditions!");
    // We check whether constraints are achievable. If coef in constraint is
    // already known, add to const term on right side. If all coefs are removed,
    // delete constraint. Issue warning if constraint is impossible.
    try {
      removeKnownCoefs(constraints);
    } catch (...) {
      THROW("Failed updating side constraints.");
    }

    // If any side constraints were removed, we must update size of matrix.
    if (int(constraints.size()) != knconstraint_) {
	int new_knconstraint = (int)constraints.size();
	int new_kncond = kncond_ - (knconstraint_ - new_knconstraint);
	// For ease of algorithm, we copy matrices to new matrices.
	vector<double> new_gmat(kdim_*kdim_*new_kncond*new_kncond);
	vector<double> new_gright(idim_*new_kncond);
	// @@sbr If we use normal conditions, these copies will be inadequate!
	//     for (i = 0; i < kdim_; ++i) // We treat one dimension at the time.
// 	int i = 0;
	for (int j = 0; j < new_kncond; ++j)
// 	    copy(gmat_.begin() + i*kdim_*kdim_*kncond_*kncond_ + j*kdim_*kncond_,
// 		 gmat_.begin() + i*kdim_*kdim_*kncond_*kncond_ + j*kdim_*kncond_ + new_kncond_,
// 		 new_gmat.begin() + j*new_kncond_);
	    copy(gmat_.begin() + j*kdim_*kncond_,
		 gmat_.begin() + j*kdim_*kncond_ + new_kncond,
		 new_gmat.begin() + j*new_kncond);
	for (int i = 0; i < idim_; ++i)
	    copy(gright_.begin() + i*kncond_,
		 gright_.begin() + i*kncond_ + new_kncond,
		 new_gright.begin() + i*new_kncond);
	knconstraint_ = new_knconstraint;
	kncond_ = new_kncond;
	gmat_ = new_gmat;
	gright_ = new_gright;
    }

    // @@ We're still assuming kdim_ == 1...
    // We update values in matries as described in the above equations.
    int nmb_free_coefs = kncond_ - knconstraint_; // Dimension of vector x (# coefs).
    // We start by updating gmat_ by adding new elements given by side constraints.
    for (size_t i = 0; i < constraints.size(); ++i)
	for (size_t j = 0; j < constraints[i].factor_.size(); ++j) { // We start with gmat_.
	    int surf_ind = constraints[i].factor_[j].first.first;
	    // We have made  sure that all elements in constraints[i] are free.
	    gmat_[kncond_*(nmb_free_coefs+i)+
		 pivot_[surf_ind][constraints[i].factor_[j].first.second]] =
		constraints[i].factor_[j].second;
	    gmat_[kncond_*pivot_[surf_ind][constraints[i].factor_[j].first.second]+
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

void
SmoothSurfSet::setApproxSideConstraints(std::vector<sideConstraintSet>& constraints,
					double weight)
//--------------------------------------------------------------------------
//     Purpose : Use least squares to minimize the error given by the constraints.
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF Oslo,  06.02.
//--------------------------------------------------------------------------
{
  if (constraints.size() == 0)
    return;

  // First we update the constraints with coef_known.
  try {
    removeKnownCoefs(constraints);
  } catch (...) {
    THROW("Failed updating side constraints.");
  }

  // As we would like to minimize ||Cx-d||, where x is the unknown coefs and C & d
  // are given by the constraints, we seek a solution to Ax = (C^T)Cx = (C^T)d = e.
  // We add values to both sides of the system. We multiply both terms by weight,
  // allowing easy adjustment of factor.
  int dim = constraints[0].dim_;
  int nmb_constraints = (int)constraints.size();
  int nn = kncond_ - knconstraint_;
  //Matrix C(nmb_constraints, nn);
  vector<vector<double> > Cmat(nmb_constraints);
  for (int ii = 0; ii < nmb_constraints; Cmat[ii++]= vector<double>(nn, 0)) {}
  
//   // We set all elements to 0.0;
//   for (int ki = 0; ki < nmb_constraints; ++ki) {
//     for (int kj = 0; kj < nn; ++kj)
//       C.element(ki, kj) = 0.0;
//   }

  // We add elements given by constraints to left side of equation.
  for (size_t ki = 0; ki < constraints.size(); ++ki) {
    // For each factor in constraint, we compute the element in the matrix C.
    for (size_t kj = 0; kj < constraints[ki].factor_.size(); ++kj) {
      int surf = constraints[ki].factor_[kj].first.first;
      int coef_idx = constraints[ki].factor_[kj].first.second;
      //C.element(ki, pivot_[surf][coef_idx]) = constraints[ki].factor_[kj].second;
      Cmat[ki][pivot_[surf][coef_idx]] = constraints[ki].factor_[kj].second;
    }
  }
  //Matrix A(nn, nn);
  vector<vector<double> > Amat(nn);
  for (int ii = 0; ii < nn; Amat[ii++].resize(nn)) {}
  //A = C.t()*C;
  for (int row = 0; row < nn; ++row) {
      for (int col = 0; col < nn; ++col) {
	  for (int k = 0; k < nmb_constraints; ++k) {
	      Amat[row][col] += Cmat[k][row] * Cmat[k][col];
	  }
      }
  }
  
  for (int ki = 0; ki < nn; ++ki)
    for (int kj = 0; kj < nn; ++kj)
	//gmat_[ki*kncond_+kj] += A.element(ki, kj)*weight;
	gmat_[ki*kncond_+kj] += Amat[ki][kj]*weight;

  // We add elements given by constraints to right side of equation.
  for (int ki = 0; ki < dim; ++ki) {
      //ColumnVector d(nmb_constraints);
      vector<double> d(nmb_constraints);
      for (size_t kj = 0; kj < constraints.size(); ++kj) {
	  //d.element(kj) = constraints[kj].constant_term_[ki];
	  d[kj] = constraints[kj].constant_term_[ki];
      }
      //ColumnVector e(nn);
      vector<double> e(nn, 0);
      //e = C.t()*d;
      for (int e_ix = 0; e_ix < nn; ++e_ix) {
	  for (int c_k = 0; c_k < nmb_constraints; ++c_k) {
	      e[e_ix] += Cmat[c_k][e_ix] * d[c_k];
	  }
      }

      for (int kj = 0; kj < nn; ++kj) {
	  //gright_[kncond_*ki+kj] += e.element(kj)*weight;
	  gright_[kncond_*ki+kj] += e[kj]*weight;
      }
  }
}


/****************************************************************************/

int
SmoothSurfSet::equationSolve(std::vector<shared_ptr<SplineSurface> >& surfaces) 
                                              // Resulting surfaces,
                                              // NB! the array must be empty.
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
   double weight;
   int idxsf;
   int numsfs = (int)srfs_.size();  // Number of surfaces in the surface set.
   vector<double>::iterator  scoef; // Pointer to surface coefficients. 
   int kn1, kn2, kn12;


//    FILE *fp = 0;
//    fp = fopen("fA.m", "w");
//    fprintf(fp,"A=[ ");
//    for (kj=0; kj<kncond_; kj++)
//      {
//        for (ki=0; ki<kncond_; ki++)
// 	 fprintf(fp, "%18.7f", gmat_[kj*kncond_+ki]);
//        if (kj<kncond_-1) fprintf(fp,"\n");
//      }
//    fprintf(fp," ]; \n");
//    fclose(fp);

   // Solve the equation system by Conjugate Gradient Method.

   vector<double> eb(idim_*kncond_, 0.0);
   for (ki=0; ki<idim_*kncond_; ki++)
     eb[ki] = gright_[ki];

       // Copy coefficients to array of unknowns

   for (idxsf=0; idxsf<numsfs; idxsf++)
     {
       bool israt = srfs_[idxsf]->rational();
       if (copy_coefs_)
	 scoef = coef_array_[idxsf].begin();
       else
	 scoef = (israt) ? srfs_[idxsf]->rcoefs_begin() 
	   : srfs_[idxsf]->coefs_begin();

       kn1 = srfs_[idxsf]->numCoefs_u();
       kn2 = srfs_[idxsf]->numCoefs_v();
       kn12 = kn1*kn2;
       for (ki=0; ki<kn12; ki+=kn1)
	 for (kj=0; kj<kn1; kj++)
	   {
	     if (coef_known_[idxsf][ki+kj] == 1 || 
		 coef_known_[idxsf][ki+kj] == 2)
	       continue;
	   
	     kl1 = (coef_known_[idxsf][ki+kj] > 2) ?
	       pivot_[idxsf][coef_known_[idxsf][ki+kj]] : 
	       pivot_[idxsf][ki+kj];

	     for (kk=0; kk<idim_; kk++)
	       gright_[kk*kncond_+kl1] = scoef[(ki+kj)*idim1_+kk]; 
	   }
     }
       
   // Set up CG-object

   SolveCG solveCg;

   // Create sparse matrix.

   solveCg.attachMatrix(&gmat_[0], kdim_*kncond_);

   // Attach parameters.

   solveCg.setTolerance(0.00000001);
   // @@sbr It seems precondRILU() assumes diagonal elements of matrix are non-zero.
   //       If constraints exist, these may be non-zero on the diagonal.
   int precond = (knconstraint_ > 0) ? 0 : 1;
   // When system has been preconditioned, we expect to solve it in kdim_*kncond_
   // iterations. However, when this is not the case, we may need more iterations.
   int nmb_iter = precond ? kdim_*kncond_ : 100*kdim_*kncond_;
   solveCg.setMaxIterations(std::min(nmb_iter, 10000));
       // Preconditioning

   //        printf("Precondintioning (0/1) ? ");
   //        scanf("%d",&precond);
   if (precond)
       {
	   double omega = 0.1;
	   // 	   printf("Omega = ");
	   // 	   scanf("%lf",&omega);
	   solveCg.precondRILU(omega);
       }

       // Solve equation systems.
       
   if (kdim_ == 3)
     {
       kstat = solveCg.solve(&gright_[0], &eb[0], kdim_*kncond_);
       //	   printf("solveCg.solve status %d \n", kstat);
       if (kstat < 0) 
	 return kstat;
       if (kstat == 1)
	 // @@ Handle this from the outside, try with twice the number of iterations?
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
	     THROW("Failed solving system (within tolerance)!");
	 }
     }

   // Copy result to output array. 
   // @@@ VSK. Is it necessary to take the Lagrange multipliers
   // into acount?

   // For each surface fetch the new coefficient vector.
   
   int kk1, kk2;  // Order of current surface in both parameter directions
   vector<double>::const_iterator st1;  // Knot vector in 1. parameter direction
   vector<double>::const_iterator st2;  // Knot vector in 2. parameter direction
   for (idxsf=0; idxsf<numsfs; idxsf++)
     {
       bool israt = srfs_[idxsf]->rational();
       if (copy_coefs_)
	 scoef = coef_array_[idxsf].begin();
       else
	 scoef = (israt) ? srfs_[idxsf]->rcoefs_begin() 
	   : srfs_[idxsf]->coefs_begin();
       kn1 = srfs_[idxsf]->numCoefs_u();
       kn2 = srfs_[idxsf]->numCoefs_v();
       kn12 = kn1*kn2;
       kk1 = srfs_[idxsf]->order_u();
       kk2 = srfs_[idxsf]->order_v();
       st1 = srfs_[idxsf]->basis_u().begin();
       st2 = srfs_[idxsf]->basis_v().begin();

       for (ki=0; ki<kn12; ki+=kn1)
	 for (kj=0; kj<kn1; kj++)
	   {
	     if (coef_known_[idxsf][ki+kj] == 1 || 
		 coef_known_[idxsf][ki+kj] == 2)
	       continue;
	   
	     kl1 = (coef_known_[idxsf][ki+kj] > 2) ?
	       pivot_[idxsf][coef_known_[idxsf][ki+kj]] : 
	       pivot_[idxsf][ki+kj];

	     for (kk=0; kk<idim_; kk++)
	       scoef[(ki+kj)*idim1_+kk] = gright_[kk*kncond_+kl1]; 
	   }

       // Multiply the coefficients with weights if necessary.
       if (idim1_ > idim_)
	 {
	   for (ki=0; ki<kn12; ki+=kn1)
	     for (kj=0; kj<kn1; kj++)
	       {
		 weight = scoef[(ki+kj)*idim1_+idim_];

		 for (kk=0; kk<idim_; kk++)
		   scoef[(ki+kj)*idim1_+kk] *= weight;
	       }
	 }

       // Create SplineSurface. 

       shared_ptr<SplineSurface> surf;
       surf = shared_ptr<SplineSurface>(new SplineSurface(kn1,kn2,kk1,kk2,
							      st1,st2,scoef,
							      idim_, israt));
       surfaces.push_back(surf);
    }
   return 0;
}


/****************************************************************************/

void
SmoothSurfSet::removeKnownCoefs(std::vector<sideConstraintSet>& constraints) const
//--------------------------------------------------------------------------
//     Purpose : Insert values for known coefs, remove constraints when all are known
//
//     Calls   : 
//
//     Written by : Sverre Briseid,  SINTEF Oslo,  06.02.
{
    double equality_eps = 1e-06;
    for (size_t i = 0; i < constraints.size(); ++i) {
	int dim = constraints[i].dim_;
	ALWAYS_ERROR_IF(dim != idim_,
			"Constraint not member of the same geometric space as surface!");

	for (size_t j = 0; j < constraints[i].factor_.size(); ++j) {
	    int surf_ind = constraints[i].factor_[j].first.first;
	    int coef_ind = constraints[i].factor_[j].first.second;
	    vector<double>::iterator scoef = srfs_[surf_ind]->coefs_begin();
	    if (coef_known_[surf_ind][coef_ind] != 0) {
		// As coef was already known, we add to const term on right side.
		for (int k = 0; k < dim; ++k)
		    constraints[i].constant_term_[k] -=
			scoef[dim*coef_ind+k]*constraints[i].factor_[j].second;
		constraints[i].factor_.erase(constraints[i].factor_.begin() + j);
		--j;
	    }
	}
	// If all coefs in constraint were known, we remove constraint.
	if (constraints[i].factor_.size() == 0) {
	    double sum = 0.0;
	    for (int j = 0; j < dim; ++j)
		sum += constraints[i].constant_term_[j];
	    if (fabs(sum) > equality_eps) // Check whether constraint was fulfilled.
		MESSAGE("All coefs already known, requiring: 0.0 = " << sum);
	    constraints.erase(constraints.begin() + i); // We remove constraint.
	    --i;
	}
    }
}
