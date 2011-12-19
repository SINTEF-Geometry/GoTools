// ----------------------------------------------------------------
//       Implementation file for class SmoothCurve.
// ----------------------------------------------------------------
//
// This file contains the following member functions:
//             a. SmoothCurve()
//             b. ~SmoothCurve()
//             c. attach()
//             e. setOptim()
//             f. setLeastSquares()
//             g. equationSolve
// ----------------------------------------------------------------

#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/creators/Integrate.h"
//#include "newmat.h"
#include "GoTools/utils/LUDecomp.h"

#include "GoTools/utils/Values.h"
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
#include "GoTools/geometry/Utils.h"     // make std::min and std::max work (redefined in boost/smart_ptr.hpp)
#endif

#include <math.h>
#include <algorithm>


using namespace Go;
using std::vector;
using std::max;
using std::min;

SmoothCurve::SmoothCurve(int dim)
  //--------------------------------------------------------------------------
  //     Constructor for class SmoothCurve.
  //
  //     Purpose : Initialize class variables
  //
  //     Calls   :
  //
  //     Written by : Vibeke Skytt,  SINTEF,  04.98.
  //--------------------------------------------------------------------------
{
  // Set values of variables global to the class.  

  ider_ = 3;
  idim_ = kdim_ = dim;
  cont_seam_ = 0;

  kn_ = kk_ = 0;
}

//***************************************************************************


SmoothCurve::~SmoothCurve()
  //--------------------------------------------------------------------------
  //     Destructor for class SmoothCurve.
  //
  //     Written by : Vibeke Skytt,  SINTEF SI,  09.93.
  //--------------------------------------------------------------------------
{
  releaseScratch();
  // if (qcurve_) freeCurve(qcurve_);
}

void SmoothCurve::releaseScratch()
  //--------------------------------------------------------------------------
  //     Free all memory allocated for the class members of SmoothCurve.
  //
  //     Written by : Sverre Briseid,  SINTEF,  06.02
  //--------------------------------------------------------------------------
{
  int ki, kj;

  for (ki=0; ki<=ider_; ki++)
    {
      for (kj=0; kj<kn_; kj++)
	if (integral_[ki][kj]) delete [] integral_[ki][kj];

      if (integral_[ki]) delete [] integral_[ki];
    }
  if (integral_) delete [] integral_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  SmoothCurve::attach
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

int
SmoothCurve::attach(const shared_ptr<SplineCurve>& incurve,       
		      // Input curve representing the spline space
		      // 		      int fixatend[],
		      int coef_known[],
		      int numSideConstraints)
  //--------------------------------------------------------------------------
  //
  //     Purpose : Initialize class variables related to the spline space.
  //
  //     Calls   :
  //
  //     Written by : Vibeke Skytt,  SINTEF,  04.98.
  //     Revised by :
  //--------------------------------------------------------------------------
{
  int ki, kj;
  int kstat;

  // Fetch data from B-spline curve.
  rational_ = incurve->rational();
  st_ = incurve->basis().begin();
  if (rational_)
    {
      scoef_ = incurve->rcoefs_begin();
      kdim_ = idim_ + 1;
      int in_coefs = incurve->numCoefs();
      int in_ord = incurve->order();
      vector<double> bspl_coefs(in_coefs*2, 0.0);
      for (int pos_in = idim_, pos_out = 1;
	   pos_in < kdim_ * in_coefs;
	   pos_in += kdim_, pos_out +=2)
	bspl_coefs[pos_out] = scoef_[pos_in];
      bspline_curve_ = shared_ptr<SplineCurve>(new SplineCurve(in_coefs, in_ord, st_, bspl_coefs.begin(), 1, true));
    }
  else
    {
      scoef_ = incurve->coefs_begin();
      kdim_ = idim_;
    }

  //    for (ki=0; ki<2; ki++)
  //    cont_bound[ki] = fixatend[ki];


  if (incurve.get() != 0 && incurve->order() == kk_ && incurve->numCoefs() == kn_ &&
      numSideConstraints == kconstraint_)
    {
      // It is assumed that the new spline curve (incurve) lies in the same
      // spline space that the old curve.
      // If these conditions are broken, memory faults and/or wrong
      // results will occur. Including higher order derivatives in the
      // smoothing functional than in previous runs, will have no effect.

      // Zero out the arrays of the equation system.
      for (ki = 0; ki < kcond_; ++ki)
	for (kj = 0; kj < kcond_; ++kj)
	    //gmat_.element(ki, kj) = 0.0;
	    gmat_[ki][kj] = 0.0;
      for (ki = 0; ki < idim_*kcond_; ++ki)
	  //gright_.element(ki) = 0.0;
	  gright_[ki] = 0.0;

      //        integralset = 1;
      qcurve_ = incurve;
    }
  else
    {
      kk_ = incurve->order();
      kn_ = incurve->numCoefs();

      if (qcurve_.get() != 0)
	{
	  // Memory already allocated for the previous iteration of
	  // curve editing and smoothing. Free this memory.
	  releaseScratch(); 
	}


      //        integralset = 0;
      qcurve_ = incurve;

      // Allocate scratch for pivot_ array.
      pivot_.resize(kn_);

      // Update coefknown_ according to information about continuity over
      // a seam. Also update coefficients at the seam if the required 
      // information exists.
   
      coefknown_ = coef_known;

      //        if (seam[0] > 0 || seam[1] > 0)
      // 	 preparePeriodicity(seam);

      // Store information about free, fixed and irrelevant coefficients
      // and create pivot_ array. Count number of free coefficients.
      kcond_ = 0;
      for (ki = 0; ki < kn_; ++ki)
	if (coefknown_[ki] == 0) {
	  pivot_[ki] = kcond_;
	  ++kcond_;
	} else
	  pivot_[ki] = -MAXINT;

      // Add side constraints
      kconstraint_ = numSideConstraints;
      kcond_ += kconstraint_;

      // Allocate scratch for arrays of integrals of inner product of 
      // B-splines.
      kstat = prepareIntegral();
      if (kstat < 0)
	return kstat;

      // Allocate scratch for arrays in the equation system. 
      //gmat_.ReSize(kcond_,kcond_);
      gmat_.resize(kcond_);
      for (int ii = 0; ii < kcond_; gmat_[ii++].resize(kcond_)) {}
      for (ki = 0; ki < kcond_; ++ki)
	for (kj = 0; kj < kcond_; ++kj)
	    //gmat_.element(ki, kj) = 0.0;
	    gmat_[ki][kj] = 0.0;
      //gright_.ReSize(idim_*kcond_);
      gright_.resize(idim_ * kcond_);
      for (ki = 0; ki < idim_*kcond_; ++ki)
	  //gright_.element(ki) = 0.0;
	  gright_[ki] = 0.0;
    }

  return 0;
}

//**************************************************************************

int SmoothCurve::prepareIntegral()

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
  int ki, kj, kr;

  // Allocate scratch for arrays of integrals of inner product of B-splines.
  integral_ = new double**[ider_+1];
  if (integral_ == NULL)
    return -101;

  for (ki=0; ki<=ider_; ki++)
    {
      integral_[ki] = NULL;
      integral_[ki] = new double*[kn_];
      if (integral_[ki] == NULL)
	return -101;
      for (kj=0; kj<kn_; kj++)
	{
	  integral_[ki][kj] = NULL;
	  integral_[ki][kj] = new double[kn_];
	  if (integral_[ki][kj] == NULL)
	    return -101;
	  for (kr=0; kr<kn_; kr++)
	    integral_[ki][kj][kr] = (double)0;
	}
    }

  return 0;
}


void
SmoothCurve::setOptim(const double wgt1,  // Weight of 1. order term 
			const double wgt2,  // Weight of 2. order term 
			const double wgt3)  // Weight of 3. order term 
  //--------------------------------------------------------------------------
  //     Purpose : Compute the contribution to the equation system
  //		    from the smoothness functional.
  //
  //     Calls   : Integrate::GaussQuad
  //               
  //
  //     Written by : Vibeke Skytt,  SINTEF,  28-04.98.
  //--------------------------------------------------------------------------
{
    //int grstat = 0;   // Initialize status variable.
  int ki, kp, kr;
  int kl1, kl2;
  double tval;      // Contribution to equation system from smoothness term.

  // Set parameter area.

  double ta = st_[kk_-1];   // Start of par. interval.
  double tb = st_[kn_];     // End of par. interval.

  double *sc;  // Pointer into coefficient array of the original curve.


  // Compute all integrals of inner product of B-splines.

  //Integrate GaussQuad;

  //grstat = 
  //GaussQuad.
  if (rational_)
    GaussQuadInnerRational(qcurve_->basis(), ider_, ta, tb, bspline_curve_, integral_);
  else
    GaussQuadInner(qcurve_->basis(), ider_, ta, tb, integral_);
//   if (grstat != 0)
//     return grstat;


  // Traverse all B-splines and set up matrices of equation system.

  for (ki=0; ki<kn_; ki++)
    for (kp=0; kp<kn_; kp++)
      {
	// @@sbr What about coef_known[ki] == 2?
	if ((kp > ki && coefknown_[ki] == 0) || (coefknown_[kp] == 1))
	  continue;
// 	kl1 = ki - cont_bound[0];
// 	kl2 = kp - cont_bound[0];

// 	if (kl2 > kl1 &&
// 	    !(ki < cont_bound[0] || kn_-ki-1 < cont_bound[1]))
// 	  continue;

// 	if (kp < cont_bound[0] || kn_-kp-1 < cont_bound[1])
// 	  continue;

	kl1 = pivot_[ki];
	kl2 = pivot_[kp];


	// Compute an element in the left side matrix.
	// Complete term of the smoothness function.

	tval = wgt1*integral_[1][ki][kp] + 
	  wgt2*integral_[2][ki][kp] + wgt3*integral_[3][ki][kp];

// 	if (ki < cont_bound[0] || kn_-ki-1 < cont_bound[1])
 	if (coefknown_[ki] == 1)
	  {
	    // The contribution of this term is added to the right
	    // side of the equation system. First fetch the known
	    // coefficient.

	    sc = &*scoef_ + ki*kdim_;
	    for (kr=0; kr<idim_; kr++)
		//gright_.element(kr*kcond_+kl2) -= sc[kr]*tval;
		gright_[kr*kcond_+kl2] -= sc[kr]*tval;
	  }
	else
	  {
	    // The contribution of this term is added to the left
	    //  side of the equation system.

	    //gmat_.element(kl1, kl2) += tval;
	    gmat_[kl1][kl2] += tval;
	    if (kl2 < kl1)
		//gmat_.element(kl2, kl1) += tval;
		gmat_[kl2][kl1] += tval;
	  }
      }

  return;
}

//****************************************************************************

void
SmoothCurve::setLeastSquares(std::vector<double>& pnts, // Data points.
			       std::vector<double>& param_pnts, // Par.
			       std::vector<double>& pnt_weights,
			       double wgt)  // Weight of current term.    
  //--------------------------------------------------------------------------
  //     Purpose : Compute the contribution to the equation system
  //		 from the approximation of data points.
  //
  //     Calls   : BsplineBasis::computeBasisValues
  //
  //     Written by : Vibeke Skytt,  SINTEF,  04.98.
  //--------------------------------------------------------------------------
{
    int num_points = (int)param_pnts.size();

  int k1, k2, k3, k4, kj;
  int kl1, kl2;
  int kleft=0;  // Parameter used in s1220 for positioning in the knot vector.
  double tz;     // Help variable.  
  double *sc;    // Pointer into the coefficient array of the original curve.

  // Allocate scratch for B-spline basis functions. 
  std::vector<double> sbasis(kk_);
  double *sb = &sbasis[0];

  // Traverse all points in the pointset. 

//   double *pnt = pnts;
  for (int kr =0; kr < num_points; ++kr)
    {

      // Fetch B-spline basis functions different from zero. 

      if (rational_)
	{
	  Point p(1);
	  kleft = bspline_curve_->basis().knotInterval(param_pnts[kr]);
	  vector<double>::iterator bspl_it = bspline_curve_->rcoefs_begin();
	  for (int i = kleft - kk_ + 1, j = 0; i <= kleft; ++i, ++j)
	   {
	     bspl_it[i*2] = 1.0;
	     bspline_curve_->point(p, param_pnts[kr]);
	     sb[j] = p[0];
	     bspl_it[i*2] = 0.0;
	   }
	}
      else
	{
	  qcurve_->basis().computeBasisValues(param_pnts[kr], sb, 0);
	  kleft = qcurve_->basis().lastKnotInterval();
	}

      for (k1=kleft-kk_+1, k2=0; k1<=kleft; k1++, k2++)
	{
	  // Test if the B-spline basis is too close to the boundary, i.e.
	  // the coefficient is set by the boundary conditions.           

	  if (coefknown_[k1] == 1)
	    continue;

// 	  if (k1 < cont_bound[0] || kn_-k1-1 < cont_bound[1])
// 	    continue;

// 	  kl1 = k1-cont_bound[0];
	  kl1 = pivot_[k1];

	  tz = wgt*pnt_weights[kr]*sb[k2];

	  // Add contribution to right hand side.

	  for (kj=0; kj<idim_; kj++)
	      //gright_.element(kj*kcond_+kl1) += tz*pnts[kr*idim_+kj];
	      gright_[kj*kcond_+kl1] += tz*pnts[kr*idim_+kj];

	  for (k3=kleft-kk_+1, k4=0; k3<=kleft; k3++, k4++)
	    {
	      // Test if the B-spline basis is too close to the boundary.

// 	      if (k3 < cont_bound[0] || kn-k3-1 < cont_bound[1])
	      if (coefknown_[k3] == 1)
		{
		  // Adjust for boundary conditions. 

		  sc = &*scoef_ + k3*kdim_;
		  for (kj=0; kj<idim_; kj++)
		    //gright_.element(kj*kcond_+kl1) -= sc[kj]*tz*sb[k4];
		    gright_[kj*kcond_+kl1] -= sc[kj]*tz*sb[k4];
		}
	      else
		{
		  // The term gives a contribution on the left hand side. 

// 		  kl2 = k3-cont_bound[0];
		  kl2 = pivot_[k3];
		  if (kl2 > kl1)
		    continue;

		  //gmat_.element(kl1, kl2) += tz*sb[k4];
		  gmat_[kl1][kl2] += tz*sb[k4];
		  if (kl2 < kl1)
		      //gmat_.element(kl2, kl1) += tz*sb[k4];
		      gmat_[kl2][kl1] += tz*sb[k4];
		}
	    }
	}
    }

  return;
}


/****************************************************************************/

void
SmoothCurve::setPeriodicity(int cont, double weight)
  //--------------------------------------------------------------------------
  //     Purpose : Set periodicity constraints
  //
  //     Calls   : 
  //
  //     Written by : Vibeke Skytt,  SINTEF Oslo,  04.98.
  //--------------------------------------------------------------------------
{
    int nmbconstraint = cont + 1;

    //  if (coefknown_[0] > 0 && coefknown_[kn_-1] > 0)
//   if (cont_bound[0] > 0 && cont_bound[1] > 0)
//    return;                   // Everything fixed already

  // Store continuity requirement
  cont_seam_ = std::min(nmbconstraint, std::min(kk_,2)); // At most C1-continuity currently

  // Set constraints on C0-continuity
  if (nmbconstraint > 0 && coefknown_[0] == 0 && coefknown_[kn_-1] == 0)
    setC0AtSeam(weight);
        
  // Set constraints on C1-continuity
  if (nmbconstraint > 1)
    setC1AtSeam(weight);
  
  return;
}


/****************************************************************************/

void SmoothCurve::setSideConstraints(std::vector<sideConstraint>& constraints)
//--------------------------------------------------------------------------
//     Purpose : Set linear side constraints to the minimization problem,
//               using Lagrange multiplier.
//
//     Calls   :
//
//     Written by : Sverre Briseid,  SINTEF Oslo,  06.02.
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

    // We check whether constraints are achievable. If coef in constraint is
    // already known, add to const term on right side. If all coefs are removed,
    // delete constraint. Issue warning if constraint is impossible.
    for (size_t i = 0; i < constraints.size(); ++i) {
      int dim = constraints[i].dim_;
      ALWAYS_ERROR_IF(dim != idim_,
		      "Constraint not member of the same geometric space as curve!");

      for (size_t j = 0; j < constraints[i].factor_.size(); ++j) {
	int coef_ind = constraints[i].factor_[j].first;
	if (coefknown_[coef_ind] != 0) {
	  // As coef was already known, we add to const term on right side.
	  for (int k = 0; k < dim; ++k)
	    constraints[i].constant_term_[k] -=
	      scoef_[kdim_*coef_ind+k]*constraints[i].factor_[j].second;
	  constraints[i].factor_.erase(constraints[i].factor_.begin() + j);
	  --j;
	}
      }
      // If all coefs in constraint were known, we remove constraint.
      if (constraints[i].factor_.size() == 0) {
	double sum = 0.0;
	for (int j = 0; j < dim; ++j)
	  sum += constraints[i].constant_term_[j];
	if (sum != 0.0) // We perform check whether constraint was fulfilled.
	  MESSAGE("All coefs already known, side constraint impossible.");
	constraints.erase(constraints.begin() + i); // We remove constraint.
	--i;
      }
    }

    // Make sure that dimension of equation system is valid.
  if (int(constraints.size()) != kconstraint_) {
      int new_kconstraint = (int)constraints.size();
    int new_kcond = kcond_ - (kconstraint_ - new_kconstraint);
    // For ease of algorithm, we copy matrices to new matrices.
    //Matrix new_gmat(new_kcond, new_kcond);
    //ColumnVector new_gright(idim_*new_kcond);
    vector<vector<double> > new_gmat(new_kcond);
    for (int ii = 0; ii < new_kcond; new_gmat[ii++].resize(kcond_)) {}
    vector<double> new_gright(idim_*new_kcond);

    // @@sbr If we use normal conditions, these copies will be inadequate!
    //     for (i = 0; i < kdim; ++i) // We treat one dimension at the time.
    for (int i = 0; i < new_kcond; ++i)
      for (int j = 0; j < new_kcond; ++j)
	  //new_gmat.element(i, j) = gmat_.element(i, j); @@sbr - possible bug here 
	  new_gmat[i][j] = gmat_[i][j]; // @@sbr - what if nb. of constr. is lower?
    for (int i = 0; i < idim_*new_kcond; ++i)
	//new_gright.element(i) = gright_.element(i);
	new_gright[i] = gright_[i];

    //gmat_ = new_gmat;
    //gright_ = new_gright;
    gmat_.swap(new_gmat);
    gright_.swap(new_gright);
    kconstraint_ = new_kconstraint;
    kcond_ = new_kcond;
  }

  // We update values in matries as described in the above equations.
  int nmb_free_coefs = kcond_ - kconstraint_; // Dimension of vector x (# coefs).
  // We start by updating gmat_ by adding new elements given by side constraints.
  for (size_t i = 0; i < constraints.size(); ++i)
    for (size_t j = 0; j < constraints[i].factor_.size(); ++j) { // We start with gmat_
      // We have made  sure that all elements in constraints[i] are free.
//       gmat_.element(nmb_free_coefs + i, pivot_[constraints[i].factor_[j].first]) =
// 	constraints[i].factor_[j].second;
//       gmat_.element(pivot_[constraints[i].factor_[j].first], nmb_free_coefs+i) =
// 	constraints[i].factor_[j].second;
      gmat_[nmb_free_coefs + i][pivot_[constraints[i].factor_[j].first]] =
	  constraints[i].factor_[j].second;
      gmat_[pivot_[constraints[i].factor_[j].first]][nmb_free_coefs+i] =
	  constraints[i].factor_[j].second;
    }

  // We next update gright_ by adding const values given by side constraints.
  for (int i = 0; i < idim_; ++i)
    for (int j = 0; j < kconstraint_; ++j)
//       gright_.element(i*kcond_+nmb_free_coefs+j) =
// 	constraints[j].constant_term_[i];
	gright_[i*kcond_+nmb_free_coefs+j] = constraints[j].constant_term_[i];
}


/****************************************************************************/

void
SmoothCurve::equationSolve(shared_ptr<SplineCurve>& curve)  // Resulting curve.
  //--------------------------------------------------------------------------
  //     Purpose : Solve linear equation system.
  //
  //     Calls   : s6lufacp, s6lusolp
  //
  //     Written by : Vibeke Skytt,  SINTEF Oslo,  10.95.
  //--------------------------------------------------------------------------
{
  int ki, kr;
  std::vector<int> nlvec(kcond_);   // Pivot_ array used in equation solver.  

  // Solve the equation system.  
  vector<vector<double> > mat_LU = gmat_;
  vector<int> permutation(gmat_.size());
  bool parity;
  LUDecomp(mat_LU, (int)gmat_.size(), &permutation[0], parity);

  //vector<double> sol(kcond_);
  vector<double> sol(gmat_.size());

  for (kr=0; kr<idim_; kr++) {
      for (ki = 0; ki < kcond_; ++ki) {
	  sol[ki] = gright_[kr * kcond_ + permutation[ki]];
      }
      forwardSubstitution(mat_LU, &sol[0], (int)gmat_.size());
      backwardSubstitution(mat_LU,&sol[0], (int)gmat_.size());

      for (ki = 0; ki < kn_; ++ki) {
	  if (coefknown_[ki] == 0) {
	      scoef_[ki*kdim_+kr] = sol[pivot_[ki]];
	  }
      }
  }

  // Make sure that continuity requirements at a seam are satisfied.
  if (cont_seam_ > 0)
    adjustAtSeam();

  // Create curve.
  curve = shared_ptr<SplineCurve>(new SplineCurve(kn_, kk_, st_, scoef_, idim_, rational_));
  return;
}

// --------below code is newmat dependent, so it has been replaced-------
// void
// SmoothCurve::equationSolve(shared_ptr<SplineCurve>& curve)  // Resulting curve.
//   //--------------------------------------------------------------------------
//   //     Purpose : Solve linear equation system.
//   //
//   //     Calls   : s6lufacp, s6lusolp
//   //
//   //     Written by : Vibeke Skytt,  SINTEF Oslo,  10.95.
//   //--------------------------------------------------------------------------
// {
//   int ki, kr;
//   std::vector<int> nlvec(kcond_);   // Pivot_ array used in equation solver.  

//   // Solve the equation system.  
//   CroutMatrix ALUfact = gmat_;
//   ALWAYS_ERROR_IF(ALUfact.IsSingular(), 
// 		  "Matrix is singular! This should never happen!");


//   for (kr=0; kr<idim_; kr++)
//     {
//       ColumnVector sb(kcond_);
//       ColumnVector sx(kcond_);
//       for (ki=0; ki<kcond_; ki++)
//         sb.element(ki) = gright_.element(kr*kcond_+ki);
//       sx = ALUfact.i()*sb;

// //       for (ki=cont_bound[0]; ki<kn-cont_bound[1]; ki++)
//       for (ki = 0; ki < kn_; ++ki)
// 	if (coefknown_[ki] == 0)
// 	  scoef_[ki*kdim_+kr] = sx.element(pivot_[ki]);
// // 	  scoef_[ki*kdim_+kr] = sx.element(ki-cont_bound[0]);
//     }

//   // Make sure that continuity requirements at a seam are satisfied.
//   if (cont_seam_ > 0)
//     adjustAtSeam();

//   // Create curve.

//   curve = shared_ptr<SplineCurve>(new SplineCurve(kn_, kk_, st_, scoef_, idim_, rational_));

//   return;
// }


void SmoothCurve::adjustAtSeam()
//--------------------------------------------------------------------------
//     Purpose : Ensure that the expected continuity at the seam
//               is satisfied.
//
//     Calls   : 
//
//     Written by : Vibeke Skytt,  SINTEF Oslo,  04.99.
//--------------------------------------------------------------------------
{
  //MESSAGE("Not tested after quick port...");

//   double tmean;
//   double t1 = (cont_bound[0]>0) ? 1.0 : 
//     ((cont_bound[1]<=0) ? 0.5 : 0.0);
//   //       double t1 = (cont_bound[0]>0) ? 1.0 : 
//   // 	((cont_bound[1]<=0) ? 0.5 : 0.0);
//   double t2 = 1.0 - t1;

//   for (kr=0; kr<idim_; kr++) {
//       tmean = t1*scoef_[kr] + t2*scoef_[(kn_-1)*kdim_+kr];
//       scoef_[kr] = scoef_[(kn_-1)*kdim_+kr] = tmean;
//   }

//   if (cont_seam_ > 1) {
//     t1 = (cont_bound[0]>1) ? 1.0 : ((cont_bound[1]<=1) ? 0.5 : 0.0);
//     t2 = 1.0 - t1;
//     double del1 = 1.0/(st_[kk_] - st_[1]);
//     double del2 = 1.0/(st_[kn_+kk_-2] - st_[kn_-1]);

//     for (kr=0; kr<idim_; kr++)
//       {
// 	tmean = t1*(scoef_[kdim_+kr] - scoef_[kr])*del1 +
// 	  t2*(scoef_[(kn_-1)*kdim_+kr] - scoef_[(kn_-2)*kdim_+kr])*del2;
// 	scoef_[kdim_+kr] = (tmean+del1*scoef_[kr])/del1;
// 	scoef_[(kn_-2)*kdim_+kr] = 
// 	  (del2*scoef_[(kn_-1)*kdim_+kr]-tmean)/del2;
//       }
//   }
  return;
}

/****************************************************************************/

void
SmoothCurve::setC0AtSeam(double weight)
  //--------------------------------------------------------------------------
  //     Purpose : Set C0-continuity constraints.
  //
  //     Calls   : 
  //
  //     Written by : Vibeke Skytt,  SINTEF Oslo,  04.98.
  //--------------------------------------------------------------------------
{
  // Adding term weight*(C(left_param) - C(right_param))^2

  bool left_free = coefknown_[0] == 0;
  bool right_free = coefknown_[kn_-1] == 0;

  if (!left_free && !right_free)
    return;

  double inv_w_left, inv_w_right;
  if (rational_)
    {
      inv_w_left = 1.0 / scoef_[idim_];
      inv_w_right = 1.0 / scoef_[kn_ * kdim_ - 1];
    }
  else
    inv_w_left = inv_w_right = 1.0;

  int piv_left = pivot_[0];
  int piv_right = pivot_[kn_-1];

  if (left_free) gmat_[piv_left][piv_left] += weight*inv_w_left*inv_w_left;
  if (right_free) gmat_[piv_right][piv_right] += weight*inv_w_right*inv_w_right;

  if (left_free && right_free)
    {
      gmat_[piv_left][piv_right] -= weight*inv_w_left*inv_w_right;
      gmat_[piv_right][piv_left] -= weight*inv_w_left*inv_w_right;
    }

  if (!left_free)
    for (int i = 0; i < idim_; ++i)
      gright_[i*kcond_ + piv_right] += weight*inv_w_left*inv_w_right*scoef_[i];
  if (!right_free)
    for (int i = 0; i < idim_; ++i)
      gright_[i*kcond_ + piv_left] += weight*inv_w_left*inv_w_right*scoef_[(kn_-1)*kdim_ + i];
}

/****************************************************************************/

void
SmoothCurve::setC1AtSeam(double weight)
  //--------------------------------------------------------------------------
  //     Purpose : Set C1-continuity constraints.
  //
  //     Calls   : 
  //
  //     Written by : Vibeke Skytt,  SINTEF Oslo,  04.98.
  //--------------------------------------------------------------------------
{
  // Adding term weight*(C'(left_param) - C'(right_param))^2

  int coef_pos[4];
  bool is_free[4];
  int piv_pos[4];
  double term_coefs[4];

  coef_pos[0] = 0;
  coef_pos[1] = 1;
  coef_pos[2] = kn_-2;
  coef_pos[3] = kn_-1;

  bool any_free = false;
  for (int i = 0; i < 4; ++i)
    any_free |= (is_free[i] = (coefknown_[coef_pos[i]] == 0));
  if (!any_free)
    return;

  for (int i = 0; i < 4; ++i)
    piv_pos[i] = pivot_[coef_pos[i]];

  double del[2];
  del[0] = double(kk_-1)/(st_[kk_] - st_[1]);
  del[1] = double(kk_-1)/(st_[kn_+kk_-2] - st_[kn_-1]);

  term_coefs[0] = -del[0];
  term_coefs[1] = del[0];
  term_coefs[2] = del[1];
  term_coefs[3] = -del[1];

  if (rational_)
    {
      double wl0 = scoef_[idim_];
      double wl1 = scoef_[idim_+kdim_];
      double wr1 = scoef_[(kn_-1) * kdim_ - 1];
      double wr0 = scoef_[kn_ * kdim_ - 1];
      term_coefs[0] *= wl1 / (wl0*wl0);
      term_coefs[1] /= wl0;
      term_coefs[2] /= wr0;
      term_coefs[3] *= wr1 / (wr0*wr0);
    }

  for (int i = 0; i < 4; ++i)
    if (is_free[i])
      {
	double term_i = weight * term_coefs[i];
	for (int j = 0; j < 4; ++j)
	  {
	    if (is_free[j])
	      gmat_[piv_pos[i]][piv_pos[j]] += term_i * term_coefs[j];
	    else
	      for (int k = 0; k < idim_; ++k)
		gright_[k*kcond_ + piv_pos[i]] -= term_i * term_coefs[j] * scoef_[coef_pos[j]*kdim_ + k];
	  }
      }
}
