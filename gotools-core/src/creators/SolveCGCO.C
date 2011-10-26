

#include "GoTools/creators/SolveCGCO.h"
#include <math.h>

using namespace Go;

SolveCGCO::SolveCGCO(int m, int n)
  : SolveCG(), m_(m), n_(n)
{
}


SolveCGCO::~SolveCGCO()
{
}

/****************************************************************************/

void SolveCGCO::precondRILU(double relaxfac)
//--------------------------------------------------------------------------
//
//     Purpose : Prepare for preconditioning.
//
//     Calls   : Setting preconditioner for the constrained optimization
//               problem.
//
//     Written by : Vibeke Skytt,  SINTEF, 10.99
//--------------------------------------------------------------------------
{
  omega_ = relaxfac;

  // Allocate storage for the preconditioning matrix.

  // We're using a diagonal block structure:
  //       ( M  0 )
  // M_ =  (      ),
  //       ( 0  N )
  // where M is the preconditioner for the unknown coefs matrix A (as in
  // SolveCG) and N is a diagonal matrix corresponding to the
  // constraints.  We want N to be appr equal to
  // (CA^{-1}C^T)^{-1}.  But currently we settle for the diagonal
  // matrix, with a suitable scaling.

  M_.reserve(np_+n_); // The last n_ elements are used for the
		      // diagonal elements of the lower right block.
  int kr, kp;
  for (kr=0; kr<np_; kr++) {
    M_.push_back(0.0);
  }
  int ki, kj;
  // We start by constructing the preconditioning matrix for the
  // matrix A (avoid including elements in A_ corresponding
  // to the constraints).
  for (ki = 0; ki < m_; ++ki)
    for (kj = 0; kj < m_; ++kj)
      {
	int ind = getIndex(ki, kj);
	if (ind > -1)
	  M_[ind] = A_[ind];
      }

  // Create vector of indexes along the diagonal of A_ and M_.
  diagonal_.reserve(nn_);
  for (kr=0; kr<nn_; kr++)
    diagonal_.push_back(getIndex(kr, kr));
  diagset_ = 1;

  // Factorize the M_ matrix.
  // nn_ is the size of the system.
  // We should make sure that all entries in M_ deriving from lower left and
  // upper right blocks are set to 0.0.
  int k1, k2;
  int rr, ir, ii, ij;
  int kstop;
  double diag, elem;
  for (kr=0; kr<nn_ - 1; kr++)
    {
      rr = getIndex(kr, kr);
      if (rr < 0)
	continue; // The diagonal in lower right matrix consists of zeros.
      diag = M_[rr];
      kstop = irow_[kr+1];
      for (k1=rr+1; k1<kstop; k1++)
	{
	  ki = jcol_[k1];
	  ir = getIndex(ki, kr);

	  if (ir < 0)
	    continue;   // Not a symmetric matrix. Unpredictable result.

	  if (fabs(M_[ir])> 1.0e-12)
	    {
	      elem = M_[ir]/diag;
	      M_[ir] = elem;
	      ii = getIndex(ki, ki);
	      for (k2=rr+1; k2<kstop; k2++)
		{
		  kj = jcol_[k2];
		  if (fabs(M_[k2])>1.0e-12)
		    {
		      ij = getIndex(ki, kj);
		      if (ij >= 0)
			M_[ij] -= elem*M_[k2];
		      else
			M_[ii] -= elem*omega_*M_[k2];
		    }
		}
	    }
	  else
	    M_[ir] = 0.0;
	}
    }

  // We then set elements in lower right block.  As we want the matrix
  // N to be appr the inverse of the above mentioned product, we must
  // consider the values in the original problem when choosing our
  // diagonal elments.  We scale the diagonal matrix N using the
  // sup-norm of the product elements (which only yields an upper limit).
//   double l2_sum_A = 0.0;
//   double l2_sum_C = 0.0;
//   for (ki = 0; ki < irow_.size() - 1; ++ki)
//     {
//       for (kj = irow_[ki]; kj < irow_[ki+1]; ++kj)
// 	if (ki < m_ && jcol_[kj] < m_)
// 	  l2_sum_A += A_[kj]*A_[kj];
// 	else
// 	  l2_sum_C += A_[kj]*A_[kj];
//     }
  double sup_A = 0.0;
  double sup_C = 0.0;
  for (ki = 0; ki < (int)irow_.size() - 1; ++ki)
    {
      for (kj = irow_[ki]; kj < irow_[ki+1]; ++kj)
	if (ki < m_ && jcol_[kj] < m_)
	  {
	    if (fabs(A_[kj]) > sup_A)
	      sup_A = fabs(A_[kj]);
	  }
	else
	  {
	    if (fabs(A_[kj]) > sup_C)
	      sup_C = fabs(A_[kj]);
	  }
    }
//  double appr_N_l2_norm = sqrt(l2_sum_A)/(l2_sum_C); // Upper limit.
//  double diag_scale = sqrt(n_)*appr_N_l2_norm;
  // The 5.0 scale necessary for 2 specific examples ...
  // Hmm, running system on a bunch of examples the identiy matrix yielded
  // the overall best result (i.e. the largest 
  double diag_scale = 1.0;//5.0*appr_N_sup_norm;
//   if (diag_scale < 1e-02)
//     diag_scale = 1e-02;
//   else if (diag_scale > 1e02) // Too large value seems to lead to bad
// 			      // accuracy for the constraints.
//     diag_scale =1e02;
  for (kr = 0; kr < n_; ++kr)
    { // We add elements to M_ (for the diagonal matrix).
      // The elements of M_ expected to be ordered!
      int ju = irow_[m_+kr+1];
      // As A_ and M_ are expected to share indexing we must add
      // diagonal zeros to A_.
      M_.insert(M_.begin() + ju, diag_scale);//1.0);
      A_.insert(A_.begin() + ju, 0.0);
      for (kp = kr; kp < n_; ++kp)
	irow_[m_+kp+1] += 1;
      jcol_.insert(jcol_.begin() + ju, m_ + kr);
      ju = irow_[m_+kr+1]; // New value for ju.
      diagonal_[m_+kr] = ju - 1;
    }

  np_ = (int)A_.size();

//   printPrecond();
}
