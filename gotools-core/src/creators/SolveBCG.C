//===========================================================================
//                                                                           
// File: solveBCG.C                                                        
//                                                                           
// Created: Wed Sep 15 09:48:45 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: 
//                                                                           
// Description: 
//                                                                           
//===========================================================================


#include "GoTools/creators/SolveBCG.h"
#include "GoTools/geometry/Utils.h"

using std::vector;
using namespace Go;

SolveBCG::SolveBCG(int conv_type, bool symm)
  : SolveCG(), conv_type_(conv_type), symm_(symm)
{
  ASSERT(conv_type_ > 0 && conv_type_ < 5);
}


/****************************************************************************/

SolveBCG::~SolveBCG()
{
}


/****************************************************************************/

void SolveBCG::precond(double relaxfac)
//--------------------------------------------------------------------------
//
//     Purpose : Prepare for preconditioning.
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF, 10.99
//--------------------------------------------------------------------------
{
  precondRILU(relaxfac);
}


/****************************************************************************/

int SolveBCG::solve(double *x, double *b, int nn)
//--------------------------------------------------------------------------
//
//     Purpose : Solve the equation system by conjugate gradient method.
//
//     Input   : x   -  Guess on the unknowns (i.e. start pt fotr iteration).
//               b   -  Right side of the equation system.
//               nn  -  Number of unknowns.
//
//     Output  : solve - Status.
//                        1  -  No convergence within the given number
//                              of iterations.
//                        0  -  Equation system solved, OK.
//                     -106  -  Conflicting dimension of arrays.
//               x         - The solution to the equation system.
//
//     Calls   :
//
//     Written by : Sverre Briseid,  SINTEF, 09.04
//--------------------------------------------------------------------------
{
  int kstat =0;

  int ki;
  double inner_dir, inner_dir_frac, inner_res, inner_res_frac;
  double prev_inner_res = 1.0;  // Initialize for the compiler
  double x_l, der_x_l, b_l, prev_step_l, dir_l;
  double step_l = 0.0;
  vector<double> res(nn, 0.0); // Residual
  vector<double> cres(nn, 0.0); // The conjugate to res
  vector<double> dir(nn, 0.0); // Direction vector
  vector<double> cdir(nn, 0.0); // The conjugate to dir
  vector<double> step(nn, 0.0); // Step vector
  vector<double> cstep(nn, 0.0); // The conjugate to step
  // If the matrix is symmetric the conjugates should be equal to their counterparts.
  // If the matrix is also positive definite it cannot break down (in theory, w/inf prec).
  int iter = 0;
  double error = -1.0;
  const double epsilon = 1.0e-14;
  bool precond = (M_.size() > 0);
  double tol = sqrt((double)nn)*tolerance_; // Equivivalent to tolerance used in conjugate gradient method.

  vector<double> A_diag; // The diagonal of A, a simple preconditioner.
  A_diag.reserve(nn);
  for (ki = 0; ki < nn; ++ki)
    {
      int ind = getIndex(ki, ki);
      if (ind > -1)
	A_diag.push_back(A_[ind]);
      else
	A_diag.push_back(0.0);
    }

  // We must first initialize our initial residual.
  matrixProduct(x, &res[0]);
  for (ki = 0; ki < nn; ++ki)
    {
      res[ki] = b[ki] - res[ki];
      cres[ki] = res[ki];
    }

  if (symm_) // Then method specializes the minimal residual method. No longer guaranteed conv.
    matrixProduct(&res[0], &cres[0]);

  if (conv_type_ == 1)
    {
      b_l = sqrt(inner(b, b+nn, b));
      if (kstat < 0)
	{
	  return kstat;
	}
      if (precond)
	{ // Then we should use the preconditioner.
	  forwBack(&res[0], &step[0]);
	}
      else
	{
	  for (ki = 0; ki < nn; ++ki)
	    step[ki] = (A_diag[ki] != 0.0) ? res[ki]/A_diag[ki] : res[ki];
	}
    }
  else if (conv_type_ == 2)
    {
      if (precond)
	{ // Then we should use the preconditioner.
	  forwBack(b, &step[0]);
	}
      else
	{
	  for (ki = 0; ki < nn; ++ki)
	    step[ki] = (A_diag[ki] != 0.0) ? b[ki]/A_diag[ki] : b[ki];
	}
      b_l = sqrt(inner(&step[0], &step[0]+nn, &step[0]));
      if (kstat < 0)
	{
	  return kstat;
	}
      if (precond)
	{ // Then we should use the preconditioner.
	  forwBack(&res[0], &step[0]);
	}
      else
	{
	  for (ki = 0; ki < nn; ++ki)
	    step[ki] = (A_diag[ki] != 0.0) ? res[ki]/A_diag[ki] : res[ki];
	}
    }
  else if (conv_type_ == 3 || conv_type_ == 4)
    {
      if (precond)
	{ // Then we should use the preconditioner.
	  forwBack(b, &step[0]);
	}
      else
	{
	  for (ki = 0; ki < nn; ++ki)
	    step[ki] = (A_diag[ki] != 0.0) ? b[ki]/A_diag[ki] : b[ki];
	}
      b_l = (conv_type_ == 3) ?
	sqrt(inner(&step[0], &step[0]+nn, &step[0])) : max_abs_element(&step[0], nn);
      if (precond)
	{ // Then we should use the preconditioner.
	  forwBack(&res[0], &step[0]);
	}
      else
	{
	  for (ki = 0; ki < nn; ++ki)
	    step[ki] = (A_diag[ki] != 0.0) ? res[ki]/A_diag[ki] : res[ki];
	}
      b_l = (conv_type_ == 3) ?
	sqrt(inner(&step[0], &step[0]+nn, &step[0])) : max_abs_element(&step[0], nn);
    }
  else
    { // Unexpected conv_type_
      return -1;
    }

  // We're now ready for our main loop.
  while (iter < max_iterations_)
    {
      ++iter;
      if (precond)
	{ // Then we should use the preconditioner.
	  forwBack(&cres[0], &cstep[0]);
// 	  transposedForwBack(cres.begin(), cstep.begin()); @@sbr Implement!
	}
      else
	{
	  for (ki = 0; ki < nn; ++ki)
	    cstep[ki] = (A_diag[ki] != 0.0) ? cres[ki]/A_diag[ki] : cres[ki];
	}

      inner_res = 0.0;
      for (ki = 0; ki < nn; ++ki)
	{
	  inner_res += step[ki]*cres[ki];
	}

      // We then proceed to the direction vectors.
      if (iter == 1)
	{
	  for (ki = 0; ki < nn; ++ki)
	    {
	      dir[ki] = step[ki];
	      cdir[ki] = cstep[ki];
	    }
	}
      else
	{
	  inner_res_frac = inner_res/prev_inner_res;
	  for (ki = 0; ki < nn; ++ki)
	    {
	      dir[ki] = inner_res_frac*dir[ki] + step[ki];
	      cdir[ki] = inner_res_frac*cdir[ki] + cstep[ki];
	    }
	}

      // We then proceed to calculate our new residuals res & cres, and update
      // solution vector x.
      prev_inner_res = inner_res;
      matrixProduct(&dir[0], &step[0]);
      inner_dir = 0.0;
      for (ki = 0; ki < nn; ++ki)
	{
	  inner_dir += step[ki]*cdir[ki];
	}
      inner_dir_frac = inner_res/inner_dir;
      transposedMatrixProduct(&cdir[0], &cstep[0]);
      for (ki = 0; ki < nn; ++ki)
	{
	  x[ki] += inner_dir_frac*dir[ki];
	  res[ki] -= inner_dir_frac*step[ki];
	  cres[ki] -= inner_dir_frac*cstep[ki];
	}

      if (precond)
	{ // Then we should use the preconditioner.
	  forwBack(&res[0], &step[0]);
	}
      else
	{
	  for (ki = 0; ki < nn; ++ki)
	    step[ki] = (A_diag[ki] != 0.0) ? res[ki]/A_diag[ki] : res[ki];
	}

      // We then perform our error estimate to see if we are done.
      if (conv_type_ == 1)
	{
	  error = sqrt(inner(&res[0], &res[0]+nn, &res[0]))/b_l;
	  if (kstat < 0)
	    {
	      return kstat;
	    }
	}
      else if (conv_type_ == 2)
	{
	  error = sqrt(inner(&step[0], &step[0]+nn, &step[0]))/b_l;
	  if (kstat < 0)
	    {
	      return kstat;
	    }
	}
      else if (conv_type_ == 3 || conv_type_ == 4)
	{
	  prev_step_l = step_l;
	  step_l = (conv_type_ == 3) ?
	    sqrt(inner(&step[0], &step[0]+nn, &step[0])) : 
	    max_abs_element(&step[0], nn);
	  if (fabs(prev_step_l - step_l) > epsilon*step_l)
	    {
	      dir_l = (conv_type_ == 3) ?
		sqrt(inner(&dir[0], &step[0]+nn, &dir[0])) : 
		max_abs_element(&dir[0], nn);
	      if (kstat < 0)
		{
		  return kstat;
		}
	      der_x_l = fabs(inner_dir_frac)*dir_l;
	      error = step_l/(fabs(prev_step_l - step_l)*der_x_l);
	    }
	  else
	    { // Suspecting estimate not quite accurate, one more iteration.
	      error = step_l/b_l;
// 	      continue;
	    }

	  // Finally we compute our progress.
	  x_l = (conv_type_ == 3) ?
	    sqrt(inner(&x[0], &x[0]+nn, &x[0])) : max_abs_element(&x[0], nn);
	  if (error <= 0.5*x_l)
	    {
	      error /= x_l;
	    }
	  else
	    {
	      error = step_l/b_l;
// 	      continue;
	    }
	}

      if (error < tol)
	{
	  return 0; // We're done.
	}
    }



  return 1;
}


/****************************************************************************/

void SolveBCG::precondRILU(double relaxfac)
//--------------------------------------------------------------------------
//
//     Purpose : Prepare for preconditioning.
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF, 10.99
//--------------------------------------------------------------------------
{
    // @@sbr Currently using preconditioner for a symm AND indef system...

  omega_ = relaxfac;

  // Allocate storage for the preconditioning matrix.

  M_.reserve(np_);
  int kr;
  for (kr=0; kr<np_; kr++)
    M_.push_back(A_[kr]);

  // Create vector of indexes along the diagonal of A_ and M_.
  diagonal_.reserve(nn_);
  for (kr=0; kr<nn_; kr++)
    diagonal_.push_back(getIndex(kr, kr));
  diagset_ = 1;

  // Factorize the M_ matrix.

  int k1, k2, ki, kj;
  int rr, ir, ii, ij;
  int kstop;
  double diag, elem;
  int nn1 = nn_ - 1;
  for (kr=0; kr<nn1; kr++)
    {
      rr = getIndex(kr, kr);
      diag = M_[rr];
      kstop = irow_[kr+1];
      for (k1=rr+1; k1<kstop; k1++)
	{
	  ki = jcol_[k1];
	  ir = getIndex(ki, kr);

	  if (ir < 0)
	    continue;   // Not a symmetric matrix. Unpredictable result.

// 	  if ((HPSLpt_ DNEQUAL(M_[ir], 0.0)) && (diag != 0.0))
 	  if (fabs(M_[ir])>1.0e-12)
	    {
	      elem = M_[ir]/diag;
	      M_[ir] = elem;
	      ii = getIndex(ki, ki);
	      for (k2=rr+1; k2<kstop; k2++)
		{
		  kj = jcol_[k2];
		  if (fabs(M_[k2]) > 1.0e-12)
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
//  printPrecond();
}


/****************************************************************************/

double SolveBCG::max_abs_element(double* x, int nn)
{
  int ki;
  ASSERT(nn > 0);
  double max_elem = fabs(x[0]);
  for (ki = 1; ki < nn; ++ki)
    if (fabs(x[ki]) > max_elem)
      max_elem = fabs(x[ki]);

  return max_elem;
}
