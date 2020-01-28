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

#include "GoTools/creators/SolveCG.h"
#include "GoTools/utils/errormacros.h"

#include <stdio.h>
#include <math.h>
#include <iostream>


using namespace Go;


// afr: I added this function to avoid using s6scpr from SISL
namespace {
  inline double scalar_product(double* v1, double* v2, int n)
  {
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
      res += v1[i]*v2[i];
    }
    return res;
  }
}

SolveCG::SolveCG()
//--------------------------------------------------------------------------
//     Constructor for class SolveCG
//
//     Purpose : Initialize class variables
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF, 09.99
//--------------------------------------------------------------------------
{
  nn_ = np_ = 0;
  tolerance_ = 1.0e-6;
  max_iterations_ = 0;
  diagset_ = 0;
}

/****************************************************************************/

SolveCG::~SolveCG()
//--------------------------------------------------------------------------
//     Destructor for class SolveCG
//
//     Purpose :
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF, 09.99
//--------------------------------------------------------------------------
{
}

/****************************************************************************/

int SolveCG::getIndex(int ki, int kj)
//--------------------------------------------------------------------------
//
//     Purpose :  Given an index in the full equation system, get the
//                index in A_ (and M_).
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF, 10.99
//--------------------------------------------------------------------------
{
  if (ki == kj && diagset_)
    return diagonal_[ki];

// find index by the bisection method, the algorithm assumes that the
// column indices in a row are given in increasing order

  int jl  = irow_[ki];
  int ju  = irow_[ki+1];

  int jm,jam;
  while( (ju - jl) >= 1 )
    {
      jm = (ju+jl) >> 1;
      jam = jcol_[jm];
      if( kj == jam )
	return jm;

      if( kj > jam )
	jl = jm+1;
      else
	ju = jm;
  }
  // If the loop is terminated the index kj  was not found:
  return -1;
}


/****************************************************************************/

void SolveCG::attachMatrix(double *gmat, int nn)
//--------------------------------------------------------------------------
//
//     Purpose : Attach the left side of the equation system to the current
//               object and represent the matrix as a sparse matrix. No test
//               is applied on whether the matrix really is symmetric and
//               positive definite.
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF, 09.99
//--------------------------------------------------------------------------
{
  nn_ = nn;

  // Count the number of non-zero elements in the input matrix.

  int ki, kj, idx;
  np_ = 0;
  for (ki=0; ki<nn*nn; ki++)
    if (gmat[ki] != 0.0)
      np_++;

  // Reserve the required scratch for the matrix arrays.

  A_.reserve(np_);
  jcol_.reserve(np_);
  irow_.reserve(nn_ + 1);

  // Fill in the non-zero elements of the input matrix.

  for (idx=0, kj=0; kj<nn; kj++)
    {
      irow_.push_back(idx);
      for (ki=0; ki<nn; ki++)
	if (gmat[kj*nn+ki] != 0.0)
	  {
	    A_.push_back(gmat[kj*nn+ki]);
	    jcol_.push_back(ki);
	    idx++;
	  }
    }
  if (irow_.size() > 0 && irow_[irow_.size()-1] == idx)
    THROW("Singular equation system");
  irow_.push_back(idx);
}

/****************************************************************************/

void SolveCG::precondRILU(double relaxfac)
//--------------------------------------------------------------------------
//
//     Purpose : Prepare for preconditioning.
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF, 10.99
//--------------------------------------------------------------------------
{
    omega_ = relaxfac;

    // PrecondRILU() assumes diagonal elements of matrix are non-zero.

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
    for (kr=0; kr<nn1; kr++) {
	rr = getIndex(kr, kr);
#ifdef DEBUG
	if (rr < 0)
	  std::cout << "SolveCG. Error in left hand side matrix" << std::endl;
#endif
	diag = M_[rr];
	kstop = irow_[kr+1];
	for (k1=rr+1; k1<kstop; k1++) {
	    ki = jcol_[k1];
	    ir = getIndex(ki, kr);

	    if (ir < 0)
		continue;   // Not a symmetric matrix. Unpredictable result.

	    if (M_[ir] != 0.0) {
		elem = M_[ir]/diag;
		M_[ir] = elem;
		ii = getIndex(ki, ki);
		for (k2=rr+1; k2<kstop; k2++) {
		    kj = jcol_[k2];
		    if (M_[k2] != 0.0) {
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

void SolveCG::forwBack(double *r, double *s)
//--------------------------------------------------------------------------
//
//     Purpose : Solve the equation system M_*s = r, where M_ stores an
//               LU-factorized matrix. Forward - backward substitution
//               is used.
//
//     Calls   :
//
//     Written by : Vibeke Skytt,  SINTEF, 10.99
//--------------------------------------------------------------------------
{
  int ki, kj, kd, kstop;
  double tmp;

  for (ki=0; ki<nn_; ki++)
    s[ki] = r[ki];
  for (ki=0; ki<nn_; ki++)
    {
      tmp = 0.0;
      for (kj=irow_[ki]; jcol_[kj]<ki; kj++)
	tmp += M_[kj]*s[jcol_[kj]];

      s[ki] -= tmp;
    }

  s[nn_-1] /= M_[getIndex(nn_-1, nn_-1)];
  for (ki=nn_-2; ki>=0; ki--)
    {
      tmp = 0.0;
      kstop = irow_[ki+1];
      kd = getIndex(ki, ki);
#ifdef DEBUG
      if (kd < 0)
	std::cout << "SolveCG. Error in left hand side matrix" << std::endl;
#endif
      for (kj=kd+1; kj<kstop; kj++)
	tmp += M_[kj]*s[jcol_[kj]];

      s[ki] = (s[ki] - tmp)/M_[kd];
    }
}


/****************************************************************************/

void SolveCG::transposedMatrixProduct(double *sx, double *sy)
//--------------------------------------------------------------------------
//
//     Purpose : Compute sy = A_^T * sx.
//
//     Calls   :
//
//     Written by : Sverre Briseid,  SINTEF, Sep 2004
//--------------------------------------------------------------------------
{
  int kj, ki;
  for (kj = 0; kj < nn_; ++kj)
    sy[kj] = 0.0;
  for(kj=0; kj<nn_; kj++)
  {
    for(ki=irow_[kj]; ki<irow_[kj+1]; ki++)
    {
      sy[jcol_[ki]] += A_[ki] * sx[kj];
    }
  }
}


/****************************************************************************/

int SolveCG::solve(double *x, double *b, int nn)
//--------------------------------------------------------------------------
//
//     Purpose : Solve the equation system by conjugate gradient method.
//
//     Input   : x   -  Guess on the unknowns.
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
//     Written by : Vibeke Skytt,  SINTEF, 09.99
//--------------------------------------------------------------------------
{
  if (M_.size() > 0)
    return solveRILU(x, b, nn);
  else
    return solveStd(x, b, nn);
}


/****************************************************************************/

int SolveCG::solveStd(double *x, double *b, int nn)
//--------------------------------------------------------------------------
//
//     Purpose : Solve the equation system by conjugate gradient method.
//
//     Input   : x   -  Guess on the unknowns.
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
//     Written by : Vibeke Skytt,  SINTEF, 09.99
//--------------------------------------------------------------------------
{
  double tol = nn * tolerance_ * tolerance_;

  if (nn != nn_)
    return -106;   // Conflicting dimensions of equation system.

  int kj;

  std::vector<double> r(nn, 0.0);
  //r = b - Ax
  matrixProduct(x, r.begin());
  for(kj=0; kj<nn; kj++)
    r[kj] = b[kj] - r[kj];

  std::vector<double> p(nn, 0.0);
  for(kj=0; kj<nn; kj++)
    p[kj] = r[kj];
  double alpha, beta, rnorm, rnorm2, rnorm0;
  rnorm0 = rnorm = scalar_product(&r[0], &r[0], nn);

  if (fabs(rnorm) < tol)
    return 0;

  std::vector<double> q(nn, 0.0);

  for (int ki=0; ki< max_iterations_; ki++)
  {
    matrixProduct(p.begin(), q.begin());
    alpha = rnorm / scalar_product(&p[0], &q[0], nn);

    //r := r - alpha * A p
    for(kj=0; kj<nn; kj++)
      r[kj] -= alpha * q[kj];

    rnorm2 = scalar_product(&r[0], &r[0], nn);
    beta = rnorm2 / rnorm;

    //x := x + alpha p
    for(kj=0; kj<nn; kj++)
      x[kj] += alpha * p[kj];

    //p = r + beta * p
    for(kj=0; kj<nn; kj++)
      p[kj] = r[kj] + beta * p[kj];

    if (rnorm2 < tol && rnorm2/rnorm0 < tolerance_)
      {
// 	printf("No. of iterations %d, error %7.13f %7.13f \n", ki, rnorm2,
// 	       rnorm2/rnorm0);
	return 0;
      }

//     if (ki % 500 == 0)
// 	printf("No. of iterations %d, error %7.13f %7.13f \n", ki, rnorm2,
// 	       rnorm2/rnorm0);

    rnorm = rnorm2;
  }

//   printf("No. of iterations %d, error %7.13f %7.13f \n", ki, rnorm2,
// 	 rnorm2/rnorm0);
  return 1;
}

/****************************************************************************/

int SolveCG::solveRILU(double *x, double *b, int nn)
//--------------------------------------------------------------------------
//
//     Purpose : Solve the equation system by conjugate gradient method
//               using an already given RILU-preconditioner.
//
//     Input   : x   -  Guess on the unknowns.
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
//     Written by : Vibeke Skytt,  SINTEF, 10.99
//--------------------------------------------------------------------------
{
  double tol = nn * tolerance_ * tolerance_;

  if (nn != nn_)
    return -106;   // Conflicting dimensions of equation system.

  int kj;

  std::vector<double> r(nn, 0.0);
  //r = b - Ax
  matrixProduct(x, r.begin());
  for(kj=0; kj<nn; kj++)
    r[kj] = b[kj] - r[kj];

  std::vector<double> p(nn, 0.0);
  forwBack(&r[0], &p[0]);

  double alpha, beta, rnorm, rnorm2, rnorm0;
  rnorm0 = rnorm = scalar_product(&p[0], &r[0], nn);

  if (fabs(rnorm) < tol)
    return 0;

  std::vector<double> q(nn, 0.0);
  std::vector<double> s(nn, 0.0);

  for (int ki=0; ki< max_iterations_; ki++)
  {
    matrixProduct(p.begin(), q.begin());
    alpha = rnorm / scalar_product(&p[0], &q[0], nn);

    //r := r - alpha * A p
    for(kj=0; kj<nn; kj++)
      r[kj] -= alpha * q[kj];

    //x := x + alpha p
    for(kj=0; kj<nn; kj++)
      x[kj] += alpha * p[kj];

    forwBack(&r[0], &s[0]);

    rnorm2 = scalar_product(&s[0], &r[0], nn);
    beta = rnorm2 / rnorm;

    //p = r + beta * p
    for(kj=0; kj<nn; kj++)
      p[kj] = s[kj] + beta * p[kj];

//     if (rnorm2 < 0.0)
// 	printf("No. of iterations %d, error %7.13f %7.13f \n", ki, rnorm2,
// 	       rnorm2/rnorm0);

    if (fabs(rnorm2) < tol && fabs(rnorm2/rnorm0) < tolerance_)
      {
// 	printf("No. of iterations %d, error %7.13f %7.13f \n", ki, rnorm2,
// 	       rnorm2/rnorm0);
	return 0;
      }

//     if (ki % 100 == 0)
// 	printf("No. of iterations %d, error %7.13f %7.13f \n", ki, rnorm2,
// 	       rnorm2/rnorm0);

    rnorm = rnorm2;
  }

//   printf("No. of iterations %d, error %7.13f %7.13f \n", ki, rnorm2,
// 	 rnorm2/rnorm0);
  return 1;
}


void SolveCG::printPrecond()
{
  FILE* fp = NULL;
  fp = fopen("fM.m", "w");
  fprintf(fp," M=[ ");
  int i,j,ij;
  for (i=0; i< nn_; i++)
  {
    for (j=0; j< nn_; j++)
    {
      ij=getIndex(i,j);
      if (ij == -1)
	fprintf(fp," %lf",0.0);
      else
	fprintf(fp," %lf",M_[ij]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"];\n");

  fclose(fp);
}
