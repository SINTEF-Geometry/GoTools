/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrCG.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Dec 98
 DESCRIPTION : Implementation of methods in the class PrCG.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrCG.h"
#include "GoTools/utils/CPUclock.h"

using namespace Go;

//-----------------------------------------------------------------------------
PrCG::PrCG()
//-----------------------------------------------------------------------------
{
  tolerance_ = 1.0e-6;
  max_iterations_ = 0;
  it_count_ = 0;
  cpu_time_ = 0.0;
  converged_ = 0;
}


//-----------------------------------------------------------------------------
void PrCG::solve(const PrMatrix& A, PrVec& x, const PrVec& b)
//-----------------------------------------------------------------------------
{
  CPUclock rolex;
  double time0 = rolex.getTime();
  double tol = tolerance_ * tolerance_;

  int n = x.size();
  int j;

  PrVec r(n);
  //r = b - Ax
  A.prod(x,r);
  for(j=0; j<n; j++) r(j) = b(j) - r(j);

  PrVec p(n);
  for(j=0; j<n; j++) p(j) = r(j);
  double alpha, beta, rnorm, rnorm2;
  rnorm = r.inner(r);

  if(rnorm < tol)
  {
    it_count_ = 0;
    cpu_time_ = 0.0;
    converged_ = true;
    return;
  }

  PrVec rhat(n);
  PrVec q(n);

  for(int i=1; i<= max_iterations_; i++)
  {
    //s_o << "i = " << i << endl;

    A.prod(p,q);
    alpha = rnorm / (p.inner(q));

    //r := r - alpha * A p
    for(j=0; j<n; j++) r(j) -= alpha * q(j);

    //x := x + alpha p
    for(j=0; j<n; j++) x(j) += alpha * p(j);

    rnorm2 = r.inner(r);
    beta = rnorm2 / rnorm;

    //p = r + beta * p
    for(j=0; j<n; j++) p(j) = r(j) + beta * p(j);

    if(rnorm2 < tol)
    {
      it_count_ = i;
      cpu_time_ = rolex.getTime() - time0;
      converged_ = true;
      return;
    }

    rnorm = rnorm2;
  }

  it_count_ = max_iterations_;
  cpu_time_ = rolex.getTime() - time0;
  converged_ = false;

}

