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

