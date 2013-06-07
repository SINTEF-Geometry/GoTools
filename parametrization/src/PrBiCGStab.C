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

#include "GoTools/parametrization/PrBiCGStab.h"
#include "GoTools/utils/timeutils.h"

//-----------------------------------------------------------------------------
PrBiCGStab::PrBiCGStab()
//-----------------------------------------------------------------------------
{
  tolerance_ = 1.0e-6;
  max_iterations_ = 0;
  it_count_ = 0;
  cpu_time_ = 0.0;
  converged_ = 0;
}

//-----------------------------------------------------------------------------
void PrBiCGStab::solve(const PrMatrix& A, PrVec& x, const PrVec& b)
//-----------------------------------------------------------------------------
{
  double time0 = Go::getCurrentTime();
  double tol = tolerance_ * tolerance_;

  int n = x.size();
  int j;

  PrVec r(n);
  //r = b - Ax
  A.prod(x,r);
  for(j=0; j<n; j++) r(j) = b(j) - r(j);

  if(r.inner(r) < tol)
  {
    it_count_ = 0;
    cpu_time_ = 0.0;
    converged_ = true;
    return;
  }

  PrVec rhat(n);
  for(j=0; j<n; j++) rhat(j) = r(j);
  double rho0 = 1.0, alpha = 1.0, omega0 = 1.0;

  double rho1,omega1,beta;
  PrVec s(n);
  PrVec t(n);

  // All these must be set to zero. Changed PrVec so it actually happens!
  PrVec v0(n);
  PrVec v1(n);
  PrVec p0(n);
  PrVec p1(n);
  double snorm;

  for(int i=1; i<= max_iterations_; i++)
  {
    //s_o << "i = " << i << endl;

    rho1 = rhat.inner(r);
    beta = (rho1 / rho0) * (alpha / omega0);

    //p1 = r + beta * (p0 - omega0 * v0)
    for(j=0; j<n; j++) p1(j) = r(j) + beta * (p0(j) - omega0 * v0(j));

    //v1 = A * p1
    A.prod(p1,v1);

    alpha = rho1 / rhat.inner(v1);

    //s = r - alpha * v1
    for(j=0; j<n; j++) s(j) = r(j) - alpha * v1(j);

    snorm = s.inner(s);
    //s_o << "snorm = " << snorm << endl;

    if(snorm < tol)
    {
      //x = x + alpha * p1
      for(j=0; j<n; j++) x(j) += alpha * p1(j);

      it_count_ = i;
      cpu_time_ = Go::getCurrentTime() - time0;
      converged_ = true;
      return;
    }

    //t = A * s
    A.prod(s,t);

    omega1 = t.inner(s) / t.inner(t);

    //x = x + alpha * p1 + omega1 * s
    for(j=0; j<n; j++) x(j) += alpha * p1(j) + omega1 * s(j);

    //r = s - omega1 * t
    for(j=0; j<n; j++) r(j) = s(j) - omega1 * t(j);

    rho0 = rho1;
    omega0 = omega1;

    //v0 = v1
    for(j=0; j<n; j++) v0(j) = v1(j);
    //p0 = p1
    for(j=0; j<n; j++) p0(j) = p1(j);
  }

  it_count_ = max_iterations_;
  cpu_time_ = Go::getCurrentTime() - time0;
  converged_ = false;

 
}

