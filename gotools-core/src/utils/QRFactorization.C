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

#include "GoTools/utils/QRFactorization.h"
#include "GoTools/utils/errormacros.h"
#include <cmath>
#include <vector>

using namespace Go;
using std::vector;

//==============================================================================
void QRFactorization::QRDecomp(vector<double>& A, int n, int m,
			       vector<double>& Q, vector<double>& R)
//==============================================================================
{
  if ((int)A.size() != n*m)
    THROW("QRFactorization::QRDecomp: Expecting matrix size n*m");
  
  double eps = 1.0e-12;
  int kt = std::min(n,m);

  R = A;
  Q.resize(m*m, 0.0);
  for (int kr=0; kr<m; ++kr)
    Q[kr*m+kr] = 1.0;

  for (int kr=0; kr<kt; ++kr)
    {
      vector<double> v0(m-kr, 0.0);
      for (int ki=kr; ki<m; ++ki)
	v0[ki-kr] = R[ki*n+kr];

      double normv0 = 0.0;
      for (size_t i=0; i<v0.size(); ++i)
	normv0 += v0[i]*v0[i];
      normv0 = sqrt(normv0);
      //double sgn = (v0[0] >= 0.0) ? 1.0 : -1.0;
      double sgn = (fabs(v0[0]-normv0) < 0.1*normv0) ? 1.0 : -1.0;

      vector<double> v(v0.begin(), v0.end());
      v[0] += sgn*normv0;
      double normv = 0.0;
      for (size_t i=0; i<v.size(); ++i)
	normv += v[i]*v[i];
      normv = sqrt(normv);
      if (normv > eps)
	{
	  for (size_t i=0; i<v.size(); ++i)
	    v[i] /= normv;
	}

      for (int kj=kr; kj<n; ++kj)
	{
	  double dot = 0.0;
	  for (int ki=kr; ki<m; ++ki)
	    dot += v[ki-kr]*R[ki*n+kj];
	  for (int ki=kr; ki<m; ++ki)
	    R[ki*n+kj] -= 2.0*v[ki-kr]*dot;
	}

      for (int kj=0; kj<m; ++kj)
	{
	  double dot = 0.0;
	  for (int ki=kr; ki<m; ++ki)
	    dot += v[ki-kr]*Q[ki*m+kj];
	  for (int ki=kr; ki<m; ++ki)
	    Q[ki*m+kj] -= 2.0*v[ki-kr]*dot;
	}
    }
}

//==============================================================================
void QRFactorization::QRSolve(vector<double>& Q, vector<double>& R, int n, int m,
			      vector<double>& b, int dim, vector<double>& x)
//==============================================================================
{
  if ((int)Q.size() != m*m)
    THROW("QRFactorization::QRSolve: Expecting matrix size m*m");
  if ((int)b.size() != dim*m)
    THROW("QRFactorization::QRSolve: Expecting left hand size dim*m");
    
  x.resize(dim*n, 0.0);
  for (int kb=0; kb<dim; ++kb)
    {
      vector<double> b2(n, 0.0);
      for (int kj=0; kj<n; ++kj)
	for (int ki=0; ki<m; ++ki)
	  b2[kj] += Q[kj*m+ki]*b[kb*m+ki];

      for (int kr=n-1; kr>=0; --kr)
	{
	  x[kb*n+kr] = b2[kr];
	  double val = 0.0;
	  for (int ki=kr+1; ki<n; ++ki)
	    val += R[kr*n+ki]*x[kb*n+ki];
	  x[kb*n+kr] -= val;
	  if (fabs(R[kr*n+kr]) < 1.0e-18)
	    {
	      MESSAGE("QRFactorization::QRSolve: Division with zero");
	      x[kb*n+kr] = 0.0;
	    }
	  else
	    x[kb*n+kr] /= R[kr*n+kr];
	}
    }
}

