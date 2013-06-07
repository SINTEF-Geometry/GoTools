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

#include "GoTools/parametrization/PrPrewavelet_F.h"
#include "GoTools/parametrization/PrWaveletUtil.h"

#include "GoTools/parametrization/PrCG.h"
#include "GoTools/parametrization/PrMatSparse.h"
#include "GoTools/parametrization/PrVec.h"


//-----------------------------------------------------------------------------
PrPrewavelet_F::~PrPrewavelet_F()
//-----------------------------------------------------------------------------
{
}


// PRIVATE METHODS

//-----------------------------------------------------------------------------
void
PrPrewavelet_F::makeSparseMatrix(int jlev, int dim)
//-----------------------------------------------------------------------------
// The nodes in the triangulation belonging to V^j
// form a piecewise linear function f^j in the space S^j.
// There is a unique decomposition f^j = f^{j-1} + g^{j-1}
// where f^j \in S^{j-1} and g^j \in W^{j-1}.
// This routine constructs the matrix Q2 - P2 Q1 which
// can be used to decompose f^j into f^{j-1} and g^{j-1}.
//
// Specifically we have to solve
//
//    c + Q1 d = c_old              (i)
//    P2 c + Q2 d = c_new           (ii)
//
// where [c_old c_new]^T are the given coefficients of f^j at level j
// and c and d are the coefficients of f^{j-1} and g^{j-1} respectively.
// Subtracting P2 times (i) from (ii) we replace (ii) by
//
//   (Q2 - P2 Q1) d = c_new - P2 c_old.    (iii)
//
// This routine finds the matrix A := Q2 - P2 Q1 and the right hand side
// b := c_new - P2 c_old.
// Once d has been found, it can be substituted into (i) so that
//
//  c = c_old - Q1 d.
//
{
  if(dim != 1 && dim != 3) return;
  if(jlev < 1 || jlev > t_->getFinestLevel()) return;

// Scale all prewavelets by 2 to the power (jlev-1).
  int powerOfTwo = 1 << (jlev - 1);
  int n = t_->getNumNodes(jlev);
  int nprev = t_->getNumNodes(jlev-1);
  int nnew = n - nprev;
//    cout << "n = " << n << endl;
//    cout << "nprev = " << nprev << endl;
//    cout << "nnew = " << nnew << endl;

  int i,j,k;
  // initialize right hand side to c_new.
  if(dim == 3)
  {
    b1_.redim(nnew);
    for(i=nprev; i<n; i++) b1_[i-nprev] = t_->getX(i);
    b2_.redim(nnew);
    for(i=nprev; i<n; i++) b2_[i-nprev] = t_->getY(i);
  }
  b3_.redim(nnew);
  for(i=nprev; i<n; i++) b3_[i-nprev] = t_->getZ(i);

  // The matrix A will have the structure A = Q2 - P2 Q1
  // where P1 is nprev * nprev
  // where P2 is nnew * nprev
  // where Q1 is nprev * nnew
  // where Q2 is nnew * nnew
  // where also nnew = n - nprev.
  // Nodes 1,...,nprev as the "old" nodes, i.e. nodes in V^{j-1}
  // and nodes nprev+1,...,n as the new nodes, i.e. nodes in V^j \ V^{j-1}.

  // Store the cumulative number of nonzeros on each row in the
  // array sparseDS.irow.
  // The number of nonzeros of the row of A corresponding
  // to a vertex u is d(v1) + d(v2) - 1 where d(v) is the
  // degree of the vertex v and v1 and v2 are the two coarse
  // neighbours of u.
  // If the number of nonzeros in the rows are a, b, c, ..., k then
  // sparseDS.irow will contain 0, a, a+b, ..., a+b+c...+k.

  vector<int> offset(nnew+1);
  vector<int> nonZeroesPerLine(nnew, -1);

//  for(i=1; i<= nnew; i++) offset[i] = -1;

  int t,jj,kk;
  for(i=0; i< nprev; i++)
  {
    t_->getNeighbours(i,jlev,neighbours_);
    for(j=0; j< int(neighbours_.size()); j++)
    {
      jj = neighbours_[j] - nprev;
//      offset[jj+1] += neighbours_.size();
      nonZeroesPerLine[jj] += (int)neighbours_.size();
    }
  }
  // Now nonZeroesPerLine contains ?,a,b,c,...,k,
//  // Now offset contains ?,a,b,c,...,k,

  // Accumulate
  offset[0] = 0;
//  for(i=1; i<= nnew; i++) offset[i] += offset[i-1];
  for(i=1; i<= nnew; i++) offset[i] = offset[i-1] + nonZeroesPerLine[i-1];

  // Total number of nozeros in A is:
  int numNonZeros = offset[nnew];
  A_.redim(nnew,nnew,numNonZeros);

  for(i=0; i<= nnew; i++) A_.irow(i) = offset[i];

  // Initialize A to zero because we'll ADD terms later.

  for(i=0; i< numNonZeros; i++) A_(i) = 0.0;

  double semiprewavelet;
  for(i=0; i< nprev; i++)
  {
    t_->getNeighbours(i,jlev,neighbours_);
    if(t_->isBoundary(i)) t = (int)neighbours_.size() - 1;
    else                  t = (int)neighbours_.size();

    for(j=0; j< int(neighbours_.size()); j++)
    {
      jj = neighbours_[j] - nprev;
      if(dim == 3)
      {
        b1_[jj] -= 0.5 * t_->getX(i);
        b2_[jj] -= 0.5 * t_->getY(i);
      }
      b3_[jj] -= 0.5 * t_->getZ(i);

      for(k=0; k< int(neighbours_.size()); k++)
      {
        kk = neighbours_[k] - nprev;
        // The diagonal will be the sum of two terms.
        semiprewavelet = powerOfTwo * (6.0 / (double)(7 * t)
				       + theta(j,k,i,(int)neighbours_.size(),t_->isBoundary(i)));

//	cout << "spw (" << i << "," << j << "," << k << ") = " << semiprewavelet << endl;

        if(k == j)
        {
          A_(A_.irow(jj)) += semiprewavelet;
          A_.jcol(A_.irow(jj)) = kk;
        }
        else
        {
          offset[jj]++;
          A_(offset[jj]) += semiprewavelet;
          A_.jcol(offset[jj]) = kk;
//          offset[jj]++;
        }
      }
    }
  }

//  A_.print(cerr);
//  b1_.print(cerr);
//  b2_.print(cerr);
//  b3_.print(cerr);
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
void
PrPrewavelet_F::decompose(int jlev, int dim)
//-----------------------------------------------------------------------------
// The nodes in the triangulation belonging to V^j
// form a piecewise linear function f^j in the space S^j.
// This routines decomposes f^j
// into f^{j-1} + g^{j-1} where
// f^{j-1} \in S^{j-1} and g^{j-1} \in W^{j-1}.
{
  if(dim != 1 && dim != 3) return;
  if(jlev < 1 || jlev > t_->getFinestLevel()) return;

// Make the matrix A = Q2 - P2 Q1 and right hand side c_new - P2 c_old.
  makeSparseMatrix(jlev,dim);

// Scale all prewavelets by 2 to the power (jlev-1).
  int powerOfTwo = 1 << (jlev - 1);
  int n = t_->getNumNodes(jlev);
  int nprev = t_->getNumNodes(jlev-1);
  int nnew = n - nprev;
  if(dim == 3)
  {
    d1_.redim(nnew);
    c1_.redim(nprev);
    d2_.redim(nnew);
    c2_.redim(nprev);
  }
  d3_.redim(nnew);
  c3_.redim(nprev);
  int i,j;

  // Choose start vector. After testing, zero was found
  // to be just as good as any (of those commented out).

  if(dim == 3)
  {
    for(i=0; i<nnew; i++) d1_[i] = 0.0;
    for(i=0; i<nnew; i++) d2_[i] = 0.0;
  }
  for(i=0; i<nnew; i++) d3_[i] = 0.0;

  //for(i=1; i<=nnew; i++) (*dp_)(i) = (*bp_)(i);
  //for(i=1; i<=nnew; i++) (*dp_)(i) = 0.5 * (*bp_)(i);
  //for(i=1; i<=nnew; i++) (*dp_)(i) = (*bp_)(i) / (*Ap_)(Ap_->irow(i));

  PrCG solver;
  solver.setMaxIterations(nnew);
  solver.setTolerance(tolerance_);

  if(dim == 3)
  {
    solver.solve(A_,d1_,b1_);
    solver.solve(A_,d2_,b2_);
  }
  solver.solve(A_,d3_,b3_);
   
  double sum;
  int t;
  for(i=0; i< nprev; i++)
  {
    t_->getNeighbours(i,jlev,neighbours_);
    if(t_->isBoundary(i)) t = (int)neighbours_.size() - 1;
    else                  t = (int)neighbours_.size();

    if(dim == 3)
    {
      sum = 0.0;
      for(j=0; j< int(neighbours_.size()); j++)
              sum += d1_[neighbours_[j]-nprev];
      c1_[i] = t_->getX(i) + powerOfTwo * (1.5 / (double)t) * sum;
      sum = 0.0;
      for(j=0; j< int(neighbours_.size()); j++)
              sum += d2_[neighbours_[j]-nprev];
      c2_[i] = t_->getY(i) + powerOfTwo * (1.5 / (double)t) * sum;
    }
    sum = 0.0;
    for(j=0; j< int(neighbours_.size()); j++)
            sum += d3_[neighbours_[j]-nprev];
    c3_[i] = t_->getZ(i) + powerOfTwo * (1.5 / (double)t) * sum;
  }

  // Copy the result to the triangulation.
  // This could be changed if we don't want to overwrite.

  if(dim == 3)
  {
    for(i=0; i<nprev; i++) t_->setX(i,c1_[i]);
    for(i=nprev; i<n; i++) t_->setX(i,d1_[i-nprev]);
    for(i=0; i<nprev; i++) t_->setY(i,c2_[i]);
    for(i=nprev; i<n; i++) t_->setY(i,d2_[i-nprev]);
  }
  for(i=0; i<nprev; i++) t_->setZ(i,c3_[i]);
  for(i=nprev; i<n; i++) t_->setZ(i,d3_[i-nprev]);
}

//-----------------------------------------------------------------------------
void
PrPrewavelet_F::compose(int jlev, int dim)
//-----------------------------------------------------------------------------
// The nodes in the triangulation belonging to V^j
// form a piecewise linear function f^j in the space S^j.
// This routines composes f^j
// from f^{j-1} + g^{j-1} where
// f^{j-1} \in S^{j-1} and g^{j-1} \in W^{j-1}.
{
  if(dim != 1 && dim != 3) return;
  if(jlev < 1 || jlev > t_->getFinestLevel()) return;

// Scale all prewavelets by 2 to the power (jlev-1).
  int powerOfTwo = 1 << (jlev - 1);
  int n = t_->getNumNodes(jlev);
  int nprev = t_->getNumNodes(jlev-1);

  if(dim == 3)
  {
    b1_.redim(n);
    b2_.redim(n);
  }
  b3_.redim(n);

// Apply the matrix A = [P Q] to generate b = Ax, where
// x is the combined coefficient vector of f^{j-1} and g^{j-1}.

  if(dim == 3)
  {
    for(int i=0; i< n; i++) b1_[i] = 0.0;
    for(int i=0; i< n; i++) b2_[i] = 0.0;
  }
  for(int i=0; i< n; i++) b3_[i] = 0.0;
  
  double sum, sum2;
  int t;
  for(int i=0; i< nprev; i++)
  {
    t_->getNeighbours(i,jlev,neighbours_);
    if(t_->isBoundary(i)) t = (int)neighbours_.size() - 1;
    else                  t = (int)neighbours_.size();

    if(dim == 3)
    {
      sum = 0.0;
      for(size_t j=0; j< neighbours_.size(); j++)
      {
        sum += t_->getX(neighbours_[j]);
        sum2 = 0.0;
        for(size_t k=0; k< neighbours_.size(); k++)
        {
          sum2 += (3.0 / (double)(28 * t)
		   + theta((int)j,(int)k,(int)i,(int)neighbours_.size(),t_->isBoundary(i)) )
                        * t_->getX(neighbours_[k]);
        }
        b1_[neighbours_[j]] += 0.5 * t_->getX(i) + powerOfTwo * sum2;
      }
      b1_[i] = t_->getX(i) - powerOfTwo * (1.5 / (double)t) * sum;

      sum = 0.0;
      for(size_t j=0; j< neighbours_.size(); j++)
      {
        sum += t_->getY(neighbours_[j]);
        sum2 = 0.0;
        for(size_t k=0; k< neighbours_.size(); k++)
        {
          sum2 += (3.0 / (double)(28 * t)
		   + theta((int)j,(int)k,(int)i,(int)neighbours_.size(),t_->isBoundary(i)) )
                        * t_->getY(neighbours_[k]);
        }
        b2_[neighbours_[j]] += 0.5 * t_->getY(i) + powerOfTwo * sum2;
      }
      b2_[i] = t_->getY(i) - powerOfTwo * (1.5 / (double)t) * sum;
    }
      sum = 0.0;
      for(size_t j=0; j< neighbours_.size(); j++)
      {
        sum += t_->getZ(neighbours_[j]);
        sum2 = 0.0;
        for(size_t k=0; k< neighbours_.size(); k++)
        {
          sum2 += (3.0 / (double)(28 * t)
		   + theta((int)j,(int)k,(int)i,(int)neighbours_.size(),t_->isBoundary(i)) )
                        * t_->getZ(neighbours_[k]);
        }
        b3_[neighbours_[j]] += 0.5 * t_->getZ(i) + powerOfTwo * sum2;
      }
      b3_[i] = t_->getZ(i) - powerOfTwo * (1.5 / (double)t) * sum;
  }

  // Copy the result to the triangulation.
  // This could be changed if we want to retain the wavelet
  // decomposition.

  if(dim == 3)
  {
    for(int i=0; i<n; i++) t_->setX(i,b1_[i]);
    for(int i=0; i<n; i++) t_->setY(i,b2_[i]);
  }
  for(int i=0; i<n; i++) t_->setZ(i,b3_[i]);
}

