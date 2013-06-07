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

#ifndef PRMATSPARSE_H
#define PRMATSPARSE_H

#include "GoTools/parametrization/PrMatrix.h"
#include "GoTools/utils/errormacros.h"
#include <vector>

/*<PrMatSparse-syntax: */

class PrMatSparse : public PrMatrix
{
private:
  int m_, n_; // matrix dimensions
  int p_;     // number of non-zeros
  std::vector<int> irow_; // the indexes in a_ and jcol_ of the
                     // first non-zeros of the m rows
  std::vector<int> jcol_; // the p indexes j of the non-zero elements
  std::vector<double> a_;   // the p non-zero elements
  

public:

  /** @name Derived from base class */
  //@{
  virtual int rows() const;
  virtual int colmns() const;
  virtual double operator () (int i, int j) const;
  /// Find y = Ax
  virtual void prod(const PrVec& x, PrVec& y) const;
  virtual void print(std::ostream& os);
  virtual void read(std::istream& is);
  /// virtual destructor
  virtual ~PrMatSparse();
  //@}

  /** @name Other functions */
  //@{
  /// find C = A*B.
  /// Multiplies two sparse matrices, "A=(this)" times "B" and
  /// stores the result in "C"
  void matProd(PrMatSparse& B, PrMatSparse& C) const;
  /// multiply all elements with 'd'
  void scalProd(double d);
  /// print the matrix to stream in 'full' (not sparse) format
  void printFull(std::ostream& os);

  /// constructor for a sparse matrix of zero size.
  PrMatSparse() : m_(0), n_(0), p_(0) {}
  /// constructor for an (m x n) matrix with 'num_nonzero' nonzero elements
  PrMatSparse(int m, int n, int num_nonzero);
  /// Set to sparse matrix given by irow, jcol and data.
  PrMatSparse(int m, int n, int num_nonzero,
	      const int* irow, const int* jcol, const double* data);
  /// resize matrix to (m x n), with 'num_nonzero' nonzero elements
  void redim(int m, int n, int num_nonzero);
  /// Set equal to another matrix, sparsify if entries <= tol.
  void setToMatrix(const PrMatrix& m, double tol = 0.0);

    /// < 0 <= k <= m_ (last element is past-the-end index of a_)
    int& irow(int k);        
    /// < 0 <= k <= m_ (last element is past-the-end index of a_)
    const int& irow(int k) const;  
    /// < 0 <= k <= p_-1
    int& jcol(int k);        
    /// < 0 <= k <= p_-1
    const int& jcol(int k) const; 
    /// < 0 <= k <= p_-1
    double& operator () (int k);
    /// < 0 <= k <= p_-1
    const double& operator () (int k) const;
  //@}
};


/*>PrMatSparse-syntax: */

//-----------------------------------------------------------------------------
inline void PrMatSparse::scalProd(double d)
//-----------------------------------------------------------------------------
{
    for (std::vector<double>::iterator it = a_.begin(); it != a_.end(); ++it) {
	(*it) *= d;
    }
}


//-----------------------------------------------------------------------------
inline const int& PrMatSparse::irow(int k) const
//-----------------------------------------------------------------------------
{
  return irow_[k];
}

//-----------------------------------------------------------------------------
inline int& PrMatSparse::irow(int k)
//-----------------------------------------------------------------------------
{
  return irow_[k];
}

//-----------------------------------------------------------------------------
inline const int& PrMatSparse::jcol(int k) const
//-----------------------------------------------------------------------------
{
  return jcol_[k];
}

//-----------------------------------------------------------------------------
inline int& PrMatSparse::jcol(int k)
//-----------------------------------------------------------------------------
{
  return jcol_[k];
}

//-----------------------------------------------------------------------------
inline const double& PrMatSparse::operator() (int k) const
//-----------------------------------------------------------------------------
{
  return a_[k];
}

//-----------------------------------------------------------------------------
inline double& PrMatSparse::operator() (int k)
//-----------------------------------------------------------------------------
{
  return a_[k];
}

/*Class:PrMatSparse

Name:              PrMatSparse
Syntax:	           @PrMatSparse-syntax
Keywords:
Description:       This class implements a matrix
Member functions:

Constructors:
Files:
Example:

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Dec. 98
*/

#endif // PRMATSPARSE_H





