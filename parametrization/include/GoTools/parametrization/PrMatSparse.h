/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

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





