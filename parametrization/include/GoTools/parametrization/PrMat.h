/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRMAT_H
#define PRMAT_H

#include <vector>
#ifndef __BORLANDC__
class ostream;
#endif

#include "GoTools/parametrization/PrMatrix.h"

/*<PrMat-syntax: */

/// This class implements a matrix
class PrMat : public PrMatrix
{
protected:
  // elements stored column-wise
  std::vector<double> a_;
  int m_,n_;

public:
  /** @name Derived from base class */
  //@{
  virtual int rows() const {return m_;}
  virtual int colmns() const {return n_;}
  virtual double operator () (int i, int j) const;
  /// Find y = Ax
  virtual void prod(const PrVec& x, PrVec& y) const;
  virtual void read(std::istream& is);

  virtual ~PrMat();
  //@}

  /** @name Other functions */
  //@{
  /// Constructor
  /// \param m number of columns
  /// \param n number of rows
  PrMat(int m, int n, double val = 0.0);
  /// Empty default constructor
  PrMat() {}
  /// access element (i, j)
  double& operator () (int i, int j);
  /// change size of matrix
  void redim(int m, int n);
  //@}
};

/*>PrMat-syntax: */

//-----------------------------------------------------------------------------
inline double& PrMat::operator () (int i, int j)
//-----------------------------------------------------------------------------
{
  return a_[i*n_+j];
}

/*Class:PrMat

Name:              PrMat
Syntax:	           @PrMat-syntax
Keywords:
Description:
Member functions:

Constructors:
Files:
Example:

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Dec. 98
*/

#endif // PRMAT_H
