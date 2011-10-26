/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRMATRIX_H
#define PRMATRIX_H


#include "GoTools/parametrization/PrVec.h"

/*<PrMatrix-syntax: */

/// This class implements a matrix
class PrMatrix
{

public:
  // pure virtual functions...
  /// Number of rows
  virtual int rows() const  = 0;
  /// Number of columns
  virtual int colmns() const = 0;
  /// value of matrix element
  virtual double operator () (int i, int j) const = 0;
  /// Multiply matrix with 'x' and return result in 'y'.  (y = Ax)
  virtual void prod(const PrVec& x, PrVec& y) const = 0;
  /// Virtual destructor
  virtual ~PrMatrix();

  /// print contents to stream
  virtual void print(std::ostream& os);
  /// read contents from stream
  virtual void read(std::istream& is) = 0;
};

/*>PrMatrix-syntax: */

/*Class:PrMatrix

Name:              PrMatrix
Syntax:	           @PrMatrix-syntax
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

#endif // PRMATRIX_H
