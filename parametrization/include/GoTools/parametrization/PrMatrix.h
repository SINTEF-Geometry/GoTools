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
