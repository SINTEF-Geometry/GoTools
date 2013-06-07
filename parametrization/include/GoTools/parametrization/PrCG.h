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

#ifndef PRCG_H
#define PRCG_H

#include "GoTools/parametrization/PrMatrix.h"
#include "GoTools/parametrization/PrVec.h"

/*<PrCG-syntax: */

/** PrCG - This class implements the (unpreconditioned) CG method
 * for solving sparse symmetric positive definite linear systems.
 */
class PrCG
{
protected:

  double tolerance_;
  int max_iterations_;

  int it_count_;
  double cpu_time_;
  bool converged_;

public:
  /// Constructor
  PrCG();
  /// Destructor
  ~PrCG() {}

  /// Set the tolerance for the residual.
  void setTolerance(double tolerance = 1.0e-6) {tolerance_ = tolerance;}
  /// Set the maximum number of iterations.
  void setMaxIterations(int max_iterations)
           {max_iterations_ = max_iterations;}
  ///Solve the linear system, replacing the start vector with the solution.
  void solve(const PrMatrix& A, PrVec& x, const PrVec& b);

  /// Get the number of iterations spent for the last call of 'solve()'.
  int getItCount() {return it_count_; }

  /// Get the CPU time spent for the last call of 'solve()'.
  double getCPUTime() {return cpu_time_; }

  /// Check if the last call of 'solve()' managed to converge to a solution.
  bool converged() {return converged_; }
};

/*>PrCG-syntax: */

/*Class:PrCG

Name:              PrCG
Syntax:	           @PrCG-syntax
Keywords:
Description:       This class implements the (unpreconditioned) CG method
                   for solving sparse symmetric positive definite
                   linear systems.
Member functions:
                   "setTolerance()" --\\
                   Set the tolerance for the residual.

                   "setMaxIterations()" --\\
                   Set the maximum number of iterations.

                   "solve(const PrMatrix& A, PrVec& x, const PrVec& b)" --\\
                   Solve the linear system, replacing the start vector
                   with the solution.

Constructors:
Files:
Example:

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Dec. 98
*/

#endif // PRCG_H
