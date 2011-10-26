/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRBICGSTAB_H
#define PRBICGSTAB_H

#include "GoTools/parametrization/PrMatrix.h"
#include "GoTools/parametrization/PrVec.h"

/*<PrBiCGStab-syntax: */

/** PrBiCGStab - Implements the BiCGStab method for solving sparse
 * linear systems.
 */
class PrBiCGStab
{
protected:

  double    tolerance_;
  int     max_iterations_;

  int it_count_;
  double cpu_time_;
  bool converged_;

public:
  /// Constructor
  PrBiCGStab();
  /// Destructor
  ~PrBiCGStab() {}

  /// Set the tolerance for the residual.
  void setTolerance(double tolerance = 1.0e-6) {tolerance_ = tolerance;}

  /// Set the maximum number of iterations.
  void setMaxIterations(int max_iterations)
           {max_iterations_ = max_iterations;}
  /// Solve the linear system, replacing the start vector with the solution.
  void solve(const PrMatrix& A, PrVec& x, const PrVec& b);

  /// Get the number of iterations spent for the last call of 'solve()'.
  int getItCount() {return it_count_; }

  /// Get the CPU time spent for the last call of 'solve()'.
  double getCPUTime() {return cpu_time_; }
 
  /// Check if the last call of 'solve()' managed to converge to a solution.
  bool converged() {return converged_; }
};

/*>PrBiCGStab-syntax: */

/*Class:PrBiCGStab

Name:              PrBiCGStab
Syntax:	           @PrBiCGStab-syntax
Keywords:
Description:       This class implements the BiCGStab method
                   for solving sparse linear systems.
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
Date:              Nov. 98
*/

#endif // PRBICGSTAB_H
