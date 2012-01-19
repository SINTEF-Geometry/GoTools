//===========================================================================
//                                                                           
// File: SolveBCG.h                                                        
//                                                                           
// Created: Wed Sep 15 09:39:53 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: 
//                                                                           
// Description: Using the biconjugate gradient method to solve sparse and
//              positive definite matrix systems.
//              If matrix is symmetric (but indefinite) we switch to the
//              minimal residual method. Dependent on a good preconditioner.
//                                                                           
//===========================================================================

#ifndef _SOLVEBCG_H
#define _SOLVEBCG_H


#include "GoTools/creators/SolveCG.h"

namespace Go
{

  /// Using the biconjugate gradient method to solve sparse and
  /// positive definite matrix systems.

class SolveBCG : public SolveCG
{

public:

  /// Constructor.
  /// \param conv_type in the range 1 to 4. Convergence estimate. The selection
  /// of preconditioner is based on this number
  SolveBCG(int conv_type, bool symm);

  /// Destructor.
  virtual ~SolveBCG();

  /// Prepare for preconditioning.
  virtual void precond(double relaxfac);

  /// Solve the equation system.
  virtual int solve(double *ex, double *eb, int nn);


private:

  double omega_;        // Relaxation parameter.
  int conv_type_; // We operate with 4 different types of convergence estimates.
  bool symm_;

  void precondRILU(double relaxfac); // @@sbr Not sure RILU is what we will be using...

  // Helper funtion, return largest element in array (i.e. the sup-norm).
  double max_abs_element(double* x, int nn);

};

} // end namespace Go


#endif // _SOLVEBCG_H

