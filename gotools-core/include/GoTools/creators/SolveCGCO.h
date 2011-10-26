#ifndef _SOLVECGCO_H_
#define _SOLVECGCO_H_


#include "GoTools/creators/SolveCG.h"


//   -----------------------------------------------------------------------
//      Interface file for class SolveCGCO
//   -----------------------------------------------------------------------
//
//       Solve the equation system Ax=b where A is a symmetric
//       positive definite matrix using Conjugate Gradient Method.
//       Assuming A describes a constrained optimization problem, i.e.
//           (B^TB C^T)
//       A = (        ),
//           (C     0 )
//       We create a block preconditioner using the CG preconditioner for the
//       upper left block, and a (currently) diagonal preconditioner to the
//       lower right.
//
//   -----------------------------------------------------------------------
//    Written by: Vibeke Skytt                            09.99
//    Based on  : PrCGCO.h written by Mike Floater
//   -----------------------------------------------------------------------

namespace Go
{

class SolveCGCO : public SolveCG
{

public:

  // Constructor.
  // m: size of original system.
  // n: number of constraints.
  SolveCGCO(int m, int n);

  // Prepare for preconditioning.
  virtual void precondRILU(double relaxfac);

  // Destructor.
  virtual ~SolveCGCO();


private:


  int m_;
  int n_;

};

} // end namespace Go

#endif // _SOLVECGCO_H_
