/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRPREWAVELET_F_H
#define PRPREWAVELET_F_H

#include "GoTools/parametrization/PrFilterbank.h"
#include "GoTools/parametrization/PrMatSparse.h"
#include "GoTools/parametrization/PrVec.h"


/*<PrPrewavelet_F-syntax: */

/** PrPrewavelet _  This class implements algorithms for decomposing
 * and decomposing piecewise linear functions over nested
 * triangulations. The complemenet spaces are taken
 * here to be wavelet spaces with the prewavelet basis
 * constructed in
 * 
 * M. S. Floater and E. G. Quak,
 * Piecewise Linear Prewavelets on Arbitrary Triangulations,
 * Numer. Math. \b 82 (1999), 221--252.
 */
class PrPrewavelet_F : public PrFilterbank
{
private:

  PrMatSparse A_;
  PrVec b1_,b2_,b3_;
  PrVec c1_,c2_,c3_;
  PrVec d1_,d2_,d3_;
  double tolerance_;
  vector<int> neighbours_;

// PRIVATE METHODS

  void makeSparseMatrix(int jlev, int dim = 1);

public:
  /// Default constructor
  PrPrewavelet_F() {tolerance_ = 1.0e-6; }
  /// Empty destructor
  virtual ~PrPrewavelet_F();

  /// Set the numerical tolerance used for internal computations
  /// with the conjugated gradient algorithm.
  void    setCGTolerance(double tolerance = 1.0e-6)
                  {tolerance_ = tolerance;}


  /** The nodes in the triangulation belonging to \f$V^j\f$
   * form a piecewise linear function \f$f^j\f$ in the space \f$S^j\f$
   * (where j = level).
   * This routines decomposes \f$f^j\f$
   * into \f$f^{j-1} + g^{j-1}\f$ where
   * \f$f^{j-1} \in S^{j-1}\f$ and \f$g^{j-1} \in W^{j-1}\f$.
   * Use dim = 1 for explicit triangulations and
   * dim = 3 for non-explicit. Default is dim = 1.
   */
  virtual void decompose(int jlev, int dim = 1);
                  // 1 <= jlev <= finest level
  /** The nodes in the triangulation belonging to \f$V^j\f$
   * form a piecewise linear function \f$f^j\f$ in the space \f$S^j\f$.
   * This routines composes \f$f^j\f$
   * from \f$f^{j-1}\f$ + \f$g^{j-1}\f$ where
   * \f$f^{j-1} \in S^{j-1}\f$ and \f$g^{j-1} \in W^{j-1}\f$.
   * Use dim = 1 for explicit triangulations and
   * dim = 3 for non-explicit. Default is dim = 1.
   */
  virtual void compose(int jlev, int dim = 1);
                  // 1 <= jlev <= finest level

};

/*>PrPrewavelet_F-syntax: */

/*Class:PrPrewavelet_F

Name:              PrPrewavelet_F
Syntax:	           @PrPrewavelet_F-syntax
Keywords:
Description:       This class implements algorithms for decomposing
                   and decomposing piecewise linear functions over nested
                   triangulations. The complemenet spaces are taken
                   here to be wavelet spaces with the prewavelet basis
                   constructed in
       
                   M. S. Floater and E. G. Quak,
                   Piecewise Linear Prewavelets on Arbitrary Triangulations,
                   Numer. Math. {\bf 82} (1999), 221--252.

Member functions:
                   "decompose(int level, int dim)" --\\
                   The nodes in the triangulation belonging to V^j
                   form a piecewise linear function f^j in the space S^j
                   (where j = level).
                   This routines decomposes f^j
                   into f^{j-1} + g^{j-1} where
                   f^{j-1} \in S^{j-1} and g^{j-1} \in W^{j-1}.
                   Use dim = 1 for explicit triangulations and
                   dim = 3 for non-explicit. Default is dim = 1.
		   We use the Conjugated Gradient method with an l2
		   norm absolute residual stop criterion w.r.t. the tolernace.
                   
                   "compose(int level, int dim)" --\\
                   The nodes in the triangulation belonging to V^j
                   form a piecewise linear function f^j in the space S^j.
                   This routines composes f^j
                   from f^{j-1} + g^{j-1} where
                   f^{j-1} \in S^{j-1} and g^{j-1} \in W^{j-1}.
                   Use dim = 1 for explicit triangulations and
                   dim = 3 for non-explicit. Default is dim = 1.
                   
Constructors:
Files:
Example:

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Feb. 98
*/

#endif // PRPREWAVELET_F_H
