/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRFILTERBANK_H
#define PRFILTERBANK_H

#include "GoTools/parametrization/PrNestedTriangulation.h"
#include "GoTools/utils/config.h"
#include <memory>
#include <iostream>
using std::istream;
using std::ostream;
using std::cout;
using std::endl;

/*<PrFilterbank-syntax: */

/** PrFilterbank - implements an algorithm for decomposing
 * a piecewise linear function \f$f^k\f$ over a nested
 * triangulation \f$T^k\f$ into the sum
 * \f$f^0 + g^0 + g^1 + ... + g^{k-1}\f$
 * where \f$f^0 \in S^0\f$ and \f$g^0 \in W^0, ... g^{k-1} \in W^{k-1}\f$.
 * Here the set of nodes \f$V^k\f$ of \f$T^k\f$ satisfy
 * \f$V^0 \subset V^1 \subset ... \subset V^k\f$
 * and the corresponding spaces of piecewise linear 
 * \f$S^0 \subset S^1 \subset ... \subset S^k\f$.
 * The space \f$W^{j-1}\f$ is a complement space of \f$S^{j-1} \in S^j\f$.
 * Examples are prewavelet and Faber decomposition.
 */
class PrFilterbank
{
protected:

  shared_ptr<PrNestedTriangulation> t_;
  vector<int>            neighbours_;

public:
  // pure virtual functions

  /** The nodes in the triangulation belonging to \f$V^j\f$
   * form a piecewise linear function \f$f^j\f$ in the space \f$S^j\f$
   * (where j = level).
   * This routines decomposes \f$f^j\f$
   * into \f$f^{j-1} + g^{j-1}\f$ where
   * \f$f^{j-1} \in S^{j-1}\f$ and \f$g^{j-1} \in W^{j-1}\f$.
   * Use dim = 1 for explicit triangulations and
   * dim = 3 for non-explicit. Default is dim = 1.
  */
  virtual void   decompose(int jlev, int dim = 1) = 0;

  /** The nodes in the triangulation belonging to \f$V^j\f$
   * form a piecewise linear function \f$f^j\f$ in the space \f$S^j\f$.
   * This routines composes \f$f^j\f$
   * from \f$f^{j-1} + g^{j-1}\f$ where
   * \f$f^{j-1} \in S^{j-1}\f$ and \f$g^{j-1} \in W^{j-1}\f$.
   * Use dim = 1 for explicit triangulations 
   * dim = 3 for non-explicit. Default is dim = 1.
   */
  virtual void   compose(int jlev, int dim = 1) = 0;

  // other functions

  /// Destructor
  virtual ~PrFilterbank();
  /// Set the graph.
  void       attach(shared_ptr<PrNestedTriangulation> t);
  /** The triangulation is a piecewise linear function
   * \f$f^k\f$ in the space \f$S^k\f$. This routines decomposes \f$f^k\f$
   * into \f$f^0 + g^0 + g^1 + ... + g^{k-1}\f$ where
   * \f$f^0 \in S^0, g^0 \in W^0, ... g^{k-1} \in W^{k-1}\f$.
   * Use dim = 1 for explicit triangulations and
   * dim = 3 for non-explicit. Default is dim = 1.
   */
  void       decomposeAll(int dim = 1);
  void       decomposeFrom(int jlev, int dim = 1);
               // 1 <= jlev <= finest level

  /** The coarsest triangulation is a piecewise linear function
   * \f$f^0\f$ in the space \f$S^0\f$. This routines composes \f$f^k\f$
   * from \f$f^0 + g^0 + g^1 + ... + g^{k-1}\f$ where
   * \f$f^k \in S^k, g^0 \in W^0, ... g^{k-1} \in W^{k-1}\f$.
   * Use dim = 1 for explicit triangulations and
   * dim = 3 for non-explicit. Default is dim = 1.
   */
  void       composeAll(int dim = 1);
  void       composeUpTo(int jlev, int dim = 1);
               // 1 <= jlev <= finest level
};

/*>PrFilterbank-syntax: */

/*Class:PrFilterbank

Name:              PrFilterbank
Syntax:	           @PrFilterbank-syntax
Keywords:
Description:       This class implements an algorithm for decomposing
                   a piecewise linear function f^k over a nested
                   triangulation T^k into the sum
                   
                       f^0 + g^0 + g^1 + ... + g^{k-1}

                   where f^0 \in S^0 and g^0 \in W^0, ... g^{k-1} \in W^{k-1}.
                   Here the set of nodes V^k of T^k satisfy

                       V^0 \subset V^1 \subset ... \subset V^k

                   and the corresponding spaces of piecewise linear 
                   functions are

                       S^0 \subset S^1 \subset ... \subset S^k

                   The space W^{j-1} is a complement space of S^{j-1} in S^j.
                   Examples are prewavelet and Faber decomposition.

Member functions:
                   "attach(PrNestedTriangulation& t)" --\\
                   Set the graph.

                   "decompose(int level, int dim)" --\\
                   The nodes in the triangulation belonging to V^j
                   form a piecewise linear function f^j in the space S^j
                   (where j = level).
                   This routines decomposes f^j
                   into f^{j-1} + g^{j-1} where
                   f^{j-1} \in S^{j-1} and g^{j-1} \in W^{j-1}.
                   Use dim = 1 for explicit triangulations and
                   dim = 3 for non-explicit. Default is dim = 1.
                   
                   "decomposeAll(int dim)" --\\
                   The triangulation is a piecewise linear function
                   f^k in the space S^k. This routines decomposes f^k
                   into f^0 + g^0 + g^1 + ... + g^{k-1} where
                   f^0 \in S^0, g^0 \in W^0, ... g^{k-1} \in W^{k-1}.
                   Use dim = 1 for explicit triangulations and
                   dim = 3 for non-explicit. Default is dim = 1.

                   "compose(int level, int dim)" --\\
                   The nodes in the triangulation belonging to V^j
                   form a piecewise linear function f^j in the space S^j.
                   This routines composes f^j
                   from f^{j-1} + g^{j-1} where
                   f^{j-1} \in S^{j-1} and g^{j-1} \in W^{j-1}.
                   Use dim = 1 for explicit triangulations and
                   dim = 3 for non-explicit. Default is dim = 1.
                   
                   "composeAll(int dim)" --\\
                   The coarsest triangulation is a piecewise linear function
                   f^0 in the space S^0. This routines composes f^k
                   from f^0 + g^0 + g^1 + ... + g^{k-1} where
                   f^k \in S^k, g^0 \in W^0, ... g^{k-1} \in W^{k-1}.
                   Use dim = 1 for explicit triangulations and
                   dim = 3 for non-explicit. Default is dim = 1.

Constructors:
Files:
Example:

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Feb. 98
Revised:           Dec. 98 (changed to abstract base class).
Revised:           Dec. 2000 (remove BasicTools dependencies).
*/

#endif // PRFILTERBANK_H
