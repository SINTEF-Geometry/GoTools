/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRNESTEDTRIANGULATION_H
#define PRNESTEDTRIANGULATION_H

#include <vector>
using std::vector;

/*<PrNestedTriangulation-syntax: */

/** This class implements a set of nested triangulations.
 * The methods in this class are used by PrFilterBank.
 */
class PrNestedTriangulation
{
public:
  /** @name Pure virtual functions... */
  //@{

  /// Return the finest level in the hieararchy, i.e. k if the vertices 
  /// are arranged \f$V^0 \subset V^1 \subset ... \subset V^k\f$.
  virtual int      getFinestLevel() = 0;

  /// return the number of nodes in \f$V^j, j=0,1,...,k\f$.
  virtual int      getNumNodes(int jlev) = 0;

  /// return the x coeff. of the i-th node in the triangulation.
  virtual double     getX(int i) = 0;
  /// return the y coeff. of the i-th node in the triangulation.
  virtual double     getY(int i) = 0;
  /// return the z coeff. of the i-th node in the triangulation.
  virtual double     getZ(int i) = 0;

  /// set the x coeff. of the i-th node in the triangulation.
  virtual void     setX(int i, const double& x) = 0;
  /// set the y coeff. of the i-th node in the triangulation.
  virtual void     setY(int i, const double& y) = 0;
  /// set the z coeff. of the i-th node in the triangulation.
  virtual void     setZ(int i, const double& z) = 0;

  /** Return the indices of the neighbours of the i-th node
   * with respect to \f$V^j\f$. This routine assumes that
   * node i also belongs to \f$V^j\f$. The order of the neighbours
   * returned will be
   * 1. any anticlockwise order if i is an interior node
   * 2. the unique anticlockwise order if i is a boundary node.
   */
  virtual void     getNeighbours(int i, int jlev, vector<int>& neighbours) = 0;
  /// is the i-th node a boundary node or interior node?
  virtual bool   isBoundary(int i) = 0;
  //@}

  /** @name Other functions */
  //@{
  /// Empty virtual default destructor.
  virtual ~PrNestedTriangulation() {}

  /// Assuming that node i belongs to V^j\V^{j-1} and that
  /// jlev >= 1, find i's two parent nodes in \f$V^{j-1}\f$ p1 and p2.
  void getParents(int i, int jlev, int& p1, int& p2);

  // print routines
  //@}
};

/*>PrNestedTriangulation-syntax: */

/*Class:PrNestedTriangulation

Name:              "PrNestedTriangulation"
Syntax:	           @PrNestedTriangulation-syntax
Keywords:
Description:       This class implements a set of nested triangulations.
                   The methods in this class are used by PrFilterBank.
Member functions:
                   "getFinestLevel()" --\\
                   return the finest level in the hieararchy, i.e. k
                   if the vertices are arranged
                   V^0 \subset V^1 \subset ... \subset V^k.
                 
                   "getNumNodes(int jlev)" --\\
                   return the number of nodes in V^j, j=0,1,...,k.
                 
                   "getX(int i)" --\\
                   return the x coeff. of the i-th node in the triangulation.

                   "getY(int i)" --\\
                   return the y coeff. of the i-th node in the triangulation.

                   "getZ(int i)" --\\
                   return the z coeff. of the i-th node in the triangulation.

                   "setX(int i, double& x)" --\\
                   set the x coefficient of the i-th node in the triangulation.

                   "setY(int i, double& y)" --\\
                   set the y coefficient of the i-th node in the triangulation.

                   "setZ(int i, double& z)" --\\
                   set the z coefficient of the i-th node in the triangulation.

                   "getNeighbours(int i, int jlev, vector<int>& neighbours)"
                   --\\
                   Return the indices of the neighbours of the i-th node
                   with respect to V^j. This routine assumes that
                   node i also belongs to V^j. The order of the neighbours
                   returned will be
                   1. any anticlockwise order if i is an interior node
                   2. the unique anticlockwise order if i is a boundary node.

                   "isBoundary(int i)" --\\
                   is the i-th node a boundary node or interior node?

Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Feb. 98
*/

#endif // PRNESTEDTRIANGULATION_H
