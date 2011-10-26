/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 2000 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRTHIN_H
#define PRTHIN_H

#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrHeap.h"

/*<PrThin-syntax: */

/** PrThin - This class implements an algorithm for thinning a triangulation.
 */ 
class PrThin
{
protected:

  PrTriangulation_OP* trn_;
  std::vector<PrNode>* nlist_ptr; // the nodes in trn_
  std::vector<PrTriangle>* tlist_ptr; // the triangles in trn_

  PrHeap* heap_;

  double error_;
  int steps_;

  void makeHeap();
  double findError(int i);

public:
  /// Default constructor
  PrThin() {error_ = 0.0; }
  /// Empty destructor
  virtual ~PrThin();

  /// Set the triangulation.
  void attach(PrTriangulation_OP* trn);

  /// Set the tolerated error
  void setError(double error) {error_ = error; }

  /// Set number of steps used for thinning
  void setStepNumber(int n) {steps_ = n; }

  /** Find nearest neighbour j to an interior node i (in 3D) and
   * return the triangle tr which is left of the directed edge (i,j).
   * Note that j is then tr.getAnticlockwiseNode(i).
   */
  void findNearestNeighbour(int i, int& tr);

  /// Remove the node i and retriangulate by collapsing one of its edges.
  void halfEdgeCollapse(int i);

  /// Remove the node i from the triangulation and update heap.
  /// It is assumed that i is not in the heap
  void removePoint(int i);

  /// Thin the points, until the error of further removals is above the threshold error_.
  void thin();

  void printHeap();
};

/*>PrThin-syntax: */

/*Class:PrThin

Name:              PrThin
Syntax:	           @PrThin-syntax
Keywords:
Description:       This class implements an algorithm for thinning
                   a triangulation.
Member functions:
                   "attach(PrPrTriangulation_OP& tr)" --\\
                   Set the triangulation.

                   "halfEdgeCollapse(int i)" --\\
                   Remove the node i and retriangulate by collapsing
                   one of its edges.

Constructors:
Files:
Example:


See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Sep. 2000
*/

#endif // PRTHIN_H
