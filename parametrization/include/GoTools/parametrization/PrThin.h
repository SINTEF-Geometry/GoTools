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
