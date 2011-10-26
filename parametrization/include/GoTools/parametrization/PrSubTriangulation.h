/*****************************************************************************/
/*                                                                           */
/* (c) Copyright 2000 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRSUBTRIANGULATION_H
#define PRSUBTRIANGULATION_H

#include "GoTools/parametrization/PrExplicitConnectivity.h"
#include "GoTools/parametrization/PrNode.h"
#include "GoTools/parametrization/PrTriangle.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"

#include <vector>
using std::vector;
#include <map>
using std::map;

/*<PrSubTriangulation-syntax: */

/** PrSubTriangulation -  This class represents a sub-triangulation of another
 * triangulation, using a permutation vector to access
 * the nodes and triangles of the other triangulation
 * It implements the virtual functions in PrOrganizedPoints.
 * This kind of triangulation can be parametrized by PrParametrize.
 * 
 * The internal indices of the vertices are [0..numBdy_-1]
 * for the boundary nodes and [numBdy_..numBdy_+numInt_-1]
 * for the interior nodes. They are referred to as "new"
 * indices, whereas the indices in the associated
 * PrTriangulation_OP are called "old" indices. Therefore,
 * "new2old_" provides a mapping from [0..numBdy_+numInt_-1]
 * to the "old" indices and "old2new_" is its inverse.
 * 
 * We assume the boundary vertices to be ordered anticlockwise.
 */
class PrSubTriangulation : public PrExplicitConnectivity
{
private:
  int numBdy_;                  // number of boundary vertices
  int numInt_;                  // number of interior vertices
  vector<int> new2old_;         // mapping from the "new" to the "old" indices
  map<int, int> old2new_;       // mapping from the "old" to the "new" indices
  PrTriangulation_OP* g_;       // the data itself

public:
  /// Default constructor
  PrSubTriangulation() {}
  /** Constructor. Construct the PrSubTriangulation from an array with the indices
   * of the "boundary" nodes and another carrying the indices of the
   * "interior" nodes. All indices are relative to the indexing in
   * the associated PrTriangulation_OP "graph"
   */
  PrSubTriangulation(PrTriangulation_OP* graph, 
		     vector<int> boundary,
		     vector<int> interior);
  /// Empty destructor
  ~PrSubTriangulation() {}

  /// Set the top-level triangulation to 'graph'
  void attach(PrTriangulation_OP* graph) {g_ = graph;}
    
  /// Initialize the PrSubTriangulation.  The attach() function should already
  /// have been called in order to define the PrTriangulation_OP that constitutes
  /// the base triangulation.  The 'boundary' vector specifies the indexes of the
  /// boundary nodes, and the 'interior' vector specifies the indexes of the 
  /// interior nodes.  All indices are relative to the indexing in the associated 
  /// PrTriangulation_OP.
  void initialize(vector<int> boundary, vector<int> interior);

  /** @name Derived from base class */
  //@{
  virtual int getNumNodes() const {return numBdy_ + numInt_;}
  virtual int findNumBdyNodes() const {return numBdy_;}
  virtual Vector3D get3dNode(int i) const
    {return g_->get3dNode(new2old_[i]);}
  virtual void set3dNode(int i, const Vector3D& p)
    { g_->set3dNode(new2old_[i], p); };
  virtual void getNeighbours(int i, vector<int>& neighbours) const;
  virtual bool isBoundary(int i) const {return (i < numBdy_);}

  virtual double getU(int i) const {return g_->getU(new2old_[i]);}
  virtual double getV(int i) const {return g_->getV(new2old_[i]);}
  virtual void setU(int i, double u) {g_->setU(new2old_[i], u);}
  virtual void setV(int i, double v) {g_->setV(new2old_[i], v);}
  //@}

  /** @name Other routines */
  //@{
  /// get the global index of the node with local index 'i'.
  int getGlobalIndex(int i) const {return new2old_[i];}
  /// get the local index of the node with global index 'i'.
  int getLocalIndex(int i) const {
    return (*(old2new_.find(i))).second;
  }

  /// Check if the triangle indexed 'i' in the global triangulation also
  /// exists in this sub-triangulation.  This is true if all three nodes of the
  /// triangle can be found in the sub-triangulation.
  bool triangleInThis(int i) const;

  /// Get all the neighbour triangles in the global triangulation of a node 'i'.
  /// (i is global index).
  void getNeighbourTriangles(int i, vector<int>& neighbours) const;

  /// Get all indices of triangles made up of interior points in the global 
  /// triangulation, as well as those made up of boundary/interior points and
  /// equally belong to the sub-triangulation.
  void getTriangleIndices(vector<int>& indices) const;
 
  //@}
};

/*>PrSubTriangulation-syntax: */

/*Class:PrSubTriangulation

Name:              PrSubTriangulation
Syntax:	           @PrSubTriangulation-syntax
Keywords:
Description:       This class represents a sub-triangulation of another
                   triangulation, using a permutation vector to access
		   the nodes and triangles of the other triangulation
                   It implements the virtual functions in PrOrganizedPoints.
                   This kind of triangulation can be parametrized
                   by PrParametrize.

		   The internal indices of the vertices are [0..numBdy_-1]
		   for the boundary nodes and [numBdy_..numBdy_+numInt_-1]
		   for the interior nodes. They are referred to as "new"
		   indices, whereas the indices in the associated
		   PrTriangulation_OP are called "old" indices. Therefore,
		   "new2old_" provides a mapping from [0..numBdy_+numInt_-1]
		   to the "old" indices and "old2new_" is its inverse.
		   
		   We assume the boundary vertices to be ordered
		   anticlockwise.

Member functions:
                   "getNumNodes()" --\\
                   Return the number of nodes in this subtriangulation.
                 
                   "get3dNode(int i)" --\\
                   Return the i-th node of this subtriangulation.

                   "getNeighbours(int i, vector<int>& neighbours)" --\\
                   Return the indices of the neighbours of the i-th node in:
                   1. any anticlockwise order if i is an interior node
                   2. the unique anticlockwise order if i is a boundary node.

                   "isBoundary(int i)" --\\
                   Is i a boundary node or interior node?

                   "getU(int i)" --\\
                   return the u parameter value of the i-th node.

                   "getV(int i)" --\\
                   return the v parameter value of the i-th node.

                   "setU(int i, double u)" --\\
                   reset the u parameter value of the i-th node.

                   "setV(int i, double v)" --\\
                   reset the v parameter value of the i-th node.

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Kai Hormann, SINTEF
Date:              Oct. 2000
*/

#endif // PRSUBTRIANGULATION_H
