/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRPLANARGRAPH_OP_H
#define PRPLANARGRAPH_OP_H

#include "GoTools/parametrization/PrExplicitConnectivity.h"
#include "GoTools/parametrization/PrPGNode.h"

/*<PrPlanarGraph_OP-syntax: */

/** PrPlanarGraph_OP - This class represents a planar graph
 * of points in \f$R^3\f$ with or without parameter points in \f$R^2\f$.
 * It implements the virtual functions in PrOrganizedPoints.
 * This kind of planar graph can be parametrized by PrParametrize.
 */
class PrPlanarGraph_OP : public PrExplicitConnectivity
{
private:
  std::vector<PrPGNode> node_;
      // array of nodes indexed from 0 to node_.size() - 1.
  std::vector<int> adj_;
     // adjacency list, giving indices of neighbours,
     // indexed from 0 to node_(node_.size()-1).end().

public:
  /// Empty default constructor
  PrPlanarGraph_OP() {}

  /** Constructor.
   * Construct the PrPlanarGraph_OP from an array of nodes and
   * an array of adjacency lists.
   * The length of the arrays xyz_nodes and uv_nodes should be 3*npts
   * and 2*npts respectively.
   * The length of the array end should be npts.
   * The length of the array adj should be end[npts].
   * If node 1 is an interior node, its neighbours are
   * adj[1],...,adj[end[1]] in any anticlockwise sequence.
   * If node 1 is a boundary node, its neighbours are
   * adj[1],...,adj[end[1]-1] in the unique anticlockwise sequence.
   * and adj[end[1]] = 0 to indicate that node 1 is on the boundary.
   * For i=2,...,n:
   * If node i is an interior node, its neighbours are
   * adj[end[i-1]],...,adj[end[i]] in any anticlockwise sequence.
   * If node i is a boundary node, its neighbours are
   * adj[end[i-1]],...,adj[end[i]-1] in the unique anticlockwise sequence.
   * and adj[end[i]] = 0 to indicate that node i is on the boundary.
   * This is the data structure proposed by Cline and Renka.
   */
  PrPlanarGraph_OP(int npts, double* xyz_nodes, double* uv_nodes, 
                   int* end, int* adj);

  /** Constructor.
   * Construct the PrPlanarGraph_OP from an array of nodes and
   * an array of adjacency lists. Like the previous constructor,
   * only we set the uv's to zero.
   */
  PrPlanarGraph_OP(int npts, double* xyz_nodes, int* end, int* adj);

  /// Constructor. Construct the PrPlanarGraph_OP from a PrOrganizedPoints class
  PrPlanarGraph_OP(PrOrganizedPoints& op);

  /// Empty destructor
  ~PrPlanarGraph_OP() {}

  // derived from base class

  /// Return the number of nodes in the graph.
    virtual int       getNumNodes() const {return (int)node_.size(); }

  /// Return the i-th node in the graph.
  virtual Vector3D get3dNode(int i) const {return node_[i].pnt(); }

  /// Doesn't really do anything. But needed anyway.
  virtual void      set3dNode(int i, const Vector3D& p) {;}

  /** Return the indices of the neighbours of the i-th node in:
   * 1. any anticlockwise order if i is an interior node
   * 2. the unique anticlockwise order if i is a boundary node.
   */
  virtual void      getNeighbours(int i, vector<int>& neighbours) const;

  /// Is i a boundary node or interior node?
  virtual bool      isBoundary(int i) const;

  /// Return the u parameter value of the i-th node.
  virtual double getU(int i) const {return node_[i].u(); };

  /// Return the v parameter value of the i-th node
  virtual double getV(int i) const {return node_[i].v(); };

  /// Reset the u parameter value of the i-th node.
  virtual void setU(int i, double u) {node_[i].u() = u; };
 
  /// Reset the v parameter value of the i-th node.
  virtual void setV(int i, double v) {node_[i].v() = v; };

  //print and scan routines
  void print(std::ostream& os);
  void scan(std::istream& is);
  /// alternative scan, reading the end array after points
  void scan2(std::istream& is);
};

/*>PrPlanarGraph_OP-syntax: */

/*Class:PrPlanarGraph_OP

Name:              PrPlanarGraph_OP
Syntax:	           @PrPlanarGraph_OP-syntax
Keywords:
Description:       This class represents a planar graph
                   of points in R^3 with or without parameter points in R^2.
                   It implements the virtual functions in PrOrganizedPoints.
                   This kind of planar graph can be parametrized
                   by PrParametrizer.
Member functions:
                   "getNumNodes()" --\\
                   Return the number of nodes in the graph.
                 
                   "get3dNode(int i)" --\\
                   Return the i-th node in the graph.

                   "getNeighbours(int i, PrListInt& neighbours)" --\\
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

Constructors:
Files:
Example:

| main()
| {
|   s_o << "Here is a planar graph:\n";
|   s_o << "\n";
|   s_o << "       7---6 \n";
|   s_o << "      /|   | \n";
|   s_o << "     / /   | \n";
|   s_o << "    / |    | \n";
|   s_o << "   2  /    5 \n";
|   s_o << "   | |    /  \n";
|   s_o << "   | 4   /   \n";
|   s_o << "   |/ \\ /    \n";
|   s_o << "   1---3     \n";
| 
|   double *xyz_nodes;
|   int  *end, *adj;
|   int npts=7;
|   xyz_nodes = new double[3*npts];
|   xyz_nodes[0] = 0.0; xyz_nodes[1] = 0.0; xyz_nodes[2] = 4.0; 
|   xyz_nodes[3] = 0.0; xyz_nodes[4] = 1.0; xyz_nodes[5] = 0.0; 
|   xyz_nodes[6] = 1.0; xyz_nodes[7] = 0.0; xyz_nodes[8] = 0.0; 
|   xyz_nodes[9] = 0.5; xyz_nodes[10] = 0.5; xyz_nodes[11] = -1.0; 
|   xyz_nodes[12] = 2.0; xyz_nodes[13] = 1.0; xyz_nodes[14] = 0.0; 
|   xyz_nodes[15] = 2.0; xyz_nodes[16] = 2.0; xyz_nodes[17] = 0.0; 
|   xyz_nodes[18] = 1.0; xyz_nodes[19] = 2.0; xyz_nodes[20] = 0.0; 
|   end = new int[npts];
|   end[0] = 4; end[1] = 7; end[2] = 11; end[3] = 14;
|   end[4] = 17; end[5] = 20; end[6] = 24;
|   adj = new int[end[npts-1]];
|   adj[0] = 3; adj[1] = 4; adj[2] = 2; adj[3] = 0;
|   adj[4] = 1; adj[5] = 7; adj[6] = 0;
|   adj[7] = 5; adj[8] = 4; adj[9] = 1; adj[10] = 0;
|   adj[11] = 3; adj[12] = 7; adj[13] = 1;
|   adj[14] = 6; adj[15] = 3; adj[16] = 0;
|   adj[17] = 7; adj[18] = 5; adj[19] = 0;
|   adj[20] = 2; adj[21] = 4; adj[22] = 6; adj[23] = 0;
| 
|   // make a planar graph
|   PrPlanarGraph_OP pl_graph(npts,xyz_nodes,end,adj);
|   delete [] xyz_nodes; delete [] end; delete [] adj; 
|   s_o << "The data of the planar graph are \n";
|   pl_graph.print(s_o);
|   s_o << "The information concerning the graph is \n";
|   pl_graph.printInfo(s_o);
| 
|   s_o << "The nodes of the graph are \n";
|   pl_graph.printXYZNodes(s_o);
|   s_o << "The edges of the graph are \n";
|   pl_graph.printXYZEdges(s_o);
|   s_o << "The faces of the graph are \n";
|   pl_graph.printXYZFaces(s_o);
| 
|   return 0;
| }

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Aug. 97
*/

#endif // PRPLANARGRAPH_OP_H
