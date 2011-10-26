/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 2000 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRLEVELTRIANGULATION_OP_H
#define PRLEVELTRIANGULATION_OP_H

#include "GoTools/parametrization/PrExplicitConnectivity.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrNestedNode.h"
#include "GoTools/parametrization/PrTriangle.h"

/*<PrLevelTriangulation_OP-syntax: */

/** PrLevelTriangulation_OP - This class represents a triangulation
 * of points in \f$R^3\f$ with or without parameter points in \f$R^2\f$.
 * It implements the virtual functions in PrOrganizedPoints.
 * This kind of triangulation can be parametrized by PrParametrize.
 */
class PrLevelTriangulation_OP : public PrExplicitConnectivity
{
private:
  vector<PrNestedNode>* const node_;
  vector<PrTriangle> triangle_;
  const int level_;
  const int numNodes_;  // 0 <= numNodes_ < node->size()

  //int getNghrTriangle(int n1, int n2, vector<int>& tlist);
  //void buildTopology();

public:
  /// Empty default constructor.
  PrLevelTriangulation_OP()
      : node_(NULL), level_(0), numNodes_(0)
    {}

  /// Constructor.
  /// \param node pointer to a vector of nodes.  The vector will not be
  ///             copied internally.
  /// \param nt vector or triangles.  Will be copied internally
  /// \param level specify the level represented by this triangulation
  PrLevelTriangulation_OP(vector<PrNestedNode>* node,
                          vector<PrTriangle>& nt, int level = 0);

  /// Empty default destructor.
  ~PrLevelTriangulation_OP() {}

  /** @name Derived from base class */
  //@{

  /// Return the number of nodes in the graph.
  virtual int        getNumNodes() const {return numNodes_; }

  /// Return the i-th node in the graph.
  virtual Vector3D get3dNode(int i) const {return ((*node_)[i]).pnt(); }
  virtual void       set3dNode(int i, const Vector3D& p) 
     {(*node_)[i].x() = p.x(); (*node_)[i].y() = p.y();
      (*node_)[i].z() = p.z();}
  /** Return the indices of the neighbours of the i-th node in:
   *      1. any anticlockwise order if i is an interior node
   *      2. the unique anticlockwise order if i is a boundary node.
   */
  virtual void       getNeighbours(int i, vector<int>& neighbours) const;

  /// Is i a boundary node or interior node?
  virtual bool       isBoundary (int i) const;

  /// return the u parameter value of the i-th node.
  virtual double getU(int i) const {return (*node_)[i].u(); };

  /// return the v parameter value of the i-th node.
  virtual double getV(int i) const {return (*node_)[i].v(); };

  /// reset the u parameter value of the i-th node. 
  virtual void setU(int i, double u) {(*node_)[i].u() = u; };

  /// reset the v parameter value of the i-th node. 
  virtual void setV(int i, double v) {(*node_)[i].v() = v; };

    virtual int  findNumFaces() const {return (int)triangle_.size(); }
  //@}


  /** @name Other functions */
  //@{

  /// Grab internal structure for doing thinning
  inline PrNestedNode& getPrNestedNode(int i) {return (*node_)[i];}
  /// Grab internal structure for doing thinning
  inline PrTriangle& getPrTriangle(int i) {return triangle_[i];}

  /// Print out the triangles of the graph, useful for plotting.
  /// If num = 1, print first the number of triangles.
  void printXYZTriangles(ostream& os, bool num = false);
  /// Print out the uv triangles of the graph, useful for plotting.
  /// If num = 1, print first the number of triangles.
  void printUVTriangles(ostream& os, bool num = false);
  /// print whole class
  void print(ostream& os);
  /// print without neighbourhood information
  void printRawData(ostream& os);
  /// for converting uv to xyz with z=0
  void printUV(ostream& os);
  //@}
};

/*>PrLevelTriangulation_OP-syntax: */

/*Class:PrLevelTriangulation_OP

Name:              PrLevelTriangulation_OP
Syntax:	           @PrLevelTriangulation_OP-syntax
Keywords:
Description:       This class represents a triangulation
                   of points in R^3 with or without parameter points in R^2.
                   It implements the virtual functions in PrOrganizedPoints.
                   This kind of triangulation can be parametrized
                   by PrParametrize.
Member functions:
                   "getNumNodes()" --\\
                   Return the number of nodes in the graph.
                 
                   "get3dNode(int i)" --\\
                   Return the i-th node in the graph.

                   "getNeighbours(int i, PrIntList& neighbours)" --\\
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

| #include "GoTools/parametrization/PrLevelTriangulation_OP.h"
| #include "GoTools/parametrization/PrPlanarGraph_OP.h"
| 
| main()
| {
|   s_o << "Here is an explicit triangulation:\n";
|   s_o << "\n";
|   s_o << "   3-------4 \n";
|   s_o << "   |\\      | \n";
|   s_o << "   |\\\\     | \n";
|   s_o << "   | |\\    | \n";
|   s_o << "   | | \\   | \n";
|   s_o << "   | |  \\  | \n";
|   s_o << "   | 5-  \\ | \n";
|   s_o << "   |/  \\--\\| \n";
|   s_o << "   1-------2 \n";
| 
|   double *points;
|   int  *triangles;
|   int numpnts=5, numtrs=4;
|   points = new double[3*numpnts];
|   points[0] = 0.0; points[1] = 0.0; points[2] = 0.0; 
|   points[3] = 1.0; points[4] = 0.0; points[5] = 0.0; 
|   points[6] = 0.0; points[7] = 1.0; points[8] = 0.0; 
|   points[9] = 1.0; points[10] = 1.0; points[11] = 1.0; 
|   points[12] = 0.25; points[13] = 0.25; points[14] = 0.5; 
|   triangles = new int[3*numtrs];
|   triangles[0] = 1; triangles[1] = 2; triangles[2] = 5;
|   triangles[3] = 5; triangles[4] = 3; triangles[5] = 1;
|   triangles[6] = 2; triangles[7] = 4; triangles[8] = 3;
|   triangles[9] = 2; triangles[10] = 3; triangles[11] = 5;
| 
|   PrLevelTriangulation_OP pr_triang(points,numpnts,triangles,numtrs);
|   delete [] points; delete [] triangles;
| 
|   s_o << "The data of pr_triang is \n";
|   pr_triang.print(s_o);
|   s_o << "The information concerning pr_triang is \n";
|   pr_triang.printInfo(s_o);
| 
|   s_o << "The nodes of pr_triang are \n";
|   pr_triang.printXYZNodes(s_o);
|   s_o << "The edges of pr_triang are \n";
|   pr_triang.printXYZEdges(s_o);
|   s_o << "The triangles of pr_triang are \n";
|   pr_triang.printXYZFaces(s_o);
| 
|   // make a planar graph from the triangulation
|   PrPlanarGraph_OP pl_graph(pr_triang);
|   s_o << "The data of pl_graph is \n";
|   pl_graph.print(s_o);
| 
|   s_o << "Here is a terahedron (a closed triangulation):\n";
|   s_o << "\n";
|   s_o << "   3-------4 \n";
|   s_o << "   |\\     /| \n";
|   s_o << "   | \\   / | \n";
|   s_o << "   |  \\ /  | \n";
|   s_o << "   |   \\   | \n";
|   s_o << "   |  / \\  | \n";
|   s_o << "   | /   \\ | \n";
|   s_o << "   |/     \\| \n";
|   s_o << "   1-------2 \n";
| 
|   numpnts=4; numtrs=4;
|   points = new double[3*numpnts];
|   points[0] = 0.0; points[1] = 0.0; points[2] = 0.0; 
|   points[3] = 1.0; points[4] = 0.0; points[5] = 0.0; 
|   points[6] = 0.0; points[7] = 1.0; points[8] = 0.0; 
|   points[9] = 0.0; points[10] = 0.0; points[11] = 1.0; 
|   triangles = new int[3*numtrs];
|   triangles[0] = 2; triangles[1] = 3; triangles[2] = 4;
|   triangles[3] = 3; triangles[4] = 1; triangles[5] = 4;
|   triangles[6] = 1; triangles[7] = 2; triangles[8] = 4;
|   triangles[9] = 1; triangles[10] = 3; triangles[11] = 2;
| 
|   PrLevelTriangulation_OP pr_triang2(points,numpnts,triangles,numtrs);
|   delete [] points; delete [] triangles;
| 
|   s_o << "The data of pr_triang2 is \n";
|   pr_triang.print(s_o);
|   s_o << "The information concerning pr_triang2 is \n";
|   pr_triang.printInfo(s_o);
| 
|   return 0;
| }

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Nov. 2000
*/

#endif // PRLEVELTRIANGULATION_OP_H
