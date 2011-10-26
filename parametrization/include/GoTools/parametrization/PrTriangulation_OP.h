/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRTRIANGULATION_OP_H
#define PRTRIANGULATION_OP_H

#include "GoTools/parametrization/PrExplicitConnectivity.h"
#include "GoTools/parametrization/PrNode.h"
#include "GoTools/parametrization/PrTriangle.h"
#include "GoTools/utils/ScratchVect.h"

/*<PrTriangulation_OP-syntax: */

/** PrTriangulation_OP -  This class represents a triangulation
 * of points in \f$R^3\f$ with or without parameter points in \f$R^2\f$.
 * It implements the virtual functions in PrOrganizedPoints.
 * This kind of triangulation can be parametrized by PrParametrize.
 */
class PrTriangulation_OP : public PrExplicitConnectivity
{
private:
  std::vector<PrNode> node_;
  std::vector<PrTriangle> triangle_;

  int getNghrTriangle(int n1, int n2, std::vector<int>& tlist);
  void buildTopology();

public:
  /// Default constructor.
  PrTriangulation_OP() {}

  /** Constructor.Construct the PrTriangulation_OP from an array of nodes and
   * an array of triangles.
   * The nodes array should contain x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,...
   * and the triangle array i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,....
   * where the i-th node is the 3D point (xi,yi,zi) and the
   * r-th triangle has the three nodes indexed ir, jr, and kr in
   * the node array. These nodes must be ordered anticlockwise (consistently
   * throughout the triangulation).
   * The length of the array "xyz_points" is 3 * np and
   * the length of the array "triangles" is 3 * nt.
   * The uv points for the nodes are set to zero.
   *
   * The function buildTopology is called which finds the three
   * neighbouring triangles of each triangle and the first triangle of each node.
   */
  PrTriangulation_OP(const double *xyz_points, int np, const int *triangles, int nt);

  /** Constructor.Construct the PrTriangulation_OP from an array of nodes and
   * an array of triangles.
   * The xyz_points array should contain x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,...
   * The uv_points array should contain u0,v0,u1,v1,u2,v2,u3,...
   * and the triangle array i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,....
   * where the i-th node is the 3D point (xi,yi,zi) and the
   * r-th triangle has the three nodes indexed ir, jr, and kr in
   * the node array. These nodes must be ordered anticlockwise (consistently
   * throughout the triangulation).
   * The length of the array "xyzpoints" is 3 * np and
   * The length of the array "uv_points" is 2 * np and
   * the length of the array "triangles" is 3 * nt.
   * 
   * The function buildTopology is called which finds the three
   * neighbouring triangles of each triangle and the first triangle of each node.
   */
  PrTriangulation_OP(const double *xyz_points, const double *uv_points,
                     int np, const int *triangles, int nt);

  /** Constructor. Construct the PrTriangulation_OP from a PrExplicitConnectivity class
   * ASSUMING the PrOrganizedPoints class is a triangulation (every
   * face has three vertices).
   */
  PrTriangulation_OP(PrExplicitConnectivity& op);
    ~PrTriangulation_OP() {}

  /** @name Derived from base class */
  //@{
    virtual int        getNumNodes() const  {return (int)node_.size(); }
  virtual Vector3D get3dNode(int i) const {return node_[i].point(); }
  virtual void       set3dNode(int i, const Vector3D& p) 
           {node_[i].x() = p[0]; node_[i].y() = p[1]; node_[i].z() = p[2];}

  /** Return the indices of the neighbours of the i-th node in:
   * 1. any anticlockwise order if i is an interior node
   * 2. the unique anticlockwise order if i is a boundary node.
   */
  virtual void       getNeighbours(int i, std::vector<int>& neighbours) const;
  virtual bool       isBoundary(int i) const ;
    virtual int  findNumFaces() const {return (int)triangle_.size(); }
  virtual double getU(int i) const {return node_[i].u(); };
  virtual double getV(int i) const {return node_[i].v(); };
  virtual void setU(int i, double u) {node_[i].u() = u; };
  virtual void setV(int i, double v) {node_[i].v() = v; };
  //@}

  /** @name Other functions */
  //@{
  /// Get all triangles that are incident with node 'i'.
  virtual void getTriangles(int i, Go::ScratchVect<int, 20> &triangles) const;

  // mesh modifications
  /// swap triangle t1 and t2
  bool swapTriangles(int t1, int t2);
  /// swap node n1 and n2
  void swapNodes(int n1, int n2);
  /// splits the common edge of triangles t1 and t2 and inserts the new
  /// vertex v there.
  bool splitTriangles (int t1, int t2, Vector3D& v);
    
  /// Split the boundary vertex indexed 'i' into two, and replace the edge
  /// between 'i' and 'j' with two boundary edges.  The requirements are that
  /// node 'i' should be on the boundary, and 'j' should share an edge with 'i'.
  bool splitVertex(int i, int j, std::vector<int> &new_nodes);

  // Grab internal structure for doing thinning
  /// Get reference to node indexed 'i'.
  inline PrNode& getPrNode(int i) {return node_[i];}
  inline const PrNode& getPrNode(int i) const {return node_[i];}
  /// Get reference to triangle indexed 'i'.
  inline const PrTriangle& getPrTriangle(int i) const {return triangle_[i];}
  /// Get reference to triangle indexed 'i'.    
  inline PrTriangle& getPrTriangle(int i) {return triangle_[i];}
  /// Get reference to all triangles in the triangulation.
  inline std::vector<PrTriangle>& getTriangleArray() {return triangle_;}
  /// Get reference to all nodes in the triangulation
  inline std::vector<PrNode>& getNodeArray() {return node_;}

  //print and scan routines

  /// Print out the triangles of the graph, useful for plotting
  /// If num = 1, print first the number of triangles.
  void printXYZTriangles(std::ostream& os, bool num = false) const;

  /// Print out the uv triangles of the graph, useful for plotting
  /// If num = 1, print first the number of triangles.
  void printUVTriangles(std::ostream& os, bool num = false) const;

  /// print whole class
  void print(std::ostream& os) const;
  /// scan whole class
  void scan(std::istream& is);
  /// print minimum data
  void printRawData(std::ostream& os) const;
  /// scan minimum data
  void scanRawData(std::istream& is);
  /// for converting uv to xyz with z=0
  void printUV(std::ostream& os) const;
  //@}
};

/*>PrTriangulation_OP-syntax: */

/*Class:PrTriangulation_OP

Name:              PrTriangulation_OP
Syntax:	           @PrTriangulation_OP-syntax
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

| #include "GoTools/parametrization/PrTriangulation_OP.h"
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
|   PrTriangulation_OP pr_triang(points,numpnts,triangles,numtrs);
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
|   PrTriangulation_OP pr_triang2(points,numpnts,triangles,numtrs);
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
Date:              Aug. 97
*/

#endif // PRTRIANGULATION_OP_H
