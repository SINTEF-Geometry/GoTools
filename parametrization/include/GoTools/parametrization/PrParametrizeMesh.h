#ifndef PRPARAMETRIZEMESH_H
#define PRPARAMETRIZEMESH_H

#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrSubTriangulation.h"
#include "GoTools/parametrization/PrParamTriangulation.h"
#include "GoTools/parametrization/PrParametrizeInt.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
//#include "GoTools/parametrization/PrDijkstra.h"

#include <vector>
using std::vector;
#include <set>
using std::set;
#include <queue>

using std::shared_ptr;

/** PrParametrizeMesh -  Short description.
 * Detailed description.
 */

class PrParametrizeMesh {

  PrTriangulation_OP* mesh_;
  PrTriangulation_OP* basemesh_;

  PrParametrizeInt*  pi_;
  PrParametrizeBdy*  pb_;
//Dijkstra*          distance_finder_;

  double convTolerance;

public:

  /// Default constructor.
  PrParametrizeMesh();
  /// Default destructor.
  virtual ~PrParametrizeMesh();

  void attach(PrTriangulation_OP* mesh, PrTriangulation_OP* basemesh) 
    {mesh_ = mesh; basemesh_ = basemesh;}

  PrParamTriangulation* parametrize();

  void makeSubTriangulationFromTriangle(int i, 
					PrSubTriangulation& subtri,
					vector<int>& corners);

  /** Given a "triangulation" and a (ordered) set of "nodes",
   * this routine returns a "polygon", consisting of those
   * nodes (ordered) who form the shortest paths between
   * the given nodes in the triangulation
   *
   * i.e.: given a triangle (or arbitrary polygon) in the
   * coarse mesh, this algorithm finds the corresponding
   * triangle (polygon) in the original (fine) mesh
   */
  bool makePolygon( PrTriangulation_OP& triangulation, 
		    const vector<int>&   nodes, 
		    vector<int>&         polygon);
  bool makePath( PrTriangulation_OP& triangulation, 
		 int                 n1, 
		 int                 n2,
		 vector<int>&   path);

  /** Make a subtriangulation of triang along the polygon induced by the 
   * shortest paths between the nodes in nodes. Use the same nodeset ?
   */
  bool makeSubTriangulation ( PrTriangulation_OP&       triangulation, 
			      const vector<int>&   nodes, 
			      PrSubTriangulation&       subtriangulation);

  /// Make a subtriangulation of triang inside polygon. Assume polygon to be a 
  /// closed ring of neighbouring nodes
  bool makeSubTriangulationInsidePolygon ( PrTriangulation_OP& triangulation, 
					   const vector<int>&  polygon, 
					   PrSubTriangulation& subtriangulation);
 
  /** Assume the vertices in polygon is a closed oriented polygon on
   * triangulation. Returns the indices of the vertices interior to
   * polygon. If the orientation is reversed it returns the exterior
   * If the polygon is open and does not split triangulation all 
   * nodes will be found
   */
  bool getConnectedNodes( PrTriangulation_OP&       triangulation, 
			  const vector<int>&   polygon, 
			  set<int>&            nodes);

  /// Find the third node of the triangle with nodes n1, n2 in anti-clockwise 
  /// order (i.e. the triangle to the "left")
  int getLeftNode (PrTriangulation_OP& triangulation,
		   int n1,
		   int n2);

  void parametrizeSubTriangulation(std::shared_ptr<PrSubTriangulation> sub_tri,
				   vector<int>& corners);
  void getInteriorNeighbours( int                v,
			      const vector<int>& neighbours,
			      const vector<int>& boundary,
			      vector<int>&       new_nbrs );
  bool edgeInPolygon( int n1, int n2, const vector<int>& polygon );

  /** Casts a ray from the vertex "vs" into the direction of 
   * "angle". angle = 0 is the direction to the ccw node in
   * vs's leading triangle (remark: that is also the first
   * node in the result of the "getNeighbours" routine).
   * The ray is cast as long as its length does not exceed "dMax" 
   * and no triangle in "T" is visited. The triangles visited are 
   * stored in "t_path". If the path does not hit a triangle in "T", 
   * the return value is -1, otherwise it is the triangle index
   * "v_path" contains the vertices of the path.
   */
  int castRay( PrTriangulation_OP& triangulation, 
	       int vs, double angle, double dMax, const set<int>& T, 
	       vector<int>& t_path, vector<Vector3D>& v_path);

  /// Determine barycentric coordinates (r,s,t) of p4 w.r.t (p1,p2,p3)
  /// after flattening
  bool getNextV( const Vector3D& p1, const Vector3D& p2, 
		 const Vector3D& p3, const Vector3D& p4, 
		 double a, double b, bool f, double& c);

  /** runs the Dijkstra algorithm from source point "vs" to
   * destination point "vd" and determines the iso-distance-line of vd.
   * The vertices of that isoline are stored in "v_path", whereas
   * the indices of the triangles visited by that line are
   * stored in "t_path".
   * Return value is the distance from "vs" to "vd".
   */
   static double getIsoline(const PrTriangulation_OP& triangulation, 
			    int vs, 
			    int vd, 
			    vector<int>& t_path, 
			    vector<Vector3D>& v_path);

  /** Converts the direction indicated by "theta" around the node
   * "v" into a (double,int) pair. The "int" value is the index of
   * the triangle around "v" in direction "theta", the "double"
   * value is the barycentric coordinate of the intersection
   * point of a ray in direction "theta" with the edge opposite
   * to "v" in the triangle w.r.t to the vertices of that edge.
   */
   static void convertAngle(const PrTriangulation_OP& triangulation, 
			    double theta, int v, double& bc, int& t) ;
};


#endif // PRPARAMETRIZEMESH_H

