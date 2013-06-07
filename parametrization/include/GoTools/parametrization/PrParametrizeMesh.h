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

#ifndef PRPARAMETRIZEMESH_H
#define PRPARAMETRIZEMESH_H

#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrSubTriangulation.h"
#include "GoTools/parametrization/PrParamTriangulation.h"
#include "GoTools/parametrization/PrParametrizeInt.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
//#include "GoTools/parametrization/PrDijkstra.h"

#include <vector>
#include <set>
#include <queue>

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
					std::vector<int>& corners);

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
		    const std::vector<int>&   nodes, 
		    std::vector<int>&         polygon);
  bool makePath( PrTriangulation_OP& triangulation, 
		 int                 n1, 
		 int                 n2,
		 std::vector<int>&   path);

  /** Make a subtriangulation of triang along the polygon induced by the 
   * shortest paths between the nodes in nodes. Use the same nodeset ?
   */
  bool makeSubTriangulation ( PrTriangulation_OP&       triangulation, 
			      const vector<int>&   nodes, 
			      PrSubTriangulation&       subtriangulation);

  /// Make a subtriangulation of triang inside polygon. Assume polygon to be a 
  /// closed ring of neighbouring nodes
  bool makeSubTriangulationInsidePolygon ( PrTriangulation_OP& triangulation, 
					   const std::vector<int>&  polygon, 
					   PrSubTriangulation& subtriangulation);
 
  /** Assume the vertices in polygon is a closed oriented polygon on
   * triangulation. Returns the indices of the vertices interior to
   * polygon. If the orientation is reversed it returns the exterior
   * If the polygon is open and does not split triangulation all 
   * nodes will be found
   */
  bool getConnectedNodes( PrTriangulation_OP&       triangulation, 
			  const std::vector<int>&   polygon, 
			  std::set<int>&            nodes);

  /// Find the third node of the triangle with nodes n1, n2 in anti-clockwise 
  /// order (i.e. the triangle to the "left")
  int getLeftNode (PrTriangulation_OP& triangulation,
		   int n1,
		   int n2);

  void parametrizeSubTriangulation(shared_ptr<PrSubTriangulation> sub_tri,
				   std::vector<int>& corners);
  void getInteriorNeighbours( int                v,
			      const std::vector<int>& neighbours,
			      const std::vector<int>& boundary,
			      std::vector<int>&       new_nbrs );
  bool edgeInPolygon( int n1, int n2, const std::vector<int>& polygon );

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
	       int vs, double angle, double dMax, const std::set<int>& T, 
	       std::vector<int>& t_path, std::vector<Vector3D>& v_path);

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
			    std::vector<int>& t_path, 
			    std::vector<Vector3D>& v_path);

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

