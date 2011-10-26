/********************************************************************
 FILENAME    : PrPathTriangleSeq.h
 AUTHOR      : Valerie PHAM-TRONG, SINTEF
 DATE        : Mai 2002
 DESCRIPTION : Shortest path in a 3D sequence of triangles, 
               using unfolding tools.
 CHANGE LOG  :
*********************************************************************/

#ifndef _PRPATHTRIANGLESEQ_H
#define _PRPATHTRIANGLESEQ_H

#include <fstream>
#include <iostream>
#include <list>

#include "GoTools/parametrization/PrOrganizedPoints.h"
#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrGeodesics.h"
#include "GoTools/parametrization/PrInterface.h"

#include "GoTools/parametrization/PrTriangle.h"
//#include "GoTools/parametrization/PrDijkstra.h"
#include "GoTools/utils/Array.h"

/** EdgeType -  Short description.
* Detailed description.
*/
class EdgeType
{
 public:
    int vertex_[2];

    void initEdgeType();
    bool isVertex(int node);

};

/** UnfNodeType -  Short description.
* Detailed description.
*/
class UnfNodeType
{
 public:
    Go::Vector2D n_;
    int origine_;

    void printNode();
};

void printNode(vector<UnfNodeType>& v);

/** PathType -  Short description.
* Detailed description.
*/
class PathType
{
 public:
    std::vector<int> l_path_;
    std::vector<int> r_path_;
    std::vector<int> funnel_;

    void initPathType();

  /// Computes the funnel for the funnel algorithm
    void funnelOfPath(PathType prec);

  /// Modification of the left path in the funnel algorithm
    void modification_l_path_(int new_pt, const vector<UnfNodeType>& nodes_unf);

  /// Modification of the right path in the funnel algorithm
    void modification_r_path_(int new_pt, const vector<UnfNodeType>& nodes_unf);
    void print();
};


/**
 *  input: triangulation and triangle sequence. 
 *  output: 
 *  - shortest path in that sequence = polygonal line given by its 3D vertices.
 *  - list of pivot vertices, ie vertices of the triangulation 
 *        belonging to this shortest path.
 *  - length of this shortest path.
*/
vector<Vector3D> shortest_path_triangle_sequence(
    const PrTriangulation_OP& t, 
    int sce_vertex, int dest_vertex,
    const vector<int> &tr_seq3d, list<int>& list_pivots, 
    double& length);

/**
 *  input: triangulation and triangle sequence. 
 *  output: 
 *  - shortest path in that sequence = polygonal line given by its 3D vertices.
 *  - list of pivot vertices, ie vertices of the triangulation. 
 *        belonging to this shortest path.
 *  - length of this shortest path.
*/
void shortest_path_triangle_sequence(
    const PrTriangulation_OP& t, 
    const int& sce_vertex, const int& dest_vertex,
    vector<int> &tr_seq3d, list<int>& list_pivots, 
    std::vector<EdgeType> &edge_seq3d, std::vector<double> &ratio_on_edge,
    double& length);

/// Computes the 3D points of the path from their position (ratio) on
/// the edges of the sequence.
void path_3Dpts_from_ratio_on_edge(
    const PrTriangulation_OP& t, std::vector<Go::Vector3D>& sh_path,
    const std::vector<double>& ratio_on_edge, const std::vector<EdgeType>& edge_seq3d);

/// Computes a 3D point from its position (ratio) on an edge
Go::Vector3D pt3D_from_ratio_on_edge(double r, Go::Vector3D A, Go::Vector3D B);

/// Computes the 2D points of the path from their position (ratio) on
/// the edges of the sequence
void path_2Dpts_from_ratio_on_edge(
    const vector<UnfNodeType> &nodes_unf, vector<Vector2D>& path,
    const vector<double>& ratio_on_edge, const vector<EdgeType>& edge_seq2d);

/// Computes a 2D point from its position (ratio) on an edge
Vector2D pt2D_from_ratio_on_edge(double r, Vector2D a, Vector2D b);

/**  computes the 3d pivot list corresponding to the pivot in the unfolded sequence.
 *  The pivots are the vertices of the shortest path.
 *  The pivots in the unfolded sequence contain all the vertices of the path (source and dest. too).
 *  The 3d pivots do not contain the source and dest. vertices.
 */
list<int> pivots3D_from_pivots_in_unf_seq(
    const list<int> &list_pivots_unf, const vector<UnfNodeType> &nodes_unf);

/** At each interior pivot, the path is deviated from the straight line.
 * This function computes the deviation angles.
 */
list<double> deviation_from_pivots(
    const list<int> &list_pivots_unf, const vector<UnfNodeType>& nodes_unf);

/// absolute value
double abs_val(double a);

/// computes the shortes path in an ordered 2d seuqence of triangles, 
/// corresponding to a 3d seuqence that was unfolded
vector<double> shortest_path_triangle_sequence_2D(
    const PrTriangulation_OP &t, 
    const int sce_vertex, const int dest_vertex,
    const int sce_vertexUnf, int dest_vertexUnf,
    const std::vector<int>& tr_seq3d, std::vector<PrTriangle>& tr_seq_unf, 
    std::list<int>& list_pivots_unf, 
    std::vector<UnfNodeType> nodes_unf, double& length);

void print_list_int(std::list<int>& list_pivots_unf);
void print_deviation(std::list<double>& list_deviation);

double length_polygonal_path(
    const list<int> &list_pivots_unf, vector<UnfNodeType>& nodes_unf);

/// computes the intersection of the path, given by the pivots, with the
/// edge sequence, returns the ratio/position on each edge
std::vector<double> path_ratio_on_edge_from_pivots(
    std::list<int>& list_pivots_unf, std::vector<EdgeType>& edge_seq_unf, 
    const std::vector<UnfNodeType> &nodes_unf);

/** intersection between [ab] and [cd].
 * returns the position of the intersection on the oriented segment [cd] given by a ratio
 */
double segment_intersection(Go::Vector2D a, Go::Vector2D b, 
			    Go::Vector2D c, Go::Vector2D d);

/** computes i, the intersection between (ab) and (cd)
 * returns 1 if there is a unique solution,
 * 2 if the lines are the same, 0 if they are parallel
 */
int line_intersection( Go::Vector2D a, Go::Vector2D b, 
		       Go::Vector2D c, Go::Vector2D d, 
		       Go::Vector2D& i);

/** Computes the path from a source point to an edge of the 2d sequence:.
 * This path is given by:. 
 *   - a l_path_, polygonal chain of vertices to the left vertex of the edge.
 *   - a r_path_, polygonal chain of vertices to the right vertex of the edge.
 *   - funnel_, common part of the path. 
 */
void path_to_edge( 
    const int sce_vertex_unf, int i,
    std::vector<PathType>& p, const std::vector<UnfNodeType>& nodes_unf, 
    const std::vector<PrTriangle>& tr_seq_unf, std::vector<EdgeType>& edge_seq_unf);

/// given a triangle sequence, computes an edge sequence consisting of 
/// the common edges to two successiv triangles in the sequence
std::vector<EdgeType> edge_sequence(
    const std::vector<int> &tr_seq, const PrTriangulation_OP& t);

/// given a triangle sequence, computes an edge sequence consisting of 
/// the common edges to two successiv triangles in the sequence
std::vector<EdgeType> edge_sequence(std::vector<PrTriangle> tr_seq);

/// common edge between two triangles, with the nodes given in the order of t1
EdgeType common_edge(PrTriangle t1, PrTriangle t2);

/** unfolding of a sequence of triangle 3D from a triangulation.
 * returns a set of unfolded nodes, a sequence of unfolded triangles and 
 * a sequence of unfolded edges
 */
void unfolding_triangle_sequence(
    const PrTriangulation_OP& t, 
    const vector<int>& tr_seq3d, vector<PrTriangle>& tr_seq_unf,
    const vector<EdgeType>& edge_seq3d, vector<EdgeType>& edge_seq_unf,
    vector<UnfNodeType>& nodes_unf);             

/// first triangle unfolded with 1st vertex = (0,0), 2nd on the x-axis,
/// 3rd in the y>0 half plane 
void unfolding_first_triangle( const PrTriangulation_OP& t, 
			       const vector<int>& tr_seq, 
			       vector<PrTriangle>& tr_seq_unf,
			       vector<UnfNodeType>& nodes_unf);

/// associates the unfolded node to a vertex of the sequence
/// (for example for the source and destination vertices)
void unfolding_vertex(
    const PrTriangulation_OP& t, 
    const vector<int>& tr_seq3d, vector<PrTriangle>& tr_seq_unf,
    const int vertex, int& vertex_unf);

/// local frame of triangle ABC
void local_frame( 
    Go::Vector3D a, Go::Vector3D b, Go::Vector3D c,  
    Go::Vector3D& U, Go::Vector3D& V, Go::Vector3D& W);

/// axis and normal vector in the current triangle plane for the local frame
Go::Vector3D vector_W(
    Go::Vector3D a, Go::Vector3D b, Go::Vector3D c);

/// first axis/vector in the current triangle plane for the local frame
Go::Vector3D vector_U(const Go::Vector3D W);

/// second axis/vector in the current triangle plane for the local frame
Go::Vector3D vector_V(const Go::Vector3D W, const Go::Vector3D U);

/// local coordinates of a vector A in an axis frame given by (U, V, W)
Go::Vector3D local_coordinates(const Go::Vector3D a, 
			     const Go::Vector3D& U, const Go::Vector3D& V, const Go::Vector3D& W);

/// angle between two vectors AB and CD in the plane
double angle(const Go::Vector2D& a, const Go::Vector2D& b, 
	     const Go::Vector2D& c, const Go::Vector2D& d);

void rotation2D(const Go::Vector2D a, const double theta, 
		const Go::Vector2D c, Go::Vector2D& vect);

void axial_symmetry(const Go::Vector2D c, const Go::Vector2D a, 
		    const Go::Vector2D b, Go::Vector2D& v);

/// third vertex pf triangle t not belonging to edge e
int third_vertex(PrTriangle Tr, EdgeType e);

/// check if pt1 and pt2 are both on the same side of [ptA ptB]
int same_side(const Go::Vector2D pt1, const Go::Vector2D pt2, 
	      const Go::Vector2D ptA, const Go::Vector2D ptB);

double vect_prod2D(Go::Vector2D a, Go::Vector2D b, Go::Vector2D c, Go::Vector2D d);

void printUnfoldedSequence(std::ofstream& os, 
			   std::vector<UnfNodeType>& nodes_unf, std::vector<PrTriangle>& tr_seq_unf);

void printUnfoldedPath(std::ofstream& os, 
		       std::vector<UnfNodeType>& nodes_unf,     
		       std::list<int>& list_pivots_unf);

void printUnfoldedPath(std::ofstream& os, std::vector<Go::Vector2D> path);

void printUnfoldedObjects(std::ofstream& os, 
			  std::vector<UnfNodeType>& nodes_unf, std::vector<PrTriangle>& tr_seq_unf,     
			  std::list<int>& list_pivots_unf);

/// print the lengths of the edges of the triangles in the unfolded sequence
void print_edges_lengths(std::vector<PrTriangle> tr, std::vector<UnfNodeType> nodes_unf);
void print_triangle_sequence(std::vector<PrTriangle> tr, std::vector<UnfNodeType> nodes_unf);
void print_tr_vertices(std::vector<int>& tr, PrTriangulation_OP& t);
#endif

