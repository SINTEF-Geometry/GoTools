//===========================================================================
//                                                                           
// File: EdgeVertex.h
//                                                                           
// Created: October 20, 2009
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _EDGEVERTEX_H
#define _EDGEVERTEX_H

#include <vector>
#include <memory>
#include "GoTools/utils/config.h"

namespace Go
{

    class ftEdge;
    class ftSurface;
    class Body;
    class Vertex;

//===========================================================================
/** Represents the radial edge structure in a non-manifold.
 *  Edges in a volume model setting. This entity corresponds to Vertex
 *  in the surface model setting.
 *
 */
//===========================================================================

class EdgeVertex
{
 public:
  /// Constructor given a number of adjacent half edges
  EdgeVertex(std::vector<ftEdge*> edges);

  /// Constructor given one half edge
  EdgeVertex(ftEdge* edge);

  /// Destructor
  ~EdgeVertex();

  /// Add another half edge to this radial edge. Used in topology
  /// build for volume models
  void addEdge(ftEdge *edge);
  
  /// Remove a half edge from this radial edge. Used in topology
  /// modifications for volume models
  void removeEdge(ftEdge *edge);

  /// Join the content of two edge vertex instances
  void addEdgeVertex(EdgeVertex* other);

  /// Returns all edges meeting in this vertex including twins
  std::vector<ftEdge*> allEdges() const;

  /// Returns all edges belonging to a given body
  /// meeting in this vertex including twins
  std::vector<ftEdge*> allEdges(Body *bd);

  /// Returns all geometrically unique edges meeting in this
  /// vertex. One edge is return for a pair of twins.
  std::vector<ftEdge*> uniqueEdges();

  /// Returns all geometrically unique edges belonging to a given body
  /// meeting in this vertex. One edge is return for a pair of twins.
  std::vector<ftEdge*> uniqueEdges(Body *bd);

  /// Number of unique edges meeting in this radial edge, twin information
  /// is counted only once. The number coincides with the number of
  /// bodies meeting in this edge
  int nmbUniqueEdges() 
  {
    return (int)edges_.size();
  }

  /// Number of unique edges of one body meeting in this radial edge, 
  /// twin information is counted only once. 
  int nmbUniqueEdges(Body *bd) const;

  /// Fetch a specified half edge, twins are counted once and an arbitrary
  /// member of a twin pair is returned
  ftEdge* getEdge(int idx)
  {
    return edges_[idx].first;
  }

  /// Check if a half edge is collected in this radial edge
  bool hasEdge(ftEdge *edge) const;

  /// Check if a half edge is collected in this radial edge and has no
  /// twin information
  bool hasEdgeSingle(ftEdge *edge) const;

  /// Get all adjacent faces
  std::vector<ftSurface*> getAdjacentFaces() const;

  /// Get all adjacent faces belonging to a given body
  std::vector<ftSurface*> getAdjacentFaces(Body *bd) const;

  /// Get all adjacent bodies
  std::vector<Body*> getAdjacentBodies() const;

  /// Reset twin organization
  void organizeTwins();

  /// Reorganize edges to get better pairs of twins
  void reOrganize();

  void disconnectTwin(ftEdge* e1, ftEdge *e2);

  /// Split edge vertex
  void splitAtVertex(shared_ptr<Vertex> v1,
		     shared_ptr<Vertex> v2, 
		     shared_ptr<Vertex> split_vx);

  /// Average coefficients of all spline surfaces meeting in this
  /// radial edge
  void averageSplineEdges(double eps);

  bool checkRadialEdgeTopology();

 private:
  /// Half edges collected in this radial edge, twins are represented in 
  /// pairs
  std::vector<std::pair<ftEdge*,ftEdge*> > edges_;
};

} // namespace Go


#endif // _EDGEVERTEX_H_


