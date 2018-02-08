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

  /// Disconnect twin edges collected into this edge vertex, also in this
  /// data structure. Called from ftEdge.
  void disconnectTwin(ftEdge* e1, ftEdge *e2);

  /// Update twin information for given edge
  void updateEdgeInfo(ftEdge* edge);

  /// Split edge vertex
  void splitAtVertex(shared_ptr<Vertex> v1,
		     shared_ptr<Vertex> v2, 
		     shared_ptr<Vertex> split_vx);

  /// Average coefficients of all spline surfaces meeting in this
  /// radial edge
  void averageSplineEdges(double eps);

  /// Debug functionality
  bool checkRadialEdgeTopology();

 private:
  /// Half edges collected in this radial edge, twins are represented in 
  /// pairs
  std::vector<std::pair<ftEdge*,ftEdge*> > edges_;
};

} // namespace Go


#endif // _EDGEVERTEX_H_


