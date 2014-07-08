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

#ifndef _VERTEX_H
#define _VERTEX_H

#include "GoTools/utils/Point.h"
#include <vector>
#include <memory>

namespace Go
{

    class ftEdge;
    class ftSurface;
    class Body;

/// \brief The vertex class represents the vertex entity in
/// a boundary represented solid or face set

class Vertex
{
 public:
  /// Constructor. Give the geometric position of the vertex
    Vertex(Point vertex_point);

    /// Constructor. Give the geometric position of the vertex and
    /// associated edges
    Vertex(Point vertex_point, std::vector<ftEdge*> edges);

    /// Constructor. Give the geometric position of the vertex and
    /// one associated edge
    Vertex(Point vertex_point, ftEdge* edges);

    /// Constructor. Give one associated edge and an indication on
    /// which of the two vertices belonging to this edge should be
    /// constructed
     Vertex(ftEdge* edge, bool at_start);

     /// Destructor
    ~Vertex();

    /// Vertices belonging to two adjacent edges are represented
    /// as one entity. Used in topology build
    void joinVertex(shared_ptr<Vertex> other);

    /// Add a new edge to this vertex. Used in topology build
    void addEdge(ftEdge* edge);

    /// Remove an edge from this vertex. Used in connetion with topology 
    /// changes in the associated model
    void removeEdge(ftEdge* edge);

    /// Given an edge connected to this vertex, remove twin information
    /// about the edge in the vertex. Used in connetion with topology 
    /// changes in the associated model
    void disconnectTwin(ftEdge* edge);

    /// Returns all edges meeting in this vertex including twins
    std::vector<ftEdge*> allEdges() const;

    /// Returns all geometrically unique edges meeting in this
    /// vertex. One edge is return for a pair of twins.
    std::vector<ftEdge*> uniqueEdges();

    /// Returns all edges belonging to a given body. Twin edges are
    /// represented only as one edge
    std::vector<ftEdge*> uniqueEdges(Body *bd);

    /// Number of unique edges meeting in this vertex, twin edges
    /// are counted only once
    int nmbUniqueEdges()
    {
      return (int)edges_.size();
    }

    int nmbUniqueEdges(Body *bd);

    /// Get edges which are not associated a face
    std::vector<ftEdge*> freeEdges();

    /// Get the specified edge, twin edges are counted once and it is
    /// arbitrary which twin edge is returned
    ftEdge* getEdge(int idx)
    {
      return edges_[idx].first;
    }

    /// Get the geometrical position associated to this vertex
    Point getVertexPoint()
	{
	    return vertex_point_;
	}

    /// Set the geometrical position associated to this vertex
    void setVertexPoint(Point vertex_point)
    {
      vertex_point_ = vertex_point;
    }
      
    /// Get all faces meeting in this vertex
    std::vector<ftSurface*> faces() const;

    /// Get all faces belonging to a given body meeting in this vertex
    std::vector<ftSurface*> faces(Body *bd) const;

    /// Get all faces meeting in this vertex and the parameter value in
    /// the face corresponding to the vertex
    std::vector<std::pair<ftSurface*, Point> > getFaces();

    /// Get all faces belonging to a given body meeting in this vertex 
    /// and the parameter value in the face corresponding to the vertex
   std::vector<std::pair<ftSurface*, Point> > getFaces(Body *bd);

    /// Average corners of spline surfaces corresponding to this vertex
    void averageVertexPos();

    /// Get parameter of associated face corresponding to vertex
    Point getFacePar(ftSurface* face);

    /// Get all boides meeting in this vertex
    std::vector<Body*> getBodies();

    /// The distance between this vertex an another vertex
    double getDist(shared_ptr<Vertex> other_point)
	{
	    return vertex_point_.dist(other_point->getVertexPoint());
	}

    /// Check if the vertex is connected to the given edge
    bool hasEdge(ftEdge *edge) const;

    /// Check if the vertex is adjacent to the given face
    bool hasFace(ftSurface *face) const;

    /// Check if the vertex is connected to the given edge, and this edge
    /// is represented in the vertex with no twin
    bool hasEdgeSingle(ftEdge *edge) const;

    /// Check if two given edges meet in this vertex
    bool meetInVertex(ftEdge *e1, ftEdge *e2) const;

    /// Check if this vertex lies at a model boundary, i.e. is connected to
    /// edges with no twin
    bool isBoundaryVertex() const;

    /// Check if this vertex and the other vertex belongs to the same edge
    bool sameEdge(Vertex* other) const;


    /// Check if this vertex and the other vertex belongs to the same face
    bool sameFace(Vertex* other) const;

    /// Check if this vertex and the other vertex belongs to the same 
    /// underlying surface (one surface can give rise to several faces)
    bool sameUnderlyingSurface(Vertex* other) const;

    /// Check if this vertex and the other vertex are connected to the same
    /// vertex
    bool connectedToSameVertex(Vertex* other) const;

    /// Fetch the vertex connected (through an edge) to both this and
    /// the other vertex, if any
    Vertex* getCommonVertex(Vertex* other) const;

    /// Get the edge associated with two vertices, if any
    ftEdge* getCommonEdge(Vertex* other) const;

    /// Fetch the edge, if any, joining this vertex with the vertex other
    /// belonging to the specified face
    ftEdge* getCommonEdgeInFace(Vertex* other,
				ftSurface* face) const;

    /// Get the faces associated with two vertices, if any
    std::vector<ftSurface*> getCommonFaces(Vertex* other) const;

    /// Get the edges meeting in this vertex associated with a given face
    std::vector<ftEdge*> getFaceEdges(ftSurface *face) const;

    /// Check if the vertex is a corner in the given face
    bool isCornerInFace(ftSurface *face, double tol);

    /// Fetch the next vertices in the given face
    std::vector<shared_ptr<Vertex> > getNextVertex(ftSurface* face) const;

    /// Collect attached edges where the distance between the endpoints
    /// are larger than the specified tolerance or where the curves meet 
    /// with an angle that are more than the kink tolerance, but less than
    /// the corner tolerance
    void getEdgeDiscontinuities(std::vector<std::pair<ftEdge*, ftEdge*> >& gaps, double tol, 
				std::vector<std::pair<ftEdge*, ftEdge*> >& kinks, double angtol,
				double angmax) const;

     /// Reorganize edges according to twin information
    void reOrganize();

   /// Check the consistency of the edge information in this vertex
    // Debug functionality
    bool checkVertexTopology();

 private:
    /// The spacial position of the vertex
    Point vertex_point_;  

    /// Edges meeting in this vertex. Twin edges are collected in pairs.
    std::vector<std::pair<ftEdge*, ftEdge*> > edges_;
};

} // namespace Go


#endif // _VERTEX_H_
