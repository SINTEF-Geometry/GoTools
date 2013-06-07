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

#ifndef _TRIANG_H_
#define _TRIANG_H_



#include "GoTools/compositemodel/ttlPoint.h"
#include "ttl/ttl.h"
#include "ttl/ttl_util.h"
#include <memory>
#include "GoTools/compositemodel/ftPointSet.h"
#include <fstream>
#include <list>
#include <vector>

// The half-edge data structure
namespace hetriang {


  /// \b Node class in the half-edge data structure

  // Added pointer to specific object in class Node. This should be
  // handled by templates. (The user should define adequate node structure,
  // then pass it on to createDelaunay.)

  class Node {
    //static int id_count;
    double x_, y_, z_;

    Go::PointIter pnt_iter_;

  public:
    //int id_;
    Node(){init(0,0,0);}
    Node(double x, double y, double z = 0.0) {init(x,y,z);} 

  Node(Go::PointIter iter, double x, double y, double z = 0)
    : pnt_iter_(iter) {init(x,y,z);}

    ~Node(){}
    void init(double x, double y, double z) {x_ = x; y_ = y; z_ = z;} 
    
    double x() const {return x_;}
    double y() const {return y_;}
    double z() const {return z_;}
    const Go::PointIter& pointIter() const {return pnt_iter_;}
  };
  
  
  // ---------------------------------------------------------------------
  /// \b Edge class in the half-edge data structure
  class Edge {
    shared_ptr<Node> sourceNode_;
    Edge*        twinEdge_;
    Edge*        nextEdgeInFace_;
    bool         isLeadingEdge_;
    
  public:
  Edge() : sourceNode_(), twinEdge_(NULL), nextEdgeInFace_(NULL), isLeadingEdge_(false)
      {}
    ~Edge() {if(twinEdge_) twinEdge_->setTwinEdge(NULL);}
    
    void setSourceNode(shared_ptr<Node> node) {sourceNode_ = node;}
    void setNextEdgeInFace(Edge* edge) {nextEdgeInFace_ = edge;}
    void setTwinEdge(Edge* edge) {twinEdge_ = edge;}
    void setAsLeadingEdge(bool val=true){isLeadingEdge_ = val;}
    bool isLeadingEdge() const {return isLeadingEdge_;}
    Edge* getTwinEdge() const {return twinEdge_;};
    Edge* getNextEdgeInFace() const {return nextEdgeInFace_;}
    shared_ptr<Node> getSourceNode() {return sourceNode_;}
    shared_ptr<Node> getTargetNode() {return getNextEdgeInFace()->getSourceNode();}
  };

  // ---------------------------------------------------------------
  class Dart;       // Forward declaration (class in this namespace)

  /// \b Triangulation class for the half-edge data structure with adaption to TTL.
  class Triangulation {
  protected:  
    std::list<Edge*> leadingEdges_; // one half-edge for each arc
    void addLeadingEdge(Edge* edge) {edge->setAsLeadingEdge(); leadingEdges_.push_front(edge);}
    bool removeLeadingEdgeFromList(Edge* leadingEdge);
    void cleanAll();
    
  public:
    Triangulation() {}
    ~Triangulation() {cleanAll();}
    
    /// Create Delaunay triangulation from a set of points
    void createDelaunay(std::vector<Go::ttlPoint*>::iterator first,
                        std::vector<Go::ttlPoint*>::iterator last);

    /// When using rectangular boundary - loop through all points and expand.
    /// (Called from createDelaunay(...) when starting)
    Edge* initTwoEnclosingTriangles(std::vector<Go::ttlPoint*>::iterator first,
                                    std::vector<Go::ttlPoint*>::iterator last);

    // These two functions are required by TTL for Delaunay triangulation
    void swapEdge(Edge& diagonal);
    Edge* splitTriangle(Edge& edge, const Go::ttlPoint& point);

    // Functions required by TTL for removing nodes in a Delaunay triangulation
    void removeTriangle(Edge& edge); // boundary triangle required    
    void reverse_splitTriangle(Edge& edge);
    
    /// Create an arbitrary CCW dart
    Dart createDart();

    /// Get a list of "triangles" (one leading half-edge per triangle)
    const std::list<Edge*>& getLeadingEdges() const {return leadingEdges_;}
      int noTriangles() const {return (int)leadingEdges_.size();}

    
    /// One half-edge for each arc
    std::list<Edge*>* getEdges(bool skip_boundary_edges = false) const;

    /// Swap edges until the triangulation is Delaunay
    void optimizeDelaunay();

    /// Returns the start and stop-edge on a constrained edge with spesified nodes
    Edge* getConstrainedEdges(Go::ttlPoint point) const;

    /// Check if this triangulation is Delaunay
    bool checkDelaunay() const;    

    /// Get an arbitrary interior node (as the source node of the returned edge)
    Edge* getInteriorNode() const;
    
    /// Get an arbitrary boundary edge
    Edge* getBoundaryEdge() const;

    /// Print edges for plotting with, e.g., gnuplot
    void printEdges(std::ofstream& os) const;
  };
}; //end of half_edge_triang namespace

#endif
