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

#ifndef _Dijkstra_H
#define _Dijkstra_H

#include <vector>
#include <queue>
#include <functional>

#define GraphType PrOrganizedPoints
#include "GoTools/parametrization/PrOrganizedPoints.h"

//#define Point Vector3D
 
using std::vector;

class HeapNode {
    
public:
  int     idx_;
  double  key_;

  /// Constructor    
  HeapNode() {}
  /// Constructor
  HeapNode( int idx, double key) { idx_ = idx; key_ = key;}
  /// Destructor
  ~HeapNode(){}
    
  bool operator >  (const HeapNode& x) const {return (this->key_ >   x.key_);}
  bool operator >= (const HeapNode& x) const {return (this->key_ >=  x.key_);}
};

typedef std::priority_queue< HeapNode, vector<HeapNode>, std::greater<HeapNode> > HeapType;


/** Dijkstra 
 * Class implementing the Dijkstra algorithm
 */
class Dijkstra {

protected:

  const GraphType*      graph_;
  HeapType*       queue_;
  vector<double>* distances_;
  vector<int>*    back_trace_;
  vector<int>*    flags_;
  double          large_distance_;

  void    setFinished(int node_idx);
  int     isFinished(int node_idx);
  void    insertInCandidates(int node_idx); 

public:
  /// Constructor
  Dijkstra();
  /// Destructor
  ~Dijkstra();
  
  /// Set the graph to run the algorithm on.  Must be run before 'initialize()'.
  void    setGraph(const GraphType*  tri) {graph_ = tri;}
  /// Initialize internal data structures.  Must be done after 'setGraph()' has
  /// been called.
  void    initialize();
  /// Set a 'source' for the algorithm.  This is a node for which the associated
  /// minimum distance is known.  At least one source must be set for the algorithm
  /// to work properly.
  /// \param node_idx the index of the node to be designated as a source
  /// \param dist the distance at this node (typically 0).
  void    setSource(int node_idx, double dist = 0);
  /// Run the algorithm.  Compute the associated distance for all nodes in the
  /// graph.
  void    run();
  /// Run the algorithm.  Compute the associated distance for all nodes in the
  /// graph.  Edges longer than 'radius' will be ignored.
  void    run(double radius);
  /// Run the algorithm.  Compute the associated distance for nodes in the graph
  /// until the distance for the 'target' node has been found.  Then skip further
  /// computations.
  void    run(int target);
  /// Get the distance associated with the node 'node_ix'.  Does not give a valid
  /// answer until the Dijkstra algorithm has been run().
  double  getDistance(int node_idx) {return (*distances_)[node_idx];}
  /// Compute the euclidean distance between two nodes, supposedly neighbours.
  /// If they are boundary nodes, their mutual distance will be set to a very
  /// large number.
  double  getDistance(int node_idx, int neighbour);

  /// Get the index of the closest neighbour to the node 'node_idx'.  Does not give
  /// a valid answer until the Dijkstra algorithm has been run().
  int closestNeighbour(int node_idx);
};


//-----------------------------------------------------------------------------
inline void Dijkstra::setFinished( int node_idx) 
//-----------------------------------------------------------------------------
{
  (*flags_)[node_idx] = 1;
}

//-----------------------------------------------------------------------------
inline int Dijkstra::isFinished( int node_idx) 
//-----------------------------------------------------------------------------
{
  return (*flags_)[node_idx];
}


//-----------------------------------------------------------------------------
inline void Dijkstra::insertInCandidates( int node_idx) 
//-----------------------------------------------------------------------------
{
  queue_->push(HeapNode(node_idx,getDistance(node_idx)));
}


//-----------------------------------------------------------------------------
inline double Dijkstra::getDistance( int n1, int n2) 
//-----------------------------------------------------------------------------
{
  if (graph_->isBoundary(n1) && graph_->isBoundary(n2))
    return 1e200;
  // assume n1 & n2 are neighbours in graph_
  Vector3D  p1 = graph_->get3dNode(n1);
  Vector3D  p2 = graph_->get3dNode(n2);
  return p1.dist(p2);
}


#endif 
