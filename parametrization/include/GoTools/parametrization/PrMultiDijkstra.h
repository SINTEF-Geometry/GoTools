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

#ifndef _MultiDijkstra_H
#define _MultiDijkstra_H

#include <vector>
using std::vector;
#include <queue>

#define GraphType PrOrganizedPoints
#include "GoTools/parametrization/PrOrganizedPoints.h"

#define Point Vector3D


/** HeapNode2 -  Short description.
 * Detailed description.
 */
class HeapNode2 {
    
public:
  int     idx_;
  double  key_;
    
  HeapNode2() {}
  HeapNode2( int idx, double key) { idx_ = idx; key_ = key;}
  ~HeapNode2(){}
    
  bool operator >  (const HeapNode2& x) const {return (this->key_ >   x.key_);}
  bool operator >= (const HeapNode2& x) const {return (this->key_ >=  x.key_);}
};

typedef std::priority_queue< HeapNode2, vector<HeapNode2>, std::greater<HeapNode2> > HeapType2;

/** MultiDijkstra 
 * Class implementing the Dijkstra algorithm, where each node is labeled according to which
 * "candidate node" its path ends up in.
 */
class MultiDijkstra {

protected:

  GraphType*      graph_;
  HeapType2*      queue_;
  vector<double>* distances_;
  vector<int>*    back_trace_;
  vector<int>*    label_;
  vector<int>*    flags_;
  double          large_distance_;

  void    setFinished(int node_idx);
  int     isFinished(int node_idx);
  void    insertInCandidates(int node_idx); 

public:
  /// Constructor
  MultiDijkstra();
  /// Destructor
  ~MultiDijkstra();
  
  /// Set the graph to run the algorithm on.  Must be run before 'initialize()'.
  void    setGraph(GraphType*  tri) {graph_ = tri;}
  /// Initialize internal data structures.  Must be done after 'setGraph()' has
  /// been called.
  void    initialize();
  /// Set a 'source' for the algorithm.  This is a node for which the associated
  /// minimum distance is known.  At least one source must be set for the algorithm
  /// to work properly.
  /// \param node_idx the index of the node to be designated as a source
  /// \param node_label set a label for this source (makes it possibly to identify
  ///                   which nodes in the graph whose solution depends on this node.
  /// \param dist the distance at this node (typically 0).
  void    setSource(int node_idx, int node_label, double dist = 0);

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
  /// After the algorithm has been run (using the run() function _without_
  /// arguments), this function can be used to find which label each node has been
  /// assigned.
  int getLabel(int node_idx) {return (*label_)[node_idx];}

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
inline void MultiDijkstra::setFinished( int node_idx) 
//-----------------------------------------------------------------------------
{
  (*flags_)[node_idx] = 1;
}

//-----------------------------------------------------------------------------
inline int MultiDijkstra::isFinished( int node_idx) 
//-----------------------------------------------------------------------------
{
  return (*flags_)[node_idx];
}


//-----------------------------------------------------------------------------
inline void MultiDijkstra::insertInCandidates( int node_idx) 
//-----------------------------------------------------------------------------
{
  queue_->push(HeapNode2(node_idx,getDistance(node_idx)));
}


//-----------------------------------------------------------------------------
inline double MultiDijkstra::getDistance( int n1, int n2) 
//-----------------------------------------------------------------------------
{
  // assume n1 & n2 are neighbours in graph_
  Point  p1 = graph_->get3dNode(n1);
  Point  p2 = graph_->get3dNode(n2);
  return p1.dist(p2);
}


#endif 
