#ifndef _Dijkstra_H
#define _Dijkstra_H

#include <vector>
#include <queue>

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
