#include "GoTools/parametrization/PrDijkstra.h"
#include <assert.h>
//#include <limits>

//-----------------------------------------------------------------------------
Dijkstra::Dijkstra() 
//-----------------------------------------------------------------------------
{
  large_distance_ = 1e100; //(double)std::numeric_limits<double>::infinity();
  queue_      = new HeapType;
  distances_  = new vector<double>;
  back_trace_ = new vector<int>;
  flags_      = new vector<int>;
  graph_      = NULL;
}


//-----------------------------------------------------------------------------
Dijkstra::~Dijkstra() 
//-----------------------------------------------------------------------------
{
  if(queue_)      delete queue_;
  if(distances_)  delete distances_;
  if(flags_)      delete flags_;
  if(back_trace_) delete back_trace_;
}



//-----------------------------------------------------------------------------
void Dijkstra::initialize() 
//-----------------------------------------------------------------------------
{
  int i=0;

  if (graph_ == NULL) {
    assert(graph_ != NULL);
    return;
  }

  int size = graph_->getNumNodes();
  // queue_->reserve(size);
  // queue_->empty();
  if (queue_) delete queue_;
  queue_ = new HeapType;
  
  distances_->reserve(size);
  back_trace_->reserve(size);
  flags_->reserve(size);

  
  for (i=0; i<size; i++) {
    (*distances_)[i]  = large_distance_;
    (*flags_)[i]      = 0;
    (*back_trace_)[i] = -1;
  }

}



//-----------------------------------------------------------------------------
void Dijkstra::setSource( int node_idx, double dist)
//-----------------------------------------------------------------------------
{
  (*distances_)[node_idx] = dist;
  // setFinished(node_idx);
  insertInCandidates(node_idx);
}




//-----------------------------------------------------------------------------
void Dijkstra::run()
//-----------------------------------------------------------------------------
{
  vector<int>  neighbours;
  double dist,currdist,nextdist;
  int curr,next;
  int i;
  
  
  while (!queue_->empty()) {
    curr     = queue_->top().idx_; 
    queue_->pop();
    while(isFinished(curr) && !queue_->empty()) {
      curr     = queue_->top().idx_; 
      queue_->pop();
    }

    // Return if node is finished, i.e. queue is empty 
    if (isFinished(curr)) {
      return ;
    }

    currdist = (*distances_)[curr];


    setFinished(curr);
    graph_->getNeighbours(curr,neighbours);
    
    int s = (int)neighbours.size();
    for ( i=0; i<s; i++) {
      next = neighbours[i];
      if (!isFinished(next)){
        dist = getDistance(curr,next);
        nextdist = currdist+dist;
        if (nextdist < (*distances_)[next]) {
          (*distances_)[next] = nextdist;
          insertInCandidates(next);
        } // if
      } // if
    } // for
  } // while
}


//-----------------------------------------------------------------------------
void Dijkstra::run( double radius)
//-----------------------------------------------------------------------------
{
  vector<int>  neighbours;
  double dist,currdist,nextdist;
  int curr,next;
  int i;
  
  
  while (!queue_->empty()) {
    curr     = queue_->top().idx_; 
    queue_->pop();
    while(isFinished(curr) && !queue_->empty()) {
      curr     = queue_->top().idx_; 
      queue_->pop();
    }

    // Return if node is finished, i.e. queue is empty 
    if (isFinished(curr)) {
      return ;
    }

    currdist = (*distances_)[curr];
    setFinished(curr);
    if (currdist>=radius)
      continue;
    
    graph_->getNeighbours(curr,neighbours);
    
    int s = (int)neighbours.size();
    for ( i=0; i<s; i++) {
      next = neighbours[i];
      if (!isFinished(next)){
        dist = getDistance(curr,next);
        nextdist = currdist+dist;
        if (nextdist < (*distances_)[next]) {
          (*distances_)[next] = nextdist;
          insertInCandidates(next);
        } // if
      } // if
    } // for
  } // while
}

//-----------------------------------------------------------------------------
void    Dijkstra::run( int  target)
//-----------------------------------------------------------------------------
{
  vector<int>  neighbours;
  double dist,currdist,nextdist;
  int curr,next;
  int i;

  
  
  while (!queue_->empty()) {
    curr     = queue_->top().idx_; 
    queue_->pop();
    while(isFinished(curr) && !queue_->empty()) {
      curr     = queue_->top().idx_; 
      queue_->pop();
    }

    // Return if node is finished, i.e. queue is empty 
    if (isFinished(curr)) {
      return ;
    }

    currdist = (*distances_)[curr];
    setFinished(curr);

    // Finished
    if (curr == target)
      return;
    
    graph_->getNeighbours(curr,neighbours);
    
    int s = (int)neighbours.size();
    for ( i=0; i<s; i++) {
      next = neighbours[i];
      if (!isFinished(next)){
        dist = getDistance(curr,next);
        nextdist = currdist+dist;
        if (nextdist < (*distances_)[next]) {
          (*distances_)[next] = nextdist;
	  (*back_trace_)[next] = curr;
          insertInCandidates(next);
        } // if
      } // if
    } // for
  } // while
}


//-----------------------------------------------------------------------------
int Dijkstra::closestNeighbour(int node_idx)
//-----------------------------------------------------------------------------
{
  if ((*back_trace_)[node_idx] == -1) {
    std::cerr << "!!! PrDijkstra: backtrace didn't work !!!" << std::endl;
    return (0);
  }

  return (*back_trace_)[node_idx];

  // bad code !!!
/*
  vector<int> neighbours;
  graph_->getNeighbours(node_idx,neighbours);

  if (neighbours.size()==0)
    return -1;

  double closest_dist = getDistance(neighbours[0]);
  int    closest_node = neighbours[0];
  for (int i=1; i<neighbours.size();i++) {
    int    node = neighbours[i];
    double dist = getDistance(node);
    
    if (fabs (dist-closest_dist) < 0.0000000001)
      std::cerr << "!!! PrDijkstr: two paths with same length !!!" << std::endl;

    if (dist<closest_dist) {
      closest_node = node;
      closest_dist = dist;
    }
  }
  return closest_node;
*/
}




