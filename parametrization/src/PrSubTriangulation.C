/*****************************************************************************/
/*                                                                           */
/* (c) Copyright 2000 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrSubTriangulation.C
 AUTHOR      : Kai Hormann, SINTEF
 DATE        : Oct. 2000
 DESCRIPTION : Implementation of methods in the class PrSubTriangulation.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrSubTriangulation.h"
#include <set>
using std::set;
using std::cerr;
using std::cout;
using std::endl;

// PUBLIC MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
PrSubTriangulation::PrSubTriangulation(PrTriangulation_OP* graph, 
				       vector<int> boundary,
				       vector<int> interior)
//-----------------------------------------------------------------------------
//   Construct the PrSubTriangulation from an array with the indices
//   of the "boundary" nodes and another carrying the indices of the
//   "interior" nodes. All indices are relative to the indexing in
//   the associated PrTriangulation_OP "graph"
{
  // make the data of the PrTriangulation_OP available
  attach(graph);

  // initialize the index mappings
  initialize(boundary, interior);
}


//----------------------------------------------------------------------------
void
PrSubTriangulation::initialize(vector<int> boundary, vector<int> interior)
//-----------------------------------------------------------------------------
{
  size_t i;
  // set size of boundary and interior nodes
  numBdy_ = (int)boundary.size();
  numInt_ = (int)interior.size();

  // build the mapping "new2old"
  new2old_ = boundary;
  for (i=0; i<interior.size(); i++)
    new2old_.push_back(interior[i]);

  // build the mapping "old2new", needed esp. by getNeighbours
  old2new_.clear();
  for (i=0; i<new2old_.size(); i++)
      old2new_[new2old_[i]] = (int)i;

  cerr << " * A SubTriangulation with " << numBdy_ << " boundary and ";
  cerr << numInt_ << " interior nodes has been created!" << endl;
/* DEBUG
  cerr << "   new2old: ";
  for (i=0; i<new2old_.size(); i++)
    cerr << i << "->" << new2old_[i] << ", ";
  cerr << "\n   old2new: ";
  typedef map<int, int>::const_iterator MI;
  for(MI m=old2new_.begin(); m!=old2new_.end(); ++m)
    cerr << m->first << "->" << m->second << ", ";
  cerr << endl;
*/
}


//----------------------------------------------------------------------------
void
PrSubTriangulation::getNeighbours(int k, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the neighbours of the k-th node in:
//   1. any anticlockwise order if the k-th node is an interior node
//   2. the unique anticlockwise order if the k-th node is a boundary node.
{
  // get the "old" indices of the neighbours from the original triangulation
  vector<int> oldNeighbours;
  g_->getNeighbours(new2old_[k], oldNeighbours);

  // build the array of "new" neighbour indices
  neighbours.clear();
  if (isBoundary(k)) {
    // since we assume the boundary vertices to be ordered anticlockwise,
    // it's easy to find the neighbouring boundary vertices:
    int next = (k+1)%numBdy_;
    int last = (k+numBdy_-1)%numBdy_;

    // the first neighbour shall be the next (ccw) boundary node
    neighbours.push_back(next);

    // search this node in the "old" neighbours
    int i=0;
    while (oldNeighbours[i] != new2old_[next])
      i++;

    // in case of a boundary node with the same two neighbours
    // (degeneracy!) add that neighbour twice (the while-loop
    // will be skipped)
    if (new2old_[next] == new2old_[last])
      neighbours.push_back(last);
      
    // transfer the "old" indices in anticlockwise order 
    // up to the "last" node
    while (new2old_[next] != new2old_[last]) {
      i++;
      if (old2new_.find(oldNeighbours[i%oldNeighbours.size()])
	  == old2new_.end())
	std::cerr << "OOPS anticlockwise ordering possible????\n";
      else
      {
	next = (*(old2new_.find(oldNeighbours[i%oldNeighbours.size()]
				))).second;
      neighbours.push_back(next);
    }
  }
  }
  else {
    // for interior nodes, copy all neighbour indices
    for (size_t i=0; i<oldNeighbours.size(); i++)
      neighbours.push_back((*(old2new_.find(oldNeighbours[i]))).second);
  }
}

//----------------------------------------------------------------------------
void
PrSubTriangulation::getTriangleIndices(vector<int>& indices) const
//-----------------------------------------------------------------------------
{
  set<int> tri_indices;
  vector<int> nbrs(20);

  // find the triangles around interior vertices
  for (int i=numBdy_; i<numBdy_+numInt_; i++) {
    getNeighbourTriangles(i, nbrs);
    for (size_t j=0; j<nbrs.size(); j++) {
      tri_indices.insert(nbrs[j]);
    }
  }

  // find the triangles around boundary vertices 
  // REMARK: take only those with all three vertices in this SubTriangulation!
  for (int i=0; i<numBdy_; i++) {
    getNeighbourTriangles(i, nbrs);
    for (size_t j=0; j<nbrs.size(); j++) {
      if (triangleInThis(nbrs[j]))
	tri_indices.insert(nbrs[j]);
    }
  }

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
  indices.clear();
  std::copy(tri_indices.begin(), tri_indices.end(), indices.begin());
#else
  indices.assign (tri_indices.begin(), tri_indices.end());
#endif

}

//----------------------------------------------------------------------------
void
PrSubTriangulation::getNeighbourTriangles(int i, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
{
  int trEnd  = g_->getPrNode(new2old_[i]).tr();
  int trNext = g_->getPrTriangle(trEnd).getLeftTriangle(new2old_[i]);

  neighbours.clear();
  neighbours.push_back(trEnd);

  while ((trNext != trEnd) && (trNext != -1)) {
    neighbours.push_back(trNext);    
    trNext = g_->getPrTriangle(trNext).getLeftTriangle(new2old_[i]);
  }
}

//----------------------------------------------------------------------------
bool
PrSubTriangulation::triangleInThis(int i) const
//-----------------------------------------------------------------------------
{
  PrTriangle tri = g_->getPrTriangle(i);

  if (old2new_.find(tri.n1()) == old2new_.end())
    return false;
  if (old2new_.find(tri.n2()) == old2new_.end())
    return false;
  if (old2new_.find(tri.n3()) == old2new_.end())
    return false;

  return true;
}
