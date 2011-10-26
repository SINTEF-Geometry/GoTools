/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 2000 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrLevelTriangulation_OP.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Nov. 2000
 DESCRIPTION : Implementation of methods in the class PrLevelTriangulation_OP.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrLevelTriangulation_OP.h"

// PUBLIC MEMBER FUNCTIONS

//----------------------------------------------------------------------------
PrLevelTriangulation_OP::
PrLevelTriangulation_OP(vector<PrNestedNode>* node, vector<PrTriangle>& nt, int level) 
    : node_(node), triangle_(nt), level_(level), numNodes_((int)node->size())
//-----------------------------------------------------------------------------
{
//   node_ = node;
//   //triangle_.resize(nt.size());
//   //int j;
//   //for(j=0; j<t.findNumFaces(); j++)
//   //{
//     //triangle_[j] = nt[j];
//   //}
//   triangle_ = nt; //copy constructor for STL vectors

//   level_ = level;
//   numNodes_ = node->size();
}

//----------------------------------------------------------------------------
void
PrLevelTriangulation_OP::getNeighbours(int k, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the neighbours of the k-th node in:
//   1. any anticlockwise order if the k-th node is an interior node
//   2. the unique anticlockwise order if the k-th node is a boundary node.
{
  int tr1 = (*node_)[k].tr(level_);
  neighbours.clear();
  if (tr1>=0)
  {
     neighbours.push_back(triangle_[tr1].getAnticlockwiseNode(k));
     int tr,trNext;
     for(tr = tr1, trNext = triangle_[tr].getLeftTriangle(k);
	 trNext > -1 && trNext != tr1;
	 tr = trNext, trNext = triangle_[tr].getLeftTriangle(k))
     {
	neighbours.push_back(triangle_[tr].getClockwiseNode(k));
     }
     if(trNext == -1) neighbours.push_back(triangle_[tr].getClockwiseNode(k));
     // k is a boundary node
  }
}

//----------------------------------------------------------------------------
bool PrLevelTriangulation_OP::isBoundary(int k) const
//-----------------------------------------------------------------------------
//   Given a node and its leading edge,
//   return 1 if the node is a boundary node and 0 otherwise.
{
  int tr = (*node_)[k].tr(level_);
  if(triangle_[tr].getRightTriangle(k) == -1) return true;
  else return false;
}

//-----------------------------------------------------------------------------
void PrLevelTriangulation_OP::printXYZTriangles(ostream& os, bool num)
//-----------------------------------------------------------------------------
// Print out the triangles of the graph, useful for plotting
// If num = 1, print first the number of triangles.
{
  if(num) os << triangle_.size() << '\n';
  for(size_t i=0; i<triangle_.size(); i++)
  {
    (*node_)[triangle_[i].n1()].printXYZ(os);
    (*node_)[triangle_[i].n2()].printXYZ(os);
    (*node_)[triangle_[i].n3()].printXYZ(os);
    os << "\n";
  }
}

//-----------------------------------------------------------------------------
void PrLevelTriangulation_OP::printUVTriangles(ostream& os, bool num)
//-----------------------------------------------------------------------------
// Print out the uv triangles of the graph, useful for plotting
// If num = 1, print first the number of triangles.
{
  if(num) os << triangle_.size() << '\n';
  for(size_t i=0; i<triangle_.size(); i++)
  {
    (*node_)[triangle_[i].n1()].printUV(os);
    (*node_)[triangle_[i].n2()].printUV(os);
    (*node_)[triangle_[i].n3()].printUV(os);
    os << "\n";
  }
}

//-----------------------------------------------------------------------------
void PrLevelTriangulation_OP::print(ostream& os)
//-----------------------------------------------------------------------------
{
  os << numNodes_ << ' ' <<triangle_.size() << '\n';
  int i;
  for(i=0; i<numNodes_; i++) (*node_)[i].print(os);
  os << "\n";
  int numTri = (int)triangle_.size();
  for(i=0; i<numTri; i++) triangle_[i].print(os);
}

//-----------------------------------------------------------------------------
void PrLevelTriangulation_OP::printRawData(ostream& os)
//-----------------------------------------------------------------------------
{
    os << numNodes_ << ' ' << triangle_.size() << '\n';
  int i;
  for(i=0; i<numNodes_; i++)
    os << (*node_)[i].x() << ' ' << (*node_)[i].y()
       << ' ' << (*node_)[i].z() << '\n';
  os << "\n";
  int numTri = (int)triangle_.size();
  for(i=0; i<numTri; i++)
    os << triangle_[i].n1() << ' '
       << triangle_[i].n2() << ' '
       << triangle_[i].n3() << '\n';
}

//-----------------------------------------------------------------------------
void PrLevelTriangulation_OP::printUV(ostream& os)
//-----------------------------------------------------------------------------
{
  os << numNodes_ << ' ' << triangle_.size() << '\n';
  int i;
  for(i=0; i<numNodes_; i++)
    os << (*node_)[i].u() << ' ' << (*node_)[i].v() << " 0.0\n";
  os << "\n";
  int numTri = (int)triangle_.size();
  for(i=0; i<numTri; i++)
    os << triangle_[i].n1() << ' '
       << triangle_[i].n2() << ' '
       << triangle_[i].n3() << std::endl;
}

