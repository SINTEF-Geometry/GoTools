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

