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

#include "GoTools/parametrization/PrThin.h"
using std::cerr;
using std::cout;
using std::endl;

// abbreviations:

#define Node trn_->getPrNode
#define Triangle trn_->getPrTriangle

// PRIVATE METHODS

//-----------------------------------------------------------------------------
void PrThin::makeHeap()
//-----------------------------------------------------------------------------
{
  //  cout << "Constructing the heap..." << endl;
  //double time,time2;
  //CPUclock rolex;
  //time = rolex.getTime();

  //heap_.redim(trn_->getNumNodes());
 
  heap_ = new PrHeap(trn_->getNumNodes());

  for(int i=0; i<trn_->getNumNodes(); i++)
    heap_->push(findError(i), i+1);

  //time2 = rolex.getTime() - time;
  //cout << "Heap constructed in " << time2 << " units CPU time" << endl;
}

// PUBLIC METHODS

//-----------------------------------------------------------------------------
PrThin::~PrThin()
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void PrThin::attach(PrTriangulation_OP* trn)
//-----------------------------------------------------------------------------
{
  trn_ = trn;

  makeHeap();
}

//-----------------------------------------------------------------------------
double PrThin::findError(int i)
//-----------------------------------------------------------------------------
{
  int tr;

  if(trn_->isBoundary(i))
    return 1.0e30;

  findNearestNeighbour(i,tr);
  int j = Triangle(tr).getAnticlockwiseNode(i);
  Vector3D diff = Node(j).point() - Node(i).point();
  return diff.length();
}

//-----------------------------------------------------------------------------
void PrThin::findNearestNeighbour(int i, int& tr)
//-----------------------------------------------------------------------------
// Find nearest neighbour j to an interior node i (in 3D) and
// return the triangle tr which is left of the directed edge (i,j).
// Note that j is then tr.getAnticlockwiseNode(i).
{
  if(trn_->isBoundary(i))
  {
    tr = 0;
    return;
  }

  double dist,mindist;
  Vector3D diff;

  int tr1 = Node(i).tr();
  int j = Triangle(tr1).getAnticlockwiseNode(i);
  diff = Node(j).point() - Node(i).point();
  mindist = diff.length();
  tr = tr1;

  int t = Triangle(tr1).getLeftTriangle(i);
  while(t != tr1)
  {
    j = Triangle(t).getAnticlockwiseNode(i);
    diff = Node(j).point() - Node(i).point();
    dist = diff.length();
    if(dist < mindist)
    {
      mindist = dist;
      tr = t;
    }
    t = Triangle(t).getLeftTriangle(i);
  }
}

//-----------------------------------------------------------------------------
void PrThin::halfEdgeCollapse(int i)
//-----------------------------------------------------------------------------
// Remove node i. 
{
  if(trn_->isBoundary(i)) return;

  // First find nearest neighbour j.
  int u0;
  findNearestNeighbour(i,u0);

  int j = Triangle(u0).getAnticlockwiseNode(i);

  // Collapse cell to node j.

  // Set triangle pointer for i to -1
  Node(i).tr() = -1;

  // Find triangle to left of directed edge (i,j).

  int s0 = Triangle(u0).getOppositeTriangle(i);
  int t0 = Triangle(u0).getOppositeTriangle(j);

  if (s0 >= 0)
    Triangle(s0).replaceTriangle(u0,t0);
  Triangle(t0).replaceTriangle(u0,s0);

  int u1 = Triangle(u0).getRightTriangle(i);
  int s1 = Triangle(u1).getOppositeTriangle(i);
  int t1 = Triangle(u1).getOppositeTriangle(j);

  if (s1 >= 0)
    Triangle(s1).replaceTriangle(u1,t1);
  Triangle(t1).replaceTriangle(u1,s1);

  // updating the vertex->triangle pointers

  if(Node(j).tr() == u0) Node(j).tr() = t0;
  if(Node(j).tr() == u1) Node(j).tr() = t1;

  int k1 = Triangle(u1).getAnticlockwiseNode(i);
  if(Node(k1).tr() == u1) Node(k1).tr() = t1;

  int k0 = Triangle(u0).getAnticlockwiseNode(j);
  if(Node(k0).tr() == u0) Node(k0).tr() = t0;

  //for triangles t0 to t1 anticlockwise around i,
  //replace the pointer to node i by a pointer to node j.

  Triangle(t0).replaceNode(i,j);

  int t = t0;
  while(t != t1)
  {
    t = Triangle(t).getLeftTriangle(j);
    Triangle(t).replaceNode(i,j);
  }

  Triangle(u0).n1() = -1;
  Triangle(u1).n1() = -1;
}

//-----------------------------------------------------------------------------
void PrThin::removePoint(int i)
//-----------------------------------------------------------------------------
// Remove the node i from the triangulation and update heap.
// It is assumed that i is not in the heap.
{
  //  cout << "removing point " << i << endl;

  // Store the neighbours. Their priorities will be updated
  // after the cell has been retriangulated.

  vector<int> neighbours;
  trn_->getNeighbours(i,neighbours);

  // Delete the node from the triangulation.

  halfEdgeCollapse(i);

  // Update the priorities of the neighbours.

  int k;
  double error;
  for(size_t j=0; j< neighbours.size(); j++)
  {
    k = neighbours[j];
    if(trn_->isBoundary(k)) continue;
    error = findError(k);
    heap_->modify(error,k+1);
  }
}

//-----------------------------------------------------------------------------
void PrThin::thin()
//-----------------------------------------------------------------------------
// Thin the points, until the error of further removals is above the
// threshold error_;
{
  //  cout << "Thinning started..." << endl;
  //CPUclock rolex;
  //double time = rolex.getTime();

  // Pop the top node (with least error) from the heap.
  static double error=0.0;
  int i;
  int nos = steps_;

  heap_->pop(error,i); i--;

  for (int z=0; z<nos; z++)
    if (error <= error_)
    {
      int tr;
      findNearestNeighbour(i,tr);
      int j = Triangle(tr).getAnticlockwiseNode(i);

      vector<int> i_nbrs;
      vector<int> j_nbrs;

      trn_->getNeighbours(i, i_nbrs);
      trn_->getNeighbours(j, j_nbrs);

      int cnt = 0;

      for (size_t k=0; k<i_nbrs.size(); k++)
	for (size_t l=0; l<j_nbrs.size(); l++)
          if (i_nbrs[k] == j_nbrs[l])
            cnt++;

      if ((cnt == 2) && (!trn_->isBoundary(j)))
        removePoint(i);
      else
      {
	//	cout << "collapsing " << i << " to " << j << " causes problems ";
	//cout << "because they have " << cnt << " common neighbours!" << endl;
        heap_->push(1.0e30, i+1);
      }

      if (z < nos-1)
        heap_->pop(error,i); i--;
    }

//  while(error <= error_)
//  {
//    removePoint(i);
//    heap_->pop(error,i); i--;
//  }

  //double time2 = rolex.getTime() - time;
  //cout << "Thinning done after " << time2 << " units CPU time" << endl;
}



//-----------------------------------------------------------------------------
void PrThin::printHeap()
//-----------------------------------------------------------------------------
{
  //cout << "max. allowed error:" << error_ << endl;
  heap_->print(cout);
}

