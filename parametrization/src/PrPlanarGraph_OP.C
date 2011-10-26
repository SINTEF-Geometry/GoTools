/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1997 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/


/********************************************************************
 FILENAME    : PrPlanarGraph_OP.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Aug. 97
 DESCRIPTION : Implementation of methods in the class PrPlanarGraph_OP.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrPlanarGraph_OP.h"

// PUBLIC MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
PrPlanarGraph_OP::PrPlanarGraph_OP(int npts, double* xyz_nodes, double* uv_nodes, 
                                   int* end, int* adj)
//-----------------------------------------------------------------------------
//   Construct the PrPlanarGraph_OP from an array of nodes and
//   an array of adjacency lists.
//   The length of the arrays xyz_nodes and uv_nodes should be 3*npts
//   and 2*npts respectively.
//   The length of the array end should be npts.
//   The length of the array adj should be end[npts].
//   If node 1 is an interior node, its neighbours are
//   adj[1],...,adj[end[1]] in any anticlockwise sequence.
//   If node 1 is a boundary node, its neighbours are
//   adj[1],...,adj[end[1]-1] in the unique anticlockwise sequence.
//   and adj[end[1]] = 0 to indicate that node 1 is on the boundary.
//   For i=2,...,n:
//   If node i is an interior node, its neighbours are
//   adj[end[i-1]],...,adj[end[i]] in any anticlockwise sequence.
//   If node i is a boundary node, its neighbours are
//   adj[end[i-1]],...,adj[end[i]-1] in the unique anticlockwise sequence.
//   and adj[end[i]] = 0 to indicate that node i is on the boundary.
//   This is the data structure proposed by Cline and Renka.
    : node_(npts)
{
  int j;
  for(j=0; j< npts; j++)
  {
    node_[j].x() = xyz_nodes[3*j];
    node_[j].y() = xyz_nodes[3*j+1];
    node_[j].z() = xyz_nodes[3*j+2];
    node_[j].u() = uv_nodes[2*j];
    node_[j].v() = uv_nodes[2*j+1];
    node_[j].end() = end[j];
  }

  adj_.resize(node_[npts-1].end() + 1);
  for(j=0; j< (node_[npts-1].end() + 1); j++) adj_[j] = adj[j];
}

//-----------------------------------------------------------------------------
PrPlanarGraph_OP::PrPlanarGraph_OP(int npts, double* xyz_nodes,
                                   int* end, int* adj)
//-----------------------------------------------------------------------------
//   Construct the PrPlanarGraph_OP from an array of nodes and
//   an array of adjacency lists. Like the previous constructor,
//   only we set the uv's to zero.
    : node_(npts)
{
  int j;
  for(j=0; j< npts; j++)
  {
    node_[j].x() = xyz_nodes[3*j];
    node_[j].y() = xyz_nodes[3*j+1];
    node_[j].z() = xyz_nodes[3*j+2];
    node_[j].u() = 0.0;
    node_[j].v() = 0.0;
    node_[j].end() = end[j];
  }

  adj_.resize(node_[npts-1].end() + 1);
  for(j=0; j< (node_[npts-1].end() + 1); j++) adj_[j] = adj[j];
}

//-----------------------------------------------------------------------------
PrPlanarGraph_OP::PrPlanarGraph_OP(PrOrganizedPoints& op)
//-----------------------------------------------------------------------------
//   Construct the PrPlanarGraph_OP from a PrOrganizedPoints class
{
  int npts = op.getNumNodes();
  node_.resize(npts);
  Vector3D nd;
  int j, end = -1;
  vector<int> neighbours;
  for(j=0; j< npts; j++)
  {
    nd = op.get3dNode(j); 
    node_[j].x() = nd.x();
    node_[j].y() = nd.y();
    node_[j].z() = nd.z();
    node_[j].u() = op.getU(j);
    node_[j].v() = op.getV(j);
    op.getNeighbours(j,neighbours);
    end += (int)neighbours.size();
    if(op.isBoundary(j)) end++;
    node_[j].end() = end;
  }

  adj_.resize(node_[npts-1].end() + 1);
  int end_prev;
  int deg;
  for(j=0; j< npts; j++)
  {
    end_prev = (j == 1 ? -1 : node_[j-1].end() );
    op.getNeighbours(j,neighbours);
    deg = (int)neighbours.size();
    for(int i=0; i<deg; i++)
    {
      adj_[end_prev+i] = neighbours[i];
    }
    if(op.isBoundary(j)) adj_[end_prev+deg+1] = 0;
  }
}

//----------------------------------------------------------------------------
void
PrPlanarGraph_OP::getNeighbours(int k, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the neighbours of the k-th node in:
//   1. any anticlockwise order if the k-th node is an interior node
//   2. the unique anticlockwise order if the k-th node is a boundary node.
{
  neighbours.clear();
  int start = (k == 0 ? 0 : node_[k-1].end() + 1);
  int end   = (isBoundary(k) ? node_[k].end() - 1 : node_[k].end());
  for(int i=start; i<= end; i++) neighbours.push_back(adj_[i]);
}

//----------------------------------------------------------------------------
bool PrPlanarGraph_OP::isBoundary(int k) const
//-----------------------------------------------------------------------------
//   Given a node
//   return 1 if the node is a boundary node and 0 otherwise.
{
  if(adj_[node_[k].end()] == -1) return true;
  else return false;
}

//-----------------------------------------------------------------------------
void PrPlanarGraph_OP::print(std::ostream& os)
//-----------------------------------------------------------------------------
{
    os << node_.size() << '\n';
    for(size_t i=0; i<node_.size(); i++) node_[i].print(os);
    for(size_t i=0; i<node_.size(); i++)
	{
	    int start = (i == 0 ? 0 : node_[i-1].end() + 1);
	    for(int j=start; j<= node_[i].end(); j++) os << adj_[j] << '\n';
	    os << std::endl;
	}
}

//-----------------------------------------------------------------------------
void PrPlanarGraph_OP::scan(std::istream& is)
//-----------------------------------------------------------------------------
{
  int numpnts;
  is >> numpnts;
  node_.resize(numpnts);
  int i;
  for(i=0; i<numpnts; i++) node_[i].scan(is);

  adj_.resize(node_[numpnts-1].end() + 1);
  int as = (int)adj_.size();
  for(i=0; i<as; i++) is >> adj_[i];
}

//-----------------------------------------------------------------------------
void PrPlanarGraph_OP::scan2(std::istream& is)
//-----------------------------------------------------------------------------
{
  int numpnts = 0;
  is >> numpnts;
  node_.resize(numpnts);
  int i;
  for(i=0; i<numpnts; i++)
  {
    is >> node_[i].x() >> node_[i].y() >> node_[i].z();
    is >> node_[i].u() >> node_[i].v(); 
  }
  for(i=0; i<numpnts; i++)
  {
    is >> node_[i].end(); 
  }

  adj_.resize(node_[numpnts-1].end() + 1);
  int as = (int)adj_.size();
  for(i=0; i<as; i++) is >> adj_[i];
}

