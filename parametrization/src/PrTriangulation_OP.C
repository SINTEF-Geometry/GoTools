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

#include "GoTools/parametrization/PrTriangulation_OP.h"

// PRIVATE MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
int PrTriangulation_OP::getNghrTriangle(int n1, int n2,
                        vector<int>& tlist)
//-----------------------------------------------------------------------------
{
// Given a list (tlist) of triangles, all of which contain node n1,
// find (if it exists) the triangle which lies
// to the right of edge (n1,n2) as seen from node n1.
  for(size_t i = 0; i < tlist.size(); ++i)
  {
    if(triangle_[tlist[i]].getClockwiseNode(n1) == n2) return tlist[i];
  }
  return -1;
}

//-----------------------------------------------------------------------------
void PrTriangulation_OP::buildTopology()
//-----------------------------------------------------------------------------
//   This algorithm finds the three neighbouring triangles (any of
//   which may be -1) of each triangle. If the number of neighbours
//   of each node is O(1) then this algorithm takes O(N) steps.
//   The algorithm also finds the "first" triangle of each node (also O(N)).
{
  int np = (int)node_.size();
  int nt = (int)triangle_.size();

  vector< vector<int> > incident_triang(np);
  int j;
  for(j=0; j< nt; j++)
  {
    incident_triang[triangle_[j].n1()].push_back(j);
    incident_triang[triangle_[j].n2()].push_back(j);
    incident_triang[triangle_[j].n3()].push_back(j);
  }

  for(j=0; j< np; j++)
  {
    if(incident_triang[j].size() == 0)
    {
      THROW("No triangle contains node in constructing PrTriangulation_OP.");
       // std::cerr << "Warning. No triangle contains node " << j;
       // std::cerr << " in constructing PrTriangulation_OP. ";
       // std::cerr << "Aborting from buildTopology()." << std::endl;
       // For the moment we allow this because of the wavelet 
       // application
       //return;
    }
    // @afr: I do not understand what the purpose of the test below
    // is. Anyway it didn't abort either! So I commented it all out.
//     if(incident_triang[j].size() == 1)
//     {
//        std::cerr << "Warning. Boundary node " << j;
//        std::cerr << " in constructing PrTriangulation_OP. ";
//        std::cerr << "Aborting from buildTopology()." << std::endl;
//     }
  }

  int n1,n2,n3;
  for(j=0; j< nt; j++)
  {
    n1 = triangle_[j].n1();
    n2 = triangle_[j].n2();
    n3 = triangle_[j].n3();

    triangle_[j].t3() = getNghrTriangle(n1,n2,incident_triang[n1]);
    triangle_[j].t1() = getNghrTriangle(n2,n3,incident_triang[n2]);
    triangle_[j].t2() = getNghrTriangle(n3,n1,incident_triang[n3]);
  }

  for(j=0; j< np; j++)
  {
    if (incident_triang[j].size()==0)
    {
      node_[j].tr()=-1;
      continue;
    }
    // store the "first" triangle -- there should be just one
    int ifirst = 0;
    int pfirst = -1;
    for(size_t k = 0; k < incident_triang[j].size(); ++k)
    {
      if(triangle_[incident_triang[j][k]].getRightTriangle(j) == -1)
      {
        pfirst = incident_triang[j][k];
        ifirst++;
      }
    }
    if(ifirst == 0)
    {
       // j is an interior node
       node_[j].tr() = incident_triang[j][0];
    }
    else if(ifirst == 1) 
    {
        // j is a boundary node
       node_[j].tr() = pfirst;
    }
    else
    {
       std::cerr << "Error. Node" << j;
       std::cerr << " has more than one first triangle";
       std::cerr << " in constructing PrTriangulation_OP" << std::endl;
       node_[j].tr() = pfirst;
    }
  }
}

// PUBLIC MEMBER FUNCTIONS

//-----------------------------------------------------------------------------
PrTriangulation_OP::PrTriangulation_OP(const double *xyz_points, int np,
                                       const int *triangles, int nt)
//-----------------------------------------------------------------------------
//   Construct the PrTriangulation_OP from an array of nodes and
//   an array of triangles.
//   The nodes array should contain x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,...
//   and the triangle array i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,....
//   where the i-th node is the 3D point (xi,yi,zi) and the
//   r-th triangle has the three nodes indexed ir, jr, and kr in
//   the node array. These nodes must be ordered anticlockwise (consistently
//   throughout the triangulation).
//   The length of the array "xyz_points" is 3 * np and
//   the length of the array "triangles" is 3 * nt.
//   The uv points for the nodes are set to zero.
//
//   The function buildTopology is called which finds the three
//   neighbouring triangles of each triangle and the first triangle
//   of each node.
{
  node_.resize(np);
  int j;
  for(j=0; j< np; j++)
  {
    node_[j].x() = xyz_points[3*j];
    node_[j].y() = xyz_points[3*j+1];
    node_[j].z() = xyz_points[3*j+2];
    node_[j].u() = 0.0;
    node_[j].v() = 0.0;
  }

  triangle_.resize(nt);
  int n1,n2,n3;
  for(j=0; j< nt; j++)
  {
    n1 = triangles[3*j];
    n2 = triangles[3*j+1];
    n3 = triangles[3*j+2];
    if(n1 < 0 || n1 >= np || n2 < 0 || n2 >= np || n3 < 0 || n3 >= np)
    {
       std::cerr << "Node indices out of bounds ";
       std::cerr << "in constructing PrTriangulation_OP" << std::endl;
       return;
    }
    if(n1 == n2 || n1 == n3 || n2 == n3)
    {
       std::cerr << "Non-distinct nodes ";
       std::cerr << "in constructing PrTriangulation_OP, triangle " 
		 << j << std::endl;
       std::cerr << n1 << "  " << n2 << "  " << n3 << std::endl;
       return;
    }

    triangle_[j].n1() = n1;
    triangle_[j].n2() = n2;
    triangle_[j].n3() = n3;
  }

  buildTopology();
}

//-----------------------------------------------------------------------------
PrTriangulation_OP::PrTriangulation_OP(const double *xyz_points, const double *uv_points,
                                       int np, const int *triangles, int nt)
//-----------------------------------------------------------------------------
//   Construct the PrTriangulation_OP from an array of nodes and
//   an array of triangles.
//   The xyz_points array should contain x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,...
//   The uv_points array should contain u0,v0,u1,v1,u2,v2,u3,...
//   and the triangle array i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,....
//   where the i-th node is the 3D point (xi,yi,zi) and the
//   r-th triangle has the three nodes indexed ir, jr, and kr in
//   the node array. These nodes must be ordered anticlockwise (consistently
//   throughout the triangulation).
//   The length of the array "xyzpoints" is 3 * np and
//   The length of the array "uv_points" is 2 * np and
//   the length of the array "triangles" is 3 * nt.
//
//   The function buildTopology is called which finds the three
//   neighbouring triangles of each triangle and the first triangle
//   of each node.
{
  node_.resize(np);
  int j;
  for(j=0; j< np; j++)
  {
    node_[j].x() = xyz_points[3*j];
    node_[j].y() = xyz_points[3*j+1];
    node_[j].z() = xyz_points[3*j+2];
    node_[j].u() = uv_points[2*j];
    node_[j].v() = uv_points[2*j+1];
  }

  triangle_.resize(nt);
  int n1,n2,n3;
  for(j=0; j< nt; j++)
  {
    n1 = triangles[3*j];
    n2 = triangles[3*j+1];
    n3 = triangles[3*j+2];
    if(n1 < 0 || n1 >= np || n2 < 0 || n2 >= np || n3 < 0 || n3 >= np)
    {
       std::cerr << "Node indices out of bounds ";
       std::cerr << "in constructing PrTriangulation_OP" << std::endl;
       return;
    }
    if(n1 == n2 || n1 == n3 || n2 == n3)
    {
       std::cerr << "Non-distinct nodes ";
       std::cerr << "in constructing PrTriangulation_OP" << std::endl;
       return;
    }

    triangle_[j].n1() = n1;
    triangle_[j].n2() = n2;
    triangle_[j].n3() = n3;
  }

  buildTopology();

}

//-----------------------------------------------------------------------------
PrTriangulation_OP::PrTriangulation_OP(PrExplicitConnectivity& op)
//-----------------------------------------------------------------------------
//   Construct the PrTriangulation_OP from a PrOrganizedPoints class
//   ASSUMING the PrOrganizedPoints class is a triangulation (every
//   face has three vertices).
{
  int npts = op.getNumNodes();
  node_.resize(npts);
  Vector3D nd;
  int j;
  for(j=0; j< npts; j++)
  {
    nd = op.get3dNode(j); 
    node_[j].x() = nd.x();
    node_[j].y() = nd.y();
    node_[j].z() = nd.z();
    node_[j].u() = op.getU(j);
    node_[j].v() = op.getV(j);
  }

  int nt = op.findNumFaces();
  triangle_.resize(nt);
  vector<int> neighbours;

  int i,deg,it = 0;
  for(j=0; j< npts; j++)
  {
    op.getNeighbours(j,neighbours); 
    deg = (int)neighbours.size();
    for(i=0; i<(deg-1); i++)
    {
      if(j < neighbours[i] && j < neighbours[i+1])
      {
        triangle_[it].n1() = j;
        triangle_[it].n2() = neighbours[i];
        triangle_[it].n3() = neighbours[i+1];
        it++;
      } 
    }
    if(!op.isBoundary(j))
    {
      if(j < neighbours[deg] && j < neighbours[0])
      {
        triangle_[it].n1() = j;
        triangle_[it].n2() = neighbours[deg];
        triangle_[it].n3() = neighbours[0];
        it++;
      } 
    }
  }

  buildTopology();
}


//----------------------------------------------------------------------------
void
PrTriangulation_OP::getNeighbours(int k, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the neighbours of the k-th node in:
//   1. any anticlockwise order if the k-th node is an interior node
//   2. the unique anticlockwise order if the k-th node is a boundary node.
{
  int tr1 = node_[k].tr();
  neighbours.clear();
  neighbours.push_back(triangle_[tr1].getAnticlockwiseNode(k));
  int tr,trNext;
  for(tr = node_[k].tr(), trNext = triangle_[tr].getLeftTriangle(k);
      trNext > -1 && trNext != tr1;
      tr = trNext, trNext = triangle_[tr].getLeftTriangle(k))
  {
    neighbours.push_back(triangle_[tr].getClockwiseNode(k));
  }
  if(trNext == -1) neighbours.push_back(triangle_[tr].getClockwiseNode(k));
       // k is a boundary node
  return;
}


//----------------------------------------------------------------------------
void PrTriangulation_OP::getTriangles(int k, 
				      Go::ScratchVect<int, 20>& triangles) const
//----------------------------------------------------------------------------
{
  int tr1 = node_[k].tr();
  triangles.clear();
  if (tr1>=0)
  {
    triangles.push_back(tr1);

    int trNext;
    for(trNext = triangle_[tr1].getLeftTriangle(k);
	trNext > -1 && trNext != tr1;
	trNext = triangle_[trNext].getLeftTriangle(k))
    {
      triangles.push_back(trNext);
    }
#if 0
    for (int i=0; i<triangle_.size(); i++)
      {
	if (triangle_[i].n1()==k ||
	    triangle_[i].n2()==k ||
	    triangle_[i].n2()==k)
	  {
	    int fnd=0;
	    for (int j=0; j<triangles.size(); j++)
	      {
		if(triangles[j]==i)
		  fnd=1;
	      }
	    if (!fnd)
	      {
		printf("PROBLEM \n");
	      }
	  }
      }
#endif
  }

}

//----------------------------------------------------------------------------
bool PrTriangulation_OP::isBoundary(int k) const
//-----------------------------------------------------------------------------
//   Given a node and its leading edge,
//   return 1 if the node is a boundary node and 0 otherwise.
{
  if(triangle_[node_[k].tr()].getRightTriangle(k) == -1) return true;
  else return false;
}

//-----------------------------------------------------------------------------
bool PrTriangulation_OP::swapTriangles(int t1, int t2)
//-----------------------------------------------------------------------------
//
//       t6                       t6
// n1 -------- n4           n1 -------- n4
//    |\     |                 |     /|
//    | \ t2 |                 | t1 / |
// t3 |  \   | t5    --->   t3 |   /  | t5
//    |   \  |                 |  /   |
//    | t1 \ |                 | / t2 |
//    |     \|                 |/     |
// n2 -------- n3           n2 -------- n3
//       t4                       t4
//
{
  int n1,n2,n3,n4;
  int t3,t4,t5,t6;
  if(triangle_[t1].t1() == t2)
  {
    t3 = triangle_[t1].t2();
    t4 = triangle_[t1].t3();
    n2 = triangle_[t1].n1();
    n3 = triangle_[t1].n2();
    n1 = triangle_[t1].n3();
  }
  else if(triangle_[t1].t2() == t2)
  {
    t3 = triangle_[t1].t3();
    t4 = triangle_[t1].t1();
    n2 = triangle_[t1].n2();
    n3 = triangle_[t1].n3();
    n1 = triangle_[t1].n1();
  }
  else if(triangle_[t1].t3() == t2)
  {
    t3 = triangle_[t1].t1();
    t4 = triangle_[t1].t2();
    n2 = triangle_[t1].n3();
    n3 = triangle_[t1].n1();
    n1 = triangle_[t1].n2();
  }
  else return false;

  if(triangle_[t2].t1() == t1)
  {
    t5 = triangle_[t2].t2();
    t6 = triangle_[t2].t3();
    n4 = triangle_[t2].n1();
  }
  else if(triangle_[t2].t2() == t1)
  {
    t5 = triangle_[t2].t3();
    t6 = triangle_[t2].t1();
    n4 = triangle_[t2].n2();
  }
  else //triangle_[t2].t3() == t1
  {
    t5 = triangle_[t2].t1();
    t6 = triangle_[t2].t2();
    n4 = triangle_[t2].n3();
  }

  triangle_[t1].init(n4,n1,n2,t3,t2,t6);
  triangle_[t2].init(n4,n2,n3,t4,t5,t1);

  if(t4 != -1)
  {
    if(triangle_[t4].t1() == t1) triangle_[t4].t1() = t2;
    else if(triangle_[t4].t2() == t1) triangle_[t4].t2() = t2;
    else triangle_[t4].t3() = t2;
  }
  else node_[n2].tr() = t2;

  if(t6 != -1)
  {
    if(triangle_[t6].t1() == t2) triangle_[t6].t1() = t1;
    else if(triangle_[t6].t2() == t2) triangle_[t6].t2() = t1;
    else triangle_[t6].t3() = t1;
  }
  else node_[n4].tr() = t1;
  return true;
}


//-----------------------------------------------------------------------------
void PrTriangulation_OP::swapNodes(int n1, int n2)
//-----------------------------------------------------------------------------
{
  Go::ScratchVect<int, 20> trv1, trv2;

  getTriangles(n1, trv1);
  getTriangles(n2, trv2);

  for (int *it=trv1.begin(); it!=trv1.end(); it++)
  {
    triangle_[*it].replaceNode(n1, -2);
  }
  for (int *it=trv2.begin(); it!=trv2.end(); it++)
  {
    triangle_[*it].replaceNode(n2, n1);
  }
  for (int *it=trv1.begin(); it!=trv1.end(); it++)
  {
    triangle_[*it].replaceNode(-2, n2);
  }
  std::swap(node_[n1], node_[n2]);
}

//-----------------------------------------------------------------------------
bool PrTriangulation_OP::splitTriangles (int t1, int t2, Vector3D& v)
//-----------------------------------------------------------------------------
// splits the common edge of triangles t1 and t2 and inserts the new
// vertex v there. Takes care of the mesh consistency
// local configuration looks as follows:
//
// v0-------------v2     v0------------v2     
//  |            /|       |\    tn    /|      
//  |  t1      /  |       |  \      /  |      
//  |        /    |       |    \  /    |      
//  |      /      |  ==>  | t1  vn  t2 |      
//  |    /        |       |    /  \    |      
//  |  /     t2   |       |  /      \  |      
//  |/            |       |/   tn+1   \|      
// v1-------------v3     v1------------v3     
//
{
  int tn = (int)triangle_.size();
  int vn = (int)node_.size();

  // assign v0..v3 such that v1,v2 define the common edge of t1 and t2,
  // v0 is opposite of t2 in t1 and v3 is opposite of t1 in t2
  int v0, v1, v2, v3;
  if (triangle_[t1].t1() == t2) {
    v0 = triangle_[t1].n1();
    v1 = triangle_[t1].n2();
    v2 = triangle_[t1].n3();
  }
  else
    if (triangle_[t1].t2() == t2) {
      v0 = triangle_[t1].n2();
      v1 = triangle_[t1].n3();
      v2 = triangle_[t1].n1();
    }
    else
      if (triangle_[t1].t3() == t2) {
	v0 = triangle_[t1].n3();
	v1 = triangle_[t1].n1();
	v2 = triangle_[t1].n2();
      }
      else
	return (false);
  v3 = triangle_[t2].getClockwiseNode(v2);

  // insert new vertex in node list
  node_.push_back( PrNode( v.x(), v.y(), v.z(), 0.0, 0.0, t1));

  // split t1 and insert new triangle tn
  int t1n = triangle_[t1].getOppositeTriangle(v1);
  triangle_.push_back( PrTriangle( v0, vn, v2,  t2, t1n, t1 ));
  triangle_[t1].replaceTriangle(t1n, tn);
  triangle_[t1].replaceTriangle(t2, tn+1);
  triangle_[t1].replaceNode(v2, vn);
  if (t1n > -1)
    triangle_[t1n].replaceTriangle(t1, tn);
  if (node_[v2].tr() == t1)
    node_[v2].tr() = tn;

  // split t2 and insert new triangle tn+1
  int t2n = triangle_[t2].getOppositeTriangle(v2);
  triangle_.push_back( PrTriangle( v3, vn, v1,  t1, t2n, t2 ));
  triangle_[t2].replaceTriangle(t2n, tn+1);
  triangle_[t2].replaceTriangle(t1, tn);
  triangle_[t2].replaceNode(v1, vn);
  if (t2n > -1)
    triangle_[t2n].replaceTriangle(t2, tn+1);
  if (node_[v1].tr() == t2)
    node_[v1].tr() = tn+1;

  return (true);
}

//-----------------------------------------------------------------------------
bool PrTriangulation_OP::splitVertex(int i, int j, vector<int> &new_nodes)
//-----------------------------------------------------------------------------
{
  new_nodes.clear();
  if(!isBoundary(i))
    return false;

  bool j_is_boundary=isBoundary(j);

  int curr=node_[i].tr();
  while (triangle_[curr].getClockwiseNode(i)!=j)
    {
      curr=triangle_[curr].getLeftTriangle(i);
    }
  
  int next=triangle_[curr].getLeftTriangle(i);
  int last=curr;
  int first=next;
  triangle_[last].replaceTriangle(first, -1);
  triangle_[first].replaceTriangle(last, -1);

  PrNode pr_node(node_[i].x(), node_[i].y(), node_[i].z(), 0.0, 0.0, next);
  int new_ind=(int)node_.size();
  new_nodes.push_back(new_ind);
  node_.push_back(pr_node);

  curr=first;
  while(curr!=-1)
    {
      triangle_[curr].replaceNode(i, new_ind);
      curr=triangle_[curr].getLeftTriangle(new_ind);
    }

  node_[j].tr()=last;
  {
    if (j_is_boundary)
      {
	  new_ind=(int)node_.size();
	int prev = -1;
	curr=first;
	while(curr!=-1)
	  {
	    prev=curr;
	    curr=triangle_[curr].getRightTriangle(j);
	    triangle_[prev].replaceNode(j,new_ind);
	  }

	PrNode new_node(node_[j].x(), node_[j].y(), node_[j].z(), 
			0.0, 0.0, prev);
	node_.push_back(new_node);
	new_nodes.push_back(new_ind);
      }
  }
  return(true);
}


//-----------------------------------------------------------------------------
void PrTriangulation_OP::printXYZTriangles(std::ostream& os, bool num) const
//-----------------------------------------------------------------------------
// Print out the triangles of the graph, useful for plotting
// If num = 1, print first the number of triangles.
{
  if(num) os << triangle_.size() << '\n';
  for(size_t i=0; i<triangle_.size(); i++)
  {
    node_[triangle_[i].n1()].printXYZ(os);
    node_[triangle_[i].n2()].printXYZ(os);
    node_[triangle_[i].n3()].printXYZ(os);
    os << "\n";
  }
}

//-----------------------------------------------------------------------------
void PrTriangulation_OP::printUVTriangles(std::ostream& os, bool num) const
//-----------------------------------------------------------------------------
// Print out the uv triangles of the graph, useful for plotting
// If num = 1, print first the number of triangles.
{
  if(num) os << triangle_.size() << '\n';
  for(size_t i=0; i<triangle_.size(); i++)
  {
    node_[triangle_[i].n1()].printUV(os);
    node_[triangle_[i].n2()].printUV(os);
    node_[triangle_[i].n3()].printUV(os);
    os << "\n";
  }
}

//-----------------------------------------------------------------------------
void PrTriangulation_OP::print(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  os << node_.size() << ' ' <<triangle_.size() << '\n';
  size_t i;
  for(i=0; i<node_.size(); i++) node_[i].print(os);
  os << "\n";
  for(i=0; i<triangle_.size(); i++) triangle_[i].print(os);
}

//-----------------------------------------------------------------------------
void PrTriangulation_OP::scan(std::istream& is)
//-----------------------------------------------------------------------------
{
  int numpnts,numtrs;
  is >> numpnts;
  is >> numtrs;
  node_.resize(numpnts);
  triangle_.resize(numtrs);
  int i;
  for(i=0; i<numpnts; i++) node_[i].scan(is);
  for(i=0; i<numtrs; i++) triangle_[i].scan(is);
}

//-----------------------------------------------------------------------------
void PrTriangulation_OP::printRawData(std::ostream& os) const
//-----------------------------------------------------------------------------
{
    os << node_.size() << ' ' << triangle_.size() << '\n';
  size_t i;
  for(i=0; i<node_.size(); i++)
    os << node_[i].x() << ' ' << node_[i].y() << ' ' << node_[i].z() << '\n';
  os << "\n";
  for(i=0; i<triangle_.size(); i++)
    os << triangle_[i].n1() << ' '
       << triangle_[i].n2() << ' '
       << triangle_[i].n3() << '\n';
}

//-----------------------------------------------------------------------------
void PrTriangulation_OP::scanRawData(std::istream& is)
//-----------------------------------------------------------------------------
{
  int numpnts = 0;
  int numtrs = 0;
  is >> numpnts;
  is >> numtrs;
  node_.resize(numpnts);
  triangle_.resize(numtrs);
  int i;
  for(i=0; i<numpnts; i++)
    is >> node_[i].x() >> node_[i].y() >> node_[i].z();
  int n1,n2,n3;
  for(i=0; i<numtrs; i++)
  {
    is >> n1 >> n2 >> n3;
    if(n1 < 0 || n1 >= numpnts || n2 < 0 ||
       n2 >= numpnts || n3 < 0 || n3 >= numpnts)
    {
       std::cerr << "Node indices out of bounds ";
       std::cerr << "in reading PrTriangulation_OP" << std::endl;
       return;
    }
    if(n1 == n2 || n1 == n3 || n2 == n3)
    {
       std::cerr << "Non-distinct nodes ";
       std::cerr << "in reading PrTriangulation_OP" << std::endl;
       return;
    }
    triangle_[i].n1() = n1;
    triangle_[i].n2() = n2;
    triangle_[i].n3() = n3;
  }
  buildTopology();
}

//-----------------------------------------------------------------------------
void PrTriangulation_OP::printUV(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  os << "320 1 0 0\n";
  os << node_.size() << ' ' << triangle_.size() << '\n';
  size_t i;
  for(i=0; i<node_.size(); i++)
    os << node_[i].u() << ' ' << node_[i].v() << " 0.0\n";
  os << "\n";
  for(i=0; i<triangle_.size(); i++)
      os << "3 " << triangle_[i].n1() << ' '
       << triangle_[i].n2() << ' '
       << triangle_[i].n3() << std::endl;
}

