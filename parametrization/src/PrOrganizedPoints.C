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

#include "GoTools/parametrization/PrOrganizedPoints.h"
//#include "GoTools/parametrization/PrDijkstra.h"
#include <stack>

//----------------------------------------------------------------------------
PrOrganizedPoints::~PrOrganizedPoints()
//-----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
int PrOrganizedPoints::findNumBdyNodes() const
//-----------------------------------------------------------------------------
{
  int nBdy = 0;
  for(int i=0; i<getNumNodes(); i++)
       if(isBoundary(i)) nBdy++;
  return nBdy;
}

//----------------------------------------------------------------------------
int PrOrganizedPoints::findNumEdges() const
//-----------------------------------------------------------------------------
{
  int n = 0;
  vector<int> neighbours;
  for(int i=0; i<getNumNodes(); i++)
  {
    getNeighbours(i,neighbours);
    n+= (int)neighbours.size();
  }
  return n/2;
}


//----------------------------------------------------------------------------
bool PrOrganizedPoints::isMinimum(int i, vector<int>& face) 
//-----------------------------------------------------------------------------
{
  for(size_t j=0; j< face.size(); j++) if(i > face[j]) return false;
  return true;
}

//----------------------------------------------------------------------------
int PrOrganizedPoints::findNumComponents() const
//-----------------------------------------------------------------------------
{
  int n = getNumNodes();
  vector<int> component(n, 0);
  int i;
  int ic = 0;
  for(i=0; i<n; i++)
  {
    if(component[i] == 0)
    {
      ic++;
      labelNode(i,ic,component);
    }
  }
  return ic;
}

//----------------------------------------------------------------------------
void PrOrganizedPoints::labelNode(int i, int ic, vector<int>& component) const
//-----------------------------------------------------------------------------
{
   std::stack<int> ns;
   ns.push(i);
   int n, nj;
   vector<int> neighbours;
   while (ns.empty()==false)
   {
      n=ns.top();
      ns.pop();
      component[n]=ic;
      getNeighbours(n,neighbours);
      for(size_t j=0; j<neighbours.size(); j++)
      {
	 nj = neighbours[j];
	 if(component[nj] != ic)
	 {
	    ns.push(nj);
	 }
      }
   }
      
/*  int nj;
  component[i] = ic;
  vector<int> neighbours;
  getNeighbours(i,neighbours);
  for(int j=0; j<neighbours.size(); j++)
  {
    nj = neighbours[j];
    if(component[nj] != ic)
    {
      labelNode(nj,ic,component);
    }
  }
   */
}

//----------------------------------------------------------------------------
int PrOrganizedPoints::findNumBdyComponents() const
//-----------------------------------------------------------------------------
{
  int n = getNumNodes();
  vector<int> component(n, 0);
  int i;
  int ic = 0;
  for(i=0; i<n; i++)
  {
    if(isBoundary(i) && component[i] == 0)
    {
      ic++;
      labelBdyNode(i,ic,component);
    }
  }
  return ic;
}

//----------------------------------------------------------------------------
void PrOrganizedPoints::labelBdyNode(int i, int ic, 
				     vector<int>& component) const
//-----------------------------------------------------------------------------
{
  int nj;
  component[i] = ic;
  vector<int> neighbours;
  getNeighbours(i,neighbours);
  nj = neighbours[0];
  if(component[nj] != ic)
  {
    labelBdyNode(nj,ic,component);
  }
}

//----------------------------------------------------------------------------
int PrOrganizedPoints::indexComponents(vector<int>& component,
                                       vector<int>& newIndex) const
//-----------------------------------------------------------------------------
{
  int n = getNumNodes();
  component.resize(n);
  newIndex.resize(n);
  std::fill(component.begin(), component.end(), 0);
  int ic = 0;
  int index, i;
  for(i=0; i<n; i++)
  {
    if(component[i] == 0)
    {
      ic++;
      index = 0;
      labelNode(i,ic,index,component,newIndex);
    }
  }
  return ic;
}

//----------------------------------------------------------------------------
void PrOrganizedPoints::labelNode(int i, int ic, int& index,
                         vector<int>& component,
                         vector<int>& newIndex) const
//-----------------------------------------------------------------------------
{
  int nj;
  component[i] = ic;
  index++;
  newIndex[i] = index;
  vector<int> neighbours;
  getNeighbours(i,neighbours);
  for(size_t j=0; j<neighbours.size(); j++)
  {
    nj = neighbours[j];
    if(component[nj] != ic)
    {
      labelNode(nj,ic,index,component,newIndex);
    }
  }
}

//----------------------------------------------------------------------------
int PrOrganizedPoints::findIndex(Go::Vector3D& point) const
//-----------------------------------------------------------------------------
//   Given a point p in 3D find the index of the node
//   in the graph whose xyz point equals p.
{
  for(int i=0; i<getNumNodes(); i++)
  {
    if(get3dNode(i).dist(point) < 1e-15) return i;
  }
  return 0;
}


//----------------------------------------------------------------------------
void PrOrganizedPoints::topologicalDistToBdy(vector<int>& label) const
//-----------------------------------------------------------------------------
{
  int np = getNumNodes();
  label.resize(np);

  int i,j;
  for(i=0; i<np; i++)
  {
    if(!isBoundary(i)) label[i] = -1;
    else label[i] = 0;
  }

  int changes = 1;
  int distance = 0;
  vector<int> neighbours_;

  while(changes == 1)
  {
    changes = 0;
    for(i=0; i<np; i++)
      if (label[i] < 0)
      {
        getNeighbours(i,neighbours_);

        for(j=0; j< int(neighbours_.size()); j++)
        {
          if(label[neighbours_[j]] == distance)
          {
            label[i] = distance+1;
            changes = 1;
          }
        }
      }
    distance++;
  }

//  for(i=0; i<np; i++)
//    if(label[i] == -1) std::cout << "label = -1" << std::endl;

}

// //----------------------------------------------------------------------------
// void PrOrganizedPoints::geometricalDistToBdy(vector<double>& label) const
// //-----------------------------------------------------------------------------
// {
//   int np = getNumNodes();
//   label.resize(np);

//   Dijkstra* dijkstra = new Dijkstra;

//   dijkstra->setGraph(this);
//   dijkstra->initialize();

//   int i;
//   for(i=0; i<np; i++)
//     if(isBoundary(i))
//       dijkstra->setSource(i);

//   dijkstra->run();

//   double max_dist = 0.0;
//   for(i=0; i<np; i++)
//   {
//     label[i] = dijkstra->getDistance(i);  

//     if (label[i] > max_dist)
//       max_dist = label[i];
//   }

//   std::cout << "max geometrical distance = " << max_dist << std::endl;

//   delete dijkstra;
// }

//----------------------------------------------------------------------------
void
PrOrganizedPoints::getCommonNeighbours(int j, int k, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the common neighbours of the j-th and k-th node in some order
{
  neighbours.clear();

  vector<int> nj;
  vector<int> nk;

  getNeighbours (j, nj);
  getNeighbours (k, nk);

  for (size_t ij=0; ij<nj.size(); ij++)
    for (size_t ik=0; ik<nk.size(); ik++)
      if (nj[ij] == nk[ik])
	neighbours.push_back(nj[ij]);
    
  return;
}

//----------------------------------------------------------------------------
void
PrOrganizedPoints::get2Neighbours(int i, vector<int>& neighbours) const
//-----------------------------------------------------------------------------
//   Return the common neighbours of the j-th and k-th node in some order
//   including (!) boundary vertices
{
  neighbours.clear();

  vector<int> neighbours1;
  vector<int> neighbours2;

  getNeighbours(i,neighbours1);

  for (size_t j=0; j<neighbours1.size(); j++)
  {
    getNeighbours(neighbours1[j], neighbours2);
    
    for (size_t k=0; k<neighbours2.size(); k++)
    {
      bool newVertex = true;

      if (neighbours2[k] == i)
	newVertex = false;
      else
      {
	for (size_t l=0; l<neighbours.size(); l++)
	  if (neighbours[l] == neighbours2[k])
	  {
	    newVertex = false;
	    break;
	  }
      
	if (newVertex)
	  for (size_t l=0; l<neighbours1.size(); l++)
	    if (neighbours1[l] == neighbours2[k])
	    {
	      newVertex = false;
	      break;
	    }
      }

      if(newVertex)
	neighbours.push_back(neighbours2[k]);
    }
  }

  return;
}

//-----------------------------------------------------------------------------
void PrOrganizedPoints::printXYZNodes(std::ostream& os, bool num) const
//-----------------------------------------------------------------------------
// Print out the XYZ nodes of the graph.
// If num = true, print first the number of nodes.
{
    if(num) os << getNumNodes() << std::endl;
    Vector3D node;
    for(int i=0; i<getNumNodes(); i++)
	{
	    node = get3dNode(i);
//	    os << node;
	    node.write(os);
	}
}

//-----------------------------------------------------------------------------
void PrOrganizedPoints::printUVNodes(std::ostream& os, bool num) const
//-----------------------------------------------------------------------------
// Print out the UV nodes of the graph.
// If num = 1, print first the number of nodes.
{
    if(num) os << getNumNodes() << std::endl;
    for(int i=0; i<getNumNodes(); i++)
	{
	    os << getU(i) << ' ' << getV(i) << std::endl;
	}
}

//-----------------------------------------------------------------------------
void PrOrganizedPoints::printUVXYZNodes(std::ostream& os, bool num) const
//-----------------------------------------------------------------------------
// Print out the nodes of the graph.
// If num = 1, print first the number of nodes.
{
    if(num) os << getNumNodes() << std::endl;
    Vector3D node;
    for(int i=0; i<getNumNodes(); i++)
	{
	    node = get3dNode(i);
	    os << getU(i) << ' ' << getV(i) << std::endl;
//	    os << node;
	    node.write(os);
	}
}

//-----------------------------------------------------------------------------
void PrOrganizedPoints::printXYZEdges(std::ostream& os) const
//-----------------------------------------------------------------------------
// Print out the edges of the graph.
{
  vector<int> neighbours;
  int i,j,k,n;
  Vector3D node;

  for(i=0; i<getNumNodes(); i++)
  {
    getNeighbours(i,neighbours);
    n = (int)neighbours.size();
    for(k=0; k<n; k++)
    {
      j = neighbours[k];
      if(j > i)
      {
        node = get3dNode(i);
//        os << node;
        node.write(os);
	os << "  ";
        node = get3dNode(j);
//        os << node;
        node.write(os);
        os << std::endl;
        os << std::endl; // For some reason (my version of) gnuplot demands two separator lines.
      }
    }
  }
}

//-----------------------------------------------------------------------------
void PrOrganizedPoints::printUVEdges(std::ostream& os) const
//-----------------------------------------------------------------------------
// Print out the edges of the parametrization.
{
  vector<int> neighbours;
  int i,j,k,n;
  for(i=0; i<getNumNodes(); i++)
  {
    getNeighbours(i,neighbours);
    n = (int)neighbours.size();
    for(k=0; k<n; k++)
    {
      j = neighbours[k];
      if(j > i)
      {
	os << getU(i) << ' ' << getV(i) << std::endl;
	os << getU(j) << ' ' << getV(j) << std::endl;
        os << std::endl;
        os << std::endl;
      }
    }
  }
}

//----------------------------------------------------------------------------
void PrOrganizedPoints::printInfo(std::ostream& os) const
//-----------------------------------------------------------------------------
{
  os << "Number of nodes = " << getNumNodes() << std::endl;
  os << "Number of edges = " << findNumEdges() << std::endl;
  os << "Number of boundary nodes = " << findNumBdyNodes() << std::endl;
  os << "Number of connected components = ";
  os << findNumComponents() << std::endl;
  os << "Number of connected boundary components = ";
  os << findNumBdyComponents() << std::endl;
}

