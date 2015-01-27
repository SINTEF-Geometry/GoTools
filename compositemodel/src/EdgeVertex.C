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

#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/SplineSurface.h"

using std::vector;
using std::set;
using std::pair;
using std::make_pair;

namespace Go
{


//===========================================================================
    EdgeVertex::EdgeVertex(std::vector<ftEdge*> edges)
//===========================================================================
    {
	ftEdge *dummy = 0;
	// Check input
	for (size_t ki=0; ki<edges.size(); ++ki)
	  for (size_t kj=ki+1; kj<edges.size(); )
	    {
	      if (edges[ki] == edges[kj])
		edges.erase(edges.begin()+kj);
	      else
		kj++;
	    }

	for (size_t ki=0; ki<edges.size(); ++ki)
	{
	    size_t kj;
	    for (kj=0; kj<edges_.size(); kj++)
	      {
		ftSurface *face;
		face = edges_[kj].first->face()->asFtSurface();

		if (/*(face && face->twin() == edges[ki]->face() && 
		      edges_[kj].second == 0) ||*/
		    (edges_[kj].first->twin() == edges[ki] &&
		     edges[ki]->twin() == edges_[kj].first &&
		     edges_[kj].second == 0))
		{
		    edges_[kj].second = edges[ki];
		    break;
		}
	      }
	    if (kj < edges_.size())
		continue;
	    edges_.push_back(make_pair(edges[ki],dummy));
	}
    }

//===========================================================================
    EdgeVertex::EdgeVertex(ftEdge* edge)
//===========================================================================
    {
	ftEdge *dummy = 0;
	edges_.push_back(make_pair(edge,dummy));
    }

//===========================================================================
    EdgeVertex::~EdgeVertex()
//===========================================================================
    {
    }

//===========================================================================
    void EdgeVertex::addEdgeVertex(EdgeVertex* other)
//===========================================================================
    {
#ifdef DEBUG
      if (!checkRadialEdgeTopology())
	std::cout << "addEdgeVertex1" << std::endl;
#endif
      vector<ftEdge*> other_edges = other->allEdges();
      for (size_t ki=0; ki<other_edges.size(); ++ki)
	addEdge(other_edges[ki]);
#ifdef DEBUG
      if (!checkRadialEdgeTopology())
	std::cout << "addEdgeVertex2" << std::endl;
#endif
    }

//===========================================================================
  void EdgeVertex::addEdge(ftEdge* edge)
//===========================================================================
    {

	ftEdge *dummy = 0;
	size_t kj;

	for (kj=0; kj<edges_.size(); kj++)
	{
	  //	  ftSurface *face = edges_[kj].first->face()->asFtSurface();

	  if (edges_[kj].first == edge || edges_[kj].second == edge)
	    break;
	}

	if (kj == edges_.size())
	  {

	    for (kj=0; kj<edges_.size(); kj++)
	      {
		if  (edge->twin() == edges_[kj].first && 
		     edges_[kj].second == 0)
		  {
		    edges_[kj].second = edge;
		    break;
		  }
	      }
	  }

	if (kj == edges_.size())
	  edges_.push_back(make_pair(edge,dummy));  
    }

//===========================================================================
    void EdgeVertex::removeEdge(ftEdge* edge)
//===========================================================================
    {
	    
	size_t kj;

	for (kj=0; kj<edges_.size(); )
	{
	    if (edges_[kj].first == edge)
	    {
		
		edges_[kj].first = edges_[kj].second;
		edges_[kj].second = 0;
		if (edges_[kj].first == 0)
		    edges_.erase(edges_.begin()+kj);
		else
		  kj++;
		//		break;
	    }
	    else if (edges_[kj].second == edge)
	    {
		edges_[kj].second = 0;
		kj++;
		//		break;
	    }
	    else
	      kj++;
	}
    }

//===========================================================================
    vector<ftEdge*> EdgeVertex::allEdges() const
//===========================================================================
    {
	vector<ftEdge*> edges;
	for (size_t kj=0; kj<edges_.size(); kj++)
	{
	    edges.push_back(edges_[kj].first);
	    if (edges_[kj].second)
		edges.push_back(edges_[kj].second);
	}
	return edges;
    }

//===========================================================================
    vector<ftEdge*> EdgeVertex::allEdges(Body *bd) 
//===========================================================================
    {
	vector<ftEdge*> edges;
	for (size_t kj=0; kj<edges_.size(); kj++)
	{
	    if (edges_[kj].first->face())
	      {
		Body *bd2 = edges_[kj].first->face()->asFtSurface()->getBody();
		if (bd2 == bd)
		  edges.push_back(edges_[kj].first);
	      }
		
	    if (edges_[kj].second)
	      {
		if (edges_[kj].second->face())
		  {
		    Body *bd2 = edges_[kj].second->face()->asFtSurface()->getBody();
		    if (bd2 == bd)
		      edges.push_back(edges_[kj].second);
		  }
	      }
	}
	return edges;
    }

//===========================================================================
    vector<ftEdge*> EdgeVertex::uniqueEdges()
//===========================================================================
    {
	vector<ftEdge*> edges;
	for (size_t kj=0; kj<edges_.size(); kj++)
	{
	    edges.push_back(edges_[kj].first);
	}
	return edges;
    }

//===========================================================================
    vector<ftEdge*> EdgeVertex::uniqueEdges(Body *bd)
//===========================================================================
    {
	vector<ftEdge*> edges;
	for (size_t kj=0; kj<edges_.size(); kj++)
	{
	  if (!edges_[kj].first->face())
	    continue;

	  ftSurface *f1 = edges_[kj].first->face()->asFtSurface();
	  if (!f1)
	    continue;

	  if (f1->getBody() != bd)
	    continue;

	  edges.push_back(edges_[kj].first);
	}
	return edges;
    }

//===========================================================================
    int EdgeVertex::nmbUniqueEdges(Body *bd) const
//===========================================================================
    {
      int nmb = 0;
      for (size_t kj=0; kj<edges_.size(); kj++)
	{
	  if (!edges_[kj].first->face())
	    continue;
	  
	  ftSurface *f1 = edges_[kj].first->face()->asFtSurface();
	  if (!f1)
	    continue;
	  
	  if (f1->getBody() != bd)
	    continue;

	  nmb++;
	}
      return nmb;
    }
//===========================================================================
    bool EdgeVertex::hasEdge(ftEdge *edge) const
//===========================================================================
    {
	for (size_t kj=0; kj<edges_.size(); kj++)
	  if (edges_[kj].first == edge || edges_[kj].second == edge)
	    return true;

	return false;
    }

//===========================================================================
    bool EdgeVertex::hasEdgeSingle(ftEdge *edge) const
//===========================================================================
    {
	for (size_t kj=0; kj<edges_.size(); kj++)
	  if (edges_[kj].first == edge && edges_[kj].second == 0)
	    return true;

	return false;
    }


//===========================================================================
  vector<ftSurface*> EdgeVertex::getAdjacentFaces() const
//===========================================================================
  {
    // First fetch all edges
    vector<ftEdge*> edges = allEdges();
    
    // Fetch faces
    vector<ftSurface*> faces;
    faces.reserve(edges.size());
    for (size_t ki=0; ki<edges.size(); ++ki)
      {
	ftSurface* tmp = edges[ki]->face()->asFtSurface();
	// size_t kr;
	// for (kr=0; kr<faces.size(); ++kr)
	//   if (faces[kr] == tmp)
	//     break;

	// if (kr >= faces.size())
	  faces.push_back(tmp);
      }

    return faces;
  }


//===========================================================================
  vector<ftSurface*> EdgeVertex::getAdjacentFaces(Body *bd) const
//===========================================================================
  {
    // First fetch all edges
    vector<ftEdge*> edges = allEdges();
    
    // Fetch faces
    vector<ftSurface*> faces;
    for (size_t ki=0; ki<edges.size(); ++ki)
      {
	ftSurface* tmp = edges[ki]->face()->asFtSurface();
	if (tmp->getBody() == bd)
	  faces.push_back(tmp);
      }

    return faces;
  }

//===========================================================================
  vector<Body*> EdgeVertex::getAdjacentBodies() const
//===========================================================================
  {
    // First fetch all faces
    vector<ftSurface*> faces = getAdjacentFaces();
    set<Body*> all_bodies;
    for (size_t ki=0; ki<faces.size(); ++ki)
      {
	Body* body = faces[ki]->getBody();
	if (body)
	  all_bodies.insert(body);
      }

    vector<Body*> adj_bodies;
    adj_bodies.insert(adj_bodies.end(), all_bodies.begin(), all_bodies.end());
    return adj_bodies;
  }

//===========================================================================
  void EdgeVertex::disconnectTwin(ftEdge *e1, ftEdge *e2)
//===========================================================================
  {
    ftEdge *dummy = 0;
    for (size_t ki=0; ki<edges_.size(); ++ki)
      {
	if ((edges_[ki].first == e1 && edges_[ki].second == e2) ||
	    (edges_[ki].first == e2 && edges_[ki].second == e1))
	  {
	    edges_.push_back(make_pair(edges_[ki].second, dummy));
	    edges_[ki].second = dummy;
	    break;
	  }
      }
  }

//===========================================================================
  void EdgeVertex::organizeTwins()
//===========================================================================
  {
    // Fetch all edges
    ftEdge *dummy = 0;
    vector<ftEdge*> edges = allEdges();

    edges_.clear();
    for (size_t ki=0; ki<edges.size(); ki++)
      {
	size_t kj;
	for (kj=0; kj<edges_.size(); kj++)
	  {
	    ftSurface *face;
	    face = edges_[kj].first->face()->asFtSurface();

	    if (edges_[kj].first->twin() == edges[ki] &&
		 edges[ki]->twin() == edges_[kj].first &&
		 edges_[kj].second == 0)
	      {
		edges_[kj].second = edges[ki];
		break;
	      }
	  }
	if (kj < edges_.size())
	  continue;
	edges_.push_back(make_pair(edges[ki],dummy));
      }
  }

//===========================================================================
  void EdgeVertex::reOrganize()
//===========================================================================
  {
    int stat = 0;

    // Fetch all edges
    vector<ftEdge*> edges = allEdges();

    // Reorganization is limited within one Body or for edges that
    // do not belong to a body
    size_t ki, kj, kr;
    for (ki=0; ki<edges.size();)
      {
	ftSurface *curr_face = edges[ki]->face()->asFtSurface();
	Body *curr = (curr_face) ? curr_face->getBody() : NULL;

	// Find all edges pointing to this body
	vector<ftEdge*> curr_edges;
	curr_edges.push_back(edges[ki]);
	for (kj=ki+1; kj<edges.size();)
	  {
	    curr_face = edges[kj]->face()->asFtSurface();
	    Body *curr2 = (curr_face) ? curr_face->getBody() : NULL;
	    if (curr == curr2)
	      {
		curr_edges.push_back(edges[kj]);
		edges.erase(edges.begin()+kj);
	      }
	    else
	      kj++;
	  }
	edges.erase(edges.begin()+ki);

	// Reorganize current edges
	if (curr_edges.size() <= 2)
	  continue;  // Leave it as it is

	// Sort according to orientation of edge
	// Store also tangent and normal vectors
	vector<pair<Point,Point> > tnvec(curr_edges.size());
	kj = 0;
	double t1 = 0.5*(curr_edges[kj]->tMin() + curr_edges[kj]->tMax());
	vector<Point> der1(2);
	curr_edges[kj]->point(t1, 1, der1);
	Point norm = curr_edges[kj]->normal(t1);
	tnvec[kj] = make_pair(der1[1], norm);
	for (kr=1; kr<curr_edges.size(); ++kr)
	  {
	    double t2 = 0.5*(curr_edges[kr]->tMin() + 
			     curr_edges[kr]->tMax());
	    Point clo_pt;
	    double d2;
	    curr_edges[kr]->closestPoint(der1[0], t2, clo_pt, d2, &t2);
		
	    vector<Point> der2(2);
	    curr_edges[kr]->point(t2, 1, der2);
	    norm = curr_edges[kr]->normal(t2);
	    tnvec[kr] = make_pair(der2[1], norm);

	    if (der1[1]*der2[1] >= 0.0)
	      {
		// Same orientation
		std::swap(curr_edges[kj+1], curr_edges[kr]);
		std::swap(tnvec[kj+1], tnvec[kr]);
		kj++;
	      }
	  }

	size_t k2 = kj+1;  // Index to first edge of opposite orientation

	// Select one edge from the smallest group of edges with the same
	// orientation
	while (k2 >= 1 && curr_edges.size()-k2 >= 1)
	  {
	    size_t last;
	    if (k2 < curr_edges.size() - k2)
	      {
		kj = 0;
		kr = k2;
		last = curr_edges.size();
	      }
	    else
	      {
		kj = last = k2;
		kr = 0;
	      }

	    // Select the corresponding edge from the other group in such a 
	    // way that the angle between the normal vector of the other edge
	    // and the binormal vector of this edge is the smallest possible
	    // VSK. 102013. I Believe this is wrong, but I cannot remember
	    // how I thought when I implemented it. I make another try
	    // Point bivec = tnvec[kj].first.cross(tnvec[kj].second);
	    // size_t kmin = kr;
	    // double min_ang = bivec.angle2(tnvec[kr].second);
	    // for (kr++; kr<last; ++kr)
	    //   {
	    // 	double ang = bivec.angle2(tnvec[kr].second);
	    // 	if (ang < min_ang)
	    // 	  {
	    // 	    min_ang = ang;
	    // 	    kmin = kr;
	    // 	  }
	    //   }
	    // It might be that this will fail when more than
	    // 4 half edges meet in the same edge
	    Point norm1 = tnvec[kj].second;
	    double max_ang = norm1.angle(tnvec[kr].second);
	    max_ang = std::min(max_ang, fabs(M_PI - max_ang));
	    size_t kmin = kr;
	    for (kr++; kr<last; ++kr)
	      {
		double ang = norm1.angle(tnvec[kr].second);
		ang = std::min(ang, fabs(M_PI - ang));
		if (ang > max_ang)
		  {
		    max_ang = ang;
		    kmin = kr;
		  }
	      }

	    // A pair is found. Check if it is stored already
	    if (curr_edges[kj]->twin() != curr_edges[kmin])
	      {
		// Set twin pointers
		ftEdgeBase* twin1 = curr_edges[kj]->twin();
		ftEdgeBase* twin2 = curr_edges[kmin]->twin();
		if (twin1)
		  curr_edges[kj]->ftEdgeBase::disconnectTwin();
		if (twin2)
		  curr_edges[kmin]->ftEdgeBase::disconnectTwin();
		curr_edges[kj]->ftEdgeBase::connectTwin(curr_edges[kmin], stat);
		if (twin1 && twin2)
		  twin1->ftEdgeBase::connectTwin(twin2, stat);
	      }

	    // Update array. First remove information
	    for (kr=0; kr<edges_.size(); kr++)
	      {
		if (edges_[kr].first == curr_edges[kj])
		  {
		    if (edges_[kr].second)
		      {
			edges_[kr].first = edges_[kr].second;
			edges_[kr].second = NULL;
		      }
		    else
		      edges_.erase(edges_.begin()+kr);
		    break;
		  }
		else if (edges_[kr].second == curr_edges[kj])
		  {
		    edges_[kr].second = NULL;
		    break;
		  }
	      }

	    for (kr=0; kr<edges_.size(); kr++)
	      {
		if (edges_[kr].first == curr_edges[kmin])
		  {
		    if (edges_[kr].second)
		      {
			edges_[kr].first = edges_[kr].second;
			edges_[kr].second = NULL;
		      }
		    else
		      edges_.erase(edges_.begin()+kr);
		    break;
		  }
		else if (edges_[kr].second == curr_edges[kmin])
		  {
		    edges_[kr].second = NULL;
		    break;
		  }
	      }

	    // Add new pair
	    edges_.push_back(make_pair(curr_edges[kj], curr_edges[kmin]));

	    // Remove candidates
	    curr_edges.erase(curr_edges.begin()+std::max(kj, kmin));
	    curr_edges.erase(curr_edges.begin()+std::min(kj, kmin));
	    k2--;
	  }
      }

    // Make sure that edges with twins are placed first
    for (ki=0; ki<edges_.size(); ++ki)
      {
	if (edges_[ki].second)
	  continue;
	for (kj=ki+1; kj<edges_.size(); ++kj)
	  if (edges_[ki].second == NULL && edges_[kj].second != NULL)
	    {
	      std::swap(edges_[ki], edges_[kj]);
	      break;
	    }
      }

    // Reorganize vertices accordingly
    shared_ptr<Vertex> v1, v2;
    edges_[0].first->getVertices(v1, v2);
    v1->reOrganize();
    v2->reOrganize();
  }



//===========================================================================
  void EdgeVertex::splitAtVertex(shared_ptr<Vertex> v1,
				 shared_ptr<Vertex> v2, 
				 shared_ptr<Vertex> split_vx)
//===========================================================================
  {
    // Fetch edges. All edges are expected to end either at the original
    // vertcies v1 and v2 or at the new split vertex.
    vector<ftEdge*> edges = uniqueEdges();

    // Sort edges according to the position regarding the split vertex.
    // Edges are split if necessary
    vector<ftEdge*> edges1, edges2;
    int status = 1;
    for (size_t ki=0; ki<edges.size(); ++ki)
      {
	shared_ptr<Vertex> v3, v4;
	ftEdge *e2 = NULL;
	ftEdge *e1_1 = NULL, *e1_2 = NULL;
	edges[ki]->getVertices(v3, v4);
	if (edges[ki]->twin())
	  {
	    e2 = edges[ki]->twin()->geomEdge();
	    edges[ki]->ftEdgeBase::disconnectTwin();
	  }

	if ((v3.get() == v1.get() && v4.get() == split_vx.get()) ||
	    (v3.get() == split_vx.get() && v4.get() == v1.get()))
	  {
	    edges1.push_back(edges[ki]);
	    e1_1 = edges[ki];
	  }
	else if ((v3.get() == v2.get() && v4.get() == split_vx.get()) ||
		 (v3.get() == split_vx.get() && v4.get() == v2.get()))
	  {
	    edges2.push_back(edges[ki]);
	    e1_2 = edges[ki];
	  }
	else
	  {
	    shared_ptr<ftEdge> e1 = edges[ki]->splitAtVertex(split_vx);
	    shared_ptr<Vertex> other_vx = e1->getOtherVertex(split_vx.get());
	    if (other_vx.get() == v1.get())
	      {
		edges1.push_back(e1.get());
		edges2.push_back(edges[ki]);
		e1_1 = e1.get();
		e1_2 = edges[ki];
	      }
	    else
	      {
		edges1.push_back(edges[ki]);
		edges2.push_back(e1.get());
		e1_1 = edges[ki];
		e1_2 = e1.get();
	      }
	  }

	if (e2)
	  {
	    shared_ptr<Vertex> v5, v6;
	    e2->getVertices(v5, v6);
	    if ((v5.get() == v1.get() && v6.get() == split_vx.get()) ||
		(v5.get() == split_vx.get() && v6.get() == v1.get()))
	      {
		edges1.push_back(e2);
		if (e1_1)
		  e1_1->ftEdgeBase::connectTwin(e2, status);
	      }
	    else if ((v5.get() == v2.get() && v6.get() == split_vx.get()) ||
		     (v5.get() == split_vx.get() && v6.get() == v2.get()))
	      {
		edges2.push_back(e2);
		if (e1_2)
		  e1_2->ftEdgeBase::connectTwin(e2, status);
	      }
	    else
	      {	 
		shared_ptr<ftEdge> e3 = e2->splitAtVertex(split_vx);
		shared_ptr<Vertex> other_vx = e3->getOtherVertex(split_vx.get());
		if (other_vx.get() == v1.get())
		  {
		    edges1.push_back(e3.get());
		    edges2.push_back(e2);
		  }
		else
		  {
		    edges1.push_back(e2);
		    edges2.push_back(e3.get());
		  }
		if (e1_1)
		  e1_1->ftEdgeBase::connectTwin(edges1[edges1.size()-1],status);
		if (e1_2)
		  e1_2->ftEdgeBase::connectTwin(edges2[edges2.size()-1], status);
	      }
	  }
      }

    // Make new radial edges and connect all associated edges to those
    // This edge vertex will go out of scope
    if (edges1.size() > 0)
      {
	shared_ptr<EdgeVertex> edgevx = 
	  shared_ptr<EdgeVertex>(new EdgeVertex(edges1));
	for (size_t ki=0; ki<edges1.size(); ++ki)
	  edges1[ki]->setEdgeVertex(edgevx);
      }
	
    if (edges2.size() > 0)
      {
	shared_ptr<EdgeVertex> edgevx = 
	  shared_ptr<EdgeVertex>(new EdgeVertex(edges2));
	for (size_t ki=0; ki<edges2.size(); ++ki)
	  edges2[ki]->setEdgeVertex(edgevx);
      }
  }

 //===========================================================================
  void EdgeVertex::averageSplineEdges(double eps)
//===========================================================================
  {
    vector<ftEdge*> edges = allEdges();
    
    // Fetch first spline surface
    size_t ki;
    vector<shared_ptr<SplineSurface> > sfs;
    shared_ptr<CurveOnSurface> sfcv;
    vector<vector<pair<int,int> > > correspondance;
    vector<bool> same_orientation;
    for (ki=0; ki<edges.size(); ++ki)
      {
	shared_ptr<ParamSurface> sf1 = edges[ki]->face()->surface();
	shared_ptr<SplineSurface> sf2 = 
	  dynamic_pointer_cast<SplineSurface, ParamSurface>(sf1);
	shared_ptr<ParamCurve> bdcv = edges[ki]->geomCurve();
	shared_ptr<CurveOnSurface> sfcv = 
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv);

	if (sf2.get() && sfcv.get())
	  {
	    sfs.push_back(sf2);
	    break;
	  }
      }

    if (ki >= edges.size()-1)
      return;  // At most one spline surface. Nothing to average

    for (ki++; ki<edges.size(); ++ki)
      {
	shared_ptr<ParamSurface> sf1 = edges[ki]->face()->surface();
	shared_ptr<SplineSurface> sf2 = 
	  dynamic_pointer_cast<SplineSurface, ParamSurface>(sf1);
	shared_ptr<ParamCurve> bdcv = edges[ki]->geomCurve();
	shared_ptr<CurveOnSurface> sfcv2 = 
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv);

	if (sf2.get() && sfcv2.get())
	  {
	    int bd1, bd2;
	    bool same;
	    bool found = SurfaceTools::getSfAdjacencyInfo(sfs[0], sfcv, sf2, sfcv2, eps,
					    bd1, bd2, same);
	    if (!found)
	      continue;
	    
	    vector<pair<int,int> > enumeration;
	    bool pair = SurfaceTools::getCorrCoefEnum(sfs[0], sf2, bd1, bd2, same, 
					enumeration);
	    if (pair)
	      {
		sfs.push_back(sf2);
		correspondance.push_back(enumeration);
		same_orientation.push_back(same);
	      }
	  }
      }

    int dim = sfs[0]->dimension();
    for (size_t kj=0; kj<correspondance[0].size(); ++kj)
      {
	int ix1 = correspondance[0][kj].first;
	Point coef(sfs[0]->coefs_begin()+ix1*dim, 
		   sfs[0]->coefs_begin()+(ix1+1)*dim);
	for (ki=1; ki<sfs.size(); ++ki)
	  {
	    bool same = (same_orientation[ki-1] == same_orientation[0]);
	    int idx = (same) ? (int)kj : (int)correspondance[ki-1].size() - (int)kj - 1;
	    int ix2 = correspondance[ki-1][idx].second;
	    Point tmp(sfs[ki]->coefs_begin()+ix2*dim, 
		      sfs[ki]->coefs_begin()+(ix2+1)*dim);

	    coef += tmp;
	  }

	coef /= (double)(sfs.size());
	sfs[0]->replaceCoefficient(ix1, coef);

	for (ki=1; ki<sfs.size(); ++ki)
	  {
	    bool same = (same_orientation[ki-1] == same_orientation[0]);
	    int idx = (same) ? (int)kj : (int)correspondance[ki-1].size() - (int)kj - 1;
	    int ix2 = correspondance[ki-1][idx].second;
	    sfs[ki]->replaceCoefficient(ix2, coef);
	  }
      }
	    
  }

 //===========================================================================
  bool EdgeVertex::checkRadialEdgeTopology()
 //===========================================================================
  {
    bool isOK = true;
    for (size_t ki=0; ki<edges_.size(); ++ki)
      {
	if (edges_[ki].first->twin() &&
	    !(edges_[ki].first->twin() == edges_[ki].second))
	  {
	    std::cout << "Error in radial edge configuration = " << this;
	    std::cout << ", edges: " << edges_[ki].first << ", " << edges_[ki].second;
	    std::cout << std::endl;
	    isOK = false;
	  }

	bool edgeOK = edges_[ki].first->checkEdgeTopology();
	if (!edgeOK)
	  isOK = false;

	shared_ptr<EdgeVertex> radedg = edges_[ki].first->getEdgeMultiplicityInstance();
	if (radedg.get() != this)
	  {
	    std::cout << "Inconsistence in edgevertex pointer. Edge: " << edges_[ki].first;
	    std::cout << ", edgevertex: " << radedg.get() << std::endl;
	    isOK = false;
	  }

	if (edges_[ki].second &&
	    !(edges_[ki].second->twin() == edges_[ki].first))
	  {
	    std::cout << "Error in radial edge configuration = " << this;
	    std::cout << ", edges: " << edges_[ki].first << ", " << edges_[ki].second;
	    std::cout << std::endl;
	    isOK = false;
	  }

	if (edges_[ki].second)
	  {
	    edgeOK = edges_[ki].second->checkEdgeTopology();
	    if (!edgeOK)
	      isOK = false;

	    radedg = edges_[ki].second->getEdgeMultiplicityInstance();
	    if (radedg.get() != this)
	      {
		std::cout << "Inconsistence in edgevertex pointer. Edge: " << edges_[ki].second;
		std::cout << ", edgevertex: " << radedg.get() << std::endl;
		isOK = false;
	      }
	  }
      }
    return isOK;
  }

} // namespace Go
