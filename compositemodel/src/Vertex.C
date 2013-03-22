//===========================================================================
//                                                                           
// File: Vertex.C
//                                                                           
// Created: August 25, 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/geometry/GapRemoval.h"
#include "GoTools/geometry/SplineSurface.h"

using std::vector;
using std::make_pair;

namespace Go
{


//===========================================================================
    Vertex::Vertex(Point vertex_point)
//===========================================================================
	: vertex_point_(vertex_point)
    {
    }

//===========================================================================
    Vertex::Vertex(Point vertex_point, std::vector<ftEdge*> edges)
//===========================================================================
	: vertex_point_(vertex_point)
    {
	ftEdge *dummy = 0;
	for (size_t ki=0; ki<edges.size(); ki++)
	{
	    size_t kj;
	    for (kj=0; kj<ki; kj++)
		if (edges_[kj].first->twin() == edges[ki] && edges_[kj].second == 0)
		{
		    edges_[kj].second = edges[ki];
		    break;
		}
	    if (kj < ki)
		continue;
	    edges_.push_back(make_pair(edges[ki],dummy));
	}
    }

//===========================================================================
    Vertex::Vertex(Point vertex_point, ftEdge* edge)
//===========================================================================
	: vertex_point_(vertex_point)
    {
	ftEdge *dummy = 0;
	edges_.push_back(make_pair(edge,dummy));
    }

//===========================================================================
    Vertex::Vertex(ftEdge* edge, bool at_start)
//===========================================================================
    {
	ftEdge *dummy = 0;
	vertex_point_ = edge->point(at_start ? edge->tMin() : edge->tMax());
	edges_.push_back(make_pair(edge,dummy));
    }

//===========================================================================
    Vertex::~Vertex()
//===========================================================================
    {
    }

  
//===========================================================================
    void Vertex::joinVertex(shared_ptr<Vertex> other)
//===========================================================================
    {
	vector<ftEdge*> other_edges = other->allEdges();
	size_t ki, kj;
	ftEdge *dummy = 0;
	bool is_found;
	for (ki=0; ki<other_edges.size(); ki++)
	{
	  is_found = false;
	    for (kj=0; kj<edges_.size(); kj++)
	    {
		if (edges_[kj].first == other_edges[ki] || 
		    edges_[kj].second == other_edges[ki])
		  {
		    is_found = true;
		    break;
		  }
		else if  (edges_[kj].first->twin()  == other_edges[ki] && 
			  edges_[kj].second == 0)
		{
		    edges_[kj].second = other_edges[ki];
		    is_found = true;
		    break;
		}
	    }
	    for (kj++; kj<edges_.size(); ++kj)
	      {
		if (is_found && edges_[kj].first == other_edges[ki] &&
		    edges_[kj].second == 0)
		  {
		    edges_.erase(edges_.begin()+kj);
		    kj--;
		  }
	      }

	    if (is_found)
		continue;

	    edges_.push_back(make_pair(other_edges[ki],dummy));
	}
	
	double dist = vertex_point_.dist(other->getVertexPoint());
	if (dist > 0.01)
	  {
	    MESSAGE("Large vertex distance " << dist);
	  }
	vertex_point_ = 0.5*(vertex_point_ + other->getVertexPoint());
    }

//===========================================================================
    void Vertex::addEdge(ftEdge* edge)
//===========================================================================
    {

	ftEdge *dummy = 0;
	size_t kj;

	for (kj=0; kj<edges_.size(); kj++)
	{
	    if (edges_[kj].first == edge || edges_[kj].second == edge)
		break;
	}

	if (kj == edges_.size())
	  {
	    for (kj=0; kj<edges_.size(); kj++)
	      {
		if  (edges_[kj].first->twin()  == edge && edges_[kj].second == 0)
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
    void Vertex::removeEdge(ftEdge* edge)
//===========================================================================
    {
	    
	size_t kj;

	for (kj=0; kj<edges_.size();)
	{
	    if (edges_[kj].first == edge)
	    {
		
		edges_[kj].first = edges_[kj].second;
		edges_[kj].second = 0;
		if (edges_[kj].first == 0)
		    edges_.erase(edges_.begin()+kj);
		else kj++;
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
    bool Vertex::hasEdge(ftEdge *edge) const
//===========================================================================
    {
	for (size_t kj=0; kj<edges_.size(); kj++)
	  if (edges_[kj].first == edge || edges_[kj].second == edge)
	    return true;

	return false;
    }

//===========================================================================
    bool Vertex::hasFace(ftSurface *face) const
//===========================================================================
    {
	for (size_t kj=0; kj<edges_.size(); kj++)
	  if (edges_[kj].first->face() == face || 
	      (edges_[kj].second && edges_[kj].second->face() == face))
	    return true;

	return false;
    }

//===========================================================================
    bool Vertex::hasEdgeSingle(ftEdge *edge) const
//===========================================================================
    {
	for (size_t kj=0; kj<edges_.size(); kj++)
	  if (edges_[kj].first == edge && edges_[kj].second == 0)
	    return true;

	return false;
    }

//===========================================================================
    bool Vertex::meetInVertex(ftEdge *e1, ftEdge *e2) const
//===========================================================================
    {
      bool found1=false, found2=false;
      for (size_t kj=0; kj<edges_.size(); kj++)
	{
	  if (edges_[kj].first == e1 || edges_[kj].second == e1)
	    found1 = true;
	  if (edges_[kj].first == e2 || edges_[kj].second == e2)
	    found2 = true;
	}
      return (found1 && found2);
    }

  //===========================================================================
  bool Vertex::isBoundaryVertex() const
  //===========================================================================
  {
    for (size_t ki=0; ki<edges_.size(); ++ki)
      if (edges_[ki].first == 0 || edges_[ki].second == 0)
	return true;

    return false;
  }

  //===========================================================================
  bool Vertex::sameEdge(Vertex* other) const
  //===========================================================================
  {
    for (size_t ki=0; ki<edges_.size(); ++ki)
      {
	shared_ptr<Vertex> vx = edges_[ki].first->getOtherVertex(other);
	if (vx.get() == this)
	  return true;

	if (edges_[ki].second)
	  {
	     vx = edges_[ki].second->getOtherVertex(other);
	     if (vx.get() == this)
	       return true;
	  }
      }
    return false;
  }

  //===========================================================================
  bool Vertex::sameFace(Vertex* other) const
  //===========================================================================
  {
    vector<ftSurface*> faces1 = this->faces();
    vector<ftSurface*> faces2 = other->faces();
    for (size_t ki=0; ki<faces1.size(); ++ki)
      for (size_t kj=0; kj<faces2.size(); ++kj)
	if (faces1[ki] == faces2[kj])
	  return true;

    return false;
  }

  //===========================================================================
  bool Vertex::connectedToSameVertex(Vertex* other) const
  //===========================================================================
  {
    vector<ftEdge*> edges = other->uniqueEdges();
    for (size_t ki=0; ki<edges.size(); ++ki)
      {
	shared_ptr<Vertex> vx = edges[ki]->getOtherVertex(other);
	if (sameEdge(vx.get()))
	  return true;
      }
    return false;
  }
  //===========================================================================
  ftEdge* Vertex::getCommonEdge(Vertex* other) const
  //===========================================================================
  {
    vector<ftEdge*> edges = other->allEdges();
    for (size_t ki=0; ki<edges.size(); ++ki)
      {
	shared_ptr<Vertex> vx = edges[ki]->getOtherVertex(other);
	if (vx.get() == this)
	  return edges[ki];
      }
    return 0;
  }

  //===========================================================================
  ftEdge* Vertex::getCommonEdgeInFace(Vertex* other,
				      ftSurface* face) const
  //===========================================================================
  {
    vector<ftEdge*> edges = other->allEdges();
    for (size_t ki=0; ki<edges.size(); ++ki)
      {
	shared_ptr<Vertex> vx = edges[ki]->getOtherVertex(other);
	if (vx.get() == this && edges[ki]->face() == face)
	  return edges[ki];
      }
    return 0;
  }

  //===========================================================================
  vector<ftSurface*> Vertex::getCommonFaces(Vertex* other) const
  //===========================================================================
  {
    vector<ftSurface*> faces;
    vector<ftSurface*> faces1 = this->faces();
    vector<ftSurface*> faces2 = other->faces();
    for (size_t ki=0; ki<faces1.size(); ++ki)
      {
	for (size_t kj=0; kj<faces2.size(); ++kj)
	  if (faces1[ki] == faces2[kj])
	    {
	      faces.push_back(faces1[ki]);
	      break;
	    }
      }

    return faces;
  }

  //===========================================================================
  vector<ftEdge*> Vertex::getFaceEdges(ftSurface *face) const
  //===========================================================================
  {
    vector<ftEdge*> all_edges = allEdges();
    vector<ftEdge*> edges;
    for (size_t ki=0; ki<all_edges.size(); ++ki)
      {
	ftFaceBase *curr = all_edges[ki]->face();
	if (curr == face)
	  edges.push_back(all_edges[ki]);
      }
    return edges;
  }

  //===========================================================================
  bool  Vertex::isCornerInFace(ftSurface *face, double tol) 
  //===========================================================================
  {
    vector<ftEdge*> edges = getFaceEdges(face);
    
    if (edges.size() != 2)
      return true;

    // Check angle
    double t1 = edges[0]->parAtVertex(this);
    double t2 = edges[1]->parAtVertex(this);
    Point tan1 = edges[0]->tangent(t1);
    Point tan2 = edges[1]->tangent(t2);
    double tang = tan1.angle(tan2);
    
    return (tang > tol);
  }
 
  //===========================================================================
  vector<shared_ptr<Vertex> > Vertex::getNextVertex(ftSurface* face) const
  //===========================================================================
  {
    vector<shared_ptr<Vertex> > vxs;
    vector<ftEdge*> edges = getFaceEdges(face);
    for (size_t ki=0; ki<edges.size(); ++ki)
      {
	shared_ptr<Vertex> curr = edges[ki]->getOtherVertex(this);
	vxs.push_back(curr);
      }
    return vxs;
  }

  //===========================================================================
  void Vertex::disconnectTwin(ftEdge* edge)
  //===========================================================================
  {
    size_t kj;
    ftEdge *dummy = 0;

    for (kj=0; kj<edges_.size(); kj++)
      if (edges_[kj].first == edge || edges_[kj].second == edge)
	{
	  if (edges_[kj].second)
	    edges_.push_back(make_pair(edges_[kj].second,dummy));  
	  edges_[kj].second = 0;
	  break;
	}
  }

    
//===========================================================================
    vector<ftEdge*> Vertex::allEdges() const
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
    vector<ftEdge*> Vertex::uniqueEdges()
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
    vector<ftEdge*> Vertex::uniqueEdges(Body *bd)
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
	}
	return edges;
    }

//===========================================================================
    vector<ftEdge*> Vertex::freeEdges()
//===========================================================================
    {
	vector<ftEdge*> edges;
	for (size_t kj=0; kj<edges_.size(); kj++)
	{
	  if ((!edges_[kj].second) && (!edges_[kj].first->face()))
	    edges.push_back(edges_[kj].first);
	}
	return edges;
    }

//===========================================================================
  vector<pair<ftSurface*, Point> > Vertex::getFaces()
//===========================================================================
    {
      vector<ftEdge*> edges = allEdges();
      vector<pair<ftSurface*, Point> > faces;
      for (size_t ki=0; ki<edges.size(); ++ki)
	{
	  ftSurface* curr_face = edges[ki]->face()->asFtSurface();
	  size_t kj;
	  for (kj=0; kj<faces.size(); ++kj)
	    if (faces[kj].first == curr_face)
	      break;
	  if (kj < faces.size())
	    continue;
	  Point param = edges[ki]->faceParameter(edges[ki]->parAtVertex(this));
	  faces.push_back(make_pair(curr_face, param));
	}

      return faces;
    }

//===========================================================================
  vector<pair<ftSurface*, Point> > Vertex::getFaces(Body *bd)
//===========================================================================
    {
      vector<ftEdge*> edges = allEdges();
      vector<pair<ftSurface*, Point> > faces;
      for (size_t ki=0; ki<edges.size(); ++ki)
	{
	  ftSurface* curr_face = edges[ki]->face()->asFtSurface();
	  if (curr_face->getBody() != bd)
	    continue;

	  size_t kj;
	  for (kj=0; kj<faces.size(); ++kj)
	    if (faces[kj].first == curr_face)
	      break;
	  if (kj < faces.size())
	    continue;
	  Point param = edges[ki]->faceParameter(edges[ki]->parAtVertex(this));
	  faces.push_back(make_pair(curr_face, param));
	}

      return faces;
    }

//===========================================================================
  void Vertex::averageVertexPos()
//===========================================================================
    {
      vector<pair<ftSurface*, Point> > face_at_vx = getFaces();
      vector<shared_ptr<SplineSurface> > sfs;
      Point pos(0.0, 0.0, 0.0);
      vector<Point> param;
      for (size_t ki=0; ki<face_at_vx.size(); ++ki)
	{
	  shared_ptr<SplineSurface> tmp_sf =
	    dynamic_pointer_cast<SplineSurface,ParamSurface>(face_at_vx[ki].first->surface());
	  if (tmp_sf.get())
	    {
	      sfs.push_back(tmp_sf);
	      Point par = face_at_vx[ki].second;
	      Point tmp = sfs[ki]->ParamSurface::point(par[0], par[1]);
	      param.push_back(par);
	      pos += tmp;
	    }
	}
      if (sfs.size() == 0)
	return;

      pos /= (double)(sfs.size());

      double eps = 1.0e-15;  // Minor
      for (size_t ki=0; ki<sfs.size(); ++ki)
	GapRemoval::modifyAtVertex(sfs[ki], param[ki], pos, eps);
    }

//===========================================================================
  vector<ftSurface*> Vertex::faces() const
//===========================================================================
    {
      vector<ftSurface*> faces;
      for (size_t ki=0; ki<edges_.size(); ++ki)
	{
	  ftSurface* curr_face = edges_[ki].first->face()->asFtSurface();
	  size_t kj;
	  for (kj=0; kj<faces.size(); ++kj)
	    if (faces[kj] == curr_face)
	      break;
	  if (kj >= faces.size())
	    faces.push_back(curr_face);
	  
	  if (edges_[ki].second)
	    {
	      curr_face = edges_[ki].second->face()->asFtSurface();
	      for (kj=0; kj<faces.size(); ++kj)
		if (faces[kj] == curr_face)
		  break;
	      if (kj >= faces.size())
		faces.push_back(curr_face);
	    }
	}

      return faces;
    }

//===========================================================================
  vector<ftSurface*> Vertex::faces(Body *bd) const
//===========================================================================
    {
      vector<ftSurface*> faces;
      for (size_t ki=0; ki<edges_.size(); ++ki)
	{
	  ftSurface* curr_face = edges_[ki].first->face()->asFtSurface();
	  if (curr_face->getBody() != bd)
	    continue;

	  size_t kj;
	  for (kj=0; kj<faces.size(); ++kj)
	    if (faces[kj] == curr_face)
	      break;
	  if (kj >= faces.size())
	    faces.push_back(curr_face);
	  
	  if (edges_[ki].second)
	    {
	      curr_face = edges_[ki].second->face()->asFtSurface();
	      for (kj=0; kj<faces.size(); ++kj)
		if (faces[kj] == curr_face)
		  break;
	      if (kj >= faces.size())
		faces.push_back(curr_face);
	    }
	}

      return faces;
    }

//===========================================================================
  Point Vertex::getFacePar(ftSurface* face)
//===========================================================================
    {
      Point param;
      vector<ftEdge*> edges = allEdges();
      size_t ki;
      for (ki=0; ki<edges.size(); ++ki)
	if (edges[ki]->face() == face)
	  break;

      if (ki < edges.size())
	param = edges[ki]->faceParameter(edges[ki]->parAtVertex(this));

      return param;
    }

//===========================================================================
  vector<Body*> Vertex::getBodies()
//===========================================================================
    {
      // Fetch all faces
      vector<pair<ftSurface*, Point> > faces = getFaces();
      std::set<Body*> all_bodies;   // All bodies only once
      for (size_t ki=0; ki<faces.size(); ++ki)
	{
	  Body *curr = faces[ki].first->getBody();
	  if (curr)
	    all_bodies.insert(all_bodies.end(), curr);
	}
      vector<Body*> bodies(all_bodies.begin(), all_bodies.end());
      return bodies;
    }

//===========================================================================
    // Collect attached edges where the distance between the endpoints
    // are larger than the specified tolerance or where the curves meet with an angle that
    // are more than the kink tolerance, but less than the corner tolerance
    void Vertex::getEdgeDiscontinuities(vector<pair<ftEdge*, ftEdge*> >& gaps, double tol, 
					vector<pair<ftEdge*, ftEdge*> >& kinks, double angtol,
					double angmax) const
//===========================================================================
    {
	size_t ki, kj;
	ftEdge *e1, *e2;
	double dist, dist2;
	double ang, ang2;
	for (ki=0; ki<edges_.size(); ++ki)
	{
	    double t1 = edges_[ki].first->parAtVertex(this);
	    vector<Point> deriv1, deriv1_2;
	    edges_[ki].first->point(t1, 1, deriv1);
	    double t1_2 = -1;
	    if (edges_[ki].second)
	    {
		t1_2 = edges_[ki].second->parAtVertex(this);
		edges_[ki].second->point(t1_2, 1, deriv1_2);
	    }
	    for (kj=ki+1; kj<edges_.size(); ++kj)
	    {
		double t2 = edges_[kj].first->parAtVertex(this);
		vector<Point> deriv2, deriv2_2;
		edges_[kj].first->point(t2, 1, deriv2);
		double t2_2 = -1;
		if (edges_[kj].second)
		{
		    t2_2 = edges_[kj].second->parAtVertex(this);
		    edges_[kj].second->point(t2_2, 1, deriv2_2);
		}
		
		// Check distance
		e1 = edges_[ki].first;
		e2 = edges_[kj].first;
		dist = deriv1[0].dist(deriv2[0]);
		if (edges_[ki].second)
		{
		    dist2 = deriv1_2[0].dist(deriv2[0]);
		    if (dist2 > dist)
		    {
			dist = dist2;
			e1 = edges_[ki].second;
		    }

		    if (edges_[kj].second)
		    {
			dist2 = deriv1_2[0].dist(deriv2_2[0]);
			if (dist2 > dist)
			{
			    dist = dist2;
			    e2 = edges_[kj].second;
			}
		    }
		}
		if (edges_[kj].second)
		{
		    dist2 = deriv1[0].dist(deriv2_2[0]);
		    if (dist2 > dist)
		    {
			dist = dist2;
			e2 = edges_[kj].second;
		    }
		}

		if (dist > tol)
		    gaps.push_back(make_pair(e1, e2));
	
		// Check angle between tangents
		e1 = edges_[ki].first;
		e2 = edges_[kj].first;
		ang = deriv1[1].angle(deriv2[1]);
		if ((t1 == e1->tMin() && t2 == e2->tMin()) ||
		    (t1 == e1->tMax() && t2 == e2->tMax()))
		    ang = M_PI - ang;

		if (edges_[ki].second)
		{
		    ang2 = deriv1_2[1].angle(deriv2[1]);
		    if ((t1_2 == edges_[ki].second->tMin() && t2 == e2->tMin()) ||
			(t1_2 == edges_[ki].second->tMax() && t2 == e2->tMax()))
			ang2 = M_PI - ang2;
		    if (ang2 > ang && ang2 < angmax)
		    {
			ang = ang2;
			e1 = edges_[ki].second;
		    }

		    if (edges_[kj].second)
		    {
			ang2 = deriv1_2[1].angle(deriv2_2[1]);
		    if ((t1_2 == edges_[ki].second->tMin() && t2_2 == edges_[kj].second->tMin()) ||
			(t1_2 == edges_[ki].second->tMax() && t2_2 == edges_[kj].second->tMax()))
			ang2 = M_PI - ang2;
			if (ang2 > ang && ang2 <angmax)
			{
			    ang = ang2;
			    e2 = edges_[kj].second;
			}
		    }
		}
		if (edges_[kj].second)
		{
		    ang2 = deriv1[1].angle(deriv2_2[1]);
		    if ((t1 == edges_[ki].first->tMin() && t2_2 == edges_[kj].second->tMin()) ||
			(t1 == edges_[ki].first->tMax() && t2_2 == edges_[kj].second->tMax()))
			ang2 = M_PI - ang2;
		    if (ang2 > ang && ang2 < angmax)
		    {
			ang = ang2;
			e2 = edges_[kj].second;
		    }
		}

		if (ang > angtol && ang < angmax)
		    kinks.push_back(make_pair(e1, e2));
	    }
	}
    }

//===========================================================================
  void Vertex::reOrganize()
//===========================================================================
    {
      // Fetch all edges
      vector<ftEdge*> edges = allEdges();

      // Remove existing radial edge information
      edges_.clear();

      // Rebuild radial edge information
      size_t ki, kj;
      for (ki=0; ki<edges.size();)
	{
	  ftEdge *curr1 = edges[ki];
	  ftEdge *curr2 = NULL;
	  for (kj=ki+1; kj<edges.size(); ++kj)
	    if (curr1->twin() == edges[kj])
	      {
		curr2 = edges[kj];
		edges.erase(edges.begin()+kj);
		break;
	      }

	  edges_.push_back(make_pair(curr1, curr2));
	  edges.erase(edges.begin()+ki);
	}
    }

//===========================================================================
  bool Vertex::checkVertexTopology()
//===========================================================================
    {
      bool isOK = true;
      for (size_t ki=0; ki<edges_.size(); ++ki)
	{
	  shared_ptr<Vertex> v1;
	  shared_ptr<Vertex> v2;
	  edges_[ki].first->getVertices(v1, v2);
	  if (v1.get() != this && v2.get() != this)
	    {
	      std::cout << "Edge-vertex inconsistence, edge = " << edges_[ki].first;
	      std::cout << ", vertex = " << this << std::endl;
	      isOK = false;
	    }
	  bool edgeOK = edges_[ki].first->checkEdgeTopology();
	  if (!edgeOK)
	    isOK = false;
	  if (edges_[ki].second)
	    {
	      edges_[ki].second->getVertices(v1, v2);
	      if (v1.get() != this && v2.get() != this)
		{
		  std::cout << "Edge-vertex inconsistence, edge = " << edges_[ki].second;
		  std::cout << ", vertex = " << this << std::endl;
		  isOK = false;
		}
	      if (edges_[ki].first->twin() != edges_[ki].second || 
		  edges_[ki].second->twin() != edges_[ki].first)
		{
		  std::cout << "Error in vertex twin configuration, vertex = " << this;
		  std::cout <<", edges: " << edges_[ki].first << ", " << edges_[ki].second;
		  std::cout << std::endl;
		  isOK = false;
		}
	      edgeOK = edges_[ki].second->checkEdgeTopology();
	      if (!edgeOK)
		isOK = false;
	    }
	  else if (edges_[ki].first->twin())
	    {
	      std::cout << "Twin edge not represented in vertex, vertex = " << this;
	      std::cout <<", edges: " << edges_[ki].first << ", " << edges_[ki].first->twin();
	      std::cout << std::endl;
	      isOK = false;
	    }
	      
	}
      return isOK;
    }

} // namespace Go
