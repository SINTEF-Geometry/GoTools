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

#include "GoTools/compositemodel/ftSurface.h"

#include "GoTools/compositemodel/ftSmoothSurf.h"
#include "GoTools/compositemodel/AdaptSurface.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/compositemodel/PointOnEdge.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#include "GoTools/creators/CoonsPatchGen.h"

#include "GoTools/geometry/LineCloud.h"
#include <fstream>

#include "GoTools/geometry/GapRemoval.h"

using std::vector;
using std::set;
using std::make_pair;

namespace Go
{

//---------------------------------------------------------------------------
ftSurface::ftSurface(shared_ptr<ParamSurface> sf, int id)
  : ftFaceBase(id), surf_(sf), prio_type_(ftNoType), twin_(0), body_(0)
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
ftSurface::ftSurface(shared_ptr<ParamSurface> sf, shared_ptr<Loop> loop, int id)
  : ftFaceBase(id), surf_(sf), prio_type_(ftNoType), twin_(0), body_(0)
//---------------------------------------------------------------------------
{
    addOuterBoundaryLoop(loop);
}

//---------------------------------------------------------------------------
ftSurface::~ftSurface()
//---------------------------------------------------------------------------
{
#ifdef DEBUG
  std::cout << "Delete face: " << this << std::endl;
  int stop_break = 1;
#endif
}


//===========================================================================
ftSurface* ftSurface::asFtSurface()
//===========================================================================
{
  return this;
}


//---------------------------------------------------------------------------
// Reset loop information
void ftSurface::clearInitialEdges()
//---------------------------------------------------------------------------
{
    boundary_loops_.clear();
}

//---------------------------------------------------------------------------
vector<shared_ptr<ftEdgeBase> > ftSurface::createInitialEdges(double degenerate_epsilon,
							      double kink, 
							      bool no_split)
//---------------------------------------------------------------------------
{
    vector<shared_ptr<ftEdgeBase> > edges;
    edges.reserve(4);
    
    // Store the tolerances in case this operation must be repeated from
    // inside this class
    degenerate_eps_ = degenerate_epsilon;
    kink_ = kink;

    if (boundary_loops_.size() == 0)
    {
      vector<CurveLoop> loops = SurfaceTools::allBoundarySfLoops(surf_, degenerate_epsilon);
	boundary_loops_.reserve(loops.size());

	// For every loop
	for (size_t ki = 0; ki < loops.size(); ++ki) 
	{
	  shared_ptr<Loop> curr_loop = shared_ptr<Loop>(new Loop(this, loops[ki], kink, true, no_split));
	    boundary_loops_.push_back(curr_loop);
	    vector<shared_ptr<ftEdgeBase> > curr_edges = curr_loop->getEdges();
	    edges.insert(edges.end(), curr_edges.begin(), curr_edges.end());
	}
    }
    else
    {
	for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) 
	{
	    shared_ptr<Loop> curr_loop = boundary_loops_[ki];
	    vector<shared_ptr<ftEdgeBase> > curr_edges = curr_loop->getEdges();
	    edges.insert(edges.end(), curr_edges.begin(), curr_edges.end());
	}
    }

    return edges;
}


// //---------------------------------------------------------------------------
// vector<shared_ptr<ftEdgeBase> > ftSurface::setOrientation(double degenerate_epsilon)
// //---------------------------------------------------------------------------
// {
//   if (is_turned_)
//     {
//       surf_->turnOrientation();
//       is_turned_ = false;
//     }
//   return createInitialEdges(degenerate_epsilon);
// }

//---------------------------------------------------------------------------
vector<shared_ptr<ftEdgeBase> > ftSurface::startEdges()
//---------------------------------------------------------------------------
{
    vector<shared_ptr<ftEdgeBase> > first_edges(boundary_loops_.size());
    for (size_t ki=0; ki<boundary_loops_.size(); ki++)
	first_edges[ki] = boundary_loops_[ki]->getEdge(0);

    return first_edges;
}


//---------------------------------------------------------------------------
  vector<Body*> ftSurface::getAdjacentBodies() const
//---------------------------------------------------------------------------
  {
    vector<Body*> bodies;
    if (body_)
      bodies.push_back(body_);
    if (twin_ && twin_->body_ && twin_->body_ != body_)
      bodies.push_back(twin_->body_);

    return bodies;
  }

//---------------------------------------------------------------------------
  bool ftSurface::connectTwin(ftSurface* newtwin, double tol, bool no_snap)
//---------------------------------------------------------------------------
  {
    // Assumes the same number of boundary loops
    if (nmbBoundaryLoops() != newtwin->nmbBoundaryLoops())
      return false;   // No connection performed

    ftEdgeBase *first1, *first2;
    bool same_dir;
    vector<ftEdgeBase*> loop_first1, loop_first2;
    vector<int> forward;
    for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
      {
	// Ensure correspondence between edges in the boundary loops
	  shared_ptr<Loop> other_loop = newtwin->getBoundaryLoop((int)ki);
	bool correspond = 
	  boundary_loops_[ki]->correspondingEdges(other_loop,
						  tol, first1, first2,
						  same_dir, no_snap);
	if (!correspond)
	  return false;
	loop_first1.push_back(first1);
	loop_first2.push_back(first2);
	forward.push_back(same_dir ? 1 : 0);
      }

    connectTwin(newtwin, loop_first1, loop_first2, forward);
    return true;
  }


//---------------------------------------------------------------------------
  void ftSurface::setTwin(ftSurface* newtwin)
//---------------------------------------------------------------------------
  {
    // Check if connection is possible
    if (twin_ && twin_ != newtwin)
      return;  // Twin already connected

    // Connect
    twin_ = newtwin;
    newtwin->twin_ = this;
  }

  //---------------------------------------------------------------------------
  void ftSurface::connectTwin(ftSurface* newtwin, vector<ftEdgeBase*> first1, 
			      vector<ftEdgeBase*> first2, vector<int> forward)
//---------------------------------------------------------------------------
  {
    // Check if connection is possible
    if (twin_ && twin_ != newtwin)
      return;  // Twin already connected

    // Connect
    twin_ = newtwin;
    newtwin->twin_ = this;

    for (size_t ki=0; ki<first1.size(); ++ki)
      {
	// Update edge information
	ftEdge *e1 = first1[ki]->geomEdge();
	ftEdge *e2 = first2[ki]->geomEdge();
	do
	  {
#ifdef DEBUG
	    if (e1->hasEdgeMultiplicity() && 
		!e1->getEdgeMultiplicityInstance()->hasEdge(e1))
	      std::cout << "Twin face 1_1. Radial edge missing" << std::endl;
	    if (e2->hasEdgeMultiplicity() && 
		!e2->getEdgeMultiplicityInstance()->hasEdge(e2))
	      std::cout << "Twin face 1_2. Radial edge missing" << std::endl;
#endif

	    e1->addEdgeMultiplicityInstance(e2);
	    e1->joinVertices(e2);

#ifdef DEBUG
	    if (e1->hasEdgeMultiplicity() && 
		!e1->getEdgeMultiplicityInstance()->hasEdge(e1))
	      std::cout << "Twin face 2_1. Radial edge missing" << std::endl;
	    if (e2->hasEdgeMultiplicity() && 
		!e2->getEdgeMultiplicityInstance()->hasEdge(e2))
	      std::cout << "Twin face 2_2. Radial edge missing" << std::endl;
#endif

	    e1 = e1->next()->geomEdge();
	    e2 = (forward[ki]) ? e2->next()->geomEdge() : 
	      e2->prev()->geomEdge();
	  }
	while (e1!=first1[ki] && e2!=first2[ki]);
      }
  }


//---------------------------------------------------------------------------
  void ftSurface::disconnectTwin()
//---------------------------------------------------------------------------
  {
    if (twin_)
      twin_->twin_ = NULL;
    twin_ = NULL;
  }

//---------------------------------------------------------------------------
  bool ftSurface::isSpline() const
//---------------------------------------------------------------------------
  {
    return surf_->isSpline();
  }

//---------------------------------------------------------------------------
/// Check if two faces are adjacent, and return information about
/// the edge along which this acjacency ocurrs
bool ftSurface::areNeighbours(ftSurface *other, shared_ptr<ftEdge>& edge1, 
			      shared_ptr<ftEdge>& edge2, int adj_idx) const
//---------------------------------------------------------------------------
{
  vector<pair<shared_ptr<ftEdgeBase>, shared_ptr<ftEdgeBase> > > found_edges;

  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      size_t nmb_edges = boundary_loops_[ki]->size();
      for (size_t kj=0; kj<nmb_edges; ++kj)
	{
	  shared_ptr<ftEdgeBase> curr = boundary_loops_[ki]->getEdge(kj);

	  vector<ftEdge*> edges;
	  if (curr->geomEdge()->hasEdgeMultiplicity())
	    edges = curr->geomEdge()->getEdgeMultiplicityInstance()->allEdges();
	    
	  if (!curr->twin() && edges.size()==0)
	    continue;

	  for (size_t kr=0; kr<other->boundary_loops_.size(); ++kr)
	    {
	      size_t nmb_edges2 = other->boundary_loops_[kr]->size();
	      for (size_t kh=0; kh<nmb_edges2; ++kh)
		{
		  shared_ptr<ftEdgeBase> curr2 = 
		    other->boundary_loops_[kr]->getEdge(kh);
		  // if (!curr2->twin())
		  //   continue;

		  if (curr->twin() == curr2.get() && curr2->twin() == curr.get())
		    {
		      // Decrement local adjecency index, return when desired index is found
		      if (adj_idx-- > 0)
			{
			  found_edges.push_back(make_pair(curr,curr2));
			  continue;
			}
		      edge1 = dynamic_pointer_cast<ftEdge,ftEdgeBase>(curr);
		      edge2 = dynamic_pointer_cast<ftEdge,ftEdgeBase>(curr2);
		      return true;
		    }
		}

	      for (size_t kr=0; kr<edges.size(); ++kr)
		{
		  if (edges[kr] == curr.get())
		    continue;
		  if (edges[kr]->face() == other &&
		      curr->geomEdge()->hasCommonVertices(edges[kr]))
		    {
		      // It is necessary to check for common vertices even if the two
		      // curves have the same radial edge and the face pointers are
		      // correct since the radial edge may have some inconsistencies
		      // cause by splitting of edges. The radial edge is not updated in 
		      // that situation to avoid inconsistencies in tpTopologyTable

		      size_t kh;
		      for (kh=0; kh<found_edges.size(); ++kh)
			if (found_edges[kh].first.get() == curr.get() ||
			    found_edges[kh].second.get() == curr.get())
			  break;

		      if (kh < found_edges.size())
			continue;

		      if (adj_idx-- > 0)
			continue;
		      
		      edge1 = 
			dynamic_pointer_cast<ftEdge,ftEdgeBase>(curr);
			      
		      // Find shared pointer
		      vector<shared_ptr<ftEdge> > f_edgs = other->getAllEdges();
		      for (size_t kv=0; kv<f_edgs.size(); ++kv)
			if (edges[kr] == f_edgs[kv].get())
			  {
			    for (kh=0; kh<found_edges.size(); ++kh)
			      if (found_edges[kh].first.get() == f_edgs[kv].get() ||
				  found_edges[kh].second.get() == f_edgs[kv].get())
				break;
			    if (kh < found_edges.size())
			      continue;
			    
			    edge2 = f_edgs[kv];
			    return true;
			  }
		    }
		}
	    }
	}
    }
  return false;
}
		       
//---------------------------------------------------------------------------
/// Number of edges along which the two faces are adjacent
int ftSurface::nmbAdjacencies(ftSurface *other) const
//---------------------------------------------------------------------------
{
  int nmb = 0;
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      size_t nmb_edges = boundary_loops_[ki]->size();
      for (size_t kj=0; kj<nmb_edges; ++kj)
	{
	  shared_ptr<ftEdgeBase> curr = boundary_loops_[ki]->getEdge(kj);

	  vector<ftEdge*> edges;
	  if (curr->geomEdge()->hasEdgeMultiplicity())
	    edges = curr->geomEdge()->getEdgeMultiplicityInstance()->allEdges();
	    
	  if (!curr->twin() && edges.size()==0)
	    continue;

	  for (size_t kr=0; kr<other->boundary_loops_.size(); ++kr)
	    {
	      size_t nmb_edges2 = other->boundary_loops_[kr]->size();
	      size_t kh;
	      for (kh=0; kh<nmb_edges2; ++kh)
		{
		  shared_ptr<ftEdgeBase> curr2 = 
		    other->boundary_loops_[kr]->getEdge(kh);
		  // if (!curr2->twin())
		  //   continue;

		  if (curr->twin() == curr2.get() && curr2->twin() == curr.get())
		    {
		      // Adjacency found
		      nmb++;
		    }
		}

	      for (kh=0; kh<edges.size(); ++kh)
		{
		  if (edges[kh] == curr.get())
		    continue;
		  if (edges[kh]->face() == other &&
		      curr->geomEdge()->hasCommonVertices(edges[kh]) &&
		      curr->twin() != edges[kh])
		    {
		      // Adjacency found
		      nmb++;
		    }
		}
	    }
	}
    }
  return nmb;
}

//---------------------------------------------------------------------------
// Number neighbouring faces common to the two given faces
int ftSurface::nmbNextNeighbours(ftSurface *other) const
//---------------------------------------------------------------------------
{
  int nmb = 0;
  // Fetch all neighbours of both faces
  vector<ftSurface*> neighbour1;
  getAdjacentFaces(neighbour1);

  vector<ftSurface*> neighbour2;
  other->getAdjacentFaces(neighbour2);

  // Compute the number of faces belonging to both groups of neighbours
  size_t ki, kj;
  for (ki=0; ki<neighbour1.size(); ++ki)
    {
      for (kj=0; kj<neighbour2.size(); ++kj)
	if (neighbour1[ki] == neighbour2[kj])
	  nmb++;
    }

  return nmb;
}

//---------------------------------------------------------------------------
Point ftSurface::point(double u, double v) const
//---------------------------------------------------------------------------
{
  return surf_->point(u,v);
}

//---------------------------------------------------------------------------
Point ftSurface::normal(double u, double v) const
//---------------------------------------------------------------------------
{
    Point pt(surf_->dimension());
    surf_->normal(pt, u, v);
//     if (is_turned_)
// 	pt *= -1;
    return pt;
}


//---------------------------------------------------------------------------
BoundingBox ftSurface::boundingBox()
//---------------------------------------------------------------------------
{
    return surf_->boundingBox();
}



//---------------------------------------------------------------------------
void ftSurface::setId(int id)
//---------------------------------------------------------------------------
{
    id_ = id;
}


//---------------------------------------------------------------------------
int ftSurface::getId()
//---------------------------------------------------------------------------
{
    return id_;
}

//---------------------------------------------------------------------------
ftTangPriority ftSurface::getPrioType() const
//---------------------------------------------------------------------------
{
    return prio_type_;
}

//---------------------------------------------------------------------------
void ftSurface::getError(double& max_error, double& mean_error)
//---------------------------------------------------------------------------
{
    max_error = mean_error = 0.0;
}

// //---------------------------------------------------------------------------
// void ftSurface::turnOrientation()
// //---------------------------------------------------------------------------
// {
//    is_turned_ = (is_turned_) ? false : true;
// }


// //---------------------------------------------------------------------------
// bool ftSurface::getOrientation()
// //---------------------------------------------------------------------------
// {
//     return is_turned_;
// }


//---------------------------------------------------------------------------
void ftSurface::updateBoundaryLoops(shared_ptr<ftEdgeBase> new_edge)
//---------------------------------------------------------------------------
{
    for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
	boundary_loops_[ki]->updateLoop(new_edge);
}

//---------------------------------------------------------------------------
  void ftSurface::isolateFace()
//---------------------------------------------------------------------------
{
  // Disconnect edges
  vector<shared_ptr<ftEdgeBase> > tmp_edges;
  tmp_edges = createInitialEdges();

  size_t ki, kj;
  for (ki=0; ki<tmp_edges.size(); ++ki)
    {
      if (tmp_edges[ki]->twin())
	tmp_edges[ki]->disconnectTwin();

      // Disconnect radial edge
      shared_ptr<EdgeVertex> radial_edge = 
	tmp_edges[ki]->geomEdge()->getEdgeMultiplicityInstance();
      if (radial_edge.get())
	{
	  radial_edge->removeEdge(tmp_edges[ki]->geomEdge());
	  tmp_edges[ki]->geomEdge()->removeEdgeVertex();
	  //#ifdef DEBUG
	  // TEST
	  // if (!radial_edge->checkTwins())
	  //   std::cout << "Isolate face " << ki <<". Radial edge inconsistency" << std::endl;
	  //#endif
	}
    }
	  

  // Disconnect vertices
  vector<shared_ptr<Vertex> > vx = vertices();
  for (ki=0; ki<vx.size(); ++ki)
    {
      vector<ftEdge*> curr_edges = vx[ki]->getFaceEdges(this);
      for (kj=0; kj<curr_edges.size(); ++kj)
	vx[ki]->removeEdge(curr_edges[kj]);

      shared_ptr<Vertex> new_vx = 
	shared_ptr<Vertex>(new Vertex(vx[ki]->getVertexPoint(), curr_edges));
      for (kj=0; kj<curr_edges.size(); ++kj)
	curr_edges[kj]->replaceVertex(vx[ki], new_vx);
    }

}

//---------------------------------------------------------------------------
ftMessage ftSurface::createSurf(double& max_error, double& mean_error)
//---------------------------------------------------------------------------
{
  // VSK, 0108. Nothing currently
  ftMessage status;
  max_error = mean_error = 0.0;
  return status;
}


//---------------------------------------------------------------------------
void ftSurface::addBoundaryLoops(vector<shared_ptr<Loop> >& bd_loops)
{
    // Check input
    size_t ki;
    for (ki=0; ki<bd_loops.size(); ki++)
	ALWAYS_ERROR_IF(!bd_loops[ki]->isFaceConsistent(), 
			"Inconsistency in face pointers");

    // Here a function checking the consistency of the loops and 
    // consistency between the loops and the ParamSurface surf_, is
    // required

    // Existing boundary loops are removed
    boundary_loops_.clear();

    for (ki=0; ki<bd_loops.size(); ki++)
	boundary_loops_.push_back(bd_loops[ki]);

}


//---------------------------------------------------------------------------
void ftSurface::addOuterBoundaryLoop(shared_ptr<Loop> outer_loop)
{
    // Check input
    ALWAYS_ERROR_IF(!outer_loop->isFaceConsistent(), 
		    "Inconsistency in face pointers");

    const ftFaceBase* face = outer_loop->getFace();
    ALWAYS_ERROR_IF(face != 0 && face != this, "Inconsistency in face pointers");

    if (!face)
	outer_loop->setFace(this);

    // Here a function checking the consistency of the loop and 
    // consistency between the loop and the ParamSurface surf_, is
    // required

    // Existing boundary loops are removed
    boundary_loops_.clear();

    boundary_loops_.push_back(outer_loop);
}

//---------------------------------------------------------------------------
  int ftSurface::nmbEdges() const
//---------------------------------------------------------------------------
  {
    int nmb = 0;
    for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
	nmb += (int)boundary_loops_[ki]->size();
    return nmb;
  }

//---------------------------------------------------------------------------
  vector<shared_ptr<ftEdge> > ftSurface::getAllEdges() const
//---------------------------------------------------------------------------
  {
    vector<shared_ptr<ftEdge> > edges;
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) 
      {
	shared_ptr<Loop> curr_loop = boundary_loops_[ki];
	vector<shared_ptr<ftEdgeBase> > curr_edges = curr_loop->getEdges();
	for (size_t kj=0; kj<curr_edges.size(); ++kj)
	  {
	    shared_ptr<ftEdge> curr = 
	      dynamic_pointer_cast<ftEdge,ftEdgeBase>(curr_edges[kj]);
	    if (curr.get())
	      edges.push_back(curr);  // Should always be the case
	  }
      }
    return edges;
  }

//---------------------------------------------------------------------------
  vector<shared_ptr<ftEdge> > ftSurface::getAllEdges(int loop_idx) const
//---------------------------------------------------------------------------
  {
    vector<shared_ptr<ftEdge> > edges;
    if (loop_idx < 0 || loop_idx >= (int)(boundary_loops_.size()))
      return edges;

    shared_ptr<Loop> curr_loop = boundary_loops_[loop_idx];
    vector<shared_ptr<ftEdgeBase> > curr_edges = curr_loop->getEdges();
    for (size_t kj=0; kj<curr_edges.size(); ++kj)
      {
	shared_ptr<ftEdge> curr = 
	  dynamic_pointer_cast<ftEdge,ftEdgeBase>(curr_edges[kj]);
	if (curr.get())
	  edges.push_back(curr);  // Should always be the case
      }

    return edges;
  }

//---------------------------------------------------------------------------
  vector<ftEdge*> ftSurface::getAllEdgePtrs() const
//---------------------------------------------------------------------------
  {
    vector<shared_ptr<ftEdge> > edges = getAllEdges();
    vector<ftEdge*> edge_ptrs(edges.size());
    for (size_t ki=0; ki<edges.size(); ++ki)
      edge_ptrs[ki] = edges[ki].get();
    return edge_ptrs;
  }

//---------------------------------------------------------------------------
  vector<ftEdge*> ftSurface::getAllEdgePtrs(int loop_idx) const
//---------------------------------------------------------------------------
  {
    vector<shared_ptr<ftEdge> > edges = getAllEdges(loop_idx);
    vector<ftEdge*> edge_ptrs(edges.size());
    for (size_t ki=0; ki<edges.size(); ++ki)
      edge_ptrs[ki] = edges[ki].get();
    return edge_ptrs;
  }

///---------------------------------------------------------------------------
int ftSurface::nmbOuterBdCrvs(double gap, double neighbour, double angtol) const
//---------------------------------------------------------------------------
  {
    vector<shared_ptr<ftEdgeBase> > edges = boundary_loops_[0]->getEdges();
    
    int ki, kj;
    int nmbedge = (int)edges.size();
    int nmb_bds = 0;
    for (ki=nmbedge-1, kj=0; kj<nmbedge; ki=(ki+1)%nmbedge, kj++)
	{
	  Point p1 = edges[ki]->point(edges[ki]->tMax());
	  Point p2 = edges[kj]->point(edges[kj]->tMin());
	  Point tan1 = edges[ki]->tangent(edges[ki]->tMax());
	  Point tan2 = edges[kj]->tangent(edges[kj]->tMin());

	  double dist = p1.dist(p2);
	  double ang = tan1.angle(tan2);
	  if (dist > neighbour || ang > angtol) 
	    nmb_bds++;
	}
    return nmb_bds;
  }

//---------------------------------------------------------------------------
shared_ptr<ParamSurface>
ftSurface::getUntrimmed(double gap, double neighbour, double angtol, bool only_corner)
//---------------------------------------------------------------------------
{
  shared_ptr<ParamSurface> surf2;  // Output surface

  // Check if the face is really regular
  int nmb = nmbOuterBdCrvs(gap, neighbour, angtol);
  if (boundary_loops_.size() > 1 || nmb > 4)
    return surf2;  // Not possible to replace the surface with a 
                   // non-trimmed one

#ifdef DEBUG
  std::ofstream of("trimmed.g2");
  surf_->writeStandardHeader(of);
  surf_->write(of);
#endif

  // Check if the surface has to be replaced at all
  if (surf_->instanceType() == Class_SplineSurface)
    return surf_;  // Nothing to do

  // Fetch boundary curve information
  vector<pair<shared_ptr<ParamCurve>,shared_ptr<ParamCurve> > > cvs;
  getBoundaryCurves(angtol, cvs, only_corner);
  if (cvs.size()!=4)
    return surf2;

  // For each parameter direction, approximate the curves in the same
  // or refined spline space
  vector<shared_ptr<SplineCurve> > bd_cvs;
  for (int ki=0; ki<2; ++ki)
    getApproxCurves(cvs.begin()+2*ki, 2, bd_cvs, gap);

#ifdef DEBUG
  std::ofstream of3("approx_cvs.g2");
  for (size_t k3=0; k3<bd_cvs.size(); ++k3)
    {
       bd_cvs[k3]->writeStandardHeader(of3);
      bd_cvs[k3]->write(of3);
    }
#endif
  // Make surface, first a Coons patch approximating the boundary curves
  // Prepare orientation
  Point pt1 = bd_cvs[0]->ParamCurve::point(bd_cvs[0]->endparam());
  for (int ki=1; ki<4; ++ki)
    {
      Point pt2 = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->startparam());
      Point pt3 = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->endparam());
      if (pt1.dist(pt3) < pt1.dist(pt2))
	bd_cvs[ki]->reverseParameterDirection();
      pt1 = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->endparam());
    }
	  
  vector<shared_ptr<ParamCurve> > tmp_cvs(bd_cvs.begin(), bd_cvs.end());
  CurveLoop boundary(tmp_cvs, gap);
  shared_ptr<SplineSurface> init_sf(CoonsPatchGen::createCoonsPatch(boundary));

  // Approximate the original surface
  surf2 = AdaptSurface::adaptSurface(surf_, init_sf, gap);

  #ifdef DEBUG
  std::ofstream of2("untrimmed.g2");
  surf2->writeStandardHeader(of2);
  surf2->write(of2);
#endif

  return surf2;
}

//---------------------------------------------------------------------------
void 
ftSurface::getApproxCurves(vector<pair<shared_ptr<ParamCurve>,
			   shared_ptr<ParamCurve> > >::iterator cvs_in,
			   int nmb_cvs, 
			   vector<shared_ptr<SplineCurve> >& cvs_out,
			   double tol)
//---------------------------------------------------------------------------
{
  vector<shared_ptr<SplineCurve> > tmp_cvs(nmb_cvs);
  vector<shared_ptr<ParamCurve> > init_cvs;
  vector<BsplineBasis> crv_basis;
  vector<int> bb_idx;
  int ki;
  int min_cont = 2;
  
  for (ki=0; ki<nmb_cvs; ++ki)
    {
      shared_ptr<SplineCurve> spcv1, spcv2;

      // Check if the first curve can be used as it is
       if (cvs_in[ki].first.get())
	{
	    spcv1 = dynamic_pointer_cast<SplineCurve,ParamCurve>(cvs_in[ki].first);
	  if (!spcv1.get())
	    {
	      shared_ptr<CurveOnSurface> sfcv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cvs_in[ki].first);
	      if (sfcv.get() && sfcv->isConstantCurve())
		{
		  spcv1 = dynamic_pointer_cast<SplineCurve,ParamCurve>(sfcv->spaceCurve());
		}
	    }
	}

      // Check continuity
      if (spcv1.get() && (spcv1->basis().getMinContinuity() < min_cont /*spcv1->order()-2*/ ||
			  spcv1->rational()))
	spcv1.reset();

      // Check if the second curve can be used as it is
     if (cvs_in[ki].second.get())
	{
	  spcv2 = dynamic_pointer_cast<SplineCurve,ParamCurve>(cvs_in[ki].second);
	  if (!spcv2.get())
	    {
	      shared_ptr<CurveOnSurface> sfcv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cvs_in[ki].second);
	      if (sfcv.get() && sfcv->isConstantCurve())
		{
		  spcv2 = dynamic_pointer_cast<SplineCurve,ParamCurve>(sfcv->spaceCurve());
		}
	    }
	}

      // Check continuity
      if (spcv2.get() && (spcv2->basis().getMinContinuity() < min_cont /*spcv2->order()-2*/ ||
			  spcv2->rational()))
	spcv2.reset();


      // Fetch current information
      if ((spcv1.get() && spcv2.get() && spcv1->numCoefs() <= spcv2->numCoefs())
	  || (spcv1.get() && !spcv2.get()))
	{
	  tmp_cvs[ki] = spcv1;
	  crv_basis.push_back(spcv1->basis());
	  bb_idx.push_back(ki);
	}
      else if (spcv2.get())
	{
	  tmp_cvs[ki] = spcv2;
	  crv_basis.push_back(spcv2->basis());
	  bb_idx.push_back(ki);
	}

      // Always store initial curves. They will be removed later if they
      // are kept without approximation
      if (cvs_in[ki].first.get() && !(spcv2.get() &&
				      tmp_cvs[ki].get() == spcv2.get()))
	init_cvs.push_back(cvs_in[ki].first);
      else
	init_cvs.push_back(cvs_in[ki].second);
    }

  // Approximate curves in the same spline space up to possible
  // refinements
  vector<shared_ptr<SplineCurve> > app_cvs;
  if (crv_basis.size() > 0)
    {
      // Define initial knot vector as the union of the existing
      // knot vectors
      double start = crv_basis[0].startparam();
      double end = crv_basis[0].endparam();
      int order = crv_basis[0].order();
      for (ki=1; ki<(int)crv_basis.size(); ++ki)
	{
	  order = std::max(order, crv_basis[ki].order());
	  start += crv_basis[ki].startparam();
	  end += crv_basis[ki].endparam();
	}
      start /= (double)(crv_basis.size());
      end /= (double)(crv_basis.size());
      for (ki=0; ki<(int)crv_basis.size(); ++ki)
	{
	  crv_basis[ki].rescale(start, end);
	  if (crv_basis[ki].order() < order)
	    crv_basis[ki].increaseOrder(order);
	}

      vector<double> knots;
      GeometryTools::makeUnionKnots(crv_basis, tol, knots);

      // Check the distribution of knotw in the union knot vector
      int nmb_basis = (int)knots.size()-order;
      double tdel = (knots[nmb_basis] - knots[order-1])/(double)nmb_basis;
      tdel /= 5.0;
      double prev = knots[order-1];
      int kj;
      for (kj=order; kj<=nmb_basis; ++kj)
	if (knots[kj] > prev)
	  {
	    if (knots[kj] - prev < tdel)
	      break;
	    else
	      prev = knots[kj];
	  }
      if (kj <= nmb_basis)
	{
	  // Bad knot distribution. Use at most one kept curve
	  size_t kj;
	  if (bb_idx.size() == 1)
	    {
	      knots.erase(knots.begin()+order, knots.begin()+nmb_basis-1);
	      tmp_cvs[bb_idx[0]].reset();
	    }
	  else
	    {
	      knots.clear();
	      knots.insert(knots.begin(), tmp_cvs[bb_idx[0]]->basis().begin(),
			   tmp_cvs[bb_idx[0]]->basis().end());
	      for (kj=0; kj<init_cvs.size(); ++kj)
		if (init_cvs[kj].get() == cvs_in[bb_idx[0]].first.get() ||
		    init_cvs[kj].get() == cvs_in[bb_idx[0]].second.get())
		  {
		    init_cvs.erase(init_cvs.begin()+kj);
		    break;
		  }
	      for (kj=1; kj<bb_idx.size(); ++kj)
		tmp_cvs[bb_idx[kj]].reset();
	    }
	}
      else
	for (size_t kj=0; kj<bb_idx.size(); ++kj)
	  {
	    size_t kh;
	    for (kh=0; kh<init_cvs.size(); ++kh)
	      if (init_cvs[kh].get() == cvs_in[bb_idx[kj]].first.get() ||
		  init_cvs[kh].get() == cvs_in[bb_idx[kj]].second.get())
		{
		  init_cvs.erase(init_cvs.begin()+kh);
		  break;
		}
	  }
    
      
      BsplineBasis init_basis((int)knots.size()-order, order, knots.begin());
	  
      app_cvs = AdaptSurface::curveApprox(&init_cvs[0], 
					  (int)init_cvs.size(),
					  init_basis, tol);
    }
  else
    app_cvs = AdaptSurface::curveApprox(&init_cvs[0], 
					(int)init_cvs.size(), tol);
  // Collect final curves
  int kj;
  for (ki=0, kj=0; ki<nmb_cvs; ++ki)
    {
      if (!tmp_cvs[ki].get())
	cvs_out.push_back(app_cvs[kj++]);
      else
	cvs_out.push_back(tmp_cvs[ki]);
    }
}

//---------------------------------------------------------------------------
void 
ftSurface::getBoundaryCurves(double kink,
			     vector<pair<shared_ptr<ParamCurve>,
					 shared_ptr<ParamCurve> > >& cvs, 
				      bool only_corner)
//---------------------------------------------------------------------------
{
  double eps = degenerate_eps_;

  // Fetch outer edge loop
  size_t nmb_edges = boundary_loops_[0]->size();
  vector<ftEdge*> edges(nmb_edges);
  for (size_t ki=0; ki<nmb_edges; ++ ki)
    edges[ki] = (boundary_loops_[0]->getEdge(ki))->geomEdge();

  // Organize the edge loop such that it starts with a corner or a 
  // change in the adjacent surface
  ftEdge *curr = edges[0];
  ftEdge *prev = edges[edges.size()-1];
  ftFaceBase *twin1 = NULL, *twin2 = NULL;
  if (curr->twin())
    twin1 = curr->twin()->geomEdge()->face();
  Point tan1 = curr->tangent(curr->tMin());
  while (prev != curr)
    {
      // Check continuity
      Point tan2 = prev->tangent(prev->tMax());
      if (tan1.angle(tan2) > kink)
	break; // A corner is found

      if (!only_corner)
	{
	  // Check adjacent face
	  if (prev->twin())
	    twin2 = prev->twin()->geomEdge()->face();
	  if (twin1 != twin2)
	    break;  // Different adjacent faces
	}
      
      // Reorganize
      edges.insert(edges.begin(), prev);
      edges.pop_back();
      twin1 = twin2;
      tan1 = prev->tangent(prev->tMin());
      prev = edges[edges.size()-1];
    }

  // Collect curve information for each edge
  size_t idx;
  bool turned = false;
  shared_ptr<ParamCurve> dummy;
  for (idx=0; idx<edges.size(); )
    {
      vector<double> parmin1, parmin2, parmax1, parmax2;
      vector<shared_ptr<ParamCurve> > crvs1, crvs2;
      ftEdge *e1 = edges[idx];
      ftEdge *e2 = (e1->twin()) ? e1->twin()->geomEdge() : NULL;
      shared_ptr<ParamCurve> cv1 = e1->geomCurve();
      shared_ptr<ParamCurve> cv2;
      crvs1.push_back(cv1);
      parmin1.push_back(e1->tMin());
      parmax1.push_back(e1->tMax());
      if (e2)
	{
	  cv2 = e2->geomCurve();
	  crvs2.push_back(cv2);
	  parmin2.push_back(e2->tMin());
	  parmax2.push_back(e2->tMax());
	}

      while (true)
	{
	  ftEdge *e3=NULL, *e4=NULL;
	  e3 = e1->next()->geomEdge();
	  if (e3 == edges[0])
	    {
	      idx = edges.size();
	      break;  // The loop is finshed
	    }

	  e4 = (e3->twin()) ? e3->twin()->geomEdge() : NULL;
	  if (!only_corner)
	    {
	      if ((e2 && e4 && e2->face() != e4->face()) ||
		  (e2!=NULL && e4==NULL) || (e4!=NULL && e2==NULL))
		{
		  idx++;
		  break;  // A new adjacent surface is found
		}
	    }

	  // if (e2 == NULL && e4 == NULL)
	  //   {
	      // Check if we have reached a corner
	      Point tan3 = e1->tangent(e1->tMax());
	      Point tan4 = e3->tangent(e3->tMin());
	      if (tan3.angle(tan4) > kink)
		{
		  idx++;
		  break;
		}
	    // }

	  // Check if a new curve is found
	  shared_ptr<ParamCurve> cv3 = e3->geomCurve();
	  shared_ptr<ParamCurve> cv4;
	  if (cv3.get() != cv1.get() || fabs(e3->tMin()-e1->tMax()) > eps)
	    {
	      crvs1.push_back(cv3);
	      parmin1.push_back(e3->tMin());
	      parmax1.push_back(e3->tMax());
	    }
	  else
	    {
	      size_t last = parmin1.size() - 1;
	      parmin1[last] = std::min(parmin1[last], e3->tMin());
	      parmax1[last] = std::max(parmax1[last], e3->tMax());
	    }

	  if (e4)
	    {
	      cv4 = e4->geomCurve();
	      if (cv4.get() != cv2.get() || !e2 ||
		  fabs(e4->tMin()-e2->tMax()) > eps)
		{
		  crvs2.push_back(cv4);
		  parmin2.push_back(e4->tMin());
		  parmax2.push_back(e4->tMax());
		}
	      else
		{
		  size_t last = parmin2.size() - 1;
		  parmin2[last] = std::min(parmin2[last], e4->tMin());
		  parmax2[last] = std::max(parmax2[last], e4->tMax());
		}
	    }

	  // Update
	  e1 = e3;
	  e2 = e4;
	  cv1 = cv3;
	  cv2 = cv4;
	  idx++;
	}

      // Collect information for the current boundary curve
      int use_curve = 0;
      if (crvs1.size() == 1)
	{
	  cv1 = shared_ptr<ParamCurve>(crvs1[0]->subCurve(parmin1[0], parmax1[0]));
	  if (crvs2.size() > 1 || crvs2.size() == 0)
	    use_curve = 1;
	}
      if (crvs2.size() == 1)
	{
	  cv2 = shared_ptr<ParamCurve>(crvs2[0]->subCurve(parmin2[0], parmax2[0]));
	  if (crvs1.size() > 1 || crvs2.size() == 0)
	    use_curve = 2;
	}
      
      if (use_curve == 0 && crvs1.size() > 1)
	{
	  cv1 = shared_ptr<ParamCurve>(crvs1[0]->subCurve(parmin1[0], parmax1[0]));
	  double dist;
	  for (size_t kr=1; kr<crvs1.size(); ++kr)
	    {
	      Point pt1 = cv1->point(cv1->startparam());
	      Point pt2 = cv1->point(cv1->endparam());
	      shared_ptr<ParamCurve> tmp = 
		shared_ptr<ParamCurve>(crvs1[kr]->subCurve(parmin1[kr], parmax1[kr]));
	      Point pt3 = tmp->point(tmp->startparam());
	      Point pt4 = tmp->point(tmp->endparam());
	      if (std::min(pt2.dist(pt3), pt2.dist(pt4)) > 
		  std::min(pt1.dist(pt3), pt1.dist(pt4)))
		{
		  cv1->reverseParameterDirection();
		  std::swap(pt1, pt2);
		}
	      if (pt2.dist(pt4) < pt2.dist(pt3))
		tmp->reverseParameterDirection();
	  
	      // Perform reparametrization and join with C0 continuity
	      vector<Point> pts1(2);
	      cv1->point(pts1, cv1->endparam(), 1);
	      vector<Point> pts2(2);
	      tmp->point(pts2, tmp->startparam(), 1);
	      double fac = pts2[1].length()/pts1[1].length();
	      
	      double len1 = cv1->estimatedCurveLength();
	      double len2 = tmp->estimatedCurveLength();
	      // TESTING
	      //fac = 1;
	      // END TESTING
	      
	      double s1 = cv1->startparam();
	      double s2 = cv1->endparam();
	      double t1 = tmp->startparam();
	      double t2 = tmp->endparam();

	      fac = len2*(s2-s1)/(len1*(t2-t1));

	      tmp->setParameterInterval(t1, t1+fac*(t2-t1));
#ifdef DEBUG
	      std::ofstream of("pre_append.g2");
	      shared_ptr<CurveOnSurface> bcv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
	      shared_ptr<ParamCurve> dcv = (bcv.get()) ? bcv->spaceCurve() : cv1;
	      dcv->writeStandardHeader(of);
		dcv->write(of);
	       bcv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      dcv = (bcv.get()) ? bcv->spaceCurve() : tmp;
	      dcv->writeStandardHeader(of),
		dcv->write(of);
#endif
	      cv1->appendCurve(tmp.get(), 0, dist, false);
#ifdef DEBUG
	      std::ofstream of2("post_append.g2");
	      bcv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
	      dcv = (bcv.get()) ? bcv->spaceCurve() : cv1;
	      dcv->writeStandardHeader(of2);
	      dcv->write(of2);
#endif
	    }
	  if (crvs2.size() == 0)
	    use_curve = 1;
	}

      if (use_curve == 0 && crvs2.size() > 1)
	{
	  cv2 = shared_ptr<ParamCurve>(crvs2[0]->subCurve(parmin2[0], parmax2[0]));
	  double dist;
	  for (size_t kr=1; kr<crvs2.size(); ++kr)
	    {
	      Point pt1 = cv2->point(cv2->startparam());
	      Point pt2 = cv2->point(cv2->endparam());
	      shared_ptr<ParamCurve> tmp = 
		shared_ptr<ParamCurve>(crvs2[kr]->subCurve(parmin2[kr], parmax2[kr]));
	      Point pt3 = tmp->point(tmp->startparam());
	      Point pt4 = tmp->point(tmp->endparam());
	      if (std::min(pt2.dist(pt3), pt2.dist(pt4)) > 
		  std::min(pt1.dist(pt3), pt1.dist(pt4)))
		{
		  cv2->reverseParameterDirection();
		  std::swap(pt1, pt2);
		}
	      if (pt2.dist(pt4) < pt2.dist(pt3))
		tmp->reverseParameterDirection();
	      //cv2->appendCurve(tmp.get(), 0, dist, false);
	      // try {
	      //   cv2->appendCurve(tmp.get(), 1, dist, true);
	      // }
	      // catch (...)
	      //   {
	      // Perform reparametrization and join with C0 continuity
	      vector<Point> pts1(2);
	      cv2->point(pts1, cv2->endparam(), 1);
	      vector<Point> pts2(2);
	      tmp->point(pts2, tmp->startparam(), 1);
	      double fac = pts2[1].length()/pts1[1].length();

	      double len1 = cv2->estimatedCurveLength();
	      double len2 = tmp->estimatedCurveLength();
	      
	      double s1 = cv2->startparam();
	      double s2 = cv2->endparam();
	      double t1 = tmp->startparam();
	      double t2 = tmp->endparam();
	      fac = len2*(s2-s1)/(len1*(t2-t1));

	      tmp->setParameterInterval(t1, t1+fac*(t2-t1));
#ifdef DEBUG
	      std::ofstream of("pre_append.g2");
	      shared_ptr<CurveOnSurface> bcv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
	      shared_ptr<ParamCurve> dcv = (bcv.get()) ? bcv->spaceCurve() : cv2;
	      dcv->writeStandardHeader(of);
		dcv->write(of);
	       bcv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	      dcv = (bcv.get()) ? bcv->spaceCurve() : tmp;
	      dcv->writeStandardHeader(of),
		dcv->write(of);
#endif
	      cv2->appendCurve(tmp.get(), 0, dist, false);
#ifdef DEBUG
	      std::ofstream of2("post_append.g2");
	      bcv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
	      dcv = (bcv.get()) ? bcv->spaceCurve() : cv2;
	      dcv->writeStandardHeader(of2);
	      dcv->write(of2);
#endif
	      // }
	    } 
	  if (crvs1.size() == 0)
	    use_curve = 2;
	}
	  
      // Check orientation
      bool opposite = false;
      if (cv2.get())
	{
	  Point pt10 = cv1->point(cv1->startparam());
	  Point pt20 = cv1->point(cv1->endparam());
	  Point pt30 = cv2->point(cv2->startparam());
	  Point pt40 = cv2->point(cv2->endparam());
	  if (pt10.dist(pt40) + pt20.dist(pt30) < pt10.dist(pt30) + pt20.dist(pt40))
	    opposite = true;
	}

// #ifdef DEBUG
//       shared_ptr<CurveOnSurface> sfcv1 = 
// 	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv1);
//       shared_ptr<CurveOnSurface> sfcv2 = 
// 	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv2);
//       if (sfcv1.get() && sfcv2.get()) {
// 	sfcv1->spaceCurve()->writeStandardHeader(of);
// 	sfcv1->spaceCurve()->write(of);
// 	sfcv2->spaceCurve()->writeStandardHeader(of);
// 	sfcv2->spaceCurve()->write(of);
//       }
// #endif

      // Store curves. Mind the orientation
      if (use_curve != 2 && turned)
	cv1->reverseParameterDirection();
      if ((turned && !opposite) ||
	  (!turned && opposite))
	{
	  if (use_curve != 1)
	    cv2->reverseParameterDirection();
	  turned = true;
	}
      else
	turned = false;
      
      if (use_curve == 0)
	cvs.push_back(make_pair(cv1, cv2));
      else if (use_curve == 1)
	cvs.push_back(make_pair(cv1, dummy));
      else
	cvs.push_back(make_pair(dummy, cv2));
      
    }
  
}

//---------------------------------------------------------------------------
void ftSurface::closestPoint(const Point& pt,
			     double&  clo_u,
			     double&  clo_v, 
			     Point& clo_pt,
			     double&  clo_dist,
			     double   epsilon) const
//---------------------------------------------------------------------------
{
    surf_->setIterator(Iterator_geometric);
    surf_->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
}

//---------------------------------------------------------------------------
void ftSurface::closestPoint(const Point& pt,
			     double&  clo_u,
			     double&  clo_v, 
			     Point& clo_pt,
			     double&  clo_dist,
			     double   epsilon,
			     const RectDomain* domain_of_interest,
			     double   *seed) const
//---------------------------------------------------------------------------
{
    surf_->setIterator(Iterator_geometric);
    surf_->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
			domain_of_interest, seed);
}

//---------------------------------------------------------------------------
ftEdgeBase* ftSurface::closestBoundaryPoint(const Point& pt,
					    const Point& in_vec,
					    double&  clo_u,
					    double&  clo_v, 
					    Point& clo_pt,
					    double&  clo_dist,
					    double& clo_par) const
//---------------------------------------------------------------------------
{
    ftEdgeBase* e;
    clo_dist = 1e05; // @@sbr Should be large enough.
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
      vector<shared_ptr<ftEdgeBase> > edges = boundary_loops_[ki]->getEdges();
      for (size_t kj=0; kj<edges.size(); ++kj)
	{
	  double clo_t, clo_dist2;
	  Point clo_pt2;
	  edges[kj]->closestPoint(pt, clo_t, clo_pt2, clo_dist2);
	if (in_vec.dimension() == pt.dimension() && 
	    (clo_pt2-pt)*in_vec < 0.0)
	  continue;

	  if (clo_dist2 < clo_dist) {
	    clo_dist = clo_dist2;
	    e = edges[kj].get();
	    clo_par = clo_t;
	    clo_pt = clo_pt2;
	  }
	}
    }
    Point face_par = e->geomEdge()->faceParameter(clo_par);
    clo_u = face_par[0];
    clo_v = face_par[1];

    return e;
}

//---------------------------------------------------------------------------
ftEdgeBase* ftSurface::closestOuterBoundaryPoint(const Point& pt,
						 const Point& in_vec,
						 double&  clo_u,
						 double&  clo_v, 
						 Point& clo_pt,
						 double&  clo_dist,
						 double& clo_par) const
//---------------------------------------------------------------------------
{
    ftEdgeBase* e;
    clo_dist = 1e05; // @@sbr Should be large enough.
    vector<shared_ptr<ftEdgeBase> > edges = boundary_loops_[0]->getEdges();
    for (size_t kj=0; kj<edges.size(); ++kj)
      {
	double clo_t, clo_dist2;
	Point clo_pt2;
	edges[kj]->closestPoint(pt, clo_t, clo_pt2, clo_dist2);
	if (in_vec.dimension() == pt.dimension() && 
	    (clo_pt2-pt)*in_vec < 0.0)
	  continue;

	if (clo_dist2 < clo_dist) {
	  clo_dist = clo_dist2;
	  e = edges[kj].get();
	  clo_par = clo_t;
	  clo_pt = clo_pt2;
	}
      }
    Point face_par = e->geomEdge()->faceParameter(clo_par);
    clo_u = face_par[0];
    clo_v = face_par[1];
    
    return e;
}

// In parameter domain.
//---------------------------------------------------------------------------
ftEdgeBase* ftSurface::edgeClosestToPoint(double u, double v)
//---------------------------------------------------------------------------
{
    ftEdgeBase* e;

    // If a spline surface has a neighbour along only a part of an
    // edge, the edge will be split. Hence we may not assume that a
    // surface has 4 topological edges.

//     if (surf_->instanceType() == Class_SplineSurface) {
// 	// @@ Assumes that the domain is rectangular
// 	const RectDomain& dom
// 	    = dynamic_cast<const RectDomain&>(surf_->parameterDomain());
// 	// The diagonals divide the domain in 4 triangles. Find out
// 	// which we're in.
// 	Vector2D pt = Vector2D(u,v) - 0.5*(dom.lowerLeft() + dom.upperRight());
// 	Vector2D diag1norm(dom.vmin() - dom.vmax(), dom.umax() - dom.umin());
// 	Vector2D diag2norm(dom.vmin() - dom.vmax(), dom.umin() - dom.umax());
// 	bool left1 = pt*diag2norm > 0.0;
// 	bool left2 = pt*diag1norm < 0.0;
// 	int edgenum;// = left1 ? 2 : 0;
// // 	int l2mod = left1 ? 1 : -1;
// // 	edgenum += left2 ? l2mod : -l2mod;
// 	if (left1)
// 	  edgenum = (left2) ? 0 : 3;
// 	else
// 	  edgenum = (left2) ? 1 : 2;

// 	// Get the right edge.
// 	e = first_edges_[0].get();
// 	for (int i = 0; i < edgenum; ++i) {
// 	    e = e->next();
// 	}
//     } else {
	// This operation costs more, but needed for handling of trimmed sfs.
	Point space_pt = point(u, v);
	double global_clo_dist = 1e05; // @@sbr Should be large enough.
	ftEdgeBase* global_clo_edge;
	for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) {
	    ftEdgeBase* first_edge = boundary_loops_[ki]->getEdge(0).get();
	    double clo_t, clo_dist;
	    Point clo_pt;
	    first_edge->closestPoint(space_pt, clo_t, clo_pt, clo_dist);
	    if (clo_dist < global_clo_dist) {
		global_clo_dist = clo_dist;
		global_clo_edge = first_edge;
	    }
	    ftEdgeBase* curr_edge = first_edge->next();
	    while(curr_edge != first_edge) {
		curr_edge->closestPoint(space_pt, clo_t, clo_pt, clo_dist);
		if (clo_dist < global_clo_dist) {
		    global_clo_dist = clo_dist;
		    global_clo_edge = curr_edge;
		}
		curr_edge = curr_edge->next();
	    }
	}
	e = global_clo_edge;
//     }

    return e;
}

//---------------------------------------------------------------------------
void ftSurface::write(std::ostream& os)
//---------------------------------------------------------------------------
{
    try {
	surf_->writeStandardHeader(os);
	surf_->write(os);
    }
    catch(...) {
	std::cerr << "Surface didn't want to be written!" << std::endl;
    }
}

//===========================================================================
ftMessage ftSurface::smoothOutFace(int edge_cont, double approx_orig_tol, 
				   double deg_tol, double& maxerr,
				   double& meanerr, double approx_weight)
//===========================================================================
{
    ftMessage status;

    // We must put surface in a ftSmoothSurf, and make a call to inner smoothing.
    shared_ptr<SplineSurface> spline_sf =
	dynamic_pointer_cast<SplineSurface, ParamSurface>(surf_);
    ALWAYS_ERROR_IF(spline_sf.get() == 0, "Surface is not a spline surface.");

    double approx_tol = -1.0; // We will only perform smoothing wrt original surface.
    vector<int> edge_derivs(4, edge_cont);
    int max_update_iter = 8; // How many iterations will we perform?

    // The current surface must be reparameterized to get a reasonable
    // reliability of the smoothing process.
    double umin1 = spline_sf->startparam_u();
    double vmin1 = spline_sf->startparam_v();
    double umax1 = spline_sf->endparam_u();
    double vmax1 = spline_sf->endparam_v();

    double length_u, length_v;
    GeometryTools::estimateSurfaceSize(*spline_sf, length_u, length_v);
    spline_sf->setParameterDomain(umin1, umin1+length_u,
				  vmin1, vmin1+length_v);

    // Perform smoothing
    ftSmoothSurf smoothsrf(spline_sf, approx_tol, approx_orig_tol, edge_derivs, max_update_iter);
    smoothsrf.setApproxWeight(approx_weight);
    ftPointSet points; // Empty set.
    bool isOK = false;
    try {
	isOK = smoothsrf.update(points, deg_tol, false);
    } catch(...) {
	status.setError(FT_ERROR_IN_SURFACE_CREATION);
	return status;
    }

    // Fetch distance between original and modified surface
    smoothsrf.getError(maxerr, meanerr);

    if (!isOK)
	status.addWarning(FT_NO_SMOOTHING);

    // Return to initial parametrization
    spline_sf->setParameterDomain(umin1, umax1, vmin1, vmax1);
    return status;
}

//===========================================================================
shared_ptr<SplineCurve>  
ftSurface::getBoundaryPiece(Point& pt1, Point& pt2, double eps)
//===========================================================================
{
  shared_ptr<BoundedSurface> bd_surf =
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf_);
  if (bd_surf.get())
    {
      vector<shared_ptr<CurveOnSurface> > cvs1;
      bd_surf->getBoundaryInfo(pt1, pt2, cvs1);

      // We don't know the sequence of the input points so we
      // try both

      vector<shared_ptr<CurveOnSurface> > cvs2;
      bd_surf->getBoundaryInfo(pt2, pt1, cvs2);

      
      if (cvs1.size() > 0 && cvs1.size() <= cvs2.size())
	{
	  shared_ptr<SplineCurve> cv(cvs1[0]->geometryCurve());
	  int cont = 1;
	  double dist;
	  for (size_t ki=1; ki<cvs1.size(); ++ki)
	    cv->appendCurve(cvs1[ki]->geometryCurve(), cont, dist);
  
	  return cv;
	}
      else
	{
 	  shared_ptr<SplineCurve> cv(cvs2[0]->geometryCurve());
	  int cont = 1;
	  double dist;
	  for (size_t ki=1; ki<cvs2.size(); ++ki)
	    cv->appendCurve(cvs2[ki]->geometryCurve(), cont, dist);
	  cv->reverseParameterDirection();
	  return cv;
	}
    }
  else
    {
      SplineCurve *cv, *crosscv;
      surf_->getBoundaryInfo(pt1, pt2, eps, cv, crosscv);
      if (crosscv)
	delete crosscv;
      shared_ptr<SplineCurve> cv2(cv);

      return cv2;
    }
}

 //===========================================================================
bool ftSurface::getSurfaceDisconts(double tol,  vector<double>& disc_u,
				   vector<double>& disc_v)
//===========================================================================
{
  shared_ptr<SplineSurface> surf = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(surf_);

  if (surf.get())
  {
      // Spline surface
      // The last parameters tells that we are not to compute geometric
      // discontinuities
      GeometryTools::surfaceKinks(*surf, tol, disc_u, disc_v, false);
      return (disc_u.size() > 0 || disc_v.size() > 0);
  }
  else
  {
      shared_ptr<BoundedSurface> bd_surf = 
	  dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf_);
      if (bd_surf.get())
      {
	  // Bounded surface
	  BoundedUtils::trimSurfaceKinks(*bd_surf, tol, disc_u, disc_v, false);
	  return (disc_u.size() > 0 || disc_v.size() > 0);
      }
  }

      // Other surfaces will (nov 2008) be elementary and have no kinks
    return false;
}

//===========================================================================
bool ftSurface::getSurfaceKinks(double angtol, vector<double>& g1_disc_u,
				vector<double>& g1_disc_v)
//===========================================================================
{
  shared_ptr<SplineSurface> surf = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(surf_);

  if (surf.get())
  {
      // Spline surface
      GeometryTools::surfaceKinks(*surf, angtol, g1_disc_u, g1_disc_v);
      return (g1_disc_u.size() > 0 || g1_disc_v.size() > 0);
  }
  else
  {
      shared_ptr<BoundedSurface> bd_surf = 
	  dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf_);
      if (bd_surf.get())
      {
	  // Bounded surface
	  BoundedUtils::trimSurfaceKinks(*bd_surf, angtol, g1_disc_u, g1_disc_v);
	  return (g1_disc_u.size() > 0 || g1_disc_v.size() > 0);
      }

      // Other surfaces will (nov 2008) be elementary and have no kinks
      return false;
  }
}

//===========================================================================
vector<shared_ptr<ftSurface> > ftSurface::splitAlongKinks(double angtol)
    // Split a surface along constant parameter line kinks. Currently
    // only spline surface is implemented.
//===========================================================================
{
  vector<shared_ptr<ftSurface> > subsfs;
  
  shared_ptr<SplineSurface> surf = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(surf_);
  if (surf.get() == 0)
    return subsfs;

  vector<double> g1_disc_u;
  vector<double> g1_disc_v;

  GeometryTools::surfaceKinks(*surf, angtol, g1_disc_u, g1_disc_v);
  vector<shared_ptr<SplineSurface> > subspline;
  subspline = GeometryTools::splitInKinks(*surf, g1_disc_u, g1_disc_v);

  for (size_t ki=0; ki<subspline.size(); ki++)
      subsfs.push_back(shared_ptr<ftSurface>(new ftSurface(subspline[ki], (int)ki)));

  return subsfs;
}

//===========================================================================
  ftMessage ftSurface::removeGap(ftEdgeBase* e1, ftEdgeBase* e2, ftFaceBase *other,
				 double epsge)
//===========================================================================
{
    ftMessage status;

    // Check input
    if (e1->twin() != e2)
      {
	status.addWarning(FT_SURFACE_NOT_MODIFIED);
	return status;  // Does not make sense
      }

    // Fetch the geometry entities.
    shared_ptr<ParamSurface> srf1 = surface();
    shared_ptr<ParamSurface> srf2 = other->surface();

    ftEdge* e1g = e1->geomEdge();
    ftEdge* e2g = e2->geomEdge();
    if (!e1g || !e2g)
      {
	status.addWarning(FT_SURFACE_NOT_MODIFIED);
	return status;  // Does not make sense
      }

    shared_ptr<ParamCurve> bdcv1 = e1g->geomCurve();
    shared_ptr<ParamCurve> bdcv2 = e2g->geomCurve();

    Point vertex1 = e1g->getVertex(true)->getVertexPoint();
    Point vertex2 = e1g->getVertex(false)->getVertexPoint();

    double start1, end1, start2, end2;
    start1 = e1->tMin();
    end1 = e1->tMax();
    start2 = e2->tMin();
    end2 = e2->tMax();

    // Case distinction
    shared_ptr<SplineSurface> splsf1 = 
      dynamic_pointer_cast<SplineSurface, ParamSurface>(srf1);
    shared_ptr<SplineSurface> splsf2 = 
      dynamic_pointer_cast<SplineSurface, ParamSurface>(srf2);
    shared_ptr<BoundedSurface> bdsf1 = 
      dynamic_pointer_cast<BoundedSurface, ParamSurface>(srf1);
    shared_ptr<BoundedSurface> bdsf2 = 
      dynamic_pointer_cast<BoundedSurface, ParamSurface>(srf2);

    shared_ptr<CurveOnSurface> sfcv1 = 
      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv1);
    shared_ptr<CurveOnSurface> sfcv2 = 
      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv2);


    if (splsf1 && splsf2 && sfcv1 && sfcv2)
      GapRemoval::removeGapSpline(splsf1, sfcv1, start1, end1,
				  splsf2, sfcv2, start2, end2,
				  vertex1, vertex2, epsge);
    else if (bdsf1 && bdsf2 && sfcv1 && sfcv2)
      GapRemoval::removeGapTrim(sfcv1, start1, end1,
				sfcv2, start2, end2,
				vertex1, vertex2, epsge);
//     else if (splsf1 && bdsf2 && sfcv1 && sfcv2) 
//       GapRemoval::removeGapSplineTrim(splsf1, sfcv1, start1, end1,
// 				      bdsf2, sfcv2, start2, end2,
// 				      vertex1, vertex2, epsge);
//     else if (splsf2 && bdsf1 && sfcv1 && sfcv2) 
//       GapRemoval::removeGapSplineTrim(splsf2, sfcv2, start2, end2,
// 				      bdsf1, sfcv1, start1, end1,
// 				      vertex1, vertex2, epsge);
    else
      {
	status.addWarning(FT_SURFACE_NOT_MODIFIED);
	return status;  // Not implemented
      }

    // Update topology information

    return FT_OK;
}

//===========================================================================
bool ftSurface::isClose(shared_ptr<Vertex> v, double tol) const
//===========================================================================
  {
    Point clo_pt;
    double clo_u, clo_v;
    double clo_dist;

    closestPoint(v->getVertexPoint(),
		 clo_u, clo_v,
		 clo_pt,
		 clo_dist,
		 1.0e-9);
    return (clo_dist <= tol);

  }

//===========================================================================
vector<shared_ptr<Vertex> > ftSurface::vertices() const
//===========================================================================
{
  vector<shared_ptr<Vertex> > result;

  for (size_t i = 0; i < boundary_loops_.size(); ++i)
    {
      shared_ptr<Loop> aLoop = boundary_loops_[i];
      for (size_t j = 0; j < aLoop->size(); ++j)
	{
	  shared_ptr<ftEdgeBase> edge = aLoop->getEdge(j);
	  for (int k = 0; k < 2; ++k)
	    {
	      shared_ptr<Vertex> vert;
	      vert = edge->geomEdge()->getVertex(k == 0);
	      size_t l = 0;
	      for (; l< result.size(); ++l)
		if (result[l].get() == vert.get())
		  break;
	      if (l == result.size())
		result.push_back(vert);
	    }
	}
    }

  return result;
}


//===========================================================================
vector<shared_ptr<Vertex> > ftSurface::getNonCornerVertices(double kink) const
//===========================================================================
{
  vector<shared_ptr<Vertex> > result;

  for (size_t i = 0; i < boundary_loops_.size(); ++i)
    {
      shared_ptr<Loop> aLoop = boundary_loops_[i];

      for (size_t j = 0; j < aLoop->size(); ++j)
	{
	  shared_ptr<ftEdgeBase> edge1 = aLoop->getEdge(j);
	  ftEdgeBase* edge2 = edge1->next();
	    
	  shared_ptr<Vertex> vert = edge1->geomEdge()->getVertex(false);

	  // Check angle
	  Point tan1 = edge1->tangent(edge1->tMax());
	  Point tan2 = edge2->tangent(edge2->tMin());
	  double tang = tan1.angle(tan2);
	  if (tang <= kink)
	    result.push_back(vert);
	}
    }
  return result;
}

//===========================================================================
vector<shared_ptr<Vertex> > ftSurface::getNonCornerVertices(double kink,
							    int loop_idx) const
//===========================================================================
{
  vector<shared_ptr<Vertex> > result;

  shared_ptr<Loop> aLoop = boundary_loops_[loop_idx];

  for (size_t j = 0; j < aLoop->size(); ++j)
    {
      shared_ptr<ftEdgeBase> edge1 = aLoop->getEdge(j);
      ftEdgeBase* edge2 = edge1->next();
	    
      shared_ptr<Vertex> vert = edge1->geomEdge()->getVertex(false);

      // Check angle
      Point tan1 = edge1->tangent(edge1->tMax());
      Point tan2 = edge2->tangent(edge2->tMin());
      double tang = tan1.angle(tan2);
      if (tang <= kink)
	result.push_back(vert);
    }
  return result;
}

//===========================================================================
vector<shared_ptr<Vertex> > ftSurface::getCornerVertices(double kink) const
//===========================================================================
{
  vector<shared_ptr<Vertex> > result;

  for (size_t i = 0; i < boundary_loops_.size(); ++i)
    {
      shared_ptr<Loop> aLoop = boundary_loops_[i];

      for (size_t j = 0; j < aLoop->size(); ++j)
	{
	  shared_ptr<ftEdgeBase> edge1 = aLoop->getEdge(j);
	  ftEdgeBase* edge2 = edge1->next();
	    
	  shared_ptr<Vertex> vert = edge1->geomEdge()->getVertex(false);

	  // Check angle
	  Point tan1 = edge1->tangent(edge1->tMax());
	  Point tan2 = edge2->tangent(edge2->tMin());
	  double tang = tan1.angle(tan2);
	  if (tang > kink)
	    result.push_back(vert);
	}
    }
  return result;
}

//===========================================================================
vector<shared_ptr<Vertex> > ftSurface::getCornerVertices(double kink,
							 int loop_idx) const
//===========================================================================
{
  vector<shared_ptr<Vertex> > result;

  if (loop_idx < 0 || loop_idx >= (int)boundary_loops_.size())
    return result;

  shared_ptr<Loop> aLoop = boundary_loops_[loop_idx];

  for (size_t j = 0; j < aLoop->size(); ++j)
    {
      shared_ptr<ftEdgeBase> edge1 = aLoop->getEdge(j);
      ftEdgeBase* edge2 = edge1->next();
      
      shared_ptr<Vertex> vert = edge1->geomEdge()->getVertex(false);

      // Check angle
      Point tan1 = edge1->tangent(edge1->tMax());
      Point tan2 = edge2->tangent(edge2->tMin());
      double tang = tan1.angle(tan2);
      if (tang > kink)
	result.push_back(vert);
    }

  return result;
}

//===========================================================================
vector<shared_ptr<Vertex> > ftSurface::getCommonVertices(ftSurface *other) const
//===========================================================================
{
  vector<shared_ptr<Vertex> > vx1 = vertices();
  vector<shared_ptr<Vertex> > vx2 = other->vertices();
  vector<shared_ptr<Vertex> > vx3;
  for (size_t ki=0; ki<vx1.size(); ++ki)
    for (size_t kj=0; kj<vx2.size(); ++kj)
      {
	if (vx1[ki].get() == vx2[kj].get())
	  {
	    vx3.push_back(vx1[ki]);
	    break;
	  }
      }
  return vx3;
}

//===========================================================================
vector<shared_ptr<ftEdge> > ftSurface::getCommonEdges(ftSurface *other) const
//===========================================================================
{
    vector<shared_ptr<ftEdge> > edges;
    for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) 
      {
	shared_ptr<Loop> curr_loop = boundary_loops_[ki];
	vector<shared_ptr<ftEdgeBase> > curr_edges = curr_loop->getEdges();
	for (size_t kj=0; kj<curr_edges.size(); ++kj)
	  {
	    if (!curr_edges[kj]->twin())
	      continue;
	    ftEdge* twin = curr_edges[kj]->twin()->geomEdge();
	    if (!twin || twin->face() != other)
	      continue;
	    shared_ptr<ftEdge> curr = 
	      dynamic_pointer_cast<ftEdge,ftEdgeBase>(curr_edges[kj]);
	    if (curr.get())
	      edges.push_back(curr);  // Should always be the case
	  }
      }
    return edges;
}

//===========================================================================
shared_ptr<Vertex> ftSurface::getClosestVertex(const Point& pnt) const
//===========================================================================
{
  shared_ptr<Vertex> dummy;

  // First fetch all vertices
  vector<shared_ptr<Vertex> > vertices = this->vertices();
  if (vertices.size() == 0)
    return dummy;  // No vertices

  size_t min_ind = 0;
  double min_dist = vertices[0]->getVertexPoint().dist(pnt);
  for (size_t ki=1; ki<vertices.size(); ++ki)
    {
      double dist = vertices[ki]->getVertexPoint().dist(pnt);
      if (dist < min_dist)
	{
	  min_ind = ki;
	  min_dist = dist;
	}
    }
  return vertices[min_ind];
}

//===========================================================================
void ftSurface::getBadDistance(vector<pair<ftSurface*, ftEdge*> >& badPairs,
			       double tol) const
//===========================================================================
{
  RectDomain domain = surf_->containingDomain();

  for (size_t i = 0; i < boundary_loops_.size(); ++i)
    boundary_loops_[i]->getBadDistance(badPairs, &domain, tol);
}


//===========================================================================
void ftSurface::getBadDistance(vector<pair<ftEdge*, shared_ptr<Vertex> > >& badPairs,
			       double tol) const
//===========================================================================
{
  for (size_t i = 0; i < boundary_loops_.size(); ++i)
    boundary_loops_[i]->getBadDistance(badPairs, tol);
}



//===========================================================================
void ftSurface::getBadDistance(vector<pair<ftSurface*, shared_ptr<Vertex> > >& badPairs,
			       double tol)
//===========================================================================
{
  vector<shared_ptr<Vertex> > vert_vector;
  vert_vector = vertices();

  for (size_t i = 0; i < vert_vector.size(); ++i)
    if (!isClose(vert_vector[i], tol))
      badPairs.push_back(make_pair(this, vert_vector[i]));
}


//===========================================================================
void ftSurface::getPosTangentSurfaceDiscont(vector<ftEdge*>& badPos,
					    vector<ftEdge*>& badTangent,
					    double tol, double kink, double bend, int leastSurfIndex,
					    shared_ptr<SurfaceModel> sm) const
//===========================================================================
{
  for (size_t i = 0; i < boundary_loops_.size(); ++i)
    boundary_loops_[i]->getPosTangentSurfaceDiscont(badPos, badTangent, tol, kink, bend, leastSurfIndex, sm);
}

//===========================================================================
bool ftSurface::checkLoopOrientation(vector<shared_ptr<Loop> >& inconsistent_loops) const
//===========================================================================
{
    bool loop_OK, is_OK = true;

    // For all loops
    size_t ki;
    for (ki=0; ki<boundary_loops_.size(); ++ki)
    {
	// Check first for consistency in orientation between edges and curves
	loop_OK = boundary_loops_[ki]->checkConsistency();
	if (!loop_OK)
	{
	    inconsistent_loops.push_back(boundary_loops_[ki]);
	    is_OK = false;
	    continue;
	}
    }

    // Check surface type
    if (surf_->instanceType() != Class_BoundedSurface)  
	return (is_OK) ? true : false;   // No more test required for non-trimmed surfaces

    shared_ptr<BoundedSurface> srf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf_);

    // For bounded surfaces, check the consistence between this loop and the
    // curve loop corresponding to the surface

    double eps = boundary_loops_[0]->getTol();
    vector<CurveLoop> crv_loops = srf->allBoundaryLoops(eps); // Get loops as when building this face
    for (size_t kj=0; kj<crv_loops.size(); ++kj)
    {
	int nmb_crvs = crv_loops[kj].size();
	if (nmb_crvs != (int)boundary_loops_[kj]->size())
	    return false;

	vector< shared_ptr<ParamCurve> >::const_iterator crvs = crv_loops[kj].begin();
	for (int ki=0; ki<nmb_crvs; ++ki)
	{
	    shared_ptr<ParamCurve> curr = boundary_loops_[kj]->getEdge(ki)->geomEdge()->geomCurve();
	    loop_OK = (curr.get() == crvs[ki].get());
	    if (!loop_OK)
	    {
		inconsistent_loops.push_back(boundary_loops_[kj]);
		is_OK = false;
		break;
	    }
	}
    }
    if (!is_OK)
	return false;
    
    vector<CurveLoop> bd_loops =
	srf->absolutelyAllBoundaryLoops(); // Return actual loops, including deg cvs.
    double loop_tol = bd_loops[0].getSpaceEpsilon();
    double int_tol = 1.0e-6;
    vector<vector<shared_ptr<CurveOnSurface> > > all_bd_cvs;
    for (ki=0; ki<bd_loops.size(); ++ki)
    {
	bool CCW;
	// Check for loop orientation
	// First represent the loop with surface curves
	vector<shared_ptr<CurveOnSurface> > cvs_on_sf;
	vector<shared_ptr<ParamCurve> > bd_cvs = bd_loops[ki].getCurves();
	try {
	    LoopUtils::representAsSurfaceCurves(bd_cvs, srf, cvs_on_sf);
	}
	catch(...)
	{
	    inconsistent_loops.push_back(boundary_loops_[ki]);
	    is_OK = false;
	    continue;
	}

	all_bd_cvs.push_back(cvs_on_sf);
	loop_tol = std::max(loop_tol, bd_loops[ki].getSpaceEpsilon());

	// Check if the loop has the correct orientation
	try {
	    CCW = LoopUtils::paramIsCCW(cvs_on_sf, int_tol, int_tol);
	    //CCW = LoopUtils::paramIsCCW(cvs_on_sf, loop_tol, int_tol);
	}
	catch(...)
	{
	    inconsistent_loops.push_back(boundary_loops_[ki]);
	    is_OK = false;
	    continue;
	}
	if ((ki==0 && !CCW) || (ki>0 && CCW))
	{
	    inconsistent_loops.push_back(boundary_loops_[ki]);
	    is_OK = false;
	    continue;
	}
    }

    if (!is_OK)
	return false;

    if (bd_loops.size() > 1)
    {
	// Check if the inner loops intersects the outer one
	for (ki=1; ki<bd_loops.size(); ++ki)
	{
	    loop_OK = LoopUtils::firstLoopInsideSecond(all_bd_cvs[ki], all_bd_cvs[0],
						       loop_tol, int_tol);
	    if (!loop_OK)
	    {
		inconsistent_loops.push_back(boundary_loops_[ki]);
		is_OK = false;
		continue;
	    }
	}
    }
    return (is_OK) ? true : false;  
}

//===========================================================================
bool ftSurface::hasAcuteAngle(ftEdge* along_edge, double angtol) const
//===========================================================================
{
    // Get the other face
    ftEdge *twin = along_edge->twin()->geomEdge();
    if (!twin)
	return false;  // No neighbour => no acute angle

    ftSurface* other = twin->face()->asFtSurface();
    if (!other)
	return false;  // Not a ftSurface

    // Estimate number of sample points
    double t1 = along_edge->tMin();
    double t2 = along_edge->tMax();
    double len = along_edge->geomCurve()->estimatedCurveLength(t1, t2);
    double fac = 0.1;
    int nmb_sample = (int)(len/fac);
    nmb_sample = std::min(nmb_sample, 50);
    nmb_sample = std::max(nmb_sample, 4);

    double par;
    double tint = (t2 - t1)/(double)(nmb_sample-1);
    int ki;
    double seed = twin->tMax();
    double del = std::min(0.001*len, 0.1);
    double eps = 1.0e-6;
    for (ki=0, par=t1; ki<nmb_sample; ++ki, par+=tint)
    {
	Point pnt1 = along_edge->point(par);
	Point tan1 = along_edge->tangent(par);
	Point norm1 = along_edge->normal(par);
	Point dir1 = norm1.cross(tan1);
	dir1.normalize();
	Point pos1 = pnt1 + del*dir1;

	double u1, v1, d1;
	Point clo1;
	closestPoint(pos1, u1, v1, clo1, d1, eps);

	Point pnt2;
	double clo_t, clo_dist;
	twin->closestPoint(pnt1, clo_t, pnt2, clo_dist, &seed);
	seed = clo_t;

	Point tan2 = twin->tangent(clo_t);
	Point norm2 = twin->normal(clo_t);
	Point dir2 = norm2.cross(tan2);
	dir2.normalize();
	Point pos2 = pnt2 + del*dir2;

	double u2, v2, d2;
	Point clo2;
	other->closestPoint(pos2, u2, v2, clo2, d2, eps);

	Point vec1 = pnt1 - clo1;
	Point vec2 = clo2 - pnt2;
	double ang = vec1.angle(vec2);
	if (fabs(M_PI - ang) < angtol)
	    return true;
    }

    return false;
}


//===========================================================================
//  Get pairs of close boundary loop points
void ftSurface::getNarrowRegion(double gap_tol, double tol, 
				vector<pair<shared_ptr<PointOnEdge>, 
				shared_ptr<PointOnEdge> > >& narrow_pt) 
//===========================================================================
{
    double fac = 1.0e-4;  // In removal of trivial intersections
    double p_eps = 1.0e-8;  // To remove trivial SurfaceCurves

    // Get all edges belonging to the boundary loops
    vector<shared_ptr<ftEdgeBase> > edges = createInitialEdges(gap_tol);

    // Compute closest point on edges
    size_t ki, kj;
    for (ki=0; ki<edges.size(); ++ki)
    {
	ftEdge* e1 = edges[ki]->geomEdge();
	if (!e1)
	    continue;  // Unable to do closest point
	shared_ptr<ParamCurve> crv1 = e1->geomCurve();
	double t1_1 = e1->tMin();
	double t1_2 = e1->tMax();
	double del1 = fac*(t1_2 - t1_1);
	for (kj=ki+1; kj<edges.size(); ++kj)
	{
	    ftEdge *e2 = edges[kj]->geomEdge();
	    if (!e2)
		continue;
	    shared_ptr<ParamCurve> crv2 = e2->geomCurve();
	    double t2_1 = e2->tMin();
	    double t2_2 = e2->tMax();
	    double del2 = fac*(t2_2 - t2_1);

	    double par1, par2, dist;
	    double seed1 = 0.5*(t1_1 + t1_2); // For the time being
	    double seed2 = 0.5*(t2_1 + t2_2); // For the time being
	    Point pnt1, pnt2;
	    ClosestPoint::closestPtCurves(crv1.get(), crv2.get(), t1_1, t1_2,
			    t2_1, t2_2, seed1, seed2, par1, par2, 
			    dist, pnt1, pnt2);

	    if (dist <= tol)
	    {
		// Check for the trivial intersection between adjacent edges
		if (e1->next() == e2 && fabs(par1-t1_2) < del1 && fabs(par2-t2_1) < del2)
		    continue;  // Try next edges
		// Check for the trivial intersection between adjacent edges
		if (e2->next() == e1 && fabs(par1-t1_1) < del1 && fabs(par2-t2_2) < del2)
		    continue;  // Try next edges

		// Make an estimate of the curve length between the boundary points
		// First make CurveOnSurfaceCurve
		// Get surface parameter corresponding to the edge parameter
		Point facepar1 = e1->faceParameter(par1);
		Point facepar2 = e2->faceParameter(par2);
		if (facepar1.dist(facepar2) < p_eps)
		    continue;  // Same point in face

		shared_ptr<SplineCurve> par_crv = shared_ptr<SplineCurve>(new SplineCurve(facepar1,
											  facepar2));
		shared_ptr<CurveOnSurface> sf_crv = 
		    shared_ptr<CurveOnSurface>(new CurveOnSurface(surf_, par_crv, true));

		double len = sf_crv->estimatedCurveLength();
		if (len <= tol)
		{
		    shared_ptr<PointOnEdge> pt_e1 = 
			shared_ptr<PointOnEdge>(new PointOnEdge(e1, par1));
		    shared_ptr<PointOnEdge> pt_e2 = 
			shared_ptr<PointOnEdge>(new PointOnEdge(e2, par2));
		    narrow_pt.push_back(make_pair(pt_e1, pt_e2));
		}
		
	    }
	}
    }
}

//===========================================================================
//  Check and fix orientation of boundary loops of trimmed surfaces
bool ftSurface::checkAndFixBoundaries()
//===========================================================================
{
  shared_ptr<BoundedSurface> bd_surf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf_);
  if (bd_surf.get())
    {
      // A trimmed surface is found
      int fix = BoundedUtils::checkAndFixLoopOrientation(bd_surf);
      if (fix == 2)
	{
	  // The curve loop is changed in this step. Must update
	  // edge loop
	  clearInitialEdges();
	  (void)createInitialEdges(degenerate_eps_, kink_);
	}
      return (fix) ? true : false;
    }
  else 
    return false;
}


//===========================================================================
// 
double ftSurface::area(double tol) const
//===========================================================================
{
  double tol2d = 1.0e-4; // For the time being

  shared_ptr<BoundedSurface> bd_surf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf_);
  if (bd_surf.get())
    {
      // A trimmed surface is found
      // Get underlying surface 
      shared_ptr<ParamSurface> sf = bd_surf->underlyingSurface();
      shared_ptr<SplineSurface> sf2 = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);

      // Check if the surface is trimmed along constant parameter
      // curves
      bool iso_trim = bd_surf->isIsoTrimmed(tol2d);
      iso_trim = false;  // For test formaal
      if (iso_trim)
	{
	  // Iso trimmed surface. Get right size surface
	  RectDomain domain = bd_surf->containingDomain();

	  RectDomain dom2 = sf->containingDomain();  // To avoid problems due to numerics
	  double umin = std::max(domain.umin(), dom2.umin());
	  double umax = std::min(domain.umax(), dom2.umax());
	  double vmin = std::max(domain.vmin(), dom2.vmin());
	  double vmax = std::min(domain.vmax(), dom2.vmax());
	  vector<shared_ptr<ParamSurface> > sfs = 
	    sf->subSurfaces(umin, vmin, umax, vmax);

	  sf = sfs[0];
	}

      if (iso_trim && sf2.get())
	{
	  // Iso trimmed spline surface. Compute area of rectangular
	  // domain inside trimming
	  return sf->area(tol);
	}
      
      // Triangulate surface. First estimate size of underlying surface
      double len_u, len_v;
      GeometryTools::estimateSurfaceSize(*sf, len_u, len_v);

      // Define mesh size
      // int total_nmb = 400;  // This number should somehow be dependent
      //                       // on the tolerance
      double fac = len_u/len_v;
      int nmb1 = std::max(2, (int)(fac*20.0));
      int nmb2 = std::max(2, (int)(20.0/fac));
      
      // Make mesh
      shared_ptr<GeneralMesh> mesh;
      if (iso_trim)
	{
	  RectangularSurfaceTesselator tesselator(*sf.get(), nmb1, nmb2, false);
	  tesselator.tesselate();
	  mesh = tesselator.getMesh();

	}
      else
	{
	  ParametricSurfaceTesselator tesselator(*bd_surf.get());
	  //tesselator.changeRes(nmb2, nmb1);
	  tesselator.tesselate();
	  mesh = tesselator.getMesh();
	}

      // Compute area as the total size of all triangles
      double *nodes = mesh->vertexArray();
      int num_triang =  mesh->numTriangles();
      int num_triang_nodes = 3*num_triang;
      unsigned int *tri = mesh->triangleIndexArray();

      std::ofstream out("triangles.g2");

      double area = 0;
      for (int ki=0; ki<num_triang_nodes; ki+=3)
	{
	  Point pnt1 = Point(&nodes[3*tri[ki]],&nodes[3*(tri[ki]+1)]);
	  Point pnt2 = Point(&nodes[3*tri[ki+1]],&nodes[3*(tri[ki+1]+1)]);
	  Point pnt3 = Point(&nodes[3*tri[ki+2]],&nodes[3*(tri[ki+2]+1)]);

	std::vector<double> triang;
	triang.insert(triang.end(), pnt1.begin(), pnt1.end());
	triang.insert(triang.end(), pnt2.begin(), pnt2.end());
	triang.insert(triang.end(), pnt2.begin(), pnt2.end());
	triang.insert(triang.end(), pnt3.begin(), pnt3.end());
	triang.insert(triang.end(), pnt3.begin(), pnt3.end());
	triang.insert(triang.end(), pnt1.begin(), pnt1.end());
	  
	LineCloud polyline(&triang[0], 3);
	polyline.writeStandardHeader(out);
	polyline.write(out);

	  Point vec1 = pnt1 - pnt2;
	  Point vec2 = pnt3 - pnt2;
	  double l1 = vec1.length();
	  double l2 = vec2.length();
	  double ang = vec1.angle(vec2);

	  double curr_area = 0.5*sin(ang)*l1*l2;
	  area += curr_area;
	}

      return area;
      
    }
  else
    return surf_->area(tol);
}

//===========================================================================
// 
// Get neighbouring faces
    void ftSurface::getAdjacentFaces(std::vector<ftSurface*>& neighbours) const
//===========================================================================
{
    std::set<ftSurface*> all_faces;
    for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
	// Get edges
	vector<shared_ptr<ftEdgeBase> > edges = boundary_loops_[ki]->getEdges();

	// Get faces associated to twins
	for (size_t kj=0; kj<edges.size(); ++kj)
	{
	    if (!edges[kj]->twin())
		continue;

	    ftSurface* face = edges[kj]->twin()->face()->asFtSurface();
	    if (face)
		all_faces.insert(face);
	}
    }
    neighbours.insert(neighbours.end(), all_faces.begin(), all_faces.end());
}

//===========================================================================
// 
// 
bool ftSurface::commonSplineSpace(ftSurface *other, double tol) 
//===========================================================================
{
  if (!isSpline() || !other->isSpline())
    return false;

  shared_ptr<ftEdge> edge1, edge2;
  bool neighbours = areNeighbours(other, edge1, edge2);
  if (!neighbours)
    return true;  // Not applicable

  // A common spline space requires corner-to-corner configuration
  if (!isCornerToCorner(other, tol))
    return false;

  // Fetch geometry, check surface type and edge type
  shared_ptr<ParamSurface> srf1 = surface();
  shared_ptr<ParamSurface> srf2 = other->surface();
  shared_ptr<ParamCurve> bdcv1 = edge1->geomCurve();
  shared_ptr<ParamCurve> bdcv2 = edge2->geomCurve();
  shared_ptr<SplineSurface> splsf1 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf1);
  shared_ptr<SplineSurface> splsf2 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf2);
  
  if (!splsf1.get() || !splsf2.get())
    {
      MESSAGE("Check data structure. Inconsistency of surface types");
      return false;  // This should not occur. Inconsistency
    }

  // Check that either none or both surfaces are rational
  if ((splsf1->rational() && !splsf2->rational()) ||
      (!splsf1->rational() && splsf2->rational()))
    return false;

  AdjacencyInfo adj_info = getAdjacencyInfo(other, DEFAULT_SPACE_EPSILON);
  if (!adj_info.adjacency_found_)
    return false;
  bool same_orient = adj_info.same_orient_;
  int bd1 = adj_info.bd_idx_1_; // Index of common boundary curve, 
  // 0=umin, 1=umax, 2=vmin, 3=vmax  
  int bd2 = adj_info.bd_idx_2_;

  // Fetch spline spaces
  BsplineBasis basis1 = (bd1 == 0 || bd1 == 1) ? splsf1->basis_v() : 
    splsf1->basis_u();
  BsplineBasis basis2 = (bd2 == 0 || bd2 == 1) ? splsf2->basis_v() : 
    splsf2->basis_u();

  // Make copy of basis 2
  BsplineBasis basis2_2 = basis2;
  
  // Transform basis 2 to match basis 1 with regard to orientation and domain
  if (!same_orient)
    basis2_2.reverseParameterDirection();
  basis2_2.rescale(basis1.startparam(), basis1.endparam());

  // Check equality of spline space
  bool same = basis1.sameSplineSpace(basis2_2, tol);
  return same;
}

//===========================================================================
// 
// 
void ftSurface::makeCommonSplineSpace(ftSurface *other) 
//===========================================================================
{
  shared_ptr<ftEdge> edge1, edge2;
  bool neighbours = areNeighbours(other, edge1, edge2);
  if (!neighbours)
    return;

  // Fetch geometry, check edge type
  shared_ptr<ParamSurface> srf1 = surface();
  shared_ptr<ParamSurface> srf2 = other->surface();
  shared_ptr<ParamCurve> bdcv1 = edge1->geomCurve();
  shared_ptr<ParamCurve> bdcv2 = edge2->geomCurve();
  shared_ptr<CurveOnSurface> sfcv1 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv1);
  shared_ptr<CurveOnSurface> sfcv2 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv2);
  shared_ptr<SplineSurface> splsf1 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf1);
  shared_ptr<SplineSurface> splsf2 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf2);
  
  if (!splsf1.get() || !splsf2.get())
    return;  // Not spline surfaces
  
  if (!sfcv1.get() || !sfcv2.get())
    {
      MESSAGE("Check data structure. Expecting curve on surface");
      return;
    }
  
  double start1, end1, start2, end2;
  start1 = edge1->tMin();
  end1 = edge1->tMax();
  start2 = edge2->tMin();
  end2 = edge2->tMax();
  Point vertex1 = edge1->getVertex(true)->getVertexPoint();
  Point vertex2 = edge1->getVertex(false)->getVertexPoint();

  AdjacencyInfo adj_info = getAdjacencyInfo(other, DEFAULT_SPACE_EPSILON);
  if (adj_info.adjacency_found_)
    {
      bool same_orient = adj_info.same_orient_;
      //int bd1 = adj_info.bd_idx_1_; // Index of common boundary curve, 
      // 0=umin, 1=umax, 2=vmin, 3=vmax  
      //int bd2 = adj_info.bd_idx_2_;

      double tol = 1.0e-6;  // Not used
      GapRemoval::removeGapSpline(splsf1, sfcv1, start1, end1, 
				  splsf2, sfcv2, start2, end2,
				  vertex1, vertex2, tol, &same_orient);
    }
}

 //===========================================================================
// 
// 
bool ftSurface::isCornerToCorner(ftSurface *other, double tol, int adj_idx) 
//===========================================================================
{
  shared_ptr<ftEdge> edge1, edge2;
  bool neighbours = areNeighbours(other, edge1, edge2, adj_idx);
  if (!neighbours)
    return true;  // Not applicable

  // Fetch geometry, check edge type
  shared_ptr<ParamSurface> srf1 = surface();
  shared_ptr<ParamSurface> srf2 = other->surface();
  shared_ptr<ParamCurve> bdcv1 = edge1->geomCurve();
  shared_ptr<ParamCurve> bdcv2 = edge2->geomCurve();
  shared_ptr<CurveOnSurface> sfcv1 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv1);
  shared_ptr<CurveOnSurface> sfcv2 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv2);
  
  if (!sfcv1.get() || !sfcv2.get())
    {
      MESSAGE("Check data structure. Expecting curve on surface");
      return false;
    }

  bool at_corners = SurfaceTools::cornerToCornerSfs(srf1, sfcv1, srf2, sfcv2, tol);
  return at_corners;
}

//===========================================================================
// 
void ftSurface::splitAtInternalCorner(ftSurface* other,
				      vector<shared_ptr<ftSurface> >& new_face1,
				      vector<shared_ptr<ftSurface> >& new_face2,
				      double tol)
//===========================================================================
{
  // Fetch adjacency information
  int bd1, bd2;
  bool same_orient;
  bool found = getAdjacencyInfo(other, tol, bd1, bd2, same_orient);
  if (!found)
    return;  // Nothing to do

  shared_ptr<ftEdge> edge1, edge2;
  bool neighbours = areNeighbours(other, edge1, edge2);
  if (!neighbours)
    return;

  // Fetch geometry, check edge and surface type
  shared_ptr<ParamSurface> srf1 = surface();
  shared_ptr<ParamSurface> srf2 = other->surface();
  shared_ptr<ParamCurve> bdcv1 = edge1->geomCurve();
  shared_ptr<ParamCurve> bdcv2 = edge2->geomCurve();
  shared_ptr<CurveOnSurface> sfcv1 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv1);
  shared_ptr<CurveOnSurface> sfcv2 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv2);
  shared_ptr<SplineSurface> splsf1 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf1);
  shared_ptr<SplineSurface> splsf2 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf2);
  
  if (!splsf1.get() || !splsf2.get())
    return;  // Not spline surfaces
  
  
  if (!sfcv1.get() || !sfcv2.get())
    {
      MESSAGE("Check data structure. Expecting curve on surface");
      return;
    }

  // Fetch endparameters in face of common boundary
  Point par1 = sfcv1->faceParameter(edge1->tMin());
  Point par2 = sfcv1->faceParameter(edge1->tMax());
  Point par3 = sfcv2->faceParameter(edge2->tMin());
  Point par4 = sfcv2->faceParameter(edge2->tMax());

  // Check if the first surface must be split
  double ptol = 1.0e-6;  // Not a lasting solution
  vector<double> splitpar1;
  splitpar1.push_back(splsf1->startparam_u());
  if (bd1 == 2 || bd1 == 3)
    {
      double tmin = std::min(par1[0], par2[0]);
      double tmax = std::max(par1[0], par2[0]);
      if (tmin > splitpar1[splitpar1.size()-1]+ptol)
	{
	  splsf1->basis_u().knotIntervalFuzzy(tmin, ptol);
	  splitpar1.push_back(tmin);
	}
      if (tmax < splsf1->endparam_u()-ptol && 
	  tmax > splitpar1[splitpar1.size()-1]+ptol)
	{
	  splsf1->basis_u().knotIntervalFuzzy(tmax, ptol);
	  splitpar1.push_back(tmax);
	}
    }
  splitpar1.push_back(splsf1->endparam_u());
  
  vector<double> splitpar2;
  splitpar2.push_back(splsf1->startparam_v());
  if (bd1 == 0 || bd1 == 1)
    {
      double tmin = std::min(par1[1], par2[1]);
      double tmax = std::max(par1[1], par2[1]);
      if (tmin > splitpar2[splitpar2.size()-1]+ptol)
	{
	  splsf1->basis_v().knotIntervalFuzzy(tmin, ptol);
	  splitpar2.push_back(tmin);
	}
      if (tmax < splsf1->endparam_v()-ptol && 
	  tmax > splitpar2[splitpar2.size()-1]+ptol)
	{
	  splsf1->basis_v().knotIntervalFuzzy(tmax, ptol);
	  splitpar2.push_back(tmax);
	}
    }
  splitpar2.push_back(splsf1->endparam_v());
  
  if (splitpar1.size() > 2 || splitpar2.size() > 2)
    {
      for (size_t ki=1; ki<splitpar1.size(); ++ki)
	for (size_t kj=1; kj<splitpar2.size(); ++kj)
	  {
	    shared_ptr<ParamSurface> newsf =
	      shared_ptr<ParamSurface>(splsf1->subSurface(splitpar1[ki-1],
							   splitpar2[kj-1],
							   splitpar1[ki],
							   splitpar2[kj]));
	    shared_ptr<ftSurface> newface =
	      shared_ptr<ftSurface>(new ftSurface(newsf, -1));
	    new_face1.push_back(newface);
	  }
    }

  // Check if the second surface must be split
  vector<double> splitpar3;
  splitpar3.push_back(splsf2->startparam_u());
  if (bd2 == 2 || bd2 == 3)
    {
      double tmin = std::min(par3[0], par4[0]);
      double tmax = std::max(par3[0], par4[0]);
      if (tmin > splitpar3[splitpar3.size()-1]+ptol)
	{
	  splsf2->basis_u().knotIntervalFuzzy(tmin, ptol);
	  splitpar3.push_back(tmin);
	}
      if (tmax < splsf2->endparam_u()-ptol && 
	  tmax > splitpar3[splitpar3.size()-1]+ptol)
	{
	  splsf2->basis_u().knotIntervalFuzzy(tmax, ptol);
	  splitpar3.push_back(tmax);
	}
    }
  splitpar3.push_back(splsf2->endparam_u());
  
  vector<double> splitpar4;
  splitpar4.push_back(splsf2->startparam_v());
  if (bd2 == 0 || bd2 == 1)
    {
      double tmin = std::min(par3[1], par4[1]);
      double tmax = std::max(par3[1], par4[1]);
      if (tmin > splitpar4[splitpar4.size()-1]+ptol)
	{
	  splsf2->basis_v().knotIntervalFuzzy(tmin, ptol);
	  splitpar4.push_back(tmin);
	}
      if (tmax < splsf2->endparam_v()-ptol && 
	  tmax > splitpar4[splitpar4.size()-1]+ptol)
	{
	  splsf2->basis_v().knotIntervalFuzzy(tmax, ptol);
	  splitpar4.push_back(tmax);
	}
    }
  splitpar4.push_back(splsf2->endparam_v());
  
  if (splitpar3.size() > 2 || splitpar4.size() > 2)
    {
      for (size_t ki=1; ki<splitpar3.size(); ++ki)
	for (size_t kj=1; kj<splitpar4.size(); ++kj)
	  {
	    shared_ptr<ParamSurface> newsf =
	      shared_ptr<ParamSurface>(splsf2->subSurface(splitpar3[ki-1],
							   splitpar4[kj-1],
							   splitpar3[ki],
							   splitpar4[kj]));
	    shared_ptr<ftSurface> newface =
	      shared_ptr<ftSurface>(new ftSurface(newsf, -1));
	    new_face2.push_back(newface);
	  }
    }
}

//===========================================================================
// 
// 
bool ftSurface::getAdjacencyInfo(ftSurface *other, double tol,
				 int& bd1, int& bd2, bool& same_orient)
//===========================================================================
{
  MESSAGE("Using deprecated getAdjacencyInfo() method. Use other instead");

  AdjacencyInfo adj_info = getAdjacencyInfo(other, tol);
  if (adj_info.adjacency_found_)
    {
      bd1 = adj_info.bd_idx_1_;
      bd2 = adj_info.bd_idx_2_;
      same_orient = adj_info.same_orient_;
    }
  return adj_info.adjacency_found_;
}



//===========================================================================
// 
// 
AdjacencyInfo ftSurface::getAdjacencyInfo(ftSurface *other, double tol,
					  int adj_idx, bool test_corner)
//===========================================================================
{
  AdjacencyInfo adj_info;

  shared_ptr<ftEdge> edge1, edge2;
  bool neighbours = areNeighbours(other, edge1, edge2, adj_idx);

  shared_ptr<ParamSurface> srf1 = surface();
  shared_ptr<ParamSurface> srf2 = other->surface();
  shared_ptr<ParamCurve> bdcv1, bdcv2;
  if (!neighbours)
    {
      vector<shared_ptr<Vertex> > vx = getCommonVertices(other);

      if (vx.size() == 0)
	{
	  adj_info.adjacency_found_ = false;  // Not applicable
	  return adj_info;
	}

      // Check if there is a degenerate situation
      size_t kj;
      for (kj=0; kj<vx.size(); ++kj)
	{
	  bool deg_found = checkDegAdjacency(other, vx[kj], tol,
					     bdcv1, bdcv2);
	  if (deg_found)
	    break;
	}
      if (kj == vx.size())
	{
	  // No degenerate case found
	  adj_info.adjacency_found_ = false;  // Not applicable
	  return adj_info;
	}	  
    }
  else
    {

      // Fetch geometry, check edge type
      bdcv1 = edge1->geomCurve();
      bdcv2 = edge2->geomCurve();
    }

  // Neighbourhood established. Fetch adjacency details
  shared_ptr<CurveOnSurface> sfcv1 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv1);
  shared_ptr<CurveOnSurface> sfcv2 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv2);
  
  if (!sfcv1.get() || !sfcv2.get())
    {
      MESSAGE("Check data structure. Expecting curve on surface");
      adj_info.adjacency_found_ = false;
      return adj_info;
    }

  int bd1, bd2;
  bool same_orient;
  adj_info.adjacency_found_ = SurfaceTools::getSfAdjacencyInfo(srf1, sfcv1, srf2, sfcv2, 
						 tol, bd1, bd2, same_orient);
  adj_info.bd_idx_1_ = bd1;
  adj_info.bd_idx_2_ = bd2;
  adj_info.same_orient_ = same_orient;

  if (test_corner && adj_info.adjacency_found_)
    adj_info.corner_failed_ = !SurfaceTools::cornerToCornerSfs(srf1, sfcv1, srf2, sfcv2, tol);

  return adj_info;
}

//===========================================================================
// 
// 
AdjacencyInfo ftSurface::getAdjacencyInfo(ftEdge *edge, ftSurface *other,
					  double tol)
//===========================================================================
{
  AdjacencyInfo adj_info;
  shared_ptr<ParamSurface> srf1 = surface();
  shared_ptr<ParamSurface> srf2 = other->surface();
  ftEdge *edge2 = edge->twin()->geomEdge();
  if (!edge2 || edge2->face() != other)
    {
      adj_info.adjacency_found_ = false;
      return adj_info;
    }
  shared_ptr<ParamCurve> bdcv1 = edge->geomCurve();
  shared_ptr<ParamCurve> bdcv2 = edge2->geomCurve();

  // Fetch adjacency details
  shared_ptr<CurveOnSurface> sfcv1 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv1);
  shared_ptr<CurveOnSurface> sfcv2 = 
    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv2);
  
  if (!sfcv1.get() || !sfcv2.get())
    {
      MESSAGE("Check data structure. Expecting curve on surface");
      adj_info.adjacency_found_ = false;
      return adj_info;
    }

  int bd1, bd2;
  bool same_orient;
  adj_info.adjacency_found_ = SurfaceTools::getSfAdjacencyInfo(srf1, sfcv1, srf2, sfcv2, 
						 tol, bd1, bd2, same_orient);
  adj_info.bd_idx_1_ = bd1;
  adj_info.bd_idx_2_ = bd2;
  adj_info.same_orient_ = same_orient;

  return adj_info;
}

//===========================================================================
// 
// 
bool ftSurface::checkDegAdjacency(ftSurface *other, shared_ptr<Vertex> vx,
				  double tol,
				  shared_ptr<ParamCurve>& bdcv1, 
				  shared_ptr<ParamCurve>& bdcv2)
//===========================================================================
{
  // Check if both surfaces are degenerate
  shared_ptr<ParamSurface> srf1 = surface();
  shared_ptr<ParamSurface> srf2 = other->surface();
  bool deg1, deg2;
  bool b1, b2, r1, r2, t1, t2, l1, l2;
  deg1 = srf1->isDegenerate(b1, r1, t1, l1, tol);
  deg2 = srf2->isDegenerate(b2, r2, t2, l2, tol);
  if (!(deg1 && deg2))
    return false;

  // Check if the degeneracy corresponds with the common vertex
  // Get points at degenerate edges for the first surface
  Point vx_pt = vx->getVertexPoint();
  bool found1 = false;
  RectDomain dom1 = srf1->containingDomain();
  if (b1)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	srf1->constParamCurves(dom1.vmin(), true);

      if (cvs.size() > 0)
	{
	  Point mid = cvs[0]->point(0.5*(cvs[0]->startparam()+
					 cvs[0]->endparam()));
	  if (vx_pt.dist(mid) < tol)
	    {
	      // Boundary found
	      found1 = true;

	      // Make corresponding curve-on-surface
	      if (cvs[0]->instanceType() == Class_CurveOnSurface)
		bdcv1 = cvs[0];
	      else
		{
		  shared_ptr<ParamCurve> pcv = 
		    shared_ptr<SplineCurve>
		    (new SplineCurve(Point(dom1.umin(),dom1.vmin()),
				     cvs[0]->startparam(),
				     Point(dom1.umax(), dom1.vmin()),
				     cvs[0]->endparam()));
		  bdcv1 = shared_ptr<ParamCurve>
		    (new CurveOnSurface(srf1, pcv, cvs[0], false, 3));
		}
	    }
	}
    }

  if (r1 && !found1)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	srf1->constParamCurves(dom1.umax(), false);

      if (cvs.size() > 0)
	{
	  Point mid = cvs[0]->point(0.5*(cvs[0]->startparam()+
					 cvs[0]->endparam()));
	  if (vx_pt.dist(mid) < tol)
	    {
	      // Boundary found
	      found1 = true;

	      // Make corresponding curve-on-surface
	      if (cvs[0]->instanceType() == Class_CurveOnSurface)
		bdcv1 = cvs[0];
	      else
		{
		  shared_ptr<ParamCurve> pcv = 
		    shared_ptr<SplineCurve>
		    (new SplineCurve(Point(dom1.umax(),dom1.vmin()),
				     cvs[0]->startparam(),
				     Point(dom1.umax(), dom1.vmax()),
				     cvs[0]->endparam()));
		  bdcv1 = shared_ptr<ParamCurve>
		    (new CurveOnSurface(srf1, pcv, cvs[0], false, 3));
		}
	    }
	}
    }
      
  if (t1 && !found1)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	srf1->constParamCurves(dom1.vmax(), true);

      if (cvs.size() > 0)
	{
	  Point mid = cvs[0]->point(0.5*(cvs[0]->startparam()+
					 cvs[0]->endparam()));
	  if (vx_pt.dist(mid) < tol)
	    {
	      // Boundary found
	      found1 = true;

	      // Make corresponding curve-on-surface
	      if (cvs[0]->instanceType() == Class_CurveOnSurface)
		bdcv1 = cvs[0];
	      else
		{
		  shared_ptr<ParamCurve> pcv = 
		    shared_ptr<SplineCurve>
		    (new SplineCurve(Point(dom1.umin(),dom1.vmax()),
				     cvs[0]->startparam(),
				     Point(dom1.umax(), dom1.vmax()),
				     cvs[0]->endparam()));
		  bdcv1 = shared_ptr<ParamCurve>
		    (new CurveOnSurface(srf1, pcv, cvs[0], false, 3));
		}
	    }
	}
    }

  if (l1 && !found1)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	srf1->constParamCurves(dom1.umin(), false);

      if (cvs.size() > 0)
	{
	  Point mid = cvs[0]->point(0.5*(cvs[0]->startparam()+
					 cvs[0]->endparam()));
	  if (vx_pt.dist(mid) < tol)
	    {
	      // Boundary found
	      found1 = true;

	      // Make corresponding curve-on-surface
	      if (cvs[0]->instanceType() == Class_CurveOnSurface)
		bdcv1 = cvs[0];
	      else
		{
		  shared_ptr<ParamCurve> pcv = 
		    shared_ptr<SplineCurve>
		    (new SplineCurve(Point(dom1.umin(),dom1.vmin()),
				     cvs[0]->startparam(),
				     Point(dom1.umin(), dom1.vmax()),
				     cvs[0]->endparam()));
		  bdcv1 = shared_ptr<ParamCurve>
		    (new CurveOnSurface(srf1, pcv, cvs[0], false, 3));
		}
	    }
	}
    }
 
  // Get points at degenerate edges for the second surface
  bool found2 = false;
  RectDomain dom2 = srf1->containingDomain();
  if (b2)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	srf2->constParamCurves(dom2.vmin(), true);

      if (cvs.size() > 0)
	{
	  Point mid = cvs[0]->point(0.5*(cvs[0]->startparam()+
					 cvs[0]->endparam()));
	  if (vx_pt.dist(mid) < tol)
	    {
	      // Boundary found
	      found2 = true;

	      // Make corresponding curve-on-surface
	      if (cvs[0]->instanceType() == Class_CurveOnSurface)
		bdcv2 = cvs[0];
	      else
		{
		  shared_ptr<ParamCurve> pcv = 
		    shared_ptr<SplineCurve>
		    (new SplineCurve(Point(dom2.umin(),dom2.vmin()),
				     cvs[0]->startparam(),
				     Point(dom2.umax(), dom2.vmin()),
				     cvs[0]->endparam()));
		  bdcv2 = shared_ptr<ParamCurve>
		    (new CurveOnSurface(srf1, pcv, cvs[0], false, 3));
		}
	    }
	}
    }

  if (r2 && !found2)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	srf2->constParamCurves(dom2.umax(), false);

      if (cvs.size() > 0)
	{
	  Point mid = cvs[0]->point(0.5*(cvs[0]->startparam()+
					 cvs[0]->endparam()));
	  if (vx_pt.dist(mid) < tol)
	    {
	      // Boundary found
	      found2 = true;

	      // Make corresponding curve-on-surface
	      if (cvs[0]->instanceType() == Class_CurveOnSurface)
		bdcv2 = cvs[0];
	      else
		{
		  shared_ptr<ParamCurve> pcv = 
		    shared_ptr<SplineCurve>
		    (new SplineCurve(Point(dom2.umax(),dom2.vmin()),
				     cvs[0]->startparam(),
				     Point(dom2.umax(), dom2.vmax()),
				     cvs[0]->endparam()));
		  bdcv2 = shared_ptr<ParamCurve>
		    (new CurveOnSurface(srf1, pcv, cvs[0], false, 3));
		}
	    }
	}
    }
      
  if (t2 && !found2)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	srf2->constParamCurves(dom2.vmax(), true);

      if (cvs.size() > 0)
	{
	  Point mid = cvs[0]->point(0.5*(cvs[0]->startparam()+
					 cvs[0]->endparam()));
	  if (vx_pt.dist(mid) < tol)
	    {
	      // Boundary found
	      found2 = true;

	      // Make corresponding curve-on-surface
	      if (cvs[0]->instanceType() == Class_CurveOnSurface)
		bdcv2 = cvs[0];
	      else
		{
		  shared_ptr<ParamCurve> pcv = 
		    shared_ptr<SplineCurve>
		    (new SplineCurve(Point(dom2.umin(),dom2.vmax()),
				     cvs[0]->startparam(),
				     Point(dom2.umax(), dom2.vmax()),
				     cvs[0]->endparam()));
		  bdcv2 = shared_ptr<ParamCurve>
		    (new CurveOnSurface(srf1, pcv, cvs[0], false, 3));
		}
	    }
	}
    }

  if (l2 && !found2)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	srf2->constParamCurves(dom2.umin(), false);

      if (cvs.size() > 0)
	{
	  Point mid = cvs[0]->point(0.5*(cvs[0]->startparam()+
					 cvs[0]->endparam()));
	  if (vx_pt.dist(mid) < tol)
	    {
	      // Boundary found
	      found2 = true;

	      // Make corresponding curve-on-surface
	      if (cvs[0]->instanceType() == Class_CurveOnSurface)
		bdcv2 = cvs[0];
	      else
		{
		  shared_ptr<ParamCurve> pcv = 
		    shared_ptr<SplineCurve>
		    (new SplineCurve(Point(dom2.umin(),dom2.vmin()),
				     cvs[0]->startparam(),
				     Point(dom2.umin(), dom2.vmax()),
				     cvs[0]->endparam()));
		  bdcv2 = shared_ptr<ParamCurve>
		    (new CurveOnSurface(srf1, pcv, cvs[0], false, 3));
		}
	    }
	}
    }
  return (found1 && found2);
}

//===========================================================================
// 
// 
bool ftSurface::getCorrCoefEnumeration(ftSurface *other, double tol,
				       vector<pair<int,int> >& enumeration) 
//===========================================================================
{
  // Fetch adjacency information
  int bd1, bd2;
  bool same_orient;
  bool found = getAdjacencyInfo(other, tol, bd1, bd2, same_orient);
  if (!found)
    return false;

  // Check that the surfaces have corresponding spline spaces
  bool common = commonSplineSpace(other, tol);
  if (!common)
    return false; // No corresponding coefficient enumeration

  // Fetch surface geometry
  shared_ptr<ParamSurface> srf1 = surface();
  shared_ptr<ParamSurface> srf2 = other->surface();
  shared_ptr<SplineSurface> splsf1 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf1);
  shared_ptr<SplineSurface> splsf2 = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf2);
  if (!splsf1.get() || !splsf2.get())
    return false;

  bool pairwise = SurfaceTools::getCorrCoefEnum(splsf1, splsf2,
				  bd1, bd2, same_orient,
				  enumeration);
  return pairwise;
}

//===========================================================================
// 
// 
bool ftSurface::getFreeBoundaryInfo(double tol, 
				    vector<int>& free_boundaries) 
//===========================================================================
{
  bool all_at_boundaries = true;
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      // Get edges
      vector<shared_ptr<ftEdgeBase> > edges = boundary_loops_[ki]->getEdges();

      // Get no twin edges
      for (size_t kj=0; kj<edges.size(); ++kj)
	{
	  if (edges[kj]->twin())
	    continue;

	  // Get associated boundary curve and check for surface boundary
	  shared_ptr<ParamCurve> bdcv = edges[kj]->geomEdge()->geomCurve();
	  shared_ptr<CurveOnSurface> sfcv = 
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv);
	  if (sfcv.get())
	    {
	      bool same_orientation;
	      int bd = sfcv->whichBoundary(tol, same_orientation);
	      if (bd < 0)
		all_at_boundaries = false;
	      else
		free_boundaries.push_back(bd);
	    }
	  else
	    all_at_boundaries = false;
	}
    }
  return all_at_boundaries;
}

//===========================================================================
// 
// 
bool ftSurface::getBoundaryCoefEnumeration(int bd, 
					   std::vector<int>& enumeration) 
//===========================================================================
{
  // Fetch surface geometry
  shared_ptr<ParamSurface> srf = surface();
  shared_ptr<SplineSurface> splsf = 
    dynamic_pointer_cast<SplineSurface, ParamSurface>(srf);
  if (!splsf.get() || bd < 0 || bd > 3)
    return false;

  bool found = SurfaceTools::getCoefEnumeration(splsf, bd, enumeration);
  return found;
}

//===========================================================================
vector<shared_ptr<EdgeVertex> > ftSurface::getRadialEdges() const
//===========================================================================
{
  vector<shared_ptr<EdgeVertex> > result;

  set<shared_ptr<EdgeVertex> > tmp_result;  // To get the radial edges
  // only once
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      shared_ptr<Loop> aLoop = boundary_loops_[ki];
      for (size_t kj = 0; kj < aLoop->size(); ++kj)
	{
	  ftEdge* edge = aLoop->getEdge(kj)->geomEdge();
	  if (edge->hasEdgeMultiplicity())
	    {
	      tmp_result.insert(edge->getEdgeMultiplicityInstance());
	    }
	}
    }
  result.insert(result.end(), tmp_result.begin(), tmp_result.end());

  return result;
}

//===========================================================================
bool ftSurface::hasRadialEdges() const
//===========================================================================
{
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      shared_ptr<Loop> aLoop = boundary_loops_[ki];
      if (aLoop->hasRadialEdges())
	return true;
    }
  return false;
}

//===========================================================================
bool ftSurface::hasRealRadialEdges() const
//===========================================================================
{
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      shared_ptr<Loop> aLoop = boundary_loops_[ki];
      size_t nmb = aLoop->size();
      for (size_t kj=0; kj<nmb; ++kj)
	{
	  ftEdge *edge = aLoop->getEdge(kj)->geomEdge();
	  if (edge->hasEdgeMultiplicity() && 
	      edge->getEdgeMultiplicityInstance()->nmbUniqueEdges() > 1)
	    return true;
	}
    }
  return false;
}

///===========================================================================
bool ftSurface::allRadialEdges() const
//===========================================================================
{
  for (size_t ki=0; ki<boundary_loops_.size(); ++ki)
    {
      shared_ptr<Loop> aLoop = boundary_loops_[ki];
      if (!aLoop->hasRadialEdges())
	return false;
    }
  return true;
}

//===========================================================================
vector<ftSurface*> ftSurface::fetchCorrespondingFaces() const
//===========================================================================
{
  vector<ftSurface*> result;
//   if (!allRadialEdges())
//     return result;   // Only relevant if all edges have radial edges

  // Get radial edges
  vector<shared_ptr<EdgeVertex> > radial_edges = getRadialEdges();
  if (radial_edges.size() == 0)
    return result;

  // Find faces with relation to all radial edges
  result = radial_edges[0]->getAdjacentFaces();
  size_t ki;
  for (ki=1; ki<radial_edges.size(); ++ki)
    {
      vector<ftSurface*> curr = radial_edges[ki]->getAdjacentFaces();
      
      // Remove faces from result that do not also exist in curr
      size_t kj, kr;
      for (kj=0; kj<result.size();)
	{
	  for (kr=0; kr<curr.size(); ++kr)
	    if (curr[kr] == result[kj])
	      break;

	  if (kr >= curr.size())
	    {
	      // Face not found
	      result.erase(result.begin()+kj);
	    }
	  else
	    kj++;
	}
    }
  // Remove this face
  size_t nmb = result.size();
  for (ki=0; ki<nmb;)
    {
      if (result[ki] == this)
	{
	  result.erase(result.begin()+ki);
	  nmb--;
	}
      else
	ki++;
    }

  return result;
}

//===========================================================================
bool ftSurface::checkFaceTopology()
//===========================================================================
{
  bool isOK = true;
  size_t ki, kj;
  for (ki=0; ki<boundary_loops_.size(); ++ki)
    {
      bool loopface = boundary_loops_[ki]->isFaceConsistent();
      if (!loopface)
	{
	  std::cout << "Loop face inconsistence, loop = " << boundary_loops_[ki] << std::endl;
	  isOK = false;
	}
      else
	{
	  ftFaceBase *face = boundary_loops_[ki]->getEdge(0)->geomEdge()->face();
	  if (face != this)
	    {
	      std::cout << "Inconsistence in face back pointer, face = " << this;
	      std::cout << ", edge = " << boundary_loops_[ki]->getEdge(0);
	      std::cout << ", back pointer face = " << face << std::endl;
	      isOK = false;
	    }
	}
    }

  if (twin_)
    {
      if (twin_->twin() != this)
	{
	  std::cout << "Surface twin inconsistencies, face1 = " << this;
	  std::cout <<", face2 = " << twin_ << std::endl;
	  isOK = false;
	}

      int nmb1 = nmbBoundaryLoops();
      int nmb2 = twin_->nmbBoundaryLoops();
      if (nmb1 != nmb2)
	{
	  std::cout << "Twin inconsistence. Different number of boundary loops. ";
	  std::cout << "Face1 = " << this << ", face2 = " << twin_ << std::cout;
	  isOK = false;
	}
      else
	{
	  // Check consistence of outer loop
	  shared_ptr<Loop> loop1 = getBoundaryLoop(0);
	  shared_ptr<Loop> loop2 = twin_->getBoundaryLoop(0);
	  if (loop1->size() != loop2->size())
	    {
	      std::cout << "Twin inconsistence. Different number of edges in loops. ";
	      std::cout << "Face1 = " << this << ", loop1 = " << loop1;
	      std::cout << ", face2 = " << twin_ << ", loop2 = " << loop2 << std::endl;
	      isOK = false;
	    }
	  else
	    {
	      vector<shared_ptr<ftEdgeBase> > edges1 = loop1->getEdges();
	      vector<shared_ptr<ftEdgeBase> > edges2 = loop2->getEdges();
		  for (ki=0; ki<edges1.size(); ++ki)
		    {
		      for (kj=0; kj<edges2.size(); ++kj)
			if (edges1[ki]->geomEdge()->hasCommonRadialEdge(edges2[kj]->geomEdge()))
			  break;
		      if (kj == edges2.size())
			{
			  std::cout << "Twin inconsistence. Non matching edges in loops. ";
			  std::cout << "Face1 = " << this << ", loop1 = " << loop1;
			  std::cout << ", face2 = " << twin_ << ", loop2 = " << loop2 << std::endl;
			  isOK = false;
			}
		    }
	    }
	}
      // Other loops are not checked for the time being
    }
  return isOK;
}

} // namespace Go
