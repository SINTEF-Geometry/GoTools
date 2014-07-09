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

#include "GoTools/compositemodel/CompleteEdgeNet.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/compositemodel/Path.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/geometry/BoundedSurface.h"
#include <fstream>

//#define DEBUG

using std::vector;
using std::make_pair;
using std::pair;
using namespace Go;

//===========================================================================
CompleteEdgeNet::CompleteEdgeNet(shared_ptr<SurfaceModel> sfmodel,
				 bool perform_step2, bool smooth_connections)
//===========================================================================
  : model_(sfmodel), perform_step2_(perform_step2), 
    smooth_connections_(smooth_connections)
{
}

//===========================================================================
CompleteEdgeNet::~CompleteEdgeNet()
//===========================================================================
{
}

//===========================================================================
bool CompleteEdgeNet::perform(vector<pair<Point,Point> >& corr_vx_pts)
//===========================================================================
{
  // Check if the surface model represents a solid
  int nmb_bd = model_->nmbBoundaries();
  if (nmb_bd > 0)
    return false;  // Not a solid

  // VSK 06.14. This is already expected to be done
  // // Prepare the model for edge completion
  // RegularizeFaceSet regularize(model_, true);
  // model_ = regularize.getRegularModel();

  // Fetch already identified missing edges
  // if (corr_vx_pts.size() == 0)
  //   corr_vx_pts = regularize.fetchVxPntCorr();
  if (corr_vx_pts.size() > 0)
    addIdentifiedEdges(corr_vx_pts);

#ifdef DEBUG
  std::ofstream of("missing_edges0.g2");
  of << "410 1 0 4 255 0 0 255" << std::endl;
  of << missing_edges_.size() << std::endl;
  for (size_t ki=0; ki<missing_edges_.size(); ++ki)
    {
      of << missing_edges_[ki].first->getVertexPoint() << "  ";
      of << missing_edges_[ki].second->getVertexPoint() << std::endl;
    }
#endif

  // Search for and add other missing edges
  addMissingEdges();
  return true;
}


//===========================================================================
void CompleteEdgeNet::addIdentifiedEdges(vector<pair<Point,Point> >& corr_vx_pts)
//===========================================================================
{
  double eps = model_->getTolerances().gap;  // Equality tolerance

  // Fetch all vertices
  vector<shared_ptr<Vertex> > vx;
  model_->getAllVertices(vx);

  for (size_t ki=0; ki<corr_vx_pts.size(); ++ki)
    {
      Point p1 = corr_vx_pts[ki].first;
      Point p2 = corr_vx_pts[ki].second;

      // Find associated vertices
      size_t kj, kr;
      for (kj=0; kj<vx.size(); ++kj)
	if (p1.dist(vx[kj]->getVertexPoint()) < eps)
	  break;

      for (kr=0; kr<vx.size(); ++kr)
	if (p2.dist(vx[kr]->getVertexPoint()) < eps)
	  break;

      if (kj < vx.size() && kr < vx.size())
	{
	  if (vx[kj]->sameEdge(vx[kr].get()))
	    continue;   // Vertices already connected

	  // Look for a vertex situated inbetween the identified one where
	  // the missing edge should connect
	  // @@@ VSK 012014. This is a case specific test which will fail if the
	  // model gets complex enough. Wait for a failure.
	  int idx1, idx2;
	  identifyVertexConnection(vx, kj, kr, idx1, idx2);
	  // @@@ VSK 022014. TEST MORE
	  // if (idx1 < 0 || idx2 < 0)
	  if (idx1 < 0 || idx2 < 0 || idx1 >= vx.size() || idx2 >= vx.size())
	      continue;

	  missing_edges_.push_back(make_pair(vx[idx1], vx[idx2]));
	}
    }
}

//===========================================================================
void CompleteEdgeNet::identifyVertexConnection(vector<shared_ptr<Vertex> > vxs,
					       size_t ki1, size_t ki2, 
					       int& ix1, int& ix2)
//===========================================================================
{
  Point pt1 = vxs[ki1]->getVertexPoint();
  Point pt2 = vxs[ki2]->getVertexPoint();
  Point vec = pt2 - pt1;
  double len = vec.normalize_checked();
  double fac = 0.1;
  double min_ang = 0.05*M_PI;
  double max_ang = 0.25*M_PI;

  // Check the neighbours of the given vertices
  for (size_t kj=0; kj<vxs.size(); ++kj)
    {
      if (kj == ki1 || kj == ki2)
	continue;
      if (vxs[kj]->sameEdge(vxs[ki1].get()) || vxs[kj]->sameEdge(vxs[ki2].get()))
	{
	  // Check distance to the cord between the given vertices
	  Point pt3 = vxs[kj]->getVertexPoint();
	  double t = vec*(pt3 - pt1);
	  if (t <= 0.0 || t >= len)
	    continue;  // Outside of limit

	  Point q = pt1 + t*vec;
	  double dist = pt3.dist(q);
	  if (dist > fac*len)
	    continue;  // Too far distant

	  // Check angles for all edges meeting in the candidate vertex
	  vector<ftEdge*> edges = vxs[kj]->uniqueEdges();
	  size_t kr;
	  for (kr=0; kr<edges.size(); ++kr)
	    {
	      double par = edges[kr]->parAtVertex(vxs[kj].get());
	      Point tan = edges[kr]->tangent(par);
	      double ang = vec.angle(tan);
	      ang = std::min(ang, M_PI-ang);
	      if (edges[kr]->hasVertex(vxs[ki1].get()) || 
		  edges[kr]->hasVertex(vxs[ki2].get()))
		{
		  if (ang > min_ang)
		    break;
		  else
		    min_ang = ang;
		}
	      else
		{
		  if (ang < max_ang)
		    break;
		}
	    }
	  if (kr < edges.size())
	    continue;

	  // A new vertex is found
	  if (vxs[kj]->sameEdge(vxs[ki1].get()) && vxs[kj]->sameEdge(vxs[ki2].get()))
	    {
	      // No connection
	      ix1 = -1;
	      ix2 = -1;
	    }
	  else if (vxs[kj]->sameEdge(vxs[ki1].get()))
	    {
	      ix1 = (int)kj;
	      ix2 = (int)ki2;
	    }
	  else
	    {
	      ix1 = (int)ki1;
	      ix2 =(int) kj;
	    }
	}
    }
}

//===========================================================================
void CompleteEdgeNet::addMissingEdges()
//===========================================================================
{
  // Get start edges
  vector<ftEdge*> start_edges = getStartEdges();
#ifdef DEBUG
  std::ofstream of0("start_edges0.g2");
  for (size_t kh=0; kh<start_edges.size(); ++kh)
    {
      shared_ptr<ParamCurve> cv = start_edges[kh]->geomCurve();
      shared_ptr<ParamCurve> cv2 = 
	shared_ptr<ParamCurve>(cv->subCurve(start_edges[kh]->tMin(),
					    start_edges[kh]->tMax()));
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(of0);
	  sfcv->spaceCurve()->write(of0);
	}
      else
	{
	  cv2->writeStandardHeader(of0);
	  cv2->write(of0);
	}
    }
#endif
  
  // // Fetch all edges in the surface model, twins are represented twice
  // vector<ftEdge*> edges;
  // int nmb_faces = model_->nmbEntities();
  // int ki;
  // for (ki=0; ki<nmb_faces; ++ki)
  //   {
  //     vector<ftEdge*> curr_edges = model_->getFace(ki)->getAllEdgePtrs();
  //     edges.insert(edges.end(), curr_edges.begin(), curr_edges.end());
  //   }

  // Traverse edges to create loops
  while (start_edges.size() > 0)
    {
      ftEdge *curr_edge = start_edges[0];
      shared_ptr<Vertex> vx = curr_edge->getVertex(false);

      vector<ftEdge*> curr_path;
      vector<int> curr_idx;
      traverseEdges(start_edges, curr_path, curr_idx, curr_edge, vx, false);

      ftEdge *twin = curr_edge->twin()->geomEdge();
      vector<ftEdge*>::iterator e1 = std::find(start_edges.begin(), 
					       start_edges.end(),
					       curr_edge);
      if (e1 != start_edges.end())
	start_edges.erase(e1);
      e1 = std::find(start_edges.begin(), start_edges.end(), twin);
      if (e1 != start_edges.end())
	start_edges.erase(e1);

#ifdef DEBUG
      std::ofstream of("curr_missing_edges.g2");
      of << "410 1 0 4 255 0 0 255" << std::endl;
      of << missing_edges_.size() << std::endl;
      for (size_t ki=0; ki<missing_edges_.size(); ++ki)
	{
	  of << missing_edges_[ki].first->getVertexPoint() << "  ";
	  of << missing_edges_[ki].second->getVertexPoint() << std::endl;
	}
#endif

    }

  // This operation has a risk for wrong connections. It is recommended to
  // perform it after a first division of the solid is performed
  if (perform_step2_)
    addRemainingEdges();
}

//===========================================================================
void CompleteEdgeNet::traverseEdges(vector<ftEdge*>& edges,
				    vector<ftEdge*>& curr_path,
				    vector<int>& curr_idx,
				    ftEdge *curr_edge,
				    shared_ptr<Vertex> vx,
				    bool search_end)
//===========================================================================
{
  int next_idx = 0;
  curr_path.push_back(curr_edge);
  curr_idx.push_back(next_idx);

  writePath(curr_path, vx);

#ifdef DEBUG
  std::ofstream out("remainingedges.g2");
  for (size_t kh=0; kh<edges.size(); ++kh)
    {
     shared_ptr<ParamCurve> cv = edges[kh]->geomCurve();
      shared_ptr<ParamCurve> cv2 = 
	shared_ptr<ParamCurve>(cv->subCurve(edges[kh]->tMin(),
					    edges[kh]->tMax()));
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
      if (sfcv.get())
	{
	  sfcv->spaceCurve()->writeStandardHeader(out);
	  sfcv->spaceCurve()->write(out);
	}
      else
	{
	  cv2->writeStandardHeader(out);
	  cv2->write(out);
	}
    }
#endif

  // Check if we have found a loop
  vector<ftEdge*> loop = Path::identifyLoop(curr_path, vx);
  if (loop.size() > 0)
    {
      // Check number of significant vertices and add missing edges if any
      bool found = regularizeEdgeLoop(loop);

      // if (found)
      // 	{
	  // Remove loop from path
	  size_t e1;
	  for (e1=0; e1<curr_path.size(); ++e1)
	    if (curr_path[e1] == loop[0])
	      break;
	  curr_path.erase(curr_path.begin()+e1, curr_path.begin()+e1+loop.size());
	  curr_idx.erase(curr_idx.begin()+e1, curr_idx.begin()+e1+loop.size());
      // 	}
      // else
      // 	{
      // 	  curr_path.pop_back();  // Remove the last and search for a new
      // 	  // loop configuration
      // 	  curr_idx.pop_back();
      // 	}

      // Remove loop from start edges
      for (size_t kh=0; kh<loop.size(); ++kh)
	{
	  ftEdge *twin = loop[kh]->twin()->geomEdge();
	  vector<ftEdge*>::iterator e1 = std::find(edges.begin(), 
						   edges.end(),
						   loop[kh]);
	  if (e1 != edges.end())
	    edges.erase(e1);
	  e1 = std::find(edges.begin(), edges.end(), twin);
	  if (e1 != edges.end())
	    edges.erase(e1);
	}
    
      if (curr_path.size() > 0)
	{
	  curr_edge = curr_path[curr_path.size()-1];
	  vx = curr_edge->getVertex(search_end);
	  if (curr_path.size() > 1)
	    {
	      shared_ptr<Vertex> tmp = 
		curr_path[curr_path.size()-2]->hasVertex(vx.get()) ? vx :
		curr_edge->getVertex(!search_end);
	      vx = curr_edge->getOtherVertex(tmp.get());
	    }
	}
      else
	vx.reset();
    }
  
  if (vx.get())
    {
      // Check if the path ends in an identified missing edge
      size_t ki, kj;
      for (ki=0; ki<missing_edges_.size(); ++ki)
	if (missing_edges_[ki].first.get() == vx.get() || 
	    missing_edges_[ki].second.get() == vx.get())
	  break;

      if (ki < missing_edges_.size())
	{
	  // Check if the start vertex ends in a missing edge as well
	  //shared_ptr<Vertex> vx2 = curr_path[0]->getVertex(!search_end);
	  shared_ptr<Vertex> vx2 = curr_path[0]->getVertex(true);
	  shared_ptr<Vertex> vx3 = curr_path[0]->getVertex(false);
	  // Both ends are already traversed
	  if (search_end && (vx2.get() == missing_edges_[ki].first.get() ||
			     vx2.get() == missing_edges_[ki].second.get() || 
			     vx3.get() == missing_edges_[ki].first.get() ||
			     vx3.get() == missing_edges_[ki].second.get()))
	    kj = ki;
	  else
	    {
	      for (kj=0; kj<missing_edges_.size(); ++kj)
		if (missing_edges_[kj].first.get() == vx2.get() || 
		    missing_edges_[kj].second.get() == vx2.get() ||
		    missing_edges_[kj].first.get() == vx3.get() || 
		    missing_edges_[kj].second.get() == vx3.get())
		  break;
	    }

#ifdef DEBUG
	      std::ofstream out2("curr_missing_edges.g2");
	      out2 << "410 1 0 4 155 100 0 255 " << std::endl;
	      out2 << missing_edges_.size() << std::endl;
	      for (size_t km=0; km<missing_edges_.size(); ++km)
		{
		  out2 << missing_edges_[km].first->getVertexPoint() << " "; 
		  out2 << missing_edges_[km].second->getVertexPoint() << std::endl;
		}
#endif
		
	  if (ki == kj)
	    {
	      // The path meets the same missing edge in both ends
	      // A loop is found
	      bool found = regularizeEdgeLoop(curr_path);
	      //	      if (found)
		curr_path.clear();
	      // else
	      // 	{
	      // 	  curr_path.pop_back();
	      // 	}

		// Remove path from start edges
		for (size_t kh=0; kh<curr_path.size(); ++kh)
		  {
		    ftEdge *twin = curr_path[kh]->twin()->geomEdge();
		    vector<ftEdge*>::iterator e1 = std::find(edges.begin(), 
							     edges.end(),
							     curr_path[kh]);
		    if (e1 != edges.end())
		      edges.erase(e1);
		    e1 = std::find(edges.begin(), edges.end(), twin);
		    if (e1 != edges.end())
		      edges.erase(e1);
		  }

	      curr_edge = 0;
	      vx.reset();
	    }
	  else if (search_end || kj < missing_edges_.size())
	    {
	      // Will not perform a legal loop. Step one edge back
	      // Try to continue along another edge
	      writePath(curr_path, vx);

	      int stop_break0 = 1;
	      if (kj < missing_edges_.size())
		{
		  if ((vx.get() == missing_edges_[ki].first.get() &&
		      (missing_edges_[ki].second.get() == missing_edges_[kj].first.get() ||
		       missing_edges_[ki].second.get() == missing_edges_[kj].second.get())) ||
		      (vx.get() == missing_edges_[ki].second.get() &&
		      (missing_edges_[ki].first.get() == missing_edges_[kj].first.get() ||
		       missing_edges_[ki].first.get() == missing_edges_[kj].second.get())))
		    {
		      curr_path.pop_back();
		      curr_idx.pop_back();
		      if (curr_path.size() > 0)
			{
			  curr_edge = curr_path[curr_path.size()-1];
			  vx = curr_edge->getVertex(search_end);
			  if (curr_path.size() > 1)
			    {
			      shared_ptr<Vertex> tmp = 
				curr_path[curr_path.size()-2]->hasVertex(vx.get()) ? vx :
				curr_edge->getVertex(!search_end);
			      vx = curr_edge->getOtherVertex(tmp.get());
			    }
			}
		      else
			vx.reset();
		    }
		}
	    }
	  else
	    {
	      // Turn path and continue search from the other end
	      writePath(curr_path, vx);
	      std::reverse(curr_path.begin(), curr_path.end());
	      search_end = !search_end;
	      std::reverse(curr_idx.begin(), curr_idx.end());
	      curr_edge = curr_path[curr_path.size()-1]; 
	      vx = curr_edge->getVertex(search_end);
	      if (curr_path.size() > 1)
		{
		  shared_ptr<Vertex> tmp = 
		    curr_path[curr_path.size()-2]->hasVertex(vx.get()) ? vx :
		    curr_edge->getVertex(!search_end);
		  vx = curr_edge->getOtherVertex(tmp.get());
		}
	    }
	}
    }

  ftEdge *next_edge = 0;
  while (vx.get() && !next_edge)
    {
      // Fetch next legal edge
      next_idx = curr_idx[curr_idx.size()-1];
      next_edge = fetchNextEdge(curr_edge, vx, next_idx);
      curr_idx[curr_idx.size()-1] = next_idx;
      if (!next_edge)
	{
	  writePath(curr_path, vx);
	  curr_path.pop_back();
	  curr_idx.pop_back();
	  if (curr_path.size() > 0)
	    {
	      curr_edge = curr_path[curr_path.size()-1];
	      vx = curr_edge->getVertex(search_end);
	      if (curr_path.size() > 1)
		{
		  shared_ptr<Vertex> tmp = 
		    curr_path[curr_path.size()-2]->hasVertex(vx.get()) ? vx :
		    curr_edge->getVertex(!search_end);
		  vx = curr_edge->getOtherVertex(tmp.get());
		}
	    }
	  else
	    vx.reset();
	}
    }

  if (next_edge)
    {
      vx = next_edge->getOtherVertex(vx.get());
      traverseEdges(edges, curr_path, curr_idx, next_edge, vx, search_end);
    }
}

//===========================================================================
ftEdge* CompleteEdgeNet::fetchNextEdge(ftEdge *curr_edge,
				       shared_ptr<Vertex> vx,
				       int& next_idx)
//===========================================================================
{
  ftEdge *next = NULL;

  ftEdge *twin = curr_edge->twin()->geomEdge();
  if (!twin)
    return next;  // Should not happen

  // Fetch faces associate to the current edge
  vector<ftSurface*> faces1;
  Body *bd = NULL;
  if (curr_edge->face())
    bd = curr_edge->face()->asFtSurface()->getBody();
  if (curr_edge->hasEdgeMultiplicity())
    faces1 = curr_edge->getEdgeMultiplicityInstance()->getAdjacentFaces(bd);
  else
    {
      faces1.push_back(curr_edge->face()->asFtSurface());
      faces1.push_back(twin->face()->asFtSurface());
    }

  // Fetch vertex edges
  vector<ftEdge*> vx_edges = vx->uniqueEdges(bd);
  int ki;
  for (ki=next_idx; ki<(int)vx_edges.size(); ++ki)
    {
      if (vx_edges[ki] == curr_edge || vx_edges[ki] == twin)
	continue;  // Current edge

      // Fetch assocated faces and check that they are different from
      // the current ones
      vector<ftSurface*> faces2;
      if (vx_edges[ki]->hasEdgeMultiplicity())
      	faces2 = vx_edges[ki]->getEdgeMultiplicityInstance()->getAdjacentFaces(bd);
      else
      	{
      	  faces2.push_back(vx_edges[ki]->face()->asFtSurface());
      	  faces2.push_back(vx_edges[ki]->twin()->face()->asFtSurface());
      	}
      size_t kj, kr;
      for (kj=0; kj<faces1.size(); ++kj)
      	{
      	  for (kr=0; kr<faces2.size(); ++kr)
      	    if (faces1[kj] == faces2[kr])
      	      break;
      	  if (kr < faces2.size())
      	    break;
      	}
      if (kj < faces1.size())
      	continue;
      
      break;
    }
     
  if (ki < (int)vx_edges.size())
    next = vx_edges[ki];

  next_idx = ki + 1;
  return next;
}

//===========================================================================
bool CompleteEdgeNet::regularizeEdgeLoop(vector<ftEdge*>& edges)
//===========================================================================
{

  bool to_add_edges = false;
  bool found = false;

  // Fetch vertices
  vector<shared_ptr<Vertex> > vxs;
  int ki;
  //int ki, kj, kr, kh;
  shared_ptr<Vertex> tmp_vx = edges[edges.size()-1]->getVertex(false);
  if (edges[0]->hasVertex(edges[edges.size()-1]->getVertex(true).get()))
    tmp_vx = edges[edges.size()-1]->getVertex(true);
  for (ki=0; ki<(int)edges.size(); ++ki)
    {
      if (!edges[ki]->hasVertex(tmp_vx.get()))
	{
	  tmp_vx = edges[ki]->getVertex(true);
	  if (ki+1 < (int)edges.size() && edges[ki+1]->hasVertex(tmp_vx.get()))
	    tmp_vx = edges[ki]->getVertex(false);
	}
	vxs.push_back(tmp_vx);
	tmp_vx = edges[ki]->getOtherVertex(tmp_vx.get());
    }
  if (tmp_vx.get() != vxs[0].get())
    {
      vxs.push_back(tmp_vx);
      to_add_edges = true;  // Does already contain a missing edge,
      // check if more are missing
    }

#ifdef DEBUG
  std::ofstream of("edge_loop.g2");
  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      shared_ptr<ParamCurve> curr_crv = edges[ki]->geomCurve();
      shared_ptr<ParamCurve> curr_crv2 = 
	shared_ptr<ParamCurve>(curr_crv->subCurve(edges[ki]->tMin(),
						  edges[ki]->tMax()));
      shared_ptr<CurveOnSurface> sf_cv = 
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curr_crv2);
      if (sf_cv.get())
	{
	  sf_cv->spaceCurve()->writeStandardHeader(of);
	  sf_cv->spaceCurve()->write(of);
	}
      else
	{
	  curr_crv2->writeStandardHeader(of);
	  curr_crv2->write(of);
	}

      Point pnt = vxs[ki]->getVertexPoint();
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << "1" << std::endl;
      of << pnt << std::endl;
    }

  if (vxs.size() > edges.size())
    {
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << "1" << std::endl;
      of << vxs[vxs.size()-1]->getVertexPoint() << std::endl;
    }
#endif

  if (vxs.size() <= 4)
    return found;  // No need for extra edges


  // // Count the number of the same underlying surface on both sides of the
  // // edge
  // int nmb_same = 0;
  // for (ki=0; ki<edges.size(); ++ki)
  //   {
  //     ftSurface *f1 = edges[ki]->face()->asFtSurface();
  //     ftSurface *f2 = edges[ki]->twin()->geomEdge()->face()->asFtSurface();
  //     shared_ptr<BoundedSurface> bd_sf1 = 
  // 	dynamic_pointer_cast<BoundedSurface,ParamSurface>(f1->surface());
  //     shared_ptr<BoundedSurface> bd_sf2 = 
  // 	dynamic_pointer_cast<BoundedSurface,ParamSurface>(f2->surface());
  //     if (bd_sf1.get() && bd_sf2.get() && 
  // 	  bd_sf1->underlyingSurface().get() == bd_sf2->underlyingSurface().get())
  // 	nmb_same++;
  //   }

  // std::cout << "Number of edges: " << edges.size();
  // std::cout <<", number of same surface: " << nmb_same << std::endl;

  // Split the edge loop where there exists additional edges between
  // vertices in the loop 
  vector<vector<ftEdge*> > split_loops;
  vector<vector<shared_ptr<Vertex> > > split_vxs;
  vector<bool> add_edges_split;
  splitLoop(edges, vxs, to_add_edges, split_loops, split_vxs, 
	    add_edges_split);

  // Treat each edges separately
  for (ki=0; ki<(int)split_loops.size(); ++ki)
    {
      bool found2 = regularizeCurrLoop(split_loops[ki], split_vxs[ki], 
				       add_edges_split[ki]);
      if (found2)
	found = true;
    }
  return found;
}

//===========================================================================
void CompleteEdgeNet::splitLoop(vector<ftEdge*>& edges,
				vector<shared_ptr<Vertex> >& vxs,
				bool to_add_edges,
				vector<vector<ftEdge*> >& split_loops,
				vector<vector<shared_ptr<Vertex> > >& split_vxs,
				vector<bool>& add_edges_split)
//===========================================================================
{
  size_t idx1 = 0;
  size_t idx2 = vxs.size()-1;
  size_t ki, kj, kr;

  // Traverse the edge loop and extract sub loops
  for (ki=idx1; ki<=idx2; ++ki)
    {
      for (kj=ki+2; kj<=idx2; ++kj)
	{
	  // Check if the two vertices are connected
	  if (vxs[ki]->sameEdge(vxs[kj].get()) && (ki > idx1 || kj < idx2))
	    {
	      // Extract local loop
	      // First vertices
	      vector<shared_ptr<Vertex> > curr_vxs;
	      for (kr=idx1; kr<=ki; ++kr)
		curr_vxs.push_back(vxs[kr]);
	      for (kr=kj; kr<=idx2; ++kr)
		curr_vxs.push_back(vxs[kr]);
	      split_vxs.push_back(curr_vxs);

	      // Edges
	      vector<ftEdge*> curr_edges;
	      for (kr=idx1; kr<ki; ++kr)
		curr_edges.push_back(edges[kr]);
	      curr_edges.push_back(vxs[ki]->getCommonEdge(vxs[kj].get()));
	      for (kr=kj; kr<idx2; ++kr)
		curr_edges.push_back(edges[kr]);
	      if (split_loops.size() > 0)
		curr_edges.push_back(curr_vxs[0]->getCommonEdge(curr_vxs[curr_vxs.size()-1].get()));
	      split_loops.push_back(curr_edges);

	      // Joint parameter
	      add_edges_split.push_back(add_edges_split.size() == 0 ? 
					to_add_edges : false);

	      // Set indices
	      idx1 = ki;
	      idx2 = kj;
	    }
	}
    }
      
  // The last loop
  // First vertices
  vector<shared_ptr<Vertex> > curr_vxs;
  for (kr=idx1; kr<=idx2; ++kr)
    curr_vxs.push_back(vxs[kr]);
  split_vxs.push_back(curr_vxs);

  // Edges
  vector<ftEdge*> curr_edges;
  for (kr=idx1; kr<idx2; ++kr)
    curr_edges.push_back(edges[kr]);
  if (split_loops.size() > 0 || (!to_add_edges))
    curr_edges.push_back(curr_vxs[0]->getCommonEdge(curr_vxs[curr_vxs.size()-1].get()));
  split_loops.push_back(curr_edges);

  // Joint parameter
  add_edges_split.push_back(add_edges_split.size() == 0 ? 
			    to_add_edges : false);
  
}

//===========================================================================
bool CompleteEdgeNet::regularizeCurrLoop(vector<ftEdge*>& edges,
					 vector<shared_ptr<Vertex> >& vxs,
					 bool to_add_edges)
//===========================================================================
{
  bool found = false;
  if (vxs.size() <= 4)
    return found;  // Not required to add edges

  Body *bd = model_->getBody();

  // Check for repeated vertices
  int ki, kj, kr, kh;
  for (ki=0; ki<(int)vxs.size(); ++ki)
    for (kj=ki+1; kj<(int)vxs.size(); ++kj)
      if (vxs[ki].get() == vxs[kj].get())
	return found;

#ifdef DEBUG
  std::ofstream ofc("curr_edge_loop.g2");
  for (ki=0; ki<(int)edges.size(); ++ki)
    {
      shared_ptr<ParamCurve> curr_crv = edges[ki]->geomCurve();
      shared_ptr<ParamCurve> curr_crv2 = 
	shared_ptr<ParamCurve>(curr_crv->subCurve(edges[ki]->tMin(),
						  edges[ki]->tMax()));
      shared_ptr<CurveOnSurface> sf_cv = 
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curr_crv2);
      if (sf_cv.get())
	{
	  sf_cv->spaceCurve()->writeStandardHeader(ofc);
	  sf_cv->spaceCurve()->write(ofc);
	}
      else
	{
	  curr_crv2->writeStandardHeader(ofc);
	  curr_crv2->write(ofc);
	}

      Point pnt = vxs[ki]->getVertexPoint();
      ofc << "400 1 0 4 0 255 0 255" << std::endl;
      ofc << "1" << std::endl;
      ofc << pnt << std::endl;
    }

  if (vxs.size() > edges.size())
    {
      ofc << "400 1 0 4 0 255 0 255" << std::endl;
      ofc << "1" << std::endl;
      ofc << vxs[vxs.size()-1]->getVertexPoint() << std::endl;
    }
#endif

  // Check if the plane(s) defined by the edge loop are significantly
  // different from the tangent planes of the assiciated faces
  // To check if the loop should be prosessed further
  double bend_tol = model_->getTolerances().bend;
  if (!to_add_edges)
    {
      int nmb_plane = 0;
      int nmb_check = 0;
      for (ki=0; ki<(int)edges.size(); ++ki)
	{
	    kj = (ki == 0) ? (int)edges.size()-1 : ki-1; 
	  //kr = (ki == edges.size()-1) ? 0 : ki+1;
	  Point tan1 = edges[kj]->tangent(edges[kj]->tMin());
	  Point tan2 = edges[ki]->tangent(edges[ki]->tMax());
	  if (tan1.angle(tan2) < bend_tol)
	    continue;  // Not a clear plane

	  Point norm = tan1.cross(tan2);
	  if (norm.length() < model_->getTolerances().gap)
	    continue;

	  ftSurface *f1 = edges[ki]->face()->asFtSurface();
	  ftSurface *f2 = edges[ki]->twin()->geomEdge()->face()->asFtSurface();
	  Point par1 = vxs[ki]->getFacePar(f1);
	  Point par2 = vxs[ki]->getFacePar(f2);
	  Point norm1 = f1->normal(par1[0], par1[1]);
	  Point norm2 = f2->normal(par2[0], par2[1]);
	  double ang1 = std::min(norm.angle(norm1), norm.angle(-norm1));
	  double ang2 = std::min(norm.angle(norm2), norm.angle(-norm2));
	  double ang = std::min(ang1, ang2);
	  nmb_check++;
	  if (ang < bend_tol)
	    nmb_plane++;
	}

#ifdef DEBUG
      std::cout << "Nmb check: " << nmb_check << ", nmb plane: " << nmb_plane << std::endl;
#endif

      if (nmb_plane < (int)(0.75*nmb_check + 1))
	//if (nmb_plane < (int)(0.5*nmb_check + 1))
	//if (nmb_plane < (int)(0.9*nmb_check + 1))
	to_add_edges = true;
    }

  if (!to_add_edges)
    return found;

  // Select vertices between which to add edges
  // First project all vertices into the best plane given by the vertices
  vector<Point> vxs2(vxs.size());
  Point norm(0.0, 0.0, 0.0);
  Point pos(0.0, 0.0, 0.0);
  int nmb_vec = 0;
  for (ki=0; ki<(int)vxs.size(); ++ki)
    {
      pos += vxs[ki]->getVertexPoint();
      for (kj=ki+1; kj<(int)vxs.size(); ++kj)
	{
	  Point vec1 = vxs[kj]->getVertexPoint() - vxs[ki]->getVertexPoint();
	  for (kr=ki+1; kr<(int)vxs.size(); ++kr)
	    for (kh=kr+1; kh<(int)vxs.size(); ++kh)
	      {
		Point vec2 = vxs[kh]->getVertexPoint() - vxs[kr]->getVertexPoint();
		Point tmp = vec1.cross(vec2);
		if (tmp.length() > model_->getTolerances().gap)
		  {
		    tmp.normalize();
		    norm += tmp;
		    nmb_vec++;
		  }
	      }
	}
    }
  norm /= (double)nmb_vec;
  norm.normalize();
  pos /= (double)(vxs.size());

  for (ki=0; ki<(int)vxs.size(); ++ki)
    {
      Point vxpt = vxs[ki]->getVertexPoint();
      vxs2[ki] = vxpt - ((vxpt-pos)*norm)*norm;
    }

#ifdef DEBUG
  std::ofstream of2d("proj_edges.g2");
  for (ki=1; ki<(int)vxs.size(); ++ki)
    {
      of2d << "410 1 0 4 0 255 0 255" << std::endl;
      of2d << "1" << std::endl;
      of2d << vxs2[ki-1] << "  " << vxs2[ki] << std::endl;
    }
  of2d << "410 1 0 4 0 255 0 255" << std::endl;
  of2d << "1" << std::endl;
  of2d << vxs2[vxs2.size()-1] << "  " << vxs2[0] << std::endl;
#endif

  // Check feasability of loop. All non-corner vertices cannot
  // lie between the same two corners
  int last_corner = -1;
  int first_corner = -1;
  int last_idx = -1;
  int first_idx = -1;
  int nmb_corner = 0;
  double corner_ang = 0.1*M_PI;
  for (ki=0; ki<(int)vxs2.size(); ++ki)
    {
      kr = (ki==0) ? vxs2.size()-1 : ki-1;
      kj = (ki+1)%(vxs2.size());
      Point vec1 = vxs2[kj] - vxs2[ki];
      Point vec2 = vxs2[ki] - vxs2[kr];
      double ang = vec1.angle(vec2);
      if (ang > corner_ang) 
	{
	  last_corner = ki;
	  if (first_corner < 0)
	    first_corner = ki;
	  nmb_corner++;
	}
      else
	{
	  if (last_idx < 0)
	    first_idx = last_idx = ki;
	  else if (last_idx < first_corner)
	    last_idx = ki;
	  else if (last_corner > last_idx && last_corner > first_corner)
	    break;
	}
    }
  if (ki == vxs2.size() && nmb_corner <= 4 && last_idx >= 0 &&
      ((first_idx < first_corner && last_idx > last_corner) || first_idx == last_idx))
    return found;  // All non corner vertices are placed between two corners

  // Select start vertex for the construction of missing edges
  int vx_idx = -1;
  int vx_idx2 = -1;
  double min_ang = MAXDOUBLE;
  double min_dist = MAXDOUBLE;
  double ang_tol = model_->getTolerances().kink;
  double dist_fac = 10.0;
  if (false /*vxs.size() > edges.size()*/)
    {
      // The edge loop start and ends in an already added edge,
      // continue to add in the same pattern
      vx_idx = (int)vxs.size() - 2;
      vx_idx2 = 1;

      // Check
      for (size_t kn=0; kn<missing_edges_.size(); ++kn)
	{
	  if (missing_edges_[kn].first.get() == vxs[vx_idx].get() || 
	      missing_edges_[kn].first.get() == vxs[vx_idx2].get() ||
	      missing_edges_[kn].second.get() == vxs[vx_idx].get() || 
	      missing_edges_[kn].second.get() == vxs[vx_idx2].get())
	    {
	      vx_idx = -1;
	      break;
	    }
	}
    }
  else
    {
	for (ki=0, kj=ki+1, kr=kj+1, kh=kr+1; ki<(int)vxs2.size(); 
	   ++ki, ++kj, ++kr, ++kh)
	{
	    kj = kj % (int)vxs2.size();
	    kr = kr % (int)vxs2.size();
	    kh = kh % (int)vxs2.size();

	  Point vec1 = vxs2[kh] - vxs2[ki];
	  Point vec2 = vxs2[kr] - vxs2[kj];
	  double ang = vec1.angle(vec2);

	  // Check if the two intermediate vertices intersects the candidate
	  double v1len = vec1.length();
	  double v12 = vec1*vec1;
	  double mfac = 0.001;
	  Point v2 = vxs2[kj] - vxs2[ki];
	  double td = (v2*vec1)/v12;
	  double dd = vxs2[kj].dist(vxs2[ki]+td*vec1);
	  if (td >= 0.0 && td <= 1.0 && dd <= mfac*v1len)
	    continue;
	  v2 = vxs2[kr] - vxs2[ki];
	  td = (v2*vec1)/v12;
	  dd = vxs2[kr].dist(vxs2[ki]+td*vec1);
	  if (td >= 0.0 && td <= 1.0 && dd <= mfac*v1len)
	    continue;

	  int ka;
	  for (ka=(kh+1)%((int)vxs2.size()); ka!=ki; ka=(ka+1)%((int)vxs2.size()))
	    {
	      v2 = vxs2[ka] - vxs2[ki];
	      td = (v2*vec1)/v12;
	      dd = vxs2[ka].dist(vxs2[ki]+td*vec1);
	      if (td >= 0.0 && td <= 1.0 && dd <= mfac*v1len)
		break;
	    }
	  if (ka != ki)
	    continue;

	  // Check if the candidate crosses the polygon of projected
	  // vertices
	  size_t kk;
	  double tol = 1.0e-8;
	  Point n1 = vec1.cross(norm);
	  Point p1 = 0.5*(vxs2[ki] + vxs2[kh]);
	  for (kk=1; kk<vxs2.size(); ++kk)
	    {
	      Point vecn = vxs2[kk] - vxs2[kk-1];
	      Point n2 = vecn.cross(norm);
	      Point p2 = 0.5*(vxs2[kk-1] + vxs2[kk]);
	      double angn = vec1.angle(vecn);
	      double t1, s1;
	      if (angn < ang_tol)
		{
		  //t1 = s1 = (vxs2[ki] - p2)*(p2 - vxs2[kh]);
		  Point tmp = vec1;
		  tmp.normalize();
		  s1 = (p2 - vxs2[ki])*tmp;
		  t1 = (p2 - vxs2[ki] - s1*tmp).length();
		  if (t1 < tol)
		    break;  // Coincidence (possibly in the extension)
		}
	      else
		{
		  t1 = ((p2 - vxs2[ki])*n2)/(vec1*n2);
		  s1 = ((p1 - vxs2[kk-1])*n1)/((vxs2[kk]-vxs2[kk-1])*n1);
		  if (t1 > tol && t1 < 1.0-tol && s1 > tol && s1 < 1.0-tol)
		    break;
		}
	    }
	  
	  if (kk <vxs2.size())
	    continue;

	  Point vecn = vxs2[0] - vxs2[vxs2.size()-1];
	  Point n2 = vecn.cross(norm);
	  Point p2 = 0.5*(vxs2[vxs2.size()-1] + vxs2[0]);
	  double angn = vec1.angle(vecn);
	  double t1, s1;
	  if (angn < ang_tol)
	    {
	      t1 = s1 = (vxs2[ki] - p2)*(p2 - vxs2[kh]);
	    }
	  else
	    {
	      t1 = ((p2 - vxs2[ki])*n2)/(vec1*n2);	
	      s1 = ((p1 - vxs2[vxs2.size()-1])*n1)/((vxs2[0]-vxs2[vxs2.size()-1])*n1);
	    }
	  if (t1 > tol && t1 < 1.0-tol && s1 < 1.0-tol)
	    continue;
	  
	  Point vec3 = vxs2[kh] - vxs2[kr];
	  Point vec4 = vxs2[kj] - vxs2[ki];
	  double ang2 = vec1.angle(vec3);
	  double ang3 = vec2.angle(vec4);
	  double ang4 = vec3.angle(vec4);
	  double dist = vxs[kh]->getDist(vxs[ki]);
	  double fac = 1.5;
	  if (ki == 0 || 
	      (ang < min_ang-ang_tol && ang2 > fac*ang && 
	       ang*dist < fac*min_ang*min_dist) ||
	      (dist < dist_fac*min_dist && ang*dist < fac*min_ang*min_dist))
	    {
	      min_ang = ang;
	      min_dist = dist;
	      vx_idx = ki;
	    }
	  else if (fabs(ang - min_ang) < ang_tol && ang2 > fac*ang && 
		   ang*dist < fac*min_ang*min_dist &&
		   dist < min_dist)
	    {
	      min_ang = ang;
	      min_dist = dist;
	      vx_idx = ki;
	    }
	}
    }

  // @@@ Preferably, we should compare the missing edge direction with
  // the directions of both existing corresponding edges. This is currently
  // not done

  if (vx_idx < 0)
    return found;  // No combination is found

#ifdef DEBUG
  std::ofstream of("edge_loop2.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << vxs[vx_idx]->getVertexPoint() << std::endl;
#endif

  // Reduced set of edges and vertices
  vector<ftEdge*> red_edges;
  vector<shared_ptr<Vertex> > red_vxs;
  vector<shared_ptr<ftEdge> > extra_edges;
  size_t curr_nmb_edges = missing_edges_.size();
  bool found_already = false;

  // Combine vertices
  Point prev_proj_vec;
  if (vx_idx2 < 0)
      vx_idx2 = (vx_idx+3) % (int)vxs2.size();
  for (ki=vx_idx, kj=vx_idx2; abs(int(kj-ki))>1;
       ki=(ki>0)?ki-1:(int)vxs2.size()-1, kj=(kj+1)%(int)vxs2.size())
    {
	if ((ki==0 && kj==(int)vxs2.size()-1) ||
	    (kj==0 && ki==(int)vxs2.size()) ||
	  abs(int(ki-kj)) <= 1)
	break;

#ifdef DEBUG
	std::ofstream ofcc("conn_cand.g2");
	ofcc << "400 1 0 4 0 255 0 255" << std::endl;
	ofcc << "1" << std::endl;
	ofcc << vxs[ki]->getVertexPoint() << std::endl;
	ofcc << "400 1 0 4 155 100 0 255" << std::endl;
	ofcc << "1" << std::endl;
	ofcc << vxs[kj]->getVertexPoint() << std::endl;
#endif

	
	// Check if both vertices are T-joints
	if (vxs[ki]->nmbUniqueEdges(bd) <= 2)
	  ki=(ki>0)?ki-1:(int)vxs2.size()-1;
	if (vxs[kj]->nmbUniqueEdges(bd) <= 2)
	  kj=(kj+1)%(int)vxs2.size();
	if ((ki==0 && kj==(int)vxs2.size()-1) ||
	    (kj==0 && ki==(int)vxs2.size()) ||
	  abs(int(ki-kj)) <= 1)
	break;

      if (vxs[ki]->sameEdge(vxs[kj].get()))
	  continue;  // Already connected

      if (vxs[ki]->sameFace(vxs[kj].get()))
	  continue;  // Don't make connections across an existing face

      if (vxs[ki]->sameUnderlyingSurface(vxs[kj].get()))
	  continue;  // Don't make connections across an existing surface

      // Check if the two vertices are joined to the same vertex
      // through edges
      // @@@ VSK. This test is probably not sufficient for all cases,
      // but I wait until a more complex case turns up before I try
      // to refine it
      if (vxs[ki]->connectedToSameVertex(vxs[kj].get()))
	continue;

      // Check if any of the current vertices are connected to a
      // non-neighbouring vertex in the loop already
      for (kr=0; kr<(int)vxs.size(); ++kr)
	{
	  if (kr == vx_idx || abs(vx_idx-kr) == 1 ||
	      abs(kr-vx_idx) == (int)vxs.size()-1)
	    continue;
	    if (vxs[kr]->sameEdge(vxs[vx_idx].get()))
	      break;
	  }
      if (kr < (int)vxs.size())
	continue;  // No not connect
      
       for (kr=0; kr<(int)vxs.size(); ++kr)
	{
	  if (kr == vx_idx2 || abs(vx_idx2-kr) == 1 ||
	      abs(kr-vx_idx2) == (int)vxs.size()-1)
	    continue;
	    if (vxs[kr]->sameEdge(vxs[vx_idx2].get()))
	      break;
	  }
      if (kr < (int)vxs.size())
	continue;  // No not connect

      if (!smooth_connections_ && vxs[ki]->nmbUniqueEdges() == 4 &&
	  vxs[kj]->nmbUniqueEdges() == 4)
	{
	  // Do not make a connection between two vertices that both
	  // are surrounded with faces defining the same plane in the
	  // vertex
	  vector<pair<ftSurface*, Point> > faces1 = vxs[ki]->getFaces();
	  vector<pair<ftSurface*, Point> > faces2 = vxs[kj]->getFaces();
	  
	  Point norm1 = faces1[0].first->normal(faces1[0].second[0],
						faces1[0].second[1]);
	  for (kr=1; kr<faces1.size(); ++kr)
	    {
	      Point norm = faces1[kr].first->normal(faces1[kr].second[0],
						    faces1[kr].second[1]);
	      if (norm1.angle(norm) > bend_tol)
		break;
	    }

	  Point norm2 = faces2[0].first->normal(faces2[0].second[0],
						faces2[0].second[1]);
	  for (kh=1; kh<faces2.size(); ++kh)
	    {
	      Point norm = faces2[kh].first->normal(faces2[kh].second[0],
						    faces2[kh].second[1]);
	      if (norm2.angle(norm) > bend_tol)
		break;
	    }

	  if (kr == faces1.size() && kh == faces2.size())
	    continue;  // Do not connect at the current stage
	}

      // // Check if the current candidate connection bypasses a better connection
      // // in the same face
      // if (betterConnectionInFace(vxs[ki], bd, 
      // 				 edges[(ki<edges.size())?ki:(int)edges.size()-1], 
      // 				 edges[(ki>0)?ki-1:(int)edges.size()-1], vxs[kj]))
      // 	continue;
      // if (betterConnectionInFace(vxs[kj], bd, 
      // 				 edges[(kj<edges.size())?kj:(int)edges.size()-1],
      // 				 edges[(kj>0)?kj-1:(int)edges.size()-1], vxs[ki]))
      // 	continue;
      
      // Check if the new connection intersects another vertex in the loop
      Point v1 = vxs2[kj] - vxs2[ki];
      double v1len = v1.length();
      double v12 = v1*v1;
      double mfac = 0.1;
      for (kr=0; kr<(int)vxs.size(); ++kr)
	{
	  if (kr==ki || kr==kj)
	    continue;

	  Point v2 =  vxs2[kr] - vxs2[ki];
	  double td = (v2*v1)/v12;
	  double dd = vxs2[kr].dist(vxs2[ki] + td*v1);
	  if (td >= 0.0 && td <= 1.0 && dd <= mfac*v1len)
	    break;
	}

      if (found && kr==(int)vxs.size())
	{
	  int ki2 = (ki>0)?ki-1:(int)vxs2.size()-1;
	  int kj2 = (kj+1)%(int)vxs2.size();
	  Point v2 = vxs2[kj2] - vxs2[ki2];
	  double proj_angle = prev_proj_vec.angle(v1) + v1.angle(v2);
	  if (proj_angle > 0.2*M_PI)
	    kr = 0;   // Too large angular difference between adjacent
	  // connections. Treat loop recursively
	}

      prev_proj_vec = v1;

      if (kr < (int)vxs.size())
	{
	  if (found)
	    {
	      // Create reduced set of edges and vertices
	      kh = ki;
	      while (true)
		{
		  red_vxs.push_back(vxs[kh]);
		  if (kh < edges.size())
		    red_edges.push_back(edges[kh]);
		  kr = (kh + 1)%(int)(vxs.size());
		  size_t kn;
		  for (kn=curr_nmb_edges; kn<missing_edges_.size(); ++kn)
		    {
		      if (missing_edges_[kn].first.get() == vxs[kr].get() || 
			  missing_edges_[kn].second.get() == vxs[kr].get())
			break;
		    }
		  if (kn < missing_edges_.size())
		    {
		      shared_ptr<Vertex> vxn1 = missing_edges_[kn].first;
		      shared_ptr<Vertex> vxn2 = missing_edges_[kn].second;
		      shared_ptr<ParamCurve> cv = 
			shared_ptr<ParamCurve>(new SplineCurve(vxn1->getVertexPoint(),
							       vxn2->getVertexPoint()));
		      
		      shared_ptr<ftEdge> edgen = 
			shared_ptr<ftEdge>(new ftEdge(cv, cv->startparam(), cv->endparam()));
		      red_vxs.push_back(vxs[kr]);
		      extra_edges.push_back(edgen);
		      red_edges.push_back(edgen.get());
		      if (vxn2.get() == vxs[kr].get())
			std::swap(vxn1, vxn2);
		      for (kh=0; kh<(int)vxs.size(); ++kh)
			if (vxs[kh].get() == vxn2.get())
			  break;
		    }
		  else
		    {
		      kh = kr;
		    }

		  if (kh == ki)
		    break;
		}
	    }
	  break;
	}

      found = true;

      // Check if the missing edge exists already
      size_t kn;
      for (kn=0; kn<missing_edges_.size(); ++kn)
	if ((missing_edges_[kn].first.get() == vxs[ki].get() && 
	     missing_edges_[kn].second.get() == vxs[kj].get()) ||
	    (missing_edges_[kn].first.get() == vxs[kj].get() && 
	     missing_edges_[kn].second.get() == vxs[ki].get()))
	  {
	    found_already = true;
	    break;
	  }

      if (kn == missing_edges_.size())
	missing_edges_.push_back(make_pair(vxs[ki], vxs[kj]));

#ifdef DEBUG      
      of << "410 1 0 4 255 0 0 255" << std::endl;
      of << "1" << std::endl;
      of << vxs[ki]->getVertexPoint() << " " << vxs[kj]->getVertexPoint() << std::endl;
#endif
    }
      
  if (red_edges.size() > 0 && red_edges.size() < edges.size() && !found_already)
    {
      bool found2 = regularizeCurrLoop(red_edges, red_vxs, true);
      if (found2)
	found = true;
      for (ki=0; ki<(int)extra_edges.size(); ++ki)
	{
	  shared_ptr<ftEdge> curr = extra_edges[ki];
	  shared_ptr<Vertex> v1, v2;
	  curr->getVertices(v1, v2);

	  // Remove edge from associated vertices
	  v1->removeEdge(curr.get());
	  v2->removeEdge(curr.get());
	}
    }      

  int stop_break;
  stop_break = 1;

  return found;
}

//===========================================================================
double CompleteEdgeNet::getVertexAngle(ftEdge *edge1, ftEdge *edge2)
//===========================================================================
{
  Point tan1 = edge1->tangent(edge1->tMax());
  tan1.normalize();
  Point tan2 = edge2->tangent(edge2->tMin());
  tan2.normalize();

  // Project into tangent plane
  Point norm = tan1.cross(tan2);
  double tol = 0.1;
  double ang = tan1.angle(tan2); 
  if (norm.length() > tol)
    {
      norm.normalize();

      Point vec = norm.cross(tan1);
      if (tan2*vec >= 0.0)
	ang = M_PI - ang;
      else
	ang = M_PI + ang;
    }
  else
    {
      if (tan1*tan2 < 0.0)
	ang = M_PI - ang;
      else
	ang = M_PI + ang;
    }
     
  return ang;
}

//===========================================================================
bool compare_angle(pair<shared_ptr<Vertex>,pair<Point,double> > f1, 
		     pair<shared_ptr<Vertex>,pair<Point,double> > f2)
{
  return (f1.second.second < f2.second.second);
}
//===========================================================================
void CompleteEdgeNet::addRemainingEdges()
//===========================================================================
{
  double ptol = 1.0e-4;

  // Fetch all vertices
  vector<shared_ptr<Vertex> > vx;
  model_->getAllVertices(vx);
  Body *bd = model_->getFace(0)->getBody();
  shared_ptr<Body> bd2;
  if (!bd)
    {
      bd2 = shared_ptr<Body>(new Body(model_));
      bd = bd2.get();
    }

  // Pick verticesnot lying in a corner, and compute the opening
  // angle of the normal vector in corners
  size_t ki, kj, kh;
  vector<pair<shared_ptr<Vertex>,pair<Point,double> > > corners;
  double ang;
  Point centre;
#ifdef DEBUG
  std::ofstream of2("corners_init.g2");
#endif
  for (ki=0; ki<vx.size(); ++ki)
    {
#ifdef DEBUG
      of2 << "400 1 0 4 0 255 0 255" << std::endl;
      of2 << 1 << std::endl;
      of2 << vx[ki]->getVertexPoint() << std::endl;
#endif

      bool in_corner = vertexInfo(vx[ki], ang, centre);
      if (in_corner)
	corners.push_back(make_pair(vx[ki],make_pair(centre,ang)));
    }

  if (corners.size() <= 8)
    return;  // No further splitting is required

  // Remove vertices where a missing edge already are identified
  int kr;
  for (kr=0; kr<(int)corners.size(); ++kr)
    {
      for (kj=0; kj<missing_edges_.size(); ++kj)
	{
	  if (missing_edges_[kj].first.get() == corners[kr].first.get() ||
	      missing_edges_[kj].second.get() == corners[kr].first.get())
	    {
	      corners.erase(corners.begin()+kr);
	      kr--;
	      break;
	    }
	}
    }

  if (corners.size() <= 8)
    return;  // No further splitting is required (or don't know how
  // to handle this)

  // Sort the corner vertices according to opening angle and
  // select to split in the most convex vertices, leaving at least 8 vertices
  std::sort(corners.begin(), corners.end(), compare_angle);
  
#ifdef DEBUG
  std::ofstream of0("corners0.g2");
  of0 << "400 1 0 4 255 0 0 255" << std::endl;
  of0 << corners.size() << std::endl;
  for (ki=0; ki<corners.size(); ++ki)
  of0 << corners[ki].first->getVertexPoint() << std::endl;
#endif

  // Count number of convex vertices
  for (ki=0; ki<corners.size(); ++ki)
    if (corners[ki].second.second >= M_PI)
      break;
  int nmb = std::max((int)(corners.size()-ki), 8);
  corners.erase(corners.begin()+corners.size()-nmb, corners.end());

  // Remove the selected corner vertices and the vertices next to them
  // from the vertex pool
  // Remove also vertices with only two connected edges
  for (ki=0; ki<corners.size(); ++ki)
    {
      for (kj=0; kj<vx.size();)
	{
	  if (corners[ki].first.get() == vx[kj].get() ||
	      corners[ki].first->sameEdge(vx[kj].get()) ||
	      vx[kj]->nmbUniqueEdges() < 3)
	    vx.erase(vx.begin()+kj);
	  else
	    kj++;
	}
    }
	  
#ifdef DEBUG
  std::ofstream of("corners.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << corners.size() << std::endl;
  for (ki=0; ki<corners.size(); ++ki)
  of << corners[ki].first->getVertexPoint() << std::endl;
#endif

  // For each remaining corner, define a missing edge between this 
  // corner and the closest vertex in the pool. One pool vertex can
  // only be part of one missing edge
  // !!! This may be a too simple solution in the longer run
  
  // Select first corner to connect to
  vector<double> acc_dist(corners.size());
  for (kh=0; kh<corners.size(); ++kh)
    {
      vector<pair<shared_ptr<Vertex>,pair<Point,double> > > corners2(corners.begin(),
								     corners.end());
      vector<shared_ptr<Vertex> > vx2(vx.begin(), vx.end());
      vector<double> curr_dist(corners.size());
      acc_dist[kh] = 0.0;
      for (ki=0; ki<corners2.size(); ++ki)
	{
	  size_t ki2 = (kh + ki)%corners2.size();
	  Point pnt1 = corners2[ki2].first->getVertexPoint();
	  double mindist = HUGE;
	  int minind = -1;
	  for (kj=0; kj<vx2.size(); ++kj)
	    {
	      double dist = corners2[ki2].first->getDist(vx2[kj]);
	      if (dist < mindist)
		{
		  // Check if the new edge goes into the material
		  Point pnt2 = vx2[kj]->getVertexPoint();
		  Point pnt3 = pnt1 + 0.1*(pnt2 - pnt1);
		  Point pnt4 = pnt1 - 0.05*(pnt2 - pnt1);
		  
		  if (bd->isInside(pnt3) && (!bd->isInside(pnt4)))
		    {
		      mindist = dist;
		      minind = (int)kj;
		    }
		}
	    }
	  if (minind < 0)
	    {
	      std::cout << "Negative index in add missing edges" << std::endl;
	    }
	  acc_dist[kh] += mindist;
	  if (minind >= 0)
	    vx2.erase(vx2.begin()+minind);
	}
    }

  // Do the actual connection
  double mind = HUGE;
  int mincorner = -1;
  for (kh=0; kh<acc_dist.size(); ++kh)
    {
      if (acc_dist[kh] < mind)
	{
	  mind = acc_dist[kh];
	  mincorner = kh;
	}
    }

  for (kh=0; kh<corners.size(); ++kh)
    {
      ki = (kh + mincorner)%corners.size();

      double mindist = HUGE;
      int minind = -1;
      for (kj=0; kj<vx.size(); ++kj)
	{
	  double dist = corners[ki].first->getDist(vx[kj]);
	  if (dist < mindist)
	    {
	      // Check if the new edge goes into the material
	      Point pnt1 = corners[ki].first->getVertexPoint();
	      Point pnt2 = vx[kj]->getVertexPoint();
	      Point pnt3 = pnt1 + 0.1*(pnt2 - pnt1);
	      Point pnt4 = pnt1 - 0.05*(pnt2 - pnt1);
		  
	      if (bd->isInside(pnt3) && (!bd->isInside(pnt4)))
		{
		  mindist = dist;
		  minind = (int)kj;
		}
	    }
	}
      // Connect
      // Make sure not to create a 3-sided surface
      bool connected = corners[ki].first->connectedToSameVertex(vx[minind].get());
      if (!connected)
	missing_edges_.push_back(make_pair(corners[ki].first, vx[minind]));
    
#ifdef DEBUG
      of << "410 1 0 4 100 55 100 255 " << std::endl;
      of << "1" << std::endl;
      of << corners[ki].first->getVertexPoint() << "  ";
      of << vx[minind]->getVertexPoint() << std::endl;
#endif
      vx.erase(vx.begin()+minind);
    }
      
//   while (corners.size() > 0)
//     {
//       double mindist = HUGE;
//       int minind = -1;
//       int min_corner = -1;
//       for (ki=0; ki<corners.size(); ++ki)
// 	{
// 	  for (kj=0; kj<vx.size(); ++kj)
// 	    {
// 	      double dist = corners[ki].first->getDist(vx[kj]);
// 	      if (dist < mindist)
// 		{
// 		  // Check if the new edge goes into the material
// 		  Point pnt1 = corners[ki].first->getVertexPoint();
// 		  Point pnt2 = vx[kj]->getVertexPoint();
// 		  Point pnt3 = pnt1 + 0.1*(pnt2 - pnt1);
// 		  Point pnt4 = pnt1 - 0.05*(pnt2 - pnt1);
		
// 		  if (bd->isInside(pnt3) && (!bd->isInside(pnt4)))
// 		    {
// 		      mindist = dist;
// 		      minind = (int)kj;
// 		      min_corner = (int)ki;
// 		    }
// 		}
// 	    }
// 	}

//       // Connect
//       missing_edges_.push_back(make_pair(corners[min_corner].first, vx[minind]));
    
// #ifdef DEBUG
//       of << "410 1 0 4 100 55 100 255 " << std::endl;
//       of << "1" << std::endl;
//       of << corners[min_corner].first->getVertexPoint() << "  ";
//       of << vx[minind]->getVertexPoint() << std::endl;
// #endif
//       vx.erase(vx.begin()+minind);
//       corners.erase(corners.begin()+min_corner);
//     }
}	  

//===========================================================================
bool CompleteEdgeNet::vertexInfo(shared_ptr<Vertex> vx, double& angle,
				 Point& centre)
//===========================================================================
{
  double tol = 1.0e-6; //1.0e-9;

  // Fetch all faces meeting in this vertex and belonging to this body
  Body *bd = model_->getBody();
  vector<pair<ftSurface*,Point> > faces = vx->getFaces(bd);
  if (faces.size() < 3)
    return false;

  // Compute the surface normal corresponding to all faces
  vector<Point> norm(faces.size());
  size_t ki, kj;
  for (ki=0; ki<faces.size(); ++ki)
    norm[ki] = faces[ki].first->normal(faces[ki].second[0],
				       faces[ki].second[1]);

  // Compute the vector cone corresponging to these normals
  DirectionCone cone(norm[0]);
  for (size_t ki=1; ki<norm.size(); ++ki)
    cone.addUnionWith(norm[ki]);

  // Check if the normal vector span a volume
  Point vec = norm[0].cross(norm[1]);
  double ang = norm[0].angle(norm[1]);
  double angtol = model_->getTolerances().bend;
  double eps = model_->getTolerances().gap;
  for (ki=2; ki<norm.size(); ++ki)
    {
      if (ang > angtol)
	break;
      vec = norm[0].cross(norm[ki]);
      ang = norm[0].angle(norm[ki]);
    }
  ki--;

  for (kj=ki+1; kj<norm.size(); ++kj)
    {
      double ang1 = norm[0].angle(norm[kj]);
      double ang2 = norm[ki].angle(norm[kj]);
      if (ang1 > angtol && ang2 > angtol)
	break;
    }

  // if (kj == norm.size())
  //   {
  //     angle = M_PI;
  //     return false;  // The normal vectors do not span a volume
  //   }
      
  centre = cone.centre();
  if (cone.greaterThanPi())
    {
      angle = 1.5*M_PI;
      return true;
    }
  else
    angle = cone.angle();
  
  // Check if the corner is convex or concave
  // For each associated surface, compute the partial derivatives in
  // the vertex and project the cone centre into the tangent plane
  Point tan_pt;  // Point in the area to be used in indeterminate situations
  int sgnpluss = 0, sgnminus = 0;
  for (ki=0; ki<norm.size(); ++ki)
    {
      vector<ftEdge*> edges = vx->getFaceEdges(faces[ki].first->asFtSurface());
      if (edges.size() != 2)
	continue;
      double t1 = edges[0]->parAtVertex(vx.get());
      double t2 = edges[1]->parAtVertex(vx.get());
      Point tan1 = edges[0]->tangent(t1);
      if (fabs(edges[0]->tMax()-t1) < fabs(t1-edges[0]->tMin()))
	tan1 *= -1.0;
      Point tan2 = edges[1]->tangent(t2);
      if (fabs(edges[1]->tMax()-t2) < fabs(t2-edges[1]->tMin()))
	tan2 *= -1.0;
      tan1.normalize();
      tan2.normalize();
      vec = centre - (centre*norm[ki])*norm[ki];
      Point vec2 = 0.5*(tan1+tan2);
      double scpr = vec*vec2;
      if (scpr > tol)
	sgnpluss++;
      else if (scpr < -tol)
	sgnminus++;

      if (tan_pt.dimension() == 0 && vec2.length() > eps)
	{
	  vec2.normalize();
	  double len = edges[0]->estimatedCurveLength();
	  double dist_fac = 0.01;
	  tan_pt = vx->getVertexPoint() + dist_fac*len*vec2;
	}
    }

  if (sgnpluss > 0 && sgnminus > 0)
    angle = M_PI;  // A saddle point
  else if (sgnminus > 0)
    angle = 2*M_PI - angle;  // A convex corner
  else if (sgnminus == 0 && sgnpluss == 0)
    {
      // Extra test in indeterminate situations. Check if the
      // point in the tangent plane lies inside the model. In
      // that case it is a concave corner
      if (bd == 0 || tan_pt.dimension() == 0)
	angle = 2*M_PI - angle;  // Doesn't really make sense
      else
	{
	  bool inside = bd->isInside(tan_pt);
	  if (!inside)
	    angle = 2*M_PI - angle;  // A convex corner
	}
    }
      

  return true;
}

//===========================================================================
// 
// 
vector<ftEdge*> CompleteEdgeNet::getStartEdges()
//===========================================================================
{
  // Fetch all unique edges. First fetch shells
  size_t ki;
  vector<shared_ptr<ftEdge> > edges = 
    model_->getUniqueInnerEdges();

  // Keep the edges that split one underlying surface
  vector<ftEdge*> start_edges;
  for (ki=0; ki<edges.size(); ++ki)
    {
      if (!edges[ki]->twin())
	continue; // Unexpected for a solid, but definitely not a splitting edge

      if (!edges[ki]->face())
	continue; // Newly constructed edge, treated elsewhere

      shared_ptr<ParamSurface> sf1 = 
	edges[ki]->face()->asFtSurface()->surface();
      shared_ptr<ParamSurface> sf2 = 
	edges[ki]->twin()->geomEdge()->face()->asFtSurface()->surface();

      // Check for trimmed surfaces
      shared_ptr<BoundedSurface> bd_sf1 = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf1);
      shared_ptr<BoundedSurface> bd_sf2 = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf2);
      if (!(bd_sf1.get() && bd_sf2.get()))
	continue;
      if (bd_sf1->underlyingSurface().get() ==
	  bd_sf2->underlyingSurface().get())
	start_edges.push_back(edges[ki].get());  // A start edge is found
      else
	{
	  // An edge may be relevant for a loop around a missing surface
	  // even if the underlying surface on both sides are different.
	  // This can be due to merging of surfaces across a seam
	  // Count the number of bodies meeting in this edge, and the number of
	  // faces for each body
	  // First check continuity in the current body
	  shared_ptr<Vertex> vx1 = edges[ki]->getVertex(true);
	  shared_ptr<Vertex> vx2 = edges[ki]->getVertex(false);
	  Point norm1 = edges[ki]->normal(edges[ki]->tMin());
	  double t1 = edges[ki]->twin()->geomEdge()->parAtVertex(vx1.get());
	  Point norm2 = edges[ki]->twin()->normal(t1);
	  Point norm3 = edges[ki]->normal(edges[ki]->tMax());
	  double t2 = edges[ki]->twin()->geomEdge()->parAtVertex(vx2.get());
	  Point norm4 = edges[ki]->twin()->normal(t2);
	  double angtol = model_->getTolerances().bend;
	  if (norm1.angle(norm2) < angtol && norm3.angle(norm4) < angtol)
	    {
	      // Sufficient continuity. Check configuration
	      vector<ftSurface*> faces = edges[ki]->getAllAdjacentFaces();

	      // Check equality of faces
	      size_t kr, kh;
	      for (kr=0; kr<faces.size(); ++kr)
		  for (kh=kr+1; kh<faces.size(); )
		    {
		      if (faces[kr] == faces[kh])
			faces.erase(faces.begin()+kh);
		      else
			kh++;
		    }

	      vector<Body*> bodies;
	      vector<int> nmb_faces;
	      for (size_t kj=0; kj<faces.size(); ++kj)
		{
		  Body *bd = faces[kj]->getBody();
		  for (kr=0; kr<bodies.size(); ++kr)
		    {
		      if (bd == bodies[kr])
			break;
		    }
		  if (kr == bodies.size())
		    {
		      bodies.push_back(bd);
		      nmb_faces.push_back(1);
		    }
		  else 
		    nmb_faces[kr]++;
		}
	  
	      // Check if the number of associated faces for eah body is always two
	      for (kr=0; kr<nmb_faces.size(); ++kr)
		if (nmb_faces[kr] != 2)
		  break;

	      if ((bodies.size() == 1 || bodies.size() == 2) &&
		  kr == nmb_faces.size())
		start_edges.push_back(edges[ki].get());  // A start edge is found
	    }
	}
    }

  return start_edges;
}

//===========================================================================
bool CompleteEdgeNet::betterConnectionInFace(shared_ptr<Vertex> source,
					     Body *bd,
					     ftEdge* edg1, ftEdge *edg2,
					     shared_ptr<Vertex> dest)
//===========================================================================
{
  // Compute tangent at source
  double par1 = edg1->parAtVertex(source.get());
  double par2 = edg2->parAtVertex(source.get());
  Point tan1 = edg1->tangent(par1);
  Point tan2 = edg2->tangent(par2);
  if (tan1*tan2 < 0.0)
    tan2 *= -1;
  Point tan = 0.5*(tan1 + tan2);

  // Fetch all vertices in all faces belonging to the destination vertex
  vector<ftSurface*> faces = dest->faces(bd);
  std::set<shared_ptr<Vertex> > all_vertices;  // All vertices represented once
  size_t ki;
  for (ki=0; ki<faces.size(); ++ki)
    {
      vector<shared_ptr<Vertex> > curr_vertices = 
	faces[ki]->vertices();
      all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
    }
  vector<shared_ptr<Vertex> > vertices(all_vertices.begin(), all_vertices.end());

  // Parameters to compare with
  Point vx1 = source->getVertexPoint();
  Point vx2 = dest->getVertexPoint();
  Point vec = vx2 - vx1;
  double angle = tan.angle(vec);
  double len = vec.length();

  int idx = -1;
  double fac = 0.9;
  double pihalf = 0.5*M_PI;
  for (ki=0; ki<vertices.size(); ++ki)
    {
      Point v1 = vertices[ki]->getVertexPoint();
      double ang1 = tan.angle(v1 - vx1);
      double l1 = vx1.dist(v1);
      if (fabs(pihalf - ang1) < fac*fabs(pihalf - angle) && l1 < fac*len)
	idx = (int)ki;
    }

  if (idx >= 0)
    return true;
  else
    return false;
}

//===========================================================================
void CompleteEdgeNet::writePath(vector<ftEdge*>& edges,
				shared_ptr<Vertex> vx)
//===========================================================================
{
#ifdef DEBUG
  std::ofstream of("path.g2");
  Point pt1 = edges[0]->point(edges[0]->tMin());
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pt1 << std::endl;

  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      shared_ptr<ParamCurve> curr_crv = edges[ki]->geomCurve();
      shared_ptr<ParamCurve> curr_crv2 = 
	shared_ptr<ParamCurve>(curr_crv->subCurve(edges[ki]->tMin(),
						  edges[ki]->tMax()));
      shared_ptr<CurveOnSurface> sf_cv = 
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curr_crv2);
      if (sf_cv.get())
	{
	  sf_cv->spaceCurve()->writeStandardHeader(of);
	  sf_cv->spaceCurve()->write(of);
	}
      else
	{
	  curr_crv2->writeStandardHeader(of);
	  curr_crv2->write(of);
	}

      Point pnt = edges[ki]->point(edges[ki]->tMax());
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << "1" << std::endl;
      of << pnt << std::endl;
    }

  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << vx->getVertexPoint() << std::endl;
#endif
}

