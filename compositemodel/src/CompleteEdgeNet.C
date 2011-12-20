//===========================================================================
//                                                                           
// File: CompleteEdgeNet.C                                                   
//                                                                           
// Created: May. 2011                                         
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/CompleteEdgeNet.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/compositemodel/Path.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/geometry/BoundedSurface.h"
#include <fstream>

#define DEBUG

using std::vector;
using std::make_pair;
using namespace Go;

//===========================================================================
CompleteEdgeNet::CompleteEdgeNet(shared_ptr<SurfaceModel> sfmodel)
//===========================================================================
  : model_(sfmodel)
{
}

//===========================================================================
CompleteEdgeNet::~CompleteEdgeNet()
//===========================================================================
{
}

//===========================================================================
bool CompleteEdgeNet::perform()
//===========================================================================
{
  // Check if the surface model represents a solid
  int nmb_bd = model_->nmbBoundaries();
  if (nmb_bd > 0)
    return false;  // Not a solid

  // Prepare the model for edge completion
  RegularizeFaceSet regularize(model_);
  model_ = regularize.getRegularModel();

  addMissingEdges();
  return true;
}


//===========================================================================
void CompleteEdgeNet::addMissingEdges()
//===========================================================================
{
  // Fetch all edges in the surface model, twins are represented twice
  vector<ftEdge*> edges;
  int nmb_faces = model_->nmbEntities();
  int ki;
  for (ki=0; ki<nmb_faces; ++ki)
    {
      vector<ftEdge*> curr_edges = model_->getFace(ki)->getAllEdgePtrs();
      edges.insert(edges.end(), curr_edges.begin(), curr_edges.end());
    }

  // Traverse edges to create loops
  while (edges.size() > 0)
    {
      ftEdge *curr_edge = edges[0];
      shared_ptr<Vertex> vx = curr_edge->getVertex(false);

      vector<ftEdge*> curr_path;
      traverseEdges(edges, curr_path, curr_edge, vx, false);
    }

  addRemainingEdges();
}

//===========================================================================
void CompleteEdgeNet::traverseEdges(vector<ftEdge*>& edges,
				    vector<ftEdge*>& curr_path,
				    ftEdge *curr_edge,
				    shared_ptr<Vertex> vx,
				    bool search_end)
//===========================================================================
{
  // Move edge (and twin) from list of allowed edges to current path
  curr_path.push_back(curr_edge);
  ftEdge *twin = curr_edge->twin()->geomEdge();

  writePath(curr_path, vx);

  vector<ftEdge*>::iterator e1 = std::find(edges.begin(), edges.end(),
					   curr_edge);
  if (e1 != edges.end())
    edges.erase(e1);
  // e1 = std::find(edges.begin(), edges.end(), twin);
  // if (e1 != edges.end())
  //   edges.erase(e1);

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
      regularizeEdgeLoop(loop);

      // Remove loop from path
      e1 = std::find(curr_path.begin(), curr_path.end(), loop[0]);
      vector<ftEdge*>::iterator e2 = 
	std::find(curr_path.begin(), curr_path.end(), loop[loop.size()-1]);
      e2++;
      curr_path.erase(e1, e2);
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
	      regularizeEdgeLoop(curr_path);
	      curr_path.clear();
	      curr_edge = 0;
	      vx.reset();
	    }
	  else if (search_end || kj < missing_edges_.size())
	    {
	      // Will not perform a legal loop. Step one edge back
	      // Try to continue along anothe edge
	      writePath(curr_path, vx);

	      int stop_break0 = 1.0;
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
      next_edge = fetchNextEdge(curr_edge, vx, edges);
      if (!next_edge)
	{
	  writePath(curr_path, vx);
	  curr_path.pop_back();
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
      traverseEdges(edges, curr_path, next_edge, vx, search_end);
    }
}

//===========================================================================
ftEdge* CompleteEdgeNet::fetchNextEdge(ftEdge *curr_edge,
				       shared_ptr<Vertex> vx,
				       vector<ftEdge*> edges)
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
  vector<ftEdge*> vx_edges = vx->uniqueEdges();
  size_t ki;
  for (ki=0; ki<vx_edges.size(); ++ki)
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

      // Check if the edge is legal
      ftEdge *twin = vx_edges[ki]->twin()->geomEdge();
      vector<ftEdge*>::iterator e1 = std::find(edges.begin(), edges.end(),
					       vx_edges[ki]);      
      if (e1 != edges.end())
	break; // An edge is found

      vector<ftEdge*>::iterator e2 = std::find(edges.begin(), edges.end(),
					       twin);      
      
      if (e2 != edges.end())
	{
	  vx_edges[ki] = twin;
	  break; // An edge is found
	}
    }
     
  if (ki < vx_edges.size())
    next = vx_edges[ki];

  return next;
}

//===========================================================================
void CompleteEdgeNet::regularizeEdgeLoop(vector<ftEdge*> edges)
//===========================================================================
{

  bool to_add_edges = false;

  // Fetch vertices
  vector<shared_ptr<Vertex> > vxs;
  int ki, kj, kr, kh;
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

  // Set index pointers to significant vertices (more than two edges meet)
  vector<int> sign;
  for (ki=0; ki<(int)vxs.size(); ++ki)
    if (vxs[ki]->nmbUniqueEdges() > 2)
      sign.push_back((int)ki);

  if (vxs.size() <= 4)
    return;  // No need for extra edges

  // Check for repeated vertices
  for (ki=0; ki<(int)vxs.size(); ++ki)
    for (kj=ki+1; kj<(int)vxs.size(); ++kj)
      if (vxs[ki].get() == vxs[kj].get())
	return;

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

  // Check if the plane(s) defined by the edge loop are significantly
  // different from the tangent planes of the assiciated faces
  // To check if the loop should be prosessed further
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
	  if (tan1.angle(tan2) < model_->getTolerances().bend)
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
	  if (ang < model_->getTolerances().bend)
	    nmb_plane++;
	}

      std::cout << "Nmb check: " << nmb_check << ", nmb plane: " << nmb_plane << std::endl;

      //if (nmb_plane < (int)(0.5*nmb_check + 1))
      if (nmb_plane < (int)(0.9*nmb_check + 1))
	to_add_edges = true;
    }

  if (!to_add_edges)
    return;

  // Select vertices between which to add edges
  // First project all vertices into the best plane given by the vertices
  vector<Point> vxs2(vxs.size());
  Point norm(0.0, 0.0, 0.0);
  int nmb_vec = 0;
  for (ki=0; ki<(int)vxs.size(); ++ki)
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
  norm /= (double)nmb_vec;

  for (ki=0; ki<(int)vxs.size(); ++ki)
    {
      Point vxpt = vxs[ki]->getVertexPoint();
      vxs2[ki] = vxpt - (vxpt*norm)*norm;
    }

  // Select start vertex for the construction of missing edges
  int vx_idx = -1;
  int vx_idx2 = -1;
  double min_ang = MAXDOUBLE;
  double min_dist = MAXDOUBLE;
  double ang_tol = model_->getTolerances().kink;
  double dist_fac = 10.0;
  if (vxs.size() > edges.size())
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

	  // Check if the candidate crosses the polygon of projected
	  // vertices
	  size_t kk;
	  double tol = 1.0e-8;
	  Point n1 = vec1.cross(norm);
	  Point p1 = 0.5*(vxs2[ki] + vxs2[kh]);
	  for (kk=1; kk<vxs2.size(); ++kk)
	    {
	      Point n2 = (vxs2[kk] - vxs2[kk-1]).cross(norm);
	      Point p2 = 0.5*(vxs2[kk-1] + vxs2[kk]);
	      double t1 = ((p2 - vxs2[ki])*n2)/(vec1*n2);	
	      double s1 = ((p1 - vxs2[kk-1])*n1)/((vxs2[kk]-vxs2[kk-1])*n1);
	      if (t1 > tol && t1 < 1.0-tol && s1 < 1.0-tol)
		break;
	    }
	  
	  if (kk <vxs2.size())
	    continue;

	  Point n2 = (vxs2[0] - vxs2[vxs2.size()-1]).cross(norm);
	  Point p2 = 0.5*(vxs2[vxs2.size()-1] + vxs2[0]);
	  double t1 = ((p2 - vxs2[ki])*n2)/(vec1*n2);	
	  double s1 = ((p1 - vxs2[vxs2.size()-1])*n1)/((vxs2[0]-vxs2[vxs2.size()-1])*n1);
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
	      dist < dist_fac*min_dist && ang*dist < fac*min_ang*min_dist)
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
    return;  // No combination is found

#ifdef DEBUG
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << vxs[vx_idx]->getVertexPoint() << std::endl;
#endif

  // Combine vertices
  if (vx_idx2 < 0)
      vx_idx2 = (vx_idx+3) % (int)vxs2.size();
  for (ki=vx_idx, kj=vx_idx2; abs(int(kj-ki))>1;
       ki=(ki>0)?ki-1:(int)vxs2.size()-1, kj=(kj+1)%(int)vxs2.size())
    {
	if ((ki==0 && kj==(int)vxs2.size()-1) ||
	    (kj==0 && ki==(int)vxs2.size()) ||
	  abs(int(ki-kj)) <= 1)
	break;

      if (vxs[ki]->sameEdge(vxs[kj].get()))
	  continue;  // Already connected

      // Check if the missing edge exists already
      size_t kn;
      for (kn=0; kn<missing_edges_.size(); ++kn)
	if ((missing_edges_[kn].first.get() == vxs[ki].get() && 
	     missing_edges_[kn].second.get() == vxs[kj].get()) ||
	    (missing_edges_[kn].first.get() == vxs[kj].get() && 
	     missing_edges_[kn].second.get() == vxs[ki].get()))
	  break;

      if (kn == missing_edges_.size())
	missing_edges_.push_back(make_pair(vxs[ki], vxs[kj]));
#ifdef DEBUG      
      of << "410 1 0 4 255 0 0 255" << std::endl;
      of << "1" << std::endl;
      of << vxs[ki]->getVertexPoint() << " " << vxs[kj]->getVertexPoint() << std::endl;
#endif
    }
      

  // // Set index pointers to the four most convex corners
  // vector<int> corner;
  // vector<double> angle;
  // for (ki=0; ki<sign.size(); ++ki)
  //   {
  //     // Compute angle
  //     int idx = (sign[ki] == 0) ? sign[sign.size()-1] : sign[ki] - 1;
  //     double ang = getVertexAngle(edges[idx], edges[sign[ki]]);
  //     if (corner.size() < 4)
  // 	{
  // 	  angle.push_back(ang);
  // 	  corner.push_back(sign[ki]);
  // 	}
  //     else
  // 	{
  // 	  // Replace least significant corner
  // 	  size_t max_idx = angle.size();
  // 	  size_t max_ang = ang;
  // 	  for (size_t kr=0; kr<angle.size(); ++kr)
  // 	    if (angle[kr] < max_ang)
  // 	      {
  // 		max_ang = angle[kr];
  // 		max_idx = kr;
  // 	      }

  // 	  if (max_idx < angle.size())
  // 	    {
  // 	      angle.erase(angle.begin()+max_idx);
  // 	      corner.erase(corner.begin()+max_idx);
  // 	      angle.push_back(ang);
  // 	      corner.push_back(sign[ki]);
  // 	    }
  // 	}
  //   }	  

  // for (ki=0; ki<corner.size(); ++ki)
  //   {
  //     of << "400 1 0 4 255 0 0 255" << std::endl;
  //     of << "1" << std::endl;
  //     of << vxs[corner[ki]]->getVertexPoint() << std::endl;
  //   }
  int stop_break;
  stop_break = 1;
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
bool compare_angle(pair<shared_ptr<Vertex>,double> f1, 
		   pair<shared_ptr<Vertex>,double> f2)
{
  return (f1.second < f2.second);
}
//===========================================================================
void CompleteEdgeNet::addRemainingEdges()
//===========================================================================
{
  // Fetch all vertices
  vector<shared_ptr<Vertex> > vx;
  model_->getAllVertices(vx);

  // Pick verticesnot lying in a corner, and compute the opening
  // angle of the normal vector in corners
  size_t ki, kj;
  vector<pair<shared_ptr<Vertex>,double> > corners;
  double ang;
  for (ki=0; ki<vx.size(); ++ki)
    {
      bool in_corner = vertexInfo(vx[ki], ang);
      if (in_corner)
	corners.push_back(make_pair(vx[ki],ang));
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
  
  // Count number of convex vertices
  for (ki=0; ki<corners.size(); ++ki)
    if (corners[ki].second >= M_PI)
      break;
  int nmb = std::max((int)(corners.size()-ki), 8);
  corners.erase(corners.begin()+corners.size()-nmb, corners.end());

  // Remove the selected corner vertices and the vertices next to them
  // from the vertex pool
  for (ki=0; ki<corners.size(); ++ki)
    {
      for (kj=0; kj<vx.size();)
	{
	  if (corners[ki].first.get() == vx[kj].get() ||
	      corners[ki].first->sameEdge(vx[kj].get()))
	    vx.erase(vx.begin()+kj);
	  else
	    kj++;
	}
    }
	  
  // For each remaining corner, define a missing edge between this 
  // corner and the closest vertex in the pool. One pool vertex can
  // only be part of one missing edge
  // !!! This may be a too simple solution in the longer run
  for (ki=0; ki<corners.size(); ++ki)
    {
      double mindist = HUGE;
      int minind = -1;
      for (kj=0; kj<vx.size(); ++kj)
	{
	  double dist = corners[ki].first->getDist(vx[kj]);
	  if (dist < mindist)
	    {
	      mindist = dist;
	      minind = (int)kj;
	    }
	}

      // Connect
      missing_edges_.push_back(make_pair(corners[ki].first, vx[minind]));
      vx.erase(vx.begin()+minind);
    }
}	  

//===========================================================================
bool CompleteEdgeNet::vertexInfo(shared_ptr<Vertex> vx, double& angle)
//===========================================================================
{
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

  if (kj == norm.size())
    {
      angle = M_PI;
      return false;  // The normal vectors do not span a volume
    }
      
  if (cone.greaterThanPi())
    angle = 1.5*M_PI;
  else
    angle = cone.angle();
  
  // Check if the corner is convex or concave
  // For each associated surface, compute the partial derivatives in
  // the vertex and project the cone centre into the tangent plane
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
      Point centre = cone.centre();
      vec = centre - (centre*norm[ki])*norm[ki];
      Point vec2 = 0.5*(tan1+tan2);
      double scpr = vec*vec2;
      if (scpr > 0.0)
	sgnpluss++;
      else if (scpr < 0.0)
	sgnminus++;
    }

  if (sgnpluss > 0 && sgnminus > 0)
    angle = M_PI;  // A saddle point
  else if (sgnminus > 0)
    angle = 2*M_PI - angle;  // A convex corner
      

  return true;
}

//===========================================================================
void CompleteEdgeNet::writePath(vector<ftEdge*> edges,
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

