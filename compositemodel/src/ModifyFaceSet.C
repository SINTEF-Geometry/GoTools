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

#include "GoTools/compositemodel/RegularizeFace.h"
#include "GoTools/compositemodel/ModifyFaceSet.h"
#include "GoTools/compositemodel/RegularizeUtils.h"
#include "GoTools/compositemodel/ftEdge.h"
#include <fstream>
#include <cstdlib>

//#define DEBUG

using std::vector;
using std::set;
using std::make_pair;

using namespace Go;

//==========================================================================
    ModifyFaceSet::ModifyFaceSet(shared_ptr<SurfaceModel> model)
      : model_(model)
//==========================================================================
{
}

//==========================================================================
ModifyFaceSet::~ModifyFaceSet()
//==========================================================================
{

}

//==========================================================================
shared_ptr<SurfaceModel> ModifyFaceSet::getModifiedModel(int& nmb)
//==========================================================================
{
  nmb = divide();
  return model_;
}

  // Collect all faces
//==========================================================================
int ModifyFaceSet::divide()
//==========================================================================
{
  // Collect all faces
  vector<shared_ptr<ftSurface> > faces = model_->allFaces();
  int nmb_faces = (int)faces.size();

  tpTolerances tptol = model_->getTolerances();

  // Fetch all sharp edges
  vector<ftEdge*> sharp_edges = fetchSharpEdges();

  int nmb = 0;
  vector<shared_ptr<ftSurface> > split_faces;
  for (size_t ki=0; ki<sharp_edges.size(); ++ki)
    {
      nmb++;
      vector<shared_ptr<Vertex> > vxs;
      shared_ptr<Vertex> vx1, vx2;
      sharp_edges[ki]->getVertices(vx1, vx2);
      vxs.push_back(vx1);
      vxs.push_back(vx2);
      vector<ftEdge*> curr_edges;;
      curr_edges.push_back(sharp_edges[ki]);
      curr_edges.push_back(sharp_edges[ki]);
      vector<ftSurface*> adj_face;

      // Split faces until a full edge loop including the sharp edge
      // is obtained. The computatition is not expected to involve all
      // faces, but the loop is limited to avoid to go infinite
      int nmb_split = 0;
      bool finish = false;
      vector<ftEdge*> next_edge;
      vector<double> angle;
      for (size_t kj=0; kj<curr_edges.size(); ++kj)
	{
	  double ang2;
	  ftEdge* edge2;
	  ftSurface* face2 = fetchNextFace(curr_edges[kj], vxs[kj].get(),
					  tptol.bend, edge2, ang2);
	  next_edge.push_back(edge2);
	  adj_face.push_back(face2);
	  angle.push_back(ang2);
	}

      if (next_edge.size() < 2 ||
	  (next_edge[0] == NULL && next_edge[1] == NULL && nmb == 1))
	break;  // The situation is probably simple enough to ommit
      // splitting in concave edges

	  // bool is_concave[2];
	  // for (int kj=0; kj<2; ++kj)
	  //   is_concave[kj] = vxs[kj]->isConcave(adj_face[kj], tptol.bend);

      if (adj_face[0] == adj_face[1])
	{
	  adj_face.pop_back();
	  next_edge.pop_back();
	  angle.pop_back();
	  finish = true;
	}
	  // if (adj_face[0] == adj_face[1])
	  //   {
	  //     if (is_concave[0])
	  // 	adj_face[0] = NULL;
	  //     else
	  // 	adj_face[1] = NULL;   // Same face, split only once
	  //     finish = true;
	  //   }

	  // bool divided = false;
	  // int nmb_try = 0;
	  // while (!divided && nmb_try < 2)
	  //   {
	  //     nmb_try++;
	  //for (int kj=0; kj<2; ++kj)
      for (size_t kj=0; kj<adj_face.size();)
	{
	  ++nmb_split;
	  if (nmb_split == nmb_faces)
	    break;  // Probably an infinite loop

	  // if (/*nmb_try == 1 && */is_concave[kj] && (!is_concave[1-kj]))
	  //   continue;  // Ambigous situation, postpone if possible
	  
	  int next_ix;
	  if (adj_face[kj] && (!next_edge[kj])/*M_PI-angle[kj] > tptol.bend*/)
	    {
	      // Check if the current face is the result of a previous split.
	      // In that case the split counter must be decreased to set
	      // the correct expected edges in the edge loop of the splitting
	      // surface to be created corresponding to the concave edge
	      for (size_t kr=0; kr<split_faces.size(); ++kr)
		if (split_faces[kr].get() == adj_face[kj])
		  {
		    nmb--;
		    break;
		  }

	      // Regularize face
	      shared_ptr<ftSurface> curr = 
		model_->fetchAsSharedPtr(adj_face[kj]);
	      if (!curr.get())
		{
		  finish = true;
		  continue;
		}
	      RegularizeFace regularize(curr, model_);
	      vector<shared_ptr<Vertex> > vx_pri;
	      vx_pri.push_back(vxs[kj]);

	      shared_ptr<Vertex> vx1_curr = 
		curr->hasVertexPoint(vx1->getVertexPoint(), tptol.neighbour);
	      shared_ptr<Vertex> vx2_curr = 
		curr->hasVertexPoint(vx2->getVertexPoint(), tptol.neighbour);
	      if (vx1_curr && vx1_curr.get() != vxs[kj].get())
		vx_pri.push_back(vx1_curr);
	      if (vx2_curr && vx2_curr.get() != vxs[kj].get())
		vx_pri.push_back(vx2_curr);
	      regularize.setPriVx(vx_pri, false);
	      
	      vector<shared_ptr<Vertex> > corner = 
		curr->getCornerVertices(tptol.bend);
	      if (corner.size() > 4)
		regularize.setDivideInT(false);
	      if (corner.size() < 4 || vx_pri.size() > 1/*<= 4*/)
		regularize.setDegenFlag(true);
	      regularize.setMaxRec(1); // Only one split
	      vector<shared_ptr<ftSurface> > faces2 = 
		regularize.getRegularFaces();
#ifdef DEBUG
	      std::ofstream of2("post_regface.g2");
	      for (size_t kr=0; kr<faces2.size(); ++kr)
		{
		  faces2[kr]->surface()->writeStandardHeader(of2);
		  faces2[kr]->surface()->write(of2);
		}
#endif

	      // Collect information of created faces
	      split_faces.insert(split_faces.end(), 
				 faces2.begin(), faces2.end());
		  // if (faces2.size() > 1)
		  // 	divided = true;

		  // Identify edge continuation. First fetch candidates
	      vector<shared_ptr<ftEdge> > cand_edg;
	      for (size_t kr=0; kr<faces2.size(); ++kr)
		for (size_t kh=kr+1; kh<faces2.size(); ++kh)
		  {
		    vector<shared_ptr<ftEdge> > tmp_edg = 
		      faces2[kr]->getCommonEdges(faces2[kh].get());
		    for (size_t kb=0; kb<tmp_edg.size(); ++kb)
		      {
			if (curr_edges[kj]->commonVertex(tmp_edg[kb].get()))
			  cand_edg.push_back(tmp_edg[kb]);
		      }
		  }
	      if (cand_edg.size() == 1)
		next_edge[kj] = cand_edg[0].get();

	      next_ix = (int)adj_face.size();
	    }
	  else
	    next_ix = (int)kj+1;
	  
	  // Fetch the next vertex in the chain
	  if (next_edge[kj])
	    {
	      shared_ptr<Vertex> common_vx = 
		next_edge[kj]->getCommonVertex(curr_edges[kj]);
	      shared_ptr<Vertex> curr_vx =
		next_edge[kj]->getOtherVertex(common_vx.get());
	      if (curr_vx->getDist(vx1) < tptol.neighbour ||
		  curr_vx->getDist(vx2) < tptol.neighbour)
		break;
	      double ang2;
	      ftEdge *edge2;
	      ftSurface *curr_adj =
		fetchNextFace(next_edge[kj], curr_vx.get(),
			      tptol.bend, edge2, ang2);
	      size_t kr;
	      for (kr=0; kr<adj_face.size(); ++kr)
		if (curr_adj == adj_face[kr])
		  break;
	      if (kr == adj_face.size())
		{
		  vxs.insert(vxs.begin()+next_ix, curr_vx);
		  adj_face.insert(adj_face.begin()+next_ix, curr_adj);
		  next_edge.insert(next_edge.begin()+next_ix, edge2);
		  angle.insert(angle.begin()+next_ix, ang2);
		  curr_edges.insert(curr_edges.begin()+next_ix, next_edge[kj]);

		  vxs.erase(vxs.begin()+kj);
		  adj_face.erase(adj_face.begin()+kj);
		  next_edge.erase(next_edge.begin()+kj);
		  angle.erase(angle.begin()+kj);
		  curr_edges.erase(curr_edges.begin()+kj);
		}
	      else
		{
		  ++kj;
		  //finish = true;
		}
	    }
	  else
	    ++kj;
	  
	  nmb++;
	}
	  //	    }

      // if (adj_face[0] == adj_face[1])
      // 	nmb--;
      // if (finish || 
      // 	  vxs[0]->getVertexPoint().dist(vxs[1]->getVertexPoint()) < tptol.neighbour)
      // 	break;
    }

  return nmb;
}

//==========================================================================
vector<ftEdge*>  ModifyFaceSet::fetchSharpEdges()
//==========================================================================
{
  vector<ftEdge*> sharp_edges;
  model_->getCorners(sharp_edges);

  // Remove concex edges
  for (size_t ki=0; ki<sharp_edges.size(); )
    {
      if (!sharp_edges[ki]->twin())
	{
	  sharp_edges.erase(sharp_edges.begin()+ki);
	  continue;
	}
      int nmb_samples = 10;
      int nmb_concave = 0;
      int nmb_convex = 0;
      double t1 = sharp_edges[ki]->tMin();
      double t2 = sharp_edges[ki]->tMax();
      double del = (t2 - t1)/(double)(nmb_samples-1);
      double par;
      int kr;
      for (kr=0, par=t1; kr<nmb_samples; ++kr, par+=del)
	{
	  Point pos1 = sharp_edges[ki]->point(par);
	  Point norm1 = sharp_edges[ki]->normal(par);
	  Point tan1 = sharp_edges[ki]->tangent(par);
	  Point vec1 = tan1.cross(norm1);

	  double clo_t, clo_dist;
	  Point clo_pt;
	  sharp_edges[ki]->twin()->closestPoint(pos1, clo_t, clo_pt,
					       clo_dist);
	  Point norm2 = sharp_edges[ki]->twin()->normal(clo_t);
	  Point tan2 = sharp_edges[ki]->twin()->tangent(clo_t);
	  Point vec2 = tan2.cross(norm2);
	  double ang2 = vec1.angle(norm2);
	  double ang3 = vec2.angle(norm1);
	  if (ang2 + ang3 > M_PI)
	    nmb_concave++;
	  else
	    nmb_convex++;
	  int stop_break;
	}
      if (nmb_convex >= nmb_concave)
	  sharp_edges.erase(sharp_edges.begin()+ki);
      else
	++ki;
    }

#ifdef DEBUG
  std::ofstream of("convex_crvs.g2");
  for (size_t ki=0; ki<sharp_edges.size(); ++ki)
    {
      shared_ptr<ParamCurve> cv = 
	shared_ptr<ParamCurve>(sharp_edges[ki]->geomCurve()->subCurve(sharp_edges[ki]->tMin(),
								      sharp_edges[ki]->tMax()));
      cv->geometryCurve()->writeStandardHeader(of);
      cv->geometryCurve()->write(of);
    }
#endif
  return sharp_edges;
}

//==========================================================================
ftSurface*  ModifyFaceSet::fetchNextFace(ftEdge* edge, Vertex* vx, double angtol,
					 ftEdge*& next_edge, double& angle)
//==========================================================================
{
  next_edge = NULL;  // Initially
  ftFaceBase *edg_faces[2];
  edg_faces[0] = edge->face();
  edg_faces[1] = edge->twin()->face();
  vector<ftSurface*> vx_faces = vx->faces();
  Point vxpoint = vx->getVertexPoint();
  shared_ptr<Vertex> vx2 = edge->getOtherVertex(vx);
  Point vxdir = vxpoint - vx2->getVertexPoint();
  angle = M_PI;
  for (size_t kr=0; kr<vx_faces.size();)
    {
      int kh;
      for (kh=0; kh<2; ++kh)
	if (vx_faces[kr] == edg_faces[kh])
	  break;
      if (kh < 2)
	vx_faces.erase(vx_faces.begin()+kr);
      else
	++kr;
    }
  if (vx_faces.size() == 1)
    {
#ifdef DEBUG
      std::ofstream of1("adj_faces.g2");
      vx_faces[0]->surface()->writeStandardHeader(of1);
      vx_faces[0]->surface()->write(of1);
#endif

      // Check angles between associated edges
      Point tan1 = edge->tangent(edge->parAtVertex(vx));
      vector<ftEdge*> edgs = vx->getEdges(vx_faces[0]);
      vector<double> scpr(edgs.size(), 0.0);
      for (size_t kr=0; kr<edgs.size(); ++kr)
	{
	  Point tan2 = edgs[kr]->tangent(edgs[kr]->parAtVertex(vx));
	  double ang = tan1.angle(tan2);
	  angle = std::min(angle, ang);
	  if (ang < angtol || M_PI-ang < angtol)
	    {
	      angle = ang;   // Continue along edge chain without face split

	      // Find direction of edge continuation
	      shared_ptr<Vertex> vx3 = edgs[kr]->getOtherVertex(vx);
	      Point vxdir3 = vx3->getVertexPoint() - vxpoint;
	      scpr[kr] = vxdir*vxdir3;
	    }
	}
      double maxpr = 0.0;
      int ix = -1;
      for (size_t kr=0; kr<scpr.size(); ++kr)
	{
	  if (scpr[kr] > maxpr)
	    {
	      maxpr = scpr[kr];
	      ix = (int)kr;
	    }
	}
      if (ix >= 0)
	next_edge = edgs[ix];
      return vx_faces[0];
    }
  return NULL;
}


