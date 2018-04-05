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
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
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
void
ModifyFaceSet::getSplittingSurface(vector<shared_ptr<ParamSurface> >& split_sfs,
				   vector<ftSurface*>& corr_faces,
				   vector<vector<ftEdge*> >& edges)
//==========================================================================
{
  tpTolerances tptol = model_->getTolerances();

  if (model_->nmbEntities() < 7)
    return;  // No point in looking for splitting surfaces

  // Fetch all sharp edges
  vector<ftEdge*> sharp_edges = fetchSharpEdges();
#ifdef DEBUG
  std::cout << "Nmb concave edges: " << sharp_edges.size() << std::endl;
  for (size_t ka=0; ka<sharp_edges.size(); ++ka)
    std::cout << sharp_edges[ka] << " ";
  std::cout << std::endl;
#endif

  if (sharp_edges.size() == 0)
    return;

  // Sort concave edges into chains
  vector<vector<ftEdge* > > edg_chain;
  vector<ftEdge*> curr_chain(1, sharp_edges[0]);
  edg_chain.push_back(curr_chain);
  for (size_t ki=1; ki<sharp_edges.size(); ++ki)
    {
      size_t kj, kr;
      for (kj=0; kj<edg_chain.size(); ++kj)
	{
	  for (kr=0; kr<edg_chain[kj].size(); ++kr)
	    {
	      if (sharp_edges[ki]->commonVertex(edg_chain[kj][kr]))
		{
		  edg_chain[kj].push_back(sharp_edges[ki]);
		  break;
		}
	    }
	  if (kr < edg_chain[kj].size())
	    break;
	}
      if (kj == edg_chain.size())
	{
	  vector<ftEdge*> curr_chain2(1, sharp_edges[ki]);
	  edg_chain.push_back(curr_chain2);
	}
    }
	     
  // For each edge chain, extract end vertices and check if they
  // end up in a T-joint vertex. In that case other model split
  // methods are more adapted
  for (size_t kr=0; kr<edg_chain.size(); )
    {
      size_t ix1=0, ix2=0;
      shared_ptr<Vertex> v1, v2;
      edg_chain[kr][0]->getVertices(v1, v2);
      for (size_t ki=1; ki<edg_chain[kr].size(); ++ki)
	{
	  size_t kj;
	  for (kj=0; kj<ki; ++kj)
	    {
	      shared_ptr<Vertex> common_vx = 
		edg_chain[kr][kj]->getCommonVertex(edg_chain[kr][ki]);
	      if (!common_vx.get())
		continue;

	      if (common_vx.get() == v1.get())
		{
		  v1 = edg_chain[kr][ki]->getOtherVertex(common_vx.get());
		  ix1 = ki;
		}
	      else if (common_vx.get() == v2.get())
		{
		  v2 = edg_chain[kr][ki]->getOtherVertex(common_vx.get());
		  ix2 = ki;
		}
	    }
	}
#ifdef DEBUG
      std::cout << "ix1 = " << ix1 << ", ix2 = " << ix2 << std::endl;
#endif
      if (!(edg_chain[kr].size() > 1 && ix1 == ix2))
	{
	  // Check if the edge corresponds to Tjoint vertices in the next
	  // faces
	  double ang1, ang2;
	  ftEdge *edge1, *edge2;
	  ftSurface *face1 = fetchNextFace(edg_chain[kr][ix1], v1.get(),
					   tptol.bend, edge1, ang1);
	  ftSurface *face2 = fetchNextFace(edg_chain[kr][ix2], v2.get(),
					   tptol.bend, edge2, ang2);

#ifdef DEBUG
	  std::cout << "edge1: " << edge1 << ", edge2: " << edge2 << std::endl;
#endif
	  if (edge1 != NULL || edge2 != NULL)
	    {
	      // An edge continuation to the concave edge is found
	      edg_chain.erase(edg_chain.begin()+kr);  	 
	    }
	  else
	    ++kr;
	}
      else
	{
	  edg_chain.erase(edg_chain.begin()+kr);  	 
	}
    }

  if (edg_chain.size() == 0)
    return;

  // Collect adjacent faces appropriate for use as a splitting surface
  vector<ftSurface*> cand_faces;
  vector<ElementarySurface*> elem_sfs;
  vector<vector<ftEdge*> > corr_edgs;
  vector<ClassType> sfs_type;
  for (size_t ki=0; ki<sharp_edges.size(); ++ki)
    {
      ftSurface *f1 = sharp_edges[ki]->face()->asFtSurface();
      ftSurface *f2 = (sharp_edges[ki]->twin()) ? 
	sharp_edges[ki]->twin()->face()->asFtSurface() : NULL;

      // Check if face f1 is found already
      size_t kr;
      if (f1)
	{
	  for (kr=0; kr<cand_faces.size(); ++kr)
	    if (cand_faces[kr] == f1)
	      {
		corr_edgs[kr].push_back(sharp_edges[ki]);
		break;
	      }
	  if (kr == cand_faces.size())
	    {
	      shared_ptr<ParamSurface> parent1 = f1->surface()->getParentSurface();
	      ElementarySurface *elem1 = f1->surface()->elementarySurface();
	      if ((!elem1) && parent1.get())
		{
		  shared_ptr<ElementarySurface> elem01 = 
		    dynamic_pointer_cast<ElementarySurface,ParamSurface>(parent1);
		  elem1 = elem01.get();
		}
	      if (elem1)
		{
		  cand_faces.push_back(f1);
		  elem_sfs.push_back(elem1);
		  sfs_type.push_back(elem1->instanceType());
		  vector<ftEdge*> curr;
		  curr.push_back(sharp_edges[ki]);
		  corr_edgs.push_back(curr);
		}
	    }
	}

      if (f2)
	{
	  // Check if face f2 is found already
	  for (kr=0; kr<cand_faces.size(); ++kr)
	    if (cand_faces[kr] == f2)
	      {
		corr_edgs[kr].push_back(sharp_edges[ki]);
		break;
	      }
	  if (kr == cand_faces.size())
	    {
	      shared_ptr<ParamSurface> parent2 = 
		f2->surface()->getParentSurface();
	      ElementarySurface *elem2 = f2->surface()->elementarySurface();
	      if ((!elem2) && parent2.get())
		{
		  shared_ptr<ElementarySurface> elem02 = 
		    dynamic_pointer_cast<ElementarySurface,ParamSurface>(parent2);
		  elem2 = elem02.get();
		}
	      if (elem2)
		{
		  cand_faces.push_back(f2);
		  elem_sfs.push_back(elem2);
		  sfs_type.push_back(elem2->instanceType());
		  vector<ftEdge*> curr;
		  curr.push_back(sharp_edges[ki]);
		  corr_edgs.push_back(curr);
		}
	    }
	}
    }

  if (cand_faces.size() == 0)
    return;

  // Sort candidate faces according to importance. First compute some info
  vector<BoundingBox> bbox(cand_faces.size());
  //vector<DirectionCone> dcone(cand_faces.size());
  vector<double> bsize(cand_faces.size());
  for (size_t kr=0; kr<cand_faces.size(); ++kr)
    {
      shared_ptr<ParamSurface> surf = cand_faces[kr]->surface();
      bbox[kr] = surf->boundingBox();
      //dcone[ki] = surf->normalCone();
      bsize[kr] = bbox[kr].low().dist(bbox[kr].high());
    }

  for (size_t ki=0; ki<cand_faces.size(); ++ki)
    for (size_t kr=ki+1; kr<cand_faces.size(); ++kr)
      {
	if (corr_edgs[kr].size() > corr_edgs[ki].size() ||
	    (corr_edgs[kr].size() == corr_edgs[ki].size() && 
	     sfs_type[kr] < sfs_type[ki]) || 
	    (corr_edgs[kr].size() == corr_edgs[ki].size() && 
	     sfs_type[kr] == sfs_type[ki] && bsize[kr] > bsize[ki]))
	  {
	    std::swap(cand_faces[ki], cand_faces[kr]);
	    std::swap(elem_sfs[ki], elem_sfs[kr]);
	    std::swap(sfs_type[ki], sfs_type[kr]);
	    std::swap(corr_edgs[ki], corr_edgs[kr]);
	    std::swap(bbox[ki], bbox[kr]);
	    //std::swap(dcone[ki], dcone[kr]);
	    std::swap(bsize[ki], bsize[kr]);
	  }
      }

  // Check if any candidates overlap
  for (size_t ki=1; ki<cand_faces.size();)
    {
      size_t kr;
      for (kr=0; kr<corr_edgs[ki].size(); ++kr)
	{
	  ftEdge *curr = corr_edgs[ki][kr];
	  size_t kj;
	  for (kj=0; kj<ki; ++kj)
	    {
	      size_t kh;
	      for (kh=0; kh<corr_edgs[kj].size(); ++kh)
		{
		  if (corr_edgs[kj][kh] == curr)
		    break;
		}
	      if (kh < corr_edgs[kj].size())
		break;
	    }
	  if (kj < ki)
	    break;
	}
      if (kr < corr_edgs[ki].size())
	{
	  // Match found
	  cand_faces.erase(cand_faces.begin()+ki);
	  elem_sfs.erase(elem_sfs.begin()+ki);
	  corr_edgs.erase(corr_edgs.begin()+ki);
	}
      else
	++ki;
    }

  // Check if several instances of found elementary surfaces really
  // represents the same surface
  for (size_t ki=0; ki<elem_sfs.size(); ++ki)
    {
      size_t kj;
      for (kj=ki+1; kj<elem_sfs.size();)
	{
	  bool same = SurfaceModelUtils::sameElementarySurface(elem_sfs[ki],
							       elem_sfs[kj],
							       tptol.gap,
							       tptol.kink);
	  if (same)
	    {
	      elem_sfs.erase(elem_sfs.begin()+kj);
	      corr_edgs[ki].insert(corr_edgs[ki].end(), corr_edgs[kj].begin(),
				   corr_edgs[kj].end());
	      corr_edgs.erase(corr_edgs.begin()+kj);
	    }
	  else
	    ++kj;
	}
    }


  BoundingBox box = model_->boundingBox();
  double len = box.high().dist(box.low());
  for (size_t ki=0; ki<elem_sfs.size(); ++ki)
    {
      // Extend elementary surface to ensure intersection
      shared_ptr<ElementarySurface> tmp_elem(elem_sfs[ki]->clone());
      tmp_elem->enlarge(len, len, len, len);
      shared_ptr<SplineSurface> tmp_sf(tmp_elem->createSplineSurface());
      tmp_sf->setElementarySurface(tmp_elem);
      split_sfs.push_back(tmp_sf);
      // for (size_t kj=0; kj<corr_edgs[ki].size(); ++kj)
      // 	{
      // 	  shared_ptr<ParamCurve> tmp_cv = corr_edgs[ki][kj]->geomCurve();
      // 	  shared_ptr<CurveOnSurface> tmp_sfcv = 
      // 	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_cv);
      // 	  if (!tmp_sfcv.get())
      // 	    continue;
      // 	  tmp_sfcv->setUnderlyingSurface(tmp_sf);
      // 	  if (elem_sfs[ki]->instanceType() > Class_Plane)
      // 	    {
      // 	      // The edges corresponding to the splitting surface may
      // 	      // have inaccurate parameter curves
      // 	      shared_ptr<ParamCurve> tmp_par = tmp_sfcv->parameterCurve();
      // 	      tmp_sfcv->unsetParameterCurve();
      // 	      bool found = tmp_sfcv->ensureParCrvExistence(tptol.gap);
      // 	      if (!found)
      // 		tmp_sfcv->setParameterCurve(tmp_par);
      // 	    }
      // 	}

    }
#ifdef DEBUG
  std::ofstream of("split_sfs.g2");
  for (size_t ki=0; ki<split_sfs.size(); ++ki)
    {
      elem_sfs[ki]->writeStandardHeader(of);
      elem_sfs[ki]->write(of);
      split_sfs[ki]->writeStandardHeader(of);
      split_sfs[ki]->write(of);      
    }
  for (size_t ki=0; ki<corr_edgs.size(); ++ki)
    {
      for (size_t kr=0; kr<corr_edgs[ki].size(); ++kr)
	{
	  shared_ptr<SplineCurve> tmp( 
				      corr_edgs[ki][kr]->geomCurve()->geometryCurve());
	  tmp->writeStandardHeader(of);
	  tmp->write(of);
	}
    }
#endif

  corr_faces = cand_faces;
  edges = corr_edgs;
  
  return;
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
	  
	  bool concave = vxs[kj]->isConcave(adj_face[kj], tptol.bend);
	  if (concave && kj < adj_face.size() - 1)
	    {
	      // Postpone
	      size_t last = adj_face.size()-1;
	      std::swap(adj_face[kj], adj_face[last]);
	      std::swap(next_edge[kj], next_edge[last]);
	      std::swap(vxs[kj], vxs[last]);
	      std::swap(angle[kj], angle[last]);
	      std::swap(curr_edges[kj], curr_edges[last]);
	    }

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
	      vector<shared_ptr<Vertex> > corner = 
		curr->getCornerVertices(tptol.bend);
	      bool degen_flag = (corner.size() < 4 || vx_pri.size() > 1/*<= 4*/);

	      if (vx_pri.size() <= 1 && kj < vxs.size()-1)
		{
		  addPrioritizedVertex(curr, vxs[kj+1], vx_pri);
		}

	      regularize.setPriVx(vx_pri, false);
	      
	      if (corner.size() > 4)
		regularize.setDivideInT(false);
	      if (degen_flag)
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
		  nmb++;
		}
	      else
		{
		  ++kj;
		  nmb++;
		  //finish = true;
		}
	    }
	  else
	    ++kj;
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


//==========================================================================
void ModifyFaceSet::addPrioritizedVertex(shared_ptr<ftSurface> face,
					 shared_ptr<Vertex> vx,
					 vector<shared_ptr<Vertex> >& vx_pri)
//==========================================================================
{
  // Compute closest boundary point to current vertex
  Point vx_pos = vx->getVertexPoint();
  shared_ptr<Loop> loop = face->getBoundaryLoop(0);

  double par, dist;
  int ind;
  Point close;
  loop->closestPoint(vx_pos, ind, par, close, dist);

  shared_ptr<ftEdgeBase> edg = loop->getEdge(ind);
  double t1 = edg->tMin();
  double t2 = edg->tMax();
  double tdel = (t2 - t1)/10.0;
  ftEdge *edg2 = edg->geomEdge();
  if (edg2)
    {
      // Check if the edge already is connected to a prioritized vertex
      size_t ki;
      for (ki=0; ki<vx_pri.size(); ++ki)
	if (edg2->hasVertex(vx_pri[ki].get()))
	  break;
      if (ki < vx_pri.size())
	edg2 = NULL;
    }
	    
  if (edg2 && par > t1+tdel && par < t2-tdel)
    {
      // Create split vertex
      shared_ptr<ftEdge> newedge = edg2->split2(par);
      shared_ptr<Vertex> tmp_vx = edg2->getCommonVertex(newedge.get());
      vx_pri.push_back(tmp_vx);
    }

}
