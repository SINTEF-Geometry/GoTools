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
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/compositemodel/RegularizeUtils.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/compositemodel/Path.h"
//#include "GoTools/topology/tpTopologyTable.h"
#include "GoTools/geometry/SurfaceTools.h"
#include <fstream>
#include <cstdlib>

//#define DEBUG_REG;

using std::vector;
using std::set;
using std::make_pair;

namespace Go {

//==========================================================================
  RegularizeFaceSet::RegularizeFaceSet(vector<shared_ptr<ftSurface> > faces, 
				       double epsge, double angtol,
				       bool split_in_cand)
    : split_mode_(1), split_in_cand_(split_in_cand)
//==========================================================================
{
  model_ = shared_ptr<SurfaceModel>(new SurfaceModel(epsge, epsge, 10.0*epsge,
						     angtol, 10.0*angtol,
						     faces));
  cand_split_.resize(faces.size());
}

//==========================================================================
  RegularizeFaceSet::RegularizeFaceSet(vector<shared_ptr<ftSurface> > faces, 
				       double gap, double neighbour, 
				       double kink, double bend, 
				       bool split_in_cand)
    : split_mode_(1), split_in_cand_(split_in_cand)
//==========================================================================
{
  model_ = shared_ptr<SurfaceModel>(new SurfaceModel(gap, gap, neighbour,
						     kink, bend, faces));
  cand_split_.resize(faces.size());
}

//==========================================================================
    RegularizeFaceSet::RegularizeFaceSet(shared_ptr<SurfaceModel> model,
					 bool split_in_cand)
      : split_mode_(1), split_in_cand_(split_in_cand)
//==========================================================================
{
  model_ = model;
  cand_split_.resize(model->nmbEntities());
}

//==========================================================================
RegularizeFaceSet::~RegularizeFaceSet()
//==========================================================================
{

}

//==========================================================================
void RegularizeFaceSet::setFaceCorrespondance(int idx1, int idx2)
//==========================================================================
{
  corr_faces_.push_back(std::make_pair(idx1, idx2));
}

//==========================================================================
vector<shared_ptr<ftSurface> > 
  RegularizeFaceSet::getRegularFaces(bool reverse_sequence)
//==========================================================================
{
  divide(reverse_sequence);
  vector<shared_ptr<ftSurface> > faces;
  int nmb_faces = model_->nmbEntities();
  faces.reserve(nmb_faces);
  for (int ki=0; ki<nmb_faces; ++ki)
    faces.push_back(model_->getFace(ki));
  return faces;
}

//==========================================================================
shared_ptr<SurfaceModel> RegularizeFaceSet::getRegularModel(bool reverse_sequence)
//==========================================================================
{
  divide(reverse_sequence);
  return model_;
}

//==========================================================================
  void RegularizeFaceSet::divide(bool reverse_sequence)
//==========================================================================
{
  // Divide the faces one by one. First collect all faces
  vector<shared_ptr<ftSurface> > faces = model_->allFaces();
  int nmb_faces = (int)faces.size();

#ifdef DEBUG_REG
  // Check all vertices
  vector<shared_ptr<Vertex> > curr_vx;
  model_->getAllVertices(curr_vx);
  for (size_t kf=0; kf<curr_vx.size(); ++kf)
    if (!curr_vx[kf]->checkVertexTopology())
      {
	std::ofstream vx_of("error_vx.g2");
	vx_of << "400 1 0 4 255 0 0 255 " << std::endl;
	vx_of << "1" << std::endl;
	vx_of << curr_vx[kf]->getVertexPoint() << std::endl;
	std::cout << " Error in vertex topology " << std::endl;
      }
#endif

  vector<shared_ptr<ftSurface> > remaining_faces;
  remaining_faces.insert(remaining_faces.end(), faces.begin(), faces.end());

  // Set face correspondance
  computeFaceCorrespondance(faces);

  // Prioritize division of faces to avoid premature decisions
  vector<int> perm(nmb_faces);
  int kj;
  for (kj=0; kj<nmb_faces; ++kj)
    perm[kj] = kj;

  // Perform sorting
  prioritizeFaces(faces, perm);

  if (reverse_sequence)
    {
      for (int ka=0; ka<nmb_faces/2; ++ka)
	std::swap(perm[ka], perm[nmb_faces-1-ka]);
    }

#ifdef DEBUG_REG
      std::ofstream of0("all_post_regface.g2");
#endif
  for (int kj=0; kj<nmb_faces; ++kj)
    {
      vector<shared_ptr<Vertex> > pre_vx1;
      model_->getAllVertices(pre_vx1);

      shared_ptr<ftSurface> curr = faces[perm[kj]];
#ifdef DEBUG_REG
      std::ofstream of1("pre_regface.g2");
      curr->surface()->writeStandardHeader(of1);
      curr->surface()->write(of1);
#endif
#ifdef DEBUG_REG2
      std::ofstream t1("missing_twin1.g2");
      vector<shared_ptr<ftEdge> > edges = curr->getAllEdges();
      for (size_t ix=0; ix<edges.size(); ++ix)
	if (!edges[ix]->twin())
	  {
	    std::cout << "Missing twin pointer. RegularizeFaceSet. initial divide" << std::endl;
	    t1 << "410 1 0 4 155 55 0 255" << std::endl;
	    t1 << "1" << std::endl;
	    t1 << edges[ix]->point(edges[ix]->tMin()) << " " << edges[ix]->point(edges[ix]->tMax()) << std::endl;
	  }
#endif

      // Check configuration
      if (curr->twin() && curr->twin()->getBody() == curr->getBody())
	continue;  // This operation is risky in the current configuration.
      // Skip it and hope the situation is resolved at a later stage

      // Check if the topology update implies that some radial edge
      // information may get lost. In that case, store a pointer
      // to the relevant EdgeVertex instances
      // Fetch this information before the model is updated
      vector<shared_ptr<EdgeVertex> > edgevx;
      vector<std::pair<Point,Point> > endpts;
      getSeamRadialEdge(curr.get(), edgevx, endpts);

      ftSurface *twin = curr->twin();

      bool split_in_cand = split_in_cand_;
      if (cand_split_[perm[kj]].size() >  0 && split_mode_ > 1)
	{
	  size_t kr;
	  for (kr=0; kr<corr_faces_.size(); ++kr)
	    if (corr_faces_[kr].first == perm[kj] ||
		corr_faces_[kr].second == perm[kj])
	      break;
	  if (kr < corr_faces_.size())
	    {
	      size_t idx = (corr_faces_[kr].first == perm[kj]) ?
		corr_faces_[kr].second : corr_faces_[kr].first;
	      shared_ptr<ftSurface> other = faces[idx];
	      if (curr->nmbBoundaryLoops() < other->nmbBoundaryLoops())
		split_in_cand = true;
	    }
	}

      RegularizeFace regularize(curr, model_, split_in_cand);
      regularize.setSplitMode(split_mode_);
      // if (kj == nmb_faces-1)
      // 	regularize.setDivideInT(false);  // The last T-joint division is better done here
      if (cand_split_[perm[kj]].size() >  0)
	regularize.setCandSplit(cand_split_[perm[kj]]);
      regularize.classifyVertices();

      // Provide regularize with information about the faces not treated yet.
      // First remove current face
      for (size_t kh=0; kh<remaining_faces.size(); ++kh)
	{
	  if (remaining_faces[kh].get() == curr.get())
	    {
	      remaining_faces.erase(remaining_faces.begin()+kh);
	      break;
	    }
	}
      regularize.setNonTjointFaces(remaining_faces);

      vector<shared_ptr<ftSurface> > faces2 = 
	regularize.getRegularFaces();

      // Store info about vertex point correspondance
      vector<pair<Point,Point> > corr_vx_pts = regularize.fetchVxPntCorr();
      if (corr_vx_pts.size() > 0)
	corr_vx_pts_.insert(corr_vx_pts_.end(), corr_vx_pts.begin(), 
			    corr_vx_pts.end());

#ifdef DEBUG_REG
      std::ofstream of2("post_regface.g2");
      for (size_t kr=0; kr<faces2.size(); ++kr)
	{
	  faces2[kr]->surface()->writeStandardHeader(of2);
	  faces2[kr]->surface()->write(of2);
	}
#endif
#ifdef DEBUG_REG
      for (size_t kr=0; kr<faces2.size(); ++kr)
	{
	  faces2[kr]->surface()->writeStandardHeader(of0);
	  faces2[kr]->surface()->write(of0);
	}
#endif
      // Update topology information
      if (faces2.size() > 1)
	{
	  // Set candidate pairs of split parameters
	  size_t kr;
	  for (kr=0; kr<corr_faces_.size(); ++kr)
	    if (corr_faces_[kr].first == perm[kj] ||
		corr_faces_[kr].second == perm[kj])
	      break;
	  if (kr < corr_faces_.size())
	    {
	      size_t idx = (corr_faces_[kr].first == perm[kj]) ?
		corr_faces_[kr].second : corr_faces_[kr].first;
	      vector<pair<pair<Point,int>, std::pair<Point,int> > > end_split =
		getEndSplit(curr, faces2);
	      cand_split_[idx] = end_split;
	    }

	  // Fetch info about removed seams
	  vector<Point> seam_joints = regularize.getSeamJointInfo();
	  if (seam_joints.size() > 0)
	    seam_joints_.insert(seam_joints_.end(), seam_joints.begin(),
				seam_joints.end());

	  for (kr=0; kr<faces2.size(); ++kr)
	    {
	      attachRadialEdge(faces2[kr].get(), edgevx, endpts, 
			       model_->getTolerances().neighbour);
	    }
#ifdef DEBUG_REG2
      std::ofstream t2("missing_twin2.g2");
	  for (kr=0; kr<faces2.size(); ++kr)
	    {
	      vector<shared_ptr<ftEdge> > edges = 
		faces2[kr]->getAllEdges();
	      for (size_t ix=0; ix<edges.size(); ++ix)
		if (!edges[ix]->twin())
		  {
		    std::cout << "Missing twin pointer. RegularizeFaceSet. divide" << std::endl;
		    t2 << "410 1 0 4 155 55 0 255" << std::endl;
		    t2 << "1" << std::endl;
		    t2 << edges[ix]->point(edges[ix]->tMin()) << " " << edges[ix]->point(edges[ix]->tMax()) << std::endl;
		  }
	    }
#endif

	  if (twin)
	    {
	      Body *bd = twin->getBody();
	      if (bd)
		{
		  shared_ptr<SurfaceModel> shell = bd->getShell(twin);
		  if (shell.get())
		    {
		      shell->regularizeTwin(twin, faces2);
		      size_t kr;
		      for (kr=0; kr<modified_models_.size(); ++kr)
			if (modified_models_[kr] == shell.get())
			  break;
		      if (kr == modified_models_.size())
			modified_models_.push_back(shell.get());
		    }
		      
		}
	    }
	}
      vector<shared_ptr<Vertex> > post_vx1;
      model_->getAllVertices(post_vx1);
    }

#ifdef DEBUG_REG
    std::ofstream debug("split_faces.g2");
    for (int kv=0; kv<model_->nmbEntities(); ++kv)
      {
	model_->getSurface(kv)->writeStandardHeader(debug);
	model_->getSurface(kv)->write(debug);
      }
#endif
  
  // Split in T-joint like situations that ruins the corner-to-corner
  // configuaration with 4-sided surfaces that we want
    vector<shared_ptr<Vertex> > pre_vx;
    model_->getAllVertices(pre_vx);

    // Add still not removed seams to the seam repository
    seamClassification();
    
    // Perform split
    splitInTJoints();

#ifdef DEBUG_REG
    std::ofstream of("post_T_split.g2");
    for (int kr=0; kr<model_->nmbEntities(); ++kr)
      {
	shared_ptr<ParamSurface> sf = model_->getSurface(kr);
	sf->writeStandardHeader(of);
	sf->write(of);
      }
#endif

    vector<shared_ptr<Vertex> > post_vx;
    model_->getAllVertices(post_vx);

}

//==========================================================================
void RegularizeFaceSet::splitInTJoints()
//==========================================================================
{
  int nmb_faces = model_->nmbEntities();
  //double angtol = model_->getTolerances().kink;
  double angtol = model_->getTolerances().bend;

  // Postpone the treatment of faces without corners
  int last_face = nmb_faces - 1;
  for (int ki=0; ki<nmb_faces; ++ki)
    {
      if (last_face < 0)
	break;  // Only faces without corners

      vector<shared_ptr<Vertex> > corner = 
	model_->getFace(ki)->getCornerVertices(angtol, 0);
      if (corner.size() == 0)
	{
	  model_->swapFaces(ki, last_face);
	  ki--;
	  last_face--;
	}
    }
  
#ifdef DEBUG_REG
  std::ofstream of("T_face.g2");
  for (int kr=0; kr<nmb_faces; ++kr)
    {
      shared_ptr<ParamSurface> sf = model_->getSurface(kr);
      sf->writeStandardHeader(of);
      sf->write(of);
    }

  bool shell_ok = model_->checkShellTopology();
  if (!shell_ok)
    std::cout << "Corrupt shell topology" << std::endl;

  vector<shared_ptr<Vertex> > vxs;
  model_->getAllVertices(vxs);
  std::ofstream of2("Error_vx.g2");
  for (size_t kr=0; kr<vxs.size(); ++kr)
    if (!vxs[kr]->checkVertexTopology())
      {
	of2 << "400 1 0 4 200 0 55 255" << std::endl;
	of2 << "1 " << std::endl;
	of2 << vxs[kr]->getVertexPoint() << std::endl;
      }
#endif

  bool changed = true;
  while (changed)
    {
      nmb_faces = model_->nmbEntities();

#ifdef DEBUG_REG
      std::ofstream debug("T_faces.g2");
      for (int kv=0; kv<nmb_faces; ++kv)
	{
	  model_->getSurface(kv)->writeStandardHeader(debug);
	  model_->getSurface(kv)->write(debug);
	}
#endif
      changed = false;
      for (int ki=0; ki<nmb_faces; ++ki)
	{
	  shared_ptr<ftSurface> curr = model_->getFace(ki);
	  Body *curr_body = curr->getBody();

	  // Find potensial T-joints
	  if (curr->nmbBoundaryLoops() == 0)
	    continue;  // Can't split

#ifdef DEBUG_REG
	  std::ofstream debug("curr_T_face.g2");
	  curr->surface()->writeStandardHeader(debug);
	  curr->surface()->write(debug);
#endif

	  // Get corner vertices
	  vector<shared_ptr<Vertex> > corner = 
	    curr->getCornerVertices(angtol, 0);

	  int kj;
	  if (corner.size() == 0)
	    {
	      // Postpone the treatment of faces without corners if a face with 
	      // corners can be treated first
	      for (kj=ki+1; kj<nmb_faces; ++kj)
		{
		  shared_ptr<ftSurface> curr2 = model_->getFace(kj);
		  if (curr2->nmbBoundaryLoops() == 0)
		    continue;  // Can't split
		  
		  vector<shared_ptr<Vertex> > corner2 = 
		    curr2->getCornerVertices(angtol, 0);
		  if (corner2.size() > 0)
		    {
		      model_->swapFaces(ki, kj);
		      curr = curr2;
		      curr_body = curr->getBody();
		      corner = corner2;
#ifdef DEBUG_REG
		      curr->surface()->writeStandardHeader(debug);
		      curr->surface()->write(debug);
#endif
		      break;
		    }
		}
	    }
	      

	  // Get potensial T-joints
	  vector<shared_ptr<Vertex> > Tvx = 
	    curr->getNonCornerVertices(angtol, 0);
	  removeInsignificantVertices(Tvx);

	  // Check if the vertex really indicates a T-joint
	  for (kj=0; kj<(int)Tvx.size(); )
	    {
	      vector<ftSurface*> vx_faces = Tvx[kj]->faces();
	      for (size_t kr=0; kr<vx_faces.size();)
		{
		  if (vx_faces[kr]->getBody() != curr_body)
		    vx_faces.erase(vx_faces.begin()+kr);
		  else
		    kr++;
		}
	      if (vx_faces.size() > 2)
		kj++;
	      else if (vx_faces.size() == 1)
		{
		  // Only this face. Not a T-joint
		  Tvx.erase(Tvx.begin()+kj);
		}
	      else
		{
		  // One adjacent face. Check if the vertex lies in a 
		  // corner of this face
		  ftSurface *other = (vx_faces[0] == curr.get()) ?
		    vx_faces[1] : vx_faces[0];
		  vector<shared_ptr<Vertex> > other_corner = 
		    other->getCornerVertices(angtol);
		  vector<shared_ptr<Vertex> >::iterator vxp =
		    std::find(other_corner.begin(), other_corner.end(), 
			      Tvx[kj]);
		  if (vxp == other_corner.end())
		    Tvx.erase(Tvx.begin() + kj); // Not a corner
		  else
		    {
		      // Check if the vertex lies at a seam of the
		      // adjacent surface. In that case, skip this vertex
		      vector<ftEdge*> edg = Tvx[kj]->getFaceEdges(other);
		      bool found_seam = false;
		      for (size_t k1=0; k1<edg.size(); ++k1)
			for (size_t k2=k1+1; k2<edg.size(); ++k2)
			  if (edg[k1]->twin() == edg[k2])
			    found_seam = true;
		      if (found_seam)
			Tvx.erase(Tvx.begin() + kj);  // Seam
		      else
			kj++;
		    }
		}
	    }

	  // Check configuration
	  if (curr->twin() && curr->twin()->getBody() == curr->getBody())
	    continue;  // This operation is risky in the current configuration.
	  // Skip it and hope the situation is resolved at a later stage

	  if (corner.size() + Tvx.size() > 4)
	    {
	      // Check current number of faces
	      int nmb_faces = model_->nmbEntities();

 	      // Number of corners and T-joints exceeds 4. Not
	      // appropriate for a 4-sided surface divide
	      vector<shared_ptr<ftSurface> > faces = 
		divideInTjoint(curr, Tvx, corner, changed);
#ifdef DEBUG_REG
		  vector<shared_ptr<Vertex> > vxs;
		  model_->getAllVertices(vxs);
		  std::ofstream of3("Error_vx4.g2");
		  for (size_t kr=0; kr<vxs.size(); ++kr)
		    if (!vxs[kr]->checkVertexTopology())
		      {
			of3 << "400 1 0 4 200 0 55 255" << std::endl;
			of3 << "1 " << std::endl;
			of3 << vxs[kr]->getVertexPoint() << std::endl;
		      }
#endif
	      if (faces.size() > 0)
		{
		  // Check if the topology update implies that some radial edge
		  // information may get lost. In that case, store a pointer
		  // to the relevant EdgeVertex instances
		  vector<shared_ptr<EdgeVertex> > edgevx;
		  vector<std::pair<Point,Point> > endpts;
		  getSeamRadialEdge(curr.get(), edgevx, endpts);

		  ftSurface *twin = curr->twin();
		  model_->removeFace(curr);
		  for (size_t kr=0; kr<faces.size(); ++kr)
		    {
		      model_->append(faces[kr]);
		      attachRadialEdge(faces[kr].get(), edgevx, endpts, 
				       model_->getTolerances().neighbour);
		    }
		  for (size_t kr=0; kr<faces.size(); ++kr)
		    {
		      vector<shared_ptr<ftEdge> > edges = 
			faces[kr]->getAllEdges();
#ifdef DEBUG_REG2
		      for (size_t ix=0; ix<edges.size(); ++ix)
			if (!edges[ix]->twin())
			  std::cout << "Missing twin pointer. RegularizeFaceSet. splitInTjoints" << std::endl;
#endif
		    }

		  if (twin)
		    {
		      Body *bd = twin->getBody();
		      if (bd)
			{
			  shared_ptr<SurfaceModel> shell = 
			    bd->getShell(twin);
			  if (shell.get())
			    shell->regularizeTwin(twin, faces);
			}
		    }

		  for (kj=0; kj<(int)faces.size(); ++kj)
		    {
		      // Identify surfaces that can be joined
		      // across a seam
		      int pardir, status, status2;
		      ftSurface *other =
			identifySeamFaces(faces[kj], pardir, status);
		      if (other)
			{
#ifdef DEBUG_REG
			  std::ofstream seam("seam_surfs.g2");
			  faces[kj]->surface()->writeStandardHeader(seam);
			  faces[kj]->surface()->write(seam);
			  other->surface()->writeStandardHeader(seam);
			  other->surface()->write(seam);
#endif

			  // Check twins
			  ftSurface *twin1 = faces[kj]->twin();
			  ftSurface *twin2 = other->twin();
			  shared_ptr<ftSurface> merge1;
			  shared_ptr<ftSurface> merge2;
			  bool update = true;
			  if ((twin1 && !twin2) || (twin2 && !twin1))
			    update = false;
			  if (twin1 && twin2)
			    {
#ifdef DEBUG_REG
			      twin1->surface()->writeStandardHeader(seam);
			      twin1->surface()->write(seam);
			      twin2->surface()->writeStandardHeader(seam);
			      twin2->surface()->write(seam);
#endif

			      shared_ptr<ftEdge> edge1, edge2;
			      // bool neighbour = 
			      // twin1->areNeighbours(twin2,
			      // edge1, edge2, 0);
			      int dir2;
			      status2 = 
				getParameterDirection(twin1, twin2,
						      edge1, edge2,
						      model_->getTolerances().neighbour,
						      dir2);
			      if (!status2)
				update = false;
			      else if (twin1->getBody() != twin2->getBody())
				update = false;
			      
			      if (update)
				{
				  shared_ptr<SurfaceModel> shell =
				    twin1->getBody()->getShell(twin1);
				  vector<Point> seam_joints2;
				  if (shell.get())
				    {
				      if (status2 == 1)
					merge2 = shell->mergeSeamFaces(twin1, twin2,
								       dir2, 
								       seam_joints2);
				      else
					merge2 = shell->mergeSeamCrvFaces(twin1, twin2,
									  seam_joints2);
				    }
					
				  if (!merge2.get())
				    update = false;
#ifdef DEBUG_REG
				  if (merge2.get())
				    {
				      merge2->surface()->writeStandardHeader(seam);
				      merge2->surface()->write(seam);
				    }
#endif
				}
			    }
			  if (update)
			    {
			      vector<Point> seam_joints;
			      if (status == 1)
				merge1 = model_->mergeSeamFaces(faces[kj].get(), 
								other, 
								pardir,
								seam_joints);
			      else
				merge1 = model_->mergeSeamCrvFaces(faces[kj].get(), 
								   other, 
								   seam_joints);

			      if (merge1.get())
				{
#ifdef DEBUG_REG
				  merge1->surface()->writeStandardHeader(seam);
				  merge1->surface()->write(seam);
#endif

				  // Add new removed seams to the seam repository
				  seam_joints_.insert(seam_joints_.end(),
						      seam_joints.begin(),
						      seam_joints.end());
				}
			      if (merge1.get() && merge2.get())
				merge1->connectTwin(merge2.get(), 
						    model_->getTolerances().neighbour);
			    }
			}
		    }

		  changed = true;
#ifdef DEBUG_REG
		  vxs.clear();
		  model_->getAllVertices(vxs);
		  for (size_t kr=0; kr<vxs.size(); ++kr)
		    if (!vxs[kr]->checkVertexTopology())
		      {
			of3 << "400 1 0 4 0 200 55 255" << std::endl;
			of3 << "1 " << std::endl;
			of3 << vxs[kr]->getVertexPoint() << std::endl;
		      }
#endif

		}
	      if (model_->nmbEntities() != nmb_faces)
		changed = true;  // A modification has been applied
	    }
	  if (changed)
	    break;
	}
    }
}

//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeFaceSet::divideInTjoint(shared_ptr<ftSurface>& face,
				  vector<shared_ptr<Vertex> >& Tvx,
				  vector<shared_ptr<Vertex> >& corner,
				  bool& changed)
//==========================================================================
{
  vector<shared_ptr<ftSurface> > faces;
  if (Tvx.size() + corner.size() <= 4)
    return faces;  // No point in dividing

  // Select the T-joint to divide in
  shared_ptr<Vertex> curr;
  double min_dist = 0.0;
  size_t curr_idx;
  size_t ki, kr;
  while (Tvx.size() > 0)
    {
      min_dist = 0.0;
      for (ki=0; ki<Tvx.size(); ++ki)
	{
	  // Compute distance between T vertex and corners on both sides
	  vector<ftEdge*> edges = Tvx[ki]->getFaceEdges(face.get());
	  if (edges.size() != 2)
	    continue;  // Does not make sense

	  double len = 0.0;
	  for (size_t kj=0; kj<edges.size(); ++kj)
	    {
	      ftEdge* curr_edge = edges[kj];
	      shared_ptr<Vertex> curr_vx = Tvx[ki];
	      double tp = edges[kj]->parAtVertex(Tvx[ki].get());
	      bool next = (tp - edges[kj]->tMin() < edges[kj]->tMax() - tp);
	      shared_ptr<Vertex> other;
	      if (corner.size() > 0)
		{
		  while (true)
		    {
		      other = curr_edge->getOtherVertex(curr_vx.get());
		      if (!other.get())
			{
			  other = Tvx[ki];
			  break;
			}
		      vector<shared_ptr<Vertex> >::iterator vxp =
			std::find(corner.begin(), corner.end(), other);
		      if (vxp != corner.end())
			break;

		      curr_vx = other;
		      curr_edge = (next) ? curr_edge->next()->geomEdge() :
			curr_edge->prev()->geomEdge();
		      
		      // To avoid an infinite loop
		      if (curr_vx.get() == Tvx[ki].get())
			break;
		    }
		}
	      else
		other = curr_edge->getOtherVertex(Tvx[ki].get());
	      len += Tvx[ki]->getVertexPoint().dist(other->getVertexPoint());
	    }

	  if (len > min_dist)
	    {
	      min_dist = len;
	      curr = Tvx[ki];
	      curr_idx = ki;
	    }
	}

      if (!curr.get())
	break;  // Cannot divide

      // Some testing
      vector<ftSurface*> vx_faces = curr->faces();
      if (vx_faces.size() == 3)
	{
	  int idx1, idx2;
	  if (vx_faces[0] == face.get())
	    {
	      idx1 = 1;
	      idx2 = 2;
	    }
	  else if (vx_faces[1] == face.get())
	    {
	      idx1 = 0;
	      idx2 = 2;
	    }
	  else
	    {
	      idx1 = 0;
	      idx2 = 1;
	    }

	  shared_ptr<BoundedSurface> srf1 =
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(vx_faces[idx1]->surface());
	  shared_ptr<BoundedSurface> srf2 =
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(vx_faces[idx2]->surface());

	  if (srf1.get() && srf2.get())
	    {
	      if (srf1->underlyingSurface().get() == 
		  srf2->underlyingSurface().get())
		{
		  shared_ptr<ftEdge> edg1, edg2;
		  bool neighbour = vx_faces[idx1]->areNeighbours(vx_faces[idx2],
								 edg1, edg2);
		  if (neighbour)
		    {
		      shared_ptr<CurveOnSurface> bd_cv =
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>(edg1->geomCurve());
		      bool same;
		      // int bd = bd_cv->whichBoundary(1.0e-6, same);
		      bd_cv->whichBoundary(1.0e-6, same);

#ifdef DEBUG_REG
		      std::ofstream of0("adj_faces.g2");
		      srf1->writeStandardHeader(of0);
		      srf1->write(of0);
		      srf2->writeStandardHeader(of0);
		      srf2->write(of0);
		      int stop_brek;
		      stop_brek = 1;
#endif
		    }
		}
	    }
	}

#ifdef DEBUG_REG
      std::ofstream of("curr_T_vx.g2");
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << "1" << std::endl;
      of << curr->getVertexPoint() << std::endl;
#endif

      // Check if merge is preferred before divide in this situation
      ftSurface *merge1 = NULL;
      ftSurface *merge2 = NULL;
      shared_ptr<Vertex> other_vx;
      int pardir1, pardir2;
      double parval1, parval2;
      bool atstart1, atstart2;
      pair<Point, Point> co_par1;
      pair<Point, Point> co_par2;
      int found_merge = mergeSituation(face.get(), curr, merge1, pardir1, 
				       parval1, atstart1, merge2, pardir2,
				       parval2, atstart2, other_vx, co_par1, co_par2);

#ifdef DEBUG_REG
    if (!curr->checkVertexTopology())
      {
	std::ofstream of2("Error_vx2.g2");
	of2 << "400 1 0 4 200 0 55 255" << std::endl;
	of2 << "1 " << std::endl;
	of2 << curr->getVertexPoint() << std::endl;
      }
#endif

    if (found_merge > 0 && merge1->nmbAdjacencies(merge2) > 1)
      {
	// Check if the edges are adjacent
	vector<shared_ptr<ftEdge> > edges;
	int idx=0;
	shared_ptr<ftEdge> edg1, edg2;
	while (merge1->areNeighbours(merge2, edg1, edg2, idx))
	  {
	    edges.push_back(edg1);
	    idx++;
	  }
	for (size_t kv=1; kv<edges.size(); ++kv)
	  {
	    shared_ptr<Vertex> commonv = edges[kv-1]->getCommonVertex(edges[kv].get());
	    if (!commonv.get())
	    {
	      found_merge = 3;
	      return faces;
	    }
	  }
      }
      
      if (found_merge == 2)
	{
#ifdef DEBUG_REG
	  std::ofstream of1("merge_situation.g2");
	  merge1->surface()->writeStandardHeader(of1);
	  merge1->surface()->write(of1);
	  merge2->surface()->writeStandardHeader(of1);
	  merge2->surface()->write(of1);
#endif
	  if (merge1->twin() || merge2->twin())
	    MESSAGE("RegularizeFaceSet::divideInTjoint. Twin info not maintained in merge");

	  vector<Point> seam_joints;
	  shared_ptr<ftSurface> merged_face = model_->mergeFaces(merge1, pardir1, parval1, 
								 atstart1, merge2, pardir2, 
								 parval2, atstart2, co_par1, 
								 co_par2, seam_joints);

	  if (merged_face.get())
	    break;   // T-joint resolved
	  else
	    found_merge = 1;  // Try another approach
	    // Tvx.erase(Tvx.begin() + curr_idx);  // Look for another T-joint

	  // // Return empty face vector to avoid further topology updates
	  // vector<shared_ptr<ftSurface> > faces;
	  // return faces;
	}

      if (found_merge == 1)
	{
	  // The faces meet with requested continuity, but the common boundary
	  // is not a constant parameter curve in both faces. If possible,
	  // approximate the face with a non-trimmed spline surface
	  int nmb_bd1 = merge1->nmbOuterBdCrvs(model_->getTolerances().gap,
					       model_->getTolerances().neighbour,
					       model_->getTolerances().bend); 
	  int nmb_bd2 = merge2->nmbOuterBdCrvs(model_->getTolerances().gap,
					       model_->getTolerances().neighbour,
					       model_->getTolerances().bend); 
	  if (nmb_bd1 == 4 && nmb_bd2 == 4)
	    {
	      // Both faces have 4 corners
	      shared_ptr<ftSurface> m1 = model_->replaceRegularSurface(merge1, true);
	      shared_ptr<ftSurface> m2 = model_->replaceRegularSurface(merge2, true);
	      if (m1.get() || m2.get())
	      	{
		  // The T-vertex has changed. Thus, we cannot continue the merge
		  // at this stage. Return and continue the search for T-joints
		  changed = true;
		  // Return empty face vector to avoid further topology updates
		  vector<shared_ptr<ftSurface> > faces;
#ifdef DEBUG_REG
		  vector<shared_ptr<Vertex> > vxs;
		  model_->getAllVertices(vxs);
		  std::ofstream of3("Error_vx3.g2");
		  for (size_t kr=0; kr<vxs.size(); ++kr)
		    if (!vxs[kr]->checkVertexTopology())
		      {
			of3 << "400 1 0 4 200 0 55 255" << std::endl;
			of3 << "1 " << std::endl;
			of3 << vxs[kr]->getVertexPoint() << std::endl;
		      }
#endif
		  return faces;
	      // 	  // Check if merging is possible with the new situation
	      // 	  found_merge = mergeSituation(face.get(), curr, merge1, pardir1, 
	      // 				       parval1, atstart1, merge2, pardir2,
	      // 				       parval2, atstart2, other_vx, 
	      // 				       co_par1, co_par2);
	      	}
	      //Tvx.erase(Tvx.begin() + curr_idx);  // Look for another T-joint
	    }
	  //Tvx.erase(Tvx.begin() + curr_idx);  // Look for another T-joint
	}

      if (!found_merge)
	{
	  // A T-joint vertex is found
	  // Get candidate vertices
	  vector<shared_ptr<Vertex> > vx = face->vertices();
	  vector<shared_ptr<Vertex> > cand_vx;
	  ftEdge* cand_edge = 0;
	  selectCandidateSplit(face, curr, vx, cand_vx, cand_edge);
 
	  // Divide
	  //cand_edge = 0;
	  vector<shared_ptr<Vertex> > non_corner;
	  Point dummy;
	  vector<shared_ptr<Vertex> > Tvx2 = Tvx;
	  Tvx2.erase(Tvx2.begin() + curr_idx);
	  for (kr=0; kr<Tvx2.size(); ++kr)
	    {
	      int kh;
	      for (kh=(int)cand_vx.size()-1; kh>=0; --kh)
		if (Tvx2[kr].get() == cand_vx[kh].get())
		  {
		    // Vertex exists in both samples
		    cand_vx.erase(cand_vx.begin()+kh);
		  }
	    }

	  faces = RegularizeUtils::divideVertex(face, curr, cand_vx, 
						cand_edge, Tvx2,
						model_->getTolerances().gap,
						model_->getTolerances().neighbour,
						model_->getTolerances().kink,
						model_->getTolerances().bend,
						non_corner, dummy, dummy);
	  if (faces.size() > 1)
	    break;  // T-joints resolved
	  else
	    Tvx.erase(Tvx.begin() + curr_idx);  // Look for another T-joint
	}
    }
  return faces;
	    
}

//==========================================================================
void 
RegularizeFaceSet::selectCandidateSplit(shared_ptr<ftSurface> face,
					shared_ptr<Vertex> select_vx,
					vector<shared_ptr<Vertex> >& vx,
					vector<shared_ptr<Vertex> >& cand_vx,
					ftEdge*& cand_edge)
//==========================================================================
{
  cand_edge = 0;
  vector<shared_ptr<Vertex> > other_vx;
  vector<ftEdge*> edges = select_vx->getFaceEdges(face.get());
  for (size_t ki=0; ki<edges.size(); ++ki)
    other_vx.push_back(edges[ki]->getOtherVertex(select_vx.get()));

  for (size_t ki=0; ki<vx.size(); ++ki)
    {
      if (vx[ki].get() == select_vx.get())
	continue;   // Not a candiate
      if (vx[ki]->sameEdge(select_vx.get()))
	continue; 
      size_t kj;
      for (kj=0; kj<other_vx.size(); ++kj)
	if (vx[ki]->sameEdge(other_vx[kj].get()))
	  break;

      if (kj >= other_vx.size())
	cand_vx.push_back(vx[ki]);
    }
  removeInsignificantVertices(cand_vx);

  if (cand_vx.size() == 0)
    {
      // Find candiate edge
      if (edges.size() > 0)
	{
	  if (edges[0]->parAtVertex(select_vx.get()) == edges[0]->tMin())
	    cand_edge = edges[0]->next()->next()->geomEdge();
	  else
	    cand_edge = edges[0]->prev()->prev()->geomEdge();
	}	    
    }

}					   

//==========================================================================
double RegularizeFaceSet::getSegmentAngle(shared_ptr<ftSurface> face,
					  shared_ptr<Vertex> vx1,
					  shared_ptr<Vertex> vx2,
					  Point& pnt, Point& normal)
//==========================================================================
{
  // vx1 is the source vertex
  Point vec = vx2->getVertexPoint() - vx1->getVertexPoint();
  double ang = vec.angle(normal);
  ang = 0.5*M_PI - ang;
  ang = fabs(ang);

  // Check direction
  vector<ftEdge*> edges = vx1->getFaceEdges(face.get());
  if (edges.size() < 1)
      ang += M_PI;
  else
    {
      double t1 = edges[0]->parAtVertex(vx1.get());
      double t2 = edges[1]->parAtVertex(vx1.get());
      Point tan1 = edges[0]->tangent(t1);
      Point tan2 = edges[1]->tangent(t2);
      if (edges[0]->tMax() - t1 < t1 - edges[0]->tMin())
	tan1 *= -1;
      else 
	tan2 *= -1;
      tan1.normalize();
      tan2.normalize();
      Point vec2 = 0.5*(tan1 + tan2);
      Point par = vx1->getFacePar(face.get());
      Point norm1 = face->normal(par[0], par[1]);
      Point norm2 = tan2.cross(tan1);
      if (norm1*norm2 < 0.0)
	vec2 *= -1;
      if (vec*vec2 < 0.0)
	ang += M_PI;
    }
      
  return ang;
}

//==========================================================================
  vector<pair<pair<Point,int>, pair<Point,int> > > 
  RegularizeFaceSet::getEndSplit(shared_ptr<ftSurface> prev_face,
				 vector<shared_ptr<ftSurface> >& faces)
//==========================================================================
{
  vector<pair<pair<Point,int>, pair<Point,int> > > params;
  double tol = model_->getTolerances().neighbour;

  // Fetch all edges once
  vector<shared_ptr<ftEdge> > edges;
  size_t ki, kj, kr;
  for (ki=0; ki<faces.size(); ++ki)
    {
      vector<shared_ptr<ftEdge> > curr = faces[ki]->getAllEdges();

      for (kj=0; kj<curr.size(); kj++)
	{
	  // Check if the edge or a twin edge exists already
	  for (kr=0; kr<edges.size(); ++kr)
	    if (edges[kr].get() == curr[kj].get() ||
		edges[kr]->twin() == curr[kj].get())
	      break;
	  if (kr == edges.size())
	    edges.push_back(curr[kj]);
	}
    }

  // Find end parameters corresponding to edges
  shared_ptr<Loop> outerbd = prev_face->getBoundaryLoop(0);
  for (ki=0; ki<edges.size(); ++ki)
    {
      // Fetch edge points in chain
      shared_ptr<Vertex> v1, v2;
      vector<ftEdge*> chain = Path::edgeChain(edges[ki].get(), 
					      /*model_->getTolerances().bend,*/
					      model_->getTolerances().kink,
					      v1, v2);
      Point pos1 = v1->getVertexPoint();
      Point pos2 = v2->getVertexPoint();

      // Check if the vertex lies at a face boundary
      Point dummy_vec;
      double upar1, vpar1, dist1, par1, upar2, vpar2, dist2, par2;
      Point clo1, clo2;
      ftEdgeBase *e1 = prev_face->closestBoundaryPoint(pos1, dummy_vec, upar1, vpar1,
						       clo1, dist1, par1);
      ftEdgeBase *e2 = prev_face->closestBoundaryPoint(pos2, dummy_vec, upar2, vpar2,
						       clo2, dist2, par2);
      int bd1 = (dist1 < tol) ? 
	((outerbd->isInLoop(e1)) ? 2 : 1) : 0;
      int bd2 = (dist2 < tol) ? 
	 ((outerbd->isInLoop(e2)) ? 2 : 1)  : 0;
      
      params.push_back(std::make_pair(make_pair(pos1,bd1), make_pair(pos2,bd2)));

      // Remove used edges from list
      for (kj=0; kj<chain.size(); ++kj)
	for (kr=ki+1; kr<edges.size(); ++kr)
	  {
	    if (chain[kj] == edges[kr].get() || chain[kj]->twin() == edges[kr].get())
	      edges.erase(edges.begin()+kr);
	  }
    }

  return params;
}

//==========================================================================
ftSurface*
RegularizeFaceSet::identifySeamFaces(shared_ptr<ftSurface> face1,
				     int& pardir, int& status)
//==========================================================================
{
  status = 0;
  ftSurface* face2;
  double epsge = model_->getTolerances().gap;

  // Fetch all edges
  vector<shared_ptr<ftEdge> > edges = face1->getAllEdges();

  // Look for transitions where the underlying surface is
  // the same for both faces
  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      shared_ptr<EdgeVertex> radial_edge = edges[ki]->getEdgeMultiplicityInstance();

      if (!edges[ki]->twin() && !radial_edge.get())
	continue;  // No adjacent face

      // Fetch faces adjacent to this edge
      vector<ftSurface*> faces;
      if (radial_edge.get())
	faces = radial_edge->getAdjacentFaces();
      else
	faces.push_back(edges[ki]->twin()->geomEdge()->face()->asFtSurface());

      // Remove faces belonging to other bodies
      Body *bd = face1->getBody();
      size_t kj;
      for (kj=0; kj<faces.size(); )
	{
	  if (faces[kj]->getBody() != bd)
	    faces.erase(faces.begin()+kj);
	  else 
	    kj++;
	}
      
      for (kj=0; kj<faces.size(); ++kj)
	{
	  face2 = faces[kj];

	  // Test if the faces have the same underlying surface
	  shared_ptr<ParamSurface> surf1 = face1->surface();
	  shared_ptr<ParamSurface> surf2 = face2->surface();
	  shared_ptr<BoundedSurface> bd_sf1 = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
	  shared_ptr<BoundedSurface> bd_sf2 = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
	  if (!(bd_sf1.get() && bd_sf2.get()))
	    continue;  // At least one surface is not trimmed

	  shared_ptr<ParamSurface> base1 = bd_sf1->underlyingSurface();
	  shared_ptr<ParamSurface> base2 = bd_sf2->underlyingSurface();
	  if (base1.get() != base2.get())
	    continue;  // Different underlying surfaces

#ifdef DEBUG_REG
	  std::ofstream of("seam_cand.g2");
	  surf1->writeStandardHeader(of);
	  surf1->write(of);
	  surf2->writeStandardHeader(of);
	  surf2->write(of);
#endif
	  
	  // Count number of adjacencies
	  // NB! This approach may be too strict as it does not take 
	  // into account that the same boundary curves may lead to more
	  // than one edge adjacency
	  int nmb_adjacency = face1->nmbAdjacencies(face2);
	  if (nmb_adjacency != 1)
	    continue;

	  shared_ptr<ftEdge> edge1, edge2;
	  bool neighbours = face1->areNeighbours(face2, edge1, edge2);
	  if (!neighbours)
	    continue;  // Something mysterious

	  int at_seam = getParameterDirection(face1.get(), face2, edge1,
					       edge2, epsge, pardir);
	  if (!at_seam)
	    continue;

	  // Must check for approximate C1 continuity
	  int nsample = 15;
	  double t1 = edge1->tMin();
	  double t2 = edge1->tMax();
	  double tdel = (t2 - t1)/(double)(nsample-1);
	  double tpar, tpar2;
	  int kh;
	  double eps = model_->getTolerances().neighbour;
	  double angtol = model_->getTolerances().kink;
	  double *seed = NULL;
	  for (kh=0, tpar=t1; kh<nsample; kh++, tpar+=tdel)
	    {
	      Point pos1 = edge1->point(tpar);
	      Point pos2;
	      double dist;
	      edge2->closestPoint(pos1, tpar2, pos2, dist, seed);
	      if  (dist > eps)
		break;
	      Point norm1 = edge1->normal(tpar);
	      Point norm2 = edge2->normal(tpar2);
	      if (norm1.angle(norm2) > angtol)
		break;
	      seed = &tpar2;
	    }
	  if (kh < nsample)
	    continue;  // Not sufficent continuity
	  
	  // It remains to check if the underlying surface is closed
	  // in the relevant parameter direction

	  status = at_seam;
	  return face2;
	}
    }

  return 0;  // No match found

}

//==========================================================================
int
RegularizeFaceSet::getParameterDirection(ftSurface* face1, ftSurface* face2,
					 shared_ptr<ftEdge> edge1,
					 shared_ptr<ftEdge> edge2,
					 double eps, int& pardir)
//==========================================================================
{
  // Check for constant parameter curves
  // Fetch information about trimming curves
  int status = 0;
  AdjacencyInfo info = face1->getAdjacencyInfo(face2, eps);

  // Check for adjacent seam surfaces
  ftSurface *f1=NULL, *f2=NULL;
  bool at_seam = false;
  f1 = edge1->twin()->geomEdge()->face()->asFtSurface();
  f2 = edge2->twin()->geomEdge()->face()->asFtSurface();
  if (f1 && f2 && f1->twin()==f2 && f2->twin()==f1)
    at_seam = true;
  else if (edge1->hasEdgeMultiplicity())
    {
      vector<ftSurface*> faces = 
  	edge1->getEdgeMultiplicityInstance()->getAdjacentFaces();
      Body *bd = face1->getBody();
      for (size_t ki=0; ki<faces.size(); ++ki)
  	for (size_t kj=ki+1; kj<faces.size(); ++kj)
  	  if (faces[ki]->twin() == faces[kj] &&
  	      faces[kj]->twin() == faces[ki] &&
  	      faces[ki]->getBody() == bd &&
  	      faces[kj]->getBody() == bd)
  	    at_seam = true;
    }

  if (info.bd_idx_1_ >= 0 && info.bd_idx_2_ >= 0 &&
      /* fabs(info.bd_idx_1_ - info.bd_idx_2_) == 1 && */
      (info.bd_idx_1_/2) == (info.bd_idx_2_/2))
    {
      pardir = (info.bd_idx_1_ <= 1) ? 0 : 1;
      status = 1;
    }
  else if (at_seam)
    {
      status = 2;

      // Check parameter values of associated vertices
      shared_ptr<Vertex> v1, v2, v3, v4;
      edge1->getVertices(v1, v2);
      edge2->getVertices(v3, v4);
      Point p1 = v1->getFacePar(face1); 
      Point p2 = v2->getFacePar(face1);
      Point p3 = v3->getFacePar(face2);
      Point p4 = v4->getFacePar(face2);
      if (fabs(p1[0]-p2[0]) < eps && fabs(p3[0]-p4[0]) < eps)
  	pardir = 1;
      else if (fabs(p1[1]-p2[1]) < eps && fabs(p3[1]-p4[1]) < eps)
  	pardir = 0;
      else
  	return 0;
      
    }
  else
    return 0;

  return status;
}

 //==========================================================================
void
RegularizeFaceSet::removeInsignificantVertices(vector<shared_ptr<Vertex> >& vx)
//==========================================================================
{
  for (size_t ki=0; ki<vx.size(); )
    {
      vector<ftEdge*> edges = vx[ki]->uniqueEdges();

      // Sort edges with respect to associated body
      size_t kj, kr;
      for (kj=0; kj<edges.size(); ++kj)
	{
	  Body *body = edges[kj]->face()->asFtSurface()->getBody();
	  for (kr=kj+1; kr<edges.size(); ++kr)
	    {
	      if (edges[kr]->face()->asFtSurface()->getBody() == body)
		std::swap(edges[kj+1], edges[kr]);
	    }
	}
      
      bool non_sign = true;
      bool seam = true;
      for (size_t kj=0; kj<edges.size(); kj=kr)
	{
	  Body *body = edges[kj]->face()->asFtSurface()->getBody();
	  for (kr=kj+1; kr<edges.size(); ++kr)
	    if (edges[kr]->face()->asFtSurface()->getBody() != body)
	      break;

	  // Check if the vertex is insignificant in this body
	  if (kr - kj == 2)
	    {
	      // Check if the underlying geometry curve passes this
	      // vertex
	      shared_ptr<ParamCurve> cv1 = edges[kj]->geomCurve();
	      shared_ptr<ParamCurve> cv2;
	      if (edges[kj]->twin())
		cv2 = edges[kj]->twin()->geomEdge()->geomCurve();
	      shared_ptr<ParamCurve> cv3 = edges[kj+1]->geomCurve();
	      shared_ptr<ParamCurve> cv4;
	      if (edges[kj+1]->twin())
		cv4 = edges[kj+1]->twin()->geomEdge()->geomCurve();
	      if (!(cv2.get() && cv4.get()))
		{
		  if (cv1.get() != cv3.get())
		    {
		      // Check if the vertex is a corner
		      double t1 = edges[kj]->parAtVertex(vx[ki].get());
		      double t2 = edges[kj+1]->parAtVertex(vx[ki].get());
		      Point tan1 = edges[kj]->tangent(t1);
		      Point tan2 = edges[kj+1]->tangent(t2);
		      if (tan1.angle(tan2) >= model_->getTolerances().bend)
			non_sign = false;
		    }
		}
	      else if (cv2.get() && cv4.get())
		{
		  if (!((cv1.get() == cv3.get() && cv2.get() == cv4.get()) ||
			(cv1.get() == cv4.get() && cv2.get() == cv3.get())))
		    non_sign = false;
		}
	      else
		non_sign = false;
	    }
	  else
	    non_sign = false;

	  if (non_sign)
	    continue;   // Check next body

	  // Check if the vertex lies at a seam in this body
	  // Fetch associated faces in this body
	  set<ftSurface*> faces;
	  for (size_t kh=kj; kh<kr; ++kh)
	    {
	      faces.insert(edges[kh]->face()->asFtSurface());
	      if (edges[kh]->twin())
		faces.insert(edges[kh]->twin()->geomEdge()->face()->asFtSurface());
	    }
	  if (faces.size() != 2)
	    seam = false;

	  if (!seam)
	    {
	      // Check if this vertex is removed previously
	      size_t kr;
	      Point vx_pos = vx[ki]->getVertexPoint();
	      for (kr=0; kr<seam_joints_.size(); ++kr)
		if (seam_joints_[kr].dist(vx_pos) < model_->getTolerances().gap)
		  break;
	      if (kr < seam_joints_.size())
		seam = true;
	      else
		break;
	    }
	}

      if (non_sign || seam)
	vx.erase(vx.begin()+ki);
      else 
	ki++;
    }
}

//==========================================================================
void
RegularizeFaceSet::seamClassification()
//==========================================================================
{
  double angtol = model_->getTolerances().kink;
  std::set<shared_ptr<Vertex> > vxs;
  int nmb = model_->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      vector<shared_ptr<Vertex> > curr_vx =
	model_->getFace(ki)->getNonCornerVertices(angtol);
      vxs.insert(curr_vx.begin(), curr_vx.end());
    }
  vector<shared_ptr<Vertex> > vx;
  vx.insert(vx.end(), vxs.begin(), vxs.end());

  for (size_t ki=0; ki<vx.size(); ++ki)
    {
      vector<ftEdge*> edges = vx[ki]->uniqueEdges();

      // Sort edges with respect to associated body
      size_t kj, kr;
      for (kj=0; kj<edges.size(); ++kj)
	{
	  Body *body = edges[kj]->face()->asFtSurface()->getBody();
	  for (kr=kj+1; kr<edges.size(); ++kr)
	    {
	      if (edges[kr]->face()->asFtSurface()->getBody() == body)
		std::swap(edges[kj+1], edges[kr]);
	    }
	}
      
      bool non_sign = true;
      bool seam = true;
      for (size_t kj=0; kj<edges.size(); kj=kr)
	{
	  Body *body = edges[kj]->face()->asFtSurface()->getBody();
	  for (kr=kj+1; kr<edges.size(); ++kr)
	    if (edges[kr]->face()->asFtSurface()->getBody() != body)
	      break;

	  // Check if the vertex is insignificant in this body
	  if (kr - kj == 2)
	    {
	      // Check if the underlying geometry curve passes this
	      // vertex
	      shared_ptr<ParamCurve> cv1 = edges[kj]->geomCurve();
	      shared_ptr<ParamCurve> cv2;
	      if (edges[kj]->twin())
		cv2 = edges[kj]->twin()->geomEdge()->geomCurve();
	      shared_ptr<ParamCurve> cv3 = edges[kj+1]->geomCurve();
	      shared_ptr<ParamCurve> cv4;
	      if (edges[kj+1]->twin())
		cv4 = edges[kj+1]->twin()->geomEdge()->geomCurve();
	      if (!(cv2.get() && cv4.get()))
		{
		  if (cv1.get() != cv3.get())
		    non_sign = false;
		}
	      else if (cv2.get() && cv4.get())
		{
		  if (!((cv1.get() == cv3.get() && cv2.get() == cv4.get()) ||
			(cv1.get() == cv4.get() && cv2.get() == cv3.get())))
		    non_sign = false;
		}
	      else
		non_sign = false;
	    }
	  else
	    non_sign = false;

	  if (non_sign)
	    continue;   // Check next body

	  // Check if the vertex lies at a seam in this body
	  // Fetch associated faces in this body
	  set<ftSurface*> faces;
	  size_t kh, kw;
	  for (kh=kj; kh<kr; ++kh)
	    {
	      faces.insert(edges[kh]->face()->asFtSurface());
	      if (edges[kh]->twin())
		faces.insert(edges[kh]->twin()->geomEdge()->face()->asFtSurface());
	    }
	  vector<ftSurface*> faces2;
	  faces2.insert(faces2.end(), faces.begin(), faces.end());
	  // Remove twins
	  for (kh=0; kh<faces2.size(); )
	    {
	      for (kw=kh+1; kw<faces2.size(); ++kw)
		if (faces2[kh]->twin() == faces2[kw] &&
		    faces2[kw]->twin() == faces2[kh])
		  break;
	      if (kw < faces2.size())
		{
		  faces2.erase(faces2.begin()+kw);
		  faces2.erase(faces2.begin()+kh);
		}
	      else
		kh++;
	    }
	  if (faces2.size() != 2)
	    seam = false;

	  if (!seam)
	    break;
	}
      
      if (non_sign)
	seam = false;
      if (seam)
	seam_joints_.push_back(vx[ki]->getVertexPoint());
    }
}

//==========================================================================
void
RegularizeFaceSet::getSeamRadialEdge(ftSurface* face, 
				     vector<shared_ptr<EdgeVertex> >& edgevx, 
				     vector<std::pair<Point,Point> >& endpts)
//==========================================================================
{
  // Fetch all radial edges
  vector<shared_ptr<EdgeVertex> > radedges = face->getRadialEdges();

  // Current body
  Body *bd = face->getBody();

  for (size_t ki=0; ki<radedges.size(); ++ki)
    {
      vector<ftEdge*> edges = radedges[ki]->allEdges();
      size_t kj;
      for (kj=0; kj<edges.size(); ++kj)
	{
	  ftSurface *curr_face = edges[kj]->geomEdge()->face()->asFtSurface();
	  Body *curr_bd = curr_face->getBody();
	  if (curr_bd == bd && curr_face != face)
	    break;  // More faces are meeting in the current edge
	}

      if (kj == edges.size())
	{
	  // Due to a possible inconsistency in the radial edge, the
	  // edge must be limited by end points
	  for (kj=0; kj<edges.size(); ++kj)
	    if (edges[kj]->geomEdge()->face() == face)
	      break;

	  if (kj < edges.size())
	    {
	      // A candidate radial edge to get lost
	      edgevx.push_back(radedges[ki]);
	      
	      Point pt1 = edges[kj]->point(edges[kj]->tMin());
	      Point pt2 = edges[kj]->point(edges[kj]->tMax());
	      endpts.push_back(make_pair(pt1,pt2));
	    }
	}
    }
}	  
      
//==========================================================================
void
RegularizeFaceSet::attachRadialEdge(ftSurface* face, 
				    vector<shared_ptr<EdgeVertex> >& edgevx, 
				    vector<std::pair<Point,Point> >& endpts,
				    double tol)
//==========================================================================
{
  // For each face edge that does not have a radial edge, check if
  // it fits with a given radial edge
  vector<shared_ptr<ftEdge> > edges = face->getAllEdges();
  size_t ki, kj, kh;
  for (kj=0; kj<edges.size(); ++kj)
    {
      if (edges[kj]->hasEdgeMultiplicity())
	continue;

      Point pt1 = edges[kj]->point(edges[kj]->tMin());
      Point pt2 = edges[kj]->point(edges[kj]->tMax());
      for (ki=0; ki<edgevx.size(); ++ki)
	{
	  if (std::min(pt1.dist(endpts[ki].first), pt1.dist(endpts[ki].second)) > tol ||
	      std::min(pt2.dist(endpts[ki].first), pt2.dist(endpts[ki].second)) > tol)
	    continue;  // Not corresponding end points

	  // Check if the radial edge contains any new information
	  vector<ftEdge*> edges2 = edgevx[ki]->allEdges();
	  for (kh=0; kh<edges2.size(); ++kh)
	    if (edges2[kh]->face() != face)
	      break;
	  if (kh < edges2.size())
	    {
	      // Connect to radial edge
	      edgevx[ki]->addEdge(edges[kj].get());
	      edges[kj]->setEdgeVertex(edgevx[ki]);
	      break;
	    }
	}
    }
}

//==========================================================================
int
RegularizeFaceSet::mergeSituation(ftSurface* face, 
				  shared_ptr<Vertex> vx,
				  ftSurface*& merge1,
				  int& dir1, double& val1, bool& atstart1,
				  ftSurface*& merge2,
				  int& dir2, double& val2, bool& atstart2,
				  shared_ptr<Vertex>& other_vx,
				  pair<Point, Point>& co_par1, 
				  pair<Point, Point>& co_par2)
//==========================================================================
{
  double angtol = model_->getTolerances().bend;
  int found_merge_cand = RegularizeUtils::noExtension(vx, face,
						      other_vx, co_par1,
						      co_par2, dir1,
						      dir2, val1, val2,
						      angtol, true);
  // bool found_merge_cand = RegularizeUtils::noExtension(vx, face,
  // 						       other_vx, co_par1,
  // 						       co_par2, dir1,
  // 						       dir2, val1, val2,
  // 						       angtol, false);
  if (found_merge_cand == 0)
    return found_merge_cand;  // No merge candiates are found 

  // Fetch faces

  ftEdge *edge = vx->getCommonEdge(other_vx.get());;
  merge1 = edge->face()->asFtSurface();
  merge2 = edge->twin()->face()->asFtSurface();

  if (found_merge_cand < 2)
    return found_merge_cand;  // The common boundary is not a
  //  constant parameter curve in both faces

  // Check whether the merge is at start- or end of the surfaces
  double u1, u2, v1, v2;
  merge1->surface()->getInternalPoint(u1, v1);
  merge2->surface()->getInternalPoint(u2, v2);
  atstart1 = (dir1 == 0) ? (u1 > val1) : (v1 > val1);
  atstart2 = (dir2 == 0) ? (u2 > val2) : (v2 > val2);
      
  return 2;
}

//==========================================================================
void
RegularizeFaceSet::prioritizeFaces(vector<shared_ptr<ftSurface> >& faces,
				   vector<int>& perm)
//==========================================================================
{
  // We start by a simple selection where faces with holes, but no corners in
  // the outer loop is handled last
  size_t ki, kj, kr;
  double kink = model_->getTolerances().bend;
  for (ki=0, kj=faces.size()-1; ki<kj; )
    {
      int nmb_loops = faces[perm[ki]]->nmbBoundaryLoops();
      if (nmb_loops <= 1)
  	{
  	  ki++;
  	  continue;
  	}
      vector<shared_ptr<Vertex> > corners = 
  	faces[perm[ki]]->getCornerVertices(kink, 0);
      if (nmb_loops == 2 && corners.size() == 0)
  	{
  	  // Postpone splitting of this face
	  std::swap(perm[ki], perm[kj]);
  	  kj--;
  	}
      else
  	ki++;
    }

  // The second last face type to handle is the ones without holes and
  // with 4 corners
  for (ki=0; ki<kj; )
    {
      int nmb_loops = faces[perm[ki]]->nmbBoundaryLoops();
      if (nmb_loops > 1)
	{
	  ki++;
	  continue;
	}
      vector<shared_ptr<Vertex> > corners = 
	faces[perm[ki]]->getCornerVertices(kink, 0);
      if (nmb_loops == 1 && corners.size() == 4)
	{
	  // Postpone splitting of this face
	  std::swap(perm[ki], perm[kj]);
	  kj--;
	}
      else
	ki++;
    }

  // Treat adjacent faces to faces with hole directly after their
  // neighbour
  int kh = 1;
  for (ki=0; ki<faces.size(); ki+=kh)
    {
      kh = 1;
      int nmb_loops = faces[perm[ki]]->nmbBoundaryLoops();
      if (nmb_loops > 1)
	{
	  // Check if there is a corresponding face and this face has
	  // more loops than the current one
	  for (kr=0; kr<corr_faces_.size(); ++kr)
	    if (corr_faces_[kr].first == perm[ki] ||
		corr_faces_[kr].second == perm[ki])
	      break;
	  if (kr < corr_faces_.size())
	    {
	      size_t ix = (corr_faces_[kr].first == perm[ki]) ? corr_faces_[kr].second :
		corr_faces_[kr].first;
	      for (kr=0; kr<perm.size(); ++kr)
		if (perm[kr] == (int)ix)
		  break;
	      if (kr >= perm.size())
		kr = ki;
	      int nmb_loops2 = faces[ix]->nmbBoundaryLoops();
	      if (kr > ki && nmb_loops2 > nmb_loops)
		{
		  std::swap(perm[ki], perm[kr]);
		  nmb_loops = nmb_loops2;
		}
	    }
	  
	  vector<ftSurface*> neighbours;
	  faces[perm[ki]]->getAdjacentFaces(neighbours);
	  for (kr=0; kr<neighbours.size(); ++kr)
	    {
	      for (kj=ki+kh; kj<faces.size(); ++kj)
		if (faces[perm[kj]].get() == neighbours[kr])
		  break;
	      if (kj < faces.size())
		{
		  int tmp = perm[kj];
		  perm.insert(perm.begin()+ki+1, tmp);
		  perm.erase(perm.begin()+kj+1);
		  kh++;
		}
	    }
	}
    }
}

//==========================================================================
void
RegularizeFaceSet::computeFaceCorrespondance(vector<shared_ptr<ftSurface> >& faces)
//==========================================================================
{
  // For each combination of faces, check if they are next neighbours through
  // a number of (at least three) other faces. This will not catch all occurances
  // of faces correspondance
  int ki, kj;
  size_t kr, kh;

  double ang_tol = 0.1;
  vector<double> corr_dist;
  vector<double> len_frac;

  // Collect information
  vector<DirectionCone> normalcone(faces.size());
  vector<DirectionCone> tangcone1(faces.size());
  vector<DirectionCone> tangcone2(faces.size());
  vector<BoundingBox> box(faces.size());
  vector<Point> mid(faces.size());
  for (kr=0; kr<faces.size(); ++kr)
    {
      shared_ptr<ParamSurface> sf = faces[kr]->surface();
      normalcone[kr] = sf->normalCone();
      tangcone1[kr] = sf->tangentCone(true);
      tangcone2[kr] = sf->tangentCone(false);
      box[kr] = sf->boundingBox();
      mid[kr] = 0.5*(box[kr].low() + box[kr].high());
    }

  for (ki=0; ki<(int)faces.size(); ++ki)
    {
      // Check if the face already occurs in a face correspondance pair
      for (kr=0; kr<corr_faces_.size(); ++kr)
	if (corr_faces_[kr].first == ki || corr_faces_[kr].second == ki)
	  break;
      // if (kr < corr_faces_.size())
      // 	continue;

      for (kj=ki+1; kj<(int)faces.size(); ++kj)
	{
	  // Check if the face already occurs in a face correspondance pair
	  for (kh=0; kh<corr_faces_.size(); ++kh)
	    if (corr_faces_[kh].first == kj || corr_faces_[kh].second == kj)
	      break;
	  // if (kr < corr_faces_.size())
	  //   continue;

#ifdef DEBUG_REG
	  std::ofstream of("curr_corr_face.g2");
	  faces[ki]->surface()->writeStandardHeader(of);
	  faces[ki]->surface()->write(of);
	  faces[kj]->surface()->writeStandardHeader(of);
	  faces[kj]->surface()->write(of);
#endif
	
	  // Check if the faces are neighbours
	  if (faces[ki]->nmbAdjacencies(faces[kj].get()) > 0)
	    continue;  // No correspondance

	  // Check if the faces may be opposite to each other
	  Point axis1 = normalcone[ki].centre();
	  Point axis2 = normalcone[kj].centre();
	  double cone_angle = axis1.angle(axis2);
	  cone_angle = M_PI - cone_angle;
	  // if (axis1*axis2 > -1.0e-6)
	  //   cone_angle = 2*M_PI-cone_angle;
	  if (cone_angle < ang_tol)
	    {
	      // Near parallel and opposite normal cone axes.
	      // Possibly corresponding faces
	      // Make distance estimate
	      double dist = mid[ki].dist(mid[kj]);
	      double len1 = box[ki].high().dist(box[ki].low());
	      double len2 = box[kj].high().dist(box[kj].low());
	      double frac = std::min(len1,len2)/std::max(len1,len2);
	      if ((kr == corr_faces_.size() || len_frac[kr]*dist < frac*corr_dist[kr]) &&
		  (kh == corr_faces_.size() || len_frac[kh]*dist < frac*corr_dist[kh]))
		{
		  // A match is found
		  if (kr < corr_faces_.size() && len_frac[kr]*dist < frac*corr_dist[kr])
		    {
		      corr_faces_[kr] = make_pair(ki, kj);
		      corr_dist[kr] = dist;
		      len_frac[kr] = frac;
		    }
		  else if (kh < corr_faces_.size() && len_frac[kh]*dist < frac*corr_dist[kh])
		    {
		      corr_faces_[kh] = make_pair(ki, kj);
		      corr_dist[kh] = dist;
		      len_frac[kh] = frac;
		    }
		  else 
		    {
		       setFaceCorrespondance(ki, kj);
		       corr_dist.push_back(dist);
		       len_frac.push_back(frac);
		    }
		}
	    }
	  else
	    {
	      // Alternative approach
	      int nmb_next = faces[ki]->nmbNextNeighbours(faces[kj].get());
	      if (nmb_next >= 3)
		{
		  // Check if a closest point internal to one face hits the other
		  Point pos1, pos2;
		  Point norm1, norm2, vec;
		  double u1, v1, u2, v2, dist;
		  pos1 = faces[ki]->surface()->getInternalPoint(u1, v1);
		  norm1 = faces[ki]->normal(u1, v1);

		  faces[kj]->closestPoint(pos1, u2, v2, pos2, dist, 
					  model_->getTolerances().gap);
		  norm2 = faces[kj]->normal(u2, v2);
		  vec = pos2 - pos1;

		  if (vec.length() > model_->getTolerances().neighbour)
		    {
		      double ang1 = vec.angle(norm1);
		      double ang2 = vec.angle(norm2);
		      double ang_tol = 0.01*M_PI;
		      if ((ang1 < ang_tol || M_PI-ang1 < ang_tol) && 
			  (ang2 < ang_tol || M_PI-ang2 < ang_tol))
			{
			  // A correspondance is found
			  double dist = mid[ki].dist(mid[kj]);
			  double len1 = box[ki].high().dist(box[ki].low());
			  double len2 = box[kj].high().dist(box[kj].low());
			  double frac = std::min(len1,len2)/std::max(len1,len2);
			  if (kr == corr_faces_.size() || dist < corr_dist[kr])
			    {
			      if (kr < corr_faces_.size() && dist < corr_dist[kr])
				{
				  corr_faces_[kr] = make_pair(ki, kj);
				  corr_dist[kr] = dist;
				  len_frac[kr] = frac;
				}
			      else if (kh < corr_faces_.size() && dist < corr_dist[kh])
				{
				  corr_faces_[kh] = make_pair(ki, kj);
				  corr_dist[kh] = dist;
				  len_frac[kh] = frac;
				}
			      else 
				{
				  setFaceCorrespondance(ki, kj);
				  corr_dist.push_back(dist);
				  len_frac.push_back(frac);
				}
			    }
			}
		    }
		}
	    }
	}
    } 
}

}  // namespace Go
