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

#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/VolumeAdjacency.h"
#include "GoTools/tesselator/GeneralMesh.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariate/ElementaryVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/trivariate/VolumeTools.h"
#include <fstream>

//#define DEBUG
//#define DEBUG_VOL2

using namespace Go;
using std::vector;

//===========================================================================
VolumeModel::VolumeModel(std::vector<shared_ptr<ftVolume> >& volumes,
			 double space_epsilon,
			 double kink)  // Kink between adjacent surfaces 
  //===========================================================================
  : CompositeModel(space_epsilon, 10.0*space_epsilon, kink, 10.0*kink),
    approxtol_(space_epsilon)
{
  if (volumes.size() == 0)
    return;

  bodies_.reserve(volumes.size());
  for (size_t ki=0; ki<volumes.size(); ++ki)
    bodies_.push_back(volumes[ki]);

  buildTopology();
}

 
//===========================================================================
VolumeModel::VolumeModel(std::vector<shared_ptr<ftVolume> >& volumes,
			 double space_epsilon, double neighbour, 
			 double kink, double bend,  // Kink between adjacent surfaces 
			 bool adjacency_set)
  //===========================================================================
  : CompositeModel(space_epsilon, neighbour, kink, bend),
    approxtol_(space_epsilon)
{
  if (volumes.size() == 0)
    return;

  bodies_.reserve(volumes.size());
  for (size_t ki=0; ki<volumes.size(); ++ki)
    bodies_.push_back(volumes[ki]);

  if (adjacency_set)
    {
      setBoundarySfs();
      setVertexIdentity();
    }
  else
    buildTopology();
}

 
//===========================================================================
VolumeModel::VolumeModel(double space_epsilon, double neighbour, 
			 double kink, double bend)  // Kink between adjacent surfaces 
  //===========================================================================
  : CompositeModel(space_epsilon, neighbour, kink, bend),
    approxtol_(space_epsilon)
{
}

//===========================================================================
VolumeModel::VolumeModel(const VolumeModel& vm)
  //===========================================================================
  : CompositeModel(vm.toptol_.gap, vm.toptol_.neighbour, vm.toptol_.kink, 
		   vm.toptol_.bend),
    approxtol_(vm.approxtol_)
{
}

//===========================================================================
VolumeModel::~VolumeModel()
//===========================================================================
{
}


//===========================================================================
// Bounding box of the entire volume model
BoundingBox VolumeModel::boundingBox()
//===========================================================================
{
  BoundingBox box;
  if (bodies_.size() == 0)
    return box;

  box = bodies_[0]->boundingBox();
  for (size_t ki=1; ki<bodies_.size(); ++ki)
    box.addUnionWith(bodies_[ki]->boundingBox());
  return box;
}


//===========================================================================
// Bounding box of one volume
BoundingBox VolumeModel::boundingBox(int idx) const
//===========================================================================
{
  ASSERT(idx < (int)bodies_.size());

  BoundingBox box = bodies_[0]->boundingBox();
  for (size_t ki=1; ki<bodies_.size(); ++ki)
    {
      box.addUnionWith( bodies_[ki]->boundingBox());
    }
  return box;
}


//===========================================================================
int VolumeModel::nmbEntities() const
//===========================================================================
{
  return (int)bodies_.size();
}

//===========================================================================
shared_ptr<ftVolume> VolumeModel::getBody(int idx) const
//===========================================================================
{
  shared_ptr<ftVolume> result;
  if (idx >= 0 && idx < (int)bodies_.size())
    result = static_pointer_cast<ftVolume>(bodies_[idx]);
  return result;
}


//===========================================================================
shared_ptr<ParamVolume> VolumeModel::getVolume(int idx) const
//===========================================================================
{
  shared_ptr<ftVolume> body = getBody(idx);
  if (body.get() == 0)
    {
      shared_ptr<ParamVolume> dummy;
      return dummy;
    }
  return body -> getVolume();
}


//===========================================================================
shared_ptr<SplineVolume> VolumeModel::getSplineVolume(int idx) const
//===========================================================================
{
  shared_ptr<ParamVolume> par_vol = getVolume(idx);
  shared_ptr<SplineVolume> spl_vol;
  if (par_vol.get() == 0)
    return spl_vol;

  spl_vol = dynamic_pointer_cast<SplineVolume, ParamVolume>(par_vol);
  if (spl_vol.get() != 0)
    return spl_vol;

  shared_ptr<ElementaryVolume> elem_vol = dynamic_pointer_cast<ElementaryVolume, ParamVolume>(par_vol);
  if (elem_vol.get() == 0)
    return spl_vol;

  return shared_ptr<SplineVolume>(elem_vol->geometryVolume());
}


//===========================================================================
int VolumeModel::getIndex(shared_ptr<ftVolume> body) const
//===========================================================================
{
  return getIndex(body.get());
}


//===========================================================================
int VolumeModel::getIndex(ftVolume* body) const
//===========================================================================
{
  for (size_t i = 0; i < bodies_.size(); ++i)
    if (bodies_[i].get() == body)
	return (int)i;

  return -1;
}


//===========================================================================
shared_ptr<ftVolume> VolumeModel::fetchAsSharedPtr(Body *body) const
//===========================================================================
{
  shared_ptr<ftVolume> result;
  for (size_t i = 0; i < bodies_.size(); ++i)
    if (bodies_[i].get() == body)
      {
	result = static_pointer_cast<ftVolume>(bodies_[i]);
	break;
      }
  return result;
}


//===========================================================================
void VolumeModel::evaluate(int idx,      // Index of volume
			   double par[], // Parameter value
			   Point& pnt) const
//===========================================================================
{
  bodies_[idx]->getVolume()->point(pnt, par[0], par[1], par[2]);
}

//===========================================================================
void VolumeModel:: evaluate(int idx,      // Index
			    double par[], // Parameter value
			    int nder,     // Number of derivatives to compute, 0=only position
			    std::vector<Point>& der) const
//===========================================================================
{
  bodies_[idx]->getVolume()->point(der, par[0], par[1], par[2], nder);
}

//===========================================================================
void VolumeModel::closestPoint(Point& pnt,     // Input point
		 Point& clo_pnt, // Found closest point
		 int& idx,           // Index of volume where the closest point is found
		 double clo_par[],   // Parameter value corrsponding to the closest point
		 double& dist)  // Distance between input point and found closest point
//===========================================================================
{
  // Not implemented
}

//===========================================================================
shared_ptr<IntResultsModel> VolumeModel::intersect(const ftLine& line)
//===========================================================================
{
  // Not implemented
    MESSAGE("Not implemented");
    return shared_ptr<IntResultsModel>();
}

//===========================================================================
shared_ptr<IntResultsModel> VolumeModel::intersect_plane(const ftPlane& plane)
//===========================================================================
{
  // Not implemented
    MESSAGE("Not implemented");
    return shared_ptr<IntResultsModel>();
}

//===========================================================================
void VolumeModel::extremalPoint(Point& dir,     // Direction
		  Point& clo_pnt, // Found closest point
		  int& idx,           // Index of volume where the closest point is found
		  double ext_par[])   // Parameter value of extremal point
//===========================================================================
{
  // Not implemented
}

//===========================================================================
bool VolumeModel::isDegenerate(int idx) const
//===========================================================================
{
    // Not implemented
    MESSAGE("Not implemented - returning arbitrary value");
    return false;
}

//===========================================================================
double VolumeModel::curvature(int idx, // Index of entity
			      double *par) const  // Parameter value at which to compute curvature
//===========================================================================
{
  // Not implemented
    MESSAGE("Not implemented - returning arbitrary value");
    return 0.0;
}

//===========================================================================
void VolumeModel:: turn(int idx)  // Turn parameter directions of one entity
//===========================================================================
{
  // Not implemented
}

//===========================================================================
void VolumeModel::turn()
//===========================================================================
{
  // Not implemented
}

//===========================================================================
void VolumeModel::append(shared_ptr<ftVolume> volume)
//===========================================================================
{
// #ifdef DEBUG
//   bool isOK = checkModelTopology();
//   if (!isOK)
//     std::cout << "VolumeTopology, append (before). Topology inconsistencies" << std::endl;
// #endif

  bodies_.push_back(volume);
  buildTopology(volume);

  boundary_shells_.clear();
  setBoundarySfs();

#ifdef DEBUG
  bool isOK = checkModelTopology();
  if (!isOK)
    std::cout << "VolumeTopology, append (after). Topology inconsistencies" << std::endl;
#endif

 }

//===========================================================================
void VolumeModel::append(vector<shared_ptr<ftVolume> > volumes)
//===========================================================================
{
  for (vector<shared_ptr<ftVolume> >::const_iterator it = volumes.begin();
       it != volumes.end();
       ++it)
    append(*it);
}

//===========================================================================
void VolumeModel::append(shared_ptr<VolumeModel> anotherModel)
  //===========================================================================
{
  append(anotherModel->bodies_);
}

//===========================================================================
void VolumeModel::removeSolid(shared_ptr<ftVolume> vol)
  //===========================================================================
{
#ifdef DEBUG
  bool isOK = checkModelTopology();
  if (!isOK)
    std::cout << "VolumeTopology, removeSolid (before). Topology inconsistencies" << std::endl;
#endif

  // Find index
  int idx = getIndex(vol.get());
  if (idx < 0)
    return; // Does not exist

  // Remove the boundary faces one by one to erase the data structure
  int nmb_shells = vol->nmbOfShells();
  for (int ki=0; ki<nmb_shells; ki++)
    {
      shared_ptr<SurfaceModel> sfmodel = vol->getShell(ki);
      int nmb_faces = sfmodel->nmbEntities();
      for (int kj=nmb_faces-1; kj>=0; --kj)
	{
	  shared_ptr<ftSurface> face = sfmodel->getFace(kj);
	  face->disconnectTwin();
	  sfmodel->removeFace(face);
	}
    }
  bodies_.erase(bodies_.begin() + idx);

  // Regenerate model boundaries
  boundary_shells_.clear();
  setBoundarySfs();
  
#ifdef DEBUG
  isOK = checkModelTopology();
  if (!isOK)
    std::cout << "VolumeTopology, removeSolid (after). Topology inconsistencies" << std::endl;
#endif

}

//===========================================================================
void VolumeModel::tesselate(vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  // Not implemented
}

//===========================================================================
void VolumeModel::tesselate(int resolution[],
			    vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  // Not implemented
}

//===========================================================================
void VolumeModel::tesselate(double density,
			    vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  // Not implemented
}
//===========================================================================
void VolumeModel::tesselatedCtrPolygon(vector<shared_ptr<LineCloud> >& ctr_pol) const
//===========================================================================
{
  // Not implemented
}


//===========================================================================
void VolumeModel::buildTopology()
//===========================================================================
{
  VolumeAdjacency computeTop(toptol_.gap, toptol_.neighbour);
  vector<shared_ptr<Body> > solids(bodies_.begin(), bodies_.end());
  computeTop.setAdjacency(solids);

  setBoundarySfs();

  // Add information about faces at the boundary meeting only in radial edges
  setVertexIdentity();
}


//===========================================================================
void VolumeModel::buildTopology(shared_ptr<ftVolume> body)
//===========================================================================
{
  int body_idx;
  for (body_idx = 0; body_idx < (int)bodies_.size(); ++body_idx)
    if (bodies_[body_idx] == body)
      break;

  if (body_idx < (int)bodies_.size())
    {
      VolumeAdjacency computeTop(toptol_.gap, toptol_.neighbour);
      vector<shared_ptr<Body> > solids(bodies_.begin(), bodies_.end());
      computeTop.setAdjacency(solids, body_idx);
    }

  // Add information about faces at the boundary meeting only in radial edges
  setVertexIdentity();
}

  //===========================================================================
  void VolumeModel::setVertexIdentity()
  //---------------------------------------------------------------------------
  //
  // Purpose: Add information about faces at the boundary meeting only 
  //          in radial edges
  //
  //===========================================================================
  {
    // Fetch all faces lying at outer boundaries, i.e. all faces with no twin
    vector<shared_ptr<ftSurface> > bd_faces = getBoundaryFaces();

    // Fetch all boundary vertices
    size_t ki, kj;
    std::set<shared_ptr<Vertex> > bd_vx0;
    for (ki=0; ki<bd_faces.size(); ++ki)
      {
	vector<shared_ptr<Vertex> > curr_vx = bd_faces[ki]->vertices();
	bd_vx0.insert(curr_vx.begin(), curr_vx.end());
      }
    vector<shared_ptr<Vertex> > bd_vx;
    bd_vx.insert(bd_vx.end(), bd_vx0.begin(), bd_vx0.end());
    
    // Find identitical vertices
    vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > cand_vx;
    for (ki=0; ki<bd_vx.size(); ++ki)
      for (kj=ki+1; kj<bd_vx.size(); ++kj)
	{
	  double dist = 
	    bd_vx[ki]->getVertexPoint().dist(bd_vx[kj]->getVertexPoint());
	  if (dist < toptol_.neighbour)
	    {
	      pair<shared_ptr<Vertex>, shared_ptr<Vertex> > ident = 
			std::make_pair(bd_vx[ki], bd_vx[kj]);
	      cand_vx.push_back(ident);
	    }
	}

    // Find pairs of identical vertices that could be joined with 
    // an edge vertex
    for (ki=0; ki<cand_vx.size(); ++ki)
      {
	for (kj=ki+1; kj<cand_vx.size(); ++kj)
	  {
	    // Find joining edges, if any
	    ftEdge *e1 = 
	      cand_vx[ki].first->getCommonEdge(cand_vx[kj].first.get());
	    ftEdge *e2 = 
	      cand_vx[ki].first->getCommonEdge(cand_vx[kj].second.get());
	    ftEdge *e3 = 
	      cand_vx[ki].second->getCommonEdge(cand_vx[kj].first.get());
	    ftEdge *e4 = 
	      cand_vx[ki].second->getCommonEdge(cand_vx[kj].second.get());
	    if (e1 && e3)
	      {
		e1->addEdgeMultiplicityInstance(e3);
	      }
	    else if (e1 && e4)
	      {
		e1->addEdgeMultiplicityInstance(e4);
	      }
	    else if (e2 && e3)
	      {
		e2->addEdgeMultiplicityInstance(e3);
	      }
	    else if (e2 && e4)
	      {
		e2->addEdgeMultiplicityInstance(e4);
	      }
	    else
	      continue;  // Cannot join with radial edge
	    
	  }
      }
  }

//===========================================================================
void VolumeModel::setBoundarySfs()
//===========================================================================
{
  // Fetch all faces lying at outer boundaries, i.e. all faces with no twin
  vector<shared_ptr<ftSurface> > bd_faces = getBoundaryFaces();

#ifdef DEBUG_VOL2
   std::ofstream of("bd_sfs.g2");
#endif

  // Make copy of all faces to avoid destroying existing topology 
  // information in the shells of the solids
  vector<shared_ptr<ftSurface> > bd_faces2(bd_faces.size());
  for (size_t ki=0; ki<bd_faces.size(); ++ki)
    {
      bd_faces2[ki] = shared_ptr<ftSurface>(new ftSurface(bd_faces[ki]->surface(),
							  (int)ki));
      bd_faces2[ki]->setBody(bd_faces[ki]->getBody());

 #ifdef DEBUG_VOL2
     shared_ptr<ParamSurface> surf = bd_faces[ki]->surface();
      shared_ptr<SurfaceOnVolume> vsurf = 
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf);
      if (vsurf.get())
	{
	  vsurf->spaceSurface()->writeStandardHeader(of);
	  vsurf->spaceSurface()->write(of);
	}
      else
	{
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
#endif
    }

  // Represent the outer boundaries as a SurfaceModel
  shared_ptr<SurfaceModel> sfmodel =
    shared_ptr<SurfaceModel>(new SurfaceModel(toptol_.gap, toptol_.gap,
					      toptol_.neighbour, 
					      toptol_.kink, toptol_.bend,
					      bd_faces2));

  // Fetch all connected models
  vector<shared_ptr<SurfaceModel> > models = sfmodel->getConnectedModels();

  if (models.size() == 1)
    boundary_shells_.push_back(models);
  else
    {
      // Sort models into connected submodels and inner and outer shells
      // First get maximum size of model
      BoundingBox box = boundingBox();
      
      // Classify shells
      for (size_t ki=0; ki<models.size(); ++ki)
	{
	  // Check if the current model is already classified
	  size_t kj, kr;
	  bool found = findBoundaryShell(models[ki], kj, kr);
	  if (found)
	    continue;  // Shell already classified

	  if (models[ki]->nmbEntities() == 0)
	    continue;  // Empty model
	  
	  // Make linear curve through current model
	  shared_ptr<ftSurface> face = models[ki]->getFace(0);
	  
	  // Fetch a point on this face
	  Point pnt, norm;
	  shared_ptr<ParamSurface> srf = face->surface();
	  RectDomain dom = srf->containingDomain();
	  double upar = 0.5*(dom.umin()+dom.umax());
	  double vpar = 0.5*(dom.vmin()+dom.vmax());
	  Point pnt0 = srf->point(upar, vpar);
	  double dist;
	  face->closestPoint(pnt0, upar, vpar, pnt, dist, toptol_.gap);
	  norm = face->normal(upar, vpar);

	  // Make large enough curve
	  double len = box.low().dist(box.high()); // Estimated box size
	  norm.normalize();
	  Point pt1 = pnt - len*norm;
	  Point pt2 = pnt + len*norm;
	  shared_ptr<SplineCurve> crv = 
	    shared_ptr<SplineCurve>(new SplineCurve(pt1, pt2));

	  // Get model intersections
	  vector<intersection_point> int_pts = getIntSfModelsCrv(models, crv);

	  // Sort intersection points according to curve parameter
	  std::sort(int_pts.begin(), int_pts.end(), par_compare);
	  size_t i1, i2, i3, i4;
	  for (kj=0; kj<int_pts.size(); ++kj)
	    {
	      // Find the next intersection with the same model
	      for (kr=kj+1; kr<int_pts.size(); ++kr)
		if (int_pts[kj].shell_idx == int_pts[kr].shell_idx)
		  break;
	      if (kr == int_pts.size())
		kr = kj;

	      found = findBoundaryShell(models[int_pts[kj].shell_idx], 
					i1, i2);
	      if (!found)
		{
		  vector<shared_ptr<SurfaceModel> > curr;
		  curr.push_back(models[int_pts[kj].shell_idx]);
		  boundary_shells_.push_back(curr);
		  i1 = boundary_shells_.size() - 1;
		}

	      // Add inner boundary shells
	      for (size_t kh=kj+1; kh<=kr; ++kh)
		{
		  found = findBoundaryShell(models[int_pts[kh].shell_idx], 
					    i3, i4);
		  if (!found)
		    boundary_shells_[i1].push_back(models[int_pts[kh].shell_idx]);
		}
	    }
	}
    }
}

//===========================================================================
bool VolumeModel::allSplines() const
//===========================================================================
{
  for (size_t ki=0; ki<bodies_.size(); ++ki)
    {
      if (!bodies_[ki]->isSpline())
	return false;
    }

  return true;
}

//===========================================================================
void VolumeModel::getAllVertices(vector<shared_ptr<Vertex> >& vertices) const
//===========================================================================
{
  // Collect all vertices
  std::set<shared_ptr<Vertex> > all_vertices;  // All vertices in the model represented once

  for (size_t ki=0; ki<bodies_.size(); ++ki)
    {
      vector<shared_ptr<Vertex> > curr_vertices = 
	bodies_[ki]->vertices();
      all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
    }

  vertices.clear();
  vertices.insert(vertices.end(), all_vertices.begin(), all_vertices.end());
}

//===========================================================================
void 
VolumeModel::getRadialEdges(vector<shared_ptr<EdgeVertex> >& rad_edges) const
//===========================================================================
{
  // Collect all radial edges
  std::set<shared_ptr<EdgeVertex> > radialedges;  // All vertices in the model represented once

  for (size_t ki=0; ki<bodies_.size(); ++ki)
    {
      vector<shared_ptr<EdgeVertex> > curr_edges = 
	bodies_[ki]->radialEdges();
      radialedges.insert(curr_edges.begin(), curr_edges.end());
    }

  rad_edges.clear();
  rad_edges.insert(rad_edges.end(), radialedges.begin(), radialedges.end());
}

//===========================================================================
void 
VolumeModel::uniqueNonRadialEdges(vector<shared_ptr<ftEdge> >& edges) const
//===========================================================================
{
  edges.clear();
  for (size_t ki=0; ki<bodies_.size(); ++ki)
    {
      vector<shared_ptr<ftEdge> > curr_edges = 
	bodies_[ki]->uniqueNonRadialEdges();
      edges.insert(edges.end(), curr_edges.begin(), curr_edges.end());
    }
}

//===========================================================================
bool VolumeModel::isCornerToCorner(double tol) const
//===========================================================================
{
//  MESSAGE("VolumeModel::isCornerToCorner. Not implemented");
  size_t nmb_bodies = bodies_.size();
  size_t ki, kj;
  shared_ptr<ftSurface> face1;
  shared_ptr<ftSurface> face2;

  for (ki=0; ki<nmb_bodies; ++ki)
    {
      for (kj=0; kj<nmb_bodies; ++kj)
	{
	  if (!bodies_[ki]->areNeighbours(bodies_[kj].get(),
					  face1, face2))
	    continue;

	  if (!bodies_[ki]->isCornerToCorner(bodies_[kj], tol))
	    return false;
	}
    }
  return true;
}

//===========================================================================
void VolumeModel::makeCornerToCorner(double tol)
//===========================================================================
{
//  MESSAGE("VolumeModel::makeCornerToCorner. Not implemented");
  VolumeAdjacency computeTop(toptol_.gap, toptol_.neighbour);
  bool changed = true;
  while (changed)
    {
      // As long as there is a new or unknown configuration, continue
      // to check for corner mismatches
      
      changed = false;  // No modifications performed yet

      size_t ki, kj, kr;
      shared_ptr<ftSurface> face1;
      shared_ptr<ftSurface> face2;

      for (ki=0; ki<bodies_.size(); ++ki)
	{
	  for (kj=0; kj<bodies_.size(); ++kj)
	    {
	      if (!bodies_[ki]->areNeighbours(bodies_[kj].get(),
					      face1, face2))
		continue;

	      if (!bodies_[ki]->isCornerToCorner(bodies_[kj], tol))
		{
		  changed = true;

		  vector<shared_ptr<ftVolume> > nbodies1;
		  vector<shared_ptr<ftVolume> > nbodies2;
		  bodies_[ki]->splitAtInternalCorner(bodies_[kj].get(), 
						     nbodies1, nbodies2,
						     tol);
		  if (nbodies1.size() > 0)
		    {
		      // The first volume has been split
		      // Update topology and volume lists
		      //computeTop.removeSolid(bodies_, bodies_[ki]);

		      int idx = getIndex(bodies_[ki]);
		      bodies_.erase(bodies_.begin() + idx);
		    }

		  if (nbodies2.size() > 0)
		    {
		      // The second volume has been split
		      // Update topology and volume lists
		      //computeTop.removeSolid(bodies_, bodies_[kj]);

		      int idx = getIndex(bodies_[kj]);
		      bodies_.erase(bodies_.begin() + idx);
		    }

		  for (kr=0; kr<nbodies1.size(); ++ kr)
		    {
		      //computeTop.addSolid(bodies_, nbodies1[kr]);
		      bodies_.push_back(nbodies1[kr]);
		    }
      
		  for (kr=0; kr<nbodies2.size(); ++ kr)
		    {
		      //computeTop.addSolid(bodies_, nbodies2[kr]);
		      bodies_.push_back(nbodies2[kr]);
		    }

		  break;
		}

	      if (changed)
		break;
	    }
	}
    }
}

//===========================================================================
void VolumeModel::makeCommonSplineSpaces()
//===========================================================================
{
  bool changed = true;
  while (changed)
    {
      // As long as there is a new or unknown configuration, continue
      // to check for spline space mismatches
      
      changed = false;  // No modifications performed yet

      size_t ki, kj;
      shared_ptr<ftSurface> face1;
      shared_ptr<ftSurface> face2;

      for (ki=0; ki<bodies_.size(); ++ki)
	{
	  for (kj=ki+1; kj<bodies_.size(); ++kj)
	    {
	      if (!bodies_[ki]->areNeighbours(bodies_[kj].get(), face1, face2))
		continue;

	      if (!bodies_[ki]->commonSplineSpace(bodies_[kj].get(), 
						  toptol_.gap))
		{
		  changed = 
		    bodies_[ki]->makeCommonSplineSpace(bodies_[kj].get());
		  break;
		}
	      if (changed)
		break;
	    }
	}
    }
}

//===========================================================================
void VolumeModel::averageCorrespondingCoefs()
//===========================================================================
{
  // First average coefficients at vertices
  vector<shared_ptr<Vertex> > vx;
  getAllVertices(vx);
  size_t ki, kj;
  for (ki=0; ki<vx.size(); ++ki)
    {
      averageVolCorner(vx[ki].get());
    }


  // Then average along radial edges
  vector<shared_ptr<EdgeVertex> > radial;
  getRadialEdges(radial);
  for (size_t ki=0; ki<radial.size(); ++ki)
    averageVolBoundaries(radial[ki].get());

  // Finally average inner coefficients
  shared_ptr<ftSurface> face1;
  shared_ptr<ftSurface> face2;
  for (ki=0; ki<bodies_.size(); ++ki)
    {
      for (kj=ki+1; kj<bodies_.size(); ++kj)
	{
	  if (!bodies_[ki]->areNeighbours(bodies_[kj].get(), face1, face2))
	    continue;

	  if (!bodies_[ki]->commonSplineSpace(bodies_[kj].get(), 
					      toptol_.gap))
	    continue;

	  
	  shared_ptr<SplineVolume> vol1 = 
	    dynamic_pointer_cast<SplineVolume, ParamVolume>(bodies_[ki]->getVolume());
	  shared_ptr<SplineVolume> vol2 = 
	    dynamic_pointer_cast<SplineVolume, ParamVolume>(bodies_[kj]->getVolume());
	  if (!(vol1.get() && vol2.get()))
	    continue;
	      
	  int dim = vol1->dimension();
	  vector<pair<int,int> > coef_corr;
	  bool found;
	  found = bodies_[ki]->getCorrCoefEnumeration(bodies_[kj].get(), 
						      toptol_.gap, 
						      coef_corr);
	  for (size_t kr=0; kr<coef_corr.size(); ++kr)
	    {
	      Point coef1(vol1->coefs_begin()+dim*coef_corr[kr].first,
			  vol1->coefs_begin()+dim*(coef_corr[kr].first+1));
	      Point coef2(vol2->coefs_begin()+dim*coef_corr[kr].second,
			  vol2->coefs_begin()+dim*(coef_corr[kr].second+1));
	      Point coef = 0.5*(coef1 + coef2);
	      vol1->replaceCoefficient(coef_corr[kr].first, coef);
	      vol2->replaceCoefficient(coef_corr[kr].second, coef);
	    }

	  (void)vol1->getBoundarySurfaces(true);
	  (void)vol2->getBoundarySurfaces(true);

	  #ifdef DEBUG_VOL2
	  std::ofstream of("av_vols.g2");
	  vol1->writeStandardHeader(of);
	  vol1->write(of);
	  vol2->writeStandardHeader(of);
	  vol2->write(of);
	  int stop_break = 1;
	  #endif
	}
    }

  for (ki=0; ki<bodies_.size(); ++ki)
    bodies_[ki]->updateBoundaryInfo();
}

//===========================================================================
int VolumeModel::nmbBoundaries() const
//===========================================================================
{
  int nmb = 0;
  for (size_t ki=0; ki<boundary_shells_.size(); ++ki)
      nmb += (int)boundary_shells_[ki].size();
  return nmb;
}

//===========================================================================
shared_ptr<SurfaceModel> VolumeModel::getOuterBoundary(int idx) const
//===========================================================================
{
  size_t ki, kj;
  int idx2;
  for (ki=0, idx2=0; ki<boundary_shells_.size(); ++ki)
    for (kj=0; kj<boundary_shells_[ki].size(); ++kj, ++idx2)
      if (idx2 == idx)
	return boundary_shells_[ki][kj];

  return shared_ptr<SurfaceModel>();
}


//===========================================================================
vector<shared_ptr<VolumeModel> >  VolumeModel::getConnectedModels()
//===========================================================================
{
    vector<shared_ptr<VolumeModel> > models;
    vector<shared_ptr<ftVolume> > curr_set;
    vector<shared_ptr<ftVolume> > all_sets;

    // Find start face for collecting a compact set
    for (size_t ki=0; ki<bodies_.size(); ++ki)
      {
	curr_set.clear();

	// Check if the body is used already
	size_t kj;
	for (kj=0; kj<all_sets.size(); ++kj)
	  if (bodies_[ki].get() == all_sets[kj].get())
	    break;

	if (kj < all_sets.size())
	  continue;  // The face already belongs to a set

	shared_ptr<ftVolume> vol = bodies_[ki];

	// Store first face and all connected faces
	getCurrConnectedModel(vol, curr_set, all_sets);

	// Make surface model
	shared_ptr<VolumeModel> curr_model = 
	  shared_ptr<VolumeModel>(new VolumeModel(curr_set,
						  toptol_.gap,
						  toptol_.neighbour,
						  toptol_.kink,
						  toptol_.bend));
	models.push_back(curr_model);
      }
    return models;
  }

  
  //===========================================================================
  void 
  VolumeModel::getCurrConnectedModel(shared_ptr<ftVolume>& vol,
				    vector<shared_ptr<ftVolume> >& curr_set,
				    vector<shared_ptr<ftVolume> >& all_sets) const
  //===========================================================================
  {
    // Store current volume
    curr_set.push_back(vol);
    all_sets.push_back(vol);

    // Fetch all neighbours
    vector<ftVolume*> neighbours;
    vol->getAdjacentBodies(neighbours);
    
    // For all neighbours, store the neighbour and all its neighbours
    // as long as they are not found already
    for (size_t ki=0; ki<neighbours.size(); ki++)
      {
	size_t kj;
	for (kj=0; kj<all_sets.size(); ++kj)
	  if (neighbours[ki] == all_sets[kj].get())
	    break;

	if (kj < all_sets.size())
	  continue;  // Vol categorized already
	
	// Handle neighbours to neighbour
	shared_ptr<ftVolume> curr_vol = fetchAsSharedPtr(neighbours[ki]);
	getCurrConnectedModel(curr_vol, curr_set, all_sets);
      }
  }


//===========================================================================
vector<shared_ptr<ftSurface> > VolumeModel::getBoundaryFaces() const
//===========================================================================
{
  // Fetch all faces lying at outer boundaries, i.e. all faces with no twin
  vector<shared_ptr<ftSurface> > bd_faces;
  for (size_t ki=0; ki<bodies_.size(); ++ki)
    {
      int nmb_shells = bodies_[ki]->nmbOfShells();
      for (int kj=0; kj<nmb_shells; ++kj)
	{
	  shared_ptr<SurfaceModel> shell = bodies_[ki]->getShell(kj);

	  // For each face, check if it has a twin
	  int nmb_faces = shell->nmbEntities();
	  for (int kr=0; kr<nmb_faces; ++kr)
	    {
	      shared_ptr<ftSurface> face = shell->getFace(kr);
	      if (!face->twin())
		bd_faces.push_back(face);
	    }
	}
    }
  return bd_faces;
}


//===========================================================================
vector<shared_ptr<ftSurface> > VolumeModel::getBoundaryFaces(int boundary_idx) const
//===========================================================================
{
  vector<shared_ptr<ftSurface> > bd_faces;

  if (boundary_idx < 0)
    return bd_faces;

  for (int i = 0; i < (int)boundary_shells_.size(); ++i)
    {
      if (boundary_idx < (int)boundary_shells_[i].size())
	{
	  shared_ptr<SurfaceModel> curr = boundary_shells_[i][boundary_idx];
	  int nmb = curr->nmbEntities();
	  for (int j = 0; j < nmb; ++j)
	    {
	      shared_ptr<ftFaceBase> face = curr->getFace(j);
	      shared_ptr<ftSurface> face2 = 
		dynamic_pointer_cast<ftSurface, ftFaceBase>(face);
	      if (face2.get())
		bd_faces.push_back(face2);  // The face should always be an ftSurface
	    }
	  return bd_faces;
	}
      boundary_idx -= (int)boundary_shells_[i].size();
    }

  return bd_faces;
}


//===========================================================================
vector<shared_ptr<ftSurface> > VolumeModel::getUniqueInnerFaces() const
//===========================================================================
{
  vector<shared_ptr<ftSurface> > faces;
  for (size_t ki=0; ki<bodies_.size(); ++ki)
    {
      int nmb_shells = bodies_[ki]->nmbOfShells();
      for (int kj=0; kj<nmb_shells; ++kj)
	{
	  shared_ptr<SurfaceModel> shell = bodies_[ki]->getShell(kj);

	  // For each face, check if it has a twin
	  int nmb_faces = shell->nmbEntities();
	  for (int kr=0; kr<nmb_faces; ++kr)
	    {
	      shared_ptr<ftSurface> curr_face = shell->getFace(kr);
	      if (curr_face->twin())
		{
		  // Check if this face or its twin is collected already
		  size_t kh;
		  for (kh=0; kh<faces.size(); ++kh)
		    if (faces[kh].get() == curr_face.get() ||
			faces[kh].get() == curr_face->twin())
		      break;

		  if (kh == faces.size())
		    faces.push_back(curr_face);
		}
	    }
	}
    }
  return faces;
}

//===========================================================================
vector<VolumeModel::intersection_point> 
 VolumeModel::getIntSfModelsCrv(vector<shared_ptr<SurfaceModel> >& models, 
				shared_ptr<SplineCurve> crv) const
//===========================================================================
{
  // Intersect all surface models with the curve
  vector<intersection_point> result;
  for (size_t ki=0; ki<models.size(); ++ki)
    {
      vector<bool> seg;
      vector<pair<ftPoint, double> > curr_result =
	models[ki]->intersect(crv, seg);
      for (size_t kj=0; kj<curr_result.size(); ++kj)
	{
	  result.push_back(intersection_point(curr_result[kj].second, (int)ki));
	}
    }
  return result;
}

//===========================================================================
 bool VolumeModel::findBoundaryShell(shared_ptr<SurfaceModel> model,
				     size_t& idx1, size_t& idx2) const
//===========================================================================
{
  for (idx1=0; idx1<boundary_shells_.size(); ++idx1)
    {
      for (idx2=0; idx2<boundary_shells_[idx1].size(); ++idx2)
	if (boundary_shells_[idx1][idx2].get() == model.get())
	  break;
      if (idx2 < boundary_shells_[idx1].size())
	break;
    }
  if (idx1 < boundary_shells_.size())
    return true;

  return false;
}
//===========================================================================
 void VolumeModel::regularizeBdShells()
//===========================================================================
{
  bool modified = true;
  bool changed = false;
  vector<pair<Point,Point> > dummy;
  vector<SurfaceModel*> modified_adjacent;
  int ki, kj;

  // Sort volumes according to an increasing number of holes in the
  // boundary surfaces
  int nmb_vols = bodies_.size();
  vector<int> perm(nmb_vols);
  for (ki=0; ki<nmb_vols; ++ki)
    perm[ki] = ki;
						
  // Count number of holes
  vector<int> nmb_holes(nmb_vols, 0);
  for (ki=0; ki<nmb_vols; ++ki)
    {
      shared_ptr<SurfaceModel> shell = bodies_[perm[ki]]->getOuterShell();
      int nmb = shell->nmbEntities();
      int curr_nmb_holes = 0;
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> face = shell->getFace(kj);
	  curr_nmb_holes += (face->nmbBoundaryLoops() - 1);
	}
      nmb_holes[ki] = curr_nmb_holes;
    }
  
  // Do the sorting
  for (ki=0; ki<nmb_vols; ++ki)
    for (kj=ki+1; kj<nmb_vols; ++kj)
      if (nmb_holes[perm[ki]] > nmb_holes[perm[kj]])
	std::swap(perm[ki], perm[kj]);

  while (modified)
    {
      // As long as one connected volumes are modified, proceed
      modified = false;

      for (size_t ki=0; ki<bodies_.size(); ++ki)
	{
	  bool mod2 = bodies_[perm[ki]]->regularizeBdShells(dummy, 
						      modified_adjacent);
	  if (mod2)
	    {
	      modified = true;
	      changed = true;
	    }
	}
    }

  if (changed)
    {
      // Update boundary information
      boundary_shells_.clear();
      setBoundarySfs();
    }
}

//===========================================================================
void VolumeModel::replaceNonRegVolumes(int degree, int split_mode)
//===========================================================================
{
  bool pattern_split = true; //false;
  int nmb_vols = nmbEntities();
  vector<SurfaceModel*> modified_ajacent;
  int ki, kj;

  // Sort volumes according to an increasing number of holes in the
  // boundary surfaces
  vector<int> perm(nmb_vols);
  for (ki=0; ki<nmb_vols; ++ki)
    perm[ki] = ki;
						
  // Count number of holes
  vector<int> nmb_holes(nmb_vols, 0);
  for (ki=0; ki<nmb_vols; ++ki)
    {
      shared_ptr<SurfaceModel> shell = bodies_[perm[ki]]->getOuterShell();
      int nmb = shell->nmbEntities();
      int curr_nmb_holes = 0;
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> face = shell->getFace(kj);
	  curr_nmb_holes += (face->nmbBoundaryLoops() - 1);
	}
      nmb_holes[ki] = curr_nmb_holes;
    }
  
  // Do the sorting
  for (ki=0; ki<nmb_vols; ++ki)
    for (kj=ki+1; kj<nmb_vols; ++kj)
      if (nmb_holes[perm[ki]] > nmb_holes[perm[kj]])
	std::swap(perm[ki], perm[kj]);
  
  // Create regular volumes
  bool changed = true;
  while (changed)
    {
      changed = false;
      for (ki=0; ki<nmb_vols; ++ki)
	{
	  if (!bodies_[perm[ki]]->isRegularized())
	    {
	      vector<shared_ptr<ftVolume> > regvols =
		bodies_[perm[ki]]->replaceWithRegVolumes(degree, modified_ajacent,
							 false, split_mode, 
							 pattern_split);
	      
	      if (regvols.size() > 0)
		{
		  bodies_.erase(bodies_.begin() + perm[ki]);
		  for (kj=ki+1; kj<nmb_vols; ++kj)
		    perm[kj] -= 1;
		  perm.erase(perm.begin() + ki);
		  ki--;
		  nmb_vols--;
		  append(regvols);
		  changed = true;
		}
	      pattern_split = true;
	    }
	}
    }

  // Update boundary information
  boundary_shells_.clear();
  setBoundarySfs();
 }

//===========================================================================
 void VolumeModel::averageVolCorner(Vertex* vx) 
//===========================================================================
{
  vector<Body*> bds = vx->getBodies();
  Point vx_point = vx->getVertexPoint();

  if (bds.size() <= 1)
    return;  // Nothing to average

  // Collect spline volumes
  vector<shared_ptr<SplineVolume> > vols;
  vector<int> idxs;
  Point pos(0.0, 0.0, 0.0);
  for (size_t ki=0; ki<bds.size(); ++ki)
    {
      ftVolume* vol1 = dynamic_cast<ftVolume*>(bds[ki]);
      if (!vol1)
	continue;

      shared_ptr<SplineVolume> vol2 =
	dynamic_pointer_cast<SplineVolume,ParamVolume>(vol1->getVolume());
      if (!vol2)
	continue;

      double upar, vpar, wpar, dist;
      Point corner;
      int idx;
      idx = vol2->closestCorner(vx_point, upar, vpar, wpar, corner, dist);
      if (dist < toptol_.neighbour)
	{
	  vols.push_back(vol2);
	  idxs.push_back(idx);
	  pos += corner;
	}
    }

#ifdef DEBUG_VOL2
  std::ofstream of1("av_vols1.g2");
  for (size_t kh=0; kh<vols.size(); ++kh)
    {
      vols[kh]->writeStandardHeader(of1);
      vols[kh]->write(of1);
    }
  int stop_break = 1;
#endif
  // Replace corner coefficients
  pos /= (double)(vols.size());
  for (size_t ki=0; ki<vols.size(); ++ki)
    vols[ki]->replaceCoefficient(idxs[ki], pos);

#ifdef DEBUG_VOL2
  std::ofstream of2("av_vols2.g2");
  for (size_t kh=0; kh<vols.size(); ++kh)
    {
      vols[kh]->writeStandardHeader(of2);
      vols[kh]->write(of2);
    }
  stop_break = 1;
#endif
}

//===========================================================================
 void VolumeModel::averageVolBoundaries(EdgeVertex* edge) 
//===========================================================================
 {
  // Get all adjacent volumes
  vector<Body*> bodies = edge->getAdjacentBodies();
  if (bodies.size() <= 1)
    return;  // Nothing to average

  vector<ftVolume*> ftvols;
  vector<shared_ptr<SplineVolume> > vols;
  vector<vector<int> > coef_enum;
  vector<bool> same_orientation;
  size_t ki;

  // Get first spline volume
  for (ki=0; ki<bodies.size(); ++ki)
    {
      ftVolume *vol1 = dynamic_cast<ftVolume*>(bodies[ki]);
      if (!vol1)
	continue;
      shared_ptr<SplineVolume> vol0 = 
	dynamic_pointer_cast<SplineVolume,ParamVolume>(vol1->getVolume());
      if (vol0.get())
	{
	  vols.push_back(vol0);
	  ftvols.push_back(vol1);
	  break;
	}
    }

  if (ki >= bodies.size())
    return;  // At most one spline volume. Nothing to average

  for (ki++; ki<bodies.size(); ++ki)
    {
      ftVolume *vol1 = dynamic_cast<ftVolume*>(bodies[ki]);
      if (!vol1)
	continue;
      shared_ptr<SplineVolume> vol0 = 
	dynamic_pointer_cast<SplineVolume,ParamVolume>(vol1->getVolume());
      if (!vol0.get())
	continue;

      // Adjacency info
      VolumeAdjacencyInfo info = ftvols[0]->getCornerAdjacencyInfo(vol1, edge,
								   toptol_.neighbour);
      if (!info.corner_adjacency_)
	continue;

      if (ftvols.size() == 1)
	{
	  vector<int> coef_idx0;
	  bool found = VolumeTools::getVolBdCoefEnumeration(vols[0], info.bd_idx_1_,
					       info.edg_idx_1_, coef_idx0);
	  if (!found)
	    return;  
	  coef_enum.push_back(coef_idx0);
	}

      vector<int> coef_idx;
      bool found = VolumeTools::getVolBdCoefEnumeration(vol0, info.bd_idx_2_,
					   info.edg_idx_2_, coef_idx);
      if (!found)
	continue;
      ftvols.push_back(vol1);
      vols.push_back(vol0);
      coef_enum.push_back(coef_idx);
      same_orientation.push_back(info.same_orient_edge_);
    }
	 
  #ifdef DEBUG_VOL2
  std::ofstream of("rad_vol.g2");
  for (size_t kj=0; kj<vols.size(); ++kj)
    {
      vols[kj]->writeStandardHeader(of);
      vols[kj]->write(of);
    }
  #endif

  if (vols.size() <= 1)
    return; // Nothing to average

  size_t nmb = coef_enum[0].size();
  int dim = vols[0]->dimension();
  for (size_t kj=0; kj<coef_enum[0].size(); ++kj)
    {
      int ix1 = coef_enum[0][kj];
      Point coef(vols[0]->coefs_begin()+ix1*dim, 
		 vols[0]->coefs_begin()+(ix1+1)*dim);
      for (ki=1; ki<vols.size(); ++ki)
	{
	  // Check 
	  if (coef_enum[ki].size() != nmb)
	    continue;

	  int idx = (same_orientation[ki-1]) ? (int)kj : (int)(nmb - kj - 1);
	  int ix2 = coef_enum[ki][idx];
	  Point tmp(vols[ki]->coefs_begin()+ix2*dim, 
		    vols[ki]->coefs_begin()+(ix2+1)*dim);

	  coef += tmp;
	}

      coef /= (double)(vols.size());
      vols[0]->replaceCoefficient(ix1, coef);

      for (ki=1; ki<vols.size(); ++ki)
	{
	  // Check 
	  if (coef_enum[ki].size() != nmb)
	    continue;

	  int idx = (same_orientation[ki-1]) ? (int)kj : (int)(nmb - kj - 1);
	  int ix2 = coef_enum[ki][idx];
	  vols[ki]->replaceCoefficient(ix2, coef);
	}
    }
      
  #ifdef DEBUG_VOL2
  std::ofstream of2("rad_vol2.g2");
  for (size_t kj=0; kj<vols.size(); ++kj)
    {
      vols[kj]->writeStandardHeader(of2);
      vols[kj]->write(of2);
    }
  int stop_break = 1;
  #endif
 }


//===========================================================================
bool VolumeModel::checkModelTopology()
//===========================================================================
 {
   bool isOK = true;
   size_t ki, kh;
   for (ki=0; ki<bodies_.size(); ++ki)
     {
       bool bodyOK = bodies_[ki]->checkBodyTopology();
       if (!bodyOK)
	 isOK = false;
     }

   for (ki=0; ki<boundary_shells_.size(); ++ki)
     {
       for (kh=0; kh<boundary_shells_[ki].size(); ++kh)
	 {
	   // Check back pointers
	   int nmb = boundary_shells_[ki][kh]->nmbEntities();
	   for (int kj=0; kj<nmb; ++kj)
	     {
	       Body *bd = boundary_shells_[ki][kh]->getFace(kj)->getBody();
	       size_t kr;
	       for (kr=0; kr<bodies_.size(); ++kr)
		 if (bd == bodies_[kr].get())
		   break;
	       if (kr == bodies_.size())
		 {
		   std::cout << "Boundary shell back pointer inconsistency, face = ";
		   std::cout << boundary_shells_[ki][kh]->getFace(kj) << std::endl;
		   isOK = false;
		 }
	     }
	 }
     }

   vector<shared_ptr<EdgeVertex> > rad;
   getRadialEdges(rad);
   for (ki=0; ki<rad.size(); ++ki)
     {
       bool radOK = rad[ki]->checkRadialEdgeTopology();
       if (!radOK)
	 isOK = false;
     }

   return isOK;
 }
