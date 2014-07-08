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
#include "GoTools/compositemodel/RegularizeUtils.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/Path.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/HermiteInterpolator.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/topology/FaceAdjacency.h"
#include "GoTools/topology/FaceConnectivityUtils.h"
#include "GoTools/geometry/GoIntersections.h"

#include <fstream>
#include <cstdlib>

//#define DEBUG_REG

using std::vector;
using std::set;
using std::make_pair;
using std::pair;


namespace Go {

//==========================================================================
RegularizeFace::RegularizeFace(shared_ptr<ftSurface> face, 
			       double epsge, double angtol, 
			       double tol2, bool split_in_cand)
//==========================================================================
  : epsge_(epsge), angtol_(angtol), tol2_(tol2), bend_(5.0*angtol), face_(face),
    split_in_cand_(split_in_cand), split_mode_(1), divideInT_(2), 
    top_level_(true), isolate_fac_(0.6)
{
  vector<shared_ptr<ftSurface> > faces(1);
  faces[0] = face;
  model_ = shared_ptr<SurfaceModel>(new SurfaceModel(epsge, epsge, tol2,
						     angtol, 5.0*angtol,
						     faces));
}

//==========================================================================
RegularizeFace::RegularizeFace(shared_ptr<ftSurface> face, 
			       double epsge, double angtol, 
			       double tol2, double bend, 
			       bool split_in_cand)
//==========================================================================
  : epsge_(epsge), angtol_(angtol), tol2_(tol2), bend_(bend), face_(face),
    split_in_cand_(split_in_cand), split_mode_(1), divideInT_(2), 
    top_level_(true), isolate_fac_(0.6)
{
  vector<shared_ptr<ftSurface> > faces(1);
  faces[0] = face;
  model_ = shared_ptr<SurfaceModel>(new SurfaceModel(epsge, epsge, tol2,
						     angtol, bend,
						     faces));
}
//==========================================================================
RegularizeFace::RegularizeFace(shared_ptr<ftSurface> face, 
			       shared_ptr<SurfaceModel> model,
			       bool split_in_cand)
//==========================================================================
  : face_(face), split_in_cand_(split_in_cand), split_mode_(1), 
    divideInT_(1), top_level_(true), isolate_fac_(0.6)
{
  epsge_ = model->getTolerances().gap;
  tol2_ = model->getTolerances().neighbour;
  angtol_ = model->getTolerances().kink;
  bend_ = model->getTolerances().bend;
  if (model->getIndex(face) >= 0)
    model_ = model;
  else
    {
      vector<shared_ptr<ftSurface> > faces(1);
      faces[0] = face;
      model_ = shared_ptr<SurfaceModel>(new SurfaceModel(epsge_, epsge_, tol2_,
							 angtol_, bend_,
							 faces));
    }
}

///==========================================================================
RegularizeFace::~RegularizeFace()
//==========================================================================
{

}

//==========================================================================
void RegularizeFace::setAxis(Point& centre, Point& axis)
//==========================================================================
{
  // Should we also store/check whether the axis information is 
  // already set?
  centre_ = centre;
  axis_ = axis;
}


//==========================================================================
void RegularizeFace::unsetAxis()
//==========================================================================
{
  centre_.resize(0);
  axis_.resize(0);
}

//==========================================================================
vector<shared_ptr<ftSurface> > RegularizeFace::getRegularFaces()
//==========================================================================
{
#ifdef DEBUG_REG
  std::ofstream of("reg_face.g2");
  face_->surface()->writeStandardHeader(of);
  face_->surface()->write(of);

  // INFO
  vector<shared_ptr<Vertex> > non_corner = face_->getNonCornerVertices(bend_);
  std::ofstream ofvx("cand_vx.g2");
  ofvx << "400 1 0 4 155 100 0 255" << std::endl;
  ofvx << non_corner.size() << std::endl;
  for (size_t kh=0; kh<non_corner.size(); ++kh)
    {
      Point vx_pt = non_corner[kh]->getVertexPoint();
      ofvx << vx_pt << std::endl;
    }
#endif

  divide();
  return sub_faces_;
}

//==========================================================================
  void RegularizeFace::divide()
//==========================================================================
{
  // Check whether a division is required. If the face has no inner trimming
  // loops and 4 corners it is not.
  // An exception to this is if the shape of the face is very awkward, but
  // this is not handled currently
  corners_ = face_->getCornerVertices(bend_);

#ifdef DEBUG_REG
  std::ofstream of("corners.g2");
  of << "400 1 0 4 155 100 0 255" << std::endl;
  of << corners_.size() << std::endl;
  for (size_t kh=0; kh<corners_.size(); ++kh)
    {
      Point vx_pt = corners_[kh]->getVertexPoint();
      of << vx_pt << std::endl;
    }
#endif

  int nmb_loops = face_->nmbBoundaryLoops();
  if (nmb_loops <= 1 && corners_.size() <= 4 &&
      !(cand_split_.size() > 1 && split_in_cand_))
    {
      sub_faces_.push_back(face_);
      // return;  // Less than 4 corners may be feasible
    }
  else
   {
     // All vertices
     vx_ = face_->vertices();

     // Identify incomplete holes along the outer boundary of the face
     vector<vector<ftEdge*> > half_holes;
     half_holes = getHalfHoles();
     int nmb_half_hole = (int)half_holes.size();

     bool is_split = false;
     if (cand_split_.size() > 1 && split_in_cand_ /*corners_.size() <= 4*/)
       {
	 // Extract extra internal loops intended to match patterns first
	 is_split = splitWithPatternLoop();
       }

     if (!is_split)
       {
	 // Check if any holes must be separated
	 if (nmb_half_hole + nmb_loops > 2)
	   {
	     // Separate holes
	     faceWithHoles(half_holes);

	   }

	 else if (nmb_half_hole + nmb_loops == 2)
	   {
	     // Divide according to the found hole
	     faceOneHole(half_holes);
	   }
	 else if (corners_.size() > 4)
	   {
	     // An outer boundary a different number of corners than 4 is found
	     // Divide accordingly
	     faceOuterBd(half_holes);
	   }
       }
   }

  // Split in T-joint like situations that ruins the corner-to-corner
  // configuaration with 4-sided surfaces that we want
#ifdef DEBUG_REG
  std::ofstream of1("pre_T_split.g2");
  for (size_t kr=0; kr<sub_faces_.size(); ++kr)
    {
      shared_ptr<ParamSurface> sf = sub_faces_[kr]->surface();
      sf->writeStandardHeader(of1);
      sf->write(of1);
    }
#endif

  std::set<shared_ptr<Vertex> > all_vx1;  // All vertices in the model re  
  for (size_t kr=0; kr<sub_faces_.size(); ++kr)
    {
      vector<shared_ptr<Vertex> > curr_vertices = 
	sub_faces_[kr]->vertices();
      all_vx1.insert(curr_vertices.begin(), curr_vertices.end());
    }
  vector<shared_ptr<Vertex> > pre_vx;
  pre_vx.insert(pre_vx.end(), all_vx1.begin(), all_vx1.end());
      
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

  if (axis_.dimension() > 0 && centre_.dimension() == 0)
    unsetAxis();  // Do not need this information any more. It can
  // mess up the division
  if (top_level_)
    {
      unsetAxis();  // Can do harm in this context
      splitInTJoints();

      std::set<shared_ptr<Vertex> > all_vx2;  // All vertices in the model re  
      for (size_t kr=0; kr<sub_faces_.size(); ++kr)
	{
	  vector<shared_ptr<Vertex> > curr_vertices = 
	    sub_faces_[kr]->vertices();
	  all_vx2.insert(curr_vertices.begin(), curr_vertices.end());
	}
      vector<shared_ptr<Vertex> > post_vx;
      post_vx.insert(post_vx.end(), all_vx1.begin(), all_vx1.end());

#ifdef DEBUG_REG
      std::ofstream of2("post_T_split.g2");
      for (size_t kr=0; kr<sub_faces_.size(); ++kr)
	{
	  shared_ptr<ParamSurface> sf = sub_faces_[kr]->surface();
	  sf->writeStandardHeader(of2);
	  sf->write(of2);
	}
#endif

      int break_val = 1;
   }
}

//==========================================================================
void RegularizeFace::splitInTJoints()
//==========================================================================
{
#ifdef DEBUG_REG
  std::ofstream of("T_face.g2");
  for (size_t kr=0; kr<sub_faces_.size(); ++kr)
    {
      shared_ptr<ParamSurface> sf = sub_faces_[kr]->surface();
      sf->writeStandardHeader(of);
      sf->write(of);
    }
#endif

  bool changed = true;
  while (changed)
    {
      changed = false;
      for (size_t ki=0; ki<sub_faces_.size(); ++ki)
	{
	  shared_ptr<ftSurface> curr = sub_faces_[ki];

	  // Find potensial T-joints
	  if (curr->nmbBoundaryLoops() == 0)
	    continue;  // Can't split

	  // Get potensial T-joints
	  vector<shared_ptr<Vertex> > Tvx = 
	    curr->getNonCornerVertices(bend_, 0);

#ifdef DEBUG_REG
	  std::ofstream of2("curr_T_face.g2");
	  curr->surface()->writeStandardHeader(of2);
	  curr->surface()->write(of2);
	  of2 << "400 1 0 4 100 155 0 255" << std::endl;
	  of2 << Tvx.size() << std::endl;
	  for (size_t ka=0; ka<Tvx.size(); ++ka)
	    of2 << Tvx[ka]->getVertexPoint() << std::endl;
#endif
	  // Remove the insignificant ones and the ones representing
	  // a seam
	  if (divideInT_ != 2)
	    removeInsignificantVertices(Tvx, false, curr.get());
#ifdef DEBUG_REG
	  of2 << "400 1 0 4 155 0 100 255" << std::endl;
	  of2 << Tvx.size() << std::endl;
	  for (size_t ka=0; ka<Tvx.size(); ++ka)
	    of2 << Tvx[ka]->getVertexPoint() << std::endl;
#endif	  
	  // Check if the vertex really indicates a T-joint
	  // Get corner vertices
	  //size_t kj;
	  vector<shared_ptr<Vertex> > corner = 
	    curr->getCornerVertices(bend_, 0);

	  if (corner.size() + Tvx.size() > 4 && divideInT_)
	    {
	      // Number of corners and T-joints exceeds 4. Not
	      // appropriate for a 4-sided surface divide
	      vector<shared_ptr<ftSurface> > faces = 
		divideInTjoint(curr, Tvx, corner);
	      if (faces.size() > 0 &&
		  !(faces.size() == 1 && faces[0].get() == face_.get()))
		{
		  model_->removeFace(face_);
		  sub_faces_.erase(sub_faces_.begin()+ki);
		  model_->append(faces, false, false);
		  sub_faces_.insert(sub_faces_.end(), faces.begin(), faces.end());
		  changed = true;
		}
	    }
	  if (changed)
	    break;
	}
    }
}

//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeFace::divideInTjoint(shared_ptr<ftSurface> face,
			       vector<shared_ptr<Vertex> >& Tvx,
			       vector<shared_ptr<Vertex> >& corner)
//==========================================================================
{
  vector<shared_ptr<ftSurface> > faces;
  if (Tvx.size() + corner.size() <= 4)
    return faces;  // No point in dividing

  // Select the T-joint to divide in
  shared_ptr<Vertex> curr;
  double min_dist = 0.0;
  size_t ki, curr_idx;
  for (ki=0; ki<Tvx.size(); ++ki)
    {
      // Compute distance between T vertex and corners on both sides
      vector<ftEdge*> edges = Tvx[ki]->getFaceEdges(face.get());
      if (edges.size() != 2)
	continue;  // Does not make sense

      if (Tvx[ki]->nmbUniqueEdges() == 2)
	continue;  // Not a T-joint

      // Check if the T-joint is allowed for split
      vector<ftSurface*> vx_faces = Tvx[ki]->faces();
      size_t kj, kr;
      for (kj=0; kj<vx_faces.size(); ++kj)
	{
	  for (kr=0; kr<nonTjoint_faces_.size(); ++kr)
	    if (vx_faces[kj] == nonTjoint_faces_[kr].get())
	      break;
	  if (kr < nonTjoint_faces_.size())
	    break;
	}
      if (kj < vx_faces.size())
	continue;  // Not a feasible T-joint
	
      double len = 0.0;
      for (kj=0; kj<edges.size(); ++kj)
	{
	  ftEdge* curr_edge = edges[kj];
	  shared_ptr<Vertex> curr_vx = Tvx[ki];
	  double tp = edges[kj]->parAtVertex(Tvx[ki].get());
	  bool next = (tp - edges[kj]->tMin() < edges[kj]->tMax() - tp);
	  shared_ptr<Vertex> other;
	  while (true)
	    {
	      other = curr_edge->getOtherVertex(curr_vx.get());
	      vector<shared_ptr<Vertex> >::iterator vxp =
		std::find(corner.begin(), corner.end(), other);
	      if (vxp != corner.end())
		break;

	      curr_vx = other;
	      curr_edge = (next) ? curr_edge->next()->geomEdge() :
		curr_edge->prev()->geomEdge();
	    }
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
    return faces;  // Cannot divide

#ifdef DEBUG_REG
  std::ofstream of("curr_T_vx.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << curr->getVertexPoint() << std::endl;
#endif

  // A T-joint vertex is found. Set current face
  face_ = face;
  
  // Get candidate vertices
  vector<shared_ptr<Vertex> > vx = face->vertices();
  vector<shared_ptr<Vertex> > cand_vx;
  ftEdge* cand_edge = 0;
  selectCandidateSplit(curr, vx, cand_vx, cand_edge, false);
 
  // Fetch info about vertices not belonging to the corners of the
  // initial face
  vector<shared_ptr<Vertex> > non_corner = face_->getNonCornerVertices(bend_);
  //removeInsignificantVertices(non_corner);

  // The remaining T-joint vertices are highly prioritized vertices
  // for splitting. Remove the current vertex
  Tvx.erase(Tvx.begin()+curr_idx);
  /*for (ki=0; ki<Tvx.size(); ++ki)
    {
      int kh;
      for (kh=(int)cand_vx.size()-1; kh>=0; --kh)
      if (Tvx[ki].get() == cand_vx[kh].get())
	{
	  // Vertex exists in both samples
	  cand_vx.erase(cand_vx.begin()+kh);
	}
	}*/

  // Divide
  //cand_edge = 0;
  vector<shared_ptr<Vertex> > dummy;
  faces = RegularizeUtils::divideVertex(face_, curr, cand_vx, cand_edge,
					/*Tvx*/  dummy, epsge_, tol2_, angtol_, 
					bend_, non_corner, centre_, axis_);
  return faces;
	    
}

//==========================================================================
void RegularizeFace::faceWithHoles(vector<vector<ftEdge*> >& half_holes)
//==========================================================================
{
  // For each hole and half hole, estimate a centre, an axis and a radius
  int nmb_loops = face_->nmbBoundaryLoops();
  vector<hole_info> holes(half_holes.size() + nmb_loops - 1);
  Point mid, axis;
  double rad, ang;
  double max_rad = 0.0;
  int ki;
  for (ki=1; ki<nmb_loops; ++ki)
    {
      shared_ptr<Loop> loop = face_->getBoundaryLoop(ki);  // Inner loop
      size_t nmb_edges = loop->size();
      vector<ftEdge*> edges(nmb_edges);
      for (size_t kj=0; kj<nmb_edges; ++kj)
	edges[kj] = loop->getEdge(kj)->geomEdge();

      bool done = Path::estimateHoleInfo(edges, mid, axis, rad, ang);
      if (done)
	{
	  holes[ki-1].setInfo(mid, axis, rad);
	  max_rad = std::max(max_rad, rad);
	}
      else
	{
	  holes.erase(holes.begin()+ki-1);
	  nmb_loops--;
	  ki--;
	}
    }

  double half_hole_fac = 1.6; //1.5; //1.75;
  for (ki=0; ki<(int)half_holes.size(); ++ki)
    {
      bool done = Path::estimateHoleInfo(half_holes[ki], mid, axis, rad, ang);
      if (done && rad < half_hole_fac*max_rad && split_mode_ == 1)
	{
	  // Large half holes are better treaded as a part of the outer
	  // boundary
	  holes[nmb_loops+ki-1].setInfo(mid, axis, rad);
	}
      else
	{
	  holes.erase(holes.begin()+nmb_loops+ki-1);
	  half_holes.erase(half_holes.begin()+ki);
	  ki--;
	}
    }
  for (ki=1; ki<(int)holes.size(); ++ki)
    if (holes[ki-1].hole_axis_*holes[ki].hole_axis_ < 0.0)
      holes[ki].hole_axis_ *= -1;

  // Check if there is any holes registered
  if (holes.size() == 0)
    {
      // Consider the outer boundary of this face
      faceOuterBd(half_holes);
      return;
    }

  // Check if there is more one hole to consider
  if (holes.size() == 1)
    {
      faceOneHole(half_holes);
      return;
    }

  // Compute the weight point and a cone surronding the hole axes
  // corresponding to all holes
  Point wgt_pt(0.0, 0.0, 0.0);
  DirectionCone axis_cone(holes[0].hole_axis_);
  for (ki=0; ki<(int)holes.size(); ++ki)
    {
      wgt_pt += holes[ki].hole_centre_;
      axis_cone.addUnionWith(holes[ki].hole_axis_);
    }
  wgt_pt /= (double)(holes.size());

  // Check the location of this weight point,
  // whether it is outside the outer trimming loop, inside an
  // inner trimming loop or inside the trimmed surface
  // First project the point onto the underlying surface
  double upar, vpar, dist;
  Point wgt_proj;
  shared_ptr<ParamSurface> surf = face_->surface();
  shared_ptr<BoundedSurface> bd_surf = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
  shared_ptr<ParamSurface> surf2;
  if (bd_surf.get())
    surf2 = bd_surf->underlyingSurface();
  else
    surf2 = surf;
  surf2->closestPoint(wgt_pt, upar, vpar, wgt_proj, dist, epsge_);
  Point wgt_par(upar, vpar);
  Point wgt_norm;
  surf2->normal(wgt_norm, upar, vpar);

  int position;
  if (bd_surf.get())
    {
      // Check if the closest point is internal to the surface
      CurveBoundedDomain dom = bd_surf->parameterDomain();
      Vector2D param(upar,vpar);
      if (!dom.isOnBoundary(param, epsge_))
	dist = 0.0;
    }
  if (dist > epsge_)
    position = -1;
  else
    position = positionWeigthPoint(wgt_par);
   // Check if the holes can be sorted along a line
  // In that case divide according to this chord
  Point pnt, dir;
  vector<double> parval;
  vector<int> perm;
  vector<shared_ptr<ftSurface> > faces;
  bool sorted = false;
//  if (position >= 0)
    sorted = sortAlongLine(holes, pnt, dir, parval, perm);
  if (sorted && position >= 0)
    {
#ifdef DEBUG_REG
      std::ofstream hole_out("hole_out.g2");
      for (size_t k4=0; k4<holes.size(); ++k4)
	{
	  hole_out << "400 1 0 4 255 0 0 255" << std::endl;
	  hole_out << " 1 " << std::endl;
	  hole_out << holes[perm[k4]].hole_centre_ << std::endl;
	}
	  
#endif
    faces = divideAcrossLine(half_holes, holes, pnt, dir, perm);
    }
  
  if (faces.size() == 0 /*&& position >= 0*/)
    {
      // The weight point lies inside a hole
      // Divide according to this hole
      // First remove hole with weight point from the list
      Point centre, axis;
      if (position > 0)
	{
	  centre = holes[position-1].hole_centre_;
	  axis = holes[position-1].hole_axis_;
	  holes.erase(holes.begin() + position - 1);

	  // Store axis information
	  setAxis(centre, axis);
	}
      else
	{
	  centre = wgt_pt;
	  Point tmp_axis(0.0, 0.0, 0.0);
	  for (size_t kr=0; kr<holes.size(); ++kr)
	    tmp_axis += holes[kr].hole_axis_;
	  tmp_axis /= (double)(holes.size());
	  axis = tmp_axis;
	}
      vector<double> angles;
      double ang_limit = 0.25*M_PI;
      if (!sorted)
	{
	  if (axis_cone.greaterThanPi() || axis_cone.angle() > ang_limit)
	    sorted = false;
	  else
	    sorted = sortRadially(holes, wgt_pt, axis, angles, perm);
	}
      if (sorted)
	{
	  // Split between surrounding holes with respect to the
	  // mid hole
	  if (position > 0 && holes.size() == 1)
	    faces = isolateOneHoleRadially(centre, axis, holes[0]);
	  else /*if (position > 0)*/
	    faces = isolateHolesRadially(half_holes,
					 centre,
					 axis,
					 position - 1,
					 holes, perm);
	  // else
	  //   faces = isolateHolesRadially2(half_holes,
	  // 				  centre,
	  // 				  axis,
	  // 				  position - 1,
	  // 				  holes, perm);
	    
	}
    }      

  if (faces.size() == 0 && position == -1)
    {
      faces = faceOuterBdFaces(half_holes);
    }

  if (faces.size() == 0 && holes.size() > 0)
    {
      // Check if the holes may be sorted along a line in the parameter domain
      identifyParLine(holes, dir);
      faces = isolateHolesParametrically(half_holes, dir);
    }
      
  if (faces.size() == 0)
    {
      if (position == -1)
	{
	  // The weight point lies outside the face. Split according to outer
	  // boundary
	  faceOuterBd(half_holes);
	  return;
	}

      // Sort the holes along a chord or by angle around the weight point
      // Divide according the weight point
      vector<double> angles;
      sorted = sortRadially(holes, wgt_pt, wgt_norm, angles, perm);

      if (sorted)
	{
	  // Select the two holes with minimum distance
	  double min_dist = 1.0e6;
	  int idx1 = -1, idx2 = -1;
	  for (size_t kj=0; kj<holes.size(); ++kj)
	    for (size_t kr=kj+1; kr<holes.size(); ++kr)
	      {
		double dist = 
		  holes[perm[kj]].hole_centre_.dist(holes[perm[kr]].hole_centre_);
		if (dist < min_dist)
		  {
		    idx1 = perm[kj];
		    idx2 = perm[kr];
		    min_dist = dist;
		  }
	      }

	  if (idx1 < 0)
	    {
	      // Don't know how to handle this situation
	      // Write problem surface and quit
#ifdef DEBUG_REG
	      std::ofstream of("unhandled.g2");
	      face_->surface()->writeStandardHeader(of);
	      face_->surface()->write(of);
#endif
	      THROW("Cannot sort holes of face");
	    }

	  vector<hole_info> holes2(2);
	  holes2[0] = holes[idx1];
	  holes2[1] = holes[idx2];
	  vector<int> perm2(2);
	  perm2[0] = 0;
	  perm2[1] = 1;
	  Point pnt2 = 0.5*(holes2[0].hole_centre_ + holes2[1].hole_centre_);
	  Point dir2 = holes2[1].hole_centre_ - holes2[0].hole_centre_;
	  dir2.normalize();
	  faces = divideAcrossLine(half_holes, holes2, pnt2, dir2, perm2);
	}
      else
	{
	  // Don't know how to handle this situation
	  // Write problem surface and quit
#ifdef DEBUG_REG
	  std::ofstream of("unhandled.g2");
	  face_->surface()->writeStandardHeader(of);
	  face_->surface()->write(of);
#endif
	  THROW("Cannot sort holes of face");
	}
    }

  // Update topology and treat each sub face
  int nmb_faces = (int)faces.size();
  if (nmb_faces > 1)
    {
      model_->removeFace(face_);
      model_->append(faces, false, false);
      for (int kj=0; kj<nmb_faces; )
	{
	  RegularizeFace regularize(faces[kj], model_);
	  if (axis_.dimension() > 0)
	    regularize.setAxis(centre_, axis_);
	  regularize.setDivideInT(divideInT_);
	  regularize.setCandSplit(cand_split_);
	  regularize.setSplitMode(split_mode_);
	  regularize.unsetTopLevel();
	  vector<shared_ptr<ftSurface> > faces2 = 
	    regularize.getRegularFaces();

	  if (faces2.size() > 1)
	    {
	      faces.erase(faces.begin()+kj);
	      nmb_faces--;
	      faces.insert(faces.end(), faces2.begin(), faces2.end());

	      // Check if any new faces may be joined across the seam
	      mergeSeams(faces, nmb_faces, faces2);
	    }
	  else
	    kj++;
	}
    }
  sub_faces_.insert(sub_faces_.end(), faces.begin(), faces.end());
}

//==========================================================================
void RegularizeFace::faceOneHole(vector<vector<ftEdge*> >& half_holes)
//==========================================================================
{
  int nmb_loops = face_->nmbBoundaryLoops();
  if (nmb_loops + half_holes.size() != 2)
    return;  // Wrong function

  // Make an estimate of the size of the face (largest distance
  // between vertices)
  double face_size = 0.0;
  for (size_t ki=0; ki<vx_.size(); ++ki)
    for (size_t kj=ki+1; kj<vx_.size(); ++kj)
      face_size = std::max(face_size, vx_[ki]->getDist(vx_[kj]));

  // Check if the hole must be isolated
  vector<shared_ptr<ftSurface> > faces = initIsolateHole(half_holes);
  if (faces.size() == 0)
    {
      // Estimate mid point, axis and radius of hole/half hole
      vector<ftEdge*> edges;
      if (half_holes.size() == 1)
	edges = half_holes[0];
      else
	{
	  shared_ptr<Loop> loop = face_->getBoundaryLoop(1);  // Inner loop
	  size_t nmb_edges = loop->size();
	  edges.resize(nmb_edges);
	  for (size_t ki=0; ki<nmb_edges; ++ki)
	    edges[ki] = loop->getEdge(ki)->geomEdge();
	}
      double halfhole_ang;
      bool done = Path::estimateHoleInfo(edges, centre_, axis_, radius_, halfhole_ang);
      //double size_fac = 0.25; //0.1;
      if (!done || (half_holes.size() > 0 && split_mode_ > 1 && halfhole_ang <= M_PI) 
	  /*|| radius_ > size_fac*face_size*/)
	{
	  // Missing hole information. Divide according to the outer boundary
	 vector<vector<ftEdge*> > dummy;
	 faceOuterBd(dummy);
	 return;
	}

      // Make parameter box around hole
      vector<Point> parpt(edges.size());
      for (size_t ki = 0; ki<edges.size(); ++ki)
	{
	  ftEdge *curr = edges[ki]->geomEdge();
	  parpt[ki] = curr->faceParameter(curr->tMin());
	}
      BoundingBox parbox;
      parbox.setFromPoints(parpt);
      Point low = parbox.low();
      Point high = parbox.high();

      // Compute the points in the face corresponding to the corners of
      // the parameter box
      vector<Point> holebox(4);
      holebox[0] = face_->point(low[0],low[1]);
      holebox[1] = face_->point(high[0],low[1]);
      holebox[2] = face_->point(low[0],high[1]);
      holebox[3] = face_->point(high[0],high[1]);
  
      // For each corner in the outer boundary loop, make splitting curves
      // connecting to the hole. The curves would intersect the mid point
      // in their extension.
      vector<shared_ptr<Vertex> > corner = face_->getCornerVertices(bend_, 0);
      vector<shared_ptr<Vertex> > hole_vx;
      if (half_holes.size() == 1)
	{
	  // Fetch vertices belonging to the half hole
	  // Remove corner vertices belonging to the half hole and those
	  // next to it
	  shared_ptr<Vertex> curr_vx;
	  vector<shared_ptr<Vertex> >::iterator vx;
	  for (size_t ki=0; ki<half_holes[0].size()-1; ++ki)
	    {
	      curr_vx = half_holes[0][ki]->getVertex(false);
	      hole_vx.push_back(curr_vx);
	      vx = std::find(corner.begin(), corner.end(), curr_vx);
	      if (vx != corner.end())
		corner.erase(vx);
	    }

	  curr_vx = half_holes[0][0]->getVertex(true);
	  vx = std::find(corner.begin(), corner.end(), curr_vx);
	  if (vx != corner.end())
	    {
	      if (vx != corner.begin())
		{
		  vx--;
		  corner.erase(vx, vx+2);
		}
	      else
		{
		  corner.erase(vx);
		  corner.erase(corner.end()-1);
		}
	    }

	  curr_vx = half_holes[0][half_holes[0].size()-1]->getVertex(false);
	  vx = std::find(corner.begin(), corner.end(), curr_vx);
	  if (vx != corner.end())
	    {
	      if (vx+1 != corner.end())
		corner.erase(vx, vx+2);
	      else
		{
		  corner.erase(vx);
		  corner.erase(corner.begin());
		}
	    }
	}
      else
	// Fetch vertices belonging to the hole
	hole_vx = face_->getBoundaryLoop(1)->getVertices();

      // Remove insignificant vertices, but keep the information to avoid
      // an intersection very close to an existing vertex
      vector<shared_ptr<Vertex> > hole_vx2(hole_vx.begin(), hole_vx.end());
      removeInsignificantVertices(hole_vx);
      for (size_t kr=0; kr<hole_vx.size(); ++kr)
	{
	  vector<shared_ptr<Vertex> >::iterator vxp =
	    std::find(hole_vx2.begin(), hole_vx2.end(), hole_vx[kr]);
	  if (vxp != hole_vx2.end())
	    hole_vx2.erase(vxp);
	}

      // Check if any corners are found in the outer loop
      if (corner.size() == 0)
	{
	  // No corner. Use non-corner vertices in the split
	  corner = face_->getBoundaryLoop(0)->getVertices();
	  removeInsignificantVertices(corner);
	}

      bool OKseg = true;

      // Number of corners belonging to the hole
      vector<shared_ptr<Vertex> > hole_corners = face_->getCornerVertices(bend_, 1);

      // Check if the same underlying surface lies on both sides of an edge in the 
      // inner loop
      size_t kh;
      for (kh=0; kh<edges.size(); ++kh)
	{
	  shared_ptr<ParamSurface> surf1 = face_->surface();
	if (surf1->instanceType() == Class_BoundedSurface)
	  {
	    shared_ptr<BoundedSurface> bdsf =
	      dynamic_pointer_cast<BoundedSurface, GeomObject>(surf1);
	    surf1 = bdsf->underlyingSurface();
	  }
	  shared_ptr<ParamSurface> surf2;
	  if (edges[kh]->twin())
	    surf2 = edges[kh]->twin()->geomEdge()->face()->surface();
	  if (surf2.get() && surf2->instanceType() == Class_BoundedSurface)
	      {
		shared_ptr<BoundedSurface> bdsf =
		  dynamic_pointer_cast<BoundedSurface, GeomObject>(surf2);
		surf2 = bdsf->underlyingSurface();
	      }

	  if (surf2.get() && surf1.get() == surf2.get())
	    break;
	}

      if ((corner.size() == 0 && half_holes.size() == 0) ||
	  (kh < edges.size() && hole_corners.size() >= corner.size()))
	{
	  // Two loops. No significant vertices in the outer loop
	  faceOneHole2();
	}
      else
	{
	  // Check if the sequence of chords from the hole midpoint to
	  // the corners have increasing angles
	  size_t ki;
	  double min_corner_dist = MAXDOUBLE;
	  double max_corner_dist = 0.0;
	  if (corner.size() > 2)
	    {
	      double ang_tol = 0.1*M_PI;
	      Point vec1 = corner[0]->getVertexPoint() - centre_;
	      Point vec2 = corner[1]->getVertexPoint() - centre_;
	      min_corner_dist = std::min(vec1.length(), vec2.length());
	      max_corner_dist = std::max(vec1.length(), vec2.length());
	      double ang1 = vec1.angle(vec2);
	      for (ki=2; ki<corner.size(); ++ki)
		{
		  Point vec3 = corner[ki]->getVertexPoint() - centre_;
		  min_corner_dist = std::min(min_corner_dist, vec3.length());
		  max_corner_dist = std::max(max_corner_dist, vec3.length());
		  double ang2 = vec1.angle(vec3);
		  double ang3 = vec2.angle(vec3);
		  if (ang1 + ang3 > M_PI)
		    ang2 = 2*M_PI - ang2;
		  if (ang2 < ang1+ang_tol)
		    {
		      OKseg = false;
		      break;
		    }
		  vec1 = vec2;
		  ang1 = vec1.angle(vec3);
		}
	    }

	  if (!OKseg && min_corner_dist < MAXDOUBLE)
	    {
	      // Check if a subset of the corner points should be considered
	      vector<shared_ptr<Vertex> > sub_corner;
	      double frac = 0.1;
	      for (ki=0; ki<corner.size(); ++ki)
		{
		  double dist = corner[ki]->getVertexPoint().dist(centre_);
		  if (dist < min_corner_dist + frac*(max_corner_dist - min_corner_dist))
		    sub_corner.push_back(corner[ki]);
		}
	      if (sub_corner.size() > 2)
		{
		  corner = sub_corner;
		  OKseg = true;
		}
	    }

	  vector<shared_ptr<CurveOnSurface> > segments;
	  Point segment_point;
	  shared_ptr<BoundedSurface> bd_sf;
	  shared_ptr<ParamSurface> surf = face_->surface();
	  double len_fac = 6.0;
	  double max_edge_len = 2.0*len_fac*radius_; //4.0*radius_; //12.0*radius_; // 8.0*radius_; // 5.0*radius_;
	  for (ki=0; ki<corner.size(); ++ki)
	    {
	      Point pnt = corner[ki]->getVertexPoint();
	      Point vec = pnt - centre_;
	      vec.normalize();

	      // Compute closest point on hole boundary to the given vertex
	      int close_idx;
	      double close_par, close_dist;
	      Point close_pnt;
	      Path::closestPoint(edges, pnt, close_idx, close_par,
				 close_pnt, close_dist);

	      // Compute splitting curve
	      shared_ptr<CurveOnSurface> trim_segment = 
		computeCornerSplit(corner[ki], hole_vx, hole_vx2, bd_sf,
				   edges[close_idx]->faceParameter(close_par));
	      if (!trim_segment.get())
		OKseg = false;

	      // Check the length of the splitting curves compared to the hole
	      // radius
	      double len = 0.0;
	      double min_dist = 100.0*radius_;
	      if (trim_segment.get())
		{
		  len = trim_segment->estimatedCurveLength();
		  Point tmp1 = 
		    trim_segment->ParamCurve::point(trim_segment->startparam());
		  Point tmp2 = 
		    trim_segment->ParamCurve::point(trim_segment->endparam());
		  min_dist = std::min(pnt.dist(tmp1), pnt.dist(tmp2));
		}
	      else
		break;  // Could not compute curve. This problem is not
	      // likely to disappear. Try to reduse the faze size

	      // If the length is small
	      if (len < max_edge_len && !(min_dist > 0.1*radius_))
		segments.push_back(trim_segment);
	      else if (trim_segment.get())
		{
		  // Select point on trim segment a reasonable distance
		  // from the hole
		  double edge_len = 1.2*radius_;
		  Point p0 = pnt;
		  Point p0tan = vec;
		  vector<Point> p1(2); 
		  trim_segment->point(p1, trim_segment->startparam(), 1);
		  vector<Point> p2(2);
		  trim_segment->point(p2, trim_segment->endparam(), 1);
		  if (p1[0].dist(centre_) < p0.dist(centre_))
		    {
		      p0 = p1[0];
		      p0tan = p1[1];
		    }
		  if (p2[0].dist(centre_) < p0.dist(centre_))
		    {
		      p0 = p2[0];
		      p0tan = p2[1];
		    }

		  // Find initial point
		  //Point dir = pnt - p0;
		  Point dir = p0 - centre_;
		  if (dir*p0tan < 0.0)
		    p0tan *= -1.0;
		  p0tan.normalize();
		  Point seg_pnt;
		  if (p0.dist(centre_) < 3.0*radius_)
		    seg_pnt = p0 + edge_len*p0tan;
		  else
		    {
		      // Select box corner closest to the current corner
		      Point b0 = holebox[0];
		      double d0 = b0.dist(pnt);
		      for (int kb=1; kb<4; kb++)
			{
			  double d1 = holebox[kb].dist(pnt);
			  if (d1 < d0)
			    {
			      d0 = d1;
			      b0 = holebox[kb];
			    }
			}
		      dir = b0 - centre_;
		      dir.normalize();
		      seg_pnt = b0 + edge_len*dir;
		    }

		  double init_dist = seg_pnt.dist(centre_);
		  segment_point = seg_pnt;

		  // Find closest point on segment to initial point
		  double par, dist;
		  Point seg_pnt2;
		  trim_segment->closestPoint(seg_pnt, 
					     trim_segment->startparam(),
					     trim_segment->endparam(),
					     par, seg_pnt2, dist);
		  if (dist < init_dist)
			{
			  init_dist = dist;
			  segment_point = seg_pnt2;
			}
		  break;
		}
	    }

#ifdef DEBUG_REG
	  std::ofstream out_file("split_segments.g2");
	  for (size_t kj=0; kj<segments.size(); ++kj)
	    {
	      shared_ptr<ParamCurve> cv = segments[kj]->spaceCurve();
	      cv->writeStandardHeader(out_file);
	      cv->write(out_file);
	    }
#endif

	  if (ki < corner.size())
	    {
	      // Isolate hole
	      if (segment_point.size() > 0)
		faces = isolateHole(segment_point, corner[ki], half_holes); 
	      else
		{
		  while (faces.size() == 0 && isolate_fac_ > 0.3)
		    {
		      isolate_fac_ *= 0.75;
		      faces = initIsolateHole(half_holes);
		    }
		}
	    }
	  else if (!OKseg && corner.size() > 4)
	    {
	      // Try to divided according to the outer boundary before
	      // the hole is removed
	      unsetAxis();
	      faceOuterBd(half_holes);
	    }
	  else
	    {
	      // Divide out hole
	      // Define faces
	      vector<shared_ptr<BoundedSurface> > sub_sfs =
		BoundedUtils::splitWithTrimSegments(bd_sf, segments, epsge_);

#ifdef DEBUG_REG
	      std::ofstream of("split_surf.g2");
	      for (size_t kr=0; kr<sub_sfs.size(); ++kr)
		{
		  sub_sfs[kr]->writeStandardHeader(of);
		  sub_sfs[kr]->write(of);
		}
#endif
      
	      // Fetch info about vertices not belonging to the corners of the
	      // initial face
	      vector<shared_ptr<Vertex> > non_corner = 
		face_->getNonCornerVertices(bend_);
	      removeInsignificantVertices(non_corner);

	      faces = RegularizeUtils::createFaces(sub_sfs, face_, epsge_, tol2_,
						   angtol_, non_corner);
	    }
	}
    }

  // Update topology and treat each sub face
  int nmb_faces = (int)faces.size();
  if (nmb_faces > 1)
    {
      model_->removeFace(face_);
      model_->append(faces, false, false);
      for (int ki=0; ki<nmb_faces; )
	{
	  RegularizeFace regularize(faces[ki], model_);
	  if (axis_.dimension() > 0)
	    regularize.setAxis(centre_, axis_);
	  regularize.setDivideInT(divideInT_);
	  regularize.setCandSplit(cand_split_);
	  regularize.setSplitMode(split_mode_);
	  regularize.unsetTopLevel();
	  if (cand_split_.size() >  0)
	    regularize.setCandSplit(cand_split_);

	  vector<shared_ptr<ftSurface> > faces2 = 
	    regularize.getRegularFaces();

	  if (faces2.size() > 1)
	    {
	      faces.erase(faces.begin()+ki);
	      nmb_faces--;
	      faces.insert(faces.end(), faces2.begin(), faces2.end());

	      // Check if any new faces may be joined across the seam
	      mergeSeams(faces, nmb_faces, faces2);
	    }
	  else
	    ki++;
	}
      sub_faces_.insert(sub_faces_.end(), faces.begin(), faces.end());
    }
}


//==========================================================================
shared_ptr<CurveOnSurface> 
RegularizeFace::computeCornerSplit(shared_ptr<Vertex> corner,
				   vector<shared_ptr<Vertex> >& hole_vx,
				   vector<shared_ptr<Vertex> >& hole_vx2,
				   shared_ptr<BoundedSurface>& bd_sf,
				   const Point& close_par, bool outer_vx)
//==========================================================================
{
  shared_ptr<CurveOnSurface> dummy;
  shared_ptr<ParamSurface> surf = face_->surface();
  

  // Define plane
  double plane_ang = M_PI/4.0;
  Point pnt = corner->getVertexPoint();
  Point norm;
  if (centre_.dimension() > 0)
    {
      Point vec = pnt - centre_;
      vec.normalize();
      double angle = vec.angle(axis_);
      norm = (pnt - centre_).cross(axis_);
      norm.normalize();
      if (angle < plane_ang)
	{
	  norm = vec.cross(norm);
	  norm.normalize();
	}
    }
  else
    {
      norm = axis_;
      norm.normalize();
    }

  // Check if an inner vertex approximately defines the same plane
  double level_ang = M_PI/8.0; //M_PI/15.0; //M_PI/10.0;
  double level_ang2 = M_PI/75.0; //M_PI/75.0; // M_PI/100.0; //M_PI/50.0; //M_PI/10.0;
  int min_idx = -1, min_idx2 = -1;
  double min_ang = 2.0*level_ang, min_ang2 = 2.0*level_ang;
  double fac = 0.5;

  if (hole_vx.size() > 0)
    {
      for (size_t kj=0; kj<hole_vx.size(); ++kj)
	{
	  double ang = getSegmentAngle(corner, hole_vx[kj], 
				       pnt, norm);
	  ang = std::min(ang, fabs(M_PI-ang));
	  Point vx_pt = hole_vx[kj]->getVertexPoint();
	  double side = (vx_pt - centre_)*(pnt - centre_);
	  double dist = fabs((vx_pt - pnt)*norm);
	  if (ang < min_ang && side > 0.0 && dist < fac*radius_)
	    {
	      min_ang = ang;
	      min_idx = (int)kj;
	    }
	}
    }

  if (hole_vx2.size() > 0)
    {
      for (size_t kj=0; kj<hole_vx2.size(); ++kj)
	{
	  double ang = getSegmentAngle(corner, hole_vx2[kj], 
				       pnt, norm);
	  ang = std::min(ang, fabs(M_PI-ang));
	  Point vx_pt = hole_vx2[kj]->getVertexPoint();
	  double side = (vx_pt - centre_)*(pnt - centre_);
	  double dist = fabs((vx_pt - pnt)*norm);
	  if (ang < min_ang2 && side > 0.0 && dist < fac*radius_)
	    {
	      min_ang2 = ang;
	      min_idx2 = (int)kj;
	    }
	}
    }

  // Check if there exists a pattern of division from previous surfaces
  bool found = false;
  Point parval1, parval2;
  if (cand_split_.size() > 0)
    {
      // Fetch an appropriate split curve
      Point vx_point = corner->getVertexPoint();
      found = fetchPatternSplit(vx_point, parval1, parval2);
    }

  vector<shared_ptr<CurveOnSurface> > trim_segments;
  Point pnt2;
  bool OKseg = true;
  shared_ptr<CurveOnSurface> save_segment;
  if (trim_segments.size() == 0 && min_idx >= 0 && min_ang < level_ang)
    {
      // A corresponding vertex if found. Split between vertices
      parval1 = corner->getFacePar(face_.get());
      pnt2 = corner->getVertexPoint();
      parval2 = hole_vx[min_idx]->getFacePar(face_.get());
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
						     parval2, epsge_,
						     bd_sf);

      // Remove vertex from candidate list so it will not be used 
      // again
      hole_vx.erase(hole_vx.begin() + min_idx);

#ifdef DEBUG_REG
      std::ofstream out_file_1("trim_segments.g2");
      for (size_t kv=0; kv<trim_segments.size(); ++kv)
	{
	  shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
	  cv->writeStandardHeader(out_file_1);
	  cv->write(out_file_1);
	}
#endif
      OKseg = checkTrimSegments(trim_segments, corner, pnt, bd_sf, outer_vx);
      if (!OKseg && !save_segment.get() && trim_segments.size() == 1)
	{
	  save_segment = trim_segments[0];
	  trim_segments.clear();
	}
    }

  if (trim_segments.size() == 0 && min_idx2 >= 0 && min_ang2 < level_ang2)
    {
      // A corresponding vertex if found. Split between vertices
      parval1 = corner->getFacePar(face_.get());
      pnt2 = corner->getVertexPoint();
      parval2 = hole_vx2[min_idx2]->getFacePar(face_.get());
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
						     parval2, epsge_,
						     bd_sf);

      // Remove vertex from candidate list so it will not be used 
      // again
      hole_vx2.erase(hole_vx2.begin() + min_idx2);

#ifdef DEBUG_REG
      std::ofstream out_file_1("trim_segments.g2");
      for (size_t kv=0; kv<trim_segments.size(); ++kv)
	{
	  shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
	  cv->writeStandardHeader(out_file_1);
	  cv->write(out_file_1);
	}
#endif
      OKseg = checkTrimSegments(trim_segments, corner, pnt, 
				      bd_sf, outer_vx);
      if (OKseg && save_segment.get() && trim_segments.size() == 1)
	save_segment.reset();
      else if (!OKseg && !save_segment.get() && trim_segments.size() == 1)
	{
	  save_segment = trim_segments[0];
	  trim_segments.clear();
	}
	  
    }

  if (trim_segments.size() == 0 && close_par.dimension() > 0)
    {
      // Connect to closest point
      parval1 = corner->getFacePar(face_.get());
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
						     close_par, epsge_,
						     bd_sf);
#ifdef DEBUG_REG
      std::ofstream out_file_1("trim_segments.g2");
      for (size_t kv=0; kv<trim_segments.size(); ++kv)
	{
	  shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
	  cv->writeStandardHeader(out_file_1);
	  cv->write(out_file_1);
	}
#endif
      OKseg = checkTrimSegments(trim_segments, corner, pnt, bd_sf, outer_vx);
      if (OKseg && save_segment.get() && trim_segments.size() == 1)
	save_segment.reset();
      else if (!OKseg && !save_segment.get() && trim_segments.size() == 1)
	{
	  save_segment = trim_segments[0];
	  trim_segments.clear();
	}

      if (trim_segments.size() == 0)
	{
	  // Try to define a curve in the parameter domain of the surface
	  Point close_pos = face_->point(close_par[0], close_par[1]);
	  shared_ptr<ParamCurve> pcrv = 
	    RegularizeUtils::checkStrightParCv(face_, corner, close_pos, epsge_);
	  if (pcrv.get())
	    {
	      // Create corresponding curve-on-surface curve
	      trim_segments = BoundedUtils::getTrimCrvsPcrv(surf, pcrv, epsge_, bd_sf);
	      if (trim_segments.size() > 1)
		trim_segments.clear();   // No valid splitting curve
#ifdef DEBUG_REG
	      std::ofstream out_file_2("trim_segments0.g2");
	      for (size_t kv=0; kv<trim_segments.size(); ++kv)
		{
		  shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
		  cv->writeStandardHeader(out_file_2);
		  cv->write(out_file_2);
		}
#endif
	      OKseg = checkTrimSegments(trim_segments, corner, pnt, bd_sf, 
					outer_vx);
	      if (OKseg && save_segment.get() && trim_segments.size() == 1)
		save_segment.reset();
	      else if (!OKseg && !save_segment.get() && trim_segments.size() == 1)
		{
		  save_segment = trim_segments[0];
		  trim_segments.clear();
		}
	    }
	}
    }

  if (trim_segments.size() == 0 && found)
    {
      // Use split pattern
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
						     parval2, epsge_,
						     bd_sf);
	  
#ifdef DEBUG_REG
      std::ofstream out_file_1("trim_segments.g2");
      for (size_t kv=0; kv<trim_segments.size(); ++kv)
	{
	  shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
	  cv->writeStandardHeader(out_file_1);
	  cv->write(out_file_1);
	}
#endif
      OKseg = checkTrimSegments(trim_segments, corner, pnt, bd_sf, outer_vx);
      if (OKseg && save_segment.get() && trim_segments.size() == 1)
	save_segment.reset();
      else if (!OKseg && !save_segment.get() && trim_segments.size() == 1)
	{
	  save_segment = trim_segments[0];
	  trim_segments.clear();
	}
    }

    if (trim_segments.size() == 0)
    {
      int kt=0;
      while (trim_segments.size() == 0 && kt<2)
	{
	  kt++;

	  // Split by plane intersection
	  pnt2 = pnt;
	  trim_segments = BoundedUtils::getPlaneIntersections(surf, pnt,
							      norm, epsge_,
							      bd_sf);
	  OKseg = checkTrimSegments(trim_segments, corner, pnt, bd_sf, 
				    outer_vx);
	  if (OKseg && save_segment.get() && trim_segments.size() == 1)
	    save_segment.reset();
	  else if (!OKseg && !save_segment.get() && trim_segments.size() == 1)
	    {
	      save_segment = trim_segments[0];
	      trim_segments.clear();
	    }

#ifdef DEBUG_REG
	  std::ofstream out_file_1("trim_segments.g2");
	  for (size_t kv=0; kv<trim_segments.size(); ++kv)
	    {
	      shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
	      cv->writeStandardHeader(out_file_1);
	      cv->write(out_file_1);
	    }
#endif

	  norm = (pnt - centre_).cross(norm);
	  norm.normalize();
	}
    }

    if (trim_segments.size() == 0)
      {
	// Try to define a curve in the parameter domain of the surface
	shared_ptr<ParamCurve> pcrv = 
	  RegularizeUtils::checkStrightParCv(face_, corner, centre_, epsge_);
	if (pcrv.get())
	  {
	    // Create corresponding curve-on-surface curve
	    trim_segments = BoundedUtils::getTrimCrvsPcrv(surf, pcrv, epsge_, bd_sf);
	    if (trim_segments.size() > 1)
	      trim_segments.clear();   // No valid splitting curve
#ifdef DEBUG_REG
	  std::ofstream out_file_2("trim_segments0.g2");
	  for (size_t kv=0; kv<trim_segments.size(); ++kv)
	    {
	      shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
	      cv->writeStandardHeader(out_file_2);
	      cv->write(out_file_2);
	    }
#endif
	  OKseg = checkTrimSegments(trim_segments, corner, pnt, bd_sf, 
				    outer_vx);
	  if (OKseg && save_segment.get() && trim_segments.size() == 1)
	    save_segment.reset();
	  else if (!OKseg && !save_segment.get() && trim_segments.size() == 1)
	    {
	      save_segment = trim_segments[0];
	      trim_segments.clear();
	    }
#ifdef DEBUG_REG
	  std::ofstream out_file_3("trim_segments.g2");
	  for (size_t kv=0; kv<trim_segments.size(); ++kv)
	    {
	      shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
	      cv->writeStandardHeader(out_file_3);
	      cv->write(out_file_3);
	    }
#endif
	  }
      }

  if (trim_segments.size() != 1)
    {
      if (save_segment.get())
	return save_segment;
      else
	return dummy;  // Either no legal splitting curve is found
  // or it is divided into several pieces. Thus, another split
  // prior to this one might be favorable
    }
  else
    return trim_segments[0];
}

//==========================================================================
  bool 
    RegularizeFace::checkTrimSegments(vector<shared_ptr<CurveOnSurface> >& trim_segments,
				      shared_ptr<Vertex> corner, Point pnt,
				      shared_ptr<BoundedSurface>& bd_sf,
				      bool outer_vx)
//==========================================================================
{
  shared_ptr<ParamSurface> surf = face_->surface();
  bool OKseg = true;

  // Remove intersections not connected with the initial point and
  // intersections going in the opposite direction than to the centre
  int kr;
  for (kr=0; kr<(int)trim_segments.size(); )
    {
      Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
      Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());
      double d1 = pos1.dist(pnt);
      double d2 = pos2.dist(pnt);
      double d3 = 0, d4 = 0;
      if (centre_.dimension() > 0)
	{
	  d3 = pos1.dist(centre_) - 2.0*radius_;
	  d4 = pos2.dist(centre_) - 2.0*radius_;
	}
      if (d1 > epsge_ && d2 > epsge_)
	trim_segments.erase(trim_segments.begin()+kr);
      // else if (d3 > epsge_ && d4 > epsge_)
      // 	trim_segments.erase(trim_segments.begin()+kr);
      else if (outer_vx &&
	       ((d1 < epsge_ && pos2.dist(centre_) > pos1.dist(centre_)) ||
		(d2 < epsge_ && pos1.dist(centre_) > pos2.dist(centre_))))
	trim_segments.erase(trim_segments.begin()+kr);
      else
	{
	  // Check the angle between the splitting curve and the boundary curves
	  // in the current corner
	  vector<Point> der(2);
	  if (d1 < d2)
	    trim_segments[kr]->point(der, trim_segments[kr]->startparam(), 1);
	  else
	    trim_segments[kr]->point(der, trim_segments[kr]->endparam(), 1);

	  vector<ftEdge*> edges = corner->getFaceEdges(face_.get());
	  if (edges.size() == 2)
	    {
	      double t1 = edges[0]->parAtVertex(corner.get());
	      double t2 = edges[1]->parAtVertex(corner.get());
	      Point tan1 = edges[0]->tangent(t1);
	      Point tan2 = edges[1]->tangent(t2);
	      double ang1 = der[1].angle(tan1);
	      if (fabs(M_PI-ang1) < ang1)
		ang1 = fabs(M_PI-ang1);
	      double ang2 = der[1].angle(tan2);
	      if (fabs(M_PI-ang2) < ang2)
		ang2 = fabs(M_PI-ang2);
	      if (ang1 < bend_ || ang2 < bend_)
		{
		  // Replace the splitting curve with a parameter based curve
		  Point param = 
		    trim_segments[kr]->faceParameter((d1<d2) ? trim_segments[kr]->endparam()
						     : trim_segments[kr]->startparam());
		  Point parval1 = corner->getFacePar(face_.get());
		  vector<shared_ptr<CurveOnSurface> > tmp_segments = 
		    BoundedUtils::getTrimCrvsParam(surf, parval1, param, epsge_, bd_sf);

		  if (tmp_segments.size() > 0)
		    {
		      // Check if the alternative is any better
		      vector<Point> der2(2);
		      tmp_segments[0]->point(der2, tmp_segments[0]->startparam(), 1);
		      double ang3 = der2[1].angle(tan1);
		      if (fabs(M_PI-ang3) < ang3)
			ang3 = fabs(M_PI-ang3);
		      double ang4 = der2[1].angle(tan2);
		      if (fabs(M_PI-ang4) < ang4)
			ang4 = fabs(M_PI-ang4);
		      if (std::min(ang1, ang2) < std::min(ang3, ang4) &&
			  ang3 >= bend_ && ang4 >= bend_)
			trim_segments[kr] = tmp_segments[0];
		      else 
			{
			  OKseg = false;
			  // trim_segments.erase(trim_segments.begin()+kr);
			  // kr--;
			}
		    }
		} 
	    }
	  kr++;
	}
    }
  return OKseg;
}

//==========================================================================
void RegularizeFace::faceOneHole2()
//==========================================================================
{
  ASSERT(face_->nmbBoundaryLoops() == 2);

  // // Search for half holes in the inner loop
  // vector<vector<ftEdge*> > half_holes = getHalfHoles(1);

  // if (half_holes.size() == 1)
  //   {
  //   }
  // else if (half_holes.size() > 1)
  //   {
  //     faceWithHoles(half_holes);
  //   }
  // else
  //   {
      // Fetch significant vertices in the inner loop
      vector<shared_ptr<Vertex> > corner = face_->getCornerVertices(bend_, 1);
      vector<shared_ptr<Vertex> > vxs = face_->getBoundaryLoop(0)->getVertices();
      removeInsignificantVertices(vxs);
      if (corner.size() > 0)
	{
	  // Split in corners
	  vector<shared_ptr<CurveOnSurface> > segments;
	  shared_ptr<BoundedSurface> bd_sf;
	  shared_ptr<ParamSurface> surf = face_->surface();
	  vector<shared_ptr<Vertex> > dummy_vx;  // No vertices in outer loop
	  for (size_t ki=0; ki<corner.size(); ++ki)
	    {
	      Point pnt = corner[ki]->getVertexPoint();
	      Point vec = pnt - centre_;
	      vec.normalize();
	      
	      // Fetch a vector in the given vertex pointing into the surface
	      Point in_vec = RegularizeUtils::getInVec(corner[ki], face_);

	      // Compute closest point on the outer boundary
	      //int close_idx;
	      double close_par, close_dist, upar, vpar;
	      Point close_pnt;
	      (void)face_->closestOuterBoundaryPoint(pnt, in_vec, upar, vpar, 
						     close_pnt, close_dist, close_par);
	      Point face_par(upar, vpar);
	      
	      double fac = 0.8;
	      double lim = fac*pnt.dist(close_pnt);
	      snapToVertex(close_pnt, face_par, vxs, lim);

#ifdef DEBUG_REG
	      std::ofstream of0("vx_close.g2");
	      of0 << "400 1 0 4 255 0 0 255" << std::endl;
	      of0 << "1" << std::endl;
	      of0 << pnt << std::endl;
	      of0 << "400 1 0 4 0 255 0 255" << std::endl;
	      of0 << "1" << std::endl;
	      of0 << close_pnt << std::endl;
#endif

	      // Compute splitting curve
	      Point dummy;
	      shared_ptr<CurveOnSurface> trim_segment = 
		computeCornerSplit(corner[ki], dummy_vx, dummy_vx, bd_sf,
				   face_par, false);
	      if (trim_segment.get())
		segments.push_back(trim_segment);
	    }

	  if (segments.size() > 0)
	    {
	      // Divide out hole
	      // Define faces
	      vector<shared_ptr<BoundedSurface> > sub_sfs =
		BoundedUtils::splitWithTrimSegments(bd_sf, segments, epsge_);

#ifdef DEBUG_REG
	      std::ofstream of("split_surf.g2");
	      for (size_t kr=0; kr<sub_sfs.size(); ++kr)
		{
		  sub_sfs[kr]->writeStandardHeader(of);
		  sub_sfs[kr]->write(of);
		}
#endif

	      vector<shared_ptr<Vertex> > dummy_vx;
	      vector<shared_ptr<ftSurface> > faces = 
		RegularizeUtils::createFaces(sub_sfs, face_, epsge_, tol2_,
					     angtol_, dummy_vx);

	      // Update topology and treat each sub face
	      int nmb_faces = (int)faces.size();
	      if (nmb_faces > 1)
		{
		  model_->removeFace(face_);
		  model_->append(faces, false, false);
		  for (int ki=0; ki<nmb_faces; )
		    {
		      RegularizeFace regularize(faces[ki], model_);
		      if (centre_.dimension() > 0)
			regularize.setAxis(centre_, axis_);
		      regularize.setDivideInT(divideInT_);
		      regularize.setCandSplit(cand_split_);
		      regularize.setSplitMode(split_mode_);
		      regularize.unsetTopLevel();
		      if (cand_split_.size() >  0)
			regularize.setCandSplit(cand_split_);
		  
		      vector<shared_ptr<ftSurface> > faces2 = 
			regularize.getRegularFaces();
		  
		      if (faces2.size() > 1)
			{
			  faces.erase(faces.begin()+ki);
			  nmb_faces--;
			  faces.insert(faces.end(), faces2.begin(), faces2.end());
			  // Check if any new faces may be joined across the seam
			  mergeSeams(faces, nmb_faces, faces2);
			}
		      else
			ki++;
		    }
		}
	      sub_faces_.insert(sub_faces_.end(), faces.begin(), faces.end());
	    }
	}
      else
	{
	  // Fetch other vertices in the inner loop
	  vector<shared_ptr<Vertex> > vx = face_->getNonCornerVertices(bend_, 1);
	  if (vx.size() > 0)
	    {
	      // Split in vertice(s)
	    }
	  else
	    {
	      // Split in one arbitrary point on the outer boundary
	    }
	}
      //    }
    
}


//==========================================================================
void RegularizeFace::snapToVertex(Point& pos, Point& par, 
				  vector<shared_ptr<Vertex> >& vxs,
				  double lim)
//==========================================================================
{
  size_t kr;
  int ix = -1;
  double mindist = HUGE;
  for (kr=0; kr<vxs.size(); ++kr)
    {
      double dist = vxs[kr]->getVertexPoint().dist(pos);
      if (dist < mindist && dist <= lim)
	{
	  ix = (int)kr;
	  mindist = dist;
	}
    }

  if (ix >= 0)
    {
      pos = vxs[ix]->getVertexPoint();
      par = vxs[ix]->getFacePar(face_.get());
    }
}


//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::faceOuterBdFaces(vector<vector<ftEdge*> >& half_holes)
//==========================================================================
{
  vector<shared_ptr<ftSurface> > subfaces;
  vector<shared_ptr<Vertex> > non_corner;

  // The face is to be treated according to its outer boundary
  int nmb_loops = face_->nmbBoundaryLoops();
  if (nmb_loops == 0)
    return subfaces;  // Nothing to do

  // Find the corner vertex with largest angle
  // First find corners of the outer boundary
  vector<shared_ptr<Vertex> > corners = 
    face_->getCornerVertices(bend_, 0);

  // Dismiss vertices belonging to the half holes
  size_t ki;
  for (ki=0; ki<half_holes.size(); ++ki)
    {
      shared_ptr<Vertex> curr_vx = half_holes[ki][0]->getVertex(true);
      vector<shared_ptr<Vertex> >::iterator vxp =
	std::find(corners.begin(), corners.end(), curr_vx);
      if (vxp != corners.end())
	corners.erase(vxp);
      for (size_t kj=0; kj<half_holes[ki].size(); ++kj)
	{
	  curr_vx = half_holes[ki][kj]->getVertex(false);
	  vxp = std::find(corners.begin(), corners.end(), curr_vx);
	  if (vxp != corners.end())
	    corners.erase(vxp);
	}
    }
  
  // Fetch all concave corner vertices
  vector<shared_ptr<Vertex> > concave_corners;
  getConcaveCorners(corners, concave_corners);
  if (concave_corners.size() > 0)
    {
#ifdef DEBUG_REG
      std::ofstream of01("concave_vx.g2");
      of01 << "400 1 0 4 255 0 0 255 " << std::endl;
      of01 << concave_corners.size() << std::endl;
      for (size_t k3=0; k3<concave_corners.size(); ++k3)
	of01 << concave_corners[k3]->getVertexPoint() << std::endl;
#endif
    }
  
  if (concave_corners.size() > 0)
    {
      // Make possible regular blocks between vertices
      subfaces = chopOffRegBlocks(concave_corners);
      if (subfaces.size() > 1)
	return subfaces;
    }

  if (concave_corners.size() > 0)
    {
      // Split by making connections between concave corners and
      // other vertices
      subfaces = connectToVertex(concave_corners);
      if (subfaces.size() > 1)
	return subfaces;
    }

   if (corners.size() == 0)
    return subfaces;
   //  shared_ptr<Vertex> split_vx = getSignificantVertex(corners);
   // Prioritize corners according to opening angle
   vector<shared_ptr<Vertex> > sorted_corners = prioritizeCornerVx(corners);

  vector<shared_ptr<CurveOnSurface> > trim_segments;
  shared_ptr<BoundedSurface> bd_sf;
  for (ki=0; ki<sorted_corners.size(); ++ki)
     {
       shared_ptr<Vertex> split_vx = sorted_corners[ki];

       // Divide the face according to a plane/parameter line trough this
       // vertex, dividing the angle (approximately) in two equal pieces

       // Prefer to divide between two vertices, but divide to an edge
       // if this seems to be more adequate

       // Traverse all vertices to find preferred candidates
       // corners or all?
       vector<shared_ptr<Vertex> > cand_vx;
       ftEdge* cand_edge = 0;
       vector<shared_ptr<Vertex> > vx = face_->getBoundaryLoop(0)->getVertices();
       // Dismiss vertices belonging to the half holes
       for (size_t kr=0; kr<half_holes.size(); ++kr)
	 {
	   shared_ptr<Vertex> curr_vx = half_holes[kr][0]->getVertex(true);
	   vector<shared_ptr<Vertex> >::iterator vxp =
	     std::find(vx.begin(), vx.end(), curr_vx);
	   if (vxp != vx.end())
	     vx.erase(vxp);
	   for (size_t kj=0; kj<half_holes[kr].size(); ++kj)
	     {
	       curr_vx = half_holes[kr][kj]->getVertex(false);
	       vxp = std::find(vx.begin(), vx.end(), curr_vx);
	       if (vxp != vx.end())
		 vx.erase(vxp);
	     }
	 }

       selectCandidateSplit(split_vx, vx, cand_vx, cand_edge);

#ifdef DEBUG_REG
       std::ofstream of0("split_vx.g2");
       of0 << "400 1 0 4 255 0 0 255 " << std::endl;
       of0 << "1" << std::endl;
       of0 << split_vx->getVertexPoint() << std::endl;
#endif

       // Fetch info about vertices not belonging to the corners of the
       // initial face
       non_corner = face_->getNonCornerVertices(bend_);
       removeInsignificantVertices(non_corner);

       // Perform split
       vector<shared_ptr<Vertex> > prio;
       for (size_t kr=0; kr<corners_.size(); ++kr)
	 {
	   size_t kh;
	   for (kh=0; kh<cand_vx.size(); ++kh)
	     if (cand_vx[kh].get() == corners_[kr].get())
	       {
		 prio.push_back(cand_vx[kh]);
		 cand_vx.erase(cand_vx.begin()+kh);
		 break;
	       }
	 }
       trim_segments = RegularizeUtils::findVertexSplit(face_, split_vx, cand_vx, 
							cand_edge, prio, epsge_, 
							tol2_, angtol_, bend_, 
							non_corner, centre_, axis_,
							bd_sf);

       // Check split configuration to avoid 3-sided surfaces
       RegularizeUtils::checkTrimConfig(face_, trim_segments, split_vx,
					corners, tol2_);
       
       if (trim_segments.size() > 0)
	 break;
     }

  if (trim_segments.size() == 0)
    return subfaces;   // No splitting is performed

  // Define faces
   vector<shared_ptr<BoundedSurface> > sub_sfs =
     BoundedUtils::splitWithTrimSegments(bd_sf, trim_segments, epsge_);

#ifdef DEBUG_REG
   std::ofstream of("split_surf.g2");
   for (size_t kr=0; kr<sub_sfs.size(); ++kr)
     {
       sub_sfs[kr]->writeStandardHeader(of);
       sub_sfs[kr]->write(of);
     }
#endif

   // Create faces
   subfaces = RegularizeUtils::createFaces(sub_sfs, face_, epsge_,
					   tol2_, angtol_, non_corner);
   return subfaces;
}

//==========================================================================
void RegularizeFace::faceOuterBd(vector<vector<ftEdge*> >& half_holes)
//==========================================================================
{
  vector<shared_ptr<ftSurface> > subfaces = faceOuterBdFaces(half_holes);

#ifdef DEBUG_REG
  std::ofstream of("split_face.g2");
  shared_ptr<ParamSurface> surf = face_->surface();
  surf->writeStandardHeader(of);
  surf->write(of);
  for (size_t kr=0; kr<subfaces.size(); ++kr)
    {
      surf = subfaces[kr]->surface();
      surf->writeStandardHeader(of);
      surf->write(of);
    }
#endif

  // Update topology and treat each sub face
  int nmb_faces = (int)subfaces.size();
  if (nmb_faces > 1)
    {
      model_->removeFace(face_);
      model_->append(subfaces, false, false);
      for (int ki=0; ki<nmb_faces; )
	{
	  RegularizeFace regularize(subfaces[ki], model_);
	  if (axis_.dimension() > 0)
	    regularize.setAxis(centre_, axis_);
	  regularize.setDivideInT(divideInT_);
	  regularize.setCandSplit(cand_split_);
	  regularize.setSplitMode(split_mode_);
	  regularize.unsetTopLevel();
	  vector<shared_ptr<ftSurface> > faces = 
	    regularize.getRegularFaces();

	  if (faces.size() > 1)
	    {
	      subfaces.erase(subfaces.begin()+ki);
	      nmb_faces--;
	      subfaces.insert(subfaces.end(), faces.begin(), faces.end());

	      // Check if any new faces may be joined across the seam
	      mergeSeams(subfaces, nmb_faces, faces);
	    }
	  else
	    ki++;
	}
    }
  sub_faces_.insert(sub_faces_.end(), subfaces.begin(), subfaces.end());
}

//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeFace::divideVertex(ftEdge* edge, double par,
			     vector<shared_ptr<Vertex> > cand_vx,
			     ftEdge* cand_edge)
//==========================================================================
{
  // Get the plane with which to divide the current face to get subdivision
  // information
  Point pnt = edge->point(par);
  Point normal = edge->tangent(par);
  Point norm = edge->normal(par);
  Point vec2 = normal.cross(norm);

  // Traverse all candidate vertices and check if any is feasible for 
  // division
  // Compute angle between segment between vertices and found plane
  // Select vertex with minimum angle
  double level_ang = M_PI/8.0;
  int min_idx = -1;
  double min_ang = 2.0*level_ang;

for (int ki=0; ki<(int)cand_vx.size(); ++ki)
  {
    Point vec = cand_vx[ki]->getVertexPoint() - pnt;
    double ang = vec.angle(vec2);
    ang = std::min(ang, M_PI-ang);
    if (ang < min_ang)
      {
	min_ang = ang;
	min_idx = ki;
      }
  }

  vector<shared_ptr<CurveOnSurface> > trim_segments;
  shared_ptr<BoundedSurface> bd_sf;
  shared_ptr<ParamSurface> surf = face_->surface();
  if (min_ang < level_ang)
    {
      // Find division curve between vertices
      Point parval1 = edge->faceParameter(par);
      Point parval2 = cand_vx[min_idx]->getFacePar(face_.get());
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
						     parval2, epsge_,
						     bd_sf);
    }
  else if (cand_edge)
    {
      // Let the division curve end at the midpoint of the given edge
      Point parval1 = edge->faceParameter(par);
      double tmid = 0.5*(cand_edge->tMin() + cand_edge->tMax());
      Point parval2 = cand_edge->faceParameter(tmid);
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
						     parval2, epsge_,
						     bd_sf);
    }
  else
    {
      // Check for a split pattern
      bool found = false;
      if (cand_split_.size() > 0)
	{
	  // Fetch an appropriate split curve
	  Point parval1, parval2;
	  found = fetchPatternSplit(pnt, parval1, parval2, false);

	  if (found)
	    {
	      // Check if both endpoints lie at the outer boundary of the face
	      Point pt1 = face_->point(parval1[0], parval1[1]);
	      Point pt2 = face_->point(parval2[0], parval2[1]);
	      shared_ptr<Loop> bd = face_->getBoundaryLoop(0);
	      int idx1, idx2;
	      double par1, par2, d1, d2;
	      Point clo1, clo2;
	      bd->closestPoint(pt1, idx1, par1, clo1, d1);
	      bd->closestPoint(pt2, idx2, par2, clo2, d2);
	      if (d1 > epsge_ || d2 > epsge_)
		{
		  if (d1 <= epsge_)
		    {
		      // Reposition plane to fit with pattern from opposite
		      // face
		      pnt = face_->point(parval1[0], parval1[1]);
		    }
		  found = false;
		}
	    }
	  if (found)
	    trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
							   parval2, epsge_,
							   bd_sf);
	  if (trim_segments.size() == 0)
	    found = false;
	    
	}
      if (!found)
	// Find intersections between the face and this plane
	trim_segments = BoundedUtils::getPlaneIntersections(surf, pnt,
							    normal, epsge_,
							    bd_sf);
    }

  // Check feasability of intersections
  //MESSAGE("Some code is missing to check feasability of intersection");

#ifdef DEBUG_REG
  std::ofstream out_file("split_segments.g2");
  for (size_t kj=0; kj<trim_segments.size(); ++kj)
    {
      shared_ptr<ParamCurve> cv = trim_segments[kj]->spaceCurve();
      cv->writeStandardHeader(out_file);
      cv->write(out_file);
    }
#endif

  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_sf, trim_segments, epsge_);

  // Fetch info about vertices not belonging to the corners of the
  // initial face
  vector<shared_ptr<Vertex> > non_corner = 
    face_->getNonCornerVertices(bend_);
  removeInsignificantVertices(non_corner);

  vector<shared_ptr<ftSurface> > faces = 
    RegularizeUtils::createFaces(sub_sfs, face_, epsge_,  tol2_, angtol_, 
				 non_corner);

  return faces;
}

//==========================================================================
double RegularizeFace::getSegmentAngle(shared_ptr<Vertex> vx1,
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
  vector<ftEdge*> edges = vx1->getFaceEdges(face_.get());
  if (edges.size() < 1)
      ang += M_PI;
  else
    {
      int ix0 = 0, ix1 = 1;
      if (edges[0]->prev() == edges[1])
	std::swap(ix0, ix1);
      double t1 = edges[ix0]->parAtVertex(vx1.get());
      double t2 = edges[ix1]->parAtVertex(vx1.get());
      Point tan1 = edges[ix0]->tangent(t1);
      Point tan2 = edges[ix1]->tangent(t2);
      if (edges[ix0]->tMax() - t1 < t1 - edges[ix0]->tMin())
	tan1 *= -1;
      else 
	tan2 *= -1;
      tan1.normalize();
      tan2.normalize();
      Point vec2 = 0.5*(tan1 + tan2);
      Point par = vx1->getFacePar(face_.get());
      Point norm1 = face_->normal(par[0], par[1]);
      Point norm2 = tan2.cross(tan1);
      if (norm1*norm2 < 0.0)
	vec2 *= -1;
      if (vec*vec2 < 0.0)
	ang += M_PI;
    }
      
  return ang;
}

//==========================================================================
shared_ptr<Vertex> 
RegularizeFace::getSignificantVertex(vector<shared_ptr<Vertex> > cand_vx)
//==========================================================================
{
  vector<shared_ptr<Vertex> > cand_vx2;
  if (cand_split_.size() > 0)
    {
      // Compute minimum distance between candidate vertices
      size_t ki, kj;
      double min_dist = 1.0e8;  // A huge number
      for (ki=0; ki<cand_vx.size(); ++ki)
	{
	  Point vxp1 = cand_vx[ki]->getVertexPoint();
	  for (kj=ki+1; kj<cand_vx.size(); ++kj)
	    min_dist = std::min(min_dist, 
				vxp1.dist(cand_vx[kj]->getVertexPoint()));
	}
      double level_dist = 0.2*min_dist;

      // Select vertices which corresponds to T-joints in the corresponding 
      // face
      vector<pair<Point,int> > corr_pts;
      for (ki=0; ki<cand_split_.size(); ++ki)
	{
	  for (kj=0; kj<corr_pts.size(); ++kj)
	    if (corr_pts[kj].first.dist(cand_split_[ki].first.first) < epsge_)
	      break;
	  if (kj < corr_pts.size())
	    corr_pts[kj].second++;
	  else
	    corr_pts.push_back(make_pair(cand_split_[ki].first.first, 1));

	  for (kj=0; kj<corr_pts.size(); ++kj)
	    if (corr_pts[kj].first.dist(cand_split_[ki].second.first) < epsge_)
	      break;
	  if (kj < corr_pts.size())
	    corr_pts[kj].second++;
	  else
	    corr_pts.push_back(make_pair(cand_split_[ki].second.first, 1));
	}

      for (ki=0; ki<corr_pts.size(); ++ki)
	{
	  if (corr_pts[ki].second < 3)
	    continue;   // Not a T-joint vertex

	  // Fetch corresponding point in this face
	  Point pos1;
	  double u1, v1, d1;
	  face_->closestPoint(corr_pts[ki].first, u1, v1, pos1, d1, epsge_);

	  // Find closest candidate vertex
	  double min_dist2 = min_dist;
	  size_t min_ind = 0;
	  for (kj=0; kj<cand_vx.size(); ++kj)
	    {
	      double dist = pos1.dist(cand_vx[kj]->getVertexPoint());
	      if (dist < min_dist2)
		{
		  min_dist2 = dist;
		  min_ind = kj;
		}
	    }
	  
	  if (min_dist2 < level_dist)
	    cand_vx2.push_back(cand_vx[min_ind]);
	}
      if (cand_vx2.size() == 0)
	cand_vx2.insert(cand_vx2.end(), cand_vx.begin(), cand_vx.end());
    }
  else
    cand_vx2.insert(cand_vx2.end(), cand_vx.begin(), cand_vx.end());

  // Compute angles between in and outcoming edge in the candidate vertices
  // Use the angle internally to the face
  vector<double> angle(cand_vx2.size());
  size_t ki;
  for (ki=0; ki<cand_vx2.size(); ++ki)
    {
      vector<ftEdge*> edges = cand_vx2[ki]->getFaceEdges(face_.get());
      if (edges.size() < 2)
	{
	  angle[ki] = 0.0;  // Does not really make sense
	  continue;
	}

      // Don't expect more than two edges
      double t1 = edges[0]->parAtVertex(cand_vx2[ki].get());
      double t2 = edges[1]->parAtVertex(cand_vx2[ki].get());
      Point tan1 = edges[0]->tangent(t1);
      Point tan2 = edges[1]->tangent(t2);
//       if (edges[0]->tMax() - t1 < t1 - edges[0]->tMin())
// 	tan1 *= -1;
//       else 
// 	tan2 *= -1;
      if (edges[0]->tMax() - t1 > t1 - edges[0]->tMin())
	std::swap(tan1, tan2);
      tan1 *= -1;
      double ang = tan1.angle(tan2);

      Point par = cand_vx2[ki]->getFacePar(face_.get());
      Point norm1 = face_->normal(par[0], par[1]);
      Point norm2 = tan2.cross(tan1);
      if (norm1*norm2 < 0.0)
	ang = 2.0*M_PI - ang;
      angle[ki] = ang;
    }

  size_t max_idx = 0;
  double max_ang = angle[0];
  double prev_ang = angle[0];
  for (ki=1; ki<cand_vx2.size(); ++ki)
    {
      if (angle[ki] > max_ang)
	{
	  max_idx = ki;
	  prev_ang = max_ang;
	  max_ang = angle[ki];
	}
    }
      
  if (max_ang-prev_ang < 0.1*max_ang)
    std::cout << "Significant vertex insignificant" << std::endl;

  return cand_vx2[max_idx];
}

//==========================================================================
vector<shared_ptr<Vertex> >
RegularizeFace::prioritizeCornerVx(vector<shared_ptr<Vertex> > cand_vx)
//==========================================================================
{
  vector<shared_ptr<Vertex> > cand_vx2(cand_vx.begin(), cand_vx.end());
  int nmb_T_split = 0;
  if (cand_split_.size() > 0)
    {
      // Compute minimum distance between candidate vertices
      size_t ki, kj;
      double min_dist = 1.0e8;  // A huge number
      for (ki=0; ki<cand_vx.size(); ++ki)
	{
	  Point vxp1 = cand_vx[ki]->getVertexPoint();
	  for (kj=ki+1; kj<cand_vx.size(); ++kj)
	    min_dist = std::min(min_dist, 
				vxp1.dist(cand_vx[kj]->getVertexPoint()));
	}
      double level_dist = 0.2*min_dist;

      // Select vertices which corresponds to T-joints in the corresponding 
      // face
      vector<pair<Point,int> > corr_pts;
      for (ki=0; ki<cand_split_.size(); ++ki)
	{
	  for (kj=0; kj<corr_pts.size(); ++kj)
	    if (corr_pts[kj].first.dist(cand_split_[ki].first.first) < epsge_)
	      break;
	  if (kj < corr_pts.size())
	    corr_pts[kj].second++;
	  else
	    corr_pts.push_back(make_pair(cand_split_[ki].first.first, 1));

	  for (kj=0; kj<corr_pts.size(); ++kj)
	    if (corr_pts[kj].first.dist(cand_split_[ki].second.first) < epsge_)
	      break;
	  if (kj < corr_pts.size())
	    corr_pts[kj].second++;
	  else
	    corr_pts.push_back(make_pair(cand_split_[ki].second.first, 1));
	}

      for (ki=0; ki<corr_pts.size(); ++ki)
	{
	  if (corr_pts[ki].second < 3)
	    continue;   // Not a T-joint vertex

	  // Fetch corresponding point in this face
	  Point pos1;
	  double u1, v1, d1;
	  face_->closestPoint(corr_pts[ki].first, u1, v1, pos1, d1, epsge_);

	  // Find closest candidate vertex
	  double min_dist2 = min_dist;
	  size_t min_ind = 0;
	  for (kj=nmb_T_split; kj<cand_vx2.size(); ++kj)
	    {
	      double dist = pos1.dist(cand_vx2[kj]->getVertexPoint());
	      if (dist < min_dist2)
		{
		  min_dist2 = dist;
		  min_ind = kj;
		}
	    }
	  
	  if (min_dist2 < level_dist)
	    {
	      std::swap(cand_vx2[nmb_T_split], cand_vx2[min_ind]);
	      ++nmb_T_split;
	    }
	}
      if (cand_vx2.size() == 0)
	cand_vx2.insert(cand_vx2.end(), cand_vx.begin(), cand_vx.end());
    }

  // Compute angles between in and outcoming edge in the candidate vertices
  // Use the angle internally to the face
  vector<double> angle(cand_vx2.size());
  size_t ki;
  for (ki=0; ki<cand_vx2.size(); ++ki)
    {
#ifdef DEBUG_REG
      std::ofstream of("prio_vx.g2");
      of << "400 1 0 4 100 0 155 255" << std::endl;
      of << "1" << std::endl;
      of << cand_vx2[ki]->getVertexPoint() << std::endl;
#endif

      vector<ftEdge*> edges = cand_vx2[ki]->getFaceEdges(face_.get());
      if (edges.size() < 2)
	{
	  angle[ki] = 0.0;  // Does not really make sense
	  continue;
	}

      // Don't expect more than two edges
      double t1 = edges[0]->parAtVertex(cand_vx2[ki].get());
      double t2 = edges[1]->parAtVertex(cand_vx2[ki].get());
      Point tan1 = edges[0]->tangent(t1);
      Point tan2 = edges[1]->tangent(t2);
//       if (edges[0]->tMax() - t1 < t1 - edges[0]->tMin())
// 	tan1 *= -1;
//       else 
// 	tan2 *= -1;
      if (edges[0]->tMax() - t1 > t1 - edges[0]->tMin())
	std::swap(tan1, tan2);
      tan1 *= -1;
      double ang = tan1.angle(tan2);

      Point par = cand_vx2[ki]->getFacePar(face_.get());
      Point norm1 = face_->normal(par[0], par[1]);
      Point norm2 = tan2.cross(tan1);
      if (norm1*norm2 < 0.0)
	ang = 2.0*M_PI - ang;
      angle[ki] = ang;
    }

  // Sort vector according to decreasing angle, but prioritize T-joints
  int kr, kh;
  for (kr=0; kr<nmb_T_split; ++kr)
    for (kh=kr+1; kh<nmb_T_split; ++kh)
      if (angle[kh] > angle[kr])
	{
	  std::swap(cand_vx2[kr], cand_vx2[kh]);
	  std::swap(angle[kr], angle[kh]);
	}
      
  for (kr=nmb_T_split; kr<(int)cand_vx2.size(); ++kr)
    for (kh=kr+1; kh<(int)cand_vx2.size(); ++kh)
      if (angle[kh] > angle[kr])
	{
	  std::swap(cand_vx2[kr], cand_vx2[kh]);
	  std::swap(angle[kr], angle[kh]);
	}

  return cand_vx2;
}

//==========================================================================
void RegularizeFace::getConcaveCorners(vector<shared_ptr<Vertex> >& corners, 
				       vector<shared_ptr<Vertex> >& concave_corners)
//==========================================================================
{
  for (size_t ki=0; ki<corners.size(); ++ki)
    {
      vector<ftEdge*> edges = corners[ki]->getFaceEdges(face_.get());
      if (edges.size() < 2)
	continue;

      double t1 = edges[0]->parAtVertex(corners[ki].get());
      double t2 = edges[1]->parAtVertex(corners[ki].get());
      Point tan1 = edges[0]->tangent(t1);
      Point tan2 = edges[1]->tangent(t2);
      if (edges[0]->tMax() - t1 > t1 - edges[0]->tMin())
	std::swap(tan1, tan2);
      tan1 *= -1;
      double ang = tan1.angle(tan2);

      Point par = corners[ki]->getFacePar(face_.get());
      Point norm1 = face_->normal(par[0], par[1]);
      Point norm2 = tan2.cross(tan1);
      if (norm1*norm2 < 0.0)
	ang = 2.0*M_PI - ang;
      
      if (ang > M_PI)
	concave_corners.push_back(corners[ki]);
    }
}

//==========================================================================
bool RegularizeFace::getVertexProperties(shared_ptr<Vertex> vx, Point& parpnt,
					 double& ang, bool& T_joint)
//==========================================================================
{
  vector<ftEdge*> edges = vx->getFaceEdges(face_.get());
  if (edges.size() < 2)
      ang = 0.0;  // Does not really make sense

  // Don't expect more than two edges
  double t1 = edges[0]->parAtVertex(vx.get());
  double t2 = edges[1]->parAtVertex(vx.get());
  Point tan1 = edges[0]->tangent(t1);
  Point tan2 = edges[1]->tangent(t2);
  if (edges[0]->tMax() - t1 > t1 - edges[0]->tMin())
    std::swap(tan1, tan2);
  tan1 *= -1;
  ang = tan1.angle(tan2);

  Point par = vx->getFacePar(face_.get());
  Point norm1 = face_->normal(par[0], par[1]);
  Point norm2 = tan2.cross(tan1);
  if (norm1*norm2 < 0.0)
    ang = 2.0*M_PI - ang;

  int nmb_edges = vx->nmbUniqueEdges();
  T_joint = (nmb_edges > 2);

  // Check if the chord from the vertex to the given parameter position
  // lies inside the surface domain
  Point vx_par = vx->getFacePar(face_.get());
  shared_ptr<ParamSurface> sf = face_->surface();
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf);
  if (bd_sf.get())
    {
      CurveBoundedDomain dom = bd_sf->parameterDomain();
      shared_ptr<SplineCurve> cv = shared_ptr<SplineCurve>(new SplineCurve(parpnt, 
									   vx_par));
      vector<double> paramint;
      dom.findPcurveInsideSegments(*cv, epsge_, paramint);
      if (paramint.size() == 1 && paramint[0] > cv->endparam()-epsge_)
	return true;
      else
	return false;
    }
  else
    return true;
 }

//==========================================================================
vector<vector<ftEdge*> > RegularizeFace::getHalfHoles(int idx)
//==========================================================================
{
  vector<vector<ftEdge*> > half_holes;

  if (face_->nmbBoundaryLoops() == 0)
    return half_holes;  // Double closed surface, no outer boundary

  // Info about face curvature
  shared_ptr<ParamSurface> surf = face_->surface();
  DirectionCone cone1 = surf->tangentCone(true);
  DirectionCone cone2 = surf->tangentCone(false);
  
  // Fetch outer loop
  shared_ptr<Loop> loop = face_->getBoundaryLoop(idx);

  // Get edges
  vector<shared_ptr<ftEdgeBase> > edges = loop->getEdges();

  size_t ki, kj;
  vector<int> corner_ind;
  for (ki=0; ki<edges.size(); ++ki)
    {
      Point tan1 = edges[ki]->tangent(edges[ki]->tMax());
      Point tan2 = edges[ki]->next()->tangent(edges[ki]->next()->tMin());
      double ang = tan1.angle(tan2);
      if (ang > bend_)
	  corner_ind.push_back((int)ki);
    }

  if (corner_ind.size() <= 5 && face_->nmbBoundaryLoops() < 2)
    return half_holes;  // No possibility for half edges


  // Compute "opening angle" of candidate half hole
  vector<double> cone_ang;
  vector<pair<double,double> > opening_ang;
  vector<pair<double,double> > axis_ang;
  for (ki=0; ki<corner_ind.size(); ++ki)
    {
      kj = ki+1;
      if (kj == corner_ind.size())
	kj = 0;

      ftEdge *e1 = edges[corner_ind[ki]]->geomEdge();
      ftEdge *e2 = edges[corner_ind[kj]]->next()->geomEdge();
      ftEdge *e3 = edges[corner_ind[ki]]->next()->geomEdge();
      ftEdge *e4 = edges[corner_ind[kj]]->geomEdge();
      Point pos1 = e1->point(e1->tMax());
      Point tan1 = e1->tangent(e1->tMax());
      Point pos2 = e2->point(e2->tMin());
      Point tan2 = e2->tangent(e2->tMin());

      Point tan3 = e3->tangent(e3->tMin());
      Point tan4 = e4->tangent(e4->tMax());

      // Make curve bypassing the candidate half edge
      shared_ptr<SplineCurve> arc = 
	shared_ptr<SplineCurve>(new SplineCurve());
      HermiteInterpolator interpolator;
      vector<double> parval(4);
      vector<double> data;
      data.reserve(4*pos1.dimension());
      double dist = pos1.dist(pos2);
      parval[0] = parval[1] = 0.0;
      //parval[2] = parval[3] = dist;
      parval[2] = parval[3] = 1.0;
      tan1.normalize();
      tan1 *= 0.33*dist;
      tan2.normalize();
      tan2 *= 0.33*dist;
      data.insert(data.end(), pos1.begin(), pos1.end());
      data.insert(data.end(), tan1.begin(), tan1.end());
      data.insert(data.end(), pos2.begin(), pos2.end());
      data.insert(data.end(), tan2.begin(), tan2.end());
      arc->interpolate(interpolator, 4, pos1.dimension(), &parval[0],
		       &data[0]);
      DirectionCone cone = arc->directionCone();
      double angle;
      if (cone.greaterThanPi())
	angle = 1.5*M_PI;
      else
	  angle = cone.angle();
      cone_ang.push_back(angle);

#ifdef DEBUG_REG
      std::ofstream of0("arc.g2");
      arc->writeStandardHeader(of0);
      arc->write(of0);
#endif

      Point chord = pos2 - pos1;
      double ang3 = chord.angle(tan3);
      double ang4 = chord.angle(tan4);
      // if (tan3*tan4 < 0.0)
      // 	{
      // 	  // Adjacent tangents point in opposite directions. Not a
      // 	  // likely half hole
      // 	  opening_ang.push_back(make_pair(0.0, 0.0));
      // 	}
      // else
	opening_ang.push_back(make_pair(ang3, ang4));
	axis_ang.push_back(make_pair(cone.centre().angle(cone1.centre()),
				     cone.centre().angle(cone2.centre())));
    }

  double ang_limit = M_PI/2.0 - angtol_; //M_PI/4.0;
  double opening_limit = M_PI/4.0; 
  double fac = 0.5;
  // int nmb = (idx == 0) ? (int)corner_ind.size() - 4 :
  //   (int)corner_ind.size();  // Max number of half holes
  // int kr;
  // for (kr=0; kr<nmb; kr+=2)
  //   {
      // Find the candidate with the smallest opening angle
      //size_t min_idx;
      //double min_ang = 2*M_PI;
      for (ki=0; ki<cone_ang.size(); ++ki)
	{
	  double edge_ang = 0.5*(opening_ang[ki].first+opening_ang[ki].second);
      // 	  if (cone_ang[ki] < min_ang && edge_ang > opening_limit &&
      // 	      /* VSK 1208 I am not sure about this, but try it to avoid false
      // 		 half holes */ fabs(M_PI-edge_ang) > opening_limit)
      // 	    {
      // 	      min_idx = ki;
      // 	      min_ang = cone_ang[ki];
      // 	    }
      // 	}

      // if (min_ang < ang_limit)
      // 	{
      	  if (cone_ang[ki] < ang_limit && edge_ang > opening_limit &&
      	      /* VSK 1208 I am not sure about this, but try it to avoid false
      		 half holes */ fabs(M_PI-edge_ang) > 0.5*opening_limit)
      	    {
	      // Check if the candiate half hole may be the surface itself
	      if ((cone1.angle() >= ang_limit || cone1.greaterThanPi()) &&
		  axis_ang[ki].first > fac*opening_limit)
		continue;
	      if ((cone2.angle() >= ang_limit || cone2.greaterThanPi()) &&
		  axis_ang[ki].second > fac*opening_limit)
		continue;

	      // A half hole is found
	      vector<ftEdge*> curr_hole;
	      //ftEdge* curr_edge = edges[corner_ind[min_idx]]->next()->geomEdge();
	      ftEdge* curr_edge = edges[corner_ind[ki]]->next()->geomEdge();
	      curr_hole.push_back(curr_edge);
	      //kj = min_idx + 1;
	      kj = ki + 1;
	      if (kj == corner_ind.size())
		kj = 0;
	      while (true)
		{
		  if (curr_edge == edges[corner_ind[kj]].get())
		    break;
		  curr_edge = curr_edge->next()->geomEdge();
		  curr_hole.push_back(curr_edge);
		}
	      half_holes.push_back(curr_hole);
	  
	      // // Mark the half hole as used
	      // cone_ang[min_idx] = 2*M_PI;
	      // cone_ang[min_idx+1] = 2*M_PI;
	    }
	  //   else
	  // 	break;  // No more half holes
	  // }
	}

#ifdef DEBUG_REG
      std::ofstream of("half_holes.g2");
      for (size_t k2=0; k2<half_holes.size(); ++k2)
	{
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << 2*half_holes[k2].size() << std::endl;
	  for (size_t k3=0; k3<half_holes[k2].size(); ++k3)
	    {
	      of << half_holes[k2][k3]->getVertex(true)->getVertexPoint() << std::endl;
	      of << half_holes[k2][k3]->getVertex(false)->getVertexPoint() << std::endl;
	    }
	}
#endif

  return half_holes;
}

//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeFace::initIsolateHole(vector<vector<ftEdge*> >& half_holes)
//==========================================================================
{
  vector<shared_ptr<ftSurface> > dummy_vec;

  // First check if the face is too curved to behave well in the algorithm
  // in isolateHole. In that case, make an initial split
  // Pick the minium undarlying surface
  shared_ptr<ParamSurface> surf = face_->surface();
  shared_ptr<BoundedSurface> bd_surf =
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
  if (!bd_surf.get())
    return dummy_vec;

  RectDomain domain = bd_surf->containingDomain();
  shared_ptr<ParamSurface> under = bd_surf->underlyingSurface();
  SplineSurface* srf = under->asSplineSurface();

  if (!srf)
    return dummy_vec;
  
  RectDomain dom2 = srf->containingDomain();  // To avoid problems due to numerics
  double umin = std::max(domain.umin(), dom2.umin());
  double umax = std::min(domain.umax(), dom2.umax());
  double vmin = std::max(domain.vmin(), dom2.vmin());
  double vmax = std::min(domain.vmax(), dom2.vmax());

  shared_ptr<SplineSurface> srf2 = 
    shared_ptr<SplineSurface>(srf->subSurface(umin, vmin, umax, vmax));

  // Make parameter box around hole
  vector<ftEdge*> edges0;
  if (half_holes.size() == 1)
    edges0 = half_holes[0];
  else
    {
      shared_ptr<Loop> loop = face_->getBoundaryLoop(1);  // Inner loop
      size_t nmb_edges = loop->size();
      edges0.resize(nmb_edges);
      for (size_t ki=0; ki<nmb_edges; ++ki)
	edges0[ki] = loop->getEdge(ki)->geomEdge();
    }

  vector<Point> parpt0(edges0.size());
  for (size_t ki = 0; ki<edges0.size(); ++ki)
    {
      ftEdge *curr = edges0[ki]->geomEdge();
      parpt0[ki] = curr->faceParameter(curr->tMin());
    }
  BoundingBox parbox0;
  parbox0.setFromPoints(parpt0);
  Point low0 = parbox0.low();
  Point high0 = parbox0.high();

  // Compute the points in the face corresponding to the corners of
  // the parameter box
  vector<Point> holebox(4);
  holebox[0] = face_->point(low0[0],low0[1]);
  holebox[1] = face_->point(high0[0],low0[1]);
  holebox[2] = face_->point(low0[0],high0[1]);
  holebox[3] = face_->point(high0[0],high0[1]);

  
   // Compute tangent cones
  DirectionCone cone1 = srf2->tangentCone(true);
  DirectionCone cone2 = srf2->tangentCone(false);
  double fac = isolate_fac_;
  int gtpi1 = cone1.greaterThanPi();
  int gtpi2 = cone2.greaterThanPi();
  if (cone1.angle() < fac*M_PI && gtpi1 <= 0 && 
      cone2.angle() < fac*M_PI && gtpi2 <= 0)
    {
      // The surface is not too curved to be handled by the isolation algorithm
      return dummy_vec;
    }

  int dir = (cone1.angle() >= cone2.angle()) ? 0 : 1;  // Parameter direction of split

  // Make parameter box around hole
  vector<shared_ptr<ftEdgeBase> > edges = face_->getBoundaryLoop(1)->getEdges();
  vector<Point> parpt(4*edges.size());
  for (size_t ki = 0; ki<edges.size(); ++ki)
    {
      ftEdge *curr = edges[ki]->geomEdge();
      double t1 = curr->tMin();
      double t2 = curr->tMax();
      double tdel = (t2 - t1)/(double)4;
      for (int kr=0; kr<4; ++kr, t1+=tdel)
	parpt[4*ki+kr] = curr->faceParameter(t1);
    }
  BoundingBox parbox;
  parbox.setFromPoints(parpt);
  Point low = parbox.low();
  Point high = parbox.high();

  // Define split parameter
  double splitpar;
  double fac1 = 0.25;
  double fac2 = 0.75;
  if (srf2->endparam(dir) - high[dir] > low[dir] - srf2->startparam(dir))
    splitpar = std::min(high[dir] + fac1*(high[dir]-low[dir]), 
			fac1*high[dir]+fac2*srf2->endparam(dir));
  else
    splitpar = std::max(low[dir] - fac1*(high[dir]-low[dir]), 
			fac1*low[dir]+fac2*srf2->startparam(dir));

  Point par1, par2;
  if (dir == 0)
    {
      par1 = Point(splitpar, srf2->startparam(1));
      par2 = Point(splitpar, srf2->endparam(1));
    }
  else
    {
      par1 = Point(srf2->startparam(0), splitpar);
      par2 = Point(srf2->endparam(0), splitpar);
    }

  // Find closest point in domain
  CurveBoundedDomain dom = bd_surf->parameterDomain();
  Array<double,2> p1(par1[0],par1[1]);
  Array<double,2> p2(par2[0],par2[1]);
  Array<double,2> c1, c2;
  dom.closestOnBoundary(p1, c1, epsge_);
  dom.closestOnBoundary(p2, c2, epsge_);
  Point parpt1(c1[0],c1[1]);
  Point parpt2(c2[0],c2[1]);

  // Get candidate split vertices (outer loop)
  double maxdist = 0.5*low.dist(high);
  vector<shared_ptr<Vertex> > vxs = face_->getNonCornerVertices(bend_, 0);
  removeInsignificantVertices(vxs);
  double d1min = 1.0e8, d2min = 1.0e8;
  Point tmp1, tmp2;
  for (size_t ki=0; ki<vxs.size(); ++ki)
    {
      Point vxpar = vxs[ki]->getFacePar(face_.get());
      double d1 = vxpar.dist(parpt1);
      double d2 = vxpar.dist(parpt2);
      if (d1 < d1min)
	{
	  d1min = d1;
	  tmp1 = vxpar;
	}
      if (d2 < d2min)
	{
	  d2min = d2;
	  tmp2 = vxpar;
	}
    }

  if (d1min < maxdist && 
      !(d2min < d1min && tmp1.dist(tmp2) < 0.1*epsge_))
    parpt1 = tmp1;

  if (d2min < maxdist && 
      !(d1min < d2min && tmp1.dist(tmp2) < 0.1*epsge_))
    parpt2 = tmp2;

  // Perform split
  vector<shared_ptr<CurveOnSurface> > trim_segments;
  trim_segments = BoundedUtils::getTrimCrvsParam(surf, parpt1,
						 parpt2, epsge_,
						 bd_surf);
#ifdef DEBUG_REG
  std::ofstream out_file("split_segments.g2");
  for (size_t kj=0; kj<trim_segments.size(); ++kj)
    {
      shared_ptr<ParamCurve> cv = trim_segments[kj]->spaceCurve();
      cv->writeStandardHeader(out_file);
      cv->write(out_file);
    }
#endif

  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_surf, trim_segments, epsge_);

  // Fetch info about vertices not belonging to the corners of the
  // initial face
  vector<shared_ptr<Vertex> > non_corner = 
    face_->getNonCornerVertices(bend_);
  removeInsignificantVertices(non_corner);

  vector<shared_ptr<ftSurface> > faces = 
    RegularizeUtils::createFaces(sub_sfs, face_, epsge_,  tol2_, angtol_, 
				 non_corner);

  return faces;
}


//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeFace::isolateHole(const Point& seg_pnt, shared_ptr<Vertex> vx,
			    vector<vector<ftEdge*> >& half_holes)
//==========================================================================
{
  vector<shared_ptr<ftSurface> > faces;

  // Find the edges adjacent to the current vertex
  vector<ftEdge*> edges = vx->getFaceEdges(face_.get());

  // INFO
  vector<shared_ptr<Vertex> > non_corner = face_->getNonCornerVertices(bend_,
								       0);
#ifdef DEBUG_REG
  std::ofstream ofvx("cand_vx.g2");
  ofvx << "400 1 0 4 155 100 0 255" << std::endl;
  ofvx << non_corner.size() << std::endl;
  for (size_t kh=0; kh<non_corner.size(); ++kh)
    {
      Point vx_pt = non_corner[kh]->getVertexPoint();
      ofvx << vx_pt << std::endl;
    }
#endif

  // For each edges project the given point onto that edge
  // Select the edge with the closest projected point
  double bd_dist = 1.0e5*radius_;
  ftEdge *curr_edge = 0;
  ftEdge *found_edge = 0;
  shared_ptr<Vertex> curr_vx;
  double bd_par = 0.0;
  size_t ki;
  for (ki=0; ki<edges.size(); ++ki)
    {
      double par, dist;
      Point pnt;
      edges[ki]->closestPoint(seg_pnt, par, pnt, dist);
      if (dist < bd_dist)
	{
	  found_edge = edges[ki];
	  bd_dist = dist;
	  bd_par = par;
	}

      curr_vx = vx;
      curr_edge = edges[ki];
      while (true)
	{
	  shared_ptr<Vertex> other = curr_edge->getOtherVertex(curr_vx.get());
	  vector<shared_ptr<Vertex> >::iterator vxp =
	    std::find(corners_.begin(), corners_.end(), other);	  
	  if (vxp != corners_.end())
	    break;

	  ftEdge *next_edge;
	  double tp = curr_edge->parAtVertex(curr_vx.get());
	  if (tp - curr_edge->tMin() < curr_edge->tMax() - tp)
	    next_edge = curr_edge->next()->geomEdge();
	  else
	    next_edge = curr_edge->prev()->geomEdge();

	  next_edge->closestPoint(seg_pnt, par, pnt, dist);
	  if (dist < bd_dist)
	    {
	      found_edge = next_edge;
	      bd_dist = dist;
	      bd_par = par;
	    }
	  curr_edge = next_edge;
	  curr_vx = other;
	}

	  
    }

  // Search for a vertex nearby the boundary point
  shared_ptr<Vertex> select_vx;
  shared_ptr<Vertex> v1, v2;
  found_edge->getVertices(v1, v2);
  double t1 = found_edge->parAtVertex(v1.get());
  double t2 = found_edge->parAtVertex(v2.get());
  double len1 = found_edge->estimatedCurveLength(std::min(t1, bd_par),
						    std::max(t1, bd_par));
  double len2 = found_edge->estimatedCurveLength(std::min(t2, bd_par),
						    std::max(t2, bd_par));

  // Compare also distance to the hole centre
  Point bd_pnt = found_edge->point(bd_par);
  double bd_pnt_dist = bd_pnt.dist(centre_);
  double vx_dist1 = v1->getVertexPoint().dist(centre_);
  double vx_dist2 = v2->getVertexPoint().dist(centre_);
  
  double fac = 3.0;
  bool v1_ok = true;
  vector<shared_ptr<Vertex> > v1_end;
  vector<shared_ptr<Vertex> >::iterator vx1 =
    std::find(corners_.begin(), corners_.end(), v1);	  
  if (vx1 != corners_.end())
    v1_ok = false;
  else
    {
      // Search for candidate end vertices in a chain
      vector<shared_ptr<Vertex> > met_already;
      v1_end = RegularizeUtils::endVxInChain(face_, face_.get(), 
					     NULL, v1, v1, v1, met_already);
      for (ki=0; ki<v1_end.size(); )
	{
	  vector<shared_ptr<Vertex> >::iterator vx_it =
	    std::find(corners_.begin(), corners_.end(), v1_end[ki]);
	  vector<shared_ptr<Vertex> >::iterator vx_it2 =
	    std::find(non_corner.begin(), non_corner.end(), v1_end[ki]);
	  if (vx_it == corners_.end() && vx_it2 == non_corner.end())
	    v1_end.erase(v1_end.begin()+ki);
	  else
	    ki++;
	}
    }

  bool v2_ok = true;
  vector<shared_ptr<Vertex> > v2_end;
  vector<shared_ptr<Vertex> >::iterator vx2 =
    std::find(corners_.begin(), corners_.end(), v2);	  
  if (vx2 != corners_.end())
    v2_ok = false;
  else
    {
       // Search for candidate end vertices in a chain
      vector<shared_ptr<Vertex> > met_already;
      v2_end = RegularizeUtils::endVxInChain(face_, face_.get(), 
					     NULL, v2, v2, v2, met_already);
      for (ki=0; ki<v2_end.size(); )
	{
	  vector<shared_ptr<Vertex> >::iterator vx_it =
	    std::find(corners_.begin(), corners_.end(), v2_end[ki]);
	  vector<shared_ptr<Vertex> >::iterator vx_it2 =
	    std::find(non_corner.begin(), non_corner.end(), v2_end[ki]);
	  if (vx_it == corners_.end() && vx_it2 == non_corner.end())
	    v2_end.erase(v2_end.begin()+ki);
	  else
	    ki++;
	}
    }

  vector<shared_ptr<Vertex> > select_end;
  if (len1 <= len2 && len1 < fac*radius_ && 
      (vx_dist1 >= bd_pnt_dist || v1_end.size() > 0) &&
      v1.get() != vx.get() && v1_ok)
    {
      select_vx = v1;
      select_end = v1_end;
      bd_par = t1;
    }
  else if (len2 < fac*radius_ && 
	   (vx_dist2 >= bd_pnt_dist || v2_end.size() > 0) &&
	   v2.get() != vx.get() && v2_ok)
    {
      select_vx = v2;
      select_end = v2_end;
      bd_par = t2;
    }

  // Fetch all vertices belong to half holes and those next to it
  vector<shared_ptr<Vertex> > hole_vx;
  for (ki=0; ki<half_holes.size(); ++ki)
    {
      shared_ptr<Vertex> curr_vx = half_holes[ki][0]->getVertex(true);
      shared_ptr<Vertex> other_vx = 
	half_holes[ki][0]->prev()->geomEdge()->getOtherVertex(curr_vx.get());
      hole_vx.push_back(other_vx);
      hole_vx.push_back(curr_vx);
      for (size_t kj=0; kj<half_holes[ki].size(); ++kj)
	{
	  curr_vx = half_holes[ki][kj]->getVertex(false);
	  hole_vx.push_back(curr_vx);
	}
      other_vx = 
	half_holes[ki][0]->next()->geomEdge()->getOtherVertex(curr_vx.get());
      hole_vx.push_back(other_vx);
    }
  //removeInsignificantVertices(hole_vx);

  // Select also candidate vertices for splitting and eventual candidate edge
  // Perform splitting
  vector<shared_ptr<Vertex> > cand_vx;
  ftEdge* cand_edge = 0;
  vector<shared_ptr<Vertex> > vxs = face_->getBoundaryLoop(0)->getVertices();
  removeInsignificantVertices(vxs);
  if (select_vx.get())
    {
      // Avoid making 3-sided faces
      for (ki=0; ki<select_end.size(); )
	{
	  if (select_vx->sameEdge(select_end[ki].get()) ||
	      select_vx->connectedToSameVertex(select_end[ki].get()))
	    select_end.erase(select_end.begin()+ki);
	  else
	    ++ki;
	}

      bool strong = false;
      if (select_end.size() == 0)
	selectCandidateSplit(select_vx, vxs, cand_vx, cand_edge);
      else
	{
	  cand_vx = select_end;
	  strong = true;
	}

      // Remove vertices belonging to half holes and those next to it
      // from the candidate vertices 
      for (size_t kj=0; kj<hole_vx.size(); ++kj)
	{
	  vector<shared_ptr<Vertex> >::iterator vx =
	    std::find(cand_vx.begin(), cand_vx.end(), hole_vx[kj]);
	  if (vx != cand_vx.end())
	    cand_vx.erase(vx);
	}

      // Fetch info about vertices not belonging to the corners of the
      // initial face
      vector<shared_ptr<Vertex> > non_corner = 
	face_->getNonCornerVertices(bend_, 0);
      removeInsignificantVertices(non_corner);

      // Perform splitting
      // Unset axis information to avoid a bad choice of a splitting plane
      unsetAxis();
      cand_edge = 0;
      vector<shared_ptr<Vertex> > dummy_prio;
      faces = RegularizeUtils::divideVertex(face_, select_vx, 
					    cand_vx, cand_edge, dummy_prio,
					    epsge_, tol2_, angtol_, bend_,
					    non_corner, centre_, axis_,
					    strong);
    }
  else
    {
      selectCandidateSplit(found_edge, vxs, cand_vx, cand_edge);

      // Remove vertices belonging to half holes and those next to it
      // from the candidate vertices
      for (size_t kj=0; kj<hole_vx.size(); ++kj)
	{
	  vector<shared_ptr<Vertex> >::iterator vx =
	    std::find(cand_vx.begin(), cand_vx.end(), hole_vx[kj]);
	  if (vx != cand_vx.end())
	    cand_vx.erase(vx);
	}

      // Perform splitting
      cand_edge = 0;
      faces = divideVertex(found_edge, bd_par, cand_vx, cand_edge);
    }


  return faces;
}

//==========================================================================
void 
RegularizeFace::selectCandidateSplit(shared_ptr<Vertex> select_vx,
				      vector<shared_ptr<Vertex> >& vx,
				      vector<shared_ptr<Vertex> >& cand_vx,
				     ftEdge*& cand_edge, bool keep_T_joints)
//==========================================================================
{
  cand_edge = 0;
  for (size_t ki=0; ki<vx.size(); ++ki)
    {
      if (vx[ki].get() == select_vx.get())
	continue;   // Not a candiate
      ftEdge *common_edge = vx[ki]->getCommonEdgeInFace(select_vx.get(),
							face_.get());
      if (common_edge)
	continue;  // Already connected in this face
      vector<ftEdge*> edges = vx[ki]->getFaceEdges(face_.get());
      size_t kj;
      for (kj=0; kj<edges.size(); ++kj)
	{
	  shared_ptr<Vertex> other = 
	    edges[kj]->getOtherVertex(vx[ki].get());
	  if (other->sameEdge(select_vx.get()))
	    break;

	  if (select_vx->connectedToSameVertex(other.get()))
	    {
	      if (other->nmbUniqueEdges() == 2)
		{
		  vector<ftEdge*> edg2 = other->uniqueEdges();
		  double t1 = edg2[0]->parAtVertex(other.get());
		  double t2 = edg2[0]->parAtVertex(other.get());
		  Point tan1 = edg2[0]->tangent(t1);
		  Point tan2 = edg2[1]->tangent(t2);
		  if (tan1.angle(tan2) < bend_)
		    break;  // Not a significant vertex
		}

	      Vertex *common = select_vx->getCommonVertex(other.get());
	      if (common->nmbUniqueEdges() == 2)
		{
		  vector<ftEdge*> edg2 = common->uniqueEdges();
		  double t1 = edg2[0]->parAtVertex(common);
		  double t2 = edg2[0]->parAtVertex(common);
		  Point tan1 = edg2[0]->tangent(t1);
		  Point tan2 = edg2[1]->tangent(t2);
		  if (tan1.angle(tan2) < bend_)
		    break;  // Not a significant vertex
		}
	    }
	}	
      if (kj < edges.size())
	continue;
      
      cand_vx.push_back(vx[ki]);
    }
  removeInsignificantVertices(cand_vx, keep_T_joints); // Keep T-joints

  if (cand_vx.size() == 0)
    {
      // Find candiate edge
      vector<ftEdge*> edges = select_vx->getFaceEdges(face_.get());
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
void 
RegularizeFace::selectCandidateSplit(ftEdge* edge,
				      vector<shared_ptr<Vertex> >& vx,
				      vector<shared_ptr<Vertex> >& cand_vx,
				      ftEdge*& cand_edge)
//==========================================================================
{
  cand_edge = 0;
  ftEdge *e1 = edge->next()->geomEdge();
  ftEdge *e2 = edge->prev()->geomEdge();
  shared_ptr<Vertex> v1, v2, v3, v4;
  e1->getVertices(v1, v2);
  e2->getVertices(v3, v4);
  for (size_t ki=0; ki<vx.size(); ++ki)
    {
      if (vx[ki].get() == v1.get() || vx[ki].get() == v2.get())
	continue;   // Not a candiate
      if (vx[ki].get() == v3.get() || vx[ki].get() == v4.get())
	continue;   // Not a candiate

      cand_vx.push_back(vx[ki]);
    }
  removeInsignificantVertices(cand_vx);

  if (cand_vx.size() == 0)
    {
      // Find candiate edge
      cand_edge = edge->next()->next()->geomEdge();
    }

}

//==========================================================================
bool
RegularizeFace::sortAlongLine(vector<hole_info>& holes, Point& pnt,
			      Point& dir, vector<double>& parvals,
			      vector<int>& perm)
//==========================================================================
{
  // Find the most distant holes and make a line throuhg the hole centra
  // Fetch also the smallest and largest radius
  int max_ind1=-1, max_ind2=-1;
  double max_dist = 0.0;
  double min_rad = 1.0e6, max_rad = 0.0;
  int ki, kj;
  for (ki=0; ki<(int)holes.size(); ++ki)
    {
      min_rad = std::min(min_rad, holes[ki].hole_radius_);
      max_rad = std::max(max_rad, holes[ki].hole_radius_);

      for (kj=ki+1; kj<(int)holes.size(); ++kj)
	{
	  double dist = holes[ki].hole_centre_.dist(holes[kj].hole_centre_);
	  if (dist > max_dist)
	    {
	      max_ind1 = ki;
	      max_ind2 = kj;
	      max_dist = dist;
	    }
	}
    }
  if (max_ind1 < 0)
    return false;  // Less than two points

  pnt = holes[max_ind1].hole_centre_;
  dir = holes[max_ind2].hole_centre_ - holes[max_ind1].hole_centre_;
  dir.normalize();

  // Project the hole centra on to the line
  double level_rad1 = 2.0*max_rad;
  for (ki=0; ki<(int)holes.size(); ++ki)
    {
      double t1 = ( holes[ki].hole_centre_ - pnt)*dir;

      // Check distance
      Point tmp = pnt + t1*dir;
      double dist = tmp.dist(holes[ki].hole_centre_);
      if (dist > level_rad1)
	return false;  // Not a good case for sorting

      parvals.push_back(t1);
    }

  // Make permutation array
  perm.resize(holes.size());
  for (ki=0; ki<(int)holes.size(); ++ki)
    perm[ki] = ki;

  for (ki=0; ki<(int)holes.size(); ++ki)
    for (kj=ki+1; kj<(int)holes.size(); ++kj)
      if (parvals[perm[kj]] < parvals[perm[ki]])
	std::swap(perm[ki], perm[kj]);

  // Check distance between the points on the line corresponding
  // to the holes
  // @@@ VSK. 1113. After all, the holes are separate. Skips this test
  // until cases that should have been excluded by this test turns up
  //double fac = 2.0;
  double fac = 0.1; //0.5;
  for (ki=1; ki<(int)holes.size(); ++ki)
    {
      Point tmp1 = pnt + parvals[perm[ki-1]]*dir;
      Point tmp2 = pnt + parvals[perm[ki]]*dir;
      double dist = tmp1.dist(tmp2);
      if (dist < fac*(holes[perm[ki-1]].hole_radius_ +
  		      holes[perm[ki]].hole_radius_))
  	return false;
    }

  return true;  // Sorting performed
}

//==========================================================================
void
RegularizeFace::identifyParLine(vector<hole_info>& holes, Point& dir)
//==========================================================================
{
  // First project the hole centra into the parameter domain of the surface
  shared_ptr<ParamSurface> srf = face_->surface();
  shared_ptr<BoundedSurface> bd_srf = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(srf);
  if (bd_srf.get())
    srf = bd_srf->underlyingSurface();
  vector<Point> par_centre(holes.size());
  size_t ki;
  for (ki=0; ki<holes.size(); ++ki)
    {
      double u1, v1, d1;
      Point pos;
      srf->closestPoint(holes[ki].hole_centre_, u1, v1,
			  pos, d1, epsge_);
      par_centre[ki] = Point(u1, v1);
    }

  // Find a line through these parameter points
  RectDomain dom = face_->surface()->containingDomain();

  // First the first and the last point in each parameter direction
  int ix1=-1, ix2=-1, ix3=-1, ix4=-1;
  double minp = dom.umax(), maxp=dom.umin();
  for (ki=0; ki<par_centre.size(); ++ki)
    {
      if (par_centre[ki][0] < minp)
	{
	  minp = par_centre[ki][0];
	  ix1 = (int)ki;
	}
      if (par_centre[ki][0] > maxp)
	{
	  maxp = par_centre[ki][0];
	  ix2 = (int)ki;
	}
    }
  
  minp = dom.vmax();
  maxp=dom.vmin();
  for (ki=0; ki<par_centre.size(); ++ki)
    {
      if (par_centre[ki][1] < minp)
	{
	  minp = par_centre[ki][1];
	  ix3 = (int)ki;
	}
      if (par_centre[ki][1] > maxp)
	{
	  maxp = par_centre[ki][1];
	  ix4 = (int)ki;
	}
    }

  if ((dom.vmax() - dom.vmin())*(par_centre[ix2][0] - par_centre[ix1][0]) >
      (dom.umax() - dom.umin())*(par_centre[ix4][1] - par_centre[ix3][1]))
    dir = par_centre[ix2] - par_centre[ix1];
  else
    dir = par_centre[ix4] - par_centre[ix3];
 }

//==========================================================================
bool
RegularizeFace::sortRadially(vector<hole_info>& holes, const Point& wgt_pnt,
			     const Point& norm, vector<double>& angles, 
			     vector<int>& perm)
//==========================================================================
{
  if (holes.size() == 0)
    return false;  // No holes

  int nmb = (int)holes.size();
  angles.resize(nmb);

  // Make initial chord and project into the weight point plane
  int ki, kj;
  Point vec1 = wgt_pnt - holes[0].hole_centre_;
  vec1 -= (vec1*norm)*norm;
  Point dir = norm%vec1;
  angles[0] = 0.0;

  for (ki=1; ki<nmb; ++ki)
    {
      Point vec2 = wgt_pnt - holes[ki].hole_centre_;
      vec2 -= (vec2*norm)*norm;
      double ang = vec1.angle(vec2);
      if (dir*vec2 < 0.0)
	ang = 2.0*M_PI - ang;
      angles[ki] = ang;
    }
  // Make permutation array
  perm.resize(nmb);
  for (ki=0; ki<nmb; ++ki)
    perm[ki] = ki;

  for (ki=0; ki<nmb; ++ki)
    for (kj=ki+1; kj<nmb; ++kj)
      if (angles[perm[kj]] < angles[perm[ki]])
	std::swap(perm[ki], perm[kj]);

  // Make sure that the angles are distinct
  double tol = 0.01;
  for (ki=1; ki<nmb; ++ki)
    if (!(angles[perm[ki]] > angles[perm[ki-1]] + tol))
      return false;

  // Shift the holes such that the largest distance is between the first 
  // and the last hole
  double max_dist = holes[perm[0]].hole_centre_.dist(holes[perm[nmb-1]].hole_centre_);
  int ix = nmb;
  for (ki=0; ki<nmb-1; ++ki)
    {
      double dist = 
	holes[perm[ki]].hole_centre_.dist(holes[perm[ki+1]].hole_centre_);
      if (dist > max_dist)
	{
	  max_dist = dist;
	  ix = ki;
	}
    }
  if (ix < nmb)
    {
      // Shift holes
      perm.insert(perm.end(), perm.begin(), perm.begin()+ix+1);
      perm.erase(perm.begin(), perm.begin()+ix+1);
    } 

  return true;
}

//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::divideAcrossLine(vector<vector<ftEdge*> >& half_holes,
				 vector<hole_info>& holes, Point& pnt,
				 Point& dir, vector<int>& perm)

//==========================================================================
{
  // For each pair of consequetive holes, make one or two splitting curves
  // between the holes to separate them
  // In a later version, we may also consider the distance to the outer
  // boundary
  shared_ptr<ParamSurface> surf = face_->surface();
  vector<Point> split_pnt;
  vector<Point> split_norm;
  size_t ki;
  double min_rad = holes[0].hole_radius_;
  vector<double> min_dist;
  double len_fac = 0.2;
  vector<pair<int, int> > hole_idx;
  vector<double> seg_lengths;
  vector<pair<Point, Point> > seg_endpt;
  double len;
  for (ki=1; ki<holes.size(); ++ki)
    {
      // Make a plane following a chord between holes midpoints,
      // approximately being perpendicular to the hole axis
      Point tmp1 = holes[perm[ki-1]].hole_centre_;
      Point tmp2 = holes[perm[ki]].hole_centre_;
      Point mid = 0.5*(tmp1 + tmp2);
      Point axis = 0.5*(holes[perm[ki-1]].hole_axis_ +
			holes[perm[ki]].hole_axis_);
      Point normal = dir.cross(axis);
      normal = (tmp2 - tmp1).cross(axis);
      normal.normalize();

      // Compute intersections between the face and this plane and pick
      // the appropriate piece(s)
      shared_ptr<BoundedSurface> bd_surf;
      vector<shared_ptr<CurveOnSurface> > trim_seg =
	BoundedUtils::getPlaneIntersections(surf, mid, normal, 
					    epsge_, bd_surf);

      // Remove the segments lying outside the area between the two holes
      size_t kj;
      for (kj=0; kj<trim_seg.size(); )
	{
	  double tmp_par = 0.5*(trim_seg[kj]->startparam() +
				trim_seg[kj]->endparam());
	  Point tmp3 = trim_seg[kj]->ParamCurve::point(tmp_par);
	  if ((tmp3 - tmp1)*(tmp3 - tmp2) >= 0.0)
	    trim_seg.erase(trim_seg.begin()+kj);
	  else
	    kj++;
	}
      
#ifdef DEBUG_REG
  std::ofstream out_file_1("trim_segments.g2");
  for (kj=0; kj<trim_seg.size(); ++kj)
    {
      shared_ptr<ParamCurve> cv = trim_seg[kj]->spaceCurve();
      cv->writeStandardHeader(out_file_1);
      cv->write(out_file_1);
    }
#endif
       // Estimate length of remaining segments
      len = 0.0;
      Point p1 = tmp2, p2 = tmp1;
      for (kj=0; kj<trim_seg.size(); ++kj)
	{
	  len += trim_seg[kj]->estimatedCurveLength();
	  Point pos1 = 
	    trim_seg[kj]->ParamCurve::point(trim_seg[kj]->startparam());
	  Point pos2 = 
	    trim_seg[kj]->ParamCurve::point(trim_seg[kj]->endparam());
	  if (tmp1.dist(pos1) < tmp1.dist(p1))
	    p1 = pos1;
	  if (tmp1.dist(pos2) < tmp1.dist(p1))
	    p1 = pos2;
	  if (tmp2.dist(pos1) < tmp2.dist(p2))
	    p2 = pos1;
	  if (tmp2.dist(pos2) < tmp2.dist(p2))
	    p2 = pos2;
	}

      // Define split points
      double r1 = holes[perm[ki-1]].hole_radius_;
      double r2 = holes[perm[ki]].hole_radius_;
      min_rad = std::min(min_rad, r2);

      int nmb_cand_split = 0;
      if (cand_split_.size() > 0)
	nmb_cand_split = nmbSplitPattern(p1, p2);
	
      double max_edge_len = 3.0*(r1 + r2); //6.0*(r1 + r2); //4.0*(r1 + r2);
      //double min_edge_len = 0.2*std::min(r1, r2);
      double av_rad = 0.5*(r1 + r2);
      vector<Point> seg_pnt;
      // if (len < min_edge_len && nmb_cand_split == 0)
      //if (len < 2.0*av_rad && split_mode_ > 1)
      if (len < 10.0*av_rad && split_mode_ > 1)
      	{
      	  // Perform alternative split
	  hole_idx.push_back(make_pair(perm[ki-1], perm[ki]));
	  seg_lengths.push_back(len);
	  seg_endpt.push_back(make_pair(p1, p2));
      	}
      else
	{
	  size_t curr_seg_size = seg_pnt.size();
	  if ((len > max_edge_len && nmb_cand_split != 1) || 
	      nmb_cand_split > 1)
	    {
	      // Define two split points
	      double d1 = p1.dist(p2);
	      if (d1 < 2.0*(r1 + r2))
		d1 = 2.0*(r1 + r2);
	      double t1 = std::min(0.33, 2.0*r1/d1);
	      Point tmp0 = (1.0-t1)*p1 + t1*p2;
	      seg_pnt.push_back(tmp0);
	      double t2 = std::min(0.33, 2.0*r2/d1);
	      tmp0 = t2*p1 + (1.0-t2)*p2;
	      seg_pnt.push_back(tmp0);
	    }
	  else
	    {
	      // One split point
	      // double t1 = r1/(r1 + r2);
	      double t1 = r2/(r1 + r2);
	      Point tmp0 = t1*p1 + (1.0-t1)*p2;
	      seg_pnt.push_back(tmp0);
	    }

	  // Fetch split point lying on the trim segments
	  for (kj=curr_seg_size; kj<seg_pnt.size(); kj++)
	    {
	      double curr_par;
	      double init_dist = 2.0*len;
	      int curr_idx = -1;
	      for (size_t kr=0; kr<trim_seg.size(); ++kr)
		{
		  // Find closest point on segment to initial point
		  double par, dist;
		  Point seg_pnt2;
		  trim_seg[kr]->closestPoint(seg_pnt[kj], 
					     trim_seg[kr]->startparam(),
					     trim_seg[kr]->endparam(),
					     par, seg_pnt2, dist);
		  if (dist < init_dist)
		    {
		      init_dist= dist;
		      curr_par = par;
		      curr_idx = (int)kr;
		    }
		}

	      if (curr_idx >= 0)
		{
		  // Evaluate point and tangent corresponding to the trimming curve
		  vector<Point> der(2);
		  trim_seg[curr_idx]->point(der, curr_par, 1);
		  split_pnt.push_back(der[0]);
		  der[1].normalize();
		  split_norm.push_back(der[1]);
		  min_dist.push_back(std::max(min_rad, len_fac*len));
		}
	    }
	}
    }

  if (split_mode_ == 3)
    {
      // Check if splitting from the first and/or last hole to the consequtive
      // boundary should be performed
      int kj;
      for (ki=0, kj=1; ki<holes.size(); 
	   ki+=(holes.size()-1), kj+=((int)holes.size()-3))
	{
	  // Construct the chord between holes midpoints,
	  Point mid = holes[perm[ki]].hole_centre_;
	  Point tmp = holes[perm[kj]].hole_centre_;
	  Point chord = tmp - mid;

	  // Compute the closest points from the hole centre to the outer
	  // boundaries of the face and select the one being close to the
	  // centre and most parallel to the chord
	  vector<shared_ptr<ftEdge> > bd_edges = face_->getAllEdges(0);
	  vector<Point> clo_pts;
	  vector<double> clo_dist;
	  for (size_t kr=0; kr<bd_edges.size(); ++kr)
	    {
	      double par, dist;
	      Point pos;
	      bd_edges[kr]->closestPoint(mid, par, pos, dist);
	      clo_pts.push_back(pos);
	      clo_dist.push_back(dist);
	    }

	  // Select closest point on edge
	  int min_idx = -1;
	  double len = HUGE;
	  double angle = HUGE;
	  for (size_t kr=0; kr<clo_pts.size(); ++kr)
	    {
	      Point vec = mid - clo_pts[kr];
	      double ang = chord.angle(vec);
	      ang = std::min(M_PI-ang, ang);
	      if ((clo_dist[kr] < 1.5*len && ang < angle) ||
		  (clo_dist[kr] < len && ang < angle + bend_))
		{
		  min_idx = (int)kr;
		  len = clo_dist[kr];
		  angle = ang;
		}
	    }

	  if (min_idx >= 0)
	    {
	      // Compute closest point on hole boundary
	      shared_ptr<Loop> loop = face_->getBoundaryLoop(perm[ki]+1);
	      double par, dist;
	      Point pos;
	      int ind;
	      loop->closestPoint(clo_pts[min_idx], ind, par, pos, dist);
 	      if (dist < 2.0*holes[perm[ki]].hole_radius_)
		{
		  
		  hole_idx.push_back(make_pair(perm[ki], -1));
		  seg_lengths.push_back(dist);
		  seg_endpt.push_back(make_pair(pos, clo_pts[min_idx]));
		}
	    }
	}
    }

  vector<shared_ptr<ftSurface> > dummy_faces;
  if (split_pnt.size() > 0)
    {
      // Divide face with respect to the specified planes, but make modifications
      // for vertices at the outer boundary.
#ifdef DEBUG_REG
      std::ofstream split_out("split_out.g2");
      for (size_t k5=0; k5<split_pnt.size(); ++k5)
	{
	  split_out << "400 1 0 4 0 255 0 255" << std::endl;
	  split_out << " 1 " << std::endl;
	  split_out << split_pnt[k5] << std::endl;
	}
#endif
	  
    return divideByPlanes(split_pnt, split_norm, half_holes, min_dist);
    }
  else if (hole_idx.size() > 0)
    {
      // Make two splits from hole to hole to combine them into one hole
      // instead of separating the holes
      vector<shared_ptr<ftSurface> > faces =
	holeToHoleSplit(half_holes, holes, hole_idx, seg_lengths, seg_endpt, -1);
      if (faces.size() > 0)
	return faces;
    }
  return dummy_faces;
}

//===========================================================================
bool compare_first(pair<double,double> f1, pair<double,double> f2)
{
  return (f1.first < f2.first);
}

//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::isolateHolesParametrically(vector<vector<ftEdge*> >& half_holes,
					   Point& dir)

//==========================================================================
{
  vector<shared_ptr<ftSurface> > faces;

  // Set parameter direciton in which to split
  int dir_idx = (fabs(dir[0]) > fabs(dir[1])) ? 0 : 1;

  // For each hole and half hole, make parameter boxes and extract the
  // relevant limiting parameters
  int nmb = face_->nmbBoundaryLoops();
  vector<pair<double,double> > limit(half_holes.size()+nmb-1);
  size_t ki, kr;
  for (ki=0, kr=0; ki<half_holes.size(); ++ki, ++kr)
    {
      vector<ftEdge*> edges = half_holes[ki];
      // Make parameter box around hole
      vector<Point> parpt(2*edges.size()+1);
      size_t kj = 0;
      ftEdge *curr;
      for (kj=0; kj<edges.size(); ++kj)
	{
	  curr = edges[kj]->geomEdge();
	  parpt[2*kj] = curr->faceParameter(curr->tMin());
	  parpt[2*kj+1] = curr->faceParameter(0.5*(curr->tMin()+curr->tMax()));
	}
	parpt[2*kj] = curr->faceParameter(curr->tMax());
	
	BoundingBox parbox;
	parbox.setFromPoints(parpt);
	Point low = parbox.low();
	Point high = parbox.high();
	
	limit[ki].first = low[dir_idx];
	limit[ki].second = high[dir_idx];
      }

  for (ki=1; ki<(size_t)nmb; ++ki, ++kr)
    {
	shared_ptr<Loop> loop = face_->getBoundaryLoop((int)ki);
      vector<shared_ptr<ftEdgeBase> > edges = loop->getEdges();
      // Make parameter box around hole
      vector<Point> parpt(2*edges.size()+1);
      size_t kj = 0;
      ftEdge *curr = edges[0]->geomEdge();
      for (kj=0; kj<edges.size(); ++kj)
	{
	  curr = edges[kj]->geomEdge();
	  parpt[2*kj] = curr->faceParameter(curr->tMin());
	  parpt[2*kj+1] = curr->faceParameter(0.5*(curr->tMin()+curr->tMax()));
	}
	parpt[2*kj] = curr->faceParameter(curr->tMax());
	
	BoundingBox parbox;
	parbox.setFromPoints(parpt);
	Point low = parbox.low();
	Point high = parbox.high();
	
	limit[kr].first = low[dir_idx];
	limit[kr].second = high[dir_idx];
      }

  if (limit.size() < 2)
    return faces;
	  
  // Sort parameter limits
  std::sort(limit.begin(), limit.end(), compare_first);
      
  // // Check overlaps
  // for (ki=1; ki<limit.size(); ++ki)
  //   if (limit[ki].first < limit[ki-1].second + tol2_)
  //     return faces;

  // Join parameter intervals in case of overlaps
  vector<pair<double, double> > max_size(limit.size());
  for (ki=1; ki<limit.size(); )
    {
      if (limit[ki].first < limit[ki-1].second + tol2_)
	{
	  max_size[ki-1] = make_pair(limit[ki-1].second-limit[ki-1].first,
				     limit[ki].second-limit[ki].first);
	  max_size.erase(max_size.begin()+ki);
	  limit[ki-1].first = std::min(limit[ki-1].first, limit[ki].first);
	  limit[ki-1].second = std::max(limit[ki-1].second, limit[ki].second);
	  limit.erase(limit.begin()+ki);
	}
      else
	{
	  ki++;
	  max_size[ki-1] = make_pair(limit[ki-1].second-limit[ki-1].first,
				     limit[ki-1].second-limit[ki-1].first);
	}
    }
  max_size[ki-1] = make_pair(limit[ki-1].second-limit[ki-1].first,
			     limit[ki-1].second-limit[ki-1].first);


  if (limit.size() > 3)
    {
      // Check if the distance between the parameter intervals
      // varies largely. Do not split between relatively close
      // intervals. First compute mean distance
      double av_dist = 0.0;
      double max_dist = 0.0;
      for (ki=1; ki<limit.size(); ++ki)
	{
	  av_dist += (limit[ki].first - limit[ki-1].second);
	  max_dist = std::max(max_dist, limit[ki].first - limit[ki-1].second);
	}
      av_dist /= (double)(limit.size()-1);

      // Merge too dense parameter intervals
      double frac = (av_dist > 0.5*max_dist) ? 0.9 : 0.5*max_dist/av_dist;
      for (ki=1; ki<limit.size(); )
	{
	  if (limit[ki].first - limit[ki-1].second < frac*av_dist)
	    {
	      max_size[ki-1] = make_pair(max_size[ki-1].first, max_size[ki].second);
	      max_size.erase(max_size.begin()+ki);
	      limit[ki-1].first = std::min(limit[ki-1].first, limit[ki].first);
	      limit[ki-1].second = std::max(limit[ki-1].second, limit[ki].second);
	      limit.erase(limit.begin()+ki);
	    }
	  else
	    ki++;
	}
    }
      

  shared_ptr<ParamSurface> surf = face_->surface();
  RectDomain dom = surf->containingDomain();

  // Compute split curves
  vector<shared_ptr<CurveOnSurface> > trim_segments;
  shared_ptr<BoundedSurface> bd_surf;
  for (ki=1; ki<limit.size(); ++ki)
    {
      // Check if one or two split curves are feasible
      double fac = 4.0;
      // int nmb_split =  (limit[ki].first - limit[ki-1].second >
      // 			fac*std::max(limit[ki-1].second-limit[ki-1].first,
      // 				     limit[ki].second-limit[ki].first)) ? 2 : 1;
      int nmb_split =  (limit[ki].first - limit[ki-1].second >
			fac*std::max(max_size[ki-1].second, max_size[ki].first)) ? 2 : 1;
      double par;
      double del = (limit[ki].first - limit[ki-1].second)/(double)(nmb_split+1);
      double del1 = (nmb_split == 1) ? del : std::min(del, 2.0*max_size[ki-1].second);
      double del2 = (nmb_split == 1) ? del : std::min(del, 2.0*max_size[ki].first);
      int kj;
      for (kj=0, par = limit[ki-1].second+del1; kj<nmb_split; ++kj, par=limit[ki].first-del2)
	{
	  Point parpt1(2), parpt2(2);
	  parpt1[dir_idx] = parpt2[dir_idx] = par;
	  parpt1[1-dir_idx] = (dir_idx == 0) ? dom.vmin() : dom.umin();
	  parpt2[1-dir_idx] = (dir_idx == 0) ? dom.vmax() : dom.umax();
	  vector<shared_ptr<CurveOnSurface> > segments;
	  segments = BoundedUtils::getTrimCrvsParam(surf, parpt1,
						    parpt2, epsge_,
						    bd_surf);
	  trim_segments.insert(trim_segments.end(), segments.begin(), segments.end());
	}
    }
  
#ifdef DEBUG_REG
  std::ofstream out_file_1("trim_segments.g2");
  for (size_t kv=0; kv<trim_segments.size(); ++kv)
    {
      shared_ptr<ParamCurve> cv = trim_segments[kv]->spaceCurve();
      cv->writeStandardHeader(out_file_1);
      cv->write(out_file_1);
    }
#endif
  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_surf, trim_segments, epsge_);

#ifdef DEBUG_REG
	  std::ofstream of("split_surf.g2");
	  for (size_t kr=0; kr<sub_sfs.size(); ++kr)
	    {
	      sub_sfs[kr]->writeStandardHeader(of);
	      sub_sfs[kr]->write(of);
	    }
#endif
      

  // Fetch info about vertices not belonging to the corners of the
  // initial face
  vector<shared_ptr<Vertex> > non_corner = 
    face_->getNonCornerVertices(bend_);
  removeInsignificantVertices(non_corner);

  faces = 
    RegularizeUtils::createFaces(sub_sfs, face_, epsge_,  tol2_, angtol_, 
				 non_corner);

  return faces;
}
//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::isolateHolesRadially(vector<vector<ftEdge*> >& half_holes,
				     const Point& mid, const Point& axis,
				     int loop_idx,
				     vector<hole_info>& holes, 
				     vector<int>& perm)
//==========================================================================
{
  // For each pair of consequetive holes, make one or two splitting curves
  // between the holes to separate them
  // In a later version, we may also consider the distance to the outer
  // boundary
  shared_ptr<ParamSurface> surf = face_->surface();
  vector<Point> seg_pnt;
  vector<Point> seg_norm;
  size_t ki, kj, kh, kr;
  double min_rad = holes[0].hole_radius_;
  size_t nmb_holes = (loop_idx < 0 && split_mode_ == 1) ? 1 : holes.size();
  vector<double> min_dist;
  vector<pair<int, int> > hole_idx;
  vector<double> seg_lengths;
  vector<pair<Point, Point> > seg_endpt;

  for (ki=0; ki<nmb_holes/*holes.size()*/; ++ki)
    {
      kh = ki+1;
      if (kh == holes.size())
	kh = 0;

      // Check if a cylinder or a plane intersection is appropriate
      double len_fac = 3.0; //6.0;
      double d1 = mid.dist(holes[perm[ki]].hole_centre_);
      double d2 = mid.dist(holes[perm[kh]].hole_centre_);
      double min_rad = std::min(holes[perm[ki]].hole_radius_,
				holes[perm[kh]].hole_radius_);
      double cylinder_fac = 0.5;
      shared_ptr<BoundedSurface> bd_surf;
      vector<shared_ptr<CurveOnSurface> > trim_seg;
      Point tmp1 = holes[perm[ki]].hole_centre_;
      Point tmp2 = holes[perm[kh]].hole_centre_;
      for (kr=0; kr<2; ++kr)
	{
	  // Construct a curve between the two holes to collect
	  // information for the splitting between the holes
	  // Try first a cylinder intersection
	  if (kr == 0 && fabs(d1 - d2) < cylinder_fac*min_rad)
	    {
	      // Intersect with cylinder
	      double cyl_rad = 0.5*(d1 + d2);
	      trim_seg = BoundedUtils::getCylinderIntersections(surf, mid, 
								axis, cyl_rad,
								epsge_, bd_surf);
	      //len_fac = 3.0;
	      
	    }
	  else if (kr == 0)
	    continue;
	  else
	    {
	      // Make a plane following a chord between holes midpoints,
	      // approximately being perpendicular to the hole axis
	      Point mid2 = 0.5*(tmp1 + tmp2);
	      Point dir = tmp2 - tmp1;
	      Point axis2 = 0.5*(holes[perm[ki]].hole_axis_ +
				 holes[perm[kh]].hole_axis_);
	      Point normal = dir.cross(axis2);
	      normal.normalize();
	      
	      // Compute intersections between the face and this plane and pick
	      // the appropriate piece(s)
	      trim_seg = BoundedUtils::getPlaneIntersections(surf, mid2, normal, 
							     epsge_, bd_surf);
	    }

#ifdef DEBUG_REG
	  std::ofstream out_file0("split_segments0.g2");
	  for (kj=0; kj<trim_seg.size(); ++kj)
	    {
	      shared_ptr<ParamCurve> cv = trim_seg[kj]->spaceCurve();
	      cv->writeStandardHeader(out_file0);
	      cv->write(out_file0);
	    }
#endif

	  // Remove the segments lying outside the area between the two holes
	  for (kj=0; kj<trim_seg.size(); )
	    {
	      double tmp_par = 0.5*(trim_seg[kj]->startparam() +
				    trim_seg[kj]->startparam());
	      Point tmp3 = trim_seg[kj]->ParamCurve::point(tmp_par);
	      if ((tmp3 - tmp1)*(tmp3 - tmp2) >= 0.0)
		trim_seg.erase(trim_seg.begin()+kj);
	      else
		kj++;
	    }
      
	  // TESTING
	  if (trim_seg.size() > 1 && (nmb_holes == 1 || kr == 0))
	    trim_seg.clear();

	  if (trim_seg.size() > 0)
	    break;
	}

      if (trim_seg.size() == 0 && nmb_holes == 1)
	{
	  // No valid information so far. Try another approach
	  getSegPntRadialSplit(holes[perm[ki]], holes[perm[kh]],
			       len_fac, seg_pnt, seg_norm, min_dist);
	}
      else if (trim_seg.size() == 0)
	continue;
      else
	{
	  // Estimate length of remaining segments
	  double len = 0.0;
	  Point p1 = tmp2, p2 = tmp1;
	  for (kj=0; kj<trim_seg.size(); ++kj)
	    {
	      len += trim_seg[kj]->estimatedCurveLength();
	      Point pos1 = 
		trim_seg[kj]->ParamCurve::point(trim_seg[kj]->startparam());
	      Point pos2 = 
		trim_seg[kj]->ParamCurve::point(trim_seg[kj]->endparam());
	      if (tmp1.dist(pos1) < tmp1.dist(p1))
		p1 = pos1;
	      if (tmp1.dist(pos2) < tmp1.dist(p1))
		p1 = pos2;
	      if (tmp2.dist(pos1) < tmp2.dist(p2))
		p2 = pos1;
	      if (tmp2.dist(pos2) < tmp2.dist(p2))
		p2 = pos2;
	    }
	  len = p1.dist(p2);

#ifdef DEBUG_REG
	  std::ofstream out_file("split_segments.g2");
	  for (kj=0; kj<trim_seg.size(); ++kj)
	    {
	      shared_ptr<ParamCurve> cv = trim_seg[kj]->spaceCurve();
	      cv->writeStandardHeader(out_file);
	      cv->write(out_file);
	    }
	  out_file << "400 1 0 4 255 0 0 255" << std::endl;
	  out_file << "2" << std::endl;
	  out_file << tmp1 << std::endl;
	  out_file << tmp2 << std::endl;
	  out_file << "400 1 0 4 0 255 0 255" << std::endl;
	  out_file << "2" << std::endl;
	  out_file << p1 << std::endl;
	  out_file << p2 << std::endl;
#endif

	  // Define split points
	  double r1 = tmp1.dist(p1); //holes[perm[ki]].hole_radius_;
	  double r2 = tmp2.dist(p2); //holes[perm[kj]].hole_radius_;
	  double max_edge_len = len_fac*(r1 + r2); //4.0*(r1 + r2);
	  //min_rad = std::min(min_rad, r2);
	  double dist_fac = 0.5; //0.25;
	  min_dist.push_back(dist_fac*p1.dist(p2));

	  int nmb_cand_split = 0;
	  if (cand_split_.size() > 0)
	    nmb_cand_split = nmbSplitPattern(p1, p2);

	  double fac = 1.2;
	  double av_rad = 0.5*(r1 + r2);
	  //if (len < 3.0*av_rad && split_mode_ > 1)
	  if (len < 10.0*av_rad && split_mode_ > 1)
	    {
	      // Perform alternative split
	      hole_idx.push_back(make_pair(perm[ki], perm[kh]));
	      seg_lengths.push_back(len);
	      seg_endpt.push_back(make_pair(p1, p2));
	    }
	  else if (trim_seg.size() > 1 && split_mode_ > 1)
	    continue;   // No splitting related to these holes
	  else if ((len > max_edge_len && nmb_cand_split != 1) || 
		   nmb_cand_split > 1)
	    {
	      // Define two split points
	      double d1 = p1.dist(p2);
	      if (d1 < 2.0*(r1 + r2))
		d1 = 2.0*(r1 + r2);
	      double t1 = fac*r1/d1;
	      Point tmp0 = (1.0-t1)*p1 + t1*p2;
	      seg_pnt.push_back(tmp0);
	      Point tmpnorm = p2 - p1;
	      tmpnorm.normalize();
	      seg_norm.push_back(tmpnorm);
	      double t2 = fac*r2/d1;
	      tmp0 = t2*p1 + (1.0-t2)*p2;
	      seg_pnt.push_back(tmp0);
	      seg_norm.push_back(tmpnorm);
	      min_dist[min_dist.size()-1] *= 0.5;
	      min_dist.push_back(min_dist[min_dist.size()-1]);
	    }
	  else
	    {
	      // One split point
	      //double t1 = r1/(r1 + r2);
	      double t1 = r2/(r1 + r2);
	      Point tmp0 = t1*p1 + (1.0-t1)*p2;
	      seg_pnt.push_back(tmp0);
	      Point tmpnorm = p2 - p1;
	      tmpnorm.normalize();
	      seg_norm.push_back(tmpnorm);
	    }
	}
    }

  if (split_mode_ == 3)
    {
      // Check if splitting from the first and/or last hole to the consequtive
      // boundary should be performed
      int kj;
      for (ki=0, kj=1; ki<holes.size(); 
	   ki+=(holes.size()-1), kj+=((int)holes.size()-3))
	{
	  // Construct the chord between holes midpoints,
	  Point mid = holes[perm[ki]].hole_centre_;
	  double rad = holes[perm[ki]].hole_radius_;
	  Point tmp = holes[perm[kj]].hole_centre_;
	  Point chord = tmp - mid;
	  chord.normalize();
	  Point mid2 = mid - rad*chord;

	  // Compute the closest points from the hole centre to the outer
	  // boundaries of the face and select the one being close to the
	  // centre and most parallel to the chord
	  vector<shared_ptr<ftEdge> > bd_edges = face_->getAllEdges(0);
	  vector<Point> clo_pts;
	  vector<double> clo_dist;
	  for (size_t kr=0; kr<bd_edges.size(); ++kr)
	    {
	      double par, dist;
	      Point pos;
	      bd_edges[kr]->closestPoint(mid2, par, pos, dist);
	      clo_pts.push_back(pos);
	      clo_dist.push_back(dist);
	    }

	  // Select closest point on edge
	  int min_idx = -1;
	  double len = HUGE;
	  double angle = HUGE;
	  for (size_t kr=0; kr<clo_pts.size(); ++kr)
	    {
	      Point vec = mid - clo_pts[kr];
	      double ang = chord.angle(vec);
	      ang = std::min(M_PI-ang, ang);
	      if ((clo_dist[kr] < 1.5*len && ang < angle) ||
		  (clo_dist[kr] < len && ang < angle + bend_))
		{
		  min_idx = (int)kr;
		  len = clo_dist[kr];
		  angle = ang;
		}
	    }

	  if (min_idx >= 0)
	    {
	      // Project point onto chord and find a new position at the boundary
	      Point pq = mid + ((clo_pts[min_idx] - mid)*chord)*chord;
	      shared_ptr<Loop> loop = face_->getBoundaryLoop(0);
	      double par, dist;
	      Point pos;
	      int ind;
	      loop->closestPoint(pq, ind, par, pos, dist);
	      pq = pos;
 
	      // Compute closest point on hole boundary
	      loop = face_->getBoundaryLoop(perm[ki]+1);
	      loop->closestPoint(pq, ind, par, pos, dist);
 	      if (dist < 3.0*holes[perm[ki]].hole_radius_)
		{
		  
		  hole_idx.push_back(make_pair(perm[ki], -1));
		  seg_lengths.push_back(dist);
		  seg_endpt.push_back(make_pair(pos, pq));
		}
	    }
	}
    }

  vector<shared_ptr<ftSurface> > dummy_faces;
  if (seg_pnt.size() > 0)
    {
      // Divide face with respect to planes passing through the specified
      // points and the specified axis
      if (loop_idx >= 0)
	{
	  double min_dist2;
	  if (min_dist.size() == 0)
	    min_dist2 = 1.5*min_rad;
	  else
	    {
	      min_dist2 = min_dist[0];
	      for (size_t km=1; km<min_dist.size(); ++km)
		min_dist2 = std::min(min_dist2, min_dist[km]);
	    }
	  return divideByPlanes(seg_pnt, mid, axis, loop_idx, half_holes, 
				holes, min_dist2);
	}
      else
	return divideByPlanes(seg_pnt, seg_norm, half_holes, min_dist);
    }
  else if (hole_idx.size() > 0)
    {
      // Make two splits from hole to hole to combine them into one hole
      // instead of separating the holes
      vector<shared_ptr<ftSurface> > faces =
	holeToHoleSplit(half_holes, holes, hole_idx, seg_lengths, seg_endpt, loop_idx);
      if (faces.size() > 0)
	return faces;
    }
  return dummy_faces;
}

//==========================================================================
void RegularizeFace::getSegPntRadialSplit(hole_info& hole1, hole_info& hole2,
					  double len_fac,
					  vector<Point>& seg_pnt,
					  vector<Point>& seg_norm,
					  vector<double>& min_dist)
//==========================================================================
{
  // Make a plane orthogonal to the chord between the holes midpoints
  size_t ki;
  Point tmp1 = hole1.hole_centre_;
  Point tmp2 = hole2.hole_centre_;
  Point mid = 0.5*(tmp1 + tmp2);
  Point normal = tmp2 - tmp1;
  normal.normalize();
  
  // Intersect with plane
  shared_ptr<ParamSurface> surf = face_->surface();
  shared_ptr<BoundedSurface> bd_surf;
  vector<shared_ptr<CurveOnSurface> > trim_seg;
  trim_seg = BoundedUtils::getPlaneIntersections(surf, mid, normal, 
						 epsge_, bd_surf);

  if (trim_seg.size() == 0)
    return;   // Nothing to do

#ifdef DEBUG_REG
  std::ofstream of("seg_split.g2");
  for (ki=0; ki<trim_seg.size(); ++ki)
    {
      shared_ptr<ParamCurve> cv = trim_seg[ki]->spaceCurve();
      cv->writeStandardHeader(of);
      cv->write(of);
    }
#endif
      

  // Select a point on the intersection curve(s)
  int idx = ((int)trim_seg.size() / 2);
  double par = 0.5*(trim_seg[idx]->startparam() + trim_seg[idx]->endparam());
  Point pos = trim_seg[idx]->ParamCurve::point(par);

#ifdef DEBUG_REG
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pos << std::endl;
#endif

  // Intersect with planes through the midpoint and each hole to find
  // points at the hole boundaries
  Point dir1 = pos - tmp1;
  Point normal1 = dir1.cross(hole1.hole_axis_);
  normal1.normalize();

  // Compute intersections between the face and this plane and pick
  // the appropriate piece(s)
  vector<shared_ptr<CurveOnSurface> > trim_seg1 = 
    BoundedUtils::getPlaneIntersections(surf, pos, normal1, 
					epsge_, bd_surf);

  // Keep only the closest intersection
  if (trim_seg1.size() > 0)
    {
      int min_idx = -1;
      double min_dist = HUGE;
      for (ki=0; ki<trim_seg1.size(); ++ki)
	{
	  double par, dist;
	  Point close;
	  trim_seg1[ki]->closestPoint(pos, trim_seg1[ki]->startparam(), 
				      trim_seg1[ki]->endparam(), par,
				      close, dist);
	  if (dist < min_dist)
	    {
	      min_dist = dist;
	      min_idx = (int)ki;
	    }
	}
      std::swap(trim_seg1[0], trim_seg1[min_idx]);
      trim_seg1.erase(trim_seg1.begin()+1, trim_seg1.end());
    }
	  

 #ifdef DEBUG_REG
  for (ki=0; ki<trim_seg1.size(); ++ki)
    {
      shared_ptr<ParamCurve> cv = trim_seg1[ki]->spaceCurve();
      cv->writeStandardHeader(of);
      cv->write(of);
    }
#endif
      
  Point p1 = pos;
  for (ki=0; ki<trim_seg1.size(); ++ki)
    {
      Point pos1 = 
	trim_seg1[ki]->ParamCurve::point(trim_seg1[ki]->startparam());
      Point pos2 = 
	trim_seg1[ki]->ParamCurve::point(trim_seg1[ki]->endparam());
      if (tmp1.dist(pos1) < tmp1.dist(p1))
	p1 = pos1;
      if (tmp1.dist(pos2) < tmp1.dist(p1))
	p1 = pos2;
    }

  Point dir2 = pos - tmp2;
  Point normal2 = dir2.cross(hole2.hole_axis_);
  normal2.normalize();

  vector<shared_ptr<CurveOnSurface> > trim_seg2 = 
    BoundedUtils::getPlaneIntersections(surf, pos, normal2, 
					epsge_, bd_surf);

  // Keep only the closest intersection
  if (trim_seg2.size() > 0)
    {
      int min_idx = -1;
      double min_dist = HUGE;
      for (ki=0; ki<trim_seg2.size(); ++ki)
	{
	  double par, dist;
	  Point close;
	  trim_seg2[ki]->closestPoint(pos, trim_seg2[ki]->startparam(), 
				      trim_seg2[ki]->endparam(), par,
				      close, dist);
	  if (dist < min_dist)
	    {
	      min_dist = dist;
	      min_idx = (int)ki;
	    }
	}
      std::swap(trim_seg2[0], trim_seg2[min_idx]);
      trim_seg2.erase(trim_seg2.begin()+1, trim_seg2.end());
    }
	  

#ifdef DEBUG_REG
  for (ki=0; ki<trim_seg2.size(); ++ki)
    {
      shared_ptr<ParamCurve> cv = trim_seg2[ki]->spaceCurve();
      cv->writeStandardHeader(of);
      cv->write(of);
    }
#endif

  Point p2 = pos;
  for (ki=0; ki<trim_seg2.size(); ++ki)
    {
      Point pos1 = 
	trim_seg2[ki]->ParamCurve::point(trim_seg2[ki]->startparam());
      Point pos2 = 
	trim_seg2[ki]->ParamCurve::point(trim_seg2[ki]->endparam());
      if (tmp2.dist(pos1) < tmp2.dist(p2))
	p2 = pos1;
      if (tmp2.dist(pos2) < tmp2.dist(p2))
	p2 = pos2;
    }

#ifdef DEBUG_REG
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << "2" << std::endl;
  of << p1 << std::endl;
  of << p2 << std::endl;
#endif

  // Compute distances and split points
  double d1 = p1.dist(pos);
  double d2 = p2.dist(pos);
  double r1 = tmp1.dist(p1); 
  double r2 = tmp2.dist(p2); 
  double dist_fac = 0.5; //0.25;
  double max_edge_len = len_fac*(r1 + r2); 

  int nmb_cand_split = 0;
  if (cand_split_.size() > 0)
    nmb_cand_split = nmbSplitPattern(p1, p2);

  double fac = 1.2;
  if ((d1+d2 > max_edge_len && nmb_cand_split != 1) || 
      nmb_cand_split > 1)
    {
      // Define two split points
      if (d1 < 2.0*r1)
	d1 = 2.0*r1;
      double t1 = fac*r1/d1;
      Point tmp0 = (1.0-t1)*p1 + t1*pos;
      seg_pnt.push_back(tmp0);
      Point tmpnorm = pos - p1;
      tmpnorm.normalize();
      seg_norm.push_back(tmpnorm);
      min_dist.push_back(dist_fac*d1);
      
      if (d2 < 2.0*r2)
	d2 = 2.0*r2;
      double t2 = fac*r2/d2;
      tmp0 = (1.0-t2)*p2 + t2*pos;
      seg_pnt.push_back(tmp0);
      tmpnorm = pos - p2;
      tmpnorm.normalize();
      seg_norm.push_back(tmpnorm);
      min_dist.push_back(dist_fac*d2);
    }
  else
    {
      seg_pnt.push_back(pos);
      seg_norm.push_back(normal);
      min_dist.push_back(dist_fac*(d1+d2));
    }
}

//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::isolateHolesRadially2(vector<vector<ftEdge*> >& half_holes,
				      const Point& mid, const Point& axis,
				      int loop_idx,
				      vector<hole_info>& holes, 
				      vector<int>& perm)
//==========================================================================
{
  // Fetch all vertices on the outer loop
  vector<shared_ptr<Vertex> > vx = face_->getBoundaryLoop(0)->getVertices();
  removeInsignificantVertices(vx);  

  // For each pair of consequetive holes, make one or two splitting curves
  // between the holes to separate them
  // In a later version, we may also consider the distance to the outer
  // boundary
  shared_ptr<ParamSurface> surf = face_->surface();
  vector<Point> seg_pnt;
  vector<Point> seg_norm;
  size_t ki, kh;
  double min_rad = holes[0].hole_radius_;
  size_t nmb_holes =  holes.size();
  vector<double> min_dist;

#ifdef DEBUG_REG
  std::ofstream of0("sel_vx.g2");
#endif

  // Estimate average chord length between holes (flattened)
  double av_chord = 0.0;
  double chord_fac = 0.2;
  for (ki=0; ki<nmb_holes; ++ki)
    {
      kh = ki+1;
      if (kh == holes.size())
	kh = 0;
      Point chord = holes[perm[kh]].hole_centre_ - holes[perm[ki]].hole_centre_;
      av_chord += (chord.length() - holes[perm[ki]].hole_radius_ -
		   holes[perm[kh]].hole_radius_);
    }
  av_chord /= (double)(nmb_holes);

  vector<int> used_idx;
  int nmb_vx_split = 0;
  for (ki=0; ki<nmb_holes; ++ki)
    {
      kh = ki+1;
      if (kh == holes.size())
	kh = 0;

      // Check if a cylinder or a plane intersection is appropriate
      double len_fac = 6.0;
      double d1 = mid.dist(holes[perm[ki]].hole_centre_);
      double d2 = mid.dist(holes[perm[kh]].hole_centre_);
      double min_rad = std::min(holes[perm[ki]].hole_radius_,
				holes[perm[kh]].hole_radius_);
      double cylinder_fac = 0.5;
      shared_ptr<BoundedSurface> bd_surf;
      vector<shared_ptr<CurveOnSurface> > trim_seg;
      Point tmp1 = holes[perm[ki]].hole_centre_;
      Point tmp2 = holes[perm[kh]].hole_centre_;
      if (fabs(d1 - d2) < cylinder_fac*min_rad)
	{
	  // Intersect with cylinder
	  double cyl_rad = 0.5*(d1 + d2);
	  trim_seg = BoundedUtils::getCylinderIntersections(surf, mid, 
							    axis, cyl_rad,
							    epsge_, bd_surf);
	  //len_fac = 3.0;
	  
	}
      else
	{
	  // Make a plane following a chord between holes midpoints,
	  // approximately being perpendicular to the hole axis
	  Point mid2 = 0.5*(tmp1 + tmp2);
	  Point dir = tmp2 - tmp1;
	  Point axis2 = 0.5*(holes[perm[ki]].hole_axis_ +
			     holes[perm[kh]].hole_axis_);
	  Point normal = dir.cross(axis2);
	  normal.normalize();

	  // Compute intersections between the face and this plane and pick
	  // the appropriate piece(s)
	  trim_seg = BoundedUtils::getPlaneIntersections(surf, mid2, normal, 
							 epsge_, bd_surf);
	}

      // Remove the segments lying outside the area between the two holes
      size_t kj;
      for (kj=0; kj<trim_seg.size(); )
	{
	  double tmp_par = 0.5*(trim_seg[kj]->startparam() +
				trim_seg[kj]->startparam());
	  Point tmp3 = trim_seg[kj]->ParamCurve::point(tmp_par);
	  if ((tmp3 - tmp1)*(tmp3 - tmp2) >= 0.0)
	    trim_seg.erase(trim_seg.begin()+kj);
	  else
	    kj++;
	}
      
      // Estimate length of remaining segments
      double len = 0.0;
      Point p1 = tmp2, p2 = tmp1;
      for (kj=0; kj<trim_seg.size(); ++kj)
	{
	  len += trim_seg[kj]->estimatedCurveLength();
	  Point pos1 = 
	    trim_seg[kj]->ParamCurve::point(trim_seg[kj]->startparam());
	  Point pos2 = 
	    trim_seg[kj]->ParamCurve::point(trim_seg[kj]->endparam());
	  if (tmp1.dist(pos1) < tmp1.dist(p1))
	    p1 = pos1;
	  if (tmp1.dist(pos2) < tmp1.dist(p1))
	    p1 = pos2;
	  if (tmp2.dist(pos1) < tmp2.dist(p2))
	    p2 = pos1;
	  if (tmp2.dist(pos2) < tmp2.dist(p2))
	    p2 = pos2;
	}
      // len = p1.dist(p2); // ?
      if (p1.dist(p2) < chord_fac*av_chord)
	continue;  // Do not split betwen very close holes at this stage

#ifdef DEBUG_REG
      std::ofstream out_file("split_segments.g2");
      for (kj=0; kj<trim_seg.size(); ++kj)
	{
	  shared_ptr<ParamCurve> cv = trim_seg[kj]->spaceCurve();
	  cv->writeStandardHeader(out_file);
	  cv->write(out_file);
	}
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "2" << std::endl;
      out_file << tmp1 << std::endl;
      out_file << tmp2 << std::endl;
      out_file << "400 1 0 4 0 255 0 255" << std::endl;
      out_file << "2" << std::endl;
      out_file << p1 << std::endl;
      out_file << p2 << std::endl;
#endif

      if (trim_seg.size() != 1)
	continue;  // Try to simplify the situation before handling this case

      // Identify vertices from the outer boundary which are relevant 
      // for splitting between the current holes
      size_t kr;
      double length_fac = 0.1;
      double ang_fac = 0.1*M_PI;
      vector<pair<Point, Point> > cand_seg;
      vector<pair<pair<double, bool>,pair<double,double> > > vx_info;
      vector<int> vx_idx;
      for (kj=0; kj<vx.size(); ++kj)
	{
	  // Check if the vertex is used in a split already
	  for (kr=0; kr<used_idx.size(); ++kr)
	    if ((int)kj == used_idx[kr])
	      break;
	  if (kr < used_idx.size())
	    continue;  // Do not reuse vertex at this recursion level

	  Point vx_pnt = vx[kj]->getVertexPoint();
	  for (kr=0; kr<trim_seg.size(); ++kr)
	    {
	      double close_par, close_dist;
	      Point close_point, close_tan;
	      Point close_param;
	      vector<Point> pts(2);
	      double tmin = trim_seg[kr]->startparam();
	      double tmax = trim_seg[kr]->endparam();
	      trim_seg[kr]->closestPoint(vx_pnt, tmin, tmax, close_par,
					 close_point, close_dist);
	      if (close_par < tmin+length_fac*(tmax-tmin))
		{
		  trim_seg[kr]->point(pts, tmin+length_fac*(tmax-tmin), 1);
		  close_point = pts[0];
		  close_tan = pts[1];
		  close_param = 
		    trim_seg[kr]->parameterCurve()->point(tmin+length_fac*(tmax-tmin));
		}
	      else if (close_par > tmax-length_fac*(tmax-tmin))
		{
		  trim_seg[kr]->point(pts, tmax-length_fac*(tmax-tmin), 1);
		  close_point = pts[0];
		  close_tan = pts[1];
		  close_param = 
		    trim_seg[kr]->parameterCurve()->point(tmax-length_fac*(tmax-tmin));
		}
	      else
		{
		  trim_seg[kr]->point(pts, close_par, 1);
		  close_tan = pts[1];
		  close_param = trim_seg[kr]->parameterCurve()->point(close_par);
		}
	      double ang = close_tan.angle(close_point - vx_pnt);
	      if (fabs(0.5*M_PI-ang) < ang_fac)
		{
		  double ang2;
		  bool T_joint;
		  bool connected = getVertexProperties(vx[kj], close_param, 
						       ang2, T_joint);
		  if (connected)
		    {
		      Point tmp_vec = close_tan.cross(vx_pnt - close_point);
		      Point norm = tmp_vec.cross(vx_pnt - close_point);
		      norm.normalize();
		      cand_seg.push_back(make_pair(close_point, norm));
		      vx_info.push_back(make_pair(make_pair(ang2, T_joint),
						  make_pair(ang, vx_pnt.dist(close_point))));
		      vx_idx.push_back((int)kj);
		    }
		}
	    }
	}

      // Sort candidate vertices according to importance
      // Dismiss vertices far away from the chord
      int kr2;
      double min_d = 1.0e8, max_d = 0.0;
      double d_fac = 0.3;
      for (kr2=0; kr2<(int)vx_info.size(); ++kr2)
	{
	  min_d = std::min(min_d, vx_info[kr2].second.second);
	  max_d = std::max(max_d, vx_info[kr2].second.second);
	}
      double d_lim = std::max(d_fac*(max_d - min_d), 0.5*min_d);
      int nmb_cand = (int)cand_seg.size();
      for (kr2=0; kr2<nmb_cand; ++kr2)
	if (vx_info[kr2].second.second > min_d + d_lim)
	  {
	    std::swap(vx_info[kr2], vx_info[nmb_cand-1]);
	    std::swap(cand_seg[kr2], cand_seg[nmb_cand-1]);
	    std::swap(vx_idx[kr2], vx_idx[nmb_cand-1]);
	    nmb_cand--;
	    kr2--;
	  }

      // Dismiss very convex corners
      double min_a = 2.0*M_PI, max_a = 0.0;
      for (kr2=0; kr2<nmb_cand; ++kr2)
	{
	  min_a = std::min(min_a, vx_info[kr2].first.first);
	  max_a = std::max(max_a, vx_info[kr2].first.first);
	}
      double a_lim = std::max(d_fac*(max_a - min_a), M_PI/6.0);
       for (kr2=0; kr2<nmb_cand; ++kr2)
	if (vx_info[kr2].first.first < max_a - a_lim)
	  {
	    std::swap(vx_info[kr2], vx_info[nmb_cand-1]);
	    std::swap(cand_seg[kr2], cand_seg[nmb_cand-1]);
	    std::swap(vx_idx[kr2], vx_idx[nmb_cand-1]);
	    nmb_cand--;
	    kr2--;
	  }
      
       // Sort according to tangent angle
       for (int kr1=0; kr1<nmb_cand; ++kr1)
	 {
	   double ang1 = fabs(0.5*M_PI-vx_info[kr1].second.first); 
	 for (kr2=kr1+1; kr2<nmb_cand; ++kr2)
	   {
	     double ang2 = fabs(0.5*M_PI-vx_info[kr2].second.first);
	     if (ang2 < ang1)
	       {
		 std::swap(vx_info[kr2], vx_info[kr1]);
		 std::swap(cand_seg[kr2], cand_seg[kr1]);
		 std::swap(vx_idx[kr2], vx_idx[kr1]);
	       }
	   }
	 }

      // Define split points
      double r1 = tmp1.dist(p1); //holes[perm[ki]].hole_radius_;
      double r2 = tmp2.dist(p2); //holes[perm[kj]].hole_radius_;
       double max_edge_len = len_fac*(r1 + r2); //4.0*(r1 + r2);
     //min_rad = std::min(min_rad, r2);
      double dist_fac = 0.25;
      min_dist.push_back(dist_fac*p1.dist(p2));

      int nmb_cand_split = 0;
      vector<pair<Point,Point> > pattern;
      vector<pair<int,int> > pattern_vx;
      if (cand_split_.size() > 0)
	{
	  vector<double> ang_info(vx_info.size());
	  for (kr2=0; kr2<(int)vx_info.size(); ++kr2)
	    ang_info[kr2] = vx_info[kr2].second.first;
	  nmb_cand_split = fetchSplitPattern(p1, p2, vx, vx_idx, nmb_cand,
					     pattern, pattern_vx);
	  //nmb_cand_split = nmbSplitPattern(p1, p2);
	}

      double fac = 2.5; //1.2;
      if (cand_split_.size() > 0 && pattern.size() == 0)
	continue;   // A split here will not fit with the opposite face
      else if (((len > max_edge_len && nmb_cand_split != 1) || 
		nmb_cand_split > 1) &&
	       pattern.size() != 1)
	{
	  // Define two split points
	  double d1 = p1.dist(p2);
	  if (d1 < 2.0*(r1 + r2))
	    d1 = 2.0*(r1 + r2);
	  double t1 = fac*r1/d1;
	  Point tmp0 = (1.0-t1)*p1 + t1*p2;
	  Point tmpnorm = p2 - p1;
	  tmpnorm.normalize();
	  int idx = -1; 
	  if (nmb_cand > 0)
	    {
	      double dd = 0.5*d1;
	      for (kr2=0; kr2<nmb_cand; ++kr2)
		{
		  double dd2 = cand_seg[kr2].first.dist(tmp0);
		  for (size_t kr3=0; kr3<pattern_vx.size(); ++kr3)
		    if (pattern_vx[kr3].first == vx_idx[kr2] ||
			pattern_vx[kr3].second == vx_idx[kr2])
		      dd2 = 0.0;   // Select this vertex
		  if (dd2 < dd)
		    {
		      dd = dd2;
		      idx = kr2;
		    }
		}
	      if (idx >= 0)
		{
		  tmp0 = cand_seg[idx].first;
		  seg_pnt.push_back(cand_seg[idx].first);
		  seg_norm.push_back(cand_seg[idx].second);
		  nmb_vx_split++;
		  used_idx.push_back(vx_idx[idx]);

#ifdef DEBUG_REG
		  of0 << "400 1 0 4 0 155 100 255" << std::endl;
		  of0 << "1" << std::endl;
		  of0 << vx[vx_idx[idx]]->getVertexPoint() << std::endl;
#endif
		}
	    }
	  if (idx < 0)
	    {
	      seg_pnt.push_back(tmp0);
	      seg_norm.push_back(tmpnorm);
	    }
	  Point tmp1 = (p1.dist(tmp0) > p2.dist(tmp0)) ? p1 + 0.5*(tmp0-p1) :
	    p2 - 0.5*(p2 - tmp0);
	  int idx2 = -1;
	  if (nmb_cand > 0)
	    {
	      double dd = 1.0e8;
	      for (kr2=0; kr2<nmb_cand; ++kr2)
		{
		  if (kr2 == idx)
		    continue;
		  double dd2 = cand_seg[kr2].first.dist(tmp1);
		  for (size_t kr3=0; kr3<pattern_vx.size(); ++kr3)
		    if (pattern_vx[kr3].first == vx_idx[kr2] ||
			pattern_vx[kr3].second == vx_idx[kr2])
		      dd2 = 0.0;   // Select this vertex
		  if (dd2 < dd)
		    {
		      dd = cand_seg[kr2].first.dist(tmp1);
		      idx2 = kr2;
		    }
		}
	      if (idx2 >= 0)
		{
		  seg_pnt.push_back(cand_seg[idx2].first);
		  seg_norm.push_back(cand_seg[idx2].second);
		  nmb_vx_split++;
		  used_idx.push_back(vx_idx[idx2]);
#ifdef DEBUG_REG
		  of0 << "400 1 0 4 0 155 100 255" << std::endl;
		  of0 << "1" << std::endl;
		  of0 << vx[vx_idx[idx2]]->getVertexPoint() << std::endl;
#endif
		}
	    }
	  if (idx2 < 0 && !(idx >= 0 && nmb_vx_split > 0))
	    {
	      seg_pnt.push_back(tmp1);
	      seg_norm.push_back(tmpnorm);
	    }
	  min_dist[min_dist.size()-1] *= 0.5;
	  min_dist.push_back(min_dist[min_dist.size()-1]);
	}
      else
	{
	  double t1 = r1/(r1 + r2);
	  Point tmp0 = t1*p1 + (1.0-t1)*p2;

	  // One split point
	  if (nmb_cand > 0)
	    {
	      int idx = 0;
	      double dd = 0.75*d1;
	      for (kr2=0; kr2<nmb_cand; ++kr2)
		{
		  double dd2 = cand_seg[kr2].first.dist(tmp0);
		  for (size_t kr3=0; kr3<pattern_vx.size(); ++kr3)
		    if (pattern_vx[kr3].first == vx_idx[kr2] ||
			pattern_vx[kr3].second == vx_idx[kr2])
		      dd2 = 0.0;   // Select this vertex
		  if (dd2 < dd)
		    {
		      dd = dd2;
		      idx = kr2;
		    }
		}
	      seg_pnt.push_back(cand_seg[idx].first);
	      seg_norm.push_back(cand_seg[idx].second);
	      nmb_vx_split++;
	      used_idx.push_back(vx_idx[idx]);
#ifdef DEBUG_REG
	      of0 << "400 1 0 4 0 155 100 255" << std::endl;
	      of0 << "1" << std::endl;
	      of0 << vx[vx_idx[0]]->getVertexPoint() << std::endl;
#endif
	    }
	  else
	    {
	      seg_pnt.push_back(tmp0);
	      Point tmpnorm = p2 - p1;
	      tmpnorm.normalize();
	      seg_norm.push_back(tmpnorm);
	    }
	}
    }

  // Divide face with respect to planes passing through the specified
  // points and the specified axis
  if (loop_idx >= 0)
    return divideByPlanes(seg_pnt, mid, axis, loop_idx, half_holes, 
			  holes, min_rad);
  else
    return divideByPlanes(seg_pnt, seg_norm, half_holes, min_dist);
  }

//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::isolateOneHoleRadially(const Point& mid, const Point& axis,
				      hole_info& hole)
//==========================================================================
{
  // Intersect with cylinder with centre in the axis passing through
  // the hole centre
  shared_ptr<ParamSurface> surf = face_->surface();
  BoundingBox box = surf->boundingBox();
  double box_size = box.high().dist(box.low());
  shared_ptr<BoundedSurface> bd_surf;
  vector<shared_ptr<CurveOnSurface> > trim_seg;
  Point tmp1 = hole.hole_centre_;
  double cyl_rad = mid.dist(tmp1);
  trim_seg = BoundedUtils::getCylinderIntersections(surf, mid, 
						    axis, cyl_rad,
						    epsge_, bd_surf);
  
#ifdef DEBUG_REG
      std::ofstream out_file("split_segments.g2");
      for (size_t kj=0; kj<trim_seg.size(); ++kj)
	{
	  shared_ptr<ParamCurve> cv = trim_seg[kj]->spaceCurve();
	  cv->writeStandardHeader(out_file);
	  cv->write(out_file);
	}
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "1" << std::endl;
      out_file << tmp1 << std::endl;
#endif

      // Select points at the trimming segments at both sides of the hole
      // to govern the split. Create cylinder surface with which to intersect
      // the trim segments
      double fac = 3.0;
      double cyl_rad2 = fac*hole.hole_radius_;
      Point vec = hole.hole_centre_ - mid;
      Cylinder cyl(cyl_rad2, hole.hole_centre_, hole.hole_axis_, vec);
      cyl.setParameterBounds(0.0, -box_size, 2*M_PI, box_size);
      shared_ptr<SplineSurface> cyl_sf(cyl.createSplineSurface());

      vector<Point> seg_pnt;
      vector<Point> seg_norm;
      for (size_t ki=0; ki<trim_seg.size(); ++ki)
	{
	  shared_ptr<SplineCurve> cv(trim_seg[ki]->geometryCurve());

	  vector<pair<double, Point> > int_pts;
	  vector<int> pretop;
	  vector<pair<pair<double,Point>, pair<double,Point> > > int_crvs;
	  intersectCurveSurf(cv.get(), cyl_sf.get(), epsge_, int_pts,
			     pretop, int_crvs);
	  for (size_t kr=0; kr<int_pts.size(); ++kr)
	    {
	      vector<Point> pts(2);
	      cv->point(pts, int_pts[kr].first, 1);
	      size_t kh;
	      for (kh=0; kh<seg_pnt.size(); ++kh)
		{
		  if (pts[0].dist(seg_pnt[kh]) <= epsge_)
		    break;
		}
	      if (kh == seg_pnt.size())
		{
		  seg_pnt.push_back(pts[0]);
		  seg_norm.push_back(pts[1]);
#ifdef DEBUG_REG
		  out_file << "400 1 0 4 0 255 0 255" << std::endl;
		  out_file << "1" << std::endl;
		  out_file << pts[0] << std::endl;
#endif
		}
	    }
	}

      vector<vector<ftEdge*> > half_holes_dummy;
      vector<double> min_dist(seg_pnt.size(), hole.hole_radius_);
      return divideByPlanes(seg_pnt, seg_norm, half_holes_dummy, min_dist);
}

//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::divideByPlanes(vector<Point>& pnts, vector<Point>& normals,
			       vector<vector<ftEdge*> >& half_holes,
			       const vector<double>& level_dist)
//==========================================================================
{
  // Fetch all vertices on the outer loop
  vector<shared_ptr<Vertex> > vx = face_->getBoundaryLoop(0)->getVertices();

  // Make copy
  vector<shared_ptr<Vertex> > vx2(vx.begin(), vx.end());

  removeInsignificantVertices(vx);  

  // Remove vertices belonging to half holes from the candidate vertices
  removeHalfHoleVx(vx, half_holes);

  size_t ki, kj;
  for (ki=0; ki<vx2.size(); )
    {
      for (kj=0; kj<vx.size(); ++kj)
	if (vx[kj].get() == vx2[ki].get())
	  break;
      if (kj < vx.size())
	vx2.erase(vx2.begin()+ki);
      else
	ki++;
    }
  
  // For each plane, perform intersection
  vector<shared_ptr<CurveOnSurface> > segments;
  shared_ptr<BoundedSurface> bd_sf;
  shared_ptr<ParamSurface> surf = face_->surface();
  for (ki=0; ki<pnts.size(); ++ki)
    {
      Point curr_pnt = pnts[ki];

      // Find the closest vertex with distance less than the level value
      shared_ptr<Vertex> cand, cand2;
      shared_ptr<Vertex> cand3;
      double min_dist = level_dist[ki];
      double prev_dist3 = 1.0e8;
      for (size_t kj=0; kj<vx.size(); ++kj)
	{
	  Point tmp = vx[kj]->getVertexPoint();
	  double dist = fabs((tmp - pnts[ki])*normals[ki]);
	  double dist3 = tmp.dist(curr_pnt);
	  if (dist < min_dist)
	    {
	      min_dist = dist;
	      cand = vx[kj];
	    }
	  if (dist3 < prev_dist3)
	    {
	      prev_dist3 = dist3;
	      cand3 = vx[kj];
	    }
	}

      double level_ang = 0.1*M_PI;
      Point parval1, parval2;
      bool found = false;
      if (cand.get())
	{
	  // Check if a splitting curve in an opposite face can be mimiced
	  if (cand_split_.size() > 0)
	    {
	      Point vx_point = cand->getVertexPoint();
	      found = fetchPatternSplit2(vx_point, parval1, parval2, level_dist[ki], 
					 level_ang, normals[ki], true);
	    }

	  min_dist = level_dist[ki];
	  int idx = -1;
	  double idx_dist = HUGE;
	  if (found)
	    {
	      // Search for a vertex close to the found parameter value
	      Point pos = face_->point(parval2[0], parval2[1]);
	      size_t kj;
	      for (kj=0; kj<vx.size(); ++kj)
		{
		  if (vx[kj].get() == cand.get())
		    continue;  // Vertex already selected
		  
		  double dist = vx[kj]->getVertexPoint().dist(pos);
		  if (dist < idx_dist)
		    {
		      idx = (int)kj;
		      idx_dist = dist;
		    }
		}
	      if (idx_dist < level_dist[ki])
		{
		  // Check if this vertex is the end of a chain
		  vector<shared_ptr<Vertex> > met_already;
		  vector<shared_ptr<Vertex> > v1_end;
		  v1_end = RegularizeUtils::endVxInChain(face_, face_.get(), 
							 NULL, cand, cand,
							 cand, met_already);
		  for (kj=0; kj<v1_end.size(); ++kj)
		    {
		      if (v1_end[kj].get() == vx[idx].get())
			break;
		    }

		  if (kj<v1_end.size() || idx_dist < tol2_)
		    {
		      // This is a very strict tolerance and should be replaced
		      // after getting some more experience
		      cand2 = vx[idx];
		    }
		}
	    }
	  else
	    {
	      // Search for an alternative vertex at an edge not adjacent to
	      // the current one which also lies within the limit distance
	      // from the specified plane
	      for (size_t kj=0; kj<vx.size(); ++kj)
		{
		  if (vx[kj].get() == cand.get())
		    continue;  // Vertex already selected

		  // Compute angle between this segment and the candidate plane 
		  // intersection
		  Point tmp_pnt = cand->getVertexPoint();
		  double angle = getSegmentAngle(cand, vx[kj], tmp_pnt, 
						 normals[ki]);

		  if (vx[kj]->sameEdge(cand.get()) && angle >= level_ang)
		    continue;  // Edge adjacent to already found vertex

		  Point tmp = vx[kj]->getVertexPoint();
		  double dist = fabs((tmp - pnts[ki])*normals[ki]);
		  if ((tmp-pnts[ki])*(tmp_pnt-pnts[ki]) > 0.0)
		    continue;  // The two vertices is on the same side of the
		  // chord between the holes

		  if ((dist < min_dist && angle < 2.0*level_ang) || 
		      (angle < level_ang && dist < 2.0*min_dist))
		    {
		      min_dist = dist;
		      cand2 = vx[kj];
		    }
		}
	    }
	}

      vector<shared_ptr<CurveOnSurface> > trim_segments;
      if (cand.get() && (found || cand2.get()))
	{
	  // Find division curve between vertices
	  parval1 = cand->getFacePar(face_.get());
	  if (cand2)
	    parval2 = cand2->getFacePar(face_.get());
	  trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
							 parval2, epsge_,
							 bd_sf);

#ifdef DEBUG_REG
	  std::ofstream of10("trim_crvs0.g2");
	  for (size_t ki2=0; ki2<trim_segments.size(); ++ki2)
	    {
	      trim_segments[ki2]->spaceCurve()->writeStandardHeader(of10);
	      trim_segments[ki2]->spaceCurve()->write(of10);
	    }
#endif
	  // Remove trim segments very distant from the initial point
	  size_t kj;
	  for (kj=0; kj<trim_segments.size();)
	    {
	      Point close_pt;
	      double close_par, close_dist;
	      // trim_segments[kj]->closestPoint(curr_pnt, trim_segments[kj]->startparam(),
	      // 				  trim_segments[kj]->endparam(), close_par,
	      // 				  close_pt, close_dist);
	      // if (close_dist > tol2_/*level_dist[ki]*/)
	      trim_segments[kj]->closestPoint(pnts[ki], trim_segments[kj]->startparam(),
					      trim_segments[kj]->endparam(), close_par,
					      close_pt, close_dist);
	      if (close_dist > level_dist[ki])
		trim_segments.erase(trim_segments.begin()+kj);
	      else
		kj++;
	    }
	}
      if (trim_segments.size() == 0 && cand.get())
	{
	  // Move plane to pass through the candidate vertex
	  curr_pnt = cand->getVertexPoint();
	  trim_segments = BoundedUtils::getPlaneIntersections(surf, curr_pnt,
							      normals[ki], 
							      epsge_,
							      bd_sf);
#ifdef DEBUG_REG
	  std::ofstream of10("trim_crvs0.g2");
	  for (size_t ki2=0; ki2<trim_segments.size(); ++ki2)
	    {
	      trim_segments[ki2]->spaceCurve()->writeStandardHeader(of10);
	      trim_segments[ki2]->spaceCurve()->write(of10);
	    }
#endif
	  // Remove trim segments very distant from the initial point
	  for (kj=0; kj<trim_segments.size();)
	    {
	      Point close_pt;
	      double close_par, close_dist;
	      // trim_segments[kj]->closestPoint(curr_pnt, trim_segments[kj]->startparam(),
	      // 				  trim_segments[kj]->endparam(), close_par,
	      // 				  close_pt, close_dist);
	      // if (close_dist > tol2_/*level_dist[ki]*/)
	      trim_segments[kj]->closestPoint(pnts[ki], trim_segments[kj]->startparam(),
					      trim_segments[kj]->endparam(), close_par,
					      close_pt, close_dist);
	      if (close_dist > level_dist[ki])
		trim_segments.erase(trim_segments.begin()+kj);
	      else
		kj++;
	    }
	}
      if (trim_segments.size() == 0)
	{
	  // Split with plane
	  trim_segments = BoundedUtils::getPlaneIntersections(surf, pnts[ki],
							      normals[ki], 
							      epsge_,
							      bd_sf);
#ifdef DEBUG_REG
	  std::ofstream of10("trim_crvs0.g2");
	  for (size_t ki2=0; ki2<trim_segments.size(); ++ki2)
	    {
	      trim_segments[ki2]->spaceCurve()->writeStandardHeader(of10);
	      trim_segments[ki2]->spaceCurve()->write(of10);
	    }
#endif
	  // Remove trim segments very distant from the initial point
	  size_t kj;
	  for (kj=0; kj<trim_segments.size();)
	    {
	      Point close_pt;
	      double close_par, close_dist;
	      // trim_segments[kj]->closestPoint(curr_pnt, trim_segments[kj]->startparam(),
	      // 				  trim_segments[kj]->endparam(), close_par,
	      // 				  close_pt, close_dist);
	      // if (close_dist > tol2_/*level_dist[ki]*/)
	      trim_segments[kj]->closestPoint(pnts[ki], trim_segments[kj]->startparam(),
					      trim_segments[kj]->endparam(), close_par,
					      close_pt, close_dist);
	      if (close_dist > level_dist[ki])
		trim_segments.erase(trim_segments.begin()+kj);
	      else
		kj++;
	    }
	}

      // Modify trimming curves if an endpoint is very close to an insignificant vertex
      Point dummy_vec;
      for (size_t kr=0; kr<trim_segments.size(); kr++)
	{
	  Point pt1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
	  Point pt2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());

	  Point trim_parval1(2), trim_parval2(2);
	  Point close_pt;
	  double close_dist, close_par;
	  (void)face_->closestBoundaryPoint(pt1, dummy_vec, trim_parval1[0], 
					    trim_parval1[1], close_pt,
					    close_dist, close_par);
	  (void)face_->closestBoundaryPoint(pt2, dummy_vec, trim_parval2[0], 
					    trim_parval2[1], close_pt,
					    close_dist, close_par);
	  bool modify = false;
	  for (kj=0; kj<vx2.size(); ++kj)
	    {
	      Point vx_pt = vx2[kj]->getVertexPoint();
	      if (vx_pt.dist(pt1) < 10.0*tol2_)
		{
		  trim_parval1 = vx2[kj]->getFacePar(face_.get());
		  modify = true;
		}
	      if (vx_pt.dist(pt2) < 10.0*tol2_)
		{
		  trim_parval2 = vx2[kj]->getFacePar(face_.get());
		  modify = true;
		}
	    }
	  if (modify)
	    {
	      vector<shared_ptr<CurveOnSurface> > trim_segments2 =
		BoundedUtils::getTrimCrvsParam(surf, trim_parval1, trim_parval2,
					       epsge_, bd_sf);
	      if (trim_segments2.size() == 1)
		trim_segments[kr] = trim_segments2[0];
	    }
	}
	  
      segments.insert(segments.end(), trim_segments.begin(), 
		      trim_segments.end());
    }

#ifdef DEBUG_REG
  std::ofstream of0("trim_crvs.g2");
  for (ki=0; ki<segments.size(); ++ki)
    {
      segments[ki]->spaceCurve()->writeStandardHeader(of0);
      segments[ki]->spaceCurve()->write(of0);
    }
#endif
  
  // The trim segments may intersect each other, remove the
  // last one in case of an intersection
  splitTrimSegments(segments);

  if (segments.size() == 0)
    {
      // No split defined. Try to split according to the outer
      // boundary
       vector<shared_ptr<ftSurface> > faces = faceOuterBdFaces(half_holes);
       return faces;
    }

  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_sf, segments, epsge_);


#ifdef DEBUG_REG
  std::ofstream of("split_surf.g2");
  for (size_t kr=0; kr<sub_sfs.size(); ++kr)
    {
      sub_sfs[kr]->writeStandardHeader(of);
      sub_sfs[kr]->write(of);
    }
#endif
  // Fetch info about vertices not belonging to the corners of the
  // initial face
  vector<shared_ptr<Vertex> > non_corner = 
    face_->getNonCornerVertices(bend_);
  removeInsignificantVertices(non_corner);

  vector<shared_ptr<ftSurface> > faces = 
    RegularizeUtils::createFaces(sub_sfs, face_, epsge_, tol2_, angtol_, 
				 non_corner);

  return faces;
}

//==========================================================================
void 
RegularizeFace::splitTrimSegments(vector<shared_ptr<CurveOnSurface> >& segments)
//==========================================================================
{
  int ki, kj;
  for (ki=0; ki<(int)segments.size(); ++ki)
    {
      for (kj=ki+1; kj<(int)segments.size(); ++kj)
	{
	  double par1, par2, dist;
	  Point pt1, pt2;
	  ClosestPoint::closestPtCurves(segments[ki].get(), segments[kj].get(),
			  par1, par2, dist, pt1, pt2);
	  if (dist < epsge_)
	    {
	      Point p1 = segments[ki]->ParamCurve::point(segments[ki]->startparam());
	      Point p2 = segments[ki]->ParamCurve::point(segments[ki]->endparam());
	      Point p3 = segments[kj]->ParamCurve::point(segments[kj]->startparam());
	      Point p4 = segments[kj]->ParamCurve::point(segments[kj]->endparam());

	      if (std::min(p1.dist(p3), p1.dist(p4)) < epsge_ &&
		  std::min(p2.dist(p3), p2.dist(p4)) < epsge_)
		continue;

	      segments.erase(segments.begin()+kj);
	      kj--;
	      // // Split curves
	      // if (p1.dist(pt1) > epsge_ && p2.dist(pt1) > epsge_)
	      // 	{
	      // 	  shared_ptr<CurveOnSurface> cv1 = 
	      // 	    shared_ptr<CurveOnSurface>(segments[ki]->subCurve(segments[ki]->startparam(),
	      // 							      par1));
	      // 	  shared_ptr<CurveOnSurface> cv2 = 
	      // 	    shared_ptr<CurveOnSurface>(segments[ki]->subCurve(par1,
	      // 							      segments[ki]->endparam()));
	      // 	  segments[ki] = cv1;
	      // 	  segments.push_back(cv2);
	      // 	}

	      // if (p3.dist(pt2) > epsge_ && p4.dist(pt2) > epsge_)
	      // 	{
	      // 	  shared_ptr<CurveOnSurface> cv1 = 
	      // 	    shared_ptr<CurveOnSurface>(segments[kj]->subCurve(segments[kj]->startparam(),
	      // 							      par2));
	      // 	  shared_ptr<CurveOnSurface> cv2 = 
	      // 	    shared_ptr<CurveOnSurface>(segments[kj]->subCurve(par2,
	      // 							      segments[kj]->endparam()));
	      // 	  segments[kj] = cv1;
	      // 	  segments.push_back(cv2);
	      // 	}
	    }
	}
    }
}

//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::divideByPlanes(vector<Point>& pnts, 
			       const Point& mid, const Point& axis,
			       int loop_idx,
			       vector<vector<ftEdge*> >& half_holes,
			       vector<hole_info>& holes, 
			       double level_dist)
//==========================================================================
{
#ifdef DEBUG_REG
  std::ofstream of0("seg_pnts.g2");
  of0 << "400 1 0 4 0 255 0 255" << std::endl;
  of0 << pnts.size() << std::endl;
  for (size_t k2=0; k2<pnts.size(); ++k2)
    of0 << pnts[k2] << std::endl;
#endif

  // Fetch all vertices on the outer loop
  vector<shared_ptr<Vertex> > notsignvx;  // Storate of non-significant vertices
  vector<shared_ptr<Vertex> > vx = face_->getBoundaryLoop(0)->getVertices();
  notsignvx.insert(notsignvx.end(), vx.begin(), vx.end());
  removeInsignificantVertices(vx);

  // Remove vertices belonging to half holes from the candidate vertices
  removeHalfHoleVx(vx, half_holes);

  // Fetch vertices belonging to the hole
  vector<shared_ptr<Vertex> > hole_vx = 
    face_->getBoundaryLoop(loop_idx+1)->getVertices();
  notsignvx.insert(notsignvx.end(), hole_vx.begin(), hole_vx.end());
  removeInsignificantVertices(hole_vx);

   // Extract non-significant vertices
   for (size_t kr=0; kr<vx.size(); ++kr)
     {
       vector<shared_ptr<Vertex> >::iterator vxp =
	 std::find(notsignvx.begin(), notsignvx.end(), vx[kr]);
       if (vxp != notsignvx.end())
	 notsignvx.erase(vxp);
     }
   for (size_t kr=0; kr<hole_vx.size(); ++kr)
     {
       vector<shared_ptr<Vertex> >::iterator vxp =
	 std::find(notsignvx.begin(), notsignvx.end(), hole_vx[kr]);
       if (vxp != notsignvx.end())
	 notsignvx.erase(vxp);
     }
   

  // For each plane, perform intersection
  vector<shared_ptr<CurveOnSurface> > segments;
  shared_ptr<BoundedSurface> bd_sf;
  shared_ptr<ParamSurface> surf = face_->surface();
  for (size_t ki=0; ki<pnts.size(); ++ki)
    {
      Point curr_pt = pnts[ki];
      Point norm = (curr_pt - mid).cross(axis);
      norm.normalize();

      // Find the closest vertex with distance less than the level value
      shared_ptr<Vertex> cand, cand2;
      double min_dist = level_dist;
      for (size_t kj=0; kj<vx.size(); ++kj)
	{
	  Point tmp = vx[kj]->getVertexPoint();

	  // Make sure that the vertex is at the same side of the axis
	  // as the initial point
	  if ((mid-curr_pt)*(mid-tmp) < 0.0)
	    continue;

	  double dist = fabs((tmp - curr_pt)*norm);
	  if (dist < min_dist)
	    {
	      min_dist = dist;
	      cand = vx[kj];
	    }
	}

      if (cand.get())
	{
	  // Check that the candidate vertex will not be in conflict with
	  // any holes
	  size_t kj;
	  double rfac = 1.2;
	  Point vec2 = cand->getVertexPoint() - curr_pt;
	  double dist2 = fabs(vec2*norm);
	  for (kj=0; kj<holes.size(); ++kj)
	    {
	      Point vec1 = holes[kj].hole_centre_ - curr_pt;
	      double dist1 = 
		fabs((holes[kj].hole_centre_ - cand->getVertexPoint())*norm);
	      if (vec1*vec2 > 0.0 && dist1 < rfac*holes[kj].hole_radius_ &&
		  dist2 <  level_dist)
		break;
	    }
	  if (kj < holes.size())
	    cand.reset();
	}

      if (cand.get())
	{
	  // Remove vertex from candidate list so it will not be used 
	  // again
	  vector<shared_ptr<Vertex> >::iterator vxp =
	    std::find(vx.begin(), vx.end(), cand);
	  if (vxp != vx.end())
	    vx.erase(vxp);
	}


       if (cand.get())
	{
	  // Reset plane
	  curr_pt = cand->getVertexPoint();
	  norm = (curr_pt - mid).cross(axis);
	}

      // Check if an inner vertex approximately defines the same plane
      // double level_ang = M_PI/10.0;
      min_dist = level_dist;

      for (size_t kj=0; kj<hole_vx.size(); ++kj)
	{
	  Point vx_pt = hole_vx[kj]->getVertexPoint();
	  double dist = fabs((vx_pt - curr_pt)*norm);
	  double side = (vx_pt - mid)*(curr_pt - mid);
	  if (dist < min_dist && side > 0.0)
	    {
	      min_dist = dist;
	      cand2 = hole_vx[kj];
	    }
	}

      if (cand2.get())
	{
	  // Remove vertex from candidate list so it will not be used 
	  // again
	  vector<shared_ptr<Vertex> >::iterator vxp =
	    std::find(hole_vx.begin(), hole_vx.end(), cand2);
	  if (vxp != hole_vx.end())
	    hole_vx.erase(vxp);
	}


      vector<shared_ptr<CurveOnSurface> > trim_segments;
      if (cand.get() && cand2.get())
	{
	  // Find division curve between vertices
	  Point parval1 = cand->getFacePar(face_.get());
	  Point parval2 = cand2->getFacePar(face_.get());
	  trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
							 parval2, epsge_,
							 bd_sf);
	}
      else
	{
	  // Split with plane
	  if (cand2.get())
	    {
	      // Modify plane
	      curr_pt = cand2->getVertexPoint();
	      norm = (curr_pt - mid).cross(axis);
	    }
	      
	  trim_segments = BoundedUtils::getPlaneIntersections(surf, 
							      curr_pt,
							      norm, 
							      epsge_,
							      bd_sf);
	  // Keep only the segments lying at the correct side of the
	  // axis
	  size_t kj;
	  size_t nmb_seg = trim_segments.size();
	  for (kj=0; kj<nmb_seg; )
	    {
	      Point tmp1 = curr_pt - mid;
	      tmp1 -= ((curr_pt - mid)*axis)*axis;
	      Point tmp2;
	      trim_segments[kj]->point(tmp2, trim_segments[kj]->startparam());
	      Point tmp3 = tmp2 - mid;
	      tmp3 -=((tmp2 - mid)*axis)*axis;
	      if (tmp1*tmp3 < 0.0)
		{
		  trim_segments.erase(trim_segments.begin()+kj);
		  nmb_seg--;
		}
	      else
		{
		  // Control trim segments towards the non-significant vertices
		  size_t kr;
		  Point tmp4;
		  trim_segments[kj]->point(tmp4, trim_segments[kj]->endparam());
		  Point parval1, parval2;
		  for (kr=0; kr<notsignvx.size(); ++kr)
		    {
		      Point tmp = notsignvx[kr]->getVertexPoint();
		      if (tmp.dist(tmp2) < tol2_)
			parval1 = notsignvx[kr]->getFacePar(face_.get());
		      else if (tmp.dist(tmp4) < tol2_)
			parval2 = notsignvx[kr]->getFacePar(face_.get());
		    }
		  if (parval1.dimension() > 0 || parval2.dimension() > 0)
		    {
		      // Modify curve to fit with vertex
		      if (parval1.dimension() == 0)
			{
			  double u, v, dist;
			  Point close;
			  face_->closestPoint(tmp2, u, v, close, dist, epsge_);
			  parval1 = Point(u,v);
			}
		      if (parval2.dimension() == 0)
			{
			  double u, v, dist;
			  Point close;
			  face_->closestPoint(tmp4, u, v, close, dist, epsge_);
			  parval2 = Point(u,v);
			}
		      vector<shared_ptr<CurveOnSurface> > trim_segments2 = 
			BoundedUtils::getTrimCrvsParam(surf, parval1,
						       parval2, epsge_,
						       bd_sf);
		      trim_segments.erase(trim_segments.begin()+kj);
		      nmb_seg--;
		      trim_segments.insert(trim_segments.end(), 
					   trim_segments2.begin(),
					   trim_segments2.end());
		    }
		  else
		    kj++;
		}
	    }

	}

      // Remove trim segments very distant from the initial point
      size_t kj;
      for (kj=0; kj<trim_segments.size(); ++kj)
	{
	  Point close_pt;
	  double close_par, close_dist;
	  trim_segments[kj]->closestPoint(curr_pt, trim_segments[kj]->startparam(),
					  trim_segments[kj]->endparam(), close_par,
					  close_pt, close_dist);
	  if (close_dist > tol2_/*level_dist[ki]*/)
	    trim_segments.erase(trim_segments.begin()+kj);
	  else
	    kj++;
	}

      segments.insert(segments.end(), trim_segments.begin(), 
		      trim_segments.end());
    }
	  
#ifdef DEBUG_REG
	  std::ofstream out_file("split_segments.g2");
	  for (size_t kj=0; kj<segments.size(); ++kj)
	    {
	      shared_ptr<ParamCurve> cv = segments[kj]->spaceCurve();
	      cv->writeStandardHeader(out_file);
	      cv->write(out_file);
	    }
#endif
      
  
  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_sf, segments, epsge_);

#ifdef DEBUG_REG
  std::ofstream of("split_surf.g2");
  for (size_t kr=0; kr<sub_sfs.size(); ++kr)
    {
      sub_sfs[kr]->writeStandardHeader(of);
      sub_sfs[kr]->write(of);
    }
#endif
  // Fetch info about vertices not belonging to the corners of the
  // initial face
  vector<shared_ptr<Vertex> > non_corner = 
    face_->getNonCornerVertices(bend_);
  removeInsignificantVertices(non_corner);

  vector<shared_ptr<ftSurface> > faces = 
    RegularizeUtils::createFaces(sub_sfs, face_, epsge_, tol2_,
				 angtol_, non_corner);

  return faces;
}

//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::holeToHoleSplit(vector<vector<ftEdge*> >& half_holes,
				vector<hole_info>& holes, 
				vector<pair<int,int> >& hole_idx,
				vector<double>& seg_lengths,
				vector<pair<Point, Point> >& seg_endpt,
				int loop_idx)
//==========================================================================
{
  vector<shared_ptr<CurveOnSurface> > trim_seg;
  shared_ptr<BoundedSurface> bd_sf;
  shared_ptr<ParamSurface> surf = face_->surface();
  for (size_t kh=0; kh<hole_idx.size(); ++kh)
    {
      int nmb_seg = 0;
      int idx1 = hole_idx[kh].first;
      int idx2 = hole_idx[kh].second;
      double len = seg_lengths[kh];

#ifdef DEBUG_REG
	  std::ofstream of0("split_segments.g2");
#endif

      if (idx1 >=0 && idx2 >=0)
	{
	  // Split between two holes
	  // For each hole, collect all corner and non-corner vertices
	  vector<shared_ptr<Vertex> > corner1;
	  vector<shared_ptr<Vertex> > corner2;
	  vector<shared_ptr<Vertex> > non_corner1;
	  vector<shared_ptr<Vertex> > non_corner2;
	  vector<ftEdge*> edges1;
	  vector<ftEdge*> edges2;
	  int nmb_holes = face_->nmbBoundaryLoops() - 1;
	  if (idx1 < nmb_holes)
	    edges1 = face_->getAllEdgePtrs(idx1+1);
	  else
	    edges1 = half_holes[idx1-nmb_holes];

	  Path::classifyCorners(edges1,  bend_,
				corner1, non_corner1);

	  if (idx2 < nmb_holes)
	    edges2 = face_->getAllEdgePtrs(idx2+1);
	  else
	    edges2 = half_holes[idx2-nmb_holes];

	  Path::classifyCorners(edges2, bend_, 
				corner2, non_corner2);

	  // Sort possible combinations of vertices according to the distance
	  // between them and whether or not one or both vertices are corners
	  vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > vx_pair;
	  vector<double> vx_dist;


	  size_t ki, kj;
	  for (ki=0; ki<corner1.size(); ++ki)
	    {
	      for (kj=0; kj<corner2.size(); ++kj)
		{
		  double dist = corner1[ki]->getDist(corner2[kj]);
		  vx_pair.push_back(make_pair(corner1[ki], corner2[kj]));
		  vx_dist.push_back(dist);
		}
	    }
	  // Sort
	  for (ki=0; ki<vx_pair.size(); ++ki)
	    for (kj=ki+1; kj<vx_pair.size(); ++kj)
	      if (vx_dist[kj] < vx_dist[ki])
		{
		  std::swap(vx_dist[ki], vx_dist[kj]);
		  std::swap(vx_pair[ki], vx_pair[kj]);
		}

	  // For each candidate pair until two candidates are selected, 
	  // check if a splitting curve between the vertices is legal
	  for (ki=0; ki<vx_pair.size(); ++ki)
	    {
	      Point parval1 = vx_pair[ki].first->getFacePar(face_.get());
	      Point parval2 = vx_pair[ki].second->getFacePar(face_.get());
	      vector<shared_ptr<CurveOnSurface> > local_seg 
		= BoundedUtils::getTrimCrvsParam(surf, parval1,
						 parval2, epsge_,
						 bd_sf);

	      // Remove segments not connected with the original vertices
	      RegularizeUtils::checkTrimSeg3(local_seg, parval1, parval2, epsge_);

	      if (local_seg.size() == 1)
		{
		  bool modified = adjustTrimSeg(local_seg[0], vx_pair[ki].first, 
						edges1, vx_pair[ki].second, 
						edges2, len);
#ifdef DEBUG_REG
		  shared_ptr<ParamCurve> cv = local_seg[0]->spaceCurve();
		  cv->writeStandardHeader(of0);
		  cv->write(of0);
#endif
	  
		  nmb_seg++;
		  trim_seg.push_back(local_seg[0]);

		  if (!modified)
		    {
		      // Remove candidates involving the current vertices
		      for (kj=ki+1; kj<vx_pair.size(); )
			{
			  if (vx_pair[kj].first.get() == vx_pair[ki].first.get() ||
			      vx_pair[kj].second.get() == vx_pair[ki].second.get())
			    {
			      vx_pair.erase(vx_pair.begin()+kj);
			      vx_dist.erase(vx_dist.begin()+kj);
			    }
			  else
			    ++kj;
			}
		      for (kj=0; kj<corner1.size(); )
			{
			  if (vx_pair[ki].first.get() == corner1[kj].get())
			    corner1.erase(corner1.begin()+kj);
			  else
			    ++kj;
			}
		      for (kj=0; kj<corner2.size(); )
			{
			  if (vx_pair[ki].second.get() == corner2[kj].get())
			    corner2.erase(corner2.begin()+kj);
			  else
			    ++kj;
			}
		    }
		}
	      if (nmb_seg == 2)
		break; 
	    }
	}

      if (nmb_seg < 2)
	{
	  // Not enough splitting curves. Try another approach
	  // Intersect each hole with a plane perpendicular to the
	  // chord between the holes a little distance from the
	  // points at the hole boundaries
	  Point p01 = seg_endpt[kh].first;
	  Point p02 = seg_endpt[kh].second;

	  // Modify points
	  int ix1 = hole_idx[kh].first;
	  int ix2 = hole_idx[kh].second;
	  int l_ix1 = (loop_idx < 0) ? ix1+1 : 
	    ((ix1 < loop_idx) ? ix1+1 : ix1+2);
	  int l_ix2 = (loop_idx < 0) ? ix2+1 : 
	    ((ix2 < loop_idx) ? ix2+1 : ix2+2);
	  shared_ptr<Loop> loop1 = face_->getBoundaryLoop(l_ix1);
	  shared_ptr<Loop> loop2 = face_->getBoundaryLoop(l_ix2);
	  int ind1, ind2;
	  double par1, par2, dist1, dist2;
	  Point p1, p2;
	  loop2->closestPoint(p01, ind2, par2, p2, dist2);
	  loop1->closestPoint(p02, ind1, par1, p1, dist1);
	  
	  // Make a selection
	  int ind3, ind4;
	  double par3, par4, dist3, dist4;
	  Point p3, p4;
	  loop1->closestPoint(p01, ind3, par3, p3, dist3);
	  loop2->closestPoint(p02, ind4, par4, p4, dist4);
	  
	  Point vec0 = p02 - p01;
	  Point vec = p2 - p1;
	  Point tan1 = loop1->getEdge(ind1)->tangent(par1);
	  Point tan2 = loop2->getEdge(ind2)->tangent(par2);
	  Point tan3 = loop1->getEdge(ind3)->tangent(par3);
	  Point tan4 = loop2->getEdge(ind4)->tangent(par4);
	  double ang1 = vec.angle(tan1);
	  double ang2 = vec.angle(tan2);
	  double ang3 = vec0.angle(tan3);
	  double ang4 = vec0.angle(tan4);
	  ang1 = fabs(0.5*M_PI - ang1);
	  ang2 = fabs(0.5*M_PI - ang2);
	  ang3 = fabs(0.5*M_PI - ang3);
	  ang4 = fabs(0.5*M_PI - ang4);
	  if (ang1 + ang2 > ang3 + ang4)
	    {
	      std::swap(p1, p01);
	      std::swap(p2, p02);
	    }

	  double r1 = (ix1 >= 0) ? holes[ix1].hole_radius_ : 0.0;
	  double r2 = (ix2 >= 0) ? holes[ix2].hole_radius_ : 0.0;
	  vec = p2 - p1;
	  double len = vec.length();
	  vec.normalize();
	  double fac = 0.15;
	  Point q1 = p1 - fac*std::min(r1, len)*vec;
	  Point q2 = p2 + fac*std::min(r2, len)*vec;
	  Point o1 = (ix1 >= 0) ? holes[ix1].hole_centre_ : p1;
	  Point o2 = (ix2 >= 0) ? holes[ix2].hole_centre_ : p2;
	  //q1 = 0.5*(q1 + o1 + (1.0-fac)*r1*vec);
	  //q2 = 0.5*(q2 + o2 - (1.0-fac)*r2*vec);
	  //Point vec1 = (q1 - o1);
	  Point vec1 = q1 - p1;
	  vec1.normalize_checked();
	  //Point vec2 = (q2 - o2);
	  Point vec2 = q2 - p2;
	  vec2.normalize_checked();
	  Point vec3 = o1 - p1; //o1 - q1;
	  vec3.normalize_checked();
	  Point vec4 = o2 - p2; //o2 - q2;
	  vec4.normalize_checked();
	  vec1 = vec3; //0.75*vec3 + 0.25*vec1; //0.25*vec3 + 0.75*vec1;
	  vec2 = vec4; //0.75*vec4 + 0.25*vec2; //0.25*vec4 + 0.75*vec2;
	  
#ifdef DEBUG_REG
	  std::ofstream of3("split_segments0.g2");
	  of3 << "400 1 0 4 155 0 100 255 " << std::endl;
	  of3 << "2 " << std::endl;
	  of3 << p01 << std::endl;
	  of3 << p02 << std::endl;
	  of3 << "400 1 0 4 255 0 0 255 " << std::endl;
	  of3 << "2 " << std::endl;
	  of3 << p1 << std::endl;
	  of3 << p2 << std::endl;
	  of3 << "400 1 0 4 0 255 0 255 " << std::endl;
	  of3 << "2 " << std::endl;
	  of3 << q1 << std::endl;
	  of3 << q2 << std::endl;
	  of3 << "400 1 0 4 0 155 100 255 " << std::endl;
	  of3 << "2 " << std::endl;
	  of3 << o1 << std::endl;
	  of3 << o2 << std::endl;
#endif
	  // Intersecting 
	  vector<shared_ptr<CurveOnSurface> > trim_seg1 =
	    BoundedUtils::getPlaneIntersections(surf, q1, vec1, epsge_, bd_sf);
	  vector<shared_ptr<CurveOnSurface> > trim_seg2 =
	    BoundedUtils::getPlaneIntersections(surf, q2, vec2, epsge_, bd_sf);

#ifdef DEBUG_REG
	  for (size_t kv=0; kv<trim_seg1.size(); ++kv)
	    {
	      shared_ptr<ParamCurve> cv = trim_seg1[kv]->spaceCurve();
	      cv->writeStandardHeader(of3);
	      cv->write(of3);
	    }
	  for (size_t kv=0; kv<trim_seg2.size(); ++kv)
	    {
	      shared_ptr<ParamCurve> cv = trim_seg2[kv]->spaceCurve();
	      cv->writeStandardHeader(of3);
	      cv->write(of3);
	    }

	      
#endif
	  // Pick appropriate endpoints of the trimming curves
	  Point cand_pt1[2], cand_pt2[2];
	  Point cand_par1[2], cand_par2[2];
	  if (trim_seg1.size() > 1)
	    extractCandPt(q1, l_ix1, trim_seg1, cand_pt1, cand_par1);
	  if (trim_seg2.size() > 1)
	    extractCandPt(q2, l_ix2, trim_seg2, cand_pt2, cand_par2);
	  if (trim_seg1.size() == 0 && trim_seg2.size() == 0)
	    continue;  // Not enough information to continue
	  Point dummy_vec;
	  if (trim_seg1.size() <= 1)
	    {
	      Point m1 = p1 + cand_pt2[0] - q2;
	      Point m2 = p1 + cand_pt2[1] - q2;
	      double upar, vpar, dist, par;
	      (void)face_->closestBoundaryPoint(m1, dummy_vec, upar, vpar, cand_pt1[0],
						dist, par);
	      cand_par1[0] = Point(upar, vpar);
	      (void)face_->closestBoundaryPoint(m2, dummy_vec, upar, vpar, cand_pt1[1],
						dist, par);
	      cand_par1[1] = Point(upar, vpar);
	    }
	  if (trim_seg2.size() <= 1)
	    {
	      Point m1 = p2 + cand_pt1[0] - q1;
	      Point m2 = p2 + cand_pt1[1] - q1;
	      double upar, vpar, dist, par;
	      (void)face_->closestBoundaryPoint(m1, dummy_vec, upar, vpar, cand_pt2[0],
						dist, par);
	      cand_par2[0] = Point(upar, vpar);
	      (void)face_->closestBoundaryPoint(m2,dummy_vec,  upar, vpar, cand_pt2[1],
						dist, par);
	      cand_par2[1] = Point(upar, vpar);
	    }

	  // Combine endpoints into pairs
	  if (cand_pt1[0].dist(cand_pt2[1]) + cand_pt1[1].dist(cand_pt2[0]) <
	      cand_pt1[0].dist(cand_pt2[0]) + cand_pt1[1].dist(cand_pt2[1]))
	    {
	      std::swap(cand_pt2[0], cand_pt2[1]);
	      std::swap(cand_par2[0], cand_par2[1]);
	    }

	  // Define one or two trimming curves
	  if (nmb_seg == 1)
	    {
	      // Select the one of the "other side" of the initial point
	    }
	  else
	    {
	      double fac = 0.3; //0.2;
	      double lim1 = fac*cand_pt1[0].dist(cand_pt1[1]);
	      double lim2 = fac*cand_pt2[0].dist(cand_pt2[1]);

	      // Snap to nearby vertices if any
	      size_t kr;
	      int idx1 = -1, idx2 = -1, idx3 = -1, idx4 = -1;
	      double mindist1 = HUGE, mindist2 = HUGE, mindist3 = HUGE, mindist4 = HUGE;
	      vector<shared_ptr<Vertex> > vx1 = loop1->getVertices();
	      removeInsignificantVertices(vx1);
	      vector<shared_ptr<Vertex> > vx2 = loop2->getVertices();
	      removeInsignificantVertices(vx2);
	      for (kr=0; kr<vx1.size(); ++kr)
		{
		  double dist1 = vx1[kr]->getVertexPoint().dist(cand_pt1[0]);
		  double dist3 = vx1[kr]->getVertexPoint().dist(cand_pt1[1]);
		  if (dist1 < mindist1 && dist1 <= lim1)
		    {
		      idx1 = (int)kr;
		      mindist1 = dist1;
		    }
		  if (dist3 < mindist3 && dist3 <= lim1)
		    {
		      idx3 = (int)kr;
		      mindist3 = dist3;
		    }
		}
	      for (kr=0; kr<vx2.size(); ++kr)
		{
		  double dist2 = vx2[kr]->getVertexPoint().dist(cand_pt2[0]);
		  double dist4 = vx2[kr]->getVertexPoint().dist(cand_pt2[1]);
		  if (dist2 < mindist2 && dist2 <= lim2)
		    {
		      idx2 = (int)kr;
		      mindist2 = dist2;
		    }
		  if (dist4 < mindist4 && dist4 <= lim2)
		    {
		      idx4 = (int)kr;
		      mindist4 = dist4;
		    }
		}

	      if (idx1 >= 0)
		{
		  cand_pt1[0] = vx1[idx1]->getVertexPoint();
		  cand_par1[0] = vx1[idx1]->getFacePar(face_.get());
		}
	      if (idx2 >= 0)
		{
		  cand_pt2[0] = vx2[idx2]->getVertexPoint();
		  cand_par2[0] = vx2[idx2]->getFacePar(face_.get());
		}
	      if (idx3 >= 0)
		{
		  cand_pt1[1] = vx1[idx3]->getVertexPoint();
		  cand_par1[1] = vx1[idx3]->getFacePar(face_.get());
		}
	      if (idx4 >= 0)
		{
		  cand_pt2[1] = vx2[idx4]->getVertexPoint();
		  cand_par2[1] = vx2[idx4]->getFacePar(face_.get());
		}
	      
	      
	      vector<shared_ptr<CurveOnSurface> > local_seg1 
		= BoundedUtils::getTrimCrvsParam(surf, cand_par1[0],
						 cand_par2[0], epsge_, bd_sf);

#ifdef DEBUG_REG
	      for (size_t kv=0; kv<local_seg1.size(); ++kv)
		{
		  shared_ptr<ParamCurve> cv = local_seg1[kv]->spaceCurve();
		  cv->writeStandardHeader(of3);
		  cv->write(of3);
		}
#endif
	      // Update segments not connected with the original vertices
	      updateTrimSeg(local_seg1, cand_pt1[0], cand_par1[0], (idx1>=0),
			    l_ix1, cand_pt2[0], cand_par2[0], (idx2>=0), l_ix2);
	      if (local_seg1.size() == 1)
		{
#ifdef DEBUG_REG
		shared_ptr<ParamCurve> cv = local_seg1[0]->spaceCurve();
		cv->writeStandardHeader(of0);
		cv->write(of0);
#endif
		
		nmb_seg++;
		trim_seg.push_back(local_seg1[0]);
		}

	      vector<shared_ptr<CurveOnSurface> > local_seg2 
		= BoundedUtils::getTrimCrvsParam(surf, cand_par1[1],
						 cand_par2[1], epsge_, bd_sf);

#ifdef DEBUG_REG
	      for (size_t kv=0; kv<local_seg2.size(); ++kv)
		{
		  shared_ptr<ParamCurve> cv = local_seg2[kv]->spaceCurve();
		  cv->writeStandardHeader(of3);
		  cv->write(of3);
		}
#endif
	      // Update segments not connected with the original vertices
	      updateTrimSeg(local_seg2, cand_pt1[1], cand_par1[1], (idx3>=0),
			    l_ix1, cand_pt2[1], cand_par2[1], (idx4>=0), l_ix2);
	      if (local_seg2.size() == 1)
		{
#ifdef DEBUG_REG
		  shared_ptr<ParamCurve> cv = local_seg2[0]->spaceCurve();
		  cv->writeStandardHeader(of0);
		  cv->write(of0);
#endif
		
		  nmb_seg++;
		  trim_seg.push_back(local_seg2[0]);
		}
	    }
	}
    }
  
  vector<shared_ptr<ftSurface> > faces;
  if (trim_seg.size() < 2)
    return faces;   // Trimming curves will not split
      
  // Divide surface
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_sf, trim_seg, epsge_);

#ifdef DEBUG_REG
  std::ofstream of("split_surf.g2");
  for (size_t kr=0; kr<sub_sfs.size(); ++kr)
    {
      sub_sfs[kr]->writeStandardHeader(of);
      sub_sfs[kr]->write(of);
    }
#endif
      
  // Fetch info about vertices not belonging to the corners of the
  // initial face
  vector<shared_ptr<Vertex> > non_corner = 
    face_->getNonCornerVertices(bend_);
  removeInsignificantVertices(non_corner);

  faces = RegularizeUtils::createFaces(sub_sfs, face_, epsge_, tol2_,
				       angtol_, non_corner);
 
  return faces;
}

//==========================================================================
void RegularizeFace::updateTrimSeg(std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
				   const Point& pnt1, const Point& param1, 
				   bool at_vx1, int loop_idx1, 
				   const Point& pnt2, const Point& param2, 
				   bool at_vx2, int loop_idx2)
//==========================================================================
{
  shared_ptr<ParamSurface> surf = face_->surface();
  shared_ptr<Loop> loop1 = face_->getBoundaryLoop(loop_idx1);
  shared_ptr<Loop> loop2 = face_->getBoundaryLoop(loop_idx2);
  shared_ptr<BoundedSurface> bd_sf;
  vector<shared_ptr<Vertex> > vx1 = loop1->getVertices();
  vector<shared_ptr<Vertex> > vx1_2 = vx1;
  removeInsignificantVertices(vx1);
  vector<shared_ptr<Vertex> > vx2 = loop2->getVertices();
  vector<shared_ptr<Vertex> > vx2_2 = vx2;
  removeInsignificantVertices(vx2);
  double fac = 0.3; //0.4;
  for (size_t kr=0; kr<trim_segments.size(); ++kr)
    {
      Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
      Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());
      Point parval1, parval2;
      if (pnt1.dist(pos1) > epsge_ && (!at_vx1))
	{
	  int ind1, ind2;
	  double par1, par2, dist1, dist2;
	  Point p1, p2;
	  loop1->closestPoint(pnt1, ind1, par1, p1, dist1);
	  loop1->closestPoint(pos1, ind2, par2, p2, dist2);
	  if (ind1 == ind2)
	    {
	      double par = fac*par1 + (1-fac)*par2;
	      parval1 = loop1->getEdge(ind1)->geomEdge()->faceParameter(par);
	      pos1 = face_->point(parval1[0], parval1[1]);
	    }
	  else if (fabs(ind2-ind1) == 1)
	    {
	      ftEdge* e1 = loop1->getEdge(ind1)->geomEdge();
	      ftEdge* e2 = loop1->getEdge(ind2)->geomEdge();
	      shared_ptr<Vertex> v1 = e1->getCommonVertex(e2);
	      if (loop1->size() == 2)
		{
		  shared_ptr<Vertex> v2 = e1->getOtherVertex(v1.get());
		  if (v2->getVertexPoint().dist(pos1) <
		      v1->getVertexPoint().dist(pos1))
		    v1 = v2;
		}
	      pos1 = v1->getVertexPoint();
	      parval1 = v1->getFacePar(face_.get());
	    }
	  else
	    {
	      int ind3;
	      double par3, dist3;
	      Point p3;
	      Point mid = fac*p1 + (1-fac)*p2;
	      loop1->closestPoint(mid, ind3, par3, p3, dist3);
	      parval1 = loop1->getEdge(ind3)->geomEdge()->faceParameter(par3);
	      pos1 = p3;
	    }
	}
      else
	{
	  pos1 = pnt1;
	  parval1 = param1;
	}
	  
      if (pnt2.dist(pos2) > epsge_ && (!at_vx2))
	{
	  int ind1, ind2;
	  double par1, par2, dist1, dist2;
	  Point p1, p2;
	  loop2->closestPoint(pnt2, ind1, par1, p1, dist1);
	  loop2->closestPoint(pos2, ind2, par2, p2, dist2);
	  if (ind1 == ind2)
	    {
	      double par = fac*par1 + (1-fac)*par2;
	      parval2 = loop2->getEdge(ind1)->geomEdge()->faceParameter(par);
	      pos2 = face_->point(parval2[0], parval2[1]);
	    }
	  else if (fabs(ind2-ind1) == 1)
	    {
	      ftEdge* e1 = loop2->getEdge(ind1)->geomEdge();
	      ftEdge* e2 = loop2->getEdge(ind2)->geomEdge();
	      shared_ptr<Vertex> v1 = e1->getCommonVertex(e2);
	      if (loop2->size() == 2)
		{
		  shared_ptr<Vertex> v2 = e2->getOtherVertex(v1.get());
		  if (v2->getVertexPoint().dist(pos2) <
		      v1->getVertexPoint().dist(pos2))
		    v1 = v2;
		}
	      pos2 = v1->getVertexPoint();
	      parval2 = v1->getFacePar(face_.get());
	    }
	  else
	    {
	      int ind3;
	      double par3, dist3;
	      Point p3;
	      Point mid = fac*p1 + (1-fac)*p2;
	      loop2->closestPoint(mid, ind3, par3, p3, dist3);
	      parval2 = loop2->getEdge(ind3)->geomEdge()->faceParameter(par3);
	      pos2 = p3;
	    }
	}
      else
	{
	  pos2 = pnt2;
	  parval2 = param2;
	}

      if (parval1.dimension() > 0 || parval2.dimension() > 0)
	{
	  // Update
	  int idx1 = -1, idx2 = -1;
	  double lim = 0.3*pos1.dist(pos2);
	  double lim2 = 2.0*model_->getTolerances().neighbour;
	  if (parval1.dimension() == 0 || at_vx1)
	    parval1 = param1;
	  else
	    {
	      // Check with existing vertices
	      double mindist = HUGE;
	      for (size_t kh=0; kh<vx1.size(); ++kh)
		{
		  double dist = vx1[kh]->getVertexPoint().dist(pos1);
		  if (dist < mindist && dist <= lim)
		    {
		      idx1 = (int)kh;
		      mindist = dist;
		    }
		}
	      if (idx1 < 0)
		{
		  // Check all vertices
		  for (size_t kh=0; kh<vx1_2.size(); ++kh)
		    {
		      double dist = vx1_2[kh]->getVertexPoint().dist(pos1);
		      if (dist < mindist && dist <= lim2)
			{
			  idx1 = (int)kh;
			  mindist = dist;
			}
		    }
		  if (idx1 >= 0)
		    {
		      pos1 = vx1_2[idx1]->getVertexPoint();
		      parval1 = vx1_2[idx1]->getFacePar(face_.get());
		      idx1 = -1;
		    }
		}
	      if (idx1 >= 0)
		parval1 = vx1[idx1]->getFacePar(face_.get());
	    }
	  if (parval2.dimension() == 0 || at_vx2)
	    parval2 = param2;
	  else
	    {
	      // Check with existing vertices
	      double mindist = HUGE;
	      for (size_t kh=0; kh<vx2.size(); ++kh)
		{
		  double dist = vx2[kh]->getVertexPoint().dist(pos2);
		  if (dist < mindist && dist <= lim)
		    {
		      idx2 = (int)kh;
		      mindist = dist;
		    }
		}
	      if (idx2 < 0)
		{
		  // Check all vertices
		  for (size_t kh=0; kh<vx2.size(); ++kh)
		    {
		      double dist = vx2[kh]->getVertexPoint().dist(pos2);
		      if (dist < mindist && dist <= lim2)
			{
			  idx2 = (int)kh;
			  mindist = dist;
			}
		    }
		  if (idx2 >= 0)
		    {
		      pos2 = vx2_2[idx1]->getVertexPoint();
		      parval2 = vx2_2[idx2]->getFacePar(face_.get());
		      idx2 = -1;
		    }
		}
	      if (idx2 >= 0)
		parval2 = vx2[idx2]->getFacePar(face_.get());
	    }
	  
	  shared_ptr<ParamCurve> pcrv;
	  if (idx1 >=0 && idx2 >= 0)
	    pcrv = RegularizeUtils::checkStrightParCv(face_, vx1[idx1], vx2[idx2],
						      epsge_);
	  else if (idx1 >= 0)
	    pcrv = RegularizeUtils::checkStrightParCv(face_, vx1[idx1], pnt2,
						      epsge_);
	  else if (idx2 >= 0)
	    pcrv = RegularizeUtils::checkStrightParCv(face_, vx2[idx2], pnt1,
						      epsge_);
	  else
	    pcrv = RegularizeUtils::checkStrightParCv(face_, pos1, pos2,
						      epsge_);

	  vector<shared_ptr<CurveOnSurface> > local_seg;
	  if (pcrv.get())
	    local_seg = BoundedUtils::getTrimCrvsPcrv(surf, pcrv, epsge_, bd_sf);
	  else
	    local_seg = BoundedUtils::getTrimCrvsParam(surf, parval1,
						       parval2, epsge_, bd_sf);
	  if (local_seg.size() == 1)
	    trim_segments[kr] = local_seg[0];
	}
    }
}

//==========================================================================
void RegularizeFace::extractCandPt(Point mid, int hole_ix,
				   vector<shared_ptr<CurveOnSurface> >& seg,
				   Point cand_pt[], Point cand_par[])
//==========================================================================
{
  // For each segment, select the endpoint closest to the current hole boundary
  shared_ptr<Loop> loop = face_->getBoundaryLoop(hole_ix);
  vector<Point> seg_pt(seg.size());
  vector<Point> seg_par(seg.size());
  for (size_t ki=0; ki<seg.size(); ++ki)
    {
      Point pos1 = seg[ki]->ParamCurve::point(seg[ki]->startparam());
      Point pos2 = seg[ki]->ParamCurve::point(seg[ki]->endparam());

      int ix1, ix2;
      double par1, par2, dist1, dist2;
      Point clo1, clo2;
      loop->closestPoint(pos1, ix1, par1, clo1, dist1);
      loop->closestPoint(pos2, ix2, par2, clo2, dist2);
      seg_pt[ki] = (dist1 <= dist2) ? pos1 : pos2;
      seg_par[ki] = (dist1 <= dist2) ?
	seg[ki]->parameterCurve()->point(seg[ki]->startparam()) :
	seg[ki]->parameterCurve()->point(seg[ki]->endparam());
    }

  if (seg.size() == 2)
    {
      for (int kr=0; kr<2; ++kr)
	{
	  cand_pt[kr] = seg_pt[kr];
	  cand_par[kr] = seg_par[kr];
	}
    }
  else
    {
      // Select the two closest point to the given mid point
      // First sort
      for (size_t ki=0; ki<seg_pt.size(); ++ki)
	for (size_t kj=ki+1; kj<seg_pt.size(); ++kj)
	  {
	    if (seg_pt[kj].dist(mid) < seg_pt[ki].dist(mid))
	      {
		std::swap(seg_pt[ki], seg_pt[kj]);
		std::swap(seg_par[ki], seg_par[kj]);
	      }
	  }
      for (int kr=0; kr<2; ++kr)
	{
	  cand_pt[kr] = seg_pt[kr];
	  cand_par[kr] = seg_par[kr];
	}
    }
}

//==========================================================================
bool RegularizeFace::adjustTrimSeg(shared_ptr<CurveOnSurface>& trim_seg,
				   shared_ptr<Vertex> vx1,
				   const vector<ftEdge*>& edges1,
				   shared_ptr<Vertex> vx2,
				   const vector<ftEdge*>& edges2,
				   double len)
//==========================================================================
{
  // Compute minimum angle between new segment and existing edges in
  // both endpoints
  double min_ang1, min_ang2;
  RegularizeUtils::angleInEndpoints(trim_seg, vx1, vx2, face_,
				    min_ang1, min_ang2);
  Point pos1 = vx1->getVertexPoint();
  Point pos2 = vx2->getVertexPoint();

  // Compute segment length
  double seg_len = trim_seg->estimatedCurveLength();

  // Check if the curve segment should be modified
  double fac = 10.0;
  if (min_ang1 < bend_ || min_ang2 < bend_ ||
      seg_len > fac*len)
    {
      // Compute closest points to the other edge sequence from both 
      // vertices and also from the found closest point to this edge
      // sequence
      int ind1, ind2, ind3, ind4;
      double par1, par2, par3, par4, dist1, dist2, dist3, dist4;
      Point pt1, pt2, pt3, pt4;
      Path::closestPoint(edges2, pos1, ind1, par1, pt1, dist1);
      Path::closestPoint(edges1, pt1, ind2, par2, pt2, dist2);
      Path::closestPoint(edges1, pos2, ind3, par3, pt3, dist3);
      Path::closestPoint(edges2, pt3, ind4, par4, pt4, dist4);

      // Select candidate (pt1-pt2 or pt3-pt4)
      Point vec1 = pt2 - pt1;
      Point vec2 = pt4 - pt3;
      Point tan1 = edges2[ind1]->tangent(par1);
      Point tan2 = edges1[ind2]->tangent(par2);
      Point tan3 = edges1[ind3]->tangent(par3);
      Point tan4 = edges2[ind4]->tangent(par4);
      double ang1 = vec1.angle(tan1);
      ang1 = std::min(ang1, fabs(M_PI-ang1));
      double ang2 = vec1.angle(tan2);
      ang2 = std::min(ang2, fabs(M_PI-ang2));
      double ang3 = vec2.angle(tan3);
      ang3 = std::min(ang3, fabs(M_PI-ang3));
      double ang4 = vec2.angle(tan4);
      ang4 = std::min(ang4, fabs(M_PI-ang4));

      // Select the candidate where the new distance is less than
      // the threshold (or smallest) and if both, the one that enables
      // another trimming curve to the same point
      int config = 0;
      if (dist2 < fac*len && dist4 < fac*len)
	{
	  if (std::min(ang1, ang2) <= bend_)
	    config = 2;
	  else if (std::min(ang3, ang4) <= bend_)
	    config = 1;
	  else
	    {
	      double ang12 = std::max(std::min(ang1, fabs(0.5*M_PI-ang1)),
				      std::min(ang2, fabs(0.5*M_PI-ang2)));
	      double ang34 = std::max(std::min(ang3, fabs(0.5*M_PI-ang3)),
				      std::min(ang4, fabs(0.5*M_PI-ang4)));
	      if (ang12 > ang34)
		config = 1;
	      else
		config = 2;
	    }
	}
      else if (dist2 < fac*len || dist2 < dist4)
	config = 1;
      else 
	config = 2;
				  
      if (config == 0)
	return false;
      
      Point parval1, parval2;
      if (config == 1)
	{
	  parval1 = edges2[ind1]->faceParameter(par1);
	  parval2 = edges1[ind2]->faceParameter(par2);
	}
      else
	{
	  parval1 = edges1[ind3]->faceParameter(par3);
	  parval2 = edges2[ind4]->faceParameter(par4);
	}

      shared_ptr<BoundedSurface> bd_sf;
      shared_ptr<ParamSurface> surf = face_->surface();
      vector<shared_ptr<CurveOnSurface> > new_seg =
	 BoundedUtils::getTrimCrvsParam(surf, parval1, parval2, epsge_,
					bd_sf);
      RegularizeUtils::checkTrimSeg3(new_seg, parval1, parval2, epsge_);
      if (new_seg.size() == 1)
	{
	  trim_seg = new_seg[0];
	  return true;
	}
      else
	return false;
    }
  return false;
}

//==========================================================================
int
RegularizeFace::positionWeigthPoint(const Point& wgt_par)
//==========================================================================
{
  // Get surface
  shared_ptr<ParamSurface> surf = face_->surface();
  shared_ptr<BoundedSurface> bd_surf =
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
  if (!bd_surf.get())
    {
      bool inside = surf->inDomain(wgt_par[0], wgt_par[1]);
      if (inside)
	return 0;
      else 
	return -1;
    }
  else
    {
      // Fetch curve bounded domain corresponding to the surface
      CurveBoundedDomain domain = bd_surf->parameterDomain();

      // Check if the point lies at the boundary
      Array<double,2> par(wgt_par[0], wgt_par[1]);
      if (domain.isOnBoundary(par, epsge_))
	return 0;  // The boundary is defined to be inside
      else
	{
	  int domain_pos;
	  int ki;
	  for (ki=0; ki<3; ++ki)
	    {
	      domain_pos = -2;
	      try {
		domain_pos = domain.positionPointInDomain(ki, wgt_par[0],
							  wgt_par[1], epsge_);
	      }
	      catch (...)
		{
		}

	      if (domain_pos > -2)
		break;
	    }
	  return domain_pos;
	}
    }
}

//==========================================================================
void RegularizeFace::removeHalfHoleVx(vector<shared_ptr<Vertex> >& vx,
				      vector<vector<ftEdge*> >& half_holes)
//==========================================================================
{
  size_t ki;
  for (ki=0; ki<half_holes.size(); ++ki)
    {
      shared_ptr<Vertex> curr_vx = half_holes[ki][0]->getVertex(true);
      vector<shared_ptr<Vertex> >::iterator vxp =
	std::find(vx.begin(), vx.end(), curr_vx);
      if (vxp != vx.end())
	vx.erase(vxp);
      for (size_t kj=0; kj<half_holes[ki].size(); ++kj)
	{
	  curr_vx = half_holes[ki][kj]->getVertex(false);
	  vxp = std::find(vx.begin(), vx.end(), curr_vx);
	  if (vxp != vx.end())
	    vx.erase(vxp);
	}
    }
}

//==========================================================================
ftSurface*
RegularizeFace::identifySeamFaces(shared_ptr<ftSurface> face1,
				  vector<shared_ptr<ftSurface> > faces,
				  int& pardir)
//==========================================================================
{
  ftSurface* face2;

  // Fetch all edges
  vector<shared_ptr<ftEdge> > edges = face1->getAllEdges();

  // Fetch info about all non-smooth transitions
  FaceConnectivityUtils<ftEdgeBase,ftFaceBase> connectivity;
  vector<ftEdgeBase*> kinks;
  vector<shared_ptr<ftFaceBase> > tmp(faces.begin(), faces.end());
  connectivity.cornersAndKinks(tmp, kinks);

  // Look for smooth transitions where the underlying surface is
  // the same for both faces
  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      if (!edges[ki]->twin())
	continue;  // No adjacent face

      size_t kj;
      for (kj=0; kj<kinks.size(); ++kj)
	if (edges[ki].get() == kinks[kj] ||
	    edges[ki].get() == kinks[kj]->twin())
	  break;
      if (kj < kinks.size())
	continue;  // Not a smooth transition

      face2 = edges[ki]->twin()->geomEdge()->face()->asFtSurface();
      if (!face2)
	continue;  
      
      // Fetch information about trimming curves
      AdjacencyInfo info = face1->getAdjacencyInfo(face2, epsge_);
      if (info.bd_idx_1_ < 0 || info.bd_idx_2_ < 0)
	continue;   // Both edges do not follow a boundary curve

      if (abs(info.bd_idx_1_ - info.bd_idx_2_) != 1)
	continue;  // Not min and max parameter
      if ((info.bd_idx_1_/2) != (info.bd_idx_2_/2))
	continue;  // Not the same parameter direction

      // Check if more than one adjacency instance exists for the same two faces
      shared_ptr<ftEdge> edg1, edg2;
      bool neighbours = face1->areNeighbours(face2, edg1, edg2, 1);
      if (neighbours)
	continue;

      pardir = (info.bd_idx_1_ <= 1) ? 0 : 1;

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

      // Check if the next or previous edge follow the same constant boundary
      int dir1 = (info.bd_idx_1_ <= 1) ? 1 : 2;
      int dir2 = (info.bd_idx_2_ <= 1) ? 1 : 2;
      int dir;
      double val;
      shared_ptr<ParamCurve> tmp_crv = edges[ki]->next()->geomEdge()->geomCurve();
      shared_ptr<CurveOnSurface> sf_crv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
      if (sf_crv.get())
	sf_crv->isConstantCurve(epsge_, dir, val);
      else
	dir = -1;
      if (dir == dir1)
	continue;
      tmp_crv = edges[ki]->prev()->geomEdge()->geomCurve();
      sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
     if (sf_crv.get())
	sf_crv->isConstantCurve(epsge_, dir, val);
      else
	dir = -1;
      if (dir == dir1)
	continue;
      tmp_crv = edges[ki]->twin()->next()->geomEdge()->geomCurve();
      sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
     if (sf_crv.get())
	sf_crv->isConstantCurve(epsge_, dir, val);
      else
	dir = -1;
      if (dir == dir2)
	continue;
      tmp_crv = edges[ki]->twin()->prev()->geomEdge()->geomCurve();
      sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
     if (sf_crv.get())
	sf_crv->isConstantCurve(epsge_, dir, val);
      else
	dir = -1;
      if (dir == dir2)
	continue;

      // Check if the joining curves are smooth
      shared_ptr<Vertex> v1, v2;
      edges[ki]->getVertices(v1, v2);
      ftEdge *e1, *e2, *e3, *e4;
      if (v1->hasEdge(edges[ki]->next()->geomEdge()))
	{
	  e1 = edges[ki]->next()->geomEdge();
	  e3 = edges[ki]->prev()->geomEdge();
	}
      else
	{
	  e1 = edges[ki]->prev()->geomEdge();
	  e3 = edges[ki]->next()->geomEdge();
	}
      if (v1->hasEdge(edges[ki]->twin()->next()->geomEdge()))
 	{
	  e2 = edges[ki]->twin()->next()->geomEdge();
	  e4 = edges[ki]->twin()->prev()->geomEdge();
	}
      else
	{
	  e2 = edges[ki]->twin()->prev()->geomEdge();
	  e4 = edges[ki]->twin()->next()->geomEdge();
	}
      double t1 = e1->parAtVertex(v1.get());
      double t2 = e2->parAtVertex(v1.get());
      Point tan1 = e1->tangent(t1);
      Point tan2 = e2->tangent(t2);
      if (tan1.angle(tan2) > angtol_)
	continue;  // Not a smooth transition

      t1 = e3->parAtVertex(v2.get());
      t2 = e4->parAtVertex(v2.get());
      tan1 = e3->tangent(t1);
      tan2 = e4->tangent(t2);
      if (tan1.angle(tan2) > angtol_)
	continue;
      

      // It remains to check if the underlying surface is closed
      // in the relevant parameter direction

      return face2;
    }

  return 0;  // No match found

}

//==========================================================================
shared_ptr<ftSurface>
RegularizeFace::mergeSeamFaces(ftSurface* face1, ftSurface* face2, int pardir)
//==========================================================================
{
  // Get surfaces and check consistency
  shared_ptr<ftSurface> dummy;
  shared_ptr<ParamSurface> surf1 = face1->surface();
  shared_ptr<ParamSurface> surf2 = face2->surface();
  shared_ptr<BoundedSurface> bd_sf1 = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
  shared_ptr<BoundedSurface> bd_sf2 = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
  if (!(bd_sf1.get() && bd_sf2.get()))
    return dummy;

  shared_ptr<ParamSurface> base = bd_sf1->underlyingSurface();
  if (base.get() != bd_sf2->underlyingSurface().get())
    return dummy;  // Different underlying surfaces

  // Get domains
  RectDomain dom1 = bd_sf1->containingDomain();
  RectDomain dom2 = bd_sf2->containingDomain();
  RectDomain dom3 = base->containingDomain();

  // Set common domain, including periodicity
  double umin1, umin2, umax1, umax2, vmin1, vmin2, vmax1, vmax2;
  bool move_first;
  double udel, vdel;
  if (pardir == 0)
    {
      umin1 = std::max(dom3.umin(),dom1.umin());
      umin2 = std::max(dom3.umin(),dom2.umin());
      umax1 = std::min(dom3.umax(),dom1.umax());
      umax2 = std::min(dom3.umax(),dom2.umax());
      vmin1 = vmin2 = dom3.vmin();
      vmax1 = vmax2 = dom3.vmax();
      move_first = (dom1.umin() < dom2.umin());
      udel = (move_first) ? umax2 - umin1 : umax1 - umin2;
      vdel = 0.0;
    }
  else
    {
      umin1 = umin2 = dom3.umin();
      umax1 = umax2 = dom3.umax();
      vmin1 = std::max(dom3.vmin(),dom1.vmin());
      vmin2 = std::max(dom3.vmin(),dom2.vmin());
      vmax1 = std::min(dom3.vmax(),dom1.vmax());
      vmax2 = std::min(dom3.vmax(),dom2.vmax());
      move_first = (dom1.vmin() < dom2.vmin());
      udel = 0.0;
      vdel = (move_first) ? vmax2 - vmin1 : vmax1 - vmin2;
    }

  // Make new underlying surface
  SplineSurface *base2 = base->asSplineSurface();
  if (!base2)
    return dummy;
  shared_ptr<SplineSurface> sub1 = 
    shared_ptr<SplineSurface>(base2->subSurface(umin1, vmin1, umax1, vmax1));
  shared_ptr<SplineSurface> sub2 = 
    shared_ptr<SplineSurface>(base2->subSurface(umin2, vmin2, umax2, vmax2));
  if (!move_first)
    std::swap(sub1, sub2);

#ifdef DEBUG_REG
  std::ofstream of("merge_sfs.g2");
  sub1->writeStandardHeader(of);
  sub1->write(of);
  sub2->writeStandardHeader(of);
  sub2->write(of);
#endif

  double dist;
  sub2->appendSurface(sub1.get(), pardir+1, 1, dist, false);

#ifdef DEBUG_REG
  sub2->writeStandardHeader(of);
  sub2->write(of);
#endif
  if (dist > epsge_)
    return dummy;

  // Make boundary curves of merged surface.
  // First adjust parameter curves to fit with the modified domain

  // Fetch boundary loops
  CurveLoop loop1 = bd_sf1->outerBoundaryLoop();
    bd_sf1->outerBoundaryLoop();
  CurveLoop loop2 = bd_sf2->outerBoundaryLoop();

  // Remove coincident curves
  // DEBUG. Draw
#ifdef DEBUG_REG
  std::ofstream space("space.g2");
  sub2->writeStandardHeader(space);
  sub2->write(space);
#endif
  int nmb1 = loop1.size();
  int nmb2 = loop2.size();
  int ki, kj;
  vector<shared_ptr<CurveOnSurface> > bd_cvs;
  Point vec(udel, vdel);
  for (ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamCurve> cv = loop1[ki];
      shared_ptr<CurveOnSurface> sf_cv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
      if (!sf_cv.get())
	continue;
#ifdef DEBUG_REG
      sf_cv->spaceCurve()->writeStandardHeader(space);
      sf_cv->spaceCurve()->write(space);
#endif
      shared_ptr<CurveOnSurface> tmp = 
	shared_ptr<CurveOnSurface>(sf_cv->clone());
      tmp->setUnderlyingSurface(sub2);
      if (move_first)
	{
	  bool changed = tmp->translateParameterCurve(vec);
	  if (!changed)
	    {
	      shared_ptr<ParamCurve> par_cv = tmp->parameterCurve();
	      Point par1 = par_cv->point(par_cv->startparam());
	      par1 += vec;
	      Point par2 = par_cv->point(par_cv->endparam());
	      par2 += vec;
	      tmp->makeParameterCurve(epsge_, par1, par2);
	    }
	}
      bd_cvs.push_back(tmp);
    }

  for (ki=0; ki<nmb2; ++ki)
    {
      shared_ptr<ParamCurve> cv = loop2[ki];
      shared_ptr<CurveOnSurface> sf_cv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
      if (!sf_cv.get())
	continue;
#ifdef DEBUG_REG
      sf_cv->spaceCurve()->writeStandardHeader(space);
      sf_cv->spaceCurve()->write(space);
#endif
      shared_ptr<CurveOnSurface> tmp = 
	shared_ptr<CurveOnSurface>(sf_cv->clone());
      tmp->setUnderlyingSurface(sub2);
      if (!move_first)
	{
	  bool changed = tmp->translateParameterCurve(vec);
	  if (!changed)
	    {
	      shared_ptr<ParamCurve> par_cv = tmp->parameterCurve();
	      Point par1 = par_cv->point(par_cv->startparam());
	      par1 += vec;
	      Point par2 = par_cv->point(par_cv->endparam());
	      par2 += vec;
	      tmp->makeParameterCurve(epsge_, par1, par2);
	    }
	}
      bd_cvs.push_back(tmp);
    }

  // Combine loops, first remove duplicates
  bool found = false;
  for (ki=0; ki<nmb1;)
    {
      Point pos1 = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->startparam());
      Point pos2 = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->endparam());
      for (kj=nmb1; kj<nmb1+nmb2;)
	{
	  found = false;
	  Point pos3 = bd_cvs[kj]->ParamCurve::point(bd_cvs[kj]->startparam());
	  Point pos4 = bd_cvs[kj]->ParamCurve::point(bd_cvs[kj]->endparam());
	  if ((pos1.dist(pos4) < epsge_ && pos2.dist(pos3) < epsge_) ||
	      (pos1.dist(pos3) < epsge_ && pos2.dist(pos4) < epsge_))
	    {
	      // Maybe some more checking?
	      bd_cvs.erase(bd_cvs.begin()+kj);
	      nmb2--;
	      bd_cvs.erase(bd_cvs.begin()+ki);
	      nmb1--;
	      found = true;
	      break;
	    }
	  kj++;
	}
      if (!found)
	ki++;
    }

  // Sort curves
  for (ki=0; ki<(int)bd_cvs.size()-1; ++ki)
    {
      Point pos = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->endparam());
      Point pos1 = bd_cvs[ki+1]->ParamCurve::point(bd_cvs[ki+1]->startparam());
      Point pos2 = bd_cvs[ki+1]->ParamCurve::point(bd_cvs[ki+1]->endparam());
      double dist1 = pos.dist(pos1);
      double dist2 = pos.dist(pos2);
      for (kj=ki+2; kj<(int)bd_cvs.size(); ++kj)
	{
	  Point pos3 = bd_cvs[kj]->ParamCurve::point(bd_cvs[kj]->startparam());
	  Point pos4 = bd_cvs[kj]->ParamCurve::point(bd_cvs[kj]->endparam());
	  double dist3 = pos.dist(pos3);
	  double dist4 = pos.dist(pos4);
	  if (std::min(dist3,dist4) < std::min(dist1,dist2))
	    {
	      if (dist4 < dist3 && dist4  < tol2_)
		{
		  bd_cvs[kj]->reverseParameterDirection();
		  std::swap(dist3, dist4);
		}
	      std::swap(bd_cvs[ki+1], bd_cvs[kj]);
	      dist1 = dist3;
	      dist2 = dist4;
	    }
	}
    }

#ifdef DEBUG_REG
  std::ofstream space2("space2.g2");
  std::ofstream par2("par2.g2");
  for (ki=0; ki<(int)bd_cvs.size(); ++ki)
    {
      bd_cvs[ki]->spaceCurve()->writeStandardHeader(space2);
      bd_cvs[ki]->spaceCurve()->write(space2);
      bd_cvs[ki]->parameterCurve()->writeStandardHeader(par2);
      bd_cvs[ki]->parameterCurve()->write(par2);
    }
#endif

  // Fetch joint positions in boundary loop
  vector<Point> joints(bd_cvs.size());
  for (ki=0; ki<(int)bd_cvs.size(); ++ki)
    joints[ki] = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->startparam());

  // Make new bounded surface
  shared_ptr<BoundedSurface> merged = 
    shared_ptr<BoundedSurface>(new BoundedSurface(sub2, bd_cvs, epsge_));
  merged->analyzeLoops();
  double merge_dist;
  bool success;
  success = merged->simplifyBdLoops(epsge_, angtol_, merge_dist);
  merged->analyzeLoops();

#ifdef DEBUG_REG
  std::ofstream merge("merge_sf.g2");
  merged->writeStandardHeader(merge);
  merged->write(merge);
#endif

  // Make face
  shared_ptr<ftSurface> merged_face = 
    shared_ptr<ftSurface>(new ftSurface(merged, -1));
  merged_face->setBody(face1->getBody());
  (void)merged_face->createInitialEdges(epsge_, angtol_);

  // Check if any joints has been removed. Fetch face vertices
  vector<shared_ptr<Vertex> > vx = merged_face->vertices();
  for (ki=0; ki<(int)joints.size();)
    {
      for (kj=0; kj<(int)vx.size(); ++kj)
	if (joints[ki].dist(vx[kj]->getVertexPoint()) < epsge_)
	  break;
      if (kj < (int)vx.size())
	joints.erase(joints.begin()+ki);
      else
	ki++;
    }
  // Store info about removed joints
  seam_joints_.insert(seam_joints_.end(), joints.begin(), joints.end());

  return merged_face;
  
}

//==========================================================================
bool
RegularizeFace::fetchPatternSplit(const Point& corner,
				  Point& parval1, Point& parval2,
				  bool use_input_point)
//==========================================================================
{
  bool found = false;

  // Compute parameter value in corner
  double vx_upar, vx_vpar, vx_dist, vx_t;
  Point vx_pos;
  Point dummy_vec;
  ftEdgeBase *tmp_edge = face_->closestBoundaryPoint(corner, dummy_vec, vx_upar, vx_vpar,
						     vx_pos, vx_dist, vx_t);
  parval1 = Point(vx_upar, vx_vpar);

  // First check if there is an endpoint of previous splits close
  // to the given corner
  size_t ki;
  double fac = 0.1;
  double fac2 = 0.5;
  int idx = -1;
  size_t min_idx = -1;
  double min_dist = 1.0e8;
  double cv_len = 0.0;
  shared_ptr<ParamSurface> surf = face_->surface();
  Point cv_pos, other_pos;
  //const Domain& dom = surf->parameterDomain();
  for (ki=0; ki<cand_split_.size(); ++ki)
    {
      // Fetch endpoints in this face
      Point pos1, pos2;
      double u1, v1, u2, v2, d1, d2;
      face_->closestPoint(cand_split_[ki].first.first, u1, v1, pos1, d1, epsge_);
      face_->closestPoint(cand_split_[ki].second.first, u2, v2, pos2, d2, epsge_);

      double len = pos1.dist(pos2);
      double dist1 = corner.dist(pos1);
      double dist2 = corner.dist(pos2);
      ftEdgeBase *edge1 = face_->edgeClosestToPoint(u1, v1);
      ftEdgeBase *edge2 = face_->edgeClosestToPoint(u2, v2);
      if (edge1 == edge2 || edge1->next() == edge2 ||
	  edge1->prev() == edge2)
	continue;  // Candidate split curve ends up in same or adjacent edge
      shared_ptr<Vertex> close1 = face_->getClosestVertex(pos1);
      shared_ptr<Vertex> close2 = face_->getClosestVertex(pos2);
      double vx_d1 = close1->getVertexPoint().dist(pos1);
      double vx_d2 = close2->getVertexPoint().dist(pos2);
      if (close1->sameEdge(close2.get()) && vx_d1 < fac*len && vx_d2 < fac*len)
	  continue;  // Candidate split curve ends up in same edge

      if ((dist1 < fac*len || (dist1 < fac2*len && tmp_edge == edge1)) && 
	  dist1 < min_dist)
	{
	  idx = 1;
	  min_dist = dist1;
	  min_idx = ki;
	  cv_len = len;
	  cv_pos = pos2;
	  other_pos = pos1;
	}
      if ((dist2 < fac*len || (dist2 < fac2*len && tmp_edge == edge2)) && 
	  dist2 < min_dist)
	{
	  idx = 2;
	  min_dist = dist2;
	  min_idx = ki;
	  cv_len = len;
	  cv_pos = pos1;
	  other_pos = pos2;
	}
    }

  if (idx < 0)
    return false;  // No candidat split is found
      
      
  // Find closest boundary point to opposite end of selected split
  double par_u, par_v, dist, par_t;
  Point bd_pos;
  ftEdgeBase *edge = face_->closestBoundaryPoint(cv_pos, dummy_vec, par_u, par_v,
						 bd_pos, dist, par_t);
  if (bd_pos.dist(cv_pos) < fac*cv_len)
    {
      parval2 = Point(par_u, par_v);

      // // Extend the curve to ensure a proper split
      // Point vec = parval2 - parval1;
      // parval2 += 0.05*vec;
      
      double u1, v1, d1;
      Point face_pos;
      face_->closestPoint(bd_pos, u1, v1, face_pos, d1, epsge_);
      Point tmp1 = face_->point(par_u, par_v);
      Point tmp2 = edge->point(par_t);
      Point seed = parval2;
      parval2 = edge->geomEdge()->faceParameter(par_t, seed.begin());
      found = true;
    }

  if (!use_input_point)
    {
      parval2 = Point(par_u, par_v);
      (void)face_->closestBoundaryPoint(other_pos, dummy_vec, par_u, par_v,
					bd_pos, dist, par_t);
      if (bd_pos.dist(other_pos) < fac*cv_len)
	{
	  parval1 = Point(par_u, par_v);
	  found = true;
	}
      else
	found = false;
    }
  return found;
}

//==========================================================================
bool
RegularizeFace::fetchPatternSplit2(const Point& cand_pt,
				   Point& parval1, Point& parval2,
				   double level_dist, double level_ang,
				   const Point& normal, bool only_outer_bd)
//==========================================================================
{
  // Input point
  double vx_upar, vx_vpar, vx_dist, vx_t;
  Point vx_pos;
  Point dummy_vec;
  ftEdgeBase *vx_edg = face_->closestBoundaryPoint(cand_pt, dummy_vec, vx_upar, vx_vpar,
						   vx_pos, vx_dist, vx_t);
  parval1 = Point(vx_upar, vx_vpar);

  // For each pattern split, project the endpoints onto the current face
  // and select the endpoint mostly in line with the projected input points.
  // If both endpoints are cleary distant, both points are rejected
  vector<Point> endpts;
  int idx = -1;
  double min_dist = HUGE;
  for (size_t ki=0; ki<cand_split_.size(); ++ki)
    {
      // Fetch endpoints in this face
      Point pos1, pos2;
      double u1, v1, u2, v2, d1, d2;
      face_->closestPoint(cand_split_[ki].first.first, u1, v1, pos1, d1, epsge_);
      face_->closestPoint(cand_split_[ki].second.first, u2, v2, pos2, d2, epsge_);

      // Get closest boundary point
      double bd_upar1, bd_vpar1, bd_dist1, bd_t1;
      double bd_upar2, bd_vpar2, bd_dist2, bd_t2;
      Point bd_pos1, bd_pos2;
      ftEdgeBase *edg1, *edg2;
      if (only_outer_bd)
	{
	  edg1 = face_->closestOuterBoundaryPoint(pos1, dummy_vec, bd_upar1, bd_vpar1,
						  bd_pos1, bd_dist1, bd_t1);
	  edg2 = face_->closestOuterBoundaryPoint(pos2, dummy_vec, bd_upar2, bd_vpar2,
						  bd_pos2, bd_dist2, bd_t2);
	}
      else
	{
	  edg1 = face_->closestBoundaryPoint(pos1, dummy_vec, bd_upar1, bd_vpar1,
					     bd_pos1, bd_dist1, bd_t1);
	  edg2 = face_->closestBoundaryPoint(pos2, dummy_vec, bd_upar2, bd_vpar2,
					     bd_pos2, bd_dist2, bd_t2);
	}

      // Check if the current split is reasonably close (same edge or 
      // neighbouring edge)
      if (!(edg1 == vx_edg || edg1->next() == vx_edg || edg1->prev() == vx_edg ||
	    edg2 == vx_edg || edg2->next() == vx_edg || edg2->prev() == vx_edg))
	continue;

      // Check angle of candidate split
      double d3 = vx_pos.dist(bd_pos1);
      double d4 = vx_pos.dist(bd_pos2);
      Point other = (d3 < d4) ? bd_pos2 : bd_pos1;
      // double ang1 = normal.angle(bd_pos1 - cand_pt);
      // double ang2 = normal.angle(bd_pos2 - cand_pt);
      double ang = normal.angle(other - cand_pt);
      ang = fabs(0.5*M_PI - ang);
      if (ang > level_ang)
	continue;

      // Check if this candidate matches better than the previously found ones
      if (d3 < min_dist || d4 < min_dist)
	{
	  idx = (int)ki;
	  min_dist = std::min(d3, d4);
	  if (d4 < d3)
	    {
	      parval2 = Point(bd_upar1, bd_vpar1);
	    }
	  else
	    {
	      parval2 = Point(bd_upar2, bd_vpar2);
	    }
	}
    }

  // Make a first check on whether the found pattern split is appropriate
  if (min_dist > level_dist || parval1.dist(parval2) < epsge_)
    return false;

  return true;
}

//==========================================================================
int
RegularizeFace::nmbSplitPattern(const Point& p1, const Point& p2)
//==========================================================================
{
  double level_ang = M_PI/6.0;
  double level_dist = p1.dist(p2);
  level_dist *= 0.1;

  // Find the outer boundary points closest to the given input points
  int ind1, ind2;
  double par1, par2, dist1, dist2;
  Point close1, close2;
  face_->getBoundaryLoop(0)->closestPoint(p1, ind1, par1, close1, dist1);
  face_->getBoundaryLoop(0)->closestPoint(p2, ind2, par2, close2, dist2);

  // For each pattern split, project the endpoints onto the current face
  // and select the endpoint mostly in line with the projected input points.
  // If both endpoints are cleary distant, both points are rejected
  vector<Point> endpts;
  for (size_t ki=0; ki<cand_split_.size(); ++ki)
    {
      // Fetch endpoints in this face
      Point pos1, pos2;
      double u1, v1, u2, v2, d1, d2;
      face_->closestPoint(cand_split_[ki].first.first, u1, v1, pos1, d1, epsge_);
      face_->closestPoint(cand_split_[ki].second.first, u2, v2, pos2, d2, epsge_);

      size_t kj;
      Point vec1 = pos1 - close1;
      Point vec2 = close2 - pos1;
      double ang1 = vec1.angle(vec2);
      double scp1 = vec1*vec2;
      Point vec3 = pos2 - close1;
      Point vec4 = close2 - pos2;
      double ang2 = vec3.angle(vec4);
      double scp2 = vec3*vec4;
      if (fabs(d1-d2) < tol2_ && fabs(ang1-ang2) < angtol_)
	continue;  // Both endpoints seem to correspond to the outer boundary
      else if (ang1 < ang2 && scp1 > 0.0 && ang1 < level_ang &&
	  pos1.dist(close1) > level_dist && pos1.dist(close2) > level_dist)
	{
	  // Check if the point is found already
	  for (kj=0; kj<endpts.size(); ++kj)
	    //if (pos1.dist(endpts[kj]) < epsge_)
	    if (pos1.dist(endpts[kj]) < level_dist)
	      break;
	  if (kj == endpts.size())
	    endpts.push_back(pos1);
	  // for (kj=0; kj<endpts.size(); ++kj)
	  //   if (cand_split_[ki].first.first.dist(endpts[kj]) < level_dist)//epsge_)
	  //     break;
	  // if (kj == endpts.size())
	  //   endpts.push_back(cand_split_[ki].first.first);
	}
      else if (ang2 < ang1 && scp2 > 0.0 &&  ang2 < level_ang &&
	       pos2.dist(close1) > level_dist && pos2.dist(close2) > level_dist)
	{
	  // Check if the point is found already
	  for (kj=0; kj<endpts.size(); ++kj)
	    //if (pos2.dist(endpts[kj]) < epsge_)
	    if (pos2.dist(endpts[kj]) < level_dist)
	      break;
	  if (kj == endpts.size())
	    endpts.push_back(pos2);
	  // for (kj=0; kj<endpts.size(); ++kj)
	  //   if (cand_split_[ki].second.first.dist(endpts[kj]) < level_dist)//epsge_)
	  //     break;
	  // if (kj == endpts.size())
	  //   endpts.push_back(cand_split_[ki].second.first);
	}
    }
  return (int)endpts.size();
}

//==========================================================================
int
RegularizeFace::fetchSplitPattern(const Point& p1, const Point& p2,
				  vector<shared_ptr<Vertex> >& vx,
				  vector<int>& vx_idx, int nmb_idx,
				  vector<pair<Point,Point> >& pattern,
				  vector<pair<int,int> >& pattern_vx)
//==========================================================================
{
  if (cand_split_.size() == 0)
    return 0;  // No appropriate pattern information

#ifdef DEBUG_REG
  std::ofstream of("cross.g2");
#endif

  int nmb_split = 0;

  shared_ptr<Loop> loop = face_->getBoundaryLoop(0);
  Point vec1 = p2 - p1;
  double ang_fac = 0.1*M_PI;

  // Fetch parameter values to endpoints of chord
  Point pos1, pos2;
  double u1, v1, u2, v2, d1, d2;
  face_->closestPoint(p1, u1, v1, pos1, d1, epsge_);
  face_->closestPoint(p2, u2, v2, pos2, d2, epsge_);
  Point a(u1, v1);
  Point b(u2, v2);

  double mean_bd_dist = 0.0;
  double del = 0.05;
  for (size_t ki=0; ki<cand_split_.size(); ++ki)
    {
      // Fetch endpoints in this face
      face_->closestPoint(cand_split_[ki].first.first, u1, v1, pos1, d1, epsge_);
      face_->closestPoint(cand_split_[ki].second.first, u2, v2, pos2, d2, epsge_);

      // Check if this curve crosses the given chord
      Point c(u1, v1);
      Point d(u2, v2);
      Point w1 = b - a;
      Point w2 = d - c;
      Point w3 = a - c;
      double r1 = w3*w1;
      double r2 = w3*w2;
      double r3 = w1*w2;
      double r4 = w1*w1;
      double r5 = w2*w2;
      double det = r3*r3 - r4*r5;
      double s = (r1*r5 - r2*r3)/det;
      double t = (r1*r3 - r2*r4)/det;

#ifdef DEBUG_REG
      of << "410 1 0 4 255 0 0 255" << std::endl;
      of << "1" << std::endl;
      of << pos1 << " " << pos2 << std::endl;
#endif

      if (s <= del || s >= 1.0-del || t <= del || t >= 1.0-del)
	continue;   // Chords don't cross

      nmb_split++;

      // Check if the endpoints correspond to points on the outer boundary
      Point dummy_vec;
      double u3, v3, d3, t3, u4, v4, d4, t4;
      Point pos3, pos4;
      ftEdgeBase *edge1 = face_->closestBoundaryPoint(pos1, dummy_vec, 
						      u3, v3, pos3, d3, t3);
      if (!loop->isInLoop(edge1))
	continue;
      ftEdgeBase *edge2 = face_->closestBoundaryPoint(pos2, dummy_vec, 
						      u4, v4, pos4, d4, t4);
      if (!loop->isInLoop(edge2))
	continue;

      // Check if the pattern chord describes the current split curve well
      Point vec2 = pos4 - pos3;
      double ang = vec1.angle(vec2);
      if (fabs(0.5*M_PI-ang) >= ang_fac)
	continue;
      
      // Parameter pair corresponding to endpoints of pattern split curve
      pair<Point,Point> endpt = make_pair(Point(u3,v3), Point(u4,v4));
      pattern.push_back(endpt);
      pattern_vx.push_back(make_pair(-1, -1));

      mean_bd_dist += d3;
      mean_bd_dist += d4;
    }
  if (pattern.size() == 0)
    return nmb_split;

  mean_bd_dist /= (double)(2*pattern.size());
  double level_dist = 3.0*mean_bd_dist;

  // Search for correspondance between endpoints and candidate vertices
  // For each candidate vertex, find the closest endpoint
  vector<double> curr_dist(pattern.size(), 1.0e8);
  for (int kr=0; kr<nmb_idx; ++kr)
    {
      Point vx_pt = vx[vx_idx[kr]]->getVertexPoint();
      double mindist = 1.0e8;
      int min_ind = -1;
      for (size_t kj=0; kj<pattern.size(); ++kj)
	{
	  Point pt1 = face_->point(pattern[kj].first[0], pattern[kj].first[1]); 
	  Point pt2 = face_->point(pattern[kj].second[0], pattern[kj].second[1]);
	  double dist1 = vx_pt.dist(pt1);
	  double dist2 = vx_pt.dist(pt2);
	  if (dist1 < mindist)
	    {
	      mindist = dist1;
	      min_ind = 2*(int)kj;
	    }
	  if (dist2 < mindist)
	    {
	      mindist = dist2;
	      min_ind = 2*(int)kj + 1;
	    }
	}
      int id1 = min_ind/2;
      if (mindist < level_dist && mindist < curr_dist[id1])
	{
	  if (min_ind % 2 == 0)
	    pattern_vx[id1].first = vx_idx[kr];
	  else
	    pattern_vx[id1].second = vx_idx[kr];
	  curr_dist[id1] = mindist;
	}
    }
	  
  return nmb_split;
}

//==========================================================================
bool
RegularizeFace::splitWithPatternLoop()
//==========================================================================
{
  if (cand_split_.size() <= 1)
    return false;  // No appropriate pattern information

  // Fetch associated surface
  shared_ptr<ParamSurface> surf = face_->surface();

#ifdef DEBUG_REG
  std::ofstream of0("pattern0.g2");
#endif

  // Represent all split information as CurveOnSurface curves related to
  // the current face. Remember point correspondance.
  size_t ki;
  vector<shared_ptr<CurveOnSurface> > cand_cvs;
  shared_ptr<BoundedSurface> bd_sf;
  vector<pair<Point,Point> > corr_pts;
  Point dummy_vec;
  for (ki=0; ki<cand_split_.size(); ++ki)
    {
      // Consider only curves where both endpoints lie in the inner of the surface
      if (cand_split_[ki].first.second > 1 || cand_split_[ki].second.second > 1)
	continue;

      // Fetch endpoints in this face
      Point pos1, pos2;
      double u1, v1, u2, v2, d1, d2;
      face_->closestPoint(cand_split_[ki].first.first, u1, v1, pos1, d1, epsge_);
      face_->closestPoint(cand_split_[ki].second.first, u2, v2, pos2, d2, epsge_);

#ifdef DEBUG_REG
      of0 << "400 1 0 4 255 0 0 255" << std::endl;
      of0 << "2" << std::endl;
      of0 << cand_split_[ki].first.first << std::endl;
      of0 << pos1 << std::endl;
      of0 << "400 1 0 4 0 255 0 255" << std::endl;
      of0 << "2" << std::endl;
      of0 << cand_split_[ki].second.first << std::endl;
      of0 << pos2 << std::endl;
#endif

      // Avoid curves reaching boundaries in this face. 
      double bdu1, bdv1, bdd1, bdp1, bdu2, bdv2, bdd2, bdp2;
      Point bdpos1, bdpos2;
      /*ftEdgeBase *e1 =*/ 
      (void)face_->closestBoundaryPoint(pos1, dummy_vec, bdu1, bdv1,
					bdpos1, bdd1, bdp1);
      /*ftEdgeBase *e2 =*/ 
      (void)face_->closestBoundaryPoint(pos2, dummy_vec, bdu2, bdv2,
					bdpos2, bdd2, bdp2);
      if (bdd1 < tol2_ && bdd2 < tol2_)
	continue;
      
      if (pos1.dist(pos2) < tol2_)
      	continue;  // Not a valid curve segment

      corr_pts.push_back(make_pair(cand_split_[ki].first.first, pos1));
      corr_pts.push_back(make_pair(cand_split_[ki].second.first, pos2));

      // Represent split as CurveOnSurface curve
      vector<shared_ptr<CurveOnSurface> > tmp_cvs =
      	BoundedUtils::getTrimCrvsParam(surf, Point(u1,v1), Point(u2,v2),
      				       epsge_, bd_sf);
      cand_cvs.insert(cand_cvs.end(), tmp_cvs.begin(), tmp_cvs.end());
    }
#ifdef DEBUG_REG
  std::ofstream of1("pattern1.g2");
  for (size_t kv=0; kv<cand_cvs.size(); ++kv)
    {
      cand_cvs[kv]->spaceCurve()->writeStandardHeader(of1);
      cand_cvs[kv]->spaceCurve()->write(of1);
    }
#endif

//   // Remove the outer loop of candidate curves and those curves directly
//   // connected to them
//   removeOuterCands(cand_cvs);
// #ifdef DEBUG_REG
//   std::ofstream of2("pattern2.g2");
//   for (size_t kv=0; kv<cand_cvs.size(); ++kv)
//     {
//       cand_cvs[kv]->spaceCurve()->writeStandardHeader(of2);
//       cand_cvs[kv]->spaceCurve()->write(of2);
//     }
// #endif

  // Store remaining point correspondances
  size_t nmb_curr = corr_vx_pts_.size();
  for (ki=0; ki<corr_pts.size(); ++ki)
    {
      // Check if the registered point correspondance is still valid
      size_t kj, kr;
      for (kj=0; kj<cand_cvs.size(); ++kj)
	{
	  Point p1 = cand_cvs[kj]->ParamCurve::point(cand_cvs[kj]->startparam());
	  Point p2 = cand_cvs[kj]->ParamCurve::point(cand_cvs[kj]->endparam());
	  if (corr_pts[ki].second.dist(p1) < epsge_ ||
	      corr_pts[ki].second.dist(p2) < epsge_)
	    {
	      // Check if the correspondance is stored already
	      for (kr=nmb_curr; kr<corr_vx_pts_.size(); ++kr)
		{
		  if (corr_pts[ki].second.dist(corr_vx_pts_[kr].first) < epsge_ ||
		      corr_pts[ki].second.dist(corr_vx_pts_[kr].second) < epsge_)
		    break;
		}
	      if (kr == corr_vx_pts_.size())
		{
		  corr_vx_pts_.push_back(corr_pts[ki]);
		  break;
		}
	    }
	}
    }

  if (cand_cvs.size() > 0)
    {
      // Split the current surface according to the remaining pattern curves
      vector<shared_ptr<BoundedSurface> > sub_sfs =
	BoundedUtils::splitWithTrimSegments(bd_sf, cand_cvs, epsge_);

#ifdef DEBUG_REG
      std::ofstream of("split_surf.g2");
      for (size_t kr=0; kr<sub_sfs.size(); ++kr)
	{
	  sub_sfs[kr]->writeStandardHeader(of);
	  sub_sfs[kr]->write(of);
	}
#endif

      // Fetch info about vertices not belonging to the corners of the
      // initial face
      vector<shared_ptr<Vertex> > non_corner = 
	face_->getNonCornerVertices(bend_);
      removeInsignificantVertices(non_corner);
      
      // Create associated faces
      vector<shared_ptr<ftSurface> > faces = 
	RegularizeUtils::createFaces(sub_sfs, face_, epsge_, tol2_,
				     angtol_, non_corner);

      // Update topology and treat each sub face
      int nmb_faces = (int)faces.size();
      if (nmb_faces > 1)
	{
	  model_->removeFace(face_);
	  model_->append(faces, false, false);
	  for (int kj=0; kj<nmb_faces; )
	    {
	      RegularizeFace regularize(faces[kj], model_);
	      if (axis_.dimension() > 0)
		regularize.setAxis(centre_, axis_);
	      regularize.setDivideInT(divideInT_);
	      regularize.setCandSplit(cand_split_);
	      regularize.setSplitMode(split_mode_);
	      regularize.unsetTopLevel();
	      if (cand_split_.size() >  0)
		regularize.setCandSplit(cand_split_);

	      vector<shared_ptr<ftSurface> > faces2 = 
		regularize.getRegularFaces();

	      if (faces2.size() > 1)
		{
		  faces.erase(faces.begin()+kj);
		  nmb_faces--;
		  faces.insert(faces.end(), faces2.begin(), faces2.end());

		  // Check if any new faces may be joined across the seam
		  mergeSeams(faces, nmb_faces, faces2);
		}
	      else
		kj++;
	    }
	}
      sub_faces_.insert(sub_faces_.end(), faces.begin(), faces.end());
      return true;
    }
  return false;
}

//==========================================================================
void
RegularizeFace::removeOuterCands(vector<shared_ptr<CurveOnSurface> >& cand_cvs)
//==========================================================================
{
  if (cand_cvs.size() == 0)
    return; // Nothing to do

  // Get parameter box
  RectDomain dom = cand_cvs[0]->underlyingSurface()->containingDomain();
  Point corner1(dom.umin(), dom.vmin());
  Point corner2(dom.umax(), dom.vmax());

  // Find indices of the vectors to remove
  size_t ki, kj, kr;
  vector<int> idx;
  for (ki=0; ki<cand_cvs.size(); ++ki)
    {
      // Define a linear parameter domain curve passing through the 
      // midpoint of this curve being perpendicular to the tangent vector
      // of this parameter curve and large enough to cover the entire domain
      shared_ptr<ParamCurve> pcv = cand_cvs[ki]->parameterCurve();
      vector<Point> res = pcv->point(0.5*(pcv->startparam()+ pcv->endparam()), 1);
      Point dir(res[1][1], -res[1][0]);
      dir.normalize();
      double d1 = res[0].dist(corner1);
      double d2 = res[0].dist(corner2);
      Point p1 = res[0] - std::max(d1, d2)*dir;
      Point p2 = res[0] + std::max(d1, d2)*dir;

      shared_ptr<SplineCurve> lin_cv(new SplineCurve(p1, p2));
      
      // Compute the intersection between the current curve and this line
      vector<pair<double,double> > intersection_par;
      vector<int> pretopology;
      vector<pair<pair<double,double>, pair<double,double> > > int_crvs;
      intersect2Dcurves(lin_cv.get(), pcv.get(), epsge_, intersection_par,
			pretopology, int_crvs);

      if (intersection_par.size() == 0)
	continue;  // No intersection is found. Should not happen

      // Only one intersection is expected since the parameter curves are
      // all linear. Remember the parameter
      double par1 = intersection_par[0].first;

      // For all other curves, compute intersections with the line and remember
      // the smallest and larges intersection parameter
      double min_par = lin_cv->endparam();
      double max_par = lin_cv->startparam();
      for (kj=0; kj<cand_cvs.size(); ++kj)
	{
	  if (kj == ki)
	    continue;

	  // Compute the intersection between the current curve and the line
	  shared_ptr<ParamCurve> pcv2 = cand_cvs[kj]->parameterCurve();
	  vector<pair<double,double> > intersection_par2;
	  vector<int> pretopology2;
	  vector<pair<pair<double,double>, pair<double,double> > > int_crvs2;
	  intersect2Dcurves(lin_cv.get(), pcv2.get(), epsge_, intersection_par2,
			    pretopology2, int_crvs2);

	  if (intersection_par2.size() == 0)
	    continue;  // No intersection is found

	  // Only one intersection is expected since the parameter curves are
	  // all linear. Remember the parameter
	  double par2 = intersection_par2[0].first;
	  min_par = std::min(min_par, par2);
	  max_par = std::max(max_par, par2);
	}
	 
      // Check if the current curve is the first or the last intersection
      // with the line. In that case, it should be removed
      if (par1 <= min_par || par1 >= max_par)
	idx.push_back((int)ki);
    }

#ifdef DEBUG_REG
  std::ofstream of1("remove_cvs1.g2");
  for (size_t i1=0; i1<idx.size(); ++i1)
    {
      cand_cvs[idx[i1]]->spaceCurve()->writeStandardHeader(of1);
      cand_cvs[idx[i1]]->spaceCurve()->write(of1);
    }
#endif // DEBUG_REG

  // Remove also curves connected to the already indentified curves
  size_t nmb = idx.size();
  for (ki=0; ki<nmb; ++ki)
    {
      Point p1 = 
	cand_cvs[idx[ki]]->ParamCurve::point(cand_cvs[idx[ki]]->startparam());
      Point p2 = 
	cand_cvs[idx[ki]]->ParamCurve::point(cand_cvs[idx[ki]]->endparam());

      for (kj=0; kj<cand_cvs.size(); ++kj)
	{
	  // Check if the curve is already identified
	  for (kr=0; kr<idx.size(); ++kr)
	    if (idx[kr] == (int)kj)
	      break;
	  if (kr < idx.size())
	    continue;
	  
	  Point p3 = cand_cvs[kj]->ParamCurve::point(cand_cvs[kj]->startparam());
	  Point p4 = cand_cvs[kj]->ParamCurve::point(cand_cvs[kj]->endparam());
	  if (p1.dist(p3) < epsge_ || p1.dist(p4) < epsge_ ||
	      p2.dist(p3) < epsge_ || p2.dist(p4) < epsge_)
	    idx.push_back((int)kj);
	}
    }

#ifdef DEBUG_REG
  std::ofstream of2("remove_cvs2.g2");
  for (size_t i1=0; i1<idx.size(); ++i1)
    {
      cand_cvs[idx[i1]]->spaceCurve()->writeStandardHeader(of2);
      cand_cvs[idx[i1]]->spaceCurve()->write(of2);
    }
#endif // DEBUG_REG

  // Remove curves indicated by index vector. First sort indices
  std::sort(idx.begin(), idx.end());

  // Start from the last to remove to avoid messing up the indices
  for (int kh=(int)idx.size()-1; kh>=0; --kh)
    cand_cvs.erase(cand_cvs.begin()+idx[kh]);
}

//==========================================================================
void
RegularizeFace::classifyVertices()
//==========================================================================
{
  non_sign_vx_.clear();
  seam_vx_.clear();

  vector<shared_ptr<Vertex> > vx = face_->getNonCornerVertices(bend_);
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
	      if (edges[kj]->twin())
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
	  for (size_t kh=kj; kh<kr; ++kh)
	    {
	      faces.insert(edges[kh]->face()->asFtSurface());
	      if (edges[kh]->twin())
		faces.insert(edges[kh]->twin()->geomEdge()->face()->asFtSurface());
	    }
	  if (faces.size() != 2)
	    seam = false;

	  if (!seam)
	    break;
	}
      if (non_sign)
	non_sign_vx_.push_back(vx[ki]);
      else if (seam)
	seam_vx_.push_back(vx[ki]);
    }
}

//==========================================================================
void
RegularizeFace::removeInsignificantVertices(vector<shared_ptr<Vertex> >& vx, 
					    bool keep_T_joints, 
					    ftSurface *face)
//==========================================================================
{
  size_t kj, kr;
  Body *bd = model_->getBody();
  for (kj=0; kj<vx.size(); )
    {
      for (kr=0; kr<non_sign_vx_.size(); ++kr)
	if (vx[kj].get() == non_sign_vx_[kr].get())
	  break;
      if (kr < non_sign_vx_.size())
	{
	  vx.erase(vx.begin()+kj);
	  continue;
	}
      for (kr=0; kr<seam_vx_.size(); ++kr)
	{
	  if (vx[kj].get() == seam_vx_[kr].get())
	    break;
	  
	  // Check also distance
	  double td = vx[kj]->getDist(seam_vx_[kr]);
	  if (td < epsge_)
	    break;
	}
      if (kr < seam_vx_.size())
	{
	  vx.erase(vx.begin()+kj);
	  continue;
	}

      //vector<ftEdge*> edges = vx[kj]->allEdges();
      vector<ftEdge*> edges = vx[kj]->uniqueEdges(bd);
     if (edges.size() == 1 ||
	 (edges.size() == 2 /*&& !edges[0]->twin() && !edges[1]->twin()*/))
	{
	  // Closed loop, no neighbour. Check if the vertex is a corner
	  if (edges.size() == 2)
	    {
	      double t1 = edges[0]->parAtVertex(vx[kj].get());
	      double t2 = edges[1]->parAtVertex(vx[kj].get());
	      Point tan1 = edges[0]->tangent(t1);
	      Point tan2 = edges[1]->tangent(t2);
	      if (tan1.angle(tan2) < bend_)
		{
		  vx.erase(vx.begin()+kj);
		  continue;
		}
	    }
	  else
	    {
	      vx.erase(vx.begin()+kj);
	      continue;
	    }
	}

      if (keep_T_joints && vx[kj]->nmbUniqueEdges() > 2)
	{
	  kj++;
	  continue;  // Keep this vertex
	}

      // Check if the vertex belongs to an edge that is not continuing
      // past the first face
      shared_ptr<Vertex> vx2;
      pair<Point, Point> co_par1; 
      pair<Point, Point> co_par2;
      int dir1, dir2;
      double val1, val2;
      int not_extend = 
	RegularizeUtils::noExtension(vx[kj], 
				     (face!=NULL)?face:face_.get(), vx2,
				     co_par1, co_par2, dir1,
				     dir2, val1, val2,
				     bend_, false);
      if (not_extend > 0)
	{
	  vx.erase(vx.begin()+kj);
	  continue;
	}

      // Fetch the faces meeting in this vertex
      vector<ftSurface*> adjacent_faces = vx[kj]->faces();
      if (adjacent_faces.size() == 2)
	{
	  // Count number of twins
	  int no_twin = 0;
	  for (kr=0; kr<edges.size(); ++kr)
	    if (!edges[kr]->twin())
	      no_twin++;

	  // @@@ VSK. The check on number of twins may need some extra
	  // tuning. 
	  if (edges.size() == 3 && no_twin != 2)
	    {
	      // A seam. Remove vertex
	      vx.erase(vx.begin()+kj);
	      continue;
	    }
	  else
	    {
	      // If there is 4 edges and two of them has no twin,
	      // it might be a not notified seam
	      if (edges.size() == 4 && no_twin == 2)
	    {
	      // A seam. Remove vertex
	      vx.erase(vx.begin()+kj);
	      continue;
	    }
	    }
	}
	  

      kj++;
    }
}

//==========================================================================
void RegularizeFace::mergeSeams(vector<shared_ptr<ftSurface> >& faces, 
				int& nmb_faces,
				vector<shared_ptr<ftSurface> >& faces2)
//==========================================================================
{
  for (int kr=0; kr<(int)faces2.size(); ++kr)
    {
      // Identify surfaces that can be joined across a seam
      int pardir;
      ftSurface *other =
	identifySeamFaces(faces2[kr], faces, pardir);
      if (other)
	{
#ifdef DEBUG_REG
	  std::ofstream seam("seam_surfs.g2");
	  faces2[kr]->surface()->writeStandardHeader(seam);
	  faces2[kr]->surface()->write(seam);
	  other->surface()->writeStandardHeader(seam);
	  other->surface()->write(seam);
#endif

	  shared_ptr<ftSurface> joined = mergeSeamFaces(faces2[kr].get(), other, pardir);
	  if (joined.get())
	    {
	      size_t kh;
	      for (kh=0; kh<faces.size(); ++kh)
		if (faces[kh].get() == other)
		  break;
	      
	      if (kh < faces.size())
		{
		  model_->removeFace(faces[kh]);
		  faces.erase(faces.begin()+kh);
		  if ((int)kh < nmb_faces)
		    nmb_faces--;
		}

	      for (kh=0; kh<faces.size(); ++kh)
		if (faces[kh].get() == faces2[kr].get())
		  break;
	      
	      if (kh < faces.size())
		{
		  model_->removeFace(faces2[kr]);
		  faces.erase(faces.begin()+kh);
		  faces2.erase(faces2.begin()+kr);
		  kr--;
		}
	      
	      for (kh=0; kh<sub_faces_.size(); ++kh)
		if (sub_faces_[kh].get() == other)
		  break;
	      
	      if (kh < sub_faces_.size())
		sub_faces_.erase(sub_faces_.begin()+kh);

	      model_->append(joined, false, false);
	      faces.push_back(joined);
	    }
	}
    }
}



//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::chopOffRegBlocks(vector<shared_ptr<Vertex> >& concave_corners)
//==========================================================================
{
  vector<shared_ptr<ftSurface> > subfaces;

  // Collect significant vertices in the outer loop
  vector<shared_ptr<Vertex> > vxs = face_->getBoundaryLoop(0)->getSeqVertices();
  removeInsignificantVertices(vxs, true);
#ifdef DEBUG_REG
   std::ofstream of00("cand_vx.g2");
   of00 << "400 1 0 4 255 0 0 255" << std::endl;
   of00 << vxs.size()  << std::endl;
   for (size_t kh=0; kh<vxs.size(); kh++)
     of00 << vxs[kh]->getVertexPoint() << std::endl;
#endif
  if (vxs.size() < 6)
    return subfaces;  // Not a good candidate for this type of split

  if (concave_corners.size() == 0)
    return subfaces;  // No candidate split vertices

  size_t ki;
  vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > split_vxs;
  for (ki=0; ki<concave_corners.size(); ++ki)
    {
      // Look for candiate split from this corner. First localize the
      // candidate in the complete corner array
      int idc;
      shared_ptr<Vertex> vx = concave_corners[ki];
      for (idc=0; idc<(int)vxs.size(); ++idc)
	if (vxs[idc].get() == vx.get())
	  break;

      // Two possible blocks starts from this vertex. Check the
      // feasability of both
      vector<shared_ptr<Vertex> > cand1(4), cand2(4);
      int kj, kr;
      for (kj=0; kj<4; ++kj)
	{
	  kr = (idc + kj)%((int)vxs.size());
	  cand1[kj] = vxs[kr];
	}
      for (kj=0; kj<4; ++kj)
	{
	  kr = (idc - kj);
	  if (kr < 0)
	    kr += (int)vxs.size();
	  cand2[kj] = vxs[kr];
	}

#ifdef DEBUG_REG
      std::ofstream of0("cand_reg1_vx.g2");
      of0 << "400 1 0 4 0 255 0 255" << std::endl;
      of0 << "4 " << std::endl;
      for (kj=0; kj<4; kj++)
	of0 << cand1[kj]->getVertexPoint() << std::endl;
      
      std::ofstream of1("cand_reg2_vx.g2");
      of1 << "400 1 0 4 0 255 0 255" << std::endl;
      of1 << "4 " << std::endl;
      for (kj=0; kj<4; kj++)
	of1 << cand2[kj]->getVertexPoint() << std::endl;
#endif
      
      bool OKreg1 = false, OKreg2 = false;
      // Check if the end verticex has been involved in a split previously
      size_t kv, kh;
      vector<ftSurface*> faces = cand1[cand1.size()-1]->faces();
      for (kv=0; kv<faces.size(); ++kv)
	{
	  if (faces[kv] == face_.get())
	    continue;

	  int nmb = -1;
	  for (kh=0; kh<nonTjoint_faces_.size(); ++kh)
	    if (faces[kv] == nonTjoint_faces_[kh].get())
	      nmb++;
	  if (nmb < (int)nonTjoint_faces_.size())
	    break;
	}
      if (kv < faces.size())
	OKreg1 = false;
      else
	OKreg1 = RegularizeUtils::checkRegularity(cand1, face_);

      faces = cand2[cand2.size()-1]->faces();
      for (kv=0; kv<faces.size(); ++kv)
	{
	  if (faces[kv] == face_.get())
	    continue;

	  int nmb = -1;
	  for (kh=0; kh<nonTjoint_faces_.size(); ++kh)
	    if (faces[kv] == nonTjoint_faces_[kh].get())
	      nmb++;
	  if (nmb < (int)nonTjoint_faces_.size())
	    break;
	}
      if (kv < faces.size())
	OKreg2 = false;
      else
	OKreg2 = RegularizeUtils::checkRegularity(cand2, face_);
      pair<shared_ptr<Vertex>, shared_ptr<Vertex> > cand_split;
      bool found = false;
      if (OKreg1 && OKreg2)
	{
	  cand_split = make_pair(cand1[0], cand1[3]);
	  found = true;
	}
      else if (vxs.size() > 6 && OKreg1)
	{
	  cand_split = make_pair(cand1[0], cand1[3]);
	  found = true;
	}
      else if (vxs.size() > 6 && OKreg2)
	{
	  cand_split = make_pair(cand2[0], cand2[3]);
	  found = true;
	}
   
      if (found)
	{
	  // Make sure that the candidate is not found already
	  for (kr=0; kr<(int)split_vxs.size(); ++kr)
	    // if ((split_vxs[kr].first == cand_split.first &&
	    // 	 split_vxs[kr].second == cand_split.second) ||
	    // 	(split_vxs[kr].first == cand_split.second &&
	    // 	 split_vxs[kr].second == cand_split.first))
	    if (split_vxs[kr].first == cand_split.first ||
		split_vxs[kr].second == cand_split.second ||
		split_vxs[kr].first == cand_split.second ||
		split_vxs[kr].second == cand_split.first)
	      break;
	  if (kr == (int)split_vxs.size())
	    split_vxs.push_back(cand_split);
	}
    }

  if (split_vxs.size() == 0)
    return subfaces;  // No regular blocks are found

  // Perform splitting
  shared_ptr<ParamSurface> surf = face_->surface();
  vector<shared_ptr<CurveOnSurface> > trim_segments;
  shared_ptr<BoundedSurface> bd_sf;
  for (ki=0; ki<split_vxs.size(); ++ki)
    {
      // Find division curve between vertices
      Point parval1 = split_vxs[ki].first->getFacePar(face_.get());
      Point parval2 = split_vxs[ki].second->getFacePar(face_.get());
      vector<shared_ptr<CurveOnSurface> > curr_seg = 
	BoundedUtils::getTrimCrvsParam(surf, parval1, parval2, epsge_, bd_sf);

      // Check current segments
      vector<shared_ptr<Vertex> > next_vxs = 
	split_vxs[ki].first->getNextVertex(face_.get());
      vector<shared_ptr<Vertex> > next2 =
	split_vxs[ki].second->getNextVertex(face_.get());
      next_vxs.insert(next_vxs.end(), next2.begin(), next2.end());
      
      Point vx_point = split_vxs[ki].first->getVertexPoint();
      double curr_tol = 
	0.05*vx_point.dist(split_vxs[ki].second->getVertexPoint());

      // In this case, the result will not be a regular block if more than one
      // trimming curve is found. Dismiss the result in other cases
      if (curr_seg.size() > 1)
	curr_seg.clear();

      if (curr_seg.size() == 0 && centre_.dimension() > 0 &&
	  axis_.dimension() > 0)
	{
	  // No valid trim curve and an axis is defined. Check if a cylinder
	  // intersection is feasible
	  double d1 = split_vxs[ki].first->getVertexPoint().dist(centre_);
	  double d2 = split_vxs[ki].second->getVertexPoint().dist(centre_);
	  if (fabs(d1 - d2) < epsge_)
	    {
	      double cyl_rad = 0.5*(d1 + d2);
	      curr_seg = BoundedUtils::getCylinderIntersections(surf, centre_, 
								axis_, cyl_rad,
								epsge_, bd_sf);
	      if (curr_seg.size() > 1)
		curr_seg.clear();
	    }
	}

      if (curr_seg.size() > 0)
	RegularizeUtils::checkTrimSeg(curr_seg, next_vxs, 
				      split_vxs[ki].first->getVertexPoint(), 
				      split_vxs[ki].second->getVertexPoint(), 
				      curr_tol);

       if (curr_seg.size() > 0)
	trim_segments.insert(trim_segments.end(), curr_seg.begin(), curr_seg.end());
   }
      
  if (trim_segments.size() == 0)
    return subfaces;  // Segments not found or dismissed

 #ifdef DEBUG_REG
  std::ofstream out_file("split_segments.g2");
  for (size_t kj=0; kj<trim_segments.size(); ++kj)
    {
      shared_ptr<ParamCurve> cv = trim_segments[kj]->spaceCurve();
      cv->writeStandardHeader(out_file);
      cv->write(out_file);
    }
#endif

  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_sf, trim_segments, epsge_);

#ifdef DEBUG_REG
  std::ofstream of("split_surf.g2");
  for (size_t kr=0; kr<sub_sfs.size(); ++kr)
    {
      sub_sfs[kr]->writeStandardHeader(of);
      sub_sfs[kr]->write(of);
    }
#endif

  // Create faces
  vector<shared_ptr<Vertex> > non_corner = face_->getNonCornerVertices(bend_);
  subfaces = RegularizeUtils::createFaces(sub_sfs, face_, epsge_, 
					  tol2_, angtol_, non_corner);
  return subfaces;
}

//==========================================================================
vector<shared_ptr<ftSurface> >
RegularizeFace::connectToVertex(vector<shared_ptr<Vertex> >& concave_corners)
//==========================================================================
{
  vector<shared_ptr<ftSurface> > subfaces;

  // Collect significant vertices in the outer loop
  vector<shared_ptr<Vertex> > vxs = face_->getBoundaryLoop(0)->getVertices();
  removeInsignificantVertices(vxs);
  if (vxs.size() == 0)
    return subfaces;  // No vertices to connect to

  // Face properties
  shared_ptr<ParamSurface> surf = face_->surface();
  RectDomain dom = surf->containingDomain();

  size_t ki;
  vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > split_vxs;
  for (ki=0; ki<concave_corners.size(); ++ki)
    {
      // Look for candiate split from this corner. 
      shared_ptr<Vertex> curr_vx = concave_corners[ki];
      Point vx_point = curr_vx->getVertexPoint();

      // Fetch a vector in the given vertex pointing into the surface
      Point in_vec = RegularizeUtils::getInVec(curr_vx, face_);

      // Get the plane with which to divide the current face to get subdivision
      // information
      Point pnt;
      Point normal;
      if (centre_.dimension() > 0)
	{
	  pnt = centre_;
	  normal = (vx_point - centre_).cross(axis_);
	}
      else if (axis_.dimension() > 0)
	{
	  pnt = vx_point;
	  normal = axis_;
	}
      else
	RegularizeUtils::getDivisionPlane(face_, curr_vx, epsge_, pnt, normal);

      // Fetch boundary curve information
      size_t kr, kh;
      vector<ftEdge*> vx_edg = curr_vx->getFaceEdges(face_.get());
      vector<shared_ptr<ParamCurve> > vx_cvs;
      for (kr=0; kr<vx_edg.size(); ++kr)
	{
	  shared_ptr<ParamCurve> tmp = vx_edg[kr]->geomCurve();
	  for (kh=0; kh<vx_cvs.size(); ++kh)
	    if (vx_cvs[kh].get() == tmp.get())
	      break;
	  if (kh == vx_cvs.size())
	    vx_cvs.push_back(tmp);
	}

      // Compute distance to closest point on the boundary (except nearby edges)
      int close_idx;
      double close_dist;
      Point close_par;
      RegularizeUtils::getClosestBoundaryPar(face_, curr_vx, vx_cvs, vx_point, epsge_, 
					     close_idx, close_dist, close_par);
      Point close_pt = face_->point(close_par[0], close_par[1]);

      double cyl_rad = -1.0;
      int min_idx = RegularizeUtils::selectCandVx(face_, curr_vx, in_vec, 
						  vxs, dom, epsge_, angtol_,
						  centre_, normal, vx_cvs, 
						  close_dist, close_pt, cyl_rad);

      if (min_idx >= 0)
	{
	  pair<shared_ptr<Vertex>, shared_ptr<Vertex> > cand_split =
	    make_pair(curr_vx, vxs[min_idx]);

	  // Make sure that the candidate is not found already
	  for (kr=0; kr<split_vxs.size(); ++kr)
	    // if ((split_vxs[kr].first == cand_split.first &&
	    // 	 split_vxs[kr].second == cand_split.second) ||
	    // 	(split_vxs[kr].first == cand_split.second &&
	    // 	 split_vxs[kr].second == cand_split.first))
	    if (split_vxs[kr].first == cand_split.first ||
		split_vxs[kr].second == cand_split.second ||
		split_vxs[kr].first == cand_split.second ||
		split_vxs[kr].second == cand_split.first)
	      break;
	  if (kr == split_vxs.size())
	    split_vxs.push_back(cand_split);
	}
    }
   
  // Perform splitting
  vector<shared_ptr<CurveOnSurface> > trim_segments;
  shared_ptr<BoundedSurface> bd_sf;
  for (ki=0; ki<split_vxs.size(); ++ki)
    {
      // Find division curve between vertices
      Point parval1 = split_vxs[ki].first->getFacePar(face_.get());
      Point parval2 = split_vxs[ki].second->getFacePar(face_.get());
      vector<shared_ptr<CurveOnSurface> > curr_seg = 
	BoundedUtils::getTrimCrvsParam(surf, parval1, parval2, epsge_, bd_sf);

      // Check current segments
      vector<shared_ptr<Vertex> > next_vxs = 
	split_vxs[ki].first->getNextVertex(face_.get());
      vector<shared_ptr<Vertex> > next2 =
	split_vxs[ki].second->getNextVertex(face_.get());
      next_vxs.insert(next_vxs.end(), next2.begin(), next2.end());
      
      Point vx_point = split_vxs[ki].first->getVertexPoint();
      double curr_tol = 
	0.05*vx_point.dist(split_vxs[ki].second->getVertexPoint());
      RegularizeUtils::checkTrimSeg(curr_seg, next_vxs, 
				    split_vxs[ki].first->getVertexPoint(), 
				    split_vxs[ki].second->getVertexPoint(), 
				    curr_tol);

      // Check also angles between current edge loop and the new trimming curve
      // in the endpoints of this curve
      size_t kj;
      for (kj=0; kj<curr_seg.size(); )
	{
	  double min_ang1, min_ang2;
	  RegularizeUtils::angleInEndpoints(curr_seg[kj], split_vxs[ki].first,
					    split_vxs[ki].second, face_,
					     min_ang1, min_ang2);
	  if (min_ang1 < bend_ || min_ang2 < bend_)
	    {
	      // To small angle
	      curr_seg.erase(curr_seg.begin()+kj);
	    }
	  else
	    kj++;
	}

       if (curr_seg.size() > 0)
	trim_segments.insert(trim_segments.end(), curr_seg.begin(), curr_seg.end());
    }
      
  if (trim_segments.size() == 0)
    return subfaces;  // Segments not found or dismissed

 #ifdef DEBUG_REG
  std::ofstream out_file("split_segments.g2");
  for (size_t kj=0; kj<trim_segments.size(); ++kj)
    {
      shared_ptr<ParamCurve> cv = trim_segments[kj]->spaceCurve();
      cv->writeStandardHeader(out_file);
      cv->write(out_file);
    }
#endif

  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_sf, trim_segments, epsge_);

#ifdef DEBUG_REG
  std::ofstream of("split_surf.g2");
  for (size_t kr=0; kr<sub_sfs.size(); ++kr)
    {
      sub_sfs[kr]->writeStandardHeader(of);
      sub_sfs[kr]->write(of);
    }
#endif

  // Create faces
  vector<shared_ptr<Vertex> > non_corner = face_->getNonCornerVertices(bend_);
  subfaces = RegularizeUtils::createFaces(sub_sfs, face_, epsge_, 
					  tol2_, angtol_, non_corner);
  return subfaces;
}

 }
