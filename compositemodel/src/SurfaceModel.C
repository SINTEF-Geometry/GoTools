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

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/compositemodel/Path.h"
#include "GoTools/compositemodel/AdaptSurface.h"
#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/SISLCurveInterface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/compositemodel/cmUtils.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#include "GoTools/tesselator/TesselatorUtils.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/compositemodel/ftSurfaceSetPoint.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/ElementaryUtils.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/topology/FaceAdjacency.h"
#include "GoTools/topology/FaceConnectivityUtils.h"

//#define DEBUG
//#define DEBUG_REG

using std::vector;
using std::make_pair;
using std::max;
using std::min;
using std::ofstream;

namespace Go
{

//===========================================================================
  SurfaceModel::SurfaceModel(std::vector<shared_ptr<ftSurface> >& faces,
			     double space_epsilon,
			     double kink,  // Kink between adjacent surfaces 
			     bool adjacency_set) // Input faces
    //===========================================================================
    : CompositeModel(space_epsilon, 10.0*space_epsilon, kink, 10.0*kink),
      approxtol_(space_epsilon), tol2d_(1.0e-4)
  {
      if (faces.empty())
	  return;

    faces_.reserve(faces.size());
    for (size_t i = 0; i < faces.size(); ++i)
      faces_.push_back(faces[i]);
    initializeCelldiv();
    if (adjacency_set)
	setTopology();
    else
	buildTopology();
  }




  //===========================================================================
  SurfaceModel::SurfaceModel(double approxtol,
			     double gap,   // Gap between adjacent surfaces
			     double neighbour,  // Threshold for whether surfaces are adjacent
			     double kink,  // Kink between adjacent surfaces 
			     double bend, // Intended G1 discontinuity between adjacent surfaces
			     std::vector<shared_ptr<ftSurface> >& faces,
			     bool adjacency_set) // Input faces
    //===========================================================================
    : CompositeModel(gap, neighbour, kink, bend),
      approxtol_(approxtol), tol2d_(1.0e-4)
  {
      if (faces.empty())
	  return;

    faces_.reserve(faces.size());
    for (size_t i = 0; i < faces.size(); ++i)
      faces_.push_back(faces[i]);
    initializeCelldiv();
    if (adjacency_set)
	setTopology();
    else
	buildTopology();
  }




  //===========================================================================
  SurfaceModel::SurfaceModel(double approxtol,
			     double gap,   // Gap between adjacent surfaces
			     double neighbour,  // Threshold for whether surfaces are adjacent
			     double kink,  // Kink between adjacent surfaces 
			     double bend, // Intended G1 discontinuity between adjacent surfaces
			     std::vector<shared_ptr<ParamSurface> >& surfaces) // Input surfaces
    //===========================================================================
    : CompositeModel(gap, neighbour, kink, bend),
      approxtol_(approxtol), tol2d_(1.0e-4)
  {
      if (surfaces.empty())
	  return;

    faces_.reserve(surfaces.size());
    for (size_t i = 0; i < surfaces.size(); ++i)
      {
	// Check if the surface is likely to be closed. In that case split the surface
	vector<shared_ptr<ParamSurface> > surface_pieces = 
	  SurfaceModelUtils::checkClosedFaces(surfaces[i], neighbour);
	for (size_t j=0; j<surface_pieces.size(); ++j)
	  {
	    shared_ptr<ftSurface> newSurf;
	    newSurf.reset(new ftSurface(surface_pieces[j], (int)i));
	    faces_.push_back(newSurf);
	  }
      }
    initializeCelldiv();
    buildTopology();
  }




  //===========================================================================
  SurfaceModel::SurfaceModel(double approxtol,
			     double gap,   // Gap between adjacent surfaces
			     double neighbour,  // Threshold for whether surfaces are adjacent
			     double kink,  // Kink between adjacent surfaces 
			     double bend) // Intended G1 discontinuity between adjacent surfaces
    //===========================================================================
    : CompositeModel(gap, neighbour, kink, bend),
      approxtol_(approxtol), tol2d_(1.0e-4)
  {

  }





  //===========================================================================
  SurfaceModel::SurfaceModel(const SurfaceModel& sm)
    //===========================================================================
    : CompositeModel(sm),
      approxtol_(sm.approxtol_),
      tol2d_(sm.tol2d_),
      face_checked_(sm.face_checked_),
      limit_box_(sm.limit_box_)
  {
    // Rebuild faces based on ParamSurface. Edges between surfaces will be created
    // in buildTopology()

    closest_idx_ = sm.closest_idx_;
    faces_.reserve(sm.faces_.size());
    for (size_t i = 0; i < sm.faces_.size(); ++i)
      {
	shared_ptr<ParamSurface> newParamSurf;
	newParamSurf.reset(sm.faces_[i]->surface()->clone());
	shared_ptr<ftSurface> newSurf;
	newSurf.reset(new ftSurface(newParamSurf, sm.faces_[i]->getId()));
	faces_.push_back(newSurf);
      }

    initializeCelldiv();
    buildTopology();   // Sets boundary_curves_ and connectivity between surfaces
  }




  //===========================================================================
  SurfaceModel::~SurfaceModel()
  //===========================================================================
  {
  }


  //===========================================================================
  void SurfaceModel::setTolerances(double approxtol, double gap, double kink)
  //===========================================================================
  {
    setTolerances(approxtol, gap, 10.0*gap, kink, 10.0*kink);
  }


  //===========================================================================
  void SurfaceModel::setTolerances(double approxtol, double gap, double neighbour,
			      double kink, double bend)
  //===========================================================================
  {
    CompositeModel::setTolerances(gap, neighbour, kink, bend);
    approxtol_ = approxtol;
    //initializeCelldiv();
    buildTopology();
  }



  //===========================================================================
  void SurfaceModel::setBoundaryCurves()
  //===========================================================================
  {
    int i, j, k;
    FaceConnectivityUtils<ftEdgeBase,ftFaceBase> face_connectivity;

    vector< vector<shared_ptr<ftEdgeBase> > > loopvec;
    face_connectivity.BoundaryLoops(faces_, loopvec);
    vector<vector<ftFaceBase*> > grouped_faces;
    face_connectivity.disjointObjects(faces_, grouped_faces);

    // Separate the boundary curves, according to grouped_faces.
    vector<vector<shared_ptr<Loop> > > bd_loops(grouped_faces.size());
    // We save info about (appr) curve lengths.
    vector<vector<double> > lengths(grouped_faces.size());
    shared_ptr<ftEdgeBase> curr_edge;
    vector<shared_ptr<ftEdgeBase> > loop_edges;
    for (i = 0; i < (int)loopvec.size(); ++i) {
      curr_edge = loopvec[i][0];
      for (j = 0; j < (int)grouped_faces.size(); ++j) {
	for (k = 0; k < (int)grouped_faces[j].size(); ++k)
	  if (curr_edge->face() == grouped_faces[j][k])
	    break;
	if (k < (int)grouped_faces[j].size())
	  break;
      }
      loop_edges.clear();
      double length = 0;
      for (k = 0; k < (int)loopvec[i].size(); ++k) {
	length += cmUtils::estimatedCurveLength(loopvec[i][k].get()); //->geomEdge()->estimatedCurveLength();
	loop_edges.push_back(loopvec[i][k]);
      }
      lengths[j].push_back(length);
      bd_loops[j].push_back(shared_ptr<Loop>(new Loop(loop_edges, 
						      toptol_.neighbour)));
    }

    // We next push the longest curve up front, believing it is the outer curve.
    for (i = 0; i < (int)bd_loops.size(); ++i) {
      if (bd_loops[i].size() == 0)
	continue;

      std::vector<double>::const_iterator max_iter = std::max_element(lengths[i].begin(),
								      lengths[i].end());
      int index = (int)(max_iter - lengths[i].begin());
      std::swap(bd_loops[i][0], bd_loops[i][index]);
      std::swap(lengths[i][0], lengths[i][index]);
    }

    // We should be done.
    boundary_curves_ = bd_loops;
  }


  //===========================================================================
  ftMessage SurfaceModel::buildTopology(int first_idx, bool set_twin_face_info)
  //---------------------------------------------------------------------------
  //
  // Purpose: Find adjacency between faces and build a topology table 
  //          representing this adjacency.
  //
  //===========================================================================
  {
    ftMessage status;

    // Perform adjacency analysis
    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    adjacency.computeAdjacency(faces_, inconsistent_orientation_, first_idx);

    setBoundaryCurves();

    if (set_twin_face_info)
      {
	// Add information about faces at the boundary meeting only in vertices
	setVertexIdentity();
	
	// Add information about twin faces
	setTwinFaceInfo();
      }
    
    return status;
  }



  //===========================================================================
  ftMessage SurfaceModel::setTopology()
  //---------------------------------------------------------------------------
  //
  // Purpose: Fetch existing adjacency information between faces and build a 
  //          topology table representing this adjacency.
  //
  //===========================================================================
  {
    ftMessage status;

    // Compute connectivity information
    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    adjacency.setConnectivity(faces_);
    
    setBoundaryCurves();

    if (false)
      {
    // Add information about faces at the boundary meeting only in vertices
    setVertexIdentity();
   
    // Add information about twin faces
    setTwinFaceInfo();
      }

    return status;
  }



  //===========================================================================
  void SurfaceModel::setVertexIdentity()
  //---------------------------------------------------------------------------
  //
  // Purpose: Add information about faces at the boundary meeting only 
  //          in vertices
  //
  //===========================================================================
  {
    // Fetch all boundary vertices
    vector<shared_ptr<Vertex> > bd_vx;
    getBoundaryVertices(bd_vx);

    // Find identitical vertices
    vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > cand_vx;
    size_t ki, kj;
    for (ki=0; ki<bd_vx.size(); ++ki)
      for (kj=ki+1; kj<bd_vx.size(); ++kj)
	{
	  double dist = 
	    bd_vx[ki]->getVertexPoint().dist(bd_vx[kj]->getVertexPoint());
	  if (dist < toptol_.neighbour)
	    {
	      pair<shared_ptr<Vertex>, shared_ptr<Vertex> > ident = 
		make_pair(bd_vx[ki], bd_vx[kj]);
	      cand_vx.push_back(ident);
	    }
	}

    // Ensure one representation of identical vertices
    vector<shared_ptr<Vertex> > removed;
    for (ki=0; ki<cand_vx.size(); ++ki)
      {
	// Check if the first vertex in the pair is removed
	for (kj=0; kj<removed.size(); ++kj)
	  if (cand_vx[ki].first.get() == removed[kj].get())
	    break;

	if (kj < removed.size())
	  {
	    // Join vertices. Keep the second
	    ftEdge* edge = cand_vx[ki].second->getEdge(0);
	    edge->joinVertex(cand_vx[ki].first, cand_vx[ki].second);
	    removed.push_back(cand_vx[ki].first);
	  }
	else
	  {
	    // Join vertices. Keep the first
	    ftEdge* edge = cand_vx[ki].first->getEdge(0);
	    edge->joinVertex(cand_vx[ki].second, cand_vx[ki].first);
	    removed.push_back(cand_vx[ki].second);
	  }
       }
  }

  //===========================================================================
  void SurfaceModel::setTwinFaceInfo()
  //---------------------------------------------------------------------------
  //
  // Purpose: Add information about twin faces
  //
  //===========================================================================
  {
    for (size_t ki=0; ki<faces_.size(); ++ki)
      {
	ftSurface* curr = faces_[ki]->asFtSurface();
	if (!curr)
	  continue;

	if (curr->twin())
	  continue;  // Face already contains twin information
	
// 	if (!curr->allRadialEdges())
// 	  continue;  // Not surrounded by radial edges

	if (!curr->hasRadialEdges())
	  continue;  

	// Fetch candidate twin faces
	vector<ftSurface*> twin_cand = curr->fetchCorrespondingFaces();

	// Check for coincidence
	int idx = -1;
	int nmb_found = 0;
	Identity ident;
	shared_ptr<ParamSurface> srf1 = curr->surface();
	for (size_t kr=0; kr<twin_cand.size(); ++kr)
	  {
	    shared_ptr<ParamSurface> srf2 = twin_cand[kr]->surface();
#ifdef DEBUG_REG
	    std::ofstream of("twin_cand.g2");
	    srf1->writeStandardHeader(of);
	    srf1->write(of);
	    srf2->writeStandardHeader(of);
	    srf2->write(of);
#endif
	    
	    int res = ident.identicalSfs(srf1, srf2, toptol_.neighbour);
	    if (res == 1)
	      {
		// Coincidence
		idx = (int)kr;
		nmb_found++;
	      }
	  }

	if (nmb_found == 1)
	  {
	    // One coincident face found
	    curr->setTwin(twin_cand[idx]);
	  }
      }

  }

  //===========================================================================
  void SurfaceModel::limitVolume(double xmin, double xmax,
				 double ymin, double ymax,
				 double zmin, double zmax)
  //===========================================================================
  {
    BoundingBox b(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));
    limit_box_ = b;
  }


  //===========================================================================
  bool SurfaceModel::pointWithinLimits(const ftPoint& point)
  //===========================================================================
  {
    return pointWithinLimits(point.position());
  }


  //===========================================================================
  bool SurfaceModel::pointWithinLimits(const Point& point)
  //===========================================================================
  {
    bool res;
    try {
      res = limit_box_.containsPoint(point);
    }
    //catch (UnInitialized) {
    catch(std::exception) {
      MESSAGE("Limits not set. Call LimitVolume(...)");
      res = false;
    }
    return res;
  }


  //===========================================================================
  // Bounding box of the entire surface model
  BoundingBox SurfaceModel::boundingBox()
  //===========================================================================
  {
    return celldiv_ -> big_box();
  }


  //===========================================================================
  // Bounding box of one face
  BoundingBox SurfaceModel::boundingBox(int idx) const
  //===========================================================================
  {
    ASSERT(idx < (int)faces_.size());

    return faces_[idx]->boundingBox();
  }



  //===========================================================================
  bool SurfaceModel::isDegenerate(int idx) const
  //===========================================================================
  {
    ASSERT(idx < (int)faces_.size());
    ASSERT(idx >= 0);

    bool dummy[4];
    shared_ptr<ParamSurface> face = getSurface(idx);
    return face -> isDegenerate(dummy[0], dummy[1],
				dummy[2], dummy[3],
				toptol_.neighbour);
    
  }


  //===========================================================================
  void SurfaceModel::evaluate(int idx,      // Index of surface
			      double par[], // Parameter value
			      Point& pnt) const  // Result
  //===========================================================================
  {
    shared_ptr<ParamSurface> surf = getSurface(idx);
    if (surf.get() == 0) return;
    surf -> point(pnt, par[0], par[1]);
  }



  //===========================================================================
  void SurfaceModel::evaluate(int idx,      // Index
			      double par[], // Parameter value
			      int nder,     // Number of derivatives to compute, 0=only position
			      vector<Point>& der) const  // Result
  //===========================================================================
  {
    shared_ptr<ParamSurface> surf = getSurface(idx);
    if (surf.get() == 0) return;
    surf -> point(der, par[0], par[1], nder);
  }




  //===========================================================================
  void SurfaceModel::append(shared_ptr<ftSurface> face, bool set_twin,
			    bool adjacency_set, bool remove_twins, int idx)
  //===========================================================================
  {
#ifdef DEBUG
  bool isOK = checkShellTopology();
  if (!isOK)
    std::cout << "Shell, append (before). Topology inconsistencies" << std::endl;
#endif

    // Compute connectivity information related to the new face
    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    if (adjacency_set)
      {
	// Fetch all neighbours
	vector<ftSurface*> curr_faces;
	face->getAdjacentFaces(curr_faces);
	vector<ftFaceBase*> curr_faces2(curr_faces.begin(), curr_faces.end());
	
	// Add current
	curr_faces2.push_back(face.get());

	adjacency.setConnectivity(curr_faces2);
      }
    else
      {
	vector<pair<ftFaceBase*,ftFaceBase*> > orientation_inconsist;
	adjacency.computeFaceAdjacency(faces_, face, orientation_inconsist);
	if (orientation_inconsist.size() > 0)
	  inconsistent_orientation_.insert(inconsistent_orientation_.end(),
					   orientation_inconsist.begin(),
					   orientation_inconsist.end());
      }

    if (idx < 0)
      faces_.push_back(face);
    else
      faces_.insert(faces_.begin()+idx, face);

    if (boundary_curves_.size() > 0)
      boundary_curves_.erase(boundary_curves_.begin(), boundary_curves_.end());
    setBoundaryCurves();

    initializeCelldiv();

    // Add twin info for new face
    if (set_twin && !face->twin() /*&& face->allRadialEdges()*/)
      {
	// Fetch candidate twin faces
	vector<ftSurface*> twin_cand;
	if (face->allRadialEdges())
	  twin_cand = face->fetchCorrespondingFaces();
	else
	  face->getAdjacentFaces(twin_cand); 
	
	// Check for coincidence
	int idx = -1;
	int nmb_found = 0;
	Identity ident;
	shared_ptr<ParamSurface> srf1 = face->surface();
	for (size_t kr=0; kr<twin_cand.size(); ++kr)
	  {
	    shared_ptr<ParamSurface> srf2 = twin_cand[kr]->surface();
	    int res = ident.identicalSfs(srf1, srf2, toptol_.neighbour);
	    if (res == 1)
	      {
		// Coincidence
		idx = (int)kr;
		nmb_found++;
	      }
	  }

	if (nmb_found == 1)
	  {
	    if (remove_twins)
	      {
		removeFace(face);
		shared_ptr<ftSurface> other_face = 
		  fetchAsSharedPtr(twin_cand[idx]);

		if (other_face.get())
		  {
		    // Check if the twin face is a part of this model
		    double u1, v1;
		    Point pos1 = face->surface()->getInternalPoint(u1, v1);
		    double u2, v2, dist;
		    Point pos2;
		    other_face->closestPoint(pos1, u2, v2, pos2, dist, toptol_.gap);
		    Point norm = other_face->normal(u2, v2);
		    norm.normalize();

		    // Check if the model intersects a beam from close to the
		    // double face and inwards. If not, remove the face
		    norm *= -1;
		    Point pos = pos2 + 2*toptol_.neighbour*norm;
		    ftPoint result(0,0,0);
		    bool found = hit(pos, norm, result);
		    if (!found)
		      removeFace(other_face);
		  }
		
	      }
	    else
	      {
		// One coincident face found
		face->setTwin(twin_cand[idx]);
	      }
	  }
      }
#ifdef DEBUG
  isOK = checkShellTopology();
  if (!isOK)
    std::cout << "Shell, append (after). Topology inconsistencies" << std::endl;
#endif

  }



  //===========================================================================
  void SurfaceModel::append(std::vector<shared_ptr<ftSurface> > faces,
			    bool adjacency_set, bool set_twin)
  //===========================================================================
  {
#ifdef DEBUG
  bool isOK = checkShellTopology();
  if (!isOK)
    std::cout << "Shell, append (before). Topology inconsistencies" << std::endl;
#endif

  int nmb_faces = (int)faces_.size();
    for (size_t i = 0; i < faces.size(); ++i)
      faces_.push_back(faces[i]);
    initializeCelldiv();
    if (adjacency_set)
      setTopology();
    else
      buildTopology(nmb_faces, set_twin);

#ifdef DEBUG
  isOK = checkShellTopology();
  if (!isOK)
    std::cout << "Shell, append (after). Topology inconsistencies" << std::endl;
#endif
  }


  //===========================================================================
  void SurfaceModel::append(shared_ptr<SurfaceModel> anotherModel)
  //===========================================================================
  {
#ifdef DEBUG
  bool isOK = checkShellTopology();
  if (!isOK)
    std::cout << "Shell, append (before). Topology inconsistencies" << std::endl;
#endif

    for (size_t i = 0; i < anotherModel->faces_.size(); ++i)
      faces_.push_back(anotherModel->faces_[i]);
    initializeCelldiv();
    buildTopology();

#ifdef DEBUG
  isOK = checkShellTopology();
  if (!isOK)
    std::cout << "Shell, append (after). Topology inconsistencies" << std::endl;
#endif
  }

  //===========================================================================
  void SurfaceModel::closestBoundaryPoint(Point& pnt,     // Input point
					  Point& clo_pnt, // Found closest point
					  int& idx,           // Index of surface where the closest point is found
					  double clo_par[],   // Parameter value corresponding to the closest point
					  double& dist)       // Distance between input point and found closest point
  //===========================================================================
  {
    idx = -1;
    int loop_ix1 = -1, loop_ix2 = -1;
    size_t cv_ix;
    dist = std::numeric_limits<double>::max();
    double cv_par;
    for (size_t ki=0; ki<boundary_curves_.size(); ++ki)
      {
	for (size_t kj=0; kj<boundary_curves_[ki].size(); ++kj)
	  {
	    int cv_ind;
	    double tpar, tdist;
	    Point cv_close;
	    boundary_curves_[ki][kj]->closestPoint(pnt, cv_ind, tpar,
						   cv_close, tdist);
	    if (tdist < dist)
	      {
		cv_ix = cv_ind;
		loop_ix1 = (int)ki;
		loop_ix2 = (int)kj;
		cv_par = tpar;
		dist = tdist;
		clo_pnt = cv_close;
	      }
	  }
      }

    if (loop_ix1 >= 0)
      {
	shared_ptr<ftEdgeBase> edg0 =
	  boundary_curves_[loop_ix1][loop_ix2]->getEdge(cv_ix);
	ftEdge* edg = edg0->geomEdge();
	ftFaceBase *face0 = edg->face();
	ftSurface *face = face0->asFtSurface();
	idx = getIndex(face);
	Point fpar = edg->faceParameter(cv_par);
	clo_par[0] = fpar[0];
	clo_par[1] = fpar[1];
      }
  }



  //===========================================================================
  void SurfaceModel::closestPoint(Point& pnt,     // Input point
				  Point& clo_pnt, // Found closest point
				  int& idx,           // Index of surface where the closest point is found
				  double clo_par[],   // Parameter value corresponding to the closest point
				  double& dist)       // Distance between input point and found closest point
  //===========================================================================
  {
    ftPoint closest = closestPoint(pnt);
    clo_pnt = closest.position();
    clo_par[0] = closest.u();
    clo_par[1] = closest.v();
    idx = getIndex(closest.face());
    dist = clo_pnt.dist(pnt);
    closest_idx_ = idx;
  }


  //===========================================================================
  ftPoint SurfaceModel::closestPoint(const Point& point)
  //===========================================================================
  {
    // Before start, make sure that the face_checked vector is cleared
    fill(face_checked_.begin(), face_checked_.end(), false);

    // First, we locate the cell we're in @@@ now no longer necessary
    //int ix, iy, iz;
    //CellContaining(point, ix, iy, iz);
    //const ftCell& c = cells_[ix + iy*ncellsx_ + iz*ncellsx_*ncellsy_];

    // cout << "In cell " << ix + iy*ncellsx_ + iz*ncellsx_*ncellsy_ << endl;

    // Compute minimal distance from point to all cells
    double d;
    vector<ftCellInfo> cell_info(celldiv_ -> numCells());
    for (int i = 0; i < celldiv_ -> numCells(); ++i)
      {
	d = boxVecDist(celldiv_ -> getCell(i).box(), point);
	cell_info[i] = ftCellInfo(i, d);
      }
    // Sort the cell info vector
    sort(cell_info.begin(), cell_info.end());

    //for (int i=0; i<cells_.size(); ++i)
    //cout << "Cell " << cell_info[i].index_ 
    //     << " distance " << cell_info[i].dist_ << endl;
    //cout << "Closest cell : " << cell_info.front().dist_ << endl;
    //cout << "Farthest cell: " << cell_info.back().dist_  << endl;

    // Traverse the faces of that cell, computing the distance
    // to those faces' bounding boxes.
    int nmb_test = 0;
    Point cp;
    double dist;
    Point bestcp;
    double bestu = 0.0, bestv = 0.0;
    double bestdist = 1e100; // A gogool should be enough
    ftSurface* bestface = 0;
    int id;
    for (int j = 0; j< celldiv_ -> numCells(); ++j) {
      if (cell_info[j].dist_ > bestdist) {
	//      	    cout << cell_info[j].dist_ << " > " << bestdist << endl;
	break;
      }
      const ftCell& c = celldiv_ -> getCell(cell_info[j].index_);
      //      	cout << "Cell number: " << cell_info[j].index_ << endl;
      if (closest_idx_ >= 0)
      {
	  int i;
	  for (i=0; i<c.num_faces(); ++i) {
	    id = getIndex(c.face(i));
	      if (id == closest_idx_)
		  break;
	  }
	  if (i < c.num_faces() && !face_checked_[id])
	      if (!face_checked_[id]) {
		  ftPoint ret = closestPointLocal(ftPoint(point, c.face(i)));
		  nmb_test++;
		  if (ret.face() != 0) { // That is, a new point was found
		      cp = ret.position();
		      dist = point.dist(cp);
		      if (dist < bestdist) {
			  bestdist = dist;
			  bestcp = cp;
			  bestu = ret.u();
			  bestv = ret.v();
			  bestface = ret.face();
		      }
		  }
	      }
      }

      for (int i=0; i<c.num_faces(); ++i) {
	id = getIndex(c.face(i));
	if (!face_checked_[id]) {
	  ftPoint ret = closestPointLocal(ftPoint(point, c.face(i)));
	  nmb_test++;
	  if (ret.face() != 0) { // That is, a new point was found
	    cp = ret.position();
	    dist = point.dist(cp);
	    if (dist < bestdist) {
	      bestdist = dist;
	      bestcp = cp;
	      bestu = ret.u();
	      bestv = ret.v();
	      bestface = ret.face();
	    }
	  }
	}
      }
    }

    // Clear the face_checked vector
    fill(face_checked_.begin(), face_checked_.end(), false);

#ifdef DEBUG_SFMOD
    std::cout << "Number of faces checked: " << nmb_test << std::endl;
#endif
    return ftPoint(bestcp, bestface, bestu, bestv);
  }

  //===========================================================================
  void SurfaceModel::closestPoint(Point& point, int seed_ix, double seed[],
				  Point& clo_pt, int& idx, double clo_par[],
				  double& dist)
  //===========================================================================
  {
    fill(face_checked_.begin(), face_checked_.end(), false);
   ftPoint inpt;
    bool use_seed = (seed_ix >= 0);
    int ix = (seed_ix >= 0) ? seed_ix : 0;
    ftSurface *first = static_pointer_cast<ftSurface>(faces_[ix]).get();
    if (seed_ix >= 0)
      inpt = ftPoint(point, first, seed[0], seed[1]);
    else
      inpt = ftPoint(point, first);
    ftPoint close = closestPointLocal(inpt, use_seed);
    if (close.face() != 0)
      {
	clo_pt = close.position();
	idx = getIndex(close.face());
	clo_par[0] = close.u();
	clo_par[1] = close.v();
	dist = point.dist(clo_pt);
      }
    else
      closestPoint(point, clo_pt, idx, clo_par, dist);
  }
 

  //===========================================================================
  int SurfaceModel::nmbEntities() const
  //===========================================================================
  {
      return (int)faces_.size();
  }


  //===========================================================================
  shared_ptr<ftSurface> SurfaceModel::getFace(int idx) const
  //===========================================================================
  {
    shared_ptr<ftSurface> result;
    if (idx >= 0 && idx < (int)faces_.size())
      result = static_pointer_cast<ftSurface>(faces_[idx]);
    return result;
  }


  //===========================================================================
  vector<shared_ptr<ftSurface> > SurfaceModel::allFaces() const
  //===========================================================================
  {
    vector<shared_ptr<ftSurface> > result;
    for (size_t ki=0; ki<faces_.size(); ++ki)
      result.push_back(static_pointer_cast<ftSurface>(faces_[ki]));
    return result;
  }


  //===========================================================================
  shared_ptr<ParamSurface> SurfaceModel::getSurface(int idx) const
  //===========================================================================
  {
    shared_ptr<ftSurface> face = getFace(idx);
    if (face.get() == 0)
      {
	shared_ptr<ParamSurface> dummy;
	return dummy;
      }
    return face -> surface();
  }


  //===========================================================================
  shared_ptr<SplineSurface> SurfaceModel::getSplineSurface(int idx) const
  //===========================================================================
  {
    shared_ptr<ParamSurface> par_surf = getSurface(idx);
    shared_ptr<SplineSurface> spl_surf;
    if (par_surf.get() == 0)
      return spl_surf;

    spl_surf = dynamic_pointer_cast<SplineSurface, ParamSurface>(par_surf);
    if (spl_surf.get() != 0)
      return spl_surf;

    shared_ptr<ElementarySurface> elem_surf = dynamic_pointer_cast<ElementarySurface, ParamSurface>(par_surf);
    if (elem_surf.get() == 0)
      return spl_surf;

    return shared_ptr<SplineSurface>(elem_surf->geometrySurface());
  }

  //===========================================================================
  bool SurfaceModel::hasFace(ftSurface* face) const
  //===========================================================================
  {
    for (size_t i = 0; i < faces_.size(); ++i)
      if (faces_[i].get() == face)
	return true;

    return false;
  }

  //===========================================================================
  int SurfaceModel::getIndex(shared_ptr<ftSurface> face) const
  //===========================================================================
  {
    return getIndex(face.get());
  }


  //===========================================================================
  int SurfaceModel::getIndex(ftSurface* face) const
  //===========================================================================
  {
    for (size_t i = 0; i < faces_.size(); ++i)
      if (faces_[i].get() == face)
	  return (int)i;

    return -1;
  }


  //===========================================================================
  int SurfaceModel::getIndex(ParamSurface* surf) const
  //===========================================================================
  {
    for (size_t i = 0; i < faces_.size(); ++i)
      if (faces_[i]->surface().get() == surf)
	  return (int)i;

    return -1;
  }

//===========================================================================
shared_ptr<ftSurface> SurfaceModel::fetchAsSharedPtr(ftFaceBase *face) const
//===========================================================================
{
  shared_ptr<ftSurface> result;
  for (size_t i = 0; i < faces_.size(); ++i)
    if (faces_[i].get() == face)
      {
	result = static_pointer_cast<ftSurface>(faces_[i]);
	break;
      }
  return result;
}

//===========================================================================
void SurfaceModel::swapFaces(int idx1, int idx2)
//===========================================================================
{
  if (idx1 < 0 || idx1 >= (int)faces_.size() || idx2 < 0 ||
      idx2 >= (int)faces_.size())
    return;  // Do nothing

  // Swap
  std::swap(faces_[idx1], faces_[idx2]);
}

  //===========================================================================
  void SurfaceModel::turn(int idx)
  //===========================================================================
  {
    shared_ptr<ftFaceBase> curr = faces_[idx];
    shared_ptr<ParamSurface> srf = getSurface(idx);
    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    adjacency.releaseFaceAdjacency(curr);
    faces_.erase(faces_.begin()+idx);

    for (size_t ki=0; ki<inconsistent_orientation_.size(); )
      {
	if (inconsistent_orientation_[ki].first == curr.get() ||
	    inconsistent_orientation_[ki].second == curr.get())
	  inconsistent_orientation_.erase(inconsistent_orientation_.begin()+ki);
	else
	  ki++;
      }

    curr->clearInitialEdges();
    srf->swapParameterDirection();

    vector<pair<ftFaceBase*,ftFaceBase*> > orientation_inconsist;
    adjacency.computeFaceAdjacency(faces_, curr, orientation_inconsist);
    faces_.insert(faces_.begin()+idx, curr);
    if (orientation_inconsist.size() > 0)
      inconsistent_orientation_.insert(inconsistent_orientation_.end(),
				       orientation_inconsist.begin(),
				       orientation_inconsist.end());
  }

  //===========================================================================
  void SurfaceModel::turn()
  //===========================================================================
  {
    // First turn all surfaces, removing loop information
    for (size_t ki=0; ki<faces_.size(); ++ki)
      {
	faces_[ki]->clearInitialEdges();
	shared_ptr<ParamSurface> srf = getSurface((int)ki);
	srf->turnOrientation();
      }

    // Recompute topology information
    buildTopology();
  }

   //===========================================================================
  void SurfaceModel::initializeCelldiv()
  //===========================================================================
  {

      // Check if there are any faces. @jbt
      if (faces_.empty()) {
	  MESSAGE("No faces - return empty CellDivision object.");
	  celldiv_ = shared_ptr<CellDivision>();
	  return;
      }

      int nf = (int)faces_.size();
    vector<ftSurface*> surfaces;
    for (size_t i = 0; i < faces_.size(); ++i)
      {
	ftSurface* asSurf = faces_[i] -> asFtSurface();
	asSurf->setId((int)i);
	if (asSurf != 0) surfaces.push_back(asSurf);
      }

    face_checked_ = vector<bool>(nf, false);

    int min_cell = 3;
    int m = max(1, min(min_cell, nf/50));
    celldiv_ = shared_ptr<CellDivision> (new CellDivision(surfaces, m, m, m));
  }


  //===========================================================================
  const ftCell& SurfaceModel::getCell(int i) const
  //===========================================================================
  {
    return celldiv_ -> getCell(i);
  }

  //===========================================================================
  bool SurfaceModel::isClosed() const
  //===========================================================================
  {
      return nmbBoundaries() == 0;
  }


  //===========================================================================
  int SurfaceModel::nmbBoundaries() const
  //---------------------------------------------------------------------------
  //
  // Purpose: Return the number of boundaries (including holes).
  //
  //===========================================================================
  {
    int nmb_boundaries = 0;
    for (size_t i = 0; i < boundary_curves_.size(); ++i)
	nmb_boundaries += (int)boundary_curves_[i].size();

    return nmb_boundaries;
  }


  //===========================================================================
  ftCurve SurfaceModel::getBoundary(int idx)
  //---------------------------------------------------------------------------
  //
  // Purpose: Return the specified boundary (may be a hole).
  //
  //===========================================================================
  {
    ALWAYS_ERROR_IF(idx + 1 > nmbBoundaries(),
		    "There aren't that many boundaries.");

    int i = 0;
    int counter = 0;
    while ((int)boundary_curves_[i].size() - 1 + counter < idx) {
	counter += (int)boundary_curves_[i].size();
      ++i;
    }

    ftCurve edge_curve(CURVE_EDGE);

    // Get all boundary edges belonging to the current loop
    vector<shared_ptr<ftEdgeBase> > edges = 
      boundary_curves_[i][idx - counter]->getEdges();

    // Join edges into an ftCurve
    for (size_t kj=0; kj<edges.size(); ++kj)
      addSegment(edge_curve, edges[kj].get(), CURVE_EDGE);

    return edge_curve;
  }

//   //===========================================================================
//   ftCurve SurfaceModel::getBoundary(int idx)
//   //---------------------------------------------------------------------------
//   //
//   // Purpose: Return the specified boundary (may be a hole).
//   //
//   //===========================================================================
//   {
//     ALWAYS_ERROR_IF(idx + 1 > nmbBoundaries(),
// 		    "There aren't that many boundaries.");

//     int i = 0, j;
//     int counter = 0;
//     while ((int)boundary_curves_[i].size() - 1 + counter < idx) {
//       counter += boundary_curves_[i].size();
//       ++i;
//     }

//     ftEdgeBase* first_edge = boundary_curves_[i][idx - counter];

//     ftCurve edge_curve(CURVE_EDGE);
//     vector< vector<ftEdgeBase*> > loopvec;
//     top_table_.BoundaryLoops(loopvec);

//     for (i = 0; i < (int)loopvec.size(); ++i) {
//       for (j = 0; j < (int)loopvec[i].size(); ++j)
// 	if (loopvec[i][j] == first_edge)
// 	  break;
//       if (j < (int)loopvec[i].size())
// 	break;
//     }
//     int index = i;

//     for (j = 0; j < (int)loopvec[index].size(); ++j) {
//       ftEdgeBase* edge = loopvec[index][j];
//       // 	addSegment(edge_curve, loopvec[index][j], CURVE_EDGE);
//       addSegment(edge_curve, edge, CURVE_EDGE);
//     }

//     edge_curve.joinSegments(toptol_.gap, toptol_.neighbour, 
// 			    toptol_.kink, toptol_.bend);
//     return edge_curve;
//   }


  //===========================================================================
  ftCurve SurfaceModel::getGaps()
  //===========================================================================
  {
    ftCurve curve(CURVE_GAP);
    getCurveofType(CURVE_GAP, curve);
    return curve;
  }


  //===========================================================================
  void SurfaceModel::getGaps(vector<ftEdge*>& gaps)
  //===========================================================================
  {
    FaceConnectivityUtils<ftEdgeBase,ftFaceBase> connectivity;
    vector<ftEdgeBase*> vec;
    connectivity.cornersAndKinks(faces_, vec);

    for (size_t ki=0; ki<vec.size(); ki++)
    {
	ftEdgeBase* e0 = vec[ki];
	ftEdgeBase* e1 = e0->twin();
	ALWAYS_ERROR_IF(e1 == 0, "Unexpected edge type.");
	shared_ptr<FaceConnectivity<ftEdgeBase> > info = 
	  e0->getConnectivityInfo();
	for (size_t kj=0; kj<info->status_.size(); kj++)
	{
	    if (info->status_[kj] == 3)
	    {
		// A gap is found. Collect information
		ftEdge* curr_edge = e0->geomEdge();
		if (curr_edge)
		    gaps.push_back(curr_edge);

		break;
	    }
	}
    }	
  }

  //===========================================================================
  ftCurve SurfaceModel::getKinks()
  //===========================================================================
  {
    ftCurve curve(CURVE_KINK);
    getCurveofType(CURVE_KINK, curve);
    return curve;
  }


  //===========================================================================
  void SurfaceModel::getKinks(vector<ftEdge*>& kinks)
  //===========================================================================
  {
    FaceConnectivityUtils<ftEdgeBase,ftFaceBase> connectivity;
    vector<ftEdgeBase*> vec;
    connectivity.cornersAndKinks(faces_, vec);

    for (size_t ki=0; ki<vec.size(); ki++)
    {
	ftEdgeBase* e0 = vec[ki];
	ftEdgeBase* e1 = e0->twin();
	ALWAYS_ERROR_IF(e1 == 0, "Unexpected edge type.");
	shared_ptr<FaceConnectivity<ftEdgeBase> > info = 
	  e0->getConnectivityInfo();
	for (size_t kj=0; kj<info->status_.size(); kj++)
	{
	    if (info->status_[kj] == 1)
	    {
		// A gap is found. Collect information
		ftEdge* curr_edge = e0->geomEdge();
		if (curr_edge)
		    kinks.push_back(curr_edge);

		break;
	    }
	}
    }	
  }

  //===========================================================================
  ftCurve SurfaceModel::getG1Disconts()
  //===========================================================================
  {
    ftCurve curve(CURVE_CORNER);
    getCurveofType(CURVE_CORNER, curve);
    return curve;
  }

  //===========================================================================
  void SurfaceModel::getCorners(vector<ftEdge*>& corners)
  //===========================================================================
  {
    FaceConnectivityUtils<ftEdgeBase,ftFaceBase> connectivity;
    vector<ftEdgeBase*> vec;
    connectivity.cornersAndKinks(faces_, vec);

    for (size_t ki=0; ki<vec.size(); ki++)
    {
	ftEdgeBase* e0 = vec[ki];
	ftEdgeBase* e1 = e0->twin();
	ALWAYS_ERROR_IF(e1 == 0, "Unexpected edge type.");
	shared_ptr<FaceConnectivity<ftEdgeBase> > info = 
	  e0->getConnectivityInfo();
	for (size_t kj=0; kj<info->status_.size(); kj++)
	{
	    if (info->status_[kj] == 2 || info->status_[kj] == 3)
	    {
		// A corner or is found. Collect information
		ftEdge* curr_edge = e0->geomEdge();
		if (curr_edge)
		    corners.push_back(curr_edge);

		break;
	    }
	}
    }	
  }


  //===========================================================================
  ftSurface* SurfaceModel::getSurface2(int index) const
  //===========================================================================
  {
    ftFaceBase* aFace = faces_[index].get();
    if (aFace == 0) return 0;
    return aFace -> asFtSurface();
  }


  //===========================================================================
  vector<shared_ptr<SurfaceModel> > SurfaceModel::getConnectedModels() const
  //===========================================================================
  {
    vector<shared_ptr<SurfaceModel> > models;
    vector<shared_ptr<ftSurface> > curr_set;
    vector<shared_ptr<ftSurface> > all_sets;

    // Find start face for collecting a compact set
    for (size_t ki=0; ki<faces_.size(); ++ki)
      {
	curr_set.clear();

	// Check if the face is used already
	size_t kj;
	for (kj=0; kj<all_sets.size(); ++kj)
	  if (faces_[ki].get() == all_sets[kj].get())
	    break;

	if (kj < all_sets.size())
	  continue;  // The face already belongs to a set

	shared_ptr<ftSurface> face = 
	  dynamic_pointer_cast<ftSurface, ftFaceBase>(faces_[ki]);
	if (face.get() == 0)
	  continue;  // Unexpected surface type

	// Store first face and all connected faces
	getCurrConnectedModel(face, curr_set, all_sets);

	// Make surface model
	shared_ptr<SurfaceModel> curr_model = 
	  shared_ptr<SurfaceModel>(new SurfaceModel(approxtol_, 
						    toptol_.gap,
						    toptol_.neighbour,
						    toptol_.kink,
						    toptol_.bend,
						    curr_set, true));
	models.push_back(curr_model);
      }
    return models;
  }

  
  //===========================================================================
  void 
  SurfaceModel::getCurrConnectedModel(shared_ptr<ftSurface>& face,
				    vector<shared_ptr<ftSurface> >& curr_set,
				    vector<shared_ptr<ftSurface> >& all_sets) const
  //===========================================================================
  {
    // Store current face
    curr_set.push_back(face);
    all_sets.push_back(face);

    // Fetch all neighbours
    vector<ftSurface*> neighbours;
    face->getAdjacentFaces(neighbours);
    
    // For all neighbours, store the neighbour and all its neighbours
    // as long as they are not found already
    for (size_t ki=0; ki<neighbours.size(); ki++)
      {
	size_t kj;
	for (kj=0; kj<all_sets.size(); ++kj)
	  if (neighbours[ki] == all_sets[kj].get())
	    break;

	if (kj < all_sets.size())
	  continue;  // Face categorized already
	
	// Handle neighbours to neighbour
	shared_ptr<ftSurface> curr_face = fetchAsSharedPtr(neighbours[ki]);
	if (!curr_face.get())
	  continue;
	getCurrConnectedModel(curr_face, curr_set, all_sets);
      }
  }
  
  //===========================================================================
  ftPoint SurfaceModel::closestPointLocal(const ftPoint& point, bool use_seed) const
  //===========================================================================
  {
    double tol = toptol_.gap;  // Should be a parameter tolerance
    const Point& pt = point.position();
    int id = getIndex(point.face());
    ftSurface* curface = 0;
    Point cp;
    double u, v;
    double dist;
    bool finished = false;
    Point bestcp;
    double bestu, bestv;
    double bestdist = 1e100; // A gogool should be enough
    double closestpt_epsilon = toptol_.neighbour; // Maybe gap instead?
    ftSurface* bestface = 0;
    int nmb_checked = 0;
    double seedval[2];
    double *seed = use_seed ? seedval : 0;
    seedval[0] = point.u();
    seedval[1] = point.v();
    while (!finished) {
      if (!face_checked_[id]) {
	curface = dynamic_cast<ftSurface*>(faces_[id].get());
	face_checked_[id] = true;
	//  	    cout << "Face: " << id << endl;
	ASSERT(curface != 0);
	const Domain& dom = curface->surface()->parameterDomain();
	if (seed)
	  {
	    bool in_domain = false;
	    try {
	      Vector2D seed2(seed[0], seed[1]);
	      in_domain = dom.isInDomain(seed2, tol);
             }
	    catch (...)
	      {
		in_domain = false;
	      }
	    if (!in_domain)
	      seed = 0;
	  }
	curface->closestPoint(pt, u, v, cp, dist, closestpt_epsilon, NULL,
			      seed);
	nmb_checked++;
	if (dist < bestdist) {
	  bestdist = dist;
	  bestcp = cp;
	  bestu = u;
	  bestv = v;
	  bestface = curface;
	}
	// Check if the point was on the boundary
	bool on_boundary = dom.isOnBoundary(Vector2D(u, v),
					    toptol_.neighbour);
	if (on_boundary) {
	  //  		cout << "Face " << id << endl;
	  ftEdgeBase* boundary_edge
	    = curface->edgeClosestToPoint(u, v);
	  ftEdgeBase* twin = boundary_edge->twin();
	  if (!twin) // That is, there is no neighbour
	    {
	      bool corner = dom.isOnCorner(Vector2D(u, v),
					   toptol_.neighbour);
	      if (corner)
		{
		  shared_ptr<Vertex> vx1, vx2;
		  boundary_edge->geomEdge()->getVertices(vx1, vx2);
		  Point pt1 = vx1->getVertexPoint();
		  Point pt2 = vx2->getVertexPoint();
		  shared_ptr<Vertex> vx = (pt1.dist(pt) <= pt2.dist(pt)) ? vx1 : vx2;
		  vector<ftEdge*> edgs = vx->getEdges(curface);
		  for (size_t kr=0; kr<edgs.size(); ++kr)
		    {
		      if (edgs[kr] != boundary_edge)
			{
			  twin = edgs[kr]->twin();
			  if (twin)
			    break;
			}
		    }
		}
	      if (!twin)
		finished = true;
	    }
	  if (twin)
	    {
	    //  		    cout << "We're crossing a boundary!" << endl;
	    int id2 = getIndex(twin->face()->asFtSurface());
	    if (id2 == id)
	      {
		shared_ptr<ParamSurface> surf = getSurface(id);
		SurfaceTools::surface_seedfind(pt, *surf, 0,
					       seedval[0], seedval[1]);
	      }
	    else
	      {
		id = id2;
		if (use_seed)
		  {
		    ftEdge *twin2 = twin->geomEdge();
		    if (twin2)
		      {
			double tpar, tdist;
			Point close;
			twin2->closestPoint(bestcp, tpar, close, tdist);
			Point fpar = twin2->faceParameter(tpar);
			seedval[0] = fpar[0];
			seedval[1] = fpar[1];
			seed = seedval;
		      }
		  }
	      }
	    }
	} else // point was in the interior
	  finished = true;
      } else { // if face_checked_[id]
	//  	    cout << "That face was already checked" << endl;
	if (!bestface)
	  return ftPoint(point.position(), 0);
	else
	  finished = true;
      }
    }
#ifdef DEBUG_SFMOD
    std::cout << nmb_checked << "   ";
#endif

    return ftPoint(bestcp, curface, u, v);
  }



  //===========================================================================
  void SurfaceModel::addSegment(ftCurve& cv, ftEdgeBase* edge, ftCurveType ty)
  //===========================================================================
  {
    //   ftEdge* geomedge = edge->geomEdge();

    // The geometry of the edge curve is relevant, thus we
    // need the ftEdge. 

    // As curve may have been split, we extract relevant part.
    double tmin = edge->tMin();
    double tmax = edge->tMax();
    double knot_diff_tol = 1e-05;

    ftFaceBase* face_base_1 = edge->face();
    shared_ptr<ParamCurve> sub_curve_1 = shared_ptr<ParamCurve>
      (edge->geomEdge()->geomCurve()->subCurve(tmin, tmax, knot_diff_tol));
  
    ftFaceBase* face_base_2 = 0;
    shared_ptr<ParamCurve> sub_curve_2;
    if (edge->twin() != 0) {
      double tmin_twin = edge->twin()->tMin();
      double tmax_twin = edge->twin()->tMax();
      face_base_2 = edge->twin()->face();
      sub_curve_2 = shared_ptr<ParamCurve>
	(edge->twin()->geomEdge()->geomCurve()->subCurve(tmin_twin, 
							 tmax_twin,
							 knot_diff_tol));
    }

    shared_ptr<ParamCurve> param_curve_1;
    shared_ptr<ParamCurve> param_curve_2;
    shared_ptr<ParamCurve> common_space_curve;
  
    // looking for predefined spacecurve
    if (sub_curve_1->instanceType() == Class_SplineCurve) {
      common_space_curve = sub_curve_1;
    } else if (sub_curve_2.get() != 0 && 
	       sub_curve_2->instanceType() == Class_SplineCurve) {
      common_space_curve = sub_curve_2;
    }

    if (sub_curve_1->instanceType() == Class_CurveOnSurface) {
      shared_ptr<CurveOnSurface> curonsurf = 
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_curve_1);
      ALWAYS_ERROR_IF(curonsurf.get() == 0,
		      "Logically impossible error encountered in "
		      "ftFairingToolbox::addSegment()");

      param_curve_1 = curonsurf->parameterCurve();
      if (common_space_curve.get() == 0) 
	  common_space_curve = curonsurf->spaceCurve();
      if (common_space_curve.get() == 0) {
	shared_ptr<ParamCurve> param_curve = curonsurf->parameterCurve();
	shared_ptr<ParamSurface> base_surface = 
	  curonsurf->underlyingSurface();
	common_space_curve = 
	  shared_ptr<ParamCurve>(CurveCreators::liftParameterCurve(param_curve,
								   base_surface,
								   approxtol_));
      }
    }
    if (sub_curve_2.get() != 0 &&
	sub_curve_2->instanceType() == Class_CurveOnSurface) {
      shared_ptr<CurveOnSurface> curonsurf = 
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_curve_2);
      ALWAYS_ERROR_IF(curonsurf.get() == 0,
		      "Logically impossible error encountered in "
		      "ftFairingToolbox::addSegment()");
      param_curve_2 = curonsurf->parameterCurve();
      if (common_space_curve.get() == 0) 
	  common_space_curve = curonsurf->spaceCurve();
      if (common_space_curve.get() == 0) {
	shared_ptr<ParamCurve> param_curve = curonsurf->parameterCurve();
	shared_ptr<ParamSurface> base_surface = 
	  curonsurf->underlyingSurface();
	common_space_curve = 
	  shared_ptr<ParamCurve>(CurveCreators::liftParameterCurve(param_curve,
								   base_surface,
								   approxtol_));
      }
    }
  
    cv.appendSegment(ftCurveSegment(ty, JOINT_DISC,
				    face_base_1, 
				    face_base_2,
				    param_curve_1,
				    param_curve_2,
				    common_space_curve));
  }


  //===========================================================================
  void SurfaceModel::getCurveofType(ftCurveType type, ftCurve& curve)
  //===========================================================================
  {
    FaceConnectivityUtils<ftEdgeBase,ftFaceBase> connectivity;
    vector<ftEdgeBase*> vec;
    connectivity.cornersAndKinks(faces_, vec);

    vector<ftEdgeBase*> ft_vec;
    int i;
    for (i = 0; i < (int)vec.size(); ++i) {
      ftEdgeBase* edge = vec[i];
      ft_vec.push_back(edge);
    }

    int info_type = 0;
    if (type == CURVE_GAP)
      info_type = 3;
    else if (type == CURVE_KINK)
      info_type = 1;
    else if (type == CURVE_CORNER)
      info_type = 2;
    // The status is:
    // 0 : edges join smoothly. G1.
    // 1 : edges join, but normals are slightly discontinous. A kink.
    // 2 : edges join, but the normals are discontinous. G0.
    // 3 : edges almost join. A gap.
    // 4 : edges are totally discontinous.

    //     // cvtypes is indexed on status
    //     ftCurveType cvtypes[] = { CURVE_NOTYPE, CURVE_KINK, CURVE_CORNER,
    // 			      CURVE_GAP, CURVE_NOTYPE };

    // Now we have lots of ftEdgeBase pointers, we have to extract all the
    // 'bad' parts (status > 0) and create curve segments
    Point ps, pm, pe;
    for (i = 0; i < (int)ft_vec.size(); ++i)
      {
	ftEdgeBase* e0 = ft_vec[i];
	ftEdgeBase* e1 = e0->twin();
	//ALWAYS_ERROR_IF(e1 == 0, "Unexpected edge type.");
	if (e1 == 0)
	  continue; // Unexpected edge type
	shared_ptr<FaceConnectivity<ftEdgeBase> > info = 
	  e0->getConnectivityInfo();
	if (e0 == info->e2_ || e1 == info->e1_)
	  std::swap(e0, e1);
	  
	int num_subint = (int)info->status_.size();
	ALWAYS_ERROR_IF(num_subint<1, "Corrupt info!");

	for (int j=0; j<num_subint; ++j)
	  {
	    if (info->status_[j] != info_type)
	      continue;        // Not the correct type of incident

	    double par[2][2];
	    par[0][0] = info->parameters_[j].first;
	    par[0][1] = info->parameters_[j+1].first;
	    par[1][0] = info->parameters_[j].second;
	    par[1][1] = info->parameters_[j+1].second;
	    // If the parameter interval is tiny, skip to next segment.
	    if (fabs(par[0][1] - par[0][0]) < 1e-10) break; 

	    // Get curve piece of space curve
	    ftEdge* geomedge = e0->geomEdge();
	    ParamCurve* sc 
	      = geomedge->geomCurve()->subCurve(par[0][0],par[0][1]);
        
	    // If the curve's end and midpoint are separated by a
	    // distance smaller than the neighbourhood tolerance,
	    // delete it and skip to next segment.
	    // Is this alway right, or only if it is a boundary
	    // curve?
	    ps = sc->point(par[0][0]);
	    pm = sc->point(0.5*par[0][0] + 0.5*par[0][1]);
	    pe = sc->point(par[0][1]);
	    double dist =  ps.dist(pm) + pe.dist(pm);
	    if (dist < toptol_.neighbour) {
	      delete sc;
	      break;
	    }
        
	    // Get curve piece of parameter curves
	    SplineCurve* pc1 = 0;
	    //	pc1  = geomedge->parCurve()->subCurve(par[0][0],par[0][1]);
	    SplineCurve* pc2 = 0;
	    //	pc2  = e1->geomEdge()->parCurve()->subCurve(par[1][0],par[1][0]);
	    ftFaceBase* fb0 = e0->face();
	    ftFaceBase* fb1 = e1->face();
	    curve.appendSegment(ftCurveSegment(type, JOINT_DISC,
					       fb0, fb1,
					       shared_ptr<ParamCurve>(pc1), 
					       shared_ptr<ParamCurve>(pc2), 
					       shared_ptr<ParamCurve>(sc)));
	  }
      }
  }

  //===========================================================================
  bool SurfaceModel::removeFace(shared_ptr<ftSurface> face)
  //===========================================================================
  {
#ifdef DEBUG
  bool isOK = checkShellTopology();
  if (!isOK)
    std::cout << "Shell, remove face (before). Topology inconsistencies" << std::endl;
#endif

    int idx = getIndex(face);
    if (idx < 0 || idx >= (int)faces_.size())
      return false;

    if (face->twin())
      face->disconnectTwin();

    faces_.erase(faces_.begin()+idx);
    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    adjacency.releaseFaceAdjacency(face);

    if (boundary_curves_.size() > 0)
      boundary_curves_.erase(boundary_curves_.begin(), boundary_curves_.end());
    setBoundaryCurves();

    if (faces_.size() > 0)
      initializeCelldiv();

#ifdef DEBUG
    isOK = checkShellTopology();
    if (!isOK)
      std::cout << "Shell, remove face (after). Topology inconsistencies" << std::endl;
#endif

    return true;
  }

#if 0
  // Functionality not through quality assurance. There are still problems
  //===========================================================================
  bool SurfaceModel::replaceFace(shared_ptr<ftSurface> old_face,
				shared_ptr<ftSurface> face)
  //===========================================================================
  {
    int idx = getIndex(old_face);
    if (idx < 0 || idx >= (int)faces_.size())
      return false;

    ftSurface *twin = NULL;
    if (face->twin())
      {
	twin = face->twin();
	face->disconnectTwin();
      }
    
    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    bool replaced = adjacency.replaceFace(old_face, face);
    if (replaced)
      {
	faces_.erase(faces_.begin()+idx);
	faces_.push_back(face);
      }

    if (twin)
      face->setTwin(twin);

    if (replaced && faces_.size() > 0)
      initializeCelldiv();

    return replaced;
  }
#endif

  //===========================================================================
  void SurfaceModel::updateFaceTopology(shared_ptr<ftSurface> face)
  //===========================================================================
  {
    int idx = getIndex(face);
    if (idx < 0 || idx >= (int)faces_.size())
      return;

    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    adjacency.releaseFaceAdjacency(face);
    faces_.erase(faces_.begin()+idx);
    
    for (size_t ki=0; ki<inconsistent_orientation_.size(); )
      {
	if (inconsistent_orientation_[ki].first == face.get() ||
	    inconsistent_orientation_[ki].second == face.get())
	  inconsistent_orientation_.erase(inconsistent_orientation_.begin()+ki);
	else
	  ki++;
      }

    vector<pair<ftFaceBase*,ftFaceBase*> > orientation_inconsist;
    adjacency.computeFaceAdjacency(faces_, face, orientation_inconsist);
    faces_.insert(faces_.begin()+idx, face);
    if (orientation_inconsist.size() > 0)
      inconsistent_orientation_.insert(inconsistent_orientation_.end(),
				       orientation_inconsist.begin(),
				       orientation_inconsist.end());
  }

  //===========================================================================
  void SurfaceModel::tesselate(vector<shared_ptr<GeneralMesh> >& meshes) const
  //===========================================================================
  {
    int res[2];
    res[0] = res[1] = 20;
    tesselate(res, meshes);
  }

  //===========================================================================
  void SurfaceModel::tesselate(int uv_res, 
			       vector<shared_ptr<GeneralMesh> >& meshes) const
  //===========================================================================
  {
    tesselate(faces_, uv_res, meshes);
  }

  //===========================================================================
  void SurfaceModel::tesselate(const vector<shared_ptr<ftFaceBase> >& faces,
			       int uv_res, 
			       vector<shared_ptr<GeneralMesh> >& meshes) const
  //===========================================================================
  {
    meshes.clear();
    int u_res, v_res;
    for (size_t ki=0; ki<faces.size(); ki++)
    {
	// Make sure that boundary loops are oriented correctly
	bool fix;
	fix = faces[ki]->asFtSurface()->checkAndFixBoundaries();

	shared_ptr<ParamSurface> surf = faces[ki]->surface();

	TesselatorUtils::getResolution(surf.get(), u_res, v_res, uv_res);

	shared_ptr<GeneralMesh> mesh;
	try {
	  SurfaceModelUtils::tesselateOneSrf(surf, mesh, tol2d_, u_res, v_res);
	}
	catch (...)
	  {
	    // Don't get a mesh here
	    continue;
	  }
	meshes.push_back(mesh);
    }
  }

  //===========================================================================
  void SurfaceModel::tesselate(int resolution[],
			       vector<shared_ptr<GeneralMesh> >& meshes) const
  //===========================================================================
  {
    tesselate(faces_, resolution, meshes);
  }

  //===========================================================================
  void SurfaceModel::tesselate(const vector<shared_ptr<ftFaceBase> >& faces,
			       int resolution[],
			       vector<shared_ptr<GeneralMesh> >& meshes) const
  //===========================================================================
  {
    meshes.clear();
    for (size_t ki=0; ki<faces.size(); ki++)
    {
	// Make sure that boundary loops are oriented correctly
	bool fix;
	fix = faces[ki]->asFtSurface()->checkAndFixBoundaries();

	shared_ptr<ParamSurface> surf = faces[ki]->surface();

	shared_ptr<GeneralMesh> mesh;
	try {
	  SurfaceModelUtils::tesselateOneSrf(surf, mesh, tol2d_, 
					     resolution[0], resolution[1]);
	}
	catch (...)
	  {
	    // Don't get a mesh here
	    continue;
	  }
	meshes.push_back(mesh);
    }
  }

  //===========================================================================
  void SurfaceModel::tesselate(double density,
			       vector<shared_ptr<GeneralMesh> >& meshes) const
  //===========================================================================
  {
    tesselate(faces_, density, meshes);
  }

  //===========================================================================
  void SurfaceModel::tesselate(const vector<shared_ptr<ftFaceBase> >& faces,
			       double density,
			       vector<shared_ptr<GeneralMesh> >& meshes) const
  //===========================================================================
  {
    meshes.clear();

    int min_nmb = 3;
    int max_nmb = (int)(sqrt(1000000.0/(int)faces.size()));
    int u_res = 8; //20;
    int v_res = 8; //20;

    for (size_t ki=0; ki<faces.size(); ki++)
    {
	// Make sure that boundary loops are oriented correctly
	bool fix;
	fix = faces[ki]->asFtSurface()->checkAndFixBoundaries();

	shared_ptr<ParamSurface> surf = faces[ki]->surface();

	// Get resolution
	SurfaceModelUtils::setResolutionFromDensity(surf, density, min_nmb, 
						    max_nmb, tol2d_, 
						    u_res, v_res);

	shared_ptr<GeneralMesh> mesh;
	try {
	  SurfaceModelUtils::tesselateOneSrf(surf, mesh, tol2d_, u_res, v_res);
	}
	catch (...)
	  {
	    // Don't get a mesh here
	    continue;
	  }
	meshes.push_back(mesh);
    }
  }

  //===========================================================================
  shared_ptr<ftPointSet>  SurfaceModel::triangulate(double density) const
  //===========================================================================
  {
      int min_nmb = 3;
      int max_nmb = (int)(sqrt(1000000.0/(int)faces_.size()));
      int n = 8; //20;
      int m = 8; //20;
    shared_ptr<ftPointSet> triang = shared_ptr<ftPointSet>(new ftPointSet());

    vector<shared_ptr<ftSurface> > faces;
    vector<pair<int, int> > pnt_range;

    int nmb_pnt = 0;
    bool node_identity = false;
    for (size_t ki=0; ki<faces_.size(); ki++)
    {
	// Make sure that boundary loops are oriented correctly
	shared_ptr<ftSurface> curr_face = getFace((int)ki);
	bool fix;
	fix = curr_face->checkAndFixBoundaries();

	shared_ptr<ParamSurface> surf = curr_face->surface();

	// Get resolution
	SurfaceModelUtils::setResolutionFromDensity(surf, density, tol2d_,
						    min_nmb, max_nmb, n, m);

	// Make a mesh describing the surface
	shared_ptr<GeneralMesh> mesh;
	try {
	  SurfaceModelUtils::tesselateOneSrf(surf, mesh, tol2d_, n, m);
	}
	catch (...)
	  {
	    // Don't get a mesh here
	    continue;
	  }

	// Check mesh type
	node_identity = false;
	GenericTriMesh *trimesh = mesh->asGenericTriMesh();
	if (trimesh)
	    node_identity = true;
	
	// Make local triangulation
	shared_ptr<ftPointSet> local_triang = shared_ptr<ftPointSet>(new ftPointSet());
	meshToTriang(curr_face, mesh, n, m, local_triang, node_identity);

	// Add the surface triangulation to the current triangulation
	triang->append(local_triang);

#ifdef DEBUG
	// Debug output
	std::ofstream pointsout("pointsdump.g2");
	//std::ofstream edgesout0("data/triangedges.dat");
	std::ofstream edgessout("triangedges.g2");
	//triang->printXYZEdges(edgesout0);
	//triang->printXYZNodes(pointsout0, true);
	triang->write(edgessout);
	vector<Vector3D> bd_nodes;
	vector<Vector3D> inner_nodes;
	int k2;
	for (k2=0; k2<(int)triang->size(); ++k2)
	{
	    if ((*triang)[k2]->isOnBoundary())
		bd_nodes.push_back((*triang)[k2]->getPoint());
	    else
		inner_nodes.push_back((*triang)[k2]->getPoint());
	}
		
	pointsout << "400 1 0 4 255 0 0 255" << std::endl;
	pointsout << bd_nodes.size() << std::endl;
	for (k2=0; k2<(int)bd_nodes.size(); ++k2)
	    pointsout << bd_nodes[k2][0] << " " << bd_nodes[k2][1] << " " << bd_nodes[k2][2] << std::endl;
	pointsout << "400 1 0 4 0 255 0 255" << std::endl;
	pointsout << inner_nodes.size() << std::endl;
	for (k2=0; k2<(int)inner_nodes.size(); ++k2)
	    pointsout << inner_nodes[k2][0] << " " << inner_nodes[k2][1] << " " << inner_nodes[k2][2] << std::endl;
#endif

	// Handle common boundaries
	// First find pairs of faces meeting at a common boundary
	vector<ftSurface*> neighbours;
	curr_face->getAdjacentFaces(neighbours);
	for (size_t kj=0; kj<neighbours.size(); ++kj)
	{
	    // Check if this face is meshed already
	    size_t kr;
	    for (kr=0; kr<faces.size(); ++kr)
		if (faces[kr].get() == neighbours[kj])
		{
		    // A common boundary is found
		    triang->mergeBoundary(faces[kr], pnt_range[kr].first, 
					  pnt_range[kr].second, curr_face,
					  nmb_pnt, triang->size(), toptol_.gap);
		}
	}
	
	// Set range information
	faces.push_back(curr_face);
	pnt_range.push_back(make_pair(nmb_pnt, triang->size()));
	nmb_pnt = triang->size();

    }

#ifdef DEBUG_SFMOD
//     std::ofstream pointsout2("data/pointsdump2.g2");
     std::ofstream edgesout2("data/triangedges2.g2");
     std::ofstream edgesout3("data/triangedges3.dat");
//     triang->printXYZNodes(pointsout2, true);
    
     triang->write(edgesout2);
     triang->printXYZEdges(edgesout3);
#endif
	    

    return triang;
  }

  //===========================================================================
  void SurfaceModel::tesselatedCtrPolygon(vector<shared_ptr<LineCloud> >& ctr_pol) const
  //===========================================================================
  {
    tesselatedCtrPolygon(faces_, ctr_pol);
  }

  //===========================================================================
  void SurfaceModel::tesselatedCtrPolygon(const vector<shared_ptr<ftFaceBase> >& faces,
					  vector<shared_ptr<LineCloud> >& ctr_pol) const
  //===========================================================================
  {
    for (size_t ki=0; ki<faces.size(); ++ki)
      {
	shared_ptr<LineCloud> curr_pol = TesselatorUtils::getCtrPol(faces[ki]->surface().get());
	ctr_pol.push_back(curr_pol);
      }
  }


  //===========================================================================
  void SurfaceModel::meshToTriang(shared_ptr<ftSurface> face,
				  shared_ptr<GeneralMesh> mesh,
				  int n, int m, shared_ptr<ftPointSet> triang,
				  bool check_endpoint_identity) const
  //===========================================================================
  {
      int num_nodes = mesh->numVertices();

      // Insert all nodes in the triangulation
      int ki;
      double *nodes = mesh->vertexArray();
      double *param = mesh->paramArray();
      for (ki=0; ki<num_nodes; ++ki)
      {
	  Vector3D pos(nodes[3*ki], nodes[3*ki+1],nodes[3*ki+2]);
	  Vector2D par(param[2*ki], param[2*ki+1]);
	  int bnd = mesh->atBoundary(ki);
	  shared_ptr<ftFaceBase> facebase = face;
	  shared_ptr<ftSurfaceSetPoint> sfpnt = 
	      shared_ptr<ftSurfaceSetPoint>(new ftSurfaceSetPoint(pos, bnd, facebase, par));
	  triang->addEntry(sfpnt);
      }

      // Set neighbourhood information to complete the triangulation
      int num_triang = mesh->numTriangles();
      int num_triang_nodes = 3*num_triang;
      unsigned int *triang_idx = mesh->triangleIndexArray();
      for (ki=0; ki<num_triang_nodes; ki+=3)
      {
	  ftSamplePoint* pt1 = (*triang)[triang_idx[ki]];
	  ftSamplePoint* pt2 = (*triang)[triang_idx[ki+1]];
	  ftSamplePoint* pt3 = (*triang)[triang_idx[ki+2]];
	  pt1->addNeighbour(pt2);
	  pt2->addNeighbour(pt1);
	  pt2->addNeighbour(pt3);
	  pt3->addNeighbour(pt2);
	  pt3->addNeighbour(pt1);
	  pt1->addNeighbour(pt3);
      }
	  
      // In the case of a generic tri-mesh, some nodes may be unused.

      // Remove identitical boundary nodes
      if (check_endpoint_identity)
	  triang->cleanNodeIdentity(toptol_.gap);
  }

  //===========================================================================
  void 
  SurfaceModel::fetchSamplePoints(double density,
				  vector<SamplePointData>& sample_points) const
  //===========================================================================
  {
    sample_points.clear();
    int min_nmb = 3;
    int max_nmb = (int)(sqrt(1000000.0/(int)faces_.size()));

    // For each face, estimate the number of sample points and compute points
    for (size_t ki=0; ki<faces_.size(); ++ki)
      {
	ftSurface *curr = faces_[ki]->asFtSurface();
	if (!curr)
	  continue;  // Unexpected situation

	// Fetch number of sampling points
	int nmb_u, nmb_v;
	SurfaceModelUtils::setResolutionFromDensity(curr->surface(), density, 
						    tol2d_, min_nmb, max_nmb, 
						    nmb_u, nmb_v);

	// Sample face boundaries
	FaceUtilities::getBoundaryData(curr, 2*(nmb_u+nmb_v),
				       sample_points);

	// Sample the inner of the current face
	FaceUtilities::getInnerData(curr, nmb_u, nmb_v, sample_points);
      }
	
  }

  //===========================================================================
  void 
  SurfaceModel::getAllVertices(vector<shared_ptr<Vertex> >& vertices) const
  //===========================================================================
  {
    // Collect all vertices
    std::set<shared_ptr<Vertex> > all_vertices;  // All vertices in the model represented once

    for (size_t ki=0; ki<faces_.size(); ++ki)
      {
	vector<shared_ptr<Vertex> > curr_vertices = 
	  faces_[ki]->asFtSurface()->vertices();
	all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
      }

    vertices.clear();
    vertices.insert(vertices.end(), all_vertices.begin(), all_vertices.end());
  }

  //===========================================================================
  void 
  SurfaceModel::getBoundaryVertices(vector<shared_ptr<Vertex> >& vertices) const
  //===========================================================================
  {
    vertices.clear();
    for (size_t ki=0; ki<boundary_curves_.size(); ++ki)
      {
	for (size_t kj=0; kj<boundary_curves_[ki].size(); ++kj)
	  {
	    vector<shared_ptr<ftEdgeBase> > edges =  
	      boundary_curves_[ki][kj]->getEdges();
	    for (size_t kr=0; kr<edges.size(); ++kr)
	      {
		ftEdge* curr = edges[kr]->geomEdge();
		vertices.push_back(curr->getVertex(true));
	      }
	  }
      }
  }


///===========================================================================
vector<shared_ptr<ftEdge> > SurfaceModel::getBoundaryEdges() const
//===========================================================================
{
  // Fetch all faces lying at outer boundaries, i.e. all faces with no twin
  vector<shared_ptr<ftEdge> > bd_edges;
  size_t ki, kj, kr;
  for (ki=0; ki<boundary_curves_.size(); ++ki)
    for (kj=0; kj<boundary_curves_[ki].size(); ++kj)
      {
	shared_ptr<Loop> curr = boundary_curves_[ki][kj];
	size_t nmb = curr->size();
	for (kr=0; kr<nmb; ++kr)
	  {
	    shared_ptr<ftEdgeBase> edge = curr->getEdge(kr);
	    shared_ptr<ftEdge> edge2 = 
	      dynamic_pointer_cast<ftEdge, ftEdgeBase>(edge);
	    if (edge2.get())
	      bd_edges.push_back(edge2);  // The edge should always be an ftEdge
	  }
      }
  return bd_edges;
}

//===========================================================================
vector<shared_ptr<ftEdge> > SurfaceModel::getBoundaryEdges(int boundary_idx) const
//===========================================================================
{
  vector<shared_ptr<ftEdge> > bd_edges;

  if (boundary_idx < 0)
    return bd_edges;

  for (int i = 0; i < (int)boundary_curves_.size(); ++i)
    {
      if (boundary_idx < (int)boundary_curves_[i].size())
	{
	  shared_ptr<Loop> curr = boundary_curves_[i][boundary_idx];
	  int nmb = (int)curr->size();
	  for (int j = 0; j < nmb; ++j)
	    {
	      shared_ptr<ftEdgeBase> edge = curr->getEdge(j);
	      shared_ptr<ftEdge> edge2 = 
		dynamic_pointer_cast<ftEdge, ftEdgeBase>(edge);
	      if (edge2.get())
		bd_edges.push_back(edge2);  // The edge should always be an ftEdge
	    }
	  return bd_edges;
	}
      boundary_idx -= (int)boundary_curves_[i].size();
    }

  return bd_edges;
}

//===========================================================================
vector<shared_ptr<ftEdge> > SurfaceModel::getUniqueInnerEdges() const
//===========================================================================
{
  vector<shared_ptr<ftEdge> > edges;
  for (size_t ki=0; ki<faces_.size(); ++ki)
    {
      vector<shared_ptr<ftEdge> > curr_edges = 
	faces_[ki]->asFtSurface()->getAllEdges();
      
      // For each edge, check if it has a twin
      for (size_t kr=0; kr<curr_edges.size(); ++kr)
	{
	  if (curr_edges[kr]->twin())
	    {
	      // Check if this face or its twin is collected already
	      size_t kh;
	      for (kh=0; kh<edges.size(); ++kh)
		if (edges[kh].get() == curr_edges[kr].get() ||
		    edges[kh].get() == curr_edges[kr]->twin())
		  break;

		  if (kh == edges.size())
		    edges.push_back(curr_edges[kr]);
	    }
	}
    }
  return edges;
}

//===========================================================================
Body* SurfaceModel::getBody()
//===========================================================================
{
  for (size_t ki=0; ki<faces_.size(); ++ki)
    {
      Body *bd = faces_[ki]->asFtSurface()->getBody();
      if (bd)
	return bd;
    }
  return NULL;
}

   //===========================================================================
  bool SurfaceModel::simplifyTrimLoops(double& max_dist)
  //===========================================================================
  {
      int nmb_faces = (int)faces_.size();
    bool modified = false;
    max_dist = 0.0;
    int edge_pr_loop = 1; //4;

    vector<int> mod_faces;
    int ki;
    for (ki=0; ki<nmb_faces; ++ ki)
      {
	shared_ptr<ParamSurface> srf = getSurface(ki);
	shared_ptr<BoundedSurface> bd_srf = 
	  dynamic_pointer_cast<BoundedSurface, ParamSurface>(srf);
	if (!bd_srf.get())
	  continue;  // Not a trimmed surfaces, no reason to simplify loop

	ftSurface *face = faces_[ki]->asFtSurface();
	if (!face)
	  continue;

	int nmb_loop = face->nmbBoundaryLoops();
	int nmb_edge = face->nmbEdges();
	if (nmb_edge <= edge_pr_loop*nmb_loop)
	  continue;  // Loops not fragmented

	// Try to simplify loop
	double dist;
	bool mod = bd_srf->simplifyBdLoops(toptol_.gap, toptol_.kink, dist);

	if (mod)
	  {
	    mod_faces.push_back((int)ki);
	    modified = true;
	    max_dist = std::max(dist, max_dist);
	  }
      }

    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    for (ki=0; ki < (int)mod_faces.size(); ++ki)
      {
	// The boundary loops of this surface is changed. The edge loops
	// of the face must be updated accordingly
	// Remove topology information regarding this face
	shared_ptr<ftFaceBase> face = faces_[mod_faces[ki]];
	adjacency.releaseFaceAdjacency(face);

	// Remove outdated edges
	face->clearInitialEdges();
	faces_.erase(faces_.begin()+mod_faces[ki]);

	for (size_t kr=0; kr<inconsistent_orientation_.size(); )
	  {
	    if (inconsistent_orientation_[kr].first == face.get() ||
		inconsistent_orientation_[kr].second == face.get())
	      inconsistent_orientation_.erase(inconsistent_orientation_.begin()+kr);
	    else
	      kr++;
	  }

	// Compute new edges and update adjacency information
	vector<pair<ftFaceBase*,ftFaceBase*> > orientation_inconsist;
	adjacency.computeFaceAdjacency(faces_, face, orientation_inconsist);
	faces_.insert(faces_.begin()+mod_faces[ki], face);
	if (orientation_inconsist.size() > 0)
	  inconsistent_orientation_.insert(inconsistent_orientation_.end(),
					   orientation_inconsist.begin(),
					   orientation_inconsist.end());
      }

    if (modified)
      setBoundaryCurves();

    return modified;
  }

//===========================================================================
bool SurfaceModel::isAxisRotational(Point& centre, Point& axis, 
				    Point& vec, double& angle,
				    double& min_ang)
//===========================================================================
{
  // Only applicable for closed surface models
  if (nmbBoundaries() > 0)
    return false;

  // For each face, check if it is rotational. In that case extract the rotational
  // information and check if it is consistent. Non-rotational faces are collected for 
  // processing at a later stage
  centre.resize(0);  // Just to be sure
  angle = 0.0;
  vector<ftSurface*> not_rotational;
  vector<pair<Point,double> > slices;
  min_ang = 2.0*M_PI;
  for (size_t ki=0; ki<faces_.size(); ++ki)
    {
      Point curr_mid, curr_axis, curr_vec;
      double curr_ang;
      bool rotational = 
	faces_[ki]->asFtSurface()->surface()->isAxisRotational(curr_mid, curr_axis, 
							       curr_vec, curr_ang);
      min_ang = std::min(min_ang, curr_ang);
      if (rotational)
	{
	  if (centre.dimension() == 0)
	    {
	      centre = curr_mid;
	      axis = curr_axis;
	      vec = curr_vec;
	      angle = curr_ang;
	      slices.push_back(make_pair(curr_vec, curr_ang));
	    }
	  else
	    {
	      // Check consistence
	      double tmp_ang = axis.angle(curr_axis);
	      if (axis.angle(curr_axis) > toptol_.kink && 
		  fabs(M_PI-tmp_ang) > toptol_.kink)
		return false;
	      if (tmp_ang > toptol_.kink)
		{
		  // The axes are oppositely oriented. Adjust the start vector
		  Array<double,3> tmp_vec(curr_vec[0], curr_vec[1], curr_vec[2]);
		  MatrixXD<double, 3> mat;
		  mat.setToRotation(curr_ang, curr_axis[0], 
				    curr_axis[1], curr_axis[2]);  // Rotate the 
		  // start vector the angle curr_ang around curr_axis
		  Array<double,3> tmp_vec2 = mat*tmp_vec;
		  curr_vec = Point(tmp_vec2[0], tmp_vec2[1], tmp_vec2[2]);
		}

	      Point tmp_vec = curr_mid - centre;
	      tmp_ang = axis.angle(tmp_vec);
	      if (tmp_vec.length() > toptol_.gap && tmp_ang > toptol_.kink &&
		  fabs(M_PI-tmp_ang) > toptol_.kink)
		return false;

	      if (true /*vec.angle(curr_vec) > toptol_.kink || 
			 fabs(angle-curr_ang) > toptol_.gap*/)
		{
		  slices.push_back(make_pair(curr_vec, curr_ang));
		}
	    }
	}
      else
	not_rotational.push_back(faces_[ki]->asFtSurface());
    }
  
  // Check for a candidate axis;
  if (axis.dimension() == 0)
    return false;  // Not a rotational model

  // Join slices to find the complete sector of rotation.
  // Find candidate
  slices.erase(slices.begin());  // Already in current estimate
  while (true)
    {
      size_t ki;
      for (ki=0; ki<slices.size(); ++ki)
	{
	  if (fabs(vec.angle(slices[ki].first) - angle) < toptol_.kink &&
	      angle*slices[ki].second > 0.0 && 
	      angle+slices[ki].second < 2.0*M_PI+toptol_.kink)
	    {
	      angle += slices[ki].second;
	      slices.erase(slices.begin()+ki);
	      break;
	    }
	  else if (false /*More configurations later*/)
	    {
	      break;
	    }
	}
      if (angle > 2.0*M_PI - toptol_.kink)
	break;
      if (ki == slices.size())
	break;
    }

  // Check if all other slices fit into the same radial section
  while (true)
    {
      size_t ki;
      double angle2;
      // First remove slices that are equal to the complete sector
      for (ki=0; ki<slices.size();)
	{
	  if (fabs(angle - slices[ki].second) < toptol_.kink)
	    slices.erase(slices.begin()+ki);
	  else
	    ki++;
	}

      // Find a slice that fits with the start of the sector
      for (ki=0; ki<slices.size(); ++ki)
	{
	  if (vec.angle(slices[ki].first)< toptol_.kink)
	    {
	      angle2 = slices[ki].second;
	      break;
	    }
	}
      if (ki == slices.size())
	break;

      slices.erase(slices.begin()+ki);

      // Complete the current sector
      for (ki=0; ki<slices.size(); ++ki)
	{
	  if (fabs(vec.angle(slices[ki].first) - angle2) < toptol_.kink &&
	      angle2*slices[ki].second > 0.0 && 
	      angle2+slices[ki].second < angle+toptol_.kink)
	    {
	      angle2 += slices[ki].second;
	      slices.erase(slices.begin()+ki);
	      break;
	    }
	  else if (false /*More configurations later*/)
	    {
	      break;
	    }
	}
      if (angle2 > angle - toptol_.kink)
	angle2 = 0.0;
      else if (angle2 < angle - toptol_.kink)
	break;  // Not possible to complete the sector
    }

  if (slices.size() > 0)
    return false;   // All slices do not fit into the same sector

  // Check if the non-rotational faces agree with the identified
  // radial section
  Point norm1 = axis.cross(vec);
  norm1.normalize();
  Array<double,3> tmp_norm1(norm1[0], norm1[1], norm1[2]);
  MatrixXD<double, 3> mat;
  mat.setToRotation(angle, axis[0], axis[1], axis[2]);  // Rotate the normal
  // vector angle around axis
  Array<double,3> tmp_norm2 = mat*tmp_norm1;
  Point norm2(tmp_norm2[0], tmp_norm2[1], tmp_norm2[2]);

  // For each non-rotational face, check if it is planar and has a normal
  // vector parallel to one of the specified vectors
  for (size_t ki=0; ki<not_rotational.size(); ++ki)
    {
      Point face_norm;
      if (!not_rotational[ki]->surface()->isPlanar(face_norm, toptol_.kink))
	return false;

      if (face_norm.angle(norm1) > toptol_.kink && 
	  face_norm.angle(norm2) > toptol_.kink)
	return false;
    }
  
  return true;
}

//===========================================================================
bool SurfaceModel::isLinearSwept(Point& pnt, Point& axis, double& len)
//===========================================================================
{
  // Only applicable for closed surface models
  if (nmbBoundaries() > 0)
    return false;

  // To be valid, all faces in this model lie either in one out of two
  // planes or they are linear and the linear direction coincide and
  // coincide with both plane normals.
  // This is a restriction to the class of linearily swept models, but
  // will cover quite a number of cases.

  pnt.resize(0);
  axis.resize(0);
  len = 0.0;
  double eps = toptol_.gap;
  double angtol = toptol_.kink;
  Point ax;
  Point pt;
  double len2;
  double u1, v1;
  for (size_t ki=0; ki<faces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = faces_[ki]->asFtSurface()->surface();

      // Check if the surface is planar
      Point normal;
      bool planar = sf->isPlanar(normal, eps);

      // Check if the surface is linear
      Point dir, dir2;
      bool linear = sf->isLinear(dir, dir2, eps);

      if (!(planar || linear))
	return false;  // Not a linear sweep model

      // Check if the configuration allows for linear sweep
      if (axis.dimension() == 0)
	{
	  if (planar && linear)
	    {
	      if (ax.dimension() == 0)
		{
		  ax = normal;
		  pt = sf->getInternalPoint(u1, v1);
		}
	      else
		{
		  double ang = normal.angle(ax);
		  if (ang > angtol && fabs(M_PI-ang) > angtol && 
		      fabs(0.5*M_PI-ang) > angtol)
		    return false;
		  if (ang <= angtol || fabs(M_PI-ang) <= angtol)
		    {
		      Point pos = sf->getInternalPoint(u1, v1);
		      double dist = (pos - pnt)*axis;
		      if (fabs(dist) > eps)
			len2 = dist;
		    }
		}
	    }
	  else if (planar)
	    {
	      if (ax.dimension())
		{
		  double ang = normal.angle(ax);
		  if (ang > angtol && fabs(M_PI-ang) > angtol && 
		      fabs(0.5*M_PI-ang) > angtol)
		    return false;
		}
	      axis = normal;
	      pnt = sf->getInternalPoint(u1, v1);
	      axis.normalize();
	    }
	  else
	    {
	      axis = dir;
	      if (ax.dimension() && fabs(0.5*M_PI-dir.angle(ax)) < angtol)
		pnt = pt;
	      axis.normalize();
	    }
	}
      else
	{
	  if (planar && linear)
	    {
	      double ang = axis.angle(normal);
	      if (ang > angtol && fabs(M_PI-ang) > angtol && 
		  fabs(0.5*M_PI-ang) > angtol)
		return false;
	    }
	  else if (planar)
	    {
	      // This surface is situated at the start or end of
	      // the sweep. Check axis
	      double ang = axis.angle(normal);
	      if (ang > angtol && fabs(M_PI-ang) > angtol)
		return false;

	      // Check the position of the plane. First fetch a point
	      // in the plane
	      Point pos = sf->getInternalPoint(u1, v1);
	      if (pnt.dimension() == 0)
		pnt = pos;
	      else
		{
		  double dist = (pos - pnt)*axis;
		  if (fabs(dist) > eps)
		    {
		      if (fabs(len) <= eps)
			len = dist;
		      else if (fabs(len-dist) > eps)
			return false;
		    }
		}
	    }
	  else
	    {
	      // Check axis
	      double ang = axis.angle(dir);
	      if (ang > angtol && fabs(M_PI-ang) > angtol)
		return false;
	    }
	}
    } 
  if (axis.dimension() == 0)
    {
      if (ax.dimension() == 0)
	return false;
      else
	{
	  axis = ax;
	  pnt = pt;
	  len = len2;
	}
    }

  return true;
}

//===========================================================================
vector<shared_ptr<ftSurface> >  SurfaceModel::facesInPlane(Point& pnt, Point& axis)
//===========================================================================
{
  vector<shared_ptr<ftSurface> > faces;
  double eps = toptol_.gap;
  double angtol = toptol_.kink;
  double u1, v1;
  for (size_t ki=0; ki<faces_.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = faces_[ki]->asFtSurface()->surface();

      // Check if the surface is planar
      Point normal;
      (void)sf->isPlanar(normal, eps);
      double ang = axis.angle(normal);
      if (ang <= angtol || fabs(M_PI-ang) <= angtol)
	{
	  // The face normal coincides. Check the position of the plane
	  Point pos = sf->getInternalPoint(u1, v1);
	  double dist = (pos - pnt)*axis;
	  if (fabs(dist) <= eps)
	    faces.push_back(static_pointer_cast<ftSurface>(faces_[ki]));
	}
    }
  return faces;
}

//===========================================================================
bool SurfaceModel::allSplines() const
//===========================================================================
{
  for (size_t ki=0; ki<faces_.size(); ++ki)
    {
      if (!faces_[ki]->asFtSurface()->isSpline())
	return false;
    }

  return true;
}

//===========================================================================
bool SurfaceModel::isCornerToCorner() const
//===========================================================================
{
  size_t nmb_faces = faces_.size();
  size_t ki, kj;
  shared_ptr<ftEdge> edge1;
  shared_ptr<ftEdge> edge2;

  for (ki=0; ki<nmb_faces; ++ki)
    {
      for (kj=ki+1; kj<nmb_faces; ++kj)
	{
	  ftSurface* curr = faces_[kj]->asFtSurface();
	  if (!faces_[ki]->asFtSurface()->areNeighbours(curr, edge1, edge2))
	    continue;

	  if (!faces_[ki]->asFtSurface()->isCornerToCorner(curr, 
							   toptol_.gap))
	    return false;
	}
    }
  return true;
}

//===========================================================================
void SurfaceModel::makeCornerToCorner()
//===========================================================================
{
  bool changed = true;
  while (changed)
    {
      // As long as there is a new or unknown configuration, continue
      // to check for corner mismatches
      
      changed = false;  // No modifications performed yet

      size_t ki, kj;
      shared_ptr<ftEdge> edge1;
      shared_ptr<ftEdge> edge2;

      for (ki=0; ki<faces_.size(); ++ki)
	{
	  ftSurface *curr1 = faces_[ki]->asFtSurface();
	  for (kj=ki+1; kj<faces_.size(); ++kj)
	    {
	      ofstream of("corner_split.g2");

	      ftSurface *curr2 = faces_[kj]->asFtSurface();
	      if (!curr1->areNeighbours(curr2, edge1, edge2))
		continue;

	      if (!curr1->isCornerToCorner(curr2, toptol_.gap))
		{
		  changed = true;

		  curr1->surface()->writeStandardHeader(of);
		  curr1->surface()->write(of);
		  curr2->surface()->writeStandardHeader(of);
		  curr2->surface()->write(of);

		  vector<shared_ptr<ftSurface> > nfaces1;
		  vector<shared_ptr<ftSurface> > nfaces2;
		  curr1->splitAtInternalCorner(curr2, nfaces1, nfaces2,
					       toptol_.gap);
		  if (nfaces1.size() > 0)
		    {
		      for (size_t ix=0; ix<nfaces1.size(); ++ix)
			{
			  nfaces1[ix]->surface()->writeStandardHeader(of);
			  nfaces1[ix]->surface()->write(of);
			}

		      // The first face has been split
		      // Remove from model

		      
		      bool removed; 
		      removed = 
			removeFace(dynamic_pointer_cast<ftSurface, 
				   ftFaceBase>(faces_[ki]));
		    }

		  if (nfaces2.size() > 0)
		    {
		      // The first face has been split
		      for (size_t ix=0; ix<nfaces2.size(); ++ix)
			{
			  nfaces2[ix]->surface()->writeStandardHeader(of);
			  nfaces2[ix]->surface()->write(of);
			}

		      // Remove from model

		      
		      bool removed;
		      removed = 
			removeFace(dynamic_pointer_cast<ftSurface, 
				   ftFaceBase>(faces_[kj]));		    
		    }

		  if (nfaces1.size() > 0)
		    append(nfaces1);

		  if (nfaces2.size() > 0)
		    append(nfaces2);



		  break;
		}

	      if (changed)
		break;
	    }
	}
    }
}


//===========================================================================
void SurfaceModel::makeCommonSplineSpaces()
//===========================================================================
{
  bool changed = true;
  while (changed)
    {
      // As long as there is a new or unknown configuration, continue
      // to check for spline space mismatches
      
      changed = false;  // No modifications performed yet

      size_t ki, kj;
      shared_ptr<ftEdge> edge1;
      shared_ptr<ftEdge> edge2;

      for (ki=0; ki<faces_.size(); ++ki)
	{
	  ftSurface *curr1 = faces_[ki]->asFtSurface();
	  for (kj=ki+1; kj<faces_.size(); ++kj)
	    {
	      ftSurface *curr2 = faces_[kj]->asFtSurface();
	      if (!curr1->areNeighbours(curr2, edge1, edge2))
		continue;

	      if (!curr1->commonSplineSpace(curr2, toptol_.gap))
		{
		  changed = true;

#ifdef DEBUG
		  ofstream fp("common_tmp.g2");
		  shared_ptr<ParamSurface> sf1 = curr1->surface();
		  shared_ptr<ParamSurface> sf2 = curr2->surface();
		  sf1->writeStandardHeader(fp);
		  sf1->write(fp);
		  sf2->writeStandardHeader(fp);
		  sf2->write(fp);
#endif
		  curr1->makeCommonSplineSpace(curr2);
#ifdef DEBUG
		  ofstream fp2("common_tmp2.g2");
		  sf1 = curr1->surface();
		  sf2 = curr2->surface();
		  sf1->writeStandardHeader(fp2);
		  sf1->write(fp2);
		  sf2->writeStandardHeader(fp2);
		  sf2->write(fp2);
#endif
		  
		  break;
		}
	      if (changed)
		break;
	    }
	}
    }
}

//===========================================================================
void SurfaceModel::enforceCoLinearCoefs()
//===========================================================================
{
  double tol = toptol_.neighbour;
  double ang_tol = toptol_.bend;

  // Fetch all edges in the model (twins are represented once)
  vector<shared_ptr<ftEdge> > edges = getUniqueInnerEdges();
  size_t ki;
  for (ki=0; ki<edges.size(); ++ki)
    {
      if (!edges[ki]->twin())
	continue;  // No adjacent face with wich to enforce colinearity

      ftSurface *face1 = edges[ki]->face()->asFtSurface();
      ftEdge *twin = edges[ki]->twin()->geomEdge();
      ftSurface *face2 = twin->face()->asFtSurface();
      if (!(face1 && face2))
	continue;  // Not enough faces

      if ((!face1->isSpline()) || (!face2->isSpline()))
	continue;  // Associated surfaces are not spline surfaces, cannot
                   // modify coefficients

      // Enforce colinearity
      (void)FaceUtilities::enforceCoLinearity(face1, edges[ki].get(),
					      face2, tol, ang_tol);
    }

  // Ensure co linearity at vertices.
  // First fetch the vertices
  vector<shared_ptr<Vertex> > vxs;
  getAllVertices(vxs);
  for (ki=0; ki<vxs.size(); ++ki)
    {
      (void)FaceUtilities::enforceVxCoLinearity(vxs[ki], tol, ang_tol);
    }
}


//===========================================================================
void SurfaceModel::regularizeTwin(ftSurface *face, 
				  vector<shared_ptr<ftSurface> >& twinset)
//===========================================================================
{
  // Fetch surface to regularize 
  shared_ptr<ParamSurface> srf = face->surface();

#ifdef DEBUG_REG
  std::ofstream of("twin_reg.g2");
  srf->writeStandardHeader(of);
  srf->write(of);
#endif

  // Domain
  RectDomain dom = srf->containingDomain();

  // Get geometry tolerance
  double space_epsilon = toptol_.gap;

  // Check if the surface is trimmed
  shared_ptr<BoundedSurface> bd_srf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(srf);
  if (bd_srf.get())
    {
      srf = bd_srf->underlyingSurface();
      space_epsilon = bd_srf->outerBoundaryLoop().getSpaceEpsilon();

#ifdef DEBUG_REG
      std::ofstream bd("init_twin.g2");
      bd_srf->writeStandardHeader(bd);
      bd_srf->write(bd);
#endif
    }


  // Fetch debug data
  vector<shared_ptr<CurveOnSurface> > bd_crvs;
  if (bd_srf.get())
    {
      vector<CurveLoop> loops = bd_srf->allBoundaryLoops();
      vector<shared_ptr<CurveOnSurface> > bd_crvs;
      for (size_t k1=0; k1<loops.size(); ++k1)
	for (int k2=0; k2<loops[k1].size(); ++k2)
	  bd_crvs.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(loops[k1][k2]));
    }
    
  // For each twin face, fetch the boundary loop(s) and make corresponding
  // loops with respect to this surface
  double ptol = 10.0*toptol_.gap; // To limit the search
  size_t ki;
  vector<vector<shared_ptr<CurveOnSurface> > > twin_loops;
  for (ki=0; ki<twinset.size(); ++ki)
    {
      int nmb_loops = twinset[ki]->nmbBoundaryLoops();
      for (int kj=0; kj<nmb_loops; ++kj)
	{
	  vector<shared_ptr<CurveOnSurface> > loop_cvs;
	  shared_ptr<Loop> loop = twinset[ki]->getBoundaryLoop(kj);
	  
	  // Fetch all edges
	  size_t nmb_edges = loop->size();
	  size_t kr, km;
	  int missing = 0;
	  for (kr=0; kr<nmb_edges; kr=km)
	    {
	      shared_ptr<ftEdge> edge = 
		dynamic_pointer_cast<ftEdge,ftEdgeBase>(loop->getEdge(kr));
	      if (!edge.get())
		continue;  // Unexpected edge type

#ifdef DEBUG_REG
	      shared_ptr<CurveOnSurface> tmp_cv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(edge->geomCurve());
	      if (tmp_cv.get())
		{
		  tmp_cv->spaceCurve()->writeStandardHeader(of);
		  tmp_cv->spaceCurve()->write(of);
		}
#endif

	      // Check if a corresponding edge is assocated the face face
	      shared_ptr<ParamCurve> other_cv;
	      double t1, t2;
	      if (edge->hasEdgeMultiplicity())
		{
		  Point pos1 = edge->point(edge->tMin());
		  Point pos2 = edge->point(edge->tMax());
		  Point pos3 = edge->point(0.5*(edge->tMin()+edge->tMax()));
		  shared_ptr<Vertex> v1 =  edge->getVertex(true);
		  shared_ptr<Vertex> v2 =  edge->getVertex(false);
		  vector<ftEdge*> r_edges = 
		    edge->getEdgeMultiplicityInstance()->allEdges();
		  size_t kh;
		  for (kh=0; kh<r_edges.size(); ++kh)
		    if (r_edges[kh]->face() == face)
		      {
			other_cv = r_edges[kh]->geomCurve();
			double d1, d2, d3, t3;
			Point close1, close2, close3;
			shared_ptr<Vertex> v3 =  r_edges[kh]->getVertex(true);
			shared_ptr<Vertex> v4 =  r_edges[kh]->getVertex(false);
			double seed;
			if (v1.get() == v3.get())
			  seed = r_edges[kh]->parAtVertex(v3.get());
			else if (v1.get() == v4.get())
			  seed = r_edges[kh]->parAtVertex(v4.get());
			else
			  seed = 0.5*(r_edges[kh]->parAtVertex(v3.get()) +
				      r_edges[kh]->parAtVertex(v4.get()));
			other_cv->closestPoint(pos1, other_cv->startparam(), 
					       other_cv->endparam(), t1,  
					       close1, d1, &seed);

			if (v2.get() == v3.get())
			  seed = r_edges[kh]->parAtVertex(v3.get());
			else if (v2.get() == v4.get())
			  seed = r_edges[kh]->parAtVertex(v4.get());
			else
			  seed = 0.5*(r_edges[kh]->parAtVertex(v3.get()) +
				      r_edges[kh]->parAtVertex(v4.get()));
			other_cv->closestPoint(pos2, other_cv->startparam(), 
					       other_cv->endparam(), t2,  
					       close2, d2, &seed);
			seed = 0.5*(t1 + t2);
			other_cv->closestPoint(pos3, other_cv->startparam(), 
					       other_cv->endparam(), t3, close3, d3, &seed);
			if (d1 < toptol_.neighbour && d2 < toptol_.neighbour &&
			    d3 < toptol_.neighbour)
			  break;
		      }
		  if (kh == r_edges.size())
		    other_cv = shared_ptr<ParamCurve>();
		}
	      else 
		missing++;

	      for (km=kr+1; km<nmb_edges; ++km)
		{
		  if (!other_cv.get())
		    break;

		  // Check if the next edge is associated the same curve
		  shared_ptr<ftEdge> edge2 = 
		    dynamic_pointer_cast<ftEdge,ftEdgeBase>(loop->getEdge(km));
		  if (!edge2.get() || 
		      edge->geomCurve().get() != edge2->geomCurve().get())
		    break;

		  // Check if a corresponding edge is assocated the face face
		  double t3, t4;
		  if (!edge2->hasEdgeMultiplicity())
		    break;

		  Point pos1 = edge2->point(edge2->tMin());
		  Point pos2 = edge2->point(edge2->tMax());
		  Point pos3 = edge2->point(0.5*(edge2->tMin()+edge2->tMax()));
		  shared_ptr<Vertex> v1 =  edge2->getVertex(true);
		  shared_ptr<Vertex> v2 =  edge2->getVertex(false);
		  vector<ftEdge*> r_edges = 
		    edge2->getEdgeMultiplicityInstance()->allEdges();
		  size_t kh;
		  for (kh=0; kh<r_edges.size(); ++kh)
		    if (r_edges[kh]->face() == face &&
			r_edges[kh]->geomCurve().get() == other_cv.get())
		      {
			double d1, d2, d3, t5;
			Point close1, close2, close3;
			shared_ptr<Vertex> v3 =  r_edges[kh]->getVertex(true);
			shared_ptr<Vertex> v4 =  r_edges[kh]->getVertex(false);
			double seed;
			if (v1.get() == v3.get())
			  seed = r_edges[kh]->parAtVertex(v3.get());
			else if (v1.get() == v4.get())
			  seed = r_edges[kh]->parAtVertex(v4.get());
			else
			  seed = 0.5*(r_edges[kh]->parAtVertex(v3.get()) +
				      r_edges[kh]->parAtVertex(v4.get()));
			other_cv->closestPoint(pos1, other_cv->startparam(), 
					       other_cv->endparam(), t3,  
					       close1, d1, &seed);
			if (v2.get() == v3.get())
			  seed = r_edges[kh]->parAtVertex(v3.get());
			else if (v2.get() == v4.get())
			  seed = r_edges[kh]->parAtVertex(v4.get());
			else
			  seed = 0.5*(r_edges[kh]->parAtVertex(v3.get()) +
				      r_edges[kh]->parAtVertex(v4.get()));
			other_cv->closestPoint(pos2, other_cv->startparam(), 
					       other_cv->endparam(), t4,  
					       close2, d2, &seed);
			seed = 0.5*(t3 + t4);
			other_cv->closestPoint(pos3, other_cv->startparam(), 
					       other_cv->endparam(), t5, close3, d3, &seed);
			if (d1 < toptol_.neighbour && d2 < toptol_.neighbour &&
			    d3 < toptol_.neighbour)
			  break;
		      }
		  if (kh < r_edges.size())
		    t2 = t4;
		  else
		    break;
		}

	      shared_ptr<CurveOnSurface> sf_crv;
	      if (other_cv.get())
		{
		  // Store already exising boundary curve of face

		  // @@@ VSK 1012. Due to possible inconsistencies between
		  // ftEdge and EdgeVertex, the end parameters of other
		  // is not necessarily correct. Perform closest point to find
		  // end parameters
		  if (t1 > t2)
		    std::swap(t1, t2);
		  if (t2 - t1 > toptol_.gap)
		    {
		      // Avoid subdividing close to a knot
		      //double fuzzy_tol = std::min(1.0e-5, 1.e-3*(t2-t1));
		      //double fuzzy_tol = std::min(1.0e-4, 1.e-3*(t2-t1));
		      double fuzzy_tol = std::min(1.0e-7, 1.e-3*(t2-t1));
		      shared_ptr<SplineCurve> tmp_crv =
			shared_ptr<SplineCurve>(other_cv->geometryCurve());
		      if (tmp_crv.get())
			{
			  (void)tmp_crv->basis().knotIntervalFuzzy(t1, fuzzy_tol);
			  (void)tmp_crv->basis().knotIntervalFuzzy(t2, fuzzy_tol);
			}

		      shared_ptr<ParamCurve> crv = 
			shared_ptr<ParamCurve>(other_cv->subCurve(t1, t2));
		      sf_crv = 
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>(crv);
		      if (!sf_crv)
			sf_crv = 
			  shared_ptr<CurveOnSurface>(new CurveOnSurface(srf,
									crv,
									false));
		    }
		  else
		    other_cv = shared_ptr<ParamCurve>();
		}
	      if (!other_cv.get())
		{
		  // Make new trimming curve
		  shared_ptr<ParamCurve> crv = 
		    shared_ptr<ParamCurve>(edge->geomCurve()->subCurve(edge->tMin(), edge->tMax()));
		  shared_ptr<CurveOnSurface> sf_crv0 = 
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(crv);
		  if (sf_crv0)
		    sf_crv =
		      shared_ptr<CurveOnSurface>(new CurveOnSurface(srf,
								    sf_crv0->spaceCurve(),
								    false));
		  else
		    sf_crv = 
		      shared_ptr<CurveOnSurface>(new CurveOnSurface(srf,
								    crv,
								    false));
		}
	      bool success = sf_crv->ensureParCrvExistence(toptol_.gap, &dom);
	      if (!success)
		{
		    int stop_break;
		    stop_break = 1;
		}
	      loop_cvs.push_back(sf_crv);
	    }
#ifdef DEBUG_REG2
	  std::cout << "Missing radial edge (" << ki <<"): " << missing << std::endl;
#endif

	  // Ensure a consistent loop orientation
	  // A check for CCW orientation is performed later
	  Point pos0 = loop_cvs[0]->ParamCurve::point(loop_cvs[0]->startparam());
	  Point pos1 = loop_cvs[0]->ParamCurve::point(loop_cvs[0]->endparam());
	  for (kr=1; kr<loop_cvs.size(); ++kr)
	    {
	      Point pos2 = loop_cvs[kr]->ParamCurve::point(loop_cvs[kr]->startparam());
	      Point pos3 = loop_cvs[kr]->ParamCurve::point(loop_cvs[kr]->endparam());
	      if (kr == 1)
		{
		  double d0 = std::min(pos0.dist(pos2), pos0.dist(pos3));
		  double d1 = std::min(pos1.dist(pos2), pos1.dist(pos3));
		  if (d0 < d1)
		    {
		      // First curve has inconsistent orientation
		      loop_cvs[0]->reverseParameterDirection();
		      std::swap(pos0, pos1);
		    }
		}
	      
	      if (pos1.dist(pos2) < pos1.dist(pos3))
		pos1 = pos3;
	      else
		{
		  loop_cvs[kr]->reverseParameterDirection();
		  pos1 = pos2;
		}
	    }

	  // Check for the parameter curve of seam curves. Such curves may
	  // be inconsistent in the parameter domain
	  Point par0 = 
	    loop_cvs[loop_cvs.size()-1]->parameterCurve()->point(loop_cvs[loop_cvs.size()-1]->endparam());
	  for (kr=0; kr<loop_cvs.size(); ++kr)
	    {
	      Point par1 = loop_cvs[kr]->parameterCurve()->point(loop_cvs[kr]->startparam());
	      Point par2 = loop_cvs[kr]->parameterCurve()->point(loop_cvs[kr]->endparam());
	      if (par1.dist(par0) > ptol)
		{
		  // Try to regenerate parameter curve
		  int idx = (kr == loop_cvs.size()-1) ? 0 : (int)kr+1;
		  Point par3 = loop_cvs[idx]->parameterCurve()->point(loop_cvs[idx]->startparam());
		  bool success =
		    loop_cvs[kr]->makeParameterCurve(toptol_.gap, par0, par3);
		  if (success)
		    par2 = par3;
		}
	      par0 = par2;
	    }
	  twin_loops.push_back(loop_cvs);
	}
    }

#ifdef DEBUG_REG2
  for (ki=0; ki<twin_loops.size(); ++ki)
    for (size_t kr=0; kr<twin_loops[ki].size(); ++kr)
      {
	bool ok1 = twin_loops[ki][kr]->sameParameterDomain();
	bool ok2 = twin_loops[ki][kr]->sameOrientation();
	bool ok3 = twin_loops[ki][kr]->sameTrace(toptol_.neighbour);
	bool ok4 = twin_loops[ki][kr]->sameCurve(toptol_.neighbour);
	if (!(ok1 && ok2 && ok3 && ok4))
	  std::cout << " Loop curve inconsistence, crv nr: " << kr << std::endl;
      }
#endif

  for (ki=0; ki<twin_loops.size(); ++ki)
    if (!LoopUtils::paramIsCCW(twin_loops[ki], space_epsilon, 10.0*toptol_.gap))
      {
	int kr;
	int nmb = (int)twin_loops[ki].size();
	for (kr=0; kr<nmb; ++kr)
	  twin_loops[ki][kr]-> reverseParameterDirection();
	for (kr=0; kr<nmb/2; ++kr)
	  std::swap(twin_loops[ki][kr], twin_loops[ki][nmb-kr-1]);
      }

#ifdef DEBUG_REG
  std::ofstream of3("twin_loops.g2");
  std::ofstream of4("twin_par.g2");
  for (ki=0; ki<twin_loops.size(); ++ki)
    {
      for (size_t kr=0; kr<twin_loops[ki].size(); ++kr)
	{
	  twin_loops[ki][kr]->spaceCurve()->writeStandardHeader(of3);
	  twin_loops[ki][kr]->spaceCurve()->write(of3);
	  of3 << "400 1 0 4 255 0 0 255" << std::endl;
	  of3 << "1 " << std::endl;
	  of3 << twin_loops[ki][kr]->ParamCurve::point(twin_loops[ki][kr]->startparam()) << std::endl;
	  twin_loops[ki][kr]->parameterCurve()->writeStandardHeader(of4);
	  twin_loops[ki][kr]->parameterCurve()->write(of4);
	}
    }
#endif

  // Create trimmed surfaces
  vector<shared_ptr<BoundedSurface> > bd_srfs = 
    BoundedUtils::createTrimmedSurfs(twin_loops, srf, toptol_.gap);
  
#ifdef DEBUG_REG
  std::ofstream of2("twin_reg2.g2");
  for (ki=0; ki<bd_srfs.size(); ++ki)
    {
      bd_srfs[ki]->writeStandardHeader(of2);
      bd_srfs[ki]->write(of2);

      // Test
      bd_srfs[ki]->analyzeLoops();
      //BoundedUtils::fixInvalidBoundedSurface(bd_srfs[ki]);
      int stop_break;
      stop_break = 1;
    }
#endif

  vector<shared_ptr<ftSurface> > twin_faces(bd_srfs.size());
  for (ki=0; ki<bd_srfs.size(); ++ki)
    {
      shared_ptr<ftSurface> curr =
	shared_ptr<ftSurface>(new ftSurface(bd_srfs[ki], -1));
      curr->setBody(face->getBody());
      (void)curr->createInitialEdges(toptol_.gap, toptol_.kink);
      twin_faces[ki] = curr;
    }
  
  shared_ptr<ftSurface> face2 = fetchAsSharedPtr(face);
  (void)removeFace(face2);
  for (ki=0; ki<twin_faces.size(); ++ki)
    {
      append(twin_faces[ki], false, false);
    }

#ifdef DEBUG_REG2
  for (ki=0; ki<twin_faces.size(); ++ki)
    {
      vector<shared_ptr<ftEdge> > edges = twin_faces[ki]->getAllEdges();
      for (size_t kr=0; kr<edges.size(); ++kr)
	std::cout << kr << " " << edges[kr].get() << " " << edges[kr]->twin() << std::endl;
      std::cout << std::endl;
    }
#endif

  // Set twin pointers
  for (ki=0; ki<twin_faces.size(); ++ki)
    {
      // Check for coincidence
      // Quick test checking in one internal point
      int idx = -1;
      int nmb_found = 0;
      // Identity ident;
      shared_ptr<ParamSurface> srf1 = twin_faces[ki]->surface();
      vector<double> sf_dist(twinset.size());
      for (size_t kr=0; kr<twinset.size(); ++kr)
	  {
	    shared_ptr<ParamSurface> srf2 = twinset[kr]->surface();
#ifdef DEBUG_REG
	    std::ofstream of("twin_cand.g2");
	    srf1->writeStandardHeader(of);
	    srf1->write(of);
	    srf2->writeStandardHeader(of);
	    srf2->write(of);
#endif
	    
	    double upar1, vpar1;
	    Point pos1 = srf1->getInternalPoint(upar1, vpar1);

	    // Compute closest point in the other surface
	    double upar2, vpar2, dist;
	    Point pos2;
	    srf2->closestPoint(pos1, upar2, vpar2, pos2, dist,
			       toptol_.gap);
	    // int res = ident.identicalSfs(srf1, srf2, toptol_.neighbour);
	    // if (res > 0)
	    sf_dist[kr] = pos1.dist(pos2);
	    if (sf_dist[kr] < toptol_.neighbour)
	      {
		// Coincidence
		idx = (int)kr;
		nmb_found++;
	      }
	  }

	if (nmb_found == 1)
	  {
	    // One coincident face found
	    twin_faces[ki]->connectTwin(twinset[idx].get(), toptol_.neighbour);
#ifdef DEBUG_REG2
	    vector<shared_ptr<ftEdge> > edges = twin_faces[ki]->getAllEdges();
	    for (size_t kk=0; kk<edges.size(); ++kk)
	      std::cout << kk << " " << edges[kk].get() << " " << edges[kk]->twin() << std::endl;
	    std::cout << std::endl;

	    edges = twinset[idx]->getAllEdges();
	    for (size_t kk=0; kk<edges.size(); ++kk)
	      std::cout << kk << " " << edges[kk].get() << " " << edges[kk]->twin() << std::endl;
	    std::cout << std::endl;

#endif
	  }
	int stop_break;
	stop_break = 1;
    }
     
#ifdef DEBUG_REG
  std::ofstream of5("regularized_model.g2");
  for (size_t kv=0; kv<faces_.size(); ++kv)
    {
      shared_ptr<ParamSurface> sf = getSurface(kv);
      sf->writeStandardHeader(of5);
      sf->write(of5);
      }
#endif
      
}
	  
//===========================================================================
shared_ptr<ParamSurface> 
SurfaceModel::mergeAllSurfaces()
//===========================================================================
{
  shared_ptr<ParamSurface> result;
  // Try to merge all faces by performing append between appropriate
  // underlying surfaces
  // First fetch smooth joints
  double eps = getTolerances().gap;
  double bend = getTolerances().bend;
  bool changed = true;
  while (changed)
    {
#ifdef DEBUG
      std::ofstream of1("pre_merge_all.g2");
      for (size_t ki=0; ki<faces_.size(); ++ki)
	{
	  faces_[ki]->surface()->writeStandardHeader(of1);
	  faces_[ki]->surface()->write(of1);
	}
#endif
      changed = false;

      FaceConnectivityUtils<ftEdgeBase,ftSurface> connectivity;
      vector<ftEdgeBase*> vec;
      vector<shared_ptr<ftSurface> > faces = allFaces();
      connectivity.smoothEdges(faces, vec, bend);
      
      for (size_t ki=0; ki<vec.size(); ++ki)
	{
	  ftEdge *edg = vec[ki]->geomEdge();
	  if (edg == NULL)
	    continue;
	  vector<ftSurface*> adj_faces = edg->getAdjacentFaces();

	  if (adj_faces.size() == 2)
	    {
	      // Check if a merge is feasible
	      vector<shared_ptr<ftEdge> > common_edgs = 
		adj_faces[0]->getCommonEdges(adj_faces[1]);
	      shared_ptr<Vertex> vx1, vx2;
	      edg->getVertices(vx1, vx2);
	      int dir1, dir2;
	      double val1, val2;
	      bool at_start1, at_start2;
	      pair<Point,Point> co_par1, co_par2;
	      int found;  
	      if (common_edgs.size() > 1)
		found = 0;
	      else
		found = SurfaceModelUtils::mergeSituation(adj_faces[0], 
							  adj_faces[1],
							  vx1, vx2, 
							  dir1, val1, at_start1,
							  dir2, val2, at_start2,
							  co_par1, co_par2, eps);

	      if (found == 2)
		{
		  vector<Point> seam_joints;
		  shared_ptr<ftSurface> merged_face;
		  try {
		    merged_face = 
		      mergeFaces(adj_faces[0], dir1, val1, at_start1,
				 adj_faces[1], dir2, val2, at_start2, 
				 co_par1, co_par2, seam_joints);
		  }
		  catch (...)
		    {
		      merged_face.reset();
		    }
		  if (merged_face.get())
		    {
#ifdef DEBUG
		      std::ofstream of2("post_merge_all.g2");
		      for (size_t ki=0; ki<faces_.size(); ++ki)
			{
			  faces_[ki]->surface()->writeStandardHeader(of2);
			  faces_[ki]->surface()->write(of2);
			}
#endif
		      changed = true;
		      break;
		    }
		}
	    }
	  
	}
    }
  if (faces_.size() == 1)
    {
      result = faces_[0]->surface();
      shared_ptr<BoundedSurface> bd_res =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(result);
      if (bd_res.get())
	{
	  double max_dist = 0.0;
	  bd_res->simplifyBdLoops(eps, bend, max_dist);
	  int stop_break = 1;
	}
    }
  return result;
}

//===========================================================================
shared_ptr<ftSurface> 
SurfaceModel::mergeFaces(ftSurface* face1, int pardir1, double parval1,
			 bool atstart1, ftSurface* face2, int pardir2, 
			 double parval2, bool atstart2,
			 pair<Point,Point> co_par1, pair<Point,Point> co_par2,
			 vector<Point>& seam_joints)
//===========================================================================
{
  double eps = toptol_.gap;

  // Get surfaces and check consistency
  // The second surface is copied to avoid changing the input in case a merge
  // cannot be completed
  shared_ptr<ftSurface> dummy;
  shared_ptr<ParamSurface> surf1 = face1->surface();
  shared_ptr<ParamSurface> surf2 = shared_ptr<ParamSurface>(face2->surface()->clone());

  shared_ptr<BoundedSurface> bd_sf1 = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
  if (!bd_sf1.get())
    {
      // @@@ VSK, 022012. Creating bounded surfaces is not an ideal solution. Better would
      // be to merge spline surfaces directly, but this is a quick experiment to
      // be able to use existing code
      shared_ptr<SplineSurface> spl = 
	dynamic_pointer_cast<SplineSurface,ParamSurface>(surf1);
      if (spl.get())
	bd_sf1 = shared_ptr<BoundedSurface>(BoundedUtils::convertToBoundedSurface(*spl, eps));
    }
  shared_ptr<BoundedSurface> bd_sf2_0 = 
    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
  if (!bd_sf2_0.get())
    {
      shared_ptr<SplineSurface> spl = 
	dynamic_pointer_cast<SplineSurface,ParamSurface>(surf2);
      if (spl.get())
	bd_sf2_0 = shared_ptr<BoundedSurface>(BoundedUtils::convertToBoundedSurface(*spl, eps));
    }
   if (!(bd_sf1.get() && bd_sf2_0.get()))
    return dummy;

  // Ensure consistent parameter curves across joint by changing the parameter domain
  // of surface two. Copy initial surface
  shared_ptr<BoundedSurface> bd_sf2 = 
    shared_ptr<BoundedSurface>(bd_sf2_0->clone());
  double a1 = (pardir1 == 0) ? co_par1.first[1] : co_par1.first[0];
  double a2 = (pardir1 == 0) ? co_par2.first[1] : co_par2.first[0];
  double b1 = (pardir2 == 0) ? co_par1.second[1] : co_par1.second[0];
  double b2 = (pardir2 == 0) ? co_par2.second[1] : co_par2.second[0];

  // TEST
  if (false /*(a2-a1)*(b2-b1) < 0.0*/)
    {
#ifdef DEBUG_REG
      std::ofstream pts("pts.g2");
      pts << "400 1 0 4 255 0 0 255" << std::endl;
      pts << "1 " << std::endl;
      pts << face1->point(co_par1.first[0], co_par1.first[1]) << std::endl;
      pts << "400 1 0 4 0 255 0 255" << std::endl;
      pts << "1 " << std::endl;
      pts << face1->point(co_par2.first[0], co_par2.first[1]) << std::endl;
      pts << "400 1 0 4 155 100 0 255" << std::endl;
      pts << "1 " << std::endl;
      pts << face2->point(co_par1.second[0], co_par1.second[1]) << std::endl;
      pts << "400 1 0 4 0 100 155  255" << std::endl;
      pts << "1 " << std::endl;
      pts << face2->point(co_par2.second[0], co_par2.second[1]) << std::endl;
#endif
      // Opposite direction of surfaces
      if (a2 < a1)
	{
	  RectDomain doma = bd_sf1->underlyingSurface()->containingDomain();
	  double a3 = (pardir1 == 0) ? doma.vmin() : doma.umin();
	  double a4 = (pardir1 == 0) ? doma.vmax() : doma.umax();
	  bd_sf1->reverseParameterDirection(pardir1 == 1);
	  a1 = a3 + (a4 - a1);
	  a2 = a3 + (a4 - a2);
	}
      if (b2 < b1)
	{
	  RectDomain domb = bd_sf2->underlyingSurface()->containingDomain();
	  double b3 = (pardir2 == 0) ? domb.vmin() : domb.umin();
	  double b4 = (pardir2 == 0) ? domb.vmax() : domb.umax();
	  bd_sf2->reverseParameterDirection(pardir2 == 1);
	  b1 = b3 + (b4 - b1);
	  b2 = b3 + (b4 - b2);
	}
    }
  // END TEST

  double c1, c2, d1, d2;
  RectDomain dom0 = bd_sf2->underlyingSurface()->containingDomain();
  if (pardir2 == 0)
    {
      c1 = dom0.umin();
      c2 = dom0.umax();
      d1 = a1 - (b1 - dom0.vmin())*(a2 - a1)/(b2 - b1);
      d2 = a1 - (b1 - dom0.vmax())*(a2 - a1)/(b2 - b1);
      if (d1 > d2)
	std::swap(d1, d2);
    }
  else
    {
      c1 = a1 - (b1 - dom0.umin())*(a2 - a1)/(b2 - b1);
      c2 = a1 - (b1 - dom0.umax())*(a2 - a1)/(b2 - b1);
      d1 = dom0.vmin();
      d2 = dom0.vmax();
      if (c1 > c2)
	std::swap(c1, c2);
    }
  bd_sf2->setParameterDomain(c1, c2, d1, d2);
  RectDomain bd_dom1 = bd_sf1->containingDomain();
  RectDomain bd_dom2 = bd_sf2->containingDomain();
 
  shared_ptr<ParamSurface> base1 = bd_sf1->underlyingSurface();
  shared_ptr<ParamSurface> base2 = bd_sf2->underlyingSurface();
  // shared_ptr<SplineSurface> base1 = 
  //   dynamic_pointer_cast<SplineSurface,ParamSurface>(bd_sf1->underlyingSurface());
  // shared_ptr<SplineSurface> base2 =  
  //   dynamic_pointer_cast<SplineSurface,ParamSurface>(bd_sf2->underlyingSurface());
  if (!(base1.get() && base2.get()))
      return dummy;

  // Get domains
  RectDomain dom1 = base1->containingDomain();
  RectDomain dom2 = base2->containingDomain();

  // Define sub domains
  // @@@ VSK, The underlying surfaces does not necessarily have the
  // same domain. It might be necessary to adjust
  double umin1, umin2, umax1, umax2, vmin1, vmin2, vmax1, vmax2;
  if (pardir1 == 0)
    {
      umin1 = (atstart1) ? parval1 : 
	bd_dom1.umin() - 0.1*(bd_dom1.umin()-dom1.umin());
      umax1 = (atstart1) ? bd_dom1.umax() + 0.1*(dom1.umax()-bd_dom1.umax()) : 
	parval1;
      vmin1 = dom1.vmin();
      vmax1 = dom1.vmax();
    }
  else
    {
      umin1 = dom1.umin();
      umax1 = dom1.umax();
      vmin1 = (atstart1) ? parval1 : 
	bd_dom1.vmin() - 0.1*(bd_dom1.vmin()-dom1.vmin());
      vmax1 = (atstart1) ? bd_dom1.vmax() + 0.1*(dom1.vmax()-bd_dom1.vmax()) : 
	parval1;
    }

  if (pardir2 == 0)
    {
      umin2 = (atstart2) ? parval2 : 
	bd_dom2.umin() - 0.1*(bd_dom2.umin() - dom2.umin());
      umax2 = (atstart2) ? bd_dom2.umax() + 0.1*(dom2.umax() - bd_dom2.umax()) : 
	parval2;
      vmin2 = dom2.vmin();
      vmax2 = dom2.vmax();
      if (pardir1 == 0)
	{
	  vmin1 = vmin2 = std::max(vmin1, vmin2);
	  vmax1 = vmax2 = std::min(vmax1, vmax2);
	}
      else
	{
	  umin1 = vmin2 = std::max(umin1, vmin2);
	  umax1 = vmax2 = std::min(umax1, vmax2);
	}
    }
 else
    {
      umin2 = dom2.umin();
      umax2 = dom2.umax();
      vmin2 = (atstart2) ? parval2 : 
	bd_dom2.vmin() - 0.1*(bd_dom2.vmin()-dom2.vmin());
      vmax2 = (atstart2) ? bd_dom2.vmax() + 0.1*(dom2.vmax()-bd_dom2.vmax()) : 
	parval2;
      if (pardir1 == 0)
	{
	  vmin1 = umin2 = std::max(vmin1, umin2);
	  vmax1 = umax2 = std::min(vmax1, umax2);
	}
      else
	{
	  umin1 = umin2 = std::max(umin1, umin2);
	  umax1 = umax2 = std::min(umax1, umax2);
	}
     }

  // Check if the trimmed surfaces are contained within the domain of the
  // sub surfaces of the base surfaces
  if (umin1 > bd_dom1.umin()+eps || umax1 < bd_dom1.umax()-eps || 
      vmin1 > bd_dom1.vmin()+eps || vmax1 < bd_dom1.vmax()-eps)
    return dummy;
  if (umin2 > bd_dom2.umin()+eps || umax2 < bd_dom2.umax()-eps || 
					    vmin2 > bd_dom2.vmin()+eps || vmax2 < bd_dom2.vmax()-eps)
    return dummy;

  // Make sub surfaces of the base surfaces
  // shared_ptr<SplineSurface> sub_base1 =
  //   shared_ptr<SplineSurface>(base1->subSurface(umin1, vmin1, umax1, vmax1));
  // shared_ptr<SplineSurface> sub_base2 =
  //   shared_ptr<SplineSurface>(base2->subSurface(umin2, vmin2, umax2, vmax2));
  vector<shared_ptr<ParamSurface> > tmp_sub1 = 
    base1->subSurfaces(umin1, vmin1, umax1, vmax1);
  vector<shared_ptr<ParamSurface> > tmp_sub2 = 
    base2->subSurfaces(umin2, vmin2, umax2, vmax2);
  if (tmp_sub1.size() != 1 || tmp_sub2.size() != 1)
    return dummy;
  shared_ptr<ParamSurface> sub_base1 = tmp_sub1[0];
  shared_ptr<ParamSurface> sub_base2 = tmp_sub2[0];

  // Make bounded surfaces, reusing the boundary loops
  vector<CurveLoop> loops1 = bd_sf1->allBoundaryLoops();
  vector<CurveLoop> loops2 = bd_sf2->allBoundaryLoops();
  vector<vector<shared_ptr<CurveOnSurface> > > vec1(loops1.size());
  vector<vector<shared_ptr<CurveOnSurface> > > vec2(loops2.size());
  int ki, kj;
  for (ki=0; ki<(int)loops1.size(); ++ki)
    {
      for (kj=0; kj<loops1[ki].size(); ++kj)
	{
	  shared_ptr<ParamCurve> tmp = 
	    shared_ptr<ParamCurve>(loops1[ki][kj]->clone());
	  shared_ptr<CurveOnSurface> tmp2 = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	  tmp2->setUnderlyingSurface(sub_base1);
	  vec1[ki].push_back(tmp2);
	}
    }
  for (ki=0; ki<(int)loops2.size(); ++ki)
    {
      for (kj=0; kj<loops2[ki].size(); ++kj)
	{
	  shared_ptr<ParamCurve> tmp = 
	    shared_ptr<ParamCurve>(loops2[ki][kj]->clone());
	  shared_ptr<CurveOnSurface> tmp2 = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp);
	  tmp2->setUnderlyingSurface(sub_base2);
	  vec2[ki].push_back(tmp2);
	}
    }

  shared_ptr<BoundedSurface> sub1 = 
    shared_ptr<BoundedSurface>(new BoundedSurface(sub_base1, vec1, toptol_.gap));
  shared_ptr<BoundedSurface> sub2 = 
     shared_ptr<BoundedSurface>(new BoundedSurface(sub_base2, vec2, toptol_.gap));
  if (!sub1.get() || !sub2.get())
    return dummy;

  // Prepare for append
  int reverse = 0;
  if (atstart1 && (!atstart2))
    {
      std::swap(sub1, sub2);
      std::swap(atstart1, atstart2);
      std::swap(pardir1, pardir2);
    }

  if (atstart1)
    {
      sub1->reverseParameterDirection(pardir1 == 0);
      reverse = pardir1 + 1;
    }
  if (!atstart2)
    sub2->reverseParameterDirection(pardir2 == 0);
  if (pardir1 != pardir2)
    sub2->swapParameterDirection();
  RectDomain dom3 = sub1->underlyingSurface()->containingDomain();
  RectDomain dom4 = sub2->underlyingSurface()->containingDomain();
  if (pardir1 == 0)
    sub2->setParameterDomain(dom3.umax(),  dom3.umax()+dom4.umax()-dom4.umin(), 
			     dom4.vmin(), dom4.vmax());
  else
    sub2->setParameterDomain(dom4.umin(), dom4.umax(),
			     dom3.vmax(), dom3.vmax()+dom4.vmax()-dom4.vmin());
	
  // Finally, check direction in joint
  // double u1 = (pardir1 == 0) ? umax1 : 0.5*(umin1 + umax1);
  // double v1 = (pardir1 == 0) ? 0.5*(vmin1 + vmax1) : vmax1;
  // double u2 = (pardir2 == 0) ? umin2 : 0.5*(umin2 + umax2);
  // double v2 = (pardir2 == 0) ? 0.5*(vmin2 + vmax2) : vmin2;
  // if (pardir1 != pardir2)
  //   std::swap(u2, v2);
  // vector<Point> pts1 = sub1->ParamSurface::point(u1, v1, 1);
  // vector<Point> pts2 = sub2->ParamSurface::point(u2, v2, 1);
  double sc = (a2-a1)*(b2-b1);  //(pardir1 == 0) ? pts1[2]*pts2[2] : pts1[1]*pts2[1];
  if (sc < 0)
    sub2->reverseParameterDirection(pardir1 == 1);
  
#ifdef DEBUG_REG
  std::ofstream of0("merge_before.g2");
  sub1->writeStandardHeader(of0);
  sub1->write(of0);
  sub2->writeStandardHeader(of0);
  sub2->write(of0);
#endif


  // Append
  // shared_ptr<SplineSurface> base3 = 
  //   dynamic_pointer_cast<SplineSurface, ParamSurface>(sub1->underlyingSurface());
  // shared_ptr<SplineSurface> base4 = 
  //   dynamic_pointer_cast<SplineSurface, ParamSurface>(sub2->underlyingSurface());
  shared_ptr<ParamSurface> base3 = sub1->underlyingSurface();
  shared_ptr<ParamSurface> base4 = sub2->underlyingSurface();
  if (!base3.get() || !base4.get())
    return dummy;

#ifdef DEBUG_REG
  std::ofstream of("base_merge.g2");
  base3->writeStandardHeader(of);
  base3->write(of);
  base4->writeStandardHeader(of);
  base4->write(of);
#endif

  // Fetch boundary loops
  // loop1 = sub1->outerBoundaryLoop();
  // loop2 = sub2->outerBoundaryLoop();
  vector<CurveLoop> bd_loops1 = sub1->allBoundaryLoops();
  vector<CurveLoop> bd_loops2 = sub2->allBoundaryLoops();
  // int nmb1 = loop1.size();
  int nmb2 = (int)bd_loops2.size();

  // Scale the second underlying surface to get approximately the
  // same parameterization of the two surfaces
  double tanlen1 = SurfaceTools::estimateTangentLength(base3.get(), 
						       pardir1+1, false);
  double tanlen2 = SurfaceTools::estimateTangentLength(base4.get(), 
						       pardir1+1, true);
  double fac = tanlen2/tanlen1;
  double umin, umax, vmin, vmax;
  double dist;

  RectDomain base4_dom = base4->containingDomain();
  umin = base4_dom.umin(); //base4->startparam_u();
  umax = base4_dom.umax(); //base4->endparam_u();
  vmin = base4_dom.vmin(); //base4->startparam_v();
  vmax = base4_dom.vmax(); //base4->endparam_v();
  if (pardir1 == 0)
    umax = umin + fac*(umax - umin);
  else
    vmax = vmin + fac*(vmax - vmin);

  if (fabs(fac - 1.0) > 0.1)
    {
  bool scaled = true;
  for (kj=0; kj<nmb2; ++kj)
    {
      int nmb2_2 = bd_loops2[kj].size();
      vector<shared_ptr<ParamCurve> > scaled_crvs(nmb2_2);
      for (ki=0; ki<nmb2_2; ++ki)
	{
	  shared_ptr<CurveOnSurface> sf_cv = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(bd_loops2[kj][ki]);
	  if (!sf_cv.get())
	    {
	      scaled = false;
	      break;
	    }
	  shared_ptr<CurveOnSurface> tmp = 
	    shared_ptr<CurveOnSurface>(sf_cv->clone());
	  scaled = tmp->setDomainParCrv(umin, umax, vmin, vmax,
					base4_dom.umin(), base4_dom.umax(),
					base4_dom.vmin(), base4_dom.vmax());
					// base4->startparam_u(), base4->endparam_u(), 
					// base4->startparam_v(), base4->endparam_v());
	  if (!scaled)
	    break;
	  scaled_crvs[ki] = tmp;
	}
      if (scaled)
	{
	  CurveLoop tmp_loop(scaled_crvs, bd_loops2[kj].getSpaceEpsilon());
	  bd_loops2[kj].swap(tmp_loop);
	}
    }
      
  // Scale the parameter domain accordingly
  if (scaled)
    base4->setParameterDomain(umin, umax, vmin, vmax);
    }

  RectDomain dom5 = base4->containingDomain();
  int continuity = 1; //(base3->rational() || base4->rational()) ? 0 : 1;

  // Check for rational spline surfaces
  shared_ptr<SplineSurface> base3_2 =
    dynamic_pointer_cast<SplineSurface,ParamSurface>(base3);
  shared_ptr<SplineSurface> base4_2 =
    dynamic_pointer_cast<SplineSurface,ParamSurface>(base4);
  if ((base3_2.get() && base3_2->rational()) ||
      (base4_2.get() && base4_2->rational()))
      continuity = 0;
  // // Make copy
  // shared_ptr<ParamSurface> base3_tmp(base3->clone());
  // shared_ptr<ParamSurface> base4_tmp(base4->clone());
  // double adjust_wgt = base3->appendSurface(base4.get(), pardir1+1, 
  // 					   continuity, dist, false);
  // if (adjust_wgt > toptol_.neighbour)
  //   return dummy;   // The tolerance here is quite arbitrary

  shared_ptr<ParamSurface> joined_sf;
  try {
    joined_sf = 
      base3->getAppendSurface(base4.get(), pardir1+1, continuity, dist, false);
  }
  catch (...)
    {
      return dummy;
    }
  if (!joined_sf.get())
    return dummy;

#ifdef DEBUG_REG
  joined_sf->writeStandardHeader(of);
  joined_sf->write(of);
#endif

  if (dist > toptol_.gap)
    //if (dist > toptol_.neighbour)
    {
      if (continuity == 1)
	{
	  // Try with positional continuity
	  shared_ptr<ParamSurface> joined_sf = 
	    base3->getAppendSurface(base4.get(), pardir1+1, 0, dist, false);
	  // base3_tmp->appendSurface(base4_tmp.get(), pardir1+1, 0, dist, false);
	  if (dist > toptol_.neighbour)
	    return dummy;
	  // base3 = base3_tmp;
	}
      else if (dist > toptol_.neighbour)
	return dummy;
    }
  sub1->replaceSurf(joined_sf);
  
  // Check continuity
  int cont = 1; //base3->basis(pardir1).getMinContinuity();
  
  shared_ptr<ftSurface> merged_face = 
    performMergeFace(sub1->underlyingSurface(),
		     bd_loops1, bd_loops2, 
		     face1->getBody(), seam_joints, reverse,
		     cont);

  if (merged_face.get())
    {
      // Update topology structure
      shared_ptr<ftSurface> shrface1 = fetchAsSharedPtr(face1);
      shared_ptr<ftSurface> shrface2 = fetchAsSharedPtr(face2);
      removeFace(shrface1);
      removeFace(shrface2);
      append(merged_face);

#ifdef DEBUG_REG2
      vector<shared_ptr<ftEdge> > edges = merged_face->getAllEdges();
      for (size_t kr=0; kr<edges.size(); ++kr)
	std::cout << kr << " " << edges[kr].get() << " " << edges[kr]->twin() << std::endl;
      std::cout << std::endl;
#endif
    }

#ifdef DEBUG_REG
  std::ofstream mod("post_merge.g2");
  for (size_t k2=0; k2<faces_.size(); ++k2)
    {
      shared_ptr<ParamSurface> sf = faces_[k2]->surface();
      sf->writeStandardHeader(mod);
      sf->write(mod);
    }
#endif

  return merged_face;
}

//===========================================================================
shared_ptr<ftSurface> 
SurfaceModel::mergeSeamFaces(ftSurface* face1, ftSurface* face2, int pardir,
			     vector<Point>& seam_joints)
//===========================================================================
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
  shared_ptr<SplineSurface> base2 = 
    dynamic_pointer_cast<SplineSurface,ParamSurface>(base);
  if (!base2.get())
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

  // Fetch boundary loops
  CurveLoop loop1 = bd_sf1->outerBoundaryLoop();
    bd_sf1->outerBoundaryLoop();
  CurveLoop loop2 = bd_sf2->outerBoundaryLoop();
  int nmb1 = loop1.size();
  int nmb2 = loop2.size();

  // // Scale the second underlying surface to get approximately the
  // // same parameterization of the two surfaces
   int ki, kj;
  // double tanlen1 = estimateTangentLength(sub1.get(), pardir+1, true);
  // double tanlen2 = estimateTangentLength(sub2.get(), pardir+1, false);
  // double fac = tanlen2/tanlen1;
  // double u1, u2, v1, v2;
  // if (move_first)
  //   {
  //     u1 = umin1;
  //     u2 = umax1;
  //     v1 = vmin1;
  //     v2 = vmax1;
  //   }
  // else
  //   {
  //     u1 = umin2;
  //     u2 = umax2;
  //     v1 = vmin2;
  //     v2 = vmax2;
  //   }
  // if (pardir == 0)
  //   u2 = u1 + fac*(u2 - u1);
  // else
  //   v2 = v1 + fac*(v2 - v1);

  // // First scale the assiciated parameter curves 
  // bool scaled = true;
  // vector<shared_ptr<ParamCurve> > scaled_crvs(move_first ? nmb1 : nmb2);
  // if (move_first)
  //   {
  //     for (ki=0; ki<nmb1; ++ki)
  // 	{
  // 	  shared_ptr<CurveOnSurface> sf_cv = 
  // 	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(loop1[ki]);
  // 	  if (!sf_cv.get())
  // 	    {
  // 	      scaled = false;
  // 	      break;
  // 	    }
  // 	  shared_ptr<CurveOnSurface> tmp = 
  // 	    shared_ptr<CurveOnSurface>(sf_cv->clone());
  // 	  scaled = tmp->setDomainParCrv(u1, u2, v1, v2,
  // 					umin1, umax1, vmin1, vmax1);
  // 	  if (!scaled)
  // 	    break;
  // 	  scaled_crvs[ki] = tmp;
  // 	}
  //     if (scaled)
  // 	{
  // 	  CurveLoop tmp_loop(scaled_crvs, loop1.getSpaceEpsilon());
  // 	  loop1.swap(tmp_loop);
  // 	}
  //   }
  // else
  //   {
  //     for (ki=0; ki<nmb2; ++ki)
  // 	{
  // 	  shared_ptr<CurveOnSurface> sf_cv = 
  // 	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(loop2[ki]);
  // 	  if (!sf_cv.get())
  // 	    {
  // 	      scaled = false;
  // 	      break;
  // 	    }
  // 	  shared_ptr<CurveOnSurface> tmp = 
  // 	    shared_ptr<CurveOnSurface>(sf_cv->clone());
  // 	  scaled = tmp->setDomainParCrv(u1, u2, v1, v2,
  // 					umin2, umax2, vmin2, vmax2);
  // 	  if (!scaled)
  // 	    break;
  // 	  scaled_crvs[ki] = tmp;
  // 	}
  //     if (scaled)
  // 	{
  // 	  CurveLoop tmp_loop(scaled_crvs, loop2.getSpaceEpsilon());
  // 	  loop2.swap(tmp_loop);
  // 	}
  //    }
      
  // // Scale the parameter domain accordingly
  // if (scaled)
  //   sub2->setParameterDomain(u1, u2, v1, v2);

  double dist;
  sub2->appendSurface(sub1.get(), pardir+1, 1, dist, false);
#ifdef DEBUG_REG
  sub2->writeStandardHeader(of);
  sub2->write(of);
#endif

  if (dist > toptol_.gap)
    return dummy;

  // Make boundary curves of merged surface.
  // First adjust parameter curves to fit with the modified domain

  // Remove coincident curves
#ifdef DEBUG_REG
  // DEBUG. Draw
  std::ofstream space("space.g2");
  sub2->writeStandardHeader(space);
  sub2->write(space);
#endif
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
	      tmp->makeParameterCurve(toptol_.gap, par1, par2);
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
	      tmp->makeParameterCurve(toptol_.gap, par1, par2);
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
	  if ((pos1.dist(pos4) < toptol_.gap && pos2.dist(pos3) < toptol_.gap) ||
	      (pos1.dist(pos3) < toptol_.gap && pos2.dist(pos4) < toptol_.gap))
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
	      if (dist4 < dist3 && dist4  < toptol_.neighbour)
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
    shared_ptr<BoundedSurface>(new BoundedSurface(sub2, bd_cvs, toptol_.gap));
  merged->analyzeLoops();
   double merge_dist;
  bool success;
  success = merged->simplifyBdLoops(toptol_.gap, toptol_.kink, merge_dist);
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
  (void)merged_face->createInitialEdges(toptol_.gap, toptol_.kink);

  // Check if any joints has been removed. Fetch face vertices
  vector<shared_ptr<Vertex> > vx = merged_face->vertices();
  for (ki=0; ki<(int)joints.size();)
    {
      for (kj=0; kj<(int)vx.size(); ++kj)
	if (joints[ki].dist(vx[kj]->getVertexPoint()) < toptol_.gap)
	  break;
      if (kj < (int)vx.size())
	joints.erase(joints.begin()+ki);
      else
	ki++;
    }
  seam_joints = joints;

  shared_ptr<ftSurface> shrface1 = fetchAsSharedPtr(face1);
  shared_ptr<ftSurface> shrface2 = fetchAsSharedPtr(face2);
  removeFace(shrface1);
  removeFace(shrface2);
  append(merged_face);
  
  vector<shared_ptr<ftEdge> > edges = merged_face->getAllEdges();
#ifdef DEBUG_REG2
  for (size_t kr=0; kr<edges.size(); ++kr)
    std::cout << kr << " " << edges[kr].get() << " " << edges[kr]->twin() << std::endl;
  std::cout << std::endl;
#endif


#ifdef DEBUG_REG
  std::ofstream mod("post_merge.g2");
  for (size_t k2=0; k2<faces_.size(); ++k2)
    {
      shared_ptr<ParamSurface> sf = faces_[k2]->surface();
      sf->writeStandardHeader(mod);
      sf->write(mod);
    }
#endif

  return merged_face;
}		  

//===========================================================================
shared_ptr<ftSurface> 
SurfaceModel::mergeSeamCrvFaces(ftSurface* face1, ftSurface* face2, 
				vector<Point>& seam_joints)
//===========================================================================
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
  shared_ptr<ParamSurface> base2 = bd_sf2->underlyingSurface();
  shared_ptr<BoundedSurface> bd_sf1_copy;
  shared_ptr<BoundedSurface> bd_sf2_copy;
  // Make copies to rewind in case of now success
  bd_sf1_copy = shared_ptr<BoundedSurface>(bd_sf1->clone());
  bd_sf2_copy = shared_ptr<BoundedSurface>(bd_sf2->clone());

  if (base.get() != base2.get())
    {
      // Check for identical underlying elementary surfaces.
      bool same = 
	ElementaryUtils::sameElementarySurface(base.get(), base2.get(),
					       toptol_.gap, toptol_.kink);
      if (same)
	{
	  // Make sure that the underlying surface is large enough
	  BoundingBox box2 = bd_sf2->boundingBox();
	  double len = box2.low().dist(box2.high());

	  // Surface extension still missing

	  // Release all parameter curves belonging to surface 2
	  int ix1, ix2;
	  vector<CurveLoop> bd_loops = bd_sf2->allBoundaryLoops();
	  for (ix1=0; ix1<(int)bd_loops.size(); ++ix1)
	    {
	      int nmb_cvs = bd_loops[ix1].size();
	      for (ix2=0; ix2<nmb_cvs; ++ix2)
		{
		  shared_ptr<CurveOnSurface> sf_cv =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bd_loops[ix1][ix2]);
		  if (sf_cv)
		    sf_cv->unsetParameterCurve();
		}
	    }

	  // Replace underlying surface
	  bd_sf2->replaceSurf(base);

	  // Make new parameter curves
	  for (ix1=0; ix1<(int)bd_loops.size(); ++ix1)
	    {
	      int nmb_cvs = bd_loops[ix1].size();
	      for (ix2=0; ix2<nmb_cvs; ++ix2)
		{
		  shared_ptr<CurveOnSurface> sf_cv =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bd_loops[ix1][ix2]);		 
		  if (sf_cv)
		    {
		      bool OK = 
			sf_cv->ensureParCrvExistence(toptol_.gap);
		      if (!OK)
			break;
		    }
		}
	      if (ix2 < nmb_cvs)
		break;
	    }
	  if (ix1 < (int)bd_loops.size())
	    {
	      // Operation failed. Reset
	      bd_sf2 = bd_sf2_copy;
	      return dummy;
	    }
#ifdef DEBUG_REG
	  std::ofstream of1("geom_cvs.g2");
	  std::ofstream of2 ("par_cvs.g2");
	  vector<CurveLoop> bdl1 = bd_sf1->allBoundaryLoops();
	  for (ix1=0; ix1<(int)bdl1.size(); ++ix1)
	    {
	      int nmb_cvs = bdl1[ix1].size();
	      for (ix2=0; ix2<nmb_cvs; ++ix2)
		{
		  shared_ptr<CurveOnSurface> sf_cv =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdl1[ix1][ix2]);
		  sf_cv->spaceCurve()->writeStandardHeader(of1);
		  sf_cv->spaceCurve()->write(of1);
		  sf_cv->parameterCurve()->writeStandardHeader(of2);
		  sf_cv->parameterCurve()->write(of2);
		}
	    }
	  
	  vector<CurveLoop> bdl2 = bd_sf2->allBoundaryLoops();
	  for (ix1=0; ix1<(int)bdl2.size(); ++ix1)
	    {
	      int nmb_cvs = bdl2[ix1].size();
	      for (ix2=0; ix2<nmb_cvs; ++ix2)
		{
		  shared_ptr<CurveOnSurface> sf_cv =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdl2[ix1][ix2]);
		  sf_cv->spaceCurve()->writeStandardHeader(of1);
		  sf_cv->spaceCurve()->write(of1);
		  sf_cv->parameterCurve()->writeStandardHeader(of2);
		  sf_cv->parameterCurve()->write(of2);
		}
	    }
#endif
	}
      else
	return dummy;  // Different underlying surfaces
    }

  // Make boundary curves of merged surface.
  // Fetch boundary loops
  // CurveLoop loop1 = bd_sf1->outerBoundaryLoop();
  // CurveLoop loop2 = bd_sf2->outerBoundaryLoop();
  vector<CurveLoop> loops1 = bd_sf1->allBoundaryLoops();
  vector<CurveLoop> loops2 = bd_sf2->allBoundaryLoops();

  shared_ptr<ftSurface> merged_face = 
    performMergeFace(base, loops1, loops2, face1->getBody(), seam_joints, 0);
  if (!merged_face.get())
    {
      bd_sf1 = bd_sf1_copy;
      bd_sf2 = bd_sf2_copy;
      return merged_face;
    }

  shared_ptr<ftSurface> shrface1 = fetchAsSharedPtr(face1);
  shared_ptr<ftSurface> shrface2 = fetchAsSharedPtr(face2);
  removeFace(shrface1);
  removeFace(shrface2);
  append(merged_face);
  
#ifdef DEBUG_REG2
  vector<shared_ptr<ftEdge> > edges = merged_face->getAllEdges();
  for (size_t kr=0; kr<edges.size(); ++kr)
    std::cout << kr << " " << edges[kr].get() << " " << edges[kr]->twin() << std::endl;
  std::cout << std::endl;
#endif


#ifdef DEBUG_REG
  std::ofstream mod("post_merge.g2");
  for (size_t k2=0; k2<faces_.size(); ++k2)
    {
      shared_ptr<ParamSurface> sf = faces_[k2]->surface();
      sf->writeStandardHeader(mod);
      sf->write(mod);
    }
#endif

  return merged_face;
}		  

//===========================================================================
shared_ptr<ftSurface> 
SurfaceModel::performMergeFace(shared_ptr<ParamSurface> base,
			       vector<CurveLoop>& loops1,
			       vector<CurveLoop>& loops2,
			       Body *bd,
			       vector<Point>& seam_joints,
			       int reverse, int cont)
//===========================================================================
{
  shared_ptr<ftSurface> merged_face;
  Point par_eps = SurfaceTools::getParEpsilon(*base, toptol_.gap);
  double ptol = 0.5*(par_eps[0]+par_eps[1]);

  // Remove coincident curves
#ifdef DEBUG_REG
  // DEBUG. Draw
  std::ofstream space("space.g2");
  base->writeStandardHeader(space);
  base->write(space);
#endif
  int nmb1 = loops1[0].size();
  int nmb2 = loops2[0].size();
  int ki, kj, kr;
  vector<shared_ptr<CurveOnSurface> > bd_cvs;
  for (ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ParamCurve> cv = loops1[0][ki];
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

      // Make sure that the back pointer to the base surface is correct
      tmp->setUnderlyingSurface(base);

      // Check parameter curve
      bool samecrv = tmp->sameCurve(toptol_.gap);

      bd_cvs.push_back(tmp);
    }

  for (ki=0; ki<nmb2; ++ki)
    {
      shared_ptr<ParamCurve> cv = loops2[0][ki];
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

      // Make sure that the back pointer to the base surface is correct
      tmp->setUnderlyingSurface(base);

      // Check parameter curve
      bool samecrv = tmp->sameCurve(toptol_.gap);

      bd_cvs.push_back(tmp);
    }

  // Fetch estimate of accuracy
  double eps = toptol_.neighbour;
  eps = std::max(eps, loops1[0].getSpaceEpsilon());
  eps = std::max(eps, loops2[0].getSpaceEpsilon());

  // Combine loops, first remove duplicates
  bool found2 = false;
  bool found = false;
  for (ki=0; ki<nmb1;)
    {
      if (!bd_cvs[ki]->hasParameterCurve())
	break;
      Point pos1 = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->startparam());
      Point pos2 = bd_cvs[ki]->ParamCurve::point(bd_cvs[ki]->endparam());
      for (kj=nmb1; kj<nmb1+nmb2;)
	{
	  found = false;
	  if (!bd_cvs[kj]->hasParameterCurve())
	    break;
	  Point pos3 = bd_cvs[kj]->ParamCurve::point(bd_cvs[kj]->startparam());
	  Point pos4 = bd_cvs[kj]->ParamCurve::point(bd_cvs[kj]->endparam());
	  if ((pos1.dist(pos4) < toptol_.gap && pos2.dist(pos3) < toptol_.gap) ||
	      (pos1.dist(pos3) < toptol_.gap && pos2.dist(pos4) < toptol_.gap))
	    {
	      // Check also distances in the parameter domain
	      Point parpt1 = 
		bd_cvs[ki]->parameterCurve()->point(bd_cvs[ki]->startparam());
	      Point parpt2 = 
		bd_cvs[ki]->parameterCurve()->point(bd_cvs[ki]->endparam());
	      Point parpt3 = 
		bd_cvs[kj]->parameterCurve()->point(bd_cvs[kj]->startparam());
	      Point parpt4 = 
		bd_cvs[kj]->parameterCurve()->point(bd_cvs[kj]->endparam());
	      if (!((parpt1.dist(parpt4) < ptol && parpt2.dist(parpt3) < ptol) ||
		    (parpt1.dist(parpt3) < ptol && parpt2.dist(parpt4))))
		{
		  ++kj;
		  continue;
		}
	      
	      // Maybe some more checking?
	      bd_cvs.erase(bd_cvs.begin()+kj);
	      nmb2--;
	      bd_cvs.erase(bd_cvs.begin()+ki);
	      nmb1--;
	      found2 = found = true;
	      break;
	    }
	  kj++;
	}
      if (!found)
	ki++;
    }

  if (!found2)
    return merged_face;

  // Sort curves
  vector<vector<shared_ptr<CurveOnSurface> > > bd_cvs_vec;
  bd_cvs_vec.push_back(bd_cvs);
  for (kr=0; kr<(int)bd_cvs_vec.size(); ++kr)
    {
      vector<shared_ptr<CurveOnSurface> > bd_cvs2;
      for (ki=0; ki<(int)bd_cvs_vec[kr].size()-1; ++ki)
	{
	  Point pos = bd_cvs_vec[kr][ki]->ParamCurve::point(bd_cvs_vec[kr][ki]->endparam());
	  Point pos1 = bd_cvs_vec[kr][ki+1]->ParamCurve::point(bd_cvs_vec[kr][ki+1]->startparam());
	  Point pos2 = bd_cvs_vec[kr][ki+1]->ParamCurve::point(bd_cvs_vec[kr][ki+1]->endparam());
	  double dist1 = pos.dist(pos1);
	  double dist2 = pos.dist(pos2);
	  for (kj=ki+2; kj<(int)bd_cvs_vec[kr].size(); ++kj)
	    {
	      Point pos3 = bd_cvs_vec[kr][kj]->ParamCurve::point(bd_cvs_vec[kr][kj]->startparam());
	      Point pos4 = bd_cvs_vec[kr][kj]->ParamCurve::point(bd_cvs_vec[kr][kj]->endparam());
	      double dist3 = pos.dist(pos3);
	      double dist4 = pos.dist(pos4);
	      if (std::min(dist3,dist4) < std::min(dist1,dist2))
		{
		  if (dist4 < dist3 && dist4  < toptol_.neighbour)
		    {
		      bd_cvs_vec[kr][kj]->reverseParameterDirection();
		      std::swap(dist3, dist4);
		    }
		  std::swap(bd_cvs_vec[kr][ki+1], bd_cvs_vec[kr][kj]);
		  dist1 = dist3;
		  dist2 = dist4;
		}
	    }
	  if (dist1 > eps && dist2 > eps) 
	    {
	      if  (loops1.size() == 1 && loops2.size() == 1)
		{
		}
	      else
		{
		  bd_cvs2.push_back(bd_cvs_vec[kr][ki+1]);
		  bd_cvs_vec[kr].erase(bd_cvs_vec[kr].begin()+ki+1);
		  ki--;
		}
	    }
	}
      if (bd_cvs2.size() > 0)
	bd_cvs_vec.push_back(bd_cvs2);
    }

  if (bd_cvs_vec.size() > 1)
    {
      // Make sure that the outer loop is the first. Check orientation
      for (ki=0; ki<(int)bd_cvs_vec.size(); ++ki)
	if (LoopUtils::paramIsCCW(bd_cvs_vec[ki], eps, eps))
	  break;
      if (ki != 0)
	std::swap(bd_cvs_vec[0], bd_cvs_vec[ki]);
    }

  // Add inner trimming loops
  for (kr=1; kr<(int)loops1.size(); ++kr)
    {
      vector<shared_ptr<CurveOnSurface> > bd_cvs2;
      int size = loops1[kr].size();
      for (ki=0; ki<size; ++ki)
	{
	  shared_ptr<ParamCurve> tmp_cv = loops1[kr][ki];
	  shared_ptr<CurveOnSurface> tmp_bd =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(tmp_cv);
	  if (tmp_bd.get())
	    {
	      tmp_bd->setUnderlyingSurface(base);
	      bd_cvs2.push_back(tmp_bd);
	    }
	}
      bd_cvs_vec.push_back(bd_cvs2);
    }

  for (kr=1; kr<(int)loops2.size(); ++kr)
    {
      vector<shared_ptr<CurveOnSurface> > bd_cvs2;
      int size = loops2[kr].size();
      for (ki=0; ki<size; ++ki)
	{
	  shared_ptr<ParamCurve> tmp_cv = loops2[kr][ki];
	  shared_ptr<CurveOnSurface> tmp_bd =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(tmp_cv);
	  if (tmp_bd.get())
	    {
	      tmp_bd->setUnderlyingSurface(base);
	      bd_cvs2.push_back(tmp_bd);
	    }
	}
      bd_cvs_vec.push_back(bd_cvs2);
    }

#ifdef DEBUG_REG
  std::ofstream space2("space2.g2");
  std::ofstream par2("par2.g2");
  for (kr=0; kr<(int)bd_cvs_vec.size(); ++kr)
    for (ki=0; ki<(int)bd_cvs_vec[kr].size(); ++ki)
      {
	bd_cvs_vec[kr][ki]->spaceCurve()->writeStandardHeader(space2);
	bd_cvs_vec[kr][ki]->spaceCurve()->write(space2);
	bd_cvs_vec[kr][ki]->parameterCurve()->writeStandardHeader(par2);
	bd_cvs_vec[kr][ki]->parameterCurve()->write(par2);
      }
#endif

  // Fetch joint positions in outer boundary loop
  vector<Point> joints(bd_cvs_vec[0].size());
  for (ki=0; ki<(int)bd_cvs_vec[0].size(); ++ki)
    joints[ki] = bd_cvs_vec[0][ki]->ParamCurve::point(bd_cvs_vec[0][ki]->startparam());

  // Make new bounded surface
  shared_ptr<BoundedSurface> merged;
  merged = 
    shared_ptr<BoundedSurface>(new BoundedSurface(base, bd_cvs_vec, toptol_.gap));
  double merge_dist;
  bool success = false;
  if (cont >= 1)
    try {
      success = merged->simplifyBdLoops(toptol_.gap, toptol_.kink, merge_dist);
    }
    catch (...)
      {
	success = false;
      }
  if (reverse)
    merged->reverseParameterDirection(reverse == 1);

#ifdef DEBUG_REG
  std::ofstream merge("merge_sf.g2");
  merged->writeStandardHeader(merge);
  merged->write(merge);
#endif

  // Make face
  merged_face = shared_ptr<ftSurface>(new ftSurface(merged, -1));
  merged_face->setBody(bd);
  (void)merged_face->createInitialEdges(toptol_.gap, toptol_.kink);

  // Check if any joints has been removed. Fetch face vertices
  vector<shared_ptr<Vertex> > vx = merged_face->vertices();
  for (ki=0; ki<(int)joints.size();)
    {
      for (kj=0; kj<(int)vx.size(); ++kj)
	if (joints[ki].dist(vx[kj]->getVertexPoint()) < toptol_.gap)
	  break;
      if (kj < (int)vx.size())
	joints.erase(joints.begin()+ki);
      else
	ki++;
    }
  seam_joints = joints;

  return merged_face;
}

//===========================================================================
void
SurfaceModel::replaceRegularSurfaces()
//===========================================================================
{
  for (int ki=0; ki<(int)faces_.size(); ++ki)
    {
      shared_ptr<ftSurface> face = getFace(ki);
      Body *bd = face->getBody();

      shared_ptr<ParamSurface> surf = face->getUntrimmed(toptol_.gap,
							 toptol_.neighbour,
							 toptol_.bend);
      if (!surf.get())
	continue;  // Not a regular surface

      if (surf.get() == face->surface().get())
	continue;  // Surface not changed

      // Remove the face from the model
      bool performed = removeFace(face);
      if (!performed)
	continue;

      // Add the new surface to the model. First make face
      shared_ptr<ftSurface> face2 = 
	shared_ptr<ftSurface>(new ftSurface(surf, -1));
      face2->setBody(bd);
      append(face2);
      ki--;
    }
}

//===========================================================================
shared_ptr<ftSurface>
SurfaceModel::replaceRegularSurface(ftSurface *face, bool only_corner)
//===========================================================================
{
  shared_ptr<ftSurface> face3;

  Body *bd = face->getBody();
  shared_ptr<ftSurface> face2 = fetchAsSharedPtr(face); // Need the shared pointer
  if (!face2.get())
    return face3;  // Not in this surface model

  shared_ptr<ParamSurface> surf = face->getUntrimmed(toptol_.gap,
						     toptol_.neighbour,
						     toptol_.bend,
						     only_corner);
  if (!surf.get())
    return face3;  // Not a regular surface

  if (surf.get() == face->surface().get())
    return face3;  // Surface not changed

  // Remove the face from the model
  bool performed = removeFace(face2);
  if (!performed)
    return face3;

  // Add the new surface to the model. First make face
  face3 =  shared_ptr<ftSurface>(new ftSurface(surf, -1));
  face3->setBody(bd);
  append(face3);

  return face3;
}

//===========================================================================
void SurfaceModel::simplifyShell()
//===========================================================================
{
  // Merge two surfaces and update topology before the next instance is found
  bool changed = true;

  while (changed)
    {
      changed = false;

      // Collect edges across which the continuity is g1
      FaceConnectivityUtils<ftEdgeBase,ftFaceBase> connectivity;
      vector<ftEdgeBase*> vec;
      connectivity.smoothEdges(faces_, vec, toptol_.bend);

#ifdef DEBUG_REG
      if (vec.size() > 0)
	{
	  std::ofstream edgof("smoothedges.g2");
	  for (size_t ki=0; ki<vec.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> cv = vec[ki]->geomEdge()->geomCurve();
	      shared_ptr<CurveOnSurface> sf_cv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	      if (sf_cv.get())
		cv = sf_cv->spaceCurve();
	      cv->writeStandardHeader(edgof);
	      cv->write(edgof);
	    }
	}
#endif
      for (size_t ki=0; ki<vec.size(); ++ki)
	{
	  // Extract candidate face
	  ftSurface *face1 = vec[ki]->geomEdge()->face()->asFtSurface();
	  ftSurface *face2 = NULL;
	  if (vec[ki]->twin())
	    face2 = vec[ki]->twin()->geomEdge()->face()->asFtSurface();

#ifdef DEBUG_REG
	  std::ofstream of("merge_faces.g2");
	  face1->surface()->writeStandardHeader(of);
	  face1->surface()->write(of);	      
	  face2->surface()->writeStandardHeader(of);
	  face2->surface()->write(of);	      
#endif

	  // Check if a merge is possible
	  int dir1, dir2;
	  double val1, val2;
	  bool atstart1, atstart2;
	  pair<Point, Point> co_par1;
	  pair<Point, Point> co_par2;
	  int merge_cand = mergeSituation(face1,dir1, val1, atstart1,
					    face2, dir2, val2, atstart2, 
					   co_par1, co_par2);
	  if (merge_cand == 2)
	    {
	      vector<Point> seam_joints;
	      shared_ptr<ftSurface> merged = 
		mergeFaces(face1, dir1, val1, atstart1, face2, dir2, val2, 
			   atstart2, co_par1, co_par2, seam_joints);
	      if (merged.get())
		{
#ifdef DEBUG_REG
		  merged->surface()->writeStandardHeader(of);
		  merged->surface()->write(of);	      
#endif
		  changed = true;
		  break;
		}
	      else
		merge_cand = 1;
	    }
	  if (false/*merge_cand == 1*/)
	    {
	      // Try an alternative approach
	      // Make intermediate surface model
	      vector<shared_ptr<ParamSurface> > sfs(2);
	      sfs[0] = shared_ptr<ParamSurface>(face1->surface()->clone());
	      sfs[1] = shared_ptr<ParamSurface>(face2->surface()->clone());
	      shared_ptr<SurfaceModel> 
		tmp_model(new SurfaceModel(approxtol_,
					   toptol_.gap,
					   toptol_.neighbour,
					   toptol_.kink,
					   toptol_.bend,
					   sfs));

	      // Perform approximation
	      double error;
	      shared_ptr<ParamSurface> approx_surf = tmp_model->approxFaceSet(error);
	      if (approx_surf.get())
		{
		  // Replace the two original faces in this model with the new
		  // surface
		  shared_ptr<ftSurface> merged = 
		    shared_ptr<ftSurface>(new ftSurface(approx_surf, -1));
		  (void)merged->createInitialEdges(toptol_.gap, 
						   toptol_.kink);
		  // Update topology structure
		  shared_ptr<ftSurface> shrface1 = fetchAsSharedPtr(face1);
		  shared_ptr<ftSurface> shrface2 = fetchAsSharedPtr(face2);
		  removeFace(shrface1);
		  removeFace(shrface2);
		  append(merged);
#ifdef DEBUG_REG
		  merged->surface()->writeStandardHeader(of);
		  merged->surface()->write(of);	      
#endif
		  changed = true;
		  break;
		}
	    }
 	}
    }
}

//===========================================================================
int
SurfaceModel::mergeSituation(ftSurface* face1, int& dir1, double& val1, bool& atstart1,
			     ftSurface* face2, int& dir2, double& val2, bool& atstart2,
			     pair<Point, Point>& co_par1, pair<Point, Point>& co_par2)
//===========================================================================
{
  int can_merge = 2;
  double eps = toptol_.gap;

  // Fetch common edges
  vector<shared_ptr<ftEdge> > edges = face1->getCommonEdges(face2);
  if (edges.size() == 0)
    return 0;

  // Move startpoint of edge chain if necessary
  if (edges.size() > 1)
    {
      while (true)
	{
	  Point p1 = edges[0]->point(edges[0]->tMin());
	  Point p2 = edges[0]->point(edges[0]->tMax());
	  Point p3 = edges[edges.size()-1]->point(edges[edges.size()-1]->tMin());
	  Point p4 = edges[edges.size()-1]->point(edges[edges.size()-1]->tMax());
	  if (p1.dist(p4) < p2.dist(p3))
	    {
	      edges.insert(edges.begin(), edges[edges.size()-1]);
	      edges.pop_back();
	    }
	  else
	    break;
	}
    }
  // Check if the edges are connected and whether the edges are isoparametric in the
  // associated surface
  shared_ptr<ParamCurve> cv1 = edges[0]->geomCurve();
  shared_ptr<ParamCurve> cv2 = edges[0]->twin()->geomEdge()->geomCurve();
  shared_ptr<CurveOnSurface> sf_cv1 = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv1);
  shared_ptr<CurveOnSurface> sf_cv2 = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv2);
  if (!(sf_cv1.get() && sf_cv2.get()))
    return 0;

  bool iso1 = sf_cv1->isConstantCurve(eps, dir1, val1);
  bool iso2 = sf_cv2->isConstantCurve(eps, dir2, val2);
  if (!(iso1 && iso2))
    can_merge = 1;
  
  for (size_t ki=1; ki<edges.size(); ++ki)
    {
      tpJointType join = edges[ki-1]->checkContinuity(edges[ki].get(), toptol_.neighbour,
						      toptol_.gap, toptol_.bend,
						      toptol_.kink);
      if (join >= JOINT_DISC || join == JOINT_G0)
	return 0;

      shared_ptr<ParamCurve> cv3 = edges[ki]->geomCurve();
      shared_ptr<ParamCurve> cv4 = edges[ki]->twin()->geomEdge()->geomCurve();
      shared_ptr<CurveOnSurface> sf_cv3 = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv3);
      shared_ptr<CurveOnSurface> sf_cv4 = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv4);
      if (!(sf_cv4.get() && sf_cv4.get()))
	return 0;
      
      int dir3, dir4;
      double val3, val4;
      iso1 = sf_cv3->isConstantCurve(eps, dir3, val3);
      iso2 = sf_cv4->isConstantCurve(eps, dir4, val4);
      if (!(iso1 && iso2))
	can_merge = 1;
     
      else if (dir1 != dir3 || dir2 != dir4)
	return can_merge = 1; // Not the same iso parameter
      else if (fabs(val3-val1) > eps || fabs(val4-val2) > eps)
	return can_merge = 1;  // Not the same parameter value
    }

  // Check the tangency of adjacent boundary edges
  ftEdgeBase *next1 = edges[edges.size()-1]->next();
  ftEdgeBase *prev1 = edges[edges.size()-1]->twin()->prev();
  tpJointType join = prev1->checkContinuity(next1, toptol_.neighbour, toptol_.gap, 
					    toptol_.bend, toptol_.kink);
  if (join >= JOINT_DISC || join == JOINT_G0)
    return 0;
  join = edges[edges.size()-1]->checkContinuity(next1, toptol_.neighbour, 
						toptol_.gap, toptol_.bend, 
						toptol_.kink);
  if (join == JOINT_G1 || join == JOINT_KINK)
    return 0;
  join = prev1->checkContinuity(edges[edges.size()-1]->twin(), 
				toptol_.neighbour, toptol_.gap, toptol_.bend, 
				toptol_.kink);
  if (join == JOINT_G1 || join == JOINT_KINK)
    return 0;

  ftEdgeBase *prev2 = edges[0]->prev();
  ftEdgeBase *next2 = edges[0]->twin()->next();
  join = prev2->checkContinuity(next2, toptol_.neighbour, toptol_.gap, 
				toptol_.bend, toptol_.kink);
  if (join >= JOINT_DISC || join == JOINT_G0)
    return 0;
  join = edges[0]->twin()->checkContinuity(next2, toptol_.neighbour, 
					   toptol_.gap, toptol_.bend, 
					   toptol_.kink);
  if (join == JOINT_G1 || join == JOINT_KINK)
    return 0;
  join = prev2->checkContinuity(edges[0].get(), toptol_.neighbour, toptol_.gap, 
				toptol_.bend, toptol_.kink);
  if (join == JOINT_G1 || join == JOINT_KINK)
    return 0;
  
  if (can_merge == 2)
    {
      dir1--;
      dir2--;
      
      // Check whether the merge is at start- or end of the surfaces
      double u1, u2, v1, v2;
      face1->surface()->getInternalPoint(u1, v1);
      face2->surface()->getInternalPoint(u2, v2);
      atstart1 = (dir1 == 0) ? (u1 > val1) : (v1 > val1);
      atstart2 = (dir2 == 0) ? (u2 > val2) : (v2 > val2);

      // Extract information on co parameterization
      shared_ptr<Vertex> vx1 = edges[0]->getVertex(true);
      shared_ptr<Vertex> vx2 = edges[edges.size()-1]->getVertex(false);
      Point par1_1 = vx1->getFacePar(face1);
      Point par1_2 = vx1->getFacePar(face2);
      Point par2_1 = vx2->getFacePar(face1);
      Point par2_2 = vx2->getFacePar(face2);
      co_par1 = make_pair(par1_1, par1_2);
      co_par2 = make_pair(par2_1, par2_2);
    }
  
  return can_merge;
}

//===========================================================================
shared_ptr<SplineSurface> SurfaceModel::approxFaceSet(double& error, int degree)
//===========================================================================
{
  // Approximate the current face set by one non-trimmed spline surface
  // if possible
  shared_ptr<SplineSurface> dummy;  // Default result equals no result
  error = 2.0*toptol_.neighbour;

  // Triangulate surface set
  // First set the density of the triangulation
  // BoundingBox box = boundingBox();
  // double len = box.low().dist(box.high());
  // double fac = 100.0;
  // double density = 2.0*toptol_.neighbour; // len/fac;
  // shared_ptr<ftPointSet> triang = triangulate(density);

  // // Make sure that the triangulation corners are appropriate for
  // // parameterization
  // triang->checkAndUpdateTriangCorners();
  shared_ptr<ftPointSet> triang = shared_ptr<ftPointSet>(new ftPointSet()); 
  vector<shared_ptr<ftSurface> > faces = allFaces();
  SurfaceModelUtils::triangulateFaces(faces, triang, toptol_.gap);

  // Fetch outer boundary
  vector<shared_ptr<ftEdge> > edges = getBoundaryEdges();
  vector<ftEdge*> edg(edges.size());
  for (size_t kj=0; kj<edges.size(); ++kj)
    edg[kj] = edges[kj].get();

  // Collect edges into curve bounding the final surface
  vector<shared_ptr<ParamCurve> > curves1;
  vector<Point> joint_pts;
  Path::getEdgeCurves(edg, curves1, joint_pts, toptol_.gap,
		      toptol_.bend, false);
  if (curves1.size() != 4)
    return dummy;

  // Create initial surface as a Coons patch
  // Approximate boundary curves if necessary
  // Check and fix orientation
  vector<shared_ptr<ParamCurve> > curves2;
  makeCoonsBdCvs(curves1, toptol_.gap, degree, curves2);

#ifdef DEBUG
  std::ofstream of1("bd_cvs.g2");
  for (size_t kr=0; kr<curves2.size(); ++kr)
    {
      curves2[kr]->writeStandardHeader(of1);
      curves2[kr]->write(of1);
    }
#endif

  // Create surface
  CurveLoop boundary(curves2, toptol_.gap);
  shared_ptr<SplineSurface> init_surf(CoonsPatchGen::createCoonsPatch(boundary));

#ifdef DEBUG
  std::ofstream of2("init_coons.g2");
  init_surf->writeStandardHeader(of2);
  init_surf->write(of2);
#endif

  // // Identify corners
  // vector<int> corner;
  // vector<Point> corner_pnts(4);
  // for (int kii=0; kii<4; ++kii)
  //   corner_pnts[kii] = curves2[kii]->point(curves2[kii]->startparam());
  // triang->identifyBdPnts(corner_pnts, corner);

  // // Define seed point at the boundary for triangulation purposes
  // PointIter first = (*triang)[corner[0]];
  // triang->setFirst(first);
  // vector<PointIter> next = first->getNeighbours();
  // for (size_t kj=0; kj<next.size(); ++kj)
  //   if (next[kj]->isOnBoundary())
  //     {
  // 	triang->setSecond(next[kj]);
  // 	break;
  //     }

  // try {
  //   // Parameterize sample points
  //   error = AdaptSurface::parameterizePoints(init_surf, triang, corner);
  // } catch (...)
  //   {
      // Parameterization of points did not succed. Try to parameterize by
      // projection on the initial surface
      error = AdaptSurface::projectPoints(init_surf, triang);
    // }
  if (error < approxtol_)
    return init_surf;

  // Approximate sample points
  int maxiter = 3; //2; //5;
  double max_error2, mean_error;
  shared_ptr<SplineSurface> result = AdaptSurface::doApprox(init_surf, maxiter, 
							    triang, approxtol_, 
							    max_error2, mean_error);
  if (max_error2 < error)
    {
      error = max_error2;
      return result;
    }
  else
      return init_surf;
}

//===========================================================================
shared_ptr<ParamSurface> SurfaceModel::representAsOneSurface(double& dist, 
							     int degree)
//===========================================================================
{
  shared_ptr<ParamSurface> result;
  double ang_lim = 10.0*getTolerances().bend;
  double a_tol = 1.0e-8;

  vector<shared_ptr<ftEdge> > bd_edgs = getBoundaryEdges(0);
#ifdef DEBUG
  std::ofstream of("bd_edgs.g2");
  for (size_t ki=0; ki<bd_edgs.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = bd_edgs[ki]->geomCurve();
      shared_ptr<SplineCurve> tmp_spl(tmp_cv->geometryCurve());
      tmp_spl->writeStandardHeader(of);
      tmp_spl->write(of);
    }
#endif

  vector<ftEdge*> edg(bd_edgs.size());
  for (size_t kj=0; kj<bd_edgs.size(); ++kj)
    edg[kj] = bd_edgs[kj].get();

  // Collect edges into curve bounding the final surface
  vector<shared_ptr<ParamCurve> > curves1;
  vector<Point> joint_pts;
  Path::getEdgeCurves(edg, curves1, joint_pts, toptol_.gap,
		      toptol_.bend, false);

  // Make sure that all curves are pure spline curves
  vector<shared_ptr<SplineCurve> > curves2(curves1.size());
  for (size_t kj=0; kj<curves1.size(); ++kj)
    curves2[kj] = shared_ptr<SplineCurve>(curves1[kj]->geometryCurve());
  
  // Compute curve lengths
  vector<double> len(curves2.size());
  for (size_t kj=0; kj<curves2.size(); ++kj)
    len[kj] = curves2[kj]->estimatedCurveLength();

  // Copy current state
  vector<shared_ptr<ParamCurve> > trim_cvs(curves2.size());
  for (size_t kj=0; kj<trim_cvs.size(); ++kj)
    trim_cvs[kj] = shared_ptr<ParamCurve>(curves2[kj]->clone());
  
#ifdef DEBUG
  std::ofstream of3_0("bd_cvs0.g2");
  for (size_t ki=0; ki<curves2.size(); ++ki)
    {
      curves2[ki]->writeStandardHeader(of3_0);
      curves2[ki]->write(of3_0);
    }
#endif

  // Reduce the number of boundary curves to 4.
  vector<shared_ptr<ParamCurve> > removed_cvs;
  while (curves2.size() > 4)
    {
      // Compute average and minimum curve length
      // Recompute also curve lengths
      double av_len = 0.0;
      double min_len = std::numeric_limits<double>::max();
      len.resize(curves2.size());
      for (size_t kj=0; kj<curves2.size(); ++kj)
	{
	  len[kj] = curves2[kj]->estimatedCurveLength();
	  av_len += len[kj];
	  min_len = std::min(min_len, len[kj]);
	}
      av_len /= (int)curves2.size();
       
      // First try to remove curves by connecting adjacent curves
      bool removed = false;

      double fac_len = 0.2;
      //if (min_len > fac_len*av_len)
	removed = reduceCrvNmb(curves2, degree, removed_cvs);
      
      if (!removed)
	{
	  // For each candidate curve to remove, compute the potential 
	  // new corner and compute the distance to the candidate curve
	  int ix1, ix2;
	  vector<double> dist2cv(curves2.size(), 0.0);
	  for (ix1=0, ix2=2; ix1<(int)curves2.size(); 
	       ++ix1, ix2=(ix2+1)%((int)curves2.size()))
	    {
	      // Construct a new corner in the extension of the previous and
	      // next curve to the current curve
	      shared_ptr<ParamCurve> prev = curves2[ix1];
	      shared_ptr<ParamCurve> next = curves2[ix2];

	      vector<Point> der1(2), der2(2);
	      prev->point(der1, prev->endparam(), 1);
	      next->point(der2, next->startparam(), 1);
	      der1[1].normalize();
	      der2[1].normalize();
	      der2[1] *= -1;
	      
	      int ix3 = (ix1 + 1)%(int)curves2.size();
	      double det = (der1[1]*der1[1])*(der2[1]*der2[1]) - 
		(der1[1]*der2[1])*(der1[1]*der2[1]);
	      if (det > 1.0e-8)
		{
		  double s = ((der2[0]*der1[1]-der1[0]*der1[1])*(der1[1]*der2[1]) -
			      (der2[0]*der2[1] - der1[0]*der2[1])*(der1[1]*der1[1]))/det;

		  Point pos = der2[0]+s*der2[1];

		  // Compute distance
		  double tpar, dist;
		  Point clo_pt;
		  curves2[ix3]->closestPoint(pos, curves2[ix3]->startparam(),
					     curves2[ix3]->endparam(), tpar,
					     clo_pt, dist);
		  dist2cv[ix3] = dist;
		}
	      else
		dist2cv[ix3] = std::numeric_limits<double>::max();
	    }

	  // Identify a short curve with a small distance to the new
	  // corner
	  double min_len = len[0] + dist2cv[0];
	  int min_ix = 0;
	  for (size_t kj=1; kj<len.size(); ++kj)
	    {
	      if (len[kj] + dist2cv[kj] < min_len)
		{
		  min_len = len[kj] + dist2cv[kj];
		  min_ix = (int)kj;
		}
	    }

	  ix1 = (min_ix>0) ? min_ix-1 : (int)curves2.size()-1;
	  ix2 = (min_ix+1) % (int)curves2.size();
	  
	  // Construct a new corner in the extension of the previous and next
	  // curve to the curves identified for removal
	  shared_ptr<ParamCurve> prev = curves2[ix1];
	  shared_ptr<ParamCurve> next = curves2[ix2];
      
	  vector<Point> der1(3), der2(3);
	  prev->point(der1, prev->endparam(), 2);
	  next->point(der2, next->startparam(), 2);
	  der1[1].normalize();
	  der2[1].normalize();
	  der2[1] *= -1;
	  double d1len = der1[2].length();
	  double d2len = der2[2].length();
      
	  double s = ((der2[0]*der1[1]-der1[0]*der1[1])*(der1[1]*der2[1]) -
		      (der2[0]*der2[1] - der1[0]*der2[1])*(der1[1]*der1[1]))/
	    ((der1[1]*der1[1])*(der2[1]*der2[1]) - (der1[1]*der2[1])*(der1[1]*der2[1]));
	  double t = ((der2[0]*der1[1]-der1[0]*der1[1]) + 
		      s*(der1[1]*der2[1]))/(der1[1]*der1[1]);
	  Point pos1 = der2[0]+s*der2[1];
	  Point pos2 = der1[0]+t*der1[1];
	  Point pos = (d1len*pos1 + d2len*pos2)/(d1len+d2len);
#ifdef DEBUG
	  std::ofstream of2("corner.g2");
	  of2 << "400 1 0 4 255 0 0 255" << std::endl;
	  of2 << "1" << std::endl;
	  of2 << pos << std::endl;
#endif

	  // Create curves interpolating the endpoints and tangents of the adjacent
	  // curves and the new corner point
	  der2[1] *= -1;

	  // First curve
	  double dd1 = der1[0].dist(pos);
	  der1[1] *= (dd1/3.0);
	  vector<Point> points1(3);
	  points1[0] = der1[0];
	  points1[1] = der1[1];
	  points1[2] = pos;
	  vector<int> type1(3);
	  type1[0] = type1[2] = 1;
	  type1[1] = 2;
	  shared_ptr<SplineCurve> crv1 = 
	    SISLCurveInterface::interpolate(points1, type1, degree);

	  // Second curve
	  double dd2 = der2[0].dist(pos);
	  der2[1] *= (dd2/3.0);
	  vector<Point> points2(3);
	  points2[0] = pos;
	  points2[1] = der2[1];
	  points2[2] = der2[0];
	  vector<int> type2(3);
	  type2[0] = type2[2] = 1;
	  type2[1] = 3;
	  shared_ptr<SplineCurve> crv2 = 
	    SISLCurveInterface::interpolate(points2, type2, degree);
#ifdef DEBUG
	  std::ofstream of3("curve_pieces.g2");
	  crv1->writeStandardHeader(of3);
	  crv1->write(of3);
	  crv2->writeStandardHeader(of3);
	  crv2->write(of3);
#endif

	  // Merge curves
	  double dist1, dist2;
	  ParamCurve *tmp1 = crv1.get();
	  ParamCurve *tmp2 = next.get();
	  try {
	    prev->appendCurve(tmp1, 1, dist1);
	    crv2->appendCurve(tmp2, 1, dist2);
	  }
	  catch (...)
	    {
	      return result;  // No success
	    }
      
	  curves2[ix2] = crv2;
	  if (ix2 > ix1)
	    {
	      removed_cvs.insert(removed_cvs.end(), curves2.begin()+ix1+1,
				 curves2.begin()+ix2);
	      curves2.erase(curves2.begin()+ix1+1, curves2.begin()+ix2);
	    }
	  else
	    {
	      removed_cvs.insert(removed_cvs.end(), curves2.begin()+ix1+1,
				 curves2.end());
	      removed_cvs.insert(removed_cvs.end(), curves2.begin(),
				 curves2.begin()+ix2);
	      curves2.erase(curves2.begin()+ix1+1, curves2.end());
	      curves2.erase(curves2.begin(), curves2.begin()+ix2);
	    }
	}
#ifdef DEBUG
	  std::ofstream of3_2("bd_cvs.g2");
	  for (size_t ki=0; ki<curves2.size(); ++ki)
	    {
	      curves2[ki]->writeStandardHeader(of3_2);
	      curves2[ki]->write(of3_2);
	    }
#endif

	  int stop_break = 1;
	}

  //     // Identify concave corners in reduced curve sequence
  //     concave.clear();
  //     cvder1 = 
  // 	curves2[curves2.size()-1]->ParamCurve::point(curves2[curves2.size()-1]->endparam(), 1);
  //     for (size_t kj=0; kj<curves2.size(); ++kj)
  // 	{
  // 	  vector<Point> cvder2 = curves2[kj]->ParamCurve::point(curves2[kj]->startparam(), 1);
  // 	  //cvder2[1] *= -1;
  // 	  double angle = cvder1[1].angle(cvder2[1]);

  // 	  // Point close;
  // 	  // double par[2];
  // 	  // double dist;
  // 	  // int idx;
  // 	  // closestPoint(cvder1[0], close, idx, par, dist);
  // 	  // Point norm1 = faces_[idx]->normal(par[0],par[1]);
  // 	  // Point norm2 = cvder2[1].cross(cvder1[1]);
  // 	  // if (norm1*norm2 < 0.0)
  // 	  //   angle = 2.0*M_PI-angle;
  // 	  // if (angle > M_PI)
  // 	  if (fabs(0.5*M_PI-angle) > ang_lim &&
  // 	      cvder1[1]*cvder2[1] < -a_tol)
  // 	    concave.push_back(make_pair((kj==0)?(int)curves2.size()-1:(int)kj-1, (int)kj));
	  
  // 	  cvder1  = curves2[kj]->ParamCurve::point(curves2[kj]->endparam(), 1);
  // 	}

  //   }
  // if (concave.size() > 0)
  //   return result;  // Configuration not suitable for approximation by one
  // // surface
#ifdef DEBUG
  std::ofstream of_cvs("coons_bd_cvs.g2");
  for (size_t kj=0; kj<curves2.size(); ++kj)
    {
      curves2[kj]->writeStandardHeader(of_cvs);
      curves2[kj]->write(of_cvs);
    }
#endif

  // Create initial surface
  vector<shared_ptr<ParamCurve> > curves3(curves2.begin(), curves2.end());
  vector<shared_ptr<ParamCurve> > curves4;
  double maxbddist = makeCoonsBdCvs(curves3, toptol_.gap, degree, curves4);
  CurveLoop boundary(curves4, toptol_.gap);
  shared_ptr<SplineSurface> init_surf(CoonsPatchGen::createCoonsPatch(boundary));

#ifdef DEBUG
  std::ofstream of4("init_coons.g2");
  init_surf->writeStandardHeader(of4);
  init_surf->write(of4);
#endif

  // Update initial surface
  // First create triangulated sample points
  shared_ptr<ftPointSet> triang = shared_ptr<ftPointSet>(new ftPointSet()); 
  vector<shared_ptr<ftSurface> > faces = allFaces();
  SurfaceModelUtils::triangulateFaces(faces, triang, toptol_.gap);

  // Erase boundary curve flag for points on removed boundaries
  triang->removeOuterBoundaryFlag(removed_cvs, toptol_.neighbour);

  // Parameterize by projection on the initial surface
  double error = AdaptSurface::projectPoints(init_surf, triang);

  shared_ptr<SplineSurface> surf;
  if (error > approxtol_)
    {
      // Approximate sample points
      int maxiter = 3; //2; //5;
      double max_error2, mean_error;
      surf = AdaptSurface::doApprox(init_surf, maxiter, triang, approxtol_, 
				    max_error2, mean_error);
      if (max_error2 > error)
	{
	  surf = init_surf;
	}
      else
	{
	  error = max_error2;
	  if (mean_error > approxtol_ || 
	      (max_error2 > 10.0*approxtol_ && max_error2 > maxbddist))
	    surf.reset();  // Don't accept huge errors
	}
    }
  else 
    surf = init_surf;

  if (!surf.get())
    return result;

#ifdef DEBUG
  std::ofstream of5("updated_sf.g2");
  surf->writeStandardHeader(of5);
  surf->write(of5);
#endif

  if (removed_cvs.size() > 0)
    {
      // Make trimmed surface
      // Maybe better to use updated curve, trim with the removed ones and set
      // boundary information
      vector<shared_ptr<CurveOnSurface> > trim_cvs2(trim_cvs.size());
      for (size_t kj=0; kj<trim_cvs.size(); ++kj)
	{
	  trim_cvs2[kj] = 
	    shared_ptr<CurveOnSurface>(new CurveOnSurface(surf, trim_cvs[kj], false));
	  bool found = trim_cvs2[kj]->ensureParCrvExistence(toptol_.gap);
	  if (!found)
	    found = trim_cvs2[kj]->ensureParCrvExistence(toptol_.neighbour);
	}
  
      result = shared_ptr<ParamSurface>(new BoundedSurface(surf, trim_cvs2, 
							   toptol_.gap, false));
    }
  else
    result = surf;

#ifdef DEBUG
  std::ofstream of6("result_sf.g2");
  result->writeStandardHeader(of6);
  result->write(of6);
#endif

  return result;
}

//===========================================================================
  bool SurfaceModel::reduceCrvNmb(vector<shared_ptr<SplineCurve> >& curves,
				  int degree,
				  vector<shared_ptr<ParamCurve> >& removed_curves)
  //===========================================================================
  {
    double ang_lim1 = 0.35*M_PI;  // Treshhold for candidate edge to remove
    double ang_lim2 = 0.1;        // Treshhold for sharp edge
    double eps = getTolerances().gap;

    bool removed = true;
    while (removed && curves.size() > 4)
      {
	removed = false;
	// Compute angle of tangent cone for each curve
	vector<double> cone_ang(curves.size());
	int nmb_large = 0;
	for (size_t ki=0; ki<curves.size(); ++ki)
	  {
	    DirectionCone cone = curves[ki]->directionCone();
	    if (cone.greaterThanPi())
	      cone_ang[ki] = 1.5*M_PI;
	    else
	      cone_ang[ki] = cone.angle();
	    if (cone_ang[ki] >= ang_lim1)
	      nmb_large++;
	  }

	double max_ang = cone_ang[0];
	double prev_ang = max_ang;
	size_t ix = 0;
	while (!removed && ix < cone_ang.size())
	  {
	    // Find best candidate
	    for (size_t ki=ix+1; ki<cone_ang.size(); ++ki)
	      {
		if (cone_ang[ki] > max_ang)
		  {
		    prev_ang = max_ang;
		    max_ang = cone_ang[ki];
		    ix = ki;
		  }
		else if (cone_ang[ki] > prev_ang)
		  prev_ang = cone_ang[ki];
	      }
	    
	    if (max_ang < ang_lim1)
	      break;

	    if (max_ang - prev_ang < 0.2*max_ang); //0.1*max_ang)
	      break;  // No unique candidate

	    // Check if the curve can be bypassed
	    // Find indices of adjacent curves
	    size_t ix1 = (ix > 0) ? ix-1 : curves.size()-1;
	    size_t ix2 = (ix < curves.size()-1) ? ix+1 : 0;

	    // Create candidate joining curve
	    vector<Point> der1(2), der2(2), der3(2), der4(2);
	    curves[ix1]->point(der1, curves[ix1]->endparam(), 1);
	    curves[ix2]->point(der2, curves[ix2]->startparam(), 1);
	    curves[ix]->point(der3, curves[ix]->startparam(), 1);
	    curves[ix]->point(der4, curves[ix]->endparam(), 1);
	    
	    double ang1 = der1[1].angle(der3[1]);
	    double ang2 = der2[1].angle(der4[1]);
	    ang1 = std::min(M_PI-ang1,ang1);
	    ang2 = std::min(M_PI-ang2,ang2);

	    if (ang1 < ang_lim2 && ang2 < ang_lim2)
	      {
		// Bypassing the curve will not lead to removal
		++ix;
		continue;
	      }

	    if (curves.size() == 5 && ang1 < ang2)
	      ang1 = 0;
	    else if (curves.size() == 5)
	      ang2 = 0;

	    // Create curve arc bypassing the current curve
	    vector<Point> points;
	    vector<int> type;
	    points.push_back(der1[0]);
	    type.push_back(1);
	    if (ang1 >= ang_lim2)
	      {
		points.push_back(der1[1]);
		type.push_back(2);
	      }
	    if (ang2 >= ang_lim2)
	      {
		points.push_back(der2[1]);
		type.push_back(3);
	      }
	    points.push_back(der2[0]);
	    type.push_back(1);
	    
	    shared_ptr<SplineCurve> arc = 
	      SISLCurveInterface::interpolate(points, type, 
					      degree);
	      
#ifdef DEBUG
	    std::ofstream of("arc.g2");
	    arc->writeStandardHeader(of);
	    arc->write(of);
#endif

	    // Check curve
	    DirectionCone cone = arc->directionCone();
	    if (cone.greaterThanPi() || cone.angle() > max_ang)
	      {
		++ix;
		continue; // Bypassing curve more complex
	      }

	    // Check if the curve lies outside the current face set
	    // Define point between the curves
	    Point pos1, pos2;
	    curves[ix]->point(pos1, 0.5*(curves[ix]->startparam() +
					  curves[ix]->endparam()));
	    arc->point(pos2, 0.5*(0.5*(arc->startparam() +
				       arc->endparam())));
	    Point mid = 0.5*(pos1 + pos2);

	    double mindist1 = std::numeric_limits<double>::max();
	    double mindist2 = std::numeric_limits<double>::max();
	    for (size_t ki=0; ki<faces_.size(); ++ki)
	      {
		double upar1, upar2, vpar1, vpar2, dist1, dist2, par2;
		Point close1, close2, invec;

		shared_ptr<ftSurface> face = getFace(ki);
		face->closestPoint(mid, upar1, vpar1, close1, dist1, eps);
		ftEdgeBase* edg =
		  face->closestBoundaryPoint(mid, invec, upar2, vpar2, 
					     close2, dist2, par2);
		mindist1 = std::min(mindist1, dist1);
		mindist2 = std::min(mindist2, dist2);
	      }
	    if (mindist1 < mindist2-eps)
	      {
		++ix;
		continue;  // Arc inside surface set
	      }

	    // Merge curves
	    removed = true;
	    double d1, d2;
	    if (ang1 >= ang_lim2)
	      {
		curves[ix1]->appendCurve(arc.get(), 1, d1);
		if (ang2 >= ang_lim2)
		  {
		    curves[ix1]->appendCurve(curves[ix2].get(), 1, d2);
		    curves.erase(curves.begin()+ix2);
		  }
	      }
	    else
	      {
		arc->appendCurve(curves[ix2].get(), 1, d2);
		curves[ix2] = arc;
	      }
	    removed_curves.push_back(curves[ix]);
	    curves.erase(curves.begin()+ix);
	  }
      }
    return removed;
  }

// //===========================================================================
// void SurfaceModel::makeCoonsBdCvs(vector<shared_ptr<ParamCurve> >& cvs1,
// 				  double tol,
// 				  vector<shared_ptr<ParamCurve> >& cvs2)
// //===========================================================================
// {
//   // Check input
//   if (cvs1.size() != 4)
//     return;

//   // The curves are expected to be given in a head-to-tail orientation.
//   // Change direction of the last two curves
//   cvs1[2]->reverseParameterDirection();
//   cvs1[3]->reverseParameterDirection();

//   cvs2.resize(cvs1.size());

//   // For each pair of curves, create approximative curves on the same knot vector
//   static int min_cont = 2; //1;
//   for (int ki=0; ki<2; ++ki)
//     {
//       vector<shared_ptr<ParamCurve> > init_cvs;
//       vector<BsplineBasis> crv_basis;
//       int kj;
//       for (kj=0; kj<=2; kj+=2)
// 	{
// 	  shared_ptr<SplineCurve> spcv;
// 	  spcv = dynamic_pointer_cast<SplineCurve,ParamCurve>(cvs1[ki+kj]);
// 	  if (!spcv.get())
// 	    {
// 	      shared_ptr<CurveOnSurface> sfcv =
// 		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cvs1[ki+kj]);
// 	      if (sfcv.get() && sfcv->isConstantCurve())
// 		spcv = dynamic_pointer_cast<SplineCurve,ParamCurve>(sfcv->spaceCurve());
// 	    }
// 	  if (spcv.get() && (spcv->basis().getMinContinuity() < min_cont ||
// 			     spcv->rational()))
// 	    spcv.reset();

// 	  init_cvs.push_back(cvs1[ki+kj]);
// 	  if (spcv.get())
// 	    {
// 	      crv_basis.push_back(spcv->basis());
// 	      cvs2[ki+kj] = spcv;
// 	    }
// 	}

//       // Approximate curves in the same spline space up to possible
//       // refinements
//       vector<shared_ptr<SplineCurve> > app_cvs;
//       vector<double> knots;
//       int order = 0;

//       if (crv_basis.size() > 0/* && max_basis < max_coef*/)
// 	{
// 	  // Define initial knot vector as the union of the existing
// 	  // knot vectors
// 	  double start = crv_basis[0].startparam();
// 	  double end = crv_basis[0].endparam();
// 	  order = crv_basis[0].order();
// 	  for (kj=1; kj<(int)crv_basis.size(); ++kj)
// 	    {
// 	      order = std::max(order, crv_basis[kj].order());
// 	      start += crv_basis[kj].startparam();
// 	      end += crv_basis[kj].endparam();
// 	    }
// 	  start /= (double)(crv_basis.size());
// 	  end /= (double)(crv_basis.size());
// 	  for (kj=0; kj<(int)crv_basis.size(); ++kj)
// 	    {
// 	      crv_basis[kj].rescale(start, end);
// 	      if (crv_basis[kj].order() < order)
// 		crv_basis[kj].increaseOrder(order);
// 	    }

// 	  GeometryTools::makeUnionKnots(crv_basis, tol, knots);
	  
// 	  // Check the distribution of knots in the union knot vector
// 	  int nmb_basis = (int)knots.size()-order;
// 	  double tdel = (knots[nmb_basis] - knots[order-1])/(double)nmb_basis;
// 	  tdel /= 5.0;
// 	  double prev = knots[order-1];
// 	  for (kj=order; kj<=nmb_basis; ++kj)
// 	    if (knots[kj] > prev)
// 	      {
// 		if (knots[kj] - prev < tdel)
// 		  break;
// 		else
// 		  prev = knots[kj];
// 	      }
// 	  if (kj <= nmb_basis)
// 	    {
// 	      // Bad knot distribution. Use at most one kept curve
// 	      if (crv_basis.size() == 1)
// 		{
// 		  crv_basis.clear();
// 		  knots.clear();
// 		}
// 	      else
// 		{
// 		  crv_basis.erase(crv_basis.begin()+1, crv_basis.end());
// 		  knots = crv_basis[0].getKnots();
// 		}
// 	    }
// 	}

//       if (knots.size() >  1)
// 	{
// 	  BsplineBasis init_basis((int)knots.size()-order, order, knots.begin());
// 	  app_cvs = AdaptSurface::curveApprox(&init_cvs[0], 
// 					      (int)init_cvs.size(),
// 					      init_basis, tol);
// 	}
//       else
// 	app_cvs = AdaptSurface::curveApprox(&init_cvs[0], 
// 					    (int)init_cvs.size(), tol);

	  
//       // Collect final curves
//       for (kj=0; kj<2; ++kj)
// 	{
// 	  if (!cvs2[ki+2*kj].get())
// 	    cvs2[ki+2*kj] = app_cvs[kj];
// 	}
//     }

//   // The curves are expected to be given in a head-to-tail orientation.
//   // Change direction of the last two curves
//   cvs2[2]->reverseParameterDirection();
//   cvs2[3]->reverseParameterDirection();
// }

//===========================================================================
double SurfaceModel::makeCoonsBdCvs(vector<shared_ptr<ParamCurve> >& cvs1,
				  double tol, int degree,
				  vector<shared_ptr<ParamCurve> >& cvs2)
//===========================================================================
{
  double maxdist = 0.0;

  // Check input
  if (cvs1.size() != 4)
    return maxdist;

  // The curves are expected to be given in a head-to-tail orientation.
  // Change direction of the last two curves
  cvs1[2]->reverseParameterDirection();
  cvs1[3]->reverseParameterDirection();

  cvs2.resize(cvs1.size());

  for (int ki=0; ki<4; ++ki)
    {
      // Evaluate start- and end points and derivatives
      //vector<Point> startpt(2), endpt(2);
      vector<Point> startpt(1), endpt(1);
      cvs1[ki]->point(startpt, cvs1[ki]->startparam(), 0);
      cvs1[ki]->point(endpt, cvs1[ki]->endparam(), 0);

      // double len = cvs1[ki]->estimatedCurveLength(10);
      // startpt[1].normalize();
      // startpt[1] *= len/3.0;
      // endpt[1].normalize();
      // endpt[1] *= len/3.0;

      // Approximate
      vector<shared_ptr<ParamCurve> > curr_cvs(1, cvs1[ki]);

      int max_iter = 5;
      double max_dist;
      shared_ptr<SplineCurve> appr_cv(CurveCreators::approxCurves(&(curr_cvs[0]),
								  &(curr_cvs[1]),
								  startpt, endpt,
								  tol, max_dist, 
								  max_iter, degree));
      if (max_dist > tol) {
	MESSAGE("Failed approximating within tolerance (" << tol <<
		"), using cv anyway. Dist: " << max_dist);
	maxdist = std::max(maxdist, max_dist);
      }
      cvs2[ki] = appr_cv;
    }
	  
  // The curves are expected to be given in a head-to-tail orientation.
  // Change direction of the last two curves
  cvs2[2]->reverseParameterDirection();
  cvs2[3]->reverseParameterDirection();

  return maxdist;
}

//===========================================================================
bool
SurfaceModel::checkShellTopology()
//===========================================================================
{
  bool isOK = true;
  size_t ki, kj, kr, kh;
  for (ki=0; ki<faces_.size(); ++ki)
    {
      bool faceOK = faces_[ki]->asFtSurface()->checkFaceTopology();
      if (!faceOK)
	isOK = false;
    }

  // if (!isOK)
  //   return isOK;

  for (ki=0; ki<boundary_curves_.size(); ++ki)
    {
      for (kh=0; kh<boundary_curves_[ki].size(); ++kh)
	{
	  if (boundary_curves_[ki][kh]->getFace())
	    {
	      std::cout << "Loop face set. Loop = " << boundary_curves_[ki][kh] << std::endl;
	      isOK = false;
	    }

	  vector<shared_ptr<ftEdgeBase> > edges = boundary_curves_[ki][kh]->getEdges();
	  for (kj=0; kj<edges.size(); ++kj)
	    {
	      ftFaceBase *curr = edges[kj]->geomEdge()->face();
	      for (kr=0; kr<faces_.size(); ++kr)
		if (curr == faces_[kr].get())
		  break;
	      if (kr >= faces_.size())
		{
		  std::cout << "Boundary loop inconsistency, edge = " << edges[ki];
		  std::cout << ", face = " << curr << std::endl;
		  isOK = false;
		}
	    }
	}
    }

  vector<shared_ptr<Vertex> > vx;
  getAllVertices(vx);
  for (ki=0; ki<vx.size(); ++ki)
  {
      const double epsgeo = toptol_.gap;
      bool vxOK = vx[ki]->checkVertexTopology(epsgeo);
      if (!vxOK)
	isOK = false;
  }

  return isOK;
}

 } // namespace Go
