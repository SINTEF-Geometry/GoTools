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

#include "GoTools/qualitymodule/FaceSetRepair.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/qualitymodule/QualityUtils.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftMessage.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/GapRemoval.h"


using std::make_pair;

namespace Go
{
    using namespace qualityUtils;


  //===========================================================================
  FaceSetRepair::FaceSetRepair(shared_ptr<SurfaceModel> sfmodel)
  //===========================================================================
    : ModelRepair(), vertex_update_(false), edges_update_(false)
  {
    sfmodel_ = sfmodel;
    quality_ = 
      shared_ptr<FaceSetQuality>(new FaceSetQuality(sfmodel->getTolerances(),
						    sfmodel->getApproximationTol()));

    results_ = quality_->getResults();
  }

  //===========================================================================
  FaceSetRepair::FaceSetRepair(shared_ptr<FaceSetQuality> quality)
  //===========================================================================
    : ModelRepair(quality->getResults()), quality_(quality), 
      vertex_update_(false), edges_update_(false)
  {
    sfmodel_ = quality->getAssociatedSfModel();
  }

  //===========================================================================
  FaceSetRepair::~FaceSetRepair()
  //===========================================================================
  {
  }

  //===========================================================================
  void FaceSetRepair::mendGaps()  
  //===========================================================================
  {
    // Check if quality results exist already
    ASSERT(quality_.get());
    ASSERT(results_.get());

    // Make sure that the gaps are computed and that the tolerance matches
    double epsge = 0.5*sfmodel_->getTolerances().gap;
    // double tol = -1.0;
    vector<pair<ftEdge*, ftEdge*> > pos_discont;
    quality_->facePositionDiscontinuity(pos_discont);

    // Make sure that the vertex positions are updated
    if (!vertex_update_)
      {
	std::set<shared_ptr<Vertex> > vertices;
	for (size_t ki=0; ki<pos_discont.size(); ++ki)
	  {
	    shared_ptr<Vertex>  vx1, vx2;
	    pos_discont[ki].first->getVertices(vx1, vx2);
	    vertices.insert(vertices.end(), vx1);
	    vertices.insert(vertices.end(), vx2);
	  }

	std::set<shared_ptr<Vertex> >::iterator it = vertices.begin();
	while (it != vertices.end())
	  {
	    improveVertexPos(*it, epsge);
	    ++it;
	  }
	
      }

	
    // First try to remove gaps by improving the trimming curves
    gapTrimming(pos_discont, epsge, false);

    // Search for cases where two B-spline surfaces meet in a common
    // boundary
    updateMeetingSplines(pos_discont, epsge);

    // One B-spline surfaces meet another surface along its boundary
    updateOneSpline(pos_discont, epsge);

    // Two B-spline surfaces
    modifyTwoSplines(pos_discont, epsge);

    // One B-spline surface and one surface that is not a B-spline surface
    modifyOneSpline(pos_discont, epsge);

    // Search for B-spline surfaces where the edge is an iso parameter
    // and this is not the case for the opposite surface
    updateOneSpline(pos_discont, epsge);

    //mendEdgeDistance();

    // Update data in face set
    sfmodel_->setTopology();

    // Update quality results
    pos_discont.clear();
    results_->reset(FACE_POSITION_DISCONT);
    quality_->facePositionDiscontinuity(pos_discont);

    // Average two B-spline surfaces meeting in a common
    // boundary
    averageMeetingSplines(pos_discont, epsge);

    // Try to remove remaining gaps by improving the trimming curves
    gapTrimming(pos_discont, epsge, true);

    // Update data in face set
    sfmodel_->setTopology();

    // Update quality results
    pos_discont.clear();
    results_->reset(FACE_POSITION_DISCONT);
    quality_->facePositionDiscontinuity(pos_discont);

 
  }

  //===========================================================================
  void FaceSetRepair::identicalAndEmbeddedFaces()  
  //===========================================================================
  {
    // Check if quality results exist already
    ASSERT(quality_.get());
    ASSERT(results_.get());

    // Fetch information about identical and embedded faces
    // Compute if necessary
    vector<pair<shared_ptr<ftSurface>,shared_ptr<ftSurface> > > identical;
    vector<pair<shared_ptr<ftSurface>,shared_ptr<ftSurface> > > embedded;

    quality_->identicalOrEmbeddedFaces(identical, embedded);

    // First remove embedded faces
    size_t ki;
    vector<shared_ptr<ftSurface> > other;
    for (ki=0; ki<embedded.size(); ++ki)
      {
	sfmodel_->removeFace(embedded[ki].second);
	other.push_back(embedded[ki].first);
      }

    // Identical faces. Remove the one with less number of neighbours
    for (ki=0; ki<identical.size(); ++ki)
      {
	vector<ftSurface*> neighbour1, neighbour2;
	identical[ki].first->getAdjacentFaces(neighbour1);
	identical[ki].second->getAdjacentFaces(neighbour2);
	if (neighbour1.size() < neighbour2.size())
	  {
	    sfmodel_->removeFace(identical[ki].first);
	    other.push_back(identical[ki].second);
	  }
	else
	  {
	    sfmodel_->removeFace(identical[ki].second);
	    other.push_back(identical[ki].first);
	  }
      }

    // Update topology related to remaining faces
    for (ki=0; ki<other.size(); ++ki)
      {
	sfmodel_->updateFaceTopology(other[ki]);
      }

    if (other.size() > 0)
      {
	// Update quality results
	identical.clear();
	embedded.clear();
	results_->reset(IDENTICAL_FACES);
	results_->reset(EMBEDDED_FACES);
	quality_->identicalOrEmbeddedFaces(identical, embedded);
      }
  }

  //===========================================================================
  void FaceSetRepair::identicalVertices()  
  //===========================================================================
  {
    // Check if quality results exist already
    ASSERT(quality_.get());
    ASSERT(results_.get());

    // Fetch information about identical vertices
    vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > vertex_pair;
    quality_->identicalVertices(vertex_pair);

    // Join info stored in identical vertices
    vector<shared_ptr<Vertex> > removed;
    for (size_t ki=0; ki<vertex_pair.size(); ++ki)
      {
	// Check if the first vertex in the pair is removed
	size_t kj;
	for (kj=0; kj<removed.size(); ++kj)
	  if (vertex_pair[ki].first.get() == removed[kj].get())
	    break;

	if (kj < removed.size())
	  {
	    // Join vertices. Keep the second
	    ftEdge* edge = vertex_pair[ki].second->getEdge(0);
	    edge->joinVertex(vertex_pair[ki].first, vertex_pair[ki].second);
	    removed.push_back(vertex_pair[ki].first);
	  }
	else
	  {
	    // Join vertices. Keep the first
	    ftEdge* edge = vertex_pair[ki].first->getEdge(0);
	    edge->joinVertex(vertex_pair[ki].second, vertex_pair[ki].first);
	    removed.push_back(vertex_pair[ki].second);
	  }
      }

    if (vertex_pair.size() > 0)
      {
	// Recompute vertex identity
	vertex_pair.clear();
	results_->reset(IDENTICAL_VERTICES);
	quality_->identicalVertices(vertex_pair);
      }
  }

  //===========================================================================
  void FaceSetRepair::consistentFaceNormal()  
  //===========================================================================
  {
    // Check if quality results exist already
    ASSERT(quality_.get());
    ASSERT(results_.get());

    // Fetch information about face normal inconsistency
    vector<shared_ptr<ftSurface> > faces;
    quality_->faceNormalConsistency(faces);
    std::cout << "Number of inconsistent faces: " << faces.size() << std::endl;

    int nmb_inconsistent = (int)faces.size();
    ftFaceBase* last_turned;
    vector<ftFaceBase*> turned;
    vector<ftFaceBase*> all_turned;
    while (faces.size() > 0)
      {
	// 	for (size_t ki=0; ki<faces.size(); ++ki)
	// 	  {
	// 	    int idx = sfmodel_->getIndex(faces[ki].get());
	// 	    sfmodel_->turn(idx);
	// 	  }
	size_t kj;
	for (kj=0; kj<faces.size(); ++kj)
	  {
	    size_t ki;
	    for (ki=0; ki<turned.size(); ++ki)
	      if (turned[ki] == faces[kj].get())
		break;
	    if (ki == turned.size())
	      break;
	  }

	ftFaceBase* curr_face = 0;
	if (kj == faces.size())
	  {
	    std::cout << "Cannot turn face " << std::endl;
	    break;
	    //	    kj = 0;
	    size_t ki, kh;
	    vector<pair<ftFaceBase*, ftFaceBase*> > candidate_faces;
	    sfmodel_->getInconsistentFacePairs(candidate_faces);
	    for (ki=0; ki<all_turned.size(); ++ki)
	      {
		for (kh=0; kh<candidate_faces.size(); ++kh)
		  if (candidate_faces[kh].first == all_turned[ki])
		    {
		      curr_face = candidate_faces[kh].second;
		      break;
		    }
		  else if (candidate_faces[kh].second == all_turned[ki])
		    {
		      curr_face = candidate_faces[kh].first;
		      break;
		    }
		if (kh < candidate_faces.size())
		  {
		    all_turned.erase(all_turned.begin()+ki);
		    break;
		  }
	      }
	  }

	int idx;
	if (curr_face)
	  {
	    idx = sfmodel_->getIndex(faces[kj].get());
	    sfmodel_->turn(idx);
	    last_turned = faces[kj].get();
// 	    turned.push_back(last_turned);
// 	    all_turned.push_back(last_turned);
	  }
	else
	  {
	    idx = sfmodel_->getIndex(faces[kj].get());
	    sfmodel_->turn(idx);
	    last_turned = faces[kj].get();
	    turned.push_back(last_turned);
	    all_turned.push_back(last_turned);
	  }

// 	std::ofstream out_file("curr_sfs.g2");
// 	int nmb_faces = sfmodel_->nmbEntities();
// 	for (int kr=0; kr<nmb_faces; ++kr)
// 	  {
// 	    shared_ptr<ParamSurface> sf = sfmodel_->getSurface(kr);
// 	    sf->writeStandardHeader(out_file);
// 	    sf->write(out_file);
// 	  }


	faces.clear();
	results_->reset(FACE_ORIENTATION);
	quality_->faceNormalConsistency(faces);
	std::cout << "Number of inconsistent faces: " << faces.size() << std::endl;
	if ((int)faces.size() >= nmb_inconsistent && !curr_face)
	  //if (faces.size() > nmb_inconsistent)
	  {
	    all_turned.erase(all_turned.begin() + all_turned.size() - 1);

	    // The last turn did not improve the situation
	    // Redo
	    sfmodel_->turn(idx);

	    // Find the face corresponding to the last one turned
	    vector<pair<ftFaceBase*, ftFaceBase*> > candidate_faces;
	    sfmodel_->getInconsistentFacePairs(candidate_faces);

	    for (kj=0; kj<candidate_faces.size(); ++kj)
	      if (candidate_faces[kj].first == last_turned ||
		  candidate_faces[kj].second == last_turned)
		break;
	    if (kj < candidate_faces.size())
	      {
		if (candidate_faces[kj].first == last_turned)
		  {
		    int idx2 = 
		      sfmodel_->getIndex(candidate_faces[kj].second->asFtSurface());
		    sfmodel_->turn(idx2);
		    last_turned = candidate_faces[kj].second;
		  }
		else
		  {
		    int idx2 = 
		      sfmodel_->getIndex(candidate_faces[kj].first->asFtSurface());
		    sfmodel_->turn(idx2);
		    last_turned = candidate_faces[kj].first;
		  }
		turned.push_back(last_turned);
		all_turned.push_back(last_turned);

// 		std::ofstream out_file("curr_sfs.g2");
// 		int nmb_faces = sfmodel_->nmbEntities();
// 		for (int kr=0; kr<nmb_faces; ++kr)
// 		  {
// 		    shared_ptr<ParamSurface> sf = sfmodel_->getSurface(kr);
// 		    sf->writeStandardHeader(out_file);
// 		    sf->write(out_file);
// 		  }

	      }

	faces.clear();
	results_->reset(FACE_ORIENTATION);
	quality_->faceNormalConsistency(faces);
	std::cout << "Number of inconsistent faces (2): " << faces.size() << std::endl;
	  }
	else
	  turned.clear();

	nmb_inconsistent = (int)faces.size();   
	
      }
  }

  //===========================================================================
  void FaceSetRepair::optimizeVertexPosition()
  //===========================================================================
  {
    // All vertices are iterated to get the best possible position with
    // regard to the faces meeting in the vertex
    // First fetch all vertices
    double fac = 0.1;
    double epsge = fac*sfmodel_->getTolerances().gap;
    vector<shared_ptr<Vertex> > vertices;
    sfmodel_->getAllVertices(vertices);

    // Modify position
    for (size_t ki=0; ki<vertices.size(); ++ki)
      {
	improveVertexPos(vertices[ki], epsge);

      }
    vertex_update_ = true;
  }

  //===========================================================================
  void FaceSetRepair::mendEdgeDistance()
  //===========================================================================
  {
    // Check if quality results exist already
    ASSERT(quality_.get());
    ASSERT(results_.get());
    double epsge = 0.5*sfmodel_->getTolerances().gap;

    // Fetch information about too distance face-edge relationships
    vector<pair<ftSurface*, ftEdge*> > edges;
    quality_->faceEdgeDistance(edges);

    for (size_t ki=0; ki<edges.size(); ++ki)
      {
	// Fetch associated vertices
	shared_ptr<Vertex> vx1, vx2;
	ftEdge *curr = edges[ki].second;
	curr->getVertices(vx1, vx2);

	// Check if the vertices are updated already
	if (!vertex_update_)
	  {
	    // Improve position of vertices
	    // 1. vertex
	    improveVertexPos(vx1, epsge);

	    // 2. vertex
	    improveVertexPos(vx2, epsge);
	  }
	
	// Regenerate edge
	shared_ptr<ParamCurve> bdcv = curr->geomCurve();
	shared_ptr<CurveOnSurface> sfcv = 
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bdcv);
	if (sfcv.get())
	  {
	    sfcv->updateCurves(vx1->getVertexPoint(), vx2->getVertexPoint(),
			       0.9*epsge);
	  }

      }
    edges_update_ = true;

    edges.clear();
    results_->reset(FACE_EDGE_DISTANCE);
    quality_->faceEdgeDistance(edges);
  }

  //===========================================================================
  void 
  FaceSetRepair::gapTrimming(vector<pair<ftEdge*, ftEdge*> >& pos_discont,
			     double epsge, bool update_iso)
  //===========================================================================
  {
    // Try to remove gaps by improving the trimming curves
    for (size_t ki=0; ki<pos_discont.size(); )
      {
	// Fetch the related faces
	ftSurface *face1 = pos_discont[ki].first->face()->asFtSurface();
	ftSurface *face2 = pos_discont[ki].second->face()->asFtSurface();

	if (!(face1 && face2))
	  continue;   // This should not happen. Cannot mend gaps

	// Check edge curves
	ftEdge* e1g = pos_discont[ki].first->geomEdge();
	ftEdge* e2g = pos_discont[ki].second->geomEdge();
	 shared_ptr<CurveOnSurface> sfcv1 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e1g->geomCurve());
	 shared_ptr<CurveOnSurface> sfcv2 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e2g->geomCurve());
	 if (!(sfcv1.get() && sfcv2.get()))
	     continue;

// 	// Update faces
// 	ftMessage status = face1->removeGap(pos_discont[ki].first,
// 					    pos_discont[ki].second,
// 					    face2, epsge);

	// Check if both faces are trimmed. In that case a better positioning
	// of the trim curve may remove the gap
	shared_ptr<BoundedSurface> bd1 = 
	  dynamic_pointer_cast<BoundedSurface, ParamSurface>(face1->surface());
	shared_ptr<BoundedSurface> bd2 = 
	  dynamic_pointer_cast<BoundedSurface, ParamSurface>(face2->surface());

	if (!(bd1.get() && bd2.get()))
	  {
	    ki++;
	    continue;   // Not handled in this attempt
	  }

	if (!update_iso)
	  {
	    // Check if the edge corresponds to a boundary curve for
	    // any of the surfaces
	    bool same1, same2;
	    int idx1 = sfcv1->whichBoundary(epsge, same1); 
	    int idx2 = sfcv2->whichBoundary(epsge, same2); 
	    if (idx1 >= 0 || idx2 >= 0)
	      {
		ki++;
		continue;  // Do not update
	      }
	  }
	    

	Point vertex1 = e1g->getVertex(true)->getVertexPoint();
	Point vertex2 = e1g->getVertex(false)->getVertexPoint();
	double max_gap = 
	  GapRemoval::removeGapTrim(sfcv1, e1g->tMin(), e1g->tMax(),
				    sfcv2, e2g->tMin(), e2g->tMax(),
				    vertex1, vertex2, epsge);
	if (max_gap < epsge)
	  {
	    // Gap removed.
	    pos_discont.erase(pos_discont.begin()+ki);
	  }
	else
	  ki++;
    }

  }

  //===========================================================================
  void 
    FaceSetRepair::improveVertexPos(shared_ptr<Vertex> vx, double epsge)
  //===========================================================================
  {
    // Get faces meeting in this vertex
    vector<pair<ftSurface*, Point> > faces = vx->getFaces();

	// Extract geometrical information
    vector<pair<shared_ptr<ParamSurface>, Point> > sfs;
    size_t kj;
    for (kj=0; kj<faces.size(); ++kj)
      sfs.push_back(make_pair(faces[kj].first->surface(),
			      faces[kj].second));

    // Modify vertex position. Find also improved associated
    // parameter values
    Point vertex_pos = vx->getVertexPoint();
    SurfaceTools::iterateCornerPos(vertex_pos, sfs, epsge);
    vx->setVertexPoint(vertex_pos);

    for (kj=0; kj<sfs.size(); ++kj)
      {
	shared_ptr<ParamSurface> srf = sfs[kj].first;
	shared_ptr<SplineSurface> s1 = 
	  dynamic_pointer_cast<SplineSurface, ParamSurface>(srf);
	shared_ptr<BoundedSurface> b1 = 
	  dynamic_pointer_cast<BoundedSurface, ParamSurface>(srf);
	if (b1.get())
	  {
	    bool is_spline;
	    is_spline = b1->hasUnderlyingSpline(s1);
	  }
	if (s1.get())
	  {
	    bool mod; 
	    mod = GapRemoval::modifyAtVertex(s1, sfs[kj].second, 
						  vertex_pos, epsge);
	  }
      }

  }

  //===========================================================================
  void 
  FaceSetRepair::updateOneSpline(vector<pair<ftEdge*, ftEdge*> >& pos_discont,
				 double epsge)
  //===========================================================================
  {
    // Search for B-spline surfaces where the edge is an iso parameter
    // and this is not the case for the opposite surface
    for (size_t ki=0; ki<pos_discont.size(); )
      {
	// Fetch the related faces
	ftSurface *face1 = pos_discont[ki].first->face()->asFtSurface();
	ftSurface *face2 = pos_discont[ki].second->face()->asFtSurface();

	if (!(face1 && face2))
	  continue;   // This should not happen. Cannot mend gaps

	// Check edge curves
	ftEdge* e1g = pos_discont[ki].first->geomEdge();
	ftEdge* e2g = pos_discont[ki].second->geomEdge();
	 shared_ptr<CurveOnSurface> sfcv1 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e1g->geomCurve());
	 shared_ptr<CurveOnSurface> sfcv2 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e2g->geomCurve());
	 if (!(sfcv1.get() && sfcv2.get()))
	     continue;

	 shared_ptr<BoundedSurface> bd1 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(face1->surface());
	 shared_ptr<BoundedSurface> bd2 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(face2->surface());
	 shared_ptr<SplineSurface> s1 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(face1->surface());
	 shared_ptr<SplineSurface> s2 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(face2->surface());

	 bool is_spline1 = false, is_spline2 = false;
	 if (bd1.get())
	   is_spline1 = bd1->hasUnderlyingSpline(s1);
	     
	 if (bd2.get())
	   is_spline2 = bd2->hasUnderlyingSpline(s2);

	 // Check if the edge follows the surface boundary
	 bool same1, same2;
	 int idx1 = sfcv1->whichBoundary(epsge, same1); 
	 int idx2 = sfcv2->whichBoundary(epsge, same2); 
	 shared_ptr<CurveOnSurface> dummy;

	 vector<ftEdge*> bd_edges;
	 bool modify_first;
	 if ((s1.get() && idx1 != -1) && !(s2.get() && idx2 != -1))
	   {
	     // We have found a candidate for a B-spline boundary update
	     // Fetch all edges following the surface boundary
	     getBoundaryEdges(e1g, bd_edges, epsge);
	     modify_first = true;
	   }
	 else if ((s2.get() && idx2 != -1) && !(s1.get() && idx1 != -1))
	   {
	     // We have found a candidate for a B-spline boundary update
	     // Fetch all edges following the surface boundary
	     getBoundaryEdges(e2g, bd_edges, epsge);
	     modify_first = false;
	   }

	 // Extract geometry information
	 vector<shared_ptr<CurveOnSurface> > bd_cvs1;
	 vector<shared_ptr<CurveOnSurface> > bd_cvs2;
	 vector<double> start1, end1, start2, end2;
	 size_t kj;
	 for (kj=0; kj<bd_edges.size(); ++kj)
	   {
	     shared_ptr<ParamCurve> curr1 = bd_edges[kj]->geomCurve();
	     bd_cvs1.push_back(dynamic_pointer_cast<CurveOnSurface,
			       ParamCurve>(curr1));
	     start1.push_back(bd_edges[kj]->tMin());
	     end1.push_back(bd_edges[kj]->tMax());
	     if (bd_edges[kj]->twin())
	       {
		 ftEdge *twin = bd_edges[kj]->twin()->geomEdge();
		 shared_ptr<ParamCurve> curr2 = twin->geomCurve();
		 bd_cvs2.push_back(dynamic_pointer_cast<CurveOnSurface,
				   ParamCurve>(curr2));
		 start2.push_back(twin->tMin());
		 end2.push_back(twin->tMax());
	       }
	     else
	       {
		 bd_cvs2.push_back(dummy);
		 start2.push_back(-1.0);
		 end2.push_back(-1.0);
	       }
	   }

	 bool modified = false;
	 if (bd_edges.size() > 0)
	   {
	     Point vertex1 = bd_edges[0]->getVertex(true)->getVertexPoint();
	     Point vertex2 = 
	       bd_edges[bd_edges.size()-1]->getVertex(false)->getVertexPoint();
	     if (modify_first)
	       modified = GapRemoval::removeGapSplineTrim(s1, bd_cvs1, start1, end1, 
							  bd_cvs2, start2, end2, 
							  vertex1, vertex2, epsge);
	     else
	       modified = GapRemoval::removeGapSplineTrim(s2, bd_cvs1, start1, end1, 
							  bd_cvs2, start2, end2, 
							  vertex1, vertex2, epsge);
	   }

	 if (modified)
	   {
	     for (kj=0; kj<bd_edges.size(); ++kj)
	       {
		 // Check if this edge belongs to the collection of
		 // gap edges
		 size_t kr;
		 for (kr=0; kr<pos_discont.size(); ++kr)
		   if (pos_discont[kr].first == bd_edges[kj] ||
		       pos_discont[kr].second == bd_edges[kj])
		     pos_discont.erase(pos_discont.begin()+kr);
	       }
	   }
	 else ki++;

      }
  }

  //===========================================================================
  void 
  FaceSetRepair::updateMeetingSplines(vector<pair<ftEdge*, ftEdge*> >& pos_discont,
				      double epsge)
  //===========================================================================
  {
    // Search for cases where two B-spline surfaces meet in a common
    // boundary
    // Search for B-spline surfaces where the edge is an iso parameter
    // and this is not the case for the opposite surface
    for (size_t ki=0; ki<pos_discont.size(); )
      {
	// Fetch the related faces
	ftSurface *face1 = pos_discont[ki].first->face()->asFtSurface();
	ftSurface *face2 = pos_discont[ki].second->face()->asFtSurface();

	if (!(face1 && face2))
	  continue;   // This should not happen. Cannot mend gaps

	// Check edge curves
	ftEdge* e1g = pos_discont[ki].first->geomEdge();
	ftEdge* e2g = pos_discont[ki].second->geomEdge();
	 shared_ptr<CurveOnSurface> sfcv1 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e1g->geomCurve());
	 shared_ptr<CurveOnSurface> sfcv2 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e2g->geomCurve());
	 if (!(sfcv1.get() && sfcv2.get()))
	     continue;

	 shared_ptr<BoundedSurface> bd1 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(face1->surface());
	 shared_ptr<BoundedSurface> bd2 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(face2->surface());
	 shared_ptr<SplineSurface> s1 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(face1->surface());
	 shared_ptr<SplineSurface> s2 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(face2->surface());

	 bool is_spline1 = false, is_spline2 = false;
	 if (bd1.get())
	   is_spline1 = bd1->hasUnderlyingSpline(s1);
	     
	 if (bd2.get())
	   is_spline2 = bd2->hasUnderlyingSpline(s2);
	 
	 // Check if the edge follows the surface boundary
	 bool same1, same2;
	 int idx1 = sfcv1->whichBoundary(epsge, same1); 
	 int idx2 = sfcv2->whichBoundary(epsge, same2); 

	 vector<pair<ftEdge*,ftEdge*> > bd_edges;
	 if (s1.get() && idx1 != -1 && s2.get() && idx2 != -1)
	   {
	     // Two spline surfaces meet in a common boundary curve
	     // Fetch neighbouring surfaces with the same configuration
	     // and where the composed boundary curve is smooth
	     getBoundaryEdges(e1g, e2g, bd_edges, epsge);

	     vector<shared_ptr<CurveOnSurface> > bd_cvs1;
	     vector<shared_ptr<CurveOnSurface> > bd_cvs2;
	     vector<double> start1, end1, start2, end2;
	     size_t kj;
	     for (kj=0; kj<bd_edges.size(); ++kj)
	       {
		 shared_ptr<ParamCurve> curr1 = bd_edges[kj].first->geomCurve();
		 bd_cvs1.push_back(dynamic_pointer_cast<CurveOnSurface,
				   ParamCurve>(curr1));
		 start1.push_back(bd_edges[kj].first->tMin());
		 end1.push_back(bd_edges[kj].first->tMax());
		 ftEdge *twin = bd_edges[kj].second->geomEdge();
		 shared_ptr<ParamCurve> curr2 = twin->geomCurve();
		 bd_cvs2.push_back(dynamic_pointer_cast<CurveOnSurface,
				   ParamCurve>(curr2));
		 start2.push_back(twin->tMin());
		 end2.push_back(twin->tMax());
	       }

	     vector<Point> vertex;
	     for (size_t kr=0; kr<bd_edges.size(); ++kr)
	       vertex.push_back(bd_edges[kr].first->getVertex(true)->
				getVertexPoint());
	     vertex.push_back(bd_edges[bd_edges.size()-1].first->
			      getVertex(false)->getVertexPoint());
	     GapRemoval::removeGapSpline2(bd_cvs1, start1, end1, 
					  bd_cvs2, start2, end2, 
					  vertex, epsge);
	   }

	 for (size_t kj=0; kj<bd_edges.size(); ++kj)
	   {
	     // Check if this edge belongs to the collection of
	     // gap edges
	     size_t kr;
	     for (kr=0; kr<pos_discont.size(); ++kr)
	       if (pos_discont[kr].first == bd_edges[kj].first ||
		   pos_discont[kr].second == bd_edges[kj].first ||
		   pos_discont[kr].first == bd_edges[kj].second ||
		   pos_discont[kr].second == bd_edges[kj].second)
		 pos_discont.erase(pos_discont.begin()+kr);
	   }

	 if (bd_edges.size() == 0)
	   ki++;
      }

  }

  //===========================================================================
  void 
  FaceSetRepair::averageMeetingSplines(vector<pair<ftEdge*, ftEdge*> >& pos_discont,
				      double epsge)
  //===========================================================================
  {
    // Search for cases where two B-spline surfaces meet in a common
    // boundary
    // Search for B-spline surfaces where the edge is an iso parameter
    // and this is not the case for the opposite surface
    for (size_t ki=0; ki<pos_discont.size(); )
      {
	// Fetch the related faces
	ftSurface *face1 = pos_discont[ki].first->face()->asFtSurface();
	ftSurface *face2 = pos_discont[ki].second->face()->asFtSurface();

	if (!(face1 && face2))
	  continue;   // This should not happen. Cannot mend gaps

	// Check edge curves
	ftEdge* e1g = pos_discont[ki].first->geomEdge();
	ftEdge* e2g = pos_discont[ki].second->geomEdge();
	 shared_ptr<CurveOnSurface> sfcv1 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e1g->geomCurve());
	 shared_ptr<CurveOnSurface> sfcv2 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e2g->geomCurve());
	 if (!(sfcv1.get() && sfcv2.get()))
	     continue;

	 shared_ptr<BoundedSurface> bd1 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(face1->surface());
	 shared_ptr<BoundedSurface> bd2 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(face2->surface());
	 shared_ptr<SplineSurface> s1 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(face1->surface());
	 shared_ptr<SplineSurface> s2 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(face2->surface());

	 bool is_spline1 = false, is_spline2 = false;
	 if (bd1.get())
	   is_spline1 = bd1->hasUnderlyingSpline(s1);
	     
	 if (bd2.get())
	   is_spline2 = bd2->hasUnderlyingSpline(s2);
	 
	 // Check if the edge follows the surface boundary
	 bool same1, same2;
	 int idx1 = sfcv1->whichBoundary(epsge, same1); 
	 int idx2 = sfcv2->whichBoundary(epsge, same2); 

	 if (s1.get() && idx1 != -1 && s2.get() && idx2 != -1)
	   {
	     // Two spline surfaces meet in a common boundary curve
	     // Average coefficients along common boundary
	     double start1, end1, start2, end2;
	     start1 = e1g->tMin();
	     end1 = e1g->tMax();
	     start2 = e2g->tMin();
	     end2 = e2g->tMax();
	     Point vertex1 = e1g->getVertex(true)->getVertexPoint();
	     Point vertex2 = e1g->getVertex(false)->getVertexPoint();
	     GapRemoval::removeGapSpline(s1, sfcv1, start1, end1, 
					 s2, sfcv2, start2, end2, 
					 vertex1, vertex2, epsge);

	     // Remove this edge belongs to the collection of gap edges
	     pos_discont.erase(pos_discont.begin()+ki);
	   }
	 else
	   ki++;
      }

  }

    //===========================================================================
  void 
  FaceSetRepair::modifyOneSpline(vector<pair<ftEdge*, ftEdge*> >& pos_discont,
				 double epsge)
  //===========================================================================
  {
    // Search for configurations with one B-spline surface and one not
    for (size_t ki=0; ki<pos_discont.size(); )
      {
	// Fetch the related faces
	ftSurface *face1 = pos_discont[ki].first->face()->asFtSurface();
	ftSurface *face2 = pos_discont[ki].second->face()->asFtSurface();

	if (!(face1 && face2))
	  continue;   // This should not happen. Cannot mend gaps

	// Check edge curves
	ftEdge* e1g = pos_discont[ki].first->geomEdge();
	ftEdge* e2g = pos_discont[ki].second->geomEdge();
	 shared_ptr<CurveOnSurface> sfcv1 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e1g->geomCurve());
	 shared_ptr<CurveOnSurface> sfcv2 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e2g->geomCurve());
	 if (!(sfcv1.get() && sfcv2.get()))
	     continue;

	 shared_ptr<ParamSurface> srf1 = face1->surface();
	 shared_ptr<ParamSurface> srf2 = face2->surface();
	 shared_ptr<BoundedSurface> bd1 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(srf1);
	 shared_ptr<BoundedSurface> bd2 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(srf2);
	 shared_ptr<SplineSurface> s1 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(srf1);
	 shared_ptr<SplineSurface> s2 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(srf2);

	 bool is_spline1 = false, is_spline2 = false;
	 if (bd1.get())
	   is_spline1 = bd1->hasUnderlyingSpline(s1);
	     
	 if (bd2.get())
	   is_spline2 = bd2->hasUnderlyingSpline(s2);

	 vector<ftEdge*> bd_edges;
	 bool modify_first;
	 if (s1.get() && !s2.get())
	   {
	     // We have found a candidate for a B-spline update
	     // Fetch all edges in this surface
	     getEdges(e1g, bd_edges, epsge);
	     modify_first = true;
	   }
	 else if (s2.get() && !s1.get())
	   {
	     // We have found a candidate for a B-splineupdate
	     // Fetch all edges in this surface
	     getEdges(e2g, bd_edges, epsge);
	     modify_first = false;
	   }

	 // Extract geometry information
	 vector<shared_ptr<CurveOnSurface> > bd_cvs1;
	 vector<shared_ptr<CurveOnSurface> > bd_cvs2;
	 vector<double> start1, end1, start2, end2;
	 shared_ptr<CurveOnSurface> dummy;
	 size_t kj;
	 for (kj=0; kj<bd_edges.size(); ++kj)
	   {
	     shared_ptr<ParamCurve> curr1 = bd_edges[kj]->geomCurve();
	     bd_cvs1.push_back(dynamic_pointer_cast<CurveOnSurface,
			       ParamCurve>(curr1));
	     start1.push_back(bd_edges[kj]->tMin());
	     end1.push_back(bd_edges[kj]->tMax());
	     if (bd_edges[kj]->twin())
	       {
		 ftEdge *twin = bd_edges[kj]->twin()->geomEdge();
		 shared_ptr<ParamCurve> curr2 = twin->geomCurve();
		 bd_cvs2.push_back(dynamic_pointer_cast<CurveOnSurface,
				   ParamCurve>(curr2));
		 start2.push_back(twin->tMin());
		 end2.push_back(twin->tMax());
	       }
	     else
	       {
		 bd_cvs2.push_back(dummy);
		 start2.push_back(-1.0);
		 end2.push_back(-1.0);
	       }
	   }

	 if (bd_edges.size() > 0)
	   {
	     if (modify_first)
	       GapRemoval::modifySplineSf(srf1, bd_cvs1, start1, end1, 
					  srf2, bd_cvs2, 
					  start2, end2, epsge);
	     else
	       GapRemoval::modifySplineSf(srf2, bd_cvs1, start1, end1, 
					  srf1, bd_cvs2, 
					  start2, end2, epsge);
	   }
	 else ki++;

	 for (kj=0; kj<bd_edges.size(); ++kj)
	   {
	     // Check if this edge belongs to the collection of
	     // gap edges
	     size_t kr;
	     for (kr=0; kr<pos_discont.size(); ++kr)
	       if (pos_discont[kr].first == bd_edges[kj] ||
		   pos_discont[kr].second == bd_edges[kj])
		 pos_discont.erase(pos_discont.begin()+kr);
	   }

      }
  }

  //===========================================================================
  void 
  FaceSetRepair::modifyTwoSplines(vector<pair<ftEdge*, ftEdge*> >& pos_discont,
				  double epsge)
  //===========================================================================
  {
    // Search for configurations with two B-spline surfaces 
    for (size_t ki=0; ki<pos_discont.size(); )
      {
	// Fetch the related faces
	ftSurface *face1 = pos_discont[ki].first->face()->asFtSurface();
	ftSurface *face2 = pos_discont[ki].second->face()->asFtSurface();

	if (!(face1 && face2))
	  continue;   // This should not happen. Cannot mend gaps

	// Check edge curves
	ftEdge* e1g = pos_discont[ki].first->geomEdge();
	ftEdge* e2g = pos_discont[ki].second->geomEdge();
	 shared_ptr<CurveOnSurface> sfcv1 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e1g->geomCurve());
	 shared_ptr<CurveOnSurface> sfcv2 =
	   dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	   (e2g->geomCurve());
	 if (!(sfcv1.get() && sfcv2.get()))
	     continue;

	 shared_ptr<ParamSurface> psurf1 = face1->surface();
	 shared_ptr<ParamSurface> psurf2 = face2->surface();
	 shared_ptr<BoundedSurface> bd1 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(psurf1);
	 shared_ptr<BoundedSurface> bd2 = 
	   dynamic_pointer_cast<BoundedSurface, ParamSurface>(psurf2);
	 shared_ptr<SplineSurface> s1 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(psurf1);
	 shared_ptr<SplineSurface> s2 =
	   dynamic_pointer_cast<SplineSurface, ParamSurface>(psurf2);

	 bool is_spline1 = false, is_spline2 = false;
	 if (bd1.get())
	   is_spline1 = bd1->hasUnderlyingSpline(s1);
	     
	 if (bd2.get())
	   is_spline2 = bd2->hasUnderlyingSpline(s2);

	 vector<pair<ftEdge*,ftEdge*> > bd_edges;
	 if (s1.get() && s2.get())
	   {
	     // We have found a candidate for a B-spline update
	     // Fetch all edges in this surface
	     getEdges(e1g,e2g, bd_edges, epsge);
	   }

	 // Extract geometry information
	 vector<shared_ptr<CurveOnSurface> > bd_cvs1;
	 vector<shared_ptr<CurveOnSurface> > bd_cvs2;
	 vector<double> start1;
	 vector<double> end1;
	 vector<double> start2; 
	 vector<double> end2;
	 size_t kj;
	 vector<Point> vertex;
	 for (kj=0; kj<bd_edges.size(); ++kj)
	   {
	     shared_ptr<ParamCurve> curr1 = bd_edges[kj].first->geomCurve();
	     bd_cvs1.push_back(dynamic_pointer_cast<CurveOnSurface,
			       ParamCurve>(curr1));
	     start1.push_back(bd_edges[kj].first->tMin());
	     end1.push_back(bd_edges[kj].first->tMax());
	     shared_ptr<ParamCurve> curr2 = bd_edges[kj].second->geomCurve();
	     bd_cvs2.push_back(dynamic_pointer_cast<CurveOnSurface,
			       ParamCurve>(curr2));
	     start2.push_back(bd_edges[kj].second->tMin());
	     end2.push_back(bd_edges[kj].second->tMax());
	     vertex.push_back(bd_edges[kj].first->getVertex(true)->
			      getVertexPoint());
	   }

	 if (bd_edges.size() > 0)
	   {
	     vertex.push_back(bd_edges[bd_edges.size()-1].first->
			      getVertex(false)->getVertexPoint());
	     GapRemoval::modifySplines(psurf1, bd_cvs1, start1, end1, 
				       psurf2, bd_cvs2, start2, end2, 
				       vertex, epsge);
	   }
	 else ki++;

	 for (kj=0; kj<bd_edges.size(); ++kj)
	   {
	     // Check if this edge belongs to the collection of
	     // gap edges
	     size_t kr;
	     for (kr=0; kr<pos_discont.size(); ++kr)
	       if (pos_discont[kr].first == bd_edges[kj].first ||
		   pos_discont[kr].second == bd_edges[kj].first ||
		   pos_discont[kr].first == bd_edges[kj].second ||
		   pos_discont[kr].second == bd_edges[kj].second)
		 pos_discont.erase(pos_discont.begin()+kr);
	   }

      }
  }

  //===========================================================================
  void 
  FaceSetRepair::getBoundaryEdges(ftEdge *edge1, vector<ftEdge*>& edges,
				  double epsge)
  //===========================================================================
  {
    edges.clear();

    ftEdge *e1 = edge1;

    shared_ptr<CurveOnSurface> sfcv =
      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());

    if (!sfcv.get())
      return;

    bool same;
    int idx = sfcv->whichBoundary(epsge, same);
    if (idx == -1)
      return;

    shared_ptr<ParamSurface> sf = sfcv->underlyingSurface();

    edges.push_back(e1);
    while (true)
      {
	e1 = e1->next()->geomEdge();
	if (e1 == edges[0])
	  break;

	sfcv = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
	shared_ptr<ParamSurface> sf2 = sfcv->underlyingSurface();
	if (sf.get() != sf2.get())
	  break;

	int idx2 = sfcv->whichBoundary(epsge, same);
	if (idx2 != idx)
	  break;
	edges.push_back(e1);

      }

    e1 = edge1;
    while (true)
      {
	e1 = e1->prev()->geomEdge();
	if (e1 == edges[edges.size()-1])
	  break;

	sfcv = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
	shared_ptr<ParamSurface> sf2 = sfcv->underlyingSurface();
	if (sf.get() != sf2.get())
	  break;

	int idx2 = sfcv->whichBoundary(epsge, same);
	if (idx2 != idx)
	  break;
	edges.insert(edges.begin(), e1);

      }

    
    
  }

  //===========================================================================
  void 
  FaceSetRepair::getEdges(ftEdge *edge1, vector<ftEdge*>& edges, double epsge)
  //===========================================================================
  {
    edges.clear();

    ftEdge *e1 = edge1;

    shared_ptr<CurveOnSurface> sfcv =
      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());

    if (!sfcv.get())
      return;

    shared_ptr<ParamSurface> sf = sfcv->underlyingSurface();

    edges.push_back(e1);
    while (true)
      {
	e1 = e1->next()->geomEdge();
	if (e1 == edges[0])
	  break;

	if (!e1->twin())
	  continue;

	sfcv = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
	shared_ptr<ParamSurface> sf2 = sfcv->underlyingSurface();
	if (sf.get() != sf2.get())
	  break;

	edges.push_back(e1);

      }

    e1 = edge1;
    while (true)
      {
	e1 = e1->prev()->geomEdge();
	if (e1 == edges[edges.size()-1])
	  break;

	if (!e1->twin())
	  continue;

	sfcv = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
	shared_ptr<ParamSurface> sf2 = sfcv->underlyingSurface();
	if (sf.get() != sf2.get())
	  break;

	edges.insert(edges.begin(), e1);

      }

  }

  //===========================================================================
  void 
  FaceSetRepair::getEdges(ftEdge *edge1, ftEdge *edge2,
			  vector<pair<ftEdge*, ftEdge*> >& edges, double epsge)
  //===========================================================================
  {
    edges.clear();

    ftEdge *e1 = edge1;
    ftEdge *e2 = edge2;

    shared_ptr<CurveOnSurface> sfcv1 =
      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
    shared_ptr<CurveOnSurface> sfcv2 =
      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e2->geomCurve());

    if (!sfcv1.get() || !sfcv2.get())
      return;

    shared_ptr<ParamSurface> sf1 = sfcv1->underlyingSurface();
    shared_ptr<ParamSurface> sf2 = sfcv2->underlyingSurface();

    edges.push_back(make_pair(e1, e2));
    while (true)
      {
	e1 = e1->next()->geomEdge();
	if (e1 == edges[0].first)
	  break;

	sfcv1 = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
	shared_ptr<ParamSurface> sf1_2 = sfcv1->underlyingSurface();
	if (sf1.get() != sf1_2.get())
	  break;

	if (!e1->twin())
	  break;
	e2 = e1->twin()->geomEdge();
	if (!e2)
	  break;

	sfcv2 = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e2->geomCurve());
	shared_ptr<ParamSurface> sf2_2 = sfcv2->underlyingSurface();
	if (sf2.get() != sf2_2.get())
	  break;

	edges.push_back(make_pair(e1, e2));

      }

    e1 = edge1;
    while (true)
      {
	e1 = e1->prev()->geomEdge();
	if (e1 == edges[edges.size()-1].first)
	  break;

	sfcv1 = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
	shared_ptr<ParamSurface> sf1_2 = sfcv1->underlyingSurface();
	if (sf1.get() != sf1_2.get())
	  break;

	if (!e1->twin())
	  break;
	e2 = e1->twin()->geomEdge();
	if (!e2)
	  break;

	sfcv2 = dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e2->geomCurve());
	shared_ptr<ParamSurface> sf2_2 = sfcv2->underlyingSurface();
	if (sf2.get() != sf2_2.get())
	  break;

	edges.insert(edges.begin(), make_pair(e1, e2));

      }

  }

  //===========================================================================
  void 
  FaceSetRepair::getBoundaryEdges(ftEdge *edge1, ftEdge *edge2, 
				  vector<pair<ftEdge*,ftEdge*> >& edges,
				  double epsge)
  //===========================================================================
  {
    edges.clear();

    ftEdge* e1 = edge1;
    ftEdge* e2 = edge2;

    shared_ptr<CurveOnSurface> sfcv1 =
      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());

    shared_ptr<CurveOnSurface> sfcv2 =
      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e2->geomCurve());

    if (!(sfcv1.get() && sfcv2.get()))
      return;

    bool same1;
    int idx1 = sfcv1->whichBoundary(epsge, same1);
    if (idx1 == -1)
      return;

    bool same2;
    int idx2 = sfcv2->whichBoundary(epsge, same2);
    if (idx2 == -1)
      return;

    edges.push_back(make_pair(e1, e2));
    while (true)
      {
	ftEdge *e3 = e1->next()->geomEdge();
	if (e3 == edges[0].first)
	  break;

	ftEdge *e4 = e2->prev()->geomEdge();
	if (e4 == edges[0].second)
	  break;

	// Check continuity
	tpJointType joint1 = e1->checkContinuity(e3, 
						sfmodel_->getTolerances().neighbour,
						sfmodel_->getTolerances().gap,
						sfmodel_->getTolerances().bend,
						sfmodel_->getTolerances().kink);
	tpJointType joint2 = e4->checkContinuity(e2, 
						sfmodel_->getTolerances().neighbour,
						sfmodel_->getTolerances().gap,
						sfmodel_->getTolerances().bend,
						sfmodel_->getTolerances().kink);
	if (joint1 >= JOINT_KINK && joint2 >= JOINT_KINK)
	  break;
	
	if (joint1 < JOINT_KINK)
	  {
	    e1 = e3;
	    if (!e1->twin())
	      break;
	    e2 = e1->twin()->geomEdge();
	    if (!e2)
	      break;
	  }
	else
	  {
	    e2 = e4;
	    if (!e2->twin())
	      break;
	    e1 = e2->twin()->geomEdge();
	    if (!e1)
	      break;
	  }
	
	shared_ptr<CurveOnSurface> sfcv =
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
	if (!sfcv.get())
	  break;

	shared_ptr<ParamSurface> sf = sfcv->underlyingSurface();
	shared_ptr<BoundedSurface> bd1 = 
	  dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
	shared_ptr<SplineSurface> s1 =
	  dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);
	if (!(s1.get() || (bd1.get() &&  bd1->hasUnderlyingSpline(s1))))
	  break;  // Not a spline

	bool same;
	int idx = sfcv->whichBoundary(epsge, same);
	if (idx == -1)
	  break;  // Not a boundary trim

	sfcv =
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e2->geomCurve());
	if (!sfcv.get())
	  break;

	sf = sfcv->underlyingSurface();
	bd1 = dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
	s1 = dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);
	if (!(s1.get() || (bd1.get() &&  bd1->hasUnderlyingSpline(s1))))
	  break;  // Not a spline

	idx = sfcv->whichBoundary(epsge, same);
	if (idx == -1)
	  break; // Not a boundary trim

	// A continuous continuation of the surface boundary is found
	edges.push_back(make_pair(e1,e2));
      }

    e1 = edge1;
    e2 = edge2;
    while (true)
      {
	ftEdge *e3 = edge1->prev()->geomEdge();
	if (e3 == edges[edges.size()-1].first)
	  break;

	ftEdge *e4 = e2->next()->geomEdge();
	if (e4 == edges[edges.size()-1].second)
	  break;

	// Check continuity
	tpJointType joint1 = e3->checkContinuity(e1, 
						 sfmodel_->getTolerances().neighbour,
						 sfmodel_->getTolerances().gap,
						 sfmodel_->getTolerances().bend,
						 sfmodel_->getTolerances().kink);
	tpJointType joint2 = e2->checkContinuity(e4, 
						sfmodel_->getTolerances().neighbour,
						sfmodel_->getTolerances().gap,
						sfmodel_->getTolerances().bend,
						sfmodel_->getTolerances().kink);
	if (joint1 >= JOINT_KINK && joint2 >= JOINT_KINK)
	  break;
	
	if (joint1 < JOINT_KINK)
	  {
	    e1 = e3;
	    if (!e1->twin())
	      break;
	    e2 = e1->twin()->geomEdge();
	    if (!e2)
	      break;
	  }
	else
	  {
	    e2 = e4;
	    if (!e2->twin())
	      break;
	    e1 = e2->twin()->geomEdge();
	    if (!e1)
	      break;
	  }
	
	shared_ptr<CurveOnSurface> sfcv =
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e1->geomCurve());
	if (!sfcv.get())
	  break;

	shared_ptr<ParamSurface> sf = sfcv->underlyingSurface();
	shared_ptr<BoundedSurface> bd1 = 
	  dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
	shared_ptr<SplineSurface> s1 =
	  dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);
	if (!(s1.get() || (bd1.get() &&  bd1->hasUnderlyingSpline(s1))))
	  break;  // Not a spline

	bool same;
	int idx = sfcv->whichBoundary(epsge, same);
	if (idx == -1)
	  break; // Not a boundary trim

	sfcv =
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(e2->geomCurve());
	if (!sfcv.get())
	  break;

	sf = sfcv->underlyingSurface();
	bd1 = dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
	s1 = dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);
	if (!(s1.get() || (bd1.get() &&  bd1->hasUnderlyingSpline(s1))))
	  break;  // Not a spline

	idx = sfcv->whichBoundary(epsge, same);
	if (idx == -1)
	  break; // Not a boundary trim


	edges.insert(edges.begin(), make_pair(e1,e2));

      }
  }

}  // namespace Go
