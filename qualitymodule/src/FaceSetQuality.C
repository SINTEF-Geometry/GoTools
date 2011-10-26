//===========================================================================
//                                                                           
// File: FaceSetQuality.C
//                                                                           
// Created: April 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/intersections/Singular.h"
#include "GoTools/qualitymodule/QualityUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/CurvatureAnalysis.h"
#include "GoTools/geometry/Curvature.h"
#include "GoTools/geometry/PointOnCurve.h"
#include <fstream>

using std::set;
using std::make_pair;
using std::shared_ptr;
using std::dynamic_pointer_cast;

namespace Go
{
    using namespace qualityUtils;

  using namespace qualityUtils;

  //===========================================================================
  FaceSetQuality::FaceSetQuality(double gap,   // Gap between adjacent surfaces
			     double kink,  // Kink between adjacent surfaces 
			     double approx)
  //===========================================================================
      : ModelQuality(gap, kink, approx)
  {
  }


  //===========================================================================
  FaceSetQuality::FaceSetQuality(const tpTolerances& toptol, 
				 double approx)
  //===========================================================================
      : ModelQuality(toptol, approx)
  {
  }


  //===========================================================================
  FaceSetQuality::FaceSetQuality(shared_ptr<SurfaceModel> sfmodel)
  //===========================================================================
    : ModelQuality(sfmodel->getTolerances(), sfmodel->getApproximationTol())
  {
      model_ = sfmodel;
   }


  //===========================================================================
  FaceSetQuality::~FaceSetQuality()
  //===========================================================================
  {
  }

  //===========================================================================
  void FaceSetQuality::attach(shared_ptr<SurfaceModel> sfmodel)
  //===========================================================================
  {
      model_ = sfmodel;
  }

  //===========================================================================
  void FaceSetQuality::degenSurfaces(vector<shared_ptr<ParamSurface> >& deg_sfs)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(DEGEN_SRF_BD, tol) && tol == toptol_.neighbour)
      {
	  vector<shared_ptr<ftSurface> > deg_faces;
	  deg_faces = results_->getDegSfs();
	  deg_sfs.resize(deg_faces.size());
	  for (size_t ki=0; ki<deg_faces.size(); ++ki)
	      deg_sfs[ki] = deg_faces[ki]->surface();
	  return;
      }

      deg_sfs.clear();

      results_->reset(DEGEN_SRF_BD);
      results_->performtest(DEGEN_SRF_BD, toptol_.neighbour);

      int nmb_sfs = model_->nmbEntities();
      for (int ki=0; ki<nmb_sfs; ki++)
      {
	  bool is_degen = model_->isDegenerate(ki);
	  if (is_degen)
	  {
	      // Store in result container
	      results_->addDegSf(model_->getFace(ki));

	      // Return surface
	      deg_sfs.push_back(model_->getSurface(ki));
	  }
      }
  }



  //===========================================================================
  void FaceSetQuality::degenerateSfCorners(vector<shared_ptr<ftPoint> >& deg_corners)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(DEGEN_SRF_CORNER, tol) && tol == toptol_.kink)
      {
	  deg_corners = results_->getDegCorners();
	  return;
      }

      deg_corners.clear();

      results_->reset(DEGEN_SRF_CORNER);
      results_->performtest(DEGEN_SRF_CORNER, toptol_.kink);

      int nmb_sfs = model_->nmbEntities();
      for (int ki=0; ki<nmb_sfs; ++ki)
      {
	  shared_ptr<ftSurface> face = model_->getFace(ki);
	  shared_ptr<ParamSurface> surf = face->surface();
	  vector<Point> tmp_corners;
	  surf->getDegenerateCorners(tmp_corners, toptol_.kink);

	  for (size_t kj=0; kj<tmp_corners.size(); ++kj)
	  {
	      Point pnt = surf->point(tmp_corners[kj][0], tmp_corners[kj][1]);
	      shared_ptr<ftPoint> curr = shared_ptr<ftPoint>(new ftPoint(pnt, face.get(), 
									 tmp_corners[kj][0], 
									 tmp_corners[kj][1]));
	      deg_corners.push_back(curr);
	      results_->addDegenerateSfCorner(curr);
	  }
      }
  }

  //===========================================================================
  void FaceSetQuality::identicalVertices(vector<pair<shared_ptr<Vertex>,
					 shared_ptr<Vertex> > >& identical_vertices)

  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(IDENTICAL_VERTICES, tol) && tol == toptol_.neighbour)
      {
	  identical_vertices = results_->getIdenticalVertices();
	  return;
      }

      identical_vertices.clear();
      results_->reset(IDENTICAL_VERTICES);
      results_->performtest(IDENTICAL_VERTICES, toptol_.neighbour);
      int nmb_sfs = model_->nmbEntities();
      int ki;
      std::set<shared_ptr<Vertex> > all_vertices;  // All vertices in the model represented once

      // Collect all vertices
      for (ki=0;  ki<nmb_sfs; ki++)
      {
	  vector<shared_ptr<Vertex> > curr_vertices = model_->getFace(ki)->vertices();
	  all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
      }

      // Check distance between all pairs of vertices
      set<shared_ptr<Vertex> >::iterator iter1 = all_vertices.begin();
      set<shared_ptr<Vertex> >::iterator iter2;
      set<shared_ptr<Vertex> >::iterator last = all_vertices.end();
      for (; iter1 != last; ++iter1)
	  for (iter2=iter1, ++iter2; iter2 != last; ++iter2)
	  {
	      double dist = (*iter1)->getVertexPoint().dist((*iter2)->getVertexPoint());
	      if (dist < toptol_.neighbour)
	      {
		  pair<shared_ptr<Vertex>, shared_ptr<Vertex> > identical =
		      make_pair(*iter1, *iter2);
		  identical_vertices.push_back(identical);
		  results_->addIdenticalVertices(identical);
	      }
	  }
		  
  }


  //===========================================================================
  void FaceSetQuality::identicalOrEmbeddedEdges(vector<pair<shared_ptr<ftEdge>,
						shared_ptr<ftEdge> > >& identical_edges,
						vector<pair<shared_ptr<ftEdge>,
						shared_ptr<ftEdge> > >& embedded_edges)

  //===========================================================================
  {
      // Check if the SurfaceModel is undefined. @jbt
      if (model_->nmbEntities() == 0) {
	  MESSAGE("SurfaceModel is undefined - returning empty results");
	  identical_edges.clear();
	  embedded_edges.clear();
	  return;
      }

      // Check if the test is performed already
      double tol1, tol2;
      if (results_->testPerformed(IDENTICAL_EDGES, tol1) && tol1 == toptol_.neighbour &&
	  results_->testPerformed(EMBEDDED_EDGES, tol2) && tol2 == toptol_.neighbour)
      {
	  identical_edges = results_->getIdenticalEdges();
	  embedded_edges = results_->getEmbeddedEdges();
	  return;
      }

      identical_edges.clear();
      embedded_edges.clear();

      results_->reset(IDENTICAL_EDGES);
      results_->performtest(IDENTICAL_EDGES, toptol_.neighbour);

      results_->reset(EMBEDDED_EDGES);
      results_->performtest(EMBEDDED_EDGES, toptol_.neighbour);

      
      std::ofstream file("id_edge_out.txt");   // Debug output
//      int nmb_sfs = model_->nmbEntities();
//      int ki;
      Identity ident;
      int coincidence;

//       // Collect all edges
//       std::set<shared_ptr<ftEdgeBase> > all_edges;  // All edges in the model represented once
//       for (ki=0;  ki<nmb_sfs; ki++)
//       {
// 	  // The function returns existing edges if there are any
// 	  vector<shared_ptr<ftEdgeBase> > curr_edges = 
// 	      model_->getFace(ki)->createInitialEdges();
// 	  all_edges.insert(curr_edges.begin(), curr_edges.end());
//       }

//        // Check coincidence between all pairs of edges (not twins)
//       set<shared_ptr<ftEdgeBase> >::iterator iter1 = all_edges.begin();
//       set<shared_ptr<ftEdgeBase> >::iterator iter2;
//       set<shared_ptr<ftEdgeBase> >::iterator last = all_edges.end();
//       for (; iter1 != last; ++iter1)
// 	  for (iter2=iter1, ++iter2; iter2 != last; ++iter2)
// 	  {
// 	      if ((*iter1)->twin() == (*iter2).get())
// 		  continue;

      // Get candidate edges
      vector<pair<shared_ptr<ftEdgeBase>, shared_ptr<ftEdgeBase> > > candidates;
      model_->getOverlappingEdges(toptol_.neighbour, candidates);
      for (size_t kj=0; kj<candidates.size(); ++kj)
      {
	  coincidence = ident.identicalCvs(candidates[kj].first->geomEdge()->geomCurve(), 
					   candidates[kj].first->tMin(), 
					   candidates[kj].first->tMax(), 
					   candidates[kj].second->geomEdge()->geomCurve(),
					   candidates[kj].second->tMin(), 
					   candidates[kj].second->tMax(), 
					   toptol_.neighbour);
	  if (coincidence > 0)
	  {
	      file << kj << "edge1: " << candidates[kj].first.get() << ", twin: ";
	      file << candidates[kj].first.get()->twin() << ", face: ";
	      file << candidates[kj].first.get()->face();
	      file << ", edge2: " << candidates[kj].second.get() << ", twin: ";
	      file << candidates[kj].second.get()->twin() <<  ", face: ";
	      file << candidates[kj].second.get()->face() << std::endl;
	  }

	  if (coincidence == 1)
	  {
	      pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > coinc_crvs = 
		  make_pair(dynamic_pointer_cast<ftEdge, ftEdgeBase>(candidates[kj].first),
			    dynamic_pointer_cast<ftEdge, ftEdgeBase>(candidates[kj].second));
	      identical_edges.push_back(coinc_crvs);
	      results_->addIdenticalEdges(coinc_crvs);
	  }
	  else if (coincidence == 2)
	  {
	      pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > embedded_crvs = 
		  make_pair(dynamic_pointer_cast<ftEdge, ftEdgeBase>(candidates[kj].second),
			    dynamic_pointer_cast<ftEdge, ftEdgeBase>(candidates[kj].first));
	      embedded_edges.push_back(embedded_crvs);
	      results_->addEmbeddedEdges(embedded_crvs);
	  }
	  else if (coincidence == 3)
	  {
	      pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > embedded_crvs = 
		  make_pair(dynamic_pointer_cast<ftEdge, ftEdgeBase>(candidates[kj].first),
			    dynamic_pointer_cast<ftEdge, ftEdgeBase>(candidates[kj].second));
	      embedded_edges.push_back(embedded_crvs);
	      results_->addEmbeddedEdges(embedded_crvs);
	  }
      }
		  
 }

  //===========================================================================
  void FaceSetQuality::identicalOrEmbeddedFaces(std::vector<std::pair<std::shared_ptr<ftSurface>,
					      std::shared_ptr<ftSurface> > >& identical_faces,
					      std::vector<std::pair<std::shared_ptr<ftSurface>,
					      std::shared_ptr<ftSurface> > >& embedded_faces)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol1, tol2;
      if (results_->testPerformed(IDENTICAL_FACES, tol1) && tol1 == toptol_.neighbour &&
	  results_->testPerformed(EMBEDDED_FACES, tol2) && tol2 == toptol_.neighbour)
      {
	  identical_faces = results_->getIdenticalFaces();
	  embedded_faces = results_->getEmbeddedFaces();
	  return;
      }

      identical_faces.clear();
      embedded_faces.clear();

      results_->reset(IDENTICAL_FACES);
      results_->performtest(IDENTICAL_FACES, toptol_.neighbour);

      results_->reset(EMBEDDED_FACES);
      results_->performtest(EMBEDDED_FACES, toptol_.neighbour);

      Identity ident;
      int coincidence;
//      int nmb_sfs = model_->nmbEntities();
//       int ki, kj;
//       for (ki=0;  ki<nmb_sfs; ki++)
//       {
// 	  shared_ptr<ParamSurface> surf1 = model_->getSurface(ki);
// 	  for (kj=ki+1; kj<nmb_sfs; kj++)
// 	  {
// 	      shared_ptr<ParamSurface> surf2 = model_->getSurface(kj);
      vector<pair<ftSurface*, ftSurface*> > candidates;
      model_->getOverlappingFaces(toptol_.neighbour, candidates);
      for (size_t kj=0; kj<candidates.size(); ++kj)
      {
	  shared_ptr<ParamSurface> surf1 = candidates[kj].first->surface();
	  shared_ptr<ParamSurface> surf2 = candidates[kj].second->surface();
	  coincidence = ident.identicalSfs(surf1, surf2, toptol_.neighbour);

	  if (coincidence > 0)
	  {
	      int idx1 = model_->getIndex(candidates[kj].first);
	      int idx2 = model_->getIndex(candidates[kj].second);
	      pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > hit = 
		  make_pair(model_->getFace(idx1), model_->getFace(idx2));
	      
	      if (coincidence == 1)
	      {
		  identical_faces.push_back(hit);
		  results_->addIdenticalFaces(hit);
	      }
	      else if (coincidence == 2)
	      {
		  embedded_faces.push_back(hit);
		  results_->addEmbeddedFaces(hit);
	      }
	      else if (coincidence == 3)
	      {
		  embedded_faces.push_back(hit);
		  results_->addEmbeddedFaces(hit);
	      }
	  }
      }
  }    


//===========================================================================
  void FaceSetQuality::miniEdges(vector<shared_ptr<ftEdge> >& mini_edges)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(MINI_EDGE, tol) && tol == small_size_)
      {
	  mini_edges = results_->getMiniEdges();
	  return;
      }

      mini_edges.clear();

      results_->reset(MINI_EDGE);
      results_->performtest(MINI_EDGE, small_size_);

      int nmb_sfs = model_->nmbEntities();
      int ki;
      // Collect all edges
      std::set<shared_ptr<ftEdgeBase> > all_edges;  // All edges in the model represented once
      for (ki=0;  ki<nmb_sfs; ki++)
      {
	  // The function returns existing edges if there are any
	  vector<shared_ptr<ftEdgeBase> > curr_edges = 
	      model_->getFace(ki)->createInitialEdges();
	  all_edges.insert(curr_edges.begin(), curr_edges.end());
      }

      // Check edge length
      set<shared_ptr<ftEdgeBase> >::iterator iter = all_edges.begin();
      set<shared_ptr<ftEdgeBase> >::iterator last = all_edges.end();
      for (; iter != last; ++iter)
      {
	  // Quick check on three points
	  shared_ptr<ParamCurve> crv = (*iter)->geomEdge()->geomCurve();
	  Point pt1, pt2, pt3;
	  double parmin = (*iter)->tMin();
	  double parmax = (*iter)->tMax();
	  pt1 = crv->point(parmin);
	  pt2 = crv->point(0.5*(parmin+parmax));
	  pt3 = crv->point(parmax);
	  double len = pt1.dist(pt2) + pt2.dist(pt3);
	  if (len >= small_size_)
	      continue;

	  // A more exact length computation
	  len = crv->length(toptol_.gap, parmin, parmax);
	  if (len < small_size_)
	  {
	      mini_edges.push_back(dynamic_pointer_cast<ftEdge, ftEdgeBase>(*iter));
	      results_->addMiniEdge(dynamic_pointer_cast<ftEdge, ftEdgeBase>(*iter));
	  }
      }
      
      return;
  }

  //===========================================================================
  void FaceSetQuality::miniSurfaces(vector<shared_ptr<ftSurface> >& mini_surfaces)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(MINI_FACE, tol) && tol == small_size_*small_size_)
      {
	  mini_surfaces = results_->getMiniFaces();
	  return;
      }

      mini_surfaces.clear();

      results_->reset(MINI_SURFACE);
      results_->reset(MINI_FACE);
      results_->performtest(MINI_SURFACE, small_size_*small_size_);
      results_->performtest(MINI_FACE, small_size_*small_size_);

      // Check surface area
      int nmb_sfs = model_->nmbEntities();
      int ki;
      double size_fac = 10.0;
      double small_size2 = small_size_*small_size_;
      for (ki=0; ki<nmb_sfs; ++ki)
      {
	  // The function needs that bounded surfaces have correctly oriented
	  // boundary loops. Make sure that this is the case
	  bool turned;
	  turned = model_->getFace(ki)->checkAndFixBoundaries();

	  // A pre check to find out if a proper area calculation is needed
	  shared_ptr<ParamSurface> surf = model_->getSurface(ki);
	  double area_estimate = qualityUtils::estimateArea(surf);
	  if (area_estimate > size_fac*small_size2)
	    continue;

	  // Compute area
	  double area = model_->getFace(ki)->area(toptol_.neighbour);
	  if (area < small_size2)
	  {
	      mini_surfaces.push_back(model_->getFace(ki));
	      results_->addMiniSurface(model_->getSurface(ki));
	      results_->addMiniFace(model_->getFace(ki));
	  }
      }
  }

  //===========================================================================
  void FaceSetQuality::vanishingSurfaceNormal(vector<shared_ptr<ftPoint> >& singular_points,
					      vector<shared_ptr<ftCurve> >& singular_curves)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(VANISHING_NORMAL, tol) && tol == toptol_.gap)
      {
	  singular_points = results_->getSingPnts();
	  singular_curves = results_->getSingCrvs();
	  return;
      }

      singular_points.clear();
      singular_curves.clear();

      results_->reset(VANISHING_NORMAL);
      results_->performtest(VANISHING_NORMAL, toptol_.gap);

      int nmb_sfs = model_->nmbEntities();
      for (int ki=0; ki<nmb_sfs; ki++)
      {
	  shared_ptr<ftSurface> face = model_->getFace(ki);
	  shared_ptr<ParamSurface> surf = model_->getSurface(ki);
	  vector<Point> singular_pts;
	  vector<vector<Point> > singular_sequences;

	  Singular::vanishingNormal(surf, toptol_.gap, singular_pts, singular_sequences);

	  size_t kj, kr;
	  for (kj=0; kj<singular_pts.size(); kj++)
	  {
	      double u = singular_pts[kj][0];
	      double v = singular_pts[kj][1];
	      Point pos = surf->point(u, v);
	      shared_ptr<ftPoint> sing_pt = shared_ptr<ftPoint>(new ftPoint(pos, face.get(), u, v));
	      singular_points.push_back(sing_pt);

	      results_->addSingPnt(sing_pt);
	  }

	  for (kj=0; kj<singular_sequences.size(); kj++)
	  {
	      shared_ptr<ftCurve> curr_crv = shared_ptr<ftCurve>(new ftCurve(CURVE_SINGULAR));
	      for (kr=1; kr<singular_sequences[kj].size(); kr++)
	      {
		  Point pt1 = singular_sequences[kj][kr-1];
		  Point pt2 = singular_sequences[kj][kr];

		  shared_ptr<ParamCurve> paramcurve = shared_ptr<ParamCurve>(new SplineCurve(pt1,pt2));
		  shared_ptr<ParamCurve> dummycrv;
		  ftCurveSegment curr_seg(CURVE_SINGULAR, JOINT_G0, face.get(), 0, paramcurve, dummycrv,
					  dummycrv, toptol_.gap);
		  curr_crv->appendSegment(curr_seg);
	      }
	      singular_curves.push_back(curr_crv);
	      results_->addSingCurve(curr_crv);
	  }
		
      }
  }


  //===========================================================================
  void 
  FaceSetQuality::vanishingCurveTangent(vector<shared_ptr<PointOnCurve> >& sing_points,
					vector<pair<shared_ptr<PointOnCurve>, 
					shared_ptr<PointOnCurve> > >& sing_curves)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(VANISHING_TANGENT, tol) && tol == toptol_.gap)
      {
	  sing_points = results_->getSingCrvPnts();
	  sing_curves = results_->getSingCrvCrvs();
	  return;
      }

      sing_points.clear();
      sing_curves.clear();

      results_->reset(VANISHING_TANGENT);
      results_->performtest(VANISHING_TANGENT, toptol_.gap);

      int nmb_sfs = model_->nmbEntities();
      int ki;

      // Collect all edges
      std::set<shared_ptr<ftEdgeBase> > all_edges;  // All edges in the model represented once
      for (ki=0;  ki<nmb_sfs; ki++)
      {
	  // The function returns existing edges if there are any
	  vector<shared_ptr<ftEdgeBase> > curr_edges = 
	      model_->getFace(ki)->createInitialEdges();
	  all_edges.insert(curr_edges.begin(), curr_edges.end());
      }

      // Check curve singularities
      set<shared_ptr<ftEdgeBase> >::iterator iter = all_edges.begin();
      set<shared_ptr<ftEdgeBase> >::iterator last = all_edges.end();
      for (; iter != last; ++iter)
      {
	  shared_ptr<ParamCurve> crv = (*iter)->geomEdge()->geomCurve();
	  double parmin = (*iter)->tMin();
	  double parmax = (*iter)->tMax();
	  vector<double> singular_pts;
	  vector<vector<double> > singular_sequences;

	  Singular::vanishingTangent(crv, parmin, parmax, toptol_.gap, 
				     singular_pts, singular_sequences);

	  size_t kj;
	  for (kj=0; kj<singular_pts.size(); kj++)
	  {
	      double u = singular_pts[kj];
	      shared_ptr<PointOnCurve> sing_pt = 
		  shared_ptr<PointOnCurve>(new PointOnCurve(crv, u));
	      sing_points.push_back(sing_pt);

	      results_->addSingPntCrv(sing_pt);
	  }

	  for (kj=0; kj<singular_sequences.size(); kj++)
	  {
	      double u1 = singular_sequences[kj][0];
	      double u2 = singular_sequences[kj][singular_sequences[kj].size()-1];
	      shared_ptr<PointOnCurve> sing_pt1 = 
		  shared_ptr<PointOnCurve>(new PointOnCurve(crv, u1));
	      shared_ptr<PointOnCurve> sing_pt2 = 
		  shared_ptr<PointOnCurve>(new PointOnCurve(crv, u2));
	      pair<shared_ptr<PointOnCurve>, shared_ptr<PointOnCurve> > sing_crv =
		  make_pair(sing_pt1, sing_pt2);
	      sing_curves.push_back(sing_crv);
	      results_->addSingCurveCrv(sing_crv);
	  }
		
      }
  }


  //===========================================================================
  void FaceSetQuality::narrowRegion(vector<pair<shared_ptr<PointOnEdge>, 
				  shared_ptr<PointOnEdge> > >& narrow_regions)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(NARROW_REGION, tol) && tol == toptol_.neighbour)
      {
	  narrow_regions = results_->getNarrowRegion();
	  return;
      }

      narrow_regions.clear();

    results_->reset(NARROW_REGION);
    results_->performtest(NARROW_REGION, toptol_.neighbour);

    int nmb_sfs = model_->nmbEntities();
    for (int ki=0; ki<nmb_sfs; ki++)
      {
	  vector<pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > regions;
	  model_->getFace(ki)->getNarrowRegion(toptol_.gap, toptol_.neighbour,
					       regions);
	  for (size_t kr=0; kr<regions.size(); ++kr)
	  {
	      narrow_regions.push_back(regions[kr]);
	      results_->addNarrowRegion(regions[kr]);
	  }
      }

  }

    //===========================================================================
  void FaceSetQuality::sliverSurfaces(vector<shared_ptr<ParamSurface> >& sliver_sfs,
				      double thickness,
				      double factor)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(SLIVER_FACE, tol) && tol == thickness)
      {
	  vector<shared_ptr<ftSurface> > sliver =  results_->getSliverSfs();
	  sliver_sfs.resize(sliver.size());
	  for (size_t ki=0; ki<sliver.size(); ++ki)
	      sliver_sfs[ki] = sliver[ki]->surface();
	  return;
      }

    sliver_sfs.clear();

    results_->reset(SLIVER_FACE);
    results_->performtest(SLIVER_FACE, thickness);

    int nmb_sfs = model_->nmbEntities();
    for (int ki=0; ki<nmb_sfs; ki++)
      {
	bool is_sliver = isSliverFace(model_->getSurface(ki), thickness, factor);
	if (is_sliver)
	  {
	    // Store in result container
	    results_->addSliverSf(model_->getFace(ki));

	    // Return surface
	    sliver_sfs.push_back(model_->getSurface(ki));
	  }
      }
  }


  //===========================================================================
  void FaceSetQuality::edgeVertexDistance(vector<pair<ftEdge*,
					  shared_ptr<Vertex> > >& edge_vertices)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(EDGE_VERTEX_DISTANCE, tol) && tol == toptol_.gap)
      {
	  edge_vertices = results_->getDistantEdgeVertex();
	  return;
      }

    edge_vertices.clear();

    results_->reset(EDGE_VERTEX_DISTANCE);
    results_->performtest(EDGE_VERTEX_DISTANCE, toptol_.gap);

    vector<pair<ftEdge*, shared_ptr<Vertex> > > result;

    int nmb_sfs = model_->nmbEntities();
    for (int i = 0; i < nmb_sfs; ++i)
      model_->getFace(i)->getBadDistance(result, toptol_.gap);

    for (size_t i = 0; i < result.size(); ++i)
      {
	results_->addDistantEdgeVertex(result[i]);
	edge_vertices.push_back(result[i]);
      }

  }


  //===========================================================================
  void FaceSetQuality::faceVertexDistance(vector<pair<ftSurface*,
					  shared_ptr<Vertex> > >& face_vertices)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(FACE_VERTEX_DISTANCE, tol) && tol == toptol_.gap)
      {
	  face_vertices = results_->getDistantFaceVertex();
	  return;
      }

    face_vertices.clear();
    results_->reset(FACE_VERTEX_DISTANCE);
    results_->performtest(FACE_VERTEX_DISTANCE, toptol_.gap);

    vector<pair<ftSurface*, shared_ptr<Vertex> > > result;

    int nmb_sfs = model_->nmbEntities();
    for (int i = 0; i < nmb_sfs; ++i)
      model_->getFace(i)->getBadDistance(result, toptol_.gap);

    for (size_t i = 0; i < result.size(); ++i)
      {
	results_->addDistantFaceVertex(result[i]);
	face_vertices.push_back(result[i]);
      }

  }

  //===========================================================================
  void FaceSetQuality::faceEdgeDistance(vector<pair<ftSurface*,
					ftEdge*> >& face_edges)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(FACE_EDGE_DISTANCE, tol) && tol == toptol_.gap)
      {
	  face_edges = results_->getDistantFaceEdge();
	  return;
      }

    face_edges.clear();

    results_->reset(FACE_EDGE_DISTANCE);
    results_->performtest(FACE_EDGE_DISTANCE, toptol_.gap);

    vector<pair<ftSurface*, ftEdge*> > result;

    int nmb_sfs = model_->nmbEntities();
    for (int i = 0; i < nmb_sfs; ++i)
      model_->getFace(i)->getBadDistance(result, toptol_.gap);

    for (size_t i = 0; i < result.size(); ++i)
      {
	results_->addDistantFaceEdge(result[i]);
	face_edges.push_back(result[i]);
      }

  }


  //===========================================================================
  void FaceSetQuality::edgePosAndTangDiscontinuity(vector<pair<ftEdge*, ftEdge*> >& pos_disconts,
						   vector<pair<ftEdge*, ftEdge*> >& tang_disconts)
  //===========================================================================
  {
       // Check if the test is performed already
      double tol1, tol2;
      if (results_->testPerformed(EDGE_POSITION_DISCONT, tol1) && tol1 == toptol_.gap &&
	  results_->testPerformed(EDGE_TANGENTIAL_DISCONT, tol2) && tol2 == toptol_.kink)
      {
	  pos_disconts = results_->getEdgePosDiscont();
	  tang_disconts = results_->getEdgeTangDiscont();
	  return;
      }

      pos_disconts.clear();
    tang_disconts.clear();

    results_->reset(EDGE_POSITION_DISCONT);
    results_->performtest(EDGE_POSITION_DISCONT, toptol_.gap);

    results_->reset(EDGE_TANGENTIAL_DISCONT);
    results_->performtest(EDGE_TANGENTIAL_DISCONT, toptol_.kink);

    int nmb_sfs = model_->nmbEntities();
    int ki;
    std::set<shared_ptr<Vertex> > all_vertices;  // All vertices in the model represented once

    // Collect all vertices
    for (ki=0;  ki<nmb_sfs; ki++)
    {
	vector<shared_ptr<Vertex> > curr_vertices = model_->getFace(ki)->vertices();
	all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
    }

    // For all vertices, collect attached edges where the distance between the endpoints
    // are larger than the specified tolerance or where the curves meet with an angle that
    // are more than the kink tolerance, but less than the corner tolerance
    set<shared_ptr<Vertex> >::iterator iter = all_vertices.begin();
    set<shared_ptr<Vertex> >::iterator last = all_vertices.end();
    size_t kj;
    for (; iter != last; ++iter)
    {
	vector<pair<ftEdge*, ftEdge*> > gaps;
	vector<pair<ftEdge*, ftEdge*> > kinks;
	(*iter)->getEdgeDiscontinuities(gaps, toptol_.gap, 
					kinks, toptol_.kink, toptol_.bend);
	for (kj=0; kj<gaps.size(); ++kj)
	{
	    pos_disconts.push_back(gaps[kj]);
	    results_->addEdgePositionDiscont(gaps[kj]);
	}
	for (kj=0; kj<kinks.size(); ++kj)
 	{
	    tang_disconts.push_back(kinks[kj]);
	    results_->addEdgeTangentDiscont(kinks[kj]);
	}
   }
  }

  //===========================================================================
  void FaceSetQuality::facePositionDiscontinuity(vector<pair<ftEdge*, ftEdge*> >& pos_disconts)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(FACE_POSITION_DISCONT, tol) && tol == toptol_.gap)
      {
	  pos_disconts = results_->getFacePosDiscont();
	  return;
      }

    pos_disconts.clear();

    results_->reset(FACE_POSITION_DISCONT);
    results_->performtest(FACE_POSITION_DISCONT, toptol_.gap);

    vector<ftEdge*> disconts;
    model_->getGaps(disconts);
    for (size_t ki=0; ki<disconts.size(); ++ki)
    {
	pair<ftEdge*, ftEdge*> gap = make_pair(disconts[ki],disconts[ki]->twin()->geomEdge());
	pos_disconts.push_back(gap);
	results_->addFacePositionDiscont(gap);
    }
  }

  //===========================================================================
  void FaceSetQuality::faceTangentDiscontinuity(vector<pair<ftEdge*, ftEdge*> >& tangent_disconts)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(FACE_TANGENTIAL_DISCONT, tol) && tol == toptol_.kink)
      {
	  tangent_disconts = results_->getFaceTangDiscont();
	  return;
      }

    tangent_disconts.clear();
    results_->reset(FACE_TANGENTIAL_DISCONT);
    results_->performtest(FACE_TANGENTIAL_DISCONT, toptol_.kink);

    vector<ftEdge*> disconts;
    model_->getKinks(disconts);
    for (size_t ki=0; ki<disconts.size(); ++ki)
    {
	pair<ftEdge*, ftEdge*> kink = make_pair(disconts[ki],disconts[ki]->twin()->geomEdge());
	tangent_disconts.push_back(kink);
	results_->addFaceTangentDiscont(kink);
    }
  }


  //===========================================================================
  bool compare_number(pair<ftFaceBase*, int> f1, pair<ftFaceBase*, int> f2)
  //===========================================================================
  {
      return (f1.second >= f2.second);
  }


  //===========================================================================
  void FaceSetQuality::loopOrientationConsistency(vector<shared_ptr<Loop> >& inconsistent_loops)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(LOOP_ORIENTATION, tol))
      {
	  inconsistent_loops = results_->getInconsistLoop();
	  return;
      }

    inconsistent_loops.clear();

    results_->reset(LOOP_ORIENTATION);
    results_->performtest(LOOP_ORIENTATION, 0.0);  // No use of tolerance

    // Check each face
    int nmb_sfs = model_->nmbEntities();
    int nmb_loops = 0;
    for (int kj=0; kj<nmb_sfs; ++kj)
    {
	bool is_consistent = model_->getFace(kj)->checkLoopOrientation(inconsistent_loops);
	if (!is_consistent)
	{
	    for (int ki=nmb_loops; ki<(int)inconsistent_loops.size(); ++ki)
	    {
		results_->addInconsistentBdLoop(inconsistent_loops[ki]);
	    }
	    nmb_loops = (int)inconsistent_loops.size();
	}
    }
    
    return;
  }

  //===========================================================================
  void FaceSetQuality::faceNormalConsistency(vector<shared_ptr<ftSurface> >& inconsistent_faces)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(FACE_ORIENTATION, tol))
      {
	  inconsistent_faces = results_->getInconsistFace();
	  return;
      }

    inconsistent_faces.clear();

    results_->reset(FACE_ORIENTATION);
    results_->performtest(FACE_ORIENTATION, 0.0);  // No use of tolerance

    // Fetch candidate faces found by the topology analysis
    vector<pair<ftFaceBase*, ftFaceBase*> > candidate_faces;
    model_->getInconsistentFacePairs(candidate_faces);

    // Remove pair with an acute angle
    size_t ki, kj;
    for (ki=0; ki<candidate_faces.size(); )
      {
	ftSurface* curr1 = candidate_faces[ki].first->asFtSurface();
	ftSurface* curr2 = candidate_faces[ki].second->asFtSurface();
	if (!curr1 || !curr2)
	  candidate_faces.erase(candidate_faces.begin()+ki);
	else
	  {
// 		std::ofstream debug("orient.g2");
// 		ParamSurface *s1 = curr1->surface().get();
// 		ParamSurface *s2 = curr2->surface().get();
// 		s1->writeStandardHeader(debug);
// 		s1->write(debug);
// 		s2->writeStandardHeader(debug);
// 		s2->write(debug);

	    shared_ptr<ftEdge> edge1, edge2;
	    bool neighbours = curr1->areNeighbours(curr2, edge1, edge2);
	    if (!neighbours)
	      candidate_faces.erase(candidate_faces.begin()+ki);
	    else
	      {
		bool acute_angle = 
		  curr1->hasAcuteAngle(edge1.get(), toptol_.kink);
		if (acute_angle)
		  candidate_faces.erase(candidate_faces.begin()+ki);
		else
		  ki++;
	      }
	  }
      }

//     std::ofstream out_file("candiate_sfs.g2");
//     size_t nmb_faces = candidate_faces.size();
//     for (size_t kr=0; kr<nmb_faces; ++kr)
//       {
// 	shared_ptr<ParamSurface> sf = 
// 	  candidate_faces[kr].first->asFtSurface()->surface();
// 	sf->writeStandardHeader(out_file);
// 	sf->write(out_file);
// 	sf = candidate_faces[kr].second->asFtSurface()->surface();
// 	sf->writeStandardHeader(out_file);
// 	sf->write(out_file);
//       }

   // Fetch the faces participating in most entries
    vector<ftFaceBase*> faces;
    while (candidate_faces.size() > 0)
    {
	vector<pair<ftFaceBase*, int> > face_occurances;
	for (ki=0; ki<candidate_faces.size(); ++ki)
	{
	    bool f1 = false, f2 = false;
	    for (kj=0; kj<face_occurances.size(); ++kj)
	    {
		if (candidate_faces[ki].first == face_occurances[kj].first)
		{
		    face_occurances[kj].second++;
		    f1 = true;
		}

		if (candidate_faces[ki].second == face_occurances[kj].first)
		{
		    face_occurances[kj].second++;
		    f2 = true;
		}
	    }
	    if (!f1)
		face_occurances.push_back(make_pair(candidate_faces[ki].first, 1));
	    if (!f2)
		face_occurances.push_back(make_pair(candidate_faces[ki].second, 1));
	}

	//	std::sort(face_occurances.begin(), face_occurances.end(), compare_number);

	// Switch meaning of number to number of consistent adjacence
	for (ki=0; ki<face_occurances.size(); ++ki)
	  {
	    vector<ftSurface*> adj;
	    face_occurances[ki].first->asFtSurface()->getAdjacentFaces(adj);
	    face_occurances[ki].second = (int)adj.size() - face_occurances[ki].second;
	  }

	for (ki=0; ki<face_occurances.size(); ++ki)
	  for (kj=ki+1; kj<face_occurances.size(); ++kj)
	    {
	      if (face_occurances[ki].second > face_occurances[kj].second)
		std::swap(face_occurances[ki], face_occurances[kj]);
	    }

	// Save the face with most occurances
	ftFaceBase* curr_face = face_occurances[0].first;
	faces.push_back(curr_face);

	// Remove all entries containing the current face
	for (ki=0; ki<candidate_faces.size(); )
	{
	    if (candidate_faces[ki].first == curr_face || 
		candidate_faces[ki].second == curr_face)
		candidate_faces.erase(candidate_faces.begin() + ki);
	    else
		ki++;
	}
    }

    // Fetch the faces corresponding to the found occurances and store them
    int nmb_sfs = model_->nmbEntities();
    int kr;
    for (ki=0; ki<faces.size(); ++ki)
    {
	for (kr=0; kr<nmb_sfs; ++kr)
	    if (model_->getFace(kr).get() == faces[ki])
	    {
		results_->addInconsistentFaceInSet(model_->getFace(kr));
		inconsistent_faces.push_back(model_->getFace(kr));
	    }
    }
  }

  //===========================================================================
  void FaceSetQuality::sfG1Discontinuity(vector<shared_ptr<ftSurface> >& discont_sfs)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(SF_G1DISCONT, tol) && tol == toptol_.kink)
      {
	  discont_sfs = results_->getG1DiscontSfs();
	  return;
      }

    discont_sfs.clear();

    results_->reset(SF_G1DISCONT);
    results_->performtest(SF_G1DISCONT, toptol_.kink);

    int nmb_sfs = model_->nmbEntities();
    for (int ki = 0; ki < nmb_sfs; ++ki)
    {
	vector<double> g1_disc_u, g1_disc_v;
	shared_ptr<ftSurface> curr_face = model_->getFace(ki);
	bool has_kinks = curr_face->getSurfaceKinks(toptol_.kink, g1_disc_u, g1_disc_v);
	if (has_kinks)
	{
	    discont_sfs.push_back(curr_face);
	    results_->addG1DiscontSf(curr_face);
	}
    }
  }
    
  //===========================================================================
  void FaceSetQuality::sfC1Discontinuity(vector<shared_ptr<ftSurface> >& discont_sfs)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(SF_C1DISCONT, tol) && tol == toptol_.gap)
      {
	  discont_sfs = results_->getC1DiscontSfs();
	  return;
      }

    discont_sfs.clear();

    results_->reset(SF_C1DISCONT);
    results_->performtest(SF_C1DISCONT, toptol_.gap);

    int nmb_sfs = model_->nmbEntities();
    for (int ki = 0; ki < nmb_sfs; ++ki)
    {
	vector<double> g1_disc_u, g1_disc_v;
	shared_ptr<ftSurface> curr_face = model_->getFace(ki);
	bool has_kinks = curr_face->getSurfaceDisconts(toptol_.gap, g1_disc_u, 
						       g1_disc_v);
	if (has_kinks)
	{
	    discont_sfs.push_back(curr_face);
	    results_->addC1DiscontSf(curr_face);
	}
    }
  }
    
    //===========================================================================
  void FaceSetQuality::cvC1G1Discontinuity(vector<shared_ptr<ParamCurve> >& c1_discont,
					 vector<shared_ptr<ParamCurve> >& g1_discont)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol1, tol2;
      if (results_->testPerformed(CV_C1DISCONT, tol1) && tol1 == toptol_.gap &&
	  results_->testPerformed(CV_G1DISCONT, tol2) && tol2 == toptol_.kink)
      {
	  c1_discont = results_->getC1DiscontCvs();
	  g1_discont = results_->getG1DiscontCvs();
	  return;
      }
      c1_discont.clear();
      g1_discont.clear();

    results_->reset(CV_C1DISCONT);
    results_->performtest(CV_C1DISCONT, toptol_.gap);
    results_->reset(CV_G1DISCONT);
    results_->performtest(CV_G1DISCONT, toptol_.kink);

    int nmb_sfs = model_->nmbEntities();
    int ki;
    size_t kr;
    // Collect all edges
    std::set<shared_ptr<ParamCurve> > all_crvs;  // All curves in the model represented once
    for (ki=0;  ki<nmb_sfs; ki++)
    {
	// The function returns existing edges if there are any
	vector<shared_ptr<ftEdgeBase> > curr_edges = 
	    model_->getFace(ki)->createInitialEdges();
	for (kr=0; kr<curr_edges.size(); ++kr)
	    all_crvs.insert(curr_edges[kr]->geomEdge()->geomCurve());
    }

    // Check edge length
    set<shared_ptr<ParamCurve> >::iterator iter = all_crvs.begin();
    set<shared_ptr<ParamCurve> >::iterator last = all_crvs.end();
    for (; iter != last; ++iter)
    {
	// Get curve
	SplineCurve *spline = (*iter)->geometryCurve();
	if (!spline)
	    continue;

	vector<double> c1disconts;
	vector<double> g1disconts;
	curveKinks(*spline, toptol_.gap, toptol_.kink, c1disconts, g1disconts);
	if (c1disconts.size() > 0)
	{
	    c1_discont.push_back(*iter);
	    results_->addC1DiscontCv(*iter);
	}
	if (g1disconts.size() > 0)
	{
	    g1_discont.push_back(*iter);
	    results_->addG1DiscontCv(*iter);
	}
    }   

  }

  //===========================================================================
  void FaceSetQuality::cvCurvatureRadius(vector<pair<shared_ptr<PointOnCurve>, double> >& small_curv_rad,
					 pair<shared_ptr<PointOnCurve>, double>& minimum_curv_rad)
  //===========================================================================
  {
      // Check if the test is performed already
      double tol;
      if (results_->testPerformed(CV_CURVATURE_RADIUS, tol) && tol == curvature_radius_)
      {
	  small_curv_rad = results_->getSmallCvCurvatureR();
	  minimum_curv_rad = results_->getMinCvCurvatureR();
	  return;
      }

      shared_ptr<PointOnCurve> min_pos;
      small_curv_rad.clear();
      double min_rad = MAXDOUBLE;

      results_->reset(CV_CURVATURE_RADIUS);
      results_->performtest(CV_CURVATURE_RADIUS,curvature_radius_); 

    // Collect all edges
    int nmb_sfs = model_->nmbEntities();
    std::set<shared_ptr<ParamCurve> > all_crvs;  // All curves in the model represented once
    for (int ki=0;  ki<nmb_sfs; ki++)
    {
	// The function returns existing edges if there are any
	vector<shared_ptr<ftEdgeBase> > curr_edges = 
	    model_->getFace(ki)->createInitialEdges();
	for (size_t kr=0; kr<curr_edges.size(); ++kr)
	    all_crvs.insert(curr_edges[kr]->geomEdge()->geomCurve());
    }


    double mincurv, param;
    set<shared_ptr<ParamCurve> >::iterator iter = all_crvs.begin();
    set<shared_ptr<ParamCurve> >::iterator last = all_crvs.end();
    for (; iter != last; ++iter)
    {
	// Get curve
	SplineCurve *spline = (*iter)->geometryCurve();
	if (!spline)
	    continue;

	minimalCurvatureRadius(*spline, mincurv, param);
	if (mincurv < min_rad)
	{
	    min_rad = mincurv;
	    min_pos = shared_ptr<PointOnCurve>(new PointOnCurve(*iter, param));
	}
	if (mincurv < curvature_radius_)
	{
	    shared_ptr<PointOnCurve> curr_point 
	    = shared_ptr<PointOnCurve>(new PointOnCurve(*iter, param));
	    pair<shared_ptr<PointOnCurve>, double> curr_rad = make_pair(curr_point, mincurv);
	    small_curv_rad.push_back(curr_rad);
	    results_->smallCvCurvRad(curr_rad);
	}
    }
    minimum_curv_rad = make_pair(min_pos, min_rad);
    results_->setMinimumCvCurvatureRadius(minimum_curv_rad);
  }

  //===========================================================================
  void FaceSetQuality::sfCurvatureRadius(vector<pair<shared_ptr<ftPoint>, double> >& small_curv_rad,
					 pair<shared_ptr<ftPoint>, double>& minimum_curv_rad)
  //===========================================================================
  {
      double tol;
      if (results_->testPerformed(SF_CURVATURE_RADIUS, tol) && tol == curvature_radius_)
      {
	  small_curv_rad = results_->getSmallSfCurvatureR();
	  minimum_curv_rad = results_->getMinSfCurvatureR();
	  return;
      }

      shared_ptr<ftPoint> min_pos;
      small_curv_rad.clear();
      double min_rad = MAXDOUBLE;

      results_->reset(SF_CURVATURE_RADIUS);
      results_->performtest(SF_CURVATURE_RADIUS, curvature_radius_);

    int nmb_sfs = model_->nmbEntities();
    double mincurv, par_u, par_v;
    for (int ki = 0; ki < nmb_sfs; ++ki)
    {
	shared_ptr<ftSurface> face = model_->getFace(ki);
	shared_ptr<ParamSurface> surf = face->surface();

	minimalCurvatureRadius(*surf, curvature_radius_, mincurv, par_u, par_v, toptol_.gap);
	if (mincurv < min_rad)
	{
	    min_rad = mincurv;
	    Point pos = surf->point(par_u, par_v);
	    min_pos = shared_ptr<ftPoint>(new ftPoint(pos, face.get(), par_u, par_v));
	}
	if (mincurv < curvature_radius_)
	{
	    Point pos = surf->point(par_u, par_v);
	    shared_ptr<ftPoint> curr_ftpoint 
		= shared_ptr<ftPoint>(new ftPoint(pos, face.get(), par_u, par_v));
	    pair<shared_ptr<ftPoint>, double> curr_rad = make_pair(curr_ftpoint, mincurv);
	    small_curv_rad.push_back(curr_rad);
	    results_->smallSfCurvRad(curr_rad);
	}
    }

    minimum_curv_rad = make_pair(min_pos, min_rad);
    results_->setMinimumCurvatureRadius(minimum_curv_rad);
  }

  //===========================================================================
    void FaceSetQuality::acuteEdgeAngle(vector<pair<ftEdge*, ftEdge*> >& edge_acute)
  //===========================================================================
    {
      double tol;
      if (results_->testPerformed(EDGE_ACUTE_ANGLE, tol) && tol == toptol_.kink)
      {
	  edge_acute = results_->getEdgeAcuteAngle();
	  return;
      }

	edge_acute.clear();

	results_->reset(EDGE_ACUTE_ANGLE);
	results_->performtest(EDGE_ACUTE_ANGLE, toptol_.kink);

	int nmb_sfs = model_->nmbEntities();
	int ki, kj;

	for (ki=0; ki<nmb_sfs; ++ki)
	{
	    shared_ptr<ftSurface> curr_face = model_->getFace(ki);
	    int nmb_loop = curr_face->nmbBoundaryLoops();
	    for (kj=0; kj<nmb_loop; ++kj)
	    {
		vector<pair<ftEdge*, ftEdge*> > acute_edges;
		curr_face->getBoundaryLoop(kj)->getAcuteEdges(acute_edges, toptol_.kink);

		for (size_t kr=0; kr<acute_edges.size(); ++kr)
		{
		    edge_acute.push_back(acute_edges[kr]);
		    results_->addEdgeAcuteAngle(acute_edges[kr]);
		}
	    }
	}

    }

  //===========================================================================
    void FaceSetQuality::acuteFaceAngle(vector<pair<ftSurface*, ftSurface* > >& face_acute)
  //===========================================================================
    {
      double tol;
      if (results_->testPerformed(FACE_ACUTE_ANGLE, tol) && tol == toptol_.kink)
      {
	  face_acute = results_->getFaceAcuteAngle();
	  return;
      }

	face_acute.clear();

	results_->reset(FACE_ACUTE_ANGLE);
	results_->performtest(FACE_ACUTE_ANGLE, toptol_.kink);

	// The candiate acute angles are along corners in the model
	vector<ftEdge*> corners;
	model_->getCorners(corners);

	for (size_t ki=0; ki<corners.size(); ++ki)
	{
	    ftSurface* curr_face = corners[ki]->face()->asFtSurface();
	    if (!curr_face)
		continue;   // Not handled

	    // Check current face
	    bool acute_angle = curr_face->hasAcuteAngle(corners[ki], toptol_.kink);
	    if (acute_angle)
	    {
		pair<ftSurface*, ftSurface*> acute_faces = 
		    make_pair(curr_face, corners[ki]->twin()->face()->asFtSurface());
		face_acute.push_back(acute_faces);
		results_->addFaceAcuteAngle(acute_faces);
	    }
	}
    }

    //===========================================================================
    void FaceSetQuality::loopIntersection(vector<pair<shared_ptr<PointOnEdge>, 
					shared_ptr<PointOnEdge> > >& loop_intersection)
  //===========================================================================
    {
      double tol;
      if (results_->testPerformed(LOOP_INTERSECTION, tol) && tol == toptol_.gap)
      {
	  loop_intersection = results_->getIntersectingBdLoops();
	  return;
      }

	loop_intersection.clear();

	results_->reset(LOOP_INTERSECTION);
	results_->performtest(LOOP_INTERSECTION, toptol_.gap);

	int nmb_sfs = model_->nmbEntities();
	int ki, kj, kh;

	for (ki=0; ki<nmb_sfs; ++ki)
	{
	    shared_ptr<ftSurface> curr_face = model_->getFace(ki);
	    int nmb_loop = curr_face->nmbBoundaryLoops();
	    for (kj=0; kj<nmb_loop; ++kj)
	    {
		shared_ptr<Loop> loop1 = curr_face->getBoundaryLoop(kj);
		for (kh=kj+1; kh<nmb_loop; ++kh)
		{
		    vector<pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > int_pt;
		    shared_ptr<Loop> loop2 = curr_face->getBoundaryLoop(kh);

		    // DEBUG. Draw loop
// 		    std::ofstream out_file("int_loops.g2");
// 		    for (size_t i2=0; i2<loop1->size(); ++i2)
// 		      {
// 			shared_ptr<ftEdgeBase> e1 = loop1->getEdge(i2);
// 			shared_ptr<ParamCurve> crv = e1->geomEdge()->geomCurve();
// 			SplineCurve* spline = crv->geometryCurve();
// 			if (spline)
// 			  {
// 			    spline->writeStandardHeader(out_file);
// 			    spline->write(out_file);
// 			  }
// 		      }
// 		    for (size_t i2=0; i2<loop2->size(); ++i2)
// 		      {
// 			shared_ptr<ftEdgeBase> e1 = loop2->getEdge(i2);
// 			shared_ptr<ParamCurve> crv = e1->geomEdge()->geomCurve();
// 			SplineCurve* spline = crv->geometryCurve();
// 			if (spline)
// 			  {
// 			    spline->writeStandardHeader(out_file);
// 			    spline->write(out_file);
// 			  }
// 		      }

		    loop1->getLoopIntersections(loop2, toptol_.gap, int_pt);

		    for (size_t kr=0; kr<int_pt.size(); ++kr)
		    {
			loop_intersection.push_back(int_pt[kr]);
			results_->addLoopIntersection(int_pt[kr]);
		    }
		}
	    }
	}
    }

    //===========================================================================
    void FaceSetQuality::loopSelfIntersection(vector<pair<shared_ptr<PointOnEdge>, 
					      shared_ptr<PointOnEdge> > >& loop_self_intersection)
  //===========================================================================
    {
      double tol;
      if (results_->testPerformed(LOOP_SELF_INTERSECTION, tol) && tol == toptol_.gap)
      {
	  loop_self_intersection = results_->getSelfIntersectingBdLoops();
	  return;
      }

	loop_self_intersection.clear();


	results_->reset(LOOP_SELF_INTERSECTION);
	results_->performtest(LOOP_SELF_INTERSECTION, toptol_.gap);

	int nmb_sfs = model_->nmbEntities();
	int ki, kj;

	for (ki=0; ki<nmb_sfs; ++ki)
	{
	    shared_ptr<ftSurface> curr_face = model_->getFace(ki);
	    int nmb_loop = curr_face->nmbBoundaryLoops();
	    for (kj=0; kj<nmb_loop; ++kj)
	    {
		shared_ptr<Loop> loop1 = curr_face->getBoundaryLoop(kj);
		vector<pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > int_pt;
		loop1->getLoopSelfIntersections(toptol_.gap, int_pt);
		for (size_t kr=0; kr<int_pt.size(); ++kr)
		{
		    loop_self_intersection.push_back(int_pt[kr]);
		    results_->addLoopSelfIntersection(int_pt[kr]);
		}
	    }
	}
    }


  //===========================================================================
    void FaceSetQuality::indistinctKnots(vector<shared_ptr<ParamCurve> >& cv_knots,
				       vector<shared_ptr<ParamSurface> >& sf_knots,
				       double tol)
  //===========================================================================
    {
      double tol2;
      if (results_->testPerformed(INDISTINCT_KNOTS, tol2) && tol2 == tol)
      {
	  cv_knots = results_->getCvIndistinctKnots();
	  sf_knots = results_->getSfIndistinctKnots();
	  return;
      }
	cv_knots.clear();
	sf_knots.clear();

	results_->reset(INDISTINCT_KNOTS);
	results_->performtest(INDISTINCT_KNOTS, tol);

	int nmb_sfs = model_->nmbEntities();
	int ki;

	for (ki=0; ki<nmb_sfs; ++ki)
	{
	    shared_ptr<ParamSurface> surf = model_->getSurface(ki);
	    vector<shared_ptr<ParamCurve> > trim_knots;
	    bool has_indistinct_knots = qualityUtils::hasIndistinctKnots(surf, tol, trim_knots);
	    if (has_indistinct_knots)
	    {
		sf_knots.push_back(surf);
		results_->addSfIndistinctKnot(surf);
	    }
	    for (size_t kr=0; kr<trim_knots.size(); ++kr)
	    {
		cv_knots.push_back(trim_knots[kr]);
		results_->addCvIndistinctKnot(trim_knots[kr]);
	    }
	}
    }

} // namespace Go
