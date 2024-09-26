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

#include "GoTools/compositemodel/ftSurfaceSet.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftSSfEdge.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/cmUtils.h"
#include "GoTools/creators/SurfaceCreators.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/compositemodel/ftSmoothSurf.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/topology/FaceAdjacency.h"

#include "GoTools/compositemodel/cmUtils.h"

//#include "GoTools/model_toolbox/ftUtils.h"

#include <cstdio> // for debugging

#ifndef NDEBUG
// #define DEBUG
// #define FANTASTIC_DEBUG
#endif

using std::vector;
using std::min;
using std::max;

namespace Go
{

//===========================================================================
ftSurfaceSet::ftSurfaceSet(const vector<shared_ptr<ParamSurface> >& surfaces,
			   tpTolerances& topeps, double approxeps)
  : ftFaceBase(-1), toptol_(topeps),
    approxtol_(approxeps), approxweight_(1.0-1.0e-6)
//===========================================================================
{
  ftFaceBase* f = 0;
  for (size_t i=0; i<surfaces.size(); ++i)
  {
      f = new ftSurface(surfaces[i], (int)i);
    faces_.push_back(shared_ptr<ftFaceBase>(f));
  }
  buildTopology();
}

//===========================================================================
ftSurfaceSet::~ftSurfaceSet()
//===========================================================================
{
}

//===========================================================================
void
ftSurfaceSet::getSurfaces(vector<shared_ptr<ParamSurface> >& surfaces)
//---------------------------------------------------------------------------
//
// Purpose: Return the surfaces of the surface set
//
//===========================================================================
{
  int no_faces = (int)faces_.size();
  surfaces.resize(no_faces);
  for (int ki=0; ki<no_faces; ki++)
    surfaces[ki] = faces_[ki]->surface();
}

// //---------------------------------------------------------------------------
// vector<shared_ptr<ftEdgeBase> >
// ftSurfaceSet::setOrientation(double degenerate_epsilon)
// //---------------------------------------------------------------------------
// {
//   if (is_turned_)
//     {
//       if (surf_.get() != 0)
// 	surf_->turnOrientation();
//       for (size_t ki=0; ki<faces_.size(); ki++)
// 	  faces_[ki]->turnOrientation();
//       buildTopology();
//       is_turned_ = false;
//     }
//   return createInitialEdges(degenerate_epsilon);
// }

//---------------------------------------------------------------------------
vector<shared_ptr<ftEdgeBase> > ftSurfaceSet::startEdges()
//---------------------------------------------------------------------------
{
    vector<shared_ptr<ftEdgeBase> > first_edges(boundary_loops_.size());
    for (size_t ki=0; ki<boundary_loops_.size(); ki++)
	first_edges[ki] = boundary_loops_[ki]->getEdge(0);

    return first_edges;
}

//---------------------------------------------------------------------------
BoundingBox ftSurfaceSet::boundingBox()
//---------------------------------------------------------------------------
{
  if (surf_.get() != 0)
    return surf_->boundingBox();
  else
    {
	int nface = (int)faces_.size();
      BoundingBox box = faces_[0]->boundingBox();
      for (int ki=1; ki<nface; ki++)
	box.addUnionWith(faces_[ki]->boundingBox());
      return box;
    }
}

//===========================================================================
void ftSurfaceSet::setId(int id)
//===========================================================================
{
    id_ = id;
}

//===========================================================================
int ftSurfaceSet::getId()
//===========================================================================
{
    return id_;
}

// //---------------------------------------------------------------------------
// void ftSurfaceSet::turnOrientation()
// //---------------------------------------------------------------------------
// {
//     is_turned_ = (is_turned_) ? false : true;
// }


// //---------------------------------------------------------------------------
// bool ftSurfaceSet::getOrientation()
// //---------------------------------------------------------------------------
// {
//     return is_turned_;
// }

//===========================================================================
ftMessage ftSurfaceSet::trimWithPlane(const ftPlane& plane)
//===========================================================================
{
    ftMessage status, local_status;

    double epsge = toptol_.neighbour;
    int ki, kj;
    for (ki = 0; ki < (int)faces_.size(); ++ki) {
	shared_ptr<ParamSurface> surf = faces_[ki]->surface();

	bool boundingbox_intersected = plane.intersectsBox(surf->boundingBox());
	bool surface_trimmed = false;

	if (boundingbox_intersected) {
	    // As we are about to trim surface, we convert to a BoundedSurface.

	    // We intersect surf with plane. Routine expects a SplineSurface.
	    vector<shared_ptr<BoundedSurface> > trimmed_surfs;
	    try {
		// @@sbr Geom epsilon (0.01*epsge) must be given a suitable value.
		// Should match up with epsilon of boundary loop of bounded surface.
		trimmed_surfs = BoundedUtils::trimWithPlane(surf, plane.point(),
							      plane.normal(), 0.01*epsge);
	    } catch (...) {
		status.addWarning(FT_ERROR_IN_SURFACE_TRIMMING);
	    }
	    if (trimmed_surfs.size() != 0) { // surf was not intersected by plane.
		for (kj = 0; kj < (int)trimmed_surfs.size(); ++kj)
		    faces_.push_back(shared_ptr<ftFaceBase>
				    (new ftSurface(trimmed_surfs[kj], faces_[ki]->getId())));
		faces_.erase(faces_.begin() + ki);
		surface_trimmed = true;
	    }
	}

	// If plane did not intersect surf, we check on which side surf lies.
	// We compute inner product of normal and vector given by a plane point
	// and a surface_pt. Positive inner product <=> surf is above plane.
	if (!surface_trimmed) {
	    Vector2D pt(0.0, 0.0); // We choose a random point.
	    Vector2D domain_pt;
	    surf->parameterDomain().closestInDomain(pt, domain_pt, epsge);
	    Point surf_pt = surf->point(domain_pt[0], domain_pt[1]);
	    Point normal = plane.normal();
	    if (cmUtils::abovePlane(surf_pt, plane.point(), plane.normal())) {
		faces_.erase(faces_.begin() + ki);
		--ki;
	    }
	}
    }

//     // For debugging
//     Point nrm;
//     std::ofstream fileout("data/output/trimmed_boundaries.g2");
//     int j;
//     shared_ptr<ParamSurface> face;
//     for (i = 0; i < faces_.size(); ++i) {
//  	face = faces_[i]->Surface();

// 	CurveLoop loop;
// 	    try {
// 		loop = face->outerBoundaryLoop(epsge);
// 	    } catch (...) {
// 		MESSAGE("Failed extracting outer loop! Something is very wrong!");
// 		continue;
// 	    }

// 	for (j = 0; j < loop.size(); ++j) {
// 	    shared_ptr<ParamCurve> spacecurve =
// 		((face.do_dynamic_cast<BoundedSurface>()).get() != 0) ?
// 		loop[j].do_dynamic_cast<CurveOnSurface>()->geomCurve() :
// 		loop[j];
//  	    spacecurve->writeStandardHeader(fileout);
//  	    spacecurve->write(fileout);
// 	}
//     }
//     // end of debugging

    // Set face-id corresponding to the array index
    for (ki=0; ki<(int)faces_.size(); ki++)
      faces_[ki]->setId(ki);

    buildTopology();

    return status;
}

//===========================================================================
ftMessage
ftSurfaceSet::getBoundaryConditions(const vector<ftEdgeBase*>& edgeloop, 
				    vector<int>& corner, 
				    vector< shared_ptr<SplineCurve> >& 
				    bd_curves, 
				    vector< shared_ptr<SplineCurve> >&
				    cross_curves,
				    bool compute_cross_curves,
				    vector<BoundaryPiece>& bdpiece,
				    bool require_cord_length_param)
//===========================================================================
{
  ftMessage status;

  // Do also make an estimate of the curve length of each boundary curve.
  int ki, kj, kh, kk;
  double length, minlength = 10000000000.0;
  int idxmin = 0;
  ftFaceBase *adjsurf;
//   double tmin, tmax;
  Point ptmin, ptmax;
  SplineCurve *cv = 0, *crosscv = 0; 
  shared_ptr<SplineCurve> cv2;
  ftSurfaceSet *adjssf;
  ftSurface *adjface;
  ftTangPriority adj_prio = ftNoType;
  shared_ptr<SplineCurve> dummycrv;  // Empty pointer to curve
  double par1, par2;
  int bdindex;
  vector<bool> local_at_end(2*(corner.size()-1), false);

  // Copy edgeloop array for local use
  int nmbedges = (int)edgeloop.size();
  vector<ftEdgeBase*> edges;
  edges.insert(edges.end(), edgeloop.begin(), edgeloop.end());
  if (nmbedges < corner[corner.size()-1])
    {
      // Extend the edge array.
      edges.insert(edges.end(), edgeloop.begin(),
		   edgeloop.begin()+corner[corner.size()-1]-nmbedges);
    }

  vector<vector<shared_ptr<SplineCurve> > > all_bd_curves;
  vector<vector<shared_ptr<SplineCurve> > > all_cross_curves;
  vector<vector<bool> > curve_from_neighbour;

  // For each boundary of the current surface ...
  for (ki=0; ki<(int)corner.size()-1; ki++)
    {
      // Expecting that all parts of edge belong to created twin faces, or none.
        bool local_edge;
        local_edge = true;
      // Boundary and cross tangent curves belonging to the current
      // boundary of the current surface. The pieces are stored in 
      // local vectors.
      vector<shared_ptr<SplineCurve> > curr_bd;
      vector<shared_ptr<SplineCurve> > curr_cross;
      vector<bool> crv_from_neighbour;

      // Fetch boundary information from the various pieces along 
      // this boundary. Try to avoid unnecessary fragmentation of
      // one curve.

//       ftSSfEdge* s_edge;
//       ftEdge* g_edge;
      kj = corner[ki];
      length = 0.0;

      adjsurf = 0;

      vector<ftEdgeBase*> local_edges;
      if (corner[ki+1] < corner[ki]) {
	  local_edges.insert(local_edges.end(), edges.begin() + corner[ki], edges.end());
	  local_edges.insert(local_edges.end(), edges.begin(), edges.begin() + corner[ki+1]);
      } else {
	  local_edges.insert(local_edges.end(),
			     edges.begin() + corner[ki], edges.begin() + corner[ki+1]);
      }
      std::vector<ftEdgeBase*>::iterator first_edge = local_edges.begin();
// 	  = edges.begin() + corner[ki];
      std::vector<ftEdgeBase*>::iterator last_edge = local_edges.end();
// 	  = edges.begin() + corner[ki+1];


//      std::vector<ftEdgeBase*>::iterator first_edge = edges.begin() + corner[ki];
//      std::vector<ftEdgeBase*>::iterator last_edge = edges.begin() + corner[ki+1];

      while (first_edge != last_edge)
	{
	  // As long as the boundary conditions are uniform (no neighbouring 
	  // surface or the same neighbouring surface), compute the size 
	  // of the current type of condition, and return the condition.
	  double len = 0.0;
	  status = getNextEdgePieceInfo(first_edge, last_edge, adjsurf,
					cv, ptmin, ptmax, len);
	  length += len;

	  // Fetch the current boundary information
	  if (adjsurf)
	    {
	      crosscv = cv = 0;
	      local_edge = false;
	      // Fetch the boundary and cross boundary curve
	      // from the neighbouring surface.
	      SplineSurface *g_surf = 
		  dynamic_cast<SplineSurface*>(adjsurf->surface().get());
	      if (g_surf != 0)
		{
//  		  g_surf->getBoundaryIdx(ptmin, ptmax, toptol_.gap,
//  					 bdindex, par1, par2);

#ifdef FANTASTIC_DEBUG
		    std::ofstream debug("data/debug.g2");
		    g_surf->writeStandardHeader(debug);
		    g_surf->write(debug);
		    SplineCurve debug_crv(ptmin, ptmax);
		    debug_crv.writeStandardHeader(debug);
		    debug_crv.write(debug);
#endif // FANTASTIC_DEBUG

		  g_surf->getBoundaryIdx(ptmin, ptmax, 10.0*approxtol_,
					 bdindex, par1, par2);
		  ALWAYS_ERROR_IF(bdindex < 0, "No boundary information");

		  g_surf->getBoundaryInfo(par1, par2, bdindex,
					  cv, crosscv);
		}

	      // Check if the neighbouring face is of type 
	      // ftSuperSurface. In that case fetch the common
	      // boundary curve from this face.
	      adjssf = dynamic_cast<ftSurfaceSet*>(adjsurf);
	      adjface = dynamic_cast<ftSurface*>(adjsurf);
	      if ((adjssf != 0) || (adjface != 0))
		{
		  if (adjssf != 0)
		    {
		      cv2 = adjssf->getBoundaryPiece(ptmin, ptmax, 
						     toptol_.neighbour);
		      adj_prio = adjssf->getPrioType();
		    }
		  else
		    {
		      cv2 = adjface->getBoundaryPiece(ptmin, ptmax, 
						     toptol_.neighbour);
		      adj_prio = adjface->getPrioType();
		    }

		  double minpar, maxpar;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
		  minpar = min(par1,par2);
		  maxpar = max(par1,par2);
#else
		  minpar = std::min(par1,par2);
		  maxpar = std::max(par1,par2);		  
#endif
		  
		  bdpiece.push_back(std::make_pair(std::make_pair(minpar, maxpar),
						   std::make_pair(adjsurf, bdindex)));
		}

	      if (cv2.get() != 0)
		{
		  curr_bd.push_back(cv2);
		  if (crosscv != 0)
		    crosscv->setParameterInterval(cv2->startparam(),
						  cv2->endparam());
		}
	      else if (cv != 0)
		curr_bd.push_back(shared_ptr<SplineCurve>(cv));
	      if (compute_cross_curves && crosscv != 0 && 
		  getPrioType() == ftSlave && adj_prio == ftMaster)
		curr_cross.push_back(shared_ptr<SplineCurve>(crosscv));
	      else 
		curr_cross.push_back(dummycrv);
	      crv_from_neighbour.push_back(true);
	    }
	  else
	    {
	      // Fetch the boundary information from the current face.
	      // The edge may have been split.
// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
// 	      cv = dynamic_cast<SplineCurve*>(cv->subCurve(tmin, tmax));
// #else
// 	      cv = cv->subCurve(tmin, tmax);
// #endif

// 	      ALWAYS_ERROR_IF(cv==0, "Curve is of wrong kind!", CorruptData());

	      if (curr_bd.size() == 0)
		local_at_end[2*ki] = true;
	      if (first_edge == last_edge)
		local_at_end[2*ki+1] = true;
	      curr_bd.push_back(shared_ptr<SplineCurve>(cv));
	      curr_cross.push_back(dummycrv);
	      crv_from_neighbour.push_back(false);
	    }

	}  // end while

      // Join the various curves making up this boundary into one curve.
      // Boundary curves are to be joined with C1 continuity and
      // Cross boundary curves with C(-1) continuity.
      ALWAYS_ERROR_IF(curr_bd.size() == 0, "Missing boundary curve");


      // Look for minimum edge length.
      if (length < minlength)
	{
	  minlength = length;
	  idxmin = ki;
	}

//       // Make sure that already created curves (approximations) taken from neighbour
//       // have cord length parametrization. Maybe they should be stored that way?
//       for (kj = 0; kj < curr_bd.size(); ++kj)
// 	if (crv_from_neighbour[kj])
// 	  curr_bd[kj]->setParameterInterval(curr_bd[kj]->startparam(),
// 					    curr_bd[kj]->startparam() + curr_bd[kj]->
// 					    ParamCurve::estimatedCurveLength());

      all_bd_curves.push_back(curr_bd);
      all_cross_curves.push_back(curr_cross);
      curve_from_neighbour.push_back(crv_from_neighbour);
    }

//   // @@ At this instance we're ready to split edge cvs if additional_corner_pts_ exist.
//   // Should be done before edges are approximated.
//   double gap = toptol_.gap;
//   ftUtils::splitBdCvs(all_bd_curves, all_cross_curves, curve_from_neighbour,
// 		      additional_corner_pts_, corner, gap);

  int nmb_sides = (int)all_bd_curves.size();
  double cross_tol = 1e-05; // If less than this value, do not use end tangent.
  // We then loop through all curves and join segments.
  for (ki = 0; ki < (int)all_bd_curves.size(); ++ki) {
      // We estimate length of sum of segments.
      double sum_segments = 0.0;
      for (kk = 0; kk < (int)all_bd_curves[ki].size(); ++kk) {
	  sum_segments +=
	      all_bd_curves[ki][kk]->ParamCurve::estimatedCurveLength();
      }

      // Create missing pieces of the boundary curve by approximation
      for (kj=0; kj<(int)all_bd_curves[ki].size(); kj++) {
	  kh=kj+1;
	  if (curve_from_neighbour[ki][kj]) {
	      continue;   // Already a final curve.
	  }

	  for (; kh<(int)all_bd_curves[ki].size() && !curve_from_neighbour[ki][kh]; kh++);

	  // Replace the current boundary pieces by an approximation.
	  // Continuity towards neighbours should be C1.
	  vector<Point> start_point, end_point;
	  if (kj > 0) { //  // Inner segment, collect start point from prev segment.
	      vector<Point> pt = all_bd_curves[ki][kj-1]->ParamCurve::point
		  (all_bd_curves[ki][kj-1]->endparam(), 1);
	      start_point.push_back(pt[0]);
	      start_point.push_back(pt[1]);
	  } else {
	      if (ki > 0) { // First segment, fetch from prev crv.
		  int ci = (int)all_bd_curves[ki-1].size() - 1;
		  Point pt = all_bd_curves[ki-1][ci]->ParamCurve::point
		      (all_bd_curves[ki-1][ci]->endparam());
		  start_point.push_back(pt);
	      }
	      int ci1 = (ki + nmb_sides - 1) % nmb_sides;
	      int ci2 = (int)all_cross_curves[ci1].size() - 1;
	      while ((ci2 > -1) && (all_cross_curves[ci1][ci2].get() == 0))
		  --ci2; // Looking for last cross curve not 0.
	      if (ci2 > -1) {
		  Point cross_pt = all_cross_curves[ci1][ci2]->ParamCurve::point
		      (all_cross_curves[ci1][ci2]->endparam());
		  if (cross_pt.length() > cross_tol) {
		      vector<Point> bnd_pt = all_bd_curves[ci1][ci2]->ParamCurve::point
			  (all_bd_curves[ci1][ci2]->endparam(), 1);
		      vector<Point> curr_pt = all_bd_curves[ki][kj]->ParamCurve::point
			  (all_bd_curves[ki][kj]->startparam(), 1);
		      double coef1, coef2;
		      CoonsPatchGen::blendcoef(bnd_pt[1].begin(), cross_pt.begin(),
						 curr_pt[1].begin(), 3, 1, &coef1, &coef2);
		      // We would like to use proj of current start point into parameter
		      // plane given by bd_curve and cross_curve of neighbouring face.
		      Point new_start_tan(bnd_pt[1]*coef1 + cross_pt*coef2);
		      new_start_tan.normalize();
		      start_point.push_back(new_start_tan);
		  }
	      } else { // Fetch local tangent.
		  vector<Point> curr_pt = all_bd_curves[ki][kj]->ParamCurve::point
		      (all_bd_curves[ki][kj]->startparam(), 1);
		  Point new_start_tan(curr_pt[1]);
		  new_start_tan.normalize();
		  //start_point[1] = new_start_tan;
		  // 	    start_point[1] = shared_ptr<Point>(new Point(new_start_tan));
	      }
	  }

	  // @@sbr We should also use points from neighbour faces.
	  // 	  else if (

	  if (kh < (int)all_bd_curves[ki].size()) { // Inner segment.
	      vector<Point> pt = all_bd_curves[ki][kh]->ParamCurve::point
		  (all_bd_curves[ki][kh]->startparam(), 1);
	      end_point.push_back(pt[0]);
	      end_point.push_back(pt[1]);
	  } else {
	      if (ki == (int)all_bd_curves.size() - 1) {
		  // Last segment, fetch from first crv.
		  Point pt = all_bd_curves[0][0]->ParamCurve::point
		      (all_bd_curves[0][0]->startparam());
		  end_point.push_back(pt);
	      } else if (curve_from_neighbour[ki+1][0]) {
		  Point pt = all_bd_curves[ki+1][0]->ParamCurve::point
		      (all_bd_curves[ki+1][0]->startparam());
		  end_point.push_back(pt);
	      }
	      int ci1 = (ki + 1) % nmb_sides;
	      int ci2 = 0;
	      while ((ci2 < (int)all_cross_curves[ci1].size()) &&
		     (all_cross_curves[ci1][ci2].get() == 0))
		  ++ci2; // Looking for last cross curve not 0.
	      if (ci2 < (int)all_cross_curves[ci1].size()) {
		  Point cross_pt = all_cross_curves[ci1][ci2]->ParamCurve::point
		      (all_cross_curves[ci1][ci2]->startparam());
		  if (cross_pt.length() > cross_tol) {
		      vector<Point> bnd_pt = all_bd_curves[ci1][ci2]->ParamCurve::point
			  (all_bd_curves[ci1][ci2]->startparam(), 1);
		      vector<Point> curr_pt = all_bd_curves[ki][kh-1]->ParamCurve::point
			  (all_bd_curves[ki][kh-1]->endparam(), 1);
		      double coef1, coef2;
		      CoonsPatchGen::blendcoef(bnd_pt[1].begin(), cross_pt.begin(),
						 curr_pt[1].begin(), 3, 1, &coef1, &coef2);
		      // We would like to use proj of current start point into parameter
		      // plane given by bd_curve and cross_curve of neighbouring face.
		      Point new_end_tan(bnd_pt[1]*coef1 + cross_pt*coef2);
		      new_end_tan.normalize();
		      end_point.push_back(new_end_tan);
		      // 	Point pt = all_cross_curves[ci1][ci2]->ParamCurve::point
		      //          (all_cross_curves[ci1][ci2]->startparam());
		      //        end_point[1] = shared_ptr<Point>(new Point(pt));
		  }
	      } else { // Fetch local tangent.
		  vector<Point> curr_pt = all_bd_curves[ki][kh-1]->ParamCurve::point
		      (all_bd_curves[ki][kh-1]->startparam(), 1);
		  Point new_end_tan(curr_pt[1]);
		  new_end_tan.normalize();
		  //end_point[1] = new_end_tan;
		  //end_point[1] = shared_ptr<Point>(new Point(new_end_tan));
	      }
	  }

	  int max_iter = 7;
	  try {
	      double max_dist;
	      vector<shared_ptr<ParamCurve> > tmp_bd_curves(all_bd_curves[ki].begin(),
							     all_bd_curves[ki].end());
	      // shared_ptr<SplineCurve> appr_cv
	      // 	  (CurveCreators::approxCurves(&(all_bd_curves[ki][0]) + kj,
	      // 					 &(all_bd_curves[ki][0]) + kh,
	      // 					 start_point, end_point,
	      // 					 approxtol_, max_dist, max_iter));
	      shared_ptr<SplineCurve> appr_cv
		  (CurveCreators::approxCurves(&(tmp_bd_curves[0]) + kj,
						 &(tmp_bd_curves[0]) + kh,
						 start_point, end_point,
						 approxtol_, max_dist, max_iter));
	      if ((max_dist < approxtol_) || require_cord_length_param) {
		  all_bd_curves[ki][kj] = appr_cv;
		  if (max_dist > approxtol_) {
		      MESSAGE("Failed approximating within tolerance (" << approxtol_ <<
				 "), using cv anyway. Dist: " << max_dist);
		  }
	      } else {
		  MESSAGE("Failed approximating, using input segment.");
		  break;
	      }
	  } catch (...) {
	      MESSAGE("Failed approximating, using input segment.");
	      break;
	  }
	  // 	approxBoundaryPiece(&all_bd_curves[ki][kj], kh-kj, 
	  // 			    start_point, end_point, approxtol_);
	  // @@sbr Essential that the approxtol_ < toptol_.neighbour.
	  // 			      toptol_.neighbour);

	  // Now only the first curve of the sequence is valid. Remove
	  // the remaining curves and corresponding cross tangent curves.
	  if (kh-kj > 1) {
	      all_bd_curves[ki].erase(all_bd_curves[ki].begin()+kj+1,
				      all_bd_curves[ki].begin()+kh);
	      // These cross curves are all dummy curves.
	      all_cross_curves[ki].erase(all_cross_curves[ki].begin()+kj+1,
					 all_cross_curves[ki].begin()+kh);
	      curve_from_neighbour[ki].erase(curve_from_neighbour[ki].begin()+kj+1, 
					     curve_from_neighbour[ki].begin()+kh);
	  }
      }

      // Join boundary curves.
      all_bd_curves[ki][0]->makeKnotStartRegular();
      for (kj=1; kj<(int)all_bd_curves[ki].size(); kj++)
	  {
              // double dummy_dist;
              // const int debug_cont = 0; // @@sbr101706 Calculate the input continuity!
	      all_bd_curves[ki][kj]->makeKnotEndRegular();
	      all_bd_curves[ki][0]->appendCurve(all_bd_curves[ki][kj].get());//, debug_cont, dummy_dist);
	  }
      bd_curves.push_back(all_bd_curves[ki][0]);

      // Make sure that the boundary curves and the cross boundary
      // curves still have matching parameter intervals.
      for (kj=1; kj<(int)all_bd_curves[ki].size(); kj++)
	  if (all_cross_curves[ki][kj].get() != 0)
	      all_cross_curves[ki][kj]->
		  setParameterInterval(all_bd_curves[ki][kj]->startparam(),
				       all_bd_curves[ki][kj]->endparam());

      // Join cross boundary curves. Create missing curves.
      if (compute_cross_curves) {
	  cross_curves.push_back
	      (joinCrossCurves(all_cross_curves[ki], 
			       all_bd_curves[ki][0]->startparam(), 
			       all_bd_curves[ki][all_bd_curves[ki].size()-1]->
			       endparam()));
      } else {
	  cross_curves.push_back(shared_ptr<SplineCurve>());
      }
  } // End for each boundary


    // Ensure C0 continuity in the corners. Try also to ensure G1
    // continuity when there is cross boundary information along one
    // of the meeting edges
    // @@sbr This should not be necessary!!!
//   adjustEndPoints(bd_curves, cross_curves, local_at_end);

  // We must make sure that cross_curves share parameter interval with boundary cvs.
  for (ki = 0; ki < (int)cross_curves.size(); ++ki)
    if (cross_curves[ki].get() != 0)
      cross_curves[ki]->setParameterInterval(bd_curves[ki]->startparam(),
					     bd_curves[ki]->endparam());

  if (bd_curves.size() == 2 || bd_curves.size() == 3)
    {
       if (bd_curves.size() != cross_curves.size()) {
	  MESSAGE("This should never happen, fix! Continuing without cross tangent curves.");
	  cross_curves.resize(bd_curves.size());
	  for (kj = 0; kj < (int)cross_curves.size(); ++kj) {
	     cross_curves[kj] = shared_ptr<SplineCurve>();
	  }
       }

      // Two or three boundary curves exist. Extend the
      // set of boundaries with degenerate curves to be
      // able to make a surface with degenerate edge(s).
      vector<shared_ptr<ParamCurve> > boundary;
      vector<shared_ptr<ParamCurve> > crossboundary;
      for (kj=0; kj<(int)bd_curves.size(); kj++)
	{
	  boundary.push_back(bd_curves[kj]);
	  crossboundary.push_back(cross_curves[kj]);
	}
      cmUtils::extendWithDegBd(corner, boundary, crossboundary, idxmin);
      status.addWarning(FT_DEGENERATE_SURFACE);

      for (kj=0; kj<(int)boundary.size(); kj++)
	if (kj >= (int)bd_curves.size() || 
	    boundary[kj].get() != bd_curves[kj].get())
	  {
	    bd_curves.insert(bd_curves.begin()+kj, 1, 
			     dynamic_pointer_cast<SplineCurve,
			     ParamCurve>(boundary[kj]));
	    cross_curves.insert(cross_curves.begin()+kj, 1, dummycrv);
	  }
    }
  

  return status;
}

//===========================================================================
shared_ptr<SplineCurve>  
ftSurfaceSet::getBoundaryPiece(Point& ptmin, Point& ptmax, double eps)
//---------------------------------------------------------------------------
//
// Purpose: Return a part of the boundary of the current surface if the
//          boundary is stored separately.
//
//===========================================================================
{ 
  double knot_diff_tol = 1e-05;
  double mindist = 100.0*eps;
  int minind = -1;
  int ki;
  shared_ptr<SplineCurve> crvpiece;  // Empty curve so far
  SplineCurve *currcv;
  Point clo_pt;
  double t1, t2, tmin, tmax, cdist1, cdist2;
  bool turncurve;

  // Traverse all boundary edges, looking for the right one.
  for (ki=0; ki<(int)approx_bd_curves_.size(); ki++)
    {
      currcv = approx_bd_curves_[ki].get();
      if (currcv->isDegenerate(toptol_.gap))
	{
	  // We are not interested in this curve.
	  cdist1 = cdist2 = mindist;
	}
      else
	{
	  currcv->closestPoint(ptmin, currcv->startparam(), currcv->endparam(),
			       t1, clo_pt, cdist1);
	  currcv->closestPoint(ptmax, currcv->startparam(), currcv->endparam(),
			       t2, clo_pt, cdist2);
	}
      if (cdist1 + cdist2 < mindist)
	{
	  mindist = cdist1 + cdist2;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	  tmin = min(t1,t2);
	  tmax = max(t1,t2);
#else
	  tmin = std::min(t1,t2);
	  tmax = std::max(t1,t2);
#endif
	  turncurve = (t1 > t2) ? true : false;
	  minind = ki;
	}
    }

  if (mindist > eps)
    return crvpiece;   // No relevant surface boundary

  // Pick the relevant part of the surface boundary
  currcv = approx_bd_curves_[minind].get();

//   // We do not want to insert knots too close to existing knots.
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
  try {
      crvpiece = shared_ptr<SplineCurve>
	  (dynamic_cast<SplineCurve*>(currcv->subCurve(tmin, tmax, knot_diff_tol)));
  } catch (...) {
      THROW("Failed extrating subcurve!");
  }
#else
//   crvpiece = shared_ptr<SplineCurve>(currcv->subCurve(tmin, tmax));
  try {
      crvpiece = shared_ptr<SplineCurve>(currcv->subCurve(tmin, tmax, knot_diff_tol));
  } catch (...) {
      THROW("Failed extrating subcurve!");
  }
#endif
  
  if (turncurve)
    crvpiece->reverseParameterDirection();
  return crvpiece;
}

//===========================================================================
ftMessage
ftSurfaceSet::getNextEdgePieceInfo(std::vector<ftEdgeBase*>::iterator& first_edge,
				   std::vector<ftEdgeBase*>::iterator& last_edge,
				   ftFaceBase*& adjsurf, SplineCurve*& cv,
				   Point& ptmin, Point& ptmax,
// 				     double& tmin, double& tmax, 
				   double& length)
//-----------------------------------------------------------------------
// As long as the boundary conditions are uniform (no neighbouring 
// surface or the same neighbouring surface), compute the size 
// of the current type of condition, and return the condition.
//===========================================================================
{
  double tmin, tmax;
  double knot_diff_tol = 1e-05;
  ftMessage status;
  ftFaceBase *firstadj = 0;
  adjsurf = 0;
//   cv = 0;
  length = 0.0;
  ftSSfEdge* s_edge;
  ftEdge* g_edge;
  shared_ptr<ParamCurve> edge_cv;

  while (first_edge != last_edge)
    {
      s_edge = dynamic_cast<ftSSfEdge*>(first_edge[0]);      
      if (s_edge == 0)
	{
	  status.setError(FT_UNEXPECTED_DATA_TYPE);
	  return status;
	}

      g_edge = s_edge->geomEdge();
      if (g_edge == 0)
	{
	  status.setError(FT_UNEXPECTED_DATA_TYPE);
	  return status;
	}

      ftEdgeBase *other_edge = s_edge->twin();
      if (other_edge)
	{
	  // A neighbouring face exist. Check if the surface is 
	  // created.
	  adjsurf = other_edge->face();
	  if (adjsurf->surface().get() == 0)
	    adjsurf = 0;   // The surface is not created.
	}
      else
	adjsurf = 0;
		  
      if (edge_cv && adjsurf)
	// A neighboruing face exist, but did not previously. Pick
	// the edge so far before picking from the neighbouring face.
	break;
      else if (adjsurf && firstadj == 0)
	{
	  // A neighbouring face exist, and this is the first face
	  firstadj = adjsurf;
	  ptmin = s_edge->point(s_edge->tMin());
	}
      else if (adjsurf != firstadj)
	// A new neighbouring face along the current boundary of the
	// current surface
	break;
      else if (edge_cv && !adjsurf)
	{
	  // No neighbouring surface exist. Fetch the boundary
	  // curve from the current edge, and check if we are 
	  // still along the same geometry curve as for the
	  // previous edge.
	    if (g_edge->geomCurve() != edge_cv)
		break;
	    tmax = g_edge->tMax();
	}
      else if (adjsurf && firstadj && adjsurf == firstadj);
      else
	{
	  // Fetch the boundary curve from the local edge.
	  // No previous geometry curve exist.
	  edge_cv = g_edge->geomCurve();
	  tmin = g_edge->tMin();
	  tmax = g_edge->tMax();
	}

      ptmax = s_edge->point(s_edge->tMax());
      length += g_edge->estimatedCurveLength();

      first_edge++;
    }
  adjsurf = firstadj;

  // We extract geometry curve, this to handle CurveOnSurface.
  if (edge_cv) {
      shared_ptr<ParamCurve> sub_cv;
      try {
	  sub_cv = shared_ptr<ParamCurve>(edge_cv->subCurve(tmin, tmax, knot_diff_tol));
      } catch (...) {
	  THROW("Failed extrating subcurve!");
      }
      try {
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	  cv = (SplineCurve*)(sub_cv->geometryCurve()->clone());
#else
	  cv = sub_cv->geometryCurve()->clone();
#endif
      } catch (...) {
	  if (sub_cv->instanceType() == Class_CurveOnSurface) {
	      shared_ptr<CurveOnSurface> cv_on_sf =
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(sub_cv);
	      shared_ptr<ParamCurve> pcv = cv_on_sf->parameterCurve();
	      ASSERT(pcv.get() != 0);
	      shared_ptr<ParamSurface> srf = cv_on_sf->underlyingSurface();
	      double epsge = toptol_.neighbour;
	      try {
		  cv = CurveCreators::liftParameterCurve(pcv, srf, epsge);
	      } catch (...) {
		  THROW("Emergency method failed!");
	      }
	  } else {
	      THROW("Unexpected incident occured!");
	  }
      }
  }

  return status;
}

//===========================================================================
shared_ptr<SplineCurve> 
ftSurfaceSet::joinCrossCurves(vector< shared_ptr<SplineCurve> >& cross_curves,
			      double startpar, double endpar)
//---------------------------------------------------------------------------
//
// Purpose: Join cross boundary curves along an edge into one curve, and
//          create missing curves. If no curve exist return a dummy 
//          shared pointer.
//
//===========================================================================
{
  int ki;
  shared_ptr<SplineCurve> curr;
  SplineCurve *prev = 0;  // Pointer to previous curve
  double dist;
  Point pnt1, pnt2;

  // Traverse the cross-boundary curves and append the curves into one
  curr = cross_curves[0];
  prev = curr.get();
  for (ki=1; ki<(int)cross_curves.size(); ki++)
    {
      if (curr.get() == 0 && cross_curves[ki].get() == 0)
	continue;   // No cross tangent curve found yet.

      else if (curr.get() != 0 && prev != 0 && cross_curves[ki].get() != 0)
	{
	  // Append this curve to the current cross boundary curve
	  curr->appendCurve(cross_curves[ki].get(), -1, dist);
	  prev = cross_curves[ki].get();
	}

      else if (curr.get() != 0 && cross_curves[ki].get() != 0)
	{
	  // The previous cross boundary curve does not exist. Make a
	  // linear interpolation between the last boundary curve and
	  // the current one. Append the linear blend and this curve
	  // to the current boundary curve.
	  curr->point(pnt1, curr->endparam());
	  cross_curves[ki]->point(pnt2, cross_curves[ki]->startparam());
	  SplineCurve midcrv(pnt1, curr->endparam(), 
			       pnt2, cross_curves[ki]->startparam());

	  curr->appendCurve(&midcrv, -1, dist);
	  curr->appendCurve(cross_curves[ki].get(), -1, dist);
	  prev = cross_curves[ki].get();
	}

      else if (cross_curves[ki].get() != 0)
	{
	  // No previous curve exist. Construct the current boundary
	  // up to now and append this curve.
	  cross_curves[ki]->point(pnt2, cross_curves[ki]->startparam());
	  curr = shared_ptr<SplineCurve>
	    (new SplineCurve(pnt2, startpar, 
			       pnt2, cross_curves[ki]->startparam()));
	  curr->appendCurve(cross_curves[ki].get(), -1, dist);
	  prev = cross_curves[ki].get();
	}
      
      else if (cross_curves[ki].get() == 0)
	{
	  // This cross tangent curve does not exist.
	  prev = cross_curves[ki].get();
	}
    }

  // Check if the last piece of the cross tangent curve exist.
  if (prev == 0 && curr.get() != 0)
    {
      // Construct the last piece of the cross tangent curve.
      curr->point(pnt1, curr->endparam());
      shared_ptr<SplineCurve> lastcrv = 
	shared_ptr<SplineCurve>
	    (new SplineCurve(pnt1, curr->endparam(), 
			       pnt1, endpar));
	  curr->appendCurve(lastcrv.get(), -1, dist);
    }

  return curr;
}

  //===========================================================================
ftMessage
ftSurfaceSet::fetchSamplePoints(const vector<ftEdgeBase*>& edgeloop,
				vector<int>& corner, vector<int>& cn,
				ftPointSet& points)
//-----------------------------------------------------------------------
// Sample data points representing the shape of the current
// ftSuperSurface.
//===========================================================================
{
  ftMessage status;
  int ki;

#ifndef NDEBUG
  {
      std::ofstream debug_out("tmp/edgeloop.g2");
      for (size_t kj = 0; kj < edgeloop.size(); ++kj)
      {
          ftEdge* ft_edge = edgeloop[kj]->geomEdge();
          shared_ptr<ParamCurve> geom_cv = ft_edge->geomCurve();
          if (geom_cv->instanceType() == Class_CurveOnSurface)
          {
              shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(geom_cv);
              shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
              if (space_cv)
              {
                  space_cv->writeStandardHeader(debug_out);
                  space_cv->write(debug_out);
              }
          }
      }
  }
#endif
  
  // Fetch edges starting at corners
  vector<ftEdgeBase*> edgc;
  edgc.reserve(4);
  for (ki=0; ki<(int)corner.size()-1; ki++)
    edgc.push_back(dynamic_cast<ftSSfEdge*>(edgeloop[corner[ki]])->geomEdge());
  if (edgc.size() != 4)
    {
      status.setError(FT_NON_4_SIDED_SURF);
      return status;
    }

  // Fetch data points
  cn.resize(edgc.size());
  std::fill(cn.begin(), cn.end(), -1);
  //  points = shared_ptr<ftPointSet>(new ftPointSet());
  // faces_points[i] intended to contain PointIters to elements in the above 
  // introduced points. The indexing is the same as in faces_.
  //  vector<vector<ttlPoint*> > faces_points(faces_.size());

  // We sample the boundary points.
  //  getInitBndData(edgc, *points, faces_points, &cn[0]);
  try
    {
      getInitBndData(edgc, points, &cn[0]);
    }
  catch (...)
    {
      THROW("Something went wrong!");
    }

  if (!sampleOnlyEdges()) {
    // Maximum number of inner pointer for each patch
    int total_max = 200; //50;// 100; //1000; // Number of pts to sample in one direction.
    int max_sample = (int)(total_max/(sqrt(double(faces_.size()))));
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    max_sample = min(max_sample, 100);
    //max_sample = max(5,max_sample);
    //max_sample = max(50,max_sample);
#else
    max_sample = std::min(max_sample, 100);
    //max_sample = std::max(5,max_sample);
    //max_sample = std::max(50,max_sample);
#endif

    // @@sbr Remove when done debugging
//     max_sample = 3;

    // We sample the interior points (topology is not updated).
  try
    {
    getInitInnerData(points, max_sample);
    }
  catch (...)
    {
      THROW("Something went wrong!");
    }
  }

//   // For debugging
//   ofstream of("data/output/sampled_points.g2");
//   SplineCurve out_curve;
//   Point origo(0.0, 0.0, 0.0);
//   for (ki = 0; ki < points_->size(); ++ki) {
//     Vector3D sampled_pt = (*points_)[ki]->getPoint();
//     Point pt(sampled_pt[0], sampled_pt[1], sampled_pt[2]);
//     Point pt2 = 0.1*origo+0.9*pt;
//     out_curve = SplineCurve(pt2, pt);
//     out_curve.writeStandardHeader(of);
//     out_curve.write(of);
//   }
//   // end of debugging

  return status;
}

//===========================================================================
ftMessage
ftSurfaceSet::modifySurface(int maxiter, int max_update_iter, 
			      vector<int>& edge_derivs,
			      ftPointSet& points,
			      double& max_error, double& mean_error)
//-----------------------------------------------------------------------
// Modify the current surface. The surface may be representing the
// spline space without containing coefficient data. Parameter iteration,
// iteration on the weight on approximation compared to smoothness, and
// refinement of the spline space is included.
//===========================================================================
{
  ftMessage status;
  bool isOK;
  double prevmax, prevmean;
  double prevapprox;

  max_error = points.getMaxDist();
  mean_error = points.getMeanDist();

  // Make instance for surface modification and smoothing.
  double approx_orig_tol = -1.0; // We will not perform inner smoothing.
  ftSmoothSurf smoothsrf(surf_, approxtol_, approx_orig_tol, edge_derivs, 
			 max_update_iter);
  smoothsrf.setApproxWeight(approxweight_);

  try {
      bool reparam = true; //false; //@@sbr Should be true, fix when done debugging!!!
    isOK = smoothsrf.update(points, toptol_.gap, reparam);
  }
  catch(...)
    {
	status.addWarning(FT_COULD_NOT_SMOOTH_SURFACE);
//       status.setError(FT_ERROR_IN_SURFACE_CREATION);
	return status;
    }

//   // debugging
//   // We loop through points checking accuracy.
//   int nmb_pts = points->size();
//   double comp_dist;
//   Vector3D vec;
//   Point pt;
//   double clo_u, clo_v, clo_dist;
//   Point clo_pt;
//   double tol = 1e-10;
//   bool bnd_pt = true;
//   double seed[2];
//   int i;
//   for (i = 0; i < nmb_pts; ++i)
//     if ((*points)[i]->isOnBoundary() && ((*points)[i]->getDist() > approxtol_)) {
//       comp_dist = (*points)[i]->getDist();
//       vec = points->get3dNode(i);
//       pt = Point(vec[0], vec[1], vec[2]);
//       seed[0] = points->getU(i);
//       seed[1] = points->getV(i);
// //       surf_->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, tol, bnd_pt, seed);
//       RectDomain *rd = 0;
//       surf_->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist, rd, seed);
// //       std::cout << "Dist: " << comp_dist << ", dist: " << clo_dist
// // 		<< std::endl;
//     }
//   // end of debugging

//   // debugging
//        char file_name3[300];
//        sprintf(file_name3, "data/output/surf_%d.g2",1);
//        ofstream of1(file_name3);
//        surf_->writeStandardHeader(of1);
//        of1 << *surf_;
//   // end debugging


  // Iterate to make sure that the surface is accurate enough.
//    max_error = prevmax = points.getMaxDist();
//    mean_error = prevmean = points.getMeanDist();
  smoothsrf.getError(max_error, mean_error);
  prevmax = max_error;
  prevmean = mean_error;
  int iter=0;

  // Make a copy of the previous version of the surface.
  shared_ptr<SplineSurface> prev_surf;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
  prev_surf =
    shared_ptr<SplineSurface>(dynamic_cast<SplineSurface*>(surf_->clone()));
  ALWAYS_ERROR_IF(prev_surf.get()==NULL,"Bad cast");

#else
  prev_surf = shared_ptr<SplineSurface>(surf_->clone());
#endif

  approxweight_ = smoothsrf.getApproxWeight();
  while (!isOK && iter < maxiter) {
    iter++;

    // Parameterize the extended set of data points.
    //     param_instance.reparametrize(points, PROJECT_SURF, surf);

    // Add more degrees of freedom to the surface where the
    // error is too large.
    try {
      smoothsrf.refineSurf(points);
    }
    catch(...)
      {
	status.addWarning(FT_COULD_NOT_REFINE_SURFACE);
	break;
      }

    // Update the surface with the new points and refined data
    // size.
    try {
	bool reparam = true; //false; //@@sbr Should be true, fix when done debugging!!!
      isOK = smoothsrf.update(points, toptol_.gap, reparam);
    }
    catch(...)
      {
	status.addWarning(FT_COULD_NOT_SMOOTH_SURFACE);
	return status;
      }

    // debugging
    //        char file_name[300];
    //        sprintf(file_name, "data/output/surf_%d.g2",iter+1);
    //        ofstream of3(file_name);
    //        surf_->writeStandardHeader(of3);
    //        of3 << *surf_;
    // end debugging

    prevapprox = approxweight_;
    approxweight_ = smoothsrf.getApproxWeight();
//      max_error = points.getMaxDist();
//      mean_error = points.getMeanDist();
  smoothsrf.getError(max_error, mean_error);
    if (!(max_error < 0.95*prevmax || mean_error < 0.95*prevmean))
      break;
    prevmax = max_error;
    prevmean = mean_error;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    prev_surf  = shared_ptr<SplineSurface>
      (dynamic_cast<SplineSurface*>(surf_->clone()));
    ALWAYS_ERROR_IF(prev_surf.get()==NULL, "Bad cast");

#else
    prev_surf  = shared_ptr<SplineSurface>(surf_->clone());
#endif
  }

  if (max_error > prevmax)
    {
      surf_ = prev_surf;  // Revert to the previous version
      approxweight_ = prevapprox;
    }
  
  //   // debugging
  //   std::ofstream fileout("data/output/super_surf.g2");
  //   surf_->writeStandardHeader(fileout);
  //   surf_->write(fileout);
  //   // end of debugging

  if (!isOK)
    status.addWarning(FT_APPROX_ERROR_TOO_LARGE);

  return status;
}


//===========================================================================
ftMessage ftSurfaceSet::getOuterLoop(vector<ftEdgeBase*>& outer_loop,
				     vector<int>& corners,
				     double bend_tol, double kink_tol)
//===========================================================================
{
    ftMessage status;
    corners.clear();

    // First fetch boundary information (number of boundary loops,
    // number of curves for the outer loop).
    vector< vector<ftEdgeBase*> > loopvec;
    try {
	BoundaryLoops(loopvec);
    } catch (...) {
	status.setError(FT_NO_DATA); // It appears first_edges_ was not set.
	return status;
    }

    //    if (loopvec.size() > 0)
    //    {
    //          std::ofstream dump("debugdump.g2");
    //        for (int ki=0; ki<1; ki++)
    //          for (int kj=0; kj<loopvec[ki].size(); kj++)
    //    	{
    //    	  ftEdge *geomedge = loopvec[ki][kj]->geomEdge();
    //    	  geomedge->SpaceCurve()->writeStandardHeader(dump);
    //    	  geomedge->SpaceCurve()->write(dump);
    //    	}
    //    }
    if (loopvec.size() > 1) // We do not handle "trimmed" sets.
	{
	    status.setError(FT_WRONG_NO_OF_BOUNDARIES);
	    return status;
	}

    if (loopvec.size() == 0)
	{
	    status.setError(FT_NO_DATA);
	    return status;
	}

    if (loopvec[0].size() == 1)
	{
	    status.setError(FT_WRONG_NO_OF_BOUNDARIES);
	    return status;
	}

#ifdef FANTASTIC_DEBUG
    std::ofstream debug("data/debug.g2");
    for (size_t ki = 0; ki < loopvec[0].size(); ++ki) {
	ftEdge* edge1 = loopvec[0][ki]->geomEdge();
	double t0 = edge1->tMin();
	double t1 = edge1->tMax();
	double knot_diff_tol = 1e-05;
	try {
	    shared_ptr<ParamCurve> scp(edge1->geomCurve()->subCurve(t0, t1, knot_diff_tol));
	    if (scp->instanceType() == Class_SplineCurve) {
		shared_ptr<SplineCurve> sc1 =
		    dynamic_pointer_cast<SplineCurve, ParamCurve>(scp);
		sc1->writeStandardHeader(debug);
		sc1->write(debug);
	    } else if (scp->instanceType() == Class_CurveOnSurface) {
		shared_ptr<CurveOnSurface> cv_on_sf =
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(scp);
		if (cv_on_sf->spaceCurve().get() != 0) {
		    cv_on_sf->spaceCurve()->writeStandardHeader(debug);
		    cv_on_sf->spaceCurve()->write(debug);
		}
	    }
	} catch (...) {
	    MESSAGE("Failed extrating subcurve!");
	}
    }
#endif // FANTASTIC_DEBUG
  
    int ki, kj;
    int nmbedge = (int)loopvec[0].size();
    for (ki=nmbedge-1, kj=0; kj<nmbedge; ki=(ki+1)%nmbedge, kj++)
	{
	    tpJointType cont;
	    //    std::cout << "ki=" << ki << ", kj=" << kj << std::endl;
	    cont = loopvec[0][ki]->checkContinuity(loopvec[0][kj],
						   toptol_.neighbour, toptol_.gap,
						   bend_tol, kink_tol);
	    /* if (cont >= JOINT_GAP)
	       {
	       message = shared_ptr<ftMessage>(new ftError(FT_DISCONTINUITY));
	       return message;
	       }
	       else */ if (cont >= JOINT_G0) {
		   corners.push_back(kj);
#ifdef FANTASTIC_DEBUG
		   std::ofstream debug("data/debug.g2");
		   ftEdge* edge1 = loopvec[0][ki]->geomEdge();
		   double knot_diff_tol = 1e-05;
		   try {
		       shared_ptr<SplineCurve> sc1
			   (dynamic_cast<SplineCurve*>
			    (edge1->geomCurve()->subCurve(edge1->tMin(), edge1->tMax(), knot_diff_tol)));
		       ftEdge* edge2 = loopvec[0][kj]->geomEdge();
		       shared_ptr<SplineCurve> sc2
			   (dynamic_cast<SplineCurve*>
			    (edge2->geomCurve()->subCurve(edge2->tMin(), edge2->tMax(), knot_diff_tol)));
		       if (sc1.get() != 0) {
			   sc1->writeStandardHeader(debug);
			   sc1->write(debug);
		       }
		       if (sc2.get() != 0) {
			   sc2->writeStandardHeader(debug);
			   sc2->write(debug);
		       }
		   } catch (...) {
		       MESSAGE("Failed extrating subcurve!");
		   }
#endif // FANTASTIC_DEBUG
	       }
	}

    // We then run through additional_corner_pts_.
    // Assuming edges have already been split.
    cmUtils::updateWithNewCorners(loopvec[0], corners, additional_corner_pts_, toptol_.gap);

    // If the number of boundaries was not correct we see if we may still consider the surface set as
    // having a rectangular domain.
    if (corners.size() > 4)
    {
        vector<int> new_corners = cmUtils::removeInnerCorners(loopvec[0], corners);
        if (new_corners.size() == 4)
        {
            MESSAGE("DEBUG: Successfully reduced the number of corners to 4!");
            corners = new_corners;
        }
    }

    if (corners.size() < 2 || corners.size() > 4)
    {
        status.setError(FT_WRONG_NO_OF_BOUNDARIES);
        return status;
    }

    int one_past_last = (corners.front() > 0) ? corners[0] : corners[0] + (int)loopvec[0].size();
    corners.push_back(one_past_last);
    outer_loop = loopvec[0];

    return status;
}

//===========================================================================
void ftSurfaceSet::getInitBndData(vector<ftEdgeBase*>& edgc, ftPointSet& points,
				  int cn[])
//---------------------------------------------------------------------------
//
// Purpose: Get initial data points from the input surface set.
//
//===========================================================================
{
  int bnd = -1;    // Boundary type of current edge. 1 == outer bnd, 2 == inner bnd.
  int csidx;  // Current_edge surface-index (in faces_).

  // Set up vectors for handling points which are not properly
  // connected to its neighbour.
  vector<ftEdgeBase*> free_edges;
  vector<bool> is_start;
  vector<PointIter> vertex_points; // We store every vertex point.

  ftEdgeBase *curr_edge = 0;
  // For ease of algorithm, we start by picking an outer boundary edge.
  for (csidx = 0; csidx < (int)faces_.size(); ++csidx) {
      ftEdgeBase* start_edge = faces_[csidx]->startEdges()[0].get();
      curr_edge = start_edge;
      if (curr_edge->onBoundary()) {
	  bnd = 1;
	  break;
      }
      curr_edge = start_edge->next();
      while ((!(curr_edge->onBoundary())) && (curr_edge != start_edge))
	  curr_edge = curr_edge->next();
      if (curr_edge->onBoundary()) {
	  bnd = 1;
	  break;
      }
  }
  ASSERT(bnd != -1);
  // Fetch all edges meeting in this crossing (start of curr_edge), not counting edg1
  vector<ftEdgeBase*> edg1;
  vector<bool> start1;
  curr_edge->adjacentEdges(true, edg1, start1);

  // Evaluate the start edge in the start point.
  double parmin = curr_edge->tMin();
  double parmax = curr_edge->tMax();
  Point pnt = curr_edge->point(parmin);
  Vector3D pnt3D(pnt[0], pnt[1], pnt[2]);
  shared_ptr<ftSurfaceSetPoint> ftpnt(new ftSurfaceSetPoint(pnt3D, bnd));
  double* seed = NULL;
  shared_ptr<Point> par_pt;
  double tpar = parmin;
  getParPoint(curr_edge, tpar, par_pt);
  if (par_pt.get() != 0) {
      seed = (*par_pt).begin();
  }
  addFaceRef(*ftpnt, faces_[csidx], seed); // edg1 not returned by AdjacentEdges().
  int ki, kj;
  for (ki = 0; ki < (int)edg1.size(); ++ki)
      if (start1[ki])
        {
	  int fi = 0; // faces_-index.
	  while (edg1[ki]->face() != faces_[fi].get())
	      ++fi;
	  par_pt.reset();
	  tpar = edg1[ki]->tMin();
	  getParPoint(edg1[ki], tpar, par_pt);
	  if (par_pt.get() != 0) {
	      seed = (*par_pt).begin();
	  } else {
	      seed = NULL;
	  }
	  addFaceRef(*ftpnt, faces_[fi], seed);
	}

  // Store the point in the pointset
  points.addEntry(ftpnt);
  PointIter prevpt = points.lastAdded();
  points.setFirst(prevpt);
  bool set_second = true;

  // Gather information about possible corner point.
  int cidx = points.size() - 1;
  int kcn;
  for (kcn=0; kcn<(int)edgc.size(); ++kcn)
      if (edgc[kcn] == curr_edge)
	  cn[kcn] = cidx;
    
  // Remember all edges meeting in this corner point.
  for (ki=0; ki<(int)edg1.size(); ki++)
    {
      free_edges.push_back(edg1[ki]);
      if (curr_edge->twin() == edg1[ki])
	is_start.push_back(false); // To avoid duplicate of point.
      else
	is_start.push_back(start1[ki]);
      vertex_points.push_back(prevpt);
    }

  // Traverse all edges
  while (true)
    {
	// Get number of points to evaluate along the edge (including end points).
        ftEdge* ft_edge = dynamic_cast<ftEdge*>(curr_edge);
	int nmb_eval = nmbToEval(ft_edge,
                                 curr_edge->tMin(), curr_edge->tMax()); // >= 2
        //MESSAGE("nmb_eval: " << nmb_eval);
        double cv_length = ft_edge->estimatedCurveLength();
        //MESSAGE("DEBUG: cv_length: " << cv_length);
	const int min_samples = 5; // With too few samples we may fail picking the correct sub-face.
	const int max_samples = 80;
        nmb_eval = std::max(min_samples, std::min(nmb_eval, max_samples));
	getEdgeInnerData(curr_edge, prevpt, points, edgc, cn, set_second, nmb_eval);
	if (nmb_eval != 2)
	    set_second = false;

	// Check if the last point of curr_edge exist already
	for (kj=0; kj<(int)free_edges.size(); kj++)
	  if (free_edges[kj] == curr_edge && !is_start[kj])
	    break;

	if (kj < (int)free_edges.size())
	  {
            // As a bnd point is to have 2 bnd neighbours only, we may have to
	    // sample mid-point on edge.
	      if ((nmb_eval == 2) && (curr_edge->twin() != 0) &&
		  (prevpt->isOnBoundary()) && (vertex_points[kj]->isOnBoundary()))
		  getEdgeInnerData(curr_edge, prevpt, points, edgc, cn,
				   set_second, 3);

	    // Connect to the already existing corner point
            prevpt->addNeighbour(vertex_points[kj]);
            vertex_points[kj]->addNeighbour(prevpt);

	    // Remove free_edge instance
	    free_edges.erase(free_edges.begin()+kj, free_edges.begin()+kj+1);
	    is_start.erase(is_start.begin()+kj, is_start.begin()+kj+1);
	    vertex_points.erase(vertex_points.begin()+kj,
				vertex_points.begin()+kj+1);
	    // To avoid duplicate of points we run through free_edges to see
	    // if curr_edge.twin() is stored as a start edge.
	    int ki2;
	    for (ki2 = 0; ki2 < (int)free_edges.size(); ++ki2)
		if (free_edges[ki2]->twin() == curr_edge) {
		    is_start[ki2] = false;
		    //	    break; // We're done
		}
	  }
	else
	  {
	    // Check if the point lies at an outer boundary.
	    bnd = curr_edge->onBoundary() ? 1 : 2;

	    // Fetch edges meeting in this corner point (end of curr_edge).
	    curr_edge->adjacentEdges(false, edg1, start1);
	    for (ki=0; ki<(int)edg1.size(); ki++)
	      if (edg1[ki]->onBoundary())
		  bnd = 1;

	    // Evaluate and store the endpoint.
	    pnt = curr_edge->point(parmax);

	    pnt3D = Vector3D(pnt[0], pnt[1], pnt[2]);
	    ftpnt = shared_ptr<ftSurfaceSetPoint>(new ftSurfaceSetPoint(pnt3D, bnd));
	    // Store reference in faces_points.
	    //	    addFaceRef(*ftpnt, faces_[csidx].get());
	    for (ki = 0; ki < (int)edg1.size(); ++ki)
		if (start1[ki])
		    {
			int fi = 0; // faces_ indeks		
			while (edg1[ki]->face() != faces_[fi].get())
			    ++fi;
			par_pt.reset();
			tpar = edg1[ki]->tMin();
			getParPoint(edg1[ki], tpar, par_pt);
			if (par_pt.get() != 0) {
			    seed = (*par_pt).begin();
			} else {
			    seed = NULL;
			}
			addFaceRef(*ftpnt, faces_[fi], seed);
		    }

            // As a bnd point is to have 2 bnd neighbours only, we may have to
	    // sample mid-point on edge.
	    if ((nmb_eval == 2) && (curr_edge->twin() != 0) &&
		(prevpt->isOnBoundary()) && (ftpnt->isOnBoundary()))
		getEdgeInnerData(curr_edge, prevpt, points, edgc, cn, set_second, 3);

	    addBndPoint(ftpnt, prevpt, points);
	    if (set_second) {
		points.setSecond(prevpt);
		set_second = false;
	    }

	    // Remember all edges meeting in this corner point.
	    for (ki=0; ki<(int)edg1.size(); ki++)
	      {
		free_edges.push_back(edg1[ki]);
		if (curr_edge->twin() == edg1[ki])
		  is_start.push_back(false); // To avoid duplicate of point.
		else
		  is_start.push_back(start1[ki]);
		vertex_points.push_back(prevpt);
	      }
	  }

	// Get the next edge to traverse
	for (kj=0; kj<(int)free_edges.size(); kj++)
	  if (is_start[kj])
	    break;

	if (kj == (int)free_edges.size())
	  break;

	curr_edge = free_edges[kj];
	prevpt = vertex_points[kj];
	parmin = curr_edge->tMin();
	parmax = curr_edge->tMax();
      
	// Gather information about possible corner point.
	for (kcn=0; kcn<(int)edgc.size(); ++kcn)
	  if (edgc[kcn] == curr_edge) 
	    cn[kcn] = prevpt->getIndex();

	// Remove free_edge instance
	free_edges.erase(free_edges.begin()+kj, free_edges.begin()+kj+1);
	is_start.erase(is_start.begin()+kj, is_start.begin()+kj+1);
	vertex_points.erase(vertex_points.begin()+kj,
			    vertex_points.begin()+kj+1);

    }
}

//===========================================================================
void ftSurfaceSet::addBndPoint(shared_ptr<ftSurfaceSetPoint>& ftpnt,
			       PointIter& prevpt, ftPointSet& points)
//===========================================================================
{
    PointIter latestpt = points.addEntry(ftpnt);
#ifndef NDEBUG
    { // Many of these errors seem to be consecutive in order.
        const double dist = prevpt->pntDist(latestpt);
        const double num_tol = 1.0e-14;
        if (dist < num_tol)
        {
            MESSAGE("Something wrong going on, adding the same point! " << ftpnt->getPoint());
        }
    }
#endif
    prevpt->addNeighbour(latestpt);
    latestpt->addNeighbour(prevpt);
    prevpt = latestpt;
}

//===========================================================================
void ftSurfaceSet::addFaceRef(ftSurfaceSetPoint& point, shared_ptr<ftFaceBase>& face,
			      double* seed)
//===========================================================================
{
    double clo_u, clo_v;
    Point clo_pt;
    double clo_dist;
    double epsilon = 1e-6; // @@@ Change value? Initilalize from what?
    
    ParamSurface* curr_surf =
	dynamic_cast<ParamSurface*>(face->surface().get());
    // @@ Weirdness follows -- got dimension 6 until we split statement
    Vector3D p = point.getPoint();
    //    cout << "p = " << p << "ptrs: " << p.begin() << "   " << p.end() << endl;
    Point pt(p.begin(), p.end(), false);
    //    Point pt(point.getPoint().begin(), point.getPoint().end(), false);
    RectDomain rect_domain = curr_surf->containingDomain();
    if (point.isOnSubSurfaceBoundary()) {
	curr_surf->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist,
					epsilon, &rect_domain, seed);
    } else {
	curr_surf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist,
				epsilon, &rect_domain, seed);
    }

#ifdef FANTASTIC_DEBUG
    if (clo_dist > approxtol_*10) {
	std::cout << "closest dist in subsurface: " << clo_dist << std::endl;
    }
#endif // FANTASTIC_DEBUG

    Vector2D par_pt(clo_u, clo_v);
    point.addPair(face, par_pt);
}

//===========================================================================
void ftSurfaceSet::getEdgeInnerData(ftEdgeBase* curr_edge, PointIter& prevpt,
				    ftPointSet& points, vector<ftEdgeBase*>& edgc,
				    int cn[], bool set_second, int nmb_eval)
//---------------------------------------------------------------------------
//
// Purpose: Sample points from interior of edge. Called by getInitBndData().
//
//===========================================================================
{
  int bnd = curr_edge->onBoundary() ? 1 : 2;
  // Get index of curr_edge surface and curr_edge->Twin() surface
  int csidx = 0;
  const double num_tol = 1.0e-14;
  while ((curr_edge->face()) != faces_[csidx].get()) ++csidx;
  int ctidx = -1; // To denote that curr_edge->Twin() == 0
  if (!(curr_edge->onBoundary())) {
      ctidx = 0;
      while ((curr_edge->twin()->face()) != faces_[ctidx].get())
	  ++ctidx;
  }

  double parmin = curr_edge->tMin();
  double parmax = curr_edge->tMax();
  double parinc = (parmax - parmin)/(nmb_eval-1);

  // Initiate for evaluate of start point on edge
  double tpar=parmin;
  int kj=0, kcn;
  for (kcn=0; kcn<(int)edgc.size(); kcn++)
    if (edgc[kcn] == curr_edge)
      break;

  double tpar2;
  double tol = 1.0e-5*(parmax-parmin);
  shared_ptr<ParamCurve> curr_curve = curr_edge->geomEdge()->geomCurve();
  if (!(kcn+1 < (int)edgc.size() && edgc[kcn+1] == edgc[kcn])) {
      // The start point is to be evaluated only at degenerate 
      // outer edges. Update parameter to avoid evaluating in
      // the start point in other cases.
      kj++;
      tpar2 = curr_curve->nextSegmentVal(tpar, true, tol);
      //tpar = std::min(tpar + parinc, tpar2);
      tpar += parinc;
    }

  //  for (; kj<nmb_eval-1; kj++, tpar+=parinc) {
  while (tpar < parmax - num_tol) // We subtract num_tol to avoid sampling in end point.
    {
      if (kj == 0)
	// A degenerate edge is expected. Update corner incides.
	cn[kcn+1] = prevpt->getIndex() + 1;

      // Evaluate and store the point in the pointset
      Point pnt = curr_edge->point(tpar);
      Vector3D pnt3D(pnt[0], pnt[1], pnt[2]);
      //  ftSamplePoint ftpnt(pnt3D, bnd);
      shared_ptr<ftSurfaceSetPoint> ftpnt(new ftSurfaceSetPoint(pnt3D, bnd));
      // Store reference in faces_points.
      shared_ptr<Point> par_pt;
      par_pt.reset();
      getParPoint(curr_edge, tpar, par_pt);
      double* seed = NULL;
      if (par_pt.get() != 0) {
	  seed = (*par_pt).begin();
      }
      addFaceRef(*ftpnt, faces_[csidx], seed);
      if (ctidx != -1) {
	  par_pt.reset();
// 	  tpar = curr_edge->twin()->tMin();
	  double twin_dec = (curr_edge->twin()->tMax() - curr_edge->twin()->tMin())/(nmb_eval - 1);
	  double twin_tpar = curr_edge->twin()->tMax() - kj*twin_dec;
	  getParPoint(curr_edge->twin(), twin_tpar, par_pt);
	  if (par_pt.get() != 0) {
	      seed = (*par_pt).begin();
	  } else {
	      seed = NULL;
	  }
	  addFaceRef(*ftpnt, faces_[ctidx], seed);
      }
      // Add sampled point to points.
      addBndPoint(ftpnt, prevpt, points);

      if (set_second) {
	points.setSecond(prevpt);
	set_second = false;
      }
      kj++;
      //tpar = std::min(tpar+parinc, curr_curve->nextSegmentVal(tpar, true, tol));
      tpar += parinc;
    }
}


//===========================================================================
ftMessage ftSurfaceSet::getInitInnerData(ftPointSet& points, int max_sample)
//---------------------------------------------------------------------------
//
//Purpose: Get initial data points from the the interior of a set of surfaces.
//
//===========================================================================
{
  // If surf is degenerate, a uniform sampling approach will result in a dense cloud
  // of points near degenerate edge. As this is far from ideal for the
  // parametrization method, we choose another approach: if the shortest edge is
  // in the u-direction, we sample along u-iso-lines, choosing a new number for each
  // line (numbers may be equal).

  ftMessage status;

  // The topological structure of points is set outside this fuction, afterwards.
  //    double step = 20.0*toptol_.neighbour;
  double step = toptol_.neighbour; //5.0*toptol_.neighbour;
//   int csidx = 0;  // Current surface-index (as given in faces_)

  int i, j, m;
  for (i = 0; i < (int)faces_.size(); ++i) {   

      // Index of current surface.
//       csidx = faces_[i]->getId();

      shared_ptr<ParamSurface> curr_surface = faces_[i]->surface();
      shared_ptr<SplineSurface> cont_surf(createContainingSurface(curr_surface.get()));
#ifdef DEBUG
      std::ofstream of("containing_sf.g2");
      cont_surf->writeStandardHeader(of);
      cont_surf->write(of);
#endif

      //     shared_ptr<SplineSurface> under_surf; // The underlying surface.
      //     if (faces_[i]->surface()->instanceType() == Class_SplineSurface)
      //       under_surf = std::dynamic_pointer_cast<SplineSurface, ParamSurface>
      // 	(faces_[i]->surface());
      //     else if (faces_[i]->surface()->instanceType() == Class_BoundedSurface)
      //       under_surf = std::dynamic_pointer_cast<SplineSurface, ParamSurface>
      // 	((std::dynamic_pointer_cast<BoundedSurface, ParamSurface>
      // 	  (faces_[i]->surface()))->underlyingSurface());
      //     else {
      //       status.setError(FT_UNEXPECTED_DATA_TYPE);
      //       return status;
      //     }

      vector<double> lengths; // vmin, vmax, umin, umax.
      try {
	  lengths = estimateSurfSides(cont_surf.get()); //under_surf.get());
      } catch (...) {
	  status.setError(FT_NO_DATA); // Should never happen!
	  return status;
      }

      double par_step_u =
	step*(cont_surf->endparam_u() - cont_surf->startparam_u())/min(lengths[0], lengths[1]);
      double par_step_v =
	step*(cont_surf->endparam_v() - cont_surf->startparam_v())/min(lengths[2], lengths[3]);
      double par_step = 0.5*min(par_step_u, par_step_v);

      bool sample_along_u_lines = // true = sample points along u-iso-crvs (dir2).
	  min(lengths[2], lengths[3]) < min(lengths[0], lengths[1]);

      double length1 = (sample_along_u_lines) ? 0.5*(lengths[0] + lengths[1]) :
	  0.5*(lengths[2] + lengths[3]);

      int min_samples = 3; //1;

      // @@@220302 max_samples should be set based on geometry. The sampling
      // really should make sure the triangles are quite similar, as it seems
      // to be a requirement of the parametrization.
      //int max_samples = 2;
      //int max_samples = 40;
      int nmb_eval_1, nmb_eval_2;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
      nmb_eval_1 = max(min_samples, min(max_sample, (int)(length1/step)));
#else
      nmb_eval_1 =
	  std::max(min_samples, std::min(max_sample, (int)(length1/step)));
#endif
      // Containingdomain of a trimmed surface is that of the underlying surface.
      const Domain& par_dom = curr_surface->parameterDomain();
      RectDomain rect_domain = cont_surf->containingDomain(); //curr_surface->containingDomain();
      RectDomain tmp_domain = curr_surface->containingDomain();
      double tmin_1 = sample_along_u_lines ? rect_domain.umin() : rect_domain.vmin();
      double tmax_1 = sample_along_u_lines ? rect_domain.umax() : rect_domain.vmax();
      double tmin_2 = sample_along_u_lines ? rect_domain.vmin() : rect_domain.umin();
      double tmax_2 = sample_along_u_lines ? rect_domain.vmax() : rect_domain.umax();
      int dir =  sample_along_u_lines ? 0 : 1;

      double parinc_1 = (tmax_1 - tmin_1)/(nmb_eval_1+1);

      Point pnt;
      Vector3D pnt3D;
      int at_edge = 0;
//       double tol = 1e-10;
      double par1, length2, parinc_2, par2, par_u, par_v;

      vector<Vector2D> bd_train = getBdTrain(points, faces_[i]);

#ifdef FANTASTIC_DEBUG
      std::ofstream debug("data/debug.g2");
      int ki;
      int nmb_segs = (int)bd_train.size();
      vector<double> pts(nmb_segs*6);
      for (ki = 0; ki < nmb_segs; ++ki) {
	  int fi = ki;
	  int si = (ki + 1)%nmb_segs;
	  Vector3D from(bd_train[fi][0], bd_train[fi][1], 0.0);
	  Vector3D to(bd_train[si][0], bd_train[si][1], 0.0);
	  copy(from.begin(), from.end(), pts.begin() + 6*ki);
	  copy(to.begin(), to.end(), pts.begin() + 6*ki + 3);
      }
      LineCloud lc(pts.begin(), nmb_segs);
      lc.writeStandardHeader(debug);
      lc.write(debug);
#endif // FANTASTIC_DEBUG


      double tol1 = 1.0e-5*(tmax_1 - tmin_1);
      double tol2 = 1.0e-5*(tmax_2 - tmin_2);
      // par1 = std::min(tmin_1 + parinc_1, 
      // 		      cont_surf->nextSegmentVal(dir, tmin_1, true, tol1));
      par1 = tmin_1 + parinc_1;
      // for (j = 1; j < nmb_eval_1 + 1; ++j) {
      // 	  par1 = tmin_1 + j*parinc_1;
      while (par1 < tmax_1)
	{
	  shared_ptr<SplineCurve> iso_curve(cont_surf->constParamCurve //under_surf->constParamCurve
					    (par1, !sample_along_u_lines));
	  length2 = iso_curve->ParamCurve::estimatedCurveLength();
	  length2 = (sample_along_u_lines) ? 0.5*(lengths[2] + lengths[3]) :
	      0.5*(lengths[0] + lengths[1]); // @@sbr We may sample uniformly after all...
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	  nmb_eval_2 = max(min_samples, min(max_sample, (int)(length2/step)));
#else
	  nmb_eval_2= 
	      std::max(min_samples, std::min(max_sample, (int)(length2/step)));
#endif
	  parinc_2 = (tmax_2 - tmin_2)/(nmb_eval_2+1);
	  // par2 = std::min(tmin_2 + parinc_2, 
	  // 		  cont_surf->nextSegmentVal(1-dir, tmin_2, true, tol2));
	  par2 = tmin_2 + parinc_2;
	  // for (m = 1; m < nmb_eval_2 + 1; ++m) {
	  //     par2 = tmin_2 + m*parinc_2;
	  while (par2 < tmax_2)
	    {
	      par_u = sample_along_u_lines ? par1 : par2;
	      par_v = sample_along_u_lines ? par2 : par1;
	      Vector2D par_pt(par_u, par_v);
	      // 	if (!((curr_surface->parameterDomain()).isInDomain(par_pt, tol)))
	      // 	  continue; // We do not sample points outside a trimmed surface.
	      // @@sbr This may lead to sampling of pts not on actual sf (defined by containing domain).
	      // But as they stay inside bd_train it shouldn't really mather.
	      // We must futhermore check that we're inside chain defined by bd_pts.
	      if (cmUtils::insideBdTrain(par_pt, bd_train)) {

#ifdef FANTASTIC_DEBUG
		  std::ofstream debug("data/debug.g2");
		  Point par_space_pt(par_pt[0], par_pt[1], 0.0);
		  vector<double> pts(3);
		  copy(par_space_pt.begin(), par_space_pt.end(), pts.begin());
		  PointCloud3D pc(pts.begin(), 1);
		  pc.writeStandardHeader(debug);
		  pc.write(debug);
#endif // FANTASTIC_DEBUG

		  double tol = 1e-06;
		  bool is_in_dom;
		  try {
		      is_in_dom = par_dom.isInDomain(par_pt, tol);
		      if (is_in_dom) {
			  Vector2D clo_bd_pt;
			  par_dom.closestOnBoundary(par_pt, clo_bd_pt, tol);
			  // @@ Step seldom refers to empirical sampling
			  //if (par_pt.dist(clo_bd_pt) > par_step) {
			      curr_surface->point(pnt, par_u, par_v);
			      pnt3D = Vector3D(pnt[0], pnt[1], pnt[2]);
			      shared_ptr<ftSurfaceSetPoint>
				  ftpnt(new ftSurfaceSetPoint(pnt3D, at_edge,
// 						  faces_[csidx], par_pt));
							      faces_[i], par_pt));
			      points.addEntry(ftpnt);
// 		      } else {
// 			  MESSAGE("Avoided using sample pt close to bd.");
			  }
// 		  } else {
// 		      MESSAGE("Avoided using sample pt outside domain.");
		      //}
		  } catch (...) {
		      MESSAGE("Failed inDomain-test, moving on to next inner pt.");
		  }
	      }
	      // par2 = std::min(par2+parinc_2, 
	      // 		      cont_surf->nextSegmentVal(1-dir, par2, true, tol2));
	      par2 += parinc_2;
	    }
	  // par1 = std::min(par1+parinc_1, 
	  // 		  cont_surf->nextSegmentVal(dir, par1, true, tol1));
	  par1 += parinc_1;
	}

      // Add points on inner boundaries
      vector<shared_ptr<ftEdgeBase> > start_edges = faces_[i]->startEdges();
      ftEdgeBase *inner_edge = 0;
      int max_samples = 80;
      for (size_t kr=1; kr<start_edges.size();  kr++)
      {
	  inner_edge = start_edges[kr].get();
	  while (true)
	  {
	      double parmin = inner_edge->tMin();
	      double parmax = inner_edge->tMax();
	      int nmb_eval = min(max_samples, nmbToEval(dynamic_cast<ftEdge*>(inner_edge),
							parmin, parmax));
	      double tpar;
	      double tint = (parmax-parmin)/(int)(nmb_eval-1);
	      int kh;
	      for (kh=0, tpar=parmin; kh<nmb_eval; kh++, tpar+=tint)
	      {
		  Point pnt = inner_edge->point(tpar);
		  Vector3D pnt3D(pnt[0], pnt[1], pnt[2]);
		  shared_ptr<ftSurfaceSetPoint> ftpnt(new ftSurfaceSetPoint(pnt3D, 2));
		  double* seed = NULL;
		  shared_ptr<Point> par_pt;
		  double tpar2 = tpar;
		  getParPoint(inner_edge, tpar2, par_pt);
		  if (par_pt.get() != 0) {
		      seed = (*par_pt).begin();
		  }
		  addFaceRef(*ftpnt, faces_[i], seed); // edg1 not returned by AdjacentEdges().
		  points.addEntry(ftpnt);
	      }

	      inner_edge = inner_edge->next();
	      if (inner_edge == 0 || inner_edge == start_edges[kr].get())
		  break;
	  }
      }
  }

  return status;
}

//===========================================================================
ftMessage ftSurfaceSet::getInitInnerData2(ftPointSet& points)
//---------------------------------------------------------------------------
//
// Purpose: Get initial data points from the the interior of a set of surfaces.
//
//===========================================================================
{
    ftMessage status;

    // The topological structure of points is set outside this fuction, afterwards.
    //    double step = 20.0*toptol_.neighbour;
    double step = 5.0*toptol_.neighbour;
    int csidx = 0;  // Current surface-index (as given in faces_)

    int i, j, m;
    for (i = 0; i < (int)faces_.size(); ++i) {
	// Index of current surface.
	csidx = faces_[i]->getId();

	shared_ptr<ParamSurface> curr_surface = faces_[i]->surface();
	RectDomain rect_domain = curr_surface->containingDomain();

	double tmin_u = rect_domain.umin();
	double tmax_u = rect_domain.umax();
	double tmin_v = rect_domain.vmin();
	double tmax_v = rect_domain.vmax();

	shared_ptr<SplineSurface> under_surf; // The underlying surface.
	if (faces_[i]->surface()->instanceType() == Class_SplineSurface)
	    under_surf = dynamic_pointer_cast<SplineSurface, ParamSurface>
		(faces_[i]->surface());
	else if (faces_[i]->surface()->instanceType() == Class_BoundedSurface)
	    under_surf = dynamic_pointer_cast<SplineSurface, ParamSurface>
		((dynamic_pointer_cast<BoundedSurface, ParamSurface>
		  (faces_[i]->surface()))->underlyingSurface());
	else {
	    status.setError(FT_UNEXPECTED_DATA_TYPE);
	    return status;
	}
	double length1, length2;
	try {
 	    vector<double> lengths = estimateSurfSides(under_surf.get());
	    length1 = 0.5*(lengths[0] + lengths[1]);
	    length2 = 0.5*(lengths[2] + lengths[3]);
	} catch (...) {
	    status.setError(FT_NO_DATA); // Should never happen!
	    return status;
	}
	//int min_samples = 1;
	int min_samples = 3;
	// The sampling really should make sure the triangles are quite similar,
	// as it seems to be a requirement of the parametrization.
		//int max_samples = 2;
	int max_samples = 40;

	int nmb_eval_u, nmb_eval_v;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	nmb_eval_u = max(min_samples, min(max_samples, (int)(length1/step)));
	nmb_eval_v = max(min_samples, min(max_samples, (int)(length2/step)));
#else
 	nmb_eval_u =
	    std::max(min_samples, std::min(max_samples, (int)(length1/step)));
	nmb_eval_v = 
	    std::max(min_samples, std::min(max_samples, (int)(length2/step)));
#endif
 	double parinc_u = (tmax_u - tmin_u)/(nmb_eval_u+1);
 	double parinc_v = (tmax_v - tmin_v)/(nmb_eval_v+1);

	Point pnt;
	Vector3D pnt3D;
	double par_u, par_v;
	int at_edge = 0;
// 	double tol = 1e-10;
	vector<Vector2D> bd_train = getBdTrain(points, faces_[i]);

	for (j = 1; j < nmb_eval_u + 1; ++j)
	    for (m = 1; m < nmb_eval_v + 1; ++m) {
		par_u = tmin_u + parinc_u * j;
		par_v = tmin_v + parinc_v * m;
		Vector2D par_pt(par_u, par_v);
// 		if (!((curr_surface->parameterDomain()).isInDomain(par_pt, tol)))
// 		    continue; // We do not sample points outside a trimmed surface.
		// @@sbr This may lead to sampling of pts not on actual sf (defined by containing domain).
		// But as they stay inside bd_train it shouldn't really mather.
		// We must futhermore check that we're inside chain defined by bd_pts.
		if (cmUtils::insideBdTrain(par_pt, bd_train)) {
		   curr_surface->point(pnt, par_u, par_v);
		   pnt3D = Vector3D(pnt[0], pnt[1], pnt[2]);
		   shared_ptr<ftSurfaceSetPoint>
		      ftpnt(new ftSurfaceSetPoint(pnt3D, at_edge,
						  faces_[csidx], par_pt));
		   points.addEntry(ftpnt);
		}
	    }
    }

    return status;
}


//===========================================================================
ftMessage ftSurfaceSet::updatePointTopology(ftPointSet& points)
//===========================================================================
{
    ftMessage status;
    int i, j, m;
    vector<vector<ttlPoint*> > faces_points(faces_.size());
    for (i = 0; i < (int)faces_.size(); ++i) {
	// For each surface we rescale parameter domain of spline/underlying surface.
	// Parameter values are then mapped accordingly.
	RectDomain rect_dom = faces_[i]->surface()->containingDomain();
	RectDomain new_dom = cmUtils::geometricParamDomain(faces_[i]->surface().get());
	for (j = 0; j < points.size(); ++j) {
	    ftSurfaceSetPoint* sspnt = dynamic_cast<ftSurfaceSetPoint*>(points[j]);
	    if (sspnt == 0) {
		status.setError(FT_UNEXPECTED_DATA_TYPE);
		return status;
	    }
	    for (m = 0; m < sspnt->nmbFaces(); ++m) {
		if (faces_[i].get() == sspnt->face(m).get()) {
		    double u = sspnt->parValue(m)[0];
		    double v = sspnt->parValue(m)[1];
		    double new_u = (new_dom.umax() - new_dom.umin())*(u - rect_dom.umin())/
			(rect_dom.umax() - rect_dom.umin()) + new_dom.umin();
		    double new_v = (new_dom.vmax() - new_dom.vmin())*(v - rect_dom.vmin())/
			(rect_dom.vmax() - rect_dom.vmin()) + new_dom.vmin();
		    faces_points[i].push_back(new ttlPoint(sspnt, new_u, new_v));
		}
	    }
	}
    }

    vector<hetriang::Triangulation> triang(faces_points.size());
    bool missing_face_points = false;
    for (i = 0; i < (int)faces_.size(); ++i) {
        if (faces_points[i].size() == 0)
        {
            missing_face_points = true;
            break;
        }
	triang[i].createDelaunay(faces_points[i].begin(),
				 faces_points[i].end());

#ifdef FANTASTIC_DEBUG
	std::ofstream debug("data/debug.dat");
	triang[i].printEdges(debug);
#endif // FANTASTIC_DEBUG
    }

    if (missing_face_points)
    {
        for (i = 0; i < (int)faces_.size(); ++i) {
            for (j = 0; j < (int)faces_points[i].size(); ++j)
            {
                if (faces_points[i][j])
                {
                    delete faces_points[i][j]; // No more need for object.
                }
            }
        }
        status.setError(FT_TOPOLOGY_PROBLEM);
        return status;
    }


    // We run through the vector, updating structure for each face.
//     double epsge = 1e-05;
     for (i = 0; i < (int)triang.size(); ++i) {
	const list<hetriang::Edge*>& l_edges = triang[i].getLeadingEdges();
	// For triangulation of current face, we run through all triangles.
	list<hetriang::Edge*>::const_iterator leading_edge_it = l_edges.begin();
	while (leading_edge_it != l_edges.end()) {
	    hetriang::Edge* curr_edge = *leading_edge_it;
	    shared_ptr<hetriang::Node> source_node = curr_edge->getSourceNode();
	    shared_ptr<hetriang::Node> target_node = curr_edge->getTargetNode();
	    for (m = 0; m < 3; ++m) {
		std::vector<PointIter> neighbours =
		    source_node->pointIter()->getNeighbours();
		for (j = 0; j < (int)neighbours.size(); ++j)
		    if (target_node->pointIter() == neighbours[j])
			break;
		// If break was executed, connection already exists.
		// If both nodes are on boundary we dont't make the points they
		// are referring to neighbours (as all boundary points already
		// have got their maximum of two boundary neighbours).
		// We could have allowed two points on a subsurfaceboundary
		// to be neighbours, but it would reault in a conflict when
		// topology is to be used in the context of a graph.
// 		double source_u = source_node->x();
// 		double source_v = source_node->y();
// 		double target_u = target_node->x();
// 		double target_v = target_node->y();
		if (j == (int)neighbours.size())
		  // 	 if (!(source_node->pointIter()->isOnBoundary()
		  // 	  && target_node->pointIter()->isOnBoundary()))
		  if (!((source_node->pointIter()->isOnSubSurfaceBoundary()
			 && target_node->pointIter()->isOnSubSurfaceBoundary())))
// 			  // We expect boundary params to be exact.
// 		      && (fabs(source_u - target_u) > epsge &&
// 			  fabs(source_v - target_v) > epsge))
		    // Add neighbour (in one direction at the time).
		    (source_node->pointIter())->
		      addNeighbour(target_node->pointIter());
		if (m == 2)
		    break;
		curr_edge = curr_edge->getNextEdgeInFace();
		source_node = curr_edge->getSourceNode();
		target_node = curr_edge->getTargetNode();
	    }
	    ++leading_edge_it; // We iterate to next triangle.
	}
    }

    return status;
}

//===========================================================================
ftMessage ftSurfaceSet::merge(double& max_error, double& mean_error,
			      double bend_tol, double kink_tol)
//---------------------------------------------------------------------------
//
// Purpose: Merge the current surface patches to produce one
//          SplineSurface
//
//===========================================================================
{
  int ki;
  ftMessage status, local_status;
  // Check if the surface is created already
  if (surf_.get() != 0)
    {
      // Nothing to do
      status.addWarning(FT_SURFACE_ALREADY_CREATED);
      return status;
    }

  vector<ftEdgeBase*> outer_loop;
  vector<int> corner;
  local_status = getOuterLoop(outer_loop, corner, bend_tol, kink_tol);

#ifdef FANTASTIC_DEBUG
  std::ofstream debug("data/debug.g2");
  for (ki = 0; ki < (int)outer_loop.size(); ++ki) {
      shared_ptr<ParamCurve> cv = outer_loop[ki]->geomEdge()->geomCurve();
      if (cv->instanceType() == Class_SplineCurve) {
	  cv->writeStandardHeader(debug);
	  cv->write(debug);
      } else if (cv->instanceType() == Class_CurveOnSurface) {
	  shared_ptr<CurveOnSurface> cv_on_sf =
	      dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
	  shared_ptr<SplineCurve> par_cv =
	      dynamic_pointer_cast<SplineCurve, ParamCurve>(cv_on_sf->parameterCurve());
	  if (par_cv.get() != 0) {
	      SplineDebugUtils::writeSpaceParamCurve(*par_cv, debug);
	  }
	  if (cv_on_sf->spaceCurve().get() != 0) {
	      cv_on_sf->spaceCurve()->writeStandardHeader(debug);
	      cv_on_sf->spaceCurve()->write(debug);
	  }
      }
  }
#endif // FANTASTIC_DEBUG
  
  int no_warn = local_status.noOfWarnings();
  for (ki=0; ki<no_warn; ki++)
      status.addWarning(local_status.getWarning(ki));
  if (!local_status.isOK()) {
      status.setError(local_status.getMessage());
      return status;
  }

  local_status = makeSurface(outer_loop, corner, max_error, mean_error);

  if (!local_status.isOK())
    status.setError(local_status.getMessage());
  no_warn = local_status.noOfWarnings();
  for (ki=0; ki<no_warn; ki++)
    status.addWarning(local_status.getWarning(ki));

  return status;
}

//===========================================================================
vector<double> ftSurfaceSet::estimateSurfSides(SplineSurface* surf)
//===========================================================================
{
    ALWAYS_ERROR_IF(surf == 0,"Surface does not exist!");

    vector<double> surf_sides(4); // Lengths of surf sides.

    // Sampling in nmb1 and nmb2 points along each boundary curve
    int nmb1 = 4, nmb2 = 4;
    int ki;

    // 1. parameter direction
    double upar, vpar;
    double length1 = 0.0;
    vpar = surf->startparam_v();
    upar = surf->startparam_u();
    double udel = (surf->endparam_u() - upar)/(double)(nmb1-1);
    double vdel = (surf->endparam_v() - vpar)/(double)(nmb2-1);
    Point pt1, pt2;
    surf->point(pt1, upar, vpar);
    for (ki=1, upar+= udel; ki<nmb1; ki++, upar+=udel)
	{
	    surf->point(pt2, upar, vpar);
	    length1 += pt1.dist(pt2);
	    pt1 = pt2;
	}
    surf_sides[0] = length1;

    length1 = 0.0;
    vpar = surf->endparam_v();
    upar = surf->startparam_u();
    surf->point(pt1, upar, vpar);
    for (ki=1, upar+= udel; ki<nmb1; ki++, upar+=udel)
	{
	    surf->point(pt2, upar, vpar);
	    length1 += pt1.dist(pt2);
	    pt1 = pt2;
	}
    surf_sides[1] = length1;

    // 2. parameter direction
    double length2 = 0.0;
    vpar = surf->startparam_v();
    upar = surf->startparam_u();
    surf->point(pt1, upar, vpar);
    for (ki=1, vpar+=vdel; ki<nmb2; ki++, vpar+=vdel)
	{
	    surf->point(pt2, upar, vpar);
	    length2 += pt1.dist(pt2);
	    pt1 = pt2;
	}
    surf_sides[2] = length2;

    length2 = 0.0;
    vpar = surf->startparam_v();
    upar = surf->endparam_u();
    surf->point(pt1, upar, vpar);
    for (ki=1, vpar+=vdel; ki<nmb2; ki++, vpar+=vdel)
	{
	    surf->point(pt2, upar, vpar);
	    length2 += pt1.dist(pt2);
	    pt1 = pt2;
	}
    surf_sides[3] = length2;

    return surf_sides;
}


//===========================================================================
ftMessage ftSurfaceSet::reparametrizeSurf(const vector<ftEdgeBase*>& edgeloop,
					  const vector<int>& corners,
					  ftPointSet& points, double umin, double umax,
					  double vmin, double vmax)
//===========================================================================
{
    ftMessage status;
    double epsge = 1e-10;
    double minlen = 0.01;

    bool compute = (umax-umin<minlen || vmax-vmin<minlen);

    if (corners.size() != 5) {
	status.setError(FT_NON_4_SIDED_SURF);
	return status;
    }

    // We estimate length of each of the four boundary curves.
    // @@sbr Not that smart when surf is trimmed in two ends.
    vector<double> lengths(corners.size() - 1, 0.0);
    int nmb = (int)edgeloop.size();
//     int side = 0;
    int i;
    //    lengths[side] = edgeloop[0]->geomEdge()->estimatedCurveLength();
    for (i = 0; i < (int)corners.size() - 1; ++i) {
      int j = corners[i];
      while (j < corners[i+1]) {
	  lengths[i] += cmUtils::estimatedCurveLength(edgeloop[j%nmb]);//->geomEdge()->estimatedCurveLength();	
	  ++j;
      }
    }

    // If surface is trimmed in both ends, we need info on mid lengths.
    double mid_length_u = surf_->constParamCurve(0.5*(surf_->startparam_v()+ surf_->endparam_v()), true)
	->ParamCurve::estimatedCurveLength();
    double mid_length_v = surf_->constParamCurve(0.5*(surf_->startparam_u()+ surf_->endparam_u()), false)
	->ParamCurve::estimatedCurveLength();
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    double length_u = max(minlen, max(mid_length_u, max(lengths[0], lengths[2])));
    double length_v = max(minlen, max(mid_length_v, max(lengths[1], lengths[3])));
#else
    double length_u = std::max(minlen, max(mid_length_u, max(lengths[0], lengths[2])));
    double length_v = std::max(minlen, max(mid_length_v, max(lengths[1], lengths[3])));
#endif

    double umax_old = surf_->endparam_u();
    double vmax_old = surf_->endparam_v();
    double umin2 = surf_->startparam_u();
    double vmin2 = surf_->startparam_v();
    if (compute)
      surf_->setParameterDomain(umin2, umin2 + length_u, vmin2, vmin2 + length_v);
    else
      surf_->setParameterDomain(umin, umax, vmin, vmax);

     umin2 = surf_->startparam_u();
     vmin2 = surf_->startparam_v();
    double umax2 = surf_->endparam_u();
    double vmax2 = surf_->endparam_v();
    length_u = umax2 - umin2;
    length_v = vmax2 - vmin2;

    // We next must update the parameter values of the sampled points according
    // to the new parameter domain.
    double u_frac = length_u / (umax_old - umin2);
    double v_frac = length_v / (vmax_old - vmin2);
    for (i = 0; i < points.size(); ++i) {
      double u = umin2 + (points.getU(i) - umin2) * u_frac;
      double v = vmin2 + (points.getV(i) - vmin2) * v_frac;
      if (points[i]->isOnBoundary()) { // If bnd point, make sure it stays on bnd.
	if (fabs(u - umin2) < epsge)
	  u = umin2;
	else if (fabs(u - umax2) < epsge)
	  u = umax2;
	if (fabs(v - vmin2) < epsge)
	  v = vmin2;
	else if (fabs(v - vmax2) < epsge)
	  v = vmax2;
      }
      points.setU(i, u);
      points.setV(i, v);
    }

    return status;
}

//===========================================================================
void ftSurfaceSet::buildTopology()
//---------------------------------------------------------------------------
//
// Purpose: Find adjacency between faces and build a topology table 
//          representing this adjacency.
//
//===========================================================================
{
    // Perform adjacency analysis
    FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(toptol_);
    adjacency.computeAdjacency(faces_, 0);
}

//===========================================================================
void ftSurfaceSet::BoundaryLoops(vector< vector<ftEdgeBase*> >& loopvec)
//---------------------------------------------------------------------------
//
// Purpose: Fetch all boundary loops of the current ftSuperSurface. Start
//          from first_edges_.
//
//===========================================================================
{
  // If first_edges_ has not been set we give an error. Object must be handled
  // by a tpTopologyTable (which sets and stores first_edges_)!
  if (boundary_loops_.size() == 0)
      THROW("boundary_loops_ was not set.");

  loopvec.resize(boundary_loops_.size());
  for (size_t ki=0; ki<boundary_loops_.size(); ki++)
    {
      for (size_t kj=0; kj<boundary_loops_[ki]->size(); ++kj)
	loopvec[ki].push_back(boundary_loops_[ki]->getEdge(kj).get());
    }
}

//===========================================================================
vector<Vector2D> ftSurfaceSet::getBdTrain(ftPointSet& pt_set, shared_ptr<ftFaceBase>& face)
//===========================================================================
{
   vector<Vector2D> bd_train;
   // We run through pt_set locating element with the input idx.
   int ki, kj;
   ftSurfaceSetPoint* first_pt = 0;
   for (ki = 0; ki < pt_set.size(); ++ki) {
      ftSurfaceSetPoint* pt = dynamic_cast<ftSurfaceSetPoint*>(pt_set[ki]);
      ASSERT(pt != 0);
      for (kj = 0; kj < pt->nmbFaces(); ++kj)
	 if (pt->face(kj) == face) {
	    first_pt = pt;
	    bd_train.push_back(first_pt->parValue(kj));
	    break;
	 }
      if (kj < pt->nmbFaces()) {
	 break;
      }
   }
   if (first_pt == 0) {
      return bd_train;
   }

   // We then locate neighbour pt with the same idx.
   ftSurfaceSetPoint* curr_pt = NULL;
   vector<PointIter> neighbours = first_pt->getNeighbours();
   for (ki = 0; ki < (int)neighbours.size(); ++ki) {
      ftSurfaceSetPoint* pt = dynamic_cast<ftSurfaceSetPoint*>(neighbours[ki]);
      ASSERT(pt != 0);
      for (kj = 0; kj < pt->nmbFaces(); ++kj)
	 if (pt->face(kj) == face) {
	    curr_pt = pt;
	    bd_train.push_back(curr_pt->parValue(kj));
	    break;
	 }
      if (kj < pt->nmbFaces()) {
	 break;
      }
      if (curr_pt == pt) {
	 break;
      }
   }
   ASSERT(curr_pt != 0);

   ftSurfaceSetPoint* prev_pt = first_pt;
   while (curr_pt != first_pt) {
      neighbours = curr_pt->getNeighbours();
      for (ki = 0; ki < (int)neighbours.size(); ++ki) {
	 ftSurfaceSetPoint* pt = dynamic_cast<ftSurfaceSetPoint*>(neighbours[ki]);
	 ASSERT(pt != 0);
	 for (kj = 0; kj < pt->nmbFaces(); ++kj)
	    if ((pt->face(kj) == face) && (neighbours[ki] != prev_pt)) {
	       prev_pt = curr_pt;
	       curr_pt = pt;
	       if (curr_pt != first_pt) {
		  bd_train.push_back(curr_pt->parValue(kj));
	       }	
	       break;
	    }
	 if (kj < pt->nmbFaces()) {
	    break;
	 }
      }
   }

   return bd_train;
}


//===========================================================================
void ftSurfaceSet::getParPoint(ftEdgeBase* edge, double tpar,
			       shared_ptr<Point>& par_pt)
//===========================================================================
{
  if (edge->geomEdge()->geomCurve()->instanceType() == Class_CurveOnSurface) {
      shared_ptr<CurveOnSurface> cv_on_sf =
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(edge->geomEdge()->geomCurve());
      if (cv_on_sf->parPref()) {
	  par_pt = shared_ptr<Point>(new Point(cv_on_sf->parameterCurve()->point(tpar)));
      }
  }
}

//===========================================================================
SplineSurface* ftSurfaceSet::createContainingSurface(ParamSurface* sf)
//===========================================================================
{
    // If trim cvs are rational, then what?
    double knot_diff_tol = 1e-05;
    SplineSurface* cont_sf = NULL;
    if (sf->instanceType() == Class_SplineSurface) {
	cont_sf = dynamic_cast<SplineSurface*>(sf->clone());
    } else if (sf->instanceType() == Class_BoundedSurface) {
	BoundedSurface* bd_sf = dynamic_cast<BoundedSurface*>(sf);
#if _MSC_VER > 0 && _MSC_VER < 1300
	const CurveBoundedDomain& bd_dom = 
	  dynamic_cast<const CurveBoundedDomain&>(bd_sf->parameterDomain());
#else
	CurveBoundedDomain bd_dom = bd_sf->parameterDomain();
#endif
	RectDomain cont_dom = bd_dom.containingDomain();
	shared_ptr<SplineSurface> under_sf =
	    dynamic_pointer_cast<SplineSurface, ParamSurface>(bd_sf->underlyingSurface());
	ASSERT(under_sf.get() != 0);
#if _MSC_VER > 0 && _MSC_VER < 1300
	const RectDomain& under_dom =  dynamic_cast<const RectDomain&>(under_sf->parameterDomain());
#else
	RectDomain under_dom = under_sf->parameterDomain();
#endif
	double umin = max(cont_dom.umin(), under_dom.umin());
	double vmin = max(cont_dom.vmin(), under_dom.vmin());
	double umax = min(cont_dom.umax(), under_dom.umax());
	double vmax = min(cont_dom.vmax(), under_dom.vmax());
	cont_sf = under_sf->subSurface(umin, vmin, umax, vmax, knot_diff_tol);
    }

    return cont_sf;
}

} // namespace Go
