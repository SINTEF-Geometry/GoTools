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

#include "GoTools/compositemodel/ftChartSurface.h"

#include "GoTools/parametrization/PrPlanarGraph_OP.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
#include "GoTools/parametrization/PrPrmUniform.h"
#include "GoTools/parametrization/PrOrganizedPoints.h"
#include "GoTools/compositemodel/ftSmoothSurf.h"
#include "GoTools/compositemodel/ftSSfEdge.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/utils/CurvatureUtils.h"
#include "GoTools/geometry/HermiteInterpolator.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/topology/FaceConnectivityUtils.h"

#include "GoTools/compositemodel/cmUtils.h"

//#include "GoTools/model_toolbox/ftUtils.h"

#include <cstdio>

#define FANTASTIC_DEBUG

using std::vector;
using std::make_pair;
using std::min;
using std::max;
using std::swap;

namespace Go
{

//===========================================================================
ftChartSurface::ftChartSurface(const vector<shared_ptr<ParamSurface> >& surfaces,
                               tpTolerances& topeps, double approxeps,
			       double curvature_tol, int m, int n)
    : ftSurfaceSet(surfaces, topeps, approxeps),
      curvature_tol_(curvature_tol), maxerror_(0.0), meanerror_(0.0), m_(m), n_(n),
      symm_distr_functions_(false)
//===========================================================================
{
   grid_distr_functions_.resize(4);
   secn_distr_.resize(4);
}

//===========================================================================
ftChartSurface::~ftChartSurface()
//===========================================================================
{
}

//---------------------------------------------------------------------------
vector<shared_ptr<ftEdgeBase> >
ftChartSurface::createInitialEdges(double degenerate_epsilon,
				   double kink, bool no_split)
//---------------------------------------------------------------------------
{
  vector<shared_ptr<ftEdgeBase> > return_edges;
  if (boundary_loops_.size() == 0)
    {
      FaceConnectivityUtils<ftEdgeBase,ftFaceBase> connectivity;
      vector< vector<ftEdgeBase*> > loopvec;
      connectivity.BoundaryLoops(faces_, loopvec);
	
      // For every loop
      ftEdgeBase* e;
      ftEdge* geomedge;
      for (size_t i = 0; i < loopvec.size(); ++i) 
	{
	  vector<shared_ptr<ftEdgeBase> > edges;
	  int n = (int)loopvec[i].size();
	  edges.reserve(n);
	  for (int j = 0; j < n; ++j) 
	    {
	      geomedge = dynamic_cast<ftEdge*>(loopvec[i][j]);

	      e = new ftSSfEdge(this, geomedge);
	      edges.push_back(shared_ptr<ftEdgeBase> (e));
	      if (j > 0) {
		int ne = (int)edges.size();
		edges[ne - 1]->connectAfter(edges[ne - 2].get());
	      }
	    }
	  int ne = (int)edges.size();
	  edges[0]->closeLoop(edges[ne - 1].get());
	  shared_ptr<Loop> curr_loop = shared_ptr<Loop>(new Loop(this, edges, degenerate_epsilon));
	  boundary_loops_.push_back(curr_loop);

	  // We make sure that our chosen corner pts do not lie in the interior of an edge.
	  cmUtils::splitEdgesInCorners(edges, additional_corner_pts_, degenerate_epsilon);

	  return_edges.insert(return_edges.begin(), edges.begin(), edges.end());
	}
    }
  else
    {
      for (size_t ki = 0; ki < boundary_loops_.size(); ++ki) 
	{
	  shared_ptr<Loop> curr_loop = boundary_loops_[ki];
	  vector<shared_ptr<ftEdgeBase> > curr_edges = curr_loop->getEdges();
	  return_edges.insert(return_edges.end(), curr_edges.begin(), curr_edges.end());
	}
    }
  return return_edges;
}

//===========================================================================
Point ftChartSurface::point(double u, double v) const
//===========================================================================
{
    // We find local parameters and face in graph.
    double new_u = u;
    double new_v = v;
    shared_ptr<ftFaceBase> dummy_face;
    Point clo_pt = point(new_u, new_v, dummy_face);

    return clo_pt;
}

//===========================================================================
Point ftChartSurface::point(double& u, double& v, shared_ptr<ftFaceBase>& face,
			    double* seed, bool use_input_face) const
//===========================================================================
{
    Point space_pt = surf_->ParamSurface::point(u, v);
    if (!use_input_face) {
	graph_.getLocalParameters(u, v, face);
#if 0
        // If the boundary curve is not approximated strict enough the closest point evaluation approach
        // may fail due to multiple surface points projecting to the same parameter point (if boundary of
        // surf_ extends outside the surface set). For this scenario it may be beneficial to rely on a
        // good parametrization mapping instead of relying on the space pt in the approximating surf_.
        space_pt = face->point(u, v);
        return space_pt;
 #endif
    } else {
        MESSAGE("Using input face! Why is this not an option for the normal() function?");
	ASSERT(face.get() != 0 && seed != NULL);
	u = seed[0];
	v = seed[1];
    }
    
#ifndef NDEBUG
    if (0)
    {
        Point debug_local_pt = face->point(u,v);
        double debug_dist = space_pt.dist(debug_local_pt);
        if (debug_dist > 0.07) // This value matches current case (fanta_ro2_sub.g2) ...
        {
            MESSAGE("DEBUG: dist from global to local pt: " << debug_dist);
        }
    }
#endif
    
    // Local parameters are used as seed in a closest point iteration on found face.
    double clo_u, clo_v, clo_dist;
    double clo_u_bd, clo_v_bd, clo_dist_bd;
    Point clo_pt, clo_pt_bd;
    Vector2D par_pt(u, v);
    double bd_tol = 1e-6; // @@sbr Hardcoded value!
    double knot_tol = 1e-12;
    if (seed == 0) {
	seed = par_pt.begin();
    }
    bool bd_pt = false;
    try {
	bd_pt = face->surface()->parameterDomain().isOnBoundary(par_pt, bd_tol);//knot_tol);
    } catch (...) {
	// 
    }
    // @@sbr201706 If the boundary point is along an inner edge we should check the distance for the adjacent surface.
//    if (bd_pt) {
	face->surface()->closestBoundaryPoint(space_pt, clo_u_bd, clo_v_bd, clo_pt_bd, clo_dist_bd,
					      bd_tol, NULL, seed);
        // We find the topological edge.
        
        //MESSAGE("Evaluating in a boundary point: clo_u: " << clo_u << ", clo_v: " << clo_v);
        //  } else {
	face->surface()->closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist,
				      bd_tol, NULL, seed);
        //}
        if (clo_dist_bd < clo_dist)
        {
            clo_u = clo_u_bd;
            clo_v = clo_v_bd;
            clo_pt = clo_pt_bd;
            clo_dist = clo_dist_bd;
        }
#if 1
        {
            if (clo_dist > 1.0e-02)
            {
//    if ((clo_u == 0.0 && u > 0.0) || (clo_v == 0.0 && v > 0.0)) {
                MESSAGE("u: " << u << ", v: " << v << ", clo_u: " << clo_u << ", clo_v: " << clo_v <<
                        ", clo_dist: " << clo_dist);
            }
        }
#endif
    
#ifdef FANTASTIC_DEBUG
    std::ofstream debug("tmp/debug.g2");
    vector<double> pts(6);
    copy(space_pt.begin(), space_pt.end(), pts.begin());
    copy(clo_pt.begin(), clo_pt.end(), pts.begin() + 3);
    LineCloud lc(pts.begin(), 1);
    lc.writeStandardHeader(debug);
    lc.write(debug);
    Point seed_space_pt = face->surface()->point(u, v);
    copy(space_pt.begin(), space_pt.begin() + 3, pts.begin());
    copy(seed_space_pt.begin(), seed_space_pt.end(), pts.begin() + 3);
    LineCloud lc2(pts.begin(), 1);
    lc2.writeStandardHeader(debug);
    lc2.write(debug);
#endif // FANTASTIC_DEBUG

    u = clo_u;
    v = clo_v;
    return clo_pt;
}

//===========================================================================
Point ftChartSurface::normal(double u, double v) const
//===========================================================================
{
    Point normal;
    if (surf_.get() == 0) {
        THROW("Graph not created yet!");
    }
    // We find local parameters and face in graph.
    double new_u = u;
    double new_v = v;
    shared_ptr<ftFaceBase> face;
    graph_.getLocalParameters(new_u, new_v, face);

    // Local parameters are used as seed in a closest point iteration on found face.
    Point space_pt = surf_->ParamSurface::point(u, v);
    double clo_u, clo_v, clo_dist;
    double clo_u_bd, clo_v_bd, clo_dist_bd;
    Point clo_pt, clo_pt_bd;
    Vector2D par_pt(u, v);
    double bd_tol = 1e-6; // @@sbr Hardcoded value!
    double knot_tol = 1e-12;
    double seed[2];
    seed[0] = new_u;
    seed[1] = new_v;
    bool bd_pt = false;
    try {
	bd_pt = face->surface()->parameterDomain().isOnBoundary(par_pt, bd_tol);//knot_tol);
    } catch (...) {
	// 
    }
//    if (bd_pt) {
	// @@sbr Make sure found edge is without twin?
	face->surface()->closestBoundaryPoint(space_pt, clo_u_bd, clo_v_bd, clo_pt_bd, clo_dist_bd,
					      bd_tol, NULL, seed);
//    } else {
	face->surface()->closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist,
				      bd_tol, NULL, seed);
        //  }
        if (clo_dist_bd < clo_dist)
        {
            clo_u = clo_u_bd;
            clo_v = clo_v_bd;
        }
#if 1
        {
            if (clo_dist > 1.0e-02)
            {
//    if ((clo_u == 0.0 && u > 0.0) || (clo_v == 0.0 && v > 0.0)) {
                MESSAGE("u: " << u << ", v: " << v << ", clo_u: " << clo_u << ", clo_v: " << clo_v <<
                        ", clo_dist: " << clo_dist);
            }
        }
#endif

    normal = face->normal(clo_u, clo_v);
    return normal;
}

//===========================================================================
ftMessage ftChartSurface::createSurf(double& max_error, double& mean_error)
//---------------------------------------------------------------------------
//
// Purpose: Merge the current surface patches to produce one
//          SplineSurface
//
//===========================================================================
{
  ftMessage status;
  double bend_tol = toptol_.bend;
  double kink_tol = toptol_.kink;
  if (surf_.get() == 0) {
      try {
	  status = merge(max_error, mean_error, bend_tol, kink_tol);
      } catch (...) {
	  status.setError(FT_ERROR_IN_SURFACE_CREATION);
	  return status;
      }
      maxerror_ = max_error;
      meanerror_ = mean_error;
  } else
    status.addWarning(FT_SURFACE_ALREADY_CREATED);

  return status;
}

//---------------------------------------------------------------------------
void ftChartSurface::getError(double& max_error, double& mean_error)
//---------------------------------------------------------------------------
{
  max_error = maxerror_;
  mean_error = meanerror_;
}

//===========================================================================
ftTangPriority ftChartSurface::getPrioType() const
//===========================================================================
{
    return ftNoType;
}

//===========================================================================
void ftChartSurface::setPrioType(ftTangPriority type)
//===========================================================================
{
    prio_type_ = type;
}

//===========================================================================
void ftChartSurface::setGridResolution(int m, int n)
//===========================================================================
{
   m_ = m;
   n_ = n;
}

//===========================================================================
bool ftChartSurface::gridCreated(int& m, int& n) const
//===========================================================================
{
   m = m_;
   n = n_;

   // Should be enough to conclude whether a grid has been made.
   return (m_*n_ > 0 && (int) grid_pts_.size() == m_ && (int) grid_pts_[0].size() == n_);
}

//===========================================================================
void ftChartSurface::setEdgeScales(vector<pair<Point, double> >& edge_scales)
//===========================================================================
{
    edge_scales_ = edge_scales;
}

//===========================================================================
void ftChartSurface::setSymmDistrFunctions(bool symm_distr_functions)
//===========================================================================
{
    symm_distr_functions_ = symm_distr_functions;

}

//===========================================================================
ftMessage ftChartSurface::prepareGrid(RotationInfo* rot_info, ftCurve* outer_bd)
//===========================================================================
{
    ftMessage status;

    int ki;
    // We try to make neighbouring blocks as equal as possible.
    ftMessage local_status;
    if (rot_info != 0) { // Applies to rot_info only.
	local_status = modifyGridDistrFunctions(rot_info, outer_bd);
	for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
	    status.addWarning(local_status.getWarning(ki));
	if (!local_status.isOK()) {
	    status.setError(local_status.getMessage());
	    return status;
	}
    }

    // We make sure that distr functions along common bd match.
    local_status = updateGridDistribution(rot_info, outer_bd);
    for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
	status.addWarning(local_status.getWarning(ki));
    if (!local_status.isOK()) {
	status.setError(local_status.getMessage());
	return status;
    }

    return status;
}

//===========================================================================
ftMessage ftChartSurface::createGrid(RotationInfo* rot_info, ftCurve* outer_bd)
//===========================================================================
{
   ftMessage status;

   ASSERT(m_ > 1 && n_ > 1);

   // We set elements of grid_pts_.
   sampleGridPts();

   // We should then see if any neighbouring sfs of the same type already exist.
   // If one exists along an edge, with the same nmb of samples, we replace sampled
   // pts along edge with those from the neighbour.
   vector<ftEdgeBase*> outer_loop;
   vector<int> corners;
   double bend_tol = toptol_.bend;
   double kink_tol = toptol_.kink;
   ftMessage local_status = getOuterLoop(outer_loop, corners, bend_tol, kink_tol);
   if (!local_status.isOK()) {
      status.setError(local_status.getMessage());
      return status;
   }
   // If corners do not define a rectangular surface, regular gridding will be impossible to achieve
   // (with our current strategy of a m_*n_-mesh).
   if (corners.size() != 5) {
      status.setError(FT_NON_4_SIDED_SURF);
      return status;
   }

   // For each of the four edges we see whether one other surface has corners matching edge.
   vector<pair<int, int> > matching_edges; // We store index of those edges denoting a matching edge.
   vector<pair<int, int> > matching_grid_res; // We also store sampling res along edges.
   vector<ftChartSurface*> matching_faces;
   vector<double> rot_angles; // 0.0 for edges which match without rotation.
   int ki, kj;
   // We extract info about matching edges aorund created sf.
   local_status = getMatchingEdges(outer_loop, corners, matching_edges, matching_grid_res,
				   matching_faces, rot_angles);
   for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
     status.addWarning(local_status.getWarning(ki));
   if (!local_status.isOK()) {
     status.setError(local_status.getMessage());
     return status;
   }

   for (ki = 0; ki < (int) matching_edges.size(); ++ki) {
       int other_res_m, other_res_n;
       if (!(matching_faces[ki]->gridCreated(other_res_m, other_res_n))) {
	 matching_edges.erase(matching_edges.begin() + ki);
	 matching_faces.erase(matching_faces.begin() + ki);
	 rot_angles.erase(rot_angles.begin() + ki);
	 --ki;
       }
   }

   // @@@@@@ VSK @@@@@@ 030704
   // For rotational surfaces, check correspondence along outer edges
   vector<pair<int, int> > new_matching_edges;
   vector<ftChartSurface*> new_matching_faces;
   vector<double> new_rot_angles;
   local_status = getMatchingEdges(outer_loop, corners, rot_info, outer_bd,
				   new_matching_edges, new_matching_faces, new_rot_angles);
   for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
     status.addWarning(local_status.getWarning(ki));
   if (!local_status.isOK()) {
     status.setError(local_status.getMessage());
     return status;
   }

   // If grid is not created for the other sf we remove match.
   for (ki = 0; ki < (int) new_matching_edges.size(); ++ki) {
       int other_res_m, other_res_n;
       if (!(new_matching_faces[ki]->gridCreated(other_res_m, other_res_n))) {
	   new_matching_edges.erase(new_matching_edges.begin() + ki);
	   new_matching_faces.erase(new_matching_faces.begin() + ki);
	   new_rot_angles.erase(new_rot_angles.begin() + ki);
	   --ki;
       }
   }
   matching_edges.insert(matching_edges.end(), new_matching_edges.begin(), new_matching_edges.end());
   matching_faces.insert(matching_faces.end(), new_matching_faces.begin(), new_matching_faces.end());
   rot_angles.insert(rot_angles.end(), new_rot_angles.begin(), new_rot_angles.end());

   // We replace those samples corresponding to matching edges.
   for (ki = 0; ki < (int) matching_edges.size(); ++ki) {
       int face_edge = matching_edges[ki].first;
       int twin_face_edge = matching_edges[ki].second;
       ftChartSurface* twin_face = matching_faces[ki];       
       vector<shared_ptr<ftSamplePoint> > edge_grid = getEdgeGrid(face_edge);
       vector<shared_ptr<ftSamplePoint> > twin_edge_grid = twin_face->getEdgeGrid(twin_face_edge);
       if (edge_grid.size() == twin_edge_grid.size()) {
	   for (kj = 0; kj < (int) edge_grid.size(); ++kj) {
	       Vector3D twin_space_pt = twin_edge_grid[twin_edge_grid.size()-1-kj]->getPoint();
	       if (rot_angles[ki] != 0.0) {
		   Vector3D transvec(rot_info->center_pt_.begin());
		   twin_space_pt -= transvec;
		   GeometryTools::rotatePoint(rot_info->rot_axis_, rot_angles[ki], twin_space_pt.begin());
		   twin_space_pt += transvec;
	       }
	       edge_grid[kj]->setPoint(twin_space_pt);
	       // We then update parameter value.
	       Point space_pt(twin_space_pt[0], twin_space_pt[1], twin_space_pt[2]);
	       double clo_u, clo_v, clo_dist;
	       Point clo_pt;
	       Vector2D curr_par_pt = edge_grid[kj]->getPar();
	       double space_epsilon = 1e-08;
	       surf_->closestBoundaryPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist,
					   space_epsilon, NULL, curr_par_pt.begin());
	       // Not replacing space_pt as par_pt is to be used as a seed in closest point
	       // (not used in surface evaluation).
	       Vector2D new_par_pt(clo_u, clo_v);
	       edge_grid[kj]->setPar(new_par_pt);
	   }
       } else {
	   MESSAGE("Size of meeting grids did not match, using local grid.");	   
       }
   }

   return status;
}


//===========================================================================
vector<shared_ptr<ftSamplePoint> > ftChartSurface::getEdgeGrid(int edge_ind)
//===========================================================================
{
    vector<shared_ptr<ftSamplePoint> > edge_grid;
    int nmb_edge_pts = (edge_ind < 2) ? m_ : n_;
    int start_m = (edge_ind == 3) ? (m_ - 1) : 0;
    int start_n = (edge_ind == 1) ? (n_ - 1) : 0;
    int step_m = (edge_ind < 2) ? 1 : 0;
    int step_n = (edge_ind < 2) ? 0 : 1;
    int ki;
    for (ki = 0; ki < nmb_edge_pts; ++ki) {
	edge_grid.push_back(grid_pts_[start_m+ki*step_m][start_n+ki*step_n]);
    }
    if ((edge_ind == 1) || (edge_ind == 2)) {
	reverse(edge_grid.begin(), edge_grid.end());
    }

    return edge_grid;
}

//===========================================================================
void ftChartSurface::writeGrid(std::ostream& os) const
//===========================================================================
{
   int dummy;
   if (!gridCreated(dummy, dummy)) {
       MESSAGE("Not writing grid as it was not created!");
       return;
   }

   os << n_ << " " << m_ << " " << 1 << std::endl;
   int ki, kj, kk;
   int dim = 3;
   for (ki = 0; ki < dim; ++ki) {
       for (kj = 0; kj < (int) grid_pts_.size(); ++kj) {
	   for (kk = 0; kk < (int) grid_pts_[ki].size(); ++kk) {
	       Vector3D space_pt = grid_pts_[kj][kk]->getPoint();
	       os << space_pt[ki] << " ";
	   }
       }
       os << std::endl;
   }
}

//===========================================================================
void ftChartSurface::writeDebugGrid(std::ostream& os) const
//===========================================================================
{
   int dummy;
   if (!gridCreated(dummy, dummy)) {
       MESSAGE("Not writing grid as it was not created!");
       return;
   }

   vector<double> pts;
   int ki, kj;
   for (ki = 0; ki < m_; ++ki)
      for (kj = 0; kj < n_; ++kj) {
	 if (ki < m_ - 1) {
	    Vector3D from = grid_pts_[ki][kj]->getPoint();
	    Vector3D to = grid_pts_[ki+1][kj]->getPoint();
	    pts.insert(pts.end(), from.begin(), from.end());
	    pts.insert(pts.end(), to.begin(), to.end());
	 }
	 if (kj < n_ - 1) {
	    Vector3D from = grid_pts_[ki][kj]->getPoint();
	    Vector3D to = grid_pts_[ki][kj+1]->getPoint();
	    pts.insert(pts.end(), from.begin(), from.end());
	    pts.insert(pts.end(), to.begin(), to.end());
	 }
      }


   int nmb_lines = (int)pts.size()/6;
   LineCloud lc(pts.begin(), nmb_lines);
   lc.writeStandardHeader(os);
   lc.write(os);
}


//===========================================================================
vector<shared_ptr<FaceConnectivity<ftEdgeBase> > > ftChartSurface::getInnerEdgeCont() const
//===========================================================================
{
    MESSAGE("Under construction!");
    vector<shared_ptr<FaceConnectivity<ftEdgeBase> > > inner_edge_conn;

    // The boundary_loops_ are wrt the surface set. We need boundary loops for each of the surfaces in
    // the set.
    for (size_t ki = 0; ki < faces_.size(); ++ki)
    {
        shared_ptr<ftFaceBase> face = faces_[ki];
        std::vector<shared_ptr<ftEdgeBase> > start_edges = face->startEdges();
        for (size_t kj = 0; kj < start_edges.size(); ++kj)
        {
            ftEdgeBase* first_edge = start_edges[kj].get();
            ftEdgeBase* curr_edge = first_edge;
            //ftEdgeBase* next_edge = first_edge->next();
            while (true)
            {
                // We only store edges with a twin.
                if (curr_edge->twin() != NULL)
                {
                    shared_ptr<FaceConnectivity<ftEdgeBase> > conn_info = curr_edge->getConnectivityInfo();
                    inner_edge_conn.push_back(conn_info);
                }

                curr_edge = curr_edge->next();
                // next_edge = curr_edge->next();
                if (curr_edge == first_edge) {
                    break;
                }
            }
        }
    }

    return inner_edge_conn;
}


//===========================================================================
bool ftChartSurface::sampleOnlyEdges() const
//===========================================================================
{
    return false;
}

//===========================================================================
ftMessage
ftChartSurface::makeSurface(const vector<ftEdgeBase*>& edgeloop,
			    vector<int>& corner,
			    double& max_error, double& mean_error)
//---------------------------------------------------------------------------
//
// Purpose: Approximate surface set with rectangular surface.
//
//===========================================================================
{
  int ki;
  ftMessage status, local_status;
  
      // Fetch the boundary curves as SplineCurves
  vector< shared_ptr<SplineCurve> > bd_curves;
  vector< shared_ptr<SplineCurve> > cross_curves;
  bd_curves.reserve(4);
  cross_curves.reserve(4);  

      // Fetch the curves and cross tangents from the boundary edges.
  vector<BoundaryPiece> bdpiece;
  // We require return cvs to be cord length parametrized.
  status = getBoundaryConditions(edgeloop, corner, bd_curves, cross_curves,
				 false, bdpiece, true);

  if (cross_curves.size() != bd_curves.size()) {
     MESSAGE("Size of bd and cross curves do not match!");
     cross_curves.resize(4);
  }

  // We make a copy of the approximated curves.
  approx_bd_curves_.clear();
  for (ki = 0; ki < (int) bd_curves.size(); ++ki)
    {
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
       approx_bd_curves_.push_back
	 (shared_ptr<SplineCurve>
	  (dynamic_cast<SplineCurve*>(bd_curves[ki]->clone())));
#else    
       approx_bd_curves_.push_back
	 (shared_ptr<SplineCurve>(bd_curves[ki]->clone()));
#endif
    }

  // @@sbr Highly experimental routine to handle bd_cvs with high curvature!
  bool debug = false;
  if (debug)
    {
      try {
	vector<shared_ptr<SplineCurve> > rep_bd_cvs;
	cmUtils::reparametrizeBdCvs2(bd_curves, approxtol_, rep_bd_cvs);
	bd_curves = rep_bd_cvs;
      }
      catch (...) {
	MESSAGE("Method failed, using existing curves.");
      }
    }

#ifndef NDEBUG
  {
      std::ofstream fileout_dbg("tmp/bd_cvs.g2");
      for (size_t ki = 0; ki < bd_curves.size(); ++ki)
      {
          bd_curves[ki]->writeStandardHeader(fileout_dbg);
          bd_curves[ki]->write(fileout_dbg);
      }
  }
#endif
  
  // We construct the initial surface by interpolating boundary curves.
  // @@sbr Should approximate boundary curves and store these prior to interpolation.
  vector<int> edge_derivs(4, 1); // We keep the boundary fixed.
  try {
    // Create a Coons patch from the boundary curves.

    vector<shared_ptr<ParamCurve> > boundary;
    for (ki=0; ki<(int)bd_curves.size(); ki++)
      boundary.push_back(bd_curves[ki]);
    CurveLoop loop(boundary, toptol_.neighbour);
    surf_ = shared_ptr<SplineSurface>
      (CoonsPatchGen::createCoonsPatch(loop));
    
  } catch(...) {
    status.setError(FT_ERROR_IN_SURFACE_CREATION);
    return status;
  }

  // Sample data points representing the shape of the current ftChartSurface.
  vector<int> cn;
  shared_ptr<ftPointSet> points = shared_ptr<ftPointSet>(new ftPointSet());
  try
    {
      local_status = fetchSamplePoints(edgeloop, corner, cn, *points);
    }
  catch (...)
    {
        status.setError(FT_ERROR_IN_SURFACE_CREATION);
        return status;
    }
  for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
    status.addWarning(local_status.getWarning(ki));
  if (!local_status.isOK()) {
    status.setError(local_status.getMessage());
    return status;
  }
  vector<ftSamplePoint*> graph_nodes;
  for (ki = 0; ki < points->size(); ++ki) {
      if ((*points)[ki]->isOnSubSurfaceBoundary()) {
	  graph_nodes.push_back((*points)[ki]);
      }
  }
  // Before we calculate an approximating surf_, we add some more points to 
  // the outer boundary. These additional pts will not be used in the graph.
//   // As the boundary points are to be used in the 
//   // approximation (and not in the graph), we do not mark them as boundary 
//   // points (thus they will not become nodes), and do not set/update any 
//   // neighbourhood information. Number of new points for an edge depends on 
//   // distance between sampled points.
  // The new points must have only two bd neighbours (in local face).
  addOuterBoundaryPoints(*points);

  // If the method somehow failed fetching points we exit.
  if (points->size() == 0)
  {
      status.setError(FT_ERROR_IN_SURFACE_CREATION);
      return status;
  }

  // We set up topology for inner points (both inner-inner & inner-outer, using TTL).
  local_status = updatePointTopology(*points);
  for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
    status.addWarning(local_status.getWarning(ki));
  if (!local_status.isOK()) {
    status.setError(local_status.getMessage());
    return status;
  }

  // Check if the surface is degenerate using a sensible tolerance
  // before this is done with an arbitrary value.
  bool b, r, t, l;
  surf_->isDegenerate(b, r, t, l, toptol_.gap);

  // We must make sure that the ftPointSet has the neighbour
  // structure the way the parametrization code expects it.
  points->orderNeighbours();

  // Parametrize boundary points.
  // The parametrization routine seems to work best with the unit domain.
  surf_->setParameterDomain(0.0, 1.0, 0.0, 1.0);
  // Parameterize
  //   ftParameterize param_instance(points, PLANAR_GRAPH, domain);
  //PrPrmShpPres par;
  PrPrmUniform par;
  PrParametrizeBdy bdy;
  shared_ptr<PrOrganizedPoints> op = shared_ptr<PrOrganizedPoints>(points);

#ifdef FANTASTIC_DEBUG
  std::ofstream pointsout("tmp/pointsdump.dat");
  std::ofstream edgessout("tmp/triangedges.dat");
  points->printXYZNodes(pointsout, true);
  points->printXYZEdges(edgessout);
#endif // FANTASTIC_DEBUG

  try {
    bdy.attach(op);
    //    bdy.setParamKind(PrUNIFBDY);
    double umin = surf_->startparam_u();
    double umax = surf_->endparam_u();
    double vmin = surf_->startparam_v();
    double vmax = surf_->endparam_v();
    if (true)
      bdy.parametrize(cn[0], cn[1], cn[2], cn[3], umin, umax, vmin, vmax);
    else
      points->reparBdy(surf_);

#ifdef FANTASTIC_DEBUG
    points->orderNeighbours();
    std::ofstream debug("tmp/debug.g2");
    vector<double> pts;
    PointIter first_iter = (*points)[cn[0]];
    Vector3D space_pt = first_iter->getPoint();
    PointIter curr_iter = first_iter;
    PointIter next_iter = first_iter->getFirstNeighbour();
    Vector3D curr_pt = curr_iter->getPoint();
    copy(curr_pt.begin(), curr_pt.end(), std::back_inserter(pts));
    Vector3D next_pt = next_iter->getPoint();
    copy(next_pt.begin(), next_pt.end(), std::back_inserter(pts));
    curr_iter = next_iter;
    next_iter = curr_iter->getFirstNeighbour();
    while (curr_iter != first_iter) {
	curr_pt = curr_iter->getPoint();
	copy(curr_pt.begin(), curr_pt.end(), std::back_inserter(pts));
	next_pt = next_iter->getPoint();
	copy(next_pt.begin(), next_pt.end(), std::back_inserter(pts));
	curr_iter = next_iter;
	next_iter = curr_iter->getFirstNeighbour();
    }
    LineCloud lc(pts.begin(), (int)pts.size()/6);
    lc.writeStandardHeader(debug);
    lc.write(debug);
#endif // FANTASTIC_DEBUG

    par.attach(op);
    par.setBiCGTolerance(0.00001);
    par.parametrize();
  } catch(...) {
    status.setError(FT_ERROR_IN_PARAMETERIZE);
    return status;
  }

#ifdef FANTASTIC_DEBUG
  std::ofstream parout("tmp/pardump.dat");
  std::ofstream paredgesout("tmp/paredges.dat");
  points->printUVNodes(parout, true);
  points->printUVEdges(paredgesout);
#endif // FANTASTIC_DEBUG

  // Rescale surf_ based on lengths of outer boundary loop.
  // Parameter values of points are also updated.
  if (true) {
    local_status = reparametrizeSurf(edgeloop, corner, *points);
    for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
      status.addWarning(local_status.getWarning(ki));
    if (!local_status.isOK()) {
      status.setError(local_status.getMessage());
      return status;
    }
  }
  if (true) {
    points->reparBdy(surf_, true); // Second parametrization, we use existing uv.
    points->reparInnerPoints(surf_, true);
  }

  //  std::cout << "Number of data points : " << points->size() << std::endl;

  // Make surface. Parameter iteration and iteration on the
  // weight on approximation compared to smoothness is included.
  int maxiter = 5;
  int max_update_iter = 1;
  local_status = modifySurface(maxiter, max_update_iter, edge_derivs, 
			       *points, max_error, mean_error);
  for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
    status.addWarning(local_status.getWarning(ki));
  if (!local_status.isOK()) {
    status.setError(local_status.getMessage());
    return status;
  }

  // We make sure that corners lie in corners of parameter domain.
  for (ki = 0; ki < (int) cn.size(); ++ki) {
      double umin = surf_->startparam_u();
      double umax = surf_->endparam_u();
      double vmin = surf_->startparam_v();
      double vmax = surf_->endparam_v();
      Vector2D cn_pt = (*points)[cn[ki]]->getPar();
      if (cn_pt[0] != umin && cn_pt[0] != umax) {
	  if (fabs(cn_pt[0] - umin) < fabs(umax - cn_pt[0])) {
	      cn_pt[0] = umin;
	  } else {
	      cn_pt[0] = umax;
	  }
      }
      if (cn_pt[1] != vmin && cn_pt[1] != vmax) {
	  if (fabs(cn_pt[1] - vmin) < fabs(vmax - cn_pt[1])) {
	      cn_pt[1] = vmin;
	  } else {
	      cn_pt[1] = vmax;
	  }
      }
      ((*points)[cn[ki]])->setPar(cn_pt);
  }

  // We have created our surface surf_.
  // We make sure that m_ & n_ correspond to longest/shortest dir.
  if (m_ != n_) { // We may need to swap values.
      vector<double> appr_edge_lengths(4);
      int nmb_samples = 10;
      for (ki = 0; ki < 4; ++ki) {
	  shared_ptr<SplineCurve> ccw_edge_cv(surf_->edgeCurve(ki));
	  appr_edge_lengths[ki] = ccw_edge_cv->ParamCurve::estimatedCurveLength(nmb_samples);
      }
      double avg_length_1 = 0.5*(appr_edge_lengths[0] + appr_edge_lengths[2]);
      double avg_length_2 = 0.5*(appr_edge_lengths[1] + appr_edge_lengths[3]);
      if (avg_length_2 > avg_length_1) {
	  swap(m_, n_);
      }
  }

   // If there exists a neighbour sf which was created prior to this we check grid resolutions
   // along matching edges. If it differs we use that of neighbour. If there exist for both
   // edges in u-/v-direction and they differ we exit (impossible input).
   // @@sbr It does not fix problem with grid sizes for a rotational model.
   vector<pair<int, int> > matching_edges, matching_grid_res;
   vector<ftChartSurface*> matching_faces;
   vector<double> rot_alphas;
   local_status = getMatchingEdges(edgeloop, corner, matching_edges, matching_grid_res,
				   matching_faces, rot_alphas);
   for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
     status.addWarning(local_status.getWarning(ki));
   if (!local_status.isOK()) {
     status.setError(local_status.getMessage());
     return status;
   }

   vector<int> edge_res(4);
   edge_res[0] = edge_res[1] = m_;
   edge_res[2] = edge_res[3] = n_;
   for (ki = 0; ki < (int) matching_edges.size(); ++ki) {
     edge_res[matching_edges[ki].first] = matching_grid_res[ki].second;
   }
   if ((edge_res[0] != edge_res[1]) || (edge_res[0] != m_)) {
     // We must surely alter m_.
     MESSAGE("Altering grid res to match that of neighbour sf...");
     if (edge_res[0] == edge_res[1]) {
       m_ = edge_res[0];
     } else if (edge_res[0] != m_ && edge_res[1] != m_) {
       status.setError(FT_NEIGHBOUR_GRID_SIZE_DIFFER);
     } else {
       m_ = (edge_res[0] != m_) ? edge_res[0] : edge_res[1];
     }
   }
   if ((edge_res[2] != edge_res[3]) || (edge_res[2] != n_)) {
     // We must surely alter n_.
     MESSAGE("Altering grid res to match that of neighbour sf...");
     if (edge_res[2] == edge_res[3]) {
       n_ = edge_res[2];
     } else if (edge_res[2] != n_ && edge_res[3] != n_) {
       status.setError(FT_NEIGHBOUR_GRID_SIZE_DIFFER);
     } else {
       n_ = (edge_res[2] != n_) ? edge_res[2] : edge_res[3];
     }
   }

   // We next construct the graph_.
   points->orderNeighbours(); // As we've picked from points graph_nodes will have consistent topology.
   try {
     createGraph(graph_nodes);
   } catch (...) {
     status.setError(FT_FAILED_CREATING_GRAPH);
     return status;
   }

#if 1
    std::ofstream fileout("tmp/surf_.g2");
    surf_->writeStandardHeader(fileout);
    surf_->write(fileout);
#endif

  // Create the grid distribution functions using local information
  // for this surface
  createGridDistrFunctions();

  // Modify grid distribution functions with respect to neighbouring
  // surfaces
  // @@ As we modify gdf's for each new sf, routine is called more times than needed.
  // @@ If moved outside class we'll fix that. The extra time spent is hardly critical.
  local_status = modifyGridDistrFunctions(bdpiece);
  for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
    status.addWarning(local_status.getWarning(ki));
  if (!local_status.isOK()) {
    status.setError(local_status.getMessage());
    return status;
  }

  return status;
}

//===========================================================================
void ftChartSurface::constructBasises(ftPointSet& pointset,
				      const vector<int>& corners, double epsge,
				      BsplineBasis& basis_u,
				      BsplineBasis& basis_v)
//===========================================================================
{
    pointset.orderNeighbours();

    // For each edge we collect sampled points.
    int dim = 3;
    int nmb_edges = (int)corners.size();
    vector<vector<Vector3D> > space_pts(nmb_edges);
    vector<vector<Vector2D> > par_pts(nmb_edges);
    int i;
    
    for (i = 0; i < nmb_edges; ++i) {
	PointIter pt_iter = pointset[corners[i]];
	space_pts[i].push_back(pt_iter->getPoint());
	par_pts[i].push_back(pt_iter->getPar());
	while (pt_iter != pointset[corners[(i+1)%nmb_edges]]) {
	    pt_iter = pt_iter->getFirstNeighbour();
	    space_pts[i].push_back(pt_iter->getPoint());
	    par_pts[i].push_back(pt_iter->getPar());
	}
    }

    // We must reverse ordering for edges 2 & 3.
    for (i = 2; i < 4; ++i) {
	reverse(space_pts[i].begin(), space_pts[i].end());
	reverse(par_pts[i].begin(), par_pts[i].end());
    }
    // Extract values to vectors.
    vector<vector<double> > points(nmb_edges);
    vector<vector<double> > parvals(nmb_edges);
    for (i = 0; i < nmb_edges; ++i)
	for (size_t j = 0; j < space_pts[i].size(); ++j) {
	    points[i].insert(points[i].end(),
			     space_pts[i][j].begin(), space_pts[i][j].end());
	    // The points are sampled from edges defining a rectangle.
	    parvals[i].push_back(par_pts[i][j][i%2]);
	}

    // Points are collected, time for spline approximation of edge 0 & 1.
    vector<BsplineBasis> basises;
    for (i = 0; i < 2; ++i) {
	ApproxCurve approx_curve(points[i], parvals[i], dim, epsge);
	double maxdist, avdist;
	// We use inital value of 5 iterations.
	shared_ptr<SplineCurve> approx_crv = 
	  approx_curve.getApproxCurve(maxdist, avdist);
	basises.push_back(approx_crv->basis());
    }

    // We then use the basis as a basis for approximating edge 2 & 3.
    for (i = 2; i < 4; ++i) {
	vector<double> knots(basises[i-2].begin(), basises[i-2].end());
	ApproxCurve approx_curve(points[i], parvals[i], dim, epsge,
				   basises[i-2].numCoefs(), basises[i-2].order(),
				   knots);
	double maxdist, avdist;
	shared_ptr<SplineCurve> approx_crv = 
	  approx_curve.getApproxCurve(maxdist, avdist);
	basises[i-2] = approx_crv->basis();
    }

    basis_u = basises[0];
    basis_v = basises[1];
}


//===========================================================================
ftMessage ftChartSurface::createGraph(vector<ftSamplePoint*>& graph_nodes)//ftPointSet& points)
//===========================================================================
{
    ftMessage status;

    try {
	graph_.setGraph(graph_nodes); //points); // @@ Should return ftMessage.
    } catch (...) {
	status.setError(FT_NON_4_SIDED_SURF);
    }

    return status;
}

//===========================================================================
int ftChartSurface::nmbToEval(ftEdge* edge, double tmin, double tmax)
//===========================================================================
{
    int i;
    double max_dist = std::min(curvature_tol_, toptol_.neighbour*0.1); // @@sbr Dividing
    // Based on curvature of edge (distance from sampled points to straight line
    // between end points) we return number of points to be evaluated (>=2).
    // Return value on form 2^n + 1 (we split in two when not within max_dist).
    
    // We test line segment in 20 interior points. If inside, we approve.
    // Otherwise we split edge in equal (in parameter domain) halfs and test again.
    Point start_pt = edge->point(tmin);
    Point end_pt = edge->point(tmax);
    int nmb_test_pts = 40;//20;
    double length = start_pt.dist(end_pt);
    double step = (tmax - tmin) / (nmb_test_pts + 1);
    for (i = 1; i < nmb_test_pts + 1; ++i) {
	Point pt = edge->point(tmin + i*step);
	double length1 = pt.dist(start_pt);
	double length2 = pt.dist(end_pt);	
	// We use Herons formula: A = sqrt(s(s-a)(s-b)(s-c)), s=(a+b+c)/2.
	double s = 0.5*(length + length1 + length2);
	double height = 2*sqrt(s*(s-length)*(s-length1)*(s-length2))/length;

	if (height > max_dist)
	    break;
    }

    if (i < nmb_test_pts + 1) {
	// # of segments = nmbToEval - 1 (we evaluate in both end points).
	int nmb_1 = nmbToEval(edge, tmin, tmin + 0.5*(tmax - tmin)) - 1;
	int nmb_2 = nmbToEval(edge, tmin + 0.5*(tmax - tmin), tmax) - 1;
	return (2*max(nmb_1, nmb_2) + 1); // @@sbr Remove when done debugging!
    } else
	return 2; // No need to (possibly further) split edge.
}

//===========================================================================
void ftChartSurface::addOuterBoundaryPoints(ftPointSet& points)
//===========================================================================
{
    //    return; // @@sbr Remove when done debugging!!!
    points.orderNeighbours();
    // We start by finding a sampled point on the outer boundary.
    ftSamplePoint* first_pt = NULL;
    int ki, kj;
    for (ki = 0; ki < points.size(); ++ki) {
	if (points[ki]->isOnBoundary()) {
	    first_pt = points[ki];
	    break;
	}
    }

    ftSurfaceSetPoint* last_pt = dynamic_cast<ftSurfaceSetPoint*>(first_pt);
    ftSurfaceSetPoint* next_pt = dynamic_cast<ftSurfaceSetPoint*>(first_pt->getFirstNeighbour());
    ASSERT(last_pt != 0 && next_pt != 0);
    shared_ptr<ftSurfaceSetPoint> new_pt;
    Vector3D median_space_pt;
    Vector2D median_par_pt;

    int min_samples = 3;
    int max_samples = 18; // Sample no less than 3 points between two samples points.
    // According to the current routine only outer (four) edges are traversed.
    const int at_bd = 1;
    PointIter second_pt = NULL;
    while (next_pt != first_pt) {
	next_pt->removeNeighbour(last_pt);
	last_pt->removeNeighbour(next_pt);
	ftSurfaceSetPoint* next_next_pt = // In order to avoid unnecessary sorting.
	    dynamic_cast<ftSurfaceSetPoint*>(next_pt->getFirstNeighbour());

	vector<shared_ptr<ftFaceBase> > faces;
	for (ki = 0; ki < last_pt->nmbFaces(); ++ki) {
	    for (kj = 0; kj < next_pt->nmbFaces(); ++kj) {
		if (last_pt->face(ki) == next_pt->face(kj)) {
		    faces.push_back(last_pt->face(ki));
		}
	    }
	}
	ASSERT(faces.size() > 0);

	double length = last_pt->getPoint().dist(next_pt->getPoint());
	double incr = 20.0*toptol_.neighbour;
	
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	int nmb_to_sample = max(min_samples, min(max_samples,
						 (int)(0.25*length/incr)));
	
#else
	int nmb_to_sample =  std::max(min_samples, 
				      std::min(max_samples, (int)(0.25*length/incr)));
#endif

	double weight = 1.0 / (nmb_to_sample + 1);

	for (ki = 1; ki < nmb_to_sample + 1; ++ki) {
	    double step = weight*ki;
	    vector<Vector2D> par_pts;
	    for (kj = 0; kj < (int) faces.size(); ++kj) {
		median_par_pt =
		    (1 - step)*last_pt->getPar(faces[kj].get()) + step*next_pt->getPar(faces[kj].get());
		Vector2D clo_in_dom;
		double tolerance = 1e-08;
		const Domain& face_domain = faces[kj]->surface()->parameterDomain();
		face_domain.closestInDomain(median_par_pt, clo_in_dom, tolerance);
		par_pts.push_back(clo_in_dom);
	    }

	    Point space_pt = faces[0]->surface()->point(par_pts[0][0], par_pts[0][1]);
	    median_space_pt.setValue(space_pt.begin());
	    new_pt = shared_ptr<ftSurfaceSetPoint>
		(new ftSurfaceSetPoint(median_space_pt, at_bd));
	    if (second_pt == 0) {
		second_pt = new_pt.get();
	    }
	    for (kj = 0; kj < (int) faces.size(); ++kj) {
		new_pt->addPair(faces[kj], par_pts[kj]);
	    }
// 	    new_pt->setPar(median_par_pt);
	    if (ki == 1) { // @@ Important that nmb_to_sample > 1!!!
		last_pt->addNeighbour(new_pt.get());
		new_pt->addNeighbour(last_pt);
	    } else if (ki == nmb_to_sample) {
		points.lastAdded()->addNeighbour(new_pt.get());
		new_pt->addNeighbour(points.lastAdded());       
		next_pt->addNeighbour(new_pt.get());
		new_pt->addNeighbour(next_pt);
	    } else {
		points.lastAdded()->addNeighbour(new_pt.get());
		new_pt->addNeighbour(points.lastAdded());
	    }

#ifdef FANTASTIC_DEBUG
	    std::ofstream debug("tmp/debug.g2");
	    vector<double> pts(6);
	    Vector3D from = new_pt->getPoint();
	    copy(from.begin(), from.end(), pts.begin());
	    Vector3D to = new_pt->getFirstNeighbour()->getPoint();
	    copy(to.begin(), to.end(), pts.begin() + 3);
	    LineCloud lc(pts.begin(), 1);
	    lc.writeStandardHeader(debug);
	    lc.write(debug);
#endif // FANTASTIC_DEBUG

	    points.addEntry(new_pt);
	}
	last_pt = next_pt;
	next_pt = next_next_pt;
	ASSERT(next_pt != 0);
    }

    vector<shared_ptr<ftFaceBase> > faces;
    for (ki = 0; ki < last_pt->nmbFaces(); ++ki) {
	for (kj = 0; kj < next_pt->nmbFaces(); ++kj) {
	    if (last_pt->face(ki) == next_pt->face(kj)) {
		faces.push_back(last_pt->face(ki));
	    }
	}
    }
    ASSERT(faces.size() > 0);

    // We add points between the last two points.
    double length = last_pt->getPoint().dist(next_pt->getPoint());
    double incr = 20.0*toptol_.neighbour;
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    int nmb_to_sample= max(min_samples, min(max_samples, (int)(0.25*length/incr)));
#else
    int nmb_to_sample= std::max(min_samples, std::min(max_samples, (int)(0.25*length/incr)));
#endif
    double weight = 1.0 / (nmb_to_sample + 1);
    next_pt->removeNeighbour(last_pt);
    last_pt->removeNeighbour(next_pt);
    for (ki = 1; ki < nmb_to_sample + 1; ++ki) {
	double step = weight*ki;
	vector<Vector2D> par_pts;
	for (kj = 0; kj < (int) faces.size(); ++kj) {
	    median_par_pt =
		(1 - step)*last_pt->getPar(faces[kj].get()) + step*next_pt->getPar(faces[kj].get());
	    Vector2D clo_in_dom;
	    double tolerance = 1e-08;
	    const Domain& face_domain = faces[kj]->surface()->parameterDomain();
	    face_domain.closestInDomain(median_par_pt, clo_in_dom, tolerance);
	    par_pts.push_back(clo_in_dom);
	}

	Point space_pt = faces[0]->surface()->point(par_pts[0][0], par_pts[0][1]);
	median_space_pt.setValue(space_pt.begin());
	new_pt = shared_ptr<ftSurfaceSetPoint>
	    (new ftSurfaceSetPoint(median_space_pt, at_bd));
	for (kj = 0; kj < (int) faces.size(); ++kj) {
	    new_pt->addPair(faces[kj], par_pts[kj]);
	}
	// 	    new_pt->setPar(median_par_pt);
	if (ki == 1) { // @@ Important that nmb_to_sample > 1!!!
	    last_pt->addNeighbour(new_pt.get());
	    new_pt->addNeighbour(last_pt);
	} else if (ki == nmb_to_sample) {
	    points.lastAdded()->addNeighbour(new_pt.get());
	    new_pt->addNeighbour(points.lastAdded());       
	    next_pt->addNeighbour(new_pt.get());
	    new_pt->addNeighbour(next_pt);
	} else {
	    points.lastAdded()->addNeighbour(new_pt.get());
	    new_pt->addNeighbour(points.lastAdded());
	}
	points.addEntry(new_pt);
    }

    points.setFirst(first_pt);
    points.setSecond(second_pt);
}

   
//===========================================================================
void ftChartSurface::createGridDistrFunctions()
//===========================================================================
{
   ASSERT(surf_.get() != 0);
   ASSERT(grid_distr_functions_.size() == 4);

   // As a first implementation we let the functions be the identity mappings.
   // Later on they will probably be hermite interpolants.
   // f: [tmin, tmax] -> [tmin, tmax].

   // Fetch endparameters
   double uvpars[4], par[4];
   uvpars[0] = surf_->startparam_u();
   uvpars[1] = surf_->endparam_u();
   uvpars[2] = surf_->startparam_v();
   uvpars[3] = surf_->endparam_v();
   int ki, kj, kh;
   for (ki=0; ki<2; ki++)
     {
       for (kj=0; kj<2; kj++)
	 {
	   par[1-ki] = par[3-ki] = uvpars[2*(1-ki)+kj];
	   for (kh=0; kh<2; kh++)
	     par[ki+2*kh] = uvpars[2*ki+kh];
	   grid_distr_functions_[2*ki+kj] = 
	     createOneDistrFunction(par, par+2, ki==0);
	   secn_distr_[2*ki+kj] = true;
	 }
     }

   if (symm_distr_functions_)
     {
       // Make opposite distribution functions equal
       adaptGridDistrFunc(grid_distr_functions_[0].get(),
			  grid_distr_functions_[1].get(), false);
       adaptGridDistrFunc(grid_distr_functions_[2].get(),
			  grid_distr_functions_[3].get(), false);
     }
       

   // Modify the distribution functions according to edge scale
   vector<double> scale = getEdgeScales();

   if (scale[0] != 1.0)
     {
       // Bottom side. Modify start of left and right distribution function
       // at the start
       modifyOneDistrFunction(grid_distr_functions_[2], scale[0], true, n_-1);
       modifyOneDistrFunction(grid_distr_functions_[3], scale[0], true, n_-1);
     }

   if (scale[1] != 1.0)
     {
       // Top side. Modify start of left and right distribution function
       // at the end
       modifyOneDistrFunction(grid_distr_functions_[2], scale[1], false, n_-1);
       modifyOneDistrFunction(grid_distr_functions_[3], scale[1], false, n_-1);
     }

   if (scale[2] != 1.0)
     {
       // Left side. Modify start of bottom and top distribution function
       // at the start
       modifyOneDistrFunction(grid_distr_functions_[0], scale[2], true, m_-1);
       modifyOneDistrFunction(grid_distr_functions_[1], scale[2], true, m_-1);
     }

   if (scale[3] != 1.0)
     {
       // Right side. Modify start of bottom and top distribution function
       // at the end
       modifyOneDistrFunction(grid_distr_functions_[0], scale[3], false, m_-1);
       modifyOneDistrFunction(grid_distr_functions_[1], scale[3], false, m_-1);
     }

   
}

//===========================================================================
shared_ptr<SplineCurve> 
ftChartSurface::createOneDistrFunction(double par1[], double par2[],
				       bool isudir)
//===========================================================================
{
  // Evaluate the surface in the endpoints of the current boundary curve
  vector<Point> bd1 = surf_->ParamSurface::point(par1[0], par1[1], 2);
  vector<Point> bd2 = surf_->ParamSurface::point(par2[0], par2[1], 2);
  
  // Fetch the derivatives along the current boundary
  int idx = (isudir) ? 0 : 1;
  vector<Point> der1(3), der2(3);
  der1[0] = bd1[0];
  der1[1] =  bd1[1+idx];
  der1[2] = bd1[3+2*idx];
  der2[0] = bd2[0];
  der2[1] =  bd2[1+idx];
  der2[2] = bd2[3+2*idx];

  // Get information about tangent lengts and parameter interval
  double parint, len1, len2;
  getHermiteData(der1, der2, parint, len1, len2);

  // Modify the tangent lengths in such aq way that
  // len1 + len2 = 2*parint/3
  double fac = 3.0*(len1 + len2)/(2.0*parint);
  len1 /= fac;
  len2 /= fac;

  // Modify the Hermite information to make them fit to the
  // current parameter interval, and make Hermite curve segment
//    double parvals[4], int_cond[4];
//    parvals[0] = parvals[1] = par1[idx];
//    parvals[2] = parvals[3] = par2[idx];
//    int_cond[0] = par1[idx];
//    int_cond[1] = len1*(par2[idx]-par1[idx])/parint;
//    int_cond[2] = par2[idx];
//    int_cond[3] = len2*(par2[idx]-par1[idx])/parint;

//    HermiteInterpolator hermite;
//    hermite.interpolate(4, 1, parvals, int_cond, coefs);

  // Make spline curve
  vector<double> knots(8);
  int ki;
  for (ki=0; ki<4; ki++)
    {
      knots[ki] = par1[idx];
      knots[4+ki] = par2[idx];
    }
  vector<double> coefs(4);
  coefs[0] = par1[idx];
  coefs[1] = par1[idx] + len1*(par2[idx]-par1[idx])/parint;
  coefs[2] = par2[idx] - len2*(par2[idx]-par1[idx])/parint;
  coefs[3] = par2[idx];
  shared_ptr<SplineCurve> bdcv = shared_ptr<SplineCurve>
    (new SplineCurve(4, 4, &knots[0], &coefs[0], 1));

//    // Make degree 5 curve
//    bdcv->raiseOrder(2);

  return bdcv;
}

//===========================================================================
void ftChartSurface::modifyOneDistrFunction(shared_ptr<SplineCurve> func,
					    double scale, bool atStart, 
					    int nbints)
//===========================================================================
{
  // Define interpolation problems
  double aint = func->endparam() - func->startparam();
  double ta1 = func->startparam() + aint/(double)nbints;
  double ta2 = func->endparam() - aint/(double)nbints;

  Point pt1, pt2;
  Point p1, p2;
  p1 = func->ParamCurve::point(func->startparam());
  p2 = func->ParamCurve::point(func->endparam());

  pt1 = func->ParamCurve::point(ta1);
  pt2 = func->ParamCurve::point(ta2);

  vector<double> param;
  vector<double> data;
    
  SplineInterpolator interpolate;
  interpolate.setFreeConditions();

  // Start condition
  param.push_back(func->startparam());
  data.push_back(p1[0]);

  // Second condition
  param.push_back(ta1);
  data.push_back(atStart ? p1[0]+scale*(pt1[0]-p1[0]) : pt1[0]);

  // Third condition
  param.push_back(ta2);
  data.push_back(atStart ? pt2[0] : p2[0]-scale*(p2[0]-pt2[0]));

  // Last condition
  param.push_back(func->endparam());
  data.push_back(p2[0]);

  // Solve interpolation problems
  vector<double> coefs;
  interpolate.interpolate((int)param.size(), 1, &param[0], &data[0], coefs);

  // Copy result to distribution function
  int ki;
  vector<double>::iterator c1 = func->coefs_begin();
  for (ki=0; ki<(int)coefs.size(); ki++)
    c1[ki] = coefs[ki];
 }


// //===========================================================================
// void ftChartSurface::testLocateInGraph(ftPointSet& points)
// //===========================================================================
// {
//     vector<shared_ptr<ftFaceBase> > faces; // For debugging

//     for (int i = 0; i < points.size(); ++i) {
// 	ftSamplePoint* sampled_pt = points[i];
// 	Vector3D vector = sampled_pt->getPoint();
// 	Point pt(vector[0], vector[1], vector[2]);
// 	// We do a closest point on created surface.
// 	double clo_u, clo_v, clo_dist;
// 	Point clo_pt;
// 	double epsilon = 1e-10;
// 	surf_->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon);

// 	// We then locate corresponding patch and estimated local parameter values.
// 	shared_ptr<ftFaceBase> face;
// 	graph_.getLocalParameters(clo_u, clo_v, face);
// 	Point guessed_pt = face->Surface()->ParamSurface::point(clo_u, clo_v);

// 	// for debugging
// 	faces.push_back(face);
// 	// end of debugging

// // 	Vector2D old_par =
// // 	    (static_cast<ftSurfaceSetPoint*>(sampled_pt))->getPar(face);
// // 	Vector2D new_par(clo_u, clo_v);

// 	double guessed_dist = pt.dist(guessed_pt); //old_par.dist(new_par);
// // 	cout << "Param dist for point number " << i << ": " << guessed_dist << endl;
//     }
// }

//===========================================================================
ftMessage 
ftChartSurface::getMatchingEdges(Point from_space_pt, Point to_space_pt,
				 ftFaceBase* face1, ftFaceBase* face2,
				 int& edge1, double& par1, 
				 double& par2, int& edge2,
				 double& par3, double& par4,
				 double epsgeo)
//===========================================================================
{
   ftMessage status;

   shared_ptr<SplineSurface> sf1 =
      dynamic_pointer_cast<SplineSurface, ParamSurface>(face1->surface());
   shared_ptr<SplineSurface> sf2 =
      dynamic_pointer_cast<SplineSurface, ParamSurface>(face2->surface());
   if ((sf1.get() == 0) || (sf2.get() == 0)) {
       status.setError(FT_NOT_SPLINE_SURF);
       return status;
   }

//    double space_eps = 1e-05;
   //   double par1, par2;
   sf1->getBoundaryIdx(from_space_pt, to_space_pt,
		       epsgeo, edge1, par1, par2);
   sf2->getBoundaryIdx(from_space_pt, to_space_pt,
		       epsgeo, edge2, par3, par4);

   return status;
}

//===========================================================================
ftMessage 
ftChartSurface::modifyGridDistrFunctions(vector<BoundaryPiece>& bdpiece)
  // Modify the grid distribution functions to get smooth
  // transitions across block boundaries and correspondence
  // along block boundaries
//===========================================================================
{
  ftMessage status;
  int ki, kj;
  double par1, par2;
  double param1[2], param2[2]; //, param3[2], param4[2];
  Point pt1, pt2;
  shared_ptr<SplineSurface> currsf; 
  ftSurfaceSet *currsfset;
  ftChartSurface *currchart;
  ftFaceBase *curr_face;
//   bool oppisite;
  int curridx, idx;
  for (ki=0; ki<(int)bdpiece.size(); ki++) {
    curr_face = bdpiece[ki].second.first;
    currsfset = dynamic_cast<ftSurfaceSet*>(curr_face);
    if (currsfset)
      currsf = dynamic_pointer_cast<SplineSurface, ParamSurface>
	(currsfset->surface());
      if (currsf.get() == 0) {
	continue;  // Not possible to make a boundary match condition.
      }
      currchart = dynamic_cast<ftChartSurface*>(currsfset);
      if (currchart == 0) {
	continue;
      }

      // Evaluate endpoints of the matching interval on the other
      // surface.
      curridx = (bdpiece[ki]).second.second;
      idx = (curridx > 1) ? 1 : 0;
      param1[idx] = (bdpiece[ki]).first.first;
      param1[1-idx] = ((bdpiece[ki]).second.second % 2 == 0) ?
	currsf->basis(1-idx).startparam() : currsf->basis(1-idx).endparam();
      param2[idx] = (bdpiece[ki]).first.second;
      param2[1-idx] = param1[1-idx];
      currsf->point(pt1, param1[0], param1[1]);
      currsf->point(pt2, param2[0], param2[1]);

      // Find the matching interval on this surface.
      double epsgeo = toptol_.gap;
      int bdidx = -1;
      surf_->getBoundaryIdx(pt1, pt2, epsgeo,
			    bdidx, par1, par2);
      if (bdidx == -1) {
	  status.setError(FT_ERROR_IN_SURFACE_GRIDDING);
	  return status;
      }

      // To avoid long recursions to adapt grid distribution functions
      // just adapt this function to that of the adjacent surface.
      SplineCurve *dfunc1 = 
	getGridDistrFunc(bdidx, par1, par2, epsgeo);
      SplineCurve *dfunc2 =
	currchart->getGridDistrFunc(curridx, param1[idx], param2[idx], epsgeo);
      if (dfunc1 == 0 || dfunc2 == 0) {
	continue;  // Only total match is considered
      }
      bool opposite = ((bdidx == curridx) ||
		       (abs(bdidx - curridx) == 3) ||
		       (bdidx == 1 && curridx == 2) ||
		       (bdidx == 2 && curridx == 1));
      adaptGridDistrFunc(dfunc2, dfunc1, opposite);

//       int nmb = 2;
//       vector<bool> start(nmb, false); //, start2;
//       vector<double> length(nmb); //, length2;
//       vector<bool> secn(nmb, false); //, secn2;
//       vector<int> nbints(nmb); //, nbints2;
//       vector<SplineCurve*> adjfunc(nmb);
//       bool local_bool1, local_bool2;
//       adjfunc[0] = getPriorDistrFunc(bdidx, opposite, local_bool1, length[0], local_bool2, nbints[0]);
//       start[0] = local_bool1;
//       secn[0] = local_bool2;
//       adjfunc[1] = currchart->getPriorDistrFunc(curridx, opposite,
// 						local_bool1, length[1], local_bool2, nbints[1]);
//       start[1] = local_bool1;
//       secn[1] = local_bool2;
//       averageGridDistrFunc(adjfunc, start, length, secn, nbints);
// // 			   adjfunc2, start2, length2, secn2, nbints2);

//       adjfunc[0] = getPriorDistrFunc(bdidx, !opposite, local_bool1, length[0], local_bool2, nbints[0]);
//       start[0] = local_bool1;
//       secn[0] = local_bool2;
//       adjfunc[1] = currchart->getPriorDistrFunc(curridx, !opposite,
// 						local_bool1, length[1], local_bool2, nbints[1]);
//       start[1] = local_bool1;
//       secn[1] = local_bool2;
//       averageGridDistrFunc(adjfunc, start, length, secn, nbints);
// // 			   adjfunc2, start2, length2, secn2, nbints2);
  }

  // We then try to make sure that size of blocks in corner match.
  // Given input piece, we must extract top edge, and then the edges surrounding
  // If not all faces have been created, we do not bother to rescale cvs.
  ftEdgeBase* first_edge = boundary_loops_[0]->getEdge(0).get();
  ftEdgeBase* prev_edge = first_edge;
  ftEdgeBase* curr_edge = prev_edge->next();
  bool finished = false;
  while (!finished) { // We continue if in a corner.
      if (curr_edge == first_edge) {
	  finished = true; // This will be our last lap.
      }
      // If start pt of edge is in additional_corner_pts_ we're in a corner.
      Point start_pt = curr_edge->point(curr_edge->tMin());
      for (ki = 0; ki < (int)additional_corner_pts_.size(); ++ki) {
	double dist = start_pt.dist(additional_corner_pts_[ki]);
	if (dist < toptol_.gap) {
	  break;
	}
      }
      if (ki == (int)additional_corner_pts_.size()) { // start_pt not in additional_corner_pts_.
	tpJointType cont = prev_edge->checkContinuity(curr_edge,
						      toptol_.neighbour, toptol_.gap,
						      toptol_.bend, toptol_.kink);
	if (cont < JOINT_G0) {
	  prev_edge = curr_edge;
	  curr_edge = curr_edge->next();
	  continue;
	}
      }

      vector<ftEdgeBase*> edg;
      vector<bool> edg_start;
      curr_edge->adjacentEdges(true, edg, edg_start);
      edg.insert(edg.begin(), curr_edge);
      edg_start.insert(edg_start.begin(), true);
      // We run through edg, moving on to next corner if not all sfs were constructed.
      for (ki = 0; ki < (int)edg.size(); ++ki) {
	  shared_ptr<ParamSurface> sf(edg[ki]->face()->surface());
	  if (sf.get() == 0) {
	      break;
	  }
      }
      if (ki < (int)edg.size()) {
	  prev_edge = curr_edge;
	  curr_edge = curr_edge->next();
	  continue;
      }
      // We make sure that edges are given in the cw direction.
//       double angle_tol = toptol_.kink;
//       cmUtils::cwOrientation(edg, edg_start, angle_tol);
      cmUtils::cwOrientation2(edg, edg_start);
      if (edg.size() == 2) {
	  prev_edge = curr_edge;
	  curr_edge = curr_edge->next();
	  continue;
      }

      // If 2*odd #edges meet we must divide into 2 groups.
      vector<bool> start;
      vector<double> length;
      vector<bool> secn;
      vector<int> nbints;
      vector<SplineCurve*> gdfs; // Grid distribution functions.
      for (ki = 0; ki < (int)edg.size(); ++ki) {
	  ftChartSurface* other_chart = dynamic_cast<ftChartSurface*>(edg[ki]->face());
	  ASSERT(other_chart != 0);
	  // Using end pts we decide along which edge (b, t, l, r) we are.
	  shared_ptr<SplineSurface> face_sf =
	      dynamic_pointer_cast<SplineSurface, ParamSurface>(other_chart->surface());
	  if (face_sf.get() == 0) {
	      break;
	  }
	  Point mid_pt = edg[ki]->point(0.5*(edg[ki]->tMin() + edg[ki]->tMax()));
	  double epsge = 10*toptol_.gap; // As it should be in the middle of an edge of spline sf
	                                 // (not deg) we may use large epsilon.
	  int bd_idx; // (b, t, l, r)
	  face_sf->getBoundaryIdx(mid_pt, epsge, bd_idx);
	  if (bd_idx == -1) {
#ifdef FANTASTIC_DEBUG
	      std::ofstream debug("tmp/debug.g2");
	      SplineCurve* space_cv = edg[ki]->geomEdge()->geomCurve()->geometryCurve();
	      if (space_cv != 0) {
		  space_cv->writeStandardHeader(debug);
		  space_cv->write(debug);
	      }
	      face_sf->writeStandardHeader(debug);
	      face_sf->write(debug);
#endif // FANTASTIC_DEBUG
	      THROW("Failed finding edge within tol!");
	  }
	  double local_length;
	  bool local_start, local_secn;
	  int local_nbints;
	  // Input bd_idx refers to that of a neighbour edge.
	  int other_bdidx;
	  bool opposite;
	  if (bd_idx == 0) {
	      other_bdidx = (edg_start[ki]) ? 2 : 3;
	      opposite = false;
	  } else if (bd_idx == 1) {
	      other_bdidx = (edg_start[ki]) ? 3 : 2;
	      opposite = true;
	  } else if (bd_idx == 2) {
	      other_bdidx = (edg_start[ki]) ? 1 : 0;
	      opposite = false;
	  } else if (bd_idx == 3) {
	      other_bdidx = (edg_start[ki]) ? 0 : 1;
	      opposite = true;
	  }
	  gdfs.push_back(other_chart->getPriorDistrFunc(other_bdidx, opposite, local_start,
							local_length, local_secn, local_nbints));
	  start.push_back(local_start);
	  length.push_back(local_length);
	  secn.push_back(local_secn);
	  nbints.push_back(local_nbints);
      }
      if (ki == (int)edg.size()) { // If all faces meeting in corner were made we continue.
	  // We must then fix the gdfs.
	  // @@sbr We may need more advanced algorithm (what if common pt lies on an outer edge?)
	  if ((gdfs.size())%4 == 0) {
	      // Based on tangent in start_pt we decide how to split cvs into two groups.
	      for (ki = 0; ki < 2; ++ki) {
		  vector<bool> start1;
		  vector<double> length1;
		  vector<bool> secn1;
		  vector<int> nbints1;
		  vector<SplineCurve*> gdfs1;
		  for (kj = 2*ki; kj < (int)gdfs.size(); kj += 4) {
		      gdfs1.push_back(gdfs[kj]); // Using the fact that edg is sorted based on angles.
		      start1.push_back(start[kj]);
		      length1.push_back(length[kj]);
		      secn1.push_back(secn[kj]);
		      nbints1.push_back(nbints[kj]);
		      gdfs1.push_back(gdfs[kj+1]); // Using the fact that edg is sorted based on angles.
		      start1.push_back(start[kj+1]);
		      length1.push_back(length[kj+1]);
		      secn1.push_back(secn[kj+1]);
		      nbints1.push_back(nbints[kj+1]);
		  }
 		  averageGridDistrFunc(gdfs1, start1, length1, secn1, nbints1);
// 		  averageGridDistrFunc(gdfs1[0], start1[0], length1[0], secn1[0], nbints1[0],
// 				       gdfs1[1], start1[1], length1[1], secn1[1], nbints1[1]);
	      }
	  } else {
	      averageGridDistrFunc(gdfs, start, length, secn, nbints);
	  }
      }

      prev_edge = curr_edge;
      curr_edge = curr_edge->next();
  }

  return status;
}

//===========================================================================
ftMessage
ftChartSurface::modifyGridDistrFunctions(RotationInfo* rot_info, ftCurve* total_outer_bd)
//===========================================================================
{
  int ki;
    ftMessage status;
    if (rot_info == 0) {
	return status;
    }
    ASSERT(total_outer_bd != 0);

   vector<ftEdgeBase*> outer_loop;
   vector<int> corners;
   double bend_tol = toptol_.bend;
   double kink_tol = toptol_.kink;
   ftMessage local_status = 
     getOuterLoop(outer_loop, corners, bend_tol, kink_tol);
   if (!local_status.isOK()) {
      status.setError(local_status.getMessage());
      return status;
   }

    // We must extract bd cvs from outer_bd. If rotation results in matching
    // edges, we must then modify corresponding bd cvs in both end pts of common edge.
    vector<pair<int, int> > matching_edges;
    vector<ftChartSurface*> matching_faces;
    vector<double> rot_alphas;
    local_status = getMatchingEdges(outer_loop, corners, rot_info, total_outer_bd,
				    matching_edges, matching_faces, rot_alphas);
    for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
      status.addWarning(local_status.getWarning(ki));
    if (!local_status.isOK()) {
      status.setError(local_status.getMessage());
      return status;
    }
    double epsgeo = toptol_.gap;
    for (ki = 0; ki < (int)matching_edges.size(); ++ki) {
	// To avoid long recursions to adapt grid distribution functions
	// just adapt this function to that of the adjacent surface.
	double ta = (matching_edges[ki].first < 2) ? surf_->startparam_u() : surf_->startparam_v();
	double tb = (matching_edges[ki].first < 2) ? surf_->endparam_u() : surf_->endparam_v();
	SplineCurve *dfunc1 =
	    getGridDistrFunc(matching_edges[ki].first, ta, tb, epsgeo);
	shared_ptr<SplineSurface> other_sf =
	    dynamic_pointer_cast<SplineSurface, ParamSurface>(matching_faces[ki]->surface());
	ASSERT(other_sf.get() != 0);
	ta = (matching_edges[ki].second < 2) ? other_sf->startparam_u() : other_sf->startparam_v();
	tb = (matching_edges[ki].second < 2) ? other_sf->endparam_u() : other_sf->endparam_v();
	SplineCurve *dfunc2 =
	    matching_faces[ki]->getGridDistrFunc(matching_edges[ki].second, ta, tb, epsgeo);
	if (dfunc1 == 0 || dfunc2 == 0)
	    continue;  // Only total match is considered

	bool opposite = ((matching_edges[ki].first == matching_edges[ki].second) ||
			 (abs(matching_edges[ki].first - matching_edges[ki].second) == 3) ||
			 (matching_edges[ki].first == 1 && matching_edges[ki].second == 2) ||
			 (matching_edges[ki].first == 2 && matching_edges[ki].second == 1));
	adaptGridDistrFunc(dfunc2, dfunc1, opposite);

	bool start1, start2;
	double length1, length2;
	bool secn1, secn2;
	int nbints1, nbints2;
	SplineCurve *adjfunc1 = getPriorDistrFunc(matching_edges[ki].first, 
						  opposite, start1, length1,
						  secn1, nbints1);
	SplineCurve *adjfunc2 =
	    matching_faces[ki]->getPriorDistrFunc(matching_edges[ki].second, 
						  opposite, start2, length2,
						  secn2, nbints2);
	averageGridDistrFunc(adjfunc1, start1, length1, secn1, nbints1,
			     adjfunc2, start2, length2, secn2, nbints2);

	adjfunc1 = getPriorDistrFunc(matching_edges[ki].first, !opposite, 
				     start1, length1, secn1, nbints1);
	adjfunc2 = 
	  matching_faces[ki]->getPriorDistrFunc(matching_edges[ki].second, 
						!opposite, start2, length2,
						secn2, nbints2);
	averageGridDistrFunc(adjfunc1, start1, length1, secn1, nbints1,
			     adjfunc2, start2, length2, secn2, nbints2);
    }

    return status;
}

//===========================================================================
ftMessage 
ftChartSurface::updateGridDistribution(RotationInfo* rot_info, ftCurve* total_outer_bd)
  // Last update of the grid distribution functions to get correspondance
  // along common boundaries
//===========================================================================
{
  ftMessage status;
  int ki, kj;

  // We should then see if any neighbouring sfs of the same type already exist.
  // In that case average common grid distribution function
   vector<ftEdgeBase*> outer_loop;
   vector<int> corners;
   double bend_tol = toptol_.bend;
   double kink_tol = toptol_.kink;
   ftMessage local_status = 
     getOuterLoop(outer_loop, corners, bend_tol, kink_tol);
   if (!local_status.isOK()) {
      status.setError(local_status.getMessage());
      return status;
   }

   // For each of the edges we see whether one other surface has corners 
   // matching edge.
   // We store index of matching edges.
   vector<pair<SplineCurve*, SplineCurve*> > matching_spline_cvs;
//    vector<double> rot_angles;
   vector<bool> opposites;
   for (ki = 0; ki < (int)corners.size() - 1; ++ki) {
      ftEdgeBase* curr_edge = outer_loop[corners[ki]];
      ftEdgeBase* twin_edge = curr_edge->twin();
      ftChartSurface* neighbour_face = (twin_edge != 0) 
	? dynamic_cast<ftChartSurface*>(twin_edge->face()) : 0;
      if (neighbour_face == 0 || neighbour_face->surface().get() == 0)
	 continue;

      ftEdgeBase* last_edge = outer_loop[(corners[ki+1]-1)%(corners.size())];
      ftSSfEdge* s_edge = dynamic_cast<ftSSfEdge*>(curr_edge);      
      if (s_edge == 0)
	{
	  status.setError(FT_UNEXPECTED_DATA_TYPE);
	  return status;
	}

      double tmin = s_edge->tMin();
      Point ptmin = s_edge->point(tmin);

      s_edge = dynamic_cast<ftSSfEdge*>(last_edge);      
      if (s_edge == 0)
	{
	  status.setError(FT_UNEXPECTED_DATA_TYPE);
	  return status;
	}

      double tmax = s_edge->tMax();
      Point ptmax = s_edge->point(tmax); 

      double par1, par2, par3, par4;
      int edge1 = -1, edge2 = -1;
      double epsgeo = toptol_.gap;
      local_status = getMatchingEdges(ptmin, ptmax, this, neighbour_face,
				      edge1, par1, par2,
				      edge2, par3, par4, epsgeo);
      for (kj = 0; kj < local_status.noOfWarnings(); ++kj)
	  status.addWarning(local_status.getWarning(kj));
      if (!local_status.isOK()) {
	  status.setError(local_status.getMessage());
	  return status;
      }

      if (edge1 != -1 && edge2 != -1) {
	  bool opposite = ((par2 - par1)*(par4 - par3) < 0.0);
	  opposites.push_back(opposite);
	  SplineCurve *dfunc1 = 
	      getGridDistrFunc(edge1, par1, par2, epsgeo); // Assuming par domain reflects geometry.
	  SplineCurve *dfunc2 =
	      neighbour_face->getGridDistrFunc(edge2, par3, par4, epsgeo);
	  if (dfunc1 == 0 || dfunc2 == 0)
	      continue;  // Only total match is considered
	  matching_spline_cvs.push_back(make_pair(dfunc1, dfunc2));
	  //       rot_angles.push_back(0.0);
      } else {
	  MESSAGE("Should these edges match?");
      }
   }

   if (rot_info != 0) { // Must check to see if rotation results in further match along edges.
       double epsgeo = toptol_.gap;
       ASSERT(total_outer_bd != 0);
       vector<pair<int, int> > matching_edges;
       vector<ftChartSurface*> matching_faces;
       vector<double> rot_alphas;
       local_status = getMatchingEdges(outer_loop, corners, rot_info, total_outer_bd,
				       matching_edges, matching_faces, rot_alphas);
       for (ki = 0; ki < local_status.noOfWarnings(); ++ki)
	 status.addWarning(local_status.getWarning(ki));
       if (!local_status.isOK()) {
	 status.setError(local_status.getMessage());
	 return status;
       }

       for (ki = 0; ki < (int)matching_edges.size(); ++ki) {
	   double ta = (matching_edges[ki].first < 2) ? surf_->startparam_u() : surf_->startparam_v();
	   double tb = (matching_edges[ki].first < 2) ? surf_->endparam_u() : surf_->endparam_v();
	   SplineCurve *dfunc1 = 
	       getGridDistrFunc(matching_edges[ki].first, ta, tb, epsgeo);
	   shared_ptr<SplineSurface> other_sf =
	       dynamic_pointer_cast<SplineSurface, ParamSurface>(matching_faces[ki]->surface());
	   ASSERT(other_sf.get() != 0);
	   ta = (matching_edges[ki].second < 2) ? other_sf->startparam_u() : other_sf->startparam_v();
	   tb = (matching_edges[ki].second < 2) ? other_sf->endparam_u() : other_sf->endparam_v();
	   SplineCurve *dfunc2 =
	       matching_faces[ki]->getGridDistrFunc(matching_edges[ki].second, ta, tb, epsgeo);
	   if (dfunc1 == 0 || dfunc2 == 0)
	       continue;  // Only total match is considered
	   matching_spline_cvs.push_back(make_pair(dfunc1, dfunc2));
// 	   rot_angles.push_back(rot_alphas[ki]);
	   bool opposite = ((matching_edges[ki].first == matching_edges[ki].second) ||
			    (abs(matching_edges[ki].first - matching_edges[ki].second) == 3) ||
			    (matching_edges[ki].first == 1 && matching_edges[ki].second == 2) ||
			    (matching_edges[ki].first == 2 && matching_edges[ki].second == 1));
	   opposites.push_back(opposite);
       }
   }

   for (ki = 0; ki < (int)matching_spline_cvs.size(); ++ki) {
       averageAdjGridDistrFunc(matching_spline_cvs[ki].second,
			       matching_spline_cvs[ki].first, opposites[ki]);
   }

   return status;
}      
      
//===========================================================================
void ftChartSurface::adaptGridDistrFunc(SplineCurve *from_func,
					SplineCurve *to_func, bool opposite)
//===========================================================================
{
  // Adjusted for parameter interval and orientation, the grid distribution
  // to_func is set to be identical to the function from_func. We modify
  // the coefficients of the spline curve directly to change the grid
  // distribution function itself.

  // First make a copy of from_func and represent this with the same
  // orientation and length as to_func.
#ifdef _MSC_VER
  SplineCurve *from_func2 = dynamic_cast<SplineCurve*>(from_func->clone());
#else
  SplineCurve *from_func2 = from_func->clone();
#endif
  vector<double>::iterator c1 = from_func2->coefs_begin();
  vector<double>::iterator c2 = from_func2->coefs_end();
  if (opposite)
    {
      double del1 = c1[1] - c1[0];
      double del2 = c2[-1] - c2[-2];
      c1[1] = c1[0] + del2;
      c2[-2] = c2[-1] - del1;
    }
 
  // Exchange the coefficients of to_func
  vector<double>::iterator c3 = to_func->coefs_begin();
  vector<double>::iterator c4 = to_func->coefs_end();
  double ta1 = c1[0];
  double ta3 = c3[0];
  double frac = (c4[-1]-c3[0])/(c2[-1]-c1[0]);
  for (; c1<c2; c1++, c3++)
    *c3 = ta3 + (*c1 - ta1)*frac;
}

//===========================================================================
void ftChartSurface::averageAdjGridDistrFunc(SplineCurve *func1,
					     SplineCurve *func2, bool opposite)
//===========================================================================
{
  // Adjusted for parameter interval and orientation, the grid distribution
  // func1 and func2 is set equal.

  // First take the orientation of the correponding surface boundaries
  // into account
  vector<double>::iterator c1 = func1->coefs_begin();
  vector<double>::iterator c2 = func1->coefs_end();
  if (opposite)
    {
      double del[5];  // At most degree 5 curve
      int ki;
      int nn = func1->numCoefs();
      for (ki=1; ki<nn; ki++)
	del[ki-1] = c1[ki] - c1[ki-1];
      for (ki=nn-2; ki>0; ki--)
	c1[ki] = c1[ki+1] - del[nn-ki-2];
    }
 
  // Make average coefficient vector
  vector<double> coef(func1->numCoefs());

  vector<double>::iterator c3 = func2->coefs_begin();
  vector<double>::iterator c4 = func2->coefs_end();
  vector<double>::iterator c5 = coef.begin();

  for (; c1<c2; c1++, c3++, c5++)
    *c5 = 0.5*(*c1 + *c3);

  // Exchange the coefficients 
  c1 = func1->coefs_begin();
  c3 = func2->coefs_begin();
  c5 = coef.begin();

  double ta1 = c1[0];
  double ta3 = c3[0];
  double ta5 = c5[0];
  double frac1 = (c2[-1]-c1[0])/(c5[coef.size()-1]-c5[0]);
  double frac2 = (c4[-1]-c3[0])/(c5[coef.size()-1]-c5[0]);
  for (; c1<c2; c1++, c3++, c5++)
    {
      *c1 = ta1 + (*c5 - ta5)*frac1;
      *c3 = ta3 + (*c5 - ta5)*frac2;
    }

  c1 = func1->coefs_begin();
  if (opposite)
    {
      double del[5];  // At most degree 5 curve
      int ki;
      int nn = func1->numCoefs();
      for (ki=1; ki<nn; ki++)
	del[ki-1] = c1[ki] - c1[ki-1];
      for (ki=nn-2; ki>0; ki--)
	c1[ki] = c1[ki+1] - del[nn-ki-2];
    }
}

//  //===========================================================================
//  void ftChartSurface::averageGridDistrFunc(SplineCurve *func1, bool start1,
//  					  double l1, bool secn1, int nbints1,
//  					  SplineCurve *func2, bool start2,
//  					  double l2, bool secn2, int nbints2)
//  //===========================================================================
//  {
//      // co[ki] +/- del must stay inside end values! Should also be monotone.
//      double length1 = fabs(func1->coefs_begin()[0] - func1->coefs_end()[-1]);
//      double length2 = fabs(func2->coefs_begin()[0] - func2->coefs_end()[-1]);
//      double frac1 = l1/length1; //length1/l1; //;
//      double frac2 = l2/length2; //length2/l2; //;
//      double min_length = min(length1/frac1, length2/frac2);

//    double tmp1, tmp2;
//    if (start1 && start2)
//      {
//        vector<double>::iterator c1 = func1->coefs_begin();
//        vector<double>::iterator c2 = func2->coefs_begin();
//        double del = 0.5*((c1[1]-c1[0])/frac1 + (c2[1]-c2[0])/frac2);
//        if (del > min_length*0.8) {
//  	  del = min_length*0.8;
//        }
//        double del2 = 0.5*((c1[2]-c1[0]) + (c2[2]-c2[0]));
//        double d1 = 0.1*(c1[3]-c1[2]);
//        double d2 = 0.1*(c2[3]-c2[2]);
//        c1[1] = c1[0] + del*frac1;
//        c2[1] = c2[0] + del*frac2;
//        if (secn1)
//  	{
//  	  tmp1 = 0.7*c1[2] + 0.3*(c1[0]+del2*frac1);
//  	  c1[2] = (tmp1 < c1[3]-d1) ? tmp1 : c1[3]-d1;
//  	}
//        if (secn2)
//  	{
//  	  tmp2 = 0.7*c2[2] + 0.3*(c2[0]+del2*frac2);
//  	  c2[2] = (tmp2 < c2[3]-d2) ? tmp2 : c2[3]-d2;
//  	}
//      }
//    else if (start1)
//      {
//        vector<double>::iterator c1 = func1->coefs_begin();
//        vector<double>::iterator c2 = func2->coefs_end()-2;
//        double del = 0.5*((c1[1]-c1[0])/frac1 + (c2[1]-c2[0])/frac2);
//        if (del > min_length*0.8) {
//  	  del = min_length*0.8;
//        }
//        double del2 = 0.5*((c1[2]-c1[0]) + (c2[1]-c2[-1]));
//        double d1 = 0.1*(c1[3]-c1[2]);
//        double d2 = 0.1*(c2[-1]-c2[-2]);
//        c1[1] = c1[0] + del*frac1;
//        c2[0] = c2[1] - del*frac2;
//        if (secn1)
//  	{
//  	  tmp1 = 0.7*c1[2] + 0.3*(c1[0]+del2*frac1);
//  	  c1[2] = (tmp1 < c1[3]-d1) ? tmp1 : c1[3]-d1;
//  	}
//        if (secn2)
//  	{
//  	  tmp2 = 0.7*c2[-1] + 0.3*(c2[1]-del2*frac2);
//  	  c2[-1] = (tmp2 > c2[-2]+d2) ? tmp2 : c2[-2]+d2;
//  	}
//      }
//    else if (start2)
//      {
//        vector<double>::iterator c1 = func1->coefs_end()-2;
//        vector<double>::iterator c2 = func2->coefs_begin();
//        double del = 0.5*((c1[1]-c1[0]) + (c2[1]-c2[0]));
//        if (del > min_length*0.8) {
//  	  del = min_length*0.8;
//        }
//        double del2 = 0.5*((c1[1]-c1[-1])/frac1 + (c2[2]-c2[0])/frac2);
//        double d1 = 0.1*(c1[-1]-c1[-2]);
//        double d2 = 0.1*(c2[3]-c2[2]);
//        c1[0] = c1[1] - del*frac1;
//        c2[1] = c2[0] + del*frac2;
//        if (secn1)
//  	{
//  	  tmp1 = 0.7*c1[-1] + 0.3*(c1[1]-del2*frac1);
//  	  c1[-1] = (tmp1 > c1[-2]+d1) ? tmp1 : c1[-2]+d1;
//  	}
//        if (secn2)
//  	{
//  	  tmp2 = 0.7*c2[2] + 0.3*(c2[0]+del2*frac2);
//  	  c2[2] = (tmp2 < c2[3]-d2) ? tmp2 : c2[3]-d2;
//  	}
//      }
//    else
//      {
//        vector<double>::iterator c1 = func1->coefs_end()-2;
//        vector<double>::iterator c2 = func2->coefs_end()-2;
//        double del = 0.5*((c1[1]-c1[0]) + (c2[1]-c2[0]));
//        if (del > min_length*0.8) {
//  	  del = min_length*0.8;
//        }
//        double del2 = 0.5*((c1[1]-c1[-1])/frac1 + (c2[1]-c2[-1])/frac2);
//        double d1 = 0.1*(c1[-1]-c1[-2]);
//        double d2 = 0.1*(c2[-1]-c2[-2]);
//        c1[0] = c1[1] - del*frac1;
//        c2[0] = c2[1] - del*frac2;
//        if (secn1)
//  	{
//  	  tmp1 = 0.7*c1[-1] + 0.3*(c1[1]-del2*frac1);
//  	  c1[-1] = (tmp1 > c1[-2]+d1) ? tmp1 : c1[-2]+d1;
//  	}
//        if (secn2)
//  	{
//  	  tmp2 = 0.7*c2[-1] + 0.3*(c2[1]-del2*frac2);
//  	  c2[-1] = (tmp2 > c2[-2]+d2) ? tmp2 : c2[-2]+d2;
//  	}
//      }
//  }

//===========================================================================
void ftChartSurface::averageGridDistrFunc(SplineCurve *func1, bool start1,
					  double l1, bool secn1, int nbints1,
					  SplineCurve *func2, bool start2,
					  double l2, bool secn2, int nbints2)
//===========================================================================
{
    MESSAGE("Not up to date, use other version (taking a vector)!");

    // co[ki] +/- del must stay inside end values! Should also be monotone.
    double length1 = fabs(func1->coefs_begin()[0] - func1->coefs_end()[-1]);
    double length2 = fabs(func2->coefs_begin()[0] - func2->coefs_end()[-1]);
    double lint1 = /*l1/(double)nbints1; */length1/(double)nbints1;
    double lint2 = /*l2/(double)nbints2; */length2/(double)nbints2;
    double frac1 = l1/length1; //length1/l1; //;
    double frac2 = l2/length2; //length2/l2; //;
    double min_length = min(length1, length2);


    // Define interpolation problems
    double aint = func1->endparam() - func1->startparam();
    aint *= frac1; // @@sbr Not stable yet. Not ever I suspect... Not working unless domain
//                    // reflects surface size.
    double ta1 = func1->startparam() + aint/(double)nbints1;
    double ta2 = func1->endparam() - aint/(double)nbints1;
    double bint = func2->endparam() - func2->startparam();
    bint *= frac2; // @@sbr Not stable yet.
    double tb1 = func2->startparam() + bint/(double)nbints2;
    double tb2 = func2->endparam() - bint/(double)nbints2;
    Point pt1, pt2;
    Point p1, p2, p3, p4;
    p1 = func1->ParamCurve::point(func1->startparam());
    p2 = func1->ParamCurve::point(func1->endparam());
    p3 = func2->ParamCurve::point(func2->startparam());
    p4 = func2->ParamCurve::point(func2->endparam());

    pt1 = (start1) ? func1->ParamCurve::point(ta1) : func1->ParamCurve::point(ta2);
    pt2 = (start2) ? func2->ParamCurve::point(tb1) : func2->ParamCurve::point(tb2);
    double t1 = (start1) ? pt1[0] - p1[0] : p2[0] - pt1[0];
    double t2 = (start2) ? pt2[0] - p3[0] : p4[0] - pt2[0];
    double gamma = (lint2*t1 + lint1*t2)/(lint1+lint2);
    if (gamma > 0.8*min_length) {
      gamma = 0.8*min_length;
    }

    vector<double> param1, param2;
    vector<double> data1, data2;
    
    SplineInterpolator interpolate1;
    interpolate1.setFreeConditions();
    SplineInterpolator interpolate2;
    interpolate2.setFreeConditions();

    // First curve

    // Start condition
    param1.push_back(func1->startparam());
    data1.push_back(p1[0]);

    // Second condition
    if (start1 || !secn1)
      {
	param1.push_back(ta1);
	data1.push_back(start1 ? p1[0]+gamma :
			(func1->ParamCurve::point(ta1))[0]);
      }
    else
      interpolate1.setNaturalStartCondition();

    // Third condition
    if (!start1 || !secn1)
      {
	param1.push_back(ta2);
	data1.push_back(!start1 ? p2[0]-gamma :
			(func1->ParamCurve::point(ta2))[0]);
      }
    else
      interpolate1.setNaturalEndCondition();
    
    // Last condition
    param1.push_back(func1->endparam());
    data1.push_back(p2[0]);

    // Second curve
    // Start condition
    param2.push_back(func2->startparam());
    data2.push_back(p3[0]);

    // Second condition
    if (start2 || !secn2)
      {
	param2.push_back(tb1);
	data2.push_back(start2 ? p3[0]+gamma :
			(func2->ParamCurve::point(tb1))[0]);
      }
    else
      interpolate2.setNaturalStartCondition();

    // Third condition
    if (!start2 || !secn2)
      {
	param2.push_back(tb2);
	data2.push_back(!start2 ? p4[0]-gamma :
			(func2->ParamCurve::point(tb2))[0]);
      }
    else
      interpolate2.setNaturalEndCondition();
    
    // Last condition
    param2.push_back(func2->endparam());
    data2.push_back(p4[0]);
	
    // Solve interpolation problems
    vector<double> coefs1, coefs2;
    interpolate1.interpolate((int)param1.size(), 1, &param1[0], &data1[0], coefs1);
    interpolate2.interpolate((int)param2.size(), 1, &param2[0], &data2[0], coefs2);

    // Copy result to distribution functions
    int ki;
    vector<double>::iterator c1 = func1->coefs_begin();
    vector<double>::iterator c2 = func2->coefs_begin();
    for (ki=0; ki<(int)coefs1.size(); ki++)
      c1[ki] = coefs1[ki];
    for (ki=0; ki<(int)coefs2.size(); ki++)
      c2[ki] = coefs2[ki];
}

//===========================================================================
void ftChartSurface::averageGridDistrFunc(vector<SplineCurve*>& funcs, vector<bool> start,
					  vector<double> in_lengths, vector<bool> secn,
					  const vector<int>& nbints)
//===========================================================================
{
    // co[ki] +/- del must stay inside end values! Should also be monotone.
    int nmb = (int)funcs.size();
    vector<double> lengths(nmb), lints(nmb), fracs(nmb);
    int ki;
    for (ki = 0; ki < nmb; ++ki) {
	lengths[ki] = fabs(funcs[ki]->coefs_begin()[0] - funcs[ki]->coefs_end()[-1]);
	lints[ki] = /*in_*/lengths[ki]/(double)nbints[ki]; //lengths[ki]/(double)nbints[ki];
	fracs[ki] = in_lengths[ki]/lengths[ki]; //length1/l1; //;
    }
    vector<double>::const_iterator min_length_iter = min_element(in_lengths.begin(), in_lengths.end());
    double min_length = *min_length_iter;
    int min_ind = (int)(min_length_iter - in_lengths.begin());

    // Define interpolation problems
    vector<double> ints(nmb), t1s(nmb), t2s(nmb), ts(nmb);
    vector<Point> p1s(nmb), p2s(nmb), p3(nmb), p4(nmb);
    for (ki = 0; ki < nmb; ++ki) {
	ints[ki] = funcs[ki]->endparam() - funcs[ki]->startparam();
//  	ints[ki] *= fracs[ki]; // @@sbr Not stable yet.
	t1s[ki] = funcs[ki]->startparam() + ints[ki]/(double)nbints[ki];
	t2s[ki] = funcs[ki]->endparam() - ints[ki]/(double)nbints[ki];
	p1s[ki] = funcs[ki]->ParamCurve::point(funcs[ki]->startparam());
	p2s[ki] = funcs[ki]->ParamCurve::point(funcs[ki]->endparam());
	p3[ki] = funcs[ki]->ParamCurve::point(t1s[ki]);
	p4[ki] = funcs[ki]->ParamCurve::point(t2s[ki]);
	ts[ki] = (start[ki]) ? p3[ki][0] - p1s[ki][0] : p2s[ki][0] - p4[ki][0];
	ts[ki] *= fracs[ki]; // We redefine length as the one in space.
    }

    double gamma = 0.0;
    // New method: taking the sum over all ts*lints for largest/smallest elements.
    vector<double> ts_sorted(ts.begin(), ts.end());
    sort(ts_sorted.begin(), ts_sorted.end());
    vector<double> lints_sorted(lints.begin(), lints.end());
    sort(lints_sorted.begin(), lints_sorted.end());
    reverse(lints_sorted.begin(), lints_sorted.end());
    double sum_prod = 0.0;
    double sum_lints = 0.0;
    for (ki = 0; ki < (int)ts_sorted.size(); ++ki) {
	sum_prod += ts_sorted[ki]*lints_sorted[ki];
	sum_lints += lints_sorted[ki];
    }
    gamma = sum_prod/sum_lints;
    // We allow first cell to be 2.0 times as long as the even distribution.
    double length_frac = 2.0/((double)nbints[min_ind]);
    if (gamma > length_frac*min_length) {
	MESSAGE("gamma too large: " << gamma << ", new gamma: " << length_frac*min_length);
	gamma = length_frac*min_length;
    }

    for (ki = 0; ki < nmb; ++ki) {
	vector<double> params, data;
	SplineInterpolator interpolate;
	interpolate.setFreeConditions();

	// Start condition
	params.push_back(funcs[ki]->startparam());
	data.push_back(p1s[ki][0]);

	// Second condition
	if (start[ki] || !secn[ki]) {
	    params.push_back(t1s[ki]);
	    data.push_back(start[ki] ? p1s[ki][0]+gamma/fracs[ki] :
			   (funcs[ki]->ParamCurve::point(t1s[ki]))[0]);
	} else {
	    if (false) { // @@sbr Should be replaced by curvature test!
		// Since we got here there must be an interpolation cond in end of cv.
		// Hence we adjust accordingly in the start of cv.
		params.push_back(t1s[ki]);
		// We compute the new distance from start pt.
// 	    double new_s = (p3[ki][0] - p1s[ki][0])*(p2s[ki][0] - gamma/fracs[ki] - p1s[ki][0]) /
// 		(p4[ki][0] - p1s[ki][0]);
		// This second version seems smarte: comparing lengths towards end pt.
		double new_s = (p3[ki][0] - p1s[ki][0])*(p2s[ki][0] - p4[ki][0])/
		    (gamma/fracs[ki]);
		data.push_back(p1s[ki][0] + new_s);
	    } else {
		interpolate.setNaturalStartCondition();
	    }
	}

	// Third condition
	if (!start[ki] || !secn[ki]) {
	    params.push_back(t2s[ki]);
	    data.push_back(!start[ki] ? p2s[ki][0]-gamma/fracs[ki] :
			    (funcs[ki]->ParamCurve::point(t2s[ki]))[0]);
	} else {
	    if (false) { // @@sbr Should be replaced by curvature test!
		// Since we got here there must be an interpolation cond in start of cv.
		// Hence we adjust accordingly in the end of cv.
		params.push_back(t2s[ki]);
		// We compute the new distance from end pt.
// 	    double new_e = (p2s[ki][0] - p4[ki][0])*(p2s[ki][0] - gamma/fracs[ki] - p1s[ki][0]) /
// 		(p2s[ki][0] - p3[ki][0]);
		double new_e = (p2s[ki][0] - p4[ki][0])*(p3[ki][0] - p1s[ki][0])/
		    (gamma/fracs[ki]);
		data.push_back(p2s[ki][0] - new_e);
	    } else {
		interpolate.setNaturalEndCondition();
	    }
	}

	// Last condition
	params.push_back(funcs[ki]->endparam());
	data.push_back(p2s[ki][0]);

	// Solve interpolation problems
	vector<double> new_coefs;
	interpolate.interpolate((int)params.size(), 1, &params[0], &data[0], new_coefs);

#ifdef FANTASTIC_DEBUG
	// We perform a debug test to see if the coefs are strictly increasing.
	int kj;
	for (kj = 1; kj < (int)new_coefs.size() - 1; ++kj) // Currently only two steps.
	    if ((new_coefs[kj] < new_coefs[0]) || (new_coefs[kj] > new_coefs[new_coefs.size()-1])) {
		MESSAGE("Grid distr coefs not strictly increasing!!! Expect failure.");
	    }
#endif // FANTASTIC_DEBUG

	// Copy result to distribution function
	copy(new_coefs.begin(), new_coefs.end(), funcs[ki]->coefs_begin());
    }
}

//===========================================================================
SplineCurve* ftChartSurface::getGridDistrFunc(int bdidx, double par1,
					      double par2, double par_tol)
//===========================================================================
{
  double ta = grid_distr_functions_[bdidx]->startparam();
  double tb = grid_distr_functions_[bdidx]->endparam();
  if (!((fabs(ta-par1)<par_tol || fabs(ta-par2)<par_tol) &&
	(fabs(tb-par1)<par_tol || fabs(tb-par2)<par_tol))) {
    return 0;
  } else {
    return grid_distr_functions_[bdidx].get();
  }
}

//===========================================================================
SplineCurve* ftChartSurface::getPriorDistrFunc(int bdidx, bool opposite,
					       bool& start, double& length,
					       bool& secn, int& nbints)
//===========================================================================
{
  int previdx;

  if (opposite)
    {
      switch(bdidx)
	{
	case 0:
	  previdx = 3;
	  start = true;
	  break;
	case 1:
	  previdx = 3;
	  start = false;
	  break;
	case 2:
	  previdx = 1;
	  start = true;
	  break;
	case 3:
	  previdx = 1;
	  start = false;
	  break;
	}
    }
  else
    switch(bdidx)
      {
      case 0:
	previdx = 2;
	start = true;
	break;
      case 1:
	previdx = 2;
	start = false;
	break;
      case 2:
	previdx = 0;
	start = true;
	break;
      case 3:
	previdx = 0;
	start = false;
	break;
      }

  bool dir_u = (previdx < 2) ? true : false; // : true;
  double par;
  switch(previdx)
    {
    case 0:
      par = surf_->startparam_v();
      break;
    case 1:
      par = surf_->endparam_v();
      break;
    case 2:
      par = surf_->startparam_u();
      break;
    case 3:
      par = surf_->endparam_u();
      break;
    } 
  GeometryTools::estimateIsoCurveLength(*surf_, dir_u, par, length);
  secn = secn_distr_[previdx];
  secn_distr_[previdx] = false;
  nbints = (previdx <= 1) ? m_ - 1 : n_ - 1;
  return grid_distr_functions_[previdx].get(); 
}

//===========================================================================
void ftChartSurface::sampleGridPts()
//===========================================================================
{
   // The sampling is performed as the tensor product
   // p(u, v) = ((vmax - v)*f1(u)/(vmax - vmin) + v*f2(u)/(vmax - vmin),
   //            (umax - u)*g1(v)/(umax - umin) + u*g2(v)/(umax - umin)), where
   // fi: [umin, umax] -> [umin, umax] && gi: [vmin, vmax] - > [vmin, vmax]
   // are the grid_distr_functions_.

    // We sample points on surface according to the distribution functions and grid size.
    double umin = surf_->startparam_u();
    double umax = surf_->endparam_u();
    double vmin = surf_->startparam_v();
    double vmax = surf_->endparam_v();
    double step_u = (umax - umin)/(m_ - 1);
    double step_v = (vmax - vmin)/(n_ - 1);
    grid_pts_.clear();
    int ki, kj, kk, km, kn;
    for (ki = 0; ki < m_; ++ki) {
	vector<shared_ptr<ftSurfaceSetPoint> > iso_grid_pts;
	double u = umin + ki*step_u;
	Point f1u = grid_distr_functions_[0]->ParamCurve::point(u);
	Point f2u = grid_distr_functions_[1]->ParamCurve::point(u);
	double u2 = 0.5*(f1u[0]+f2u[0]);
	for (kj = 0; kj < n_; ++kj) {
	    double v = vmin + kj*step_v;
	    Point g1v = grid_distr_functions_[2]->ParamCurve::point(v);
	    Point g2v = grid_distr_functions_[3]->ParamCurve::point(v);
//    	    double new_u = (vmax - v)*f1u[0]/(vmax - vmin) + 
//  	      v*f2u[0]/(vmax - vmin);
//    	    double new_v = (umax - u)*g1v[0]/(umax - umin) + 
//  	      u*g2v[0]/(umax - umin);
	    double v2 = 0.5*(g1v[0]+g2v[0]);
	    double new_u = (vmax - v2)*f1u[0]/(vmax - vmin) + 
	      v2*f2u[0]/(vmax - vmin);
	    double new_v = (umax - u2)*g1v[0]/(umax - umin) + 
	      u2*g2v[0]/(umax - umin);
	    u2 = new_u;
	    v2 = new_v;
	    new_u = (vmax - v2)*f1u[0]/(vmax - vmin) + 
	      v2*f2u[0]/(vmax - vmin);
	    new_v = (umax - u2)*g1v[0]/(umax - umin) + 
	      u2*g2v[0]/(umax - umin);
	    if (new_u < umin || new_u > umax || new_v < vmin || new_v > vmax) {
		double knot_diff_tol = 1e-10; // Issue warning if this is not of numerical nature.
		if (new_u + knot_diff_tol < umin || new_u - knot_diff_tol > umax ||
		    new_v + knot_diff_tol < vmin || new_v - knot_diff_tol > vmax) {
		    MESSAGE("Something wrong with distr functions. "
			       "(new_u, new_v) outside domain, moving it inside.");
		}
		new_u = max(umin, min(umax, new_u));
		new_v = max(vmin, min(vmax, new_v));
		if (new_u - umin < knot_diff_tol)
		  new_u = umin;
		if (umax - new_u < knot_diff_tol)
		  new_u = umax;
		if (new_v - vmin < knot_diff_tol)
		  new_v = vmin;
		if (vmax - new_v < knot_diff_tol)
		  new_v = vmax;
	    }
	    Vector2D par_pt(new_u, new_v);
	    Point appr_pt = surf_->ParamSurface::point(new_u, new_v);
	    shared_ptr<ftFaceBase> face;
	    double local_u = new_u;
	    double local_v = new_v;
	    Point sf_pt;
	    double dist = 1e10; // Initializing dist (as it is used when comparing vs appr_pt).
	    try {
		sf_pt = point(local_u, local_v, face);
		dist = appr_pt.dist(sf_pt);
	    } catch (...) { // This should be stabilized in a better manner, but not now...
		if ((ki == 0) && (kj == 0)) {
		    local_u = umin;
		    local_v = vmin;
		    sf_pt = surf_->ParamSurface::point(local_u, local_v);
		    dist = appr_pt.dist(sf_pt);
		} else if ((ki == m_ - 1) && (kj == n_ - 1)) {
		    local_u = umax;
		    local_v = vmax;
		    sf_pt = surf_->ParamSurface::point(local_u, local_v);
		    dist = appr_pt.dist(sf_pt);
		} else {
		    // @@sbr sf_pt not initialized, but hardly within tolerance.
		    MESSAGE("Not stable yet... Fix!"); //, UnknownError());
		}
	    }

#ifdef FANTASTIC_DEBUG
	    std::ofstream debug("tmp/debug.g2");
	    vector<double> pts(6);
	    copy(appr_pt.begin(), appr_pt.end(), pts.begin());
	    copy(sf_pt.begin(), sf_pt.end(), pts.begin() + 3);
	    LineCloud lc(pts.begin(), 1);
	    lc.writeStandardHeader(debug);
	    lc.write(debug);
#endif // FANTASTIC_DEBUG

	    if (dist > approxtol_) {
// 	       MESSAGE("Failed finding point on chart_sf within approxtol_, trying new method.");
		// We try to use par_pt from previous evaluation as seed.
		// The dir index closest to the boundary edge is to be held locked.
		// If we fail dist test for chosen face, we try with neighbour face.
		vector<shared_ptr<ftFaceBase> > faces;
		faces.push_back(face);
		shared_ptr<ftSurface> ft_sf = dynamic_pointer_cast<ftSurface, ftFaceBase>(face);
		if (ft_sf.get() != 0) {
		    ftEdgeBase* clo_edge = ft_sf->edgeClosestToPoint(local_u, local_v);
		    ftEdgeBase* twin_edge = clo_edge->twin();
		    if (twin_edge != 0) {
			ftFaceBase* twin_face = twin_edge->face();
			for (kk = 0; kk < (int)faces_.size(); ++kk) {
			    if (faces_[kk].get() == twin_face) {
				faces.push_back(faces_[kk]);
				break;
			    }
			}
		    }

		    for (kk = 0; kk < (int)faces.size(); ++kk) {
			face = faces[kk];

			// If necessary we try 4 neighbour pts (if they exist) to find matching face.
			vector<shared_ptr<ftSurfaceSetPoint> > neighbour_sf_set_pts;
			if (kj > 0) {
			    neighbour_sf_set_pts.push_back(iso_grid_pts[kj-1]);
			}
			if (ki > 0) {
			    neighbour_sf_set_pts.push_back(grid_pts_[ki-1][kj]);
			    if (kj > 0) {
				neighbour_sf_set_pts.push_back(grid_pts_[ki-1][kj-1]);
			    }
			    if (kj < n_ - 1) {
				neighbour_sf_set_pts.push_back(grid_pts_[ki-1][kj+1]);
			    }
			}
			vector<shared_ptr<ftSurfaceSetPoint> > sf_set_pts;
			for (km = 0; km < (int)neighbour_sf_set_pts.size(); ++km) {
			    for (kn = 0; kn < neighbour_sf_set_pts[km]->nmbFaces(); ++kn) {
				if (neighbour_sf_set_pts[km]->face(kn).get() == face.get()) {
				    sf_set_pts.push_back(neighbour_sf_set_pts[km]);
 				    break;
				}
			    }
// 			    if (kn < neighbour_sf_set_pts[km]->nmbFaces()) {
// 				break;
// 			    }
			}
// 			if (sf_set_pt.get() != 0) {
			for (km = 0; km < (int)sf_set_pts.size(); ++km) {
			    Vector2D prev_par_pt;
			    try {
				prev_par_pt = sf_set_pts[km]->getPar(face.get());
				double new_local_u = new_u;
				double new_local_v = new_v;
				bool use_input_face = true; //(kk == 1);
				Point new_sf_pt = point(new_local_u, new_local_v, face,
							prev_par_pt.begin(), use_input_face);
				double new_dist = new_sf_pt.dist(appr_pt);
				if (new_dist < dist) {
				    dist = new_dist;
// 				    if (new_dist > approxtol_) {
//  // 			       MESSAGE("Still not there yet.. Using new pt anyhow.");
// 				    }
				    sf_pt = new_sf_pt;
				    local_u = new_local_u;
				    local_v = new_local_v;
				} else {
// 			   MESSAGE("Failed finding closesr pt, using previous attempt!");
				}
			    } catch (...) {
				MESSAGE("Not stable yet...  Picked the wrong face?");
			    }
			}
			if (dist < approxtol_) {
			    break;
			}
		    }
		}
	    }

	    ASSERT(sf_pt.size() != 0);
	    Vector3D space_pt(sf_pt.begin());
	    int bd = ((ki == 0) || (ki == m_ - 1) || (kj == 0) || (kj == n_ - 1)) ? 1 : 0;
	    Vector2D local_par(local_u, local_v);
	    shared_ptr<ftSurfaceSetPoint> sample_pt(new ftSurfaceSetPoint(space_pt, bd,
									  face, local_par));
	    sample_pt->setPar(par_pt);

	    if (dist > 30*approxtol_) {
		// If we failed creating appr within tol many samples will be outside approxtol_.
		MESSAGE("Failed finding point on chart_sf within 30*approxtol_, "
			   "using found pt anyway.");

#ifdef FANTASTIC_DEBUG
		std::ofstream debug("tmp/debug.g2");
		vector<double> pts(6);
		copy(appr_pt.begin(), appr_pt.end(), pts.begin());
		copy(space_pt.begin(), space_pt.end(), pts.begin() + 3);
		LineCloud lc(pts.begin(), 1);
		lc.writeStandardHeader(debug);
		lc.write(debug);
#endif // FANTASTIC_DEBUG

	    }

	    iso_grid_pts.push_back(sample_pt);
	}
	grid_pts_.push_back(iso_grid_pts);
    }
}

//===========================================================================
ftMessage ftChartSurface::getMatchingEdges(vector<ftEdgeBase*> local_outer_loop, vector<int> corners,
					   RotationInfo* rot_info, ftCurve* total_outer_bd,
					   vector<pair<int, int> >& matching_edges,
					   vector<ftChartSurface*>& matching_faces,
					   vector<double>& rot_angles)
//===========================================================================
{
    ftMessage status;
    int ki, kj;
   if (rot_info != NULL) {
       ASSERT(total_outer_bd != 0); // Either both are NULL or neither.
       vector<pair<shared_ptr<ParamCurve>, ftFaceBase*> > total_outer_loop_cvs =
	   cmUtils::getG1FaceCurves(*total_outer_bd);

       vector<shared_ptr<ParamCurve> > local_outer_loop_cvs;
       for (ki = 0; ki < (int)corners.size() - 1; ++ki) {
	   ftEdgeBase* curr_edge = local_outer_loop[corners[ki]];
	   ftEdgeBase* twin_edge = curr_edge->twin();
	   if (twin_edge != 0) {
	       continue;   // Not an outer edge
	   }
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	   shared_ptr<ParamCurve> first_space_cv
	       ((ParamCurve*)(curr_edge->geomEdge()->geomCurve()->geometryCurve()->clone()));
#else
	   shared_ptr<ParamCurve> first_space_cv
	       (curr_edge->geomEdge()->geomCurve()->geometryCurve()->clone());
#endif
	   for (kj = corners[ki] + 1; kj < corners[ki+1]; ++kj) {
	       if (local_outer_loop[kj]->twin() != 0) {
		   break;
	       }
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	       shared_ptr<ParamCurve> next_space_cv
		   ((ParamCurve*)(local_outer_loop[kj]->geomEdge()->geomCurve()->geometryCurve()->clone()));
#else
	       shared_ptr<ParamCurve> next_space_cv
		   (local_outer_loop[kj]->geomEdge()->geomCurve()->geometryCurve()->clone());
#endif
	       first_space_cv->appendCurve(next_space_cv.get());
	   }
	   if (kj == corners[ki+1]) {
// 	       shared_ptr<ParamCurve> space_cv;
// 	       if (edge_space_cv->instanceType() == Class_SplineCurve) {
// 		   space_cv = edge_space_cv;
// 	       } else if (edge_space_cv->instanceType() == Class_CurveOnSurface) {
// 		   space_cv = (dynamic_pointer_cast<CurveOnSurface, ParamCurve>(edge_space_cv))->
// 		       spaceCurve();
// 	       } else {
// 		   GO_ERROR("Unexpected curve type!", InputError());
// 	       }
	       local_outer_loop_cvs.push_back(first_space_cv);
	   }
       }

       for (ki = 0; ki < (int)local_outer_loop_cvs.size(); ++ki) {
	   // We rotate edge cv in both positivie and negative direction.
	   shared_ptr<SplineCurve> curr_cv_pos(dynamic_cast<SplineCurve*>(local_outer_loop_cvs[ki]->clone()));
	   ASSERT(curr_cv_pos.get() != 0);
	   curr_cv_pos->reverseParameterDirection();
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	   shared_ptr<SplineCurve> curr_cv_neg((SplineCurve*)(curr_cv_pos->clone()));
#else
	   shared_ptr<SplineCurve> curr_cv_neg(curr_cv_pos->clone());
#endif
	   GeometryTools::translateSplineCurve(-(rot_info->center_pt_), *curr_cv_pos);
	   GeometryTools::rotateSplineCurve(rot_info->rot_axis_, rot_info->rot_angle_, *curr_cv_pos);
	   GeometryTools::translateSplineCurve(rot_info->center_pt_, *curr_cv_pos);
	   GeometryTools::translateSplineCurve(-(rot_info->center_pt_), *curr_cv_neg);
	   GeometryTools::rotateSplineCurve(rot_info->rot_axis_, -(rot_info->rot_angle_), *curr_cv_neg);
	   GeometryTools::translateSplineCurve(rot_info->center_pt_, *curr_cv_neg);

	   for (kj = 0; kj < (int)total_outer_loop_cvs.size(); ++kj) {
	       ftChartSurface* other_chart =
		   dynamic_cast<ftChartSurface*>(total_outer_loop_cvs[kj].second);
	       if (other_chart == 0) {
		   continue;
	       }

	       double epsgeo = toptol_.gap;
	       bool match_pos = GeometryTools::isCoincident(*curr_cv_pos,
					     *total_outer_loop_cvs[kj].first, epsgeo);
	       bool match_neg = GeometryTools::isCoincident(*curr_cv_neg,
					     *total_outer_loop_cvs[kj].first, epsgeo);
	       //   This routine could also return the endpoints of the matching curves
	       //   in order to be able to call getMatchingEdges.
	       if (match_pos || match_neg) {
		   // Need more info to find edge1, edge2 and neighbour_face.
		   // Need also a flag to tell that the grid points of
		   // neighbour_face must be rotated similar to cv2.
		   int edge1 = -1;
		   int edge2 = -1;
		   double par1, par2, par3, par4;
		   double epsgeo = toptol_.gap;
		   shared_ptr<SplineSurface> sf1 =
		       dynamic_pointer_cast<SplineSurface, ParamSurface>(this->surface());
		   shared_ptr<SplineSurface> sf2 = dynamic_pointer_cast<SplineSurface, ParamSurface>
		       (total_outer_loop_cvs[kj].second->surface());
		   if ((sf1.get() == 0) || (sf2.get() == 0)) {
		       status.setError(FT_NOT_SPLINE_SURF);
		       return status;
		   }
		   Point from_space_pt1 = local_outer_loop_cvs[ki]->ParamCurve::point
		       (local_outer_loop_cvs[ki]->startparam());
		   Point to_space_pt1 = local_outer_loop_cvs[ki]->ParamCurve::point
		       (local_outer_loop_cvs[ki]->endparam());
		   sf1->getBoundaryIdx(from_space_pt1, to_space_pt1,
				       epsgeo, edge1, par1, par2);
		   Point from_space_pt2 = total_outer_loop_cvs[kj].first->ParamCurve::point
		       (total_outer_loop_cvs[kj].first->startparam());
		   Point to_space_pt2 = total_outer_loop_cvs[kj].first->ParamCurve::point
		       (total_outer_loop_cvs[kj].first->endparam());
		   sf2->getBoundaryIdx(from_space_pt2, to_space_pt2,
				       epsgeo, edge2, par3, par4);
		   if ((edge1 == -1) || (edge2 == -1)) {
		       status.setError(FT_ERROR_IN_SURFACE_GRIDDING);
		       return status;
		   }
		   matching_edges.push_back(make_pair(edge1, edge2));
		   matching_faces.push_back(other_chart);
		   double alpha = match_pos ? -rot_info->rot_angle_ : rot_info->rot_angle_;
		   rot_angles.push_back(alpha);
	       }
	   }
       }
   }

   return status;
}

//===========================================================================
ftMessage ftChartSurface::getMatchingEdges(vector<ftEdgeBase*> outer_loop, vector<int> corners,
					   vector<pair<int, int> >& matching_edges,
					   vector<pair<int, int> >& matching_grid_res,
					   vector<ftChartSurface*>& matching_faces,
					   vector<double>& rot_angles)
//===========================================================================
{
  ftMessage status;

   int ki;
   for (ki = 0; ki < (int)corners.size() - 1; ++ki) {
      ftEdgeBase* curr_edge = outer_loop[corners[ki]];
      ftEdgeBase* twin_edge = curr_edge->twin();
      ftChartSurface* neighbour_face =
	 (twin_edge != 0) ? dynamic_cast<ftChartSurface*>(twin_edge->face()) : 0;
      if (neighbour_face == 0 || neighbour_face->surface().get() == 0) {
	 continue;
      }
      int other_res_m, other_res_n;
      neighbour_face->gridCreated(other_res_m, other_res_n);

      while ((curr_edge != outer_loop[corners[(ki+1)%4]]) && (twin_edge != 0)) {
	 curr_edge = curr_edge->next();
	 twin_edge = curr_edge->twin();
	 if ((twin_edge != 0) && (twin_edge->face() != neighbour_face)) {
	    twin_edge = 0;
	 }
      }
      if (curr_edge == outer_loop[corners[(ki+1)%4]]) {
	  int edge1 = -1;
	  int edge2 = -1;
	  double par1, par2, par3, par4;
	  // We then find which of the face's edging that match.
	  ftEdgeBase* first_edge = outer_loop[corners[ki]];
	  ftEdgeBase* last_edge = outer_loop[corners[(ki+1)%4]]; // Beginning of next side.
	  Point from_pt = first_edge->point(first_edge->tMin());
	  Point to_pt = last_edge->point(last_edge->tMin());
	  double epsgeo = toptol_.gap;
	  ftMessage local_status =
	      getMatchingEdges(from_pt, to_pt, this, neighbour_face, edge1, 
			       par1, par2, edge2, par3, par4, epsgeo);
	  if (!local_status.isOK()) {
	    status.setError(local_status.getMessage());
	    return status;
	  }
	  if ((edge1 != -1) && (edge2 != -1)) {
	      matching_edges.push_back(make_pair(edge1, edge2));
	      matching_faces.push_back(neighbour_face);
	      rot_angles.push_back(0.0);
	      int local_res = (edge1 < 2) ? m_ : n_;
	      int other_res = (edge2 < 2) ? other_res_m : other_res_n;
	      matching_grid_res.push_back(make_pair(local_res, other_res));
	  } else {
	      // May occur if approxtol_ allows sf to move outside gap_.
	      MESSAGE("Something not right, failed finding matching edge indices! "
			 "Continuing nonetheless.");
	  }
      }
   }

   return status;
}

//===========================================================================
vector<double> ftChartSurface::getEdgeScales()
//===========================================================================
{
    int ki;
    vector<double> edge_scales(4, 1.0);
    if (surf_.get() == 0) {
	MESSAGE("sf_ was not created... Returning default values.");
	return edge_scales;
    }

    // For each element we must decide which edge it corresponds to.
    for (ki = 0; ki < (int)edge_scales_.size(); ++ki) {
	int bd_idx; // b, t, l, r
	double epsgeo = toptol_.gap;
	surf_->getBoundaryIdx(edge_scales_[ki].first, epsgeo, bd_idx);

	edge_scales[bd_idx] = edge_scales_[ki].second;
    }

    return edge_scales;
}

} // namespace Go
