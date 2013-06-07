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

#include "GoTools/qualitymodule/ModelQuality.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/PointOnCurve.h"

using std::vector;
using std::make_pair;

namespace Go
{

  //===========================================================================
  ModelQuality::ModelQuality(double gap,   // Gap between adjacent surfaces
			     double kink,  // Kink between adjacent surfaces 
			     double approx)
  //===========================================================================
      : toptol_(tpTolerances(gap, 10.0*gap, kink, 10.0*kink)),
	approx_(approx), small_size_(20.0*gap), curvature_radius_(100.0*gap)
  {
    results_ = shared_ptr<QualityResults>(new QualityResults());
    // results_.reset(new QualityResults());
  }


  //===========================================================================
  ModelQuality::ModelQuality(const tpTolerances& toptol, 
			     double approx)
  //===========================================================================
      : toptol_(tpTolerances(toptol.gap, toptol.neighbour, toptol.kink, toptol.bend)),
	approx_(approx), small_size_(20.0*toptol.gap), curvature_radius_(100.0*toptol.gap)
  {
    results_ = shared_ptr<QualityResults>(new QualityResults());
    // results_.reset(new QualityResults());
  }


  //===========================================================================
  ModelQuality::~ModelQuality()
  //===========================================================================
  {
  }


  //===========================================================================
  void ModelQuality::degenSurfaces(vector<shared_ptr<ParamSurface> >& deg_sfs)
  //===========================================================================
  {
      deg_sfs.clear();

      // No surfaces 
      return;
  }


  //===========================================================================
  void ModelQuality::degenerateSfCorners(vector<shared_ptr<ftPoint> >& deg_corners)
  //===========================================================================
  {
      deg_corners.clear();

      // No surfaces 
      return;
  }

  //===========================================================================
  void ModelQuality::identicalVertices(vector<pair<shared_ptr<Vertex>,
				       shared_ptr<Vertex> > >& identical_vertices)
  {
      identical_vertices.clear();

      // No vertices 
      return;
  }


   //===========================================================================
  void ModelQuality::identicalOrEmbeddedEdges(vector<pair<shared_ptr<ftEdge>,
						shared_ptr<ftEdge> > >& identical_edges,
						vector<pair<shared_ptr<ftEdge>,
						shared_ptr<ftEdge> > >& embedded_edges)

  //===========================================================================
  {
      identical_edges.clear();
      embedded_edges.clear();

      // No edges 
      return;
  }

 //===========================================================================
  void ModelQuality::identicalOrEmbeddedFaces(std::vector<std::pair<shared_ptr<ftSurface>,
					      shared_ptr<ftSurface> > >& identical_faces,
					      std::vector<std::pair<shared_ptr<ftSurface>,
					      shared_ptr<ftSurface> > >& embedded_faces)
  //===========================================================================
  {
      identical_faces.clear();
      embedded_faces.clear();

      // No surfaces 
      return;
  }


  //===========================================================================
  void ModelQuality::miniEdges(vector<shared_ptr<ftEdge> >& mini_edges)
  //===========================================================================
  {
      mini_edges.clear();

      return;
  }

  //===========================================================================
  void ModelQuality::miniSurfaces(vector<shared_ptr<ftSurface> >& mini_surfaces)
  //===========================================================================
  {
      mini_surfaces.clear();

      return;
  }

  //===========================================================================
  void ModelQuality::vanishingSurfaceNormal(vector<shared_ptr<ftPoint> >& singular_points,
					    vector<shared_ptr<ftCurve> >& singular_curves)
  //===========================================================================
  {
      singular_points.clear();
      singular_curves.clear();

      // No surfaces 
      return;
  }


  //===========================================================================
  void 
  ModelQuality::vanishingCurveTangent(vector<shared_ptr<PointOnCurve> >& sing_points,
				      vector<pair<shared_ptr<PointOnCurve>, 
				      shared_ptr<PointOnCurve> > >& sing_curves)
  //===========================================================================
  {
      sing_points.clear();
      sing_curves.clear();

      // No surfaces 
      return;
  }


  //===========================================================================
  void ModelQuality::sliverSurfaces(vector<shared_ptr<ParamSurface> >& sliver_sfs,
				    double thickness,
				    double factor)
  //===========================================================================
  {
      sliver_sfs.clear();

      // No surfaces 
      return;
  }


  //===========================================================================
  void ModelQuality::narrowRegion(vector<pair<shared_ptr<PointOnEdge>, 
				  shared_ptr<PointOnEdge> > >& narrow_regions)
  //===========================================================================
  {
      narrow_regions.clear();

      return;
  }

  //===========================================================================
  void ModelQuality::edgeVertexDistance(vector<pair<ftEdge*,
					shared_ptr<Vertex> > >& edge_vertices)
  //===========================================================================
  {
      edge_vertices.clear();

      // No surfaces 
      return;
  }


  //===========================================================================
  void ModelQuality::faceVertexDistance(vector<pair<ftSurface*,
					shared_ptr<Vertex> > >& face_vertices)
  //===========================================================================
  {
      face_vertices.clear();

      // No surfaces 
      return;
  }

  //===========================================================================
  void ModelQuality::faceEdgeDistance(vector<pair<ftSurface*, ftEdge*> >& face_edges)
  //===========================================================================
  {
    face_edges.clear();

    // No surfaces 
    return;
  }

  //===========================================================================
  void ModelQuality::edgePosAndTangDiscontinuity(vector<pair<ftEdge*, ftEdge*> >& pos_disconts,
						 vector<pair<ftEdge*, ftEdge*> >& tang_disconts)
  //===========================================================================
  {
    pos_disconts.clear();
    tang_disconts.clear();
    return;
  }

  //===========================================================================
  void ModelQuality::facePositionDiscontinuity(vector<pair<ftEdge*, ftEdge*> >& pos_disconts)
  //===========================================================================
  {
    pos_disconts.clear();

    // No surfaces 
    return;
  }

  //===========================================================================
  void ModelQuality::faceTangentDiscontinuity(vector<pair<ftEdge*, ftEdge*> >& tangent_disconts)
  //===========================================================================
  {
    tangent_disconts.clear();

    // No surfaces 
    return;
  }

  //===========================================================================
  void ModelQuality::loopOrientationConsistency(vector<shared_ptr<Loop> >& inconsistent_loops)
  //===========================================================================
  {
    inconsistent_loops.clear();
    // No surfaces 
    return;
  }

  //===========================================================================
  void ModelQuality::faceNormalConsistency(vector<shared_ptr<ftSurface> >& inconsistent_faces)
  //===========================================================================
  {
    inconsistent_faces.clear();
    // No surfaces 
    return;
  }


  //===========================================================================
  void ModelQuality::sfG1Discontinuity(vector<shared_ptr<ftSurface> >& discont_sfs)
  //===========================================================================
  {
    discont_sfs.clear();

    // No surfaces 
    return;
  }

  //===========================================================================
  void ModelQuality::sfC1Discontinuity(vector<shared_ptr<ftSurface> >& discont_sfs)
  //===========================================================================
  {
    discont_sfs.clear();

    // No surfaces 
    return;
  }

  //===========================================================================
  void ModelQuality::cvC1G1Discontinuity(vector<shared_ptr<ParamCurve> >& c1_discont,
					 vector<shared_ptr<ParamCurve> >& g1_discont)
  //===========================================================================
  {
      c1_discont.clear();
      g1_discont.clear();

    return;
  }

  //===========================================================================
  void ModelQuality::cvCurvatureRadius(std::vector<std::pair<shared_ptr<PointOnCurve>, double> >& small_curv_rad,
				       std::pair<shared_ptr<PointOnCurve>, double>& minimum_curv_rad)
  //===========================================================================
  {
      shared_ptr<PointOnCurve> dummy;
      small_curv_rad.clear();
      minimum_curv_rad = make_pair(dummy, MAXDOUBLE);

      return;
  }

  //===========================================================================
  void ModelQuality::sfCurvatureRadius(std::vector<std::pair<shared_ptr<ftPoint>, double> >& small_curv_rad,
				       std::pair<shared_ptr<ftPoint>, double>& minimum_curv_rad)
  //===========================================================================
  {
      shared_ptr<ftPoint> dummy;
      small_curv_rad.clear();
      minimum_curv_rad = make_pair(dummy, MAXDOUBLE);

      return;
  }

  //===========================================================================
    void ModelQuality::acuteEdgeAngle(vector<pair<ftEdge*, ftEdge*> >& edge_acute)
  //===========================================================================
    {
	edge_acute.clear();

	return;
    }

  //===========================================================================
    void ModelQuality::acuteFaceAngle(vector<pair<ftSurface*, ftSurface*> >& face_acute)
  //===========================================================================
    {
	face_acute.clear();

	return;
    }

  //===========================================================================
    void ModelQuality::loopIntersection(vector<pair<shared_ptr<PointOnEdge>, 
					shared_ptr<PointOnEdge> > >& loop_intersection)
  //===========================================================================
    {
	loop_intersection.clear();

	return;
    }

  //===========================================================================
    void ModelQuality::loopSelfIntersection(vector<pair<shared_ptr<PointOnEdge>, 
					shared_ptr<PointOnEdge> > >& loop_self_intersection)
  //===========================================================================
    {
	loop_self_intersection.clear();

	return;
    }

  //===========================================================================
    void ModelQuality::indistinctKnots(vector<shared_ptr<ParamCurve> >& cv_knots,
				       vector<shared_ptr<ParamSurface> >& sf_knots,
				       double tol)
  //===========================================================================
    {
	cv_knots.clear();
	sf_knots.clear();

	return;
    }

} // namespace Go
