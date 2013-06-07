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

#ifndef _MODELQUALITY_H
#define _MODELQUALITY_H

#include "GoTools/topology/tpTolerances.h"
#include "GoTools/qualitymodule/QualityResults.h"
#include "GoTools/qualitymodule/testSuite.h"
#include <vector>

namespace Go
{
    class SurfaceModel;
    class PointOnCurve;

    class ModelQuality
	{
	public:
	    // Constructor
	    ModelQuality(double gap, double kink, double approx);

	    ModelQuality(const tpTolerances& toptol, double approx);

	    // Destructor
	    virtual
	    ~ModelQuality();

	    void setMiniElementSize(double element_size)
		{
		    small_size_ = element_size;
		}

	    void setCurvatureRadius(double radius)
		{
		    curvature_radius_ = radius;
		}

	    /// Does not generate new geometry
	    virtual
	    void degenSurfaces(std::vector<shared_ptr<ParamSurface> >& deg_sfs);

	    /// Generates new points, pointers to existing surfaces
	    virtual
		void degenerateSfCorners(std::vector<shared_ptr<ftPoint> >& deg_corners);

	    /// Uses existing topological entity
	    virtual
		void identicalVertices(std::vector<std::pair<shared_ptr<Vertex>,
					      shared_ptr<Vertex> > >& identical_vertices);

	    /// Uses existing topological entities. May contain generated curves
	    virtual
		void identicalOrEmbeddedEdges(std::vector<std::pair<shared_ptr<ftEdge>,
					      shared_ptr<ftEdge> > >& identical_edges,
					      std::vector<std::pair<shared_ptr<ftEdge>,
					      shared_ptr<ftEdge> > >& embedded_edges);

	    /// Uses existing topological entities. Points to existing surfaces
	    virtual
		void identicalOrEmbeddedFaces(std::vector<std::pair<shared_ptr<ftSurface>,
					      shared_ptr<ftSurface> > >& identical_faces,
					      std::vector<std::pair<shared_ptr<ftSurface>,
					      shared_ptr<ftSurface> > >& embedded_faces);

	    /// Uses existing topological entities
	    virtual 
		void miniEdges(std::vector<shared_ptr<ftEdge> >& mini_edges);


	    /// Uses existing topological entities. Points to existing surfaces
	    virtual 
		void miniSurfaces(std::vector<shared_ptr<ftSurface> >& mini_surfaces);

	    /// Returns ftPoints and ftCurves that will contain generated geometry
	    virtual
		void vanishingSurfaceNormal(std::vector<shared_ptr<ftPoint> >& singular_points,
					    std::vector<shared_ptr<ftCurve> >& singular_curves);

	    /// Contains generated points and curves fetched from topological entities
	    virtual
		void vanishingCurveTangent(std::vector<shared_ptr<PointOnCurve> >& sing_points,
					   std::vector<std::pair<shared_ptr<PointOnCurve>, 
					   shared_ptr<PointOnCurve> > >& sing_curves);

	    /// Uses existing geometry
	    virtual
	      void sliverSurfaces(std::vector<shared_ptr<ParamSurface> >& sliver_sfs,
				  double thickness,
				  double factor = 2.0);

	    /// Uses existing topological entities. Contains generated geometry
	    virtual
		void narrowRegion(std::vector<std::pair<shared_ptr<PointOnEdge>, 
				  shared_ptr<PointOnEdge> > >& narrow_regions);

	    /// Uses existing topological entities
	    virtual
	      void edgeVertexDistance(std::vector<std::pair<ftEdge*,
				    shared_ptr<Vertex> > >& edge_vertices);

	    /// Existing topology. Partly generated geometry, existing surfaces
	    virtual
	      void faceVertexDistance(std::vector<std::pair<ftSurface*,
				    shared_ptr<Vertex> > >& face_vertices);

	    /// Existing topology. Partly generated geometry, existing surfaces
	    virtual
	      void faceEdgeDistance(std::vector<std::pair<ftSurface*, ftEdge*> >& face_edges);

	    /// Uses existing topological entities. 
	    virtual
	      void edgePosAndTangDiscontinuity(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_disconts,
					       std::vector<std::pair<ftEdge*, ftEdge*> >& tang_disconts);

	    /// Uses existing topological entities. Existing surfaces available
	    virtual
	      void facePositionDiscontinuity(std::vector<std::pair<ftEdge*, ftEdge*> >& pos_disconts);

	    /// Uses existing topological entities. Existing surfaces available
	    virtual
	      void faceTangentDiscontinuity(std::vector<std::pair<ftEdge*, ftEdge*> >& tangent_disconts);

	    /// Uses existing topological entities. 
	    virtual
		void loopOrientationConsistency(std::vector<shared_ptr<Loop> >& inconsistent_loops);

	    /// Uses existing topological entities. Existing surfaces available
	    virtual
		void faceNormalConsistency(std::vector<shared_ptr<ftSurface> >& inconsistent_faces);

	    /// Uses existing topological entities. Existing surfaces available
	    virtual
		void sfG1Discontinuity(std::vector<shared_ptr<ftSurface> >& discont_sfs);

	    /// Uses existing topological entities. Existing surfaces available
	    virtual
		void sfC1Discontinuity(std::vector<shared_ptr<ftSurface> >& discont_sfs);

	    /// Contains curves fetched from topological entities
	    virtual
		void cvC1G1Discontinuity(std::vector<shared_ptr<ParamCurve> >& c1_discont,
					 std::vector<shared_ptr<ParamCurve> >& g1_discont);

	    /// Contains generated points and curves fetched from topological entities
	    virtual
		void cvCurvatureRadius(std::vector<std::pair<shared_ptr<PointOnCurve>, double> >& small_curv_rad,
				       std::pair<shared_ptr<PointOnCurve>, double>& minimum_curv_rad);

	    /// Generates new points, pointers to existing surfaces
	    virtual
		void sfCurvatureRadius(std::vector<std::pair<shared_ptr<ftPoint>, double> >& small_curv_rad,
				       std::pair<shared_ptr<ftPoint>, double>& minimum_curv_rad);
	    /// Uses existing topological entities. 
	    virtual
		void acuteEdgeAngle(std::vector<std::pair<ftEdge*, ftEdge*> >& edge_acute);

	    /// Uses existing topological entities. Existing surfaces available
	    virtual
		void acuteFaceAngle(std::vector<std::pair<ftSurface*, ftSurface*> >& face_acute);

	    /// Uses existing topological entities. Contains generated geometry
	    virtual
		void loopIntersection(std::vector<std::pair<shared_ptr<PointOnEdge>, 
				      shared_ptr<PointOnEdge> > >& loop_intersection);

	    /// Uses existing topological entities. Contains generated geometry
	    virtual
		void loopSelfIntersection(std::vector<std::pair<shared_ptr<PointOnEdge>, 
					  shared_ptr<PointOnEdge> > >& loop_self_intersection);

	    /// Uses existing surface and possibly generated curves
	    virtual
		void indistinctKnots(std::vector<shared_ptr<ParamCurve> >& cv_knots,
				     std::vector<shared_ptr<ParamSurface> >& sf_knots,
				     double tol = 1.0e-8);

	    // Get results of quality check
	    shared_ptr<QualityResults> getResults()
	      {
		return results_;
	      }

	protected:
	    tpTolerances toptol_;
	    double approx_;
	    double small_size_;
	    double curvature_radius_;
	    shared_ptr<QualityResults> results_;
	};

} // namespace Go

#endif // _MODELQUALITY_H

