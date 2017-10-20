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

#ifndef LRSURFSTITCH_H
#define LRSURFSTITCH_H

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Direction2D.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"

namespace Go
{
  /// This class stitches a collection of LR B-spline surfaces to
  /// ensure C0 or C1 continuity between adjacent surfaces. The surfaces
  /// must be organized in a regular grid. Not all surfaces in the grid
  /// must exist.
  /// The parameter domain of the surfaces must join such that adjacent
  /// surfaces have adjacent, none-overlapping parameter domains without
  /// any gaps.
  /// The surfaces must be oriented consistently
  /// C1 continuity is implemented only for 1D surfaces
  /// Due to the restrictions on the surface configurations only the
  /// organizing function is made public. Otherwise, general functionality
  /// is made private as they do not check if the configuration 
  /// constraints are satisfied.

  class LRSurfStitch
  {
  public:
    /// Perform stitching. The surfaces are organized from bottom to
    /// top and from left to right
    /// Only LR B-spline surfaces are modified, but the surfaces may
    /// be trimmed
    void stitchRegSfs(std::vector<shared_ptr<ParamSurface> >& sfs,
		      int nmb_u, int nmb_v, double eps,
		      int cont);

    /// Perform stitching. The surfaces are organized from bottom to
    /// top and from left to right
    void stitchRegSfs(std::vector<shared_ptr<LRSplineSurface> >& sfs,
		      int nmb_u, int nmb_v, double eps,
		      int cont);

    /// We calculate the max distance (cont = 0) / angle (cont = 1) between corresponding edges.
    std::vector<double> analyzeContinuity(std::vector<shared_ptr<ParamSurface> >& sfs,
					  int nmb_u, int nmb_v, int cont,
					  int num_edge_samples = 100);

  private:
    void consistentSplineSpaces(std::vector<shared_ptr<LRSplineSurface> >& sfs,
				int nmb_u, int nmb_v, double eps,
				int cont);

    // Make sure that the refinements along all edges have a full tensor product structure for
    // the first 'element_width' elements.
    // For corners additional knots are inserted.
    void tensorStructure(shared_ptr<LRSplineSurface> surf, int element_width,
			 bool edges[4]);

    bool matchSplineSpace(shared_ptr<LRSplineSurface> surf1, int edge1,
			  shared_ptr<LRSplineSurface> surf2, int edge2, 
			  int element_width, double tol);

    void checkCornerMatch(std::vector<std::pair<shared_ptr<LRSplineSurface>,int> >& sfs,
			 double tol, std::vector<int>& nmb_match,
			 std::vector<Point>& corner_val, 
			 std::vector<Point>& param);

    int averageCorner(std::vector<std::pair<shared_ptr<ParamSurface>,int> >& sfs,
		      double tol);

   int averageCorner(std::vector<std::pair<shared_ptr<LRSplineSurface>,int> >& sfs,
		      double tol);

   int makeCornerC1(std::vector<std::pair<shared_ptr<LRSplineSurface>,int> >& sfs,
		    double tol);

   bool averageEdge(shared_ptr<ParamSurface> surf1, int edge1,
		     shared_ptr<ParamSurface> surf2, int edge2, double tol);

    bool averageEdge(shared_ptr<LRSplineSurface> surf1, int edge1,
		     shared_ptr<LRSplineSurface> surf2, int edge2, 
		     int cont, double tol);

    void makeLineC1(LRBSpline2D* bsp[4], Direction2D dir);

    void fetchEdgeCorners(shared_ptr<LRSplineSurface> surf, int edge,
			  double& u1, double& v1, double& u2, double& v2);

    void extractMissingKnots(std::vector<double>& union_vec, 
			     std::vector<double>& vec,
			     double tol, int order,
			     std::vector<double>& resvec);

    void defineRefinements(const Mesh2D& mesh, Direction2D dir,
			   int edge, int ix, const std::vector<double>& knot_vals, 
			   int element_width,
			   std::vector<LRSplineSurface::Refinement2D>& refs); 

    void extractBoundaryBsplines(shared_ptr<LRSplineSurface> surf,
				 int edge,
				 std::vector<std::vector<LRBSpline2D*> >& bsplines);

    // Utility function for debugging.
    // Compute distance, tangent and normal differences in input parameter points.
    void surfaceDifference(const LRSplineSurface* sf1, const Point& param1,
			   const LRSplineSurface* sf2, const Point& param2,
//			   bool along_u_dir,
			   double& sf_dist, double& tang1_ang, double& tang2_ang);

};
};

#endif
