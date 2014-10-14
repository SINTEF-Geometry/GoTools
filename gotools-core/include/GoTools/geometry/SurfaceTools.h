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

#ifndef _SURFACETOOLS_H
#define _SURFACETOOLS_H

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/RectDomain.h"


namespace Go
{

/// Free functions operating on parametric surfaces
namespace SurfaceTools
{

  /// Fetch the outer boundary loop of a parametric surface
  CurveLoop outerBoundarySfLoop(shared_ptr<ParamSurface> surf,
				double degenerate_epsilon);

  /// Fetch all boundary loops of a parametric surface
  std::vector<CurveLoop> allBoundarySfLoops(shared_ptr<ParamSurface> surf,
					    double degenerate_epsilon); 

  /// Fetch all boundary loops of a parametric surface including degenerate
  /// curves
  std::vector<CurveLoop> 
    absolutelyAllBoundarySfLoops(shared_ptr<ParamSurface> surf,
				 double degenerate_epsilon); 
  /// Compute the corner where a number of surfaces meet. This function is
  /// relevant if the surfaces almost interpolate the same corner. 
  /// If the surface is elementary (cone, sphere, ...), then this surface
  /// has more influence on the the corner position than spline surfaces.
  /// \param vertex the position of the corner. This point is moved during
  /// the computation to get closer to a point shared by all surfaces
  /// \param sfs the surfaces meeting in this corner and the surface parameter
  /// associated with the corner. Note that the parameter does not necessarily
  /// lie at a corner of the underlying non-trimmed surface and it might also
  /// represent a vertex in a surface model and may lie internally at a
  /// surface boundary
  void 
    iterateCornerPos(Point& vertex, 
		     std::vector<std::pair<shared_ptr<ParamSurface>, Point> > sfs,
		     double tol);

  /// Check if two neighbouring surfaces in a surface set have a 
  /// corner-to-corner configuration, i.e. no T-joints
  /// \param sf1 one surface
  /// \param sf_cv1 the boundary curve/trimming curve corresponding
  /// to sf1 at the boundary common with sf2
  /// \param sf2 the other surface
  /// \param sf_cv2 the boundary curve/trimming curve corresponding
  /// to sf2 at the boundary common with sf1
  /// \param return parameter true of the configuration is corner-to-corner
  bool cornerToCornerSfs(shared_ptr<ParamSurface> sf1,
			 shared_ptr<CurveOnSurface> sf_cv1,
			 shared_ptr<ParamSurface> sf2,
			 shared_ptr<CurveOnSurface> sf_cv2,
			 double tol);

  /// Given two neighbouring surfaces, sf1 and sf2, and information
  /// about the common boundary, fetch information about the
  /// adjacency configuration
  /// \param sf1 one surface
  /// \param sf_cv1 the boundary curve/trimming curve corresponding
  /// to sf1 at the boundary common with sf2
  /// \param sf2 the other surface
  /// \param sf_cv2 the boundary curve/trimming curve corresponding
  /// to sf2 at the boundary common with sf1
  /// \param tol adjacency tolerance
  /// \param return parameter true if the two surfaces share the same boundary and
  /// that boundary is boundary trimmed in both surfaces
  /// \param bd1 boundary curve of first surface which follows the common 
  /// boundary, -1=no such boundary curve, 0=umin, 1=umax, 2=vmin, 3=vmax
  /// \param bd2 boundary curve of second surface which follows the common 
  /// boundary
  /// \param same_orient whether or not the surface boundaries have got the
  /// same orientation along the common boundary
  bool getSfAdjacencyInfo(shared_ptr<ParamSurface> sf1,
			  shared_ptr<CurveOnSurface> sf_cv1,
			  shared_ptr<ParamSurface> sf2,
			  shared_ptr<CurveOnSurface> sf_cv2,
			  double tol,
			  int& bd1, int& bd2, bool& same_orient);

  /// Given two spline surfaces and information about a common boundary
  /// between these surfaces, compute the enumeration of corresponding 
  /// coefficients along this boundary provided that the surfaces have
  /// got the same spline space along the boundary
  /// \param sf1 first surface
  /// \param sf2 second surface
  /// \param bd1 boundary curve of first surface which follows the common 
  /// boundary, -1=no such boundary curve, 0=umin, 1=umax, 2=vmin, 3=vmax
  /// \param bd2 boundary curve of second surface which follows the common 
  /// boundary
  /// \param same_orient whether or not the surface boundaries have got the
  /// same orientation along the common boundary
  /// \param return value whether or not the surfaces have got the same
  /// spline space
  /// \param enumeration pairwise enumeration of the corresponding coefficients
  /// in the two surfaces
  bool getCorrCoefEnum(shared_ptr<SplineSurface> sf1,
		       shared_ptr<SplineSurface> sf2,
		       int bd1, int bd2, bool same_orient,
		       std::vector<std::pair<int,int> >& enumeration);

  /// Fetch the local enumeration of the surface coefficients along a
  /// specified boundary (bd=0 - umin, bd=1 - umax, bd=2 - vmin, bd=3 - vmax)
  bool getCoefEnumeration(shared_ptr<SplineSurface> sf, int bd,
			  std::vector<int>& enumeration);

  /// Return both the boundary coefficients as well as coefficient row
  /// number 2 when counting from the edge.
  bool getCoefEnumeration(shared_ptr<SplineSurface> sf, int bd,
			  std::vector<int>& enumeration_bd,
			  std::vector<int>& enumeration_bd2);

  /// Fetch the coefficient enumeration of the corner of a spline surface
  /// specified by the two boundaries meeting in this corner
  bool getCornerCoefEnum(shared_ptr<SplineSurface> sf, int bd1, int bd2,
			  int& enumeration);

  /// Find coefficient enumeration of possible colinear coefficients
  /// at a common edge and check whether the coefficients are indeed
  /// colinear (return value: 0 = no correspondance, 1 = almost colinear, 
  /// 2 = conditions defined)
  int checkCoefCoLinearity(shared_ptr<SplineSurface> sf1,
			    shared_ptr<SplineSurface> sf2,
			    int bd1, int bd2, bool same_orient,
			    double tol, double ang_tol,
			    std::vector<std::vector<int> >& enumeration);

  /// Estimate a start point to a closest point iteration between a given
  /// point, pt, and a given surface, sf, possibly restriced to the parameter 
  /// domain rd.
  void surface_seedfind(const Point& pt, 
			const ParamSurface& sf, 
			const RectDomain* rd,
			double& u,
			double& v);

  /// Parameterize a point set by projecting the points onto a given base 
  /// surface
  void parameterizeByBaseSurf(const  ParamSurface& sf, 
			      const std::vector<double>& points,
			      std::vector<double>& parvals);

  /// Estimate the average length of the cross tangent of a spline
  /// surface at a given boundary
  /// \param surf the surface
  /// \param pardir = 1 : the cross tangent curve has constant u-parameter
  /// pardir = 0 : the cross tangent curve has constant v-parameter
  /// \param at_start the cross tangent curve is evaluated at the start
  /// of the parameter domain in the given parameter direction
  double estimateTangentLength(SplineSurface *surf, int pardir, 
			       bool at_start);

    /// Check if the surface is closed in one or both paramater
    /// directions
    void GO_API checkSurfaceClosed(const Go::ParamSurface& sf,
		        bool& closed_dir_u, bool& closed_dir_v,
		        double closed_tol=1e-06);

    /// Check if the surface is closed in one or both paramater
    /// directions
    void GO_API surfaceClosed(const Go::SplineSurface& sf,
		        bool& closed_dir_u, bool& closed_dir_v,
		        double closed_tol=1e-06);

    // Given a geometric epsilon, we calculate the corresponding global lengths in the parameter domain.
    // The pareps is returned as a 2D point.
    Point getParEpsilon(const ParamSurface& sf, double epsgeo);

} // namespace SurfaceTools
} // namespace Go





#endif // _SURFACETOOLS_H

