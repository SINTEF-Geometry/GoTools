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

#ifndef _TRIMCRVUTILS_H
#define _TRIMCRVUTILS_H


#include "GoTools/geometry/GeomObject.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Domain.h"
#include "GoTools/utils/Point.h"

#include <vector>


namespace Go
{
  /// Utility functionality related to trimming of LR B-spline surfaces with
  /// respect to a corresponding point cloud. Turn polygons in the parameter domain
  /// of the surface into spline curves appropriate for trimming.
  namespace TrimCrvUtils
  {
    /// File input
    std::vector<double> readTrimPoints(std::ifstream& filein,
				       Point& translate_vec);

    /// Split polygon in corners and ensure that the polygons make a closed loop
    std::vector<std::vector<double> > 
      splitTrimPoints(std::vector<double>& pts_2d,
		      double epsgeo, double kink_tol);

    /// Extract constant parameter sub sequences from polygon
    std::vector<std::vector<double> > 
      extractConstParSeqs(std::vector<double>& pts_2d, 
			  int parix, double par,
			  int nmb_match,
			  double tol1, 
			  double tol2);

    SplineCurve* 
      createVariableOffsetDistFunction(const SplineCurve& par_cv, const std::vector<double>& pts_2d);

    SplineCurve* 
      createOffsetTrimCurve(const SplineCurve& par_cv, const SplineCurve& offset_dist);

    /// Approximate point sequence by a spline curve
    // Assuming the 3D trim points correspond to a 2.5D surface given by (u, v, f(u,v)), we discard the z parameter when
    // approximating the points with a curve in the parameter domain. Assuming the curve should form a simple loop.
    shared_ptr<SplineCurve> 
      approximateTrimPts(std::vector<double> trim_pts_3d,
			 int dim, double epsgeo, int max_iter);

    /// We make sure that all the coefs in the parameter curves lie inside the domain given by the surface.
    // Assuming that the surface domain is rectangular.
    void moveCurveCoefsInsideSurfDomain(const ParamSurface* sf,
					std::vector<shared_ptr<SplineCurve> >& par_cvs);

    /// We approximate the 2D curve with a set of linear segments.
    // Currently we compare the distance between the spline coefs, which is an upper bound.
    std::vector<shared_ptr<Line> > 
      approximateCurve(const SplineCurve& par_curve, double epspar);

    void largestLineCoefDistance(const SplineCurve& curr_segment,
				 int& max_dist_ind, double& max_coef_dist, double& max_coef_par);

    /// Make sure that the trimming curve lies in the parameter domain of the given surface
    shared_ptr<SplineCurve> 
      clipToDomain(const SplineCurve& par_cv, const Domain& domain);

    // We split the trim_pts_2d in kinks. End pts are repeated to ensure continuity between segments.
    // We furthermore expect the input to be closed (i.e. form a loop with equal end points).
    /// Split polygon in corners
    std::vector<std::vector<double> > 
      splitCurvePointsInKinks(const std::vector<double>& trim_pts_2d,
			      double kink_tol);

    /// Ensure that the polygon defines a closed loop
    void makeConnectedLoop(std::vector<std::vector<double> >& trim_seqs_2d, double epsgeo);

    /// Translation of points, compute vector to origin
    void translateToOrigin(GeomObject& go_object, Point& translate_vec);

    /// Do the actual translation given vector of translation
    void translateObject(GeomObject& go_object, const Point& translate_vec);

    /// Translation of surface
    void translateSurfaceDomain(ParamSurface* sf, const Point& translate_vec);

    // We rescale the value of the z dimension. Typically relevant for 2.5D surfaces.
    void scaleZ(ParamSurface& sf, double scale_factor);

    // from and to referr to index of a 2-tuple in pts_2d. We allow negative indices as
    // well as indices larger than the number of 2-tuples, with the modulo operation
    // defining the correct index.
    Point averageDirection(const std::vector<double>& pts_2d, int from, int to);

    /// Tolerance in parameter space
    double getEpsgeo(const std::vector<double>& pts_2d);


  }
}

#endif // _IQMULUSUTILS_H

