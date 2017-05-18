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

#ifndef __CMUTILS_H
#define __CMUTILS_H


#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/compositemodel/ftCurve.h"



namespace Go
{

    
struct RotationInfo
{
    Go::Point center_pt_;
    Go::Point rot_axis_;
    double rot_angle_; ///< Measured in radians.
};

class ParamSurface;

  /// Various utility functions for the compositemodel module
namespace cmUtils
{

  // Return 2-dimensional point representing parameter values of input point.
  // Point faceParameter(ftEdgeBase* edge, double t);

  /// Estimate the curve length corresponding to an edge
  double estimatedCurveLength(ftEdgeBase* edge, int nmb_samples = 4);

  /// Domain of surface is rescaled according to surface lengths in geometry space.
  RectDomain geometricParamDomain(ParamSurface* sf);

  /// Return true if pt lies above plane as defined
  bool abovePlane(Go::Point pt, Go::Point plane_pt, Go::Point normal);

    /// Extend the set of boundary curves with degenerate curves if
    /// a degenerate surface is to be created. Handles 2 or 3 bnd_curves.
    void extendWithDegBd(std::vector<int>& corner, 
			 std::vector< shared_ptr<Go::ParamCurve> >& bd_curves, 
			 std::vector< shared_ptr<Go::ParamCurve> >& cross_curves, 
			 int idxmin);

    /// Assuming input corner pts lie inside gap of outer_loop, and not in the middle of a segment.
    void updateWithNewCorners(const std::vector<Go::ftEdgeBase*>& outer_loop, std::vector<int>& corners,
			      const std::vector<Go::Point>& add_corner_pts, double gap);


        /// We decide whether sample_pt lies inside polygonial domain defined by bd_train.
    /// bd_train must be a simple loop, direction may be either.
    bool insideBdTrain(const Go::Vector2D& sample_pt, const std::vector<Go::Vector2D>& bd_train);

    // Based on curvature the bd_cvs are reparametrized.
    void reparametrizeBdCvs(const std::vector<shared_ptr<Go::SplineCurve> >& bd_cvs,
			    double appr_tol,
			    std::vector<shared_ptr<Go::SplineCurve> >& new_bd_cvs);

    // Almost the same, slightly altered parametrization (avg curv of neighbour pts).
    void reparametrizeBdCvs2(const std::vector<shared_ptr<Go::SplineCurve> >& bd_cvs,
			     double appr_tol,
			     std::vector<shared_ptr<Go::SplineCurve> >& new_bd_cvs);

        /// We make sure that all corner_pts (within gap from an edge) do not lie in the
    /// interior of an edge. To be called prior to connection with twin.
    void splitEdgesInCorners(std::vector<shared_ptr<Go::ftEdgeBase> >& edges,
			     const std::vector<Go::Point>& corner_pts, double gap);

    
    /// Given input of edges, meeting in a common vertex (as given by start), make sure that the edges
    /// after the first element is sorted based on angle, in cw direction.
    /// In addition we also make sure that edges with common sf_id are given in sequence.
    /// Note: Expecting that all sfs have consistent orientation (which they'll have if
    /// topology has been built)!!!
    void cwOrientation(std::vector<Go::ftEdgeBase*>& meeting_edges, std::vector<bool>& start,
		       double angle_tol = 1e-05);

    /// Method does not compute tangents and angles but is based on next_, prev_ & twin_.
    /// We're assuming that at most two input edges are without twin.
    void cwOrientation2(std::vector<Go::ftEdgeBase*>& meeting_edges, std::vector<bool>& start);

        /// Divide the input curve into G1 segments lying on one surface only (input expected to lie on bd).
  std::vector<std::pair<shared_ptr<Go::ParamCurve>, Go::ftFaceBase*> >
    getG1FaceCurves(Go::ftCurve& bd_curve);

        /// Return angle from first_leg to second_leg, in ccw direction. Dimension is 2 or 3. For dimension
    /// 3 the input normal is used to define plane. Returned angle lies in [0, 2*pi).
    /// normal assumed to exist if dim of legs equals 3, need not be normal to legs. Only defines ccw.
    double ccwAngle(Go::Point first_leg, Go::Point second_leg, Go::Point* normal = NULL);

    std::vector<int> removeInnerCorners(const std::vector<Go::ftEdgeBase*>& outer_loop, std::vector<int>& corners);
    
}    // Namespace cmUtils

} // namespace Go

#endif    // #ifndef __CMUTILS_H
