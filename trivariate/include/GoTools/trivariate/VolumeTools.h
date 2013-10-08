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

#ifndef _VOLUMETOOLS_H
#define _VOLUMETOOLS_H

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/utils/Array.h"
#include <memory>
#include <vector>



namespace Go
{

  class SurfaceOnVolume;

/// This namespace contains free functions operating on parametric volumes

namespace VolumeTools {

    /// Analyze periodicity of volume based on number of repeating
    /// knots and control points. The return value is -1 if the volume
    /// edges are disjoint, otherwise k if sf is C^k continuous across the
    /// seam. These are sufficient but not necessary conditions for periodicity,
    /// so it is possible that a call to analyzePeriodicityDerivs() will yield a
    /// higher degree of periodicity.
    /// The current implementation is quite slow, and not optimized for speed.
    /// \param sf reference to the SplineVolume to be analyzed
    /// \param direction specify 'direction' to be '0' to check for periodicity in
    ///                  the first parameter direction, or '1' to check
    ///                  the second parameter direction, or 2 for the third.
    /// \param knot_tol the tolerance used when comparing knot intervals
    /// \return -1 if the volume edges are disjoint, otherwise k if the 
    ///         
    int GO_API
    analyzePeriodicity(const SplineVolume& sf, int direction,
                       double knot_tol = 1e-12);

  /// Describe a volume as a high-dimensional curve in a given direction.
  /// If the volume is rational, the curve will be non-rational
  /// and living in the homogenous space.
  /// \param volume the volume to express as a curve
  /// \param cv_dir the parameter direction that will be kept when defining 
  ///               the curve (the other two will disappear, as the control
  ///               points in this direction will be lumped together and expressed
  ///               as single control points in a higher-dimensional space.
  ///               'cv_dir' takes the values 0, 1 or 2 for keeping the first,
  ///               second or third parameter direction respectively
  /// \return shared pointer to a new SplineCurve, expressing the volume
  ///         as a curve in a high-dimensional space.
  shared_ptr<SplineCurve>
    representVolumeAsCurve(const SplineVolume& volume,
			   int cv_dir);

  /// Describe a curve as a lower-dimensional volume in a given direction.
  /// \param curve the curve that we want to express as a volume
  /// \param cv_dir If this variable is set to 0, then the curve's parameter
  ///               will become the \em first parameter in the generated volume.  
  ///               If it is set to 1, the curve's parameter will become the
  ///               \em second parameter in the generated volume.
  ///               If it is set to 2, the curve's parameter will become the
  ///               \em third parameter in the generated volume.  Other values
  ///               are illegal.
  /// \param other_bas1 the first BsplineBasis for the additional parameter directions.
  ///                   If cv_dir is 0, this will be the second parameter direction
  ///                   on the volume. Otherwise, it wil be the first parameter direction
  /// \param other_bas2 the second BsplineBasis for the additional parameter directions.
  ///                   If cv_dir is 2, this will be the second parameter direction
  ///                   on the volume. Otherwise, it wil be the third parameter direction
  /// \param rational define whether the generated volume shall be specified as 
  ///                 \em rational or not.
  /// \return a shared pointer to a new SplineVolume, expressing the curve
  ///         in a space of lower dimensionality.
  shared_ptr<SplineVolume>
    representCurveAsVolume(const SplineCurve& curve,
			   int cv_dir,
			   const BsplineBasis& other_bas1,
			   const BsplineBasis& other_bas2,
			   bool rational);

  /// Describe a volume as a high-dimensional surface in two given directions.
  /// If the volume is rational, the surface will be non-rational
  /// and living in the homogenous space.
  /// \param volume the volume to express as a surface
  /// \param sf_dir1 the volume parameter direction to become the first parameter
  ///                direction of the surface, either 0, 1 or 2.
  /// \param sf_dir2 the volume parameter direction to become the second parameter
  ///                direction of the surface, either 0, 1 or 2, and different from
  ///                sf_dir1. The control points in the direction different from
  ///                sf_dir1 and sf_dir2 will be lumped together and expressed as
  ///                single control points in a higher-dimensional space.
  /// \return shared pointer to a new SplineSurface, expressing the volume
  ///         as a surface in a high-dimensional space.
  shared_ptr<SplineSurface>
    representVolumeAsSurface(const SplineVolume& volume,
			     int sf_dir1,
			     int sf_dir2);

  /// Describe a surface as a lower-dimensional volume in given directions.
  /// \param surface the surface that we want to express as a volume
  /// \param sf_dir1 The volume parameter direction from the first surface parameter
  ///                direction. The value of sf_dir1 must be either 0, 1 or 2. Other values
  ///                are illegal.
  /// \param sf_dir2 The volume parameter direction from the second surface parameter
  ///                direction. The value of sf_dir2 must be either 0, 1 or 2, and different
  ///                from sf_dir1. Other values are illegal.
  /// \param other_bas the BsplineBasis for the additional volume parameter direction, different
  ///                  sf_dir1 and sfdir_2.
  /// \param rational define whether the generated volume shall be specified as 
  ///                 \em rational or not.
  /// \return a shared pointer to a new SplineVolume, expressing the surface
  ///         in a space of lower dimensionality.
  shared_ptr<SplineVolume>
    representSurfaceAsVolume(const SplineSurface& surface,
			     int sf_dir1,
			     int sf_dir2,
			     const BsplineBasis& other_bas,
			     bool rational);


  /// Check if two neighbouring volumes in a volume set have a 
  /// corner-to-corner configuration, i.e. no T-joints
  /// \param vol1 one volume
  /// \param vol_sf1 the boundary face corresponding
  /// to vol1 at the boundary common with vol2
  /// \param vol2 the other volume
  /// \param vol_sf1 the boundary face corresponding
  /// to vol2 at the boundary common with vol1
  /// \param return parameter true of the configuration is corner-to-corner
  bool cornerToCornerVols(shared_ptr<ParamVolume> vol1,
			  shared_ptr<SurfaceOnVolume> vol_sf1,
			  shared_ptr<ParamVolume> vol2,
			  shared_ptr<SurfaceOnVolume> vol_sf2,
			  double tol);

  /// Given two neighbouring volumes, vol1 and vol2, and information
  /// about the common boundary, fetch information about the
  /// adjacency configuration
  /// \param vol1 the first volume
  /// \param vol_sf1 the boundary surface/trimming surface corresponding
  /// to vol1 at the boundary common with vol2
  /// \param vol2 the other volume
  /// \param vol_sf2 the boundary surface/trimming surface corresponding
  /// to vol2 at the boundary common with vol1
  /// \param tol adjacency tolerance
  /// \param return parameter true if the two volumes share the same boundary and
  /// that boundary is boundary trimmed in both volumes
  /// \param bd1 boundary surface of first volume which follows the common 
  /// boundary, -1=no such boundary surface, 0=umin, 1=umax, 2=vmin, 3=vmax,
  /// 4=wmin, 5=wmax
  /// \param bd2 boundary surface of second volume which follows the common 
  /// boundary
  /// \param orientation =0: The surfaces at the common boundary have the 
  /// same orientation in both parameter directions, =1: the orientation is
  /// opposite in the first parameter direction, =2: the orientation is
  /// opposite in the second parameter direction, =3: the orientation is
  /// opposite in both parameter directions
  /// \param same_seq tells if 1. parameter direction of the first boundary
  /// surface corresponds to the 1. parameter direction of the second boundary
  /// surface (true) or if the parameter directions are swapped (false)
  bool getVolAdjacencyInfo(shared_ptr<ParamVolume> vol1,
			   shared_ptr<SurfaceOnVolume> vol_sf1,
			   shared_ptr<ParamVolume> vol2,
			   shared_ptr<SurfaceOnVolume> vol_sf2,
			   double tol,
			   int& bd1, int& bd2, int& orientation,
			   bool& same_seq);

  /// Given information about the adjacency relationship between two spline
  /// volumes (the enumeration of the common boundary and the correspondance
  /// in parameterization, see the function getVolAdjacencyInfo), fetch the
  /// enumeration of pairwise corresponding coefficients. The function returns
  /// false if the spline spaces of the two surfaces at the common boundary
  /// are not the same.
  bool getCorrCoefVolEnum(shared_ptr<SplineVolume> vol1,
			  shared_ptr<SplineVolume> vol2,
			  int bd1, int bd2, int orientation,
			  bool same_seq, 
			  std::vector<std::pair<int, int> >& enumeration);

  /// Given a spline volume and a boundary enumeration (bd=0:umin, 
  /// bd=1:umax, bd=2:vmin, bd=3:vmax, bd=4:wmin, bd=5:wmax), fetch
  /// the enumeration of the volume coefficents at that boundary surface
  bool getVolCoefEnumeration(shared_ptr<SplineVolume> vol, int bd,
			     std::vector<int>& enumeration);

  /// Given a spline volume and a boundary enumeration (bd=0:umin, 
  /// bd=1:umax, bd=2:vmin, bd=3:vmax, bd=4:wmin, bd=5:wmax), fetch
  /// the enumeration of the volume coefficents at that boundary surface
  bool getVolCoefEnumeration(shared_ptr<SplineVolume> vol, int bd,
			     std::vector<int>& enumeration_bd,
			     std::vector<int>& enumeration_bd2);

  /// Given a spline volume, a boundary surface enumeration (bd=0:umin, 
  /// bd=1:umax, bd=2:vmin, bd=3:vmax, bd=4:wmin, bd=5:wmax) and a
  /// boundary curve enumeration corresponding to the boundary surface
  /// (bd_cv=0:umin, bd_cv=1:umax, bd_cv=2:vmin, bd_cv=3:vmax), fetch
  /// the enumeration of the volume coefficents along that boundary curve
  bool getVolBdCoefEnumeration(shared_ptr<SplineVolume> vol, int bd,
			       int bd_cv, std::vector<int>& enumeration);

  /// Given a parametric volume, fetch all boundary surfaces
 std::vector<shared_ptr<ParamSurface> > 
    getBoundarySurfaces(shared_ptr<ParamVolume> vol);

 /// Given a parametric volume, fetch all boundary surfaces and ensure
 /// a constistent normal vector direction of these surfaces
 std::vector<shared_ptr<ParamSurface> > 
    getOrientedBoundarySurfaces(shared_ptr<ParamVolume> vol);

  /// Given a spline volume, fetch all boundary surfaces
 shared_ptr<SurfaceOnVolume> 
   getBoundarySurface(shared_ptr<SplineVolume> vol, int idx);

 /// Given a spline volume, fetch all boundary surfaces and ensure
 /// a constistent normal vector direction of these surfaces
 shared_ptr<SurfaceOnVolume> 
   getOrientedBoundarySurface(shared_ptr<SplineVolume> vol, int idx);

#if 0
  // Currently removed as this is done somewhat differently than the 2D case.
  void
  averageBoundaryCoefs(shared_ptr<SplineVolume>& vol1, int bd1, bool keep_first,
		       shared_ptr<SplineVolume>& vol2, int bd2, bool keep_second,
		       std::vector<bool> found_corners, std::vector<Point> corners,
		       int orientation);
#endif

 /// Given two spline volumes sharing a common boundary and information 
 /// about the configuration of this boundary in the spline volumes
 /// (see the explanation of getVolAdjacencyInfo), ensure that the
 /// two volumes share the same spline space at the common boundary
 void volCommonSplineSpace(shared_ptr<SplineVolume> vol1, int bd1,
			   shared_ptr<SplineVolume> vol2, int bd2,
			   int orientation, bool same_seq);

 /// Approximate a parameter curve in the parameter space of a given volume
 /// by a space curve
 shared_ptr<SplineCurve> 
   liftVolParamCurve(shared_ptr<ParamCurve> pcurve, 
		     shared_ptr<ParamVolume> vol,
		     double tol);

 /// Approximate a given space curve by a curve in the parameter space of a 
 /// given volume
 shared_ptr<SplineCurve> 
   projectVolParamCurve(shared_ptr<ParamCurve> spacecurve, 
		     shared_ptr<ParamVolume> vol,
		     double tol);

 shared_ptr<SplineCurve> 
   approxVolParamCurve(shared_ptr<ParamCurve> spacecurve, 
		       shared_ptr<ParamVolume> vol,
		       double tol, int max_iter, double& maxdist);

} // namespace VolumeTools

} // namespace Go





#endif // _VOLUMETOOLS_H

