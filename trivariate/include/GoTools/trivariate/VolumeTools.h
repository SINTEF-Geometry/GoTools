//===========================================================================
//
// File : VolumeTools.h
//
// Created: Tue Nov 25 10:58:13 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: VolumeTools.h,v 1.2 2008-11-27 12:59:19 kfp Exp $
//
// Description:
//
//===========================================================================


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
  std::shared_ptr<SplineCurve>
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
  std::shared_ptr<SplineVolume>
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
  std::shared_ptr<SplineSurface>
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
  std::shared_ptr<SplineVolume>
    representSurfaceAsVolume(const SplineSurface& surface,
			     int sf_dir1,
			     int sf_dir2,
			     const BsplineBasis& other_bas,
			     bool rational);

  bool getVolAdjacencyInfo(std::shared_ptr<ParamVolume> vol1,
			   std::shared_ptr<SurfaceOnVolume> vol_sf1,
			   std::shared_ptr<ParamVolume> vol2,
			   std::shared_ptr<SurfaceOnVolume> vol_sf2,
			   double tol,
			   int& bd1, int& bd2, int& orientation,
			   bool& same_seq);

  bool getCorrCoefVolEnum(std::shared_ptr<SplineVolume> vol1,
			  std::shared_ptr<SplineVolume> vol2,
			  int bd1, int bd2, int orientation,
			  bool same_seq, 
			  std::vector<std::pair<int, int> >& enumeration);

  bool getVolCoefEnumeration(std::shared_ptr<SplineVolume> vol, int bd,
			     std::vector<int>& enumeration);

  bool getVolBdCoefEnumeration(std::shared_ptr<SplineVolume> vol, int bd,
			       int bd_cv, std::vector<int>& enumeration);

 std::vector<std::shared_ptr<ParamSurface> > 
    getBoundarySurfaces(std::shared_ptr<ParamVolume> vol);

 std::vector<std::shared_ptr<ParamSurface> > 
    getOrientedBoundarySurfaces(std::shared_ptr<ParamVolume> vol);

 std::shared_ptr<SurfaceOnVolume> 
   getBoundarySurface(std::shared_ptr<SplineVolume> vol, int idx);

 std::shared_ptr<SurfaceOnVolume> 
   getOrientedBoundarySurface(std::shared_ptr<SplineVolume> vol, int idx);

 void volCommonSplineSpace(std::shared_ptr<SplineVolume> vol1, int bd1,
			   std::shared_ptr<SplineVolume> vol2, int bd2,
			   int orientation, bool same_seq);

 // Approximate a parameter curve in the parameter space of a given volume
 // by a space curve
 std::shared_ptr<SplineCurve> 
   liftVolParamCurve(std::shared_ptr<ParamCurve> pcurve, 
		     std::shared_ptr<ParamVolume> vol,
		     double tol);

 // Approximate a space curve by a curve in the parameter space of a 
 // given volume
  
 std::shared_ptr<SplineCurve> 
   projectVolParamCurve(std::shared_ptr<ParamCurve> spacecurve, 
		     std::shared_ptr<ParamVolume> vol,
		     double tol);
} // namespace Go





#endif // _VOLUMETOOLS_H

