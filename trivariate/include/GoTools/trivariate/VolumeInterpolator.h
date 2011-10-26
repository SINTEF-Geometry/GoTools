//===========================================================================
//
// File : VolumeInterpolator
//
// Created: December 2010
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description: Interpolation of regular data points
//
//===========================================================================


#ifndef __VOLUMEINTERPOLATOR_H
#define __VOLUMEINTERPOLATOR_H

#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/utils/Point.h"

namespace Go
{
  class SplineVolume;

/// This namespace contains functions used to interpolate a set of points
  namespace VolumeInterpolator
  {
    /// Interpolate a set of regular, parameterized interpolation
    /// points. The parameterization is assumed to correspond to the
    /// given B-spline basises (select interpolation points in the
    /// Greville points). The function throws if the input parameters
    /// are inconsistent
    SplineVolume* regularInterpolation(const BsplineBasis& basis_u,
				       const BsplineBasis& basis_v,
				       const BsplineBasis& basis_w,
				       std::vector<double>& par_u,
				       std::vector<double>& par_v,
				       std::vector<double>& par_w,
				       std::vector<double>& points,
				       int dimension,
				       bool rational,
				       std::vector<double>& weights);

  };    // namespace VolumeInterpolator


} // namespace Go


#endif    // #ifndef __VOLUMEINTPEROLATOR_H

