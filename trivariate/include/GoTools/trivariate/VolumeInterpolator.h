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
    /// \param basis_u spline basis in the first parameter direction
    /// \param basis_v spline basis in the second parameter direction
    /// \param basis_w spline basis in the third parameter direction
    /// \param par_u parameter values in 1. parameter direction corresponding
    /// to point
    /// \param par_v parameter values in 2. parameter direction corresponding
    /// to point
    /// \param par_w parameter values in 2. parameter direction corresponding
    /// to point
    /// \param points the regular point set
    /// \param dimension dimension of geometry space
    /// \param rational whether or not a rational surface is expected
    /// \param weights the weights of the rational volume, used only if 
    /// rational==true
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

