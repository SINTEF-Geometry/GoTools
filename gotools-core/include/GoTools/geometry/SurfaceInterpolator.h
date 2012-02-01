//===========================================================================
//
// File : SurfaceInterpolator
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


#ifndef __SURFACEINTERPOLATOR_H
#define __SURFACEINTERPOLATOR_H

#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/utils/Point.h"

namespace Go
{
  class SplineSurface;

  /// \brief Functionality for creating interpolating surfaces.
  namespace SurfaceInterpolator
  {
    /// Interpolate a set of regular, parameterized interpolation
    /// points. The parameterizetion is assumed to correspond to the
    /// given B-spline basises (select interpolation points in the
    /// Greville points). The function throws if the input parameters
    /// are inconsistent
    SplineSurface* regularInterpolation(const BsplineBasis& basis_u,
					const BsplineBasis& basis_v,
					std::vector<double>& par_u,
					std::vector<double>& par_v,
					std::vector<double>& points,
					int dimension,
					bool rational,
					std::vector<double>& weights);

  };    // namespace SurfaceInterpolator


} // namespace Go


#endif    // #ifndef __SURFACEINTPEROLATOR_H

