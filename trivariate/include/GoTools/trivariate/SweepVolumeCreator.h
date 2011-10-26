//===========================================================================
//
// File : SweepVolumeCreator.h
//
// Created: Thu Dec  4 07:55:18 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: SweepVolumeCreator.h,v 1.3 2008-12-11 09:10:23 kfp Exp $
//
// Description: Class with static methods for volume creation by sweep methods
//
//===========================================================================




#ifndef __SWEEPVOLUMECREATOR_H
#define __SWEEPVOLUMECREATOR_H


#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"

namespace Go
{

  /// Class with static methods for volume creation by sweep methods
  class SweepVolumeCreator
  {

  public:

    // Destructor
    virtual ~SweepVolumeCreator() { }

    /// Create a linearly swept volume as a tensor product of a surface and a curve. If \f$P\f$ is the point to
    /// be swept, the surface has  B-spline coefficients \f$S_{i,j}\f$ and the curve has B-spline coefficients
    /// \f$C_k\f$ then the volume has B-spline coefficients \f$S_{i,j}+C_k-P\f$. If the point lies on the curve,
    /// the curve will be swept along the surface. If the point lies on the surface, the surface will
    /// be swept along the curve.
    /// \param surface the surface
    /// \param curve the curve
    /// \param pt the point to be swept
    /// \return a pointer to the swept volume
    static SplineVolume* linearSweptVolume(const SplineSurface &surface,
					   const SplineCurve &curve,
					   const Point &pt);

    /// Create a volume as the rotation of a surface around an axis
    /// \param surface the surface
    /// \param angle the rotation angle (in radians). If more than \f$2\pi\f$, it is truncated down to \f$2\pi\f$.
    ///              If less than \f$-2\pi\f$, it is truncated up to \f$-2\pi\f$.
    /// \param pt a point on the rotation axis
    /// \param axis a vector describing the direction of the rotation axis
    /// \return a pointer to the rotated volume
    static SplineVolume* rotationalSweptVolume(const SplineSurface &surface,
					       double angle,
					       const Point &pt,
					       const Point &axis);

  };    // Class SweepVolumeCreator


} // namespace Go

#endif    // #ifndef __SWEEPVOLUMECREATOR_H

