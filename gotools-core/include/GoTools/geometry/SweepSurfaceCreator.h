//===========================================================================
//
// File : SweepSurfaceCreator.h
//
// Created: Wed Dec  3 11:07:54 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: SweepSurfaceCreator.h,v 1.3 2008-12-11 09:11:28 kfp Exp $
//
// Description: Class with static methods for surface creation by sweep methods
//
//===========================================================================


#ifndef __SWEEPSURFACECREATOR_H
#define __SWEEPSURFACECREATOR_H

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"

namespace Go
{

  /// Functionality to create a spline surface by linear or rotational sweep
  class SweepSurfaceCreator
  {

  public:

    // Destructor
    virtual ~SweepSurfaceCreator() { }

    /// Create a linearly swept surface as a tensor product of two curves. If \f$P\f$ is the point to
    /// be swept, and the curves have B-spline coefficients \f$C_i\f$ and \f$D_j\f$, then the surface has
    /// B-spline coefficients \f$C_i+D_j-P\f$. If the point lies on the first curve,
    /// the first curve will be swept along the second curve. If the point lies on the second curve,
    /// the second curve will be swept along the first curve.
    /// \param curv1 the first curve
    /// \param curv2 the second curve
    /// \param pt the point to be swept
    /// \return a pointer to the swept surface
    static SplineSurface* linearSweptSurface(const SplineCurve &curv1,
					     const SplineCurve &curv2,
					     const Point &pt);

    /// Create a surface as the rotation of a curve around an axis
    /// \param curve the curve
    /// \param angle the rotation angle (in radians). If less than \f$2\pi\f$, it is truncated down to \f$2\pi\f$.
    ///              If more than \f$-2\pi\f$, it is truncated up to \f$-2\pi\f$.
    /// \param pt a point on the rotation axis
    /// \param axis a vector describing the direction of the rotation axis
    /// \return a pointer to the rotated surface
    static SplineSurface* rotationalSweptSurface(const SplineCurve &curve,
						 double angle,
						 const Point &pt,
						 const Point &axis);


  };    // Class SweepSurfaceCreator


} // namespace Go


#endif    // #ifndef __SWEEPSURFACECREATOR_H
