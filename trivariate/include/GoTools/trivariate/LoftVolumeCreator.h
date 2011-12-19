//===========================================================================
//
// File : LoftVolumeCreator.h
//
// Created: Wed Dec 17 12:39:36 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: LoftVolumeCreator.h,v 1.1 2008-12-19 09:54:29 kfp Exp $
//
// Description: Methods for volume creation
//
//===========================================================================

#ifndef _LOFTVOLUMECREATOR_H
#define _LOFTVOLUMECREATOR_H


#include "GoTools/geometry/CurveLoop.h"


namespace Go {


class SplineVolume;
class SplineSurface;
class SplineCurve;


/// This namespace contains functions used to create lofted volumes

namespace LoftVolumeCreator {

  /// Create a lofting volume based on the input surfaces. The surfaces are
  /// not changed during the lofting process. The surfaces must all lie in the
  /// same space
  /// \param first_surface iterator to first iso-surface in the lofted volume.
  /// \param nmb_srfs the number of surfaces referred to by first_surface.
  /// \return pointer to the created lofting volume.
  SplineVolume* loftVolume(std::vector<shared_ptr<SplineSurface> >::iterator
			   first_surface, int nmb_srfs);


  /// Create a lofting volume interpolating the input surfaces in the input
  /// parameters. The surfaces are not changed during the lofting process.
  /// The surfaces must all lie in the same space
  /// \param first_surface iterator to first iso-surface in the lofted volume.
  /// \param first_param iso parameter to corresponding surface referred to by first_surface.
  /// \param nmb_srfs the number of surfaces referred to by first_surface.
  /// \return pointer to the created lofting volume.
  SplineVolume* loftVolume(std::vector<shared_ptr<SplineSurface> >::iterator
			   first_surface,
			   std::vector<double>::iterator first_param,
			   int nmb_srfs);


  /// Create a vector of surfaces holding a copy of the input surfaces, but where
  /// the spline bases have been changed (reparametrized, knot inserted and order raised)
  /// so that the B-spline spaces in the u-direction for the new surfaces all are identical,
  /// same with v-direction, and where all surfaces are now rational if at least on
  /// of the input surfaces were rational. The input surfaces are not changed during the process.
  /// The surfaces must all lie in the same space.
  /// \param first_surface iterator to first input surface.
  /// \param nmb_srfs the number of surfaces referred to by first_surface.
  /// \return vector holding the unified surfaces.
  std::vector<shared_ptr<SplineSurface> >
    unifiedSurfacesCopy(std::vector<shared_ptr<SplineSurface> >::iterator first_surface,
			int nmb_srfs);


  /// Calculate iso parameters for the input surfaces. The surfaces are expected to be
  /// ordered, i.e. corresponding to increasing iso parameters.
  /// Surfaces are given iso-parameters in the range 0.0 to param_length.
  /// All the paramter differences between two neighbouring surfaces will be propotional
  /// to the quadratic mean of the distance between corresponding control points, i.e.
  /// if, for two surfaces S=first_surface[i] and T=first_surface[i+1], Q is the quadratic
  /// mean of the distance between corresponding control points of S and T, and P is the
  /// difference between the calculated parameter values for S and T, then the ratio
  /// Q/P is the same for all i.
  /// \param first_surface iterator to first iso surface.
  /// \param nmb_crvs the number of input surfaces.
  /// \param param_length length of parameter domain.
  /// \param params the computed iso parameters for the input surfaces.
  void makeLoftParams(std::vector<shared_ptr<SplineSurface> >::const_iterator first_surface,
		      int nmb_crvs, double param_length, std::vector<double>& params);



  /// Create a lofting volume interpolating the input surfaces in the input
  /// parameters. For each parameter direction (u or v), the input surfaces are expected
  /// to have identical B-spline space. Either all or none of the surfaces are expected to be rational.
  /// The returned lofted volume will be rational if and only if the surfaces are rational.
  /// \param first_surface iterator to first iso-surface in the lofted volume.
  /// \param first_param iso parameter to corresponding surface referred to by first_surface.
  /// \param nmb_srfs the number of surfaces referred to by first_surface.
  /// \return pointer to the created lofting volume.
  SplineVolume* loftVolumeFromUnifiedSurfaces(std::vector<shared_ptr<SplineSurface> >::iterator first_surface,
					      std::vector<double>::iterator first_param,
					      int nmb_srfs);



  /// Create a lofting volume interpolating the input surfaces in the input
  /// parameters. For each parameter direction (u or v), the input surfaces are expected
  /// to have identical B-spline space. Either all or none of the surfaces are expected
  /// to be rational. In case the surfaces are rational, we also loft
  /// the denominator function. In this way, we might risk that the resulting volume has
  /// negative controll points, this must be handled by the calling function.
  /// The returned lofted volume will non-rational if the surfaces are non-rational.
  /// \param first_surface iterator to first iso-surface in the lofted volume.
  /// \param first_param iso parameter to corresponding surface referred to by first_surface.
  /// \param nmb_srfs the number of surfaces referred to by first_surface.
  /// \return pointer to the created lofting volume.
  SplineVolume* loftNonrationalVolume(std::vector<shared_ptr<SplineSurface> >::iterator first_surface,
				      std::vector<double>::iterator first_param,
				      int nmb_srfs);


  /// Create a lofting volume interpolating the input rational surfaces in the input
  /// parameters by Hermite interpolation, to ensure positive denominator function.
  /// For each parameter direction (u or v), the input surfaces are expected
  /// to have identical B-spline space. All of the surfaces are expected to be rational.
  /// The returned lofted volume will also be rational.
  /// \param first_surface iterator to first iso-surface in the lofted volume.
  /// \param first_param iso parameter to corresponding surface referred to by first_surface.
  /// \param nmb_srfs the number of surfaces referred to by first_surface.
  /// \return pointer to the created lofting volume.
  SplineVolume* loftRationalVolume(std::vector<shared_ptr<SplineSurface> >::iterator first_surface,
				   std::vector<double>::iterator first_param,
				   int nmb_srfs);
  

} // namespace LoftVolumeCreator


} // namespace Go


#endif // _LOFTVOLUMECREATOR_H
