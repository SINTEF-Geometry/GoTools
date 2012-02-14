//===========================================================================
//
// File : LoftSurfaceCreator.h
//
//===========================================================================

#ifndef _LOFTSURFACECREATOR_H
#define _LOFTSURFACECREATOR_H


#include "GoTools/geometry/CurveLoop.h"


namespace Go {


class SplineSurface;
class SplineCurve;


/// This namespace contains functions used to create lofted surfaces

namespace LoftSurfaceCreator {

  /// Create a lofting surface based on the input curves. The curves are
  /// not changed during the lofting process. The curves must all lie in the
  /// same space
  /// \param first_curve iterator to first iso-curve in the lofted surface.
  /// \param nmb_crvs the number of curves referred to by first_curve.
  /// \return pointer to the created lofting surface.
  SplineSurface* loftSurface(std::vector<shared_ptr<SplineCurve> >::iterator
			   first_curve, int nmb_crvs);


  /// Create a lofting surface interpolating the input curves in the input
  /// parameters. The curves are not changed during the lofting process.
  /// The curves must all lie in the same space
  /// \param first_curve iterator to first iso-curve in the lofted surface.
  /// \param first_param iso parameter to corresponding curve referred to by first_curve.
  /// \param nmb_crvs the number of curves referred to by first_curve.
  /// \return pointer to the created lofting surface.
  SplineSurface* loftSurface(std::vector<shared_ptr<SplineCurve> >::iterator
			   first_curve,
			   std::vector<double>::iterator first_param,
			   int nmb_crvs);


  /// Create a lofting surface based on the input curves. The curves are
  /// changed during the lofting process.
  /// The curves must all lie in the same space.
  /// \param loft_curves the input curves.
  /// \param nmb_crvs the number of input curves.
  /// \return pointer to the created lofting surface.
  SplineSurface* loftSurface(std::vector<shared_ptr<SplineCurve> >& loft_curves,
			     int nmb_crvs);


  /// Create a vector of curves holding a copy of the input curves, but where
  /// the spline bases have been changed (reparametrized, knot inserted and order raised)
  /// so that the B-spline spaces for the new curves all are identical, and where all curves
  /// are now rational if at least on of the input curves were rational.
  /// The input curves are not changed during the process.
  /// The curves must all lie in the same space.
  /// \param first_curve iterator to first input curve.
  /// \param nmb_crvs the number of curves referred to by first_curve.
  /// \return vector holding the unified curves.
  std::vector<shared_ptr<SplineCurve> >
    unifiedCurvesCopy(std::vector<shared_ptr<SplineCurve> >::iterator first_curve,
		      int nmb_crvs);


  /// Calculate iso parameters for the input curves. The curves are expected to be
  /// ordered, i.e. corresponding to increasing iso parameters.
  /// Curves are given iso-parameters in the range 0.0 to param_length.
  /// All the paramter differences between two neighbouring curves will be propotional
  /// to the quadratic mean of the distance between corresponding control points, i.e.
  /// if, for two curves S=first_curve[i] and T=first_curve[i+1], Q is the quadratic
  /// mean of the distance between corresponding control points of S and T, and P is the
  /// difference between the calculated parameter values for S and T, then the ratio
  /// Q/P is the same for all i.
  /// \param first_curve iterator to first iso curve.
  /// \param nmb_crvs the number of input curves.
  /// \param param_length length of parameter domain.
  /// \param params the computed iso parameters for the input curves.
  void makeLoftParams(std::vector<shared_ptr<SplineCurve> >::const_iterator first_curve,
		      int nmb_crvs, double param_length, std::vector<double>& params);



  /// Create a lofting surface interpolating the input curves in the input
  /// parameters. The input curves are expected to have identical B-spline space.
  /// Either all or none of the curves are expected to be rational.   ???
  /// The returned lofted surface will be rational if and only if the curves are rational.
  /// \param first_curve iterator to first iso-curve in the lofted surface.
  /// \param first_param iso parameter to corresponding curve referred to by first_curve.
  /// \param nmb_crvs the number of curves referred to by first_curve.
  /// \return pointer to the created lofting surface.
  SplineSurface* loftSurfaceFromUnifiedCurves(std::vector<shared_ptr<SplineCurve> >::iterator first_curve,
					      std::vector<double>::iterator first_param,
					      int nmb_crvs);



  /// Create a lofting surface interpolating the input curves in the input
  /// parameters. The input curves are expected to have identical B-spline space.
  /// Either all or none of the curves are expected to be rational. In case the curves are rational, we also loft
  /// the denominator function. In this way, we might risk that the resulting surface has
  /// negative controll points, this must be handled by the calling function.
  /// The returned lofted surface will be non-rational if the curves are non-rational.
  /// \param first_curve iterator to first iso-curve in the lofted surface.
  /// \param first_param iso parameter to corresponding curve referred to by first_curve.
  /// \param nmb_crvs the number of curves referred to by first_curve.
  /// \return pointer to the created lofting surface.
  SplineSurface* loftNonrationalSurface(std::vector<shared_ptr<SplineCurve> >::iterator first_curve,
				      std::vector<double>::iterator first_param,
				      int nmb_crvs);


  /// Create a lofting surface interpolating the input rational curves in the input
  /// parameters by Hermite interpolation, to ensure positive denominator function.
  /// The input curves are expected to have identical B-spline space.
  /// All of the curves are expected to be rational.
  /// The returned lofted surface will also be rational.
  /// \param first_curve iterator to first iso-curve in the lofted surface.
  /// \param first_param iso parameter to corresponding curve referred to by first_curve.
  /// \param nmb_crvs the number of curves referred to by first_curve.
  /// \return pointer to the created lofting surface.
  SplineSurface* loftRationalSurface(std::vector<shared_ptr<SplineCurve> >::iterator first_curve,
				   std::vector<double>::iterator first_param,
				   int nmb_crvs);
  

} // namespace LoftSurfaceCreator


} // namespace Go


#endif // _LOFTSURFACECREATOR_H
