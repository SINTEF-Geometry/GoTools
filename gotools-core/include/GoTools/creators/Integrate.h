#ifndef _INTEGRATE_H_
#define _INTEGRATE_H_
#include <memory>
#include "GoTools/geometry/SplineCurve.h"

namespace Go
{

  class BsplineBasis;

  /// Functions used to compute integrals of inner products of B-splines

  /** Store parameters and weights used for numerical integration
   *  by Gauss quadrature. All parameter values for all Bezier segments
   *  are stored in parameters, while the weights within an interval is
   *  stored only once in weights. Thus, the lenght of weights gives
   *  the number of samples for each Bezier segment.
   * @param basis B-spline basis
   * @param parameters Parameter values for all points
   * @param par_weights Weight for each parameter
   */

  void GaussQuadValues(const BsplineBasis& basis,
		       std::vector<double>& parameters,
		       std::vector<double>& par_weights);


  /** Compute all definite integrals of inner products of
   *  derivatives of B-splines up to a given order where the
   *  differentiation is of the same order for both B-splines.
   *  The interval of integration are equal to the parameter
   *  intervals of the surface in the current par. dir.
   * @param basis B-spline basis.
   * @param ider Number of derivatives to compute.
   * @param lim1 Start of parameter interval.
   * @param lim2 End of parameter interval.
   * @param integral Computed integrals.
   */
    void GaussQuadInner(const BsplineBasis& basis, int ider, double lim1,
			double lim2, double*** integral);


  /** Compute all definite integrals of inner products of
   *  derivatives of B-splines up to a given order where the
   *  gap between the order of differentiation on the first and
   *  second B-spline is constant.
   *  The interval of integration are equal to the parameter
   *  intervals of the surface in the current par. dir.
   * @param basis B-spline basis.
   * @param derivs Number of derivatives to compute.
   * @param gap Difference between derivation order.
   * @param start_der First derivative to compute.
   * @param lim1 Start of parameter interval.
   * @param lim2 End of parameter interval.
   * @param integral Computed integrals.
   */
    void GaussQuadInnerFlat(const BsplineBasis& basis, int derivs, int start_der, int gap,
			    double lim1, double lim2, std::vector<double>& integral);

  /** Compute all definite integrals of inner products of
   *  derivatives of B-splines up to a given order where the
   *  differentiation is of the same order for both B-splines.
   *  The interval of integration are equal to the parameter
   *  intervals of the surface in the current par. dir.
   * @param basis B-spline basis.
   * @param ider Number of derivatives to compute.
   * @param lim1 Start of parameter interval.
   * @param lim2 End of parameter interval.
   * @param integral Computed integrals.
   */
    void GaussQuadInner2(const BsplineBasis& basis, int ider, double lim1,
			 double lim2, double** integral);

  /** Compute all definite integrals of inner products of
   *  derivatives of rational B-splines up to a given order where the
   *  differentiation is of the same order for both B-splines.
   *  The interval of integration are equal to the parameter
   *  intervals of the surface in the current par. dir.
   * @param basis B-spline basis.
   * @param ider Number of derivatives to compute.
   * @param lim1 Start of parameter interval.
   * @param lim2 End of parameter interval.
   * @param bspline_curve 1-dim rational 0-function defining weights
   * @param coefs Wieghts for denominator function
   * @param integral Computed integrals.
   */
    void GaussQuadInnerRational(const BsplineBasis& basis, int ider, double lim1,
				double lim2, std::shared_ptr<SplineCurve> bspline_curve,
				double*** integral);


};
#endif
