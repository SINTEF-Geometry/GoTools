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

#ifndef _INTEGRATE_H_
#define _INTEGRATE_H_
#include <memory>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/config.h"

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
				double lim2, shared_ptr<SplineCurve> bspline_curve,
				double*** integral);


};
#endif
