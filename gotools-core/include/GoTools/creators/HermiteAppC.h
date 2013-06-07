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

#ifndef _HERMITEAPPC_H_
#define _HERMITEAPPC_H_

#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/creators/HermiteGrid1D.h"


namespace Go
{

/// This class is used to generate a SplineCurve from a EvalCurve
/// using Hermite interpolation.  The generated curve will approximate the
/// EvalCurve within specified tolerances.

class HermiteAppC
{
public:

    /// Constructor where the tolerances and the curve to approximate are specified.
    /// \param crv the curve that we want to generate a Hermite approximation of.
    ///            The curve is \em not copied, only pointed to by the HermiteAppC.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalCurves
    HermiteAppC(EvalCurve* crv, double tolerance1, double tolerance2);

    /// Constructor where the tolerances and the curve to approximate are specified,
    /// as well as the parameters for which we will sample the input curve before
    /// starting the approximating process.
    /// \param crv the curve that we want to generate a Hermite approximation of.
    ///            The curve is \em not copied, only pointed to by the HermiteAppC.
    /// \param initpars pointer to the array of parameter values for which we will
    ///        sample the input curve.
    /// \param n number of parameter values in the array 'initpars'.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalCurves
    HermiteAppC(EvalCurve* crv, double initpars[],int n, 
		  double tolerance1, double tolerance2);

    /// Empty destructor
    ~HermiteAppC(){}

    /// Refine the internal sampling of the curve to approximate such that
    /// the Hermite interpolated curve of this sampling approximates parametrically
    /// the original curve within the specified tolerance.
    void refineApproximation();	// Refine Hermite Grid.

    /// Return the cubic spline curve Hermite interpolating the grid.
    /// \return the cubic spline curve that Hermite interpolates the sampled 
    ///         points of the EvalCurve specified in the constructor.
    shared_ptr<SplineCurve> getCurve();	

 private:
    EvalCurve* curve_;	// Pointer to original curve existing outside *this.
    const double tol1_;
    const double tol2_;
    double min_interval_;	// Smaller intervals are not refined
    HermiteGrid1D grid_;
//     shared_ptr<SplineCurve> curve_approx_; // Spline representation of approximation

    // Distance to evaluator ok (with current grid)? Return value: 1 =
    // ok, 0 = insert new_knot, -1 = failed.
    int testSegment(int j, double& new_knot);
    int bisectSegment(int);
    bool method_failed_;


};

} // namespace Go;

#endif
