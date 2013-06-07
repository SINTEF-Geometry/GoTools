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

#ifndef _HERMITEAPPS_H_
#define _HERMITEAPPS_H_

#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/creators/HermiteGrid1DMulti.h"

namespace Go
{

/// This class is used to generate a set of SplineCurves from a EvalCurveSet (which
/// itself represents a set of related curves) using Hermite interpolation.  The generated
/// curves will approximate those defined by the EvalCurveSet within specified 
/// tolerances.  This class is really a generalization of HermiteAppC.

class HermiteAppS
{
public:
    /// Constructor where the tolerances and the curves to approximate are 
    /// specified.
    /// \param surface the curve set that we want to approximate by Hermite 
    ///                interpolation of sampled values.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalCurveSets.
    /// \param dims vector that specifies the spatial dimensions of the curves 
    ///             contained in the EvalCurveSet.  The size of the vector should
    ///             be equal to the total number of curves in 'surf', ie. the return
    ///             value of its EvalCurveSet::nmbCvs() function.
    HermiteAppS(EvalCurveSet* surface, 
		  double tolerance1, 
		  double tolerance2, 
		  std::vector<int> dims);

    /// Constructor where the tolerances and the curve set to approximate are 
    /// specified, as well as the parameters where they should be sampled prior to
    /// the Hermite interpolation.
    /// \param surface the curve set that we want to approximate by Hermite
    ///                interpolation of sampled values.
    /// \param initpars pointer to the array of parameter values for which we will
    ///                 sample the input curves.
    /// \param n number of parameter values in the array 'initpars[]'.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalCurveSets.
    /// \param dims vector that specifies the spatial dimensions of the curves 
    ///             contained in the EvalCurveSet.  The size of the vector should
    ///             be equal to the total number of curves in 'surf', ie. the return
    ///             value of its EvalCurveSet::nmbCvs() function.
    HermiteAppS(EvalCurveSet* surface, double initpars[], int n,
		  double tolerance1, double tolerance2, std::vector<int> dims);

    /// Empty destructor
    ~HermiteAppS(){}

    /// Refine the internal sampling of the curves such that the Hermite 
    /// interpolated curves parametrically approximates the original curve within
    /// a specified tolerance.
    void refineApproximation();	// Refine Hermite Grid.

    /// Return the cubic spline curves intepolating the grid (ie. approximating
    /// the original curve set). 
    /// \return a vector containing shared pointers to the newly created 
    ///         spline curves that Hermite interpolate the sampled points
    ///         of the EvalCurveSet (curve set) specified in the constructor.
    std::vector<shared_ptr<SplineCurve> > getCurves();

private:
    EvalCurveSet* surface_;     // Pointer to original surface existing outside *this.
    HermiteGrid1DMulti grid_; // We define a grid for all curves, separating point values.
    const double tol1_;
    const double tol2_;
    const double min_interval_;	// Smaller intervals are not refined
    std::vector<shared_ptr<SplineCurve> > curve_approx_; // Spline representations of approximation.
    /*   shared_ptr<SplineSurface> surface_approx_; // Spline representation of approximation */

    bool testSegment(int j, double& new_knot);	// Distance to _original
    int bisectSegment(int);


};


} // namespace Go

#endif // _HERMITEAPPS_H_
