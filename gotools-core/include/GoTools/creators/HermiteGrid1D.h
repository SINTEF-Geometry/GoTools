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

#ifndef _HERMITEGRID1D_
#define _HERMITEGRID1D_

#include <vector>
#include "GoTools/utils/Point.h"

namespace Go
{

class EvalCurve;

/// The type "GoHemiteGrid1D" holds a one dimensional grid containing sampled
/// points and derivatives from a curve.  It can be used to generate bezier
/// curve segments obtained by Hermite interpolation of intervals between
/// the sampled parameter values.

class HermiteGrid1D
{
public:

    /// Construct a HermiteGrid1D from a curve and a given interval.
    /// The grid will only contain the sampled values (position, derivative)
    /// for the first and last value of the interval.
    /// \param crv curve to sample
    /// \param start start of interval
    /// \param end end of interval
    HermiteGrid1D(const EvalCurve& crv, double start, double end);

    /// Construct a HermiteGrid1D from a curve and a set of parameter
    /// values.
    /// \param crv curve to sample
    /// \param param array of strictly increasing parameters contained in
    ///              the parameter domain of 'crv'.
    /// \param n number of elements in 'param[]'.
    HermiteGrid1D(const EvalCurve& crv, double param[], int n);
    
    /// Default destructor.
    ~HermiteGrid1D();

    /// Add another sample (parameter, position, tangent) to the grid.
    /// Returns the index of the new knot (parameter value) in the sorted
    /// knot vector after insertion.  
    /// \param crv curve to evaluate
    /// \param knot the new sample value (parameter value, knot)
    int addKnot(const EvalCurve& crv, double knot);
  
    /// Calculate Bezier coefficients of the cubic curve interpolating 
    /// the point and tangent values at grid nodes with indices "left" 
    /// and "right".
    /// \param left  indicating grid node for start of curve segment
    /// \param right indicating grid node for end of curve segment
    /// \param spar  start parameter of segment
    /// \param epar  end parameter of segment
    /// \param bezcoef array of cubic Bezier coefficients
    void getSegment(int left, int right, double& spar,
		    double& epar, Point bezcoef[4]);
    
    /// Return the grid parameters
    std::vector<double> getKnots() { return knots_; }

    /// Return the sample values (positions and first derivatives)
    std::vector<Point> getData() { return array_; }

    /// Return the spatial dimension
    int dim(){ return dim_; }

    /// Return the number of samples in the grid
    int size() {return MM_;}

private:
  std::vector<double> knots_;     // Sorted array of DISTINCT parameters of curve
  std::vector<Point> array_; 	// Array holding position, and
  				// directional derivative
  int dim_;        		// Spatial dimension of position,
  int MM_;         		// Number of grid points
  int elem_size_;		// Number of Point stored for each
  				// grid point (typically 2: pt, der)
  int index_;                   // Index into knot array


  int getPosition(double knot);

};

} // namespace Go

#endif
