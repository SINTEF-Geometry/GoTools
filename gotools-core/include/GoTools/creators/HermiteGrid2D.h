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

#ifndef _HERMITEGRID2D_H
#define _HERMITEGRID2D_H


#include <vector>
#include "GoTools/utils/Point.h"

namespace Go
{

class EvalSurface;

/// The type "GoHemiteGrid1D" holds a one dimensional grid containing sampled
/// points and derivatives from a curve.  It can be used to generate bezier
/// curve segments obtained by Hermite interpolation of intervals between
/// the sampled parameter values.

class HermiteGrid2D
{
public:

    /// Construct a HermiteGrid1D from a curve and a given interval.
    /// The grid will only contain the sampled values (position, derivative)
    /// for the first and last value of the interval.
    /// \param crv curve to sample
    /// \param start start of interval
    /// \param end end of interval
    HermiteGrid2D(const EvalSurface& sf,
                  double u1, double u2, double v1, double v2);

    /// Construct a HermiteGrid1D from a curve and a set of parameter
    /// values.
    /// \param crv curve to sample
    /// \param param array of strictly increasing parameters contained in
    ///              the parameter domain of 'crv'.
    /// \param n number of elements in 'param[]'.
    HermiteGrid2D(const EvalSurface& sf,
                  double param_u[], double param_v[], int mm, int nn);
    
    /// Default destructor.
    ~HermiteGrid2D();

    /// Add another sample (parameter, position, tangent) to the grid.
    /// Returns the index of the new knot (parameter value) in the sorted
    /// knot vector after insertion.  
    /// \param crv curve to evaluate
    /// \param knot the new sample value (parameter value, knot)
    int addKnot(const EvalSurface& sf, double knot, bool dir_is_u);
  
    /// Calculate Bezier coefficients of the cubic curve interpolating 
    /// the point and tangent values at grid nodes with indices "left" 
    /// and "right".
    /// \param left  indicating grid node for start of curve segment
    /// \param right indicating grid node for end of curve segment
    /// \param spar  start parameter of segment
    /// \param epar  end parameter of segment
    /// \param bezcoef array of cubic Bezier coefficients
    void getSegment(int left1, int right1,
                    int left2, int right2,
                    double& spar1, double& epar1,
                    double& spar2, double& epar2,
                    Point bezcoef[4]);
    
    /// Return the grid parameters
    std::vector<double> getKnots(bool dir_is_u) const { return (dir_is_u ? knots_u_ : knots_v_); }

    /// Return the sample values (positions and first derivatives)
    std::vector<Point> getData() const { return array_; }

    /// Return the spatial dimension
    int dim() const { return dim_; }

    /// Return the number of samples in the grid
    int size1() const {return MM_;}

    int size2() const {return NN_;}

    void removeGridLines(const std::vector<int>& grid_lines_u,
                         const std::vector<int>& grid_lines_v);

    // Get the no split status for the element.
    int getNoSplitStatus(int ind_u, int ind_v);

    // The split status is replaced (not added).
    void setNoSplitStatus(int ind_u, int ind_v, int no_split_status);

private:
    std::vector<double> knots_u_;     // Sorted array of DISTINCT parameters of sf in u-dir.
    std::vector<double> knots_v_;     // Sorted array of DISTINCT parameters of sf in v-dir.
    std::vector<Point> array_; 	// Array holding position, and
  				// directional derivative, including twist.
    int dim_;        		// Spatial dimension of position,
    int MM_;         		// Number of grid points in u dir.
    int NN_;         		// Number of grid points in v dir.
    int elem_size_;		// Number of Point stored for each
  				// grid point (typically 3: pt, der_u, der_v)
    int index_u_;                   // Index into knot array
    int index_v_;                   // Index into knot array

    std::vector<int> removed_grid_u_;
    std::vector<int> removed_grid_v_;

    std::vector<int> no_split_status_; // Size MM*NN. 0 => no restrictions, 1 => not in dir 1 (u),
                                       // 2 => not in dir 2 (v), 3 => no split.
    
    int getPosition(double knot, bool dir_is_u);

};

} // namespace Go


#endif // _HERMITEGRID2D_H

