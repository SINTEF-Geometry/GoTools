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

#ifndef _HERMITEAPPREVALSURF_H
#define _HERMITEAPPREVALSURF_H


#include "GoTools/creators/EvalSurface.h"
#include "GoTools/creators/HermiteGrid2D.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"


namespace Go
{

/// This class is used to generate a SplineSurface from a EvalSurface
/// using Hermite interpolation.  The generated surface will approximate the
/// EvalSurface within specified tolerances.

class HermiteApprEvalSurf
{
public:

    /// Constructor where the tolerances and the surface to approximate are specified.
    /// \param sf the surface that we want to generate a Hermite approximation of.
    ///            The surface is \em not copied, only pointed to by the HermiteApprEvalSurf.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalSurfaces
    HermiteApprEvalSurf(EvalSurface* sf, double tolerance1, double tolerance2);

    /// Constructor where the tolerances and the surface to approximate are specified,
    /// as well as the parameters for which we will sample the input surface before
    /// starting the approximating process.
    /// \param sf the surface that we want to generate a Hermite approximation of.
    ///            The surface is \em not copied, only pointed to by the HermiteApprEvalSurf.
    /// \param initpars pointer to the array of parameter values for which we will
    ///        sample the input surface.
    /// \param n number of parameter values in the array 'initpars'.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalSurfaces
    HermiteApprEvalSurf(EvalSurface* sf,
                        double initpars_u[],
                        int mm,
                        double initpars_v[],
                        int nn,
                        double tolerance1, double tolerance2);

    /// Empty destructor
    ~HermiteApprEvalSurf(){}

    /// Refine the internal sampling of the surface to approximate such that
    /// the Hermite interpolated surface of this sampling approximates parametrically
    /// the original surface within the specified tolerance.
    void refineApproximation();	// Refine Hermite Grid.

    /// Return the cubic spline surface Hermite interpolating the grid.
    /// \return the cubic spline surface that Hermite interpolates the sampled 
    ///         points of the EvalSurface specified in the constructor.
    shared_ptr<SplineSurface> getSurface(bool& method_failed);

    const HermiteGrid2D& getGrid() const;

    // To be called after refinement is completed.
    void removeGridLines(const std::vector<double>& grid_lines_u,
                         const std::vector<double>& grid_lines_v);

    // We do not split closer than no_split_dist to a curve in 2d_cvs.
    void setNoSplit(const std::vector<shared_ptr<SplineCurve> >& no_split_cvs_2d,
                    double no_split_dist_2d);
    
 private:
    EvalSurface* surface_;	// Pointer to original surface existing outside *this.
    const double tol1_; // Used by surface_ in approximationOK().
    const double tol2_; // Used by surface_ in approximationOK().
    double min_interval_;	// Smaller intervals are not refined
    HermiteGrid2D grid_;
//     shared_ptr<SplineSurface> surface_approx_; // Spline representation of approximation
    
#if 1
    // The knot values of grid elements which are not to be used. This approach is easier to extend and
    // maintain than using the grid index (whihc varies as the grid is extended).
    std::vector<double> removed_grid_u_;
    std::vector<double> removed_grid_v_;
#else
    std::vector<int> removed_grid_u_;
    std::vector<int> removed_grid_v_;
#endif

    // We may restrict the split parameters in the grid.

    // @@sbr 201705 Should this be part of the HermiteGrid2D? Along with the removed_grid_u_ &
    // removed_grid_v_? Seems like the better choice for an easier data structure.
    std::vector<shared_ptr<SplineCurve> > no_split_cvs_2d_;
    std::vector<BoundingBox> no_split_bd_box_; // We're representing the cvs with an axis linear bounding box.
    double no_split_2d_tol_;

        // Distance to evaluator ok (with current grid)?
    // Return value: 1 = ok, 0 = insert new_knot, -1 = failed.
    int testSegment(int left1, int left2, double& new_knot, bool& dir_is_u);
    int bisectSegment(int left1, int left2, bool& dir_is_u);
    bool method_failed_;

    int splitDomain(double spar1, double epar1, double spar2, double epar2, Point bezcoef[16],
                    bool& dir_is_u, double& new_knot);
    
};

} // namespace Go;


#endif // _HERMITEAPPREVALSURF_H

