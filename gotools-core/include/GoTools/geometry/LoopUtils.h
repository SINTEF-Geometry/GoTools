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

#ifndef _LOOPUTILS_H
#define _LOOPUTILS_H


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"

#include <vector>
#include <memory>

namespace Go{

/// Functions for checking the orientation of loops (closed curves), and 
/// whether one loop on a surface encloses another.
namespace LoopUtils {
    /// Represent a vector of ParamCurve curves as a vector of CurveOnSurface curves
    /// The surface is given as additional input
    /// Note that the function throws if the surface information is inconsistent
    void representAsSurfaceCurves(const std::vector< shared_ptr<ParamCurve> >& curves,
				  shared_ptr<BoundedSurface> surf,
				  std::vector<shared_ptr<CurveOnSurface> >& cvs_on_sf);

    /// Check if a closed 2D-loop is oriented counterclockwise or not.
    /// \param simple_par_loop a sequence of 2D curves that are joined start-to-end and that form
    ///                        a closed loop in the plane.
    /// \param int_tol (geometric) tolerance used for internal computations (intersection detections)
    /// \return 'true' if the loop was found to be oriented CCW, otherwise 'false'.
    bool loopIsCCW(const std::vector<shared_ptr<Go::ParamCurve> >& simple_par_loop, 
		   double space_epsilon, double int_tol);

    /// To support previous interface
    bool loopIsCCW(const std::vector<shared_ptr<Go::SplineCurve> >& simple_par_loop, 
		   double space_epsilon, double int_tol);

    /// Check if a loop defined by CurveOnSurface s is oriented counterclockwise in the surface's 
    /// parametric domain.
    /// \param loop a sequence of CurveOnSurface s that are jointed start-to-end and that form
    ///             a closed loop on the surface.  
    /// \param int_tol (geometric) tolerance used for internal computations (intersection detections)
    /// \return 'true' if the loop was found to be oriented CCW, otherwise 'false'.
    bool
    paramIsCCW(const std::vector< shared_ptr<Go::CurveOnSurface> >& loop,
	       double space_epsilon, double int_tol);

    /// Check if a closed 2D-loop is oriented counterclockwise or not. The
    /// loop is given as a CurveLoop.
    bool loopIsCCW(const CurveLoop& loop, double int_tol);
    
    /// Loops expected to be disjoint, except possibly share part of boundary.

    /// Test whether one loop lies entirely within another.  This function does not work in the 
    /// general case; it makes the assumption that the loops do NOT intersect each other 
    /// transversally (their boundaries are allowed to tangentially touch though).  The algorithm
    /// works by testing a single point on the first loop for being inside the second loop, so 
    /// if the first loop lay partially inside, partially outside the second, the answer would be
    /// arbitrary.
    /// \param first_loop the first loop
    /// \param second_loop the second loop
    /// \param loop_tol the tolerance for defining coincidence between start/endpoints on the 
    ///                 consecutive curve segments consituting a loop.
    /// \param int_tol tolerance used for intersection calculations
    /// \return 'true' if 'first_loop' was found to be located inside 'second_loop' (given the 
    ///         assumptions above).  'false' otherwise.
    bool firstLoopInsideSecond(const std::vector<shared_ptr<Go::CurveOnSurface> >& first_loop,
			       const std::vector<shared_ptr<Go::CurveOnSurface> >& second_loop,
			       double loop_tol, double int_tol);

    /// Reorganize curves to create a ccw loop
    /// The curves are assumed to have correct sequence, but possible
    /// wrong parameter direction
    /// Return value: false = reorganization not possible
    bool makeLoopCCW(std::vector<shared_ptr<ParamCurve> >& loop_cvs,
		     double tol);

} // end namespace Go
} // end namespace LoopUtils

#endif // _LOOPUTILS_H

