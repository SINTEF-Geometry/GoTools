//===========================================================================
//                                                                           
// File: LoopUtils.h                                                         
//                                                                           
// Created: Tue Jun  3 15:07:01 2003                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: LoopUtils.h,v 1.11 2008-12-01 14:00:48 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
    void representAsSurfaceCurves(std::vector< shared_ptr<ParamCurve> >& curves,
				  shared_ptr<BoundedSurface> surf,
				  std::vector<shared_ptr<CurveOnSurface> >& cvs_on_sf);

    /// Check if a closed 2D-loop is oriented counterclockwise or not.
    /// \param simple_par_loop a sequence of 2D curves that are joined start-to-end and that form
    ///                        a closed loop in the plane.
    /// \param int_tol (geometric) tolerance used for internal computations (intersection detections)
    /// \return 'true' if the loop was found to be oriented CCW, otherwise 'false'.
    bool loopIsCCW(const std::vector<shared_ptr<Go::SplineCurve> >&
		   simple_par_loop, 
		   double int_tol);

    /// Check if a loop defined by CurveOnSurface s is oriented counterclockwise in the surface's 
    /// parametric domain.
    /// \param loop a sequence of CurveOnSurface s that are jointed start-to-end and that form
    ///             a closed loop on the surface.  
    /// \param int_tol (geometric) tolerance used for internal computations (intersection detections)
    /// \return 'true' if the loop was found to be oriented CCW, otherwise 'false'.
    bool
    paramIsCCW(const std::vector< shared_ptr<Go::CurveOnSurface> >& loop,
	       double int_tol);

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

} // end namespace Go
} // end namespace LoopUtils

#endif // _LOOPUTILS_H

