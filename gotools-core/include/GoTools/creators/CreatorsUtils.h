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

#ifndef _CREATORSUTILS_H
#define _CREATORSUTILS_H


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/BoundedSurface.h"
#include <memory>
#include "GoTools/utils/config.h"

namespace Go {

/// Related to the generation of cross tangent curves.
namespace CreatorsUtils
{

    /// Helper function. Returns the parametric representation of the input
    /// curve.
    /// The returned curve is created (i.e. it should be be handled by a
    /// smart ptr)!
    /// This function takes as argument a vector of (shared pointers to)
    /// CurveOnSurfaces, which
    /// are assumed to lie on the same surface and to represent consecutive
    /// parts of a larger
    /// curve, and return the parametric representation of this larger curve.
    /// \param cv a vector of (shared pointers to) CurveOnSurfaces.
    ///           We suppose that these curves
    ///           are all lying on the same surface and that they constitute
    ///           consecutive parts
    ///           of a global curve.  
    /// \return if successful, a pointer to the newly generated SplineCurve,
    ///         which represents
    ///         the global curve that is the union of the input curves.
    ///         This SplineCurve is 
    ///         represented in the parameter plane of the surface, and is
    ///         therefore a 2D curve.
    ///         If it could not be constructed (due a failure of any of
    ///         the provided CurveOnSurfaces
    ///         to provide their internal parameter curve), a null pointer
    ///         is returned instead.
    ///         The user assumes ownership of the SplineCurve.
    SplineCurve GO_API *getParametricCurve
    (const std::vector<shared_ptr<const CurveOnSurface> >& cv);

    /// Generate the inwards pointing cross-tangent curve along the input trim
    /// curve.
    /// \param cv the trim curve.
    /// \return the inwards cross tangent curve.
    shared_ptr<Go::SplineCurve> GO_API
    createCrossTangent(const Go::CurveOnSurface& cv);

    /// The cross tangent cv along input cv is created. Length given by length
    /// of derivs in
    /// sf along cv or by input_cross_cv if it exists. Direction of created
    /// cross tangent cv
    /// is to the left of cv (i.e. it should point into the surface).
    /// If input_cross_cv exists so does basis_space_cv and they must share
    /// parametrization.
    /// If input_cross_cv exists, but not basis_space_cv, it shares
    /// parametrization with cv.
    /// \param cv the trim curve along which we are to compute the cross
    ///        tangent cv.
    /// \param basis_space_cv if != NULL we use the basis and parametrization
    ///                       of basis_space_cv.
    /// \param cross_cv_ref if != NULL the curve defines the angle between
    ///                     the boundary curve and the new cross tangent
    ///                     curve.
    ///                     Additionally it defines the length of the new
    ///                     cross tangent curve.
    /// \param appr_offset_cv whether the method should approximate the offset
    ///                       curve
    ///                       (as opposed to the cross tangent curve).
    /// \return the cross tangent (or offset) curve.
    shared_ptr<Go::SplineCurve> GO_API
    createCrossTangent(const Go::CurveOnSurface& cv,
		       shared_ptr<Go::SplineCurve> basis_space_cv,
		       const Go::SplineCurve* cross_cv_ref,
		       bool appr_offset_cv = true);

    /// Project a point in space onto a surface. If the surface is
    /// closed, and the point is on the seam, then the set of possible
    /// candidates is returned.
    /// \param sf the surface the point is projected onto
    /// \param closed_dir_u boolean that is \c true if the surface is closed
    /// in the \f$u\f$-direction
    /// \param closed_dir_v boolean that is \c true if the surface is closed
    /// in the \f$v\f$-direction
    /// \param space_pt the point to be projected
    /// \return vector of candidates of projected points in the form
    /// of parameter pairs
    std::vector<Go::Point> GO_API
    projectPoint(const Go::ParamSurface* sf,
		 bool closed_dir_u, bool closed_dir_v,
		 const Go::Point& space_pt, double epsgeo = 1e-04);

    shared_ptr<Go::Point>
    projectCurvePoint(const ParamSurface* sf,
		      bool closed_dir_u, bool closed_dir_v,
		      const Go::ParamCurve* space_cv, double cv_par, double epsgeo = 1e-04);

    /// Project a point on a space curve onto a surface. If the
    /// surface is closed, and the curve follows the seam at the point
    /// in question, then the point corresponding to a
    /// counterclockwise curve is chosen. If the programmer knows
    /// that this is not the case, then this must be handled.
    /// \param sf the surface the point is projected onto
    /// \param closed_dir_u boolean that is \c true if the surface is closed
    /// in the \f$u\f$-direction
    /// \param closed_dir_v boolean that is \c true if the surface is closed
    /// in the \f$v\f$-direction
    /// \param space_cv the curve to be projected
    /// \param cv_par the parameter of the curve point to be projected
    /// \return the projected point in the form of a parameter pair
    shared_ptr<Go::Point>
    projectCurvePoint(const SplineSurface& sf,
		      bool closed_dir_u, bool closed_dir_v,
		      const Go::ParamCurve* space_cv, double cv_par, double epsgeo = 1e-04);

    /// Repair erranous seem curves in bounded surfaces with
    /// closed underlying surfaces
    void GO_API
      fixSeemCurves(shared_ptr<BoundedSurface> bd_sf, 
		    std::vector<shared_ptr<CurveOnSurface> >& loop_cvs,
		    bool closed_dir_u, bool closed_dir_v,
		    double tol);

    /// Repair erranous trimming curves in bounded surfaces and compute
    /// missing parameter curves in the corresponding CurveOnSurface curves
    // The user is given the option to overrun the default tolerance epsgeo
    // from the boundary_loops with epsgeo*epsgeo_frac.
    void GO_API
    fixTrimCurves(shared_ptr<Go::BoundedSurface> bd_sf,
		  double epsgeo_frac = 1.0, double tol = 1.0e-3,
		  double tol2 = 1.0e-2, double ang_tol = 1.0e-2);

} // of namespace CreatorsUtils.

}; // end namespace Go



#endif // _CREATORSUTILS_H
