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

#ifndef _SURFACEGEN_H
#define _SURFACEGEN_H

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/config.h"


namespace Go{

/// This namespace contains functions used to create a number of surface covering
/// a whole defined by of a number of boundary curves (3<=nmb_cvs<=6). Hahns
/// method is used in the construction.
/// Typically used for patching a model in which one or more surfaces are missing.
namespace HahnsSurfaceGen{


    /// The routine which actually creates vector of blending surfaces, with given
    /// curves as their total bnd curve. # of return surfaces equals # of curves.
    /// All bnd curves share orientation (clockwise), and all cross curves point
    /// inwards. A missing cross curve is indicated by a 0 pointer.
    /// bnd_curves must form a loop. As we approximate modified cross curves,
    /// we include parameters to denote the exactness of the approximation.
    /// All curves expected to live in 3-dimensional space.
    /// \param bnd_curves edge curves for the Hahns Surface.
    /// \param cross_curves corresponding cross tangent curves.
    /// \param neighbour_tol allowed distance between corresponding end points.
    /// \param kink_tol allowed angle between corresponding tangents.
    /// \param knot_diff_tol parametric tolerance for equality of knots.
    /// \return vector containing the Hahns Surface.
    std::vector<shared_ptr<Go::ParamSurface> >
    constructPolygonialSurface(std::vector<shared_ptr<Go::ParamCurve> >& bnd_curves,
			       std::vector<shared_ptr<Go::ParamCurve> >& cross_curves,
			       double neighbour_tol,
			       double kink_tol,
			       double knot_diff_tol);

    /// Given input of bnd_curves and cross_tangent curves, create the corresponding
    /// Hahns Surfaces.  Assuming input curves fulfill corner conditions.
    /// All curves expected to live in 3-dimensional space.
    /// \param bnd_curves edge curves for the Hahns Surface.
    /// \param mod_cross_curves corresponding cross tangent curves.
    /// \param neighbour_tol allowed distance between corresponding end points.
    /// \param kink_tol allowed angle between corresponding tangents.
    /// \param knot_diff_tol parametric tolerance for equality of knots.
    /// \return vector containing the Hahns Surface.
    std::vector<shared_ptr<Go::ParamSurface> >
    constructHahnsSurface(std::vector<shared_ptr<Go::SplineCurve> >& bnd_curves,
			  std::vector<shared_ptr<Go::SplineCurve> >& mod_cross_curves,
			  double neighbour_tol,
			  double kink_tol,
			  double knot_diff_tol);

} // end of namespace

} // end of namespace


#endif // _SURFACEGEN_H
