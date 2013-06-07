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

#ifndef _SURFACECREATORS_H
#define _SURFACECREATORS_H


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"


namespace Go
{

/// Various functions for generating SplineSurface s by approximation, blending,
/// etc.
namespace SurfaceCreators
{

    /// Given input of two parametric surfaces, which intersect along an edge,
    /// create a smooth surface by moving boundary of surfs to inner part of
    /// current surfs. Normals of input surfs are assumed to be consistent.
    /// \param surfs the input surfaces, size of vector is 2.
    /// \param int_cvs intersection curves between the two surfaces.  Vector has size 2.
    ///                Both parameter curves should exist, the space curve should be the same.
    /// \param dist_0 the offset space distance in surfs[0].
    /// \param dist_1 the offset space distance in surfs[1].
    /// \param epsge the geometrical tolerance when offsetting.
    /// \param trim_crvs the offset trim curves: par_cv0, space_cv0, par_cv1, space_cv1.
    shared_ptr<SplineSurface>
    createSmoothTransition(const std::vector<shared_ptr<const ParamSurface> >& surfs,
			   const std::vector<shared_ptr<const CurveOnSurface> >& int_cvs,
			   double dist_0, double dist_1, double epsge,
			   std::vector<shared_ptr<SplineCurve> >& trim_crvs);

    /// Return the product of the two spline surfaces. Expecting them to be
    /// non-rational.
    /// \param sf1 the first surface.
    /// \param sf2 the second surface.
    /// \return the surface product.
    shared_ptr<SplineSurface> mult1DSurfaces(const SplineSurface& sf1,
						    const SplineSurface& sf2);

    /// Return the product of the two Bezier patches. Expecting them to be
    /// non-rational.
    /// \param patch1 SplineSurface of Bezier type (i.e. no inner knots).
    /// \param patch2 SplineSurface of Bezier type (i.e. no inner knots).
    /// \return the surface product.
    shared_ptr<SplineSurface> mult1DBezierPatches(const SplineSurface& patch1,
							 const SplineSurface& patch2);

    /// Return the rational surface 'nom_sf/den_sf'.
    /// \param nom_sf the nominator surface.
    /// \param den_sf the denominator surface.
    /// \param weights_in_first true if the coefficients of the nom_sf have been multiplied
    ///                         by the corresonding rational coefficients in den_sf.
    /// \return the rational surface.
    shared_ptr<Go::SplineSurface> mergeRationalParts(const Go::SplineSurface& nom_sf,
							    const Go::SplineSurface& den_sf,
							    bool weights_in_first = false);

    /// Given input of 1d-sf, return the 3d visualization (u, v, f(u, v)).
    /// \param sf_1d 1-dimensional surface.
    /// \return the 3-dimensional surface.
    shared_ptr<Go::SplineSurface> insertParamDomain(const Go::SplineSurface& sf_1d);

    /// Assuming the sf is rational, separate the geometric space from the homogenuous.
    /// \param sf the input spline surface.
    /// \return the non-rational parts of sf.
    std::vector<shared_ptr<Go::SplineSurface> >
    separateRationalParts(const Go::SplineSurface& sf);

} // end of namespace SurfaceCreators

} // end of namespace Go


#endif // _SURFACECREATORS_H

