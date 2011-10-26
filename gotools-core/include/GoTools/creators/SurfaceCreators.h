//===========================================================================
//                                                                           
// File: SurfaceCreators.h                                                 
//                                                                           
// Created: Thu Feb 21 09:32:50 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: SurfaceCreators.h,v 1.7 2007-12-04 16:11:39 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
    std::shared_ptr<SplineSurface>
    createSmoothTransition(const std::vector<std::shared_ptr<const ParamSurface> >& surfs,
			   const std::vector<std::shared_ptr<const CurveOnSurface> >& int_cvs,
			   double dist_0, double dist_1, double epsge,
			   std::vector<std::shared_ptr<SplineCurve> >& trim_crvs);

    /// Return the product of the two spline surfaces. Expecting them to be
    /// non-rational.
    /// \param sf1 the first surface.
    /// \param sf2 the second surface.
    /// \return the surface product.
    std::shared_ptr<SplineSurface> mult1DSurfaces(const SplineSurface& sf1,
						    const SplineSurface& sf2);

    /// Return the product of the two Bezier patches. Expecting them to be
    /// non-rational.
    /// \param patch1 SplineSurface of Bezier type (i.e. no inner knots).
    /// \param patch2 SplineSurface of Bezier type (i.e. no inner knots).
    /// \return the surface product.
    std::shared_ptr<SplineSurface> mult1DBezierPatches(const SplineSurface& patch1,
							 const SplineSurface& patch2);

    /// Return the rational surface 'nom_sf/den_sf'.
    /// \param nom_sf the nominator surface.
    /// \param den_sf the denominator surface.
    /// \param weights_in_first true if the coefficients of the nom_sf have been multiplied
    ///                         by the corresonding rational coefficients in den_sf.
    /// \return the rational surface.
    std::shared_ptr<Go::SplineSurface> mergeRationalParts(const Go::SplineSurface& nom_sf,
							    const Go::SplineSurface& den_sf,
							    bool weights_in_first = false);

    /// Given input of 1d-sf, return the 3d visualization (u, v, f(u, v)).
    /// \param sf_1d 1-dimensional surface.
    /// \return the 3-dimensional surface.
    std::shared_ptr<Go::SplineSurface> insertParamDomain(const Go::SplineSurface& sf_1d);

    /// Assuming the sf is rational, separate the geometric space from the homogenuous.
    /// \param sf the input spline surface.
    /// \return the non-rational parts of sf.
    std::vector<std::shared_ptr<Go::SplineSurface> >
    separateRationalParts(const Go::SplineSurface& sf);

} // end of namespace SurfaceCreators

} // end of namespace Go


#endif // _SURFACECREATORS_H

