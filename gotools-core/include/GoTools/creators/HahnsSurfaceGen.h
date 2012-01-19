//===========================================================================
//                                                                           
// File: HahnsSurfaceGen.h                                                      
//                                                                           
// Created: Wed Dec 12 12:26:29 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: HahnsSurfaceGen.h,v 1.5 2007-12-04 16:11:39 jbt Exp $
//                                                                           
// Description: Utility functions for generating a surface when given bnd curves.
//              Cross tangents may be given as well. Assumes 3 <= #curves <= 6.
//              Code ported from SISL (written by vsk). Uses Hahn's method:
//              Joerg Hahn: "Filling Polygonial Holes with Rectangular Patches",
//              Theory and Practice of Geometric Modelling, Blackburn Oct 1988.
//                                                                           
//===========================================================================

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
