//==========================================================================
//                                                                          
// File: BernsteinUtils.h                                                    
//                                                                          
// Created: Wed Feb 26 13:13:18 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: BernsteinUtils.h,v 1.12 2006-04-07 13:47:53 afr Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _BERNSTEINUTILS_H
#define _BERNSTEINUTILS_H


#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/implicitization/BernsteinMulti.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <vector>


namespace Go {


class SplineCurve;
class SplineSurface;


/// Takes a segment (defined as having numCoefs() == order() and an
/// order()-regular knot vector) and returns a bernstein polynomial
/// equal to the curve's coordinate number dd. Allowed values for dd
/// are in [0, dimension()-1] for nonrational segments, and in [0,
/// dimension()] for rational (NURBS) segments. The polynomials
/// returned in the rational case are the weighted coordinate
/// functions for dd in [0, dimension()-1] and to the weight
/// polynomial for dd equal to dimension().
void spline_to_bernstein(const SplineCurve& seg, int dd,
			 BernsteinPoly& bp);

/// Converts the dd-component of a surface patch on SplineSurface form
/// to a BernsteinMulti
void spline_to_bernstein(const SplineSurface& pat, int dd,
			 BernsteinMulti& bm);

/// Converts a curve segment on SplineCurve form to a vector of
/// BernsteinPolys
void spline_to_bernstein(const SplineCurve& seg,
			 std::vector<BernsteinPoly>& seg_bp);


/// Converts a surface patch on SplineSurface form to a vector of
/// BernsteinMultis
void spline_to_bernstein(const SplineSurface& pat,
			 std::vector<BernsteinMulti>& pat_bm);


/// Converts a Bezier curve segment on SplineCurve form to an Array of
/// \a Ndim BernsteinPolys. If the SplineCurve is rational, \a Ndim is
/// equal to \a dim+1, where \a dim is the dimension of
/// space. Furthermore, in \a curve_bp, the weights are included in
/// the space part of the coefficients, i.e. they have the form \a
/// wP1, \a wP2, ..., \a wPdim, \a w.
/// \param segment the Bezier curve segment on SplineCurve form
/// \param curve_bp the resulting Array of BernsteinPolys
template <int Ndim>
void splineToBernstein(const SplineCurve& segment,
		       Array<BernsteinPoly, Ndim>& curve_bp)
{
    int order = segment.order();
    ALWAYS_ERROR_IF(order != segment.numCoefs(),
		    "The curve is not Bezier");
    int dim = segment.dimension();
    bool rational = segment.rational();
    int ndim = rational ? dim+1 : dim;
    ALWAYS_ERROR_IF(ndim != Ndim,
		    "The dimension of curve_bp does not match the "
		    "dimension of seg");
    std::vector<double> coefs(order);
    for (int i = 0; i < Ndim; ++i) {
	curve_bp[i] = BernsteinPoly(coefs);
    }
    std::vector<double>::const_iterator segment_coefs;
    if (rational) {
	segment_coefs = segment.rcoefs_begin();
    } else {
	segment_coefs = segment.coefs_begin();
    }
    for (int i = 0; i < order; ++i) {
	for (int dd = 0; dd < Ndim; ++dd) {
	    curve_bp[dd][i] = segment_coefs[Ndim*i + dd];
	}
    }

    return;
}


/// Converts a Bezier surface patch on SplineSurface form to an Array
/// of \a Ndim BernsteinMultis. If the SplineSurface is rational, \a
/// Ndim is equal to \a dim+1, where \a dim is the dimension of
/// space. Furthermore, in \a surface_bm, the weights are included in
/// the space part of the coefficients, i.e. they have the form \a
/// wP1, \a wP2, ..., \a wPdim, \a w.
/// \param patch the Bezier surface on SplineSurface form
/// \param surface_bm the resulting Array of BernsteinMultis
template <int Ndim>
void splineToBernstein(const SplineSurface& patch,
		       Array<BernsteinMulti, Ndim>& surface_bm)
{
    int orderu = patch.order_u();
    int orderv = patch.order_v();
    ALWAYS_ERROR_IF((orderu != patch.numCoefs_u())
		    || (orderv != patch.numCoefs_v()),
		    "The surface is not Bezier");
    int dim = patch.dimension();
    bool rational = patch.rational();
    int ndim = rational ? dim+1 : dim;
    ALWAYS_ERROR_IF(ndim != Ndim,
		    "The dimension of surface_bm does not match the "
		    "dimension of pat");
    std::vector<double> coefs(orderu * orderv);
    for (int i = 0; i < Ndim; ++i) {
	surface_bm[i] = BernsteinMulti(orderu-1, orderv-1, coefs);
    }
    std::vector<double>::const_iterator patch_coefs;
    if (rational) {
	patch_coefs = patch.rcoefs_begin();
    } else {
	patch_coefs = patch.coefs_begin();
    }
    for (int j = 0; j < orderv; ++j) {
	for (int i = 0; i < orderu; ++i) {
	    for (int dd = 0; dd < Ndim; ++dd) {
		surface_bm[dd][orderu*j + i]
		    = patch_coefs[Ndim*(orderu*j + i) + dd];
	    }
	}
    }

    return;
}


/// Converts an Array of \a Ndim BernsteinPolys to a SplineCurve
/// segment. If the wanted SplineCurve is rational, \a Ndim is equal
/// to \a dim+1, where \a dim is the dimension of space. Furthermore,
/// in \a curve_bp, the weights are included in the space part of the
/// coefficients, i.e. they have the form \a wP1, \a wP2, ..., \a
/// wPdim, \a w.
/// \param curve_bp an Array of \a Ndim BernsteinPolys
/// \param rational flag to indicate if the objects are rational
/// \param segment the resulting spline curve segment
template <int Ndim>
void bernsteinToSpline(const Array<BernsteinPoly, Ndim>& curve_bp,
		       bool rational,
		       SplineCurve& segment)
{
    int npoints = curve_bp[0].degree() + 1;
    int order = npoints;
    int dim = (rational ? Ndim - 1 : Ndim);
    std::vector<double> knots(2 * order, 0.0);
    for (int i = order; i < 2 * order; ++i)
	knots[i] = 1.0;
    std::vector<double> coefs(Ndim * npoints);
    for (int j = 0; j < npoints; ++j) {
	for (int dd = 0; dd < dim; ++dd) {
	    coefs[j*Ndim + dd] = curve_bp[dd][j];
	}
	if (rational) {
	    coefs[j*Ndim + dim] = curve_bp[Ndim-1][j];
	}
    }

    segment = SplineCurve(npoints, order,
			  knots.begin(), coefs.begin(), dim, rational);
    return;
}


/// Converts an Array of \a Ndim BernsteinMultis to a SplineSurface
/// patch. If the wanted SplineSurface is rational, \a Ndim is equal
/// to \a dim+1, where \a dim is the dimension of space. Furthermore,
/// in \a surface_bm, the weights are included in the space part of
/// the coefficients, i.e. they have the form \a wP1, \a wP2, ..., \a
/// wPdim, \a w.
/// \param surface_bm an Array of \a Ndim BernsteinMultis
/// \param rational flag to indicate if the objects are rational
/// \param patch the resulting spline surface patch
template <int Ndim>
void bernsteinToSpline(const Array<BernsteinMulti, Ndim>& surface_bm,
		       bool rational,
		       SplineSurface& patch)
{
    int npointsu = surface_bm[0].degreeU() + 1;
    int npointsv = surface_bm[0].degreeV() + 1;
    int orderu = npointsu;
    int orderv = npointsv;
    int dim = (rational ? Ndim - 1 : Ndim);
    std::vector<double> knotsu(2 * orderu, 0.0);
    for (int i = orderu; i < 2 * orderu; ++i)
	knotsu[i] = 1.0;
    std::vector<double> knotsv(2 * orderv, 0.0);
    for (int i = orderv; i < 2 * orderv; ++i)
	knotsv[i] = 1.0;
    std::vector<double> coefs(Ndim * npointsu * npointsv);
    for (int j = 0; j < npointsv; ++j) {
	for (int i = 0; i < npointsu; ++i) {
	    for (int dd = 0; dd < dim; ++dd) {
		coefs[(npointsu*j + i) * Ndim + dd]
		    = surface_bm[dd][npointsu*j + i];
	    }
	    if (rational) {
		coefs[(npointsu*j + i) * Ndim + dim]
		    = surface_bm[Ndim-1][npointsu*j + i];
	    }
	}
    }

    patch = SplineSurface(npointsu, npointsv, orderu, orderv,
			  knotsu.begin(), knotsv.begin(),
			  coefs.begin(), dim, rational);
    return;
}


} // namespace Go


#endif // _BERNSTEINUTILS_H

