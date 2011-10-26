//===========================================================================
//                                                                           
// File: GSSinsertKnot.C                                                     
//                                                                           
// Created: Wed Apr  4 15:09:53 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: GSSinsertKnot.C,v 1.5 2003-05-08 15:37:13 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"

using namespace Go;


// @@ These functions are quite slow and may be optimized a lot!
//
// The basic idea is to treat the surface as a curve (may have to
// turn the surface, though)
//

//===========================================================================
void SplineSurface::insertKnot_v(double apar)
//===========================================================================
{
    int kdim = rational_ ? dim_+1 : dim_;
    // Make a hypercurve from this surface
    SplineCurve cv(numCoefs_v(), order_v(), basis_v_.begin(),
		     activeCoefs().begin(), kdim*numCoefs_u(), false);
    // Insert the knot in the curve
    cv.insertKnot(apar);
    // Insert knot into basis
    basis_v_.insertKnot(apar);
    // Copy back the coefficients
    if (!rational_) {
	coefs_.resize(numCoefs_u()*numCoefs_v()*dim_);
	std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
    } else {
	rcoefs_.resize(numCoefs_u()*numCoefs_v()*(dim_+1));
	std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
	updateCoefsFromRcoefs();
    }
}

//===========================================================================
void SplineSurface::insertKnot_v(const std::vector<double>& new_knots)
//===========================================================================
{
    int kdim = rational_ ? dim_+1 : dim_;
    // Make a hypercurve from this surface
    SplineCurve cv(numCoefs_v(), order_v(), basis_v_.begin(),
		     activeCoefs().begin(), kdim*numCoefs_u(), false);
    // Insert the knot in the curve
    cv.insertKnot(new_knots);
    // Insert knot into basis
    basis_v_.insertKnot(new_knots);
    // Copy back the data
    if (!rational_) {
	coefs_.resize(numCoefs_u()*numCoefs_v()*dim_);
	std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
    } else {
	rcoefs_.resize(numCoefs_u()*numCoefs_v()*(dim_+1));
	std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
	updateCoefsFromRcoefs();
    }

}

//===========================================================================
void SplineSurface::insertKnot_u(double apar)
//===========================================================================
{
    swapParameterDirection();
    insertKnot_v(apar);
    swapParameterDirection();
}

//===========================================================================
void SplineSurface::insertKnot_u(const std::vector<double>& new_knots)
//===========================================================================
{
    swapParameterDirection();
    insertKnot_v(new_knots);
    swapParameterDirection();
}
