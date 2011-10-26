//===========================================================================
//                                                                           
// File: GSSraiseOrder.C                                                     
//                                                                           
// Created: Thu Aug 16 18:39:29 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: GSSraiseOrder.C,v 1.5 2005-10-18 15:11:05 sbr Exp $
//                                                                           
// Description: Raise order of spline surface. Method: treat surface as a curve.
//                                                                           
//===========================================================================

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"


using namespace Go;


//===========================================================================
void SplineSurface::raiseOrder(int raise_u, int raise_v)
//===========================================================================
{
    ALWAYS_ERROR_IF(raise_u < 0 || raise_v < 0,
		    "Order to raise by must be positive!");

    // We're raising in the u-direction first.
    bool rat = rational_;
    std::vector<double>::iterator coef_iter = (rat) ? rcoefs_.begin() : coefs_.begin();
    int rdim = (rat) ? dimension() + 1 : dimension();
    if (raise_u > 0) {
	swapParameterDirection();
	// This makes for a really easy routine. We now raise in v-direction.
	SplineCurve huge_curve(numCoefs_v(), order_v(),
			       basis_v().begin(), coef_iter,
			       rdim*numCoefs_u(), false);

	huge_curve.raiseOrder(raise_u);
	basis_v_ = huge_curve.basis();
	coefs_.empty();
	if (rat) {
	    rcoefs_.assign(huge_curve.coefs_begin(), huge_curve.coefs_end());
	} else {
	    coefs_.assign(huge_curve.coefs_begin(), huge_curve.coefs_end());
	}
	swapParameterDirection();
    }

    // We then raise in the v-direction.
    coef_iter = (rat) ? rcoefs_.begin() : coefs_.begin();
    if (raise_v > 0) {
	SplineCurve huge_curve(numCoefs_v(), order_v(),
			       basis_v().begin(), coef_iter,
			       rdim*numCoefs_u(), false);

	huge_curve.raiseOrder(raise_v);
	basis_v_ = huge_curve.basis();
	coefs_.empty();
	if (rat) {
	    rcoefs_.assign(huge_curve.coefs_begin(), huge_curve.coefs_end());
	} else {
	    coefs_.assign(huge_curve.coefs_begin(), huge_curve.coefs_end());
	}
    }

    if (rat) {
	updateCoefsFromRcoefs();
    }
}
