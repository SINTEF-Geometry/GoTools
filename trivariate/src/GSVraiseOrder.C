//===========================================================================
//
// File : GSVraiseOrder.C
//
// Created: Mon Nov 24 08:46:08 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: GSVraiseOrder.C,v 1.1 2008-11-24 10:36:11 kfp Exp $
//
// Description: Raise order of spline volume. Method: treat volume as a curve.
//
//===========================================================================


#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"


using namespace Go;


//===========================================================================
void SplineVolume::raiseOrder(int raise_u, int raise_v, int raise_w)
//===========================================================================
{
    ALWAYS_ERROR_IF(raise_u < 0 || raise_v < 0 || raise_w < 0,
		    "Order to raise by must be non-negative!");

    // We're raising in the u-direction first.
    if (raise_u > 0)
      {
	swapParameterDirection(0, 2);
	raiseOrder_wdir(raise_u);
	swapParameterDirection(0, 2);
      }

    // We then raise in the v-direction.
    if (raise_v > 0)
      {
	swapParameterDirection(1, 2);
	raiseOrder_wdir(raise_v);
	swapParameterDirection(1, 2);
      }

    // Finally we raise in the w-direction.
    if (raise_w > 0)
      raiseOrder_wdir(raise_w);

    // If neccessary, we rebuild the coefs from the rational coefs
    if (rational_ && (raise_u > 0 || raise_v > 0 || raise_w > 0))
      updateCoefsFromRcoefs();

}





//===========================================================================
void SplineVolume::raiseOrder_wdir(int raise)
//===========================================================================
{
  std::vector<double>::iterator coef_iter = (rational_) ? rcoefs_.begin() : coefs_.begin();
  int rdim = (rational_) ? dimension() + 1 : dimension();

  SplineCurve huge_curve(numCoefs(2), order(2),
			 basis(2).begin(), coef_iter,
			 rdim*numCoefs(0)*numCoefs(1), false);

  huge_curve.raiseOrder(raise);
  basis_w_ = huge_curve.basis();
  coefs_.empty();
  if (rational_) {
    rcoefs_.assign(huge_curve.coefs_begin(), huge_curve.coefs_end());
  } else {
    coefs_.assign(huge_curve.coefs_begin(), huge_curve.coefs_end());
  }

}
