//===========================================================================
//
// File : GSVinsertKnot.C
//
// Created: Fri Nov 21 10:47:47 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: GSVinsertKnot.C,v 1.1 2008-11-21 11:33:29 kfp Exp $
//
// Description:
//
//===========================================================================


#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"

using namespace Go;

//===========================================================================
void SplineVolume::insertKnot(int pardir, double apar)
//===========================================================================
{
  int kdim = rational_ ? dim_+1 : dim_;

  if (pardir == 0)
    swapParameterDirection(0,2);
  else if (pardir == 1)
    swapParameterDirection(1,2);

  // Make a hypercurve from this volume
  SplineCurve cv(numCoefs(2), order(2), basis_w_.begin(),
		 activeCoefs().begin(), kdim*numCoefs(1)*numCoefs(0), false);

  // Insert the knot in the curve
  cv.insertKnot(apar);

  // Insert knot into basis
  basis_w_.insertKnot(apar);

  // Copy back the data
  if (!rational_)
    {
      coefs_.resize(numCoefs(0)*numCoefs(1)*numCoefs(2)*dim_);
      std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
    }
  else
    {
      rcoefs_.resize(numCoefs(0)*numCoefs(1)*numCoefs(2)*(dim_+1));
      std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
      updateCoefsFromRcoefs();
    }

  if (pardir == 0)
    swapParameterDirection(0,2);
  else if (pardir == 1)
    swapParameterDirection(1,2);

}



//===========================================================================
void SplineVolume::insertKnot(int pardir, const std::vector<double>& new_knots)
//===========================================================================
{
  int kdim = rational_ ? dim_+1 : dim_;

  if (pardir == 0)
    swapParameterDirection(0,2);
  else if (pardir == 1)
    swapParameterDirection(1,2);

  // Make a hypercurve from this volume
  SplineCurve cv(numCoefs(2), order(2), basis_w_.begin(),
		 activeCoefs().begin(), kdim*numCoefs(1)*numCoefs(0), false);

  // Insert the knots in the curve
  cv.insertKnot(new_knots);

  // Insert knots into basis
  basis_w_.insertKnot(new_knots);

  // Copy back the data
  if (!rational_)
    {
      coefs_.resize(numCoefs(0)*numCoefs(1)*numCoefs(2)*dim_);
      std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
    }
  else
    {
      rcoefs_.resize(numCoefs(0)*numCoefs(1)*numCoefs(2)*(dim_+1));
      std::copy(cv.coefs_begin(), cv.coefs_end(), activeCoefs().begin());
      updateCoefsFromRcoefs();
    }

  if (pardir == 0)
    swapParameterDirection(0,2);
  else if (pardir == 1)
    swapParameterDirection(1,2);

}

