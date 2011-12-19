//===========================================================================
//
// File : GSVderivVolume.C
//
// Created: Tue Nov 25 14:16:17 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: GSVderivVolume.C,v 1.2 2008-11-27 11:41:35 vsk Exp $
//
// Description:
//
//===========================================================================

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/trivariate/VolumeTools.h"
#include <memory>
#include <algorithm>


namespace Go {

//==========================================================================
SplineVolume* SplineVolume::derivVolume(int ider1, int ider2, int ider3) const
//==========================================================================
{
  // Currently works for non-rational volumes only

  ALWAYS_ERROR_IF(ider1 < 0 || ider2 < 0 || ider3 < 0,
		  "Trying to take a negative number of derivatives!");

  ALWAYS_ERROR_IF(rational_,
		  "Trying to get a derivative volume for rational volume. Not implemented yet");

  // Create curves representing the volume as a curve in the
  // arameter direction, make derivatives of these curves

  shared_ptr<SplineVolume> vol;
  shared_ptr<SplineCurve> temp_cv;
  shared_ptr<SplineCurve> temp_cv_deriv;

  // First in the u direction.
  if (ider1 > 0)
  {
      temp_cv = representVolumeAsCurve(*this, 0);
      temp_cv_deriv.reset(temp_cv->derivCurve(ider1));
      vol  = representCurveAsVolume(*temp_cv_deriv, 0, basis_v_, basis_w_, 0);
  }
  else
      vol = shared_ptr<SplineVolume>(clone());

  BsplineBasis new_ubasis = vol->basis(0);

  // Then in the v direction.
  if (ider2 > 0)
  {
      temp_cv = representVolumeAsCurve(*vol, 1);
      temp_cv_deriv.reset(temp_cv->derivCurve(ider2));
      vol = representCurveAsVolume(*temp_cv_deriv, 1, new_ubasis, basis_w_, 0);
      // vol = representCurveAsVolume(*temp_cv, 1, new_ubasis, basis_w_, 0);
  }
  BsplineBasis new_vbasis = vol->basis(1);

  // Finally in the w direction.
  if (ider3 > 0)
  {
      temp_cv = representVolumeAsCurve(*vol, 2);
      temp_cv_deriv.reset(temp_cv->derivCurve(ider3));
      vol = representCurveAsVolume(*temp_cv_deriv, 2, new_ubasis, new_vbasis, 0);
  }

  // We have to clone the return value because we cannot return a
  // shared pointer from this function.
  return vol->clone();

}


} // namespace Go;
