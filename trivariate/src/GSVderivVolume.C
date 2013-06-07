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
      temp_cv = VolumeTools::representVolumeAsCurve(*this, 0);
      temp_cv_deriv.reset(temp_cv->derivCurve(ider1));
      vol  = VolumeTools::representCurveAsVolume(*temp_cv_deriv, 0, basis_v_, basis_w_, 0);
  }
  else
      vol = shared_ptr<SplineVolume>(clone());

  BsplineBasis new_ubasis = vol->basis(0);

  // Then in the v direction.
  if (ider2 > 0)
  {
      temp_cv = VolumeTools::representVolumeAsCurve(*vol, 1);
      temp_cv_deriv.reset(temp_cv->derivCurve(ider2));
      vol = VolumeTools::representCurveAsVolume(*temp_cv_deriv, 1, new_ubasis, basis_w_, 0);
      // vol = representCurveAsVolume(*temp_cv, 1, new_ubasis, basis_w_, 0);
  }
  BsplineBasis new_vbasis = vol->basis(1);

  // Finally in the w direction.
  if (ider3 > 0)
  {
      temp_cv = VolumeTools::representVolumeAsCurve(*vol, 2);
      temp_cv_deriv.reset(temp_cv->derivCurve(ider3));
      vol = VolumeTools::representCurveAsVolume(*temp_cv_deriv, 2, new_ubasis, new_vbasis, 0);
  }

  // We have to clone the return value because we cannot return a
  // shared pointer from this function.
  return vol->clone();

}


} // namespace Go;
