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
