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
