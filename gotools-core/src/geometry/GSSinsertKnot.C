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
