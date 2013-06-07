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

using namespace std;

namespace Go {


//==========================================================================
void SplineSurface::makeBernsteinKnotsU()
//==========================================================================
{
    // @@ WARNING: Comparing floating point numbers for equality.

    vector<double> new_knots;
    SplineSurface the_surface = *this;
    int k = order_u(); // Order of spline surface in u-direction

    // Loop through u-knot vector
    vector<double>::const_iterator iter = the_surface.basis_u().begin();
    while (iter != the_surface.basis_u().end()) {
	// Finds how many equal knots there are
	int counter = 1;
	while (iter+counter != the_surface.basis_u().end()
	       && *iter == *(iter+counter))
	    ++counter;

	ALWAYS_ERROR_IF(counter > k, "More equal knots than order.");

	// Fill up knots
	new_knots.insert(new_knots.end(), k - counter, *iter);

	iter += counter;
    }

    the_surface.insertKnot_u(new_knots);
    *this = the_surface;
    return;
}

//==========================================================================
void SplineSurface::makeBernsteinKnotsV()
//==========================================================================
{
    // @@ WARNING: Comparing floating point numbers for equality.

    vector<double> new_knots;
    SplineSurface the_surface = *this;
    int k = order_v(); // Order of spline surface in v-direction

    // Loop through v-knot vector
    vector<double>::const_iterator iter = the_surface.basis_v().begin();
    while (iter != the_surface.basis_v().end()) {
	// Finds how many equal knots there are
	int counter = 1;
	while (iter+counter != the_surface.basis_v().end()
	       && *iter == *(iter + counter))
	    ++counter;

	ALWAYS_ERROR_IF(counter > k, "More equal knots than order.");


	// Fill up knots
	new_knots.insert(new_knots.end(), k - counter, *iter);

	iter += counter;
    }

    the_surface.insertKnot_v(new_knots);
    *this = the_surface;
    return;
}


//==========================================================================
int SplineSurface::numberOfPatches_u() const
//==========================================================================
{
    // @@ WARNING: Comparing floating point numbers for equality.

    // Loop through u-knot vector
    vector<double>::const_iterator iter = basis_u().begin()+order_u()-1;
    int num = 0;
    while (iter != basis_u().begin()+numCoefs_u()) {
	int counter = 1;
	// Finds how many equal knots there are
	while (*iter == *(iter+counter)
	       && iter+counter != basis_u().begin()+numCoefs_u())
	    ++counter;
	++num;
	iter += counter;
    }

    return num;
}


//==========================================================================
int SplineSurface::numberOfPatches_v() const
//==========================================================================
{
    // @@ WARNING: Comparing floating point numbers for equality.

    // Loop through u-knot vector
    vector<double>::const_iterator iter = basis_v().begin()+order_v()-1;
    int num = 0;
    while (iter != basis_v().begin()+numCoefs_v()) {
	int counter = 1;
	// Finds how many equal knots there are
	while (*iter == *(iter+counter)
	       && iter+counter != basis_v().begin()+numCoefs_v())
	    ++counter;
	++num;
	iter += counter;
    }

    return num;
}


} // namespace Go;
