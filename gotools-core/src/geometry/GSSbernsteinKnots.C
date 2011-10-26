//===========================================================================
//                                                                           
// File: GSSbernsteinKnots.C                                                 
//                                                                           
// Created: Tue Jun 19 17:33:47 2001                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: GSSbernsteinKnots.C,v 1.8 2005-06-27 15:32:59 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
