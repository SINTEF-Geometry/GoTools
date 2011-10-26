//===========================================================================
//
// File : GSVbernsteinKnots.C
//
// Created: Thu Oct 30 13:12:07 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: GSVbernsteinKnots.C,v 1.1 2008-10-31 13:48:15 vsk Exp $
//
// Description:
//
//===========================================================================


#include "GoTools/trivariate/SplineVolume.h"

using namespace std;

namespace Go {


//==========================================================================
void SplineVolume::makeBernsteinKnots(int pardir)
//==========================================================================
{
    // @@ WARNING: Comparing floating point numbers for equality.

    vector<double> new_knots;
    SplineVolume the_volume = *this;
    int k = order(pardir); // Order of spline volume in specified direction

    // Loop through knot vector
    vector<double>::const_iterator iter = the_volume.basis(pardir).begin();
    while (iter != the_volume.basis(pardir).end()) {
	// Finds how many equal knots there are
	int counter = 1;
	while (iter+counter != the_volume.basis(pardir).end()
	       && *iter == *(iter+counter))
	    ++counter;

	ALWAYS_ERROR_IF(counter > k, "More equal knots than order.");

	// Fill up knots
	new_knots.insert(new_knots.end(), k - counter, *iter);

	iter += counter;
    }

    the_volume.insertKnot(pardir,new_knots);
    *this = the_volume;
    return;
}

//==========================================================================
int SplineVolume::numberOfPatches(int pardir) const
//==========================================================================
{
    // @@ WARNING: Comparing floating point numbers for equality.

    // Loop through knot vector
    vector<double>::const_iterator iter = basis(pardir).begin()+order(pardir)-1;
    int num = 0;
    while (iter != basis(pardir).begin()+numCoefs(pardir)) {
	int counter = 1;
	// Finds how many equal knots there are
	while (*iter == *(iter+counter)
	       && iter+counter != basis(pardir).begin()+numCoefs(pardir))
	    ++counter;
	++num;
	iter += counter;
    }

    return num;
}



} // namespace Go;
