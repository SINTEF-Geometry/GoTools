//===========================================================================
//                                                                           
// File: GSCbernsteinKnots.C                                                 
//                                                                           
// Created: Wed May  2 16:17:54 2001                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: GSCbernsteinKnots.C,v 1.5 2005-06-09 07:15:52 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineCurve.h"


using namespace std;


namespace Go {


//==========================================================================
void SplineCurve::makeBernsteinKnots()
//==========================================================================
{
    // @@ WARNING: Comparing floating point numbers for equality.

    vector<double> new_knots;
    SplineCurve the_curve = *this;
    int k = order(); // Order of spline curve

    // Loop through knot vector
    vector<double>::const_iterator iter = the_curve.basis().begin();
    while (iter < the_curve.basis().end()) {
	// Finds how many equal knots there are
	int counter = 1;
	while (iter+counter < the_curve.basis().end()
	       && *iter == *(iter + counter))
	    ++counter;

	ALWAYS_ERROR_IF(counter > k, "More equal knots than order.");


	// Fill up knots
	new_knots.insert(new_knots.end(), k - counter, *iter);

	iter += counter;
    }

    the_curve.insertKnot(new_knots);
    *this = the_curve;
    return;
}


} // namespace Go;
