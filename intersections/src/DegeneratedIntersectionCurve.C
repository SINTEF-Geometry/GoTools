//===========================================================================
//                                                                           
// File: DegeneratedIntersectionCurve.C                                      
//                                                                           
// Created: Wed Apr 20 17:03:38 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision:
// $Id: DegeneratedIntersectionCurve.C,v 1.3 2006-06-20 10:32:23 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include <stdexcept>


using std::logic_error;
using std::shared_ptr;


namespace Go {


//===========================================================================
DegeneratedIntersectionCurve::~DegeneratedIntersectionCurve() {}
//===========================================================================


//===========================================================================
shared_ptr<ParamCurve> DegeneratedIntersectionCurve::getCurve() const 
//===========================================================================
{
    // Curve is assumed to be collapsed into a single point.  Creating 
    // degenerated spline curve.

    const int dim = 3; //@ can we always assume this?
    const double knot[] = {0, 0, 1, 1};
    double coef[dim * 2];
    const Point temp = ipoints_.front()->getPoint();
    ASSERT(temp.dimension() == dim);
    for (int i = 0; i < dim; ++i) {
	coef[i] = coef[i + dim] = temp[i];
    }

    MESSAGE("Converting a degenerated IntersectionCurve into a ParamCurve!");
    return shared_ptr<ParamCurve>(new SplineCurve(2, 2, knot, coef,
						  dim, false));
}


//===========================================================================
shared_ptr<ParamCurve>
DegeneratedIntersectionCurve::getParamCurve(int obj_nmb) const 
//===========================================================================
{
    shared_ptr<ParamCurve> result;
    switch (obj_nmb) {
    case 1:
	// implement this
	// result = something
	break;
    case 2:
	// implement this
	//result = something
	break;
    default:
	throw logic_error("Argument to getParamCurve() should be 1 or 2.");
    }
    if (result.get() == 0) {
	MESSAGE("Warning;  Returned isocurve is a zero pointer.\n"
		"It should have been precalculated, but this functionality\n"
		"is not yet implemented. ");
    }

    MESSAGE("Converting a degenerated IntersectionCurve into a ParamCurve!");
    return result;
}


//===========================================================================


}; // end namespace Go
