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

#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include <stdexcept>


using std::logic_error;


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
