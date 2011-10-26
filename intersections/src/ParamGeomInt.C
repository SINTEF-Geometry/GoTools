//===========================================================================
//                                                                           
// File: ParamGeomInt.C 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: ParamGeomInt.C,v 1.5 2007-01-15 10:12:39 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/ParamGeomInt.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/CompositeBox.h"


namespace Go {


//===========================================================================
DirectionCone ParamGeomInt::reducedDirectionCone(bool reduce_at_bd[4],
						 double epsge) const
//===========================================================================
{
    // Default implementation
    return directionCone();
}


//===========================================================================
bool ParamGeomInt::isLinear(double epsge)
//===========================================================================
{
    // First make angular condition based on geometry size and tolerance
    CompositeBox box = compositeBox();
    Point low = box.low();
    Point high = box.high();

    int dim = low.dimension();
    double bsize = 0.0;
    for (int kj=0; kj<dim; kj++)
	bsize = std::max(bsize, high[kj]-low[kj]);

    double ang_tol = epsge/(2.0*bsize);
    DirectionCone cone;
    try {
	cone = directionCone();
    }
    catch (...)
    {
	return false;
    }

    if (cone.greaterThanPi() || cone.angle() >= ang_tol)
	return false;  // The current object is not linear
    else
	return true;
    
}

//===========================================================================
bool ParamGeomInt::coneLargerThanPi()
//===========================================================================
{
    DirectionCone cone;
    try {
	cone = directionCone();
    }
    catch (...)
    {
	return true;
    }

    if (cone.greaterThanPi())
	return true;

    return false;
}

//===========================================================================


} // namespace Go

