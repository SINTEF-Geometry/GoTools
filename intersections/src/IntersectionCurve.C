//===========================================================================
//                                                                           
// File: IntersectionCurve.C                                                 
//                                                                           
// Created: Wed Apr 20 16:51:58 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: IntersectionCurve.C,v 1.30 2006-06-20 10:49:32 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/intersections/IntersectionCurve.h"
#include <list>
#include <memory>

using std::list;
using std::ostream;
using std::shared_ptr;


namespace Go {


//===========================================================================
shared_ptr<IntersectionPoint>
IntersectionCurve::getGuidePoint(int index) const 
//===========================================================================
{
    DEBUG_ERROR_IF(index < 0 || index >= int(ipoints_.size()),
	     "index out of range in IntersectionCurve::getGuidePoint()");

    list<shared_ptr<IntersectionPoint> >::const_iterator it
	= ipoints_.begin();
    for (int i = 0; i < index; ++i) {
	++it;
    }
    return *it;
} 


//===========================================================================
bool IntersectionCurve::
getGuidePointTangent(shared_ptr<IntersectionPoint> pt, 
		     Point& tan, 
		     int type) const
//===========================================================================
{
    list<shared_ptr<IntersectionPoint> >::const_iterator it;
    it = find(ipoints_.begin(), ipoints_.end(), pt);
    if (it != ipoints_.end()) {
	switch(type) {
	case 1:
	    tan = (*it)->getPar1Dir(); // bool second_branch = false ??
	    break;
	case 2:
	    tan = (*it)->getPar2Dir(); // bool second_branch = false ??
	    break;
	default:
	    tan = (*it)->getTangent();
	}
	return true;
    } 
    // could not find this point
    return false;
}

//===========================================================================
void IntersectionCurve::writeIPointsToStream(ostream& os) const
//===========================================================================
{
    size_t num_pts = ipoints_.size();
    os << num_pts << '\n';

    for (list<shared_ptr<IntersectionPoint> >::const_iterator it
	     = ipoints_.begin(); it != ipoints_.end(); ++it) {
	(*it)->write(os);
	os << '\n';
    }
}


//===========================================================================


}; // end namespace Go


