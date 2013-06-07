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
#include <list>
#include <memory>

using std::list;
using std::ostream;


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


