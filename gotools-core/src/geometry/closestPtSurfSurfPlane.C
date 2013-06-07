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

#include "GoTools/geometry/ClosestPoint.h"

namespace Go {

//===========================================================================
void  ClosestPoint::closestPtSurfSurfPlane(const std::vector<Point>& epoint,
			    const std::vector<Point>& epnt1,
			    const std::vector<Point>& epnt2,
			    const Point& epar1,
			    const Point& epar2,
			    const ParamSurface* psurf1,
			    const ParamSurface* psurf2,
			    double aepsge,
			    std::vector<Point>& gpnt1,
			    std::vector<Point>& gpnt2,
			    Point& gpar1, 
			    Point& gpar2, 
			    int& jstat,
			    AlgorithmChoice algo)
//===========================================================================
{
    if (algo == GEOMETRICAL) {
	ClosestPoint::closestPtSurfSurfPlaneGeometrical(epoint, 
					  epnt1, 
					  epnt2, 
					  epar1, 
					  epar2, 
					  psurf1, 
					  psurf2, 
					  aepsge, 
					  gpnt1, 
					  gpnt2, 
					  gpar1, 
					  gpar2, 
					  jstat);
    } else {
	// algo == FUNCTIONAL
	ClosestPoint::closestPtSurfSurfPlaneFunctional(epoint, 
					  epnt1, 
					  epnt2, 
					  epar1, 
					  epar2, 
					  psurf1, 
					  psurf2, 
					  aepsge, 
					  gpnt1, 
					  gpnt2, 
					  gpar1, 
					  gpar2, 
					  jstat);
    }
}

}; // end namespace Go;

