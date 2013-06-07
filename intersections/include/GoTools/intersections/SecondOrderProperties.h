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

#ifndef _SECONDORDERPROPERTIES_H
#define _SECONDORDERPROPERTIES_H

#include <vector>
#include "GoTools/utils/Point.h"
#include "GoTools/intersections/SingularityType.h"

namespace Go {

// this structure encapsulates all information about a IntersectionPoint that 
// is or might be influenced by second order properties of the objects it is
// lying on. The need for this class arises because a geometrical object may
// carry a discontinuity in its second derivative at an IntersectionPoint, and
// its differential properties may therefore vary according to from where one
// look at the point.
// An IntersectionPoint lying on such a discontinuity (or more!), will therefore
// contain multiple SecondOrderProperties, all of which are correct depending
// on which point of view one takes.

struct SecondOrderProperties
{
    SecondOrderProperties() : 
	singularity_info_is_cached(false), tangent_is_oriented(false), tangent_2d_is_cached(false){}

    void clear() 
    {
	singularity_info_is_cached = false;
	tangent_is_oriented = false;
	tangent_2d_is_cached = false;
    }

    // caching of singularity type.  If the intersecting objects have 3D tangents (ie. they
    // are not parametrical functions), the 3D tangents will be cached simultaneously with
    // the SingularitytType.
    bool singularity_info_is_cached;
    Point tangent_3d[2]; // space for two tangents, in case of branch points
    SingularityType singularity_type;

    bool tangent_is_oriented;

    // caching of parametric tangent and related info
    bool tangent_2d_is_cached;
    Point tangent_2d_1[2]; // space for two tangents in case of branching point
    Point tangent_2d_2[2]; // space for two tangents in case of branching point

};

}; // end namespace Go

#endif // _SECONDORDERPROPERTIES_H

