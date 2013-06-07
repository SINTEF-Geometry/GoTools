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

#ifndef _INTERSECTIONPOINTUTILS_H
#define _INTERSECTIONPOINTUTILS_H


#include "GoTools/intersections/SingularityType.h"


namespace Go {


/// Enumeration used to classify the direction of an intersection curve
/// passing through an intersection point with respect to surface
/// boundaries

enum IntPtClassification {
    /// Not defined
    DIR_UNDEF = 0,
    /// Tangent of the intersection curve pointing into the domain
    DIR_IN,
    /// Tangent of the intersection curve pointing out from the comain
    DIR_OUT,
    /// Tangent of the intersection curve parallel to the domain
    DIR_PARALLEL,
    /// Tangent of the intersection curve perpendicular to the domain
    DIR_PERPENDICULAR,
    /// No defined tangent
    DIR_HIGHLY_SINGULAR,
    /// Tangent touching the domain at a corner
    DIR_TOUCH
};


/// Enumeration used to classify the location of an intersection point
/// in the parameter domain of the object.
// @jbt: The values are modeled after the sisl routine
// sh6sislinside().

enum IntPtLocation {
    /// Point is outside domain for both objects
    LOC_OUTSIDE_BOTH = 0,
    /// Point is inside domains (inner point) for both objects
    LOC_INSIDE_BOTH,
    /// Point is on the edge of at least one object
    LOC_EDGE_ONE,
    /// Point is on the corner of at least one object
    LOC_CORNER_ONE,
    /// Point is on the corner of both objects
    LOC_CORNER_BOTH,
    /// Point is on the edge of both objects
    LOC_EDGE_BOTH
};


/// Helper struct for saving bracketed bounds of influence areas

struct CachedInterval {

    CachedInterval() : inside(0), outside(0), cached(false) {}
	
    double inside;
    double outside;
    bool cached;
};


/// Struct that holds diagnostic information about an intersection
/// point. Used by IntersectionPool::checkIntersectionPoints().

struct IntPtInfo {

    IntPtInfo() : is_ok(false) {}

    int nneighbours;
    SingularityType singularity_type;
    IntPtClassification direction;
    IntPtLocation location;
    bool is_ok;
    
};


} // namespace Go


#endif // _INTERSECTIONPOINTUTILS_H

