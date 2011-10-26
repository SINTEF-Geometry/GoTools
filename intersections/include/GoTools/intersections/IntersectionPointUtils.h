//==========================================================================
//                                                                          
// File: IntersectionPointUtils.h                                            
//                                                                          
// Created: Wed Jun 28 16:23:40 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: IntersectionPointUtils.h,v 1.3 2006-09-05 09:09:56 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

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

