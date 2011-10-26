//===========================================================================
//                                                                           
// File: SecondOrderProperties.h                                             
//                                                                           
// Created: Thu Jan  6 13:19:17 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: SecondOrderProperties.h,v 1.2 2006-08-22 09:22:42 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

