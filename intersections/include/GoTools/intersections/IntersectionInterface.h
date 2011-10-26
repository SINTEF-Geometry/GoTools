//===========================================================================
//                                                                           
// File: IntersectionInterface.h
//                                                                           
// Created: Dec. 08
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
//===========================================================================

#ifndef _INTERSECTIONINTERFACE_H
#define _INTERSECTIONINTERFACE_H

#include "GoTools/geometry/ParamCurve.h"
#include <vector>

// This collection of functions provides an interface to the GoTools intersection
// functionality

namespace Go
{
    /// Intersection between two parametric curves
    void intersectCurves(std::shared_ptr<ParamCurve> crv1, std::shared_ptr<ParamCurve> crv2,
			 double tol, std::vector<std::pair<double, double> >& intersection_points);

} // namespace Go

#endif // _INTERSECTIONINTERFACE_H

