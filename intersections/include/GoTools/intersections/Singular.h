//===========================================================================
//                                                                           
// File: Singular.h
//                                                                           
// Created: April 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision:
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SINGULAR_H
#define _SINGULAR_H

#include "GoTools/utils/Point.h"
#include <vector>
#include <memory>

namespace Go {

/// Compute vanishing curve tangents and surface normals

    class ParamSurface;
    class ParamCurve;

namespace Singular
{

    void vanishingNormal(std::shared_ptr<ParamSurface> srf, double tol,
			 std::vector<Point>& singular_pts,  // Singular points in the parameter domain 
			 std::vector<std::vector<Point> >& singular_sequences);  // Sequences of parameter points
                                                                        	// making a singular curve

    void vanishingTangent(std::shared_ptr<ParamCurve> crv, 
			  double start, double end, double tol,
			 std::vector<double>& singular_pts,  // Singular points in the parameter domain 
			 std::vector<std::vector<double> >& singular_sequences);  // Sequences of parameter points
                                                                        	// making a singular curve

} // of namespace Singular

} // end namespace Go



#endif // _SINGULAR_H
