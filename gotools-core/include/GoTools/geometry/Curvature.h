//===========================================================================
//
// File : Curvature.h
//
// Created: Fri Aug 22 16:58:07 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: Curvature.h,v 1.1 2008-08-27 10:48:19 kfp Exp $
//
// Description:
//
//===========================================================================


#ifndef __CURVATURE_H
#define __CURVATURE_H


#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SplineCurve.h"
#include <vector>





namespace Go
{

  /// Curvature analysis related to curves

  /// Get the points on a SplineCurve with a given curvature radius
  /// Result is stored in a vector as parameter values
  void curvatureRadiusPoints(const SplineCurve& curve,
			     double curveRad,
			     std::vector<double>& pos);

  /// Get the minimal curvature radius, and the parameter value of the point
  /// with the minimal curvature radius
  void minimalCurvatureRadius(const SplineCurve& curve,
			      double& mincurv,
			      double& pos);

} // namespace Go



#endif    // #ifndef __CURVATURE_H

