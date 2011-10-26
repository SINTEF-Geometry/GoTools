//===========================================================================
//
// File : CurveInterpolator
//
// Created: May 2011
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description: Interpolation of regular data points
//
//===========================================================================


#ifndef __CURVEINTERPOLATOR_H
#define __CURVEINTERPOLATOR_H

#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/utils/Point.h"

namespace Go
{
  class SplineCurve;

  namespace CurveInterpolator
  {
    /// Interpolate a set of regular, parameterized interpolation
    /// points. The parameterizetion is assumed to correspond to the
    /// given B-spline basis (select interpolation points in the
    /// Greville points). The function throws if the input parameters
    /// are inconsistent
    SplineCurve* regularInterpolation(const BsplineBasis& basis,
					std::vector<double>& par,
					std::vector<double>& points,
					int dimension,
					bool rational,
					std::vector<double>& weights);

  };    // namespace CurveInterpolator


} // namespace Go


#endif    // #ifndef __CURVEINTPEROLATOR_H

