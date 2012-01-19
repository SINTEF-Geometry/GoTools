//===========================================================================
//                                                                           
// File: ModifySurf.h
//                                                                           
// Created: Mars 2010
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision:
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _MODIFYSURF_H
#define _MODIFYSURF_H

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

namespace Go {

  /// Functionality used to modify one or more spline surfaces with
  /// respect to given conditions
namespace ModifySurf
{
  /// Replace one boundary curve of a SplineSurface
  /// Approximate the initial surface
  /// bd_idx = 0: umin
  ///          1: umax
  ///          2: vmin
  // /         3: vmax
  void replaceBoundary(shared_ptr<SplineSurface> surf,
		       shared_ptr<SplineCurve> curve,
		       int bd_idx, double tol);

  /// Enfore colinearity between coeffients at the common boundary between
  /// two spline surfaces. The two outer rows of coefficients are involved
  /// in the linearity constraints, but a few additional rows are modified
  /// for reasons of smoothness. All surface boundaries, except the affected
  /// one, are kept fixed.
  bool enforceCoefCoLinearity(shared_ptr<SplineSurface> sf1, int bd1, 
			      shared_ptr<SplineSurface> sf2, int bd2, 
			      double tol, 
			      std::vector<std::vector<int> >& enumeration);

} // of namespace ModifySurf

}; // end namespace Go

#endif // _MODIFYSURF_H
