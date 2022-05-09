#ifndef _LRTRACEISOCONTOURS_H
#define _LRTRACEISOCONTOURS_H


#include <vector>
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/TraceContoursTypedefs.h"
#include "GoTools/geometry/SplineCurve.h"

namespace Go
{
  /// Compute the level-set curves of the LR spline function 'lrs'.
  /// The input surface is split into a number of tensor-product spline patches,
  /// using the function LRSplineSurface::subdivideIntoSimpler,
  /// level curves are computed for every patch and the resulting
  /// curves are merged across patch boundaries.
  // In other words, compute the intersections between a set of horizontal planes and
  // 'lrs', interpreted as a 2 1/2 D surface. The height of the planes
  // (level-set values) are provided in the vector 'isovals'.  The output is a
  // vector with one entry per entry in 'isovals', containing the curves that
  // represents this level-set.  These are 2D curves in the parameter plane of
  // 'lrs'.  If 'include_3D_curves' is set to 'true', the corresponding 3D
  // curves will also be returned.  These will have constant z-value (equal to
  // the corresponding entry in 'isovals'), but are useful for visualization
  // purposes in a 3D viewer.  By default, the function employs the function
  // 'traceIsovals' to march out the curves.  If desired, SISL routine s1314 can
  // be used instead by setting 'use_sisl_marching' to true.  The SISL routine
  // is slower.
  std::vector<CurveVec> LRTraceIsocontours(const shared_ptr<ParamSurface>& surf,
					   const std::vector<double>& isovals,
					   const int threshold_missing,
					   const double tol); 

  std::vector<CurveVec> LRTraceIsocontours(const LRSplineSurface& lrs,
					   const std::vector<double>& isovals,
					   const int threshold_missing,
					   const double tol = 1e-6, 
					   const bool include_3D_curves = false,
					   const bool use_sisl_marching = false,
					   const CurveBoundedDomain* domain = NULL);
}; // end namespace Go;

#endif
