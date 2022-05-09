#ifndef _SSURFTRACEISOCONTOURS_H
#define _SSURFTRACEISOCONTOURS_H


#include <vector>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/TraceContoursTypedefs.h"

namespace Go
{
  /// Compute the level-set curves of the tensor product spline function 'ss'.
  // In other words, compute the intersections between a set of horizontal
  // planes and 'ss', interpreted as a 2 1/2 D surface. The height of the
  // planes (level-set values) are provided in the vector 'isovals'.  The output
  // is a vector with one entry per entry in 'isovals', containing the curves
  // that represents this level-set.  These are 2D curves in the parameter plane
  // of 'ss'.  If 'include_3D_curves' is set to 'true', the corresponding 3D
  // curves will also be returned.  These will have constant z-value (equal to
  // the corresponding entry in 'isovals'), but are useful for visualization
  // purposes in a 3D viewer.  By default, the function employs the function
  // 'traceIsovals' to march out the curves.  If desired, SISL routine s1314 can
  // be used instead by setting 'use_sisl_marching' to true.  The SISL routine
  // is slower.
  std::vector<CurveVec> SSurfTraceIsocontours(const SplineSurface& ss,
					      const std::vector<double>& isovals,
					      const double tol = 1e-6,
					      const bool include_3D_curves = false,
					      const bool use_sisl_marching = false);
}; // end namespace Go;

#endif
