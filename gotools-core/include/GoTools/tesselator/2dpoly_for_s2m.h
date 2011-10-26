#ifndef TDPOLY_FOR_S2M_H_INCLUDED
#define TDPOLY_FOR_S2M_H_INCLUDED






//#include <functional>
//#include <algorithm>
#include <vector>

// Hvilken av alle disse definerer navnerommet Go?!

//#include "DefaultDataHandler.h"

// #include "gvCurveTesselator.h"
// #include "gvRectangularSurfaceTesselator.h"
// #include "gvNoopTesselator.h"
// #include "gvLineCloudTesselator.h"
// #include "gvRectGridTesselator.h"

// #include "gvCurvePaintable.h"
// #include "gvRectangularSurfacePaintable.h"
// #include "gvPointCloudPaintable.h"
// #include "gvLineCloudPaintable.h"
// #include "gvQuadsPaintable.h"

// #include "gvParametricSurfacePaintable.h"
// #include "gvParametricSurfaceTesselator.h"
// #include "ParametricSurfacePropertySheet.h"

// #include "gvPropertySheet.h"
// #include "GoTools/geometry/SplineCurvePropertySheet.h"
// #include "RectangularSurfacePropertySheet.h"
// #include "PointCloudPropertySheet.h"

// #include "GoTools/geometry/SplineSurface.h"
// #include "BoundedSurface.h"
// #include "GoTools/geometry/SplineCurve.h"
// #include "PointCloud.h"
// #include "GoTools/geometry/LineCloud.h"
// #include "GoTools/geometry/ParamSurface.h"
// #include "RectGrid.h"
// #include "ClassType.h"

// #include "Factory.h"

#include "GoTools/utils/Array.h" // for Vector3D

#include "GoTools/utils/errormacros.h"






using namespace Go;

// using std::vector; // 100213: Not a good idea to use in a header file?




typedef std::pair<short, std::vector<short> *> short_list;
typedef std::pair<short, std::vector<short_list> *> short_list_short_list;


namespace Go
{
  
  
  
  bool point_inside_contour(const double x0, const double y0,
			    const double * const vertices,
			    const std::vector<int> &contour
			    
			    // 090129: Usage of this has not been implemented so far. I am not sure it is a good
			    //         idea to do it either. Assumptions upon which its usefulness was once
			    //         based may not be present.  Commenting it out in order to avoid computing
			    //         it...
			    //vector<short_list_short_list> &sorted_segments
    );
  
  // 090115:
  // bool point_on_contour_corner(const double x0, const double y0,
  // 			     const double * const vertices, const vector<int> &contour);
  
  bool segment_contour_intersection_for_s2m(const double x0, const double y0,
					    const double x1, const double y1,
					    const double * const vertices,
					    const std::vector<int> &contour,
					    
					    // 090129: Usage of this has not been implemented so far. I am not
					    //         sure it is a good idea to do it either. Assumptions upon
					    //         which its usefulness was once based may not be present.
					    //         Commenting it out in order to avoid computing it...
					    //const vector<short_list_short_list> &sorted_segments,
					    
					    double &x, double &y, double &s,
					    const bool snap_ends = false
#ifdef DBG
					    , const bool dbg = false
#endif
    );
  
  std::vector< short_list_short_list > sort_2dpoly_segments(const double * const vertices,
						       const std::vector<int> &contour,
						       const bool transposed=false);
  
  int is_inside(const std::vector< Vector3D > &trim_curve_p, const std::vector<int> &contour,
		const double u, const double v
#ifdef DBG
		, const bool dbg=false
#endif
    );
  
  
  // 090115: This must be (re)checked before being used...
  bool is_on_corner(const std::vector< Vector3D > &trim_curve_p, const std::vector<int> &contour,
		    const double u, const double v);
  
  // 090117:
  int is_on_contour(const std::vector< Vector3D > &trim_curve_p, const std::vector<int> &contour,
		    const double u, const double v
#ifdef DBG
		    , const bool dbg=false
#endif
    );
  
  
  // 090204:
  bool degenerate_triangle(const Vector2D &c1, const Vector2D &c2, const Vector2D &c3
#ifdef DBG
			   , const bool dbg = false
#endif
    );
  
  
  



} // namespace Go






#endif
