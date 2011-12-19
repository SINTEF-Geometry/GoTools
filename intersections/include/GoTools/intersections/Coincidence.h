//===========================================================================
//                                                                           
// File: coincidence.h
//                                                                           
// Created: September 2004
//                                                                           
// Author: Vibeke Skytt
//          
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================


#ifndef _COINCIDENCE_H
#define _COINCIDENCE_H

#include <memory>
#include "GoTools/utils/Point.h"
#include "GoTools/intersections/ParamObjectInt.h"
#include "GoTools/intersections/GeoTol.h"
#include "GoTools/intersections/ParamFunctionInt.h"

#include <vector>

// Coincidence testing

namespace Go {

  class ParamCurveInt;
  class ParamSurfaceInt;
  class SplineSurfaceInt;
    class Param1FunctionInt;
    class Param0FunctionInt;

    /// Check if two curves are coinciding between parameters start and
    ///  end on this curve, and other_start and other_end on the other curve.
    int checkCoincide(ParamCurveInt *curve,
		      double start, double end, shared_ptr<GeoTol> tol,
		      ParamCurveInt *other,
		      double other_start, double other_end);

    /// Check if this curve is coinciding between parameters start and
    ///  end, with the ParamSurface between parameters su_start and su_end.
    int checkCoincide(ParamCurveInt *curve, double start, double end,
		      ParamSurfaceInt *surf,
		      Point su_start,
		      Point su_end,
		      shared_ptr<GeoTol> tol);

    /// Check if this curve is coinciding between parameters start and
    ///  end, with the SplineSurface between parameters su_start and su_end.
    int checkCoincide(ParamCurveInt *curve, double start, double end,
		      SplineSurfaceInt *surf,
		      const Point& su_start,
		      const Point& su_end,
		      shared_ptr<GeoTol> tol);

    /// Check if this function is coinciding between parameters start and
    /// end, with the Point.
    int checkCoincide(Param1FunctionInt *func1, double start, double end,
		      Param0FunctionInt *C,
		      shared_ptr<GeoTol> tol);

     /// Check if the two given surfaces coincide within the loop given
    /// by the parameter values of the intersection points representing it
    int checkCoincide(ParamSurfaceInt *surf1, ParamSurfaceInt *surf2,
		      std::vector<double>& par_loop, 
		      const shared_ptr<GeoTol> tol);

    // Computes the next step along a curve when we are testing if two
    // curves coincide. 
    int stepLength(const std::vector<Point>& ft, const std::vector<Point>& gs, 
		   double delta, bool forward, double& delta_t);


     // Computes the next step along a curve when we are testing if a
    // curve and a surface coincide. 
    double stepLength(const std::vector<Point>& ft, const std::vector<Point>& gs, 
		      bool forward, std::vector<Point>& cvder,
		      double min_step, double max_step, double aepsge); 

 
    // Evaluate the curve representing the projection of a spatial curve onto a surface
    int evalProjectedCurve(const std::vector<Point>& ft, 
			   const std::vector<Point>& gs, 
			   std::vector<Point>& res);

    void evalDistCurve(const std::vector<Point>& marching_curve, 
		       const std::vector<Point>& other_curve, 
		       std::vector<Point>& dist_curve);

    // Start in an intersection point and march along one parameter direction
    // as long as there is "coincidence" (within a given tolerance).  Then return the
    // parameter value for the last point found that was closer to the other object
    // than the tolerance, and the parameter value for the first point found that is
    // more distant than the tolerance.  
    // 'obj1'                          - the first object, could be a curve or surface
    // 'obj2'                          - the second object, could be a curve or surface 
    // 'geotol'                        - the geometric tolerance for determining 
    //                                   coincidence.
    // 'current_params'                - a pointer to the parameter array that specify
    //                                   the position of the intersection point on the
    //                                   two objects.
    // 'dir'                           - specify along which parameter direction to march
    // 'forward'                       - specify whether we will be marching forwards or
    //                                   backwards.
    // 'last_param_val_inside'         - returns the last value found along parameter
    //                                   'dir' that corresponded to a position still
    //                                   within the given tolerance 'geotol'.
    // 'first_param_val_outside'       - returns the last value found along parameter
    //                                   'dir' that corresponded to a position outside
    //                                   the given tolerance 'geotol'.
    // returns: 'true' if it managed to bracket the interval properly, 'false'
    //          if it reached the end of the marching curve before the outside
    //          bracket was found.
    bool determineCoincidenceRegion(const ParamObjectInt* obj1,
				    const ParamObjectInt* obj2,
				    shared_ptr<const GeoTol> tol,
				    const double* current_params,
				    int dir,
				    bool forward,
				    double& last_param_val_inside,
				    double& first_param_val_outside);

    // Special version for 0-dim 2nd object (1st object is 1- or 2-par of dim 1).
    bool determineCoincidenceRegion(const ParamFunctionInt* obj_1d,
				    double C,
				    shared_ptr<const GeoTol> tol,
				    const double* current_params,
				    int dir,
				    bool forward,
				    double& last_param_val_inside,
				    double& first_param_val_outside);

//     void determine_confidence_interval(const ParamCurveInt* const marching_curve,
// 				       const ParamCurveInt* const other_curve,
// 				       double& inside_param_marching,
// 				       double& inside_param_other,
// 				       double& outside_param_marching,
// 				       double& outside_param_other,
// 				       double& inside_val,
// 				       double& outside_val,
// 				       double center_value,
// 				       double bracket,
// 				       bool& succeeded);

//     // Start in an intersection point and march along two curves in a given
//     // direction as long as they coincide.  Then return the two parameter values
//     // for the last points found that are closer to each other than the tolerance, 
//     // and the two parameter values for the first points found that are more 
//     // distant.
//     // 'marching_curve'                - the first curve, along which we will march
//     // 'other_curve'                   - the second curve
//     // 'isect_param_marching'          - the marching curve's parameter for the
//     //                                   departing intersection point
//     // 'isect_param_other'             - the second curve's parameter for the
//     //                                   departing intersection point
//     // 'step_forward_marching'         - set to 'false' if we want to march in the
//     //                                   inverse parameter direction for the 
//     //                                   marching curve.
//     // 'aepsge'                        - contains the geometric tolerance (aepsge) 
//     //                                   and the bracket interval size (bracket_size)
//     // 'last_param_inside_marching'    - return the last parameter found on the
//     //                                   marching curve that corresponds with a
//     //                                   geometric point on the other curve with
//     //                                   a distance inferior to the geometric 
//     //                                   tolerance. It should be inside 
//     //                                   [(1-bracket_size) * aepsge, aepsge]
//     // 'first_param_outside_marching'  - return the first parameter found on the
//     //                                   marching curve that corresponds with a 
//     //                                   geometric point on the other curve with
//     //                                   a distance superior to the geometric 
//     //                                   tolerance, but still within the interval
//     //                                   [aepsge, (1+bracket_size) * aepsge]
//     // 'last_param_inside_other'       - return the parameter on the second curve
//     //                                   that represents the closest point found
//     //                                   to the point on the marching curve 
//     //                                   represented by 'last_param_inside_marching'
//     // 'first_param_outside_other'     - return the parameter on the second curve
//     //                                   that represents the closest point found to 
//     //                                   the point on the marching curve represented
//     //                                   by 'first_param_outside_marching'.
//     // returns: 'true' if it managed to bracket the interval properly, 'false'
//     //          if it reached the end of the marching curve before the outside
//     //          bracket was found.
//     bool measureCoincidenceRegion(const ParamCurveInt* marching_curve,
// 				  const ParamCurveInt* other_curve,
// 				  double isect_param_marching,
// 				  double isect_param_other,
// 				  bool step_forward_marching,
// 				  const shared_ptr<const GeoTol> aepsge,
// 				  double& last_param_inside_marching,
// 				  double& first_param_outside_marching,
// 				  double& last_param_inside_other,
// 				  double& first_param_outside_other);

//     // Helper functions for 'measureCoincidenceRegion(...)'
//     void expand_coincidence_region(const ParamCurveInt* const marching_curve,
// 				   const ParamCurveInt* const secondry_curve,
// 				   bool marching_fwd,
// 				   double& marching_par,
// 				   double& secondry_par,
// 				   double& dist,
// 				   double lower_limit, // smaller than eps
// 				   double eps,
// 				   bool& passed_lower_limit);


 
} // namespace Go

#endif // _COINCIDENCE_H
