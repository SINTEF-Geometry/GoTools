/**
\example adapt_curve adapt_curve.C
\verbatim
\endverbatim
This program demonstrates the use of the class AdaptCurve.
The class can generate a B-spline curve that approximates
an evaluator based curve to satisfy a given accuracy.

The program reads a 2D parameter curve and a surface from file, and
lifts the curve onto the surface to make an evaluator based curve.
This evaluator based curve are then used as input to the AdaptCurve
constructor.

\example append_curve append_curve.C
\verbatim
\endverbatim
This program demonstrates the use of the function
SplineCurve::append_curve declared in SplineCurve.h.
The function joins two SplineCurves by appending the start of the second
curve to the end of the first curve. The two curves must be of the same type

\example approx_curve approx_curve.C
\verbatim
\endverbatim
This program demonstrates the use of the class ApproxCurve.
The class can generate a B-spline curve that approximates a set
of parametrized points for a given accuracy by inserting new knots in
the curve until the required accuracy is reached.

\example approx_surface approx_surface.C
\verbatim
\endverbatim
This program demonstrates the use of the class ApproxSurf.
The class can generate a new tensor product B-spline surface with four boundary
curves that approximates a set of parametrized points for a given accuracy, or
modify an old surface by a set of parametrized points.

\example circle circle.C
\verbatim
\endverbatim
This program demonstrates the use of the class Circle.
It is a subclass of ElementaryCurve.
The space dimension of a circle is either 2 or 3.
The default parametrization is an angle from 0 to 2*PI.

\example closestpoint_curve closestpoint_curve.C
\verbatim
\endverbatim
This program demonstrates the use of the two functions
ParamCurve::closestPoint(const Point& pt, double& clo_t, Point& clo_pt,
double& clo_dist)  and
SplineCurve::closestPoint(const Point& pt, double tmin, double tmax,
double& clo_t, Point& clo_pt, double& clo_dist, double const *seed = 0).
They compute the closest point on a curve from a specified point.

\example closestpoint_degenerate_sf closestpoint_degenerate_sf.C
\verbatim
\endverbatim
This program demonstrates the use of the function
SplineSurface::closestPoint(const Point& pt, double& clo_u, double& clo_v, 
Point& clo_pt, double& clo_dist, double epsilon,
const RectDomain* domain_of_interest = NULL, double *seed = 0)
The declaration of the function is in 'ParamSurface.h'.
The function compute the closest point on a surface from a specified point. 
It reads a spline surface file in Go-format and a file in plain ASCII format
with the xyz-coordinates of the points we want to find the closest point on
the surface to.

\example closestpoint_surface closestpoint_surface.C
\verbatim
\endverbatim
This program demonstrates the use of the function
'SplineSurface::closestPoint(const Point& pt, double& clo_u, double& clo_v, 
Point& clo_pt, double& clo_dist, double epsilon,
const RectDomain* domain_of_interest = NULL, double *seed = 0)'.
The declaration of the function is in 'ParamSurface.h'.
The function compute the closest point on a surface from a specified point. 
It reads a spline surface file in Go-format and a file in plain ASCII format
with the xyz-coordinates of the points we want to find the closest point on
the surface to.

\example cone cone.C
\verbatim
\endverbatim
This program demonstrates the use of the class Cone.
It is a subclass of ElementarySurface.
The space dimension of a cone is 3.
The default parametrization is the angle u from 0 to 2*PI and the
distance v from minus infinity to plus infinity.

\example const_param_curves const_param_curves.C
\verbatim
\endverbatim
The program demonstrates the use of the function
SplineSurface::constParamCurve(double parameter, bool pardir_is_u).
The declaration of the function is in SplineSurface.h.

\example coons_patch_gen coons_patch_gen.C
\verbatim
\endverbatim
This program demonstrates the use of some of the functions in namespace
CoonsPatchGen.
The functions can be used to create a Coons Patch or a Gordon Surface.
The functions returns a SplineSurface pointer to the created surface.

\example cylinder cylinder.C
\verbatim
\endverbatim
This program demonstrates the use of the class Cylinder.
It is a subclass of ElementarySurface.
The space dimension of a cylinder is 3.
The default parametrization is the angle u from 0 to 2*PI and the
distance v from minus infinity to plus infinity.

\example ellipse ellipse.C
\verbatim
\endverbatim
This program demonstrates the use of the class Ellipse.
It is a subclass of ElementaryCurve.
The space dimension of an ellipse is either 2 or 3.
The default parametrization is an angle from 0 to 2*PI.

\example interpol_curve_free interpol_curve_free.C
\verbatim
\endverbatim
This program reads a point data set from a file, interpolates a spline curve
through the points and writes an output file with the spline curve and the
input points.

\example interpol_curve_hermite interpol_curve_hermite.C
\verbatim
\endverbatim
This program reads a point data set from a file, interpolates a spline curve
through the points and write two output files: One file with spline curve
data and one file with tangent vectors from the input points.

\example linear_swept_surface linear_swept_surface.C
\verbatim
\endverbatim
This program demonstrates the use of the static function 'linearSweptSurface'
in the class 'SweepSurfaceCreator'.
The function can generate a B-spline surface by sweeping one curve along
another. A given point on the sweeping curve will be swept along the other
curve.

\example project_curve project_curve.C
\verbatim
\endverbatim
This program demonstrates the use of the function
void CurveCreators::projectCurve(shared_ptr<ParamCurve>& space_cv,
                                 shared_ptr<ParamSurface>& surf,
  			            double epsge,
			            shared_ptr<SplineCurve>& proj_cv,
			            shared_ptr<SplineCurve>& par_cv)

The function generates a cubic spline curve(order four) which lies on a given
SplineSurface and is the projection of a given space curve (in 3D space) onto
that surface, within a given tolerance. The given space curve should be close
to the surface.

\example rotational_swept_surface rotational_swept_surface.C
\verbatim
\endverbatim
This program demonstrates the use of the static function
\em rotationalSweptSurface in the class \em SweepSurfaceCreator.
The function can generate a B-spline surface by rotating a curve around 
an axis.

\example sphere sphere.C
\verbatim
\endverbatim
This program demonstrates the use of the class Sphere.
It is a subclass of ElementarySurface.
The space dimension of a sphere is 3.
The default parametrization is the angles u from 0 to 2*PI and v from
-PI/2 to +PI/2.

\example surface_of_revolution surface_of_revolution.C
\verbatim
\endverbatim
This program demonstrates the use of the class SurfaceOfRevolution.
SurfaceOfRevolution is swept out by a SplineCurve that is rotated
around an axis with a complete revolution, and is thereby a
parametric surface. The space dimension is 3.
The curve must be such that it doesn't lead to a self-intersecting surface.

\example torus torus.C
\verbatim
\endverbatim
This program demonstrates the use of the class Torus.
It is a subclass of ElementarySurface.
The space dimension of a torus is 3.
The default parametrization is the angles u along the major circle and v
along the minor circle from 0 to 2*PI.

*/
