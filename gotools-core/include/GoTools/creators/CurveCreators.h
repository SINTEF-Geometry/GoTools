//===========================================================================
//                                                                           
// File: CurveCreators.h                                                    
//                                                                           
// Created: Fri Jun 15 14:14:59 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: CurveCreators.h,v 1.5 2009-05-13 07:30:28 vsk Exp $
//                                                                           
// Description: Various SplineCurve creators.
//                                                                           
//===========================================================================

#ifndef _CURVECREATORS_H
#define _CURVECREATORS_H

#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/SplineSurface.h"


namespace Go
{


class SplineCurve;

/// Various functions for generating SplineCurve s by approximation, blending,
/// etc.
namespace CurveCreators
{

    /// Given a spline curve f_i on a space of dimension n, and a (1-dim)
    /// spline function alpha_i (not necessarily a Bezier curve) defined on
    /// the same space, return the product.

    /// Construct a SplineCurve which is the 'product' of another SplineCurve
    /// and a function (which is given as a 1-dimensional SplineCurve).  
    /// \b NB: The two SplineCurves given as argument must be defined on exactly
    /// the same parameter interval.  Their knotvectors and orders do not have
    /// to be equal, though.
    /// \param alpha the multiplier function, given here as a 1-dimensional
    ///              SplineCurve
    /// \param f the SplineCurve to be multiplied with
    /// \return A raw pointer to the newly generated SplineCurve, a product
    ///         between 'f' and the 'alpha' function.  The user assumes ownership.
    SplineCurve GO_API *multCurveWithFunction(const SplineCurve& alpha,
                                              const SplineCurve& f);

    /// Given blend functions alpha_1 & alpha_2 (represented as 1-dimensional
    /// SplineCurves), and the spline curves f_1 and f_2, return the blended
    /// expression 'alpha_1 * f_1 + alpha_2 * f_2'.
    /// \b NB: All the curves given as input must e defined on the same parameter
    /// interval.  Their knotvectors and orders do not have to be equal, though.
    /// \param alpha_1 the multiplier function for the first curve, 'f_1'.  
    ///                It is represented as an 1-dimensional SplineCurve.
    /// \param f_1 The first SplineCurve to be blended.
    /// \param alpha_2 the multiplier function for the second curve, 'f_2'.
    ///                It is represented as an 1-dimensional SplineCurve
    /// \param f_2 The second SplineCurve to be blended.
    /// \return A raw pointer to the newly generated SplineCurve, a blend
    ///         between 'f_1' and 'f_2'.  The user assumes ownership.
    SplineCurve GO_API *blend(const SplineCurve& alpha_1,
                              const SplineCurve& f_1,
                              const SplineCurve& alpha_2,
                              const SplineCurve& f_2);

    /// Given input curves, we approximate with one cubic spline curve,
    /// fulfilling end requirements.

    /// This function is used generate one SplineCurve which approximates
    /// a sequence of (supposedly end-to-end continuous) input SplineCurves.
    /// In addition, the user can define explicitly the start and end point and
    /// tangent of the generated curve.
    /// \param first_crv pointer to the start of a range of (shared pointers to)
    ///                  SplineCurves, supposedly stored so that the end point 
    ///                  of each curve coincides with the start point of the 
    ///                  succeeding one.
    /// \param last_crv  pointer to one-past-last of the range of SplineCurves.
    /// \param start_pt  The user can choose to explicitly specify the start point
    ///                  and tangent of the curve to be generated.  If so, then 
    ///                  this vector must contain one or two elements: the start 
    ///                  point and optionally the start tangent.
    /// \param end_pt    The user can choose to explicitly specify the end point
    ///                  and tangent of the curve to be generated.  If so, then
    ///                  this vector must contain one or two elements: the end 
    ///                  point and optionally the end tangent.
    /// \param approxtol The routine will try to generate a curve that fits the
    ///                  given curves within a tolerance of 'approxtol', but is 
    ///                  not guaranteed to succeed.  If unsuccessful, a message
    ///                  will be written to standard output.
    /// \param maxdist   Reports the maximum found distance between the generated
    ///                  curve and the input curves.
    /// \retval max_iter Specify the maximum number of iterations that is allowed
    ///                  in order to converge to a solution.
    /// \return a raw pointer to the generated SplineCurve.  User assumes ownership.
    SplineCurve GO_API *approxCurves(std::shared_ptr<SplineCurve>* first_crv,
                                     std::shared_ptr<SplineCurve>* last_crv,
                                     const std::vector<Point>& start_pt, 
                                     const std::vector<Point>& end_pt,
                                     double approxtol, 
                                     double& maxdist, 
                                     int max_iter = 5);

/*     SplineCurve GO_API * */
/*     projectSpaceCurve(std::shared_ptr<ParamCurve>& space_cv, */
/*                       std::shared_ptr<ParamSurface>& surf, */
/*                       std::shared_ptr<Point>& start_par_pt, */
/*                       std::shared_ptr<Point>& end_par_pt, */
/*                       double epsge, */
/*                       const RectDomain* domain_of_interest = NULL); */

    void
    projectCurve(std::shared_ptr<ParamCurve>& space_cv,
		 std::shared_ptr<ParamSurface>& surf,
		 double epsge,
		 std::shared_ptr<SplineCurve>& proj_cv,
		 std::shared_ptr<SplineCurve>& par_cv);

    /// Project the space_cv into parameter domain given by surf.
    /// start_par_pt (& end*_ not needed, but useful to ensure correct evaluation.

    /// Generate a SplineCurve which lies on a given (part of a) SplineSurface and 
    /// is the projection og a given SplineCurve onto that surface.
    /// \param space_cv the SplineCurve (in 3D space) that we want to project onto
    ///                 the surface.
    /// \param surf The surface we will project the curve onto.  The resulting curve
    ///             will therefore lie in this surface.
    /// \param start_par_pt if the user wants to explicitly define where the start
    ///                     point of the curve should lie, 'start_par_pt' could be
    ///                     set to point to this point.  It can also be a zero pointer,
    ///                     in which case the start point is determined just like any
    ///                     other point
    /// \param end_par_pt if the user wants to explicitly define where the end point
    ///                   of the curve should lie, 'end_par_pt' could be set to point
    ///                   to this point.  It can also be a zero pointer, in which case
    ///                   the end point is determined just like any other point.
    /// \param epsge The geometrical tolerance to use when generating the projected 
    ///               curve. Max dist in surface from exact projection.
    /// \param domain_of_interest If the user do not want a projection onto the whole
    ///                           'surf', he/she can specify the part of the surface
    ///                           to use by limiting the parametric domain to use.
    ///                           If this pointer is left as a null pointer, the whole
    ///                           surface is used.
    /// \return a raw pointer to the generated SplineCurve.  User assumes ownership.
    SplineCurve GO_API *
    projectSpaceCurve(std::shared_ptr<ParamCurve>& space_cv,
                      std::shared_ptr<ParamSurface>& surf,
                      std::shared_ptr<Point>& start_par_pt,
                      std::shared_ptr<Point>& end_par_pt,
                      double epsge,
                      const RectDomain* domain_of_interest = NULL);
    
    /// Lift the parameter_cv onto surf.

    /// 'Lift' a 2D parameter curve onto a surface.  (This means generate a space curve
    /// lying \em on the surface, defined as the curved obtained when evaluating 
    /// the surface at the parameter points given by the 2D curve).  
    /// \param parameter_cv the 2D parameter curve.  Its values should be kept within
    ///                    the parameter domain of 'surf'.
    /// \param surf the surface on which the obtained 3D curve will lie.
    /// \param epsge geometrical tolerance used when generating the curve
    /// \return a raw pointer to the generated spatial SplineCurve.  User assumes 
    ///         ownership.
    SplineCurve GO_API *
    liftParameterCurve(std::shared_ptr<ParamCurve>& parameter_cv,
                       std::shared_ptr<ParamSurface>& surf,
                       double epsge);

    // Assuming dim = 3 (or 2?)
    // Axis defines the start-/end-point of the curve.
    /// Create a circle.

    /// Generate a 3D SplineCurve representing a circle. Assuming dimension is 3.
    /// \param center the center point of the circle,
    /// \param axis a vector pointing from the center point of the circle and towards
    ///             start and end point of the circle.
    /// \param normal the normal to the plane in which the circle lie.
    /// \param radius the radius of the circle to generate.
    /// \return a raw pointer to the generated SplineCurve, representing the
    ///         specified circle.  User assumes ownership.
    SplineCurve GO_API *createCircle(Point center, Point axis, Point normal,
                                     double radius);

    /// Given input of 1-dimensional curve, return the 2-dimensional visualization (x, f(x)).
    std::shared_ptr<Go::SplineCurve> GO_API
    insertParamDomain(const Go::SplineCurve& cv_1d, double knot_tol = 1e-08);

    /// Return the spline curve given by cv(t) + offset_val.
    SplineCurve GO_API *offsetCurve(const SplineCurve& cv, Point offset_val);


}


} // namespace go


#endif // _CURVECREATORS_H

