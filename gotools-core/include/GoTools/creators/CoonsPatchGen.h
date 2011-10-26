//===========================================================================
//                                                                           
// File: CoonsPatchGen.h                                                   
//                                                                           
// Created: Tue Apr  3 13:52:02 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CoonsPatchGen.h,v 1.5 2007-12-04 16:11:39 jbt Exp $
//                                                                           
// Description: Methods for generating surfaces which interpolate curves.
//                                                                           
//===========================================================================

#ifndef _COONSPATCHGEN_H
#define _COONSPATCHGEN_H


#include "GoTools/geometry/CurveLoop.h"


namespace Go {


class SplineSurface;
class SplineCurve;
class ParamCurve;


/// This namespace contains functions used to create a Coons Patch or
/// a Gordon Surface.

namespace CoonsPatchGen {

    // We introduce some exception classes.

    /// Exception class
    class UnKnownError{};

    /// Create a new SplineSurface representing the coons patch as
    /// defined by a loop of four SplineCurves.
    /// \param boundary the CurveLoop describing the boundary of the
    /// coons patch to be generated.  There must be exactly four
    /// curves, they must all be of type 'SplineCurve', and they must
    /// be nonrational.
    /// \return a pointer to a newly created SplineSurface
    /// representing the coons patch.  The user assumes ownership of
    /// the object.
    SplineSurface* createCoonsPatch(const CurveLoop& boundary);


    // bd_curves must contain four curves cross_curves must contain
    // four curves (may be 0 curves) bd_curves and cross_curves must
    // be SplineCurves bd_curves and cross_curves must lie in the same
    // dimension and be nonrational bd curves will be reparameterized
    // to have reasonable tangent values at endpoints and
    // corresponding curves will share parameter interval.
    // neighbour_tol - hvor mye en kantkurve kan flyttes

    /// bd_curves form a loop, all cross curves point into the surface.
    /// cross_curves need not fulfill twist and tangent requirements.
    /// Create a new SplineSurface representing the coons patch
    /// defined by a loop of four boundary curves, and their
    /// respective cross-tangent curves.
    /// \param bd_curves ccw loop defining the boundary of the patch,
    /// size = 4.
    /// \param cross_curves corresponding cross tangent curves, size =
    /// 4, an element may be a NULL pointer.
    /// \param epsge allowed distance between corresponding end points
    /// of bd_curves.
    /// \param kink_tol allowed angle between end tangents of
    /// bd_curves and corresponding end points of cross_curves.
    /// \return pointer to the created Gordon Surface.
    SplineSurface*
    createCoonsPatch(std::vector<std::shared_ptr<ParamCurve> >& bd_curves,
		     std::vector<std::shared_ptr<ParamCurve> >& cross_curves,
		     double epsge, double kink_tol);

    /// Create a Gordon Surface from input curves.  bd_curves forms a
    /// loop, all cross curves point into the surface.  Assumes curves
    /// have been preprocessed, i.e. cross_curves fulfill twist and
    /// tangent requirements.
    /// \param bd_curves ccw loop defining the boundary of the patch,
    /// size = 4.
    /// \param cross_curves corresponding cross tangent curves, size =
    /// 4, an element may be a NULL pointer.
    /// \return pointer to the created Gordon Surface.
    SplineSurface*
    createCoonsPatch(std::vector<std::shared_ptr<SplineCurve> >& bd_curves,
		     std::vector<std::shared_ptr<SplineCurve> >&
		     cross_curves);

    // None of the Gordon Surfaces perform manipulation of the cross
    // tangs.

    /// Create a Gordon Surface from input curves and parameters.  We
    /// require that both boundary curves in one direction (at least)
    /// are given.  If use_param_values == true, we require the curves
    /// to be ordered (with parameters as given by params), with
    /// nmb_u_crvs u-curves up front.  If use_param_values == false,
    /// mesh_curves is reordered as explained in splitMeshCurves() and
    /// sortMeshCurves(). Intersection parameters of curves are found
    /// and set in params, and nmb_u_crvs is also set.
    /// \param mesh_curves the iso-curves defining the Gordon Surface.
    /// \param params the corresponding iso parameters.  Updated if
    /// use_param_values == false.
    /// \param nmb_u_crvs the number of curves parametrized in the
    /// u-direction.  Updated if use_param_values == false.
    /// \param use_param_values whether to use the input parameters.
    /// \return pointer to the created Gordon Surface.
    // @@sbr Include new version which sets up parameters.
    SplineSurface*
    createGordonSurface(std::vector<std::shared_ptr<SplineCurve> >&
			mesh_curves,
			std::vector<double>& params, int& nmb_u_crvs,
			bool use_param_values);

    /// Create a Gordon Surface from input curves and parameters.
    /// cross_index refers to indexing of mesh_curves, which means
    /// that the boundary curves to which it is a tangent must be
    /// included.  All cross_curves are assumed to fulfill twist and
    /// tangent requirements.
    /// \param mesh_curves the iso-curves defining the Gordon Surface.
    /// \param params the corresponding iso parameters.  Updated if
    /// use_param_values == false.
    /// \param nmb_u_crvs the number of curves parametrized in the
    /// u-direction.  Updated if use_param_values == false.
    /// \param cross_curves cross tangent curves for the Gordon
    /// Surface, corresponding to element in mesh_curves.
    /// \param cross_index referring to index element in mesh_curves.
    /// \param use_param_values whether to use the input parameters.
    /// \return pointer to the created Gordon Surface.
    SplineSurface*
    createGordonSurface(std::vector<std::shared_ptr<SplineCurve> >&
			mesh_curves,
			std::vector<double>& params, int& nmb_u_crvs,
			std::vector<std::shared_ptr<SplineCurve> >&
			cross_curves,
			std::vector<int>& cross_index,
			bool use_param_values = true);

    /// Create a Gordon Surface interpolating the input curves in the input
    /// parameters.
    /// mesh_curves assumed to be ordered (i.e. u-curves first, with ascending
    /// parameter values), and missing boundary curve(s) in at most one direction.
    /// Expecting end conditions are satisfied.
    /// \param mesh_curves iso-curves for Gordon Surface.
    /// \param params iso-parameters for the mesh_curves.
    /// \param nmb_u_crvs the number of curves parametrized in the u-direction.
    ///                   May alter as the parameter directions may swap.
    /// \param cross_curves the cross tangent curves along iso-curves for the Gordon Surface.
    /// \param cross_index index in mesh_curves of corresponding boundary curve.
    /// \return pointer to the created Gordon Surface.
    SplineSurface*
    doCreateSurface(std::vector<std::shared_ptr<SplineCurve> >& mesh_curves,
		    std::vector<double>& params, int& nmb_u_crvs,
		    std::vector<std::shared_ptr<SplineCurve> >& cross_curves,
		    std::vector<int>& cross_index);

    /// Create a lofting surface based on the input curves.
    /// \param first_curve iterator to first iso-curve in the lofted surface.
    /// \param nmb_crvs the number of curves referred to by first_curve.
    /// \return pointer to the created lofting surface.
    SplineSurface* loftSurface(std::vector<std::shared_ptr<SplineCurve> >::iterator
			       first_curve, int nmb_crvs);

    /// Create a lofting surface interpolating the input curves in the input
    /// parameters.
    /// \param first_curve iterator to first iso-curve in the lofted surface.
    /// \param first_param iso parameter to corresponding curve referred to by first_curve.
    /// \param nmb_crvs the number of curves referred to by first_curve.
    /// \return pointer to the created lofting surface.
    SplineSurface* loftSurface(std::vector<std::shared_ptr<SplineCurve> >::iterator
			       first_curve,
			       std::vector<double>::iterator first_param,
			       int nmb_crvs);

    /// Create a lofting surface interpolating the input curves and cross tangent curves
    /// in the input parameters. All curves assumed to share spline space.
    /// \param first_curve iterator to first iso-curve in the lofted surface.
    /// \param first_param iso parameter to corresponding curve referred to by first_curve.
    /// \param nmb_crvs the number of curves referred to by first_curve.
    /// \param first_cross_curve iterator to first cross tangent curve in the lofted surface.
    /// \param cross_index referring to index of corresponding boundary curve.
    /// \return pointer to the created lofting surface.
    SplineSurface*
      loftSurface(std::vector<std::shared_ptr<SplineCurve> >::iterator first_curve,
		  std::vector<double>::iterator first_param,
		  int nmb_crvs,
		  std::vector<std::shared_ptr<SplineCurve> >::iterator first_cross_curve,
		  std::vector<int>& cross_index);

    /// Make tensor product surface which interpolates given grid points.
    /// All curves assumed to share spline space.
    /// \param mesh_curves curves to be interpolated in the Gordon Surface, both u- and
    ///                    v-curves.
    /// \param params the iso parameters for the curves.
    /// \param nmb_u_crvs the number of curves parametrized in the u-direction.
    /// \param cross_curves cross tangent curves for the Gordon Surface.
    /// \param cross_index referring to index of corresponding boundary curve.
    /// \return pointer the created lofting surface.
    SplineSurface* tpSurface(const std::vector<std::shared_ptr<SplineCurve> >& mesh_curves,
			      std::vector<double> params, int nmb_u_crvs,
			      const std::vector<std::shared_ptr<SplineCurve> >& cross_curves,
			      std::vector<int>& cross_index);

    /// Given input iso-curves, the curves are analyzed and the ordering altered such
    /// that the nmb_u_crvs first elements are u-curves, and the rest are v-curves.
    /// Calculated iso parameter for the curves is returned in params.
    /// \param mesh_curves iso-curves for a surface, sorted inside function.
    /// \param params calculated parameters for the iso-curves.
    /// \param nmb_u_crvs the number of curves parametrized in the u-direction.
    /// \param cross_index elements referring to mesh_curves, updated inside function.
    /// \param epsgeo geometrical tolerance defining intersections between mesh_curves.
    void splitMeshCurves(std::vector<std::shared_ptr<SplineCurve> >& mesh_curves,
			 std::vector<double>& params, int& nmb_u_crvs,
			 std::vector<int>& cross_index, double epsgeo);

    /// Given that mesh_curves has been regrouped with u-curves up front,
    /// we sort the vector according to values in params.
    void sortMeshCurves(std::vector<std::shared_ptr<SplineCurve> >& mesh_curves,
			std::vector<double>& params, int nmb_u_crvs,
			std::vector<int>& cross_index);


    /// Prepare cross tangents for surface creation.
    /// curves contains boundary_curve, cross_curve; one edge at a time
    /// (iedge*2 curves). boundary curves form a loop. All curves share
    /// spline space, and all cross curves point inwards = into the surface.
    /// A missing cross_curve is indicated by a null pointer.
    /// Output mod_cross_curves fulfills tangent and twist conditions.
    void getCrossTangs(const std::vector<std::shared_ptr<SplineCurve> >& curves,
		       std::vector<std::shared_ptr<SplineCurve> >& mod_cross_curves,
		       double tol1, double tol2);

    /// Generates missing cross boundary curves. In case of missing cross_crv,
    /// we're expecting a NULL pointer.
    /// boundary_crvs form a loop. cross_crvs share spline space, and all cross
    /// curves point inwards = into the surface. A missing cross_curve is
    /// indicated by a null pointer. The added cross curve will also point inwards,
    /// following parametrization given by boundary curves.
    void addMissingCrossCurves(const std::vector<std::shared_ptr<SplineCurve> >& bnd_curves,
			       std::vector<std::shared_ptr<SplineCurve> >& cross_crvs);

    /// Find blending functions used to blend two derivative
    /// along some boundary curves corresponding to a surface,
    /// into a cross derivative curve pr edge.
    /// Get curves for blending of two tangent curves into one cross tangent curve.
    /// Size of curves = iedge * 3 (i.e. bd_curve, cross_curve, tangent curve).
    /// Boundary curves form a loop, and all curves share orientation and
    ///line space. cross-curves point inwards.
    void getTangBlends(std::vector<std::shared_ptr<SplineCurve> >& curves, int iedge,
		       std::vector<std::shared_ptr<SplineCurve> >& blend_functions);

    /// Project the etang vector onto the plane defined by evecu and evecv.
    /// The projected vector will be represented as (*coef1)*evecu + (*coef2)*evecv.
    /// \param evecu one of the vectors spanning the plane.
    /// \param evecv one of the vectors spanning the plane.
    /// \param etang vector to be projected.
    /// \param idim the geometric dimension of the space.
    /// \param isign sign with wich etang is to be multiplied.
    /// \param coef1 scalar corresponding to evecu.
    /// \param coef2 scalar corresponding to evecv.
    void blendcoef(double evecu[],double evecv[],double etang[],
		   int idim,int isign,double *coef1, double *coef2);

    /// Hermite interpolate the input points.
    /// \param econd array containing the points and derivatives to be interpolated.
    ///              The ordering is <start_pt, (start_tan,) end_der, end_pt>
    ///              The array will contain the coefficients of the spline curve
    ///              fulfilling the input requirements.
    /// \param icond number of interpolation conditions.
    /// \param hasder1 Whether derivative is given in the start of the interval
    /// \param astart start parameter of the produced curve.
    /// \param aend end parameter of the produced curve.
    /// \param idim dimension of the geometric space.
    void hermit(double econd[], int icond, bool hasder1,
		double astart, double aend, int idim);

    /// Make sure that the end points of the existing cross_tangent_curves match
    /// the end derivatives of the corresponding bd_curves.
    /// Curves are oriented CCW.
    /// \param bd_curves curves along the boundary, ccw orientated. Will not be altered.
    /// \param cross_curves the corresponding cross tangent curves (same orientation).
    ///                     The curve objects referred to by the pointers may be altered.
    void fixCrossEndPts(const std::vector<std::shared_ptr<SplineCurve> >& bd_curves,
			const std::vector<std::shared_ptr<SplineCurve> >& cross_curves);


    /// Calculate iso parameters for the input curves. The curves are expected to be
    /// ordered, i.e. corresponding to increasing iso parameters.
    /// Curves are given iso-parameters in the range 0.0 to param_length.
    /// \param first_curve iterator to first iso curve.
    /// \param nmb_crvs the number of input curves.
    /// \param param_length length of parameter domain.
    /// \param params the computed iso parameters for the input curves.
    void
    makeLoftParams(std::vector<std::shared_ptr<SplineCurve> >::const_iterator first_curve,
		   int nmb_crvs, double param_length, std::vector<double>& params);

    /// Check length of tangent vectors at the endpoints of a curve compared to
    /// the size of the curve. If the vectors are too long, reparametrize the curve
    /// and all the other curves in vector.
    /// All curves expected to share parametrization.
    /// First curve is a bnd curve. Optional additional curves may be cross tangent
    /// curve. Function reparametrizes bnd curve so that end tangents have length
    /// equal to aconst. Typical value is 1/3 of edge-length.
    /// \param curves input and output curves.
    /// \param aconst scalar defining end tangent length.
    void
    reparamBoundaryCurve(std::vector<std::shared_ptr<SplineCurve> >& curves, double aconst);

} // namespace CoonsPatchGen


} // namespace Go


#endif // _COONSPATCHGEN_H
