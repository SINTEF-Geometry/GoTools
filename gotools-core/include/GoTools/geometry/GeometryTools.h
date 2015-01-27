/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#ifndef _GEOMETRYTOOLS_H
#define _GEOMETRYTOOLS_H

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/utils/Array.h"
#include <memory>
#include <vector>
#include "GoTools/utils/config.h"


namespace Go
{

/// This namespace contains free functions operating on objects from the
/// geometry module. These functions are in general not called from
/// the classes of the module, but are considered to be part of the
/// user's programming interface.
namespace GeometryTools
{

    //-----------------------------------------------------------------------
    //
    //  Tool functions in alphabetical order.
    //
    //-----------------------------------------------------------------------


    /// Analyze periodicity of curve based on number of repeating
    /// knots and control points. The return value is -1 if the curve
    /// ends are disjoint, otherwise k if cv is C^k continuous. These are
    /// sufficient but not necessary conditions for periodicity, so it is
    /// possible that a call to analyzePeriodicityDerivs() will yield a
    /// higher degree of periodicity.
    /// \param cv the curve to analyse
    /// \param knot_tol the tolerance used when comparing knot intervals.
    /// \return -1 if the curve ends are disjoint, or k if the curve is 
    ///         proven to be C^k continuous.
    int GO_API
    analyzePeriodicity(const SplineCurve& cv, double knot_tol = 1e-12);

    /// Analyze periodicity of surface based on number of repeating
    /// knots and control points. The return value is -1 if the surface
    /// edges are disjoint, otherwise k if sf is C^k continuous across the
    /// seam. These are sufficient but not necessary conditions for periodicity,
    /// so it is possible that a call to analyzePeriodicityDerivs() will yield a
    /// higher degree of periodicity.
    /// The current implementation is quite slow, and not optimized for speed.
    /// \param sf reference to the SplineSurface to be analyzed
    /// \param direction specify 'direction' to be '0' to check for periodicity in
    ///                  the first parameter direction, or '1' to check
    ///                  the second parameter direction.
    /// \param knot_tol the tolerance used when comparing knot intervals
    /// \return -1 if the surface edges are disjoint, otherwise k if the 
    ///         
    int GO_API
    analyzePeriodicity(const SplineSurface& sf, int direction,
                       double knot_tol = 1e-12);

    /// Analyze periodicity of basis based on number of repeating
    /// knots. The return value is -1 if there is no repeating (or 
    /// order-tuple end knots), up to k if there are (order-1-k)-tuple
    /// end knots and sufficient repeats for a C^k curve. This is
    /// mostly a helper function for the curve and surface analysis
    /// functions.
    /// \param basis reference to the BsplineBasis to analyze
    /// \param knot_tol the tolerance used when comparing knots and knot intervals
    /// \return -1 if there is no repeating (or order-tuple end knots), up to
    ///         k if there are (order-1-k)-tuple end knots and sufficient
    ///         repeats for a C^k curve.
    int GO_API analyzePeriodicity(const BsplineBasis& basis,
                                  double knot_tol = 1e-12);

    /// Analyze periodicity of curve based on evaluating derivatives at both endpoints.
    /// The return value is -1 if the curve
    /// ends are disjoint, otherwise k if cv is C^k continuous.
    /// A maximum of max_derivs derivatives are computed, so the analysis
    /// cannot yield a higher return value then max_derivs.
    /// The tolerance is the amount the point and derivatives are allowed to be
    /// different at the two ends. The default is quite tight!
    /// \param cv the curve to be analyzed for periodicity
    /// \param max_derivs the maximum of computed derivatives (the analysis can
    ///                   not detect continuities beyond this value).
    /// \param tol the tolerance used when comparing points and derivatives for 
    ///            approximate equality
    /// \return -1 if the curve ends are disjoint, k if the curve is proven 
    ///         to be C^k continuous.
    int GO_API analyzePeriodicityDerivs(const ParamCurve& cv,
                                        int max_derivs,
                                        double tol = 1e-14);

    /// Analyze periodicity of surface based on evaluating derivatives
    /// at opposing ends of the surface. The return value is -1 if the surface
    /// edges are disjoint, otherwise k if sf is C^k continuous across the
    /// seam. These are sufficient but not necessary conditions for periodicity,
    /// so it is possible that a call to analyzePeriodicityDerivs() will yield a
    /// higher degree of periodicity.
    /// The current implementation is quite slow, and not optimized for speed.
    /// \param sf the SplineSurface to analyze for periodicity
    /// \param direction the parameter direction to check for periodicity (0 to
    ///                  check in the first parameter direction, 1 to check in the
    ///                  second.
    /// \param max_derivs the maximum of computed derivatives (the analysis cannot
    ///                   detect continuities beyond this value
    /// \param tol the tolerance used when comparing points and derivatives for 
    ///            approximate equality.
    /// \return -1 if the surface edges are disjoint, or k if the surface is 
    ///         proven to be C^k continuous across the seam.
    int GO_API analyzePeriodicityDerivs(const SplineSurface& sf,
                                        int direction,
                                        int max_derivs,
                                        double tol = 1e-14);

    /// Addition of two signed SplineCurves, i.e. this function can also
    /// be used for subtraction. The curves is assumed to live on the same
    /// parameter domain, but may have different knot vectors.
    /// resulting curve = fac1 * crv1 + fac2 * crv2
    /// \param crv1 the first of the curves to be added
    /// \param fac1 multiplicative factor of the first curve
    /// \param crv2 the second of the curves to be  added
    /// \param fac2 multiplicative factor of the second curve
    /// \param num_tol tolerance used when unifying the spline spaces
    ///                in which 'crv1' and 'crv2' lie.
    /// \return shared pointed to a newly created SplineCurve which is the 
    ///         sum of 'crv1' and 'crv2'.
    shared_ptr<SplineCurve> GO_API
    curveSum(const SplineCurve& crv1, double fac1,
	     const SplineCurve& crv2, double fac2,
	     double num_tol = 1e-05);

    // The same requirements as in curveSum().
    shared_ptr<SplineSurface> GO_API
    surfaceSum(const SplineSurface& sf1, double fac1,
	       const SplineSurface& sf2, double fac2,
	       double num_tol = 1e-05);

    /// Rough estimate of the size of a parametric surface
    /// \param srf the surface we want to estimate the size of
    /// \param length_u the estimated average length of the surface, moving along
    ///                 the first parameter
    /// \param length_v the estimated average length of the surface, moving along
    ///                 the second parameter
    void GO_API estimateSurfaceSize(const ParamSurface& srf, double& length_u, 
                                    double& length_v, double *area=NULL);

    /// estimate the length of an iso-curve on the surface
    /// \param srf the surface containing the iso-curve
    /// \param dir_u 'true' if the iso-curve is along the first parameter (second
    ///        parameter fixed), 'false' if it is the first parameter that is fixed.
    /// \param par parameter value for the fixed parameter
    /// \param length returns the estimated length of the iso-curve
    void GO_API estimateIsoCurveLength(const SplineSurface& srf, bool dir_u, 
                                       double par, double& length);

    /// Check if a given spline surface degnerates to a curve within a 
    /// given tolerance
    bool GO_API degenerateToCurve(const SplineSurface& srf, bool dir_u,
                                  double tol);

    /// Make a specified surface boundary exactly degenerate
    void GO_API makeBdDegenerate(SplineSurface& srf, int bd_idx);   // left, right, bottom, top

    /// Compute the union of a set of knot vectors
    void GO_API makeUnionKnots(std::vector<BsplineBasis>& bbasis,
			       double tol, std::vector<double>& union_knots);

    void GO_API makeUnionKnots(std::vector<std::vector<double> >& knots,
			       double tol, std::vector<double>& union_knots);

    /// Check if a curve coefficient is equal to a constant in a specified dimension
    /// provided it already lies close
    bool GO_API checkConstantCoef(SplineCurve& cv, int idx, double val,
                                  double max_dist, double tol);

    /// Modify surface along specified boundary to match a specific constant in 
    /// one direction. Currently only non-rational surfaces
    void GO_API setSfBdCoefToConst(SplineSurface& srf, int bd_idx, int idx_d,
                                   double val, double deg_tol);

    /// Finds the "dominant u- and v-vectors" for a surface, defined to be
    /// the sum of all the vectors pointing from one control point to the
    /// next in the u- and v-directions
    /// \param surface the surface we want to analyze
    /// \param dominant_u returns the dominant u-vector for the surface
    /// \param dominant_v returns the dominant v-vector for the surface
    void GO_API findDominant(const SplineSurface& surface,
                             Vector3D& dominant_u, 
                             Vector3D& dominant_v);

    /// Partition the given curve into segments where each segment is at least
    /// G^n continuous (currently supporting up to G2 continuity).
    /// \param curve the curve to analyze 
    /// \param cont the tolerance used for defining continuity.  cont[0] is the
    ///        tolerance used for defining the G^0 continuity, cont[1] is used
    ///        for checking G^1 continuity, etc.  The length of this vector also
    ///        determines 'n'  ('n' is equal to cont.size() - 1).
    /// \retval gn_joints returns the parameter values where one segment stops
    ///                  and the next begins.  The vector will at least contain
    ///                  two values (the start and end parameter of the curve).
    void GO_API getGnJoints(const ParamCurve& curve,
                            const std::vector<double>& cont,
                            std::vector<double>& gn_joints);

    /// Partition a given CurveLoop into segments where each segment is at least
    /// G^n continuous (currently supporting up to G2 continuity).
    /// \param loop the CurveLoop to analyze
    /// \param cont the tolerance used for defining continuity.  cont[0] is the
    ///             tolerance used for defining the G^0 continuity, cont[1] is used
    ///             for checking G^1 continuity, etc.  The length of this vector also
    ///             determines 'n' ('n' is equal to cont.size() - 1).
    /// \param gn_joints a vector of a vector, reporting the result of the analysis.  
    ///                  The outer vector contains one entry per curve in the 
    ///                  CurveLoop.  This entry contains the parameters for which 
    ///                  that curve must be split in order to obtain G^n segments.
    ///                  NB: the start- and end parameters are NOT included here, 
    ///                  as opposed to the result obtained from the other
    ///                  \ref getGnJoints() function.  The exception to this is when
    ///                  there is an G^n discontinuity at the transition between
    ///                  the curve 'i-1' and curve 'i' in the CurveLoop.  To indicate
    ///                  this, the start parameter value for curve 'i' will be 
    ///                  included in the corresponding entry in 'gn_joints'.
    void GO_API getGnJoints(const CurveLoop& loop,
                            const std::vector<double>& cont,
                            std::vector<std::vector<double> >& gn_joints);

    /// Returns true if the curves are approximately coincident (that
    /// is, they have the same parameter space and overlap in
    /// space). The function uses a (possibly high) number of closest
    /// point calculations to check for spatial coincidence.
    /// \param cv1 the first curve to check for coincidence
    /// \param cv2 the second curve to check for coincidence
    /// \param epsge the tolerance for specifying what we mean by 
    ///              \em approximately coincident
    /// \return 'true' if the curves are found to be approximately coincident,
    ///         'false' otherwise.
    bool GO_API isCoincident(const ParamCurve& cv1, const ParamCurve& cv2,
                             double epsge);

    /// Returns true if any vector difference between neighboring control
    /// points in the u- or v-directions has negative projection on the
    /// given reference vector.
    /// \param surface the SplineSurface containing the control points
    /// \param refvector the vector onto which we carry out the projection
    /// \param eps toleance used when defining  whether a vector has a negative
    ///            projection onto 'refvector' (the scalar product must be
    ///            less than -eps).
    /// \return 'true' if we have found a difference between neighboring control
    ///         points that has a negative projection on the reference vector.
    bool GO_API negativeProj(const SplineSurface& surface,
                             const Array<Vector3D, 2>& refvector,
                             const double eps = 0.0);

    /// Project a 3D parametric curve into a given plane. The curve can be
    /// returned either as a 3D curve lying in the plane or as a 2D
    /// curve. In the latter case, the coordinate system is rotated such
    /// that the given plane normal coincides with the z-axis.
    /// \param incurve the curve to project
    /// \param normal the normal to the plane of projection
    /// \param planar 'true' if we want the returned curve to be a 2D curve.
    /// \return a shared pointer to a newly constructed, planar SplineCurve
    ///         represented as a parametric curve expressing the projection 
    /// of 'incurve' onto the given plane.
     shared_ptr<ParamCurve> GO_API
    projectCurve(shared_ptr<ParamCurve> incurve,
                 const Point& normal,
                 bool planar);

    /// Project a 3D SplineCurve into a given plane. The curve can be
    /// returned either as a 3D curve lying in the plane or as a 2D
    /// curve. In the latter case, the coordinate system is rotated such
    /// that the given plane normal coincides with the z-axis.
    /// \param incurve the curve to project
    /// \param normal the normal to the plane of projection
    /// \param planar 'true' if we want the returned curve to be a 2D curve.
    /// \return a shared pointer to a newly constructed, planar SplineCurve
    ///         expressing the projection of 'incurve' onto the given plane.
    shared_ptr<SplineCurve> GO_API
    projectCurve(const SplineCurve& incurve,
                 const Point& normal,
                 bool planar);

//     /// Project the input pt onto space defined by leg1 & leg2.
//     /// The projected pt is given as: proj_pt = coef1*leg1 + coef2*leg2.
//     void projectPoint(Point pt, Point leg1, Point leg2,
// 		      double& coef1, double& coef2);

    /// Describe a surface as a high-dimensional curve in a given direction.
    /// If the surface is rational, the curve will be non-rational
    /// and living in the homogenous space.
    /// \param surface the surface to express as a curve
    /// \param cv_dir the parameter direction that will be kept when defining 
    ///               the curve (the other one will disappear, as the control
    ///               points in this direction will be lumped together and expressed
    ///               as single control points in a higher-dimensional space.
    ///               'cv_dir' takes either the value '1' (keep the first parameter
    ///               direction) or '2' (keep the second parameter direction).
    /// \return shared pointer to a new SplineCurve, expressing the surface
    ///         as a curve in a high-dimensional space.
    shared_ptr<SplineCurve> GO_API
    representSurfaceAsCurve(const SplineSurface& surface,
			    int cv_dir);

    /// Describe a curve as surface in a given direction. The surface
    /// output is rational if the argument rational is set to 'true'. In that
    /// case the curve is supposed to live in homogenous space.

    /// Describe a curve as a lower-dimensional surface in a given direction.
    /// \param curve the curve that we want to express as a surface
    /// \param cv_dir If this variable is set to '1', then the curve's parameter
    ///               will become the \em first parameter in the generated surface.  
    ///               If it is set to '2', the curve's parameter will become the
    ///               \em second parameter in the generated surface.  Other values
    ///               are illegal.
    /// \param other_bas the BsplineBasis for the additional parameter direction.
    /// \param rational define whether the generated surface shall be specified as 
    ///                 \em rational or not.
    /// \return a shared pointer to a new SplineSurface, expressing the curve
    ///         in a space of lower dimensionality.
    shared_ptr<SplineSurface> GO_API
    representCurveAsSurface(const SplineCurve& curve,
			    int cv_dir,
			    const BsplineBasis& other_bas,
			    bool rational);

    /// Compute the elements of the matrix describing a rotation or a given 
    /// number of radians around a given axis going through the origin.
    /// \param unit_axis_dir this vector defines the axis of rotation.  It must be of
    ///                      unit length, although this is not checked by the function.
    /// \param alpha the angle of rotation, in radians
    /// \return a vector containing the matrix elements of the corresponding rotation
    ///         matrix, stored row-wise.
    std::vector<double> GO_API getRotationMatrix(const Point& unit_axis_dir,
                                                 double alpha);

    /// Rotate the given SplineSurface a certain angle around a given axis.
    /// \param rot_axis the axis of rotation.  It does not have to be normalized, but
    ///                 must of course be nonzero.
    /// \param alpha the angle of rotation, given in radians
    /// \param sf reference to the surface that is to be rotated.
    void GO_API rotateSplineSurf(Point rot_axis, double alpha,
                                 SplineSurface& sf);

    /// Rotate the given SplineCurve a certain angle around a given axis.
    /// \param rot_axis the axis of rotation.  It does not have to be normalized, but
    ///                 must of course be nonzero.
    /// \param alpha the angle of rotation, given in radians
    /// \param cv reference to the curve that is to be rotated
    void GO_API rotateSplineCurve(Point rot_axis, double alpha,
                                  SplineCurve& cv);

    /// Rotate the given LineCloud a certain angle around a given axis.
    /// \param rot_axis the axis of rotation.  It does not have to be normalized, but
    ///                 must of course be nonzero.
    /// \param alpha the angle of rotation, given in radians
    /// \param lc reference to the LineCloud that is to be rotated
    void GO_API rotateLineCloud(Point rot_axis, double alpha, LineCloud& lc);

    /// Rotate the given 3D point a certain angle around a certain axis.
    /// \param rot_axis the axis of rotation.  It does not have to be normalized, but
    ///                 must of course be nonzero.
    /// \param alpha the angle of rotation, given in radians
    /// \param space_pt pointer to the memory location where the 3D coordinates of 
    ///                 the point are stored.  These will be overwritten with the
    ///                 rotated coordinates.
    void GO_API rotatePoint(Point rot_axis, double alpha, double* space_pt);
    /// Rotate the given 3D point a certain angle around a certain axis.
    /// \param rot_axis the axis of rotation.  It does not have to be normalized, but
    ///                 must of course be nonzero.
    /// \param alpha the angle of rotation, given in radians
    /// \param space_pt reference to teh point to be rotated.  This will be overwritten with the
    ///                 rotated coordinates.
    void GO_API rotatePoint(Point rot_axis, double alpha, Point& space_pt);


    /// Split a spline curve into Bezier segments
    void GO_API splitCurveIntoSegments(const SplineCurve& cv,
                                       std::vector<SplineCurve>& seg);


    /// Extract sub patches from the surface given by input parameters.
    std::vector<shared_ptr<SplineSurface> > GO_API
    splitInKinks(const SplineSurface& sf,
		 const std::vector<double>& u_kinks,
		 const std::vector<double>& v_kinks);


    /// Splits a spline surface into Bezier patches
    void GO_API splitSurfaceIntoPatches(const SplineSurface& sf,
                                        std::vector<SplineSurface>& pat);


    /// Surface assumed to be continuous. Return parameter values
    /// failing to achieve G1-continuity.
    void GO_API surfaceKinks(const SplineSurface& sf, double max_normal_angle,
                             std::vector<double>& g1_disc_u, 
                             std::vector<double>& g1_disc_v,
			     bool compute_g1_disc = true);

    /// Find parameter values where a curve is G1- or C1-discontinuous
    void GO_API curveKinks(const SplineCurve& cv, double tol, double ang_tol,
			   std::vector<double>& c1_disconts,
			   std::vector<double>& g1_disconts);

    /// Translate the given SplineSurface by trans_vec.
    void GO_API translateSplineSurf(const Point& trans_vec, SplineSurface& sf);

    /// Translate the given SplineCurve by trans_vec.
    void GO_API translateSplineCurve(const Point& trans_vec, SplineCurve& cv);

    /// Translate the given LineCloud by trans_vec.
    void GO_API translateLineCloud(const Point& trans_vec, LineCloud& lc);

    /// Average specified boundary coefficients between two spline surfaces
    /// to ensure a C0 transition
    void GO_API
    averageBoundaryCoefs(shared_ptr<SplineSurface>& srf1, int bd1,
                         bool keep_first,
                         shared_ptr<SplineSurface>& srf2, int bd2,
                         bool keep_second, bool found_corner1, Point corner1,
                         bool found_corner2, Point corner2, bool opposite);

    /// Make sure that a set of curves live on the same knot vector
    /// tol-equal knots are set equal (i.e. if they differ within tol).
    void GO_API
    unifyCurveSplineSpace(std::vector<shared_ptr<SplineCurve> >& curves,
			  double tol);

    /// Make sure that a set of surfaces live on the same knot vectors
    /// tol-equal knots are set equal (i.e. if they differ within tol).
    /// dir 0 means both, 1 is u, 2 is v
    void GO_API
    unifySurfaceSplineSpace(std::vector<shared_ptr<SplineSurface> >& surfaces,
			    double tol, int dir = 0);

    /// Make sure that a set of surfaces live on the same knot vectors in one
    /// parameter direction.
    /// tol-equal knots are set equal (i.e. if they differ within tol).
    /// Nothing is changed in the other parameter direction (as opposed to
    /// unifySurfaceSplineSpace() where orders are raised to same value for
    /// all surfaces in both directions)
    /// Warning! Objects being pointed to may be recreated inside function.
    /// Remember to update other shared pointers if needed.
    /// \param surfaces Surfaces to end up with common knot vectors
    /// \param tol tolerance for identifying equal knot values
    /// \param unify_u_dir if 'true', unify bases in first parameter direction
    ///        if 'false', unify bases in second parameter direction
    void GO_API
    unifySurfaceSplineSpaceOneDir(std::vector<shared_ptr<SplineSurface> >& surfaces,
				  double tol, bool unify_u_dir);

    /// Join patches with continuity according to input basis_u & basis_v.
    /// The patches are assumed to be in Bezier form, and form a continuous
    /// set. There should be num_u x num_v patches, where num_u = numcoefs_u/order_u
    /// and num_v = numcoefs_v/order_v (as given by basis_u and basis_v).
    shared_ptr<SplineSurface> GO_API
    joinPatches(const std::vector<shared_ptr<SplineSurface> >& patches,
		const SplineSurface& spline_space);

    /// Insert num_knots new knots in the num_knots largest knot
    /// intervals in the B-spline basis basis. The knots are inserted
    /// one by one. Thus, more than one knot may be inserted in one
    /// knot interval of the initial basis. And we do not end up with
    /// the most even distribution of knots possible.
    void GO_API insertKnotsEvenly(BsplineBasis& basis, int num_knots);

    /// Insert num_knots new knots in the num_knots largest knot
    /// intervals in the B-spline basis basis in the parameter
    /// interval [tmin,tmax].  The knots are inserted one by
    /// one. Thus, more than one knot may be inserted in one knot
    /// interval of the initial basis.  And we do not end up with the
    /// most even distribution of knots possible.
    void GO_API insertKnotsEvenly(BsplineBasis& basis, double tmin, double tmax,
                                  int num_knots, double knot_diff_tol = 1e-05);

    /// Return the mid parameter of the largest
    /// parameter interval in a given B-spline basis. If several knot
    /// intervals have the same size, the first is returned.
    double GO_API getKnotAtLargestInterval(const BsplineBasis& basis);

    /// Return the start and end parameter corresponding to the largest
    /// parameter interval in a given B-spline basis. If several knot
    /// intervals have the same size, the first is returned.
    std::pair<double, double> GO_API
    getLargestParameterInterval(const BsplineBasis& basis);

} //namespace GeometryTools

} // namespace Go





#endif // _GEOMETRYTOOLS_H

