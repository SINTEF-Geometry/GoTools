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

#ifndef _BOUNDEDUTILS_H
#define _BOUNDEDUTILS_H


#include "BoundedSurface.h"
#include "SplineSurface.h"
#include "CurveOnSurface.h"
#include "LoopUtils.h"
#include "Plane.h"

#include <memory>


namespace Go {

/// Functions related to the trimming of surfaces, etc.  Also contains functions
/// for spatial transformations of such surfaces (rotation, translation, etc.)
namespace BoundedUtils {

    /// Extract those parts of a given CurveOnSurface where its parameter curve
    /// lies inside the parameter domain of a BoundedSurface.
    /// \param curve the CurveOnSurface that we want to check
    /// \param bounded_surf the BoundedSurface whose parameter domain we will
    ///                     check against
    /// \param epsge geometric tolerance used in calculations
    /// \return a vector containing those segments of 'curve' that have parameter
    ///         descriptions inside the parameter domain of the 'bounded_surf'.
    std::vector<shared_ptr<CurveOnSurface> >
      intersectWithSurface(CurveOnSurface& curve,
			   BoundedSurface& bounded_surf, double epsge,
			   bool only_inner_curves=false);

    /// We have two set of CurveOnSurface s, 'curves1' and 'curves2', and two
    /// BoundedSurface s, 'bd_sf1' and 'bd_sf2'.  We extract those segments of
    /// curves in 'curves1' and 'curves2' that have parameter curves in the 
    /// parameter domains of respective 'bd_sf1' and 'bd_sf2'.  Then we compare
    /// the resulting segments from the two sets against each others, and keep 
    /// those that overlap spatially.  Segments are split in two if necessary
    /// in order to make start and end points from one set coincide with those from
    /// the other set.  'curves1' and 'curves2' are then cleared and filled with the
    /// resulting curve segments.
    /// \param curves1 see above, we suppose that the curves are NOT self-intersecting
    /// \param bd_sf1 see above, we suppose that the underlying surface is the same
    ///               as the one refered to by the curves in 'curves1'.
    /// \param curves2 see above, we suppose that the curves are NOT self-intersecting
    /// \param bd_sf2 see above, we suppose that the underlying surface is the same
    ///               as the one refered to by the curves in 'curves2'.
    /// \param epsge geometric epsilon used in closest point calculations, checks for
    ///              coincidence, etc.
    void intersectWithSurfaces(std::vector<shared_ptr<CurveOnSurface> >& curves1,
			       shared_ptr<BoundedSurface>& bd_sf1,
			       std::vector<shared_ptr<CurveOnSurface> >& curves2,
			       shared_ptr<BoundedSurface>& bd_sf2,
			       double epsge, bool only_inner_curves=false);

    void 
      splitIntersectingCurves(std::vector<shared_ptr<CurveOnSurface> >& cvs1, 
			      std::vector<shared_ptr<CurveOnSurface> >& cvs2, 
			      double epsge, double knot_diff_tol);

    /// Intersect a parametric surface with a plane and fetch the
    /// intersections curves represented as curve on surface 
    std::vector<shared_ptr<CurveOnSurface> >
      getPlaneIntersections(const shared_ptr<ParamSurface>& surf,
			    Point point, Point normal, double epsge,
			    shared_ptr<BoundedSurface>& bounded_sf);

    /// Intersect a parametric surface with a cylinder and fetch the
    /// intersections curves represented as curve on surface 
    std::vector<shared_ptr<CurveOnSurface> >
      getCylinderIntersections(const shared_ptr<ParamSurface>& surf,
			       Point point, Point axis, double radius,  double epsge,
			       shared_ptr<BoundedSurface>& bounded_sf);

    /// Intersect two parametric surfaces and fetch the resulting
    /// intersection curves represented as curve on surface
    void
      getSurfaceIntersections(const shared_ptr<ParamSurface>& surf1,
			      const shared_ptr<ParamSurface>& surf2,
			      double epsge,
			      std::vector<shared_ptr<CurveOnSurface> >& int_cv1,
			      shared_ptr<BoundedSurface>& bounded_sf1,
			      std::vector<shared_ptr<CurveOnSurface> >& int_cv2,
			      shared_ptr<BoundedSurface>& bounded_sf2,
			      bool only_inner_curves=false);

    /// Split a parametric surface by intersecting it with a plane and split along
    /// intersection curves
    std::vector<shared_ptr<BoundedSurface> >
      splitWithPlane(const shared_ptr<ParamSurface>& surf,
		    Point point, Point normal, double epsge);

    /// Split a parametric surface betweem specified parameter values
    std::vector<shared_ptr<BoundedSurface> >
      splitBetweenParams(const shared_ptr<ParamSurface>& surf,
			 Point parval1, Point parval2, double epsge);

    std::vector<shared_ptr<BoundedSurface> >
      splitBetweenParPairs(const shared_ptr<ParamSurface>& surf,
			   std::vector<std::pair<Point,Point> > parvals, 
			   double epsge);

    /// Get the split curves between specified parameter values
    std::vector<shared_ptr<CurveOnSurface> >
      getTrimCrvsParam(const shared_ptr<ParamSurface>& surf,
		       Point parval1, Point parval2, double epsge,
		       shared_ptr<BoundedSurface>& bounded_sf);

    // Given a curve in the parameter domain of the surface, construct
    // the corresponding trimming curve in the surface
    std::vector<shared_ptr<CurveOnSurface> >
      getTrimCrvsPcrv(const shared_ptr<ParamSurface>& surf,
		      shared_ptr<ParamCurve>& pcurve, double epsge,
		      shared_ptr<BoundedSurface>& bounded_sf);

    /// We intersect a parametric surface with a plane, and return the surface(s)
    /// consisting only of the part(s) of the surface that were located on the 
    /// positive side of the intersection.  If there was no intersection, an empty
    /// stl-vector is returned.  The plane is defined by its normal and a point 
    /// located on it.
    /// \param surf the parametric surface.  It must be either a BoundedSurface or 
    ///             a SplineSurface.
    /// \param point a point on the plane to intersect against
    /// \param normal normal of the plane to intersect against
    /// \param epsge geometric tolerance
    /// \return the surface(s) consisting of the part(s) of 'surf' that were located
    ///         on the positive side of the intersection.
    std::vector<shared_ptr<BoundedSurface> >
      trimWithPlane(const shared_ptr<ParamSurface>& surf,
		    Point point, Point normal, double epsge);

    /// must be BoundedSurfaces or SplineSurface
    /// underlying surfaces must be of type SplineSurface
    
    /// If the argument surfaces intersect, and if the intersection curves result in
    /// new boundary loops being defined, then the new surface parts defined within these
    /// domains will be returned.  Otherwise, the return vector will be empty.
    /// \param sf1 the first surface to participate in the intersection
    /// \param sf2 the second surface to participate in the intersection
    /// \param epsge geometrical tolerance to be used in computations
    /// \return a vector with the surfaces representing the parts of the original surfaces
    ///         enclosed by new parametrical loops arising when combining existing loops
    ///         with the curves defined by the intersection.
    std::vector<shared_ptr<BoundedSurface> >
    trimSurfWithSurf(const shared_ptr<ParamSurface>& sf1,
		     const shared_ptr<ParamSurface>& sf2, double epsge);
    

    std::vector<std::vector<shared_ptr<BoundedSurface> > >
	trimSurfsWithSurfs(const std::vector<shared_ptr<ParamSurface> >& sfs1,
			   const std::vector<shared_ptr<ParamSurface> >& sfs2, double epsge);


    std::vector<shared_ptr<BoundedSurface> > 
	trimSurfWithSurfs(shared_ptr<ParamSurface>&  sf,
			  const std::vector<shared_ptr<ParamSurface> >& sfs2, 
			  double epsge);


    /// If surf already is a BoundedSurface, return clone. If SplineSurface,
    /// convert. Otherwise, Error. Return surface is created inside function.

    /// Convert a SplineSurface to a BoundedSurface.  All information is copied, nothing
    /// is shared. 
    /// \param surf the SplineSurface to convert
    /// \param space_epsilon the tolerance assigned to the newly created BoundedSurface
    /// \return (pointer to) a BoundedSurface that represent the same surface as 'surf'.  
    ///         The user assumes ownership of the object.
    BoundedSurface* convertToBoundedSurface(const SplineSurface& surf,
					      double space_epsilon);

    shared_ptr<BoundedSurface> convertToBoundedSurface(shared_ptr<ParamSurface> surf,
						       double space_epsilon);

    /// Given input of partial boundary curves, extract parts of boundary making it a
    /// boundary loop (or more). Input segments expected to be ordered, going in the same direction.
    /// part_bd_cvs should lie on sf (as oppsed to only parts of cvs).

    /// This function tries to complete "partial" boundary loops by filling out
    /// the missing parts using fragments from the domain boundary of a BoundedSurface.
    /// \param sf the surface whose domain boundary will be used
    /// \param part_bnd_cvs a vector of (shared pointers to) curve segments that represent
    ///                     incomplete loops.  Upon function return, this vector will be emptied.
    /// \return a vector contained the loops that the function was able to completely
    ///         close using curve segments from 'part_bnd_cvs' and the domain boundaries of 'sf'.
    std::vector< std::vector< shared_ptr< CurveOnSurface > > >
    getBoundaryLoops(const BoundedSurface& sf, 
		     std::vector< shared_ptr< CurveOnSurface > >&
		     part_bnd_cvs, double eps, int last_split=-1);

    /// Help function for getBoundaryLoops (few sample)
    int checkCurveCoinc(shared_ptr<ParamCurve> cv1, 
			shared_ptr<ParamCurve> cv2, double tol);
      
    /// Help function for getBoundaryLoops (few sample)
    int checkCurveCoincidence(shared_ptr<CurveOnSurface> cv1, 
			      shared_ptr<CurveOnSurface> cv2, 
			      double tol, bool same_orient);
      
   /// All input loops are expected to be simple, lying on surface. They are sorted based
   /// on orientation. No pair of curves with the same orientation may lie inside/outside eachother.
   /// It they do the outer/inner (ccw/cw) loop(s) will be erased.
   /// Furthermore assuming no pair of loops intersect (may touch tangentially).


    /// Define the surfaces that result from trimming a given SplineSurface with a set of 
    /// boundary loops.  Counterclockwise loops define the interior of a parameter domain, while
    /// clockwise loops define holes in the domain. Method assumes all parameter curves exist.
    /// \param loops each entry in the outermost vector represent a vector of curves that together
    ///              specify a loop in 2D parametrical space of the surface 'under_sf'.  These
    ///              are the trim curves.  The curves are expected to be simple, lying on the 
    ///              surface 'under_sf'.  No pair of curve loops with the same orientation may lie
    ///              inside/outside of each other; in that case, the irrelevant loop will be erased.
    ///              Furthermore, we assume that no pair of loops intersect transversally (they are
    ///              still allowed to touch tangentially).
    /// \param under_sf the underlying SplineSurface that we are going to trim with the curves
    ///                 given in 'loops'.
    /// \param epsgeo geometrical tolerance used in computations
    /// \return a vector containing BoundedSurface s that each represent a trimmed part of the
    ///         'under_sf' surface.
    std::vector<shared_ptr<BoundedSurface> >
     createTrimmedSurfs(std::vector<std::vector<shared_ptr<CurveOnSurface> > >& loops,
			shared_ptr<ParamSurface> under_sf, 
			double epsgeo);

    /// Split a given bounded surface according to given trimming curves in
    /// this surface
    /// Notice that the function expects the given curves to actually split
    /// the surface
    std::vector<shared_ptr<BoundedSurface> >
      splitWithTrimSegments(shared_ptr<BoundedSurface> surf,
			    std::vector< shared_ptr< CurveOnSurface > >& bnd_cvs,
			    double eps);

    /// Subtract the part of a trimmed surface corresponding to a given boundary
    /// loop from the surface and return the remaining pieces
    std::vector<shared_ptr<BoundedSurface> >
      subtractSfPart(shared_ptr<BoundedSurface> surf,
		     std::vector< shared_ptr< CurveOnSurface > >& bnd_cvs,
		     double eps);

    /// Find the intersection curve(s) between a parametric surface and a given plane.
    /// The plane is defined by its normal and a point on the plane.
    /// \param surf the surface to intersect with the plane
    /// \param pnt a point lying on the plane
    /// \param normal the normal of the plane
    /// \param geom_tol geometrical tolerance to be used for intersection computations
    /// \return a vector with (shared pointers to) CurveOnSurface s, which represent
    ///         the intersection curves found.
    std::vector<shared_ptr<CurveOnSurface> >
      intersectWithPlane(shared_ptr<ParamSurface>& surf,
			 Point pnt, Point normal, double geom_tol);

    /// Find the intersection curve(s) between a parametric surface and a given plane.
    /// The line is defined by its direction and a point on the line.
    /// \param surf the surface to intersect with the line
    /// \param pnt a point lying on the line
    /// \param dir the direction of the line
    /// \param geom_tol geometrical tolerance to be used for intersection computations
    /// \return a vector with surface parameter and position of intersection.
    std::vector<std::pair<Point, Point> >
      intersectWithLine(shared_ptr<ParamSurface>& surf,
			Point pnt, Point dir, double geom_tol);

    /// Find the intersection curve(s) between a SplineSurface and a given cylinder.
    /// The cylinder is defined by a point on the axis, the axis direction and the radius
    /// \param surf the SplineSurface to intersect with the cylinder
    /// \param pnt a point lying on the cylinder axis
    /// \param vec is the direction of the cylinder axis
    /// \radius is the cylinder radius
    /// \param geom_tol geometrical tolerance to be used for intersection computations
    /// \return a vector with (shared pointers to) CurveOnSurface s, which represent
    ///         the intersection curves found.
    std::vector<shared_ptr<CurveOnSurface> >
      intersectWithCylinder(shared_ptr<ParamSurface>& surf,
			    Point pnt, Point vec, double radius, double geom_tol);

    /// Find the intersction curve(s) between two spline surfaces
    /// \param sf1 the first surface to intersect
    /// \param sf2 the second surface to intersect
    /// \retval int_segments1 a vector of CurveOnSurface s, representing the intersection 
    ///                       curves as lying on 'sf1'.
    /// \retval int_segments2 a vector of CurveOnSurface s, representing the intersection
    ///                       curves as lying on 'sf2'.
    /// \param epsge geometrical tolerance to be used for intersection computations
    void getIntersectionCurve(shared_ptr<ParamSurface>& sf1,
			      shared_ptr<ParamSurface>& sf2,
			      std::vector<shared_ptr<CurveOnSurface> >& int_segments1,
			      std::vector<shared_ptr<CurveOnSurface> >& int_segments2,
			      double epsge);


    /// Translate a given BoundedSurface 
    /// \param trans_vec the vector specifying the translation to apply to the surface
    /// \param bd_sf the surface to translate
    /// \param deg_eps an epsilon value used when determining degenerate boundary loops
    void translateBoundedSurf(Point trans_vec, BoundedSurface& bd_sf,
			      double deg_eps);

    // Rotate a given BoundedSurface
    /// \param rot_axis a vector specifying the axis of rotation
    /// \param alpha the angle of rotation (in radians)
    /// \param bf_sf the surface to rotate
    /// \param deg_eps an epsilon value used when determining degenerate boundary loops
    void rotateBoundedSurf(Point rot_axis, double alpha,
			   BoundedSurface& bf_sf, double deg_eps);

    /// Surface assumed to be continuous. Return parameter values
    /// failing to achieve G1-continuity.
    void trimSurfaceKinks(const BoundedSurface& sf, double max_normal_angle,
			  std::vector<double>& g1_disc_u, 
			  std::vector<double>& g1_disc_v,
			  bool compute_g1_disc = true);

    /// Check loop orientation and fix if necessary. 
    /// NB! It is assumed that the loops has got the correct sequence
    /// \param surf the BoundedSurface to check
    /// \return \c true if the loop orientation was fixed, \c false otherwise
    int checkAndFixLoopOrientation(shared_ptr<BoundedSurface> surf);

    /// Create spline surface description of plane. The plane should
    /// be bounded according to the space_crvs (which should lie in
    /// the plane).
    /// \param plane
    /// \param space_crvs Curves lying in the plane,.
    /// \return Spline surface description of a bounded plane.
    shared_ptr<Go::SplineSurface>
    makeTrimmedPlane(shared_ptr<Go::Plane>& plane,
		     std::vector<shared_ptr<Go::ParamCurve> >&
		     space_crvs);

    /// Change the position of a plane to adapt it to curves supposed to lie
    /// in this plane. Geometry fix related to external models
    void translatePlaneToCurves(shared_ptr<Go::Plane>& plane,
				std::vector<shared_ptr<Go::ParamCurve> >&
				space_crvs);

    /// Fix the boundary loops of a surface in case of inconsistencies.
    void fixInvalidBoundedSurface(shared_ptr<Go::BoundedSurface>& bd_sf,
				  double max_tol_mult = 1.0);

    /// Check if the distance between two curves in the loop is less than
    /// the tolerance epsgeo. Common endpoints between curves are not counted.
    bool loopIsDegenerate(std::vector<shared_ptr<CurveOnSurface> >& loop,
			  double epsgeo);

    bool createMissingParCvs(Go::BoundedSurface& bd_sf);

    bool createMissingParCvs(std::vector<Go::CurveLoop>& bd_loops);

    // The bd_loop should consist of CurveOnSurface's.
    std::vector<std::pair<shared_ptr<Go::Point>, shared_ptr<Go::Point> > >
    getEndParamPoints(const Go::CurveLoop& bd_loop, bool ccw_loop);


} // namespace Go
} // namespace BoundedUtils

#endif // _BOUNDEDUTILS_H
