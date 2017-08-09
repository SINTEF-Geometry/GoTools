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

#ifndef _BOUNDEDSURFACE_H
#define _BOUNDEDSURFACE_H


#include <memory>
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/utils/config.h"


namespace Go
{

  class SplineSurface;

/// The class representing trimmed surfaces in Go. Surface should be
/// connected.
class GO_API BoundedSurface : public ParamSurface
{
public:

    /// Create an empty BoundedSurface that can be assigned or read() into
    BoundedSurface();

    /// Create a BoundedSurface by specifying the underlying surface and a loop
    /// of curves that specifies the trimming of the surface.
    /// \param surf the created BoundedSurface will represent a trimmed version of 
    ///             this surface
    /// \param loop a vector of CurveOnSurface s that together describe the boundary
    ///             that defines the trimming of the  surface.  The curves in this
    ///             vector should all lie on the surface in question, and when placed
    ///             head-to-tail they should form a closed loop with counterclockwise
    ///             orientation.
    /// \param space_epsilon geometrical tolerance used when treating the loops.
    /// \param fix_trim_cvs the constructor may alter loop curves if they are not
    ///                     consistent.
    BoundedSurface(shared_ptr<ParamSurface> surf,
		   std::vector<shared_ptr<CurveOnSurface> > loop,
		   double space_epsilon,
		   bool fix_trim_cvs = true);


    /// Create a BoundedSurface by specifying the underlying surface and a number of 
    /// loops of curves that specify the trimming of the surface.
    /// \param surf the created BoundedSurface will represent a trimmed version of 
    ///             this surface
    /// \param loops each entry in 'loop is a vector of CurveOnSurface s that describe 
    ///             a closed loop forming a part of the trimmed surface's boundary.  (Since
    ///             the surface may have internal holes, more than one loop might be 
    ///             required to describe its boundary).  The curve loops should all
    ///             lie on the surface in question.  The first entry in this vector
    ///             describes the outermost boundary, which should be oriented 
    ///             counterclockwise.  The other entries represent holes, and should
    ///             be oriented clockwise.
    /// \param space_epsilon geometrical tolerance used when treating the loops.
    /// \param fix_trim_cvs the constructor may alter loop curves if they are not
    ///                     consistent.
    BoundedSurface(shared_ptr<ParamSurface> surf,
		   std::vector<std::vector<shared_ptr<CurveOnSurface> > > loops,
		   double space_epsilon,
		   bool fix_trim_cvs = true);

    /// Create a BoundedSurface by specifying the underlying surface
    /// and a number of loops of curves that specify the trimming of
    /// the surface.
    /// \param surf the created BoundedSurface will represent a
    ///             trimmed version of this surface
    /// \param loops each entry in 'loops' is a vector of CurveOnSurface
    ///             s that describe a closed loop forming a part of
    ///             the trimmed surface's boundary.  (Since the
    ///             surface may have internal holes, more than one
    ///             loop might be required to describe its boundary).
    ///             The curve loops should all lie on the surface in
    ///             question.  The first entry in this vector
    ///             describes the outermost boundary, which should be
    ///             oriented counterclockwise.  The other entries
    ///             represent holes, and should be oriented clockwise.
    /// \param space_epsilons geometrical tolerances used when treating
    ///             the loops.
    /// \param fix_trim_cvs the constructor may alter loop curves if they are not
    ///                     consistent.
    BoundedSurface(shared_ptr<ParamSurface> surf,
		   std::vector<std::vector<shared_ptr<CurveOnSurface> > > loops,
		   std::vector<double> space_epsilons,
		   bool fix_trim_cvs = true);

    /// Create a bounded surface from a non-trimmed one
    BoundedSurface(shared_ptr<ParamSurface> surf,
		   double space_epsilon);

    /// Create a bounded surface from a non-trimmed one
    BoundedSurface(shared_ptr<ParamSurface> surf,
		   std::vector<CurveLoop>& loops);

    /// Create a bounded surface from a non-trimmed one
    BoundedSurface(shared_ptr<ParamSurface> surf,
		   std::vector<shared_ptr<CurveLoop> >& loops);

    /// Virtual destructor ensures safe inheritance
    virtual ~BoundedSurface();


    // From Streamable

    // @afr: These should not be called!
    /// read this BoundedSurface from a stream
    virtual void read (std::istream& is);

    void read (std::istream& is,
	       bool fix_trim_cvs);

    /// write this BoundedSurface to a stream
    virtual void write (std::ostream& os) const;

    // From GeomObject

    /// Return the object's bounding box
    virtual BoundingBox boundingBox() const;

    /// Return the dimension of the space in which the object lies (usually 2 or 3)    
    virtual int dimension() const;

    /// Return the class type identifier of type BoundedSurface
    virtual ClassType instanceType() const;

    /// Return the class type identifier of type BoundedSurface
    static ClassType classType()
    { return Class_BoundedSurface; }

    /// clone this BoundedSurface and return a pointer to the clone
    virtual BoundedSurface* clone() const;

    /// Return the spline surface represented by this surface, if any
    virtual SplineSurface* asSplineSurface() 
    {
      return surface_->asSplineSurface();
    }

    // From ParamSurface

    /// Creates a DirectionCone covering all normals to this surface.
    /// \return a DirectionCone (not necessarily the smallest) containing all normals 
    ///         to this surface.
    virtual DirectionCone normalCone() const;

    /// Creates a DirectionCone covering all tangents to 
    /// this surface along a given parameter direction.
    /// \param pardir_is_u if 'true', then the DirectionCone will be defined on basis 
    ///        of the surface's tangents along the first parameter direction.  Otherwise
    ///        the second parameter direction will be used.
    /// \return a DirectionCone (not necessarily the smallest) containing all tangents
    ///         to this surface along the specified parameter direction.
    virtual DirectionCone tangentCone(bool pardir_is_u) const;

    /// Return the parameter domain of the surface.  This may be a simple
    /// rectangular domain (\ref RectDomain) or any other subclass of
    /// \ref Domain (such as CurveBoundedDomain, found in the
    /// \c sisl_dependent module).
    /// \return a Domain object describing the parametric domain of the surface
    virtual const CurveBoundedDomain& parameterDomain() const;

    /// Get a rectangular parameter domain that is guaranteed to contain the
    /// surface's \ref parameterDomain().  It may be the same.  There is no
    /// guarantee that this is the smallest domain containing the actual domain.
    /// \return a RectDomain that is guaranteed to include the surface's total
    ///         parameter domain.
    virtual RectDomain containingDomain() const;

    /// Check if a parameter pair lies inside the domain of this surface
    virtual bool inDomain(double u, double v, double eps=1.0e-4) const;

    /// Check if a parameter pair lies inside the domain of this surface
    /// return value = 0: outside
    ///              = 1: internal
    ///              = 2: at the boundary
    virtual int inDomain2(double u, double v, double eps=1.0e-4) const;

    /// Check if a parameter pair lies at the boundary of this surface
    virtual bool onBoundary(double u, double v, double eps=1.0e-4) const;


    virtual Point closestInDomain(double u, double v) const;

    /// Returns the anticlockwise, outer boundary loop of the surface.
    /// \param degenerate_epsilon edges whose length is smaller than this value
    ///        are ignored.
    /// \return a CurveLoop describing the anticlockwise, outer boundary loop of 
    ///         the surface.
    virtual CurveLoop outerBoundaryLoop(double degenerate_epsilon
					  = DEFAULT_SPACE_EPSILON) const;

    /// Returns the anticlockwise outer boundary loop of the surface, together with 
    /// clockwise loops of any interior boundaries, such that the surface always is
    /// 'to the left of' the loops.
    /// \param degenerate_epsilon edges whose length is smaller than this value are
    ///                           ignored.
    /// \return a vector containing CurveLoops.  The first of these describe the 
    ///         outer boundary of the surface (counterclockwise), whereas the others 
    ///         describe boundaries of interior holes (clockwise).
    virtual std::vector<CurveLoop> allBoundaryLoops(double degenerate_epsilon
						      = DEFAULT_SPACE_EPSILON) const;

    /// Returns the anticlockwise outer boundary loop of the surface, together with 
    /// clockwise loops of any interior boundaries, such that the surface always is
    /// 'to the left of' the loops.  This function works like \ref allBoundaryLoops(),
    /// except that it includes degenerate edges.
    /// \return vector containing CurveLoops.  The first of these describe the
    ///         outer boundary of the surface (counterclockwise) whereas the others
    ///         describe boundaries of interior holes (clockwise).
    std::vector<CurveLoop> absolutelyAllBoundaryLoops() const;


    /// Evaluates the surface's position for a given parameter pair.
    /// \param pt the result of the evaluation is written here 
    /// \param upar the first parameter
    /// \param vpar the second parameter
    virtual void point(Point& pt, double upar, double vpar) const;


    /// Evaluates the surface's position and a certain number of derivatives
    /// for a given parameter pair.
    /// \param pts the vector containing the evaluated values.  Its size must be 
    ///            set by the user prior to calling this function, and should be
    ///            equal to (derivs+1) * (derivs+2) / 2.  Upon completion of the
    ///            function, its first entry is the surface's position at the given
    ///            parameter pair.  Then, if 'derivs' > 0, the two next entries will
    ///            be the surface tangents along the first and second parameter 
    ///            direction.  The next three entries are the second- and cross 
    ///            derivatives, in the order (du2, dudv, dv2), and similar for 
    ///            even higher derivatives.
    /// \param upar the first parameter 
    /// \param vpar the second parameter
    /// \param derivs number of requested derivatives
    /// \param u_from_right specify whether derivatives along the first parameter are
    ///                     to be calculated from the right ('true', default) or from
    ///                     the left ('false')
    /// \param v_from_right specify whether derivatives along the second parameter are
    ///                     to be calculated from the right ('true', default) or from
    ///                     the left ('false')
    /// \param resolution tolerance used when determining whether parameters are located 
    ///                   at special values of the parameter domain (in particualar; knot
    ///                   values in case of spline objects.
    virtual void point(std::vector<Point>& pts, 
		       double upar, double vpar,
		       int derivs,
		       bool u_from_right = true,
		       bool v_from_right = true,
		       double resolution = 1.0e-12) const;

    // /// Evaluates the surface's position and a certain number of derivatives for 
    // /// a given parameter pair.
    // /// \param pts the vector containing the evaluated values.  Its size must be 
    // ///            set by the user prior to calling this function, and should be
    // ///            equal to (derivs+1) * (derivs+2) / 2.  Upon completion of the
    // ///            function, its first entry is the surface's position at the given
    // ///            parameter pair.  Then, if 'derivs' > 0, the two next entries will
    // ///            be the surface tangents along the first and second parameter 
    // ///            direction.  The next three entries are the second- and cross 
    // ///            derivatives, in the order (du2, dudv, dv2), and similar for 
    // ///            even higher derivatives.
    // /// \param upar the first parameter 
    // /// \param vpar the second parameter
    // /// \param derivs number of requested derivatives
    // virtual void point(std::vector<Point>& pts, 
    // 		       double upar, double vpar,
    // 		       int derivs) const;

    //    using ParamSurface::point;
    /// Evaluates the surface normal for a given parameter pair
    /// \param n the computed normal will be written to this variable
    /// \param upar the first parameter
    /// \param vpar the second parameter
    virtual void normal(Point& n, double upar, double vpar) const;

    /// Evaluate points in a grid
    /// The nodata value is applicable for bounded surfaces
    /// and grid points outside the trimming loop(s) will
    /// get this value
    virtual void evalGrid(int num_u, int num_v, 
			  double umin, double umax, 
			  double vmin, double vmax,
			  std::vector<double>& points,
			  double nodata_val = -9999) const;

    /// Fetch an arbitrary internal point in the surface
    /// Used for localization purposes
    virtual Point getInternalPoint(double& u, double& v) const;

    /// Get the curve(s) obtained by intersecting the surface with one of its constant
    /// parameter curves.  For surfaces without holes, this will be the parameter curve
    /// itself; for surfaces with interior holes this may be a collection of several, 
    /// disjoint curves.  
    /// \param parameter parameter value for the constant parameter (either u or v)
    /// \param pardir_is_u specify whether the \em moving parameter (as opposed to the 
    ///                    \em constant parameter) is the first ('true') or the second
    ///                    ('false') one.
    /// \return a vector containing shared pointers to the obtained, newly constructed
    ///          constant-parameter curves.
    virtual std::vector<shared_ptr<ParamCurve> >
    constParamCurves(double parameter, bool pardir_is_u) const;

    /// Get the surface(s) obtained by cropping the parameter domain of this surface
    /// between given values for the first and second parameter.  In general, for 
    /// surfaces with no interior holes, the result will be \em one surface; however,
    /// for surfaces with interior holes, the result might be \em several \em disjoint
    /// surfaces.
    /// \param from_upar lower value for the first parameter in the subdomain
    /// \param from_vpar lower value for the second parameter in the subdomain
    /// \param to_upar upper value for the first parameter in the subdomain
    /// \param to_vpar upper value for the second parameter in the subdomain
    /// \param fuzzy tolerance used when determining intersection with interior 
    ///        boundaries
    /// \return a vector contained shared pointers to the obtained, newly constructed
    ///         sub-surfaces.
    virtual std::vector<shared_ptr<ParamSurface> >
    subSurfaces(double from_upar, double from_vpar,
		double to_upar, double to_vpar,
		double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// Determine the parameter value of the start of the 'next
    /// segment' from a parameter value, along a given parameter
    /// direction.  A 'segment' is here defined as a parameter
    /// interval in which there will be no discontinuities in
    /// derivatives or other artifacts.  For spline objects, a segment
    /// will typically be the interval between two consecutive,
    /// non-coincident knots.
    /// \param dir the parameter direction in which we search for the
    /// next segment (0 or 1)
    /// \param par the parameter value starting from which we search
    /// for the start value of the next segment
    /// \param forward define whether we shall move forward ('true')
    /// or backwards when searching along this parameter
    /// \param tol tolerance used for determining whether the 'par' is
    /// already located \em on the next segment value
    /// \return the value of the start value of the next segment (or
    /// the end of the previous segment, if we are moving
    /// backwards...)
    virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const;

    /// Iterates to the closest point to pt on the surface. 
    /// \param pt the point to find the closest point to
    /// \param clo_u u parameter of the closest point
    /// \param clo_v v parameter of the closest point
    /// \param clo_pt the geometric position of the closest point
    /// \param clo_dist the distance between pt and clo_pt
    /// \param epsilon parameter tolerance (will in any case not be higher than
    ///                sqrt(machine_precision) x magnitude of solution
    /// \param domain_of_interest pointer to parameter domain in which to search for 
    ///                           closest point. If a NULL pointer is used, the entire
    ///                           surface is searched.
    /// \param seed pointer to parameter values where iteration starts.
    virtual void closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v, 
			      Point&       clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      const RectDomain* domain_of_interest = NULL,
			      double   *seed = 0) const;


    /// Iterates to the closest point to pt on the boundary of the surface.
    /// \see closestPoint()
    virtual void closestBoundaryPoint(const Point& pt,
				      double&        clo_u,
				      double&        clo_v, 
				      Point&       clo_pt,
				      double&        clo_dist,
				      double epsilon,
				      const RectDomain* domain_of_interest = NULL,
				      double *seed = 0) const;

    /// Get the boundary curve segment between two points on the boundary, as 
    /// well as the cross-tangent curve.  If the given points are not positioned 
    /// on the same boundary (within a certain tolerance), no curves will be created.
    /// \b NB: This function has not yet been implemented!
    /// \param pt1 the first point on the boundary, given by the user
    /// \param pt2 the second point on the boundary, given by the user
    /// \param epsilon the tolerance used when determining whether the given points
    ///                are lying on a boundary, and if they do, whether they both lie
    ///                on the \em same boundary.
    /// \param cv upon return, this will point to a newly created curve representing
    ///           the boundary curve between 'pt1' and 'pt2'.  The user assumes ownership
    ///           of this object and is responsible for its deletion.  No curve is created
    ///           if the given points are not found to lie on the same boundary.
    /// \param crosscv upon return, this will point to a newly created curve representing
    ///                the cross-boundary curve between 'pt1' and 'pt2'  The user assumes
    ///                ownership of this object and is responsible for its deletion.
    ///                The direction is outwards from the surface.
    ///                No curve is created if the given points are not found to lie on the
    ///                same boundary.
    /// \param knot_tol tolerance used when working with the knot-vector, to specify how
    ///                 close a parameter value must be to a knot in order to be considered
    ///                 'on' the knot.
    virtual void getBoundaryInfo(Point& pt1, Point& pt2, 
				 double epsilon, SplineCurve*& cv,
				 SplineCurve*& crosscv, double knot_tol = 1e-05) const;


    /// Get the boundary curve segment between two points on the same
    /// boundary loop.  If the given points are not positioned on the
    /// same boundary loop (within a certain tolerance), no curves
    /// will be retuned.
    /// \param pt1 the first point on the boundary, given by the user
    /// \param pt2 the second point on the boundary, given by the user
    /// \retval bd_cvs upon return, this will contain shared pointers to curves that, 
    ///                taken consecutively, describe the requested boundary segment in 
    ///                its entirety.
    void getBoundaryInfo(Point& pt1, Point& pt2,
			 std::vector<shared_ptr<CurveOnSurface> >& bd_cvs) const;

    /// Turns the direction of the normal of the surface.
    virtual void turnOrientation();

    /// Reverses the direction of the basis in input direction.
    /// \b NB: This function has not yet been implemented!
    /// \param direction_is_u if 'true', the first parameter direction
    ///                       will be reversed, otherwise, the second
    ///                       parameter direction will be reversed
    virtual void reverseParameterDirection(bool direction_is_u);

    // If a segment in a loop in boundary_loops_ is not G1, curve is split.
    // Function to be called prior to a topology builder relying on smooth segments.

    /// This function processes all the curves that participate in
    /// defining the surface's (trimmed) boundary.  Those curves that
    /// are not G1 within a certain tolerance are split into
    /// severalcurves, so that all G1-discontinuities will end up \em
    /// between consecutive curve segments.
    /// \param kink the tolerance to use for checking G1 continuity
    void makeBoundaryCurvesG1(double kink);

    /// This function processes all the curves that participate in
    /// defining the surface's (trimmed) boundary.  Those curves that
    /// are smaller than a given tolerance is merged with one of the
    /// adjacent curves if possible
    /// \param gap positional tolerance used to check for possibility of merge
    /// \param neighbour tolerance used for checking if a curve is small
    /// \param kink angular tolerance used to check for possibility of merge
    void removeSmallBoundaryCurves(double gap, double neighbour,
				   double kink);

     /// Swaps the two parameter directions
    virtual void swapParameterDirection();

    /// Compute the total area of this surface up to some tolerance
    /// \param tol the relative tolerance when approximating the area, i.e.
    ///            stop iteration when error becomes smaller than
    ///            tol/(surface area)
    /// \return the area calculated
    /// NB! Intermediate solution with lower accuracy
    virtual double area(double tol) const;

     /// Change the parameter domain of the underlying surface, and modify the
    /// boundary loops with respect to this change
    /// \param u1 new start value of first parameter
    /// \param u2 new end value of first parameter
    /// \param v1 new start value of second parameter
    /// \param v2 new end value of second parameter
    virtual void setParameterDomain(double u1, double u2, double v1, double v2);

    // If a boundary loop is represented as a single curve, it is split into 3 parts.
    // Handy tue to current limitations in topology analysator.

    /// Split all boundary loops defined by only \em curve into three parts.
    /// (This somewhat exotic member function is included due to its handiness with the
    /// GoTools topology analysator).
    void splitSingleLoops();

    // Access functions

    /// Get a pointer to the underlying surface
    /// \return shared pointer to the underlying surface
    shared_ptr<ParamSurface> underlyingSurface()
    { return surface_; }

    /// Get a pointer to the underlying surface
    /// \return shared pointer to the underlying surface
    shared_ptr<const ParamSurface> underlyingSurface() const
    { return surface_; }

    /// Check if the final underlying surface is a spline surface
    /// and in that case return this surface
    bool hasUnderlyingSpline(shared_ptr<SplineSurface>& srf);

    /// Get the number of boundary loops that describe the trimmed surface.
    int numberOfLoops() const
    { return (int)boundary_loops_.size(); }

    /// Get a shared pointer to a specific boundary loop
    shared_ptr<CurveLoop> loop(int idx)
      { return boundary_loops_[idx]; }

    /// Get the space-curve resulting from fixing one of the surface's
    /// parameters and moving the other along its allowed range
    /// (inside the trimmed domain).  If this results in several
    /// disjoint curves, an exception is thrown.
    /// \param parameter the parameter value of the fixed parameter
    /// \param direction_is_u if 'true' then the "free" parameter will be the first one,
    ///                       and the second parameter will be fixed.  If 'false', it is
    ///                       the other way around.
    /// \return a newly created SplineCurve representing the requested space-curve.
    ///         The ownership is assumed by the user.
    SplineCurve* constParamCurve(double parameter, bool direction_is_u) const;

    /// Query whether any of the four boundaries of the \em underlying \em surface
    /// are degenerate (zero length) within a certain tolerance.  In the below, we refer
    /// to 'u' as the first parameter and 'v' as the second.
    /// \param b 'true' upon return of function if the boundary (v = v_min) is degenerate
    /// \param r 'true' upon return of function if the boundary (v = v_max) is degenerate
    /// \param t 'true' upon return of function if the boundary (u = u_min) is degenerate
    /// \param l 'true' upon return of function if the boundary (u = u_max) is degenerate
    /// \param tolerance boundaries are considered degenerate if their length is shorter
    ///        than this value, given by the user
    /// \return 'true' if at least one boundary curve was found to be degenerate, 'false'
    ///         otherwise.
    virtual bool isDegenerate(bool& b, bool& r,
			      bool& t, bool& l, double tolerance) const;


    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    /// Return surface corners, i.e joints betwee trimming curves, 
    /// geometric and parametric points in that sequence
    virtual void 
      getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;

    virtual void setIterator(IteratorType type)
	{
	    iterator_ = type;
	    surface_->setIterator(type);
	}

    /// Check if the current surface is trimmed along constant parameter curves
    virtual bool isIsoTrimmed(double tol) const;

    /// Fetch
    shared_ptr<ParamSurface> getIsoTrimSurface(double tol) const;

    /// Check if the loop orientation is tested and corrected
    bool orientationIsSet()
	{
	    return (loop_fixed_.size() == boundary_loops_.size());
	}

    /// Check the status of the loop orientation. If a fix has been performed
    /// or the loop is not OK, the return value will be false
   bool orientationOK()
	{
	    if (loop_fixed_.size() != boundary_loops_.size())
		return false;
	    for (size_t ki=0; ki<loop_fixed_.size(); ++ki)
		if (loop_fixed_[ki])
		    return false;
	    return true;
	}

   /// Remove history of loop orientation fixes
    void setOrientationOK()
	{
	    loop_fixed_.resize(boundary_loops_.size());
	    std::fill(loop_fixed_.begin(), loop_fixed_.end(), 0);
	}

    /// Turn orientation of specified loop, and remember turning it
    void turnLoopOrientation(int idx);

    /// The boundary loops may be outside the loop tolerance, or the
    /// loops do not fullfill the loop requirements (1 outer loop
    /// which is ccw, cw loops inside, loops should be simple and
    /// disjoint).
    bool isValid(int& valid_state) const;

    /// We try to fix the invalid loops. Return value: true if the
    /// loops are valid. Returned max_gap is the largest gap between
    /// end segments in the loops.  Assuming that max_tol_mult >=
    /// 1.0. The routine is allowed to alter the eps to
    /// eps*max_tol_mult.
    bool fixInvalidSurface(double& max_loop_gap, double max_tol_mult = 1.0);

    /// Checking all boundary_loops_ to see it they fulfill
    /// requirements. Routine sets valid_state_. Must be called after
    /// surface is altered.
    void analyzeLoops();

    /// If both parameter and space curve are given for a segment, and
    /// they do not match, one of them is removed, unless we may alter
    /// the tolerance slightly (at most to epsgeo*max_tol_mult).
    void removeMismatchCurves(double max_tol_mult);

    /// We measure the largest distance from loop to the surface. If
    /// the loops are defined by curves in the parameter domain, then it
    /// is trivially 0.0.
    /// Useful for testing whether the tolerance makes any sense.
    double maxLoopSfDist(int loop_ind, int nmb_seg_samples = 100);

    /// We measure the largest distance between loop segments.
    double maxLoopGap();

    /// Given a parameter value corresponding to on specified curve in
    /// a specified boundary loop, return the corresponding surface
    /// parameter
    Point getSurfaceParameter(int loop_idx, int cv_idx, double bd_par) const;

    /// Simplify boundary loops by reducing the number of curves if possible
    bool simplifyBdLoops(double tol, double ang_tol, double& max_dist);

    /// Test if underlying surface is spline, if not try to convert it to spline
    /// /return wether underlying surface is spline afterwards
    bool makeUnderlyingSpline();

    /// Test if underlying surface and all loop curves are CurveOnSurface
    /// where the space curves are SplineCurve.
    bool allIsSpline() const;

    /// Return copy where underlying surface and all loop curves are CurveOnSurface
    /// where the space curves are SplineCurve.
    BoundedSurface* allSplineCopy() const;

    // This surface is axis rotational if the underlying surface is and
    // all the trimming curves are consistent with being rotational
    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle);

    /// This surface is planar if the underlying surface is
    virtual bool isPlanar(Point& normal, double tol);

    /// Check if the surface is linear in one or both parameter directions
    virtual bool isLinear(Point& dir1, Point& dir2, double tol);

    /// Run through the boundary loops, returning the smallest epsgeo.
    double getEpsGeo() const;

    // Check if any of the input params are close to the end parameters of surface_.
    // True if inside (umax - min)*domain_fraction of an end parameter.
    bool closeToUnderlyingBoundary(double upar, double vpar,
				   double domain_fraction = 0.01) const;

  /// Check if a polynomial element (for spline surfaces) intersects the
    /// (trimming) boundaries of this surface
    /// \param elem_ix: Element index counted according to distinct knot
    /// values. Sequence of coordinates: x runs fastest, then y
    /// \param eps: Intersection tolerance
    /// \return -1: Not a spline surface or element index out of range
    ///          0: Not on boundary or touching a boundary curve
    ///          1: On boundary (intersection with boundary found)
    /// Note that a touch with the boundaries of the underlying surfaces
    /// is not consdered a boundary intersection while touching a trimming
    /// curve is seen as an intersection
    virtual int ElementOnBoundary(int elem_ix, double eps);

   /// Check if a polynomial element (for spline surfaces) intersects the
    /// (trimming) boundaries of this ftSurface, is inside or outside
    /// \param elem_ix: Element index counted according to distinct knot
    /// values. Sequence of coordinates: x runs fastest, then y
    /// \param eps: Intersection tolerance
    /// \return -1: Not a spline surface or element index out of range
    ///          0: Outside trimmed volume
    ///          1: On boundary (intersection with boundary found)
    ///          2: Internal to trimmed surfaces
    /// Note that a touch with the boundaries of the underlying surface
    /// is not consdered a boundary intersection while touching a trimming
    /// curve is seen as an intersection
    virtual int ElementBoundaryStatus(int elem_ix, double eps);

    friend void 
      GeometryTools::setParameterDomain(std::vector<shared_ptr<BoundedSurface> >& sfs,
				       double u1, double u2, 
				       double v1, double v2);
private:
    /// The underlying surface
    shared_ptr<ParamSurface> surface_;

    /// The curves describing the boundaries of the surface.  First
    /// element is the outer boundary loop (ordering is done by
    /// constructor).
    std::vector<shared_ptr<CurveLoop> > boundary_loops_;

    /// Indicates if the boundary loop has been fixed with respect to
    /// orientation
    std::vector<int> loop_fixed_;

    mutable CurveBoundedDomain domain_;

    mutable bool iso_trim_;
    mutable double iso_trim_tol_;

    mutable BoundingBox box_;

    // The trim curves should be valid loops. Additionally the first
    // element should be the outer ccw loop, all other loops should be
    // cw loops lying inside the ccw loop.
    int valid_state_; //  0 = not validated (analyze not performed / failed).
                      //  1 = valid.
                      // -1 = par & space cv mismatch.
                      // -2 = par cv(s) missing (required).
                      // -4 = loop(s) not closed (dir of segments, order, gaps).
                      // -8 = loops not ordered or direction wrong.
                      // -15 = -1 -2 -4 -8, i.e. all artifacts/features.

    /// Helper function. When called with analyze = true no fixing is
    /// performed. Otherwise (i.e. if analyze = false) the routine
    /// tries to fix gap(s). If unsuccessful, nothing is changed.
    bool fixLoopGaps(double& max_loop_gap, bool analyze);

    /// Helper functions Order boundary_loops_ w/outer boundary loop
    /// first. When called with analyze=true the loops are not
    /// rearranged (even if they were found to be invalid). Otherwise
    /// the routine tries to fix loops. Expecting consistent loops
    /// (closed and simple). Routine also fixes the direction of the
    /// loops if not consistent (the outer should be ccw and lie first
    /// in loop vector, all other cw).
    bool orderBoundaryLoops(bool analyze,
			    double degenerate_epsilon = DEFAULT_SPACE_EPSILON);

    /// We then look for missing par cvs.
    bool parameterCurveMissing();

    // We see if the par cv and the space cv match (i.e. if trace
    // is the same, as well as direction).
    // By calling this routine the eps tmay be increased to
    // max_tol_mult*eps. Assuming that max_tol_mult >= 1.0.
    // @@sbr072009 Mismatch in domain and orientation should be handled
    // from this class. Missing par cv and mismatch between par &
    // space cv should be handled from the outside as it requires more
    // machinery. To be implemented!
    bool fixParSpaceMismatch(bool analyze, double max_tol_mult,
			     int nmb_seg_samples);

    // We want the boundary curve to be at least c1. To be called from
    // public function.
    std::vector<shared_ptr<CurveOnSurface> >
    splitIntoC1Curves(shared_ptr<CurveOnSurface>& curve,
		      double space_epsilon, double kink);

    // Used to avoid code duplication in two nearly equal
    // constructors.
    void
    constructor_implementation(shared_ptr<ParamSurface> surf,
			       std::vector<std::vector<shared_ptr<CurveOnSurface> > >
			       loops,
			       std::vector<double> space_epsilons,
			       bool fix_trim_cvs);

    bool checkParCrvsAtSeam();

    void setParameterDomainBdLoops(double u1, double u2, 
				   double v1, double v2);
};


} // namespace Go


#endif // _BOUNDEDSURFACE_H

