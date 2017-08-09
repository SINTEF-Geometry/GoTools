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

#ifndef _INTERSECTIONPOOL_H
#define _INTERSECTIONPOOL_H


#include "GoTools/intersections/ParamObjectInt.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/utils/Point.h"
#include <memory>
#include <vector>
#include <set>
#include <ostream>


namespace Go {


class IntersectionPool;
struct BoundaryIntersectionData;


/// Class providing access to the intersection points and surfaces
/// related to a particular subproblem.

class IntersectionPool {
public:

    /// Constructor
    /// \param obj1 the first object of the (possible) intersection
    /// \param obj2 the second object of the (possible) intersection
    /// \param parent pointer to parent pool (null pointer if no
    /// parent)
    /// \param missing_dir index of the missing parameter (negative if
    /// no missing parameter)
    /// \param missing_value value of the missing parameter
    IntersectionPool(shared_ptr<ParamObjectInt> obj1, 
		     shared_ptr<ParamObjectInt> obj2,
		     shared_ptr<IntersectionPool> parent
		     = shared_ptr<IntersectionPool>(),
		     int missing_dir = -1,
		     double missing_value = 0);

    /// Destructor
    virtual ~IntersectionPool() {}

    /// Add a number of IntersectionPoint s to the pool by specifying
    /// their parameter values.
    /// \param nmb_int_pts the number of IntersectionPoint s that
    /// shall be added to the pool.
    /// \param pointpar1 pointer to an array containing the
    /// consecutive parameter values in the first object for the
    /// IntersectionPoint s to add.
    /// \param pointpar2 pointer to an array containing the
    /// consecutive parameter values in the second object for the
    /// IntersectionPoint s to add.
    /// \param epsge shared pointer to the object specifying the
    /// tolerances to use for the added IntersectionPoint s.
    void addParamPoints(int nmb_int_pts, 
			double* pointpar1, 
			double* pointpar2,
			shared_ptr<GeoTol> epsge);

    /// Add a number of IntersectionCurves to the pool by specifying
    /// iterators to the IntersectionPoint making up the curves.
    /// \param nmb_int_cvs number of curves to add
    /// \param start_iterators pointer to a range of iterators to
    /// IntersectionPoint.  The range should have length equal to the
    /// number of curves to add.  Each iterator represents the start
    /// of a range of IntersectionPoints that defines an
    /// IntersectionCurve to add to the pool. (The corresponding end
    /// of the range is specified by the corresponding iterator in \a
    /// end_iterators.
    /// \param end_iterators pointer to a range of iterators to
    /// IntersectionPoint.  The range should have length equal to the
    /// number of curves to add.  Each iterator represents the
    /// one-past-end of a range of IntersectionPoint s that defines an
    /// IntersectionCurve to add to the pool. (The corresponding start
    /// of the range is specified by the corresponding iterator in \a
    /// start_iterators)
    template<class ip_iterator>
    void addCurves(int nmb_int_cvs, 
		   ip_iterator* start_iterators, 
		   ip_iterator* end_iterators);

    /// Add one IntersectionCurve to the pool, by specifying a range
    /// of IntersectionPoints that make up the curve.
    /// \param start iterator to the start of the range of
    /// IntersectionPoints.
    /// \param end iterator to the end of the range of
    /// IntersectionPoints.
    template<class ip_iterator>
    void addCurve(ip_iterator start, ip_iterator end);

    /// Get a vector containing shared pointers to all
    /// IntersectionPoints in the pool that is located on the
    /// boundary of one of the pool's objects.
    /// \retval bd_ints vector of shared pointers to boundary
    /// intersection points
    void getBoundaryIntersections(std::vector<shared_ptr<IntersectionPoint> >& bd_ints);

    // @@@ VSK, Change name when implementing.
    // This function is not currently in use.
    bool  checkIfBothPointsLieOnOneEndpointAndNotOfTheSameCurve();


    /// Get a vector containing shared pointers to all the
    /// IntersectionPoints in the pool.  The entries in the vector
    /// will be sorted according to parameter values.
    /// \retval int_pts vector containing the pool's
    /// IntersectionPoints in a sorted manned.
    void getSortedIntersections(std::vector<shared_ptr<IntersectionPoint> >& int_pts);

    /// Get a vector containing the parametric values along a specific
    /// parameter direction for all the pool's IntersectionPoints
    /// that do not lie 'on' or 'close to' the boundary along this
    /// parameter.  The returned values will be sorted according to
    /// increasing value.
    /// \param pardir the concerned parameter direction
    /// \return the vector containing the requested parameter values.
    std::vector<double> getSortedInnerInts(int pardir);

    /// Get a vector containing shared pointers to the pool's
    /// IntersectionPoints, sorted along a given spatial or
    /// parametrical direction. NOTE: The name of this function
    /// contains 'Bd' for 'boundary' for historical reasons. In fact,
    /// we get \em all intersections.
    /// \param vec the direction along which to sort the
    /// IntersectionPoint.  If the number of vector elements is 3, it
    /// is interpreted as a 'spatial' vector; otherwise it is
    /// interpreted as a vector in the combined parametric domain of
    /// the two objects.
    /// \retval result vector of shared pointers to
    /// IntersectionPoints
    void getSortedBdInts(const Point& vec,
			 std::vector<shared_ptr<IntersectionPoint> >& result,
			 int obj_nmb = -1);

    /// Check the consistency between the direction of the intersection
    /// curve indicated by vec and the direction computed in the
    /// existing intersection points. Turn the direction of vec if
    /// necessary.
    /// \param vec the direction of the intersection curve. This mey
    /// be turned if necessary.
    /// \param sorting_obj the object - 0 or 1 - used for sorting
    /// \retval 'true' if ???, 'false' if ???
    // @@@jbt - NB: Find out about the return value...
    bool checkSortingDir(Point& vec, int sorting_obj);

    // This function is not yet implemented
    void getSortedInts(/* AROUND_BOUNDARIES ,*/
		       std::vector<shared_ptr<IntersectionPoint> >& result)
    { }

    // This function is not yet implemented
    void getAllIntersections(std::vector<shared_ptr<IntersectionPoint> >& result)
    { }

    /// Get a vector containing all those IntersectionPoints in the
    /// pool that are branchpoints.
    /// \retval pts vector of shared pointers to IntersectionPointss
    void getBranchPoints(std::vector<shared_ptr<IntersectionPoint> > & pts);

    /// Compute the middle parameter value of all intersection points
    /// of this intersection pool
    /// \retval mid points to an array of \c doubles where the
    /// result will be written.  The number of elements in the array
    /// is equal to the total number of parameter directions for this
    /// pool.
    void getMidParameter(double *mid);

    // @@@ VSK. Let addIntersection point return the new intersection
    // point. It is sometimes useful.

    /// Construct a new IntersectionPoint and add it to this pool, as
    /// well as the pool's parents (and older ancestors).  NB: The
    /// objects given as input should be equal to - or sub-objects of
    /// the objects that are already pointed to by this pool.
    /// \param obj_int1_ The first of the intersecting objects of the
    /// new IntersectionPoint
    /// \param obj_int2_ The second of the intersecting objects of the
    /// new IntersectionPoint
    /// \param epsge the tolerances used for the new IntesectionPoint
    /// \param par1 pointer to the IntersectionPoint's parameters in
    /// the first object
    /// \param par2 pointer to the IntersectionPoint's parameters in
    /// the second object
    shared_ptr<IntersectionPoint>
    addIntersectionPoint(shared_ptr<ParamObjectInt> obj_int1_, 
			 shared_ptr<ParamObjectInt> obj_int2_,
			 shared_ptr<GeoTol> epsge,
			 double *par1, double *par2);

    /// Include the IntersectionPoints in the lower_order_pool into
    /// 'this' pool, generating new representations of the points
    /// which includes the missing parameter.  It is assumed that the
    /// IntersectionPoints in \a lower_order_pool have no
    /// IntersectionLink with points outside the pool.
    /// \param lower_order_pool a pointer to the lower-order pool from
    /// which we want to include its IntersectionPoints.
    void includeReducedInts(shared_ptr<IntersectionPool>
			    lower_order_pool);

    /// Reorganize self-intersection parameters.  A 'twin point' is a
    /// concept when working with self-intersection objects.  It
    /// represent an IntersectionPoint that already exist in the
    /// parent pool, but whose ordering of parameters had to be
    /// switched in order to conform with its objects (which are
    /// really two parts of the same, global surface).  If a pool
    /// contains twin points, we must assure that the curves where
    /// these points are included have a consistent ordering of the
    /// parameters.  This function takes care of making this ordering
    /// consistent.  Moreover, since the twin points are really just
    /// another representation of an already-existing point, they will
    /// go out of scope and disappear when the sub-pool in which they
    /// lie is destroyed.  Therefore, we must make sure that points
    /// linked to by the twin points will now be linked to by the
    /// original points that the twin-points are representing.
    /// 
    /// \b NB: This is a very special function that is only used with
    /// the following preconditions:
    /// \li this pool should represent a self-intersection, ie. its
    /// first and second object should be the same.
    /// \li the 'sub_pool' (given as argument) should be a child of
    /// 'this' pool.
    /// \li 'this' pool and the one pointed to by \a sub_pool should
    /// have the same number of parameters.
    ///
    /// \param sub_pool shared pointer to the sub-pool, as described
    /// above.
    void 
    selfIntersectParamReorganise(shared_ptr<IntersectionPool> sub_pool);

    // This function is to be used only in self intersection context, when
    // boundary intersections are computed.
    void mirrorIntersectionPoints(size_t first);


    /// Check if this IntersectionPool contains IntersectionPoints
    /// that are not "boundary points" with respect to the parameter
    /// \a pardir.
    /// \param pardir the parameter direction we want to check for
    /// \return 'true' if there was found at least one
    /// IntersectionPoint whose \a pardir parameter did not lie on the
    /// extremal values of its interval.  'false' otherwise.
    bool hasPointsInInner(int pardir);

    /// Check if this IntersectionPool contains any
    /// IntersectionPoints at all.
    /// \return 'true' if 'this' pool contains at least one
    /// IntersectionPoint.
    bool hasIntersectionPoints() const
    { return (int_points_.size() > 0); }
    
    /// Query the number of IntersectionPoints contained in this
    /// IntersectionPool.
    /// \return the number of IntersectionPoints in this pool.
    int numIntersectionPoints()	const
    { return (int)int_points_.size(); }

    /// \brief Get the number of IntersectionPoints in this pool who
    /// are not classified as ORDINARY_POINTs (see \ref
    /// SingularityType for more info), i.e. points that lie on
    /// non-transversal intersections.
    ///
    /// \return number of singular IntersectionPoints.
    int numSingularIntersectionPoints();

    /// Check if a given IntersectionPoint is in this
    /// IntersectionPool.
    /// \param pnt pointer to the IntersectionPoint that we want to
    /// check.
    /// \return 'true' if this IntersectionPoint was found in this
    /// IntersectionPool, 'false' otherwise.
    bool isInDomain(IntersectionPoint *pnt) const;

    /// Locate the IntersectionPoint in this IntersectionPool whose
    /// parameter values are closest (in the Euclidean norm on the
    /// parameter space) to a given set of parameters.
    /// \param param[] pointer to an array containing the parameter
    /// values that we want to find the closest IntersectionPoint to.
    /// \retval pnt upon function return, if a closest
    /// IntersectionPoint was found, \a pnt will be set to that
    /// IntersectionPoint.
    /// \return 'true' if a closest IntersectionPoint was found,
    /// 'false' otherwise.
    bool closestInDomain(double param[],
			 shared_ptr<IntersectionPoint>& pnt) const;

    /// Locate the IntersectionPoint in this IntersectionPool whose
    /// parameter values are closest (in the Euclidean norm on the
    /// parameter space) to a given set of parameters.
    /// \param param[] pointer to an array containing the parameter
    /// values that we want to find the closest IntersectionPoint to.
    /// \retval pnt upon function return, if a closest
    /// IntersectionPoint was found, \a pnt will be set to that
    /// IntersectionPoint.
    /// \return 'true' if a closest IntersectionPoint was found,
    /// 'false' otherwise.
    bool closestInDomain(const double param[],
			 shared_ptr<IntersectionPoint>& pnt) const;

    /// Checks the existence of any intersection points with the given
    /// parameter in the given direction.
    /// \param dir specifies the concerned parameter direction.
    /// \param par specifies the value of this parameter.
    /// \return 'true' if the pool contains an IntersectionPoint whose
    /// parameter \a dir is approximately equal to \a par.
    bool existIntersectionPoint(int dir, double par);

    /// Check if a given parameter value in a given parameter
    /// direction lies in the influence area of any
    /// IntersectionPoints in the pool.
    /// \param pardir the concerned parameter direction
    /// \param par the parameter value to check against
    /// \param first_outside if this is set to 'true' then the outer
    /// brackets of the influence intervals will be used (which in
    /// practice means that if the function returns 'false', then \a
    /// par is \em guaranteed to be outside the influence area of all
    /// points).  If this argument is set to 'false', then the inner
    /// brackets of the influence interval will be used (which in
    /// practice means that if the function returns 'true', then \a
    /// par is \em guaranteed to be inside the influence area of at
    /// least one point).
    /// \return
    /// \li 0 if no IntersectionPoint with an influence area covering
    /// the specified value was found.
    /// \li 1 if we found an IntersectionPoint whose influence area
    /// covered the specified value, but where the specified value did
    /// not hit \em exactly on the IntersectionPoint's own parameter.
    /// \li 2 if we found an IntersectionPoint whose parameter value
    /// in the specified direction coincide with the specified value.
    int inInfluenceArea(int pardir, double par, bool first_outside = true);

    /// Check if a given parameter value in a given parameter
    /// direction lies in the influence area of any
    /// IntersectionPoints in the specified vector.
    /// \param pardir the concerned parameter direction
    /// \param par the parameter value to check against
    /// \param int_pts the vector of IntersectionPoints to check
    /// against.
    /// \param first_outside if this is set to 'true' then the outer
    /// brackets of the influence intervals will be used (which in
    /// practice means that if the function returns 'false', then \a
    /// par is \em guaranteed to be outside the influence area of all
    /// points).  If this argument is set to 'false', then the inner
    /// brackets of the influence interval will be used (which in
    /// practice means that if the function returns 'true', then \a
    /// par is \em guaranteed to be inside the influence area of at
    /// least one point).
    /// \return
    /// \li 0 if no IntersectionPoint with an influence area covering
    /// the specified value was found.
    /// \li 1 if we found an IntersectionPoint whose influence area
    /// covered the specified value, but where the specified value did
    /// not hit \em exactly on the IntersectionPoint's own parameter.
    /// \li 2 if we found an IntersectionPoint whose parameter value
    /// in the specified direction coincide with the specified value.
    int inInfluenceArea(int pardir, 
			double par, 
			std::vector<shared_ptr<IntersectionPoint> >& int_pts, 
			bool first_outside = true);

    /// Fill the vector given as argument with shared pointers to all
    /// the IntersectionPoints in 'this' IntersectionPool.
    /// \retval int_pts vector of shared pointers to all of the pool's
    /// IntersectionPoints.
    void getIntersectionPoints(std::vector<shared_ptr<IntersectionPoint> >& int_pts) const;
    
    /// Fill the argument vector with shared pointers to all the
    /// IntersectionCurves in 'this' IntersectionPool.
    /// \param int_curves vector of shared pointers to all of the
    /// pool's IntersectionCurves.
    void getIntersectionCurves(std::vector<shared_ptr<IntersectionCurve> >&
			       int_curves) const;
    
    /// Get a reference to the vector containing (shared pointers to)
    /// 'this' IntersectionPool's IntersectionPoints.
    /// \retval reference to vector of shared pointers to
    /// IntersectionPoints
    std::vector<shared_ptr<IntersectionPoint> >&
    getIntersectionPoints();

    /// Get a reference to the vector containing (shared pointers to)
    /// 'this' IntersectionPools IntersectionPoint s.
    /// \retval reference to vector of shared pointers to
    /// IntersectionPoints
    const std::vector<shared_ptr<IntersectionPoint> >&
    getIntersectionPoints() const;

    /// Get a reference to the vector containing (shared pointers to)
    /// 'this' IntersectionPools IntersectionCurves.
    /// \retval reference to vector of shared pointers to
    /// IntersectionCurves
    const std::vector<shared_ptr<IntersectionCurve> >&
    getIntersectionCurves() const;

    /// Fetch the IntersectionPoints with a specified parameter value
    /// in a specified direction
    /// \param dir the specified parameter direction
    /// \param par the specified parameter value
    /// \retval result the vector containing those of the pool's
    /// IntersectionPoints that (approximately) had the specified
    /// parameter value in the specified parameter direction.
    void getIntersectionPoints(int dir, double par,
			       std::vector<shared_ptr<IntersectionPoint> >& result) const;

    /// Check if the IntersectionPool contains at least one
    /// IntersectionPoint with the specified parameter value in the
    /// specified direction.
    /// \param dir the specified parameter direction
    /// \param par the specified parameter value
    /// \return 'true' if at least one IntersectionPoint with the
    /// specified parameter value in the specified parameter direction
    /// was found in the pool.  'false' otherwise.
    bool hasIntersectionPoints(int dir, double par) const;

    /// Check if the intersection points are OK. If they are not, give
    /// further diagnostics.
    /// \param int_pt_info vector with information about the
    /// intersection points. This information has the type of
    /// struct IntPtInfo, which indicates if the point is ok, the number
    /// of neighbours, the singularity type, and the location of the
    /// point is the parameter domains of the objects
    /// \return 'true' if all intersection points are OK, 'false' if
    /// at least one point is not OK
    bool checkIntersectionPoints(std::vector<IntPtInfo>& int_pt_info) const;

    /// Check if the intersection links are OK. If they are not, give
    /// further diagnostics.
    bool checkIntersectionLinks() const;

    /// Remove double points by merging into one single point. Double
    /// points are points that are the same within the tolerance. The
    /// new points that results from this will inherit all the links
    /// from the old points.
    void removeDoublePoints();

    // Remove ordinary intersection points with only one neighbour when these
    // points lie reasonably close, their tangents point roughly in the same
    // direction, and the neighbour has significantly better accuracy.
    void removeLooseEnds();

    // Check ordinary intersection points to check if the links they belong to
    // are consistent
    void removeFalseCorners();

    /// Remove defect intersection links. Whether or not a link is
    /// defect depends on the return value of
    /// verifyIntersectionLink().
    void removeDefectLinks();

    void getAllLinks(std::set<shared_ptr<IntersectionLink> >& links);

    /// Get all isolated IntersectionPoints in the pool (those
    /// with no neighbours), as well as all IntersectionCurves.
    /// \retval int_points upon function return, this vector will
    /// contain shared pointers to all the IntersectionPoints in
    /// this pool that had no neighbours.
    /// \retval int_curves upon function return, this vector will
    /// contain shared pointers to all the IntersectionCurves
    /// contained in this pool.
    void getResult(std::vector<shared_ptr<IntersectionPoint> >& int_points, 
		   std::vector<shared_ptr<IntersectionCurve> >& int_curves);

    /// Get the index of the missing parameter (-1 if none missing)
    /// \return the index of the missing parameter
    int lackingParameter() const;

    /// Get value of missing parameter (if it exists)
    /// \return value of the missing parameter
    double lackingParameterValue() const;

    /// Check if two IntersectionPoints lie on the same boundary of
    /// the parametric domain of this IntersectionPools objects
    /// (\a obj1_ and \a obj2_).
    /// \param pt1 the first IntersectionPoint
    /// \param pt2 the second IntersectionPoint
    /// \return 'true' if \a pt1 and \a pt2 both were found to lie on
    /// the same (parametric) boundary.  'false' otherwise.
    bool atSameBoundary(shared_ptr<IntersectionPoint> pt1,
			shared_ptr<IntersectionPoint> pt2);

    /// Check if two IntersectionPoints can be found on two
    /// different boundaries of the parametric domain of this
    /// IntersectionPools objects (\a obj1_ and \a obj2_)
    /// \param pt1 the first IntersectionPoint
    /// \param pt2 the second IntersectionPoint
    /// \return 'true' if \a pt1 and \a pt2 could be found to lie on
    /// two distinct borders (i.e. \a pt1 lay on one border and \a pt2
    /// lay on another).
    bool atDifferentBoundary(shared_ptr<IntersectionPoint> pt1,
			     shared_ptr<IntersectionPoint> pt2);

    /// Check if the two specified IntersectionPoints represent a
    /// boundary intersection in this pool.  In order to be
    /// interpreted as such, they must both lie on the same
    /// (parametric) boundary, and be connected with an iso-parametric
    /// IntersectionLink
    /// \param pt1 the first IntersectionPoint
    /// \param pt2 the second IntersectionPoint
    /// \return 'true' if \a pt1 and \a pt2 satisfy the criteria
    /// specified above, 'false' otherwise.
    bool isBoundaryIntersection(shared_ptr<IntersectionPoint> pt1,
				shared_ptr<IntersectionPoint> pt2);

    /// Check if a point lies at a boundary in the current domain(s)
    bool isBoundaryPoint(shared_ptr<IntersectionPoint> pt1)
	{
	    return isBoundaryPoint(pt1.get());
	}

    bool isBoundaryPoint(IntersectionPoint* pt1);

    // @@@ VSK. Yet some new functions. I just make a dummy implementation

    /// Get all the intersection loops contained in this
    /// IntersectionPool.  This function assumes that the topologies
    /// of the loops in this pool are not too 'warped'.  Notably, it
    /// might fail if two loops share some edges, since in that case
    /// there will be ambiguous choices for what constitutes the
    /// loops.  (This is a fact from graph theory, and if we want to
    /// work around it, we might have to consider additional
    /// information of a non-topological nature, like looking at the
    /// geometric/parametric position of the points).
    /// \retval loop_ints upon function return, this vector will
    /// contain all the loops that were found in this IntersectionPool
    /// (each loop is represented as a vector of
    /// IntersectionPoints).
    /// \return 'true' if at least one loop was found. 'false'
    /// otherwise.
    bool fetchLoops(std::vector<std::vector<shared_ptr<IntersectionPoint> > >& loop_ints);

    /// Check if a given loop delimits a partial coincidence area
    /// (PAC).  To determine this the meta-information in the
    /// IntersectionLink s is examined.
    /// \param int_loop the loop to check
    /// \return 'true' if the loop delimits a PAC, 'false' otherwise.
    bool isPAC(std::vector<shared_ptr<IntersectionPoint> >& int_loop);

    /// Set the specified loop to be a partial coincidence area (PAC).
    /// This is written into the loop's IntersectionLinks as
    /// meta-information.
    /// \param int_loop the loop in question
    void setCoincidence(std::vector<shared_ptr<IntersectionPoint> >& int_loop);

    // Check if there is a connection between the given points through
    // the points in the intersection pool. For the time being

    /// Checks if there is a path of connected IntersectionPoints
    /// going between the two IntersectionPoints specified.  The
    /// specified points are assumed to be in this IntersectionPool,
    /// and all IntersectionPoints of the path must also be in the
    /// pool.
    /// \param pt1 the first IntersectionPoint.  Must already be in
    /// the pool.  Otherwise, an error will occur.
    /// \param pt2 the second IntersectionPoint.  Must already be in
    /// the pool.  Otherwise, an error will occur.
    /// \return 'true' if a path was found connecting \a pt1 with \a
    /// pt2 using only IntersectionPoints in the pool.  'false'
    /// otherwise.
    bool isConnectedInside(shared_ptr<IntersectionPoint> pt1,
			   shared_ptr<IntersectionPoint> pt2);

    /// Prepare output by making IntersectionCurves. Points that
    /// have one or three or more IntersectionLinks are considered
    /// endpoints to curves, points with two IntersectionLinks lie
    /// in the inner of a curve.
    void makeIntersectionCurves();

    /// Get the "original" intersection points of this intersection
    /// pool. If this pool has no parents with the same number of
    /// parameters as itself, then return the IntersectionPoints of
    /// this pool.  Otherwise, return the IntersectionPoints of the
    /// parent pool.
    /// \retval int_pts the returned IntersectionPoints.
    void getOrigPoints(std::vector<shared_ptr<IntersectionPoint> >& int_pts) const;

    /// Write various  debug information
    /// \param singular use 1 for singular case, 0 is default 
    void writeDebug(int singular = 0);

    /// Among all IntersectionLinks between IntersectionPoints in
    /// this pool, split those who cross a specified parameter
    /// direction at a specified parameter value.  Split means to
    /// insert a new IntersectionPoint at this place.
    /// \param fixed_dir the specified parameter direction
    /// \param fixed_val the specified parameter value
    void splitIntersectionLinks(int fixed_dir, double fixed_val);

    /// Get the combined number of parameter directions in the two
    /// objects.
    /// \return the number of parameters in the first object plus the
    /// number of parameters in the second object.
    int numParams() const
    { return obj1_->numParams() + obj2_->numParams(); }

    /// Get the start of the parameter interval for the specified
    /// parameter direction
    /// \param dir the specified parameter direction
    /// \return the start of the parameter interval for the \a dir
    /// parameter direction.
    double startParam(int dir) const;

    /// Get the end of the parameter interval for the specified
    /// parameter direction
    /// \param dir the specified parameter direction
    /// \return the end of the parameter interval for the specified
    /// parameter direction.
    double endParam(int dir) const;

    /// Check all the neighbours of the IntersectionPoints in this
    /// pool, and add those of them that are not already in the pool,
    /// but whose parameters are covered by the pool's parameter
    /// domain.
    void includeCoveredNeighbourPoints();

    /// Create a new IntersectionPoint and add it to the the
    /// IntersectionPool.  The new IntersectionPoint (A) will lie
    /// inside the influence interval of another of the pool's points;
    /// the one (B) indicated by \a pt_in_pool (supposedly in the pool
    /// already).  The parameter direction in which (A) lies inside
    /// the influence interval of (B) is indicated by \a pardir.  The
    /// parameter values of the new point is pointed to by \a parvals
    /// (an array which also contains the parameter value in \a
    /// pardir.  There will be inserted an IntersectionLink between
    /// (A) and (B), which will be isoparametric in the other
    /// parameter direction if \a pardir is on a two-parametric object
    /// (surface or function).
    /// \param pt_in_pool the IntersectionPoint (supposedly already in
    /// the pool) whose influence area the new IntersectionPoint is
    /// going to lie.
    /// \param parvals pointer to an array specifying the parameter
    /// values for the new IntersectionPoint to create.
    /// \param pardir the parameter direction that will be
    /// iso-parametric in the newly established IntersectionLink (see
    /// above).
    void 
    insertInInfluenceInterval(shared_ptr<IntersectionPoint> pt_in_pool,
                              double *parvals, int pardir);

    /// Synchronize pool with its parent. We check that the pool is
    /// "valid" by making sure that the intersection points on the
    /// current level is also present at previous levels. If there are
    /// intersection points that are redundant in this sense, we
    /// remove them. This function us useful when several sibling
    /// subintersectors are around, and running compute on one of them
    /// has removed intersection points that are also present in
    /// someof the others.
    void synchronizePool();

    /// Remove redundant IntersectionPoints from the pool, starting
    /// from the point indexed \a first_idx and upwards.  An
    /// IntersectionPoint is considered redundant if:
    /// \li it has exactly two neighbours \em AND
    /// \li it lies on an isoparametric intersection curve \em AND
    /// \li both of its (two) neighbours are in the pool
    ///
    /// \param first_idx index of the first point
    void cleanUpPool(int first_idx = 0, double epsge=1.0e-15);

    /// Remove an IntersectionPoint from the pool. Should only be
    /// called if the number of neighbours are less than or equal to
    /// two.

    /// \param int_point the IntersectionPoint to remove from the pool
    /// and its ancestors (as long as the ancestors have the same
    /// number of parameter directions).  IntersectionLinks will be
    /// broken up and the removed IntersectionPoint's neighbours will
    /// be reconnected.
    /// \param int_point shared pointer to the IntersectionPoint to
    /// remove.
    void removeIntPoint(shared_ptr<IntersectionPoint> int_point);
    
    /// Remove an IntersectionPoint from the pool.
    /// \param int_point the IntersectionPoint to remove from the pool
    /// and its ancestors (as long as the ancestors have the same
    /// number of parameter directions).  IntersectionLinks will be
    /// broken up.
    /// \param int_point shared pointer to the IntersectionPoint to
    /// remove.
    void removeIntPoint2(shared_ptr<IntersectionPoint> int_point);
    
    /// Remove all intersection points from the pool whose parameters
    /// lie in the box defined by \a frompar and \a
    /// topar. IntersectionLinks will be broken up.
    /// \param frompar the array of parameters defining the lower
    /// corner of the box in question.
    /// \param topar the array of parameters defining the upper
    /// corner of the box in question.
    void removeIntPoints(double* frompar, double* topar, bool only_inner = false);
    
    /// \brief Check if two surfaces intersect along a common boundary
    /// (i.e. a linked intersection list exist along one boundary in
    /// both objects).
    ///
    /// Intersections whose parametric extent is smaller than a
    /// specified fraction (\a frac) of the surface's parameter
    /// intervals are ignored.  The data on the found boundary
    /// intersections are returned as entries into a vector of
    /// BoundaryIntersectionData.
    /// \param frac the fraction of a parameter's total span that must
    /// be covered in order to consider this intersection "long
    /// enough" to be considered. (A value between 0 and 1, typically
    /// small).
    /// \retval isects vector with objects containing information on
    /// the found boundary intersections.
    void intersectAlongCommonBoundary(double frac,
				      std::vector<BoundaryIntersectionData>&
				      isects);

    /// Removes IntersectionLinks positioned along the boundary of
    /// both objects, The concerned IntersectionPoints are also
    /// removed, except when they also links to points which are NOT
    /// on the boundary.
    void removeBoundaryIntersections(bool single_pt=true);

    /// Verify that the given intersection link is in the pool and
    /// represents a connected piece of the intersection curve.
    /// \param link shared pointer to the intersection link we wish to
    /// verify
    /// \retval 'true' if the link represents a connected piece of the
    /// intersection curve, 'false' otherwise
    bool
    verifyIntersectionLink(const shared_ptr<IntersectionLink>& link,
			   int recursion_limit = 10)
	const;

    /// Check consistence of chain of intersection links with respect to
    /// number of neighbours
    bool 
	checkIntersectionChain(IntersectionPoint *pnt,
			       IntersectionPoint *first, 
			       IntersectionPoint *prev = 0);

    /// Checks if the pool is "valid", that is, if the intersection
    /// points on the present subdivision level exist also on previous
    /// levels.
    /// \return \a true if the pool is valid, \a false otherwise
    bool validate() const;

    // Remove unnecessary or harmful intersection points (called
    // before making IntersectionCurves.
    void weedOutClutterPoints() ; 

    /// Writes out a list of the current intersection points in the
    /// pool. Writes to cout.
    void writeIntersectionPoints() const;

    /// Writes out a list of links beween intersection points in the
    /// pool. Writes to cout.
    void writeIntersectionLinks() const;

    friend class SfSfIntersector;

private:
    // Data members

    // The intersecting objects
    shared_ptr<ParamObjectInt> obj1_;
    shared_ptr<ParamObjectInt> obj2_;

    // Intersection points (both isolated and on curves)
    std::vector<shared_ptr<IntersectionPoint> > int_points_;

    // Intersection curves
    std::vector<shared_ptr<IntersectionCurve> > int_curves_;

    int missing_param_index_;
    double missing_param_value_;

    shared_ptr<IntersectionPool> prev_pool_;
//     static const double REL_PAR_RES_; // a tolerance - @move this
// 				      // somewhere else??

    // Functions

    void 
    add_point_and_propagate_upwards(shared_ptr<IntersectionPoint> point);

    void iterateToSplitPoint(shared_ptr<IntersectionLink> link,
			     int fixed_dir, double fixed_value,
			     double par[], double& dist);

    void fetch_relevant_points(const std::vector<shared_ptr<IntersectionPoint> >& ipoints,
			       int missing_dir,
			       double missing_value);

    void constructor_implementation(shared_ptr<ParamObjectInt> obj1, 
				    shared_ptr<ParamObjectInt> obj2,
				    shared_ptr<IntersectionPool> parent,
				    int missing_dir,
				    double missing_value);

    void associate_parent_points(std::vector<shared_ptr<IntersectionPoint> >& children,
				 int lacking_ix,
				 double lacking_val);

    void transfer_links_to_parent_points(std::vector<shared_ptr<IntersectionPoint> >&
					 children, int lacking_ix);

    void get_param_limits(int pardir,
			  double& min_param, double& max_param) const;

    void
    set_correct_differentiating_domain(shared_ptr<IntersectionPoint>
                                       ip) const;

    void separate_isoparametric_curves(std::vector<std::vector<int> >&
				       curve_indices) const;

    void separate_at_problematic_points(std::vector<std::vector<int> >&
					curve_indices) const;

    int endof_noniso_curve(const std::vector<int>& cur_cv,
			   int start_ix) const;

    int endof_isoparametric_curve(const std::vector<int>& cur_cv,
				  int start_ix) const;

    // Remove curves which are judged unnecessary or just 'clutter'.
    void weedOutClutterCurves(std::vector<std::vector<int> >& curves) const;

    void locateEdgepoints(Array<std::vector<shared_ptr<IntersectionPoint> >, 4>& obj1_edges,
			  Array<std::vector<shared_ptr<IntersectionPoint> >, 4>& obj2_edges);

    void determine_free_dir_parameter(double* par, 
				      shared_ptr<IntersectionLink> link,
				      int free_dir,
				      int fixed_dir,
				      double fx_val);
};


//===========================================================================
//                  IMPLEMENTATION OF INLINE METHODS
//===========================================================================


//===========================================================================
inline void 
IntersectionPool::
getIntersectionPoints(std::vector<shared_ptr<IntersectionPoint> >& int_pts) const
//===========================================================================
{
    int_pts = int_points_;
}


//===========================================================================
inline void 
IntersectionPool::
getIntersectionCurves(std::vector<shared_ptr<IntersectionCurve> >& int_curves) const
//===========================================================================
{
    int_curves = int_curves_;
}


//===========================================================================
inline std::vector<shared_ptr<IntersectionPoint> >&
IntersectionPool::getIntersectionPoints()
//===========================================================================
{
    return int_points_;
}


//===========================================================================
inline const std::vector<shared_ptr<IntersectionPoint> >& 
IntersectionPool::getIntersectionPoints() const
//===========================================================================
{
    return int_points_;
}


//===========================================================================
inline const std::vector<shared_ptr<IntersectionCurve> >&
IntersectionPool::getIntersectionCurves() const 
//===========================================================================
{
    return int_curves_;
}


//===========================================================================
template<class ip_iterator>
inline void IntersectionPool::addCurves(int nmb_int_cvs, 
					ip_iterator* start_iterators, 
					ip_iterator* end_iterators)
//===========================================================================
{
    for (int i = 0; i < nmb_int_cvs; ++i) {
	addCurve(start_iterators[i], end_iterators[i]);
    }
}


//===========================================================================
template<class ip_iterator>
inline void IntersectionPool::addCurve(ip_iterator start, ip_iterator end)
//===========================================================================
{
    shared_ptr<IntersectionCurve> temp = std::make_shared<IntersectionCurve>(start, end);
    int_curves_.push_back(temp);
}


//===========================================================================
inline int IntersectionPool::lackingParameter() const
//===========================================================================
{
    return missing_param_index_;
}


//===========================================================================
inline double IntersectionPool::lackingParameterValue() const
//===========================================================================
{
    return missing_param_value_;
}


/// Helper struct to be used with
/// IntersectionPool::intersectAlongCommonBoundary()

//===========================================================================
struct BoundaryIntersectionData {
//===========================================================================

    /// Running parameter direction in object 1 and 2
    int dir[2];

    /// \brief Parameter value for the fixed parameter in object 1 and
    /// 2 (-1 if the object does not have a fixed parameter, i.e. when
    /// the object is a curve).
    double par[2];

    /// Chain of IntersectionPoints making up the intersection
    std::vector<shared_ptr<IntersectionPoint> > pts;

};


} // namespace Go


#endif  // _INTERSECTIONPOOL_H
