//===========================================================================
//                                                                           
// File: IntersectionPoint.h
//                                                                           
// Created: 
//                                                                           
// Author: oan
//                                                                           
// Revision: $Id: IntersectionPoint.h,v 1.76 2007-11-01 14:31:37 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _INTERSECTIONPOINT_H
#define _INTERSECTIONPOINT_H


#include "GoTools/intersections/IntersectionPointUtils.h"
#include "GoTools/intersections/LinkType.h"
#include "GoTools/intersections/GeoTol.h"
#include "GoTools/intersections/SecondOrderProperties.h"
#include <memory>
#include <vector>
#include <map>


namespace Go {


class ParamObjectInt;
class ParamGeomInt;
class ParamSurface;
class IntersectionLink;


/// Object describing a point located on the intersection of two
/// geometrical objects.

class IntersectionPoint {
public:
    /// Constructor. This constructor specifies an IntersectionPoint
    /// from all its defining components.
    /// \param obj1 the first of the two intersecting objects on which
    /// the IntersectionPoint lies.
    /// \param obj2 the second of the two intersecting objects on
    /// which the IntersectionPoint lies.
    /// \param epsge an object defining various tolerances
    /// \param obj_1_params parameters of the intersection point for
    /// the first object
    /// \param obj_2_params parameters of the intersection point for
    /// the second object
    IntersectionPoint(const ParamObjectInt* obj1, 
		      const ParamObjectInt* obj2,
		      const std::shared_ptr<GeoTol> epsge, 		      
		      const double* obj1_params,
		      const double* obj2_params);

    /// Constructor. This constructor is for creating an
    /// IntersectionPoint that is a replica of an existing one, but
    /// where one of the objects has one parameter less (i.e. is
    /// 'picked' from one of the objects of the existing
    /// IntersectionPoint by fixing one of its parameters).
    /// \param obj1 the first of the two intersecting objects on which
    /// the IntersectionPoint lies.  It is either the same as the
    /// first object defining the point \a ip, or picked from it (one
    /// parameter less).
    /// \param obj2 the second of the two intersecting objects on
    /// which the IntersectionPoint lies.  It is either the same as
    /// the second object defining the point \a ip, or picked from it
    /// (one parameter less).
    /// \param ip The already-existing IntersectionPoint that 'this'
    /// IntersectionPoint will represent.
    /// \param missing_param indicates which parameter of the
    /// point \a ip is 'missing' for 'this' IntersectionPoint.
    /// (Parameters are numbered from 0 and up to the total number of
    /// parameters in 'obj1' and 'obj2').
    IntersectionPoint(const ParamObjectInt* obj1,
		      const ParamObjectInt* obj2,
		      const std::shared_ptr<IntersectionPoint> ip,
		      int missing_param);

    /// Destructor
    ~IntersectionPoint();

    /// Default constructor
    IntersectionPoint() {} // create an undefined point which can be
			   // assigned or read() into.
    
    /// Write IntersectionPoint to stream (NB: topological and parent
    /// information will be lost)
    /// \param os output stream
    void write(std::ostream& os) const;

    /// Read IntersectionPoint from stream (NB: no topological or
    /// parent information)
    /// \param is input stream
    // @@This function is for debugging.  Do not
    // try to do anything serious with it.  It must create 'new'
    // objects for its pointers, and thus causes memory leaks!
    void read(std::istream& is);

    /// Write the parameters of the intersection point to stream
    /// \param os output stream
    void writeParams(std::ostream& os) const;

    /// Write available info about an IntersectionPoint to standard
    /// output. Makes use of the struct IntPtInfo.
    void writeInfo() const;

    /// Replace the parameters defining this IntersectionPoint and
    /// recalculate internal information (the underlying geometric
    /// objects are kept).
    /// \param param pointer to an array of \c double s, defining the
    /// new values for the point's parameters.  The array should of
    /// course have length equal to the total number of parameters
    /// defining this IntersectionPoint.
    void replaceParameter(double *param);
    
    /// Get average of the IntersectionPoint's position in space as
    /// evaluated in the two underlying objects.
    /// \return the position of the IntersectionPoint in space
    Point getPoint() const;

    /// Get the IntersectionPoint's position in space as evaluated in
    /// the first underlying object.
    /// \return the position of the IntersectionPoint in space
    Point getPoint1() const
    { return point1_; }

    /// Get the IntersectionPoint's position in space as evaluated in
    /// the second underlying object.
    /// \return the position of the IntersectionPoint in space
    Point getPoint2() const
    { return point2_; }

    /// This function currently only works (and makes sense) for 3D
    /// geometrical objects (not functions).  Given a vector 'dir', it
    /// projects this vector onto the two objects at the
    /// IntersectionPoint, and calculates this direction in the two
    /// parameter spaces.  These normalized directions are returned in
    /// 'par_1_dir' and 'par_2_dir'.
    /// \param dir the direction vector that the user wants to project
    /// on the objects
    /// \retval par_1_dir the direction of the projected vector in the
    /// parameter space of the first object.
    /// \retval par_2_dir the direction of the projected vector in the
    /// parameter space of the second object.
    void projectToParamPlanes(const Point& dir, 
			      Point& par_1_dir, Point& par_2_dir) const;

    /// Check if there is a unique tangent direction for the
    /// intersection at the position of this IntersectionPoint.  This
    /// is the case for ORDINARY_POINTs and TANGENTIAL_POINTs, and
    /// false in all other cases (see \ref SingularityType).  See also
    /// \ref tangentIsOriented().
    /// \return 'true' if the IntersectionPoint has a unique tangent
    /// direction (not necessarily uniquely oriented).  'false'
    /// otherwise.
    bool hasUniqueTangentDirection() const; 

    /// Returns true if this IntersectionPoint is sure about the
    /// orientation of its tangent.  It is sure about the orientation
    /// if the intersection is transversal, or if the intersection is
    /// singular and the orientation has been set manually through the
    /// fixTangentOrientation() command.  See also \ref
    /// hasUniqueTangentDirection().
    /// \return 'true' if the \em tangent \em orientation is clear (a
    /// prerequisite is of course that the \em tangent \em direction
    /// is clear). 'false' otherwise.
    bool tangentIsOriented() const;

    /// Distance between the position of the IntersectionPoint as
    /// described in the first object and as described in the second
    /// object.
    /// \return the distance between the two points describing the
    /// IntersectionPoint in each of the objects.
    double getDist() const
    { return dist_; }

    /// Get a vector containing the IntersectionPoint's parameter
    /// values.
    const std::vector<double>& getPar() const
    { return par_; }

    /// Get a specified parameter value of the IntersectionPoint
    /// \param idx the number of ther parameter for which we seek the
    /// value
    /// \return the specified parameter value
    double getPar(int idx) const
    { return par_[idx]; }

    /// Get a pointer to the parameter values in the first object
    /// \return a pointer to the array where the parameter values for
    /// the first object are (consecutively) stored.
    const double* getPar1() const;

    /// Get a pointer to the parameter values in the second object
    /// \return a pointer to the array where the parameter values for
    /// the second object are (consecutively) stored.
    const double* getPar2() const;

    /// Get a Point representing the IntersectionPoint in the
    /// parameter domain of the first object
    Point getPar1Point() const;
    /// Get a Point representing the IntersectionPoint in the
    /// parameter domain of the second object
    Point getPar2Point() const;

    /// Get number of parameter values in object number 1
    int numParams1() const;

    /// Get number of parameter values in object number 2
    int numParams2() const;

    /// Get the tangent of the intersection at the point.  If the
    /// point is a branch point (several intersection curves meet at
    /// this point), then the user has to specify the branch that she
    /// wants the tangent to.
    /// \param second_branch usually set to 'false' (default).  If the
    /// IntersectionPoint lies at a branchpoint, then the user can set
    /// this parameter to 'true' in order to get the tangent of the
    /// \em second of the two branches that supposedly meet in this
    /// point.  However, setting this parameter to 'true' \em will
    /// cause an error if the IntersectionPoint is \em not a branch
    /// point.
    /// \return the tangent of intersection curve at this point.
    Point getTangent(bool second_branch = false) const;

    /// Get the parameter direction along intersection in first object
    /// at this IntersectionPoint.  If this point is a branch point
    /// (several intersection curves meet at this point), then the
    /// user has to specify the branch that she wants the tangent to.
    /// \param second_branch usually set to 'false' (default).  If the
    /// IntersectionPoint lies at a branchpoint, then the user can set
    /// this parameter to 'true' in order to get the parameter
    /// direction for the \em second of the two branches that
    /// supposedly meet in this point. However, setting this parameter
    /// to 'true' \em will cause an error if the IntersectionPoint is
    /// \em not a branch point.
    /// \return the parameter direction along the intersection in the
    /// first object at this point.
    Point getPar1Dir(bool second_branch = false) const;

    /// Get the parameter direction along intersection in second
    /// object at this IntersectionPoint.  If this point is a branch
    /// point (several intersection curves meet at this point), then
    /// the user has to specify the branch that she wants the tangent
    /// to.
    /// \param second_branch usually set to 'false' (default).  If the
    /// IntersectionPoint lies at a branchpoint, then the user can set
    /// this parameter to 'true' in order to get the parameter
    /// direction for the \em second of the two branches that
    /// supposedly meet in this point. However, setting this parameter
    /// to 'true' \em will cause an error if the IntersectionPoint is
    /// \em not a branch point.
    /// \return the parameter direction along the intersection in the
    /// second object at this point.
    Point getPar2Dir(bool second_branch = false) const;

    /// Get pointer to first object
    const ParamObjectInt* getObj1() const
    { return obj1_; }

    /// Get pointer to second object
    const ParamObjectInt* getObj2() const
    { return obj2_; }

    /// Get the start value of the parameter interval for the
    /// specified parameter.
    /// \param pardir the parameter for which we seek the start value.
    /// \return the start value of the parameter interval for the
    /// parameter 'pardir'.
    double startParam(int pardir);

    /// Get the end value of the parameter interval for the specified
    /// parameter.
    /// \param pardir the parameter for which we seek the end value.
    /// \return the end value of the parameter interval for the
    /// parameter 'pardir'.
    double endParam(int pardir);

    /// Return 'true' if the point was found to be singular, ie. not
    /// on a normal, transversal intersection (see \ref
    /// SingularityType).
    bool pointIsSingular() const
    { return getSingularityType() != ORDINARY_POINT; }

    /// Return 'true' if the point was found to be near-singular, that
    /// is, the tangent calculated in this point is very small.
    bool isNearSingular() const;

    /// Return the tolerance for the parameter in the direction
    /// 'pardir'. This is equal to GeoTol::rel_par_res_ multiplied by
    /// the span of the parameter.
    /// \param pardir the parameter for which we want to find the
    /// parameter tolerance
    /// \return the tolerance for the specified parameter.
    double parameterTolerance(int pardir); // pardir == 0 || pardir == 1

    /// If the tangent is not oriented (because of singular
    /// intersection), it can be set manually with this command.  If
    /// the argument is "true", the tangent's orientation will be
    /// flipped.  In the opposite case, it will be kept as it is.
    /// Subsequent calls to the \ref tangentIsOriented() function will
    /// now report that the tangent is oriented.  Not to be called for
    /// an ORDINARY_POINT (which should have been handled already).
    /// See \ref SingularityType for a definition of what
    /// ORDINARY_POINT means.
    /// \param flip if 'true', then the tangent's orientation will be
    /// flipped, otherwise it will be kept as it is.  In any case, the
    /// IntersectionPoint will now consider the orientation of its
    /// tangent to be 'fixed'.
    void fixTangentOrientation(bool flip); 

    /// The first and second argument are set to 'true' or 'false'
    /// according to whether the point is lying on the boundary of
    /// respectively the first and second object.
    /// \retval first 'true' if this IntersectionPoint is located on
    /// the boundary of the first object, 'false' otherwise.
    /// \retval second 'true' if this IntersectionPoint is located on
    /// the boundary of the second object, 'false' otherwise.
    void isBoundaryPoint(bool& first, bool& second) const;

    /// The first and second argument are set to 'true' or 'false'
    /// according to whether the tangent direction in this point is
    /// clearly oriented toward the inside of the domain for
    /// respectively the ffirst and second surface.  It can only be
    /// pointing outwards if the point is lying on the boundary.  If
    /// the point is on the boundary and its tangent is close to
    /// parallel with the boundary, 'first_along_boundary' and/or
    /// 'second_along_boundary' will be set to true, otherwise they
    /// are false.  Of course, points whose tangents are directed
    /// along the boundary will return 'first' and/or 'second' false.
    /// \retval first 'true' if the tangent direction of the
    /// intersection at this point is clearly oriented toward the
    /// inside of the domain of the first object. 'false' if not.
    /// \retval second 'true' if the tangent direction of the
    /// intersection at this point is clearly oriented toward the
    /// inside of the domain of the second object. 'false' if not.
    /// \retval first_along_boundary 'true' if the point is on the
    /// boundary of the first object, and with a tangent that is close
    /// to parallel with the boundary. 'false' otherwise.
    /// \retval second_along_boundary 'true' if the point is on the
    /// boundary of the second object, and with at tangent that is
    /// close to parallel with the boundary. 'false' otherwise.
    void tangentPointingInwards(bool& first, 
				bool& second,
				bool& first_along_boundary,
				bool& second_along_boundary) const;

    /// \brief Check if this point is topologically connected with
    /// another IntersectionPoint (i.e. that there is an
    /// IntersectionLink between them).
    ///
    /// \param p a shared pointer to the point we are checking
    /// against.
    /// \return 'true' if there is an IntersectionLink between 'this'
    /// point and \a p.  'false' otherwise.
    bool isConnectedTo(std::shared_ptr<IntersectionPoint> p) const
    { return isConnectedTo(p.get()); }

    /// \brief Check if this point is topologically connected with
    /// another IntersectionPoint (i.e. that there is an
    /// IntersectionLink between them).
    ///
    /// \param point a pointer to the point we are checking against.
    /// \return 'true' if there is an IntersectionLink between 'this'
    /// point and \a point.  'false' otherwise.
    bool isConnectedTo(const IntersectionPoint* point) const;

    /// This function topologically "connects" this point with another
    /// IntersectionPoint.  It generates a new IntersectionLink
    /// combining these points.  This IntersectionLink can inherit the
    /// meta-information from \a model_link. This is typically useful
    /// when an IntersectionPoint is "inserted" between two
    /// already-connected points, and the two new links should contain
    /// the same information as the older, split-up link.  Another
    /// example where information from \a model_link needs to be
    /// copied is when parents of two connected points needs to be
    /// connected.  In that case, it is important that we know which
    /// parameter direction has been introduced when moving to the
    /// parent, so that we can "shift" the meta-information pertaining
    /// to individiual parameters accordingly.  If no new parameter
    /// information has been introduced, \a added_parameter_dir should
    /// be left at its default value of -1.
    ///  \param point the IntersectionPoint we want to connect to \c
    /// this point
    /// \param type the type of the link in terms of a LinkType
    /// \param model_link if this parameter is specified, the
    /// information associated with \a model_link will be copied by
    /// the newly created IntersectionLink.
    /// \param added_parameter_dir relevant parameter when \a
    /// model_link contains one parameter direction less than the new
    /// IntersectionLink does.  In this case, \a added_parameter_dir
    /// will specify the index of the parameter that is lacking.
    /// \return shared pointer to the generated link
    std::shared_ptr<IntersectionLink>
    connectTo(IntersectionPoint *const point, 
	      LinkType type,
	      std::shared_ptr<IntersectionLink> model_link
	      = std::shared_ptr<IntersectionLink>(),
	      int added_parameter_dir = -1);

    /// This function topologically "connects" this point with another
    /// IntersectionPoint.  It generates a new IntersectionLink
    /// combining these points.  This IntersectionLink can inherit the
    /// meta-information from \a model_link. This is typically useful
    /// when an IntersectionPoint is "inserted" between two
    /// already-connected points, and the two new links should contain
    /// the same information as the older, split-up link.  Another
    /// example where information from \a model_link needs to be
    /// copied is when parents of two connected points needs to be
    /// connected.  In that case, it is important that we know which
    /// parameter direction has been introduced when moving to the
    /// parent, so that we can "shift" the meta-information pertaining
    /// to individiual parameters accordingly.  If no new parameter
    /// information has been introduced, \a added_parameter_dir should
    /// be left at its default value of -1.
    /// \param point the IntersectionPoint we want to connect to \c
    /// this point
    /// \param type the type of the link in terms of a LinkType
    /// \param model_link if this parameter is specified, the
    /// information associated with \a model_link will be copied by
    /// the newly created IntersectionLink.
    /// \param added_parameter_dir relevant parameter when \a
    /// model_link contains one parameter direction less than the new
    /// IntersectionLink does.  In this case, \a added_parameter_dir
    /// will specify the index of the parameter that is lacking.
    /// \return shared pointer to the generated link
    std::shared_ptr<IntersectionLink>
    connectTo(std::shared_ptr<IntersectionPoint> point, 
	      LinkType type,
	      std::shared_ptr<IntersectionLink> model_link
	      = std::shared_ptr<IntersectionLink>(),
	      int added_parameter_dir = -1)
    { return connectTo(point.get(), type,
		       model_link, added_parameter_dir); }

    /// Topologically disconnect 'this' IntersectionPoint from another
    /// IntersectionPoint (which means removing whatever
    /// IntersectionLink there might be between them).
    /// \param point the IntersectionPoint we want to disconnect from
    /// 'this' IntersectionPoint.
    void disconnectFrom(IntersectionPoint* point);
    
    /// Returns how many 'tangential branches' emanates from this
    /// point (usually 1, but 2 for branch points, and 0 for isolated-
    /// or higher-order points).  See \ref SingularityType for a
    /// classification of IntersectionPoint s.
    int numBranches() const;

    /// Classifies an IntersectionPoint according to the direction of
    /// the corresponding intersection curve with respect to the
    /// parameter domain of two surfaces.  \a branch_num indicates
    /// which tangential branch to use (usually there's just one
    /// (\a branch_num should be 0), but for branchpoints there are two
    /// (\a branch_num can be 0 or 1).
    /// \param surf1 the first surface
    /// \param surf2 the second surface
    /// \param branch_num indicates which tangential branch to use
    /// \return the classification of this IntersectionPoint according
    /// to the domain of the given surfaces. See also \ref
    /// IntPtClassification.
    IntPtClassification getClassification(const ParamSurface *surf1, 
					  const ParamSurface *surf2,
					  int branch_num = 0) const;

    /// Classifies an IntersectionPoint according to the direction of
    /// the corresponding intersection curve with respect to the
    /// parameter domain of a \em function.  \a branch_num indicates
    /// which tangential branch to use (usually there's just one
    /// (\a branch_num should be 0), but for branchpoints there are two
    /// (\a branch_num can be 0 or 1).
    /// \param par_func the function, expressed as a 1D ParamSurface
    /// \param C the iso-value that defines the intersection on this
    /// function.
    /// \param branch_num indicates which tangential branch to use
    /// \return the classification of this IntersectionPoint according
    /// to the domain of the given function. See also \ref
    /// IntPtClassification.
    IntPtClassification getClassification(const ParamSurface *par_func,
					  double C,
					  int branch_num = 0) const;

    /// Get the SingularityType of this IntersectionPoint.
    SingularityType getSingularityType() const;

    /// Get the object containing the \em tolerances used by this
    /// IntersectionPoint
    std::shared_ptr<GeoTol> getTolerance() 
    { return epsge_; }

    /// Get the object containing the \em tolerances used by this
    /// IntersectionPoint
    std::shared_ptr<const GeoTol> getTolerance() const
    { return epsge_; }

    /// Get a shared pointer to the parent point of this
    /// IntersectionPoint.  If it has no parent, the shared pointer
    /// will be a null-pointer.
    std::shared_ptr<IntersectionPoint> parentPoint()
    { return parent_point_; }

    /// Set a given IntersectionPoint to be the parent point of this
    /// IntersectionPoint.
    /// \param p the point that shall be the parent point of this
    /// IntersectionPoint.
    void setParentPoint(std::shared_ptr<IntersectionPoint> p)
    {
	//ASSERT(!parent_point_); // should only be used with orphans
	parent_point_ = p;
    }

    bool hasParentPoint()
	{
	    return (parent_point_.get() != 0);
	}

    /// Flip the two intersecting objects (useful sometimes for
    /// self-intersections, should be avoided otherwise).
    void flipObjects();

    /// Returns the parameter value corresponding to the influence
    /// area of the intersection point in a specified direction and
    /// orientation.  Since the bounds of the influence area are
    /// bracketed, we do not have the exact value of its limits.  If
    /// 'first_outside' is 'true', we will return the first parameter
    /// value detected outside the influence area.  Otherwise, we will
    /// return the last parameter detected inside the influence area.
    /// \param dir the parameter direction for which we seek the
    /// influence area.
    /// \param forward 'true' if we want to locate the upper end of
    /// the point's influence area, 'false' if we want to locate the
    /// lower end.
    /// \param first_outside 'true' if we want the returned value to
    /// be the first parameter value \em outside the influence area
    /// that was found by the function. 'false' if we seek the last
    /// detected value \em inside the influence area.
    /// \return a value that (according to the function arguments)
    /// either represents the last found parameter value \em inside or
    /// first found parameter value \em outside the influence area of
    /// this IntersectionPoint along the specified parameter
    /// direction.
    double getInfluenceArea(int dir, bool forward, bool first_outside,
			    double aeps=0.0) const;

    /// Returns 0 if the given parameter value \a par was not in the
    /// point's influence area for the point's parameter \a pardir.
    /// Returns 1 if it is inside the point's influence area, but not
    /// exactly on the point.  Returns 2 if it is exactly on the
    /// point.  Since the bounds of the influence area are bracketed,
    /// we do not have the exact value of its limits.  If \a
    /// first_outside is 'true', we will check against the first
    /// parameter value detected outside the influence area.
    /// Otherwise, we will check against the last parameter detected
    /// inside the influence area.
    /// \param pardir the concerned parameter directon
    /// \param par the concerned parameter value
    /// \param first_outside specify whether we want to check against
    /// the first detected value \em outside the influence area
    /// ('true'), or the last detected value \em inside it ('false').
    /// \return 0, 1 or 2, according to the above description.
    int inInfluenceArea(int pardir, double par, bool first_outside) const;

    /// This function is mostly for optimizing reasons.  This point is
    /// supposed to lie inside the influence area of \a other_pt, and
    /// by this function, this influence interval is split up and
    /// shared between this point and the other point.  By deducing
    /// this interval directly from \a other_pt, we avoid having to
    /// march it out again.  The concerned parameter direction is
    /// \a pardir.
    /// \param other_pt refers to the IntersectionPoint that we assume
    /// share influence area with \c this IntersectionPoint.
    /// \param pardir the concerned parameter directon
    void shareInfluenceAreaWith(std::shared_ptr<IntersectionPoint> other_pt,
				int pardir);

    // @@@ VSK. For the time being, but needs implementing when the
    // CachedInterval::inside is OK
    /// This function is not properly defined at the moment.
    /// \param pardir the concerned parameter directon
    /// \param parval the concerned parameter value
    bool inInfluenceAreaBracket(int pardir, double parval)
    { return false; }

    /// Return the intersection link connecting \c this point with
    /// \a point.  If no such link exists, a null pointer will be
    /// returned.
    /// \param point the point for which we seek the IntersectionLink
    /// to this point.
    /// \return the IntersectionLink in question.
    std::shared_ptr<IntersectionLink>
    getIntersectionLink(const IntersectionPoint* point);
    
    /// Get a vector containing pointers ot all IntersectionPoints
    /// topologically connected with \c this
    /// IntersectionPoint. IntersectionPoints are topologically
    /// connected if they share an IntersectionLink.
    /// \retval pnts upon function return, this will be the vector
    /// sought for.
    void getNeighbours(std::vector<IntersectionPoint*>& pnts) const;

    /// Get a vector containing all the IntersectionLinks that
    /// defines connections from \c this IntersectionPoint.
    /// \retval links upon function return, this will be the vector
    /// sought for.
    void getNeighbourLinks(std::vector<std::
			   shared_ptr<IntersectionLink> >& links) const;

    /// Query the number of IntersectionPoints ('neighbours') that
    /// are connected to \c this IntersectionPoint via an
    /// IntersectionLink.
    int numNeighbours()const { return int(intersection_links_.size()); }

    /// \a it should point to a p-array of \c bool, where \a p is the
    /// total number of parameters for this IntersectionPoint.  The
    /// value at \a it[i] determines whether second-order differentiation
    /// for the parameter \a i should be carried out from the left,
    /// rather than from the right (default).  This sets an internal
    /// state variable that will not be changed until the next call of
    /// setDifferentiateFromLeft().  (The reason a template was used
    /// here is to be able to use iterators to vectors of bool, which
    /// are special in that they cannot be accessed in the usual
    /// pointer way).
    /// \param it an iterator (pointer, etc.) to the aforementioned
    /// array of \c bool.
    template<class bool_iterator>
    void setDifferentiateFromLeft(bool_iterator it) const  
    {
	int n = numParams1() + numParams2();
	ASSERT(n == int(diff_from_left_.size()));
	for (int i = 0; i < n; ++i, ++it) {
	    diff_from_left_[i] = g2_discontinuous_params_[i] && (*it);
	}
	int key = bool2int(diff_from_left_.begin(), diff_from_left_.end());
	cur_active_2nd_order_properties_ = &sec_order_properties_[key];
    }

    /// Returns whether the internal state for parameter 'pardir' is
    /// set so that it differentiates from left.
    /// \param pardir the concerned parameter direction
    /// \return 'true' if differentiation is set to be carried out
    /// 'from the left' for this parameter, 'false' otherwise.
    bool differentiatesFromLeft(int pardir) const
    {
	ASSERT(pardir >= 0 && pardir < numParams1() + numParams2());
	return diff_from_left_[pardir];
    }

    /// Returns the number of parameter directions for which this
    /// point is connected with links that are isoparametric.  The
    /// indices of the isoparametric parameter directions are written
    /// to an array pointed to by 'iso_par_dirs' (allocated by the
    /// user).
    /// \param iso_par_dirs pointer to the array where the indices
    /// will be written.
    /// \return the number of parameter directions for which this
    /// point is connected with links that are isoparametric.
    int hasIsoLinks(int* const iso_par_dirs) const;

    /// Returns the number of parameter directions for which all
    /// IntersectionLink s of this IntersectionPoint are
    /// isoparametric. The indices of these parameter directions are
    /// written to the array pointed to by 'iso_par_dirs' (allocated
    /// by the user).
    /// \param iso_par_dirs pointer to the array where the indices
    /// will be written.
    /// \return the number of parameter directions for which all
    /// IntersectionLink s of this IntersectionPoint are
    /// isoparametric.
    int commonIsoLinks(int* const iso_par_dirs) const;

    /// Check if there is a G2-discontinuity for one of the parameters
    /// at this IntersectionPoint.
    /// \return 'true' if at least one G2-discontinuity was found,
    /// 'false' otherwise.
    bool containsG2Discontinuity() const 
    {
	int ki;
	for (ki=0; ki < int(g2_discontinuous_params_.size()); ki++)
	    if (g2_discontinuous_params_[ki])
		return true;
	return false;
    }

    /// Check if the point lies on a degenerate edge (or other
    /// degenerate part) of the surface.
    /// \return 'true' if degenerate, 'false' otherwise.
    bool isDegenerate() const;

    /// Check if the given intersection point represents the same
    /// point as 'this' point within the tolerance.
    /// \param point pointer to the intersection point we wish to test
    /// on
    /// \retval 'true' if the points are the same, 'false' otherwise
    bool isSamePoint(const IntersectionPoint* point) const;

    /// Check if the intersection point is OK, and get some info about
    /// the point.
    /// \param int_pt_info struct of type IntPtInfo with info about
    /// the point. Indicates if the point is ok, the number of
    /// neighbours, the singularity type, and the location of the
    /// point is the parameter domains of the objects
    /// \return 'true' if the intersection points is OK, 'false'
    /// otherwise
    bool checkIntersectionPoint(IntPtInfo& int_pt_info) const;

    /// Check if the point is in the domain specified by a parameter
    /// box.
    /// \param frompar the array of parameters defining the lower
    /// corner of the box in question.
    /// \param topar the array of parameters defining the upper corner
    /// of the box in question.
    /// \retval \a true if inside domain, \a false otherwise.
    bool isInDomain(double* frompar, double* topar) const;

    /// Global tolerance for tangent length (considered singular if
    /// length is less than this value.
    static const double tangent_tol;

private:

    // ----- Data members -----

    // NB: the IntersectionPoint has no ownership to the objects
    // pointed to by obj1_ and obj2_
    const ParamObjectInt* obj1_; // Object number one
    const ParamObjectInt* obj2_; // Object number two

    Point point1_; // Intersection point in space in object number one
    Point point2_; // Intersection point in space in object number two
    double dist_;  // Distance between the points in the two objects

    // Parameter values of this intersection point
    std::vector<double> par_;

    // If obj1_ or obj2_ has parents, then parent_point_ is the
    // corresponding intersection point
    std::shared_ptr<IntersectionPoint> parent_point_;

    // Contains link objects to other IntersectionPoints. The links
    // are unique.
    std::vector<std::shared_ptr<IntersectionLink> > intersection_links_;
    
    // Tolerance-related stuff
    std::shared_ptr<GeoTol> epsge_;

    // Caching of influence area
    mutable std::vector<CachedInterval> cached_influence_area_forwards_;
    mutable std::vector<CachedInterval> cached_influence_area_backwards_;

    // Information (that might be) influenced by the second derivative
    // of the surface.
    std::vector<bool> g2_discontinuous_params_;
    mutable std::map<int, SecondOrderProperties> sec_order_properties_;
    mutable SecondOrderProperties* cur_active_2nd_order_properties_;
    mutable std::vector<bool> diff_from_left_;

    // ------ member functions -----

    void constructor_implementation(const ParamObjectInt* obj1, 
				    const ParamObjectInt* obj2,
				    std::shared_ptr<GeoTol> epsge,
				    const double* obj_1_params,
				    const double* obj_2_params);

    inline bool tangent_2d_is_cached() const
    { return cur_active_2nd_order_properties_->tangent_2d_is_cached; }

    inline bool singularity_info_is_cached() const
    { return cur_active_2nd_order_properties_->singularity_info_is_cached; }

    // Expecting object(s) to be 2-parametric.
    void compute_parametrical_tangents() const; 
    // The two following function are called (only) by
    // compute_parametrical_tangents())
    void compute_surface_parametrical_tangents_safe() const;
    void compute_function_parametrical_tangents() const; 

    // Calculates the tangent at a discovered singular point, and
    // return the point's SingularityType
    SingularityType calculate_tangent_at_singular_point() const;

    // Check each parameter for a second order discontinuity in the
    // IntersectionPoint
    //std::shared_ptr<ParamObjectInt> obj2,
    static std::vector<bool> 
    detect_2nd_order_discontinuities(const ParamObjectInt* o1,
				     const ParamObjectInt* o2,
				     const std::shared_ptr<GeoTol> gtol,
				     const double* obj1_par,
				     const double* obj2_par);

    void calculate_singularity_type() const;
    void calculate_singularity_type_function() const;
    void calculate_singularity_type_surface() const; // Assuming space is 3-d
    void calculate_singularity_type_curves() const ;
    void calculate_singularity_type_point() const;
    void calculate_singularity_type_curve_surface() const;

    static double determinant(const Point&, const Point&, const Point&);
    
    // We suppose that 'T' lies in the plane spanned by 'b1' and 'b2'.
    // We want to calculate the coefficients in order to express 'T'
    // as a linear combination of 'b1' and 'b2'.

    // Calculate the vector 'T', supposedly in the linear span of 'b1'
    // and 'b2', as a linear combination of the two.  The appropriate
    // coefficients are returned in 'b1_coef' and 'b2_coef'.  If 'T'
    // is _not_ in the plane, the projected vector is used.
    static void decompose(const Point& T, 
			  const Point& b1, 
			  const Point& b2, 
			  double& b1_coef,
			  double& b2_coef);

    // Project a vector to the tangent space of a geometric object,
    // and find the parametric direction in the object's parameter
    // space that correspond to this projection.  'dir' is the vector
    // to be projected.  The object is given by 'go', and the
    // parametrization for the point where the tangent space is
    // defined, is given by 'par'.  Whether the tangent space should
    // be computed used derivation "from left" or "from right" in each
    // parameter direction is dictated by 'from_left'.  The resulting
    // direction in parametric space is returned by 'par_dir'.
    static void proj_2_pplane(const Point& dir, 
			      const double* par,
			      std::vector<bool>::const_iterator from_left,
			      const ParamGeomInt* go, 
			      Point& par_dir);

    // Takes a range of bools and converts them to an int. The
    // resulting int has the binary representation that corresponds to
    // the range of bools.
    template<class iterator>
    static int bool2int(iterator begin_range, iterator end_range)
    {
	int result = 0;
	for (int pos = 0; begin_range != end_range; ++pos, ++begin_range) {
	    // If the bool at the current position is true, then set the
	    // bit at position 'pos' from the right in integer 'result'
	    if (*begin_range) {
		result |= (1 << pos); 
	    }
	}
	return result;
    }
    
    // Takes an int and converts it to a range of bools. The resulting
    // bools correspond to the binary representation of the int.
    template<class iterator>
    static void int2bool(int number, iterator cur_pos, iterator end_range)
    {
	for (int pos = 0; cur_pos != end_range; ++pos, ++cur_pos) {
	    // Set the bool at the current position to true if 'number'
	    // has the bit at position 'pos' from the right turned on
            *cur_pos = (number & (1 << pos)) != 0;
	}
    }

    std::map<int, SecondOrderProperties> 
    generate_2nd_order_property_map(const std::vector<bool>& disc_params);

    inline bool exactly_on_point(int pardir, double par) const
    { return (fabs(par_[pardir] - par) < epsge_->getRelParRes()); }

    // Should be called whenever the IntersectionPoint is changed
    void clearCache() const
    {
	std::map<int, SecondOrderProperties>::iterator it;
	for (it = sec_order_properties_.begin();
	     it != sec_order_properties_.end(); ++it) {
	    (it->second).clear();
	}
	cached_influence_area_forwards_
	    = cached_influence_area_backwards_
	    = std::vector<CachedInterval>(par_.size());
	diff_from_left_ = std::vector<bool>(par_.size(), false);
	cur_active_2nd_order_properties_ = &sec_order_properties_[0];
    }

};


}; // namespace Go


#endif  // _INTERSECTIONPOINT_H
