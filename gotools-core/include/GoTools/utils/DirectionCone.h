//===========================================================================
//                                                                           
// File: DirectionCone.h                                                     
//                                                                           
// Created: Wed Jan 30 19:06:42 2002                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: DirectionCone.h,v 1.13 2009-05-13 07:33:08 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _DIRECTIONCONE_H
#define _DIRECTIONCONE_H

#include "GoTools/utils/Point.h"
#include "GoTools/utils/config.h"

namespace Go {

    /** Class representing a direction cone in Euclidian N-space
     *  The direction cone is characterised by a central direction
     *  and an angle 'alpha'.  All directions whose angle with the
     *  central direction is less than 'alpha' is considered to be
     *  covered by the cone.  If ever 'alpha' equals or exceeds
     *  180 degrees, any possible direction will end up being
     *  covered by the cone.
     */

class GO_API DirectionCone {
public:
    ///  Constructor making an undefined direction cone
    DirectionCone() : greater_than_pi_(-1), zero_tol_(1.0e-10) {}

    ///  Constructor making a direction cone whose central direction
    ///  is given by 'pt', and whose angle is zero.  The length of 'pt'
    ///  should be assured to be greater than zero, if not, then 
    ///  the DirectionCone will not be well defined.
    DirectionCone(const Point& pt)
	: centre_(pt), angle_(0.0), greater_than_pi_(0), zero_tol_(1.0e-10)
    {
	if (pt.length() > 0.0)
	    centre_.normalize();
    }

    /// Constructor making a direction cone whose central direction
    /// is given by 'centre', and whose angle (in radians) is 'angle'.
    /// The length of 'centre' should be assured to be greater than
    /// zero; if not, then the DirectionCone will not be well defined.
    DirectionCone(const Point& centre, double angle)
	: centre_(centre), angle_(angle), zero_tol_(1.0e-10)
    {
	if (centre_.length() > 0.0)
	    centre_.normalize();
	check_angle();
    }
    
    /// Do not inherit from this class -- nonvirtual destructor.
    ~DirectionCone() { }

    /// Re-define an already-existing DirectionCone from an array of 
    /// 'dim'-dimensional vectors stored in the memory area between
    /// 'start' and 'end'.  The first vector in this range will define
    /// the central direction of the DirectionCone, and cone's angle
    /// will be set wide enough to contain all the other directions given.
    /// NB: the length of the range (end - start) must be an exact 
    /// multiplum of 'dim'.
    void setFromArray(const double* start, const double* end, int dim);

    /// Return the dimension of the Euclidean space in which this 
    /// DirectionCone is defined
    int dimension() const { return centre_.size(); }

    /// Return the central direction of the DirectionCone
    const Point& centre() const { return centre_; }

    /// Return the angle of the DirectionCone
    double angle() const { return angle_; }

    /// Inform if the DirectionCone's angle is greater then PI. 
    /// If this is the case, then naturally ANY spatial direction
    /// is included in the cone.
    int greaterThanPi() const { return greater_than_pi_; }

    /// Return 'true' if there exists at least one direction that is 
    /// covered by both 'this' and the 'cone' DirectionCones.
    bool overlaps(const DirectionCone& cone) const;

    /// Return 'true' if there exist a direction in one cone that is
    /// perpendicular to a direction in the other cone.
    bool perpendicularOverlaps(const DirectionCone& cone) const;

    /// Return 'true' if the direction specified by 'pt' is contained
    /// within the DirectionCone.
    bool containsDirection(const Point& pt, double tol = 0.0) const;

    /// If necessary, increase the angle of the DirectionCone so that
    /// it covers the direction represented by 'pt'.
    void addUnionWith(const Point& pt);

    /// If necessary, increase the angle of the DirectionCone so that
    /// it covers all directions covered by 'cone'.
    void addUnionWith(const DirectionCone& cone);

    /// Read DirectionCone from stream
    void read(std::istream& is);

    /// Write DirectionCone to stream
    void write(std::ostream& os) const;

//     void setZeroTol(double zero_tol)
//     {
// 	zero_tol_ = zero_tol;
//     }

private:

    void check_angle() const;

    Point centre_;
    double angle_;
    mutable int greater_than_pi_;

    double zero_tol_; // By default the zero-tol is 0.0. Used in greaterThanPi().

};


} // namespace Go


#endif // _DIRECTIONCONE_H

