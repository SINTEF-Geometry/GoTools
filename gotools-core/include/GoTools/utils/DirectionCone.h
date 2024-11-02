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

    /// Test run. Which angle would an added vector result in. Angles larger
/// than M_PI is reported as M_PI.
double unionAngle(const Point& pt) const;

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

