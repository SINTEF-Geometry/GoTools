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

#include "GoTools/utils/DirectionCone.h"
#include <cmath>
#include <math.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace Go;


//===========================================================================
void
DirectionCone::setFromArray(const double* start, const double* end, int dim)
//===========================================================================
{
    ASSERT(end-start >= dim);
    ASSERT((end-start) % dim == 0);

    Point current(start, start + dim);
    const double* it = start;
    while (current.length() == 0.0) {
      it += dim;
      if (it == end) {
	MESSAGE("Setting a DirectionCone from an array of zero vectors");
	centre_ = current;
	angle_ = 0.0;
	greater_than_pi_ = 1;
	return;
      }
      current.setValue(it);
    }

//     if (current.length() == 0.0) {
// 	THROW("Cannot set cone from array containing zero vectors");
// 	return;
//     }

    centre_ = current;
    centre_.normalize();
    angle_ = 0.0;
    greater_than_pi_ = 0;

    for (const double* it = start + dim; it != end; it += dim) {
	current.setValue(it);
	if (current.length() > 0.0) {
	    addUnionWith(current);
	    if (greater_than_pi_ > 0)
		return;
	}
    }
    return;
}

//===========================================================================
bool DirectionCone::overlaps(const DirectionCone& cone) const
//===========================================================================
{
    if (centre_.length() == 0.0 || cone.centre_.length() == 0.0) {
	return false;   // Degenerate cone
    }
    if (greater_than_pi_ || cone.greater_than_pi_) {
	return true; // all directions covered
    }
    
    double ang = centre_.angle(cone.centre_);

    if (ang + angle_ + cone.angle_ < M_PI && angle_ + cone.angle_ < ang) {
	return false;
    }
    return true;
}


//===========================================================================
bool DirectionCone::perpendicularOverlaps(const DirectionCone& cone) const
//===========================================================================
{
    const double pihalf = 0.5*M_PI;
    if (centre_.length() == 0.0 || cone.centre_.length() == 0.0) {
	return false;   // Degenerate cone
    }

    if (greater_than_pi_ || cone.greater_than_pi_) {
	return true;
    }

    double ang = centre_.angle(cone.centre_);
    if (ang + angle_ < pihalf - cone.angle_ || 
	ang - pihalf - angle_ > cone.angle_) {
	return false;
    }
    return true;
}


//===========================================================================
bool DirectionCone::containsDirection(const Point& pt, double tol) const
//===========================================================================
{
    ALWAYS_ERROR_IF(greater_than_pi_ < 0, "Not initialized - cannot call.");
    if (greater_than_pi_) {
	return true; // surely covered
    }

    double ang = (centre_ * pt) / pt.length();
    ang = acos(ang);
    if (ang > angle_ + tol)
	return false;
    return true;
}


//===========================================================================
void DirectionCone::addUnionWith(const Point& pt)
//===========================================================================
{
    ALWAYS_ERROR_IF(greater_than_pi_ < 0,
		    "Must initialize DirectionCone");
    ALWAYS_ERROR_IF(centre_.dimension() != pt.dimension(),
		    "Dimension mismatch");

    Point t = pt;
    double tlength = t.length();
    if (tlength == 0.0) {
	MESSAGE("Ignoring vector of zero length");
	return;
    }
    t /= tlength;
    double theta = centre_ * t;
    if (theta < -1.0) // This test and the following is needed for
	theta = -1.0; // some reason... @jbt
    if (theta > 1.0)
	theta = 1.0;
    theta = acos(theta);
    if (theta <= angle_)
	return;
    if (theta + angle_ >= M_PI + zero_tol_) {
	greater_than_pi_ = 1;
	return;
    }

    // New code
    double t1 = sin(0.5 * (theta + angle_));
    double t2 = sin(0.5 * (theta - angle_));
    centre_ = centre_ * t1 + t * t2;
    if (centre_.length() == 0.0) {
	// The cone spans 180 degrees.
	//  Happens for instance if first 2 vectors point in opposite directions.
	greater_than_pi_ = 1;
	return;
// 	THROW("Can't normalize zero vector!");
    } else {
	centre_.normalize();
    }

    angle_ = 0.5 * (theta + angle_);
    check_angle();

    return;
}


//===========================================================================
void DirectionCone::addUnionWith(const DirectionCone& cone)
//===========================================================================
{
    ALWAYS_ERROR_IF(greater_than_pi_ < 0 || cone.greater_than_pi_ < 0,
		    "Must initialize DirectionCones.");
    ALWAYS_ERROR_IF(centre_.dimension() != cone.centre_.dimension(),
		    "Dimension mismatch.");

    if (cone.greater_than_pi_) {
	greater_than_pi_ = 1;
	return;
    }
    double theta = centre_ * cone.centre_;
    theta = acos(theta);
    if (theta + angle_ + cone.angle_ >= M_PI + zero_tol_) {
	greater_than_pi_ = 1;
	return;
    } else if (theta + cone.angle_ > angle_) {
	double t1 = sin(0.5 * (theta + angle_ - cone.angle_));
	double t2 = sin(0.5 * (theta - angle_ + cone.angle_));
	centre_ = centre_ * t1 + cone.centre_ * t2;
	centre_.normalize();

	angle_ = 0.5 * (theta + angle_ + cone.angle_);
    }
    check_angle();
    return;
}


//===========================================================================
void DirectionCone::check_angle() const
//===========================================================================
{
    ALWAYS_ERROR_IF(angle_ < 0.0, "Check failed -- negative angle.");

    if (angle_ <= M_PI + zero_tol_) 
	greater_than_pi_ = 0;
    else
	greater_than_pi_ = 1;

    return;
}


//===========================================================================
void DirectionCone::read(std::istream& is)
//===========================================================================
{
    ALWAYS_ERROR_IF(centre_.dimension() == 0,
		    "DirectionCone has no set dimension yet - cannot read.");
    is >> centre_ >> angle_;
    check_angle();
    return;
}


//===========================================================================
void DirectionCone::write(std::ostream& os) const
//===========================================================================
{
    ALWAYS_ERROR_IF(greater_than_pi_ < 0,
		    "Not initialized - cannot write.");

    streamsize prev = os.precision(15);
    os << centre_ << endl << angle_;
    os.precision(prev);   // Reset precision to it's previous value
    return;
}
