//===========================================================================
//                                                                           
// File: DirectionCone.C                                                     
//                                                                           
// Created: Wed Jan 30 19:40:00 2002                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                           
// Revision: $Id: DirectionCone.C,v 1.22 2005-09-22 17:17:14 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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

	angle_ = 0.5 * (theta + angle_);
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

    os << centre_ << endl << angle_;
    return;
}
