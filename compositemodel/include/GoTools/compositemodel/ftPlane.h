//===========================================================================
//                                                                           
// File: ftPlane.h                                                           
//                                                                           
// Created: Thu Mar 23 11:34:23 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ftPlane.h,v 1.3 2008-04-02 07:18:27 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _FTPLANE_H
#define _FTPLANE_H

#include "GoTools/utils/Point.h"

namespace Go
{
    class BoundingBox;


//===========================================================================
/** ftPlane -  Representing a plane
 * 
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * 
 * 
 */
//===========================================================================

class ftPlane
{
private:

protected:

    Point normal_;  
    Point point_;

public:

    ftPlane()
      {}

    /// Constructor.
    /// First parameter is a vector normal to the plane,
    /// second parameter is a point in the plane.
    ftPlane(const Point& n, const Point& p)
	: normal_(n), point_(p)
	{}


    /// Constructor.
    /// Parameters are three points in the plane
    ftPlane(const Point& x, const Point& y, const Point& z)
      : point_(x)
    {
      normal_ = (z-x) % (y-x);
      normal_ /= normal_.length();
    }

    /// Copy constructor
    ftPlane(const ftPlane& plane);

    /// Destructor.
    ~ftPlane() {}

    /// Get plane normal
    const Point& normal() const { return normal_; }
    /// Get point in plane
    const Point& point()  const { return point_;  }
    /// Check if the plane intersects a bounding box
    bool intersectsBox(const BoundingBox& box) const;
};

} // namespace Go


#endif // _FTPLANE_H

