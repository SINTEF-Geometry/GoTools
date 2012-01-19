//===========================================================================
//                                                                           
// File: ftLine.h                                                           
//                                                                           
// Created: Feb. 2007
//                                                                           
// Author: Vibeke Skytt, SINTEF
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _FTLINE_H
#define _FTLINE_H

#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/ftPlane.h"



#include "GoTools/utils/Point.h"

namespace Go
{

//===========================================================================
/** ftLine -  A line
 * 
//===========================================================================
*/

class GO_API ftLine
{
private:

protected:

    Point dir_;  
    Point point_;

public:
    /// Default constructor
    ftLine()
      {}

    /// Constructor.
    /// First parameter is a vector normal to the plane,
    /// second parameter is a point in the plane.
    ftLine(const Point& dir, const Point& pnt)
	: dir_(dir), point_(pnt)
	{}

    /// Destructor.
    ~ftLine();

    /// Find two planes intersecting in this line
    void getTwoPlanes(ftPlane& plane1, ftPlane& plane2) const;

    /// Determine if this line comes close enough to a given point after
    /// a rescaling of pt's coordinates
    bool closeToScaledPoint(const Point& pt,      // Input point before scaling
			    const Point& scale,   // Scale factor
			    double dist2) const;  // Square of maximum distance

    /// Return a new line corresponding to this line after a rescaling of the coordinates
    ftLine scaled(const Point& scale) const;

    /// Determine if the line intersects the smallest eliposide containing a box
    /// Is about 2 times faster than intersectsBox, but returns true
    /// about 50% more often than intersectsBox
    bool intersectsSphereOfBox(const BoundingBox& box) const;

    /// Determine if the line intersects a box
    bool intersectsBox(const BoundingBox& box) const;

    /// Determine if two surfaces of the line intersect a box
    /// Returns true about 45% more often than intersectsBox
    bool planesIntersectBox(const BoundingBox& box) const;

    /// The direction of this line
    const Point& direction() const { return dir_; }
    /// A point on this line
    const Point& point()  const { return point_;  }

};

} // namespace Go

#endif // _FTPLANE_H

