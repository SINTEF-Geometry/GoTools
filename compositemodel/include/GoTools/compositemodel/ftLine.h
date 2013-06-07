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

