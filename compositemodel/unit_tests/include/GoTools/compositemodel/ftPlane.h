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

