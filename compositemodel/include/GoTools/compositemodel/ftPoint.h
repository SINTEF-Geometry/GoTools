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

#ifndef _FTPOINT_H
#define _FTPOINT_H


#include "GoTools/utils/Point.h"

namespace Go
{

class ftSurface;

//===========================================================================
/** ftPoint - represents a point, possibly lying on a surface
 * ftPoint contains representations for both a 3D point in space
 * and possibly also parameters for the point on a specific surface.
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 */
//===========================================================================
class ftPoint
{
protected:
    Point pt_;              // The point
    ftSurface* surface_;         // The surface on which it lies (if nonzero)
    double u_, v_;            // Coordinates on that surface

public:
    /// Default constructor
    ftPoint() {}
    /// Constructor
    /// \param pt Point
    /// \param sf Associated surface
    /// \param u First parameter in surface corresponding to point
    /// \param v Second parameter in surface corresponding to point
    ftPoint(Point pt, ftSurface* sf = 0, double u = 0, double v = 0)
	: pt_(pt), surface_(sf), u_(u), v_(v) {}
    /// Constructor. Point in 3d given by space coordinates
    ftPoint(double x, double y, double z)
	: pt_(x, y, z), surface_(0), u_(0), v_(0) {}
    /// Position of point
    const Point& position() const { return pt_; }
    /// Associated face/surface
    ftSurface* face() const { return surface_; }
    /// 1. parameter in surface
    double u() const { return u_; }
    /// 2. parameter in surface
    double v() const { return v_; }
    /// Normal of surface in point
    Point normal() const;
};

} // namespace Go


#endif // _FTPOINT_H

