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

#ifndef _SPHEREINT_H
#define _SPHEREINT_H


#include "GoTools/intersections/AlgObj3DInt.h"
#include "GoTools/utils/Point.h"


namespace Go {


class SplineSurface;


/// Class representing spherical algebraic intersection objects.

class SphereInt : public AlgObj3DInt {
public:
    /// Constructor.
    /// Used when reading from file.
    SphereInt();

    /// Constructor.
    /// \param center the center point of the sphere.
    /// \param radius the radius of the sphere.
    SphereInt(Point center, double radius);

    /// Destructor.
    virtual ~SphereInt();

    /// Read a sphere description from file.
    /// \param is the stream containing the sphere description.
    void read(std::istream& is);

    /// Get center point of the sphere.
    /// \return The center point of the sphere.
    Point center() const;

    /// Get the radius of the sphere.
    /// \return The radius of the sphere.
    double radius() const;

    /// Get a SplineSuface representing the sphere. This can be used to
    /// visualize the object.
    /// \return The spline surface of the sphere.
    SplineSurface* surface() const;

private:

    Point center_; // 3D ref point (x_0, y_0, z_0).
    double radius_; // (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = radius_^2

};


} // namespace Go


#endif // _SPHEREINT_H

