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

#ifndef _PLANEINT_H
#define _PLANEINT_H


#include "GoTools/intersections/AlgObj3DInt.h"
#include "GoTools/utils/Point.h"
#include <memory>


namespace Go {


class SplineSurface;


/// Class representing planar algebraic intersection objects.

class PlaneInt : public AlgObj3DInt
{
public:
    /// Constructor.
    /// Used when reading from file.
    PlaneInt();

    /// Constructor.
    /// \param point reference point in the plane.
    /// \param normal normal to the plane.
    PlaneInt(Point point, Point normal);

    /// Constructor.
    /// The plane is described by the expression \f$ax + by + cz + d = 0\f$.
    /// \param a the x multiplicator.
    /// \param b the y multiplicator.
    /// \param c the z multiplicator.
    /// \param d the constant.
    PlaneInt(double a, double b, double c, double d);

    /// Constructor.
    virtual ~PlaneInt();

    /// Read a plane description from file.
    /// \param is the stream containing the plane description.
    void read(std::istream& is);

    /// Get the x multiplicator.
    /// \return a The x multiplicator.
    double a() const;
    /// Get the y multiplicator.
    /// \return b The y multiplicator.
    double b() const;
    /// Get the z multiplicator.
    /// \return c The z multiplicator.
    double c() const;
    /// Get the constant.
    /// \return d The constant.
    double d() const;

    /// Get a SplineSuface representing the plane. This can be used to
    /// visualize the object. Since the object is to be tesselated we
    /// make sure it is bounded.
    /// \param mid_pt the intersection of the diagonals in a rectangle
    /// in the plane.
    /// \param length_x the length of the rectangle in the x direction.
    /// \param length_y the length of the rectangle in the y direction.
    /// \return The spline surface of the plane.
    shared_ptr<SplineSurface>
    surface(Point mid_pt, double length_x, double length_y) const;


private:

    // (n_x, n_y, n_z)*(x - x_0, y - y_0, z - z_0) = 0.
    Point point_; // 3D ref point (x_0, y_0, z_0).
    Point normal_; // 3D normal direction (n_x, n_y, n_z).

    // Compute point_ & normal_.
    void computePoints();

    // Compute a, b, c & d based on point_ & normal_.
    void computeConstants(Point point, Point normal,
			  double& a, double& b, double& c, double& d) const;

    // We return the projection of the 3D input space_pt.
    Point projectPoint(const Point& space_pt) const;

};


} // namespace Go


#endif // _PLANEINT_H

