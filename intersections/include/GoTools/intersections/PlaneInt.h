//===========================================================================
//                                                                           
// File: PlaneInt.h                                                          
//                                                                           
// Created: Mon Jan 24 15:25:54 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: PlaneInt.h,v 1.6 2006-04-20 10:27:44 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
    std::shared_ptr<SplineSurface>
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

