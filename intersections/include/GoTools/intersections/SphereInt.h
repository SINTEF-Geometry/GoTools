//===========================================================================
//                                                                           
// File: SphereInt.h                                                         
//                                                                           
// Created: Mon Jan 31 13:19:38 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: SphereInt.h,v 1.7 2006-03-08 09:31:19 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

