//===========================================================================
//                                                                           
// File: CylinderInt.h                                                       
//                                                                           
// Created: Wed Apr  6 14:01:08 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: CylinderInt.h,v 1.5 2006-03-08 09:31:19 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#ifndef _CYLINDERINT_H
#define _CYLINDERINT_H


#include "GoTools/intersections/AlgObj3DInt.h"
#include "GoTools/utils/Point.h"


namespace Go {


class SplineSurface;


/// Class for cylindrical algebraic intersection objects.

class CylinderInt : public AlgObj3DInt
{
public:
    /// Constructor.
    /// Used when reading from file.
    CylinderInt();

    /// Constructor.
    /// \param ax_pt reference point on the center axis of the
    /// cylinder.
    /// \param ax_dir direction vector of the center axis of the
    /// cylinder.
    /// \param radius the radius of the cylinder.
    CylinderInt(Point ax_pt, Point ax_dir, double radius);

    /// Destructor.
    virtual ~CylinderInt();

    /// Read a cylinder description from file.
    /// \param is the stream containing the cylinder description.
    void read(std::istream& is);

    /// Get the axis reference point.
    /// \return the axis reference point.
    Point ax_pt() const;
    /// Get the axis direction vector..
    /// \return the axis direction vector.
    Point ax_dir() const;
    /// Get the radius of the cylinder.
    /// \return the radius of the cylinder.
    double radius() const;

    /// Get a SplineSurface representing the cylinder. This can be
    /// used to visualize the object.  Since the object is to be
    /// tesselated we make sure it is bounded.
    /// \param bottom_pos the lowest point on the direction axis.
    /// \param height the height of the cylinder.
    /// \return The spline surface of the cylinder.
    SplineSurface* surface(Point bottom_pos, double height) const;


private:

    Point ax_pt_; // 3D ref point (x_0, y_0, z_0).
    Point ax_dir_; // 3D ref point (x_0, y_0, z_0).
    double radius_; // (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = radius_^2

    // Setting the values for the implicit surface representation
    // using ac_pt_, ax_dir_ & radius_.
    void setImplicitValues();

};


} // namespace Go


#endif // _CYLINDERINT_H

