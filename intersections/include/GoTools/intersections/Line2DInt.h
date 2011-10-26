//===========================================================================
//                                                                           
// File: Line2DInt.h                                                           
//                                                                           
// Created: Mon Jan 24 11:24:54 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Line2DInt.h,v 1.3 2006-02-22 14:52:04 jbt Exp $
//                                                                           
// Description: Implementation of (algebraic) 2d-line.
//                                                                           
//===========================================================================


#ifndef _LINEINT_H
#define _LINEINT_H


#include "GoTools/intersections/AlgObj2DInt.h"
#include "GoTools/utils/Point.h"


namespace Go {


/// Class representing an algebraic line in 2-dimensional space.

class Line2DInt : public AlgObj2DInt {
public:
    /// Constructor.
    /// \param point reference point which lies on the line.
    /// \param dir direction of the line.
    Line2DInt(Point point, Point dir);

    /// Constructor.
    /// The line is given by the expression \f$ax + by + c = 0\f$.
    /// \param a the x multiplicator.
    /// \param b the y multiplicator.
    /// \param c the constant
    Line2DInt(double a, double b, double c);

    /// Destructor.
    virtual ~Line2DInt(){};

    /// Get the x multiplicator.
    /// \return The x multiplicator.
    double a();

    /// Get the y multiplicator.
    /// \return The y multiplicator.
    double b();

    /// Get the constant.
    /// \return The constant.
    double c();

private:

    Point point_; // 2D ref point.
    Point dir_; // 2D line dir.

    // Compute point_ & dir_ based on a, b & c (as given by factors_).
    void computePoints();

    // Compute a, b & c based on point_ & dir_.
    void computeConstants(Point point, Point dir,
			  double& a, double& b, double& c);

};


} // namespace Go


#endif // _LINEINT_H
