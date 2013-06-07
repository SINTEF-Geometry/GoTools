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
