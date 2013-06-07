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

#ifndef _COORDINATESYSTEM_H
#define _COORDINATESYSTEM_H

#include "GoTools/utils/Array.h"
#include "GoTools/utils/MatrixXD.h"

namespace Go
{

    /// Defines a Cartesian coordinate system

template <int Dim>
class CoordinateSystem
{
public:
    typedef MatrixXD<double, Dim> Matrix;
    typedef Array<double, Dim> Vector;

    /// The default constructor creates an object with 
    /// an identity rotation matrix and zero translation vector.
    CoordinateSystem()
    {
	rotation_.identity();
	translation_.zero();
    }

    /// Creates a CoordinateSystem with the specified 
    /// rotation matrix and translation vector.
    CoordinateSystem(const Matrix& rot, const Vector& trans)
	: rotation_(rot), translation_(trans)
    {}

    /// Get the 'global' coordinates of a vector 'v' expressed
    /// in the CoordinateSystem.
    Vector operator* (const Vector& v)
    {
	return rotation_*v + translation_;
    }

    /// Create a CoordinateSystem that is a composition of 'this' 
    /// CoordinateSystem and 'c'.
    CoordinateSystem operator* (const CoordinateSystem& c)
    {
	Matrix newrot = rotation_*c.rotation_;
	Vector newtrans = translation_ + rotation_*c.translation_;
	return CoordinateSystem(newrot, newtrans);
    }

    /// Get the rotation matrix.
    Matrix& rot() { return rotation_; }

    /// Get the rotation matrix. Const version.
    const Matrix& rot() const { return rotation_; }

    /// Get the translation vector.
    Vector& tr() { return translation_; }

    /// Get the translation vector. Const version.
    const Vector& tr() const { return translation_; }

private:
    Matrix rotation_;
    Vector translation_;
};


} // namespace Go



#endif // _COORDINATESYSTEM_H

