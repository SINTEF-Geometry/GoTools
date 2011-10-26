//===========================================================================
//                                                                           
// File: CoordinateSystem.h                                                  
//                                                                           
// Created: Tue May 21 18:28:46 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CoordinateSystem.h,v 1.4 2005-06-06 09:32:01 oan Exp $
//                                                                           
// Description: Defines a Cartesian coordinate system
//                                                                           
//===========================================================================

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

