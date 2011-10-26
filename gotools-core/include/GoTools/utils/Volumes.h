//==========================================================================
//                                                                          
// File: Volumes.h                                                           
//                                                                          
// Created: Sun May 25 20:54:43 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Volumes.h,v 1.9 2006-04-19 09:16:24 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _VOLUMES_H
#define _VOLUMES_H

#include "GoTools/utils/Factorial.h"
#include "GoTools/utils/Array.h"
#include <vector>


namespace Go {


/// Calculates the determinant of a 2 x 2 matrix, represented in memory as a
/// sequence of 2 Array s of length 2.  Same function also exists for 3 x 3
/// matrices.
template< typename T >
inline T determinantOf(const Array<T, 2>* a) 
{
    return a[0][0] * a[1][1] - a[1][0] * a[0][1];
};


/// Calculates the determinant of a 3 x 3 matrix, represented in memory as a
/// sequence of 3 Array s of length 3.  Same function also exists for 2 x 2
/// matrices.
template<typename T>
inline T determinantOf(const Array<T, 3>* a) 
{
    return 
	a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) -
	a[0][1] * (a[1][0] * a[2][2] - a[2][0] * a[1][2]) +
	a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
};


/// Computes the volume of a simplex consisting of (Dim+1) vertices embedded
/// in Euclidean space of dimension (Dim)
template<typename T, int Dim>
inline T simplex_volume(const Array<T, Dim>* a)
{
    Array<T, Dim> tmp[Dim];
    for (int i = 0; i < Dim; ++i) {
	tmp[i] = a[i] - a[i+1];
    }
    return determinantOf(tmp) * InverseFactorial<double, Dim>::val();
    // determinant / factorial
}


/// Computes the area of a 2-dimensional triangle.  Input is an array of
/// corner points.  Same function also exists for 3-dimensional triangles.
template <typename T>
inline T area(const Array<T, 2>* c)
{ return simplex_volume(c); }


/// Computes the area of a 2-dimensional triangle.  Input is an array of
/// corner points.  Same function also exists for 2-dimensional triangles.
template < typename T >
inline T area(const Array<T, 3>* c)
{
    // Using the one-half cross product rule
    Array<T, 3> d0 = c[1] - c[0];
    Array<T, 3> d1 = c[2] - c[0];
    Array<T, 3> crossprod = d0.cross(d1);
    return 0.5 * crossprod.length();
}


/// Computes the volume of a 3D simplex (embedded i 3D space).
template <typename T>
inline T volume(const Array<T, 3>* c)
{ return simplex_volume(c); }


/// Computes the signed area of a triangle embedded in 3D space. Input is an
/// array of corner points and a normal to determine the sign.
template <typename T>
T signed_area(const Array<T, 3>* c, const Array<T, 3>& normal)
{
    // Using the one-half cross product rule
    Array<T, 3> d0 = c[1] - c[0];
    Array<T, 3> d1 = c[2] - c[0];
    Array<T, 3> crossprod = d0.cross(d1);
    if (crossprod*normal > 0) {
	return 0.5 * crossprod.length();
    } else {
	return -0.5 * crossprod.length();
    }
}


} // namespace Go


#endif // _VOLUMES_H

