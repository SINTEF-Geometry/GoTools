//===========================================================================
//                                                                           
// File: BaryCoordSystemTriangle3D.h                                         
//                                                                           
// Created: Mon Sep 27 13:57:02 2004                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: BaryCoordSystemTriangle3D.h,v 1.3 2005-06-06 09:32:01 oan Exp $
//                                                                           
//===========================================================================

#ifndef _BARYCOORDSYSTEMTRIANGLE3D_H
#define _BARYCOORDSYSTEMTRIANGLE3D_H


#include "GoTools/utils/Array.h"
#include "GoTools/utils/Volumes.h"
#include <iostream>
#include <algorithm>


namespace Go {

    /** A barycentric coordinate system for a triangle (2-manifold) embedded in 3D.
     *  Note that this differs from what can be expressed using the 
     *  BaryCoordSystem template in that for the latter, the dimension of the simplex
     *  and the dimension of the space is equal.
     */

class BaryCoordSystemTriangle3D {
public:
    /// Empty default constructor
    BaryCoordSystemTriangle3D() { }

    /// Constructor. Takes an array of points that will become
    /// the corners of the coordinate simplex.
    BaryCoordSystemTriangle3D(const Array<double, 3>* corners)
    {
	std::copy(corners, corners+3, corners_);
	total_area_ = area(corners);
	normal_ = (corners[1]-corners[0]) % (corners[2]-corners[0]);
	normal_.normalize();
    }

    /// Input is a barycentric point, output is the corresponding
    /// cartesian point.
    template <typename T>
    Array<T, 3> baryToCart(const Array<T, 3>& bary_pt) const
    {
        Array<T, 3> cart_pt(T(0.0), T(0.0), T(0.0));
	for (int i = 0; i < 3; ++i) {
	    cart_pt[0] += corners_[i][0] * bary_pt[i];
	    cart_pt[1] += corners_[i][1] * bary_pt[i];
	    cart_pt[2] += corners_[i][2] * bary_pt[i];
	}
	return cart_pt;
    }

    /// Input is a cartesian point, output is the corresponding
    /// barycentric point.
    template <typename T>
    Array<T, 3> cartToBary(const Array<T, 3>& cart_pt) const
    {
        static Array<T, 3> subtriangle[3];
	int i;
	for (i = 1; i < 3; ++i) {
	    subtriangle[i] = corners_[i];
	}

	Array<T, 3> bary_pt;
	for (i = 0; i < 3; ++i) {
	    subtriangle[i] = cart_pt;
	    bary_pt[i] = signed_area(subtriangle, normal_);
	    subtriangle[i] = corners_[i];
	}
	return bary_pt / total_area_;
    }

private:
    Array<double, 3> corners_[3];
    Array<double, 3> normal_;
    double total_area_;
};


} // namespace Go

#endif // _BARYCOORDSYSTEMTRIANGLE3D_H

