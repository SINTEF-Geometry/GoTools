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

