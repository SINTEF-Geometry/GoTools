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

#include "GoTools/trivariate/RegularVolMesh.h"
#include "GoTools/utils/errormacros.h"

namespace Go
{


//===========================================================================
    RegularVolMesh::RegularVolMesh(int m,
			     bool use_normals,
			     bool use_texcoords)
	: use_norm_(use_normals),
	  use_texc_(use_texcoords)
    //===========================================================================
    {
	resize(m);
    }


    RegularVolMesh::~RegularVolMesh()
    {
    }

    RegularMesh* RegularVolMesh::asRegularMesh()
	{
	    return NULL;
	}

    void RegularVolMesh::resize(int m)
    {
//    	std::cout << m << ' ' << n << std::endl;
	num_strips_ = 6*(m-1);
	strip_length_ = 2*m;
	int N = 6*m*m; // total number of vertices
	vert_.resize(N*3);
	param_.resize(N*2);
	if (use_norm_)
	    norm_.resize(N*3);
	if (use_texc_)
	    texc_.resize(N*2);
	calculateIndices();
    }

    int RegularVolMesh::atBoundary(int idx)
    {
	THROW("Not implemented yet.");

	int m = num_strips_ + 1;
	int n = strip_length_/2;
	if (idx < n)
	    return 1;  // Lower boundary
	if (idx >= n*(m-1))
	    return 1;  // Upper boundary
	if (idx%n == 0)
	    return 1;  // Left boundary
	if (idx%n == n-1)
	    return 1;  // Right boundary

	return 0;  // The node lies in the inner of the surface
    }


    void RegularVolMesh::calculateIndices()
    {
	///
	/// Ok, remember that 'strip_length_' is the number of vertices used by
	/// a strip. 
	///
	strips_.resize(num_strips_*strip_length_);
	///
	/// Compute the indices (triangle strips) on a rectangular grid.
	///
	/// ('strips_' contains the indices (three of them) making up each
	/// of the triangles, this doesn't have much to do with strips, except
	/// that they are later grouped into 'num_strips_' equally large
	/// chunks, which then make up 'num_strips_' strips.
	///
	/// Since we are dealing with surfaces we must split the vertices according to sf.
	int m = strip_length_/2; // number of vertices along strip
	int i, j, ki;
	for (ki = 0; ki < 6; ++ki)
	{
	    int last = ki*m*(m-1)*2;
	    for ( i = 0; i < num_strips_/6; ++i ) {
		for ( j = 0; j < m; ++j ) {
		    strips_[last++] = (i+1)*m+j +ki*m*m;
		    strips_[last++] = i*m+j + ki*m*m;
		}
	    }
	}

	///
	/// 'triangles_' doesn't have anything to do with strips, it is rather
	/// an array containing 'Triangles'.
	///
	triangles_.reserve(num_strips_*(strip_length_-2));
	triangles_.clear();
	triangle_index_.reserve(3*(num_strips_*(strip_length_-2)));
	for (i = 0; i < num_strips_; ++i) {
	    for (j = 0; j < strip_length_-2; ++j) {
		triangles_.push_back(Triangle(i*strip_length_ + j));
		for (int h=0; h<3; ++h)
		    triangle_index_.push_back(strips_[i*strip_length_ + j + h]);
	    }
	}
    }

} // namespace Go

