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

#include "GoTools/tesselator/RegularMesh.h"

namespace Go
{


//===========================================================================
    RegularMesh::RegularMesh(int m, int n,
			     bool use_normals,
			     bool use_texcoords)
	: use_norm_(use_normals),
	  use_texc_(use_texcoords)
    //===========================================================================
    {
	resize(m, n);
    }


    RegularMesh::~RegularMesh()
    {
    }

    RegularMesh* RegularMesh::asRegularMesh()
	{
	    return this;
	}

    void RegularMesh::resize(int m, int n)
    {
//    	std::cout << m << ' ' << n << std::endl;
	num_strips_ = n-1;
	strip_length_ = 2*m;
	int N = m*n; // total number of vertices
	vert_.resize(N*3);
	param_.resize(N*2);
	if (use_norm_)
	    norm_.resize(N*3);
	if (use_texc_)
	    texc_.resize(N*2);
	calculateIndices();
    }

    int RegularMesh::atBoundary(int idx)
    {
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

} // namespace Go

