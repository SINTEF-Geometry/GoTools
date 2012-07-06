//===========================================================================
//                                                                           
// File: RegularVolMesh.C                                                    
//                                                                           
// Created: Thu Jul  5 15:52:04 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

