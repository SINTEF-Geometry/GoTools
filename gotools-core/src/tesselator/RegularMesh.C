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

