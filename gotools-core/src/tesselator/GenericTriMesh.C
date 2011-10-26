#include "GoTools/tesselator/GenericTriMesh.h"
#include "GoTools/utils/errormacros.h"

namespace Go
{


//===========================================================================
GenericTriMesh::GenericTriMesh(int num_vert, int num_tri,
			       bool use_normals,
			       bool use_texcoords)
    : use_norm_(use_normals),
      use_texc_(use_texcoords)
//===========================================================================
{
    resize(num_vert, num_tri);
}


//===========================================================================
GenericTriMesh::~GenericTriMesh()
//===========================================================================
{
}


//===========================================================================
GenericTriMesh* GenericTriMesh::asGenericTriMesh()
//===========================================================================
{
    return this;
}


//===========================================================================
void GenericTriMesh::resize(int num_vert, int num_tri)
//===========================================================================
{
    vert_.resize(num_vert*3);
    param_.resize(num_vert*2);
    bd_.resize(num_vert);
    if (use_norm_)
	norm_.resize(num_vert*3);
    if (use_texc_)
	texc_.resize(num_vert*2);
    triangles_.resize(num_tri*3);
}

//===========================================================================
double* GenericTriMesh::vertexArray()
//===========================================================================
{
    if (vert_.empty()) {
	MESSAGE("Trying to get pointer to empty vector - returning NULL");
	return NULL;
    }
    return &vert_[0];
}


//===========================================================================
double* GenericTriMesh::paramArray()
//===========================================================================
{
    if (param_.empty()) {
	MESSAGE("Trying to get pointer to empty vector - returning NULL");
	return NULL;
    }
    return &param_[0];
}


//===========================================================================
int GenericTriMesh::atBoundary(int idx)
//===========================================================================
{
    if (bd_.empty()) {
	MESSAGE("Trying to get int element of empty vector - returning 0");
	return 0;
    }
    return bd_[idx];
}


//===========================================================================
int* GenericTriMesh::boundaryArray()
//===========================================================================
{
    if (bd_.empty()) {
	MESSAGE("Trying to get pointer to empty vector - returning NULL");
	return NULL;
    }
    return &bd_[0];
}


//===========================================================================
double* GenericTriMesh::normalArray()
//===========================================================================
{
    if (norm_.empty()) {
	MESSAGE("Trying to get pointer to empty vector - returning NULL");
	return NULL;
    }
    return &norm_[0];
}


//===========================================================================
double* GenericTriMesh::texcoordArray()
//===========================================================================
{
    if (texc_.empty()) {
	MESSAGE("Trying to get pointer to empty vector - returning NULL");
	return NULL;
    }
    return &texc_[0];
}


//===========================================================================
unsigned int* GenericTriMesh::triangleIndexArray()
//===========================================================================
{
    if (triangles_.empty()) {
	MESSAGE("Trying to get pointer to empty vector - returning NULL");
	return NULL;
    }
    return &triangles_[0];
}



} // namespace Go

