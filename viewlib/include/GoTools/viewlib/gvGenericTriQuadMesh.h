//===========================================================================
//                                                                           
// File: gvGenericTriQuadMesh.h                                              
//                                                                           
// Created: Fri Mar  1 12:01:18 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvGenericTriQuadMesh.h,v 1.1 2007-04-17 12:25:36 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVGENERICTRIQUADMESH_H
#define _GVGENERICTRIQUADMESH_H


#include <vector>


/** Structure for ...
 */

template <typename FloatType>
class gvGenericTriQuadMesh
{
public:
    gvGenericTriQuadMesh(int num_vert, int num_tri, int num_quads,
			 bool use_normals = true,
			 bool use_texcoords = false)
	: use_norm_(use_normals),
	  use_texc_(use_texcoords)
    {
	resize(num_vert, num_tri, num_quads);
    }

    bool useNormals()
    { return use_norm_; }

    bool useTexCoords() 
    { return use_texc_; }

    /// The total number of triangles
    int numTriangles()
    { return triangles_.size()/3; }
    /// The total number of quads
    int numQuads()
    { return quads_.size()/4; }

    int numVertices()
    { return vert_.size()/3; }

    void resize(int num_vert, int num_tri, int num_quads)
    {
	vert_.resize(num_vert*3);
	if (use_norm_)
	    norm_.resize(num_vert*3);
	if (use_texc_)
	    texc_.resize(num_vert*3);
	triangles_.resize(num_tri*3);
	quads_.resize(num_quads*4);
    }

    FloatType* vertexArray() { return &vert_[0]; }
    FloatType* normalArray() { return &norm_[0]; }
    FloatType* texcoordArray() { return &texc_[0]; }
    int* triangleIndexArray() { return &triangles_[0]; }
    int* quadIndexArray() { return &triangles_[0]; }

private:
    bool use_norm_;
    bool use_texc_;
  
    std::vector<FloatType> vert_;
    std::vector<FloatType> norm_;
    std::vector<FloatType> texc_;
    std::vector<int> triangles_;
    std::vector<int> quads_;
  
};





#endif // _GVGENERICTRIQUADMESH_H

