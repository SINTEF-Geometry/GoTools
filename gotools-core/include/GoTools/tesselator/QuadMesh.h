//===========================================================================
//                                                                           
// File: QuadMesh.h                                                        
//                                                                           
// Created: Fri Jan 21 14:09:46 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: 
//                                                                           
//===========================================================================

#ifndef _QUADMESH_H
#define _QUADMESH_H

#include "GoTools/tesselator/GeneralMesh.h"
#include <vector>

namespace Go
{

/** Based on gvGenericTriMesh,
 *  designed to store data for GL visualizing.
 */
class QuadMesh : public GeneralMesh
{
public:
  /// Constructor. Give mesh size. Set whether or not normals
  /// and texture coordinates should be computed
    QuadMesh(int num_vert = 0, int num_quads = 0,
	       bool use_normals = true,
	       bool use_texcoords = false)
	: use_norm_(use_normals),
	  use_texc_(use_texcoords)
    {
	resize(num_vert, num_quads);
    }

    /// Check if normals will be computed
     bool useNormals()
    { return use_norm_; }

    /// Check if texture coordinates will be computed
    bool useTexCoords() 
    { return use_texc_; }

    /// The total number of quads
    int numQuads()
    { return (int)quads_.size()/4; }

     /// The total number of triangles. Not used in this context
   virtual int numTriangles()
    { return (int)quads_.size()/4; }

    virtual int numVertices()
    { return (int)vert_.size()/3; }

    /// Change mesh size
    void resize(int num_vert, int num_quads)
    {
	vert_.resize(num_vert*3);
	if (use_norm_)
	    norm_.resize(num_vert*3);
	if (use_texc_)
	    texc_.resize(num_vert*2);
	quads_.resize(num_quads*4);
    }

    virtual double* vertexArray() { return &vert_[0]; }
    /// Normal information
    double* normalArray() { return &norm_[0]; }
    /// Texture coordinate information
    double* texcoordArray() { return &texc_[0]; }
    /// Indices to quads
    unsigned int* quadIndexArray() { return &quads_[0]; }
    /// Indices to triangle. Not used in this context
    virtual unsigned int* triangleIndexArray() { return &quads_[0]; }
    virtual double* paramArray() {return 0;}
    virtual int atBoundary(int idx) {return 0;}


private:
    bool use_norm_;
    bool use_texc_;
  
    std::vector<double> vert_;
    std::vector<double> norm_;
    std::vector<double> texc_;
    std::vector<unsigned int> quads_;
  
};

} // namespace Go



#endif // _GVQUADMESH_H

