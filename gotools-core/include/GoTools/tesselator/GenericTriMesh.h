//===========================================================================
//                                                                           
// File: GenericTriMesh.h                                                  
//                                                                           
// Created: Fri Nov 30 13:33:46 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: GenericTriMesh.h,v 1.5 2009-05-13 07:32:50 vsk Exp $
//                                                                           
// Description: Based on RegularMesh, with lots of stuff simplified,
//              designed to store data for GL visualizing (regular triangles,
//              not triangle strips).
//                                                                           
//===========================================================================

#ifndef GENERICTRIMESH_H
#define GENERICTRIMESH_H


#include "GoTools/tesselator/GeneralMesh.h"
#include <vector>
/* #ifdef _MSC_VER */
/* #include <windows.h> */
/* #endif */
/* #include <GL/gl.h> */

namespace Go
{

/** Based on RegularMesh, 
 * designed to store data for GL visualizing (regular triangles,
 * not triangle strips).
 */
class GO_API GenericTriMesh : public GeneralMesh
{
public:
  /// Constructor. Give mesh size. Set whether or not normals
  /// and texture coordinates should be computed
    GenericTriMesh(int num_vert = 0, int num_tri = 0,
		   bool use_normals = true,
		   bool use_texcoords = false);

    /// Destructor
    ~GenericTriMesh();

    /// Check if normals will be computed
    bool useNormals()
    { return use_norm_; }

    /// Check if texture coordinates will be computed
    bool useTexCoords() 
    { return use_texc_; }

    /// The total number of triangles
    virtual int numTriangles()
    { return (int)triangles_.size()/3; }

    virtual int numVertices()
    { return (int)vert_.size()/3; }

    /// Change mesh size
    void resize(int num_vert, int num_tri);

    virtual double* vertexArray();
    virtual double* paramArray();
    virtual int atBoundary(int idx);

    /// Indices to boundary nodes
    int* boundaryArray();
    /// Normal information
    double* normalArray();
    /// Texture coordinate information
    double* texcoordArray();
    virtual unsigned int* triangleIndexArray();

    /// Casting. Return as generic tri mesh 
     virtual GenericTriMesh* asGenericTriMesh();

private:
    bool use_norm_;
    bool use_texc_;
  
    std::vector<double> vert_;
    std::vector<double> param_;
    std::vector<int> bd_;
    std::vector<double> norm_;
    std::vector<double> texc_;
    std::vector<unsigned int> triangles_;
  
};

} // namespace Go



#endif // GENERICTRIMESH_H

