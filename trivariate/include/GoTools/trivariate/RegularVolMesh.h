//===========================================================================
//                                                                           
// File: RegularVolMesh.h                                                    
//                                                                           
// Created: Thu Jul  5 15:40:17 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _REGULARVOLMESH_H
#define _REGULARVOLMESH_H

#include "GoTools/tesselator/RegularMesh.h"


namespace Go
{


/** Simple triangulation class, designed to store data for
 * GL visualizing with triangle strips.
 * Vertices are stored in x1 y1 z1 x2 y2 z2... format.
 * Triangles are stored as startindices into the stripe array.
 * Based on earlier work.
 */
class GO_API RegularVolMesh : public GeneralMesh
{
public:
  /// Constructor. Give mesh size. Set whether or not normals
  /// and texture coordinates should be computed
    ///
    /// 010613: To me it seems that 'm' and 'n' is the dimension of the
    ///         (rectangular) grid of vertices, i.e., that in total there
    ///         will be (m-1)*(n-1)*2 triangles. Trying to fix things using
    ///         this assumption... (jon)
    ///
    RegularVolMesh(int m = 20,
		   bool use_normals = true,
		   bool use_texcoords = false);

    /// Destructor
    virtual ~RegularVolMesh();

    /// Check if normals will be computed
    bool useNormals()
    { return use_norm_; }

    /// Check if texture coordinates will be computed
   bool useTexCoords() 
    { return use_texc_; }
  
    /// Number of strips
    int numStrips()
    { return num_strips_; }

    /// The number of vertices in a strip
    /// (i.e the number of triangles in a strip+ 2).
    int stripLength()
    { return strip_length_; }

    /// The total number of triangles
    virtual int numTriangles()
    { return (int)triangles_.size(); }

    /// Number of vertices 
    virtual int numVertices()
    { return (int)vert_.size()/3; }

    /// Change mesh size
    void resize(int m);


    virtual double* vertexArray() { return &vert_[0]; }
    virtual double* paramArray() { return &param_[0]; }
    virtual int atBoundary(int idx);
    /// Normal information
    double* normalArray() { return &norm_[0]; }
    /// Texture coordinate information
    double* texcoordArray() { return &texc_[0]; }
    unsigned int* stripArray() { return &strips_[0]; }
    Triangle* triangleArray() { return &triangles_[0]; }
    virtual unsigned int* triangleIndexArray() { return &triangle_index_[0]; }
  
    /// Casting. Return as regular mesh
    virtual RegularMesh* asRegularMesh();


private:
    bool use_norm_;
    bool use_texc_;
    int num_strips_;
    int strip_length_;
  
    typedef std::vector<double> Vd;
    Vd vert_;
    Vd param_;
    Vd norm_;
    Vd texc_;
    std::vector<unsigned int> strips_;
    std::vector<Triangle> triangles_;
    std::vector<unsigned int> triangle_index_;
  
    ///
    /// Creates strips and triangles
    ///
    void calculateIndices();
  
};

} // namespace Go

#endif // _REGULARVOLMESH_H


