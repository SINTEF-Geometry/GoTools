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

#ifndef _REGULARMESH_H
#define _REGULARMESH_H

#include "GoTools/tesselator/GeneralMesh.h"
#include <vector>

namespace Go
{


/** Documentation ...
 */
class Triangle
{
public:
    explicit Triangle(int index) : index_(index) {}
    void setIndex(int i) { index_ = i; }
    int index() const { return index_; }
private:
    int index_;
};

/** Simple triangulation class, designed to store data for
 * GL visualizing with triangle strips.
 * Vertices are stored in x1 y1 z1 x2 y2 z2... format.
 * Triangles are stored as startindices into the stripe array.
 * Based on earlier work.
 */
class GO_API RegularMesh : public GeneralMesh
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
     RegularMesh(int m = 20, int n = 20,
		  bool use_normals = true,
		bool use_texcoords = false);

    /// Destructor
    virtual ~RegularMesh();

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
    void resize(int m, int n);


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

     /// Translate all vertices by vert_translation, wrt the geometry.
     void translate(const std::vector<double>& vert_translation);

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

    /// The 3D-translation wrt the geometry.
    std::vector<double> vert_translation_;
  
    ///
    /// Creates strips and triangles
    ///
    void calculateIndices()
    {
	///
	/// Ok, remember that 'strip_length_' is the number of vertices used by
	/// a strip. 
	///
	strips_.resize(num_strips_*strip_length_);
	int last = 0;
	///
	/// Compute the indices (triangle strips) on a rectangular grid.
	///
	/// ('strips_' contains the indices (three of them) making up each
	/// of the triangles, this doesn't have much to do with strips, except
	/// that they are later grouped into 'num_strips_' equally large
	/// chunks, which then make up 'num_strips_' strips.
	///
	int m = strip_length_/2; // number of vertices along strip
	int i;
	for ( i = 0; i < num_strips_; ++i ) {
	    for ( int j = 0; j < m; ++j ) {
		strips_[last++] = (i+1)*m+j ;
		strips_[last++] = i*m+j;
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
	    for (int j = 0; j < strip_length_-2; ++j) {
		triangles_.push_back(Triangle(i*strip_length_ + j));
		for (int h=0; h<3; ++h)
		    triangle_index_.push_back(strips_[i*strip_length_ + j + h]);
	    }
	}
    }
  
};

} // namespace Go


#endif // _REGULARMESH_H

