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

