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

     // /// Translate all vertices by vert_translation, wrt the geometry.
     // void translate(const std::vector<double>& vert_translation);

private:
    bool use_norm_;
    bool use_texc_;
  
    std::vector<double> vert_;
    std::vector<double> param_;
    std::vector<int> bd_;
    std::vector<double> norm_;
    std::vector<double> texc_;
    std::vector<unsigned int> triangles_;

    // /// The 3D-translation wrt the geometry.
    // std::vector<double> vert_translation_;
  
};

} // namespace Go



#endif // GENERICTRIMESH_H

