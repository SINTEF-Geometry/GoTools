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

#ifndef _GENERALMESH_H
#define _GENERALMESH_H

#include "GoTools/utils/config.h"

#include <vector>

namespace Go
{

class RegularMesh;
 class LineStrip;
 class GenericTriMesh;

/** Super class for meshes. Used in transfer of mesh.
*/


class GO_API GeneralMesh
{
public:

    GeneralMesh();

    /// Destructor
    virtual ~GeneralMesh();

    /// Number of nodes in mesh
    virtual int numVertices() = 0;

    /// Fetch all nodes 
    virtual double* vertexArray() = 0;

    /// Fetch parameter values corresponding to the nodes
    virtual double* paramArray() = 0;

    /// Check if a given node lies at the boundary
    virtual int atBoundary(int idx) = 0;
    
    /// The number of triangles represented by this mesh
    virtual int numTriangles()
    {
      return 0;
    }

    /// Indices of nodes belonging to each triangle
    virtual unsigned int* triangleIndexArray() = 0;

    /// Casting. Return as regular mesh if possible
    virtual RegularMesh* asRegularMesh();

    /// Casting. Return as line strip if possible
    virtual LineStrip* asLineStrip();

    /// Casting. Return as generic tri mesh if possible
    virtual GenericTriMesh* asGenericTriMesh();

    /// Translate all vertices by vert_translation, wrt the geometry, discarding any previus translation.
    void translate(const std::vector<double>& vert_translation);

    /// The 3D-translation wrt the geometry.
    std::vector<double> vert_translation_;


};

} // namespace Go




#endif // _GENERALMESH_H
