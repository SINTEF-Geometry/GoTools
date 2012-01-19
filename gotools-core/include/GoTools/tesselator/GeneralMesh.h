//===========================================================================
//                                                                           
// File: GeneralMesh.h                                                      
//                                                                           
// Created: Mars 08
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================


#ifndef _GENERALMESH_H
#define _GENERALMESH_H

#include "GoTools/utils/config.h"

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

};

} // namespace Go




#endif // _GENERALMESH_H
