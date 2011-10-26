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

/** Super class for meshes. Used in transfer of mesh.
*/

class RegularMesh;
 class LineStrip;
 class GenericTriMesh;


class GO_API GeneralMesh
{
public:
    virtual ~GeneralMesh();

    virtual int numVertices() = 0;

    virtual double* vertexArray() = 0;

    virtual double* paramArray() = 0;

    virtual int atBoundary(int idx) = 0;
    
    virtual int numTriangles()
    {
      return 0;
    }

    virtual unsigned int* triangleIndexArray() = 0;

    virtual RegularMesh* asRegularMesh();

    virtual LineStrip* asLineStrip();

    virtual GenericTriMesh* asGenericTriMesh();

};

} // namespace Go




#endif // _GENERALMESH_H
