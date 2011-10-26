//===========================================================================
//                                                                           
// File: gvGenericTriQuadPaintable.h                                         
//                                                                           
// Created: Fri Mar  1 12:22:22 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvGenericTriQuadPaintable.h,v 1.1 2007-04-17 12:25:36 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVGENERICTRIQUADPAINTABLE_H
#define _GVGENERICTRIQUADPAINTABLE_H


#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/viewlib/gvGenericTriMesh.h"
#include "GoTools/viewlib/gvGenericTriQuadMesh.h"
#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif
#include <GL/gl.h>

/** Documentation ...
    etc
 */

template <typename FloatType>
class gvGenericTriQuadPaintable : public gvPaintable
{
public:
    gvGenericTriQuadPaintable(gvGenericTriQuadMesh<FloatType>& mesh,
			  const gvColor& ncolor,
			  const gvColor& scolor,
			  int id)
	: gvPaintable(ncolor, scolor, id),
	  mesh_(mesh)
    {}
    gvGenericTriQuadPaintable(gvGenericTriQuadMesh<FloatType>& mesh,
			  const gvColor& ncolor,
			  int id)
	: gvPaintable(ncolor, id),
	  mesh_(mesh)
    {}
    virtual ~gvGenericTriQuadPaintable()
    {}

    virtual void paint()
    {
	GLenum fptype = GLfptype(FloatType());
	// Draw surfaces a little bit behind other things 
	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1.0, 1.0);

	// Deal with visibility and selection state
	if (!visible_) return;
	if (selected_) {
	    glMaterialfv(GL_FRONT_AND_BACK,
			 GL_AMBIENT_AND_DIFFUSE,
			 selected_color_.rgba);
	} else {
	    glMaterialfv(GL_FRONT_AND_BACK,
			 GL_AMBIENT_AND_DIFFUSE,
			 normal_color_.rgba);
	}
	GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);

    // Set up the vertex and normal arrays
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glVertexPointer(3, fptype, 0, mesh_.vertexArray());
	//      std::cout << mesh_.vertexArray()[300] << ' '
	//  	      << mesh_.vertexArray()[301] << ' '
	//  	      << mesh_.vertexArray()[302] << std::endl;
	// @@@ We should check for the use of normals and textures in mesh_.
	glNormalPointer(fptype, 0, mesh_.normalArray());

    // Draw the triangles
	glDrawElements(GL_TRIANGLES, mesh_.numTriangles()*3,
		       GL_UNSIGNED_INT, mesh_.triangleIndexArray());
	glDrawElements(GL_QUADS, mesh_.numQuads()*3,
		       GL_UNSIGNED_INT, mesh_.quadIndexArray());

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisable(GL_POLYGON_OFFSET_FILL);
    }


protected:
    gvGenericTriQuadMesh<FloatType>& mesh_;
};





#endif // _GVGENERICTRIQUADPAINTABLE_H

