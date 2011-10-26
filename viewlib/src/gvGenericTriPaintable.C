//===========================================================================
//                                                                           
// File: gvGenericTriPaintable.C                                             
//                                                                           
// Created: Fri Nov 30 14:26:43 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvGenericTriPaintable.C,v 1.1 2007-04-17 12:25:51 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================
#if 0

#include "GoTools/viewlib/gvGenericTriPaintable.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif
#include <GL/gl.h>


//===========================================================================
gvGenericTriPaintable::~gvGenericTriPaintable()
//===========================================================================
{
}

//===========================================================================
void gvGenericTriPaintable::paint()
//===========================================================================
{
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
    glVertexPointer(3, GL_DOUBLE, 0, tri_.vertexArray());
//      std::cout << tri_.vertexArray()[300] << ' '
//  	      << tri_.vertexArray()[301] << ' '
//  	      << tri_.vertexArray()[302] << std::endl;
    // @@@ We should check for the use of normals and textures in tri_.
    glNormalPointer(GL_DOUBLE, 0, tri_.normalArray());

    // Draw the triangles
    glDrawElements(GL_TRIANGLES, tri_.numTriangles()*3,
		   GL_UNSIGNED_INT, tri_.triangleIndexArray());

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_POLYGON_OFFSET_FILL);
}


#endif
