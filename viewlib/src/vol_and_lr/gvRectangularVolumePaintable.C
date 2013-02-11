//===========================================================================
//                                                                           
// File: gvRectangularVolumePaintable.C                                      
//                                                                           
// Created: Thu Jul  5 15:24:10 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================





#include "GoTools/viewlib/vol_and_lr/gvRectangularVolumePaintable.h"

#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/viewlib/gvTexture.h"
#ifdef _MSC_VER
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif


//===========================================================================
gvRectangularVolumePaintable::~gvRectangularVolumePaintable()
//===========================================================================
{
}

//===========================================================================
void gvRectangularVolumePaintable::paint(gvTexture* texture)
//===========================================================================
{
    if (!visible_) return;

    // Draw surfaces a little bit behind other things 
    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1.0, 1.0);

    // Deal with visibility and selection state
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
    GLfloat grey[] = { 0.8f, 0.8f, 0.8f, 1.0 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);

    // Set up the vertex and normal arrays
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, tri_.vertexArray());
    glNormalPointer(GL_DOUBLE, 0, tri_.normalArray());
    bool use_textures = (texture != 0) && tri_.useTexCoords();
    if (use_textures) {
	glEnable(GL_TEXTURE_2D);
	texture->bind();
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glTexCoordPointer(2, GL_DOUBLE, 0, tri_.texcoordArray());
	if (selected_) {
	   glMaterialfv(GL_FRONT_AND_BACK,
			GL_AMBIENT_AND_DIFFUSE,
			grey);
	} else
	{
	   glMaterialfv(GL_FRONT_AND_BACK,
			GL_AMBIENT_AND_DIFFUSE,
			white);
	}
    }

    // Draw the triangle strips
    for (int i = 0; i < tri_.numStrips(); ++i) {
	glDrawElements(GL_TRIANGLE_STRIP, tri_.stripLength(),
		       GL_UNSIGNED_INT, tri_.stripArray()+i*tri_.stripLength());
    }
    if (use_textures) {
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisable(GL_TEXTURE_2D);
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_POLYGON_OFFSET_FILL);
}
