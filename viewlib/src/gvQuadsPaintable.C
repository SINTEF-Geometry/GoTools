//===========================================================================
//                                                                           
// File: gvQuadsPaintable.C                                                  
//                                                                           
// Created: Fri Jan 21 14:40:07 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvQuadsPaintable.C,v 1.1 2007-04-17 12:25:55 sbr Exp $
//                                                                           
//===========================================================================


#include "GoTools/viewlib/gvQuadsPaintable.h"
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

#include <iostream>
using namespace std;

//===========================================================================
gvQuadsPaintable::~gvQuadsPaintable()
//===========================================================================
{
}

//===========================================================================
void gvQuadsPaintable::paint(gvTexture* texture)
//===========================================================================
{
    // Deal with visibility and selection state
    if (!visible_) return;
    glDisable(GL_LIGHTING);
    if (selected_) {
	glColor4fv(selected_color_.rgba);
    } else {
	glColor4fv(normal_color_.rgba);
    }
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, quads_.vertexArray());
    int mode[2];
    glGetIntegerv(GL_POLYGON_MODE, mode);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawElements(GL_QUADS, quads_.numQuads()*4,
		   GL_UNSIGNED_INT, quads_.quadIndexArray());
    glPolygonMode(GL_FRONT, mode[0]);
    glPolygonMode(GL_BACK, mode[1]);
    glDisableClientState(GL_VERTEX_ARRAY);
    glEnable(GL_LIGHTING);
}


