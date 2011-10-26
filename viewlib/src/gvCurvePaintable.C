//===========================================================================
//                                                                           
// File: gvCurvePaintable.C                                                  
//                                                                           
// Created: Thu Jun 21 16:29:41 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvCurvePaintable.C,v 1.1 2007-04-17 12:25:50 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/gvCurvePaintable.h"
#include "GoTools/tesselator/LineStrip.h"
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
gvCurvePaintable::~gvCurvePaintable()
//===========================================================================
{
}

//===========================================================================
void gvCurvePaintable::paint(gvTexture*)
//===========================================================================
{
    if (!visible_) return;
    glEnableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_LIGHTING);
    if (selected_) {
	glColor4fv(selected_color_.rgba);
    } else {
	glColor4fv(normal_color_.rgba);
    }
    glVertexPointer(3, GL_DOUBLE, 0, line_.vertexArray());
    glDrawElements(GL_LINE_STRIP, line_.numVertices(),
		   GL_UNSIGNED_INT, line_.stripArray());
    if(draw_pts_) {
	// Inverting the colors for the points.
	if (!selected_) {
	    glColor4fv(selected_color_.rgba);
	} else {
	    glColor4fv(normal_color_.rgba);
	}
	double orig_sz;
	glGetDoublev(GL_POINT_SIZE, &orig_sz);
	glPointSize(5);
	glVertexPointer(3, GL_DOUBLE, 0, line_.vertexArray());
	glDrawArrays(GL_POINTS, 0, line_.numVertices());
	glPointSize((float)orig_sz);
    }
    glEnable(GL_LIGHTING);
    glDisableClientState(GL_VERTEX_ARRAY);
}
