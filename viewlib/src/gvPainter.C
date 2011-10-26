//===========================================================================
//                                                                           
// File: gvPainter.C                                                         
//                                                                           
// Created: Wed Jun 20 11:37:10 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvPainter.C,v 1.1 2007-04-17 12:25:54 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/viewlib/gvPainter.h"
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
gvPainter::~gvPainter()
//===========================================================================
{
}

//===========================================================================
void gvPainter::drawScene(QGLWidget *w)
//===========================================================================
{
    ASSERT(paintables_.size() == textures_.size());
    for (size_t i = 0; i < paintables_.size(); ++i) {
	if (paintables_[i].get() == 0)
	    continue;
	glLoadName(paintables_[i]->id()); // A no-op when not in selection mode
	if (w!=NULL)
	  paintables_[i]->paintGL(w, textures_[i].get());
	else
	paintables_[i]->paint(textures_[i].get());
    }
}

