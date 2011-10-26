//===========================================================================
//                                                                           
// File: gvLineCloudPaintable.h                                              
//                                                                           
// Created: Tue Oct 29 09:41:34 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvLineCloudPaintable.h,v 1.1 2007-04-17 12:25:37 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVLINECLOUDPAINTABLE_H
#define _GVLINECLOUDPAINTABLE_H



#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/geometry/LineCloud.h"
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
#include <cmath>

/** Documentation ...
    etc
 */

class gvLineCloudPaintable : public gvPaintable
{
public:
    gvLineCloudPaintable(const Go::LineCloud& lc,
			 const gvColor& ncolor,
			 const gvColor& scolor,
			 int id)
	: gvPaintable(ncolor, scolor, id),
	  lc_(lc), fractionrendered_(1.0)
    {}
    gvLineCloudPaintable(const Go::LineCloud& lc,
			  const gvColor& ncolor,
			  int id)
	: gvPaintable(ncolor, id),
	  lc_(lc), fractionrendered_(1.0)
    {}
    virtual ~gvLineCloudPaintable()
    {}

    virtual void paint(gvTexture*)
    {
	if (!visible_) return;
	glEnableClientState(GL_VERTEX_ARRAY);
	glDisable(GL_LIGHTING);
	if (selected_) {
	    glColor4fv(selected_color_.rgba);
	} else {
	    glColor4fv(normal_color_.rgba);
	}
	glVertexPointer(3, GL_DOUBLE, 0, lc_.point(0).begin());
	int numrendered = int(floor(lc_.numLines() * fractionrendered_) + 0.5);
	glDrawArrays(GL_LINES, 0, numrendered*2);
	glEnable(GL_LIGHTING);
	glDisableClientState(GL_VERTEX_ARRAY);
    }

    void setFractionRendered(double f)
    {
	fractionrendered_ = f;
    }

    double fractionRendered()
    {
	return fractionrendered_;
    }

protected:
    const Go::LineCloud& lc_;
    double fractionrendered_;
};






#endif // _GVLINECLOUDPAINTABLE_H

