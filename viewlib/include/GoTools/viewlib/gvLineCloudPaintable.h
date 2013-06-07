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

