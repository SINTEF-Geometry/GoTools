//===========================================================================
//                                                                           
// File: gvCurvePaintable.h                                                  
//                                                                           
// Created: Thu Jun 21 16:26:34 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvCurvePaintable.h,v 1.1 2007-04-17 12:25:34 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVCURVEPAINTABLE_H
#define _GVCURVEPAINTABLE_H


#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/tesselator/LineStrip.h"

//class LineStrip;

/** OpenGL calls for a parametric curve (represented by a
    LineStrip).
*/

class gvCurvePaintable : public gvPaintable
{
public:
    gvCurvePaintable(Go::LineStrip& line,
		     const gvColor& ncolor,
		     const gvColor& scolor,
		     int id)
	: gvPaintable(ncolor, scolor, id),
	  line_(line), draw_pts_(false)
    {}
    gvCurvePaintable(Go::LineStrip& line,
		     const gvColor& ncolor,
		     int id)
	: gvPaintable(ncolor, id),
	  line_(line), draw_pts_(false)
    {}
    virtual ~gvCurvePaintable();

    virtual void paint(gvTexture* texture);

    void setDrawPoints(bool dp)
    {
	draw_pts_ = dp;
    }
    bool drawPoints()
    {
	return draw_pts_;
    }


protected:
    Go::LineStrip& line_;
    bool draw_pts_;
};



#endif // _GVCURVEPAINTABLE_H

