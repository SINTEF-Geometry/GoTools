//===========================================================================
//                                                                           
// File: gvRectangularSurfacePaintable.h                                                
//                                                                           
// Created: Thu Jun 21 16:46:34 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvRectangularSurfacePaintable.h,v 1.1 2007-04-17 12:25:43 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVRECTANGULARSURFACEPAINTABLE_H
#define _GVRECTANGULARSURFACEPAINTABLE_H


#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/tesselator/RegularMesh.h"
//class RegularMesh;

/** gvRectangularSurfacePaintable: OpenGL calls for a parametric
    surface on a rectangular domain.
*/

class gvRectangularSurfacePaintable : public gvPaintable
{
public:
    gvRectangularSurfacePaintable(Go::RegularMesh& tri,
		       const gvColor& ncolor,
		       const gvColor& scolor,
		       int id)
	: gvPaintable(ncolor, scolor, id),
	  tri_(tri)
    {}
    gvRectangularSurfacePaintable(Go::RegularMesh& tri,
		       const gvColor& ncolor,
		       int id)
	: gvPaintable(ncolor, id),
	  tri_(tri)
    {}
    virtual ~gvRectangularSurfacePaintable();

    virtual void paint(gvTexture* texture);


protected:
    Go::RegularMesh& tri_;
};


#endif // _GVRECTANGULARSURFACEPAINTABLE_H

