//===========================================================================
//                                                                           
// File: gvQuadsPaintable.h                                                  
//                                                                           
// Created: Fri Jan 21 14:06:46 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvQuadsPaintable.h,v 1.1 2007-04-17 12:25:42 sbr Exp $
//                                                                           
//===========================================================================

#ifndef _GVQUADSPAINTABLE_H
#define _GVQUADSPAINTABLE_H

    /** Brief description. 
     *  Detailed description.
     */



#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/tesselator/QuadMesh.h"

class gvQuadsPaintable : public gvPaintable
{
public:
 gvQuadsPaintable(Go::QuadMesh& quads,
		     const gvColor& ncolor,
		     const gvColor& scolor,
		     int id)
	: gvPaintable(ncolor, scolor, id),
	  quads_(quads)
    {}
 gvQuadsPaintable(Go::QuadMesh& quads,
		     const gvColor& ncolor,
		     int id)
	: gvPaintable(ncolor, id),
	  quads_(quads)
    {}
    virtual ~gvQuadsPaintable();

    virtual void paint(gvTexture* texture);


protected:
    Go::QuadMesh& quads_;
};



#endif // _GVQUADSPAINTABLE_H

