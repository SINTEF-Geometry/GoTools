//===========================================================================
//                                                                           
// File: gvNoopPaintable.h                                                   
//                                                                           
// Created: Fri Apr  4 07:07:30 2003                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: gvNoopPaintable.h,v 1.1 2007-04-17 12:25:38 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#ifndef _GVNOOPPAINTABLE_H
#define _GVNOOPPAINTABLE_H


#include "GoTools/viewlib/gvColor.h"
#include "GoTools/viewlib/gvPaintable.h"
class gvTexture;

/** Documentation ...
 */

class gvNoopPaintable : public gvPaintable
{
public:
    gvNoopPaintable(const gvColor& ncolor,
		    const gvColor& scolor,
		    int id)
	: gvPaintable(ncolor, scolor, id)
    {}
    gvNoopPaintable(const gvColor& ncolor,
		    int id)
	: gvPaintable(ncolor, id)
    {}

    virtual ~gvNoopPaintable();
    virtual void paint(gvTexture* texture);
};


#endif // _GVNOOPPAINTABLE_H

