//===========================================================================
//                                                                           
// File: gvRectangularVolumePaintable.h                                      
//                                                                           
// Created: Thu Jul  5 15:24:24 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVRECTANGULARVOLUMEPAINTABLE_H
#define _GVRECTANGULARVOLUMEPAINTABLE_H



#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/trivariate/RegularVolMesh.h"
//class RegularMesh;

/** gvRectangularSurfacePaintable: OpenGL calls for a parametric
    surface on a rectangular domain.
*/

class gvRectangularVolumePaintable : public gvPaintable
{
public:
    gvRectangularVolumePaintable(Go::RegularVolMesh& tri,
				 const gvColor& ncolor,
				 const gvColor& scolor,
				 int id)
	: gvPaintable(ncolor, scolor, id),
	  tri_(tri)
    {}
    gvRectangularVolumePaintable(Go::RegularVolMesh& tri,
				 const gvColor& ncolor,
				 int id)
	: gvPaintable(ncolor, id),
	  tri_(tri)
    {}
    virtual ~gvRectangularVolumePaintable();

    virtual void paint(gvTexture* texture);


protected:
    Go::RegularVolMesh& tri_;
};

#endif // _GVRECTANGULARVOLUMEPAINTABLE_H

