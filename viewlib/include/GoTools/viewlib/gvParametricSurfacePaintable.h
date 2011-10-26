//===========================================================================
//                                                                           
// File: gvParametricSurfacePaintable.h                                                
//                                                                           
// Created: Thu Jun 21 16:46:34 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvParametricSurfacePaintable.h,v 1.1 2007-04-17 12:25:40 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVPARAMETRICSURFACEPAINTABLE_H
#define _GVPARAMETRICSURFACEPAINTABLE_H

#ifdef _MSC_VER
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif
#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#include "GoTools/geometry/BoundedSurface.h"

typedef Go::GenericTriMesh genMesh;

/** gvParametricSurfacePaintable: OpenGL calls for a parametric surface.
*/

class gvParametricSurfacePaintable : public gvPaintable
{
public:
    gvParametricSurfacePaintable(genMesh& tri,
				 const gvColor& ncolor,
				 const gvColor& scolor,
				 int id)
	: gvPaintable(ncolor, scolor, id),
	  tri_(tri)
    {}
    gvParametricSurfacePaintable(genMesh& tri,
				 const gvColor& ncolor,
				 int id)
	: gvPaintable(ncolor, id),
	  tri_(tri)
    {}
/*     // Testing whether call to GLU may generate a good triangulation. */
/*     gvParametricSurfacePaintable(std::shared_ptr<Go::BoundedSurface> surf, */
/* 				 const gvColor& ncolor, */
/* 				 int id) */
/* 	: gvPaintable(ncolor, id), */
/*       surf_(surf) */
/*     { */
/*       createSurface(); */
/*     } */
    virtual ~gvParametricSurfacePaintable();

    virtual void paint(gvTexture* texture);


protected:
    genMesh& tri_;

/*     std::shared_ptr<Go::BoundedSurface> surf_; */
/*     GLUnurbsObj* nurbSurface_; // Remove when done? */

    void createSurface();

    void drawSurface();

};


#endif // _GVPARAMETRICSURFACEPAINTABLE_H

