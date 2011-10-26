//===========================================================================
//                                                                           
// File: gvGenericTriPaintable.h                                             
//                                                                           
// Created: Fri Nov 30 14:25:16 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvGenericTriPaintable.h,v 1.1 2007-04-17 12:25:35 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVGENERICTRIPAINTABLE_H
#define _GVGENERICTRIPAINTABLE_H


#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif
#include <GL/gl.h>


namespace
{
    inline GLenum GLfptype(float)
    {
	return GL_FLOAT;
    }
    inline GLenum GLfptype(double)
    {
	return GL_DOUBLE;
    }
}

/** Documentation ...
    etc
 */

template <typename FloatType>
class gvGenericTriPaintable : public gvPaintable
{
public:
    gvGenericTriPaintable(GenericTriMesh<FloatType>& tri,
			  const gvColor& ncolor,
			  const gvColor& scolor,
			  int id)
	: gvPaintable(ncolor, scolor, id),
	  tri_(tri)
    {}
    gvGenericTriPaintable(GenericTriMesh<FloatType>& tri,
			  const gvColor& ncolor,
			  int id)
	: gvPaintable(ncolor, id),
	  tri_(tri)
    {}
    virtual ~gvGenericTriPaintable()
    {}

    virtual void paint(gvTexture*)
    {
	GLenum fptype = GLfptype(FloatType());
	// Draw surfaces a little bit behind other things 
	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1.0, 1.0);

	// Deal with visibility and selection state
	if (!visible_) return;
	if (selected_) {
	    glMaterialfv(GL_FRONT_AND_BACK,
			 GL_AMBIENT_AND_DIFFUSE,
			 selected_color_.rgba);
	} else {
	    glMaterialfv(GL_FRONT_AND_BACK,
			 GL_AMBIENT_AND_DIFFUSE,
			 normal_color_.rgba);
	}
	GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);

    /// Set up the vertex arrays
	glEnableClientState(GL_VERTEX_ARRAY);
    /// Set up the normal arrays
	glEnableClientState(GL_NORMAL_ARRAY);
	glVertexPointer(3, fptype, 0, tri_.vertexArray());
	//      std::cout << tri_.vertexArray()[300] << ' '
	//  	      << tri_.vertexArray()[301] << ' '
	//  	      << tri_.vertexArray()[302] << std::endl;
	// @@@ We should check for the use of normals and textures in tri_.
	glNormalPointer(fptype, 0, tri_.normalArray());

    /// Draw the triangles
	glDrawElements(GL_TRIANGLES, tri_.numTriangles()*3,
		       GL_UNSIGNED_INT, tri_.triangleIndexArray());

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisable(GL_POLYGON_OFFSET_FILL);
    }


protected:
    GenericTriMesh<FloatType>& tri_;
};




#endif // _GVGENERICTRIPAINTABLE_H

