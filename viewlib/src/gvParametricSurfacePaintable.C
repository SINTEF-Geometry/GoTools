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

#include "GoTools/viewlib/gvParametricSurfacePaintable.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/viewlib/gvTexture.h"

#ifdef _MSC_VER
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif


using namespace Go;


//===========================================================================
gvParametricSurfacePaintable::~gvParametricSurfacePaintable()
//===========================================================================
{
//     gluDeleteNurbsRenderer(nurbSurface_);
}

//===========================================================================
void gvParametricSurfacePaintable::paint(gvTexture* texture)
//===========================================================================
{
    // Deal with visibility and selection state
    if (!visible_)
       return;

    // Draw surfaces a little bit behind other things 
    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1.0, 1.0);

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

    // Set up the vertex and normal arrays
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, tri_.vertexArray());
    //      std::cout << tri_.vertexArray()[300] << ' '
    //  	      << tri_.vertexArray()[301] << ' '
    //  	      << tri_.vertexArray()[302] << std::endl;
    // @@@ We should check for the use of normals and textures in tri_.
    glNormalPointer(GL_DOUBLE, 0, tri_.normalArray());

    bool use_textures = (texture != 0) && tri_.useTexCoords();
    if (use_textures) {
	glEnable(GL_TEXTURE_2D);
	texture->bind();
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glTexCoordPointer(2, GL_DOUBLE, 0, tri_.texcoordArray());
    }


    // Draw the triangle strips
    glDrawElements(GL_TRIANGLES, tri_.numTriangles()*3,
		   GL_UNSIGNED_INT, tri_.triangleIndexArray());

    if (use_textures) {
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisable(GL_TEXTURE_2D);
    }

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_POLYGON_OFFSET_FILL);
}




#if 0
//===========================================================================
void gvParametricSurfacePaintable::createSurface()
//===========================================================================
{
// setup nurb 
   nurbSurface_ = gluNewNurbsRenderer();
   gluNurbsProperty( nurbSurface_, GLU_SAMPLING_TOLERANCE, 25.0 );
   gluNurbsProperty( nurbSurface_, GLU_DISPLAY_MODE, GLU_FILL );
}

//=========================================================================== 
void gvParametricSurfacePaintable::drawSurface()
//===========================================================================
{
    //shared_ptr<BoundedSurface> trimmed_sf(surf_->clone());
    shared_ptr<BoundedSurface> trimmed_sf(surf_);
    trimmed_sf->setParameterDomain(0.0, 1.0, 0.0, 1.0);

    double loop_tol = 0.001;
    vector<CurveLoop> loops = trimmed_sf->allBoundaryLoops(0.001);
    shared_ptr<SplineSurface> under_sf(dynamic_pointer_cast<SplineSurface, ParamSurface>
					 (trimmed_sf->underlyingSurface()));
    ALWAYS_ERROR_IF(under_sf.get() == 0,
		"Member surface does not exist!", InputError());
    int dim = under_sf->dimension();
    double umin = under_sf->startparam_u();
    double umax = under_sf->endparam_u();
    double vmin = under_sf->startparam_v();
    double vmax = under_sf->endparam_v();

    vector<GLfloat> knots_u(under_sf->basis_u().begin(), under_sf->basis_u().end());
    vector<GLfloat> knots_v(under_sf->basis_v().begin(), under_sf->basis_v().end());
    vector<GLfloat> surf_coefs(under_sf->coefs_begin(), under_sf->coefs_end());

    glColor3f( 1.0, 1.0, 1.0 );
    gluBeginSurface( nurbSurface_ );
     gluNurbsSurface( nurbSurface_, knots_u.size(), &knots_u[0], knots_v.size(), &knots_v[0],
 		     dim, under_sf->numCoefs_u()*dim, &surf_coefs[0],
 		     under_sf->order_u(), under_sf->order_v(), GL_MAP2_VERTEX_3 );
    for (int ki = 0; ki < loops.size(); ++ki) {
	gluBeginTrim( nurbSurface_ );
	for (int kj = 0; kj < loops[ki].size(); ++kj) {
	    shared_ptr<CurveOnSurface> cv =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(loops[ki][kj]);
	    ALWAYS_ERROR_IF(cv.get() == 0,
			"Unexpected curve type.", InputError());
	    shared_ptr<SplineCurve> trim_cv =
		dynamic_pointer_cast<SplineCurve, ParamCurve>(cv->parameterCurve());
	    ALWAYS_ERROR_IF(trim_cv.get() == 0,
			"Unexpected curve type.", InputError());
	    vector<GLfloat> knots(trim_cv->basis().begin(), trim_cv->basis().end());
	    vector<GLfloat> cv_coefs(trim_cv->coefs_begin(), trim_cv->coefs_end());
	    gluNurbsCurve( nurbSurface_, knots.size(), &knots[0], 2,
			   &cv_coefs[0], trim_cv->order(), GLU_MAP1_TRIM_2);
	}
	gluEndTrim( nurbSurface_ );
    }
    gluEndSurface( nurbSurface_ );
}
#endif
