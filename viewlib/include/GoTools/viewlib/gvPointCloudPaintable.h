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

#ifndef _GVPOINTCLOUDPAINTABLE_H
#define _GVPOINTCLOUDPAINTABLE_H



#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/geometry/PointCloud.h"
#ifdef _MSC_VER
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define NOGDI
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <cmath>

#include <QString>
#include <QtOpenGL>

/** Documentation ...
    etc
 */

class gvPointCloudPaintable : public gvPaintable
{
public:
    gvPointCloudPaintable(const Go::PointCloud3D& pc,
			  const gvColor& ncolor,
			  const gvColor& scolor,
			  int id, bool paintId=false)
	: gvPaintable(ncolor, scolor, id),
	  pc_(pc), fractionrendered_(1.0), pointsize_(5.0), paintId_(paintId),
	  vert_translation_(3, 0.0)
    {}
    gvPointCloudPaintable(const Go::PointCloud3D& pc,
			  const gvColor& ncolor,
			  int id, bool paintId=false)
	: gvPaintable(ncolor, id),
	  pc_(pc), fractionrendered_(1.0), pointsize_(5.0), paintId_(paintId),
	  vert_translation_(3, 0.0)
    {}
    virtual ~gvPointCloudPaintable()
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
	glVertexPointer(3, GL_DOUBLE, 0, pc_.point(0).begin());
	int numrendered = int(floor(pc_.numPoints() * fractionrendered_)
			      + 0.5);

	double orig_sz;
	glGetDoublev(GL_POINT_SIZE, &orig_sz);
	glPointSize((GLfloat)pointsize_);
	glDrawArrays(GL_POINTS, 0, numrendered);
	glPointSize((GLfloat)orig_sz);
	glEnable(GL_LIGHTING);
	glDisableClientState(GL_VERTEX_ARRAY);
    }

    virtual void paintGL(QGLWidget *w, gvTexture* texture)
    {
       if (!visible_) return;
       gvPaintable::paintGL(w, texture);
       glDisable(GL_LIGHTING);
       if (paintId_ && w)
       {
	  for (int i=0; i<pc_.numPoints(); i++)
	  {
	     Go::Vector3D p=pc_.point(i);
	     if (selected_) {
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, 
			     selected_color_.rgba);
	     } else {
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, 
			     normal_color_.rgba);
	     }
#if QT_VERSION >= 0x030100
	     w->renderText(p[0], p[1], p[2], QString::number(id_));
#endif
	  }
       }
       glEnable(GL_LIGHTING);
    }

    void setFractionRendered(double f)
    {
	fractionrendered_ = f;
    }

    double fractionRendered()
    {
	return fractionrendered_;
    }

    void setPointSize(double sz)
    {
	pointsize_ = sz;
    }

    double pointSize()
    {
	return pointsize_;
    }

    bool getPaintId() const
    {
       return paintId_;
    }

    void setPaintId(bool paintId) 
    {
       paintId_=paintId;
    }

    /// Translate all vertices by vert_translation, wrt the geometry, discarding any previus translation.
    void translate(const std::vector<double>& vert_translation)
    {
	MESSAGE("Not implemented yet!");
    }

protected:
    const Go::PointCloud3D& pc_;
    double fractionrendered_;
    double pointsize_;
    bool paintId_;

    /// The 3D-translation wrt the geometry.
    std::vector<double> vert_translation_;

};



#endif // _GVPOINTCLOUDPAINTABLE_H






