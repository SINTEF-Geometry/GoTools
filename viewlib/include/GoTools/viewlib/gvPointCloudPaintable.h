//===========================================================================
//                                                                           
// File: gvPointCloudPaintable.h                                             
//                                                                           
// Created: Wed Mar 13 15:38:25 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvPointCloudPaintable.h,v 1.3 2008-03-07 13:03:52 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
      pc_(pc), fractionrendered_(1.0), pointsize_(5.0), paintId_(paintId)
    {}
    gvPointCloudPaintable(const Go::PointCloud3D& pc,
			  const gvColor& ncolor,
			  int id, bool paintId=false)
	: gvPaintable(ncolor, id),
      pc_(pc), fractionrendered_(1.0), pointsize_(5.0), paintId_(paintId)
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

protected:
    const Go::PointCloud3D& pc_;
    double fractionrendered_;
    double pointsize_;
    bool paintId_;
};



#endif // _GVPOINTCLOUDPAINTABLE_H






