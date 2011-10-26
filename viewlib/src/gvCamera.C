//===========================================================================
//                                                                           
// File: gvCamera.C                                                          
//                                                                           
// Created: Wed Jun 20 14:20:53 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvCamera.C,v 1.1 2007-04-17 12:25:50 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#include "GoTools/viewlib/gvCamera.h"
#include "GoTools/viewlib/gvUtilities.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/utils/Values.h"
#include <math.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <algorithm>

//===========================================================================
gvCamera::gvCamera(const Vector3D& fpoint, int x, int y, int w, int h)
//===========================================================================
    : focal_point_(fpoint),
      distance_(10.0),
      axissize_(0.0),
      viewx_(x), viewy_(y), viewwidth_(w), viewheight_(h),
      viewfield_(20.0),
      winz_(0.0),
      perspective_mode_(true)
{
}

//===========================================================================
gvCamera::~gvCamera()
//===========================================================================
{
}

//===========================================================================
void gvCamera::initializeGL()
//===========================================================================
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glGetDoublev(GL_MODELVIEW_MATRIX, mv_);
}

//===========================================================================
Vector3D gvCamera::getWorldRay()  const
//===========================================================================
{
  Vector3D wr(mv_[2], 
		mv_[2+4],
		mv_[2+8]);
  wr.normalize();
  return(wr);
}

//===========================================================================
Vector3D gvCamera::getWorldUp() const
//===========================================================================
{
  Vector3D up(mv_[1], 
		mv_[1+4],
		mv_[1+8]);
  up.normalize();
  return(up);
}


//===========================================================================
Vector3D gvCamera::getWorldSide() const
//===========================================================================
{
  Vector3D side(mv_[0], 
		  mv_[0+4],
		  mv_[0+8]);
  side.normalize();
  return(side);
}

//===========================================================================
void gvCamera::rotate(double axisx, double axisy, double angle)
//===========================================================================
{
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glRotated(angle, axisx, axisy, 0.0);
    glMultMatrixd(mv_);
    glGetDoublev(GL_MODELVIEW_MATRIX, mv_);    
    glPopMatrix();
}


//===========================================================================
void gvCamera::rotateTransversal(double angle)
//===========================================================================
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRotated(angle, 0.0, 0.0, 1.0);
    glMultMatrixd(mv_);
    glGetDoublev(GL_MODELVIEW_MATRIX, mv_);    
    glPopMatrix();
//      std::cout << angle << std::endl;
//      for (int i = 0; i < 16; ++i) {
//  	std::cout << mv_[i] << ((i+1)%4==0 ? '\n' : ' ');
//      }
}


//===========================================================================
void gvCamera::resetRotation()
//===========================================================================
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glGetDoublev(GL_MODELVIEW_MATRIX, mv_);    
    glPopMatrix();
}


//===========================================================================
void gvCamera::moveFocalPointRelative(const Vector3D& eyeptrel)
//===========================================================================
{
    // Compute the window z coordinate of the focal point
    int viewport[4];
    double projection[16];
    double modelview[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    Vector3D ppp;
    gluProject(focal_point_[0], focal_point_[1], focal_point_[2],
	       modelview, projection, viewport,
	       &ppp[0], &ppp[1], &ppp[2]);
    double winz = ppp[2];
    //    std::cout << "PPP = " << ppp << std::endl;

    // Interpret the third value of eyeptrel as a relative
    // amount to move the viewpoint
//     double objz = eyeptrel[2];

    Vector3D pmiddle(viewx_ + viewwidth_/2,
		       viewy_ + viewheight_/2,
		       winz);//(distance_-near_)/(far_-near_));
    Vector3D pnew = pmiddle + eyeptrel;
    Vector3D trpm = eyeCoordsToObjectCoords(pmiddle);
    Vector3D trpn = eyeCoordsToObjectCoords(pnew);
//      std::cout << "------------------\n" << pmiddle << pnew << trpm << trpn
//  	      << trpn - trpm;
    focal_point_ -= trpn;
    focal_point_ += trpm;
}


//===========================================================================
Vector3D gvCamera::eyeCoordsToObjectCoords(const Vector3D& eyept) const
//===========================================================================
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    use();

    int viewport[4];
    double projection[16];
    double modelview[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    Vector3D pt;
    int success = gluUnProject(eyept[0], eyept[1], eyept[2],
			       modelview, projection, viewport,
			       &pt[0], &pt[1], &pt[2]);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    ALWAYS_ERROR_IF(!success, "gluUnProject() failed!");
    return pt;
}

//===========================================================================
Vector3D gvCamera::objectCoordsToWindowCoords(const Vector3D& pt) const
//===========================================================================
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    use();

    int viewport[4];
    double projection[16];
    double modelview[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    Vector3D winpt;
    int success = gluProject(pt[0], pt[1], pt[2],
			       modelview, projection, viewport,
			       &winpt[0], &winpt[1], &winpt[2]);

    // It seems that the y value is referring to the opposite edge. @@sbr Fix this!
    winpt[1] = viewheight_ - winpt[1];

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    ALWAYS_ERROR_IF(!success, "gluProject() failed!");
    return winpt;
}


//===========================================================================
void gvCamera::use() const
//===========================================================================
{
    // Set up the OpenGL projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (perspective_mode_) {
	gluPerspective( /* field of view in degree */ viewfield_,
			/* aspect ratio*/ double(viewwidth_)/double(viewheight_),
			/* Z near */ distance_/10.0, /* Z far */ distance_*10);
    } else {
	double right = distance_*10 * tan(viewfield_*2*M_PI/360.0);
	double h = 2*right*double(viewheight_)/double(viewwidth_);
	glOrtho(-right, right, -h/2.0, h/2.0, distance_/10.0, distance_*10);
    }


    glViewport(viewx_, viewy_, viewwidth_, viewheight_);


    // Set up the modelview for the model
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(0.0, 0.0, -distance_);
    glMultMatrixd(mv_);
    if (axissize_ > 0.0)
	draw_gl_axes(axissize_);
    glTranslated(-focal_point_[0], -focal_point_[1], -focal_point_[2]);


//      for (int i = 0; i < 16; ++i) {
//  	std::cout << mv_[i] << ((i+1)%4==0 ? '\n' : ' ');
//      }

}

//===========================================================================
void gvCamera::useModelView() const
//===========================================================================
{
    // Set up the modelview for the model
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(0.0, 0.0, -distance_);
    glMultMatrixd(mv_);
    glTranslated(-focal_point_[0], -focal_point_[1], -focal_point_[2]);
}

//===========================================================================
void gvCamera::pick(int x, int y, int w, int h) const
//===========================================================================
{
    // Set up the OpenGL projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    int vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    gluPickMatrix(x, vp[3]-y-1, w, h, vp);
    gluPerspective( /* field of view in degree */ viewfield_,
		    /* aspect ratio*/  double(viewwidth_)/double(viewheight_),
		    /* Z near */ distance_/10.0, /* Z far */ distance_*10);

    //    glViewport(viewx_, viewy_, viewwidth_, viewheight_);


    // Set up the modelview for the model
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(0.0, 0.0, -distance_);
    glMultMatrixd(mv_);
    glTranslated(-focal_point_[0], -focal_point_[1], -focal_point_[2]);
}
