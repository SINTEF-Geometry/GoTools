//===========================================================================
//                                                                           
// File: gvCamera.h                                                          
//                                                                           
// Created: Wed Jun 20 12:42:23 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvCamera.h,v 1.1 2007-04-17 12:25:34 sbr Exp $
//                                                                           
// Description: A gvCamera deals with perspective and modelview
// transformations, with GL viewports, and with drawing 3D axes.
//                                                                           
//===========================================================================

#ifndef _GVCAMERA_H
#define _GVCAMERA_H


#include <utility>
#include "GoTools/utils/Array.h"
using Go::Vector3D;

/** gvCamera deals with perspective and modelview
 * transformations, with GL viewports, and with drawing 3D axes.
 */

class gvCamera
{
public:
    /// This constructor takes the initial focus point and the
    /// viewport as parameters. Default distance (from camera
    /// to focal point) is 10.0. Default axissize is zero.
    gvCamera(const Vector3D& fpoint, int x, int y, int w, int h);
    ~gvCamera();

    /// This is called to initialize the camera. You need a valid
    /// OpenGL context to call this.
    void initializeGL();

    /// Sets up the projection, viewport and modelview. Draws axis
    /// if axissize is nonzero.
    void use() const;

    /// Sets up the modelview only.
    void useModelView() const;

    /// Sets up the projection, viewport and modelview, for picking.
    void pick(int x, int y, int w, int h) const;

    /// Sets the focal point.
    /// The focal point is the point towards which the camera is
    /// looking. Also, the rotate() call rotates the camera
    /// about this point.
    void setFocalPoint(const Vector3D& fpoint)
    { focal_point_ = fpoint; }
    /// Gets the current focal point.
    void getFocalPoint(Vector3D& fpoint) const
    { fpoint = focal_point_; }

    /// Sets the distance.
    /// The distance is the modelspace distance from the camera
    /// to the focal point.
    void setDistance(double d)
    { distance_ = d; }
    /// Gets the current distance.
    void getDistance(double& d) const
    { d = distance_; }

    /// Sets the size (in modelspace units) of the axis.
    /// Setting the size to zero turns off the axis.
    void setAxisSize(double size)
    { axissize_ = size; }
    /// Gets the current axissize.
    void getAxisSize(double& size) const
    { size = axissize_; }

    Vector3D getWorldRay() const; 
    Vector3D getWorldUp() const; 
    Vector3D getWorldSide() const; 

    /// The axis of rotation is in screen coordinates, and will
    /// be translated to the correct coordinate system by the
    /// function. To give an example, if you detect a mouse drag
    /// in the direction (mx, my) you may want to rotate
    /// like this: camera.rotate(-my, mx, 0.1*sqrt(mx*mx + my*my));
    void rotate(double axisx, double axisy, double angle);
    /// Rotates the camera about the current (screen coordinate)
    /// z-axis.
    void rotateTransversal(double angle);
    /// Sets the internal modelview matrix to identity.
    /// The effect will be visible the next time use() is called.
    void resetRotation();

    /// Moves the focal point relative to screen coordinates
    void moveFocalPointRelative(const Vector3D& eyeptrel);
    /// Computes a ray in the model space corresponding to
    /// the screen point indicated.
    Vector3D eyeCoordsToObjectCoords(const Vector3D& eyept) const;

    Vector3D objectCoordsToWindowCoords(const Vector3D& pt) const;

    /// Sets the viewport.
    /// The viewport should be set whenever the window size changes.
    void setViewPort(int x, int y, int w, int h)
    { viewx_ = x; viewy_ = y; viewwidth_ = w; viewheight_ = h; }
    /// Gets the viewport.
    void getViewPort(int& x, int& y, int& w, int& h) const
    { x = viewx_; y = viewy_; w = viewwidth_; h = viewheight_; }

    /// Sets perspective mode. True means perspective, false
    /// means orthographic projection mode.
    void setPerspectiveMode(bool perspective_mode)
    { perspective_mode_ = perspective_mode; }
    /// Gets the current perspective mode.
    bool getPerspectiveMode() const
    { return perspective_mode_; }

    

protected:
    Vector3D focal_point_;
    double distance_;

    /// mv_ is the rotation matrix
    double mv_[16];

    double axissize_;

    int viewx_;
    int viewy_;
    int viewwidth_;
    int viewheight_;

    double viewfield_;

    double winz_;

    bool perspective_mode_;

    void computeWinz();
};


#endif // _GVCAMERA_H

