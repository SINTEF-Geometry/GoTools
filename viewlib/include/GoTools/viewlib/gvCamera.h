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

