//===========================================================================
//                                                                           
// File: gvView.h                                                            
//                                                                           
// Created: Fri Apr 20 16:26:18 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvView.h,v 1.3 2007-05-02 14:39:25 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVVIEW_H
#define _GVVIEW_H


#include "GoTools/viewlib/gvCamera.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/viewlib/gvObserver.h"
#include "GoTools/viewlib/gvTexture.h"
#include "GoTools/geometry/ParamSurface.h"

// Qt includes
#include <QtOpenGL>
//Added by qt3to4:
#include <QtGui/QMouseEvent>
#include <QtGui/QKeyEvent>

// Standard includes
#include <iostream>
#include <vector>


class gvData;


/** Documentation ...
 *   etc
 */

class gvView : public QGLWidget, public gvObserver
{

Q_OBJECT

public:
    /// The gvData& argument is the gvData object this object
    /// should observe. The other arguments are standard arguments
    /// for the QGlWidget superclass.
    gvView(gvData& data,
	   QWidget* parent=0, const char* name=0,
	   const QGLWidget * shareWidget = 0, Qt::WFlags f=0);
    virtual ~gvView();

    /// This function does the actual drawing into the widget's area.
    /// After clearing the area, it will call use() on the camera,
    /// and ask the gvPainter owned by the observed gvData object to
    /// do a repaint.
    virtual void paintGL();

    virtual void saveSnapshot(int w, int h, const QString &filename);
    /// Initialize GL state (shademodel, lights, etc.)
    virtual void initializeGL();
    /// Called when the widget is resized.
    /// It will set the viewport of the camera and then
    /// do a repaint.
    virtual void resizeGL( int w, int h );

    /// Searching among visible and selected objects.
    virtual void getObjAndParam(int mousex, int mousey, 
				shared_ptr<const Go::ParamSurface> &obj,
				double &tex_u, double &tex_v);

    /// Pick the geometrical object at a mouse coordinate.
    void pick(int mousex, int mousey);
    /// Pick the geometrical objects in a region of the window.
    void pickRegion(int mousex, int mousey, int w, int h);

    // Given input of 3D point in world coordinates, transform to window coordinates.
    void getWindowCoords(const Vector3D& pt, int& mousex, int& mousey) const;

    /// @@ Not sure if I need this
    virtual QSize sizeHint() const;

    /// This function is called whenever the observed gvData
    /// object is changed.
    virtual void observedChanged();

    /// Put the center of the bounding box of the model in the
    /// middle of the screen, at a (hopefully) useful distance
    /// from the camera.
    void focusOnBox();
    void focusOnVisible();

    /// Are we in wireframe mode?
    bool wireframe() { return wireframe_; }
    /// Are we in selection mode?
    bool selectionmode() { return selection_mode_; }
    /// Are we in paint-axis mode?
    bool axis() { return axis_; }
    /// Are we in bacface-culling mode?
    bool backCull() { return backcull_; }
    /// Are we in specular highlight mode?
    bool specular() { return specular_; }
    /// Are we in perspective projection mode?
    bool perspective() { return camera_.getPerspectiveMode(); }
    /// Are we in feedback mode?
    bool feedbackmode() const { return feedback_mode_;}
    /// Should the alpha value be used for blending?
    bool blendingmode() const { return blending_mode_;}


  virtual QSize minimumSizeHint() const
  {
    return QSize(50, 50);
  }

public slots:
    /// Set wireframe state.
    void setWireframe(bool mode);
    /// Set selection mode.
    void setSelectionmode(bool mode);
    /// Set paint-axis state.
    void setAxis(bool mode);
    /// Set backface-culling state.
    void setBackCull(bool mode);
    /// Set specular highlighting state.
    void setSpecular(bool mode);
    /// Set perspective projection state.
    void setPerspective(bool mode);
    /// Set feedback state.
    void setFeedbackmode(bool mode);
    /// Set blending state
    void setBlendingmode(bool mode);
    /// Set get click mode
    void setGetClickmode(bool mode);
    /// Centers view on point at x,y in back buffer
    void setCenter(int x, int y);

signals:
    /// An object has been picked
    void objectPicked(unsigned int name);
    /// An object has been picked
    void objectsPicked(unsigned int* names, int numnames);
    /// In feedback mode, indicates mouse position relative
    /// to start of drag.
    void feedback(int xrel, int yrel);

protected:
    gvView(const QGLFormat &format, gvData& data,
	   QWidget* parent, const char* name,
	   const QGLWidget * shareWidget, Qt::WFlags f);
    /// Overridden in order to make mouse moves translate,
    /// zoom and rotate the camera.
    virtual void mousePressEvent(QMouseEvent* e);
    /// Overridden in order to make mouse moves translate,
    /// zoom and rotate the camera.
    virtual void mouseReleaseEvent(QMouseEvent* e);
    /// Overridden in order to make mouse moves translate,
    /// zoom and rotate the camera.
    virtual void mouseMoveEvent(QMouseEvent* e);

//     /// Overridden.
//     virtual void keyPressEvent(QKeyEvent* e);
//     /// Overridden.
//     virtual void keyReleaseEvent(QKeyEvent* e);

    virtual bool get3Dpoint(int mousex, int mousey, Vector3D &objpt);

    void drawOverlay();

    gvTexture* makeFineCheckImage();
    gvTexture* makeCoarseCheckImage();
    static void texCoordInt(int num_objs, int ind, 
			    double &min_s, double &max_s,
			    double &min_t, double &max_t);


protected:
    bool gl_initialized_;
    bool no_data_;
    bool selection_mode_;
    bool selecting_;
    bool feedback_mode_;
    bool get_click_mode_;
    /// True if mouse was pressed in widget area and is not yet released.
    bool mouse_is_active_; 
//     bool keyboard_is_active_; 
    bool blending_mode_;

    gvData& data_;

    Go::BoundingBox box_;
    gvCamera camera_;
    gvCamera lights_camera_;
    double base_axis_size_;

//     ButtonState mouse_button_;
    Qt::MouseButton mouse_button_;
//     Qt::KeyboardModifier keyboard_modifier_;
//     QKeyEvent* keyboard_event_;
    QPoint last_mouse_pos_;
    QPoint starting_mouse_pos_;
    double draglength_;
    double unitx_;
    double unity_;

    bool wireframe_;
    bool axis_;
    bool backcull_;
    bool specular_;
    gvTexture *coarseTex_;
    gvTexture *fineTex_;

  QPainter* painter_;
};



#endif // _GVVIEW_H

