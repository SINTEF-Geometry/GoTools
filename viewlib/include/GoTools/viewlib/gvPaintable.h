//===========================================================================
//                                                                           
// File: gvPaintable.h                                                       
//                                                                           
// Created: Wed Jun 20 11:42:37 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvPaintable.h,v 1.1 2007-04-17 12:25:39 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVPAINTABLE_H
#define _GVPAINTABLE_H


#include "GoTools/viewlib/gvColor.h"
class gvTexture;

class QGLWidget;

/** gvPaintable: Super class for OpenGL calls to geometric objects.
*/

class gvPaintable
{
public:
    /// The arguments of this constructor (which will be called
    /// only by subclasses) are the colors to use in normal and
    /// selected states.
    gvPaintable(const gvColor& ncolor,
		const gvColor& scolor,
		int id)
	: visible_(true), selected_(false),
	  normal_color_(ncolor),
	  selected_color_(scolor),
	  id_(id)
    {}
    /// If the selected-mode color is not specifyed, a lighter
    /// shade of the normal color will be used.
    gvPaintable(const gvColor& ncolor, int id)
	: visible_(true), selected_(false),
	  normal_color_(ncolor),
	  id_(id)
    {
	// Default behaviour:
	// Set the selected_color_ to a lighter shade of the
	// normal_color_.
	selected_color_.rgba[0] = 0.5f + 0.5f*normal_color_.rgba[0];
	selected_color_.rgba[1] = 0.5f + 0.5f*normal_color_.rgba[1];
	selected_color_.rgba[2] = 0.5f + 0.5f*normal_color_.rgba[2];
	selected_color_.rgba[3] = 0.5f + 0.5f*normal_color_.rgba[3];
    }
    virtual ~gvPaintable();

    /// The function that does the painting. Must be overridden
    /// by all subclasses. If your subclass does not use textures,
    /// just ignore the texture argument.
    virtual void paint(gvTexture* texture) = 0;

    virtual void paintGL(QGLWidget *w, gvTexture* texture)
    {
       paint(texture);
    }

    void setVisible(bool state)
    {
	visible_ = state;
    }
    bool visible() const
    {
	return visible_;
    }
    void setSelected(bool state)
    {
	selected_ = state;
    }
    bool selected() const
    {
	return selected_;
    }
    void setColor(const gvColor& ncolor,
		  const gvColor& scolor)
    {
	normal_color_ = ncolor;
	selected_color_ = scolor;
    }
    void setColor(const gvColor& ncolor)
    {
	normal_color_ = ncolor;
	// Set the selected_color_ to a lighter shade of the
	// normal_color_.
	selected_color_.rgba[0] = 0.5f + 0.5f*normal_color_.rgba[0];
	selected_color_.rgba[1] = 0.5f + 0.5f*normal_color_.rgba[1];
	selected_color_.rgba[2] = 0.5f + 0.5f*normal_color_.rgba[2];
	selected_color_.rgba[3] = 0.5f + 0.5f*normal_color_.rgba[3];
    }
    gvColor getNormalColor()
    {
	return normal_color_;
    }
    int id() const
    {
	return id_;
    }
    void setId(int id)
    {
	id_ = id;
    }

protected:
    bool visible_;
    bool selected_;
    gvColor normal_color_;
    gvColor selected_color_;
    int id_;
};



#endif // _GVPAINTABLE_H

