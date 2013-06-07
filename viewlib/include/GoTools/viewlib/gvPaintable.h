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

