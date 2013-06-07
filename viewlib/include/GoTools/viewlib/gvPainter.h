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

#ifndef GV_PAINTER_H_INCLUDED
#define GV_PAINTER_H_INCLUDED

#include <vector>
#include "GoTools/utils/config.h"
#include "GoTools/viewlib/gvPaintable.h"
#include "GoTools/viewlib/gvCamera.h"
#include "GoTools/viewlib/gvTexture.h"
#include "GoTools/utils/errormacros.h"

class QGLWidget;

/// The gvPainter is responsible for painting the 3D scene.
class gvPainter
{
public:
    gvPainter()
    {}
    virtual ~gvPainter();

    /// Draws the 3D scene by calling the paint() virtual method of
    /// all its contained gvPaintables.
    virtual void drawScene(QGLWidget *w=NULL);

    /// Add a new paintable. gvPaintable is the abstract bas class
    /// of all paintables.
    void addPaintable(shared_ptr<gvPaintable> pa)
    {
	paintables_.push_back(pa);
	textures_.push_back(shared_ptr<gvTexture>());
    }
//     void addPaintable(shared_ptr<gvPaintable> pa, int id)
//     {
// 	ASSERT(id < paintables_.size() + 1);
// 	paintables_.insert(paintables_.begin() + id, pa);
// 	textures_.insert(textures_.begin() + id, shared_ptr<gvTexture>());
//     }
    /// Removes all paintables, if you want to paint something
    /// you have to addPaintable() again.
    void removeAllPaintables()
    {
	paintables_.clear();
	textures_.clear();
    }

    void removePaintable(int id)
    {
      paintables_[id].reset();
      textures_[id].reset();
    }

    void removeLastPaintable()
    {
	paintables_.erase(paintables_.end() - 1);
	textures_.erase(textures_.end() - 1);
    }
    /// Call this to access the individual paintables, in order to
    /// manipulate selection state, visibility state or color.
    gvPaintable& getPaintable(int index)
    {
	return *(paintables_[index]);
    }

    void setTexture(int index, shared_ptr<gvTexture> tex)
    {
	ASSERT(index < int(textures_.size()));
	textures_[index] = tex;
    }

    shared_ptr<gvTexture> getTexture(int index)
    {
      ASSERT(index < int(textures_.size()));
      return(textures_[index]);
    }

    shared_ptr<const gvTexture> setTexture(int index)
    {
	ASSERT(index < int(textures_.size()));
	return textures_[index];
    }

protected:
    std::vector< shared_ptr<gvPaintable> > paintables_;
    std::vector< shared_ptr<gvTexture> > textures_;
    
};


#endif // of #ifndef GV_PAINTER_H_INCLUDED
