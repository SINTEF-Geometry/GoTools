//===========================================================================
//                                                                           
// File: gvPainter.h                                                         
//                                                                           
// Created: Mon Apr 30 13:03:20 2001                                         
//                                                                           
// Author: Jens Olav Nygaard <jnygaard@sintef.math.no>
//                                                                           
// Revision: $Id: gvPainter.h,v 1.1 2007-04-17 12:25:40 sbr Exp $
//                                                                           
// Description: The gvPainter is responsible for painting the 3D scene.
//                                                                           
//===========================================================================


#ifndef GV_PAINTER_H_INCLUDED
#define GV_PAINTER_H_INCLUDED

#include <vector>
#include <memory>
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
    void addPaintable(std::shared_ptr<gvPaintable> pa)
    {
	paintables_.push_back(pa);
	textures_.push_back(std::shared_ptr<gvTexture>());
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

    void setTexture(int index, std::shared_ptr<gvTexture> tex)
    {
	ASSERT(index < int(textures_.size()));
	textures_[index] = tex;
    }

    std::shared_ptr<gvTexture> getTexture(int index)
    {
      ASSERT(index < int(textures_.size()));
      return(textures_[index]);
    }

    std::shared_ptr<const gvTexture> setTexture(int index)
    {
	ASSERT(index < int(textures_.size()));
	return textures_[index];
    }

protected:
    std::vector< std::shared_ptr<gvPaintable> > paintables_;
    std::vector< std::shared_ptr<gvTexture> > textures_;
    
};


#endif // of #ifndef GV_PAINTER_H_INCLUDED
