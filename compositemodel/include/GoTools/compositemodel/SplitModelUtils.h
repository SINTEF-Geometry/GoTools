//===========================================================================
//
// File : SplitModelUtils.h
//
// Created: Mars 2010
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description:
//
//===========================================================================


#ifndef __SPLITMODELUTILS_H
#define __SPLITMODELUTILS_H

#include "GoTools/compositemodel/SurfaceModel.h"

namespace Go
{
  /// Utility functionality for splitting of surface models. 
  namespace SplitModelUtils
  {
    /// The model is split from corners vertices belonging only to one face
    /// towards the axis defined by pnt and axis
    void splitInFreeCorners(shared_ptr<SurfaceModel> sfmodel,
			    const Point& pnt, const Point& axis);

    /// The model is split from vertices that do not constitute a corner for the
    /// current face towards the axis defined by pnt and axis
     void splitInNonCorners(shared_ptr<SurfaceModel> sfmodel,
			    const Point& pnt, const Point& axis);

     /// The model is split from vertices in outer loops towards the axis 
     /// defined by pnt and axis or towards vertices in inner loops
     void splitInOuterVertices(shared_ptr<SurfaceModel> sfmodel,
			       shared_ptr<ftSurface> face,
			       const Point& pnt, const Point& axis);
  }
}

#endif    // #ifndef __SPLITMODELUTILS_H
