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
  /// Utility functionality for splitting SurfaceModels
  namespace SplitModelUtils
  {
    void splitInFreeCorners(shared_ptr<SurfaceModel> sfmodel,
			    const Point& pnt, const Point& axis);

     void splitInNonCorners(shared_ptr<SurfaceModel> sfmodel,
			    const Point& pnt, const Point& axis);

     void splitInOuterVertices(shared_ptr<SurfaceModel> sfmodel,
			       shared_ptr<ftSurface> face,
			       const Point& pnt, const Point& axis);
  }
}

#endif    // #ifndef __SPLITMODELUTILS_H
