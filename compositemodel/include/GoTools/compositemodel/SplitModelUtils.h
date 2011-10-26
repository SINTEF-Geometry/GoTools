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
    void splitInFreeCorners(std::shared_ptr<SurfaceModel> sfmodel,
			    const Point& pnt, const Point& axis);

     void splitInNonCorners(std::shared_ptr<SurfaceModel> sfmodel,
			    const Point& pnt, const Point& axis);

     void splitInOuterVertices(std::shared_ptr<SurfaceModel> sfmodel,
			       std::shared_ptr<ftSurface> face,
			       const Point& pnt, const Point& axis);
  }
}

#endif    // #ifndef __SPLITMODELUTILS_H
