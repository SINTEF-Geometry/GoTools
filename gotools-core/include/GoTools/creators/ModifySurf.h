//===========================================================================
//                                                                           
// File: ModifySurf.h
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

#ifndef _MODIFYSURF_H
#define _MODIFYSURF_H

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

namespace Go {

namespace ModifySurf
{
  // Replace one boundary curve of a SplineSurface
  // Approximate the initial surface
  // bd_idx = 0: umin
  //          1: umax
  //          2: vmin
  //          3: vmax
  void replaceBoundary(std::shared_ptr<SplineSurface> surf,
		       std::shared_ptr<SplineCurve> curve,
		       int bd_idx, double tol);

} // of namespace ModifySurf

}; // end namespace Go

#endif // _MODIFYSURF_H
