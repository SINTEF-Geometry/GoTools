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
  void replaceBoundary(shared_ptr<SplineSurface> surf,
		       shared_ptr<SplineCurve> curve,
		       int bd_idx, double tol);

  bool enforceCoefCoLinearity(shared_ptr<SplineSurface> sf1, int bd1, 
			      shared_ptr<SplineSurface> sf2, int bd2, 
			      double tol, 
			      std::vector<std::vector<int> >& enumeration);

} // of namespace ModifySurf

}; // end namespace Go

#endif // _MODIFYSURF_H
