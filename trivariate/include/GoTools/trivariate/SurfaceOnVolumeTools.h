//===========================================================================
//                                                                           
// File: SurfaceOnVolumeTools.h
//                                                                           
// Created: Jan. 2010
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SURFACEONVOLUMETOOLS_H
#define _SURFACEONVOLUMETOOLS_H

#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/CurveLoop.h"

namespace Go {

/// This namespace contains functions used to get boundary curves on volumes
  namespace SurfaceOnVolumeTools {
    
    CurveLoop getOuterBoundaryLoop(shared_ptr<SurfaceOnVolume> sf,
				   double eps);
  }
}
#endif // _SURFACEONVOLUMETOOLS_H
