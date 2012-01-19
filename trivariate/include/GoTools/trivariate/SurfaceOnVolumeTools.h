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

/// This namespace contains functions used to fetch boundary curves from
/// surface on volumes
  namespace SurfaceOnVolumeTools {
    
    /// Fetch the outer loop of the surface. This function is implemented
    /// outside of SurfaceOnVolume to allow the resulting curves to be
    /// of type CurveOnSurface and possess all the information related to
    /// this curve type.
    CurveLoop getOuterBoundaryLoop(shared_ptr<SurfaceOnVolume> sf,
				   double eps);
  }
}
#endif // _SURFACEONVOLUMETOOLS_H
