//===========================================================================
//                                                                           
// File: SurfaceOnVolumeTools.C
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

#include "GoTools/trivariate/SurfaceOnVolumeTools.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/SurfaceTools.h"

namespace Go {

//===========================================================================
  CurveLoop 
  SurfaceOnVolumeTools::getOuterBoundaryLoop(shared_ptr<SurfaceOnVolume> sf,
					     double eps)
//===========================================================================
  {
    CurveLoop spacecrvs = SurfaceTools::outerBoundarySfLoop(sf->spaceSurface(), eps);

    int nmb = spacecrvs.size();
    for (int ki=0; ki<nmb; ++ki)
      {
	shared_ptr<ParamCurve> cv = spacecrvs[ki];
	shared_ptr<CurveOnSurface> sfcrv = 
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
	if (sfcrv.get())
	  {
	    sfcrv->setUnderlyingSurface(sf);
	    sfcrv->ensureParCrvExistence(eps);
	  }
      }
    return spacecrvs;
  }
}
