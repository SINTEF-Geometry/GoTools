//===========================================================================
//                                                                           
// File: GapRemovalVolume.h                                                  
//                                                                           
// Created: Wed Nov  2 16:12:32 2011                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GAPREMOVALVOLUME_H
#define _GAPREMOVALVOLUME_H


#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"



namespace Go {

  // Removal of gaps between two adjacent surfaces of various types

namespace GapRemoval
{


  // We average the volumes along the matching faces on the
  // rectangular domain given by vertex lower left and upper right
  // (assuming that such a domain is well defined).
  void
  removeGapSpline(shared_ptr<SplineVolume>& vol1, 
		  shared_ptr<SurfaceOnVolume>& bd_sf1,
		  double sf1_start1, double sf1_end1,
		  double sf1_start2, double sf1_end2,
		  shared_ptr<SplineVolume>& vol2, 
		  shared_ptr<SurfaceOnVolume>& bd_sf2,
		  double sf2_start1, double sf2_end1,
		  double sf2_start2, double sf2_end2,
		  Point vertex_ll, Point vertex_ur,
		  double epsge, int orientation);

}

}

#endif // _GAPREMOVALVOLUME_H

