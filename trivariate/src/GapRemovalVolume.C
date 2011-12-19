//===========================================================================
//                                                                           
// File: GapRemovalVolume.C                                                  
//                                                                           
// Created: Wed Nov  2 16:12:41 2011                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/trivariate/GapRemovalVolume.h"


namespace Go
{


//===========================================================================
void
GapRemoval::removeGapSpline(shared_ptr<SplineVolume>& vol1, 
			    shared_ptr<SurfaceOnVolume>& bd_sf1,
			    double sf1_start1, double sf1_end1,
			    double sf1_start2, double sf1_end2,
			    shared_ptr<SplineVolume>& vol2, 
			    shared_ptr<SurfaceOnVolume>& bd_sf2,
			    double sf2_start1, double sf2_end1,
			    double sf2_start2, double sf2_end2,
			    Point vertex_ll, Point vertex_ur,
			    double epsge, int orientation)
//===========================================================================
{
  MESSAGE("removeGapSpline() not yet implemented.");

  // @@sbr201111 I guess a crv for domain to be altered may be more
  // appropriate.  Otherwise we expect the domain to correspond to
  // iso-curves.

}


}
