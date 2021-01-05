//===========================================================================
//                                                                           
// File: LRBenchmarkUtils3D.h                                                
//                                                                           
// Created: Thu Apr 17 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LRBENCHMARKUTILS3D_H
#define _LRBENCHMARKUTILS3D_H


#include "GoTools/lrsplines3D/LRSplineVolume.h"


namespace Go
{

    double benchmarkVolRefinement(LRSplineVolume& lr_vol,
				  const std::vector<LRSplineVolume::Refinement3D>& refs,
				  bool single_insertions = false);

}

#endif // _LRBENCHMARKUTILS3D_H

