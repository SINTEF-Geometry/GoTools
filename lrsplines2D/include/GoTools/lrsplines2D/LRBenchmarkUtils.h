//===========================================================================
//                                                                           
// File: LRBenchmarkUtils.h                                                  
//                                                                           
// Created: Thu Nov 22 15:41:59 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LRBENCHMARKUTILS_H
#define _LRBENCHMARKUTILS_H


#include "GoTools/lrsplines2D/LRSplineSurface.h"


namespace Go
{

    double benchmarkSfRefinement(LRSplineSurface& lr_sf,
				 const std::vector<LRSplineSurface::Refinement2D>& refs,
				 bool single_insertions = false);

}

#endif // _LRBENCHMARKUTILS_H

