//===========================================================================
//                                                                           
// File: LRBenchmarkUtils.C                                                  
//                                                                           
// Created: Thu Nov 22 15:48:28 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/lrsplines2D/LRBenchmarkUtils.h"
#include "GoTools/utils/timeutils.h"
#include <omp.h>


using std::vector;


namespace Go
{

double benchmarkSfRefinement(LRSplineSurface& lr_sf,
			     const vector<LRSplineSurface::Refinement2D>& refs,
			     bool single_insertions)
{
    // double time0 = omp_get_wtime();
    double time0 = getCurrentTime();

    if (single_insertions)
    { // On refinement at the time.
	for (size_t ki = 0; ki < refs.size(); ++ki)
	{
#ifndef NDEBUG
	    std::cout << "DEBUG: ref num: ki = " << ki << std::endl;
#endif
	    lr_sf.refine(refs[ki]);
	}
    }
    else
    {
	lr_sf.refine(refs);
    }

    // double time1 = omp_get_wtime(); // Measured in seconds.
    double time1 = getCurrentTime();
    double time_spent = time1 - time0;

    return time_spent;
}

}
