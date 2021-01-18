//===========================================================================
//                                                                           
// File: LRBenchmarkUtils3D.C                                                
//                                                                           
// Created: Wed Apr 17 17:49:51 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================



#include "GoTools/lrsplines3D/LRBenchmarkUtils3D.h"
#include "GoTools/utils/timeutils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;


namespace Go
{

double benchmarkVolRefinement(LRSplineVolume& lr_vol,
			      const vector<LRSplineVolume::Refinement3D>& refs,
			      bool single_insertions)
{
  if (refs.size() == 0)
    return 0.0;

    // double time0 = omp_get_wtime();
    double time0 = getCurrentTime();

    if (single_insertions)
      { // On refinement at the time.
	for (size_t ki = 0; ki < refs.size(); ++ki)
	  {
#if 0//ndef NDEBUG
	    std::cout << "DEBUG: ref num: ki = " << ki << std::endl;
#endif
	    lr_vol.refine(refs[ki]);
	}
    }
    else
    {
	lr_vol.refine(refs);
    }

    // double time1 = omp_get_wtime(); // Measured in seconds.
    double time1 = getCurrentTime();
    double time_spent = time1 - time0;

    return time_spent;
}

}
