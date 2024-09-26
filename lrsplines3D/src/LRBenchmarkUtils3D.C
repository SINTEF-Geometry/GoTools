/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
