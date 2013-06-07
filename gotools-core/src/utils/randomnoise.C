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

#include <stdexcept>
#include <math.h>
#include <cstdlib>
#include "GoTools/utils/randomnoise.h"

using namespace std;

namespace Go {

//===========================================================================
void normalNoise(double* res, double mean_err, int num_samples)
//===========================================================================
{
    double scale_factor = double(2) / double(RAND_MAX);
    double v1, v2, r2, fact;
    for (int i = 0; i < num_samples; i+=2) {
	do {
	    v1 = double(rand()) * scale_factor - 1;
	    v2 = double(rand()) * scale_factor - 1;
	    r2 = v1 * v1 + v2 * v2;
	} while (r2 > 1 || r2 == double(0) || false);
	// we now know that the point (v1, v2) is within the unit circle, away from 0
	fact = sqrt(- 2 * log(r2)/r2);
	fact *= mean_err;
	res[i] = v1 * fact;
	if (i+1 < num_samples) {
	    res[i+1] = v2 * fact;
	}
    }
}

//===========================================================================
void uniformNoise(double* res, double lval, double uval, int num_samples)
//===========================================================================
{
    if (uval <= lval) {
	throw runtime_error("uniformNoise(...) : erroneous range.");
    }
    double range = uval - lval;
    double scale_factor = range / double(RAND_MAX);
    for (int i = 0; i < num_samples; ++i) {
	res[i] = double(rand()) * scale_factor + lval;
    }
}
 
}; // namespace Go
