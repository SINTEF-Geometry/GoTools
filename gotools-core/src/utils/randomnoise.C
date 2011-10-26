//===========================================================================
//                                                                           
// File: randomnoise.C                                                       
//                                                                           
// Created: Thu Feb 12 09:55:27 2004                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: randomnoise.C,v 1.3 2008-12-05 11:37:54 jnygaard Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
