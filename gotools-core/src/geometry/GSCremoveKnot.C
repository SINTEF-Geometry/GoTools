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

#include "GoTools/geometry/SplineCurve.h"
#include <algorithm>
#include <iostream>

//#ifdef __BORLANDC__
#include <iterator> // For back_inserter.  This one should be required by VC++ and GCC as well...
//#endif

using std::back_inserter;

namespace Go{


//===========================================================================
void SplineCurve::removeKnot(double tpar)
//===========================================================================
{
    std::vector<double>::const_iterator ki = basis().begin();
    std::vector<double>::const_iterator kend = basis().end();
    std::vector<double>::const_iterator t_iter = std::find(ki, kend, tpar);
    ALWAYS_ERROR_IF(t_iter == kend,
		    "Attempting to remove unknown knot vector.");

    //ALWAYS_ERROR_IF(rational_,	"Not implemented for rational case!");

    int mt = 1; // left multiplicity of tpar
    // moving to last occurence of t_par
    while (t_iter[1] == tpar) {
	++t_iter;
	++mt;
    }
    int ti = (int)(t_iter - ki); // ti = index of t_iter
    int k = order();

    int n = k - mt + 2; // size of matrix

    // arrays for our computations
    std::vector<double> lambda;
    std::vector<double> mu;
    std::vector<double> sign;
    std::vector<double> coefs;
    std::vector<double> x;

    int dim = (rational_) ? dimension() + 1 : dimension();
    std::vector<double>::iterator sc = (rational_) ? rcoefs_begin() : 
      coefs_begin();

    // initialize arrays
    lambda.push_back(0.0);
    mu.push_back(1.0);
    sign.push_back((n & 1) ? 1.0 : (-1.0)); // (-1) if divisible by 2
    int i, j;
    for (j = 0; j < dim; ++j) {
	coefs.push_back(sc[(ti - k) * dim + j]);
	x.push_back(0.0);
    }
    for (i = 1; i < (n - 1); ++i) {
	lambda.push_back((ki[ti + i] - tpar) / (ki[ti + i] - ki[ti - k + i]));
	mu.push_back((tpar - ki[ti - k + i]) / (ki[ti + i] - ki[ti - k + i]));
	sign.push_back(sign[i - 1] * (-1.0));
	for (j = 0; j < dim; ++j) {
	    coefs.push_back(sc[(ti - k + i) * dim + j]);
	    x.push_back(0.0);
	}
    }
    lambda.push_back(1.0);
    mu.push_back(0.0);
    sign.push_back(1.0);
    for (j= 0; j < dim; ++j) {
	coefs.push_back(sc[(ti - mt +1) * dim + j]);
	x.push_back(0.0);
    }

    // Algorithm to solve our system; Gaussian elimination.

    //  1.
    int p = 0;
    while (mu[p + 1] >= 0.5) ++p;

    // 2a. Gaussian elimination.
    for (i = 1; i < p + 1; ++i) {
	sign[i] += - sign[i - 1] * lambda[i] / mu[i - 1];
	for (j = 0; j < dim; ++j)
	    coefs[i * dim + j] += 
		- coefs[(i - 1) * dim + j] * lambda[i] / mu[i - 1];
    }

    // 2b. Backward Gaussian elimination.
    for (i = n - 2; i > p; --i) {
	sign[i] += - sign[i + 1] * mu[i] / lambda[i + 1];
	for (j = 0; j < dim; ++j)
	    coefs[i * dim + j] += 
		- coefs[(i + 1) * dim + j] * mu[i] / lambda[i + 1];
    }

    //  3. Solve for x(p) and x(n-1).
    sign[p] += - sign[p + 1] * mu[p] / lambda[p + 1];
    for (j = 0; j < dim; ++j) {
	coefs[p * dim + j] += 
	    - coefs[(p + 1) * dim + j] * mu[p] / lambda[p + 1];
	x[(n - 1) * dim + j] = coefs[p * dim + j] / sign[p];
	x[p * dim + j] = (coefs[(p + 1) * dim + j] 
				  - sign[p + 1] * x[(n - 1) * dim + j]) 
	    / lambda[p + 1];
    }

    for (j = 0; j < dim; ++j) {
	// 4a. Backward determination of x(i).
	for (i = p - 1; i > -1; --i)
	    x[i * dim + j] = 
		(coefs[i * dim + j] 
		 - sign[i] * x[(n - 1) * dim + j]) / mu[i];
	
	// 4b. Forward determination of x(i-1).
	for (i = (p + 2); i < n; ++i)
	    x[(i - 1) * dim + j] = 
		(coefs[i * dim + j] 
		 - sign[i] * x[(n - 1) * dim + j]) / lambda[i];
    }

    // Algorithm finished.

    // Change our sc.
    for (j = 0; j < dim; ++j) {
	for (i = ti - k; i < ti - mt + 1; ++i)
 	    sc[i * dim + j] = x[(k - ti +i) * dim + j];
	for (i = ti - mt + 1; i < numCoefs() - 1; ++i)
	    sc[i * dim + j] = sc[(i + 1) * dim + j];
    }	

    if (rational_)
      rcoefs_.erase(rcoefs_.begin() + (numCoefs() - 1) * dim,
		    rcoefs_.begin() + numCoefs() * dim);
    else
      coefs_.erase(coefs_.begin() + (numCoefs() - 1) * dim,
		   coefs_.begin() + numCoefs() * dim);

    if (rational_)
      updateCoefsFromRcoefs();

    // Update basis_.
    std::vector<double> new_knots;
    std::copy(ki, ki + ti, std::back_inserter(new_knots));
    std::copy(ki + ti + 1, kend, std::back_inserter(new_knots));
    double *ks;
    ks = &new_knots[0];
    basis_ = BsplineBasis(numCoefs() - 1, order(), ks);


}


} // namespace Go;
