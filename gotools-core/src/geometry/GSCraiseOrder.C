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

namespace Go{

//===========================================================================
void SplineCurve::raiseOrder(int raise)
//===========================================================================

/*********************************************************************
*
* PURPOSE    : 	To raise the description of a B-spline curve one order.
*
*
* INPUT &    : 	basis_.begin()	- Description of knot vector of original description
* VARIABLES	coefs_	- Coefficients of original description
*		numCoefs()	- Number of vertices of original description
*		order()	- Order of original description
*		dim_	- The dimension of the space in which the curve lies
*		new_knots	- Knot vector of the raised basis
*		num_newCoefs	- Number of vertices in the raised curve
*		ecc	- Array for internal use only
*		ecw	-        ---- " ----
*               new_coefs	- Coefs of the raised curve
*
*
*
*
* METHOD     : 	The order raising algorithm of Cohen, Lyche and Schumaker
*		is used.
*
*
* REFERENCES :	Fortran version:
*		T.Dokken, SI, 1984-06
*
*
* CALLS      :  s6err.
*
*
* WRITTEN BY : 	Christophe R. Birkeland, SI, 1991-07
* REWRITTEN BY : Sverre Briseid, SINTEF, 2001-05 (s1753() ported to C++)
* REVISED BY :
*
*********************************************************************
*/
{
    ALWAYS_ERROR_IF(raise < 0, "Raise must be positive!");

    bool rat = rational_;
    int kdim = (rat) ? dim_ + 1 : dim_;
    int j;
    // This may not be the most efficient method when raise is big.
    for (j = 0; j < raise; ++j) {

	std::vector<double>::iterator coef_iter =
	    (rat) ? rcoefs_.begin() : coefs_.begin();
	// lots of dummy variables
	int ki, kj, kk, kl, kr, kstop;/* Loop control variables 		*/
	int kjmid, ikmid;		/* kjmid=(kj-1)*kdim  ikmid=(order()-1)*kdim */
	double ty1, ty2, tyi, tyik;	/* Parameters used in Main Loop		*/
	double dummy;
	double tden;
	std::vector<double> new_knots;  // vector to hold knots of raised basis
	std::vector<double> new_coefs; // vector to hold coefs acc. to raised basis
	int new_order = order() + 1;

	// Initialize new_knots.
	for (int i = 0; i < numCoefs() + order(); ++i) {
	    new_knots.push_back(basis_.begin()[i]);
	    if ((i < numCoefs() + order() - 1) &&
		(basis_.begin()[i] < basis_.begin()[i + 1]))
		new_knots.push_back(basis_.begin()[i]);
	}
	new_knots.push_back(endparam());

	int num_newCoefs = (int)new_knots.size() - new_order;
	double* ecc = new double[num_newCoefs * kdim];
	double* ecw = new double[num_newCoefs * kdim];

	/* Initiate local variables. */
	kr = 1;

	for (kj = 1; kj <= num_newCoefs; kj++)
	    {
		/* Find kr, such that 
		   basis_.begin()[kr-1]<=new_knots[kj-1]<basis_.begin()[kr] */
		for (kr--; basis_.begin()[kr] <= new_knots[kj - 1]; kr++) ;


		/* Set ecc and ecw to zero. */
		for (ki = 0; ki < order() * kdim; ki++)
		    {
			ecc[ki] = (double) 0.0;
			ecw[ki] = (double) 0.0;
		    }

		/* Initialize the remaining ecc and ecw entries. */
		kstop = std::min(order(), numCoefs() + order() - kr);
		for (ki = std::max(0, order() - kr); ki < kstop; ki++)
		    for (kl = 0; kl < kdim; kl++)
			{
			    dummy = coef_iter[(ki + kr - order()) * kdim + kl];
			    ecc[ki * kdim + kl] = dummy;
			    ecw[ki * kdim + kl] = dummy;
			}

		/* MAIN LOOP. */
		for (kk = order() - 1; kk > 0; kk--)
		    {
			ty1 = new_knots[kj + kk - 1];
			ty2 = new_knots[kj + kk];
			kstop = std::max(order() - kk, order() - kr);

			for (ki = std::min(order() - 1, numCoefs() +2*order() - kk - kr - 1);
			     ki >= kstop; ki--)
			    {
				tyi = basis_.begin()[kr + ki - order()];
				tyik = basis_.begin()[kr + ki + kk - order()];
				tden = tyik - tyi;

				for (kl = 0; kl < kdim; kl++)
				    {
					ecc[ki * kdim + kl] =
					    ((ty2 - tyi) * ecc[ki * kdim + kl] +
					     (tyik - ty2) * ecc[(ki - 1) * kdim + kl]) / tden;
					ecw[ki * kdim + kl] =
					    ((ty1 - tyi) * ecw[ki * kdim + kl] +
					     (tyik - ty1) * ecw[(ki - 1) * kdim + kl]) / tden +
					    ecc[ki * kdim + kl];
				    }
			    }
		    }
		kjmid = (kj - 1) * kdim;
		ikmid = (order() - 1) * kdim;

		for (kl = 0; kl < kdim; kl++)
		    new_coefs.push_back(ecw[ikmid + kl] / order());
	    }

	// Update our basis_ object according to new_coefs & new_knots.
	double *ks;
	ks = &new_knots[0];
	basis_ = BsplineBasis(num_newCoefs, new_order, ks);

	if (rat) {
	    std::swap(rcoefs_, new_coefs);
	    updateCoefsFromRcoefs();
	} else {
	    std::swap(coefs_, new_coefs);
	}
	// Cleaning up.
	delete[] ecc;
	delete[] ecw;

    }

}
} // namespace Go;
