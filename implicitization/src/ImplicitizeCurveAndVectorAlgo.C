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

#include "GoTools/implicitization/ImplicitizeCurveAndVectorAlgo.h"
#include "GoTools/implicitization/ImplicitizeSurfaceAlgo.h"
// #include "GoTools/implicitization/ImplicitUtils.h"
// #include "GoTools/geometry/GeometryTools.h"


using namespace std;


namespace Go {


//==========================================================================
int ImplicitizeCurveAndVectorAlgo::perform()
//==========================================================================
{
    // Make a ruled surface

    // Number of coefficients
    int ncoefsu = crv_.numCoefs();
    int ncoefsv = 2;

    // Orders
    int ordu = crv_.order();
    int ordv = 2;

    // Knot vectors
    typedef vector<double>::iterator iter;
    iter knotsu_begin = crv_.basis().begin();
    double knotsv_begin[] = { 0.0, 0.0, 1.0, 1.0 };

    // Rationality and dimension
    bool rational = crv_.rational();
    int dim = 3;
    int effdim = (rational ? dim+1 : dim);

    // Now the tricky part - the coefficients

    // First we need a scale. (We normalize the vector pt_.)
    BoundingBox box = crv_.boundingBox();
    Point diagonal = box.high() - box.low();
    double scale = diagonal.length();
    pt_.normalize();

    vector<double> coefs(effdim * ncoefsu * ncoefsv);
    iter curr = (rational) ? crv_.rcoefs_begin() : crv_.coefs_begin();
    iter down = coefs.begin(); // First row of coefs ("down")
    iter up = coefs.begin() + effdim * ncoefsu; // Second row of coefs
						// ("up")
    for (int i = 0; i < ncoefsu; ++i) {
	double w = (rational ? *(curr + dim) : 1.0);
	Point tmppt(curr, curr+dim);
	// Set a down-point
	tmppt -= w * 0.5 * scale * pt_;
	copy(tmppt.begin(), tmppt.end(), down);
	if (rational) {
  	    *(down + dim) = w;
	}
	down += effdim;
	// Set an up-point
	tmppt += w * scale * pt_;
	copy(tmppt.begin(), tmppt.end(), up);
	if (rational) {
	    *(up + dim) = w;
	}
	up += effdim;

	curr += effdim;
    }

    SplineSurface surf(ncoefsu, ncoefsv, ordu, ordv,
		       knotsu_begin, knotsv_begin, coefs.begin(), dim,
		       rational);

    // Run surface algorithm
    ImplicitizeSurfaceAlgo algo(surf, deg_);
    algo.setTolerance(tol_);
    algo.perform();
    algo.getResultData(implicit_, bc_, sigma_min_);

    return 0;
}


//==========================================================================


} // namespace Go
