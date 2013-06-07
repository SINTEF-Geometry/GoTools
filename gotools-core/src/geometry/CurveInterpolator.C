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

#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/SplineCurve.h"

using std::vector;

namespace Go
{
  SplineCurve* 
  CurveInterpolator::regularInterpolation(const BsplineBasis& basis,
					    vector<double>& par,
					    vector<double>& points,
					    int dimension,
					    bool rational,
					    vector<double>& weights)
  {
    // Check input
    ASSERT(par.size() == points.size()/dimension);
    ASSERT(basis.numCoefs() == (int)par.size());

    vector<double> points2;
    if (rational)
      {
	ASSERT(weights.size() == points.size()/dimension);

	// Include weight information
	// First get weights in interpolation points
	shared_ptr<SplineCurve> denom = 
	  shared_ptr<SplineCurve>(new SplineCurve(basis, 
						  weights.begin(), 1,
						  false));
	vector<double> wgtval;
	denom->gridEvaluator(wgtval, par);
	size_t nmb_pnt = par.size();
	points2.reserve(nmb_pnt*(dimension+1));
	for (size_t kr=0; kr<nmb_pnt; ++kr)
	  {
	    for (int kh=0; kh<dimension; kh++)
	      points2.push_back(points[kr*dimension+kh]*wgtval[kr]);
	    points2.push_back(wgtval[kr]);
	  }
      }
    else
      points2 = points;

    // Interpolate curve
    vector<int> tg_idx;
    vector<double> tg_pnt;
    vector<double> coefs;
    SplineInterpolator interpolator;
    interpolator.setBasis(basis);
    interpolator.interpolate(par, points2, tg_idx, tg_pnt, coefs);

    // Make curve
    SplineCurve* crv = new SplineCurve(basis, coefs.begin(), dimension,
				       rational);
    return crv;
  }

} // namespace Go
