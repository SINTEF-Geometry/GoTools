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

#include "GoTools/geometry/GeometryTools.h"
#include <memory>

//***************************************************************************
//
// Implementation file of the free function curveSum defined in
// GeometryTools.h/
//
//***************************************************************************

using std::vector;
using std::max;
using std::min;

namespace Go
{

  shared_ptr<SplineCurve>
  GeometryTools::curveSum(const SplineCurve& crv1, double fac1,
	   const SplineCurve& crv2, double fac2, double num_tol)

    //********************************************************************
    // Addition of two signed SplineCurves, i.e. this function can
    // also be used for subtraction. The curves is assumed to live on
    // the same parameter domain, but may have different knot vectors.
    //********************************************************************
  {
    // Check input
    ALWAYS_ERROR_IF(fabs(crv1.startparam() - crv2.startparam()) > num_tol ||
		fabs(crv1.endparam() - crv2.endparam()) > num_tol,
		"Inconsistent parameter domain.");

    // For the time being
    if (crv1.rational() || crv2.rational()) {
	THROW("Sum of rational curves is not impelemented");
    }

    // Make copy of curves
    vector<shared_ptr<SplineCurve> > curves;
    curves.reserve(2);
    shared_ptr<SplineCurve> cv;
// #ifdef _MSC_VER
//     cv = shared_ptr<SplineCurve>(dynamic_cast<SplineCurve*>(crv1.clone()));
// #else
    cv = shared_ptr<SplineCurve>(crv1.clone());
// #endif
    curves.push_back(cv);
// #ifdef _MSC_VER
//     cv = shared_ptr<SplineCurve>(dynamic_cast<SplineCurve*>(crv2.clone()));
// #else
    cv = shared_ptr<SplineCurve>(crv2.clone());
// #endif
    curves.push_back(cv);

    // Make sure that the curves live on the same knot vector
//     double tol = 0.00001;
    try {
	unifyCurveSplineSpace(curves, num_tol);
    } catch (...) {
	THROW("Failed unifying spline spaces!");
    }

    // Add signed coefficients
    vector<double> coefs;
    int nmb_coefs = curves[0]->numCoefs();
    int dim = curves[0]->dimension();
    coefs.resize(dim*nmb_coefs);
    int ki;
    std::vector<double>::iterator c1 = curves[0]->coefs_begin();
    std::vector<double>::iterator c2 = curves[1]->coefs_begin();
    for (ki=0; ki<dim*nmb_coefs; ki++)
      coefs[ki] = fac1*c1[ki] + fac2*c2[ki];

    // Create output curve
    shared_ptr<SplineCurve>
	curvesum(new SplineCurve(nmb_coefs, 
				 curves[0]->order(),
				 curves[0]->basis().begin(),
				 &coefs[0], 
				 dim,
				 false));
    return curvesum;
  }
}

