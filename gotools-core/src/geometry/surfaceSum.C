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

using std::vector;
using std::setprecision;
using std::endl;
using std::pair;
using std::make_pair;

namespace Go
{

  shared_ptr<SplineSurface>
  GeometryTools::surfaceSum(const SplineSurface& sf1, double fac1,
	     const SplineSurface& sf2, double fac2, double num_tol)

    //********************************************************************
    // Addition of two signed SplineSurfaces, i.e. this function can
    // also be used for subtraction. The surfaces is assumed to live on
    // the same parameter domain, but may have different knot vectors.
    //********************************************************************
  {
    // Check input
    ALWAYS_ERROR_IF(fabs(sf1.startparam_u() - sf2.startparam_u()) > num_tol ||
		    fabs(sf1.endparam_u() - sf2.endparam_u()) > num_tol ||
		    fabs(sf1.startparam_v() - sf2.startparam_v()) > num_tol ||
		    fabs(sf1.endparam_v() - sf2.endparam_v()) > num_tol,
		    "Inconsistent parameter domain.");

    // For the time being
    if (sf1.rational() || sf2.rational()) {
	THROW("Sum of rational surfaces is not implemented");
    }

    // Make copy of surfaces
    vector<shared_ptr<SplineSurface> > surfaces;
    surfaces.reserve(2);
    shared_ptr<SplineSurface> sf;
// #ifdef _MSC_VER
//     sf = shared_ptr<SplineSurface>(dynamic_cast<SplineSurface*>(sf1.clone()));
// #else
    sf = shared_ptr<SplineSurface>(sf1.clone());
// #endif
    surfaces.push_back(sf);
// #ifdef _MSC_VER
//     sf = shared_ptr<SplineSurface>(dynamic_cast<SplineSurface*>(sf2.clone()));
// #else
    sf = shared_ptr<SplineSurface>(sf2.clone());
// #endif
    surfaces.push_back(sf);

    // Make sure that the surfaces live on the same knot vector
    GeometryTools::unifySurfaceSplineSpace(surfaces, num_tol);

    // Add signed coefficients
    vector<double> coefs;
    int nmb_coefs_u = surfaces[0]->numCoefs_u();
    int nmb_coefs_v = surfaces[0]->numCoefs_v();
    int dim = surfaces[0]->dimension();
    coefs.resize(dim*nmb_coefs_u*nmb_coefs_v);
    int ki;
    std::vector<double>::iterator s1 = surfaces[0]->coefs_begin();
    std::vector<double>::iterator s2 = surfaces[1]->coefs_begin();
    for (ki=0; ki<dim*nmb_coefs_u*nmb_coefs_v; ki++)
      coefs[ki] = fac1*s1[ki] + fac2*s2[ki];

    // Create output curve
    shared_ptr<SplineSurface>
	surfacesum(new SplineSurface(nmb_coefs_u, nmb_coefs_v,
				     surfaces[0]->order_u(),
				     surfaces[0]->order_v(),
				     surfaces[0]->basis_u().begin(),
				     surfaces[0]->basis_v().begin(),
				     &coefs[0], dim, false));

    return surfacesum;
  }
}
