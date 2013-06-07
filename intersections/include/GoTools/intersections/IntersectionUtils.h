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

#ifndef _INTERSECTIONUTILS_H
#define _INTERSECTIONUTILS_H


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/intersections/AlgObj2DInt.h"
#include "GoTools/intersections/AlgObj3DInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/implicitization/BernsteinMulti.h"
#include <vector>


namespace Go {


/// Various functions related to the intersection algorithms
namespace IntersectionUtils {


// Extract 1-dim part of the cv (dim_id == 0 => implies extract
// x-coord part etc).
shared_ptr<SplineCurve>
create1DSplineCurve(const SplineCurve& cv, int dim_id);

shared_ptr<SplineSurface>
create1DSplineSurface(const SplineSurface& sf, int dim_id);

// Return the splinecurve as 1d product of it's nth degree
// dimensional parts (x^(nx)*y^(ny)*z^(nz))
shared_ptr<SplineCurve>
splineCurveProduct(std::vector<shared_ptr<SplineCurve> >& cv,
		   Alg2DElem term);

shared_ptr<SplineSurface>
splineSurfaceProduct(std::vector<shared_ptr<SplineSurface> >& sf,
		     Alg3DElem term);

shared_ptr<SplineCurve>
insertCvInAlgcv(const SplineCurve& cv, AlgObj2DInt* alg_obj2d_int);

// In this first version we transform to BernsteinMulti, for which
// multiplication is believed to be faster and more accurate.
shared_ptr<SplineSurface>
insertSfInAlgsf(const SplineSurface& sf, AlgObj3DInt* alg_obj3d_int);

// Outdated version using multiplication between SplineSurface
// objects.  Probably to be removed.
shared_ptr<SplineSurface>
insertSfInAlgsf2(const SplineSurface& sf, AlgObj3DInt* alg_obj3d_int);

// In order to insert the spline-sf into the equation we first put
// it on the input barycentric coordinate system.  The return sf
// is the 1D-surface resulting from converting (x,y,z) to
// barycentric coordinated and evaluating in the impl function.
shared_ptr<SplineSurface>
insertSfInImplObj(const SplineSurface& spline_sf,
		  const BernsteinTetrahedralPoly& impl,
		  const BaryCoordSystem3D& bc);

// Compute the distance between the original component
// representation to the composite representation. @@sbr Remove
// when stable.  comp_1d_sf = impl(spline_sf_x, spline_sf_y,
// spline_sf_z);
double
distImplRepresentationCompFunction(const SplineSurface& spline_sf,
				   const BernsteinTetrahedralPoly& impl,
				   const BaryCoordSystem3D& bc,
				   const SplineSurface& comp_1d_sf,
				   double upar, double vpar);
    

} // end namespace IntersectionUtils


} // end namespace Go


#endif // _INTERSECTIONUTILS_H

