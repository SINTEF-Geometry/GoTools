//===========================================================================
//                                                                           
// File: IntersectionUtils.h                                                 
//                                                                           
// Created: Mon Jan 31 17:28:49 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: IntersectionUtils.h,v 1.10 2006-02-23 14:46:04 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

