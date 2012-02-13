//===========================================================================
//                                                                           
// File: IntersectionUtils.C                                                 
//                                                                           
// Created: Mon Jan 31 17:31:37 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: IntersectionUtils.C,v 1.23 2006-04-20 10:27:46 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/IntersectionUtils.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/SurfaceCreators.h"
#include "GoTools/implicitization/BernsteinMulti.h"
#include "GoTools/implicitization/BernsteinUtils.h"
#include "GoTools/implicitization/ImplicitUtils.h"
#include "GoTools/creators/CreatorsUtils.h"

#include <fstream>


using std::vector;

namespace Go {


//===========================================================================
shared_ptr<SplineCurve>
IntersectionUtils::create1DSplineCurve(const SplineCurve& cv, int dim_id)
//===========================================================================
{
    int space_dim = cv.dimension();
    int rat = cv.rational();
    int total_dim = cv.dimension() + rat;
    ASSERT(dim_id < space_dim);

    int new_total_dim = 1 + rat;
    vector<double> coefs(cv.numCoefs()*new_total_dim);
    for (size_t ki = 0; ki < coefs.size(); ++ki) {
	coefs[ki*new_total_dim]
	    = (rat) ? cv.rcoefs_begin()[ki*total_dim+dim_id] :
	    cv.coefs_begin()[ki*total_dim+dim_id];;
	if (rat) {
	    coefs[ki*new_total_dim+1]
		= cv.rcoefs_begin()[ki*total_dim+space_dim];
	}
    }

    shared_ptr<SplineCurve>
	return_cv(new SplineCurve(cv.numCoefs(), cv.order(),
				  cv.basis().begin(), coefs.begin(), 1, (rat!=0)));

    return return_cv;
}


//===========================================================================
shared_ptr<SplineSurface>
IntersectionUtils::create1DSplineSurface(const SplineSurface& sf, int dim_id)
//===========================================================================
{
    int space_dim = sf.dimension();
    int rat = sf.rational();
    int total_dim = sf.dimension() + rat;
    ASSERT(dim_id < space_dim);

    int new_total_dim = 1 + rat;
    int num_coefs = sf.numCoefs_u()*sf.numCoefs_v();
    vector<double> coefs(sf.numCoefs_u()*sf.numCoefs_v()*new_total_dim);
    for (int ki = 0; ki < num_coefs; ++ki) {
	coefs[ki*new_total_dim]
	    = (rat) ? sf.rcoefs_begin()[ki*total_dim+dim_id] :
	    sf.coefs_begin()[ki*total_dim+dim_id];
	if (rat) {
	    coefs[ki*new_total_dim+1]
		= sf.rcoefs_begin()[ki*total_dim+space_dim];
	}
    }

    shared_ptr<SplineSurface>
	return_sf(new SplineSurface(sf.numCoefs_u(), sf.numCoefs_v(),
				    sf.order_u(), sf.order_v(),
				    sf.basis_u().begin(),
				    sf.basis_v().begin(),
				    coefs.begin(), 1, (rat!=0)));

    return return_sf;
}


//===========================================================================
shared_ptr<SplineCurve>
IntersectionUtils::splineCurveProduct(vector<shared_ptr<SplineCurve> >& cv,
				      Alg2DElem term)
//===========================================================================
{
    ASSERT(cv.size() > 0);
    vector<double> coefs(cv[0]->numCoefs(), 1.0); // We start by
						  // defining the cv
						  // as the identity.
    shared_ptr<SplineCurve> return_cv;
    int ki, kj;
    for (ki = 0; ki < 2; ++ki) {
	int degree = term.degrees_[ki];
	for (kj = 0; kj < degree; ++kj) {
	    if (return_cv.get() == 0) {
		return_cv = shared_ptr<SplineCurve>(cv[ki]->clone());
	    } else {
		shared_ptr<SplineCurve>
		    mult_cv(CurveCreators::
			    multCurveWithFunction(*return_cv, *cv[ki]));
		*return_cv = *mult_cv;
	    }
	}
    }

    if (return_cv.get() == 0)
	return_cv = shared_ptr<SplineCurve>
	    (new SplineCurve(cv[0]->numCoefs(), cv[0]->order(),
			     cv[0]->basis().begin(), coefs.begin(), 1));

    for (ki = 0; ki < return_cv->numCoefs(); ++ki)
	return_cv->coefs_begin()[ki] *= term.factor_;

    return return_cv;
}


//===========================================================================
shared_ptr<SplineSurface> IntersectionUtils::
splineSurfaceProduct(vector<shared_ptr<SplineSurface> >& sf, Alg3DElem term)
//===========================================================================
{
    ASSERT(sf.size() > 0);
    shared_ptr<SplineSurface> return_sf;
    int ki, kj;
    for (ki = 0; ki < 3; ++ki) {
	int degree = term.degrees_[ki];
	for (kj = 0; kj < degree; ++kj) {
	    if (return_sf.get() == 0) {
		return_sf = shared_ptr<SplineSurface>(sf[ki]->clone());
	    } else {
		shared_ptr<SplineSurface> mult_sf
		    (SurfaceCreators::mult1DSurfaces(*return_sf, *sf[ki]));
		*return_sf = *mult_sf;
	    }
	}
    }

    if (return_sf.get() == 0) { // If the degree is 0 in all directions
			        // we return a constant.
	vector<double> coefs(sf[0]->numCoefs_u()*sf[0]->numCoefs_v(), 1.0);
	return_sf = shared_ptr<SplineSurface>
	    (new SplineSurface(sf[0]->numCoefs_u(), sf[0]->numCoefs_v(),
			       sf[0]->order_u(), sf[0]->order_v(),
			       sf[0]->basis_u().begin(),
			       sf[0]->basis_v().begin(),
			       coefs.begin(), 1));
    }

    vector<double>::iterator iter = return_sf->coefs_begin();
    while (iter != return_sf->coefs_end()) {
	*iter *= term.factor_;
	++iter;
    }

    return return_sf;
}


//===========================================================================
shared_ptr<SplineCurve>
IntersectionUtils::insertCvInAlgcv(const SplineCurve& cv,
				   AlgObj2DInt* alg_obj2d_int)
//===========================================================================
{
    vector<shared_ptr<SplineCurve> > spline_cv_parts(2);
    int ki;
    for (ki = 0; ki < 2; ++ki) // We extract x- and y-part of the cv.
	spline_cv_parts[ki]
	    = shared_ptr<SplineCurve>(create1DSplineCurve(cv, ki));

    shared_ptr<SplineCurve> sum_cv;
    for (ki = 0; ki < alg_obj2d_int->numTerms(); ++ki) {
	shared_ptr<SplineCurve> part_cv
	    (splineCurveProduct(spline_cv_parts, alg_obj2d_int->term(ki)));
	try {
	    sum_cv = (ki == 0) ? part_cv : 
		shared_ptr<SplineCurve>(GeometryTools::curveSum(*sum_cv, 1.0,
						 *part_cv, 1.0));
	} catch (...) {
	    THROW("Failed adding curves");
	}
    }

    return shared_ptr<SplineCurve>(sum_cv->clone());
}


//===========================================================================
shared_ptr<SplineSurface>
IntersectionUtils::insertSfInAlgsf(const SplineSurface& sf,
				   AlgObj3DInt* alg_obj3d_int)
//===========================================================================
{
    int dim = sf.dimension();
    ASSERT(dim == 3);
    if (!alg_obj3d_int->usingPowerBasis()) {
	BernsteinTetrahedralPoly impl;
	BaryCoordSystem3D bc;
	alg_obj3d_int->getImplicit(impl, bc);
	return insertSfInImplObj(sf, impl, bc);
    } else {
	// We give the sf a bezier representation.
	vector<SplineSurface> patches;
	vector<shared_ptr<SplineSurface> > mult_patches;
	GeometryTools::splitSurfaceIntoPatches(sf, patches);
	for (size_t ki = 0; ki < patches.size(); ++ki) {
	    vector<BernsteinMulti> bern_mult(dim);
	    spline_to_bernstein(patches[ki], bern_mult);
	    vector<double> init_bern(1, 1.0); // Initialized to 1.0.
	    BernsteinMulti final_sum;
	    for (int kj = 0; kj < alg_obj3d_int->numTerms(); ++kj) {
		Alg3DElem alg_3d_elem = alg_obj3d_int->term(kj);
		BernsteinMulti final_mult(0, 0, init_bern);
		for (int kg = 0; kg < dim; ++kg) {
		    int order = alg_3d_elem.degrees_[kg];
		    for (int kh = 0; kh < order; ++kh) {
			final_mult *= bern_mult[kg];
		    }
		}
		final_mult *= alg_3d_elem.factor_;
		if (kj == 0) {
		    final_sum = final_mult;
		} else {
		    final_sum += final_mult;
		}
	    }

	    // We then must create the final bezier spline surface.
	    int nmb_u = final_sum.degreeU() + 1;
	    int nmb_v = final_sum.degreeV() + 1;
	    vector<double> new_knots_u(nmb_u, patches[ki].startparam_u());
	    new_knots_u.insert(new_knots_u.end(),
			       nmb_u, patches[ki].endparam_u());
	    vector<double> new_knots_v(nmb_v, patches[ki].startparam_v());
	    new_knots_v.insert(new_knots_v.end(),
			       nmb_v, patches[ki].endparam_v());
	    mult_patches.push_back
		(shared_ptr<SplineSurface>
		 (new SplineSurface(nmb_u, nmb_v, nmb_u, nmb_v,
				    new_knots_u.begin(), new_knots_v.begin(),
				    final_sum.coefsBegin(), 1,
				    patches[ki].rational())));
	}

	// We then must recreate the B-spline surface.  The inner
	// continuity is taken care of by input basises (not the basis
	// being used, only used to calc the correct 'order -
	// inner_mult').
	return GeometryTools::joinPatches(mult_patches, sf);
    }
}


//===========================================================================
shared_ptr<SplineSurface>
IntersectionUtils::insertSfInAlgsf2(const SplineSurface& sf,
				    AlgObj3DInt* alg_obj3d_int)
//===========================================================================
{
    if (!alg_obj3d_int->usingPowerBasis()) {
	BernsteinTetrahedralPoly impl;
	BaryCoordSystem3D bc;
	alg_obj3d_int->getImplicit(impl, bc);
	return insertSfInImplObj(sf, impl, bc);
    } else {
	vector<shared_ptr<SplineSurface> > spline_sf_parts(3);
	int ki;
	for (ki = 0; ki < 3; ++ki) // We extract x- and y-part of the sf.
	    spline_sf_parts[ki]
		= shared_ptr<SplineSurface>(create1DSplineSurface(sf, ki));

	shared_ptr<SplineSurface> sum_sf;
	for (ki = 0; ki < alg_obj3d_int->numTerms(); ++ki) {
	    shared_ptr<SplineSurface> part_sf
		(splineSurfaceProduct(spline_sf_parts,
				      alg_obj3d_int->term(ki)));

	    if (ki == 0) {
		sum_sf = part_sf;
	    } else {
		// If the sfs do not live in the same spline space we
		// must raise order.
		int order_diff_u = sum_sf->order_u() - part_sf->order_u();
		int order_diff_v = sum_sf->order_v() - part_sf->order_v();
		if (order_diff_u > 0 || order_diff_v > 0) {
		    part_sf->raiseOrder(std::max(order_diff_u, 0),
					std::max(order_diff_v, 0));
		}
		if (order_diff_u < 0 || order_diff_v < 0) {
		    sum_sf->raiseOrder(std::max(-order_diff_u, 0),
				       std::max(-order_diff_v, 0));
		}
		// Spline spaces are now equal, allowing us to add coefs.
		vector<double>::iterator sum_sf_iter = sum_sf->coefs_begin();
		vector<double>::const_iterator part_sf_iter
		    = part_sf->coefs_begin();
		while (sum_sf_iter != sum_sf->coefs_end()) {
		    *sum_sf_iter += *part_sf_iter;
		    ++sum_sf_iter;
		    ++part_sf_iter;
		}
	    }
	}

	return shared_ptr<SplineSurface>(sum_sf->clone());
    }
}


//===========================================================================
shared_ptr<SplineSurface>
IntersectionUtils::insertSfInImplObj(const SplineSurface& spline_sf,
				     const BernsteinTetrahedralPoly& impl,
				     const BaryCoordSystem3D& bc)
//===========================================================================
{
    // First convert the spline surface to barycentric coordinates
    SplineSurface surf_bc;
    cart_to_bary(spline_sf, bc, surf_bc);

    // Split into patches
    vector<SplineSurface> patches;
    GeometryTools::splitSurfaceIntoPatches(surf_bc, patches);
    int npatches = (int)patches.size();

    // Loop through all patches
    vector<Array<BernsteinMulti, 4> > beta(npatches);
    vector<BernsteinMulti> betatmp;
    BernsteinMulti eval_func;
    Array<BernsteinMulti, 1> eval_array;
    vector<shared_ptr<SplineSurface> > scalar_patches(npatches);
    SplineSurface scalar_tmp;
    for (int i = 0; i < npatches; ++i) {

	// Convert to BernsteinMultis. If spline_sf is rational, only
	// the first four components are used.
	spline_to_bernstein(patches[i], betatmp);
	for (int j = 0; j < 4; ++j) {
	    beta[i][j] = betatmp[j];
	}

	// Evaluate the implicit function on the current patch. We are
	// taking advantage of the fact that an evaluation operator is
	// defined for BernsteinTetrhedralPolys (i.e., operator()),
	// and takes a template argument of the type Array<T, 4>. Here
	// the typename T equals BernsteinMulti.
	eval_func = impl(beta[i]);
	eval_array[0] = eval_func;

	// Convert back to a 1D spline. We must pretend that
	// eval_array represents a rational spline since we ignored
	// denominators in the evaluation above. (If rational==true,
	// the last component of eval_array is interpreted as a
	// denominator.)
	bool rational = false;
	bernsteinToSpline(eval_array, rational, scalar_tmp);
	double u1 = patches[i].startparam_u();
	double u2 = patches[i].endparam_u();
	double v1 = patches[i].startparam_v();
	double v2 = patches[i].endparam_v();
	scalar_tmp.setParameterDomain(u1, u2, v1, v2);
	scalar_patches[i] = shared_ptr<SplineSurface>(scalar_tmp.clone());

    }

    // Join the patches again. spline_sf is used to define the spline
    // space.
    shared_ptr<SplineSurface> joined_sf
	= GeometryTools::joinPatches(scalar_patches, spline_sf);

    return joined_sf;

}


//===========================================================================
double IntersectionUtils::
distImplRepresentationCompFunction(const SplineSurface& spline_sf,
				   const BernsteinTetrahedralPoly& impl,
				   const BaryCoordSystem3D& bc,
				   const SplineSurface& comp_1d_sf,
				   double upar, double vpar)
//===========================================================================
{
    Point spline_pt = spline_sf.ParamSurface::point(upar, vpar);
    Vector3D c_pt(spline_pt[0], spline_pt[1], spline_pt[2]);
    Vector4D bc_pt = bc.cartToBary(c_pt);
    double impl_val = impl(bc_pt);
    Point comp_pt = comp_1d_sf.ParamSurface::point(upar, vpar);
    double dist = fabs(impl_val - comp_pt[0]);
//     double scale = std::max(fabs(impl_val), fabs(comp_pt[0]));
//     double rel_dist = (scale == 0.0) ? fabs(impl_val - comp_pt[0]) :
// 	fabs(impl_val - comp_pt[0])/scale;

    return dist;
}


//===========================================================================


} // end namespace Go
