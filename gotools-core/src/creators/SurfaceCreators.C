//===========================================================================
//                                                                           
// File: SurfaceCreators.C                                                 
//                                                                           
// Created: Thu Feb 21 09:33:00 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: SurfaceCreators.C,v 1.11 2009-05-13 07:30:32 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/creators/SurfaceCreators.h"

#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineUtils.h"
//#include "sisl.h"
#include "GoTools/creators/HermiteAppC.h"
#include "GoTools/creators/SmoothTransition.h"
#include "GoTools/creators/HermiteAppS.h"
#include "GoTools/creators/CoonsPatchGen.h"
//#include "newmat.h"
#include "GoTools/utils/LUDecomp.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineDebugUtils.h"

#include <fstream>


using namespace Go;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::vector;
using std::max;
using std::min;



namespace {

/// We make sure that the parameter and space curves form loops
/// (assuming that they nearly do).
/// \param param_bd_curves parametric curves.
/// \param space_bd_curves the space curves corresponding to param_bd_curves.
/// \param sfs the parametric surfaces on which the input curves lie.
void transformToLoops(std::vector<std::shared_ptr<SplineCurve> >& param_bd_curves,
		      std::vector<std::shared_ptr<SplineCurve> >& space_bd_curves,
		      const std::vector<std::shared_ptr<const ParamSurface> >& sfs);

};// end anonymous namespace


//===========================================================================
std::shared_ptr<SplineSurface>
SurfaceCreators::createSmoothTransition(const vector<shared_ptr<const ParamSurface> >& surfs,
					const vector<shared_ptr<const CurveOnSurface> >& int_cvs,
					double dist_0, double dist_1, double epsge,
					vector<shared_ptr<SplineCurve> >& trim_crvs)
//===========================================================================
{
    ALWAYS_ERROR_IF(surfs.size() != 2,
		    "Expecting input of two surfaces.");
    ALWAYS_ERROR_IF(int_cvs.size() != 2,
		    "Expecting input of two boundary curves.");
    shared_ptr<const SplineCurve> space_cv
	(dynamic_pointer_cast<const SplineCurve, const ParamCurve>(int_cvs[0]->spaceCurve()));
    shared_ptr<const SplineCurve> param_cv1
	(dynamic_pointer_cast<const SplineCurve, const ParamCurve>(int_cvs[0]->parameterCurve()));
    shared_ptr<const SplineCurve> param_cv2
	(dynamic_pointer_cast<const SplineCurve, const ParamCurve>(int_cvs[1]->parameterCurve()));
    // @@sbr Check direction of curves! Suspecting the other should be turned...
    ALWAYS_ERROR_IF(param_cv1.get() == 0 || param_cv2.get() == 0 || space_cv.get() == 0,
		    "Missing parameter or space curve.");

    SmoothTransition smooth_trans(space_cv, param_cv1, param_cv2, surfs[0], surfs[1],
				  dist_0, dist_1, epsge);
    vector<double> initpars;
    int order = space_cv->order();
    int nb_coef = space_cv->numCoefs();
    vector<double>::const_iterator knots = space_cv->basis().begin();
    initpars.push_back(knots[order-1]);
    for (int kj = order; kj <= nb_coef; ++kj)
	if (knots[kj] > initpars[initpars.size()-1])
	    initpars.push_back(knots[kj]);
    vector<int> dims(6);
    dims[0] = dims[2] = dims[3] = dims[5] = 3;
    dims[1] = dims[4] = 2;
    HermiteAppS approximator(&smooth_trans, &initpars[0],
			     (int)initpars.size(), epsge, epsge, dims);

    approximator.refineApproximation();
    vector<shared_ptr<SplineCurve> > smooth_trans_cvs = approximator.getCurves();

    // Finally we use lofting to create a smooth transition between surfs.
    SplineSurface* smooth_trans_sf;
    // @@sbr We ought to choose parameter values based on geometry.
    vector<shared_ptr<SplineCurve> > param_bd_cvs, space_bd_cvs, cross_cvs;
    vector<double> params;
    vector<int> cross_index;
    for (int ki = 0; ki < 2; ++ki) {
	param_bd_cvs.push_back(smooth_trans_cvs[3*ki+1]);
	space_bd_cvs.push_back(smooth_trans_cvs[3*ki]);
	cross_cvs.push_back(smooth_trans_cvs[3*ki+2]);
	params.push_back((double)ki);
	cross_index.push_back(ki);
    }
    // If input space curves are loops, we make sure that the boundary curves form also are.
    double loop_tol = epsge;
    if ((space_cv->ParamCurve::point(space_cv->startparam())).dist
	((space_cv->ParamCurve::point(space_cv->endparam()))) < loop_tol) {
	// Input curves expected to share basis.
 	transformToLoops(param_bd_cvs, space_bd_cvs, surfs);

	// We next must make sure that the end tangents correspond. We do this by using Coon's Patch.
	vector<double> pts(6);
	Point start_pt = 0.5*(space_bd_cvs[0]->ParamCurve::point(space_bd_cvs[0]->startparam()) +
			      space_bd_cvs[0]->ParamCurve::point(space_bd_cvs[0]->endparam()));
	Point end_pt = 0.5*(space_bd_cvs[1]->ParamCurve::point(space_bd_cvs[1]->startparam()) +
			    space_bd_cvs[1]->ParamCurve::point(space_bd_cvs[1]->endparam()));
	copy(start_pt.begin(), start_pt.end(), pts.begin());
	copy(end_pt.begin(), end_pt.end(), pts.begin() + 3);
	Point start_tangent = 0.5*(cross_cvs[0]->ParamCurve::point(cross_cvs[0]->startparam()) +
				   cross_cvs[0]->ParamCurve::point(cross_cvs[0]->endparam()));
	Point end_tangent = 0.5*(cross_cvs[1]->ParamCurve::point(cross_cvs[1]->startparam()) +
				 cross_cvs[1]->ParamCurve::point(cross_cvs[1]->endparam()));
	end_tangent *= -1.0;
	// We create a curve describing the closed edge.
	SplineInterpolator int_pol;
	int_pol.setHermiteConditions(start_tangent, end_tangent);
	vector<double> coefs;
	int_pol.interpolate(2, 3, &params[0], &pts[0], coefs);
	vector<double> knots(8);
	knots[0] = knots[1] = knots[2] = knots[3] = params[0];
	knots[4] = knots[5] = knots[6] = knots[7] = params[1];
	vector<shared_ptr<SplineCurve> > common_bd_cv(2);
	common_bd_cv[0] = shared_ptr<SplineCurve>
	    (new SplineCurve(4, 4, &knots[0], &coefs[0], 3));
	common_bd_cv[1] = shared_ptr<SplineCurve>
	    (new SplineCurve(4, 4, &knots[0], &coefs[0], 3));

	vector<shared_ptr<SplineCurve> > all_bd_cvs(4), all_cross_cvs(4);
	for (int ki = 0; ki < 2; ++ki) {
	    all_bd_cvs[2*ki] = space_bd_cvs[ki];
	    all_bd_cvs[2*ki+1] = common_bd_cv[ki];
	    all_cross_cvs[2*ki] = cross_cvs[ki];
	    all_cross_cvs[2*ki+1] = shared_ptr<SplineCurve>();
	}
	all_bd_cvs[2]->reverseParameterDirection();
	all_bd_cvs[3]->reverseParameterDirection();
	all_cross_cvs[2]->reverseParameterDirection();
	CoonsPatchGen::fixCrossEndPts(all_bd_cvs, all_cross_cvs);
	all_bd_cvs[2]->reverseParameterDirection(); // Last bd_cv is not used as we are lofting.
	all_bd_cvs[3]->reverseParameterDirection(); // @@sbr Not necessary.
	all_cross_cvs[2]->reverseParameterDirection();

	vector<double>::iterator iter; // We must reverse direction of cross tangent.
	for (iter = cross_cvs[1]->coefs_begin(); iter != cross_cvs[1]->coefs_end(); ++iter)
	    *iter *= -1.0;
    }

    smooth_trans_sf =
	CoonsPatchGen::loftSurface(space_bd_cvs.begin(), params.begin(), 2,
				   cross_cvs.begin(), cross_index);

    // Transition is represented as a SplineSurface.
    trim_crvs.push_back(param_bd_cvs[0]);
    trim_crvs.push_back(space_bd_cvs[0]);
    space_bd_cvs[1]->reverseParameterDirection();
    param_bd_cvs[1]->reverseParameterDirection();
    trim_crvs.push_back(param_bd_cvs[1]);
    trim_crvs.push_back(space_bd_cvs[1]);

    return std::shared_ptr<SplineSurface>(smooth_trans_sf);
}


//===========================================================================
std::shared_ptr<SplineSurface>
SurfaceCreators::mult1DSurfaces(const SplineSurface& sf1,
				const SplineSurface& sf2)
//===========================================================================
{
    // Expecting the sfs to share spline space.
    vector<SplineSurface> sf1_patches, sf2_patches;
    splitSurfaceIntoPatches(sf1, sf1_patches);
    splitSurfaceIntoPatches(sf2, sf2_patches);

    vector<shared_ptr<SplineSurface> > mult_patches;
    for (size_t ki = 0; ki < sf1_patches.size(); ++ki) {
	mult_patches.push_back
	    (shared_ptr<SplineSurface>(mult1DBezierPatches(sf1_patches[ki], sf2_patches[ki])));
    }

#ifdef CREATORS_DEBUG
    double upar = 0.935*sf1_patches[0].startparam_u() + (1-0.935)*sf1_patches[0].endparam_u();
    double vpar = 0.915*sf1_patches[0].startparam_v() + (1-0.915)*sf1_patches[0].endparam_v();
    Point sf1_pt = sf1_patches[0].ParamSurface::point(upar, vpar);
    Point sf2_pt = sf2_patches[0].ParamSurface::point(upar, vpar);
    double mult = sf1_pt*sf2_pt;
    Point mult_pt = mult_patches[0]->ParamSurface::point(upar, vpar);
    std::cout << "Correct product: " << mult << ", calulcated roduct: " << mult_pt[0] << std::endl;
#endif // CREATORS_DEBUG

    // We then must recreate the B-spline surface.
    std::shared_ptr<SplineSurface> joined_patches = joinPatches(mult_patches, sf1);

    return joined_patches;
}

#ifdef __BORLANDC__
// 'pascal' is a keyword in Borland C++ (declares a variable or a function using
// a Pascal-style naming convention).
#define pascal pascal__
#endif

//===========================================================================
std::shared_ptr<SplineSurface>
SurfaceCreators::mult1DBezierPatches(const SplineSurface& patch1,
				     const SplineSurface& patch2)
//===========================================================================
{
    // @@sbr This should be fixed shortly. Nothing more than separating the
    // spatial and rational components.
    ASSERT(!patch1.rational() && !patch2.rational());

    // We should of course also check the actual knots, but why bother (trusting the user).
    //     ASSERT(basis1_u.numCoefs() == basis2_u.numCoefs() &&
    // 	   basis1_v.numCoefs() == basis2_v.numCoefs() &&
    // 	   basis1_u.order() == basis2_u.order() &&
    // 	   basis1_v.order() == basis2_v.order());
    ASSERT((patch1.dimension() == 1) && (patch2.dimension() == 1));
    // @@sbr Suppose we could allow for differing orders (but equal parameter domain).

    // Ported from SISL routine s6multsfs().
    int order = max(2*(patch1.order_u() - 1) + 1, 2*(patch1.order_v() - 1) + 1);;
    vector<double> pascal((order+1)*(order+2)/2, 0.0); // Binomial coefficients (Pascal's triangle)
    int ki, kj;
    vector<double>::iterator psl1;     /* Pointer used in Pascals triangle */
    vector<double>::iterator psl2;     /* Pointer used in Pascals triangle */
    for(ki = 0, psl2 = pascal.begin(); ki <= order ; ki++, psl1 = psl2, psl2 += ki) {
	psl2[0] = 1.0;

	for(kj = 1; kj < ki; kj++)
	    psl2[kj] = psl1[kj-1] + psl1[kj];

	psl2[ki] = 1.0;
    }

    int order11 = patch1.order_u();
    int order12 = patch1.order_v();
    int order21 = patch2.order_u();
    int order22 = patch2.order_v();

    vector<double>::const_iterator c1 = patch1.coefs_begin();
    vector<double>::const_iterator c2 = patch2.coefs_begin();

    //     vector<double> mult_coefs(patch1.numCoefs_u()*patch1.numCoefs_v()*patch1.dimension(), 0.0);
    int p1,p2,r1,r2;
    int kgrad11 = order11-1;
    int kgrad12 = order12-1;
    int kgrad21 = order21-1;
    int kgrad22 = order22-1;
    int kgrad1  =  kgrad11 + kgrad21; // Degree of mult basis functions in 1st dir.
    int kgrad2  =  kgrad12 + kgrad22;
    int kstop2 = order12+order22-1;
    int kstop1 = order11+order21-1;
    vector<double> mult_coefs(kstop1*kstop2, 0.0);
    vector<double>::const_iterator psl_kgrad11 = pascal.begin()+kgrad11*(kgrad11+1)/2;
    vector<double>::const_iterator psl_kgrad12 = pascal.begin()+kgrad12*(kgrad12+1)/2;
    vector<double>::const_iterator psl_kgrad21 = pascal.begin()+kgrad21*(kgrad21+1)/2;
    vector<double>::const_iterator psl_kgrad22 = pascal.begin()+kgrad22*(kgrad22+1)/2;
    vector<double>::const_iterator psl_kgrad1  = pascal.begin()+kgrad1 *(kgrad1 +1)/2;
    vector<double>::const_iterator psl_kgrad2  = pascal.begin()+kgrad2 *(kgrad2 +1)/2;
    vector<double>::const_iterator qsc1, qsc2;
    double tsum, sumi;
    vector<double>::iterator temp = mult_coefs.begin();
    double tdiv, t2;
    int kstop3, kstop4;

    for (p2 = 0; p2 < kstop2; ++p2)
	for (p1 = 0; p1 <kstop1; p1++, temp++) {
	    tdiv  =  psl_kgrad1[p1]*psl_kgrad2[p2];
	    kstop4  =  min(p2,kgrad12);
	    for (r2 = max(0,p2-kgrad22),tsum = 0.0; r2 <= kstop4; r2++) {
		t2  =  psl_kgrad12[r2]*psl_kgrad22[p2-r2];
		kstop3  =  min(p1,kgrad11);
		for (r1 = max(0,p1-kgrad21),sumi = 0.0,
			 qsc1 = c1+r2*order11,qsc2 = c2+(p2-r2)*order21;
		     r1 <= kstop3; r1++)
		    sumi +=  psl_kgrad11[r1]*psl_kgrad21[p1-r1]*qsc1[r1]*qsc2[p1-r1];
		tsum +=  t2*sumi;
	    }
	    tsum /=  tdiv;
	    *temp  =  tsum;
	}
    //     *order_newsurf1 = kstop1;
    //     *order_newsurf2 = kstop2; 

    // We must add knots to input basises according to new order.
    vector<double> new_knots_u, new_knots_v;
    new_knots_u.insert(new_knots_u.begin(), kstop1, patch1.startparam_u());
    new_knots_u.insert(new_knots_u.end(), kstop1, patch1.endparam_u());
    new_knots_v.insert(new_knots_v.begin(), kstop2, patch1.startparam_v());
    new_knots_v.insert(new_knots_v.end(), kstop2, patch1.endparam_v());

    // Finally we create the spline sf with the multiplied coefs.
    std::shared_ptr<SplineSurface> mult_sf(new SplineSurface(kstop1, kstop2, kstop1, kstop2,
							       new_knots_u.begin(), new_knots_v.begin(),
							       mult_coefs.begin(), patch1.dimension(),
							       patch1.rational()));
    //     SplineSurface* mult_sf = new SplineSurface(patch1.numCoefs_u(), patch1.numCoefs_v(),
    // 					       patch1.order_u(), patch1.order_v(),
    // 					       patch1.basis_u().begin(), patch1.basis_v().begin(),
    // 					       mult_coefs.begin(), patch1.dimension(),
    // 					       patch1.rational());

    return mult_sf;
}


//===========================================================================
shared_ptr<SplineSurface> SurfaceCreators::mergeRationalParts(const SplineSurface& nom_sf,
							      const SplineSurface& den_sf,
							      bool weights_in_first)
//===========================================================================
{
    ASSERT((!nom_sf.rational()) && (!den_sf.rational()));
    ASSERT(den_sf.dimension() == 1);

    int dim = nom_sf.dimension();

    // We first make sure they share spline space.
    vector<shared_ptr<SplineSurface> > sfs;
    sfs.push_back(shared_ptr<SplineSurface>(nom_sf.clone()));
    sfs.push_back(shared_ptr<SplineSurface>(den_sf.clone()));
    double knot_diff_tol = 1e-06;
    unifySurfaceSplineSpace(sfs, knot_diff_tol);

    vector<double> rcoefs;
    vector<double>::const_iterator iter = sfs[0]->coefs_begin();
    vector<double>::const_iterator riter = sfs[1]->coefs_begin();
    int num_coefs = sfs[0]->numCoefs_u()*sfs[0]->numCoefs_v();
    for (int ki = 0; ki < num_coefs; ++ki) {
	for (int kj = 0; kj < dim; ++kj) {
	    if (weights_in_first) {
		rcoefs.push_back(iter[ki*dim+kj]);
	    } else {
		rcoefs.push_back(iter[ki*dim+kj]*riter[ki]);
	    }
	}
	rcoefs.push_back(riter[ki]);
    }

    shared_ptr<SplineSurface> rat_sf(new SplineSurface
				     (sfs[0]->numCoefs_u(), sfs[0]->numCoefs_v(),
				      sfs[0]->order_u(), sfs[0]->order_v(),
				      sfs[0]->basis_u().begin(), sfs[0]->basis_v().begin(),
				      rcoefs.begin(), dim, true));

    return rat_sf;
}


//===========================================================================
shared_ptr<SplineSurface> SurfaceCreators::insertParamDomain(const SplineSurface& sf_1d)
//===========================================================================
{
    shared_ptr<SplineSurface> sf_1d_cp(sf_1d.clone());
    int dim = sf_1d.dimension();
    ASSERT(dim == 1);

    bool rat = sf_1d_cp->rational();
    // The returned object should be linear in the first two directions.
    // We create an additional 1d-sf describing the linear param space.
    vector<double> lin_knots_u(4, sf_1d_cp->startparam_u());
    lin_knots_u[2] = lin_knots_u[3] = sf_1d_cp->endparam_u();
    vector<double> lin_knots_v(4, sf_1d_cp->startparam_v());
    lin_knots_v[2] = lin_knots_v[3] = sf_1d_cp->endparam_v();
    int rdim = (rat) ? dim + 1 : dim;
    vector<double> lin_coefs_u(4, 1.0);
    lin_coefs_u[0] = lin_coefs_u[2] = lin_knots_u[0];
    lin_coefs_u[1] = lin_coefs_u[3] = lin_knots_u[2];
    shared_ptr<SplineSurface> lin_sf_u(new SplineSurface(2, 2, 2, 2,
							 lin_knots_u.begin(), lin_knots_v.begin(),
							 lin_coefs_u.begin(), 1));
    vector<double> lin_coefs_v(4*rdim, 1.0);
    lin_coefs_v[0] = lin_coefs_v[1] = lin_knots_v[0];
    lin_coefs_v[2] = lin_coefs_v[3] = lin_knots_v[2];
    shared_ptr<SplineSurface> lin_sf_v(new SplineSurface(2, 2, 2, 2,
							 lin_knots_u.begin(), lin_knots_v.begin(),
							 lin_coefs_v.begin(), 1));

    if (rat) {
	// We extract the rational part (i.e. the denominator sf) and mult it the linear parts.
	vector<shared_ptr<SplineSurface> > rat_parts = separateRationalParts(*sf_1d_cp);
	lin_sf_u = SurfaceCreators::mult1DSurfaces(*lin_sf_u, *rat_parts[1]);
	lin_sf_v = SurfaceCreators::mult1DSurfaces(*lin_sf_v, *rat_parts[1]);

	// We must then raise the order of sf_1d_cp by 1.
	rat_parts[0]->raiseOrder(1, 1);
	rat_parts[1]->raiseOrder(1, 1);
	sf_1d_cp = mergeRationalParts(*rat_parts[0], *rat_parts[1], false);
    } else {
	int raise_u = sf_1d_cp->order_u() - 2;
	int raise_v = sf_1d_cp->order_v() - 2;
	lin_sf_u->raiseOrder(raise_u, raise_v);
	lin_sf_v->raiseOrder(raise_u, raise_v);
    }

    // If not bezier we must also refine the space.
    int ik1 = sf_1d_cp->order_u();
    int ik2 = sf_1d_cp->order_v();
    int in1 = sf_1d_cp->numCoefs_u();
    int in2 = sf_1d_cp->numCoefs_v();
    if (ik1 < in1 || ik2 < in2)
      {
	vector<double> new_knots_u(sf_1d_cp->basis_u().begin() + ik1,
				   sf_1d_cp->basis_u().begin() + in1);
	vector<double> new_knots_v(sf_1d_cp->basis_v().begin() + ik2,
				   sf_1d_cp->basis_v().begin() + in2);
	lin_sf_u->insertKnot_u(new_knots_u);
	lin_sf_u->insertKnot_v(new_knots_v);
	lin_sf_v->insertKnot_u(new_knots_u);
	lin_sf_v->insertKnot_v(new_knots_v);
      }

    // Finally we create our space sf (i.e. living in a 3-dimensional env).
    vector<double> all_coefs;
    int coefs_size = sf_1d_cp->numCoefs_u()*sf_1d_cp->numCoefs_v();
    for (int ki = 0; ki < coefs_size; ++ki) {
	if (rat) {
	    all_coefs.push_back(lin_sf_u->coefs_begin()[ki*dim]);
	    all_coefs.push_back(lin_sf_v->coefs_begin()[ki*dim]);
	    all_coefs.push_back(sf_1d_cp->rcoefs_begin()[ki*rdim]);
	    all_coefs.push_back(sf_1d_cp->rcoefs_begin()[ki*rdim+1]);
	} else {
	    all_coefs.push_back(lin_sf_u->coefs_begin()[ki*dim]);
	    all_coefs.push_back(lin_sf_v->coefs_begin()[ki*dim]);
	    all_coefs.push_back(sf_1d_cp->coefs_begin()[ki*dim]);
	}
    }
    shared_ptr<SplineSurface> return_sf(new SplineSurface(sf_1d_cp->numCoefs_u(), sf_1d_cp->numCoefs_v(),
							  sf_1d_cp->order_u(), sf_1d_cp->order_v(),
							  sf_1d_cp->basis_u().begin(),
							  sf_1d_cp->basis_v().begin(),
							  all_coefs.begin(), 3, rat));

    return return_sf;
}


//===========================================================================
vector<shared_ptr<SplineSurface> >
SurfaceCreators::separateRationalParts(const SplineSurface& sf)
//===========================================================================
{
    bool rat = sf.rational();
    ASSERT(rat);

    int dim= sf.dimension();
    int rdim = dim + 1;
    vector<shared_ptr<SplineSurface> > sep_sfs;
    vector<double> coefs(sf.coefs_begin(), sf.coefs_end());
    int nmb1 = sf.numCoefs_u();
    int nmb2 = sf.numCoefs_v();
    vector<double> rcoefs;
    int num_coefs = nmb1*nmb2;
    vector<double>::const_iterator rcoef_iter = sf.rcoefs_begin();
    for (int ki = 0; ki < num_coefs; ++ki) {
	rcoefs.push_back(rcoef_iter[ki*rdim+1]);
	for (int kj = 0; kj < dim; ++kj) {
	    coefs[ki*dim+kj] /= (rcoefs.back());
	}
    }
    sep_sfs.push_back(shared_ptr<SplineSurface>
		      (new SplineSurface(nmb1, nmb2, sf.order_u(), sf.order_v(),
					 sf.basis_u().begin(), sf.basis_v().begin(),
					 coefs.begin(), dim)));
    sep_sfs.push_back(shared_ptr<SplineSurface>
		      (new SplineSurface(nmb1, nmb2, sf.order_u(), sf.order_v(),
					 sf.basis_u().begin(), sf.basis_v().begin(),
					 rcoefs.begin(), 1)));

    return sep_sfs;
}


namespace {

//===========================================================================
void transformToLoops(vector<shared_ptr<SplineCurve> >& param_bd_cvs,
		      vector<shared_ptr<SplineCurve> >& space_bd_cvs,
		      const vector<shared_ptr<const ParamSurface> >& sfs)
//===========================================================================
{
    int nmb_bd_cvs = (int)space_bd_cvs.size();
    if (nmb_bd_cvs < 2)
	return;

    // We run through the curves, adding suitable segments.
    for (size_t ki = 0; ki < param_bd_cvs.size(); ++ki) {
// 	ALWAYS_ERROR_IF(self_int_params[ki].size() != 0,
// 		    "Self intersecting curves currently not supported.", InputError());
	RectDomain rd = sfs[ki]->containingDomain();
	double min_box_side = min(rd.umax() - rd.umin(), rd.vmax() - rd.vmin());
	vector<Point> from_pt = param_bd_cvs[ki]->ParamCurve::point(param_bd_cvs[ki]->startparam(), 1);
	vector<Point> to_pt = param_bd_cvs[ki]->ParamCurve::point(param_bd_cvs[ki]->endparam(), 1);
	// If surface is closed, param cv corr to space loop must be treated in a different manner.
	if (from_pt[0].dist(to_pt[0]) > min_box_side/2.0) 
	  {
	    // First find closest point on the boundary loop
	    int ind;
	    double cl_par, cl_dist;
	    Point cl_from, cl_to;
	    const CurveLoop& boundary = sfs[ki]->outerBoundaryLoop();
	    boundary.closestParPoint(from_pt[0], ind, cl_par, cl_from, cl_dist);
	    ALWAYS_ERROR_IF(ind<0, "Missing parameter loop");
	    boundary.closestParPoint(to_pt[0], ind, cl_par, cl_to, cl_dist);
	    ALWAYS_ERROR_IF(ind<0, "Missing parameter loop");

	    from_pt[0] = cl_from;
	    to_pt[0] = cl_to;
	    
	    // Then change the almost identical parameter values (average).
	    if (fabs(from_pt[0][0] - to_pt[0][0]) < (fabs(from_pt[0][1] - to_pt[0][1]))) {
		from_pt[0][0] = to_pt[0][0] = 0.5*(from_pt[0][0] + to_pt[0][0]);
		//param_bd_cvs[ki]->coefs_begin()[0] = param_bd_cvs[ki]->coefs_end()[-2] = from_pt[0][0];
	    } else {
		from_pt[0][0] = to_pt[0][0] = 0.5*(from_pt[0][0] + to_pt[0][0]);
		//param_bd_cvs[ki]->coefs_begin()[1] = param_bd_cvs[ki]->coefs_end()[-1] = from_pt[0][1];
	    }		
	    Point new_pt_start = sfs[ki]->ParamSurface::point(from_pt[0][0], from_pt[0][1]);
	    Point new_pt_end =  sfs[ki]->ParamSurface::point(to_pt[0][0], to_pt[0][1]);
	    Point mid_pt = 0.5*(new_pt_start + new_pt_end);
	    copy(mid_pt.begin(), mid_pt.end(), space_bd_cvs[ki]->coefs_begin());
	    copy(mid_pt.begin(), mid_pt.end(), space_bd_cvs[ki]->coefs_end() - 3);
	    copy(from_pt[0].begin(), from_pt[0].end(), param_bd_cvs[ki]->coefs_begin());
	    copy(to_pt[0].begin(), to_pt[0].end(), param_bd_cvs[ki]->coefs_end() - 2);

	} else {
	    double dl = (rd.upperRight() - rd.lowerLeft()).length(); // Diagonal length.
	    // We make sure that the lines start and end outside sf.
	    Point from1 = from_pt[0] - (dl + 1.0)*from_pt[1]/(from_pt[1].length());
	    Point to1 = from_pt[0] + (dl + 1.0)*from_pt[1]/(from_pt[1].length());
	    Point from2 = to_pt[0] - (dl + 1.0)*to_pt[1]/(to_pt[1].length());
	    Point to2 = to_pt[0] + (dl + 1.0)*to_pt[1]/(to_pt[1].length());
	    // We solve: s*from1 + (1 - s)*to1 = t*from2 + (1 - t)*to2
	    //Matrix A(2, 2);
	    vector<vector<double> > A(2);
	    A[0].resize(2); A[1].resize(2);
	    //ColumnVector b(2);
	    vector<double> b(2);
	    for (int kj = 0; kj < 2; ++kj) {
// 		A.element(kj, 0) = from1[kj] - to1[kj];
// 		A.element(kj, 1) = to2[kj] - from2[kj];
// 		b.element(kj) = to2[kj] - to1[kj];
		A[kj][0] = from1[kj] - to1[kj];
		A[kj][1] = to2[kj] - from2[kj];
		b[kj] = to2[kj] - to1[kj];
	    }
// 	    ColumnVector x(2);
// 	    x = A.i()*b;
	    LUsolveSystem(A, (int)A.size(), &b);

	    // We then calculate param value and corr space value.
	    //Point param_int_pt = x.element(0)*from1 + (1 - x.element(0))*to1;
	    Point param_int_pt = b[0]*from1 + (1 - b[0]) * to1;

	    Point space_int_pt = sfs[ki]->ParamSurface::point(param_int_pt[0], param_int_pt[1]);
	    // Finally we change values of end ctrl pts.
	    copy(param_int_pt.begin(), param_int_pt.end(), param_bd_cvs[ki]->coefs_begin());
	    copy(param_int_pt.begin(), param_int_pt.end(), param_bd_cvs[ki]->coefs_end() - 2);
	    copy(space_int_pt.begin(), space_int_pt.end(), space_bd_cvs[ki]->coefs_begin());
	    copy(space_int_pt.begin(), space_int_pt.end(), space_bd_cvs[ki]->coefs_end() - 3);
	}
    }
}

}; // end anonymous namespace
