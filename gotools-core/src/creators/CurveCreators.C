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

#include "GoTools/creators/CurveCreators.h"

#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/creators/ApproxCurve.h"
#include "GoTools/creators/ProjectCurve.h"
#include "GoTools/creators/HermiteAppC.h"
#include "GoTools/creators/HermiteAppS.h"
#include "GoTools/creators/LiftCurve.h"
#include "GoTools/creators/TrimCurve.h"
#include "GoTools/creators/AdaptCurve.h"
#include "GoTools/creators/EvalParamCurve.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/creators/OffsetCurveNormalDir.h"

#include <vector>
#include <algorithm>
#include <fstream>
#include <iterator>

#ifdef __BORLANDC__
#include <iterator>
#endif

// #ifndef NDEBUG
// #ifndef DEBUG
//#define DEBUG
// #endif // DEBUG
// #endif // NDEBUG

#include "GoTools/geometry/Plane.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;

// Returning parameters to be used by approximation. For a spline curve the unique knots are returned
// (possibly removing knot(s) which are too close to an adjacent knot). For a general curve the end
// parameters are returned.
vector<double> getInitParams(ParamCurve* cv);

//===========================================================================
SplineCurve*
CurveCreators::multCurveWithFunction(const SplineCurve& alpha,
				       const SplineCurve& f)
//===========================================================================
{
    ALWAYS_ERROR_IF(alpha.dimension() != 1,
		    "Dimension of function is different from 1.");

    // Testing parameter interval.
    double start_par = alpha.startparam();
    double end_par = alpha.endparam();
    ALWAYS_ERROR_IF(f.startparam() != start_par || f.endparam() != end_par,
		    "Parameter intervals of the curves do not coincide.");


    // We start by constructing the needed knot vector.
    vector<double> new_alpha_knots, new_f_knots, knots;
    int i, j;
    std::vector<double>::const_iterator iter;
    for (iter = alpha.basis().begin(); iter < alpha.basis().end(); ++iter) {
	new_alpha_knots.push_back(iter[0]);
	// We must insert multiple knots for cont. reasons.
	if ((iter[0] < iter[1]) || (alpha.basis().end() - iter == 1))
	    for (i = 0; i < f.order() - 1; ++i)
		new_alpha_knots.push_back(iter[0]);
    }
    for (iter = f.basis().begin(); iter < f.basis().end(); ++iter) {
	new_f_knots.push_back(iter[0]);
	// We must insert multiple knots for cont. reasons.
	if ((iter[0] < iter[1]) || (f.basis().end() - iter == 1))
	    for (i = 0; i < alpha.order() - 1; ++i)
		new_f_knots.push_back(iter[0]);
    }
    std::set_union(new_alpha_knots.begin(), new_alpha_knots.end(),
		   new_f_knots.begin(), new_f_knots.end(),
		   std::back_inserter(knots));

    int order = alpha.order() + f.order() - 1;
    int num_coefs = (int)knots.size() - order;
    vector<double> coefs_par; // Parameter values for the coefs (Greville).
    vector<double> coefs;

    for (i = 0; i < num_coefs; ++i) {
	double tpar = 0.0;
	for (j = 1; j < order; ++j)
	    tpar += knots[i + j];
	tpar /= (order - 1);
	coefs_par.push_back(tpar);

	// We now must evaluate the product of the curves in tpar.
	Point alpha_pt, f_pt;
	alpha.point(alpha_pt, tpar);
	f.point(f_pt, tpar);
	for (j = 0; j < f.dimension(); ++j)
	    coefs.push_back(alpha_pt[0] * f_pt[j]);
    }

    // We interpolate the points. As we use up all our degrees of freedom
    // we're assured to end up with the wanted spline product.
    vector<double> new_coefs;
    vector<double> dummy_tangents;
    vector<int> dummy_index;
    BsplineBasis basis(num_coefs, order, knots.begin());
    SplineInterpolator interpolator;
    interpolator.setBasis(basis);
    interpolator.interpolate(coefs_par, coefs, dummy_index,
			     dummy_tangents, new_coefs);

    //    num_coefs = (new_coefs.size()/f.dimension());
    SplineCurve* return_curve =
	new SplineCurve(num_coefs, order,
			  interpolator.basis().begin(),
			  new_coefs.begin(), f.dimension());

    return return_curve;
}


//===========================================================================
SplineCurve* CurveCreators::blend(const SplineCurve& alpha_1,
				    const SplineCurve& f_1,
				    const SplineCurve& alpha_2,
				    const SplineCurve& f_2)
//===========================================================================
{
    // Testing dimensions.
    ALWAYS_ERROR_IF(alpha_1.dimension() != 1 || alpha_2.dimension() != 1,
		    "Dimension of (at least one) blend function differs from 1.");

    ALWAYS_ERROR_IF(f_1.dimension() != f_2.dimension(),
		    "Dimension mismatch between the spline curves to be blended.");

    // Testing parameter interval.
    double start_par = alpha_1.startparam();
    double end_par = alpha_1.endparam();
    ALWAYS_ERROR_IF(alpha_2.startparam() != start_par || alpha_2.endparam() != end_par
		|| f_1.startparam() != start_par || f_1.endparam() != end_par ||
		f_2.startparam() != start_par || f_2.endparam() != end_par,
		    "Parameter intervals of the curves do not coincide.");


    SplineCurve* first_curve = multCurveWithFunction(alpha_1, f_1);
    SplineCurve* second_curve = multCurveWithFunction(alpha_2, f_2);

    // Make sure the two curves have the same order.
    int i;
    int order_1 = first_curve->order();
    int order_2 = second_curve->order();
    for (i = 0; i < abs(order_1 - order_2); ++i)
	if (order_1 < order_2) first_curve->raiseOrder();
	else second_curve->raiseOrder();

    // We next put the two curves on the same knot vector.
    vector<double> all_knots, new_knots_first, new_knots_second;
    std::set_union(first_curve->basis().begin(), first_curve->basis().end(),
		   second_curve->basis().begin(), second_curve->basis().end(),
		   std::back_inserter(all_knots));
    std::set_difference(all_knots.begin(), all_knots.end(),
			first_curve->basis().begin(), first_curve->basis().end(),
			std::back_inserter(new_knots_first));
    std::set_difference(all_knots.begin(), all_knots.end(),
			second_curve->basis().begin(), second_curve->basis().end(),
			std::back_inserter(new_knots_second));
    first_curve->insertKnot(new_knots_first);
    second_curve->insertKnot(new_knots_second);

    // We add the spline coefficients of second_curve to first_curve.
    for (i = 0; i < first_curve->numCoefs() * first_curve->dimension(); ++i)
	first_curve->coefs_begin()[i] += second_curve->coefs_begin()[i];

    delete second_curve;

    return first_curve;
}


//===========================================================================
SplineCurve* CurveCreators::approxCurves(shared_ptr<ParamCurve>* first_crv,
					 shared_ptr<ParamCurve>* last_crv,
					 const vector<Point>& start_pt,
					 const vector<Point>& end_pt, 
					 double approxtol,
					 double& maxdist, 
					 int max_iter, int degree)
//---------------------------------------------------------------------------
//
// Purpose: Replace the current boundary pieces by an approximation.
//          Continuity towards neighbours should be C1.
//
//===========================================================================
{
    vector<Point> start_pt_cpy = start_pt;
    vector<Point> end_pt_cpy = end_pt;
    int nmbcrvs = (int)(last_crv - first_crv);

  int dim = -1;
  double len;
  Point pt1, pt2;
  double t1, t2, tint, tpar;
  int nmbsample;
  vector<double> points;
  vector<double> params;
  for (int ki=0; ki<nmbcrvs; ki++)
    {
      dim = first_crv[ki]->dimension();
      // First make an estimate of the length of the curve piece to replace.
      len = first_crv[ki]->estimatedCurveLength();

      // Compute number of sample points.
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
     nmbsample = max(5, min(1000, (int)(len/approxtol)));
#else
    nmbsample = std::max(5, std::min(1000, (int)(len/approxtol)));
#endif

      t1 = first_crv[ki]->startparam(); 
      t2 = first_crv[ki]->endparam();

      // Evaluate the curve in the sample points and make a centripetal
      // length parameterization of the points. Remember the index
      // of the last point.
      if (ki == 0)
	{
	  first_crv[ki]->point(pt1, t1);
	  points.insert(points.end(), pt1.begin(), pt1.end());
	  params.push_back(0.0);
	}
	  
      tint = (t2 - t1)/(double)(nmbsample-1);
      int kj;
      for (kj=1, tpar=t1+tint; kj<nmbsample; kj++, tpar+=tint)
	{
            if (kj == nmbsample - 1)
                tpar = t2; // We want exact match in the point.
	  first_crv[ki]->point(pt2, tpar);
	  points.insert(points.end(), pt2.begin(), pt2.end());
	  params.push_back(params[params.size()-1] + sqrt(pt1.dist(pt2)));
	  //params.push_back(params[params.size()-1] + pt1.dist(pt2));

	  pt1 = pt2;
	}
    }

  // If start_pt or end_pt is missing the point, method will approximate while
  // keeping start/end point of input curve fixed.
  // As approximation method uses cord length parametrization, we make sure input
  // tangents are of length 1.
  int nmb_derivatives = 0;
  if (start_pt_cpy.size() > 1) {
    start_pt_cpy[1].normalize();
    ++nmb_derivatives;
  }
  if (end_pt_cpy.size() > 1) {
    end_pt_cpy[1].normalize();
    ++nmb_derivatives;
  }

#ifdef DEBUG
  std::ofstream of("point_sequence.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  int nmbpt = points.size()/dim;
  of << nmbpt << std::endl;
  for (int ki=0; ki<nmbpt; ++ki)
    {
      for (int kj=0; kj<dim; ++kj)
	of << points[ki*dim+kj] << "  ";
      of << std::endl;
    }
#endif
  
  // Create a curve approximating the points.
  // If max_iter is too large, we risk ending up with spline curve with dense inner knot spacing.
  double avdist;
  int order = degree + 1;
  ApproxCurve approx_curve(points, params, dim, approxtol,
			     order + nmb_derivatives, order);
  approx_curve.setEndPoints(start_pt_cpy, end_pt_cpy);
  approx_curve.setSmooth(1.0e-6);
  //approx_curve.unsetSmooth();
  shared_ptr<SplineCurve> crv = approx_curve.getApproxCurve(maxdist, avdist,
							      max_iter);

  if (maxdist > approxtol) {
    // Either the number of iterations was too small or we didn't make any progress.
      MESSAGE("Failed iterating towards curve pieces! User must decide if close enough.");
  }

  return dynamic_cast<SplineCurve*>(crv->clone());
}


//===========================================================================
void
CurveCreators::projectCurve(shared_ptr<ParamCurve>& space_cv,
			    shared_ptr<ParamSurface>& surf,
			    double epsge,
			    shared_ptr<SplineCurve>& proj_cv,
			    shared_ptr<SplineCurve>& par_cv)
//===========================================================================
{
    ASSERT(space_cv->dimension() == 3);


    // Make surface curve
    shared_ptr<CurveOnSurface> sf_cv = 
      shared_ptr<CurveOnSurface>(new CurveOnSurface(surf, space_cv, false));

    // We must first construct a EvalCurveSet for use in GoHermitAppS.
    shared_ptr<TrimCurve> proj_crv(new TrimCurve(sf_cv.get()));

    // Approximate
    vector<double> initpars;
    if (space_cv->instanceType() == Class_SplineCurve) {
	shared_ptr<SplineCurve> spline_cv =
	    dynamic_pointer_cast<SplineCurve>(space_cv);
	int order = spline_cv->order();
	int nb_coef = spline_cv->numCoefs();
	vector<double>::const_iterator knots = spline_cv->basis().begin();
	initpars.push_back(knots[order-1]);
	for (int kj = order; kj <= nb_coef; ++kj)
	    if (knots[kj] > initpars[initpars.size()-1])
		initpars.push_back(knots[kj]);
    } else {
	initpars.push_back(space_cv->startparam());
	initpars.push_back(space_cv->endparam());
    }
    double tol2 = epsge;
    vector<int> dims(2);
    dims[0] = 3;
    dims[1] = 2;
    HermiteAppS approximator(proj_crv.get(),
			     &initpars[0], (int)initpars.size(),
			     epsge, tol2, dims);
    approximator.refineApproximation();
    vector<shared_ptr<SplineCurve> > proj_crvs = approximator.getCurves();

    proj_cv = proj_crvs[0];
    par_cv = proj_crvs[1];
}


//===========================================================================
vector<shared_ptr<SplineCurve> > 
CurveCreators::curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
                           const BsplineBasis& init_basis, double tol,
                           int nmb_init_sample_pr_seg)
//===========================================================================
{
    vector<shared_ptr<SplineCurve> > result;

    // Get basis information
    int order = init_basis.order();
    int nmb_coef = init_basis.numCoefs();
    vector<double> knots(order+nmb_coef);
    vector<double>::const_iterator start = init_basis.begin();
    vector<double>::const_iterator end = init_basis.end();
    int ki, kj;
    for (ki=0; start<end; ++start, ++ki)
        knots[ki] = *start;

    // Ensure at least cubic degree
    if (order < 4)
    {
	int nn = 4 - order;
	vector<double> knotval;
	init_basis.knotsSimple(knotval);
	vector<double> newknots(knotval.size()*nn);
	size_t kr;
	for (kr=0, kj=0; kr<knotval.size(); ++kr)
            for (ki=0; ki<nn; ++ki)
                newknots[kj++] = knotval[kr];

	vector<double> knots2(knots.size() + newknots.size());;
	std::merge(knots.begin(), knots.end(),
		   newknots.begin(), newknots.end(),
		   knots2.begin());
	order = 4;
	knots = knots2;
	nmb_coef = (int)knots.size() - order;
    }

    double ta = knots[order-1];
    double tb = knots[nmb_coef];

    for (ki=0; ki<nmb_cvs; ++ki)
    {
	shared_ptr<EvalParamCurve> eval_crv(new EvalParamCurve(cvs[ki]));

	// Adapt knot vector to the parameter interval of the curve
	double tc = eval_crv->start();
	double td = eval_crv->end();
	double t1 = knots[order-1];
	double t2 = knots[nmb_coef];
	for (size_t kr=0; kr<knots.size(); ++kr)
            knots[kr] = tc + (knots[kr]-t1)*(td-tc)/(t2-t1);

	// Approximate
	AdaptCurve adapt(eval_crv.get(), tol, nmb_coef, order, knots);
	adapt.unsetSmooth();
	adapt.setNmbInitSamples(nmb_init_sample_pr_seg);
 	//adapt.setSmooth(smoothw);
	// int performed = adapt.approximate();
	int max_iter = 3;
	adapt.approximate(max_iter);
	double maxdist, avdist;
	shared_ptr<SplineCurve> res = adapt.getAdaptCurve(maxdist, avdist);

	res->setParameterInterval(ta, tb);
	result.push_back(res);

#ifdef DEBUG
	std::cout << "Curve accuracy(1): " << maxdist << " " << avdist;
	std::cout << ", nmb coefs = " << res->numCoefs();
	std::cout << ", initial = " << nmb_coef << std::endl;
#endif

	// We could update the initial spline space here to make it fit for
	// the next aproximation, but since AdaptCurve performes refinement
	// in the mid parameter of knot intervals, this is not crucial
    }
    return result;
}


//===========================================================================
vector<shared_ptr<SplineCurve> > 
CurveCreators::curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
                           double tol, int degree,
                           int nmb_init_sample_pr_seg)
//===========================================================================
{
    vector<shared_ptr<SplineCurve> > result;

    int order = degree + 1;
    int nmb_coef = order;  // Initially the spline space is cubic Bezier

    // Approximate first curve
    shared_ptr<EvalParamCurve> eval_crv(new EvalParamCurve(cvs[0]));
    static double fac_tol = 10.0;
    AdaptCurve adapt0(eval_crv.get(), fac_tol*tol, nmb_coef, order);
    adapt0.unsetSmooth();
    adapt0.setNmbInitSamples(nmb_init_sample_pr_seg);
    static int maxiter = 4;
    int performed = adapt0.approximate(maxiter);
    double maxdist, avdist;
    shared_ptr<SplineCurve> res = adapt0.getAdaptCurve(maxdist, avdist);
#ifdef DEBUG_ADAPT
    std::cout << "Curve accuracy(2): " << maxdist << " " << avdist;
    std::cout << ", nmb coefs = " << res->numCoefs() << std::endl;
#endif
    
    // Fetch knot vector
    vector<double> knots;
    knots.insert(knots.end(), res->knotsBegin(), res->knotsEnd());
    order = res->order();
    nmb_coef = res->numCoefs(); // Update info
    double ta = knots[order-1];
    double tb = knots[nmb_coef];
    
    //result.push_back(res);

    // The remaining curves
    // double smoothw = 0.01;
    //for (int ki=1; ki<nmb_cvs; ++ki)
    for (int ki=0; ki<nmb_cvs; ++ki)
    {
	eval_crv = shared_ptr<EvalParamCurve>(new EvalParamCurve(cvs[ki]));

	// Adapt knot vector to the parameter interval of the curve
	double tc = eval_crv->start();
	double td = eval_crv->end();
	double t1 = knots[order-1];
	double t2 = knots[nmb_coef];
	for (size_t kr=0; kr<knots.size(); ++kr)
            knots[kr] = tc + (knots[kr]-t1)*(td-tc)/(t2-t1);

	// Approximate
	AdaptCurve adapt(eval_crv.get(), tol, nmb_coef, order, knots);
	//adapt.setSmooth(smoothw);
	adapt.unsetSmooth();
	adapt.setNmbInitSamples(nmb_init_sample_pr_seg);
 	performed = adapt.approximate();
	res = adapt.getAdaptCurve(maxdist, avdist);

	res->setParameterInterval(ta, tb);
	result.push_back(res);

#ifdef DEBUG_ADAPT
	std::cout << "Curve accuracy(2): " << maxdist << " " << avdist;
	std::cout << ", nmb coefs = " << res->numCoefs();
	std::cout << ", initial = " << nmb_coef << std::endl;
#endif

	// We could update the initial spline space here to make it fit for
	// the next aproximation, but since AdaptCurve performes refinement
	// in the mid parameter of knot intervals, this is not crucial
    }
    return result;
}


//===========================================================================
SplineCurve*
CurveCreators::projectSpaceCurve(shared_ptr<ParamCurve>& space_cv,
				 shared_ptr<ParamSurface>& surf,
				 shared_ptr<Point>& start_par_pt,
				 shared_ptr<Point>& end_par_pt,
				 double epsge,
				 const RectDomain* domain_of_interest)
//===========================================================================
{
    ASSERT(space_cv->dimension() == 3);

    shared_ptr<ParamSurface> surf2 = shared_ptr<ParamSurface>(surf->clone());

    if (surf2->instanceType() == Class_SplineSurface)
    {
        shared_ptr<SplineSurface> tmp_srf = 
            dynamic_pointer_cast<SplineSurface,ParamSurface>(surf2);
	// We make sure the surf is k-regular.
	surf2 = shared_ptr<ParamSurface>
            (tmp_srf->subSurface(tmp_srf->startparam_u(), tmp_srf->startparam_v(),
                                 tmp_srf->endparam_u(), tmp_srf->endparam_v()));
    }

    // If the input curve is degenerate there is no need to go through the whole projection process.
    double deg_eps = std::min(epsge, 1.0e-6);;
    double space_cv_length = space_cv->estimatedCurveLength();
    if (space_cv_length < deg_eps)
    {
        if (start_par_pt && end_par_pt)
        {
            double tmin = space_cv->startparam();
            double tmax = space_cv->endparam();
            SplineCurve* deg_line = new SplineCurve(*start_par_pt, tmin, *end_par_pt, tmax);
            return deg_line;
        }
        else
        {
            return nullptr;
        }
    }

    // We must first construct a EvalCurve for use in GoHermitAppC.
//     shared_ptr<ParamCurve> cv = space_cv;
    shared_ptr<ProjectCurve> proj_crv(new ProjectCurve(space_cv, surf2,
						       start_par_pt, end_par_pt,
						       epsge,
// 						       epsge,
						       domain_of_interest));

#if 0//ndef NDEBUG
    // Check the evaluator base curve in the endpoints. If they don't satisfy the
    // accuracy requirement, then the approximation will never succeed
    Point proj1, proj2, space1, space2, sf_pt1, sf_pt2;
    proj1 = proj_crv->eval(proj_crv->start());
    proj2 = proj_crv->eval(proj_crv->end());
    space_cv->point(space1, space_cv->startparam());
    space_cv->point(space2, space_cv->endparam());
    surf2->point(sf_pt1, proj1[0], proj1[1]);
    surf2->point(sf_pt2, proj2[0], proj2[1]);
    const double dist1 = sf_pt1.dist(space1);
    const double dist2 = sf_pt2.dist(space2);
    // Even though end distance is larger than tolerance, the routine
    // actually compares distance from projection to approximation
    // ... Which makes sense. Hence no need to demand curve to be
    // within tolerance in end points.
    const double max_dist = std::max(dist1, dist2);
    if (max_dist > epsge)
    {
      MESSAGE("Distance in space curve end points: max_dist = "
	      << std::max(dist1, dist2) << ", epsge = " << epsge);
    }
#endif
    
#ifdef DEBUG
      std::ofstream out("project.g2");
      surf->writeStandardHeader(out);
      surf->write(out);
      space_cv->writeStandardHeader(out);
      space_cv->write(out);
#endif

    // Approximate
    vector<double> initpars = getInitParams(space_cv.get());
    double kink_tol = 1e-02; // We only want it to be reasonably smooth.
    HermiteAppC approximator(proj_crv.get(),
			     &initpars[0], (int)initpars.size(),
			       epsge, kink_tol);
    approximator.refineApproximation();
    shared_ptr<SplineCurve> projection_cv = approximator.getCurve();

#if _MSC_VER > 0 && _MSC_VER < 1300
    return dynamic_cast<SplineCurve*>(projection_cv->clone());
#else
    return (projection_cv.get() != NULL) ? projection_cv->clone() : NULL;
#endif
}


//===========================================================================
SplineCurve*
CurveCreators::liftParameterCurve(shared_ptr<ParamCurve>& parameter_cv,
				  shared_ptr<ParamSurface>& surf,
				  double epsge)
//===========================================================================
{
    ASSERT(parameter_cv->dimension() == 2);

    // We must first construct a EvalCurve for use in GoHermitAppC.
    shared_ptr<LiftCurve> lift_crv(new LiftCurve(parameter_cv, surf, epsge));

    // Approximate
    vector<double> initpars = getInitParams(parameter_cv.get());
    HermiteAppC approximator(lift_crv.get(),
			     &initpars[0], (int)initpars.size(),
			       epsge, epsge); // Using input epsge for both geom and kink tol.
    approximator.refineApproximation();
    shared_ptr<SplineCurve> lifted_cv = approximator.getCurve();

    if (!lifted_cv)
    {
        return nullptr;
    }

#if _MSC_VER > 0 && _MSC_VER < 1300
    return dynamic_cast<SplineCurve*>(lifted_cv->clone());
#else
    return lifted_cv->clone();
#endif
}


//===========================================================================
shared_ptr<SplineCurve>
CurveCreators::offsetCurveNormalDir(shared_ptr<ParamCurve>& parameter_cv,
                                    shared_ptr<ParamCurve>& space_cv,
                                    shared_ptr<ParamSurface>& surf,
                                    double epsge,
                                    double offset_dist)
//===========================================================================
{
    ASSERT((parameter_cv.get() == nullptr) || (parameter_cv->dimension() == 2));

    // We must first construct a EvalCurve for use in GoHermitAppC.
    shared_ptr<OffsetCurveNormalDir> offset_crv
        (new OffsetCurveNormalDir(parameter_cv, space_cv, surf, epsge, offset_dist));

    // Approximate
    shared_ptr<ParamCurve> cv = (parameter_cv.get() != nullptr) ? parameter_cv : space_cv;
    vector<double> initpars = getInitParams(cv.get());
    HermiteAppC approximator(offset_crv.get(),
			     &initpars[0], (int)initpars.size(),
                             epsge, epsge); // Using iput epsge for both geom and kink tol.
    approximator.refineApproximation();
    shared_ptr<SplineCurve> offset_cv = approximator.getCurve();

    return offset_cv;
}

#if 0
//===========================================================================
SplineCurve*
CurveCreators::offsetParameterCurve(shared_ptr<ParamCurve>& parameter_cv,
                                    shared_ptr<ParamSurface>& surf,
                                    double offset_dist, double epsge)
//===========================================================================
{
    ASSERT(parameter_cv->dimension() == 2);

    // We must first construct a EvalCurve for use in GoHermitAppC.
    shared_ptr<OffsetCurve> offset_cv(new OffsetCurve(parameter_cv, surf, offset_dist, epsge));

    // Approximate
    vector<double> initpars = getInitParams(parameter_cv.get());
    HermiteAppC approximator(offset_cv.get(),
			     &initpars[0], (int)initpars.size(),
                             epsge, epsge); // Using iput epsge for both geom and kink tol.
    approximator.refineApproximation();
    shared_ptr<SplineCurve> appr_cv = approximator.getCurve();

#if _MSC_VER > 0 && _MSC_VER < 1300
    return dynamic_cast<SplineCurve*>(appr_cv->clone());
#else
    return appr_cv->clone();
#endif
}
#endif


//===========================================================================
SplineCurve* CurveCreators::createCircle(Point center, 
					   Point axis, 
					   Point normal, 
					   double radius)
//===========================================================================
{
    ALWAYS_ERROR_IF(center.dimension() != 3,
		"Input pts must be of dimension 3!"); // We may accept 2 also, but not now...

    // Just to be on the safe side...
    axis.normalize();
    axis *= radius;
    normal.normalize();
    int ki;
    int dim = center.dimension();
    int num_coefs_circ = 9;
    std::vector<double> circle_coefs((dim + 1)*num_coefs_circ); // Rational cv.
    // We're starting in the yz-plane. Circle is defined as unit circle in input plane.
    double weight = 1/sqrt(2.0);
    Point cross = normal % axis;
    for (ki = 0; ki < 3; ++ki) {
	circle_coefs[ki] = center[ki] + axis[ki];
	circle_coefs[ki+4] = weight*(center[ki] + axis[ki] + cross[ki]);
	circle_coefs[ki+8] = center[ki] + cross[ki];
	circle_coefs[ki+12] = weight*(center[ki] - axis[ki] + cross[ki]);
	circle_coefs[ki+16] = center[ki] - axis[ki];
	circle_coefs[ki+20] = weight*(center[ki] - axis[ki] - cross[ki]);
	circle_coefs[ki+24] = center[ki] - cross[ki];
	circle_coefs[ki+28] = weight*(center[ki] + axis[ki] - cross[ki]);
	circle_coefs[ki+32] = center[ki] + axis[ki];
    }
    for (ki = 0; ki < 9; ++ki) {
	circle_coefs[ki*4+3] = (ki % 2 == 0) ? 1.0 : weight;
    }

    // We try to make a rational spline curve from coefs.
    int order_circ = 3; // Quadratic.
    std::vector<double> knots(num_coefs_circ+order_circ);
    int kk;
    for (kk = 0; kk < order_circ; ++kk) {
	knots[kk] = 0.0;
    }
    int val = 1;
    for (; kk < num_coefs_circ; kk+= 2, ++val) {
	knots[kk] = knots[kk+1] = (double)(val);
    }
    for (; kk < num_coefs_circ + order_circ; ++kk) {
	knots[kk] = (double)(val);
    }

    SplineCurve* circle_cv = new SplineCurve(num_coefs_circ, order_circ, knots.begin(),
					     circle_coefs.begin(), dim, true);

    return circle_cv;
}



//===========================================================================
shared_ptr<SplineCurve> CurveCreators::insertParamDomain(const SplineCurve& cv_1d,
							 double knot_tol)
//===========================================================================
{
    if (cv_1d.rational()) {
	MESSAGE("Rational case not supported yet! Using non-rational coefs.");
    }

    // The returned object should be linear in the first direction.
    // We create an additional 1d-cv describing the linear param space.
    vector<double> lin_knots(4, cv_1d.startparam());
    lin_knots[2] = lin_knots[3] = cv_1d.endparam();
    vector<double> lin_coefs(2);
    lin_coefs[0] = lin_knots[0];
    lin_coefs[1] = lin_knots[2];
    shared_ptr<SplineCurve> lin_cv(new SplineCurve(2, 2, lin_knots.begin(), lin_coefs.begin(), 1));

    // We then make sure the domains are equal (i.e. we may need to insert knots into basis of lin_cv).
    vector<shared_ptr<SplineCurve> > cvs;
    cvs.push_back(lin_cv);
    cvs.push_back(shared_ptr<SplineCurve>(cv_1d.clone()));
    GeometryTools::unifyCurveSplineSpace(cvs, knot_tol);

    // We then create our param cv (i.e. living a 2-dimensional domain).
    vector<double> all_coefs;
    int nmb_ctl_pts = cvs[0]->numCoefs();
    for (int ki = 0; ki < nmb_ctl_pts; ++ki) {
	all_coefs.push_back(cvs[0]->coefs_begin()[ki]);
	all_coefs.push_back(cvs[1]->coefs_begin()[ki]);
    }
    shared_ptr<SplineCurve> return_cv(new SplineCurve(cvs[0]->numCoefs(), cvs[0]->order(),
						      cvs[0]->basis().begin(), all_coefs.begin(), 2));

    return return_cv;
}

//===========================================================================
SplineCurve* CurveCreators::offsetCurve(const SplineCurve& base_cv, Point offset_val)
//===========================================================================
{

    SplineCurve* offset_cv(new SplineCurve(base_cv));

    // We then add offset_val to all coefs (should handle rational cvs as well).
    std::vector<double>::iterator coefs_iter = offset_cv->coefs_begin();
    int dim = offset_cv->dimension();
    ASSERT(dim == offset_val.dimension());
    while (coefs_iter < offset_cv->coefs_end())
    {
	for (int ki = 0; ki < dim; ++ki)
	    coefs_iter[ki] += offset_val[ki];

	coefs_iter += dim;
    }

    return offset_cv;
}


//===========================================================================
vector<double> getInitParams(ParamCurve* cv)
//===========================================================================
{
    vector<double> init_pars;

    if (cv->instanceType() == Class_SplineCurve)
    {
        SplineCurve* spline_cv = dynamic_cast<SplineCurve*>(cv);
        int order = spline_cv->order();
        int nb_coef = spline_cv->numCoefs();
        vector<double>::const_iterator knots = spline_cv->basis().begin();
        init_pars.push_back(knots[order-1]);
        for (int kj = order; kj <= nb_coef; ++kj)
            if (knots[kj] > init_pars[init_pars.size()-1])
                init_pars.push_back(knots[kj]);
    }
    else
    {
        init_pars.push_back(cv->startparam());
        init_pars.push_back(cv->endparam());
    }

    // // The allowed knot distance is scaled according to the curve length (i.e. as if the curve was curve
    // // length parametrized).
    // double cv_length = cv->estimatedCurveLength(std::max((int)init_pars.size(), 4));
    // double scale = cv_length/(cv->endparam() - cv->startparam());

    // Check knot suggestions for distribution
    double avdist = (init_pars[init_pars.size()-1] - init_pars[0])/(double)(init_pars.size()-1);
    double mindist = 0.1*avdist;///scale;
    for (int kj=1; kj<(int)init_pars.size(); ++kj)
    {
        if (init_pars.size() <= 2)
            break;
        if (init_pars[kj]-init_pars[kj-1] < mindist)
        {
            // Remove initial knot
            int ix;
            if (kj == 1)
                ix = kj;
            else if (kj == (int)init_pars.size()-1)
                ix = kj-1;
            else
            {
                // If the length of the tangents are highly varying on each side of the parameter we
                // remove the knot with the least variance.
                vector<Point> left1 = cv->point(init_pars[kj-2], 1);
                double left1_length = left1[1].length();
                vector<Point> left2 = cv->point(init_pars[kj-1], 1);
                double left2_length = left2[1].length();
                vector<Point> right1 = cv->point(init_pars[kj], 1);
                double right1_length = right1[1].length();
                vector<Point> right2 = cv->point(init_pars[kj+1], 1);
                double right2_length = right2[1].length();
                // We compute the change of tangent length if we were to remove either knot.
                double frac_left = (left1_length > right1_length) ?
                    (left1_length/right1_length) : (right1_length/left1_length);
                double frac_right = (left2_length > right2_length) ?
                    (left2_length/right2_length) : (right2_length/left2_length);
                double min_frac = min(frac_left, frac_right);
                double max_frac = max(frac_left, frac_right);
                if (max_frac < 2.0) // 201905 Note that this value may need some tuning.
                {
                    ix = (init_pars[kj-1]-init_pars[kj-2] < init_pars[kj+1]-init_pars[kj])
                        ? kj-1 : kj;
                }
                else
                {
                    if (min_frac > 2.0)
                    {
                        // The change of tangents is too large, we keep both knots.
                        continue;
                    }
                    else
                    {
                        ix = (frac_left < frac_right) ? kj -1 : kj;
                    }
                }
            }
            init_pars.erase(init_pars.begin()+ix);
            if (ix == kj-1)
                --kj;
        }
    }

    return init_pars;
}
