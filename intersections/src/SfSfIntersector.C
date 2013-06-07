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

#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/intersections/SfSelfIntersector.h"
#include "GoTools/intersections/SfCvIntersector.h"
#include "GoTools/geometry/extremalPtSurfSurf.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/utils/RotatedBox.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/intersections/IntersectionLink.h"
#include "GoTools/utils/Values.h"
#include "GoTools/geometry/ObjectHeader.h" // for debugging
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/IntersectionUtils.h"
#include "GoTools/intersections/Spline2FunctionInt.h"
#include "GoTools/implicitization/ImplicitizeCurveAndVectorAlgo.h"
#include "GoTools/intersections/Par2FuncIntersector.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/geometry/Utils.h"
#include <limits>
#include <stdio.h> // for debugging
#include <iostream>
#include <fstream> // For debugging


using std::vector;
using std::swap;
using std::cout;
using std::endl;
using std::max;
using std::pair;
using std::make_pair;
using std::not1;
using std::set;
using std::numeric_limits;


namespace { // Anonymous namespace


// predicate for STL function
class CrossesValue {
public:
    CrossesValue(int par_dir, double value, double eps) 
	: par_dir_(par_dir), value_(value), eps_(eps) {}
    bool operator()(const shared_ptr<Go::IntersectionLink>& l) const 
    {
	Go::IntersectionPoint *p1, *p2;
	l->getIntersectionPoints(p1, p2);
	double t1 = p1->getPar()[par_dir_];
	double t2 = p2->getPar()[par_dir_];
	if (t1 > t2) {
	    swap(t1, t2);
	}
	return (value_ - t1) > eps_ && (t2 - value_) > eps_;
    }
    typedef const shared_ptr<Go::IntersectionLink> argument_type;
    typedef bool result_type;
private:
    int par_dir_;
    double value_;
    double eps_;
};


// predicate for STL function
class TestInDomain {
public:
    TestInDomain(const Go::ParamObjectInt* obj1, 
		 const Go::ParamObjectInt* obj2,
		 shared_ptr<Go::IntersectionPoint> ref_point)
    {
	lower_limit_.reserve(4);
	upper_limit_.reserve(4);
	epsilon_.reserve(4);
	int obj1_params = obj1->numParams();
	int obj2_params = obj2->numParams();
	int i;
	for (i = 0; i < obj1_params; ++i) {
	    lower_limit_.push_back(obj1->startParam(i));
	    upper_limit_.push_back(obj1->endParam(i));
	    epsilon_.push_back(ref_point->parameterTolerance(i));
	}
	for (i = 0; i < obj2_params; ++i) {
	    lower_limit_.push_back(obj2->startParam(i));
	    upper_limit_.push_back(obj2->endParam(i));
	    epsilon_.push_back(ref_point->parameterTolerance(obj1_params + i));
	}
    }

    bool operator()(const shared_ptr<Go::IntersectionLink>& l) const
    {
	Go::IntersectionPoint *p1, *p2;
	l->getIntersectionPoints(p1, p2);

	int num_param = (int)lower_limit_.size();
	for (int i = 0; i < num_param; ++i) {
	    if (p1->getPar(i) < lower_limit_[i] - epsilon_[i] ||
		p1->getPar(i) > upper_limit_[i] + epsilon_[i] ||
		p2->getPar(i) < lower_limit_[i] - epsilon_[i] ||
		p2->getPar(i) > upper_limit_[i] + epsilon_[i]) {
		// p1 or p2 for this link was found to lie outside of
		// the domain
		return false;
	    }
	}
	// everything was inside the domain
	return true;
    }
    typedef const shared_ptr<Go::IntersectionLink> argument_type;
    typedef bool result_type;
private:
    vector<double> lower_limit_;
    vector<double> upper_limit_;
    vector<double> epsilon_;
};


} // End anonymous namespace


namespace Go {


//===========================================================================
SfSfIntersector::SfSfIntersector(shared_ptr<ParamGeomInt> obj1, 
				 shared_ptr<ParamGeomInt> obj2,
				 double epsge,
				 Intersector* prev)
    : Intersector2Obj(obj1, obj2, epsge, prev),
      approx_implicit_(2), approx_implicit_err_(2, -1.0),
      approx_implicit_gradsize_(2, -1.0), approx_implicit_gradvar_(2, -1.0)
//===========================================================================
{
    approx_implicit_[0] = vector<shared_ptr<Param2FunctionInt> >(2);
    approx_implicit_[1] = vector<shared_ptr<Param2FunctionInt> >(2);
    prev_implicit_[0] = prev_implicit_[1] = false;
    if (prev) {
	setApproxImplicitFromPrev();
    }
    sorting_obj_[0] = -1;
    sorting_obj_[1] = -1;
}


//===========================================================================
SfSfIntersector::SfSfIntersector(shared_ptr<ParamGeomInt> obj1, 
				 shared_ptr<ParamGeomInt> obj2,
				 shared_ptr<GeoTol> epsge, 
				 Intersector* prev)
    : Intersector2Obj(obj1, obj2, epsge, prev), 
      approx_implicit_(2), approx_implicit_err_(2, -1.0),
      approx_implicit_gradsize_(2, -1.0), approx_implicit_gradvar_(2, -1.0)
//===========================================================================
{
    approx_implicit_[0] = vector<shared_ptr<Param2FunctionInt> >(2);
    approx_implicit_[1] = vector<shared_ptr<Param2FunctionInt> >(2);
    prev_implicit_[0] = prev_implicit_[1] = false;
    if (prev) {
	setApproxImplicitFromPrev();
    }
    sorting_obj_[0] = -1;
    sorting_obj_[1] = -1;
}


//===========================================================================
SfSfIntersector::~SfSfIntersector()
//===========================================================================
{
  // Currently empty
}


//===========================================================================
shared_ptr<Intersector> 
SfSfIntersector::lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
				       shared_ptr<ParamGeomInt> obj2, 
				       Intersector* prev,
				       int eliminated_parameter,
				       double eliminated_value)
//===========================================================================
{
    // Necessarily a SfCvIntersector

    shared_ptr<SfCvIntersector>
	intersector(new SfCvIntersector(obj1, obj2, epsge_, prev,
					eliminated_parameter,
					eliminated_value));

    return intersector;
}


//===========================================================================
int SfSfIntersector::performRotatedBoxTest(double eps1, double eps2)
//===========================================================================
{
    // Purpose: Perform a rotated box test between two surfaces

    double ang_tol = epsge_->getAngleTol();

    ParamSurfaceInt *surf[2];
    surf[0] = obj_int_[0]->getParamSurfaceInt();
    surf[1] = obj_int_[1]->getParamSurfaceInt();
    DEBUG_ERROR_IF(surf[0] == 0 || surf[1] == 0, "Error in data structure");

    //Only done for simple surfaces
    if (!surf[0]->isSimple() ||
	!surf[1]->isSimple())
	return 1;

    // Define coordinate system for the rotation. First check if an
    // intersection point between the surfaces is found
    vector<Point> axis(2), der(3);
    Point norm;
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);
    int nmb_rec = nmbRecursions();
    int sf_ind = (nmb_rec % 2 == 0) ? 0 : 1;
    if (int_pts.size() > 0)
    {
	// Make coordinate system from the partial derivatives of one surface
	// in the (first) intersection point
	const double *param = (sf_ind == 0) ? int_pts[0]->getPar1() : 
	    int_pts[0]->getPar2();
	obj_int_[sf_ind]->point(der, param, 1);
	if (der[1].angle(der[2]) <= ang_tol)
	{
	    param = (sf_ind == 1) ? int_pts[0]->getPar1() : 
		int_pts[0]->getPar2();
	    obj_int_[1-sf_ind]->point(der, param, 1);
	}
	axis[0] = der[1];
	norm = der[1].cross(der[2]);
	axis[1] = norm.cross(der[1]);
    }
    else
    {
	// Make axis from surface corners
	surf[sf_ind]->axisFromCorners(der[1], der[2]);
	axis[0] = der[1];
	norm = der[1].cross(der[2]);
	axis[1] = norm.cross(der[1]);
    }

    // Make rotated boxes
    if (axis[0].length() < epsge_->getNumericalTol() ||
	axis[1].length() < epsge_->getNumericalTol())
	return 1;
    axis[0].normalize();
    axis[1].normalize();
    if (axis[0].angle_smallest(axis[1]) < epsge_->getNumericalTol())
	return 1;

    RotatedBox box1 = surf[0]->getRotatedBox(axis);
    RotatedBox box2 = surf[1]->getRotatedBox(axis);

    bool overlap = box1.overlaps(box2, eps1, eps2);
    return (overlap) ? 1 : 0;
}


//===========================================================================
bool SfSfIntersector::foundIntersectionNearBoundary()
//===========================================================================
{
    //  Purpose: Check for an intersection configuration that
    //  contradicts the conclusion no intersections in the inner

    // Fetch ParamSurface information
    ParamSurfaceInt* surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt* surf2 = obj_int_[1]->getParamSurfaceInt();
    ASSERT(surf1 != 0 && surf2 != 0);

    ParamSurface *srf1 = surf1->getParamSurface().get();
    ParamSurface *srf2 = surf2->getParamSurface().get();

    // Fetch existing intersection points
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    for (size_t ki=0; ki<int_pts.size(); ki++)
    {
	IntPtClassification type = int_pts[ki]->getClassification(srf1, srf2);
	if (type != DIR_TOUCH && type != DIR_PARALLEL)
	    return true;

	// Check if another intersection branch exist
	int nmb_branch = int_pts[ki]->numBranches();
	if (nmb_branch == 2)
	{
	    type = int_pts[ki]->getClassification(srf1, srf2, 1);
	    if (type != DIR_TOUCH && type != DIR_PARALLEL)
		return true;
	}
    }
    return false;
}


//===========================================================================
int SfSfIntersector::performInterceptionByImplicitization()
//===========================================================================
{
    // Purpose: Use implicitization to check if an intersection
    // between the current surfaces is possible.

    // return 0 -> Guaranteed no intersection
    // return 1 -> Maybe intersection

    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1') {
	int_results_->writeDebug();
	cout << "Attempt to do interception by implicitization..." << endl;
    }

    // VSK, 0612. For the time being...
    /*if (obj_int_[0]->coneLargerThanPi() ||
	obj_int_[1]->coneLargerThanPi())
	return 1;  // No point in trying to implicitize */

    // This number represents a reasonable, but not too small, error
    // for the implicitization. NOTE: It corresponds to the threshold
    // for finding the nullspace in the implicitization process.
    double approx_implicit_err = 1.0e-12;

    // Create the implicit objects - if they don't already exist
    createApproxImplicit(approx_implicit_err);

    bool maybe_intersection = true;
    bool first_implicit_exists
	= (approx_implicit_[0][0].get() != 0)
	&& (approx_implicit_[0][1].get() != 0);
    if (first_implicit_exists) {
	CompositeBox box1 = approx_implicit_[0][0]->compositeBox();
	CompositeBox box2 = approx_implicit_[0][1]->compositeBox();
	double overlap;
	box1.getOverlap(box2, overlap);
	double maxerr = max(approx_implicit_err_[0], approx_implicit_err);
	if (overlap < -maxerr) {
// 	if (overlap < -approx_implicit_err_[0]) {
// 	if (overlap < 0) {
	    maybe_intersection = false;
	}
    }
    bool second_implicit_exists
	= (approx_implicit_[1][0].get() != 0)
	&& (approx_implicit_[1][1].get() != 0);
    if (second_implicit_exists) {
	CompositeBox box1 = approx_implicit_[1][0]->compositeBox();
	CompositeBox box2 = approx_implicit_[1][1]->compositeBox();
	double overlap;
	box1.getOverlap(box2, overlap);
	double maxerr = max(approx_implicit_err_[1], approx_implicit_err);
	if (overlap < -maxerr) {
// 	if (overlap < -approx_implicit_err_[1]) {
// 	if (overlap < 0) {
	    maybe_intersection = false;
	}
    }

    if (maybe_intersection) {
	if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1')
	    cout << "...No success." << endl;
	return 1;
    } else {
	if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1')
	    cout << "...Success." << endl;
	return 0;
    }

    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1')
	cout << "...No success." << endl;
    return 1;   //  0; @@@ VSK, Check this up
}


//===========================================================================
int SfSfIntersector::interceptionBySeparationSurface()
//===========================================================================
{
    // Purpose: Use implicitization to create a separation surface in
    // order to check if an intersection between the current surfaces
    // is possible. Mainly of interest for the self-intersection
    // case.

    // return 0 -> Guaranteed no intersection
    // return 1 -> Maybe intersection

    // Check if the surfaces are splines
    shared_ptr<SplineSurfaceInt> spsf_int1 =
	dynamic_pointer_cast<SplineSurfaceInt, ParamObjectInt>(obj_int_[0]);
    shared_ptr<SplineSurfaceInt> spsf_int2 =
	dynamic_pointer_cast<SplineSurfaceInt, ParamObjectInt>(obj_int_[1]);
    if(spsf_int1.get() == 0 || spsf_int2.get() == 0) {
	// Not splines - we can't guarantee anything
	return 1; 
    }

    // Check existence of intersection along a common boundary
    double frac = 0.1;
    vector<BoundaryIntersectionData> bd_ints;
    int_results_->intersectAlongCommonBoundary(frac, bd_ints);
    if (bd_ints.size() == 0) {
	// No intersection on boundary - we can't guarantee anything
	return 1;
    }

    // Make the average centre vector of the two surface normal cones to
    // use for the creation of the separation surface
    DirectionCone cone1 = spsf_int1->directionCone();
    DirectionCone cone2 = spsf_int2->directionCone();
    if (cone1.greaterThanPi() && cone2.greaterThanPi()) {
	// Cones are too big - anything can happen
	return 1;
    }

    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1') {
	int_results_->writeDebug();
	cout << "Attempt to do interception by implicit separation "
	     << "surface..." << endl;
    }
    
    Point mid_vec = (cone1.greaterThanPi()) ? cone2.centre() :
	((cone2.greaterThanPi()) ? cone1.centre() : 0.5*(cone1.centre()
							 + cone2.centre()));

    // What tolerance should we use? Well, 'tol' is related to the
    // implicit representation, so we are using the "universal
    // implicitization tolerance" 1.0e-12. @jbt
    double tol = 1.0e-12;
//     double tol = 0.5*(epsge_->getEpsge());

    // Make separation surface
    // First fetch boundary curve 
    for (size_t ki=0; ki<bd_ints.size(); ki++) {

	// Find the parameter span covered by the intersection
	//int idx1 = 1 - bd_ints[ki].dir[0];
	//int idx2 = 5 - bd_ints[ki].dir[1];  // Surface number 2
	const int idx1 = bd_ints[ki].dir[0];
	const int idx2 = bd_ints[ki].dir[1]; // surface number 2
	
	size_t npt = bd_ints[ki].pts.size();
	double frac1 = fabs(bd_ints[ki].pts[0]->getPar(idx1) -
			    bd_ints[ki].pts[npt-1]->getPar(idx1));
	frac1 /= (obj_int_[0]->endParam(idx1) -
		  obj_int_[0]->startParam(idx1));
	double frac2 = fabs(bd_ints[ki].pts[0]->getPar(idx2) -
			    bd_ints[ki].pts[npt-1]->getPar(idx2));
	frac2 /= (obj_int_[1]->endParam(idx2-2) -
		  obj_int_[1]->startParam(idx2-2));

	// Choose the smallest fraction
	shared_ptr<SplineCurve> bd_crv;
	double par1 = bd_ints[ki].par[0];
	double par2 = bd_ints[ki].par[1];
	if (frac1 < frac2) {
// 	    double tmppar = bd_ints[ki].par[bd_ints[ki].dir[0]];
// 	    bd_crv = shared_ptr<SplineCurve>
// 	      (spsf_int1->constParamCurve(tmppar, bd_ints[ki].dir[0] == 1));
	    bd_crv = shared_ptr<SplineCurve>
		(spsf_int1->constParamCurve(par1, bd_ints[ki].dir[0] == 0));
	} else {
// 	    double tmppar = bd_ints[ki].par[bd_ints[ki].dir[1]];
// 	    bd_crv = shared_ptr<SplineCurve>
// 		(spsf_int2->constParamCurve(tmppar, bd_ints[ki].dir[1] == 3));
	    bd_crv = shared_ptr<SplineCurve>
		(spsf_int2->constParamCurve(par2, bd_ints[ki].dir[1] == 2));
	}

	if (bd_crv->numCoefs() > 2*bd_crv->order())
	    return 1;

	int min_deg = 1, max_deg = 5;
	int deg;  // Implicit degree

	for (deg=min_deg; deg<=max_deg; deg++) {

	    // Make an implicit surface along the given boundary curve
	    // ruled by the average centre vector of the two surface
	    // normal cones.
	    ImplicitizeCurveAndVectorAlgo genSurf(*(bd_crv.get()),
						  mid_vec, deg);
	    //genSurf.setTolerance(tol);
	    genSurf.perform();

	    BernsteinTetrahedralPoly implicit;
	    BaryCoordSystem3D bc;
	    double sigma_min;
	    genSurf.getResultData(implicit, bc, sigma_min);

	    // Check if the two spline surfaces lie on either side of
	    // the implicit separation surface
	    shared_ptr<SplineSurface> spline_sf_1d_1
		= IntersectionUtils::insertSfInImplObj
		(*(spsf_int1->getSplineSurface()), implicit, bc);

	    // Check if the surface coefficients have one sign (using
	    // tolerance). Avoid checking the common boundary curve
	    typedef vector<double>::const_iterator const_iter;
	    const_iter first = spline_sf_1d_1->coefs_begin();
	    const_iter last = spline_sf_1d_1->coefs_end();
	    const_iter it;
	    int sgn = 0;
	    int kn1 = spline_sf_1d_1->numCoefs_u();
	    int kn2 = spline_sf_1d_1->numCoefs_v();
	    int ki, kj;
	    int ki_bd = -1, kj_bd = -1;
	    if (idx1 == 1) {
		ki_bd = (par1-spsf_int1->startParam(0)
			 < epsge_->getRelParRes()) ? 0 : kn1-1;
	    } else {
		kj_bd = (par1-spsf_int1->startParam(1)
			 < epsge_->getRelParRes()) ? 0 : kn2-1;
	    }
	    for (it=first, kj=0; kj<kn2; kj++) {
		for (ki=0; ki<kn1; ki++, it++) {
		    if (ki == ki_bd || kj == kj_bd) {
			continue;
		    }
		    if ((*it) < -tol) {
			if (sgn > 0) {
			    break;
			} else { 
			    sgn = -1;
			}
		    }
		    if ((*it) > tol) {
			if (sgn < 0) {
			    break;
			} else { 
			    sgn = 1;
			}
		    }
		}
		if (ki < kn1) {
		    break;
		}
	    }
	    if (it < last) {
		continue;   // This boundary did not separate
	    }

	    shared_ptr<SplineSurface> spline_sf_1d_2 =
		IntersectionUtils::insertSfInImplObj
		(*(spsf_int2->getSplineSurface()), implicit, bc);

	    first = spline_sf_1d_2->coefs_begin();
	    last = spline_sf_1d_2->coefs_end();
	    kn1 = spline_sf_1d_2->numCoefs_u();
	    kn2 = spline_sf_1d_2->numCoefs_v();
	    sgn *= -1;  // The other spline surface is supposed to lie
			// on the other side of the implicit surface
	    ki_bd = kj_bd = -1;

	    if (idx2 == 3) {
		ki_bd = (par2-spsf_int2->startParam(0) 
			 < epsge_->getRelParRes()) ? 0 : kn1-1;
	    } else {
		kj_bd = (par2-spsf_int2->startParam(1) 
			 < epsge_->getRelParRes()) ? 0 : kn2-1;
	    }
	    for (it=first, kj=0; kj<kn2; kj++) {
		for (ki=0; ki<kn1; ki++, it++) {
		    if (ki == ki_bd || kj == kj_bd) {
			continue;
		    }
		    if ((*it) < -tol) {
			if (sgn > 0) {
			    break;
			} else {
			    sgn = -1;
			}
		    }
		    if ((*it) > tol) {
			if (sgn < 0) {
			    break;
			} else {
			    sgn = 1;
			}
		    }
		}
		if (ki < kn1) {
		    break;
		}
	    }
	    if (it < last) {
		continue;   // This boundary did not separate
	    }
	    break;  // Separation succeeded
	}
	if (deg <= max_deg) {
	    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1')
		cout << "...Success." << endl;
	    return 0;  // Separation succeeded
	}
    }

    // The separation surface does not split the two spline surfaces
    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1')
	cout << "...No success." << endl;
    return 1;
}


//===========================================================================
int SfSfIntersector::simpleCase2(Point& axis1, Point& axis2)
//===========================================================================
{
    //  Purpose: Do extra simple case test between two parametric
    //  surfaces Based on the sisl function s1795.

    // The function only works for spline surface. Check if we have the
    // correct input.
    int kr;
    for (kr=0; kr<2; kr++)
    {
	if (!(obj_int_[kr]->isSpline()))
	    return 0;
    }

    // We are making a cone surrounding the orientating surface
    // on the unit sphere. The cone is representated with centre
    // coordinates and an angle. The orientation is computed
    // from aproximation of the normal to the surface.
    double angle = axis1.angle_smallest(axis2);
    if (angle > 0.5*M_PI)
	angle = M_PI - angle;

    double tangle[2];
    double tlen = axis1*axis2;
    Point scen1, scen2;

    // Compute the cone angles
    for (kr=0, scen1=axis1, scen2=axis2-tlen*axis1; 
	 kr<2; kr++, scen1=axis2, scen2=axis1-tlen*axis2)
    {
	scen2.normalize();
	tangle[kr] = obj_int_[kr]->getOptimizedConeAngle(scen1, scen2); 
    }
    if (tangle[0] + tangle[1] < angle + epsge_->getRelParRes())
	return 1;  // Simple case
    else
	return 0;
}


//===========================================================================
int SfSfIntersector::simpleCaseByImplicitization()
//===========================================================================
{
    //  Purpose: Use implicitization to check if an inner closed
    //  intersection loop is possible for the current surfaces

    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1') {
	int_results_->writeDebug();
	cout << "Attempt to get simple case by implicitization..." << endl;
    }

/*    // VSK, 0612. For the time being...
    if (obj_int_[0]->coneLargerThanPi() ||
	obj_int_[1]->coneLargerThanPi())
	return 0;  // No point in trying to implicitize */

    const double approx_implicit_err = 1.0e-12;
    createApproxImplicit(approx_implicit_err);

    // Check for cylinder-plane intersection. (Actually, we can
    // replace 'cylinder' with something that is linear in one
    // direction.)
    int simple_case = 0; // Not cylinder or plane
    if (approx_implicit_[0][1].get() != 0
	&& approx_implicit_[1][0].get() != 0) {
	shared_ptr<ParamSurface>
	    parfunc1(approx_implicit_[0][1]->getParamSurface());
	shared_ptr<ParamSurface>
	    parfunc2(approx_implicit_[1][0]->getParamSurface());
	shared_ptr<SplineSurface> splinefunc1
	    = dynamic_pointer_cast<SplineSurface, ParamSurface>(parfunc1);
	shared_ptr<SplineSurface> splinefunc2
	    = dynamic_pointer_cast<SplineSurface, ParamSurface>(parfunc2);
	int deg[4];
	deg[0] = splinefunc1->order_u() - 1;
	deg[1] = splinefunc1->order_v() - 1;
	deg[2] = splinefunc2->order_u() - 1;
	deg[3] = splinefunc2->order_v() - 1;
	int ndeg1 = 0;
	int ndeg2 = 0;
	for (int i = 0; i < 4; ++i) {
	    if (deg[i] == 1) {
		++ndeg1;
	    }
	    if (deg[i] >= 2) {
		++ndeg2;
	    }
	}
	// Cylinder-plane happens when one of the degrees equals 1, while
	// the rest equals 2.
	if (ndeg1 == 1 && ndeg2 == 3) {
	    simple_case = 1;
	}
    }

    // Check for monotonicity
    if (simple_case == 0) {
	Point mon_dir;
	double zero_tol = 0.1*epsge_->getEpsge();
	bool is_monotone;
	is_monotone = approx_implicit_[0][1].get() != 0
	    && approx_implicit_err_[0] <= approx_implicit_err
	    && approx_implicit_[0][1]->monotone(mon_dir, zero_tol);
	if (is_monotone) {
	    simple_case = 1; // Monotone
	    sorting_obj_[0] = 1;
	    sorting_dir_[0].setValue(-mon_dir[1],mon_dir[0]);
	    //sorting_dir_.setValue(mon_dir[0],mon_dir[1]);
	    bool success = int_results_->checkSortingDir(sorting_dir_[0],
							 sorting_obj_[0]);
	    sorting_pardir_[0] = success;
	    if (!success) {
		cout << "Inconsistent tangent directions  " << sorting_obj_
		     << endl;
		//simple_case = 0;
		//sorting_obj_ = -1;
	    }
	}
	is_monotone = approx_implicit_[1][0].get() != 0
	    && approx_implicit_err_[1] <= approx_implicit_err
	    && approx_implicit_[1][0]->monotone(mon_dir, zero_tol);
	if (is_monotone) {
	    int idx = (sorting_obj_[0] < 0) ? 0 : 1;
	    if (idx == 1 && sorting_pardir_[0] == false)
	    {
		sorting_obj_[1] = sorting_obj_[0];
		sorting_dir_[1] = sorting_dir_[0];
		sorting_pardir_[1] = sorting_pardir_[0];
		idx = 0;
	    }
	    simple_case = 1; // Monotone
	    sorting_obj_[idx] = 0;
	    sorting_dir_[idx].setValue(-mon_dir[1],mon_dir[0]);
	    //sorting_dir_.setValue(mon_dir[0],mon_dir[1]);
	    bool success = int_results_->checkSortingDir(sorting_dir_[idx],
							 sorting_obj_[idx]);
	    sorting_pardir_[idx] = success;
	    if (!success) {
		cout << "Inconsistent tangent directions  " << sorting_obj_
		     << endl;
		//simple_case = 0;
		//sorting_obj_ = -1;
	    }
	}
    }

    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1') {
	if (simple_case) {
	    cout << "...Success." << endl;
	}
	else
	    cout << "...No success." << endl;
    }
    return simple_case;
}


//===========================================================================
int SfSfIntersector::checkCoincidence()
//===========================================================================
{
    // Intersection pool checks for closed loops of intersection points
    // along the surface boundaries
    vector<vector<shared_ptr<IntersectionPoint> > > loop_ints;
    bool found_loop = int_results_->fetchLoops(loop_ints);
    if (!found_loop)
	return 0;  // No closed loop of intersections along the boundaries

    double min_len = obj_int_[0]->endParam(0)-obj_int_[0]->startParam(0);
    for (int ki=1; ki<obj_int_[0]->numParams(); ki++)
    {
	double td = obj_int_[0]->endParam(ki)-obj_int_[0]->startParam(ki);
	if (td < min_len)
	    min_len = td;
    }
    for (int ki=0; ki<obj_int_[1]->numParams(); ki++)
    {
	double td = obj_int_[1]->endParam(ki)-obj_int_[1]->startParam(ki);
	if (td < min_len)
	    min_len = td;
    }

    ParamSurfaceInt *surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt *surf2 = obj_int_[1]->getParamSurfaceInt();
    DEBUG_ERROR_IF(surf1 == 0 || surf2 == 0, "Error in data structure");
    for (int ki=0; ki<int(loop_ints.size()); ki++)
    {
	// For each boundary loop, check if the surface areas surrounded by
	// the loop represents a PAC area.
	// First fetch if coincidence is already marked
	int is_coincident = int_results_->isPAC(loop_ints[ki]);

	if (!is_coincident)
	{
	    // VSK, 0609. Check if the loop is degenerate or 
	    // (0610) cover only a minor part of the surfaces
	    double ptol = std::max(0.1*min_len, epsge_->getRelParRes());
	   
	    vector<double> par0, par1;
	    par0 = loop_ints[ki][0]->getPar();
	    int k1, k2;
	    for (k1=0; k1<int(loop_ints[ki].size()); k1++)
	    {
		par1 = loop_ints[ki][k1]->getPar();
		for (k2=0; k2<int(par0.size()); k2++)
		    if (fabs(par0[k2]-par1[k2]) > ptol)
			break;
		if (k2 < int(par0.size()))
		    break;
	    }
	    if (k1 == int(loop_ints[ki].size()))
		return 0;  // Degenerate loop. Do not set coincidence

	    // Fetch the parameter values of the boundary intersections
	    vector<double> loop_par;
	    for (int kj=0; kj<int(loop_ints[ki].size()); kj++)
	    {
		vector<double> parval = loop_ints[ki][kj]->getPar();
		loop_par.insert(loop_par.end(), parval.begin(), parval.end());
	    }

			
	    // Check coincidence
	    is_coincident = checkCoincide(surf1, surf2, loop_par, epsge_);
	    if (!is_coincident)
		return 0;  // Not coincidence

	    // Otherwise mark coincidence
	    int_results_->setCoincidence(loop_ints[ki]);
	}
    }
    return 1;  // Coincidence between the two surfaces
}


//===========================================================================
void SfSfIntersector::microCase()
//===========================================================================
{
    // Purpose: Given two parametric surfaces, try to set a
    // coincident intersection result when no solution is found and
    // the surfaces are too small to subdivide any further

    // Fetch all intersections at the boundaries sorted around the
    // boundaries (using polar coordinates and sort according to the
    // angle).
    vector<shared_ptr<IntersectionPoint> > bd_ints;
    int_results_->getSortedInts(/* AROUND_BOUNDARIES, */ bd_ints);

    if (bd_ints.size() <= 1)
	return;   // At most one intersection point, nothing to do

    if (bd_ints.size() == 2)
    {
	// Make sure that the points are connected
	bd_ints[0]->connectTo(bd_ints[1], MICRO_SFSF);
	return;
    }

    // Check if the points are connected already
    if (isConnected(bd_ints))
	return;

    // Connect all points to a branch point in the middle of the surfaces
    // @@@ VSK. What to do if some points are connected already and some
    // are not. I don't want to make small artificial PAC areas.
    // Everything is by default very tiny so I don't care so much about
    // if the midpoint is the best possble intersection point or that the
    // best intersection points at the boundaries are connected to the
    // midpoint.
    double mid_par[4];
    int_results_->getMidParameter(mid_par);
    shared_ptr<IntersectionPoint> mid_int =
	int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1],
					   epsge_, mid_par, mid_par+2);
    bd_ints[0]->connectTo(mid_int, MICRO_SFSF);

    int ki;
    for (ki=1; ki<int(bd_ints.size()); ki++)
    {
	if (!(int_results_->isConnectedInside(bd_ints[ki], mid_int)))
	    bd_ints[ki]->connectTo(mid_int, MICRO_SFSF);
    }
    
    return;
}


//===========================================================================
bool SfSfIntersector::degTriangleSimple()
//===========================================================================
{
    // Purpose: Given two parametric surfaces where one of them is a
    // small degenerate triangle. Make connections if possible.

    // Fetch ParamSurface information
    ParamSurfaceInt* surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt* surf2 = obj_int_[1]->getParamSurfaceInt();
    ASSERT(surf1 != 0 && surf2 != 0);

    ParamSurface *srf1 = surf1->getParamSurface().get();
    ParamSurface *srf2 = surf2->getParamSurface().get();

    bool triang1 = surf1->getDegTriang();
    bool triang2 = surf2->getDegTriang();
    if (triang1 == false && triang2 == false)
	return false;    // No degenerate tiny triangles

    if (triang1 && !triang2 && surf2->isDegenerate(epsge_->getEpsge()))
	return false;    // Complex situation, no rules

    if (triang2 && !triang1 && surf1->isDegenerate(epsge_->getEpsge()))
	return false;    // Complex situation, no rules

    // Fetch intersection points at all degenerate boundaries
    vector<shared_ptr<IntersectionPoint> > ints_at_deg;
    getPointsAtDegEdges(ints_at_deg);
    if (ints_at_deg.size() == 0)
	return false;    // No points to connect to at degenerate edges.
                         // No need of special treatment

    // Fetch intersection points at the opposite boundaries
    vector<pair<shared_ptr<IntersectionPoint>,int> > other_ints;
    getPointsOppositeDegEdges(other_ints);
    if (other_ints.size() == 0)
	return false;  // TEST

    for (size_t ki=0; ki<other_ints.size(); ki++)
    {
	// Check if the current intersection point is a corner
	// touch. In that case no connection should be performed.
	IntPtClassification type, type2;
	shared_ptr<IntersectionPoint> curr = other_ints[ki].first;
	type = curr->getClassification(srf1, srf2);
	int nmb_branch = curr->numBranches();
	if (type == DIR_TOUCH && nmb_branch == 2)
	{
	    type2 = curr->getClassification(srf1, srf2, 1);
	    if (type == DIR_TOUCH)
		type = type2;
	}
	if (type == DIR_TOUCH)
	    continue;       // No connection

	// Connect to a point at some degenerate boundary
	// First compute boundary specification
	int bd_idx = other_ints[ki].second;
	int obj_idx = (bd_idx < 4) ? 0 : 1;
	int par_idx = (bd_idx < 4) ? bd_idx/2 : (bd_idx-4)/2;
	int idx = obj_idx*2 + (1 - par_idx);
	double par =  (bd_idx % 2 == 0) ? obj_int_[obj_idx]->startParam(par_idx) : 
		obj_int_[obj_idx]->endParam(par_idx);
	vector<shared_ptr<IntersectionPoint> > curr_deg;
	int_results_->getIntersectionPoints(2*obj_idx+par_idx, par, curr_deg);

	// Dismiss the point if already connected
	size_t kr;
	for (kr=0; kr<curr_deg.size(); kr++)
	    if (int_results_->isConnectedInside(curr, curr_deg[kr]))
		break;
	if (kr < curr_deg.size())
	    continue;   // No connection required

	vector<shared_ptr<IntersectionLink> > links;
	double par2 = curr->getPar(idx);
	getLinksAtDegEdges(curr_deg, idx, curr, links);

	if (links.size() > 0)
	{
	    // Fetch most relevant link
	    int lnk_idx = -1;
	    double t1f=-1.0, t2f=-1.0;
	    IntersectionPoint *p1, *p2;
	    double param[4], ppt[4];
	    double tx;
	    double min_dist2 = 1.0e10;

	    // Compute the parameter value of the point at the degenerate edge found
	    // by continuing the direction of the intersection curve in the point at
	    // the opposit edge
	    Point ptang1 = curr->getPar1Dir();
	    Point ptang2 = curr->getPar2Dir();
	    if (obj_idx == 0)
	    {
		tx = (fabs(ptang1[par_idx]) < epsge_->getRelParRes()) ? 1.0 : 
		    (par-curr->getPar(par_idx))/ptang1[par_idx];
		ppt[par_idx] = par;
		ppt[1-par_idx] = curr->getPar(1-par_idx) + tx*ptang1[1-par_idx];
		ppt[2] = curr->getPar(2) + tx*ptang2[0];
		ppt[3] = curr->getPar(3) + tx*ptang2[1];
	    }
	    else
	    {
		tx = (fabs(ptang2[par_idx]) < epsge_->getRelParRes()) ? 1.0 : 
		    (par-curr->getPar(2+par_idx))/ptang2[par_idx];
		ppt[2+par_idx] = par;
		ppt[3-par_idx] = curr->getPar(3-par_idx) + tx*ptang2[1-par_idx];
		ppt[0] = curr->getPar(0) + tx*ptang1[0];
		ppt[1] = curr->getPar(1) + tx*ptang1[1];
	    }

	    for (kr=0; kr<links.size(); kr++)
	    {
		// Relevant link between intersection points at the degenerate edge is found
		// Get endpoints
		links[kr]->getIntersectionPoints(p1, p2);

		// Compute parameter of new intersection point
		double t1, t2;
		double ta = p1->getPar(idx);
		double tb = p2->getPar(idx);
		t1 = (par2 - ta)/(tb - ta);
		t2 = (tb - par2)/(tb - ta);
		for (int kh=0; kh<4; kh++)
		{
		    param[kh] = t2*p1->getPar(kh) + t1*p2->getPar(kh);
		}
		double  dist2 = Utils::distance_squared(param, param+4, ppt);
		if (lnk_idx < 0 || dist2<min_dist2)
		{
		    lnk_idx = (int)kr;
		    min_dist2 = dist2;
		    t1f = t1;
		    t2f = t2;
		}
	    }
	    
	    links[lnk_idx]->getIntersectionPoints(p1, p2);
	    if (fabs(t1f) < curr->parameterTolerance(idx))
	    {
		curr->connectTo(p1, DEG_TRIANGLE);
	    }
	    else if (fabs(t2f) < curr->parameterTolerance(idx))
	    {
		curr->connectTo(p2, DEG_TRIANGLE);
	    }
	    else
	    {
		for (int kh=0; kh<4; kh++)
		{
		    double ta = p1->getPar(kh);
		    double tb = p2->getPar(kh);
		    param[kh] = ppt[kh];
		    if (param[kh] < std::min(ta,tb))
			param[kh] = std::min(ta,tb);
		    if (param[kh] > std::max(ta,tb))
			param[kh] = std::max(ta,tb);
		}
		shared_ptr<IntersectionPoint> tmp = 
		    int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1], 
						       getTolerance(),
						       param, param+2);
		tmp->connectTo(p1, DEG_TRIANGLE);
		tmp->connectTo(p2, DEG_TRIANGLE);
		p1->disconnectFrom(p2);
		curr->connectTo(tmp, DEG_TRIANGLE);
	    }
	}
	else
	{
	    // Connect to the closest intersection point
	    double dist;
	    double mindist = 1.0e6;  // A large number
	    int minind = -1;
	    vector<double> par1 = curr->getPar();
	    for (kr=0; kr<curr_deg.size(); kr++)
	    {
		vector<double> par2 = curr_deg[kr]->getPar();
		dist = 0.0;
		for (int kh=0; kh<(int)(par1.size()); kh++)
		    dist += (par1[kh] - par2[kh])*(par1[kh] - par2[kh]);
		if (dist < mindist)
		{
		    mindist = dist;
		    minind = (int)kr;
		}
	    }
	    if (minind >=0)
		curr->connectTo(curr_deg[minind], DEG_TRIANGLE);
	}
    }
    return true;  // Connections performed
}
    

//===========================================================================
void SfSfIntersector::
getLinksAtDegEdges(vector<shared_ptr<IntersectionPoint> >& deg_pnts, 
		   int idx, shared_ptr<IntersectionPoint> curr,
		   vector<shared_ptr<IntersectionLink> >& links)
//===========================================================================
{
    size_t ki, kj;
    double par = curr->getPar(idx);
    double ptol = curr->parameterTolerance(idx);
    vector<shared_ptr<IntersectionLink> > local;
    vector<shared_ptr<IntersectionLink> >::iterator end;
    CrossesValue crosses(idx, par, ptol);
    TestInDomain ends_in_domain(obj_int_[0].get(), obj_int_[1].get(), curr);
    for (size_t kr=0; kr<deg_pnts.size(); kr++)
    {
	local.clear();
	deg_pnts[kr]->getNeighbourLinks(local);
	end = remove_if(local.begin(), end, not1(crosses));
	end = remove_if(local.begin(), end, not1(ends_in_domain));
	for (ki=0; ki<local.size(); ki++)
	{
	    for (kj=0; kj<links.size(); kj++)
		if (links[kj].get() == local[ki].get())
		    break;
	    if (kj >= links.size())
		links.push_back(local[ki]);
	}
    }
}

//===========================================================================
void SfSfIntersector::
getPointsAtDegEdges(vector<shared_ptr<IntersectionPoint> >& result)
//===========================================================================
{
    // For each degenerate boundary, fetch points
    for (int ki=0; ki<2; ki++)
    {
	for (int dir=0; dir<2; dir++)
	{
	    int deg_edge = obj_int_[ki]->isDegenerate(epsge_->getEpsge(), dir);
	    if (deg_edge == 1 || deg_edge == 3)
	    {
		vector<shared_ptr<IntersectionPoint> > local;
		double par = obj_int_[ki]->startParam(dir);
		int_results_->getIntersectionPoints(2*ki+dir, par, local);
		result.insert(result.end(), local.begin(), local.end());
	    }
	    if (deg_edge == 2 || deg_edge == 3)
	    {
		vector<shared_ptr<IntersectionPoint> > local;
		double par = obj_int_[ki]->endParam(dir);
		int_results_->getIntersectionPoints(2*ki+dir, par, local);
		result.insert(result.end(), local.begin(), local.end());
	    }
	}
    }
}

    
//===========================================================================
void SfSfIntersector::getPointsOppositeDegEdges(
    vector<pair<shared_ptr<IntersectionPoint>, int> >& result)
//===========================================================================
{
    // For each degenerate boundary, fetch points at the opposite boundary
    int bd_idx = 0;
    size_t kr;
    for (int ki=0; ki<2; ki++)
    {
	for (int dir=0; dir<2; dir++, bd_idx+=2)
	{
	    int deg_edge = obj_int_[ki]->isDegenerate(epsge_->getEpsge(), dir);
	    if (deg_edge == 1)
	    {
		vector<shared_ptr<IntersectionPoint> > local;
		double par = obj_int_[ki]->endParam(dir);
		int_results_->getIntersectionPoints(2*ki+dir, par, local);
		for (size_t kj=0; kj<local.size(); kj++)
		{
		    for (kr=0; kr<result.size(); kr++)
			if ((result[kr].first).get() == (local[kj]).get())
			    break;
		    if (kr < result.size())
			continue;  // Point already registered
		    result.push_back(make_pair(local[kj],bd_idx));
		}
	    }
	    if (deg_edge == 2)
	    {
		vector<shared_ptr<IntersectionPoint> > local;
		double par = obj_int_[ki]->startParam(dir);
		int_results_->getIntersectionPoints(2*ki+dir, par, local);
		for (size_t kj=0; kj<local.size(); kj++)
		{
		    for (kr=0; kr<result.size(); kr++)
			if ((result[kr].first).get() == (local[kj]).get())
			    break;
		    if (kr < result.size())
			continue;  // Point already registered
		    result.push_back(make_pair(local[kj],bd_idx+1));
		}
	    }
	}
    }
}

    
//===========================================================================
int SfSfIntersector::updateIntersections()
//===========================================================================
{
    // Purpose: Given two parametric surfaces in a simple case
    // situation, connect intersection points at the surface
    // boundaries into tracks to represent the topology of the
    // intersection problem
    // Updated by VSK, 0806. Correct sorting direction if simple case
    // is found by implicitization

    // Fetch ParamSurface information
    ParamSurfaceInt* surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt* surf2 = obj_int_[1]->getParamSurfaceInt();
    ASSERT(surf1 != 0 && surf2 != 0);

    ParamSurface *srf1 = surf1->getParamSurface().get();
    ParamSurface *srf2 = surf2->getParamSurface().get();

    // Simple case is detected by checking overlap in the normal cones
    // of the surface. Sort the intersection points at the boundaries
    // of this surface with respect to the normal vector of the plane
    // spanned by the cone axis.
    Point vec;
    if (sorting_obj_[0] < 0)
    {
	DirectionCone cone1, cone2;
	try {
	    cone1 = obj_int_[0]->directionCone();
	    cone2 = obj_int_[1]->directionCone();
	} catch (...) {
	// What to do here
	}
	Point axis1 = cone1.centre();
	Point axis2 = cone2.centre();
	vec = axis1 % axis2;
    }
    else
    {
	vec = sorting_dir_[0];

	// The orientation of the sorting direction may be wrong. Turn
	// the direction of the vector if necessary
	bool is_turned;
	is_turned = int_results_->checkSortingDir(vec, sorting_obj_[0]);
    }

    // Fetch the sorted intersection points
    vector<shared_ptr<IntersectionPoint> > bd_ints;
    int_results_->getSortedBdInts(vec, bd_ints, sorting_obj_[0]);


    // Classify the intersection points
    int ki;
    int nmb_nottouch = 0;
    int nmb_ints;
    int nmb_in = 0;
    int nmb_out = 0;
    int nmb_perpendicular = 0;
    int nmb_branch;
    IntPtClassification type, type2;
    vector<pair<shared_ptr<IntersectionPoint>,
	IntPtClassification> > bd_ints_typed;
    nmb_ints = (int)bd_ints.size();
    if (nmb_ints == 0)
	return 1;  // No points to connect. Finished.

    for (ki=0; ki<nmb_ints; ki++) {
	type = bd_ints[ki]->getClassification(srf1, srf2);

	// Check if another intersection branch exist
	nmb_branch = bd_ints[ki]->numBranches();
	if (nmb_branch == 2) {
	    type2 = bd_ints[ki]->getClassification(srf1, srf2, 1);
	    if (type == DIR_TOUCH) {
		type = type2;
	    } else if (type2 != DIR_TOUCH) {
		bd_ints_typed.push_back(make_pair(bd_ints[ki],type2));
		if (type2 != DIR_TOUCH) {
		    nmb_nottouch++;
		}
		if (type2 == DIR_IN) {
		    nmb_in++;
		}
		if (type2 == DIR_OUT) {
		    nmb_out++;
		}
	    }
	}
	bd_ints_typed.push_back(make_pair(bd_ints[ki],type));
	if (type != DIR_TOUCH) {
	    nmb_nottouch++;
	}
	if (type == DIR_IN) {
	    nmb_in++;
	}
	if (type == DIR_OUT) {
	    nmb_out++;
	}
	if (type == DIR_PERPENDICULAR) {
	    nmb_perpendicular++;
	}
    }

    // Move corner points to the end of the array
    if (nmb_nottouch < int(bd_ints_typed.size())) {
	for (ki=0; ki<nmb_nottouch; ki++) {
	    if (bd_ints_typed[ki].second == DIR_TOUCH) {
		pair<shared_ptr<IntersectionPoint>, IntPtClassification> curr
		    = bd_ints_typed[ki];
		bd_ints_typed.erase(bd_ints_typed.begin()+ki);
		bd_ints_typed.push_back(curr);
		ki--;
	    }
	}
    }

    if (selfint_case_)
    {
	// The surface-surface intersection is performed in a self intersection
	// setting. If the two surfaces are neighbours, superfluous intersection
	// points/curves will occur along the common boundary between the two
	// surfaces. Remove these extra intersections from the vector of boundary
	// intersections as they disturbe the configuration.
	vector<BoundaryIntersectionData> common_bd;
	int_results_->intersectAlongCommonBoundary(0.0, common_bd);
	for (size_t h1=0; h1<common_bd.size(); h1++)
	    for (size_t h2=0; h2<common_bd[h1].pts.size(); h2++)
	    {
		// Look for the point on the common boundary in the vector
		// of boundary intersections
		size_t h3;
		for (h3=0; h3<bd_ints_typed.size(); h3++)
		    if (bd_ints_typed[h3].first.get() == common_bd[h1].pts[h2].get())
			break;
		if (h3 < bd_ints_typed.size() && bd_ints_typed[h3].second ==  DIR_HIGHLY_SINGULAR)
		{
		    // False intersection, remove
		    bd_ints_typed.erase(bd_ints_typed.begin()+h3);
		    bd_ints.erase(bd_ints.begin()+h3);
		    nmb_nottouch--;
		    nmb_ints--;
		}
	    }
    }

    // Perform connection of intersection points into track according
    // to the current situation.
    // First check if more than one relevant point exist
    if ((nmb_ints == 1)
	|| ((nmb_nottouch == 1)
	    && !((nmb_ints == 2)
		 && (nmb_in == 1 || nmb_out == 1
		     || nmb_perpendicular == 1)))) {
	return 1;   // Finished, not necessary to connect
    }

    // Try to make connections between the intersection points
    // according to their sequence. The classification must be
    // consistent with the sequence
    if (connectDirected(bd_ints_typed, nmb_nottouch))
	return 1;    // Connections performed

    // Then check if the intersection points are connected already
    if (isConnected(bd_ints_typed, nmb_nottouch))
	return 1;   // Nothing more to do

    // VSK, 0709. The most believed connection rules are tryed without
    // success. Check if an alternative sorting direction exist
    if (sorting_obj_[0] >= 0 && sorting_obj_[1] >= 0)
    {
	vec = sorting_dir_[1];

	// The orientation of the sorting direction may be wrong. Turn
	// the direction of the vector if necessary
	bool is_turned;
	is_turned = int_results_->checkSortingDir(vec, sorting_obj_[1]);

	bd_ints.clear();
	int_results_->getSortedBdInts(vec, bd_ints, sorting_obj_[1]);
	nmb_nottouch = nmb_in = nmb_out = nmb_perpendicular = 0;

	bd_ints_typed.clear();

	for (ki=0; ki<nmb_ints; ki++) {
	    type = bd_ints[ki]->getClassification(srf1, srf2);

	    // Check if another intersection branch exist
	    nmb_branch = bd_ints[ki]->numBranches();
	    if (nmb_branch == 2) {
		type2 = bd_ints[ki]->getClassification(srf1, srf2, 1);
		if (type == DIR_TOUCH) {
		    type = type2;
		} else if (type2 != DIR_TOUCH) {
		    bd_ints_typed.push_back(make_pair(bd_ints[ki],type2));
		    if (type2 != DIR_TOUCH) {
			nmb_nottouch++;
		    }
		    if (type2 == DIR_IN) {
			nmb_in++;
		    }
		    if (type2 == DIR_OUT) {
			nmb_out++;
		    }
		}
	    }
	    bd_ints_typed.push_back(make_pair(bd_ints[ki],type));
	    if (type != DIR_TOUCH) {
		nmb_nottouch++;
	    }
	    if (type == DIR_IN) {
		nmb_in++;
	    }
	    if (type == DIR_OUT) {
		nmb_out++;
	    }
	    if (type == DIR_PERPENDICULAR) {
		nmb_perpendicular++;
	    }
	}

	// Move corner points to the end of the array
	if (nmb_nottouch < int(bd_ints_typed.size())) {
	    for (ki=0; ki<nmb_nottouch; ki++) {
		if (bd_ints_typed[ki].second == DIR_TOUCH) {
		    pair<shared_ptr<IntersectionPoint>, IntPtClassification> curr
			= bd_ints_typed[ki];
		    bd_ints_typed.erase(bd_ints_typed.begin()+ki);
		    bd_ints_typed.push_back(curr);
		    ki--;
		}
	    }
	}

	if (selfint_case_)
	{
	    // The surface-surface intersection is performed in a self intersection
	    // setting. If the two surfaces are neighbours, superfluous intersection
	    // points/curves will occur along the common boundary between the two
	    // surfaces. Remove these extra intersections from the vector of boundary
	    // intersections as they disturbe the configuration.
	    vector<BoundaryIntersectionData> common_bd;
	    int_results_->intersectAlongCommonBoundary(0.0, common_bd);
	    for (size_t h1=0; h1<common_bd.size(); h1++)
		for (size_t h2=0; h2<common_bd[h1].pts.size(); h2++)
		{
		    // Look for the point on the common boundary in the vector
		    // of boundary intersections
		    size_t h3;
		    for (h3=0; h3<bd_ints_typed.size(); h3++)
			if (bd_ints_typed[h3].first.get() == common_bd[h1].pts[h2].get())
			    break;
		    if (h3 < bd_ints_typed.size() && bd_ints_typed[h3].second ==  DIR_HIGHLY_SINGULAR)
		    {
			// False intersection, remove
			bd_ints_typed.erase(bd_ints_typed.begin()+h3);
			bd_ints.erase(bd_ints.begin()+h3);
			nmb_nottouch--;
			nmb_ints--;
		    }
		}
	}

	// Perform connection of intersection points into track according
	// to the current situation.
	// First check if more than one relevant point exist
	if ((nmb_ints == 1)
	    || ((nmb_nottouch == 1)
		&& !((nmb_ints == 2)
		     && (nmb_in == 1 || nmb_out == 1
			 || nmb_perpendicular == 1)))) {
	    return 1;   // Finished, not necessary to connect
	}

	// Try to make connections between the intersection points
	// according to their sequence. The classification must be
	// consistent with the sequence
	if (connectDirected(bd_ints_typed, nmb_nottouch))
	    return 1;    // Connections performed

	// Then check if the intersection points are connected already
	if (isConnected(bd_ints_typed, nmb_nottouch))
	    return 1;   // Nothing more to do

    }	

    // If only two relevant (not corner touching) points exist and one
    // points in or one points out, the two points must be
    // connected. A corner touch is overruled if one point points in
    // or out
    if ((nmb_ints == 2)
	&& ((nmb_in + nmb_out == 2) || nmb_in == 1 || nmb_out == 1)) {
	// Connect the two non-touching points
	int idx = std::max(nmb_nottouch-1,1);
	bd_ints_typed[0].first->connectTo(bd_ints_typed[idx].first,
					  SIMPLE_TWO_POINTS);
	return 1;
    }

    if ((nmb_ints == 2)
	&& ((nmb_in + nmb_out + nmb_perpendicular == 2)
	    || nmb_in == 1 || nmb_out == 1 || nmb_perpendicular == 1)) {
	// Check if a connection between the two non-touch points is
	// correct
	int idx = std::max(nmb_nottouch-1,1);
	if (canConnect(bd_ints_typed[0].first, bd_ints_typed[idx].first)) {
	    bd_ints_typed[0].first->connectTo(bd_ints_typed[idx].first,
					      SIMPLE_TWO_POINTS);
	}
	return 1;
    }

    // VSK, 070529
    // Check if connected points in sequence both point in or out.
    // In that case remove the less accurate one
    // We now try to iterate on the intersection points to see if they
    // move away from the edge.
    int num_before = int_results_->numIntersectionPoints();
    for (ki=1; ki < (int)bd_ints.size(); ki++)
    {
	if ((bd_ints_typed[ki-1].second == DIR_IN && bd_ints_typed[ki].second == DIR_IN) ||
	     (bd_ints_typed[ki-1].second == DIR_OUT && bd_ints_typed[ki].second == DIR_OUT))
	{
	    if (bd_ints_typed[ki-1].first->isConnectedTo(bd_ints_typed[ki].first) ||
		canConnect(bd_ints_typed[ki-1].first, bd_ints_typed[ki].first))
	    {
		if (bd_ints_typed[ki-1].first->getDist() > bd_ints_typed[ki].first->getDist())
		{
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfSfIntersector::updateIntersections 1. Removed point at: ";
			std::cout << bd_ints[ki-1]->getPar(0) << " ";
			std::cout << bd_ints[ki-1]->getPar(1) << " ";
			std::cout << bd_ints[ki-1]->getPar(2) << " ";
			std::cout << bd_ints[ki-1]->getPar(3);
			std::cout << std::endl;
		    }
		    int_results_->removeIntPoint(bd_ints_typed[ki-1].first);
		    bd_ints_typed.erase(bd_ints_typed.begin()+ki-1);
		    bd_ints.erase(bd_ints.begin()+ki-1);
		}
		else
		{
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfSfIntersector::updateIntersections 1. Removed point at: ";
			std::cout << bd_ints[ki]->getPar(0) << " ";
			std::cout << bd_ints[ki]->getPar(1) << " ";
			std::cout << bd_ints[ki]->getPar(2) << " ";
			std::cout << bd_ints[ki]->getPar(3);
			std::cout << std::endl;
		    }
		    int_results_->removeIntPoint(bd_ints_typed[ki].first);
		    bd_ints_typed.erase(bd_ints_typed.begin()+ki);
		    bd_ints.erase(bd_ints.begin()+ki);
		}
		ki = 0;
	    }
	}
    }
    nmb_ints = (int)bd_ints.size();

    int nmb_orig = 0;
    /*Point axis_x(1.0, 0.0);
    Point axis_y(0.0, 1.0);
    Point tangent2d;
    double maxcosangle = 0.0;
    int dir = 0;
    for (int i = 0; i < nmb_ints; ++i) {
	if (!bd_ints[i]->hasUniqueTangentDirection()) {
	    continue;
	}
	tangent2d = bd_ints[i]->getPar1Dir();
	if (fabs(axis_x.cosAngle(tangent2d)) > maxcosangle) {
	    maxcosangle = fabs(axis_x.cosAngle(tangent2d));
	    dir = 0;
	}

	if (fabs(axis_y.cosAngle(tangent2d)) > maxcosangle) {
	    maxcosangle = fabs(axis_y.cosAngle(tangent2d));
	    dir = 1;
	}
	tangent2d = bd_ints[i]->getPar2Dir();
	if (fabs(axis_x.cosAngle(tangent2d)) > maxcosangle) {
	    maxcosangle = fabs(axis_x.cosAngle(tangent2d));
	    dir = 2;
	}

	if (fabs(axis_y.cosAngle(tangent2d)) > maxcosangle) {
	    maxcosangle = fabs(axis_y.cosAngle(tangent2d));
	    dir = 3;
	}
    }
    if (maxcosangle > 0.0) {*/
    for (int dir=0; dir<4; dir++)
    {
	cout << "dir = " << dir << endl;
	double par = 0.0; // Doesn't matter since along == false
	bool along = false;
	bool only_isolated = false;
	postIterate2(nmb_orig, dir, par, along, only_isolated);
	// Points that have moved to the inner need to be removed
	vector<shared_ptr<IntersectionPoint> > int_pts;
	int_results_->getIntersectionPoints(int_pts);
	int num_ints = (int)int_pts.size();
	double tol = epsge_->getRelParRes();
	for (int i = 0; i < num_ints; ++i) {
	    if (!int_pts[i]->hasUniqueTangentDirection())
		continue;
	    const double* par1 = int_pts[i]->getPar1();
	    const double* par2 = int_pts[i]->getPar2();
	    if (!obj_int_[0]->boundaryPoint(par1, tol)
		&& !obj_int_[1]->boundaryPoint(par2, tol)) {
		if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		{
		    std::cout << "SfSfIntersector::updateIntersections 2. Removed point at: ";
		    std::cout << int_pts[i]->getPar(0) << " ";
		    std::cout << int_pts[i]->getPar(1) << " ";
		    std::cout << int_pts[i]->getPar(2) << " ";
		    std::cout << int_pts[i]->getPar(3);
		    std::cout << std::endl;
		}
		if (int_pts[i]->numNeighbours() <= 2) {
		    int_results_->removeIntPoint(int_pts[i]);
		}
		else {
		    int_results_->removeIntPoint2(int_pts[i]);
		}
	    }
	}
    }
	int num_after = int_results_->numIntersectionPoints();
	if (num_after < num_before) {
	    // Recursive call!
	    cout << "Recursive call to update " << endl;
	    int stat = updateIntersections();
	    return stat;
	}
	//return 1;

    // This situation is not handled. Write debug output.
    if ((getenv("DEBUG") && (*getenv("DEBUG"))=='1') ||
	(getenv("DEBUG_FINISH") && (*getenv("DEBUG_FINISH"))=='1'))
    {
	cout << "Inconsistent conditions for connections found " << endl;
	writeDebugConnect(bd_ints_typed);
    }
    handleComplexity();

    /*if (selfint_case_)
    {
	// Another try
	selfintComplex2();
	}*/

    return 0;
}


//===========================================================================
int SfSfIntersector::repairIntersections()
//===========================================================================
{
    // Purpose: This function may be called at the end of the
    // intersection algorithm to repair incorrect holes, branch
    // points, etc., in the intersection curves.

     iterateOnIntersectionPoints();

     removeIsolatedPoints();

     int_results_->removeDoublePoints();

     repairSingularityBox();
 
//     repairFalseBranchPoints();

     int_results_->removeLooseEnds();

     int_results_->removeDefectLinks();

    repairFalseBranchPoints();

     int_results_->removeFalseCorners();

    fixCrossingLinks();

    repairMissingLinks();

    removeIsolatedPoints();

   return 0;
}


/// Utility in surface-surface intersection
struct intersection_pair
{
    shared_ptr<IntersectionPoint> pt1;
    shared_ptr<IntersectionPoint> pt2;
    double dist;

    intersection_pair(shared_ptr<IntersectionPoint> p1, 
		      shared_ptr<IntersectionPoint> p2, 
		      double d)
	{
	    pt1 = p1;
	    pt2 = p2;
	    dist = d;
	}
};

bool compare_pairs(const intersection_pair& p1, const intersection_pair& p2)
{
    return (p1.dist < p2.dist);
}

//===========================================================================
void SfSfIntersector::repairSingularityBox()
//===========================================================================
{
    // Cleans up areas around branch points. The method is: 1)
    // Identify a branch point (the singularity), 2) Define a suitable
    // box around the branch point, 3) Run intersection algorithms on
    // the boundary curves, 4) Check if boundary intersections are
    // consistent, and if so 5) Connect up with the branch point.

    vector<shared_ptr<IntersectionPoint> > int_pts
	= int_results_->getIntersectionPoints();
    int int_sz = int(int_pts.size());

    // Get vector of branch points. Also include near-singularities.
    vector<shared_ptr<IntersectionPoint> > branch_pts;
    for (int i = 0; i < int_sz; ++i) {
	if (int_pts[i]->getSingularityType() == BRANCH_POINT
	    || int_pts[i]->isNearSingular()) {
	    branch_pts.push_back(int_pts[i]);
	}
    }
    int nbranchpts = (int)branch_pts.size();

    typedef vector<shared_ptr<IntersectionPoint> >::iterator iter;

    // Loop through branch points
    double ptol = 10.0*epsge_->getRelParRes();
    for (int i = 0; i < nbranchpts; ++i) {
	// First define the box around the singularity
	double frompar[4];
	double topar[4];
	getSingularityBox(branch_pts[i], frompar, topar);

	// Check box
	int kr;
	for (kr=0; kr<4; kr++)
	    if (topar[kr]-frompar[kr] < ptol)
		break;
	if (kr < 4)
	    continue;  // Illegal singularity box found

	// Check if other branch points are in the singularity box. If
	// so, remove them.
	iter it = branch_pts.begin() + i + 1;
	while (it != branch_pts.end()) {
	    if ((*it)->isInDomain(frompar, topar)) {
		it = branch_pts.erase(it);
		--nbranchpts;
	    }
	    else {
		++it;
	    }
	}

	// Check if the singularity box actually needs fixing
	int numinside = 0; // At least the branch point is inside
	int_pts = int_results_->getIntersectionPoints();
	int_sz = int(int_pts.size());
	for (int j = 0; j < int_sz; ++j) {
	    if (int_pts[j]->isInDomain(frompar, topar)) {
		++numinside;
	    }
	}
	if (numinside == 1) {
	    // Only the branch point is inside - continue
	    continue;
	}

	// Set up the subintersector and get intersections on the
	// boundaries.
	shared_ptr<SfSfIntersector> intersector
	    = getSubIntersector(frompar, topar);
	intersector->getBoundaryIntersections();
	intersector->iterateOnIntersectionPoints();

	// If, for some reason, there are more boundary points on the
	// same branch, we must make sure that only one gets connected
	// to the branch point. We do that by removing all but one
	// boundary point on the branch. Which one is left is
	// arbitrary. Also: If, for some reason, we get a branch point
	// among the boundary points, we assume that it is another
	// version of the original branch point, and therefore remove
	// it.
	vector<shared_ptr<IntersectionPoint> > bdints
	    = intersector->getIntPool()->getIntersectionPoints();
	it = bdints.begin();
	while (it != bdints.end()) {
	    // Remove if branch point
	    if ((*it)->getSingularityType() == BRANCH_POINT || 
		(*it)->getSingularityType() == HIGHER_ORDER_POINT ||
		(*it)->getSingularityType() == ISOLATED_POINT) {
		it = bdints.erase(it);
		continue;
	    }
	    // Check for more than one point on the same branch
	    Point tan1 = (*it)->getTangent();
	    Point diff1 = (*it)->getPoint() - branch_pts[i]->getPoint();
	    iter jt = it + 1;
	    while (jt != bdints.end()) {
		if ((*jt)->getSingularityType() == BRANCH_POINT || 
		    (*jt)->getSingularityType() == HIGHER_ORDER_POINT ||
		    (*jt)->getSingularityType() == ISOLATED_POINT) 
		{
		    ++jt;  // No tangent exist
		}
		else
		{
		    Point tan2 = (*jt)->getTangent();
		    double cos = tan1.cosAngle(tan2);

		    // Check if the two intersection points lies at the
		    // same boundary. If not, it can be dangerous to remove
		    // one of them
		    bool diff_bd = intersector->getIntPool()->atDifferentBoundary(*it, *jt);

		    // cos = 0.97 is roughly 15 degrees
		    if (cos > 0.97 && !diff_bd) {
			jt = bdints.erase(jt);
		    }
		    else {
			++jt;
		}
		}
	    }
	    ++it;
	}

	// Store the intersection points outside the subintersector's
	// domain that links to points on the inside.
	vector<IntersectionPoint*> outside_nbours;
	int_pts = int_results_->getIntersectionPoints();
	int_sz = int(int_pts.size());
	for (int j = 0; j < int_sz; ++j) {
	    if (int_pts[j]->isInDomain(frompar, topar)) {
		vector<IntersectionPoint*> neighbours;
		int_pts[j]->getNeighbours(neighbours);
		int nneighbours = (int)neighbours.size();
		for (int k = 0; k < nneighbours; ++k) {
		    if (!neighbours[k]->isInDomain(frompar, topar)) {
			outside_nbours.push_back(neighbours[k]);
		    }
		}
	    }
	}
	int noutside = (int)outside_nbours.size();

	// Remove all intersection points in 'this' pool that lies
	// within the subintersector's domain, replace the branch
	// point (since it was removed). Add the new boundary
	// intersections. Connect them to the branch point, and try to
	// connect them to the previously stored outside points.
	// VSK, 0607. Keep boundary points
	int_results_->removeIntPoints(frompar, topar/*, true*/);
	int_results_->add_point_and_propagate_upwards(branch_pts[i]);

	int nbdints = (int)bdints.size();
	vector<intersection_pair> int_pairs;
	for (int j = 0; j < nbdints; ++j) {
	    vector<double> par = bdints[j]->getPar();
	    shared_ptr<IntersectionPoint> newpt
		= int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1],
						     epsge_, &par[0], &par[2]);
	    if (newpt->getSingularityType() == HIGHER_ORDER_POINT ||
		newpt->getSingularityType() == ISOLATED_POINT)
		continue;

	    newpt->connectTo(branch_pts[i], BRANCH_CONNECTION);


	    // Sort candidate connections to the outside points according to
	    // consistency in tangent direction and position
	    Point curr1 = newpt->getPoint();
	    Point tangent1 = newpt->getTangent();
	    for (int k = 0; k < noutside; ++k) 
	    {
		int kh;
		double angle;
		for (kh=0; kh<int_sz; kh++)
		    if (int_pts[kh].get() == outside_nbours[k])
			break;
		if (kh == int_sz)
		    continue;

		if (int_pts[kh]->getSingularityType() == HIGHER_ORDER_POINT ||
		    int_pts[kh]->getSingularityType() == ISOLATED_POINT)
		    angle = M_PI/2.0;
		else
		{
		    Point curr2 = int_pts[kh]->getPoint();
		    Point tangent2 = int_pts[kh]->getTangent();
		    Point diff = curr1 - curr2;
		    diff.normalize();
		    double ta1 = 1.0 - tangent1*tangent2;
		    double ta2 = 1.0 - 0.5*(fabs(diff*tangent1) + fabs(diff*tangent2));
		    angle = ta1 + ta2;
		}
		intersection_pair curr_pair(newpt, int_pts[kh], angle);
		int_pairs.push_back(curr_pair);
	    }
	}

	std::sort(int_pairs.begin(), int_pairs.end(), compare_pairs);
		
	for (int j = 0; j < nbdints; ++j) {
	    for (int k = 0; k < (int)int_pairs.size(); ++k) {
		if (int_pairs[k].pt1->numNeighbours() >= 2 ||
		    int_pairs[k].pt2->numNeighbours() >= 2)
		    continue; // Sufficient connections already performed

		shared_ptr<IntersectionLink> link
		    = int_pairs[k].pt1->connectTo(int_pairs[k].pt2,
				       INSIDE_OUTSIDE_SINGULARITY_BOX);
		bool link_ok = int_results_->verifyIntersectionLink(link);
		if (!link_ok) {
		    int_pairs[k].pt1->disconnectFrom(int_pairs[k].pt2.get());
		}
	    }
	}

    }

    return;
}


//===========================================================================
void SfSfIntersector::repairFalseBranchPoints()
//===========================================================================
{
    // Repair false branch points. The method is: 1) Identify ordinary
    // intersection points with more than two neighbours that are in
    // the inner and not on a degenerate edge, 2) Find the tangent to
    // the intersection curve, 3) Divide the neighbours into those in
    // the "forward" and "backward" directions, 4) Keep only the link
    // to the best forward and backward neighbour, based on the angles
    // of the links with the tangents.

    // Get intersection point info
    vector<IntPtInfo> int_pt_info;
//     double tol = epsge_->getRelParRes();
    bool is_ok = int_results_->checkIntersectionPoints(int_pt_info);
    if (is_ok)
	return;

    // Main loop
    vector<shared_ptr<IntersectionPoint> > int_pts
	= int_results_->getIntersectionPoints();
    int int_sz = int(int_pts.size());
    for (int i = 0; i < int_sz; ++i) {
	bool false_branch_point = !int_pt_info[i].is_ok
	    && int_pt_info[i].nneighbours > 2
	    && (int_pt_info[i].singularity_type == ORDINARY_POINT
		|| int_pt_info[i].singularity_type == TANGENTIAL_POINT)
	    && int_pt_info[i].location == LOC_INSIDE_BOTH
	    && !int_pts[i]->isNearSingular()
	    && !int_pts[i]->isDegenerate();
	if (false_branch_point) {
	    // Get tangent in space
	    Point tangent = int_pts[i]->getTangent();
	    tangent.normalize();
	    // Loop over neighbours to find the best forward and
	    // backward
	    vector<IntersectionPoint*> neighbours;
	    int_pts[i]->getNeighbours(neighbours);
	    int nneigh = int(neighbours.size());
	    int best_forward = -1; // Initialize to "none"
	    int best_backward = -1; // Initialize to "none"
	    double best_forward_angles = 0.0;
	    double best_backward_angles = 0.0;
	    int j;
	    for (j = 0; j < nneigh; ++j) {
		Point point = neighbours[j]->getPoint()
		    - int_pts[i]->getPoint();
		if (point.length() == 0.0)
		    continue;
		point.normalize();
		Point other_tangent;
		try {
		other_tangent= neighbours[j]->getTangent();
		}
		catch (...) {
		    break;
		}

		if (other_tangent.length() != 0.0)
		    other_tangent.normalize();
		double cos = tangent * point;
		double other_cos = other_tangent * point;
		bool second_branch = (neighbours[j]->numBranches() == 2);
		if (second_branch) {
		    other_tangent = neighbours[j]->getTangent(second_branch);
		    if (other_tangent.length() != 0.0)
			other_tangent.normalize();
		    double tmp = other_tangent * point;
		    if (fabs(other_cos) < fabs(tmp)) {
			other_cos = tmp;
		    }
		}
		double angles = cos * cos + other_cos * other_cos;
		// Is this neighbour "forward" or "backward"
		bool forward = (cos > 0.0) ? true : false;
		// Keep track of best forward and backward point
		if (forward) {
		    if (angles > best_forward_angles) {
			best_forward = j;
			best_forward_angles = angles;
		    }
		}
		else {
		    if (angles > best_backward_angles) {
			best_backward = j;
			best_backward_angles = angles;
		    }
		}
	    }
	    if (j < nneigh)
		continue;  // Cannot check branch point

	    // Finally, disconnect neighbours other than the best
	    // forward and backward
	    for (int j = 0; j < nneigh; ++j) {
		if (!(j == best_forward || j == best_backward)) {
		    int_pts[i]->disconnectFrom(neighbours[j]);
		}
	    }
	}
    }

    return;
}


//===========================================================================
void SfSfIntersector::removeIsolatedPoints()
//===========================================================================
{
    // Remove isolated points. Method: 1) Identify ordinary
    // intersection points with no neighbours, 2) Remove these from
    // the pool.

    // Get intersection point info
    vector<IntPtInfo> int_pt_info;
//     double tol = epsge_->getRelParRes();
    bool is_ok = int_results_->checkIntersectionPoints(int_pt_info);
    if (is_ok)
	return;

    // Main loop
    vector<shared_ptr<IntersectionPoint> > int_pts
	= int_results_->getIntersectionPoints();
    typedef vector<shared_ptr<IntersectionPoint> >::iterator iter;
    iter it = int_pts.begin();
    int i = 0;
    while (it != int_pts.end()) {
	bool remove = !int_pt_info[i].is_ok
	    && int_pt_info[i].nneighbours == 0
	    && (int_pt_info[i].singularity_type == ORDINARY_POINT
		|| int_pt_info[i].singularity_type == TANGENTIAL_POINT);

	if (remove && int_results_->isBoundaryPoint(*it))
	    remove = false;

	if (remove) {
	    int_results_->removeIntPoint(*it);
	    it = int_pts.erase(it);
	}
	else {
	    ++it;
	}
	++i;
    }

    return;
}

//===========================================================================
void SfSfIntersector::repairMissingLinks()
//===========================================================================
{
    // Repair holes. The method for this is: 1) Identify ordinary
    // intersection points in the pool with only one link, 2) Set up
    // subproblems for each pair of intersection points such that
    // these are in the corners of the parameter domains, 3) Test for
    // simple case, and 4) If simple case, update the intersections.

    // Get the intersection points but keep only non-OK, ordinary
    // points with zero or one neighbour
    vector<shared_ptr<IntersectionPoint> > int_pts
	= int_results_->getIntersectionPoints();
    IntPtInfo info;
    typedef vector<shared_ptr<IntersectionPoint> >::iterator iter;
    iter it = int_pts.begin();
    while (it != int_pts.end()) {
	(*it)->checkIntersectionPoint(info);
	bool keep = !info.is_ok
	    && info.nneighbours <= 1
	    && (info.singularity_type == ORDINARY_POINT
		|| info.singularity_type == TANGENTIAL_POINT);
	if (!keep) {
	    it = int_pts.erase(it);
	}
	else {
	    ++it;
	}
    }
    if (int_pts.empty()) {
	return;
    }

/*    // Loop through all pairs of intersection points and check
    // subproblems
    it = int_pts.begin();
    bool link_ok = false;
    while (it != int_pts.end()) {
	shared_ptr<IntersectionPoint> ipoint = *it;
	vector<double> ipar = ipoint->getPar();
	iter jt = it + 1;    
	while (jt != int_pts.end()) {
	    shared_ptr<IntersectionPoint> jpoint = *jt;
	    vector<double> jpar = jpoint->getPar();

	    shared_ptr<IntersectionLink> link
		= ipoint->connectTo(jpoint, REPAIRED_MISSING_LINK);
	    link_ok = int_results_->verifyIntersectionLink(link);
	    if (!link_ok) {
		ipoint->disconnectFrom(jpoint.get());
	    }
	    ++jt;
	}
	++it;
	}*/

    if (int_pts.size() == 1) {
	cout << "Cannot repair a single loose end in "
	     << "SfSfIntersector::repairIntersections()..." << endl;
    }

    // VSK, 0607. Try to connect loose ends
    vector<intersection_pair> int_pairs;
    size_t ki, kj;
    double d1;
    for (ki=0; ki<int_pts.size(); ki++)
	for (kj=ki+1; kj<int_pts.size(); kj++)
	{
	    Point curr = int_pts[ki]->getPoint();
	    d1 = curr.dist2(int_pts[kj]->getPoint());

	    intersection_pair curr_pair(int_pts[ki], int_pts[kj], d1);
	    int_pairs.push_back(curr_pair);
	}

    // Sort according to distance
    std::sort(int_pairs.begin(), int_pairs.end(), compare_pairs);

    // Move links containing points with bad accuracy to the end of the candidate
    // list
    double tol = 0.01*epsge_->getEpsge();
    int kr;
    int nmb_pairs = (int)int_pairs.size();
    for (kr=0; kr<nmb_pairs; kr++)
    {
	if (int_pairs[kr].pt1->getDist() > tol || 
	    int_pairs[kr].pt2->getDist() > tol)
	{
	    intersection_pair curr_pair = int_pairs[kr];
	    int_pairs.erase(int_pairs.begin()+kr);
	    int_pairs.push_back(curr_pair);
	    kr--;
	    nmb_pairs--;
	}
    }

    // For each pair of candidate intersection points for connection, check
    // if a connection is possible
    for (ki=0; ki<int_pairs.size(); ki++)
    {
	if (int_pairs[ki].pt1->numNeighbours() >= 2 ||
	    int_pairs[ki].pt2->numNeighbours() >= 2)
	    continue; // Sufficient connections already performed

	// bool connected = connectIfPossible(int_pairs[ki].pt1, int_pairs[ki].pt2);
    }

    /*size_t ki, kj, kr;
    double d1, d2;
    for (ki=0; ki<int_pts.size(); ki++)
    {
	// Fetch best point to start with
	int idx = ki;
	d1 = int_pts[idx]->getDist();
	for (kj=ki+1; kj<int_pts.size(); kj++)
	{
	    d2 = int_pts[kj]->getDist();
	    if (d2 < d1)
	    {
		idx = kj;
		d1 = d2;
	    }
	}
	if (idx != ki)
	    std::swap(int_pts[ki], int_pts[idx]);

	Point curr = int_pts[ki]->getPoint();
	// Sort remaining intersection points according to their distance to the current
	for (kj=ki+1; kj<int_pts.size(); kj++)
	{
	    d1 = curr.dist2(int_pts[kj]->getPoint());
	    for (kr=kj+1; kr<int_pts.size(); kr++)
	    {
		d2 = curr.dist2(int_pts[kr]->getPoint());
		if (d2 < d1)
		{
		    std::swap(int_pts[kj], int_pts[kr]);
		    d1 = d2;
		}
	    }
	}

	// Connect to the closest point that provides a legal connection
	for (kj=ki+1; kj<int_pts.size(); kj++)
	{
	    bool connected = connectIfPossible(int_pts[ki], int_pts[kj]);
	    if (connected)
	    {
		// Remove points without loose ends from the array of candidates
		if (int_pts[kj]->numNeighbours() >= 2)
		{
		    int_pts.erase(int_pts.begin()+kj);
		}
		if (int_pts[ki]->numNeighbours() >= 2)
		{
		    int_pts.erase(int_pts.begin()+ki);
		    ki--;
		}
		break;
	    }
	}
	}*/

    return;
}


//===========================================================================
bool SfSfIntersector::connectIfPossible(shared_ptr<IntersectionPoint> pt1,
					shared_ptr<IntersectionPoint> pt2)
//===========================================================================
{
    // First check if there exist and intersection between these two points
    if (canConnect(pt1, pt2))
    {
	// Fetch already existing links
	set<shared_ptr<IntersectionLink> > links;
	int_results_->getAllLinks(links);
	
	// Check consistency of tangent directions
	shared_ptr<IntersectionLink> link
	    = pt1->connectTo(pt2, REPAIRED_MISSING_LINK);
	bool link_ok = int_results_->verifyIntersectionLink(link, 1);
	if (!link_ok) 
	{
	    pt1->disconnectFrom(pt2.get());
	    return false;
	}

	// Check that this link is the only one representing the missing
	// piece of the intersection curve
	double ang_tol = M_PI/4.0; // M_PI/3.0;
	if (pt1->numNeighbours() >= 2)
	{
	    vector<IntersectionPoint*> neighbours;
	    pt1->getNeighbours(neighbours);
	    Point diff_1 = pt1->getPar1Point() - pt2->getPar1Point();
	    Point diff_2 = pt1->getPar2Point() - pt2->getPar2Point();
	    for (size_t ki=0; ki<neighbours.size(); ki++)
	    {
		if (neighbours[ki] == pt2.get())
		    continue;
		Point diff2_1 = pt1->getPar1Point() - neighbours[ki]->getPar1Point();
		Point diff2_2 = pt1->getPar2Point() - neighbours[ki]->getPar2Point();
		if (diff_1.angle(diff2_1) < ang_tol || diff_2.angle(diff2_2) < ang_tol)
		{
		    // Both connections point in the same direction
		    pt1->disconnectFrom(pt2.get());
		    return false;
		}
	    }
	}

	if (pt2->numNeighbours() >= 2)
	{
	    vector<IntersectionPoint*> neighbours;
	    pt2->getNeighbours(neighbours);
	    Point diff_1 = pt2->getPar1Point() - pt1->getPar1Point();
	    Point diff_2 = pt2->getPar2Point() - pt1->getPar2Point();
	    for (size_t ki=0; ki<neighbours.size(); ki++)
	    {
		if (neighbours[ki] == pt1.get())
		    continue;
		Point diff2_1 = pt2->getPar1Point() - neighbours[ki]->getPar1Point();
		Point diff2_2 = pt2->getPar2Point() - neighbours[ki]->getPar2Point();
		//if (diff_1.angle(diff2_1) < ang_tol || diff_2.angle(diff2_2) < ang_tol)
		if (diff_1.angle(diff2_1) < ang_tol && diff_2.angle(diff2_2) < ang_tol)
		{
		    // Both connections point in the same direction
		    pt1->disconnectFrom(pt2.get());
		    return false;
		}
	    }
	}

	// Check if the current link crosses one of the already existing ones
	shared_ptr<IntersectionLink> curr_link = pt1->getIntersectionLink(pt2.get());
	typedef set<shared_ptr<IntersectionLink> >::iterator iter;
	for (iter it = links.begin(); it != links.end(); ++it) 
	{
	    int crosses = (*it)->crosses(curr_link);
	    if (crosses == 2)
	    {
		// Illegal connection
		IntersectionPoint *q1, *q2;
		(*it)->getIntersectionPoints(q1, q2);
		if (!(q1->getSingularityType() == HIGHER_ORDER_POINT &&
		      q2->getSingularityType() == HIGHER_ORDER_POINT))
		{
		    // Dismiss special case originating from self intersection
		    pt1->disconnectFrom(pt2.get());
		    return false;
		}
	    }
		
	}
	
	// Connection OK
	return true;
    }
    else
	return false;
}

//===========================================================================
void SfSfIntersector::fixCrossingLinks()
//===========================================================================
{
    // Loop through all pairs of links in the pool. Check if a pair is
    // crossing. If they are, disconnect them and try to reconnect.

    // Defining an std::set with unique links. Note that we are
    // implicitly "ordering" links using operator<() in shared_ptr,
    // since this supplies us with a strict weak ordering.
    set<shared_ptr<IntersectionLink> > links;
    int_results_->getAllLinks(links);

    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    // Loop through pairs
    IntersectionPoint *p1, *p2, *q1, *q2;
    vector<shared_ptr<IntersectionPoint> > points_to_remove;
    typedef set<shared_ptr<IntersectionLink> >::iterator iter;
    for (iter it = links.begin(); it != links.end(); ++it) {
	iter begin = it;
	++begin;
	for (iter jt = begin; jt != links.end(); ++jt) {
	    int crosses = (*it)->crosses(*jt);
	    if (crosses == 2) 
	    {
		// Links are crossing. Try to disconnect.
		(*it)->getIntersectionPoints(p1, p2);
		(*jt)->getIntersectionPoints(q1, q2);

		// Make sure that all endpoints lie in the current
		// surface-surface intersection
		if (!int_results_->isInDomain(p1) ||
		    !int_results_->isInDomain(p2) ||
		    !int_results_->isInDomain(q1) ||
		    !int_results_->isInDomain(q2))
		    continue;

		if ((p1->getSingularityType() == HIGHER_ORDER_POINT &&
		     p2->getSingularityType() == HIGHER_ORDER_POINT) ||
		    (q1->getSingularityType() == HIGHER_ORDER_POINT &&
		     q2->getSingularityType() == HIGHER_ORDER_POINT))
		    continue;  // Special case related to self intersection

		// Check if both links is neighbours to the same link. In that case
		// find a point in the interior of this link and connect both outer
		// endpoints of the current links to this new point.
		bool found = false;
		if (p1->isConnectedTo(q1) || p1->isConnectedTo(q2) ||
		    p2->isConnectedTo(q1) || p2->isConnectedTo(q2))
		{
		    IntersectionPoint *inner1, *inner2, *outer1, *outer2;
		    if (p1->isConnectedTo(q1))
		    {
			inner1 = p1;
			inner2 = q1;
			outer1 = p2;
			outer2 = q2;
		    }
		    else if (p1->isConnectedTo(q2))
		    {
			inner1 = p1;
			inner2 = q2;
			outer1 = p2;
			outer2 = q1;
		    }
		    else if (p2->isConnectedTo(q1))
		    {
			inner1 = p2;
			inner2 = q1;
			outer1 = p1;
			outer2 = q2;
		    }
		    else 
		    {
			inner1 = p2;
			inner2 = q2;
			outer1 = p1;
			outer2 = q1;
		    }

		    double param[4];
		    double dist;
		    size_t ki;
		    shared_ptr<IntersectionPoint> pnt1, pnt2;
		    vector<shared_ptr<IntersectionPoint> > int_pts;
		    int_results_->getIntersectionPoints(int_pts);
		    for (ki=0; ki<int_pts.size(); ki++)
			if (int_pts[ki].get() == inner1)
			    break;
		    if (ki < int_pts.size())
			pnt1 = int_pts[ki];
		    for (ki=0; ki<int_pts.size(); ki++)
			if (int_pts[ki].get() == inner2)
			    break;
		    if (ki < int_pts.size())
			pnt2 = int_pts[ki];

		    if (pnt1.get() && pnt2.get())
			found = findMiddlePoint(pnt1, pnt2, param, dist);
		    else 
			found = false;
		    if (found)
		    {
			shared_ptr<IntersectionPoint> mid_int =
			    int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1],
							       epsge_, param, param+2);
			inner1->disconnectFrom(inner2);
			int_results_->removeIntPoint(pnt1);
			int_results_->removeIntPoint(pnt2);
			mid_int->connectTo(outer1, REPAIRED_MISSING_LINK);
			mid_int->connectTo(outer2, REPAIRED_MISSING_LINK);
		    }

		}
		else
		{
		    // Check also the neighbouring links. Maybe the links would be OK if one
		    // of the endpoints are removed.
		    size_t ki, kj, kr;
		    vector<shared_ptr<IntersectionLink> > links2;
		    p1->getNeighbourLinks(links2);
		    for (kr=0; kr<links2.size(); kr++)
		    {
			if (links2[kr].get() == (*it).get())
			    continue;
			int crosses2 = (*jt)->crosses(links2[kr]);
			if (crosses2 == 2)
			{
			    // Two links from the point p1 crosses the linke
			    // between q1 and q2. Remove p1.
			    for (ki=0; ki<int_pts.size(); ki++)
				if (int_pts[ki].get() == p1)
				    break;
			    for (kj=0; kj<points_to_remove.size(); kj++)
				if (int_pts[ki].get() == points_to_remove[kj].get())
				    break;
			    if (kj >= points_to_remove.size())
				points_to_remove.push_back(int_pts[ki]);
			    found = true;
			}
		    }

		    if (!found)
		    {
			links2.clear();
			p2->getNeighbourLinks(links2);
			for (kr=0; kr<links2.size(); kr++)
			{
			    if (links2[kr].get() == (*it).get())
				continue;
			    int crosses2 = (*jt)->crosses(links2[kr]);
			    if (crosses2 == 2)
			    {
				// Two links from the point p2 crosses the linke
				// between q1 and q2. Remove p2.
				for (ki=0; ki<int_pts.size(); ki++)
				    if (int_pts[ki].get() == p2)
					break;
				for (kj=0; kj<points_to_remove.size(); kj++)
				    if (int_pts[ki].get() == points_to_remove[kj].get())
					break;
				if (kj >= points_to_remove.size())
				    points_to_remove.push_back(int_pts[ki]);
				found = true;
			    }
			}
		    }

		    if (!found)
		    {
			links2.clear();
			q1->getNeighbourLinks(links2);
			for (kr=0; kr<links2.size(); kr++)
			{
			    if (links2[kr].get() == (*jt).get())
				continue;
			    int crosses2 = (*it)->crosses(links2[kr]);
			    if (crosses2 == 2)
			    {
				// Two links from the point q1 crosses the linke
				// between p1 and p2. Remove q1.
				for (ki=0; ki<int_pts.size(); ki++)
				    if (int_pts[ki].get() == q1)
					break;
				for (kj=0; kj<points_to_remove.size(); kj++)
				    if (int_pts[ki].get() == points_to_remove[kj].get())
					break;
				if (kj >= points_to_remove.size())
				    points_to_remove.push_back(int_pts[ki]);
				found = true;
			    }
			}
		    }
		    if (!found)
		    {
			links2.clear();
			q2->getNeighbourLinks(links2);
			for (kr=0; kr<links2.size(); kr++)
			{
			    if (links2[kr].get() == (*jt).get())
				continue;
			    int crosses2 = (*it)->crosses(links2[kr]);
			    if (crosses2 == 2)
			    {
				// Two links from the point q2 crosses the linke
				// between p1 and p2. Remove q2.
				for (ki=0; ki<int_pts.size(); ki++)
				    if (int_pts[ki].get() == q2)
					break;
				for (kj=0; kj<points_to_remove.size(); kj++)
				    if (int_pts[ki].get() == points_to_remove[kj].get())
					break;
				if (kj >= points_to_remove.size())
				    points_to_remove.push_back(int_pts[ki]);
				found = true;
			    }
			}
		    }
		}

		if (!found)
		{
		    // Try to refine the longest link and check if the links still cross
		    shared_ptr<IntersectionPoint> pnt1, pnt2;
		    vector<shared_ptr<IntersectionPoint> > int_pts;
		    int_results_->getIntersectionPoints(int_pts);
		    size_t ki;
		    double param[4];
		    double dist;
		    if (true)
			//(p1->getPoint()).dist(p2->getPoint()) > (q1->getPoint()).dist(q2->getPoint()))
		    {
			for (ki=0; ki<int_pts.size(); ki++)
			    if (int_pts[ki].get() == p1)
				break;
			if (ki < int_pts.size())
			    pnt1 = int_pts[ki];
			for (ki=0; ki<int_pts.size(); ki++)
			    if (int_pts[ki].get() == p2)
				break;
			if (ki < int_pts.size())
			    pnt2 = int_pts[ki];

			if (pnt1.get() && pnt2.get())
			    found = findMiddlePoint(pnt1, pnt2, param, dist);
			else 
			    found = false;
			if (found)
			{
			    p1->disconnectFrom(p2);
			    shared_ptr<IntersectionPoint> mid_int =
				int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1],
								   epsge_, param, param+2);
			    mid_int->connectTo(p1, REPAIRED_MISSING_LINK);
			    mid_int->connectTo(p2, REPAIRED_MISSING_LINK);
			}
		    }
		    if (true) //else
		    {
			for (ki=0; ki<int_pts.size(); ki++)
			    if (int_pts[ki].get() == q1)
				break;
			if (ki < int_pts.size())
			    pnt1 = int_pts[ki];
			for (ki=0; ki<int_pts.size(); ki++)
			    if (int_pts[ki].get() == q2)
				break;
			if (ki < int_pts.size())
			    pnt2 = int_pts[ki];

			if (pnt1.get() && pnt2.get())
			    found = findMiddlePoint(pnt1, pnt2, param, dist);
			else 
			    found = false;
			if (found)
			{
			    q1->disconnectFrom(q2);
			    shared_ptr<IntersectionPoint> mid_int =
				int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1],
								   epsge_, param, param+2);
			    mid_int->connectTo(q1, REPAIRED_MISSING_LINK);
			    mid_int->connectTo(q2, REPAIRED_MISSING_LINK);
			}
		    }
			
		}

		// At this stage, we know that all illegal links are removed.
		// Check anyway if one link is better than the other. If no such
		// decisions can be made, disconnect both and let repairMissingLinks
		// mend the holes later
		if (!found)
		{
		    if ((q1->numNeighbours() < 2 || q2->numNeighbours() < 2) &&
			p1->numNeighbours() >= 2 && p2->numNeighbours() >= 2)
			q1->disconnectFrom(q2);
		    else if ((p1->numNeighbours() < 2 || p2->numNeighbours() < 2) &&
			     q1->numNeighbours() >= 2 && q2->numNeighbours() >= 2)
			p1->disconnectFrom(p2);
		    else
		    {
			p1->disconnectFrom(p2);
			q1->disconnectFrom(q2);
		    }
		}
	    }
	    else if (crosses == 1)
	    {
		// Intersections links are not crossing, but lie close.
		// Check if they represent the same intersection
		bool remove;
		size_t ki, kj;

		(*it)->getIntersectionPoints(p1, p2);
		(*jt)->getIntersectionPoints(q1, q2);

		if ((p1->getSingularityType() == HIGHER_ORDER_POINT &&
		     p2->getSingularityType() == HIGHER_ORDER_POINT) ||
		    (q1->getSingularityType() == HIGHER_ORDER_POINT &&
		     q2->getSingularityType() == HIGHER_ORDER_POINT))
		    continue;  // Special case related to self intersection

		for (ki=0; ki<int_pts.size(); ki++)
		    if (int_pts[ki].get() == p1)
			break;
		if (ki < int_pts.size())
		{
		    remove = checkCloseEndpoint(int_pts[ki], *jt);
		    if (remove)
		    {
			for (kj=0; kj<points_to_remove.size(); kj++)
			    if (int_pts[ki].get() == points_to_remove[kj].get())
				break;
			if (kj >= points_to_remove.size())
			    points_to_remove.push_back(int_pts[ki]);
		    }
		}

		for (ki=0; ki<int_pts.size(); ki++)
		    if (int_pts[ki].get() == p2)
			break;
		if (ki < int_pts.size())
		{
		    remove = checkCloseEndpoint(int_pts[ki], *jt);
		    if (remove)
		    {
			for (kj=0; kj<points_to_remove.size(); kj++)
			    if (int_pts[ki].get() == points_to_remove[kj].get())
				break;
			if (kj >= points_to_remove.size())
			    points_to_remove.push_back(int_pts[ki]);
		    }
		}

		for (ki=0; ki<int_pts.size(); ki++)
		    if (int_pts[ki].get() == q1)
			break;
		if (ki < int_pts.size())
		{
		    remove = checkCloseEndpoint(int_pts[ki], *it);
		    if (remove)
		    {
			for (kj=0; kj<points_to_remove.size(); kj++)
			    if (int_pts[ki].get() == points_to_remove[kj].get())
				break;
			if (kj >= points_to_remove.size())
			    points_to_remove.push_back(int_pts[ki]);
		    }
		}

		for (ki=0; ki<int_pts.size(); ki++)
		    if (int_pts[ki].get() == q2)
			break;
		if (ki < int_pts.size())
		{
		    remove = checkCloseEndpoint(int_pts[ki], *it);
		    if (remove)
		    {
			for (kj=0; kj<points_to_remove.size(); kj++)
			    if (int_pts[ki].get() == points_to_remove[kj].get())
				break;
			if (kj >= points_to_remove.size())
			    points_to_remove.push_back(int_pts[ki]);
		    }
		}
	    }
	}
    }
	
    for (size_t kj=0; kj<points_to_remove.size(); kj++)
	int_results_->removeIntPoint(points_to_remove[kj]);
	
    return;
}


//===========================================================================
bool SfSfIntersector::checkCloseEndpoint(shared_ptr<IntersectionPoint> pnt, 
					 shared_ptr<IntersectionLink> link)
//===========================================================================
{
    // First check that the given point belongs to a short intersection curve
    size_t ki;
    if (int_results_->isBoundaryPoint(pnt))
	return false;

    if (pnt->numNeighbours() > 1)
    {
	// Check neighbours
	vector<IntersectionPoint*> neighbours;
	pnt->getNeighbours(neighbours);
	for (ki=0; ki<neighbours.size(); ki++)
	{
	    if (neighbours[ki]->numNeighbours() < 2 && 
		!int_results_->isBoundaryPoint(neighbours[ki]))
		break;
	}
	if (ki == neighbours.size())
	    return false;
    }

    // Check if the point lies internal to the given link by projecting the
    // given point onto the parameterized link
    IntersectionPoint *q1, *q2;
    link->getIntersectionPoints(q1, q2);
    Point a = q1->getPoint();
    Point b = q2->getPoint();
    Point c = pnt->getPoint();
    double div = (b-a)*(b-a);
    double s = (a-c)*(a-b)/div;
    if (s < 0.0 || s > 1.0)
	return false;  // Not an internal point

    // Check consistency of tangent directions in intersection points
    bool ok_dir1, ok_dir2;
    shared_ptr<IntersectionLink> ln1
	= q1->connectTo(pnt, REPAIRED_MISSING_LINK);
    ok_dir1 = int_results_->verifyIntersectionLink(ln1);
    pnt->disconnectFrom(q1);
    shared_ptr<IntersectionLink> ln2
	= q2->connectTo(pnt, REPAIRED_MISSING_LINK);
    ok_dir2 = int_results_->verifyIntersectionLink(ln2);
    pnt->disconnectFrom(q2);
    if (!(ok_dir1 && ok_dir2))
	return false;

    // Refine the intersection link at the closest point to the given intersection
    // point
    double epsge = pnt->getTolerance()->getEpsge();
    Point mid = (1.0 - s)*a + s*b;
    Point tan = b - a;
    tan.normalize();
    Point mid_par1 = (1.0 - s)*q1->getPar1Point() + s*q2->getPar1Point();
    Point mid_par2 = (1.0 - s)*q1->getPar2Point() + s*q2->getPar2Point();
    
    vector<Point> plane(2);  // Plane at midpoint of link, perpendicular to link
    plane[0] = mid;
    plane[1] = tan;

    vector<Point> input_point_1(7, Point(3)); // only entry [0], [1], [2] and [6] will be used
    vector<Point> input_point_2(7, Point(3)); // only entry [1], [1], [2] and [6] will be used
    vector<Point> result_pt_1(7, Point(3));
    vector<Point> result_pt_2(7, Point(3));
    Point surface_1_param, surface_2_param;
    
    shared_ptr<ParamSurface> s1 = obj_int_[0]->getParamSurfaceInt()->getParamSurface();
    shared_ptr<ParamSurface> s2 = obj_int_[1]->getParamSurfaceInt()->getParamSurface();
    s1->point(input_point_1, mid_par1[0], mid_par1[1], 1);
    s2->point(input_point_2, mid_par2[0], mid_par2[1], 1);
    s1->normal(input_point_1[6], mid_par1[0], mid_par1[1]);
    s2->normal(input_point_2[6], mid_par2[0], mid_par2[1]);
    int jstat=0;
    ClosestPoint::closestPtSurfSurfPlane(plane, input_point_1, input_point_2, mid_par1, mid_par2,
			   s1.get(), s2.get(), epsge, result_pt_1, result_pt_2, 
			   surface_1_param, surface_2_param, jstat);
    if (jstat != 1)
	return false;  // Iteration to point on intersection curve did not succeed

    shared_ptr<IntersectionPoint> tmp = 
	int_results_->addIntersectionPoint(obj_int_[0], obj_int_[1], pnt->getTolerance(),
					   surface_1_param.begin(), surface_2_param.begin());

    // Check if the closest point on the link represent the same intersection
    // as the input point
    bool remove = canConnect(tmp, pnt);
    int_results_->removeIntPoint(tmp);
    
    return remove;
}
	

//===========================================================================
void SfSfIntersector::iterateOnIntersectionPoints()
//===========================================================================
{
    // Purpose: We try all possible iterations of each intersection
    // point to see if we can improve their position. Results are
    // written to two dsp-files. @@@jbt - Experimental!

//     ofstream out1("iterations1.dsp");
//     ofstream out2("iterations2.dsp");
//     out1 << "pnttype: cross" << endl
// 	 << "fg: grey" << endl;
//     out2 << "pnttype: cross" << endl
// 	 << "fg: grey" << endl;

//     cout << "*** Iterating points ***" << endl;
    vector<shared_ptr<IntersectionPoint> > int_pts
	= int_results_->getIntersectionPoints();
    int int_sz = int(int_pts.size());
    const double parfac = 0.1; // 0.01;
    const double epsdist = 1.0e-12;
    for (int i = 0; i < int_sz; ++i) 
    {
	bool bd_pt[4];
	bd_pt[0] = bd_pt[1] = bd_pt[2] = bd_pt[3] = false;
	if (int_results_->isBoundaryPoint(int_pts[i]))
	{
	    // Take care upon how boundary points are iterated
	    int k1, k2, kp = 0;
	    for (k1=0, kp=0; k1<2; k1++)
		for (k2=0; k2<2; k2++, kp++)
		{
		    double start = obj_int_[k1]->startParam(k2);
		    double end = obj_int_[k1]->endParam(k2);
		    double par = int_pts[i]->getPar(kp);
		    if (fabs(par-start)<epsdist || fabs(end-par)<epsdist)
			bd_pt[kp] = true;
		}
	}
	double dist = int_pts[i]->getDist();
	if (dist < 0.1*epsdist)
	    continue;  // Don't bother
// 	cout << i << "    d = " << dist << endl;
	Point point = int_pts[i]->getPoint();
	double minmovedist = numeric_limits<double>::max();
	double minparam[4];
	double minseed[4];
	int movedir = 0;
	bool replace = false;
	vector<double> seed = int_pts[i]->getPar();
	for (int j = 0; j < 4; ++j) {
	    double param[4];
	    double newdist;
	    vector<double> limit1 = seed;
	    vector<double> limit2 = seed;
	    for (int k = 0; k < 4; ++k) {
		double start = int_results_->startParam(k);
		double end = int_results_->endParam(k);
		double delta = parfac * (start - end);
		if (!bd_pt[k])
		{
		    limit1[k] -= delta;
		    limit2[k] += delta;
		}
		if (limit1[k] < start) {
		    limit1[k] = start;
		}
		if (limit2[k] > end) {
		    limit2[k] = end;
		}
	    }
	    double dummy = 0.0;
// 	    cout << "\t" << "dir = " << j << " -> ";
	    try {
		doIterate(j, param, &limit1[0], dummy, &limit2[0], dummy,
			  newdist, &seed[0]);
	    }
	    catch (...) {
// 		cout << "Failed" << endl;
		continue;
	    }
	    // If we iterated out of the domain, we must move the point
	    bool outside = false;
	    for (int k = 0; k < 4; ++k) {
		double start = int_results_->startParam(k);
		double end = int_results_->endParam(k);
		if (param[k] < start) {
		    param[k] = start;
		    outside = true;
		}
		if (param[k] > end) {
		    param[k] = end;
		    outside = true;
		}
	    }
	    Point newpoint1, newpoint2;
	    obj_int_[0]->getParamSurfaceInt()->point(newpoint1, param);
	    obj_int_[1]->getParamSurfaceInt()->point(newpoint2, param+2);
	    if (outside) {
		Point newdiff = newpoint1 - newpoint2;
		newdist = newdiff.length();
	    }
// 	    cout << "d = " << newdist << endl;
	    if (!(newdist < dist && newdist < epsdist))
		continue;
	    // We got here - we can replace the point
	    replace = true;
	    Point diff1 = newpoint1 - point;
	    Point diff2 = newpoint2 - point;
	    double movedist = max(diff1.length(), diff2.length());
	    if (minmovedist > movedist) {
		minmovedist = movedist;
		movedir = j;
		for(int k = 0; k < 4; ++k) {
		    minseed[k] = seed[k];
		    minparam[k] = param[k];
		}
	    }
	}
	if (replace) {
	    int_pts[i]->replaceParameter(minparam);
	    cout <<  "Init par: " << seed[0] << " " << seed[1] << " ";
	    cout << seed[2] << " " << seed[3] << ", " << dist << " " << i << endl;
	    cout <<  "Moved to: " << minparam[0] << " " << minparam[1] << " ";
	    cout << minparam[2] << " " << minparam[3] << ", " << int_pts[i]->getDist() << endl;

// 	    cout << "\t" << "movedir = " << movedir << endl;

// 	    out1 << "lin:" << endl
// 		 << minseed[0]<< " " << minseed[1] << endl
// 		 << minparam[0] << " " << minparam[1] << endl
// 		 << "pnt:" << endl
// 		 << minseed[0] << " " << minseed[1] << endl;
// 	    out2 << "lin:" << endl
// 		 << minseed[2]<< " " << minseed[3] << endl
// 		 << minparam[2] << " " << minparam[3] << endl
// 		 << "pnt:" << endl
// 		 << minseed[2] << " " << minseed[3] << endl;
	}
    }

    return;
}


//===========================================================================
void SfSfIntersector::getSingularityBox(shared_ptr<IntersectionPoint> sing,
					double frompar[], double topar[])
//===========================================================================
{
    // Get a box around the singularity, i.e. branch point. The tricky
    // part is getting the size of the box. We solve this by analyzing
    // the curvature in the point. More precisely, we look at Dupin's
    // indicatrices for the two surfaces, and extract two lengths from
    // each indicatrix. The indicatrices are scaled such that it
    // corresponds to a plane parallel to the tangent plane at the
    // distance equal to the tolerance.

    shared_ptr<ParamSurfaceInt> sf1
	= dynamic_pointer_cast<ParamSurfaceInt, ParamGeomInt>(obj_int_[0]);
    shared_ptr<ParamSurfaceInt> sf2
	= dynamic_pointer_cast<ParamSurfaceInt, ParamGeomInt>(obj_int_[1]);
    ASSERT(sf1.get() != 0 && sf2.get() != 0);

    double tol = epsge_->getEpsge();
    double res = epsge_->getRelParRes();
    vector<double> singpar = sing->getPar();
    bool u_from_right, v_from_right;
    bool is_plane;
    double eps = numeric_limits<double>::min();

    // Get curvature info
    double E, F, G, L, M, N; // First and second fundamental forms
    double K, H; // Gaussian and mean curvatures
    double k1, k2; // Principal curvatures
    double d1, d2; // "Scaled sqrt(rho)'s" in Dupin's indicatrix
    double det;    // Determinant
    double mindist1 = 0.0;
    double mindist2 = 0.0;

   // First surface
    if (singpar[0] < int_results_->startParam(0) + res)
	u_from_right = true;
    else
	u_from_right = false;
    if (singpar[1] < int_results_->startParam(1) + res)
	v_from_right = true;
    else
	v_from_right = false;
    sf1->first_fund_form(singpar[0], singpar[1], u_from_right, v_from_right,
			 E, F, G);
    sf1->second_fund_form(singpar[0], singpar[1], u_from_right, v_from_right,
			  L, M, N);
    det = E * G - F * F;
    if (fabs(det) < res);
    
    else
    {
	K = (L * N - M * M) / det;
	H = 0.5 * (N * E - 2.0 * M * F + L * G) / det;
	k1 = H + sqrt(H * H - K);
	k2 = H - sqrt(H * H - K);
	is_plane = (k1 < eps && k2 < eps) ? true : false;
	if (!is_plane) {
	    d1 = sqrt(2.0 * tol / fabs(k1));
	    d2 = sqrt(2.0 * tol / fabs(k2));
	    mindist1 = (d1 < d2 ? d1 : d2);
	}
	// Second surface
	if (singpar[2] < int_results_->startParam(2) + res)
	    u_from_right = true;
	else
	    u_from_right = false;
	if (singpar[3] < int_results_->startParam(3) + res)
	    v_from_right = true;
	else
	    v_from_right = false;
	sf2->first_fund_form(singpar[2], singpar[3], u_from_right, v_from_right,
			     E, F, G);
	sf2->second_fund_form(singpar[2], singpar[3], u_from_right, v_from_right,
			      L, M, N);
	K = (L * N - M * M) / (E * G - F * F);
	H = 0.5 * (N * E - 2.0 * M * F + L * G) / (E * G - F * F);
	k1 = H + sqrt(H * H - K);
	k2 = H - sqrt(H * H - K);
	d1 = sqrt(2.0 * tol / fabs(k1));
	d2 = sqrt(2.0 * tol / fabs(k2));
	is_plane = (k1 < eps && k2 < eps) ? true : false;
	if (!is_plane) {
	    d1 = sqrt(2.0 * tol / fabs(k1));
	    d2 = sqrt(2.0 * tol / fabs(k2));
	    mindist2 = (d1 < d2 ? d1 : d2);
	}
    }

    // What are the distances in parameter space? We are not using the
    // euclidean distance in parameter space, but a scaled version:
    // ds^2 = a0^2*du0^2 + b0^2*du1^2 + a1^2*du1^2 + b1^2*dv1^2, where
    // (du0, dv0, du1, dv1) is a vector in 4-dimensional parameter
    // space, and a0 = |dS0/du0|, b0 = |dS0/dv0|, a1 = |dS1/du1|, and
    // b1 = |dS1/dv1|. In this way we are using information about the
    // tangent planes.
    vector<Point> pt1(3);
    vector<Point> pt2(3);
    obj_int_[0]->point(pt1, sing->getPar1(), 1);
    obj_int_[1]->point(pt2, sing->getPar2(), 1);
    vector<double> g(4);
    g[0] = pt1[1].length();
    g[1] = pt1[2].length();
    g[2] = pt2[1].length();
    g[3] = pt2[2].length();

    // We use the greatest of mindist1 and mindist2 as the minimum
    // distance. We can now set the box around the singularity.
    double dist = (mindist1 > mindist2 ? mindist1 : mindist2);
    for (int i = 0; i < 4; ++i) {
	double start = int_results_->startParam(i);
	double end = int_results_->endParam(i);
	double del = (fabs(g[i]) > res) ? 2.0 * dist / g[i] : 100.0*res; // Including factor 2.0
	frompar[i] = singpar[i] - del;
	if (frompar[i] < start) {
	    frompar[i] = start;
	}
	topar[i] = singpar[i] + del;
	if (topar[i] > end) {
	    topar[i] = end;
	}
	cout << i << ":  " << frompar[i] << " - " << topar[i] << endl;
    }

    return;


//     // Get a box around the singularity, i.e. branch point. The tricky
//     // part is getting the size of the box. We solve this by taking
//     // the following steps: 1) Find the normal to the surfaces in the
//     // branch point (this is common to both surfaces), 2) Set up and
//     // compute a new intersection problem by translating the second
//     // surface along the direction of the normal with the distance
//     // equal to the tolerance, 3) Set up and compute another problem
//     // by translating the second surface in the opposite direction, 4)
//     // Identify the intersection points among the new points that are
//     // closest to the original branch point, 5) Finally, get the
//     // distances from these points to the branch point in the
//     // parameter plane.

//     shared_ptr<ParamSurfaceInt> sf1
// 	= dynamic_pointer_cast<ParamSurfaceInt, ParamGeomInt>(obj_int_[0]);
//     shared_ptr<ParamSurfaceInt> sf2
// 	= dynamic_pointer_cast<ParamSurfaceInt, ParamGeomInt>(obj_int_[1]);
//     ASSERT(sf1.get() != 0 && sf2.get() != 0);

//     // Get the normal
//     Point normal1, normal2;
//     sf1->getParamSurface()->normal(normal1, sing->getPar(0), sing->getPar(1));
//     sf2->getParamSurface()->normal(normal2, sing->getPar(2), sing->getPar(3));
//     Point normal = normal1 + normal2;
//     normal.normalize();
//     double tol = epsge_->getEpsge();
//     double fac = 2.0; // What factor should we use? 2.0 seems
// 		      // OK. @@@jbt.
//     normal *= fac * tol;

//     // Get translated surfaces
//     shared_ptr<SplineSurfaceInt> ssi
// 	= dynamic_pointer_cast<SplineSurfaceInt, ParamSurfaceInt>(sf2);
//     shared_ptr<SplineSurface> spline = ssi->getSplineSurface();
//     ASSERT(spline.get() != 0);
//     shared_ptr<SplineSurface> splinetmp(spline->clone());
//     shared_ptr<SplineSurface> splineup(spline->clone());
//     shared_ptr<SplineSurface> splinedown(spline->clone());
//     int dim = spline->dimension();
//     ASSERT(dim == normal.dimension());
//     bool rational = spline->rational();
//     int ncoefsu = spline->numCoefs_u();
//     int ncoefsv = spline->numCoefs_v();
//     vector<double>::iterator coefsup, coefsdown;
//     if (!rational) {
// 	coefsup = splineup->coefs_begin();
// 	coefsdown = splinedown->coefs_begin();
//     }
//     else {
// 	coefsup = splineup->rcoefs_begin();
// 	coefsdown = splinedown->rcoefs_begin();

//     }
//     int kdim = (rational ? dim+1 : dim);
//     for (int j = 0; j < ncoefsv; ++j) {
// 	for (int i = 0; i < ncoefsu; ++i) {
// 	    if (!rational) {
// 		for (int k = 0; k < dim; ++k) {
// 		    *(coefsup+k) += normal[k];
// 		    *(coefsdown+k) -= normal[k];
// 		}
// 	    }
// 	    else {
// 		double w = *(coefsup+dim);
// 		for (int k = 0; k < dim; ++k) {
// 		    *(coefsup+k) += (w * normal[k]);
// 		    *(coefsdown+k) -= (w * normal[k]);
// 		}
// 	    }
// 	    coefsup += kdim;
// 	    coefsdown += kdim;
// 	}
//     }

// //     exit(0); // Stop before we translate

//     // Set up and run the two intersection problems
//     shared_ptr<ParamSurface> ps2up(splineup->clone());
//     shared_ptr<ParamGeomInt> obj2up(new SplineSurfaceInt(ps2up));
//     shared_ptr<SfSfIntersector> prev_dummy_up
// 	= shared_ptr<SfSfIntersector>
// 	(new SfSfIntersector(obj_int_[0], obj2up, epsge_));
//     shared_ptr<SfSfIntersector> intersectorup
// 	= shared_ptr<SfSfIntersector>
// 	(new SfSfIntersector(obj_int_[0], obj2up, epsge_,
// 			     prev_dummy_up.get()));
//     intersectorup->compute();

// //     exit(0); // Stop after upward translation

//     shared_ptr<ParamSurface> ps2down(splinedown->clone());
//     shared_ptr<ParamGeomInt> obj2down(new SplineSurfaceInt(ps2down));
//     shared_ptr<SfSfIntersector> prev_dummy_down
// 	= shared_ptr<SfSfIntersector>
// 	(new SfSfIntersector(obj_int_[0], obj2down, epsge_));
//     shared_ptr<SfSfIntersector> intersectordown
// 	= shared_ptr<SfSfIntersector>
// 	(new SfSfIntersector(obj_int_[0], obj2down, epsge_,
// 			     prev_dummy_down.get()));
//     intersectordown->compute();

// //     exit(0); // Stop after downward translation

//     // What are the closest points for the up and down cases? We are
//     // not using the euclidean distance in parameter space, but a
//     // scaled version: ds^2 = a0^2*du0^2 + b0^2*du1^2 + a1^2*du1^2 +
//     // b1^2*dv1^2, where (du0, dv0, du1, dv1) is a vector in
//     // 4-dimensional parameter space, and a0 = |dS0/du0|, b0 =
//     // |dS0/dv0|, a1 = |dS1/du1|, and b1 = |dS1/dv1|. In this way
//     // we are using information about the tangent planes.
//     vector<Point> pt1(3);
//     vector<Point> pt2(3);
//     obj_int_[0]->point(pt1, sing->getPar1(), 1);
//     obj_int_[1]->point(pt2, sing->getPar2(), 1);
//     vector<double> g(4);
//     g[0] = pt1[1].length();
//     g[1] = pt1[2].length();
//     g[2] = pt2[1].length();
//     g[3] = pt2[2].length();
//     vector<shared_ptr<IntersectionPoint> > intup
// 	= intersectorup->getIntPool()->getIntersectionPoints();
//     int nintup = intup.size();
//     shared_ptr<IntersectionPoint> closestup;
//     double mindist2up = numeric_limits<double>::max();
//     vector<double> singpar = sing->getPar();
//     vector<double> intpar;
//     for (int i = 0; i < nintup; ++i) {
// 	intpar = intup[i]->getPar();
// 	double dist2 = 0.0;
// 	for (int j = 0; j < 4; ++j) {
// 	    double d = g[j] * (singpar[j] - intpar[j]);
// 	    dist2 += d * d;
// 	}
// 	if (dist2 < mindist2up) {
// 	    mindist2up = dist2;
// 	    closestup = intup[i];
// 	}
//     }
//     vector<shared_ptr<IntersectionPoint> > intdown
// 	= intersectordown->getIntPool()->getIntersectionPoints();
//     int nintdown = intdown.size();
//     shared_ptr<IntersectionPoint> closestdown;
//     double mindist2down = numeric_limits<double>::max();
//     for (int i = 0; i < nintdown; ++i) {
// 	intpar = intdown[i]->getPar();
// 	double dist2 = 0.0;
// 	for (int j = 0; j < 4; ++j) {
// 	    double d = g[j] * (singpar[j] - intpar[j]);
// 	    dist2 += d * d;
// 	}
// 	if (dist2 < mindist2down) {
// 	    mindist2down = dist2;
// 	    closestdown = intdown[i];
// 	}
//     }

//     // We use the _farthest_ of closestup and closestdown as the
//     // closest point. We can now set the box around the singularity.
//     double mindist2
// 	= (mindist2up > mindist2down ? mindist2up : mindist2down);
//     double dist = sqrt(mindist2);
//     cout << "dist = " << dist << endl;
//     for (int i = 0; i < 4; ++i) {
// 	double start = int_results_->startParam(i);
// 	double end = int_results_->endParam(i);
// 	double del = dist / g[i];
// 	frompar[i] = singpar[i] - del;
// 	if (frompar[i] < start) {
// 	    frompar[i] = start;
// 	}
// 	topar[i] = singpar[i] + del;
// 	if (topar[i] > end) {
// 	    topar[i] = end;
// 	}
// 	cout << i << ":  " << frompar[i] << " - " << topar[i] << endl;
//     }

//     return;
}


//===========================================================================
shared_ptr<SfSfIntersector> SfSfIntersector::
getSubIntersector(double frompar[], double topar[])
//===========================================================================
{
    double fuzzy = epsge_->getRelParRes();
    shared_ptr<ParamSurfaceInt> obj1
	= dynamic_pointer_cast<ParamSurfaceInt, ParamGeomInt>
	(obj_int_[0]);
    vector<shared_ptr<ParamSurfaceInt> > sub1
	= obj1->subSurfaces(frompar[0], frompar[1],
			    topar[0], topar[1], fuzzy);
    shared_ptr<ParamSurfaceInt> obj2
	= dynamic_pointer_cast<ParamSurfaceInt, ParamGeomInt>
	(obj_int_[1]);
    vector<shared_ptr<ParamSurfaceInt> > sub2
	= obj2->subSurfaces(frompar[2], frompar[3],
			    topar[2], topar[3], fuzzy);
    // Create intersector
    return shared_ptr<SfSfIntersector>
	(new SfSfIntersector(sub1[0], sub2[0], epsge_));
}


//===========================================================================
bool SfSfIntersector::
isConnected(vector<pair<shared_ptr<IntersectionPoint>,
	    IntPtClassification> >& bd_ints,
	    int nmbbd)
//===========================================================================
{
    // Purpose: Check if all the points in a set of intersection
    // points at the surface boundaries is connected to some other
    // point in the set

    vector<bool> found_connection(nmbbd, false);
    int ki, kj;
    for (ki=0; ki<nmbbd; ki++)
    {
	for (kj=ki+1; kj<nmbbd; kj++)
	{
	    if (found_connection[ki] && found_connection[kj])
		continue;
	    if (bd_ints[ki].first == bd_ints[kj].first)
		continue;  // The same branch point

	    if (bd_ints[ki].first->isConnectedTo(bd_ints[kj].first))
		found_connection[ki] = found_connection[kj] = true;
	}
    }
    for (ki=0; ki<nmbbd; ki++)
	if (!found_connection[ki])
	    return false;

    return true;
}


//===========================================================================
bool SfSfIntersector::
isConnected(vector<shared_ptr<IntersectionPoint> >& bd_ints)
//===========================================================================
{
    // Purpose: Check if all the points in a set of intersection
    // points at the surface boundaries is connected to some other
    // point in the set

    int nmbbd = (int)bd_ints.size();
    vector<bool> found_connection(nmbbd, false);
    int ki, kj;
    for (ki=0; ki<nmbbd; ki++)
    {
	for (kj=ki+1; kj<nmbbd; kj++)
	{
	    if (found_connection[ki] && found_connection[kj])
		continue;
	    if (bd_ints[ki]->isConnectedTo(bd_ints[kj]))
		found_connection[ki] = found_connection[kj] = true;
	}
    }
    for (ki=0; ki<nmbbd; ki++)
	if (!found_connection[ki])
	    return false;

    return true;
}


//===========================================================================
bool SfSfIntersector::
connectDirected(vector<pair<shared_ptr<IntersectionPoint>,
		IntPtClassification> >& bd_ints, int nmbbd)
//===========================================================================
{
    // Purpose: Connect boundary intersection in correct sequence if
    // their type makes it possible

    if (nmbbd == 1)
	return false;  // No rule

    // @@@ VSK. Should we test on distance or influence area. Connect
    // anyway if the points are close enough.

    // First test if connection is possible
    int ki;
    double max_ang = 0.2;  
    for (ki = 1; ki < nmbbd; ki++) {
	if (bd_ints[ki-1].first == bd_ints[ki].first) {
	    // Same branch point
	} else if (bd_ints[ki-1].first->isConnectedTo(bd_ints[ki].first)) {
	    // Already connected. Pass to the next point
	} else if (bd_ints[ki-1].second == DIR_OUT) {
	    // The intersection curve cannot start outwards
	    return false;
	} else if (bd_ints[ki-1].second == DIR_IN &&
		   bd_ints[ki].second == DIR_IN) {
	    // Inconsistent directions
	    return false;
	} else if ((bd_ints[ki-1].second == DIR_IN
		    || bd_ints[ki-1].second == DIR_PERPENDICULAR)
		   && (bd_ints[ki].second == DIR_OUT
		       || bd_ints[ki].second == DIR_PERPENDICULAR)) {
	    // Check next point
	    if (ki<nmbbd-1 && bd_ints[ki+1].second == DIR_OUT)
	    {
		return false;  // Something is wrong. Maybe some other
		// connection method is better
	    }
	    // OK
	    ki++;
	} else if ((bd_ints[ki-1].second == DIR_IN
		    || bd_ints[ki-1].second == DIR_PERPENDICULAR) &&
		   bd_ints[ki].second == DIR_PARALLEL) {
	    // OK @@@jbt: Hmmm...?
	} else if (bd_ints[ki-1].second == DIR_PARALLEL
		   && (bd_ints[ki].second == DIR_OUT
		       || bd_ints[ki].second == DIR_PERPENDICULAR)) {
	    // @@@jbt: OK?
	} 
	else if (bd_ints[ki-1].second == DIR_PARALLEL &&
		 bd_ints[ki].second == DIR_PARALLEL)
	{
	    // VSK, 0611. These cases occur and they may or may not be
	    // OK. Make an extra check. First check consistence in 
	    // the tangent of the intersection curve
	    Point tang1 = bd_ints[ki-1].first->getTangent();
	    Point tang2 = bd_ints[ki].first->getTangent();
	    if (tang1.angle(tang2) > max_ang)
		return false;

	    // Check if an intersection point can be inserted between
	    // the current ones.
	    bool isOK = canConnect(bd_ints[ki-1].first, bd_ints[ki].first);
	    if (!isOK)
		return false;
	}
	else if (bd_ints[ki-1].second == DIR_HIGHLY_SINGULAR ||
		 bd_ints[ki].second == DIR_HIGHLY_SINGULAR)
	{
	    // VSK, 0611. These cases probably originate from self
	    // intersection. It is a simple case so check if 
	    // a connection is possible
	    // A tangency check does not make sense as the tangents cannot
	    // be computed, but check if an intersection point can be inserted 
	    // between the current ones.
	    bool isOK = canConnect(bd_ints[ki-1].first, bd_ints[ki].first);
	    if (!isOK)
		return false;
	}
	else if (bd_ints[ki-1].second == DIR_PARALLEL &&
		 bd_ints[ki].second == DIR_IN)
	{
	    // VSK, 0705. After all it is a simple case so its nothing wrong in the
	    // configuration. Do not connect, but accept
	    ;
	}
	else {
	    // No rule
	    return false;
	}
    }
    if (ki == nmbbd-1 && (bd_ints[ki].second == DIR_IN
			  || bd_ints[ki].second == DIR_OUT
			  || bd_ints[ki].second == DIR_PERPENDICULAR)) {
	// Hanging boundary point not possible to connect
	return false;
    }

    // Perform connections
    bool has_connected = false;
    for (ki = 1; ki < nmbbd; ki++) {
	if (bd_ints[ki-1].first == bd_ints[ki].first) {
	    // Same branch point
	} else if (bd_ints[ki-1].first->isConnectedTo(bd_ints[ki].first)) {
	    // Already connected. Pass to the next point
	} else if ((bd_ints[ki-1].second == DIR_IN
		    || bd_ints[ki-1].second == DIR_PERPENDICULAR)
		   && (bd_ints[ki].second == DIR_OUT
		       || bd_ints[ki].second == DIR_PERPENDICULAR)) {
	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first,
						 SIMPLE_CONE);
	    has_connected = true;
	    ki++; 
	} else if ((bd_ints[ki-1].second == DIR_IN
		  || bd_ints[ki-1].second == DIR_PERPENDICULAR)
		 && (bd_ints[ki].second == DIR_PARALLEL)) {
	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first,
						 SIMPLE_CONE);
	    has_connected = true;
	} else if ((bd_ints[ki-1].second == DIR_PARALLEL)
		 && (bd_ints[ki].second == DIR_OUT
		     || bd_ints[ki].second == DIR_PERPENDICULAR)) {
	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first,
						 SIMPLE_CONE);
	    has_connected = true;
	}
	else if (bd_ints[ki-1].second == DIR_PARALLEL &&
		 bd_ints[ki].second == DIR_PARALLEL)
	{
	    // VSK, 0611. This case has been subject to an extra check. 
	    // OK. Connect
	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first,
						 SIMPLE_CONE);
	    has_connected = true;
	}
	else if (bd_ints[ki-1].second == DIR_HIGHLY_SINGULAR ||
		 bd_ints[ki].second == DIR_HIGHLY_SINGULAR)
	{
	    // VSK, 0611. This case has been subject to an extra check. 
	    // OK. Connect
	    bd_ints[ki-1].first.get()->connectTo(bd_ints[ki].first,
						 SIMPLE_CONE);
	    has_connected = true;
	}
    }
    return has_connected;
}


//===========================================================================
bool SfSfIntersector::canConnect(shared_ptr<IntersectionPoint> pt1,
				 shared_ptr<IntersectionPoint> pt2)
//===========================================================================
{
    // Purpose: Check if a connection between two intersection points
    // at the surface boundaries is possible

    double param[4];
    double dist;
    return findMiddlePoint(pt1, pt2, param, dist);

}

//===========================================================================
bool SfSfIntersector::findMiddlePoint(shared_ptr<IntersectionPoint> pt1,
				      shared_ptr<IntersectionPoint> pt2,
				      double param[], double& dist)
//===========================================================================
{
    // Try to iterate to an intersection point at a constant parameter curve
    // between the two points. First define the curve.
    vector<double> par1 = pt1->getPar();
    vector<double> par2 = pt2->getPar();
    double maxdist = 0.0, pdist;
    int max_ind=-1, ki;
    int size = (int)(par1.size());
    for (ki=0; ki<size; ki++)
    {
	pdist = fabs(par1[ki]-par2[ki]);
	if (pdist > maxdist)
	{
	    maxdist = pdist;
	    max_ind = ki;
	}
    }
    if (max_ind < 0)
	return false;   // Do not connect

    ParamSurfaceInt *surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt *surf2 = obj_int_[1]->getParamSurfaceInt();
    shared_ptr<ParamCurve> crv;
    shared_ptr<ParamSurface> srf;
    RectDomain domain;  
    int cv_idx, sf_idx;
   if (max_ind < 2)
    {
	crv = surf1->getConstantParameterCurve(max_ind,
					       0.5*(par1[max_ind]
						    +par2[max_ind]));
	srf = surf2->getParentParamSurface(domain);
	cv_idx = 0;
	sf_idx = 1;
    }
   else
    {
	crv = surf2->getConstantParameterCurve(max_ind-2,
					       0.5*(par1[max_ind]
						    +par2[max_ind]));
	srf = surf1->getParentParamSurface(domain);
	cv_idx = 2;
	sf_idx = 0;
    }
	
   // Make seed
   double seed[3];
   int kj;
   for (ki=0, kj=0; ki<size; ki++)
   {
       if (ki == max_ind)
	   continue;
       seed[kj++] = 0.5*(par1[ki]+par2[ki]);
   }

   // Iterate to a closest point
   double cpos, gpos2[2];
   Point ptc, pts;
   ClosestPoint::closestPtCurveSurf(crv.get(), srf.get(), epsge_->getEpsge(),
		      crv->startparam(),
		      crv->endparam(), &domain, seed[cv_idx], seed+sf_idx,
		      cpos, gpos2, dist, ptc, pts);
   if (cv_idx == 0)
   {
       param[max_ind] = 0.5*(par1[max_ind]+par2[max_ind]);
       param[1-max_ind] = cpos;
       param[2] = gpos2[0];
       param[3] = gpos2[1];
   }
   else
   {
       param[0] = gpos2[0];
       param[1] = gpos2[1];
       param[max_ind] = 0.5*(par1[max_ind]+par2[max_ind]);
       param[5-max_ind] = cpos;
   }

   if (dist < epsge_->getEpsge())
       return true;  // An intersection point is found
   else
       return false;
}


//===========================================================================
void SfSfIntersector::
writeDebugLinear(vector<shared_ptr<IntersectionPoint> >& bd_ints)
//===========================================================================
{
    // Purpose: Write debug info to a file

    static int number = 500;
    number++;
    char buffer[20];
    sprintf(buffer, "debug/linear.out");
    sprintf(buffer, "%i",number);

    ParamSurfaceInt* surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt* surf2 = obj_int_[1]->getParamSurfaceInt();
    shared_ptr<ParamSurface> srf1 = surf1->getParamSurface();
    shared_ptr<ParamSurface> srf2 = surf2->getParamSurface();

    // Open debug output

   std::ofstream debug(buffer);
    int ki, kj;
    for (ki=0; ki<int(bd_ints.size()); ki++)
    {
	vector<double> par = bd_ints[ki]->getPar();
	for (kj=0; kj<int(par.size()); kj++)
	    debug << par[kj] << " ";
	debug << std::endl;
    }
    srf1->writeStandardHeader(debug);
    srf1->write(debug);
    srf2->writeStandardHeader(debug);
    srf2->write(debug);
    
}


//===========================================================================
void SfSfIntersector::
writeDebugConnect(vector<pair<shared_ptr<IntersectionPoint>,
		  IntPtClassification> >& bd_ints)
//===========================================================================
{
    // Purpose : Write debug info to a file

    static int number = 0;
    number++;
    char buffer[20];
    sprintf(buffer, "debug/connect.out");
    sprintf(buffer, "%i",number);

    ParamSurfaceInt* surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt* surf2 = obj_int_[1]->getParamSurfaceInt();
    shared_ptr<ParamSurface> srf1 = surf1->getParamSurface();
    shared_ptr<ParamSurface> srf2 = surf2->getParamSurface();

    // Open debug output

   std::ofstream debug(buffer);
    int ki, kj;
    for (ki=0; ki<int(bd_ints.size()); ki++)
    {
	vector<double> par = bd_ints[ki].first->getPar();
	for (kj=0; kj<int(par.size()); kj++)
	    debug << par[kj] << " ";
	debug << bd_ints[ki].second << std::endl;
    }
    srf1->writeStandardHeader(debug);
    srf1->write(debug);
    srf2->writeStandardHeader(debug);
    srf2->write(debug);
    
}


//===========================================================================
bool SfSfIntersector::complexityReduced()
//===========================================================================
{
    // Purpose: Check if the complexity of the intersection problem is
    // reduced during the last recursive subdivisions.

    // Check for divided high priority singularity
    if (singularity_info_.get() == 0 && prev_intersector_ &&
	prev_intersector_->numParams() == 4 && 
	prev_intersector_->hasSingularityInfo())
    {
	singularity_info_ = (shared_ptr<SingularityInfo>)
	    (new SingularityInfo(prev_intersector_->getSingularityInfo()));
    }
    if (singularity_info_
	&& singularity_info_->hasHighPriSing() == DIVIDED_SING) {
	// Check if the singularity lie within the current surfaces
	vector<double> sing = singularity_info_->getHighPriSing();
	double ptol = epsge_->getRelParRes();
	int ki;
	for (ki=0; ki<2; ki++)
	    if (sing[ki] < obj_int_[0]->startParam(ki) - ptol ||
		sing[ki] > obj_int_[0]->endParam(ki) + ptol)
		break;
	if (ki == 2)
	{
	    for (ki=0; ki<2; ki++)
		if (sing[2+ki] < obj_int_[1]->startParam(ki) - ptol ||
		    sing[2+ki] > obj_int_[1]->endParam(ki) + ptol)
		    break;
	    if (ki == 2)
		return false;  // This is a complex case
	}
    }
    
    if (!hasComplexityInfo())
	return true;  // Necessary prerequisite information is not fetched.
                      // Continue subdividing

     // Fetch number of intersection points
    complexity_info_->setNmbIntpts(int_results_->numIntersectionPoints());
    complexity_info_->setNmbSingpts(int_results_->numSingularIntersectionPoints());
    
   // Check sufficient recursion level. At least a few recursions are required
//    int min_rec = 5;
     int min_rec = 6;
   if (nmbRecursions() < min_rec)
	return true;   // Continue subdividing


    // Compare with the situation in the previous recursion level
    shared_ptr<ComplexityInfo> prev_complexity = 
	prev_intersector_->getComplexityInfo();

    if (prev_complexity.get() == 0)
	return true;  // Previous intersector is not surface-surface.

    double overlap_fac = 0.95;
    if (complexity_info_->isComplex())
	return false;
    else if ((complexity_info_->getBoxOverlap() > 
	overlap_fac*prev_complexity->getBoxOverlap() &&
	complexity_info_->getConeOverlap() >=
	prev_complexity->getConeOverlap()) ||
	complexity_info_->getNmbIntpts() > prev_complexity->getNmbIntpts() ||
	complexity_info_->getNmbSingpts() > prev_complexity->getNmbSingpts())
	return false;
    else
	return true;
}


//===========================================================================
int SfSfIntersector::linearCase()
//===========================================================================
{
    // Purpose: Given two linear parametric surfacess, compute the
    // intersection.

    // First check existance of any intersection points
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    if (int_pts.size() <= 1)
	return 0;   // One or less intersection point at the boundary.
                    // No further intersections are expected

    if (int_pts.size() == 2)
    {
	// Two intersection points. Make sure that they are connected
	if (!(int_pts[0]->isConnectedTo(int_pts[1])))
	    int_pts[0]->connectTo(int_pts[1], LINEAR_SFSF);
	return 1;
    }
 	
    // More than two intersection points. This is an unexpected sitution.
    // Make debug output.
    writeDebugLinear(int_pts);
    return 0;

}


//===========================================================================
int SfSfIntersector::doSubdivide()
//===========================================================================
{
    // Purpose : Given two parametric surfaces and no simple case
    // situation, subdivide the surfaces to produce more simple
    // subproblems. Perform intersection with the subdivision points.
    // Prepare for the next recursion level.
    //
    // NB! This function is almost identical to
    // SfCvIntersector::doSubdivide and also CvCv::doSubdivide can be
    // made to follow this pattern (that would make it a bit more
    // complex than necessary). Thus, the function might be moved to
    // Intersector2Obj, but this is delayed until I have a better
    // picture of the complete complexity of the function.

    // !!! DEBUG
    if (getenv("SUBDIV_SFSF") && *(getenv("SUBDIV_SFSF"))=='1') {
	if (!getenv("DEBUG") ||
	    (getenv("DEBUG") && *getenv("DEBUG") == '0')) {
	    cout << "================================================"
		 << endl;
	    cout << "Domain 1: ";
	    for (int ki = 0; ki < 2; ki++) {
		cout << obj_int_[0]->startParam(ki) << " ";
		cout << obj_int_[0]->endParam(ki) << " ";
	    }
	    cout << endl;
	    cout << "Domain 2: ";
	    for (int ki = 0; ki < 2; ki++) {
		cout << obj_int_[1]->startParam(ki) << " ";
		cout << obj_int_[1]->endParam(ki) << " ";
	    }
	    cout << endl;
	}
    }

    // Sort the parameter directions according to importance of
    // subdivision according to properties of the objects and already
    // computed intersections
    int perm[4];  // Three parameter directions that need sorting
    int deg_edge[4];  // =0, not degenerate, =1 degenerate at start,
                      // =2 degenerate at endpoint, =3 degenerate in
                      // both endpoints
    int nmb_subdiv = sortParameterDirections(perm, deg_edge);

    if (nmb_subdiv == 0) {
	// Not possible to subdivide neither the curve nor the surface
	return 0;   
    }

    // Declare intersection objects belonging to the next recursion
    // level.  First they point to the objects at this level. Then the
    // current objects are replaced as subdivision is taking place.
    vector<shared_ptr<ParamGeomInt> > sub_objects;
    sub_objects.push_back(obj_int_[0]);
    sub_objects.push_back(obj_int_[1]);
    // Currently one instance of each intersection object at this level
    int numobj[2];
    numobj[0] = 1;
    numobj[1] = 1;    

    // For each parameter direction in prioritized order, fetch an
    // appropriate subdivision parameter, and perform subdivision
    double subdiv_par;
    SubdivisionClassification found;
    vector<shared_ptr<ParamGeomInt> > subdiv_objs;
    vector<shared_ptr<ParamGeomInt> > bd_objs;
    bool use_sing = false;

    for (int ki = 0; ki < nmb_subdiv; ki++) {
	found = getSubdivisionParameter(perm[ki], deg_edge[perm[ki]],
					subdiv_par);
	if (found == DIVIDE_SING) {
	    use_sing = true;
	}

	if (found == CANNOT_SUBDIVIDE) {
	    // No parameter value is found. Move this parameter
	    // direction to the end of the permutation array, and
	    // decrease the number of parameter directions where it is
	    // possible to subdivide.
	    std::swap(perm[ki], perm[nmb_subdiv-1]);
	    nmb_subdiv--;
	    ki--;
	    continue;
	}

	// !!! DEBUG
	if (getenv("SUBDIV_SFSF") && (*getenv("SUBDIV_SFSF")) == '1') {
	    cout << "Subdivide dir = " << perm[ki]
		 << " par = " << subdiv_par
		 << " criterium = " << found << endl;
	}
	
	// Subdivide the current object (or objects - if the current
	// object is the surface then it can have been subdivided in
	// one parameter direction already and we have two subsurfaces
	// to subdivide)
	int nmbdir1 = obj_int_[0]->numParams();
	int idxobj = (perm[ki] >= nmbdir1) ? 1 : 0;
	int idx = idxobj * numobj[0];
	int pdir = perm[ki] - (idxobj * nmbdir1);
	for (int kj = 0; kj < numobj[idxobj]; kj++) {
	    subdiv_objs.clear();
	    bd_objs.clear();
	    try {
		sub_objects[idx+kj]->subdivide(pdir, subdiv_par, 
					       subdiv_objs, bd_objs);
	    } catch (...) {
		subdiv_objs.clear();
		bd_objs.clear();
	    }
	    if (subdiv_objs.size() < 1 || bd_objs.size() == 0) {
		continue;  // No new objects 
	    }

	    for (int kr = 0; kr < int(subdiv_objs.size()); kr++) {
		subdiv_objs[kr]->setParent(obj_int_[idxobj].get());
	    }

	    if (found == DIVIDE_DEG) {
		// Subdivision to create a small, degenerate triangle.
		// Set flag.
		setDegTriangle(subdiv_objs, deg_edge[perm[ki]],
			       pdir, subdiv_par);
	    }
      
	    for (int kr = 0; kr < int(bd_objs.size()); kr++) {
		// Intersect the subdivision object with the other object
		// @@@ VSK Faktorenes orden er ikke likegyldig!!
		shared_ptr<ParamGeomInt> obj1 
		    = (idxobj==0) ? bd_objs[kr] : obj_int_[0];
		shared_ptr<ParamGeomInt> obj2
		    = (idxobj==0) ? obj_int_[1] : bd_objs[kr];
		shared_ptr<Intersector> lower_intersector 
		    = lowerOrderIntersector(obj1, obj2,
					    this, perm[ki], subdiv_par);

		// Is it here relevant to fetch existing intersection
		// points and/or insert points into intersection
		// curves before computing intersection points with
		// the subdivision points. Normally, this is mostly of
		// interest when surfaces are involved, but
		// intersection intervals might exist?  These
		// computations do anyway involve the intersection
		// pool, but they might be trigged somehow. The
		// parameter direction and value are required
		// information

		if (found == DIVIDE_SING && hasSingularityInfo()) {
		    // Transfer singularity information to the child
		    // process
		    lower_intersector->setSingularityInfo
			(getSingularityInfo(), perm[ki]);
		}

		if (getenv("SUBDIV_SFSF")
		    && (*getenv("SUBDIV_SFSF")) == '1') {
		    vector<shared_ptr<IntersectionPoint> > ipoint1;
		    lower_intersector->getIntPool()->getIntersectionPoints
			(ipoint1);
		    cout << "Intersection points prior : " << endl;
		    for (size_t k1 = 0; k1 < ipoint1.size(); k1++) {
			int ctr = 0;
			for (int k2 = 0; k2 < 3; k2++) {
			    if (ctr == perm[ki]) {
				cout << subdiv_par << "  ";
			    }
			    cout << ipoint1[k1]->getPar(k2) << "  ";
			    ++ctr;
			}
			if (perm[ki] == 3) {
			    cout << subdiv_par << "  ";
			}
			cout << "  " << ipoint1[k1]->getDist();
			cout << endl;
		    }
		}
	      
		lower_intersector->compute(false);

		if (getenv("SUBDIV_SFSF")
		    && (*getenv("SUBDIV_SFSF")) == '1') {
		    vector<shared_ptr<IntersectionPoint> > ipoint1;
		    lower_intersector->getIntPool()->getIntersectionPoints
			(ipoint1);
		    cout << "Intersection points post : " << endl;
		    for (size_t k1 = 0; k1 < ipoint1.size(); k1++) {
			int ctr = 0;
			for (int k2 = 0; k2 < 3; k2++) {
			    if (ctr == perm[ki]) {
				cout << subdiv_par << "  ";
			    }
			    cout << ipoint1[k1]->getPar(k2) << "  ";
			    ++ctr;
			}
			if (perm[ki] == 3) {
			    cout << subdiv_par << "  ";
			}
			cout << "  " << ipoint1[k1]->getDist();
			cout << endl;
		    }
		}

		// Check quality of intersection points.  If the
		// quality is not sufficient, find a new subdivision
		// parameter and repeat the process.  Otherwise
		int nmb_orig = int_results_->numIntersectionPoints();
		int_results_->includeReducedInts
		    (lower_intersector->getIntPool());
		if (/*found != DIVIDE_SING && */
		    nmb_orig < int_results_->numIntersectionPoints()) {
		    postIterate2(nmb_orig, perm[ki], subdiv_par);
		    // Swap directions 0 <-> 1 or 2 <-> 3
		    int pdir = (perm[ki] < 2) ? 1-perm[ki] : 5-perm[ki];
		    postIterate2(nmb_orig, pdir, subdiv_par, false);
		    if (true /*!use_sing*/) {
			postIterate(0, perm[ki]);
		    }
		}
	    }

	    // Replace pointers to current objects with pointers to
	    // the sub objects after subdivision. The new objects can
	    // be the final sub objects or they can be subdivided once
	    // more depending on the number of parameter directions
	    // elected for subdivision in the current object.
	    sub_objects[idx+kj] = subdiv_objs[0];
	    sub_objects.insert(sub_objects.begin()+idx+kj+1,
			       subdiv_objs.begin()+1, subdiv_objs.end());
	    numobj[idxobj] += ((int)subdiv_objs.size()-1);
	    kj += ((int)subdiv_objs.size()-1);

	}
    }

    if (nmb_subdiv == 0) {
	return 0;  // Not subdivided
    }

    // Create new intersector objects
    sub_intersectors_.clear();
    for (int ki = 0; ki < numobj[0]; ki++) {
	for (int kj = 0; kj < numobj[1]; kj++) {
	    int kk = numobj[0] + kj;
	    shared_ptr<Intersector> intersector 
		= shared_ptr<Intersector>
		(new SfSfIntersector(sub_objects[ki], sub_objects[kk],
				     epsge_, this));
	    sub_intersectors_.push_back(intersector);
	}
    }

    return 1;
}


//===========================================================================
void SfSfIntersector::
getApproxImplicit(vector<vector<shared_ptr<Param2FunctionInt> > >&
		  approx_implicit,
		  vector<double>& approx_implicit_err,
		  vector<double>& approx_implicit_gradsize,
		  vector<double>& approx_implicit_gradvar)
//===========================================================================
{
    approx_implicit = approx_implicit_;
    approx_implicit_err = approx_implicit_err_;
    approx_implicit_gradsize = approx_implicit_gradsize_;
    approx_implicit_gradvar = approx_implicit_gradvar_;
}


//===========================================================================
int SfSfIntersector::sortParameterDirections(int perm[],
					     int deg_edge[]) 
//===========================================================================
{
    // deg_edge[index]=0, not degenerate, =1 degenerate at start, =2
    // degnerate at edge, =3 degenerate in both endpoints


    // Purpose: Fetch information related to the surfaces in order to
    // decide which one that can be subdivided and which parameter
    // direction is most important to subdivide.

  // @@@ VSK. This routine is looks similar to the one for curve-curve
  // and surface-curve intersections, but the ordering of criteria
  // differs. Or?

  double length[4], wiggle[4];  // Estimated length and wiggliness
  bool inner_knots[4], critical_val[4], can_divide[4];
  bool has_inner_ints[4];
  double rel_wiggle = 0.1;
  double rel_length = 0.1;
  double treshhold_wiggle = 0.5;
  int nmbdir[2];
  int ki, kj, kr;

  SingularityClassification sing_type = hasSingularityInfo() ?
      singularity_info_->hasHighPriSing() : NO_SING;

  nmbdir[0] = obj_int_[0]->numParams();
  nmbdir[1] = obj_int_[1]->numParams();

  // Fetch information from the surfaces
  obj_int_[0]->getLengthAndWiggle(length, wiggle);
  obj_int_[1]->getLengthAndWiggle(length+nmbdir[0], wiggle+nmbdir[0]);
  
  double med_len = 0.0, med_wiggle = 0.0;
  double min_len=length[0], max_len=length[0];
  for (ki=0; ki<4; ki++)
  {
      med_len += length[ki];
      med_wiggle += wiggle[ki];
      min_len = std::min(min_len, length[ki]);
      max_len = std::max(max_len, length[ki]);
  }
  med_len /= (double)4;
  med_wiggle /= (double)4;

  double min_length = std::max(0.5*med_len, epsge_->getEpsge());
  double min_wiggle = std::max(0.5*med_wiggle, 0.02);

  for (ki=0, kr=0; ki<2; ki++)
  {
      // Make sure that degeneracy information is computed
      //bool tmp_deg = 
      //obj_int_[0]->isDegenerate(epsge_->getEpsge(), 0);
      for (kj=0; kj<nmbdir[ki]; kj++, kr++)
      {
	  inner_knots[kr] = obj_int_[ki]->hasInnerKnots(kj);

	  critical_val[kr] = obj_int_[ki]->hasCriticalVals(kj);

	  has_inner_ints[kr] = int_results_->hasPointsInInner(kr);

	  can_divide[kr] = (obj_int_[ki]->canDivide(kj) && 
			    (obj_int_[ki]->canDivideTinyTriang(kj) || has_inner_ints[kr]));

	  // This is only relevant for surfaces!!! How to handle this?
	  deg_edge[kr] = obj_int_[ki]->isDegenerate(epsge_->getEpsge(), kj);
	  // Given a degenerate edge, it is necessary to check if there is any 
	  // intersection points situated at this edge. How?
	  if (deg_edge[kr])
	  {
	      double par;
	      if (deg_edge[kr] == 1 || deg_edge[kr] == 3)
	      {
		  par = obj_int_[ki]->startParam(kj);
		  if (!(int_results_->existIntersectionPoint(kr, par)))
		      deg_edge[kr] -= 1;
	      }
	      if (deg_edge[kr] == 2 || deg_edge[kr] == 3)
	      {
		  par = obj_int_[ki]->endParam(kj);
		  if (!(int_results_->existIntersectionPoint(kr, par)))
		      deg_edge[kr] -= 2;
	      }
	  }
      }
      // Avoid subdivision in a parameter direction opposite to a
      // degenerate edge unless in a degenerate triangle situation
      if (deg_edge[2*ki] && can_divide[2*ki] && !deg_edge[2*ki+1])
	  can_divide[2*ki+1] = false;
      if (!deg_edge[2*ki] && deg_edge[2*ki+1] && can_divide[2*ki+1])
	  can_divide[2*ki] = false;
  }


  // VSK, 0705. If the size of the surfaces vary very much, avoid subdividing in the 
  // parameter direction where the surface is shortest
  for (ki=0; ki<4; ki++)
  {
      if (length[ki] < rel_length*max_len && wiggle[ki] < treshhold_wiggle)
	  can_divide[ki] = false;
  }

  // @@@ VSK Need to divide out a degenerate edge must be included
  // Is this also relevant for surface-curve? Then the responsibility could
  // be left to the surface, and we can let the surface return something
  // like must_divide(direction)

  // Fetch information from the intersection pool according to inner 
  // intersection points

  // Number of parameter directions is four.
  int size = 4;

  int curr = 0;
  int min_nmb = 0;

  // Initiate permutation array
  for (ki=0; ki<size; ki++)
    perm[ki] = ki;

  // Sort according to the given values
  // First sort out the directions where subdivision is impossible
  for (ki=0; ki<size; ki++)
    {
      if (!can_divide[perm[ki]])
	{
	  if (ki < size-1)
	    std::swap(perm[ki], perm[size-1]);
	  ki--;
	  size--;
	}
    }

  // First priority is the need to divide out degenerate edges
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if (deg_edge[perm[kj]] && !deg_edge[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	}
  for (ki=curr; ki<size-1; ki++)
      if (deg_edge[perm[ki]] && !deg_edge[perm[ki+1]])
	  curr++;

  // Second priority is parameter directions with critical values
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if (critical_val[perm[kj]] && !critical_val[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	}
  for (ki=curr; ki<size-1; ki++)
      if (critical_val[perm[ki]])
	  curr++;

  // Check for critical singularity
  if (sing_type == INIT_SING)
  {
      // We need to divide out this singularity. No point in
      // further sorting
      return size;
  }

  // @@@ VSK, Maybe singular intersection points at the boundaries should 
  // be a criterium? 

  // Next criterium is inner knots
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if (inner_knots[perm[kj]] && !inner_knots[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	}
  for (ki=curr; ki<size-1; ki++)
      if (inner_knots[perm[ki]])
	  curr++;

  // Existing intersection points internal to the parameter interval
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if (has_inner_ints[perm[kj]] && !has_inner_ints[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	}
  for (ki=curr; ki<size; ki++)
      if (has_inner_ints[perm[ki]])
	  curr++;

  // Finally length and curvature
  // Sort with regard to length and wiggle first for the directions with inner knots
  for (ki=0; ki<curr; ki++)
    for (kj=ki+1; kj<size; kj++)
      if ((length[perm[kj]] > length[perm[ki]] && 
	   wiggle[perm[kj]] > wiggle[perm[ki]]) ||
	  (0.2*length[perm[kj]] > length[perm[ki]] && 
	   (wiggle[perm[kj]] > rel_wiggle*wiggle[perm[ki]] ||
	       wiggle[perm[ki]] < min_wiggle)) ||
	  (0.2*wiggle[perm[kj]] > wiggle[perm[ki]] && 
	   (length[perm[kj]] > rel_length*length[perm[ki]] ||
	    length[perm[ki]] < min_length)))
	{
	  std::swap(perm[ki], perm[kj]);
	}

  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if ((length[perm[kj]] > length[perm[ki]] && 
	   wiggle[perm[kj]] > wiggle[perm[ki]]) ||
	  (0.2*length[perm[kj]] > length[perm[ki]] && 
	   (wiggle[perm[kj]] > rel_wiggle*wiggle[perm[ki]] ||
	       wiggle[perm[ki]] < min_wiggle)) ||
	  (0.2*wiggle[perm[kj]] > wiggle[perm[ki]] && 
	   (length[perm[kj]] > rel_length*length[perm[ki]] ||
	    length[perm[ki]] < min_length)))
	{
	  std::swap(perm[ki], perm[kj]);
	}

  // Check that the minimum conditions for subdivision is satisfied
  for (ki=size-1; ki>=min_nmb; ki--, size--)
      if (length[perm[ki]] >= min_length || wiggle[perm[ki]] >= min_wiggle)
	  break;

  return size;
}


//===========================================================================
SubdivisionClassification 
SfSfIntersector::getSubdivisionParameter(int dir, int deg_edge, double& par)
//===========================================================================
{
    bool found;
    double aeps = epsge_->getEpsge();

    // @@@ VSK. For instance for degenerate edges, it could be
    // relevant to return two subdivision parameters. Should this be
    // implemented right away? Is it easier to do it now than later?
    int ki, kj, kr, sgn;
    double treshhold = 10.0*aeps;
    int nmbdir1 = obj_int_[0]->numParams();
    //int nmbdir2 = obj_int_[1]->numParams();

    SingularityClassification sing_type = hasSingularityInfo() ?
	singularity_info_->hasHighPriSing() : NO_SING;

    // Set pointer to the intersection object corresponding to the
    // parameter direction
    ParamGeomInt *obj = 
	(dir < nmbdir1) ? obj_int_[0].get() : obj_int_[1].get();
    ParamSurfaceInt *surf = obj->getParamSurfaceInt();
    ParamSurfaceInt *surf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt *surf2 = obj_int_[1]->getParamSurfaceInt();
    int pdir = (dir < nmbdir1) ? dir : dir-nmbdir1;
    double ta = obj->startParam(pdir);
    double tb = obj->endParam(pdir);
    double frac = 0.1; //0.01;
    double frac2 = 0.001;
    double frac3 = 0.02;

    // Fetch relevant intersection points
    vector<shared_ptr<IntersectionPoint > > ipoints;
    int_results_->getIntersectionPoints(ipoints);

    // First divide out degenerate edges
    if (deg_edge) {
	found = getSubdivParDegEdge(surf, dir, pdir, deg_edge,
				    treshhold, par);

	// Check that the parameter is inside the legal interval
	// @@@ VSK. It might be that we must include an epsilon to
	// avoid the parameter to be too close to the boundary
	if (found && par > ta && par < tb)
	    return DIVIDE_DEG;
    }

    // Get critical parameters
    // @@@ Critical parameters of different priority? Sorting?
    vector<double> critical_pars = obj->getCriticalVals(pdir);

    int size = (int)critical_pars.size();
    int is_critical = 0;
    if (size > 0) {
	// Check suitability of the critical parameters
	// @@@ VSK. inInfluenceArea must also return the intersection
	// points where the subdivision parameter belongs to the
	// influence area.  Then it is possible to test if the
	// complete constant parameter curve belongs to the
	// intersection. Should we subdivide here even if this
	// constant direction intersection is of low quality?  Is
	// there a threshold for when it is interesting?  How much
	// would this test cost? Can it be run for every candidate?
	// It would be nice to be able to march along a constant
	// parameter curve in one surface without having to pick the
	// curve.

	for (ki=0; ki<size; ki++) {
	    vector<shared_ptr<IntersectionPoint> > int_pts;
	    is_critical = int_results_->inInfluenceArea(dir, critical_pars[ki],
							int_pts);
	    if (int_pts.size() > 0) {
		// Check if the subdivision position is OK with
		// respect to the angle between the subdivision curve
		// and the intersection curve
		is_critical = checkSubdivParam(dir, critical_pars[ki], ta, tb,
					       int_pts);
	    }
	    if (is_critical == 0 || is_critical == 2) {
		par = critical_pars[ki];
		int is_OK;
		is_OK = splitIntResults(int_pts, dir, par, ta, tb, false);
		return DIVIDE_CRITICAL;
	    }
	}
    }

    // Check for critical singularity
    if (sing_type == INIT_SING || sing_type == KEEP_SING) {
	// Check if the singularity lie in a corner in both surfaces
	vector<double> sing_par = singularity_info_->getHighPriSing();
	if (surf1->inCorner(&sing_par[0], epsge_->getRelParRes()) &&
	    surf2->inCorner(&sing_par[2], epsge_->getRelParRes())) {
      
	    // We need to divide out this singularity.
	    // Find suitable division value.
	    // First get the parameter values of the singularity
	    double min_dist = 1000.0*epsge_->getRelParRes();
	    double tol = 10000.0*aeps;
	    double high_sing = sing_par[dir];
	    if (fabs(high_sing-ta) < epsge_->getRelParRes()) {
		par = surf->getParOffBd(pdir, true, tol);
		par = std::max(par, ta + frac3*(tb-ta));
		if (par > ta+min_dist && par < tb-frac*(tb-ta))
		    return DIVIDE_HIGH_SING;
	    }
	 
	    else if (fabs(tb-high_sing) < epsge_->getRelParRes()) {
		par = surf->getParOffBd(pdir, false, tol);
		par = std::min(par, tb - frac3*(tb-ta));
		if (par > ta+frac*(tb-ta) && par < tb-min_dist)
		    return DIVIDE_HIGH_SING;
	    }
	}
	else {
	    // Set information in singularity information
	    singularity_info_->setHighPriSingType(KEEP_SING);
	    // Divide in the singularity
	    par = sing_par[dir];
	    if (par > ta + frac2*(tb-ta) && par < tb - frac2*(tb-ta))
		return DIVIDE_HIGH_SING;
	    else 
		return CANNOT_SUBDIVIDE;
	}
    }

    // Iterate for a singular branch point in which to subdivide
    found = getSubdivAtSing(dir, ta, tb, par);
    if (found)
	return DIVIDE_SING; 

    // Look for a suitable knot in which to subdivide. First fetch
    // the knots sorted according to the distance to the mid-parameter
    vector<double> knot_vals = obj->getInnerKnotVals(pdir, true);
    size = (int)knot_vals.size();
    if (size > 0) {
	// Check suitability of the knots
	for (ki=0; ki<size; ki++) {
	    vector<shared_ptr<IntersectionPoint> > int_pts;
	    is_critical = int_results_->inInfluenceArea(dir, knot_vals[ki], 
							int_pts);
	    if (int_pts.size() > 0) {
		// Check if the subdivision position is OK with
		// respect to the angle between the subdivision curve
		// and the intersection curve
		is_critical = checkSubdivParam(dir, knot_vals[ki], ta, tb,
					       int_pts);
	    }

	    // Check for very close intersection points where the current parameter
	    // value is not in the influence area
	    for (kj=0; kj<(int)ipoints.size(); kj++)
	    {
		for (kr=0; kr<(int)int_pts.size(); kr++)
		    if (ipoints[ki] == int_pts[kr])
			break;
		if (kr == (int)int_pts.size())
		{
		    if (fabs(ipoints[kj]->getPar(dir) - knot_vals[ki]) < frac2*(tb-ta) &&
			fabs(ipoints[kj]->getPar(dir) - knot_vals[ki]) > epsge_->getRelParRes())
		    {
			is_critical = 1;
			break;  // Not a good subdivision parameter
		    }
		}
	    }

	    if (is_critical == 0 || is_critical == 2) {
		par = knot_vals[ki];
		int is_OK = splitIntResults(int_pts, dir, par, ta, tb, false);
		if (is_OK)
		    return DIVIDE_KNOT;
	    }
	}
    }
 
    // Check for intersection points internal to a boundary in which
    // to subdivide
    vector<shared_ptr<IntersectionPoint > > ip = ipoints;
    size = (int)ip.size();
//     vector<double> inner_ints = int_results_->getSortedInnerInts(dir);

//     // Insert parameter values between the intersection points
//     for (ki=2; ki<inner_ints.size()-1; ki++) {
// 	par = 0.5*(inner_ints[ki-1]+inner_ints[ki]);
// 	inner_ints.insert(inner_ints.begin()+ki, par);
// 	ki++;
//     }
//     size = inner_ints.size();

   
    // Sort according to angle between surface normal and distance
    // from midpoint
    double angtol = 1.0e-3;
    //int nsing = 0;
    vector<double> angle(size);
    double mind = 10.0*aeps;
    double maxd = 0.0;
    double dist;
    for (ki=0; ki<size; ki++) {
	// Compute surface normals and accuracies
	Point norm1, norm2;
	surf1->normal(ip[ki]->getPar1()[0], ip[ki]->getPar1()[1], norm1);
	surf2->normal(ip[ki]->getPar2()[0], ip[ki]->getPar2()[1], norm2);
	angle[ki] = norm1.angle(norm2);
	dist = ip[ki]->getDist();
	if (dist < mind)
	    mind = dist;
	if (dist > maxd)
	    maxd = dist;
    }
    for (ki=0; ki<size; ki++) {
	for (kj=ki+1; kj<size; kj++) {
	    if (angle[kj] < angle[ki]) {
		std::swap(ip[ki], ip[kj]);
		std::swap(angle[ki], angle[kj]);
	    }
	}
    }

    // Dismiss bad intersection points as a criterion for subdivision
    double tlevel;
    double fac = 100.0;
    if (maxd/aeps > fac*mind/aeps)
	tlevel = fac*aeps;
    else
	tlevel = aeps;

    for (ki=0; ki<size; ki++)
    {
	if (angle[ki] < angtol)
	    continue;
	if (ip[ki]->getDist() > tlevel)
	{
	    ip.erase(ip.begin()+ki);
	    angle.erase(angle.begin()+ki);
	    size--;
	}
    }

	
    double mid = 0.5*(ta + tb);
    for (ki=0; ki<size; ki++) {
	if (angle[ki] < angtol)
	    continue;
	for (kj=ki+1; kj<size; kj++) {
	    if (fabs(ip[kj]->getPar(dir)-mid) < fabs(ip[ki]->getPar(dir)-mid))
		std::swap(ip[ki], ip[kj]);
	}
    }
    if (size > 0) {
	// Check suitability of the intersection points
	//        for (ki=size/2, kj=1, sgn=-1; ki<size && ki>=0; 
	//  	   ki+=sgn*kj, kj++, sgn*=-1)
	for (ki=0; ki<size; ki++) {
	    par = ip[ki]->getPar(dir);
	    //	    par = inner_ints[ki];
	    if (par < ta+frac*(tb-ta) ||
		par > tb-frac*(tb-ta))
		continue;

	    vector<shared_ptr<IntersectionPoint> > int_pts;
	    is_critical = int_results_->inInfluenceArea(dir, par, int_pts);
	    if (int_pts.size() > 0) {
		// Check if the subdivision position is OK with
		// respect to the angle between the subdivision curve
		// and the intersection curve
		is_critical = checkSubdivParam(dir, par, ta, tb, int_pts);
	    }

	    // Check for very close intersection points where the current parameter
	    // value is not in the influence area
	    for (kj=0; kj<(int)ipoints.size(); kj++)
	    {
		if (ipoints[kj] == ip[ki])
		    continue;  // This point

		for (kr=0; kr<(int)int_pts.size(); kr++)
		    if (ipoints[ki] == int_pts[kr])
			break;
		if (kr == (int)int_pts.size())
		{
		    if (fabs(ipoints[kj]->getPar(dir) - par) < frac2*(tb-ta) &&
			fabs(ipoints[kj]->getPar(dir) - par) > epsge_->getRelParRes())
		    {
			is_critical = 1;
			break;  // Not a good subdivision parameter
		    }
		}
	    }

	    if (is_critical == 0 || is_critical == 2) {
		int is_OK = splitIntResults(int_pts, dir, par, ta, tb, true);
		if (is_OK)
		    return DIVIDE_INT;
	    }
	}
    }

    // Subdivide at a suitable parameter
    double divpar = 0.5*(ta+tb);
    double del = 0.1*(tb-ta);
    double tint = del;
    for (ki=0, sgn=-1; ki<9; ki++, divpar+=sgn*tint, tint+=del, sgn*=-1) 
    {
	vector<shared_ptr<IntersectionPoint> > int_pts;
	is_critical = int_results_->inInfluenceArea(dir, divpar, int_pts);
	if (int_pts.size() > 0) {
	    // Check if the subdivision position is OK with respect to the
	    // angle between the subdivision curve and the intersection curve
	    is_critical = checkSubdivParam(dir, divpar, ta, tb, int_pts);
	}

	// Check for very close intersection points where the current parameter
	// value is not in the influence area
	for (kj=0; kj<(int)ipoints.size(); kj++)
	{
	    if (ipoints[kj] == ip[ki])
		continue;  // This point

	    for (kr=0; kr<(int)int_pts.size(); kr++)
		if (ipoints[ki] == int_pts[kr])
		    break;
	    if (kr == (int)int_pts.size())
	    {
		if (fabs(ipoints[kj]->getPar(dir) - divpar) < frac2*(tb-ta) &&
		    fabs(ipoints[kj]->getPar(dir) - divpar) > epsge_->getRelParRes())
		{
		    is_critical = 1;
		    break;  // Not a good subdivision parameter
		}
	    }
	}

	if (is_critical == 0 || is_critical == 2) {
	    par = divpar;
	    // VSK, 0906. Do not move intersection points for such a weak split criterion
	    // splitIntResults(int_pts, dir, par, ta, tb);
	    return DIVIDE_PAR;
	}
    }

    return CANNOT_SUBDIVIDE;  // No subdivision parameter found
}


//===========================================================================
int SfSfIntersector::
splitIntResults(vector<shared_ptr<IntersectionPoint> >& int_pts,
		int pardir, double par, double start, double end,
		bool large_move)
//===========================================================================
{
    // Purpose: Insert intersection points in influence areas to
    // relevant intersection points at positions where one of the
    // current surfaces are subdivided. Alternatively move the already
    // existing intersection point if it is "close".

    int is_OK = 1;
    size_t nmb = int_pts.size();
    size_t ki;
    int npar = 4, kj, kh;
    double dequal = epsge_->getRelParRes();
    double par_tol = 1.0e-6;  // Should be tuned or maybe be dependent
                              // on the size of the current parameter
                              // interval.
    par_tol = std::min(par_tol, 0.01*(end-start));
    double par_tol2 = 0.01;
    par_tol2 = std::min(par_tol2, 0.1*(end-start));
    double eps_frac = 1.0e3;

    double ta[4], tb[4];
    for (ki=0, kh=0; ki<2; ki++)
	for (kj=0; kj<2; kj++, kh++)
	{
	    ta[kh] = obj_int_[ki]->startParam(kj);
	    tb[kh] = obj_int_[ki]->endParam(kj);
	}

    double param[4], seed[4], dist;
    int isolink[4];
    for (ki=0; ki<nmb; ki++) {
	// Iterate to the position of the new intersection point
	for (kj=0; kj<npar; kj++)
	    seed[kj] = int_pts[ki]->getPar(kj);

	if (fabs(seed[pardir] - par) < epsge_->getRelParRes())
	    continue;   // No need to modify point

	doIterate(pardir, par, param, dist, seed);

	// Check other boundaries
	for (kj=0; kj<4; kj++)
	{
	    if (seed[kj]-ta[kj]<=dequal && param[kj]-ta[kj]>dequal)
		param[kj] = ta[kj];
	    else if (tb[kj]-seed[kj]<=dequal && tb[kj]-param[kj]>dequal)
		param[kj] = tb[kj];
	}
	// Update distance
	Point pt1, pt2;
	obj_int_[0]->point(pt1, param);
	obj_int_[1]->point(pt2, param+2);
	dist = pt1.dist(pt2);
	
	//if (dist <= epsge_->getEpsge()) 
	// VSK, 0702. Adds a threshhold regarding the relative accuracy. A better solution would be
	// to make sure that no previous subdivision parameters are passed in any parameter directions,
	// but this requires some interface changes
	if (dist <= int_pts[ki]->getDist() ||
	    (dist <= epsge_->getEpsge() && fabs(param[pardir]-par) <= dequal &&
		fabs(seed[pardir]-par) > dequal && dist < eps_frac*int_pts[ki]->getDist())) 
	{
	    int nmb_links = int_pts[ki]->hasIsoLinks(isolink);
	    for (kj=0; kj<nmb_links; kj++)
		if (isolink[kj] != pardir)
		    break;

	    for (kh=0; kh<npar; kh++) {
		if ((fabs(seed[kh]-start) < dequal
		     || (fabs(seed[kh]-end) < dequal)) &&
		    !(fabs(param[kh]-start) < dequal
		      || (fabs(param[kh]-end) < dequal)))
		    break;
	    }

	    if (kh < npar) {
		// Check if a modified intersection point can do
		param[kh] = seed[kh];
		Point pt1, pt2;
		obj_int_[0]->point(pt1, param);
		obj_int_[1]->point(pt2, param+2);
		if (pt1.dist(pt2) < int_pts[ki]->getDist())
		    kh = npar;
	    }
		

	    bool same_bd = atSameBoundary(seed, param);
	    double pdist = fabs(seed[pardir]-par);
	    if (((pdist < par_tol || (pdist < par_tol2 && large_move)) && same_bd)
		&& kj == nmb_links && kh == npar) {
		// Move position of current intersection point
		if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		{
		std::cout << "SfSfIntersector::splitIntResults, Moved intersection point from ";
		std::cout  << seed[0] << " " << seed[1] << " " << seed[2];
		std::cout << " " << seed[2] << seed[3];
		std::cout << ",dist: " << int_pts[ki]->getDist() << std::endl;
		std::cout<< " to " << param[0] << " " ;
		std::cout << param[1] << " " << param[2] << " " << param[3] ;
		std::cout << " ,dist: " << dist <<std::endl;
		std::cout << "Pardist: " << fabs(seed[pardir]-par) << ", same boundary: " << same_bd;
		std::cout << std::endl;
		}
		int_pts[ki]->replaceParameter(param); //, npar);
	    }
	    else if (kj == nmb_links && kh == npar && same_bd) {
		// Insert an intersection point in the influcence area
 		//int_results_->insertInInfluenceInterval(int_pts[ki], param, 
 		//					pardir);
	    }
	    else if (pdist < par_tol && kj == nmb_links && kh == npar)
	    {
		// This point is a very hot candidate for a move, but it is impossible
		// since we don't know if we move it off a future boundary. When not moved,
		// the point is a candidate source of noice. Avoid the current subdivision
		is_OK = 0;
	    }
	}
    }
    return is_OK;
}


//===========================================================================
void SfSfIntersector::postIterate(int nmb_orig, int dir, bool keep_endpt)
//===========================================================================
{
    // Check the quality of all new interserction points. Iterate
    // along the specified parameter direction
    //double ptol = 1000.0*epsge_->getRelParRes();
    double ptol = 1000.0*epsge_->getRelParRes();
   double ptol2 = 1.0e-8;
   double tol2 = 10.0*epsge_->getEpsge();
    ParamSurfaceInt *srf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt *srf2 = obj_int_[1]->getParamSurfaceInt();
    ParamSurface *s1=0, *s2=0;
//     int pdir = (dir < 2) ? dir : dir-2;
    int idx;
    bool modified;
    if (dir < 2)
	s1 = srf1->getParamSurface().get();
    else
	s2 = srf2->getParamSurface().get();

    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    for (int ki=0; ki<(int)(int_pts.size()); ki++)
    {
	modified = false;

	// Keep singular points
	SingularityType sing_type = int_pts[ki]->getSingularityType();
	if (sing_type == ISOLATED_POINT || sing_type == BRANCH_POINT ||
	    sing_type == HIGHER_ORDER_POINT) 
	    continue;

	// It is not necessarily wise to move a point lying in the
	// corner of one surface. Skip modificiation of such points
// 	if (srf1->inCorner(int_pts[ki]->getPar1(), ptol) ||
// 	    srf2->inCorner(int_pts[ki]->getPar2(), ptol))
// 	    continue;
	
	// Check if the point has a link along the current direction
	vector<shared_ptr<IntersectionLink> > links;
	int_pts[ki]->getNeighbourLinks(links);
	for (size_t kj=0; kj<links.size(); kj++)
	{
	    // Get the other point
	    IntersectionPoint *p1, *p2;
	    links[kj]->getIntersectionPoints(p1, p2);

	    // Check if the points are linked along the current direction
	    if (fabs(p1->getPar(dir) - p2->getPar(dir)) > ptol)
	      continue; 

// 	    // Check for degeneracy
// 	    double par[2];
// 	    par[pdir] = 0.5*(p1->getPar(dir) + p2->getPar(dir));
// 	    par[1-pdir] = (dir < 2) ? p1->getPar(1-dir) :  p1->getPar(5-dir);
// 	    bool degen = (dir < 2) ? 
// 		srf1->isDegenerate(epsge_->getEpsge(), dir, par)
// 		: srf2->isDegenerate(epsge_->getEpsge(), dir-2, par);
// 	    if (degen)
// 		continue;

	    // Check if both points lie in the current domain
	    if (!(int_results_->isInDomain(p1)) || 
		!(int_results_->isInDomain(p2)))
		continue;

	    // If the piece of the intersection code really leaves the
	    // parameter domain of the surface, keep it to avoid
	    // loosing a piece of the intersection curve later on
	    // Classify the points
	    IntPtClassification type1 = p1->getClassification(s1, s2);
	    IntPtClassification type2 = p2->getClassification(s1, s2);
	    if (type1 == DIR_TOUCH && (type2 == DIR_IN || type2 == DIR_OUT))
		continue;
	    if (type2 == DIR_TOUCH && (type1 == DIR_IN || type1 == DIR_OUT))
		continue;

	    // Iterate to a closest point between the intersection 
	    // points
	    double param[4];
	    double dist;
	    vector<double> par1 = p1->getPar();
	    vector<double> par2 = p2->getPar();
	    double d1 = p1->getDist();
	    double d2 = p2->getDist();
	    double d3, d4;
	    bool bd1, bd2, bd_mid, same_bd;
	    vector<bool> boundaries;

	    bd1 = atBoundary(p1->getPar1(), boundaries);
	    bd2 = atBoundary(p2->getPar1(), boundaries);
	    same_bd = atSameBoundary(p1->getPar1(), p2->getPar1());
	    if (bd1 && bd2 && !same_bd)
		continue;  // Must keep endpoint of link

	    if (keep_endpt && (bd1 || bd2) && (p2->numNeighbours() > 1 && p1->numNeighbours() > 1))
		continue;  // Must keep endpoints anyway

	    double fac = (d1 + d2 < epsge_->getRelParRes()) ? 0.5 
		: d2/(d1 + d2);
	    doIterate(dir, param, &par1[0], fac, &par2[0], 1.0-fac, dist);

	    // Check if a good intersection point is found
	    if (dist > p1->getDist()-ptol2 && 
		dist > p2->getDist()-ptol2)
		continue;
	    
	    bd_mid = atBoundary(param, boundaries);

	    // Check if the result of the iteration is identical to 
	    // one of the input points
	    shared_ptr<IntersectionPoint> tmp;
	    bool remove_tmp = false;
	    IntersectionPoint *curr = NULL;
	    for (idx=0; idx<4; idx++)
		if (fabs(param[idx]-par1[idx]) > ptol)
		    break;
	    if (idx == 4)
	    {
		curr = p1;
		dist = p1->getDist();
	    }
	    else
	    {
		for (idx=0; idx<4; idx++)
		    if (fabs(param[idx]-par2[idx]) > ptol)
			break;
		if (idx == 4)
		{
		    curr = p2;
		    dist = p2->getDist();
		}
		else
		{
		    // Make a new point
		    tmp = int_results_->addIntersectionPoint(obj_int_[0], 
							     obj_int_[1],
							     epsge_, 
							     param, param+2);
		    tmp->connectTo(p1, POST_ITERATE);
		    tmp->connectTo(p2, POST_ITERATE);
		    p1->disconnectFrom(p2);

		    // Set iso-parametric information
		    shared_ptr<IntersectionLink> lnk1 = tmp->getIntersectionLink(p1);
		    shared_ptr<IntersectionLink> lnk2 = tmp->getIntersectionLink(p2);
		    lnk1->setIsoparametricIn(dir, true);
		    lnk2->setIsoparametricIn(dir, true);

		    curr = tmp.get();
		    d3 = p1->getPoint().dist(tmp->getPoint());
		    d4 = p2->getPoint().dist(tmp->getPoint());
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfSfIntersector::postIterate. New intersection point at: ";
			std::cout << param[0] << " " << param[1] << " " << param[2] << " " << param[3];
			std::cout << " Dist: " << tmp->getDist() << std::endl;
			std::cout << "Between " << par1[0] << " " << par1[1] << " " << par1[2];
			std::cout << " " << par1[3] << " " << d1 << " and ";
			std::cout << par2[0] << " " << par2[1] << " " << par2[2];
			std::cout << " " << par2[3] << " " << d2 << std::endl;
		    }

		}
	    }
	    
	    if (dist < p1->getDist() && p1->numNeighbours() <= 2 &&
		!(bd1 && !bd_mid) && !srf1->inCorner(p1->getPar1(), ptol) &&
		!srf2->inCorner(p1->getPar2(), ptol))
	    {
		// Remove p1. First find the index in the current array
		size_t kr;
		for (kr=0; kr<int_pts.size(); kr++)
		    if (int_pts[kr].get() == p1)
			break;

		if (kr < int_pts.size())
		{
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfSfIntersector::postIterate. Removed point at: ";
			std::cout << par1[0] << " " << par1[1] << " " << par1[2] << " " << par1[3];
			std::cout << std::endl;
		    }
		    int_results_->removeIntPoint(int_pts[kr]);
		    int_pts.erase(int_pts.begin()+kr);
		    ki--;
		    modified = true;
		}
	    }
	    else if (tmp.get() && p1->getPoint().dist(tmp->getPoint()) < tol2 /*epsge_->getEpsge()*/)
	    {
		// Remove the new point again to avoid noice
		remove_tmp = true;
	    }

	    if (dist < p2->getDist() && p2->numNeighbours() <= 2 &&
		!(bd2 && !bd_mid) && !srf1->inCorner(p2->getPar1(), ptol) &&
		!srf2->inCorner(p2->getPar2(), ptol))
	    {
		// Remove p1. First find the index in the current array
		size_t kr;
		for (kr=0; kr<int_pts.size(); kr++)
		    if (int_pts[kr].get() == p2)
			break;

		if (kr < int_pts.size())
		{
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfSfIntersector::postIterate. Removed point at: ";
			std::cout << par2[0] << " " << par2[1] << " " << par2[2] << " " << par2[3];
			std::cout << std::endl;
		    }
		    int_results_->removeIntPoint(int_pts[kr]);
		    int_pts.erase(int_pts.begin()+kr);
		    ki--;
		    modified = true;
		}
	    }
	    else if (tmp.get() && p2->getPoint().dist(tmp->getPoint()) < tol2 /*epsge_->getEpsge()*/)
	    {
		// Remove the new point again to avoid noice
		remove_tmp = true;
	    }
	    if (remove_tmp)
		int_results_->removeIntPoint(tmp);	

	    if (modified)
		break;
	}

	ki = std::max(ki, -1);
    }

}
 

//===========================================================================
void SfSfIntersector::postIterate2(int nmb_orig, int dir, double par, 
				   bool along, bool only_isolated)
//===========================================================================
{
    // Check the quality of all new interserction points. Iterate
    // along the specified parameter direction.
    //
    // If along == true then par is expected to correspond to the
    // parameter of the intersection points and it is legal to move
    // the intersection point along this line. Otherwise the post
    // iteration serves as a check for double representation of
    // points. The point cannot be moved.
    //
    // If only_isolated == true then only isolated intersection points
    // are allowed to be iterated. If postIterate2() is called within
    // a subdivision, only_isolated will be true, while if it is
    // called within updateIntersections(), only_isolated will be
    // false.

//     double ptol = epsge_->getRelParRes();
    ParamSurfaceInt *srf1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt *srf2 = obj_int_[1]->getParamSurfaceInt();
    double tol2 = 0.0;

    vector<double> mima1 = srf1->getMima();
    vector<double> mima2 = srf2->getMima();
    double limit1[4], limit2[4];
    for (int idx=0; idx<2; idx++) {
	limit1[idx] = mima1[2*idx];
	limit1[2+idx] = mima2[2*idx];
	limit2[idx] = mima1[2*idx+1];
	limit2[2+idx] = mima2[2*idx+1];
    }
    if (along) {
	limit1[dir] = limit2[dir] = par;
    }

    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    // Make limit on when to remove a new point if one of the initial
    // points are singular. A singular intersection point can have
    // lower accuracy than a transversal one and still be the best choice.
    for (int idx=0; idx<nmb_orig; idx++) {
	if (int_pts[idx]->pointIsSingular()) {
	    tol2 = std::max(tol2, int_pts[idx]->getDist());
	}
    }

    for (int ki=nmb_orig; ki<(int)(int_pts.size()); ki++) {
	if (only_isolated && (int_pts[ki]->numNeighbours() > 0)) {
	    continue;  // Only isolated points are handled
	}

	// Keep singular points
	SingularityType sing_type = int_pts[ki]->getSingularityType();
	if (sing_type == ISOLATED_POINT || sing_type == BRANCH_POINT ||
	    sing_type == HIGHER_ORDER_POINT) 
	    continue;

	// It is not necessarily wise to move a point lying in the
	// corner of one surface. Skip modificiation of such points
// 	if (srf1->inCorner(int_pts[ki]->getPar1(), ptol) ||
// 	    srf2->inCorner(int_pts[ki]->getPar2(), ptol))
// 	    continue;

	if (!along) {
	    limit1[dir] = limit2[dir] = int_pts[ki]->getPar(dir);
	}
	double param[4];
	double dist;
	vector<double> parval = int_pts[ki]->getPar();
	doIterate(dir, param, limit1, 0.5, limit2, 0.5, dist, &parval[0]);
	if (dist < int_pts[ki]->getDist()) {
	    // A better intersection point is found
	    // Check if it exists already
	    shared_ptr<IntersectionPoint> closest;
	    bool exist = int_results_->closestInDomain(param, closest);
	    if (closest.get() == int_pts[ki].get())
		exist = false;

	    if (int_pts[ki]->getDist() <= tol2)
		exist = false;  // Do not remove this point

	    // Check if the closest point can replace the current one
	    if (exist && along)
	    {
		double pclose = closest->getPar(dir);
		if (fabs(pclose - par) > epsge_->getRelParRes())
		    exist = false;
	    }

	    if (exist) {
		// Check influence area
		int kj;
		for (kj=0; kj<4; kj++) {
		    if (kj == dir) {
			continue;
		    }
		    if (closest->inInfluenceArea(kj, param[kj], false)) {
			break;
		    }
		}
		if (kj < 4) {
		    // The point is a double representation. Remove it.
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfSfIntersector::postIterate2. Removed point at: ";
			std::cout << parval[0] << " " << parval[1] << " " << parval[2] << " " << parval[3];
			std::cout << std::endl;
		    }
		    int_results_->removeIntPoint(int_pts[ki]);
		    int_pts.erase(int_pts.begin()+ki);
		    ki--;
		} else {
		    exist = false;
		}
	    }
	    if (along && (!exist)) {
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfSfIntersector::postIterate2. Moved point at: ";
			std::cout << parval[0] << " " << parval[1] << " " << parval[2] << " " << parval[3];
			std::cout << " " << int_pts[ki]->getDist();
			std::cout << " to " << param[0] << " " << param[1] << " " << param[2] << " ";
			std::cout << param[3] << " " << dist;
			std::cout << std::endl;
		    }
		// Modify parameter of intersection point
		int_pts[ki]->replaceParameter(param); //, 4);
		
	    } else if ((!only_isolated) && (!exist)) {
		// We also move the point if we are in the
		// "only_isolated mode"
		//int_pts[ki]->replaceParameter(param);
	    }
	}
    }
}


//===========================================================================
void SfSfIntersector::
doIterate(int pardir, double param[], double limit1[],
	  double w1, double limit2[], double w2,
	  double& dist, double *seed)
//===========================================================================
{
    // Purpose: Iterate to a closest point between a curve in one
    // surface and the other surface.

    // Make curve and fetch surface
    Point pnt;
    shared_ptr<ParamSurface> surf;
    shared_ptr<ParamCurve> crv;
    RectDomain domain;  
    double cstart, cend;
    int sf_idx, cv_idx;
    //double pt_par[2];
   double guess[4];
   double parval = 0.5*(limit1[pardir] + limit2[pardir]);
   if (pardir < 2)
    {
	crv = obj_int_[0]->getParamSurfaceInt()->
	    getConstantParameterCurve(pardir, parval);
	surf = obj_int_[1]->getParamSurfaceInt()->getParentParamSurface(domain);
	Vector2D lower(std::min(limit1[2],limit2[2]),std::min(limit1[3],limit2[3]));
	Vector2D upper(std::max(limit1[2],limit2[2]),std::max(limit1[3],limit2[3]));
	RectDomain dom(lower,upper);
	domain = dom;
	cv_idx = 1 - pardir;
	sf_idx = 2;
    }
    else
    {
	crv = obj_int_[1]->getParamSurfaceInt()->
	    getConstantParameterCurve(pardir-2, parval);
	surf = obj_int_[0]->getParamSurfaceInt()->getParentParamSurface(domain);
	Vector2D lower(std::min(limit1[0],limit2[0]),std::min(limit1[1],limit2[1]));
	Vector2D upper(std::max(limit1[0],limit2[0]),std::max(limit1[1],limit2[1]));
	RectDomain dom(lower,upper);
	domain = dom;
	cv_idx = 5 - pardir;
	sf_idx = 0;
   }
   cstart = std::min(limit1[cv_idx], limit2[cv_idx]);
   cend = std::max(limit1[cv_idx], limit2[cv_idx]);

   for (int ki=0; ki<4; ki++)
       guess[ki] = (seed) ? seed[ki] : w1*limit1[ki] + w2*limit2[ki];

   // Iterate
   Point pt_cv, pt_sf;
   param[pardir] = parval;
   ClosestPoint::closestPtCurveSurf(crv.get(), surf.get(), epsge_->getEpsge(), cstart, 
		      cend, &domain, guess[cv_idx], guess+sf_idx, 
		      param[cv_idx], param+sf_idx, dist, pt_cv, pt_sf);
		      
}


//===========================================================================
void SfSfIntersector::doIterate(int pardir, double parval, double param[], 
				double& dist, double seed[])
//===========================================================================
{
    // Purpose: Iterate to a closest point between a point in one
    // surface and the other surface.

    // Make point and fetch surface
    Point pnt;
    shared_ptr<ParamSurface> surf;
    RectDomain domain;  
    int sf_idx;
    double pt_par[2];
   double guess[2];
   if (pardir < 2)
    {
	pt_par[pardir] = parval;
	pt_par[1-pardir] = seed[1-pardir];
	obj_int_[0]->point(pnt, pt_par);
	surf = obj_int_[1]->getParamSurfaceInt()->getParentParamSurface(domain);
	sf_idx = 2;
    }
    else
    {
	pt_par[pardir-2] = parval;
	pt_par[3-pardir] = seed[5-pardir];
	obj_int_[1]->point(pnt, pt_par);
	surf = obj_int_[0]->getParamSurfaceInt()->getParentParamSurface(domain);
	sf_idx = 0;
   }

    // Modify seed 
   int ki, kj;
   for (ki=sf_idx, kj=0; kj<2; ki++)
       guess[kj++] = seed[ki];

    // Iterate
  Point ptc;
  for (ki=0; ki<4; ki++)
      param[ki] = seed[ki];
  param[pardir] = parval;
  surf->closestPoint(pnt, param[sf_idx], param[sf_idx+1], ptc, dist,
		     epsge_->getEpsge(), &domain, guess);
}


//===========================================================================
void SfSfIntersector::removeDegenerateConnections()
//===========================================================================
{
    //  Purpose: Avoid false connections in the degenerate case. If
    //  two connected intersection points have the same parameter
    //  value on one of the objects, remove the connection

    int nmbpar1 = obj_int_[0]->numParams();
    int nmbpar2 = obj_int_[1]->numParams();
    double ptol = epsge_->getRelParRes();

    // Fetch intersection points
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    for (size_t ki=0; ki<int_pts.size(); ki++)
    {
	// Get connections
	vector<shared_ptr<IntersectionLink> > links;
	int_pts[ki]->getNeighbourLinks(links);

	for (size_t kj=0; kj<links.size(); kj++)
	{
	    // Check equality of parameter values of the connected points.
	    // First fetch points.
	    IntersectionPoint *p1, *p2;
	    links[kj]->getIntersectionPoints(p1, p2);

	    const double *par1_1 = p1->getPar1();
	    const double *par1_2 = p1->getPar2();
	    const double *par2_1 = p2->getPar1();
	    const double *par2_2 = p2->getPar2();

	    int kr1, kr2;
	    for (kr1=0; kr1<nmbpar1; kr1++)
		if (fabs(par1_1[kr1]-par2_1[kr1]) > ptol)
		    break;
	    for (kr2=0; kr2<nmbpar2; kr2++)
		if (fabs(par1_2[kr2]-par2_2[kr2]) > ptol)
		    break;

	    if (kr1 >= nmbpar1 && kr2 >= nmbpar1)
	    {
		// Identical endpoints in both objects. Remove link
		p1->disconnectFrom(p2);
	    }
		
	    else if (kr1 >= nmbpar1)
	    {
		// Identical point in the first object
		// Make an exstra check for the midpoint of the link in 
		// the case of an iso parametric link in the other object
		int kh, idx;
		for (kh=0; kh<nmbpar2; kh++)
		    if (links[kj]->isIsoparametricIn(nmbpar1+kh))
			break;

		if (kh < nmbpar2)
		{
		    // Get endpoints
		    Point pos1 = p1->getPoint();
		    Point pos2 = p2->getPoint();
		    Point posmid;
		    double param[2];  // At most two parameter direction
		    for (idx=0; idx<nmbpar2; idx++)
			param[idx] = 0.5*(par1_2[idx] + par2_2[idx]);
		    obj_int_[1]->point(posmid, param);
		    //if (pos1.dist(posmid) + posmid.dist(pos2) < epsge_->getEpsge())
		    if (true)
		    {
			// Remove link
			p1->disconnectFrom(p2);
		    }
		}
		else
		{
		    // Remove link
		    p1->disconnectFrom(p2);
		}
	    }
	    else if (kr2 >= nmbpar2)
	    {
		// Identical point in the second object
		// Make an exstra check for the midpoint of the link in 
		// the case of an iso parametric link
		int kh, idx;
		for (kh=0; kh<nmbpar1; kh++)
		    if (links[kj]->isIsoparametricIn(kh))
			break;

		if (kh < nmbpar1)
		{
		    // Get endpoints
		    Point pos1 = p1->getPoint();
		    Point pos2 = p2->getPoint();
		    Point posmid;
		    double param[2];  // At most two parameter direction
		    for (idx=0; idx<nmbpar1; idx++)
			param[idx] = 0.5*(par1_1[idx] + par2_1[idx]);
		    obj_int_[0]->point(posmid, param);
		    //if (pos1.dist(posmid) + posmid.dist(pos2) < epsge_->getEpsge())
		    if (true)
		    {
			// Remove link
			p1->disconnectFrom(p2);
		    }
		}
		else
		{
		    // Remove link
		    p1->disconnectFrom(p2);
		}
	    }
	    
	}
    }
}


//===========================================================================
bool SfSfIntersector::
getSubdivParDegEdge(ParamSurfaceInt *surf, int dir, int pdir, 
		    int deg_edge, double treshhold, double& par)
//===========================================================================
{
    int ki;

    // Fetch suitable subdivision parameter from the object
    par = surf->isolateDegPar(pdir, deg_edge, treshhold);

    // Modify parameter according to boundary intersections at the
    // degenerate boundary. The parameter must lie outside the influence
    // area of the intersection points unless there is an intersection
    // curve along the adjacent boundaries.
    double bd_par = (deg_edge == 1 || deg_edge == 3) ? surf->startParam(pdir) 
	: surf->endParam(pdir);
    double tdel = surf->endParam(pdir) - surf->startParam(pdir);
    double par2=0.0, curr;
    double tol = 2.0*epsge_->getEpsge();
    bool forward = (deg_edge == 1 || deg_edge == 3);
    vector<shared_ptr<IntersectionPoint> > bd_ints;
    int_results_->getIntersectionPoints(dir, bd_par, bd_ints);
    for (ki=0; ki<int(bd_ints.size()); ki++)
    {
	curr = bd_ints[ki]->getInfluenceArea(dir, forward, true, tol);
	par2 = std::max(par2, fabs(curr-bd_par));
    }
    //par2 = 1.1*par2;  // Just to be sure to be outside the influence area
    //par2 = 5.0*par2;  // Just to be sure to be outside the influence area
    // Control size
    if (par2 > 100.0*fabs(par-bd_par))
      par2 = 100.0*(fabs(par-bd_par));  
    if (par2 > 0.1*tdel)
      par2 = 0.1*tdel;  
    par2 = std::max(par2, 1.0e-6);
    par = (deg_edge == 1 || deg_edge == 3) ? std::max(bd_par+par2, par)
	: std::min(bd_par-par2, par);

    // VSK, 0906. Check nearby intersection points to avoid dividing in a very 
    // close vicinity
    double ptol = 0.1*fabs(bd_par-par);
    vector<shared_ptr<IntersectionPoint> > all_ints;
    int_results_->getIntersectionPoints(all_ints);
    for (ki=0; ki<(int)(all_ints.size()); ki++)
    {
	if (fabs(all_ints[ki]->getPar(pdir) - par) < ptol)
	{
	    par = all_ints[ki]->getPar(pdir);
	    break;
	}
    }
    return true;
}


//===========================================================================
void SfSfIntersector::
setDegTriangle(vector<shared_ptr<ParamGeomInt> >& sub_objects,
	       int deg_edge, int pdir, double par)
//===========================================================================
{
    // deg_edge =0 not degenerate, =1 degenerate at start, = 2
    // degenerate at endpoint, =3 degenerate in both endpoints

    size_t ki;
    double ta;
    if (deg_edge == 1 || deg_edge == 3)
    {
	// The leftmost/lower sub surface is a tiny triangle.
	// Find these sub surfaces
	for (ki=0; ki<sub_objects.size(); ki++)
	{
	    ta = sub_objects[ki]->startParam(pdir);
	    if (ta < par - epsge_->getRelParRes())
	    {
		ParamSurfaceInt *tmp = (ParamSurfaceInt*)(sub_objects[ki].get());
		tmp->setDegTriang();
	    }
	}
    }
    else if (deg_edge == 2)
    {
	// The rightmost/upper sub surface is a tiny triangle.
	// Find these sub surfaces
	for (ki=0; ki<sub_objects.size(); ki++)
	{
	    ta = sub_objects[ki]->endParam(pdir);
	    if (ta > par +  epsge_->getRelParRes())
	    {
		ParamSurfaceInt *tmp = (ParamSurfaceInt*)(sub_objects[ki].get());
		tmp->setDegTriang();
	    }
	}
    }
}


//===========================================================================
int SfSfIntersector::
checkSubdivParam(int dir, double par, double ta, double tb,
		 vector<shared_ptr<IntersectionPoint> >& int_pts)
//===========================================================================
{
    // Return value : = 1 : Not a good parameter
    //                = 2 : Parameter accepted

    int pdir = (dir < 2) ? dir : dir-2;
    ParamSurfaceInt* surf_int = (dir < 2) ? obj_int_[0]->getParamSurfaceInt()
	: obj_int_[1]->getParamSurfaceInt();
    int ki;
    double pval;
    double frac = 0.2;
    double ang_tol = std::min(10.0*epsge_->getRefAng(), 0.1);
    double eps = epsge_->getEpsge();

    // Skip the intersection points at a safe distance from the
    // subdivision value
    int npts = (int)int_pts.size();
    for (ki=0; ki<npts; ki++) {
	pval = int_pts[ki]->getPar(dir);
	if (fabs(par - pval) >= frac*(tb - ta) &&
	    !(int_pts[ki]->inInfluenceAreaBracket(dir, pval))) {
	    std::swap(int_pts[ki], int_pts[npts-1]);
	    npts--;
	    ki--;
	}
    }

    //     if (npts > 2)
    // 	return 1;   // Not a good subdivision parameter

    if (npts > 1) {
	// Fetch the first and the last point
	int idx = (dir >= 2)*2 + (1-pdir);
	vector<shared_ptr<IntersectionPoint> > pts(2);
	pts[0] = pts[1] = int_pts[0];
	for (ki=1; ki<npts; ki++) {
	    if (int_pts[ki]->getPar(idx) < pts[0]->getPar(idx))
		pts[0] = int_pts[ki];
	    if (int_pts[ki]->getPar(idx) > pts[1]->getPar(idx))
		pts[1] = int_pts[ki];
	}

	// Check if the complete subdivision curve is in intersection.
	// First check if the two intersection points lie at the same
	// boundary.  In that case, the subdivision parameter cannot
	// be accepted
// 	if (int_results_->atSameBoundary(pts[0], pts[1]))
// 	    return 1;

// 	// Check the sequence of the intersectionpoint
// 	if (int_pts[0]->getPar(dir) > int_pts[1]->getPar(dir))
// 	    std::swap(int_pts[0], int_pts[1]);

	int coincide; 
	// Check the risk of a double connection
	if (pts[0]->isConnectedTo(pts[1]) && 
	    (fabs(pts[0]->getPar(dir)-par) > epsge_->getRelParRes() ||
	     fabs(pts[1]->getPar(dir)-par) > epsge_->getRelParRes()))
	    coincide = 1;
	else
	    coincide= checkIsoCurve(pdir, (pdir == dir), par, pts);
	if (coincide == 2)
	    return 2;    // Constant parameter curve OK
    }

    // For all (one or two) intersection point, check that the angle
    // between the direction of the candidate subdivision curve and
    // the intersection curve is not too acute.
    Point pos2, der_u, der_v;
    double param[2];
    param[pdir] = par;
    for (ki=0; ki<npts; ki++) {
	if (int_pts[ki]->getSingularityType() == HIGHER_ORDER_POINT ||
	    int_pts[ki]->getSingularityType() == ISOLATED_POINT)
	    continue;
	if (int_pts[ki]->getSingularityType() == BRANCH_POINT)
	    continue;  // We really want to subdivide here. TEST for the time being
	param[1-pdir] = (dir == pdir) ? int_pts[ki]->getPar(1-pdir) 
	    : int_pts[ki]->getPar(3-pdir);
	surf_int->point(pos2, param);
	surf_int->derivs(param[0], param[1], der_u, der_v);
	Point pos1 = int_pts[ki]->getPoint();
	Point tang1 = int_pts[ki]->getTangent();
	double alpha = tang1.angle((pdir == 0) ? der_v : der_u);
	if (pos1.dist(pos2) > eps) {
	    Point vec = pos1 - pos2;
	    alpha += vec.angle_smallest((pdir == 0) ? der_v : der_u);
	}
	if (alpha < ang_tol || fabs(M_PI-alpha) < ang_tol)
	    return 1;

	if (int_pts[ki]->getSingularityType() == BRANCH_POINT) {
	    tang1 = int_pts[ki]->getTangent(true);
	    alpha = tang1.angle((pdir == 0) ? der_v : der_u);
	    if (pos1.dist(pos2) > eps) {
		Point vec = pos1 - pos2;
		alpha += vec.angle_smallest((pdir == 0) ? der_v : der_u);
	    }
	    if (alpha < ang_tol || fabs(M_PI-alpha) < ang_tol)
		return 1;
	}
    }

    return 2;
}


//===========================================================================
int SfSfIntersector::
checkIsoCurve(int pdir, bool first, double par,
	      vector<shared_ptr<IntersectionPoint> >& int_pts)
//===========================================================================
{
    ParamSurfaceInt* surf_int1
	= (first)
	? obj_int_[0]->getParamSurfaceInt()
	: obj_int_[1]->getParamSurfaceInt();
    ParamSurfaceInt* surf_int2
	= (first)
	? obj_int_[1]->getParamSurfaceInt()
	: obj_int_[0]->getParamSurfaceInt();
    shared_ptr<ParamSurface> surf1 = surf_int1->getParamSurface();
    shared_ptr<ParamSurface> surf2 = surf_int2->getParamSurface();
    int ki, kj;

    // Represent the candidate subdivision curve as curve-in-surface
    // curve
    Point parval[2];
    double seed[2][2];
    for (ki=0; ki<2; ki++) {
	if (first) {
	    parval[ki].setValue(int_pts[ki]->getPar1()[0], 
				int_pts[ki]->getPar1()[1]);
	    for (kj=0; kj<2; kj++) {
		seed[ki][kj] = int_pts[ki]->getPar(2+kj);
	    }
	}
	else {
	    parval[ki].setValue(int_pts[ki]->getPar2()[0], 
				int_pts[ki]->getPar2()[1]);
	    for (kj=0; kj<2; kj++) {
		seed[ki][kj] = int_pts[ki]->getPar(kj);
	    }
	}
    }

    double ta = parval[0][1-pdir];
    double tb = parval[1][1-pdir];
    double start
	= (first)
	? obj_int_[0]->startParam(1-pdir)
	: obj_int_[1]->startParam(1-pdir);
    double end
	= (first)
	? obj_int_[0]->endParam(1-pdir)
	: obj_int_[1]->endParam(1-pdir);
    if (std::min(ta, tb) > start + epsge_->getRelParRes()
	|| std::max(ta, tb) < end - epsge_->getRelParRes()) {
	// Check if the two intersection points lie at different
	// boundaries in the other surface.
	if (!(int_results_->atDifferentBoundary(int_pts[0], int_pts[1]))) {
	    return 1;
	}
    }

    if (pdir == 0) {
	parval[0].setValue(par, ta);
	parval[1].setValue(par, tb);
    }
    else {
	parval[0].setValue(ta, par);
	parval[1].setValue(tb, par);
    }
    shared_ptr<SplineCurve> pcrv;
    if (ta < tb) {
	pcrv = shared_ptr<SplineCurve>(new SplineCurve(parval[0], ta,
						       parval[1], tb));
    }
    else {
	pcrv = shared_ptr<SplineCurve>(new SplineCurve(parval[1], tb,
						       parval[0], ta));
    }
    shared_ptr<ParamCurve> constcrv(new CurveOnSurface(surf1, pcrv, true));

    // Project the endpoints of the subdivision curve onto the other
    // surface
    Point pt1, pt2;
    constcrv->point(pt1, constcrv->startparam());
    constcrv->point(pt2, constcrv->endparam());

    double clo_u1, clo_u2, clo_v1, clo_v2, clo_d1, clo_d2;
    Point clo_pt1, clo_pt2;
    surf2->closestPoint(pt1, clo_u1, clo_v1, clo_pt1, clo_d1,
			epsge_->getEpsge(), NULL, seed[0]);
    if (clo_d1 > epsge_->getEpsge()) {
	return 1;
    }
    surf2->closestPoint(pt2, clo_u2, clo_v2, clo_pt2, clo_d2,
			epsge_->getEpsge(), NULL, seed[1]);
    if (clo_d2 > epsge_->getEpsge()) {
	return 1;
    }
    Point su_start(clo_u1, clo_v1), su_end(clo_u2, clo_v2);
    shared_ptr<ParamCurveInt> curve_int(new ParamCurveInt(constcrv));

    int coincide = checkCoincide(curve_int.get(), constcrv->startparam(),
				 constcrv->endparam(), surf_int2,
				 su_start, su_end, epsge_);
    if (coincide) {
	return 2;
    }

    return 0;
}


//===========================================================================
bool SfSfIntersector::
getSubdivAtSing(int dir, double ta, double tb, double& par)
//===========================================================================
{
    // Purpose: Search for a subdivision parameter in a singular
    // intersection point between the two surfaces, the point is
    // either already found at the surface boudnaries or an iteration
    // will be performed to look for a singularity.

    double frac = 0.01;

    // Fetch singularity information from the previous intersector if
    // necessary and if such an information exist
    if (singularity_info_.get() == 0 && prev_intersector_ &&
	prev_intersector_->numParams() == 4 && 
	prev_intersector_->hasSingularityInfo())
    {
	singularity_info_ = (shared_ptr<SingularityInfo>)
	    (new SingularityInfo(prev_intersector_->getSingularityInfo()));
    }
    else if (singularity_info_.get() == 0)
    {
	// Make empty singularity info instance
	singularity_info_ = (shared_ptr<SingularityInfo>)(new SingularityInfo());
    }

    // First fetch all singular intersection points (branch points) at the
    // surface boundaries
    bool is_critical;
    int subdiv_ok = 0;
    vector<shared_ptr<IntersectionPoint> > bd_sing;
    int_results_->getBranchPoints(bd_sing);
    int ki;
    for (ki=0; ki<int(bd_sing.size()); ki++)
    {
	// Check the current singularity to see if it is a suitable
	// subdivision point
	par = bd_sing[ki]->getPar(dir);
	if (par < ta+frac*(tb-ta) || par > tb-frac*(tb-ta))
	    continue;
	vector<shared_ptr<IntersectionPoint> > int_pts;
	is_critical = (int_results_->inInfluenceArea(dir, par, int_pts) != 0);
	if (int_pts.size() > 0)
	{
	    // Check if the subdivision position is OK with respect to the
	    // angle between the subdivision curve and the intersection curve
	    subdiv_ok = checkSubdivParam(dir, par, ta, tb, int_pts);
	}
	if (subdiv_ok == 0 || subdiv_ok == 2)
	    return true;
    }

    // No suitable singularity at the boundary is found. Check the singularity
    // information from the intersector.
    if (singularity_info_->hasPoint())
    {
	par = singularity_info_->getParam(dir);
	double ta = (dir < 2) ? obj_int_[0]->startParam(dir) 
	    : obj_int_[1]->startParam(dir-2);
 	double tb = (dir < 2) ? obj_int_[0]->endParam(dir) 
	    : obj_int_[1]->endParam(dir-2);
	if (par > ta+frac*(tb-ta) && par < tb-frac*(tb-ta))
	{
	    // Check the parameter value
	    vector<shared_ptr<IntersectionPoint> > int_pts;
	    is_critical = (int_results_->inInfluenceArea(dir, par, int_pts) != 0);
	    if (int_pts.size() > 0)
	    {
		// Check if the subdivision position is OK with respect to the
		// angle between the subdivision curve and the intersection curve
		subdiv_ok = checkSubdivParam(dir, par, ta, tb, int_pts);
	    }
	    if (subdiv_ok == 0 || subdiv_ok == 2)
	    {
		int is_OK = splitIntResults(int_pts, dir, par, ta, tb, false);
		if (is_OK)
		    return true;
	    }
	}
    }

    
    // Iterate for a singularity. First check if the current surfaces are 
    // simple enough.
    // Make surface pointers
    ParamSurfaceInt* surf_int1 = obj_int_[0]->getParamSurfaceInt();
    ParamSurfaceInt* surf_int2 = obj_int_[1]->getParamSurfaceInt();

    bool simple1 = surf_int1->isSimple();
    bool simple2 = surf_int2->isSimple();
    int bound_simple = 5;
    if ((!singularity_info_->iterationDone() && (simple1 && simple2))
	|| singularity_info_->nmbSimple1() > bound_simple ||
	singularity_info_->nmbSimple2() > bound_simple)
    {
	// Iterate for a singular intersection point
	
 	RectDomain domain1, domain2;
	shared_ptr<ParamSurface> surf1 = 
	    surf_int1->getParentParamSurface(domain1);
	shared_ptr<ParamSurface> surf2 = 
	    surf_int2->getParentParamSurface(domain2);

	double angle_tol = epsge_->getAngleTol();
	double seed[4];
	double result_par[4];
	int constraints[2];
	double constraints_par[2];
	double limit[8];
	limit[0] = domain1.umin();
	limit[1] = domain1.umax();
	limit[2] = domain1.vmin();
	limit[3] = domain1.vmax();
	limit[4] = domain2.umin();
	limit[5] = domain2.umax();
	limit[6] = domain2.vmin();
	limit[7] = domain2.vmax();

	// Get seed. If intersection points at the boundaries exist,
	// use the middle point between the intersections
	if (int_results_->hasIntersectionPoints())
	{
	    int_results_->getMidParameter(seed);
	}
	else
	{
	    // Use the midpoint
	    for (ki=0; ki<4; ki++)
		seed[ki] = 0.5*(limit[2*ki] + limit[2*ki+1]);
	}
	constraints[0] = constraints[1] = -1;
	constraints_par[0] = constraints_par[1] = 0.0;

	// Iterate
	int isfound;
	isfound = extremalPtSurfSurf(surf1.get(), surf2.get(), constraints, 
				     constraints_par, limit, seed, 
				     result_par, angle_tol);

	// Check if the found point is also an intersection point
	Point pt1, pt2;
	obj_int_[0]->point(pt1, result_par);
	obj_int_[1]->point(pt2, result_par+2);
	double dist = pt1.dist(pt2);
	if (dist >= epsge_->getEpsge())
	    isfound = false;

	// Update singularity information
	singularity_info_->addIterationCount();

	// Check output
	if (isfound)
	{
	    // Create a temporary intersection point to check the
	    // singularity type of the point. Tangential intersection
	    // points will not be stored as singularities for
	    // subdivision
	    shared_ptr<IntersectionPoint> tmp = 
		int_results_->addIntersectionPoint(obj_int_[0], 
						   obj_int_[1], 
						   epsge_,
						   result_par, 
						   result_par + obj_int_[0]->numParams());

	    // Check if the found singularity is a branch point
	    SingularityType sing_type = tmp->getSingularityType();

	    // Remove point
	    int_results_->removeIntPoint(tmp);

	    
	    if (sing_type == ISOLATED_POINT || sing_type == BRANCH_POINT ||
		sing_type == HIGHER_ORDER_POINT)
	    {
		// Store singularity info
		singularity_info_->setSingularPoint(result_par, 4);

		// Normal quality check for subdivision parameters
		par = result_par[dir];
		if (par >= ta+frac*(tb-ta) && par <= tb-frac*(tb-ta))
		{
		    vector<shared_ptr<IntersectionPoint> > int_pts;
		    is_critical = (int_results_->inInfluenceArea(dir, par, int_pts) != 0);
		    if (int_pts.size() > 0)
		    {
			// Check if the subdivision position is OK with respect to the
			// angle between the subdivision curve and the intersection curve
			subdiv_ok = checkSubdivParam(dir, par, ta, tb, int_pts);
		    }
		    if (subdiv_ok == 0 || subdiv_ok == 2)
		    {
			int is_OK = splitIntResults(int_pts, dir, par, ta, tb, false);
			if (is_OK)
			    return true;
		    }
		}
	    }
	}
    }
    singularity_info_->addSimpleCount(simple1, simple2);
	    
    return false;  // No singularity suitable for subdivision is found
}


//===========================================================================
void SfSfIntersector::createApproxImplicit(double tol, int recursion_limit)
//===========================================================================
{
    // Note: Default is recursion_limit = 10. To turn off this
    // limitation, use recursion_limit = 0.

    // Check if the surfaces are splines
    shared_ptr<SplineSurfaceInt> spsf_int1 =
	dynamic_pointer_cast<SplineSurfaceInt, ParamObjectInt>(obj_int_[0]);
    shared_ptr<SplineSurfaceInt> spsf_int2 =
	dynamic_pointer_cast<SplineSurfaceInt, ParamObjectInt>(obj_int_[1]);
    if(spsf_int1.get() == 0 || spsf_int2.get() == 0) {
	cout << "At least one surface is not a spline" << endl;
	return;
    }

    // Check if the surfaces are Bezier or the number of recursions is
    // large enough
    shared_ptr<SplineSurface> sf1 =
	dynamic_pointer_cast<SplineSurface, ParamSurface>
	(spsf_int1->getParamSurface());
    shared_ptr<SplineSurface> sf2 =
	dynamic_pointer_cast<SplineSurface, ParamSurface>
	(spsf_int2->getParamSurface());
    bool first_is_bezier = (sf1->numCoefs_u() == sf1->order_u()
			    && sf1->numCoefs_v() == sf1->order_v());
    bool second_is_bezier = (sf2->numCoefs_u() == sf2->order_u()
			     && sf2->numCoefs_v() == sf2->order_v());
    bool recursion_limit_reached = (nmbRecursions() >= recursion_limit);
    if ((!first_is_bezier || !second_is_bezier) && !recursion_limit_reached) {
	cout << "At least one surface is not Bezier, but the recursion "
	     << "limit is not reached." << endl
	     << "We don't implicitize." << endl;
	return;
    }
    if ((!first_is_bezier || !second_is_bezier) && recursion_limit_reached) {
	cout << "*** Implicitizing a spline... Recursion level = "
	     << nmbRecursions() << endl;
    }

    // Check if we need to implicitize the first surface
    bool first_exists = (approx_implicit_[0][0].get() != 0
			 && approx_implicit_[0][1].get() != 0);
    bool first_error_good_enough = (approx_implicit_err_[0] >= 0.0
				    && approx_implicit_err_[0] <= tol);
    bool implicitize_first = (!first_exists || (prev_implicit_[0] && !first_error_good_enough));

    // Check if we need to implicitize the second surface
    bool second_exists = (approx_implicit_[1][0].get() != 0
			  && approx_implicit_[1][1].get() != 0);
    bool second_error_good_enough = (approx_implicit_err_[1] >= 0.0
				     && approx_implicit_err_[1] <= tol);
    bool implicitize_second = (!second_exists || (prev_implicit_[1] && !second_error_good_enough));

//     // TESTING @@@jbt
//     implicitize_first = true;
//     implicitize_second = true;

    // Return if no implicitization is needed
    if (!implicitize_first && !implicitize_second) {
	cout << "No implicitization is needed" << endl;
	return;
    }

    AlgObj3DInt alg_obj_3d_int(-1);
    double dummy_tol = 0.0; // Not in use. @@@jbt
    double approx_implicit_err = 0.0;
    BernsteinTetrahedralPoly impl;
    BaryCoordSystem3D bc;

    // Implicitize first
    if (implicitize_first) {
	spsf_int1->implicitize(dummy_tol);
	bool success;
	success = spsf_int1->getImplicit(dummy_tol,
					 approx_implicit_err,
					 alg_obj_3d_int);
	ASSERT(success);
	if (approx_implicit_err < fabs(approx_implicit_err_[0])) {
	    // New error is smaller than the old - we replace the objects
	    alg_obj_3d_int.getImplicit(impl, bc);
	    shared_ptr<SplineSurface> par_sf1 
		= IntersectionUtils::insertSfInImplObj(*sf1, impl, bc);
	    shared_ptr<SplineSurface> par_sf2
		= IntersectionUtils::insertSfInImplObj(*sf2, impl, bc);
	    approx_implicit_[0][0]
		= shared_ptr<Param2FunctionInt>
		(new Spline2FunctionInt(par_sf1));
	    approx_implicit_[0][1]
		= shared_ptr<Param2FunctionInt>
		(new Spline2FunctionInt(par_sf2));
	    approx_implicit_err_[0] = approx_implicit_err;
	    setGradValues(1);
	    prev_implicit_[0] = false;
	}
    }

    // Implicitize second
    if (implicitize_second) {
	spsf_int2->implicitize(dummy_tol);
	bool success;
	success = spsf_int2->getImplicit(dummy_tol,
					 approx_implicit_err,
					 alg_obj_3d_int);
	ASSERT(success);
	if (approx_implicit_err < fabs(approx_implicit_err_[1])) {
	    // New error is smaller than the old - we replace the objects
	    alg_obj_3d_int.getImplicit(impl, bc);
	    shared_ptr<SplineSurface> par_sf1 
		= IntersectionUtils::insertSfInImplObj(*sf1, impl, bc);
	    shared_ptr<SplineSurface> par_sf2
		= IntersectionUtils::insertSfInImplObj(*sf2, impl, bc);
	    approx_implicit_[1][0]
		= shared_ptr<Param2FunctionInt>
		(new Spline2FunctionInt(par_sf1));
	    approx_implicit_[1][1]
		= shared_ptr<Param2FunctionInt>
		(new Spline2FunctionInt(par_sf2));
	    approx_implicit_err_[1] = approx_implicit_err;
	    setGradValues(0);
	    prev_implicit_[1] = false;
	}
    }

    return;

}


//===========================================================================
void SfSfIntersector::setApproxImplicitFromPrev()
//===========================================================================
{
    SfSelfIntersector* prev_sf_self_int
	= dynamic_cast<SfSelfIntersector*>(prev_intersector_);
    if (prev_sf_self_int != 0) {
	// If the previous intersector is a self-intersector, we do
	// not set approximate implicit data
	return;
    }
    SfSfIntersector* prev_sf_sf_int
	= dynamic_cast<SfSfIntersector*>(prev_intersector_);
    ASSERT(prev_sf_sf_int != 0);

    // Get previous objects
    vector<vector<shared_ptr<Param2FunctionInt> > > approx_implicit;
    vector<double> approx_implicit_err;
    vector<double> approx_implicit_gradsize;
    vector<double> approx_implicit_gradvar;
    prev_sf_sf_int->getApproxImplicit(approx_implicit,
				      approx_implicit_err,
				      approx_implicit_gradsize,
				      approx_implicit_gradvar);

    // @@sbr Currently we are solely dealing with rectangular
    // domains.

    shared_ptr<const ParamSurface> prev_appr_sf1;
    shared_ptr<const ParamSurface> prev_appr_sf2;
    RectDomain prev_impl_domain;	
    shared_ptr<ParamSurfaceInt> sf_int;
    shared_ptr<ParamSurface> sf;
    RectDomain curr_impl_domain;
    shared_ptr<SplineSurface> spline_sf1;
    shared_ptr<SplineSurface> spline_sf2;
    // First object
    if (approx_implicit[0][0].get() != 0
	&& approx_implicit[1][0].get() != 0) {

	// Get prev_impl_domain
	prev_appr_sf1 = approx_implicit[0][0]->getSurface();
	prev_appr_sf2 = approx_implicit[1][0]->getSurface();
	prev_impl_domain = prev_appr_sf1->containingDomain();

	// Get curr_impl_domain
	sf_int = dynamic_pointer_cast<ParamSurfaceInt, ParamGeomInt>
	    (obj_int_[0]);
	sf = sf_int->getParamSurface();
	curr_impl_domain = sf->containingDomain();

	// If the domain lies inside previous domain we extract
	// the subpart.
	bool inside_previous
	    = ((curr_impl_domain.umin() > prev_impl_domain.umin()) ||
	       (curr_impl_domain.vmin() > prev_impl_domain.vmin()) ||
	       (curr_impl_domain.umax() < prev_impl_domain.umax()) ||
	       (curr_impl_domain.vmax() < prev_impl_domain.vmax()));
	if (inside_previous) {
	    vector<shared_ptr<ParamSurface> > sub_sf1 =
		prev_appr_sf1->subSurfaces(curr_impl_domain.umin(), 
					   curr_impl_domain.vmin(),
					   curr_impl_domain.umax(),
					   curr_impl_domain.vmax());
	    vector<shared_ptr<ParamSurface> > sub_sf2 =
		prev_appr_sf2->subSurfaces(curr_impl_domain.umin(), 
					   curr_impl_domain.vmin(),
					   curr_impl_domain.umax(),
					   curr_impl_domain.vmax());

	    // Not expecting the surfaces to be trimmed.
	    ASSERT(sub_sf1.size() == 1 && sub_sf2.size() == 1 );
	    spline_sf1 = dynamic_pointer_cast<SplineSurface, ParamSurface>
		(sub_sf1[0]);
	    spline_sf2 = dynamic_pointer_cast<SplineSurface, ParamSurface>
		(sub_sf2[0]);
	    approx_implicit_[0][0] = shared_ptr<Param2FunctionInt>
		(new Spline2FunctionInt(spline_sf1));
	    approx_implicit_[1][0] = shared_ptr<Param2FunctionInt>
		(new Spline2FunctionInt(spline_sf2));

	    // Gradient
	    setGradValues(0);




	} else {
	    approx_implicit_[0][0] = approx_implicit[0][0];
	    approx_implicit_[1][0] = approx_implicit[1][0];
	    approx_implicit_gradsize_[1] = approx_implicit_gradsize[1];
	    approx_implicit_gradvar_[1] = approx_implicit_gradvar[1];
	}
	prev_implicit_[0] = true;
    }
    // Second object
    if (approx_implicit[0][1].get() != 0
	&& approx_implicit[1][1].get() != 0) {

	// Get prev_impl_domain
	prev_appr_sf1 = approx_implicit[0][1]->getSurface();
	prev_appr_sf2 = approx_implicit[1][1]->getSurface();
	prev_impl_domain = prev_appr_sf1->containingDomain();

	// Get curr_impl_domain
	sf_int = dynamic_pointer_cast<ParamSurfaceInt, ParamGeomInt>
	    (obj_int_[1]);
	sf = sf_int->getParamSurface();
	curr_impl_domain = sf->containingDomain();

	// If the domain lies inside previous domain we extract
	// the subpart.
	bool inside_previous
	    = ((curr_impl_domain.umin() > prev_impl_domain.umin()) ||
	       (curr_impl_domain.vmin() > prev_impl_domain.vmin()) ||
	       (curr_impl_domain.umax() < prev_impl_domain.umax()) ||
	       (curr_impl_domain.vmax() < prev_impl_domain.vmax()));
	if (inside_previous) {
	    vector<shared_ptr<ParamSurface> > sub_sf1 =
		prev_appr_sf1->subSurfaces(curr_impl_domain.umin(), 
					   curr_impl_domain.vmin(),
					   curr_impl_domain.umax(),
					   curr_impl_domain.vmax());
	    vector<shared_ptr<ParamSurface> > sub_sf2 =
		prev_appr_sf2->subSurfaces(curr_impl_domain.umin(), 
					   curr_impl_domain.vmin(),
					   curr_impl_domain.umax(),
					   curr_impl_domain.vmax());

	    // Not expecting the surfaces to be trimmed.
	    ASSERT(sub_sf1.size() == 1 && sub_sf2.size() == 1 );
	    spline_sf1 = dynamic_pointer_cast<SplineSurface, ParamSurface>
		(sub_sf1[0]);
	    spline_sf2 = dynamic_pointer_cast<SplineSurface, ParamSurface>
		(sub_sf2[0]);
	    approx_implicit_[0][1] = shared_ptr<Param2FunctionInt>
		(new Spline2FunctionInt(spline_sf1));
	    approx_implicit_[1][1] = shared_ptr<Param2FunctionInt>
		(new Spline2FunctionInt(spline_sf2));

	    // Gradient
	    setGradValues(1);

	} else {
	    approx_implicit_[0][1] = approx_implicit[0][1];
	    approx_implicit_[1][1] = approx_implicit[1][1];
	    approx_implicit_gradsize_[0] = approx_implicit_gradsize[0];
	    approx_implicit_gradvar_[0] = approx_implicit_gradvar[0];
	}
	prev_implicit_[1] = true;
    }

    approx_implicit_err_ = approx_implicit_err;

    return;
}


//===========================================================================
void SfSfIntersector::setGradValues(int i)
//===========================================================================
{
    int j = (i == 0) ? 1 : 0;

    shared_ptr<Spline2FunctionInt> impl
	= dynamic_pointer_cast<Spline2FunctionInt, Param2FunctionInt>
	(approx_implicit_[j][i]);

    shared_ptr<SplineSurface> grad = impl->createGradSurface();
    typedef vector<double>::iterator iter;
    double maxlen = 0.0;
    double minlen = 1.0e100;
    double avlen = 0.0;
    for (iter it = grad->coefs_begin();
	 it != grad->coefs_end();
	 it += 2) {
	double len = *it * *it + *(it+1) * *(it+1);
	avlen += len;
	if (len > maxlen)
	    maxlen = len;
	if (len < minlen)
	    minlen = len;
    }
    int ncoefs = grad->numCoefs_u() * grad->numCoefs_v();
    approx_implicit_gradsize_[i] = sqrt(avlen / ncoefs);
    approx_implicit_gradvar_[i] = sqrt(maxlen / minlen);

    return;
}
//===========================================================================


} // namespace Go
