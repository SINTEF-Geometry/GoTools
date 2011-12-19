//===========================================================================
//                                                                           
// File: InterpolatedIntersectionCurve.C                                     
//                                                                           
// Created: Wed Apr 20 17:08:09 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision:
// $Id: InterpolatedIntersectionCurve.C,v 1.26 2007-11-01 14:31:44 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/IntersectionCurve.h"
#include "sislP.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/geometry/closestPtSurfSurfPlane.h"
#include "GoTools/geometry/HermiteInterpolator.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/IntersectionLink.h"
#include "GoTools/geometry/LineCloud.h" // debug
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <stdexcept>


using namespace Go;
using std::list;
using std::map;
using std::logic_error;
using std::vector;
using std::pair;
using std::copy;
using std::ofstream;

namespace {
    const double PI = 3.1415926535897932384;
    double angle_with_segment(list<shared_ptr<IntersectionPoint> >::const_iterator p1, list<shared_ptr<IntersectionPoint> >::const_iterator p2, const Go::Point& dir);
    void debug_dump_plane(const Point& basept, const Point& normal, double scale, char* name);

    // Determine whether the two points 'p1' and 'p2', as well as the two vectors
    // 't1' and 't2' are close enough together to be considered equal.  The tolerances
    // are given in the two last arguments; 'point_tol' and 'angle_tol'.
    static bool within_tolerances(const Go::Point& p1,
				  const Go::Point& t1,
				  const Go::Point& p2,
				  const Go::Point& t2,
				  double point_tol,
				  double angle_tol);

}; // end anonymous namespace


namespace Go {

//===========================================================================
shared_ptr<ParamCurve> InterpolatedIntersectionCurve::getCurve() const
//===========================================================================
{
    if (geom_cached_) {
	return cached_geom_curve_;
    }
    HermiteInterpolator interpolator;
    const int dim = 3; //@ can we always assume this?
    const int num_pts = (int)ipoints_.size() * 2; // constraints are positions + tangents
    vector<Point> data_in(num_pts);
    vector<double> coefs(dim * num_pts); // spline coefficients

    // preparing in-data
    list<shared_ptr<IntersectionPoint> >::const_iterator cur_point = ipoints_.begin();
    list<shared_ptr<IntersectionPoint> >::const_iterator end_point = ipoints_.end();
    ASSERT((*cur_point)->getPoint().dimension() == dim);
    vector<Point>::iterator data_iter = data_in.begin();
    // For the 3d-scenario we are interested in the space curve.
    while (cur_point != end_point) {
	*data_iter++ = (*cur_point)->getPoint();
	*data_iter++ = tangent_of(cur_point);
	++cur_point;
    }

    // determining parametrization
    // @@ for now, we use a heuristic approach by approximating arc-length by distance
    // between control point and using this as a parametrization.  Is there a better way??
    vector<double> param(ipoints_.size());
    param[0] = 0;
    cur_point = ipoints_.begin();
    const double minimum = (*cur_point)->getTolerance()->getRelParRes();
    for (int i = 1; i < int(param.size()); ++i) {
	Point p1 = (*cur_point++)->getPoint();
	Point p2 = (*cur_point)->getPoint();
	double dist = p2.dist(p1);
	dist = (dist > minimum) ? dist : minimum;
	param[i] = param[i-1] + dist;
    }
    // normalizing parametrization to get it between 0 and 1. 
    // NB: the use of RelParRes above is somewhat heuristic here, as renormalization will change
    // the values used.  However, for code simplicity we leave it like this for now, and will
    // aim to change it if it creates trouble in practice
    double startpar, endpar;
    getParamSpan(startpar, endpar);
    //double length_inv = (endpar - startpar) / param.back();
    //transform(param.begin(), param.end(), param.begin(), _1 * length_inv + startpar);
    
    interpolator.interpolate(data_in, param, coefs);

    BsplineBasis b = interpolator.basis();
    
    cached_geom_curve_ = 
	shared_ptr<ParamCurve>(new SplineCurve(num_pts,b.order(),b.begin(),&coefs[0],dim,false));
    cached_geom_curve_->geometryCurve()->setParameterInterval(startpar, endpar);
    geom_cached_ = true;
    
    return cached_geom_curve_;
}

//===========================================================================
shared_ptr<ParamCurve> InterpolatedIntersectionCurve::getParamCurve(int obj_nmb) const
//===========================================================================
{
    if (obj_nmb == 1 && par1_cached_) {
	return cached_param_curve_1_;
    } else if (obj_nmb == 2 && par2_cached_) {
	return cached_param_curve_2_;
    }

    if (obj_nmb != 1 && obj_nmb != 2) {
	throw logic_error("Argument to getParamCurve() should be either 1 or 2.");
    }

    HermiteInterpolator interpolator;
    const int dim = 2;
    const int num_pts = (int)ipoints_.size() * 2; // constraints are positions + tangents
    vector<Point> data_in(num_pts);
    vector<double> coefs(dim * num_pts); // spline coefficients

    // preparing in-data
    list<shared_ptr<IntersectionPoint> >::const_iterator cur_point = ipoints_.begin();
    list<shared_ptr<IntersectionPoint> >::const_iterator end_point = ipoints_.end();
    vector<Point>::iterator data_iter = data_in.begin();

    // extra check, to make sure we are dealing with a 2-parametric object
    ASSERT((obj_nmb == 1 && (*cur_point)->numParams1() == 2) ||
	   (obj_nmb == 2 && (*cur_point)->numParams2() == 2));

    while (cur_point != end_point) {
	if (obj_nmb == 1) {
	    Point par_pt((*cur_point)->getPar1(), (*cur_point)->getPar1() + 2);
	    *data_iter++ = par_pt;
	    Point par_dir = param_tangent_of(cur_point, false);
	    *data_iter++ = par_dir;
	} else {
	    Point par_pt((*cur_point)->getPar2(), (*cur_point)->getPar2() + 2);
	    *data_iter++ = par_pt;
	    Point par_dir = param_tangent_of(cur_point, true);
	    *data_iter++ = par_dir;
	} 
	++cur_point;
    }
    // determining parametrization
    // @@ for now, we use a heuristic approach by approximating arc-length by distance
    // between control point and using this as a parametrization.  Is there a better way??
    vector<double> param(ipoints_.size());
    param[0] = 0;
    cur_point = ipoints_.begin();
    for (int i = 1; i < int(param.size()); ++i) {
	Point p1(2), p2(2);
	if (obj_nmb == 1) {
	    p1.setValue((*cur_point++)->getPar1());
	    p2.setValue((*cur_point)->getPar1());
	} else {
	    p1.setValue((*cur_point++)->getPar2());
	    p2.setValue((*cur_point)->getPar2()); 
	}
	double dist = p2.dist(p1);
	// assuring a minimum parameter span
	const double minimum = (*cur_point)->getTolerance()->getRelParRes();
	dist = (dist > minimum) ? dist : minimum;
	param[i] = param[i-1] + dist;
    }
    // normalizing parametrization to get it between 0 and 1. 
    // NB: the use of RelParRes above is somewhat heuristic here, as renormalization will change
    // the values used.  However, for code simplicity we leave it like this for now, and will
    // aim to change it if it creates trouble in practice
    double startpar, endpar;
    getParamSpan(startpar, endpar);
    //double length_inv = (endpar - startpar) / param.back();
    //transform(param.begin(), param.end(), param.begin(), _1 * length_inv + startpar);

    interpolator.interpolate(data_in, param, coefs);

    BsplineBasis b = interpolator.basis();

    if (obj_nmb == 1) {
	par1_cached_ = true;
	cached_param_curve_1_ = 
	    shared_ptr<ParamCurve>(new SplineCurve(num_pts,b.order(),b.begin(),&coefs[0],dim,false));
	cached_param_curve_1_->geometryCurve()->setParameterInterval(startpar, endpar);	
	return cached_param_curve_1_;
    } 
    // obj_nmb == 2
    par2_cached_ = true;
    cached_param_curve_2_ = 
	shared_ptr<ParamCurve>(new SplineCurve(num_pts,b.order(),b.begin(),&coefs[0],dim,false));
    cached_param_curve_2_->geometryCurve()->setParameterInterval(startpar, endpar);	
    return cached_param_curve_2_;
}

//===========================================================================
void InterpolatedIntersectionCurve::refine(const double& pos_tol, const double& angle_tol)
//===========================================================================
{
    if (pos_tol >= certified_pos_tol_ && angle_tol >= certified_angle_tol_) {
	return;
    }
    // determining intervals
    vector<list<shared_ptr<IntersectionPoint> >::iterator> start_points;
    list<shared_ptr<IntersectionPoint> >::iterator iter  = ipoints_.begin();
    list<shared_ptr<IntersectionPoint> >::iterator end_iter = ipoints_.end();
    --end_iter;
    while(iter != end_iter) {
	start_points.push_back(iter++);
    }

    // recursively refine each defined interval on the curve
    for (int i = 0; i < int(start_points.size()); ++i) {
	refine_interval_recursive(start_points[i], pos_tol, angle_tol);
    }
    certified_pos_tol_ = (certified_pos_tol_ < pos_tol) ? certified_pos_tol_ : pos_tol;
    certified_angle_tol_ = (certified_angle_tol_ < pos_tol) ? certified_angle_tol_ : pos_tol;
    geom_cached_ = par1_cached_ = par2_cached_ = false;
}

//===========================================================================
void InterpolatedIntersectionCurve::evaluateAt(double pval, Point& pos, Point& tan) 
//===========================================================================
{
    // refine the interpolated curve within the geometric tolerance
    shared_ptr<GeoTol> gtol = ipoints_.front()->getTolerance();
    const double pos_tol = gtol->getEpsge();
    const double angle_tol = gtol->getAngleTol();
    refine(pos_tol, angle_tol);
    
    // generate arguments for iterative procedure
    temp_.resize(2);
    getCurve()->point(temp_, pval, 1);
    Point midpoint_param_pos_1;
    Point midpoint_param_pos_2;
    getParamCurve(1)->point(midpoint_param_pos_1, pval);
    getParamCurve(2)->point(midpoint_param_pos_2, pval);

    const ParamSurfaceInt* psurf1 = 
	dynamic_cast<const ParamSurfaceInt*>(ipoints_.front()->getObj1());
    const ParamSurfaceInt* psurf2 = 
	dynamic_cast<const ParamSurfaceInt*>(ipoints_.front()->getObj2());
    DEBUG_ERROR_IF(!psurf1 || !psurf2, "evaluateAt() currently only works "
	     "for objects of type ParamSurfaceInt.");

    Point surface_1_param, surface_2_param;
    int jstat;

    // iterate down to a point on the real intersection curve
    bool tangent_found = eval_surf_point(temp_[0], // seed position
					 temp_[1], // seed tangent
					 psurf1,
					 midpoint_param_pos_1,
					 psurf2,
					 midpoint_param_pos_2,
					 pos,
					 tan,
					 surface_1_param,
					 surface_2_param, 
					 jstat);
    if (jstat > 2) {
	MESSAGE("WARNING!  Unable to recursively refine sub-interval.  "
		"Abandoning refinement.");
	return;
    }

    if (!tangent_found) {
	tan = temp_[1]; // @@ is this wise??
    }
}

//===========================================================================
void InterpolatedIntersectionCurve::
refine_interval_recursive(list<shared_ptr<IntersectionPoint> >::iterator startpt, const double& pos_tol, const double& angle_tol)
//===========================================================================
{
    // we suppose that startpt is the first IntersectionPoint of a defined interval,
    // so that we can be sure that (++startpt) exists.
    list<shared_ptr<IntersectionPoint> >::iterator endpt = startpt;
    ++endpt;

    // If start and endpoint coincide, refining does not make any sense, and 
    // hermitian interpolation is dangerous!
    Point diff = (*startpt)->getPoint() - (*endpt)->getPoint();
    if (diff.length2() < pos_tol * pos_tol) {
	MESSAGE("WARNING!  Start and end point on curve to be refined coincide! "
		"Abandoning refinement!");
	return;
    }

    // hermite interpolation of geometric point
    Point mid_pt, mid_tg;
    hermite_interpol(startpt, endpt, mid_pt, mid_tg, SPACECURVE);

    // hermite interpolation of parameters on first surface
    Point par_1_mid_pos, par_1_mid_tan;
    hermite_interpol(startpt, endpt, par_1_mid_pos, par_1_mid_tan, PARAMCURVE_1);

    // hermite interpolation of parameters on second surface
    Point par_2_mid_pos, par_2_mid_tan;
    hermite_interpol(startpt, endpt, par_2_mid_pos, par_2_mid_tan, PARAMCURVE_2);

    // get pointers to underlying surfaces, to use for searching for the closest
    // point on the real intersection
    const ParamSurfaceInt* psurf1 = dynamic_cast<const ParamSurfaceInt*>((*startpt)->getObj1());
    const ParamSurfaceInt* psurf2 = dynamic_cast<const ParamSurfaceInt*>((*startpt)->getObj2());
    DEBUG_ERROR_IF(!psurf1 || !psurf2, "Refine_interval_recursive currently only works "
	     "for objects of type ParamSurfaceInt.");

    // truncating parameters in case of outside domain
    for (int i = 0; i < 2; ++i) {
	if (par_1_mid_pos[i] < psurf1->startParam(i)) {
	    par_1_mid_pos[i] = psurf1->startParam(i);
	} else if (par_1_mid_pos[i] > psurf1->endParam(i)) {
	    par_1_mid_pos[i] = psurf1->endParam(i);
	}
	if (par_2_mid_pos[i] < psurf2->startParam(i)) {
	    par_2_mid_pos[i] = psurf2->startParam(i);
	} else if (par_2_mid_pos[i] > psurf2->endParam(i)) {
	    par_2_mid_pos[i] = psurf2->endParam(i);
	}
    }

    // detecting closest point on real intersection
    int jstat;
    Point sf_point(3), sf_tangent(3), sf_1_prm(2), sf_2_prm(2);

    bool tangent_found = eval_surf_point(mid_pt, mid_tg,
					 psurf1, par_1_mid_pos, 
					 psurf2, par_2_mid_pos, 
					 sf_point, sf_tangent, 
					 sf_1_prm, sf_2_prm,
					 jstat);
//   debug
//      Point p1 = (*startpt)->getPoint();
//      Point p2 = (*endpt)->getPoint();
//      Point p3 = sf_point;
//      ofstream os("debugpt.g2");
//      os << "400 1 0 4 255 255 0 255" << endl;
//      os << 3; // number of points
//      os << p1[0] << " " << p1[1] << " " << p1[2] << endl;
//      os << p2[0] << " " << p2[1] << " " << p2[2] << endl;
//      os << p3[0] << " " << p3[1] << " " << p3[2] << endl;
//      os.close();
//      ofstream ostan("debugtan.g2");
//      double tan_len = diff.length();
//      //Point tan1 = (*startpt)->getTangent();
//      Point tan1 = tangent_of(startpt);
//      tan1.normalize();
//      tan1 *= tan_len / 3;
//      //Point tan2 = (*endpt)->getTangent();
//      Point tan2 = tangent_of(endpt);
//      tan2.normalize();
//      tan2 *= tan_len / 3;
//      Point tan3 = sf_tangent;
//      if (tangent_found) 
// 	 tan3.normalize();
//      tan3 *= tan_len / 3;
//      vector<double> temp;
//      temp.push_back(p1[0]);
//      temp.push_back(p1[1]);
//      temp.push_back(p1[2]);
//      temp.push_back(p1[0] + tan1[0]);
//      temp.push_back(p1[1] + tan1[1]);
//      temp.push_back(p1[2] + tan1[2]);
//      temp.push_back(p2[0]);
//      temp.push_back(p2[1]);
//      temp.push_back(p2[2]);
//      temp.push_back(p2[0] + tan2[0]);
//      temp.push_back(p2[1] + tan2[1]);
//      temp.push_back(p2[2] + tan2[2]);
//      temp.push_back(p3[0]);
//      temp.push_back(p3[1]);
//      temp.push_back(p3[2]);
//      temp.push_back(p3[0] + tan3[0]);
//      temp.push_back(p3[1] + tan3[1]);
//      temp.push_back(p3[2] + tan3[2]);
//      LineCloud lcloud(temp.begin(), 3);
//      lcloud.writeStandardHeader(ostan);
//      lcloud.write(ostan);
//      ostan.close();

    
    if (jstat > 2) {
	MESSAGE("WARNING!  Unable to recursively refine sub-interval.  "
		"Abandoning refinement.");
	return;
    }

    if (!tangent_found) {
	sf_tangent = mid_tg; // @@ is this wise??
    }

    //Satisfying tolerance?
    if (within_tolerances(sf_point, sf_tangent, mid_pt, mid_tg, pos_tol, angle_tol)) {
	return;
    } else {
	// the hermite interpolation is not exact enough. Inserting additional guidepoint
	// and checking end intervals
	shared_ptr<IntersectionPoint> 
	    new_guidepoint(new IntersectionPoint((*startpt)->getObj1(), //obj1,
						 (*startpt)->getObj2(), //obj2, 
						 (*startpt)->getTolerance(), 
						 sf_1_prm.begin(),
						 sf_2_prm.begin()));

	list<shared_ptr<IntersectionPoint> >::iterator newpt = ipoints_.insert(endpt, new_guidepoint);

	// make sure the point is differentiated from the correct side
	choose_differentiation_side(newpt);
	    
	// making sure the curve has an internal (possibly unoriented) 
	// representation of this tangent
	establish_curve_spesific_tangent(newpt);

	if (!(*newpt)->tangentIsOriented()) {
	    // fixing orientation of tangent based on next point in sequence
	    Point this_tangent = tangent_of(newpt);
	    Point next_tangent = tangent_of(endpt);
	    bool flip = determine_flip(startpt, newpt, endpt);
	    flip_tangent(newpt, flip);
	}

	// recursively refining each of the created subintervals
	refine_interval_recursive(startpt, pos_tol, angle_tol); 
	refine_interval_recursive(newpt, pos_tol, angle_tol);
    }
    return;
}

//===========================================================================
bool InterpolatedIntersectionCurve::eval_surf_point(const Point& midpoint_pos,
						    const Point& midpoint_tan,
						    const ParamSurfaceInt* psurf1,
						    const Point& midpoint_param_pos_1,
						    const ParamSurfaceInt* psurf2,
						    const Point& midpoint_param_pos_2,
						    Point& surface_point,
						    Point& surface_tangent,
						    Point& surface_1_param,
						    Point& surface_2_param,
						    int& jstat) const 
//===========================================================================
{
    const double angle_tol_multiplier = 1.0e-4; // estimated from stabiltiy requirements,
                                                // see also IntersectionPoint::calculate_
                                                // stability_type_surface()
    const double ang_tol = ipoints_.front()->getTolerance()->getAngleTol();
    const double tan_length_tol = angle_tol_multiplier / ang_tol;
    const double num_tol = ipoints_.front()->getTolerance()->getNumericalTol();

    vector<Point> plane_definition(2);
    plane_definition[0] = midpoint_pos;
    plane_definition[1] = midpoint_tan;
    
    vector<Point> input_point_1(7, Point(3)); // only entry [0], [1], [2] and [6] will be used
    vector<Point> input_point_2(7, Point(3)); // only entry [1], [1], [2] and [6] will be used
    vector<Point> result_pt_1(7, Point(3));
    vector<Point> result_pt_2(7, Point(3));
    
    shared_ptr<const ParamSurface> s1 = psurf1->getParamSurface();
    shared_ptr<const ParamSurface> s2 = psurf2->getParamSurface();

    // fill in starting points information

    // in reality, the closestPtSurfSurfPlane routine demand that the second derivatives
    // should also be calculated, but an inspection of its inner workings reveal that it
    // never makes use of this information.  So in order to save a few cycles, we only 
    // evaluate the first derivatives and leave the other entries blank.  This can be
    // changed back by replacing the below two lines with the commented ones.
    s1->point(input_point_1, midpoint_param_pos_1[0], midpoint_param_pos_1[1], 1);
    s2->point(input_point_2, midpoint_param_pos_2[0], midpoint_param_pos_2[1], 1);
    //s1->point(input_point_1, midpoint_param_pos_1[0], midpoint_param_pos_1[1], 2);
    //s2->point(input_point_2, midpoint_param_pos_2[0], midpoint_param_pos_2[1], 2);

    // assume that the above two routines have not changed the size of the vectors
    ASSERT(input_point_1.size() == 7 && input_point_2.size() == 7); 

    // calculating normals
    s1->normal(input_point_1[6], midpoint_param_pos_1[0], midpoint_param_pos_1[1]);
    s2->normal(input_point_2[6], midpoint_param_pos_2[0], midpoint_param_pos_2[1]);
    
    // krull
    //debug_dump_plane(midpoint_pos, midpoint_tan, 10, "dump_plane.g2");

    AlgorithmChoice algo = FUNCTIONAL; //GEOMETRICAL;
    closestPtSurfSurfPlane(plane_definition, 
			   input_point_1, 
			   input_point_2,
			   midpoint_param_pos_1,
			   midpoint_param_pos_2,
			   psurf1->getParamSurface().get(),
			   psurf2->getParamSurface().get(),
			   sqrt(num_tol), // precision in minimization routine cannot
			              // be better than root of computational precision
			   result_pt_1,
			   result_pt_2,
			   surface_1_param,
			   surface_2_param,
			   jstat,
			   algo);

    surface_point = 0.5 * (result_pt_1[0] + result_pt_2[0]); // average from each object
    surface_tangent = result_pt_1[6].cross(result_pt_2[6]);
    //double eps = IntersectionPoint::tangent_tol_;
    if (surface_tangent.length2() > tan_length_tol * tan_length_tol) {
	surface_tangent.normalize();
    } else {
	// Surface tangent currently estimated to 0-vector.
	// Constructing a temporary IntersectionPoint to benefit from its robust tangent
	// calculation routines.
	
	// @making necessary GeoTol object - this is a pain here, since it will not
	// be used.  Any other solution?
	shared_ptr<GeoTol> temp_tol = ipoints_.front()->getTolerance();

	shared_ptr<IntersectionPoint> temp(new IntersectionPoint(psurf1, 
								 psurf2, 
								 temp_tol, // irrelevant here
								 surface_1_param.begin(),
								 surface_2_param.begin()));
	// if the following line fails, then we might not have any reasonable way
	// to calculate the tangent in this point
	try {
	    surface_tangent = temp->getTangent();
	} catch (...) {
	    // oh no! Not even this method could produce a usable tangent.
	    // I give up!  Return empty tangent, and report 'false'
	    surface_tangent = Point();
	    return false;
	}
    }
    return true;
}


//===========================================================================
void InterpolatedIntersectionCurve::resolve_tangents()
//===========================================================================
{
    // The IntersectionPoints may be in one of the following 3 cathegories:
    // (A) Points with a well-defined tangent AND orientation 
    // (B) Points with a well-defined tangent direction, but no orientation
    // (C) Points with no clearly defined tangent direction nor orientation
    //    (branch-points, higher-order-points)
    // We first aim to make sure that all points have a _direction_, obtainable
    // via 'tangent_of()'.  Then we make sure that they are all consistently
    // oriented.  

    // Prior to all calculations we must, however, make sure that differentiation
    // is carried out from the correct side.
    for (list<shared_ptr<IntersectionPoint> >::const_iterator it = ipoints_.begin(); it != ipoints_.end(); ++it) {
	choose_differentiation_side(it);
    }

    // Registering tangent directions of all points.  In case they are of
    // type (B) or (C), the direction will have to be estimated somehow.
    for (list<shared_ptr<IntersectionPoint> >::const_iterator it = ipoints_.begin(); it != ipoints_.end(); ++it) {
	establish_curve_spesific_tangent(it);
    }

    // Now, we want to make tangent orientations consistent.  This means
    // to attribute consistent orientations to points of type (B) and (C) based
    // on the orientation of points of type (A).
    invert_sequence_if_necessary(); 
    list<shared_ptr<IntersectionPoint> >::iterator ref_elem = choose_reference_direction();
    make_consistent_orientation(ref_elem);
}

//===========================================================================
void InterpolatedIntersectionCurve::
choose_differentiation_side(list<shared_ptr<IntersectionPoint> >::const_iterator pt) const
//===========================================================================
{
    static vector<bool> diff_from_left;
    int num_param = (*pt)->numParams1() + (*pt)->numParams2();
    diff_from_left.resize(num_param);
    list<shared_ptr<IntersectionPoint> >::const_iterator neigh_pt = pt;
    if (pt != ipoints_.begin()) {
	// adjusting differentiating side of this point according to relation with
	// previous point
	--neigh_pt;
    } else {
	// adjusting differentiating side of this point according to relation with
	// next point
	++neigh_pt;
    }

    for (int p = 0; p < num_param; ++p) {
	diff_from_left[p] = ((*pt)->getPar(p) > (*neigh_pt)->getPar(p));
    }
    (*pt)->setDifferentiateFromLeft(diff_from_left.begin());
}


//===========================================================================
void InterpolatedIntersectionCurve::invert_sequence_if_necessary()
//===========================================================================
{
    ASSERT(ipoints_.size() >= 2);

    // searching for a well-defined tangent
    for (list<shared_ptr<IntersectionPoint> >::iterator it = ipoints_.begin(); it != ipoints_.end(); ++it) {
	const shared_ptr<IntersectionPoint>& ip = *it;
	if (ip->tangentIsOriented()) {
	    const Point tangent = ip->getTangent();
	    list<shared_ptr<IntersectionPoint> >::iterator temp = it;
	    ++temp;
	    Point cord;
	    if (temp != ipoints_.end()) {
		cord = (*temp)->getPoint() - ip->getPoint();
	    } else {
		--temp;	--temp;
		cord = ip->getPoint() - (*temp)->getPoint();
	    }
	    // we have now found a well-defined tangent and will
	    // analyse it.  We suppose that all points with
	    // well-defined tangents are consistent in this (if not,
	    // it would be impossible to determine a 'correct'
	    // direction), so we consider ourselves finished with the
	    // task of determining the right direction of the
	    // sequence.
	    if (tangent.angle(cord) > 0.5 * PI) {
		// points seem to be listed in opposite order
		// (compared to found tangent). Inverse order and
		// return
		reverse(ipoints_.begin(), ipoints_.end());
	    } 
	    return;
	}
    }
}


//===========================================================================
list<shared_ptr<IntersectionPoint> >::iterator InterpolatedIntersectionCurve::choose_reference_direction()
//===========================================================================
{
    list<shared_ptr<IntersectionPoint> >::iterator ref_point = ipoints_.begin();
    while (ref_point != ipoints_.end() && !(*ref_point)->tangentIsOriented()) {
	++ref_point;
    }
    if (ref_point == ipoints_.end()) {
	// We did not succeed in finding any oriented tangent in the range.
	// Let us just pick the first one and set a reasonable orientation.
	ref_point = ipoints_.begin();
	list<shared_ptr<IntersectionPoint> >::const_iterator next_point = ref_point;
	++next_point;
	double ang = angle_with_segment(ref_point, next_point, tangent_of(ref_point));
	bool flip = (ang > PI * 0.5);
	flip_tangent(ref_point, flip);
    }
    return ref_point;
}


//===========================================================================
void InterpolatedIntersectionCurve::
make_consistent_orientation(list<shared_ptr<IntersectionPoint> >::iterator ref_point)
//===========================================================================
{
    list<shared_ptr<IntersectionPoint> >::iterator cur_point = ref_point;
    list<shared_ptr<IntersectionPoint> >::iterator old_point = ref_point;
    list<shared_ptr<IntersectionPoint> >::iterator last_point = ipoints_.end();
    --last_point;
    Point ref_tangent, cur_tangent;
    double angle_tol = (*ref_point)->getTolerance()->getAngleTol();

    // iterating forwards until last element
    while (cur_point != last_point) {
	old_point = cur_point++;
	if (!(*cur_point)->tangentIsOriented()) {
	    bool flip;
	    if (cur_point != last_point) {
		list<shared_ptr<IntersectionPoint> >::iterator next_point = cur_point;
		++next_point;
		flip =  determine_flip(old_point, cur_point, next_point);
	    } else {
		// we are at the endpoint.  Orient this point according to old point
		// @@sbr if geom space is not 3d we use the parametric tangent
		int dim = (*old_point)->getPoint().dimension();
		ref_tangent = (dim == 3) ?
		    tangent_of(old_point) : param_tangent_of(old_point, false);
		cur_tangent = (dim == 3) ?
		    tangent_of(cur_point) : param_tangent_of(cur_point, false);
		double ang = ref_tangent.angle(cur_tangent);
		if (fabs(fabs(ang) - PI * 0.5) < angle_tol) { // perpendicular tangents
		    ang = angle_with_segment(old_point, cur_point, cur_tangent);
		}
		flip = (ang > PI * 0.5);
	    }
	    flip_tangent(cur_point, flip);
	}
    }
    // iterating backwards until last element
    cur_point = ref_point;
    while (cur_point != ipoints_.begin()) {
	old_point = cur_point--;
	if (!(*cur_point)->tangentIsOriented()) {
	    bool flip;
	    if (cur_point != ipoints_.begin()) {
		list<shared_ptr<IntersectionPoint> >::iterator prec_point = cur_point;
		--prec_point;
		flip = determine_flip(prec_point, cur_point, old_point);
	    } else {
		// we are at the start point
		// orient this point according to old point
		int dim = (*old_point)->getPoint().dimension();
		ref_tangent = (dim == 3) ? 
		    tangent_of(old_point) : param_tangent_of(old_point, false);
		cur_tangent = (dim == 3) ? 
		    tangent_of(cur_point) : param_tangent_of(cur_point, false);
		double ang = ref_tangent.angle(cur_tangent);
		if (fabs(fabs(ang) - PI * 0.5) < angle_tol) { // perpendicular tangents
		    ang = angle_with_segment(cur_point, old_point, cur_tangent);
		}

		flip = (ang > PI * 0.5);
	    } 
	    flip_tangent(cur_point, flip);
	}
    }
}


//===========================================================================
Point InterpolatedIntersectionCurve::tangent_of(list<shared_ptr<IntersectionPoint> >::const_iterator pt) const
//===========================================================================
{
    map<shared_ptr<IntersectionPoint>,
	pair<Array<Point, 3>, bool> >::const_iterator it;
    it = tangents_.find(*pt);
    ASSERT(it != tangents_.end()); // this would be a 'breach of contract'  
                                   // by the constructor of this object.
    if (it->second.second) { // if tangent is "independent"
	return it->second.first[0];
    } 

    // if we got here, the point has no independent tangent.  Taking neighbours
    // into consideration
    
    vector<Point> estimates;
    Point tn;
    bool found = context_tangent_estimate(pt, GEOM, FORWARDS, tn);
    if (found) { 
	estimates.push_back(tn);
    }
    found = context_tangent_estimate(pt, GEOM, BACKWARDS, tn);
    if (found) {
	estimates.push_back(tn);
    }
    ASSERT(estimates.size() > 0); // at least one tangent should be guaranteed to be
                                  // found, since the curve is not degenerated
    if (estimates.size() == 1) {
	return estimates[0];
    }
    // two estimates were made; returning "average"
    Point average = estimates[0] + estimates[1];
    average.normalize(); // @@ this is not exactly what should be done, but it is a good
                         // approximation as long as the tangents do not vary to much.  
                         // If this solutions proves to be insufficient, try to compute 
                         // the average using quaternions, or some other technique.
    return average;
}

//===========================================================================
Point InterpolatedIntersectionCurve::param_tangent_of(list<shared_ptr<IntersectionPoint> >::const_iterator pt, bool second_obj) const 
//===========================================================================
{
    map<shared_ptr<IntersectionPoint>,
	pair<Array<Point, 3>, bool> >::const_iterator it;
    it = tangents_.find(*pt);
    ASSERT(it != tangents_.end()); // this would be a 'breach of contract'  
                                   // by the constructor of this object.
    if (it->second.second) { // if tangent is "independent"
	return second_obj ? it->second.first[2] : it->second.first[1];
    }
    
    // if we got here, the point has no independent tangent.  Taking neighbours
    // into consideration.
    vector<Point> estimates;
    TangentDomain dom  = second_obj ? PARAM2 : PARAM1;
    Point tn;
    bool found = context_tangent_estimate(pt, dom, FORWARDS, tn);
    if (found) {
	estimates.push_back(tn);
    }
    found = context_tangent_estimate(pt, dom, BACKWARDS, tn);
    if (found) {
	estimates.push_back(tn);
    }
    ASSERT(estimates.size() > 0); // at least one tangent should be guaranteed to be found,
                                  // since the curve is not degenerated
    if (estimates.size() == 1) {
	return estimates[0];
    }
    // two estimates were made; returning "average"
    Point average = estimates[0] + estimates[1];
    average.normalize(); // @@ this is not exactly what should be done, but it is a good
                         // approximation as long as the tangents do not vary to much.  
                         // If this solutions proves to be insufficient, try to compute 
                         // the average using quaternions, or some other technique.
    return average;
}

//===========================================================================
bool InterpolatedIntersectionCurve::determine_flip(list<shared_ptr<IntersectionPoint> >::const_iterator prec, 
						   list<shared_ptr<IntersectionPoint> >::const_iterator mid, 
						   list<shared_ptr<IntersectionPoint> >::const_iterator next)
//===========================================================================
{
    // determine whether the tangent of the 'mid' point should be flipped in order
    // to correspond better with its neighbour points
    const double eps = (*mid)->getTolerance()->getEpsge();
    Point mid_tangent = tangent_of(mid);
    Point dist = (*next)->getPoint() - (*prec)->getPoint();
    if (dist.length2() < eps * eps) {
	// the two neighbour points have practically the same coordinates.
	// We will rather work in the parameter plane of the first surface (arbitrary
	// choice).  @@ what if the first surface is not parametric??
	mid_tangent = param_tangent_of(mid, false);
	Point p1((*next)->getPar1(), (*next)->getPar1() + (*next)->numParams1());
	Point p2((*prec)->getPar1(), (*prec)->getPar1() + (*prec)->numParams1());
	dist = p1 - p2;
    }

    if (mid_tangent.angle(dist) < 0.5 * PI) {
	// no need for flipping
	return false;
    }

    // tangent should be flipped
    return true;
}


//===========================================================================
void InterpolatedIntersectionCurve::flip_tangent(list<shared_ptr<IntersectionPoint> >::iterator pt, bool flip)
//===========================================================================
{
    ASSERT(!(*pt)->tangentIsOriented()); 
    
    if (!flip) {
	return;
    }
    map<shared_ptr<IntersectionPoint>,
	pair<Array<Point, 3>, bool> >::iterator it;
    it = tangents_.find(*pt);
    ASSERT(it != tangents_.end()); // should be guaranteed
    if (it->second.second) { // tangent is defined
	for (int i = 0; i < 3; ++i) {
	    it->second.first[i] *= -1;
	}
    }
}


//===========================================================================
void InterpolatedIntersectionCurve::
establish_curve_spesific_tangent(list<shared_ptr<IntersectionPoint> >::const_iterator pt)
//===========================================================================
{
    // only direction, not orientation is taken care of in this
    // function

    tangents_[*pt] = pair<Array<Point, 3>, bool>(Array<Point,3>(), true);
    pair<Array<Point, 3>, bool>& cur_entry = tangents_[*pt];
    Array<Point, 3>& result = cur_entry.first;
    bool& tangent_found = cur_entry.second;

    Point& tangent_3d = result[0];
    Point& par_tang_1 = result[1];
    Point& par_tang_2 = result[2];

    // Let us have a look at an eventual tangent direction provided by the point
    if ((*pt)->hasUniqueTangentDirection()) {
	// an unique tangent direction exists, even though the point does not necessarily know 
	// which way the tangent is supposed to be oriented
	tangent_3d = (*pt)->getTangent();
	par_tang_1 = (*pt)->getPar1Dir();
	par_tang_2 = (*pt)->getPar2Dir();
	// ASSERT(par_tang_1.dimension() == 2 && par_tang_2.dimension() == 2);
	return;
    } 
    
    // If we got here, the point could not help us choose a tangent.  We must use
    // more general geometry considerations to estimate a direction.
    vector<Point> estimates;
    bool found = context_tangent_estimate(pt, GEOM, FORWARDS, tangent_3d);
    if (found) {
	estimates.push_back(tangent_3d);
    }
    found = context_tangent_estimate(pt, GEOM, BACKWARDS, tangent_3d);
    if (found) {
	estimates.push_back(tangent_3d);
    }
    ASSERT(int(estimates.size()) > 0);
    tangent_3d = estimates[0];
    if (estimates.size() == 2) {
	tangent_3d += estimates[1];
    }
    tangent_3d.normalize();
    //tangent_3d = geometric_tangent_estimate(pt); 

    // The 3D-tangent should be nozero by the time we get to this line.

    SingularityType stype = (*pt)->getSingularityType();
    if (stype == BRANCH_POINT) {
	// point has two well-defined tangent directions.  We must choose one.
	Point dir_1 = (*pt)->getTangent(false); // first branch
	Point dir_2 = (*pt)->getTangent(true);  // second branch
	double prod1 = (tangent_3d * dir_1) / dir_1.length();
	double prod2 = (tangent_3d * dir_2) / dir_2.length();
	if (fabs(prod1) > fabs(prod2)) {
	    // choose branch 1
	    tangent_3d = dir_1;
	    par_tang_1 = (*pt)->getPar1Dir(false);
	    par_tang_2 = (*pt)->getPar2Dir(false);
	} else {
	    // choose branch 2
	    tangent_3d = dir_2;
	    par_tang_1 = (*pt)->getPar1Dir(true);
	    par_tang_2 = (*pt)->getPar2Dir(true);
	}
	ASSERT(par_tang_1.dimension() == 2 && par_tang_2.dimension() == 2);
    } else {
	ASSERT(stype == ISOLATED_POINT || stype == HIGHER_ORDER_POINT);
	// point has no clue about its tangent.  
	tangent_found = false; // this point would have to have a "flexible" tangent,
	                       // which is always adjusted to conform to the closest
	                       // neighbors
	//(*pt)->projectToParamPlanes(tangent_3d, par_tang_1, par_tang_2);	
    }
}

//===========================================================================
bool InterpolatedIntersectionCurve::
context_tangent_estimate(list<shared_ptr<IntersectionPoint> >::const_iterator pt, 
			 TangentDomain tdom, 
			 EstimateDirection dir, 
			 Point& result) const
//===========================================================================
{
    Point (IntersectionPoint::*position)() const;
    Point (IntersectionPoint::*tan_fct)(bool) const;

    double tol = 1.0e20;
    int num_par_1 = (*pt)->numParams1();
    int num_par_2 = (*pt)->numParams2();
    switch (tdom) {
    case GEOM:	 
	tol = (*pt)->getTolerance()->getEpsge();
	position = &IntersectionPoint::getPoint; 
	tan_fct = &IntersectionPoint::getTangent;
	break;
    case PARAM1: 
	tol = (*pt)->parameterTolerance(0);
	if (num_par_1 == 2 && (*pt)->parameterTolerance(1) < tol) {
	    tol = (*pt)->parameterTolerance(1);
	}
	position = &IntersectionPoint::getPar1Point; 
	tan_fct = &IntersectionPoint::getPar1Dir;
	break;
    case PARAM2: 
	tol = (*pt)->parameterTolerance(num_par_1);
	if (num_par_2 == 2
	    && (*pt)->parameterTolerance(num_par_1 + 1) < tol) {
	    tol = (*pt)->parameterTolerance(num_par_1 + 1);
	}
	position = &IntersectionPoint::getPar2Point; 
	tan_fct = &IntersectionPoint::getPar2Dir;
	break;
    default: 	 
	ASSERT(false); // should never get here
    }

    bool found = false;
    list<shared_ptr<IntersectionPoint> >::const_iterator other_pt = pt;
    Point pt_pos = ((*pt).get()->*position)(); 
    Point other_pos;
    Point d;
    double d_len_2 = 0.0;
    if (dir == FORWARDS) { // searching forwards
	other_pt = pt;
	for (++other_pt; other_pt != ipoints_.end(); ++other_pt) {
	    other_pos = ((*other_pt).get()->*position)();
	    d = other_pos - pt_pos;
	    d_len_2 = d.length2();
	    if (d_len_2 > tol * tol) {
		found = true;
		break;
	    }
	}
    } else { // searching backwards
	while (other_pt != ipoints_.begin()) {
	    --other_pt;
	    other_pos = ((*other_pt).get()->*position)();
	    d = other_pos - pt_pos;
	    d_len_2 = d.length2();
	    if (d_len_2 > tol * tol) {
		found = true;
		break;
	    }
	}
    }
    if (found) {
	const Point other_tan = (*other_pt)->hasUniqueTangentDirection() ? 
	    ((*other_pt).get()->*tan_fct)(false) :  other_pos - pt_pos;
	ASSERT(d_len_2 > tol * tol); // this should have been taken
				     // care of above
	result = (2 * (other_tan * d) / d_len_2) * d - other_tan;	
    }
    return found;
}

// //===========================================================================
// Point InterpolatedIntersectionCurve::geometric_tangent_estimate(list<shared_ptr<IntersectionPoint> >::const_iterator pt) const
// //===========================================================================
// {
//     // finding other point
//     const double tol = (*pt)->getTolerance()->getEpsge();
//     list<shared_ptr<IntersectionPoint> >::const_iterator other_pt;
//     Point pt_pos = (*pt)->getPoint();
//     Point other_pos;
//     Point d;
//     double d_len_2;

//     // searching forwards
//     ASSERT(pt != ipoints_.end());
//     other_pt = pt;
//     for (++other_pt; other_pt != ipoints_.end(); ++other_pt) {
// 	other_pos = (*other_pt)->getPoint();
// 	d = other_pos - pt_pos;
// 	d_len_2 = d.length2();
// 	if (d_len_2 > tol * tol) {
// 	    break;
// 	}
//     }

//     if (other_pt == ipoints_.end()) {
// 	// could not find forward point
// 	bool found = false;
// 	// searching backwards
// 	if (pt != ipoints_.begin()) {
// 	    other_pt = pt;
// 	    do {
// 		--other_pt;
// 		other_pos = (*other_pt)->getPoint();
// 		d = other_pos - pt_pos;
// 		d_len_2 = d.length2();
// 		if (d_len_2 > tol * tol) {
// 		    found = true;
// 		    break;
// 		}
// 	    } while (other_pt != ipoints_.begin());
// 	}

// 	ASSERT(found); // if not found, it means that all points were coincident.  This
// 	               // is a 'breach of contract' from what the constructor expects.
//     }

//     Point other_tan;
    
//     if ((*other_pt)->hasUniqueTangentDirection()) {
// 	other_tan = (*other_pt)->getTangent();
//     } else {
// 	other_tan = other_pos - pt_pos;
//     }

//     ASSERT(d_len_2 > tol * tol); // this should have been taken care of above
//     // computing tangent as if it the curve was a circle segment between pt and other_pt
//     return (2 * (other_tan * d) / d_len_2) * d - other_tan;
// }


//===========================================================================
void InterpolatedIntersectionCurve::hermite_interpol(list<shared_ptr<IntersectionPoint> >::iterator start_point,
						     list<shared_ptr<IntersectionPoint> >::iterator end_point,
						     Point& mid_position,
						     Point& mid_tangent,
						     EvalKind kind) const 
//===========================================================================
{
    double startpt_info[10]; // nb, curvature and radius not calculated, but set to 0 below
    double endpt_info[10];
    int idim;
    for (int i = 0; i < 10; ++i) {startpt_info[i] = endpt_info[i] = 0;}
    Point pos_1, tng_1, pos_2, tng_2;

    // preparing input data
    switch (kind) {
    case SPACECURVE:
	idim = 3;
	pos_1 = (*start_point)->getPoint();
	tng_1 = tangent_of(start_point);
	pos_2 = (*end_point)->getPoint();
	tng_2 = tangent_of(end_point);
	break;
    case PARAMCURVE_1:
	idim = 2;
	pos_1 = Point((*start_point)->getPar1(), 
		      (*start_point)->getPar1() + (*start_point)->numParams1());
	tng_1 = param_tangent_of(start_point, false);
	pos_2 = Point((*end_point)->getPar1(), 
		      (*end_point)->getPar1() + (*end_point)->numParams1());
	tng_2 = param_tangent_of(end_point, false);
	break;
    case PARAMCURVE_2:
	idim = 2;
	pos_1 = Point((*start_point)->getPar2(), 
		      (*start_point)->getPar2() + (*start_point)->numParams2());
	tng_1 = param_tangent_of(start_point, true);
	pos_2 = Point((*end_point)->getPar2(), 
		      (*end_point)->getPar2() + (*end_point)->numParams2());
	tng_2 = param_tangent_of(end_point, true);
	break;
    default:
	THROW("Logic error in IntersectionCurve::hermite_interpol().");
    }
    ASSERT(pos_1.dimension() == idim && tng_1.dimension() == idim && 
	   pos_2.dimension() == idim && tng_2.dimension() == idim);

    for (int i = 0; i < idim; ++i) {
	startpt_info[i] = pos_1[i];
	endpt_info[i] = pos_2[i];
	startpt_info[i+idim]=tng_1[i];
	endpt_info[i+idim] = tng_2[i];
    }

    mid_position.resize(idim);
    mid_tangent.resize(idim);
    int jstat;
    s1361(startpt_info, 
	  endpt_info, 
	  idim, 
	  mid_position.begin(), 
	  mid_tangent.begin(), 
	  &jstat);

    //double eps = IntersectionPoint::tangent_tol_;
    double eps = (*start_point)->getTolerance()->getEpsge();
    eps *= eps;
    // fixing tangent if zero
    if (mid_tangent.length2() < eps) {
	mid_tangent = pos_2 - pos_1;
	if (mid_tangent.length2() < eps) {
	    // the two points are coincident
	    mid_tangent = tng_1 + tng_2;
	    if (mid_tangent.length2() < eps) {
		// the two input tangents are inversely equal
		mid_tangent = tng_1; // this is our last, desperate estimate
	    }
	} 
    }
    mid_tangent.normalize();
}


//===========================================================================
bool InterpolatedIntersectionCurve::
getGuidePointTangent(shared_ptr<IntersectionPoint> pt, 
		     Point& tan, int type) const
//===========================================================================
{
    list<shared_ptr<IntersectionPoint> >::const_iterator it = find(ipoints_.begin(), ipoints_.end(), pt);
    if (it == ipoints_.end()) {
	return false;
    } 
    switch(type) {
    case 1:
	tan = param_tangent_of(it, false);
	break;
    case 2:
	tan = param_tangent_of(it, true);
	break;
    default:
	tan = tangent_of(it);
    }
    return true;
}


}; // end namespace Go


namespace {


//===========================================================================
void debug_dump_plane(const Point& basept, const Point& normal,
		      double scale, char* name)
//===========================================================================
{
    Point vec(double(0), double(0), double(0));
    int smallest_index = 0;
    double smallest_value = normal[0];
    for (int i = 0; i < 3; ++i) {
	if (normal[i] < smallest_value) {
	    smallest_value = normal[i];
	    smallest_index = i;
	}
    }
    vec[smallest_index] = 1;

    Point dir1 = normal % vec;
    dir1.normalize();
    Point dir2 = normal % dir1;
    dir2.normalize();

    dir1 *= scale;
    dir2 *= scale;

    double corners[3 * 4];

    // filling in corners
    Point tmp = basept + dir1 + dir2;
    copy(tmp.begin(), tmp.end(), corners);
    tmp = basept + dir1 - dir2;
    copy(tmp.begin(), tmp.end(), corners + 3);
    tmp = basept - dir1 + dir2;
    copy(tmp.begin(), tmp.end(), corners + 6);
    tmp = basept - dir1 - dir2;
    copy(tmp.begin(), tmp.end(), corners + 9);

    double knots[] = {0, 0, 1, 1};
    SplineSurface s(2, 2, 2, 2, knots, knots, corners, 3, false);
    
    ofstream os(name);
    s.writeStandardHeader(os);
    s.write(os);
    os.close();
}


//===========================================================================
bool within_tolerances(const Go::Point& p1,  const Go::Point& t1, 
		       const Go::Point& p2,  const Go::Point& t2, 
		       double point_tol, double angle_tol)
//===========================================================================
{
    // checking positions
    if (p1.dist2(p2) > point_tol * point_tol) {
	return false;
    }
    
    if (t1.angle(t2) > angle_tol && t1.angle(-t2) > angle_tol) {
	return false;
    }

    // passed both tests
    return true;
}


//===========================================================================
double angle_with_segment(list<shared_ptr<IntersectionPoint> >::const_iterator p1,
			  list<shared_ptr<IntersectionPoint> >::const_iterator p2,
			  const Go::Point& dir)
//===========================================================================
{
    int dim = (*p1)->getPoint().dimension();
    Go::Point diff;
    if (dim == 3) {
	diff = (*p2)->getPoint() - (*p1)->getPoint();
    } else {
	// parametrical case
	const double* p1_par = (*p1)->getPar1();
	const double* p2_par = (*p2)->getPar1();
	const Go::Point temp1(p1_par, p1_par + (*p1)->numParams1());
	const Go::Point temp2(p2_par, p2_par + (*p1)->numParams1());
	diff = temp2 - temp1;
    }
    return dir.angle(diff);
}


}; // end anonymous namespace
