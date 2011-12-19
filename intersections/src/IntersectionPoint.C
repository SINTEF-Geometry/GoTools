//===========================================================================
//                                                                           
// File: IntersectionPoint.C 
//                                                                           
// Created: 
//                                                                           
// Author: oan
//                                                                           
// Revision: $Id: IntersectionPoint.C,v 1.115 2007-11-01 14:31:44 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/intersections/IntersectionLink.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/intersections/ParamObjectInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamPointInt.h"
#include "sislP.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include <sstream> // debug reasons
#include "GoTools/geometry/SplineCurve.h" // debug reasons
#include <fstream> // debug reasons
#include "GoTools/geometry/PointCloud.h" // debug reasons
#include "GoTools/utils/errormacros.h"
#include <algorithm>
#include <stdexcept>


using std::vector;
using std::ostream;
using std::istream;
using std::ofstream;
using std::string;
using std::stringstream;
using std::cout;
using std::swap;
using std::map;
using std::bad_cast;
using std::runtime_error;
using namespace Go;

namespace { // anonymous namespace


IntPtClassification 
check_against_domain(const Point& p, 
		     const Point& dir, 
		     const RectDomain& dom,
		     bool oriented, // 'dir' gives both direction AND
				    // orientation
		     double res,
		     double angle_tol);

inline bool has_defined_tangent(const SingularityType type);

bool collinear_vectors(const Point& p1, const Point& p2,
		       double tol, double angle_tol);


std::vector<double> 
generate_reduced_param_vec(const shared_ptr<IntersectionPoint> ip,
			   int missing_param);

// debug related function
void dump_tangents(const Point& startpt, const vector<Point>& tangents,
		   string filename);


}; // end anonymous namespace 


namespace Go {


//===========================================================================
const double IntersectionPoint::tangent_tol = 1.0e-7;
//===========================================================================


//===========================================================================
IntersectionPoint::~IntersectionPoint()
//===========================================================================
{
    // disconnect from all other IntersectionPoints
    while (intersection_links_.size() != 0) {
	// removing each element from intersection_links using the
	// disconnectFrom() function
	IntersectionPoint* other
	    = intersection_links_.back()->getOtherPoint(this);
	disconnectFrom(other);
    }
}


//===========================================================================
IntersectionPoint::IntersectionPoint(const ParamObjectInt* obj1, 
				     const ParamObjectInt* obj2,
				     const shared_ptr<GeoTol> epsge,
				     const double* obj1_params,
				     const double* obj2_params)
//===========================================================================
    : obj1_(obj1->getSameTypeAncestor()), 
      obj2_(obj2->getSameTypeAncestor()), 
      parent_point_(),
      epsge_(epsge), 
      g2_discontinuous_params_(detect_2nd_order_discontinuities
			       (obj1, obj2, epsge,
				obj1_params, obj2_params)),
      sec_order_properties_(generate_2nd_order_property_map
			    (g2_discontinuous_params_))
{
    ASSERT(obj1 && obj2 && epsge.get() != 0);
    // Calls constructor_implementation() in order to get around
    // debugging problems
    constructor_implementation(obj1, obj2, epsge,
			       obj1_params, obj2_params);
}


//===========================================================================
void IntersectionPoint::
constructor_implementation(const ParamObjectInt* obj1,
			   const ParamObjectInt* obj2,
			   shared_ptr<GeoTol> epsge,	      
			   const double* obj1_params,
			   const double* obj2_params)
//===========================================================================
{
    // This function makes debugging easier

    int nmb_par1 = obj1->numParams();
    int nmb_par2 = obj2->numParams();
    par_.resize(nmb_par1 + nmb_par2);

    obj1->point(point1_, obj1_params);
    obj2->point(point2_, obj2_params);
    dist_ = point1_.dist(point2_);

    // points in the two objects should be within geometrical
    // tolerance
    if (dist_ >= epsge->getEpsge()) {
// 	std::cout << "Illegal intersection point" << std::endl;
// 	ofstream os("illegalpoint.g2");
// 	double tmp[6];
// 	copy(point1_.begin(), point1_.end(), tmp);
// 	copy(point2_.begin(), point2_.end(), tmp + 3);
// 	PointCloud<3> pcloud(tmp, 2);
// 	pcloud.writeStandardHeader(os);
// 	pcloud.write(os);
// 	os.close();
// 	ASSERT(point1_.dist(point2_) < epsge->getEpsge());
// 	// VSK : [Commented out] to get results
    }

    if (nmb_par1 > 0) {
	par_[0] = obj1_params[0];
	if (nmb_par1 > 1)
	    par_[1] = obj1_params[1];
    }
    if (nmb_par2 > 0) {
	par_[nmb_par1] = obj2_params[0];
	if (nmb_par2 > 1)
	    par_[nmb_par1+1] = obj2_params[1];
    }

    // second-order differential properties should be defined at least
    // for right-side differentiation of all variables.  This test
    // verifies that this 'default' differential property aspect
    // (numbered 0) really exist.  (It should have been generated in
    // 'generate_2nd_order_property_map' - anything else is a bug).
    ASSERT(sec_order_properties_.find(0) != sec_order_properties_.end());

    // setting the currently active 'second order property' aspect to
    // the default one.
    cur_active_2nd_order_properties_ = &sec_order_properties_[0];

//     // additional check to see that parameters are ok (can be removed
//     // when sure everything is working)
//     for (int i = 0; i < nmb_par1; ++i) {
// 	ASSERT(par_[i] >= obj1->startParam(i));
// 	ASSERT(par_[i] <= obj1->endParam(i));
//     }
//     for (int i = 0; i < nmb_par2; ++i) {
// 	ASSERT(par_[i + nmb_par1] >= obj2->startParam(i));
// 	ASSERT(par_[i + nmb_par1] <= obj2->endParam(i));
//     }

    clearCache();
}


//===========================================================================
IntersectionPoint::IntersectionPoint(const ParamObjectInt* obj1,
				     const ParamObjectInt* obj2,
				     const shared_ptr<IntersectionPoint> ip,
				     int missing_param)
//===========================================================================
    : obj1_(obj1->getSameTypeAncestor()), 
      obj2_(obj2->getSameTypeAncestor()), 
      parent_point_(ip),
      epsge_(ip->getTolerance())
{

    par_ = generate_reduced_param_vec(ip, missing_param);
    if (par_.empty()) {
	MESSAGE("Reduced parameter vector is empty...!");
	clearCache();
	return;
    }

    g2_discontinuous_params_
	= detect_2nd_order_discontinuities(obj1, obj2, epsge_,
					   &par_[0], &par_[0] + obj1->numParams());
    sec_order_properties_
	= generate_2nd_order_property_map(g2_discontinuous_params_);



//     // points in the two objects should be within geometrical
//     // tolerance
//     ASSERT(point1_.dist(point2_) < epsge_->getEpsge()); 

    // setting the currently active 'second order property' aspect to
    // the default one.
    cur_active_2nd_order_properties_ = &sec_order_properties_[0];

//     // additional parameter check
//     for (int i = 0; i < obj1->numParams(); ++i) {
// 	ASSERT(par_[i] >= obj1->startParam(i));
// 	ASSERT(par_[i] <= obj1->endParam(i));
//     }
//     int nmb_par1 = obj1->numParams();
//     for (int i = 0; i < obj2->numParams(); ++i) {
// 	ASSERT(par_[i + nmb_par1] >= obj2->startParam(i));
// 	ASSERT(par_[i + nmb_par1] <= obj2->endParam(i));
//     }

    // obj1, obj2 is already used in initializer list so I guess
    // there's no point in calling this ASSERT - should this be obj1_
    // and obj2_ (with underscores)? @jbt
    ASSERT(obj1 && obj2 && epsge_.get() != 0);

    obj1_->point(point1_, getPar1());
    obj2_->point(point2_, getPar2());
    dist_ = point1_.dist(point2_);

    clearCache();
}


//===========================================================================
void IntersectionPoint::write(ostream& os) const
//===========================================================================
{
    const ParamSurfaceInt* s1 = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    const ParamSurfaceInt* s2 = dynamic_cast<const ParamSurfaceInt*>(obj2_);
    
    if (!s1 || !s2) {
	THROW("Could not write IntersectionPoint.  Only implemented for SplineSurfaces.");
    }
    shared_ptr<const SplineSurface> ss1 = 
	dynamic_pointer_cast<const SplineSurface, const ParamSurface>(s1->getParamSurface());
    shared_ptr<const SplineSurface> ss2 = 
	dynamic_pointer_cast<const SplineSurface, const ParamSurface>(s2->getParamSurface());
    if (!ss1.get() || !ss2.get()) {
	THROW("Not able to write IntersectionPoint.  Only implemented for SplineSurfaces.");
    }

    ss1->write(os);
    ss2->write(os);
    epsge_->write(os);
    int num_par = (int)par_.size();
    os << num_par << '\n';
    for (int i = 0; i < num_par; ++i) {
	os << par_[i] << '\n';
    }
    os << point1_.dimension() << '\n';
    point1_.write(os);
    os << '\n';
    os << point2_.dimension() << '\n';
    point2_.write(os);
    os << '\n';
}

//===========================================================================
void IntersectionPoint::read(istream& is)
//===========================================================================
{
    MESSAGE("Warning: IntersectionPoint::read() is dangerous and intended for debugging ONLY!");
    shared_ptr<SplineSurface> s1(new SplineSurface());
    shared_ptr<SplineSurface> s2(new SplineSurface());
    s1->read(is);
    s2->read(is);
    
    obj1_ = new ParamSurfaceInt(s1);
    obj2_ = new ParamSurfaceInt(s2);

    shared_ptr<GeoTol> newepsge(new GeoTol(double(0)));
    newepsge->read(is);
    epsge_.swap(newepsge);
    int num_par, dim;
    is >> num_par;
    par_.resize(num_par);
    for (int i = 0; i < num_par; ++i) {
	is >> par_[i];
    }
    is >> dim;
    point1_.resize(dim);
    point1_.read(is);
    is >> dim;
    point2_.resize(dim);
    point2_.read(is);
    dist_ = point1_.dist(point2_);


    g2_discontinuous_params_ = detect_2nd_order_discontinuities(obj1_,
								obj2_,
								epsge_,
								&par_[0],
								&par_[0] + obj1_->numParams());
    sec_order_properties_ = generate_2nd_order_property_map(g2_discontinuous_params_);
    cur_active_2nd_order_properties_ = &sec_order_properties_[0];
    clearCache();
}


//===========================================================================
void IntersectionPoint::writeParams(ostream& os) const
//===========================================================================
{
    int npar = (int)par_.size();
    int width = (int)os.precision()+3;
    for (int i = 0; i < npar-1; ++i) {
	os.width(width);
	os << par_[i] << " ";
    }
    os.width(width);
    os << par_[npar-1];
    os << " " << dist_;
}


//===========================================================================
void IntersectionPoint::writeInfo() const
//===========================================================================
{
    writeParams(cout);

    // Include more information only for surface-surface
    const ParamSurfaceInt* sf1 = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    const ParamSurfaceInt* sf2 = dynamic_cast<const ParamSurfaceInt*>(obj2_);
    if (sf1 == 0 || sf2 == 0)
	return;

    IntPtInfo info;
    checkIntersectionPoint(info);
    cout << "\t"
	 << "nneigh = " << info.nneighbours << "  "
	 << "sing = " << info.singularity_type << "  "
	 << "dir = " << info.direction << "  "
	 << "loc = " << info.location << "  "
	 << "OK = " << info.is_ok << "    "
	 << "d = " << getDist();

//     Point normal1, normal2;
//     sf1->getParamSurface()->normal(normal1, par_[0], par_[1]);
//     sf2->getParamSurface()->normal(normal2, par_[2], par_[3]);

//     Point tangent = normal1.cross(normal2);
//     double angle_tol_multiplier = 1.0e-8;
//     double tol = angle_tol_multiplier / epsge_->getAngleTol();
//     bool bd1, bd2;
//     isBoundaryPoint(bd1, bd2);
//     if (bd1 || bd2)
// 	tol *= 100.0;
//     double tangent2 = tangent.length2();
//     cout << "  tangent2 = " << tangent2
// 	 << "  tol2 = " << tol * tol;

    return;
}


//===========================================================================
std::vector<bool> IntersectionPoint::
detect_2nd_order_discontinuities(const ParamObjectInt* o1,
				 const ParamObjectInt* o2,
				 const shared_ptr<GeoTol> gtol,
				 const double* obj1_par,
				 const double* obj2_par)
//===========================================================================
{
    // checking if any of the parameters have discontinuous second derivatives at this point

    int n[2];
    n[0] = o1->numParams();
    n[1] = o2->numParams();

    vector<bool> result(n[0] + n[1], false);

    const ParamObjectInt* object[2];
    object[0] = o1;
    object[1] = o2;

    vector<Point> pt_right(6), pt_left(6);
    bool side[4];
    for (int i = 0; i < n[0] + n[1]; ++i) {
	side[i] = true;
    }

    int obj_ix = 0;
    
    const double* par_ptr = obj1_par;
    
    for (int par = 0; par < n[0] + n[1]; ++par) {
	if (par >= n[0]) {
	    obj_ix = 1; 
	    par_ptr = obj2_par;
	}
	object[obj_ix]->point(pt_right, par_ptr, 2, &side[obj_ix*n[0]], 
			      gtol->getRelParRes());
	side[par] = false; // specify that we want to evaluate this parameter from left
	object[obj_ix]->point(pt_left, par_ptr, 2, &side[obj_ix*n[0]],
			      gtol->getRelParRes());
	side[par] = true; // setting back to default

	// calculating position of second derivative
	int p = par - obj_ix * n[0];
	int der2_pos = (n[obj_ix] + 1) * (p + 1) - (p * (p + 1)) / 2;

	// if these vectors are not equal, it means that the second derivative for this
	// parameter is not continuous in this direction (C2-continuous)
	const Point& d2_right = pt_right[der2_pos];
	const Point& d2_left = pt_left[der2_pos];
	if (!collinear_vectors(d2_right, 
			       d2_left, 
			       gtol->getNumericalTol(),
			       gtol->getAngleTol())) {
	    // we found a discontinuity.  The surface's curvature varies according to 
	    // which side you look at it in this parameter direction
	    result[par] = true; 
	}
    }
    return result;
}

//===========================================================================
void IntersectionPoint::flipObjects()
//===========================================================================
{
    int shift = numParams1();
    rotate(par_.begin(), par_.begin() + shift, par_.end());
    rotate(g2_discontinuous_params_.begin(),
	   g2_discontinuous_params_.begin() + shift,
	   g2_discontinuous_params_.end());
    rotate(diff_from_left_.begin(), 
	   diff_from_left_.begin() + shift,
	   diff_from_left_.end());
    rotate(cached_influence_area_forwards_.begin(), 
	   cached_influence_area_forwards_.begin() + shift,
	   cached_influence_area_forwards_.end());
    rotate(cached_influence_area_backwards_.begin(), 
	   cached_influence_area_backwards_.begin() + shift,
	   cached_influence_area_backwards_.end());

    swap(obj1_, obj2_);
    swap(point1_, point2_);

    // reindex elements of property map
    map<int, SecondOrderProperties> reshuffled;
    map<int, SecondOrderProperties>::iterator it;
    vector<bool> tmp(numParams1() + numParams2());
    for (it = sec_order_properties_.begin();
	 it != sec_order_properties_.end(); ++it) {
	int index = it->first; 
	int2bool(index, tmp.begin(), tmp.end());
	rotate(tmp.begin(), tmp.begin() + shift, tmp.end());
	int new_index = bool2int(tmp.begin(), tmp.end());
	reshuffled[new_index] = it->second;
    }
    sec_order_properties_.swap(reshuffled);

    // set pointer to currently active 2nd order property right
    int key = bool2int(diff_from_left_.begin(), diff_from_left_.end());
    cur_active_2nd_order_properties_ = &sec_order_properties_[key];
}


//===========================================================================
void IntersectionPoint::isBoundaryPoint(bool& first, bool& second) const
//===========================================================================
{
    double tol = epsge_->getRelParRes();
    int nmb_par1 = obj1_->numParams();

    first = obj1_->boundaryPoint(&par_[0], tol);
    second = obj2_->boundaryPoint(&par_[nmb_par1], tol);
}


//===========================================================================
void IntersectionPoint::
tangentPointingInwards(bool& first, 
		       bool& second,
		       bool& first_along_boundary,
		       bool& second_along_boundary) const
//===========================================================================
{
//     SingularityType st = getSingularityType();
    DEBUG_ERROR_IF(!has_defined_tangent(getSingularityType()),
		   "Called tangentPointingInwards() on an "
		   "IntersectionPoint with no defined tangent.");

    const ParamSurfaceInt* s1 = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    const ParamSurfaceInt* s2 = dynamic_cast<const ParamSurfaceInt*>(obj2_);
    shared_ptr<const ParamSurface> sf1;
    shared_ptr<const ParamSurface> sf2;
    if (!s1) {
	const Param2FunctionInt* pf
	    = dynamic_cast<const Param2FunctionInt*>(obj1_);
      if (pf)
	sf1 = pf->getParamSurface();
    } else {
      sf1 = s1->getParamSurface();
    }
    if (s2) {
      sf2 = s2->getParamSurface();
    }
    DEBUG_ERROR_IF(!sf1 && !sf2,
		   "IntersectionPoint::tangentPointingInwards() "
		   "only implemented for parametric surface "
		   "intersections for now.");

    bool oriented = !pointIsSingular();

    const double& angle_tol = epsge_->getAngleTol();
    if (sf1) {
        RectDomain d1;
	try {
	    d1 = dynamic_cast<const RectDomain&>(sf1->parameterDomain());
	} catch(bad_cast) {
	  THROW("Too complicated boundary representation (trimmed?) "
		"for the current implementation of "
		"IntersectionPoint::tangentPointingInwards()");
	}

	Point p1(getPar1(), getPar1() + numParams1());
	Point dir1 = getPar1Dir();
	IntPtClassification cl1 = 
	    check_against_domain(p1, dir1, d1, oriented,
				 epsge_->getRelParRes(), angle_tol);

	first = (cl1 == DIR_UNDEF) || (cl1 == DIR_IN);
	first_along_boundary = (cl1 == DIR_PARALLEL);
    }

    if (sf2) {
        RectDomain d2;
	try {
	    d2 = dynamic_cast<const RectDomain&>(sf2->parameterDomain());
	} catch(bad_cast) {
	  THROW("Too complicated boundary representation (trimmed?) "
		"for the current implementation of "
		"IntersectionPoint::tangentPointingInwards()");
	}

	Point p2(getPar2(), getPar2() + numParams2());
	Point dir2 = getPar2Dir();
	IntPtClassification cl2 = 
	    check_against_domain(p2, dir2, d2, oriented,
				 epsge_->getRelParRes(), angle_tol);

	second = (cl2 == DIR_UNDEF) || (cl2 == DIR_IN);
	second_along_boundary = (cl2 == DIR_PARALLEL);
    }
}

// //===========================================================================
// void IntersectionPoint::setValues(const Point& point1, const Point& point2,
// 				  double* pointpar1, double* pointpar2)
// //===========================================================================
// {
//   point1_ = point1;
//   point2_ = point2;
 
//  // @bsp Assumes objects are surfaces with two parameter directions
//   par_[0] = pointpar1[0];
//   par_[1] = pointpar1[1];
//   par_[2] = pointpar2[0];
//   par_[3] = pointpar2[1];

//   clearCache();
// }

//===========================================================================
void IntersectionPoint::replaceParameter(double *param)
//===========================================================================
{
    int nmb_par1 = obj1_->numParams();
    int nmb_par2 = obj2_->numParams();
    //ASSERT(nmb_par1+nmb_par2 == npar);
    int npar = nmb_par1 + nmb_par2;

    int ki;
    for (ki=0; ki<npar; ki++)
	par_[ki] = param[ki];

    obj1_->point(point1_, param);
    obj2_->point(point2_, param+nmb_par1);
    dist_ = point1_.dist(point2_);
   
    g2_discontinuous_params_ = 
	detect_2nd_order_discontinuities(obj1_, 
					 obj2_, 
					 epsge_, 
					 &par_[0], 
					 &par_[obj1_->numParams()]);
    sec_order_properties_ = generate_2nd_order_property_map(g2_discontinuous_params_);
    cur_active_2nd_order_properties_ = &sec_order_properties_[0];
    clearCache();
}

//===========================================================================
Point IntersectionPoint::getPoint() const
//===========================================================================
{
  Point pnt(point1_);
  pnt += point2_;
  pnt *= 0.5;
  return pnt;
}

//===========================================================================
Point IntersectionPoint::getPar1Dir(bool second_branch) const
//===========================================================================
{
    SingularityType st;
    st = getSingularityType();
    DEBUG_ERROR_IF(!has_defined_tangent(st),
		   "Called getPar1Dir() on an IntersectionPoint"
		   " with no defined tangent.");
    DEBUG_ERROR_IF(second_branch && st != BRANCH_POINT,
	     "Error in getPar1Dir(): tried to get a second branch out of an "
	     "IntersectionPoint that is NOT a branch point.");

    if (!tangent_2d_is_cached()) {
	compute_parametrical_tangents();
    }

    ASSERT(tangent_2d_is_cached());
    return cur_active_2nd_order_properties_->tangent_2d_1[int(second_branch)];
}

//===========================================================================
Point IntersectionPoint::getPar2Dir(bool second_branch) const
//===========================================================================
{
    SingularityType st;
    st = getSingularityType();
    DEBUG_ERROR_IF(!has_defined_tangent(st),
		   "Called getPar2Dir() on an IntersectionPoint"
		   " with no defined tangent.");
    DEBUG_ERROR_IF(second_branch && st != BRANCH_POINT,
	     "Error in getPar2Dir(): tried to get a second branch out of an "
	     "IntersectionPoint that is NOT a branch point.");

    if (!tangent_2d_is_cached()) {
	compute_parametrical_tangents();
    } 
    ASSERT(tangent_2d_is_cached());
    return cur_active_2nd_order_properties_->tangent_2d_2[int(second_branch)];
}

//===========================================================================
Point IntersectionPoint::getTangent(bool second_branch) const
//===========================================================================
{
    SingularityType st = getSingularityType();
    if (!has_defined_tangent(st)) {
	throw runtime_error("Tried to run getTangent() on an IntersectionPoint "
			    "with no possible defined tangent.");
    }
    DEBUG_ERROR_IF(second_branch && st != BRANCH_POINT, 
	     "Error in getTangent(): tried to get a second branch out of an "
	     "IntersectionPoint that is NOT a branch point.");
    ASSERT(singularity_info_is_cached()); // should have been taken care of above

    // For the Param2FunctionInt case we are interested in the 2D intersection
    // point between the object and the plane.
    const Param2FunctionInt* pf = dynamic_cast<const Param2FunctionInt*>(obj1_);
    if (pf) {
	return cur_active_2nd_order_properties_->tangent_2d_1[int(second_branch)];
    } else {
	return cur_active_2nd_order_properties_->tangent_3d[int(second_branch)];
    }
}

// //===========================================================================
// Point IntersectionPoint::getDomainSpacePoint(int obj_ind) const
// //===========================================================================
// {
//     int nmb_par1 = numParams1();
//     int nmb_par2 = numParams1();
//     if (obj_ind == 0 && nmb_par1 > 0) {
// 	Point pnt(nmb_par1 + 1);
// 	for (int ki = 0; ki < nmb_par1; ++ki)
// 	    pnt[ki] = par_[ki];
// 	for (int ki = nmb_par1; ki < point1_.size(); ++ki)
// 	    pnt[nmb_par1+ki] = point1_[ki];
// 	return pnt;
//     } else if (obj_ind == 1 && nmb_par2 > 0) {
// 	Point pnt(nmb_par2 + 1);
// 	for (int ki = 0; ki < nmb_par2; ++ki)
// 	    pnt[ki] = par_[ki];
// 	for (int ki = nmb_par2; ki < point2_.size(); ++ki)
// 	    pnt[nmb_par2+ki] = point2_[ki];
// 	return pnt;
//     } else {
// 	GO_ERROR("Unexpected incident occured!", UnknownError());
//     }
// }

// //===========================================================================
// Point IntersectionPoint::getDomainSpaceTangent(int obj_ind) const
// //===========================================================================
// {
//     Point space_tangent = getTangent();
//     int nmb_par1 = numParams1();
//     int nmb_par2 = numParams1();
//     if (obj_ind == 0 && nmb_par1 > 0) {
// 	Point tangent(nmb_par1 + 1);
// 	Point par1_tangent = getPar1Dir();
// 	for (int ki = 0; ki < nmb_par1; ++ki)
// 	    tangent[ki] = par1_tangent[ki];
// 	for (int ki = nmb_par1; ki < space_tangent.size(); ++ki)
// 	    tangent[nmb_par1+ki] = space_tangent[ki];
// 	return tangent;
//     } else if (obj_ind == 1 && nmb_par2 > 0) {
// 	Point tangent(nmb_par2 + 1);
// 	Point par2_tangent = getPar2Dir();
// 	for (int ki = 0; ki < nmb_par2; ++ki)
// 	    tangent[ki] = par2_tangent[ki];
// 	for (int ki = nmb_par2; ki < space_tangent.size(); ++ki)
// 	    tangent[nmb_par2+ki] = space_tangent[ki];
// 	return tangent;
//     } else {
// 	GO_ERROR("Unexpected incident occured!", UnknownError());
//     }
// }

//===========================================================================
void IntersectionPoint::compute_parametrical_tangents() const
//===========================================================================
{
    if (dynamic_cast<const Param2FunctionInt*>(obj1_)) {
	// we are dealing with a parametric function
	compute_function_parametrical_tangents();
    } else if (dynamic_cast<const ParamSurfaceInt*>(obj1_) &&
	       dynamic_cast<const ParamSurfaceInt*>(obj2_)) {
	// we are dealing with two intersecting surfaces
	compute_surface_parametrical_tangents_safe();
    } else {
	// unknown object combination type
	THROW("calculate_function_parametrical_tangents failed due to "
	      "unrecognized type.");
    }
    cur_active_2nd_order_properties_->tangent_2d_is_cached = true; 
}


//===========================================================================
void IntersectionPoint::compute_function_parametrical_tangents() const
//===========================================================================
{
    const Param2FunctionInt* p2f
	= dynamic_cast<const Param2FunctionInt*>(obj1_);
    DEBUG_ERROR_IF(p2f == 0, "Intersection scenario not handled!");
    Point* tangent_1 = cur_active_2nd_order_properties_->tangent_2d_1;

    // The Param2FunctionInt does not have a normal, but we may look
    // at it's partial derivs.
    vector<Point> pt(6);
    shared_ptr<const ParamSurface> sf = p2f->getParamSurface();

    // It is vital that the derivatives are computed from the right
    // side!
    ASSERT(diff_from_left_.size() >= 2); // Should really be 2.
    sf->point(pt, par_[0], par_[1], 2,
	      !diff_from_left_[0], !diff_from_left_[1]);

    // Let s(u,v) = (u,v,sf(u,v))
    // s-normal = (1, 0, sf_u) x (0, 1, sf_v)
    Point sf_normal(-pt[1][0], pt[2][0], 1.0);
    if (fabs(pt[1][0]) < tangent_tol && fabs(pt[2][0]) < tangent_tol) {
	MESSAGE("Not yet stable ...");
	// We should calculate the direction based on curvature.
	// Either it is an isolated point, or the curvature in a dir
	// is 0.0.
 	// Tangential intersection.
	// We set the curvature to 0.0 and use the abc-formula to
	// solve for q.
	double tangent2_tol = tangent_tol;
	if (fabs(pt[3][0]) < tangent2_tol || fabs(pt[5][0]) < tangent2_tol) {
	    ASSERT(fabs(pt[3][0]) > tangent2_tol
		   || fabs(pt[5][0]) > tangent2_tol);
	    tangent_1[0] = (fabs(pt[3][0]) < tangent2_tol)
		? Point(1.0, 0.0) : Point(0.0, 1.0);
	} else {
	    // @@sbr Well, then we should never get here ...
	    double p = 1.0;
	    double a = pt[5][0];
	    double b = 2*p*pt[4][0];
	    double c = p*p*pt[3][0];
	    double non_neg = b*b - 4*a*c;
	    if (non_neg < 0.0) {
		// Otherwise it is an isolated sing.
		THROW("compute_function_parametrical_tangents: "
		      "Trying to compute square root of neg value)!");
	    }
	    double q1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
// 	    double q2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
	    tangent_1[0] = Point(p, q1);
	}
    } else {
	tangent_1[0] = Point(pt[2][0], -pt[1][0]);
    }
    // We must find the direction in which the derivative is <
    // tangent2_tol.
    // The correct would be to find the direciton in which the
    // curvature is 0.0.
    // D_u(f(x,y)) = p*f_x(x,y) + q*f_y(x,y), for u = <p,q>

    // D_uu(f(x,y)) = p^2*f_xx(x,y) + 2pq*f_xy(x,y) + q^2*f_yy(x,y),
    // for u = <p,q>

//     // We expect at least one of the 2nd order partial derivatives to
//     // differ from 0.0.
//     ASSERT(fabs(pt[3][0]) > tangent2_tol || fabs(pt[5][0]) > tangent2_tol);
//     if (pt[1].length2() < tangent2_tol && pt[2].length2() < tangent2_tol) {
// 	GO_WARNING("Not yet fully supported ...");
// 	// Tangential intersection.
// 	// We set the curvature to 0. and use the abc-formula to solve for q.
// 	double p = 1.0;
// 	double a = pt[5][0];
// 	double b = 2*p*pt[4][0];
// 	double c = p*p*pt[3][0];
// 	double q1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
// 	double q2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
// 	tangent_1[0] = Point(p, q1);
//     } else {
// 	// Tranversial intersection.
// 	if (pt[1].length2() < tangent2_tol || pt[2].length2() < tangent2_tol) {
// 	    tangent_1[0] = (fabs(pt[1][0]) > tangent2_tol)
// 		? Point(1.0, 0.0) : Point(0.0, 1.0);
// 	} else {
// 	    double p = 1.0;
// 	    double q = -p*pt[1][0]/pt[2][0];
// 	    tangent_1[0] = Point(p, q);
// 	}
//     }

    // Dir of int_cv is given by
    // Sf_n x C_n = (-S_u, S_v, 1) x (0, 0, 1) = (S_v, -S_u, 0)
    // Sf_n x C_n = (-S_u, -S_v, 1) x (0, 0, 1) = (-S_v, S_u, 0)
    // @@sbr Must make sure that IN-point points into domain (and the
    // opp for OUT-point).
    tangent_1[0].normalize();
    tangent_1[1] = Point(double(0), double(0));
//     cur_active_2nd_order_properties_->tangent_is_oriented = true;
// =======
//     double tangent2_tol = tangent_tol;
//     if (pt[3].length2() > tangent2_tol || pt[5].length2() > tangent2_tol) {
// 	ASSERT(pt[3][0] < tangent2_tol || pt[5][0] < tangent2_tol);
// 	// Dir of int_cv is given by
// 	// Sf_n x C_n = (-S_u, S_v, 1) x (0, 0, 1) = (S_v, -S_u, 0)
// 	// Sf_n x C_n = (-S_u, -S_v, 1) x (0, 0, 1) = (-S_v, S_u, 0) // @@sbr Or ...
// 	// @@sbr Must make sure that IN-point points into domain (and the opp for OUT-point).
// 	tangent_1[0] = (pt[3][0] > tangent2_tol) ? Point(1.0, 0.0) : Point(0.0, 1.0);
// 	tangent_1[0].normalize();
// 	tangent_1[1] = Point(double(0), double(0));
// 	cur_active_2nd_order_properties_->tangent_is_oriented = true;
// >>>>>>> 1.57
	// @@sbr The orientation is consistent, but dependent on correct ordering of pts.
	// Should be accomplished by deciding whether the par tangent points IN or OUT.
//     } else {
// 	GO_ERROR("Unexpected incident (sf degenerate).", InputError());
//     }
}

//===========================================================================
bool IntersectionPoint::tangentIsOriented() const
//===========================================================================
{
    if (!hasUniqueTangentDirection()) {
	return false; // cannot be oriented if it does not even have a
		      // direction
    }

//     SingularityType st = getSingularityType();
    DEBUG_ERROR_IF(!has_defined_tangent(getSingularityType()),
		   "Called tangentIsOriented() on an IntersectionPoint"
		   " with no defined tangent.");
    ASSERT(singularity_info_is_cached());
    return cur_active_2nd_order_properties_->tangent_is_oriented;
}
 
//===========================================================================
void IntersectionPoint::fixTangentOrientation(bool flip) 
//===========================================================================
{
    // For a 2par function it may be useful to flip direction.
    DEBUG_ERROR_IF(!hasUniqueTangentDirection(), 
	     "Can only define a tangent orientation for points with uniquely defined tangents.");
    if (flip) {
	// flipping directions of tangent
	ASSERT(singularity_info_is_cached()); // taken care of in getSingularityType()
	cur_active_2nd_order_properties_->tangent_3d[0] *= -1;
	if (tangent_2d_is_cached()) {
	    cur_active_2nd_order_properties_->tangent_2d_1[0] *= -1;
	    cur_active_2nd_order_properties_->tangent_2d_2[0] *= -1;
	}
    }

    // from now on, the tangent is oriented
    cur_active_2nd_order_properties_->tangent_is_oriented = true;
}

//===========================================================================
void IntersectionPoint::compute_surface_parametrical_tangents_safe() const 
//===========================================================================
{
    Point udir_s1, vdir_s1, udir_s2, vdir_s2;
    const ParamSurfaceInt* s1 = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    const ParamSurfaceInt* s2 = dynamic_cast<const ParamSurfaceInt*>(obj2_);

    DEBUG_ERROR_IF(!s1 || !s2,
	     "This function was called under the supposition that we were dealing"
	     " with parametric surfaces, but another object kind was found.");


    s1->derivs(par_[0], par_[1], udir_s1, vdir_s1, !diff_from_left_[0], !diff_from_left_[1]);
    s2->derivs(par_[2], par_[3], udir_s2, vdir_s2, !diff_from_left_[2], !diff_from_left_[3]);

    SingularityType st = getSingularityType(); 
    Point* tangent_1 = cur_active_2nd_order_properties_->tangent_2d_1;
    Point* tangent_2 = cur_active_2nd_order_properties_->tangent_2d_2;
    if (!has_defined_tangent(st)) {
	// we have no hope of finding any sensible tangent for this point
	tangent_1[0] = tangent_2[0] = Point(double(0), double(0));
    } else {
	int num_loops = (st == BRANCH_POINT) ? 2 : 1;
	for (int i = 0; i < num_loops; ++i) {
	    Point tang = getTangent(i != 0);
	    
	    // calculating parametrical tangent in first surface
	    tangent_1[i].resize(2);
	    tangent_2[i].resize(2);
	    decompose(tang, udir_s1, vdir_s1, tangent_1[i][0], tangent_1[i][1]);
	    decompose(tang, udir_s2, vdir_s2, tangent_2[i][0], tangent_2[i][1]);
	    
	    if (tangent_1[i].length() > epsge_->getRelParRes())
		tangent_1[i].normalize();
	    if (tangent_2[i].length() > epsge_->getRelParRes())
		tangent_2[i].normalize();
	}
    }
 
    cur_active_2nd_order_properties_->tangent_2d_is_cached = true; 
}

//===========================================================================
void IntersectionPoint::projectToParamPlanes(const Point& dir, 
					     Point& par_1_dir, 
					     Point& par_2_dir) const
//===========================================================================
{
    const ParamGeomInt* o1 = dynamic_cast<const ParamGeomInt*>(obj1_);
    const ParamGeomInt* o2 = dynamic_cast<const ParamGeomInt*>(obj2_);
    DEBUG_ERROR_IF(!o1 || !o2,
	     "Function IntersectionPoint::projectToParamPlanes() currently only works "
	     "for geometrical objects, while user tried to apply it on an "
	     "IntersectionPoint refering to objects of a different kind.");
    const double* param = &par_[0];
    vector<bool>::const_iterator from_left = diff_from_left_.begin();
    proj_2_pplane(dir, param, from_left, o1, par_1_dir);
    param += o1->numParams();
    from_left += o1->numParams();
    proj_2_pplane(dir, param, from_left, o2, par_2_dir);
}


//===========================================================================
void IntersectionPoint::decompose(const Point& T,
				  const Point& b1,
				  const Point& b2,
				  double& b1_coef,
				  double& b2_coef)
//===========================================================================
{
    // We suppose that 'T' lies in the plane spanned by 'b1' and 'b2'.  
    // We want to calculate the coefficients in order to express 'T' as a 
    // linear combination of 'b1' and 'b2'.
    // We solve the overdetermined system M * x = T, where x are the two coefficients
    // we are looking for.  We obtain x = (MtM)^{-1} MtT

    double MtM_inv[2][2];  // inverse of matrix MtM = | u.u  u.v |
                           //                         | u.v  v v |
                           // which is (using D for determinant):
                           //   1/D * |  v.v  -u.v |
                           //         | -u.v   u.u |
    double b1_norm2 = b1.length2();
    double b2_norm2 = b2.length2();
    double b1b2_scl = b1 * b2;

    double D_inv = double(1) / (b1_norm2 * b2_norm2 - (b1b2_scl * b1b2_scl));

    MtM_inv[0][0] = b2_norm2 * D_inv;    // adjunct matrix divided by determinant
    MtM_inv[0][1] = -b1b2_scl * D_inv;
    MtM_inv[1][0] = -b1b2_scl * D_inv;
    MtM_inv[1][1] = b1_norm2 * D_inv;

    double MtT[2];
    MtT[0] = b1 * T;
    MtT[1] = b2 * T;

    b1_coef = MtM_inv[0][0] * MtT[0] + MtM_inv[0][1] * MtT[1];
    b2_coef = MtM_inv[1][0] * MtT[0] + MtM_inv[1][1] * MtT[1];
    
}

//===========================================================================
void IntersectionPoint::calculate_singularity_type() const 
//===========================================================================
{
    if (dynamic_cast<const Param2FunctionInt*>(obj1_)) { 

	// we are dealing with a parametric function
	calculate_singularity_type_function();

    } else if (dynamic_cast<const ParamSurfaceInt*>(obj1_) && 
	       dynamic_cast<const ParamSurfaceInt*>(obj2_)) {

	// we are dealing with two parametric surfaces
	calculate_singularity_type_surface();

    } else if (dynamic_cast<const ParamCurveInt*>(obj1_) &&
	       dynamic_cast<const ParamCurveInt*>(obj2_)) {

	// we are dealing with two parametric curves
	calculate_singularity_type_curves();

    } else if (dynamic_cast<const ParamPointInt*>(obj1_) ||
	       dynamic_cast<const ParamPointInt*>(obj2_)) {
	
	// one of the objects was a point
	calculate_singularity_type_point();

    } else if ((dynamic_cast<const ParamSurfaceInt*>(obj1_) &&
		dynamic_cast<const ParamCurveInt*>(obj2_)) ||
	       (dynamic_cast<const ParamCurveInt*>(obj1_) &&
		dynamic_cast<const ParamSurfaceInt*>(obj2_))) {

	// we have one curve and one surface
	calculate_singularity_type_curve_surface();

    } else {
	THROW("Unrecognized object combination encountered in "
	      "calculate_singularity_type()");
    }
}

//===========================================================================
void IntersectionPoint::calculate_singularity_type_function() const  
//===========================================================================
{
    // We must calculate 2d-parametric tangent.
    if (!tangent_2d_is_cached()) {
	compute_parametrical_tangents();
    }

    // We must calculate partial derivatives in the function.
    const Param2FunctionInt* p2f_int = dynamic_cast<const Param2FunctionInt*>(obj1_);
    if (p2f_int != 0) {
	shared_ptr<const ParamSurface> sf = p2f_int->getParamSurface();
	vector<Point> pts(6);
	sf->point(pts, par_[0], par_[1], 2);
// 	double length = sqrt(pts[1].length2() + pts[2].length2());
	double tangent2_tol = tangent_tol;
// 	ASSERT(fabs(pts[3][0]) < tangent2_tol || fabs(pts[5][0]) < tangent2_tol);
	if (pts[3].length() < tangent2_tol && pts[5].length() < tangent2_tol) {
// 	if (pts[3].length2() < tangent2_tol && pts[5].length2() < tangent2_tol) {
	    // If both partial derivs are (almost-)0 the point is an isolated point
	    // (singularity).
	    cur_active_2nd_order_properties_->singularity_type = ISOLATED_POINT;
	    cur_active_2nd_order_properties_->tangent_is_oriented = false;
	} else {
	    // We then expect the point to be a an ordinary point. //tangential point.
	    // The direction of the intersection curve is then given by the cross
	    // product of the function with the plane (when looked upon in 3d).
	    cur_active_2nd_order_properties_->singularity_type = ORDINARY_POINT; 
	    cur_active_2nd_order_properties_->tangent_is_oriented = false; //true;
//	    cur_active_2nd_order_properties_->tangent_is_oriented = true;
	}
    } else { // We should also support 
	// Currently this branch handles intersection between 2-parametric function 
	// and constant.
	THROW("Not implemented yet!");
    }

    // @@sbr Are we not setting singularity_info_is_cached?
    cur_active_2nd_order_properties_->singularity_info_is_cached = true;
}

//===========================================================================
void IntersectionPoint::calculate_singularity_type_point() const
//===========================================================================
{
    // this function is called if one of the objects is a single point.
    // This means that our IntersectionPoint must necessarily be isolated.
    cur_active_2nd_order_properties_->tangent_3d[0] = 
	Point(double(0), double(0), double(0)) ;
    cur_active_2nd_order_properties_->singularity_type = ISOLATED_POINT;
    cur_active_2nd_order_properties_->tangent_is_oriented = false;
    cur_active_2nd_order_properties_->singularity_info_is_cached = true;
}

//===========================================================================
void IntersectionPoint::calculate_singularity_type_curves() const
//===========================================================================
{
    const ParamCurveInt* c1 = dynamic_cast<const ParamCurveInt*>(obj1_);
    const ParamCurveInt* c2 = dynamic_cast<const ParamCurveInt*>(obj2_);
    ASSERT(c1 && c2); // if not true, function should not have been called in the first place
    
    vector<Point> pts1(2), pts2(2); // this is a pain, but interface of ParamCurve requires it

    Point& tan3d_1 = cur_active_2nd_order_properties_->tangent_3d[0];
    Point& tan3d_2 = cur_active_2nd_order_properties_->tangent_3d[1];

    c1->getParamCurve()->point(pts1, par_[0], 1, !diff_from_left_[0]);
    c2->getParamCurve()->point(pts2, par_[1], 1, !diff_from_left_[1]);
    
    // calculate angle here
    double angle = 0;
    double tan_len1 = pts1[1].length();
    double tan_len2 = pts2[1].length();
    const double num_tol = epsge_->getNumericalTol();
    const double angle_tol = epsge_->getAngleTol();
    if (tan_len1 < num_tol || tan_len2 < num_tol) {
	// One of the curves are degenerate in this point (has practically zero tangent).
	// We can consider that they are "perpendicular", and treat this intersection
	// as transversal.  There we (artificially) set the angle larger than the
	// angle tolerance
	angle = 2 * angle_tol;
    } else {
	// normalizing tangents
	pts1[1] /= tan_len1;
	pts2[1] /= tan_len2;
	angle = pts1[1].cross(pts2[1]).length(); // OK, this is really the sine of the angle, 
	                                         // but for our purposes, that is good enough.
    }
    
    if (angle > angle_tol) {

	// the curves are intersecting each other transversally, the point is isolated
	tan3d_1 = tan3d_2 = Point(double(0), double(0), double(0)); 
	cur_active_2nd_order_properties_->singularity_type = ISOLATED_POINT;
	cur_active_2nd_order_properties_->tangent_is_oriented = false;

    } else {

	// the curves are intersecting each other tangentially
	tan3d_1 = pts1[1]; // could equally use pts2[1], and orientation is not important here.
	tan3d_2 = Point(double(0), double(0), double(0)); // (only relevant for branch points)
	cur_active_2nd_order_properties_->singularity_type = TANGENTIAL_POINT;
	cur_active_2nd_order_properties_->tangent_is_oriented = false;

    }
    cur_active_2nd_order_properties_->singularity_info_is_cached = true;
}

//===========================================================================
void IntersectionPoint::calculate_singularity_type_curve_surface() const
//===========================================================================
{
    const ParamSurfaceInt* sf_obj = 0;
    const ParamCurveInt* cv_obj = dynamic_cast<const ParamCurveInt*>(obj1_);
    bool first_is_curve = (cv_obj != 0);
    if (first_is_curve) {
	sf_obj = dynamic_cast<const ParamSurfaceInt*>(obj2_);
    } else {
	cv_obj = dynamic_cast<const ParamCurveInt*>(obj2_);
	sf_obj = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    }
    ASSERT(sf_obj && cv_obj); // if not, this function should not have been called
    const int cv_p_ix = (first_is_curve) ? 0 : 2;
    const int sf_p_ix = (first_is_curve) ? 1 : 0;

    // check if the tangent of the curve lies in the surface tangent plane
    vector<Point> pts(2);
    cv_obj->getParamCurve()->point(pts, par_[cv_p_ix], 1, !diff_from_left_[cv_p_ix]);
    
    Point sf_normal;
    const double num_tol = epsge_->getNumericalTol();
    const double angle_tol = epsge_->getAngleTol();
    bool right1 = !diff_from_left_[sf_p_ix];
    bool right2 = !diff_from_left_[sf_p_ix + 1];
    sf_obj->normal(par_[sf_p_ix], par_[sf_p_ix + 1], sf_normal, right1, right2);
    
    double sf_normal_length = sf_normal.length();
    double cv_tan_length = pts[1].length();
    if (sf_normal_length > num_tol && cv_tan_length > num_tol) {
	sf_normal /= sf_normal_length;
	pts[1] /= cv_tan_length;
	if (sf_normal * pts[1] < angle_tol) {
	    // we conclude that the curve tangent is virutally in the surface tangent plane
	    cur_active_2nd_order_properties_->tangent_3d[0] = pts[1];
	    cur_active_2nd_order_properties_->tangent_3d[1] = 
		Point(double(0), double(0), double(0));
	    cur_active_2nd_order_properties_->singularity_type = TANGENTIAL_POINT;
	    cur_active_2nd_order_properties_->tangent_is_oriented = false;
	    cur_active_2nd_order_properties_->singularity_info_is_cached = true;
	    return;
	}
    }
    // either the surface normal is degenerated, the curve tangent is degenerated
    // or the curve passes through the surface in a transversal intersection
    cur_active_2nd_order_properties_->tangent_3d[0] = 
	cur_active_2nd_order_properties_->tangent_3d[1] = 
	Point(double(0), double(0), double(0));
    cur_active_2nd_order_properties_->singularity_type = ISOLATED_POINT;
    cur_active_2nd_order_properties_->tangent_is_oriented = false;
    cur_active_2nd_order_properties_->singularity_info_is_cached = true;
    return;
}


//===========================================================================
void IntersectionPoint::calculate_singularity_type_surface() const
//===========================================================================
{
    // estimated from stability requirements, see external notes. Test
    // VSK
//     static const double angle_tol_multiplier = 1.0e-4;
//     static const double angle_tol_multiplier = 1.0e-6; 
    static const double angle_tol_multiplier = 1.0e-8; 

    const ParamSurfaceInt* s1 = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    const ParamSurfaceInt* s2 = dynamic_cast<const ParamSurfaceInt*>(obj2_);
    ASSERT(s1 && s2); // if not true, function should not have been
		      // called in the first place

    Point normal_1, normal_2;
    s1->getParamSurface()->normal(normal_1, par_[0], par_[1]);
    s2->getParamSurface()->normal(normal_2, par_[2], par_[3]);
    Point& tan3d_1 = cur_active_2nd_order_properties_->tangent_3d[0];
    Point& tan3d_2 = cur_active_2nd_order_properties_->tangent_3d[1];
    tan3d_1 = normal_1.cross(normal_2);

    double tol = angle_tol_multiplier / epsge_->getAngleTol();

    // VSK 0806. Check if the point is a global boundary point 
    // In that case the singularity criterion must be slightly relaxed
    bool bd1, bd2;
    isBoundaryPoint(bd1, bd2);
    if (bd1 || bd2)
	tol *= 100.0;

    if (tan3d_1.length2() > tol * tol) {
	// normal, transversal situation (not singular)
	tan3d_1.normalize();
	tan3d_2 = Point(double(0),double(0),double(0)); // unused
	cur_active_2nd_order_properties_->singularity_type = ORDINARY_POINT;
	cur_active_2nd_order_properties_->tangent_is_oriented = true;
    } else {
	// the computed tangent is so small that its direction cannot
	// be trusted within the required tolerance.  We must
	// determine it through other means.

	// @@ debug-related - we want to dump the various tangents
	// found for this point to a file.

	// ----------------------------- DEBUG RELATED STARTS HERE
//        // This needs testing! @jbt
//	static int report_number = 0;
//        stringstream snumber; 
//        snumber << report_number++;
//	string dump_file = "sing_pt_"
//	    + snumber.str() + ".g2";
//	cout << "Singular intersection point detected.  "
//	     << "Dumping tangents to file " 
//	     << dump_file << endl;
//	vector<Point> tangents;
//	tangents.push_back(tan3d_1);
//	if (tangents[0].length() > 1.0e-10)
//	    tangents[0].normalize();
	// ----------------------------- DEBUG RELATED ENDS HERE

	cur_active_2nd_order_properties_->singularity_type
	    = calculate_tangent_at_singular_point();
// 	    = HIGHER_ORDER_POINT; 
	cur_active_2nd_order_properties_->tangent_is_oriented = false;

	// ----------------------------- DEBUG RELATED STARTS HERE
// 	tangents.push_back(tan3d_1);
// 	if (cur_active_2nd_order_properties_->singularity_type
// 	    == BRANCH_POINT) {
// 	    tangents.push_back(tan3d_2);
// 	}
// 	dump_tangents(getPoint(), tangents, dump_file);
	// ----------------------------- DEBUG RELATED ENDS HERE

    }
    cur_active_2nd_order_properties_->singularity_info_is_cached = true;
}


//===========================================================================
SingularityType IntersectionPoint::calculate_tangent_at_singular_point() const
//===========================================================================
{
    // For theoretical explanation of the routine, refer to "Shape
    // Interrogation for Computer Aided Design and Manufacturing" by
    // Patrikalakis/Maekawa, Section 6.4.1, page 171.
    const double EPS = 1.0e-5;

    const ParamSurfaceInt* s1 = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    const ParamSurfaceInt* s2 = dynamic_cast<const ParamSurfaceInt*>(obj2_);
    DEBUG_ERROR_IF(s1 == 0 || s2 == 0,
		   "IntersectionPoint::calculate_tangent_at_singular_point() "
		   "is currently only implemented for parametric surface "
		   "intersections.");
    
    // Get normals (normalized)
    Point N1, N2;
    s1->getParamSurface()->normal(N1, par_[0], par_[1]);
    s2->getParamSurface()->normal(N2, par_[2], par_[3]);

    // We need the same orientation for the formulas below to
    // work. That is, the two normals must be equal. If they are not,
    // then the normals point in the opposite directions, and we must
    // apply some sign changes a few places. The sign changes are
    // supposed to have the same effect as reversing the u-parameter
    // direction of surface B.
    bool opposite_orientation = false;
    Point Ndiff = N1 - N2; // Equals either 0.0 or 2.0
    if (Ndiff.length() > 1.0) {
	opposite_orientation = true;
    }

    // calculating derivatives along isolines
    Point rA_u, rA_v, rB_u, rB_v;
    s1->derivs(par_[0], par_[1], rA_u, rA_v,
	       !diff_from_left_[0], !diff_from_left_[1]);
    s2->derivs(par_[2], par_[3], rB_u, rB_v,
	       !diff_from_left_[2], !diff_from_left_[3]);

    // Fundamental forms
    double E_A, E_B, F_A, F_B, G_A, G_B;
    s1->first_fund_form(par_[0], par_[1],
			!diff_from_left_[0], !diff_from_left_[1],
			E_A, F_A, G_A);
    s2->first_fund_form(par_[2], par_[3],
			!diff_from_left_[2], !diff_from_left_[3],
			E_B, F_B, G_B);
    
    double L_A, L_B, M_A, M_B, N_A, N_B;
    s1->second_fund_form(par_[0], par_[1],
			 !diff_from_left_[0], !diff_from_left_[1],
			 L_A, M_A, N_A);
    s2->second_fund_form(par_[2], par_[3],
			 !diff_from_left_[2], !diff_from_left_[3],
			 L_B, M_B, N_B);

    Point N;
    if (opposite_orientation) {
	N = 0.5 * (N1 - N2);
	L_B *= -1.0;
	M_B *= -1.0;
	N_B *= -1.0;
    }
    else {
	N = 0.5 * (N1 + N2);
    }
    
    double denom_inv = (1.0) / sqrt(E_B * G_B - F_B * F_B);

    double a_11 = determinant(rA_u, rB_v, N) * denom_inv;
    double a_12 = determinant(rA_v, rB_v, N) * denom_inv;
    double a_21 = determinant(rB_u, rA_u, N) * denom_inv;
    double a_22 = determinant(rB_u, rA_v, N) * denom_inv;
    
    double b_11
	= a_11*a_11*L_B
	+ 2.0*a_11*a_21*M_B
	+ a_21*a_21*N_B
	- L_A;
    double b_12
	= a_11*a_12*L_B
	+ (a_11*a_22+a_21*a_12)*M_B
	+ a_21*a_22*N_B
	- M_A;
    double b_22
	= a_12*a_12*L_B
	+ 2.0*a_12*a_22*M_B
	+ a_22*a_22*N_B
	- N_A;

    Point& tan3d_1 = cur_active_2nd_order_properties_->tangent_3d[0];
    Point& tan3d_2 = cur_active_2nd_order_properties_->tangent_3d[1];

    if (fabs(b_11) < EPS && fabs(b_22) < EPS && fabs(b_12) < EPS) {
	// no computable tangent.  We are at a higher order point
	tan3d_1 = Point(double(0), double(0), double(0));
	tan3d_2 = Point(double(0), double(0), double(0));
	return HIGHER_ORDER_POINT;
    }
    
    double discriminant = b_12 * b_12 - b_11 * b_22;
    
    if (discriminant < -EPS) {
	// no real roots - we are at an isolated point
	tan3d_1 = Point(double(0),double(0),double(0));
	tan3d_2 = Point(double(0),double(0),double(0));
	return ISOLATED_POINT;
    } else if (discriminant > EPS) {
	// two real roots - we are at a branch point
	if (fabs(b_11) > fabs(b_22)) {
	    double w_1 = (-b_12 + sqrt(discriminant)) / b_11;
	    double w_2 = (-b_12 - sqrt(discriminant)) / b_11;
	    tan3d_1 = w_1 * rA_u + rA_v; 
	    tan3d_2 = w_2 * rA_u + rA_v;
	} else {
	    double w_1 = (-b_12 + sqrt(discriminant)) / b_22;
	    double w_2 = (-b_12 - sqrt(discriminant)) / b_22;
	    tan3d_1 = w_1 * rA_v + rA_u; 
	    tan3d_2 = w_2 * rA_v + rA_u;
	}
	tan3d_1.normalize();
	tan3d_2.normalize();
	return BRANCH_POINT;
    } else {
	// there is a uniquely defined tangent
	double w;
	if (fabs(b_11) > fabs(b_22)) {
	    w = - b_12/b_11;
	    tan3d_1 = w * rA_u + rA_v;
	} else {
	    w = - b_12/b_22;
	    tan3d_1 = w * rA_v + rA_u;
	}
	tan3d_2 = Point(double(0), double(0), double(0));
	// tangent calculated, now normalizing
	tan3d_1.normalize();
	return TANGENTIAL_POINT;
    }
}

//===========================================================================
bool IntersectionPoint::hasUniqueTangentDirection() const
//===========================================================================
{
    SingularityType st = getSingularityType();
    return (st == ORDINARY_POINT || st == TANGENTIAL_POINT);
}

//===========================================================================
double IntersectionPoint::determinant(const Point& A, const Point& B, const Point& C)
//===========================================================================
{
    return A[0] * (B[1] * C[2] - B[2] * C[1]) -
	   B[0] * (A[1] * C[2] - A[2] * C[1]) +
	   C[0] * (A[1] * B[2] - A[2] * B[1]);
}

//===========================================================================
const double* IntersectionPoint::getPar1() const
//===========================================================================
{
    int nmb_par1 = obj1_->numParams();
    return (nmb_par1 > 0) ? &par_[0] : NULL;
}

//===========================================================================
const double* IntersectionPoint::getPar2() const
//===========================================================================
{
    int nmb_par1 = obj1_->numParams();
    int nmb_par2 = obj2_->numParams();
    return (nmb_par2 > 0) ? &par_[nmb_par1] : NULL;
}

//===========================================================================
Point IntersectionPoint::getPar1Point() const
//===========================================================================
{
    return Point(&par_[0], &par_[0] + obj1_->numParams());
}

//===========================================================================
Point IntersectionPoint::getPar2Point() const
//===========================================================================
{
    int nmb_par1 = obj1_->numParams();
    int nmb_par2 = obj2_->numParams();
    return Point(&par_[0] + nmb_par1, &par_[0] + nmb_par1 + nmb_par2);
}

//===========================================================================
int IntersectionPoint::numParams1() const
//===========================================================================
{
    return obj1_->numParams();
}

//===========================================================================
int IntersectionPoint::numBranches() const
//===========================================================================
{
    SingularityType st = getSingularityType();
    if (!has_defined_tangent(st)) {
	// no defined tangents at all
	return 0; 
    } else if (st != BRANCH_POINT) {
	// one tangent (usual case
	return 1;
    }
    // if we got here, we are at a branch point
    return 2;
}

//===========================================================================
int IntersectionPoint::numParams2() const
//===========================================================================
{
    return obj2_->numParams();
}


//===========================================================================
SingularityType IntersectionPoint::getSingularityType() const
//===========================================================================
{
    // This function has the responsibility for computing the
    // SingularityType and related information.  It will also
    // determine whether the tangent is oriented or not (set the
    // tangent_is_oriented flag).

    if (!singularity_info_is_cached()) {
	calculate_singularity_type(); 
    }
    ASSERT(singularity_info_is_cached());
    return cur_active_2nd_order_properties_->singularity_type;
}


//===========================================================================
bool IntersectionPoint::isNearSingular() const
//===========================================================================
{
    if (!hasUniqueTangentDirection())
	return false;

    const ParamSurfaceInt* sf1 = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    const ParamSurfaceInt* sf2 = dynamic_cast<const ParamSurfaceInt*>(obj2_);
    DEBUG_ERROR_IF(sf1 == 0 || sf2 == 0,
		   "IntersectionPoint::calculate_tangent_at_singular_point() "
		   "is currently only implemented for parametric surface "
		   "intersections.");
    
    Point normal1, normal2;
    sf1->getParamSurface()->normal(normal1, par_[0], par_[1]);
    sf2->getParamSurface()->normal(normal2, par_[2], par_[3]);
    Point tangent = normal1.cross(normal2);
    double angle_tol_multiplier = 1.0e-8;
    double tol = angle_tol_multiplier / epsge_->getAngleTol();
    bool bd1, bd2;
    isBoundaryPoint(bd1, bd2);
    if (bd1 || bd2)
	tol *= 100.0;
    double tangent2 = tangent.length2();
    double eps2 = 100.0 * tol * tol; // Using factor of 100.0
    if (tangent2 < eps2) {
	return true;
    }

    return false;
}


//===========================================================================
IntPtClassification IntersectionPoint::getClassification(const ParamSurface *surf1, 
							 const ParamSurface *surf2,
							 int branch_num) const 
//===========================================================================
{
    DEBUG_ERROR_IF(!surf1 && !surf2, "Zero pointer given to IntersectionPoint::"
	     "getClassification()");
    DEBUG_ERROR_IF(par_.size() != 4, "IntersectionPoint::getClassification() called "
	     "on a IntersectionPoint that is not apparent to lie on a surface-"
	     "surface intersection.");

    if (!has_defined_tangent(getSingularityType())) {
	return DIR_HIGHLY_SINGULAR;
    }

    DEBUG_ERROR_IF(branch_num < 0 || branch_num >= numBranches(), "Invalid branch number "
	     "given to getClassification() routine.");

    RectDomain s1dom, s2dom;
    //double tol = epsge_->getRelParRes();
    double tol = 1.0e-3;
    bool use_second_branch = (branch_num != 0);
    bool oriented = !pointIsSingular();

    // If this is a boundary point with iso-links, we report it as
    // parallel. This "hack" allows us to avoid problems with angular
    // tolerances. @@@jbt
    bool first, second;
    isBoundaryPoint(first, second);
    int isodir[4];
    int niso = hasIsoLinks(isodir);
    if ((first || second) && niso != 0) {
	return DIR_PARALLEL;
    }


    IntPtClassification res_s1 = DIR_UNDEF, res_s2 = DIR_UNDEF;
    if (surf1)
    {
	try {
	    s1dom = dynamic_cast<const RectDomain&>(surf1->parameterDomain());
	} catch (bad_cast) {
	    THROW("Unable to fetch rectangular domains from surfaces.  Are they trimmed?"
		  " Support for trimmed surfaces in getClassification() is not yet "
		  "implemented.");
	}
	res_s1 = 
	    check_against_domain(Point(par_[0], par_[1]), 
				 getPar1Dir(use_second_branch), s1dom, 
				 oriented, epsge_->getRelParRes(), tol);
    }

    if (surf2)
    {
	try {
	    s2dom = dynamic_cast<const RectDomain&>(surf2->parameterDomain());
	} catch (bad_cast) {
	    THROW("Unable to fetch rectangular domains from surfaces.  Are they trimmed?"
		  " Support for trimmed surfaces in getClassification() is not yet "
		  "implemented.");
	}
	res_s2 = 
	    check_against_domain(Point(par_[2], par_[3]), 
				 getPar2Dir(use_second_branch), s2dom, 
				 oriented, epsge_->getRelParRes(), tol);
    }

    if (res_s1 == DIR_UNDEF) {
	return res_s2;
    } else if (res_s2 == DIR_UNDEF) {
	return res_s1;
    } else if (res_s1 == DIR_PARALLEL) {
	return res_s2;
    } else if (res_s2 == DIR_PARALLEL) {
	return res_s1;
    } else if (res_s1 == res_s2) {
	return res_s1;
    } 
    // other cases will be classified as touch
    return DIR_TOUCH;
}


//===========================================================================
int IntersectionPoint::hasIsoLinks(int* const iso_par_dirs) const
//===========================================================================
{
    bool iso_in[4]; // can never be more than 4 parameter directions
    for (int i = 0; i < 4; ++i) {
	iso_in[i] = false;
    }
    int num_parameters = numParams1() + numParams2();

    typedef vector<shared_ptr<IntersectionLink> >::const_iterator iter;
    for (iter it = intersection_links_.begin();
	 it != intersection_links_.end(); ++it) {
	for (int i = 0; i < num_parameters; ++i) {
	    if ((!iso_in[i]) && (*it)->isIsoparametricIn(i)) {
		iso_in[i] = true;
	    }
	}
    }
    
    int num_iso = 0;
    for (int i = 0; i < num_parameters; ++i) {
	if (iso_in[i]) {
	    iso_par_dirs[num_iso++] = i;
	}
    }
    return num_iso;
}


//===========================================================================
int IntersectionPoint::commonIsoLinks(int* const iso_par_dirs) const
//===========================================================================
{
    bool iso_in[4]; // can never be more than 4 parameter directions
    for (int i = 0; i < 4; ++i) {
	iso_in[i] = true;
    }
    int num_parameters = numParams1() + numParams2();

    typedef vector<shared_ptr<IntersectionLink> >::const_iterator iter;
    for (iter it = intersection_links_.begin();
	 it != intersection_links_.end(); ++it) {
	for (int i = 0; i < num_parameters; ++i) {
	    if (!(*it)->isIsoparametricIn(i)) {
		iso_in[i] = false;
	    }
	}
    }
    
    int num_iso = 0;
    for (int i = 0; i < num_parameters; ++i) {
	if (iso_in[i]) {
	    iso_par_dirs[num_iso++] = i;
	}
    }
    return num_iso;
}

//===========================================================================
IntPtClassification
IntersectionPoint::getClassification(const ParamSurface *par1_sf,
				     double C,
				     int branch_num) const
//===========================================================================
{
    DEBUG_ERROR_IF(!par1_sf, "Zero pointer given to IntersectionPoint::"
		   "getClassification()");
    DEBUG_ERROR_IF(par1_sf->dimension() != 1, "Function IntersectionPoint::"
		   "getClassification() expecting 1d-surface!");
    DEBUG_ERROR_IF(par_.size() != 2,
		   "IntersectionPoint::getClassification() called "
		   "on a IntersectionPoint that is not apparent to lie on "
		   "a surface-plane intersection.");

    DEBUG_ERROR_IF(branch_num < 0 || branch_num >= numBranches(),
		   "Invalid branch number "
		   "given to getClassification() routine.");

    // @@sbr Currently not expecting branch pt to differ from 0.
    if (numBranches() > 1) {
	MESSAGE("Not sure if classification handles multiple branches!");
    }

//     SingularityType st = getSingularityType();
//     if (pointIsSingular()) {
// 	return DIR_SINGULAR;
//     } 
    DEBUG_ERROR_IF(!has_defined_tangent(getSingularityType()),
		   "Called getClassification() on an "
		   "IntersectionPoint with no defined tangent.");

    // We must then analyze what kind of intersection point we are
    // dealing with.  I do not think that a degenerate point is
    // possible for a function ...
//     singularity_type_ = ORDINARY_POINT;

    RectDomain sdom;
    try {
	sdom = dynamic_cast<const RectDomain&>(par1_sf->parameterDomain());
    } catch (bad_cast) {
	THROW("Unable to fetch rectangular domains from surface. "
	      "Is it trimmed? "
	      "Support for trimmed surfaces in getClassification() "
	      "is not yet implemented.");
    }

    //double tol = epsge_->getRelParRes();
    double tol = 1.0e-3;

    // The par dir of the intersection should be (S_v, -S_u, 0).
    // (S_3d = (u, v, S(u, v)), etc)
    vector<Point> pts(3);
    par1_sf->point(pts, par_[0], par_[1], 1);
    Point int_dir(pts[2][0], -pts[1][0], 0.0);
    int_dir.normalize();
    bool oriented = !pointIsSingular();
    IntPtClassification res_s = 
	check_against_domain(Point(par_[0], par_[1]), int_dir,
			     sdom, oriented,
			     epsge_->getRelParRes(), tol);

    return res_s;
}


//===========================================================================
bool IntersectionPoint::isConnectedTo(const IntersectionPoint* point) const
//===========================================================================
{
    typedef vector<shared_ptr<IntersectionLink> >::const_iterator iter;
    for (iter it = intersection_links_.begin();
	 it != intersection_links_.end(); ++it) {
	if ((*it)->isAttachedTo(point)) {
	    return true;
	}
    }
    return false;
}


//===========================================================================
shared_ptr<IntersectionLink> IntersectionPoint::
connectTo(IntersectionPoint *const point,
	  LinkType type,
	  shared_ptr<IntersectionLink> model_link,
	  int added_parameter_dir)
//===========================================================================
{
    DEBUG_ERROR_IF(!point,
		   "Null pointer given to IntersectionPoint::connectTo()");
    ASSERT(this != point);

    if (isConnectedTo(point)) {
	// already connected to the point in question
	return getIntersectionLink(point);
    }
    // if we got here, there is no present connection from 'this' to
    // 'point'.
    shared_ptr<IntersectionLink> new_link(new IntersectionLink(this, point));
    new_link->linkType() = type;
    if (model_link.get()) {
	new_link->copyMetaInformation(*model_link);
    }
    if (added_parameter_dir >= 0) {
	// shifting meta-information
	int total_num_dir = numParams1() + numParams2();
	for (int i = added_parameter_dir; i < total_num_dir - 1; ++i) {
	    new_link->setIsoparametricIn(i+1 ,
					 new_link->isIsoparametricIn(i));
	}
    }
    intersection_links_.push_back(new_link);
    point->intersection_links_.push_back(new_link);

    return new_link;
}


//===========================================================================
void IntersectionPoint::disconnectFrom(IntersectionPoint* point)
//===========================================================================
{
    DEBUG_ERROR_IF(!point, "Null pointer given to "
		   "IntersectionPoint::disconnectFrom()");

    vector<shared_ptr<IntersectionLink> >::iterator it
	= intersection_links_.begin();
    while (it != intersection_links_.end()) {
	if ((*it)->isAttachedTo(point)) {
	    it = intersection_links_.erase(it);
	    point->disconnectFrom(this);
	} else {
	    ++it;
	}
    }
    ASSERT(!point->isConnectedTo(this));
}


//===========================================================================
shared_ptr<IntersectionLink> 
IntersectionPoint::getIntersectionLink(const IntersectionPoint* point)
//===========================================================================
{
    typedef vector<shared_ptr<IntersectionLink> >::iterator iter;
    for (iter it = intersection_links_.begin();
	 it != intersection_links_.end(); ++it) {
	if ((*it)->isAttachedTo(point)) {
	    return *it;
	}
    }
    return shared_ptr<IntersectionLink>();
}


//===========================================================================
void 
IntersectionPoint::getNeighbours(vector<IntersectionPoint*>& pnts) const
//===========================================================================
{
    size_t nlinks = intersection_links_.size();
    pnts.resize(nlinks);
    int pos = 0;
    typedef vector<shared_ptr<IntersectionLink> >::const_iterator iter;
    for (iter it = intersection_links_.begin();
	 it != intersection_links_.end(); ++it) {
	pnts[pos] = (*it)->getOtherPoint(this);
	++pos;
    }
}


//===========================================================================
void IntersectionPoint::
getNeighbourLinks(vector<shared_ptr<IntersectionLink> >& links) const 
//===========================================================================
{
    links = intersection_links_;
}


//===========================================================================
double IntersectionPoint::startParam(int pardir)
//===========================================================================
{
    ASSERT(pardir >= 0
	   && pardir < obj1_->numParams() + obj2_->numParams());
    int temp = obj1_->numParams();
    return (pardir < temp)
	? obj1_->startParam(pardir) : obj2_->startParam(pardir - temp);
}


//===========================================================================
double IntersectionPoint::endParam(int pardir)
//===========================================================================
{
    ASSERT(pardir >= 0
	   && pardir < obj1_->numParams() + obj2_->numParams());
    int temp = obj1_->numParams();
    return (pardir < temp)
	? obj1_->endParam(pardir) : obj2_->endParam(pardir - temp);
}


//===========================================================================
double IntersectionPoint::parameterTolerance(int pardir)
//===========================================================================
{
    ASSERT(pardir >= 0
	   && pardir < obj1_->numParams() + obj2_->numParams());
    
    double min_param, max_param;
    if (pardir < obj1_->numParams()) {
	min_param = obj1_->startParam(pardir);
	max_param = obj1_->endParam(pardir);
    } else {
	pardir -= obj1_->numParams();
	min_param = obj2_->startParam(pardir);
	max_param = obj2_->endParam(pardir);
    }
    return epsge_->getRelParRes() * (max_param - min_param);
}

//===========================================================================
int IntersectionPoint::inInfluenceArea(int pardir, double par,
				       bool first_outside) const
//===========================================================================
{
    double low_bound = getInfluenceArea(pardir, 0, first_outside);
    if (par >= low_bound) {
	double high_bound = getInfluenceArea(pardir, 1, first_outside);
	if (par <= high_bound) {
	    // 'par' is in the influence area of this point
	    return exactly_on_point(pardir, par) ? 2 : 1;
	}
    }
    // not inside the influence area
    return 0;
}

//===========================================================================
map<int, SecondOrderProperties> 
IntersectionPoint::
generate_2nd_order_property_map(const std::vector<bool>& disc_params)
//===========================================================================
{
    // Explanation: If any of the parameters participating in the
    // definition of an IntersectionPoint presents a discontinuity at
    // the second derivative on its respective geometrical object,
    // then some methods from differential geometry are no longer
    // valid.  This pertains to the calculation of tangent vectors in
    // singular intersections, where we use second-order information
    // in order to determine the direction of the tangent.  When an
    // IntersectionPoint lies on a discontinuity, these calculations
    // might differ depending on from which 'side' of the parameter
    // you are looking at the IntersectionPoint, ie.  whether you
    // choose to use the double derivative to the 'left' or to the
    // 'right' of the discontinuity.  Since an object has 'p'
    // parameters, each which might possibly be located at a 2nd order
    // discontinuity, the potential number of different ways of
    // computing second order information becomes 2^p.  Each single
    // way can be considered correct when working on a particular part
    // of the parametric domain.  We therefore need to store
    // information related to 2nd derivatives in multiple instances -
    // here we use a map for the task.  This routine will generate all
    // instances that will be necessary for an IntersectionPoint
    // having discontinuous second derivatives for the parameters
    // whose indices are given in the 'disc_params' vector. The key to
    // the map is an unique integer computed based on which parameters
    // we want to differentiate 'on the left'.

    int num_param = (int)disc_params.size();
    
    vector<int> disc_par_ind;
    for (int i = 0; i < num_param; ++i) {
	if (disc_params[i]) {
	    disc_par_ind.push_back(i);
	}
    }

    map<int, SecondOrderProperties> result;

    if (disc_par_ind.size() == 0) {

	result[0] = SecondOrderProperties();

    } else {
    
	vector<bool> derivative_at_left(disc_par_ind.size(), false);
	
	bool all_keys_mapped = false;
	while (!all_keys_mapped) {
	    
	    vector<bool> bool_key_rep(num_param, false);
	    for (int i = 0; i < int(derivative_at_left.size()); ++i) {
		bool_key_rep[disc_par_ind[i]] = derivative_at_left[i];
	    }
	    
	    int key = bool2int(bool_key_rep.begin(), bool_key_rep.end());
	    result[key] = SecondOrderProperties();
	    
	    int flip_ix = 0;
	    bool finished = false;
	    while (!finished) {
		derivative_at_left[flip_ix] = !derivative_at_left[flip_ix];
		if (!derivative_at_left[flip_ix]) {
		    flip_ix++;
		} else {
		    finished = true;
		}
		
		if (flip_ix == int(derivative_at_left.size())) {
		    all_keys_mapped = true;
		    finished = true;
		}
	    }
	}
    }
    return result;
}

//===========================================================================
void IntersectionPoint::
shareInfluenceAreaWith(shared_ptr<IntersectionPoint> other_pt, int pardir)
//===========================================================================
{
    // NB! This function should only be called if 'other_pt' and 'this' point lie
    // in a common influence area in the direction of 'pardir'.
    CachedInterval& fwd_inter = cached_influence_area_forwards_[pardir];
    bool other_after = (getPar(pardir) < other_pt->getPar(pardir));

    fwd_inter.inside = other_after ? 
	other_pt->getPar(pardir) : 
	other_pt->getInfluenceArea(pardir, true, false);
    fwd_inter.outside = other_pt->getInfluenceArea(pardir, true, true);
    fwd_inter.cached = true;

    CachedInterval& backwd_inter = cached_influence_area_backwards_[pardir];
    backwd_inter.inside = other_after ?
	other_pt->getPar(pardir) :
	other_pt->getInfluenceArea(pardir, false, false);
    backwd_inter.outside = other_pt->getInfluenceArea(pardir, false , true);
    backwd_inter.cached = true;
}

//===========================================================================
double IntersectionPoint::getInfluenceArea(int dir, bool forward, 
					   bool first_outside, double aeps) const 
//===========================================================================
{
    ASSERT(dir >= 0 && dir < int(par_.size()));
    bool use_tmp_cach = (aeps > epsge_->getRelParRes() && 
			 fabs(aeps-epsge_->getEpsge())>epsge_->getRelParRes());
    shared_ptr<CachedInterval> tmp_cach(new CachedInterval());
    CachedInterval* cur_cached 
      = use_tmp_cach ? tmp_cach.get() :
	(forward ? &cached_influence_area_forwards_[dir] 
      : &cached_influence_area_backwards_[dir]);

    if (!cur_cached->cached) {
	// these values have not been calculated yet

	// Local tolerance class
	shared_ptr<GeoTol> tmp_tol;
	if (use_tmp_cach)
	    tmp_tol = shared_ptr<GeoTol>(new GeoTol(aeps));
	else
	    tmp_tol = epsge_;

	if (point1_.dist(point2_) >= epsge_->getEpsge())
	{
	    // Point of low accuracy. Let the interval be zero.
	    return par_[dir];
	}
	    

	// If second object is 0-dim we must make a special call.
	int nmb_par_obj2 = obj2_->numParams();
	if (nmb_par_obj2 > 0) {
	  // Wrapping this call in a try-catch block in order to
	  // prevent a crash. @@@jbt
	  try {
	    determineCoincidenceRegion(obj1_, 
				       obj2_, 
				       tmp_tol,
				       &par_[0], 
				       dir, 
				       forward, 
				       cur_cached->inside, 
				       cur_cached->outside);
	  } catch (...) {
	    return par_[dir];
	  }
	} else {
	  const ParamFunctionInt* par_func_int 
	    = dynamic_cast<const ParamFunctionInt*>(obj1_);
	  ASSERT(par_func_int);
	  Point C_pt(1);
	  obj2_->point(C_pt, NULL);
	  determineCoincidenceRegion(par_func_int,
				     C_pt[0],
				     tmp_tol,
				     &par_[0], 
				     dir,
				     forward,
				     cur_cached->inside,
				     cur_cached->outside);
	}
	cur_cached->cached = true;
    } 
    ASSERT(cur_cached->cached);
    return first_outside ? cur_cached->outside : cur_cached->inside;
}

//===========================================================================
void IntersectionPoint::proj_2_pplane(const Point& dir, 
				      const double* par,
				      vector<bool>::const_iterator from_left,
				      const ParamGeomInt* go, 
				      Point& par_dir) 
//===========================================================================
{
    switch (go->numParams()) {
    case 0:
	// isolated point
	par_dir = Point(); // empty, dimensionless point
	break;
    case 1: {
	// curve
	par_dir.resize(1);
	const ParamCurveInt* pci = dynamic_cast<const ParamCurveInt*>(go);
	ASSERT(pci); // this ought to work
	vector<Point> temp(2);
	bool from_right = !*from_left;
	pci->point(temp, par, 1, &from_right);
	par_dir[0] = (temp[1] * dir > 0) ? 1 : -1; 
    } 
	break;
    case 2: {
	// surface
	const ParamSurfaceInt* sf = dynamic_cast<const ParamSurfaceInt*>(go);
	ASSERT(sf); // this ought to work
	par_dir.resize(2);
	bool from_right[2];
	from_right[0] = !(*from_left++);
	from_right[1] = !(*from_left);
	Point udir, vdir;
	sf->derivs(par[0], par[1], udir, vdir, from_right[0], from_right[1]);
	decompose(dir, udir, vdir, par_dir[0], par_dir[1]);
	par_dir.normalize();
    }
	break;
    default:
	ASSERT(false); // should never get here, unless we start dealing with volumes
    }
    par_dir.normalize();
}


//===========================================================================
bool IntersectionPoint::isDegenerate() const
//===========================================================================
{
    // The point is degenerate if one of the surfaces has a zero
    // normal in the point.

    Point normal;
    double eps = epsge_->getEpsge();

    // Check first surface
    const ParamSurfaceInt* sf1 = dynamic_cast<const ParamSurfaceInt*>(obj1_);
    sf1->normal(par_[0], par_[1], normal);
    if (normal.length() < eps)
	return true;

    // Check second surface
    const ParamSurfaceInt* sf2 = dynamic_cast<const ParamSurfaceInt*>(obj2_);
    sf2->normal(par_[2], par_[3], normal);
    if (normal.length() < eps)
	return true;

    // We got here - must be non-degenerate
    return false;
}


//===========================================================================
bool IntersectionPoint::isSamePoint(const IntersectionPoint* point) const
//===========================================================================
{
    // Check if the points coincide in space, as computed from the
    // average of the two surfaces (i.e. the getPoint()-function).
    double epsge = epsge_->getEpsge();
    double ptol = 100.0*epsge_->getRelParRes();
    Point diff = getPoint() - point->getPoint();
    if (diff.length() > epsge)
	return false;

    // VSK, 0710
    // Check distances in the parameter plane. It is too dangerous to remove
    // points that are gometrically equal and distant in the parameter domain,
    // consider closed surfaces, self intersections, etc.
    Point p1 = getPar1Point();
    Point p2 = getPar2Point();
    Point p3 = point->getPar1Point();
    Point p4 = point->getPar2Point();
    if (p1.dist(p3) > ptol || p2.dist(p4) > ptol)
	return false;

    // Convert objects to spline surfaces
    const SplineSurfaceInt* ssi1
	= dynamic_cast<const SplineSurfaceInt*>(obj1_);
    const SplineSurfaceInt* ssi2
	= dynamic_cast<const SplineSurfaceInt*>(obj2_);
    if (ssi1 == 0 || ssi2 == 0) {
	MESSAGE("Only implemented for two surfaces!");
	return false;
    }

    // Now get the midpoints
    double midpar[4];
    for (int i = 0; i < 4; ++i) {
	midpar[i] = 0.5 * (getPar(i) + point->getPar(i));
    }
    Point midpt1, midpt2;
    ssi1->point(midpt1, midpar);
    ssi2->point(midpt2, midpar+2);

    // Test for equality
    diff = getPoint() - midpt1;
    if (diff.length() > epsge)
	return false;
    diff = point->getPoint() - midpt1;
    if (diff.length() > epsge)
	return false;
    diff = getPoint() - midpt2;
    if (diff.length() > epsge)
	return false;
    diff = point->getPoint() - midpt2;
    if (diff.length() > epsge)
	return false;

    return true;



//     // Method: First check if the two points coicide in space within
//     // the tolerance. Then extract a CurveOnSurface from 'this' to
//     // 'point' on the first surface and on the other surface, and
//     // check if they coincide. Works only for two surfaces.

//     // Check if the points coincide in space, as computed from the
//     // average of the two surfaces (i.e. the getPoint()-function).
//     double epsge = epsge_->getEpsge();
//     Point diff = getPoint() - point->getPoint();
//     if (diff.length() > epsge)
// 	return false;

//     // Convert objects to spline surfaces
//     const SplineSurfaceInt* ssi1
// 	= dynamic_cast<const SplineSurfaceInt*>(obj1_);
//     const SplineSurfaceInt* ssi2
// 	= dynamic_cast<const SplineSurfaceInt*>(obj2_);
//     if (ssi1 == 0 || ssi2 == 0) {
// 	MESSAGE("Only implemented for two surfaces!");
// 	return false;
//     }
//     // Can we get rid of this cloning business? Seems to be needed due
//     // to const-problems... @@@jbt
//     shared_ptr<SplineSurface> sf1(ssi1->splineSurface()->clone());
//     shared_ptr<SplineSurface> sf2(ssi2->splineSurface()->clone());

//     // Create spline curves in parameter planes
//     Point from, to;
//     double frompar = 0.0;
//     double topar = 1.0;
//     // Curve in first parameter plane
//     from = Point(par_[0], par_[1]);
//     to = Point(point->getPar(0), point->getPar(1));
//     shared_ptr<SplineCurve> cv1(new SplineCurve(from, frompar, to, topar));
//     // Curve in second parameter plane
//     from = Point(par_[2], par_[3]);
//     to = Point(point->getPar(2), point->getPar(3));
//     shared_ptr<SplineCurve> cv2(new SplineCurve(from, frompar, to, topar));

//     // Create curve-on-surface object
//     bool prefer_par = true;
//     shared_ptr<ParamCurve> cv_on_sf1(new CurveOnSurface(sf1, cv1,
// 							prefer_par));
//     shared_ptr<ParamCurve> cv_on_sf2(new CurveOnSurface(sf2, cv2,
// 							prefer_par));

//     shared_ptr<ParamCurveInt> curve1(new ParamCurveInt(cv_on_sf1));
//     shared_ptr<ParamCurveInt> curve2(new ParamCurveInt(cv_on_sf2));

//     // Finally, check coincidence. If coincidence, coincides = 1,
//     // otherwise coincides = 0.
//     int coincides = checkCoincide(curve1.get(), frompar, topar, epsge_,
// 				  curve2.get(), frompar, topar);
//     if (coincides == 1)
// 	return true;
//     else
// 	return false;

}


//===========================================================================
bool IntersectionPoint::checkIntersectionPoint(IntPtInfo& int_pt_info) const
//===========================================================================
{
    // Only implemented for two surfaces
    if (numParams1() != 2 || numParams2() != 2)
	return false;

    int nneighbours;
    SingularityType singularity_type;
    IntPtClassification direction;
    IntPtLocation location;
    bool is_ok;
    double tol = getTolerance()->getEpsge();

    // Get number of neighbours
    nneighbours =  numNeighbours();
    int_pt_info.nneighbours = nneighbours;

    // Get singularity type
    singularity_type = getSingularityType();
    int_pt_info.singularity_type = singularity_type;

    // Get tangent direction
    const ParamObjectInt* poi1 = getObj1();
    const ParamObjectInt* poi2 = getObj2();
    const ParamSurfaceInt* psi1 = dynamic_cast<const ParamSurfaceInt*>(poi1);
    const ParamSurfaceInt* psi2 = dynamic_cast<const ParamSurfaceInt*>(poi2);
    const ParamSurface* surf1 = psi1->getParamSurface().get();
    const ParamSurface* surf2 = psi2->getParamSurface().get();
    direction = getClassification(surf1, surf2);
    int_pt_info.direction = direction;

    // Get the location of the point (inner, edge, corner). First
    // of all, we know it is in the domain of the pool.
    location = LOC_INSIDE_BOTH;
    // Then check if it is on an edge
    bool edge_on_first, edge_on_second;
    isBoundaryPoint(edge_on_first, edge_on_second);
    if (edge_on_first || edge_on_second) {
	location = LOC_EDGE_ONE;
	if (edge_on_first && edge_on_second) {
	    location = LOC_EDGE_BOTH;
	}
	// Finally, check if it is in a corner.
	bool corner_in_first
	    = psi1->inCorner(getPar1(), tol);
	bool corner_in_second
	    = psi2->inCorner(getPar2(), tol);
	if (corner_in_first || corner_in_second) {
	    location = LOC_CORNER_ONE;
	}
	if (corner_in_first && corner_in_second) {
	    location = LOC_CORNER_BOTH;
	}
    }
    int_pt_info.location = location;

    // Now check if this information is consistent. We do this by
    // checking against a list of cases. Divide into main
    // types of cases depending on the number neighbours.
    is_ok = false;
    if (nneighbours == 0) {
	if (singularity_type == ORDINARY_POINT
	    || singularity_type == TANGENTIAL_POINT) {
	    // Must be in a corner
	    if (location == LOC_CORNER_ONE
		|| location == LOC_CORNER_BOTH) {
		is_ok = true;
	    }
	}

    }
    else if (nneighbours == 1) {
	if (singularity_type == ORDINARY_POINT
	    || singularity_type == TANGENTIAL_POINT
	    || singularity_type == BRANCH_POINT) {
	    // Must be on edge or in corner
	    if (location == LOC_EDGE_ONE
		|| location == LOC_EDGE_BOTH) {
		if (direction == DIR_IN || direction == DIR_OUT
		    || direction == DIR_PARALLEL
		    || direction == DIR_PERPENDICULAR) {
		    is_ok = true;
		}
	    }
	    if (location == LOC_CORNER_ONE
		|| location == LOC_CORNER_BOTH) {
		if (direction == DIR_IN || direction == DIR_OUT
		    || direction == DIR_PARALLEL
		    || direction == DIR_PERPENDICULAR
		    || direction == DIR_TOUCH) {
		    is_ok = true;
		}
	    }
	}
    }
    else if (nneighbours == 2) {
	if (singularity_type == ORDINARY_POINT
	    || singularity_type == TANGENTIAL_POINT) {
	    // Must be on an intersection curve, either
	    // inside or on edge with tangent along edge
	    if (location == LOC_INSIDE_BOTH) {
		is_ok = true;
	    }
	    else if (location == LOC_EDGE_ONE
		     || location == LOC_EDGE_BOTH) {
		if (direction == DIR_PARALLEL) {
		    is_ok = true;
		}
	    }
	}
	if (singularity_type == BRANCH_POINT) {
	    if (location == LOC_EDGE_ONE
		|| location == LOC_EDGE_BOTH
		|| location == LOC_CORNER_ONE
		|| location == LOC_CORNER_BOTH) {
	    is_ok = true;
	    }
	}
    }
    else if (nneighbours == 3) {
	// Could be a point on a degenerate edge
	if ((singularity_type == ORDINARY_POINT
	     || singularity_type == TANGENTIAL_POINT)
	    && (location == LOC_EDGE_ONE
		|| location == LOC_EDGE_BOTH
		|| location == LOC_CORNER_ONE
		|| location == LOC_CORNER_BOTH)
	    && (direction == DIR_IN || direction == DIR_OUT)
	    && isDegenerate()) {
	    is_ok = true;
	}
	// Otherwise it must be a branch point
	if (singularity_type == BRANCH_POINT) {
	    is_ok = true;
	}
    }
    else if (nneighbours == 4) {
	// Must be a nice branch point in the inner
	if (singularity_type == BRANCH_POINT) {
	    if (location == LOC_INSIDE_BOTH) {
		is_ok = true;
	    }
	}
    }
    int_pt_info.is_ok = is_ok;

    return is_ok;
}


//===========================================================================
bool IntersectionPoint::isInDomain(double* frompar, double* topar) const
//===========================================================================
{
    bool inside = true;

    int npar = numParams1() + numParams2();
    for (int i = 0; i < npar; ++i) {
	double intpar =  getPar(i);
	if (intpar < frompar[i] || topar[i] < intpar) {
	    inside = false;
	    break;
	}
    }

    return inside;
}


//===========================================================================


}; // namespace Go


namespace {


//===========================================================================
void dump_tangents(const Point& startpt, const vector<Point>& tangents, string filename)
//===========================================================================
{
    double knots[] = {0, 0, 1, 1};
    ofstream os(filename.c_str());
    for (int i = 0; i < int(tangents.size()); ++i) {
	// defining line
	vector<double> coefs(6);
	for (int j = 0; j < 3; ++j) {
	    coefs[j] = startpt[j];
	    coefs[3 + j] = startpt[j] + tangents[i][j];
	}
	SplineCurve tempcurve(2, 2, knots, &coefs[0], 3);
	tempcurve.writeStandardHeader(os);
	tempcurve.write(os);
    }
    os.close();
}

//===========================================================================
bool collinear_vectors(const Point& p1, const Point& p2, double num_tol, double angle_tol)
//===========================================================================
{
    double norm_p1 = p1.length();
    double norm_p2 = p2.length();
    if (norm_p1 < num_tol && norm_p2 < num_tol) {
	// both vectors are zero.  We can just as well consider them collinear
	return true;
    } else if (norm_p1 < num_tol || norm_p2 < num_tol) {
	// one of the vectors is zero and the other not.  There is no scalar K != 0
	// such that p1 = K * p2.  For this reason, we choose to not consider these
	// vectors collinear
	return false;
    } 
    // if we got here, both vectors are of nonzero length
    return p1.angle(p2) < angle_tol;
}

//===========================================================================
vector<double>
generate_reduced_param_vec(const shared_ptr<IntersectionPoint> ip,
			   int missing_param)
//===========================================================================
{
    int numpar = ip->numParams1() + ip->numParams2();
    vector<double> result;
    for (int ki=0; ki < numpar; ki++) {
	if (ki != missing_param) {
	    result.push_back(ip->getPar(ki));
	}
    }
    return result;
}

//===========================================================================
bool has_defined_tangent(const SingularityType type)
//===========================================================================
{
    // return 'true' for intersection point types that has defined tangents
    return (type == BRANCH_POINT ||
	    type == TANGENTIAL_POINT ||
	    type == ORDINARY_POINT);
}

//===========================================================================
IntPtClassification check_against_domain(const Point& p,
					 const Point& dir,
					 const RectDomain& dom,
					 bool oriented,
					 double res,
					 double angle_tol)
//===========================================================================
{
    // nb!  We suppose here that the direction vector 'dir' is of unit length!
    IntPtClassification result = DIR_UNDEF;
    
    // checking the four boundaries
    double umin = dom.umin();
    double umax = dom.umax();
    double tan_angle = dir[0];
    if (p[0] > umin - res && p[0] < umin + res) {
	// point lies on the lower edge of domain
	if (tan_angle < -angle_tol) {
	    result = DIR_OUT;
	} else if (tan_angle > angle_tol) {
	    result = DIR_IN;
	} else {
	    result = DIR_PARALLEL;
	}
    } else if (p[0] > umax - res && p[0] < umax + res) {
	// point lies on the upper edge of domain
	if (tan_angle < -angle_tol) {
	    result = DIR_IN;
	} else if (tan_angle > angle_tol) {
	    result = DIR_OUT;
	} else {
	    result = DIR_PARALLEL;
	}
    }
//     if (result == DIR_PARALLEL) {
// 	// the outcome of the v-parameter will not matter
// 	return result;
//     }
    double vmin = dom.vmin();
    double vmax = dom.vmax();
    tan_angle = dir[1];
    if (p[1] > vmin - res && p[1] < vmin + res) {
	// point lies on the left edge of domain
	if (tan_angle > angle_tol) {
	    // the vector is pointing 'inwards' seen from this edge
	    result = (result == DIR_OUT) ? DIR_TOUCH : DIR_IN;
	} else if (tan_angle < -angle_tol) {
	    // the vector is pointing 'outwards' seen from this edge
	    result = (result == DIR_IN) ?  DIR_TOUCH : DIR_OUT;
	} else if (result == DIR_UNDEF) {
	    result = DIR_PARALLEL;
	}
    } else if (p[1] > vmax - res && p[1] < vmax + res) {
	// point lies on the right edge of domain
	if (tan_angle < -angle_tol) {
	    // the vector is pointing 'inwards' seen from this edge
	    result = (result == DIR_OUT) ? DIR_TOUCH : DIR_IN;
	} else if (tan_angle > angle_tol) {
	    // the vector is pointing 'outwards' seen from this edge
	    result = (result == DIR_IN) ? DIR_TOUCH : DIR_OUT;
	} else if (result == DIR_UNDEF) {
	    result = DIR_PARALLEL;
	}
    }
    if (!oriented && (result == DIR_IN || result == DIR_OUT)) {
	// we will not keep information about 'in' or 'out', since we haven't really
	// oriented this tangent.  Instead, we will report that it is perpendicular to domain
	result = DIR_PERPENDICULAR;
    }
    return result;
}

}; // end anonymous namespace

// //===========================================================================
// void IntersectionPoint::compute_surface_parametrical_tangents() const
// //===========================================================================
// {
//     shared_ptr<ParamSurfaceInt> sf_1 = 
// 	dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(obj1_);
//     shared_ptr<ParamSurfaceInt> sf_2 = 
// 	dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(obj2_);

//     ASSERT(sf_1 && sf_2);
    
//     Point* tangent_1 = cur_active_2nd_order_properties_->tangent_2d_1;
//     Point* tangent_2 = cur_active_2nd_order_properties_->tangent_2d_2;
    
//     vector<Point> surface_1_info(7);
//     vector<Point> surface_2_info(7);

//     // NB: if this function is to be put to use again, you must make sure that 
//     // the right second-order differential information is used (derive from left/right
//     // in par_[0] ... par_[3].  Currently this information is not explicitly stored
//     // in IntersectionPoint!
//     sf_1->getSurface()->point(surface_1_info, par_[0], par_[1], 2); 
//     sf_2->getSurface()->point(surface_2_info, par_[2], par_[3], 2);
//     sf_1->getSurface()->normal(surface_1_info[6], par_[0], par_[1]);
//     sf_2->getSurface()->normal(surface_2_info[6], par_[2], par_[3]);
//     ASSERT(surface_1_info[0].dimension() == 3); // only works in 3D
//     double s1_info_array[21];
//     double s2_info_array[21];
//     for (int k = 0; k < 7; ++k) {
// 	for (int i = 0; i < 3; ++i) {
// 	    s1_info_array[3 * k + i] = surface_1_info[k][i];
// 	    s2_info_array[3 * k + i] = surface_2_info[k][i];
// 	}
//     }
    
//     double egeo3d[10]; // after s1304, this variable will contain a description of the
//     // intersection in 3D, including position, tangent, radius and 
//     // curvature.  Not used for now, but could come in handy later.
//     double egeo_param_1[7]; // after s1304, this variable will contain a description
//     // of the intersection in the parameter plane of the first
//     // surface.  It contains position, unit tangent, curvature
//     // and radius of curvature.
//     double egeo_param_2[7]; // after s1304, this variable will contain a description
//     // of the intersection in the parameter plane of the second
//     // surface.  It contains position, unit tangent, curvature
//     // and radius of curvature.
//     int jstat;
    
//     s1304(s1_info_array,  // differential information about first surface
// 	  s2_info_array,  // differential information about second surface
// 	  (double*)&par_[0],  // parameter pair for point in first surface
// 	  (double*)&par_[2],  // parameter pair for point in second surface
// 	  egeo3d,
// 	  egeo_param_1,
// 	  egeo_param_2,
// 	  &jstat);
//     DEBUG_ERROR_IF(jstat < 0, "SISL routine s1304 failed in "
// 	     "IntersectionPoint::compute_parametrical_tangents()");
    
//     // Check the orientation of the tangent vectors
//     Point tangent = surface_1_info[6] % surface_2_info[6];
//     Point tangent2(egeo3d[3],egeo3d[4],egeo3d[5]);
//     int sgn = (tangent*tangent2 < 0.0) ? -1 : 1;
    
//     tangent_1[0].resize(2);
//     tangent_2[0].resize(2);
//     tangent_1[0][0] = sgn*egeo_param_1[2]; 
//     tangent_1[0][1] = sgn*egeo_param_1[3]; 
//     tangent_2[0][0] = sgn*egeo_param_2[2]; 
//     tangent_2[0][1] = sgn*egeo_param_2[3]; 
    
//     // the following assertion should hold since we have supposed that the 
//     // IntersectionPoint is not singular.
//     ASSERT(tangent_1[0].length2() > tangent_tol * tangent_tol);
    
//     tangent_1[0].normalize();
//     tangent_1[1] = Point(double(0),double(0));
//     tangent_2[0].normalize();
//     tangent_2[1] = Point(double(0),double(0));

//     cur_active_2nd_order_properties_->tangent_is_oriented = true;
// }

