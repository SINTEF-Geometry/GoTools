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

#ifndef _INTERSECTIONPOOLUTILS_H
#define _INTERSECTIONPOOLUTILS_H


#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionLink.h"
#include "GoTools/intersections/ParamObjectInt.h"
#include "GoTools/intersections/ParamPointInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/Param1FunctionInt.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/intersections/generic_graph_algorithms.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/PointCloud.h"
#include <memory>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>


namespace Go {


//===========================================================================
void debug_write_line(const Point& p1, const Point& p2,
		      const char* fname)
//===========================================================================
{
    std::ofstream os(fname);
    double array[6];
    std::copy(p1.begin(), p1.end(), array);
    std::copy(p2.begin(), p2.end(), array+3);
    LineCloud cloud(array, 1);
    cloud.writeStandardHeader(os);
    cloud.write(os);
    os.close();
}


//===========================================================================
void debug_write_point(const Point& p1, const char* fname)
//===========================================================================
{
    std::ofstream os(fname);
    PointCloud<3> pcloud(p1.begin(), 1);
    pcloud.writeStandardHeader(os);
    pcloud.write(os);
    os.close();
}


//===========================================================================
bool no_parent(const shared_ptr<IntersectionPoint>& p) 
//===========================================================================
{
    return p->parentPoint().get() == 0;
}


//===========================================================================
IntersectionPoint* flip_intersecting_objects(IntersectionPoint* p)
//===========================================================================
{ 
    p->flipObjects(); 
    return p;
}


//===========================================================================
template<class T>
struct raw_pointer_comp 
//===========================================================================
{
    bool operator()(shared_ptr<T> A, shared_ptr<T> B) const
    {
	return A.get() < B.get();
    }
};


//===========================================================================
/// IntersectionPoolUtils. predicate for STL function
class CrossesValue  
//===========================================================================
{
public:
    CrossesValue(int par_dir, double value, double eps) 
	: par_dir_(par_dir), value_(value), eps_(eps) {}
    bool operator()(const shared_ptr<IntersectionLink>& l) const 
    {
	IntersectionPoint *p1, *p2;
	l->getIntersectionPoints(p1, p2);
	double t1 = p1->getPar()[par_dir_];
	double t2 = p2->getPar()[par_dir_];
	bool crosses;
	bool identity = true;
	if (t1 > t2) {
	    std::swap(t1, t2);
	}
	crosses = (((value_ - t1) > eps_ && (t2 - value_) > eps_));

	int numpar = p1->numParams1() + p1->numParams2();
	for (int ki=0; ki<numpar; ki++)
	{
	    if (ki == par_dir_)
		continue;
	    if (fabs(p1->getPar()[ki] - p2->getPar()[ki]) > eps_)
		identity = false;
	}
	return (crosses && !identity);
    }
    typedef const shared_ptr<IntersectionLink> argument_type;
    typedef bool result_type;
private:
    int par_dir_;
    double value_;
    double eps_;
};


//===========================================================================
/// IntersectionPoolUtils. predicate for STL function
class TestInDomain 
//===========================================================================
{
public:
    TestInDomain(const ParamObjectInt* obj1, 
		 const ParamObjectInt* obj2,
		 shared_ptr<IntersectionPoint> ref_point)
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
    
    bool operator()(const shared_ptr<IntersectionLink>& l) const
    {
	IntersectionPoint *p1, *p2;
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
    typedef const shared_ptr<IntersectionLink> argument_type;
    typedef bool result_type;
private:
    std::vector<double> lower_limit_;
    std::vector<double> upper_limit_;
    std::vector<double> epsilon_;
};


//===========================================================================
bool link_is_iso_in(shared_ptr<IntersectionLink> link, int dir)
//===========================================================================
{
    return link->isIsoparametricIn(dir);
}


//===========================================================================
bool link_is_iso(shared_ptr<IntersectionLink> link)
//===========================================================================
{
    return link->isIsoparametric();
}


//===========================================================================
bool link_is_iso_in_other_than(shared_ptr<IntersectionLink> link, int dir)
//===========================================================================
{
    int num_param = link->numParams();
    for (int i = 0; i < num_param; ++i) {
	if (i != dir && link->isIsoparametricIn(i)) {
	    return true;
	}
    }
    return false;
}


//===========================================================================
template<class T1, class T2>
bool compare_first(const std::pair<T1, T2>& x, const std::pair<T1, T2>& y)
//===========================================================================
{
    return x.first < y.first;
}


//===========================================================================
class ClosestPointCalculator
//===========================================================================
{
public:
    ClosestPointCalculator(shared_ptr<ParamObjectInt> obj)
	: point_(dynamic_pointer_cast<ParamPointInt>(obj)),
	  curve_(dynamic_pointer_cast<ParamCurveInt>(obj)),
	  surf_(dynamic_pointer_cast<ParamSurfaceInt>(obj)),
	  func0_(dynamic_pointer_cast<Param0FunctionInt>(obj)),
	  func1_(dynamic_pointer_cast<Param1FunctionInt>(obj)),
	  func2_(dynamic_pointer_cast<Param2FunctionInt>(obj))
    {
	// check that exactly one of free_curve_ and free_surface_ is
	// active
    }
    
    // calculate the closest point between the two objects at the
    // parameter value 'known_par' for the reference object.  The
    // distance at this point is returned.  The value(s) pointed to by
    // 'unknown_par' are also used as seed to the closest point
    // algorithm.
    double compute(const Point& pt, double* unknown_par, double eps)
    {
	// check that exactly one of free_curve_ and free_surface_ is
	// active
// 	int nmb = (point_ != 0) + (curve_ != 0) + (surf_ != 0)
// 	    + (func0_ != 0) + (func1_ != 0) + (func2_ != 0);
	ASSERT((point_.get() != 0) + (curve_.get() != 0) + (surf_.get() != 0)
	       + (func0_.get() != 0) + (func1_.get() != 0) 
               + (func2_.get() != 0));
	Point cl_pt;
	double dist;
	if (point_) { // the other object is a point
	    // very easy.  The closest point must be that point!
	    point_->point(cl_pt, 0);
	    dist = pt.dist(cl_pt);
	}
	if (func0_) { // the other object is a point
	    // very easy.  The closest point must be that point!
	    func0_->point(cl_pt, 0);
	    dist = pt.dist(cl_pt);
	}
	else if (curve_) { // the other object is a curve
	    shared_ptr<ParamCurve> parcurve = curve_->getParamCurve();
	    parcurve->closestPoint(pt, 
				   parcurve->startparam(), 
				   parcurve->endparam(),
				   *unknown_par,
				   cl_pt,
				   dist,
				   unknown_par);
	} else if (func1_) { // the other object is a one-valued function
	    shared_ptr<ParamCurve> parcurve = func1_->getParamCurve();
	    parcurve->closestPoint(pt,
				   parcurve->startparam(),
				   parcurve->endparam(),
				   *unknown_par,
				   cl_pt,
				   dist,
				   unknown_par);
	} else if (surf_) { // the other object is a surface
	    surf_->getParamSurface()->closestPoint(pt, 
						   unknown_par[0], 
						   unknown_par[1], 
						   cl_pt, 
						   dist, 
						   eps, 
						   0, 
						   unknown_par);
	} else { // the other object is a two-valued function
	    func2_->getParamSurface()->closestPoint(pt, 
						    unknown_par[0], 
						    unknown_par[1], 
						    cl_pt, 
						    dist, 
						    eps,
						    0,
						    unknown_par);
	}
	return dist;
    }
    
    
private:
    shared_ptr<ParamPointInt> point_;
    shared_ptr<ParamCurveInt> curve_;
    shared_ptr<ParamSurfaceInt> surf_;
    shared_ptr<Param0FunctionInt> func0_;
    shared_ptr<Param1FunctionInt> func1_;
    shared_ptr<Param2FunctionInt> func2_;
};


//===========================================================================
class LockedDirDistFunc
//===========================================================================
{
public:
    LockedDirDistFunc(const ParamObjectInt* const o1, 
		      const ParamObjectInt* const o2, 
		      int fixed_dir,
		      double fixed_value) :
	o1_(o1), 
	o2_(o2), 
	par2ix_(o1_->numParams()),
	fixed_(fixed_dir), 
	num_par_(o1->numParams() + o2->numParams() - 1),
	pvec1_(3), pvec2_(3)
    {
	ASSERT(num_par_ <= 3);
	ASSERT(fixed_dir >= 0 && fixed_dir < 4);
	int o1_pars = o1->numParams();
	for (int i = 0; i < num_par_; ++i) {
	    int ix = (i < fixed_) ? i : i+1;
	    min_par_[i] = (ix < o1_pars) 
		? o1->startParam(ix) : o2->startParam(ix - o1_pars);
	    max_par_[i] = (ix < o1_pars)
		? o1->endParam(ix) : o2->endParam(ix - o1_pars);
	}
	tmp_[fixed_dir] = fixed_value;
    }

    double operator()(const double* arg) const; 
    void grad(const double* arg, double* grad) const; 
    double minPar(int n) const { return min_par_[n];}
    double maxPar(int n) const { return max_par_[n];}
private:
    void fillParams(const double* arg) const;

    const ParamObjectInt* o1_;
    const ParamObjectInt* o2_;
    const int par2ix_;
    const int fixed_;
    const int num_par_;
    double min_par_[3];
    double max_par_[3];

    mutable double tmp_[4];

    mutable Point p1_, p2_, diff_;
    mutable std::vector<Point> pvec1_, pvec2_;
};


//===========================================================================
class ConnectionFunctor
//===========================================================================
{
public:
    ConnectionFunctor(const std::vector<
		      shared_ptr<IntersectionPoint> >& ipoints)
	: vec_(ipoints) {}
    bool operator()(int a, int b)
    { return vec_[a]->isConnectedTo(vec_[b]);}

private:
    const std::vector<shared_ptr<IntersectionPoint> >& vec_;
};


//===========================================================================


//===========================================================================
void add_reachables_from(const IntersectionPoint* pt,
			 const IntersectionPoint* last_pt,
			 std::set<IntersectionPoint*>& result)
//===========================================================================
{
    std::vector<IntersectionPoint*> n;
    pt->getNeighbours(n);
    for (int i = 0; i < int(n.size()); ++i) {
	if (n[i] != last_pt && result.find(n[i]) == result.end()) {
	    // this point must be added
	    result.insert(n[i]);
	    add_reachables_from(n[i], pt, result);
	}
    }
}


//===========================================================================
void estimate_seed_by_interpolation(const IntersectionPoint* p1,
				    const IntersectionPoint* p2,
				    int fixed_dir,
				    double fixed_value,
				    double* result)
//===========================================================================
{
    result[fixed_dir] = fixed_value;
    const double p1_fixed = p1->getPar(fixed_dir);
    const double p2_fixed = p2->getPar(fixed_dir);
    const double span = p2_fixed - p1_fixed;
    //ASSERT(fabs(span) > 0);
    const int num_param = p1->numParams1() + p1->numParams2();
    
    double ratio = (span == 0) ? 0.5 : (fixed_value - p1_fixed) / span;
    for (int i = 0; i < num_param; ++i) {
	if (i != fixed_dir) {
	    *result++ = (1 - ratio) * p1->getPar(i) + ratio * p2->getPar(i);
	}
    }
}


//===========================================================================
void determine_seed(double* par, 
		    int par_start,
		    int par_end,
		    const Point& pt,
		    const IntersectionPoint * const ip1,
		    const IntersectionPoint * const ip2)
//===========================================================================
{
    const Point& ip1_point = ip1->getPoint();
    const Point& ip2_point = ip2->getPoint();
    double dist1 = pt.dist(ip1_point);
    double dist2 = pt.dist(ip2_point);
    double total_dist = dist1 + dist2;
    double ratio = (total_dist == 0.0) ? 0 : dist1 / total_dist;

    int count = 0;
    for (int i = par_start; i < par_end; ++i) {
	par[count++] = (1 - ratio) * ip1->getPar(i) + ratio * ip2->getPar(i);
    }
}


//===========================================================================
int find_point_in(IntersectionPoint* p,
		  const std::vector<
		  shared_ptr<IntersectionPoint> >& vec)
//===========================================================================
{
    int result;
    for (result = 0; result < int(vec.size()); ++result) {
	if (vec[result].get() == p) {
	    break;
	}
    }
    return result;
}


//===========================================================================
shared_ptr<const CurveOnSurface>
make_curve_on_surface(shared_ptr<const ParamSurface> surf,
		      double u_start,
		      double v_start,
		      double u_end,
		      double v_end,
		      double p_start,
		      double p_end)
//===========================================================================
{
    ASSERT(p_start < p_end);
    std::vector<double> cpoints(4);
    cpoints[0] = u_start;
    cpoints[1] = v_start;
    cpoints[2] = u_end;
    cpoints[3] = v_end;
    std::vector<double> kvec(4);
    kvec[0] = kvec[1] = p_start;
    kvec[2] = kvec[3] = p_end;
    shared_ptr<ParamCurve> pcurve(new SplineCurve(2, 2, &kvec[0],
						  &cpoints[0], 2));
    shared_ptr<ParamSurface> temp_surf
	= const_pointer_cast<ParamSurface>(surf);
    shared_ptr<const CurveOnSurface>
	curve_on_surf(new CurveOnSurface(temp_surf, pcurve, true));
    return curve_on_surf;
}


//===========================================================================
void extract_chains(const std::vector<
		    shared_ptr<IntersectionPoint> > pts,
		    std::vector<std::vector<
		    shared_ptr<IntersectionPoint> > >& chains)
//===========================================================================
{
    chains.clear();
    // find out which IntersectionPoints are connected
    std::vector<std::vector<int> > paths, cycles;
    std::vector<int> singles;
    get_individual_paths((int)pts.size(), ConnectionFunctor(pts),
			 paths, cycles, singles);

    // generating chains of 0 length  (single points)
    for (int i = 0; i < int(singles.size()); ++i) {
	chains.push_back(std::vector<shared_ptr<IntersectionPoint> >
			 (1, pts[singles[i]]));
    }

    // generating chains of nonzero length (linked points)
    for (int p = 0; p < int(paths.size()); ++p) {
	chains.push_back(std::vector<
			 shared_ptr<IntersectionPoint> >());
	std::vector<shared_ptr<IntersectionPoint> >& cur_vec
	    = chains.back();
	for (int pp = 0; pp < int(paths[p].size()); ++pp) {
	    cur_vec.push_back(pts[paths[p][pp]]);
	}
    }
}


//===========================================================================
void LockedDirDistFunc::fillParams(const double* arg) const 
//===========================================================================
{
    for (int i = 0; i < num_par_; ++i) {
	int ix = (i < fixed_) ? i : i+1;
	tmp_[ix] = arg[i];
    }
}


//===========================================================================
double LockedDirDistFunc::operator()(const double* arg) const
//===========================================================================
{
    fillParams(arg);
    o1_->point(p1_, tmp_);
    o2_->point(p2_, tmp_ + par2ix_);
    return p1_.dist2(p2_);
}


//===========================================================================
void LockedDirDistFunc::grad(const double* arg, double* grad) const
//===========================================================================
{
    fillParams(arg);
    int o1_pars = o1_->numParams();
    int o2_pars = o2_->numParams();

    if (o1_pars > 0) {
	o1_->point(pvec1_, tmp_, 1);
    } else {
	o1_->point(pvec1_[0], tmp_);
    }
    if (o2_pars > 0) {
	o2_->point(pvec2_, tmp_ + par2ix_, 1);
    } else {
	o2_->point(pvec2_[0], tmp_ + par2ix_);
    }

    diff_ = pvec2_[0] - pvec1_[0];
    double* comp = grad;
    for (int i = 0; i < o1_pars; ++i) {
	if (i != fixed_) {
	    *comp++ = -2 * diff_ * pvec1_[i+1];
	}
    }
    for (int i = 0; i < o2_pars; ++i) {
	if (i + o1_pars != fixed_) {
	    *comp++ = 2 * diff_ * pvec2_[i+1];
	}
    }
}


// //===========================================================================
// void find_most_distant(const std::vector<shared_ptr<IntersectionPoint> > pts, 
// 		       int pdir1, 
// 		       int pdir2, 
// 		       double& dist1, 
// 		       double& dist2, 
// 		       shared_ptr<IntersectionPoint>& pt1_dir1,
// 		       shared_ptr<IntersectionPoint>& pt2_dir1,
// 		       shared_ptr<IntersectionPoint>& pt1_dir2,
// 		       shared_ptr<IntersectionPoint>& pt2_dir2)
// //===========================================================================
// {
//     // find out which IntersectionPoints are connected
//     std::vector<std::vector<int> > paths, cycles;
//     get_individual_paths(pts.size(), ConnectionFunctor(pts), paths, cycles);
//     // ASSERT(cycles.size() == 0); // this should always be true

//     // find 'longest parametric span' of a chain of IntersectionPoints in the two specified
//     // parametric directions
//     dist1 = dist2 = 0;
//     pt1_dir1 = pt2_dir1 = pt1_dir2 = pt2_dir2 = shared_ptr<IntersectionPoint>(); 
//     shared_ptr<IntersectionPoint> tmp1, tmp2;
//     for (int p = 0; p < int(paths.size()); ++p) {
// 	tmp1 = pts[paths[p].front()];
// 	tmp2 = pts[paths[p].back()];
// 	double d1_tmp = fabs(tmp1->getPar(pdir1) - tmp2->getPar(pdir1));
// 	double d2_tmp = fabs(tmp1->getPar(pdir2) - tmp2->getPar(pdir2));
// 	if (d1_tmp >= dist1) {
// 	    dist1 = d1_tmp;
// 	    pt1_dir1 = tmp1;
// 	    pt2_dir1 = tmp2;
// 	} 
// 	if (d2_tmp >= dist2) {
// 	    dist2 = d2_tmp;
// 	    pt1_dir2 = tmp1;
// 	    pt2_dir2 = tmp2;
// 	}
//     }
// }


// //===========================================================================

} // namespace Go


#endif // _INTERSECTIONPOOLUTILS_H

