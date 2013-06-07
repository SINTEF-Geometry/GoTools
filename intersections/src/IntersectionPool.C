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

#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/IntersectionPoolUtils.h"
#include "GoTools/intersections/IntersectionLink.h"
#include "GoTools/intersections/Intersector.h"
#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/intersections/ParamPointInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/Param1FunctionInt.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include <set>
#include <stdexcept>
#include <functional>
#include <iostream>
#include <iterator>

#include <fstream> // For debugging
#include <stdio.h> // for debugging
#include "GoTools/geometry/LineCloud.h" // for debugging
#include "GoTools/geometry/PointCloud.h" // for debugging


using std::vector;
using std::set;
using std::not1;
using std::cout;
using std::endl;
using std::pair;
using std::swap;
using std::ofstream;
using std::back_inserter;


namespace Go {


//const double IntersectionPool::REL_PAR_RES_ = 1.0e-12;


//===========================================================================
IntersectionPool::IntersectionPool(shared_ptr<ParamObjectInt> obj1, 
				   shared_ptr<ParamObjectInt> obj2,
				   shared_ptr<IntersectionPool> parent, 
				   int missing_dir,
				   double missing_value)
//===========================================================================
    : obj1_(obj1), 
      obj2_(obj2), 
      missing_param_index_(missing_dir),
      missing_param_value_(missing_value),
      prev_pool_(parent)
{
    // we have taken this into a separate function as a workaround for
    // certain problems with the GCC debugger, which sometimes has
    // problems setting breakpoints in constructors
    constructor_implementation(obj1, obj2, parent,
			       missing_dir, missing_value);
}


//===========================================================================
void IntersectionPool::
constructor_implementation(shared_ptr<ParamObjectInt> obj1, 
			   shared_ptr<ParamObjectInt> obj2,
			   shared_ptr<IntersectionPool> parent,
			   int missing_dir,
			   double missing_value)
//===========================================================================
{
    // If the parent pool contains IntersectionPoints, we must see if
    // we need to fetch any of them
    if (prev_pool_.get() == 0) {
	// no parent pool - nothing more to do
	return;
    } 
    if (missing_dir >= 0) {
	parent->splitIntersectionLinks(missing_dir, missing_value);
    }
    const vector<shared_ptr<IntersectionPoint> >& parent_ipoints
	= parent->getIntersectionPoints();
    if (parent_ipoints.size() == 0) {
	// no points to fetch
	return;
    }
    
    // if we got here, the parent pool exists and contains at least
    // one point
    int num_parent_pool_params;
    int num_this_pool_params;
    num_parent_pool_params
	= parent_ipoints[0]->numParams1() 
	+ parent_ipoints[0]->numParams2();
    num_this_pool_params
	= obj1_->numParams()
	+ obj2_->numParams();
    ASSERT(num_parent_pool_params == num_this_pool_params ||
	   num_parent_pool_params == num_this_pool_params + 1);

    fetch_relevant_points(parent_ipoints, missing_dir, missing_value);
}


//===========================================================================
void IntersectionPool::removeBoundaryIntersections(bool single_pt)
//===========================================================================
{
    double frac = single_pt ? 0 : 0.1;
    vector<BoundaryIntersectionData> bound_isects;
    intersectAlongCommonBoundary(frac, bound_isects);
    for (int i = 0; i < int(bound_isects.size()); ++i) {
	// the IntersectionLink itself must always be removed.  The
	// IntersectionPoints themselves will also be removed unless
	// they link to other IntersectionPoints than the ones
	// involved in this boundary intersection.
	vector<shared_ptr<IntersectionPoint> >& cur_chain
	    = bound_isects[i].pts;
	int chain_len = (int)cur_chain.size();
	// removing all links between points in chain
	for (int j = 0; j < chain_len; ++j) {
	    cur_chain[j]->disconnectFrom(cur_chain[(j+1)%chain_len].get());
	}

	// removing all points which are completely without neighbours
	for (int j = 0; j < chain_len; ++j) {
	    if (cur_chain[j]->numNeighbours() == 0) {
		removeIntPoint(cur_chain[j]);
	    }
	}
    }
 }


//===========================================================================
void IntersectionPool::
insertInInfluenceInterval(shared_ptr<IntersectionPoint> pt_in_pool, 
			  double *parvals, 
			  int pardir)
//===========================================================================
{
    // We suppose the 'pt_in_pool' point already exists in this pool
    ASSERT(find(int_points_.begin(), int_points_.end(), pt_in_pool)
	   != int_points_.end());

    // Making new point and inserting it into pool
    const int obj_1_num_par = pt_in_pool->numParams1();
    const int obj_2_num_par = pt_in_pool->numParams2();
    shared_ptr<IntersectionPoint> new_pt
	= addIntersectionPoint(obj1_, obj2_,
			       pt_in_pool->getTolerance(),
			       parvals,
			       parvals + obj_1_num_par);

    // setting link between points
    shared_ptr<IntersectionLink> link
	= new_pt->connectTo(pt_in_pool, LINK_UNDEFINED);
    
    // setting isoparametric if necessary
    int other_par = -1;
    if (obj_1_num_par == 2 && pardir < obj_1_num_par) {
	other_par = (pardir + 1) % 2;
    } else if (obj_2_num_par == 2 && pardir >= obj_1_num_par) {
	other_par = (pardir == obj_1_num_par) ? pardir + 1 : pardir - 1;
    }
    if (other_par >= 0) {
	// The object containing 'pardir' has a second parameter,
	// which we will set to iso.  The following ASSERT is to
	// validate our assumption that this is indeed an isoparameter
	ASSERT(fabs(parvals[other_par] - pt_in_pool->getPar(other_par)) 
	       < pt_in_pool->parameterTolerance(pardir));
	link->setIsoparametricIn(other_par, true);
    }

    // adjusting influence areas
    new_pt->shareInfluenceAreaWith(pt_in_pool, pardir);

}


// //===========================================================================
// bool IntersectionPool::intersectAlongCommonBoundary(double frac,
// 						    int& dir1, 
// 						    double& par1,
// 						    int& dir2,   
// 						    double& par2,
// 						    shared_ptr<IntersectionPoint>& pt1,
// 						    shared_ptr<IntersectionPoint>& pt2)
// //===========================================================================
// {
//     ASSERT(frac > 0 && frac < 1);
//     ASSERT(obj1_->numParams() == obj2_->numParams() == 2);
//     Array<vector<shared_ptr<IntersectionPoint> >, 4> obj_1_on_edges; // entry [i] is the collection of points
//                                                     // in this pool lying on edge i of obj1.
//     Array<vector<shared_ptr<IntersectionPoint> >, 4> obj_2_on_edges; // entry [i] is the collection of points
//                                                     // in this pool lying on edge i of obj2.
//     int i, j;
//     // determine minimal parametric span required for considering an intersection
//     // along a given parameter to be 'non-ignorable'.
//     Array<double, 4> min_param_span;
//     for (i = 0; i < 4; ++i) {
// 	min_param_span[i] = frac * (endParam(i) - startParam(i));
//     }

//     // detetmine on which edges the IntersectionPoints lie
//     locateEdgepoints(obj_1_on_edges, obj_2_on_edges);
    
//     // sort the vectors (necessary for using the set_intersection algorithm later)
//     for (i = 0; i < 4; ++i) {
// 	sort(obj_1_on_edges[i].begin(), obj_1_on_edges[i].end());
// 	sort(obj_2_on_edges[i].begin(), obj_2_on_edges[i].end());
//     }

//     // look for a boundary curve with a paramspan bigger than par_tol
//     vector<shared_ptr<IntersectionPoint> > scratch;
//     for (i = 0; i < 4; ++i, dir1 = i / 2) {
// 	for (j = 0; j < 4; ++j, dir2 = j / 2 + 2) {
// 	    scratch.clear();
// 	    set_intersection(obj_1_on_edges[i].begin(), obj_1_on_edges[i].end(),
// 			     obj_2_on_edges[j].begin(), obj_2_on_edges[j].end(),
// 			     back_inserter(scratch));
// 	    if (scratch.size() > 0) {
// 		shared_ptr<IntersectionPoint> tmp1, tmp2, tmp3, tmp4;
// 		double dist1, dist2;
// 		find_most_distant(scratch, dir1, dir2, dist1, dist2, tmp1, tmp2, tmp3, tmp4);
// 		bool found = false;
// 		if (dist1 >= min_param_span[dir1]) {
// 		    // a non-ignorable intersection has been found
// 		    pt1 = tmp1; pt2 = tmp2;
// 		    found = true;
// 		} else if (dist2 >= min_param_span[dir2]) {
// 		    // a non-ignorable intersection has been found
// 		    pt1 = tmp3; pt2 = tmp4;
// 		    found = true;
// 		}
// 		if (found) {
// 		    par1 = (i % 2 == 0) ? startParam(dir1) : endParam(dir1);
// 		    par2 = (j % 2 == 0) ? startParam(dir2) : endParam(dir2);
// 		    return true;
// 		}
// 	    }
// 	}
//     }
//     // could not find any intersection
//     return false;
// }


// //===========================================================================
// void IntersectionPool::intersectAlongCommonBoundary(double frac,
// 						    vector<BoundaryIntersectionData>& isects)
// //===========================================================================
// {
//     isects.clear();
//     ASSERT(frac >= 0 && frac <= 1);
//     ASSERT(obj1_->numParams() == 2 && obj2_->numParams() == 2);
//     Array<vector<shared_ptr<IntersectionPoint> >, 4> obj_1_on_edges; // entry [i] is the collection of points
//                                                     // in this pool lying on edge i of obj1.
//     Array<vector<shared_ptr<IntersectionPoint> >, 4> obj_2_on_edges; // entry [i] is the collection of points
//                                                     // in this pool lying on edge i of obj2.
//     int i, j;
//     // determine minimal parametric span required for considering an intersection
//     // along a given parameter to be 'non-ignorable'.
//     Array<double, 4> min_pspan;
//     for (i = 0; i < 4; ++i) {
// 	min_pspan[i] = frac * (endParam(i) - startParam(i));
//     }

//     // detetmine on which edges the IntersectionPoints lie
//     locateEdgepoints(obj_1_on_edges, obj_2_on_edges);
    
//     // sort the vectors (necessary for using the set_intersection algorithm later)
//     for (i = 0; i < 4; ++i) {
// 	sort(obj_1_on_edges[i].begin(), obj_1_on_edges[i].end());
// 	sort(obj_2_on_edges[i].begin(), obj_2_on_edges[i].end());
//     }

//     // look for a boundary curve with a paramspan bigger than par_tol
//     vector<shared_ptr<IntersectionPoint> > scratch;
//     vector<vector<shared_ptr<IntersectionPoint> >> chains;
//     BoundaryIntersectionData d;
//     for (i = 0; i < 4; ++i) {
// 	d.dir[0] = i / 2;
// 	for (j = 0; j < 4; ++j) {
// 	    d.dir[1] = j / 2 + 2;
// 	    scratch.clear();
// 	    set_intersection(obj_1_on_edges[i].begin(), obj_1_on_edges[i].end(),
// 			     obj_2_on_edges[j].begin(), obj_2_on_edges[j].end(),
// 			     back_inserter(scratch));
// 	    if (scratch.size() > 0) {
// 		extract_chains(scratch, chains);
// 		d.par[0] = (i % 2 == 0) ? startParam(d.dir[0]) : endParam(d.dir[0]);
// 		d.par[1] = (j % 2 == 0) ? startParam(d.dir[1]) : endParam(d.dir[1]);

// 		int run_ix_1 = (d.dir[0] + 1) % 2;
// 		int run_ix_2 = (d.dir[1] - 2 + 1) % 2 + 2;
// 		for (int c = 0; c < int(chains.size()); ++c) {
// 		    shared_ptr<IntersectionPoint> p1 = chains[c].front();
// 		    shared_ptr<IntersectionPoint> p2 = chains[c].back();
// 		    if (fabs(p1->getPar(run_ix_1) - p2->getPar(run_ix_1)) >= min_pspan[run_ix_1] ||
// 			fabs(p1->getPar(run_ix_2) - p2->getPar(run_ix_2)) >= min_pspan[run_ix_2]) {
// 			// a non-ignorable intersection has been found
// 			d.pts = chains[c];
// 			isects.push_back(d);
// 		    }
// 		}
// 	    }
// 	}
//     }
// }


//===========================================================================
void IntersectionPool::
intersectAlongCommonBoundary(double frac,
			     vector<BoundaryIntersectionData>& isects)
//===========================================================================
{
    isects.clear();
    ASSERT(frac >= 0 && frac <= 1);
    ASSERT(obj1_->numParams() >= 1 && obj2_->numParams() >= 1);
    const int npar1 = obj1_->numParams();
    const int npar2 = obj2_->numParams();
    const int tot_par = npar1 + npar2;
    
    // entry [i] is the collection of points in this pool lying on
    // edge i of obj1.
    Array<vector<shared_ptr<IntersectionPoint> >, 4> obj_1_on_edges; 
    // entry [i] is the collection of points in this pool lying on
    // edge i of obj2.
    Array<vector<shared_ptr<IntersectionPoint> >, 4> obj_2_on_edges; 
    int i, j;
    // Determine minimal parametric span required for considering an
    // intersection along a given parameter to be 'non-ignorable'.
    Array<double, 4> min_pspan;
    for (i = 0; i < tot_par; ++i) {
	min_pspan[i] = frac * (endParam(i) - startParam(i));
    }

    // detetmine on which edges the IntersectionPoints lie
    locateEdgepoints(obj_1_on_edges, obj_2_on_edges);
    
    // sort the vectors (necessary for using the set_intersection
    // algorithm later)
    for (i = 0; i < 4; ++i) {
	sort(obj_1_on_edges[i].begin(), obj_1_on_edges[i].end());
	sort(obj_2_on_edges[i].begin(), obj_2_on_edges[i].end());
    }

    // look for a boundary curve with a paramspan bigger than par_tol
    vector<shared_ptr<IntersectionPoint> > scratch;
    vector<vector<shared_ptr<IntersectionPoint> > > chains;
    BoundaryIntersectionData d;
    const int num_edges_obj1 = (npar1 == 1) ? 1 : 4; // curve -> 1
						     // edge, surface
						     // -> 4 edges
    const int num_edges_obj2 = (npar2 == 1) ? 1 : 4;

    for (i = 0; i < num_edges_obj1; ++i) {
	d.dir[0] = i / 2;
	for (j = 0; j < num_edges_obj2; ++j) {
	    d.dir[1] = j / 2 + npar1;
	    scratch.clear();
	    set_intersection(obj_1_on_edges[i].begin(),
			     obj_1_on_edges[i].end(),
			     obj_2_on_edges[j].begin(),
			     obj_2_on_edges[j].end(),
			     back_inserter(scratch));
	    if (scratch.size() > 0) {
		extract_chains(scratch, chains);
		if (npar1 == 1) {
		    d.par[0] = -1;
		} else {
		    const int fix_par = (d.dir[0] + 1) % 2;
		    d.par[0] = (i % 2 == 0)
			? startParam(fix_par) : endParam(fix_par);
		}
		if (npar2 == 1) {
		    d.par[1] = -1;
		} else {
		    const int fix_par = (d.dir[1] - npar1 + 1) % 2 + npar1;
		    d.par[1] = (j % 2 == 0)
			? startParam(fix_par) : endParam(fix_par);
		}

		for (int c = 0; c < int(chains.size()); ++c) {
		    shared_ptr<IntersectionPoint> p1 = chains[c].front();
		    shared_ptr<IntersectionPoint> p2 = chains[c].back();
		    if (p1.get() == p2.get() && p1->numNeighbours() > 0);   // Chain consisting of one point
		    // is allowed only if the point is isolated
		    else if (fabs(p1->getPar(d.dir[0]) - p2->getPar(d.dir[0]))
			>= min_pspan[d.dir[0]] ||
			fabs(p1->getPar(d.dir[1]) - p2->getPar(d.dir[1]))
			>= min_pspan[d.dir[1]]) {
			// a non-ignorable intersection has been found
			d.pts = chains[c];
			isects.push_back(d);
		    }
		}
	    }
	}
    }
}


//===========================================================================
void IntersectionPool::
locateEdgepoints(Array<vector<shared_ptr<IntersectionPoint> >, 4>& obj1_edges,
		 Array<vector<shared_ptr<IntersectionPoint> >, 4>& obj2_edges)
//===========================================================================
{
//     ASSERT(obj1_->numParams() == 2 && obj2_->numParams() == 2);

    if (int_points_.size() == 0) {
	return;
    }

    const int npar1 = obj1_->numParams();
    const int npar2 = obj2_->numParams();

    if (npar1 == 1) {
	obj1_edges[0] = int_points_; // All IntersectionPoints lay on
				     // an "edge" of a curve...
    } else {
	ASSERT(npar1 == 2);
	for (int n = 0; n < int(int_points_.size()); ++n) {
	    shared_ptr<IntersectionPoint> cur_pt = int_points_[n];
	    for (int dir = 0; dir < 2; ++dir) {
		int other_dir = (dir + 1) % 2;
		if (fabs(cur_pt->getPar(other_dir) - startParam(other_dir)) < 
		    cur_pt->parameterTolerance(other_dir)) {
		    obj1_edges[2 * dir].push_back(cur_pt);
		} else if (fabs(cur_pt->getPar(other_dir)
				- endParam(other_dir)) <
			   cur_pt->parameterTolerance(other_dir)) {
		    obj1_edges[2 * dir + 1].push_back(cur_pt);
		}
	    }
	}
    }
    if (npar2 == 1) {
	obj2_edges[0] = int_points_; // All IntersectionPoints lay on
				     // an "edge" of a curve...
    } else {
	ASSERT(npar2 == 2);
	for (int n = 0; n < int(int_points_.size()); ++n) {
	    shared_ptr<IntersectionPoint> cur_pt = int_points_[n];
	    for (int dir = 0; dir < 2; ++dir) {
		int other_dir = npar1 + ((dir + 1) % 2);
		if (fabs(cur_pt->getPar(other_dir) - startParam(other_dir)) <
		    cur_pt->parameterTolerance(other_dir)) {
		    obj2_edges[2 * dir].push_back(cur_pt);
		} else if (fabs(cur_pt->getPar(other_dir)
				- endParam(other_dir)) <
			   cur_pt->parameterTolerance(other_dir)) {
		    obj2_edges[2 * dir + 1].push_back(cur_pt);
		}
	    }
	}
    }
}


//===========================================================================
double IntersectionPool::startParam(int dir) const
//===========================================================================
{
    if (dir < obj1_->numParams()) {
	return obj1_->startParam(dir);
    } else {
	return obj2_->startParam(dir - obj1_->numParams());
    }
}


//===========================================================================
double IntersectionPool::endParam(int dir) const
//===========================================================================
{
    if (dir < obj1_->numParams()) {
	return obj1_->endParam(dir);
    } else {
	return obj2_->endParam(dir - obj1_->numParams());
    }
}


//===========================================================================
void IntersectionPool::includeCoveredNeighbourPoints()
//===========================================================================
{
    for (int i = 0; i < int(int_points_.size()); ++i) {
	shared_ptr<IntersectionPoint> ip = int_points_[i];
	int parnum = ip->numParams1() + ip->numParams2();
	// see if any of the neighbours of the current
	// IntersectionPoints need to be added to the pool
	vector<IntersectionPoint*> neighs;
	ip->getNeighbours(neighs);
	for (int n = 0; n < int(neighs.size()); ++n) {
	    // searching for this point in the pool
	    int pos = find_point_in(neighs[n], int_points_);
	    if (pos == int(int_points_.size())) {
		// point not found.  Checking if it should be included
		bool include = true;
		for (int p = 0; p < parnum; ++p) {
		    int dir = (p == 0 || p == ip->numParams1()) ? 0 : 1;
		    double tol = ip->parameterTolerance(dir);
		    double p_val = neighs[n]->getPar(p);
		    if (p_val < startParam(p)-tol || p_val > endParam(p)+tol) {
			include = false;
			break;
		    }
		}
		if (include) {
		    // this point should be included in the pool.  Fetch the
		    // shared pointer from parent!
		    ASSERT(prev_pool_.get() != 0);
		    vector<shared_ptr<IntersectionPoint> >& parent_int_points
			= prev_pool_->int_points_;
		    pos = find_point_in(neighs[n], parent_int_points);
		    if (pos < int(parent_int_points.size())) // should
							     // be
							     // guaranteed
							     // by
							     // design
			int_points_.push_back(parent_int_points[pos]);
		    else
			MESSAGE("WARNING! : Illegal intersection point in "
				"includeCoveredNeighbourPoints");

		}
	    }
	}
    }
}


//===========================================================================
void IntersectionPool::synchronizePool()
//===========================================================================
{
    // Check if the current pool contains intersection points not in
    // the parent pool. If so, remove them.

    // Get intersection points from parent pool
    vector<shared_ptr<IntersectionPoint> > orig_int_pts;
    getOrigPoints(orig_int_pts);

    // Test if the points in this pool also exist in the parent pool
    for (int i = 0; i < int(int_points_.size()); ++i) {
	if (int_points_[i]->hasParentPoint())
	    continue;  // Do not touch this point

	vector<shared_ptr<IntersectionPoint> >::iterator it
	    = find(orig_int_pts.begin(), orig_int_pts.end(),
		   int_points_[i]);
	if (it == orig_int_pts.end()) {
	    vector<IntersectionPoint*> neighbours;
	    int_points_[i]->getNeighbours(neighbours);
	    int num_neighbours = (int)neighbours.size();
	    ASSERT(num_neighbours <= 2);

	    if (num_neighbours == 2) {
		int iso_pars[4];
		int nmb_iso = int_points_[i]->commonIsoLinks(iso_pars);
		neighbours[0]->connectTo(neighbours[1], LINK_UNDEFINED);
		shared_ptr<IntersectionLink> lnk
		    = neighbours[0]->getIntersectionLink(neighbours[1]);
		for (int ki=0; ki<nmb_iso; ki++) {
		    lnk->setIsoparametricIn(iso_pars[ki], true);
		}
		int_points_[i]->disconnectFrom(neighbours[0]);
		int_points_[i]->disconnectFrom(neighbours[1]);
	    } else if (num_neighbours == 1) {
		int_points_[i]->disconnectFrom(neighbours[0]);
	    }
	    int_points_.erase(int_points_.begin()+i);
	    --i;
	}
    }

}


//===========================================================================
void IntersectionPool::cleanUpPool(int first_idx, double epsge)
//===========================================================================
{
    // Remove redundant intersection points

    int iso_pars[4];
    for (size_t ki = first_idx; ki < int_points_.size(); ki++) {

	// Check if the point has two neighbours and lies at an
	// isoparametric intersection curve
	if (int_points_[ki]->numNeighbours() == 2)
	{

	    int nmb_iso = int_points_[ki]->commonIsoLinks(iso_pars);
	    if (nmb_iso == 0 && !(obj1_->numParams() == 1 || obj2_->numParams() == 1))
		continue;  // No common iso parametric direction

	    // Do not remove corner points
	    if (obj1_->inCorner(int_points_[ki]->getPar1(), 
				int_points_[ki]->getTolerance()->getRelParRes()))
		continue;
	    if (obj2_->inCorner(int_points_[ki]->getPar2(), 
				int_points_[ki]->getTolerance()->getRelParRes()))
		continue;

	    // Check if both neighbours lie in the current pool
	    vector<IntersectionPoint*> neighbours;
	    int_points_[ki]->getNeighbours(neighbours);
	    int nNeighbours = 2; // There are two neighbours
	    bool inPool[] = { false, false };
	    for (int kj = 0; kj < nNeighbours; kj++) {
		for (size_t kr = 0; kr < int_points_.size(); kr++) {
		    if (neighbours[kj] == int_points_[kr].get()) {
			inPool[kj] = true; // Neighbour in pool
			break;
		    }
		}
		if (!inPool[kj])
		    break; // Not in pool - break and continue
	    }
	    if (!inPool[0] || !inPool[1])
		continue; // Not both neighbours in pool - continue

	    // A redundant point is identified. Remove it in all
	    // intersection pools
	    removeIntPoint(int_points_[ki]);
	    ki--;
	}
	else if (int_points_[ki]->numNeighbours() == 0)
	{
	    IntersectionPoint *curr = int_points_[ki].get();
	    for (int kj=0; kj<first_idx; kj++)
	    {
		if (int_points_[kj]->numNeighbours() == 0)
		    continue;  // Check towards existing links
		double d1 = curr->getPoint().dist(int_points_[kj]->getPoint());
		if (d1 < epsge)
		{
		    removeIntPoint(int_points_[ki]);
		    ki--;
		    break;
		}
	    }
	}
    }
}


//===========================================================================
void IntersectionPool::
removeIntPoint(shared_ptr<IntersectionPoint> int_point)
//===========================================================================
{
    // Remove one intersection point from the pool

    if (prev_pool_.get() && missing_param_index_ < 0) {
	prev_pool_->removeIntPoint(int_point);
    } else {
	// Remove the point from the lists
	vector<IntersectionPoint*> neighbours;
	int_point->getNeighbours(neighbours);
	int num_neighbours = (int)neighbours.size();
	//ASSERT(num_neighbours <= 2);

	LinkType linktype1 = LINK_UNDEFINED;
	LinkType linktype2 = LINK_UNDEFINED;
	if (num_neighbours == 2) {
	    linktype1
		= int_point->getIntersectionLink(neighbours[0])->linkType();
	    linktype2
		= int_point->getIntersectionLink(neighbours[1])->linkType();
	}
	for (int ki=0; ki<num_neighbours; ki++)
	    int_point->disconnectFrom(neighbours[ki]);

	if (num_neighbours == 2) {
	    int iso_pars[4];
	    int nmb_iso = int_point->commonIsoLinks(iso_pars);
	    LinkType newtype = MERGED_UNDEFINED;
	    if ((linktype1 == SIMPLE_CONE
		 || linktype1 == MERGED_SIMPLE_CONE)
		&& (linktype2 == SIMPLE_CONE
		    || linktype2 == MERGED_SIMPLE_CONE))
		newtype = MERGED_SIMPLE_CONE;
	    if ((linktype1 == SIMPLE_MONOTONE
		 || linktype1 == MERGED_SIMPLE_MONOTONE)
		&& (linktype2 == SIMPLE_MONOTONE
		    || linktype2 == MERGED_SIMPLE_MONOTONE))
		newtype = MERGED_SIMPLE_MONOTONE;
	    if (((linktype1 == SIMPLE_CONE
		 || linktype1 == MERGED_SIMPLE_CONE)
		 && (linktype1 == SIMPLE_MONOTONE
		 || linktype1 == MERGED_SIMPLE_MONOTONE))
		|| ((linktype1 == SIMPLE_MONOTONE
		 || linktype1 == MERGED_SIMPLE_MONOTONE)
		    && (linktype2 == SIMPLE_CONE
		    || linktype2 == MERGED_SIMPLE_CONE)))
		newtype = MERGED_SIMPLE_CONE_MONOTONE;
	    if ((linktype1 == COINCIDENCE_SFCV
		 || linktype1 == MERGED_COINCIDENCE_SFCV_SFCV)
		&& (linktype2 == COINCIDENCE_SFCV
		    || linktype2 == MERGED_COINCIDENCE_SFCV_SFCV))
		newtype = MERGED_COINCIDENCE_SFCV_SFCV;
	    neighbours[0]->connectTo(neighbours[1], newtype);
// 	    ofstream out("merged_links.txt", ios_base::app);
// 	    out << newtype << "    (" << linktype1 << " + "
// 		<< linktype2 << ")" << endl;
	    shared_ptr<IntersectionLink> lnk
		= neighbours[0]->getIntersectionLink(neighbours[1]);
	    for (int ki=0; ki<nmb_iso; ki++) {
		lnk->setIsoparametricIn(iso_pars[ki], true);
	    }
	} /*else if (num_neighbours == 1) {
	    int_point->disconnectFrom(neighbours[0]);
	    }*/
	else if (num_neighbours > 2)
	{
	    // Connect to the best point. @@@ VSK. NOT A PERMANENT SOLUTION
	    int kbest=0;
	    double mindist = neighbours[0]->getDist();
	    for (int ki=1; ki<num_neighbours; ki++)
	    {
		double dist = neighbours[ki]->getDist();
		if (dist < mindist)
		{
		    dist = mindist;
		    kbest = ki;
		}
	    }
	    for (int ki=0; ki<num_neighbours; ki++)
	    {
		LinkType newtype = MERGED_UNDEFINED;
		if (ki == kbest)
		    continue;
		neighbours[ki]->connectTo(neighbours[kbest], newtype);
	    }
	}
		
    }

    // Remove the point from the pool
    vector<shared_ptr<IntersectionPoint> >::iterator it
	= find(int_points_.begin(), int_points_.end(), int_point);
    //ASSERT(it != int_points_.end());
    // VSK, 0611. In self intersection context, a corner point of a surface will be
    // found twice for removal. Instead of ensuring that no double occurance of a 
    // point is tried to be removed, I change the ASSERT into a test here.
    if (it != int_points_.end())
	int_points_.erase(it);

    return;
}


//===========================================================================
void IntersectionPool::
removeIntPoint2(shared_ptr<IntersectionPoint> int_point)
//===========================================================================
{
    // Remove one intersection point from the pool

    if (prev_pool_.get() && missing_param_index_ < 0) {
	prev_pool_->removeIntPoint2(int_point);
    } else {
	// Remove the point from the lists
	vector<IntersectionPoint*> neighbours;
	int_point->getNeighbours(neighbours);
	int nneighbours = (int)neighbours.size();
	for (int i = 0; i < nneighbours; ++i) {
	    int_point->disconnectFrom(neighbours[i]);
	}
    }

    // Remove the point from the pool
    vector<shared_ptr<IntersectionPoint> >::iterator it
	= find(int_points_.begin(), int_points_.end(), int_point);
    ASSERT(it != int_points_.end());
    int_points_.erase(it);

    return;
}


//===========================================================================
void IntersectionPool::
removeIntPoints(double* frompar, double* topar, bool only_inner)
//===========================================================================
{
    // We remove all points that are inside the parameter box. The
    // links from these points are disconnected.

    int npar = numParams();
    vector<shared_ptr<IntersectionPoint> > intpts = getIntersectionPoints();
    int nintpts = (int)intpts.size();
    for (int i = 0; i < nintpts; ++i) {
	// Is the point inside?
	bool inside = true;
	for (int j = 0; j < npar; ++j) {
	    double intpar =  intpts[i]->getPar(j);
	    if (intpar < frompar[j] || topar[j] < intpar) {
		inside = false;
		break;
	    }
	}
	if (inside &&
	    !(only_inner && isBoundaryPoint(intpts[i]))) 
	{
	    // It is inside - we remove it
	    removeIntPoint2(intpts[i]);
	}
    }

    return;
}


// //===========================================================================
// void IntersectionPool::splitIntersectionLinks(int fixed_dir, double fixed_val)
// //===========================================================================
// {
//     if (getIntersectionPoints().size() == 0) {
// 	return;
//     }
//     const shared_ptr<IntersectionPoint>& ref_point = int_points_[0];
//     int num_param_obj_1 = ref_point->numParams1();
//     int num_param_obj_2 = ref_point->numParams2();
//     int free_dir = -1;
//     bool other_object_is_surface = false;
//     // finding the other parameter direction of the object containing the fixed dir
//     if (fixed_dir < num_param_obj_1) {
// 	if (num_param_obj_1 == 2) {
// 	    free_dir = (fixed_dir + 1) % 2;
// 	}
// 	other_object_is_surface = num_param_obj_2 == 2;
//     } else {
// 	if (num_param_obj_2 == 2) {
// 	    int temp = fixed_dir - num_param_obj_1;
// 	    free_dir = (temp + 1) % 2 + num_param_obj_1;
// 	}
// 	other_object_is_surface = num_param_obj_1 == 2;
//     }
//     set<shared_ptr<IntersectionLink>, raw_pointer_comp<IntersectionLink> > isoparam_links;
//     vector<shared_ptr<IntersectionLink>> links;
//     vector<shared_ptr<IntersectionLink>>::iterator end;
//     CrossesValue crosses(fixed_dir, fixed_val, ref_point->parameterTolerance(fixed_dir));
//     TestInDomain endpoints_in_object_domains(obj1_.get(), obj2_.get(), ref_point);
    
//     for (int i = 0; i < numIntersectionPoints(); ++i) {
// 	int_points_[i]->getNeighbourLinks(links);
// 	end = remove_if(links.begin(), links.end(), not1(crosses));
// 	end = remove_if(links.begin(), end, not1(endpoints_in_object_domains));
// 	if (free_dir != -1 && other_object_is_surface) {
// 	    // we need there to be an isoparametric direction that is not 'fixed_dir'
// 	    end = remove_if(links.begin(), end, 
// 			    not1(bind2nd(ptr_fun(link_is_iso_in_other_than), fixed_dir)));
// 	}
// 	isoparam_links.insert(links.begin(), end);
//     }

//     // isoparam_links should now contain all isoparametric links that crosses the 
//     // 'fixed_dir' parameter fixed at 'fixed_val', and whose endpoints are both in the
//     // parameter domains of the objects.  They might be isoparametric in any parameter,
//     // not necessarily in 'free_dir'.
//     set<shared_ptr<IntersectionLink>, raw_pointer_comp<IntersectionLink> >::const_iterator it;
//     bool fixed_in_obj1 = (fixed_dir < num_param_obj_1);
//     IntersectionPoint* p1;
//     IntersectionPoint* p2;
//     double par[4];
//     double* known_par =  fixed_in_obj1 ? par : par + num_param_obj_1;
//     int unknown_par_start = fixed_in_obj1 ? num_param_obj_1 : 0;
//     int unknown_par_end = fixed_in_obj1 ? num_param_obj_1 + num_param_obj_2 : num_param_obj_1;
//     double* unknown_par = par + unknown_par_start;

//     par[fixed_dir] = fixed_val;
//     shared_ptr<ParamObjectInt> fixed_obj = fixed_in_obj1 ? obj1_ : obj2_;
//     shared_ptr<ParamObjectInt> other_obj = fixed_in_obj1 ? obj2_ : obj1_;
//     ClosestPointCalculator cpcalc(other_obj);

//     for (it = isoparam_links.begin(); it != isoparam_links.end(); ++it) {
// 	// inserting one point into each of the segments represented by these links
// 	(*it)->getIntersectionPoints(p1, p2);
// 	determine_free_dir_parameter(par, *it, free_dir, fixed_dir, fixed_val);
// 	Point pt;
// 	fixed_obj->point(pt, known_par);

// 	ofstream os("lc.g2");
// 	(*it)->dumpToStream(os);
// 	os.close();
	
// 	determine_seed(unknown_par, unknown_par_start, unknown_par_end, pt, p1, p2);

// 	shared_ptr<GeoTol> tol = p1->getTolerance();
// 	double dist = cpcalc.compute(pt, unknown_par, tol->getRelParRes());

// 	DEBUG_ERROR_IF(dist > tol->getEpsge(), 
// 		 "The assumed IntersectionPoint did not lie at an intersection.");
// 	ASSERT(p1->getObj1() == p2->getObj1());
// 	ASSERT(p1->getObj2() == p2->getObj2());

// 	shared_ptr<IntersectionPoint> new_pt(new IntersectionPoint(p1->getObj1(), 
// 					    p1->getObj2(), 
// 					    tol, 
// 					    par, 
// 					    par + p1->getObj1()->numParams()));
// 	add_point_and_propagate_upwards(new_pt);
// 	new_pt->connectTo(p1, LINK_UNDEFINED, *it);
// 	new_pt->connectTo(p2, LINK_UNDEFINED, *it);
// 	p1->disconnectFrom(p2);
//     }
// }


//===========================================================================
void
IntersectionPool::splitIntersectionLinks(int fixed_dir, double fixed_val)
//===========================================================================
{
    if (getIntersectionPoints().size() == 0) {
	return;
    }
    const shared_ptr<IntersectionPoint>& ref_point = int_points_[0];
    const int num_param = ref_point->numParams1()
	+ ref_point->numParams2();
    // We are using std::set because it automatically takes care of
    // duplication (i.e. the elements are unique)
    set<shared_ptr<IntersectionLink>,
	raw_pointer_comp<IntersectionLink> > links_to_split;
    vector<shared_ptr<IntersectionLink> > links;
    vector<shared_ptr<IntersectionLink> >::iterator end;
    CrossesValue crosses(fixed_dir, fixed_val,
			 ref_point->parameterTolerance(fixed_dir));
    TestInDomain endpoints_in_object_domains(obj1_.get(), obj2_.get(),
					     ref_point);
    
    // We are using std::not1(), which is a "negater"
    for (int i = 0; i < numIntersectionPoints(); ++i) {
	int_points_[i]->getNeighbourLinks(links);
	end = remove_if(links.begin(), links.end(), not1(crosses));
	end = remove_if(links.begin(), end,
			not1(endpoints_in_object_domains));
	links_to_split.insert(links.begin(), end);
    }

    IntersectionPoint *p1, *p2;
    double par[4];
    double seed[4];
    par[fixed_dir] = fixed_val;
    typedef set<shared_ptr<IntersectionLink>,
	raw_pointer_comp<IntersectionLink> >::iterator iter;
    iter it = links_to_split.begin();
    int num_links = (int)links_to_split.size();
    for (int i = 0; i < num_links; ++i) {
	if ((num_param == 4) && !(*it)->isIsoparametric()) {
	    cout << "Trying to split a link that is not isoparametric!"
		 << endl;
	    return;
	}

 	(*it)->getIntersectionPoints(p1, p2);

// 	// determine closest point for this split
	
	double dist;
	iterateToSplitPoint((*it), fixed_dir, fixed_val, par, dist);

	if (dist > p1->getTolerance()->getEpsge())
	{
	    // Try alternative solution
	    estimate_seed_by_interpolation(p1, p2, fixed_dir, fixed_val, seed);
	    LockedDirDistFunc fun(obj1_.get(), obj2_.get(), fixed_dir, fixed_val);
	    FunctionMinimizer<LockedDirDistFunc> fmin(num_param - 1, fun, seed);

	    minimise_conjugated_gradient(fmin);
	    dist = fmin.fval();
	    if (dist > p1->getTolerance()->getEpsge()) {
		DEBUG_ERROR_IF(dist > p1->getTolerance()->getEpsge(), 
			       "The assumed IntersectionPoint did not lie at "
			       "an intersection.");
	    }

	    // copying calculated parameters
	    for (int i = 0; i < num_param; ++i) {
		if (i < fixed_dir) {
		    par[i] = fmin.getPar(i);
		} else if (i > fixed_dir) {
		    par[i] = fmin.getPar(i-1);
		} else {
		    par[i] = fixed_val; // i == fixed_dir
		}
	    }
	}

	    

	    // adding new intersection point
	int offset = p1->getObj1()->numParams();
	shared_ptr<IntersectionPoint>
	    new_pt(new IntersectionPoint(p1->getObj1(), 
					 p1->getObj2(), 
					 p1->getTolerance(),
					 par, 
					 par + offset));
	add_point_and_propagate_upwards(new_pt);
	p1->disconnectFrom(p2);
	new_pt->connectTo(p1, SPLIT_LINK, *it);
	new_pt->connectTo(p2, SPLIT_LINK, *it);
	++it;
    }
}


//===========================================================================
void IntersectionPool::
iterateToSplitPoint(shared_ptr<IntersectionLink> link,
		    int fixed_dir,
		    double fixed_val,
		    double par[], double& dist)
//===========================================================================
{
    // Iterate to a split point of an intersection link with respect to
    // a given value in a given parameter direction. Keep iso-parametric
    // directions in the link.
    IntersectionPoint *p1, *p2;  // Intersection points at each end of the link
    link->getIntersectionPoints(p1, p2);
    int np1 = p1->numParams1();
    int np2 = p1->numParams2();
    int num_param = np1 + np2;

    // Set initial parameter value
    /*for (int ki=0; ki<num_param; ki++)
	par[ki] = 0.5*(p1->getPar(ki) + p2->getPar(ki));
	par[fixed_dir] = fixed_val;*/
    double seed[4];
    int ki;
    estimate_seed_by_interpolation(p1, p2, fixed_dir, fixed_val, seed);
    for (ki=0; ki<fixed_dir; ki++)
	par[ki] = seed[ki];
    par[fixed_dir] = fixed_val;
    for (ki=fixed_dir+1; ki<num_param; ki++)
	par[ki] = seed[ki-1];

    double tol = p1->getTolerance()->getEpsge();

    Point pos1, pos2;
    shared_ptr<ParamCurve> crv1, crv2;
    shared_ptr<ParamSurface> surf1, surf2;
    double tmin1, tmin2, tmax1, tmax2;
    RectDomain domain1, domain2;
    int cv1_idx=-1, cv2_idx=-1, sf1_idx=-1, sf2_idx=-1;
    int fix_idx = (fixed_dir < np1) ? fixed_dir : fixed_dir-np1;

    // Fetch object involved in closest point iteration
    if ((fixed_dir < np1 && (np1 == 1 || link->isIsoparametricIn(1-fix_idx))) ||
	(np1 == 0))
	obj1_->point(pos1, par);
    else if (fixed_dir < np1 && np1 == 2)
    {
	crv1 = obj1_->getParamSurfaceInt()->
	    getConstantParameterCurve(fix_idx, fixed_val);
	cv1_idx = 1-fix_idx;
	tmin1 = crv1->startparam();
	tmax1 = crv1->endparam();
    }
    else if (fixed_dir >= np1 && np1 == 1)
    {
	//crv1 = obj1_->getParamCurveInt()->getParentParamCurve(tmin1, tmax1);
	crv1 = obj1_->getParamCurveInt()->getParamCurve();
	cv1_idx = 0;
	tmin1 = crv1->startparam();
	tmax1 = crv1->endparam();
    }
    else if (fixed_dir >= np1 && np1 == 2 && (link->isIsoparametricIn(0) ||
					      link->isIsoparametricIn(1)))
    {
	int idx = (link->isIsoparametricIn(0)) ? 0 : 1;
	crv1 = obj1_->getParamSurfaceInt()->
			getConstantParameterCurve(idx, par[idx]);
	cv1_idx = 1-idx;
 	tmin1 = crv1->startparam();
	tmax1 = crv1->endparam();
   }
    else if (np1 == 2)
    {
	// surf1 = obj1_->getParamSurfaceInt()->getParentParamSurface(domain1);
	surf1 = obj1_->getParamSurfaceInt()->getParamSurface();
	domain1 = obj1_->getParamSurfaceInt()->getDomain();
	sf1_idx = 0;
    }
    else 
	return;
	     
    if ((fixed_dir >= np1 && (np2 == 1 || link->isIsoparametricIn(1-fix_idx))) ||
	(np2 == 0))
	obj2_->point(pos2, par+np1);
    else if (fixed_dir >= np1 && np2 == 2)
    {
	crv2 = obj2_->getParamSurfaceInt()->
	    getConstantParameterCurve(fix_idx, fixed_val);
	cv2_idx = np1 + 1 - fix_idx;
	tmin2 = crv2->startparam();
	tmax2 = crv2->endparam();
    }
    else if (fixed_dir < np1 && np2 == 1)
    {
	//crv2 = obj2_->getParamCurveInt()->getParentParamCurve(tmin2, tmax2);
	crv2 = obj2_->getParamCurveInt()->getParamCurve();
	cv2_idx = np1;
	tmin2 = crv2->startparam();
	tmax2 = crv2->endparam();
    }
    else if (fixed_dir < np1 && np2 == 2 && (link->isIsoparametricIn(np1) ||
					      link->isIsoparametricIn(np1+1)))
    {
	int idx = (link->isIsoparametricIn(np1)) ? 0 : 1;
	crv2 = obj2_->getParamSurfaceInt()->
			getConstantParameterCurve(idx, par[idx]);
	cv2_idx = np1 + 1 - idx;
	tmin2 = crv2->startparam();
	tmax2 = crv2->endparam();
    }
    else if (np2 == 2)
    {
	// surf2 = obj2_->getParamSurfaceInt()->getParentParamSurface(domain2);
	surf2 = obj2_->getParamSurfaceInt()->getParamSurface();
	domain2 = obj2_->getParamSurfaceInt()->getDomain();
	sf2_idx = np1;
    }
    else 
	return;
	     
	
    // Iterate
    if (surf1.get() && surf2.get())
	return;  // Should not occur
    else if (surf1.get() && crv2.get())
    {
	// Curve - surface iteration
	Point pt_cv, pt_sf;
	ClosestPoint::closestPtCurveSurf(crv2.get(), surf1.get(), tol, tmin2, 
		      tmin1, &domain1, par[cv2_idx], par+sf1_idx, 
		      par[cv2_idx], par+sf1_idx, dist, pt_cv, pt_sf);
    }
    else if (crv1.get() && surf2.get())
    {
 	// Curve - surface iteration
	Point pt_cv, pt_sf;
	ClosestPoint::closestPtCurveSurf(crv1.get(), surf2.get(), tol, tmin1, 
		      tmin2, &domain2, par[cv1_idx], par+sf2_idx, 
		      par[cv1_idx], par+sf2_idx, dist, pt_cv, pt_sf);
   }
    else if (crv1.get() && crv2.get())
    {
	Point pt1, pt2;
	ClosestPoint::closestPtCurves(crv1.get(), crv2.get(), tmin1, tmax1, tmin2, tmax2,
			par[cv1_idx], par[cv2_idx], par[cv1_idx], par[cv2_idx],
			dist, pt1, pt2);
    }
    else if (surf1.get())
    {
	Point clo_pt;
	surf1->closestPoint(pos2, par[sf1_idx], par[sf1_idx+1], clo_pt, dist, tol,
			    &domain1, par+sf1_idx);
    }
    else if (surf2.get())
    {
	Point clo_pt;
	surf2->closestPoint(pos1, par[sf2_idx], par[sf2_idx+1], clo_pt, dist, tol,
			    &domain2, par+sf2_idx);
    }
    else if (crv1.get())
    {
	// Iterate to point - curve intersection
	Point clo_pt;
	crv1->closestPoint(pos2, tmin1, tmax1,
			   par[cv1_idx], clo_pt, dist, par+cv1_idx);
    }
    else if (crv2.get())
    {
	// Iterate to point - curve intersection
	Point clo_pt;
	crv2->closestPoint(pos1, tmin2, tmax2,
			   par[cv2_idx], clo_pt, dist, par+cv2_idx);
    }
    else
	return;

}

//===========================================================================
void IntersectionPool::
fetch_relevant_points(const vector<shared_ptr<IntersectionPoint> >& ipoints,
		      int missing_dir,
		      double missing_value)
//===========================================================================
{
    // This function adds relevant points among 'ipoints' to this
    // pool.  We suppose that the number of parameters of the points
    // is equal to this pool's.

    if (ipoints.size() == 0)
	return;
    ASSERT(prev_pool_.get());                               
    bool selfintersect
	= (prev_pool_->obj1_ == prev_pool_->obj2_ && missing_dir == -1);
    int num_param_1 = ipoints[0]->numParams1();
    int num_param_2 = ipoints[0]->numParams2();
    if (num_param_1 != 2 || num_param_2 != 2)
	selfintersect = false;

    // Get parameter ranges for objects in 'this' pool
    vector<double> range_start;
    vector<double> range_end;
    for (int i = 0; i < obj1_->numParams(); ++i) {
	range_start.push_back(obj1_->startParam(i));
	range_end.push_back(obj1_->endParam(i));
    }
    for (int i = 0; i < obj2_->numParams(); ++i) {
	range_start.push_back(obj2_->startParam(i));
	range_end.push_back(obj2_->endParam(i));
    }

    // Add missing value if there is a missing direction
    if (missing_dir >= 0) {
	ASSERT(missing_dir <= int(range_start.size()));
	range_start.insert(range_start.begin() + missing_dir,
			   missing_value);
	range_end.insert(range_end.begin() + missing_dir,
			 missing_value);
    }

    int numpar = num_param_1 + num_param_2;
    vector<double> pt_par(numpar);
    vector<shared_ptr<IntersectionPoint> > twin_pts;
    for (int i = 0; i < int(ipoints.size()); ++i) {

	// registering parameters from first object
	shared_ptr<IntersectionPoint> cur_point = ipoints[i];
	for (int j = 0; j < num_param_1; ++j) {
	    pt_par[j] = cur_point->getPar1()[j];
	}
	// registering parameters from second object
	for (int j = 0; j < num_param_2; ++j) {
	    pt_par[j + num_param_1] = cur_point->getPar2()[j];
	}

	// checking against parameter range and adding point
	double tol = cur_point->getTolerance()->getRelParRes();
	bool param_inside_range = true;
	for (int j = 0; j < numpar; ++j) {
	    if ((pt_par[j] < range_start[j] - tol)
		|| (pt_par[j] > range_end[j] + tol)) {
		param_inside_range = false;
		break;
	    }
	}
	if (param_inside_range) {
	    // this point should be added (it's on an intersection
	    // between 2 objects)
	    if (missing_dir == -1) {
		int_points_.push_back(cur_point);
	    } else {
		shared_ptr<IntersectionPoint>
		    temp(new IntersectionPoint(obj1_.get(), obj2_.get(), 
					       cur_point, missing_dir));
		int_points_.push_back(temp);
	    } 
	} else if (selfintersect) {
	    // test on range where two objects are switched
	    std::rotate(pt_par.begin(), pt_par.begin() + num_param_1,
			pt_par.end());
	    param_inside_range = true;
	    for (int j = 0; j < numpar; ++j) {
		if ((pt_par[j] < range_start[j] - tol)
		    || (pt_par[j] > range_end[j] + tol)) {
		    param_inside_range = false;
		    break;
		}
	    }
	    if (param_inside_range) {
		// this point should be added (it's on a self-intersection)
// 		ASSERT(missing_dir == -1); // doesn't make sense to
// 					   // have a missing dir
// 					   // here...
		shared_ptr<IntersectionPoint>
		    temp(new IntersectionPoint(obj1_.get(), obj2_.get(),
					       cur_point->getTolerance(),
					       cur_point->getPar2(),
					       cur_point->getPar1()));
		temp->setParentPoint(cur_point); // not really parent,
						 // but "twin"...
		twin_pts.push_back(temp);
	    }
	}
    }
    // Pass through the twin points than can occur in a self intersection
    // setting and make sure that they are unique before they are inserted
    // into the regular pool of intersection points
    size_t num_ints = int_points_.size();
    size_t ki, kj;
    int kr;
    for (ki=0; ki<twin_pts.size(); ki++)
    {
	shared_ptr<IntersectionPoint> cur_point = twin_pts[ki];
	double tol = cur_point->getTolerance()->getRelParRes();
	for (kj=0; kj<num_ints; kj++)
	{
	    for (kr=0; kr<numpar; kr++)
		if (fabs(cur_point->getPar(kr)-int_points_[kj]->getPar(kr)) >= tol)
		    break;   // Not identical point
	    if (kr == numpar)
		break;   // An identical point is found
	}
	if (kj == num_ints)
	{
	    // No identical point is found. Use the twin
	    int_points_.push_back(twin_pts[ki]);
	}
    }

}


//===========================================================================
bool IntersectionPool::
atDifferentBoundary(shared_ptr<IntersectionPoint> pt1,
		    shared_ptr<IntersectionPoint> pt2)
//===========================================================================
{
    int total_num_par = pt1->numParams1() + pt1->numParams2();
    ASSERT(pt2->numParams1() + pt2->numParams2() == total_num_par);

    const vector<double>& pt1_par = pt1->getPar();
    const vector<double>& pt2_par = pt2->getPar();
    double rel_par_res = pt1->getTolerance()->getRelParRes();

    bool at_bd1[8], at_bd2[8];
    int ki;
    for (ki=0; ki<total_num_par; ++ki)
    {
	at_bd1[2*ki] = at_bd1[2*ki+1] = at_bd2[2*ki] = at_bd2[2*ki+1] = false;
	double min_par_val;
	double max_par_val;
	if (ki < obj1_->numParams()) {
	    min_par_val = obj1_->startParam(ki);
	    max_par_val = obj1_->endParam(ki);
	} else {
	    min_par_val = obj2_->startParam(ki - obj1_->numParams());
	    max_par_val = obj2_->endParam(ki - obj1_->numParams());
	}

	if (fabs(min_par_val - pt1_par[ki]) < rel_par_res)
	    at_bd1[2*ki] = true;
	    
	if (fabs(max_par_val - pt1_par[ki]) < rel_par_res)
	    at_bd1[2*ki+1] = true;
	    
	if (fabs(min_par_val - pt2_par[ki]) < rel_par_res)
	    at_bd2[2*ki] = true;
	    
	if (fabs(max_par_val - pt2_par[ki]) < rel_par_res)
	    at_bd2[2*ki+1] = true;
    }

    for (ki=0; ki<2*total_num_par; ++ki)
	if ((at_bd1[ki] && !at_bd2[ki]) ||
	    (!at_bd1[ki] && at_bd2[ki]))
	    return true;

    return false;
}


//===========================================================================
bool IntersectionPool::
atSameBoundary(shared_ptr<IntersectionPoint> pt1,
	       shared_ptr<IntersectionPoint> pt2)
//===========================================================================
{
    int total_num_par = pt1->numParams1() + pt1->numParams2();
    ASSERT(pt2->numParams1() + pt2->numParams2() == total_num_par);

    const vector<double>& pt1_par = pt1->getPar();
    const vector<double>& pt2_par = pt2->getPar();
    double rel_par_res = pt1->getTolerance()->getRelParRes();
    
    for (int p = 0; p < total_num_par; ++p) {
	double min_par_val;
	double max_par_val;
	if (p < obj1_->numParams()) {
	    min_par_val = obj1_->startParam(p);
	    max_par_val = obj1_->endParam(p);
	} else {
	    min_par_val = obj2_->startParam(p - obj1_->numParams());
	    max_par_val = obj2_->endParam(p - obj1_->numParams());
	}
	if (fabs(min_par_val - pt1_par[p]) < rel_par_res &&
	    fabs(min_par_val - pt2_par[p]) < rel_par_res) {
	    // both points on lower edge of object
	    return true;
	} else if (fabs(max_par_val - pt1_par[p]) < rel_par_res &&
		   fabs(max_par_val - pt2_par[p]) < rel_par_res) {
	    // both points on upper edge of object
	    return true;
	} 
    }
    // if we got here, no common edge was found
    return false;
}


//===========================================================================
bool IntersectionPool::
isBoundaryIntersection(shared_ptr<IntersectionPoint> pt1,
		       shared_ptr<IntersectionPoint> pt2)
//===========================================================================
{
    // Check if the points lie at the same boundare
    if (!atSameBoundary(pt1, pt2))
	return false;

    // Check if the points are linked with constant parameter linkage
    // For the time being, we assume that there is not intermediate points
    shared_ptr<IntersectionLink> lnk = pt1->getIntersectionLink(pt2.get());
    if (lnk.get() != 0 && lnk->isIsoparametric())
	return true;

    return false;
}


//===========================================================================
bool IntersectionPool::isBoundaryPoint(IntersectionPoint* pt1)
//===========================================================================
{
    double tol = pt1->getTolerance()->getRelParRes();
    bool first = obj1_->boundaryPoint(pt1->getPar1(), tol);
    bool second  = obj2_->boundaryPoint(pt1->getPar2(), tol);

    return (first || second);
}

//===========================================================================
bool IntersectionPool::
fetchLoops(vector<vector<shared_ptr<IntersectionPoint> > >& loop_ints)
//===========================================================================
{
    // This function assumes that the topology of the loops in this
    // pool is not too 'warped'.  Notably, it might fail if two loops
    // share some edges, since in that case there will be ambiguous
    // choices for what constitutes the loops.  This is a fact from
    // graph theory, and if we want to work around it, we might have
    // to consider additional information of a non-topological nature,
    // like looking at the geometric/parametric position of the
    // points.

    int num_points_in_pool = (int)int_points_.size();
    if (num_points_in_pool < 3) {
	return false; // must be at least three points to make a loop
    }
    vector<vector<int> > loops;
    get_fundamental_cycle_set(num_points_in_pool,
			      ConnectionFunctor(int_points_), loops);

    loop_ints.resize(loops.size());
    for (int i = 0; i < int(loops.size()); ++i) {
	loop_ints[i].resize(loops[i].size());
	for (int j = 0; j < int(loops[i].size()); ++j) {
	    loop_ints[i][j] = int_points_[loops[i][j]];
	}
    }
    return loop_ints.size() > 0;
}


//===========================================================================
bool IntersectionPool::
isPAC(vector<shared_ptr<IntersectionPoint> >& int_loop)
//===========================================================================
{
    bool result = true;
    int loop_size = (int)int_loop.size();
    for (int i = 0; i < loop_size && result; ++i) {
	shared_ptr<IntersectionLink> link 
	    = int_loop[i]->getIntersectionLink(int_loop[(i+1)
							% loop_size].get());
	DEBUG_ERROR_IF(link.get() == 0,
		       "Error in argument given to "
		       "IntersectionPool::setCoincidence().");
	result = link->delimitsPAC();
    }
    return result;
}


//===========================================================================
void IntersectionPool::
setCoincidence(vector<shared_ptr<IntersectionPoint> >& int_loop)
//===========================================================================
{
    int loop_size = (int)int_loop.size();
    for (int i = 0; i < loop_size; ++i) {
	shared_ptr<IntersectionLink> link
	    = int_loop[i]->getIntersectionLink(int_loop[(i+1)
							% loop_size].get());
	DEBUG_ERROR_IF(link.get() == 0,
		       "Error in argument given to "
		       "IntersectionPool::setCoincidence().");
	link->setPAC(true);
    }
}


//===========================================================================
bool IntersectionPool::isConnectedInside(shared_ptr<IntersectionPoint> pt1,
					 shared_ptr<IntersectionPoint> pt2)
//===========================================================================
{
    // determining indices of pt1 and pt2 in pool's point table
    int pt1_ix = -1;
    int pt2_ix = -1;
    
    for (int i = 0; i < int(int_points_.size()); ++i) {
	if (int_points_[i] == pt1) {
	    pt1_ix = i;
	}
	if (int_points_[i] == pt2) {
	    pt2_ix = i;
	}
    }
    DEBUG_ERROR_IF(pt1_ix == -1 || pt2_ix == -1, 
		   "Points not in IntersectionPool given to "
		   "IntersectionPool::isConnectedInside()");

    return is_path_connected(pt1_ix, pt2_ix, (int)int_points_.size(),
			     ConnectionFunctor(int_points_));
}


//===========================================================================
void IntersectionPool::addParamPoints(int nmb_int_pts, 
				      double* pointpar1, 
				      double* pointpar2,
				      shared_ptr<GeoTol> epsge)
//===========================================================================
{
    int num_param_1 = obj1_->numParams();
    int num_param_2 = obj2_->numParams();

    for (int i = 0; i < nmb_int_pts; ++i) {
	shared_ptr<IntersectionPoint> 
	    temp(new IntersectionPoint(obj1_.get(), obj2_.get(), epsge,
				       pointpar1, pointpar2));
	int_points_.push_back(temp);
	pointpar1 += num_param_1;
	pointpar2 += num_param_2;
    }
}


//===========================================================================
void IntersectionPool::
getBranchPoints(vector<shared_ptr<IntersectionPoint> > & pts)
//===========================================================================
{
    pts.clear();
    vector<shared_ptr<IntersectionPoint> >::const_iterator it;
    for (it = int_points_.begin(); it != int_points_.end(); ++it) {
	// assure that we are on the right differentiating side
	if ((*it)->containsG2Discontinuity()) {
	    set_correct_differentiating_domain(*it);
	}

	// collecting all branch points
	if ((*it)->getSingularityType() == BRANCH_POINT) {
	    pts.push_back(*it);
	}
    }
}


//===========================================================================
int IntersectionPool::numSingularIntersectionPoints()
//===========================================================================
{
    int nmbsing = 0;
    vector<shared_ptr<IntersectionPoint> >::const_iterator it;
    for (it = int_points_.begin(); it != int_points_.end(); ++it) 
	if ((*it)->getSingularityType() != ORDINARY_POINT)
	    nmbsing++;
    return nmbsing;
}


//===========================================================================
void IntersectionPool::
getBoundaryIntersections(vector<shared_ptr<IntersectionPoint> > & bd_ints)
  //===========================================================================
{
//    MESSAGE("Warning: IntersectionPool::getBoundaryIntersections() unimplemented.");
}


// @@@ VSK, Change name when implementing
//===========================================================================
bool IntersectionPool::
checkIfBothPointsLieOnOneEndpointAndNotOfTheSameCurve()
  //===========================================================================
{
    // currently unimplemented
//    MESSAGE("Warning: IntersectionPool::"
//	    "checkIfBothPointsLieOnOneEndpointAndNotOfTheSameCurve() unimplemented");
    return false;
}


//===========================================================================
struct sort_point
//===========================================================================
{
  int idx1, idx2;
  sort_point(int dir1, int dir2)
  { 
      idx1 = dir1; 
      idx2 = dir2;
  }

  bool operator()(shared_ptr<IntersectionPoint> p1, 
		  shared_ptr<IntersectionPoint> p2) const
  {
      if (idx1 < 0 || idx2 < 0)
	  return true;
      else if (idx1 < 0)
      {
	  if (p1->getPar(idx2) <= p2->getPar(idx2))
	      return true;
	  else
	      return false;
      }
      else if (idx2 < 0)
      {
	  if (p1->getPar(idx1) <= p2->getPar(idx1))
	      return true;
	  else
	      return false;
      }
      else
      {
	  if (p1->getPar(idx1) < p2->getPar(idx1))
	      return true;
	  else if (p1->getPar(idx1) > p2->getPar(idx1))
	      return false;
	  else
	  {
	      if (p1->getPar(idx2) <= p2->getPar(idx2))
		  return true;
	      else
		  return false;
	  }
      }
	       
  }
};


//===========================================================================
void IntersectionPool::
getSortedIntersections(vector<shared_ptr<IntersectionPoint> >& int_pts)
  //===========================================================================
{
    int_pts.clear();
    int_pts.insert(int_pts.begin(), int_points_.begin(), int_points_.end());

    // Sort the output vector according to given parameter directions
    int dir1=-1, dir2=-1;
    int npar1 = obj1_->numParams();
    int npar2 = obj2_->numParams();
    if (npar1 == 1 && npar2 == 1)
    {
	dir1 = 0;
	dir2 = 1;
    }
    else if (npar1 == 1)
	dir1 = 0;
    else if (npar2 == 1)
	dir2 = npar1;
    sort_point compare(dir1, dir2);
    std::sort(int_pts.begin(), int_pts.end(), compare);

}


//===========================================================================
bool IntersectionPool::hasPointsInInner(int pardir)
  //===========================================================================
{
    double min_param, max_param;
    get_param_limits(pardir, min_param, max_param);

    for (int i = 0; i < int(int_points_.size()); ++i) {
	double par = int_points_[i]->getPar(pardir);
	double tol = int_points_[i]->parameterTolerance(pardir);
	if (par > min_param + tol && par < max_param - tol) {
	    return true;
	}
    }
    return false;
}


//===========================================================================
void IntersectionPool::
get_param_limits(int pardir, double& min_param, double& max_param) const 
//===========================================================================
{
    if (pardir < obj1_->numParams())
    {
	min_param = obj1_->startParam(pardir);
	max_param = obj1_->endParam(pardir);
    }
    else
    {
	min_param = obj2_->startParam(pardir - obj1_->numParams());
	max_param = obj2_->endParam(pardir - obj1_->numParams());
    }
}


//===========================================================================
vector<double> IntersectionPool::getSortedInnerInts(int pardir)
//===========================================================================
{
    vector<double> result;

    double min_param, max_param;
    get_param_limits(pardir, min_param, max_param);

    for (int i = 0; i < int(int_points_.size()); ++i) {
	double par = int_points_[i]->getPar(pardir);
	double tol = int_points_[i]->parameterTolerance(pardir);

	if (par > min_param + tol && par < max_param - tol) {
	    result.push_back(par);
	}
    }
    std::sort(result.begin(), result.end());
    return result;
}


//===========================================================================
shared_ptr<IntersectionPoint>
IntersectionPool::addIntersectionPoint(shared_ptr<ParamObjectInt> obj_int1_, 
				       shared_ptr<ParamObjectInt> obj_int2_,
				       shared_ptr<GeoTol> epsge,
				       double *par1, 
				       double *par2)
//===========================================================================
{
    shared_ptr<IntersectionPoint> 
	temp(new IntersectionPoint(obj_int1_.get(), obj_int2_.get(),
				   epsge, par1, par2));

    if (temp->getDist() >= epsge->getEpsge())
    {
	// Something is wrong
	MESSAGE("WARNING! : Illegal intersection point");
    }

    add_point_and_propagate_upwards(temp);
    return temp;
}


//===========================================================================
void IntersectionPool::
add_point_and_propagate_upwards(shared_ptr<IntersectionPoint> point)
//===========================================================================
{
    int_points_.push_back(point);
    if (lackingParameter() == -1 && prev_pool_.get() != 0) {
	// propagate point up to parent
	prev_pool_->add_point_and_propagate_upwards(point);
    }
}


//===========================================================================
void IntersectionPool::
includeReducedInts(shared_ptr<IntersectionPool> lower_order_pool)
//===========================================================================
{
    // lower_order_pool should be a child pool to 'this' pool
    ASSERT(lower_order_pool->prev_pool_.get() == this); 
    
    int lacking_ind = lower_order_pool->lackingParameter();
    double lacking_val = lower_order_pool->lackingParameterValue();
    // lower_order_pool should be a reduced pool 
    ASSERT(lacking_ind >= 0);
    
    associate_parent_points(lower_order_pool->getIntersectionPoints(),
			    lacking_ind, lacking_val);
    transfer_links_to_parent_points(lower_order_pool->getIntersectionPoints(),
				    lacking_ind);
}


//===========================================================================
void IntersectionPool::
selfIntersectParamReorganise(shared_ptr<IntersectionPool> sub_pool)
//===========================================================================
{
    ASSERT(obj1_ == obj2_); // function should only be called from
			    // selfintersection pools
    ASSERT(sub_pool->prev_pool_.get() == this); // it should be a
						// child pool to
						// 'this' pool
    ASSERT(sub_pool->lackingParameter() == -1); // number of
						// parameters should
						// be the same
    vector<shared_ptr<IntersectionPoint> > sub_pts
	= sub_pool->getIntersectionPoints();
   
    // keeping only points with twins.  If it has a 'parent', in this
    // case it is really a twin (no lacking parameter - hallmark of
    // "real" parent).
    vector<shared_ptr<IntersectionPoint> >::iterator new_end;
    new_end = remove_if(sub_pts.begin(), sub_pts.end(), no_parent);

    // finding all points that would lie on the same intersection
    // curve as a point with a twin.  As we want to make sure that
    // each such point is registered only once, a 'set' is safe to
    // use.
    set<IntersectionPoint*> points_to_flip;
    vector<shared_ptr<IntersectionPoint> >::iterator it;
    for (it = sub_pts.begin(); it != new_end; ++it) {
	add_reachables_from(it->get(), 0, points_to_flip);
    }

    // flipping points
    for_each(points_to_flip.begin(), points_to_flip.end(),
	     flip_intersecting_objects);

    // The twin points will now go out of scope.  But before that
    // happens, we will keep eventual links that they have with other
    // points.
    vector<shared_ptr<IntersectionLink> > links;
    IntersectionPoint *p1, *p2;
    for (it = sub_pts.begin(); it != new_end; ++it) {
	shared_ptr<IntersectionPoint> twin = (*it)->parentPoint();
	(*it)->getNeighbourLinks(links);
 	for (int j = 0; j < int(links.size()); ++j) {
 	    links[j]->getIntersectionPoints(p1, p2);
 	    if (it->get() == p1) {
		if (twin.get() != p2)
		    twin->connectTo(p2, LINK_UNDEFINED);
 	    } else if (it->get() == p2)
	    {
		if (twin.get() != p1)
		    twin->connectTo(p1, LINK_UNDEFINED);
 	    }
 	}
    }
}



//===========================================================================
void IntersectionPool::mirrorIntersectionPoints(size_t first)
//===========================================================================
{
    // For all intersection points, add a copy with opposite order of the
    // the parameter pairs. Copy also the connection between points
    if (obj1_->numParams() != 2 || obj2_->numParams() != 2)
	return;  // Not relevant
    if (obj1_.get() != obj2_.get())
	return;  // Not self intersection context

    size_t nmb_orig = int_points_.size();
    size_t ki, kj, kh;

    // First make a mirrored version of all the intersection points and add
    // them to the pool
    double par[4];
    for (ki=first; ki<nmb_orig; ki++)
    {
	for (int kr=0; kr<4; kr++)
	    par[kr] = int_points_[ki]->getPar(kr);
	addIntersectionPoint(obj1_, obj2_, int_points_[ki]->getTolerance(),
			     par+2, par);
    }

    // Copy links
    for (ki=first; ki<nmb_orig; ki++)
    {
	if (int_points_[ki]->numNeighbours() > 0)
	{
	    // Fetch neighbours
	    vector<IntersectionPoint*> next;
	    int_points_[ki]->getNeighbours(next);
	    for (kj=0; kj<next.size(); kj++)
	    {
		// Find the index of the neigbouring points in the vector
		for (kh=first; kh<nmb_orig; kh++)
		    if (next[kj] == int_points_[kh].get())
			break;

		if (kh == nmb_orig)
		    continue;  // Not possible to copy link

		shared_ptr<IntersectionLink> link = 
		    int_points_[ki]->getIntersectionLink(next[kj]);
		int_points_[nmb_orig+ki-first]->connectTo(int_points_[nmb_orig+kh-first],
						    link->linkType());
	    }
	}
    }	    
}

//===========================================================================
void IntersectionPool::
associate_parent_points(vector<shared_ptr<IntersectionPoint> >& children,
			int lacking_ix,
			double lacking_val)
//===========================================================================
{
    // This function will make sure that each of the
    // IntersectionPoints in the 'children' vector has exactly one
    // parent point in 'this' pool.
    vector<shared_ptr<IntersectionPoint> >& parents = int_points_;

    for (int c = 0; c < int(children.size()); ++c) {
	shared_ptr<IntersectionPoint>& child = children[c];
	shared_ptr<IntersectionPoint> cur_parent = child->parentPoint(); // could be 0

	if (cur_parent) {
	    // This child already has a parent, which we assume should
	    // already be present in the 'parents' vector.  The
	    // following loop exists only to check for this data
	    // consistency; it doesn't do anything, and could ideally
	    // be removed from the final program.
	    bool found = false;
	    for (int i = 0; i < int(parents.size()); ++i) {
		if (parents[i] == cur_parent) {
		    found = true;
		    break;
		}
	    }
	    ASSERT(found);
	} else {
	    // This IntersectionPoint has not yet defined a
	    // parent. Let us create one, and add it to the 'parents'
	    // vector.
	    
	    // preparing temporary parameter array
	    double par1[2];
	    int i, cpos;
	    for (i = 0, cpos = 0; i < obj1_->numParams(); ++i) {
		par1[i] = (i == lacking_ix) 
		    ? lacking_val : child->getPar1()[cpos++];
	    }
	    int lacking_ix2 = lacking_ix - obj1_->numParams();
	    double par2[2];
	    for (i = 0, cpos = 0; i < obj2_->numParams(); ++i) {
		par2[i] = (i == lacking_ix2) 
		    ? lacking_val : child->getPar2()[cpos++];
	    }

	    shared_ptr<IntersectionPoint>
		temp(new IntersectionPoint(obj1_.get(), 
					   obj2_.get(), 
					   child->getTolerance(), 
					   par1, 
					   par2));
	    child->setParentPoint(temp);
	    add_point_and_propagate_upwards(temp);
	}
    }
}


//===========================================================================
void IntersectionPool::
transfer_links_to_parent_points(vector<shared_ptr<IntersectionPoint> >&
				children,
				int lacking_ix)
//===========================================================================
{
    for (int c = 0; c < int(children.size()); ++c) {
	shared_ptr<IntersectionPoint> cur_child = children[c];
	vector<IntersectionPoint*> neighbours;
	cur_child->getNeighbours(neighbours);
	shared_ptr<IntersectionPoint> child_parent
	    = cur_child->parentPoint();
	ASSERT(child_parent); // parent should have been provided
			      // before calling this routine
	for (int n = 0; n < int(neighbours.size()); ++n) {
	    shared_ptr<IntersectionPoint> neigh_parent
		= neighbours[n]->parentPoint();
	    if (!neigh_parent) {
		// optional test for consistency.
		// Neighbour did not have a parent point.  This must
		// mean that it lies in a "higher" pool, in which case
		// the link should already exist.  This test verifies
		// that this is the case
		//ASSERT(child_parent->isConnectedTo(neighbours[n]));
	    } else {
		if (!child_parent->isConnectedTo(neigh_parent)) {
		    shared_ptr<IntersectionLink> link
			= cur_child->getIntersectionLink(neighbours[n]);
		    LinkType linktype = link->linkType();
		    child_parent->connectTo(neigh_parent, 
					    linktype,
					    link,
					    lacking_ix);
		}
		shared_ptr<IntersectionLink> cl
		    = child_parent->getIntersectionLink(neigh_parent.get());
		ASSERT(cl);
		cl->setIsoparametricIn(lacking_ix, true); 
	    }
	}
    }
}


//===========================================================================
int IntersectionPool::
inInfluenceArea(int pardir, double par, bool first_outside)
//===========================================================================
{
    int result = 0;
    vector<shared_ptr<IntersectionPoint> >::const_iterator it;
    for (it = int_points_.begin(); it != int_points_.end(); ++it) {
	int tmp = (*it)->inInfluenceArea(pardir, par, first_outside);
	if (tmp == 1) {
	    // in influence area, but not exactly on point
	    result = 1;
	} else if (tmp == 2) {
	    // exactly on point
	    result = (result == 1) ? 1 : 2;
	}
    }    
    return result;
}


//===========================================================================
int IntersectionPool::
inInfluenceArea(int pardir, 
		double par, 
		vector<shared_ptr<IntersectionPoint> >& int_pts, 
		bool first_outside)
//===========================================================================
{
    int result = 0;
    vector<shared_ptr<IntersectionPoint> >::const_iterator it;
    for (it = int_points_.begin(); it != int_points_.end(); ++it) {
	int tmp = (*it)->inInfluenceArea(pardir, par, first_outside);
	if (tmp > 0) {
	    int_pts.push_back(*it);
	}
	if (tmp == 1) {
	    // in influence area, but not exactly on point
	    result = 1;
	} else if (tmp == 2) {
	    // exactly on point
	    result = (result == 1) ? 1 : 2;
	}
    }    
    return result;
    
}


//===========================================================================
void IntersectionPool::getMidParameter(double *mid)
//===========================================================================
{
    int i;
    if (int_points_.size() > 0) {
	int num_params
	    = int_points_[0]->numParams1() + int_points_[0]->numParams2();
	for (i = 0; i < num_params; ++i) {
	    mid[i] = 0;
	}
	vector<shared_ptr<IntersectionPoint> >::const_iterator it;
	for (it = int_points_.begin(); it != int_points_.end(); ++it) {
	    const vector<double>& pars = (*it)->getPar();
	    for (int i = 0; i < num_params; ++i) {
		mid[i] += pars[i];
	    }
	}
	for (i = 0; i < num_params; ++i) {
	    mid[i] /= (int)int_points_.size();
	}
    } else {
	// no intersection points in pool.  Return middle value of
	// parameters in objects
	int num_p_1 = obj1_->numParams();
	int num_p_2 = obj2_->numParams();
	int num_params = num_p_1 + num_p_2;
	for (int i = 0; i < num_params; ++i) {
	    if (i < num_p_1) {
		mid[i] = (obj1_->startParam(i) + obj1_->endParam(i)) * 0.5;
	    } else {
		int j = i - num_p_1;
		mid[i] = (obj2_->startParam(j) + obj2_->endParam(j)) * 0.5;
	    }
	}
    }
}


//===========================================================================
void IntersectionPool::
getResult(vector<shared_ptr<IntersectionPoint> >& int_points,
	  vector<shared_ptr<IntersectionCurve> >& int_curves) 
//===========================================================================
{
    // @@@jbt: Note: A preprocessing step has been moved to
    // makeIntersectionCurves()

    int_points.clear();
    for (size_t ki=0; ki<int_points_.size(); ki++)
	if (int_points_[ki]->numNeighbours() == 0)
	    int_points.push_back(int_points_[ki]);

    int_curves = int_curves_;
}


//===========================================================================
void IntersectionPool::
getSortedBdInts(const Point& vec, 
		vector<shared_ptr<IntersectionPoint> >& result,
		int obj_nmb)
//===========================================================================
{
    // NOTE: We get _all_ intersection points - not just those that
    // lie on the boundary.

    result = int_points_;
    int num_points = (int)int_points_.size();
    vector<pair<double, shared_ptr<IntersectionPoint> > > proj(num_points);
    for (int i = 0; i < num_points; ++i) {
	if (vec.dimension() == 3) {
	    proj[i].first = result[i]->getPoint() * vec;
	} else if (obj_nmb < 0) {
	    Point par1(result[i]->getPar1(), result[i]->getPar1() + 2);
	    proj[i].first = par1 * vec;
	}
	else {
	    Point par1(result[i]->getPar1()+2*obj_nmb, result[i]->getPar1()+2*obj_nmb + 2);
	    proj[i].first = par1 * vec;
	}
	    
	proj[i].second = result[i];
    }

    // Sort pairs in vector proj based on the first element
    sort(proj.begin(), proj.end(),
	 compare_first<double, shared_ptr<IntersectionPoint> >);

    for (int i = 0; i < num_points; ++i) {
	result[i] = proj[i].second;
    }
}


//===========================================================================
bool IntersectionPool::
checkSortingDir(Point& vec, int sorting_obj)
//===========================================================================
{
    int no_turn = 0, turn = 0;
    if (int_points_.size() == 0)
	return true;

    for (size_t ki=0; ki<int_points_.size(); ki++)
    {
	SingularityType type = int_points_[ki]->getSingularityType();
// 	if (type == HIGHER_ORDER_POINT || type == ISOLATED_POINT)
// 	    continue;  // Can't check the direction of the intersection curve
	if (type != ORDINARY_POINT)
	    continue;   // No orientation specified

	/*// Get direction of intersection curve
	Point tang = int_points_[ki]->getTangent();
	tang.normalize();
	
	// Evaluate object
	vector<Point> der(3);
	if (sorting_obj == 0)
	    int_points_[ki]->getObj1()->point(der, int_points_[ki]->getPar1(), 1);
	else
	    int_points_[ki]->getObj2()->point(der, int_points_[ki]->getPar2(), 1); 
	Point dir = der[1]*vec[0] + der[2]*vec[1];
	dir.normalize();

	Point tang2D = (sorting_obj == 0) ? int_points_[ki]->getPar1Dir() : int_points_[ki]->getPar2Dir();

	// Check consistency
	double tdir = tang*dir;*/

	// VSK, 070516 Do the computations in the parameter domain
	Point tang2D = (sorting_obj == 0) ? int_points_[ki]->getPar1Dir() : int_points_[ki]->getPar2Dir();
	double tdir = vec*tang2D;
	double tol = int_points_[ki]->getTolerance()->getAngleTol();
	if (tdir < -tol)
	    turn++;
	else if (tdir > tol)
	    no_turn++;
    }
    if (turn > no_turn)
    {
	vec *= -1;
    }

    return ((turn > 0 && no_turn == 0) || (no_turn > 0 && turn == 0));
}

//===========================================================================
bool IntersectionPool::existIntersectionPoint(int dir, double par)
//===========================================================================
{
    vector<shared_ptr<IntersectionPoint> >::const_iterator it;
    for (it = int_points_.begin(); it != int_points_.end(); ++it) {
	double tol = (*it)->getTolerance()->getRelParRes();
	double point_par = (*it)->getPar(dir);
	if (fabs(point_par - par) < tol) {
	    return true;
	}
    }
    return false;
}


//===========================================================================
void IntersectionPool::
getIntersectionPoints(int dir, 
		      double par, 
		      vector<shared_ptr<IntersectionPoint> >& result) const
//===========================================================================
{
    result.clear();    
    vector<shared_ptr<IntersectionPoint> >::const_iterator it;
    for (it = int_points_.begin(); it != int_points_.end(); ++it) {
	double tol = (*it)->getTolerance()->getRelParRes();
	double point_par = (*it)->getPar(dir);
	if (fabs(point_par - par) < tol) {
	    result.push_back(*it);
	}
    }
}


//===========================================================================
bool IntersectionPool::hasIntersectionPoints(int dir, 
					     double par) const
//===========================================================================
{
    vector<shared_ptr<IntersectionPoint> > result;
    getIntersectionPoints(dir, par, result);
    return (result.size() > 0);
}


//===========================================================================
bool IntersectionPool::
checkIntersectionPoints(vector<IntPtInfo>& int_pt_info) const
//===========================================================================
{
    // Loop through all intersection points in the pool. For each
    // point, get the number of neighbours and the singularity type
    // and check if these are consistent. Put this information in a
    // struct IntPtInfo.

    // Only implemented for two surfaces
    shared_ptr<ParamSurfaceInt> obj1
	= dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(obj1_);
    shared_ptr<ParamSurfaceInt> obj2
	= dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(obj2_);
    if (obj1.get() == 0 || obj2.get() == 0)
	return false;

    int npoints = (int)int_points_.size();
    int_pt_info.resize(npoints);
    if (npoints == 0)
	return true;

    for (int i = 0; i < npoints; ++i) {
	int_points_[i]->checkIntersectionPoint(int_pt_info[i]);
    }

    // Return 'true' only if all points are OK
    for (int i = 0; i < npoints; ++i) {
	if (!int_pt_info[i].is_ok) {
	    return false;
	}
    }
    return true;
}


//===========================================================================
bool IntersectionPool::checkIntersectionLinks() const
//===========================================================================
{
    // Loop through all intersection links in the pool, and check if
    // they are valid

    // Only implemented for two surfaces
    shared_ptr<ParamSurfaceInt> obj1
	= dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(obj1_);
    shared_ptr<ParamSurfaceInt> obj2
	= dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(obj2_);
    if (obj1.get() == 0 || obj2.get() == 0)
	return false;

    // Defining an std::set with unique links. Note that we are
    // implicitly "ordering" links using operator<() in shared_pts,
    // since this supplies us with a strict weak ordering.
    int npoints = (int)int_points_.size();
    set<shared_ptr<IntersectionLink> > links;
    vector<shared_ptr<IntersectionLink> > neighbour_links;
    for (int i = 0; i < npoints; ++i) {
	int_points_[i]->getNeighbourLinks(neighbour_links);
	links.insert(neighbour_links.begin(), neighbour_links.end());
    }

    // Write the links
    cout << "*** IntersectionLink check ***" << endl;
    int linknum = 0;
    IntersectionPoint *p1, *p2;
    Point tangent1, tangent2;
    typedef set<shared_ptr<IntersectionLink> >::iterator iter;
    for (iter it = links.begin(); it != links.end(); ++it) {
	(*it)->getIntersectionPoints(p1, p2);
	// Check if the tangents at the ends are consistent
	try {
	    tangent1 = p1->getTangent();
	    tangent2 = p2->getTangent();
	} catch (...)
	{
	    continue;
	}
	double cos = tangent1 * tangent2;
	bool is_ok = verifyIntersectionLink(*it);
	cout << linknum << "\t";
	p1->writeParams(cout);
	cout << " -- ";
	p2->writeParams(cout);
	cout << "    cos = " << cos
	     << "    " << (is_ok ? "OK" : "FAILED") << endl;
	++linknum;
    }


    return true;
}


//===========================================================================
void IntersectionPool::removeDoublePoints()
//===========================================================================
{
    // Loop through all intersection points. Check for each point if
    // there are other points that is really the same point within the
    // tolerance. Remove these other points, and transfer the links to
    // the one surviving point.

    typedef vector<shared_ptr<IntersectionPoint> >::iterator iter;
    iter it = int_points_.begin();
    while (it != int_points_.end()) {
	// If the point lies on a degenerate edge we do nothing
	if ((*it)->isDegenerate()) {
	    ++it;
	    continue;
	}
	iter jt = it + 1;
	while (jt != int_points_.end()) {
	    bool is_same_point = (*it)->isSamePoint((*jt).get());
	    if (is_same_point) {
		// The points are the same. We first transfer the
		// links.
		int nneighbours = (*jt)->numNeighbours();
		vector<IntersectionPoint*> neighbours;
		(*jt)->getNeighbours(neighbours);
		for (int i = 0; i < nneighbours; ++i) {
		    (*jt)->disconnectFrom(neighbours[i]);
		    // Do not reconnect to itself
		    if (neighbours[i] != (*it).get()) {
			neighbours[i]->connectTo(*it, LINK_UNDEFINED);
		    }
		}
		removeIntPoint(*jt);
	    }
	    else {
		++jt;
	    }
	}
	++it;
    }

    // VSK, 0607. Make a second pass and check if neighbouring points
    // represent the same intersection
    size_t kr;
    for (kr=0; kr<int_points_.size(); kr++)
    {
	// If the point lies on a degenerate edge we do nothing
	if (int_points_[kr]->isDegenerate()) 
	    continue;

	// Check intersection links
	vector<shared_ptr<IntersectionLink> > links;
	int_points_[kr]->getNeighbourLinks(links);
	for (int j=0; j<(int)links.size(); ++j)
	{
	    bool is_ok = verifyIntersectionLink(links[j]);
	    if (!is_ok)
	    {
		// Check if the points point in the same direction and is
		// connected to the same point
		IntersectionPoint *p1, *p2;
		links[j]->getIntersectionPoints(p1, p2);
		if ((p1->getSingularityType() != ORDINARY_POINT &&
		     p1->getSingularityType() != TANGENTIAL_POINT) ||
		    (p2->getSingularityType() != ORDINARY_POINT &&
		     p2->getSingularityType() != TANGENTIAL_POINT))
		    continue;  // Link cannot be tested with this tool

		Point tangent1 = p1->getTangent();
		Point tangent2 = p2->getTangent();
		double cos = tangent1.cosAngle(tangent2);
		if (cos > 0.9)
		{
		    // The tangents point in roughly the same direction
		    // Get the neighbours of the two points
		    vector<IntersectionPoint*> npnts1;
		    vector<IntersectionPoint*> npnts2;
		    p1->getNeighbours(npnts1);
		    p2->getNeighbours(npnts2);
		    size_t k1, k2;
		    for (k1=0; k1<npnts1.size(); k1++)
		    {
			for (k2=0; k2<npnts2.size(); k2++)
			{
			    if (npnts1[k1] != p2 && npnts2[k2] != p1 &&
				npnts1[k1] == npnts2[k2])
			    {
				// Same neighbour
				break;
			    }
			}
			if (k2 < npnts2.size())
			    break;
			}

		    if (k1 < npnts1.size())
		    {
			// Double point. Remove the less exact one
			if (p1->getDist() < p2->getDist())
			{
			    size_t kh;
			    for (kh=0; kh<int_points_.size(); kh++)
				if (int_points_[kh].get() == p2)
				    break;
			    if (kh < int_points_.size())
			    {
				removeIntPoint(int_points_[kh]);
			    }
			}
			else
			{
			    size_t kh;
			    for (kh=0; kh<int_points_.size(); kh++)
				if (int_points_[kh].get() == p1)
				    break;
			    if (kh < int_points_.size())
			    {
				removeIntPoint(int_points_[kh]);
			    }
			}
		    
		    }
		}
	    }
	}
    }
    return;
}


//===========================================================================
void IntersectionPool::removeLooseEnds()
//===========================================================================
{
    // Remove ordinary intersection points with only one neighbour when these
    // points lie reasonably close, their tangents point roughly in the same
    // direction, and the neighbour has significantly better accuracy.
    int kr;
    double fac = 1.0e5;
    double dist_fac = 10.0;
    for (kr=0; kr<(int)int_points_.size(); kr++)
    {
	if (isBoundaryPoint(int_points_[kr]))
	    continue;

	if (int_points_[kr]->getSingularityType() != ORDINARY_POINT)
	    continue;

	if (int_points_[kr]->numNeighbours() != 1)
	    continue;

	// Get neighbours
	vector<IntersectionPoint*> neighbours;
	int_points_[kr]->getNeighbours(neighbours);

	Point pnt1 = int_points_[kr]->getPoint();
	Point pnt2 = neighbours[0]->getPoint();
	double dist = pnt1.dist(pnt2);

	if (dist > fac*int_points_[kr]->getTolerance()->getEpsge())
	    continue;

	// Check link
	shared_ptr<IntersectionLink> link = int_points_[kr]->getIntersectionLink(neighbours[0]);
	bool is_ok = verifyIntersectionLink(link);
	if (is_ok)
	    continue;

	double d1 = int_points_[kr]->getDist();
	double d2 = neighbours[0]->getDist();
	if (d1 > dist_fac*d2)
	{
	    // A loose end is identified. Remove
	    removeIntPoint(int_points_[kr]);
	    kr--;
	}
	else if (d2 > dist_fac*d1)
	{
	    // A loose end is identified. Remove
	    size_t kh;
	    for (kh=0; kh<int_points_.size(); kh++)
	    {
		if (int_points_[kh].get() == neighbours[0])
		    break;
	    }
	    removeIntPoint(int_points_[kh]);
	    if ((int)kh < kr)
		kr--;
	}

    }
	
}

//===========================================================================
void IntersectionPool::removeFalseCorners()
//===========================================================================
{
    // Check ordinary intersection points to check if the links they belong to
    // are consistent
    size_t kr;
    double fac = 0.2;
    for (kr=0; kr<int_points_.size(); kr++)
    {
 	/*if (isBoundaryPoint(int_points_[kr]))
	  continue;*/

	if (int_points_[kr]->getSingularityType() != ORDINARY_POINT)
	    continue;

	if (int_points_[kr]->numNeighbours() != 2)
	    continue;

	// Check consistency of the entire chain. Ordinary points and
	// tangential points are expected to have one neighbour
	bool chain_ok = checkIntersectionChain(int_points_[kr].get(), int_points_[kr].get());
	if (chain_ok)
	    continue;  

	// Get neighbours
	vector<IntersectionPoint*> neighbours;
	int_points_[kr]->getNeighbours(neighbours);

	Point curr = int_points_[kr]->getPoint();
	Point pnt1 = neighbours[0]->getPoint();
	Point pnt2 = neighbours[1]->getPoint();
	Point tan = int_points_[kr]->getTangent();

	Point diff1 = pnt1 - curr;
	Point diff2 = curr - pnt2;
	
	// If the direction of the curve changes completely, something is wrong
	if (diff1*diff2 < 0.0)
	{
	    // Remove one branch. Compute angle between the difference vector and the
	    // curve tangent
	    double ang1 = diff1.angle_smallest(tan);
	    double ang2 = diff2.angle_smallest(tan);

	    // Compute distance to next point
	    double len1 = diff1.length();
	    double len2 = diff2.length();

	    // We want to keep the link with small angle and small distance
	    bool keep_first;
	    if (ang1 < fac*ang2 && ang2 > int_points_[kr]->getTolerance()->getAngleTol())
		keep_first = true;
	    else if (ang2 < fac*ang1 && ang1 > int_points_[kr]->getTolerance()->getAngleTol())
		keep_first = false;
	    else if (len1 < len2)
		keep_first = true;
	    else
		keep_first = false;

	    if (keep_first)
		int_points_[kr]->disconnectFrom(neighbours[1]);
	    else
		int_points_[kr]->disconnectFrom(neighbours[0]);
	}
   }
}

//===========================================================================
void IntersectionPool::removeDefectLinks()
//===========================================================================
{
    // Loop through all intersection links in the pool. Check the
    // validity of each link with verifyIntersectionLink(). Remove it
    // is not valid.

    // This produces some debug writeout
    checkIntersectionLinks();

    int npoints = (int)int_points_.size();
    IntersectionPoint *p1, *p2;
    vector<shared_ptr<IntersectionLink> > neighbour_links;
    for (int i = 0; i < npoints; ++i) {
	// Check consistency of the entire chain. Ordinary points and
	// tangential points are expected to have one neighbour
	bool chain_ok = checkIntersectionChain(int_points_[i].get(), int_points_[i].get());
	if (chain_ok)
	    continue;  

	int_points_[i]->getNeighbourLinks(neighbour_links);
	int nneighbours = (int)neighbour_links.size();
	for (int j = 0; j < nneighbours; ++j) {
	    bool is_ok = verifyIntersectionLink(neighbour_links[j]);
	    if (!is_ok) {
		neighbour_links[j]->getIntersectionPoints(p1, p2);
		p1->disconnectFrom(p2);
	    }
	}
    }

    return;
}

//===========================================================================
bool IntersectionPool::checkIntersectionChain(IntersectionPoint *pnt,
					      IntersectionPoint *first,
					      IntersectionPoint *prev) 
//===========================================================================
{
    SingularityType type = pnt->getSingularityType();
    int nmb = pnt->numNeighbours();

    // Look for inconsistent configuaration
    if (type == ISOLATED_POINT && nmb > 0)
	return false;

    if ((type == ORDINARY_POINT || type == TANGENTIAL_POINT) && nmb > 2)
	return false;

    if ((type == ORDINARY_POINT || type == TANGENTIAL_POINT) && 
	(nmb == 1 && !isBoundaryPoint(pnt)))
	return false;
    
    // Look for end point of the chain
    if ((type == BRANCH_POINT || type == HIGHER_ORDER_POINT) && prev != 0)
	return true;  // Chain OK so far

    if (isBoundaryPoint(pnt) && prev != 0)
	return true;

    if (pnt == first && prev != 0)
	return true;   // Closed loop

    // Traverse the chain in all directions
    vector<IntersectionPoint*> neighbours;
    pnt->getNeighbours(neighbours);
    bool chain_ok = true;
    for (size_t ki=0; ki<neighbours.size(); ki++)
    {
	if (neighbours[ki] == prev)
	    continue;  // This part of the chain is already treated
	
	chain_ok = checkIntersectionChain(neighbours[ki], first, pnt);
	if (!chain_ok)
	    return false;
    }

    return true;  // Chain found OK
}
    

//===========================================================================
void IntersectionPool::getAllLinks(set<shared_ptr<IntersectionLink> >& links)
//===========================================================================
{
    int npoints = (int)int_points_.size();
    vector<shared_ptr<IntersectionLink> > neighbour_links;
    for (int i = 0; i < npoints; ++i) {
	int_points_[i]->getNeighbourLinks(neighbour_links);
	links.insert(neighbour_links.begin(), neighbour_links.end());
    }
}

//===========================================================================
bool IntersectionPool::isInDomain(IntersectionPoint *pnt) const
//===========================================================================
{
    for (size_t ki=0; ki<int_points_.size(); ki++)
	if (int_points_[ki].get() == pnt)
	    return true;

    return false;
}


//===========================================================================
bool IntersectionPool::
closestInDomain(double param[], 
		shared_ptr<IntersectionPoint>& pnt) const
//===========================================================================
{
    if (int_points_.size() == 0)
	return false;

    double mindist = 1.0e6;  // A large number
    double dist;
    int minind = -1;
    for (size_t ki=0; ki<int_points_.size(); ki++)
    {
	vector<double> par = int_points_[ki]->getPar();
	dist = 0.0;
	for (int kj=0; kj<(int)(par.size()); kj++)
	    dist += (par[kj] - param[kj])*(par[kj] - param[kj]);
	if (dist < mindist)
	{
	    mindist = dist;
	    minind = (int)ki;
	}
    }
    if (minind < 0)
	return false;

    pnt = int_points_[minind];
    return true;
}


//===========================================================================
bool IntersectionPool::
closestInDomain(const double param[], 
		shared_ptr<IntersectionPoint>& pnt) const
//===========================================================================
{
    double par[4]; // Up to four parameters
    for (int i = 0; i < numParams(); ++i)
	par[i] = param[i];

    return closestInDomain(par, pnt);
}


// //===========================================================================
// void IntersectionPool::weedOutClutterCurves(vector<vector<int> >& curves)  const
// //===========================================================================
// {
//     // we consider an open curve to be 'clutter' if it cannot meet the following criterion:
//     // -> the start and endpoints must EITHER:  * lie on the edge of the parameter domain
//     //                                     OR:  * be a non-transversal point    
//     // This function removes clutter curves from the 'curves' vector
//     vector<vector<int> > remaining_curves;
    
//     for (int i = 0; i < int(curves.size()); ++i) {
// 	// checking if curves[i] is clutter
// 	bool is_clutter = false;
// 	shared_ptr<IntersectionPoint> temp_point = int_points_[curves[i].front()];

// 	// checking start point of curve
// 	bool on_first_boundary, on_second_boundary;
// 	if (!temp_point->pointIsSingular()) {
// 	    // this is an ordinary point.  It must lie on the boundary in order not
// 	    // to be considered 'clutter'.
// 	    temp_point->isBoundaryPoint(on_first_boundary, on_second_boundary);
// 	    is_clutter = !(on_first_boundary || on_second_boundary);
// 	}
// 	if (!is_clutter) { // checking endpoint of curve
// 	    temp_point = int_points_[curves[i].back()];
// 	    if (!temp_point->pointIsSingular()) {
// 		// this is an ordinary point.  It must lie on the boundary in order not
// 		// to be considered 'clutter'.
// 		temp_point->isBoundaryPoint(on_first_boundary, on_second_boundary);
// 		is_clutter = !(on_first_boundary || on_second_boundary);
// 	    }
// 	}
// 	// we have checked both endpoints.  If the curve has not yet been judged to be
// 	// clutter, we will keep it
// 	if (!is_clutter) {
// 	    remaining_curves.push_back(curves[i]);
// 	}
//     }

//     curves.swap(remaining_curves);
// }


//===========================================================================
void IntersectionPool::weedOutClutterPoints() 
//===========================================================================
{
    // going through the intersection points in the pool and removing
    // those that are ordinary points with only one neighbour and that
    // does not lie on the boundary.
//    vector<shared_ptr<IntersectionPoint> > kept_int_points;
    vector<IntersectionPoint*> neigh;

    for (int i = 0; i < int(int_points_.size()); ++i) {
	bool keep_point = true;
	shared_ptr<IntersectionPoint> temp_point = int_points_[i];
// 	bool not_singular
// 	    = (temp_point->getSingularityType() == ORDINARY_POINT
// 	       || temp_point->getSingularityType() == TANGENTIAL_POINT);
// 	if (not_singular
	if (!temp_point->pointIsSingular()
	    && temp_point->numNeighbours() == 1) {
	    // this is an ordinary point with only one neighbour.
	    // Remove if it not a boundary point.
	    //bool boundary1, boundary2;
	    //temp_point->isBoundaryPoint(boundary1, boundary2);
	    //keep_point = boundary1 || boundary2;
	    keep_point = isBoundaryPoint(temp_point);
	}
	if (keep_point) {
//	    kept_int_points.push_back(temp_point);
	    ;
	} else {
	    // disconnecting this point from its neighbour, and omit
	    // it from the new vector.
	    /*temp_point->getNeighbours(neigh);
	    ASSERT(neigh.size() == 1);
	    temp_point->disconnectFrom(neigh[0]);*/
	    removeIntPoint(temp_point);
	}
    }
//    kept_int_points.swap(int_points_);
}


//===========================================================================
void IntersectionPool::makeIntersectionCurves()
//===========================================================================
{
    DEBUG_ERROR_IF(prev_pool_,
		   "makeIntersectionCurves() should only be called on "
		   "parent pool.");
    int_curves_.clear(); // just in case

    vector<vector<int> > open_curve_pt_indices;
    vector<vector<int> > closed_curve_pt_indices;
    vector<int> isolated_points; // unused here

    //weedOutClutterPoints(); // remove intersection points that are
			    // useless or harmful

    get_individual_paths((int)int_points_.size(),
			 ConnectionFunctor(int_points_),
			 open_curve_pt_indices,
			 closed_curve_pt_indices,
			 isolated_points);

    if (getenv("DEBUG_DEMO") && *(getenv("DEBUG_DEMO")) == '1')
    {
	try {
	    writeDebug();
	} catch (...) {
	    MESSAGE("Failed printing debug info, continuing.");
	}
    }

    // splitting those open curves that consists of multiple
    // isoparametric segments (optimization issue) NB: Currently
    // commented out because judged indesireable after all.
    // separate_isoparametric_curves(open_curve_pt_indices);

    // adding open curves
    vector<shared_ptr<IntersectionPoint> > tempvec;
    for (int i = 0; i < int(open_curve_pt_indices.size()); ++i) {
	int num_ipoints = (int)open_curve_pt_indices[i].size();
	tempvec.resize(num_ipoints);
	for (int j = 0; j < num_ipoints; ++j) {
	    tempvec[j] = int_points_[open_curve_pt_indices[i][j]];
	}

	int_curves_.push_back(constructIntersectionCurve(tempvec.begin(),
							 tempvec.end()));
// 	int_curves_.push_back(shared_ptr<IntersectionCurve>
// 			      (new IntersectionCurve(tempvec.begin(), 
// 						     tempvec.end())));
    }

    // @@@ NB: Add support for closed curves!  As a temporary
    // solution, we will add open curves instead.
    /*MESSAGE("Closed intersection curves not yet supported.  "
      "Making open curves instead.");*/
    for (int i = 0; i < int(closed_curve_pt_indices.size()); ++i) {
	// adding start point again
	closed_curve_pt_indices[i].push_back(closed_curve_pt_indices[i][0]);
    }

    // NB: Currently commented out because judged indesireable after
    //all.  separate_isoparametric_curves(closed_curve_pt_indices);

    // adding open curves
    for (int i = 0; i < int(closed_curve_pt_indices.size()); ++i) {
	int num_ipoints = (int)closed_curve_pt_indices[i].size();
	tempvec.resize(num_ipoints);
	for (int j = 0; j < num_ipoints; ++j) {
	    tempvec[j] = int_points_[closed_curve_pt_indices[i][j]];
	}
	int_curves_.push_back(constructIntersectionCurve(tempvec.begin(),
							 tempvec.end()));
// 	int_curves_.push_back(shared_ptr<IntersectionCurve>
// 			      (new IntersectionCurve(tempvec.begin(), 
// 						     tempvec.end())));
    }

    /*
      VSK, 0701 DONT!
    // @@@jbt: This part is put here due to changes in getResult().
    vector<shared_ptr<IntersectionPoint> > tmp_points;
    for (size_t ki = 0; ki < int_points_.size(); ki++) {
	if (int_points_[ki]->numNeighbours() == 0)
	    tmp_points.push_back(int_points_[ki]);
    }
    int_points_ = tmp_points;
    */

}


//===========================================================================
void IntersectionPool::
separate_at_problematic_points(vector<vector<int> >& curve_indices) const
//===========================================================================
{
    int orig_num_curves = (int)curve_indices.size();
    for (int c = 0; c < orig_num_curves; ++c) {
	vector<int>& cur_cv = curve_indices[c];

	int end = 0;
	vector<int> split_pt;
	split_pt.push_back(0);
	
	for (end = 1; end < int(cur_cv.size() - 1); ++end) {
	    shared_ptr<const IntersectionPoint> ip = int_points_[cur_cv[end]];
	    if (!ip->hasUniqueTangentDirection()) {
		// this is a 'problem' point
		split_pt.push_back(end);
	    }
	}
	split_pt.push_back(int(cur_cv.size()) - 1);

	// making extra curve segments (those in addition to the first
	// one)
	for (int i = 1; i < int(split_pt.size()) - 1; ++i) {
	    curve_indices.push_back(vector<int>(curve_indices[c].begin()
						+ split_pt[i],
						curve_indices[c].begin()
						+ split_pt[i+1] + 1));
	}

	// resizing first curve to represent first curve segment
	curve_indices[c].resize(split_pt[1] + 1);
    }
}


//===========================================================================
void IntersectionPool::
separate_isoparametric_curves(vector<vector<int> >& curve_indices) const 
//===========================================================================
{
    int orig_num_curves = (int)curve_indices.size();
    for (int c = 0; c < orig_num_curves; ++c) {
	vector<int>& cur_cv = curve_indices[c];

	int end = 0;
	vector<int> split_pt;
	split_pt.push_back(0); // start ix. of first point in curve
	
	for (end = 0; end < int(cur_cv.size()) - 1; ) {

	    IntersectionPoint *ip1 = int_points_[cur_cv[end]].get();
	    IntersectionPoint *ip2 = int_points_[cur_cv[end + 1]].get();
	    shared_ptr<IntersectionLink> link = ip1->getIntersectionLink(ip2);
	    ASSERT(link);
	    end = link->isIsoparametric()
		? endof_isoparametric_curve(cur_cv, end)
		: endof_noniso_curve(cur_cv, end);

	    ASSERT(end > split_pt.back());
	    split_pt.push_back(end);
	}
	ASSERT(split_pt.size() >= 2);

	// making extra curve segments (those in addition to the first
	// one)
	for (int i = 1; i < int(split_pt.size()) - 1; ++i) {
	    curve_indices.
		push_back(vector<int>(curve_indices[c].begin()
				      + split_pt[i],
				      curve_indices[c].begin()
				      + split_pt[i+1] + 1));
	}
	
	// resizing first curve to represent first curve segment
	curve_indices[c].resize(split_pt[1]+1);
    }
}	


//===========================================================================
int IntersectionPool::
endof_isoparametric_curve(const vector<int>& cur_cv, 
			  int start_ix) const 
//===========================================================================
{
    // this function can be elaborated later to take into
    // consideration all the isoparameters
    
    // locating isoparameter
    IntersectionPoint* ip1 = int_points_[cur_cv[start_ix]].get();
    IntersectionPoint* ip2 = int_points_[cur_cv[start_ix + 1]].get();
    shared_ptr<IntersectionLink> link = ip1->getIntersectionLink(ip2);
    ASSERT(link);
    
    int isopar;
    for (isopar = 0; isopar < link->numParams(); ++isopar) {
	if (link->isIsoparametricIn(isopar)) {
	    break;
	}
    }
    ASSERT(isopar < link->numParams()); // 'break' should have been
				       // called in abv. loop
    
    int end_ix;
    for (end_ix = start_ix + 1; end_ix < int(cur_cv.size()) - 1; ++end_ix) {
	ip1 = int_points_[cur_cv[end_ix]].get();
	ip2 = int_points_[cur_cv[end_ix+ 1]].get();
	link = ip1->getIntersectionLink(ip2); ASSERT(link);
	if (!link->isIsoparametricIn(isopar)) {
	    // curve no longer isoparametric in this direction
	    break;
	}
    }
    return end_ix;
}


//===========================================================================
int IntersectionPool::
endof_noniso_curve(const vector<int>& cur_cv, int start_ix) const 
//===========================================================================
{
    int end_ix;
    for (end_ix = start_ix + 1; end_ix < int(cur_cv.size()) - 1; ++end_ix) {
	// checking if we have arrived at an isoparametric segment
	IntersectionPoint* ip1 = int_points_[cur_cv[end_ix]].get();
	IntersectionPoint* ip2 = int_points_[cur_cv[end_ix + 1]].get();
	shared_ptr<IntersectionLink> link = ip1->getIntersectionLink(ip2);
	ASSERT(link);
	if (link->isIsoparametric()) {
	    // we have arrived at an isoparametric part of the curve
	    // split the curve here
	    break;
	}
    }
    return end_ix;
}


//===========================================================================
void IntersectionPool::
getOrigPoints(vector<shared_ptr<IntersectionPoint> >& int_pts) const
//===========================================================================
{
    if (lackingParameter() >= 0 || prev_pool_.get() == 0)
	getIntersectionPoints(int_pts);
    else
	prev_pool_->getOrigPoints(int_pts);
}


//===========================================================================
bool IntersectionPool::
verifyIntersectionLink(const shared_ptr<IntersectionLink>& link,
		       int recursion_limit) const
//===========================================================================
{
    // Verify that the given intersection link represents a connected
    // piece of the intersection curve.

    // Works only for surface-surface
    shared_ptr<ParamSurfaceInt> obj1
	= dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(obj1_);
    shared_ptr<ParamSurfaceInt> obj2
	= dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(obj2_);
    if (obj1.get() == 0 || obj2.get() == 0) {
	MESSAGE("Only implemented for surface-surface");
	return false;
    }

    // Get the points and check if they are in the pool
    IntersectionPoint *p1, *p2;
    link->getIntersectionPoints(p1, p2);
    if (!isInDomain(p1) || !isInDomain(p2))
	return true; // false; VSK, 0607. Don't mess around with this link

    // Check if the points are within the tolerance
    shared_ptr<GeoTol> epsge = p1->getTolerance();
    double tol = epsge->getEpsge();
    if (p1->getDist() > tol || p2->getDist() > tol)
	return false;

    // Check if it's a microlink
    Point diff = p2->getPoint() - p1->getPoint();
    if (diff.length() < tol)
	return false;

    // If one of the endpoints is an isolated singularity, we do not accept the link
    if (p1->getSingularityType() == ISOLATED_POINT || p2->getSingularityType() == ISOLATED_POINT)
	return false;

    // If one of the endpoints is a branch point, we accept the link
    if (p1->getSingularityType() == BRANCH_POINT || p1->getSingularityType() == HIGHER_ORDER_POINT)
	return true;
    if (p2->getSingularityType() == BRANCH_POINT || p2->getSingularityType() == HIGHER_ORDER_POINT)
	return true;

    // Check if the tangents at the ends are consistent
    Point tangent1 = p1->getTangent();
    Point tangent2 = p2->getTangent();
    double cos = tangent1 * tangent2;
    bool passed_false = false;
    if (cos < 0.0)
    { 
	std::cout << "Inconsistent tangent direction " << std::endl;
	std::cout << "Par 1 : " << p1->getPar(0) << " " << p1->getPar(1) << " ";
	std::cout << p1->getPar(2) << " " << p1->getPar(3) << std::endl;
	std::cout << "Par 2 : " << p2->getPar(0) << " " << p2->getPar(1) << " ";
	std::cout << p2->getPar(2) << " " << p2->getPar(3) << std::endl;
	passed_false = true;
	//return false;
    }
    diff.normalize();
    if (fabs(diff * tangent1) < 0.25 || fabs(diff * tangent2) < 0.25)
    { 
	// cos = 0.25 is roughly 75 degrees
	std::cout << "Tangent inconsistent with vector " << std::endl;
	std::cout << "Par 1 : " << p1->getPar(0) << " " << p1->getPar(1) << " ";
	std::cout << p1->getPar(2) << " " << p1->getPar(3) << std::endl;
	std::cout << "Par 2 : " << p2->getPar(0) << " " << p2->getPar(1) << " ";
	std::cout << p2->getPar(2) << " " << p2->getPar(3) << std::endl;
	passed_false = true;
	//return false;
    }
    if (cos > 0.97)
	// cos = 0.97 is roughly 15 degrees
	return true;

    // Test for simple case
    vector<double> par1 = p1->getPar();
    vector<double> par2 = p2->getPar();
    // Get the parameter boxes for subdivision
    double frompar[4];
    double topar[4];
    double fuzzy = DEFAULT_PARAMETER_EPSILON; //epsge->getRelParRes();
    for (int k = 0; k < 4; ++k) {
	frompar[k] = (par1[k] <= par2[k]) ? par1[k] : par2[k];
	topar[k] = (par1[k] > par2[k]) ? par1[k] : par2[k];
	double start = startParam(k);
	double end = endParam(k);
	if (frompar[k] > topar[k] - fuzzy) {
	    frompar[k] -= 10.0 * fuzzy;
	    if (frompar[k] < start) {
		frompar[k] = start;
	    }
	    topar[k] += 10.0 * fuzzy;
	    if (topar[k] > end) {
		topar[k] = end;
	    }
	}
    }
    // Make subsurfaces
    vector<shared_ptr<ParamSurfaceInt> > sub1
	= obj1->subSurfaces(frompar[0], frompar[1],
			    topar[0], topar[1], fuzzy);
    vector<shared_ptr<ParamSurfaceInt> > sub2
	= obj2->subSurfaces(frompar[2], frompar[3],
			    topar[2], topar[3], fuzzy);
    // Create intersector and implicitize
    shared_ptr<SfSfIntersector> prev_dummy
	= shared_ptr<SfSfIntersector>
	(new SfSfIntersector(sub1[0], sub2[0], epsge));
    shared_ptr<SfSfIntersector> intersector
	= shared_ptr<SfSfIntersector>
	(new SfSfIntersector(sub1[0], sub2[0], epsge, prev_dummy.get()));
    int simple_case = intersector->simpleCase();
    if (simple_case == 1)
    {
	if (passed_false)
	    std::cout << "Simple case by cone test" << std::endl;
	return true;
    }

    double implicit_tol = 1.0e-12;
//     int recursion_limit = 0;
    //int recursion_limit = 10;
    intersector->createApproxImplicit(implicit_tol, recursion_limit);
    // Is this a simple case?
    simple_case = intersector->complexSimpleCase();
    if (simple_case == 1)
    {
	if (passed_false)
	    std::cout << "Simple case by implicitization" << std::endl;
	return true;
    }

//     // Last chance - run compute() and see if there is a path of
//     // connected points
//     intersector->compute();
//     shared_ptr<IntersectionPool> pool = intersector->getIntPool();
//     shared_ptr<IntersectionPoint> pnt1, pnt2;
//     bool found_closest;
//     found_closest
// 	= pool->closestInDomain(&(p1->getPar())[0], pnt1);
//     if (!found_closest)
// 	return false;
//     if (!p1->isSamePoint(pnt1.get()))
// 	return false;
//     found_closest
// 	= pool->closestInDomain(&(p2->getPar())[0], pnt2);
//     if (!found_closest)
// 	return false;
//     if (!p2->isSamePoint(pnt2.get()))
// 	return false;

//     if (pool->isConnectedInside(pnt1, pnt2))
// 	return true;

    return false;
}


//===========================================================================
bool IntersectionPool::validate() const
//===========================================================================
{
    // Check if the intersection points in the present pool also
    // exists in the parent pool.

    // Get intersection points from parent pool
    vector<shared_ptr<IntersectionPoint> > orig_int_pts;
    getOrigPoints(orig_int_pts);

    // Get intersection points from this pool
    vector<shared_ptr<IntersectionPoint> > int_pts
	= getIntersectionPoints();

    // Test if the points in this pool also exist in the parent pool
    int npts = int(int_pts.size());
    for (int i = 0; i < npts; ++i) {
	vector<shared_ptr<IntersectionPoint> >::iterator it
	    = find(orig_int_pts.begin(), orig_int_pts.end(), int_pts[i]);
	if (it == orig_int_pts.end()) {
	    prev_pool_->writeIntersectionPoints();
	    writeIntersectionPoints();
	    return false;
	}
    }
    return true;
}


//===========================================================================
void IntersectionPool::writeIntersectionPoints() const
//===========================================================================
{
    // Get intersection points from this pool
    vector<shared_ptr<IntersectionPoint> > int_pts
	= getIntersectionPoints();

    // Write points
    int npts = (int)int_pts.size();
    cout << "Number of points in pool " << this << ": "
	 << npts << endl;
    for (int i = 0; i < npts; ++i) {
	cout << i << "\t";
	int_pts[i]->writeInfo();
	cout << endl;
    }

    return;
}


//===========================================================================
void IntersectionPool::writeIntersectionLinks() const
//===========================================================================
{
    int numPar1 = obj1_->numParams();
    int numPar2 = obj2_->numParams();
    if (numPar1 + numPar2 < 3)
	return;

    // Get intersection points from this pool
    vector<shared_ptr<IntersectionPoint> > int_pts
	= getIntersectionPoints();
    int npts = (int)int_pts.size();

    // Compile vector of links
    vector<shared_ptr<IntersectionLink> > all_links;
    vector<shared_ptr<IntersectionLink> > links;
    for (int i = 0; i < npts; ++i) {
	int_pts[i]->getNeighbourLinks(links);
	all_links.insert(all_links.end(), links.begin(), links.end());
    }
    int nlinks = (int)all_links.size();

    vector<bool> is_ok(nlinks);
    for (int i = 0; i < nlinks; ++i) {
	is_ok[i] = verifyIntersectionLink(all_links[i]);
    }

    // Write the links
    cout << "Number of links in pool " << this << ": "
	 << nlinks << endl;
    for (int i = 0; i < nlinks; ++i) {
	cout << i << "\t";
	all_links[i]->writeInfo();
	cout << "  OK = " << is_ok[i] << endl;
    }

    return;
}


//===========================================================================
void IntersectionPool::
set_correct_differentiating_domain(shared_ptr<IntersectionPoint> ip) const
//===========================================================================
{
    int n1 = obj1_->numParams();
    int n2 = obj2_->numParams();
    int num_par = n1 + n2;
    double min_val, max_val;
    vector<bool> deriv_at_left(num_par, false);

    for (int i = 0; i < num_par; ++i) {
	get_param_limits(i, min_val, max_val);
	double tol = ip->parameterTolerance(i);
	if (fabs(ip->getPar(i) - max_val) < tol) {
	    deriv_at_left[i] = true; 
	}
    }

    ip->setDifferentiateFromLeft(deriv_at_left.begin()); 
}


//===========================================================================
void IntersectionPool::
determine_free_dir_parameter(double* par, 
			     shared_ptr<IntersectionLink> link,
			     int free_dir, 
			     int fixed_dir,
			     double fx_val)
//===========================================================================
{
    if (free_dir < 0) {	// unspecified free parameter, nothing to do
	return;
    }
    IntersectionPoint *p1, *p2;
    link->getIntersectionPoints(p1, p2);
    if (link->isIsoparametricIn(free_dir)) {
	par[free_dir] = p1->getPar(free_dir);
	return;
    }
    // if we got here, free_dir is specified, which means that the
    // object containing the fixed_dir must be a surface.  Moreover,
    // free_dir is _NOT_ isoparametric, which means that the other
    // object must have at least one isoparametric direction.  We must
    // calculate a reasonable value for the free_dir parameter.
    const ParamObjectInt* fixed_obj = p1->getObj1(); //obj1_.get(); // p1->getObj1();
    const ParamObjectInt* other_obj = p1->getObj2(); //obj2_.get(); // p1->getObj2();
    const bool first_obj_is_fixed
	= fixed_dir < p1->getObj1()->numParams();
    const int first_index_other
	= first_obj_is_fixed ? fixed_obj->numParams() : 0;
    if (!first_obj_is_fixed) {
	swap(fixed_obj, other_obj);
    }
    ASSERT(fixed_obj->numParams() == 2); // free_dir is specified,
					     // so this must be a
					     // surface

    shared_ptr<const ParamSurface> fix_psurf;
    const ParamSurfaceInt* fix_obj_surf
	= dynamic_cast<const ParamSurfaceInt*>(fixed_obj);
    if (fix_obj_surf) {
	fix_psurf = fix_obj_surf->getParamSurface();
    } else {
	const Param2FunctionInt* fix_obj2_surf
	    = dynamic_cast<const Param2FunctionInt*>(fixed_obj);
	if (fix_obj2_surf) {
	    fix_psurf = fix_obj2_surf->getParamSurface();
	} else {
	    THROW("Unexpected object!");
	}
    }

    shared_ptr<const CurveOnSurface> fc;
    if (fixed_dir == 0 || fixed_dir == p1->getObj1()->numParams()) {
	double pmin = first_obj_is_fixed
	    ? obj1_->startParam(1) : obj2_->startParam(1);
	double pmax = first_obj_is_fixed
	    ? obj1_->endParam(1) : obj2_->endParam(1);

	fc = make_curve_on_surface(fix_psurf, fx_val, pmin, fx_val, pmax,
				   pmin, pmax);
    } else {
	double pmin = first_obj_is_fixed
	    ? obj1_->startParam(0) : obj2_->startParam(0);
	double pmax = first_obj_is_fixed
	    ? obj1_->endParam(0): obj2_->endParam(0);

	fc = make_curve_on_surface(fix_psurf, pmin, fx_val, pmax, fx_val,
				   pmin, pmax);
    }

//     RectDomain d = fix_psurf->containingDomain();
//     shared_ptr<const CurveOnSurface> fc = 
// 	(fixed_dir == 0 || fixed_dir == p1->getObj1()->numParams()) ?  
// 	make_curve_on_surface(fix_psurf, fx_val, d.vmin(), fx_val, d.vmax(), d.vmin(), d.vmax()):
// 	make_curve_on_surface(fix_psurf, d.umin(), fx_val, d.umax(), fx_val, d.umin(), d.umax());
    
    // calculating seed
    const double gamma = (p1->getPar(fixed_dir) - fx_val) /
	                 (p1->getPar(fixed_dir) - p2->getPar(fixed_dir));
    const double seed = (1 - gamma) * p1->getPar(free_dir)
	+ gamma * p2->getPar(free_dir);
    Point pt1, pt2;
    double dist, min_par_other, max_par_other;
    shared_ptr<const ParamCurve> other_curve;
    switch(other_obj->numParams()) {
    case 1: // other object is a curve 
	{
	    const ParamCurveInt* pci
		= dynamic_cast<const ParamCurveInt*>(other_obj);
	    ASSERT(pci);
	    other_curve = pci->getParamCurve();
	    min_par_other = first_obj_is_fixed
		? obj2_->startParam(0) : obj1_->startParam(0);
	    max_par_other = first_obj_is_fixed
		? obj2_->endParam(0) : obj1_->endParam(0);
	}
	break;
    case 2: // other object is a surface
	{
	    const ParamSurfaceInt* psi
		= dynamic_cast<const ParamSurfaceInt*>(other_obj);
	    ASSERT(psi); 
	    shared_ptr<const ParamSurface> other_psurf
		= psi->getParamSurface();

	    if (link->isIsoparametricIn(first_index_other)) {
		min_par_other = first_obj_is_fixed
		    ? obj2_->startParam(1) : obj1_->startParam(1);
		max_par_other = first_obj_is_fixed
		    ? obj2_->endParam(1) : obj1_->endParam(1);

		other_curve
		    = make_curve_on_surface(other_psurf, 
					    p1->getPar(first_index_other),
					    min_par_other, //  d.vmin(), 
					    p1->getPar(first_index_other),
					    max_par_other, //d.vmax(), 
					    min_par_other, //d.vmin(), 
					    max_par_other); //d.vmax());
	    } else {
		ASSERT(link->isIsoparametricIn(first_index_other + 1));
		min_par_other = first_obj_is_fixed
		    ? obj2_->startParam(0) : obj1_->startParam(0);
		max_par_other = first_obj_is_fixed
		    ? obj2_->endParam(0) : obj1_->endParam(0);

		other_curve
		    = make_curve_on_surface(other_psurf,
					    min_par_other, //d.umin(),
					    p1->getPar(first_index_other+1),
					    max_par_other, //d.umax(),
					    p1->getPar(first_index_other+1),
					    min_par_other, //d.umin(),
					    max_par_other); //d.umax());
	    }
	}
	break;
    default:
	THROW("Logial error in determine_free_dir_parameter.  "
	      "Isoparametric object was neither a surface nor a curve.");
    }
    double par2;
    ClosestPoint::closestPtCurves(fc.get(), 
		    other_curve.get(),
		    fc->startparam(),
		    fc->endparam(),
		    min_par_other, // other_curve->startparam(),
		    max_par_other, //other_curve->endparam(),
		    seed, 
		    0.5 * (min_par_other + max_par_other),
		    par[free_dir],
		    par2,
		    dist,
		    pt1,
		    pt2);
}    


//===========================================================================
void IntersectionPool::writeDebug(int singular)
//===========================================================================
{
    //ASSERT(validate());

    // Set the length of tangent vectors in 2Ddisp. This shouldn't
    // really be hardcoded, but it's much easyer. To turn off
    // tangents, use a negative number.
    const double tangent_scale_fac = 1e-3;
//     const double tangent_scale_fac = -1.0;

    vector<shared_ptr<IntersectionPoint> > orig_points;
    getOrigPoints(orig_points);

    int numPar1 = obj1_->numParams();
    int numPar2 = obj2_->numParams();

    // For the time being
    /*if (numPar1 + numPar2 < 4)
      return;*/

    // Surface-surface
    if (getenv("TOTAL_POINT") && *(getenv("TOTAL_POINT")) == '1')
    {
    
    if (numPar1 == 2 && numPar2 == 2) {
	std::cout << "Num int points in total: "
		  << orig_points.size() << std::endl;
	for (int ki = 0; ki < int(orig_points.size()); ki++) {
	    orig_points[ki]->writeParams(cout);
	    cout << endl;
	}
    }
    }

    // First object is a surface
    if (numPar1 == 2) {

	ofstream debug((singular==1) ? "sing_domain1.dsp" : "domain1.dsp");

	// Write tangent vectors first - they are then in the
	// "background"
	if (tangent_scale_fac > 0.0) {
	    debug << "pnttype: point" << endl
		  << "fg: black" << endl;
	    for (int i = 0; i < int(orig_points.size()); ++i) {
		int nbranches = orig_points[i]->numBranches();
		double fromu = orig_points[i]->getPar1()[0];
		double fromv = orig_points[i]->getPar1()[1];
		for (int j = 1; j <= nbranches; ++j) {
		    bool second_branch = (j == 2 ? true : false);
		    Point tangent;
		    try {
			tangent = orig_points[i]->getPar1Dir(second_branch);
		    } catch (...)
		    {
			continue;
		    }

		    // @@@jbt - The following is a hack to get rid of
		    // tangents with nan values
		    if (!(tangent[0] == tangent[0]))
			break;
		    if (!(tangent[1] == tangent[1]))
			break;
		    try
		    {
		    tangent.normalize();
		    } catch (...)
		    {
			continue;
		    }
		    tangent *= tangent_scale_fac;
		    double tou = fromu + tangent[0];
		    double tov = fromv + tangent[1];
		    debug << "lin:" << endl
			  << fromu << " " << fromv << endl
			  << tou << " " << tov << endl
			  << "pnt:" << endl
			  << tou << " " << tov << endl;
		}
	    }
	}

	// Write brown box to indicate current parameter domain
	double startu = obj1_->startParam(0);
	double endu = obj1_->endParam(0);
	double startv = obj1_->startParam(1);
	double endv = obj1_->endParam(1);
	debug << "fg: brown" << endl;
	debug << "lin:" << endl;
	debug << startu << " " << startv << endl;
	debug << endu << " " << startv << endl;
	debug << endu << " " << endv << endl;
	debug << startu << " " << endv << endl;
	debug << startu << " " << startv << endl;

	// Write intersection points and links to neighbours
	debug << "pnttype: star" << endl;
	for (int ki = 0; ki < int(orig_points.size()); ki++) {
	    double paru = orig_points[ki]->getPar1()[0];
	    double parv = orig_points[ki]->getPar1()[1];
	    bool in_this_pool = find(int_points_.begin(), int_points_.end(),
				     orig_points[ki]) != int_points_.end();

	    // Points
	    if (singular == 1) {
		debug << "fg: black" << endl;
	    }
	    else if (in_this_pool) {
		if (orig_points[ki]->pointIsSingular()) {
		    debug << "fg: cyan" << endl;
		}
		else {
		    debug << "fg: red" << endl;
		}
	    }
	    else {
		debug << "fg: green" << endl;
	    }
	    debug << "pnt:" << endl;
	    debug << paru << " " << parv << endl;

	    // Links
	    vector<IntersectionPoint*> neighbours;
	    orig_points[ki]->getNeighbours(neighbours);
	    for (int kj = 0; kj < int(neighbours.size()); kj++) {
		double neighu = neighbours[kj]->getPar1()[0];
		double neighv = neighbours[kj]->getPar1()[1];

		shared_ptr<IntersectionLink> link
		    = orig_points[ki]->getIntersectionLink(neighbours[kj]);
		if (singular == 1) {
		    debug << "fg: black" << endl;
		} else {
		    if (link->isIsoparametricIn(0)
			|| link->isIsoparametricIn(1)) {
			debug << "fg: magenta" << endl;
		    }
		    else {
			debug << "fg: blue" << endl;
		    }
		    if (link->linkType() == COINCIDENCE_SFCV
			|| link->linkType() == MERGED_COINCIDENCE_SFCV_SFCV) {
			debug << "fg: green" << endl;
		    }
		}
		debug << "lin:" << endl;
	        debug << paru << " " << parv << endl;
		debug << neighu << " " << neighv << endl;    
	    }

	}

    }
    
    // Second object is a surface
    if (numPar2 == 2) {

	ofstream debug((singular==1) ? "sing_domain2.dsp" : "domain2.dsp");

	// Write tangent vectors first - they are then in the
	// "background"
	if (tangent_scale_fac > 0.0) {
	    debug << "pnttype: point" << endl
		  << "fg: black" << endl;
	    for (int i = 0; i < int(orig_points.size()); ++i) {
		int nbranches = orig_points[i]->numBranches();
		double fromu = orig_points[i]->getPar2()[0];
		double fromv = orig_points[i]->getPar2()[1];
		for (int j = 1; j <= nbranches; ++j) {
		    bool second_branch = (j == 2 ? true : false);
		    Point tangent;
		    try {
			tangent = orig_points[i]->getPar2Dir(second_branch);
		    }
		    catch (...)
		    {
			continue;
		    }

		    // @@@jbt - The following is a hack to get rid of
		    // tangents with nan values
		    if (!(tangent[0] == tangent[0]))
			break;
		    if (!(tangent[1] == tangent[1]))
			break;
		    try
		    {
		    tangent.normalize();
		    } catch (...)
		    {
			continue;
		    }
		    tangent *= tangent_scale_fac;
		    double tou = fromu + tangent[0];
		    double tov = fromv + tangent[1];
		    debug << "lin:" << endl
			  << fromu << " " << fromv << endl
			  << tou << " " << tov << endl
			  << "pnt:" << endl
			  << tou << " " << tov << endl;
		}
	    }
	}

	// Write brown box to indicate current parameter domain
	double startu = obj2_->startParam(0);
	double endu = obj2_->endParam(0);
	double startv = obj2_->startParam(1);
	double endv = obj2_->endParam(1);
	debug << "fg: brown" << endl;
	debug << "lin:" << endl;
	debug << startu << " " << startv << endl;
	debug << endu << " " << startv << endl;
	debug << endu << " " << endv << endl;
	debug << startu << " " << endv << endl;
	debug << startu << " " << startv << endl;

	// Write intersection points and links to neighbours
	debug << "pnttype: star" << endl;
	for (int ki = 0; ki < int(orig_points.size()); ki++) {
	    double paru = orig_points[ki]->getPar2()[0];
	    double parv = orig_points[ki]->getPar2()[1];
	    bool in_this_pool = find(int_points_.begin(), int_points_.end(),
				     orig_points[ki]) != int_points_.end();

	    // Points
	    if (singular == 1) {
		debug << "fg: black" << endl;
	    }
	    else if (in_this_pool) {
		if (orig_points[ki]->pointIsSingular()) {
		    debug << "fg: cyan" << endl;
		}
		else {
		    debug << "fg: red" << endl;
		}
	    }
	    else {
		debug << "fg: green" << endl;
	    }
	    debug << "pnt:" << endl;
	    debug << paru << " " << parv << endl;

	    // Links
	    vector<IntersectionPoint*> neighbours;
	    orig_points[ki]->getNeighbours(neighbours);
	    for (int kj = 0; kj < int(neighbours.size()); kj++) {
		double neighu = neighbours[kj]->getPar2()[0];
		double neighv = neighbours[kj]->getPar2()[1];

		shared_ptr<IntersectionLink> link
		    = orig_points[ki]->getIntersectionLink(neighbours[kj]);
		if (singular == 1) {
		    debug << "fg: black" << endl;
		} else {
		    if (link->isIsoparametricIn(numPar1)
			|| link->isIsoparametricIn(numPar1 + 1)) {
			debug << "fg: magenta" << endl;
		    }
		    else {
			debug << "fg: blue" << endl;
		    }
		    if (link->linkType() == COINCIDENCE_SFCV
			|| link->linkType() == MERGED_COINCIDENCE_SFCV_SFCV) {
			debug << "fg: green" << endl;
		    }
		}
		debug << "lin:" << endl;
	        debug << paru << " " << parv << endl;
		debug << neighu << " " << neighv << endl;    
	    }
	}

    }

    return;
}


}; // namespace Go


// //===========================================================================
// void IntersectionPool::addPoints(int nmb_int_pts, double* pointpar1, 
// 				 double* pointpar2,
// 				 SISLSurf* sisl_sf1, SISLSurf* sisl_sf2)
// //===========================================================================
// {
//   const int dim1 = sisl_sf1->idim;
//   const int dim2 = sisl_sf2->idim;

//   int stat;
//   int der1 = 0;
//   int der2 = 0;
//   int il1_1 = 0;
//   int il2_1 = 0;
//   int il1_2 = 0;
//   int il2_2 = 0;
//   vector<double> derive1((der1+1)*dim1);  // @bsp Bare OK for der1=0
//   vector<double> derive2((der2+1)*dim2);
//   int_points_.reserve(nmb_int_pts);
 
//   for (int i=0; i<nmb_int_pts; ++i) {
//       //s1424(sisl_sf1, der1, der2, pointpar1, &il1_1, &il2_1, &derive1[0], &stat);
//       s1424(sisl_sf1, der1, der2, pointpar1 + (2 * i), &il1_1, &il2_1, &derive1[0], &stat);
//       Point point1(derive1[0],derive1[1],derive1[2]);
 
//       //s1424(sisl_sf2, der1, der2, pointpar2, &il1_2, &il2_2, &derive2[0], &stat);
//       s1424(sisl_sf2, der1, der2, pointpar2 + (2 * i), &il1_2, &il2_2, &derive2[0], &stat);
//       Point point2(derive2[0],derive2[1],derive2[2]);
      
//       IntersectionPoint int_point(obj1_,obj2_ ,sisl_sf1, sisl_sf2,
// 				  point1, point2, pointpar1, pointpar2);
//       int_points_.push_back(int_point);
//   }

// }


// //===========================================================================
// void IntersectionPool::addCurves(int nmb_int_cvs, SISLIntcurve** intcurves,
// 				 SISLSurf* sisl_sf1, SISLSurf* sisl_sf2)
// //===========================================================================
// {

//   int_curves_.reserve(nmb_int_cvs);
 
//   const int draw = 0;
//   const int makecv = 2; // Make both geometric and parametric curves.
//   int status = 0;
//   const double march_eps = 0.01;
//   const double maxstep = 0.0;

//   for (int i = 0; i < nmb_int_cvs; ++i) {
//     s1310(sisl_sf1, sisl_sf2, intcurves[i], march_eps, maxstep, makecv, draw, &status);
//     ALWAYS_ERROR_IF(status < 0,
// 		"Failed intersecting surfs.", UnknownError());
//     MESSAGE_IF(status != 0, "Returned status value: " << status);

//     IntersectionCurve int_curve(obj1_,obj2_, intcurves[i], sisl_sf1,sisl_sf2);
//     int_curves_.push_back(int_curve);
//   }
	

  
// }




// //===========================================================================
// void
// IntersectionPool::includeReducedInts(shared_ptr<IntersectionPool> lower_order_pool)
// //===========================================================================
// {
//     int lacking_ix = lower_order_pool->lackingParameter();
//     double lacking_val = lower_order_pool->lackingParameterValue();
//     DEBUG_ERROR_IF(lacking_ix < 0, 
// 	     "Argument given to includeReducedInts was not a reduced pool");

//     const vector<shared_ptr<IntersectionPoint> >& 
// 	ipoints = lower_order_pool->getIntersectionPoints();

//     for (int other = 0; other < int(ipoints.size()); ++other) {
// 	const shared_ptr<IntersectionPoint> 
// 	    other_parent_ip = ipoints[other]->parentPoint(); /// could be 0
// 	bool found = false;

// 	// checking if this IntersectionPoint is already represented in 'this' pool
// 	for (int local = 0; local < int(int_points_.size()) && !found; ++local) {
// 	    const shared_ptr<IntersectionPoint> this_ip = int_points_[local];
// 	    if (this_ip == other_parent_ip) {
// 		found = true;
// 	    }
// 	}
// 	if (!found) {
// 	    // This point does not already exist in the current pool.  Let us 
// 	    // add it!
// 	    ASSERT(other_parent_ip.get() == 0); // @ just a check.  Assumed to be true.

//  	    // preparing temporary parameter array
// 	    double par1[2]; 
// 	    int i, cpos;
// 	    for (i = 0, cpos = 0; i < obj1_->numParams(); ++i) {
// 		if (i == lacking_ix) {
// 		    par1[i] = lacking_val;
// 		} else {
// 		    par1[i] = ipoints[other]->getPar1()[cpos++];
// 		}
// 	    }
// 	    int lacking_ix2 = lacking_ix - obj1_->numParams();
// 	    double par2[2];
// 	    for (i = 0, cpos = 0; i < obj2_->numParams(); ++i) {
// 		if (i == lacking_ix2) {
// 		    par2[i] = lacking_val;
// 		} else {
// 		    par2[i] = ipoints[other]->getPar2()[cpos++];
// 		}
// 	    }
// 	    shared_ptr<IntersectionPoint> 
// 		temp(new IntersectionPoint(obj1_, 
// 					   obj2_, 
// 					   ipoints[other]->getTolerance(),
// 					   par1, 
// 					   par2));
// 	    add_point_and_propagate_upwards(temp);
// 	}
//     }
// }


//===========================================================================
