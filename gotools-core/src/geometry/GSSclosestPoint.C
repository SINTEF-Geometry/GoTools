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

#include <algorithm>
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/Utils.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::max;
using std::min;

namespace {
//===========================================================================
// squared distance function between a point and a surface.  Used by the 
// minimization algorithm initiated by SplineSurface::closestPoint
class PtSfDist2 {
//===========================================================================
public:
    PtSfDist2(const Point& pt, 
	      const SplineSurface& sf,
	      const RectDomain* rd) 
	: pt_(pt), sf_(sf), tmp_ptvec_(3) {
	if (rd) {
	    ll_[0] = rd->umin(); ll_[1] = rd->vmin();	    
	    ur_[0] = rd->umax(); ur_[1] = rd->vmax();
	} else {
	    ll_[0] = sf_.startparam_u(); ll_[1] = sf_.startparam_v();
	    ur_[0] = sf_.endparam_u();   ur_[1] = sf_.endparam_v();
	}
    }
    inline double operator()(const double* arg) const;
    inline void grad(const double* arg, double* res) const;
    inline double minPar(int pardir) const { 
	return ll_[pardir];
	//return (pardir == 0) ? sf_.startparam_u() : sf_.startparam_v();
    }
    inline double maxPar(int pardir) const {
	return ur_[pardir];
	//return (pardir == 0) ? sf_.endparam_u() : sf_.endparam_v();
    }
private:
    const Point& pt_;
    const SplineSurface sf_;
    double ll_[2]; // lower left corner of domain
    double ur_[2]; // upper right corner of domain
    mutable Point tmp_pt_;
    mutable vector<Point> tmp_ptvec_;
};

//===========================================================================
double PtSfDist2::operator()(const double* arg) const
//===========================================================================
{
    sf_.point(tmp_pt_, arg[0], arg[1]);
    return pt_.dist2(tmp_pt_);
}

//===========================================================================
void PtSfDist2::grad(const double* arg, double* res) const
//===========================================================================
{
    sf_.point(tmp_ptvec_, arg[0], arg[1], 1);
    tmp_pt_ = tmp_ptvec_[0] - pt_; // distance vector from point to surface
    res[0] = 2 * tmp_ptvec_[1] * tmp_pt_;
    res[1] = 2 * tmp_ptvec_[2] * tmp_pt_;
    //return pt_.dist2(tmp_ptvec_[0]); 
}

//===========================================================================
// find a good seed for closest point computation
void robust_seedfind(const Point& pt, 
		     const SplineSurface& sf, 
		     const RectDomain* rd,
		     double& u,
		     double& v)
//===========================================================================
{
    // Finding the closest triangle in a triangulation of the control grid.
    // We improve the performance by searching on triangulation restricted
    // to coefficients as given by 'rd' (if it exists)
    int min_ind_u = 0;
    int max_ind_u = sf.numCoefs_u() - 1;
    int min_ind_v = 0;
    int max_ind_v = sf.numCoefs_v() - 1;
    if (rd != NULL) {
	int dummy;
	sf.basis_u().coefsAffectingParam(rd->umin(), min_ind_u, dummy);
	sf.basis_u().coefsAffectingParam(rd->umax(), dummy, max_ind_u);
	sf.basis_v().coefsAffectingParam(rd->vmin(), min_ind_v, dummy);
	sf.basis_v().coefsAffectingParam(rd->vmax(), dummy, max_ind_v);
    }
    double start_u = 0.0;
    double start_v = 0.0;
    vector<double>::const_iterator coefs = sf.coefs_begin();
    SplineUtils::closest_on_rectgrid(pt.begin(), &coefs[0],
			min_ind_u, max_ind_u,
			min_ind_v, max_ind_v,
			sf.numCoefs_u(),
			start_u, start_v);
    
    // The returned u and v parameters know nothing about the
    // knot vectors and such.
    int u_ind = int(floor(start_u));
    int v_ind = int(floor(start_v));
    u_ind = std::max(0, std::min(u_ind, sf.numCoefs_u()));
    v_ind = std::max(0, std::min(v_ind, sf.numCoefs_v()));
    double u1 = sf.basis_u().grevilleParameter(u_ind);
    double v1 = sf.basis_v().grevilleParameter(v_ind);
    if (u_ind == sf.numCoefs_u() - 1) {
	u = sf.basis_u().grevilleParameter(u_ind);
    } else {
	double u2 = sf.basis_u().grevilleParameter(u_ind + 1);
	double fac = start_u - double(u_ind);
	u = u1*(1.0-fac) + fac*u2;
    }
    if (v_ind == sf.numCoefs_v() - 1) {
	v = sf.basis_v().grevilleParameter(v_ind);
    } else {
	double v2 = sf.basis_v().grevilleParameter(v_ind + 1);
	double fac = start_v - double(v_ind);
	v = v1*(1.0-fac) + fac*v2;
    }

    // Ensure that the seed is really inside the given domain
    if (rd)
      {
	u = std::max(rd->umin(), std::min(u, rd->umax()));
	v = std::max(rd->vmin(), std::min(v, rd->vmax()));
      }
}

}; // end anonymous namespace 


namespace Go {

//===========================================================================
void SplineSurface::closestPoint(const Point& pt,
				 double& clo_u,
				 double& clo_v, 
				 Point& clo_pt,
				 double& clo_dist,
				 double epsilon,
				 const RectDomain* rd,
				 double *seed) const
//===========================================================================
{
    // VSK, 0611. The conjugate gradient method is much slower than
    // the closest point iterations fetched from SISL, but it seems to
    // be more stable in some tangential cases. We need a compromise!!!
  static bool use_conjugate_gradient = (iterator_ == Iterator_parametric) ? true : false;
  //static bool use_conjugate_gradient = false;
    
    double seed_buf[2];
    if (!seed) {
	// no seed given, we must compute one
	seed = seed_buf;
	robust_seedfind(pt, *this, rd, seed[0], seed[1]);
    }

    bool at_bd = false;
    clo_dist = -1.0;
    if (use_conjugate_gradient)
    {
	PtSfDist2 dist_fun(pt, *this, rd);
	// define distance function
	FunctionMinimizer<PtSfDist2> funmin(2, dist_fun, seed, epsilon);
    
	// minimise distance function
	try {
	minimise_conjugated_gradient(funmin); //, 3) ; // number of iterations in each cycle
	}
	catch (...)
	  {
	    use_conjugate_gradient = false;
	  }

	if (use_conjugate_gradient)
	  {
	    // return results
	    clo_u = funmin.getPar(0);
	    clo_v = funmin.getPar(1);
	    clo_dist = sqrt(funmin.fval());
	    point(clo_pt, clo_u, clo_v);
            // @@sbr201710 The conjugate gradient method seems to be unstable at the boundary. We should look into this.
            // Current test case where this happens involves a rational surface (Kaplan_blade_foundry_model.stp).
            at_bd = (funmin.atMin(0) || funmin.atMax(0) || funmin.atMin(1) || funmin.atMax(1));
	  }
    }

    if ((!use_conjugate_gradient) || (at_bd))
    {
	int kstat = 0;
	double par[2], start[2], end[2];
	start[0] = (rd) ? rd->umin() : startparam_u();
	start[1] = (rd) ? rd->vmin() : startparam_v();
	end[0] = (rd) ? rd->umax() : endparam_u();
	end[1] = (rd) ? rd->vmax() : endparam_v();
	s1773(pt.begin(), epsilon, start, end, seed, par, &kstat);
        Point clo_pt2 = ParamSurface::point(par[0], par[1]);
        double clo_dist2 = pt.dist(clo_pt2);
        if ((clo_dist < 0.0) || (clo_dist2 < clo_dist))
        {
            // MESSAGE_IF(at_bd,
            //            "Using sisl closest point result. clo_dist = " << clo_dis
		       // t << ", clo_dist2 = " << clo_dist2);
            clo_u = par[0];
            clo_v = par[1];
            clo_dist = clo_dist2;
            clo_pt = clo_pt2;
        }
    }
}

// ---- OLD CODE, BUT KEPT FOR FUTURE REFERENCE ----


// //===========================================================================
// void SplineSurface::closestPoint(const Point& pt,
// 				 double& clo_u,
// 				 double& clo_v, 
// 				 Point& clo_pt,
// 				 double& clo_dist,
// 				 double epsilon,
// 				 const RectDomain* rd,
// 				 double *seed) const
// //===========================================================================
// {
//     closestPointImpl(pt, clo_u, clo_v, clo_pt, clo_dist,
// 		     epsilon, rd, true, seed);
// }


// void SplineSurface::closestPointImpl(const Point& pt,
// 				     double& clo_u,
// 				     double& clo_v, 
// 				     Point& clo_pt,
// 				     double& clo_dist,
// 				     double epsilon,
// 				     const RectDomain* rd,
// 				     bool robust_seedfind,
// 				     double *seed) const
// /*****************************************************************************
//  *
//  * PURPOSE:	Given a Nurbsurface and a point it calculates the 
//  *              closest point on the surface.
//  * IMPORTANT:	It calculates only a local solution using a Newton iteration.
//  *		It makes use of the fact that 
//  *
//  *			(pt - S(u0,v0)) * Su(u0,v0) = 0                   (*)
//  *			(pt - S(u0,v0)) * Sv(u0,v0) = 0                   
//  *
//  *		if the closest parameter (u0,v0) lies in the interior of the 
//  *		surface domain.
//  *		In principle one could now do a Newton iteration for the
//  *		zeros of the function 
//  *		                  / (pt - S(u,v)) * Su(u,v) \            //
//  *		        f(u,v) =  |                         |
//  *		                  \ (pt - S(u,v)) * Sv(u,v) /
//  *		We do something else in order to avoid second derivatives:
//  *		We taylorize  (pt - S(u,v))  and get
//  *
//  *		     (pt - S(u,v) - Su(u,v)*du - Sv(u,v)*dv)*Su(u,v) = 0
//  *		     (pt - S(u,v) - Su(u,v)*du - Sv(u,v)*dv)*Sv(u,v) = 0
//  *
//  *		for a parameter (u,v) close to (u0,v0).
//  *
//  *              Solving the equation
//  *              
//  *                      / Su*Su Su*Sv \  / du \         / Su*(pt - S) \   //
//  *                      |             |  |    |    =    |             |
//  *                      \ Su*Sv Sv*Sv /  \ dv /         \ Sv*(pt - S) /
//  *              
//  *              ie. inverting the matrix
//  *              
//  *                      / Su*Su Su*Sv \-1           /  Sv*Sv -Su*Sv \    //
//  *                      |             |   =   1/det |               |  
//  *                      \ Su*Sv Sv*Sv /             \ -Su*Sv  Su*Su /  
//  *              
//  *              and applying it to the right hand side.
//  *              
//  *                    / du \             /  Sv*Sv -Su*Sv \   / Su*(pt - S) \  //
//  *                    |    |   =  1/det  |               |   |             |
//  *                    \ dv /             \ -Su*Sv  Su*Su /   \ Sv*(pt - S) /
//  *              
//  *
//  *			/ du \                                              //
//  *              Here    |    |  is the correction term for a Newton type
//  *                      \ dv /  
//  *		iteration.
//  *				
//  *
//  * WRITTEN BY:  Hermann Kellermann, SINTEF, Oslo, 1999-04-20
//  *
//  ****************************************************************************/
// {
// #ifdef DEBUG
//     if (seed != 0) {
// 	std::cout << "Seed is " << seed[0] << ' ' << seed[1] << std::endl;
//     }
// #endif

//   // Points where we want to find the closest point on a boundary is
//   // treated separately.
//     RectDomain domain = containingDomain();
//     if (!rd)
// 	rd = &domain;

//     std::vector<Point> derivs(3, Point(0.0, 0.0, 0.0));
//     Point& S = derivs[0];
//     Point& Su = derivs[1];
//     Point& Sv = derivs[2];

//     double umin = rd->umin();
//     double umax = rd->umax();
//     double vmin = rd->vmin();
//     double vmax = rd->vmax();
    
//     // ---------------- Treating the surface boundaries -----------------------
//     //@hke  Currently not done: We expect the point to lie not too far from
//     //@hke  the surface. So the iteration for the interior will deliver a 
//     //@hke  closest point which is not thaaaaaaaat bad.
//     //@hke  But: We might do something in the future.


//     // --------- searching in the interior of the domain ---------------

//     // finding a start parameter for the Newton iteration.
//     // We use the Greville parameter for the closest control point.
//     double u, v;
//     if (seed == 0) {
// 	if (!robust_seedfind) {
// 	    // Old implementation
// 	    int best_index = SplineUtils::closest_in_array(pt.begin(), &coefs_[0],
// 					      numCoefs_u()*numCoefs_v(), dim_);
// 	    int best_index_u = best_index%numCoefs_u(); 
// 	    int best_index_v = best_index/numCoefs_u();
// 	    u = basis_u().grevilleParameter(best_index_u);
// 	    v = basis_v().grevilleParameter(best_index_v);
// 	} else { // robust_seedfind == true
// 	    // New implementation, finding the closest triangle in a
// 	    // triangulation of the control grid.

// 	    // We improve the performance by searching on triangulation
// 	    // restricted to coefs as given by rect_domain (if it exists).
// 	    int min_ind_u = 0;
// 	    int max_ind_u = numCoefs_u() - 1;
// 	    int min_ind_v = 0;
// 	    int max_ind_v = numCoefs_v() - 1;
// 	    if (rd != NULL) {
// 		int dummy;
// 		basis_u().coefsAffectingParam(rd->umin(), min_ind_u, dummy);
// 		basis_u().coefsAffectingParam(rd->umax(), dummy, max_ind_u);
// 		basis_v().coefsAffectingParam(rd->vmin(), min_ind_v, dummy);
// 		basis_v().coefsAffectingParam(rd->vmax(), dummy, max_ind_v);
// 	    }
// 	    double start_u = 0.0;
// 	    double start_v = 0.0;
// 	    SplineUtils::closest_on_rectgrid(pt.begin(), &coefs_[0],
// 				min_ind_u, max_ind_u,
// 				min_ind_v, max_ind_v,
// 				numCoefs_u(),
// 				start_u, start_v);
// 	    // The returned u and v parameters know nothing about the
// 	    // knot vectors and such.
// 	    int u_ind = int(floor(start_u));
// 	    int v_ind = int(floor(start_v));
// 	    double u1 = basis_u().grevilleParameter(u_ind);
// 	    double v1 = basis_v().grevilleParameter(v_ind);
// 	    if (u_ind == numCoefs_u() - 1) {
// 		u = basis_u().grevilleParameter(u_ind);
// 	    } else {
// 		double u2 = basis_u().grevilleParameter(u_ind + 1);
// 		double fac = start_u - double(u_ind);
// 		u = u1*(1.0-fac) + fac*u2;
// 	    }
// 	    if (v_ind == numCoefs_v() - 1) {
// 		v = basis_v().grevilleParameter(v_ind);
// 	    } else {
// 		double v2 = basis_v().grevilleParameter(v_ind + 1);
// 		double fac = start_v - double(v_ind);
// 		v = v1*(1.0-fac) + fac*v2;
// 	    }
// 	}
//     } else { // seed != 0
// 	u = seed[0];
// 	v = seed[1];
//     }

//     // If outside the allowed rectangle in the parameter domain, move it
//     // to the edge of that area
//     u = std::max(u, umin);
//     u = std::min(u, umax);
//     v = std::max(v, vmin);
//     v = std::min(v, vmax);
//     Point closest(dimension());
//     point(closest, u, v);
//     double dist2 = closest.dist2(pt);
//     double mindist2 = dist2;
//     // Our best guesses so far
//     clo_u = u;
//     clo_v = v;
//     clo_pt = closest;
//     clo_dist = sqrt(dist2);


//     // not used ?
// //     double start_u = basis_u_.begin()[basis_u_.order() - 1];
// //     double end_u   = basis_u_.begin()[basis_u_.numCoefs()];
// //     double start_v = basis_v_.begin()[basis_v_.order() - 1];
// //     double end_v   = basis_v_.begin()[basis_v_.numCoefs()];


//     // the Newton iteration
//     double space_epsilon2 = epsilon * epsilon;
//     Point diff(dimension());
//     double space_delta2 = 1e100;

//     double u_old = u;
//     double v_old = v;
//     double u_first = u;
//     double v_first = v;

//     const int max_passes = 30;
//     const int max_boundary_crosses = 5;

//     int pass = 0;
//     int boundary_crossed = 0;
//     bool finished = false;

//     while (!finished) {

// 	/*
// 	 *  We correct the parameter according to the formula
// 	 *         (S + Su*du + Sv*dv - pt)*Su = 0
// 	 *         (S + Su*du + Sv*dv - pt)*Sv = 0
// 	 *  or 
// 	 *         / Su*Su Su*Sv \  / du \         / Su*(pt - S) \       //
// 	 *         |             |  |    |    =    |             |
// 	 *         \ Su*Sv Sv*Sv /  \ dv /         \ Sv*(pt - S) / 
// 	 */

// 	point(derivs, u, v, 1);

// 	// diff = pt - S
// 	diff = pt;
// 	diff -= S;

// 	dist2 = diff.length2();

// 	double SuSu = Su*Su;
// 	double SuSv = Su*Sv;
// 	double SvSv = Sv*Sv;

// 	double det = SuSu*SvSv - SuSv*SuSv;

// 	if (det < epsilon) {
// //  	    MESSAGE("Surface has degenerate tangent plane.\nSuSu = " 
// //  		       << SuSu << "\nSvSv = " << SvSv << "\nSuSv = " << SuSv);
// 	    // If the current point is closer that the closest boundary 
// 	    // point then we return this one.

// 	    if (dist2 < mindist2) {
// 		mindist2 = dist2;
// 		clo_u = u;
// 		clo_v = v;
// 		clo_pt = S;
// 		clo_dist = sqrt(dist2);
// 	    }
// 	    return;
// 	}


// 	/* Solving the equation
// 	 *
// 	 *         / Su*Su Su*Sv \  / du \         / Su*(pt - S) \       //
// 	 *         |             |  |    |    =    |             |
// 	 *         \ Su*Sv Sv*Sv /  \ dv /         \ Sv*(pt - S) /
// 	 *
// 	 * ie. inverting the matrix
// 	 *
// 	 *         / Su*Su Su*Sv \-1           /  Sv*Sv -Su*Sv \         //
// 	 *         |             |   =   1/det |               |  
// 	 *         \ Su*Sv Sv*Sv /             \ -Su*Sv  Su*Su /  
// 	 *	
// 	 * and applying it to the right hand side.
// 	 *
// 	 *       / du \             /  Sv*Sv -Su*Sv \   / Su*(pt - S) \     //
// 	 *       |    |   =  1/det  |               |   |             |
// 	 *       \ dv /             \ -Su*Sv  Su*Su /   \ Sv*(pt - S) /
// 	 */


// 	double SuDiff = Su*diff;
// 	double SvDiff = Sv*diff;

// 	double du =  (SvSv/det)*SuDiff - (SuSv/det)*SvDiff;
// 	double dv = -(SuSv/det)*SuDiff + (SuSu/det)*SvDiff;

// 	u += du;
// 	v += dv;

// 	// @afr 29032000: Staying inside the parameter domain
// 	// has been replaced by staying inside a rectangle
// 	// [umin,umax]x[vmin,vmax] in the parameter domain
// 	//
// 	// make sure that the corrected parameters still lies in the domain

// 	//        if (u < start_u)
// 	//  	 u = start_u;
// 	//        else if (u > end_u)
// 	//  	 u = end_u;

// 	//        if (v < start_v)
// 	//  	 v = start_v;
// 	//        else if (v > end_v)
// 	//  	 v = end_v;

// 	// If outside the allowed rectangle in the parameter domain, move it
// 	// to the edge of that area
// 	// We may have reason to expect the point to lie on the boundary.
// 	// 	if (boundary_pt) {
// 	// 	  if (min(fabs(u-umin), fabs(u - umax)) <
// 	// 	      min(fabs(v-vmin), fabs(v - vmax))) {
// 	// 	    if (fabs(u - umin) < (umax - umin) / 2.0)
// 	// 	      u = umin;
// 	// 	    else if (fabs(u - umax) < (umax - umin) / 2.0)
// 	// 	      u = umax;
// 	// 	  } else {
// 	// 	    if (fabs(v - vmin) < (vmax - vmin) / 2.0)
// 	// 	      v = vmin;
// 	// 	    else if (fabs(v - vmax) < (vmax - vmin) / 2.0)
// 	// 	      v = vmax;
// 	// 	  }
// 	// 	}

// 	if (u<umin || u>umax || v<vmin || v>vmax)
// 	    ++boundary_crossed;
// 	u = std::max(u, umin);
// 	u = std::min(u, umax);
// 	v = std::max(v, vmin);
// 	v = std::min(v, vmax);
// 	//	if (u==umin || u==umax || v==vmin || v==vmax)
// 	//	    ++boundary_crossed;
// 	//        else
// 	//  	  boundary_crossed = 0;

// 	// and this is a rough estimate of the correction in space. If it's
// 	// really small then we stop the iteration. The iteration stops in
// 	// particular if we are sitting on the curve boundary trying to
// 	// get out of the curve domain.

// 	space_delta2 = (Su*(u - u_old) - Sv*(v - v_old)).length2();

// 	u_old = u;
// 	v_old = v;

// 	++pass;
// 	if (space_delta2 < space_epsilon2
// 	    || pass >= max_passes
// 	    || boundary_crossed >= max_boundary_crosses)
// 	    finished = true;
//     }

//     if (pass >= max_passes) {
// 	MESSAGE("Closest point iteration did not converge.");
// 	Point firstpt;
// 	point(firstpt, u_first, v_first);
// 	if (seed && firstpt.dist2(pt) < dist2)
// 	    {
// 		u = u_first;
// 		v = v_first;
// 		S = firstpt;
// 		dist2 = firstpt.dist2(pt);
// 	    }
//     }


//     clo_u = u;
//     clo_v = v;
//     clo_pt = S;
//     clo_dist = sqrt(dist2);
//     return;
// }

//

void SplineSurface::closestBoundaryPoint(const Point& pt,
					 double& clo_u,
					 double& clo_v, 
					 Point& clo_pt,
					 double& clo_dist,
					 double epsilon,
					 const RectDomain* rd,
					 double *seed) const
///////////////////////////////////////////////////////////////////////////////
//
// PURPOSE:	Given a Nurbsurface and a point it calculates the 
//              closest point on the boundaries of the surface.
// IMPORTANT:	It calculates only a local solution using Newton iteration
//              on the surface boundaries.
//              Very ineffective, should be changed later. 
//				
//
// WRITTEN BY:  Vibeke Skytt, SINTEF, 0205
//
///////////////////////////////////////////////////////////////////////////////
{
    RectDomain domain = containingDomain();
    if (!rd)
	rd = &domain;

    Point cpt;
    double cdist, cpar;

    // Check degeneracy
    bool b, r, t, l;
    //   double tol = 0.000001;  // Arbitrary tolerance. The information should
    //                           // be present.
    (void)isDegenerate(b, r, t, l, epsilon);

    // Checking closest point on the bottom boundary
    shared_ptr<SplineCurve> bdcrv;
    clo_dist = 1.0e10;  // Initialize with a large number
    if (b == false)
	{
	    bdcrv = shared_ptr<SplineCurve>(constParamCurve(rd->vmin(), true));
	    bdcrv->closestPoint(pt, rd->umin(), rd->umax(),
				clo_u, clo_pt, clo_dist, seed);
	    clo_v = rd->vmin();
	}

    // Checking the right boundary
    if (r == false)
	{
	    bdcrv = shared_ptr<SplineCurve>(constParamCurve(rd->umax(),
							    false));
	    bdcrv->closestPoint(pt, rd->vmin(), rd->vmax(), cpar, cpt, cdist,
				(seed == 0) ? seed : seed+1);
	    if (cdist < clo_dist)
		{
		    clo_pt = cpt;
		    clo_u = rd->umax();
		    clo_v = cpar;
		    clo_dist = cdist;
		}
	}

    // Checking the upper boundary
    if (t == false)
	{
	    bdcrv = shared_ptr<SplineCurve>(constParamCurve(rd->vmax(), true));
	    bdcrv->closestPoint(pt, rd->umin(), rd->umax(), cpar, cpt, cdist,
				seed);
	    if (cdist < clo_dist)
		{
		    clo_pt = cpt;
		    clo_u = cpar;
		    clo_v = rd->vmax();
		    clo_dist = cdist;
		}
	}

    // Checking the left boundary
    if (l == false)
	{
	    bdcrv = shared_ptr<SplineCurve>(constParamCurve(rd->umin(),
							    false));
	    bdcrv->closestPoint(pt, rd->vmin(), rd->vmax(), cpar, cpt, cdist,
				(seed == 0) ? seed : seed+1);
	    if (cdist < clo_dist)
		{
		    clo_pt = cpt;
		    clo_u = rd->umin();
		    clo_v = cpar;
		    clo_dist = cdist;
		}
	}
  
}

void 
SplineSurface::s1773(const double ppoint[],double aepsge, double estart[],double eend[],double enext[],
		     double gpos[],int *jstat) const
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a surface and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameters is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : ppoint   - The point in the closest point problem.
*              psurf    - The surface in the closest point problem.
*              aepsge   - Geometry resolution.
*              estart   - Surface parameters giving the start of the search
*                         area (umin, vmin).
*              eend     - Surface parameters giving the end of the search
*                         area (umax, vmax).
*              enext    - Surface guess parameters for the closest point
*                         iteration.
*
*
*
* OUTPUT     : gpos    - Resulting surface parameters from the iteration.
*              jstat   - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in two parameter direction.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, May 1989
* Revised by : Johannes Kaasa, SINTEF Oslo, August 1995.
*              Introduced a local copy of enext, to avoid changes.
*
*********************************************************************
*/                       
{                        
  // int kstat = 0;            /* Local status variable.                      */
  int kder=1;               /* Order of derivatives to be calulated        */
  int kdim;                 /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  int kdir;                 /* Changing direction.                         */
  int kdeg;                 /* Degenaracy flag.                            */
  double tdelta[2];         /* Parameter intervals of the surface.         */
  double tdist;             /* Distance between position and origo.        */
  double td[2],t1[2],tdn[2];/* Distances between old and new parameter
			       value in the tree parameter directions.     */
  double tprev;             /* Previous difference between the curves.     */
  vector<double> sdiff;      /* Difference between the point and the surf.  */
  double snext[2];          /* Parameter values                            */
  double guess[2];          /* Local copy of enext.                        */
  double REL_COMP_RES = 1.0e-15;
  vector<Point> pts(3);
  
  guess[0] = enext[0];
  guess[1] = enext[1];
  
  kdim = dimension();
  
  
  /* Fetch endpoints and the intervals of parameter interval of curves.  */
  
  tdelta[0] = endparam_u() - startparam_u();
  tdelta[1] = endparam_v() - startparam_v();
  
  /* Allocate local used memory */
  sdiff.resize(kdim);
  
  /* Initiate variables.  */
  
  tprev = 1.0e10;
  
  /* Evaluate 0-1.st derivatives of surface */
  /* printf("\n lin: \n %#20.20g %#20.20g",
     guess[0],guess[1]); */
  
  point(pts, guess[0], guess[1], kder);
  
  /* Compute the distanse vector and value and the new step. */
  
  s1773_s9dir(&tdist,td,td+1,&sdiff[0],ppoint,pts,
	      aepsge,kdim,&kdeg);
  
  /* Correct if we are not inside the parameter intervall. */
  
  t1[0] = td[0];
  t1[1] = td[1];
  s1773_s9corr(t1,guess[0],guess[1],estart[0],eend[0],estart[1],eend[1]);
  
  /* Iterate to find the intersection point.  */
  
  for (knbit = 0; knbit < 30; knbit++)
    {
      /* Evaluate 0-1.st derivatives of surface */
      
      snext[0] = guess[0] + t1[0];
      snext[1] = guess[1] + t1[1];
      
      point(pts, snext[0], snext[1], kder);
      
      /* Compute the distanse vector and value and the new step. */
      
      s1773_s9dir(&tdist,tdn,tdn+1,&sdiff[0],ppoint,
	    pts,aepsge,kdim,&kdeg);
      
      /* Check if the direction of the step have change. */
      
      kdir = (Utils::inner(td, td+2, tdn) >= 0.0);     /* 0 if changed. */
      
      /* Ordinary converging. */
      
      if (tdist < tprev/(double)2 || kdir)
	{
	   guess[0] += t1[0];
	   guess[1] += t1[1];
  
	  /* printf("\n %#20.20g %#20.20g",
	     guess[0],guess[1]); */
  
	  
          td[0] = t1[0] = tdn[0];
          td[1] = t1[1] = tdn[1];
	  
	  /* Correct if we are not inside the parameter intervall. */
	  
	  s1773_s9corr(t1,guess[0],guess[1],estart[0],eend[0],estart[1],eend[1]);
          tprev = tdist;

	  if ( (fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES)) break;
	}
      
      /* Not converging, adjust and try again. */
      
      else
	{
          t1[0] /= (double)2;
          t1[1] /= (double)2;
          /* knbit--;  */
	}
      if (guess[0]==guess[0]+t1[0] &&
	  guess[1]==guess[1]+t1[1]) break;
    }
  
  /* Iteration stopped, test if point founds found is within resolution */
  
  if (tdist <= aepsge)
  {
     *jstat = 1;
     /* printf("\n SUCCESS!!"); */
     
  }
  else if(kdeg)
     *jstat = 9;
  else
     *jstat = 2;
  
  gpos[0] = guess[0];
  gpos[1] = guess[1];
  

  return;


  // /* Iteration completed.  */
  
  // goto out;
  
//   /* Error in allocation */
  
//  err101: *jstat = -101;
// //  s6err("s1773",*jstat,kpos);
//   goto out;                  
  
//   /* Error in input. Conflicting dimensions.  */
  
//  err106: *jstat = -106;
// //  s6err("s1773",*jstat,kpos);
//   goto out;                  
  
//   /* Error in lower level routine.  */
  
//   error : *jstat = kstat;
// //  s6err("s1773",*jstat,kpos);
//   goto out;                  
  
//  out:    
//   return;
}

void
SplineSurface::s1773_s9corr(double gd[],double acoef1,double acoef2,
			    double astart1,double aend1,double astart2,double aend2) const
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To be sure that we are inside the boorder of the
*              parameter plan. If we are outside clipping is used
*	       to corrigate the step value.
*
*
* INPUT      : acoef1  - Coeffisient in the first direction.
*              acoef2  - Coeffisient in the second direction.
*              astart1 - The lower boorder in first direction.
*              aend1   - The higher boorder in first direction.
*              estart2 - The lower boorder in second direction.
*              eend2   - The higher boorder in second direction.
*
*
*
* INPUT/OUTPUT : gd    - Old and new step value.
*
*
* METHOD     : We are cutting a line inside a rectangle.
*	       In this case we always know that the startpoint of
*	       the line is inside the rectangel, and we may therfor
*	       use a simple kind of clipping.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
*
*********************************************************************
*/                       
{
  if (acoef1 + gd[0] < astart1)  gd[0] = astart1 - acoef1;
  else if (acoef1 + gd[0] > aend1) gd[0] = aend1 - acoef1;
  
  if (acoef2 + gd[1] < astart2)  gd[1] = astart2 - acoef2;
  else if (acoef2 + gd[1] > aend2) gd[1] = aend2 - acoef2;
}

void
SplineSurface::s1773_s9dir(double *cdist,double *cdiff1,double *cdiff2,
			   double PS[],const double *eval1,vector<Point> eval2,
			   double aepsge, int idim,int *jstat) const
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the distance vector and value beetween
*	       a point and a point on a surface.
*	       And to compute a next step on both parameter direction
*	       This is equivalent to the nearest way to the
*	       parameter plan in the tangent plan from a point in the
*	       distance surface between a point and a surface.
*
*
* INPUT      : eval1 - Value in point.
*              eval2 - Value +1 and 2. der in surface.
*	       aepsge- Geometry tolerance.
*	       idim  - Dimension of space the surface lie in.
*
*
* OUTPUT     : PS   - Array to use when computing the differens vector.
*	       cdiff1  - Relative parameter value in intersection 
*                        point in first direction.
*              cdiff2  - Relative parameter value in intersection 
*                        point in second direction.
*              cdist   - The value to the point in the distance surface.
*              jstat   - 0 OK, new No degeneracy.
*                        1 Degeneracy.
*
*
* METHOD     : The method is to compute the parameter distance to the points
*	       on both tangents which is closest to the point.
*	       The difference vector beetween these points are orthogonal
*	       to both tangents. If the distance vector beetween the point and
*	       point on the surface is "diff" and the two derivativ vectors
*	       are "der1" and "der2", and the two wanted parameter distance
*	       are "dt1" and "dt2", then we get the following system of 
*	       equations:
*		 <-dist+dt1*der1+dt2*der2,der1> = 0
*		 <-dist+dt1*der1+dt2*der2,der2> = 0
*	       This is futher:
*
*		 | <der1,der1>   <der1,der2> |  | dt1 |   | <dist,der1> |
*		 |                           |  |     | = |     	|
*		 | <der1,der2>   <der2,der2> |  | dt2 |   | <dist,der2> |
*
*	       The solution of this matrix equation is the
*	       following function.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
*
*********************************************************************
*/                       
{                        
  // int kstat=0;		          /* Local status variable.       */
  register double tdet;		  /* Determinant                  */
  register double t1,t2,t3,t4,t5; /* Variables in equation system */
  register double *S, *Su, *Sv;
  /* register double *Suv, *Suu, *Svv; */
                                  /* Pointers to surf values      */
  register double ref, ang;       /* Referance value, angle       */
  register double l1, l2;         /* Vector norm                  */
  register double tcos;
  register double min_ang=10e-11; /* Min angle                    */
  register double ptol = 1.0e-12; /* Replace DEQUAL               */
  /* ____________________________________________________________ */
  
  /* Init */
  *jstat = 0;
  *cdiff1 = 0.0;
  *cdiff2 = 0.0;
  
  /* Set pointers */
  S   = eval2[0].begin();
  Su  = eval2[1].begin();
  Sv  = eval2[2].begin();
  /* Suu = Sv  + idim;
  Suv = Suu + idim;
  Svv = Suv + idim; */

  /* Degenerate if Su=0 v Sv=0 v Su||Sv */
  l1 = sqrt(Utils::inner(Su, Su+idim, Su));
  l2 = sqrt(Utils::inner(Sv, Sv+idim, Sv));
  tcos = Utils::inner(Su, Su+idim, Sv)/(l1*l2);
  ang = acos(tcos);
  if (std::min(l1,l2) < aepsge || ang < min_ang) *jstat = 1;

  /* Computing difference vector and lenght */
  for (int ki=0; ki<idim; ki++)
      PS[ki] = eval1[ki] - S[ki];
  *cdist = sqrt(Utils::inner(PS, PS+idim, PS));
  
  if (*jstat == 1)
  {
     if (l1 < aepsge)
     {
	if (l2 > aepsge)
	   /* Su = 0 */
	   *cdiff2 = Utils::inner(PS, PS+idim, Sv)/l2*l2;
     }
     else if (l2 < aepsge)
	   /* Sv = 0 */
	   *cdiff1 = Utils::inner(PS, PS+idim, Su)/(l1*l1);
     else /* Su,Sv || */
     {
	/* Best strategy? */
	*cdiff1 = Utils::inner(PS, PS+idim, Su)/(l1*l1);
      }
	
  }
  else /* *jstat == 0 */
     
  {
     
     t1 =  Utils::inner(Su, Su+idim, Su);
     t2 =  Utils::inner(Su, Su+idim, Sv);
     t3 =  Utils::inner(Sv, Sv+idim, Sv);
     t4 =  Utils::inner(PS, PS+idim, Su);
     t5 =  Utils::inner(PS, PS+idim, Sv);
     
     ref = std::max(fabs(t1),fabs(t2));
     ref = std::max(ref,fabs(t3));
     /* Computing the determinant. */
     
     tdet = t1*t3 - t2*t2;
     
     if (fabs(tdet) < ptol)
     {
	*jstat = 1;
     }
     else 
     {
	/* Using Cramer's rule to find the solution of the system. */
	
	*cdiff1 =  (t4*t3-t5*t2)/tdet;
	*cdiff2 =  (t1*t5-t2*t4)/tdet;
     }
  }
}
} // namespace Go


