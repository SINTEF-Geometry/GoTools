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

#include "GoTools/intersections/SfCvIntersector.h"
#include "GoTools/intersections/CvCvIntersector.h"
#include "GoTools/intersections/SfPtIntersector.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/intersections/GeoTol.h"
#include "GoTools/utils/RotatedBox.h"
#include "GoTools/utils/Values.h"


using namespace Go;

using std::vector;
using std::cout;
using std::endl;


namespace { // anonymous namespace 

// distance function between two curves.  Used by the minimization algorithm
// initiated by ClosestPoint::closestPtCurves.
class CrvSrfDistFun {
public:
    CrvSrfDistFun(const ParamSurface* sf, 
		  const ParamCurve* cv,
		  const double* const minpar = 0,
		  const double* const maxpar = 0);
    
    inline double operator()(const double* arg) const;
    inline double grad(const double* arg, double* res) const;
    inline double minPar(int pardir) const;
    inline double maxPar(int pardir) const;

private:
    double minpar_[6];
    double maxpar_[6];
    const ParamSurface * const sf_;
    const ParamCurve * const cv_;
    mutable Point p1_, p2_, d_;
    mutable vector<Point> p1vec_, p2vec_;
};


} // end anonymous namespace

namespace Go {


//===========================================================================
SfCvIntersector::SfCvIntersector(shared_ptr<ParamGeomInt> obj1, 
				 shared_ptr<ParamGeomInt> obj2,
				 double epsge, Intersector* prev,
				 int eliminated_parameter,
				 double eliminated_value)
  : Intersector2Obj(obj1, obj2, epsge, prev,
		    eliminated_parameter, eliminated_value)
//===========================================================================
{
    // Set array indices to the object array. We know that one object
    // is a curve an the other is a surface, but need to find out
    // which one is which
    ParamCurveInt *curve = obj1->getParamCurveInt();
    if (curve == 0) {
	cv_idx_ = 1;
	sf_idx_ = 0;
    } else {
	cv_idx_ = 0;
	sf_idx_ = 1;
    }
}


//===========================================================================
SfCvIntersector::SfCvIntersector(shared_ptr<ParamGeomInt> obj1, 
				 shared_ptr<ParamGeomInt> obj2,
				 shared_ptr<GeoTol> epsge, Intersector* prev,
				 int eliminated_parameter,
				 double eliminated_value)
  : Intersector2Obj(obj1, obj2, epsge, prev,
		    eliminated_parameter, eliminated_value)
//===========================================================================
{
    ParamCurveInt *curve = obj1->getParamCurveInt();
    if (curve == 0) {
	cv_idx_ = 1;
	sf_idx_ = 0;
    } else {
	cv_idx_ = 0;
	sf_idx_ = 1;
    }
}


//===========================================================================
SfCvIntersector::~SfCvIntersector()
//===========================================================================
{
    // Currently empty
}


//===========================================================================
shared_ptr<Intersector> 
SfCvIntersector::lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
				       shared_ptr<ParamGeomInt> obj2, 
				       Intersector* prev,
				       int eliminated_parameter,
				       double eliminated_value)
//===========================================================================
{
  // Either a CvCvIntersector or a SfPtIntersector. 

    ASSERT(eliminated_parameter >= 0 && eliminated_parameter <= 2);

    if ((cv_idx_ == 0 && eliminated_parameter > 0) ||
	 (cv_idx_ == 1 && eliminated_parameter < 2))
    {
	shared_ptr<CvCvIntersector> curr_inter; 
	curr_inter = shared_ptr<CvCvIntersector>
	    (new CvCvIntersector(obj1, obj2, epsge_, prev,
				 eliminated_parameter, eliminated_value));

	return curr_inter;
    }
    else
    {
	shared_ptr<SfPtIntersector> curr_inter; 
	curr_inter = shared_ptr<SfPtIntersector>
	    (new SfPtIntersector(obj1, obj2, epsge_, prev,
				 eliminated_parameter, eliminated_value));

	return curr_inter;
    }

}


//===========================================================================
int SfCvIntersector::checkCoincidence()
//===========================================================================
{
    // What about linearity?

    // Can intersection point occur in the inner of the curve/surface?
    // In that case all inner intersection points can be fetched an
    // sorted according the the parameter belonging to the curve
    // Coincidence can be checked between two and two points in
    // sequence.  If all intervals are coincident, then there is total
    // coincidence.  Otherwise, connections are made, but the return
    // value is no coincidence.

    double ptol = 100.0*epsge_->getRelParRes();
    double res = epsge_->getRelParRes();
    ParamCurveInt *curve = obj_int_[cv_idx_]->getParamCurveInt();
    ParamSurfaceInt *surf = obj_int_[sf_idx_]->getParamSurfaceInt();

    // Fetch intersection points at the endpoints/boundaries
    vector<shared_ptr<IntersectionPoint> > bd_ints;
    int_results_->getIntersectionPoints(bd_ints);

    if (bd_ints.size() < 2)
	return 0;  // Not sufficient number of intersection points for
		   // coincidence

    if (bd_ints.size() > 2) {
	// Sort according to the parameter direction of the curve
	for (size_t ki=0; ki<bd_ints.size(); ki++) {
	    for (size_t kj=ki+1; kj<bd_ints.size(); kj++) {
		if (bd_ints[ki]->getPar(2*cv_idx_)
		    > bd_ints[kj]->getPar(2*cv_idx_)) {
		    std::swap(bd_ints[ki], bd_ints[kj]);
		}
	    }
	}
    }

//     // Check if the most distinct intersection points lie on endpoints
//     // in the curve or corners in the surface.
//     bool corner = surf->inCorner((sf_idx_==0) ? bd_ints[0]->getPar1() :
// 				 bd_ints[0]->getPar2(), 
// 				 epsge_->getRelParRes());
//     bool endpt = curve->boundaryPoint((cv_idx_==0) ? bd_ints[0]->getPar1() :
// 				      bd_ints[0]->getPar2(),
// 				      epsge_->getRelParRes());
// //     if (!corner && !endpt)
// // 	return 0;

//     int size = (int)(bd_ints.size());
//     corner = surf->inCorner((sf_idx_==0) ? bd_ints[size-1]->getPar1() :
// 			    bd_ints[size-1]->getPar2(), 
// 			    epsge_->getRelParRes());
//     endpt = curve->boundaryPoint((cv_idx_==0) ? bd_ints[size-1]->getPar1() :
// 				 bd_ints[size-1]->getPar2(),
// 				 epsge_->getRelParRes());
// //     if (!corner && !endpt)
// // 	return 0;

    // For all intervals between intersection points
    bool coinc = true;
    bool midpoint_coinc = true;
    bool is_degenerate;
    int is_coincident = 0; // Initialized to no coincidence
    size_t ksize = bd_ints.size();
    vector<double> start, end;

    // Check the self intersection case first
    if (selfint_case_)
    {
	start = bd_ints[0]->getPar();
	end = bd_ints[ksize-1]->getPar();

	// Possibly and intersection between a surface and its boundary
	// curves or between two neighbouring surfaces. Check if an
	// alternative and less risky test method can be used.
	double ptol = epsge_->getRelParRes();
	if (fabs(start[sf_idx_]-end[sf_idx_]) < ptol)
	{
	    is_coincident  = surf->isIsoParametric(curve, 0, start[sf_idx_],
						   ptol,
						   epsge_->getEpsge());

	    // Check if the endpoints indicate selfintersection type coincidence
	    if (fabs(bd_ints[0]->getPar(sf_idx_+1) - bd_ints[0]->getPar(2*cv_idx_)) > res ||
		fabs(bd_ints[ksize-1]->getPar(sf_idx_+1) - bd_ints[ksize-1]->getPar(2*cv_idx_)) > res)
		is_coincident = 0;
	}
	else if (fabs(start[sf_idx_+1]-end[sf_idx_+1]) < ptol)
	{
	    is_coincident  = surf->isIsoParametric(curve, 1, start[sf_idx_+1],
						   ptol,
						   epsge_->getEpsge());

	    // Check if the endpoints indicate selfintersection type coincidence
	    if (fabs(bd_ints[0]->getPar(sf_idx_) - bd_ints[0]->getPar(2*cv_idx_)) > res ||
		fabs(bd_ints[ksize-1]->getPar(sf_idx_) - bd_ints[ksize-1]->getPar(2*cv_idx_)) > res)
		is_coincident = 0;
	}
	else
	    is_coincident = 0;
    }

    if (is_coincident)
    {
	// Register the coincidence by connecting the
	// intersection points
	bd_ints[0]->connectTo(bd_ints[ksize-1], COINCIDENCE_SFCV);
	coinc = true;
	midpoint_coinc = true;
    }
    else
    {
	for (size_t ki=1; ki<ksize; ki++) {
	    is_degenerate = false;  // Initally

	    // Check if the points are connected already
	    if (bd_ints[ki-1]->isConnectedTo(bd_ints[ki]));  //@@@ Pointer
	    //or shared
	    //ptr?
	    else {
   
		// Perform marching to check coincidence
		//int cv_idx;
		Point par1(2), par2(2);
		start = bd_ints[ki-1]->getPar();
		end = bd_ints[ki]->getPar();
		//double eps = epsge_->getEpsge();

		par1.setValue(&start[sf_idx_]);
		par2.setValue(&end[sf_idx_]);
		DEBUG_ERROR_IF(curve == 0 || surf == 0,
			       "Error in data structure");

		// Check midpoint
		double cv_par1 = start[2*cv_idx_];
		double cv_par2 = end[2*cv_idx_];
		if (fabs(cv_par2-cv_par1) < epsge_->getRelParRes()) {
		    // Identical curve parameters. Check the midpoint of
		    // the constant curve between the intersection points
		    // on the surface
		    Point mid_par = 0.5*(par1+par2);
		    Point mid_pt;
		    surf->point(mid_pt, mid_par.begin());

		    if (mid_pt.dist(bd_ints[ki]->getPoint()) > epsge_->getEpsge())
		    {
			midpoint_coinc = false;   // Not a degenerate case
			continue;
		    }

		    // Iterate to the curve
		    double tmin, tmax;
		    shared_ptr<ParamCurve> cv
			= curve->getParentParamCurve(tmin, tmax);
		    tmin = cv_par1;
		    tmax = cv_par2;

		    Point clo_pt;
		    double clo_par, clo_dist;
		    cv->closestPoint(mid_pt, tmin, tmax,
				     clo_par, clo_pt, clo_dist, &cv_par1);
		    if (clo_dist > epsge_->getEpsge()) {
			midpoint_coinc = false;
			continue;
		    }
		    else
			is_degenerate = true;
		}
		else if (par1.dist(par2) < ptol) {
		    // Check the midpoint of the curve with the surface
		    double mid_par = 0.5*(cv_par1 + cv_par2);
		    Point mid_pt;
		    curve->point(mid_pt, &mid_par);

		    if (mid_pt.dist(bd_ints[ki]->getPoint()) > epsge_->getEpsge())
		    {
			midpoint_coinc = false;   // Not a degenerate case
			continue;
		    }

		    RectDomain domain;
		    shared_ptr<ParamSurface> sf
			= surf->getParentParamSurface(domain);
		    Point seed = 0.5*(par1+par2);

		    Point clo_pt;
		    double clo_par[2], clo_dist;
		    sf->closestPoint(mid_pt, clo_par[0], clo_par[1],
				     clo_pt, clo_dist,
				     epsge_->getEpsge(), &domain, seed.begin());
		    if (clo_dist > epsge_->getEpsge()) {
			midpoint_coinc = false;
			continue;
		    }
		    else
			is_degenerate = true;
		}

		is_coincident = 0; // Initialized to no coincidence
		if (is_degenerate)
		    is_coincident = 1;
		else
		{
		    if (!is_coincident)
		    {
			try {
			    is_coincident = checkCoincide(curve, start[2*cv_idx_], 
							  end[2*cv_idx_], 
							  surf, par1, par2, epsge_);
			} 
			catch (...) {
			    // May throw here. If so, keep is_coincident = 0. @@@jbt
			    is_coincident = 0;
			}
		    }
		}

		/*if (!is_coincident && curve->isSpline())
		{
		    const SplineCurve *spline = curve->getSpline();
		    if (spline->order() == 2)
		    {
			// VSK 071025. This is a hack to avoid infinite recursion
			// in self intersection context
			// This configuration doesn't make sense. Connect to avoid
			// infinite recusion
			is_coincident = 1;
		    }
		    }*/

		if (is_coincident) {
		    // Register the coincidence by connecting the
		    // intersection points
		    bd_ints[ki-1]->connectTo(bd_ints[ki], COINCIDENCE_SFCV);
		}
		else {
		    coinc = false;
		}
	    }
	}
    }

    // A coincidence interval is recognized. It is necessary to check
    // for the possibility of additional intersections off the
    // coincidence interval.  Check if the most distinct intersection
    // points lie on endpoints in the curve or corners in the surface.
    // VSK, 0611. If the curve is linear, the risk of more intersection
    // points is neglectable
    if (coinc && midpoint_coinc) {
	int size = (int)bd_ints.size();
	double ptol2 = epsge_->getRelParRes();
	double* pt_par1 = const_cast<double*>((sf_idx_==0) ? 
					      bd_ints[0]->getPar1() :
					      bd_ints[0]->getPar2());
	double* pt_par2 = const_cast<double*>((sf_idx_==0) ?
					     bd_ints[size-1]->getPar1() :
					     bd_ints[size-1]->getPar2());
	double* par2 = const_cast<double*>((cv_idx_ == 0) ? bd_ints[0]->getPar1() : bd_ints[0]->getPar2());
	bool corner1 = surf->inCorner(pt_par1, ptol2);
	bool endpt1 = curve->boundaryPoint(par2, ptol2);
	bool opposite_bd = false;
	if ((pt_par1[0]-surf->startParam(0)<ptol && surf->endParam(0)-pt_par2[0]<ptol2) ||
	    (pt_par2[0]-surf->startParam(0)<ptol && surf->endParam(0)-pt_par1[0]<ptol2))
	    opposite_bd = true;
	if ((pt_par1[1]-surf->startParam(1)<ptol && surf->endParam(1)-pt_par2[1]<ptol2) ||
	    (pt_par2[1]-surf->startParam(1)<ptol && surf->endParam(1)-pt_par1[1]<ptol2))
	    opposite_bd = true;
	bool cv_lin = curve->isLinear(epsge_->getEpsge());

	if (!corner1 && !endpt1 && !opposite_bd && !cv_lin)
	{
//	    curve->setCriticalVal(0, *par2);
	    return 0;
	}

	par2 = const_cast<double*>((cv_idx_ == 0) ? bd_ints[size-1]->getPar1() : 
				   bd_ints[size-1]->getPar2());
	bool corner2 = surf->inCorner(pt_par2, ptol2);
	bool endpt2 = curve->boundaryPoint(par2, ptol2);
	if (!corner2 && !endpt2 && !opposite_bd && !cv_lin)
	{
//	    curve->setCriticalVal(0, *par2);
	    return 0;
	}

	// Check if the routine is called in a self-intersection
	// context and the coincidence curve follows a surface
	// boundary.
	if (selfint_case_ == 1 &&
	    (corner1 || endpt1) && (corner2 || endpt2)) {
	    if (int_results_->isBoundaryIntersection(bd_ints[0],
						     bd_ints[size-1])) {

		// Check cone for possibility of self intersection
		DirectionCone cone1, cone2;
		try {
		    cone1 = obj_int_[0]->directionCone();
		    cone2 = obj_int_[1]->directionCone();
		} catch (...) {
		    return 0;  // Continue subdividing
		}
		if (cone1.greaterThanPi() || cone2.greaterThanPi())
		    return 0;  // Possibilities for self intersection
		return 1;      // Accept coincidence interval
	    }
	}

	return 1;
    }
    return (coinc && midpoint_coinc);
}


//===========================================================================
void SfCvIntersector::microCase()
//===========================================================================
{
    // Purpose : The curve and the surface lie both within the same
    // epsion ball.  It is not a simple case, and there is
    // intersection.  This is an uwanted situation. We have to make a
    // result that is as consistent as possible since we cannot
    // subdivide anymore.

    // Fetch all intersection points belonging to the current object
    // The intersection points should be sorted according to the
    // parameter of the curve
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getSortedIntersections(int_pts);

    if (int_pts.size() == 0) {
	// No intersection point exist. Construct one internal to the
	// two objects.  It is probably good enough to take the
	// midpoints since the objects are so small
	int nmbdir1 = obj_int_[0]->numParams();
	double mid_par[3];

// 	double start_par, end_par;
// 	int nmbdir2;
// 	nmbdir1 = obj_int_[0]->numParams();
// 	nmbdir2 = obj_int_[1]->numParams();
// 	for (int ki=0; ki<nmbdir1; ki++) {
// 	    start_par = obj_int_[0]->startParam(ki);
// 	    end_par = obj_int_[0]->endParam(ki);
// 	    mid_par[ki] = 0.5*(start_par + end_par);
// 	}
// 	for (int ki=0; ki<nmbdir2; ki++) {
// 	    start_par = obj_int_[1]->startParam(ki);
// 	    end_par = obj_int_[1]->endParam(ki);
// 	    mid_par[nmbdir1+ki] = 0.5*(start_par + end_par);
// 	}

	double dist;
	doIterate(mid_par, dist);
	if (dist < getTolerance()->getEpsge()) {
	    int_results_->addIntersectionPoint(obj_int_[0], 
					       obj_int_[1],
					       getTolerance(),
					       mid_par, 
					       mid_par+nmbdir1); 
	}
    } else if (int_pts.size() == 1) {
	// One intersection point. Leave it like this.
    } else {
	// More than one intersection point. Connect.
	for (int ki=1; ki<int(int_pts.size()); ki++) {
	    int_pts[ki-1]->connectTo(int_pts[ki], MICRO_SFCV);
	}
    }
}


//===========================================================================
int SfCvIntersector::updateIntersections()
//===========================================================================
{
    // Purpose : Given one parametric curve and one parametric surface
    // in a simple case situation, iterate to the intersection point,
    // if any.

    // Fetch already existing intersection points
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    // VSK, 0609. Post iterate existing intersection points
    int ki, kj;
    double param[3], guess[3];
    double dist;
    double ptol = epsge_->getRelParRes();
    double ta[3], tb[3];
    int npar1 = obj_int_[0]->numParams();
    int npar2 = obj_int_[1]->numParams();
    for (ki=0; ki<npar1; ki++)
    {
	ta[ki] = obj_int_[0]->startParam(ki);
	tb[ki] = obj_int_[0]->endParam(ki);
    }
    for (ki=0; ki<npar2; ki++)
    {
	ta[npar1+ki] = obj_int_[1]->startParam(ki);
	tb[npar1+ki] = obj_int_[1]->endParam(ki);
    }

    for (ki=0; ki<(int)(int_pts.size()); ki++)
    {
	for (kj=0; kj<3; kj++)
	{
	    guess[kj] = int_pts[ki]->getPar(kj);
	    if (fabs(guess[kj]-ta[kj])<ptol || fabs(guess[kj]-tb[kj]) < ptol)
		break;
	}
	if (kj<3)
	    continue; // Do not move boundary points

	doIterate(param, dist, guess);
	for (kj=0; kj<3; kj++)
	    if ((fabs(guess[kj]-ta[kj])<ptol || fabs(guess[kj]-tb[kj])) &&
		!(fabs(param[kj]-ta[kj])<ptol || fabs(param[kj]-tb[kj])))
		break;
	if (dist < int_pts[ki]->getDist() && kj==3)
	    int_pts[ki]->replaceParameter(param);
    }

    if (int_pts.size() > 1)
    {
	// VSK, 0610. Check if the points are connected
	for (ki=1; ki<(int)(int_pts.size()); ki++)
	    if (!int_pts[ki-1]->isConnectedTo(int_pts[ki]))
		break;
	if (ki == (int)(int_pts.size()))
	    return 0;   // Nothing more to do
	
	// VSK, 0609. Only one intersection points is expected. Identify points.
	bool at_bd1, at_bd2;
	for (kj=0; kj<3; kj++)
	{
	    double t1 = int_pts[0]->getPar(kj);
	    if (fabs(t1-ta[kj])<ptol || fabs(t1-tb[kj])<ptol)
		break;
	}
	at_bd1 = (kj<3);

	for (ki=1; ki<(int)(int_pts.size()); ki++)
	{
	    for (kj=0; kj<3; kj++)
	    {
		double t1 = int_pts[ki]->getPar(kj);
		if (fabs(t1-ta[kj])<ptol || fabs(t1-tb[kj])<ptol)
		    break;
	    }
	    at_bd2 = (kj<3);

	    if (int_pts[ki]->getDist() < int_pts[ki-1]->getDist())
	    {
		int_results_->removeIntPoint(int_pts[ki-1]);
		int_pts.erase(int_pts.begin()+ki-1);
		ki--;
	    }
	    else
	    {
		int_results_->removeIntPoint(int_pts[ki]);
		int_pts.erase(int_pts.begin()+ki);
		ki--;
		at_bd1 = at_bd2;
	    }
	}
    }
		

    if (int_pts.size() > 0) {
	// At most one intersection point is expected. One or more
	// points exist already. Leave it like this.
	
	return 0;
    }


    // Iterate to an intersection point between the curve and the
    // surface
    double par[3];
   ptol = 100.0*epsge_->getRelParRes();
   doIterate(par, dist);

    if (dist <= epsge_->getEpsge()) {

	if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
	    cout << "SfCv. Int. pt. found " << par[0] << " " << par[1];
	    cout << " " << par[2] << ", dist: " << dist;
	}
	// An intersection point is found. Represent it in the data
	// structure.
	// @@@ vsk, In sisl there is a check if the found intersection
	// point is very close to the endpoints of the curve. In that
	// case it is dismissed since the intersections at the
	// endpoints are found already. Here we have already checked
	// if there exist any intersection points. Thus, we know that
	// there does not exist any intersection point at the
	// endpoint.
	// VSK, 0906 But maybe the endpoint already is in intersection, but that
	// point is represented elsewhere at a more exact position. Make a check.
	int pc = (cv_idx_ == 0) ? 0 : 2;
	int ps = (cv_idx_ == 0) ? 1 : 0;
	if (par[pc]-obj_int_[cv_idx_]->startParam(0) < ptol || 
	    obj_int_[cv_idx_]->endParam(0)-par[pc] < ptol ||
	    par[ps]-obj_int_[sf_idx_]->startParam(0) < ptol || 
	    obj_int_[sf_idx_]->endParam(0)-par[ps] < ptol||
	    par[ps+1]-obj_int_[sf_idx_]->startParam(1) < ptol || 
	    obj_int_[sf_idx_]->endParam(1)-par[ps+1] < ptol)
	{
	    if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
		cout << " Point dismissed" << std::endl;
	    }
	    return 0;
	}

	if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
	    cout << " Point registered" << std::endl;
	    }
	shared_ptr<IntersectionPoint> tmp = 
	    int_results_->addIntersectionPoint(obj_int_[0], 
					       obj_int_[1], 
					       getTolerance(),
					       par, 
					       par + obj_int_[0]
					       ->numParams());
	return 1; // A new point is found
    }

    return 0;
}


//===========================================================================
void SfCvIntersector::doIterate(double par[], double& dist, double *guess)
//===========================================================================
{
    // First fetch the parametric curve and surface. We need both the
    // current parametric curve to compute good start values for the
    // iteration and the initial curve where no numerical noice is
    // added trough subdivision to get a numerically stable result
    // from the iteration.
    double start, end;  // The parametric limitation of the current
                        // curve in the parameter interval of the
                        // parent curve
    RectDomain domain;  // Similar for the surface
    shared_ptr<ParamCurve> curve = 
	obj_int_[cv_idx_]->getParamCurveInt()
	->getParentParamCurve(start, end);
    shared_ptr<ParamSurface> surf = 
	obj_int_[sf_idx_]->getParamSurfaceInt()
	->getParentParamSurface(domain);
    DEBUG_ERROR_IF(curve.get() == 0 || surf.get() == 0,
		   "Inconsistence in the data structure");

    // Compute a startpoint for the iteration
    double seed[3];
    bool second_order = false;
    if (guess) {
	int ki;
	for (ki=0; ki<3; ki++) {
	    seed[ki] = guess[ki];
	}
	second_order = true;
    } else {
	getSeedIteration(seed);
    }

    // Iterate
    Point ptc, pts;
    ClosestPoint::closestPtCurveSurf(curve.get(), surf.get(), epsge_->getEpsge(),
		       start, end, &domain, seed[2*cv_idx_], seed+sf_idx_, 
		       par[2*cv_idx_], par+sf_idx_, dist, ptc, pts, 
		       second_order);
}


//===========================================================================
void SfCvIntersector::doIterate2(double par[], double& dist, double guess[])
//===========================================================================
{
    // First fetch the parametric curve and surface. We need both the
    // current parametric curve to compute good start values for the
    // iteration and the initial curve where no numerical noice is
    // added trough subdivision to get a numerically stable result
    // from the iteration.
    double start, end;  // The parametric limitation of the current
                        // curve in the parameter interval of the
                        // parent curve
    RectDomain domain;  // Similar for the surface
    shared_ptr<ParamCurve> curve = 
	obj_int_[cv_idx_]->getParamCurveInt()
	->getParentParamCurve(start, end);
    shared_ptr<ParamSurface> surf = 
	obj_int_[sf_idx_]->getParamSurfaceInt()
	->getParentParamSurface(domain);
    DEBUG_ERROR_IF(curve.get() == 0 || surf.get() == 0,
		   "Inconsistence in the data structure");

    const double TOL = 1.0e-8;
    double seed[3], minpar[3], maxpar[3];
    if (cv_idx_ == 0)
    {
	seed[0] = guess[2];
	seed[1] = guess[0];
	seed[2] = guess[1];
    }
    else
    {
	for (int ki=0; ki<3; ki++)
	    seed[ki] = guess[ki];
    }
    minpar[0] = domain.umin();
    minpar[1] = domain.vmin();
    minpar[2] = start;
    maxpar[0] = domain.umax();
    maxpar[1] = domain.vmax();
    maxpar[2] = end;

    // Iterate
    CrvSrfDistFun distfun(surf.get(), curve.get(), minpar, maxpar);

    FunctionMinimizer<CrvSrfDistFun> funmin(3, distfun, seed, TOL);
    minimise_conjugated_gradient(funmin);//, 3); // number of iterations in each cycle

    if (cv_idx_ == 0)
    {
	par[0] = funmin.getPar(2);
	par[1] = funmin.getPar(0);
	par[2] = funmin.getPar(1);
     }
    else
    {
	for (int ki=0; ki<3; ki++)
	    par[ki] = funmin.getPar(ki);
   }
    dist = sqrt(funmin.fval());
}


//===========================================================================
void SfCvIntersector::postIterateBd()
//===========================================================================
{
    // Get all intersection points
    vector<shared_ptr<IntersectionPoint> > ints;
    int_results_->getIntersectionPoints(ints);
    double tol = 0.001*epsge_->getEpsge();

    double seed[3], param[3];
    double dist;
    int npar = 3;
    for (size_t ki=0; ki<ints.size(); ki++)
    {
	if (ints[ki]->getDist() < tol)
	    continue;  // No point in post iteration

	// Let the iteration seed be equal to the current parameter value of
	// the intersection point
	for (int kj=0; kj<npar; kj++)
	    seed[kj] = ints[ki]->getPar(kj);

	doIterate2(param, dist, seed);
	if (dist < ints[ki]->getDist())
	    ints[ki]->replaceParameter(param) ; //, npar);
    }
	
}

//===========================================================================
double SfCvIntersector::distInCandidatePar(double par, int dir, const double* seed)
//===========================================================================
{
    double start, end;  // The parametric limitation of the current
                        // curve in the parameter interval of the
                        // parent curve
    RectDomain domain;  // Similar for the surface
    shared_ptr<ParamCurve> curve = 
	obj_int_[cv_idx_]->getParamCurveInt()->getParentParamCurve(start, end);
    shared_ptr<ParamSurface> surf = 
	obj_int_[sf_idx_]->getParamSurfaceInt()->getParentParamSurface(domain);
    DEBUG_ERROR_IF(curve.get() == 0 || surf.get() == 0,
		   "Inconsistence in the data structure");
    double guess[2];
    int ki, kj;
    for (ki=0, kj=0; ki<3; ki++)
    {
	if (ki == dir)
	    continue;
	guess[kj++] = seed[ki];
    }

    if (dir == 2*cv_idx_)
    {
	// Evaluate the curve and iterate to the closest point in the surface
	Point pt;
	curve->point(pt, par);

	double clo_u, clo_v, clo_dist;
	Point clo_pt;
	surf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsge_->getEpsge(), &domain, guess);

	return clo_dist;
    }
    else
    {
	// Fetch the constant parameter curve corresponding to the candidate
	// subdivision parameter in the surface and iterate to the closest point
	// between the two curves
	int pdir = (cv_idx_ == 0) ? dir-1 : dir;
	shared_ptr<ParamCurve> curve2 = obj_int_[sf_idx_]->getParamSurfaceInt()->
	    getConstantParameterCurve(pdir, par);
	double start2 = (pdir == 0) ? domain.umin() : domain.vmin();
	double end2 = (pdir == 0) ? domain.umax() : domain.vmax();
	double par1, par2, dist;
	Point ptc1, ptc2;
	ClosestPoint::closestPtCurves(curve.get(), curve2.get(), start, end, start2, end2, 
			(cv_idx_ == 0) ? guess[0] : guess[1], (cv_idx_ == 0) ? guess[1] : guess[0],
			par1, par2, dist, ptc1, ptc2);

	return dist;
    }
}

//===========================================================================
int SfCvIntersector::simpleCase()
//===========================================================================
{
    // Purpose: Do simple case test between one surface and one curve
    // Currently only testing between directional cones is performed

    // Get the direction cones.
    // Should the complete cones be computed if they are larger than
    // pi?
    DirectionCone cone1, cone2;
    try {
	cone1 = obj_int_[0]->directionCone();
	cone2 = obj_int_[1]->directionCone();
    } catch (...) {
	return 0;  // Can't accept this kind of degeneracy as simple
		   // case
    }

    // Check if any of the cones are greater than pi. In that case,
    // it is no simple case.
    // NB! It might be that this test should be replaced by something
    // more strict in the future.
    if (cone1.greaterThanPi() || cone2.greaterThanPi())
	return 0;

    if (cone1.perpendicularOverlaps(cone2)) {
	double overlap;
	Point axis1 = cone1.centre();
	Point axis2 = cone2.centre();
	double angle = axis1.angle_smallest(axis2);
	double angle2 = std::min(angle, fabs(M_PI-angle));
	double sf_ang = (sf_idx_ == 0) ? cone1.angle() : cone2.angle();
	double cv_ang = (cv_idx_ == 0) ? cone1.angle() : cone2.angle();
	overlap = std::max(angle + sf_ang + cv_ang - 0.5*M_PI,
			   sf_ang + cv_ang - angle + 0.5*M_PI);


	// Remember statistical information
	if (!hasComplexityInfo())
	    complexity_info_
		= shared_ptr<ComplexityInfo>(new ComplexityInfo());
	complexity_info_->setConeOverlap(overlap);

	// Extra test if there is a small overlap between the cones. 
	if (angle < M_PI - epsge_->getAngleTol() &&
	    angle > epsge_->getAngleTol() &&
	    angle2 < 0.8*(0.5*M_PI - sf_ang) &&
	    cv_ang < 0.8*(0.5*M_PI - sf_ang)) {
	    int simple_case =  simpleCase2(axis1, axis2);
	    return simple_case;
	}
  
	return 0;
    }

    // Otherwise it is a simple case
    return 1;
}


//===========================================================================
int SfCvIntersector::performRotatedBoxTest(double eps1, double eps2)
//===========================================================================
{
    // Purpose : Perform a rotated box test between a curve and a
    // surface

//     static int nmb_sep = 0;
//     static int nmb_sep2 = 0;
//     int nmbpar = numParams();
//     int rec = nmbRecursions();
//     bool sep = false;
//     double ang_tol = epsge_->getAngleTol();

    ParamCurveInt *curve = obj_int_[cv_idx_]->getParamCurveInt();
    ParamSurfaceInt *surf = obj_int_[sf_idx_]->getParamSurfaceInt();
    DEBUG_ERROR_IF(curve == 0 || surf == 0,
		   "Inconsistence in the data structure");

    // Only done for simple surfaces
    if (!surf->isSimple())
	return 1;

    // Define coordinate system for the rotation. First check if an
    // intersection point between the geometry objects is found
    vector<Point> axis(2), der(3);
    Point norm;
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);
    //int nmb_rec = nmbRecursions();
    if (int_pts.size() > 0)
    {
	// Make coordinate system from the partial derivatives of the surface
	// in the (first) intersection point
	const double *param = (sf_idx_ == 0) ? int_pts[0]->getPar1() : 
	    int_pts[0]->getPar2();
	obj_int_[sf_idx_]->point(der, param, 1);
	axis[0] = der[1];
	norm = der[1].cross(der[2]);
//	axis[1] = norm.cross(axis[0]);
	axis[1] = norm;
	if (getenv("DEBUG_ROTATEDBOX") && *(getenv("DEBUG_ROTATEDBOX")) == '1')
	    std::cout << "Intersection point found " << std::endl;
	
    }
    else
    {
	curve->axisFromEndpts(axis[0]);
	surf->axisFromCorners(der[0], der[1]);
	norm = der[0].cross(der[1]);
//	axis[1] = norm.cross(axis[0]);
	axis[1] = norm;
	if (getenv("DEBUG_ROTATEDBOX") && *(getenv("DEBUG_ROTATEDBOX")) == '1')
	    std::cout << "No intersection point " << std::endl;
	
    }
  
    // Make rotated boxes
    if (axis[0].length() < epsge_->getNumericalTol() ||
	axis[1].length() < epsge_->getNumericalTol()) 
	return 1;
    axis[0].normalize();
    axis[1].normalize();
    if (axis[0].angle_smallest(axis[1]) < epsge_->getNumericalTol())
	return 1;

    RotatedBox box1 = surf->getRotatedBox(axis);
    RotatedBox box2 = curve->getRotatedBox(axis);
    bool overlap;
    overlap = box1.overlaps(box2, eps1, 0.0);
	// TESTING
//     if (overlap)
//     {
// 	    if (!box1.overlaps(box2, eps1, -eps1)) 
// 	    {
// 		nmb_sep++;
// 		sep = true;
// 	    }
// 	    if  (!box1.overlaps(box2, 0.0, 0.0)) 
// 	    {
// 		nmb_sep2++;
// 		sep = true;
// 	    }
// 	    if (sep)
// 	    {
// 		std::cout << " Rotated separation " << nmb_sep;
// 		std::cout << ", no fat: " << nmb_sep2 << ", recursion: ";
// 		std::cout << rec << std::endl;
// 	    }
//     }

    if (overlap)
    {
	overlap = box1.overlaps(box2, eps1, eps2);
	if (!overlap)
	{
	    if (foundIntersectionNearBoundary())
		return 1;
	}
    }
    if (overlap)
    {
	std::swap(axis[0], axis[1]);
	box1 = surf->getRotatedBox(axis);
	box2 = curve->getRotatedBox(axis);
	overlap = box1.overlaps(box2, eps1, 0.0);
//     if (overlap)
//     {
// 	    if (!box1.overlaps(box2, eps1, -eps1)) 
// 	    {
// 		nmb_sep++;
// 		sep = true;
// 	    }
// 	    if  (!box1.overlaps(box2, 0.0, 0.0)) 
// 	    {
// 		nmb_sep2++;
// 		sep = true;
// 	    }
// 	    if (sep)
// 	    {
// 		std::cout << " Rotated separation " << nmb_sep;
// 		std::cout << ", no fat: " << nmb_sep2 << ", recursion: ";
// 		std::cout << rec << std::endl;
// 	    }
//     }
	if (overlap)
	{
	    overlap = box1.overlaps(box2, eps1, eps2);
	    if (!overlap)
	    {
		if (foundIntersectionNearBoundary())
		    return 1;
	    }
	}
    }
    return (overlap) ? 1 : 0;
}


//===========================================================================
bool SfCvIntersector::foundIntersectionNearBoundary()
//===========================================================================
{
    // Purpose: Check for an intersection point close to one of the
    // endpoints of the curve provided that no intersection point is
    // found at this endpoint already

    double ptol = epsge_->getRelParRes();

    // Fetch parameter interval
    double ta, tb;
    ta = obj_int_[cv_idx_]->startParam(0);
    tb = obj_int_[cv_idx_]->endParam(0);

    // Fetch already existing intersection points
    bool has_int[2];
    has_int[0] = has_int[1] = false;
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    double par;
    for (size_t ki=0; ki<int_pts.size(); ki++)
    {
	par = int_pts[ki]->getPar(cv_idx_);
	if (fabs(par-ta) < ptol)
	    has_int[0] = true;
	if (fabs(tb-par) < ptol)
	    has_int[1] = true;
    }

    double param[3], guess[3], dist;
    guess[sf_idx_] = 0.5*(obj_int_[sf_idx_]->startParam(0) +
			    obj_int_[sf_idx_]->endParam(0));
    guess[sf_idx_+1] = 0.5*(obj_int_[sf_idx_]->startParam(1) +
			    obj_int_[sf_idx_]->endParam(1));
    if (has_int[0]);
    else
    {
	guess[2*cv_idx_] = ta;

	    // Iterate to intersection point
	    doIterate(param, dist, guess);
	    if (dist <= epsge_->getEpsge())
		return true;
	}
    if (has_int[1]);
    else
    {
	guess[2*cv_idx_] = tb;

	// Iterate to intersection point
	doIterate(param, dist, guess);
	if (dist <= epsge_->getEpsge())
	    return true;
    }

    return false;
}


//===========================================================================
int SfCvIntersector::simpleCase2(Point& axis1, Point& axis2)
//===========================================================================
{
    // Purpose: Do extra simple case test between two parametric curve
    // and one surface. Based on the sisl function s1797.

    // The function only works for spline geometries. Check if we have
    // the correct input.
    int kr;
    for (kr=0; kr<2; kr++)
    {
	if (!(obj_int_[kr]->isSpline()))
	    return 0;
    }

    double angle = axis1.angle_smallest(axis2);
    if (angle > 0.5*M_PI)
	angle = M_PI - angle;

    double tangle[2];
    double tlen = axis1*axis2;
    Point scen1, scen2;

    // Compute the cone angles
    for (kr=0,scen1=axis1, scen2=axis2-tlen*axis1; 
	 kr<2; kr++, scen1=axis2, scen2=axis1-tlen*axis2)
    {
	scen2.normalize();
	tangle[kr] = obj_int_[kr]->getOptimizedConeAngle(scen1, scen2);
    }

    // Performing a simple case check. 
    if (tangle[0] + tangle[1] < 0.5*M_PI - angle + epsge_->getRelParRes())
	return 1;  // Simple case
    else
	return 0;
}


//===========================================================================
int SfCvIntersector::linearCase()
//===========================================================================
{
    // Purpose : Given linear parametric curve and one surface,
    // compute the intersection.

    int ki, kj;

    // First check existance of any intersection points
    int dir = (cv_idx_ == 0) ? 0 : 2;  // Parameter direction of the curve
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    if (int_pts.size() == 1)
	return 1;   // One intersection point already computed

    if (int_pts.size() > 1)
    {
	shared_ptr<IntersectionPoint> dum;
	// Sort according to the parameter direction of the curve
	for (ki=0; ki<int(int_pts.size()); ki++)
	    for (kj=ki+1; kj<int(int_pts.size()); kj++)
		if (int_pts[kj]->getPar(dir) < int_pts[ki]->getPar(dir))
		{
		    dum = int_pts[kj];
		    int_pts[kj] = int_pts[ki];
		    int_pts[ki] = dum;
		}

	// Make sure that all intersection points are connected
	for (ki=1; ki<int(int_pts.size()); ki++)
	    if (!(int_pts[ki-1]->isConnectedTo(int_pts[ki])))
		int_pts[ki-1]->connectTo(int_pts[ki], LINEAR_SFCV);
	return 1;
    }

    // Check if any intersections can occur
    // First compute the "midpoint" of the surface and fetch the mid surface normal
    // from the normal cone.
    double par[2];
    for (ki=0; ki<obj_int_[sf_idx_]->numParams(); ki++)
	par[ki] = 0.5*(obj_int_[sf_idx_]->startParam(ki) + obj_int_[sf_idx_]->endParam(ki));
    Point midpt;
    obj_int_[sf_idx_]->point(midpt, par);

    DirectionCone cone = obj_int_[sf_idx_]->directionCone();
    Point normal = cone.centre();

    // Fetch the endpoints of the curve
    Point pt1, pt2;
    double ta = obj_int_[cv_idx_]->startParam(0);
    double tb = obj_int_[cv_idx_]->endParam(0);
    obj_int_[cv_idx_]->point(pt1, &ta);
    obj_int_[cv_idx_]->point(pt2, &tb);

    // Check if the endpoints lie on the same side of the surface
    //double sc1 = (midpt - pt1)*normal;
    //double sc2 = (midpt - pt2)*normal;
//     if (sc1*sc2 > 0.0)
// 	return 0;   // No intersection

    // Iterate to the intersection point
    double param[3], dist;
    doIterate(param, dist);
 
  if (dist <= epsge_->getEpsge())
    {
      // An intersection point is found. Represent it in the data
      // structure

	shared_ptr<IntersectionPoint> tmp = 
	    int_results_->addIntersectionPoint(obj_int_[0], 
					 obj_int_[1], 
					 getTolerance(),
					 param, 
					 param+obj_int_[0]->numParams());
      return 1; // A new point is found
    }

  return 0; 
}


//===========================================================================
int SfCvIntersector::doSubdivide()
//===========================================================================
{
    // Purpose: Given parametric curve and one surface and no simple
    // case situation, subdivide the curves to produce more simple sub
    // problems. Perform intersection with the subdivision points.
    // Prepare for the next recursion level.

    int perm[3];  // Three parameter directions that need sorting
    int nmb_subdiv;  // Number of parameter directions in which to subdivide
    int nmbdir1 = obj_int_[0]->numParams();
    int ki, kj, kr;

//     if (getenv("SUBDIV_SFCV") && *(getenv("SUBDIV_SFCV")) == '1') {
// 	cout << "================================================" << endl;
// 	cout << "Domain 1: ";
// 	for (ki=0; ki<obj_int_[0]->numParams(); ki++) {
// 	    cout << obj_int_[0]->startParam(ki) << " ";
// 	    cout << obj_int_[0]->endParam(ki) << " ";
// 	}
// 	cout << endl;
// 	cout << "Domain 2: ";
// 	for (ki=0; ki<obj_int_[1]->numParams(); ki++) {
// 	    cout << obj_int_[1]->startParam(ki) << " ";
// 	    cout << obj_int_[1]->endParam(ki) << " ";
// 	}
// 	cout << endl;
//     }

    // Sort the parameter directions according to importance of
    // subdivision according to properties of the objects and already
    // computed intersections
    nmb_subdiv = sortParameterDirections(perm);

    if (nmb_subdiv == 0)
	return 0;   // Not possible to subdivide neither the curve nor
                    // the surface

    // Intersection objects belonging to the next recursion level.
    // First they point to the objects at this level. Then the current
    // objects are replaced as subdivision is taking place.
    int nbobj[2];
    vector<shared_ptr<ParamGeomInt> > sub_objects;
    for (ki=0; ki<2; ki++) {
	sub_objects.push_back(obj_int_[ki]);
    }
    nbobj[0] = nbobj[1] = 1; // Currently one instance of each //
                             // intersection object at this level

    // For each parameter direction in prioritized order, fetch an
    // appropriate subdivision parameter, and perform subdivision
    double subdiv_par;
    SubdivisionClassification found;
    vector<shared_ptr<ParamGeomInt> > obj_sub;
    vector<shared_ptr<ParamGeomInt> > subdivobj;

    for (ki=0; ki<nmb_subdiv; ki++) {
	found = getSubdivisionParameter(perm[ki], subdiv_par);
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

	if (getenv("SUBDIV_SFCV") && *(getenv("SUBDIV_SFCV")) == '1') {
	    cout << "Subdivide dir = " << perm[ki]
		 << " par = " << subdiv_par;
	    cout << " criterium = " << found << endl;
	}
	// Subdivide the current object (or objects. If the current
	// object is the surface then it can have been subdivided in
	// one parameter direction already and we have two sub
	// surfaces to subdivide)
	int idxobj = (perm[ki] >= nmbdir1);
	int idx = idxobj*nbobj[0];
	for (kj=0; kj<nbobj[idxobj]; kj++) {
	    obj_sub.clear();
	    subdivobj.clear();
	    try {
		sub_objects[idx+kj]->subdivide(perm[ki]-idxobj*nmbdir1, 
					       subdiv_par, obj_sub, subdivobj);
	    } catch (...) {
		obj_sub.clear();
		subdivobj.clear();
	    }
	    if (obj_sub.size() < 1 || subdivobj.size() == 0)
		continue;  // No new objects 

	    for (kr=0; kr<int(obj_sub.size()); kr++)
		obj_sub[kr]->setParent(obj_int_[idxobj].get());
      
	    for (kr = 0; kr<int(subdivobj.size()); kr++) {
		// Intersect the subdivision object with the other object
		// @@@ VSK Faktorenes orden er ikke likegyldig!!
		shared_ptr<Intersector> subdiv_intersector = 
		    lowerOrderIntersector((idxobj==0) ? subdivobj[kr] : 
					  obj_int_[0], 
					  (idxobj==0) ? obj_int_[1] :
					  subdivobj[kr], 
					  this, perm[ki], subdiv_par);

		// Is it here relevant to fetch existing intersection
		// points and/or insert points into intersection
		// curves before computing intersection points with
		// the subdivision points. Normally, this is mostly of
		// interest when surfaces are involved, but
		// intersection intervals might exist?
		// These computations do anyway involve the
		// intersection pool, but they might be trigged
		// somehow. The parameter direction and value are
		// required information

		if (found == DIVIDE_SING && hasSingularityInfo()) {
		    // Transfer singularity information to the child
		    // process
		    subdiv_intersector->setSingularityInfo
			(getSingularityInfo(), perm[ki]);
		}

		if (getenv("SUBDIV_SFCV") && *(getenv("SUBDIV_SFCV")) == '1') {
		    vector<shared_ptr<IntersectionPoint> > ipoint1;
		    subdiv_intersector->getIntPool()->getIntersectionPoints
			(ipoint1);
		    cout << "Intersection points prior : " << endl;
		    for (size_t k1 = 0; k1 < ipoint1.size(); k1++) {
			int ctr = 0;
			for (int k2 = 0; k2 < 2; k2++) {
			    if (ctr == perm[ki]) {
				cout << subdiv_par << "  ";
			    }
			    cout << ipoint1[k1]->getPar(k2) << "  ";
			    ++ctr;
			}
			if (perm[ki] == 2) {
			    cout << subdiv_par << "  ";
			}
			cout << "  " << ipoint1[k1]->getDist();
			cout << endl;
		    }
		}

		subdiv_intersector->compute(false);

		if (getenv("SUBDIV_SFCV") && *(getenv("SUBDIV_SFCV")) == '1') {
		    vector<shared_ptr<IntersectionPoint> > ipoint1;
		    subdiv_intersector->getIntPool()->getIntersectionPoints
			(ipoint1);
		    cout << "Intersection points post : " << endl;
		    for (size_t k1 = 0; k1 < ipoint1.size(); k1++) {
			int ctr = 0;
			for (int k2 = 0; k2 < 2; k2++) {
			    if (ctr == perm[ki]) {
				cout << subdiv_par << "  ";
			    }
			    cout << ipoint1[k1]->getPar(k2) << "  ";
			    ++ctr;
			}
			if (perm[ki] == 2) {
			    cout << subdiv_par << "  ";
			}
			cout << "  " << ipoint1[k1]->getDist();
			cout << endl;
		    }
		}
		// Check quality of intersection points.
		// If the quality is not sufficient, find a new
		// subdivision parameter and repeat the process.
		// Otherwise...
		int nmb_orig = int_results_->numIntersectionPoints();
		int_results_->includeReducedInts
		    (subdiv_intersector->getIntPool());

		// Post iterate new points
		if (found != DIVIDE_SING && 
		    nmb_orig < int_results_->numIntersectionPoints()) {
		    postIterate(nmb_orig, perm[ki]);
		}
	    }

	    // Replace pointers to current objects with pointers to
	    // the sub objects after subdivision. The new objects can
	    // be the final sub objects or they can be subdivided once
	    // more depending on the number of parameter directions
	    // elected for subdivision in the current object.
	    sub_objects[idx+kj] = obj_sub[0];
	    sub_objects.insert(sub_objects.begin()+idx+kj+1,
			       obj_sub.begin()+1,
			       obj_sub.end());
	    nbobj[idxobj] += ((int)obj_sub.size()-1);
	    kj += ((int)obj_sub.size()-1);

	}
    }

    if (nmb_subdiv == 0)
	return 0;  // Not subdivided

    // Create new intersector objects
    for (ki=0; ki<nbobj[0]; ki++) {
	for (kj=0; kj<nbobj[1]; kj++) {
	    shared_ptr<Intersector> intersector = 
		shared_ptr<Intersector>
		(new SfCvIntersector(sub_objects[ki],
				     sub_objects[nbobj[0]+kj],
				     epsge_, this));
	    sub_intersectors_.push_back(intersector);
	}
    }

    return 1;
}


//===========================================================================
void SfCvIntersector::postIterate(int nmb_orig, int dir, bool keep_endpt)
//===========================================================================
{
    // Purpose : Post itererate new intersection points

    int ki, kj, kr;
    int npar = 3;
    double param[3], seed[3], dist;
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);
    double ptol = epsge_->getRelParRes();
    
    // Fetch parameter domain
    double ta[3], tb[3];
    for (ki=0, kr=0; ki<2; ki++)
    {
	int kn = obj_int_[ki]->numParams();
	for (kj=0; kj<kn; kj++, kr++)
	{
	    ta[kr] = obj_int_[ki]->startParam(kj);
	    tb[kr] = obj_int_[ki]->endParam(kj);
	}
    }

    for (ki=nmb_orig; ki < (int)int_pts.size(); ki++)
    {
	// Let the iteration seed be equal to the current parameter value of
	// the intersection point
	for (kj=0; kj<npar; kj++)
	    seed[kj] = int_pts[ki]->getPar(kj);

	// Check first if a boundary intersection may have been extended
	// along the a subdivision curve into the inner of the current objects
	// and the quality of this extension is poor
	for (kj=0; kj<npar; kj++)
	{
	    if (seed[kj]-ta[kj] < ptol || tb[kj]-seed[kj]<ptol)
		break;  // Boundary point
	}
	if (kj == npar)
	{
	    // A point in the inner is found. Check the number of neighbours
	    if (int_pts[ki]->numNeighbours() == 1)
	    {
		// Make sure that the neigbouring point is contained in the
		// current pool and check the quality of the intersection points
		vector<IntersectionPoint*> neighbour;
		int_pts[ki]->getNeighbours(neighbour);
		for (kr=0; kr<nmb_orig; kr++)
		    if (neighbour[0] == int_pts[kr].get())
			break;
		if (kr < nmb_orig)
		{
		    // Point found. Check quality
		    if (neighbour[0]->getDist() <= int_pts[ki]->getDist())
		    {
			// Remove new point
			int_results_->removeIntPoint(int_pts[ki]);
			int_pts.erase(int_pts.begin()+ki);
			ki = nmb_orig - 1;   // Check all new intersection points again
			continue;
		    }
		}
	    }
	}

	// The point is still present. Try to iterate to a better intersection point
	doIterate(param, dist, seed);
	double curr_par = (dir < 0) ? 0.0 : int_pts[ki]->getPar(dir);
	if (dist < int_pts[ki]->getDist() && 
	    (dir < 0 || fabs(param[dir]-curr_par) < epsge_->getRelParRes()))
	{
	    // Check if the parameter of the current intersection point
	    // can be replaced
	    for (kj=0; kj<npar; kj++)
		if (!(int_pts[ki]->inInfluenceArea(kj, param[kj], false)) ||
		    (fabs(seed[kj]-ta[kj]) < epsge_->getRelParRes() ||
		     fabs(tb[kj]-seed[kj]) < epsge_->getRelParRes()))
		    break;

	    if (kj == npar)
	    {
		// The new position lies within the influence area of the
		// old point and the position is improved. Check if it
		// exists already
		shared_ptr<IntersectionPoint> closest;
		bool exist = int_results_->closestInDomain(param, closest);
		if (closest.get() == int_pts[ki].get())
		    exist = false;

		if (exist)
		{
		    // Check influence area
		    for (kr=0; kr<3; kr++) {
			if (closest->inInfluenceArea(kr, param[kr], false)) 
			    break;
			}
		    if (kr == 3)
			exist = false;  
		    
		}
		
		if (exist)
		{
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfCvIntersector::postIterate. Removed point at: ";
			std::cout << seed[0] << " " << seed[1] << " " << seed[2] << ", dist: ";
			std::cout << int_pts[ki]->getDist() << std::endl;
		    }
		    int_results_->removeIntPoint(int_pts[ki]);
		}
		else
		{
		    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
		    {
			std::cout << "SfCvIntersector::postIterate. Moved point at: ";
			std::cout << seed[0] << " " << seed[1] << " " << seed[2] << ", dist: ";
			std::cout << int_pts[ki]->getDist() << " to ";
			std::cout << param[0] << " " << param[1] << " " << param[2] << ", dist: ";
			std::cout << dist << std::endl;
		    }
		    int_pts[ki]->replaceParameter(param) ; //, npar);
		}
	    }
	}
    }
}


//===========================================================================
int SfCvIntersector::sortParameterDirections(int perm[])
//===========================================================================
{
    // Purpose : Fetch information related to the objects in order to
    // decide which curves that can be subdivided and which parameter
    // direction is most important to subdivide.

  // @@@ VSK. This routine is extremely similar to the one for curve-curve
  // intersections. Should it be moved one level up in the hierarchy

  double length[3], wiggle[3];  // Estimated length and wiggliness
  bool inner_knots[3], critical_val[3], can_divide[3];
  bool has_inner_ints[3];
  int deg_edge[3];
  double rel_wiggle = 0.1;
  double rel_length = 0.1;
  int nmbdir[2];
  int ki, kj, kr;

  nmbdir[0] = obj_int_[0]->numParams();
  nmbdir[1] = obj_int_[1]->numParams();

  // Fetch information from the curve and surface
  obj_int_[0]->getLengthAndWiggle(length, wiggle);
  obj_int_[1]->getLengthAndWiggle(length+nmbdir[0], wiggle+nmbdir[0]);

  double med_len = 0.0, med_wiggle = 0.0;
  for (ki=0; ki<3; ki++)
  {
      med_len += length[ki];
      med_wiggle += wiggle[ki];
  }
  med_len /= (double)3;
  med_wiggle /= (double)3;

  double min_length = std::max(0.5*med_len, epsge_->getEpsge());
  double min_wiggle = std::max(0.5*med_wiggle, 0.02);

  for (ki=0, kr=0; ki<2; ki++)
  {
    for (kj=0; kj<nmbdir[ki]; kj++, kr++)
      {
	inner_knots[kr] = obj_int_[ki]->hasInnerKnots(kj);

	critical_val[kr] = obj_int_[ki]->hasCriticalVals(kj);

	can_divide[kr] = obj_int_[ki]->canDivide(kj);

	has_inner_ints[kr] = int_results_->hasPointsInInner(kr);

	deg_edge[kr] = obj_int_[ki]->isDegenerate(epsge_->getEpsge(), kj);

	if (deg_edge[kr])
	{
	    if (deg_edge[kr] == 1 || deg_edge[kr] == 3)
	    {
		double par = obj_int_[ki]->startParam(kj);
		if (!(int_results_->existIntersectionPoint(kr, par)))
		    deg_edge[kr] -= 1;
	    }
	    if (deg_edge[kr] == 2 || deg_edge[kr] == 3)
	    {
		double par = obj_int_[ki]->endParam(kj);
		if (!(int_results_->existIntersectionPoint(kr, par)))
		    deg_edge[kr] -= 2;
	    }
	}
      }
    // Avoid subdivison across a degenerate edge
    if (nmbdir[ki] == 2)
    {
	if (deg_edge[kr-2] && can_divide[kr-2] && !deg_edge[kr-1])
	    can_divide[kr-1] = false;
	if (deg_edge[kr-1] && can_divide[kr-1] && !deg_edge[kr-2])
	    can_divide[kr-2] = false;
    }
  }

  

  // Fetch information from the intersection pool according to inner 
  // intersection points

  // Number of parameter directions is three.
  int size = 3;

  int curr = 0;
  int min_nmb = 0;

  // Initiate permutation array
  for (ki=0; ki<size; ki++)
    perm[ki] = ki;

  // Sort according to the given values
  // First sort out the directions where subdivsion is impossible
  for (ki=0; ki<size; ki++)
    {
      if (!can_divide[perm[ki]])
	{
	  if (perm[ki] < size-1)
	    std::swap(perm[ki], perm[size-1]);
	  size--;
	}
    }

  // First priority is parameter directions with critical values
  for (ki=curr; ki<size; ki++)
    for (kj=ki+1; kj<size; kj++)
      if (critical_val[perm[kj]] && !critical_val[perm[ki]])
	{
	  std::swap(perm[ki], perm[kj]);
	}
  for (ki=curr; ki<size-1; ki++)
      if (critical_val[perm[ki]])
	  curr++;

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
  for (ki=curr; ki<size-1; ki++)
      if (has_inner_ints[perm[ki]])
	  curr++;

  // Finally length and curvature
  min_nmb = curr;
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
  for (ki=size-1; ki>=min_nmb; ki--)
      if (length[perm[ki]] < min_length && wiggle[perm[ki]] < min_wiggle)
	  size--;

  return size;
}


//===========================================================================
SubdivisionClassification SfCvIntersector::getSubdivisionParameter(int dir, double& par)
//===========================================================================
{
  // @@@ VSK. This routine is extremely similar to the one for curve-curve
  // intersections. Should it be moved one level up in the hierarchy

  int ki, kj, kr, sgn;
  double aeps = epsge_->getEpsge();
  double ptol = 100.0*epsge_->getRelParRes();
  double gtol = 100.0*aeps;
 int nmbdir1 = obj_int_[0]->numParams();
 //int nmbdir2 = obj_int_[1]->numParams();
  double pval;
  double frac = 0.2;
  double frac2 = 0.01;
  bool has_candidate = false;
  double candidate_div;
  int max_rec = 15;

  // Set pointer to the intersection object corresponding to the parameter
  // direction
  ParamGeomInt *obj = 
    (dir < nmbdir1) ? obj_int_[0].get() : obj_int_[1].get();
  int pdir = (dir < nmbdir1) ? dir : dir-nmbdir1;
  double ta = obj->startParam(pdir);
  double tb = obj->endParam(pdir);
  double ptol_init = ptol;
  ptol = std::max(0.1*(tb - ta), ptol);
  gtol = std::max(0.05, gtol);

  // Get critical parameters
  // @@@ Critical parameters of different priority? Sorting?
  vector<double> critical_pars = obj->getCriticalVals(pdir);

  int size = (int)critical_pars.size();
  int is_critical = 0;
  if (size > 0)
    {
      // Check suitability of the critical parameters
      for (ki=0; ki<size; ki++)
	{
	    par = critical_pars[ki];
	    vector<shared_ptr<IntersectionPoint> > ints_in_area;
	  is_critical = int_results_->inInfluenceArea(dir, par, ints_in_area);
	  if (is_critical == 1)
	  {
	      // Check if the candidate division points lie distant from the
	      // intersection points and not in a critical area near the boundaries
	      // of the influence area.
	      is_critical = 0;
	      for (kj=0; kj<int(ints_in_area.size()); kj++)
	      {
		  pval = ints_in_area[kj]->getPar(dir);
		  if (fabs(par - pval) < epsge_->getRelParRes());
		  else if (fabs(par - pval) < frac2*(tb - ta) && 
			   fabs(par-pval) > epsge_->getRelParRes())
		  {
		      is_critical = 1;
		      break;
		  }
		  else if (distInCandidatePar(pval, dir, 
					      ints_in_area[kj]->getPar1()) >= 0.9*epsge_->getEpsge())
		  {
		      // The subdivision parameter would give rise to a low quality
		      // intersection point
		      is_critical = 1;
		      if (!has_candidate)
		      {
			  has_candidate = true;
			  candidate_div = par;
		      }
		      break;
		  }
	      }
	      
	  }
	  if (is_critical == 0 || is_critical == 2)
	    {
	      return DIVIDE_CRITICAL;
	    }
	}
    }

  // Look for high priority singularity
  double dist = 2.0*gtol;
  if (singularity_info_.get())
  {
      SingularityClassification type = singularity_info_->hasHighPriSing();
      if (type != NO_SING)
      {
	  par = singularity_info_->getHighPriSing(dir);
	  dist = 0.0;   // Not exactly known, but use this point for subdivision
      }
  }
  if (dist < gtol && par > ta+ptol_init && par < tb-ptol_init)
  {
//       par = std::max(par, ta+ptol);
//       par = std::min(par, tb-ptol);

      // Check the parameter value
      /*vector<shared_ptr<IntersectionPoint> > int_pts;
      is_critical = int_results_->inInfluenceArea(dir, par, int_pts);
      if (is_critical == 1)
      {
	  // Check if the candidate division points lie distant from the
	  // intersection points and not in a critical area near the boundaries
	  // of the influence area.
	  is_critical = 0;
	  for (kj=0; kj<int(int_pts.size()); kj++)
	  {
	      pval = int_pts[kj]->getPar(dir);
	      if (fabs(par - pval) < frac*(tb - ta) && 
		  fabs(par-pval) > epsge_->getRelParRes())
	      {
		  is_critical = 1;
		  break;
	      }
	      else if (distInCandidatePar(pval, dir, 
					  int_pts[kj]->getPar1()) >= 0.9*epsge_->getEpsge())
	      {
		  // The subdivision parameter would give rise to a low quality
		  // intersection point
		  is_critical = 1;
		  if (!has_candidate)
		  {
		      has_candidate = true;
		      candidate_div = par;
		  }
		  break;
	      }
	      }
	      
      }
      if (is_critical == 0 || is_critical == 2)
      return DIVIDE_SING;*/
     
      return DIVIDE_SING;
  }



  // Look for a suitable knot in which to subdivide. First fetch
  // the knots sorted according to the distance to the mid-parameter
  vector<double> knot_vals = obj->getInnerKnotVals(pdir, true);
  size = (int)knot_vals.size();
   if (size > 0)
    {
      // Check suitability of the knots
      for (ki=0; ki<size; ki++)
	{
	    par = knot_vals[ki];
	    vector<shared_ptr<IntersectionPoint> > ints_in_area;
	  is_critical = int_results_->inInfluenceArea(dir, par, ints_in_area);
	  if (is_critical == 1)
	  {
	      // Check if the candidate division points lie distant from the
	      // intersection points and not in a critical area near the boundaries
	      // of the influence area.
	      is_critical = 0;
	      for (kj=0; kj<int(ints_in_area.size()); kj++)
	      {
		  pval = ints_in_area[kj]->getPar(dir);
		  if (fabs(par - pval) < frac*(tb - ta) && 
		      fabs(par-pval) > epsge_->getRelParRes())
		  {
		      is_critical = 1;
		      break;
		  }
		  else if (distInCandidatePar(pval, dir, 
					      ints_in_area[kj]->getPar1()) >= 0.9*epsge_->getEpsge())
		  {
		      // The subdivision parameter would give rise to a low quality
		      // intersection point
		      is_critical = 1;
		      if (!has_candidate)
		      {
			  has_candidate = true;
			  candidate_div = par;
		      }
		      break;
		  }
	      }
	      
	  }
	  if (is_critical == 0 || is_critical == 2)
	    {
	      return DIVIDE_KNOT;
	    }
	}
    }
 
   // Check for inner intersection points in which to subdivide
   vector<double> inner_ints = int_results_->getSortedInnerInts(dir);
   size = (int)inner_ints.size();
  if (size > 0)
    {
      // Check suitability of the intersection points
      for (ki=size/2, kj=1, sgn=-1; ki<size && ki>=0; 
	   ki+=sgn*kj, kj++, sgn*=-1)
	{
	    if (inner_ints[ki] < ta+frac2*(tb-ta) ||
		inner_ints[ki] > tb-frac2*(tb-ta))
		continue;

	    vector<shared_ptr<IntersectionPoint> > int_pts;
	    par = inner_ints[ki];
	  is_critical = int_results_->inInfluenceArea(dir, par, int_pts);
	  if (is_critical == 1)
	  {
	      // Check if the candidate division points lie distant from the
	      // intersection points and not in a critical area near the boundaries
	      // of the influence area.
	      is_critical = 0;
	      for (kr=0; kr<int(int_pts.size()); kr++)
	      {
		  pval = int_pts[kr]->getPar(dir);
		  if (fabs(par - pval) < frac*(tb - ta) && 
		      fabs(par-pval) > epsge_->getRelParRes())
		  {
		      is_critical = 1;
		      break;
		  }
		  else if (distInCandidatePar(pval, dir, 
					      int_pts[kr]->getPar1()) >= 0.9*epsge_->getEpsge())
		  {
		      // The subdivision parameter would give rise to a low quality
		      // intersection point
		      is_critical = 1;
		      if (!has_candidate)
		      {
			  has_candidate = true;
			  candidate_div = par;
		      }
		      break;
		  }
	      }
	      
	  }
	  if (is_critical == 0 || is_critical == 2)
	  {
	      return DIVIDE_INT;
	  }
	}
    }

  // Iterate for a closest point in which to subdivide
  double param[3];
  if (singularity_info_.get() == 0)
    {
	// Make empty singularity info instance
	singularity_info_ = (shared_ptr<SingularityInfo>)(new SingularityInfo());
    }
  
  // Check if a closest point exist already

  if (singularity_info_->hasPoint())
  {
      par = singularity_info_->getParam(dir);
      dist = 0.0;   // Not exactly known, but use this point for subdivision
  }
  else if (singularity_info_->iterationDone() ||
	   (prev_intersector_ && prev_intersector_->hasSingularityInfo() && 
	    prev_intersector_->getSingularityInfo()->iterationDone() == true &&
	    prev_intersector_->getSingularityInfo()->iterationSucceed() == false))
  {
      // No point in iterating again
      dist =  2.0*gtol;
  }
  else
  {
      // Iterate for a closest point
      try {
	  doIterate(param, dist);
	  singularity_info_->addIterationCount();
	  if (dist > 0.9*aeps && dist <= aeps)
	      dist = 2.0*gtol;  // Avoid subdivision in a bad intersection point
	  if (dist < gtol && param[dir] > ta+ptol && param[dir] < tb-ptol)
	      singularity_info_->setSingularPoint(param, 3);
	  par = param[dir];
      } catch (...) {
	  dist = 2.0*gtol;
      }
  }

  if (dist < gtol && par > ta+ptol && par < tb-ptol)
  {
//       par = std::max(par, ta+ptol);
//       par = std::min(par, tb-ptol);

      // Check the parameter value
      vector<shared_ptr<IntersectionPoint> > int_pts;
      is_critical = int_results_->inInfluenceArea(dir, par, int_pts);
      if (is_critical == 1)
      {
	  // Check if the candidate division points lie distant from the
	  // intersection points and not in a critical area near the boundaries
	  // of the influence area.
	  is_critical = 0;
	  for (kj=0; kj<int(int_pts.size()); kj++)
	  {
	      pval = int_pts[kj]->getPar(dir);
	      if (fabs(par - pval) < frac*(tb - ta) && 
		  fabs(par-pval) > epsge_->getRelParRes())
	      {
		  is_critical = 1;
		  break;
	      }
	      else if (distInCandidatePar(pval, dir, 
					  int_pts[kj]->getPar1()) >= 0.9*epsge_->getEpsge())
	      {
		  // The subdivision parameter would give rise to a low quality
		  // intersection point
		  is_critical = 1;
		  if (!has_candidate)
		  {
		      has_candidate = true;
		      candidate_div = par;
		  }
		  break;
	      }
	  }
	      
      }
      if (is_critical == 0 || is_critical == 2)
	  return DIVIDE_SING;
     
  }


  
  // Subdivide at a suitable parameter
  double divpar = 0.5*(ta+tb);
  double del = 0.1*(tb-ta);
  double tint = del;
  for (ki=0, sgn=-1; ki<9; ki++, divpar+=sgn*tint, tint+=del, sgn*=-1)
    {
      vector<shared_ptr<IntersectionPoint> > int_pts;
      is_critical = int_results_->inInfluenceArea(dir, divpar, int_pts);
      if (is_critical == 1)
      {
	  // Check if the candidate division points lie distant from the
	  // intersection points and not in a critical area near the boundaries
	  // of the influence area.
	  is_critical = 0;
	  for (kj=0; kj<int(int_pts.size()); kj++)
	  {
	      pval = int_pts[kj]->getPar(dir);
	      if (fabs(divpar - pval) < frac*(tb - ta) && 
		  fabs(divpar-pval) > epsge_->getRelParRes())
	      {
		  is_critical = 1;
		  break;
	      }
	      else if (distInCandidatePar(pval, dir, 
					  int_pts[kj]->getPar1()) >= 0.9*epsge_->getEpsge())
	      {
		  // The subdivision parameter would give rise to a low quality
		  // intersection point
		  is_critical = 1;
		  if (nmbRecursions() < max_rec && !has_candidate)
		  {
		      has_candidate = true;
		      candidate_div = divpar;
		  }
		  break;
	      }
	  }
	      
      }
      if (is_critical == 0 || is_critical == 2)
	{
	  par = divpar;
	  return DIVIDE_PAR;
	}
    }

  if (has_candidate)
  {
      par = candidate_div;
      return DIVIDE_PAR;
  }
  else
      return CANNOT_SUBDIVIDE;  // No subdivision parameter found
}


//===========================================================================


} // namespace Go

namespace {

//===========================================================================
CrvSrfDistFun::CrvSrfDistFun(const ParamSurface* sf, 
		       const ParamCurve* cv,
		       const double* const minpar,
		       const double* const maxpar)
//===========================================================================
    : sf_(sf), cv_(cv), p1vec_(3), p2vec_(2)
{
    RectDomain dom = sf->containingDomain();
    if (!minpar) {
	minpar_[0] = dom.umin();
	minpar_[1] = dom.vmin();
	minpar_[2] = cv_->startparam();
    } else {
	minpar_[0] = minpar[0];
	minpar_[1] = minpar[1];
	minpar_[2] = minpar[2];
    }
    if (!maxpar) {
	maxpar_[0] = dom.umax();
	maxpar_[1] = dom.vmax();
	maxpar_[2] = cv_->endparam();
    } else {
	maxpar_[0] = maxpar[0];
	maxpar_[1] = maxpar[1];
	maxpar_[2] = maxpar[2];
    }
}

//===========================================================================    
double CrvSrfDistFun::operator()(const double* arg) const
//===========================================================================
{
    sf_->point(p1_, arg[0], arg[1]);
    cv_->point(p2_, arg[2]);
    return p1_.dist2(p2_);
}

//===========================================================================
double CrvSrfDistFun::grad(const double* arg, double* res) const
//===========================================================================
{
    sf_->point(p1vec_, arg[0], arg[1], 1);
    cv_->point(p2vec_, arg[2], 1);
    d_ = p1vec_[0] - p2vec_[0];
    
    res[0] = 2 * d_ * p1vec_[1];
    res[1] = 2 * d_ * p1vec_[2];
    res[2] = -2 * d_ * p2vec_[1];
    
    return d_.length2();
}

//===========================================================================
double CrvSrfDistFun::minPar(int pardir) const
//===========================================================================
{
    //ASSERT(pardir == 0 || pardir == 1);
    return minpar_[pardir];
}

//===========================================================================
double CrvSrfDistFun:: maxPar(int pardir) const
//===========================================================================
{
    //ASSERT(pardir == 0 || pardir == 1);
    return maxpar_[pardir];
}

};
