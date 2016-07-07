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

#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GoIntersections.h"
#include <algorithm>
#include <stdexcept>
#include <fstream>

//#define DEBUG

using namespace Go;
using std::vector;
using std::pair;


//===========================================================================
CurveBoundedDomain::~CurveBoundedDomain()
//===========================================================================
{
}


//===========================================================================
CurveBoundedDomain::
CurveBoundedDomain(vector<shared_ptr<CurveLoop> > loops)
//===========================================================================
{
  size_t i;
  for (i=0; i<loops.size(); i++)
    loops_.push_back(loops[i]);
}


//===========================================================================
CurveBoundedDomain::CurveBoundedDomain(shared_ptr<CurveLoop> ccw_loop)
//===========================================================================
{
  loops_.push_back(ccw_loop);
}


//===========================================================================
int CurveBoundedDomain::isInDomain2(const Array<double, 2>& pnt,
				    double tolerance) const
//===========================================================================
{

  // Boundary points are critical. Check first if the point lies at a boundary 
  if (isOnBoundary(pnt, tolerance))
    return 2;
  else if (isInDomain(pnt, tolerance))
    return 1;
  else
    return 0;
}

//===========================================================================
bool CurveBoundedDomain::isInDomain(const Array<double, 2>& pnt,
				      double tolerance) const
//===========================================================================
{

  // Boundary points are critical. Check first if the point lies at a boundary 
  if (isOnBoundary(pnt, tolerance))
    return true;

  int nmb_catches = 0;
  try {
    // Avoid constant parameter directions
      vector<pair<double, double> > inside;
      getInsideIntervals(3, pnt[0], pnt[1], tolerance, inside);
      // Check if the pnt is inside any of the parameter intervals
      // lying inside the domain.
      int nmbint = (int)inside.size();
      int ki;
      for (ki=0; ki<nmbint; ki++) {
	  if (inside[ki].first-tolerance <= pnt[0] &&
	      inside[ki].second+tolerance >= pnt[0]) {
	      return true;
	  }
      }
  } catch (...) { // It seems there were an odd number of intersections.
      // We try intersecting in the other direction.
//       MESSAGE("Unstable method, trying in the other direction.");
      ++nmb_catches;
  }

  try {
  // Intersect the domain with a constant curve in 1. parameter
  // direction
      vector<pair<double, double> > inside;
      getInsideIntervals(1, pnt[0], pnt[1], tolerance, inside);
      // Check if the pnt is inside any of the parameter intervals
      // lying inside the domain.
      int nmbint = (int)inside.size();
      int ki;
      for (ki=0; ki<nmbint; ki++) {
	  if (inside[ki].first-tolerance <= pnt[0] &&
	      inside[ki].second+tolerance >= pnt[0]) {
	      return true;
	  }
      }
  } catch (...) { // It seems there were an odd number of intersections.
      // We try intersecting in the other direction.
//       MESSAGE("Unstable method, trying in the other direction.");
      ++nmb_catches;
  }

  try {
      vector<pair<double, double> > inside;
      getInsideIntervals(2, pnt[0], pnt[1], tolerance, inside);
      // Check if the pnt is inside any of the parameter intervals
      // lying inside the domain.
      int nmbint = (int)inside.size();
      int ki;
      for (ki=0; ki<nmbint; ki++) {
	  if (inside[ki].first-tolerance <= pnt[1] && 
	      inside[ki].second+tolerance >= pnt[1]) {
	      return true;
	  }
      }
  } catch (...) {
      if (nmb_catches == 2) { // If last attempt was a success, we're satisfied.
	// Not boundary
	THROW("Failed to decide whether point was inside boundary.");
      }
  }

  // Not inside
  return false;
}

//===========================================================================
bool CurveBoundedDomain::isOnCorner(const Array<double, 2>& point,
				    double tolerance) const
//===========================================================================
{
  for (int ki=0; ki<(int)loops_.size(); ++ki)
    {
      int nmb_crvs = loops_[ki]->size();
      for (int kj=0; kj<nmb_crvs; ++kj)
	{
	  shared_ptr<ParamCurve> crv;
	  try {
	    crv = getParameterCurve(ki, kj);
	  }
	  catch (...)
	    {
	      continue;  // Cannot check this curve
	    }

	  Point pt1 = crv->point(crv->startparam());
	  if (pt1.dist(Point(point)) < tolerance)
	    return true;
	  Point pt2 = crv->point(crv->endparam());
	  if (pt2.dist(Point(point)) < tolerance)
	    return true;
	}
	  
    }
  return false;
}

//===========================================================================
bool CurveBoundedDomain::isOnBoundary(const Array<double, 2>& point,
					double tolerance) const
//===========================================================================
{
  // Intersect the point with the curves bounding the domain (2D)
  for (int ki=0; ki<(int)loops_.size(); ++ki)
    {
      int nmb_crvs = loops_[ki]->size();
      for (int kj=0; kj<nmb_crvs; ++kj)
	{
	  shared_ptr<ParamCurve> crv;
	  try {
	    crv = getParameterCurve(ki, kj);
	  }
	  catch (...)
	    {
	      continue;  // Cannot check this curve
	    }

	  vector<double> int_pars;
	  vector<pair<double,double> > int_crvs;
	  Point pnt(point[0],point[1]);
	  intersectCurvePoint(crv.get(), pnt, tolerance, 
			      int_pars, int_crvs);
	  if (int_pars.size() > 0 || int_crvs.size() > 0)
	    return true;  // Intersection with boundary curve found
	}
    }

  // Not at boundary
  return false;
} 


//===========================================================================
void CurveBoundedDomain::closestInDomain(const Array<double, 2>& pnt,
					   Array<double, 2>& clo_pt,
					   double tolerance) const
//===========================================================================
{
//     bool is_in_domain;
    try {
	if (isInDomain(pnt, tolerance)) {
	    clo_pt = pnt;
	    return;
	}
    } catch (...) {
// 	MESSAGE("Failed deciding whether point was in domain.");
    }

    double global_clo_dist = 1.0e10;
    Point global_clo_pt;
    Point local_clo_pt;
    double clo_t, clo_dist;
    shared_ptr<ParamCurve> pcurve;
    Point ppnt(pnt[0], pnt[1]);
    for (int i = 0; i < int(loops_.size()); ++i)
	for (int j = 0; j < loops_[i]->size(); ++j) {
	    pcurve = getParameterCurve(i, j);
// 	    pcurve->closestPoint(Point(pnt.begin(), pnt.end(), false),
// 				 pcurve->startparam(), pcurve->endparam(),
// 				 clo_t, local_clo_pt, clo_dist);
	    pcurve->closestPoint(ppnt,
				 pcurve->startparam(), pcurve->endparam(),
				 clo_t, local_clo_pt, clo_dist);
	    if (clo_dist < global_clo_dist) {
		global_clo_dist = clo_dist;
		global_clo_pt = local_clo_pt;
	    }
	}

    if (global_clo_pt.size() == 0) {
        THROW("Unexpected incident occured.");
    }
    else {
        clo_pt.setValue(global_clo_pt.begin());
    }

    return;
}

//===========================================================================
void CurveBoundedDomain::closestOnBoundary(const Array<double, 2>& pnt,
					     Array<double, 2>& clo_bd_pt,
					     double tolerance) const
//===========================================================================
{
    double global_clo_dist = 1.0e10;
    Point global_clo_bd_pt;
    Point local_clo_bd_pt;
    double clo_t, clo_dist;
    shared_ptr<ParamCurve> pcurve;
    Point ppnt(pnt[0], pnt[1]);
    for (int i = 0; i < int(loops_.size()); ++i)
	for (int j = 0; j < loops_[i]->size(); ++j) {
	    pcurve = getParameterCurve(i, j);
// 	    pcurve->closestPoint(Point(pnt.begin(), pnt.end(), false),
// 				 pcurve->startparam(), pcurve->endparam(),
// 				 clo_t, local_clo_bd_pt, clo_dist);
	    pcurve->closestPoint(ppnt,
				 pcurve->startparam(), pcurve->endparam(),
				 clo_t, local_clo_bd_pt, clo_dist);
	    if (clo_dist < global_clo_dist) {
		global_clo_dist = clo_dist;
		global_clo_bd_pt = local_clo_bd_pt;
	    }
	}

    if (global_clo_bd_pt.size() == 0) {
	THROW("Unexpected incident occured.");
    } else {
	clo_bd_pt.setValue(global_clo_bd_pt.begin());
    }

    return;
}

//===========================================================================
void CurveBoundedDomain::getInternalPoint(double& upar, double& vpar) const
//===========================================================================
{
  // Make initial guess
  RectDomain dom = containingDomain();
  upar = 0.5*(dom.umin() + dom.umax());
  vpar = 0.5*(dom.vmin() + dom.vmax());

  // Adjust
  double tolerance = 1.0e-4; 
  for (int ki=1; ki<2; ++ki)
    {
      bool succeeded = true;
      vector<pair<double, double> > inside;
      try {
      getInsideIntervals(1, upar, vpar, tolerance, inside);
      }
      catch(...)
	{
	  succeeded = false;
	}

      if (succeeded)
	{
	  if (inside.size() > 0)
	    {
	      if (ki == 1) 
		upar = 0.5*(inside[0].first + inside[0].second);
	      else
		vpar = 0.5*(inside[0].first + inside[1].second);
	      return;
	    }
	}
    }
  return;
}

//===========================================================================
RectDomain CurveBoundedDomain::containingDomain() const
//===========================================================================
{
    RectDomain dom;
    if ((*loops_[0])[0]->dimension() == 2) {
	// We've got a parametric bound
	BoundingBox b = (*loops_[0])[0]->boundingBox();
	int n1 = (int)loops_.size();
	int i, j;
	for (j=0; j<n1; j++)     {
	    int n2 = loops_[j]->size();
	    for (i = (j==0)?1:0; i < n2; ++i) {
		b.addUnionWith((*loops_[j])[i]->boundingBox());
	    }
	}
	dom = RectDomain(Array<double, 2>(b.low()[0], b.low()[1]),
			   Array<double, 2>(b.high()[0], b.high()[1]));
    } else {
	// We should have a bound consisting of CurveOnSurface objects,
	// and we'll throw an exception if it's not!
	const CurveOnSurface& cv0
	    = dynamic_cast<const CurveOnSurface&>(*(*loops_[0])[0]);
	dom = cv0.containingDomain();
	int n1 = (int)loops_.size();
	int i, j;
	for (j=0; j<n1; j++)     {
	    int n2 = loops_[j]->size();
	    for (i = (j==0)?1:0; i < n2; ++i) {
		const CurveOnSurface& cv1
		    = dynamic_cast<const CurveOnSurface&>(*(*loops_[j])[i]);
		dom.addUnionWith(cv1.containingDomain());
	    }
	}
    }
    return dom;
}


//===========================================================================
void CurveBoundedDomain::clipWithDomain(int pardir, double parval, 
					  double tolerance, 
					  shared_ptr<ParamSurface> srf,
					  vector<shared_ptr<CurveOnSurface> >& trim_pieces) const
//===========================================================================
  // Fetch all intervals in one parameter direction
  // going through a specific point lying inside the 
  // bounded domain.
{
    double fuzzy_tol = 1.0e-10;

#ifdef DEBUG
    std::ofstream of("curveloop.g2");
    for (size_t k1=0; k1<loops_.size(); ++k1)
      for (int k2=0; k2<loops_[k1]->size(); ++k2)
	{
	  (*loops_[k1])[k2]->writeStandardHeader(of);
	  (*loops_[k1])[k2]->write(of);
	}
#endif

  // First find the intervals lying inside the parameter domain
  vector<pair<double, double> > insideInts;
  try {
      getInsideIntervals(pardir, parval, parval, tolerance, insideInts);
  } catch (...) {
      THROW("Method unstable, exception thrown!");
  }

  // Make CurveOnSurface pices for the intervals
  int ki;
  shared_ptr<SplineCurve> pcurve, gcurve1; 
  Point pt1(2), pt2(2);
  int idx1 = (pardir == 2) ? 0 : 1;
  int idx2 = (pardir == 2) ? 1 : 0;
  pt1[idx1] = pt2[idx1] = parval;

  vector<shared_ptr<ParamCurve> > c_crvs = 
    srf->constParamCurves(parval, pardir==1);
  for (size_t kj=0; kj<c_crvs.size(); ++kj)
    {
      gcurve1 = dynamic_pointer_cast<SplineCurve, ParamCurve>(c_crvs[kj]);
      if (gcurve1.get())
	{
	  for (ki=0; ki<int(insideInts.size()); ++ki)
	    {
	      gcurve1->basis().knotIntervalFuzzy(insideInts[ki].first, 
						 fuzzy_tol);
	      gcurve1->basis().knotIntervalFuzzy(insideInts[ki].second, 
						 fuzzy_tol);
	    }
	}
    }

  for (ki=0; ki<int(insideInts.size()); ki++)
    {
      shared_ptr<ParamCurve> gcurve2;
      double ta = insideInts[ki].first; 
      double tb = insideInts[ki].second;
      bool processed = false;
      for (size_t kj=0; kj<c_crvs.size(); ++kj)
	{
	  double tc = c_crvs[kj]->startparam();
	  double td = c_crvs[kj]->endparam();
	  if (tc < ta + tolerance && tb - tolerance < td)
	    {
	      ta = std::max(ta, tc);
	      tb = std::min(tb, td);
	      gcurve2 = 
		shared_ptr<ParamCurve>(c_crvs[kj]->subCurve(ta, tb));
	      processed = true;
	    }
	  else if (tc < ta + tolerance && ta - tolerance < td)
	    {
	      ta = std::max(ta, tc);
	      tb = td;
	      gcurve2 = 
		shared_ptr<ParamCurve>(c_crvs[kj]->subCurve(ta, tb));
	      if (kj < c_crvs.size()-1)
		insideInts[ki].first = c_crvs[kj+1]->startparam();
	      else
		processed = true;
	    }
	  else if (tc < tb + tolerance && tb - tolerance < td)
	    {
	      ta = tc;
	      tb = std::min(td, tb);
	      gcurve2 = 
		shared_ptr<ParamCurve>(c_crvs[kj]->subCurve(ta, tb));
	      processed = true;
	    }
	  if (processed)
	    break;
	}

      if (gcurve2.get())
	{
	  pt1[idx2] = ta;
	  pt2[idx2] = tb;
	  pcurve = shared_ptr<SplineCurve>
	    (new SplineCurve(pt1, gcurve2->startparam(), pt2, gcurve2->endparam()));

	  trim_pieces.push_back
	    (shared_ptr<CurveOnSurface>(new CurveOnSurface(srf, pcurve, 
							   gcurve2, true)));
	}
    }
}

//===========================================================================
void CurveBoundedDomain::getInsideIntervals(int pardir, double parval1, 
					    double parval2, double tolerance,
					      vector<pair<double, double> >&
					      insideInts) const
  // Fetch all intervals in one parameter direction
  // going through a specific point lying inside the 
  // bounded domain.
//===========================================================================
{

  // Get the rectangular domain containing this domain
  RectDomain parbox = containingDomain();
  double mult_fac = 100.0;  // The curve with which to intersect
  // should be much larger than the trimmed domain

  // Make a constant curve in the parameter domain, in the given
  // direct slightly larger than the found containing domain.
  double par1[2], par2[2], vec[2];
  double parint;
  int par_idx = 0;
  if (pardir == 2)
    {
      // Make constant parameter curve in 2. parameter direction.
      par1[0] = par2[0] = parval1;
      par1[1] = parbox.vmin();
      par2[1] = parbox.vmax();
      parint = std::max(par2[1] - par1[1], 0.1);
      par1[1] -= 0.1*mult_fac*parint;
      par2[1] += 0.1*mult_fac*parint;

      par_idx = 1;
    }
  else if (pardir == 1)
    {
      // Make constant parameter curve in 1. parameter direction.
      par1[1] = par2[1] = parval2;
      par1[0] = parbox.umin();
      par2[0] = parbox.umax();
      parint = std::max(par2[0] - par1[0], 0.1);
      par1[0] -= 0.1*mult_fac*parint;
      par2[0] += 0.1*mult_fac*parint;

      par_idx = 0;
    }
  else 
    {
      // Diagonal parameter curve
      par1[0] = par2[0] = parval1;
      par1[1] = par2[1] = parval2;
      vec[0] = parbox.umax() - parbox.umin();
      vec[1] = parbox.vmax() - parbox.vmin();
      par1[0] -= (mult_fac*vec[0]);
      par1[1] -= (mult_fac*vec[1]);
      par2[0] += (mult_fac*vec[0]);
      par2[1] += (mult_fac*vec[1]);

      par_idx = 0;  // @@@ VSK, 0309. This is not really correct
    }

  Point pnt1(par1[0], par1[1]);
  Point pnt2(par2[0], par2[1]);
  SplineCurve isopar(pnt1, pnt2);

  vector<intersection_point> intpt;

  findPcurveInsideSegments(isopar, tolerance, intpt);

  int ki;
  int nmbpoint = (int)intpt.size();
  if (nmbpoint == 1)
    {
      // Assuming non-recognized tangential intersection
      return;
    }
  else if (nmbpoint%2 != 0) {
    throw std::logic_error("Odd number of intersections."); 
    }

  for (ki=1; ki<nmbpoint; ki+=2)
    {
      Point eval1, eval2;
      isopar.point(eval1, intpt[ki-1].par1);
      isopar.point(eval2, intpt[ki].par1);
      insideInts.push_back(std::make_pair(eval1[par_idx], eval2[par_idx]));
      //      std::cout << eval1[2-pardir] <<" , " << eval2[2-pardir] << std::endl;
    }

  
}

//===========================================================================
int CurveBoundedDomain::positionPointInDomain(int pardir, double parval1, 
					      double parval2, 
					      double tolerance) const
  // Fetch all intervals in one parameter direction
  // going through a specific point lying inside the 
  // bounded domain.
//
// Return value: -1 : Outside of outer loop
//                0 : Inside domain
//                j>0 : Inside hole number j, i.e. inside loop number j
//===========================================================================
{

  // Get the rectangular domain containing this domain
  RectDomain parbox = containingDomain();

  // Make a constant curve in the parameter domain, in the given
  // direct slightly larger than the found containing domain.
  double par1[2], par2[2], vec[2];
  double parint;
  int par_idx = 0;
  if (pardir == 2)
    {
      // Make constant parameter curve in 2. parameter direction.
      par1[0] = par2[0] = parval1;
      par1[1] = parbox.vmin();
      par2[1] = parbox.vmax();
      parint = std::max(par2[1] - par1[1], 0.1);
      par1[1] -= 0.1*parint;
      par2[1] += 0.1*parint;

      par_idx = 1;
    }
  else if (pardir == 1)
    {
      // Make constant parameter curve in 1. parameter direction.
      par1[1] = par2[1] = parval2;
      par1[0] = parbox.umin();
      par2[0] = parbox.umax();
      parint = std::max(par2[0] - par1[0], 0.1);
      par1[0] -= 0.1*parint;
      par2[0] += 0.1*parint;

      par_idx = 0;
    }
  else 
    {
      // Diagonal parameter curve
      par1[0] = par2[0] = parval1;
      par1[1] = par2[1] = parval2;
      vec[0] = parbox.umax() - parbox.umin();
      vec[1] = parbox.vmax() - parbox.vmin();
      par1[0] -= vec[0];
      par1[1] -= vec[1];
      par2[0] += vec[0];
      par2[1] += vec[1];

      par_idx = 0;  // @@@ VSK, 0309. This is not really correct
    }

  Point pnt1(par1[0], par1[1]);
  Point pnt2(par2[0], par2[1]);
  SplineCurve isopar(pnt1, pnt2);

  vector<intersection_point> intpt;

  findPcurveInsideSegments(isopar, tolerance, intpt);

  int ki;
  int nmbpoint = (int)intpt.size();
  if (nmbpoint == 1)
    {
      // Assuming non-recognized tangential intersection
      return -1;
    }
  else if (nmbpoint%2 != 0) {
    throw std::logic_error("Odd number of intersections."); 
    }

  double ptol = 1.0e-8;
  for (ki=1; ki<nmbpoint; ki+=1)
    {
      // Check if the initial point lies inside the current interval
      Point eval1, eval2;
      isopar.point(eval1, intpt[ki-1].par1);
      isopar.point(eval2, intpt[ki].par1);
      if (eval1[0]-ptol < parval1 && eval2[0]+ptol > parval1 &&
	  eval1[1]-ptol < parval2 && eval2[1]+ptol > parval2)
	{
	  // The point lies inside this interval. Check if the interval
	  // corresponds to an inner loop
	  int loop1 = intpt[ki-1].loop_idx;
	  int loop2 = intpt[ki].loop_idx;
	  if (loop1 != loop2)
	    return 0;  // In the domain
	  else
	    return loop1;
	}
    }

  return -1; // The point lies outside the domain
}

//===========================================================================
void CurveBoundedDomain::
findPcurveInsideSegments(const SplineCurve& curve,
			 double tolerance,
			 vector<double>& params_start_end_interval) const
//===========================================================================
{
    params_start_end_interval.clear();
    vector<intersection_point> intpt;
    findPcurveInsideSegments(curve, tolerance, intpt);
    for (int i = 0; i < int(intpt.size()); ++i) {
	params_start_end_interval.push_back(intpt[i].par1);
    }
}

//===========================================================================
void CurveBoundedDomain::
findPcurveInsideSegments(const SplineCurve& curve,
			 double tolerance,
			 vector<double>& params_start_end_interval,
			 vector<double>& boundary_params,
			 vector<int>& boundary_loops,
			 vector<int>& boundary_curves) const
//===========================================================================
{
  params_start_end_interval.clear();
  boundary_params.clear();
  boundary_loops.clear();
  boundary_curves.clear();
  vector<intersection_point> intpt;
  findPcurveInsideSegments(curve, tolerance, intpt);
  for (int i = 0; i < int(intpt.size()); ++i)
    {
      params_start_end_interval.push_back(intpt[i].par1);
      boundary_params.push_back(intpt[i].par2);
      boundary_loops.push_back(intpt[i].loop_idx);
      boundary_curves.push_back(intpt[i].curve_idx);
    }
}

//===========================================================================
void CurveBoundedDomain::
findPcurveInsideSegments(const SplineCurve& curve, 
			 double tolerance, 
			 vector<intersection_point>& intpt) const 
//===========================================================================
{
    // Find all intersections between this spline curve and the 
    // parameter loops surrounding thi s domain.  Do also collect
    // pretopology information.
    
#ifdef DEBUG
  std::ofstream of("domain2D.g2");
  curve.writeStandardHeader(of);
  curve.write(of);
#endif
    int ki, kj;
    vector<pair<double,double> > intersection_par;
    vector<int> pretopology;
    vector<pair<pair<double,double>, pair<double,double> > > int_crvs;
    vector<pair<int,int> > curve_pos;
    vector<pair<int,int> > int_curve_pos;

    double epsge = 0.000001;
    for (ki=0; ki<int(loops_.size()); ki++) {
	for (kj=0; kj< loops_[ki]->size(); kj++) {
	    shared_ptr<ParamCurve> par_crv = getParameterCurve(ki, kj);
#ifdef DEBUG
	    par_crv->writeStandardHeader(of);
	    par_crv->write(of);
#endif

	    // Intersect
	    intersect2Dcurves(&curve, par_crv.get(), epsge, intersection_par, 
			      pretopology, int_crvs);
	    curve_pos.resize(intersection_par.size(), pair<int,int>(ki, kj));
	    int_curve_pos.resize(int_crvs.size(), pair<int,int>(ki, kj));
	}
    }
    
    // Check identity of curves
    for (ki=0; ki<(int)int_crvs.size(); )
      {
	for (kj=ki+1; kj<(int)int_crvs.size(); ++kj)
	  {
	    if (fabs(int_crvs[ki].first.first - int_crvs[kj].second.first) <= tolerance &&
		fabs(int_crvs[ki].second.first - int_crvs[kj].first.first) <= tolerance)
	      {
		// Identical and oppositely oriented intersection curves
		// Assuming a tangential situation
		int_crvs.erase(int_crvs.begin()+kj);
		int_crvs.erase(int_crvs.begin()+ki);
		break;
	      }
	  }
	if (kj == (int)int_crvs.size())
	  ki++;
      }

    // Sort the intersections with increasing parameters of the 
    // iso-parametric curve.
    
    int nmbpoint;
    int nmbcrv = (int)int_crvs.size();
    for (ki=0; ki<int(intersection_par.size()); ki++) {
      for (kj=0; kj<4; kj++)
	if (pretopology[4*ki+kj] == pretop_ON)
	  break;

	if (kj<4)
	    continue;   // Touching intersections not counted

	if ((pretopology[4*ki] == pretopology[4*ki+1] && 
	     (pretopology[4*ki] == pretop_IN || pretopology[4*ki] == pretop_OUT)) ||
	    (pretopology[4*ki+2] == pretopology[4*ki+3] && 
	     (pretopology[4*ki+2] == pretop_IN || pretopology[4*ki+2] == pretop_OUT)))
	    continue;   // Touching intersections not counted
	
	nmbpoint = (int)intpt.size();
	for (kj=0; kj<nmbpoint; kj++)
	    if (fabs(intersection_par[ki].first - 
		     intpt[kj].par1) <= epsge /*tolerance*/)
	      break;
	if (kj < nmbpoint)
	    continue;   // Do not count corner points twice

	// Check towards endpoints of curves
	for (kj=0; kj<nmbcrv; ++kj)
	  {
	    if (fabs(intersection_par[ki].first - int_crvs[kj].first.first) <= epsge /*tolerance*/ ||
		fabs(intersection_par[ki].first - int_crvs[kj].second.first) <= epsge /*tolerance*/)
	      break;
	  }
	if (kj < nmbcrv)
	    continue;   // Do not count corner points twice
	
	intpt.push_back(intersection_point(intersection_par[ki].first,
					   intersection_par[ki].second,
					   curve_pos[ki].first,
					   curve_pos[ki].second,
					   &pretopology[4*ki]));
    }
    nmbpoint = (int)intpt.size();

    std::sort(intpt.begin(), intpt.end(), par1_compare);

    // Represent the point by its endpoints
    int dummy_pretop[4];
    dummy_pretop[0] = dummy_pretop[1] = dummy_pretop[2] = dummy_pretop[3] = 0;
    for (ki=0; ki<nmbcrv; ++ki)
      {
	double tpar = std::min(int_crvs[ki].first.first,int_crvs[ki].second.first);

	for (kj=0; kj<(int)intpt.size(); ++kj)
	  if (intpt[kj].par1 >= tpar)
	    break;

	intpt.insert(intpt.begin()+kj, 
		     intersection_point(int_crvs[ki].first.first,
					int_crvs[ki].first.second,
					int_curve_pos[ki].first,
					int_curve_pos[ki].second,
					dummy_pretop));

	int ix = (int_crvs[ki].first.first < int_crvs[ki].second.first)
	  ? 1 : 0;

	intpt.insert(intpt.begin()+kj+ix, 
		     intersection_point(int_crvs[ki].second.first,
					int_crvs[ki].second.second,
					int_curve_pos[ki].first,
					int_curve_pos[ki].second,
					dummy_pretop));
      }
}

//===========================================================================
shared_ptr<ParamCurve> CurveBoundedDomain::getParameterCurve(int loop_nmb,
							       int curve_nmb) const
//===========================================================================
{
    shared_ptr<ParamCurve> curr_crv = (*loops_[loop_nmb])[curve_nmb];
    shared_ptr<ParamCurve> par_crv;
    if (curr_crv->instanceType() == Class_CurveOnSurface) {
      shared_ptr<CurveOnSurface> sf_cv = 
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curr_crv);
	par_crv = sf_cv->parameterCurve();
	if (par_crv.get() == 0) {
	    // Project the 3D curve onto the surface
	  MESSAGE("CurveBoundedDomain::getParameterCurve, arbitrary tolerance");
	  double tol = 1.0e-4;
	  sf_cv->ensureParCrvExistence(tol);
	  par_crv = sf_cv->parameterCurve();
	  if (par_crv.get() == 0)
	    THROW("CurveBoundedDomain::getParameterCurve. Missing par curve");
	}
    } else if (curr_crv->dimension() == 2) {
	par_crv = curr_crv;
    } else
	THROW("Wrong curve in loop in bounded domain.");

    return par_crv;
}
