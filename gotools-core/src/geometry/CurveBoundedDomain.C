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
  else 
    {
      // Boundary intersections are caught. Can use a small tolerance
      double tol = std::min(tolerance, 1.0e-6);
      if (isInDomain(pnt, tol))
	return 1;
      else
	return 0;
    }
}

//===========================================================================
bool CurveBoundedDomain::isInDomain(const Array<double, 2>& pnt,
				      double tolerance) const
//===========================================================================
{
  // Boundary points are critical. Check first if the point lies at a boundary 
  if (isOnBoundary(pnt, tolerance))
    return true;

  RectDomain dom = containingDomain();
  double dist1 = std::min(pnt[0]-dom.umin(), dom.umax()-pnt[0]);
  double dist2 = std::min(pnt[1]-dom.vmin(), dom.vmax()-pnt[1]);
  if (dist1 < -tolerance || dist2 < -tolerance)
    return false;
  double frac = std::min(dist1,dist2)/std::max(dist1,dist2);
  
  // Sort intersection directions to avoid coincident intersections
  // Prefer non-constant parameter directions
  int dir[3];
  int ix = (dist1 > dist2) ? 0 : 1;
  dir[0] = (frac<0.1) ? 2-ix /*ix+1*/ : 3;
  dir[1] = (frac<0.1) ? 3 : 2-ix /*ix+1*/;
  dir[2] = ix+1 /*2-ix*/;
      
  int nmb_catches = 0;
  for (int kj=0; kj<3; ++kj)
    {
      try {
	// Insert the domain with a constant parameter curve
	vector<pair<double, double> > inside;
	getInsideIntervals(dir[kj], pnt[0], pnt[1], tolerance, inside);
	// Check if the pnt is inside any of the parameter intervals
	// lying inside the domain.
	int nmbint = (int)inside.size();
	ix = (dir[kj] == 2) ? 1 : 0; 
	int ki;
	for (ki=0; ki<nmbint; ki++) {
	  if (inside[ki].first-tolerance <= pnt[ix] &&
	      inside[ki].second+tolerance >= pnt[ix]) {
	    return true;
	  }
	}
	return false;
      } catch (...) { // It seems there were an odd number of intersections.
	// We try intersecting in an other direction.
	//       MESSAGE("Unstable method, trying in the other direction.");
	++nmb_catches;
      }
    }

  if (nmb_catches == 2) { // If last attempt was a success, we're satisfied.
    // Not boundary
    THROW("Failed to decide whether point was inside boundary.");
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
  for (int ki=1; ki<3; ++ki)
    {
      bool succeeded = true;
      vector<pair<double, double> > inside;
      try {
	getInsideIntervals(ki, upar, vpar, tolerance, inside, false);
      }
      catch(...)
	{
	  succeeded = false;
	}

      if (succeeded)
	{
	  if (inside.size() > 0)
	    {
	      // Find the largest interval
	      double max_len = inside[0].second - inside[0].first;
	      size_t max_ix = 0;
	      for (size_t kj=0; kj<inside.size(); ++kj)
		{
		  double len = inside[kj].second - inside[kj].first;
		  if (len > max_len)
		    {
		      max_len = len;
		      max_ix = kj;
		    }
		}
	      if (ki == 1) 
		upar = 0.5*(inside[max_ix].first + inside[max_ix].second);
	      else
		vpar = 0.5*(inside[max_ix].first + inside[max_ix].second);

	      // Check if we have found a boundary point
	      Vector2D ppnt(upar, vpar);
	      bool is_boundary = isOnBoundary(ppnt, tolerance);
	      if (!is_boundary)
		break;
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
      if (tb - ta < tolerance)
	continue;
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
					    insideInts, bool with_bd) const
  // Fetch all intervals in one parameter direction
  // going through a specific point lying inside the 
  // bounded domain.
//===========================================================================
{

  // Get the rectangular domain containing this domain
  RectDomain parbox = containingDomain();
  Point mid(0.5*(parbox.umin()+parbox.umax()), 
	    0.5*(parbox.vmin()+parbox.vmax()));
  Point parpt(parval1, parval2);
  double len1 = (parbox.umax()-parbox.umin())+(parbox.vmax()-parbox.vmin());
  double len2 = mid.dist(parpt);
  if (len1 == 0.0)
  {
      MESSAGE("Degenerate domain!");
      return;
  }

  double mult_fac = 2.0*(len1+len2)/len1;  // The curve with which to intersect
  // should be much larger than the trimmed domain

  // Make a constant curve in the parameter domain, in the given
  // direct slightly larger than the found containing domain.
  double par1[2], par2[2], vec[2];
  double parint;
  int par_idx = 0;
  double len;
  if (pardir == 2)
    {
      // Make constant parameter curve in 2. parameter direction.
      par1[0] = par2[0] = parval1;
      par1[1] = parbox.vmin();
      par2[1] = parbox.vmax();
      parint = std::max(par2[1] - par1[1], 0.1);
      len = mult_fac*parint;
      par1[1] -= len;
      par2[1] += len;

      par_idx = 1;
    }
  else if (pardir == 1)
    {
      // Make constant parameter curve in 1. parameter direction.
      par1[1] = par2[1] = parval2;
      par1[0] = parbox.umin();
      par2[0] = parbox.umax();
      parint = std::max(par2[0] - par1[0], 0.1);
      len = mult_fac*parint;
      par1[0] -= len;
      par2[0] += len;

      par_idx = 0;
    }
  else 
    {
      // Diagonal parameter curve
      par1[0] = par2[0] = parval1;
      par1[1] = par2[1] = parval2;
      vec[0] = parbox.umax() - parbox.umin();
      vec[1] = parbox.vmax() - parbox.vmin();
      len = mult_fac*sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
      par1[0] -= (mult_fac*vec[0]);
      par1[1] -= (mult_fac*vec[1]);
      par2[0] += (mult_fac*vec[0]);
      par2[1] += (mult_fac*vec[1]);

      par_idx = 0;  // @@@ VSK, 0309. This is not really correct
    }

  Point pnt1(par1[0], par1[1]);
  Point pnt2(par2[0], par2[1]);
  SplineCurve isopar(pnt1, -len, pnt2, len);

  vector<intersection_point> intpt;

  findPcurveInsideSegments(isopar, tolerance, intpt, with_bd);

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
			 vector<double>& params_start_end_interval,
			 bool with_bd) const
//===========================================================================
{
    params_start_end_interval.clear();
    vector<double> int_seg;
    vector<intersection_point> intpt;
    findPcurveInsideSegments(curve, tolerance, intpt, with_bd);
    for (int i = 0; i < int(intpt.size()); ++i) {
	int_seg.push_back(intpt[i].par1);
    }

    // Check for curve endpoints internal to the domain
    Point pos1 = curve.ParamCurve::point(curve.startparam());
    Point pos2 = curve.ParamCurve::point(curve.endparam());
    double dist1 = 2.0*tolerance;
    double dist2 = 2.0*tolerance;
    if (int_seg.size() > 0)
      {
	Point pos3 = curve.ParamCurve::point(int_seg[0]);
	Point pos4 = curve.ParamCurve::point(int_seg[int_seg.size()-1]);
	dist1 = pos1.dist(pos3);
	dist2 = pos2.dist(pos4);
      }

    bool in1 = true, in2 = true;
    if (dist1 > tolerance)
      {
	Vector2D pos1_2(pos1[0], pos1[1]);
	if (isInDomain2(pos1_2, tolerance))
	  int_seg.insert(int_seg.begin(), curve.startparam());
	else
	  in1 = false;
      }

    if (dist2 > tolerance)
      {
	Vector2D pos2_2(pos2[0], pos2[1]);
	if (isInDomain2(pos2_2, tolerance))
	  int_seg.push_back(curve.endparam());
	else
	  in2 = false;
      }

    if (dist1 < tolerance && dist2 < tolerance && int_seg.size() == 2)
      {
	// Check if the entire curve lies outside
	Point mid = curve.ParamCurve::point(0.5*(curve.startparam() + curve.endparam()));
	Vector2D mid_2(mid[0], mid[1]);
	if (!isInDomain2(mid_2, tolerance))
	  int_seg.clear();
      }
	

    if (int_seg.size() % 2 == 1)
      {
	if (!in1)
	  int_seg.erase(int_seg.begin());
	else if (!in2)
	  int_seg.pop_back();
      }

    if (int_seg.size() == 0)
      return;  // No internal segments

    // Check results. First mark segments
    double knot_diff_tol = 1.0e-8; // A small tolerance
    vector<int> int_seg_type(int_seg.size()-1, -1);
    // 0 = short curve, 1 = outside, 2 = inside, 3 = on the boundary
    for (int j = 0; j < int(int_seg.size()) - 1; ++j) {
      double from_par = int_seg[j];
      double to_par = int_seg[j+1];
      Point tmp1 = curve.ParamCurve::point(from_par);
      Point tmp2 = curve.ParamCurve::point(0.5*(from_par+to_par));
      Point tmp3 = curve.ParamCurve::point(to_par);
      double len = tmp1.dist(tmp2) + tmp2.dist(tmp3);
      if (to_par - from_par < knot_diff_tol || len < tolerance) {
	int_seg_type[j] = 0;  // Short curve
	continue;
      }
      if (from_par < curve.startparam())
	from_par = curve.startparam();
      if (to_par > curve.endparam())
	to_par = curve.endparam();
      double med_par = 0.5*(from_par + to_par);
      Point med_pt = curve.ParamCurve::point(med_par);
      int is_in_domain = -1;
      try {
	is_in_domain = isInDomain2(Vector2D(med_pt[0], med_pt[1]), tolerance);
	if (is_in_domain == 2)
	  {
	    int is_in_domain2 = 
	      isInDomain2(Vector2D(med_pt[0], med_pt[1]), knot_diff_tol);
	    if (is_in_domain2 == 1)
	      is_in_domain = 1;
	  }
      }
      catch (...)
	{
	  is_in_domain = -2;
	}
      int_seg_type[j] = is_in_domain + 1;
    }

    // Simplify segmentation by joining segments of the same type (inside/outside)
    for (int j=0; j<int(int_seg_type.size()); ++j)
      {
	int k;
	for (k=j+1; k<int(int_seg_type.size()); ++k)
	  if (int_seg_type[k] != int_seg_type[j])
	    break;
	if (k > j+1)
	  {
	    // Simplify
	    int type = int_seg_type[j];
	    if (int_seg_type[j] == 0)
	      {
		// Check if the segment is still small
		double from_par = int_seg[j];
		double to_par = int_seg[k];
		Point tmp1 = curve.ParamCurve::point(from_par);
		Point tmp2 = curve.ParamCurve::point(0.5*(from_par+to_par));
		Point tmp3 = curve.ParamCurve::point(to_par);
		double len = tmp1.dist(tmp2) + tmp2.dist(tmp3);
		if (to_par - from_par > knot_diff_tol && len > tolerance) 
		  {
		    Point med_pt = curve.ParamCurve::point(0.5*(from_par+to_par));
		    int is_in_domain = 0;
		    try {
		      is_in_domain = 
			isInDomain2(Vector2D(med_pt[0], med_pt[1]), 
					   knot_diff_tol);
		    }
		    catch (...)
		      {
			if (len > 10.0*tolerance)
			  THROW("Could not trim intersection curve with surface boundary");
			else
			  is_in_domain = 0;
		      }
		    type = (is_in_domain >= 1) ? 2 : 1;
		  }
		else
		  {
		    int_seg.erase(int_seg.begin()+j+1, int_seg.begin()+k);
		    int_seg_type.erase(int_seg_type.begin()+j+1, int_seg_type.begin()+k);
		  }
	      }
	    else
	      {
		int_seg.erase(int_seg.begin()+j+1, int_seg.begin()+k);
		int_seg_type.erase(int_seg_type.begin()+j+1, int_seg_type.begin()+k);
	      }
	    int_seg_type[j] = type;
	  }
      }

    for (int j=1; j<int(int_seg_type.size()-1); ++j)
      {
	if (int_seg_type[j] == 0)
	  {
	    if (int_seg_type[j-1] == int_seg_type[j+1])
	      {
		// Remove small segment
		int_seg.erase(int_seg.begin()+j, int_seg.begin()+j+1);
		int_seg_type.erase(int_seg_type.begin()+j, int_seg_type.begin()+j+1);
		--j;
	      }
	  }
      }

    // Translate segments to output array
    for (int j=1; j<int(int_seg.size()); ++j)
      {
	if (int_seg_type[j-1] == 2 || int_seg_type[j-1] == 3)
	  {
	    if (params_start_end_interval.size() == 0 || 
		int_seg[j-1]>params_start_end_interval[params_start_end_interval.size()-1])
	      params_start_end_interval.push_back(int_seg[j-1]);
	    params_start_end_interval.push_back(int_seg[j]);
	  }
      }
 }

//===========================================================================
void CurveBoundedDomain::
findPcurveInsideSegments(const SplineCurve& curve,
			 double tolerance,
			 vector<double>& params_start_end_interval,
			 vector<double>& boundary_params,
			 vector<int>& boundary_loops,
			 vector<int>& boundary_curves,
			 bool with_bd) const
//===========================================================================
{
  params_start_end_interval.clear();
  boundary_params.clear();
  boundary_loops.clear();
  boundary_curves.clear();
  vector<intersection_point> intpt;
  findPcurveInsideSegments(curve, tolerance, intpt, with_bd);
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
			 vector<intersection_point>& intpt,
			 bool with_bd) const 
//===========================================================================
{
    // Find all intersections between this spline curve and the 
    // parameter loops surrounding this domain.  Do also collect
    // pretopology information.
    
#ifdef DEBUG
  std::ofstream of("domain2D.g2");
  curve.writeStandardHeader(of);
  curve.write(of);
#endif
    int ki, kj;
    vector<pair<double,double> > intersection_par;
    vector<pair<int,int> > intersection_ix;
    vector<int> pretopology;
    vector<pair<pair<double,double>, pair<double,double> > > int_crvs;
    vector<pair<int,int> > curve_pos;
    vector<pair<int,int> > int_curve_pos;

    vector<vector<bool> > segment_deg(loops_.size());

    const double deg_tol = 1.0e-12;
    const double knot_diff_tol = 1.0e-8; // A small tolerance

    double epsge = tolerance; //0.000001;   // This is a potential unstability

    // Compute maximum distance between adjacent curves
    double max_dist = 0.0;
    double max_tol = 0.1;   // To avoid meaningless intersections for
    // loops with large gaps
    
    for (ki=0; ki<int(loops_.size()); ki++) 
      {
	int nmb_cvs = loops_[ki]->size();
	shared_ptr<ParamCurve> par_crv1 = getParameterCurve(ki, nmb_cvs-1);
	double par1 = par_crv1->endparam();
	Point pos1 = par_crv1->point(par1);
	for (kj=0; kj < nmb_cvs; kj++) 
	  {
	    shared_ptr<ParamCurve> par_crv2 = getParameterCurve(ki, kj);
	    double par2 = par_crv2->startparam();
	    Point pos2 = par_crv2->point(par2);
	    double dist = pos1.dist(pos2);
	    max_dist = std::max(max_dist, dist);

	    if (dist > epsge)
	      {
		// There is a risk of loosing intersections at curve joints.
		// // Make a pre check for endpoint intersections
		// vector<double> int1, int2;
		// vector<pair<double, double> > intcv1, intcv2;
		// intersectCurvePoint(&curve, pos1, std::min(dist,max_tol), 
		// 		    int1, intcv1);
		// intersectCurvePoint(&curve, pos2, std::min(dist,max_tol), 
		// 		    int2, intcv2);
		shared_ptr<ParamCurve> gap_cv(new SplineCurve(pos1, pos2));
		int curr_nmb_par = (int)intersection_par.size();
                // The pretopology object consists of 4 ints:
                // cv1 from right, cv1 from left, cv2 from right, cv2 from left.
                // The values are: enum: SI_UNDEF, SI_IN, SI_OUT, SI_ON, SI_AT.
                // SI_IN: The segment is to the right of the other curve.
                // SI_OUT: The segment is to the left of the other curve.
                // SI_ON: Intersecting in an interval.
                // SI_AT: Intersecting at the end point.
		intersect2Dcurves(&curve, gap_cv.get(), epsge, intersection_par,
				  pretopology, int_crvs);
		for (int kr=curr_nmb_par; kr<(int)intersection_par.size(); ++kr)
                {
                    intersection_ix.push_back(std::make_pair(ki, -1)); // No loop curve
                }
#ifdef DEBUG
		gap_cv->writeStandardHeader(of);
		gap_cv->write(of);
#endif
		
 	      }
	    par1 = par_crv2->endparam();
	    pos1 = par_crv2->point(par1);
	  }
      }

    double eps2 = std::min(max_dist,max_tol) + epsge;
    for (ki=0; ki<int(loops_.size()); ki++) {
	for (kj=0; kj< loops_[ki]->size(); kj++) {
	    shared_ptr<ParamCurve> par_crv = getParameterCurve(ki, kj);
#ifdef DEBUG
	    par_crv->writeStandardHeader(of);
	    par_crv->write(of);
#endif

            // Degenerate curves lack the topology needed by this routine. Any intersections will be
            // handled by the adjacent segments. And degenerate parameter curves are not needed and
            // should have been removed from the loop.
            if (par_crv->estimatedCurveLength() < deg_tol)
            {
                segment_deg[ki].push_back(true);
                continue;
            }
            else
            {
                segment_deg[ki].push_back(false);
            }
	    // Intersect
	    int curr_nmb_par = (int)intersection_par.size();
	    intersect2Dcurves(&curve, par_crv.get(), epsge, intersection_par, 
			      pretopology, int_crvs);
		curve_pos.resize(intersection_par.size(), pair<int,int>(ki, kj));
	    if (with_bd)
	      int_curve_pos.resize(int_crvs.size(), pair<int,int>(ki, kj));
	    for (int kr=curr_nmb_par; kr<(int)intersection_par.size(); ++kr)
	      intersection_ix.push_back(std::make_pair(ki,kj));
	}
    }
    
    if (!with_bd)
      {
	int_crvs.clear();
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

    // Remove double set of intersections at curveloop corners
    for (ki=0; ki<(int)intersection_par.size(); ++ki)
    {
	Point pos1 = curve.ParamCurve::point(intersection_par[ki].first);
	for (kj=ki+1; kj<(int)intersection_par.size(); ++kj)
        {
	    Point pos2 = curve.ParamCurve::point(intersection_par[kj].first);
            double dist = pos1.dist(pos2);
	    if ((dist < eps2) && (intersection_ix[ki].first == intersection_ix[kj].first))
            {
		// Potential corner, check
		int ix_min = std::min(intersection_ix[ki].second, 
				      intersection_ix[kj].second);
		int ix_max = std::max(intersection_ix[ki].second, 
				      intersection_ix[kj].second);

                // We count the number of degenerate segments between the two segments.
                int loop_ind = intersection_ix[ki].first;
                int nmb_deg = 0;
                for (int kk = ix_min + 1; kk < ix_max; ++kk)
                {
                    if (segment_deg[loop_ind][kk])
                    {
                        ++nmb_deg;
                    }
                }

                // We check the topology. If either both are SI_IN or both are SI_OUT the intersection is
                // corner tangential and can be kept as a double intersection, not influencing the
                // outcome of the odd/even number of intersections test.  I.e., for each 4-tuple
                // pretopology object we check the last two values. If the opposite pairs are both 2 or
                // both 1 the intersection is corner tangential.

                // If we are touching a corner (i.e. not crossing) we must remove both or none.
                int pre1_1 = pretopology[ki*4+2];
                int pre1_2 = pretopology[ki*4+3];
                int pre2_1 = pretopology[kj*4+2];
                int pre2_2 = pretopology[kj*4+3];
                bool both_at = (((pre1_1 == pretop_AT) || (pre1_2 == pretop_AT)) &&
                                ((pre2_1 == pretop_AT) || (pre2_2 == pretop_AT)));
                bool both_in = (((pre1_1 == pretop_IN) && (pre2_2 == pretop_IN)) ||
                                ((pre1_2 == pretop_IN) && (pre2_1 == pretop_IN)));
                bool both_out = (((pre1_1 == pretop_OUT) && (pre2_2 == pretop_OUT)) ||
                                 ((pre1_2 == pretop_OUT) && (pre2_1 == pretop_OUT)));
                bool remove_both = (both_at && (both_in || both_out));
		int remove_ix = -1;
		if (ix_min < 0 && ix_max < 0)
                    remove_ix = kj;  // Remove one
		else if (ix_min < 0)
                    remove_ix = (intersection_ix[ki].second < 0) ? ki : kj;
		
		if (remove_ix >= 0)
                { // This refers to a gap curve added to the loop.
		    // Superflous intersection found at loop gap
		    intersection_par.erase(intersection_par.begin()+remove_ix);
		    intersection_ix.erase(intersection_ix.begin()+remove_ix);
                    pretopology.erase(pretopology.begin()+remove_ix*4, pretopology.begin()+(remove_ix+1)*4);
		    --kj;
                    // If both segments lie on the same side of the input curve they should both be removed.
                    if (remove_both)
                    {
                        int other_ind = (remove_ix == ki) ? kj : ki;
                        intersection_par.erase(intersection_par.begin()+other_ind);
                        intersection_ix.erase(intersection_ix.begin()+other_ind);
                        pretopology.erase(pretopology.begin()+other_ind*4, pretopology.begin()+(other_ind+1)*4);
                        --kj;
                    }

		    if ((remove_ix == ki) || remove_both)
                    {
			--ki;
			break;
                    }
                }
                else
                {

                    if (ix_max - ix_min - nmb_deg == 1 ||
                        (ix_max == loops_[intersection_ix[ki].first]->size()-1 &&
                         ix_min == 0))
                    {
                        // Adjacent curves, check position of intersection curve
                        shared_ptr<ParamCurve> cv1 = 
			    (*loops_[intersection_ix[ki].first])[intersection_ix[ki].second];
                        shared_ptr<ParamCurve> cv2 = 
			    (*loops_[intersection_ix[kj].first])[intersection_ix[kj].second];
                        double tpar1 = intersection_par[ki].second;
                        double tpar2 = intersection_par[kj].second;
                        Point pt1, pt2, pt3, pt4;
                        pt1 = (tpar1-cv1->startparam() < cv1->endparam()-tpar1) ?
			    cv1->point(cv1->startparam()) :
			    cv1->point(cv1->endparam());
                        pt2 = (tpar2-cv2->startparam() < cv2->endparam()-tpar2) ?
			    cv2->point(cv2->startparam()) :
			    cv2->point(cv2->endparam());
                        pt3 = cv1->point(tpar1);
                        pt4 = cv2->point(tpar2);
                        // Use geometry space test to stay consistent with
                        // the computation of the intersection points
                        // if ((tpar1-cv1->startparam() < epsge ||
                        //      cv1->endparam()-tpar1 < epsge) &&
                        //     (tpar2-cv2->startparam() < epsge ||
                        //      cv2->endparam()-tpar2 < epsge))
                        if (pt1.dist(pt3) < eps2 && pt2.dist(pt4) < eps2)
                        {
                            // Corner
                            intersection_par.erase(intersection_par.begin()+kj);
                            intersection_ix.erase(intersection_ix.begin()+kj);
                            pretopology.erase(pretopology.begin()+kj*4, pretopology.begin()+(kj+1)*4);
                            --kj;
                        }
                    }
                }
            }
        }
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
	// for (kj=0; kj<nmbpoint; kj++)
	//     if (fabs(intersection_par[ki].first - 
	// 	     intpt[kj].par1) <= epsge /*tolerance*/)
	//       break;
	// if (kj < nmbpoint)
	//     continue;   // Do not count corner points twice

	// Check towards endpoints of curves
	for (kj=0; kj<nmbcrv; ++kj)
	  {
	    if (fabs(intersection_par[ki].first - int_crvs[kj].first.first) <= epsge /*tolerance*/ ||
		fabs(intersection_par[ki].first - int_crvs[kj].second.first) <= epsge /*tolerance*/)
	      {
		break;
	      }
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

    // Represent the curve depending on the configuration
    int dummy_pretop[4];
    dummy_pretop[0] = dummy_pretop[1] = dummy_pretop[2] = dummy_pretop[3] = 0;
    for (ki=0; ki<nmbcrv; ++ki)
      {
	if (int_crvs[ki].second.first < int_crvs[ki].first.first)
	  std::swap(int_crvs[ki].first, int_crvs[ki].second);
	double tpar1 = int_crvs[ki].first.first;
	double tpar2 = int_crvs[ki].second.first;

	for (kj=0; kj<(int)intpt.size(); ++kj)
	  if (intpt[kj].par1 >= tpar1)
	    break;

	int ix = 0;
	if (kj%2 == 0)
	  {
	    intpt.insert(intpt.begin()+kj, 
			 intersection_point(int_crvs[ki].first.first,
					    int_crvs[ki].first.second,
					    int_curve_pos[ki].first,
					    int_curve_pos[ki].second,
					    dummy_pretop));
	    ++ix;
	  }

	int nb;
	for (nb=(int)intpt.size(); nb>0 && intpt[nb-1].par1>=tpar2; nb--);

	if ((nb-(int)intpt.size()%2) == 0)
	  intpt.insert(intpt.begin()+kj+ix, 
		       intersection_point(int_crvs[ki].second.first,
					  int_crvs[ki].second.second,
					  int_curve_pos[ki].first,
					  int_curve_pos[ki].second,
					  dummy_pretop));
      }
}

//===========================================================================
bool CurveBoundedDomain::doIntersect(const SplineCurve& curve, 
				     double tol) const
//===========================================================================
{
  if (curve.dimension() != 2)
    THROW("Dimension of parameter curve different from 2");

    for (size_t ki=0; ki<int(loops_.size()); ki++) {
      for (size_t kj=0; kj< loops_[ki]->size(); kj++) {
	vector<pair<double,double> > intersection_par;
	vector<int> pretopology;
	vector<pair<pair<double,double>, pair<double,double> > > int_crvs;
	
	shared_ptr<ParamCurve> par_crv = getParameterCurve(ki, kj);

	intersect2Dcurves(&curve, par_crv.get(), tol, intersection_par, 
			  pretopology, int_crvs);

	if (intersection_par.size() > 0 || int_crvs.size() > 0)
	  return true;
      }
    }

    return false;
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
