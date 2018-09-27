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

#include "GoTools/intersections/Identity.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/intersections/GeoTol.h"

using std::vector;

namespace Go
{
    int Identity::identicalSfs(shared_ptr<ParamSurface> sf1,
			       shared_ptr<ParamSurface> sf2,
			       double tol)
    {
	// Initialize intersection objects
	shared_ptr<ParamSurfaceInt> intsf1 = 
	    shared_ptr<ParamSurfaceInt>(new ParamSurfaceInt(sf1));
	shared_ptr<ParamSurfaceInt> intsf2 = 
	    shared_ptr<ParamSurfaceInt>(new ParamSurfaceInt(sf2));
	shared_ptr<GeoTol> eps = shared_ptr<GeoTol>(new GeoTol(tol));

	// Fetch boundary elements
	vector<shared_ptr<BoundaryGeomInt> > bd_cvs1;
	vector<shared_ptr<BoundaryGeomInt> > bd_cvs2;

	intsf1->getBoundaryObjects(bd_cvs1);
	intsf2->getBoundaryObjects(bd_cvs2);

	// Check identity of boundary curves
	int ki, kj;
	int coincident = 0;
	double start1, start2, end1, end2;
	for (ki=0; ki<int(bd_cvs1.size()); ++ki)
	{
	    ParamCurveInt *cv1 = bd_cvs1[ki]->getObject()->getParamCurveInt();
	    start1 = cv1->startParam(0);
	    end1 = cv1->endParam(0);
	    for (kj=0; kj<int(bd_cvs2.size()); ++kj)
	    {
		ParamCurveInt *cv2 = bd_cvs2[kj]->getObject()->getParamCurveInt();
		start2 = cv2->startParam(0);
		end2 = cv2->endParam(0);

		// Check orientation
		Point pt1, pt2, pt3, pt4;
		cv1->point(pt1, &start1);
		cv1->point(pt2, &end1);
		cv2->point(pt3, &start2);
		cv2->point(pt4, &end2);
		if (!(pt1.dist(pt3) < tol || pt1.dist(pt4) < tol))
		    continue;
		if (!(pt2.dist(pt3) < tol || pt2.dist(pt4) < tol))
		    continue;
		if (pt1.dist(pt3) < pt1.dist(pt4))
		    coincident = checkCoincide(cv1, start1, end1, eps,
					       cv2, start2, end2);
		else
		    coincident = checkCoincide(cv1, start1, end1, eps,
					       cv2, end2, start2);
		if (coincident)
		    break;
	    }
	    if (kj == int(bd_cvs2.size()))
		break;  // No coincidence for this boundary curve
	}

	if (ki == int(bd_cvs1.size()))
	{
	    // Coincidence of boundary curves found. Check the inner
	    coincident = internalCoincidence(intsf1, intsf2, eps);
	    if (coincident)
		return 1;
	}

	// Check if the boundary curves of surface 1 lies in surface 2
	for (ki=0; ki<int(bd_cvs1.size()); ++ki)
	{
	    ParamCurveInt *cv1 = bd_cvs1[ki]->getObject()->getParamCurveInt();
	    start1 = cv1->startParam(0);
	    end1 = cv1->endParam(0);

	    // Project the endpoints onto surface 2
	    Point pt1, pt2, clo_pt1, clo_pt2;
	    double u1, v1, u2, v2, dist1, dist2;
	    cv1->point(pt1, &start1);
	    cv1->point(pt2, &end1);
	    sf2->closestPoint(pt1, u1, v1, clo_pt1, dist1, tol);
	    if (dist1 > tol)
		break;  // No coincidence

	    sf2->closestPoint(pt2, u2, v2, clo_pt2, dist2, tol);
	    if (dist2 > tol)
		break;  // No coincidence

	    // Check curve
	    coincident = checkCoincide(cv1, start1, end1,
				       intsf2.get(), Point(u1,v1), Point(u2,v2),
				       eps);
	    if (!coincident)
		break;
	}
	if (ki == int(bd_cvs1.size()))
	{
	    // Coincidence of boundary curves found. Check the inner
	    coincident = internalCoincidence(intsf1, intsf2, eps);
	    // if (coincident)
	    // 	return 2;
	}


	// Check if the boundary curves of surface 2 lies in surface 1
	int coincident2 = 0;
	for (ki=0; ki<int(bd_cvs2.size()); ++ki)
	{
	    ParamCurveInt *cv2 = bd_cvs2[ki]->getObject()->getParamCurveInt();
	    start2 = cv2->startParam(0);
	    end2 = cv2->endParam(0);

	    // Project the endpoints onto surface 2
	    Point pt1, pt2, clo_pt1, clo_pt2;
	    double u1, v1, u2, v2, dist1, dist2;
	    cv2->point(pt1, &start2);
	    cv2->point(pt2, &end2);
	    sf1->closestPoint(pt1, u1, v1, clo_pt1, dist1, tol);
	    if (dist1 > tol)
		break;  // No coincidence

	    sf1->closestPoint(pt2, u2, v2, clo_pt2, dist2, tol);
	    if (dist1 > tol)
		break;  // No coincidence

	    // Check curve
	    coincident = checkCoincide(cv2, start2, end2,
				       intsf1.get(), Point(u1,v1), Point(u2,v2),
				       eps);
	    if (!coincident)
		break;
	}
	if (ki == int(bd_cvs2.size()))
	{
	    // Coincidence of boundary curves found. Check the inner
	    coincident2 = internalCoincidence(intsf2, intsf1, eps);
	    if (coincident && coincident2)
	      return 1;
	    else if (coincident2)
	      return 3;
	}

	// The surfaces are neither identical nor is one embedded in the other
	return 0;
    }


    int Identity::identicalCvs(shared_ptr<ParamCurve> cv1, 
			     shared_ptr<ParamCurve> cv2, 
			     double tol)
    // Check if two curves are identical, or one is embedded in the other.
    {
      return identicalCvs(cv1, cv1->startparam(), cv1->endparam(),
			  cv2, cv2->startparam(), cv2->endparam(), tol);
    }


    int Identity::identicalCvs(shared_ptr<ParamCurve> cv1, double start1, double end1,
			       shared_ptr<ParamCurve> cv2, double start2, double end2,
			       double tol)
    // Check if two curves are identical, or one is embedded in the other.
    // The curve extension is limited by start and end parameters of each curve
    {
	// Box test
	BoundingBox box1 = cv1->boundingBox();
	BoundingBox box2 = cv2->boundingBox();
	if (!box1.overlaps(box2, tol))
	    return 0;
	
	// Initialize intersection objects
	shared_ptr<ParamCurveInt> intcv1 = 
	    shared_ptr<ParamCurveInt>(new ParamCurveInt(cv1));
	shared_ptr<ParamCurveInt> intcv2 = 
	    shared_ptr<ParamCurveInt>(new ParamCurveInt(cv2));
	shared_ptr<GeoTol> eps = shared_ptr<GeoTol>(new GeoTol(tol));

	// Check endpoints
	Point pt1, pt2, pt3, pt4;
	intcv1->point(pt1, &start1);
	intcv1->point(pt2, &end1);
	intcv2->point(pt3, &start2);
	intcv2->point(pt4, &end2);
	
	// First check coincidence
	int coincident = 0;
	if (pt1.dist(pt3) <= tol && pt2.dist(pt4) <= tol)
	    coincident = checkCoincide(intcv1.get(), start1, end1, eps,
				       intcv2.get(), start2, end2);
	else if (pt1.dist(pt4) <= tol && pt2.dist(pt3) <= tol)
	    coincident = checkCoincide(intcv1.get(), start1, end1, eps,
				       intcv2.get(), end2, start2);
	else
	{
	    // Project the endpoints on one curve onto the other curve and
	    // check for embedded curves
	    // First check if the first curve is embedded into the second
	    Point clo_pt1, clo_pt2;
	    double clo_dist1, clo_dist2;
	    double clo_par1, clo_par2;

	    cv2->closestPoint(pt1, start2, end2, clo_par1, clo_pt1, clo_dist1);
	    cv2->closestPoint(pt2, start2, end2, clo_par2, clo_pt2, clo_dist2);
	    if (clo_dist1 <= tol && clo_dist2 <= tol)
	    {
		// Posibility for embedded curve
		coincident = checkCoincide(intcv1.get(), start1, end1, eps,
					   intcv2.get(), clo_par1, clo_par2);
		if (coincident)
		    coincident = 2;
	    }

	    if (!coincident)
	    {
		// Check if curve two is embedded in curve one
		cv1->closestPoint(pt3, start1, end1, clo_par1, clo_pt1, clo_dist1);
		cv1->closestPoint(pt4, start1, end1, clo_par2, clo_pt2, clo_dist2);
		if (clo_dist1 <= tol && clo_dist2 <= tol)
		{
		    // Posibility for embedded curve
		    coincident = checkCoincide(intcv2.get(), start2, end2, eps,
					       intcv1.get(), clo_par1, clo_par2);
		    if (coincident)
			coincident = 3;
		}

	    }
	}

	return coincident;
    }

    int Identity::internalCoincidence(shared_ptr<ParamSurfaceInt>& intsf1, 
				      shared_ptr<ParamSurfaceInt>& intsf2, 
				      shared_ptr<GeoTol>& eps)
    {
	// Check if the first surface lies in the other.
	// The surface boundaries are already tested

	const RectDomain& domain = intsf1->getDomain();  // Surrounding parameter domain
	double umin = domain.umin();
	double umax = domain.umax();
	double vmin = domain.vmin();
	double vmax = domain.vmax();
	int nmb_sample_crvs = 10;
	double tint1 = (umax - umin)/(int)(nmb_sample_crvs+1);  // Only parameter values in the inner
	double tint2 = (vmax - vmin)/(int)(nmb_sample_crvs+1);  

	// Check coincidence with surface 2 along a number of constant parameter curves of 
	// surface 1 in both parameter directions
	// 1. parameter direction
	double par;
	int ki, kj;
	int coincident;
	shared_ptr<ParamSurface> surf1 = intsf1->getParamSurface();
	shared_ptr<ParamSurface> surf2 = intsf2->getParamSurface();
	double tol = eps->getEpsge();
	double rel_par_res = eps->getRelParRes();
	for (ki=0, par=umin+tint1; ki<nmb_sample_crvs; ++ki, par+=tint1)
	{
	    vector<shared_ptr<ParamCurve> > const_crvs = surf1->constParamCurves(par, false);
	    for (kj=0; kj<int(const_crvs.size()); ++kj)
	    {
		shared_ptr<ParamCurveInt> intcrv = 
		    shared_ptr<ParamCurveInt>(new ParamCurveInt(const_crvs[kj]));
		
		// Project the endpoints onto surface 2
		double start, end;
		Point pt1, pt2, clo_pt1, clo_pt2;
		double u1, v1, u2, v2, dist1, dist2;
		start = intcrv->startParam(0);
		end = intcrv->endParam(0);
		double del = std::min(0.1*(end-start), 0.1*tint2);
		start += del;
		end -= del;
		intcrv->point(pt1, &start);
		intcrv->point(pt2, &end);
		surf2->closestPoint(pt1, u1, v1, clo_pt1, dist1, tol);
		if (dist1 > tol)
		    return 0;  // No coincidence

		surf2->closestPoint(pt2, u2, v2, clo_pt2, dist2, tol);
		if (dist2 > tol)
		    return 0;  // No coincidence

		// Check curve
		Point parpt1(u1, v1);
		Point parpt2(u2, v2);
		coincident = checkCoincide(intcrv.get(), start, end,
					   intsf2.get(), parpt1, parpt2,
					   eps);
		if (!coincident)
		  {
		    // Extra check in closed configurations
		    if (pt1.dist(pt2) < tol && parpt1.dist(parpt2) < rel_par_res)
		      {
			double seed[2];
			Point pt3 = const_crvs[kj]->point(start + 0.1*(end-start));
			surf2->closestPoint(pt3, seed[0], seed[1], clo_pt1, dist1, tol);
			surf2->closestPoint(pt1, u1, v1, clo_pt1, dist1, tol, NULL, seed);

			Point pt4 = const_crvs[kj]->point(end - 0.1*(end-start));
			surf2->closestPoint(pt4, seed[0], seed[1], clo_pt2, dist2, tol);
			surf2->closestPoint(pt2, u2, v2, clo_pt2, dist2, tol, NULL, seed);
			if (dist1 > tol)
			  return 0;  // No coincidence

			// Check curve
			parpt1 = Point(u1, v1);
			parpt2 = Point(u2, v2);
			coincident = checkCoincide(intcrv.get(), start, end,
						   intsf2.get(), parpt1, parpt2,
						   eps);
			if (!coincident)
			  return 0;
		      }
		    else
		      return 0;
		  }
	    }
	}
	
	// 2. parameter direction
	for (ki=0, par=vmin+tint2; ki<nmb_sample_crvs; ++ki, par+=tint2)
	{
	    vector<shared_ptr<ParamCurve> > const_crvs = surf1->constParamCurves(par, true);
	    for (kj=0; kj<int(const_crvs.size()); ++kj)
	    {
		shared_ptr<ParamCurveInt> intcrv = 
		    shared_ptr<ParamCurveInt>(new ParamCurveInt(const_crvs[kj]));
		
		// Project the endpoints onto surface 2
		double start, end;
		Point pt1, pt2, clo_pt1, clo_pt2;
		double u1, v1, u2, v2, dist1, dist2;
		start = intcrv->startParam(0);
		end = intcrv->endParam(0);
		double del = std::min(0.1*(end-start), 0.1*tint1);
		start += del;
		end -= del;
		intcrv->point(pt1, &start);
		intcrv->point(pt2, &end);
		surf2->closestPoint(pt1, u1, v1, clo_pt1, dist1, tol);
		if (dist1 > tol)
		    return 0;  // No coincidence

		surf2->closestPoint(pt2, u2, v2, clo_pt2, dist2, tol);
		if (dist2 > tol)
		    return 0;  // No coincidence

		// Check curve
		Point parpt1(u1, v1);
		Point parpt2(u2, v2);
		coincident = checkCoincide(intcrv.get(), start, end,
					   intsf2.get(), parpt1, parpt2,
					   eps);
		if (!coincident)
		  {
		    // Extra check in closed configurations
		    if (pt1.dist(pt2) < tol && parpt1.dist(parpt2) < rel_par_res)
		      {
			double seed[2];
			Point pt3 = const_crvs[kj]->point(start + 0.1*(end-start));
			surf2->closestPoint(pt3, seed[0], seed[1], clo_pt1, dist1, tol);
			surf2->closestPoint(pt1, u1, v1, clo_pt1, dist1, tol, NULL, seed);

			Point pt4 = const_crvs[kj]->point(end - 0.1*(end-start));
			surf2->closestPoint(pt4, seed[0], seed[1], clo_pt2, dist2, tol);
			surf2->closestPoint(pt2, u2, v2, clo_pt2, dist2, tol, NULL, seed);
			if (dist1 > tol)
			  return 0;  // No coincidence

			// Check curve
			parpt1 = Point(u1, v1);
			parpt2 = Point(u2, v2);
			coincident = checkCoincide(intcrv.get(), start, end,
						   intsf2.get(), parpt1, parpt2,
						   eps);
			if (!coincident)
			  return 0;
		      }
		    else
		    return 0;
		  }
	    }
	}

	return 1;  // Coincidence
    }

}
