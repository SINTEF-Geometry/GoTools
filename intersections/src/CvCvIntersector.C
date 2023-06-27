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

#include "GoTools/intersections/CvCvIntersector.h"
#include "GoTools/intersections/CvPtIntersector.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/intersections/Coincidence.h"
#include "GoTools/intersections/GeoTol.h"
#include "GoTools/utils/RotatedBox.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Values.h"
#include <stdio.h> // for debugging
#include <iostream>
#include "GoTools/geometry/ObjectHeader.h" // for debugging
#include <fstream> // For debugging

using std::vector;
using std::cout;

namespace Go
{

//===========================================================================
CvCvIntersector::CvCvIntersector(shared_ptr<ParamGeomInt> curve1, 
				 shared_ptr<ParamGeomInt> curve2,
				 shared_ptr<GeoTol> epsge, 
				 Intersector* prev,
				 int eliminated_parameter,
				 double eliminated_value)
    : Intersector2Obj(curve1, curve2, epsge, prev, eliminated_parameter, eliminated_value)
//===========================================================================
{
  
}

//===========================================================================
CvCvIntersector::CvCvIntersector(shared_ptr<ParamGeomInt> curve1, 
				 shared_ptr<ParamGeomInt> curve2,
				 double epsge, Intersector* prev,
				 int eliminated_parameter,
				 double eliminated_value)
    : Intersector2Obj(curve1, curve2, epsge, prev, eliminated_parameter, eliminated_value)
//===========================================================================
{
  
}

//===========================================================================
CvCvIntersector::~CvCvIntersector()
//===========================================================================
{
    // Currently empty
}

//==========================================================================


//    shared_ptr<Intersector> 

//===========================================================================
shared_ptr<Intersector> 
CvCvIntersector::lowerOrderIntersector(shared_ptr<ParamGeomInt> obj1,
				       shared_ptr<ParamGeomInt> obj2, 
				       Intersector* parent, 
				       int eliminated_parameter,
				       double eliminated_value)
//===========================================================================
{
    shared_ptr<CvPtIntersector> curr_inter;
    curr_inter = shared_ptr<CvPtIntersector>(new CvPtIntersector(obj1,
								 obj2,
								 epsge_,
								 parent,
								 eliminated_parameter,
								 eliminated_value));
    return curr_inter;
}

///////////////////////////////////////////////////////////////////
//
//  Purpose    : Perform a rotated box test between two curves
//
//  Written by : Vibeke Skytt, 0105
//
/////////////////////////////////////////////////////////////////

int CvCvIntersector::performRotatedBoxTest(double eps1, double eps2)
{
    //double ang_tol = epsge_->getAngleTol();
    ParamCurveInt *curve[2];
    curve[0] = obj_int_[0]->getParamCurveInt();
    curve[1] = obj_int_[1]->getParamCurveInt();
    int dim = curve[0]->dimension();
  DEBUG_ERROR_IF(curve[0] == 0 || curve[1] == 0,
	   "Inconsistence in the data structure");

    // Define coordinate system for the rotation. First check if an
    // intersection point between the geometry objects is found
    vector<Point> axis(2), der1(2), der2(2);
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);
    int nmb_rec = nmbRecursions();
    int cv_ind = (nmb_rec % 2 == 0) ? 0 : 1;
    if (int_pts.size() > 0)
    {
	// Make coordinate system from the partial derivatives of the surface
	// in the (first) intersection point
	const double *param1 = (cv_ind == 0) ? int_pts[0]->getPar1() : 
	    int_pts[0]->getPar2();
	const double *param2 = (cv_ind == 1) ? int_pts[0]->getPar1() : 
	    int_pts[0]->getPar2();
	obj_int_[cv_ind]->point(der1, param1, 1);
	obj_int_[1-cv_ind]->point(der2, param2, 1);
	axis[0] = der1[1];
	if (dim == 3)
	  axis[1] = der1[1].cross(der2[1]);
    }
    else
    {
	curve[cv_ind]->axisFromEndpts(axis[0]);
	curve[1-cv_ind]->axisFromEndpts(der1[0]);
	if (dim== 3)
	  axis[1] = axis[0].cross(der1[0]);
    }
  
    // Make rotated boxes
    if (axis[0].length() < epsge_->getNumericalTol() ||
	(dim == 3 && axis[1].length() < epsge_->getNumericalTol()))
	return 1;
    axis[0].normalize();
    if (dim == 3)
      {
	axis[1].normalize();
	if (axis[0].angle_smallest(axis[1]) < epsge_->getNumericalTol())
	  return 1;
      }

    RotatedBox box1 = curve[0]->getRotatedBox(axis);
    RotatedBox box2 = curve[1]->getRotatedBox(axis);

    bool overlap = box1.overlaps(box2, eps1, 0.0);
    if (overlap)
    {
	overlap = box1.overlaps(box2, eps1, eps2);
	if (!overlap)
	{
	    if (foundIntersectionNearBoundary())
		return 1;
	}
    }
    return (overlap) ? 1 : 0;
}

///////////////////////////////////////////////////////////////////
//
//  Purpose    : Check for an intersection point close to one of the
//               endpoints of the curve provided that no intersection
//               point is found at this endpoint already
//
//  Written by : Vibeke Skytt, 0405
//
/////////////////////////////////////////////////////////////////

bool CvCvIntersector::foundIntersectionNearBoundary()
{
    double ptol = epsge_->getRelParRes();

    // Fetch parameter intervals
    double ta[2], tb[2];
    for (int kj=0; kj<2; kj++)
    {
	ta[kj] = obj_int_[kj]->startParam(0);
	tb[kj] = obj_int_[kj]->endParam(0);
    }

    // Fetch already existing intersection points
    bool has_int[4];
    has_int[0] = has_int[1] = has_int[2] = has_int[3] = false;
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    double par;
    for (size_t ki=0; ki<int_pts.size(); ki++)
    {
	for (int kj=0; kj<2; kj++)
	{
	    par = int_pts[ki]->getPar(kj);
	    if (fabs(par-ta[kj]) < ptol)
		has_int[2*kj] = true;
	    if (fabs(tb[kj]-par) < ptol)
		has_int[2*kj+1] = true;
	}
    }

    double param[2], guess[2], dist;
    for (int kj=0; kj<2; kj++)
    {
	if (has_int[2*kj]);
	else
	{
	    guess[kj] = ta[kj];
	    guess[1-kj] = 0.5*(ta[1-kj] + tb[1-kj]);

	    // Iterate to intersection point
	    doIterate(param[0], param[1], dist, guess);
	    if (dist <= epsge_->getEpsge())
		return true;
	}
	if (has_int[2*kj+1]);
	else
	{
	    guess[kj] = tb[kj];
	    guess[1-kj] = 0.5*(ta[1-kj] + tb[1-kj]);

	    // Iterate to intersection point
	    doIterate(param[0], param[1], dist, guess);
	    if (dist <= epsge_->getEpsge())
		return true;
	}
    }

    return false;
}

///////////////////////////////////////////////////////////////////
//
//  Purpose    : Do extra simple case test between two parametric curves
//               Based on the sisl function s1796.
//
//  Written by : Vibeke Skytt, 0105
//
/////////////////////////////////////////////////////////////////

int CvCvIntersector::simpleCase2(Point& axis1, Point& axis2)
{
    // The function only makes sense for spline curves. Check if we have the
    // correct input.
    int kr;
    for (kr=0; kr<2; kr++)
    {
	if (!(obj_int_[kr]->isSpline()))
	    return 0;
    }

    // We are making a plane spanned by the centre axis af the
    // cones surrounding the orientating surfaces
    // on the unit sphere. Then we project each tangent
    // to this plane and compute the angle beetween these
    // projection and the centre of the cone. If the new
    // projected cones do not overlap we have simple case.
    bool turned = false;
    double angle = axis1.angle_smallest(axis2);
    if (angle > 0.5*M_PI)
    {
	angle = M_PI - angle;
	turned = true;
    }

    // Compute the cone angles
    Point scen1, scen2;
    double tlen = axis1*axis2;
    double tangle[2];
    for (kr=0,scen1=axis1, scen2=axis2-tlen*axis1; 
	 kr<2; kr++, scen1=axis2, scen2=axis1-tlen*axis2)
    {
	scen2.normalize();
	if (turned)
	    scen2 *= -1;
	tangle[kr] = obj_int_[kr]->getOptimizedConeAngle(scen1, scen2);
    }

    // Performing a simple case check. 
    if (tangle[0] + tangle[1] < angle + epsge_->getRelParRes())
	return 1;  // Simple case
    else
	return 0;
}

//===========================================================================
int CvCvIntersector::checkCoincidence()
//===========================================================================
{
    // What about linearity?

    // Can intersection point occur in the inner of the curves?
    // In that case all inner intersection points can be fetched an
    // sorted according the the parameter belonging to one of the curves
    // Coincidence can be checked between two and two points in sequence.
    // If all intervals are coincident, then there is total coincidence.
    // Otherwise, connections are made, but the return value is no
    // coincidence.

    double ptol = epsge_->getRelParRes();
    ParamCurveInt *curve1 = obj_int_[0]->getParamCurveInt();
    ParamCurveInt *curve2 = obj_int_[1]->getParamCurveInt();
    DEBUG_ERROR_IF(curve1 == 0 || curve2 == 0, "Error in data structure");

  // Fetch intersection points at the endpoints of the two curves
  vector<shared_ptr<IntersectionPoint> > bd_ints;
  int_results_->getSortedIntersections(bd_ints);

  if (bd_ints.size() < 2)
    return 0;  // Intersection points in two endpoints are expected
  // for coincidence

  // Dismiss points which are identical in the parameter interval of one
  // curve, can occur in degenerate and closed cases
  int  ki; 
  size_t ksize = bd_ints.size();
  for (ki=0; ki<(int)(ksize)-1; ki++)
  {
      if (fabs(bd_ints[ki]->getPar(0)-bd_ints[ki+1]->getPar(0)) < ptol ||
	  fabs(bd_ints[ki]->getPar(1)-bd_ints[ki+1]->getPar(1)) < ptol)
      {
	  // Check midpoint
	  if (fabs(bd_ints[ki]->getPar(0)-bd_ints[ki+1]->getPar(0)) < ptol)
	  {
	      double mid_par = 0.5*(bd_ints[ki]->getPar(1)+bd_ints[ki+1]->getPar(1));
	      Point mid_pt;
	      curve2->point(mid_pt, &mid_par);

	      if (mid_pt.dist(bd_ints[ki]->getPoint()) > epsge_->getEpsge())
		  return 0;   // Not a degenerate case

	      double tmin, tmax;
	      shared_ptr<ParamCurve> cv
		  = curve1->getParentParamCurve(tmin, tmax);
	      tmin = bd_ints[ki]->getPar(0);
	      tmax = bd_ints[ki+1]->getPar(0);
	      if (tmax < tmin)
	      {
		  tmax = tmin;
		  tmin = bd_ints[ki+1]->getPar(0);
	      }

	      Point clo_pt;
	      double clo_par, clo_dist;
	      double cv_par1 = 0.5*(tmin + tmax);
	      cv->closestPoint(mid_pt, tmin, tmax,
			       clo_par, clo_pt, clo_dist, &cv_par1);
	      if (clo_dist > epsge_->getEpsge()) 
		  return 0;
	      else
		  return 1;
	  }
	  else
	  {
	      double mid_par = 0.5*(bd_ints[ki]->getPar(0)+bd_ints[ki+1]->getPar(0));
	      Point mid_pt;
	      curve1->point(mid_pt, &mid_par);

	      if (mid_pt.dist(bd_ints[ki]->getPoint()) > epsge_->getEpsge())
		  return 0;   // Not a degenerate case

	      double tmin, tmax;
	      shared_ptr<ParamCurve> cv
		  = curve2->getParentParamCurve(tmin, tmax);
	      tmin = bd_ints[ki]->getPar(1);
	      tmax = bd_ints[ki+1]->getPar(1);
	      if (tmax < tmin)
	      {
		  tmax = tmin;
		  tmin = bd_ints[ki+1]->getPar(1);
	      }

	      Point clo_pt;
	      double clo_par, clo_dist;
	      double cv_par1 = 0.5*(tmin + tmax);
	      cv->closestPoint(mid_pt, tmin, tmax,
			       clo_par, clo_pt, clo_dist, &cv_par1);
	      if (clo_dist > epsge_->getEpsge()) 
		  return 0;
	      else
		  return 1;
	      
	  }
      }
  }

  if (bd_ints.size() < 2)
    return 0;  // Intersection points in two endpoints are expected
  // for coincidence

  // For all intervals between intersection points
  bool coinc = true;
  int is_coincident = 0;
  if (selfint_case_ && (fabs(bd_ints[0]->getPar1()[0]-bd_ints[0]->getPar2()[0]) < ptol) &&
      (fabs(bd_ints[ksize-1]->getPar1()[0]-bd_ints[ksize-1]->getPar2()[0]) < ptol))
  {
      // Selfintersection case. Identical curves are expecte.
      // Avoid marching, but check a few points on the curves
      double ta = bd_ints[0]->getPar1()[0];
      double tb = bd_ints[ksize-1]->getPar1()[0];
      int nsample = 5;
      double tint = (tb - ta)/(double)(nsample+1);
      double tc;
      int kr;
      for (kr=0, tc=ta+tint; kr<nsample; kr++, tc+=tint)
      {
	  // Evaluate both curves
	  Point pt1, pt2;
	  curve1->point(pt1, &tc);
	  curve2->point(pt2, &tc);
	  if (pt1.dist(pt2) > epsge_->getEpsge()) 
	      break;
      }
      if (kr == nsample)
	  is_coincident = true;  // Coincidence
  }

  if (selfint_case_ && is_coincident)
  {
      // Register the coincidence by connecting the
      // intersection points
      bd_ints[0]->connectTo(bd_ints[ksize-1], COINCIDENCE_CVCV);
 }
  else
  {
      for (ki=1; ki<(int)(ksize); ki++)
      {
	  is_coincident = 0;
	  // Check if the points are connected already
	  if (bd_ints[ki-1]->isConnectedTo(bd_ints[ki]));  //@@@ Pointer or shared ptr?
	  else
	  {

	      // Perform marching to check coincidence
	      if (!is_coincident)
	      {
		  vector<double> start, end;
		  start = bd_ints[ki-1]->getPar();
		  end = bd_ints[ki]->getPar();
		  is_coincident = checkCoincide(curve1, start[0], end[0], epsge_,
					    curve2, start[1], end[1]);
	      }

	      if (is_coincident)
	      {
		  // Register the coincidence by connecting the
		  // intersection points
		  bd_ints[ki-1]->connectTo(bd_ints[ki], COINCIDENCE_CVCV);
	      }
	      else
		  coinc = false;
	  }
      }
  }

  // A coincidence interval is recognized. It is necessary to check 
  // for the possibility of additional intersections off the coincidence 
  // interval
  // @@@ VSK, Change name when implementing
  if (int_results_->checkIfBothPointsLieOnOneEndpointAndNotOfTheSameCurve())
    // @@@ or something of similar effect as a combination of a few less
    // specific calls
    return 0;
  else if (selfint_case_ == 1 && coinc)
  {
      // If loose points exist at this stage, the boundary has been checked before
      if (ksize > 2)
	  return 1;  // No more subdivision

      // Check cone for possibility of self intersection
      DirectionCone cone1, cone2;
      try {
	  cone1 = obj_int_[0]->directionCone();
	  cone2 = obj_int_[1]->directionCone();
      } catch (...) {
	  return 0;  // Continue subdividing
      }
      if (cone1.greaterThanPi() && cone2.greaterThanPi())
	  return 0;  // Possibilities for self intersection
      return 1;      // Accept coincidence interval
  }
  else 
    return coinc;
}

//===========================================================================
//
// Purpose : Two curves lie both within the same epsion ball. It is not a
//           simple case, and there is intersection.
//           This is an uwanted situation. We have to make a result that
//           is as consistent as possible since we cannot subdivide anymore.
//
// Written by : Vibeke Skytt, 0804.
//===========================================================================
void CvCvIntersector::microCase()
{
  // Fetch all intersection points belonging to the two curves
  // The intersection points should be sorted according to the parameter
  // of one of the curves
  vector<shared_ptr<IntersectionPoint> > int_pts;
  int_results_->getSortedIntersections(int_pts);

  int ki;
  if (int_pts.size() == 0)
    {
      // No intersection point exist. Construct one internal to the two
      // curves
      // It is probably good enough to take the midpoints since the 
      // curves are so small
      double mid_par1 = 0.5*(obj_int_[0]->startParam(0) + 
			     obj_int_[0]->endParam(0));
      double mid_par2 = 0.5*(obj_int_[1]->startParam(0) + 
			     obj_int_[1]->endParam(0));

      double guess[2];
      guess[0] = 0.5*(obj_int_[0]->startParam(0) + obj_int_[0]->endParam(0));
      guess[1] = 0.5*(obj_int_[1]->startParam(0) + obj_int_[1]->endParam(0));
      double par1, par2, dist;
      doIterate(par1, par2, dist, guess);
      if (dist < getTolerance()->getEpsge()) {
	  shared_ptr<IntersectionPoint> tmp = 
	      int_results_->addIntersectionPoint(obj_int_[0], 
						 obj_int_[1], 
						 getTolerance(),
						 &mid_par1, 
						 &mid_par2); 
      }
    }
  else if (int_pts.size() == 1)
    {
      // One intersection point. Leave it like this.
      ;
    }
  else
    {
      // More than one intersection point. Connect.
      for (ki=1; ki<int(int_pts.size()); ki++)
	  if (int_pts[ki-1]->isConnectedTo(int_pts[ki]))
	      break;

      if (ki < int(int_pts.size()))
	int_pts[0]->connectTo(int_pts[int_pts.size()-1], MICRO_CVCV);
    }
}

//===========================================================================
//
// Purpose : Given two parametric curve in a simple case situation, iterate
//           to the intersection point, if any. 
//
// Written by : Vibeke Skytt, 0804
//
//===========================================================================
int CvCvIntersector::updateIntersections()
//===========================================================================
{
    // Fetch already existing intersection points
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    // VSK, 0609. Post iterate existing intersection points
    int ki, kj;
    double param[2], guess[2];
    double dist;
    double ptol = epsge_->getRelParRes();
    double ta[2], tb[2];
    ta[0] = obj_int_[0]->startParam(0);
    ta[1] = obj_int_[1]->startParam(0);
    tb[0] = obj_int_[0]->endParam(0);
    tb[1] = obj_int_[1]->endParam(0);
    for (ki=0; ki<(int)(int_pts.size()); ki++)
    {
	guess[0] = int_pts[ki]->getPar(0);
	guess[1] = int_pts[ki]->getPar(1);
	if (fabs(guess[0]-ta[0])<ptol || fabs(guess[0]-tb[0]) < ptol ||
	    fabs(guess[1]-ta[1])<ptol || fabs(guess[1]-tb[1]) < ptol)
	    continue;  // Do not move boundary points

	doIterate(param[0], param[1], dist, guess);
	for (kj=0; kj<2; kj++)
	    if ((fabs(guess[kj]-ta[kj])<ptol || fabs(guess[kj]-tb[kj])) &&
		!(fabs(param[kj]-ta[kj])<ptol || fabs(param[kj]-tb[kj])))
		break;
	if (dist < int_pts[ki]->getDist() && kj==2)
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
	for (kj=0; kj<2; kj++)
	{
	    double t1 = int_pts[0]->getPar(kj);
	    if (fabs(t1-ta[kj])<ptol || fabs(t1-tb[kj])<ptol)
		break;
	}
	at_bd1 = (kj<2);

	for (ki=1; ki<(int)(int_pts.size()); ki++)
	{
	    for (kj=0; kj<2; kj++)
	    {
		double t1 = int_pts[ki]->getPar(kj);
		if (fabs(t1-ta[kj])<ptol || fabs(t1-tb[kj])<ptol)
		    break;
	    }
	    at_bd2 = (kj<2);

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
		

    if (int_pts.size() > 0)
	{
	    // At most one intersection point is expected. One or more points
	    // exist already. Leave it like this.
	    return 0;
	}

    // Iterate to an intersection point between the two curves
    ptol = 100.0*epsge_->getRelParRes();
    double par1, par2;
    doIterate(par1, par2, dist);
    if (dist <= epsge_->getEpsge())
	{
	if (getenv("DEBUG_ITER") && (*getenv("DEBUG_ITER"))=='1') {
	    cout << "CvCv. Int. pt. found " << par1 << " " << par2;
	    cout << ", dist: " << dist;
	}
	    // An intersection point is found. Represent it in the data
	    // structure
	    // @@@ vsk, In sisl there is a check if the found intersection
	    // point is very close to the endpoints of the curve. In that case
	    // it is dismissed since the intersections at the endpoints are found
	    // already. Here we have already checked if there exist any
	    // intersection points. Thus, we know that there does not exist any
	    // intersection point at the endpoint.
	    // VSK, 0906 But maybe the endpoint already is in intersection, but that
	    // point is represented elsewhere at a more exact position. Make a check.
	    if (par1-obj_int_[0]->startParam(0) < ptol || 
		obj_int_[0]->endParam(0)-par1 < ptol ||
		par2-obj_int_[1]->startParam(0) < ptol || 
		obj_int_[1]->endParam(0)-par2 < ptol)
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
					       &par1, 
					       &par2);
	    return 1; // A new point is found
	}

  return 0;
}

void CvCvIntersector::doIterate(double& par1, double& par2, double& dist,
				double *guess)
{
    // Iterate to a closest point between the two curves
    // First fetch the parametric curves. We need both the current
    // parametric curve to compute good start values for the iteration
    // and the initial curve where no numerical noice is added trough
    // subdivision to get a numerically stable result from the
    // iteration.
    double start1, end1;  // The parametric limitation of the current
    // curve in the parameter interval of the parent curve
    double start2, end2;  // Similar for curve 2
    ParamCurveInt* cv_int1 = obj_int_[0]->getParamCurveInt();
    ParamCurveInt* cv_int2 = obj_int_[1]->getParamCurveInt();
    shared_ptr<ParamCurve> curve1 = cv_int1->getParentParamCurve(start1, end1);
    shared_ptr<ParamCurve> curve2 = cv_int2->getParentParamCurve(start2, end2);
    shared_ptr<ParamCurve> cv1 = cv_int1->getParamCurve();
    shared_ptr<ParamCurve> cv2 = cv_int2->getParamCurve();
    DEBUG_ERROR_IF(curve1.get() == 0 || curve2.get() == 0 || 
	     cv1.get() == 0 || cv2.get() == 0,
	     "Inconsistence in the data structure");

    // Compute a startpoint for the iteration
    double seed[2];
    if (guess)
    {
	seed[0] = guess[0];
	seed[1] = guess[1];
    }
    else
	getSeedIteration(seed);

    // Iterate
    Point ptc1, ptc2;

    ClosestPoint::closestPtCurves(curve1.get(), curve2.get(), start1, end1, start2, end2, 
		    seed[0], seed[1], par1, par2, dist, ptc1, ptc2);
}

//===========================================================================
double CvCvIntersector::distInCandidatePar(double par, int dir, const double* seed)
//===========================================================================
{
    double start, end;  // The parametric limitation of the current
                          // curve in the parameter interval of the
                          // parent curve
    shared_ptr<ParamCurve> curve = 
	obj_int_[dir]->getParamCurveInt()->getParentParamCurve(start, end);
    DEBUG_ERROR_IF(curve.get() == 0,
		   "Inconsistence in the data structure");
    double guess;
    int ki;
    for (ki=0; ki<2; ki++)
    {
	if (ki == dir)
	    continue;
	guess = seed[ki];
    }

    // Evaluate one curve and iterate to the closest point in the `other curve
    Point pt;
    obj_int_[1-dir]->point(pt, &par);

    double clo_t, clo_dist;
    Point clo_pt;
    curve->closestPoint(pt, start, end, clo_t, clo_pt, clo_dist, &guess);

    return clo_dist;
}

////===========================================================================
//
// Purpose : Given two linear parametric curves, compute the intersection.
//
// Written by : Vibeke Skytt, 1204
//
//===========================================================================
int CvCvIntersector::linearCase()
//===========================================================================
{
    int ki, kj;

    // First check existance of any intersection points
    int dir = 0;
    vector<shared_ptr<IntersectionPoint> > int_pts;
    int_results_->getIntersectionPoints(int_pts);

    if (int_pts.size() == 1)
	return 1;   // One intersection point already computed

    if (int_pts.size() > 1)
    {
	shared_ptr<IntersectionPoint> dum;
	// Sort according to the parameter direction of the first curve
	for (ki = 0; ki < int(int_pts.size()); ki++)
	    for (kj = ki + 1; kj < int(int_pts.size()); kj++)
		if (int_pts[kj]->getPar(dir) < int_pts[ki]->getPar(dir))
		{
		    dum = int_pts[kj];
		    int_pts[kj] = int_pts[ki];
		    int_pts[ki] = dum;
		}

	// Make sure that all intersection points are connected
	for (ki = 1; ki < int(int_pts.size()); ki++)
	    if (!(int_pts[ki-1]->isConnectedTo(int_pts[ki])))
		int_pts[ki-1]->connectTo(int_pts[ki], LINEAR_CVCV);
	return 1;
    }

    // Try to iterate for an intersection point
    double par1, par2, dist;
    doIterate(par1, par2, dist);
 
  if (dist <= epsge_->getEpsge())
    {
      // An intersection point is found. Represent it in the data
      // structure

	shared_ptr<IntersectionPoint> tmp = 
	    int_results_->addIntersectionPoint(obj_int_[0], 
					 obj_int_[1], 
					 getTolerance(),
					 &par1, 
					 &par2);
      return 1; // A new point is found
    }

  return 0;  
}


//===========================================================================
//
// Purpose : Given two parametric curve and no simple case situation, 
//           subdivide the curves to produce more simple sub problems.
//           Perform intersection with the subdivision points.
//           Prepare for the next recursion level.
//
// Written by : Vibeke Skytt, 0804
//
//===========================================================================
int CvCvIntersector::doSubdivide()
//===========================================================================
{
  int perm[2];  // Two curves means to parameter directions that need sorting
  int nmb_subdiv;  // Number of parameter directions in which to subdivide

  // Intersection objects belonging to the next recursion level
  vector<shared_ptr<ParamGeomInt> > sub_objects1;
  vector<shared_ptr<ParamGeomInt> > sub_objects2;

  if (getenv("SUBDIV_CVCV") && *getenv("SUBDIV_CVCV") == '1')
  {
      int npar1=1, npar2=1;
	std::cout << "================================================" << std::endl;
	std::cout << "Domain 1: ";
	std::cout << obj_int_[0]->startParam(0) << " ";
	std::cout << obj_int_[0]->endParam(0) << " ";
	std::cout << std::endl;
	std::cout << "Domain 2: ";
	std::cout << obj_int_[1]->startParam(0) << " ";
	std::cout << obj_int_[1]->endParam(0) << " ";
	std::cout << std::endl;

	vector<shared_ptr<IntersectionPoint> > ipoint;
	int_results_->getIntersectionPoints(ipoint);
	std::cout << "Intersection points in this pool: "
		  << ipoint.size() << std::endl;
	for (int ki=0; ki < int(ipoint.size()); ki++) {
	    for (int kj=0; kj<npar1+npar2; kj++) {
		std::cout << ipoint[ki]->getPar(kj) << "  ";
	    }
	    Point p1 = ipoint[ki]->getPoint1();
	    Point p2 = ipoint[ki]->getPoint1();
	    for (int kj=0; kj<p1.dimension(); kj++)
		std::cout << p1[kj] << "  ";
	    for (int kj=0; kj<p2.dimension(); kj++)
		std::cout << p2[kj] << "  ";
	    std::cout << ipoint[ki]->getDist() << "  ";
	    std::cout << ipoint[ki]->numNeighbours();
	    std::cout << std::endl;
	}
  }

  // Sort the curves according to importance of subdivision according to
  // properties of the curves and already computed intersections
  nmb_subdiv = sortParameterDirections(perm);

  if (nmb_subdiv == 0)
    return 0;   // Not possible to subdivide any of the curves

  // For each parameter direction in prioritized order, fetch an
  // appropriate subdivision parameter, and perform subdivision
  int ki, kj;
  double subdiv_par;
  SubdivisionClassification found;
  vector<shared_ptr<ParamGeomInt> > curve_sub;
  vector<shared_ptr<ParamGeomInt> > subdivpt;
  for (ki=0; ki<nmb_subdiv; ki++)
    {
      found = getSubdivisionParameter(perm[ki], subdiv_par);
      if (found == CANNOT_SUBDIVIDE)
	{
	  // No parameter value is found. Move this parameter direction
	  // to the end of the permutation array, and decrease the
	  // number of parameter directions where it is possible to 
	  // subdivide.
	    std::swap(perm[ki], perm[nmb_subdiv-1]);
	    nmb_subdiv--;
	    ki--;
	    continue;
	}

    // !!! DEBUG
  if (getenv("SUBDIV_CVCV") && *getenv("SUBDIV_CVCV") == '1')
  {
      std::cout << "Subdivide dir = " << perm[ki] << " par = " << subdiv_par;
      std::cout << " criterium = " << found << std::endl;
  }
     // Subdivide the current curve
      subdivpt.clear();
      curve_sub.clear();
      try {
	  obj_int_[perm[ki]]->subdivide(0, subdiv_par, curve_sub, subdivpt);
      } catch (...)
      {
	  subdivpt.clear();
	  curve_sub.clear();
      }
	  
	  if (subdivpt.size() < 1 || curve_sub.size() == 0)
	    continue;  // No new objects 

      for (kj = 0; kj < int(subdivpt.size()); kj++)
	{
	  // Intersect the subdivision points with the other object
	  shared_ptr<CvPtIntersector> subdiv_intersector = 
	      shared_ptr<CvPtIntersector>
	      (new CvPtIntersector(perm[ki]==0 ? subdivpt[kj] : obj_int_[0],
				   perm[ki]==0 ? obj_int_[1] : subdivpt[kj],
				   epsge_, this, perm[ki], subdiv_par));


	  // Is it here relevant to fetch existing intersection points
	  // and/or insert points into intersection curves before computing
	  // intersection points with the subdivision points. Normally, this
	  // is mostly of interest when surfaces are involved, but intersection
	  // intervals might exist?
	  // These computations do anyway involve the intersection pool, but
	  // they might be trigged somehow. The parameter direction and value
	  // are required information

	  if (found == DIVIDE_SING && hasSingularityInfo())
	  {
	      // Transfer singularity information to the child
	      // process
		  subdiv_intersector->setSingularityInfo(getSingularityInfo(), 
							 perm[ki]);
	  }

	  subdiv_intersector->compute(false);

	  // Check quality of intersection points
	  // If the quality is not sufficient, find a new subdivision parameter
	  // and repeat the process.
	  // Otherwise
	  int_results_->includeReducedInts(subdiv_intersector->getIntPool());
	}

      if (perm[ki] == 0)
	sub_objects1.insert(sub_objects1.end(), curve_sub.begin(), 
			    curve_sub.end());
      else
	sub_objects2.insert(sub_objects2.end(), curve_sub.begin(), 
			    curve_sub.end());
    }

  if (nmb_subdiv == 0)
      return 0;  // Not subdivided

  if (sub_objects1.size() == 0)
      sub_objects1.push_back(obj_int_[0]);
  if (sub_objects2.size() == 0)
      sub_objects2.push_back(obj_int_[1]);

  // Create new intersector objects
  for (ki = 0; ki < int(sub_objects1.size()); ki++)
    for (kj = 0; kj < int(sub_objects2.size()); kj++)
      {
	shared_ptr<Intersector> intersector = 
	  shared_ptr<Intersector>(new CvCvIntersector(sub_objects1[ki],
						      sub_objects2[kj],
						      epsge_, 
						      this));
	//intersector->getIntPool()->setPoolInfo(int_results_);
	sub_intersectors_.push_back(intersector);
      }

  return 1;
}

//===========================================================================
//
// Purpose : Fetch information related to the curves in order to decide
//           which curves that can be subdivided and what curve is most
//           important to subdivide.
//
// Written by : Vibeke Skytt, 0804
//
//===========================================================================
int CvCvIntersector::sortParameterDirections(int perm[])
//===========================================================================
{
  double length[2], wiggle[2];  // Estimated length and wiggliness
  bool inner_knots[2], critical_val[2], can_divide[2];
  bool has_inner_ints[2];
  double rel_wiggle = 0.1;
  double rel_length = 0.1;

  // Number of parameter directions is two.
  int size = 2;

  // Fetch information from the curves
  int ki, kj;
  for (ki=0; ki<size; ki++)
    {
      obj_int_[ki]->getLengthAndWiggle(length+ki, wiggle+ki);

      inner_knots[ki] = obj_int_[ki]->hasInnerKnots(0);

      critical_val[ki] = obj_int_[ki]->hasCriticalVals(0);

      can_divide[ki] = obj_int_[ki]->canDivide(0);

      has_inner_ints[ki] = int_results_->hasPointsInInner(ki);
    }

  double med_len = 0.0, med_wiggle = 0.0;
  for (ki=0; ki<2; ki++)
  {
      med_len += length[ki];
      med_wiggle += wiggle[ki];
  }
  med_len /= (double)2;
  med_wiggle /= (double)2;

  double min_length = std::max(0.5*med_len, epsge_->getEpsge());
  double min_wiggle = std::max(0.5*med_wiggle, 0.02);

  // Fetch information from the intersection pool according to inner 
  // intersection points


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
//
// Purpose :
//
// Written by : Vibeke Skytt, 0804
//
//===========================================================================
SubdivisionClassification CvCvIntersector::getSubdivisionParameter(int dir, double& par)
//===========================================================================
{
  int ki, kj, kr, sgn;
  double aeps = epsge_->getEpsge();
  double ptol = 100.0*epsge_->getRelParRes();
  double gtol = 100.0*aeps;
  double pval;
  double frac = 0.2;
  bool has_candidate = false;
  double candidate_div;

  // Set pointer to the intersection object corresponding to the parameter
  // direction
  ParamGeomInt *obj = (dir == 0) ? obj_int_[0].get() : obj_int_[1].get();
  double ta = obj->startParam(0);
  double tb = obj->endParam(0);
  double frac2 = 0.01;
  ptol = std::max(ptol, frac2*(tb-ta));

  // Get critical parameters
  // @@@ Critical parameters of different priority? Sorting?
  vector<double> critical_pars = obj->getCriticalVals(0);

  int size = (int)critical_pars.size();
  int is_critical;
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
	      for (kj = 0; kj < int(ints_in_area.size()); kj++)
	      {
		  pval = ints_in_area[kj]->getPar(dir);
		  if (fabs(par - pval) < frac*(tb - ta) && fabs(par-pval) > epsge_->getRelParRes())
		  {
		      is_critical = 1;
		      break;
		  }
		  else if (distInCandidatePar(pval, dir, 
					      ints_in_area[kj]->getPar1()) >= 0.9*epsge_->getEpsge())
		  {
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
  // Look for a suitable knot in which to subdivide. First fetch
  // the knots sorted according to the distance to the mid-parameter
  vector<double> knot_vals = obj->getInnerKnotVals(0, true);
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
	      for (kj = 0; kj < int(ints_in_area.size()); kj++)
	      {
		  pval = ints_in_area[kj]->getPar(dir);
		  if (fabs(par - pval) < frac*(tb - ta) && fabs(par-pval) > epsge_->getRelParRes())
		  {
		      is_critical = 1;
		      break;
		  }
		  else if (distInCandidatePar(pval, dir, 
					      ints_in_area[kj]->getPar1()) >= 0.9*epsge_->getEpsge())
		  {
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
  if (singularity_info_.get() == 0)
    {
	// Make empty singularity info instance
	singularity_info_ = (shared_ptr<SingularityInfo>)(new SingularityInfo());
    }
  
  // Check if a closest point exist already
  double param[2], dist = 2.0*gtol;
  if (singularity_info_->hasPoint())
  {
      par = singularity_info_->getParam(dir);
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
      doIterate(param[0], param[1], dist);
      if (dist > 0.9*aeps && dist <= aeps)
	  dist = 2.0*gtol;  // Avoid subdivision in a bad intersection point
      if (dist < gtol && param[dir] > ta+ptol && param[dir] < tb-ptol)
	  singularity_info_->setSingularPoint(param, 2);
      par = param[dir];
  }

  if (dist < gtol && par > ta+ptol && par < tb-ptol)
  {
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

void CvCvIntersector::writeOut()
{
      ParamCurveInt *curve1 = obj_int_[0]->getParamCurveInt();
      ParamCurveInt *curve2 = obj_int_[1]->getParamCurveInt();

      SplineCurve *qc1 = curve1->getParamCurve()->geometryCurve();
      SplineCurve *qc2 = curve2->getParamCurve()->geometryCurve();
    
      std::ofstream debug("touch_case.g2");
      qc1->writeStandardHeader(debug);
      qc1->write(debug);
      qc2->writeStandardHeader(debug);
      qc2->write(debug);
}
} // namespace Go
