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

#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/SplineDebugUtils.h"

#include <fstream>


using namespace Go;
using std::vector;

// typedef struct intersection_point {
//   double par1, par2;
//   int pretop[4];

//   intersection_point(double p1, double p2, int *top)
// 	{ 
// 	  par1 = p1; par2 = p2; 
// 	  pretop[0] = top[0]; pretop[1] = top[1]; pretop[2] = top[2]; pretop[3] = top[3];
// 	}

// } intersection_point;

// Comparisement function to use in std::sort
// static bool par1_compare(const intersection_point& el1,
// 			 const intersection_point& el2)
// {
//    if (el1.par1 < el2.par1)
//       return true;
//    else
//       return false;
// }


//===========================================================================
void LoopUtils::
representAsSurfaceCurves(const std::vector< shared_ptr<ParamCurve> >& curves,
			 shared_ptr<BoundedSurface> surf,
			 std::vector<shared_ptr<CurveOnSurface> >& cvs_on_sf)
//===========================================================================
{
    cvs_on_sf.resize(curves.size());

    shared_ptr<ParamSurface> under_sf = surf->underlyingSurface();

    int sfdim = surf->dimension();
    for (size_t ki=0; ki<curves.size(); ++ki)
    {
	shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curves[ki]);
	if (sfcv.get() != 0)
	{
	    ALWAYS_ERROR_IF(sfcv->underlyingSurface() != under_sf,
			    "Inconsistent surface pointers");
	    cvs_on_sf[ki] = sfcv;
	}
	else
	{
	    int cvdim = curves[ki]->dimension();
	    if (cvdim == sfdim && sfdim > 2)
	    {
                cvs_on_sf[ki] = shared_ptr<CurveOnSurface>(new CurveOnSurface(under_sf, curves[ki], false));
	    }
	    else if (cvdim == 2)
	    {
		cvs_on_sf[ki] = shared_ptr<CurveOnSurface>(new CurveOnSurface(under_sf, curves[ki], true));
	    }
	}
    }

}

//===========================================================================
bool
LoopUtils::loopIsCCW(const vector<shared_ptr<SplineCurve> >& simple_par_loop, 
		     double space_epsilon, double int_tol)
//===========================================================================
{
  vector<shared_ptr<ParamCurve> > tmp_cvs(simple_par_loop.begin(), simple_par_loop.end());
  return loopIsCCW(tmp_cvs, space_epsilon, int_tol);
}

//===========================================================================
bool
LoopUtils::loopIsCCW(const vector<shared_ptr<ParamCurve> >& simple_par_loop, 
		     double space_epsilon, double int_tol)
//===========================================================================
{
    ALWAYS_ERROR_IF(simple_par_loop.size() == 0,
		"Empty input vector!");

    int ki;
    ALWAYS_ERROR_IF(simple_par_loop[0]->dimension() != 2,
		"Input loop must be 2-dimensional.");

    // Select longest curve
    int idx = 0;
    double max_len = simple_par_loop[idx]->estimatedCurveLength();
    for (ki=1; ki<int(simple_par_loop.size()); ki++)
      {
	double len = simple_par_loop[ki]->estimatedCurveLength();
	if (len > max_len)
	  {
	    max_len = len;
	    idx = ki;
	  }
      }

    // We choose the parameter value close to the middle on the chosen curve in the loop.
    double tpar =
      0.4*simple_par_loop[idx]->startparam() + 0.6*simple_par_loop[idx]->endparam();

    vector<Point> pnt(2);
    simple_par_loop[idx]->point(pnt, tpar, 1);

    // We compute the normal in the given point, i.e. the left normal.
    Point normal(-pnt[1][1], pnt[1][0]);

    // We locate a point which is outside loop, in normal direction from pnt[0].
    BoundingBox box1 = simple_par_loop[0]->boundingBox();
    for (ki=1; ki<int(simple_par_loop.size()); ki++)
      {
	BoundingBox box2 = simple_par_loop[ki]->boundingBox();
	box1.addUnionWith(box2);
      }

    double length = box1.low().dist(box1.high()); 
    // We make sure that we get to the outside of the loop.
    length = std::min(length+1.0, 1.5*length); 
    normal.normalize();
    Point end_pt = pnt[0] + length*normal;
    SplineCurve normal_curve = SplineCurve(pnt[0], end_pt);

    // Adjust tolerance for small loops
    double minlen = std::min(box1.high()[0]-box1.low()[0],
			     box1.high()[1]-box1.high()[1]);
    double int_tol2 = std::max(1.0e-6, std::min(int_tol, 0.01*minlen));

    // We then check for intersections between normal_curve and simple_par_loop,
    // not counting start point of normal_curve.    
    vector<double> params_interval;
    //vector<shared_ptr<ParamCurve> > par_loop;
    // for (size_t ki=0; ki<simple_par_loop.size(); ++ki)
    //   par_loop.push_back(simple_par_loop[ki]);
    shared_ptr<CurveLoop> loop = 
      shared_ptr<CurveLoop>(new CurveLoop(simple_par_loop, space_epsilon));
    CurveBoundedDomain loop_dom(loop);

    loop_dom.findPcurveInsideSegments(normal_curve, int_tol2, params_interval);

    int nmbpoint = (int)params_interval.size();

  // Remove first intersection. If odd # intersections, loop is CCW.
  if ((params_interval.size() > 0) && 
      (params_interval[0] == normal_curve.startparam())) 
    nmbpoint--;

  return ((nmbpoint % 2 == 1) ? true : false);
}


//===========================================================================
bool LoopUtils::paramIsCCW(const vector< shared_ptr<CurveOnSurface> >& loop, 
			   double space_epsilon, double int_tol)
//===========================================================================
{
    if (loop.empty())
	return true;
    vector< shared_ptr<ParamCurve> > sc;
    for (size_t i = 0; i < loop.size(); ++i) {
      loop[i]->ensureParCrvExistence(int_tol);
	ParamCurve* pcptr = loop[i]->parameterCurve().get();
	ALWAYS_ERROR_IF(pcptr == 0, "Parameter curve not found");
	shared_ptr<ParamCurve> pc(loop[i]->parameterCurve());
	// shared_ptr<SplineCurve>
	//     spc(dynamic_pointer_cast<SplineCurve, ParamCurve>(pc));
	if (!pc.get())
	  pc = loop[i]->spaceCurve();
	sc.push_back(pc);
    }

    return loopIsCCW(sc, space_epsilon, int_tol);
}


//===========================================================================
bool LoopUtils::loopIsCCW(const CurveLoop& loop, double int_tol)
//===========================================================================
{
    // We extract the 2D parameter curves.  Well, come to think of it,
    // what really matters is not that the curve is a 2D-spline, but
    // 2D ...
    vector< shared_ptr<ParamCurve> > sc;
    for (int ki = 0; ki < loop.size(); ++ki) {
	shared_ptr<ParamCurve> pcv;
	/*if (loop[ki]->instanceType() == Class_SplineCurve)
	    pcv = dynamic_pointer_cast<SplineCurve, ParamCurve>(loop[ki]);
	    else */if (loop[ki]->instanceType() == Class_CurveOnSurface) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(loop[ki]);
	    if (cv_on_sf->parameterCurve().get() == NULL)
		THROW("Method requires parameter curve!");
	    pcv = cv_on_sf->parameterCurve();
	    // if (cv_on_sf->parameterCurve()->instanceType() == Class_SplineCurve)
	    // 	pcv = dynamic_pointer_cast<SplineCurve, ParamCurve>
	    // 	    (cv_on_sf->parameterCurve());
	    // else {
	    // 	pcv = shared_ptr<SplineCurve>
	    // 	    (cv_on_sf->parameterCurve()->geometryCurve());
	    // 	if (pcv.get() == NULL)
	    // 	    THROW("Unexpected incident.");
	}
        else
	    pcv = loop[ki];
	    //THROW("Unexpected incident.");

	sc.push_back(pcv);
    }

    double space_epsilon = loop.getSpaceEpsilon();
    return loopIsCCW(sc, space_epsilon, int_tol);
}


//===========================================================================
bool
LoopUtils::firstLoopInsideSecond(const vector<shared_ptr<CurveOnSurface> >& first_loop,
				 const vector<shared_ptr<CurveOnSurface> >& second_loop,
				 double loop_tol, double int_tol)
//===========================================================================
{
   // @@sbr If first space_curve is degenerate, move on to next!
   ALWAYS_ERROR_IF(first_loop.size() == 0 || second_loop.size() == 0,
		   "zero sized loop encountered");

   // We're also expecting all cvs to lie on the same sfs. Routine not
   // depending on it.
   int ki, kj;
   double a_tol = 1.0e-8;  // To identify identity of end parameters
   // Should relate parameter tolerance to geometry tolerance

   // As we're creating bd domains we must make sure loop defines bd area!
   // We're not interested in the direction anyway.
   vector<shared_ptr<CurveOnSurface> > first_ccw_loop = first_loop;
   if (!paramIsCCW(first_ccw_loop, loop_tol, int_tol)) {
       for (ki = 0; ki < int(first_ccw_loop.size()); ++ki) {
	   first_ccw_loop[ki] =
	       shared_ptr<CurveOnSurface>(first_loop[ki]->clone());
	   first_ccw_loop[ki]->reverseParameterDirection();
       }
       std::reverse(first_ccw_loop.begin(), first_ccw_loop.end());
   }
   vector<shared_ptr<CurveOnSurface> > second_ccw_loop = second_loop;
   if (!paramIsCCW(second_ccw_loop, loop_tol, int_tol)) {
       for (ki = 0; ki < int(second_ccw_loop.size()); ++ki) {
	   second_ccw_loop[ki] =
	       shared_ptr<CurveOnSurface>(second_loop[ki]->clone());
	   second_ccw_loop[ki]->reverseParameterDirection();
       }
       std::reverse(second_ccw_loop.begin(), second_ccw_loop.end());
   }

   vector<shared_ptr<ParamCurve> > pcvs1;
   for (ki = 0; ki < int(first_ccw_loop.size()); ++ki) {
      ALWAYS_ERROR_IF(first_ccw_loop[ki]->parameterCurve().get() == 0, 
		      "parameter curve not found");
      pcvs1.push_back(first_ccw_loop[ki]->parameterCurve());
   }
   shared_ptr<CurveLoop> cv_loop1(new CurveLoop(pcvs1, loop_tol));
   CurveBoundedDomain bd_domain1(cv_loop1);
   vector<shared_ptr<ParamCurve> > pcvs2;
   for (ki = 0; ki < int(second_ccw_loop.size()); ++ki) {
      ALWAYS_ERROR_IF(second_ccw_loop[ki]->parameterCurve().get() == 0,
		      "parameter curve not found");
      pcvs2.push_back(second_ccw_loop[ki]->parameterCurve());
   }
   shared_ptr<CurveLoop> cv_loop2(new CurveLoop(pcvs2, loop_tol));
   CurveBoundedDomain bd_domain2(cv_loop2);

   // From first_ccw_loop we pick a point which is in the interior of
   // the loop.  We start by picking parameter value which is not a
   // knot, thus guaranteeing us a well defined normal along curve (in
   // the parameter domain).
   shared_ptr<SplineCurve> pcv;
   pcv = dynamic_pointer_cast<SplineCurve, ParamCurve>
       (first_ccw_loop[0]->parameterCurve());
   if (pcv.get() == 0)
     {
       pcv = shared_ptr<SplineCurve>(first_ccw_loop[0]->parameterCurve()->geometryCurve());
     }

   ALWAYS_ERROR_IF(pcv.get() == 0, "Unable to convert curve to SplineCurve");
   vector<double>::const_iterator iter = pcv->basis().begin();
   int nmb_coefs = pcv->numCoefs();
   nmb_coefs = nmb_coefs/2;
   iter += nmb_coefs;
   while (iter[0] == iter[1]) {
       ++iter;
   }
   double tpar = 0.5*(iter[0] + iter[1]);
   vector<Point> par_pt = pcv->ParamCurve::point(tpar, 1);
   Vector2D par_vec(par_pt[0][0], par_pt[0][1]);
   if (!bd_domain2.isInDomain(par_vec, int_tol)) {
       return false;
   } else {
       if (!bd_domain2.isOnBoundary(par_vec, int_tol)) {
	   return true;
       } else {
	   shared_ptr<ParamSurface> under_sf =
	     first_ccw_loop[0]->underlyingSurface();
	   vector<vector<shared_ptr<CurveOnSurface> > > trim_pieces(2);
	   par_pt[1].normalize();
	   int dir = (fabs(par_pt[1][0]) < fabs(par_pt[1][1])) ? 1 : 2;
	   int par_ind = (dir == 1) ? 1 : 0;
	   double par = par_pt[0][par_ind];
	   int other_par_ind = (dir == 1) ? 0 : 1;
	   double other_par = par_pt[0][other_par_ind];
	   try {
	     bd_domain1.clipWithDomain(dir, par, int_tol,
				       under_sf, trim_pieces[0]);
	     bd_domain2.clipWithDomain(dir, par, int_tol,
				       under_sf, trim_pieces[1]);
	   }
	   catch (...)
	     {
	       dir = 3 - dir;
	       par_ind = (dir == 1) ? 1 : 0;
	       par = par_pt[0][par_ind];
	       other_par_ind = (dir == 1) ? 0 : 1;
	       other_par = par_pt[0][other_par_ind];
	       trim_pieces[0].clear();
	       trim_pieces[1].clear();
	       bd_domain1.clipWithDomain(dir, par, int_tol,
					 under_sf, trim_pieces[0]);
	       bd_domain2.clipWithDomain(dir, par, int_tol,
					 under_sf, trim_pieces[1]);
	     }

	   // We then must locate the two pieces which contains par_pt[0].
	   vector<int> ind(2);
	   vector<double> end_other_par(2);
	   for (ki = 0; ki < int(trim_pieces.size()); ++ki) {
	       for (kj = 0; kj < int(trim_pieces[ki].size()); ++kj) {
		   Point start_pt = trim_pieces[ki][kj]->parameterCurve()->point
		       (trim_pieces[ki][kj]->parameterCurve()->startparam());
		   Point end_pt = trim_pieces[ki][kj]->parameterCurve()->point
		       (trim_pieces[ki][kj]->parameterCurve()->endparam());
		   if ((start_pt.dist(par_pt[0]) < int_tol) ||
		       (end_pt.dist(par_pt[0]) < int_tol)) {
		       ind[ki] = kj;
		       end_other_par[ki] = (start_pt.dist(par_pt[0]) <
					    (end_pt.dist(par_pt[0]))) ?
			   end_pt[other_par_ind] : start_pt[other_par_ind];
		       break;
		   }
	       }
	   }
	   ASSERT(ind[0] != -1 && ind[1] != -1);

	   if (((end_other_par[0] < other_par-a_tol) &&
		(end_other_par[0] > end_other_par[1]+a_tol)) ||
	       ((end_other_par[0] > other_par+a_tol) &&
		(end_other_par[1] > end_other_par[0]+a_tol))) {
	       return true;
	   } else {
	       return false;
	   }
       }
   }
}

//===========================================================================
bool
LoopUtils::makeLoopCCW(vector<shared_ptr<ParamCurve> >& loop_cvs,
		       double tol)
//===========================================================================
{
  // The curves are assumed to have correct sequence, but possible
  // wrong parameter direction
  // Turn parameter direction of curves if necessary
  for (size_t ki=0; ki<loop_cvs.size()-1; ++ki)
    {
      Point pos = loop_cvs[ki]->point(loop_cvs[ki]->endparam());
      Point pos1 = loop_cvs[ki+1]->point(loop_cvs[ki+1]->startparam());
      Point pos2 = loop_cvs[ki+1]->point(loop_cvs[ki+1]->endparam());
      double dist1 = pos.dist(pos1);
      double dist2 = pos.dist(pos2);
      if (ki==0 && dist1>tol && dist2>tol)
	{
	  // Check if the first curve should be turned
	  Point pos0 = loop_cvs[0]->point(loop_cvs[0]->startparam());
	  double dist3 = pos0.dist(pos1);
	  double dist4 = pos0.dist(pos2);
	  if (std::min(dist3,dist4) < std::min(dist1,dist2))
	    {
	      loop_cvs[0]->reverseParameterDirection();
	      pos = pos0;
	      dist1 = dist3;
	      dist2 = dist4;
	    }
	}
      if (std::min(dist1, dist2) > tol)
	return false;
      
      if (dist2 < dist1)
	{
	  loop_cvs[ki+1]->reverseParameterDirection();
	  pos = pos1;
	}
      else
	pos = pos2;
    }

  return true;
}

