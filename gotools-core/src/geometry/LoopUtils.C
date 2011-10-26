//===========================================================================
//                                                                           
// File: LoopUtils.C                                                         
//                                                                           
// Created: Tue Jun  3 15:09:13 2003                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: LoopUtils.C,v 1.18 2009-03-04 15:42:10 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/LoopUtils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/geometry/SplineDebugUtils.h"

#include <fstream>


using namespace Go;
using std::shared_ptr;
using std::dynamic_pointer_cast;
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
representAsSurfaceCurves(std::vector< std::shared_ptr<ParamCurve> >& curves,
			 std::shared_ptr<BoundedSurface> surf,
			 std::vector<std::shared_ptr<CurveOnSurface> >& cvs_on_sf)
//===========================================================================
{
    cvs_on_sf.resize(curves.size());

    int sfdim = surf->dimension();
    for (size_t ki=0; ki<curves.size(); ++ki)
    {
	shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curves[ki]);
	if (sfcv.get() != 0)
	{
	    ALWAYS_ERROR_IF(sfcv->underlyingSurface().get() != 
			    surf->underlyingSurface().get(),
			    "Inconsistent surface pointers");
	    cvs_on_sf[ki] = sfcv;
	}
	else
	{
	    int cvdim = curves[ki]->dimension();
	    if (cvdim == sfdim && sfdim > 2)
	    {
		sfcv 
		    = shared_ptr<CurveOnSurface>(new CurveOnSurface(surf, 
								    curves[ki],
								    false));
	    }
	    else if (cvdim == 2)
	    {
		sfcv
		    = shared_ptr<CurveOnSurface>(new CurveOnSurface(surf, 
								    curves[ki],
								    true));
	    }
	    else
		cvs_on_sf[ki] = sfcv;
	}
    }

}

//===========================================================================
bool
LoopUtils::loopIsCCW(const vector<shared_ptr<SplineCurve> >& simple_par_loop, 
		     double int_tol)
//===========================================================================
{
    ALWAYS_ERROR_IF(simple_par_loop.size() == 0,
		"Empty input vector!");

    int ki;
    int idx = (int)simple_par_loop.size()/2;
    ALWAYS_ERROR_IF(simple_par_loop[idx]->dimension() != 2,
		"Input loop must be 2-dimensional.");

    // We choose the mid parameter value on the middle curve in the loop.
    double tpar =
      0.5*(simple_par_loop[idx]->startparam() + simple_par_loop[idx]->endparam());

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

    double length = 1.0 + (box1.low()).dist(box1.high());
    normal.normalize();
    Point end_pt = pnt[0] + length*normal;
    SplineCurve normal_curve = SplineCurve(pnt[0], end_pt);

    // We then check for intersections between normal_curve and simple_par_loop,
    // not counting start point of normal_curve.    
    vector<double> params_interval;
    vector<shared_ptr<ParamCurve> > par_loop;
    for (size_t ki=0; ki<simple_par_loop.size(); ++ki)
      par_loop.push_back(simple_par_loop[ki]);
    shared_ptr<CurveLoop> loop = 
      shared_ptr<CurveLoop>(new CurveLoop(par_loop, int_tol));
    CurveBoundedDomain loop_dom(loop);

    loop_dom.findPcurveInsideSegments(normal_curve, int_tol, params_interval);

    int nmbpoint = (int)params_interval.size();

  // Remove first intersection. If odd # intersections, loop is CCW.
  if ((params_interval.size() > 0) && 
      (params_interval[0] == normal_curve.startparam())) 
    nmbpoint--;

  return ((nmbpoint % 2 == 1) ? true : false);
}


//===========================================================================
bool LoopUtils::paramIsCCW(const vector< shared_ptr<CurveOnSurface> >& loop, 
			   double int_tol)
//===========================================================================
{
    if (loop.empty())
	return true;
    vector< shared_ptr<SplineCurve> > sc;
    for (size_t i = 0; i < loop.size(); ++i) {
      loop[i]->ensureParCrvExistence(int_tol);
	ParamCurve* pcptr = loop[i]->parameterCurve().get();
	ALWAYS_ERROR_IF(pcptr == 0, "Parameter curve not found");
	shared_ptr<ParamCurve> pc(loop[i]->parameterCurve());
	shared_ptr<SplineCurve>
	    spc(std::dynamic_pointer_cast<SplineCurve, ParamCurve>(pc));
	sc.push_back(spc);
    }

    return loopIsCCW(sc, int_tol);
}


//===========================================================================
bool LoopUtils::loopIsCCW(const CurveLoop& loop, double int_tol)
//===========================================================================
{
    // We extract the 2D parameter curves.  Well, come to think of it,
    // what really matters is not that the curve is a 2D-spline, but
    // 2D ...
    vector< shared_ptr<SplineCurve> > sc;
    for (int ki = 0; ki < loop.size(); ++ki) {
	shared_ptr<SplineCurve> pcv;
	if (loop[ki]->instanceType() == Class_SplineCurve)
	    pcv = dynamic_pointer_cast<SplineCurve, ParamCurve>(loop[ki]);
	else if (loop[ki]->instanceType() == Class_CurveOnSurface) {
	    shared_ptr<CurveOnSurface> cv_on_sf =
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(loop[ki]);
	    if (cv_on_sf->parameterCurve().get() == NULL)
		THROW("Method requires parameter curve!");
	    if (cv_on_sf->parameterCurve()->instanceType() == Class_SplineCurve)
		pcv = dynamic_pointer_cast<SplineCurve, ParamCurve>
		    (cv_on_sf->parameterCurve());
	    else {
		pcv = shared_ptr<SplineCurve>
		    (cv_on_sf->parameterCurve()->geometryCurve());
		if (pcv.get() == NULL)
		    THROW("Unexpected incident.");
	    }
	} else
	    THROW("Unexpected incident.");

	sc.push_back(pcv);
    }

    return loopIsCCW(sc, int_tol);
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

   // As we're creating bd domains we must make sure loop defines bd area!
   // We're not interested in the direction anyway.
   vector<shared_ptr<CurveOnSurface> > first_ccw_loop = first_loop;
   if (!paramIsCCW(first_ccw_loop, int_tol)) {
       for (ki = 0; ki < int(first_ccw_loop.size()); ++ki) {
	   first_ccw_loop[ki] =
	       shared_ptr<CurveOnSurface>(first_loop[ki]->clone());
	   first_ccw_loop[ki]->reverseParameterDirection();
       }
       std::reverse(first_ccw_loop.begin(), first_ccw_loop.end());
   }
   vector<shared_ptr<CurveOnSurface> > second_ccw_loop = second_loop;
   if (!paramIsCCW(second_ccw_loop, int_tol)) {
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
   pcv = std::dynamic_pointer_cast<SplineCurve, ParamCurve>
       (first_ccw_loop[0]->parameterCurve());


   ALWAYS_ERROR_IF(pcv.get() == 0, "Unable to convert curve to SplineCurve");
   vector<double>::const_iterator iter = pcv->basis().begin();
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
	   bd_domain1.clipWithDomain(dir, par, int_tol,
				     under_sf, trim_pieces[0]);
	   bd_domain2.clipWithDomain(dir, par, int_tol,
				     under_sf, trim_pieces[1]);
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

	   if (((end_other_par[0] < other_par) &&
		(end_other_par[0] > end_other_par[1])) ||
	       ((end_other_par[0] > other_par) &&
		(end_other_par[1] > end_other_par[0]))) {
	       return true;
	   } else {
	       return false;
	   }
       }
   }
}
