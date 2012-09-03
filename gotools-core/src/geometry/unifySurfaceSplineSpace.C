//===========================================================================
//                                                                           
// File: unifySurfaceSplineSpace.C                                                        
//                                                                           
// Created: Thu Apr  5 17:32:06 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: unifySurfaceSplineSpace.C,v 1.15 2009-05-13 07:30:53 vsk Exp $
//                                                                           
// Description: Given input of surfaces, function puts them on a joint basis.
//              This is done by setting tol-equal-knots equal, and then insert
//              missing knots.
//              Assumes surfaces are defined over the same parameter domain.
//              Function defined in GeometryTools.h.
//===========================================================================


#include "GoTools/geometry/GeometryTools.h"
#include <algorithm>
//#ifdef __BORLANDC__
#include <iterator> // For back_inserter.  Should be required by VC++ and GCC as well.
//#endif

using std::vector;
using std::setprecision;
using std::endl;
using std::pair;
using std::make_pair;
using std::back_inserter;

namespace Go{

  void GeometryTools::unifySurfaceSplineSpace(vector<shared_ptr<SplineSurface> >& surfaces,
			       double tol, int dir) // 0 means both, 1 is u, 2 is v
    {

      if (surfaces.size() <= 1)
	return;

      int nmb_srfs = (int)surfaces.size();
      double startparam_u = surfaces[0]->startparam_u();
      double endparam_u = surfaces[0]->endparam_u();
      double startparam_v = surfaces[0]->startparam_v();
      double endparam_v = surfaces[0]->endparam_v();
      int i, j, h, g;
      // We first locate a joint parameter domain.
      for (i = 1; i < nmb_srfs; ++i) {
	startparam_u = std::max(startparam_u, surfaces[i]->startparam_u());
	endparam_u = std::min(endparam_u, surfaces[i]->endparam_u());
	startparam_v = std::max(startparam_v, surfaces[i]->startparam_v());
	endparam_v = std::min(endparam_v, surfaces[i]->endparam_v());
      }
      // We make sure the surfaces are k-regular, over tol-equal param domain.
      for (i = 0; i < nmb_srfs; ++i) {
	if (((dir != 2) &&
	   ((fabs(startparam_u - surfaces[i]->startparam_u()) > tol) ||
	     (fabs(endparam_u - surfaces[i]->endparam_u()) > tol))) ||
	    ((dir != 1) &&
	     ((fabs(startparam_v - surfaces[i]->startparam_v()) > tol) ||
	      (fabs(endparam_v - surfaces[i]->endparam_v()) > tol)))) {
	    MESSAGE("Surfaces defined over different parameter domains, reparametrizing!");
	    double umin = (dir != 2) ? startparam_u : surfaces[i]->startparam_u();
	    double umax = (dir != 2) ? endparam_u : surfaces[i]->endparam_u();
	    double vmin = (dir != 1) ? startparam_v : surfaces[i]->startparam_v();
	    double vmax = (dir != 1) ? endparam_v : surfaces[i]->endparam_v();
	    surfaces[i]->setParameterDomain(umin, umax, vmin, vmax);
	}

// #ifdef _MSC_VER
// 	surfaces[i] = shared_ptr<SplineSurface>
// 	    (dynamic_cast<SplineSurface*>(surfaces[i]->subSurface(surfaces[i]->startparam_u(),
// 								    surfaces[i]->startparam_v(),
// 								    surfaces[i]->endparam_u(),
// 								    surfaces[i]->endparam_v())));
// #else
	surfaces[i]->makeSurfaceKRegular();
// #endif
      }

      int order_u_max = surfaces[0]->order_u();
      int order_v_max = surfaces[0]->order_v();
      // Raise the order of each surface to the maximum order
      for (i = 1; i < nmb_srfs; ++i) {
	  order_u_max = std::max(order_u_max, surfaces[i]->order_u());
	  order_v_max = std::max(order_v_max, surfaces[i]->order_v());
      }
      for (i = 0; i < nmb_srfs; ++i)
	if ((order_u_max > surfaces[i]->order_u()) ||
	    (order_v_max > surfaces[i]->order_v()))
	  surfaces[i]->raiseOrder(order_u_max - surfaces[i]->order_u(),
				  order_v_max - surfaces[i]->order_v());

      // This last section is cut and paste from unifyCurveSplineSpace.C...
      // Unify knot vectors
      // First find the union of all knot vectors
      vector<double> union_knots_u, union_knots_v;
      for (g = 0; g < 2; ++g) {
	  if (((g == 0) && (dir == 2)) || ((g == 1) && (dir == 1)))
	  continue; // We will not unify in that direction.

	vector< std::vector<double>::const_iterator > c_ptr;
	vector< std::vector<double>::const_iterator > c_end;
	c_ptr.resize(nmb_srfs);
	c_end.resize(nmb_srfs);
	for (i = 0; i < nmb_srfs; ++i)
	  {
	    c_ptr[i] = (g == 0) ? surfaces[i]->basis_u().begin() :
	      surfaces[i]->basis_v().begin();
	    c_end[i] = (g == 0) ? surfaces[i]->basis_u().end() :
	      surfaces[i]->basis_v().end();
	  }

	vector<double> union_knots;
	double min_knot, knot;
	int max_mult, mult;
	while (true)
	  {
	    // More knots?
	    for (i = 0; i < nmb_srfs; ++i)
	      if (c_ptr[i] < c_end[i])
		{
		  min_knot = c_ptr[i][0];
		  break;
		}

	    if (i == nmb_srfs)
	      break;  // All knots collected

	    for (i = 1; i < nmb_srfs; ++i)
	      {
		if (c_ptr[i] < c_end[i])
		  {
		    knot = c_ptr[i][0];
		    min_knot = std::min(min_knot, knot);
		  }
	      }
       
	    mult = 1;
	    max_mult = 1;
	    for (i = 0; i < nmb_srfs; ++i)
	      {
		if (c_ptr[i] < c_end[i])
		  {
		    knot = c_ptr[i][0];
		    if (knot < min_knot + tol)
		      for (mult=0; c_ptr[i]<c_end[i] && c_ptr[i][0]==knot; 
			   mult++, c_ptr[i]++);
		    max_mult = std::max(max_mult, mult);
		  }
	      }

	    for (i = 0; i < max_mult; ++i)
	      union_knots.push_back(min_knot);
	  }

	// Extract the new knot vectors.
	if (g == 0)
	  union_knots_u = union_knots;
	else union_knots_v = union_knots;
      }

      // Insert missing knots into the surfaces
      for (i = 0; i < nmb_srfs; ++i) {
	  for (g = 0; g < 2; ++g) {
	      if (((g == 0) && (dir == 2)) || ((g == 1) && (dir == 1)))
		  continue; // We will not unify in that direction.

	      vector<double> union_knots = (g == 0) ? union_knots_u : union_knots_v;
	      double knot;
	      int order_max = (g == 0) ? order_u_max : order_v_max;
	      int num_union = (int)union_knots.size() - order_max;
	      int num_coefs = (g == 0) ? surfaces[i]->numCoefs_u() :
		  surfaces[i]->numCoefs_v();
	      vector<double> new_knots;
	      for (j = order_max, h = order_max; 
		   j < num_coefs || h < num_union;)
		  {
		      knot = (g == 0) ? surfaces[i]->basis_u().begin()[j] :
			  surfaces[i]->basis_v().begin()[j];
		      if (fabs(knot - union_knots[h]) < tol)
			  {
			      ++j;
			      ++h;
			  }
		      else if (union_knots[h] < knot)
			  {
			      new_knots.push_back(union_knots[h]);
			      ++h;
			  }
		      else
			  ++j;
		  }

	      if (g == 0)
		  surfaces[i]->insertKnot_u(new_knots);
	      else surfaces[i]->insertKnot_v(new_knots);
	  }

	  // We test whether we must make a new spline surface in order to
	  // alter tol-equal knots.
	  vector<double> difference_u, difference_v;
	  if (dir == 1)
	      union_knots_v.assign(surfaces[i]->basis_v().begin(), surfaces[i]->basis_v().end());
	  else if (dir == 2)
	      union_knots_u.assign(surfaces[i]->basis_u().begin(), surfaces[i]->basis_u().end());
	  std::set_difference(union_knots_u.begin(), union_knots_u.end(),
			      surfaces[i]->basis_u().begin(),
			      surfaces[i]->basis_u().end(),
			      std::back_inserter(difference_u));
	  std::set_difference(union_knots_v.begin(), union_knots_v.end(),
			      surfaces[i]->basis_v().begin(),
			      surfaces[i]->basis_v().end(),
			      std::back_inserter(difference_v));
	  if ((difference_u.size() != 0) || (difference_v.size() != 0))
	    {
	      SplineSurface tmp(surfaces[i]->numCoefs_u(),
				surfaces[i]->numCoefs_v(),
				surfaces[i]->order_u(),
				surfaces[i]->order_v(),
				union_knots_u.begin(),
				union_knots_v.begin(),
				surfaces[i]->coefs_begin(),
				surfaces[i]->dimension());
	      surfaces[i]->swap(tmp);
	    }
      }
    }


  void GeometryTools::unifySurfaceSplineSpaceOneDir(vector<shared_ptr<SplineSurface> >& surfaces,
				     double tol, bool unify_u_dir)
  {
    int cv_dir = unify_u_dir ? 1 : 2;
    int nmb_sfs = (int)surfaces.size();

    vector<shared_ptr<SplineCurve> > curves;
    for (int i = 0; i < nmb_sfs; ++i)
      {
	SplineSurface* sf = surfaces[i].get();
	curves.push_back(representSurfaceAsCurve(*sf, cv_dir));
      }
    GeometryTools::unifyCurveSplineSpace(curves, tol);

    if (unify_u_dir)
      for (int i = 0; i < nmb_sfs; ++i)
	{
	  SplineCurve* cv = curves[i].get();
	  surfaces[i] = representCurveAsSurface(*cv, cv_dir, surfaces[i]->basis_v(), surfaces[i]->rational());
	}
    else
      for (int i = 0; i < nmb_sfs; ++i)
	{
	  SplineCurve* cv = curves[i].get();
	  surfaces[i] = representCurveAsSurface(*cv, cv_dir, surfaces[i]->basis_u(), surfaces[i]->rational());
	}
  }



} // namespace Go;
