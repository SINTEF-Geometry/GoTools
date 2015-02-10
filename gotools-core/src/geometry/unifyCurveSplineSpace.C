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

#include "GoTools/geometry/GeometryTools.h"

//#ifdef __BORLANDC__
#include <iterator> // For back_inserter.  Should be required by VC++ and GCC as well.
//#endif
//***************************************************************************
//
// Implementation file of the free function unifyCurveSplineSpace defined in
// GeometryTools.h/
//
//***************************************************************************

using std::vector;
using std::setprecision;
using std::endl;
using std::pair;
using std::make_pair;
using std::back_inserter;

namespace Go
{

//===========================================================================
void GeometryTools::unifyCurveSplineSpace(vector<shared_ptr<SplineCurve> >& curves,
			   double tol)
//===========================================================================
//********************************************************************
// Make sure that a set of curves live on the same knot vector
//********************************************************************
{
  int ki;
  int num_curves = (int)curves.size();
  // Do not allow input curves to have inner knots with pos difference less than tol.
  for (ki = 0; ki < num_curves; ++ki)
    if (curves[ki].get() != 0) {
      vector<double>::iterator iter = curves[ki]->basis().begin();
      vector<double>::iterator next_iter = iter + 1;
      while (next_iter != curves[ki]->basis().end()) {
	if ((*iter != *next_iter) && (*next_iter - *iter < tol))
	  THROW("Inner knot difference too small!");
	iter = next_iter;
	next_iter = iter + 1;
      }
    }

  // We allow members of curves to be NULL pointers.
  vector<shared_ptr<SplineCurve> > copy_curves;
  vector<bool> curve_exists(num_curves);
  for (ki = 0; ki < num_curves; ++ki)
    if (curves[ki].get() != 0) {
      copy_curves.push_back(curves[ki]);
      curve_exists[ki] = true;
    } else
      curve_exists[ki] = false;

  if (copy_curves.size() <= 1)
    return;  // Nothing to do

  // Raise the order of each curve to the maximum order
  int nmb_crv = (int)copy_curves.size();
  int order_max = copy_curves[0]->order();
  int order;
  for (ki=1; ki<nmb_crv; ki++)
    {
      order = copy_curves[ki]->order();
      order_max = std::max(order_max, order);
    }

  for (ki=0; ki<nmb_crv; ki++)
    {
      // Also make sure that the copy_curves have multiple knots in the ends
      copy_curves[ki]->makeKnotStartRegular();
      copy_curves[ki]->makeKnotEndRegular();
      order = copy_curves[ki]->order();
      if (order < order_max)
	copy_curves[ki]->raiseOrder(order_max - order);
    }

  // Unify knot vectors
  // First find the union of all knot vector
  vector<double> union_knots;
  vector<BsplineBasis> bbasis(nmb_crv);
  for (ki=0; ki<nmb_crv; ++ki)
    bbasis[ki] = copy_curves[ki]->basis();

  GeometryTools::makeUnionKnots(bbasis, tol, union_knots);

  // Insert missing knots into the copy_curves
  int num_union = (int)union_knots.size() - order_max;
  for (ki=0; ki<nmb_crv; ki++)
    {
      vector<double> new_knots;
      double knot;
      int kj, kh;
      for (kj=order_max, kh=order_max; 
	   kj<copy_curves[ki]->numCoefs() || kh<num_union;)
	{
	  knot = copy_curves[ki]->basis().begin()[kj];
	  if (fabs(knot - union_knots[kh]) < tol)
	    {
	      kj++;
	      kh++;
	    }
	  else if (union_knots[kh] < knot)
	    {
	      new_knots.push_back(union_knots[kh]);
	      kh++;
	    }
	  else
	    kj++;
	}

      // Insert new knots
      copy_curves[ki]->insertKnot(new_knots);
      // We test whether we must make a new spline curve in order to
      // alter tol-equal knots.
      vector<double> difference;
      std::set_difference(union_knots.begin(), union_knots.end(),
			  copy_curves[ki]->basis().begin(),
			  copy_curves[ki]->basis().end(),
			  std::back_inserter(difference));
      if (difference.size() != 0)
	{
	  SplineCurve tmp(copy_curves[ki]->numCoefs(),
			  copy_curves[ki]->order(),
			  union_knots.begin(),
			  copy_curves[ki]->coefs_begin(),
			  copy_curves[ki]->dimension());
	  copy_curves[ki]->swap(tmp);
	}
    }

  // We make sure the original curve is altered.
  int kj = 0;
  for (ki = 0; ki < num_curves; ++ki)
    if (curve_exists[ki]) {
      curves[ki] = copy_curves[kj];
      ++kj;
    }
}

  void GeometryTools::makeUnionKnots(vector<BsplineBasis>& bbasis,
		      double tol, vector<double>& union_knots)
  {
    union_knots.clear();

  vector< std::vector<double>::const_iterator > c_ptr;
  vector< std::vector<double>::const_iterator > c_end;
  int nmb_crv = (int)bbasis.size();
  c_ptr.resize(nmb_crv);
  c_end.resize(nmb_crv);
  int ki;
  for (ki=0; ki<nmb_crv; ki++)
    {
      c_ptr[ki] = bbasis[ki].begin();
      c_end[ki] = bbasis[ki].end();
    }

  double min_knot, knot;
  int max_mult, mult;
  while (true)
    {
      // More knots?
      for (ki=0; ki<nmb_crv; ki++)
	if (c_ptr[ki] < c_end[ki])
	  break;

      if (ki == nmb_crv)
	break;  // All knots collected

      min_knot = c_ptr[ki][0];
      for (ki=1; ki<nmb_crv; ki++)
      {
	  if (c_ptr[ki] != c_end[ki]) {
	      knot = c_ptr[ki][0];
	      min_knot = std::min(min_knot, knot);
	  }
      }
       
      mult = 0;
      max_mult = 1;
      for (ki=0; ki<nmb_crv; ki++)
	{
	  knot = c_ptr[ki][0];
	  if (knot < min_knot + tol)
	    for (mult=0; c_ptr[ki]<c_end[ki] && c_ptr[ki][0]==knot; 
		 mult++, c_ptr[ki]++);
	  max_mult = std::max(max_mult, mult);
	}

      for (ki=0; ki<max_mult; ki++)
	union_knots.push_back(min_knot);
    }
  }



void GeometryTools::makeUnionKnots(vector<vector<double> >& knots,
				   double tol, vector<double>& union_knots)
  {
    union_knots.clear();

  vector< std::vector<double>::const_iterator > c_ptr;
  vector< std::vector<double>::const_iterator > c_end;
  int nmb_vec = (int)knots.size();
  c_ptr.resize(nmb_vec);
  c_end.resize(nmb_vec);
  int ki;
  for (ki=0; ki<nmb_vec; ki++)
    {
      c_ptr[ki] = knots[ki].begin();
      c_end[ki] = knots[ki].end();
    }

  double min_knot, knot;
  int max_mult, mult;
  while (true)
    {
      // More knots?
      for (ki=0; ki<nmb_vec; ki++)
	if (c_ptr[ki] < c_end[ki])
	  break;

      if (ki == nmb_vec)
	break;  // All knots collected

      min_knot = c_ptr[ki][0];
      for (ki=1; ki<nmb_vec; ki++)
      {
	  if (c_ptr[ki] != c_end[ki]) {
	      knot = c_ptr[ki][0];
	      min_knot = std::min(min_knot, knot);
	  }
      }
       
      mult = 0;
      max_mult = 1;
      for (ki=0; ki<nmb_vec; ki++)
	{
	  knot = c_ptr[ki][0];
	  if (knot < min_knot + tol)
	    for (mult=0; c_ptr[ki]<c_end[ki] && c_ptr[ki][0]==knot; 
		 mult++, c_ptr[ki]++);
	  max_mult = std::max(max_mult, mult);
	}

      for (ki=0; ki<max_mult; ki++)
	union_knots.push_back(min_knot);
    }
  }


}
