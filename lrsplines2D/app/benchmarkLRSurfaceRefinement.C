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

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/Direction2D.h"
#include "GoTools/lrsplines2D/LRBenchmarkUtils.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
//#include "GoTools/lrsplines2D/Element.h" // sbr/trd version (based on trondheim version).
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/utils/config.h"



#include <iostream>
#include <assert.h>
#include <fstream>
#include <string>
#include "string.h"


//using LRSpline;

using namespace Go;
using std::vector;


bool my_equal_function (double i, double j)
{
    double tol = 1e-10;//14;
    return (fabs(i - j) < tol);
}

// Assuming the domain is the same.
double maxDist(const Go::ParamSurface* param_sf,
	       const Go::LRSplineSurface& lr_spline_sf,
	       int nmb_samples_u, int nmb_samples_v);


int main(int argc, char *argv[])
{
  if (argc != 4)
  {
      std::cout << "Usage: spline_sf.g2 min_num_coefs_each_dir num_iter" << std::endl;
      return -1;
  }

  char* filein_char = argv[1];
  std::ifstream filein(argv[1]); // Input lr spline.
  // If input surface is not fine enough we refine.
  int min_num_each_dir = atoi(argv[2]);
  int num_iter = atoi(argv[3]);
  if (num_iter != 1)
      puts("Not using num_iter yet.");

  const int num_samples_u = 101; // Rather random.
  const int num_samples_v = 143;

  shared_ptr<Go::SplineSurface> spline_sf;
//  shared_ptr<Go::SplineSurface> lr_spline_sf_go;
  shared_ptr<LRSplineSurface> lr_spline_sf, lr_spline_sf_single_refs;

  int order_u, order_v, num_coefs_u, num_coefs_v, dim, num_bases=-1;
  if (strstr(filein_char, ".g2"))
  {
      spline_sf = shared_ptr<Go::SplineSurface>(new Go::SplineSurface());
      puts("Reading the Go::SplineSurface.");
      Go::ObjectHeader header;
      // header.read(filein);
      filein >> header;
      // We check that a lr spline was read.
      if (header.classType() != Class_SplineSurface)
      {
	  puts("The input g2-file was not of type Class_SplineSurface");
	  return -1;
      }
      // We read using the GO:LRSplineSurface routine.
      filein >> *spline_sf;
      assert(!spline_sf->rational());

      order_u = spline_sf->order_u();
      order_v = spline_sf->order_v();
      num_coefs_u = spline_sf->numCoefs_u();
      num_coefs_v = spline_sf->numCoefs_v();
      dim = spline_sf->dimension();

      // We may refine the input surface.
      if ((num_coefs_u < min_num_each_dir) || (num_coefs_v < min_num_each_dir))
      {
	  // We want to insert knots so that they are as evenly as
	  // possible distributed. Or not far from it at least.
	  MESSAGE("Refining the surface!");
	  int num_new_knots_u = min_num_each_dir - num_coefs_u;
	  if (num_new_knots_u > 0)
	  {
	      BsplineBasis bas_u = spline_sf->basis_u();
	      GeometryTools::insertKnotsEvenly(bas_u, num_new_knots_u);
	      std::cout << "spline_sf bas_u size: " << spline_sf->basis_u().numCoefs() <<
		  ", bas_u size: " << bas_u.numCoefs() << std::endl;
	      vector<double> new_knots_u;
	      set_symmetric_difference(spline_sf->basis_u().begin(), spline_sf->basis_u().end(),
				       bas_u.begin(), bas_u.end(),
				       back_inserter(new_knots_u));
	      spline_sf->insertKnot_u(new_knots_u);
	  }

	  int num_new_knots_v = min_num_each_dir - num_coefs_v;
	  if (num_new_knots_v > 0)
	  {
	      BsplineBasis bas_v = spline_sf->basis_v();
	      GeometryTools::insertKnotsEvenly(bas_v, num_new_knots_v);
	      vector<double> new_knots_v;
	      set_symmetric_difference(spline_sf->basis_v().begin(), spline_sf->basis_v().end(),
				       bas_v.begin(), bas_v.end(),
				       back_inserter(new_knots_v));
	      spline_sf->insertKnot_v(new_knots_v);
	  }
      }
      std::ofstream fileout_spline("tmp/spline_sf_ref.g2");
      spline_sf->writeStandardHeader(fileout_spline);
      spline_sf->write(fileout_spline);
  }
  else
  {
      MESSAGE("Input was not a g2-file!");
      return -1;
  }

  // order_u = spline_sf->order_u();
  // order_v = spline_sf->order_v();
  num_coefs_u = spline_sf->numCoefs_u();
  num_coefs_v = spline_sf->numCoefs_v();
  // We must count the number of segments. We do this the easiest way
  // by looking at the number of different knots in each direction,
  // (and subtract one).
  vector<double> all_knots_u(spline_sf->basis_u().begin(), spline_sf->basis_u().end());
  vector<double> all_knots_v(spline_sf->basis_v().begin(), spline_sf->basis_v().end());
  vector<double>::iterator last_uniq_u = unique(all_knots_u.begin(), all_knots_u.end(), my_equal_function);
  vector<double>::iterator last_uniq_v = unique(all_knots_v.begin(), all_knots_v.end(), my_equal_function);
  vector<double> unique_knots_u(all_knots_u.begin(), last_uniq_u);
  vector<double> unique_knots_v(all_knots_v.begin(), last_uniq_v);

  int num_segs_u = unique_knots_u.size() - 1;
  int num_segs_v = unique_knots_v.size() - 1;
  std::cout << "num_segs_u: " << num_segs_u << ", num_segs_v: " << num_segs_v << std::endl;

  // We create the refinement vectors.
  // First segments with a constant u parameter.
  vector<LRSplineSurface::Refinement2D> refs_u, refs_v;
  int mult = 1;
  Direction2D dir = XFIXED; // XDIRFIXED, refining in the x-dir (u).
  for (size_t ki = 0; ki < num_segs_u - 1; ++ki)
  {
      // We bisect the interval.
      double wgt = 0.4; // Value in the range [0.0, 1.0].
      double ref_par = wgt*unique_knots_u[ki] + (1.0 - wgt)*unique_knots_u[ki+1];

      // We insert line segments which cover order_v - 1 intervals.
      for (size_t kj = order_v - 1; kj < num_segs_v; kj += order_v)
      {
	  // To avoid a too regular pattern we alter start
	  // values. Maybe we should use a random function ...
	  double tmin = (ki%2 == 1) ? unique_knots_v[kj - order_v + 1] : unique_knots_v[kj - order_v + 2];
	  double tmax = (ki%2 == 1) ? unique_knots_v[kj] : unique_knots_v[kj + 1];

	  LRSplineSurface::Refinement2D ref;
	  ref.kval = ref_par;
	  ref.start = tmin;
	  ref.end = tmax;
	  ref.d = dir;
	  ref.multiplicity = mult;
	  refs_u.push_back(ref); // Refinements in the u-dir.
      }
  }

  dir = YFIXED; // YDIRFIXED, refining in the y-dir (v).
  for (size_t ki = 0; ki < num_segs_v - 1; ++ki)
  {
      // We bisect the interval.
      double wgt = 0.4;
      double ref_par = wgt*unique_knots_v[ki] + (1.0 - wgt)*unique_knots_v[ki+1];

      // We insert line segments which cover order_u - 1 intervals.
      for (size_t kj = order_u - 1; kj < num_segs_u; kj += order_u)
      {
	  // To avoid a too regular pattern we alter start
	  // values. Maybe we should use a random function ...
	  double tmin = (ki%2 == 1) ? unique_knots_u[kj - order_u + 1] : unique_knots_u[kj - order_u + 2];
	  double tmax = (ki%2 == 1) ? unique_knots_u[kj] : unique_knots_u[kj + 1];

	  LRSplineSurface::Refinement2D ref;
	  ref.kval = ref_par;
	  ref.start = tmin;
	  ref.end = tmax;
	  ref.d = dir;
	  ref.multiplicity = mult;
	  refs_v.push_back(ref); // Refinements in the v-dir.
      }
  }

  vector<LRSplineSurface::Refinement2D> all_refs(refs_u.begin(), refs_u.end());
  all_refs.insert(all_refs.end(), refs_v.begin(), refs_v.end());
  // if (all_refs.size() > 0)
  //     all_refs.pop_back();
  // if (all_refs.size() > 0)
  //     all_refs.pop_back();
  // if (all_refs.size() > 0)
  //     all_refs.pop_back();
  std::cout << "num_refs_to_insert: " << all_refs.size() << std::endl;

//  all_refs.erase(all_refs.begin() + 1, all_refs.begin() + 3);
//  all_refs.erase(all_refs.begin() + 4, all_refs.end());
//  all_refs.erase(all_refs.begin(), all_refs.begin() + 2);
  std::cout << "num_refs_to_insert: " << all_refs.size() << std::endl;


  // We then create and refine the LRSplineSurface.
  lr_spline_sf = shared_ptr<LRSplineSurface>
      // (new LRSpline<vector<double>::const_iterator, vector<double>::const_iterator> >
      (new LRSplineSurface(order_u - 1, order_v - 1,
			   num_coefs_u, num_coefs_v,
			   dim,
			   spline_sf->basis_u().begin(),
			   spline_sf->basis_v().begin(),
			   spline_sf->coefs_begin()));

  lr_spline_sf_single_refs = shared_ptr<LRSplineSurface>
      // (new LRSpline<vector<double>::const_iterator, vector<double>::const_iterator> >
      (new LRSplineSurface(order_u - 1, order_v - 1,
			   num_coefs_u, num_coefs_v,
			   dim,
			   spline_sf->basis_u().begin(),
			   spline_sf->basis_v().begin(),
			   spline_sf->coefs_begin()));

  bool refine_multi = true;
  if (refine_multi)
  {
      double time_lrspline_lr_ref = benchmarkSfRefinement(*lr_spline_sf, all_refs);
      std::cout << "Time lr refinement: " << time_lrspline_lr_ref << std::endl;
      // We write to screen the number of basis functions for both
      // versions.
      int num_elem = lr_spline_sf->numElements();
      int num_basis_funcs = lr_spline_sf->numBasisFunctions();
      std::cout << "num_elem: " << num_elem << ", num_basis_funcs: " << num_basis_funcs << std::endl;

      double max_dist = maxDist(spline_sf.get(), *lr_spline_sf, num_samples_u, num_samples_v);
      std::cout << "Max dist between input and (multi) refined surface: " << max_dist << std::endl;

      std::ofstream fileout("tmp/ref_lr_multi.g2");
      lr_spline_sf->writeStandardHeader(fileout);
      lr_spline_sf->write(fileout);

      // MESSAGE("Missing writing refined surface grid to file!");
      std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");
//  writePostscriptMesh(*lrsf);
      writePostscriptMesh(*lr_spline_sf, lrsf_grid_ps);

  }

  bool refine_single = true;
  if (refine_single)
  {
      double time_lrspline_lr_ref_single = benchmarkSfRefinement(*lr_spline_sf_single_refs, all_refs, true);
      std::cout << "Time lr refinement single refs: " << time_lrspline_lr_ref_single << std::endl;
      int num_elem_single = lr_spline_sf_single_refs->numElements();
      int num_basis_funcs_single = lr_spline_sf_single_refs->numBasisFunctions();
      std::cout << "num_elem_single_refs: " << num_elem_single <<
	  ", num_basis_funcs_single_refs: " << num_basis_funcs_single << std::endl;

#if 0//ndef NDEBUG
  vector<LRBSpline2D*> bas_funcs;
  for (auto iter = lr_spline_sf_single_refs->basisFunctionsBegin(); iter != lr_spline_sf_single_refs->basisFunctionsEnd(); ++iter)
    {
      bas_funcs.push_back((*iter).second.get());
    }
  vector<Element2D*> elems;
  for (auto iter = lr_spline_sf_single_refs->elementsBegin(); iter != lr_spline_sf_single_refs->elementsEnd(); ++iter)
  {
      elems.push_back(((*iter).second.get()));
  }
  puts("Remove when done debugging!");
#endif

      double max_dist_post_ref_multi_ref = maxDist(spline_sf.get(), *lr_spline_sf_single_refs, num_samples_u, num_samples_v);
      std::cout << "Max dist input and (single) ref surface: " << max_dist_post_ref_multi_ref << std::endl;

      std::ofstream fileout2("tmp/ref_lr_single.g2");
      lr_spline_sf_single_refs->writeStandardHeader(fileout2);
      lr_spline_sf_single_refs->write(fileout2);
  }
}



double maxDist(const Go::ParamSurface* param_sf,
	       const Go::LRSplineSurface& lr_spline_sf,
	       int num_samples_u, int num_samples_v)
{
    // Assuming the domain is the same.
    double umin = lr_spline_sf.startparam_u();
    double umax = lr_spline_sf.endparam_u();
    double vmin = lr_spline_sf.startparam_v();
    double vmax = lr_spline_sf.endparam_v();
    double ustep = (umax - umin)/((double)num_samples_u - 1);
    double vstep = (vmax - vmin)/((double)num_samples_v - 1);
    Go::Point go_pt(3), lr_pt(3);
    double max_dist = -1.0;
// #ifndef NDEBUG
    double max_dist_u = 0.0;
    double max_dist_v = 0.0;
    Go::Point max_go_pt(3), max_lr_pt(3);
// #endif
    for (int kj = 0; kj < num_samples_v; ++kj)
    {
	double vpar = vmin + kj*vstep;
	for (int ki = 0; ki < num_samples_u; ++ki)
	{
	    double upar = umin + ki*ustep;
	    param_sf->point(go_pt, upar, vpar);
	    lr_spline_sf.point(lr_pt, upar, vpar);
	    double dist = go_pt.dist(lr_pt);
	    if (dist > max_dist)
	    {
		max_dist = dist;
// #ifndef NDEBUG
		max_dist_u = upar;
		max_dist_v = vpar;
		max_go_pt = go_pt;
		max_lr_pt = lr_pt;
// #endif
	    }
	}
    }    

// #ifndef NDEBUG
    std::cout << "max_dist: u = " << max_dist_u << ", v = " << max_dist_v << std::endl;
    std::cout << "max_go_pt = " << max_go_pt << std::endl;
    std::cout << "max_lr_pt = " << max_lr_pt << std::endl;
// #endif

    return max_dist;


}
