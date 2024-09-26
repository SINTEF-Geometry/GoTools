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

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/Direction3D.h"
#include "GoTools/lrsplines3D/LRBenchmarkUtils3D.h"
#include "GoTools/lrsplines3D/LRSplinePlotUtils3D.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariate/SplineVolume.h"
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
double maxDist(const Go::ParamVolume* param_vol,
	       const Go::LRSplineVolume& lr_spline_vol,
	       int nmb_samples_u, int nmb_samples_v, int num_samples_w);

// Description: Read a SplineVolume, convert to LRSplineVolume, benchmark the Refinement function. We compare the
// one-at-a-time approach vs all at once.
int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      std::cout << "Usage: spline_vol.g2 min_num_coefs_each_dir num_iter" << std::endl;
      return -1;
    }

  char* filein_char = argv[1];
  std::ifstream filein(argv[1]); // Input lr spline.
  // If input surface is not fine enough we refine.
  int min_num_each_dir = atoi(argv[2]);
  int num_iter = atoi(argv[3]);
  if (num_iter != 1)
    puts("Not using num_iter yet.");

  const int num_samples_u = 11;//101; // Rather random.
  const int num_samples_v = 13;//143;
  const int num_samples_w = 17;//143;

  shared_ptr<Go::SplineVolume> spline_vol;
//  shared_ptr<Go::SplineVolume> lr_spline_vol_go;
  shared_ptr<LRSplineVolume> lr_spline_vol, lr_spline_vol_single_refs;

  int order_u, order_v, order_w, num_coefs_u, num_coefs_v, num_coefs_w, dim, num_bases=-1;
  if (strstr(filein_char, ".g2")) // We do not expect an uppercase suffix (.G2).
    {
      spline_vol = shared_ptr<Go::SplineVolume>(new Go::SplineVolume());
      puts("Reading the Go::SplineVolume.");
      Go::ObjectHeader header;
      // header.read(filein);
      filein >> header;
      // We check that a lr spline was read.
      if (header.classType() != Class_SplineVolume)
	{
	  puts("The input g2-file was not of type Class_SplineVolume");
	  return -1;
	}
      // We read using the GO:LRSplineVolume routine.
      filein >> *spline_vol;
      assert(!spline_vol->rational());

      order_u = spline_vol->order(0);
      order_v = spline_vol->order(1);
      order_w = spline_vol->order(2);
      num_coefs_u = spline_vol->numCoefs(0);
      num_coefs_v = spline_vol->numCoefs(1);
      num_coefs_w = spline_vol->numCoefs(2);
      dim = spline_vol->dimension();

      // We may refine the input volume.
      if ((num_coefs_u < min_num_each_dir) || (num_coefs_v < min_num_each_dir) || (num_coefs_w < min_num_each_dir))
	{
	  // We want to insert knots so that they are as evenly as
	  // possible distributed. Or not far from it at least.
	  MESSAGE("Refining the surface!");
	  int num_new_knots_u = min_num_each_dir - num_coefs_u;
	  if (num_new_knots_u > 0)
	    {
	      BsplineBasis bas_u = spline_vol->basis(0);
	      GeometryTools::insertKnotsEvenly(bas_u, num_new_knots_u);
	      std::cout << "spline_vol bas_u size: " << spline_vol->basis(0).numCoefs() <<
		", bas_u size: " << bas_u.numCoefs() << std::endl;
	      vector<double> new_knots_u;
	      set_symmetric_difference(spline_vol->basis(0).begin(), spline_vol->basis(0).end(),
				       bas_u.begin(), bas_u.end(),
				       back_inserter(new_knots_u));
	      spline_vol->insertKnot(0, new_knots_u);
	    }

	  int num_new_knots_v = min_num_each_dir - num_coefs_v;
	  if (num_new_knots_v > 0)
	    {
	      BsplineBasis bas_v = spline_vol->basis(1);
	      GeometryTools::insertKnotsEvenly(bas_v, num_new_knots_v);
	      vector<double> new_knots_v;
	      set_symmetric_difference(spline_vol->basis(1).begin(), spline_vol->basis(1).end(),
				       bas_v.begin(), bas_v.end(),
				       back_inserter(new_knots_v));
	      spline_vol->insertKnot(1, new_knots_v);
	    }

	  int num_new_knots_w = min_num_each_dir - num_coefs_w;
	  if (num_new_knots_w > 0)
	    {
	      BsplineBasis bas_w = spline_vol->basis(2);
	      GeometryTools::insertKnotsEvenly(bas_w, num_new_knots_w);
	      vector<double> new_knots_w;
	      set_symmetric_difference(spline_vol->basis(2).begin(), spline_vol->basis(2).end(),
				       bas_w.begin(), bas_w.end(),
				       back_inserter(new_knots_w));
	      spline_vol->insertKnot(2, new_knots_w);
	    }

	  std::ofstream fileout_spline("tmp/spline_vol_ref.g2");
	  spline_vol->writeStandardHeader(fileout_spline);
	  spline_vol->write(fileout_spline);
	}
    }
  else
    {
      MESSAGE("Input was not a g2-file!");
      return -1;
    }

  // order_u = spline_vol->order(0);
  // order_v = spline_vol->order_v();
  num_coefs_u = spline_vol->numCoefs(0);
  num_coefs_v = spline_vol->numCoefs(1);
  num_coefs_w = spline_vol->numCoefs(2);
  // We must count the number of segments. We do this the easiest way
  // by looking at the number of different knots in each direction,
  // (and subtract one).
  vector<double> all_knots_u(spline_vol->basis(0).begin(), spline_vol->basis(0).end());
  vector<double> all_knots_v(spline_vol->basis(1).begin(), spline_vol->basis(1).end());
  vector<double> all_knots_w(spline_vol->basis(2).begin(), spline_vol->basis(2).end());
  vector<double>::iterator last_uniq_u = unique(all_knots_u.begin(), all_knots_u.end(), my_equal_function);
  vector<double>::iterator last_uniq_v = unique(all_knots_v.begin(), all_knots_v.end(), my_equal_function);
  vector<double>::iterator last_uniq_w = unique(all_knots_w.begin(), all_knots_w.end(), my_equal_function);
  vector<double> unique_knots_u(all_knots_u.begin(), last_uniq_u);
  vector<double> unique_knots_v(all_knots_v.begin(), last_uniq_v);
  vector<double> unique_knots_w(all_knots_w.begin(), last_uniq_w);

  int num_segs_u = unique_knots_u.size() - 1;
  int num_segs_v = unique_knots_v.size() - 1;
  int num_segs_w = unique_knots_w.size() - 1;
  std::cout << "num_segs_u: " << num_segs_u << ", num_segs_v: " << num_segs_v <<  ", num_segs_w: " << num_segs_w << std::endl;

  // @@sbr201304 Test refinement rectangles which overlap along the
  // border. Should the second such rectangle be considered an
  // extension. I guess not.

  // We create the refinement vectors.
  // First segments with a constant u parameter.
  vector<LRSplineVolume::Refinement3D> refs_u, refs_v, refs_w;
  int mult = 1;
  Direction3D dir = XDIR; // XDIRFIXED, refining in the x-dir (u).
  if ((num_segs_v > order_v - 1) && (num_segs_w > order_w - 1))
  for (size_t ki = 0; ki < num_segs_u; ++ki)
    {
      // We bisect the interval.
      double wgt = 0.4; // Value in the range [0.0, 1.0].
      double ref_par = wgt*unique_knots_u[ki] + (1.0 - wgt)*unique_knots_u[ki+1];

      // We insert line segments which cover order_v - 1 intervals.
      for (size_t kk = 1; kk < num_segs_w - order_w; kk+=order_w+1)
	{
	  double tmin2 = (1) ? unique_knots_w[kk] : unique_knots_w[kk+1];
	  double tmax2 = (1) ? unique_knots_w[kk+order_w] : unique_knots_w[kk+order_w+1];
	  for (size_t kj = 1; kj < num_segs_v - order_v; kj+=order_v+1)
	    {
//	      assert(unique_knots_v.size() > 3);
	      // To avoid a too regular pattern we alter start
	      // values. Maybe we should use a random function ...
	      double tmin1 = (1) ? unique_knots_v[kj] : unique_knots_v[kj+1];
	      double tmax1 = (1) ? unique_knots_v[kj+order_v] : unique_knots_v[kj+order_v+1];

	      LRSplineVolume::Refinement3D ref;
	      ref.kval = ref_par;
	      ref.start1 = tmin1;
	      ref.end1 = tmax1;
	      ref.start2 = tmin2;
	      ref.end2 = tmax2;
	      ref.d = dir;
	      ref.multiplicity = mult;
	      refs_u.push_back(ref); // Refinements in the u-dir.
	    }
	}
    }

  dir = YDIR; // XDIRFIXED, refining in the x-dir (u).
  if ((num_segs_u > order_u - 1) && (num_segs_w > order_w - 1))
    for (size_t ki = 0; ki < num_segs_v; ++ki)
      {
	// We bisect the interval.
	double wgt = 0.4; // Value in the range [0.0, 1.0].
	double ref_par = wgt*unique_knots_v[ki] + (1.0 - wgt)*unique_knots_v[ki+1];

	// We insert line segments which cover order_v - 1 intervals.
	for (size_t kk = 1; kk < num_segs_u - order_u; kk+=order_u+1)
	  {
	    double tmin2 = (1) ? unique_knots_u[kk] : unique_knots_u[kk+1];
	    double tmax2 = (1) ? unique_knots_u[kk+order_u] : unique_knots_u[kk+order_u+1];
	    for (size_t kj = 1; kj < num_segs_w - order_w; kj+=order_w+1)
	      {
//		assert(unique_knots_w.size() > 3);
		// To avoid a too regular pattern we alter start
		// values. Maybe we should use a random function ...
		double tmin1 = (1) ? unique_knots_w[kj] : unique_knots_w[kj+1];
		double tmax1 = (1) ? unique_knots_w[kj+order_w] : unique_knots_w[kj+order_w+1];

		LRSplineVolume::Refinement3D ref;
		ref.kval = ref_par;
		ref.start1 = tmin1;
		ref.end1 = tmax1;
		ref.start2 = tmin2;
		ref.end2 = tmax2;
		ref.d = dir;
		ref.multiplicity = mult;
		refs_v.push_back(ref); // Refinements in the u-dir.
	      }
	  }
      }

  dir = ZDIR; // Refining in the Z-dir (w).
  if ((num_segs_u > order_u - 1) && (num_segs_v > order_v - 1))
    for (size_t ki = 0; ki < num_segs_w; ++ki)
      {
	// We bisect the interval.
	double wgt = 0.4; // Value in the range [0.0, 1.0].
	double ref_par = wgt*unique_knots_w[ki] + (1.0 - wgt)*unique_knots_w[ki+1];

	// We insert line segments which cover order_w - 1 intervals.
	for (size_t kk = 1; kk < num_segs_v - order_v; kk+=order_v+1)
	  {
	    double tmin2 = (1) ? unique_knots_v[kk] : unique_knots_v[kk+1];
	    double tmax2 = (1) ? unique_knots_v[kk+order_v] : unique_knots_v[kk+order_v+1];
	    for (size_t kj = 1; kj < num_segs_u - order_u; kj+=order_u+1)
	      {
//		assert(unique_knots_u.size() > 3);
		// To avoid a too regular pattern we alter start
		// values. Maybe we should use a random function ...
		double tmin1 = (1) ? unique_knots_u[kj] : unique_knots_u[kj+1];
		double tmax1 = (1) ? unique_knots_u[kj+order_u] : unique_knots_u[kj+order_u+1];

		LRSplineVolume::Refinement3D ref;
		ref.kval = ref_par;
		ref.start1 = tmin1;
		ref.end1 = tmax1;
		ref.start2 = tmin2;
		ref.end2 = tmax2;
		ref.d = dir;
		ref.multiplicity = mult;
		refs_w.push_back(ref); // Refinements in the u-dir.
	      }
	  }
      }

    vector<LRSplineVolume::Refinement3D> all_refs(refs_u.begin(), refs_u.end());
    all_refs.insert(all_refs.end(), refs_v.begin(), refs_v.end());
    all_refs.insert(all_refs.end(), refs_w.begin(), refs_w.end());
    std::cout << "num_refs_to_insert: " << all_refs.size() << std::endl;
    std::cout << "num_refs_u_to_insert: " << refs_u.size() << std::endl;
    std::cout << "num_refs_v_to_insert: " << refs_v.size() << std::endl;
    std::cout << "num_refs_w_to_insert: " << refs_w.size() << std::endl;

    // if (all_refs.size() > 2)
    // 	all_refs.erase(all_refs.begin(), all_refs.end() - 2);
    // std::cout << "num_refs_to_insert: " << all_refs.size() << std::endl;

    // We then create and refine the LRSplineVolume.
    lr_spline_vol = shared_ptr<LRSplineVolume>
      // (new LRSpline<vector<double>::const_iterator, vector<double>::const_iterator> >
      (new LRSplineVolume(order_u - 1, order_v - 1,  order_w - 1,
			  num_coefs_u, num_coefs_v, num_coefs_w,
			  dim,
			  spline_vol->basis(0).begin(),
			  spline_vol->basis(1).begin(),
			  spline_vol->basis(2).begin(),
			  spline_vol->coefs_begin()));

#if 1
    lr_spline_vol_single_refs = shared_ptr<LRSplineVolume>
      // (new LRSpline<vector<double>::const_iterator, vector<double>::const_iterator> >
      (new LRSplineVolume(order_u - 1, order_v - 1, order_w - 1,
			  num_coefs_u, num_coefs_v, num_coefs_w,
			  dim,
			  spline_vol->basis(0).begin(),
			  spline_vol->basis(1).begin(),
			  spline_vol->basis(2).begin(),
			  spline_vol->coefs_begin()));
#endif

    bool refine_multi = true;
    if (refine_multi)
      {
	int num_elem = lr_spline_vol->numElements();
	int num_basis_funcs = lr_spline_vol->numBasisFunctions();
	std::cout << "num_elem: " << num_elem << ", num_basis_funcs: " << num_basis_funcs << std::endl;
	double time_lrspline_lr_ref = benchmarkVolRefinement(*lr_spline_vol, all_refs);
	std::cout << "Time lr refinement multi refs: " << time_lrspline_lr_ref << std::endl;
	// We write to screen the number of basis functions for both
	// versions.
	num_elem = lr_spline_vol->numElements();
	num_basis_funcs = lr_spline_vol->numBasisFunctions();
	std::cout << "num_elem: " << num_elem << ", num_basis_funcs: " << num_basis_funcs << std::endl;

	std::ofstream line_cl_out("tmp/line_cl.g2");
	writeElementLineCloud(*lr_spline_vol, line_cl_out);

	double max_dist = maxDist(spline_vol.get(), *lr_spline_vol, num_samples_u, num_samples_v, num_samples_w);
	std::cout << "Max dist between input and (multi) refined surface: " << max_dist << std::endl;

	std::ofstream fileout("tmp/ref_lr_multi.g2");
	lr_spline_vol->writeStandardHeader(fileout);
	lr_spline_vol->write(fileout);

	// MESSAGE("Missing writing refined surface grid to file!");
	std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");
//  writePostscriptMesh(*lrsf);
	writePostscriptMesh(*lr_spline_vol, lrsf_grid_ps);

      }

#if 1
    bool refine_single = true;
    if (refine_single)
      {
	double time_lrspline_lr_ref_single = benchmarkVolRefinement(*lr_spline_vol_single_refs, all_refs, true);
	std::cout << "Time lr refinement single refs: " << time_lrspline_lr_ref_single << std::endl;
	int num_elem_single = lr_spline_vol_single_refs->numElements();
	int num_basis_funcs_single = lr_spline_vol_single_refs->numBasisFunctions();
	std::cout << "num_elem_single_refs: " << num_elem_single <<
	  ", num_basis_funcs_single_refs: " << num_basis_funcs_single << std::endl;

#if 1//ndef NDEBUG
	vector<LRBSpline3D*> bas_funcs;
	for (auto iter = lr_spline_vol_single_refs->basisFunctionsBegin(); iter != lr_spline_vol_single_refs->basisFunctionsEnd(); ++iter)
	  {
	    bas_funcs.push_back((*iter).second.get());
	  }
	vector<Element3D*> elems;
	for (auto iter = lr_spline_vol_single_refs->elementsBegin(); iter != lr_spline_vol_single_refs->elementsEnd(); ++iter)
	  {
	    elems.push_back(((*iter).second.get()));
	  }
	puts("Remove when done debugging!");
#endif

	double max_dist_post_ref_multi_ref =
	  maxDist(spline_vol.get(), *lr_spline_vol_single_refs, num_samples_u, num_samples_v, num_samples_w);
	std::cout << "Max dist input and (single) ref surface: " << max_dist_post_ref_multi_ref << std::endl;

	std::ofstream fileout2("tmp/ref_lr_single.g2");
	lr_spline_vol_single_refs->writeStandardHeader(fileout2);
	lr_spline_vol_single_refs->write(fileout2);
      }
#endif

}



double maxDist(const Go::ParamVolume* param_vol,
	       const Go::LRSplineVolume& lr_spline_vol,
	       int num_samples_u, int num_samples_v, int num_samples_w)
{
  // Assuming the domain is the same.
  double umin = lr_spline_vol.startparam_u();
  double umax = lr_spline_vol.endparam_u();
  double vmin = lr_spline_vol.startparam_v();
  double vmax = lr_spline_vol.endparam_v();
  double wmin = lr_spline_vol.startparam_w();
  double wmax = lr_spline_vol.endparam_w();
  double ustep = (umax - umin)/((double)num_samples_u - 1);
  double vstep = (vmax - vmin)/((double)num_samples_v - 1);
  double wstep = (wmax - wmin)/((double)num_samples_w - 1);
  Go::Point go_pt(3), lr_pt(3);
  double max_dist = -1.0;
// #ifndef NDEBUG
  double max_dist_u = 0.0;
  double max_dist_v = 0.0;
  double max_dist_w = 0.0;
  Go::Point max_go_pt(3), max_lr_pt(3);
// #endif
  for (int kk = 0; kk < num_samples_w; ++kk)
    {
      double wpar = wmin + kk*wstep;
      for (int kj = 0; kj < num_samples_v; ++kj)
	{
	  double vpar = vmin + kj*vstep;
	  for (int ki = 0; ki < num_samples_u; ++ki)
	    {
	      double upar = umin + ki*ustep;
	      param_vol->point(go_pt, upar, vpar, wpar);
	      lr_spline_vol.point(lr_pt, upar, vpar, wpar);
	      double dist = go_pt.dist(lr_pt);
	      if (dist > max_dist)
		{
		  max_dist = dist;
// #ifndef NDEBUG
		  max_dist_u = upar;
		  max_dist_v = vpar;
		  max_dist_w = wpar;
		  max_go_pt = go_pt;
		  max_lr_pt = lr_pt;
// #endif
		}
	    }
	}    
    }
// #ifndef NDEBUG
  std::cout << "max_dist: u = " << max_dist_u << ", v = " << max_dist_v <<  ", w = " << max_dist_w << std::endl;
  std::cout << "max_go_pt = " << max_go_pt << std::endl;
  std::cout << "max_lr_pt = " << max_lr_pt << std::endl;
// #endif

    return max_dist;


}
