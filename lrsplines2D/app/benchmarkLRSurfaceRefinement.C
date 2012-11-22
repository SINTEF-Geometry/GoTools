//===========================================================================
//                                                                           
// File: benchmarkLRSurfaceRefinement.C                                    
//                                                                           
// Created: Wed Nov 22 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision:
//                                                                           
// Description: Read a LRSplineSurface, benchmark the refinement function.
//              We compare the one-at-a-time approach vs all at once.
//                                                                           
//===========================================================================

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


int main(int argc, char *argv[])
{
  if (argc != 4)
  {
      std::cout << "Usage: spline_sf.g2 num_ref num_iter" << std::endl;
      return -1;
  }

  char* filein_char = argv[1];
  std::ifstream filein(argv[1]); // Input lr spline.
  // If input surface is not fine enough we refine.
  int min_num_each_dir = atoi(argv[2]);
  int num_iter = atoi(argv[3]);
  if (num_iter != 1)
      MESSAGE("Not using num_iter yet.");

  shared_ptr<Go::SplineSurface> spline_sf;
//  shared_ptr<Go::SplineSurface> lr_spline_sf_go;
  shared_ptr<LRSplineSurface> lr_spline_sf;

  int order_u, order_v, num_coefs_u, num_coefs_v, dim, num_bases=-1;
  if (strcasestr(filein_char, ".g2"))
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

#if 0
#endif
      puts("Done reading the Go::SplineSurface.");
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
  Direction2D dir = XFIXED;
  for (size_t ki = 0; ki < num_segs_u - 1; ++ki)
  {
      // We bisect the interval.
      double wgt = 0.4; // Value in the range [0.0, 1.0].
      double const_par = wgt*unique_knots_u[ki] + (1.0 - wgt)*unique_knots_u[ki+1];

      // We insert line segments which cover order_v - 1 intervals.
      for (size_t kj = order_v - 1; kj < num_segs_v; kj += order_v)
      {
	  // To avoid a too regular pattern we alter start
	  // values. Maybe we should use a random function ...
	  double tmin = (ki%2 == 1) ? unique_knots_v[kj - order_v + 1] : unique_knots_v[kj - order_v + 2];
	  double tmax = (ki%2 == 1) ? unique_knots_v[kj] : unique_knots_v[kj + 1];

	  LRSplineSurface::Refinement2D ref;
	  ref.kval = const_par;
	  ref.start = tmin;
	  ref.end = tmax;
	  ref.d = dir;
	  ref.multiplicity = mult;
	  refs_u.push_back(ref);
      }
  }

  dir = YFIXED;
  for (size_t ki = 0; ki < num_segs_v - 1; ++ki)
  {
      // We bisect the interval.
      double wgt = 0.4;
      double const_par = wgt*unique_knots_v[ki] + (1.0 - wgt)*unique_knots_v[ki+1];

      // We insert line segments which cover order_u - 1 intervals.
      for (size_t kj = order_u - 1; kj < num_segs_u; kj += order_u)
      {
	  // To avoid a too regular pattern we alter start
	  // values. Maybe we should use a random function ...
	  double tmin = (ki%2 == 1) ? unique_knots_u[kj - order_u + 1] : unique_knots_u[kj - order_u + 2];
	  double tmax = (ki%2 == 1) ? unique_knots_u[kj] : unique_knots_u[kj + 1];

	  LRSplineSurface::Refinement2D ref;
	  ref.kval = const_par;
	  ref.start = tmin;
	  ref.end = tmax;
	  ref.d = dir;
	  ref.multiplicity = mult;
	  refs_v.push_back(ref);
      }
  }

  vector<LRSplineSurface::Refinement2D> all_refs(refs_u.begin(), refs_u.end());
  all_refs.insert(all_refs.end(), refs_v.begin(), refs_v.end());
  std::cout << "num_refs_to_insert: " << all_refs.size() << std::endl;


  bool refine_lr_olr = true;
  if (refine_lr_olr)
  {
      lr_spline_sf = shared_ptr<LRSplineSurface>
	  // (new LRSpline<vector<double>::const_iterator, vector<double>::const_iterator> >
	  (new LRSplineSurface(order_u - 1, order_v - 1,
				 num_coefs_u, num_coefs_v,
				 dim,
				 spline_sf->basis_u().begin(),
				 spline_sf->basis_v().begin(),
				 spline_sf->coefs_begin()));

      double time_lrspline_olr_ref = benchmarkSfRefinement(*lr_spline_sf, all_refs);
      std::cout << "Time olr refinement: " << time_lrspline_olr_ref << std::endl;
  }


  if (lr_spline_sf.get() != NULL)
  {
      // We write to screen the number of basis functions for both
      // versions.
      int num_basis_funcs = lr_spline_sf->numBasisFunctions();
      std::cout << ", num_basis_funcs: " << num_basis_funcs << std::endl;

      int num_samples_u = min_num_each_dir + 29; // Rather random number.
      int num_samples_v = min_num_each_dir + 37; // Rather random number.
      double umin = spline_sf->startparam_u();
      double umax = spline_sf->endparam_u();
      double ustep = (umax - umin)/(num_samples_u - 1);
      double vmin = spline_sf->startparam_v();
      double vmax = spline_sf->endparam_v();
      double vstep = (vmax - vmin)/(num_samples_v - 1);
      double max_dist = -1.0;
      Point pt_go;
      Point pt_lr;
      for (size_t kj = 0; kj < num_samples_v; ++kj)
      {
	  double vpar = vmin + kj*vstep;
	  for (size_t ki = 0; ki < num_samples_u; ++ki)
	  {
	      double upar = umin + ki*ustep;
	      spline_sf->point(pt_go, upar, vpar);	  
	      pt_lr = (*lr_spline_sf)(upar, vpar);
	      double dist = 0.0;
	      for (int kk = 0; kk < dim; ++kk)
		  dist += (pt_go[kk] - pt_lr[kk])*(pt_go[kk] - pt_lr[kk]);

	      dist = sqrt(dist);
	      if (dist > max_dist)
	      {
		  max_dist = dist;
	      }
	  }
      }
      std::cout << "max_dist: " << max_dist << std::endl;
  }

 // MESSAGE("Missing writing refined surface grid to file!");
  std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");
//  writePostscriptMesh(*lrsf);
  writePostscriptMesh(*lr_spline_sf, lrsf_grid_ps);

}
