//===========================================================================
//                                                                           
// File: benchmarkLRSurfacePointEval.C                                       
//                                                                           
// Created: Tue Jan 29 13:57:16 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description: Benchmark the point evaluation (with higher derivs)
//              for a LRSplineSurface.  If input is a regular
//              SplineSurface, we convert to LRSplineSurface and
//              compare. If input is a LRSplineSurface, we compare
//              with the corresponding (refined) SplineSurface.
//                                                                           
//===========================================================================


#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/Direction2D.h"
#include "GoTools/lrsplines2D/LRBenchmarkUtils.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/utils/config.h"


#include <iostream>
#include <assert.h>
#include <fstream>
#include <string>
#include "string.h"


using Go::ObjectHeader;
using Go::SplineSurface;
using Go::LRSplineSurface;
using std::vector;
using namespace Go;


// Assuming the domain is the same.
double maxDist(const SplineSurface& spline_sf,
	       const LRSplineSurface& lr_spline_sf,
	       int num_samples_u, int num_samples_v,
	       int der_u, int der_v);

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
      std::cout << "Usage: (lr)spline_sf.g2 num_dir_samples sum_derivs num_iter" << std::endl;
      return -1;
  }

  char* filein_char = argv[1];
  std::ifstream filein(argv[1]); // Input lr spline.
  // If input surface is not fine enough we refine.
  int num_dir_samples = atoi(argv[2]);
  int sum_derivs = atoi(argv[3]);
  int num_iter = atoi(argv[4]);
  if (num_iter != 1)
      puts("Not using num_iter yet.");

  shared_ptr<SplineSurface> spline_sf;
  shared_ptr<LRSplineSurface> lr_spline_sf;
  double knot_tol = 1e-10;

  // If we refine the lr-surface in order to extract a regular
  // spline-version we must perform the operations on a copy.
  shared_ptr<LRSplineSurface> lr_spline_sf_copy;

  int order_u, order_v, num_coefs_u, num_coefs_v, dim, num_bases=-1;
  if (strcasestr(filein_char, ".g2"))
    {
      ObjectHeader header;
      filein >> header;
      if (header.classType() == Class_SplineSurface)
	{
	  std::cout << "Input was a SplineSurface, creating a LRSplineSurface." << std::endl;
	  spline_sf = shared_ptr<SplineSurface>(new SplineSurface());
	  filein >> *spline_sf;
	  // We create the lr-version.
	  lr_spline_sf = shared_ptr<LRSplineSurface>(new LRSplineSurface(spline_sf.get(), knot_tol));
	}
      else if (header.classType() == Class_LRSplineSurface)
	{
	  std::cout << "Input was a LRSplineSurface, creating a SplineSurface (refining)." << std::endl;
	  lr_spline_sf = shared_ptr<LRSplineSurface>(new LRSplineSurface());
	  filein >> *spline_sf;
	  lr_spline_sf_copy = shared_ptr<LRSplineSurface>(lr_spline_sf->clone());

	}
      else
	{
	  std::cout << "Input was not a SplineSurface or a LRSplineSurface, exiting!" << std::endl;
	}
    }
  else
    {
      MESSAGE("Input was not a g2-file!");
      return -1;
    }

  for (int kj = 0; kj < sum_derivs + 1; ++kj)
    for (int ki = 0; ki < sum_derivs + 1 - kj; ++ki)
      {
	double max_dist = maxDist(*spline_sf, *lr_spline_sf, num_dir_samples, num_dir_samples, ki, kj);

	std::cout << "Max dist with der_u=" << ki << " & der_v=" << kj << ": " << max_dist << std::endl;
      }
}


double maxDist(const SplineSurface& spline_sf,
	       const LRSplineSurface& lr_spline_sf,
	       int num_samples_u, int num_samples_v,
	       int der_u, int der_v)
{
  const int derivs = std::max(der_u, der_v);
  const int sum_derivs = der_u + der_v;
  const int num_pts = (sum_derivs + 1)*(sum_derivs + 2)/2;
  const int first_pos = (sum_derivs + 1)*sum_derivs/2;
  const int der_pos = first_pos + der_v;
  // Assuming the domain is the same.
  double umin = lr_spline_sf.startparam_u();
  double umax = lr_spline_sf.endparam_u();
  double vmin = lr_spline_sf.startparam_v();
  double vmax = lr_spline_sf.endparam_v();
  double ustep = (umax - umin)/((double)num_samples_u - 1);
  double vstep = (vmax - vmin)/((double)num_samples_v - 1);
  Point go_pt(3), lr_pt(3);
  double max_dist = -1.0;
// #ifndef NDEBUG
  double max_dist_u = 0.0;
  double max_dist_v = 0.0;
  Point max_go_pt(3), max_lr_pt(3);
  vector<Point> go_pts(num_pts), lr_pts(num_pts);
// #endif
  for (int kj = 0; kj < num_samples_v; ++kj)
    {
      double vpar = vmin + kj*vstep;
      for (int ki = 0; ki < num_samples_u; ++ki)
	{
	  double upar = umin + ki*ustep;
	  spline_sf.point(go_pts, upar, vpar, sum_derivs);
	  go_pt = go_pts[der_pos];
//	  lr_spline_sf.point(lr_pts, upar, vpar, sum_derivs);
	  lr_pt = lr_spline_sf(upar, vpar, der_u, der_v);
//lr_pts[der_pos];
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
