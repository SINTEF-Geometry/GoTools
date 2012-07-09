//===========================================================================
//                                                                           
// File: test_refinedBezierCoefsCubic.C                                      
//                                                                           
// Created: Fri Jun 22 14:55:16 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description: We test the routine refinedBezierCoefs().
//              The coefs are extracted, we create a sub-surface, and
//              compare some sample pts.
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineUtils.h"


#include <fstream>
#include <vector>

using std::vector;
using namespace Go;


int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: spline_sf (.g2) knot_ind_u_min knot_ind_v_min" << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]); // Input lr spline.
  int ind_u_min = atoi(argv[2]);
  int ind_v_min = atoi(argv[3]);

  ObjectHeader header;
  header.read(filein);
  SplineSurface spline_sf;
  spline_sf.read(filein);

  // We write to screen some info about the spline.
  int dim = spline_sf.dimension();
  int order_u = spline_sf.order_u();
  int order_v = spline_sf.order_v();
  int num_coefs_u = spline_sf.numCoefs_u();
  int num_coefs_v = spline_sf.numCoefs_v();
  std::cout << "Spline sf info:  dim = " << dim << ", order_u = " << order_u << ", order_v = " << order_v << std::endl;
  std::cout << "                 num_coefs_u = " << num_coefs_u << ", num_coefs_v = " << num_coefs_v << std::endl;


  shared_ptr<SplineSurface> bez_sf;
  double umin, umax, vmin, vmax;
  if (ind_u_min == -1 && ind_v_min == -1)
  {
      puts("We use the global bezier refinement matrix code!");
      bez_sf = SplineUtils::refineToBezier(spline_sf);
      umin = bez_sf->startparam_u();
      umax = bez_sf->endparam_u();
      vmin = bez_sf->startparam_v();
      vmax = bez_sf->endparam_v();
  }
  else
  {
      const BsplineBasis& basis_u = spline_sf.basis_u();
      vector<double>::const_iterator iter = basis_u.begin() + ind_u_min;
      umin = *iter;
      while (*iter == umin)
	  ++iter;
      umax = *iter;
      ind_u_min = iter - basis_u.begin() - 1;

      const BsplineBasis& basis_v = spline_sf.basis_v();
      iter = basis_v.begin() + ind_v_min;
      vmin = *iter;
      while (*iter == vmin)
	  ++iter;
      vmax = *iter;
      ind_v_min = iter - basis_v.begin() - 1;

      vector<double> bez_coefs, bez_coefs2;
      Go::SplineUtils::refinedBezierCoefsCubic(spline_sf,
					       ind_u_min, ind_v_min,
					       bez_coefs);

      int num_bez_coefs_u = order_u;
      int num_bez_coefs_v = order_v;
      vector<double> knots_u(order_u, umin);
      knots_u.insert(knots_u.end(), order_u, umax);
      vector<double> knots_v(order_v, vmin);
      knots_v.insert(knots_v.end(), order_v, vmax);
      bez_sf = shared_ptr<SplineSurface>
	  (new SplineSurface(num_bez_coefs_u, num_bez_coefs_v,
			     order_u, order_v,
			     knots_u.begin(), knots_v.begin(),
			     bez_coefs.begin(),
			     dim));

#if 0
      refinedBezierCoefs(spline_sf,
			 ind_u_min, ind_v_min,
			 bez_coefs2);
      SplineSurface bez_sf2(num_bez_coefs_u, num_bez_coefs_v,
			    order_u, order_v,
			    knots_u.begin(), knots_v.begin(),
			    bez_coefs2.begin(),
			    dim);
#endif
  }

  // And finally we compare the distance in a number of sample params.
  const int num_samples_u = 137;
  const int num_samples_v = 143;
  double ustep = (umax - umin)/(num_samples_u -1);
  double vstep = (vmax - vmin)/(num_samples_v -1);
  Point spline_pt(dim), bez_pt(dim), bez_pt2(dim);
  double max_dist = -1.0;
  double max_dist2 = -1.0;
  for (size_t kj = 0; kj < num_samples_v; ++kj)
  {
      double vpar = vmin + kj*vstep;
      for (size_t ki = 0; ki < num_samples_u; ++ki)
      {
	  double upar = umin + ki*ustep;
	  spline_sf.point(spline_pt, upar, vpar);
	  bez_sf->point(bez_pt, upar, vpar);
	  double dist = spline_pt.dist(bez_pt);
	  if (dist > max_dist)
	      max_dist = dist;
#if 0
	  bez_sf2.point(bez_pt2, upar, vpar);
	  double dist2 = spline_pt.dist(bez_pt2);
	  if (dist2 > max_dist2)
	      max_dist2 = dist2;
#endif

      }
  }

  std::cout << "max_dist: " << max_dist << std::endl;
#if 0
  std::cout << "max_dist2: " << max_dist2 << std::endl;
#endif

}
