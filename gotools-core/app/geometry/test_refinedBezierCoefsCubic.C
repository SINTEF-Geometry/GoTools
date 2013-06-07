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
