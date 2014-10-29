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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineEvalGrid.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace Go;
using std::vector;


typedef Element2D simpleElement;


int main(int argc, char *argv[])
{
  if (argc != 3)
  {
      std::cout << "Usage: input_(lr_)spline_sf.g2 ref_lr_spline_sf.g2" << std::endl;
      return -1;
  }

  std::ifstream filein(argv[1]); // Input LRSplineSurface.
  std::ofstream fileout(argv[2]); // Output: The sampled points.

  Go::ObjectHeader header;
  filein >> header;

  LRSplineSurface lrspline_sf;
  lrspline_sf.read(filein);

  std::cout << "Done reading lrspline_sf." << std::endl;

  if (lrspline_sf.dimension() == 1)
  {
      MESSAGE("Lifint the 1D surface to 3D.");
      lrspline_sf.to3D();
  }


  LRSplineEvalGrid eval_grid(lrspline_sf);

  std::cout << "Done creating grid eval." << std::endl;

  int order_u = eval_grid.orderU();
  int order_v = eval_grid.orderU();
  int dim = eval_grid.dim();
  int num_elem = eval_grid.numElements();

  vector<double> grid_pts;
  grid_pts.reserve(dim*order_u*order_v*num_elem);
  double umin, umax, vmin, vmax, ustep, vstep, upar, vpar;
  vector<double> res(dim);
  int cntr = 0;
  double shrink_x_y = 1e-03;
  for (auto iter = eval_grid.elements_begin(); iter != eval_grid.elements_end(); ++iter)
    {
      // std::cout << "cntr: " << cntr << std::endl;
      // ++cntr;

      eval_grid.low(*iter, umin, vmin);
      eval_grid.high(*iter, umax, vmax);
      ustep = (umax - umin)/(order_u - 1);
      vstep = (vmax - vmin)/(order_v - 1);
      for (int kj = 0; kj < order_v; ++kj)
	{
	  vpar = vmin + kj*vstep;
	  for (int ki = 0; ki < order_u; ++ki)
	    {
	      upar = umin + ki*ustep;
	      eval_grid.evaluate(*iter, upar, vpar, &res[0]);
	      if (dim == 3)
		{
		  res[0] *= shrink_x_y;
		  res[1] *= shrink_x_y;
		}
	      grid_pts.insert(grid_pts.end(), res.begin(), res.end());
	    }
	}
    }

  int num_pts = grid_pts.size()/dim;
  Go::PointCloud3D pt_cl(grid_pts.begin(), num_pts);
  pt_cl.writeStandardHeader(fileout);
  pt_cl.write(fileout);

  eval_grid.testCoefComputation();

  return 0;
}
