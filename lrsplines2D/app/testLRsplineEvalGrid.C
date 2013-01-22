//===========================================================================
//                                                                           
// File: testLRsplineEvalGrid.C                                              
//                                                                           
// Created: Tue Jan 22 10:37:09 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


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

  LRSplineEvalGrid eval_grid(lrspline_sf);

  std::cout << "Done creating grid eval." << std::endl;

  int order_u = eval_grid.orderU();
  int order_v = eval_grid.orderU();
  int dim = eval_grid.dim();
  int num_elem = eval_grid.numElements();

  vector<float> grid_pts;
  grid_pts.reserve(dim*order_u*order_v*num_elem);
  float umin, umax, vmin, vmax, ustep, vstep, upar, vpar;
  vector<float> res(dim);
  int cntr = 0;
  for (auto iter = eval_grid.elements_begin(); iter != eval_grid.elements_end(); ++iter)
    {
      std::cout << "cntr: " << cntr << std::endl;
      ++cntr;

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
	      grid_pts.insert(grid_pts.end(), res.begin(), res.end());
	    }
	}
    }

  int num_pts = grid_pts.size()/dim;
  Go::PointCloud3D pt_cl(grid_pts.begin(), num_pts);
  pt_cl.writeStandardHeader(fileout);
  pt_cl.write(fileout);

  return 0;
}
