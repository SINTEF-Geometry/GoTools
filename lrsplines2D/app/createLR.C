//===========================================================================
//                                                                           
// File: createLR
//                                                                           
// Created: 20.10.12
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description: Create LRSplineSurface from SplineSurface
//                                                                           
//===========================================================================


#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: spline (.g2) lrspline_out.g2" << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);

  ObjectHeader header;
  header.read(filein);
  SplineSurface splinesf;
  splinesf.read(filein);

  double knot_tol = 1.0e-6;

  shared_ptr<LRSplineSurface> lrsf(new LRSplineSurface(&splinesf, 
						       knot_tol));
				   
  // We write to screen various attributes.
  int num_basis_funcs = lrsf->numBasisFunctions();
  std::cout << "num_basis_funcs: " << num_basis_funcs << std::endl;

  std::cout << "New knot, give pardir (0,1), fixed value, start and end: ";
  std::cout << std::endl;
  double parval, start, end;
  int dir;
  int mult = 1;
#if 1
  std::cin >> dir;
  std::cin >> parval;
  std::cin >> start;
  std::cin >> end;
#else
  // Hardcoded values for debugging.
  dir = 0;
  parval = 0.2;
  start = 0.5;
  end = 1.0;
#endif
  lrsf->refine((dir==0) ? YFIXED : XFIXED, parval, start, end, mult);

  num_basis_funcs = lrsf->numBasisFunctions();
  std::cout << "num_basis_funcs: " << num_basis_funcs << std::endl;

  std::cout << "Evaluate: Give parameter values: ";
  double paru, parv;
  std::cin >> paru;
  std::cin >> parv;
  Point pos;
  lrsf->point(pos, paru, parv);
  std::cout << "Value (lr): " << pos << std::endl;

  Point pos2;
  splinesf.point(pos2, paru, parv);
  std::cout << "Value (spline): " << pos2 << std::endl;

  puts("Writing lr-spline to file.");
  lrsf->write(fileout);

#ifndef NDEBUG
  std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");
//  writePostscriptMesh(*lrsf);
  writePostscriptMesh(*lrsf, lrsf_grid_ps);
#endif NDEBUG

  return 0;
}
