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
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: spline (.g2) lrspline_out.g2 curve_out.g2" << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  std::ofstream fileout2(argv[3]);

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

  puts("Writing initial lr-spline to file.");
  lrsf->writeStandardHeader(fileout);
  lrsf->write(fileout);

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

  for (auto it = lrsf->basisFunctionsBegin(); it != lrsf->basisFunctionsEnd(); ++it)
    {
      std::cout << &it->second << "  " << it->second->umin() << "  ";
      std::cout << it->second->umax() << "  " << it->second->vmin() << "  ";
      std::cout << it->second->vmax() << std::endl;
    }

  for (auto it = lrsf->elementsBegin(); it != lrsf->elementsEnd(); ++it)
    {
      int nmb = it->second->nmbBasisFunctions();
      std::cout << "Nmb basis: " << nmb << "  " << it->second->umin() << "  ";
      std::cout << it->second->umax() << "  " << it->second->vmin() << "  ";
      std::cout << it->second->vmax() << std::endl;

      for (int ki=0; ki<nmb; ++ki)
	{
	  LRBSpline2D* b = it->second->supportFunction(ki);
	  std::cout << b << "  ";
	}
      std::cout << std::endl;
    }

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

  CurveLoop loop = lrsf->outerBoundaryLoop();;

  puts("Writing lr-spline to file.");
  lrsf->writeStandardHeader(fileout);
  lrsf->write(fileout);

  int nmb = loop.size();
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamCurve> cv = loop[ki];
      cv->writeStandardHeader(fileout2);
      cv->write(fileout2);
    }

#ifndef NDEBUG
  std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");
//  writePostscriptMesh(*lrsf);
  writePostscriptMesh(*lrsf, lrsf_grid_ps);
#endif NDEBUG

  return 0;
}
