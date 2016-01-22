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

  // double umin, umax, vmin, vmax;
  // std::cout << " Sub surface, give parameters: ";
  // std::cin >> umin;
  // std::cin >> umax;
  // std::cin >> vmin;
  // std::cin >> vmax;
  // shared_ptr<LRSplineSurface> lrsub(lrsf->subSurface(umin, vmin, umax, vmax,
  // 						     knot_tol));
  // lrsub->writeStandardHeader(fileout);
  // lrsub->write(fileout);

#ifndef NDEBUG
  std::ofstream lrsf_grid_ps("tmp/lrsf_grid.ps");
//  writePostscriptMesh(*lrsf);
  writePostscriptMesh(*lrsf, lrsf_grid_ps);
#endif // NDEBUG

  return 0;
}
