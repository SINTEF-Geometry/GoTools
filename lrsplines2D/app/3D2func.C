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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: input surface(.g2), output surface(.g2)" << std::endl;
    return -1;
  }

  std::ifstream input(argv[1]);
  std::ofstream output(argv[2]);

  ObjectHeader header2;
  header2.read(input);
  shared_ptr<LRSplineSurface> sf(new LRSplineSurface());
  sf->read(input);

  if (sf->dimension() == 1)
    {
      std::cout << "Already 1D" << std::endl;
        sf->writeStandardHeader(output);
	sf->write(output);
	return 1;
    }
  else if (sf->dimension() != 3)
    {
      std::cout << "Not a 3D surface" << std::endl;
      return -1;
    }

  // Pass through all B-splines and change the coefficient to
  // use only the z-value. Only for surfaces that previously has
  // been lifted to 3D
  
  LRSplineSurface::BSplineMap::iterator it1 = 
    sf->basisFunctionsBeginNonconst();
  LRSplineSurface::BSplineMap::iterator it2 = 
    sf->basisFunctionsEndNonconst();
  for (; it1!=it2; ++it1)
    {
      LRBSpline2D *bspline = it1->second.get();
      Point coef1 = bspline->Coef();
      double gamma = bspline->gamma();
      Point coef2(1);
      coef2.setValue(coef1[2]);
      bspline->setCoefAndGamma(coef2, gamma);
    }

  sf->writeStandardHeader(output);
  sf->write(output);
}
