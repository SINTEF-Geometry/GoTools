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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include <fstream>

using namespace Go;
using std::vector;

//===========================================================================
//                                                                           
/// Description:
/// 
//                                                                           
//===========================================================================

int main(int argc, char *argv[])
{

  int dim = 1;
  int in1 = 11, in2 = 9;
  double deg1 = 4, deg2 = 3;
  double et1[16] = {0,0,0,0,0,1,2,3,4,5,6,7,7,7,7,7};
  double et2[13] = {0,0,0,0,1,2,3,4,5,6,6,6,6};
  vector<double> coefs(in1*in2*dim, 1.0);
  // Create spline surface
  shared_ptr<SplineSurface> surf(new SplineSurface(in1, in2,
						   deg1+1, deg2+1,
						   et1, et2,
						   &coefs[0], dim));

  shared_ptr<LRSplineSurface> lr_surf(new LRSplineSurface(surf.get(),
							  1.0e-6));
  
  // Refine
  LRSplineSurface::Refinement2D refs[5];
  refs[0].setVal(2.5, 0.0, 4.0, XFIXED, 1);
  refs[1].setVal(3.5, 2.0, 6.0, XFIXED, 1);
  refs[2].setVal(4.5, 0.0, 4.0, XFIXED, 1);
  refs[3].setVal(2.5, 2.0, 7.0, YFIXED, 1);
  refs[4].setVal(3.5, 0.0, 5.0, YFIXED, 1);
  
  for (int ki=0; ki<5; ++ki)
    lr_surf->refine(refs[ki], true);

  std::ofstream of1("lr_mesh.eps");
  writePostscriptMesh(*lr_surf, of1);
  
  // Compute weights by nesting level
  int max_level = 10;  // Should always be enough
  for (int level=0; level<max_level; ++level)
    {
      int num_update = 0;
      for (auto bspl=lr_surf->basisFunctionsBegin();
	   bspl!=lr_surf->basisFunctionsEnd(); ++bspl)
	{
	  int blevel = bspl->second->getNestLevel();
	  if (blevel != level)
	    continue;

	  std::cout << "B-spline level " << blevel << ", knots: " << std::endl;
	  vector<int> kvec1 = bspl->second->kvec(XFIXED);
	  for (size_t kj=0; kj<kvec1.size(); ++kj)
	    std::cout << bspl->second->knotval(XFIXED, kvec1[kj]) << " ";
	  std::cout << std::endl;
	  vector<int> kvec2 = bspl->second->kvec(YFIXED);
	  for (size_t kj=0; kj<kvec2.size(); ++kj)
	    std::cout << bspl->second->knotval(YFIXED, kvec2[kj]) << " ";
	  std::cout << std::endl << std::endl;
	  
	  // Compute coefficient using projection
	  Point coef(1);
	  coef[0] = 1.0;

	  if (blevel > 0)
	    {
	      // Update coefficient with respect to lower nesting level coefficients
	      bspl->second->adaptProjCoef(coef);
	      int stop_break = 1;
	    }

	  Point coef2 = bspl->second->Coef();
	  if (blevel > 0)
	    std::cout << "level= " << blevel << ", dist= " << coef.dist(coef2) << std::endl;
	  double gamma = bspl->second->gamma();
	  //lr_surf->setCoef(coef, bspl->second.get());
	}
    }
}

