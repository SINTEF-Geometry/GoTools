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
#include "GoTools/lrsplines2D/DefineRefs2D.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/LineCloud.h"
#include <fstream>

using namespace Go;
using std::vector;


int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      std::cout << "Usage: input spline surface, number of derivatives, number of refinements" << std::endl;
      exit(-1);
    }

  std::ifstream input(argv[1]);
  int ider = atoi(argv[2]);
  int numref = atoi(argv[3]);
  
  // Read spline surface from file
   // Read header specifying the type of geometry entity
  ObjectHeader header;
  try {
    header.read(input);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(-1);
    }
  
  shared_ptr<SplineSurface> surf1(new SplineSurface());
  surf1->read(input);
  shared_ptr<LRSplineSurface> surf2(new LRSplineSurface(surf1.get(),
							1.0e-6));

  std::ofstream of("lrsurf.g2");
  surf2->writeStandardHeader(of);
  surf2->write(of);

  double upar = 0.5*(surf2->startparam_u() + surf2->endparam_u());
  double vpar = 0.5*(surf2->startparam_v() + surf2->endparam_v());

  vector<Point> pts1((ider+1)*(ider+2)/2), pts2((ider+1)*(ider+2)/2);
  surf1->point(pts1, upar, vpar, ider);
  surf2->point(pts2, upar, vpar, ider);

  std::cout << "Spline surface: " << std::endl;
  for (size_t ki=0; ki<pts1.size(); ++ki)
    std::cout << pts1[ki] << " ";
  std::cout << std::endl;

  std::cout << "LR surface: " << std::endl;
  for (size_t ki=0; ki<pts2.size(); ++ki)
    std::cout << pts2[ki] << " ";
  std::cout << std::endl;

  int kder1 = ider/2;
  int kder2 = ider - kder1;

  Point p1, p2;
  std::cout << "kder1: " << kder1 << ", kder2: " << kder2 << std::endl;
  p1 = surf2->operator()(upar, vpar, kder1, kder2);
  std::cout << "p1: " << p1 << std::endl;
  p2 = surf2->operator()(upar, vpar, kder2, kder1);
  std::cout << "p2: " << p2 << std::endl;

  double u1 = surf2->startparam_u();
  double u2 = surf2->endparam_u();
  double v1 = surf2->startparam_v();
  double v2 = surf2->endparam_v();

  std::ofstream ofel1("sf_elements.g2");
  int num = 10;
  LineCloud lines = surf2->getElementBds(num);
  lines.writeStandardHeader(ofel1);
  lines.write(ofel1);
  
  std::cout << "Span: " << u1 << " " << u2  << " " << v1  << " " << v2  << std::endl;
  vector<LRSplineSurface::Refinement2D> refs_x, refs_y;
  std::cout << "Specify " << numref << " refinements: u, v, dir" << std::endl;
  for (int ki=0; ki<numref; ++ki)
    {
      double u, v;
      int d;
      std::cin >> u >> v >> d;

      int mult = 1;
      if (d == 1 || d == 3)
	{
	  Direction2D dir = XFIXED;
	  LRSplineSurface::Refinement2D ref;
	  DefineRefs2D::refineFullSpan(*surf2, u, v, dir, mult, ref);
	  refs_x.push_back(ref);
	}
      if (d == 2 || d == 3)
	{
	  Direction2D dir = YFIXED;
	  LRSplineSurface::Refinement2D ref;
	  DefineRefs2D::refineFullSpan(*surf2, u, v, dir, mult, ref);
	  refs_y.push_back(ref);
	}
     }

  for (size_t kr=0; kr<refs_x.size(); ++kr)
    {
      surf2->refine(refs_x[kr], true);
    }

  for (size_t kr=0; kr<refs_y.size(); ++kr)
    {
      surf2->refine(refs_y[kr], true);
    }

  std::ofstream ofel2("sf_elements2.g2");
  LineCloud lines2 = surf2->getElementBds(num);
  lines2.writeStandardHeader(ofel2);
  lines2.write(ofel2);
  
  vector<Point> pts3((ider+1)*(ider+2)/2);
  surf2->point(pts3, upar, vpar, ider);
  
  std::cout << "LR surface refined: " << std::endl;
  for (size_t ki=0; ki<pts3.size(); ++ki)
    std::cout << pts3[ki] << " ";
  std::cout << std::endl;

  Point p3, p4;
  p3 = surf2->operator()(upar, vpar, kder1, kder2);
  std::cout << "p3: " << p3 << std::endl;
  p4 = surf2->operator()(upar, vpar, kder2, kder1);
  std::cout << "p4: " << p4 << std::endl;
}

  
 
