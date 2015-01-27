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
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"

#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: lrspline_in (.g2) refinement_in lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ifstream filein2(argv[2]);
  std::ofstream fileout(argv[3]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read input surface
  ObjectHeader header;
  try {
    header.read(filein);
  }
  catch (...)
    {
      std::cerr << "Exiting" << std::endl;
      exit(-1);
    }
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(filein);
  
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cerr << "Input file contains no surface" << std::endl;
      exit(-1);
    }

  shared_ptr<LRSplineSurface> lrsf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sf);
  if (!lrsf.get())
    {
      shared_ptr<SplineSurface> splsf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(sf);
      if (splsf.get())
	lrsf = shared_ptr<LRSplineSurface>(new LRSplineSurface(splsf.get(), 1.0e-6));
    }
  if (!lrsf.get())
    {
      std::cerr << "Input file contains no spline surface" << std::endl;
      exit(-1);
    }
    
  
  shared_ptr<LRSplineSurface> tmp2(lrsf->clone());
  if (tmp2->dimension() == 1)
    tmp2->to3D();

  // tmp2->writeStandardHeader(fileout);
  // tmp2->write(fileout);
  // fileout << std::endl;
  // LineCloud lines2 = tmp2->getElementBds();
  // lines2.writeStandardHeader(fileout);
  // lines2.write(fileout);
  
  int nmb_refs;
  filein2 >> nmb_refs;
  for (int ki=0; ki<nmb_refs; ++ki)
    {
      double parval, start, end;
      int dir;
      int mult;

      filein2 >> parval;
      filein2 >> start;
      filein2 >> end;
      filein2 >> dir;
      filein2 >> mult;
      //lrsf->refine((dir==0) ? XFIXED : YFIXED, parval, start, end, mult);
      std::cout << "Iteration no. " << ki << std::endl;
      lrsf->refine((dir==0) ? XFIXED : YFIXED, parval, start, end, mult, true);

      puts("Writing lr-spline to file.");
      if (lrsf->dimension() == 1)
	lrsf->to3D();
      lrsf->writeStandardHeader(fileout);
      lrsf->write(fileout);
      fileout << std::endl;
    }
  return 0;
}
