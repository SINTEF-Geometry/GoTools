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
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Factory.h"

#include <iostream>
#include <fstream>
#include <string.h>




using namespace Go;

// The app moves a LRSplineSurface to the origin.
// We also make sure the parameter domain is moved.
int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: lrspline_in.g2 lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);

  GoTools::init();
  Registrator<LRSplineSurface> r293;

  ObjectHeader header;
  header.read(filein);
  // shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
  // shared_ptr<GeomObject> obj2(Factory::createObject(Class_LRSplineSurface));
  if (header.classType() != Class_BoundedSurface)
  {
      std::cout << "Unsupported surface type!" << std::endl;
      return 1;
  }
  else
  {
      BoundedSurface bd_sf;
      bd_sf.read(filein);
      shared_ptr<LRSplineSurface> lrsf = dynamic_pointer_cast<LRSplineSurface>(bd_sf.underlyingSurface());
      assert(lrsf.get() != NULL);

      MESSAGE("Setting parameter domain to the unit square!");
      lrsf->setParameterDomain(0.0, 1.0, 0.0, 1.0);

      const int dim = lrsf->dimension();
      if (dim == 3)
      {
	  MESSAGE("Translating the surface to the origin!");
	  BoundingBox bd_box = lrsf->boundingBox();
	  Point mid_pt = 0.5*(bd_box.low() + bd_box.high());
	  Point transl_pt = -mid_pt;
	  std::cout << "transl_pt: (" << transl_pt[0] << ", " << transl_pt[1] <<
              ", " << transl_pt[2] << ")" << std::endl;
	  lrsf->translate(transl_pt);
      }
      else
      {
	  std::cout << "Dimension is " << dim << "! Only handling 3!" << std::endl;
          return 1;
      }

      bd_sf.writeStandardHeader(fileout);
      bd_sf.write(fileout);
  }

  return 0;
}
