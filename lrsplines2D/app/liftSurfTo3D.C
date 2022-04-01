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

/// Lift a LRSplineSurface or a trimmed LRSplineSurface from 1D to 3D.

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Factory.h"

#include <iostream>
#include <fstream>

using namespace Go;


int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: lr_spline_sf_1D.g2 lr_spline_sf_3D.g2" << std::endl;
    return -1;
  }

  GoTools::init();
  Registrator<LRSplineSurface> r293;

  std::ifstream filein(argv[1]); // Input lr spline sf (bivariate 1D).
  std::ofstream fileout(argv[2]); // Output lr spline sf (bivariate 3D).

  ObjectHeader header;
  header.read(filein);
  if (header.classType() == Class_LRSplineSurface)
  {
      LRSplineSurface lrspline_sf;
      lrspline_sf.read(filein);
      if (lrspline_sf.dimension() != 1)
      {
          std::cout << "Expecting input to be (bivariate) 1D, not " << lrspline_sf.dimension() << "D." << std::endl;
          return 1;
      }

      lrspline_sf.to3D();

      lrspline_sf.writeStandardHeader(fileout);
      lrspline_sf.write(fileout);
  }
  else if (header.classType() == Class_BoundedSurface)
  {
      BoundedSurface bd_sf;
      bd_sf.read(filein);
      if (bd_sf.dimension() != 1)
      {
          std::cout << "Expecting input to be (bivariate) 1D, not " << bd_sf.dimension() << "D." << std::endl;
          return 1;
      }

      shared_ptr<ParamSurface> under_sf = bd_sf.underlyingSurface();
      if (under_sf->instanceType() == Class_LRSplineSurface)
      {
	  shared_ptr<LRSplineSurface> lr_sf = dynamic_pointer_cast<LRSplineSurface>(under_sf);
	  lr_sf->to3D();
	  MESSAGE("Missing support for geometric trim curves!");
	  bd_sf.writeStandardHeader(fileout);
	  bd_sf.write(fileout);
      }
      else
      {
	  MESSAGE("Surface type " << header.classType() << " in the underlying surface is not supported.");
      }
  }
  else if (header.classType() == Class_SplineSurface)
    {
      shared_ptr<SplineSurface> sf(new SplineSurface());
      sf->read(filein);
      shared_ptr<LRSplineSurface> lr_sf(new LRSplineSurface(sf.get(), 1.0e-8));
      lr_sf->to3D();
      shared_ptr<SplineSurface> splsf(lr_sf->asSplineSurface());
      splsf->writeStandardHeader(fileout);
      splsf->write(fileout);	  
    }
  else
  {
      MESSAGE("Surface type " << header.classType() << " is not supported.");
  }
}
