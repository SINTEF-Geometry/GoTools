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

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input surface file(.g2) "  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream infile(argv[1]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;
  //Registrator<BoundedSurface> r210;

  // Read input surface
  ObjectHeader header;
  try {
    header.read(infile);
  }
  catch (...)
    {
      std::cerr << "Exiting" << std::endl;
      exit(-1);
    }
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(infile);
  
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cerr << "Input file contains no surface" << std::endl;
      exit(-1);
    }

  shared_ptr<BoundedSurface> bdsf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
  if (bdsf.get())
    sf = bdsf->underlyingSurface();

  // Fetch domain
  RectDomain domain = sf->containingDomain();

  // Fetch box
  BoundingBox bb = sf->boundingBox();
  

  // Output
  int dim = sf->dimension();
  printf("Surface domain: %13.4f %13.4f %13.4f %13.4f \n", domain.umin(), 
	 domain.umax(), domain.vmin(), domain.vmax());
  printf("Bound z: %13.7f %13.7f \n", bb.low()[dim-1], bb.high()[dim-1]);

  shared_ptr<LRSplineSurface> lrsf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sf);
  if (lrsf.get())
    {
      std::ofstream ofel("sf_elements.g2");
      LineCloud lines = lrsf->getElementBds();
      lines.writeStandardHeader(ofel);
      lines.write(ofel);
    }

  if (bdsf.get())
    {
      std::ofstream ofl("sf_loop.g2");
      SplineDebugUtils::writeBoundary(*bdsf, ofl);
    }
}

