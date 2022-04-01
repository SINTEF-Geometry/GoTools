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
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 3 && argc != 4) {
    std::cout << "Usage: input surface(.g2), output surface(.g2), distance (optional)" << std::endl;
    return -1;
  }

  std::ifstream input(argv[1]);
  std::ofstream output(argv[2]);
  double dist = 0.0;
  if (argc == 4)
    dist = atof(argv[3]);

  GoTools::init();
  Registrator<LRSplineSurface> r293;

  ObjectHeader header;
  header.read(input);
   shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(input);
  
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cerr << "Input file contains no surface" << std::endl;
      exit(-1);
    }

  shared_ptr<ParamSurface> sf2 = sf;  
  shared_ptr<BoundedSurface> bdsf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
  if (bdsf.get())
    sf2 = bdsf->underlyingSurface();

  shared_ptr<LRSplineSurface> lrsf = 
    dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sf2);
  if (!lrsf.get())
    {
      std::cerr << "Input file contains no LR B-spline surface" << std::endl;
      exit(-1);
    }

  if (dist == 0.0)
    {
      // Bounding box
      BoundingBox box = sf->boundingBox();
      Point low = box.low();
      Point high = box.high();
      Point mid = 0.5*(low + high);
      
      // Translate surface
      lrsf->translate(-mid);
    }
  else
    {
      Point vec(sf->dimension());
      vec.setValue(dist);
      lrsf->translate(vec);
    }
  
  sf->writeStandardHeader(output);
  sf->write(output);
}

