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

#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/utils/config.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: input surface(.g2), output surface(.g2), translate (0/1)" << std::endl;
    return -1;
  }

  std::ifstream input(argv[1]);
  std::ofstream output(argv[2]);
  int translate = atoi(argv[3]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;
  //Registrator<BoundedSurface> r210;

  ObjectHeader header;
  header.read(input);
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(input);
  
  shared_ptr<ParamSurface> sf = dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);

  if (sf->dimension() != 1)
    {
      std::cout << "Dimension not equal to 1" << std::endl;
      return 0;
    }

  // Translate parameter domain
  if (translate)
    {
      RectDomain dom = sf->containingDomain();
      double umin = dom.umin();
      double umax = dom.umax();
      double vmin = dom.vmin();
      double vmax = dom.vmax();
      double umid = 0.5*(umin+umax);
      double vmid = 0.5*(vmin+vmax);
  
      sf->setParameterDomain(umin - umid, umax - umid,
			     vmin - vmid, vmax - vmid);
    }

  shared_ptr<LRSplineSurface> tmp_lr;
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
  if (bd_sf.get())
    {
      shared_ptr<ParamSurface> tmp_sf = bd_sf->underlyingSurface();
      tmp_lr = 
	dynamic_pointer_cast<LRSplineSurface, ParamSurface>(tmp_sf);
    }
  else
    tmp_lr = 
	dynamic_pointer_cast<LRSplineSurface, ParamSurface>(sf);
    
  tmp_lr->to3D();

  sf->writeStandardHeader(output);
  sf->write(output);
}

