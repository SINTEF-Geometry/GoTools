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

#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 5 && argc != 6) {
      std::cout << "Input parameters : centre (x,y,z), radius (turn?)" << std::endl;
    exit(-1);
  }

#ifdef __BORLANDC__
  using Go::Point;
#endif

  Point centre(atof(argv[1]), atof(argv[2]), atof(argv[3]));
  double radius = atof(argv[4]);
    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;
  bool turn = false;
  if (argc == 6)
    turn = (atoi(argv[5]) ? true : false);

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromSphere(centre, radius);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

  if (sfmodel)
  {
      std::ofstream out_file("sphere.g2");
      int nmb = sfmodel->nmbEntities();
      for (int ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

	  if (turn)
	    surf->swapParameterDirection();

	  surf->writeStandardHeader(out_file);
	  surf->write(out_file);
      }
  }

  delete model;
}
