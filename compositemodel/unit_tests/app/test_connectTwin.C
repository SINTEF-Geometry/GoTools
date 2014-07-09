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
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/intersections/Identity.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
      std::cout << "Surface file1(g2), surface file2, tolerance" << std::endl;
    exit(-1);
  }

  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ifstream file2(argv[2]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double tol = atof(argv[3]);

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);
  CompositeModel *model1 = factory.createFromG2(file1);
  CompositeModel *model2 = factory.createFromG2(file2);
  SurfaceModel *sfmodel1 = dynamic_cast<SurfaceModel*>(model1);
  SurfaceModel *sfmodel2 = dynamic_cast<SurfaceModel*>(model2);

  int nmb1 = sfmodel1->nmbEntities();
  int nmb2 = sfmodel2->nmbEntities();
  int ki, kj;

  Identity ident;
  for (ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ftSurface> face1 = sfmodel1->getFace(ki);
      for (kj=0; kj<nmb2; ++kj)
	{
	  shared_ptr<ftSurface> face2 = sfmodel2->getFace(kj);

	  // Check coincidence
	  
	  int coincide = ident.identicalSfs(face1->surface(),
					    face2->surface(),
					    tol);

	  if (coincide)
	    {
	      bool connected = face1->connectTwin(face2.get(), tol);
	      
	      std::cout << "Connected: " << connected << std::endl;
	    }
	}
    }
}
