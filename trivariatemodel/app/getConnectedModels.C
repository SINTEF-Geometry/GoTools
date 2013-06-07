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

#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/Utils.h"

using namespace Go;


int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  vector<shared_ptr<ftVolume> > volumes;
  
  int ki;
  while (!file1.eof())
    {
      // Read volume from file
      ObjectHeader head;
      file1 >> head;
      shared_ptr<SplineVolume> vol2;
      vol2 = shared_ptr<SplineVolume>(new SplineVolume());
      vol2->read(file1);

      shared_ptr<ParamVolume> pvol
          = dynamic_pointer_cast<ParamVolume, SplineVolume>(vol2);
      volumes.push_back(shared_ptr<ftVolume>(new ftVolume(pvol)));
      // volumes.push_back(shared_ptr<ftVolume>(new ftVolume(vol2)));

      Utils::eatwhite(file1);
    }
    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.05;

  shared_ptr<VolumeModel> model = 
    shared_ptr<VolumeModel>(new VolumeModel(volumes, gap, neighbour, 
					    kink, bend));
  std::ofstream of1("vol_boundary.g2");

  vector<shared_ptr<VolumeModel> > models2 =
    model->getConnectedModels();

  std::cout << "Number of connected models: " << models2.size() << std::endl;

  
  // Outer boundaries
  int kj;
  int nmb_bd = model->nmbBoundaries();
  std::cout << "Outer boundaries: " << nmb_bd << std::endl;

  for (ki=0; ki<(int)models2.size(); ++ki)
    {
      std::ofstream of2("connected.g2");
      int nmb = models2[ki]->nmbEntities();
      for (kj=0; kj<nmb; kj++)
      {
	  shared_ptr<SplineVolume> vol = models2[ki]->getSplineVolume(kj);

	  vol->writeStandardHeader(of2);
	  vol->write(of2);
      }
    }
}
