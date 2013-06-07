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

#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include <fstream>

using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 9) {
    std::cout << "Input parameters : Input file on g2 format, output file, point in plane, normal" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ofstream file2(argv[2]);

  std::ofstream file3("loftvol.g2");

  Point pnt(atof(argv[3]), atof(argv[4]), atof(argv[5]));
  Point norm(atof(argv[6]), atof(argv[7]), atof(argv[8]));

  double gap = 0.0001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  int nmb_sfs = sfmodel->nmbEntities();
  vector<shared_ptr<ParamSurface> > mirrored;
  vector<shared_ptr<ftVolume> > blocks;
  for (int ki=0; ki<nmb_sfs; ++ki)
    {
      shared_ptr<ParamSurface> curr = sfmodel->getSurface(ki);
      shared_ptr<ParamSurface> mod = 
	shared_ptr<ParamSurface>(curr->mirrorSurface(pnt, norm));
      mirrored.push_back(mod);

      vector<shared_ptr<SplineSurface> > sfs(2);
      sfs[0] = dynamic_pointer_cast<SplineSurface,ParamSurface>(curr);
      sfs[1] = dynamic_pointer_cast<SplineSurface,ParamSurface>(mod);

      if (sfs[0].get() && sfs[1].get())
	{
	  shared_ptr<ParamVolume> vol = 
	    shared_ptr<ParamVolume>(LoftVolumeCreator::loftVolume(sfs.begin(), 2));

	  vol->writeStandardHeader(file3);
	  vol->write(file3);

	  shared_ptr<ftVolume> ftvol = 
	    shared_ptr<ftVolume>(new ftVolume(vol, gap, kink));
	  blocks.push_back(ftvol);
	}
    }

  shared_ptr<VolumeModel> volmodel = 
    shared_ptr<VolumeModel>(new VolumeModel(blocks, gap, neighbour, 
					    kink, 10.0*kink));


  int nmb_vol = volmodel->nmbEntities();
  for (int ki=0; ki<nmb_vol; ++ki)
    {
      shared_ptr<ParamVolume> vol = volmodel->getVolume(ki);
       vol->writeStandardHeader(file2);
       vol->write(file2);
     }
}
  
