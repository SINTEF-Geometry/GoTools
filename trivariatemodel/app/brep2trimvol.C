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
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/CreateTrimVolume.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 4)
    {
      cout << "Usage: " << "<infile> <file type in (1/2)> <outfile>" << endl;
      exit(-1);
    }

  std::string infile(argv[1]);

  int type_in = atoi(argv[2]);
  ofstream outfile(argv[3]);

  double gap, neighbour, kink;
  shared_ptr<SurfaceModel> sfmodel;
  int material_id = -1;
  if (type_in == 2)
    {
      VolumeModelFileHandler fileread;
      shared_ptr<Body> body = fileread.readBody(infile.c_str());
      tpTolerances top = body->getTolerances();
      gap = top.gap;
      neighbour = top.neighbour;
      kink = top.kink;

      sfmodel = body->getOuterShell();
      material_id = body->getMaterial();
    }
  else
    {
      // The tolerances must be set according to the properties of the model.
      // The neighbour tolerance must be smaller than the smallest entity in the
      // model, but larger than the largest gap.
      // The gap tolerance must be smaller than the neighbour tolerance
      double gap = 0.001; //0.0001;
      double neighbour = 0.01; //0.001;
      double kink = 0.01;
      double approxtol = 0.001;

      CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

      std::ifstream is(infile);
      CompositeModel *model = factory.createFromG2(is);

      sfmodel = 
	shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
      if (!sfmodel.get())
	{
	  std::cout << "No input model read" << std::endl;
	  exit(-1);
	}
 
      if (sfmodel->nmbBoundaries() > 0)
	{
	  std::cout << "Not a brep solid. Consider increasing the neighbour tolerance" << std::endl;
	  exit(-1);
	}
      
      bool isOK = sfmodel->checkShellTopology();
      std::cout << "Shell topology: " << isOK << std::endl;
    }
  CreateTrimVolume trim(sfmodel, material_id);

  shared_ptr<ftVolume> vol = trim.fetchOneTrimVol();

  VolumeModelFileHandler filehandler;
  filehandler.writeStart(outfile);
  filehandler.writeHeader("Trimmed volume", outfile);
  filehandler.writeVolume(vol, outfile);
  filehandler.writeEnd(outfile);
  int stop_break = 1;
}
