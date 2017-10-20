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
#include "GoTools/compositemodel/Body.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format," << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.0001;
  double neighbour = 0.001;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

  if (sfmodel)
  {
    BoundingBox box = sfmodel->boundingBox();
    std::cout << "Bounding box, low: " << box.low() << std::endl;
    std::cout << "high: " << box.high() << std::endl;
      int ki;
      std::ofstream out_file("fileg2.g2");
      int nmb = sfmodel->nmbEntities();
      vector<shared_ptr<ftSurface> > faces;
      for (ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

	  surf->writeStandardHeader(out_file);
	  surf->write(out_file);

	  shared_ptr<ftSurface> curr_face = sfmodel->getFace(ki);
	  faces.push_back(curr_face);
      }

      // shared_ptr<SurfaceModel> model2 = 
      // 	shared_ptr<SurfaceModel>(new SurfaceModel(approxtol, gap,
      // 						  neighbour, kink,
      // 						  10.0*kink, faces, true));
      vector<shared_ptr<SurfaceModel> > model2 = sfmodel->getConnectedModels();

      std::ofstream out_file2("bd.g2");
      int nmb_bd = sfmodel->nmbBoundaries();
      std::cout << "Number of boundaries: " << nmb_bd << std::endl;
      for (ki=0; ki<nmb_bd; ki++)
      {
	  ftCurve bdcv = sfmodel->getBoundary(ki);
	  bdcv.writeSpaceCurve(out_file2);
      }

      // Create Body
      shared_ptr<Body> body(new Body(model2, 1));
      std::ofstream out_file3("model.g22");
      CompositeModelFileHandler filehandler;
      filehandler.writeStart(out_file3);
      filehandler.writeHeader("Translated from g2", out_file3);
      filehandler.writeBody(body, out_file3);
      filehandler.writeEnd(out_file3);
  }

  delete model;
}
