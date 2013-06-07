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

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/qualitymodule/FaceSetRepair.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, repair?" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  int do_repair = atoi(argv[2]);

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  shared_ptr<FaceSetQuality> quality = 
    shared_ptr<FaceSetQuality>(new FaceSetQuality(gap, kink, approxtol));
  quality->attach(sfmodel);

  vector<shared_ptr<Loop> > inconsistent_loops;
  quality->loopOrientationConsistency(inconsistent_loops);

  std::cout << "Number of inconsistent loops: " << inconsistent_loops.size() << std::endl;

  std::ofstream out_file0("inconsistent_loops.g2");
  size_t ki, kj;
  for (ki=0; ki<inconsistent_loops.size(); ki++)
  {
      for (kj=0; kj<inconsistent_loops[ki]->size(); ++kj)
      {
	  shared_ptr<ParamCurve> crv = inconsistent_loops[ki]->getEdge(kj)->geomEdge()->geomCurve();
	  SplineCurve *crv2 = crv->geometryCurve();
	  if (crv2)
	  {
	      crv2->writeStandardHeader(out_file0);
	      crv2->write(out_file0);
	  }
      }
  }

  vector<shared_ptr<Loop> > inconsistent_loops2;
  quality->loopOrientationConsistency(inconsistent_loops2);

  std::ofstream out_filen0("inconsistent_loops2.g2");
  for (ki=0; ki<inconsistent_loops2.size(); ki++)
  {
      for (kj=0; kj<inconsistent_loops2[ki]->size(); ++kj)
      {
	  shared_ptr<ParamCurve> crv = inconsistent_loops2[ki]->getEdge(kj)->geomEdge()->geomCurve();
	  SplineCurve *crv2 = crv->geometryCurve();
	  if (crv2)
	  {
	      crv2->writeStandardHeader(out_filen0);
	      crv2->write(out_filen0);
	  }
      }
  }

  vector<shared_ptr<ftSurface> >   inconsistent_faces;
  quality->faceNormalConsistency(inconsistent_faces);

  std::cout << "Number of inconsistent faces: " << inconsistent_faces.size() << std::endl;

  std::ofstream out_file("inconsistent_sfs.g2");
  for (ki=0; ki<inconsistent_faces.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = inconsistent_faces[ki]->surface();
      surf1->writeStandardHeader(out_file);
      surf1->write(out_file);
  }

  if (do_repair)
    {
      shared_ptr<FaceSetRepair> repair = 
	shared_ptr<FaceSetRepair>(new FaceSetRepair(quality));

      repair->consistentFaceNormal();
    }

  vector<shared_ptr<ftSurface> >   inconsistent_faces2;
  quality->faceNormalConsistency(inconsistent_faces2);

  std::ofstream out_filen("inconsistent_sfs2.g2");
  for (ki=0; ki<inconsistent_faces2.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = inconsistent_faces2[ki]->surface();
      surf1->writeStandardHeader(out_filen);
      surf1->write(out_filen);
  }

}
