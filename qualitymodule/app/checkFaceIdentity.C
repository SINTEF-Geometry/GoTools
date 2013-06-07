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
using std::pair;
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

  vector<pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > >  coinc_faces;
  vector<pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > >  embedded_faces;
  quality->identicalOrEmbeddedFaces(coinc_faces, embedded_faces);

  std::cout << "Number of pairs of identical faces: " << coinc_faces.size() << std::endl;
  std::cout << "Number of pairs of embedded faces: " << embedded_faces.size() << std::endl;

  std::ofstream out_file("identical_sfs.g2");
  size_t ki;
  for (ki=0; ki<coinc_faces.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = coinc_faces[ki].first->surface();
      surf1->writeStandardHeader(out_file);
      surf1->write(out_file);
      shared_ptr<ParamSurface> surf2 = coinc_faces[ki].second->surface();
      surf2->writeStandardHeader(out_file);
      surf2->write(out_file);
  }
  std::ofstream out_file2("embedded_sfs.g2");
  for (ki=0; ki<embedded_faces.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = embedded_faces[ki].first->surface();
      surf1->writeStandardHeader(out_file2);
      surf1->write(out_file2);
      shared_ptr<ParamSurface> surf2 = embedded_faces[ki].second->surface();
      surf2->writeStandardHeader(out_file2);
      surf2->write(out_file2);
  }
  if (do_repair)
    {
      shared_ptr<FaceSetRepair> repair = 
	shared_ptr<FaceSetRepair>(new FaceSetRepair(quality));
 
      repair->identicalAndEmbeddedFaces();

    }

  vector<pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > >  coinc_faces2;
  vector<pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > >  embedded_faces2;
  quality->identicalOrEmbeddedFaces(coinc_faces2, embedded_faces2);


  std::ofstream out_filen("identical_sfs2.g2");
  for (ki=0; ki<coinc_faces2.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = coinc_faces2[ki].first->surface();
      surf1->writeStandardHeader(out_filen);
      surf1->write(out_filen);
      shared_ptr<ParamSurface> surf2 = coinc_faces2[ki].second->surface();
      surf2->writeStandardHeader(out_filen);
      surf2->write(out_filen);
  }
  std::ofstream out_filen2("embedded_sfs2.g2");
  for (ki=0; ki<embedded_faces2.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = embedded_faces2[ki].first->surface();
      surf1->writeStandardHeader(out_filen2);
      surf1->write(out_filen2);
      shared_ptr<ParamSurface> surf2 = embedded_faces2[ki].second->surface();
      surf2->writeStandardHeader(out_filen2);
      surf2->write(out_filen2);
  }

}
