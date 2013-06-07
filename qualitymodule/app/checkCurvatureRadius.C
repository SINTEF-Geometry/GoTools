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
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, minium curvature radius" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  double curvature_radius = atof(argv[2]);
  quality.setCurvatureRadius(curvature_radius);

  vector<pair<shared_ptr<ftPoint>, double> > curv_radius;
  pair<shared_ptr<ftPoint>, double> min_curv_radius;

  quality.sfCurvatureRadius(curv_radius, min_curv_radius);

  std::cout << "Number of small radiuses: " << curv_radius.size() << std::endl;
  std::cout << "Minimum curvature radius: " << min_curv_radius.second << std::endl;

  std::ofstream out_file("curvature_sfs.g2");
  for (size_t ki=0; ki<curv_radius.size(); ki++)
  {
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "1 " << std::endl;
      const Point& pos = curv_radius[ki].first->position();
      out_file << pos[0] << "  ";
      out_file << pos[1] << "  ";
      out_file << pos[2] << "  " << std::endl;
      out_file << std::endl;
      curv_radius[ki].first->face()->surface()->writeStandardHeader(out_file);
      curv_radius[ki].first->face()->surface()->write(out_file);
      std::cout << "Curvature radius: " << curv_radius[ki].second << std::endl;
  }
  out_file << "400 1 0 4 0 255 0 255" << std::endl;
  out_file << "1 " << std::endl;
  const Point& pos2 = min_curv_radius.first->position();
  out_file << pos2[0] << "  ";
  out_file << pos2[1] << "  ";
  out_file << pos2[2] << "  " << std::endl;
  out_file << std::endl;
  min_curv_radius.first->face()->surface()->writeStandardHeader(out_file);
  min_curv_radius.first->face()->surface()->write(out_file);

  vector<pair<shared_ptr<ftPoint>, double> > curv_radius2;
  pair<shared_ptr<ftPoint>, double> min_curv_radius2;

  quality.sfCurvatureRadius(curv_radius2, min_curv_radius2);

  std::cout << "Number of small radiuses: " << curv_radius2.size() << std::endl;
  std::cout << "Minimum curvature radius: " << min_curv_radius2.second << std::endl;

  std::ofstream out_file2("curvature_sfs2.g2");
  for (size_t ki=0; ki<curv_radius2.size(); ki++)
  {
      out_file2 << "400 1 0 4 255 0 0 255" << std::endl;
      out_file2 << "1 " << std::endl;
      const Point& pos = curv_radius2[ki].first->position();
      out_file2 << pos[0] << "  ";
      out_file2 << pos[1] << "  ";
      out_file2 << pos[2] << "  " << std::endl;
      out_file2 << std::endl;
      curv_radius2[ki].first->face()->surface()->writeStandardHeader(out_file2);
      curv_radius2[ki].first->face()->surface()->write(out_file2);
      std::cout << "Curvature radius: " << curv_radius2[ki].second << std::endl;
  }
  out_file2 << "400 1 0 4 0 255 0 255" << std::endl;
  out_file2 << "1 " << std::endl;
  const Point& pos3 = min_curv_radius2.first->position();
  out_file2 << pos3[0] << "  ";
  out_file2 << pos3[1] << "  ";
  out_file2 << pos3[2] << "  " << std::endl;
  out_file2 << std::endl;
  min_curv_radius2.first->face()->surface()->writeStandardHeader(out_file2);
  min_curv_radius2.first->face()->surface()->write(out_file2);

}
