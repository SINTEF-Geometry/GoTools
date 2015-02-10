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
    std::cout << "Input parameters : Input file on g2 format, parameter tolerance" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.0001;
  double neighbour = 0.001;
  double kink = 0.01;
  double approxtol = 0.01;

  double tol = atof(argv[2]);

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<shared_ptr<ParamCurve> > cv_knots;
  vector<shared_ptr<ParamSurface> > sf_knots;
  quality.indistinctKnots(cv_knots, sf_knots, tol);

  std::cout << "Trimming curves with indistinct knots: " << cv_knots.size() << std::endl;
  std::cout << "Surfaces with indistinct knots: " << sf_knots.size() << std::endl;

  std::ofstream out_file("indistinct_knots.g2");
  size_t ki;
  for (ki=0; ki<cv_knots.size(); ki++)
  {
      cv_knots[ki]->writeStandardHeader(out_file);
      cv_knots[ki]->write(out_file);
  }

  for (ki=0; ki<sf_knots.size(); ki++)
  {
      sf_knots[ki]->writeStandardHeader(out_file);
      sf_knots[ki]->write(out_file);
  }

  vector<shared_ptr<ParamCurve> > cv_knots2;
  vector<shared_ptr<ParamSurface> > sf_knots2;
  quality.indistinctKnots(cv_knots2, sf_knots2, tol);

  std::ofstream out_file2("indistinct_knots2.g2");
  for (ki=0; ki<cv_knots2.size(); ki++)
  {
      cv_knots2[ki]->writeStandardHeader(out_file2);
      cv_knots2[ki]->write(out_file2);
  }

  for (ki=0; ki<sf_knots2.size(); ki++)
  {
      sf_knots2[ki]->writeStandardHeader(out_file2);
      sf_knots2[ki]->write(out_file2);
  }

}
