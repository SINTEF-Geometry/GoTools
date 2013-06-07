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
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format, " << std::endl;
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

  vector<shared_ptr<ParamCurve> >  g1_crvs;
  vector<shared_ptr<ParamCurve> >  c1_crvs;
  quality.cvC1G1Discontinuity(c1_crvs, g1_crvs);

  std::cout << "Number c1 discontinuity crvs: " << c1_crvs.size() << std::endl;
  std::cout << "Number g1 discontinuity crvs: " << g1_crvs.size() << std::endl;

  std::ofstream out_file("c1dis_cvs.g2");
  size_t ki;
  for (ki=0; ki<c1_crvs.size(); ki++)
  {
      c1_crvs[ki]->writeStandardHeader(out_file);
      c1_crvs[ki]->write(out_file);
   }

  std::ofstream out_file2("g1dis_cvs.g2");
  for (ki=0; ki<g1_crvs.size(); ki++)
  {
      g1_crvs[ki]->writeStandardHeader(out_file2);
      g1_crvs[ki]->write(out_file2);
   }

  vector<shared_ptr<ParamCurve> >  g1_crvs2;
  vector<shared_ptr<ParamCurve> >  c1_crvs2;
  quality.cvC1G1Discontinuity(c1_crvs2, g1_crvs2);


  std::ofstream out_filen("c1dis_cvs2.g2");
  for (ki=0; ki<c1_crvs2.size(); ki++)
  {
      c1_crvs2[ki]->writeStandardHeader(out_filen);
      c1_crvs2[ki]->write(out_filen);
   }

  std::ofstream out_filen2("g1dis_cvs2.g2");
  for (ki=0; ki<g1_crvs2.size(); ki++)
  {
      g1_crvs2[ki]->writeStandardHeader(out_filen2);
      g1_crvs2[ki]->write(out_filen2);
   }

}
