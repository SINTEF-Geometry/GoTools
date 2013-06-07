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
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : File 1 on g2 format, File 2 on g2 format" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ifstream file2(argv[2]);
  ALWAYS_ERROR_IF(file2.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model1 = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel1 = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model1);

  shared_ptr<CompositeModel> model2 = shared_ptr<CompositeModel>(factory.createFromG2(file2));
  shared_ptr<SurfaceModel> sfmodel2 = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model2);

  if (sfmodel1.get() && sfmodel2.get())
    {
      size_t ki;
      vector<shared_ptr<SurfaceModel> > res = 
	sfmodel1->splitSurfaceModels(sfmodel2);
      vector<shared_ptr<SurfaceModel> > res2;
      for (ki=0; ki<res.size(); ++ki)
	{
	  vector<shared_ptr<SurfaceModel> >res3 = res[ki]->getConnectedModels();
	  res2.insert(res2.end(), res3.begin(), res3.end());
	}

      std::cout << "Number of models: " << res2.size() << std::endl;
      for (ki=0; ki<res2.size(); ++ki)
	{
	  std::cout << "File nr " << ki << ": " << std::endl;
	  char outfile[80];
	  std::cin >> outfile;
	  std::ofstream of1(outfile);

	  int nmb = res2[ki]->nmbEntities();
	  for (int kj=0; kj<nmb; kj++)
	    {
	      shared_ptr<ParamSurface> surf = res2[ki]->getSurface(kj);
	      
	      surf->writeStandardHeader(of1);
	      surf->write(of1);
	    }
	}
    }
}

