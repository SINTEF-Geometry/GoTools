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

#include <vcl.h>
#endif


#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/igeslib/IGESconverter.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
using std::vector;
using std::ofstream;
using namespace Go;




int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  // Test number of input arguments
  if (argc != 6)
    {
      std::cout << "Input arguments : Input file(g2), degree, open/closed(1/0), give parameterization (0/1), Give type (0/1)" << std::endl;
      exit(-1);
    }


  double gap = 0.001;
    double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.1;
  double approx = 0.001;


  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");
  int degree = atoi(argv[2]);
  int open = atoi(argv[3]);
  int give_par = atoi(argv[4]);
  int give_type = atoi(argv[5]);

  IGESconverter conv;
  conv.readgo(file1);
  // Read curves and parameter values
  vector<shared_ptr<GeomObject> > gogeom = conv.getGoGeom();
  vector<shared_ptr<SplineCurve> > crvs;
  size_t ki;
  for (ki=0; ki<gogeom.size(); ++ki)
    {
      if (gogeom[ki]->instanceType() == Class_SplineCurve)
	{
	  shared_ptr<GeomObject> lg = gogeom[ki];
	  shared_ptr<SplineCurve> gocv =
	    dynamic_pointer_cast<SplineCurve, GeomObject>(lg);
	  crvs.push_back(gocv);	  
	}
    }

  vector<int> type(crvs.size());
  if (give_type)
    {
      int tp;
      for (int ki=0; ki < (int)crvs.size(); ++ki)
	{
	  std::cout << "Type of curve nr " << ki << ": " << std::endl;
	  std::cin >> tp;
	  type[ki] = tp;
	}
    }

  vector<double> param;
  if (give_par)
    {
      double par;
      for (int ki=0; ki < (int)crvs.size(); ++ki)
	{
	  if (give_par && type[ki] != 1)
	    continue;
	  std::cout << "Parameter value for curve nr " << ki << ": " << std::endl;
	  std::cin >> par;
	  param.push_back(par);
	}
      if (open <= 0)
	{
	  std::cout << "Parameter value for last curve:" << std::endl;
	  std::cin >> par;
	  param.push_back(par);
	}
    }

  CompositeModelFactory factory(approx, gap, neighbour, kink, bend);
  SurfaceModel *model1;
  if (give_type)
    {
      if (give_par)
	model1 = factory.interpolateCurves2(crvs, type, param, open, degree);
      else
	model1 = factory.interpolateCurves(crvs, type, param, open, degree);
    }
  else
    {
      if (give_par)
	model1 = factory.interpolateCurves2(crvs, param, open, degree);
      else
	model1 = factory.interpolateCurves(crvs, param, open, degree);
    }

  ofstream out_file("interpolate_surf.g2");
  if (model1)
    {
      int nmb = model1->nmbEntities();
      for (int kj=0; kj<nmb; kj++)
      {
	  shared_ptr<ParamSurface> surf = model1->getSurface(kj);

	  // Just to be sure
	  shared_ptr<SplineSurface> splinesf = 
	    dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);

	  SplineSurface* surf2 = 
	    splinesf->subSurface(splinesf->startparam_u(),
				 splinesf->startparam_v(),
				 splinesf->endparam_u(),
				 splinesf->endparam_v());
	  surf2->writeStandardHeader(out_file);
	  surf2->write(out_file);

	  delete surf2;
      }
    }
  delete model1;
}

      
  
