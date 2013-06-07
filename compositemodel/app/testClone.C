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
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
using std::cout;
using std::endl;
using namespace Go;




void dump(int modelNumber, SurfaceModel* sm)
{
  cout << "\n\n";
  cout << "Dumping model " << modelNumber << endl;

  printf("Model = %p\n",sm);

  cout << "#Faces = " << sm->nmbEntities() << endl;

  int top = sm->nmbEntities();

  for (int i=0; i<top; ++i)
    {
      shared_ptr<ftSurface> f;
      shared_ptr<ParamSurface> f2;
      f = sm -> getFace(i);
      f2 = sm -> getSurface(i);
      printf("Face %i = %p, param = %p\n", i, f.get(), f2.get());
    }
}


void dump2(int modelNumber, CompositeCurve* cv)
{
  cout << "\n\n";
  cout << "Dumping model " << modelNumber << endl;

  printf("Model = %p\n",cv);

  cout << "#Curves = " << cv->nmbEntities() << endl;

  int top = cv->nmbEntities();

  for (int i=0; i<top; ++i)
    {
      shared_ptr<ParamCurve> f2;
      f2 = cv -> getCurve(i);
      printf("Curve %i = %p\n", i, f2.get());
    }

  double start, end;
  cv->parameterRange(start, end);
  std::cout << "Parameter range: " << start << ", " << end << std::endl;
  
  double par = 0.2*start + 0.8*end;
  Point pt;
  cv->evaluateCurve(par, pt);
  std::cout << "Par: " << par << " Pos: " << pt[0] << " " << pt[1];
  std::cout << " " << pt[2] << std::endl;
}










int main( int argc, char* argv[] )
{
  // Test number of input arguments
  if (argc != 2)
    {
      std::cout << "Input arguments : Input file on IGES format" << std::endl;
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


  CompositeModelFactory factory(approx, gap, neighbour, kink, bend);
  CompositeModel* cm = factory.createFromIges(file1);

  CompositeModel *cm2 = cm->clone();
  SurfaceModel* sf1 = dynamic_cast<SurfaceModel*>(cm);
  SurfaceModel* sf2 = dynamic_cast<SurfaceModel*>(cm2);

  CompositeCurve *cv1 = dynamic_cast<CompositeCurve*>(cm);
  CompositeCurve *cv2 = dynamic_cast<CompositeCurve*>(cm2);
  if (sf1)
    {
      dump(1,sf1);
      dump(2,sf2);
    }
  else if (cv1)
    {
      dump2(1,cv1);
      dump2(2,cv2);
    }

}

















