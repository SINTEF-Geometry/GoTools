//===========================================================================
//
// File : test_evaluate.C
//
// Created: Thu Apr  3 10:35:14 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_evaluate.C,v 1.2 2009-05-13 07:29:53 vsk Exp $
//
// Description:
//
//===========================================================================


#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
using std::shared_ptr;
using namespace std;
using namespace Go;




int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  // Test number of input arguments
  if ((argc%3) != 2)
    {
      std::cout << "Input arguments : Input file on IGES format,";
      std::cout << " surface index_1, param_u_1, param_v_1, ..., param_v_N" << std::endl;
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
  SurfaceModel* sm = (SurfaceModel*) factory.createFromIges(file1);

  std::ofstream pointstr("points.g2");

  for (int i = 2; i < argc; )
    {
      int surf_ind = atoi(argv[i++]);
      double param[2];
      param[0] = atof(argv[i++]);
      param[1] = atof(argv[i++]);
      Point p;
      sm -> evaluate(surf_ind, param, p);
      
      pointstr << "400 1 0 4 255 0 0 255" << endl;
      pointstr << "1" << endl;
      pointstr << p << endl;
    }

}
