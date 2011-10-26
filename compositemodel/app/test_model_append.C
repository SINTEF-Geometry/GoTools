//===========================================================================
//
// File : test_append.C
//
// Created: Thu Apr  3 10:42:30 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_model_append.C,v 1.2 2009-05-13 07:29:54 vsk Exp $
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
  if (argc != 12)
    {
      std::cout << "Input arguments : Input file on IGES format, Input file on IGES format, ";
      std::cout << "x-value of first point on plane, y-value f first point, ... , z-value of third point" << std::endl;
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
  std::ifstream file2(argv[2]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  CompositeModelFactory factory(approx, gap, neighbour, kink, bend);

  std::shared_ptr<SurfaceModel> sm1, sm2;
  sm1.reset((SurfaceModel*) factory.createFromIges(file1));
  sm2.reset((SurfaceModel*) factory.createFromIges(file2));

  sm1 -> append(sm2);

  Point
    p_x (atof(argv[3]), atof(argv[4]), atof(argv[5])),
    p_y (atof(argv[6]), atof(argv[7]), atof(argv[8])),
    p_z (atof(argv[9]), atof(argv[10]), atof(argv[11]));

  ftPlane pl(p_x, p_y, p_z);

  ftCurve cv = sm1 -> intersect(pl);

  ofstream osCrv ("dumpCurve.g2");
  ofstream osPl ("dumpPlane.g2");

  cv.writeSpaceCurve(osCrv);

  osPl << "200 1 0 0" << std::endl;
  osPl << "3 0" << std:: endl;
  osPl << "2 2" << std::endl;
  osPl << "0 0 1 1" << std::endl;
  osPl << "2 2" << std::endl;
  osPl << "0 0 1 1" << std::endl;
  osPl << p_x << std::endl;
  osPl << p_x << std::endl;
  osPl << p_y << std::endl;
  osPl << p_z << std::endl;

}
