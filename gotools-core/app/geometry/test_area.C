//===========================================================================
//
// File : test_area.C
//
// Created: Thu Jul  3 08:23:37 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_area.C,v 1.1 2008-07-08 10:59:41 kfp Exp $
//
// Description:
//
//===========================================================================



#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "time.h"
#include <fstream>

using namespace Go;
using namespace std;

int main(int argc, char** argv)
{

  ASSERT(argc == 2 || argc == 3);

  ifstream file(argv[1]);
  double tol = 1e-6;
  if (argc == 3) tol = atof(argv[2]);

  ObjectHeader head;
  file >> head;
  if (head.classType() != SplineSurface::classType()) {
    THROW("Not a spline surface");
  }

  SplineSurface cv;
  file >> cv;

  double ti = -(double)clock();
  double area = cv.area(tol);
  ti += (double)clock();

  cout << "Time used is " << double(ti)/1000.0 << " ms" << endl;
  cout << setprecision(20);
  cout << "Estimated area is " << area << endl;

}
