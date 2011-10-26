//===========================================================================
//
// File : test_curvatureSurface.C
//
// Created: Mon Aug 11 14:31:47 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_curvatureCurve.C,v 1.1 2008-08-27 10:53:57 kfp Exp $
//
// Description:
//
//===========================================================================



#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Curvature.h"
#include <fstream>

using namespace Go;
using namespace std;

int main(int argc, char** argv)
{

  ASSERT(argc == 2);

  ifstream file(argv[1]);

  ObjectHeader head;
  file >> head;
  if (head.classType() != SplineCurve::classType()) {
    THROW("Not a spline curve");
  }

  SplineCurve sc;
  file >> sc;

  double mincurv, pos;
  minimalCurvatureRadius(sc, mincurv, pos);

  Point pt;
  sc.point(pt, pos);

  cout << "Minimal curvature radius: " << mincurv << endl;

  ofstream fout("point.g2");
  fout << "400 1 0 4 255 0 0 255" << endl;
  fout << "1" << endl;
  fout << pt << endl;

}
