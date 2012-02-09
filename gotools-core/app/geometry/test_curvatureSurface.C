//===========================================================================
//
// File : test_curvatureSurface.C
//
// Created: Mon Aug 11 14:31:47 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_curvatureSurface.C,v 1.2 2009-01-13 12:44:51 kfp Exp $
//
// Description:
//
//===========================================================================



#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/CurvatureAnalysis.h"
#include <fstream>

using namespace Go;
using namespace std;

int main(int argc, char** argv)
{

  ASSERT(argc == 2 || argc == 3 || argc == 4);

  ifstream file(argv[1]);
  double tol = 0.1;
  double gap = 0.001;
  if (argc >= 3) tol = atof(argv[2]);
  if (argc == 4) gap = atof(argv[3]);

  ObjectHeader head;
  file >> head;
  if (head.classType() != SplineSurface::classType()) {
    THROW("Not a spline surface");
  }

  SplineSurface sf;
  file >> sf;

  double mincurv, pos_u, pos_v;
  CurvatureAnalysis::minimalCurvatureRadius(sf, tol, mincurv, pos_u, pos_v, gap);

  Point pt;
  sf.point(pt, pos_u, pos_v);

  cout << "Minimal curvature radius: " << mincurv << endl;

  ofstream fout("point.g2");
  fout << "400 1 0 4 255 0 0 255" << endl;
  fout << "1" << endl;
  fout << pt << endl;

}
