//===========================================================================
//
// File : test_curvatureCurvePts.C
//
// Created: Tue Aug 26 11:36:13 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_curvatureCurvePts.C,v 1.1 2008-08-27 10:53:57 kfp Exp $
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

  ASSERT(argc == 3);

  ifstream file(argv[1]);
  double curveRad = atof(argv[2]);

  ObjectHeader head;
  file >> head;
  if (head.classType() != SplineCurve::classType()) {
    THROW("Not a spline curve");
  }

  SplineCurve sc;
  file >> sc;

  vector<double> positions;

  curvatureRadiusPoints(sc, curveRad, positions);

  if (positions.size() == 0)
    {
      cout << "No points found" << endl;
      return 0;
    }

  ofstream fout("point.g2");
  for (size_t i = 0; i < positions.size(); ++i)
    {
      Point pt;
      sc.point(pt, positions[i]);
      fout << "400 1 0 4 255 0 0 255" << endl;
      fout << "1" << endl;
      fout << pt << endl;
      cout << positions[i] << endl;
    }

}
