//===========================================================================
//
// File : test_uniformEval.C
//
// Created: Fri Jul 10 12:24:50 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_uniformEval.C,v 1.1 2009/07/10 10:38:15 kfp Exp $
//
// Description:
//
//===========================================================================



#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/Point.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char** argv)
{

  if (argc != 4)
    {
      cout << "Usage: " << argv[0] << " curveinfile pointsoutfile nmb_pts" << endl;
      exit(-1);
    }

  ifstream filein(argv[1]);
  ALWAYS_ERROR_IF(filein.bad(), "Bad or no curvee input filename");
  ObjectHeader head;
  filein >> head;
  if (head.classType() != SplineCurve::classType()) {
    THROW("Not a spline curve");
  }
  SplineCurve sc;
  filein >> sc;

  ofstream fileout(argv[2]);
  ALWAYS_ERROR_IF(fileout.bad(), "Bad points output filename");

  vector<Point> pts;
  vector<double> pars;
  int nmbpts = atoi(argv[3]);

  sc.uniformEvaluator(nmbpts, pts, pars);

  fileout << "400 1 0 4 255 255 0 255" << endl;
  fileout << nmbpts << endl;
  for (int i = 0; i < nmbpts; ++i)
    fileout << pts[i] << endl;
}

