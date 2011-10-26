//===========================================================================
//
// File : get_isocurve.C
//
// Created: Tue Dec 16 13:42:48 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace Go;
using namespace std;

int main(int argc, char** argv)
{

  if (argc != 5)
    {
      cout << "Usage: " << argv[0] << " surfaceinfile isocurveoutfile direction(0,1) parametervalue" << endl;
      exit(-1);
    }

  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");


  // Read volume from file
  ObjectHeader head;
  SplineSurface surf;
  is >> head >> surf;

  // Open outfile
  ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  SplineCurve* curve = surf.constParamCurve(atof(argv[4]), atoi(argv[3])==0);

  curve->writeStandardHeader(os);
  curve->write(os);
}

