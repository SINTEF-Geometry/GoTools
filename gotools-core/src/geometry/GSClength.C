//===========================================================================
//
// File : GSClength.C
//
// Created: Tue Jul  1 12:09:14 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: GSClength.C,v 1.1 2008-07-08 10:58:46 kfp Exp $
//
// Description:
//
//===========================================================================


#include "GoTools/geometry/SplineCurve.h"

using namespace std;
using namespace Go;



namespace Go
{


//===========================================================================
double SplineCurve::length(double tol)
//===========================================================================
{
  int deg = order() - 1;
  int num_spans = numCoefs() - deg;

  double result = 0.0;

  for (int i = 0; i < num_spans; ++i)
    {
      double start = knotsBegin()[deg + i];
      double end = knotsBegin()[deg + i + 1];
      result += ParamCurve::length(tol, start, end);
    }

  return result;
}


} // namspace Go
