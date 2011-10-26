#include <math.h>
#include <stdlib.h>
#include "GoTools/parametrization/PrWaveletUtil.h"

//-----------------------------------------------------------------------------
double theta(int j, int k, int i, int deg, bool isBoundary)
//-----------------------------------------------------------------------------
// Given a coarse vertex v of degree deg whose index is i
// and two neighbouring fine vertices u1 and u2 which are the (j+1)-th and
// (k+1)-th neighbours resp., in an anticlockwise direction,
// (i.e. j and k begin theie numbering at 0)
// (which may possibly be equal), return the function
// theta(u1,u2,v), described in a paper by Floater and Quak.
{
  double lambda = 0.5 * (sqrt(21.0) - 5.0);
  if(isBoundary)
  {
    int a1 = (j <= k ? j : k);             // a1 = min(j,k)
    int a2 = deg - (j >= k ? j : k) - 1;   // a2 = deg - max(j,k) - 1
    return   ( pow(lambda,(double)a1) + pow(lambda,(double)(-a1)) )
           * ( pow(lambda,(double)a2) + pow(lambda,(double)(-a2)) )
           / ( pow(lambda,(double)(-deg+1)) - pow(lambda,(double)(deg-1)) )
           / sqrt(21.0);
  }
  else
  {
    int a = abs(k-j);
    return ( pow(lambda,(double)a) + pow(lambda,(double)(deg - a)) )
           / ( 1.0 - pow(lambda,(double)(deg)) )
           / sqrt(21.0);
  }
}

