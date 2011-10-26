/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

/********************************************************************
 FILENAME    : PrThreshold.C
 AUTHOR      : Michael Floater, SINTEF
 DATE        : Sep. 98
 DESCRIPTION : Implementation of methods in the class PrThreshold.
 CHANGE LOG  :
*********************************************************************/

#include "GoTools/parametrization/PrParamUtil.h"
#include "GoTools/parametrization/PrThreshold.h"
#include "GoTools/parametrization/PrHeap.h"

#ifdef __BORLANDC__
using std::sqrt;
#endif

using std::cerr;
using std::cout;
using std::endl;

//-----------------------------------------------------------------------------
PrThreshold::PrThreshold()
//-----------------------------------------------------------------------------
{
  threshold_ = 0.0;
}

//-----------------------------------------------------------------------------
void PrThreshold::attach(std::shared_ptr<PrNestedTriangulation> t)
//-----------------------------------------------------------------------------
{
  t_ = t;
}

//-----------------------------------------------------------------------------
void
PrThreshold::threshold(double& comp_rate, double& wavelet_comp_rate,
                       int error_out)
//-----------------------------------------------------------------------------
// Threshold the wavelet coefficients by setting those
// whose absolute value is less than the tolerance to 0.
// Return the compression rate: ratio of no. of coarse hats plus
// non-zero wavelets over no. of coarse hats plus all wavelets.
// Return the wavelet compression rate: ratio of no. of non-zero
// over no. of wavelets.
//
// If error_out = 1, do the opposite of thresholding,
// i.e, retain only those wavelet coefficients which are smaller
// in absolute value than the tolerance.
// The resulting piecewise linear function represents the error after
// thresholding and can be used to find the l2 error, max error, etc.

{
  int t = t_->getFinestLevel();
  int no_left = 0;
  int num_coeff;

  if(error_out == 1)
  {
    for(int i = 0; i < t_->getNumNodes(0); i++)
    {
      t_->setX(i,0.0);
      t_->setY(i,0.0);
      t_->setZ(i,0.0);
    }
  }

  for(int j = 1; j <= t; j++)
  {
    num_coeff = 0;
    for(int i = t_->getNumNodes(j-1); i < t_->getNumNodes(j); i++)
    {
      Vector3D vec(t_->getX(i),t_->getY(i),t_->getZ(i));
      if(error_out == 0)
      {
        if(vec.length() < threshold_)
        {
          t_->setX(i,0.0);
          t_->setY(i,0.0);
          t_->setZ(i,0.0);
        }
        else num_coeff++;
      }
      else
      {
        if(vec.length() >= threshold_)
        {
          t_->setX(i,0.0);
          t_->setY(i,0.0);
          t_->setZ(i,0.0);
          num_coeff++;
        }
      }
    }

    no_left += num_coeff;
    
//      cout << "Level = " << j << endl;
//      cout << "number of wavelets = " <<
//                 t_->getNumNodes(j) - t_->getNumNodes(j-1) << endl;
//      cout << "number wavelets retained = " << num_coeff << endl;
      comp_rate = (double)num_coeff /
                     (double)(t_->getNumNodes(j) - t_->getNumNodes(j-1));
//      cout << "compression rate = " << comp_rate << endl;
  }

  cout << endl;
  comp_rate = (double)(t_->getNumNodes(0) + no_left)
              / (double)(t_->getNumNodes(t));
  wavelet_comp_rate = (double)no_left
              / (double)(t_->getNumNodes(t) - t_->getNumNodes(0));

//    cout << "total number of wavelets = " <<
//                    t_->getNumNodes(t) - t_->getNumNodes(0) << endl;
//    cout << "total number of retained wavelet coefficients = "
//                   << no_left << endl;
}

//-----------------------------------------------------------------------------
double PrThreshold::findThreshold(double& wavelet_comp_rate)
//-----------------------------------------------------------------------------
// Find a threshold which will give the given compressiosn rate.
// E.g. if wavelet_comp_rate = 0.05 and there are n wavelets then
// we want a threshold so that m = int(0.05 * n) are retained.
{
  int i,jlast;
  jlast = t_->getFinestLevel();
  int nmax = t_->getNumNodes(jlast);
  int nmin = t_->getNumNodes(0);
  int n = nmax - nmin; // num wavelets
  int m = (int)(wavelet_comp_rate * (double)n); // required num retained wavs.

  PrHeap heap(n);
  for(i = nmin; i < nmax; i++)
  {
    Vector3D vec(t_->getX(i),t_->getY(i),t_->getZ(i));
    heap.push(vec.length(),i-nmin+1);
  }

  double key;
  int index;
  for(i = 0; i < n - m; i++)
  {
    heap.pop(key,index);
        // numbers are popped in increasing order
    // s_o << key << endl;
  }
  return key;
}

//-----------------------------------------------------------------------------
void PrThreshold::thresholdByCompRate(double& wavelet_comp_rate,
                       int error_out)
//-----------------------------------------------------------------------------
// Threshold by given compression rate, e.g. 0.05.
{
  int i,jlast;
  jlast = t_->getFinestLevel();
  int nmax = t_->getNumNodes(jlast);
  int nmin = t_->getNumNodes(0);
  int n = nmax - nmin; // num wavelets
  int m = (int)(wavelet_comp_rate * (double)n); // required num retained wavs.

//    cout << "nmin = " << nmin << endl;
//    cout << "n = " << n << endl;
//    cout << "m = " << m << endl;

  if(error_out == 1)
  {
    for(int i = 0; i < nmin; i++)
    {
      t_->setX(i,0.0);
      t_->setY(i,0.0);
      t_->setZ(i,0.0);
    }
  }

  PrHeap heap(n);
  for(i = nmin; i < nmax; i++)
  {
    Vector3D vec(t_->getX(i),t_->getY(i),t_->getZ(i));
    heap.push(vec.length(),i-nmin+1);
  }

  double key;
  int index;
  if(error_out == 0)
  {
    for(i = 0; i < n - m; i++)
    {
      heap.pop(key,index);
          // numbers are popped in increasing order
      // s_o << key << endl;
      t_->setX(index+nmin-1,0.0);
      t_->setY(index+nmin-1,0.0);
      t_->setZ(index+nmin-1,0.0);
    }
  }
  else
  {
    for(i = 0; i < n - m; i++)
    {
      heap.pop(key,index);
          // numbers are popped in increasing order
      // s_o << key << endl;
    }
    for(i = 0; i < m; i++)
    {
      heap.pop(key,index);
          // numbers are popped in increasing order
      // s_o << key << endl;
      t_->setX(index+nmin-1,0.0);
      t_->setY(index+nmin-1,0.0);
      t_->setZ(index+nmin-1,0.0);
    }
  }

//    cout << "total number of wavelets = " <<
//                    t_->getNumNodes(jlast) - t_->getNumNodes(0) << endl;
//    cout << "total number of retained wavelet coefficients = " << m << endl;
}

//-----------------------------------------------------------------------------
void PrThreshold::truncateLevel(int jlev)
//-----------------------------------------------------------------------------
// Set all coeffs at level jlev to zero (0 <= jlev <= getFinestLevel).
{
  int istart = (jlev == 0 ? 0 : t_->getNumNodes(jlev-1));
  for(int i = istart; i < t_->getNumNodes(jlev); i++)
  {
    t_->setX(i,0.0);
    t_->setY(i,0.0);
    t_->setZ(i,0.0);
  }
}

//-----------------------------------------------------------------------------
void PrThreshold::leaveLevel(int jlev)
//-----------------------------------------------------------------------------
// Set all level zero coeffs to zero and all
// wavelet coeffs except those of level j.
{
  int jlast = t_->getFinestLevel();
  int j;
  for(j = 0; j < jlev; j++) truncateLevel(j);
  for(j = jlev+1; j <= jlast; j++) truncateLevel(j);
}

//-----------------------------------------------------------------------------
double PrThreshold::maxNorm()
//-----------------------------------------------------------------------------
// Calculate the max norm.
{
  double maxerror = 0.0;
  double error;
  int i;
  for(i = 0; i < t_->getNumNodes(t_->getFinestLevel()); i++)
  {
    Vector3D vec(t_->getX(i),t_->getY(i),t_->getZ(i));
    error = vec.length();
    if(error > maxerror) maxerror = error;
  }
  return maxerror;
}

//-----------------------------------------------------------------------------
double PrThreshold::averageAbsValue()
//-----------------------------------------------------------------------------
// Calculate the average absolute value.
{
  double average = 0.0;
  double absval;
  int i;
  for(i = 0; i < t_->getNumNodes(t_->getFinestLevel()); i++)
  {
    Vector3D vec(t_->getX(i),t_->getY(i),t_->getZ(i));
    absval = vec.length();
    average += absval;
  }
  average /= t_->getNumNodes(t_->getFinestLevel());
  return average;
}

//----------------------------------------------------------------------------
double PrThreshold::weightedL2Norm()
//-----------------------------------------------------------------------------
{
  int jlast = t_->getFinestLevel();
  int powerOfTwo = 1 << jlast; // 2 to the power jlast.

  int npts = t_->getNumNodes(jlast);
  vector<int> nghbrs;

  int i,j,deg;
  double sum = 0.0;

  for(i=0; i< npts; i++)
  {
    t_->getNeighbours(i,jlast,nghbrs); 
    deg = (int)nghbrs.size();
    for(j=0; j<deg-1; j++)
    {
      if(i < nghbrs[j] && i < nghbrs[j+1])
      {
        sum += triangleNorm(i,nghbrs[j],nghbrs[j+1]);
      }
    }
    if(!t_->isBoundary(i))
    {
      if(i < nghbrs[deg-1] && i < nghbrs[0])
      {
        sum += triangleNorm(i,nghbrs[deg-1],nghbrs[0]);
      }
    }
  }

  return sqrt(sum / 6.0) / powerOfTwo;
}

//----------------------------------------------------------------------------
double PrThreshold::triangleNorm(int i, int j, int k)
//-----------------------------------------------------------------------------
{
  return (  t_->getX(i) * t_->getX(i)
          + t_->getX(j) * t_->getX(j)
          + t_->getX(k) * t_->getX(k)
          + t_->getX(i) * t_->getX(j)
          + t_->getX(j) * t_->getX(k)
          + t_->getX(k) * t_->getX(i) )
       + (  t_->getY(i) * t_->getY(i)
          + t_->getY(j) * t_->getY(j)
          + t_->getY(k) * t_->getY(k)
          + t_->getY(i) * t_->getY(j)
          + t_->getY(j) * t_->getY(k)
          + t_->getY(k) * t_->getY(i) )
       + (  t_->getZ(i) * t_->getZ(i)
          + t_->getZ(j) * t_->getZ(j)
          + t_->getZ(k) * t_->getZ(k)
          + t_->getZ(i) * t_->getZ(j)
          + t_->getZ(j) * t_->getZ(k)
          + t_->getZ(k) * t_->getZ(i) );
}

/*
//----------------------------------------------------------------------------
double PrThreshold::L2Norm()
//-----------------------------------------------------------------------------
{
  int jlast = t_->getFinestLevel();
  int npts = t_->getNumNodes(jlast);
  vector<int> nghbrs;

  int i,j,deg,it = 0;
  double sum = 0.0;

  for(i=0; i< npts; i++)
  {
    t_->getNeighbours(i,jlast,nghbrs); 
    deg = nghbrs.size();
    for(j=0; j<-deg; j++)
    {
      if(i < nghbrs[j] && i < nghbrs[j+1])
      {
        sum += triangleAreaNorm(i,nghbrs[j],nghbrs[j+1]);
      }
    }
    if(!t_->isBoundary(i))
    {
      if(i < nghbrs[deg-1] && i < nghbrs[0])
      {
        sum += triangleAreaNorm(i,nghbrs[deg-1],nghbrs[0]);
      }
    }
  }

  return sqrt(sum / 6.0);
}

//----------------------------------------------------------------------------
double PrThreshold::triangleAreaNorm(int i, int j, int k)
//-----------------------------------------------------------------------------
{
  Vector3D v1(t_->getU(i),t_->getV(i),t_->getW(i));
  Vector3D v2(t_->getU(j),t_->getV(j),t_->getW(j));
  Vector3D v3(t_->getU(k),t_->getV(k),t_->getW(k));

  Vector3D d1 = v2 - v1;
  Vector3D d2 = v3 - v1;

  double triangle_area = 0.5 * (d1.crossProd(d2)).length();
  return triangle_area * triangleNorm(i,j,k);
}
*/

