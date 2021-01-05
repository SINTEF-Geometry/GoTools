/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/lrsplines2D/BSplineUniLR.h"
// #include "GoTools/utils/checks.h"
// #include "GoTools/utils/StreamUtils.h"
#include "GoTools/geometry/BsplineBasis.h"
#include <math.h>

//#define DEBUG

using namespace std;

namespace Go
{


//------------------------------------------------------------------------------
namespace
//------------------------------------------------------------------------------
{
// Since some static buffers (provided for efficiency reasons) need to know the maximum degree
// used at compile time, the following constant, MAX_DEGREE, is here defined.
const int MAX_DEGREE = 20;
  const int MAX_DER = 3;
  const int MAX_DIM = 3;

//------------------------------------------------------------------------------
double B(int deg, double t, const int* knot_ix, const double* kvals, bool at_end)
//------------------------------------------------------------------------------
{
  // a POD rather than a stl vector used below due to the limitations of thread_local as currently
  // defined (see #defines at the top of this file).  A practical consequence is that 
  // MAX_DEGREE must be known at compile time.
  //static double thread_local tmp[MAX_DEGREE+2];
  double tmp[MAX_DEGREE+2];

  // only evaluate if within support
  if ((t < kvals[knot_ix[0]]) || (t > kvals[knot_ix[deg+1]])) 
    return 0;

  assert(deg <= MAX_DEGREE);
  fill (tmp, tmp+deg+1, 0);

  // computing lowest-degree B-spline components (all zero except one)
  int nonzero_ix = 0;
  if (at_end)  
    while (kvals[knot_ix[nonzero_ix+1]] <  t) 
      ++nonzero_ix;
  else         
    while (nonzero_ix <= deg && kvals[knot_ix[nonzero_ix+1]] <= t) 
      ++nonzero_ix;

  if (nonzero_ix > deg+1)
    return 0.0; // Basis function defined to be 0.0 for value outside the support.
//  assert(nonzero_ix <= deg);

  tmp[nonzero_ix] = 1;

  // accumulating to attain correct degree

  double alpha, beta;
  double tt1, tt2, tt3, tt4, td1, td2;
  for (int d = 1; d != deg+1; ++d) {
    const int lbound = max (0, nonzero_ix - d);
    const int ubound = min (nonzero_ix, deg - d);
    tt1 = kvals[knot_ix[lbound]];
    tt3 = kvals[knot_ix[lbound+d]];
    td1 = tt3 - tt1;
    if (d <= nonzero_ix && lbound <= ubound)
      {
	tt2 = kvals[knot_ix[lbound+1]];
	tt4 = kvals[knot_ix[lbound+d+1]];
	td2 = tt4 - tt2;
	beta = (tt4 == tt2) ? 0.0  : (tt4 - t)/td2;
	tmp[lbound] = beta*tmp[lbound+1];
	tt1 = tt2;
	tt3 = tt4;
	td1 = td2; 
      }
    for (int i = lbound+(d <=nonzero_ix); i <= ubound; ++i) 
      {
	tt2 = kvals[knot_ix[i+1]];
	tt4 = kvals[knot_ix[i+d+1]];
	td2 = tt4 - tt2;
	alpha = (tt3 == tt1) ? 0.0 : (t - tt1)/td1;
	beta = (tt2 == tt4) ? 0.0 : (tt4 - t)/td2;
	tmp[i] = alpha * tmp[i] + beta * tmp[i+1];
	tt1 = tt2;
	tt3 = tt4;
	td1 = td2; 
    }
  }

  return tmp[0];
}

//------------------------------------------------------------------------------
// B-spline derivative evaluation
double dB(int deg, double t, const int* knot_ix, const double* kvals, bool at_end, int der=1)
//------------------------------------------------------------------------------
{
  const double k0   = kvals[knot_ix[0]];
  const double k1   = kvals[knot_ix[1]];
  const double kdeg = kvals[knot_ix[deg]];
  const double kdp1 = kvals[knot_ix[deg+1]];

  assert(der >  0); //@@ we should perhaps check that derivative <=
		    //   degree - multiplicity
  if (deg == 0) 
    return 0;
  double fac1 = (kdeg > k0) ? ( deg) / (kdeg - k0) : 0;
  double fac2 = (kdp1 > k1) ? (-deg) / (kdp1 - k1) : 0;

  double part1 = (fac1 != 0) ? 
    ( fac1 * ( (der>1) ? 
	       dB(deg-1, t, knot_ix, kvals, at_end, der-1)
	       : B(deg-1, t, knot_ix, kvals, at_end) ) ) : 0.0;

  double part2 = (fac2 != 0) ? 
      ( fac2 * ( (der>1) ? 
	       dB(deg-1, t, knot_ix+1, kvals, at_end, der-1) 
		 : B(deg-1, t, knot_ix+1, kvals, at_end) ) ) : 0.0;

  return part1 + part2; // The product rule.
}

//------------------------------------------------------------------------------
  void Bder(const int& deg, const double& t, int& nder, const int* knot_ix, 
	    const double* kvals, double der[], const bool& at_end)
//------------------------------------------------------------------------------
{
  // a POD rather than a stl vector used below due to the limitations of thread_local as currently
  // defined (see #defines at the top of this file).  A practical consequence is that 
  // MAX_DEGREE must be known at compile time.
  //static double thread_local tmp[MAX_DEGREE+2];
  double tmp[MAX_DEGREE+8]; // Assumes maximum number of derivatives equal to 3

  // Adjust derivative if too large
  nder = min(nder, deg);

  // only evaluate if within support
  if ((t < kvals[knot_ix[0]]) || (t > kvals[knot_ix[deg+1]])) 
    return;

  assert(deg <= MAX_DEGREE);
  fill (tmp, tmp+deg+1+2*nder, 0);

  // computing lowest-degree B-spline components (all zero except one)
  int nonzero_ix = 0;
  if (at_end)  
    while (kvals[knot_ix[nonzero_ix+1]] <  t) 
      ++nonzero_ix;
  else         
    while (nonzero_ix <= deg && kvals[knot_ix[nonzero_ix+1]] <= t) 
      ++nonzero_ix;

  if (nonzero_ix > deg+1)
    return; // Basis function defined to be 0.0 for value outside the support.
//  assert(nonzero_ix <= deg);

  tmp[nonzero_ix] = 1;

  // accumulating to attain correct degree

  double alpha, beta;
  double tt1, tt2, tt3, tt4, td1, td2;
  int kcurr, kprev, kcurr2;
  int i, j, k, h;
  for (int d = 1; d != deg + 1; ++d) {
    const int lbound = max (0, nonzero_ix - d);
    const int ubound = min (nonzero_ix, deg - d);
#if 1
    for (j=std::min(nder,d); j>0; --j)
      {
	kcurr = deg - d + 1 + j*(deg-d+1);
	kcurr2 = (d == 1) ? 0 : deg - d + 2 + (j-1)*(deg - d + 2);
	kcurr = std::max(kcurr, kcurr2);
	int kstop = kcurr - deg + d;
	kprev = (j == 1) ? 0 : kstop - 1;
	h = deg-d;
	tt2 = kvals[knot_ix[h+1]];
	tt4 = kvals[knot_ix[h+d+1]];
	td2 = (tt2 != tt4) ? 1.0/(tt4 - tt2) : 0.0;
	for (i=kcurr, k=kprev+deg-d; i>=kstop; --i, --k, --h)
	  {
	    tt1 = kvals[knot_ix[h]];
	    tt3 = kvals[knot_ix[h+d]];
	    td1 = (tt1 != tt3) ? 1.0/(tt3 - tt1) : 0.0;

	    tmp[i] = d*(td1*tmp[k] - td2*tmp[k+1]);
	    tt2 = tt1;
	    tt3 = tt3;
	    td2 = td1; 
	  }
      }

    // tt1 = kvals[knot_ix[0]];
    // tt3 = kvals[knot_ix[d]];
    tt1 = kvals[knot_ix[lbound]];
    tt3 = kvals[knot_ix[lbound+d]];
     td1 = (tt1 != tt3) ? 1.0/(tt3 - tt1) : 0.0;
    //for (i=0; i<=deg-d; ++i)
    for (i=lbound; i<=ubound; ++i)
      {
	tt2 = kvals[knot_ix[i+1]];
	tt4 = kvals[knot_ix[i+d+1]];
	td2 = (tt2 != tt4) ? 1.0/(tt4 - tt2) : 0.0;
	alpha = td1*(t-tt1);
	beta = td2*(tt4-t);
	tmp[i] = alpha * tmp[i] + beta * tmp[i+1];
	tt1 = tt2;
	tt3 = tt4;
	td1 = td2; 
      }
    tmp[deg-d+1] = 0;
#endif

#if 0
    tt1 = kvals[knot_ix[lbound]];
    tt3 = kvals[knot_ix[lbound+d]];
    td1 = (tt1 != tt3) ? 1.0/(tt3 - tt1) : 0.0;

    if (d <= nonzero_ix && lbound <= ubound)
      {
   	tt2 = kvals[knot_ix[lbound+1]];
   	tt4 = kvals[knot_ix[lbound+d+1]];
   	td2 = (tt2 != tt4) ? 1.0/(tt4 - tt2) : 0.0;
   	for (int j=nder; j>deg-d; --j)
   	  {
   	    int k = j + d - deg - 1;
   	    tmp[2*(k+1)+lbound] = -d*td2*tmp[lbound+2*k+1];
   	  }
   	beta = td2*(tt4-t);
   	tmp[lbound] = beta*tmp[lbound+1];
   	tt1 = tt2;
   	tt3 = tt4;
   	td1 = td2; 
      }
    for (int i = lbound+(d <=nonzero_ix); i <= ubound; ++i) 
      {
    	tt2 = kvals[knot_ix[i+1]];
    	tt4 = kvals[knot_ix[i+d+1]];
    	td2 = (tt2 != tt4) ? 1.0/(tt4 - tt2) : 0.0;
    	for (int j=nder; j>deg-d; --j)
    	  {
    	    int k = j + d - deg - 1;
    	    tmp[2*(k+1)+i] = d*(td1*tmp[i+2*k] - td2*tmp[i+2*k+1]);
    	  }

    	alpha = td1*(t-tt1);
    	beta = td2*(tt4-t);
    	tmp[i] = alpha * tmp[i] + beta * tmp[i+1];
    	tt1 = tt2;
    	tt3 = tt4;
    	td1 = td2; 
    }
    tmp[ubound+1] = 0.0;
#endif
  }

  der[0] = tmp[0];
#ifdef DEBUG
  double der2[4];
  double val = B(deg, t, knot_ix, kvals, at_end);
  if (fabs(der[0]-val) > 1.0e-6)
    std::cout << "Bspline evaluation mismatch, position: " << der[0] << ", " << val << std::endl;
#endif
  for (int i=1; i<=nder; ++i)
    {
      der[i] = tmp[i*2];
#ifdef DEBUG
      der2[i-1] = dB(deg, t, knot_ix, kvals, at_end, i);
      if (fabs(der[i]-der2[i-1]) > 1.0e-6)
	std::cout << "Bspline evaluation mismatch, der="<< i << ": " << der[i] << ", " << der2[i-1] << std::endl;
#endif
    }
  int stop_break;
}


}; // anonymous namespace


//==============================================================================
BSplineUniLR::BSplineUniLR(const BSplineUniLR& rhs)
//==============================================================================
{
  kvec_.insert(kvec_.end(), rhs.kvec_.begin(), rhs.kvec_.end());
  pardir_ = rhs.pardir_;
  mesh_ = rhs.mesh_;
  count_ = 0;
}

// //==============================================================================
// int BSplineUniLR::operator<(const BSplineUniLR& rhs) const
// //==============================================================================
// {
//   return compare_seq(kvec_.begin(), kvec_.end(), rhs.kvec_.begin(), rhs.kvec_.end());
// }

//==============================================================================
bool BSplineUniLR::operator==(const BSplineUniLR& rhs) const
//==============================================================================
{
  const int tmp = compare_seq(kvec_.begin(), kvec_.end(), 
			      rhs.kvec_.begin(), rhs.kvec_.end());
  return (tmp == 0);
}

//==============================================================================
void BSplineUniLR::write(ostream& os) const
//==============================================================================
{
  object_to_stream(os, kvec_);
}

//==============================================================================
void BSplineUniLR::read(istream& is) 
//==============================================================================
{
  object_from_stream(is, kvec_);
}

//==============================================================================
double BSplineUniLR::evalBasisFunc(double par) const
//==============================================================================
{
  bool at_end = (par >= max());
  return B(degree(), par, &kvec_[0], mesh_->knotsBegin(pardir_), at_end);
}

//==============================================================================
  double BSplineUniLR::evalBasisFunction(double par, int deriv, 
					bool at_end) const
//==============================================================================
{
  return (deriv > 0) ?
    dB(degree(), par, &kvec_[0], mesh_->knotsBegin(pardir_), at_end, deriv) :
    B(degree(), par, &kvec_[0], mesh_->knotsBegin(pardir_), at_end);
}

//==============================================================================
  void BSplineUniLR::evalBasisFunctions(double par, int deriv, 
				       double der[], bool at_end) const
//==============================================================================
{
  Bder(degree(), par, deriv, &kvec_[0], mesh_->knotsBegin(pardir_),
       der, at_end);
}


//==============================================================================
  int BSplineUniLR::endmult(bool atstart) const
//==============================================================================
{
  int idx;
  size_t ki;
  if (atstart)
    {
      idx = 1;
      for (ki=1; ki<kvec_.size(); ++ki)
	{
	  if (kvec_[ki] != kvec_[ki-1])
	    break;
	  idx++;
	}
    }
  else
    {
      idx = 1;
      for (ki=kvec_.size()-2; ki>=0; --ki)
	{
	  if (kvec_[ki] != kvec_[ki+1])
	    break;
	  idx++;
	}
     }
  return idx;
}

//==============================================================================
  double BSplineUniLR::getGrevilleParameter() const
//==============================================================================
{
  int nmb = (int)kvec_.size()-1;
  int ki;
  double par = 0;
  for (ki=1; ki<nmb; ++ki)
    par += mesh_->kval(pardir_, kvec_[ki]);
  par /= (double)(nmb-1);

  return par;
}

//==============================================================================
  bool BSplineUniLR::overlaps(double pmin, double pmax) const
//==============================================================================
{
  // Does it make sense to include equality?
  if (pmin >= max())
    return false;
  if (pmax <= min())
    return false;

  return true;
}

//==============================================================================
  void BSplineUniLR::subtractKnotIdx(int del)
//==============================================================================
{
  for (size_t kj=0; kj<kvec_.size(); ++kj)
    kvec_[kj] -= del;
}

//==============================================================================
  void BSplineUniLR::reverseParameterDirection()
//==============================================================================
{
  int num_unique_knots = mesh_->numDistinctKnots(pardir_);
  auto iter_beg = kvec_.begin();
  auto iter_end = kvec_.end();
  
  auto iter = iter_beg;
  while (iter != iter_end)
    {
      *iter = num_unique_knots - 1 - *iter;
      ++iter;
    }

  std::reverse(iter_beg, iter_end);
}

}
