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

#include <array>
#if (defined(__GNUC__) || _MSC_VER > 1600) // No <thread> in VS2010
#include <thread>
#endif
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/utils/checks.h"
#include "GoTools/utils/StreamUtils.h"

// The following is a workaround since 'thread_local' is not well supported by compilers yet
#if defined(__GNUC__)
#define thread_local __thread
#elif _MSC_VER > 1600  //defined(_WIN32)
#define thread_local __declspec( thread )
#else
#define thread_local // _MSC_VER == 1600, i.e. VS2010
#endif

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
  if ((t < kvals[knot_ix[0]]) || (t > kvals[knot_ix[deg+1]])) return 0;

  assert(deg <= MAX_DEGREE);
  fill (tmp, tmp+deg+1, 0);

  // computing lowest-degree B-spline components (all zero except one)
  int nonzero_ix = 0;
  if (at_end)  while (kvals[knot_ix[nonzero_ix+1]] <  t) ++nonzero_ix;
  else         while (kvals[knot_ix[nonzero_ix+1]] <= t) ++nonzero_ix;

  if (nonzero_ix > deg)
    return 0.0; // Basis function defined to be 0.0 for value outside the support.
//  assert(nonzero_ix <= deg);

  tmp[nonzero_ix] = 1;

  // accumulating to attain correct degree

  for (int d = 1; d != deg+1; ++d) {
    const int lbound = max (0, nonzero_ix - d);
    const int ubound = min (nonzero_ix, deg - d);
    for (int i = lbound; i <= ubound; ++i) {
      const double k_i     = kvals[knot_ix[i]];
      const double k_ip1   = kvals[knot_ix[i+1]];
      const double k_ipd   = kvals[knot_ix[i+d]]; 
      const double k_ipdp1 = kvals[knot_ix[i+d+1]];
      const double alpha =  (k_ipd == k_i) ? 0 : (t - k_i) / (k_ipd - k_i);
      const double beta  =  (k_ipdp1 == k_ip1) ? 0 : (k_ipdp1 - t) / (k_ipdp1 - k_ip1);

      tmp[i] = alpha * tmp[i] + beta * tmp[i+1];
    }
  }

  return tmp[0];
}


// //------------------------------------------------------------------------------
// // B-spline evaluation function
// double B_recursive(int deg, double t, const int* knot_ix, const double* kvals, bool at_end)
// //------------------------------------------------------------------------------
// {
//   const double k0   = kvals[knot_ix[0]];
//   const double k1   = kvals[knot_ix[1]];
//   const double kd   = kvals[knot_ix[deg]];
//   const double kdp1 = kvals[knot_ix[deg+1]];

//   assert(deg >= 0);
//   if (deg == 0) 
//     return // at_end: half-open interval at end of domain should be considered differently
//       (! at_end) ? ((k0 <= t && t <  k1) ? 1 : 0) : 
//                    ((k0 <  t && t <= k1) ? 1 : 0);

//   const double fac1 = (kd   > k0) ? (t -   k0) / (kd - k0) : 0;
//   const double fac2 = (kdp1 > k1) ? (kdp1 - t) / (kdp1 - k1) : 0;
  
//   return 
//     ( (fac1 > 0) ? fac1 * B(deg-1, t, knot_ix, kvals, at_end) : 0) +
//     ( (fac2 > 0) ? fac2 * B(deg-1, t, knot_ix+1, kvals, at_end) : 0);
// }

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
    // ( (fac1 != 0) ? 
    //   ( fac1 * ( (der>1) ? 
    //     dB(deg-1, t, knot_ix, kvals, at_end, der-1)
    // 		 : B(deg-1, t, knot_ix, kvals, at_end) ) )
    //   : 0 ) + 
    // ( (fac2 != 0) ? 
    //   ( fac2 * ( (der>1) ? 
    // 	       dB(deg-1, t, knot_ix+1, kvals, at_end, der-1) 
    // 		 : B(deg-1, t, knot_ix+1, kvals, at_end) ) )
    //   : 0 ) ;
}

//------------------------------------------------------------------------------
double compute_univariate_spline(int deg, 
				 double u, 
				 const vector<int>& k_ixes, 
				 const double* kvals, 
				 int deriv,
				 bool on_end)
//------------------------------------------------------------------------------
{
  return (deriv>0) ? 
    dB(deg, u, &k_ixes[0], kvals, on_end, deriv) : 
    B( deg, u, &k_ixes[0], kvals, on_end);
}

}; // anonymous namespace

//==============================================================================
bool LRBSpline2D::operator<(const LRBSpline2D& rhs) const
//==============================================================================
{
  const int tmp1 = compare_seq(kvec_u_.begin(), kvec_u_.end(), rhs.kvec_u_.begin(), rhs.kvec_u_.end());
  if (tmp1 != 0) return (tmp1 < 0);

  const int tmp2 = compare_seq(kvec_v_.begin(), kvec_v_.end(), rhs.kvec_v_.begin(), rhs.kvec_v_.end());
  if (tmp2 != 0) return (tmp2 < 0);

  const int tmp3 = compare_seq(coef_times_gamma_.begin(), coef_times_gamma_.end(), 
			       rhs.coef_times_gamma_.begin(), rhs.coef_times_gamma_.end());
  if (tmp3 != 0) return (tmp3 < 0);

  return gamma_ < rhs.gamma_;
}

//==============================================================================
bool LRBSpline2D::operator==(const LRBSpline2D& rhs) const
//==============================================================================
{
#if 0//ndef NDEBUG
  // double umi = umin();
  // double uma = umax();
  // double vmi = vmin();
  // double vma = vmax();
  // std::cout << "umin: " << umi << ", umax: " << uma << ", vmin: " << vmi << ", vmax: " << vma << std::endl;
  auto iter = kvec_u_.begin();
  std::cout << "DEBUG: kvec_u_: ";
  while (iter != kvec_u_.end())
    {
      std::cout << *iter << " ";
      ++iter;
    }
  std::cout << "kvec_v_: ";
  iter = kvec_v_.begin();  
  while (iter != kvec_v_.end())
    {
      std::cout << *iter << " ";
      ++iter;
    }
  std::cout << std::endl;

  iter = rhs.kvec_u_.begin();  
  std::cout << "DEBUG: rhs.kvec_u_: ";
  while (iter != rhs.kvec_u_.end())
    {
      std::cout << *iter << " ";
      ++iter;
    }

  iter = rhs.kvec_v_.begin();  
  std::cout << "rhs.kvec_v_: ";
  while (iter != rhs.kvec_v_.end())
    {
      std::cout << *iter << " ";
      ++iter;
    }
  std::cout << std::endl;
  
#endif
#ifndef NDEBUG
  int kvec_u_size = kvec_u_.size();
  int kvec_v_size = kvec_v_.size();
  int kvec_u_size2 = rhs.kvec_u_.size();
  int kvec_v_size2 = rhs.kvec_v_.size();
  if ((kvec_u_size != kvec_u_size2) || (kvec_v_size != kvec_v_size2))
      MESSAGE("DEBUG: Pairwise vectors are of different size!");
//    ;
#endif

  const int tmp1 = compare_seq(kvec_u_.begin(), kvec_u_.end(), 
			       rhs.kvec_u_.begin(), rhs.kvec_u_.end());
  if (tmp1 != 0)
    return false;

  const int tmp2 = compare_seq(kvec_v_.begin(), kvec_v_.end(), 
			       rhs.kvec_v_.begin(), rhs.kvec_v_.end());
  if (tmp2 != 0)
    return false;

  return true;
}

 //==============================================================================
void LRBSpline2D::write(ostream& os) const
//==============================================================================
{
  // @@sbr201301 We must decide on a file format for the LRBSpline2D.
  // For rational case the dimension is currently written as dim + 1.
  // It makes more sense to keep geometric dimension and set rational boolean.
  int dim = coef_times_gamma_.dimension();
  object_to_stream(os, dim);
  int rat = (rational_) ? 1 : 0;
  object_to_stream(os, rat);
  object_to_stream(os, '\n');
  object_to_stream(os, coef_times_gamma_);
  object_to_stream(os, gamma_);
  object_to_stream(os, weight_);
  object_to_stream(os, '\n');
  object_to_stream(os, kvec_u_);
  object_to_stream(os, kvec_v_);
}

//==============================================================================
void LRBSpline2D::read(istream& is)
//==============================================================================
{
  // @@sbr201301 Currently we are expecting the rational weight to be
  // included in file format, even for non-rational cases.
  int dim = -1;
  object_from_stream(is, dim);
  coef_times_gamma_.resize(dim);
  int rat = -1;
  object_from_stream(is, rat);
  rational_ = (rat == 1);
  object_from_stream(is, coef_times_gamma_);
  object_from_stream(is, gamma_);
  // if (gamma_ < 1.0)
  // {
  //     MESSAGE("DEBUGGING: Changing gamma from " << gamma_ << " to 1.0!");
  //     coef_times_gamma_ /= gamma_;
  //     gamma_ = 1.0;
  // }
  object_from_stream(is, weight_);
  object_from_stream(is, kvec_u_);
  object_from_stream(is, kvec_v_);
  coef_fixed_ = 0;
}

//==============================================================================
double LRBSpline2D::evalBasisFunction(double u, 
					  double v, 
					  int u_deriv, 
					  int v_deriv,
					  bool u_at_end,
					  bool v_at_end) const
//==============================================================================
{
  return 
    compute_univariate_spline(degree(XFIXED), u, kvec(XFIXED), mesh_->knotsBegin(XFIXED), 
			      u_deriv, u_at_end) *
    compute_univariate_spline(degree(YFIXED), v, kvec(YFIXED), mesh_->knotsBegin(YFIXED),  
			      v_deriv, v_at_end);
}


//==============================================================================
Point LRBSpline2D::getGrevilleParameter() const
{
  MESSAGE("getGrevilleParameter(): Not implemented.");
  return Point();
}

//==============================================================================
bool LRBSpline2D::overlaps(Element2D *el) const
//==============================================================================
{
  // Does it make sense to include equality?
  if (el->umin() >= umax())
    return false;
  if (el->umax() <= umin())
    return false;
  if (el->vmin() >= vmax())
    return false;
  if (el->vmax() <= vmin())
    return false;
  
  return true;
}


//==============================================================================
bool LRBSpline2D::overlaps(double domain[]) const
//==============================================================================
{
  // Does it make sense to include equality?
  if (domain[0] >= umax())
    return false;
  if (domain[1] <= umin())
    return false;
  if (domain[2] >= vmax())
    return false;
  if (domain[3] <= vmin())
    return false;
  
  return true;
}

//==============================================================================
bool LRBSpline2D::addSupport(Element2D *el)
//==============================================================================
{
  for (size_t i=0; i<support_.size(); i++) {
    if(el == support_[i]) {
      return false; 
    }
  }
  support_.push_back(el);
  return true;
}

//==============================================================================
void LRBSpline2D::removeSupport(Element2D *el)
//==============================================================================
{
  for (size_t i=0; i<support_.size(); i++) {
    if(el == support_[i]) {
      if (i < support_.size() - 1)
	{
	  support_[i] = support_.back();
	  support_[support_.size()-1] = NULL;
	}
      support_.pop_back();
      return;
    }
  }
}

//==============================================================================
std::vector<Element2D*>::iterator LRBSpline2D::supportedElementBegin()
//==============================================================================
{
  return support_.begin();
}

//==============================================================================
std::vector<Element2D*>::iterator LRBSpline2D::supportedElementEnd()
//==============================================================================
{
  return support_.end();
}

//==============================================================================
void LRBSpline2D::subtractKnotIdx(int u_del, int v_del)
//==============================================================================
{
  for (size_t kj=0; kj<kvec_u_.size(); ++kj)
    kvec_u_[kj] -= u_del;

  for (size_t kj=0; kj<kvec_v_.size(); ++kj)
    kvec_v_[kj] -= v_del;
}

//==============================================================================
void LRBSpline2D::reverseParameterDirection(bool dir_is_u)
//==============================================================================
{
  int num_unique_knots = (dir_is_u) ? mesh_->numDistinctKnots(XFIXED) : mesh_->numDistinctKnots(YFIXED);
  auto iter_beg = (dir_is_u) ? kvec_u_.begin() : kvec_v_.begin();
  auto iter_end = (dir_is_u) ? kvec_u_.end() : kvec_v_.end();
  
  auto iter = iter_beg;
  while (iter != iter_end)
    {
      *iter = num_unique_knots - 1 - *iter;
      ++iter;
    }

  std::reverse(iter_beg, iter_end);
}


//==============================================================================
void LRBSpline2D::swapParameterDirection()
//==============================================================================
{
  std::swap(kvec_u_, kvec_v_);
}


//==============================================================================
vector<double> LRBSpline2D::unitSquareBernsteinBasis(double start_u, double stop_u, double start_v, double stop_v) const
//==============================================================================
{
  vector<double> result;

  vector<double> coefs_u = unitIntervalBernsteinBasis(start_u, stop_u, XFIXED);
  vector<double> coefs_v = unitIntervalBernsteinBasis(start_v, stop_v, YFIXED);

  for (vector<double>::const_iterator it_v = coefs_v.begin(); it_v != coefs_v.end(); ++it_v)
    {
      Point coefs_point = coef_times_gamma_ * (*it_v);
      for (vector<double>::const_iterator it_u = coefs_u.begin(); it_u != coefs_u.end(); ++it_u)
	for (int i = 0; i < coefs_point.size(); ++i)
	  result.push_back(coefs_point[i] * (*it_u));
    }

  return result;
}


//==============================================================================
vector<double> LRBSpline2D::unitIntervalBernsteinBasis(double start, double stop, Direction2D d) const
//==============================================================================
{
  // Get knot vector, where the knots are translated by start -> 0.0 and stop -> 1.0
  vector<double> knots;
  vector<int> knots_int = kvec(d);

  double slope = 1.0/(stop - start);
  int deg = degree(d);

  for (int i = 0; i < deg + 2; ++i)
    knots.push_back(slope * (mesh_->kval(d,knots_int[i]) - start));

  // Get the position of the interval containing [0,1]. We assume that for
  // some k, knots[k] <= 0.0 and knots[k+1] >= 1.0, and let interval_pos be this k.
  // We use 0.5 instead of 1.0 to break the loop, in order to avoid using tolerances.
  // Any number in the open interval (0,1) would work.
  int interval_pos;
  for (interval_pos = 0; interval_pos <= deg; ++interval_pos)
    if (knots[interval_pos + 1] >= 0.5)
      break;

  // Prepare array holding the Bernstein basis coefficients.
  // After each step for each polynomial degree (value of k in outermost loop below),
  // the polynomial part on the interval [ knots[interval_pos], knots[interval_pos+1] ])
  // of the k-degree B-spline defined by knot vector knot[i],...,knot[i+k+1] is given
  // by coefficients coefs[i][0],...,coefs[i][k]. At the end, the coefficients to be
  // returned are in coefs[0]
  vector<vector<double> > coefs(deg+1);
  for (int i = 0; i <= deg; ++i)
    coefs[i].resize(deg + 1 - i);
  coefs[interval_pos][0] = 1.0;

  for (int k = 1; k <=deg; ++k)
    for (int i = 0; i <= deg - k; ++i)
      if (i >= interval_pos - k && i <= interval_pos)   // Only look at B-splines with support in interval
      {
	double coefs_i_jmin1 = 0.0;  // For caching coefs[i][j-1] in inner loop

	// Store 1/(k*(knots[i+k]-knots[i])) and same for next interval. The denominator should not be zero
	// (because knots[interval_pos] < knots[interval_pos +1]) but just in case we use the standard
	// assumption 1/0 = 0 from spline arithmetics
	double denom_0 = (double)k*(knots[i + k] - knots[i]);
	if (denom_0 != 0.0)
	  denom_0 = 1.0/denom_0;
	double denom_1 = (double)k*(knots[i + k + 1] - knots[i + 1]);
	if (denom_1 != 0.0)
	  denom_1 = 1.0/denom_1;

	// Some factors used several times
	double f0 = (1.0 - knots[i]) * denom_0;
	double f1 = (knots[i + k + 1] - 1.0) * denom_1;
	double f2 = f0 - denom_0;
	double f3 = f1 + denom_1;

	// Calculate the new coefficients
	for (int j = 0; j <= k; ++j)
	  {
	    double res = 0.0;
	    if (j > 0)
	      res += (f0 * coefs_i_jmin1 + f1 * coefs[i + 1][j - 1]) * (double)j;
	    if (j < k)
	      res += (f2 * coefs[i][j] + f3 * coefs[i + 1][j]) * (double)(k - j);
	    coefs_i_jmin1 = coefs[i][j];
	    coefs[i][j] = res;
	  }
      }

  return coefs[0];
}


}; // end namespace Go
