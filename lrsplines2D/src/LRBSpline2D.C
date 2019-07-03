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
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
#include "GoTools/utils/checks.h"
#include "GoTools/utils/StreamUtils.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/SplineUtils.h"

//#define DEBUG

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
  const int MAX_DER = 3;
  const int MAX_DIM = 3;


}; // anonymous namespace


//==============================================================================
LRBSpline2D::LRBSpline2D(const LRBSpline2D& rhs)
//==============================================================================
{
  coef_fixed_ = rhs.coef_fixed_;
  coef_times_gamma_ = rhs.coef_times_gamma_;
  gamma_ = rhs.gamma_;
  bspline_u_ = rhs.bspline_u_;
  bspline_u_->incrCount();  // Initial count is zero
  bspline_v_ = rhs.bspline_v_;
  bspline_v_->incrCount();
  rational_ = rhs.rational_;
  // don't copy the support
  weight_ = rhs.weight_;

}

//==============================================================================
bool LRBSpline2D::operator<(const LRBSpline2D& rhs) const
//==============================================================================
{
  const int tmp1 = ((*bspline_u_) < (*rhs.bspline_u_));
  if (tmp1 != 0) return (tmp1 < 0);

  const int tmp2 = ((*bspline_v_) < (*rhs.bspline_v_));
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

  const bool tmp1 = ((*bspline_u_) == (*rhs.bspline_u_));
  if (tmp1 == false)
    return false;

  const bool tmp2 = ((*bspline_v_) == (*rhs.bspline_v_));
  if (tmp2 == false)
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
  bspline_u_->write(os);
  bspline_v_->write(os);
  object_to_stream(os, '\n');
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

  // Univariate B-splines
  bspline_u_ = new BSplineUniLR();
  bspline_u_->read(is);
  bspline_v_ = new BSplineUniLR();
  bspline_v_->read(is);

  coef_fixed_ = 0;
}


//==============================================================================
  void LRBSpline2D::read(istream& is, 
			 vector<std::unique_ptr<BSplineUniLR> >& bsplineuni_u,
			 int& left1,
			 vector<std::unique_ptr<BSplineUniLR> >& bsplineuni_v,
			 int& left2)
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

  // Univariate B-splines
  BSplineUniLR *tmpu = new BSplineUniLR();
  tmpu->read(is);
  tmpu->setPardir(1);

  bool found1 = BSplineUniUtils::identify_bsplineuni(tmpu, bsplineuni_u, left1);
  if (found1)
    delete tmpu;
  else
    BSplineUniUtils::insert_univariate(bsplineuni_u, tmpu, left1);
  bspline_u_ = bsplineuni_u[left1].get();
  bspline_u_->incrCount();
  
  BSplineUniLR *tmpv = new BSplineUniLR();
  tmpv->read(is);
  tmpv->setPardir(2);

  bool found2 = BSplineUniUtils::identify_bsplineuni(tmpv, bsplineuni_v, left2);
  if (found2)
    delete tmpv;
  else
    BSplineUniUtils::insert_univariate(bsplineuni_v, tmpv, left2);
  bspline_v_ = bsplineuni_v[left2].get();
  bspline_v_->incrCount();
  
  coef_fixed_ = 0;
}

//==============================================================================
double LRBSpline2D::evalBasisFunc(double u, 
				  double v) const
//==============================================================================
{
  return 
    bspline_u_->evalBasisFunc(u)*bspline_v_->evalBasisFunc(v);
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
  double bval1 = bspline_u_->evalBasisFunction(u, u_deriv, u_at_end);
  double bval2 = bspline_v_->evalBasisFunction(v, v_deriv, v_at_end);
  return bval1*bval2;
}


//==============================================================================
  void LRBSpline2D::evalder_add(double u, double v, 
				int deriv,
				Point der[],
				bool u_at_end, bool v_at_end) const
//==============================================================================
{
  double eps = 1.0e-12;
   u_at_end = (u >= umax()-eps);
   v_at_end = (v >= vmax()-eps);

   deriv = std::min(MAX_DER, deriv);
   double dd[2*MAX_DER+2];
   double *bder1 = dd;
   double *bder2 = dd+deriv+1;
   bspline_u_->evalBasisFunctions(u, deriv, bder1, u_at_end);
   bspline_v_->evalBasisFunctions(v, deriv, bder2, v_at_end);

   int ki, kj, kr, kh;
   // vector<double> bb((deriv+1)*(deriv+2)/2);
   // for (kj=0, kr=0; kj<=deriv; ++kj)
   //   for (ki=0; ki<=kj; ++ki, ++kr)
   //     bb[kr] = evalBasisFunction(u, v, kj-ki, ki, u_at_end, v_at_end);

   if (rational_)
     {
       int dim = coef_times_gamma_.dimension();
       int nmb = (deriv+1)*(deriv+2)/2;
       double tmp[(int)((MAX_DER+1)*(MAX_DER+2)*(MAX_DIM+1))];
       double *tmpder = tmp;
       double val;
       Point tmppt(dim);
       kh = 0;
       for (ki=0; ki<=deriv; ++ki)
	 for (kj=0; kj<=ki; ++kj, ++kh)
	   {
	     val = weight_*bder1[ki-kj]*bder2[kj];
	     for (kr=0; kr<dim; ++kr)
	       tmpder[kh*(dim+1)+kr] = coef_times_gamma_[kr]*val;
	     tmpder[kh*(dim+1)+dim] = val;
	   }
       double *tmpder2 = tmpder+nmb*(dim+1);
       SplineUtils::surface_ratder(tmpder, dim, deriv, tmpder2);
       for (kh=0; kh<nmb; ++kh)
	 {
	   for (kr=0; kr<dim; ++kr)
	     tmppt[kr] = tmpder2[kh*dim+kr];
	   der[kh] = tmppt;
	 }
     }
   else
     {
       kh = 0;
       for (ki=0; ki<=deriv; ++ki)
	 for (kj=0; kj<=ki; ++kj, ++kh)
	   {
	     der[kh] += coef_times_gamma_*bder1[ki-kj]*bder2[kj];
	     // Point bb2 = eval(u, v, ki-kj, kj, u_at_end, v_at_end);
	     // int stop_break = 1;
	   }
     }


   // return 
  //   compute_univariate_spline(degree(XFIXED), u, kvec(XFIXED), mesh_->knotsBegin(XFIXED), 
  // 			      u_deriv, u_at_end) *
  //   compute_univariate_spline(degree(YFIXED), v, kvec(YFIXED), mesh_->knotsBegin(YFIXED),  
  // 			      v_deriv, v_at_end);
}


//==============================================================================
void LRBSpline2D::evalBasisGridDer(int nmb_der, const vector<double>& par1, 
				   const vector<double>& par2, 
				   vector<double>& derivs) const
//==============================================================================
{
  if (nmb_der == 0)
    return;  // No derivatives to compute
  nmb_der = std::max(nmb_der, 3); // At most third order derivatives
  // Allocate stratch
  int nmb1 = (int)(par1.size());
  int nmb2 = (int)(par2.size());
  int nmb_part_der = 2;
  if (nmb_der > 1)
    nmb_part_der += 3;
  if (nmb_der > 2)
    nmb_part_der += 4;  
  derivs.resize(nmb_part_der*nmb1*nmb2);
  vector<double> ebder1((nmb_der+1)*nmb1);
  vector<double> ebder2((nmb_der+1)*nmb2);

  // Compute derivatives of univariate basis
  int ki, kj;
  for (ki=0; ki<nmb1; ++ki)
   bspline_u_->evalBasisFunctions(par1[ki], nmb_der, 
				  &ebder1[ki*(nmb_der+1)], false);
  for (ki=0; ki<nmb2; ++ki)
   bspline_v_->evalBasisFunctions(par2[ki], nmb_der, 
				  &ebder2[ki*(nmb_der+1)], false);
  // for (ki=0; ki<nmb1; ++ki)
  //   {
  //     // For the time being. Should be made more effective
  //     for (int kii=0; kii<=nmb_der; ++kii)
  // 	{
  // 	  ebder1[ki*(nmb_der+1)+kii] = compute_univariate_spline(degree(XFIXED), 
  // 								 par1[ki], kvec(XFIXED), 
  // 								 mesh_->knotsBegin(XFIXED), 
  // 								 kii, false);
  // 	}
  //   }

  // for (ki=0; ki<nmb2; ++ki)
  //   {
  //     // For the time being. Should be made more effective
  //     for (int kii=0; kii<=nmb_der; ++kii)
  // 	{
  // 	  ebder2[ki*(nmb_der+1)+kii] = compute_univariate_spline(degree(YFIXED), 
  // 								 par2[ki], kvec(YFIXED), 
  // 								 mesh_->knotsBegin(YFIXED), 
  // 								 kii, false);
  // 	}
  //   }

  // Combine univariate results
  // NOTE that rational functions are NOT handled
  int kr;
  for (kj=0; kj<nmb2; ++kj)
    for (ki=0; ki<nmb1; ++ki)
      {
	derivs[kj*nmb1+ki] = 
	  gamma_*ebder1[ki*(nmb_der+1)+1]*ebder2[kj*(nmb_der+1)]; // du
	derivs[(nmb2+kj)*nmb1+ki] = 
	  gamma_*ebder1[ki*(nmb_der+1)]*ebder2[kj*(nmb_der+1)+1]; // dv
	if (nmb_der > 1)
	  {
	    derivs[(2*nmb2+kj)*nmb1+ki] = 
	      gamma_*ebder1[ki*(nmb_der+1)+2]*ebder2[kj*(nmb_der+1)]; // duu
	    derivs[(3*nmb2+kj)*nmb1+ki] = 
	      gamma_*ebder1[ki*(nmb_der+1)+1]*ebder2[kj*(nmb_der+1)+1]; // duv
	    derivs[(4*nmb2+kj)*nmb1+ki] = 
	      gamma_*ebder1[ki*(nmb_der+1)]*ebder2[kj*(nmb_der+1)+2]; // dvv
	    if (nmb_der > 2)
	      {
		derivs[(5*nmb2+kj)*nmb1+ki] = 
		  gamma_*ebder1[ki*(nmb_der+1)+3]*ebder2[kj*(nmb_der+1)]; // duuu
		derivs[(6*nmb2+kj)*nmb1+ki] = 
		  gamma_*ebder1[ki*(nmb_der+1)+2]*ebder2[kj*(nmb_der+1)+1]; // duuv
		derivs[(7*nmb2+kj)*nmb1+ki] = 
		  gamma_*ebder1[ki*(nmb_der+1)+1]*ebder2[kj*(nmb_der+1)+2]; // duvv
		derivs[(8*nmb2+kj)*nmb1+ki] = 
		  gamma_*ebder1[ki*(nmb_der+1)]*ebder2[kj*(nmb_der+1)+3]; // dvvv
	      }
	  }
      }
}

//==============================================================================
  void LRBSpline2D::evalBasisLineDer(int nmb_der, Direction2D d, 
				     const vector<double>& parval, 
				     vector<double>& derivs) const
//==============================================================================
{
  if (nmb_der == 0)
    return;  // No derivatives to compute
  nmb_der = std::max(nmb_der, 3); // At most third order derivatives
  // Allocate stratch
  int nmb = (int)(parval.size());
  derivs.resize(nmb_der*nmb);
  vector<double> ebder((nmb_der+1)*nmb);

  // Compute derivatives of univariate basis
  int ki;
  for (ki=0; ki<nmb; ++ki)
    {
      if (d == XFIXED)
	bspline_u_->evalBasisFunctions(parval[ki], nmb_der, 
				       &ebder[ki*(nmb_der+1)]);
      else
	bspline_v_->evalBasisFunctions(parval[ki], nmb_der, 
				       &ebder[ki*(nmb_der+1)]);
      // // For the time being. Should be made more effective
      // for (int kii=0; kii<=nmb_der; ++kii)
      // 	{
      // 	  ebder[ki*(nmb_der+1)+kii] = (d == XFIXED) ?
      // 	    bspline_u_->evalBasisFunction(parval[ki], kii, false) :
      // 	    bspline_v_->evalBasisFunction(parval[ki], kii, false);
      // 	}
    }

  // Multiply with weight
  // NOTE that rational functions are NOT handled
  int kr;
  for (ki=0; ki<nmb; ++ki)
      {
	derivs[ki] = gamma_*ebder[ki*(nmb_der+1)+1]; // dt
	if (nmb_der > 1)
	  {
	    derivs[nmb+ki] = gamma_*ebder[ki*(nmb_der+1)+2]; // dtt
	    if (nmb_der > 2)
	      {
		derivs[2*nmb+ki] = gamma_*ebder[ki*(nmb_der+1)+3]; // dttt
	      }
	  }
      }
}

//==============================================================================
int LRBSpline2D::endmult_u(bool atstart) const
//==============================================================================
{
  return bspline_u_->endmult(atstart);
}

//==============================================================================
int LRBSpline2D::endmult_v(bool atstart) const
//==============================================================================
{
  return bspline_v_->endmult(atstart);
}

//==============================================================================
Point LRBSpline2D::getGrevilleParameter() const
{
  double upar = bspline_u_->getGrevilleParameter();
  double vpar = bspline_v_->getGrevilleParameter();
  Point greville(upar, vpar);
  return greville;
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
bool LRBSpline2D::hasSupportedElement(Element2D *el) 
//==============================================================================
{
  for (size_t i=0; i<support_.size(); i++) 
    {
    if(el == support_[i]) 
      return true;
    }
  return false;
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
void LRBSpline2D::reverseParameterDirection(bool dir_is_u)
//==============================================================================
{
  if (dir_is_u)
    bspline_u_->reverseParameterDirection();
  else
    bspline_v_->reverseParameterDirection();
}


//==============================================================================
void LRBSpline2D::swapParameterDirection()
//==============================================================================
{
  std::swap(bspline_u_, bspline_v_);
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

  vector<vector<double> > coefs(deg+1);
  for (int i = 0; i < deg + 2; ++i)
    {
    knots.push_back(slope * (knotval(d,knots_int[i]) - start));

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
    }
  return coefs[0];
}


}; // end namespace Go
