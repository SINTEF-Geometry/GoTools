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
#include "sislP.h"
#include <set>
#include <algorithm>

//#define DEBUG
//#define DEBUG_PROJ

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
  nest_level_ = rhs.nest_level_; // To be computed?
  overload_ = rhs.overload_;
  visited_ = rhs.visited_;
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

  nest_level_ = -1;
  coef_fixed_ = 0;
  overload_ = false;
  visited_ = false;
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
  
  nest_level_ = -1;
  coef_fixed_ = 0;
  overload_ = false;
  visited_ = false;
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
  int dim = coef_times_gamma_.dimension();
  ALWAYS_ERROR_IF(rational_ && (der[0].dimension() != dim + 1), "Dimension mismatch for rational");
  
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
   if (rational_)
     {
       kh = 0;
       for (ki=0; ki<=deriv; ++ki)
	 for (kj=0; kj<=ki; ++kj, ++kh)
	   {
	     for (kr=0; kr<dim; ++kr)
	       der[kh][kr] +=  weight_*coef_times_gamma_[kr]*bder1[ki-kj]*bder2[kj];
	     der[kh][dim] += weight_*bder1[ki-kj]*bder2[kj];
	   }
     }
   else
     {
       kh = 0;
       for (ki=0; ki<=deriv; ++ki)
	 for (kj=0; kj<=ki; ++kj, ++kh)
	   {
	     der[kh] += coef_times_gamma_*bder1[ki-kj]*bder2[kj];
	   }
     }
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
vector<double> LRBSpline2D::kvec_val(Direction2D d) const
//==============================================================================
{
   const Mesh2D* mesh = dynamic_cast<const Mesh2D*>(bspline_u_->getMesh());
  vector<double> kval(kvec(d).size());
  for (size_t ki=0; ki<kval.size(); ++ki)
    kval[ki] = mesh->kval(d, kvec(d)[ki]);
  return kval;
}

//==============================================================================
Point LRBSpline2D::getGrevilleParameter() const
//==============================================================================
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
bool LRBSpline2D::overlaps(LRBSpline2D *bsp) const
//==============================================================================
{
  // Does it make sense to include equality?
  if (bsp->umin() >= umax())
    return false;
  if (bsp->umax() <= umin())
    return false;
  if (bsp->vmin() >= vmax())
    return false;
  if (bsp->vmax() <= vmin())
    return false;
  
  return true;
}

//==============================================================================
bool LRBSpline2D::covers(double domain[]) const
//==============================================================================
{
  if (!overlaps(domain))
    return false;
  if (domain[0] < umin())
    return false;
  if (domain[1] > umax())
    return false;
  if (domain[2] < vmin())
    return false;
  if (domain[3] > vmax())
    return false;
  
  return true;
}

//==============================================================================
bool LRBSpline2D::covers(LRBSpline2D *bsp) const
//==============================================================================
{
  if (!overlaps(bsp))
    return false;
  if (bsp->umin() < umin())
    return false;
  if (bsp->umax() > umax())
    return false;
  if (bsp->vmin() < vmin())
    return false;
  if (bsp->vmax() > vmax())
    return false;

  if (bsp->umin() == umin() && bsp->endmult_u(true) > endmult_u(true))
    return false;
  if (bsp->umax() == umax() && bsp->endmult_u(false) > endmult_u(false))
    return false;
  if (bsp->vmin() == vmin() && bsp->endmult_v(true) > endmult_v(true))
    return false;
  if (bsp->vmax() == vmax() && bsp->endmult_v(false) > endmult_v(false))
    return false;
  
  return true;
}
 
//==============================================================================
void LRBSpline2D::computeNestLevel()
//==============================================================================
{
  if (nest_level_ >= 0)
    return;   // Assumes the recorded value is correct

  visited_ = true;
  
  // Collect potential ancestors
  set<LRBSpline2D*> cand;
  for (auto el=support_.begin(); el!=support_.end(); ++el)
    {
      for (auto bsp=(*el)->supportBegin(); bsp!=(*el)->supportEnd(); ++bsp)
	if ((*bsp) != this)
	  cand.insert(*bsp);
    }

  for (auto bsp=cand.begin(); bsp!=cand.end(); ++bsp)
    {
      if ((*bsp)->visited())
	continue;
      if (!(*bsp)->hasNestLevel())
	(*bsp)->computeNestLevel();
      if ((*bsp)->covers(this))
	{
	  int count = (*bsp)->getNestLevel();
	  nest_level_ = std::max(nest_level_, count + 1);
	}
    }

  if (nest_level_ < 0)
    nest_level_ = 0;

  for (auto bsp=cand.begin(); bsp!=cand.end(); ++bsp)
    {
    if (covers(*bsp))
      {
	int level = (*bsp)->getNestLevel();
	if (level >= 0)
	  (*bsp)->setNestLevel(std::max(level, nest_level_+1));
      }
    }
      
  visited_ = false;
}

//==============================================================================
void LRBSpline2D::getOverlapping(vector<LRBSpline2D*>& overlap)
//==============================================================================
{
  set<LRBSpline2D*> cand;
  for (auto el=support_.begin(); el!=support_.end(); ++el)
    {
      for (auto bsp=(*el)->supportBegin(); bsp!=(*el)->supportEnd(); ++bsp)
	if ((*bsp) != this)
	  cand.insert(*bsp);
    }

  for (auto bsp=cand.begin(); bsp!=cand.end(); ++bsp)
    {
      if ((*bsp)->covers(this))
	overlap.push_back(*bsp);
    }
 
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
bool LRBSpline2D::adaptProjCoef(Point& coef)
//==============================================================================
{
  if (nest_level_ == 0)
    return true;

  int cdim = coef.dimension();
  
  // Collect ancestors
  set<LRBSpline2D*> ancest0;
  for (auto el=support_.begin(); el!=support_.end(); ++el)
    {
      for (auto bsp=(*el)->supportBegin(); bsp!=(*el)->supportEnd(); ++bsp)
	if ((*bsp) != this)
	  {
	    int level = (*bsp)->getNestLevel();
	    if (level < nest_level_ && (*bsp)->covers(this))
	      ancest0.insert(*bsp);
	  }
    }

#ifdef DEBUG_PROJ
   std::cout << "Nesting level: " << nest_level_ << ", scale factor: " << gamma_ << std::endl;
  std::cout << "Knots curr: [";
  vector<int> kvec_u1 = bspline_u_->kvec();
  vector<int> kvec_v1 = bspline_v_->kvec();
   for (size_t kj=0; kj<kvec_u1.size(); ++kj)
    std::cout << knotval(XFIXED, kvec_u1[kj]) << ", ";
  std::cout << "]x[";
  for (size_t kj=0; kj<kvec_v1.size(); ++kj)
    std::cout << knotval(YFIXED, kvec_v1[kj]) << ", ";
  std::cout << "]" << std::endl;
#endif
    
  vector<LRBSpline2D*> ancest(ancest0.begin(), ancest0.end());
#ifdef DEBUG_PROJ
   if (ancest.size() > 1)
    std::cout << "Number of ancestors: " << ancest.size() << std::endl;
#endif
  double tmp = 0.0;
  for (size_t ki=0; ki<ancest.size(); ++ki)
    {
      double weight = nestingWeight(ancest[ki]);
      Point coefgamma = ancest[ki]->coefTimesGamma();
      if (coefgamma.dimension() > cdim)
	{
	  Point tmp(coefgamma.begin(), coefgamma.begin()+cdim);
	  coefgamma = tmp;
	}
      coef -= weight*coefgamma;
      double gamma = ancest[ki]->gamma();
      tmp += weight*gamma;
    }
  double tmp2 = (1.0 - tmp)/gamma_;
#ifdef DEBUG_PROJ
  std::cout << "Projection coefficient adaption, nesting depth: " << nest_level_ << ", error: " << 1.0-tmp2 << std::endl;
#endif
  if (fabs(tmp2-1.0) > 1.0e-4)
    {
#ifdef DEBUG_PROJ
      std::cout << "Invariant: " << tmp2 << std::endl;
#endif
      setNestLevel(-1);
      computeNestLevel();
      return false;
    }
  coef /= gamma_;
  return true;
}

  struct knotwgt
  {
    vector<int> kvec_;
    double alpha_;

    knotwgt(vector<int> kvec, double alpha)
    {
      kvec_ = kvec;
      alpha_ = alpha;
    }
  };
  
//==============================================================================
double LRBSpline2D::nestingWeight(LRBSpline2D* other)
//==============================================================================
{
  vector<int> kvec_u1 = bspline_u_->kvec();
  vector<int> kvec_v1 = bspline_v_->kvec();
  vector<int> kvec_u2_0 = other->kvec(XFIXED);
  vector<int> kvec_v2_0 = other->kvec(YFIXED);
  vector<knotwgt> kvec_u2;
  kvec_u2.push_back(knotwgt(kvec_u2_0, 1.0));
  vector<knotwgt> kvec_v2;
  kvec_v2.push_back(knotwgt(kvec_v2_0, 1.0));

#ifdef DEBUG_PROJ
  std::cout << "Knots ancestor: [";
  for (size_t kj=0; kj<kvec_u2_0.size(); ++kj)
    std::cout << knotval(XFIXED, kvec_u2_0[kj]) << ", ";
  std::cout << "]x[";
  for (size_t kj=0; kj<kvec_v2_0.size(); ++kj)
    std::cout << knotval(YFIXED, kvec_v2_0[kj]) << ", ";
  std::cout << "]" << std::endl;
  std::cout << "Nesting level: " << other->nest_level_ << ", scale factor: " << other->gamma_ << std::endl;
#endif
  vector<int> diff1, diff2;
  std::set_difference(kvec_u1.begin(), kvec_u1.end(), kvec_u2_0.begin(),
		      kvec_u2_0.end(), std::back_inserter(diff1));
  std::set_difference(kvec_v1.begin(), kvec_v1.end(), kvec_v2_0.begin(),
		      kvec_v2_0.end(), std::back_inserter(diff2));

  const Mesh2D *mesh = getMesh();

#ifdef DEBUG_PROJ
  std::cout << "Knots in 1. parameter direction: " << diff1.size() << std::endl;
  std::cout << "Knots in 2. parameter direction: " << diff2.size() << std::endl;
 #endif

  int k1 = kvec_u1[0]; 
  int k2 = kvec_u1[kvec_u1.size()-1]; 
  for (size_t ki=0; ki<diff1.size(); ++ki)
    {
      double val = knotval(XFIXED, diff1[ki]);
      for (size_t kj=0; kj<kvec_u2.size(); ++kj)
	{
	  int ks = kvec_u2[kj].kvec_.size()-1;
	  if (k1 < kvec_u2[kj].kvec_[0] || k2 > kvec_u2[kj].kvec_[ks])
	    continue;
	  double y1 = mesh->kval(XFIXED, kvec_u2[kj].kvec_[0]);
	  double y2 = mesh->kval(XFIXED, kvec_u2[kj].kvec_[1]);
	  double y3 = mesh->kval(XFIXED, kvec_u2[kj].kvec_[ks-1]);
	  double y4 = mesh->kval(XFIXED, kvec_u2[kj].kvec_[ks]);
	  double a1 = (val >= y3) ? 1.0 : (val - y1)/(y3 - y1);
	  double a2 = (val <= y2) ? 1.0 : (y4 - val)/(y4 - y2);
	  size_t kr;
	  for (kr=1; kr<kvec_u2[kj].kvec_.size(); ++kr)
	    if (diff1[ki] > kvec_u2[kj].kvec_[kr-1] && diff1[ki] < kvec_u2[kj].kvec_[kr])
	      break;
	  if (kr == kvec_u2[kj].kvec_.size())
	    kvec_u2[kj].kvec_.push_back(diff1[ki]);
	  else
	    kvec_u2[kj].kvec_.insert(kvec_u2[kj].kvec_.begin()+kr, diff1[ki]);
	  vector<int> kv(kvec_u2[kj].kvec_.begin()+1, kvec_u2[kj].kvec_.end());
	  double alp = kvec_u2[kj].alpha_*a2;
	  kvec_u2[kj].kvec_.pop_back();
	  kvec_u2[kj].alpha_ *= a1;
	  kvec_u2.insert(kvec_u2.begin()+kj+1, knotwgt(kv, alp));
	  ++kj;
	}
      
      for (size_t kj=0; kj<kvec_u2.size();)
	{
	  int ks = kvec_u2[kj].kvec_.size()-1;
	  if (k1 < kvec_u2[kj].kvec_[0] || k2 > kvec_u2[kj].kvec_[ks])
	    kvec_u2.erase(kvec_u2.begin()+kj);
	  else
	    ++kj;
	}

      for (size_t kj=1; kj<kvec_u2.size();)
	{
	  if (std::equal(kvec_u2[kj-1].kvec_.begin(), kvec_u2[kj-1].kvec_.end(),
			 kvec_u2[kj].kvec_.begin()))
	    {
	      kvec_u2[kj-1].alpha_ += kvec_u2[kj].alpha_;
	      kvec_u2.erase(kvec_u2.begin()+kj);
	    }
	  else
	    ++kj;
	}
    }
  
  k1 = kvec_v1[0]; 
  k2 = kvec_v1[kvec_v1.size()-1]; 
  for (size_t ki=0; ki<diff2.size(); ++ki)
    {
      double val = knotval(YFIXED, diff2[ki]);
      for (size_t kj=0; kj<kvec_v2.size(); ++kj)
	{
	  int ks = kvec_v2[kj].kvec_.size()-1;
	  if (k1 < kvec_v2[kj].kvec_[0] || k2 > kvec_v2[kj].kvec_[ks])
	    continue;
	  double y1 = mesh->kval(YFIXED, kvec_v2[kj].kvec_[0]);
	  double y2 = mesh->kval(YFIXED, kvec_v2[kj].kvec_[1]);
	  double y3 = mesh->kval(YFIXED, kvec_v2[kj].kvec_[ks-1]);
	  double y4 = mesh->kval(YFIXED, kvec_v2[kj].kvec_[ks]);
	  double a1 = (val >= y3) ? 1.0 : (val - y1)/(y3 - y1);
	  double a2 = (val <= y2) ? 1.0 : (y4 - val)/(y4 - y2);
	  size_t kr;
	  for (kr=1; kr<kvec_v2[kj].kvec_.size(); ++kr)
	    if (diff2[ki] > kvec_v2[kj].kvec_[kr-1] && diff2[ki] < kvec_v2[kj].kvec_[kr])
	      break;
	  if (kr == kvec_v2[kj].kvec_.size())
	    kvec_v2[kj].kvec_.push_back(diff2[ki]);
	  else
	    kvec_v2[kj].kvec_.insert(kvec_v2[kj].kvec_.begin()+kr, diff2[ki]);
	  vector<int> kv(kvec_v2[kj].kvec_.begin()+1, kvec_v2[kj].kvec_.end());
	  double alp = kvec_v2[kj].alpha_*a2;
	  kvec_v2[kj].kvec_.pop_back();
	  kvec_v2[kj].alpha_ *= a1;
	  kvec_v2.insert(kvec_v2.begin()+kj+1, knotwgt(kv, alp));
	  ++kj;
	}
      
      for (size_t kj=0; kj<kvec_v2.size();)
	{
	  int ks = kvec_v2[kj].kvec_.size()-1;
	  if (k1 < kvec_v2[kj].kvec_[0] || k2 > kvec_v2[kj].kvec_[ks])
	    kvec_v2.erase(kvec_v2.begin()+kj);
	  else
	    ++kj;
	}

      for (size_t kj=1; kj<kvec_v2.size();)
	{
	  if (std::equal(kvec_v2[kj-1].kvec_.begin(), kvec_v2[kj-1].kvec_.end(),
			 kvec_v2[kj].kvec_.begin()))
	    {
	      kvec_v2[kj-1].alpha_ += kvec_v2[kj].alpha_;
	      kvec_v2.erase(kvec_v2.begin()+kj);
	    }
	  else
	    ++kj;
	}
    }

  double alpha = 1.0;
  for (size_t kj=0; kj<kvec_u2.size(); ++kj)
    if (std::equal(kvec_u1.begin(), kvec_u1.end(), &kvec_u2[kj].kvec_[0]))
	alpha *= kvec_u2[kj].alpha_;
  for (size_t kj=0; kj<kvec_v2.size(); ++kj)
    if (std::equal(kvec_v1.begin(), kvec_v1.end(), &kvec_v2[kj].kvec_[0]))
	alpha *= kvec_v2[kj].alpha_;

#ifdef DEBUG_PROJ
   std::cout << "Weight: " << alpha << std::endl;
#endif
  return alpha;
}


//==============================================================================
bool LRBSpline2D::checkOverload()
//==============================================================================
{
  bool overload = true;
  for (size_t ki=0; ki<support_.size(); ++ki)
    if (!support_[ki]->getOverload())
      {
       overload = false;
       break;
      }

  overload_ = overload;
  return overload;
}


//==============================================================================
void LRBSpline2D::swapParameterDirection()
//==============================================================================
{
  bspline_u_->setPardir(2);
  bspline_v_->setPardir(1);
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
    //    {
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
  //    }
  return coefs[0];
}


}; // end namespace Go
