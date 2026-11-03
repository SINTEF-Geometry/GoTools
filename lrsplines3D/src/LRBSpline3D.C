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

#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
#include "GoTools/utils/checks.h"
#include "GoTools/utils/StreamUtils.h"
#include <set>

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
LRBSpline3D::LRBSpline3D(const LRBSpline3D& rhs)
//==============================================================================
{
  coef_fixed_ = rhs.coef_fixed_;
  coef_times_gamma_ = rhs.coef_times_gamma_;
  gamma_ = rhs.gamma_;
  bspline_u_ = rhs.bspline_u_;
  bspline_u_->incrCount();  // Initial count is zero
  bspline_v_ = rhs.bspline_v_;
  bspline_v_->incrCount();
  bspline_w_ = rhs.bspline_w_;
  bspline_w_->incrCount();
   rational_ = rhs.rational_;
  // don't copy the support
  weight_ = rhs.weight_;
  nest_level_ = rhs.nest_level_; // To be computed?
  visited_ = rhs.visited_;

}
  //==============================================================================
  void LRBSpline3D::write(std::ostream& os) const
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
    bspline_w_->write(os);
  }
  
  //==============================================================================
  void LRBSpline3D::read(std::istream& is)
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
  object_from_stream(is, weight_);
  // if (gamma_ < 1.0)
  // {
  //     MESSAGE("DEBUGGING: Changing gamma from " << gamma_ << " to 1.0!");
  //     coef_times_gamma_ /= gamma_;
  //     gamma_ = 1.0;
  // }
  // Univariate B-splines
  bspline_u_ = new BSplineUniLR();
  bspline_u_->read(is);
  bspline_v_ = new BSplineUniLR();
  bspline_v_->read(is);
  bspline_w_ = new BSplineUniLR();
  bspline_w_->read(is);

  nest_level_ = -1;
  coef_fixed_ = 0;
  visited_ = false;
  }
  
//==============================================================================
  void LRBSpline3D::read(istream& is, 
			 vector<std::unique_ptr<BSplineUniLR> >& bsplineuni_u,
			 int& left1,
			 vector<std::unique_ptr<BSplineUniLR> >& bsplineuni_v,
			 int& left2,
			 vector<std::unique_ptr<BSplineUniLR> >& bsplineuni_w,
			 int& left3)
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
  
  BSplineUniLR *tmpw = new BSplineUniLR();
  tmpw->read(is);
  tmpw->setPardir(3);

  bool found3 = BSplineUniUtils::identify_bsplineuni(tmpw, bsplineuni_w, left3);
  if (found3)
    delete tmpw;
  else
    BSplineUniUtils::insert_univariate(bsplineuni_w, tmpw, left3);
  bspline_w_ = bsplineuni_w[left3].get();
  bspline_w_->incrCount();
  
  nest_level_ = -1;
  coef_fixed_ = 0;
  visited_ = false;
}

  //==============================================================================
  double LRBSpline3D::evalBasisFunc(double u,
                                    double v,
                                    double w) const
  //==============================================================================
  {
    return
      bspline_u_->evalBasisFunc(u)*bspline_v_->evalBasisFunc(v)*
      bspline_w_->evalBasisFunc(w);
  }


  //==============================================================================
  double LRBSpline3D::evalBasisFunction(double u, double v, double w,
                                        int u_deriv, int v_deriv, int w_deriv,
                                        bool u_at_end, bool v_at_end, bool w_at_end) const

  //==============================================================================
  {
    // return
    //   bspline_u_->evalBasisFunction(u, u_deriv, u_at_end)*
    //   bspline_v_->evalBasisFunction(v, v_deriv, v_at_end)*
    //   bspline_w_->evalBasisFunction(w, w_deriv, w_at_end);
    double bval1 = bspline_u_->evalBasisFunction(u, u_deriv, u_at_end);
    double bval2 = bspline_v_->evalBasisFunction(v, v_deriv, v_at_end);
    double bval3 = bspline_w_->evalBasisFunction(w, w_deriv, w_at_end);
    return bval1*bval2*bval3;
  }

//==============================================================================
  void LRBSpline3D::evalder_add(double u, double v, double w, 
				int deriv,
				Point der[],
				bool u_at_end, bool v_at_end, 
				bool w_at_end) const
//==============================================================================
{
  int dim = coef_times_gamma_.dimension();
  ALWAYS_ERROR_IF(rational_ && (der[0].dimension() != dim + 1), "Dimension mismatch for rational");
  
  double eps = 1.0e-12;
  u_at_end = (u >= umax()-eps);
  v_at_end = (v >= vmax()-eps);
  w_at_end = (w >= wmax()-eps);

   deriv = std::min(MAX_DER, deriv);
   double dd[3*MAX_DER+3];
   double *bder1 = dd;
   double *bder2 = dd+deriv+1;
   double *bder3 = bder2+deriv+1;
   bspline_u_->evalBasisFunctions(u, deriv, bder1, u_at_end);
   bspline_v_->evalBasisFunctions(v, deriv, bder2, v_at_end);
   bspline_w_->evalBasisFunctions(w, deriv, bder3, w_at_end);

   int ki, kj, kk, kr, kh;
   if (rational_)
     {
      kh = 0;
       for (ki=0; ki<=deriv; ++ki)
	 for (kj=0; kj<=ki; ++kj)
	   for (kk=0; kk<=kj; ++kk, ++kh)
	     {
	     for (kr=0; kr<dim; ++kr)
	       der[kh][kr] += weight_*coef_times_gamma_[kr]*bder1[ki-kj]*
		 bder2[kj-kk]*bder3[kk];
	     der[kh][dim] += weight_*bder1[ki-kj]*bder2[kj-kk]*bder3[kk];
	     }
     }
   else
     {
       kh = 0;
       for (ki=0; ki<=deriv; ++ki)
	 for (kj=0; kj<=ki; ++kj)
	   for (kk=0; kk<=kj; ++kk, ++kh)
	     {
	       der[kh] += coef_times_gamma_*bder1[ki-kj]*
		 bder2[kj-kk]*bder3[kk];
	     }
     }
}


//==============================================================================
int LRBSpline3D::endmult_u(bool atstart) const
//==============================================================================
{
  return bspline_u_->endmult(atstart);
}


//==============================================================================
int LRBSpline3D::endmult_v(bool atstart) const
//==============================================================================
{
  return bspline_v_->endmult(atstart);
}


//==============================================================================
int LRBSpline3D::endmult_w(bool atstart) const
//==============================================================================
{
  return bspline_w_->endmult(atstart);

}

//==============================================================================
int LRBSpline3D::endmult(Direction3D dir, bool atstart) const
//==============================================================================
{
  return getUnivariate(dir)->endmult(atstart);

}

  //==============================================================================
  Point LRBSpline3D::getGrevilleParameter() const
  //==============================================================================
  {
    double upar = bspline_u_->getGrevilleParameter();
    double vpar = bspline_v_->getGrevilleParameter();
    double wpar = bspline_w_->getGrevilleParameter();
    Point greville(upar, vpar, wpar);
    return greville;
  }


  //==============================================================================
  bool LRBSpline3D::overlaps(Element3D *el) const
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
    if (el->wmin() >= wmax())
      return false;
    if (el->wmax() <= wmin())
      return false;

    return true;
  }

  // Operations related to the support of this B-spline
  //==============================================================================
  bool LRBSpline3D::overlaps(double domain[]) const
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
    if (domain[4] >= wmax())
      return false;
    if (domain[5] <= wmin())
      return false;

    return true;
  }

//==============================================================================
bool LRBSpline3D::overlaps(LRBSpline3D *bsp) const
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
  if (bsp->wmin() >= wmax())
    return false;
  if (bsp->wmax() <= wmin())
    return false;
  
  return true;
}

//==============================================================================
bool LRBSpline3D::covers(double domain[]) const
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
  if (domain[4] < wmin())
    return false;
  if (domain[5] > wmax())
    return false;
 
  return true;
}

//==============================================================================
bool LRBSpline3D::covers(LRBSpline3D *bsp) const
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
  if (bsp->wmin() < wmin())
    return false;
  if (bsp->wmax() > wmax())
    return false;

  if (bsp->umin() == umin() && bsp->endmult_u(true) > endmult_u(true))
    return false;
  if (bsp->umax() == umax() && bsp->endmult_u(false) > endmult_u(false))
    return false;
  if (bsp->vmin() == vmin() && bsp->endmult_v(true) > endmult_v(true))
    return false;
  if (bsp->vmax() == vmax() && bsp->endmult_v(false) > endmult_v(false))
    return false;
  if (bsp->wmin() == wmin() && bsp->endmult_w(true) > endmult_w(true))
    return false;
  if (bsp->wmax() == wmax() && bsp->endmult_w(false) > endmult_w(false))
    return false;
 
  return true;
}
 
//==============================================================================
void LRBSpline3D::getOverlapping(vector<LRBSpline3D*>& overlap)
//==============================================================================
{
  set<LRBSpline3D*> cand;
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
void LRBSpline3D::computeNestLevel()
//==============================================================================
{
  if (nest_level_ >= 0)
    return;   // Assumes the recorded value is correct

  visited_ = true;
  
  // Collect potential ancestors
  set<LRBSpline3D*> cand;
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
  bool LRBSpline3D::addSupport(Element3D *el)
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
  void LRBSpline3D::removeSupport(Element3D *el)
  //==============================================================================
  {
    auto it = std::find(support_.begin(), support_.end(), el);
    if (it != support_.end())
      {
    	*it = support_.back();
    	support_.pop_back();
      }
    // for (size_t i=0; i<support_.size(); i++) {
    //   if(el == support_[i]) {
    //     if (i < support_.size() - 1)
    //       {
    //         support_[i] = support_.back();
    //         support_[support_.size()-1] = NULL;
    //       }
    //     support_.pop_back();
    //     return;
    //   }
    // }
  }

  //==============================================================================
  std::vector<Element3D*>::iterator LRBSpline3D::supportedElementBegin()
  //==============================================================================
  {
    return support_.begin();
  }

  //==============================================================================
  std::vector<Element3D*>::iterator LRBSpline3D::supportedElementEnd()
  //==============================================================================
  {
    return support_.end();
  }

//==============================================================================
bool LRBSpline3D::adaptProjCoef(Point& coef)
//==============================================================================
{
  if (nest_level_ == 0)
    return true;

  int cdim = coef.dimension();
  
  // Collect ancestors
  set<LRBSpline3D*> ancest0;
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
  vector<int> kvec_w1 = bspline_w_->kvec();
   for (size_t kj=0; kj<kvec_u1.size(); ++kj)
    std::cout << knotval(XDIR, kvec_u1[kj]) << ", ";
  std::cout << "]x[";
  for (size_t kj=0; kj<kvec_v1.size(); ++kj)
    std::cout << knotval(YDIR, kvec_v1[kj]) << ", ";
  std::cout << "]x[";
  for (size_t kj=0; kj<kvec_w1.size(); ++kj)
    std::cout << knotval(YDIR, kvec_w1[kj]) << ", ";
  std::cout << "]" << std::endl;
#endif
    
  vector<LRBSpline3D*> ancest(ancest0.begin(), ancest0.end());
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
double LRBSpline3D::nestingWeight(LRBSpline3D* other)
//==============================================================================
{
  vector<int> kvec_u1 = bspline_u_->kvec();
  vector<int> kvec_v1 = bspline_v_->kvec();
  vector<int> kvec_w1 = bspline_w_->kvec();
  vector<int> kvec_u2_0 = other->kvec(XDIR);
  vector<int> kvec_v2_0 = other->kvec(YDIR);
  vector<int> kvec_w2_0 = other->kvec(ZDIR);
  vector<knotwgt> kvec_u2;
  kvec_u2.push_back(knotwgt(kvec_u2_0, 1.0));
  vector<knotwgt> kvec_v2;
  kvec_v2.push_back(knotwgt(kvec_v2_0, 1.0));
  vector<knotwgt> kvec_w2;
  kvec_w2.push_back(knotwgt(kvec_w2_0, 1.0));

#ifdef DEBUG_PROJ
  std::cout << "Knots ancestor: [";
  for (size_t kj=0; kj<kvec_u2_0.size(); ++kj)
    std::cout << knotval(XDIR, kvec_u2_0[kj]) << ", ";
  std::cout << "]x[";
  for (size_t kj=0; kj<kvec_v2_0.size(); ++kj)
    std::cout << knotval(YDIR, kvec_v2_0[kj]) << ", ";
  std::cout << "]x[";
  for (size_t kj=0; kj<kvec_w2_0.size(); ++kj)
    std::cout << knotval(ZDIR, kvec_w2_0[kj]) << ", ";
  std::cout << "]" << std::endl;
  std::cout << "Nesting level: " << other->nest_level_ << ", scale factor: " << other->gamma_ << std::endl;
#endif
  vector<int> diff1, diff2, diff3;
  std::set_difference(kvec_u1.begin(), kvec_u1.end(), kvec_u2_0.begin(),
		      kvec_u2_0.end(), std::back_inserter(diff1));
  std::set_difference(kvec_v1.begin(), kvec_v1.end(), kvec_v2_0.begin(),
		      kvec_v2_0.end(), std::back_inserter(diff2));
  std::set_difference(kvec_w1.begin(), kvec_w1.end(), kvec_w2_0.begin(),
		      kvec_w2_0.end(), std::back_inserter(diff3));

  const Mesh3D *mesh = getMesh();

#ifdef DEBUG_PROJ
  std::cout << "Knots in 1. parameter direction: " << diff1.size() << std::endl;
  std::cout << "Knots in 2. parameter direction: " << diff2.size() << std::endl;
  std::cout << "Knots in 3. parameter direction: " << diff3.size() << std::endl;
 #endif

  int k1 = kvec_u1[0]; 
  int k2 = kvec_u1[kvec_u1.size()-1]; 
  for (size_t ki=0; ki<diff1.size(); ++ki)
    {
      double val = knotval(XDIR, diff1[ki]);
      for (size_t kj=0; kj<kvec_u2.size(); ++kj)
	{
	  int ks = kvec_u2[kj].kvec_.size()-1;
	  if (k1 < kvec_u2[kj].kvec_[0] || k2 > kvec_u2[kj].kvec_[ks])
	    continue;
	  double y1 = mesh->kval(XDIR, kvec_u2[kj].kvec_[0]);
	  double y2 = mesh->kval(XDIR, kvec_u2[kj].kvec_[1]);
	  double y3 = mesh->kval(XDIR, kvec_u2[kj].kvec_[ks-1]);
	  double y4 = mesh->kval(XDIR, kvec_u2[kj].kvec_[ks]);
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
      double val = knotval(YDIR, diff2[ki]);
      for (size_t kj=0; kj<kvec_v2.size(); ++kj)
	{
	  int ks = kvec_v2[kj].kvec_.size()-1;
	  if (k1 < kvec_v2[kj].kvec_[0] || k2 > kvec_v2[kj].kvec_[ks])
	    continue;
	  double y1 = mesh->kval(YDIR, kvec_v2[kj].kvec_[0]);
	  double y2 = mesh->kval(YDIR, kvec_v2[kj].kvec_[1]);
	  double y3 = mesh->kval(YDIR, kvec_v2[kj].kvec_[ks-1]);
	  double y4 = mesh->kval(YDIR, kvec_v2[kj].kvec_[ks]);
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

  k1 = kvec_w1[0]; 
  k2 = kvec_w1[kvec_w1.size()-1]; 
  for (size_t ki=0; ki<diff3.size(); ++ki)
    {
      double val = knotval(ZDIR, diff3[ki]);
      for (size_t kj=0; kj<kvec_w2.size(); ++kj)
	{
	  int ks = kvec_w2[kj].kvec_.size()-1;
	  if (k1 < kvec_w2[kj].kvec_[0] || k2 > kvec_w2[kj].kvec_[ks])
	    continue;
	  double y1 = mesh->kval(ZDIR, kvec_w2[kj].kvec_[0]);
	  double y2 = mesh->kval(ZDIR, kvec_w2[kj].kvec_[1]);
	  double y3 = mesh->kval(ZDIR, kvec_w2[kj].kvec_[ks-1]);
	  double y4 = mesh->kval(ZDIR, kvec_w2[kj].kvec_[ks]);
	  double a1 = (val >= y3) ? 1.0 : (val - y1)/(y3 - y1);
	  double a2 = (val <= y2) ? 1.0 : (y4 - val)/(y4 - y2);
	  size_t kr;
	  for (kr=1; kr<kvec_w2[kj].kvec_.size(); ++kr)
	    if (diff3[ki] > kvec_w2[kj].kvec_[kr-1] && diff3[ki] < kvec_w2[kj].kvec_[kr])
	      break;
	  if (kr == kvec_w2[kj].kvec_.size())
	    kvec_w2[kj].kvec_.push_back(diff3[ki]);
	  else
	    kvec_w2[kj].kvec_.insert(kvec_w2[kj].kvec_.begin()+kr, diff3[ki]);
	  vector<int> kv(kvec_w2[kj].kvec_.begin()+1, kvec_w2[kj].kvec_.end());
	  double alp = kvec_w2[kj].alpha_*a2;
	  kvec_w2[kj].kvec_.pop_back();
	  kvec_w2[kj].alpha_ *= a1;
	  kvec_w2.insert(kvec_w2.begin()+kj+1, knotwgt(kv, alp));
	  ++kj;
	}
      
      for (size_t kj=0; kj<kvec_w2.size();)
	{
	  int ks = kvec_w2[kj].kvec_.size()-1;
	  if (k1 < kvec_w2[kj].kvec_[0] || k2 > kvec_w2[kj].kvec_[ks])
	    kvec_w2.erase(kvec_w2.begin()+kj);
	  else
	    ++kj;
	}

      for (size_t kj=1; kj<kvec_w2.size();)
	{
	  if (std::equal(kvec_w2[kj-1].kvec_.begin(), kvec_w2[kj-1].kvec_.end(),
			 kvec_w2[kj].kvec_.begin()))
	    {
	      kvec_w2[kj-1].alpha_ += kvec_w2[kj].alpha_;
	      kvec_w2.erase(kvec_w2.begin()+kj);
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
  for (size_t kj=0; kj<kvec_w2.size(); ++kj)
    if (std::equal(kvec_w1.begin(), kvec_w1.end(), &kvec_w2[kj].kvec_[0]))
	alpha *= kvec_w2[kj].alpha_;

#ifdef DEBUG_PROJ
   std::cout << "Weight: " << alpha << std::endl;
#endif
  return alpha;
}

#if 0
  //==============================================================================
  std::vector<Element3D*> LRBSpline3D::getExtendedSupport()
  //==============================================================================
  {
    MESSAGE("(): Not implemented.");
    throw;
  }

  //==============================================================================
  std::vector<Element3D*> LRBSpline3D::getMinimalExtendedSupport()
  //==============================================================================
  {
    MESSAGE("(): Not implemented.");
    throw;
  }
#endif
 
  //==============================================================================
  bool LRBSpline3D::operator<(const LRBSpline3D& rhs) const
  //==============================================================================
  {
    const int tmp1 = ((*bspline_u_) < (*rhs.bspline_u_));
    if (tmp1 != 0) return (tmp1 < 0);

    const int tmp2 = ((*bspline_v_) < (*rhs.bspline_v_));
    if (tmp2 != 0) return (tmp2 < 0);
    
    const int tmp3 = ((*bspline_w_) < (*rhs.bspline_w_));
    if (tmp3 != 0) return (tmp3 < 0);
    
    const int tmp4 = compare_seq(coef_times_gamma_.begin(), coef_times_gamma_.end(), 
				 rhs.coef_times_gamma_.begin(), rhs.coef_times_gamma_.end());
    if (tmp4 != 0) return (tmp4 < 0);

    return gamma_ < rhs.gamma_;
  }


  //==============================================================================
  bool LRBSpline3D::operator==(const LRBSpline3D &rhs) const
  //==============================================================================
  {
  const bool tmp1 = ((*bspline_u_) == (*rhs.bspline_u_));
  if (tmp1 == false)
    return false;

  const bool tmp2 = ((*bspline_v_) == (*rhs.bspline_v_));
  if (tmp2 == false)
    return false;

  const bool tmp3 = ((*bspline_w_) == (*rhs.bspline_w_));
  if (tmp3 == false)
    return false;

  return true;
  }


//==============================================================================
bool LRBSpline3D::checkOverload()
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
  void LRBSpline3D::reverseParameterDirection(int pardir)
  //==============================================================================
  {
  if (pardir == 1)
    bspline_u_->reverseParameterDirection();
  else if (pardir == 2)
    bspline_v_->reverseParameterDirection();
  else
    bspline_w_->reverseParameterDirection();
  }


  //==============================================================================
  void LRBSpline3D::swapParameterDirection(int pardir1, int pardir2)
  //==============================================================================
  {
    MESSAGE("(): Not implemented.");
    throw;
  }

}; // end namespace Go
