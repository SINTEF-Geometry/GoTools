//===========================================================================
//                                                                           
// File: LRBSpline3D.C                                                       
//                                                                           
// Created: Mon Mar 11 16:17:14 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
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

  coef_fixed_ = 0;
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
  
  coef_fixed_ = 0;
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

  // //==============================================================================
  // double LRBSpline3D::evalBasisFunction(double u, double v, double w,
  //                                       const double* const kvals_u,
  //                                       const double* const kvals_v,
  //                                       const double* const kvals_w,
  // 					int u_deriv, int v_deriv, int w_deriv,
  // 					bool u_at_end, bool v_at_end, bool w_at_end) const
  
  // //==============================================================================
  // {
  //   double bval1 = bspline_u_->evalBasisFunction(u, u_deriv, u_at_end);
  //   double bval2 = bspline_v_->evalBasisFunction(v, v_deriv, v_at_end);
  //   double bval3 = bspline_w_->evalBasisFunction(w, w_deriv, w_at_end);
  //   return bval1*bval2*bval3;
  // }

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
       THROW("evalder_add for volumes and rationals not implemented");
       // int dim = coef_times_gamma_.dimension();
       // int nmb = (deriv+1)*(deriv+2)/2;
       // double tmp[(int)((MAX_DER+1)*(MAX_DER+2)*(MAX_DIM+1))];
       // double *tmpder = tmp;
       // double val;
       // Point tmppt(dim);
       // kh = 0;
       // for (ki=0; ki<=deriv; ++ki)
       // 	 for (kj=0; kj<=ki; ++kj, ++kh)
       // 	   {
       // 	     val = weight_*bder1[ki-kj]*bder2[kj];
       // 	     for (kr=0; kr<dim; ++kr)
       // 	       tmpder[kh*(dim+1)+kr] = coef_times_gamma_[kr]*val;
       // 	     tmpder[kh*(dim+1)+dim] = val;
       // 	   }
       // double *tmpder2 = tmpder+nmb*(dim+1);
       // SplineUtils::surface_ratder(tmpder, dim, deriv, tmpder2);
       // for (kh=0; kh<nmb; ++kh)
       // 	 {
       // 	   for (kr=0; kr<dim; ++kr)
       // 	     tmppt[kr] = tmpder2[kh*dim+kr];
       // 	   der[kh] = tmppt;
       // 	 }
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
