//===========================================================================
//                                                                           
// File: LRSplineVolume.C                                                    
//                                                                           
// Created: Tue Feb 26 11:28:58 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/utils/checks.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRSpline3DUtils.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include <algorithm>
#include <iostream>
#include <fstream>

using std::vector;
using std::unique_ptr;
using std::get;

//==============================================================================
namespace Go
//==============================================================================
{



// =============================================================================
LRSplineVolume::ElemKey 
  LRSplineVolume::generate_key(const double& umin, 
			       const double& vmin,
			       const double& wmin)
// =============================================================================
{
    ElemKey key = { umin, vmin, wmin};
  return key;
}




//==============================================================================
LRSplineVolume::LRSplineVolume(SplineVolume *vol, double knot_tol)
//==============================================================================
: knot_tol_(knot_tol), rational_(vol->rational()), curr_element_(NULL),
  mesh_(vol->basis(0).begin(), vol->basis(0).end(),
	vol->basis(1).begin(), vol->basis(1).end(),
	vol->basis(2).begin(), vol->basis(2).end())
{
  std::vector<int> knot_ixs_u = init_knot_indices(mesh_, XDIR);
  std::vector<int> knot_ixs_v = init_knot_indices(mesh_, YDIR);
  std::vector<int> knot_ixs_w = init_knot_indices(mesh_, ZDIR);

  int deg_u = vol->order(0) - 1;
  int deg_v = vol->order(1) - 1;
  int deg_w = vol->order(2) - 1;
  int coefs_u = vol->numCoefs(0);
  int coefs_v = vol->numCoefs(1);
  int coefs_w = vol->numCoefs(2);
  std::vector<double>::iterator rcoefs = vol->rcoefs_begin();
  std::vector<double>::iterator coefs = vol->coefs_begin();
  int dim = vol->dimension();
  int kdim = (rational_) ? dim + 1 : 0; // Setting the kdim to 0 for non-rational case as the rcoef vector will then be empty.

  // Store uni-variate B-splines
  bsplinesuni1_.resize(coefs_u);
  bsplinesuni2_.resize(coefs_v);
  bsplinesuni3_.resize(coefs_w);
  for (int u_ix = 0; u_ix != coefs_u; ++u_ix) {
    bsplinesuni1_[u_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(1, deg_u,
							       knot_ixs_u.begin() + u_ix,
							       &mesh_)));
  }

  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    bsplinesuni2_[v_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(2, deg_v,
							       knot_ixs_v.begin() + v_ix,
							       &mesh_)));
  }

  for (int w_ix = 0; w_ix != coefs_w; ++w_ix)  {
    bsplinesuni3_[w_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(3, deg_w,
							       knot_ixs_w.begin() + w_ix,
							       &mesh_)));
  }

  for (int w_ix = 0; w_ix != coefs_w; ++w_ix)  {
    for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
      for (int u_ix = 0; u_ix != coefs_u; ++u_ix, coefs+=dim, rcoefs += kdim)
	{
	  double wgt = (rational_) ? rcoefs[dim] : 1.0;
	  unique_ptr<LRBSpline3D> b(new LRBSpline3D(Point(coefs, coefs + dim),
						    wgt,
						    bsplinesuni1_[u_ix].get(),
						    bsplinesuni2_[v_ix].get(),
						    bsplinesuni3_[w_ix].get(),
						    1.0, rational_));
	  LRSplineVolume::BSKey bs_key = generate_key(*b, mesh_);
	  bsplines_.insert(std::pair<LRSplineVolume::BSKey, unique_ptr<LRBSpline3D> >(bs_key,std::move(b)));
	}
    }
  }

  emap_ = construct_element_map_(mesh_, bsplines_);
}

  // Copy constructor
//==============================================================================
  LRSplineVolume::LRSplineVolume(const LRSplineVolume& rhs)
//==============================================================================
   : knot_tol_(rhs.knot_tol_), rational_(rhs.rational_),  curr_element_(NULL),
     mesh_(rhs.mesh_)
  {
    // Clone univariate B-splines
    vector<std::unique_ptr<BSplineUniLR> >::const_iterator curruni1 = 
      rhs.bsplinesuni1_.begin();
    vector<std::unique_ptr<BSplineUniLR> >::const_iterator enduni1 = 
      rhs.bsplinesuni1_.end();
    bsplinesuni1_.resize(rhs.bsplinesuni1_.size());
    size_t ki;
    for (ki=0; curruni1 != enduni1; ++ki, ++curruni1)
      {
	unique_ptr<BSplineUniLR> b(new BSplineUniLR(*(*curruni1)));

	// Update mesh pointer
	b->setMesh(&mesh_);
	bsplinesuni1_[ki] = std::move(b);
      }
  
    vector<std::unique_ptr<BSplineUniLR> >::const_iterator curruni2 = 
      rhs.bsplinesuni2_.begin();
    vector<std::unique_ptr<BSplineUniLR> >::const_iterator enduni2 = 
      rhs.bsplinesuni2_.end();
    bsplinesuni2_.resize(rhs.bsplinesuni2_.size());
    for (ki=0; curruni2 != enduni2; ++ki, ++curruni2)
      {
	unique_ptr<BSplineUniLR> b(new BSplineUniLR(*(*curruni2)));

	// Update mesh pointer
	b->setMesh(&mesh_);
	bsplinesuni2_[ki] = std::move(b);
      }

    vector<std::unique_ptr<BSplineUniLR> >::const_iterator curruni3 = 
      rhs.bsplinesuni3_.begin();
    vector<std::unique_ptr<BSplineUniLR> >::const_iterator enduni3 = 
      rhs.bsplinesuni3_.end();
    bsplinesuni3_.resize(rhs.bsplinesuni3_.size());
    for (ki=0; curruni3 != enduni3; ++ki, ++curruni3)
      {
	unique_ptr<BSplineUniLR> b(new BSplineUniLR(*(*curruni3)));

	// Update mesh pointer
	b->setMesh(&mesh_);
	bsplinesuni3_[ki] = std::move(b);
      }

    // Clone LR B-splines
    int left1 = 0, left2 = 0, left3 = 0;
    BSplineMap::const_iterator curr = rhs.basisFunctionsBegin();
    BSplineMap::const_iterator end = rhs.basisFunctionsEnd();
    for (; curr != end; ++curr)
      {
        unique_ptr<LRBSpline3D> b(new LRBSpline3D(*curr->second));
 
	const BSplineUniLR *uni1 = b->getUnivariate(XDIR);
	bool found1 = BSplineUniUtils::identify_bsplineuni(uni1, bsplinesuni1_, left1);
	if (!found1)
	  THROW("Univariate B-spline not found");
	b->setUnivariate(XDIR, bsplinesuni1_[left1].get());

	const BSplineUniLR *uni2 = b->getUnivariate(YDIR);
	bool found2 = BSplineUniUtils::identify_bsplineuni(uni2, bsplinesuni2_, left2);
	if (!found2)
	  THROW("Univariate B-spline not found");
	b->setUnivariate(YDIR, bsplinesuni2_[left2].get());

	const BSplineUniLR *uni3 = b->getUnivariate(ZDIR);
	bool found3 = BSplineUniUtils::identify_bsplineuni(uni3, bsplinesuni3_, left3);
	if (!found3)
	  THROW("Univariate B-spline not found");
	b->setUnivariate(ZDIR, bsplinesuni3_[left3].get());

       LRSplineVolume::BSKey bs_key = generate_key(*b, mesh_);
        bsplines_.insert(std::pair<LRSplineVolume::BSKey, unique_ptr<LRBSpline3D> >(bs_key, std::move(b)));
      }

    // The ElementMap has to be generated and cannot be copied directly, since it
    // contains raw pointers.
    emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
   void  LRSplineVolume::read(std::istream& is)
//==============================================================================
{
  int rat = -1;
  object_from_stream(is, rat);
  rational_ = (rat == 1);

  // reading knot tolerances and the mesh
  object_from_stream(is, knot_tol_);
  object_from_stream(is, mesh_);

  // Reading all basis functions
  int num_bfuns;
  object_from_stream(is, num_bfuns);

  vector<unique_ptr<LRBSpline3D> > dummy_vec(num_bfuns);

  int left1 = 0, left2 = 0, left3 = 0;
  for (int i = 0; i != num_bfuns; ++i) {
    unique_ptr<LRBSpline3D> b(new LRBSpline3D());
//    LRBSpline3D* b = new LRBSpline3D();
    //object_from_stream(is, *b);
    b->read(is, bsplinesuni1_, left1, bsplinesuni2_, left2, 
	    bsplinesuni3_, left3);
    // We set the global mesh in the b basis function.
    b->setMesh(&mesh_);
    BSKey key = generate_key(*b, mesh_);
    bsplines_.insert(std::make_pair(key, std::move(b)));
  }

  // Reconstructing element map
  emap_ = construct_element_map_(mesh_, bsplines_);

  rational_ = rational_;

  auto it = bsplines_.begin();
  while (it != bsplines_.end())
    {
      it->second->setMesh(&mesh_);
      ++it;
    }

  curr_element_ = NULL;

#ifndef NDEBUG
  {
    vector<LRBSpline3D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    //puts("Remove when done debugging!");
    vector<Element3D*> elems;
    for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
    {
	elems.push_back(((*iter).second.get()));
    }
    //puts("Remove when done debugging!");
    int stop_break = 1;
  }
#endif
}

//==============================================================================
   void LRSplineVolume::write(std::ostream& os) const
 //==============================================================================
{
  std::streamsize prev = os.precision(15);

  int rat = (rational_) ? 1 : 0;
  object_to_stream(os, rat);
  object_to_stream(os, knot_tol_);
  object_to_stream(os, '\n');
  object_to_stream(os, mesh_);

  object_to_stream(os, bsplines_.size());
  object_to_stream(os, '\n');
  for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) 
    {
      object_to_stream(os, *(b->second));
      object_to_stream(os, '\n');
    }

  // NB: 'emap_' is not saved (contains raw pointers to other data).  
  // Instead, it will be regenerated the LRSplineSurface is read().
    os.precision(prev);   // Reset precision to it's previous value
}

//==============================================================================
   BoundingBox LRSplineVolume::boundingBox() const
//==============================================================================
{
    BSplineMap::const_iterator curr = basisFunctionsBegin();
    BoundingBox box(curr->second->Coef(), curr->second->Coef());
    BSplineMap::const_iterator end = basisFunctionsEnd();
    for (; curr!=end; ++curr)
      {
	const Point coef = curr->second->Coef();
	box.addUnionWith(coef);
      }
    return box;
}

//===========================================================================
   const Array<double, 6> LRSplineVolume::parameterSpan() const
//===========================================================================
{
  domain_[0] = mesh_.minParam(XDIR);
  domain_[1] = mesh_.maxParam(XDIR);
  domain_[2] = mesh_.minParam(YDIR);
  domain_[3] = mesh_.maxParam(YDIR);
  domain_[4] = mesh_.minParam(ZDIR);
  domain_[5] = mesh_.maxParam(ZDIR);
  return domain_;
}

//==============================================================================
const LRSplineVolume& LRSplineVolume::operator= (const LRSplineVolume& other)
//==============================================================================
{
  LRSplineVolume lr_spline_vol(other);
  this->swap(lr_spline_vol);
  return *this;
}

//==============================================================================
  void LRSplineVolume::swap(LRSplineVolume& rhs)
//==============================================================================
{
  std::swap(knot_tol_,    rhs.knot_tol_);
  std::swap(rational_,    rhs.rational_);
  std::swap(mesh_    ,    rhs.mesh_);
  std::swap(bsplinesuni1_,    rhs.bsplinesuni1_);
  std::swap(bsplinesuni2_,    rhs.bsplinesuni2_);
  std::swap(bsplinesuni3_,    rhs.bsplinesuni3_);
  std::swap(bsplines_,    rhs.bsplines_);
  std::swap(emap_    ,    rhs.emap_);

  // Must update mesh pointer in B-splines
  for (auto b_it = bsplines_.begin(); b_it != bsplines_.end(); ++b_it) 
    {
      b_it->second->setMesh(&mesh_);
    }
  // Must update mesh pointer in B-splines
  for (auto b_it = rhs.bsplines_.begin(); b_it != rhs.bsplines_.end(); ++b_it) 
    {
      b_it->second->setMesh(&rhs.mesh_);
    }
}

//==============================================================================
   ClassType LRSplineVolume::instanceType() const
//==============================================================================
{
  return classType();
}
    
//==============================================================================
   void LRSplineVolume::point(Point& pt, double upar, double vpar, double wpar) const
//==============================================================================
{
  if (rational_)
    pt = operator()(upar, vpar, wpar, 0, 0, 0);
  else
    {
      Element3D* elem;
      if (curr_element_ && curr_element_->contains(upar, vpar, wpar))
  	elem = curr_element_;
      else
  	{
  	  elem = coveringElement(upar, vpar, wpar);
  	  curr_element_ = (Element3D*)elem;
  	}
      
      point(pt, upar, vpar, wpar, elem);
    }
}

//===========================================================================
   void LRSplineVolume::point(Point& pt, double upar, double vpar, double wpar,
                              Element3D* elem) const
//===========================================================================
{
  if (rational_)
    pt = operator()(upar, vpar, wpar, 0, 0, 0, elem);
  else
    {
      if (elem && elem->contains(upar, vpar, wpar))
	curr_element_ = elem;
      else if (curr_element_ && curr_element_->contains(upar, vpar, wpar))
	elem = curr_element_;
      else
	{
	  elem = coveringElement(upar, vpar, wpar);
	  curr_element_ = (Element3D*)elem;
	}
      
      double eps = 1.0e-12;
      const bool u_at_end = (upar >= mesh_.maxParam(XDIR)-eps); 
      const bool v_at_end = (vpar >= mesh_.maxParam(YDIR)-eps);
      const bool w_at_end = (wpar >= mesh_.maxParam(ZDIR)-eps); 
      const vector<LRBSpline3D*>& bfunctions = curr_element_->getSupport();
      size_t bsize = bfunctions.size();
      vector<double> val(3*bsize);
      pt.resize(this->dimension());
      pt.setValue(0.0);
      for (size_t ki=0; ki<bsize; ++ki)
	  {
	    size_t kj;
	    const BSplineUniLR* uni1 =  bfunctions[ki]->getUnivariate(XDIR);
	    const BSplineUniLR* uni2 =  bfunctions[ki]->getUnivariate(YDIR);
	    const BSplineUniLR* uni3 =  bfunctions[ki]->getUnivariate(ZDIR);
	    for (kj=0; kj<ki; ++kj)
	      if (uni1 == bfunctions[kj]->getUnivariate(XDIR))
		break;
	    if (kj < ki)
	      val[ki] = val[kj];
	    else
	      val[ki] = uni1->evalBasisFunction(upar, 0, u_at_end);

	    for (kj=0; kj<ki; ++kj)
	      if (uni2 == bfunctions[kj]->getUnivariate(YDIR))
		break;
	    if (kj < ki)
	      val[bsize+ki] = val[bsize+kj];
	    else
	      val[bsize+ki] = 
		uni2->evalBasisFunction(vpar, 0, v_at_end);

	    for (kj=0; kj<ki; ++kj)
	      if (uni3 == bfunctions[kj]->getUnivariate(ZDIR))
		break;
	    if (kj < ki)
	      val[2*bsize+ki] = val[2*bsize+kj];
	    else
	      val[2*bsize+ki] = 
		uni3->evalBasisFunction(wpar, 0, w_at_end);

	    pt += val[ki]*val[bsize+ki]*val[2*bsize+ki]*
	      bfunctions[ki]->coefTimesGamma();
	  }
    }
 }

//==============================================================================
   void LRSplineVolume::point(std::vector<Point>& pts, 
		     double upar, double vpar, double wpar,
		     int derivs,
		     bool u_from_right,
		     bool v_from_right,
		     bool w_from_right,
		     double resolution) const
//==============================================================================
{
  point(pts, upar, vpar, wpar, derivs, curr_element_, u_from_right,
	v_from_right, w_from_right, resolution);
}

   //===========================================================================
  void LRSplineVolume::point(vector<Point>& pts, 
			     double upar, double vpar, double wpar,
			     int derivs,
			     Element3D* elem,
			     bool u_from_right,
			     bool v_from_right,
			     bool w_from_right,
			     double resolution) const
  //===========================================================================
  {
    int totpts = (derivs + 1)*(derivs + 2)*(derivs + 3)/6;
    DEBUG_ERROR_IF((int)pts.size() < totpts, "The vector of points must have sufficient size.");

    int dim = dimension();
    for (int ki = 0; ki < totpts; ++ki) {
	if (pts[ki].dimension() != dim) {
	    pts[ki].resize(dim);
	}
	pts[ki].setValue(0.0);
    }

  if (upar < paramMin(XDIR)) {
    //MESSAGE("upar was outside domain: " << upar << " < " << paramMin(XDIR) << ", moved inside.");
      upar = paramMin(XDIR);
  } else if (upar > paramMax(XDIR)) {
    //MESSAGE("upar was outside domain: " << upar << " > " << paramMax(XDIR) << ", moved inside.");
      upar = paramMax(XDIR);
  }

  if (vpar < paramMin(YDIR)) {
    //MESSAGE("vpar was outside domain: " << vpar << " < " << paramMin(YDIR) << ", moved inside.");
      vpar = paramMin(YDIR);
  } else if (vpar > paramMax(YDIR)) {
    //MESSAGE("vpar was outside domain: " << vpar << " > " << paramMax(YDIR) << ", moved inside.");
      vpar = paramMax(YDIR);
  }
    
  if (wpar < paramMin(ZDIR)) {
    //MESSAGE("wpar was outside domain: " << wpar << " < " << paramMin(ZDIR) << ", moved inside.");
      wpar = paramMin(ZDIR);
  } else if (wpar > paramMax(ZDIR)) {
    //MESSAGE("wpar was outside domain: " << wpar << " > " << paramMax(ZDIR) << ", moved inside.");
      wpar = paramMax(ZDIR);
  }
    
  // Check element
  if (!elem || !elem->contains(upar, vpar, wpar))
    {
      elem = coveringElement(upar, vpar, wpar);
    }
    curr_element_ = elem;

  const vector<LRBSpline3D*>& covering_B_functions = elem->getSupport();

  for (size_t kr=0; kr<covering_B_functions.size(); ++kr)
    {
      covering_B_functions[kr]->evalder_add(upar, vpar, wpar, derivs, &pts[0], 
					    u_from_right, v_from_right,
					    w_from_right);
    }
  }


//==============================================================================
void LRSplineVolume::elementGridEvaluate(Element3D *element,
					 vector<double>& upar,
					 vector<double>& vpar,
					 vector<double>& wpar,
					 vector<double>& points) const
//==============================================================================
{
  // Result of evaluation is stored sequentially from first to last parameter
  // value. The u-parameter runs fastest, then comes v and finally w
  // The function will throw if any parameter value is outside the element
  int dim = dimension();
  points.resize(upar.size()*vpar.size()*wpar.size()*dim);
  std::fill(points.begin(), points.end(), 0.0);

  double umin = element->umin();
  double umax = element->umax();
  double vmin = element->vmin();
  double vmax = element->vmax();
  double wmin = element->wmin();
  double wmax = element->wmax();

  const std::vector<LRBSpline3D*> bfunctions = element->getSupport();
  size_t bsize = bfunctions.size();

  // Evaluate univariate B-splines
  std::vector<double> val(bsize*(upar.size()+vpar.size()+wpar.size()));
  double *val1 = &val[0];
  double *val2 = val1 + bsize*upar.size();
  double *val3 = val2 + bsize*vpar.size();
  for (size_t ki=0; ki<bsize; ++ki)
    {
      size_t kj;
      const BSplineUniLR* uni1 =  bfunctions[ki]->getUnivariate(XDIR);
      const BSplineUniLR* uni2 =  bfunctions[ki]->getUnivariate(YDIR);
      const BSplineUniLR* uni3 =  bfunctions[ki]->getUnivariate(ZDIR);
      for (kj=0; kj<ki; ++kj)
	if (uni1 == bfunctions[kj]->getUnivariate(XDIR))
	  break;
      if (kj < ki)
	std::copy(val1+kj*upar.size(), val1+(kj+1)*upar.size(), val1+ki*upar.size());
      else
	{
	  for (size_t ka=0; ka<upar.size(); ++ka)
	    { 
	      if (upar[ka] < umin || upar[ka] > umax)
		THROW("u-parameter out of range");
	      const bool on_end = (upar[ka] == bfunctions[ki]->umax());
	      val1[ki*upar.size()+ka] = uni1->evalBasisFunction(upar[ka], 0, on_end);
	    }
	}
		
      for (kj=0; kj<ki; ++kj)
	if (uni2 == bfunctions[kj]->getUnivariate(YDIR))
	  break;
      if (kj < ki)
	std::copy(val2+kj*vpar.size(), val2+(kj+1)*vpar.size(), val2+ki*vpar.size());
      else
	{
	  for (size_t ka=0; ka<vpar.size(); ++ka)
	    {
	      if (vpar[ka] < vmin || vpar[ka] > vmax)
		THROW("v-parameter out of range");
	      const bool on_end = (vpar[ka] == bfunctions[ki]->vmax());
	      val2[ki*vpar.size()+ka] = uni2->evalBasisFunction(vpar[ka], 0, on_end);
	    }
	}
		
      for (kj=0; kj<ki; ++kj)
	if (uni3 == bfunctions[kj]->getUnivariate(ZDIR))
	  break;
      if (kj < ki)
	      std::copy(val3+kj*wpar.size(), val3+(kj+1)*wpar.size(), val3+ki*wpar.size());
      else
	{
	  for (size_t ka=0; ka<wpar.size(); ++ka)
	    {
	      if (wpar[ka] < wmin || wpar[ka] > wmax)
		THROW("w-parameter out of range");
	      const bool on_end = (wpar[ka] == bfunctions[ki]->wmax());
	      val3[ki*wpar.size()+ka] = uni3->evalBasisFunction(wpar[ka], 0, on_end);
	    }
	}
    }

  Point pt;
  double *p=&points[0];
  int ix2=0;
  for (size_t i=0; i<wpar.size(); ++i)
    for (size_t j=0; j<vpar.size(); ++j)
      for (size_t k=0; k<upar.size(); ++k, p+=dim)
	for (size_t ix=0; ix<bsize; ++ix)
	  {
	    pt = bfunctions[ix]->coefTimesGamma();
	    for (int ka=0; ka<dim; ++ka)
	      p[ka] +=
		pt[ka]*val1[ix*upar.size()+k]*
		val2[ix*vpar.size()+j]*
		val3[ix*wpar.size()+i];
	  }
		  
  
}

//==============================================================================
void LRSplineVolume::elementGridEvaluate(Element3D *element,
					 double *upar, int usize,
					 double *vpar, int vsize,
					 double *wpar, int wsize,
					 int deriv,
					 vector<double>& points) const
//==============================================================================
{
  // Result of evaluation is stored sequentially from first to last parameter
  // value. The u-parameter runs fastest, then comes v and finally w
  /// Derivatives are placed directly after the associated point. The sequence
  /// is F, Fu, Fv, Fw, Fuu, Fuv, Fuw, Fvv, Fvw, Fww, Fuuu, ...
  // The function will throw if any parameter value is outside the element
  int dim = dimension();
  int totder = (deriv + 1)*(deriv + 2)*(deriv + 3)/6;
  points.resize(usize*vsize*wsize*totder*dim);
  std::fill(points.begin(), points.end(), 0.0);

  double umin = element->umin();
  double umax = element->umax();
  double vmin = element->vmin();
  double vmax = element->vmax();
  double wmin = element->wmin();
  double wmax = element->wmax();

  const std::vector<LRBSpline3D*> bfunctions = element->getSupport();
  size_t bsize = bfunctions.size();

  // Evaluate univariate B-splines
  int usize2 = usize*(deriv+1);
  int vsize2 = vsize*(deriv+1);
  int wsize2 = wsize*(deriv+1);

  std::vector<double> val(bsize*(usize2+vsize2+wsize2), 0.0);
  double *val1 = &val[0];
  double *val2 = val1 + bsize*usize2;
  double *val3 = val2 + bsize*vsize2;
  for (size_t ki=0; ki<bsize; ++ki)
    {
      size_t kj;
      const BSplineUniLR* uni1 =  bfunctions[ki]->getUnivariate(XDIR);
      const BSplineUniLR* uni2 =  bfunctions[ki]->getUnivariate(YDIR);
      const BSplineUniLR* uni3 =  bfunctions[ki]->getUnivariate(ZDIR);
      for (kj=0; kj<ki; ++kj)
	if (uni1 == bfunctions[kj]->getUnivariate(XDIR))
	  break;
      if (kj < ki)
	std::copy(val1+kj*usize2, val1+(kj+1)*usize2, val1+ki*usize2);
      else
	{
	  for (int ka=0; ka<usize; ++ka)
	    { 
	      if (upar[ka] < umin || upar[ka] > umax)
		{
		  printf("u-parameter out of range \n");
		  printf("%d %d (%f, %f) - %f \n",ka, usize, umin, umax, upar[ka]);   
		THROW("u-parameter out of range");
		}
	      const bool on_end = (upar[ka] == bfunctions[ki]->umax());
	      uni1->evalBasisFunctions(upar[ka], deriv,
				       &val1[ki*usize2+ka*(deriv+1)], on_end);
	    }
	}
		
      for (kj=0; kj<ki; ++kj)
	if (uni2 == bfunctions[kj]->getUnivariate(YDIR))
	  break;
      if (kj < ki)
	std::copy(val2+kj*vsize2, val2+(kj+1)*vsize2, val2+ki*vsize2);
      else
	{
	  for (int ka=0; ka<vsize; ++ka)
	    {
	      if (vpar[ka] < vmin || vpar[ka] > vmax)
		{
		  printf("v-parameter out of range \n");
		THROW("v-parameter out of range");
		}
	      const bool on_end = (vpar[ka] == bfunctions[ki]->vmax());
	      uni2->evalBasisFunctions(vpar[ka], deriv,
				       &val2[ki*vsize2+ka*(deriv+1)], on_end);
	    }
	}
		
      for (kj=0; kj<ki; ++kj)
	if (uni3 == bfunctions[kj]->getUnivariate(ZDIR))
	  break;
      if (kj < ki)
	std::copy(val3+kj*wsize2, val3+(kj+1)*wsize2, val3+ki*wsize2);
      else
	{
	  for (int ka=0; ka<wsize; ++ka)
	    {
	      if (wpar[ka] < wmin || wpar[ka] > wmax)
		{
		  printf("w-parameter out of range \n");
		THROW("w-parameter out of range");
		}
	      const bool on_end = (wpar[ka] == bfunctions[ki]->wmax());
	      uni3->evalBasisFunctions(wpar[ka], deriv,
				       &val3[ki*wsize2+ka*(deriv+1)], on_end);
	    }
	}
    }

  Point pt;
  double *p=&points[0];
  int ix2=0;
  double der;
  int d1, d2, d3;
  int id, jd;
  int kh;
  size_t ux, vx, wx;
  for (int i=0; i<wsize; ++i)
    {
      id = i*(deriv+1);
      for (int j=0; j<vsize; ++j)
	{
	    jd = j*(deriv+1);
	    for (int k=0; k<usize; ++k, p+=(dim*totder))
	      for (size_t ix=0; ix<bsize; ++ix)
		{
		  pt = bfunctions[ix]->coefTimesGamma();
		  kh = 0;
		  wx = ix*wsize2 + id;
		  vx = ix*vsize2 + jd;
		  ux = ix*usize2 + k*(deriv+1);
		  for (d1=0; d1<=deriv; ++d1)
		    for (d2=0; d2<=d1; ++d2)
		      for (d3=0; d3<=d2; ++d3, ++kh)
			{
			  der = val1[ux+d1-d2]*val2[vx+d2-d3]*val3[wx+d3];
			  
			  // double bder1 =
			  //   bfunctions[ix]->evalBasisFunction(upar[k], vpar[j], wpar[i],
			  // 					d1-d2, d2-d3, d3);
			  for (int ka=0; ka<dim; ++ka)
			    p[kh*dim+ka] += pt[ka]*der;
			}
		}
	}
    }
		  
  
}

//==============================================================================
   double LRSplineVolume::startparam_u() const
//==============================================================================
{
    return paramMin(XDIR);
}

//==============================================================================
   double LRSplineVolume::startparam_v() const
//==============================================================================
{
    return paramMin(YDIR);
}

//==============================================================================
   double LRSplineVolume::startparam_w() const
//==============================================================================
{
    return paramMin(ZDIR);
}

//==============================================================================
   double LRSplineVolume::endparam_u() const
//==============================================================================
{
    return paramMax(XDIR);
}

//==============================================================================
   double LRSplineVolume::endparam_v() const
//==============================================================================
{
    return paramMax(YDIR);
}

//==============================================================================
   double LRSplineVolume::endparam_w() const
//==============================================================================
{
    return paramMax(ZDIR);
}


//==============================================================================
DirectionCone LRSplineVolume::tangentCone(int pardir) const
//==============================================================================
{
   MESSAGE("LRSplineVolume::tangentCone(): Not implemented yet.");
   throw;
}

//==============================================================================
  LRSplineVolume* LRSplineVolume::derivVolume(int ider1, int ider2, int ider3) const
//==============================================================================
{
   MESSAGE("LRSplineVolume:: Not implemented yet.");
   throw;
}

//==============================================================================
  LRSplineVolume*
  LRSplineVolume::subVolume(double from_upar, double from_vpar, double from_wpar,
			    double to_upar, double to_vpar, double to_wpar,
			    double fuzzy) const
//==============================================================================
{
    // Check input
    if (from_upar >= to_upar) {
	THROW("First u-parameter must be strictly less than second.");
    }
    if (from_vpar >= to_vpar) {
	THROW("First v-parameter must be strictly less than second.");
    }
    if (from_wpar >= to_wpar) {
	THROW("First w-parameter must be strictly less than second.");
    }
    if (from_upar < startparam_u()-fuzzy || from_vpar < startparam_v()-fuzzy ||
	from_wpar < startparam_w()-fuzzy) {
	THROW("Subvolume defined outside volume.");
    }

    // Periodic volumes not considered, i.e. it is not possible to pick a
    // part of a surface over a periodic seam.

    // If boundaries are close to existing knots, we snap.
    // Note that the indices refer to the vectors of unique knots.
    int ix1 = mesh_.knotIntervalFuzzy(XDIR, from_upar, fuzzy);
    int ix2 = mesh_.knotIntervalFuzzy(XDIR, to_upar, fuzzy);
    int iy1 = mesh_.knotIntervalFuzzy(YDIR, from_vpar, fuzzy);
    int iy2 = mesh_.knotIntervalFuzzy(YDIR, to_vpar, fuzzy);
    int iz1 = mesh_.knotIntervalFuzzy(ZDIR, from_wpar, fuzzy);
    int iz2 = mesh_.knotIntervalFuzzy(ZDIR, to_wpar, fuzzy);

    LRSplineVolume *vol = NULL;

    int deg1 = degree(XDIR);
    int deg2 = degree(YDIR);
    int deg3 = degree(ZDIR);
    int nmb1 = mesh_.numDistinctKnots(XDIR);
    int nmb2 = mesh_.numDistinctKnots(YDIR);
    int nmb3 = mesh_.numDistinctKnots(ZDIR);

    // Check if the sub surface is identical to the current surface
    bool identical = true;
    if (ix1 != 0 || mesh_.kval(XDIR, ix1) != from_upar ||
	mesh_.minMultInLine(XDIR, ix1) < deg1+1)
      identical = false;
    if (ix2 != mesh_.lastMeshVecIx(XDIR) || mesh_.kval(XDIR, ix2) != to_upar ||
	mesh_.minMultInLine(XDIR, ix2) < deg1+1)
      identical = false;
    if (iy1 != 0 || mesh_.kval(YDIR, iy1) != from_vpar ||
	mesh_.minMultInLine(YDIR, iy1) < deg2+1)
      identical = false;
     if (iy2 != mesh_.lastMeshVecIx(YDIR) || mesh_.kval(YDIR, iy2) != to_vpar ||
	mesh_.minMultInLine(YDIR, iy2) < deg2+1)
      identical = false;
    if (iz1 != 0 || mesh_.kval(ZDIR, iz1) != from_wpar ||
	mesh_.minMultInLine(ZDIR, iz1) < deg3+1)
      identical = false;
     if (iz2 != mesh_.lastMeshVecIx(ZDIR) || mesh_.kval(ZDIR, iz2) != to_wpar ||
	mesh_.minMultInLine(ZDIR, iz2) < deg3+1)
      identical = false;
    
     if (identical)
       {
	 // The sub volume is equal to the current. Copy.
	 vol = new LRSplineVolume(*this);
	 return vol;
       }

     // Make a copy of the current volume
     shared_ptr<LRSplineVolume> vol2(new LRSplineVolume(*this));

     // We locate all basis functions for which the support is crossed
     // by a subsurface boundary line.
     // We store the min & max index in each dir of those functions.
     // @@sbr201301 For large cases It may be faster to go through the elements.
     int umin_ind = nmb1; // Higher than max index.
     int umax_ind = 0;
     int vmin_ind = nmb2; // Higher than max index.
     int vmax_ind = 0;
     int wmin_ind = nmb3; // Higher than max index.
     int wmax_ind = 0;
     for (auto iter = basisFunctionsBegin(); iter != basisFunctionsEnd(); ++iter)
       {
	 LRBSpline3D* bas_func = iter->second.get();

	 if ((from_upar >= bas_func->umin()) && (from_upar <= bas_func->umax()))
	   umin_ind = std::min(umin_ind, bas_func->suppMin(XDIR));
	 if (to_upar >= bas_func->umin() && to_upar <= bas_func->umax())
	   umax_ind = std::max(umax_ind, bas_func->suppMax(XDIR));

	 if (from_vpar >= bas_func->vmin() && from_vpar <= bas_func->vmax())
	   vmin_ind = std::min(vmin_ind, bas_func->suppMin(YDIR));
	 if (to_vpar >= bas_func->vmin() && to_vpar <= bas_func->vmax())
	   vmax_ind = std::max(vmax_ind, bas_func->suppMax(YDIR));

	 if (from_wpar >= bas_func->wmin() && from_wpar <= bas_func->wmax())
	   wmin_ind = std::min(wmin_ind, bas_func->suppMin(ZDIR));
	 if (to_wpar >= bas_func->wmin() && to_wpar <= bas_func->wmax())
	   wmax_ind = std::max(wmax_ind, bas_func->suppMax(ZDIR));
       }

     // Define possible new knotlines
     // Note that the new knot lines must be longer than the size of the
     // sub surface to avoid LR B-splines partly overlapping the sub domain
     // @@@ VSK. Could the extension be smaller than prescribed here?
     vector<Refinement3D> refs; ///(6);
     //vector<Refinement2D> refs;
#if 0 // Old version, we need to consider the support of all the basis
      // functions.
     double umin = mesh_.kval(XDIR, std::max(ix1 - deg1, 0));
     double umax = mesh_.kval(XDIR, std::min(ix2 + deg1, nmb1-1));
     double vmin = mesh_.kval(YDIR, std::max(iy1 - deg2, 0));
     double vmax = mesh_.kval(YDIR, std::min(iy2 + deg2, nmb2-1));
#else
     double umin = mesh_.kval(XDIR, umin_ind);
     double umax = mesh_.kval(XDIR, umax_ind);
     double vmin = mesh_.kval(YDIR, vmin_ind);
     double vmax = mesh_.kval(YDIR, vmax_ind);
     double wmin = mesh_.kval(ZDIR, wmin_ind);
     double wmax = mesh_.kval(ZDIR, wmax_ind);
#endif
     if (from_upar-umin>fuzzy)
       refs.push_back(Refinement3D(from_upar, vmin, vmax, wmin, wmax,
     				   XDIR, deg1+1));
     if (umax-to_upar>fuzzy)
       refs.push_back(Refinement3D(to_upar, vmin, vmax, wmin, wmax,
     				   XDIR, deg1+1));
     if (from_vpar-vmin>fuzzy)
       refs.push_back(Refinement3D(from_vpar, wmin, wmax, umin, umax,
     				   YDIR, deg2+1));
     if (vmax-to_vpar>fuzzy)
       refs.push_back(Refinement3D(to_vpar, wmin, wmax, umin, umax,
     				   YDIR, deg2+1));
     if (from_wpar-wmin>fuzzy)
       refs.push_back(Refinement3D(from_wpar, umin, umax, vmin, vmax,
     				   ZDIR, deg3+1));
     if (wmax-to_wpar>fuzzy)
       refs.push_back(Refinement3D(to_wpar, umin, umax, vmin, vmax, 
     				   ZDIR, deg3+1));

     // refs[0].setVal(from_upar, vmin, vmax, wmin, wmax, XDIR, deg1+1);
     // refs[1].setVal(to_upar, vmin, vmax, wmin, wmax, XDIR, deg1+1);
     // refs[2].setVal(from_vpar, wmin, wmax, umin, umax, YDIR, deg2+1);
     // refs[3].setVal(to_vpar, wmin, wmax, umin, umax, YDIR, deg2+1);
     // refs[4].setVal(from_wpar, umin, umax, vmin, vmax, ZDIR, deg3+1);
     // refs[5].setVal(to_wpar, umin, umax, vmin, vmax, ZDIR, deg3+1);
     
     // Perform refinement
#ifdef DEBUG
     std::cout << "subVolume. Refine volume" << std::endl;
#endif
     vol2->refine(refs, true);

#ifdef DEBUG
     std::ofstream of("tmp_lr.g2");
     vol2->writeStandardHeader(of);
     vol2->write(of);
#endif

     // Fetch sub mesh
     // First get knot indices
     int iu1 = vol2->mesh().getKnotIdx(XDIR, from_upar, fuzzy);
     int iu2 = vol2->mesh().getKnotIdx(XDIR, to_upar, fuzzy);
     int iv1 = vol2->mesh().getKnotIdx(YDIR, from_vpar, fuzzy);
     int iv2 = vol2->mesh().getKnotIdx(YDIR, to_vpar, fuzzy);
     int iw1 = vol2->mesh().getKnotIdx(ZDIR, from_wpar, fuzzy);
     int iw2 = vol2->mesh().getKnotIdx(ZDIR, to_wpar, fuzzy);
     if (iu1 < 0 || iu2 < 0 || iv1 < 0 || iv2 < 0 || iw1 < 0 || iw2 < 0)
       THROW("LRSplineVolume::subVolume. Knot line not existing");

     shared_ptr<Mesh3D> sub_mesh = vol2->mesh().subMesh(iu1, iu2, iv1, iv2, iw1, iw2);

#ifdef DEBUG
     std::ofstream of2("sub_mesh.g2");
     object_to_stream(of2, *sub_mesh);
     of2 << std::endl;
#endif
     
     // Fetch LR B-splines living on the sub mesh
     vector<LRBSpline3D*> b_splines = vol2->collect_basis(iu1, iu2, iv1, iv2, iw1, iw2);

     //std::cout << "Generating sub volume 2" << std::endl;
     // Create sub volume
     vol = new LRSplineVolume(knot_tol_, rational_, *sub_mesh, b_splines, iu1, iv1, iw1);
#ifdef DEBUG     
     std::ofstream of3("sub_vol.g2");
     vol->writeStandardHeader(of3);
     vol->write(of3);
#endif
    return vol;

}

//==============================================================================
    std::vector<shared_ptr<SplineVolume> > 
      LRSplineVolume::split(std::vector<double>& param,
	    int pardir,
	    double fuzzy) const 
//==============================================================================
{
   MESSAGE("LRSplineVolume:: Not implemented yet.");
   throw;
}

//==============================================================================
     void LRSplineVolume::closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v, 
			      double&        clo_w, 
			      Point&         clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      double   *seed) const
//==============================================================================
{
   MESSAGE("LRSplineVolume:: Not implemented yet.");
   throw;
}

//==============================================================================
    int LRSplineVolume::closestCorner(const Point& pt,
		      double& upar, double& vpar, double& wpar,
		      Point& corner, double& dist) const
//==============================================================================
{
   MESSAGE("LRSplineVolume:: Not implemented yet.");
   throw;
}

//==============================================================================
    void LRSplineVolume::appendVolume(SplineVolume* vol, int join_dir,
		       int cont, bool repar)
//==============================================================================
{
   MESSAGE("LRSplineVolume:: Not implemented yet.");
   throw;
}

//==============================================================================
     void LRSplineVolume::reverseParameterDirection(int pardir)
//==============================================================================
{
   MESSAGE("LRSplineVolume:: Not implemented yet.");
   throw;
}

//==============================================================================
     void LRSplineVolume::swapParameterDirection(int pardir1, int pardir2)
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
void LRSplineVolume::setParameterDomain(double u1, double u2,
                                        double v1, double v2,
                                        double w1, double w2)
//==============================================================================
{
  double umin = paramMin(XDIR);
  double umax = paramMax(XDIR);
  double vmin = paramMin(YDIR);
  double vmax = paramMax(YDIR);
  double wmin = paramMin(ZDIR);
  double wmax = paramMax(ZDIR);

  // All reference to knots are fetched from the mesh_.
  mesh_.setParameterDomain(u1, u2, v1, v2, w1, w2);

  // But in addition the ElementMap emap_ also contains keys which
  // are constructed using the parameters of the elements, as well
  // as Element3D's which stores max and min values in both dirs for
  // the elements.
  // for (ElementMap::iterator iter = elementsBeginNonconst(); iter != elementsEndNonconst(); ++iter)
  // First move the elements out of the container
  vector<unique_ptr<Element3D> > all_elements;
  for (auto it=emap_.begin(); it!=emap_.end(); ++it)
    {
      unique_ptr<Element3D> ptr = std::move(it->second);
      all_elements.emplace_back(std::move(ptr));
    }

  // Empty container
  emap_.clear();

  // Update elements
  for (size_t ki=0; ki<all_elements.size(); ++ki)
    {
      //Element2D* elem = all_elements[ki].get();

      double elem_umin = all_elements[ki]->umin();
      double elem_umax = all_elements[ki]->umax();
      double elem_vmin = all_elements[ki]->vmin();
      double elem_vmax = all_elements[ki]->vmax();
      double elem_wmin = all_elements[ki]->wmin();
      double elem_wmax = all_elements[ki]->wmax();

      double elem_umin_new = (u2 - u1)/(umax - umin)*(elem_umin - umin) + u1;
      double elem_umax_new = (u2 - u1)/(umax - umin)*(elem_umax - umin) + u1;
      double elem_vmin_new = (v2 - v1)/(vmax - vmin)*(elem_vmin - vmin) + v1;
      double elem_vmax_new = (v2 - v1)/(vmax - vmin)*(elem_vmax - vmin) + v1;
      double elem_wmin_new = (w2 - w1)/(wmax - wmin)*(elem_wmin - wmin) + w1;
      double elem_wmax_new = (w2 - w1)/(wmax - wmin)*(elem_wmax - wmin) + w1;
      // We may encounter tolerance issues for the far end of the domain, snap.
      if (fabs(elem_umax_new - u2) < knot_tol_)
        elem_umax_new = u2;
      if (fabs(elem_vmax_new - v2) < knot_tol_)
        elem_vmax_new = v2;
      if (fabs(elem_wmax_new - w2) < knot_tol_)
        elem_wmax_new = w2;
      all_elements[ki]->setUmin(elem_umin_new);
      all_elements[ki]->setUmax(elem_umax_new);
      all_elements[ki]->setVmin(elem_vmin_new);
      all_elements[ki]->setVmax(elem_vmax_new);
      all_elements[ki]->setWmin(elem_wmin_new);
      all_elements[ki]->setWmax(elem_wmax_new);
      all_elements[ki]->updateApproxDataParDomain(elem_umin, elem_umax,
                                                  elem_vmin, elem_vmax,
                                                  elem_wmin, elem_wmax,
                                                  elem_umin_new, elem_umax_new,
                                                  elem_vmin_new, elem_vmax_new,
                                                  elem_wmin_new, elem_wmax_new);

      // Make new key
      ElemKey new_key;
      new_key.u_min = elem_umin_new;
      new_key.v_min = elem_vmin_new;
      new_key.w_min = elem_wmin_new;

      // Insert in container
      // emap_.insert(make_pair(new_key,
      // 		       std::move(unique_ptr<Element2D>(all_elements[ki].get()))));
      emap_.insert(make_pair(new_key, std::move(all_elements[ki])));
  }

  // Must also regenerate keys for the bsplines
  // First move the bsplines out of the container
  vector<unique_ptr<LRBSpline3D> > all_bsplines;
  for (auto it=bsplines_.begin(); it != bsplines_.end(); ++it)
    {
      unique_ptr<LRBSpline3D> ptr = std::move(it->second);
      all_bsplines.emplace_back(std::move(ptr));
    }
  
  bsplines_.clear();
  for (size_t ki=0; ki<all_bsplines.size(); ++ki)
    {
      auto key = generate_key(*all_bsplines[ki], mesh_);
      bsplines_.insert(make_pair(key, std::move(all_bsplines[ki])));
    }
  
}

//==============================================================================
   double LRSplineVolume::volume(double tol) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
   void 
  LRSplineVolume::getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
   LRSplineSurface* LRSplineVolume::constParamSurface(double parameter,
					     int pardir) const
//==============================================================================
{
    // We do it the lazy (and time consuming) way by extracting the
    // subvolume and fetching the edge surface.
    int bd_num = -1; // bd_num: 0 = umin, 1 = umax, 2 = vmin, 3 = vmax, 4 = wmin, 5 = wmax
    if (pardir == 0 && fabs(parameter - startparam_u()) < knot_tol_)
      bd_num = 0;
    else if (pardir == 0 && fabs(parameter - endparam_u()) < knot_tol_)
      bd_num = 1;
    else if (pardir == 1 && fabs(parameter - startparam_v()) < knot_tol_)
      bd_num = 2;
    else if (pardir == 1 && fabs(parameter - endparam_v()) < knot_tol_)
      bd_num = 3;
    else if (pardir == 2 && fabs(parameter - startparam_w()) < knot_tol_)
      bd_num = 4;
    else if (pardir == 2 && fabs(parameter - endparam_w()) < knot_tol_)
      bd_num = 5;
    bool par_on_edge = (bd_num != -1);

    if (par_on_edge)
      {
	return boundarySurface(bd_num);
      }
    else
      {
	// Snap if the given parameter is very close to an existing knot
	Direction3D dir = (pardir == 0) ? XDIR : ((pardir == 1) ? YDIR : ZDIR);
	int ixp = mesh_.knotIntervalFuzzy(dir, parameter, knot_tol_);

	// Define refinement
	int deg = degree(dir);
	vector<Refinement3D> refs;
	Direction3D dir1 = next(dir);
	Direction3D dir2 = prev(dir);
	refs.push_back(Refinement3D(parameter, paramMin(dir1), paramMax(dir1),
				    paramMin(dir2), paramMax(dir2), dir, deg+1));
	
	// Make a copy of the current volume
	shared_ptr<LRSplineVolume> vol2(new LRSplineVolume(*this));

	// Refine without defining element structure
	vol2->doRefine(refs, true);

	// Construct constant parameter surface
	return vol2->defineBivariate(pardir, parameter, false);
	
	// std::cout << "Generating sub volume" << std::endl;
	// shared_ptr<LRSplineVolume> sub_vol(subVolume(umin, vmin, wmin,
	// 					     umax, vmax, wmax,
	// 					     knot_tol_));
     
	// int sub_bd_num = 2*pardir + 1;
	// return sub_vol->boundarySurface(sub_bd_num);
      }
}

//==============================================================================

void LRSplineVolume::constParamSurfaces(vector<double>& param, int pardir,
					vector<shared_ptr<LRSplineSurface> >& sfs) const
//==============================================================================
{
  Direction3D dir = (pardir == 0) ? XDIR : ((pardir == 1) ? YDIR : ZDIR);
  int deg = degree(dir);
  Direction3D dir1 = next(dir);
  Direction3D dir2 = prev(dir);

  double start = paramMin(dir);
  
  // Define refinements
  vector<Refinement3D> refs;
  for (size_t ki=0; ki<param.size(); ++ki)
    {
      // Snap if the given parameter is very close to an existing knot
      int ixp = mesh_.knotIntervalFuzzy(dir, param[ki], knot_tol_);

	refs.push_back(Refinement3D(param[ki], paramMin(dir1), paramMax(dir1),
				    paramMin(dir2), paramMax(dir2), dir, deg+1));
    }
  
  // Make a copy of the current volume
  shared_ptr<LRSplineVolume> vol2(new LRSplineVolume(*this));

  std::cout << "Redy to refine" << std::endl;
  // Refine without defining element structure
  vol2->doRefine(refs, true);
  std::cout << "Refinement performed" << std::endl;
  
  for (size_t ki=0; ki<param.size(); ++ki)
    {
      // Construct constant parameter surface
      shared_ptr<LRSplineSurface> surf(vol2->defineBivariate(pardir, param[ki], param[ki] == start));
      sfs.push_back(surf);
	
	// std::cout << "Generating sub volume" << std::endl;
	// shared_ptr<LRSplineVolume> sub_vol(subVolume(umin, vmin, wmin,
	// 					     umax, vmax, wmax,
	// 					     knot_tol_));
     
	// int sub_bd_num = 2*pardir + 1;
	// return sub_vol->boundarySurface(sub_bd_num);
      }
}

//==============================================================================
LRSplineSurface* 
LRSplineVolume::boundarySurface(int bd_num) const
//==============================================================================
{
  // Assumes degree + 1 multiplicity at boundary
  int pardir = bd_num/2;
  Direction3D d = (pardir == 0) ? XDIR : ((pardir == 1) ? YDIR : ZDIR);
  double parval = (bd_num % 2 == 0) ? paramMin(d) : paramMax(d);

  return defineBivariate(pardir, parval, bd_num % 2 == 0);
}

//==============================================================================
LRSplineSurface* 
LRSplineVolume::defineBivariate(int pardir, double parval, bool at_start) const
//==============================================================================
{
  double eps = 1.0e-8;
  Direction3D d = (pardir == 0) ? XDIR : ((pardir == 1) ? YDIR : ZDIR);
  int deg = degree(d);
  int sgn = 1;

  // Identify bsplines uni with the multiplicity equal to order at the given parameter
  vector<BSplineUniLR*> uni;
  vector<std::unique_ptr<BSplineUniLR> >::const_iterator start, end;
  if (pardir == 0)
    {
      start = bsplinesuni1_.begin();
      end = bsplinesuni1_.end();
    }
  else if (pardir == 1)
    {
      start = bsplinesuni2_.begin();
      end = bsplinesuni2_.end();
    }
  else
    {
      start = bsplinesuni3_.begin();
      end = bsplinesuni3_.end();
    }
  if (at_start)
    {
      for (auto it=start; it!=end; ++it)
	{
	  double startval = (*it)->min();
	  if (startval > parval)
	    break;
	  if (startval == parval)
	    {
	      // Check multiplicity
	      int mult = (*it)->endmult(true);
	      if (mult == deg+1)
		uni.push_back((*it).get());
	    }
	}
    }
  else
    {
      for (auto it=end; it!=start; --it)
	{
	  auto it2 = it;
	  it2--;
	  double endval = (*it2)->max();
	  if (endval == parval)
	    {
	      // Check multiplicity
	      int mult = (*it2)->endmult(false);
	      if (mult == deg+1)
		uni.push_back((*it2).get());
	    }
	}
      sgn = -1;
    }

  // Construct bivariate mesh
  Direction3D d1 = next(d); //(d == XDIR) ? YDIR : ((d == YDIR) ? ZDIR : XDIR);
  Direction3D d2 = prev(d); //(d == XDIR) ? ZDIR : ((d == YDIR) ? XDIR : YDIR);
  vector<double> knots1 = mesh_.allKnots(d1);
  vector<double> knots2 = mesh_.allKnots(d2);
  int ix = mesh_.getKnotIdx(d, parval, eps);
  vector<vector<int> > mrvec1, mrvec2;
  Mesh3DUtils::sectionKnots(mesh_, d, ix, mrvec2, mrvec1);
  Mesh2D mesh2d(knots1, knots2, mrvec1, mrvec2);
  int stop_break = 1;
  
  // Identify B-splines associated to the univariate B-splines and construct
  // bivariate counter part
  vector<std::unique_ptr<BSplineUniLR> > b_uni1, b_uni2;
  vector<LRBSpline2D*> bd_bsplines;
  int left1 = 0, left2 = 0;
  for (BSplineMap::const_iterator it = basisFunctionsBegin();
       it != basisFunctionsEnd(); ++it)
    {
      BSplineUniLR *d_uni = it->second->getUnivariate(d);
      if (ix < d_uni->suppMin() || ix > d_uni->suppMax())
	continue;
      
      // Check if the current B-spline is associated with the identified
      // univariate B-splines
      size_t kr;
      for (kr=0; kr<uni.size(); ++kr)
	if (*uni[kr] == *d_uni)
	  break;

      if (kr == uni.size())
	continue;  // Not a boundary B-spline

      BSplineUniLR *d1_uni = it->second->getUnivariate(d1);
      BSplineUniLR *d2_uni = it->second->getUnivariate(d2);
     
      bool found1 = BSplineUniUtils::identify_bsplineuni(d1_uni, b_uni1, left1);
      if (!found1)
	{
	  BSplineUniLR *tmp = new BSplineUniLR(*d1_uni);
	  tmp->setPardir(1);
	  BSplineUniUtils::insert_univariate(b_uni1, tmp, left1);
	}
      b_uni1[left1]->incrCount();
      
      bool found2 = BSplineUniUtils::identify_bsplineuni(d2_uni, b_uni2, left2);
      if (!found2)
	{
	  BSplineUniLR *tmp = new BSplineUniLR(*d2_uni);
	  tmp->setPardir(2);
	  BSplineUniUtils::insert_univariate(b_uni2, tmp, left2);
	}
      b_uni2[left2]->incrCount();

      Point cfg = it->second->coefTimesGamma();
      double gamma = it->second->gamma();
      double weight = it->second->weight();
#ifdef DEBUG
      if (fabs(gamma-1.0) > 1.0e-6)
	std::cout << "Scaling factor: " << gamma << std::endl;
#endif
      LRBSpline2D *bspline = new LRBSpline2D(cfg, weight, b_uni1[left1].get(),
					     b_uni2[left2].get(), gamma,
					     it->second->rational());
      bd_bsplines.push_back(bspline);
    }

  LRSplineSurface *bd_surf = new LRSplineSurface(knot_tol_, rational_, mesh2d,
						 bd_bsplines, 0, 0);
  for (size_t ki=0; ki<bd_bsplines.size(); ++ki)
    delete bd_bsplines[ki];
  
  return bd_surf;
}

//==============================================================================
   std::vector<shared_ptr<ParamSurface> > 
  LRSplineVolume::getAllBoundarySurfaces() const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  bool LRSplineVolume::isDegenerate(int which_sf, int& type, bool& b, bool& r,
		    bool& t, bool& l, double tol) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
   void LRSplineVolume::getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
   double LRSplineVolume::nextSegmentVal(int dir, double par, bool forward, double tol) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
Point LRSplineVolume::operator()(double u, double v, double w, int u_deriv, int v_deriv, int w_deriv) const // evaluation
//==============================================================================
{
  // const bool u_on_end = (u == mesh_.maxParam(XFIXED));
  // const bool v_on_end = (v == mesh_.maxParam(YFIXED));
  // vector<LRBSpline2D*> covering_B_functions = 
  //   basisFunctionsWithSupportAt(u, v);
  Element3D* elem;
  if (curr_element_ && curr_element_->contains(u, v, w))
    elem = curr_element_;
  else
    {
      //std::cout << "Finding element for parameter value (" << u << "," << v << ")" << std::endl;
      elem = coveringElement(u, v, w);
      curr_element_ = (Element3D*)elem;
    }
  const vector<LRBSpline3D*> covering_B_functions = elem->getSupport();

  Point result(this->dimension()); 
  result.setValue(0.0); // will be initialized to 0, with the correct dimension

  // loop over LR B-spline functions
  int ki=0;
  int nmb_b = (int)covering_B_functions.size();
  double denom = (rational_) ? 0.0 : 1.0;
  double denom_pos = 0.0;
  double denom_der = 0.0;
  Point nom_pos(this->dimension());
  nom_pos.setValue(0.0);
  Point nom_der(this->dimension());
  nom_der.setValue(0.0);
  vector<double> basis_vals((u_deriv+1)*(v_deriv+1)*(w_deriv+1), 0.0); // To be used for rational cases, needed for derivs.
  double eps = 1.0e-12;
  const bool u_on_end = (u >= mesh_.maxParam(XDIR)-eps); //(u == (*b)->umax());
  const bool v_on_end = (v >= mesh_.maxParam(YDIR)-eps); // (v == (*b)->vmax());
  const bool w_on_end = (w >= mesh_.maxParam(ZDIR)-eps); 
  for (auto b = covering_B_functions.begin(); 
       b != covering_B_functions.end(); ++b, ++ki) 
    {
      // The b-function contains the coefficient.
      if (!rational_)
	{
	  result += (*b)->eval(u, 
			       v, 
			       w,
			       u_deriv, 
			       v_deriv, 
			       w_deriv, 
			       u_on_end, 
			       v_on_end, 
			       w_on_end);
	}
      else
	{
	  
//	  for (size_t ki = 0; ki < 
	  double basis_val_pos = (*b)->evalBasisFunction(u, 
							 v, 
							 w,
							 0, 
							 0, 
							 0,
							 u_on_end, 
							 v_on_end,
							 w_on_end);
	  double gamma = (*b)->gamma();
	  double weight = (*b)->weight();
	  Point coef = (*b)->coefTimesGamma();

	  // This is the nominator-position.
	  nom_pos += coef*weight*basis_val_pos;

	  denom_pos += weight*basis_val_pos;

	  if (u_deriv > 0 || v_deriv > 0 || w_deriv > 0)
	    {
	      double basis_val_der = (*b)->evalBasisFunction(u, 
							     v, 
							     w, 
							     u_deriv, 
							     v_deriv, 
							     w_deriv, 
							     u_on_end, 
							     v_on_end,
							     w_on_end);

	      // This is the nominator-deriv.
	      nom_der += coef*weight*basis_val_der;

	      denom_der += weight*basis_val_der;
	    }

#ifndef NDEBUG
	  if (u_deriv > 0 || v_deriv > 0 || w_deriv > 0)
	    {
	      ;//MESSAGE("Do not think that rational derivs are supported yet.");
//	      denom = 1.0;
	    }
//	  std::cout << "denom: " << denom << std::endl;
	  // if (rat_den == 0.0)
	  //   rat_den = 1.0;
#endif

	}
    }

  if (rational_)
    {
      if (u_deriv == 0 && v_deriv == 0 && w_deriv == 0)
	{
	  result = nom_pos/denom_pos;
	}
      else if (u_deriv + v_deriv + w_deriv == 1)
	{
	  result = (nom_der*denom_pos - denom_der*nom_pos)/(denom_pos*denom_pos);
	}
      else
	{
	  std::cout << "Not supported." << std::endl;
	}
    }

  return result;
}


//==============================================================================
  Point LRSplineVolume::operator()(double u, double v, double w,
                                   int u_deriv, int v_deriv, int w_deriv,
                                   Element3D* elem) const
//==============================================================================
{
  // Check element
  if (!elem || !elem->contains(u, v, w))
    {
      bool found = false;

      // Check neighbours
      if (elem)
        {
          vector<LRBSpline3D*> bsupp = elem->getSupport();
          std::set<Element3D*> supp_el;

          for (size_t ka=0; ka<bsupp.size(); ++ka)
            {
              vector<Element3D*> esupp = bsupp[ka]->supportedElements();
              supp_el.insert(esupp.begin(), esupp.end());
            }

          vector<Element3D*> supp_el2(supp_el.begin(), supp_el.end());

          for (size_t ka=0; ka<supp_el2.size(); ++ka)
            {
              if (supp_el2[ka]->contains(u, v, w))
                {
                  elem = supp_el2[ka];
                  found = true;
                  break;
                }
            }
        }

      if (!found)
          {
            //std::cout << "Finding element for parameter value (" << u << "," << v << ")" << std::endl;
            elem = coveringElement(u, v, w);
          }
      }

    const vector<LRBSpline3D*>& covering_B_functions = elem->getSupport();

    Point result(this->dimension());
    result.setValue(0.0); // will be initialized to 0, with the correct dimension

    // loop over LR B-spline functions
    int ki=0;

    // Distinguish between rational and non-rational to avoid
    // making temporary storage in the non-rational case
    double eps = 1.0e-12;
    const bool u_on_end = (u >= mesh_.maxParam(XDIR)-eps); //(u == (*b)->umax());
    const bool v_on_end = (v >= mesh_.maxParam(YDIR)-eps); // (v == (*b)->vmax());
    const bool w_on_end = (w >= mesh_.maxParam(ZDIR)-eps);
    if (!rational_)
      {
        for (auto b = covering_B_functions.begin();
             b != covering_B_functions.end(); ++b, ++ki)
          {
            // The b-function contains the coefficient.
            result += (*b)->eval(u,
                                 v,
                                 w,
                                 u_deriv,
                                 v_deriv,
                                 w_deriv,
                                 u_on_end,
                                 v_on_end,
                                 w_on_end);
          }
      }
    else
      {
      double denom_pos = 0.0;
      double denom_der = 0.0;
      Point nom_pos(this->dimension());
      nom_pos.setValue(0.0);
      Point nom_der(this->dimension());
      nom_der.setValue(0.0);
      for (auto b = covering_B_functions.begin();
           b != covering_B_functions.end(); ++b, ++ki)
        {
          const bool u_on_end = (u == (*b)->umax());
          const bool v_on_end = (v == (*b)->vmax());
          const bool w_on_end = (w == (*b)->wmax());

          // The b-function contains the coefficient.
          double basis_val_pos = (*b)->evalBasisFunction(u,
                                                         v,
                                                         w,
                                                         0,
                                                         0,
                                                         0,
                                                         u_on_end,
                                                         v_on_end,
                                                         w_on_end);
          //double gamma = (*b)->gamma();
          double weight = (*b)->weight();
          Point coef = (*b)->coefTimesGamma();

          // This is the nominator-position.
          nom_pos += coef*weight*basis_val_pos;

          denom_pos += weight*basis_val_pos;

          if (u_deriv > 0 || v_deriv > 0 || w_deriv > 0)
            {
              double basis_val_der = (*b)->evalBasisFunction(u,
                                                             v,
                                                             w,
                                                             u_deriv,
                                                             v_deriv,
                                                             w_deriv,
                                                             u_on_end,
                                                             v_on_end,
                                                             w_on_end);

              // This is the nominator-deriv.
              nom_der += coef*weight*basis_val_der;

              denom_der += weight*basis_val_der;
            }

  #ifndef NDEBUG
          if (u_deriv > 0 || v_deriv > 0 || w_deriv > 0)
            {
              ;//MESSAGE("Do not think that rational derivs are supported yet.");
                //	      denom = 1.0;
            }
            //	  std::cout << "denom: " << denom << std::endl;
            // if (rat_den == 0.0)
            //   rat_den = 1.0;
  #endif

        }

      if (u_deriv == 0 && v_deriv == 0 && w_deriv == 0)
        {
          result = nom_pos/denom_pos;
        }
      else if (u_deriv + v_deriv + w_deriv == 1)
        {
          result = (nom_der*denom_pos - denom_der*nom_pos)/(denom_pos*denom_pos);
        }
      else
        {
          std::cout << "Not supported." << std::endl;
        }
    }

    return result;
}


//==============================================================================
  void LRSplineVolume::computeBasis (double param_u, double param_v, double param_w,
				     BasisPtsSf& result, int iEl ) const
//==============================================================================
{
   MESSAGE("LRSplineVolume:: Not implemented yet.");
   throw;
}

//==============================================================================
  void LRSplineVolume::computeBasis (double param_u, double param_v, double param_w,
				     BasisDerivsSf& result, int iEl ) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  void LRSplineVolume::computeBasis (double param_u, double param_v, double param_w,
				     BasisDerivsSf2& result, int iEl ) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  void LRSplineVolume::computeBasis (double param_u, double param_v, double param_w,
				     std::vector<std::vector<double> >& result,
				     int derivs,
				     int iEl ) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  int LRSplineVolume::getElementContaining(double u, double v, double w) const
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
Go::Element3D* LRSplineVolume::coveringElement(double u, double v, double w) const
//==============================================================================
{
  // First check if a current element exists
  if (curr_element_)
    {
      if (curr_element_->contains(u, v, w))
	return curr_element_;
      
      // Check neighbours
      vector<LRBSpline3D*> bsupp = curr_element_->getSupport();
      std::set<Element3D*> supp_el;
      for (size_t ka=0; ka<bsupp.size(); ++ka)
	{
	  vector<Element3D*> esupp = bsupp[ka]->supportedElements();
	  for (size_t kb=0; kb<esupp.size(); ++kb)
	    {
	      if (esupp[kb]->contains(u, v, w))
		{
		  curr_element_ = esupp[kb];
		  return curr_element_;
		}
 	    }
	}
    }

  int ucorner, vcorner, wcorner;
  if (! Mesh3DUtils::identify_patch_lower_left(mesh_, u, v, w, ucorner, vcorner, wcorner) )
  {
#ifndef NDEBUG
    std::cout << "u: " << u << ", v: " << v << ", w: " << w << std::endl;
#endif
    THROW("Parameter outside domain in LRSplineVolume::coveringElement()");
  }

#if 0//ndef NDEBUG
  {
    vector<LRBSpline3D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
        bas_funcs.push_back((*iter).second.get());
      }
    //    puts("coveringElement(): Remove when done debugging!");
    vector<Element3D*> elems;
    vector<ElemKey> elem_keys;
    for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
    {
      elems.push_back(((*iter).second.get()));
      elem_keys.push_back(iter->first);
    }
//     puts("coveringElement(): Remove when done debugging!");
  }
#endif

  const LRSplineVolume::ElemKey key =
    {mesh_.knotsBegin(XDIR)[ucorner], mesh_.knotsBegin(YDIR)[vcorner], mesh_.knotsBegin(ZDIR)[wcorner]};
  const auto el = emap_.find(key);
  if (el == emap_.end())
    {
      vector<Element3D*> elems;
      vector<ElemKey> elem_keys;
      for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
	{
	  elems.push_back(((*iter).second.get()));
	  elem_keys.push_back(iter->first);
	}

      THROW("Current element not found");
    }
  curr_element_ = el->second.get();
  return el->second.get();
}


//==============================================================================
 void LRSplineVolume::constructElementMesh(vector<Element3D*>& elements) const
//==============================================================================
{
  double eps = 1.0e-12;

  std::ofstream of("element_domain.txt");
  
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XDIR);
  int nmb_knots_u = mesh_.numDistinctKnots(XDIR);

  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YDIR);
  int nmb_knots_v = mesh_.numDistinctKnots(YDIR);

  // Get all knot values in the w-direction
  const double* const wknots = mesh_.knotsBegin(ZDIR);
  int nmb_knots_w = mesh_.numDistinctKnots(ZDIR);

  // Construct mesh of element pointers
  int ki, kj, kk, kr;
  elements.resize((nmb_knots_u-1)*(nmb_knots_v-1)*(nmb_knots_w-1), NULL);
  int num = numElements();
  ElementMap::const_iterator it;
  for (it=elementsBegin(), kr=0, kj=0, ki=0, kk=0; kr<num; ++it, ++kr)
    {
      double umin = it->second->umin();
      double umax = it->second->umax();
      double vmin = it->second->vmin();
      double vmax = it->second->vmax();
      double wmin = it->second->wmin();
      double wmax = it->second->wmax();
      of << kr << ": (" << umin << "," << umax << ") x (" << vmin;
      of << "," << vmax << ") x (" << wmin << "," << wmax << std::endl;

      for (kk=0; kk<nmb_knots_w && wmin > wknots[kk]+eps; ++kk);
      for (kj=0; kj<nmb_knots_v && vmin > vknots[kj]+eps; ++kj);
      for (ki=0; ki<nmb_knots_u && umin > uknots[ki]+eps; ++ki);
      for (int kh1=kk; kh1<nmb_knots_w-1 && wmax >= wknots[kh1+1]-eps; ++kh1)
	for (int kh2=kj; kh2<nmb_knots_v-1 && vmax >= vknots[kh2+1]-eps; ++kh2)
	  for (int kh3=ki; kh3<nmb_knots_u-1 && umax >= uknots[kh3+1]-eps; ++kh3)
	  {
	    if (it->second.get() == NULL)
	      {
		std::cout << "Null element: (" << umin << "," << umax << ") x (";
		std::cout << vmin << "," << vmax << ") x (" << wmin;
		std::cout << "," << wmax << ")" << std::endl;
	      }
	    elements[(kh1*(nmb_knots_v-1)+kh2)*(nmb_knots_u-1)+kh3] = it->second.get();
	  }
      
     }
}

//==============================================================================
  std::vector<LRBSpline3D*>
  LRSplineVolume::basisFunctionsWithSupportAt(double u, double v, double w) const
//==============================================================================
{
  vector<LRBSpline3D*> support_functions;
  auto elem = coveringElement(u, v, w);
  vector<LRBSpline3D*>::const_iterator first = elem->supportBegin();
  vector<LRBSpline3D*>::const_iterator last = elem->supportEnd();
  int ki=0;
  for (; first != last; ++first, ++ki)
    {
      support_functions.push_back(*first);
    }
  
  return support_functions;
}

//==============================================================================
  bool LRSplineVolume::isFullTensorProduct() const
//==============================================================================
{
    return (LRSpline3DUtils::all_meshlines_uniform(XDIR, mesh_) &&
            LRSpline3DUtils::all_meshlines_uniform(YDIR, mesh_) &&
            LRSpline3DUtils::all_meshlines_uniform(ZDIR, mesh_));
}

//==============================================================================
  void LRSplineVolume::refine(Direction3D d, double fixed_val,
			      double start1, double end1,
			      double start2, double end2,
			      int mult, bool absolute)
//==============================================================================
{
  // vector<LRBSpline3D*> bvec0;
  // for (auto it2=bsplines_.begin(); it2!=bsplines_.end(); ++it2)
  //   {
  //     bvec0.push_back(it2->second.get());
  //   }
  // vector<Element3D*> evec0;
  // for (auto it3=emap_.begin(); it3!=emap_.end(); ++it3)
  //   {
  //     evec0.push_back(it3->second.get());
  //   }

  // Make a copy of the initial mesh
  Mesh3D mesh2 = mesh_;

  bool refined;
  const auto indices = // tuple<int, int, int, int, int, int>
  LRSpline3DUtils::refine_mesh(d, fixed_val,
                               start1, end1,
                               start2, end2,
                               mult, absolute,
                               degree(d), knot_tol_, mesh_, 
			       (d == XDIR) ? bsplinesuni1_ : 
			       (d == YDIR) ? bsplinesuni2_ : bsplinesuni3_,
			       refined);

  if (!refined)
    return;  // Mesh rectangle already existing
  
  // insert newly created elements to emap (unless refinement was on border, in which case no new element
  // could possibly be created
  const int prev_ix = get<0>(indices);  // Index of last non-larger in fixed direction
  const int fixed_ix = get<1>(indices); // Index of fixed_val in global knot vector.
  const int start_ix1 = get<2>(indices); // Index of start (of segment to insert) in global knot vector.
  const int end_ix1   = get<3>(indices); // Index of end (of segment to insert) in global knot vector.
  const int start_ix2 = get<4>(indices); // Index of start (of segment to insert) in global knot vector.
  const int end_ix2   = get<5>(indices); // Index of end (of segment to insert) in global knot vector.

  // Collect pointers to affected bsplines
  std::set<LRBSpline3D*> all_bsplines;
  double domain[6];  // Covers elements affected by the split
  domain[0] = domain[2] = domain[4] = std::numeric_limits<double>::max();
  domain[1] = domain[3] = domain[5] = std::numeric_limits<double>::lowest();
  for (int i = start_ix1; i != end_ix1; ++i)
    {
      for (int j = start_ix2; j != end_ix2; ++j)
        {
          // Check if the specified element exists in 'emap'
          int u_ix = (d == XDIR) ? prev_ix : ((d == YDIR) ? j : i);
          int v_ix = (d == YDIR) ? prev_ix : ((d == ZDIR) ? j : i);
          int w_ix = (d == ZDIR) ? prev_ix : ((d == XDIR) ? j : i);
          ElementMap::key_type key = {mesh_.kval(XDIR, u_ix),
                                      mesh_.kval(YDIR, v_ix),
                                      mesh_.kval(ZDIR, w_ix)};
          auto it = emap_.find(key);
          if (it == emap_.end())
            {
              // Element not found. The assumed start index of the element
              // is not correct. Recompute
              int u_ix2 = u_ix;
              int v_ix2 = v_ix;
              int w_ix2 = w_ix;
              if (d == XDIR)
                u_ix2 =
                    Mesh3DUtils::search_downwards_for_nonzero_multiplicity(mesh2, XDIR,
                                                                           u_ix, v_ix, w_ix);
              else if (d == YDIR)
                v_ix2 =
                    Mesh3DUtils::search_downwards_for_nonzero_multiplicity(mesh2, YDIR,
                                                                           v_ix, w_ix, u_ix);
              else
                w_ix2 =
                    Mesh3DUtils::search_downwards_for_nonzero_multiplicity(mesh2, ZDIR,
                                                                           w_ix, u_ix, v_ix);
              u_ix = u_ix2;
              v_ix = v_ix2;
              w_ix = w_ix2;

              key = {mesh2.kval(XDIR, u_ix), mesh2.kval(YDIR, v_ix), mesh2.kval(ZDIR, w_ix)};
              it = emap_.find(key);
            }

          if (it != emap_.end())
            {
              // The element exists. Collect bsplines
              Element3D* curr_el = (*it).second.get();
              all_bsplines.insert(curr_el->supportBegin(), curr_el->supportEnd());
              domain[0] = std::min(domain[0], curr_el->umin());
              domain[1] = std::max(domain[1], curr_el->umax());
              domain[2] = std::min(domain[2], curr_el->vmin());
              domain[3] = std::max(domain[3], curr_el->vmax());
              domain[4] = std::min(domain[4], curr_el->wmin());
              domain[5] = std::max(domain[5], curr_el->wmax());
            }
        }
    }
  vector<LRBSpline3D*> bsplines_affected(all_bsplines.begin(), all_bsplines.end());

  if (bsplines_affected.size() == 0)
    {
      //mesh_ = mesh2; // @@@obar should we revert the mesh to its unrefined version?
      return;  // No B-splines will be split, no knot insertion is possible
    }

    // Split univariate to prepare for bivariate split
    int last_ix = 
      BSplineUniUtils::last_overlapping_bsplineuni(fixed_ix,
						   (d == XDIR) ? bsplinesuni1_ : 
						   (d == YDIR) ? bsplinesuni2_ : bsplinesuni3_);

    if (d == XDIR)
      {
	LRSplineUtils::split_univariate(bsplinesuni1_, last_ix, fixed_ix, 
					absolute ? mult : 1);
      }
    else if (d == YDIR)
      {
	LRSplineUtils::split_univariate(bsplinesuni2_, last_ix, fixed_ix, 
					absolute ? mult : 1);
       }
    else 
      {
	LRSplineUtils::split_univariate(bsplinesuni3_, last_ix, fixed_ix, 
					absolute ? mult : 1);
       }

  // Cannot remove the bsplines from the global array at this stage since we operate
  // with pointers to it. When a bspline is split, the origin is removed from the
  // array after all pointers are updated and the the bspline is allowed to die.
  // Iteratively split affected LRBSpline3Ds
  // @@@ VSK. Will pointers to other entities in bsplines_ which are not
  // affected remain valid after removing and adding elements? If not, this
  // combination of objects and pointers will not work.
    LRSpline3DUtils::iteratively_split2(bsplines_affected, mesh_, bsplines_, domain,
					bsplinesuni1_, bsplinesuni2_, bsplinesuni3_);

    // Remove unused univariate B-splines
    if (d == XDIR)
      {
	for (int i=last_ix; i>=0; --i)
	  if (bsplinesuni1_[i]->getCount() <= 0)
	    bsplinesuni1_.erase(bsplinesuni1_.begin()+i);
      }
    else if (d == YDIR)
      {
	for (int i=last_ix; i>=0; --i)
	  if (bsplinesuni2_[i]->getCount() <= 0)
	    bsplinesuni2_.erase(bsplinesuni2_.begin()+i);
      }
    else
      {
	for (int i=last_ix; i>=0; --i)
	  if (bsplinesuni3_[i]->getCount() <= 0)
	    bsplinesuni3_.erase(bsplinesuni3_.begin()+i);
      }

  // Presumably this is updating the element map???
  if (fixed_ix > 0 && fixed_ix != mesh_.numDistinctKnots(d)-1)
    {
      for (int i = start_ix1; i != end_ix1; ++i)
        {
          for (int j = start_ix2; j != end_ix2; ++j)
            {
              if ( (mesh_.nu(next(d), i, j, j+1, fixed_ix, fixed_ix+1) > 0)         // @obar: SEEMS OK, BUT CHECK
                   && (mesh_.nu(prev(d), j, fixed_ix, fixed_ix+1, i, i+1) > 0) )    // @obar: SEEMS OK, BUT CHECK
                {
                  // this is the minimum corner of an element bordering our refinement.
                  // Check if it already exists in 'emap', if not, insert it.
                  // Also modify the current element and update bspline pointers
                  // in the elements
                  int u_ix2 = (d == XDIR) ? prev_ix : ((d == YDIR) ? j : i);
                  int v_ix2 = (d == YDIR) ? prev_ix : ((d == ZDIR) ? j : i);
                  int w_ix2 = (d == ZDIR) ? prev_ix : ((d == XDIR) ? j : i);
                  ElementMap::key_type key2 = {mesh_.kval(XDIR, u_ix2),
                                               mesh_.kval(YDIR, v_ix2),
                                               mesh_.kval(ZDIR, w_ix2)};
                  auto it2 = emap_.find(key2);

                  if (it2 == emap_.end())
                    {
                      // Element not found. The assumed start index of the element
                      // is not correct. Recompute
                      int u_ix3 = u_ix2;
                      int v_ix3 = v_ix2;
                      int w_ix3 = w_ix2;
                      if (d == XDIR)
                        u_ix3 =
                            Mesh3DUtils::search_downwards_for_nonzero_multiplicity(mesh2, XDIR,
                                                                                   u_ix2, v_ix2, w_ix2);
                      else if (d == YDIR)
                        v_ix3 =
                            Mesh3DUtils::search_downwards_for_nonzero_multiplicity(mesh2, YDIR,
                                                                                   v_ix2, w_ix2, u_ix2);
                      else
                        w_ix3 =
                            Mesh3DUtils::search_downwards_for_nonzero_multiplicity(mesh2, ZDIR,
                                                                                   w_ix2, u_ix2, v_ix2);
                      u_ix2 = u_ix3;
                      v_ix2 = v_ix3;
                      w_ix2 = w_ix3;

                      key2 = {mesh2.kval(XDIR, u_ix2), mesh2.kval(YDIR, v_ix2), mesh2.kval(ZDIR, w_ix2)};
                      it2 = emap_.find(key2);
                    }

                  int u_ix = (d == XDIR) ? fixed_ix : ((d == YDIR) ? j : i);
                  int v_ix = (d == YDIR) ? fixed_ix : ((d == ZDIR) ? j : i);
                  int w_ix = (d == ZDIR) ? fixed_ix : ((d == XDIR) ? j : i);
                  ElementMap::key_type key = {mesh_.kval(XDIR, u_ix),
                                              mesh_.kval(YDIR, v_ix),
                                              mesh_.kval(ZDIR, w_ix)};
                  auto it = emap_.find(key);

                  vector<double> data_points;
                  bool sort_in_u;
                  //double maxerr, averr, accerr;
                  int nmbout;

                  if (it2 != emap_.end())
                    {
                      // Update size of existing element
                      Mesh3DIterator m(mesh_, u_ix2, v_ix2, w_ix2);
                      // Indices below refer to the upper corner
                      it2->second->setUmax(mesh_.kval(XDIR, (*m)[3]));
                      it2->second->setVmax(mesh_.kval(YDIR, (*m)[4]));
                      it2->second->setWmax(mesh_.kval(ZDIR, (*m)[5]));

                      // Fetch scattered data from the element that no longer is
                      // inside
                      it2->second->getOutsidePoints(data_points, d, sort_in_u);
                      // it2->second->getAccuracyInfo(averr, maxerr, nmbout);
                      // accerr = it2->second->getAccumulatedError();

                      // Update accuracy statistices in element
                      it2->second->updateAccuracyInfo();

                      // Update supported LRBsplines
                      for (size_t kb=0; kb<bsplines_affected.size(); ++kb)
                        {
                          if (!bsplines_affected[kb]->overlaps(it2->second.get()))
                            {
                              it2->second->removeSupportFunction(bsplines_affected[kb]);
                              bsplines_affected[kb]->removeSupport(it2->second.get());
                            }
                          else
                            it2->second->addSupportFunction(bsplines_affected[kb]);
                        }
                    }

                  // Create new element
                  if (it == emap_.end())
                    {
                      Mesh3DIterator m(mesh_, u_ix, v_ix, w_ix);
                      unique_ptr<Element3D> elem(new Element3D(
                                                   mesh_.kval(XDIR, (*m)[0]),
                                                   mesh_.kval(YDIR, (*m)[1]),
                                                   mesh_.kval(ZDIR, (*m)[2]),
                                                   mesh_.kval(XDIR, (*m)[3]),
                                                   mesh_.kval(YDIR, (*m)[4]),
                                                   mesh_.kval(ZDIR, (*m)[5])) );

                      // Set LRBsplines
                      for (size_t kb=0; kb<bsplines_affected.size(); ++kb)
                        {
                          if (bsplines_affected[kb]->overlaps(elem.get()))
                            {
                              elem->addSupportFunction(bsplines_affected[kb]);
                              bsplines_affected[kb]->addSupport(elem.get());
                            }
                        }

                      // Store data points in the element
                      if (data_points.size() > 0)
                        elem->addDataPoints(data_points.begin(), data_points.end(),
                                            sort_in_u);
                      //if (ghost_points.size() > 0)
                      //elem->addGhostPoints(ghost_points.begin(), ghost_points.end(),
                      //sort_in_u_ghost);
                      //elem->setAccuracyInfo(accerr, averr, maxerr, nmbout);  // Not exact info as the
                      // element has been split
                      elem->updateAccuracyInfo();  // Accuracy statistic in element

                      emap_.insert(std::make_pair(key, std::move(elem)));
                      //auto it3 = emap_.find(key);

                    }
                }
            }
        }
    }

  // // Check
  // vector<LRBSpline3D*> bvec;
  // for (auto it2=bsplines_.begin(); it2!=bsplines_.end(); ++it2)
  //   {
  //     bvec.push_back(it2->second.get());
  //   }
  // vector<Element3D*> evec;
  // for (auto it3=emap_.begin(); it3!=emap_.end(); ++it3)
  //   {
  //     evec.push_back(it3->second.get());
  //   }

  // for (auto it2=bsplines_.begin(); it2!=bsplines_.end(); ++it2)
  //   {
  //     LRBSpline3D *currb = it2->second.get();
  //     vector<Element3D*> supp_el = currb->supportedElements();
  //     for (size_t ki=0; ki<supp_el.size(); ++ki)
  // 	{
  // 	  vector<LRBSpline3D*> supp_b = supp_el[ki]->getSupport();
  // 	  size_t kj;
  // 	  for (kj=0; kj<supp_b.size(); ++kj)
  // 	    if (supp_b[kj] == currb)
  // 	      break;
  // 	  if (kj == supp_b.size())
  // 	    std::cout << "Missing B-spline pointer in element" << std::endl;
  // 	}
  //   }

  // for (auto it3=emap_.begin(); it3!=emap_.end(); ++it3)
  //   {
  //     Element3D *currel = it3->second.get();
  //     vector<LRBSpline3D*> supp_b = currel->getSupport();
  //     for (size_t ki=0; ki<supp_b.size(); ++ki)
  // 	{
  // 	  vector<Element3D*> supp_el = supp_b[ki]->supportedElements();
  // 	  size_t kj;
  // 	  for (kj=0; kj<supp_el.size(); ++kj)
  // 	    if (supp_el[kj] == currel)
  // 	      break;
  // 	  if (kj == supp_el.size())
  // 	    std::cout << "Missing element pointer in B-spline " << std::endl;
  // 	}
  //   }

  curr_element_ = NULL;  // TESTING
}


//==============================================================================
  void LRSplineVolume::refine(const Refinement3D& ref, bool absolute)
//==============================================================================
{
    refine(ref.d, ref.kval,
           ref.start1, ref.end1,
           ref.start2, ref.end2,
           ref.multiplicity, absolute);
}

  struct refval
  {
    refval(int val, int m)
    {
      kval = val;
      mult = m;
    }
    int kval, mult;
  };

  int compare_refs(refval r1, refval r2)
  {
    return (r1.kval < r2.kval);
    // if (r1.kval < r2.kval)
    //   return -1;
    // else if (r2.kval < r1.kval)
    //   return 1;
    // else
    //   return 0;
  }

//==============================================================================
  void LRSplineVolume::refine(const std::vector<Refinement3D>& refs, bool absolute)
//==============================================================================
{
  doRefine(refs, absolute);
  
  //std::wcout << "Finally, reconstructing element map." << std::endl;
  emap_ = construct_element_map_(mesh_, bsplines_); // reconstructing the emap once at the end
}

//==============================================================================
  void LRSplineVolume::doRefine(const std::vector<Refinement3D>& refs, bool absolute)
//==============================================================================
{

#if 0//ndef NDEBUG
  {
    vector<LRBSpline3D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    puts("refine(): Remove when done debugging!");
  }
#endif

  // Collect B-splines affected by the refinements
  vector<LRBSpline3D*> bsplines_affected;
  LRSpline3DUtils::get_affected_bsplines(refs, emap_, knot_tol_, mesh_, 
  					 bsplines_affected);

  std::vector<refval> u_refs, v_refs, w_refs;
  for (size_t i = 0; i != refs.size(); ++i) {
    bool refined;
    const Refinement3D& r = refs[i];
    const auto indices = // tuple<int, int, int, int>
      LRSpline3DUtils::refine_mesh(r.d, 
				   r.kval, 
				   r.start1, 
				   r.end1, 
				   r.start2, 
				   r.end2, 
				   r.multiplicity, 
				   absolute,
				   degree(r.d), 
				   knot_tol_, 
				   mesh_, 
				   (r.d == XDIR) ?  bsplinesuni1_ : 
				   (r.d == YDIR) ? bsplinesuni2_ : bsplinesuni3_,
				   refined);

    if (!refined)
      continue;
    
    refval curr(r.kval, r.multiplicity);
    size_t kr;
    if (r.d == XDIR)
      {
	for (kr=0; kr<u_refs.size(); ++kr)
	  if (u_refs[kr].kval == curr.kval && 
	      u_refs[kr].mult == curr.mult)
	    break;
	if (kr == u_refs.size())
	  u_refs.push_back(curr);
      }
    else if (r.d == YDIR)
      {
	for (kr=0; kr<v_refs.size(); ++kr)
	  if (v_refs[kr].kval == curr.kval && 
	      v_refs[kr].mult == curr.mult)
	    break;
	if (kr == v_refs.size())
	  v_refs.push_back(curr);
      }
    else
      {
	for (kr=0; kr<w_refs.size(); ++kr)
	  if (w_refs[kr].kval == curr.kval && 
	      w_refs[kr].mult == curr.mult)
	    break;
	if (kr == w_refs.size())
	  w_refs.push_back(curr);
      }
  }

  if (u_refs.size() + v_refs.size() + w_refs.size() == 0)
    return; // Mesh rectangles already existing
  
  //std::cout << "Post refine mesh" << std::endl;
  std::sort(u_refs.begin(), u_refs.end(), compare_refs);
  std::sort(v_refs.begin(), v_refs.end(), compare_refs);
  std::sort(w_refs.begin(), w_refs.end(), compare_refs);

  #ifdef DEBUG10
  std::ofstream ofu0("uni_u0.txt");
  for (size_t kii=0; kii<bsplinesuni1_.size(); ++kii)
    {
      int nn = bsplinesuni1_[kii]->getCount();
      ofu0 << bsplinesuni1_[kii].get() << " " << nn << " ";
      bsplinesuni1_[kii]->write(ofu0);
      if (kii > 0)
	{
	  const vector<int> kvec1 = bsplinesuni1_[kii-1]->kvec();
	  const vector<int> kvec2 = bsplinesuni1_[kii]->kvec();
	  for (int ki2=0; ki2<nn; ++ki2)
	    {
	      if (kvec2[ki2] > kvec1[ki2])
		break;
	      else if (kvec2[ki2] < kvec1[ki2])
		{
		  ofu0 << "ERROR" << std::endl;
		  std::cout << "ERROR u0 " << kii << std::endl;
		}
	    }
	}
    }
#endif
  
  for (int ki=(int)u_refs.size()-1; ki>=0; --ki)
    {
      int kval = u_refs[ki].kval;
      int fixed_ix = Mesh3DUtils::last_nonlarger_knotvalue_ix(mesh_, XDIR, kval);

      int last_ix = 
	BSplineUniUtils::last_overlapping_bsplineuni(fixed_ix, bsplinesuni1_);
    
      LRSplineUtils::split_univariate(bsplinesuni1_, last_ix, fixed_ix, 
				      absolute? u_refs[ki].mult : 1);
    }

#ifdef DEBUG10
  std::ofstream ofu1("uni_u1.txt");
  for (size_t kii=0; kii<bsplinesuni1_.size(); ++kii)
    {
      int nn = bsplinesuni1_[kii]->getCount();
      ofu1 << bsplinesuni1_[kii].get() << " " << nn << " ";
      bsplinesuni1_[kii]->write(ofu1);
      if (kii > 0)
	{
	  const vector<int> kvec1 = bsplinesuni1_[kii-1]->kvec();
	  const vector<int> kvec2 = bsplinesuni1_[kii]->kvec();
	  for (int ki2=0; ki2<nn; ++ki2)
	    {
	      if (kvec2[ki2] > kvec1[ki2])
		break;
	      else if (kvec2[ki2] < kvec1[ki2])
		{
		  ofu1 << "ERROR" << std::endl;
		  std::cout << "ERROR u1 " << kii << std::endl;
		}
	    }
	}
    }
  
  std::ofstream ofv0("uni_v0.txt");
  for (size_t kii=0; kii<bsplinesuni2_.size(); ++kii)
    {
      int nn = bsplinesuni2_[kii]->getCount();
      ofv0 << bsplinesuni2_[kii].get() << " " << nn << " ";
      bsplinesuni2_[kii]->write(ofv0);
      if (kii > 0)
	{
	  const vector<int> kvec1 = bsplinesuni2_[kii-1]->kvec();
	  const vector<int> kvec2 = bsplinesuni2_[kii]->kvec();
	  for (int ki2=0; ki2<nn; ++ki2)
	    {
	      if (kvec2[ki2] > kvec1[ki2])
		break;
	      else if (kvec2[ki2] < kvec1[ki2])
		{
		  ofv0 << "ERROR" << std::endl;
		  std::cout << "ERROR v0 " << kii << std::endl;
		}
	    }
	}
     }
 #endif
  
  for (int ki=(int)v_refs.size()-1; ki>=0; --ki)
    {
      int kval = v_refs[ki].kval;
      int fixed_ix = Mesh3DUtils::last_nonlarger_knotvalue_ix(mesh_, YDIR, kval);

      int last_ix = 
	BSplineUniUtils::last_overlapping_bsplineuni(fixed_ix, bsplinesuni2_);
    
      LRSplineUtils::split_univariate(bsplinesuni2_, last_ix, fixed_ix, 
				      absolute? v_refs[ki].mult : 1); 
    }
  
#ifdef DEBUG10
  std::ofstream ofv1("uni_v1.txt");
  for (size_t kii=0; kii<bsplinesuni2_.size(); ++kii)
    {
      int nn = bsplinesuni2_[kii]->getCount();
      ofv1 << bsplinesuni2_[kii].get() << " " << nn << " ";
      bsplinesuni2_[kii]->write(ofv1);
      if (kii > 0)
	{
	  const vector<int> kvec1 = bsplinesuni2_[kii-1]->kvec();
	  const vector<int> kvec2 = bsplinesuni2_[kii]->kvec();
	  for (int ki2=0; ki2<nn; ++ki2)
	    {
	      if (kvec2[ki2] > kvec1[ki2])
		break;
	      else if (kvec2[ki2] < kvec1[ki2])
		{
		  ofv1 << "ERROR" << std::endl;
		  std::cout << "ERROR v1 " << kii << std::endl;
		}
	    }
	}
    }
  
  std::ofstream ofw0("uni_w0.txt");
  for (size_t kii=0; kii<bsplinesuni3_.size(); ++kii)
    {
     int nn = bsplinesuni3_[kii]->getCount();
     ofw0 << bsplinesuni3_[kii].get() << " " << nn << " ";
     bsplinesuni3_[kii]->write(ofw0);
      if (kii > 0)
	{
	  const vector<int> kvec1 = bsplinesuni3_[kii-1]->kvec();
	  const vector<int> kvec2 = bsplinesuni3_[kii]->kvec();
	  for (int ki2=0; ki2<nn; ++ki2)
	    {
	      if (kvec2[ki2] > kvec1[ki2])
		break;
	      else if (kvec2[ki2] < kvec1[ki2])
		{
		  ofw0 << "ERROR" << std::endl;
		  std::cout << "ERROR w0 " << kii << std::endl;
		}
	    }
	}
    }
#endif
  
  for (int ki=(int)w_refs.size()-1; ki>=0; --ki)
    {
      int kval = w_refs[ki].kval;
      int fixed_ix = Mesh3DUtils::last_nonlarger_knotvalue_ix(mesh_, ZDIR, kval);

      int last_ix = 
	BSplineUniUtils::last_overlapping_bsplineuni(fixed_ix, bsplinesuni3_);
    
      LRSplineUtils::split_univariate(bsplinesuni3_, last_ix, fixed_ix, 
				      absolute? w_refs[ki].mult : 1);
    }
#ifdef DEBUG10
  std::ofstream ofw1("uni_w1.txt");
  for (size_t kii=0; kii<bsplinesuni3_.size(); ++kii)
    {
     int nn = bsplinesuni3_[kii]->getCount();
      ofw1 << bsplinesuni3_[kii].get() << " " << nn << " ";
      bsplinesuni3_[kii]->write(ofw1);
       if (kii > 0)
	{
	  const vector<int> kvec1 = bsplinesuni3_[kii-1]->kvec();
	  const vector<int> kvec2 = bsplinesuni3_[kii]->kvec();
	  for (int ki2=0; ki2<nn; ++ki2)
	    {
	      if (kvec2[ki2] > kvec1[ki2])
		break;
	      else if (kvec2[ki2] < kvec1[ki2])
		{
		  ofw1 << "ERROR" << std::endl;
		  std::cout << "ERROR w1 " << kii << std::endl;
		}
	    }
	}
   }
#endif
  //std::cout << "Post split univariate" << std::endl;

  // //std::wcout << "Preparing for iterative splitting." << std::endl;
  // vector<unique_ptr<LRBSpline3D> > affected;
  // affected.reserve(bsplines_.size());
  // for (auto it = bsplines_.begin(); it!= bsplines_.end(); ++it)
  //   {
  //     // @@@ VSK. This is maybe the place to remove element information from the bsplines?
  //     unique_ptr<LRBSpline3D> ptr = std::move(it->second);
  //     affected.emplace_back(std::move(ptr));//b.second);
  //   }

  // @@@ VSK. In this case, we should not bother about splitting elements. They will
  // be regenerated later. Thus, the bsplines should NOT be updated with elements during
  // splitting
  // The bsplines should not have any pointers to elements. They will be set later
  //std::wcout << "Iteratively splitting." << std::endl;

  LRSpline3DUtils::iteratively_split2(bsplines_affected, mesh_, bsplines_, NULL,
				      bsplinesuni1_, bsplinesuni2_, 
				      bsplinesuni3_, false);
  // LRSpline3DUtils::iteratively_split(affected, mesh_, bsplinesuni1_, bsplinesuni2_,
  // 				     bsplinesuni3_);
  // bsplines_.clear();

  //std::wcout << "Splitting finished, now inserting resulting functions" << std::endl;
  // The bsplines are checked for duplicates and inserted in the global bspline map
  // for (auto it = affected.begin(); it != affected.end(); ++it)
  //   {
  //     LRSpline3DUtils::insert_basis_function(*it, mesh_, bsplines_);
  //   }

#if 0//ndef NDEBUG
  {
    vector<LRBSpline3D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    puts("refine(): Remove when done debugging!");
  }
#endif

  //affected.clear();

  //std::cout << "Post iteratively split" << std::endl;

  // Remove unused univariate B-splines
  for (int i=(int)bsplinesuni1_.size()-1; i>=0; --i)
    if (bsplinesuni1_[i]->getCount() <= 0)
      bsplinesuni1_.erase(bsplinesuni1_.begin()+i);

  for (int i=(int)bsplinesuni2_.size()-1; i>=0; --i)
    if (bsplinesuni2_[i]->getCount() <= 0)
      bsplinesuni2_.erase(bsplinesuni2_.begin()+i);

  for (int i=(int)bsplinesuni3_.size()-1; i>=0; --i)
    if (bsplinesuni3_[i]->getCount() <= 0)
      bsplinesuni3_.erase(bsplinesuni3_.begin()+i);

  //std::cout << "Pre construct element map" << std::endl;

  curr_element_ = NULL;  // Not valid any more
  
  //std::cout << "Post construct element map" << std::endl;
  //std::wcout << "Refinement now finished. " << std::endl;
#if 0//ndef NDEBUG
  {
    vector<LRBSpline3D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    puts("refine(): Remove when done debugging!");
  }
#endif


}
#if 0
//==============================================================================
void LRSplineVolume::zero_knot(Direction3D d, double knotval)
//==============================================================================
{
  LRSpline3DUtils::zero_knot(d, knotval, knot_tol_, mesh_,
			     (d == XDIR) ?  bsplinesuni1_ : 
			     (d == YDIR) ? bsplinesuni2_ : bsplinesuni3_);
}
#endif
//==============================================================================
  void LRSplineVolume::refineBasisFunction(int index)
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  void LRSplineVolume::refineBasisFunction(const std::vector<int> &indices)
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  void LRSplineVolume::refineElement(int index)
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  void LRSplineVolume::refineElement(const std::vector<int> &indices)
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================  
  void LRSplineVolume::setCoef(const Point& value, const LRBSpline3D* target)
//==============================================================================
{
  const auto it = bsplines_.find(generate_key(*target, mesh_));
  if (it == bsplines_.end())
    THROW("setCoef:: 'target' argument does not refer to member basis function.");

  if (value.dimension() != this->dimension())
    THROW("setCoef:: incorrect dimension of 'value' argument.");

  // if we got here, calling contract is fulfilled
  const double gamma = it->second->gamma();
  it->second->coefTimesGamma() = value * gamma;
}

//==============================================================================
  void LRSplineVolume::setCoef(const Point& value,
			       int umin_ix, int vmin_ix, int wmin_ix, int umax_ix, int vmax_ix, int wmax_ix,
			       int u_mult, int v_mult, int w_mult)
//==============================================================================
{
   MESSAGE("LRSplineVolume:: Not implemented yet.");
   throw;
}
  
//==============================================================================
  void LRSplineVolume::expandToFullTensorProduct() 
//==============================================================================
{
    Mesh3D tensor_mesh = mesh_;

    const vector<int> xmults = LRSpline3DUtils::set_uniform_meshlines(XDIR,
                                                                    tensor_mesh);
    const vector<int> ymults = LRSpline3DUtils::set_uniform_meshlines(YDIR,
                                                                    tensor_mesh);
    const vector<int> zmults = LRSpline3DUtils::set_uniform_meshlines(ZDIR,
                                                                    tensor_mesh);

    BSplineMap tensor_bsplines;

    ElementMap emap = LRSpline3DUtils::identify_elements_from_mesh(tensor_mesh);

    // splitting up basis functions
    for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) {
        LRSpline3DUtils::tensor_split(b->second,
                                      xmults,
                                      ymults,
                                      zmults,
                                      tensor_mesh,
				      bsplinesuni1_,
				      bsplinesuni2_,
				      bsplinesuni3_,
                                      tensor_bsplines);
      }

    // registering all the produced functions with the elements
    //std::wcout << "LRSplineVolume::ExpandToFullTensorProduct() - registering produced functions..." << std::endl;
    //std::wcout << "Number of basis functions: " << tensor_bsplines.size() << std::endl;
    //std::wcout << "Number of elements: "<< tensor_mesh.numDistinctKnots(XFIXED)-1 << " x " ;
    //std::wcout << tensor_mesh.numDistinctKnots(YFIXED)-1 << std::endl;

    // @@@ VSK. Use information in the LRB-splines or regenerate all elements ?
    for (auto b = tensor_bsplines.begin(); b != tensor_bsplines.end(); ++b)  {
        LRSpline3DUtils::update_elements_with_single_bspline(b->second.get(), emap,
                                                           tensor_mesh, false);
      }

    //std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - swapping and exiting." << std::endl;
    mesh_.swap(tensor_mesh);
    bsplines_.swap(tensor_bsplines);
    emap_.swap(emap);
}

//==============================================================================
void LRSplineVolume::addVolume(const LRSplineVolume& other_vol, double fac)
//==============================================================================
{
  double tol = 1.0e-12;  // Numeric noise
  int dim = dimension();

  if (numBasisFunctions() != other_vol.numBasisFunctions())
    THROW("Different number of basis functions in addVolume()");

  if (numElements() != other_vol.numElements())
    THROW("Different number of elements in addVolume()");

  // Update all basis functions
  BSplineMap::const_iterator it1 = basisFunctionsBegin();
  BSplineMap::const_iterator it2 = other_vol.basisFunctionsBegin();
  for (; it1 != basisFunctionsEnd(); ++it1, ++it2)
    {
      // @@@ Domain checking should be done. Postponed.
      
      // it is assumed that the gamma scaling factor is identical for both 
      // basis functions
      if (fabs(it1->second->gamma() - it2->second->gamma()) > tol)
	THROW("Not corresponding scaling factors in addVolume()");

      Point coef = 
	(it1->second->coefTimesGamma() + fac*it2->second->coefTimesGamma()) / 
	it1->second->gamma();

      setCoef(coef, it1->second.get());
    }   
}

//==============================================================================
SplineVolume* LRSplineVolume::asSplineVolume() 
//==============================================================================
{
  // Make full tensor product volume
  shared_ptr<LRSplineVolume> vol0;
  LRSplineVolume *vol;
  if (isFullTensorProduct()) 
    vol = this;
  else
    {
      vol0 = shared_ptr<LRSplineVolume>(clone());
      vol0->expandToFullTensorProduct();
      vol = vol0.get();
    }

  // Construct knot vectors
  const Mesh3D& mesh = vol->mesh();
  vector<double> knotsu, knotsv, knotsw;
  const double *knot;
  int ki;
  for (knot=mesh.knotsBegin(XDIR), ki=0; knot!= mesh.knotsEnd(XDIR); 
       ++knot, ++ki)
    {
      // Fetch knot multiplicity
      int mult = mesh.largestMultInLine(XDIR, ki);  // Constant for all intervals
      for (int kj=0; kj<mult; ++kj)
	knotsu.push_back(*knot);
    }

  for (knot=mesh.knotsBegin(YDIR), ki=0; knot!= mesh.knotsEnd(YDIR); 
       ++knot, ++ki)
    {
      // Fetch knot multiplicity
      int mult = mesh.largestMultInLine(YDIR, ki);  // Constant for all intervals
      for (int kj=0; kj<mult; ++kj)
	knotsv.push_back(*knot);
    }

  for (knot=mesh.knotsBegin(ZDIR), ki=0; knot!= mesh.knotsEnd(ZDIR); 
       ++knot, ++ki)
    {
      // Fetch knot multiplicity
      int mult = mesh.largestMultInLine(ZDIR, ki);  // Constant for all intervals
      for (int kj=0; kj<mult; ++kj)
	knotsw.push_back(*knot);
    }

  // Polynomial degree
  int deg_u = vol->degree(XDIR);
  int deg_v = vol->degree(YDIR);
  int deg_w = vol->degree(ZDIR);

  // Coefficients
  int num_u = (int)knotsu.size() - deg_u - 1;
  int num_v = (int)knotsv.size() - deg_v - 1;
  int num_w = (int)knotsw.size() - deg_w - 1;
  
  vector<double> coefs;
  for (LRSplineVolume::BSplineMap::const_iterator it=vol->basisFunctionsBegin();
       it != vol->basisFunctionsEnd(); ++it)
    {
      Point cf = it->second->Coef();
      coefs.insert(coefs.end(), cf.begin(), cf.end());
      if (rational_)
	{
	  double wgt = it->second->weight();
	  coefs.push_back(wgt);
	}
    }

  // Make spline volume
  SplineVolume *splvol = new SplineVolume(num_u, num_v, num_w, 
                                          deg_u+1, deg_v+1, deg_w+1,
                                          &knotsu[0], &knotsv[0], &knotsw[0],
                                          &coefs[0],
                                          dimension(), rational_);
  return splvol;
}


//==============================================================================
bool LRSplineVolume::rational() const
//==============================================================================
{
  return rational_;
}

//==============================================================================
void LRSplineVolume::to3D()
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
   void LRSplineVolume::translate(const Point& vec)
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  bool LRSplineVolume::isLeftHanded()
//==============================================================================
{
  MESSAGE("LRSplineVolume:: Not implemented yet.");
  throw;
}

//==============================================================================
  LRSplineVolume::LRSplineVolume(double knot_tol, bool rational,
				 Mesh3D& mesh, 
				 std::vector<LRBSpline3D*> b_splines,
				 int first_ixu, int first_ixv, int first_ixw)
//==============================================================================
  : knot_tol_(knot_tol), rational_(rational), mesh_(mesh), curr_element_(NULL)
{
  int left1=0, left2=0, left3=0;
  for (size_t ki=0; ki<b_splines.size(); ++ki)
  {
    BSplineUniLR *uni1 = new BSplineUniLR(*b_splines[ki]->getUnivariate(XDIR));
    BSplineUniLR *uni2 = new BSplineUniLR(*b_splines[ki]->getUnivariate(YDIR));
    BSplineUniLR *uni3 = new BSplineUniLR(*b_splines[ki]->getUnivariate(ZDIR));
    uni1->setMesh(&mesh_);  // Note that the input mesh will go out of scope. The mesh in the
    // surface is not the same
    uni1->subtractKnotIdx(first_ixu);
    uni2->setMesh(&mesh_);
    uni2->subtractKnotIdx(first_ixv);
    uni3->setMesh(&mesh_);
    uni3->subtractKnotIdx(first_ixw);

    // Check existence of univariate B-splines
    // Create new if necessary
    bool found1 = BSplineUniUtils::identify_bsplineuni(uni1, bsplinesuni1_, left1);
    if (found1)
      delete uni1;
    else
      {
	BSplineUniUtils::insert_univariate(bsplinesuni1_, uni1, left1);
      }

    bool found2 = BSplineUniUtils::identify_bsplineuni(uni2, bsplinesuni2_, left2);
    if (found2)
      delete uni2;
    else
      {
	BSplineUniUtils::insert_univariate(bsplinesuni2_, uni2, left2);
      }

    bool found3 = BSplineUniUtils::identify_bsplineuni(uni3, bsplinesuni3_, left3);
    if (found3)
      delete uni3;
    else
      {
	BSplineUniUtils::insert_univariate(bsplinesuni3_, uni3, left3);
      }

    LRBSpline3D* bb = new LRBSpline3D(b_splines[ki]->coefTimesGamma(),
				      b_splines[ki]->weight(), 
				      bsplinesuni1_[left1].get(),
				      bsplinesuni2_[left2].get(),
				      bsplinesuni3_[left3].get(),
				      b_splines[ki]->gamma(), b_splines[ki]->rational());

    LRSplineVolume::BSKey bs_key = generate_key(*bb);
    //b_splines[ki]->setMesh(&mesh_);
    bsplines_.insert(std::pair<LRSplineVolume::BSKey,
		     unique_ptr<LRBSpline3D> >(bs_key, unique_ptr<LRBSpline3D>(bb)));
  }

  emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
  LRSplineVolume::ElementMap LRSplineVolume::construct_element_map_(const Mesh3D& m, const BSplineMap& bmap)
//==============================================================================
{
  ElementMap emap = LRSpline3DUtils::identify_elements_from_mesh(m);

  for (auto b_it = bmap.begin(); b_it != bmap.end(); ++b_it) 
    {
      LRBSpline3D* tmp = b_it->second.get();

      // First remove all old elements in support
      tmp->removeSupportedElements();
      
      LRSpline3DUtils::update_elements_with_single_bspline(tmp, emap, 
							   m, false);
    }
  return emap;
}

//==============================================================================
std::vector<Go::LRBSpline3D *> LRSplineVolume::collect_basis(int from_u, int to_u,
				int from_v, int to_v,
				int from_w, int to_w) const
//==============================================================================
{
  vector<LRBSpline3D*> b_splines;
  for (auto it=bsplines_.begin(); it!=bsplines_.end(); ++it)
    {
      if (interval_overlap(it->second->suppMin(XDIR), 
			   it->second->suppMax(XDIR), 
			   from_u, to_u) &&
	  interval_overlap(it->second->suppMin(YDIR), 
			   it->second->suppMax(YDIR), 
			   from_v, to_v) &&
	  interval_overlap(it->second->suppMin(ZDIR), 
			   it->second->suppMax(ZDIR), 
			   from_w, to_w))
	{
	  b_splines.push_back(it->second.get());
	}
    }
  return b_splines;
}



} // end namespace Go
