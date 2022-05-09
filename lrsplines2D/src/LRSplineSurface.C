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

#include "GoTools/lrsplines2D/LRSplineSurface.h"

#include <iomanip>
#include <stdexcept>
#include <iostream> // @@ debug
#include <fstream>
#include <iterator> // @@ debug - remove
#include <string.h>
//#include <chrono>   // @@ debug
#include <set>
#include <tuple>
#include "GoTools/utils/checks.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
#include "GoTools/lrsplines2D/Mesh2DUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"
#include "GoTools/utils/StreamUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h" // @@ only for debug
#include "GoTools/geometry/Utils.h"

//#define DEBUG

using std::vector;
using std::istream;
using std::ostream;
using std::get;
using std::pair;
using std::make_pair;
using std::for_each;

using std::unique_ptr;

//==============================================================================
namespace Go
//==============================================================================
{

//==============================================================================
LRSplineSurface::ElementMap 
LRSplineSurface::construct_element_map_(const Mesh2D& m, const BSplineMap& bmap)
//==============================================================================
{
  ElementMap emap = LRSplineUtils::identify_elements_from_mesh(m);

  for (auto b_it = bmap.begin(); b_it != bmap.end(); ++b_it) 
    {
      LRBSpline2D* tmp = b_it->second.get();

      // First remove all old elements in support
      tmp->removeSupportedElements();
      
      LRSplineUtils::update_elements_with_single_bspline(tmp, emap, 
							 m, false);
    }

  return emap;
};

//==============================================================================
LRSplineSurface::LRSplineSurface(const SplineSurface* const surf, 
				 double knot_tol)
//==============================================================================
  : knot_tol_(knot_tol), rational_(surf->rational()), curr_element_(NULL),
  mesh_(surf->basis_u().begin(), surf->basis_u().end(),
	surf->basis_v().begin(), surf->basis_v().end())
{
  std::vector<int> knot_ixs_u = init_knot_indices(mesh_, XFIXED);
  std::vector<int> knot_ixs_v = init_knot_indices(mesh_, YFIXED);

  int deg_u = surf->order_u() - 1;
  int deg_v = surf->order_v() - 1;
  int coefs_u = surf->numCoefs_u();
  int coefs_v = surf->numCoefs_v();
  std::vector<double>::const_iterator rcoefs = surf->rcoefs_begin();
  std::vector<double>::const_iterator coefs = surf->coefs_begin();
  int dim = surf->dimension();
  int kdim = (rational_) ? dim + 1 : dim;

  // Store uni-variate B-splines
  bsplinesuni1_.resize(coefs_u);
  bsplinesuni2_.resize(coefs_v);
  for (int u_ix = 0; u_ix != coefs_u; ++u_ix) {
    bsplinesuni1_[u_ix] = std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(1, deg_u,
											 knot_ixs_u.begin() + u_ix,
											 &mesh_)));
  }

  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    bsplinesuni2_[v_ix] = std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(2, deg_v,
											 knot_ixs_v.begin() + v_ix,
											 &mesh_)));
  }

  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    for (int u_ix = 0; u_ix != coefs_u; ++u_ix, coefs+=dim, rcoefs += kdim) {
      double rat = (rational_) ? rcoefs[dim] : 1.0;
      unique_ptr<LRBSpline2D> b(new LRBSpline2D(Point(coefs, coefs + dim),
						rat,
						bsplinesuni1_[u_ix].get(),
						bsplinesuni2_[v_ix].get(),
						1.0, rational_));
//      bsplines_[generate_key(*b, mesh_)] = b;
      LRSplineSurface::BSKey bs_key = generate_key(*b, mesh_);
      bsplines_.insert(std::pair<LRSplineSurface::BSKey, unique_ptr<LRBSpline2D> >(bs_key, std::move(b)));
    }
  }
  emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
LRSplineSurface::LRSplineSurface(double knot_tol, bool rational,
				 Mesh2D& mesh, 
				 vector<LRBSpline2D*>& b_splines,
				 int first_ixu, int first_ixv)
//==============================================================================
  : knot_tol_(knot_tol), rational_(rational), mesh_(mesh), curr_element_(NULL)
{
  int left1=0, left2=0;
  for (size_t ki=0; ki<b_splines.size(); ++ki)
  {
    BSplineUniLR *uni1 = new BSplineUniLR(*b_splines[ki]->getUnivariate(XFIXED));
    BSplineUniLR *uni2 = new BSplineUniLR(*b_splines[ki]->getUnivariate(YFIXED));
    uni1->setMesh(&mesh_);  // Note that the input mesh will go out of scope. The mesh in the
    // surface is not the same
    uni1->subtractKnotIdx(first_ixu);
    uni2->setMesh(&mesh_);
    uni2->subtractKnotIdx(first_ixv);

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

    LRBSpline2D* bb = new LRBSpline2D(b_splines[ki]->coefTimesGamma(),
				      b_splines[ki]->weight(), 
				      bsplinesuni1_[left1].get(),
				      bsplinesuni2_[left2].get(),
				      b_splines[ki]->gamma(), b_splines[ki]->rational());

    LRSplineSurface::BSKey bs_key = generate_key(*bb);
    //b_splines[ki]->setMesh(&mesh_);
    bsplines_.insert(std::pair<LRSplineSurface::BSKey, unique_ptr<LRBSpline2D> >(bs_key, unique_ptr<LRBSpline2D>(bb)));
  }

  emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
LRSplineSurface::LRSplineSurface(const LRSplineSurface& rhs) 
//==============================================================================
  : knot_tol_(rhs.knot_tol_), rational_(rhs.rational_),  curr_element_(NULL),
    mesh_(rhs.mesh_)
{
  // Clone univariate B-splines
  vector<std::unique_ptr<BSplineUniLR> >::const_iterator curruni1 = rhs.bsplinesuni1_.begin();
  vector<std::unique_ptr<BSplineUniLR> >::const_iterator enduni1 = rhs.bsplinesuni1_.end();
  bsplinesuni1_.resize(rhs.bsplinesuni1_.size());
  size_t ki;
  for (ki=0; curruni1 != enduni1; ++ki, ++curruni1)
    {
      unique_ptr<BSplineUniLR> b(new BSplineUniLR(*(*curruni1)));

      // Update mesh pointer
      b->setMesh(&mesh_);
      bsplinesuni1_[ki] = std::move(b);
    }
  
  vector<std::unique_ptr<BSplineUniLR> >::const_iterator curruni2 = rhs.bsplinesuni2_.begin();
  vector<std::unique_ptr<BSplineUniLR> >::const_iterator enduni2 = rhs.bsplinesuni2_.end();
  bsplinesuni2_.resize(rhs.bsplinesuni2_.size());
  for (ki=0; curruni2 != enduni2; ++ki, ++curruni2)
    {
      unique_ptr<BSplineUniLR> b(new BSplineUniLR(*(*curruni2)));

      // Update mesh pointer
      b->setMesh(&mesh_);
      bsplinesuni2_[ki] = std::move(b);
    }

   // Clone LR B-splines
  int left1 = 0, left2 = 0;
  BSplineMap::const_iterator curr = rhs.basisFunctionsBegin();
  BSplineMap::const_iterator end = rhs.basisFunctionsEnd();
  for (; curr != end; ++curr)
    {
      unique_ptr<LRBSpline2D> b(new LRBSpline2D(*curr->second));

     const BSplineUniLR *uni1 = b->getUnivariate(XFIXED);
      bool found1 = BSplineUniUtils::identify_bsplineuni(uni1, bsplinesuni1_, left1);
      if (!found1)
	THROW("Univariate B-spline not found");
      b->setUnivariate(XFIXED, bsplinesuni1_[left1].get());
      //bsplinesuni1_[left1]->incrCount();

      const BSplineUniLR *uni2 = b->getUnivariate(YFIXED);
      bool found2 = BSplineUniUtils::identify_bsplineuni(uni2, bsplinesuni2_, left2);
      if (!found2)
	THROW("Univariate B-spline not found");
      b->setUnivariate(YFIXED, bsplinesuni2_[left2].get());
      //bsplinesuni2_[left2]->incrCount();
      
      // bsplines_[generate_key(*b, mesh_)] = b;
      LRSplineSurface::BSKey bs_key = generate_key(*b, mesh_);
      bsplines_.insert(std::pair<LRSplineSurface::BSKey, unique_ptr<LRBSpline2D> >(bs_key, std::move(b)));
    }

  // The ElementMap has to be generated and cannot be copied directly, since it
  // contains raw pointers.  
  emap_ = construct_element_map_(mesh_, bsplines_);
}

//===========================================================================
LRSplineSurface::~LRSplineSurface()

//===========================================================================
{
}

#if 1
//==============================================================================
const LRSplineSurface& LRSplineSurface::operator= (const LRSplineSurface& other)
//==============================================================================
{
  LRSplineSurface lr_spline_sf(other);
  this->swap(lr_spline_sf);
  
  return *this;
}
#endif

//==============================================================================
void LRSplineSurface::swap(LRSplineSurface& rhs)
//==============================================================================
{
  std::swap(knot_tol_,    rhs.knot_tol_);
  std::swap(rational_,    rhs.rational_);
  std::swap(mesh_    ,    rhs.mesh_);
  std::swap(bsplinesuni1_,    rhs.bsplinesuni1_);
  std::swap(bsplinesuni2_,    rhs.bsplinesuni2_);
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
void  LRSplineSurface::read(istream& is)
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

  vector<unique_ptr<LRBSpline2D> > dummy_vec(num_bfuns);


#if 0//ndef NDEBUG
  vector<LRBSpline2D*> debug_vec(num_bfuns);
#endif
  int left1 = 0, left2 = 0;
  for (int i = 0; i != num_bfuns; ++i) {
    unique_ptr<LRBSpline2D> b(new LRBSpline2D());
//    LRBSpline2D* b = new LRBSpline2D();
    //object_from_stream(is, *b);
    b->read(is, bsplinesuni1_, left1, bsplinesuni2_, left2);
    // We set the global mesh in the b basis function.
    b->setMesh(&mesh_);
#if 0
#if 1//ndef NDEBUG
    debug_vec[i] = b;
#endif
#else
#if 0
    /*tmp.*/bsplines_[generate_key(*b, mesh_)] = b;
#else
    BSKey key = generate_key(*b, mesh_);
    /*tmp.*/bsplines_.insert(std::make_pair(key, std::move(b)));
#endif
#endif
  }

  // Reconstructing element map
  emap_ = construct_element_map_(mesh_, bsplines_);

  rational_ = rational_;

  curr_element_ = NULL;

#if 0//ndef NDEBUG
  while (true) // @@sbr201305 Checking memory consumption.
    ;
#endif

#ifdef DEBUG
  {
    vector<LRBSpline2D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    //puts("Remove when done debugging!");
    vector<Element2D*> elems;
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
void LRSplineSurface::write(ostream& os) const
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
  int ki=0;
  for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b, ++ki) 
    {
      b->second->write(os);
      // object_to_stream(os, *(b->second));
      // object_to_stream(os, '\n');
    }
  os << std::endl;

  // NB: 'emap_' is not saved (contains raw pointers to other data).  
  // Instead, it will be regenerated the LRSplineSurface is read().
    os.precision(prev);   // Reset precision to it's previous value
}

//==============================================================================
SplineSurface* LRSplineSurface::asSplineSurface() 
//==============================================================================
{
  // Make full tensor product surface
  shared_ptr<LRSplineSurface> surf0;
  LRSplineSurface *surf;
  if (isFullTensorProduct())
    surf = this;
  else
    {
      surf0 = shared_ptr<LRSplineSurface>(clone());
      surf0->expandToFullTensorProduct();
      surf = surf0.get();
    }

  // Construct knot vectors
  const Mesh2D& mesh = surf->mesh();
  vector<double> knotsu, knotsv;
  const double *knot;
  int ki;
  for (knot=mesh.knotsBegin(XFIXED), ki=0; knot!= mesh.knotsEnd(XFIXED); 
       ++knot, ++ki)
    {
      // Fetch knot multiplicity
      int mult = mesh.largestMultInLine(XFIXED, ki);  // Constant for all intervals
      for (int kj=0; kj<mult; ++kj)
	knotsu.push_back(*knot);
    }

  for (knot=mesh.knotsBegin(YFIXED), ki=0; knot!= mesh.knotsEnd(YFIXED); 
       ++knot, ++ki)
    {
      // Fetch knot multiplicity
      int mult = mesh.largestMultInLine(YFIXED, ki);  // Constant for all intervals
      for (int kj=0; kj<mult; ++kj)
	knotsv.push_back(*knot);
    }

  // Polynomial degree
  int deg_u = surf->degree(XFIXED);
  int deg_v = surf->degree(YFIXED);

  // Coefficients
  int num_u = (int)knotsu.size() - deg_u - 1;
  int num_v = (int)knotsv.size() - deg_v - 1;
  
  vector<double> coefs;
  for (LRSplineSurface::BSplineMap::const_iterator it=surf->basisFunctionsBegin();
       it != surf->basisFunctionsEnd(); ++it)
    {
      Point cf = it->second->Coef();
      coefs.insert(coefs.end(), cf.begin(), cf.end());
      if (rational_)
	{
	  double wgt = it->second->weight();
	  coefs.push_back(wgt);
	}
    }

  // Make spline surface
  SplineSurface *splsf = new SplineSurface(num_u, num_v, deg_u+1, deg_v+1,
					   &knotsu[0], &knotsv[0], &coefs[0],
					   dimension(), rational_);
  return splsf;
 }

#if 0
//==============================================================================
void LRSplineSurface::computeBasis (double param_u, double param_v, BasisPtsSf     & result, int iEl ) const
//==============================================================================
{
  MESSAGE("LRSplineSurface::computeBasis() not implemented yet");
}

//==============================================================================
void LRSplineSurface::computeBasis (double param_u, double param_v, BasisDerivsSf  & result, int iEl ) const
//==============================================================================
{
  MESSAGE("LRSplineSurface::computeBasis() not implemented yet");
}

//==============================================================================
void LRSplineSurface::computeBasis (double param_u, double param_v, BasisDerivsSf2 & result, int iEl ) const
//==============================================================================
{
  MESSAGE("LRSplineSurface::computeBasis() not implemented yet");
}

//==============================================================================
void LRSplineSurface::computeBasis(double param_u,
				   double param_v,
				   std::vector<std::vector<double> >& result,
				   int derivs,
				   int iEl) const
//==============================================================================
{
  MESSAGE("LRSplineSurface::computeBasis() not implemented yet");
}

//==============================================================================
int LRSplineSurface::getElementContaining(double u, double v) const
//==============================================================================
{
  MESSAGE("LRSplineSurface::getElementContaining() not implemented yet");

  return -1;
}
#endif

//==============================================================================
//const LRSplineSurface::ElementMap::value_type& 
Element2D*
LRSplineSurface::coveringElement(double u, double v) const
//==============================================================================
{
  int ucorner, vcorner;
  if (! Mesh2DUtils::identify_patch_lower_left(mesh_, u, v, ucorner, vcorner) ) 
  {
#ifndef NDEBUG
      std::cout << "u: " << u << ", v: " << v << std::endl;
#endif
    THROW("Parameter outside domain in LRSplineSurface::basisFunctionsWithSupportAt()");
  }

  // //#ifndef 0 //NDEBUG
  // //  {
  //   vector<LRBSpline2D*> bas_funcs;
  //   for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
  //     {
  // 	bas_funcs.push_back((*iter).second.get());
  //     }
  //   //    puts("Remove when done debugging!");
  //   vector<Element2D*> elems;
  //   vector<ElemKey> elem_keys;
  //   for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
  //   {
  //     elems.push_back(((*iter).second.get()));
  //     elem_keys.push_back(iter->first);
  //   }
  //   //    puts("Remove when done debugging!");
  //   //  }
  //   //#endif

  const LRSplineSurface::ElemKey key = 
    {mesh_.knotsBegin(XFIXED)[ucorner], mesh_.knotsBegin(YFIXED)[vcorner]};
  const auto el = emap_.find(key);
  assert(el != emap_.end());
//  return *el;
  return el->second.get();
}


//==============================================================================
 void LRSplineSurface::constructElementMesh(vector<Element2D*>& elements) const
//==============================================================================
{
  // Get all knot values in the u-direction
  const double* const uknots = mesh_.knotsBegin(XFIXED);
  int nmb_knots_u = mesh_.numDistinctKnots(XFIXED);

  // Get all knot values in the v-direction
  const double* const vknots = mesh_.knotsBegin(YFIXED);
  int nmb_knots_v = mesh_.numDistinctKnots(YFIXED);

  // Construct mesh of element pointers
  int ki, kj, kr;
  elements.resize((nmb_knots_u-1)*(nmb_knots_v-1), NULL);
  int num = numElements();
  ElementMap::const_iterator it;
  for (it=elementsBegin(), kr=0, kj=0, ki=0; kr<num; ++it, ++kr)
    {
      double umin = it->second->umin();
      double umax = it->second->umax();
      double vmin = it->second->vmin();
      double vmax = it->second->vmax();

      for (; kj<nmb_knots_v && vmin > vknots[kj]; ++kj);
      for (ki=0; ki<nmb_knots_u && umin > uknots[ki]; ++ki);
      for (int kh1=kj; kh1<nmb_knots_v-1 && vmax >= vknots[kh1+1]; ++kh1)
	for (int kh2=ki; kh2<nmb_knots_u-1 && umax >= uknots[kh2+1]; ++kh2)
	  {
	    elements[kh1*(nmb_knots_u-1)+kh2] = it->second.get();
	  }
      
     }
}

//==============================================================================
vector<LRBSpline2D*> LRSplineSurface::basisFunctionsWithSupportAt(double u, double v) const
//==============================================================================
{
  vector<LRBSpline2D*> support_functions;
  auto elem = coveringElement(u, v);
  vector<LRBSpline2D*>::const_iterator first = elem->supportBegin();
  vector<LRBSpline2D*>::const_iterator last = elem->supportEnd();
  int ki=0;
  for (; first != last; ++first, ++ki)
    {
      support_functions.push_back(*first);
    }
  
  return support_functions;
}

// =============================================================================
LRSplineSurface::BSplineMap::iterator
LRSplineSurface::bsplineFromDomain(double start_u, double start_v, 
				   double end_u, double end_v, 
				   int startmult_u, int startmult_v,
				   int endmult_u, int endmult_v)
// =============================================================================
{
  BSKey key = {start_u, start_v, end_u, end_v, startmult_u, startmult_v, 
	       endmult_u, endmult_v};
      
  // Fetch the associated LR B-spline
  const auto bm = bsplines_.find(key);
  if (bm == bsplines_.end())
    THROW("bsplineFromDomain: There is no such basis function.");
  return bm;
}

// =============================================================================
vector<LRBSpline2D*>
LRSplineSurface::getBoundaryBsplines(Direction2D d, bool atstart)
// =============================================================================
{
  vector<LRBSpline2D*> bsplines;

  // Traverse all B-splines and check whether they have maximum multiplicity along
  // the given edge
  for (BSplineMap::iterator it=basisFunctionsBeginNonconst(); 
       it != basisFunctionsEndNonconst(); ++it)
    {
      int deg = it->second->degree(d);
      int mult = (d == XFIXED) ? it->second->endmult_u(atstart) :
	it->second->endmult_v(atstart);
      if (mult == deg+1)
	bsplines.push_back(it->second.get());
    }
  return bsplines;
 }

//==============================================================================
bool LRSplineSurface::isFullTensorProduct() const
//==============================================================================
{
  return (LRSplineUtils::all_meshlines_uniform(XFIXED, mesh_) &&  
	  LRSplineUtils::all_meshlines_uniform(YFIXED, mesh_));
}

//==============================================================================
  void LRSplineSurface::refine(const Refinement2D& ref, 
			     bool absolute)
//==============================================================================
{
  refine(ref.d, ref.kval, ref.start, ref.end, ref.multiplicity, absolute);
}

//==============================================================================
void LRSplineSurface::refine(Direction2D d, double fixed_val, double start, 
			     double end, int mult, bool absolute)
//==============================================================================
{
  #ifdef DEBUG
  // std::ofstream of("mesh0.eps");
  // writePostscriptMesh(*this, of);

  std::ofstream ofbe("basis_elems.txt");
  vector<LRBSpline2D*> bas_funcs;
  int kv = 0;
  for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
    {
      bas_funcs.push_back((*iter).second.get());
      ofbe << "Bspline " << kv << ": " << (*iter).second.get() << " ";
      ofbe << (*iter).second->umin() << " " << (*iter).second->umax() << " ";
      ofbe << (*iter).second->vmin() << " " <<  (*iter).second->vmax() << std::endl;
      kv++;

      if (false)
	{
	  checkSupport((*iter).second.get());
	}
    }
  ofbe << std::endl;
  vector<Element2D*> elems;
  vector<Element2D*> elems_affected;
  kv = 0;
  for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
    {
      Element2D* el = (*iter).second.get();
      elems.push_back(el);
      ofbe << "Element " << kv <<": " << el << " ";
      ofbe << el->umin() << " " << el->umax() << " " << el->vmin();
      ofbe << " " << el->vmax() << std::endl;
      kv++;
    }
  ofbe << std::endl;
  //puts("Remove when done debugging!");
  // printf("\n");
  int stop_break = 1;
  #endif

  // Make a copy of the initial mesh
  Mesh2D mesh2 = mesh_;

  bool refined;
  const auto indices = // tuple<int, int, int, int>
    LRSplineUtils::refine_mesh(d, fixed_val, start, end, mult, 
			       absolute, degree(d), knot_tol_, mesh_,
			       (d == XFIXED) ?  bsplinesuni1_ : bsplinesuni2_,
			       refined);

#ifdef DEBUG
  std::ofstream of2("mesh1.eps");
  writePostscriptMesh(*this, of2);
#endif

  // insert newly created elements to emap (unless refinement was on border, in which case no new element
  // could possibly be created
  const int prev_ix = get<0>(indices);
  const int fixed_ix = get<1>(indices); // Index of fixed_val in global knot vector.
  const int start_ix = get<2>(indices); // Index of start (of segment to insert) in global knot vector.
  const int end_ix   = get<3>(indices); // Index of end (of segment to insert) in global knot vector.

  // Collect pointers to affected bsplines
  std::set<LRBSpline2D*> all_bsplines;
  double domain[4];  // Covers elements affected by the split
  domain[0] = domain[2] = std::numeric_limits<double>::max();
  domain[1] = domain[3] = std::numeric_limits<double>::lowest();
  for (int i = start_ix; i != end_ix; ++i) {
    // Check if the specified element exists in 'emap'
    int u_ix = (d == XFIXED) ? prev_ix : i;
    int v_ix = (d == YFIXED) ? prev_ix : i;
    ElementMap::key_type key = {mesh_.kval(XFIXED, u_ix),
				mesh_.kval(YFIXED, v_ix)};
    auto it = emap_.find(key);
    if (it == emap_.end())
      {
	// Element not found. The assumed start index of the element
	// is not correct. Recompute
	int u_ix2 = u_ix;
	int v_ix2 = v_ix;
	if (d == XFIXED)
	  u_ix2 = 
	    Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh2, XFIXED,
								   u_ix, v_ix);
	else
	  v_ix2 = 
	    Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh2, YFIXED,
								   v_ix, u_ix);
	u_ix = u_ix2;
	v_ix = v_ix2;

#if 0
	key = {mesh2.kval(XFIXED, u_ix), mesh2.kval(YFIXED, v_ix)};
#else
	key.u_min = mesh2.kval(XFIXED, u_ix);
	key.v_min = mesh2.kval(YFIXED, v_ix);
#endif
	it = emap_.find(key);
#ifdef DEBUG
	if (it == emap_.end())
	  int stop_break = 1;
	//std::cout << "LRSplineSurface::refine : Element not found" << std::endl;
#endif
      }

    if (it != emap_.end())
      {
	// The element exists. Collect bsplines
	Element2D* curr_el = (*it).second.get();
	all_bsplines.insert(curr_el->supportBegin(), curr_el->supportEnd());
	domain[0] = std::min(domain[0], curr_el->umin());
	domain[1] = std::max(domain[1], curr_el->umax());
	domain[2] = std::min(domain[2], curr_el->vmin());
	domain[3] = std::max(domain[3], curr_el->vmax());
#ifdef DEBUG
	elems_affected.push_back(curr_el);
#endif
	// if (d == XFIXED)
	//   {
	//     double u1 = curr_el->umin();
	//     double u2 = curr_el->umax();
	//     if ((u2 - fixed_val > 0.55*(u2-u1) || fixed_val - u1 > 0.55*(u2-u1)) &&
	// 	u2-fixed_val > 0.0001 && fixed_val-u1 > 0.0001)
	//       std::cout << "Unbalanced element, 1. par: " << fixed_val << " in [" << u1 << "," << u2 << "]" << std::endl;
	//   }
	// else
	//   {
	//     double v1 = curr_el->vmin();
	//     double v2 = curr_el->vmax();
	//     if ((v2 - fixed_val > 0.55*(v2-v1) || fixed_val - v1 > 0.55*(v2-v1)) &&
	// 	v2-fixed_val > 0.0001 && fixed_val-v1 > 0.0001)
	//       std::cout << "Unbalanced element, 2. par: " << fixed_val << " in [" << v1 << "," << v2 << "]" << std::endl;
	//   }
      }
  }
  vector<LRBSpline2D*> bsplines_affected(all_bsplines.begin(), all_bsplines.end());

  #ifdef DEBUG
    bas_funcs.clear();
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	LRBSpline2D *tmpb = (*iter).second.get();
	bas_funcs.push_back(tmpb);
	for (auto eit=tmpb->supportedElementBegin();
	     eit != tmpb->supportedElementEnd(); ++eit)
	  {
	    if (!(*eit)->hasSupportFunction(tmpb))
	      std::cout << "Element " << (*eit) << " misses Bspline " << tmpb << std::endl;
	    double b_umin = tmpb->umin();
	    double b_umax = tmpb->umax();
	    double b_vmin = tmpb->vmin();
	    double b_vmax = tmpb->vmax();
	    if ((*eit)->umax() <= b_umin || (*eit)->umin() >= b_umax ||
		(*eit)->vmax() <= b_vmin || (*eit)->vmin() >= b_vmax)
	      std::cout << "Element " << (*eit) << " not in Bspline support " << tmpb << std::endl;
	  }
      }
    elems.clear();
    for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
      {
	Element2D *tmpe = (*iter).second.get();
	elems.push_back(tmpe);
	if (tmpe->nmbBasisFunctions() == 0)
	  std::cout << "Element with no support functions! " << tmpe << std::endl;
	for (auto bit=tmpe->supportBegin();
	     bit != tmpe->supportEnd(); ++bit)
	  {
	    if (!(*bit)->hasSupportedElement(tmpe))
	      std::cout << "Bspline " << (*bit) << " misses Element " << tmpe << std::endl;
	  }
       }
    //puts("Remove when done debugging!");

    // Check if found bsplines exists in the map
    for (size_t kr1=0; kr1<bsplines_affected.size(); ++kr1)
      {
	auto key = generate_key(*bsplines_affected[kr1], mesh_);
	const auto it = bsplines_.find(key);
	if (it == bsplines_.end())
	  std::cout << "Bspline not in map: " << bsplines_affected[kr1] << std::endl;
      }
    int deb = 0;
    #endif

    if (bsplines_affected.size() == 0)
      return;  // No B-splines will be split, no knot insertion is possible

    // Identify range of univariate B-splines
    // int iu1, iu2, iv1, iv2;
    // bool found1 = bsplineuni_range(bsplinesuni1_, 
    // 				   (d==XFIXED) ? fixed_ix : start_ix,
    // 				   (d==XFIXED) ? fixed_ix : end_ix, iu1, iu2);
    // bool found2 = bsplineuni_range(bsplinesuni2_, 
    // 				   (d==YFIXED) ? fixed_ix : start_ix,
    // 				   (d==YFIXED) ? fixed_ix : end_ix, iv1, iv2);

    // Split univariate to prepare for bivariate split
    int last_ix = 
      BSplineUniUtils::last_overlapping_bsplineuni(fixed_ix,
						   (d == XFIXED) ? bsplinesuni1_ : bsplinesuni2_);

    if (d == XFIXED)
      {
	// if (iu2 < (int)bsplinesuni1_.size()-1)
	//   ++iu2;
	// LRSplineUtils::split_univariate(bsplinesuni1_, iu1, iu2, fixed_ix);
	LRSplineUtils::split_univariate(bsplinesuni1_, last_ix, fixed_ix, 
					(absolute) ? mult : 1);
      }
    else
      {
	// if (iv2 < (int)bsplinesuni2_.size()-1)
	//   ++iv2;
	// LRSplineUtils::split_univariate(bsplinesuni2_, iv1, iv2, fixed_ix);
	LRSplineUtils::split_univariate(bsplinesuni2_, last_ix, fixed_ix, 
					(absolute) ? mult : 1);
       }

    // for (size_t ki=1; ki<bsplinesuni1_.size(); ++ki)
    //   {
    // 	int comp = ((*bsplinesuni1_[ki-1]) < (*bsplinesuni1_[ki]));
    // 	if (comp == 0)
    // 	  std::cout << "1. Equality of univiariate B-splines" << std::endl;
    // 	else if (comp > 0)
    // 	  std::cout << "1. Error in sequence of univariate B-splines" << std::endl;
    //   }
	  
  // Cannot remove the bsplines from the global array at this stage since we operate
  // with pointers to it. When a bspline is split, the origin is removed from the
  // array after all pointers are updated and the the bspline is allowed to die.
  // Iteratively split affected LRBSpline2Ds
  // @@@ VSK. Will pointers to other entities in bsplines_ which are not
  // affected remain valid after removing and adding elements? If not, this
  // combination of objects and pointers will not work.
    LRSplineUtils::iteratively_split2(bsplines_affected, mesh_, 
				      bsplines_, domain, 
				      bsplinesuni1_, bsplinesuni2_, true);
				      // bsplinesuni1_, iu1, iu2,
				      // bsplinesuni2_, iv1, iv2);

    // for (size_t ki=1; ki<bsplinesuni1_.size(); ++ki)
    //   {
    // 	int comp = ((*bsplinesuni1_[ki-1]) < (*bsplinesuni1_[ki]));
    // 	if (comp == 0)
    // 	  std::cout << "2. Equality of univiariate B-splines" << std::endl;
    // 	else if (comp > 0)
    // 	  std::cout << "2. Error in sequence of univariate B-splines" << std::endl;
    //   }

    // Remove unused univariate B-splines
    // if (d == XFIXED)
    //   {
	//for (int i=iu2; i>=iu1; --i)
    for (int i=(int)bsplinesuni1_.size()-1 /*last_ix*/; i>=0; --i)
	  if (bsplinesuni1_[i]->getCount() <= 0)
	    bsplinesuni1_.erase(bsplinesuni1_.begin()+i);
    //   }
    // else
    //   {
	//for (int i=iv2; i>=iv1; --i)
    for (int i=(int)bsplinesuni2_.size()-1 /*last_ix*/; i>=0; --i)
	  if (bsplinesuni2_[i]->getCount() <= 0)
	    bsplinesuni2_.erase(bsplinesuni2_.begin()+i);
      // }

    // std::ofstream ofuni("uni1.g2");
    // for (size_t ki=0; ki<bsplinesuni1_.size(); ++ki)
    //   {
    // 	for (size_t kj=0; kj<bsplinesuni1_[ki]->kvec().size(); ++kj)
    // 	  ofuni << bsplinesuni1_[ki]->kvec()[kj] << " ";
    // 	ofuni << ", count: " << bsplinesuni1_[ki]->getCount() << std::endl;
    // 	if (bsplinesuni1_[ki]->getCount() == 0)
    // 	  std::cout << "Count = 0, fixed_ix = " << fixed_ix << std::endl;
    //   }
    // for (size_t ki=1; ki<bsplinesuni1_.size(); ++ki)
    //   {
    // 	int comp = ((*bsplinesuni1_[ki-1]) < (*bsplinesuni1_[ki]));
    // 	if (comp == 0)
    // 	  std::cout << "3. Equality of univiariate B-splines" << std::endl;
    // 	else if (comp > 0)
    // 	  std::cout << "3. Error in sequence of univariate B-splines" << std::endl;
    //   }

#ifdef DEBUG
    bas_funcs.clear();
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    elems.clear();
    for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
      {
	  elems.push_back(((*iter).second.get()));
	  if ((*iter).second->nmbBasisFunctions() == 0)
	    std::cout << "Element with no support functions! " << (*iter).second.get() << std::endl;
      }
    //puts("Remove when done debugging!");
    stop_break = 1;
#endif

#ifdef DEBUG
    //std::cout << "Num elements prior: " << numElements() << std::endl;
#endif
  if (fixed_ix > 0 && fixed_ix != mesh_.numDistinctKnots(d)-1) {
    for (int i = start_ix; i != end_ix; ++i) {
      if (mesh_.nu(flip(d), i, fixed_ix, fixed_ix+1) > 0) {
	// this is the lower-left corner of an element bordering our refinement.  
	// Check if it already exists in 'emap', if not, insert it.
	// Do also modify the current element and update bspline pointers
	// in the elements
	int u_ix2 = (d == XFIXED) ? prev_ix : i;
	int v_ix2 = (d == XFIXED) ? i : prev_ix;
	ElementMap::key_type key2 = {mesh_.kval(XFIXED, u_ix2),
				    mesh_.kval(YFIXED, v_ix2)};
	auto it2 = emap_.find(key2);

	if (it2 == emap_.end())
	  {
	    // Element not found. The assumed start index of the element
	    // is not correct. Recompute
	    int u_ix3 = u_ix2;
	    int v_ix3 = v_ix2;
	    if (d == XFIXED)
	      u_ix3 = 
		Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh2, XFIXED,
								       u_ix2, v_ix2);
	    else
	      v_ix3 = 
		Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh2, YFIXED,
								       v_ix2, u_ix2);

	    u_ix2 = u_ix3;
	    v_ix2 = v_ix3;

#if 0
	    key2 = {mesh2.kval(XFIXED, u_ix2), mesh2.kval(YFIXED, v_ix2)};
#else
		key2.u_min = mesh2.kval(XFIXED, u_ix2);
	    key2.v_min = mesh2.kval(YFIXED, v_ix2);
#endif
	    it2 = emap_.find(key2);

#ifdef DEBUG
	    if (it2 == emap_.end())
	      std::cout << "LRSplineSurface::refine : Element not found" << std::endl;
#endif
	  }


	int u_ix = (d == XFIXED) ? fixed_ix : i;
	int v_ix = (d == YFIXED) ? fixed_ix : i;
	ElementMap::key_type key = {mesh_.kval(XFIXED, u_ix),
				    mesh_.kval(YFIXED, v_ix)};
	auto it = emap_.find(key);

	vector<double> data_points;
	vector<double> ghost_points;
	vector<double> significant_points;
	bool sort_in_u, sort_in_u_significant, sort_in_u_ghost;
	//double maxerr, averr, accerr;
	int nmbout;
	int pt_del;

	if (it2 != emap_.end())
	  {
	    // Update size of existing element
	    Mesh2DIterator m(mesh_, u_ix2, v_ix2);
	    it2->second->setUmax(mesh_.kval(XFIXED, (*m)[2]));
	    it2->second->setVmax(mesh_.kval(YFIXED, (*m)[3]));

	    // Fetch scattered data from the element that no longer is
	    // inside
	    it2->second->getOutsidePoints(data_points, d, sort_in_u);
	    it2->second->getOutsideSignificantPoints(significant_points, 
						     d, sort_in_u_significant);
	    it2->second->getOutsideGhostPoints(ghost_points, d, 
					       sort_in_u_ghost);
	    pt_del = it2->second->getNmbValPrPoint();

	    // Update supported LRBsplines
	    for (size_t kb=0; kb<bsplines_affected.size(); ++kb)
	      {
		if (!bsplines_affected[kb]->overlaps(it2->second.get()))
		  {
#if 0//ndef NDEBUG
//		    std::cout << "DEBUG: kb = " << kb << std::endl;
#endif
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
	    Mesh2DIterator m(mesh_, u_ix, v_ix);
	    unique_ptr<Element2D> elem(new Element2D(mesh_.kval(XFIXED, (*m)[0]),
						     mesh_.kval(YFIXED, (*m)[1]),
						     mesh_.kval(XFIXED, (*m)[2]),
						     mesh_.kval(YFIXED, (*m)[3])));

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
				  sort_in_u, pt_del);
	    if (significant_points.size() > 0)
	      elem->addSignificantPoints(significant_points.begin(), 
					 significant_points.end(),
					 sort_in_u_significant, pt_del);
	      
	    if (ghost_points.size() > 0)
	      elem->addGhostPoints(ghost_points.begin(), ghost_points.end(),
				   sort_in_u_ghost, pt_del);
	    emap_.insert(std::make_pair(key, std::move(elem)));
	    //auto it3 = emap_.find(key);

	  }

      }
    }
  }
#ifdef DEBUG
  //std::cout << "Num elements post: " << numElements() << std::endl;
  std::ofstream refsf("refine_one_sf.g2");
  this->writeStandardHeader(refsf);
  this->write(refsf);
  int stp = 1;
#endif
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

  int compare_refs_sf(refval r1, refval r2)
  {
    return (r1.kval < r2.kval);
  }


//==============================================================================
  void LRSplineSurface::refine(const vector<Refinement2D>& refs, 
			     bool absolute)
//==============================================================================
{
#if 0//ndef NDEBUG
  {
    vector<LRBSpline2D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    //puts("Remove when done debugging!");
    int stop_break = 1;
  }
#endif
  // Collect B-splines affected by the refinements
  vector<LRBSpline2D*> affected;
  LRSplineUtils::get_affected_bsplines(refs, emap_, knot_tol_, mesh_, 
  				       affected);

  vector<refval> u_refs, v_refs;
  for (size_t i = 0; i != refs.size(); ++i) {
    bool refined;
    const Refinement2D& r = refs[i];
    const auto indices = // tuple<int, int, int, int>
      LRSplineUtils::refine_mesh(r.d, 
				 r.kval, 
				 r.start, 
				 r.end, 
				 r.multiplicity, 
				 absolute,
				 degree(r.d), 
				 knot_tol_, 
				 mesh_, 
				 (r.d == XFIXED) ?  bsplinesuni1_ : bsplinesuni2_,
				 refined);

    if (!refined)
      continue;
    refval curr(r.kval, r.multiplicity);
    size_t kr;
    if (r.d == XFIXED)
      {
	for (kr=0; kr<u_refs.size(); ++kr)
	  if (u_refs[kr].kval == curr.kval && 
	      u_refs[kr].mult == curr.mult)
	    break;
	if (kr == u_refs.size())
	  u_refs.push_back(curr);
      }
    else if (r.d == YFIXED)
      {
	for (kr=0; kr<v_refs.size(); ++kr)
	  if (v_refs[kr].kval == curr.kval && 
	      v_refs[kr].mult == curr.mult)
	    break;
	if (kr == v_refs.size())
	  v_refs.push_back(curr);
      }
  }

  if (u_refs.size() + v_refs.size() == 0)
    return; // Mesh rectangles already existing
  
  //std::cout << "Post refine mesh" << std::endl;
  std::sort(u_refs.begin(), u_refs.end(), compare_refs_sf);
  std::sort(v_refs.begin(), v_refs.end(), compare_refs_sf);

  for (int ki=(int)u_refs.size()-1; ki>=0; --ki)
    {
      int kval = u_refs[ki].kval;
      int fixed_ix = Mesh2DUtils::last_nonlarger_knotvalue_ix(mesh_, XFIXED, kval);

      int last_ix = 
	BSplineUniUtils::last_overlapping_bsplineuni(fixed_ix, bsplinesuni1_);
    
      LRSplineUtils::split_univariate(bsplinesuni1_, last_ix, fixed_ix, 
				      absolute? u_refs[ki].mult : 1);
    }

   for (int ki=(int)v_refs.size()-1; ki>=0; --ki)
    {
      int kval = v_refs[ki].kval;
      int fixed_ix = Mesh2DUtils::last_nonlarger_knotvalue_ix(mesh_, YFIXED, kval);

      int last_ix = 
	BSplineUniUtils::last_overlapping_bsplineuni(fixed_ix, bsplinesuni2_);
    
      LRSplineUtils::split_univariate(bsplinesuni2_, last_ix, fixed_ix, 
				      absolute? v_refs[ki].mult : 1); 
    }
 
  //   // Not efficient
  //   int fixed_ix = get<1>(indices); // Index of fixed_val in global knot vector.
  //   int last_ix = 
  //     BSplineUniUtils::last_overlapping_bsplineuni(fixed_ix,
  // 						   (r.d == XFIXED) ? bsplinesuni1_ : bsplinesuni2_);
    
  //   if (r.d == XFIXED)
  //     LRSplineUtils::split_univariate(bsplinesuni1_, last_ix, fixed_ix, 
  // 				      (absolute) ? r.multiplicity : 1);
  //   else
  //     LRSplineUtils::split_univariate(bsplinesuni2_, last_ix, fixed_ix, 
  // 				      (absolute) ? r.multiplicity : 1);
  // }


  //std::wcout << "Preparing for iterative splitting." << std::endl;
  //vector<unique_ptr<LRBSpline2D> > affected;
//   vector<LRBSpline2D*> affected;
//   affected.reserve(bsplines_.size());
// //  for_each(bsplines_.begin(), bsplines_.end(), [&](const BSplineMap::value_type& b) {
//   for (auto it = bsplines_.begin(); it!= bsplines_.end(); ++it)
//     {
//       affected.push_back((*it).second.get());
//       // @@@ VSK. This is maybe the place to remove element information from the bsplines?
//       // unique_ptr<LRBSpline2D> ptr = std::move(it->second);
//       // affected.emplace_back(std::move(ptr));//b.second);
//     };
  
  // @@@ VSK. In this case, we should not bother about splitting elements. They will
  // be regenerated later. Thus, the bsplines should NOT be updated with elements during
  // splitting
  // The bsplines should not have any pointers to elements. They will be set later
  //std::wcout << "Iteratively splitting." << std::endl;

  LRSplineUtils::iteratively_split2(affected, mesh_, bsplines_, NULL, 
				    bsplinesuni1_, bsplinesuni2_, false);
  //bsplines_.clear();

  //std::wcout << "Splitting finished, now inserting resulting functions" << std::endl;
  // The bsplines are checked for duplicates and inserted in the global bspline map
//  for_each(affected.begin(), affected.end(), [&](unique_ptr<LRBSpline2D> b) {
// for (auto it = affected.begin(); it != affected.end(); ++it)
//   {
//     LRSplineUtils::insert_basis_function(*it, mesh_, bsplines_);
//   };

#if 0//ndef NDEBUG
  {
    vector<LRBSpline2D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    puts("Remove when done debugging!");
  }
#endif

    // Remove unused univariate B-splines
  for (int i=(int)bsplinesuni1_.size()-1; i>=0; --i)
    if (bsplinesuni1_[i]->getCount() <= 0)
      bsplinesuni1_.erase(bsplinesuni1_.begin()+i);

  for (int i=(int)bsplinesuni2_.size()-1; i>=0; --i)
    if (bsplinesuni2_[i]->getCount() <= 0)
      bsplinesuni2_.erase(bsplinesuni2_.begin()+i);

  //std::wcout << "Finally, reconstructing element map." << std::endl;
  emap_ = construct_element_map_(mesh_, bsplines_); // reconstructing the emap once at the end
  curr_element_ = NULL;  // No valid any more
  //std::wcout << "Refinement now finished. " << std::endl;
#if 0//ndef NDEBUG
  {
    vector<LRBSpline2D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    puts("Remove when done debugging!");
  }
#endif

}

//==============================================================================
  void LRSplineSurface::addSurface(const LRSplineSurface& other_sf, double fac)
//==============================================================================
{
  double tol = 1.0e-12;  // Numeric noice
  int dim = dimension();

  if (numBasisFunctions() != other_sf.numBasisFunctions())
    THROW("Different number of basis functions in addSurface()");

  if (numElements() != other_sf.numElements())
    THROW("Different number of elements in addSurface()");

  // Update all basis functions
  BSplineMap::const_iterator it1 = basisFunctionsBegin();
  BSplineMap::const_iterator it2 = other_sf.basisFunctionsBegin();
  for (; it1 != basisFunctionsEnd(); ++it1, ++it2)
    {
      // @@@ Domain checking should be done. Postponed.
      
      // it is assumed that the gamma scaling factor is identical for both 
      //basis functions (should
      if (fabs(it1->second->gamma() - it2->second->gamma()) > tol)
	THROW("Not corresponding scaling factors in addSurface()");

      Point coef = 
	(it1->second->coefTimesGamma() + fac*it2->second->coefTimesGamma()) / 
	it1->second->gamma();

      setCoef(coef, it1->second.get());
    }    
}

//==============================================================================
void LRSplineSurface::to3D()
//==============================================================================
{
  if (dimension() != 1) 
    THROW("Member method 'to3D()' only applies to one-dimensional LR-splines");
  if (degree(XFIXED) == 0 || degree(YFIXED) == 0) 
    THROW("Cannot convert a 0-degree spline to 3D.");

  //LRSplineUtils::insertParameterFunctions(this);
  for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) {
    const double x = LRSplineUtils::compute_greville(b->second->kvec(XFIXED), 
						     mesh().knotsBegin(XFIXED));
    const double y = LRSplineUtils::compute_greville(b->second->kvec(YFIXED), 
						      mesh().knotsBegin(YFIXED));
    const double z_gamma = b->second->coefTimesGamma()[0];
    const double gamma = b->second->gamma();
    b->second->coefTimesGamma() = Point(x*gamma, y*gamma, z_gamma);
//    b->second->coefTimesGamma() = Point(x, y, z_gamma);
    //wcout << b.second.coefTimesGamma() << std::endl;
    // int dim = b->second->coefTimesGamma().size();
    // double z = z_gamma/gamma;
    // std::cout << "z: " << z << std::endl;
  }
  int dim = dimension();
  // std::cout << "Global dim: " << dim << std::endl;
  // if (rational())
  // {
  //     std::cout << "Rational!" << std::endl;
  // }
}


//==============================================================================
LineCloud LRSplineSurface::getElementBds(int num_pts) const
//==============================================================================
{
  int dim = dimension();
  vector<double> pts;
  // For each element
  for (LRSplineSurface::ElementMap::const_iterator it=elementsBegin();
       it != elementsEnd(); ++it)
    {
      double umin = it->second->umin();
      double vmin = it->second->vmin();
      double umax = it->second->umax();
      double vmax = it->second->vmax();

      // Left side
      double del = (vmax - vmin)/(double)(num_pts-1);
      int ki;
      double upar = umin;
      double vpar;
      Point pos;
      for (ki=0, vpar=vmin; ki<num_pts; ++ki, vpar+=del)
	{
	  vpar = std::min(vpar, vmax);
	  point(pos, upar, vpar, it->second.get());
	  if (dim == 1)
	    {
	      pts.push_back(upar);
	      pts.push_back(vpar);
	    }
	  pts.insert(pts.end(), pos.begin(), pos.end());
	  if (ki>0 && ki<num_pts-1)
	    {
	      if (dim == 1)
		{
		  pts.push_back(upar);
		  pts.push_back(vpar);
		}
	      pts.insert(pts.end(), pos.begin(), pos.end());
	    }
	}

      // Right side
      upar = umax;
      for (ki=0, vpar=vmin; ki<num_pts; ++ki, vpar+=del)
	{
	  vpar = std::min(vpar, vmax);
	  point(pos, upar, vpar, it->second.get());
	  if (dim == 1)
	    {
	      pts.push_back(upar);
	      pts.push_back(vpar);
	    }
	  pts.insert(pts.end(), pos.begin(), pos.end());
	  if (ki>0 && ki<num_pts-1)
	    {
	      if (dim == 1)
		{
		  pts.push_back(upar);
		  pts.push_back(vpar);
		}
	      pts.insert(pts.end(), pos.begin(), pos.end());
	    }
	}
      
      // Bottom
      vpar = vmin;
      del = (umax - umin)/(double)(num_pts-1);
     for (ki=0, upar=umin; ki<num_pts; ++ki, upar+=del)
	{
	  upar = std::min(upar, umax);
	  point(pos, upar, vpar, it->second.get());
	  if (dim == 1)
	    {
	      pts.push_back(upar);
	      pts.push_back(vpar);
	    }
	  pts.insert(pts.end(), pos.begin(), pos.end());
	  if (ki>0 && ki<num_pts-1)
	    {
	      if (dim == 1)
		{
		  pts.push_back(upar);
		  pts.push_back(vpar);
		}
	      pts.insert(pts.end(), pos.begin(), pos.end());
	    }
	}

      // Top
      vpar = vmax;
     for (ki=0, upar=umin; ki<num_pts; ++ki, upar+=del)
	{
	  upar = std::min(upar, umax);
	  point(pos, upar, vpar, it->second.get());
	  if (dim == 1)
	    {
	      pts.push_back(upar);
	      pts.push_back(vpar);
	    }
	  pts.insert(pts.end(), pos.begin(), pos.end());
	  if (ki>0 && ki<num_pts-1)
	    {
	      if (dim == 1)
		{
		  pts.push_back(upar);
		  pts.push_back(vpar);
		}
	      pts.insert(pts.end(), pos.begin(), pos.end());
	    }
	}
    }
  int dim2 = (dim == 1) ? 3 : dim;
  LineCloud lines(&pts[0], (int)pts.size()/(2*dim2));
  return lines;
}

 //==============================================================================
bool LRSplineSurface::rational() const
//==============================================================================
{
  return rational_;
}


//==============================================================================
void LRSplineSurface::translate(const Point& vec)
//==============================================================================
{
    assert(vec.size() == dimension());

    // We run through all coefs and translate the coef by the given vec.
    for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) {
	const Point coef = b->second->Coef();
	Point trans_coef = coef + vec;
	const double gamma = b->second->gamma();	
	b->second->setCoefAndGamma(trans_coef, gamma);
    }
}


//==============================================================================
void LRSplineSurface::expandToFullTensorProduct()
//==============================================================================
{
  //std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - copying mesh..." << std::endl;
  Mesh2D tensor_mesh = mesh_;
  
  //std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - setting uniform meshlines..." << std::endl;
  const vector<int> xmults = LRSplineUtils::set_uniform_meshlines(XFIXED, 
								  tensor_mesh);
  const vector<int> ymults = LRSplineUtils::set_uniform_meshlines(YFIXED, 
								  tensor_mesh);


  BSplineMap tensor_bsplines;
  //std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - identify elements from mesh..." << std::endl;
  ElementMap emap = LRSplineUtils::identify_elements_from_mesh(tensor_mesh);
  //std::wcout << "Size of emap: " << emap.size() << std::endl;
  //std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - splitting up basis functions..." << std::endl;
  // splitting up basis functions
  for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) {
    LRSplineUtils::tensor_split(b->second, 
				xmults, 
				ymults, 
				tensor_mesh,
				bsplinesuni1_,
				bsplinesuni2_,
				tensor_bsplines);
  }

  // registering all the produced functions with the elements
  //std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - registering produced functions..." << std::endl;
  //std::wcout << "Number of basis functions: " << tensor_bsplines.size() << std::endl;
  //std::wcout << "Number of elements: "<< tensor_mesh.numDistinctKnots(XFIXED)-1 << " x " ;
  //std::wcout << tensor_mesh.numDistinctKnots(YFIXED)-1 << std::endl;

  // @@@ VSK. Use information in the LRB-splines or regenerate all elements ?
  for (auto b = tensor_bsplines.begin(); b != tensor_bsplines.end(); ++b)  {
    LRSplineUtils::update_elements_with_single_bspline(b->second.get(), emap, 
						       tensor_mesh, false);
  }

  //std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - swapping and exiting." << std::endl;
  mesh_.swap(tensor_mesh);
  bsplines_.swap(tensor_bsplines);
  emap_.swap(emap);
}


//==============================================================================
Point LRSplineSurface::operator()(double u, double v, int u_deriv, int v_deriv) const
//==============================================================================
{
#if 0//ndef NDEBUG
  if (rational_ && (u_deriv + v_deriv > 1))
    {
      MESSAGE("Currently the sum of the derivatives should be at most 1.");
      Point ret_pt(dimension());
      ret_pt.setValue(0.0);
      return ret_pt;
    }
#endif

  if (u < paramMin(XFIXED)) {
    //MESSAGE("u was outside domain: " << u << " < " << paramMin(XFIXED) << ", moved inside.");
      u = paramMin(XFIXED);
  } else if (u > paramMax(XFIXED)) {
    //MESSAGE("u was outside domain: " << u << " > " << paramMax(XFIXED) << ", moved inside.");
      u = paramMax(XFIXED);
  }

  if (v < paramMin(YFIXED)) {
    //MESSAGE("v was outside domain: " << v << " < " << paramMin(YFIXED) << ", moved inside.");
      v = paramMin(YFIXED);
  } else if (v > paramMax(YFIXED)) {
    //MESSAGE("v was outside domain: " << v << " > " << paramMax(YFIXED) << ", moved inside.");
      v = paramMax(YFIXED);
  }

  // const bool u_on_end = (u == mesh_.maxParam(XFIXED));
  // const bool v_on_end = (v == mesh_.maxParam(YFIXED));
  // vector<LRBSpline2D*> covering_B_functions = 
  //   basisFunctionsWithSupportAt(u, v);
  Element2D* elem;
  if (curr_element_ && curr_element_->contains(u, v))
    elem = curr_element_;
  else
    {
      bool found = false;

      // Check neighbours
      if (curr_element_)
	{
	  vector<LRBSpline2D*> bsupp = curr_element_->getSupport();
	  std::set<Element2D*> supp_el;
	  for (size_t ka=0; ka<bsupp.size(); ++ka)
	    {
	      vector<Element2D*> esupp = bsupp[ka]->supportedElements();
	      supp_el.insert(esupp.begin(), esupp.end());
	    }
	  vector<Element2D*> supp_el2(supp_el.begin(), supp_el.end());
	  for (size_t ka=0; ka<supp_el2.size(); ++ka)
	    if (supp_el2[ka]->contains(u, v))
	      {
		elem = curr_element_ = supp_el2[ka];
		found = true;
		break;
	      }
	}
      if (!found)
	{
	  //std::cout << "Finding element for parameter value (" << u << "," << v << ")" << std::endl;
	  elem = coveringElement(u, v);
	  curr_element_ = (Element2D*)elem;
	}
    }
  return operator()(u, v, u_deriv, v_deriv, elem);
}


//==============================================================================
  Point LRSplineSurface::operator()(double u, double v, int u_deriv, int v_deriv,
				    Element2D* elem) const
//==============================================================================
{
  // Check element
  if (!elem || !elem->contains(u, v))
    {
      bool found = false;

      // Check neighbours
      if (elem)
	{
	  vector<LRBSpline2D*> bsupp = elem->getSupport();
	  std::set<Element2D*> supp_el;
	  for (size_t ka=0; ka<bsupp.size(); ++ka)
	    {
	      vector<Element2D*> esupp = bsupp[ka]->supportedElements();
	      supp_el.insert(esupp.begin(), esupp.end());
	    }
	  vector<Element2D*> supp_el2(supp_el.begin(), supp_el.end());
	  for (size_t ka=0; ka<supp_el2.size(); ++ka)
	    if (supp_el2[ka]->contains(u, v))
	      {
		elem = supp_el2[ka];
		found = true;
		break;
	      }
	}
      if (!found)
	{
	  //std::cout << "Finding element for parameter value (" << u << "," << v << ")" << std::endl;
	  elem = coveringElement(u, v);
	}
    }
  
  curr_element_ = elem;
  const vector<LRBSpline2D*>& covering_B_functions = elem->getSupport();

  Point result(this->dimension()); 
  result.setValue(0.0); // will be initialized to 0, with the correct dimension

  // loop over LR B-spline functions
  int ki=0;

  // Distinguish between rational and non-rational to avoid
  // making temporary storage in the non-rational case
  double eps = 1.0e-12;
  const bool u_on_end = (u >= mesh_.maxParam(XFIXED)-eps); //(u == (*b)->umax());
  const bool v_on_end = (v >= mesh_.maxParam(YFIXED)-eps); // (v == (*b)->vmax());

  if (!rational_)
    {
      for (auto b = covering_B_functions.begin(); 
	   b != covering_B_functions.end(); ++b, ++ki) 
	{
	  // The b-function contains the coefficient.
	  result += (*b)->eval(u, 
			       v, 
			       u_deriv, 
			       v_deriv, 
			       u_on_end, 
			       v_on_end);
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
	  // const bool u_on_end = (u == (*b)->umax());
	  // const bool v_on_end = (v == (*b)->vmax());

	  // The b-function contains the coefficient.
	  double basis_val_pos = (*b)->evalBasisFunction(u, 
							 v, 
							 0, 
							 0, 
							 u_on_end, 
							 v_on_end);
	  //double gamma = (*b)->gamma();
	  double weight = (*b)->weight();
	  Point coef = (*b)->coefTimesGamma();
	  
	  // This is the nominator-position.
	  nom_pos += coef*weight*basis_val_pos;
	  
	  denom_pos += weight*basis_val_pos;
	  
	  if (u_deriv > 0 || v_deriv > 0)
	    {
	      double basis_val_der = (*b)->evalBasisFunction(u, 
							     v, 
							     u_deriv, 
							     v_deriv, 
							     u_on_end, 
							     v_on_end);

	      // This is the nominator-deriv.
	      nom_der += coef*weight*basis_val_der;
	      
	      denom_der += weight*basis_val_der;
	    }

#ifndef NDEBUG
	  if (u_deriv > 0 || v_deriv > 0)
	    {
	      ;//MESSAGE("Do not think that rational derivs are supported yet.");
	      //	      denom = 1.0;
	    }
	  //	  std::cout << "denom: " << denom << std::endl;
	  // if (rat_den == 0.0)
	  //   rat_den = 1.0;
#endif
	  
	}

      if (u_deriv == 0 && v_deriv == 0)
	{
	  result = nom_pos/denom_pos;
	}
      else if (u_deriv + v_deriv == 1)
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

//===========================================================================
void LRSplineSurface::closestPoint(const Point& pt,
				   double& clo_u,
				   double& clo_v, 
				   Point& clo_pt,
				   double& clo_dist,
				   double epsilon,
				   int maxiter,
				   Element2D* elem,
				   const RectDomain* rd,
				   double *seed) const
//===========================================================================
{
  double seed_buf[2];
  if (!seed) {
    seed = seed_buf;
    if (elem)
      {
	seed[0] = 0.5*(elem->umin()+elem->umax());
	seed[1] = 0.5*(elem->vmin()+elem->vmax());
      }
    else
      {
	// no seed given, we must compute one
	LRSplineSurface *currsf = (LRSplineSurface*)this;
	LRBSpline2D *bspline = 
	  LRSplineUtils::mostComparableBspline(currsf, pt);
	Point guesspt = bspline->getGrevilleParameter();
	seed[0] = guesspt[0];
	seed[1] = guesspt[1];
      }
  }

  // Uses closest point iteration fetched from SISL. 
  int kstat = 0;
  double par[2], start[2], end[2];
  RectDomain dom;
  if (rd)
    {
      start[0] = rd->umin();
      start[1] = rd->vmin();
      end[0] = rd->umax();
      end[1] = rd->vmax();
    }
  else
    {
      dom = containingDomain();
      start[0] = dom.umin();
      start[1] = dom.vmin();
      end[0] = dom.umax();
      end[1] = dom.vmax();
    }
  s1773(pt.begin(), epsilon, start, end, seed, par, maxiter, elem, &kstat);
  clo_u = par[0];
  clo_v = par[1];
  point(clo_pt, clo_u, clo_v, elem);
  clo_dist = pt.dist(clo_pt);
   
}
// //==============================================================================
// void LRSplineSurface::plotMesh(std::wostream& os) const 
// //==============================================================================
// {
//   plot_mesh(mesh_);
// }

// //==============================================================================
// void LRSplineSurface::plotBasisFunctionSupports(std::wostream& os) const
// //==============================================================================
// {
//   for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) 
//     plot_bspline_function(mesh_, b->second);
// }

//==============================================================================
void LRSplineSurface::setCoef(const Point& value, const LRBSpline2D* target)
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
void LRSplineSurface::setCoefTimesGamma(const Point& value, const LRBSpline2D* target)
//==============================================================================
{
  const auto it = bsplines_.find(generate_key(*target, mesh_));
  if (it == bsplines_.end()) 
    THROW("setCoef:: 'target' argument does not refer to member basis function.");

  if (value.dimension() != this->dimension())
    THROW("setCoef:: incorrect dimension of 'value' argument.");

  // if we got here, calling contract is fulfilled
  it->second->coefTimesGamma() = value;
} 

//==============================================================================
void LRSplineSurface::setCoef(const Point& value, 
		       int umin_ix, 
		       int vmin_ix, 
		       int umax_ix, 
		       int vmax_ix, 
		       int u_mult, 
		       int v_mult)
//==============================================================================
{
  const BSKey key = {mesh_.kval(XFIXED, umin_ix), 
		     mesh_.kval(YFIXED, vmin_ix), 
		     mesh_.kval(XFIXED, umax_ix),
		     mesh_.kval(YFIXED, vmax_ix),
		     u_mult, 
		     v_mult};

  const auto it = bsplines_.find(key);
	                         
  if (it == bsplines_.end())
    THROW("setCoef:: There is no such basis function.");

  if (value.dimension() != this->dimension())
    THROW("setCoef:: incorrect dimension of 'value' argument.");
  
  // if we got here, calling contract is fulfilled
  const double gamma = it->second->gamma();
  it->second->coefTimesGamma() = value * gamma;
}

//==============================================================================
//
//  Functionality inherited from GeomObject or ParamSurface
//
//==============================================================================
//===========================================================================
ClassType LRSplineSurface::instanceType() const
//===========================================================================
{
  return classType();
}

  //===========================================================================
BoundingBox LRSplineSurface::boundingBox() const
  //===========================================================================
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
const RectDomain& LRSplineSurface::parameterDomain() const
  //===========================================================================
  {
    Array<double, 2> ll(mesh_.minParam(XFIXED), mesh_.minParam(YFIXED));
    Array<double, 2> ur(mesh_.maxParam(XFIXED), mesh_.maxParam(YFIXED));
    domain_ = RectDomain(ll, ur);
    return domain_;
  }

  //===========================================================================
  RectDomain LRSplineSurface::containingDomain() const
  //===========================================================================
  {
    RectDomain rect_dom = parameterDomain();
    return rect_dom;
  }

 //===========================================================================
  bool LRSplineSurface::inDomain(double u, double v, double eps) const
  //===========================================================================
  {
    if (u < startparam_u() || u > endparam_u())
	return false;
    if (v < startparam_v() || v > endparam_v())
	return false;

    return true;
  }

//===========================================================================
  int LRSplineSurface::inDomain2(double u, double v, double eps) const
//===========================================================================
{
    if (u < startparam_u()-eps || u > endparam_u()+eps)
	return 0;
    if (v < startparam_v()-eps || v > endparam_v()+eps)
	return 0;

    if (u < startparam_u()+eps || u > endparam_u()-eps)
	return 2;
    if (v < startparam_v()+eps || v > endparam_v()-eps)
	return 2;

    return 1;
}

//===========================================================================
  bool LRSplineSurface::onBoundary(double u, double v, double eps) const
//===========================================================================
{
  if ((u > startparam_u()-eps && u < startparam_u()+eps) || 
      (u > endparam_u()-eps && u < endparam_u()+eps))
	return true;
  if ((v > startparam_v()-eps && v < startparam_v()-eps) || 
      (v > endparam_v()-eps && v < endparam_v()+eps))
	return true;

    return false;
}
  //===========================================================================
  Point LRSplineSurface::closestInDomain(double u, double v) const
  //===========================================================================
  {
    double u1 = std::min(std::max(u, startparam_u()), endparam_u());
    double v1 = std::min(std::max(v, startparam_v()), endparam_v());
    return Point(u1, v1);
  }

  //===========================================================================
  CurveLoop LRSplineSurface::outerBoundaryLoop(double degenerate_epsilon) const
  //===========================================================================
  {
    // Test for degeneracy.
    bool deg[4];
    if (true /*degenerate_epsilon < 0.0*/)  // Degeneracy test not implemented yet
      deg[0] = deg[1] = deg[2] = deg[3] = false; // All curves are wanted
    else
      isDegenerate(deg[0], deg[1], deg[2], deg[3], degenerate_epsilon);
    std::vector< shared_ptr< ParamCurve > >  vec;
    int perm[4] = {2, 1, 3, 0};
    for (int edgenum = 0; edgenum < 4; ++edgenum) {
	if (!deg[edgenum]) {
	    shared_ptr<ParamCurve> edgecurve (edgeCurve(perm[edgenum]));
	    if (perm[edgenum] == 0 || perm[edgenum] == 3)
		edgecurve->reverseParameterDirection();
	    vec.push_back(edgecurve);
	}
    }

    return CurveLoop(vec, (degenerate_epsilon < 0.0) ? DEFAULT_SPACE_EPSILON :
		     degenerate_epsilon);
  }

  //===========================================================================
  vector<CurveLoop> LRSplineSurface::allBoundaryLoops(double degenerate_epsilon) const
  //===========================================================================
  {
    // There is only one boundary loop...
    std::vector<CurveLoop> cvloopvec;
    cvloopvec.push_back(outerBoundaryLoop(degenerate_epsilon));
    return cvloopvec;
  }

  //===========================================================================
  void LRSplineSurface::point(Point& pt, double upar, double vpar) const
  //===========================================================================
  {
    //pt = operator()(upar, vpar, 0, 0, curr_element_);
    if (!curr_element_ || !curr_element_->contains(upar, vpar))
      {
      bool found = false;

      // Check neighbours
      if (curr_element_)
	{
	  vector<LRBSpline2D*> bsupp = curr_element_->getSupport();
	  std::set<Element2D*> supp_el;
	  for (size_t ka=0; ka<bsupp.size(); ++ka)
	    {
	      vector<Element2D*> esupp = bsupp[ka]->supportedElements();
	      for (size_t ka=0; ka<esupp.size(); ++ka)
		if (esupp[ka]->contains(upar, vpar))
		  {
		    curr_element_ = esupp[ka];
		    found = true;
		    break;
		  }
	    }
	}
      if (!found)
	{
	  //std::cout << "Finding element for parameter value (" << u << "," << v << ")" << std::endl;
	  curr_element_ = coveringElement(upar, vpar);
	}
      }

    if (rational_)
      pt = operator()(upar, vpar, 0, 0, curr_element_);
    else
      {
	double eps = 1.0e-12;
	const bool u_on_end = (upar >= mesh_.maxParam(XFIXED)-eps); //(u == (*b)->umax());
	const bool v_on_end = (vpar >= mesh_.maxParam(YFIXED)-eps); // (v == (*b)->vmax());
	const vector<LRBSpline2D*>& bfunctions = curr_element_->getSupport();
	size_t bsize = bfunctions.size();
	//vector<BSplineUniLR*> uni(2*bsize, NULL);
	vector<double> val(2*bsize);
	pt.resize(this->dimension());
	pt.setValue(0.0);
	for (size_t ki=0; ki<bsize; ++ki)
	  {
	    //uni[ki] = bfunctions[ki]->getUnivariate(XFIXED);
	    //uni[bsize+ki] = bfunctions[ki]->getUnivariate(YFIXED);
	    size_t kj;
	    const BSplineUniLR* uni1 =  bfunctions[ki]->getUnivariate(XFIXED);
	    const BSplineUniLR* uni2 =  bfunctions[ki]->getUnivariate(YFIXED);
	    for (kj=0; kj<ki; ++kj)
	      //if (uni[kj] == uni[ki])
	      if (uni1 == bfunctions[kj]->getUnivariate(XFIXED))
		break;
	    if (kj < ki)
	      val[ki] = val[kj];
	    else
	      //val[ki] = uni[ki]->evalBasisFunction(upar, 0, u_on_end);
	      val[ki] = uni1->evalBasisFunction(upar, 0, u_on_end);

	    for (kj=0; kj<ki; ++kj)
	      //if (uni[bsize+kj] == uni[bsize+ki])
	      if (uni2 == bfunctions[kj]->getUnivariate(YFIXED))
		break;
	    if (kj < ki)
	      val[bsize+ki] = val[bsize+kj];
	    else
	      val[bsize+ki] = 
		//uni[bsize+ki]->evalBasisFunction(vpar, 0, v_on_end);
		uni2->evalBasisFunction(vpar, 0, v_on_end);
	    pt += val[ki]*val[bsize+ki]*bfunctions[ki]->coefTimesGamma();
	  }
      }
  }

  //===========================================================================
void LRSplineSurface::point(Point& pt, double upar, double vpar,
			    Element2D* elem) const
  //===========================================================================
  {
    pt = operator()(upar, vpar, 0, 0, elem);
  }

   //===========================================================================
  void LRSplineSurface::normal(Point& pt, double upar, double vpar) const
  //===========================================================================
  {
    normal(pt, upar, vpar, curr_element_);
  }

   //===========================================================================
  void LRSplineSurface::normal(Point& pt, double upar, double vpar,
			    Element2D* elem) const
  //===========================================================================
  {
    double tol = DEFAULT_SPACE_EPSILON;

    Point pt_der1; 
    Point pt_der2;
    if (elem == NULL)
      {
	pt_der1 = operator()(upar, vpar, 1, 0);
	pt_der2 = operator()(upar, vpar, 0, 1);
      }
    else
      {
	pt_der1 = operator()(upar, vpar, 1, 0, elem);
	pt_der2 = operator()(upar, vpar, 0, 1, elem);
      }
      

    pt = pt_der1.cross(pt_der2);
    double l = pt.length();
    if (l > tol)
      pt.normalize();

    double cross_tan_ang = pt_der1.angle_smallest(pt_der2);
    cross_tan_ang = std::min(cross_tan_ang, fabs(M_PI - cross_tan_ang));
    double cross_tan_ang_tol = 1e-03;
    if (l < tol || cross_tan_ang < cross_tan_ang_tol) {
      if (cross_tan_ang < cross_tan_ang_tol)
	   MESSAGE("Too small angle between cross tangents, "
		   "degenerate point!");
	pt.setValue(0.0);
//	return false;
    }

  }

//===========================================================================
  void LRSplineSurface::evalGrid(int num_u, int num_v, 
				 double umin, double umax, 
				 double vmin, double vmax,
				 std::vector<double>& points,
				 double nodata_val) const
//===========================================================================
  {
    int dim = dimension();
    points.reserve(num_u*num_v*dim);

    // // Make intermediate tensor product spline surface to speed up the
    // // grid evaluation
    // shared_ptr<LRSplineSurface> tmp = shared_ptr<LRSplineSurface>(this->clone());;
    // shared_ptr<SplineSurface> tpsf =
    //   shared_ptr<SplineSurface>(tmp->asSplineSurface());

    // vector<double> param_u;
    // vector<double> param_v;
    // tpsf->gridEvaluator(num_u, num_v, points, param_u, param_v,
    // 			umin, umax, vmin, vmax);
  // Construct mesh of element pointers
    vector<Element2D*> elements;
    constructElementMesh(elements);
    
    // Get all knot values in the u-direction
    const double* const uknots = mesh_.knotsBegin(XFIXED);
    const double* const uknots_end = mesh_.knotsEnd(XFIXED);
    int nmb_knots_u = mesh_.numDistinctKnots(XFIXED);
    const double* knotu;
    
  // Get all knot values in the v-direction
    const double* const vknots = mesh_.knotsBegin(YFIXED);
    const double* const vknots_end = mesh_.knotsEnd(YFIXED);
    int nmb_knots_v = mesh_.numDistinctKnots(YFIXED);
    const double* knotv;

    double udel = (umax - umin)/(double)(num_u-1);
    double vdel = (vmax - vmin)/(double)(num_v-1);
    double upar = umin;
    double vpar = vmin;
    double tolu = std::max(1.0e-8, 1.0e-8*udel);
    double tolv = std::max(1.0e-8, 1.0e-8*vdel);

#ifdef DEBUG
    std::ofstream of("tmp_grid.g2");
    (void)of.precision(15);
    of << "400 1 0 4 255 0 0 255" << std::endl;
    of << num_u*num_v << std::endl;
#endif

    int ki, kj, kr, kh;
    for (kj=0, kr=0, knotv=vknots, ++knotv; knotv!=vknots_end; ++knotv, ++kj)
    {
      int lastv = (knotv+1 == vknots_end);
      for (; kr<num_v && vpar <= (*knotv)+lastv*tolv; ++kr, vpar+=vdel)
	{
	  if (lastv)
	    vpar = std::min(vpar, *knotv);
	  upar = umin;
	  for (ki=0, kh=0, knotu=uknots, ++knotu; knotu != uknots_end; 
	       ++knotu, ++ki)
	    {
	      int lastu = (knotu+1 == uknots_end);
	      Element2D *elem = elements[kj*(nmb_knots_u-1)+ki];
	      for (; kh<num_u && upar <= (*knotu)+lastu*tolu; ++kh, upar+=udel)
		{
		  if (lastu)
		    upar = std::min(upar, *knotu);
		  Point pos;
		  point(pos, upar, vpar, elem);
		  points.insert(points.end(), pos.begin(), pos.end());

#ifdef DEBUG
		  of << upar << " " << vpar << " " << pos[0] << std::endl;
#endif
		}
	    }
#ifdef DEBUG
	  int stop_break = 1;
#endif
	}
    }
  }

//===========================================================================
double LRSplineSurface::startparam_u() const
//===========================================================================
{
  return paramMin(XFIXED);
}

  //===========================================================================
double LRSplineSurface::endparam_u() const
  //===========================================================================
  {
    return paramMax(XFIXED);
  }

  //===========================================================================
double LRSplineSurface::startparam_v() const
  //===========================================================================
  {
    return paramMin(YFIXED);
  }

  //===========================================================================
double LRSplineSurface::endparam_v() const
  //===========================================================================
  {
    return paramMax(YFIXED);
  }

   //===========================================================================
  void LRSplineSurface::point(vector<Point>& pts, 
			      double upar, double vpar,
			      int derivs,
			      bool u_from_right,
			      bool v_from_right,
			      double resolution) const
  //===========================================================================
  {
    point(pts, upar, vpar, derivs, curr_element_, u_from_right,
	  v_from_right, resolution);
  }

   //===========================================================================
  void LRSplineSurface::point(vector<Point>& pts, 
			      double upar, double vpar,
			      int derivs,
			      Element2D* elem,
			      bool u_from_right,
			      bool v_from_right,
			      double resolution) const
  //===========================================================================
  {
    int totpts = (derivs + 1)*(derivs + 2)/2;
    DEBUG_ERROR_IF((int)pts.size() < totpts, "The vector of points must have sufficient size.");

    int dim = dimension();
    for (int ki = 0; ki < totpts; ++ki) {
	if (pts[ki].dimension() != dim) {
	    pts[ki].resize(dim);
	}
	pts[ki].setValue(0.0);
    }

  if (upar < paramMin(XFIXED)) {
    //MESSAGE("upar was outside domain: " << upar << " < " << paramMin(XFIXED) << ", moved inside.");
      upar = paramMin(XFIXED);
  } else if (upar > paramMax(XFIXED)) {
    //MESSAGE("upar was outside domain: " << upar << " > " << paramMax(XFIXED) << ", moved inside.");
      upar = paramMax(XFIXED);
  }

  if (vpar < paramMin(YFIXED)) {
    //MESSAGE("vpar was outside domain: " << vpar << " < " << paramMin(YFIXED) << ", moved inside.");
      vpar = paramMin(YFIXED);
  } else if (vpar > paramMax(YFIXED)) {
    //MESSAGE("vpar was outside domain: " << vpar << " > " << paramMax(YFIXED) << ", moved inside.");
      vpar = paramMax(YFIXED);
  }
    
  // Check element
  // if (!elem)
  //   elem = coveringElement(upar, vpar);
  if (!elem || !elem->contains(upar, vpar))
    {
      bool found = false;

      // Check neighbours
      if (elem)
	{
	  vector<LRBSpline2D*> bsupp = elem->getSupport();
	  std::set<Element2D*> supp_el;
	  for (size_t ka=0; ka<bsupp.size(); ++ka)
	    {
	      vector<Element2D*> esupp = bsupp[ka]->supportedElements();
	      supp_el.insert(esupp.begin(), esupp.end());
	    }
	  vector<Element2D*> supp_el2(supp_el.begin(), supp_el.end());
	  for (size_t ka=0; ka<supp_el2.size(); ++ka)
	    if (supp_el2[ka]->contains(upar, vpar))
	      {
		elem = supp_el2[ka];
		found = true;
		break;
	      }
	}
      if (!found)
	{
	  //std::cout << "Finding element for parameter value (" << u << "," << v << ")" << std::endl;
	  elem = coveringElement(upar, vpar);
	}
    }
  
  curr_element_ = elem;
  const vector<LRBSpline2D*>& covering_B_functions = elem->getSupport();

  //vector<Point> tmp(totpts, Point(dim));

  // vector<Point> pts2(totpts+1);
  // double eps = 1.0e-12;
  // const bool u_on_end = (upar >= mesh_.maxParam(XFIXED)-eps); //(u == (*b)->umax());
  // const bool v_on_end = (vpar >= mesh_.maxParam(YFIXED)-eps); // (v == (*b)->vmax());

  for (size_t kr=0; kr<covering_B_functions.size(); ++kr)
    {
      covering_B_functions[kr]->evalder_add(upar, vpar, derivs, &pts[0], 
					    u_from_right, v_from_right);
      // covering_B_functions[kr]->evalder(upar, vpar, derivs, &tmp[0], 
      // 					u_from_right, v_from_right);
      // for (size_t kh=0; kh<pts.size(); ++kh)
      // 	pts[kh] += tmp[kh];
    }
  // This is not the most efficient approach, should be faster to
  // evaluate basis functions once. Only a first implementation.
  // int cntr = 0;
  // for (int kj = 0; kj < derivs + 1; ++kj)
  //   {
  //     for (int ki = 0; ki < kj + 1; ++ki, ++cntr)
  // 	{
  //   	  pts2[cntr] = operator()(upar, vpar, kj-ki, ki, elem);
  // 	}
  //   }
  // if (pts2[0].dist(pts[0]) > 5)
  //   {
  //     int stop_break = 1;
  //   }
  }

  //===========================================================================
  DirectionCone LRSplineSurface::normalCone() const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::normalCone() not implemented yet");
    DirectionCone dc;
    return dc;
  }

  //===========================================================================
  DirectionCone LRSplineSurface::tangentCone(bool pardir_is_u) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::tangentCone() not implemented yet");
    DirectionCone dc;
    return dc;
  }

  //===========================================================================
  CompositeBox LRSplineSurface::compositeBox() const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::compositeBox() not implemented yet");
    // The CompositeBox api must be extended to set inner & outer box.
    CompositeBox cb(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0));
    return cb;
  }

  //===========================================================================
  LRSplineSurface*
    LRSplineSurface::subSurface(double from_upar, double from_vpar,
				 double to_upar, double to_vpar,
				double fuzzy) const
  //===========================================================================
  {
    // Check input
    if (from_upar >= to_upar) {
	THROW("First u-parameter must be strictly less than second.");
    }
    if (from_vpar >= to_vpar) {
	THROW("First v-parameter must be strictly less than second.");
    }
    if (from_upar < startparam_u()-fuzzy || from_vpar < startparam_v()-fuzzy) {
	THROW("Subsurface defined outside surface.");
    }

    // Periodic surfaces not considered, i.e. it is not possible to pick a
    // part of a surface over a periodic seam.

    // If boundaries are close to existing knots, we snap.
    // Note that the indices refer to the vectors of unique knots.
    int ix1 = mesh_.knotIntervalFuzzy(XFIXED, from_upar, fuzzy);
    int ix2 = mesh_.knotIntervalFuzzy(XFIXED, to_upar, fuzzy);
    int iy1 = mesh_.knotIntervalFuzzy(YFIXED, from_vpar, fuzzy);
    int iy2 = mesh_.knotIntervalFuzzy(YFIXED, to_vpar, fuzzy);

    LRSplineSurface *surf = NULL;

    int deg1 = degree(XFIXED);
    int deg2 = degree(YFIXED);
    int nmb1 = mesh_.numDistinctKnots(XFIXED);
    int nmb2 = mesh_.numDistinctKnots(YFIXED);

    // Check if the sub surface is identical to the current surface
    bool identical = true;
    if (ix1 != 0 || mesh_.kval(XFIXED, ix1) != from_upar ||
	mesh_.minMultInLine(XFIXED, ix1) < deg1+1)
      identical = false;
    if (ix2 != mesh_.lastMeshVecIx(XFIXED) || mesh_.kval(XFIXED, ix2) != to_upar ||
	mesh_.minMultInLine(XFIXED, ix2) < deg1+1)
      identical = false;
    if (iy1 != 0 || mesh_.kval(YFIXED, iy1) != from_vpar ||
	mesh_.minMultInLine(YFIXED, iy1) < deg2+1)
      identical = false;
     if (iy2 != mesh_.lastMeshVecIx(YFIXED) || mesh_.kval(YFIXED, iy2) != to_vpar ||
	mesh_.minMultInLine(YFIXED, iy2) < deg2+1)
      identical = false;
    
     if (identical)
       {
	 // The sub surface is equal to the current. Copy.
	 surf = new LRSplineSurface(*this);
	 return surf;
       }

     // Make a copy of the current surface
     shared_ptr<LRSplineSurface> sf(new LRSplineSurface(*this));

     // We locate all basis functions for which the support is crossed
     // by a subsurface boundary line.
     // We store the min & max index in each dir of those functions.
     // @@sbr201301 For large cases It may be faster to go through the elements.
     int umin_ind = nmb1; // Higher than max index.
     int umax_ind = 0;
     int vmin_ind = nmb2; // Higher than max index.
     int vmax_ind = 0;
     for (auto iter = basisFunctionsBegin(); iter != basisFunctionsEnd(); ++iter)
       {
	 LRBSpline2D* bas_func = iter->second.get();

	 if ((from_upar >= bas_func->umin()) && (from_upar <= bas_func->umax()))
	   umin_ind = std::min(umin_ind, bas_func->suppMin(XFIXED));
	 if (to_upar >= bas_func->umin() && to_upar <= bas_func->umax())
	   umax_ind = std::max(umax_ind, bas_func->suppMax(XFIXED));

	 if (from_vpar >= bas_func->vmin() && from_vpar <= bas_func->vmax())
	   vmin_ind = std::min(vmin_ind, bas_func->suppMin(YFIXED));
	 if (to_vpar >= bas_func->vmin() && to_vpar <= bas_func->vmax())
	   vmax_ind = std::max(vmax_ind, bas_func->suppMax(YFIXED));
       }

     // Define possible new knotlines
     // Note that the new knot lines must be longer than the size of the
     // sub surface to avoid LR B-splines partly overlapping the sub domain
     // @@@ VSK. Could the extension be smaller than prescribed here?
     vector<Refinement2D> refs(4);
     //vector<Refinement2D> refs;
#if 0 // Old version, we need to consider the support of all the basis
      // functions.
     double umin = mesh_.kval(XFIXED, std::max(ix1 - deg1, 0));
     double umax = mesh_.kval(XFIXED, std::min(ix2 + deg1, nmb1-1));
     double vmin = mesh_.kval(YFIXED, std::max(iy1 - deg2, 0));
     double vmax = mesh_.kval(YFIXED, std::min(iy2 + deg2, nmb2-1));
#else
     double umin = mesh_.kval(XFIXED, umin_ind);
     double umax = mesh_.kval(XFIXED, umax_ind);
     double vmin = mesh_.kval(YFIXED, vmin_ind);
     double vmax = mesh_.kval(YFIXED, vmax_ind);
#endif
     // if (umin_ind > 0)
     //   {
     // 	 Refinement2D curr;
     // 	 curr.setVal(from_upar, vmin, vmax, XFIXED, deg1+1);
     // 	 refs.push_back(curr);
     //   }
     // if (umax_ind < nmb1-1)
     //   {
     // 	 Refinement2D curr;
     // 	 curr.setVal(to_upar, vmin, vmax, XFIXED, deg1+1);
     // 	 refs.push_back(curr);
     //   }
     // if (vmin_ind > 0)
     //   {
     // 	 Refinement2D curr;
     // 	 curr.setVal(from_vpar, umin, umax, YFIXED, deg2+1);
     // 	 refs.push_back(curr);
     //   }
     // if (vmax_ind < nmb2-1)
     //   {
     // 	 Refinement2D curr;
     // 	 curr.setVal(to_vpar, umin, umax, YFIXED, deg2+1);
     // 	 refs.push_back(curr);
     //   }
     refs[0].setVal(from_upar, vmin, vmax, XFIXED, deg1+1);
     refs[1].setVal(to_upar, vmin, vmax, XFIXED, deg1+1);
     refs[2].setVal(from_vpar, umin, umax, YFIXED, deg2+1);
     refs[3].setVal(to_vpar, umin, umax, YFIXED, deg2+1);
     
     // Perform refinement
     // @@sbr201301 Remove when stable.
     bool multi_refine = false; //true; //false;
     if (multi_refine)
       {
	 sf->refine(refs, true);
       }
     else
       {
#ifndef NDEBUG
//	 puts("Debugging, remove when code is stable!");
//	 std::swap(refs[0], refs[1]);
#endif
	 for (size_t ki = 0; ki < refs.size(); ++ki)
	   {
#ifndef NDEBUG
	     MESSAGE("ki = " << ki << "\n");
#endif
	     sf->refine(refs[ki], true); // Second argument is 'true', which means that the mult is set	     
	                                 // to refs[ki].mult = deg+1.
	   }
       }

     if (false)
       {
     	 std::ofstream of("tmp_lr.g2");
     	 sf->writeStandardHeader(of);
     	 sf->write(of);
       }

     // Fetch sub mesh
     // First get knot indices
     int iu1 = sf->mesh().getKnotIdx(XFIXED, from_upar, fuzzy);
     int iu2 = sf->mesh().getKnotIdx(XFIXED, to_upar, fuzzy);
     int iv1 = sf->mesh().getKnotIdx(YFIXED, from_vpar, fuzzy);
     int iv2 = sf->mesh().getKnotIdx(YFIXED, to_vpar, fuzzy);
     if (iu1 < 0 || iu2 < 0 || iv1 < 0 || iv2 < 0)
       THROW("LRSplineSurface::subSurface. Knot line not existing");

     shared_ptr<Mesh2D> sub_mesh = sf->mesh().subMesh(iu1, iu2, iv1, iv2);

     // Fetch LR B-splines living on the sub mesh
     vector<LRBSpline2D*> b_splines = 
				     sf->collect_basis(iu1, iu2, iv1, iv2);

     // Copy LR B-splines and update indices to the new domain
     // It is not absolutely necessary to make copies since the intermediate
     // surface will die after this function, but it is done for consistency
     // vector<unique_ptr<LRBSpline2D> > b_splines2(b_splines.size());
     // for (size_t ki=0; ki<b_splines.size(); ++ki)
     //   {
     // 	 b_splines2[ki] = unique_ptr<LRBSpline2D>(new LRBSpline2D(*b_splines[ki]));
     // 	 b_splines2[ki]->setMesh(sub_mesh.get());
     // 	 b_splines2[ki]->subtractKnotIdx(iu1, iv1);
     //   }

     // Create sub surface
     //surf = new LRSplineSurface(knot_tol_, rational_, *sub_mesh, b_splines2);
     surf = new LRSplineSurface(knot_tol_, rational_, *sub_mesh, b_splines, iu1, iv1);
     
    return surf;
  }
#if 0
  //===========================================================================
  void LRSplineSurface::refineBasisFunction(int index)
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::refineBasisFunction() not implemented yet");
  }

  //===========================================================================
  void LRSplineSurface::refineBasisFunction(const std::vector<int> &indices)
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::refineBasisFunction() not implemented yet");
  }

  //===========================================================================
  void LRSplineSurface::refineElement(int index)
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::refineElement() not implemented yet");
  }

  //===========================================================================
  void LRSplineSurface::refineElement(const std::vector<int> &indices)
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::refineElement() not implemented yet");
  }
#endif
  //===========================================================================
  vector<shared_ptr<ParamSurface> >
    LRSplineSurface::subSurfaces(double from_upar, double from_vpar,
				 double to_upar, double to_vpar,
				 double fuzzy) const
  //===========================================================================
  {
    vector<shared_ptr<ParamSurface> > sub_sfs;
    sub_sfs.push_back(shared_ptr<ParamSurface>(subSurface(from_upar, from_vpar,
							 to_upar, to_vpar, fuzzy)));

    return sub_sfs;
  }

  //===========================================================================
  void LRSplineSurface::closestBoundaryPoint(const Point& pt,
					     double&        clo_u,
					     double&        clo_v, 
					     Point&         clo_pt,
					     double&        clo_dist,
					     double         epsilon,
					     const RectDomain* rd,
					     double *seed) const
  //===========================================================================
  {
    RectDomain domain = containingDomain();
    if (!rd)
	rd = &domain;

    Point cpt;
    double cdist, cpar;

    // Check degeneracy
    bool b, r, t, l;
    //   double tol = 0.000001;  // Arbitrary tolerance. The information should
    //                           // be present.
    (void)isDegenerate(b, r, t, l, epsilon);

    // Checking closest point on the bottom boundary
    shared_ptr<SplineCurve> bdcrv;
    clo_dist = 1.0e10;  // Initialize with a large number
    if (b == false)
	{
	    bdcrv = shared_ptr<SplineCurve>(constParamCurve(rd->vmin(), true));
	    bdcrv->closestPoint(pt, rd->umin(), rd->umax(),
				clo_u, clo_pt, clo_dist, seed);
	    clo_v = rd->vmin();
	}

    // Checking the right boundary
    if (r == false)
	{
	    bdcrv = shared_ptr<SplineCurve>(constParamCurve(rd->umax(),
							    false));
	    bdcrv->closestPoint(pt, rd->vmin(), rd->vmax(), cpar, cpt, cdist,
				(seed == 0) ? seed : seed+1);
	    if (cdist < clo_dist)
		{
		    clo_pt = cpt;
		    clo_u = rd->umax();
		    clo_v = cpar;
		    clo_dist = cdist;
		}
	}

    // Checking the upper boundary
    if (t == false)
	{
	    bdcrv = shared_ptr<SplineCurve>(constParamCurve(rd->vmax(), true));
	    bdcrv->closestPoint(pt, rd->umin(), rd->umax(), cpar, cpt, cdist,
				seed);
	    if (cdist < clo_dist)
		{
		    clo_pt = cpt;
		    clo_u = cpar;
		    clo_v = rd->vmax();
		    clo_dist = cdist;
		}
	}

    // Checking the left boundary
    if (l == false)
	{
	    bdcrv = shared_ptr<SplineCurve>(constParamCurve(rd->umin(),
							    false));
	    bdcrv->closestPoint(pt, rd->vmin(), rd->vmax(), cpar, cpt, cdist,
				(seed == 0) ? seed : seed+1);
	    if (cdist < clo_dist)
		{
		    clo_pt = cpt;
		    clo_u = rd->umin();
		    clo_v = cpar;
		    clo_dist = cdist;
		}
	}


  }

  //===========================================================================
  void LRSplineSurface::getBoundaryInfo(Point& pt1, Point& pt2, 
					double epsilon, SplineCurve*& cv,
					SplineCurve*& crosscv, double knot_tol) const
  //===========================================================================
  {
    double par1, par2;
    int bdidx;
    getBoundaryIdx(pt1, pt2, epsilon, bdidx, par1, par2, knot_tol);
    if (bdidx < 0)
      return;

    getBoundaryInfo(par1, par2, bdidx, cv, crosscv, knot_tol);
    return;
  }

  //===========================================================================
  void
  LRSplineSurface::getBoundaryInfo(double par1, double par2,
				   int bdindex, SplineCurve*& cv,
				   SplineCurve*& crosscv, double knot_tol) const
  //===========================================================================
  {
    // Set no output
    cv = crosscv = 0;
    // Following intuition, the curves go from pt1 to pt2.
    bool turn_curves = false;

    if (par1 > par2)
      turn_curves = true;
    double bdpar;
    switch (bdindex)
      {
      case 0:
	bdpar = startparam_v();
	break;
      case 1:
	bdpar = endparam_v();
	break;
      case 2:
	bdpar = startparam_u();
	break;
      case 3:
	bdpar = endparam_u();
	break;
      default:
	return;
      }

    // Ftech constant parameter curve in par. dir.
    SplineCurve *c1=0, *c2=0;
    constParamCurve(bdpar, bdindex<=1, c1, c2);
// #ifdef _MSC_VER
//   cv = dynamic_cast<SplineCurve*>(c1->subCurve(std::min(par1,par2),
// 						     std::max(par1,par2), knot_tol));
//   crosscv = dynamic_cast<SplineCurve*>(c2->subCurve(std::min(par1,par2),
// 							  std::max(par1,par2), knot_tol));
// #else
    cv = c1->subCurve(std::min(par1,par2), std::max(par1,par2), knot_tol);
    crosscv = c2->subCurve(std::min(par1,par2), std::max(par1,par2), knot_tol);
// #endif
    if (bdindex == 0 || bdindex == 2)
      // We must turn cross-curve so that it points outwards.
      for (std::vector<double>::iterator iter = crosscv->coefs_begin();
	   iter != crosscv->coefs_end(); ++iter)
	iter[0] *= -1.0;
    delete c1;
    delete c2;

    if (turn_curves) {
      cv->reverseParameterDirection();
      crosscv->reverseParameterDirection();
    }
    // else not a boundary curve. Return no curves.



  }

  //===========================================================================
  void
  LRSplineSurface::getBoundaryIdx(Point& pt1, Point& pt2, 
				double epsilon, int& bdindex,
				double& par1, double& par2, double knot_tol) const
  //--------------------------------------------------------------------------
  // Given two points on the surface boundary, find the number of the
  // corresponding boundary and the curve parameter of the closest points
  // on this surface boundary.
  //
  // Ordering of boundaries:
  //                       1
  //           ----------------------
  //           |                    |
  //         2 |                    | 3
  //      v    |                    |
  //      ^    ----------------------
  //      |-> u            0
  //===========================================================================
  {
    // Find parameter value between which the boundary curve passes.
    double u1, v1, u2, v2;
    Point cl1, cl2;
    double d1, d2;
    double tol = 1.0e-7;  // Tolerance in the parameter domain.
    closestBoundaryPoint(pt1, u1, v1, cl1, d1, tol);
    closestBoundaryPoint(pt2, u2, v2, cl2, d2, tol);

    bdindex = -1;
    if (d1 > epsilon || d2 > epsilon)
      return;        // Point not on surface

    // As we are seeking a boundary curve, we know that for both points at
    // least one parameter must be an endpoint in a u- or v-knot vector.
    // If within a basis knot, we snap the parameter.
    mesh_.knotIntervalFuzzy(XFIXED, u1, knot_tol);
    mesh_.knotIntervalFuzzy(XFIXED, u2, knot_tol);
    mesh_.knotIntervalFuzzy(YFIXED, v1, knot_tol);
    mesh_.knotIntervalFuzzy(YFIXED, v2, knot_tol);

    double startu = startparam_u();
    double endu = endparam_u();
    double startv = startparam_v();
    double endv = endparam_v();
    if (fabs(u1-u2) < fabs(v1-v2) && fabs(u1-u2) < 0.01*(endu-startu))
      {
	double umid = 0.5*(u1+u2);
	par1 = v1;
	par2 = v2;
	bdindex = (fabs(umid - startu) < fabs(endu - umid))
	  ? 2 : 3;
      }
    else if (fabs(v1-v2) < fabs(u1-u2) && fabs(v1-v2) < 0.01*(endv-startv))
      {
	double vmid = 0.5*(v1+v2);
	par1 = u1;
	par2 = u2;
	bdindex = (fabs(vmid - startv) < fabs(endv - vmid))
	  ? 0 : 1;
      }
    else
      {
	// No clear boundary is found. Is it possible to save the
	// situation? Check degeneracy. May assume there is only one degenerate edge.
	// In such a scenario we switch a parameter of the closest point.
	if (degen_.is_set_) {
	  if (degen_.b_) {
	    if (v1 < v2)
	      u1 = u2;
	    else
	      u2 = u1;
	  }
	  else if (degen_.t_) {
	    if (v1 > v2)
	      u1 = u2;
	    else
	      u2 = u1;
	  }
	  else if (degen_.l_) {
	    if (u1 < u2)
	      v1 = v2;
	    else
	      v2 = v1;
	  }
	  else if (degen_.r_) {
	    if (u1 > u2)
	      v1 = v2;
	    else
	      v2 = v1;
	  }

	  // We check whether new points are close enough.
	  Point new_pt1 = ParamSurface::point(u1, v1);
	  Point new_pt2 = ParamSurface::point(u2, v2);
	  double new_u1, new_u2, new_v1, new_v2;
	  closestBoundaryPoint(new_pt1, new_u1, new_v1, cl1, d1, tol);
	  closestBoundaryPoint(new_pt2, new_u2, new_v2, cl2, d2, tol);
	  if (d1 > epsilon || d2 > epsilon)
	    return;        // Point not on surface
	    
	  if (fabs(u1-u2) < fabs(v1-v2)) {
	    double umid = 0.5*(u1+u2);
	    par1 = v1;
	    par2 = v2;
	    bdindex = (fabs(umid - startu) < fabs(endu - umid))
	      ? 2 : 3;
	  }
	  else if (fabs(v1-v2) < fabs(u1-u2)) {
	    double vmid = 0.5*(v1+v2);
	    par1 = u1;
	    par2 = u2;
	    bdindex = (fabs(vmid - startv) < fabs(endv - vmid))
	      ? 0 : 1;
	  }
	}
	else {
	  return;
	}
	
      }
      
    return;
  }


  //===========================================================================
  void LRSplineSurface::turnOrientation()
  //===========================================================================
  {
    swapParameterDirection();
  }

  //===========================================================================
  void LRSplineSurface::swapParameterDirection()
  //===========================================================================
  {
    // We must update the mesh_, bsplines_, emap_ and domain_.

    // First the mesh.
    mesh_.swapParameterDirection();

    // We then update all basis functions in bsplines_ with the
    // reversed domain. It is only the knot indices which need
    // updating.
    // Since the key is const we must recreate the map.
    BSplineMap bsplines;
    auto iter = bsplines_.begin();
    while (iter != bsplines_.end())
      {
	unique_ptr<LRBSpline2D> bas_func; // Empty pointer.
	std::swap(bas_func, iter->second);
	bas_func->swapParameterDirection();

	// We create the new key.
	BSKey bs_key = iter->first;
	std::swap(bs_key.u_min, bs_key.v_min);
	std::swap(bs_key.u_max, bs_key.v_max);

	// We must remove the bas_func from bsplines_.
	bsplines.insert(std::make_pair(bs_key, std::move(bas_func)));
	++iter;
      }
    std::swap(bsplines_, bsplines);

    // Finally we update the elements_.
    // Since the key is const we must recreate the map.
    ElementMap emap;
    auto iter2 = emap_.begin();
    while (iter2 != emap_.end())
      {
	Element2D* elem = iter2->second.get();
	elem->swapParameterDirection();

	ElemKey elem_key = iter2->first;
	elem_key.u_min = elem->umin();
	elem_key.v_min = elem->vmin();
	
	iter2->second.release();
	emap.insert(make_pair(elem_key, std::move(unique_ptr<Element2D>(elem))));
	++iter2;
      }
    std::swap(emap_, emap);

  }

  //===========================================================================
  void LRSplineSurface::reverseParameterDirection(bool dir_is_u)
  //===========================================================================
  {
    // We must update the mesh_, bsplines_ and emap_.

    // We reverse the mesh grid (in the given direction).
    // It is important that this is performed first since we update keys
    // based on these values.
    mesh_.reverseParameterDirection(dir_is_u);
    MESSAGE("Done reversing the mesh dir!");

    // We then update all basis functions in bsplines_ with the
    // reversed domain. It is only the knot indices which need
    // updating.
    // Since the key is const we must recreate the map.
    BSplineMap bsplines;
    auto iter = bsplines_.begin();
    while (iter != bsplines_.end())
      {
	unique_ptr<LRBSpline2D> bas_func;
	std::swap(bas_func, iter->second);
	bas_func->reverseParameterDirection(dir_is_u);

	// We create the new key.
	BSKey bs_key = iter->first;
	if (dir_is_u)
	  {
	    bs_key.u_min = bas_func->umin();
	    bs_key.u_max = bas_func->umax();
	    std::swap(bs_key.u_mult1, bs_key.u_mult2);
	  }
	else
	  {
	    bs_key.v_min = bas_func->vmin();
	    bs_key.v_max = bas_func->vmax();
	    std::swap(bs_key.v_mult1, bs_key.v_mult2);
	  }

	bsplines.insert(make_pair(bs_key, std::move(bas_func)));
	++iter;
      }
    std::swap(bsplines_, bsplines);

    // Finally we update the elements_.
    // Since the key is const we must recreate the map.
    ElementMap emap;
    auto iter2 = emap_.begin();
    while (iter2 != emap_.end())
      {
	Element2D* elem = iter2->second.get();
	// We must update the end parameters of the element.
	if (dir_is_u)
	  {
	    double new_umin = mesh_.maxParam(XFIXED) - elem->umax();
	    // We must snap the updated knots to the knots in the mesh_.
	    mesh_.knotIntervalFuzzy(XFIXED, new_umin, knot_tol_);

	    double new_umax = new_umin + elem->umax() - elem->umin();
	    mesh_.knotIntervalFuzzy(XFIXED, new_umax, knot_tol_);

	    elem->setUmin(new_umin);
	    elem->setUmax(new_umax);
	  }
	else
	  {
	    double new_vmin = mesh_.maxParam(YFIXED) - elem->vmax();
	    // We must snap the updated knots to the knots in the mesh_.
	    mesh_.knotIntervalFuzzy(YFIXED, new_vmin, knot_tol_);

	    double new_vmax = new_vmin + elem->vmax() - elem->vmin();
	    mesh_.knotIntervalFuzzy(YFIXED, new_vmax, knot_tol_);

	    elem->setVmin(new_vmin);
	    elem->setVmax(new_vmax);
	  }

	ElemKey elem_key = iter2->first;
	if (dir_is_u)
	  elem_key.u_min = elem->umin();
	else
	  elem_key.v_min = elem->vmin();

	iter2->second.release();
	emap.insert(make_pair(elem_key, std::move(unique_ptr<Element2D>(elem))));
	++iter2;
      }
    std::swap(emap_, emap);
  }

  //===========================================================================
  void LRSplineSurface::setParameterDomain(double u1, double u2, double v1, double v2)
  {
    // @@sbr201301 Fix this I think ...
    //MESSAGE("I do think we should snap all knots to the mesh knots!");
    double umin = paramMin(XFIXED);
    double umax = paramMax(XFIXED);
    double vmin = paramMin(YFIXED);
    double vmax = paramMax(YFIXED);

    // All reference to knots are fetched from the mesh_.
    mesh_.setParameterDomain(u1, u2, v1, v2);

    // But in addition the ElementMap emap_ also contains keys which
    // are constructed using the parameters of the elements, as well
    // as Element2D's which stores max and min values in both dirs for
    // the elements.
    // for (ElementMap::iterator iter = elementsBeginNonconst(); iter != elementsEndNonconst(); ++iter)
    // First move the elements out of the container
    vector<unique_ptr<Element2D> > all_elements;
    for (auto it=emap_.begin(); it!=emap_.end(); ++it)
      {
	unique_ptr<Element2D> ptr = std::move(it->second);
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

	double elem_umin_new = (u2 - u1)/(umax - umin)*(elem_umin - umin) + u1;
	double elem_umax_new = (u2 - u1)/(umax - umin)*(elem_umax - umin) + u1;
	double elem_vmin_new = (v2 - v1)/(vmax - vmin)*(elem_vmin - vmin) + v1;
	double elem_vmax_new = (v2 - v1)/(vmax - vmin)*(elem_vmax - vmin) + v1;
	// We may encounter tolerance issues for the far end of the domain, snap.
	if (fabs(elem_umax_new - u2) < knot_tol_)
	  elem_umax_new = u2;
	if (fabs(elem_vmax_new - v2) < knot_tol_)
	  elem_vmax_new = v2;

	all_elements[ki]->setUmin(elem_umin_new);
	all_elements[ki]->setUmax(elem_umax_new);
	all_elements[ki]->setVmin(elem_vmin_new);
	all_elements[ki]->setVmax(elem_vmax_new);
	all_elements[ki]->updateLSDataParDomain(elem_umin, elem_umax, 
						elem_vmin, elem_vmax,
						elem_umin_new, elem_umax_new, 
						elem_vmin_new, elem_vmax_new);

	// Make new key
	ElemKey new_key;
	new_key.u_min = elem_umin_new;
	new_key.v_min = elem_vmin_new;
    
	// Insert in container
	// emap_.insert(make_pair(new_key, 
	// 		       std::move(unique_ptr<Element2D>(all_elements[ki].get()))));
	emap_.insert(make_pair(new_key, std::move(all_elements[ki])));
    }
   
    // Must also regenerate keys for the bsplines
    // First move the bsplines out of the container
    vector<unique_ptr<LRBSpline2D> > all_bsplines;
    for (auto it=bsplines_.begin(); it != bsplines_.end(); ++it)
      {
	unique_ptr<LRBSpline2D> ptr = std::move(it->second);
	all_bsplines.emplace_back(std::move(ptr));
      }

    bsplines_.clear();
    for (size_t ki=0; ki<all_bsplines.size(); ++ki)
      {
	auto key = generate_key(*all_bsplines[ki], mesh_);
	bsplines_.insert(make_pair(key, std::move(all_bsplines[ki])));
      }

    // ElementMap::iterator iter = emap_.begin();
    // size_t nmb_el = emap_.size();
    // //for (ElementMap::iterator iter = emap_.begin(); iter != emap_.end(); )
    // for (size_t ki=0; ki<nmb_el; ++ki)
    //   {
    // 	Element2D* elem = iter->second.get();

    // 	double elem_umin = elem->umin();
    // 	double elem_umax = elem->umax();
    // 	double elem_vmin = elem->vmin();
    // 	double elem_vmax = elem->vmax();

    // 	double elem_umin_new = (u2 - u1)/(umax - umin)*(elem_umin - umin) + u1;
    // 	double elem_umax_new = (u2 - u1)/(umax - umin)*(elem_umax - umin) + u1;
    // 	double elem_vmin_new = (v2 - v1)/(vmax - vmin)*(elem_vmin - vmin) + v1;
    // 	double elem_vmax_new = (v2 - v1)/(vmax - vmin)*(elem_vmax - vmin) + v1;
    // 	// We may encounter tolerance issues for the far end of the domain, snap.
    // 	if (fabs(elem_umax_new - u2) < knot_tol_)
    // 	  elem_umax_new = u2;
    // 	if (fabs(elem_vmax_new - v2) < knot_tol_)
    // 	  elem_vmax_new = v2;

    // 	elem->setUmin(elem_umin_new);
    // 	elem->setUmax(elem_umax_new);
    // 	elem->setVmin(elem_vmin_new);
    // 	elem->setVmax(elem_vmax_new);

    // 	// Since the key is const for a map element, we must replace the entry.
    // 	ElemKey new_key;
    // 	new_key.u_min = elem_umin_new;
    // 	new_key.v_min = elem_vmin_new;
    
    // auto nextIterator = iter;
    // nextIterator++;
    // iter->second.release();
    // emap_.erase(iter);
    // emap_.insert(make_pair(new_key, std::move(unique_ptr<Element2D>(elem))));
    // iter = nextIterator;
    // }
  }

  //===========================================================================
  double LRSplineSurface::area(double tol) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::area() not yet implemented.");

    // double area = 0.0;
    // int num_elem = (int)emap_.size();
    // int deg_u = degree(XFIXED);
    // int deg_v = degree(YFIXED);
    // // We make sure the error tolerance is fulfilled.
    // double elem_area_tol = tol/num_elem;
    // for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
    //   {
    // 	;//area += iter->second->surfaceArea(elem_area_tol);
    //   }

    return 0.0;
  }

  //===========================================================================
  bool LRSplineSurface::isDegenerate(bool& bottom, bool& right,
				     bool& top, bool& left, double epsilon) const
  //===========================================================================
  {
    if (degen_.is_set_ && fabs(degen_.tol_ - epsilon) < 1.0e-15)
      {
	bottom = degen_.b_;
	right = degen_.r_;
	top = degen_.t_;
	left = degen_.l_;
      }
    else
      {
	// This is not the fastest approach, a first approach.

	// We fetch the boundary curves and check if they are degenerate.
	for (int ki = 0; ki < 4; ++ki)
	  {
	    shared_ptr<SplineCurve> edge_cv(edgeCurve(ki));

	    bool degen = edge_cv->isDegenerate(epsilon);
	    switch (ki)
	      {
	      case 0:
	      {
		left = degen;
		break;
	      }
	      case 1:
	      {
		right = degen;
		break;
	      }
	      case 2:
	      {
		bottom = degen;
		break;
	      }
	      case 3:
	      {
		top = degen;
		break;
	      }
	      default:
	      {
		THROW("Should never happen ...");
		break;
	      }
	      }
	  }

	degen_.is_set_ = true;
	degen_.tol_ = epsilon;
	degen_.b_ = bottom;
	degen_.l_ = left;
	degen_.t_ = top;
	degen_.r_ = right;
      }

    return left || right || top || bottom;
  }

  //===========================================================================
  void LRSplineSurface::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
  //===========================================================================
  {
    // Parameter values in corners
    double param[8];
    param[0] = param[4] = startparam_u();
    param[2] = param[6] = endparam_u();
    param[1] = param[3] = startparam_v();
    param[5] = param[7] = endparam_v();

    // For all corners
    vector<Point> derivs(3);
    double ang;
    for (int ki=0; ki<4; ki++)
    {
	point(derivs, param[2*ki], param[2*ki+1], 1);
	ang = derivs[1].angle(derivs[2]);
	if (fabs(ang) < tol || fabs(M_PI-ang) < tol)
	    deg_corners.push_back(Point(param[2*ki], param[2*ki+1]));
    }

    return;
  }

  //===========================================================================
  void LRSplineSurface::getCornerPoints(vector<pair<Point,Point> >& corners) const
  //===========================================================================
  {
    corners.resize(4);

    // Parameter values in corners
    double param[8];
    param[0] = param[6] = startparam_u();
    param[2] = param[4] = endparam_u();
    param[1] = param[3] = startparam_v();
    param[5] = param[7] = endparam_v();

    // For all corners
    Point pos;
    for (int ki=0; ki<4; ki++)
      {
	point(pos, param[2*ki], param[2*ki+1]);
	corners[ki] = std::make_pair(pos, Point(param[2*ki], param[2*ki+1]));
      }

    return;
  }

//===========================================================================
SplineCurve*
LRSplineSurface::edgeCurve(int edge_num) const
//===========================================================================
{
  SplineCurve *edgcv = NULL;
  // Direction of constant parameter curve
  Direction2D d = (edge_num <= 1) ? XFIXED : YFIXED;  // Orthogonal to the curve
  Direction2D d2 = (edge_num <= 1) ? YFIXED : XFIXED; // Along the curve
  bool atstart = (edge_num == 0 || edge_num == 2);  // Whether the constant
  // parameter curve is in the start related to the opposite parameter direction
  // of the surface

  // Check if the basis has k-tupple knots in the current direction
  int ix = (atstart) ? mesh_.firstMeshVecIx(d) : 
    mesh_.lastMeshVecIx(d);  
  int mult = mesh_.minMultInLine(d, ix);
  int deg = degree(d);  // Degree orthogonal to the constant parameter curve
  if (mult < deg)
    {
      // The knot multiplicity is less than the order. Use functionality for
      // constant parameter curves
      return constParamCurve(mesh_.kval(d, ix), d == YFIXED);
    }

  // Fetch knot vector indices
  vector<int> knot_idx =  LRBSpline2DUtils::derive_knots(mesh_, d2, 
							 mesh_.firstMeshVecIx(d2),
							 mesh_.lastMeshVecIx(d2),
							 atstart ? ix : ix-1,
							 atstart ? ix+1 : ix);

  // Fetch associated coefficients
  // Fetch LRBSplines. Traverse the knot vector to make keys 
  // for the bspline map.
  int dir = d;
  int startmult[2], endmult[2];
  double startval[2], endval[2];
  int num = mesh_.numDistinctKnots(d);  // Number of knots orthogonal to the
  // constant parameter curve

  // Set parameter value at the constant parameter curve
  if (atstart)
    {
      startval[dir] = mesh_.kval(d, 0);
    }
  else
    {
      endval[dir] = mesh_.kval(d, num-1);
    }
  int deg2 = degree(d2);  // Degree along the constant parameter curve

#ifdef DEBUG
  vector<LRBSpline2D*> bfuncs;
  bfuncs.reserve(bsplines_.size());
  auto iter = bsplines_.begin();
  while (iter != bsplines_.end())
    {
      bfuncs.push_back(iter->second.get());
      ++iter;
    }
#endif

  // Since we have multiple knots in the surface boundary and we are only interested
  // in the basis functions being non-zero along the boundary, we know that these
  // basis functions must have a multiplicity equal to the order at the boundary
  // and one in the other end.
  startmult[dir] = atstart ? deg+1 : 1;
  endmult[dir] = atstart ? 1 : deg+1;
  size_t k1, k2;
  vector<double> coefs;
  for (k1=0, k2=deg2+1; k2<knot_idx.size(); ++k1, ++k2)
    {
      // Fetch domain of first minimal LR B-spline along the constant
      // parameter curve
      // First orthogonal to the curve
      if (atstart)
	{ // Assuming k-regularity at end boundary parameters! We will
	  // then only have to step to the neighbour knot to locate
	  // the basis function.
	  int k_idx = 
	    Mesh2DUtils::search_upwards_for_nonzero_multiplicity(mesh_, d, 1, 
								 knot_idx[k1], knot_idx[k2]);
	  endval[dir] = mesh_.kval(d, k_idx);
	}
      else
	{
	  int k_idx = 
	    Mesh2DUtils::search_downwards_for_nonzero_multiplicity(mesh_, d, num-2, 
								   knot_idx[k1], knot_idx[k2]);
	  startval[dir] = mesh_.kval(d, k_idx);
	}

      // Along the curve
      startval[1-dir] = mesh_.kval(d2, knot_idx[k1]);
      endval[1-dir] = mesh_.kval(d2, knot_idx[k2]);

      // Count multiplicitity along the curve
      int km = 1;
      for (; km<=deg2; ++km)
	if (knot_idx[k1+km] > knot_idx[k1])
	  break;
      startmult[1-dir] = km;

      km = 1;
      for (; km<=deg2; ++km)
	if (knot_idx[k2-km] < knot_idx[k2])
	  break;
      endmult[1-dir] = km;

      // Define key, i.e. the domain and multiplicities of the basis function
      BSKey key = {startval[0], startval[1], endval[0], endval[1], 
		   startmult[0], startmult[1], endmult[0], endmult[1]};
      
      // Fetch the associated LR B-spline
      const auto bm = bsplines_.find(key);
      if (bm == bsplines_.end())
       {
         knot_idx.erase(knot_idx.begin()+k1);
         --k1;
         --k2;
         continue;
       }
      //THROW("edgeCurve:: There is no such basis function.");

      // Fetch coefficient
      Point cf = bm->second->coefTimesGamma();
      coefs.insert(coefs.end(), cf.begin(), cf.end());
    }

  // Define spline curve
  // First compute the real knot values
  int nmbcf = (int)knot_idx.size() - deg2 - 1;
  vector<double> knots(knot_idx.size());
  for (k1=0; k1<knot_idx.size(); ++k1)
    knots[k1] = mesh_.kval(d2, knot_idx[k1]);

  edgcv = new SplineCurve(nmbcf, deg2+1, knots.begin(), coefs.begin(),
			  dimension(), rational_);

  return edgcv;
}

//===========================================================================
SplineCurve*
LRSplineSurface::constParamCurve(double parameter,
				 bool pardir_is_u) const
//===========================================================================
{
    // We do it the lazy (and time consuming) way by extracting the
    // subsurface and fetching the edge curve.
    int edge_num = -1; // edge_num: 0 = umin, 1 = umax, 2 = vmin, 3 = vmax
    if (!pardir_is_u && fabs(parameter - startparam_u()) < knot_tol_)
      edge_num = 0;
    else if (!pardir_is_u && fabs(parameter - endparam_u()) < knot_tol_)
      edge_num = 1;
    else if (pardir_is_u && fabs(parameter - startparam_v()) < knot_tol_)
      edge_num = 2;
    else if (pardir_is_u && fabs(parameter - endparam_v()) < knot_tol_)
      edge_num = 3;
    bool par_on_edge = (edge_num != -1);

    if (par_on_edge)
      {
	return edgeCurve(edge_num);
      }
    else
      {
	double umin = startparam_u();
	double vmin = startparam_v();
	double umax = (pardir_is_u) ? endparam_u() : parameter;
	double vmax = (!pardir_is_u) ? endparam_v() : parameter;

	shared_ptr<LRSplineSurface> sub_sf(subSurface(umin, vmin,
						      umax, vmax,
						      knot_tol_));
	int sub_edge_num = (pardir_is_u) ? 3 : 1;
	return sub_sf->edgeCurve(sub_edge_num);
      }
}


//===========================================================================
void LRSplineSurface::constParamCurve(double parameter, 
				      bool pardir_is_u, 
				      SplineCurve*& cv, 
				      SplineCurve*& crosscv) const
 //===========================================================================
{
  MESSAGE("Not yet implemented!");

}

  //===========================================================================
  vector< shared_ptr<ParamCurve> >
    LRSplineSurface::constParamCurves(double parameter, bool pardir_is_u) const
  //===========================================================================
  {
    vector<shared_ptr<ParamCurve> > return_cvs;
    return_cvs.push_back(shared_ptr<ParamCurve>(constParamCurve(parameter,
								pardir_is_u)));

    return return_cvs;
  }

  //===========================================================================
  double LRSplineSurface::nextSegmentVal(int dir, double par, bool forward, double tol) const
  //===========================================================================
  {
    // The given interface means we must look at the global knot vectors,
    // i.e. not the parameter corresponding to the next element.
    double startpar = (dir == 0) ? startparam_u() : startparam_v();
    double endpar = (dir == 0) ? endparam_u() : endparam_v();

    if (!forward && par <= startpar)
      return startpar;
  
    if (forward && par >= endpar)
      return endpar;

    // We locate the index for which et[ind] <= par < et[ind+1].
    Direction2D dir2d = (dir == 0) ? XFIXED : YFIXED;
    vector<double> knots(mesh_.knotsBegin(dir2d), mesh_.knotsEnd(dir2d));

    std::vector<double>::const_iterator knot;
    if (forward)
      {
	par += fabs(tol);
	knot = std::upper_bound(knots.begin(),knots.end(),par);
	if (knot == knots.end())
	  return endpar;
	else
	  return *knot;
      }
    else
      {
	par -= fabs(tol);
	for (knot=knots.end()-1; knot>knots.begin(); --knot)
	  {
	    if (*knot < par)
	      return *knot;
	  }
	return *(knots.begin());
      }

    return 0.0;
  }

//===========================================================================
LRSplineSurface* LRSplineSurface::mirrorSurface(const Point& pos, 
						const Point& norm) const
//===========================================================================
  {
    LRSplineSurface* mirrored = clone();

    Point normal = norm;
    normal.normalize();

    for (auto iter = mirrored->basisFunctionsBegin(); iter != mirrored->basisFunctionsEnd(); ++iter)
      {
	LRBSpline2D* bas_func = iter->second.get();
	Point coef = bas_func->Coef();
	double gamma = bas_func->gamma();

	Point tmp = ((coef - pos)*normal) * normal;
	coef -= 2.0*tmp;

	bas_func->setCoefAndGamma(coef, gamma);
      }

    return mirrored;
  }

//===========================================================================
vector<LRBSpline2D*> 
LRSplineSurface::collect_basis(int from_u, int to_u, 
			       int from_v, int to_v) const
//===========================================================================
{
  // vector<unique_ptr<LRBSpline2D> > b_splines;
  vector<LRBSpline2D*> b_splines;
  for (auto it=bsplines_.begin(); it!=bsplines_.end(); ++it)
    {
      if (interval_overlap(it->second->suppMin(XFIXED), 
			   it->second->suppMax(XFIXED), 
			   from_u, to_u) &&
	  interval_overlap(it->second->suppMin(YFIXED), 
			   it->second->suppMax(YFIXED), 
			   from_v, to_v))
	{
	  // unique_ptr<LRBSpline2D> ptr = std::move(it->second);
	  // b_splines.emplace_back(std::move(ptr));
	  b_splines.push_back(it->second.get());
	}
    }
  return b_splines;
}


//===========================================================================
void
LRSplineSurface::checkSupport(LRBSpline2D* basis) const
//===========================================================================
{
  // Domain
  int u1_ix = basis->suppMin(XFIXED);
  int u2_ix = basis->suppMax(XFIXED);
  int v1_ix = basis->suppMin(YFIXED);
  int v2_ix = basis->suppMax(YFIXED);

  // Collect elements
  std::set<Element2D*> elem1;
  double v1 = mesh_.kval(YFIXED, v1_ix);
  double v2, u1, u2;
  for (int kj=v1_ix+1; kj<=v2_ix; ++kj)
    {
      v2 =  mesh_.kval(YFIXED, kj);
      double vpar = 0.5*(v1 + v2);
      u1 = mesh_.kval(XFIXED, u1_ix);
      for (int ki=u1_ix+1; ki<=u2_ix; ++ki)
	{
	  u2 =  mesh_.kval(XFIXED, ki);
	  double upar = 0.5*(u1 + u2);
	  Element2D *el = coveringElement(upar, vpar);
	  elem1.insert(elem1.end(), el);
	  u1 = u2;
	}
      v1 = v2;
    }

  vector<Element2D*> elem2(elem1.begin(), elem1.end());

  // Compare with the elements in the bspline support
  for (size_t kr=0; kr<elem2.size(); ++kr)
    {
      int kn=0;
      for (auto iter = basis->supportedElementBegin(); 
	   iter != basis->supportedElementEnd(); ++iter)
	{
	  if (elem2[kr] == (*iter))
	    break;
	  ++kn;
	}
      if (kn == basis->nmbSupportedElements())
	{
	  std::cout << "Bspline " << basis << " misses element " << elem2[kr];
	  std::cout << " in its support" << std::endl;

	  if (!elem2[kr]->hasSupportFunction(basis))
	    std::cout << "Element " << elem2[kr] << " misses basis " << basis << std::endl;
	}
    }
  
	  
}

void 
LRSplineSurface::s1773(const double ppoint[],double aepsge, 
		       double estart[],double eend[],double enext[],
		       double gpos[], int maxiter,
		       Element2D* elem, int *jstat) const
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a surface and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameters is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : ppoint   - The point in the closest point problem.
*              psurf    - The surface in the closest point problem.
*              aepsge   - Geometry resolution.
*              estart   - Surface parameters giving the start of the search
*                         area (umin, vmin).
*              eend     - Surface parameters giving the end of the search
*                         area (umax, vmax).
*              enext    - Surface guess parameters for the closest point
*                         iteration.
*
*
*
* OUTPUT     : gpos    - Resulting surface parameters from the iteration.
*              jstat   - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in two parameter direction.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, May 1989
* Revised by : Johannes Kaasa, SINTEF Oslo, August 1995.
*              Introduced a local copy of enext, to avoid changes.
*
*********************************************************************
*/                       
{                        
  // int kstat = 0;            /* Local status variable.                      */
  int kder=1;               /* Order of derivatives to be calulated        */
  int kdim;                 /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  int kdir;                 /* Changing direction.                         */
  int kdeg;                 /* Degenaracy flag.                            */
  double tdelta[2];         /* Parameter intervals of the surface.         */
  double tdist;             /* Distance between position and origo.        */
  double td[2],t1[2],tdn[2];/* Distances between old and new parameter
			       value in the tree parameter directions.     */
  double tprev;             /* Previous difference between the curves.     */
  vector<double> sdiff;      /* Difference between the point and the surf.  */
  double snext[2];          /* Parameter values                            */
  double guess[2];          /* Local copy of enext.                        */
  double REL_COMP_RES = 1.0e-12; //1.0e-15;
  vector<Point> pts(3);
  double fac = 1.5;
  
  guess[0] = enext[0];
  guess[1] = enext[1];
  
  kdim = dimension();
  
  
  /* Fetch endpoints and the intervals of parameter interval of curves.  */
  
  RectDomain dom = containingDomain();
  tdelta[0] = std::max(dom.umin(), dom.umax() - dom.umin());
  tdelta[1] = std::max(dom.vmin(), dom.vmax() - dom.vmin());
  
  /* Allocate local used memory */
  sdiff.resize(kdim);
  
  /* Initiate variables.  */
  
  tprev = 1.0e10;
  
  /* Evaluate 0-1.st derivatives of surface */
  /* printf("\n lin: \n %#20.20g %#20.20g",
     guess[0],guess[1]); */
  
  point(pts, guess[0], guess[1], kder, elem);
  
  /* Compute the distanse vector and value and the new step. */
  
  s1773_s9dir(&tdist,td,td+1,&sdiff[0],ppoint,pts,
	      aepsge,kdim,&kdeg);
  
  /* Correct if we are not inside the parameter intervall. */
  
  t1[0] = td[0];
  t1[1] = td[1];
  s1773_s9corr(t1,guess[0],guess[1],estart[0],eend[0],estart[1],eend[1]);
  
  /* Iterate to find the intersection point.  */
  
  tprev = tdist;
  for (knbit = 0; knbit < maxiter; knbit++)
    {
      /* Evaluate 0-1.st derivatives of surface */
      
      snext[0] = guess[0] + t1[0];
      snext[1] = guess[1] + t1[1];
      
      point(pts, snext[0], snext[1], kder, elem);
      
      /* Compute the distanse vector and value and the new step. */
      
      s1773_s9dir(&tdist,tdn,tdn+1,&sdiff[0],ppoint,
	    pts,aepsge,kdim,&kdeg);
      
      /* Check if the direction of the step have change. */
      
      kdir = (Utils::inner(td, td+2, tdn) >= 0.0);     /* 0 if changed. */
      
      /* Ordinary converging. */
      
      if (tdist < tprev/(double)2 || (kdir && tdist < fac*tprev))
	{
	   guess[0] += t1[0];
	   guess[1] += t1[1];
  
	  /* printf("\n %#20.20g %#20.20g",
	     guess[0],guess[1]); */
  
	  
          td[0] = t1[0] = tdn[0];
          td[1] = t1[1] = tdn[1];
	  
	  /* Correct if we are not inside the parameter intervall. */
	  
	  s1773_s9corr(t1,guess[0],guess[1],estart[0],eend[0],estart[1],eend[1]);
          tprev = tdist;

	  if ( (fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES)) break;
	}
      
      /* Not converging, adjust and try again. */
      
      else
	{
          t1[0] /= (double)2;
          t1[1] /= (double)2;
          /* knbit--;  */
	}
      if (guess[0]==guess[0]+t1[0] &&
	  guess[1]==guess[1]+t1[1]) break;
    }
  
  /* Iteration stopped, test if point founds found is within resolution */
  
  if (tdist <= aepsge)
  {
     *jstat = 1;
     /* printf("\n SUCCESS!!"); */
     
  }
  else if(kdeg)
     *jstat = 9;
  else
     *jstat = 2;
  
  gpos[0] = guess[0];
  gpos[1] = guess[1];
  

  return;

}



 } // end namespace Go


