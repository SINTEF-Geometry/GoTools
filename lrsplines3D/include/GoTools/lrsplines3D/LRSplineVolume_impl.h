//===========================================================================
//                                                                           
// File: LRSplineVolume_impl.h                                               
//                                                                           
// Created: Mon Feb 25 11:29:39 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LRSPLINEVOLUME_IMPL_H
#define _LRSPLINEVOLUME_IMPL_H


#include <algorithm>
#include <iostream>    // @@ Only used for debug.  
#include "GoTools/lrsplines3D/Mesh3DUtils.h" // @@ Only used for debug.  
#include "GoTools/lrsplines3D/Mesh3D.h" // @@ Only used for debug.  


namespace Go 
{

namespace {
  // double least   (double a, double b)           { return std::min(a, b);}
  // double least   (double a, double b, double c) { return std::min(std::min(a, b), c); }
  // double greatest(double a, double b)           { return std::max(a, b);}
  // double greatest(double a, double b, double c) { return std::max(std::max(a, b), c); }

  // making a 'knot vector' (with explicit multiplicities), whose values are indices
  // into the knotvector values stored in the mesh, rather than storing the knot values
  // directly.
  std::vector<int> init_knot_indices(const Go::Mesh3D& m, Go::Direction3D d)
  {    
    std::vector<int> result;
    for (int i = 0; i != m.numDistinctKnots(d); ++i)
	result.insert(result.end(), m.nu(d, i, 0, 1, 0, 1), i);
    return result;
  }

  // /* template<typename IntIterator> */
  // /* int consecutives(IntIterator start, IntIterator end) { */
  // /*   return std::find_if(start, end, [start](int x) {return x != *start;}) - start; */

  // int consecutives(const std::vector<int>& v) {
  //   return std::find_if(v.begin(), v.end(), [&](int x) 
  // 			{return x != v.front();}) - v.begin();}

  // int predessesors(const std::vector<int>& v) {
  //   return std::find_if(v.rbegin(), v.rend(), [&](int x) 
  // 			{return x != v.back();}) - v.rbegin();}

}; // end anonymous namespace


//==============================================================================
inline int LRSplineVolume::degree(Direction3D d) const
//==============================================================================
{
  // returns zero degree if spline is invalid (uninitialized).
  return bsplines_.size() > 0 ? bsplines_.begin()->second->degree(d) : 0;
}

//==============================================================================
inline double LRSplineVolume::paramMin(Direction3D d) const
//==============================================================================
{
  // this routine will give the expected result when the mesh was initialized 
  // with a fixed number of knot multiplicities at start.  (In all normal
  // situations, including when using k-multiple knots, this is the case).
  // However, a global start multiplicity cannot be defined if the multiplicities
  // of meshrectangles vary along the border.  
  // In this case, the minimum multiplicity found would be used, as determined
  // by the nu-operator.
  const int mul = mesh_.nu(d, 0,
			   0, mesh_.numDistinctKnots(next(d)) - 1,
			   0, mesh_.numDistinctKnots(prev(d)) - 1);
  const int deg = degree(d);
  if (mul < deg + 1)
  { // Not k-regular.
      // We must find the knot with index 'deg' (couting multiplicities).
      int uniq_ind = 0;
      int tot_ind = mul - 1;
      while (tot_ind < deg)
      {
	  ++uniq_ind;
	  const int mul2 = mesh_.nu(d, uniq_ind,
				    0, mesh_.numDistinctKnots(next(d)) - 1,
				    0, mesh_.numDistinctKnots(prev(d)) - 1);
	  tot_ind += mul2;
      }
      return mesh_.kval(d, uniq_ind);
  }

  return mesh_.kval(d, std::max(deg + 1 - mul, 0));
}

//==============================================================================
inline double LRSplineVolume::paramMax(Direction3D d) const
//==============================================================================
{
  // this routine will give the expected result when the mesh was initialized 
  // with a fixed number of knot multiplicities at end.  (In all normal
  // situations, including when using k-multiple knots, this is the case).
  // However, a global end multiplicity cannot be defined if the multiplicities
  // of meshrectangles vary along the border.  
  // In this case, the minimum multiplicity found would be used, as determined
  // by the nu-operator.
  const int mul = mesh_.nu(d, mesh_.numDistinctKnots(d) - 1,
			   0, mesh_.numDistinctKnots(next(d)) - 1,
			   0, mesh_.numDistinctKnots(prev(d)) - 1);
  const int deg = degree(d);
  if (mul < deg + 1)
  { // Not k-regular.
      // We must find the knot with index 'deg' (couting multiplicities).
      int uniq_ind = mesh_.numDistinctKnots(d) - 1;
      int tot_ind = mul - 1;
      while (tot_ind < deg)
      {
	  --uniq_ind;
	  const int mul2 = mesh_.nu(d, uniq_ind,
				    0, mesh_.numDistinctKnots(next(d)) - 1,
				    0, mesh_.numDistinctKnots(prev(d)) - 1);
	  tot_ind += mul2;
      }
      return mesh_.kval(d, uniq_ind);
  }

  return mesh_.kval(d, mesh_.numDistinctKnots(d) - 1 - std::max(degree(d) + 1 - mul, 0));
}



// =============================================================================
inline LRSplineVolume::BSKey LRSplineVolume::generate_key(const LRBSpline3D& b, 
							    const Mesh3D& m)
// =============================================================================
{
  /* BSKey key = { m.kval(XFIXED, b.suppMin(XFIXED)), */
  /* 		m.kval(YFIXED, b.suppMin(YFIXED)), */
  /* 		m.kval(XFIXED, b.suppMax(XFIXED)), */
  /* 		m.kval(YFIXED, b.suppMax(YFIXED)), */
  /* 		consecutives(b.kvec(XFIXED).begin(), b.kvec(XFIXED).end()), */
  /* 		consecutives(b.kvec(YFIXED).begin(), b.kvec(YFIXED).end())}; */
    BSKey key = { m.kval(XDIR, b.suppMin(XDIR)),
		  m.kval(YDIR, b.suppMin(YDIR)),
		  m.kval(ZDIR, b.suppMin(ZDIR)),
		  m.kval(XDIR, b.suppMax(XDIR)),
		  m.kval(YDIR, b.suppMax(YDIR)),
		  m.kval(ZDIR, b.suppMax(ZDIR)),
		  consecutives(b.kvec(XDIR)), consecutives(b.kvec(YDIR)), consecutives(b.kvec(ZDIR)),
		  predessesors(b.kvec(XDIR)), predessesors(b.kvec(YDIR)), predessesors(b.kvec(ZDIR))};
#ifndef NDEBUG
  double deb_val = 0.0;
#endif
   return key;
}

// =============================================================================
inline LRSplineVolume::BSKey 
  LRSplineVolume::generate_key(const LRBSpline3D& b)
// =============================================================================
{
  /* BSKey key = { b.umin(), b.vmin(), b.umax(), b.vmax(), */
  /* 		consecutives(b.kvec(XFIXED).begin(), b.kvec(XFIXED).end()), */
  /* 		consecutives(b.kvec(YFIXED).begin(), b.kvec(YFIXED).end())}; */
  BSKey key = { b.umin(), b.vmin(), b.wmin(), b.umax(), b.vmax(), b.wmax(),
		consecutives(b.kvec(XDIR)), consecutives(b.kvec(YDIR)), consecutives(b.kvec(ZDIR)),
		predessesors(b.kvec(XDIR)), predessesors(b.kvec(YDIR)), predessesors(b.kvec(ZDIR))};
   return key;
}

// =============================================================================
inline bool LRSplineVolume::BSKey::operator<(const BSKey& rhs) const
// =============================================================================
{
  // The following way to define the '<' operator should ensure that when the
  // LRSplineVolume is expanded into full tensor-product representation, the basis
  // functions would be ordered in the usual way (stored sequentially first along
  // the u-direction (lowest stride), then along the v-direction
  return
    (w_min  < rhs.w_min ) ? true : 
    (w_min  > rhs.w_min ) ? false :
    (w_mult1 > rhs.w_mult1) ? true :
    (w_mult1 < rhs.w_mult1) ? false :
    (w_mult2 < rhs.w_mult2) ? true :
    (w_mult2 > rhs.w_mult2) ? false :
    (v_min  < rhs.v_min ) ? true : 
    (v_min  > rhs.v_min ) ? false :
    (v_mult1 > rhs.v_mult1) ? true :
    (v_mult1 < rhs.v_mult1) ? false :
    (v_mult2 < rhs.v_mult2) ? true :
    (v_mult2 > rhs.v_mult2) ? false :
    (u_min  < rhs.u_min ) ? true : 
    (u_min  > rhs.u_min ) ? false :
    (u_mult1 > rhs.u_mult1) ? true :
    (u_mult1 < rhs.u_mult1) ? false :
    (u_mult2 < rhs.u_mult2) ? true :
    (u_mult2 > rhs.u_mult2) ? false :

    (w_max  < rhs.w_max)  ? true :
    (w_max  > rhs.w_max)  ? false :
    (v_max  < rhs.v_max)  ? true :
    (v_max  > rhs.v_max)  ? false :
    (u_max  < rhs.u_max)  ? true :
    (u_max  > rhs.u_max)  ? false :
                            false; // _all_ members exactly equal
}

// =============================================================================
inline bool LRSplineVolume::ElemKey::operator<(const ElemKey& rhs) const
// =============================================================================
{
  // The following way to define the '<' operator should ensure that when the
  // LRSplineSurface is expanded into full tensor-product representation, the basis
  // functions would be ordered in the usual way (stored sequentially first along
  // the u-direction (lowest stride), then along the v-direction
  return
    (w_min  < rhs.w_min ) ? true : 
    (w_min  > rhs.w_min ) ? false :
    (v_min  < rhs.v_min ) ? true : 
    (v_min  > rhs.v_min ) ? false :
    (u_min  < rhs.u_min ) ? true : 
    (u_min  > rhs.u_min ) ? false :
                            false; // _all_ members exactly equal
}


// =============================================================================
template<typename KnotIterator, typename CoefIterator>
LRSplineVolume::LRSplineVolume(int deg_u,
			       int deg_v,
			       int deg_w,
			       int coefs_u,
			       int coefs_v,
			       int coefs_w,
			       int dimension,
			       KnotIterator knotvals_u_start,
			       KnotIterator knotvals_v_start,
			       KnotIterator knotvals_w_start,
			       CoefIterator coefs_start,
			       double knot_tol)
// =============================================================================
  : knot_tol_(knot_tol), rational_(false), curr_element_(NULL),
    mesh_(knotvals_u_start, knotvals_u_start + coefs_u + deg_u + 1,
	  knotvals_v_start, knotvals_v_start + coefs_v + deg_v + 1,
	  knotvals_w_start, knotvals_w_start + coefs_w + deg_w + 1)
{
  std::vector<int> knot_ixs_u = init_knot_indices(mesh_, XDIR);
  std::vector<int> knot_ixs_v = init_knot_indices(mesh_, YDIR);
  std::vector<int> knot_ixs_w = init_knot_indices(mesh_, ZDIR);

  bool rat = rational_;
  double rat_val = 1.0;

  // Store uni-variate B-splines
  bsplinesuni1_.resize(coefs_u);
  bsplinesuni2_.resize(coefs_v);
  bsplinesuni3_.resize(coefs_w);
  for (int u_ix = 0; u_ix != coefs_u; ++u_ix) {
    bsplinesuni1_[u_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(1, deg_u,
							       knot_ixs_u.begin() + u_ix,
							       &mesh_)));
    bsplinesuni1_[u_ix]->incrCount();
  }
  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    bsplinesuni2_[v_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(2, deg_v,
							       knot_ixs_v.begin() + v_ix,
							       &mesh_)));
    bsplinesuni2_[v_ix]->incrCount();
  }

  for (int w_ix = 0; w_ix != coefs_w; ++w_ix)  {
    bsplinesuni3_[w_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(3, deg_w,
							       knot_ixs_w.begin() + w_ix,
							       &mesh_)));
    bsplinesuni3_[w_ix]->incrCount();
  }

  for (int w_ix = 0; w_ix != coefs_w; ++w_ix)  {
      for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
	  for (int u_ix = 0; u_ix != coefs_u; ++u_ix, coefs_start += dimension) {
	    BSplineUniLR *tmp1 = bsplinesuni1_[u_ix].get();
	    BSplineUniLR *tmp2 = bsplinesuni2_[v_ix].get();
	    BSplineUniLR *tmp3 = bsplinesuni3_[w_ix].get();
	      std::unique_ptr<LRBSpline3D> b(new LRBSpline3D(Point(coefs_start, 
								   coefs_start + dimension),
							     rat_val,
							     tmp1, tmp2, tmp3,
							     1.0, rat));
	      LRSplineVolume::BSKey bs_key = generate_key(*b, mesh_);
	      //bsplines_[generate_key(*b, mesh_)] = b;
	      bsplines_.insert(std::move(std::pair<LRSplineVolume::BSKey, std::unique_ptr<LRBSpline3D> >(bs_key, std::move(b))));
	  }
      }
  }

  // Identifying all elements and mapping the basis functions to them
  emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
template<typename KnotIterator>
LRSplineVolume::LRSplineVolume(int deg_u,
			       int deg_v,
			       int deg_w,
			       int coefs_u,
			       int coefs_v,
			       int coefs_w,
			       int dimension,
			       KnotIterator knotvals_u_start,
			       KnotIterator knotvals_v_start,
			       KnotIterator knotvals_w_start,
			       double knot_tol)
//==============================================================================
  : knot_tol_(knot_tol), rational_(false), curr_element_(NULL),
    mesh_(knotvals_u_start, knotvals_u_start + coefs_u + deg_u + 1,
	  knotvals_v_start, knotvals_v_start + coefs_v + deg_v + 1,
	  knotvals_w_start, knotvals_w_start + coefs_w + deg_w + 1)
{
  std::vector<int> knot_ixs_u = init_knot_indices(mesh_, XDIR);
  std::vector<int> knot_ixs_v = init_knot_indices(mesh_, YDIR);
  std::vector<int> knot_ixs_w = init_knot_indices(mesh_, ZDIR);

  const double rat_val = 1.0;
  bool rat = rational_;
  Point p_zero(dimension);
  p_zero.setValue(0.0);

  // Store uni-variate B-splines
  bsplinesuni1_.resize(coefs_u);
  bsplinesuni2_.resize(coefs_v);
  bsplinesuni3_.resize(coefs_w);
  for (int u_ix = 0; u_ix != coefs_u; ++u_ix) {
    bsplinesuni1_[u_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(1, deg_u,
							       knot_ixs_u.begin() + u_ix,
							       &mesh_)));
    bsplinesuni1_[u_ix]->incrCount();
  }
  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    bsplinesuni2_[v_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(2, deg_v,
							       knot_ixs_v.begin() + v_ix,
							       &mesh_)));
    bsplinesuni2_[v_ix]->incrCount();
  }

  for (int w_ix = 0; w_ix != coefs_w; ++w_ix)  {
    bsplinesuni3_[w_ix] = 
      std::move(std::unique_ptr<BSplineUniLR>(new BSplineUniLR(3, deg_w,
							       knot_ixs_w.begin() + w_ix,
							       &mesh_)));
    bsplinesuni3_[w_ix]->incrCount();
  }

  for (int w_ix = 0; w_ix != coefs_w; ++w_ix)  {
    for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
      for (int u_ix = 0; u_ix != coefs_u; ++u_ix) {
	    BSplineUniLR *tmp1 = bsplinesuni1_[u_ix].get();
	    BSplineUniLR *tmp2 = bsplinesuni2_[v_ix].get();
	    BSplineUniLR *tmp3 = bsplinesuni3_[w_ix].get();
	    std::unique_ptr<LRBSpline3D> b(new LRBSpline3D(p_zero,
							   rat_val,
							   tmp1, tmp2, tmp3,
							   1.0, rat));
        LRSplineVolume::BSKey bs_key = generate_key(*b, mesh_);
        bsplines_.insert(std::move(std::pair<LRSplineVolume::BSKey, std::unique_ptr<LRBSpline3D> >(bs_key, std::move(b))));
        //	bsplines_[generate_key(*b, mesh_)] = b;
      }
    }
  }

  emap_ = construct_element_map_(mesh_, bsplines_);
}

}; // end namespace Go

#endif // _LRSPLINEVOLUME_IMPL_H

