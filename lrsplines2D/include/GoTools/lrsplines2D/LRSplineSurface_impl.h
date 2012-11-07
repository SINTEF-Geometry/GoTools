#ifndef LR_SPLINESURFACE_IMPL_H
#define LR_SPLINESURFACE_IMPL_H

#include <algorithm>
#include <iostream>    // @@ Only used for debug.  
#include "GoTools/lrsplines2D/Mesh2DUtils.h" // @@ Only used for debug.  

namespace {
  double least   (double a, double b)           { return std::min(a, b);}
  double least   (double a, double b, double c) { return std::min(std::min(a, b), c); }
  double greatest(double a, double b)           { return std::max(a, b);}
  double greatest(double a, double b, double c) { return std::max(std::max(a, b), c); }

  // making a 'knot vector' (with explicit multiplicities), whose values are indices
  // into the knotvector values stored in the mesh, rather than storing the knot values
  // directly.
  std::vector<int> init_knot_indices(const Go::Mesh2D& m, Go::Direction2D d)  
  {    
    std::vector<int> result;
    for (int i = 0; i != m.numDistinctKnots(d); ++i)
      result.insert(result.end(), m.nu(d, i, 0, 1), i);
    return result;
  }

  template<typename IntIterator>
  int consecutives(IntIterator start, IntIterator end) {
    return std::find_if(start, end, [start](int x) {return x != *start;}) - start;
  }

}; // end anonymous namespace

namespace Go 
{

// =============================================================================
inline int LRSplineSurface::degree(Direction2D d) const
// =============================================================================
{
  // returns zero degree if spline is invalid (uninitialized).
  return bsplines_.size() > 0 ? bsplines_.begin()->second.degree(d) : 0;
}

//==============================================================================
inline double LRSplineSurface::paramMin(Direction2D d) const
//==============================================================================
{
  // this routine will give the expected result when the mesh was initialized 
  // with a fixed number of knot multiplicities at start.  (In all normal
  // situations, including when using k-multiple knots, this is the case).
  // However, a global start multiplicity cannot be defined if the multiplicities
  // of meshrectangles vary along the border.  
  // In this case, the minimum multiplicity found would be used, as determined
  // by the nu-operator.
  const int mul = mesh_.nu(d, 0, 0, mesh_.numDistinctKnots(flip(d)) - 1);
  return mesh_.kval(d, std::max(degree(d) + 1 - mul, 0));
}

//==============================================================================
inline double LRSplineSurface::paramMax(Direction2D d) const
//==============================================================================
{
  // this routine will give the expected result when the mesh was initialized 
  // with a fixed number of knot multiplicities at end.  (In all normal
  // situations, including when using k-multiple knots, this is the case).
  // However, a global end multiplicity cannot be defined if the multiplicities
  // of meshrectangles vary along the border.  
  // In this case, the minimum multiplicity found would be used, as determined
  // by the nu-operator.
  const int mul = mesh_.nu(d, mesh_.numDistinctKnots(d) - 1, 0, mesh_.numDistinctKnots(flip(d)) - 1);
  return mesh_.kval(d, mesh_.numDistinctKnots(d) - 1 - std::max(degree(d) + 1 - mul, 0));
}

// =============================================================================
inline LRSplineSurface::BSKey LRSplineSurface::generate_key(const LRBSpline2D& b, 
							    const Mesh2D& m)
// =============================================================================
{
  BSKey key = { m.kval(XFIXED, b.suppMin(XFIXED)),
		m.kval(YFIXED, b.suppMin(YFIXED)),
		m.kval(XFIXED, b.suppMax(XFIXED)),
		m.kval(YFIXED, b.suppMax(YFIXED)),
		consecutives(b.kvec(XFIXED).begin(), b.kvec(XFIXED).end()),
		consecutives(b.kvec(YFIXED).begin(), b.kvec(YFIXED).end())};
  return key;
}

// =============================================================================
inline bool LRSplineSurface::BSKey::operator<(const BSKey& rhs) const
// =============================================================================
{
  // The following way to define the '<' operator should ensure that when the
  // LRSplineSurface is expanded into full tensor-product representation, the basis
  // functions would be ordered in the usual way (stored sequentially first along
  // the u-direction (lowest stride), then along the v-direction
  return
    (v_min  < rhs.v_min ) ? true : 
    (v_min  > rhs.v_min ) ? false :
    (v_mult > rhs.v_mult) ? true :
    (v_mult < rhs.v_mult) ? false :
    (u_min  < rhs.u_min ) ? true : 
    (u_min  > rhs.u_min ) ? false :
    (u_mult > rhs.u_mult) ? true :
    (u_mult < rhs.u_mult) ? false :
    (v_max  < rhs.v_max)  ? true :
    (v_max  > rhs.v_max)  ? false :
    (u_max  < rhs.u_max)  ? true :
    (u_max  > rhs.u_max)  ? false :
                            false; // _all_ members exactly equal
}

// =============================================================================
inline LRSplineSurface::ElemKey 
  LRSplineSurface::generate_key(const double& umin, 
				const double& vmin)
// =============================================================================
{
  ElemKey key = { umin, vmin};
  return key;
}

// =============================================================================
inline bool LRSplineSurface::ElemKey::operator<(const ElemKey& rhs) const
// =============================================================================
{
  // The following way to define the '<' operator should ensure that when the
  // LRSplineSurface is expanded into full tensor-product representation, the basis
  // functions would be ordered in the usual way (stored sequentially first along
  // the u-direction (lowest stride), then along the v-direction
  return
    (v_min  < rhs.v_min ) ? true : 
    (v_min  > rhs.v_min ) ? false :
    (u_min  < rhs.u_min ) ? true : 
    (u_min  > rhs.u_min ) ? false :
                            false; // _all_ members exactly equal
}

// =============================================================================
template<typename KnotIterator, typename CoefIterator>
LRSplineSurface::LRSplineSurface(int deg_u,
				 int deg_v,
				 int coefs_u,
				 int coefs_v,
				 int dimension,
				 KnotIterator knotvals_u_start,
				 KnotIterator knotvals_v_start,
				 CoefIterator coefs_start,
				 double knot_tol)
// =============================================================================
  : knot_tol_(knot_tol), rational_(false),
    mesh_(knotvals_u_start, knotvals_u_start + coefs_u + deg_u + 1,
	  knotvals_v_start, knotvals_v_start + coefs_v + deg_v + 1)
{
  std::vector<int> knot_ixs_u = init_knot_indices(mesh_, XFIXED);
  std::vector<int> knot_ixs_v = init_knot_indices(mesh_, YFIXED);

  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    for (int u_ix = 0; u_ix != coefs_u; ++u_ix, coefs_start += dimension) {
      LRBSpline2D b(Point(coefs_start, coefs_start + dimension),
		    deg_u,
		    deg_v,
		    knot_ixs_u.begin() + u_ix,
		    knot_ixs_v.begin() + v_ix,
		    1.0, &mesh_);
      bsplines_[generate_key(b, mesh_)] = b;
    }
  }
  // Identifying all elements and mapping the basis functions to them
  emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
template<typename KnotIterator>
LRSplineSurface::LRSplineSurface(int deg_u,
				 int deg_v,
				 int coefs_u,
				 int coefs_v,
				 int dimension,
				 KnotIterator knotvals_u_start,
				 KnotIterator knotvals_v_start,
				 double knot_tol)
//==============================================================================
: knot_tol_(knot_tol), rational_(false),
    mesh_(knotvals_u_start, knotvals_u_start + coefs_u + deg_u + 1,
	  knotvals_v_start, knotvals_v_start + coefs_v + deg_v + 1)
{
  std::vector<int> knot_ixs_u = init_knot_indices(mesh_, XFIXED);
  std::vector<int> knot_ixs_v = init_knot_indices(mesh_, YFIXED);

  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    for (int u_ix = 0; u_ix != coefs_u; ++u_ix) {
      LRBSpline2D b(Point(dimension),
		    deg_u,
		    deg_v,
		    knot_ixs_u.begin() + u_ix,
		    knot_ixs_v.begin() + v_ix,
		    1.0, &mesh_);
      bsplines_[generate_key(b, mesh_)] = b;
    }
  }
  emap_ = construct_element_map_(mesh_, bsplines_);
}

}; // end namespace Go

#endif
