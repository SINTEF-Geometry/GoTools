#include <assert.h>
#include <stdexcept>
//#include <iostream> // @@ debug only
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"


using namespace std;

namespace Go
{
namespace {

// Derive the knotvector resulting from moving along the u-direction (if d == XFIXED) or 
// v-direction (if d == YFIXED) from 'beg' to 'end', and  with multiplicities equal to 
// the 'nu' value of the segment in the ortogonal direction, starting at 'orto_min' and 
// extending to 'orto_max'.
//------------------------------------------------------------------------------
vector<int> 
derive_knots(const Mesh2D& m, Direction2D d, int beg, int end, int orto_min, int orto_max)
//------------------------------------------------------------------------------
{
  vector<int> result;
  for (int pos = beg; pos <= end; ++pos) 
    result.insert(result.end(), m.nu(d, pos, orto_min, orto_max), pos); // 0 multiplicities will not be inserted

  return result;
}
  
// return the multiplicity of the first knot in a knotvec
//------------------------------------------------------------------------------
inline int start_multiplicity(const vector<int>& v) 
//------------------------------------------------------------------------------
{
  assert(!v.empty());
  return find_if(v.begin(), v.end(), [&](int x) {return x != v.front();}) - v.begin();
}

// return the multiplicithy of the last knot in the knotvec
//------------------------------------------------------------------------------
inline int end_multiplicity(const vector<int>& v)
//------------------------------------------------------------------------------
{
  assert(!v.empty());
  return find_if(v.rbegin(), v.rend(), [&](int x) {return x != v.back();}) - v.rbegin();
}

// return the number of inner knots (counting multiplicities) in the knotvec
//------------------------------------------------------------------------------
int num_inner_knots(const vector<int>& v)
//------------------------------------------------------------------------------
{
  return (int)v.size() - start_multiplicity(v) - end_multiplicity(v);
}

// Find and return the index of the first knot in 'kvec1' that is not found in
// 'kvec2'.  There should be at least one such knot - otherwise an exception will
// be thrown.
//------------------------------------------------------------------------------
int find_uncovered_inner_knot(const vector<int>& kvec1, const vector<int>& kvec2)
//------------------------------------------------------------------------------
{
  const int m1 = start_multiplicity(kvec1);
  const int m2 = start_multiplicity(kvec2);
  const auto tmp = std::mismatch(kvec1.begin() + (m1 - 1), kvec1.end(), kvec2.begin() + (m2 - 1));

  if (tmp.first == kvec1.end()) throw runtime_error("No uncovered inner knot found.");

  return *tmp.first;
}
  

}; // end anonymous namespace


//==============================================================================
  void LRBSpline2DUtils::split_function(const LRBSpline2D& orig, 
					Direction2D d, 
					const double* const kvals,
					int new_knot_ix,
					LRBSpline2D& new_1,
					LRBSpline2D& new_2)
//==============================================================================
{
  assert(new_knot_ix > orig.suppMin(d) && new_knot_ix < orig.suppMax(d));

  // defining the following two constants for code readability purposes
  const int degree = orig.degree(d);
  const double new_kval = kvals[new_knot_ix];
  const double y1           = kvals[orig.suppMin(d)];
  const double y_deg_plus_2 = kvals[orig.suppMax(d)];
  const double y2           = kvals[orig.kvec(d)[1]];
  const double y_deg_plus_1 = kvals[orig.kvec(d)[degree]];

  // compute alpha_1 and alpha_2 multipliers
  const double a1 = (new_kval >= y_deg_plus_1) ? 1 : (new_kval - y1) / (y_deg_plus_1 - y1);
  const double a2 = (new_kval <= y2          ) ? 1 : (y_deg_plus_2 - new_kval) / (y_deg_plus_2 - y2);

  // compute new gamma coefficients
  const double g1 = orig.gamma() * a1;
  const double g2 = orig.gamma() * a2;
  
  //wcout << "Setting gamma to: " << g1 << " " << g2 << endl;
  
  // compute new control points
  const Point c_g1 = orig.coefTimesGamma() * a1;
  const Point c_g2 = orig.coefTimesGamma() * a2;

  // making new knotvector (copy old knotvector and insert the new knot at the correct place
  vector<int> kvec_new(orig.kvec(d));
  kvec_new.insert(std::find_if(kvec_new.begin(), 
			       kvec_new.end(), 
			       [new_knot_ix](int ix) {return ix >= new_knot_ix;}),
		  new_knot_ix);
  
  // sorting out knotvec iterators
  const vector<int>::const_iterator k1_u = (d == XFIXED) ? kvec_new.begin() : orig.kvec(XFIXED).begin();
  const vector<int>::const_iterator k1_v = (d == XFIXED) ? orig.kvec(YFIXED).begin() : kvec_new.begin();
  const vector<int>::const_iterator k2_u = (d == XFIXED) ? k1_u + 1 : k1_u;
  const vector<int>::const_iterator k2_v = (d == XFIXED) ? k1_v     : k1_v + 1;

  new_1 = LRBSpline2D(c_g1, orig.degree(XFIXED), orig.degree(YFIXED), k1_u, k1_v, g1);
  new_2 = LRBSpline2D(c_g2, orig.degree(XFIXED), orig.degree(YFIXED), k2_u, k2_v, g2);

}


//==============================================================================
void LRBSpline2DUtils::recursively_split(const LRBSpline2D& b_orig,
					 const Mesh2D& mesh,
					 std::vector<LRBSpline2D>& result)
//==============================================================================
{
  //wcout << "Enter 'recursively_split' at level: " << level << endl;
  const int umin = b_orig.suppMin(XFIXED);
  const int vmin = b_orig.suppMin(YFIXED);
  const int umax = b_orig.suppMax(XFIXED);
  const int vmax = b_orig.suppMax(YFIXED);
  const vector<int> m_kvec_u = derive_knots(mesh, XFIXED, umin, umax, vmin, vmax);
  const vector<int> m_kvec_v = derive_knots(mesh, YFIXED, vmin, vmax, umin, umax);

  // @@ The assertions below should always hold if function is called with correct 
  // argument. When code is properly debugged and tested, they can be taken away for
  // efficiency (asserts should go away anyway when compiling in optimized mode).
  // Alternatively, if it is a concern that users might call this function with wrong 
  // argument, the assertions could be replaced by exception-throwing 'if'-statements.
  assert(std::includes(m_kvec_u.begin(), m_kvec_u.end(), b_orig.kvec(XFIXED).begin(), b_orig.kvec(XFIXED).end()));
  assert(std::includes(m_kvec_v.begin(), m_kvec_v.end(), b_orig.kvec(YFIXED).begin(), b_orig.kvec(YFIXED).end()));

  if (num_inner_knots(m_kvec_u) > num_inner_knots(b_orig.kvec(XFIXED))) {
    // Since we know that m_kvec_u contains more elements than b_orig.kvec(XFIXED) and since
    // we know that the latter is included in the former, we know that there must be at least
    // one knot in 'm_kvec_u' that is not found in b_orig.kvec(XFIXED).  We can therefore call
    // the following function without risking an exception to be thrown.
    const int new_ix = find_uncovered_inner_knot(m_kvec_u, b_orig.kvec(XFIXED));

    LRBSpline2D new_fct_1, new_fct_2;
    split_function(b_orig, XFIXED, mesh.knotsBegin(XFIXED), new_ix, new_fct_1, new_fct_2);
    recursively_split(new_fct_1, mesh, result); //, level+1);
    recursively_split(new_fct_2, mesh, result); //, level+1);

  } else if (num_inner_knots(m_kvec_v) > num_inner_knots(b_orig.kvec(YFIXED))) {
    // same comment as above
    const int new_ix = find_uncovered_inner_knot(m_kvec_v, b_orig.kvec(YFIXED));
    
    LRBSpline2D new_fct_1, new_fct_2;
    split_function(b_orig, YFIXED, mesh.knotsBegin(YFIXED), new_ix, new_fct_1, new_fct_2);
    recursively_split(new_fct_1, mesh, result); // , level+1);
    recursively_split(new_fct_2, mesh, result); // , level+1);

  } else {
    // This function has minimal support in the mesh, and should not be further split
    result.push_back(b_orig);
  }
}

//==============================================================================
// if 'b' can be split at least once in the mesh 'm', split it once, and return the 
// result through 'b1' and 'b2'.  The function never carries out more than one split, 
// even when several splits are possible.
bool LRBSpline2DUtils::try_split_once(const LRBSpline2D& b, const Mesh2D& mesh, 
				      LRBSpline2D& b1, LRBSpline2D& b2)
//==============================================================================
{
  const int umin = b.suppMin(XFIXED);
  const int vmin = b.suppMin(YFIXED);
  const int umax = b.suppMax(XFIXED);
  const int vmax = b.suppMax(YFIXED);
  const vector<int> m_kvec_u = derive_knots(mesh, XFIXED, umin, umax, vmin, vmax);
  const vector<int> m_kvec_v = derive_knots(mesh, YFIXED, vmin, vmax, umin, umax);

  // @@ The assertions below should always hold if function is called with correct 
  // argument. When code is properly debugged and tested, they can be taken away for
  // efficiency (asserts should go away anyway when compiling in optimized mode).
  // Alternatively, if it is a concern that users might call this function with wrong 
  // argument, the assertions could be replaced by exception-throwing 'if'-statements.
  assert(std::includes(m_kvec_u.begin(), m_kvec_u.end(), b.kvec(XFIXED).begin(), b.kvec(XFIXED).end()));
  assert(std::includes(m_kvec_v.begin(), m_kvec_v.end(), b.kvec(YFIXED).begin(), b.kvec(YFIXED).end()));

  if (num_inner_knots(m_kvec_u) > num_inner_knots(b.kvec(XFIXED))) {
    // Since we know that m_kvec_u contains more elements than b.kvec(XFIXED) and since
    // we know that the latter is included in the former, we know that there must be at least
    // one knot in 'm_kvec_u' that is not found in b.kvec(XFIXED).  We can therefore call
    // the following function without risking an exception to be thrown.
    const int new_ix = find_uncovered_inner_knot(m_kvec_u, b.kvec(XFIXED));

    split_function(b, XFIXED, mesh.knotsBegin(XFIXED), new_ix, b1, b2);
    return true;

  } else if (num_inner_knots(m_kvec_v) > num_inner_knots(b.kvec(YFIXED))) {
    // same comment as above
    const int new_ix = find_uncovered_inner_knot(m_kvec_v, b.kvec(YFIXED));
    split_function(b, YFIXED, mesh.knotsBegin(YFIXED), new_ix, b1, b2);
    return true;
  } 
  // No splits possible
  return false;
}
}; // end namespace Go
