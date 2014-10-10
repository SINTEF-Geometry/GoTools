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

#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"
#include "GoTools/utils/checks.h"
#include "GoTools/geometry/SplineSurface.h"

//------------------------------------------------------------------------------

using std::vector;
using std::set;
using std::tuple;
using std::get;
using std::find_if;
using std::find;
using std::unique_ptr;

namespace Go
{

namespace {
  // Anonymous namespace for tools for building parameter functions into one-dimensional LR-spline surface function

  // Struct to hold necessary data for a 2D B-spline in the splitting process, see method insertParameterFunctions() below
  struct greville_lrbspline
  {
    // Sort key to make sure an LR B-spline is always split into B-splines with smaller sort key.
    // It is the sum of the keys for each parameter direction, created by greville_lrbspline_key() below
    int sort_key_;

    // The index knot vectors in each parameter direction
    vector<int> knots_u_;
    vector<int> knots_v_;

    // The gamma value of the LR B-spline (same as LRBSpline2D::gamma_)
    double gamma_;

    // The LRBSpline2D::coef_times_gamma_ values for the control points building the linear functions (u,v) |-> u and (u,v) |-> v
    double gamma_times_u_;
    double gamma_times_v_;

    greville_lrbspline()
    {
    }

    greville_lrbspline(int sort_key,
		       const vector<int>::const_iterator& knots_u, const vector<int>::const_iterator& knots_v,
		       int deg_u, int deg_v,
		       double gamma, double gamma_times_x, double gamma_times_y)
      : sort_key_(sort_key),
	gamma_(gamma), gamma_times_u_(gamma_times_x), gamma_times_v_(gamma_times_y)
    {
      knots_u_.resize(deg_u + 2);
      knots_v_.resize(deg_v + 2);
      copy(knots_u, knots_u + (deg_u + 2), knots_u_.begin());
      copy(knots_v, knots_v + (deg_v + 2), knots_v_.begin());
    }

    // Return -1 if this instance comes before the input b-spline defined by the given knot vectors, and given precalculated key
    // Return 1 if this instance comes after the input b-spline
    // Return 0 if they are equal
    int compare(int key, const vector<int>::const_iterator knots_u_start, const vector<int>::const_iterator knots_v_start) const
    {
      // Compare sort keys first. High key values come before small
      if (sort_key_ != key)
	return (sort_key_ > key) ? -1 : 1;

      // Compare knot vectors
      int len_u = (int)knots_u_.size();
      for (int i = 0 ; i < len_u; ++i)
	if (knots_u_[i] != knots_u_start[i])
	  return (knots_u_[i] < knots_u_start[i]) ? -1 : 1;

      int len_v = (int)knots_v_.size();
      for (int i = 0 ; i < len_v; ++i)
	if (knots_v_[i] != knots_v_start[i])
	  return (knots_v_[i] < knots_v_start[i]) ? -1 : 1;

      return 0;
    }
  };


  // Sort key for a univariate B-spline given by its index knot vector, the key is used to sort B-splines such that the when splitting a B-spline, the
  // new B-splines always have a smaller key
  // The key is given as (S+1)*O - (M1 + M2) where
  // - S is the index span of the knots
  // - O is the polynomial order
  // - M1 is the multiplicity of the first knot
  // - M2 is the multiplicity of the last knot
  // This is the same as S*O - M1 + P2, where P2 is the position of the last knot smaller than the end knot
  int greville_lrbspline_key(const vector<int>::const_iterator knots_start, int deg)
  {
    int first_knot = knots_start[0];
    int last_knot = knots_start[deg + 1];
    int pos_after_first = 1;
    for (;knots_start[pos_after_first] == first_knot; ++pos_after_first);
    int pos_before_last = deg + 1;
    for (;knots_start[pos_before_last] == last_knot; --pos_before_last);
    return (deg + 1) * (last_knot - first_knot) + pos_after_first - pos_before_last;
  }


  // Searches for the position of a specific 2D B-spline in the sorted list of B-splines in insertParameterFunctions(), and tells
  // whether the B-spline already existed in the list or not
  // Notice that the list, in order to save time, is already preallocated and presized. Also notice that the storage of the list might wrap,
  // i.e. the list might go to the end of the vector and continue from the top
  // - key            is the sort key of the 2D B-spline, as stored in  greville_lrbspline::sort_key_
  // - knots_u_start  is the start of the first knot vector of the LR B-spline
  // - knots_v_start  is the start of the second knot vector of the LR B-spline
  // - spline_list    is the vector holding the sorted list (possibly wrapped) of the splines to search in
  // - start_pos      is the position of the first element in the list. We must have 0 <= start_pos < spline_list.size()
  // - end_pos        is the position after the last element in the list. We must have start_pos <= end_pos <= start_pos + spline_list.size()
  //                  If end_pos > spline_list.size() the list wraps. i.e.
  //                  * If end_pos <= spline_list.size() then the list is stored in
  //                         spline_list[start_pos], ... , spline_list[end_pos-1]
  //                  * If end_pos > spline_list.size() then the list is stored in
  //                         spline_list[start_pos], ... , spline_list[spline_list.size()-1], spline_list[0], ... , spline_list[end_pos-spline_list.size()-1]
  // - exists         Ends up as true if the B-spline already existed in the list, false if not
  // returns          An integer p where start_pos <= p <= end_pos. If the B-spline existed in the list (exists == true), it is found in position p.
  //                  If the B-spline did not exist in the list, it should be inserted at position p.
  //                  Notice that p < spline_list.size() refers to position p in the vector, and
  //                  p >= spline_list.size() refers to position p-spline_list.size() in the vector, and
  int spline_list_pos(int key, const vector<int>::const_iterator knots_u_start, const vector<int>::const_iterator knots_v_start,
		      const vector<greville_lrbspline>& spline_list, int start_pos, int end_pos, bool& exists)
  {
    int spline_list_size = spline_list.size();

    int left = start_pos;    // Everything before position 'left' is known to be before the test spline
    int right = end_pos;     // Everything at or after position 'right' is known to be after the test spline

    while (left < right)
      {
	int mid = (left + right) / 2;
	int compare_result = spline_list[mid - (mid < spline_list_size ? 0 : spline_list_size)].compare(key, knots_u_start, knots_v_start);
	if (compare_result == -1)
	  left = mid + 1;
	else if (compare_result == 1)
	  right = mid;
	else
	  {
	    exists = true;
	    return mid;
	  }
      }

    exists = false;
    return left;
  }


  // Inserts a new 2D B-spline in the sorted list of B-splines in insertParameterFunctions(). It is assumed that the vector spline_list has
  // at least one free position. The vector holding the list might wrap (see comments on spline_list_pos())
  // - glb            The b-spline to be inserted
  // - spline_list    Is the vector holding the sorted list (possibly wrapped) of the splines to insert into
  // - insert_pos     The insert position
  // - end_pos        The position after the last element in the list. We must have insert_pos <= end_pos < 2 * spline_list.size()
  void spline_list_insert(const greville_lrbspline& glb, vector<greville_lrbspline>& spline_list, int insert_pos, int end_pos)
  {
    int spline_list_size = spline_list.size();

    // The insert position and/or the end position might be greater than or equal to the list size, the the actual storage position has wrapped around and started
    // from the beginning of the vector. We must handle the different cases

    // If the insert postion (and then also the end position) indicate wrapping, we just subtract both of them and handle them as if neither was wrapped
    if (insert_pos >= spline_list_size)
      {
	insert_pos -= spline_list_size;
	end_pos -= spline_list_size;
      }

    // If the end position (and then also the insert position) does not indicate wrapping, the block to be moved does not wrap, then just move it one position
    if (end_pos < spline_list_size)
      {
	for (int i = end_pos; i > insert_pos; --i)
	  spline_list[i] = spline_list[i - 1];
      }

    // Otherwise, we wrap. First move the last block (with actual storage position from the from the beginning of the vector) down one position,
    // then move the one element now being wrapped over, and finally move the first block (with actual storage position down to one before the end of the list)
    else
      {
	for (int i = end_pos - spline_list_size; i > 0; --i)
	  spline_list[i] = spline_list[i - 1];
	spline_list[0] = spline_list[spline_list_size - 1];
	for (int i = spline_list_size - 1; i > insert_pos; --i)
	  spline_list[i] = spline_list[i - 1];
      }

    // Finally we insert
    spline_list[insert_pos] = glb;
  }



  // Extend the size of the sorted list of B-splines in insertParameterFunctions().
  // The method is called because the storage vector is full, this means that the list elements are, in order, found as
  //    spline_list[top_spline_list], spline_list[spline_list.size()-1], spline_list[0], ... spline_list[top_spline_list-1]
  // - spline_list    The vector holding the sorted list (most possibly wrapped) of the splines
  // - top_spline_list  The position in the vector of the first spline in the list.
  // - extend_size      The number of elements to add to the size of the vector
  void extend_spline_list(vector<greville_lrbspline>& spline_list, int top_spline_list, int extend_size)
  {
    int old_spline_list_size = spline_list.size();
    int new_size = old_spline_list_size + extend_size;
    spline_list.reserve(new_size);
    spline_list.resize(new_size);

    // Move elements originally wrapped around to their new positions. Notice that wrapping might still take place afterwards
    // if extend_size < top_spline_list
    for (int i = 0, new_pos = old_spline_list_size; i < top_spline_list; ++i)
      {
	spline_list[new_pos] = spline_list[i];
	++new_pos;
	if (new_pos == new_size)
	  new_pos = 0;
      }
  }


}; // end anonymous namespace


//------------------------------------------------------------------------------
  LRSplineSurface::ElementMap 
  LRSplineUtils::identify_elements_from_mesh(const Mesh2D& m)
//------------------------------------------------------------------------------
{
  // The element entities are created, but no bspline information is
  // attached at this stage
  LRSplineSurface::ElementMap emap;
  for (Mesh2DIterator mit = m.begin(); mit != m.end(); ++mit)
    {
      unique_ptr<Element2D> elem(new Element2D(m.kval(XFIXED, (*mit)[0]),
					       m.kval(YFIXED, (*mit)[1]),
					       m.kval(XFIXED, (*mit)[2]),
					       m.kval(YFIXED, (*mit)[3])));
      // emap[LRSplineSurface::generate_key(elem)] = elem;
      LRSplineSurface::ElemKey e_key = LRSplineSurface::generate_key(m.kval(XFIXED, (*mit)[0]),
							       m.kval(YFIXED, (*mit)[1]));
      emap.insert(std::make_pair(e_key, std::move(elem)));
    }

  return emap;
}

//------------------------------------------------------------------------------
// Add ('remove'=false) or remove ('remove'=true) references to a particular 
// LRBSpline2D for all affected elements in 'emap'.  
void LRSplineUtils::update_elements_with_single_bspline(LRBSpline2D* b, 
							LRSplineSurface::ElementMap& emap, 
							const Mesh2D& mesh,
							bool remove)
//------------------------------------------------------------------------------
{
  const double* kvals_x = mesh.knotsBegin(XFIXED);
  const double* kvals_y = mesh.knotsBegin(YFIXED);

  //@@@ VSK. The LRBSpline2D knows wich elements are in its support. Should be
  // possible to make this more effective
  // The element list in the LR B-spline must be up to date
  for (int y = b->suppMin(YFIXED); y != b->suppMax(YFIXED); ++y) {
    for (int x = b->suppMin(XFIXED); x != b->suppMax(XFIXED); ++x) {

      // -- following two lines are meant as an optimization, but the effect seems
      // -- to be quite marginal. @@
      if (mesh.nu(XFIXED, x, y, y+1) < 1) continue; // this cannot be the lower-left -
      if (mesh.nu(YFIXED, y, x, x+1) < 1) continue; // corner of an element

      const LRSplineSurface::ElemKey key = {kvals_x[x], kvals_y[y]};
      auto it = emap.find(key);
      if (it != emap.end()) {
	// We found an element covered by the support of this basis function.
	// Update the element
	if (remove) { // we want to remove it (if it's there)
	  it->second->removeSupportFunction(b);
	} else {      // we want to insert it (if it's not there already)
	  it->second->addSupportFunction(b);
	  b->addSupport(it->second.get());  // Update bspline with respect to element
	}
      }
    }
  }
}




//------------------------------------------------------------------------------
int LRSplineUtils::locate_interval(const Mesh2D& m, Direction2D d, double value, 
				   double other_value, bool at_end)
//------------------------------------------------------------------------------
{
  const double u = (d == XFIXED) ?       value : other_value;
  const double v = (d == XFIXED) ? other_value :       value;

  int u_ix, v_ix;
  const bool found = (at_end) ?
    Mesh2DUtils::identify_patch_upper_right(m, u, v, u_ix, v_ix) : 
    Mesh2DUtils::identify_patch_lower_left (m, u, v, u_ix, v_ix);

  if (!found) 
    THROW("locate_interval() : parameter outside domain.");
  
  return (d == XFIXED) ? u_ix : v_ix;
}

//------------------------------------------------------------------------------
// increment all indices in the B-spline knotvecs in the given direction, if they
// are equal to or larger to 'from_ix' (This function is useful when new meshlines
// are inserted, causing existing meshlines to get higher indices than before.
void LRSplineUtils::increment_knotvec_indices(LRSplineSurface::BSplineMap& bmap, 
					      Direction2D d, 
					      int from_ix)
//------------------------------------------------------------------------------
{
  for (auto b = bmap.begin(); b != bmap.end(); ++b) { // b is here a <key, value> pair, where 'value' is a LRBSpline2D
    vector<int>& kvec = b->second->kvec(d);
    if (from_ix <= kvec.back()) 
      for (auto k = kvec.begin(); k != kvec.end(); ++k)
  	if (*k >= from_ix) 
  	  ++*k;
  }
}


//------------------------------------------------------------------------------
// returns a pointer to the new (or existing) function
LRBSpline2D* 
LRSplineUtils::insert_basis_function(unique_ptr<LRBSpline2D>& b, 
				     const Mesh2D& mesh, 
				     LRSplineSurface::BSplineMap& bmap)
//------------------------------------------------------------------------------
{
  // Add a bspline to the global pool of bsplines, but check first if it
  // exists already. In that case, the bspline scaling factor is updated
  auto key = LRSplineSurface::generate_key(*b, mesh);
  if (bmap.find(key) != bmap.end()) {

    // combine b with the function already present
    LRBSpline2D* target = bmap[key].get();
    target->gamma()            += b->gamma();
    target->coefTimesGamma() += b->coefTimesGamma();

    return target;
  } 
  // if we got here, there is no pre-existing basis function like b.  
  // return (bmap[key] = b).get();
  // bmap[key] = b;
//  bmap.insert(std::pair<LRSplineSurface::BSKey, unique_ptr<LRBSpline2D> >(key, b));
  // unique_ptr<LRBSpline2D> dummy_ptr;
  // std::pair<LRSplineSurface::BSKey, unique_ptr<LRBSpline2D> > key_b(key, dummy_ptr);
  // std::swap(b, key_b.second);
//  bmap.insert(key_b);//std::make_pair(key, b));
  bmap.insert(std::make_pair(key, std::move(b)));
  return b.get();
}

// For each line of the mesh in the given direcion, set the multiplicity of all meshrectangles
// to the highest multiplicity for that line.  Return a vector giving the multiplicity for
// each line.
//------------------------------------------------------------------------------
vector<int> LRSplineUtils::set_uniform_meshlines(Direction2D d, Mesh2D& mesh)
//------------------------------------------------------------------------------
{
  vector<int> mults(mesh.numDistinctKnots(d));
  const int line_len = mesh.numDistinctKnots(flip(d)) - 1;
  for (int i = 0; i != mesh.numDistinctKnots(d); ++i) {
    const int mul = mesh.largestMultInLine(d, i);
    mesh.setMult(d, i, 0, line_len, mul);
    mults[i] = mul;
  }
  return mults;
}

//------------------------------------------------------------------------------
bool LRSplineUtils::all_meshlines_uniform(Direction2D d, const Mesh2D& m)
//------------------------------------------------------------------------------
{
  const int orto_intervals = m.numDistinctKnots(flip(d)) - 1; 
  for (int i = 0; i != m.numDistinctKnots(d); ++i) 
    if (m.nu(d, i, 0, orto_intervals) != m.largestMultInLine(d, i)) return false;
  return true;
}

//------------------------------------------------------------------------------
double LRSplineUtils::compute_greville(const vector<int>& v_ixs, 
				       const double* const vals)
//------------------------------------------------------------------------------
{
  assert(v_ixs.size() > 2); // degree of spline should be at least 1.
  const int num_pts = (int)v_ixs.size() - 2;
  double result = 0;
  for (int i = 1; i != num_pts+1; ++i) {
    result += vals[v_ixs[i]];
  }
  return result / num_pts;
}

//------------------------------------------------------------------------------
  vector<double> LRSplineUtils::compute_greville(int deg, const vector<int>& k_vec_in, const double* const knotvals)
//------------------------------------------------------------------------------
  {
    int nmb_sb_splines = (int)k_vec_in.size() - (deg + 1);
    vector<double> greville_vals(nmb_sb_splines);
    for (int i = 0; i < nmb_sb_splines; ++i)
      {
	double sum = 0.0;
	for (int j = 1; j <= deg; ++j)
	  sum += knotvals[k_vec_in[i + j]];
	greville_vals[i] = sum / (double)deg;
      }
    return greville_vals;
  }


// Helper function used by 'tensor_split'.  Identifies the knots that need to 
// be inserted into some reference vector, in order to bring multiplicty of each 
// knot up to the multiplicity expressed in 'mult' (including those with a
// implicit current multiplicity of zero.  This function is used by 'tensor_split'.
//------------------------------------------------------------------------------
vector<int> LRSplineUtils::knots_to_insert(const vector<int>& ref, 
					   const vector<int>& mults)
//------------------------------------------------------------------------------
{
  assert(ref.size() > 1);
  vector<int> result;
  const int first = ref.front();
  const int last  = ref.back();
  assert(last > first);
  
  // Determining multiplicities of all internal knots (incl. those with conceptual
  // multiplicity of 0) of 'ref'
  vector<int> ref_internal_mults(last-(first+1), 0);
  for (auto k = ref.begin(); k != ref.end(); ++k) 
    if (*k > first && *k < last) 
      ++ref_internal_mults[(*k)-(first+1)];

  // determining missing knots, by comparing the current multiplicity of the knots,
  // as compared with the target multiplicities, as expressed in the 'mults' vector.
  for (size_t i = 0; i != ref_internal_mults.size(); ++i) {
    const size_t knot_ix = i + (first+1);
    const int missing_mul = mults[knot_ix] - ref_internal_mults[i];
    assert(missing_mul >= 0);
    result.insert(result.end(), missing_mul, (int)knot_ix);
  }
  return result;
}

// Compute the alpha coefficients required by function 'insert_knots'.  These are recursively
// computed, and express the coefficients of the various subdivided of an original univariate 
// bspline-function (covering the span of 'oldvec_ix'), after a number of knot insertions resulting
// in the vector 'newvec_ix'.  These are the coefficients computed by the Oslo Algorithm
//------------------------------------------------------------------------------
double LRSplineUtils::compute_alpha(int degree, 
				    const int* const oldvec_ix, 
				    const int* const newvec_ix, 
				    const double* const kvals)
//------------------------------------------------------------------------------
{
  if (degree == 0) 
    return ((newvec_ix[0] < oldvec_ix[1]) && (newvec_ix[0] >= oldvec_ix[0]) ? 1 : 0);
 
  const double nv_d   = kvals[newvec_ix[degree]];
  const double ov_0   = kvals[oldvec_ix[0]];
  const double ov_1   = kvals[oldvec_ix[1]];
  const double ov_d   = kvals[oldvec_ix[degree]];
  const double ov_dp1 = kvals[oldvec_ix[degree+1]];

  const double fac1 = (ov_d - ov_0) > 0 ? (nv_d - ov_0) / (ov_d - ov_0) : 0;
  const double fac2 = (ov_dp1 - ov_1) > 0 ? (ov_dp1 - nv_d) / (ov_dp1 - ov_1) : 0;

  return 
    ((fac1 > 0) ? fac1 * compute_alpha(degree-1, oldvec_ix, newvec_ix, kvals) : 0) +
    ((fac2 > 0) ? fac2 * compute_alpha(degree-1, oldvec_ix+1, newvec_ix, kvals) : 0);
}

// Insert knots of 'new_knots' into one of the univariate knotvectors of 'bfun' (as 
// specified by 'd'), using the Oslo Algorithm.
// The coefficients of the b-spline functions that are produced are returned as the first 
// component of the tuple.  The univariate knot vector after insertion is returned as the
// second component of the tuple.
// This function is used by 'tensor_split'.
// It is assumed that 'new_knots' is sorted.
//------------------------------------------------------------------------------
tuple<vector<double>, vector<int>> 
  LRSplineUtils::insert_knots(const vector<int>& new_knots,
			      unique_ptr<LRBSpline2D>& bfun,
			      const Direction2D d,
			      const double* const kvals)
//------------------------------------------------------------------------------
{
  // setting up structure of the returned result
  tuple<vector<double>, vector<int> > result(vector<double>(new_knots.size() + 1), 
					     bfun->kvec(d));
  vector<double>& alpha = get<0>(result);
  vector<int>& kvec_final = get<1>(result);
  
  // inserting new knots into knotvector, and ensuring that they are in the correct
  // order (nondecreasing)
  kvec_final.insert(kvec_final.end(), new_knots.begin(), new_knots.end());
  sort(kvec_final.begin(), kvec_final.end());

  // computing the alpha multiplication factors that are used to express 'bfun' as 
  // a linear sum of its (new_knots.size() + 1) subdivided parts.
  for (size_t i = 0; i != new_knots.size() + 1; ++i) {
    alpha[i] = compute_alpha(bfun->degree(d), &bfun->kvec(d)[0], &kvec_final[i], kvals);
  }
  return result;
}

// Efficiently split the function 'bfun' up according to a full tensor product mesh.
// It is assumed that 'tensor_mesh' is a tensor mesh.
//------------------------------------------------------------------------------
void LRSplineUtils::tensor_split(unique_ptr<LRBSpline2D>& bfun, 
				 const vector<int>& x_mults,
				 const vector<int>& y_mults,
				 const Mesh2D& tensor_mesh,
				  LRSplineSurface::BSplineMap& bmap)
//------------------------------------------------------------------------------
{
  const vector<int> kx = knots_to_insert(bfun->kvec(XFIXED), x_mults);
  const vector<int> ky = knots_to_insert(bfun->kvec(YFIXED), y_mults);

  const double* const x_kvals = tensor_mesh.knotsBegin(XFIXED);
  const double* const y_kvals = tensor_mesh.knotsBegin(YFIXED);
  
  const auto x_coefs_kvec = insert_knots(kx, bfun, XFIXED, x_kvals);
  const auto y_coefs_kvec = insert_knots(ky, bfun, YFIXED, y_kvals);
  const vector<double>& x_coefs = get<0>(x_coefs_kvec);
  const vector<double>& y_coefs = get<0>(y_coefs_kvec);
  const vector<int>& x_knots = get<1>(x_coefs_kvec);
  const vector<int>& y_knots = get<1>(y_coefs_kvec);

  const int deg_x = bfun->degree(XFIXED);
  const int deg_y = bfun->degree(YFIXED);
  const double gamma = bfun->gamma();
  const Point& c_g = bfun->coefTimesGamma();
  const double weight = bfun->weight();
  const bool rational = bfun->rational();

  for (int iy = 0; iy != (int)y_coefs.size(); ++iy) {
    const double yc = y_coefs[iy];
    for (int ix = 0; ix != (int)x_coefs.size(); ++ix) {
      const double xc = x_coefs[ix];
      unique_ptr<LRBSpline2D> basis(new LRBSpline2D(c_g*yc*xc,
						    weight,
						    deg_x,
						    deg_y,
						    &x_knots[ix], 
						    &y_knots[iy], 
						    yc * xc * gamma,
						    &tensor_mesh,
						    rational));
      insert_basis_function(basis, tensor_mesh, bmap);
    }
  }
}

//------------------------------------------------------------------------------
void LRSplineUtils::iteratively_split (vector<unique_ptr<LRBSpline2D> >& bfuns, 
				       const Mesh2D& mesh)
//------------------------------------------------------------------------------
{
  // The following set is used to keep track over unique b-spline functions.   
  // b-spline function is here identified by its knotvectors only, as we already
  // assume that degrees are unchanging.  Also, since we expect to find several
  // component of a given b-spline-function, and these must be added up, we will
  // not look at control points or gamma coefficients to determine uniqueness.

//  set<unique_ptr<LRBSpline2D>, LRBSpline2DUtils::support_compare> tmp_set;
//  set<LRBSpline2D*, LRBSpline2DUtils::support_compare> tmp_set;
  set<LRBSpline2D*, support_compare> tmp_set;

//  unique_ptr<LRBSpline2D> b_split_1, b_split_2;
  bool split_occurred;

  // this closure adds b_spline functions to tmp_set, or combine them if they 
  // are already in it
//  auto insert_bfun_to_set = [&tmp_set](unique_ptr<LRBSpline2D>& b)->void {
  auto insert_bfun_to_set = [&tmp_set](LRBSpline2D* b)->bool {
    auto it = tmp_set.find(b);
    if (it == tmp_set.end()) {  // not already in set
      tmp_set.insert(b);
      return true;
    } else {
    // combine b with the function already present
      bool rat = b->rational();
      if (rat)
	{ // We must alter the weight of the second basis function to match that of our reference.
	  // We multiply the coefTimesGamma.
	  double b_w = b->weight();
	  double it_w = (*it)->weight();
	  double weight = b_w + it_w;
	  // We must rescale the coefs to reflect the change in weight.
	  b->coefTimesGamma() *= b_w/weight;
	  (*it)->coefTimesGamma() *= it_w/weight;
	  b->weight() = (*it)->weight() = weight;
	  // (*it)->gamma() += b->gamma();
	  // (*it)->coefTimesGamma() += b->coefTimesGamma();
	}
      (*it)->gamma() += b->gamma();
      (*it)->coefTimesGamma() += b->coefTimesGamma();
      return false;
    }
  };

  // After a new knot is inserted, there might be bsplines that are no longer
  // minimal. Split those according to knot line information in the mesh
  // keep looping until no more basis functions were inserted
  do {
    tmp_set.clear();
    split_occurred = false;

    for (auto b = bfuns.begin(); b != bfuns.end(); ++b) {
      LRBSpline2D *b_split_1 = NULL;
      LRBSpline2D *b_split_2 = NULL;
      if (LRBSpline2DUtils::try_split_once(*(*b), mesh, b_split_1, b_split_2)) {
	// this function was splitted.  Throw it away, and keep the two splits
	bool was_inserted = insert_bfun_to_set(b_split_1);
	if (!was_inserted)
	  delete b_split_1;
	was_inserted = insert_bfun_to_set(b_split_2);
	if (!was_inserted)
	  delete b_split_2;
	split_occurred = true;
      } else {
	// this function was not split.  Keep it.
	insert_bfun_to_set(b->get());
	// We must also release the function from the unique_ptr in the bfuns vector.
	b->release();
      }
    }

    // moving the collected bsplines over to the vector
    bfuns.clear();
//    bfuns.resize(tmp_set.size());
    int cntr = 0;
    for (auto b_kv = tmp_set.begin(); b_kv != tmp_set.end(); ++b_kv) 
      {
	// bfuns.insert(bfuns.end(), std::move(*b_kv));
//	std::swap(bfuns[cntr], *b_kv);
	// unique_ptr<LRBSpline2D> ptr = std::move(*b_kv);
	// bfuns.emplace_back(std::move(ptr));
//	std::swap(bfuns[cntr], std::move(*b_kv));
	bfuns.insert(bfuns.end(), unique_ptr<LRBSpline2D>(*b_kv));
	++cntr;
      }

  } while (split_occurred);
}

//------------------------------------------------------------------------------
void LRSplineUtils::iteratively_split2 (vector<LRBSpline2D*>& bsplines,
					const Mesh2D& mesh,
					LRSplineSurface::BSplineMap& bmap,
					double domain[])
//------------------------------------------------------------------------------
{
  // The following set is used to keep track over unique b-spline functions.   
  // b-spline function is here identified by its knotvectors only, as we already
  // assume that degrees are unchanging.  Also, since we expect to find several
  // component of a given b-spline-function, and these must be added up, we will
  // not look at control points or gamma coefficients to determine uniqueness.
  // @@sbr201212 It is not suited for tmp_set since we are not guaranteed to
  // have the boundary elements, resulting in a miss if a basis function with
  // the same support is the last element.
  set<LRBSpline2D*, support_compare> tmp_set;

  bool split_occurred;

  // this closure adds b_spline functions to tmp_set, or combine them if they 
  // are already in it
  auto insert_bfun_to_set = [&tmp_set](LRBSpline2D* b,
				       LRSplineSurface::BSplineMap& bmap,
				       double domain[])->bool
    {
      auto it = tmp_set.find(b);
      bool overlap = b->overlaps(domain);
      LRBSpline2D* other = NULL;
      if (it != tmp_set.end())
	{
	  other = (*it);
	}
      else if (!overlap)
	{
	  // Search for an identical B-spline in the entire domain
	  LRSplineSurface::BSKey key = LRSplineSurface::generate_key(*b);
	  auto it3 = bmap.find(key);
	  if (it3 != bmap.end() && it3->second.get() != b)
	    other = it3->second.get();
	}
      if (it == tmp_set.end() && !other)
	{  
	  // We must check if the last element of tmp_set is equal.
	  auto it2 = tmp_set.end();
	  int set_size_pre = tmp_set.size();
	  if (set_size_pre > 0)
	    it2--;
	  bool last_elem_equal = (set_size_pre > 0 && (support_equal(*it2, b)));
	  if (last_elem_equal)
	    MESSAGE("DEBUG: Last element is equal to new element!");
	  // not already in set
	  tmp_set.insert(b);
	  int set_size_post = tmp_set.size();
	  if (set_size_pre == set_size_post)
	    MESSAGE("DEBUG: It seems we tried to insert an element already present!");
	  return true;
	}
      else
	{
	  bool rat = b->rational();
	  if (rat)
	    { // We must alter the weight of the second basis function to match that of our reference.
	      double b_w = b->weight();
	      double it_w = other->weight();
	      double weight = b_w + it_w;//0.66*b_w + 0.34*it_w;
	      // We must rescale the coefs to reflect the change in weight.
	      b->coefTimesGamma() *= b_w/weight; // c_1*w_1 = c_1*(w_1/w_n)*w_n.
	      other->coefTimesGamma() *= it_w/weight;
	      b->weight() = (*it)->weight() = weight;
	    }
	  // combine b with the function already present
	  other->gamma() += b->gamma();
	  other->coefTimesGamma() += b->coefTimesGamma();
	  // We update the support of b with its replacement.
	  std::vector<Element2D*>::iterator it2 = b->supportedElementBegin();
	  for (it2; it2 < b->supportedElementEnd(); ++it2)
	    {
	      // If there exists a support function already (such as b) it is overwritten.
	      (*it2)->addSupportFunction(other);
	      other->addSupport(*it2);
	    }

	  // Finally we remove all elements from b.
	  while (b->nmbSupportedElements() > 0)
	    {
	      auto it2 = b->supportedElementBegin();
	      b->removeSupport(*it2);
	    }

	  return false;
	}
    };

  // After a new knot is inserted, there might be bsplines that are no longer
  // minimal. Split those according to knot line information in the mesh
  // keep looping until no more basis functions were inserted
#ifndef NDEBUG
  int deb_iter = 0;
#endif

  vector<unique_ptr<LRBSpline2D> > added_basis;

  do { // Loop is run until no more splits occur.
    tmp_set.clear(); // Used to store new basis functions for each iteration.
    split_occurred = false;

#ifndef NDEBUG
//    std::cout << "deb_iter: " << deb_iter << std::endl;
      vector<LRBSpline2D*> tmp_set_vec, tmp_set_supp_supp_vec;
      vector<Element2D*> tmp_set_supp_vec;
//      for (auto iter = tmp_set.begin(); iter != tmp_set.end(); ++iter)
      for (auto iter = bsplines.begin(); iter != bsplines.end(); ++iter)
	{
	  tmp_set_vec.push_back((*iter));
	  tmp_set_supp_vec.insert(tmp_set_supp_vec.end(), (*iter)->supportedElementBegin(), (*iter)->supportedElementEnd());
	}
      for (auto iter = tmp_set_supp_vec.begin(); iter != tmp_set_supp_vec.end(); ++iter)
	{
	  tmp_set_supp_supp_vec.insert(tmp_set_supp_supp_vec.end(), (*iter)->supportBegin(), (*iter)->supportEnd());
	}      

      std::sort(tmp_set_supp_supp_vec.begin(), tmp_set_supp_supp_vec.end());
      tmp_set_supp_supp_vec.erase(std::unique(tmp_set_supp_supp_vec.begin(), tmp_set_supp_supp_vec.end()),
				  tmp_set_supp_supp_vec.end());

      std::sort(tmp_set_supp_vec.begin(), tmp_set_supp_vec.end());
      tmp_set_supp_vec.erase(std::unique(tmp_set_supp_vec.begin(), tmp_set_supp_vec.end()), tmp_set_supp_vec.end());
      // @@sbr201212 puts("Remove when done debugging!");
#endif

    int ki = 0;
    for (auto b = bsplines.begin(); b != bsplines.end(); ++b, ++ki) {

      LRBSpline2D *b_split_1 = NULL;
      LRBSpline2D *b_split_2 = NULL;
      if (LRBSpline2DUtils::try_split_once(*(*b), mesh, b_split_1, b_split_2)) {
     	// this function was splitted.  Throw it away, and keep the two splits
	// @@@ VSK. Must also update bmap and set element pointers
	// Fetch all elements
	vector<Element2D*> elements = (*b)->supportedElements();

	// Remove bspline from element
	for (size_t kr=0; kr<elements.size(); ++kr)
	  {
#ifndef NDEBUG
//	    std::cout << "DEBUG: ki = " << ki << ", kr = " << kr << ", deb_iter = " << deb_iter << std::endl;
#endif
	    elements[kr]->removeSupportFunction(*b);
	  }

	// Remove bspline from bspline map
	LRSplineSurface::BSKey key = LRSplineSurface::generate_key(*(*b));
	auto it = bmap.find(key);
	if (it != bmap.end())
	  bmap.erase(it);
	else
	  {
	    // Remove the bspline from the vector of bsplines to add
	    for (size_t kr=0; kr<added_basis.size(); ++kr)
	      if (added_basis[kr].get() == (*b))
		{
		  std::swap(added_basis[kr], added_basis[added_basis.size()-1]);
		  added_basis.pop_back();
		  break;
		}
	  }

	// // Add new bsplines to the bspline map
	// LRSplineSurface::BSKey key1 = LRSplineSurface::generate_key(*b_split_1);
	// LRSplineSurface::BSKey key2 = LRSplineSurface::generate_key(*b_split_2);
	// bmap[key1] = b_split_1;
	// bmap[key2] = b_split_2;

	// Until the elements are split, let the new bsplines store all
	// elements from their origin in their support
#if 0
	LRSplineSurface::BSKey key1 = LRSplineSurface::generate_key(*b_split_1);
	auto it1 = tmp_set.find(key1);
	if (it1 != tmp_set.end())
	  {
	    MESSAGE("b_split_1 already present in tmp_set!");
	  }

	LRSplineSurface::BSKey key2 = LRSplineSurface::generate_key(*b_split_2);
	auto it2 = tmp_set.find(key2);
	if (it2 != tmp_set.end())
	  {
	    MESSAGE("b_split_2 already present in tmp_set!");
	  }
#endif
	// LRSplineSurface::BSKey key2 = LRSplineSurface::generate_key(*b_split_1);
	// auto iter = bsplines.size();
	// Since the elements have not yet been split, the support is the same.
	b_split_1->setSupport(elements);
	b_split_2->setSupport(elements);

    	if (insert_bfun_to_set(b_split_1, bmap, domain)) // @@sbr deb_iter==0 && ki == 20. ref==4.
	  {
	    // A new LRBspline is created, remember it
	    added_basis.push_back(unique_ptr<LRBSpline2D>(b_split_1));
	    // Let the elements know about the new bsplines
	    for (size_t kr=0; kr<elements.size(); ++kr)
	      if (b_split_1->overlaps(elements[kr]))
		elements[kr]->addSupportFunction(b_split_1);
	      else
		{
#if 0//ndef NDEBUG
		  MESSAGE("No overlap!"); // @@sbr201212 This should not happen.
#endif
		  b_split_1->removeSupport(elements[kr]);
		  elements[kr]->removeSupportFunction(b_split_1);
		}
	  }
	else
	  { // Memory management.
	    delete b_split_1;
	  }

    	if (insert_bfun_to_set(b_split_2, bmap, domain))
	  {
	    // A new LRBspline is created, remember it
	    added_basis.push_back(unique_ptr<LRBSpline2D>(b_split_2));
	    // Let the elements know about the new bsplines
	    for (size_t kr=0; kr<elements.size(); ++kr)
	      if (b_split_2->overlaps(elements[kr]))
		elements[kr]->addSupportFunction(b_split_2);
	      else
		{
#if 0//ndef NDEBUG
		  MESSAGE("No overlap!"); // @@sbr201212 This should not happen.
#endif
		  b_split_2->removeSupport(elements[kr]);
		  elements[kr]->removeSupportFunction(b_split_2);
		}
	  }
	else
	  { // Memory management.
	    delete b_split_2;
	  }

    	split_occurred = true;
      } else {
     	// this function was not split.  Keep it.
     	bool was_inserted = insert_bfun_to_set(*b, bmap, domain);
	if (!was_inserted)
	  {
//	    MESSAGE("DEBUG: We should remove basis function from added_basis!");
	    // Remove the bspline from the vector of bsplines to add
	    for (size_t kr=0; kr<added_basis.size(); ++kr)
	      if (added_basis[kr].get() == (*b))
		{
		  std::swap(added_basis[kr], added_basis[added_basis.size()-1]);
		  added_basis.pop_back();
		  break;
		}
	  }
      }
     }

    // moving the collected bsplines over to the vector
    bsplines.clear();
    for (auto b_kv = tmp_set.begin(); b_kv != tmp_set.end(); ++b_kv) 
      {
	bsplines.push_back(*b_kv);
      }

#ifndef NDEBUG
    ++deb_iter;
#endif

  } while (split_occurred);

  // Add new basis functions to bmap
  for (size_t kr=0; kr<added_basis.size(); ++kr)
    {
      LRSplineSurface::BSKey key = LRSplineSurface::generate_key(*added_basis[kr]);
      auto it = bmap.find(key);
      if (it != bmap.end())
	{ // @@ I guess we handle this by adding 
	  MESSAGE("Already added to map! This is a bug. Expect core dump if not fixed ...");
	  // @@sbr201305 This will in a lost pointer and most likely a core dump!
	  LRBSpline2D* b = added_basis[kr].get();
#if 1
	  bool rat = b->rational();
	  if (rat)
	    { // We must alter the weight of the second basis function to match that of our reference.
	      double b_w = b->weight();
	      double it_w = (it->second)->weight();
	      double weight = b_w + it_w;//0.66*b_w + 0.34*it_w;
	      // We must rescale the coefs to reflect the change in weight.
	      b->coefTimesGamma() *= b_w/weight; // c_1*w_1 = c_1*(w_1/w_n)*w_n.
	      (it->second)->coefTimesGamma() *= it_w/weight;
	      b->weight() = (it->second)->weight() = weight;
	    }
	  // combine b with the function already present
	  (it->second)->gamma() += b->gamma();
	  (it->second)->coefTimesGamma() += b->coefTimesGamma();
	  // We update the support of b with its replacement.
	  std::vector<Element2D*>::iterator it2 = b->supportedElementBegin();
	  for (it2; it2 < b->supportedElementEnd(); ++it2)
	    {
	      // If there exists a support function already (such as b) it is overwritten.
	      (*it2)->addSupportFunction(it->second.get());
	      (it->second)->addSupport(*it2);
	    }
#endif

	  // We locate the corresponding pointer in bsplines.
	  for (auto iter = bsplines.begin(); iter != bsplines.end(); ++iter)
	    if (*iter == b)
	      {
		bsplines.erase(iter);
		break;
	      }
	}
      else
	{
	  bmap.insert(std::make_pair(key, std::move(added_basis[kr])));
	}
    }
      
}

//------------------------------------------------------------------------------
// Inserts a refinement into the mesh and increments indices in the BSplineMap accordingly.
// The function returns three integers:  
//   - The first is the index of the meshline on which the inserted meshrectangle is located
//     (index of the fixed parameter value in the global knot vector.  
//   - The second is the start index of the mesh rectangle in the other parameter direction.
//   - The third is the end index of the mesh rectangle in the other parameter direction.
//
// NB: Does NOT carry out any splitting of the basis functions in the BSplineMap.  It is the
//     caller's responsibility to do this afterwards!  (As such, this function is not 
//     intended for use by itself, but should only be called from one of the LRSpline::refine()
//     methods.

  tuple<int, int, int, int>
LRSplineUtils::refine_mesh(Direction2D d, double fixed_val, double start, 
			   double end, int mult, bool absolute,
			   int spline_degree, double knot_tol,
			   Mesh2D& mesh,  LRSplineSurface::BSplineMap& bmap)
//------------------------------------------------------------------------------
{
  if (mult > spline_degree + 1) 
    THROW("Cannot refine with multiplicity higher than degree+1.");

  // Determining the anchor points of the new meshline to be inserted (it must end in an existing
  // ortogonal meshline of multiplicity >= 1).
  // @@ EXPLAIN THE WEIRD USE OF SIGNS AND TOLERANCE BELOW
  double del = end - start;
  const int start_ix = locate_interval(mesh, flip(d), start + del * knot_tol, fixed_val, false);
  const int   end_ix = locate_interval(mesh, flip(d), end   - del * knot_tol, fixed_val,  true);
  // const int start_ix = locate_interval(mesh, flip(d), start + fabs(start) * knot_tol, fixed_val, false);
  // const int   end_ix = locate_interval(mesh, flip(d), end   - fabs(end)   * knot_tol, fixed_val,  true);
  //const int start_ix = locate_interval(mesh, flip(d), start * (1 + knot_tol), fixed_val, false);
  //const int   end_ix = locate_interval(mesh, flip(d), end   * (1 - knot_tol), fixed_val,  true);

  // Fetch the last nonlarger knot value index in the fixed direction
  int prev_ix = Mesh2DUtils::last_nonlarger_knotvalue_ix(mesh, d, fixed_val);

  int fixed_ix; // to be set below. Knot value index of the new knot
  // const auto existing_it = find_if(mesh.knotsBegin(d),
  // 				   mesh.knotsEnd(d), 
  // 				   [=](double x)->bool
  // 				   {return (abs(fixed_val - x) < abs(fixed_val) * knot_tol);});

  // if (existing_it != mesh.knotsEnd(d)) { // increase multiplicity of existing meshrectangles
  if (fabs(mesh.kval(d, prev_ix) - fixed_val) < knot_tol)
    {
      // increase multiplicity of existing meshrectangles
      //fixed_ix = int(existing_it - mesh.knotsBegin(d));
      fixed_ix = prev_ix;

      // check that the proposed multiplicity modification is legal
      for (int i = start_ix; i < end_ix; ++i) {
	const int cur_m = mesh.nu(d, fixed_ix, i, i+1);
	if (absolute && (cur_m > mult)) 
	  THROW("Cannot decrease multiplicity.");
	else if (!absolute && (cur_m+mult > spline_degree + 1)) 
	  THROW("Cannot increase multiplicity.");
      }
      // set or increment multiplicity
      absolute ? 
	mesh.setMult(d, fixed_ix, start_ix, end_ix, mult) :
	mesh.incrementMult(d, fixed_ix, start_ix, end_ix, mult);

    } else { 
    // insert a new line, and set relevant part to desired multiplicity
    // @@sbr Should this perhaps result in a warning? Or do we want to
    // increase grid automatically?
    fixed_ix = mesh.insertLine(d, fixed_val, 0); 
    mesh.setMult(d, fixed_ix, start_ix, end_ix, mult);

    // change index of _all_ basis functions who refer to knot values with indices >= inserted one
    increment_knotvec_indices(bmap, d, fixed_ix);
  }

  // We must also update the mesh in the basis functions.
  auto it = bmap.begin();
  while (it != bmap.end())
    {
      it->second->setMesh(&mesh);
      ++it;
    }

  // If this prev_ix corresponds to our fixed_val, we decrease the value.
  // @@sbr201212 I guess we could do this earlier, but then we need to update the working code above ...
  if (mesh.kval(d, prev_ix) == fixed_val)
    prev_ix -= mult;
  if (prev_ix < 0) // We must handle cases near the start.
    prev_ix = 0;

   return tuple<int, int, int, int>(prev_ix, fixed_ix, start_ix, end_ix);
}

bool LRSplineUtils::support_equal(const LRBSpline2D* b1, const LRBSpline2D* b2)
{
  // to compare b1 and b2, compare the x-knotvectors.  If these are identical, compare
  // the y-knotvectors instead.
  const int tmp1 = compare_seq(b1->kvec(XFIXED).begin(), b1->kvec(XFIXED).end(),
			       b2->kvec(XFIXED).begin(), b2->kvec(XFIXED).end());
  if (tmp1 != 0)
    return false;
  const int tmp2 = compare_seq(b1->kvec(YFIXED).begin(), b1->kvec(YFIXED).end(),
			       b2->kvec(YFIXED).begin(), b2->kvec(YFIXED).end());
  return (tmp2 == 0);
};

bool LRSplineUtils::elementOK(const Element2D* elem, const Mesh2D& m)
{
  MESSAGE("Under construction!");

  // We check that all functions in the support overlap.
   return true;
}

//==============================================================================
  void LRSplineUtils::insertParameterFunctions(LRSplineSurface* lr_spline_sf)
//==============================================================================
  {
    const Mesh2D& mesh = lr_spline_sf->mesh();
    int deg_u = lr_spline_sf->degree(XFIXED);
    int deg_v = lr_spline_sf->degree(YFIXED);
    int knots_u_size = mesh.numDistinctKnots(XFIXED);
    int knots_v_size = mesh.numDistinctKnots(YFIXED);
    const double* knots_begin_u = mesh.knotsBegin(XFIXED);
    const double* knots_begin_v = mesh.knotsBegin(YFIXED);

    // First we start with the biggest tensor product submesh
    vector<int> sub_tensor_knots_u;
    for (int i = 0; i < knots_u_size; ++i)
      {
	int nu = mesh.nu(XFIXED, i, 0, knots_v_size - 1);
	for (int j = 0; j < nu; ++j)
	  sub_tensor_knots_u.push_back(i);
      }
    vector<int> sub_tensor_knots_v;
    for (int i = 0; i < knots_v_size; ++i)
      {
	int nu = mesh.nu(YFIXED, i, 0, knots_u_size - 1);
	for (int j = 0; j < nu; ++j)
	  sub_tensor_knots_v.push_back(i);
      }

    // Get greville values for each direction in tensor product submesh
    vector<double> greville_u = compute_greville(deg_u, sub_tensor_knots_u, knots_begin_u);
    vector<double> greville_v = compute_greville(deg_v, sub_tensor_knots_v, knots_begin_v);
    int st_u_size = (int)greville_u.size();
    int st_v_size = (int)greville_v.size();

    // Sort keys for the origianlly inserted univariate B-splines in each direction
    vector<int> stk_key_u(st_u_size);
    vector<int> stk_key_v(st_v_size);
    vector<int>::const_iterator stk_it_u = sub_tensor_knots_u.begin();
    vector<int>::const_iterator stk_it_v = sub_tensor_knots_v.begin();
    for (int i = 0; i < st_u_size; ++i)
      stk_key_u[i] = greville_lrbspline_key(stk_it_u + i, deg_u);
    for (int i = 0; i < st_v_size; ++i)
      stk_key_v[i] = greville_lrbspline_key(stk_it_v + i, deg_v);

    // The vector holding the list of B-splines we will be working in.
    // We preallocate now to have one continous block of data for lookup
    // The list might be 'wrapped' in the vector, see comments below
    vector<greville_lrbspline> spline_list;
    int spline_list_size_step = (int)((double)(lr_spline_sf->numBasisFunctions()) * 1.2);
    int spline_list_size = spline_list_size_step;
    spline_list.reserve(spline_list_size);
    spline_list.resize(spline_list_size);

    // Positions for first element in list, and one beyond the last element in the list.
    // We always have
    //    0 <= top_spline_list < spline_list.size()
    //    top_spline_list <= end_pos_spline_list <= top_spline_list + spline_list.size()
    // If end_pos_spline_list <= spline_list.size() then the list is stored as
    //     spline_list[top_spline_list], ..., spline_list[end_pos_spline_list-1]
    // If end_pos_spline_list > spline_list.size() then the list is stored as
    //     spline_list[top_spline_list], ..., spline_list[spline_list.size()-1], spline_list[0], ..., spline_list[end_pos_spline_list-spline_list.size()-1]
    int top_spline_list = 0;
    int end_pos_spline_list = 0;

    // Insert the starting B-splines from the tensor mesh, before refinement
    for (int j = 0; j < st_v_size; ++j)
      for (int i = 0; i < st_u_size; ++i)
	{
	  int key = stk_key_u[i] + stk_key_v[j];
	  greville_lrbspline glb(key, stk_it_u + i, stk_it_v + j, deg_u, deg_v, 1.0, greville_u[i], greville_v[j]);
	  bool dummy;
	  int insert_pos = spline_list_pos(key, stk_it_u + i, stk_it_v + j, spline_list, 0, end_pos_spline_list, dummy);
	  spline_list_insert(glb, spline_list, insert_pos, end_pos_spline_list);
	  ++end_pos_spline_list;
	}

    // We now repeat the splitting algorithm down to minimal support B-splines.
    // I.e. we iterate through all B-splines in the list, in sort order. Each is either split as much as possible,
    // or recognized as a minimal support B-spline, then the control point of the corresponding B-spline in
    // the LR-spline surface function is updated
    while (top_spline_list < end_pos_spline_list)
      {

	// Fetch B-spline on top of list and update start postion
	greville_lrbspline current_spline = spline_list[top_spline_list];
	if (++top_spline_list == spline_list_size)
	  {
	    top_spline_list = 0;
	    end_pos_spline_list -= spline_list_size;
	  }

	// Minimnum and maximum knot index, and multiplicities at minimum knot index
	int min_knot_u = current_spline.knots_u_[0];
	int max_knot_u = current_spline.knots_u_[deg_u + 1];
	int min_knot_v = current_spline.knots_v_[0];
	int max_knot_v = current_spline.knots_v_[deg_v + 1];

	int mult_u_min = 0;
	int mult_v_min = 0;
	while (current_spline.knots_u_[++mult_u_min] == min_knot_u);
	while (current_spline.knots_v_[++mult_v_min] == min_knot_v);

	// Search for knots to insert, u-direction
	vector<int> new_u;
	for (int knot_idx = min_knot_u + 1, first_pos_geq_knot = mult_u_min; knot_idx < max_knot_u; ++knot_idx)
	  {
	    // Get highest possible multiplicity inside B-spline domain for given knot
	    int mu = mesh.nu(XFIXED, knot_idx, min_knot_v, max_knot_v);

	    // Subtract existing multiplicity for B-spline at given knot
	    while (current_spline.knots_u_[first_pos_geq_knot] == knot_idx)
	      {
		--mu;
		++first_pos_geq_knot;
	      }

	    // Add knots to be inserted
	    for (; mu > 0; --mu)
	      new_u.push_back(knot_idx);
	  }

	// Search for knots to insert, v-direction
	vector<int> new_v;
	for (int knot_idx = min_knot_v + 1, first_pos_geq_knot = mult_v_min; knot_idx < max_knot_v; ++knot_idx)
	  {
	    int mu = mesh.nu(YFIXED, knot_idx, min_knot_u, max_knot_u);
	    while (current_spline.knots_v_[first_pos_geq_knot] == knot_idx)
	      {
		--mu;
		++first_pos_geq_knot;
	      }
	    for (; mu > 0; --mu)
	      new_v.push_back(knot_idx);
	  }

	// Now we have the knots to be inserted. Either we found none, then the B-spline is minimal, or some were found, then we split

	if (new_u.size() == 0 && new_v.size() == 0)
	  {

	    // Current B-spline is minimal. Look up corresponding B-spline in LR-Spline surface and update control point

	    // Calculate multiplicities at maximum knot index
	    int mult_u_max = 0;
	    int mult_v_max = 0;
	    while (current_spline.knots_u_[deg_u + 1 -(++mult_u_max)] == max_knot_u);
	    while (current_spline.knots_v_[deg_v + 1 -(++mult_v_max)] == max_knot_v);

	    // Look up B-spline and update control point
	    LRSplineSurface::BSplineMap::iterator bspl_it =
	      lr_spline_sf->bsplineFromDomain(knots_begin_u[min_knot_u], knots_begin_v[min_knot_v], knots_begin_u[max_knot_u], knots_begin_v[max_knot_v],
					      mult_u_min, mult_v_min, mult_u_max, mult_v_max);
	    if (bspl_it->second->dimension() != 1)
	      THROW("Minimal support LR B-spline did not have dimension 1 as expected");

	    const double z_gamma = bspl_it->second->coefTimesGamma()[0];
	    bspl_it->second->coefTimesGamma() = Point(current_spline.gamma_times_u_, current_spline.gamma_times_v_, z_gamma);
	  }

	else
	  {
	    // Current B-spline is not of minimal support - perform splitting

	    // Insert knots and get weights and expanded knot vectors in each parameter direction
	    vector<int> full_knots_u;
	    vector<int> full_knots_v;
	    vector<double> weights_u;
	    vector<double> weights_v;
	    LRBSpline2DUtils::split_several(knots_begin_u, current_spline.knots_u_, new_u, full_knots_u, weights_u);
	    LRBSpline2DUtils::split_several(knots_begin_v, current_spline.knots_v_, new_v, full_knots_v, weights_v);

	    // Build sort keys for each parameter direction
	    int new_u_size = (int)weights_u.size();
	    int new_v_size = (int)weights_v.size();
	    vector<int> new_key_u(new_u_size);
	    vector<int> new_key_v(new_v_size);
	    vector<int>::const_iterator full_knots_u_begin = full_knots_u.begin();
	    vector<int>::const_iterator full_knots_v_begin = full_knots_v.begin();
	    for (int i = 0; i < new_u_size; ++i)
	      new_key_u[i] = greville_lrbspline_key(full_knots_u_begin + i, deg_u);
	    for (int i = 0; i < new_v_size; ++i)
	      new_key_v[i] = greville_lrbspline_key(full_knots_v_begin + i, deg_v);

	    // Run through all new LR B-splines after knot insertions, and either create or update gamma and control values
	    for (int j = 0; j < (int)weights_v.size(); ++j)
	      for (int i = 0; i < (int)weights_u.size(); ++i)
		{
		  double weight = weights_u[i] * weights_v[j];
		  int key = new_key_u[i] + new_key_v[j];
		  bool b_spline_exists;
		  int new_bspline_pos = spline_list_pos(key, full_knots_u_begin + i, full_knots_v_begin + j,
							spline_list, top_spline_list, end_pos_spline_list, b_spline_exists);

		  if (b_spline_exists)
		    {

		      // Update existing B-spline.
		      int existing_glb_idx = new_bspline_pos - (new_bspline_pos < spline_list_size ? 0 : spline_list_size);
		      spline_list[existing_glb_idx].gamma_ += weight * current_spline.gamma_;
		      spline_list[existing_glb_idx].gamma_times_u_ += weight * current_spline.gamma_times_u_;
		      spline_list[existing_glb_idx].gamma_times_v_ += weight * current_spline.gamma_times_v_;
		    }

		  else
		    {

		      // Create and insert new B-spline
		      greville_lrbspline new_glb(key, full_knots_u_begin + i, full_knots_v_begin + j, deg_u, deg_v,
						 weight * current_spline.gamma_,
						 weight * current_spline.gamma_times_u_,
						 weight * current_spline.gamma_times_v_);

		      // Before inserting, test if vector holding the splines to test is full, expand if necessary
		      if (end_pos_spline_list == top_spline_list + spline_list_size)
			{
			  extend_spline_list(spline_list, top_spline_list, spline_list_size_step);
			  spline_list_size += spline_list_size_step;
			}

		      // Insert
		      spline_list_insert(new_glb, spline_list, new_bspline_pos, end_pos_spline_list);
		      ++end_pos_spline_list;
		    }
		}
	  }
      }

    // Finally we test if all LR B-splines in the LR-spline surface have been upgraded
    LRSplineSurface::BSplineMap::const_iterator bspl_end = lr_spline_sf->basisFunctionsEnd();
    for (LRSplineSurface::BSplineMap::const_iterator bspl_it = lr_spline_sf->basisFunctionsBegin(); bspl_it != bspl_end; ++bspl_it)
      {
	if (bspl_it->second->dimension() != 3)
	  THROW("Not all minimal support LR B-spline have raised their dimension to 3");
      }
  }


SplineSurface* LRSplineUtils::fullTensorProductSurface(const LRSplineSurface& lr_spline_sf)
{
  bool is_full_tensor = (lr_spline_sf.isFullTensorProduct());
  shared_ptr<LRSplineSurface> lr_spline_sf_copy;
  const LRSplineSurface* full_tp_sf = &lr_spline_sf;
  if (!is_full_tensor)
    {
      lr_spline_sf_copy = shared_ptr<LRSplineSurface>(lr_spline_sf.clone());
      lr_spline_sf_copy->expandToFullTensorProduct();
      full_tp_sf = lr_spline_sf_copy.get();
    }

  int num_coefs_u, num_coefs_v;
  int order_u = full_tp_sf->degree(XFIXED) + 1;
  int order_v = full_tp_sf->degree(YFIXED) + 1;
  int dim = full_tp_sf->dimension();
  bool rational = full_tp_sf->rational();

  // The basis functions should be ordered with increasing u-parameter first, then the v-parameter.
  auto iter = full_tp_sf->basisFunctionsBegin();
//  int kdim = (rational) ? dim + 1 : dim;
  vector<double> sf_coefs;
  while (iter != full_tp_sf->basisFunctionsEnd())
    {
      Point coef = iter->second->coefTimesGamma();//Coef();
      sf_coefs.insert(sf_coefs.end(), coef.begin(), coef.end());
      if (rational)
	sf_coefs.insert(sf_coefs.end(), iter->second->weight());

      ++iter;
    }

  const Mesh2D& mesh = full_tp_sf->mesh();
  vector<double> knots_u = mesh.getKnots(XFIXED, 0);
  vector<double> knots_v = mesh.getKnots(YFIXED, 0);
  num_coefs_u = knots_u.size() - order_u;
  num_coefs_v = knots_v.size() - order_v;

  SplineSurface* spline_sf = new SplineSurface(num_coefs_u, num_coefs_v,
					       order_u, order_v,
					       knots_u.begin(), knots_v.begin(),
					       sf_coefs.begin(),
					       dim,
					       rational);

  return spline_sf;
}

vector<vector<double> > LRSplineUtils::elementLineClouds(const LRSplineSurface& lr_spline_sf)
{
  MESSAGE("elementLineClouds(): Not yet implemented.");

  // We make sure that the edges of the elements contain the corners of neighbouring elements.
  vector<vector<double> > dummy;
  return dummy;

}


}; // end namespace Go

