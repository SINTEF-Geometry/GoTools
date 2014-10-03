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

#include <assert.h>
#include <stdexcept>
//#include <iostream> // @@ debug only
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"


using namespace std;

namespace Go
{
namespace {

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


  // Derive the knotvector resulting from moving along the u-direction (if d == XFIXED) or 
  // v-direction (if d == YFIXED) from 'beg' to 'end', and  with multiplicities equal to 
  // the 'nu' value of the segment in the ortogonal direction, starting at 'orto_min' and 
  // extending to 'orto_max'.
  //------------------------------------------------------------------------------
  vector<int> 
  LRBSpline2DUtils::derive_knots(const Mesh2D& m, Direction2D d, int beg, int end, 
				 int orto_min, int orto_max)
  //------------------------------------------------------------------------------
  {
    vector<int> result;
    for (int pos = beg; pos <= end; ++pos) 
      {
	int start = m.nu(d, pos, orto_min, orto_max);
	result.insert(result.end(), start, pos); // 0 multiplicities will not be inserted
      }

    return result;
  }


//==============================================================================
  void LRBSpline2DUtils::split_several(const double* knotvals, const vector<int>& k_vec_in, const vector<int>& new_knots,
				       vector<int>& k_vec_out, vector<double>& b_spline_weigths)
//==============================================================================
  {
    int deg = (int)k_vec_in.size() - 2;
    int nmb_new_knots = (int)new_knots.size();

    // k_vec_out is resized to its final size, and k_vec_in is inserted at the beginning
    k_vec_out.resize(deg + 2 + nmb_new_knots);
    for (int i = 0; i <= deg + 1; ++i)
      k_vec_out[i] = k_vec_in[i];

    // b_spline_weigths is resized to its dinal size
    // 1.0 is inserted in the first position, 0.0 at the others. This is important, as the algorithm is
    // using the 0.0 after the end of the current weight list at each step in the insertions
    b_spline_weigths.resize(nmb_new_knots + 1);
    for (int i = 0; i <= nmb_new_knots; ++i)
      b_spline_weigths[i] = (i == 0) ? 1 : 0;

    // Insert the knots, one at a time
    for (int j = 0; j < nmb_new_knots; ++j)
      {
	int new_knot_idx = new_knots[j];
	double new_knot_val = knotvals[new_knot_idx];

	// first_geq_new_idx is the position in k_vec_out of the first knot index that is >= the new knot index
	int first_geq_new_idx = 0;
	for (;k_vec_out[first_geq_new_idx] < new_knot_idx; ++first_geq_new_idx);

	// insert_pos is the position in k_vec_out of the first knot index that is > the new knot index
	// This is the same as the insert position in k_vec_out of the new knot index
	int insert_pos = first_geq_new_idx;
	for (;k_vec_out[insert_pos] <= new_knot_idx; ++insert_pos);

	// Determine where to start changing the coefficients.
	// B-splines with their second-to-last knot before (or at) the inserted knot, are not affected, so they are skipped
	double prev_old_weigth = 0.0;
	int start_i = insert_pos - deg;
	if (start_i <= 0)
	  start_i = 0;
	else
	  prev_old_weigth = b_spline_weigths[start_i - 1];

	// Calculate the new B-spline coefficients
	for (int i = start_i; i <= j + 1; ++i)
	  {
	    double next_old_weigth = b_spline_weigths[i];
	    if (i >= first_geq_new_idx)
	      b_spline_weigths[i] = prev_old_weigth;
	    else
	      {
		double factor = (new_knot_val - knotvals[k_vec_out[i]]) / (knotvals[k_vec_out[i + deg]] - knotvals[k_vec_out[i]]);
		b_spline_weigths[i] = next_old_weigth * factor + prev_old_weigth * (1.0 - factor);
	      }
	    prev_old_weigth = next_old_weigth;
	  }

	// Insert the new knot
	for (int i = deg + j + 1; i >= insert_pos; --i)
	  k_vec_out[i + 1] = k_vec_out[i];
	k_vec_out[insert_pos] = new_knot_idx;
      }
  }

  
//==============================================================================
  void LRBSpline2DUtils::split_function(const LRBSpline2D& orig, 
					const Mesh2D& mesh,
					Direction2D d, 
					const double* const kvals,
					int new_knot_ix,
					LRBSpline2D*& new_1,
					LRBSpline2D*& new_2)
//==============================================================================
{
  assert(new_knot_ix > orig.suppMin(d) && new_knot_ix < orig.suppMax(d));

#ifndef NDEBUG
  if (orig.rational())
  {
      ;//MESSAGE("Rational case under implementation, expect errors!");
  }
#endif

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
  const bool rat = orig.rational();
  // If rational case we incorporate the scaling into the weight.
  // @@sbr201301 But what about rational case and gamma differerent from 1.0?
  const Point c_g1 = (rat) ? orig.coefTimesGamma() : orig.coefTimesGamma() * a1;
  const Point c_g2 = (rat) ? orig.coefTimesGamma() : orig.coefTimesGamma() * a2;
//  const Point c_g2 = orig.coefTimesGamma() * a2;

  const double weight = orig.weight();
  // compute new weights
  const double w1 = (rat) ? a1*weight : 1.0;
  const double w2 = (rat) ? a2*weight : 1.0;

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

  new_1 = new LRBSpline2D(c_g1, w1, orig.degree(XFIXED), 
			  orig.degree(YFIXED), 
			  k1_u, k1_v, g1, &mesh, rat);
  new_1->setFixCoef(orig.coefFixed());
  new_2 = new LRBSpline2D(c_g2, w2, orig.degree(XFIXED), 
			  orig.degree(YFIXED), 
			  k2_u, k2_v, g2, &mesh, rat);
  new_2->setFixCoef(orig.coefFixed());

}


//==============================================================================
// if 'b' can be split at least once in the mesh 'm', split it once, and return the 
// result through 'b1' and 'b2'.  The function never carries out more than one split, 
// even when several splits are possible.
bool LRBSpline2DUtils::try_split_once(const LRBSpline2D& b, const Mesh2D& mesh, 
				      LRBSpline2D*& b1, 
				      LRBSpline2D*& b2)
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

    // @@@ VSK. Cannot set pointers to the new bsplines before it is placed in the
    // global array. Will the position of the element change when a bspline is removed?
    split_function(b, mesh, XFIXED, mesh.knotsBegin(XFIXED), new_ix, b1, b2);
    return true;

  } else if (num_inner_knots(m_kvec_v) > num_inner_knots(b.kvec(YFIXED))) {
    // same comment as above
    const int new_ix = find_uncovered_inner_knot(m_kvec_v, b.kvec(YFIXED));
    split_function(b, mesh, YFIXED, mesh.knotsBegin(YFIXED), new_ix, b1, b2);
    return true;
  } 
  // No splits possible
  return false;
}

}; // end namespace Go
