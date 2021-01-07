//===========================================================================
//                                                                           
// File: LRBSpline3DUtils.C                                                  
//                                                                           
// Created: Thu Apr 18 08:17:21 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <assert.h>
#include <stdexcept>
//#include <iostream> // @@ debug only
#include "GoTools/lrsplines3D/LRBSpline3DUtils.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"


using namespace std;


namespace Go
{


// Anonymous namespace
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
  const auto tmp = std::mismatch(kvec2.begin() + (m2 - 1), kvec2.end(), kvec1.begin() + (m1 - 1));

  if (tmp.second == kvec1.end()) throw runtime_error("No uncovered inner knot found.");

  return *tmp.second;
}


} // end anonymous namespace



// Derive the knotvector resulting from moving along the u-direction (if d == XDIR) or
// v-direction (if d == YDIR) or w-direction (if d == ZDIR) from 'beg' to 'end', and
// with multiplicities equal to the 'nu' value of the segment in the ortogonal direction,
// starting at 'orto_min' and extending at 'orto_max'.
//------------------------------------------------------------------------------
vector<int>
LRBSpline3DUtils::derive_knots(const Mesh3D& m, Direction3D d, int beg, int end,
                               int orto1_min, int orto1_max,
                               int orto2_min, int orto2_max)
//------------------------------------------------------------------------------
{
  vector<int> result;
  for (int pos = beg; pos <= end; ++pos) {
    int start =  m.nu(d, pos, orto1_min, orto1_max, orto2_min, orto2_max);
    result.insert(result.end(), start, pos);  // 0 multiplicities will not be inserted
  }

  return result;
}



//==============================================================================
  void LRBSpline3DUtils::split_function(const LRBSpline3D& orig, 
					Direction3D d, 
					const double* const kvals,
					int new_knot_ix,
					vector<unique_ptr<BSplineUniLR> >& bspline_vec1,
					vector<unique_ptr<BSplineUniLR> >& bspline_vec2,
					vector<unique_ptr<BSplineUniLR> >& bspline_vec3,
					LRBSpline3D*& new_1,
					LRBSpline3D*& new_2)
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
  size_t kvec_size = kvec_new.size();
  kvec_new.insert(std::find_if(kvec_new.begin(), 
			       kvec_new.end(), 
			       [new_knot_ix](int ix) {return ix >= new_knot_ix;}),
		  new_knot_ix);
  
  // sorting out knotvec iterators
  const vector<int>::const_iterator k1_u = (d == XDIR) ? kvec_new.begin() : orig.kvec(XDIR).begin();
  const vector<int>::const_iterator k1_v = (d == YDIR) ? kvec_new.begin() : orig.kvec(YDIR).begin();
  const vector<int>::const_iterator k1_w = (d == ZDIR) ? kvec_new.begin() : orig.kvec(ZDIR).begin();
  const vector<int>::const_iterator k2_u = (d == XDIR) ? k1_u + 1 : k1_u;
  const vector<int>::const_iterator k2_v = (d == YDIR) ? k1_v + 1 : k1_v;
  const vector<int>::const_iterator k2_w = (d == ZDIR) ? k1_w + 1 : k1_w;

  // Identify univariate bspline
  BSplineUniLR *uni1_1=NULL, *uni1_2=NULL, *uni2_1=NULL, *uni2_2=NULL, *uni3_1=NULL, *uni3_2=NULL;
  if (d == XDIR)
    {
      int ix = 0;
      bool found1 = 
	BSplineUniUtils::identify_bsplineuni(k1_u, k1_u+kvec_size, bspline_vec1, ix);
      if (!found1)
       {
	 // Must create univariate B-spline
	 BSplineUniLR *origuni = orig.getUnivariate(XDIR);
	 BSplineUniLR *uninew = new BSplineUniLR(origuni->pardir(), origuni->degree(),
						 k1_u, origuni->getMesh());

	 // Insert in array
	 BSplineUniUtils::insert_univariate(bspline_vec1, uninew, ix);
       }
      uni1_1 = bspline_vec1[ix].get();
      bool found2 = 
	BSplineUniUtils::identify_bsplineuni(k2_u, k2_u+kvec_size, bspline_vec1, ix);
      if (!found2)
       {
	 // Must create univariate B-spline
	 BSplineUniLR *origuni = orig.getUnivariate(XDIR);
	 BSplineUniLR *uninew = new BSplineUniLR(origuni->pardir(), origuni->degree(),
						 k2_u, origuni->getMesh());

	 // Insert in array
	 BSplineUniUtils::insert_univariate(bspline_vec1, uninew, ix);
       }
      uni1_2 = bspline_vec1[ix].get();
      uni2_1 = uni2_2 = orig.getUnivariate(YDIR);
      uni3_1 = uni3_2 = orig.getUnivariate(ZDIR);
    }
  else if (d == YDIR)
    {
      int ix = 0;
      bool found1 = 
	BSplineUniUtils::identify_bsplineuni(k1_v, k1_v+kvec_size, bspline_vec2, ix);
     if (!found1)
       {
	 // Must create univariate B-spline
	 BSplineUniLR *origuni = orig.getUnivariate(YDIR);
	 BSplineUniLR *uninew = new BSplineUniLR(origuni->pardir(), origuni->degree(),
						 k1_v, origuni->getMesh());

	 // Insert in array
	 BSplineUniUtils::insert_univariate(bspline_vec2, uninew, ix);
       }
       uni2_1 = bspline_vec2[ix].get();
      bool found2 = 
	BSplineUniUtils::identify_bsplineuni(k2_v, k2_v+kvec_size, bspline_vec2, ix);
      if (!found2)
       {
	 // Must create univariate B-spline
	 BSplineUniLR *origuni = orig.getUnivariate(YDIR);
	 BSplineUniLR *uninew = new BSplineUniLR(origuni->pardir(), origuni->degree(),
						 k2_v, origuni->getMesh());

	 // Insert in array
	 BSplineUniUtils::insert_univariate(bspline_vec2, uninew, ix);
       }
      uni2_2 = bspline_vec2[ix].get();
      uni1_1 = uni1_2 = orig.getUnivariate(XDIR);
      uni3_1 = uni3_2 = orig.getUnivariate(ZDIR);
    }
  else
    {
      int ix = 0;
      bool found1 = 
	BSplineUniUtils::identify_bsplineuni(k1_w, k1_w+kvec_size, bspline_vec3, ix);
     if (!found1)
       {
	 // Must create univariate B-spline
	 BSplineUniLR *origuni = orig.getUnivariate(ZDIR);
	 BSplineUniLR *uninew = new BSplineUniLR(origuni->pardir(), origuni->degree(),
						 k1_w, origuni->getMesh());

	 // Insert in array
	 BSplineUniUtils::insert_univariate(bspline_vec3, uninew, ix);
       }
       uni3_1 = bspline_vec3[ix].get();
      bool found2 = 
	BSplineUniUtils::identify_bsplineuni(k2_w, k2_w+kvec_size, bspline_vec3, ix);
      if (!found2)
       {
	 // Must create univariate B-spline
	 BSplineUniLR *origuni = orig.getUnivariate(ZDIR);
	 BSplineUniLR *uninew = new BSplineUniLR(origuni->pardir(), origuni->degree(),
						 k2_w, origuni->getMesh());

	 // Insert in array
	 BSplineUniUtils::insert_univariate(bspline_vec3, uninew, ix);
       }
      uni3_2 = bspline_vec3[ix].get();
      uni1_1 = uni1_2 = orig.getUnivariate(XDIR);
      uni2_1 = uni2_2 = orig.getUnivariate(YDIR);
    }
  
  new_1 = new LRBSpline3D(c_g1, w1, uni1_1, uni2_1, uni3_1, g1, rat);
  //new_1->setFixCoef(orig.coefFixed());

  new_2 = new LRBSpline3D(c_g2, w2, uni1_2, uni2_2, uni3_2, g2, rat);
  //new_2->setFixCoef(orig.coefFixed());

}


//==============================================================================
// if 'b' can be split at least once in the mesh 'm', split it once, and return the 
// result through 'b1' and 'b2'.  The function never carries out more than one split, 
// even when several splits are possible.
bool LRBSpline3DUtils::try_split_once(const LRBSpline3D& b, const Mesh3D& mesh, 
				      int mult1, int mult2, int mult3,
				      vector<unique_ptr<BSplineUniLR> >& bspline_vec1,
				      vector<unique_ptr<BSplineUniLR> >& bspline_vec2,
				      vector<unique_ptr<BSplineUniLR> >& bspline_vec3,
                                      LRBSpline3D*& b1,
                                      LRBSpline3D*& b2)
//==============================================================================
{
  const int umin = b.suppMin(XDIR);
  const int vmin = b.suppMin(YDIR);
  const int wmin = b.suppMin(ZDIR);
  const int umax = b.suppMax(XDIR);
  const int vmax = b.suppMax(YDIR);
  const int wmax = b.suppMax(ZDIR);
  const int xmult = b.endmult_u(true) + b.endmult_u(false);
  const int ymult = b.endmult_v(true) + b.endmult_v(false);
  const int zmult = b.endmult_w(true) + b.endmult_w(false);

  // @@ The assertions below should always hold if function is called with correct 
  // argument. When code is properly debugged and tested, they can be taken away for
  // efficiency (asserts should go away anyway when compiling in optimized mode).
  // Alternatively, if it is a concern that users might call this function with wrong 
  // argument, the assertions could be replaced by exception-throwing 'if'-statements.

  if (umax - umin + xmult + mult1 - 3 > b.degree(XDIR)+1)
    {
      const vector<int> m_kvec_u = 
	derive_knots(mesh, XDIR, umin, umax, vmin, vmax, wmin, wmax);
      if (!std::includes(m_kvec_u.begin(), m_kvec_u.end(), 
			 b.kvec(XDIR).begin(), b.kvec(XDIR).end()))
	THROW("B-spline knot vector not included in mesh knot vector, XDIR");
      
      if (num_inner_knots(m_kvec_u) > num_inner_knots(b.kvec(XDIR))) {
	// Since we know that m_kvec_u contains more elements than b.kvec(XDIR) and since
	// we know that the latter is included in the former, we know that there must be at least
	// one knot in 'm_kvec_u' that is not found in b.kvec(XDIR).  We can therefore call
	// the following function without risking an exception to be thrown.
	const int new_ix = find_uncovered_inner_knot(m_kvec_u, b.kvec(XDIR));

	// @@@ VSK. Cannot set pointers to the new bsplines before it is placed in the
	// global array. Will the position of the element change when a bspline is removed?
	split_function(b, XDIR, mesh.knotsBegin(XDIR), new_ix, bspline_vec1,
		       bspline_vec2, bspline_vec3, b1, b2);
	return true;
      }
    }

  if (vmax - vmin + ymult + mult2 - 3 > b.degree(YDIR)+1)
    {
      const vector<int> m_kvec_v = 
	derive_knots(mesh, YDIR, vmin, vmax, wmin, wmax, umin, umax);
      if (!std::includes(m_kvec_v.begin(), m_kvec_v.end(), 
			   b.kvec(YDIR).begin(), b.kvec(YDIR).end()))
	THROW("B-spline knot vector not included in mesh knot vector, YDIR");

      if (num_inner_knots(m_kvec_v) > num_inner_knots(b.kvec(YDIR))) {
	// same comment as above
	const int new_ix = find_uncovered_inner_knot(m_kvec_v, b.kvec(YDIR));
	split_function(b, YDIR, mesh.knotsBegin(YDIR), new_ix, bspline_vec1,
		       bspline_vec2, bspline_vec3, b1, b2);
	return true;
      }
    }

  if (wmax - wmin + zmult + mult3 -3 > b.degree(ZDIR)+1)
    {
      const vector<int> m_kvec_w = 
	derive_knots(mesh, ZDIR, wmin, wmax, umin, umax, vmin, vmax);
      if (!std::includes(m_kvec_w.begin(), m_kvec_w.end(), 
			   b.kvec(ZDIR).begin(), b.kvec(ZDIR).end()))
	THROW("B-spline knot vector not included in mesh knot vector, ZDIR");

      if (num_inner_knots(m_kvec_w) > num_inner_knots(b.kvec(ZDIR))) {
	// same comment as above
	const int new_ix = find_uncovered_inner_knot(m_kvec_w, b.kvec(ZDIR));
	split_function(b, ZDIR, mesh.knotsBegin(ZDIR), new_ix, bspline_vec1,
		       bspline_vec2, bspline_vec3, b1, b2);
	return true;
      } 
    }

  // No splits possible
  return false;
}


} // end namespace Go
