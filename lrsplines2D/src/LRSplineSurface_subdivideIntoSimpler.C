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

#include <iostream> // debug
#include <algorithm>
#include <chrono> // profiling
#include <numeric>
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/geometry/CurveBoundedDomain.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h" // debug
#include "GoTools/lrsplines2D/LRSplinePlotUtils2.h" // debug
#include <fstream>

using std::vector;
using std::map;
using std::unique_ptr;
using std::pair;
using std::array;
using std::tuple;
using std::make_pair;
using std::max;
using std::min;
using std::get;
using std::cout; //debug

using namespace Go;

using IntPair = pair<int, int>;
using IntVec  = vector<int>;
using Int3Vec = vector<array<int, 3>>;

//#define DEBUG
//#define VERBOSE

namespace { // anonymous namespace

  // Structure to denote a subdomain
  struct Subdomain {
    Subdomain(const Mesh2D& m) :
      range_x({0, m.numDistinctKnots(XFIXED)-1}),
      range_y({0, m.numDistinctKnots(YFIXED)-1}) {};
    Subdomain(IntPair r1, IntPair r2) : range_x(r1), range_y(r2) {}
    
    IntPair range_x; // first and last knot along x-knotvec (ignoring multiplicities)
    IntPair range_y; // first and last knot along y-knotvec (ignoring multiplicities)
  };

    // determine the number of segments/multiplicities that are missing in order
  // to convert the specified subdomain to a full grid.  Does not consider the
  // boundary segments (first and last knot of each range) Vector contains one
  // entry per line considered.  First value is the number of missing
  // segments/multiplicities.  Second is the maximum muliplicity in the range
  // under consideration.  Third value is the number of individual segments in
  // the range under consideration.
  Int3Vec missing_tensorgrid_segments(const Mesh2D& m,
				      const Direction2D d,
				      const Subdomain& dom);

  // count total number of missing segments/multiplicities to make the specified
  // subdomain into a full grid (not counting boundaries).
  struct MissingSegInfo {pair<int,int> num; Int3Vec missing_x; Int3Vec missing_y;};
  MissingSegInfo total_missing_tensorgrid_segments(const Mesh2D&m,
						   const Subdomain& dom);



  // Suggest a split of the specified domain into two subdomains with fewer
  // total missing segments/multiplicities.  'bnd_mult' specifies the
  // multiplicity each segment of the split itself should have.
  // Return values:
  // 1) Direction (whether the split occurs for XFIXED or YFIXED)
  // 2) Index of the line along which the split should occur (or -1 if no split
  //    was found)
  // 3) The total missing tensorgrid segments before split
  // 4) The cost of the split (the number of minisegments multiplicities to
  //    raise along the split line)
  // 5) The total missing tensorgrid segments of the first new subdomain
  // 6) The total missing tensorgrid segments of the second new subdomain
  // If a split is found, then the sum of the three last returned values should
  // be strictly lower than the value of the third returned value.  If this is
  // not possible, no split is proposed, and the second returned value will be -1.
  //tuple<Direction2D, int, int, int, int, int>
  struct DomainSplitDetail {
    Direction2D d;
    int ix;
    int init_cost;
    int split_cost;
    int new_cost1;
    int new_cost2;

    DomainSplitDetail()
    {
      ix = -1;  // No split identified
      d = XFIXED; // Dummy
      init_cost = 0;
      split_cost = 0;
      new_cost1 = 0;
      new_cost2 = 0;
    }

    DomainSplitDetail(Direction2D dir, int x, int init, int split,
		      int cost1, int cost2)
    {
      d = dir;
      ix = x;
      init_cost = init;
      split_cost = split;
      new_cost1 = cost1;
      new_cost2 = cost2;
    }
  };
  DomainSplitDetail
  suggest_domain_split(const Mesh2D& m, const Subdomain& dom, const IntPair bnd_mult,
		       const int threshold_missing);

  struct ConsecutiveSplit {Direction2D d; int ix; IntPair range;};
  tuple<vector<ConsecutiveSplit>, vector<pair<Subdomain,LRSplineSurface::PatchStatus>>>
    recursive_split(const Mesh2D& m, const IntPair bnd_mult, const int threshold_missing,
		    const Subdomain dom, const CurveBoundedDomain* bddomain, 
		    double tol);

    LRSplineSurface::PatchStatus
    check_domain_trim_status(const Mesh2D& mesh, const Subdomain dom, 
			     const CurveBoundedDomain* bddomain, double tol);

  vector<LRSplineSurface::Refinement2D>
  prepare_refinements(const Mesh2D& m,
		      const vector<ConsecutiveSplit>& splits,
		      const IntPair mult);


  Mesh2D remove_unused_knots(Mesh2D& m);

  vector<int> reindex_knots(const vector<int>& ixs_old,
			    const vector<double>& kvals_new,
			    const vector<double>& kvals_old);

  unique_ptr<LRBSpline2D> 
  adapt_bspline(const LRBSpline2D* const b,
		LRSplineSurface& patch,
		std::vector<std::unique_ptr<BSplineUniLR> >& bsplinesuni1,
		std::vector<std::unique_ptr<BSplineUniLR> >& bsplinesuni2,
		const LRSplineSurface& orig_surf);
}; //end anonymous namespace

namespace Go
{
  // --------------------------------------------------------------------------
  IntVec active_knots(const Mesh2D& m, Direction2D d, IntPair r1, IntPair r2)
  // --------------------------------------------------------------------------
  {
    // The first and last knot is automatically 'active' by merit to belong to
    // the boundary of the domain under consideration
    
    if (d==YFIXED) swap(r1, r2); // ensure that r1 correspond to the fixed direction
    IntVec result {r1.first}; // include first knot (see initial comment)

    // determine which of the "interior" knots should be considered active
    for (int i = r1.first+1; i < r1.second; ++i) {
      const vector<IntPair> segs = m.segments(d, i);
      if (any_of(segs.begin(), segs.end(),
		 [r2] (IntPair p) {return (p.first < r2.second) &
		                          (p.second > r2.first);}))
	result.push_back(i);
    }

    result.push_back(r1.second); // include last knot (see inital comment)

    return result;
  }

  // --------------------------------------------------------------------------
  IntVec make_segmult_map(const vector<GPos>& segs, int start_ix, int end_ix)
  // --------------------------------------------------------------------------
  {
    // determine the knot multiplicities of minisegments starting with a given knot
    IntVec result(end_ix - start_ix, segs[0].mult);

    for (size_t i = 1; (i != segs.size()) && (segs[i].ix < end_ix); ++i) // @@ can probably be optimized
      fill(result.begin() + max(0, (segs[i].ix - start_ix)), result.end(), segs[i].mult);
    return result;
  }

  
  // --------------------------------------------------------------------------
  IntVec inisegment_multiplicities(const Mesh2D& m,
					 const Direction2D d,
					 const int seg_ix,
					 const IntVec& knots)
  // --------------------------------------------------------------------------
  {
    assert(knots.size() >= 2);

    const auto& all_segs = m.mrects(d, seg_ix); // all segments

    // segments intersecting with the intervals within the given knots

    // making map of multiplicities
    const IntVec segmult = make_segmult_map(all_segs, knots.front(), knots.back());

    // computing the multiplicities of each segment between two consecutive knots

    IntVec result(knots.size()-1);
    const auto end_it = knots.end()-1;
    const int first_ix = knots.front();
    transform(knots.begin(), end_it, result.begin(),
    	      [&] (int k) {return segmult[k - first_ix];});

    // IntVec result;
    // result.reserve(knots.size()-1);
    // const auto end_it = knots.end()-1;
    // const int first_ix = knots.front();
    // transform(knots.begin(), end_it, back_inserter(result),
    // 	      [&] (int k) {return segmult[k - first_ix];});

    
    // IntVec result(knots.size()-1);
    // auto target = result.begin();
    // int first_ix = knots.front();
    // const auto last_it = knots.end()-1;
    // for (auto it = knots.begin(); it != last_it; ++it) {
    //   *target++ = segmult[*it - first_ix];
    // }
    
    return result;
  }

  
  // ============================================================================
  vector<pair<shared_ptr<LRSplineSurface>, LRSplineSurface::PatchStatus> > 
  LRSplineSurface::subdivideIntoSimpler(const int threshold_missing, 
					double tol,
					const CurveBoundedDomain* domain) const
  // ============================================================================
  {
#ifdef VERBOSE
    std::cout << "Entering subdivision function." << std::endl;
#endif

    // Determine subdivisions
    const IntPair order {degree(XFIXED)+1, degree(YFIXED)+1};
    auto t1 = std::chrono::high_resolution_clock::now();
    const auto splits = recursive_split(mesh(), order, threshold_missing, 
					Subdomain(mesh()), domain, tol);
    auto t2 = std::chrono::high_resolution_clock::now();
#ifdef VERBOSE
    std::cout << "identifying splits took "
	       << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
	       << " milliseconds\n";
    std::cout << get<0>(splits).size() << " splits identified." << std::endl;
#endif
#ifdef DEBUG
    vector<pair<Subdomain,LRSplineSurface::PatchStatus> > tmp = get<1>(splits);
    size_t nn = tmp.size();
    for (size_t ki=0; ki<nn;  ki+=20)
      {
	for (int ka=0; ka<20 && ki+ka<nn; ++ka)
	  std::cout << tmp[ki+ka].second << "  ";
	std::cout << std::endl;
      }
#endif
    // Making working copy of the present surface, and carry out the determined
    // splits by raising internal multiplicities accordingly.
    auto lrs_copy = shared_ptr<LRSplineSurface>(this->clone());

    t1 = std::chrono::high_resolution_clock::now();
    lrs_copy->refine(prepare_refinements(mesh(), get<0>(splits), order), true);

    t2 = std::chrono::high_resolution_clock::now();
#ifdef VERBOSE
    std::cout << "refining  took "
	       << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
	       << " milliseconds\n";
    std::cout << "Finished refining." << std::endl;
    // Extract and return individual surface patches
    std::cout << "Generating " << get<1>(splits).size() << " new surfaces..." << std::endl;
#endif
    vector<pair<shared_ptr<LRSplineSurface>, LRSplineSurface::PatchStatus>> result;
    map<ElemKey, size_t> patchmap;


    t1 = std::chrono::high_resolution_clock::now();    
    // Establishing individual surface patches, empty, without Bspline-functions
    for (const auto patch : get<1>(splits)) {
      
      const auto submesh =
      	remove_unused_knots(*(lrs_copy->mesh().subMesh(patch.first.range_x.first,
       						       patch.first.range_x.second,
       						       patch.first.range_y.first,
       						       patch.first.range_y.second)));
      // const auto submesh = lrs_copy->mesh().subMesh(patch.range_x.first,
      // 						    patch.range_x.second,
      // 						    patch.range_y.first,
      // 						    patch.range_y.second);

      //plot_mesh(submesh);

      auto cur = shared_ptr<LRSplineSurface>(new LRSplineSurface());

      cur->knot_tol_ = lrs_copy->knot_tol_;
      cur->rational_ = lrs_copy->rational_;
      cur->mesh_     = submesh;

      // etablishing the element map
      cur->emap_ = LRSplineUtils::identify_elements_from_mesh(cur->mesh_);
      result.emplace_back(make_pair(cur,patch.second));

      // making the current surface easy to find when distributing the
      // Bspline-functions later
      const size_t ix = result.size()-1;
      for (auto e = cur->elementsBegin(); e != cur->elementsEnd(); ++e) {
      	assert(e->first.u_min < cur->mesh_.maxParam(XFIXED));
      	assert(e->first.v_min < cur->mesh_.maxParam(YFIXED));
      	patchmap[e->first] = ix;
      }
    }

    // Distribute Bsplines.  If surface has been correctly subdivided, each
    // Bspline function belongs to exactly one patch.
    for (auto b = lrs_copy->basisFunctionsBegin();
    	 b != lrs_copy->basisFunctionsEnd(); ++b) {
      const BSKey key = b->first;
      const size_t patch_ix = patchmap[ElemKey {key.u_min, key.v_min}];
      auto patch = result[patch_ix].first;

      // make copy of current BSpline basis function, and add it to the patch
      patch->bsplines_[key] = 
	adapt_bspline(b->second.get(), *patch, patch->bsplinesuni1_,
		      patch->bsplinesuni2_, *lrs_copy);

      auto krull1 = patch->bsplines_[key]->kvec(XFIXED); //@@@
      // auto p1 = krull1.front();
      // auto p2 = krull1.back();
      // if (p1 == p2)
      // 	cout << p1 << " " << p2 << patch_ix << std::endl;
      // update internal links
      LRSplineUtils::update_elements_with_single_bspline(patch->bsplines_[key].get(),
    							 patch->emap_,
    							 patch->mesh(), 
							 false);
    }
    t2 = std::chrono::high_resolution_clock::now();
#ifdef VERBOSE

    std::cout << "generating patches  took "
	       << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
	       << " milliseconds\n";
    std::cout << "Finished generating patches." << std::endl;
#endif
    // vector<pair<shared_ptr<LRSplineSurface>, LRSplineSurface::PatchStatus> > result2(result.size());
    // for (size_t ki=0; ki<result.size(); ++ki)
    //   result2[ki] = make_pair(result[ki], INSIDE);

    return result;
  }

    

  
}; // end namespace Go


namespace { // anonymous namespace

  // ----------------------------------------------------------------------------
  vector<int> reindex_knots(const vector<int>& ixs_old,
			    const vector<double>& kvals_new,
			    const vector<double>& kvals_old)
  // ----------------------------------------------------------------------------
  {
    vector<double> kvals; kvals.reserve(ixs_old.size());
    transform(ixs_old.begin(), ixs_old.end(), back_inserter(kvals),
	      [&kvals_old] (int ix) {return kvals_old[ix];});

    vector<int> result; result.reserve(ixs_old.size());
    transform(kvals.begin(), kvals.end(), back_inserter(result),
	      [&kvals_new] (double val) {
		return int(find(kvals_new.begin(), kvals_new.end(), val) -
			   kvals_new.begin());});
    return result;
  }
  
  // ----------------------------------------------------------------------------
  unique_ptr<LRBSpline2D> adapt_bspline(const LRBSpline2D* const b,
					LRSplineSurface& patch,
					vector<unique_ptr<BSplineUniLR> >& bsplinesuni1,
					vector<unique_ptr<BSplineUniLR> >& bsplinesuni2,
					const LRSplineSurface& orig_surf)
  // ----------------------------------------------------------------------------
  {
    // make a copy of the original bspline function
    auto result = unique_ptr<LRBSpline2D>(new LRBSpline2D(*b));

    // clear old support
    result->setSupport(vector<Element2D*> {}); 
    
    // determine knot indices according to new spline mesh
    vector<int> kvec_u = 
      reindex_knots(result->kvec(XFIXED),
		    vector<double>(patch.mesh().knotsBegin(XFIXED),
				   patch.mesh().knotsEnd(XFIXED)),
		    vector<double>(orig_surf.mesh().knotsBegin(XFIXED),
				   orig_surf.mesh().knotsEnd(XFIXED)));

    // Check if the univariate B-spline exists already
    int left1 = ((int)bsplinesuni1.size())/2;
    bool found1 = 
      BSplineUniUtils::identify_bsplineuni(kvec_u.begin(), kvec_u.end(),
					   bsplinesuni1, left1);
    if (!found1)
      {
	BSplineUniLR *uni1 = new BSplineUniLR(1, b->degree(XFIXED), kvec_u.begin(), 
					      &patch.mesh());
	BSplineUniUtils::insert_univariate(bsplinesuni1, uni1, left1);
      }
    result->setUnivariate(XFIXED, bsplinesuni1[left1].get());

    vector<int> kvec_v =
      reindex_knots(result->kvec(YFIXED),
		    vector<double>(patch.mesh().knotsBegin(YFIXED),
				   patch.mesh().knotsEnd(YFIXED)),
		    vector<double>(orig_surf.mesh().knotsBegin(YFIXED),
				   orig_surf.mesh().knotsEnd(YFIXED)));

    int left2 = ((int)bsplinesuni2.size())/2;
    bool found2 = 
      BSplineUniUtils::identify_bsplineuni(kvec_v.begin(), kvec_v.end(),
					   bsplinesuni2, left2);
    if (!found2)
      {
	BSplineUniLR *uni2 = new BSplineUniLR(2, b->degree(YFIXED), kvec_v.begin(), 
					      &patch.mesh());
	BSplineUniUtils::insert_univariate(bsplinesuni2, uni2, left2);
      }
    result->setUnivariate(YFIXED, bsplinesuni2[left2].get());

    return result;
  }

  // ==========================================================================
  Int3Vec missing_tensorgrid_segments(const Mesh2D& m,
				      const Direction2D d,
				      const Subdomain& dom)
  // ==========================================================================
  {
    // determine active knots in the opposite direction (not all zero in the segment)
    const auto d2_knots = active_knots(m, flip(d), dom.range_x, dom.range_y); 

    const IntPair r1 = (d == XFIXED) ? dom.range_x : dom.range_y; // d-direction
    const IntPair r2 = (d == XFIXED) ? dom.range_y : dom.range_x;  // other direction

    Int3Vec result;
    result.reserve(max(r1.second - r1.first - 1, 1));
    
    // we only consider interior.  Ignore first and last value of range
    for (int i = r1.first+1; i != r1.second; ++i) {

      // multiplicities for each "mini-segment" along this line, within the domain
      const auto ms_mult = inisegment_multiplicities(m, d, i, d2_knots);

      const int max_mult = *max_element(ms_mult.begin(), ms_mult.end());

      // compute number of missing minisegments, and add it to result vector
      const int missing = max_mult * (int)ms_mult.size() -
	                  std::accumulate(ms_mult.begin(), ms_mult.end(), 0);
      result.push_back( {missing, max_mult, int(ms_mult.size())});
    }
      
    return result;
  }

  // ==========================================================================
  MissingSegInfo
  total_missing_tensorgrid_segments(const Mesh2D& m, const Subdomain& dom)
  // ==========================================================================
  {
    const auto c1 = missing_tensorgrid_segments(m, XFIXED, dom);
    const auto c2 = missing_tensorgrid_segments(m, YFIXED, dom);

    const auto add_fun = [](int cur, const array<int, 3>& a) {return cur + a[0];};

    int missing_x = std::accumulate(c1.begin(), c1.end(), 0, add_fun);
    int missing_y = std::accumulate(c2.begin(), c2.end(), 0, add_fun);
    
    pair<int,int> num_missing = make_pair(missing_x, missing_y);
    return {num_missing, c1, c2};
  }

  // --------------------------------------------------------------------------
  vector<IntPair>
  nmb_t_junctions(const Mesh2D& m, Direction2D d, IntVec& cand, 
		  IntPair range, IntPair range_other)
  // --------------------------------------------------------------------------
  {
    IntPair init = make_pair(0, range.second-range.first);
    vector<IntPair> nmb_T(cand.size(), init);
    for (size_t kj=0; kj<cand.size(); ++kj)
      {
	// skip first and last value of range
	for (int i = range_other.first+1; i < range_other.second; ++i) 
	  {
	    const auto mrects = m.mrects(flip(d), i);
	    for (size_t kh=0; kh<mrects.size(); ++kh)
	      if (mrects[kh].ix == cand[kj])
		nmb_T[kj].first++;
	  }
	nmb_T[kj].second = std::min(range.second-cand[kj],
				    cand[kj]-range.first);
      }

    return nmb_T;
  }

  // --------------------------------------------------------------------------
  bool has_t_junctions(const Mesh2D& m, Direction2D d, int line_ix, IntPair range)
  // --------------------------------------------------------------------------
  {
    // skip first and last value of range
    for (int i = range.first+1; i < range.second; ++i) {
      const auto mrects = m.mrects(flip(d), i);
      if (any_of(mrects.begin(), mrects.end(), [line_ix](GPos p) {return p.ix == line_ix;}))
	return true;
    }
    return false;
  }

  // --------------------------------------------------------------------------
  int sort_splitting_candidates(IntVec& cand, vector<IntPair>& nmb_t)
  // --------------------------------------------------------------------------
  {
    // Sort candidates with respect to the number of T-joints
    for (size_t ki=0; ki<cand.size(); ++ki)
      for (size_t kj=ki+1; kj<cand.size(); ++kj)
	{
	  if (nmb_t[kj].first > nmb_t[ki].first)
	    {
	      std::swap(cand[ki], cand[kj]);
	      std::swap(nmb_t[ki], nmb_t[kj]);
	    }
	}

    int nmb = 0;
    int level = 2*nmb_t[0].first/3;
    int max_range_dist = nmb_t[0].second;
    for (; nmb<(int)nmb_t.size(); ++nmb)
      {
	max_range_dist = std::max(max_range_dist, nmb_t[nmb].second);
	if (nmb_t[nmb].first < level)
	  break;
      }

    // Sort the remaining candidates with respect to closeness to range centre
    for (int kr=nmb; kr<(int)cand.size(); ++kr)
      for (int kh=kr+1; kh<(int)cand.size(); ++kh)
	{
	  if (nmb_t[kh].second > nmb_t[kr].second)
	    {
	      std::swap(cand[kr], cand[kh]);
	      std::swap(nmb_t[kr], nmb_t[kh]);
	    }
	}

    int nmb2 = nmb;
    for (; nmb2<(int)nmb_t.size(); ++nmb2)
      if (nmb_t[nmb2].second < max_range_dist)
	break; 

    return nmb2;
  }

  struct PerfInfo {int ix; int nmb_t; int perf; int range_distr;
    int split_cost; int new_cost1; int new_cost2;};
  // --------------------------------------------------------------------------
  PerfInfo choose_splitting_candidate(const Mesh2D& m, const Subdomain& dom,
				      const MissingSegInfo& init_missing, 
				      Direction2D d, int bnd_mult)
  // --------------------------------------------------------------------------
  {
    const Int3Vec& missing_d  = (d == XFIXED) ? init_missing.missing_x :
                                                init_missing.missing_y;
    // const auto adder = [](int cur, const array<int, 3>& a) {return cur + a[0];};
    // const int num_missing = std::accumulate(missing_d.begin(), missing_d.end(), 0, adder);

    // if (num_missing == 0)
    //   return{-1, 0, 0, 0, 0}; // no split necessary. 

    const IntPair range_d     = (d == XFIXED) ? dom.range_x : dom.range_y;
    const IntPair range_other = (d == XFIXED) ? dom.range_y : dom.range_x;
    double scale = 250.0/((range_d.second-range_d.first)*(range_other.second-range_other.first));
    scale = std::min(scale, 0.1);
    
    IntVec cand(range_d.second - (range_d.first + 1), 0);
    iota(cand.begin(), cand.end(), range_d.first + 1); // consecutive values

    cand.erase(remove_if(cand.begin(), cand.end(),[&m, d, range_other](int ix)
			 {return !has_t_junctions(m, d, ix, range_other);}),
	       cand.end());
    if (cand.size() == 0) // no possible split
      return {-1, 0, 0, 0, 0, 0, 0};

    vector<IntPair> nmb_t = 
      nmb_t_junctions(m, d, cand, range_d, range_other);

    int nmb_significant = sort_splitting_candidates(cand, nmb_t);

    // search for optimal candidate
    vector<PerfInfo> perf; // assess performance of each candidate
    //transform(cand.begin(), cand.end(), back_inserter(perf), [&] (int ix) {
    //for (size_t ki=0; ki<cand.size(); ++ki)
    for (int ki=0; ki<nmb_significant; ++ki)
      {
    	int ix = cand[ki];
	const int local_ix = ix - (range_d.first + 1);
	const int cur_missing = missing_d[local_ix][0];
	const int cur_maxmult = missing_d[local_ix][1]; assert(cur_maxmult <= bnd_mult);
	const int num_miniseg = missing_d[local_ix][2];
	
	// cost to bring this line subdomain up to 'bnd_mult' multiplicity everywhere
	const int split_cost = cur_missing + (bnd_mult - cur_maxmult) * num_miniseg;
	const IntPair r1 {range_d.first, ix};
	const IntPair r2 {ix, range_d.second};
	pair<int,int> cost1 =
	  total_missing_tensorgrid_segments(m,
					    (d==XFIXED) ? Subdomain {r1, range_other} : Subdomain {range_other, r1}).num;
	const int new_cost1 = cost1.first + cost1.second;
	pair<int,int> cost2 = 
	  total_missing_tensorgrid_segments(m,
					    (d==XFIXED) ? Subdomain {r2, range_other} : Subdomain {range_other, r2}).num;
	const int new_cost2 = cost2.first + cost2.second;

	const int distr = (d==XFIXED) ? 
	  cost1.first+cost2.first+std::min(cost1.first,cost2.first) :
	  cost1.second+cost2.second+std::min(cost1.second,cost2.second); 
	// return PerfInfo {ix,
	//                  init_missing.num - (split_cost + new_cost1 + new_cost2),
	//                  split_cost, new_cost1, new_cost2};
	int num = init_missing.num.first + init_missing.num.second;
	int curr_perf = num - (split_cost + new_cost1 + new_cost2);
	if (curr_perf > 0)
	  curr_perf *= (int)(scale*(nmb_t[ki].first+nmb_t[ki].second));
      	perf.push_back(PerfInfo {ix, nmb_t[ki].first, curr_perf,
      	      /*num - (split_cost + new_cost1 + new_cost2),*/
      	      distr, split_cost, new_cost1, new_cost2});
      }
    //         });

    PerfInfo sel_perf =  *max_element(perf.begin(), perf.end(),
				      [](const PerfInfo& p1, const PerfInfo& p2)
				      {return p1.perf < p2.perf;});
    return sel_perf;
  }
  
  // ==========================================================================
  DomainSplitDetail suggest_domain_split(const Mesh2D& m,
					 const Subdomain& dom,
					 const IntPair bnd_mult,
					 const int threshold_missing)
  // ==========================================================================
  {
    DomainSplitDetail dummy;
    // Unless there are internal lines in both directions, there is no point in splitting
    if (min(dom.range_x.second - dom.range_x.first, dom.range_y.second - dom.range_y.first) < 2)
      return dummy;
      //return {XFIXED, -1, 0, 0, 0, 0}; // dummy values, except -1, which flags 'no split'
      
    // If we got here, we know that both direction has at least one internal line
    const auto init = total_missing_tensorgrid_segments(m, dom);
    int num = init.num.first+init.num.second;
    if (num < threshold_missing)
      return DomainSplitDetail(XFIXED, -1, num, 0, 0, 0);

    const PerfInfo best_split_x = choose_splitting_candidate(m, dom, init, XFIXED, bnd_mult.first);
    const PerfInfo best_split_y = choose_splitting_candidate(m, dom, init, YFIXED, bnd_mult.second);
    const Direction2D best_dir = best_split_x.ix < 0                     ? YFIXED :
                                 best_split_y.ix < 0                     ? XFIXED :
                                 (best_split_x.perf > best_split_y.perf) ? XFIXED : YFIXED;
    const PerfInfo split = (best_dir == XFIXED) ? best_split_x : best_split_y;
    
    const int split_ix = (split.perf > 0) ? split.ix :
      (split.split_cost <= max(split.new_cost1, split.new_cost2)) ? split.ix : -1;

    return DomainSplitDetail(best_dir, split_ix, num, 
			     split.split_cost, split.new_cost1, split.new_cost2);
  }

  // ==========================================================================
  tuple<vector<ConsecutiveSplit>, vector<pair<Subdomain,LRSplineSurface::PatchStatus>>>
    recursive_split(const Mesh2D& m, const IntPair bnd_mult, const int threshold_missing,
		    const Subdomain dom, const CurveBoundedDomain* bddomain,
		    double tol)
  // ==========================================================================
  {
#ifdef DEBUG
    std::ofstream mesh0("mesh_init.eps");
    writePostscriptMesh(m, dom.range_x.first, dom.range_x.second,
			dom.range_y.first, dom.range_y.second, mesh0, true);
#endif

    LRSplineSurface::PatchStatus patchtype = LRSplineSurface::INSIDE;
    if (bddomain)
      {
	// Check if the current subdomain is outside the trimming loop.
	// In that case, stop the recursive split
	patchtype = check_domain_trim_status(m, dom, bddomain, tol);
      }

    DomainSplitDetail cur;
    cur = suggest_domain_split(m, dom, bnd_mult, threshold_missing);

    if (cur.ix < 0 || patchtype == LRSplineSurface::OUTSIDE) 
      {
	// no advantageous split possible or subdomain outside trimming loop
	pair<Subdomain, LRSplineSurface::PatchStatus> subdom = 
	  make_pair(dom, patchtype);
	return make_tuple(vector<ConsecutiveSplit>{}, 
			  vector<pair<Subdomain,LRSplineSurface::PatchStatus>>{subdom});
      }

    // a split was found
    const Direction2D d = cur.d; // avoid having to spell it out every time
    const pair<Subdomain, Subdomain> new_doms = (d == XFIXED) ?
      make_pair(Subdomain {{dom.range_x.first, cur.ix}, dom.range_y},
		Subdomain {{cur.ix, dom.range_x.second}, dom.range_y}) : 
      make_pair(Subdomain {dom.range_x, {dom.range_y.first, cur.ix}},
		Subdomain {dom.range_x, {cur.ix, dom.range_y.second}});
#ifdef DEBUG
    std::ofstream mesh2("mesh_1.eps");
    std::ofstream mesh3("mesh_2.eps");
    writePostscriptMesh(m, new_doms.first.range_x.first, new_doms.first.range_x.second,
      new_doms.first.range_y.first, new_doms.first.range_y.second, mesh2, true);

    writePostscriptMesh(m, new_doms.second.range_x.first, new_doms.second.range_x.second,
      new_doms.second.range_y.first, new_doms.second.range_y.second, mesh3, true);
#endif
    const auto r1 = recursive_split(m, bnd_mult, threshold_missing, new_doms.first, 
				    (patchtype == LRSplineSurface::INTERSECT) ?
				    bddomain : NULL, tol);
    const auto r2 = recursive_split(m, bnd_mult, threshold_missing, new_doms.second, 
				    (patchtype == LRSplineSurface::INTERSECT) ?
				    bddomain : NULL, tol);

    vector<ConsecutiveSplit> splits {{d, cur.ix, (d == XFIXED) ? dom.range_y : dom.range_x} };
    splits.insert(splits.end(), get<0>(r1).begin(), get<0>(r1).end());
    splits.insert(splits.end(), get<0>(r2).begin(), get<0>(r2).end());
    
    vector<pair<Subdomain,LRSplineSurface::PatchStatus>> doms = get<1>(r1);
    doms.insert(doms.end(), get<1>(r2).begin(), get<1>(r2).end());

    return make_tuple(splits, doms);
  }

  // ==========================================================================
    LRSplineSurface::PatchStatus
    check_domain_trim_status(const Mesh2D& mesh, const Subdomain dom, 
			     const CurveBoundedDomain* bddomain, 
			     double tol)
  // ==========================================================================
    {
      if (bddomain == NULL)
	return LRSplineSurface::INSIDE;

      // Define subdomain boundaries as 2D spline curves
      double x1 = mesh.kval(XFIXED, dom.range_x.first);
      double x2 = mesh.kval(XFIXED, dom.range_x.second);
      double y1 = mesh.kval(YFIXED, dom.range_y.first);
      double y2 = mesh.kval(YFIXED, dom.range_y.second);

      vector<Point> corner(4);
      corner[0] = Point(x1, y1);
      corner[1] = Point(x2, y1);
      corner[2] = Point(x2, y2);
      corner[3] = Point(x1, y2);

      // First check for intersection
      vector<shared_ptr<ParamCurve> > bd_cvs(4);
      for (size_t ki=0; ki<corner.size(); ++ki)
	{
	  shared_ptr<SplineCurve> bd(new SplineCurve(corner[ki], 
						     corner[(ki+1)%4]));
	  bd_cvs[ki] = bd;
	  if (bddomain->doIntersect(*bd, tol))
	    return LRSplineSurface::INTERSECT;
	}

      // Check inside/outside status
      Vector2D mid(0.5*(x1+x2), 0.5*(y1+y2));
      if (bddomain->isInDomain(mid, tol))
	return LRSplineSurface::INSIDE;
      else
	{
	  // Check if the trimming domain is inside the subdomain
	  shared_ptr<CurveLoop> loop(new CurveLoop(bd_cvs, tol, false));
	  CurveBoundedDomain dom2(loop);
	  double upar, vpar;
	  bddomain->getInternalPoint(upar, vpar);
	  Vector2D mid2(upar, vpar);
	  if (dom2.isInDomain(mid2, tol))
	    return LRSplineSurface::INTERSECT;
	  else
	    return LRSplineSurface::OUTSIDE;
	}
    }

  // ==========================================================================
  vector<LRSplineSurface::Refinement2D> 
  prepare_refinements(const Mesh2D& m, const vector<ConsecutiveSplit>& splits,
		      const IntPair mult)
  // ==========================================================================
  {
    vector<LRSplineSurface::Refinement2D> result;
    transform(splits.begin(), splits.end(), back_inserter(result), [&] (const ConsecutiveSplit& s) {
	
	return LRSplineSurface::Refinement2D {
	  m.kval(s.d, s.ix),
	  m.kval(flip(s.d), s.range.first),
	  m.kval(flip(s.d), s.range.second),
	  s.d,
	  s.d==XFIXED ? mult.first : mult.second};
      });
    return result;
  }

  // ----------------------------------------------------------------------------
  Mesh2D remove_unused_knots(Mesh2D& m)
  // ----------------------------------------------------------------------------
  {
    m.removeUnusedLines(XFIXED);
    m.removeUnusedLines(YFIXED);
    return m;
  }
}; //end anonymous namespace
