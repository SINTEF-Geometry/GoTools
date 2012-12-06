#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"
#include "GoTools/utils/checks.h"

//------------------------------------------------------------------------------

using std::vector;
using std::set;
using std::tuple;
using std::get;
using std::find_if;
using std::find;

namespace Go
{

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
      shared_ptr<Element2D> elem(new Element2D(m.kval(XFIXED, (*mit)[0]),
					       m.kval(YFIXED, (*mit)[1]),
					       m.kval(XFIXED, (*mit)[2]),
					       m.kval(YFIXED, (*mit)[3])));
      // emap[LRSplineSurface::generate_key(elem)] = elem;
      emap[LRSplineSurface::generate_key(m.kval(XFIXED, (*mit)[0]),
					 m.kval(YFIXED, (*mit)[1]))] = elem;
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
LRSplineUtils::insert_basis_function(shared_ptr<LRBSpline2D> b, 
				     const Mesh2D& mesh, 
				     LRSplineSurface::BSplineMap& bmap)
//------------------------------------------------------------------------------
{
  // Add a bspline to the global pool of bsplines, but check first if it
  // exists already. In that case, the bspline scaling factor is updated
  auto key = LRSplineSurface::generate_key(*b, mesh);
  if (bmap.find(key) != bmap.end()) {

    // combine b with the function already present
    shared_ptr<LRBSpline2D> target = bmap[key];
    target->gamma()            += b->gamma();
    target->coefTimesGamma() += b->coefTimesGamma();

    return target.get();
  } 
  // if we got here, there is no pre-existing basis function like b.  
  return (bmap[key] = b).get();
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
			      shared_ptr<LRBSpline2D> bfun,
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
void LRSplineUtils::tensor_split(shared_ptr<LRBSpline2D> bfun, 
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

  for (int iy = 0; iy != (int)y_coefs.size(); ++iy) {
    const double yc = y_coefs[iy];
    for (int ix = 0; ix != (int)x_coefs.size(); ++ix) {
      const double xc = x_coefs[ix];
      shared_ptr<LRBSpline2D> basis(new LRBSpline2D(c_g*yc*xc, 
						    deg_x,
						    deg_y,
						    &x_knots[ix], 
						    &y_knots[iy], 
						    yc * xc * gamma,
						    &tensor_mesh));
      insert_basis_function(basis, tensor_mesh, bmap);
    }
  }
}

//------------------------------------------------------------------------------
void LRSplineUtils::iteratively_split (vector<shared_ptr<LRBSpline2D> >& bfuns, 
				       const Mesh2D& mesh)
//------------------------------------------------------------------------------
{
  // The following set is used to keep track over unique b-spline functions.   
  // b-spline function is here identified by its knotvectors only, as we already
  // assume that degrees are unchanging.  Also, since we expect to find several
  // component of a given b-spline-function, and these must be added up, we will
  // not look at control points or gamma coefficients to determine uniqueness.

  set<shared_ptr<LRBSpline2D>, LRBSpline2DUtils::support_compare> tmp_set;

  shared_ptr<LRBSpline2D> b_split_1, b_split_2;
  bool split_occurred;

  // this closure adds b_spline functions to tmp_set, or combine them if they 
  // are already in it
  auto insert_bfun_to_set = [&tmp_set](shared_ptr<LRBSpline2D> b)->void {
    auto it = tmp_set.find(b);
    if (it == tmp_set.end()) {  // not already in set
      tmp_set.insert(b);
    } else {
    // combine b with the function already present
      (*it)->gamma() += b->gamma();
      (*it)->coefTimesGamma() += b->coefTimesGamma();
    }
  };

  // After a new knot is inserted, there might be bsplines that are no longer
  // minimal. Split those according to knot line information in the mesh
  // keep looping until no more basis functions were inserted
  do {
    tmp_set.clear();
    split_occurred = false;

    for (auto b = bfuns.begin(); b != bfuns.end(); ++b) {
      if (LRBSpline2DUtils::try_split_once(*(*b), mesh, b_split_1, b_split_2)) {
	// this function was splitted.  Throw it away, and keep the two splits
	insert_bfun_to_set(b_split_1);
	insert_bfun_to_set(b_split_2);
	split_occurred = true;
      } else {
	// this function was not split.  Keep it.
	insert_bfun_to_set(*b);
      }
    }

    // moving the collected bsplines over to the vector
    bfuns.clear();
    for (auto b_kv = tmp_set.begin(); b_kv != tmp_set.end(); ++b_kv) 
      bfuns.push_back(*b_kv);

  } while (split_occurred);
}

//------------------------------------------------------------------------------
void LRSplineUtils::iteratively_split2 (vector<LRBSpline2D*>& bsplines,
					const Mesh2D& mesh,
					LRSplineSurface::BSplineMap& bmap)
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

  shared_ptr<LRBSpline2D> b_split_1, b_split_2;
  bool split_occurred;

  // this closure adds b_spline functions to tmp_set, or combine them if they 
  // are already in it
  auto insert_bfun_to_set = [&tmp_set](LRBSpline2D* b)->bool
    {
      auto it = tmp_set.find(b);
      if (it == tmp_set.end())
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
	  // combine b with the function already present
	  (*it)->gamma() += b->gamma();
	  (*it)->coefTimesGamma() += b->coefTimesGamma();
	  // We update the support of b with its replacement.
	  std::vector<Element2D*>::iterator it2 = b->supportedElementBegin();
	  for (it2; it2 < b->supportedElementEnd(); ++it2)
	    {
	      // If there exists a support function already (such as b) it is overwritten.
	      (*it2)->addSupportFunction(*it);
	      (*it)->addSupport(*it2);
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

  vector<shared_ptr<LRBSpline2D> > added_basis;

  do {
    tmp_set.clear(); // Used to store new basis functions for each iteration.
    split_occurred = false;

#ifndef NDEBUG
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
		  added_basis[kr] = added_basis[added_basis.size()-1];
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

    	if (insert_bfun_to_set(b_split_1.get())) // @@sbr deb_iter==0 && ki == 20. ref==4.
	  {
	    // A new LRBspline is created, remember it
	    added_basis.push_back(b_split_1);
	    // Let the elements know about the new bsplines
	    for (size_t kr=0; kr<elements.size(); ++kr)
	      if (b_split_1->overlaps(elements[kr]))
		elements[kr]->addSupportFunction(b_split_1.get());
	      else
		{
		  MESSAGE("No overlap!"); // @@sbr201212 This should not happen.
		  b_split_1->removeSupport(elements[kr]);
		  elements[kr]->removeSupportFunction(b_split_1.get());
		}
	  }

    	if (insert_bfun_to_set(b_split_2.get()))
	  {
	    // A new LRBspline is created, remember it
	    added_basis.push_back(b_split_2);
	    // Let the elements know about the new bsplines
	    for (size_t kr=0; kr<elements.size(); ++kr)
	      if (b_split_2->overlaps(elements[kr]))
		elements[kr]->addSupportFunction(b_split_2.get());
	      else
		{
		  MESSAGE("No overlap!"); // @@sbr201212 This should not happen.
		  b_split_2->removeSupport(elements[kr]);
		  elements[kr]->removeSupportFunction(b_split_2.get());
		}
	  }

    	split_occurred = true;
      } else {
     	// this function was not split.  Keep it.
     	bool was_inserted = insert_bfun_to_set(*b);
	if (!was_inserted)
	  {
//	    MESSAGE("DEBUG: We should remove basis function from added_basis!");
	    // Remove the bspline from the vector of bsplines to add
	    for (size_t kr=0; kr<added_basis.size(); ++kr)
	      if (added_basis[kr].get() == (*b))
		{
		  added_basis[kr] = added_basis[added_basis.size()-1];
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
	{
	  MESSAGE("Already added to map!");
	}
      bmap[key] = added_basis[kr];
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
  const int start_ix = locate_interval(mesh, flip(d), start + fabs(start) * knot_tol, fixed_val, false);
  const int   end_ix = locate_interval(mesh, flip(d), end   - fabs(end)   * knot_tol, fixed_val,  true);
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

  // // We must also update the mesh in the basis functions.
  // auto it = bmap.begin();
  // while (it != bmap.end())
  //   {
  //     it->second.setMesh(&mesh);
  //     ++it;
  //   }

  // If this prev_ix corresponds to our fixed_val, we decrease the value.
  // @@sbr201212 I guess we should do this earlier, but then we need to update the logic above which works ...
  if (mesh.kval(d, prev_ix) == fixed_val)
    prev_ix -= mult;

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

}


}; // end namespace Go

