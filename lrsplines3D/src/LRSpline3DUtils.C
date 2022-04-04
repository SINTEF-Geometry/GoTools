//===========================================================================
//                                                                           
// File: LRSpline3DUtils.C                                                   
//                                                                           
// Created: Wed Mar  6 17:11:23 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/lrsplines3D/LRSpline3DUtils.h"
#include "GoTools/lrsplines3D/LRBSpline3DUtils.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/utils/checks.h"
#include <fstream>
#include <iostream>
using namespace std;

namespace Go
{

//------------------------------------------------------------------------------
  LRSplineVolume::ElementMap 
  LRSpline3DUtils::identify_elements_from_mesh(const Mesh3D& m)
//------------------------------------------------------------------------------
  {
    // The element entities are created, but no bspline information is
    // attached at this stage
    LRSplineVolume::ElementMap emap;
    for (Mesh3DIterator mit = m.begin(); mit != m.end(); ++mit)
      {
        unique_ptr<Element3D> elem(new Element3D(m.kval(XDIR, (*mit)[0]),
						 m.kval(YDIR, (*mit)[1]),
						 m.kval(ZDIR, (*mit)[2]),
						 m.kval(XDIR, (*mit)[3]),
						 m.kval(YDIR, (*mit)[4]),
						 m.kval(ZDIR, (*mit)[5])));
	// emap[LRSplineVolume::generate_key(elem)] = elem;
	LRSplineVolume::ElemKey e_key = LRSplineVolume::generate_key(m.kval(XDIR, (*mit)[0]),
								     m.kval(YDIR, (*mit)[1]),
								     m.kval(ZDIR, (*mit)[2]));
	emap.insert(std::make_pair(e_key, std::move(elem)));
      }

    return emap;
  }


//------------------------------------------------------------------------------
// Add ('remove'=false) or remove ('remove'=true) references to a particular 
// LRBSpline3D for all affected elements in 'emap'.  
  void LRSpline3DUtils::update_elements_with_single_bspline(LRBSpline3D* b, 
							    LRSplineVolume::ElementMap& emap, 
							    const Mesh3D& mesh,
							    bool remove)
//------------------------------------------------------------------------------
{
    const double* kvals_x = mesh.knotsBegin(XDIR);
    const double* kvals_y = mesh.knotsBegin(YDIR);
    const double* kvals_z = mesh.knotsBegin(ZDIR);

    //@@@ VSK. The LRBSpline3D knows wich elements are in its support. Should be
    // possible to make this more effective
    // The element list in the LR B-spline must be up to date
    // LRSplineVolume::ElementMap::const_iterator it_prev = emap.end();
    
  for (int z = b->suppMin(ZDIR); z != b->suppMax(ZDIR); ++z) {
      for (int y = b->suppMin(YDIR); y != b->suppMax(YDIR); ++y) {
	for (int x = b->suppMin(XDIR); x != b->suppMax(XDIR); ++x) {


#if 0 // Removed when extending to 3D.
      // -- following two lines are meant as an optimization, but the effect seems
      // -- to be quite marginal. @@
	  if (mesh.nu(XFIXED, x, y, y+1) < 1) continue; // this cannot be the lower-left -
	  if (mesh.nu(YFIXED, y, x, x+1) < 1) continue; // corner of an element
#endif

	  const LRSplineVolume::ElemKey key = {kvals_x[x], kvals_y[y], kvals_z[z]};

	  LRSplineVolume::ElementMap::const_iterator it;
	  // if (it_prev != emap.end())
	  //   {
	  //     if (key.u_min == it_prev->first.u_min &&
	  // 	  key.v_min == it_prev->first.v_min &&
	  // 	  key.w_min == it_prev->first.w_min)
	  // 	{
	  // 	  it = it_prev;
	  // 	}
	  //     else
	  // 	it = emap.find(key);
	  //   }
	  // else
	    it = emap.find(key);
	  if (it != emap.end()) {
	    // We found an element covered by the support of this basis function.
	    // Update the element
	    if (remove) { // we want to remove it (if it's there)
	      it->second->removeSupportFunction(b);
	    } else {      // we want to insert it (if it's not there already)
	      it->second->addSupportFunction(b);
	      b->addSupport(it->second.get());  // Update bspline with respect to element
	      // it_prev = it;
	      // it_prev++;
	    }
	  }
	}
      }
    }
}

//------------------------------------------------------------------------------
int LRSpline3DUtils::locate_interval(const Mesh3D& m, Direction3D d, double value,
				     double next_value, double prev_value,
				     bool at_end)
//------------------------------------------------------------------------------
{
  double vals[3];
  vals[0] = value;
  vals[1] = next_value;
  vals[2] = prev_value;
  int u_ind = (d == XDIR) ? 0 : ((d == YDIR) ? 2 : 1);
  const double u = vals[u_ind];
  const double v = vals[(u_ind+1)%3];
  const double w = vals[(u_ind+2)%3];

  int u_ix, v_ix, w_ix;
  const bool found = (at_end) ?
    Mesh3DUtils::identify_patch_upper_right(m, u, v, w, u_ix, v_ix, w_ix) : 
    Mesh3DUtils::identify_patch_lower_left (m, u, v, w, u_ix, v_ix, w_ix);

  if (!found) 
    THROW("locate_interval() : parameter outside domain.");
  
  return (d == XDIR) ? u_ix : ((d == YDIR) ? v_ix : w_ix);
}

//------------------------------------------------------------------------------
vector<int> LRSpline3DUtils::set_uniform_meshlines(Direction3D d, Mesh3D& mesh)
//------------------------------------------------------------------------------
{
  vector<int> mults(mesh.numDistinctKnots(d));
  const int line_len1 = mesh.numDistinctKnots(next(d)) - 1;
  const int line_len2 = mesh.numDistinctKnots(prev(d)) - 1;
  for (int i = 0; i != mesh.numDistinctKnots(d); ++i) {
    const int mul = mesh.largestMultInLine(d, i);
    mesh.setMult(d, i, 0, line_len1, 0, line_len2, mul);
    mults[i] = mul;
  }
  return mults;
}

//------------------------------------------------------------------------------
bool LRSpline3DUtils::all_meshlines_uniform(Direction3D d, const Mesh3D& m)
//------------------------------------------------------------------------------
{
  const int orto_intervals1 = m.numDistinctKnots(next(d)) - 1;
  const int orto_intervals2 = m.numDistinctKnots(prev(d)) - 1;
  for (int i = 0; i != m.numDistinctKnots(d); ++i)
    if (m.nu(d, i, 0, orto_intervals1, 0, orto_intervals2) != m.largestMultInLine(d, i)) return false;
  return true;
}



//------------------------------------------------------------------------------
// returns a pointer to the new (or existing) function
LRBSpline3D* 
LRSpline3DUtils::insert_basis_function(unique_ptr<LRBSpline3D>& b,
				       const Mesh3D& mesh, 
				       LRSplineVolume::BSplineMap& bmap)
//------------------------------------------------------------------------------
{
  // Add a bspline to the global pool of bsplines, but check first if it
  // exists already. In that case, the bspline scaling factor is updated
  auto key = LRSplineVolume::generate_key(*b, mesh);
  if (bmap.find(key) != bmap.end()) {

    // combine b with the function already present
    LRBSpline3D* target = bmap[key].get();
    target->gamma()            += b->gamma();
    target->coefTimesGamma() += b->coefTimesGamma();

    return target;
  } 
  // if we got here, there is no pre-existing basis function like b.  
  bmap.insert(std::make_pair(key, std::move(b)));
  return b.get();
  //return (bmap[key] = b).get();
}

  //------------------------------------------------------------------------------
  SplineVolume* LRSpline3DUtils::fullTensorProductVolume(const LRSplineVolume& lr_spline_vol)
  //------------------------------------------------------------------------------
{
    bool is_full_tensor = (lr_spline_vol.isFullTensorProduct());
    shared_ptr<LRSplineVolume> lr_spline_vol_copy;
    const LRSplineVolume* full_tp_vol = &lr_spline_vol;
    if (!is_full_tensor)
      {
        lr_spline_vol_copy = shared_ptr<LRSplineVolume>(lr_spline_vol.clone());
        lr_spline_vol_copy->expandToFullTensorProduct();
        full_tp_vol = lr_spline_vol_copy.get();
      }

    int num_coefs_u, num_coefs_v, num_coefs_w;
    int order_u = full_tp_vol->degree(XDIR) + 1;
    int order_v = full_tp_vol->degree(YDIR) + 1;
    int order_w = full_tp_vol->degree(ZDIR) + 1;
    int dim = full_tp_vol->dimension();
    bool rational = full_tp_vol->rational();

    // The basis functions should be ordered with increasing u-parameter first,
    // then the v-parameter, then the w-parameter.
    auto iter = full_tp_vol->basisFunctionsBegin();
    //int kdim = (rational) ? dim + 1 : dim;
    vector<double> vol_coefs;
    while (iter != full_tp_vol->basisFunctionsEnd())
      {
        Point coef = iter->second->coefTimesGamma();//Coef();
        vol_coefs.insert(vol_coefs.end(), coef.begin(), coef.end());
        if (rational)
          vol_coefs.insert(vol_coefs.end(), iter->second->weight());
          ++iter;
      }

    const Mesh3D& mesh = full_tp_vol->mesh();
    vector<double> knots_u = mesh.getKnots(XDIR, 0, 0);
    vector<double> knots_v = mesh.getKnots(YDIR, 0, 0);
    vector<double> knots_w = mesh.getKnots(ZDIR, 0, 0);
    num_coefs_u = (int)knots_u.size() - order_u;
    num_coefs_v = (int)knots_v.size() - order_v;
    num_coefs_w = (int)knots_w.size() - order_w;

    SplineVolume* spline_vol = new SplineVolume(num_coefs_u, num_coefs_v, num_coefs_w,
                                                order_u, order_v, order_w,
                                                knots_u.begin(), knots_v.begin(), knots_w.begin(),
                                                vol_coefs.begin(),
                                                dim,
                                                rational);

    return spline_vol;
}



//------------------------------------------------------------------------------
void 
LRSpline3DUtils::iteratively_split (vector<std::unique_ptr<LRBSpline3D> >& bfuns,
				    const Mesh3D& mesh,
				    vector<unique_ptr<BSplineUniLR> >& bspline_vec1,
				    vector<unique_ptr<BSplineUniLR> >& bspline_vec2,
				    vector<unique_ptr<BSplineUniLR> >& bspline_vec3)
//------------------------------------------------------------------------------
{
  // The following set is used to keep track over unique b-spline functions.   
  // b-spline function is here identified by its knotvectors only, as we already
  // assume that degrees are unchanging.  Also, since we expect to find several
  // component of a given b-spline-function, and these must be added up, we will
  // not look at control points or gamma coefficients to determine uniqueness.

  set<LRBSpline3D*, support_compare> tmp_set;
  bool split_occurred;

  // this closure adds b_spline functions to tmp_set, or combine them if they 
  // are already in it
  auto insert_bfun_to_set = [&tmp_set](LRBSpline3D* b)->bool {
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
  int innermult1 = mesh.largestInnerMult(XDIR);
  int innermult2 = mesh.largestInnerMult(YDIR);
  int innermult3 = mesh.largestInnerMult(ZDIR);
  do {
    tmp_set.clear();
    split_occurred = false;
    int deb_ptr = 0; // @@sbr201304 Remove when done debugging!
    for (auto b = bfuns.begin(); b != bfuns.end(); ++b, ++deb_ptr) {
      LRBSpline3D *b_split_1 = NULL;
      LRBSpline3D *b_split_2 = NULL;
      if (LRBSpline3DUtils::try_split_once(*(*b), mesh, innermult1, innermult2,
					   innermult3, bspline_vec1, 
					   bspline_vec2, bspline_vec3, 
					   b_split_1, b_split_2)) {
          // this function was split.  Throw it away, and keep the two splits
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
    for (auto b_kv = tmp_set.begin(); b_kv != tmp_set.end(); ++b_kv)
      bfuns.insert(bfuns.end(), unique_ptr<LRBSpline3D>(*b_kv));

  } while (split_occurred);
}


//------------------------------------------------------------------------------
void 
LRSpline3DUtils::iteratively_split2 (vector<LRBSpline3D*>& bsplines,
				     const Mesh3D& mesh,
				     LRSplineVolume::BSplineMap& bmap,
				     double* domain,
				     vector<unique_ptr<BSplineUniLR> >& bspline_vec1,
				     vector<unique_ptr<BSplineUniLR> >& bspline_vec2,
				     vector<unique_ptr<BSplineUniLR> >& bspline_vec3,
				     bool support)
//------------------------------------------------------------------------------
{
  // The following set is used to keep track over unique b-spline functions.
  // b-spline function is here identified by its knotvectors only, as we already
  // assume that degrees are unchanging.  Also, since we expect to find several
  // components of a given b-spline-function, and these must be added up, we will
  // not look at control points or gamma coefficients to determine uniqueness.
  set<LRBSpline3D*, support_compare> tmp_set;

  bool split_occurred;

  // this closure adds b_spline functions to tmp_set, or combine them if they
  // are already in it
  auto insert_bfun_to_set = [&tmp_set](LRBSpline3D* b,
                                       LRSplineVolume::BSplineMap& bmap,
                                       double* domain,
				       bool support)->bool
  {
    auto it = tmp_set.find(b);
    bool overlap = (domain == NULL) ? false : b->overlaps(domain);
    LRBSpline3D* other = NULL;
    if (it != tmp_set.end())
      {
        other = (*it);
      }
    else if (!overlap)
      {
        // Search for an identical B-spline in the entire domain
        LRSplineVolume::BSKey key = LRSplineVolume::generate_key(*b);
        auto it3 = bmap.find(key);
        if (it3 != bmap.end() && it3->second.get() != b)
          other = it3->second.get();
      }
    if (it == tmp_set.end() && !other)
      {
        // We must check if the last element of tmp_set is equal.
        auto it2 = tmp_set.end();
        int set_size_pre = (int)tmp_set.size();
        if (set_size_pre > 0)
          it2--;
        bool last_elem_equal = (set_size_pre > 0 && (support_equal(*it2, b)));
        if (last_elem_equal)
	  std::cout << "DEBUG: Last element is equal to new element!" << std::endl;
	//MESSAGE("DEBUG: Last element is equal to new element!");
        // not already in set
        tmp_set.insert(b);
        int set_size_post = (int)tmp_set.size();
        if (set_size_pre == set_size_post)
          //MESSAGE("DEBUG: It seems we tried to insert an element already present!");
	std::cout << "DEBUG: It seems we tried to insert an element already present!" << std::endl;
        return true;
      }
    else
      {
        bool rat = b->rational();
        if (rat)
          { // We must alter the weight of the second basis function to match that of our reference.
            double b_w = b->weight();
            double it_w = other->weight();
            double weight = b_w + it_w;
            // We must rescale the coefs to reflect the change in weight.
            b->coefTimesGamma() *= b_w/weight;
            other->coefTimesGamma() *= it_w/weight;
            b->weight() = (*it)->weight() = weight;
          }
        // combine b with the function already present
        other->gamma() += b->gamma();
        other->coefTimesGamma() += b->coefTimesGamma();

	if (support)
	  {
	    // We update the support of b with its replacement.
	    std::vector<Element3D*>::iterator it2 = b->supportedElementBegin();
	    for (; it2 < b->supportedElementEnd(); ++it2)
	      {
		// Note that in subsequent divisions, the new bspline may point to
		// elements which are not in the support of the already existing one
		// Thus, check for overlap
		if (other->overlaps((*it2)))
		  {
		    // If there exists a support function already (such as b) it is overwritten.
		    (*it2)->addSupportFunction(other);
		    other->addSupport(*it2);
		  }
		else
		  {
		    //std::cout << "No overlap, element " << *it2 << ", bspline " << b << std::endl;
		    int stop_break = 1;
		  }
	      }

	    // Finally we remove all elements from b.
	    for (auto it4=b->supportedElementBegin(); 
		 it4!=b->supportedElementEnd(); ++it4)
	      (*it4)->removeSupportFunction(b);
	    b->removeSupportedElements();
	    // while (b->nmbSupportedElements() > 0)
	    //   {
	    //     auto it2 = b->supportedElementBegin();
	    //     (*it2)->removeSupportFunction(b);
	    //     b->removeSupport(*it2);
	    //   }
	  }

        return false;
      }
  };

  // After a new knot is inserted, there might be bsplines that are no longer
  // minimal. Split those according to knot line information in the mesh
  // keep looping until no more basis functions were inserted

  vector<unique_ptr<LRBSpline3D> > added_basis;

  int innermult1 = mesh.largestInnerMult(XDIR);
  int innermult2 = mesh.largestInnerMult(YDIR);
  int innermult3 = mesh.largestInnerMult(ZDIR);
#ifdef DEBUG10
  std::ofstream of("bsplit.tmp");
#endif
  do { // Loop is run until no more splits occur.
    tmp_set.clear(); // Used to store new basis functions for each iteration.
    split_occurred = false;

#ifdef DEBUG10
    of << "num bsplines: " << bsplines.size() << std::endl;
#endif
    int ki = 0;
    for (auto b = bsplines.begin(); b != bsplines.end(); ++b, ++ki) {

      LRBSpline3D *b_split_1 = NULL;
      LRBSpline3D *b_split_2 = NULL;

#ifdef DEBUG10
      of << *b << " " << (*b)->umin() << " " << (*b)->vmin();
      of << " " << (*b)->wmin() << " " << (*b)->umax();
      of << " " << (*b)->vmax() << " " << (*b)->wmax() << std::endl;
      of << "Uni count: " << (*b)->getUnivariate(XDIR)->getCount();
      of << ", " << (*b)->getUnivariate(YDIR)->getCount();
      of << ", " << (*b)->getUnivariate(ZDIR)->getCount() << std::endl;
      of << "Uni: " << (*b)->getUnivariate(XDIR) << " ";
      (*b)->getUnivariate(XDIR)->write(of);
      of << ", " << (*b)->getUnivariate(YDIR) << " ";
      (*b)->getUnivariate(YDIR)->write(of);
      of << ", " << (*b)->getUnivariate(ZDIR) << " ";
      (*b)->getUnivariate(ZDIR)->write(of);
      LRSplineVolume::BSKey tmp_key = LRSplineVolume::generate_key(*(*b));
#endif
      // Fetch all elements
      vector<Element3D*> elements = (*b)->supportedElements();

      if (LRBSpline3DUtils::try_split_once(*(*b), mesh, innermult1, innermult2,
					   innermult3, bspline_vec1, bspline_vec2,
					   bspline_vec3, b_split_1, b_split_2)) {
#ifdef DEBUG10
	of << "Split " << b_split_1 << " " << b_split_1->umin() << " " << b_split_1->vmin();
	of << " " << b_split_1->wmin() << " " << b_split_1->umax();
	of << " " << b_split_1->vmax() << " " << b_split_1->wmax() << std::endl;
      LRSplineVolume::BSKey tmp_key1 = LRSplineVolume::generate_key(*b_split_1);
	of << "Uni count: " << b_split_1->getUnivariate(XDIR)->getCount();
	of << ", " << b_split_1->getUnivariate(YDIR)->getCount();
	of << ", " << b_split_1->getUnivariate(ZDIR)->getCount() << std::endl;
	of << "Uni: " << b_split_1->getUnivariate(XDIR) << " ";
	b_split_1->getUnivariate(XDIR)->write(of);
	of << ", " << b_split_1->getUnivariate(YDIR) << " ";
	b_split_1->getUnivariate(YDIR)->write(of);
	of << ", " << b_split_1->getUnivariate(ZDIR) << " ";
	b_split_1->getUnivariate(ZDIR)->write(of);	
	of << "Split " << b_split_2 << " " << b_split_2->umin() << " " << b_split_2->vmin();
	of << " " << b_split_2->wmin() << " " << b_split_2->umax();
	of << " " << b_split_2->vmax() << " " << b_split_2->wmax() << std::endl;
	LRSplineVolume::BSKey tmp_key2 = LRSplineVolume::generate_key(*b_split_2);
	of << "Uni count: " << b_split_2->getUnivariate(XDIR)->getCount();
	of << ", " << b_split_2->getUnivariate(YDIR)->getCount();
	of << ", " << b_split_2->getUnivariate(ZDIR)->getCount() << std::endl;
	of << "Uni: " << b_split_2->getUnivariate(XDIR) << " ";
	b_split_2->getUnivariate(XDIR)->write(of);
	of << ", " << b_split_2->getUnivariate(YDIR) << " ";
	b_split_2->getUnivariate(YDIR)->write(of);
	of << ", " << b_split_2->getUnivariate(ZDIR) << " ";
	b_split_2->getUnivariate(ZDIR)->write(of);
#endif
        // this function was split.  Throw it away, and keep the two splits
        // @@@ VSK. Must also update bmap and set element pointers
        // Remove bspline from element
        for (size_t kr=0; kr<elements.size(); ++kr)
          {
            elements[kr]->removeSupportFunction(*b);
          }

	// Remove bspline from bspline map
	LRSplineVolume::BSKey key = LRSplineVolume::generate_key(*(*b));
	auto it = bmap.find(key);
	if (it != bmap.end())
	  {
	    bmap.erase(it);
	  }
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

	// Until the elements are split, let the new bsplines store all
	// elements from their origin in their support

	if (support)
	  {
	    // Since the elements have not yet been split, the support is the same.
	    b_split_1->setSupport(elements);
	    b_split_2->setSupport(elements);
	  }

	if (insert_bfun_to_set(b_split_1, bmap, domain, support)) // @@sbr deb_iter==0 && ki == 20. ref==4.
	  {
	    //std::cout << "deb_iter: " << deb_iter << ", ki" << ki << ", b_split_1: " << b_split_1 << std::endl;
	    // A new LRBspline is created, remember it
	    added_basis.push_back(unique_ptr<LRBSpline3D>(b_split_1));

	    if (support)
	      {
		// Let the elements know about the new bsplines
		for (size_t kr=0; kr<elements.size(); ++kr)
		  if (b_split_1->overlaps(elements[kr]))
		    elements[kr]->addSupportFunction(b_split_1);
		  else
		    {
		      b_split_1->removeSupport(elements[kr]);
		      elements[kr]->removeSupportFunction(b_split_1);
		    }
	      }
#ifdef DEBUG10
	    of << "First new B-spline inserted" << std::endl;
#endif
	  }
	else
	  { // Memory management.
	    delete b_split_1;
#ifdef DEBUG10
	    of << "First new B-spline merged with previous" << std::endl;
#endif
	  }

	if (insert_bfun_to_set(b_split_2, bmap, domain, support))
	  {
	    //std::cout << "deb_iter: " << deb_iter << ", ki" << ki << ", b_split_2: " << b_split_2 << std::endl;
	    // A new LRBspline is created, remember it
	    added_basis.push_back(unique_ptr<LRBSpline3D>(b_split_2));
	    if (support)
	      {
		// Let the elements know about the new bsplines
		for (size_t kr=0; kr<elements.size(); ++kr)
		  if (b_split_2->overlaps(elements[kr]))
		    elements[kr]->addSupportFunction(b_split_2);
		  else
		    {
		      b_split_2->removeSupport(elements[kr]);
		      elements[kr]->removeSupportFunction(b_split_2);
		    }
	      }
#ifdef DEBUG10
	    of << "Second new B-spline inserted" << std::endl;
#endif
	  }
	else
	  { // Memory management.
	    delete b_split_2;
#ifdef DEBUG10
	    of << "Second new B-spline merged with previous" << std::endl;
#endif
	  }

        split_occurred = true;
      } else {
#ifdef DEBUG10
	of << "No split" << std::endl;
#endif
        // this function was not split.  Keep it.
        bool was_inserted = insert_bfun_to_set(*b, bmap, domain, support);
        if (!was_inserted)
          {
	    if (support)
	      {
		// Remove bspline from element
		for (size_t kr=0; kr<elements.size(); ++kr)
		  {
		    elements[kr]->removeSupportFunction(*b);
		  }
	      }
            //	    MESSAGE("DEBUG: We should remove basis function from added_basis!");
            // Remove the B-spline also from the bmap if present
            // Remove the bspline from the vector of bsplines to add
            size_t kr;
            bool found = false;
            for (kr=0; kr<added_basis.size(); ++kr)
              if (added_basis[kr].get() == (*b))
                {
                  std::swap(added_basis[kr], added_basis[added_basis.size()-1]);
                  added_basis.pop_back();
                  found = true;
                  break;
                }

	    if (!found)
	      {
		LRSplineVolume::BSKey key = LRSplineVolume::generate_key(*(*b));
		auto it = bmap.find(key);
		if (it != bmap.end())
		  {
		    // Remove
		    bmap.erase(it);
		  }
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

  } while (split_occurred);
#ifdef DEBUG10
  of << "Finished splitting" << std::endl;
#endif 
#if 1//ndef NDEBUG
  {
    vector<LRBSpline3D*> bas_funcs;
    for (auto iter = bmap.begin(); iter != bmap.end(); ++iter)
      {
        bas_funcs.push_back((*iter).second.get());
      }
    //puts("Remove when done debugging!");
    int stop_break = 1;
  }
#endif

#ifdef DEBUG10
  of << "Add to bmap" << std::endl;
#endif  
  // Add new basis functions to bmap
  for (size_t kr=0; kr<added_basis.size(); ++kr)
    {
      //LRBSpline3D* tmp_b = added_basis[kr].get();
      LRSplineVolume::BSKey key = LRSplineVolume::generate_key(*added_basis[kr]);
      auto it = bmap.find(key);
      if (it != bmap.end())
        { // @@ I guess we handle this by adding
          //MESSAGE("Already added to map! This is a bug. Expect core dump if not fixed ...");
          // @@sbr201305 This will in a lost pointer and most likely a core dump!
          LRBSpline3D* b = added_basis[kr].get();
          //#if 1
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

	  if (support)
	    {
	      // We update the support of b with its replacement.
	      std::vector<Element3D*>::iterator it2 = b->supportedElementBegin();
	      for (/*it2*/; it2 < b->supportedElementEnd(); ++it2)
		{
		  // If there exists a support function already (such as b) it is overwritten.
		  (*it2)->addSupportFunction(it->second.get());
		  (it->second)->addSupport(*it2);

		  // Remove b-spline from element
		  (*it2)->removeSupportFunction(b);
		}
	    }
         //#endif

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
#ifdef DEBUG10
  of << "Finished add to bmap" << std::endl;
#endif
 }


// Helper function used by 'tensor_split'.  Identifies the knots that need to
// be inserted into some reference vector, in order to bring multiplicty of each
// knot up to the multiplicity expressed in 'mult' (including those with a
// implicit current multiplicity of zero.  This function is used by 'tensor_split'.
// This is really a univariate function so does not need reimplementing for each dim.
//------------------------------------------------------------------------------
  vector<int> LRSpline3DUtils::knots_to_insert(const vector<int>& ref,
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
  // This is really a univariate function so does not need reimplementing for each dim.
//------------------------------------------------------------------------------
  double LRSpline3DUtils::compute_alpha(int degree,
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
// This is really a univariate function (except direction)
// so does not need reimplementing for each dim.
//------------------------------------------------------------------------------
  tuple<vector<double>, vector<int> >
  LRSpline3DUtils::insert_knots(const vector<int>& new_knots,
                                unique_ptr<LRBSpline3D>& bfun,
                                const Direction3D d,
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
  void 
  LRSpline3DUtils::tensor_split(unique_ptr<LRBSpline3D>& bfun,
				const vector<int>& x_mults,
				const vector<int>& y_mults,
				const vector<int>& z_mults,
				const Mesh3D& tensor_mesh,
				vector<unique_ptr<BSplineUniLR> >& bspline_vec1,
				vector<unique_ptr<BSplineUniLR> >& bspline_vec2,
				vector<unique_ptr<BSplineUniLR> >& bspline_vec3,
				LRSplineVolume::BSplineMap& bmap)
//------------------------------------------------------------------------------
{
    const vector<int> kx = knots_to_insert(bfun->kvec(XDIR), x_mults);
    const vector<int> ky = knots_to_insert(bfun->kvec(YDIR), y_mults);
    const vector<int> kz = knots_to_insert(bfun->kvec(ZDIR), z_mults);

    const double* const x_kvals = tensor_mesh.knotsBegin(XDIR);
    const double* const y_kvals = tensor_mesh.knotsBegin(YDIR);
    const double* const z_kvals = tensor_mesh.knotsBegin(ZDIR);

    const auto x_coefs_kvec = insert_knots(kx, bfun, XDIR, x_kvals);
    const auto y_coefs_kvec = insert_knots(ky, bfun, YDIR, y_kvals);
    const auto z_coefs_kvec = insert_knots(kz, bfun, ZDIR, z_kvals);

    const vector<double>& x_coefs = get<0>(x_coefs_kvec);
    const vector<double>& y_coefs = get<0>(y_coefs_kvec);
    const vector<double>& z_coefs = get<0>(z_coefs_kvec);
    const vector<int>& x_knots = get<1>(x_coefs_kvec);
    const vector<int>& y_knots = get<1>(y_coefs_kvec);
    const vector<int>& z_knots = get<1>(z_coefs_kvec);

    const int deg_x = bfun->degree(XDIR);
    const int deg_y = bfun->degree(YDIR);
    const int deg_z = bfun->degree(ZDIR);
    const double gamma = bfun->gamma();
    const Point& c_g = bfun->coefTimesGamma();
    const double weight = bfun->weight();
    const bool rational = bfun->rational();

     // Univariate B-splines
    int left1 = 0, left2 = 0, left3 = 0;
    for (int ix = 0; ix != (int)x_coefs.size(); ++ix) {
      BSplineUniLR *tmpu = new BSplineUniLR(1, deg_x, x_knots.begin()+ix, 
					    &tensor_mesh);

      // Check if the B-spline exists already
      bool found1 = BSplineUniUtils::identify_bsplineuni(tmpu, bspline_vec1, left1);
      if (found1)
	delete tmpu;
      else
	BSplineUniUtils::insert_univariate(bspline_vec1, tmpu, left1);
    }

    for (int iy = 0; iy != (int)y_coefs.size(); ++iy) {
      BSplineUniLR *tmpv = new BSplineUniLR(2, deg_y, y_knots.begin()+iy, 
					    &tensor_mesh);

      // Check if the B-spline exists already
      bool found2 = BSplineUniUtils::identify_bsplineuni(tmpv, bspline_vec2, left2);
      if (found2)
	delete tmpv;
      else
	BSplineUniUtils::insert_univariate(bspline_vec2, tmpv, left2);
    }

   for (int iz = 0; iz != (int)z_coefs.size(); ++iz) {
      const double zc = z_coefs[iz];
      BSplineUniLR *tmpw = new BSplineUniLR(3, deg_z, z_knots.begin()+iz, 
					    &tensor_mesh);

      // Check if the B-spline exists already
      bool found3 = BSplineUniUtils::identify_bsplineuni(tmpw, bspline_vec3, left3);
      if (found3)
	delete tmpw;
      else
	BSplineUniUtils::insert_univariate(bspline_vec3, tmpw, left3);

      for (int iy = 0; iy != (int)y_coefs.size(); ++iy) {
        const double yc = y_coefs[iy];
        for (int ix = 0; ix != (int)x_coefs.size(); ++ix) {
          const double xc = x_coefs[ix];
          unique_ptr<LRBSpline3D> basis(new LRBSpline3D(c_g*zc*yc*xc,
                                                        weight,
							bspline_vec1[left1].get(),
							bspline_vec2[left2].get(),
							bspline_vec3[left3].get(),
                                                        zc * yc * xc * gamma,
                                                        rational));
        insert_basis_function(basis, tensor_mesh, bmap);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Inserts a refinement into the mesh and increments indices in the BSplineMap accordingly.
// The function returns six integers:
//   - The first is the last nonlarger knot value index in the fixed direction.
//   - The second is the index of the meshline on which the inserted meshrectangle is located
//     (index of the fixed parameter value in the global knot vector.  
//   - The third is the start index of the mesh rectangle in the next parameter direction.
//   - The fourth is the end index of the mesh rectangle in the next parameter direction.
//   - The fifth is the start index of the mesh rectangle in the next+1 parameter direction.
//   - The sixth is the end index of the mesh rectangle in the next+1 parameter direction.
//
// NB: Does NOT carry out any splitting of the basis functions in the BSplineMap.  It is the
//     caller's responsibility to do this afterwards!  (As such, this function is not 
//     intended for use by itself, but should only be called from one of the LRSpline::refine()
//     methods.
//------------------------------------------------------------------------------
  tuple<int, int, int, int, int, int>
  LRSpline3DUtils::refine_mesh(Go::Direction3D d, double fixed_val,
			       double start1, double end1,
			       double start2, double end2,
			       int mult, bool absolute,
			       int spline_degree, double knot_tol,
			       Go::Mesh3D& mesh,
			       std::vector<std::unique_ptr<BSplineUniLR> >& bsplines,
			       bool& refined)
//------------------------------------------------------------------------------
{
  if (mult > spline_degree + 1) 
    THROW("Cannot refine with multiplicity higher than degree+1.");

  // Determining the anchor points of the new meshline to be inserted (it must end in an existing
  // orthogonal meshline of multiplicity >= 1).
  // @@ EXPLAIN THE WEIRD USE OF SIGNS AND TOLERANCE BELOW
  const int start_ix1 = locate_interval(mesh, next(d),
					start1 + fabs(start1) * knot_tol,
					start2,
					fixed_val,
					false);
  const int   end_ix1 = locate_interval(mesh, next(d),
					end1   - fabs(end1)   * knot_tol,
					end2,
					fixed_val,
					true);
  const int start_ix2 = locate_interval(mesh, prev(d),
					start2 + fabs(start2) * knot_tol,
					fixed_val,
					start1,
					false);
  const int   end_ix2 = locate_interval(mesh, prev(d),
					end2   - fabs(end2)   * knot_tol,
					fixed_val,
					end1,
					true);
  //const int start_ix = locate_interval(mesh, flip(d), start * (1 + knot_tol), fixed_val, false);
  //const int   end_ix = locate_interval(mesh, flip(d), end   * (1 - knot_tol), fixed_val,  true);

  // Fetch the last nonlarger knot value index in the fixed direction
  int prev_ix = Mesh3DUtils::last_nonlarger_knotvalue_ix(mesh, d, fixed_val);

  int fixed_ix; // to be set below. Knot value index of the new knot
  // const auto existing_it = find_if(mesh.knotsBegin(d),
  // 				   mesh.knotsEnd(d), 
  // 				   [=](double x)->bool
  // 				   {return (abs(fixed_val - x) < abs(fixed_val) * knot_tol);});

  // if (existing_it != mesh.knotsEnd(d)) { // increase multiplicity of existing meshrectangles
  refined = true;
  if (fabs(mesh.kval(d, prev_ix) - fixed_val) > knot_tol &&
      fabs(mesh.kval(d, prev_ix+1) - fixed_val) < knot_tol)
    ++prev_ix;
  if (fabs(mesh.kval(d, prev_ix) - fixed_val) < knot_tol)
    {
      // increase multiplicity of existing meshrectangles
      //fixed_ix = int(existing_it - mesh.knotsBegin(d));
      fixed_ix = prev_ix;

      // check that the proposed multiplicity modification is legal
      bool perform = true;
      for (int j = start_ix2; j < end_ix2; ++j)
      	for (int i = start_ix1; i < end_ix1; ++i)
      	  {
      	    const int cur_m = mesh.nu(d, fixed_ix, i, i+1, j, j+1);
      	    if (absolute && (cur_m > mult))
	      {
		MESSAGE("Cannot decrease multiplicity.");
		perform = false;
	      }
      	    else if (!absolute && (cur_m+mult > spline_degree + 1))
	      {
		MESSAGE("Cannot increase multiplicity.");
		perform = false;
	      }
      	  }
      if (perform)
	{
	  // set or increment multiplicity
	  if (absolute)
	    refined = mesh.setMult(d, fixed_ix, start_ix1, end_ix1,
				   start_ix2, end_ix2, mult);
	  else
	    mesh.incrementMult(d, fixed_ix, start_ix1, end_ix1,
					 start_ix2, end_ix2, mult);
	}
    }
  else
    { 
      // Insert a new rectangle, and set relevant part to desired multiplicity.
      // I.e. we set 0 as the initial multiplicity.
      fixed_ix = mesh.insertRectangle(d, fixed_val, 0);
      mesh.setMult(d, fixed_ix, start_ix1, end_ix1, start_ix2, end_ix2, mult);

      // change index of _all_ basis functions who refer to knot values with indices >= inserted one
      LRSplineUtils::increment_knotvec_indices(bsplines, fixed_ix);
    }

  // // We must also update the mesh in the basis functions.
  // auto it = bmap.begin();
  // while (it != bmap.end())
  //   {
  //     it->second->setMesh(&mesh);
  //     ++it;
  //   }

  // If this prev_ix corresponds to our fixed_val, we decrease the value.
  // @@sbr201212 I guess we could do this earlier, but then we need to update the working code above ...
  if (mesh.kval(d, prev_ix) == fixed_val)
    prev_ix -= mult;

  return tuple<int, int, int, int, int, int>(prev_ix, fixed_ix, start_ix1, end_ix1, start_ix2, end_ix2);
}
#if 0
//==============================================================================
void  LRSpline3DUtils::zero_knot(Go::Direction3D d, double fixed_val,
				 double knot_tol,Go::Mesh3D& mesh,
				 std::vector<std::unique_ptr<BSplineUniLR> >& bsplines)
//==============================================================================
{
  // Fetch the last nonlarger knot value index in the fixed direction
  int prev_ix = Mesh3DUtils::last_nonlarger_knotvalue_ix(mesh, d, fixed_val);
  if (fabs(mesh.kval(d, prev_ix) - fixed_val) >= knot_tol)
    {
      // Insert a new rectangle with zero multiplicity
      int fixed_ix = mesh.insertRectangle(d, fixed_val, 0);
       // change index of _all_ basis functions who refer to knot values with indices >= inserted one
      LRSplineUtils::increment_knotvec_indices(bsplines, fixed_ix);
    }
}
#endif
//==============================================================================
void LRSpline3DUtils::get_affected_bsplines(const std::vector<LRSplineVolume::Refinement3D>& refs, 
					    const LRSplineVolume::ElementMap& emap,
					    double knot_tol, const Mesh3D& mesh,
					    std::vector<LRBSpline3D*>& affected)
//==============================================================================
{
  std::set<LRBSpline3D*> all_bsplines;
  for (size_t ki=0; ki<refs.size(); ++ki)
    {
      const LRSplineVolume::Refinement3D& r = refs[ki];
     const int start_ix1 = locate_interval(mesh, next(r.d),
					    r.start1 + fabs(r.start1) * knot_tol,
					    r.start2,
					    r.kval,
					    false);
      const int   end_ix1 = locate_interval(mesh, next(r.d),
					    r.end1   - fabs(r.end1)   * knot_tol,
					    r.end2,
					    r.kval,
					    true);
      const int start_ix2 = locate_interval(mesh, prev(r.d),
					    r.start2 + fabs(r.start2) * knot_tol,
					    r.kval,
					    r.start1,
					    false);
      const int   end_ix2 = locate_interval(mesh, prev(r.d),
					    r.end2   - fabs(r.end2)   * knot_tol,
					    r.kval,
					    r.end1,
					    true);
      int prev_ix = Mesh3DUtils::last_nonlarger_knotvalue_ix(mesh, r.d, r.kval);

      for (int i = start_ix1; i != end_ix1; ++i)
	{
	  for (int j = start_ix2; j != end_ix2; ++j)
	    {
	      int prev_ix_curr =
		Mesh3DUtils::search_downwards_for_nonzero_multiplicity(mesh, r.d,
								       prev_ix, i, j);
	      // Check if the specified element exists in 'emap'
	      int u_ix = (r.d == XDIR) ? prev_ix_curr : ((r.d == YDIR) ? j : i);
	      int v_ix = (r.d == YDIR) ? prev_ix_curr : ((r.d == ZDIR) ? j : i);
	      int w_ix = (r.d == ZDIR) ? prev_ix_curr : ((r.d == XDIR) ? j : i);
	      LRSplineVolume::ElementMap::key_type key = {mesh.kval(XDIR, u_ix),
							  mesh.kval(YDIR, v_ix),
							  mesh.kval(ZDIR, w_ix)};
	      auto it = emap.find(key);
	      if (it != emap.end())
		{
		  // The element exists. Collect bsplines
		  Element3D* curr_el = (*it).second.get();
		  all_bsplines.insert(curr_el->supportBegin(), curr_el->supportEnd());
		}
	    }
	}
    }
  affected.insert(affected.end(), all_bsplines.begin(), all_bsplines.end());
}

//==============================================================================
namespace {
int compare_u_par(const void* el1, const void* el2)
{
  if (((double*)el1)[0] < ((double*)el2)[0])
    return -1;
  else if (((double*)el1)[0] > ((double*)el2)[0])
    return 1;
  else
    return 0;
}

int compare_v_par(const void* el1, const void* el2)
{
  if (((double*)el1)[1] < ((double*)el2)[1])
    return -1;
  else if (((double*)el1)[1] > ((double*)el2)[1])
    return 1;
  else
    return 0;
}

int compare_w_par(const void* el1, const void* el2)
{
  if (((double*)el1)[2] < ((double*)el2)[2])
    return -1;
  else if (((double*)el1)[2] > ((double*)el2)[2])
    return 1;
  else
    return 0;
}
}

bool LRSpline3DUtils::support_equal(const LRBSpline3D* b1, const LRBSpline3D* b2)
{
  // to compare b1 and b2, compare the x-knotvectors.  If these are identical, compare
  // the y-knotvectors instead.
  const int tmp1 = compare_seq(b1->kvec(XDIR).begin(), b1->kvec(XDIR).end(),
                               b2->kvec(XDIR).begin(), b2->kvec(XDIR).end());
  if (tmp1 != 0) return false;
  const int tmp2 = compare_seq(b1->kvec(YDIR).begin(), b1->kvec(YDIR).end(),
                               b2->kvec(YDIR).begin(), b2->kvec(YDIR).end());
  if (tmp2 != 0) return false;
  const int tmp3 = compare_seq(b1->kvec(ZDIR).begin(), b1->kvec(ZDIR).end(),
                               b2->kvec(ZDIR).begin(), b2->kvec(ZDIR).end());
  return (tmp3 == 0);
}




//==============================================================================
void LRSpline3DUtils::distributeDataPoints(LRSplineVolume* vol, 
                                           vector<double>& points, 
                                           bool add_distance_field, 
                                           bool primary_points) 
//==============================================================================
{
  int dim = vol->dimension();
  int del = dim+3;                   // Number of entries for each point
                                     // at this stage there is no distance info...
  size_t nmb = (int)points.size()/del;  // Number of data points

  // Erase point information in the elements
  for (LRSplineVolume::ElementMap::const_iterator it = vol->elementsBegin();
       it != vol->elementsEnd(); ++it)
    it->second->eraseDataPoints();
#if 1
  // @obar: This is a simplified version of the 2D implementation, which is slower, but 
  // probably OK if the initial volume has few elements.
  // Can generalise the 2D version if a long time is spent here.
  for (size_t ix=0; ix!=nmb; ix++) {
    // Fetch associated element
    Element3D* elem = vol->coveringElement(points[del*ix],points[del*ix+1],points[del*ix+2]);
    if (add_distance_field)
      elem->addDataPoints(points.begin()+del*ix, points.begin()+del*(ix+1), del, false);
    else
      elem->addDataPoints(points.begin()+del*ix, points.begin()+del*(ix+1), false);
  }
#endif
#if 0
  // Sort the points according to the u-parameter
  qsort(&points[0], nmb, del*sizeof(double), compare_u_par);

  // Get all knot values in the u-direction
  const double* const uknots_begin = vol->mesh().knotsBegin(XDIR);
  const double* const uknots_end = vol->mesh().knotsEnd(XDIR);
  int nmb_knots_u = vol->mesh().numDistinctKnots(XDIR);
  const double* knotu;

  // Construct mesh of element pointers
  vector<Element3D*> elements;
  vol->constructElementMesh(elements);

  // Get all knot values in the v-direction
  const double* const vknots_begin = vol->mesh().knotsBegin(YDIR);
  const double* const vknots_end = vol->mesh().knotsEnd(YDIR);
  int nmb_knots_v = vol->mesh().numDistinctKnots(YDIR);
  const double* knotv;

  // Get all knot values in the w-direction
  const double* const wknots_begin = vol->mesh().knotsBegin(ZDIR);
  const double* const wknots_end = vol->mesh().knotsEnd(ZDIR);
  int nmb_knots_w = vol->mesh().numDistinctKnots(ZDIR);
  const double* knotw;

  // Traverse points and divide them according to their position in the
  // u direction
  int ki, kj, kr;
  int pp0, pp1;
  for (ki=0, pp0=0, knotu=uknots_begin, ++knotu; knotu!= uknots_end; 
       ++knotu, ++ki)
    {
      
      for (pp1=pp0; pp1<(int)points.size() && points[pp1] < (*knotu); pp1+=del);
      if (knotu+1 == uknots_end)
	pp1 = (int)points.size();

      // Sort the current sub set of points according to the v-parameter
      qsort(&points[0]+pp0, (pp1-pp0)/del, del*sizeof(double), compare_v_par);

      // Traverse points and divide them according to their position in the
      // v direction
      int pp2, pp3;
      for (kj=0, pp2=pp0, knotv=vknots_begin, ++knotv; knotv!=vknots_end; 
	   ++knotv, ++kj)
	{
	  for (pp3=pp2; pp3<pp1 && points[pp3+1] < (*knotv); pp3 += del);
	  if (knotv+1 == vknots_end)
	    pp3 = pp1;
	  
	  // Sort the current sub set of points according to the w-parameter
	  qsort(&points[0]+pp2, (pp3-pp2)/del, del*sizeof(double), compare_w_par);

	  // Traverse the relevant points and store them in the associated element
	  // Note that an extra entry will be added for each point to allow for
	  // storing the distance between the point and the surface
	  int pp4, pp5;
	  for (kr=0, pp4=pp2, knotw=wknots_begin, ++knotw; knotw!=wknots_end; 
	       ++knotw, ++kr)
	    {
	      for (pp5=pp4; pp5<pp3 && points[pp5+2] < (*knotw); pp5 += del);
	      if (knotw+1 == wknots_end)
		pp5 = pp3;
	      
	      // Fetch associated element
	      Element3D* elem = elements[(kr*(nmb_knots_v-1)+kj)*(nmb_knots_u-1)+ki];

	      if (elem == NULL)
		{
		  std::cout << "Missing element pointer" << std::endl;
		  std::cout << nmb_knots_u << " " << nmb_knots_v << " " << nmb_knots_w << std::endl;
		  std::cout << ki << " " << kj << " " << kr << std::endl;
		  int ix = 0.5*(pp4+pp5)*del;
		  elem = vol->coveringElement(points[ix], points[ix+1],
					      points[ix+2]);
		}
	      if (add_distance_field)
		elem->addDataPoints(points.begin()+pp4, points.begin()+pp5, 
				    del, false);
	      else
		elem->addDataPoints(points.begin()+pp4, points.begin()+pp5, 
				    false);
	      pp4 = pp5;
	    }
	  pp2 = pp3;
	}
      pp0 = pp1;
    }
#endif  
}

//==============================================================================
void LRSpline3DUtils::evalAllBSplines(const vector<LRBSpline3D*>& bsplines,
				      double upar, double vpar, double wpar,
				      bool u_at_end, bool v_at_end, bool w_at_end,
				      vector<double>& result)
//==============================================================================
{
  size_t bsize = bsplines.size();
  result.resize(bsize);
  vector<double> val(3*bsize);
  for (size_t ki=0; ki<bsize; ++ki)
    {
      size_t kj;
      const BSplineUniLR* uni1 =  bsplines[ki]->getUnivariate(XDIR);
      const BSplineUniLR* uni2 =  bsplines[ki]->getUnivariate(YDIR);
      const BSplineUniLR* uni3 =  bsplines[ki]->getUnivariate(ZDIR);
      for (kj=0; kj<ki; ++kj)
	if (uni1 == bsplines[kj]->getUnivariate(XDIR))
	  break;
      if (kj < ki)
	val[ki] = val[kj];
      else
	val[ki] = uni1->evalBasisFunction(upar, 0, u_at_end);

      for (kj=0; kj<ki; ++kj)
	if (uni2 == bsplines[kj]->getUnivariate(YDIR))
	  break;
      if (kj < ki)
	val[bsize+ki] = val[bsize+kj];
      else
	val[bsize+ki] = 
	  uni2->evalBasisFunction(vpar, 0, v_at_end);

      for (kj=0; kj<ki; ++kj)
	if (uni3 == bsplines[kj]->getUnivariate(ZDIR))
	  break;
      if (kj < ki)
	val[2*bsize+ki] = val[2*bsize+kj];
      else
	val[2*bsize+ki] = 
	  uni3->evalBasisFunction(wpar, 0, w_at_end);

      result[ki] = val[ki]*val[bsize+ki]*val[2*bsize+ki];
    }
 }

//==============================================================================
void LRSpline3DUtils::evalAllBSplines2(const vector<LRBSpline3D*>& bsplines,
				       double upar, double vpar, double wpar,
				       bool u_at_end, bool v_at_end, bool w_at_end,
				       double* result, double* val)
//==============================================================================
{
  size_t bsize = bsplines.size();
  // result.resize(bsize);
  // vector<double> val(3*bsize);
  for (size_t ki=0; ki<bsize; ++ki)
    {
      size_t kj;
       BSplineUniLR* uni1 =  bsplines[ki]->getUnivariate(XDIR);
       BSplineUniLR* uni2 =  bsplines[ki]->getUnivariate(YDIR);
       BSplineUniLR* uni3 =  bsplines[ki]->getUnivariate(ZDIR);
      for (kj=0; kj<ki; ++kj)
	if (uni1 == bsplines[kj]->getUnivariate(XDIR))
	  break;
      if (kj < ki)
	val[ki] = val[kj];
      else
	val[ki] = uni1->evalBasisFunction(upar, 0, u_at_end);

      for (kj=0; kj<ki; ++kj)
	if (uni2 == bsplines[kj]->getUnivariate(YDIR))
	  break;
      if (kj < ki)
	val[bsize+ki] = val[bsize+kj];
      else
	val[bsize+ki] = 
	  uni2->evalBasisFunction(vpar, 0, v_at_end);

      for (kj=0; kj<ki; ++kj)
	if (uni3 == bsplines[kj]->getUnivariate(ZDIR))
	  break;
      if (kj < ki)
	val[2*bsize+ki] = val[2*bsize+kj];
      else
	val[2*bsize+ki] = 
	  uni3->evalBasisFunction(wpar, 0, w_at_end);

      result[ki] = val[ki]*val[bsize+ki]*val[2*bsize+ki];
    }
 }

//==============================================================================
void LRSpline3DUtils::evalAllBSplinePos(const vector<LRBSpline3D*>& bsplines,
					double upar, double vpar, double wpar,
					bool u_at_end, bool v_at_end, bool w_at_end,
					vector<Point>& result)
//==============================================================================
{
  size_t bsize = bsplines.size();
  result.resize(bsize);
  vector<double> val(3*bsize);
  for (size_t ki=0; ki<bsize; ++ki)
    {
      size_t kj;
      const BSplineUniLR* uni1 =  bsplines[ki]->getUnivariate(XDIR);
      const BSplineUniLR* uni2 =  bsplines[ki]->getUnivariate(YDIR);
      const BSplineUniLR* uni3 =  bsplines[ki]->getUnivariate(ZDIR);
      for (kj=0; kj<ki; ++kj)
	if (uni1 == bsplines[kj]->getUnivariate(XDIR))
	  break;
      if (kj < ki)
	val[ki] = val[kj];
      else
	val[ki] = uni1->evalBasisFunction(upar, 0, u_at_end);

      for (kj=0; kj<ki; ++kj)
	if (uni2 == bsplines[kj]->getUnivariate(YDIR))
	  break;
      if (kj < ki)
	val[bsize+ki] = val[bsize+kj];
      else
	val[bsize+ki] = 
	  uni2->evalBasisFunction(vpar, 0, v_at_end);

      for (kj=0; kj<ki; ++kj)
	if (uni3 == bsplines[kj]->getUnivariate(ZDIR))
	  break;
      if (kj < ki)
	val[2*bsize+ki] = val[2*bsize+kj];
      else
	val[2*bsize+ki] = 
	  uni3->evalBasisFunction(wpar, 0, w_at_end);

      result[ki] = val[ki]*val[bsize+ki]*val[2*bsize+ki]*
	bsplines[ki]->coefTimesGamma();
    }
 }

}; // end namespace Go
