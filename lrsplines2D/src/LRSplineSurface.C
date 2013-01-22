#include <stdexcept>
#include <iostream> // @@ debug
#include <fstream>
#include <iterator> // @@ debug - remove
//#include <chrono>   // @@ debug
#include <set>
#include <tuple>
#include "GoTools/utils/checks.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/Mesh2DUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"
#include "GoTools/utils/StreamUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2DUtils.h"
#include "GoTools/geometry/SplineCurve.h"
//#include "GoTools/lrsplines2D/PlotUtils.h" // @@ only for debug

using std::vector;
using std::istream;
using std::ostream;
using std::get;
using std::pair;
using std::for_each;

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
      LRSplineUtils::update_elements_with_single_bspline(tmp, emap, 
							 m, false);
    }

  return emap;
};

//==============================================================================
LRSplineSurface::LRSplineSurface(SplineSurface *surf, double knot_tol)
//==============================================================================
: knot_tol_(knot_tol), rational_(surf->rational()),
  mesh_(surf->basis_u().begin(), surf->basis_u().end(),
	surf->basis_v().begin(), surf->basis_v().end())
{
  std::vector<int> knot_ixs_u = init_knot_indices(mesh_, XFIXED);
  std::vector<int> knot_ixs_v = init_knot_indices(mesh_, YFIXED);

  int deg_u = surf->order_u() - 1;
  int deg_v = surf->order_v() - 1;
  int coefs_u = surf->numCoefs_u();
  int coefs_v = surf->numCoefs_v();
  std::vector<double>::iterator rcoefs = surf->rcoefs_begin();
  std::vector<double>::iterator coefs = surf->coefs_begin();
  int dim = surf->dimension();
  int kdim = (rational_) ? dim + 1 : dim;
  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    for (int u_ix = 0; u_ix != coefs_u; ++u_ix, coefs+=dim, rcoefs += kdim) {
      double rat = (rational_) ? rcoefs[dim] : 1.0;
      shared_ptr<LRBSpline2D> b(new LRBSpline2D(Point(coefs, coefs + dim),
						rat,
						deg_u,
						deg_v,
						knot_ixs_u.begin() + u_ix,
						knot_ixs_v.begin() + v_ix,
						1.0, &mesh_, rational_));
      bsplines_[generate_key(*b, mesh_)] = b;
    }
  }
  emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
LRSplineSurface::LRSplineSurface(double knot_tol, bool rational,
				 Mesh2D& mesh, 
				 vector<shared_ptr<LRBSpline2D> > b_splines)
//==============================================================================
  : knot_tol_(knot_tol), rational_(rational), mesh_(mesh)
{
  for (size_t ki=0; ki<b_splines.size(); ++ki)
    bsplines_[generate_key(*b_splines[ki], mesh_)] = b_splines[ki];

  emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
LRSplineSurface::LRSplineSurface(const LRSplineSurface& rhs) 
//==============================================================================
  : knot_tol_(rhs.knot_tol_), rational_(rhs.rational_), 
    mesh_(rhs.mesh_)
{
  // Clone LR B-splines
  BSplineMap::const_iterator curr = rhs.basisFunctionsBegin();
  BSplineMap::const_iterator end = rhs.basisFunctionsEnd();
  for (; curr != end; ++curr)
    {
      shared_ptr<LRBSpline2D> b(new LRBSpline2D(*curr->second));
      bsplines_[generate_key(*b, mesh_)] = b;
    }

  // The ElementMap has to be generated and cannot be copied directly, since it
  // contains raw pointers.  
  emap_ = construct_element_map_(mesh_, bsplines_);
}

//==============================================================================
void LRSplineSurface::swap(LRSplineSurface& rhs)
//==============================================================================
{
  std::swap(knot_tol_,    rhs.knot_tol_);
  std::swap(rational_,    rhs.rational_);
  std::swap(mesh_    ,    rhs.mesh_);
  std::swap(bsplines_,    rhs.bsplines_);
  std::swap(emap_    ,    rhs.emap_);
}

//==============================================================================
void  LRSplineSurface::read(istream& is)
//==============================================================================
{
  LRSplineSurface tmp;
  int rat = -1;
  object_from_stream(is, rat);
  rational_ = (rat == 1);

  // reading knot tolerances and the mesh
  object_from_stream(is, tmp.knot_tol_);
  object_from_stream(is, tmp.mesh_);

  // Reading all basis functions
  int num_bfuns;
  object_from_stream(is, num_bfuns);
  for (int i = 0; i != num_bfuns; ++i) {
    shared_ptr<LRBSpline2D> b(new LRBSpline2D());
    object_from_stream(is, *b);
    // We set the global mesh in the b basis function.
    b->setMesh(&tmp.mesh_);
#if 1
    tmp.bsplines_[generate_key(*b, tmp.mesh_)] = b;
#else
    BSKey key = generate_key(*b, tmp.mesh_);
    tmp.bsplines_.insert(std::make_pair(key, b));
#endif
  }

  // Reconstructing element map
  tmp.emap_ = construct_element_map_(tmp.mesh_, tmp.bsplines_);

  tmp.rational_ = rational_;

  this->swap(tmp);

  auto it = bsplines_.begin();
  while (it != bsplines_.end())
    {
      it->second->setMesh(&mesh_);
      ++it;
    }

#if 0//ndef NDEBUG
  {
    vector<LRBSpline2D*> bas_funcs;
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    puts("Remove when done debugging!");
    vector<Element2D*> elems;
    for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
    {
	elems.push_back(((*iter).second.get()));
    }
    puts("Remove when done debugging!");
  }
#endif

}

//==============================================================================
void LRSplineSurface::write(ostream& os) const
//==============================================================================
{
  int rat = (rational_) ? 1 : 0;
  object_to_stream(os, rat);
  object_to_stream(os, knot_tol_);
  object_to_stream(os, '\n');
  object_to_stream(os, mesh_);

  object_to_stream(os, bsplines_.size());
  object_to_stream(os, '\n');
  for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) 
    {
      object_to_stream(os, *(b->second));
      object_to_stream(os, '\n');
    }

  // NB: 'emap_' is not saved (contains raw pointers to other data).  
  // Instead, it will be regenerated the LRSplineSurface is read().
}

//==============================================================================
const LRSplineSurface::ElementMap::value_type& 
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

  const LRSplineSurface::ElemKey key = 
    {mesh_.knotsBegin(XFIXED)[ucorner], mesh_.knotsBegin(YFIXED)[vcorner]};
  const auto el = emap_.find(key);
  assert(el != emap_.end());
  return *el;
}


//==============================================================================
vector<LRBSpline2D*> LRSplineSurface::basisFunctionsWithSupportAt(double u, double v) const
//==============================================================================
{
  vector<LRBSpline2D*> support_functions;
  auto it = coveringElement(u, v);
  vector<LRBSpline2D*>::const_iterator first = it.second->supportBegin();
  vector<LRBSpline2D*>::const_iterator last = it.second->supportEnd();
  int ki=0;
  for (; first != last; ++first, ++ki)
    {
      support_functions.push_back(*first);
    }
  
  return support_functions;
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
#ifndef NDEBUG
  if (rational_)
    MESSAGE("Refining a rational surface does not seem to work, yet.");
#endif

  refine(ref.d, ref.kval, ref.start, ref.end, ref.multiplicity, absolute);
}

//==============================================================================
void LRSplineSurface::refine(Direction2D d, double fixed_val, double start, 
			     double end, int mult, bool absolute)
//==============================================================================
{
#ifndef NDEBUG
  if (rational_)
    MESSAGE("Refining a rational surface does not seem to work, yet.");
#endif

#if 0//ndef NDEBUG
  vector<LRBSpline2D*> bas_funcs;
  for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
    {
      bas_funcs.push_back((*iter).second.get());
    }
  vector<Element2D*> elems;
  for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
  {
      elems.push_back(((*iter).second.get()));
  }
  puts("Remove when done debugging!");
#endif

  const auto indices = // tuple<int, int, int, int>
  LRSplineUtils::refine_mesh(d, fixed_val, start, end, mult, absolute, 
			     degree(d), knot_tol_, mesh_, bsplines_);

  // insert newly created elements to emap (unless refinement was on border, in which case no new element
  // could possibly be created
  const int prev_ix = get<0>(indices);
  const int fixed_ix = get<1>(indices); // Index of fixed_val in global knot vector.
  const int start_ix = get<2>(indices); // Index of start (of segment to insert) in global knot vector.
  const int end_ix   = get<3>(indices); // Index of end (of segment to insert) in global knot vector.

  // Collect pointers to affected bsplines
  std::set<LRBSpline2D*> all_bsplines;
  for (int i = start_ix; i != end_ix; ++i) {
    // Check if the specified element exists in 'emap'
    int u_ix = (d == XFIXED) ? prev_ix : i;
    int v_ix = (d == YFIXED) ? prev_ix : i;
    ElementMap::key_type key = {mesh_.kval(XFIXED, u_ix),
				mesh_.kval(YFIXED, v_ix)};
    auto it = emap_.find(key);
    if (it != emap_.end())
      {
	// The element exists. Collect bsplines
	all_bsplines.insert(it->second->supportBegin(), it->second->supportEnd());
      }
  }
  vector<LRBSpline2D*> bsplines_affected(all_bsplines.begin(), all_bsplines.end());

#if 0//ndef NDEBUG
    bas_funcs.clear();
    for (auto iter = bsplines_.begin(); iter != bsplines_.end(); ++iter)
      {
	bas_funcs.push_back((*iter).second.get());
      }
    elems.clear();
    for (auto iter = emap_.begin(); iter != emap_.end(); ++iter)
      {
	  elems.push_back(((*iter).second.get()));
      }
    puts("Remove when done debugging!");
#endif

  // Cannot remove the bsplines from the global array at this stage since we operate
  // with pointers to it. When a bspline is split, the origin is removed from the
  // array after all pointers are updated and the the bspline is allowed to die.
  // Iteratively split affected LRBSpline2Ds
  // @@@ VSK. Will pointers to other entities in bsplines_ which are not
  // affected remain valid after removing and adding elements? If not, this
  // combination of objects and pointers will not work.
  LRSplineUtils::iteratively_split2(bsplines_affected, mesh_, bsplines_); 

#if 0//ndef NDEBUG
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
	      MESSAGE("Element with no support functions!");
      }
    puts("Remove when done debugging!");
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

	int u_ix = (d == XFIXED) ? fixed_ix : i;
	int v_ix = (d == YFIXED) ? fixed_ix : i;
	ElementMap::key_type key = {mesh_.kval(XFIXED, u_ix),
				    mesh_.kval(YFIXED, v_ix)};
	auto it = emap_.find(key);
	if (it2 != emap_.end())
	  {
	    // Update size of existing element
	    Mesh2DIterator m(mesh_, u_ix2, v_ix2);
	    it2->second->setUmax(mesh_.kval(XFIXED, (*m)[2]));
	    it2->second->setVmax(mesh_.kval(YFIXED, (*m)[3]));

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
	    shared_ptr<Element2D> elem(new Element2D(mesh_.kval(XFIXED, (*m)[0]),
						     mesh_.kval(YFIXED, (*m)[1]),
						     mesh_.kval(XFIXED, (*m)[2]),
						     mesh_.kval(YFIXED, (*m)[3])));
	    emap_[key] = elem;
	    //auto it3 = emap_.find(key);


	    // Set LRBsplines
	    for (size_t kb=0; kb<bsplines_affected.size(); ++kb)
	    {
		if (bsplines_affected[kb]->overlaps(elem.get()))
		{
		    elem->addSupportFunction(bsplines_affected[kb]);
		    bsplines_affected[kb]->addSupport(elem.get());
		}
	    }

	  }

      }
    }
  }

}

//==============================================================================
void LRSplineSurface::refine(const vector<Refinement2D>& refs, 
			     bool absolute)
//==============================================================================
{
#ifndef NDEBUG
  if (rational_)
    MESSAGE("Refining a rational surface does not seem to work, yet.");
  std::wcout << "Inserting refinements into mesh." << std::endl;
#endif

  for (size_t i = 0; i != refs.size(); ++i) {
    const Refinement2D& r = refs[i];
       LRSplineUtils::refine_mesh(r.d, 
				  r.kval, 
				  r.start, 
				  r.end, 
				  r.multiplicity, 
				  absolute,
				  degree(r.d), 
				  knot_tol_, 
				  mesh_, 
				  bsplines_);
  }


  std::wcout << "Preparing for iterative splitting." << std::endl;
  vector<shared_ptr<LRBSpline2D> > affected;
  affected.reserve(bsplines_.size());
  for_each(bsplines_.begin(), bsplines_.end(), [&](const BSplineMap::value_type& b) {
      // @@@ VSK. This is maybe the place to remove element information from the bsplines?
      affected.push_back(b.second);
    });
  
  // @@@ VSK. In this case, we should not bother about splitting elements. They will
  // be regenerated later. Thus, the bsplines should NOT be updated with elements during
  // splitting
  // The bsplines should not have any pointers to elements. They will be set later
  std::wcout << "Iteratively splitting." << std::endl;
  LRSplineUtils::iteratively_split(affected, mesh_);
  bsplines_.clear();

  std::wcout << "Splitting finished, now inserting resulting functions" << std::endl;
  // The bsplines are checked for duplicates and inserted in the global bspline map
  for_each(affected.begin(), affected.end(), [&](shared_ptr<LRBSpline2D> b) {
      LRSplineUtils::insert_basis_function(b, mesh_, bsplines_);
    });

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

  std::wcout << "Finally, reconstructing element map." << std::endl;
  emap_ = construct_element_map_(mesh_, bsplines_); // reconstructing the emap once at the end
  std::wcout << "Refinement now finished. " << std::endl;
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
void LRSplineSurface::to3D()
//==============================================================================
{
  if (dimension() != 1) 
    THROW("Member method 'to3D()' only applies to one-dimensional LR-splines");
  if (degree(XFIXED) == 0 || degree(YFIXED) == 0) 
    THROW("Cannot convert a 0-degree spline to 3D.");

  for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) {
    const double x = LRSplineUtils::compute_greville(b->second->kvec(XFIXED), 
						     mesh().knotsBegin(XFIXED));
    const double y = LRSplineUtils::compute_greville(b->second->kvec(YFIXED), 
						      mesh().knotsBegin(YFIXED));
    const double z_gamma = b->second->coefTimesGamma()[0];
    const double gamma = b->second->gamma();
    b->second->coefTimesGamma() = Point(x*gamma, y*gamma, z_gamma);
    //wcout << b.second.coefTimesGamma() << std::endl;
  }
}


//==============================================================================
bool LRSplineSurface::rational() const
//==============================================================================
{
  return rational_;
}


//==============================================================================
void LRSplineSurface::expandToFullTensorProduct()
//==============================================================================
{
  std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - copying mesh..." << std::endl;
  Mesh2D tensor_mesh = mesh_;
  
  std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - setting uniform meshlines..." << std::endl;
  const vector<int> xmults = LRSplineUtils::set_uniform_meshlines(XFIXED, 
								  tensor_mesh);
  const vector<int> ymults = LRSplineUtils::set_uniform_meshlines(YFIXED, 
								  tensor_mesh);

  BSplineMap tensor_bsplines;
  std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - identify elements from mesh..." << std::endl;
  ElementMap emap = LRSplineUtils::identify_elements_from_mesh(tensor_mesh);
  std::wcout << "Size of emap: " << emap.size() << std::endl;
  std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - splitting up basis functions..." << std::endl;
  // splitting up basis functions
  for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) {
    LRSplineUtils::tensor_split(b->second, 
				xmults, 
				ymults, 
				tensor_mesh,
				tensor_bsplines);
  }

  // registering all the produced functions with the elements
  std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - registering produced functions..." << std::endl;
  std::wcout << "Number of basis functions: " << tensor_bsplines.size() << std::endl;
  std::wcout << "Number of elements: "<< tensor_mesh.numDistinctKnots(XFIXED)-1 << " x " ;
  std::wcout << tensor_mesh.numDistinctKnots(YFIXED)-1 << std::endl;

  // @@@ VSK. Use information in the LRB-splines or regenerate all elements ?
  for (auto b = tensor_bsplines.begin(); b != tensor_bsplines.end(); ++b)  {
    LRSplineUtils::update_elements_with_single_bspline(b->second.get(), emap, 
						       tensor_mesh, false);
  }

  std::wcout << "LRSplineSurface::ExpandToFullTensorProduct() - swapping and exiting." << std::endl;
  mesh_.swap(tensor_mesh);
  bsplines_.swap(tensor_bsplines);
  emap_.swap(emap);
}


//==============================================================================
Point LRSplineSurface::operator()(double u, double v, int u_deriv, int v_deriv) const
//==============================================================================
{
#ifndef NDEBUG
  if (u_deriv + v_deriv > 1)
    {
      MESSAGE("Currently the sum of the derivatives should be at most 1.");
    }
#endif

  // const bool u_on_end = (u == mesh_.maxParam(XFIXED));
  // const bool v_on_end = (v == mesh_.maxParam(YFIXED));
  // vector<LRBSpline2D*> covering_B_functions = 
  //   basisFunctionsWithSupportAt(u, v);
  auto it = coveringElement(u, v);
  const vector<LRBSpline2D*> covering_B_functions = it.second->getSupport();

  Point result(this->dimension()); 
  result.setValue(0.0); // will be initialized to 0, with the correct dimension

  // loop over LR B-spline functions
  int ki=0;
  int nmb_b = (int)covering_B_functions.size();
  double denom = (rational_) ? 0.0 : 1.0;
  vector<double> basis_vals((u_deriv+1)*(v_deriv+1), 0.0); // To be used for rational cases, needed for derivs.
  for (auto b = covering_B_functions.begin(); 
       b != covering_B_functions.end(); ++b, ++ki) 
    {
      const bool u_on_end = (u == (*b)->umax());
      const bool v_on_end = (v == (*b)->vmax());
      // The b-function contains the coefficient.
      if (!rational_)
	{
	  result += (*b)->eval(u, 
			       v, 
			       mesh_.knotsBegin(XFIXED), 
			       mesh_.knotsBegin(YFIXED), 
			       u_deriv, 
			       v_deriv, 
			       u_on_end, 
			       v_on_end);
	}
      else
	{
	  
//	  for (size_t ki = 0; ki < 
	  double basis_val = (*b)->evalBasisFunction(u, 
						     v, 
						     mesh_.knotsBegin(XFIXED), 
						     mesh_.knotsBegin(YFIXED), 
						     u_deriv, 
						     v_deriv, 
						     u_on_end, 
						     v_on_end);
	  double gamma = (*b)->gamma();
	  double weight = (*b)->weight();
	  Point coef = (*b)->Coef();

	  result += coef*weight*basis_val;

	  denom += weight*basis_val;

#ifndef NDEBUG
	  if (u_deriv > 0 || v_deriv > 0)
	    {
	      MESSAGE("Do not think that rational derivs are supported yet.");
	      denom = 1.0;
	    }
	  std::cout << "denom: " << denom << std::endl;
	  // if (rat_den == 0.0)
	  //   rat_den = 1.0;
#endif

	}
    }

  if (rational_ && (u_deriv + v_deriv == 1))
    denom = denom*denom;

  return result/denom;
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
  bool LRSplineSurface::inDomain(double u, double v) const
  //===========================================================================
  {
    if (u < startparam_u() || u > endparam_u())
	return false;
    if (v < startparam_v() || v > endparam_v())
	return false;

    return true;
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
    pt = operator()(upar, vpar, 0, 0);
  }

  //===========================================================================
void LRSplineSurface::normal(Point& pt, double upar, double vpar) const
  //===========================================================================
  {
    Point pt_der1 = operator()(upar, vpar, 1, 0);
    Point pt_der2 = operator()(upar, vpar, 0, 1);

    pt = pt_der1.cross(pt_der2);
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
    // @@sbr201301 Test this function!
    MESSAGE("LRSplineSurface::point() not testet yet");
    int vec_size = (derivs + 1)*(derivs + 2)/2;
    pts.resize(vec_size);

    // This is not the most efficient approach, should be faster to
    // evaluate basis functions once. Only a first implementation.
    int cntr = 0;
    for (size_t kj = 0; kj < derivs; ++kj)
	for (size_t ki = 0; ki < kj; ++ki, ++cntr)
	    pts[cntr] = operator()(upar, vpar, kj-ki, ki);
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

     // Define possible new knotlines
     // Note that the new knot lines must be longer than the size of the
     // sub surface to avoid LR B-splines partly overlapping the sub domain
     // @@@ VSK. Could the extension be smaller than prescribed here?
     vector<Refinement2D> refs(4);
     double umin = mesh_.kval(XFIXED, std::max(ix1 - deg1, 0));
     double umax = mesh_.kval(XFIXED, std::min(ix2 + deg1, nmb1-1));
     double vmin = mesh_.kval(YFIXED, std::max(iy1 - deg2, 0));
     double vmax = mesh_.kval(YFIXED, std::min(iy2 + deg2, nmb2-1));
     refs[0].setVal(from_upar, vmin, vmax, XFIXED, deg1+1);
     refs[1].setVal(to_upar, vmin, vmax, XFIXED, deg1+1);
     refs[2].setVal(from_vpar, umin, umax, YFIXED, deg2+1);
     refs[3].setVal(to_vpar, umin, umax, YFIXED, deg2+1);
     
     // Perform refinement
     sf->refine(refs, true);

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
     vector<shared_ptr<LRBSpline2D> > b_splines = 
       sf->collect_basis(iu1, iu2, iv1, iv2);

     // Copy LR B-splines and update indices to the new domain
     // It is not absolutely necessary to make copies since the intermediate
     // surface will die after this function, but it is done for consistency
     vector<shared_ptr<LRBSpline2D> > b_splines2(b_splines.size());
     for (size_t ki=0; ki<b_splines.size(); ++ki)
       {
	 b_splines2[ki] = shared_ptr<LRBSpline2D>(new LRBSpline2D(*b_splines[ki]));
	 b_splines2[ki]->setMesh(sub_mesh.get());
	 b_splines2[ki]->subtractKnotIdx(iu1, iv1);
       }
     

     // Create sub surface
     surf = new LRSplineSurface(knot_tol_, rational_, *sub_mesh, b_splines2);
     
    return surf;
  }

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
    MESSAGE("LRSplineSurface::getBoundaryInfo() not implemented yet");
  }

  //===========================================================================
  void LRSplineSurface::turnOrientation()
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::turnOrientation() not implemented yet");
  }

  //===========================================================================
  void LRSplineSurface::swapParameterDirection()
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::swapParameterDirection() not implemented yet");
  }

  //===========================================================================
  void LRSplineSurface::reverseParameterDirection(bool direction_is_u)
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::reverseParameterDirection() not implemented yet");
  }

  //===========================================================================
  double LRSplineSurface::area(double tol) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::area() not implemented yet");
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
	for (size_t ki = 0; ki < 4; ++ki)
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
      return constParamCurve(mesh_.kval(d, ix), d == XFIXED);
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
	{
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
      size_t km = 1;
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
	THROW("edgeCurve:: There is no such basis function.");

      // Fetch coefficient
      Point cf = bm->second->coefTimesGamma();
      coefs.insert(coefs.end(), cf.begin(), cf.end());
    }

  // Define spline curve
  // First compute the real knot values
  int nmbcf = knot_idx.size() - deg2 - 1;
  vector<double> knots(knot_idx.size());
  for (k1=0; k1<knot_idx.size(); ++k1)
    knots[k1] = mesh_.kval(d2, knot_idx[k1]);

  edgcv = new SplineCurve(nmbcf, deg2+1, knots.begin(), coefs.begin(),
			  dimension(), rational_);

  return edgcv;
}

//===========================================================================
SplineCurve*
LRSplineSurface::constParamCurve (double parameter,
				  bool pardir_is_u) const
//===========================================================================
{
    MESSAGE("LRSplineSurface::constParamCurves() not implemented yet");
    return NULL;
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
    MESSAGE("LRSplineSurface::nextSegmentVal() not implemented yet");
    return 0.0;
  }

//===========================================================================
LRSplineSurface* LRSplineSurface::mirrorSurface(const Point& pos, 
						const Point& norm) const
//===========================================================================
  {
    MESSAGE("LRSplineSurface::mirrorSurface() not implemented yet");
    return NULL;
  }

//===========================================================================
vector<shared_ptr<LRBSpline2D> > 
LRSplineSurface::collect_basis(int from_u, int to_u, 
			       int from_v, int to_v) const
//===========================================================================
{
  vector<shared_ptr<LRBSpline2D> > b_splines;
  for (auto it=bsplines_.begin(); it!=bsplines_.end(); ++it)
    {
      if (interval_overlap(it->second->suppMin(XFIXED), 
			   it->second->suppMax(XFIXED), 
			   from_u, to_u) &&
	  interval_overlap(it->second->suppMin(YFIXED), 
			   it->second->suppMax(YFIXED), 
			   from_v, to_v))
	b_splines.push_back(it->second);
    }
  return b_splines;
}

} // end namespace Go


