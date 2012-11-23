#include <stdexcept>
 #include <iostream> // @@ debug
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
      const LRBSpline2D* tmp = &(b_it->second);
      // The connection between the bspline and the elements in its support is set
      // @@@ VSK. Must get rid of the const cast. The split function in Element2D
      // is probably what creates the need. This function is likely to disappear
      // when the integration is completed.
      LRSplineUtils::update_elements_with_single_bspline(const_cast<LRBSpline2D*>(tmp), emap, 
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
  std::vector<double>::iterator coefs = rational_ ?
    surf->rcoefs_begin() : surf->coefs_begin();
  int dimension = surf->dimension();
  int kdim = dimension + rational_;
  for (int v_ix = 0; v_ix != coefs_v; ++v_ix)  {
    for (int u_ix = 0; u_ix != coefs_u; ++u_ix, coefs+=kdim) {
      LRBSpline2D b(Point(coefs, coefs + kdim),
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

//==============================================================================
  LRSplineSurface::LRSplineSurface(const LRSplineSurface& rhs) 
//==============================================================================
    : knot_tol_(rhs.knot_tol_), mesh_(rhs.mesh_), bsplines_(rhs.bsplines_),
      emap_(construct_element_map_(mesh_, bsplines_))
{
  // The ElementMap has to be generated and cannot be copied directly, since it
  // contains raw pointers.  Hence the call to 'construct_element_map_' in the
  // initialization
}

//==============================================================================
void LRSplineSurface::swap(LRSplineSurface& rhs)
//==============================================================================
{
  std::swap(knot_tol_,    rhs.knot_tol_);
  std::swap(mesh_    ,    rhs.mesh_);
  std::swap(bsplines_,    rhs.bsplines_);
  std::swap(emap_    ,    rhs.emap_);
}

//==============================================================================
void  LRSplineSurface::read(istream& is)
//==============================================================================
{
  LRSplineSurface tmp;

  // reading knot tolerances and the mesh
  object_from_stream(is, tmp.knot_tol_);
  object_from_stream(is, tmp.mesh_);

  // Reading all basis functions
  int num_bfuns;
  LRBSpline2D b;
  object_from_stream(is, num_bfuns);
  for (int i = 0; i != num_bfuns; ++i) {
    object_from_stream(is, b);
    // @@@ VSK. Has to add the mesh pointer to the LRBspline
    tmp.bsplines_[generate_key(b, tmp.mesh_)] = b;
  }

  // Reconstructing element map
  tmp.emap_ = construct_element_map_(tmp.mesh_, tmp.bsplines_);

  this->swap(tmp);
}

//==============================================================================
void LRSplineSurface::write(ostream& os) const
//==============================================================================
{
  object_to_stream(os, knot_tol_);
  object_to_stream(os, '\n');
  object_to_stream(os, mesh_);

  object_to_stream(os, bsplines_.size());
  object_to_stream(os, '\n');
  for (auto b = bsplines_.begin(); b != bsplines_.end(); ++b) 
    {
      object_to_stream(os, b->second);
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
  Element2D elem = it.second;
  vector<LRBSpline2D*>::const_iterator first = it.second.supportBegin();
  vector<LRBSpline2D*>::const_iterator last = it.second.supportEnd();
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
  refine(ref.d, ref.kval, ref.start, ref.end, ref.multiplicity, absolute);
}

//==============================================================================
void LRSplineSurface::refine(Direction2D d, double fixed_val, double start, 
			     double end, int mult, bool absolute)
//==============================================================================
{
  const auto indices = 
  LRSplineUtils::refine_mesh(d, fixed_val, start, end, mult, absolute, 
			     degree(d), knot_tol_, mesh_, bsplines_);

  // insert newly created elements to emap (unless refinement was on border, in which case no new element
  // could possibly be created
  const int prev_ix = get<0>(indices);
  const int fixed_ix = get<1>(indices);
  const int start_ix = get<2>(indices);
  const int end_ix   = get<3>(indices);

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
	all_bsplines.insert(it->second.supportBegin(), it->second.supportEnd());
      }
  }
  vector<LRBSpline2D*> bsplines_affected(all_bsplines.begin(), all_bsplines.end());

  // Cannot remove the bsplines from the global array at this stage since we operate
  // with pointers to it. When a bspline is split, the origin is removed from the
  // array after all pointers are updated and the the bspline is allowed to die.
  // Iteratively split affected LRBSpline2Ds
  // @@@ VSK. Will pointers to other entities in bsplines_ which are not
  // affected remain valid after removing and adding elements? If not, this
  // combination of objects and pointers will not work.
  LRSplineUtils::iteratively_split2(bsplines_affected, mesh_, bsplines_); 

  if (fixed_ix > 0 && fixed_ix != mesh_.numDistinctKnots(d)-1) {
    for (int i = start_ix; i != end_ix; ++i) {
      if (mesh_.nu(flip(d), i, fixed_ix, fixed_ix+1) > 0) {
	// this is the lower-left corner of an element bordering our refinement.  Check if it already
	// exists in 'emap', if not, insert it.
	// @@@ VSK. This piece of code must also modify the current
	// element (if prev_ix != fixed_ix) and update bspline pointers
	// in the elements
	int u_ix2 = (d == XFIXED) ? prev_ix : i;
	int v_ix2 = (d == YFIXED) ? prev_ix : i;
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
	    it2->second.setUmax(mesh_.kval(XFIXED, (*m)[2]));
	    it2->second.setVmax(mesh_.kval(YFIXED, (*m)[3]));

	    // Update supported LRBsplines
	    for (size_t kb=0; kb<bsplines_affected.size(); ++kb)
	      {
		if (!bsplines_affected[kb]->overlaps(&it2->second))
		  {
		    it2->second.removeSupportFunction(bsplines_affected[kb]);
		    bsplines_affected[kb]->removeSupport(&it2->second);
		  }
		else
		  it2->second.addSupportFunction(bsplines_affected[kb]);
	      }
	  }

	// Create new element
	if (it == emap_.end())
	  {
	    Mesh2DIterator m(mesh_, u_ix, v_ix);
	    Element2D elem(mesh_.kval(XFIXED, (*m)[0]),
			   mesh_.kval(YFIXED, (*m)[1]),
			   mesh_.kval(XFIXED, (*m)[2]),
			   mesh_.kval(YFIXED, (*m)[3]));

	    // Set LRBsplines
	    for (size_t kb=0; kb<bsplines_affected.size(); ++kb)
	      {
		if (bsplines_affected[kb]->overlaps(&elem))
		  {
		    elem.addSupportFunction(bsplines_affected[kb]);
		    bsplines_affected[kb]->addSupport(&elem);
		  }
	      }

	    emap_[key] = elem;
	  }

      }
    }
  }

  // // Decide exactly which LRBSpline2Ds were affected by this particular insertion, 
  // // remove them from 'bsplines_', and return them in a vector.
  // // @@@ VSK. Assume that the bsplines are collected. We still have to remove them from
  // // bsplines_. This involves a search using the key, but no interference testing is
  // // required.
  // // @@@ VSK. This call should be removed. The content is already performed
  // // (I hope)
  // vector<LRBSpline2D> affected = 
  // LRSplineUtils::collect_and_remove(bsplines_, d, mesh_, fixed_val, 
  // 				    start_ix, end_ix, emap_);
  
  // // Iteratively split affected LRBSpline2Ds
  // // @@@ VSK. Also this is done already
  // LRSplineUtils::iteratively_split(affected, mesh_); 
  
  // @@@ VSK. Should the bsplines be updated with the elements in the previous call
  // or should this be done here? iteratively_split may end up with splitting the 
  // bsplines several times. The elements in the support should be set when all splitting
  // of this B-spline is finished.

  // @@@ VSK. This part should be modified due to information about which
  // bsplines where modified. Only a limited set of elements can possibly
  // be affected
  // @@@ VSK. I think this is done as well, if it works out correctly
  // for_each(affected.begin(), affected.end(), [&](const LRBSpline2D& b) {
  //     const LRBSpline2D* bfun = 
  // 	LRSplineUtils::insert_basis_function(b, mesh_, bsplines_);
  //     LRSplineUtils::update_elements_with_single_bspline(const_cast<LRBSpline2D*>(bfun), emap_, 
  // 							 mesh_, false);
  //   });
}

//==============================================================================
void LRSplineSurface::refine(const vector<Refinement2D>& refs, 
			     bool absolute)
//==============================================================================
{
  std::wcout << "Inserting refinements into mesh." << std::endl;

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
  vector<LRBSpline2D> affected;
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
  for_each(affected.begin(), affected.end(), [&](const LRBSpline2D& b) {
      LRSplineUtils::insert_basis_function(b, mesh_, bsplines_);
    });

  std::wcout << "Finally, reconstructing element map." << std::endl;
  emap_ = construct_element_map_(mesh_, bsplines_); // reconstructing the emap once at the end
  std::wcout << "Refinement now finished. " << std::endl;

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
    const double x = LRSplineUtils::compute_greville(b->second.kvec(XFIXED), 
						     mesh().knotsBegin(XFIXED));
    const double y = LRSplineUtils::compute_greville(b->second.kvec(YFIXED), 
						      mesh().knotsBegin(YFIXED));
    const double z_gamma = b->second.coefTimesGamma()[0];
    const double gamma = b->second.gamma();
    b->second.coefTimesGamma() = Point(x*gamma, y*gamma, z_gamma);
    //wcout << b.second.coefTimesGamma() << std::endl;
  }
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
    LRSplineUtils::update_elements_with_single_bspline(&b->second, emap, 
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

  const bool u_on_end = (u == mesh_.maxParam(XFIXED));
  const bool v_on_end = (v == mesh_.maxParam(YFIXED));
  // vector<LRBSpline2D*> covering_B_functions = 
  //   basisFunctionsWithSupportAt(u, v);
  auto it = coveringElement(u, v);
  const vector<LRBSpline2D*> covering_B_functions = it.second.getSupport();

  Point result(this->dimension()); 
  result.setValue(0.0); // will be initialized to 0, with the correct dimension

  // loop over LR B-spline functions
  int ki=0;
  int nmb_b = (int)covering_B_functions.size();
  for (auto b = covering_B_functions.begin(); 
       b != covering_B_functions.end(); ++b, ++ki) 
    result += (*b)->eval(u, 
			 v, 
			 mesh_.knotsBegin(XFIXED), 
			 mesh_.knotsBegin(YFIXED), 
			 u_deriv, 
			 v_deriv, 
			 u_on_end, 
			 v_on_end);

  return result;
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
  const double gamma = it->second.gamma();
  it->second.coefTimesGamma() = value * gamma;
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
  const double gamma = it->second.gamma();
  it->second.coefTimesGamma() = value * gamma;
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
    BoundingBox box;
    BSplineMap::const_iterator curr = basisFunctionsBegin();
    BSplineMap::const_iterator end = basisFunctionsEnd();
    for (; curr!=end; ++curr)
      {
	const Point coef = curr->second.Coef();
	box.addUnionWith(coef);
      }
    return box;
  }

  //===========================================================================
const RectDomain& LRSplineSurface::parameterDomain() const
  //===========================================================================
  {
    Array<double, 2> ll(mesh_.minParam(YFIXED), mesh_.minParam(XFIXED));
    Array<double, 2> ur(mesh_.maxParam(YFIXED), mesh_.maxParam(XFIXED));
    return RectDomain(ll, ur);
  }

  //===========================================================================
  RectDomain LRSplineSurface::containingDomain() const
  //===========================================================================
  {
    return parameterDomain();
  }

 //===========================================================================
  bool LRSplineSurface::inDomain(double u, double v) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::inDomain() not implemented yet");
    return false;
  }

  //===========================================================================
  Point LRSplineSurface::closestInDomain(double u, double v) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::closestInDomain() not implemented yet");
    Point p;
    return p;
  }

  //===========================================================================
  CurveLoop LRSplineSurface::outerBoundaryLoop(double degenerate_epsilon) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::outerBoundaryLoop() not implemented yet");
    CurveLoop cl;
    return cl;
  }

  //===========================================================================
  vector<CurveLoop> LRSplineSurface::allBoundaryLoops(double degenerate_epsilon) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::allBoundaryLoops() not implemented yet");
    vector<CurveLoop> v_cl;
    return v_cl;
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
    MESSAGE("LRSplineSurface::normal() not implemented yet");
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
    MESSAGE("LRSplineSurface::point() not implemented yet");
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
  vector<shared_ptr<ParamSurface> >
    LRSplineSurface::subSurfaces(double from_upar, double from_vpar,
				 double to_upar, double to_vpar,
				 double fuzzy) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::subSurfaces() not implemented yet");
    vector<shared_ptr<ParamSurface> > res;
    return res;
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
    MESSAGE("LRSplineSurface::closestBoundaryPoint() not implemented yet");
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
  bool LRSplineSurface::isDegenerate(bool& b, bool& r,
				     bool& t, bool& l, double tolerance) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::isDegenerate() not implemented yet");
    return false;
  }

  //===========================================================================
  void LRSplineSurface::getDegenerateCorners(vector<Point>& deg_corners, double tol) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::getDegenerateCorners() not implemented yet");
  }

  //===========================================================================
  void LRSplineSurface::getCornerPoints(vector<pair<Point,Point> >& corners) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::getCornerPoints() not implemented yet");
  }

  //===========================================================================
  vector< shared_ptr<ParamCurve> >
    LRSplineSurface::constParamCurves(double parameter, bool pardir_is_u) const
  //===========================================================================
  {
    MESSAGE("LRSplineSurface::constParamCurves() not implemented yet");
    vector< shared_ptr<ParamCurve> > res;
    return res;
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
}; // end namespace Go


