//===========================================================================
//                                                                           
// File: GordonSurface.C                                                   
//                                                                           
// Created: Wed Jul 18 18:24:51 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: GordonSurface.C,v 1.6 2009-05-13 07:30:32 vsk Exp $
//                                                                           
// Description: 
//                                                                           
//===========================================================================

#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include <algorithm>
#include <cmath>
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include <fstream> // For debugging.
#include <iterator>

#ifdef __BORLANDC__
#include <iterator>
#endif

using namespace Go;
using std::vector;
using std::max;
using std::min;
using std::swap;
using std::back_inserter;

namespace
{

    /// Given vector of mesh curves, index of a boundary curve is returned.
    /// Returns also index of first curve in mesh_curves in the
    /// 'perpendicular' direction, and parameters of the intersection point.
    void findBndCurve(const std::vector<shared_ptr<SplineCurve> >& mesh_curves,
		      int& index_curr_crv, int& index_end_crv,
		      double& par_curr_crv, double& par_end_crv,
		      double epsgeo);

} // end anonymous namespace

namespace Go
{

// We demand a regular grid in the parameter domain, i.e. we require a u- (v-)
// curve to cross v- (u-) curves at same parameter value in the v- (u-) curve.
//  -------------  Depicted situation is a worst case scenario for the intersections
//     |    | |    of the input curves (we prefer all bndcurves to be given).
//     |    | |    splitMeshCurves() traverses vector of curves until it for a
//     |    | |    curr_crv bnd point finds a bnd-curve (which in the picture
//  ---+----+-+--  amounts to for a vertical curve finding bottom or top curve).
//     |    | |    Bnd curve is put into curves[0], while curr_crv is put into
//  -------------  curves[nmb_crvs - 1]. Respective parameters are stored in
//                 params. With curves[0] as reference, the rest of the
// curves are placed; in 2nd section in case their endpoint crosses curves[0]
// (we then also store their correspondingly parameter value), otherwise in
// first section. When done, we run through 1st section, storing intersections
// in params wrt. curves[nmb_crvs - 1]. In the case of missing boundary curves
// (in max one direction!) these are added. Our curves are now ordered and
// parameterized as curves on a surface; we then make a Gordon surface.
//===========================================================================
SplineSurface*
CoonsPatchGen::createGordonSurface(vector<shared_ptr<SplineCurve> >&
				     curves,
				     vector<double>& params, int& nmb_u_crvs,
				     bool use_param_values)
//===========================================================================
{
    vector<shared_ptr<SplineCurve> > cross_curves;
    vector<int> cross_index;
    return createGordonSurface(curves, params, nmb_u_crvs, cross_curves,
			       cross_index, use_param_values);
}


// All cross_curves are assumed to exist (no null pointers).
//===========================================================================
SplineSurface*
CoonsPatchGen::createGordonSurface(vector<shared_ptr<SplineCurve> >& mesh_curves,
				     vector<double>& params, int& nmb_u_crvs,
				     std::vector<shared_ptr<SplineCurve> >&
				     cross_curves,
				     std::vector<int>& cross_index,
				     bool use_param_values)
//===========================================================================
{
    // Checking that dimensions are the same, and that no curves are rational.
    // Also making sure the curves are all k-regular.
    int dim = mesh_curves[0]->dimension();
    for (size_t i = 0; i < mesh_curves.size(); ++i){
	ALWAYS_ERROR_IF(mesh_curves[i]->dimension() != dim, "Dimension mismatch.");
	ALWAYS_ERROR_IF(mesh_curves[i]->rational(), "Rational curves not supported.");
	mesh_curves[i]->makeKnotStartRegular();
	mesh_curves[i]->makeKnotEndRegular();
    }
    for (size_t i = 0; i < cross_curves.size(); ++i){
	ALWAYS_ERROR_IF(cross_curves[i]->dimension() != dim,
			"Dimension mismatch.");
	ALWAYS_ERROR_IF(cross_curves[i]->rational(), 
		    "Rational curves not supported.");
	cross_curves[i]->makeKnotStartRegular();
	cross_curves[i]->makeKnotEndRegular();
    }

    // In case of bad input data splitMeshCurves() may throw an exception.
    // splitMeshCurves() sets values in params and nmb_u_crvs.
    // @@ In the case of nonempty cross_index, it is currently not altered.
    double epsgeo = 1e-05;
    if (!use_param_values) {
      // We must set up parameters of the intersections.
      try {
	MESSAGE("It seems the index of cross cvs are not updated.");
	splitMeshCurves(mesh_curves, params, nmb_u_crvs, cross_index, epsgeo);
      } catch (...) {
	// @@ sbr To be verified: exception occurs when given only u-curves.
	makeLoftParams(mesh_curves.begin(), (int)mesh_curves.size(), 1.0, params);
	nmb_u_crvs = (int)mesh_curves.size();
	if (cross_curves.size() == 0)
	    return loftSurface(mesh_curves.begin(), (int)mesh_curves.size());
	else THROW("Not implemented yet!");
      }
      sortMeshCurves(mesh_curves, params, nmb_u_crvs, cross_index);
    }

    int nmb_u_cross = 0;
    // Number of cross curves with u parameter direction.
    for (size_t i = 0; i < cross_index.size(); ++i){
      if (cross_index[i] < nmb_u_crvs) {
	++nmb_u_cross;
      }
    }

    // We may need to rearrange the ordering of the mesh_curves. Index should be rising.
    for (size_t i = 0; i < cross_index.size(); ++i) {
      vector<int>::iterator iter = min_element(cross_index.begin() + i, cross_index.end());
      swap(cross_index[i], cross_index[iter - cross_index.begin()]);
      swap(cross_curves[i], cross_curves[iter - cross_index.begin()]);
    }

    // We make sure that all u- (v-) curves have the same order.
    // Also verifying that the u- (v-) curves have the same parameter domain.
    int max_u_order = mesh_curves[0]->order();
    double parmin = mesh_curves[0]->startparam();
    double parmax = mesh_curves[0]->endparam();
    int ci = 0;
    for (int i = 0; i < nmb_u_crvs; ++i) {
	int u_curve_order = mesh_curves[i]->order();
	max_u_order = std::max(max_u_order, u_curve_order);
	if ((int(cross_index.size()) > ci) && (cross_index[ci] == i)) {
	    max_u_order = std::max(cross_curves[ci]->order(), max_u_order);
	    ++ci;
	}
	ALWAYS_ERROR_IF((mesh_curves[i]->startparam() != parmin) ||
			(mesh_curves[i]->endparam() != parmax),
			"u-curves defined over different parameter domains!");
    }
    for (int i = 0; i < nmb_u_cross; ++i)
	ALWAYS_ERROR_IF((cross_curves[i]->startparam() != parmin) ||
			(cross_curves[i]->endparam() != parmax),
			"u-curves defined over different parameter domains!");
    int max_v_order = (nmb_u_crvs < int(mesh_curves.size())) ? 
	mesh_curves[nmb_u_crvs]->order() : 0;
    for (size_t i = nmb_u_crvs; i < mesh_curves.size(); ++i) {
	parmin = mesh_curves[nmb_u_crvs]->startparam();
	parmax = mesh_curves[nmb_u_crvs]->endparam();
	int v_curve_order = mesh_curves[i]->order();
	max_v_order = std::max(max_v_order, v_curve_order);
	if ((int(cross_index.size()) > ci) && (cross_index[ci] == int(i))) {
	    max_v_order = std::max(cross_curves[ci]->order(), max_v_order);
	    ++ci;
	}
	ALWAYS_ERROR_IF((mesh_curves[i]->startparam() != parmin) ||
			(mesh_curves[i]->endparam() != parmax),
			"v-curves defined over different parameter domains!");
    }
    for (size_t i = nmb_u_cross; i < cross_curves.size(); ++i)
	ALWAYS_ERROR_IF((cross_curves[i]->startparam() != parmin) ||
			(cross_curves[i]->endparam() != parmax),
			"v-curves defined over different parameter domains!");

    // We next raise the order where needed.
    for (int i = 0; i < nmb_u_crvs; ++i)
	mesh_curves[i]->raiseOrder(max_u_order - mesh_curves[i]->order());
    for (size_t i = nmb_u_crvs; i < mesh_curves.size(); ++i)
	mesh_curves[i]->raiseOrder(max_v_order - mesh_curves[i]->order());
    for (int i = 0; i < nmb_u_cross; ++i)
	cross_curves[i]->raiseOrder(max_u_order - cross_curves[i]->order());
    for (size_t i = nmb_u_cross; i < cross_curves.size(); ++i)
	cross_curves[i]->raiseOrder(max_v_order - cross_curves[i]->order());

    return doCreateSurface(mesh_curves, params, nmb_u_crvs,
			   cross_curves, cross_index);

}

// There should be no need to include a cross_curve where a mesh_curve is missing.
// If so, then these input values (with indexing w.r.t. mesh_curves) need to change.
//===========================================================================
SplineSurface*
CoonsPatchGen::doCreateSurface(vector<shared_ptr<SplineCurve> >&
				     mesh_curves,
				     vector<double>& params, int& nmb_u_crvs,
				     std::vector<shared_ptr<SplineCurve> >&
				     cross_curves,
				     std::vector<int>& cross_index)
//===========================================================================
{
    const double knot_diff_tol = 1e-05;

    vector<int>::const_iterator iter = cross_index.begin();
    while (iter != cross_index.end())
	if (*iter < nmb_u_crvs) ++iter;
	else break;
    // We gather information about cross curves.
    int nmb_u_cross = (int)(iter - cross_index.begin());

    // Checking whether we have got enough boundary curves.
    // Note that u_min (and its likes) is checked against v-curves.
    // @@@ At some level bnd-params should be rounded off, if close to endpoint.
    double umin = mesh_curves[0]->startparam();
    double umax = mesh_curves[0]->endparam();
    double vmin = mesh_curves[nmb_u_crvs]->startparam();
    double vmax = mesh_curves[nmb_u_crvs]->endparam();
    int missing_u, missing_v; // 0=0_missing, 1=1_missing, 2=2_missing
    // @@ sbr At this point the bnd-parameters should be 'fuzzed' into position.
    // Most likely they're given by the user, and therefore correct. But they
    // may have been created by the split and sort routines.
    //   double snap_tol = 1e-06;  // to be implemented @@ sbr
    if (fabs((params[nmb_u_crvs-1]-params[0]) - (vmax-vmin)) < knot_diff_tol)
		missing_u=0;
    else if ((fabs(params[nmb_u_crvs-1]- vmax) > knot_diff_tol) &&
		(fabs(params[0] - vmin) >knot_diff_tol)) missing_u=2;
    else missing_u=1;
    if ((fabs(params[mesh_curves.size()-1]- umax) < knot_diff_tol) &&
			(fabs(params[nmb_u_crvs] - umin) < knot_diff_tol))
		missing_v=0;	
    else if ((fabs(params[mesh_curves.size()-1] - umax) > knot_diff_tol) &&
	     (fabs(params[nmb_u_crvs] - umin) > knot_diff_tol))
		 missing_v=2;
    else missing_v=1;
    ALWAYS_ERROR_IF((missing_u!=0)&&(missing_v!=0),
		    "Missing boundary curves in both directions!");

    // As we are lofting the u_curves first (in v-direction), we must make
    // sure these boundary curves are given. Otherwise we swap u- with v-curves.
    if (missing_u!=0) {
	vector<shared_ptr<SplineCurve> > curves;
	std::copy(mesh_curves.begin() + nmb_u_crvs, mesh_curves.end(),
		  std::back_inserter(curves));
	std::copy(mesh_curves.begin(), mesh_curves.begin() + nmb_u_crvs,
		  std::back_inserter(curves));
	mesh_curves = curves;

	vector<double> params2;
	std::copy(params.begin() + nmb_u_crvs, params.end(),
		  std::back_inserter(params2));
	std::copy(params.begin(), params.begin() + nmb_u_crvs,
		  std::back_inserter(params2));
	params = params2;

	if (cross_curves.size() != 0) {
	    vector<shared_ptr<SplineCurve> > cross_curves2;
	    std::copy(cross_curves.begin() + nmb_u_cross, cross_curves.end(),
		      std::back_inserter(cross_curves2));
	    std::copy(cross_curves.begin(), cross_curves.begin() + nmb_u_cross,
		      std::back_inserter(cross_curves2));
	    cross_curves = cross_curves2;
	}

	nmb_u_crvs = (int)mesh_curves.size() - nmb_u_crvs;
	std::swap(umin, vmin);
	std::swap(umax, vmax);
	std::swap(missing_u, missing_v);
    }

    int nmb_v_crvs = (int)mesh_curves.size() - nmb_u_crvs;
    // Finally we must make sure that corr. curves are on the same knot vector.
    std::vector<shared_ptr<SplineCurve> > dummy_vector_u, dummy_vector_v;;
    dummy_vector_u.insert(dummy_vector_u.end(),
			  mesh_curves.begin(), mesh_curves.begin() + nmb_u_crvs);
    dummy_vector_u.insert(dummy_vector_u.end(),
			  cross_curves.begin(),
			  cross_curves.begin() + nmb_u_cross);
    GeometryTools::unifyCurveSplineSpace(dummy_vector_u, knot_diff_tol);

    dummy_vector_v.insert(dummy_vector_v.end(),
			  mesh_curves.begin() + nmb_u_crvs, mesh_curves.end());
    dummy_vector_v.insert(dummy_vector_v.end(),
			  cross_curves.begin() + nmb_u_cross, cross_curves.end());
    GeometryTools::unifyCurveSplineSpace(dummy_vector_v, knot_diff_tol);
    // As objects may have changed, we must extract the curves.
    mesh_curves.clear();
    cross_curves.clear();
    mesh_curves.insert(mesh_curves.end(), dummy_vector_u.begin(),
		       dummy_vector_u.begin() + nmb_u_crvs);
    cross_curves.insert(cross_curves.end(),
			dummy_vector_u.begin() + nmb_u_crvs, dummy_vector_u.end());
    mesh_curves.insert(mesh_curves.end(), dummy_vector_v.begin(),
		       dummy_vector_v.begin() + nmb_v_crvs);
    cross_curves.insert(cross_curves.end(),
			dummy_vector_v.begin() + nmb_v_crvs, dummy_vector_v.end());

    // We're now ready to construct our three surfaces.

    // @@ sbr This ought to be replaced when loft takes iterators.
    std::vector<shared_ptr<SplineCurve> >
    u_cross_curves(cross_curves.begin(), cross_curves.begin() + nmb_u_cross);
    std::vector<int> u_cross_index(cross_index.begin(),
				   cross_index.begin() + nmb_u_cross);
    // Lofting in v-direction. loft_u_sf indicates we're lofting u-curves.
    shared_ptr<SplineSurface> loft_u_sf =
      shared_ptr<SplineSurface>(loftSurface(mesh_curves.begin(),
					      params.begin(), nmb_u_crvs,
					      u_cross_curves.begin(),
					      u_cross_index));

    // Checking for missing boundary curves (these should be v-curves).
    // If we are adding missing curve, we must almost certainly insert knots.
    if ((nmb_u_crvs < int(params.size())) && 
			(fabs(params[nmb_u_crvs] - umin) > knot_diff_tol)) {
        // We must add missing umin-curve.
        SplineCurve* new_curve = loft_u_sf->edgeCurve(3);
	vector<double> all_knots, new_mesh_knots, new_loft_knots;
	std::set_union(new_curve->basis().begin(), new_curve->basis().end(),
		   mesh_curves[nmb_u_crvs]->basis().begin(),
		   mesh_curves[nmb_u_crvs]->basis().end(),
		   std::back_inserter(all_knots));
	std::set_difference(all_knots.begin(), all_knots.end(),
			    mesh_curves[nmb_u_crvs]->basis().begin(),
			    mesh_curves[nmb_u_crvs]->basis().end(),
			    std::back_inserter(new_mesh_knots));
	std::set_difference(all_knots.begin(), all_knots.end(),
			    new_curve->basis().begin(), new_curve->basis().end(),
			    std::back_inserter(new_mesh_knots));
	for (size_t i = nmb_u_crvs; i < mesh_curves.size(); ++i)
	    mesh_curves[i]->insertKnot(new_mesh_knots);
	new_curve->insertKnot(new_loft_knots);
	mesh_curves.insert(mesh_curves.begin() + nmb_u_crvs,
			   shared_ptr<SplineCurve>(new_curve));
	params.insert(params.begin() + nmb_u_crvs, umin); // Invalidating missing_u
    }
    if (fabs(params[params.size() - 1] - umax) > knot_diff_tol) {
        // We must add missing umax-curve.
	SplineCurve* new_curve = loft_u_sf->edgeCurve(1);
	vector<double> all_knots, new_mesh_knots, new_loft_knots;
	std::set_union(new_curve->basis().begin(), new_curve->basis().end(),
		       mesh_curves[nmb_u_crvs]->basis().begin(),
		       mesh_curves[nmb_u_crvs]->basis().end(),
		       std::back_inserter(all_knots));
	std::set_difference(all_knots.begin(), all_knots.end(),
			    mesh_curves[nmb_u_crvs]->basis().begin(),
			    mesh_curves[nmb_u_crvs]->basis().end(),
			    std::back_inserter(new_mesh_knots));
	std::set_difference(all_knots.begin(), all_knots.end(),
			    new_curve->basis().begin(), new_curve->basis().end(),
			    std::back_inserter(new_loft_knots));
	for (size_t i = nmb_u_crvs; i < mesh_curves.size(); ++i)
	    mesh_curves[i]->insertKnot(new_mesh_knots);
	new_curve->insertKnot(new_loft_knots);
	mesh_curves.insert(mesh_curves.end(),
			   shared_ptr<SplineCurve>(new_curve));
	params.insert(params.end(), umax); // Invalidating missing_u
    }

    std::vector<shared_ptr<SplineCurve> >
    v_cross_curves(cross_curves.begin() + nmb_u_cross, cross_curves.end());
    std::vector<int> v_cross_index(cross_index.begin() + nmb_u_cross,
				   cross_index.end());
    for (size_t i = 0; i < v_cross_index.size(); ++i)
	v_cross_index[i] -= nmb_u_crvs;
    // Creating the lofted surface in the u-direction, we must bear in mind
    // that the surface created needs to swap parameter directions.
    shared_ptr<SplineSurface> loft_v_sf =
      shared_ptr<SplineSurface>(loftSurface(mesh_curves.begin() + nmb_u_crvs,
					      params.begin() + nmb_u_crvs,
					    (int)mesh_curves.size() - nmb_u_crvs,
					      v_cross_curves.begin(),
					      v_cross_index));
    loft_v_sf->swapParameterDirection();

    // We then make the tp-surface.
    shared_ptr<SplineSurface> tp_sf =
      shared_ptr<SplineSurface>(tpSurface(mesh_curves, params, nmb_u_crvs,
					    cross_curves, cross_index));

    // Put the surfaces on a joint basis.
    std::vector<shared_ptr<SplineSurface> > surfaces;
    surfaces.push_back(loft_u_sf);
    surfaces.push_back(loft_v_sf);
    surfaces.push_back(tp_sf);
    GeometryTools::unifySurfaceSplineSpace(surfaces, knot_diff_tol);

    // As our push_backed elements may have changed, we use vector elements.
    // final_surface = surfaces[0] + surfaces[1] - surfaces[2]->
    for (int i = 0; i < surfaces[0]->dimension() *
	     surfaces[0]->numCoefs_u() * surfaces[0]->numCoefs_v(); ++i) {
	surfaces[0]->coefs_begin()[i] += surfaces[1]->coefs_begin()[i];
	surfaces[0]->coefs_begin()[i] -= surfaces[2]->coefs_begin()[i];
    }

    // That should be it. Returning our promised Gordon surface.
    // As surfaces[0] is controlled by a smart pointer, we must make a copy.
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    return dynamic_cast<SplineSurface*>(surfaces[0]->clone());
    
#else
    return surfaces[0]->clone();
#endif
}


//===========================================================================
SplineSurface*
CoonsPatchGen::loftSurface(std::vector<shared_ptr<SplineCurve> >::iterator
				 first_curve,
				 int nmb_crvs)
//===========================================================================
{
    vector<double> params;
    // @@ We're giving the surface a parameter domain of length 1.0 in v-direct.
    // It would be more natural to set value in accordance with length of
    // u-direction and geometric lengths, possibly take an input value.
    makeLoftParams(first_curve, nmb_crvs, 1.0, params);

    vector<shared_ptr<SplineCurve> > dummy_cross_curves;
    vector<int> dummy_index;
    return loftSurface(first_curve, params.begin(), nmb_crvs,
		       dummy_cross_curves.begin(), dummy_index);
}


//===========================================================================
SplineSurface*
CoonsPatchGen::loftSurface(std::vector<shared_ptr<SplineCurve> >::iterator
				 first_curve,
				 std::vector<double>::iterator first_param,
				 int nmb_crvs)
//===========================================================================
{
    vector<shared_ptr<SplineCurve> > dummy_cross_curves;
    vector<int> dummy_index;

    return loftSurface(first_curve, first_param, nmb_crvs,
		       dummy_cross_curves.begin(), dummy_index);

}

// bnd_cross_curves assumed to be boundary information.
// The general routine to be implemented asap.
// @@ sbr Order may be given as input; same applies for tpSurface (for consistency).
//===========================================================================
SplineSurface*
CoonsPatchGen::loftSurface(std::vector<shared_ptr<SplineCurve> >::iterator
			     first_curve,
			     std::vector<double>::iterator first_param,
			     int nmb_crvs,
			     vector<shared_ptr<SplineCurve> >::iterator
			     first_cross_curve,
			     vector<int>& cross_index)
//===========================================================================
{

    int nmb_cross_crvs = (int)cross_index.size();

  bool rational = false;
  for (int i = 0; i < nmb_crvs; ++i)
    if (first_curve[i]->rational())
      rational = true;
  if (nmb_cross_crvs != 0)
    for (int i = 0; i < nmb_crvs; ++i)
      if (first_cross_curve[i]->rational())
	rational = true;

  vector<shared_ptr<SplineCurve> >::iterator first_curve_copy;
  vector<shared_ptr<SplineCurve> >::iterator first_cross_curve_copy;
  vector<shared_ptr<SplineCurve> > curves_copy;
  vector<shared_ptr<SplineCurve> > cross_curves_copy;
  if (!rational)
    {
      first_curve_copy = first_curve;
      first_cross_curve_copy = first_cross_curve;
    }
  else
    {
      for (int i = 0; i <nmb_crvs; ++i)
	{
	  if (first_curve[i]->rational())
	    curves_copy.push_back(first_curve[i]);
	  else
	    {
	      SplineCurve* sc = first_curve[i]->clone();
	      sc->representAsRational();
	      curves_copy.push_back(shared_ptr<SplineCurve>(sc));
	    }
	}
      first_curve_copy = curves_copy.begin();
      if (nmb_cross_crvs != 0)
	{
	  for (int i = 0; i <nmb_crvs; ++i)
	    {
	      if (first_cross_curve[i]->rational())
		cross_curves_copy.push_back(first_cross_curve[i]);
	      else
		{
		  SplineCurve* sc = first_cross_curve[i]->clone();
		  sc->representAsRational();
		  cross_curves_copy.push_back(shared_ptr<SplineCurve>(sc));
		}
	    }
	  first_cross_curve_copy = cross_curves_copy.begin();
	}
      else
	first_cross_curve_copy = first_cross_curve;
    }

//     MESSAGE_IF(rational,
// 		  "Rational curves not tested yet!");

  // This may not be true, but no testing has been performed on general case.
  MESSAGE_IF(cross_index.size() > 2,
	     "Only input of boundary cross_curves has been tested!");

  // By looking at each curve as a point of dimension dim*u_coefs.size(),
  // the interpolation method gets rather easy.
  vector<double> coefs, coefs_u_loft, params(nmb_crvs);
  for (int i = 0; i < nmb_crvs; ++i) {
    if (rational)
      coefs.insert(coefs.end(),
		   first_curve_copy[i]->rcoefs_begin(), first_curve_copy[i]->rcoefs_end());
    else
      coefs.insert(coefs.end(),
		   first_curve[i]->coefs_begin(), first_curve[i]->coefs_end());
    params[i] = first_param[i];
  }

  vector<double> tangent_points;
  if (nmb_cross_crvs != 0) {
    // We must extract tangent information from the bnd_cross_curves.
    // All curves are assumed to have equal dim and basis.
    int dim = first_cross_curve[0]->dimension();
    if (rational)
      dim += 1;
    for (int i = 0; i < nmb_cross_crvs; ++i)
      if (rational)
	tangent_points.insert(tangent_points.end(),
			      first_cross_curve[i]->rcoefs_begin(),
			      first_cross_curve[i]->rcoefs_end());
      else
	tangent_points.insert(tangent_points.end(),
			      first_cross_curve[i]->coefs_begin(),
			      first_cross_curve[i]->coefs_end());
  }

  int order = std::min(4, nmb_crvs + nmb_cross_crvs);
  SplineInterpolator u_interpolator;
  u_interpolator.makeBasis(params, cross_index, order);
  u_interpolator.interpolate(params, coefs, cross_index,
			     tangent_points, coefs_u_loft);

  // We're lofting in 2nd parameter (v) direction.
  SplineSurface* surf;
  if (rational)
    surf = new SplineSurface(first_curve_copy[0]->basis(), u_interpolator.basis(),
			     coefs_u_loft.begin(), first_curve_copy[0]->dimension(), true);
  else
    surf = new SplineSurface(first_curve[0]->basis(), u_interpolator.basis(),
			     coefs_u_loft.begin(), first_curve[0]->dimension(), false);
  return surf;
}


// cross_curves assumed to be boundary information. We also assume
// u(v)-cross_curves and u(v)-mesh_curves are on the same knot vectors, with the
// same order.
// Furthermore all boundary curves (plus maybe som more) are assumed to be present.
// For both vectors we assume u-curves are given first.
//===========================================================================
SplineSurface*
CoonsPatchGen::tpSurface(const vector<shared_ptr<SplineCurve> >& mesh_curves,
			       vector<double> params, int nmb_u_crvs,
			       const vector<shared_ptr<SplineCurve> >&
			       cross_curves,
			       vector<int>& cross_index)
//===========================================================================
{
    // We currently handle cross curves on boundaries only.
    MESSAGE_IF(cross_curves.size() > 4,
		  "Only input of boundary cross_curves has been tested.");
    ALWAYS_ERROR_IF((nmb_u_crvs == 0) || (mesh_curves.size() - nmb_u_crvs == 0),
		    "tpSurface assumes input of curves in both directions! Use loft! ");


    int nmb_v_crvs = (int)mesh_curves.size() - nmb_u_crvs;

    vector<int>::iterator iter = cross_index.begin();
    while (iter != cross_index.end())
	if (*iter < nmb_u_crvs) ++iter;
	else break;

    // We gather information about cross curves.
    int nmb_u_cross = (int)(iter - cross_index.begin());
    int nmb_v_cross = (int)(cross_index.end() - iter);
    int dim = mesh_curves[0]->dimension(); // Should be the same for all curves.
    // We interpolate in the u-direction first.
    // The first thing we must do is to gather the interpolation points.
    vector<double> points; // p0x, p0y, p0z, p1x, p1y, p1z, p2x,...
    Point pnt;
    for (size_t i = nmb_u_crvs; i < mesh_curves.size(); ++i)
	for (int j = 0; j < nmb_u_crvs; ++j) {
	    mesh_curves[i]->point(pnt, params[j]);
	    for (int k = 0; k < dim; ++k)
		points.push_back(pnt[k]);
	}
    vector<double> tangent_points;
    if (nmb_v_cross != 0) {
	// We must extract tangent information from the bnd_cross_curves.
	// All curves are assumed to have equal dim and basis.
	tangent_points.resize(dim * nmb_u_crvs * nmb_v_cross);
	for (size_t i = nmb_u_cross; i < cross_index.size(); ++i)
	    for (int j = 0; j < nmb_u_crvs; ++j) {
		double fac = params[j];
		Point tangent_pt(dim);
		cross_curves[i]->point(tangent_pt, fac);
		for (int k = 0; k < dim; ++k)
		    tangent_points[j*dim + k + (i-nmb_u_cross)*dim*nmb_u_crvs] =
			tangent_pt[k];
	    }
    }

    vector<double> coefs, new_coefs;
    std::vector<double> u_params(params.begin() + nmb_u_crvs, params.end());
    std::vector<int> v_cross_index(cross_index.begin() + nmb_u_cross,
				   cross_index.end());
    for (size_t i = 0; i < v_cross_index.size(); ++i)
	v_cross_index[i] -= nmb_u_crvs;
    int u_order = std::min(4, nmb_v_crvs + nmb_v_cross);
    SplineInterpolator v_interpolator; // We're interpolating v-curves.
    v_interpolator.makeBasis(u_params, v_cross_index, u_order);
    v_interpolator.interpolate(u_params, points, v_cross_index,
			       tangent_points, coefs);

    // Transpose the coefs.
    // There should be no doubt that the basis is wrt. the u-direction...
    SplineUtils::transpose_array(dim, v_interpolator.basis().numCoefs(),
		    nmb_u_crvs, &coefs[0]);


    tangent_points.clear();
    // We interpolate the cross tangents in order to get them on the same basis.
    // We use the same interpolator as in the first lofting.
    if (nmb_u_cross != 0) {
//  	tangent_points.resize(dim * v_interpolator.basis().numCoefs() *
//  			      nmb_u_cross);
	for (int i = 0; i < nmb_u_cross; ++i) {
	    int ti = nmb_u_cross; // 
	    vector<double> cross_coefs;
	    vector<double> cross_points; //(dim*nmb_v_crvs);
	    vector<double> cross_tangent_points; //(dim*nmb_v_cross);
	    vector<int> cross_ind; // Referring to info about tangents.
	    for (size_t j = nmb_u_crvs; j < mesh_curves.size(); ++j) {
		double fac = params[j];
		vector<Point> tangent_pt(2);
		cross_curves[i]->point(tangent_pt, fac, 1);
		for (int k = 0; k < dim; ++k) {
		    cross_points.push_back(tangent_pt[0][k]);
		}
		// We must check whether tangent information is present.
		if ((ti < int(cross_index.size())) && (cross_index[ti] == int(j))) {
//  		    double fac2 = params[cross_index[ti]];
//  		    vectorPoint cross_tangent_pt(dim);
//  		    cross_curves[i]->tangent(cross_tangent_pt, fac);
		    for (int k = 0; k < dim; ++k)
			cross_tangent_points.push_back(tangent_pt[1][k]);
		    cross_ind.push_back((int)j - nmb_u_crvs);
		    ++ti;
		}
	    }
	    v_interpolator.interpolate(u_params, cross_points, cross_ind,
				       cross_tangent_points, cross_coefs);
//  	    tangent_points.insert(tangent_points.begin(),
//  				  cross_coefs.begin(), cross_coefs.end());
	    tangent_points.insert(tangent_points.end(),
				  cross_coefs.begin(), cross_coefs.end());
	}
    }

    std::vector<double> v_params(params.begin(), params.begin() + nmb_u_crvs);
    std::vector<int> u_cross_index(cross_index.begin(),
				   cross_index.begin() + nmb_u_cross);
    int v_order = std::min(4, nmb_u_crvs + nmb_u_cross);
    SplineInterpolator u_interpolator; // We're interpolating u-curves.
    u_interpolator.makeBasis(v_params, u_cross_index, v_order);
    u_interpolator.interpolate(v_params, coefs, u_cross_index,
			       tangent_points, new_coefs);

    return new SplineSurface(v_interpolator.basis(), u_interpolator.basis(),
			       new_coefs.begin(), dim);
}



// Given vector of mesh_curves, u-curves are moved to the front.
// Respective parameter values are returned in params.
//===========================================================================
void CoonsPatchGen::splitMeshCurves(vector<shared_ptr<SplineCurve> >&
				      mesh_curves,
				      vector<double>& params, int& nmb_u_crvs,
				      vector<int>& cross_index, double epsgeo)
//===========================================================================
{
    vector<int> mesh_ind;
    for (size_t ki = 0; ki < mesh_curves.size(); ++ki) {
        mesh_ind.push_back((int)ki);
    }

    int nmb_crvs = (int)mesh_curves.size();
    params.resize(nmb_crvs);
    int ibc, ipc; // Index boundary curve, index 'perpendicular' (wrt. bc) curve.
    double ppc, pbc; // Param value for perpend. curve, param value for bound curve.
    int nuc = 0; // Number of u_curves.
    int nvc = 0; // Number of v curves.
    try {
	findBndCurve(mesh_curves, ipc, ibc, ppc, pbc, epsgeo);
    }
    catch (...) {
	throw UnKnownError();
    }

    // We now order the information we have got so far.
    // The first ipc+1 elem of mesh_curves are made u-curves, last is v-curve.
    nuc = ipc;
    nvc = 1;
    // Put perpendicular curve last.
    std::swap(mesh_curves[ipc], mesh_curves[nmb_crvs - 1]);
    std::swap(mesh_ind[ipc], mesh_ind[nmb_crvs - 1]);
    if (ibc == nmb_crvs - 1) ibc = ipc; // In case bnd_crv was swapped.
    if (ibc < ipc) nuc = ipc;
    else {
	std::swap(mesh_curves[ipc], mesh_curves[ibc]);
	std::swap(mesh_ind[ipc], mesh_ind[ibc]);
	ibc = ipc;
	nuc = ipc + 1;
    }
    if (ibc != 0) { // Put bnd_curv up front.
	std::swap(mesh_curves[ipc], mesh_curves[ibc]);
	std::swap(mesh_ind[ipc], mesh_ind[ibc]);
	ibc = 0;
    }
    ipc = nmb_crvs - 1;
    params[ibc] = pbc;
    params[ipc] = ppc;
    Point bnd_pt, pt, clo_pt;
    double clo_t, clo_dist;
    shared_ptr<SplineCurve> bnd_crv = mesh_curves[ibc];
    for (int ki = nuc; ki < nmb_crvs - nvc; ++ki) {
	mesh_curves[ki]->point(bnd_pt, pbc);
	bnd_crv->closestPoint(bnd_pt, bnd_crv->startparam(),
			      bnd_crv->endparam(), clo_t, clo_pt, clo_dist);
	if (clo_dist < epsgeo) {
	    ++nvc;
	    std::swap(mesh_curves[ki], mesh_curves[nmb_crvs - nvc]);
	    std::swap(mesh_ind[ki], mesh_ind[nmb_crvs - nvc]);
	    params[nmb_crvs - nvc] = clo_t;
	    --ki;
	} else ++nuc;
    }
    // The curves are now split, all but one param-value in u dir still unknown.
    // Checking against perp_crv, we get their param values.
    double tpar = params[nmb_crvs - 1];
    shared_ptr<SplineCurve> perp_crv = mesh_curves[ipc];
    for (int ki = 1; ki < nuc; ++ki) {
	mesh_curves[ki]->point(pt, tpar);
	perp_crv->closestPoint(pt, perp_crv->startparam(), perp_crv->endparam(),
			       clo_t, clo_pt, clo_dist);
	if (clo_dist > epsgeo) {
		    // This will typically happen if curves do not form a 
		    // grid in the parameter domain.
	    MESSAGE("Suspecting bad parametrization.");
	    throw UnKnownError();
	}
	params[ki] = clo_t;
    }
    nmb_u_crvs = nuc;

    // Finally we must update cross_index vector.
    for (size_t ki = 0; ki < cross_index.size(); ++ki) {
        vector<int>::const_iterator iter = 
	    std::find(mesh_ind.begin(), mesh_ind.end(), cross_index[ki]);
	ASSERT(iter != mesh_ind.end());
	cross_index[ki] = (int)(iter - mesh_ind.begin());
    }
}


//===========================================================================
void CoonsPatchGen::sortMeshCurves(vector<shared_ptr<SplineCurve> >&
				     mesh_curves,
				     vector<double>& params,
				     int nmb_u_crvs,
				     vector<int>& cross_index)
//===========================================================================
{
    vector<double> cp_params = params; // We make a copy of params.
    // Sort values in cp_params, u values first, then v.
    std::sort(cp_params.begin(), cp_params.begin() + nmb_u_crvs);
    std::sort(cp_params.begin() + nmb_u_crvs, cp_params.end());

    vector<shared_ptr<SplineCurve> > cp_mesh_curves;
    vector<int> new_mesh_ind;
    // Sort u curves...
    for (int i = 0; i < nmb_u_crvs; ++i) {
        std::vector<double>::iterator iter =
	    std::find(params.begin(), params.begin() + nmb_u_crvs, cp_params[i]);
	cp_mesh_curves.push_back(mesh_curves[iter - params.begin()]);
	new_mesh_ind.push_back((int)(iter - params.begin()));
    }
    // ..., then v curves.
    for (size_t i = nmb_u_crvs; i < mesh_curves.size(); ++i) {
        std::vector<double>::iterator iter =
	    std::find(params.begin() + nmb_u_crvs, params.end(), cp_params[i]);
	cp_mesh_curves.push_back(mesh_curves[iter - params.begin()]);
	new_mesh_ind.push_back((int)(iter - params.begin()));
    }
    mesh_curves = cp_mesh_curves;
    params = cp_params;

    for (size_t i = 0; i < cross_index.size(); ++i) {
      vector<int>::iterator iter =
	  std::find(new_mesh_ind.begin(), new_mesh_ind.end(), cross_index[i]);
      cross_index[i] = (int)(iter - new_mesh_ind.begin());
    }
}

} // end namespace Go.


namespace
{

// Locate index of a boundary curve intersecting a curve. Get params.
// index_curr_crv is used for an efficient split of curves.
//===========================================================================
void findBndCurve(const vector<shared_ptr<SplineCurve> >&
		  mesh_curves,
		  int& index_curr_crv, int& index_bnd_crv,
		  double& par_curr_crv, double& par_bnd_crv,
		  double epsgeo)
//===========================================================================
{
    shared_ptr<SplineCurve> curr_crv;
    shared_ptr<SplineCurve> other_crv;
    Point startpt, endpt, clo_pt;
    double tpar, clo_t, clo_dist;
    index_curr_crv = -1; // Index of current curve.
    index_bnd_crv = -1; // Denotes that a boundary curve is not yet found.
    while (true) {
	++index_curr_crv;
	if (index_curr_crv == int(mesh_curves.size())) {
	    // MESSAGE("Missing all boundary curves, or only given u-curves.");
	    throw;
	}
	curr_crv = mesh_curves[index_curr_crv];
	for (size_t  i = 0; i < mesh_curves.size(); ++i)
	    if (int(i) != index_curr_crv) {
		other_crv = mesh_curves[i];
		tpar = curr_crv->startparam();
		curr_crv->point(startpt, tpar);
		other_crv->closestPoint(startpt, other_crv->startparam(),
					other_crv->endparam(), clo_t,
					clo_pt, clo_dist);
		if (clo_dist < epsgeo) {
		    index_bnd_crv = (int)i;
		    par_bnd_crv = tpar;
		    par_curr_crv = clo_t;
		    return;
		} else {
		    tpar = curr_crv->endparam();
		    curr_crv->point(endpt, tpar);
		    other_crv->closestPoint(endpt, other_crv->startparam(),
					    other_crv->endparam(), clo_t,
					    clo_pt, clo_dist);
		    if (clo_dist < epsgeo) {
			index_bnd_crv = (int)i;
			par_bnd_crv = tpar;
			par_curr_crv = clo_t;
			return;
		    }
		}
	    }
    }
}

} // end anonymous namespace
