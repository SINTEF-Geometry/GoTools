#ifndef LR_SPLINEUTILS_H
#define LR_SPLINEUTILS_H

#include <array>
#include <functional>
#include <set>
#include <map>
#include <unordered_map>
#include <tuple>
#include <iostream> // @@ debug

#include "GoTools/utils/checks.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/Element2D.h"

namespace Go
{

  namespace LRSplineUtils
  {
    LRSplineSurface::ElementMap identify_elements_from_mesh(const Mesh2D& m);

    void update_elements_with_single_bspline(LRBSpline2D* b, 
					     LRSplineSurface::ElementMap& emap, 
					     const Mesh2D& mesh,
					     bool remove);

    int locate_interval(const Mesh2D& m, Direction2D d, double value, 
			double other_value, bool at_end);

    void increment_knotvec_indices(LRSplineSurface::BSplineMap& bmap, 
				   Direction2D d, int from_ix);

    LRBSpline2D* 
      insert_basis_function(shared_ptr<LRBSpline2D> b, 
			    const Mesh2D& mesh, 
			    LRSplineSurface::BSplineMap& bmap);
    
    std::vector<int> set_uniform_meshlines(Direction2D d, Mesh2D& mesh);

    bool all_meshlines_uniform(Direction2D d, const Mesh2D& m);

    double compute_greville(const std::vector<int>& v_ixs, 
			    const double* const vals);

    std::vector<int> knots_to_insert(const std::vector<int>& ref, 
				     const std::vector<int>& mults);

    double compute_alpha(int degree, 
			 const int* const oldvec_ix, 
			 const int* const newvec_ix, 
			 const double* const kvals);

    std::tuple<std::vector<double>, std::vector<int>> 
      insert_knots(const std::vector<int>& new_knots,
		   shared_ptr<LRBSpline2D> bfun,
		   const Direction2D d,
		   const double* const kvals);

    void tensor_split(shared_ptr<LRBSpline2D> bfun, 
		      const std::vector<int>& x_mults,
		      const std::vector<int>& y_mults,
		      const Mesh2D& tensor_mesh,
		      LRSplineSurface::BSplineMap& bmap);

    void iteratively_split (std::vector<shared_ptr<LRBSpline2D> >& bfuns, 
			    const Mesh2D& mesh);

    void iteratively_split2 (std::vector<LRBSpline2D*>& bsplines,
			     const Mesh2D& mesh,
			     LRSplineSurface::BSplineMap& bmap);

    std::tuple<int, int, int, int>
      refine_mesh(Direction2D d, double fixed_val, double start, double end, 
		  int mult, bool absolute,
		  int spline_degree, double knot_tol,
		  Mesh2D& mesh, LRSplineSurface::BSplineMap& bmap);

    bool support_equal(const LRBSpline2D* b1, const LRBSpline2D* b2);

    bool elementOK(const Element2D* elem, const Mesh2D& m);

    SplineSurface* fullTensorProductSurface(const LRSplineSurface& lr_spline_sf);

    //==============================================================================
    struct support_compare
    //==============================================================================
    {
      bool operator()(const LRBSpline2D* b1, const LRBSpline2D* b2) const  {
	// to compare b1 and b2, compare the x-knotvectors.  If these are identical, compare
	// the y-knotvectors instead.
	const int tmp1 = compare_seq(b1->kvec(XFIXED).begin(), b1->kvec(XFIXED).end(),
				     b2->kvec(XFIXED).begin(), b2->kvec(XFIXED).end());
	if (tmp1 != 0) return (tmp1 < 0);
	const int tmp2 = compare_seq(b1->kvec(YFIXED).begin(), b1->kvec(YFIXED).end(),
				     b2->kvec(YFIXED).begin(), b2->kvec(YFIXED).end());
	return (tmp2 < 0);
      }
    };
      //------------------------------------------------------------------------------

  }; // end namespace LRSplineUtils

}; // end namespace Go

#endif
