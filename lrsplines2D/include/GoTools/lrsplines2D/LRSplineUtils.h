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
#include "GoTools/lrsplines2D/BSplineUniLR.h"
#include "GoTools/lrsplines2D/Element2D.h"

namespace Go
{

  /// Utility functionality used in computations involving LR B-spline surfaces
  namespace LRSplineUtils
  {
    enum PointType {
      UNDEF_POINTS = 0,
      REGULAR_POINTS = 1,
      SIGNIFICANT_POINTS = 2,
      GHOST_POINTS = 3
    };

    LRSplineSurface::ElementMap identify_elements_from_mesh(const Mesh2D& m);

    void update_elements_with_single_bspline(LRBSpline2D* b, 
					     LRSplineSurface::ElementMap& emap, 
					     const Mesh2D& mesh,
					     bool remove);

    int locate_interval(const Mesh2D& m, Direction2D d, double value, 
			double other_value, bool at_end);

    void increment_knotvec_indices(std::vector<std::unique_ptr<BSplineUniLR> >& bsplines,
				   const int& from_ix);
    
   LRBSpline2D* 
    insert_basis_function(std::unique_ptr<LRBSpline2D>& b, 
			    const Mesh2D& mesh, 
			    LRSplineSurface::BSplineMap& bmap);
    
    std::vector<int> set_uniform_meshlines(Direction2D d, Mesh2D& mesh);

    bool all_meshlines_uniform(Direction2D d, const Mesh2D& m);

    double compute_greville(const std::vector<int>& v_ixs, 
			    const double* const vals);

    // Calculate the greville values for all b-splines in a knotvector
    // - deg is the polynomial degree of the b-splines
    // - knotvals points to the begining of the list of all possible different knot values
    // - k_vec_in holds the entire knot index vector, thus the i'th B-spline has knot vector indices k_vec_in[i]...k_vec_in[i+deg+1]
    // - Returns the greville value of the b-splines, length is k_vec_in.size() - deg - 1
    std::vector<double> compute_greville(int deg, const std::vector<int>& k_vec_in, const double* const knotvals);


    std::vector<int> knots_to_insert(const std::vector<int>& ref, 
				     const std::vector<int>& mults);

    double compute_alpha(int degree, 
			 const int* const oldvec_ix, 
			 const int* const newvec_ix, 
			 const double* const kvals);

    std::tuple<std::vector<double>, std::vector<int>> 
      insert_knots(const std::vector<int>& new_knots,
		   std::unique_ptr<LRBSpline2D>& bfun,
		   const Direction2D d,
		   const double* const kvals);

    void tensor_split(std::unique_ptr<LRBSpline2D>& bfun,
		      const std::vector<int>& x_mults,
		      const std::vector<int>& y_mults,
		      const Mesh2D& tensor_mesh,
		      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
		      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,
		      LRSplineSurface::BSplineMap& bmap);

    void 
      iteratively_split(std::vector<std::unique_ptr<LRBSpline2D> >& bfuns, 
			const Mesh2D& mesh, 
			std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
			std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2);

    void 
      iteratively_split2(std::vector<LRBSpline2D*>& bsplines,
			 const Mesh2D& mesh,
			 LRSplineSurface::BSplineMap& bmap,
			 double domain[], 
			 std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
			 std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,
			 bool support);

    std::tuple<int, int, int, int>
      refine_mesh(Direction2D d, double fixed_val, double start, double end, 
		  int mult, bool absolute,
		  int spline_degree, double knot_tol,
		  Mesh2D& mesh, 
		  std::vector<std::unique_ptr<BSplineUniLR> >& bsplines,
		  bool& refined);

    bool support_equal(const LRBSpline2D* b1, const LRBSpline2D* b2);

    void split_univariate(std::vector<std::unique_ptr<BSplineUniLR> >& bsplines,
			  int& last, int fixed_ix, int mult);

    bool elementOK(const Element2D* elem, const Mesh2D& m);

    // Lifts a one-dimensional LR-spline function to a three-dimensional function by adding the
    // linear functions for the parameter values before the original function.
    // I.e. the function is changed
    //
    // from    (u,v) |-> (f(u,v))
    // to      (u,v) |-> (u, v, f(u,v))
    void insertParameterFunctions(LRSplineSurface* lr_spline_sf);

    SplineSurface* fullTensorProductSurface(const LRSplineSurface& lr_spline_sf);

    // Fetch the B-spline with the coefficient closest to pos
    // To get a seed parameter for closest point iterations
    LRBSpline2D* mostComparableBspline(LRSplineSurface* lr_spline_sf,
				       Point pos);

    std::vector<std::vector<double> > elementLineClouds(const LRSplineSurface& lr_spline_sf);

    // Distribute given data points to elements
    void distributeDataPoints(LRSplineSurface* srf, std::vector<double>& points, 
			      bool add_distance_field = false, 
			      PointType type = REGULAR_POINTS,
			      bool outlier_flag = false);

    void evalAllBSplines(const std::vector<LRBSpline2D*>& bsplines,
			 double upar, double vpar, 
			 bool u_at_end, bool v_at_end, 
			 std::vector<double>& result);

    void evalAllBSplinePos(const std::vector<LRBSpline2D*>& bsplines,
			   double upar, double vpar, 
			   bool u_at_end, bool v_at_end, 
			   std::vector<Point>& result);

    void
      get_affected_bsplines(const std::vector<LRSplineSurface::Refinement2D>& refs, 
			    const LRSplineSurface::ElementMap& emap,
			    double knot_tol, const Mesh2D& mesh,
			    std::vector<LRBSpline2D*>& affected);
    
    //==============================================================================
    struct support_compare
    //==============================================================================
    {
      bool operator()(const LRBSpline2D* b1, const LRBSpline2D* b2) const  {
	// to compare b1 and b2, compare the x-knotvectors.  If these are identical, compare
	// the y-knotvectors instead.
#if 0//ndef NDEBUG
	std::cout << "kvec_u_size1: " << b1->kvec(XFIXED).size() << std::endl;
	std::cout << "kvec_v_size1: " << b1->kvec(YFIXED).size() << std::endl;
	std::cout << "kvec_u_size2: " << b2->kvec(XFIXED).size() << std::endl;
	std::cout << "kvec_v_size2: " << b2->kvec(YFIXED).size() << std::endl;
#endif
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
