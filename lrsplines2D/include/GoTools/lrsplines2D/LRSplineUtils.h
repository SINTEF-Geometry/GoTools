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
    insert_basis_function(std::unique_ptr<LRBSpline2D>& b, 
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
		   std::unique_ptr<LRBSpline2D>& bfun,
		   const Direction2D d,
		   const double* const kvals);

    void tensor_split(std::unique_ptr<LRBSpline2D>& bfun,
		      const std::vector<int>& x_mults,
		      const std::vector<int>& y_mults,
		      const Mesh2D& tensor_mesh,
		      LRSplineSurface::BSplineMap& bmap);

    void iteratively_split (std::vector<std::unique_ptr<LRBSpline2D> >& bfuns, 
			    const Mesh2D& mesh);

    void iteratively_split2 (std::vector<LRBSpline2D*>& bsplines,
			     const Mesh2D& mesh,
			     LRSplineSurface::BSplineMap& bmap,
			     double domain[]);

    std::tuple<int, int, int, int>
      refine_mesh(Direction2D d, double fixed_val, double start, double end, 
		  int mult, bool absolute,
		  int spline_degree, double knot_tol,
		  Mesh2D& mesh, LRSplineSurface::BSplineMap& bmap);

    bool support_equal(const LRBSpline2D* b1, const LRBSpline2D* b2);

    bool elementOK(const Element2D* elem, const Mesh2D& m);

    SplineSurface* fullTensorProductSurface(const LRSplineSurface& lr_spline_sf);

    // Fetch the B-spline with the coefficient closest to pos
    // To get a seed parameter for closest point iterations
    LRBSpline2D* mostComparableBspline(LRSplineSurface* lr_spline_sf,
				       Point pos);

    std::vector<std::vector<double> > elementLineClouds(const LRSplineSurface& lr_spline_sf);

    // Distribute given data points to elements
    void distributeDataPoints(LRSplineSurface* srf, std::vector<double>& points, 
			      bool add_distance_field = false, 
			      bool primary_points = true);


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
