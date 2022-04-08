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

#ifndef LRBSPLINE2DUTILS_H
#define LRBSPLINE2DUTILS_H

#include <vector>
#include <set>
#include "GoTools/utils/checks.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"


namespace Go
{

  /// LRBSpline2D related functionality used in refinement
  
namespace LRBSpline2DUtils
{
  // Derive the knotvector resulting from moving along the u-direction (if d == XFIXED) or 
  // v-direction (if d == YFIXED) from 'beg' to 'end', and  with multiplicities equal to 
  // the 'nu' value of the segment in the ortogonal direction, starting at 'orto_min' and 
  // extending to 'orto_max'.
  std::vector<int> 
    derive_knots(const Mesh2D& m, Direction2D d, int beg, int end, int orto_min, int orto_max);


  // Splits a one-dimensional B-spline by inserting several knots, returns the new knot vector and the
  // coefficients used to express the original B-spline by the new B-splines. The knots are expressed by
  // indices pointing to a common vector of knot values
  // - knotvals points to the begining of the list of all possible different knot values
  // - k_vec_in holds the indices of the knot vector in the original b-spline
  // - new_knots holds all the indices of all the knots to be inserted. They do not have to be different from each other,
  //   or different from the knots in k_vec_in, but they must lie strictly inside the span of k_vec_in
  // - k_vec_out will end up holding the entire knot index vector after insertion, thus the i'th B-spline afterwards
  //   will have vector indices k_vec_out[i]...k_vec_out[i+degree+1].
  //   k_vec_out.size() will be k_vec_in.size() + new_knots.size()
  // - b_spline_weigths will end up holding the B-spline coefficients, thus the original B-spline defined by k_vec_in is
  //   expressed as the sum of b_spline_weigths[i] * (B-spline i)
  //   b_spline_weigths.size() will be 1 + new_knots.size()
  void split_several(const double* knotvals, const std::vector<int>& k_vec_in, const std::vector<int>& new_knots,
		     std::vector<int>& k_vec_out, std::vector<double>& b_spline_weigths);


  // Generates the two LRBSpline2Ds ('new_1', 'new_2') that results
  // from splitting the LRBSpline2D 'orig' in direction 'd', along
  // the knotvalue referenced by 'new_knot_ix'.  Since the LRBSpline2Ds
  // only contain indices to knots, not the knotvals themselves, it is 
  // necessary to supply the knotvals of the concerned direction explicitly 
  // by the two arrays pointed to by 'knotvalues'
  void split_function(const LRBSpline2D& orig, 
		      Direction2D d, 
		      const double* const knotvalues,
		      int new_knot_ix, 
		      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
		      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,
		      LRBSpline2D*& new_1,
		      LRBSpline2D*& new_2);


// if 'b' can be split at least once in the mesh 'm', split it once, and return the 
// result through 'b1' and 'b2'.  The function never carries out more than one split, 
// even when several splits are possible.
// Memory must be handled on the outside of the function.
bool try_split_once(const LRBSpline2D& b, const Mesh2D& mesh, 
		    int mult1, int mult2,
		    std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
		    std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,
		    LRBSpline2D*& b1, 
		    LRBSpline2D*& b2);


#if 0
// Comparison functor which only takes into account the support of the LRBSpline2Ds, 
// not the coefficient or gamma.  This functor is sometimes useful when generating
// STL structures that keep track of LRBSpline2Ds with unique support.
//==============================================================================
struct support_compare
//==============================================================================
{
    bool operator()(const std::unique_ptr<LRBSpline2D>& b1, const std::unique_ptr<LRBSpline2D>& b2) const
    {
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
#endif

}; // End namespace LRBSpline2DUtils
}; // End namespace Go
#endif
