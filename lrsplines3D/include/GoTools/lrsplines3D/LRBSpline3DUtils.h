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

#ifndef _LRBSPLINE3DUTILS_H
#define _LRBSPLINE3DUTILS_H



#include <vector>
#include <set>
#include "GoTools/utils/checks.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"


namespace Go
{
  /// LRBSpline3D related functionality used in refinement

namespace LRBSpline3DUtils
{


  // Derive the knotvector resulting from moving along the u-direction (if d == XFIXED) or 
  // v-direction (if d == YFIXED) from 'beg' to 'end', and  with multiplicities equal to 
  // the 'nu' value of the segment in the ortogonal direction, starting at 'orto_min' and 
  // extending to 'orto_max'.
  std::vector<int> 
    derive_knots(const Mesh3D& m, Direction3D d, int beg, int end,
		 int orto1_min, int orto1_max,
		 int orto2_min, int orto2_max);

  // Generates the two LRBSpline2Ds ('new_1', 'new_2') that results
  // from splitting the LRBSpline2D 'orig' in direction 'd', along
  // the knotvalue referenced by 'new_knot_ix'.  Since the LRBSpline2Ds
  // only contain indices to knots, not the knotvals themselves, it is 
  // necessary to supply the knotvals of the concerned direction explicitly 
  // by the two arrays pointed to by 'knotvalues'
  void split_function(const LRBSpline3D& orig, 
		      Direction3D d, 
		      const double* const kvals,
		      int new_knot_ix,
		      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
		      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,
		      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec3,
		      LRBSpline3D*& new_1,
		      LRBSpline3D*& new_2);

// if 'b' can be split at least once in the mesh 'm', split it once, and return the 
// result through 'b1' and 'b2'.  The function never carries out more than one split, 
// even when several splits are possible.
bool try_split_once(const LRBSpline3D& b, const Mesh3D& mesh, 
		    int mult1, int mult2, int mult3,
		    std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
		    std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,
		    std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec3,
                    LRBSpline3D*& b1,
                    LRBSpline3D*& b2);






} // End namespace LRBSpline3DUtils
} // End namespace Go

#endif // _LRBSPLINE3DUTILS_H

