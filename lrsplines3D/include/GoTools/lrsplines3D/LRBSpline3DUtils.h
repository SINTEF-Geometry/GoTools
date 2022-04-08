//===========================================================================
//                                                                           
// File: LRBSpline3DUtils.h                                                  
//                                                                           
// Created: Mon Feb 25 11:08:17 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

