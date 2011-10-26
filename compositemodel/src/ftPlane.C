//===========================================================================
//                                                                           
// File: ftPlane.C                                                           
//                                                                           
// Created: Wed May 31 11:52:52 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ftPlane.C,v 1.3 2008-04-03 13:43:15 kfp Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/utils/BoundingBox.h"


namespace Go
{


  //===========================================================================
  ftPlane::ftPlane(const ftPlane& plane)
    :normal_(plane.normal_), point_(plane.point_)
  //===========================================================================
  {
  }

  // A plane intersects a box if it separates its corners
  //===========================================================================
  bool ftPlane::intersectsBox(const BoundingBox& box) const
  //===========================================================================
  {
    Point corner[2];
    corner[0] = box.low() - point_;
    corner[1] = box.high() - point_;
    double firstNorm = corner[0]*normal_;
    if (firstNorm * (corner[1]*normal_) <= 0) return true;
    for (int i = 1; i < 7; ++i)
      if (firstNorm * (Point(corner[i>>2][0],
			     corner[(i>>1)&1][1],
			     corner[i&1][2])
		       *normal_)
	  <= 0)
	return true;
    return false;
  }

} // namespace Go
