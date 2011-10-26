//===========================================================================
//                                                                           
// File: SingularityType.h                                                   
//                                                                           
// Created: Fri Jan  7 12:33:17 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: SingularityType.h,v 1.2 2005-08-22 08:15:39 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SINGULARITYTYPE_H
#define _SINGULARITYTYPE_H

namespace Go {
    /// A classification of an IntersectionPoint, based on the nature of the
    /// intersection on which it lies.
  enum SingularityType
      {
	  /// IntersectionPoint lies on normal, transversal intersection
	  ORDINARY_POINT = 0, 
	  /// The two intersecting objects are tangent to each other at the
	  /// location of this IntersectionPoint.
	  TANGENTIAL_POINT, 
	  /// This is an isolated IntersectionPoint (does not lie on an 
	  /// IntersectionCurve or on a partial coincidence area).
	  ISOLATED_POINT, 
          /// This is an IntersectionPoint where flere interection curves meet.
	  BRANCH_POINT,  
	  /// The two objects intersecting at this point have same tangent
	  /// plane AND curvature.
	  HIGHER_ORDER_POINT 
      };
};

#endif // _SINGULARITYTYPE_H
