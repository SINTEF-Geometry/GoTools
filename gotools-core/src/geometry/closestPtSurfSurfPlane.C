//===========================================================================
//                                                                           
// File: closestPtSurfSurfPlane.C                                            
//                                                                           
// Created: Thu Apr 28 13:00:13 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: closestPtSurfSurfPlane.C,v 1.8 2005-06-29 13:08:49 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/closestPtSurfSurfPlane.h"

namespace Go {

//===========================================================================
void closestPtSurfSurfPlane(const std::vector<Point>& epoint,
			    const std::vector<Point>& epnt1,
			    const std::vector<Point>& epnt2,
			    const Point& epar1,
			    const Point& epar2,
			    const ParamSurface* psurf1,
			    const ParamSurface* psurf2,
			    double aepsge,
			    std::vector<Point>& gpnt1,
			    std::vector<Point>& gpnt2,
			    Point& gpar1, 
			    Point& gpar2, 
			    int& jstat,
			    AlgorithmChoice algo)
//===========================================================================
{
    if (algo == GEOMETRICAL) {
	closestPtSurfSurfPlaneGeometrical(epoint, 
					  epnt1, 
					  epnt2, 
					  epar1, 
					  epar2, 
					  psurf1, 
					  psurf2, 
					  aepsge, 
					  gpnt1, 
					  gpnt2, 
					  gpar1, 
					  gpar2, 
					  jstat);
    } else {
	// algo == FUNCTIONAL
	closestPtSurfSurfPlaneFunctional(epoint, 
					  epnt1, 
					  epnt2, 
					  epar1, 
					  epar2, 
					  psurf1, 
					  psurf2, 
					  aepsge, 
					  gpnt1, 
					  gpnt2, 
					  gpar1, 
					  gpar2, 
					  jstat);
    }
}

}; // end namespace Go;

