//===========================================================================
//                                                                           
// File: Identity.h
//                                                                           
// Created: April 2008
//                                                                           
// Author: Vibeke Skytt
//          
// Revision: 
//                                                                           
// Description: Check for identical and embedded entities
//                                                                           
//===========================================================================


#ifndef _IDENTITY_H
#define _IDENTITY_H

#include <memory>
#include "GoTools/utils/Point.h"
#include "GoTools/intersections/GeoTol.h"

#include <vector>

// Coincidence testing

namespace Go 
{

    class ParamSurfaceInt;
    class ParamSurface;
    class ParamCurve;
    
    class Identity
	{
	public:
	    /// Return value = 0 : Not coincident
	    ///              = 1 : Coincident surfaces
	    ///              = 2 : Surface one is embedded in surface two
	    ///              = 3 : Surface two is embedded in surface one
	    int identicalSfs(std::shared_ptr<ParamSurface> sf1,
			     std::shared_ptr<ParamSurface> sf2,
			     double tol);

	    /// Return value = 0 : Not coincident
	    ///              = 1 : Coincident curves
	    ///              = 2 : Curve one is embedded in curve two
	    ///              = 3 : Curve two is embedded in curve one
	    int identicalCvs(std::shared_ptr<ParamCurve> cv1, double start1, double end1,
			     std::shared_ptr<ParamCurve> cv2, double start2, double end2,
			     double tol);
	private:
	    int internalCoincidence(std::shared_ptr<ParamSurfaceInt>& intsf1, 
				    std::shared_ptr<ParamSurfaceInt>& intsf2, 
				    std::shared_ptr<GeoTol>& eps);
	};

} // namespace Go

#endif // _IDENTITY_H

