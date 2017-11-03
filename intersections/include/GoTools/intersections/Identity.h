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

#ifndef _IDENTITY_H
#define _IDENTITY_H

#include <memory>
#include "GoTools/utils/Point.h"
#include "GoTools/intersections/GeoTol.h"

#include <vector>

namespace Go 
{

    class ParamSurfaceInt;
    class ParamSurface;
    class ParamCurve;
    
    /// Check coincidence

    class Identity
	{
	public:
	    /// Return value = 0 : Not coincident
	    ///              = 1 : Coincident surfaces
	    ///              = 2 : Surface one is embedded in surface two
	    ///              = 3 : Surface two is embedded in surface one
	    int identicalSfs(shared_ptr<ParamSurface> sf1,
			     shared_ptr<ParamSurface> sf2,
			     double tol);

	    /// Return value = 0 : Not coincident
	    ///              = 1 : Coincident curves
	    ///              = 2 : Curve one is embedded in curve two
	    ///              = 3 : Curve two is embedded in curve one
	    int identicalCvs(shared_ptr<ParamCurve> cv1, double start1, double end1,
			     shared_ptr<ParamCurve> cv2, double start2, double end2,
			     double tol);

	    /// Return value = 0 : Not coincident
	    ///              = 1 : Coincident curves
	    ///              = 2 : Curve one is embedded in curve two
	    ///              = 3 : Curve two is embedded in curve one
	    int identicalCvs(shared_ptr<ParamCurve> cv1, 
			     shared_ptr<ParamCurve> cv2, 
			     double tol);
	private:
	    int internalCoincidence(shared_ptr<ParamSurfaceInt>& intsf1, 
				    shared_ptr<ParamSurfaceInt>& intsf2, 
				    shared_ptr<GeoTol>& eps);
	};

} // namespace Go

#endif // _IDENTITY_H

