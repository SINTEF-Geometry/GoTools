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

#include "GoTools/implicitization/BernsteinUtils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Array.h"
#include <iostream>


using namespace Go;
using namespace std;


int main()
{
    cout << "*** BernsteinUtils ***" << endl;
    cout << endl;

    cout << "*** splineToBernstein - curves ***" << endl;
    cout << endl;

    // Construct a spline curve
    int nctrlpoints = 3;
    int order = 3;
    double knots[] = { 0.0, 0.0, 0.0, 2.0, 2.0, 2.0};
    double ccoefs[]= { -1.0, 1.0, 1.0, 0.0, -1.0, 2.0, 1.0, 1.0, 1.0 };
    const int cdim = 2;
    bool rational = true;
    SplineCurve curve(nctrlpoints, order, knots, ccoefs, cdim, rational);
    cout << curve << endl;

    // Convert to Bernstein
    Array<BernsteinPoly, cdim+1> curve_bp;
    splineToBernstein(curve, curve_bp);
    cout << curve_bp[0] << endl 
	 << curve_bp[1] << endl
	 << curve_bp[2] << endl << endl;

    cout << "*** bernsteinToSpline - curves ***" << endl;
    cout << endl;

    // Convert again to SplineCurve
    SplineCurve curve2;
    bernsteinToSpline(curve_bp, rational, curve2);
    cout << curve2 << endl;

    cout << "*** splineToBernstein - surfaces ***" << endl;
    cout << endl;

    // Construct a spline surface
    int nctrlpointsu = 3;
    int nctrlpointsv = 3;
    int orderu = 3;
    int orderv = 3;
    double knotsu[] = { 0.0, 0.0, 0.0, 2.0, 2.0, 2.0};
    double knotsv[] = { 0.0, 0.0, 0.0, 2.0, 2.0, 2.0};
    double scoefs[] = { -1.0, -1.0, 1.0, 1.0,
			0.0, -1.0, 0.5, 2.0,
			1.0, -1.0, 1.0, 1.0,
			-1.0, 0.0, 2.0, 0.5,
			0.0, 0.0, 1.0, 1.0,
			1.0, 0.0, 1.5, 2.0,
			-1.0, 1.0, 0.5, 1.0,
			0.0, 1.0, 2.0, 0.5,
			1.0, 1.0, 1.5, 1.0 };
    const int sdim = 3;
    rational = true;
    SplineSurface surface(nctrlpointsu, nctrlpointsv, orderu, orderv,
			  knotsu, knotsv, scoefs, sdim, rational);
    cout << surface << endl;

    // Convert to Bernstein
    Array<BernsteinMulti, sdim+1> surface_bm;
    splineToBernstein(surface, surface_bm);
    cout << surface_bm[0] << endl 
	 << surface_bm[1] << endl
	 << surface_bm[2] << endl
	 << surface_bm[3] << endl << endl;

    cout << "*** bernsteinToSpline - surfaces ***" << endl;
    cout << endl;

    // Convert again to SplineSurface
    SplineSurface surface2;
    bernsteinToSpline(surface_bm, rational, surface2);
    cout << surface2 << endl;

    return 0;
}




