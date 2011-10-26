//==========================================================================
//                                                                          
// File: example_BernsteinUtils.C                                            
//                                                                          
// Created: Tue Mar 28 14:18:46 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: example_BernsteinUtils.C,v 1.1 2006-03-29 12:16:24 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


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




