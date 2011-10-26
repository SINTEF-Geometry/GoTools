//===========================================================================
//                                                                           
// File: CurveOnSurface_compare.C                                            
//                                                                           
// Created: Mon Oct 14 15:56:56 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CurveOnSurface_compare.C,v 1.2 2006-04-19 09:27:33 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>


using std::cerr;
using std::cout;
using std::endl;
using std::min;
using std::cin;
using namespace Go;
using std::shared_ptr;


int main()
{
    shared_ptr<SplineCurve> sc(new SplineCurve);
    shared_ptr<SplineCurve> pc(new SplineCurve);
    shared_ptr<SplineSurface> sf(new SplineSurface);
    ObjectHeader header;
    cin >> header >> (*sc);
    cin >> header >> (*pc);
    cin >> header >> (*sf);

    CurveOnSurface cv(sf, pc, sc, true);
    const int N = 100;
    Point pt(3);
    Point ppt(2);
    Point pt2(3);
    double t0 = cv.startparam();
    double t1 = cv.endparam();
    cout << "400 1 0 4 255 255 0 255\n" << N << endl;
    double mdiff = 1e100;
    for (int i = 0; i < N; ++i) {
	double fact = i/double(N-1);
	cv.point(pt, t0*(1.0-fact) + t1*fact);
	cout << pt;
	pc->point(ppt, t0*(1.0-fact) + t1*fact);
	sf->point(pt2, ppt[0], ppt[1]);
	mdiff = min(mdiff, pt.dist(pt2));
	//cout << pt2;
    }
    cerr << mdiff << endl;
}
