//===========================================================================
//                                                                           
// File: test_closestpoint.C                                                   
//                                                                           
// Created: Tue Mar 21 16:44:59 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: test_closestpoint.C,v 1.2 2005-06-27 15:32:49 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
    ASSERT(argc >= 2);
    ifstream file(argv[1]);
    ObjectHeader head;
    SplineSurface sf;
    SplineCurve cv;
    file >> head;
    if (head.classType() == SplineSurface::classType()) { 
	file >> sf;
    } else if (head.classType() == SplineCurve::classType()) {
	file >> cv;
    } else {
	THROW("Unrecognized object type.");
    }
    
    double x, y, z;
    cin >> x >> y >> z;
    Point p(x, y, z);
    double clo_u;
    double clo_v;
    Point clo_pt(3);
    double clo_dist;
    double epsilon = 1e-8;

    if (head.classType() == SplineSurface::classType()) {
	sf.closestPoint(p, clo_u, clo_v, clo_pt, clo_dist, epsilon);
    } else { // curve
	cv.closestPoint(p, cv.startparam(), cv.endparam(), clo_u, clo_pt, clo_dist);
    }

    cout << clo_pt << endl;
}










