//===========================================================================
//                                                                           
// File: grid_area.C                                                         
//                                                                           
// Created: Thu Feb  1 15:15:17 2007                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: grid_area.C,v 1.1 2007-03-08 10:45:17 afr Exp $
//                                                                           
//===========================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/RectGrid.h"

using namespace std;
using namespace Go;

int main(int argc, char** argv)
{
    if (argc < 2) {
	cerr << "Usage: " << argv[0] << " method (1 = quads, 2 = triangles)" << endl;
	return 1;
    }
    int method = atoi(argv[1]);

    ObjectHeader head;
    cin >> head;
    ASSERT(head.classType() == RectGrid::classType());
    RectGrid grid;
    cin >> grid;

    int numu = grid.numCoefs_u();
    int numv = grid.numCoefs_v();
    double area = 0.0;
    const double* p = grid.rawData();
    Vector3D q1;
    Vector3D q2;
    Vector3D q3;
    Vector3D q4;
    for (int i = 0; i < numu - 1; ++i) {
	for (int j = 0; j < numv - 1; ++j) {
	    Vector3D p0(p + 3*(i + j*numu));
	    Vector3D p1(p + 3*(i + j*numu + 1));
	    Vector3D p2(p + 3*(i + (j + 1)*numu + 1));
	    Vector3D p3(p + 3*(i + (j + 1)*numu));
	    switch(method) {
	    case 1:
		q1 = p2;
		q1 -= p0;
		q2 = p3;
		q2 -= p1;
		q3 = q1.cross(q2);
		area += 0.5*q3.length();
		break;
	    case 2:
		q1 = p1;
		q1 -= p0;
		q2 = p2;
		q2 -= p0;
		q3 = q1.cross(q2);
		area += 0.5*q3.length();
		q1 = p3;
		q1 -= p0;
		q3 = q2.cross(q1);
		area += 0.5*q3.length();
		break;
	    case 3:
		q1 = p1;
		q1 -= p0;
		q2 = p2;
		q2 -= p1;
		q3 = p3;
		q3 -= p2;
		q4 = p0;
		q4 -= p3;
		area += 0.25*(q1.cross(q4)).length();
		area += 0.25*(q2.cross(q3)).length();
		area += 0.25*(q1.cross(q2)).length();
		area += 0.25*(q3.cross(q4)).length();
		break;
	    case 4:
		q1 = p1;
		q1 -= p0;
		q2 = p2;
		q2 -= p1;
		q3 = p3;
		q3 -= p2;
		q4 = p0;
		q4 -= p3;
		area += 0.5*(q1.cross(q4)).length();
		area += 0.5*(q2.cross(q3)).length();
		//area += 0.25*(q1.cross(q2)).length();
		//area += 0.25*(q3.cross(q4)).length();
		break;
	    default:
		THROW("No such method: " << method);
	    }
	}
    }
    cout.precision(18);
    cout << area << endl;
}
