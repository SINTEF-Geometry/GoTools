//===========================================================================
//
// File : test_rotationalsweep.C
//
// Created: Mon Dec  8 10:42:00 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_rotationalsweep.C,v 1.1 2008-12-09 11:24:08 kfp Exp $
//
// Description:
//
//===========================================================================



#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;


int main(int argc, char* argv[] )
{

    ALWAYS_ERROR_IF(argc != 10, "Usage: " << argv[0]
		    << " curveinfile surfaceoutfile angle point_x point_y point_z axis_x axis_y axis_z" << endl);

    // Open input curve file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no curve input filename");

    // Open output surface file
    ofstream os(argv[2]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    double angle = atof(argv[3]);
    Point pt(atof(argv[4]), atof(argv[5]), atof(argv[6]));
    Point axis(atof(argv[7]), atof(argv[8]), atof(argv[9]));

    // Read curve from file
    SplineCurve curve;
    ObjectHeader head;
    is >> head >> curve;

    SplineSurface* surf = SweepSurfaceCreator::rotationalSweptSurface(curve, angle, pt, axis);

    surf->writeStandardHeader(os);
    surf->write(os);
}
