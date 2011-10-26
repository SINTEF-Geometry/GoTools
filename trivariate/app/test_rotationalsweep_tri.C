//===========================================================================
//
// File : test_rotationalsweep.C
//
// Created: Mon Dec  8 12:53:54 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_rotationalsweep.C,v 1.1 2008-12-09 11:23:07 kfp Exp $
//
// Description:
//
//===========================================================================



#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;


int main(int argc, char* argv[] )
{

    ALWAYS_ERROR_IF(argc != 10, "Usage: " << argv[0]
		    << " surfaceinfile volumeoutfile angle point_x point_y point_z axis_x axis_y axis_z" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no surface input filename");

    // Open output volume file
    ofstream os(argv[2]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    double angle = atof(argv[3]);
    Point pt(atof(argv[4]), atof(argv[5]), atof(argv[6]));
    Point axis(atof(argv[7]), atof(argv[8]), atof(argv[9]));

    // Read surface from file
    SplineSurface surface;
    ObjectHeader head;
    is >> head >> surface;

    SplineVolume* vol = SweepVolumeCreator::rotationalSweptVolume(surface, angle, pt, axis);

    vol->writeStandardHeader(os);
    vol->write(os);




    ofstream oslines("helplines.g2");
    oslines << "100 0 9 0" << endl;
    oslines << "3 0" << endl;
    oslines << "2 2" << endl;
    oslines << "0 0 1 1" << endl;
    oslines << "1 0 0" << endl;
    oslines << "1 0 3" << endl;

    oslines << "100 0 9 0" << endl;
    oslines << "3 0" << endl;
    oslines << "2 2" << endl;
    oslines << "0 0 1 1" << endl;
    oslines << cos(angle) << " " << sin(angle) << " 0" << endl;
    oslines << cos(angle) << " " << sin(angle) << " 6" << endl;

    oslines << "100 0 9 0" << endl;
    oslines << "3 0" << endl;
    oslines << "2 2" << endl;
    oslines << "0 0 1 1" << endl;
    oslines << "0.625 0 0" << endl;
    oslines << "0.625 0 3" << endl;

    oslines << "100 0 9 0" << endl;
    oslines << "3 0" << endl;
    oslines << "2 2" << endl;
    oslines << "0 0 1 1" << endl;
    oslines << cos(angle)/1.6 << " " << sin(angle)/1.6 << " 0" << endl;
    oslines << cos(angle)/1.6 << " " << sin(angle)/1.6 << " 6" << endl;
}
