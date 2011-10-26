//===========================================================================
//
// File : test_linearsweep.C
//
// Created: Thu Dec  4 09:29:10 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_linearsweep.C,v 1.1 2008-12-04 09:03:14 kfp Exp $
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

    ALWAYS_ERROR_IF(argc != 7, "Usage: " << argv[0]
		    << " surfaceinfile curveinfile volumeoutfile point_x point_y point_z" << endl);

    // Open input surface file
    ifstream is_surf(argv[1]);
    ALWAYS_ERROR_IF(is_surf.bad(), "Bad or no surface input filename");

    // Open input curve file
    ifstream is_crv(argv[2]);
    ALWAYS_ERROR_IF(is_crv.bad(), "Bad or no curve input filename");

    // Open output volume file
    ofstream os(argv[3]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    Point pt(atof(argv[4]), atof(argv[5]), atof(argv[6]));

    // Read surface from file
    SplineSurface surf;
    ObjectHeader head;
    is_surf >> head >> surf;

    // Read curve from file
    SplineCurve curve;
    is_crv >> head >> curve;

    SplineVolume* vol = SweepVolumeCreator::linearSweptVolume(surf, curve, pt);

    vol->writeStandardHeader(os);
    vol->write(os);
}
