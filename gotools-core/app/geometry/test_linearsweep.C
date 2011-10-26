//===========================================================================
//
// File : test_linearsweep.C
//
// Created: Wed Dec  3 13:48:02 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_linearsweep.C,v 1.1 2008-12-04 09:04:13 kfp Exp $
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

    ALWAYS_ERROR_IF(argc != 7, "Usage: " << argv[0]
		    << " curve1infile curve2infile surfaceoutfile point_x point_y point_z" << endl);

    // Open input curve 1 file
    ifstream is1(argv[1]);
    ALWAYS_ERROR_IF(is1.bad(), "Bad or no curve 1 input filename");

    // Open input curve 2 file
    ifstream is2(argv[2]);
    ALWAYS_ERROR_IF(is2.bad(), "Bad or no curve 2 input filename");

    // Open output surface file
    ofstream os(argv[3]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    Point pt(atof(argv[4]), atof(argv[5]), atof(argv[6]));

    // Read curves from file
    SplineCurve curv1, curv2;
    ObjectHeader head;
    is1 >> head >> curv1;
    is2 >> head >> curv2;

    SplineSurface* surf = SweepSurfaceCreator::linearSweptSurface(curv1, curv2, pt);

    surf->writeStandardHeader(os);
    surf->write(os);
}
