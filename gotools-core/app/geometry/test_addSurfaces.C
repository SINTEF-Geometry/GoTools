//===========================================================================
//
// File : test_addSurfaces.C
//
// Created: Thu May 27 12:47:04 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
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

    ALWAYS_ERROR_IF(argc != 4, "Usage: " << argv[0]
		    << " surf1_infile surf2_infile surfaceoutfile" << endl);

    // Open input surface 1 file
    ifstream is1(argv[1]);
    ALWAYS_ERROR_IF(is1.bad(), "Bad or no surface 1 input filename");

    // Open input surface 2 file
    ifstream is2(argv[2]);
    ALWAYS_ERROR_IF(is2.bad(), "Bad or no surface 2 input filename");

    // Open output surface file
    ofstream os(argv[3]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    // Read surfaces from file
    SplineSurface *surf1, *surf2;
    ObjectHeader head;

    is1 >> head;
    ALWAYS_ERROR_IF(head.classType() != SplineSurface::classType(), "First argument is not a spline surface file");
    surf1 = new SplineSurface;
    surf1->read(is1);

    is2 >> head;
    ALWAYS_ERROR_IF(head.classType() != SplineSurface::classType(), "Second argument is not a spline surface file");
    surf2 = new SplineSurface;
    surf2->read(is2);

    surf1->add(surf2);

    surf1->writeStandardHeader(os);
    surf1->write(os);
}
