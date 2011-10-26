//===========================================================================
//
// File : test_addVolumes.C
//
// Created: Thu May 27 15:14:21 2010
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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;


int main(int argc, char* argv[] )
{

    ALWAYS_ERROR_IF(argc != 4, "Usage: " << argv[0]
		    << " vol1_infile vol2_infile volumeoutfile" << endl);

    // Open input volume 1 file
    ifstream is1(argv[1]);
    ALWAYS_ERROR_IF(is1.bad(), "Bad or no volume 1 input filename");

    // Open input volume 2 file
    ifstream is2(argv[2]);
    ALWAYS_ERROR_IF(is2.bad(), "Bad or no volume 2 input filename");

    // Open output volume file
    ofstream os(argv[3]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    // Read volumes from file
    SplineVolume *vol1, *vol2;
    ObjectHeader head;

    is1 >> head;
    ALWAYS_ERROR_IF(head.classType() != SplineVolume::classType(), "First argument is not a spline volume file");
    vol1 = new SplineVolume;
    vol1->read(is1);

    is2 >> head;
    ALWAYS_ERROR_IF(head.classType() != SplineVolume::classType(), "Second argument is not a spline volume file");
    vol2 = new SplineVolume;
    vol2->read(is2);

    vol1->add(vol2);

    vol1->writeStandardHeader(os);
    vol1->write(os);
}
