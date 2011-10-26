//===========================================================================
//
// File : test_subVolume.C
//
// Created: Wed Nov 26 09:43:54 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_subVolume.C,v 1.1 2008-11-27 08:35:17 kfp Exp $
//
// Description:
//
//===========================================================================



#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;



int main(int argc, char* argv[] )
{

    ALWAYS_ERROR_IF(argc != 9, "Usage: " << argv[0]
		    << " volumeinfile volumeoutfile from_u from_v from_w to_u to_v to_w" << endl);

    // Open input volume file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Open output volume swap file 1
    ofstream os(argv[2]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    // Read volume from file
    SplineVolume vol;
    ObjectHeader head;
    is >> head >> vol;

    SplineVolume* vol2 = vol.subVolume(atof(argv[3]),
				       atof(argv[4]),
				       atof(argv[5]),
				       atof(argv[6]),
				       atof(argv[7]),
				       atof(argv[8]));

    vol2->writeStandardHeader(os);
    vol2->write(os);
}
