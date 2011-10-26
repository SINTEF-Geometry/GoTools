//===========================================================================
//
// File : test_removeKnot.C
//
// Created: Tue Dec  2 13:43:56 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_removeKnot.C,v 1.1 2008-12-02 13:21:56 kfp Exp $
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
    ALWAYS_ERROR_IF(argc != 5, "Usage: " << argv[0]
		    << " volumeinfile volumeoutfile pardir knotvalue" << endl);

    // Open input volume file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Open output volume file
    ofstream os(argv[2]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    // Read volume from file
    SplineVolume vol;
    ObjectHeader head;
    is >> head >> vol;

    int pardir = atoi(argv[3]);
    double tpar = atof(argv[4]);
    vol.removeKnot(pardir, tpar);

    vol.writeStandardHeader(os);
    vol.write(os);
}
