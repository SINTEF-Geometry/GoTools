//===========================================================================
//
// File : test_clone_swap.C
//
// Created: Fri Nov 21 12:07:04 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_clone_swap.C,v 1.1 2008-11-21 11:33:29 kfp Exp $
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
    ALWAYS_ERROR_IF(argc != 6, "Usage: " << argv[0]
		    << " volumeinfile1 volumeinfile2 outfileswap1 outfileswap2 outfileclone" << endl);

    // Open input volume file 1
    ifstream is1(argv[1]);
    ALWAYS_ERROR_IF(is1.bad(), "Bad or no input 1 filename");

    // Open input volume file 2
    ifstream is2(argv[2]);
    ALWAYS_ERROR_IF(is2.bad(), "Bad or no input 2 filename");

    // Open output volume swap file 1
    ofstream os_swap1(argv[3]);
    ALWAYS_ERROR_IF(os_swap1.bad(), "Bad output swap 1 filename");

    // Open output volume swap file 2
    ofstream os_swap2(argv[4]);
    ALWAYS_ERROR_IF(os_swap2.bad(), "Bad output swap 2 filename");

    // Open output volume clone file
    ofstream os_clone(argv[5]);
    ALWAYS_ERROR_IF(os_clone.bad(), "Bad output clone filename");

    // Read volume 1 from file
    SplineVolume vol1;
    ObjectHeader head;
    is1 >> head >> vol1;

    // Read volume 2 from file
    SplineVolume vol2;
    is2 >> head >> vol2;

    // Volume for clone
    SplineVolume* vol_clone = vol1.clone();
    vol1.swap(vol2);

    vol1.writeStandardHeader(os_swap1);
    vol1.write(os_swap1);

    vol2.writeStandardHeader(os_swap2);
    vol2.write(os_swap2);

    vol_clone->writeStandardHeader(os_clone);
    vol_clone->write(os_clone);
}
