//===========================================================================
//
// File : test_append.C
//
// Created: Tue Dec  2 10:17:40 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_append.C,v 1.1 2008-12-02 13:21:56 kfp Exp $
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
    ALWAYS_ERROR_IF(argc != 7, "Usage: " << argv[0]
		    << " volumeinfile appendvolumeinfile outfile par_dir cont(0 or 1) repar(0 or 1)" << endl);

    // Open input volume file
    ifstream is1(argv[1]);
    ALWAYS_ERROR_IF(is1.bad(), "Bad or no input filename");

    // Open input append volume file
    ifstream is2(argv[2]);
    ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename (append)");

    // Open output volume swap file 1
    ofstream os(argv[3]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    // Read volume 1 from file
    SplineVolume vol1;
    ObjectHeader head;
    is1 >> head >> vol1;

    // Read volume 2 from file
    SplineVolume *vol2 = new SplineVolume();;
    is2 >> head >> *vol2;

    int join_dir = atoi(argv[4]);
    int cont = atoi(argv[5]);
    bool repair = atoi(argv[6]) != 0;

    // Append
    vol1.appendVolume(vol2, join_dir, cont, repair);

    vol1.writeStandardHeader(os);
    vol1.write(os);
}
