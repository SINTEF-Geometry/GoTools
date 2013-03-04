//==========================================================================
//                                                                          
// File: test_BoundedSurface.C                                                   
//                                                                          
// Created: Thu Jan  7 10:11:54 2010                                         
//                                                                          
// Author: Jan B. Thomassen <jan.b.thomassen@sintef.no>
//                                                                          
// Revision: $Id:$
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include <fstream>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"


using namespace std;
using namespace Go;

int main(int argc, char** argv)
{
    GoTools::init();


    // Example input file: data/bounded_surface.g2

    if (argc != 2) {
        cout << "Usage:  " << argv[0] << " FILE" << endl;
        return 0;
    }

    // Read the surface from file
    ifstream input(argv[1]);
    if (input.bad()) {
        cout << "File error (no file or corrupt file specified)."
                  << endl;
        return 1;
    }
    ObjectHeader header;
    BoundedSurface bs;
    input >> header >> bs;

    ofstream out1("data/input.g2");
    bs.writeStandardHeader(out1);
    bs.write(out1);


    bs.swapParameterDirection();

    ofstream out2("data/output.g2");
    bs.writeStandardHeader(out2);
    bs.write(out2);




    return 0;
}
