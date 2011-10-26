//==========================================================================
//                                                                          
// File: test_makeBernsteinKnots.C                                           
//                                                                          
// Created: Fri Jun 10 10:11:52 2005                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: test_makeBernsteinKnots.C,v 1.1 2005-06-10 08:24:33 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char** argv)
{
    if (argc != 2) {
	cout << "Usage:  " << argv[0] << " FILE" << endl;
	return 0;
    }

    // Read the curve from file
    ifstream input1(argv[1]);
    if (input1.bad()) {
	cerr << "File error (no file or corrupt file specified)."
	     << endl;
	return 1;
    }
    ObjectHeader header;

    SplineCurve cv;
    input1 >> header >> cv;

    cv.makeBernsteinKnots();
    cv.writeStandardHeader(cout);
    cv.write(cout);

//     SplineSurface surface;
//     input1 >> header >> surface;

//     surface.makeBernsteinKnotsU();
//     surface.makeBernsteinKnotsV();
//     surface.writeStandardHeader(cout);
//     surface.write(cout);

    return 0;
}







