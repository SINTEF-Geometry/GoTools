//===========================================================================
//                                                                           
// File: multieval.C                                                         
//                                                                           
// Created: Thu Sep  7 15:27:04 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: multieval.C,v 1.5 2005-06-09 07:15:40 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <fstream>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc < 3, "Usage: " << argv[0]
		    << " inputsurf inputpoints" << endl);


    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read surface from file
    SplineSurface sf;
    is >> sf;

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*2);
    for (int i = 0; i < n; ++i) {
	pts >> pt[2*i] >> pt[2*i+1];
    }

    std::vector<Point> p(3, Point(sf.dimension()));
    for (int i = 0; i < 10000; ++i) {
	for (int j = 0; j < n; ++j) {
	    sf.point(p, pt[2*j], pt[2*j+1], 0);
	}
    }
//      cout << p[0] << p[1] << p[2] << (p[1] % p[2]);
}





