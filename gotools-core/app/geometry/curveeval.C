//===========================================================================
//                                                                           
// File: curveeval.C                                                         
//                                                                           
// Created: Fri Oct 13 16:34:04 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: curveeval.C,v 1.7 2005-06-09 07:15:39 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
    if (argc < 2) {
	std::cerr << "Usage: " << argv[0]
		  << " inputsurface" << std::endl;
	return 1;
    }

    // Open input curve file
    std::ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read curve from file
    ObjectHeader header;
    is >> header;
    SplineSurface sf;
    is >> sf;

    // Input point to evaluate
    double u, v;
    std::cout << "Please enter a (u, v) coordiante pair: ";
    std::cin >> u >> v;

    Point p(sf.dimension());
    for (int i = 0; i < 1000; ++i) {
	sf.point(p, u, v);
    }
}





