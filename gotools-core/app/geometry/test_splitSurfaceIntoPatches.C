//===========================================================================
//                                                                           
// File: test_splitSurfaceIntoPatches.C                                   
//                                                                           
// Created: Tue Aug 14 14:32:20 2001                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: test_splitSurfaceIntoPatches.C,v 1.7 2005-06-10 09:23:54 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>
#include <iomanip>

using namespace Go;
using namespace std;



int main(int argc, char** argv)
{
    if (argc != 2) {
	cout << "Usage:  " << argv[0] << " FILE" << endl;
	return 0;
    }

    // Read the surface from file
    std::ifstream input(argv[1]);
    if (input.bad()) {
	std::cerr << "File error (no file or corrupt file specified)."
		  << std::endl;
	return 1;
    }
    ObjectHeader header;
    SplineSurface surface;
    input >> header >> surface;

//      surface.writeStandardHeader(std::cout);
//      std::cout << surface << std::endl;

    std::vector<SplineSurface> patches;
    GeometryTools::splitSurfaceIntoPatches(surface, patches);

    cout << "There are " << patches.size() << " patches" << endl;

    // Write out n'th patch
    const int n = 0;
    patches[n].writeStandardHeader(std::cout);
    std::cout << patches[n] << std::endl;

    return 0;
}
