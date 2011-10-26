//==========================================================================
//                                                                          
// File: direction_cone_test.C                                               
//                                                                          
// Created: Fri Feb  1 12:05:36 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: direction_cone_test.C,v 1.3 2003-05-08 15:37:14 afr Exp $
//                                                                          
// Description:
//                                                                          
//===========================================================================


#include "GoTools/utils/DirectionCone.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char** argv)
{
    if (argc != 2) {
	cout << "Usage: direction_cone_test FILE" << endl;
	return 0;
    }

    // Read the spline object from file
    ifstream input(argv[1]);
    if (input.bad()) {
	cerr << "File error (no file or corrupt file specified)."
	     << endl;
	return 1;
    }
    ObjectHeader header;
//     SplineCurve curve;
//     input >> header >> curve;
    SplineSurface surf;
    input >> header >> surf;

    // Make direction cone
//     DirectionCone cone = curve.directionCone();
    DirectionCone cone = surf.normalCone();

    // Write out
    cone.write(cout);
    cout << endl;

    return 0;
}
