//==========================================================================
//                                                                          
// File: deriv_surface_test.C                                                
//                                                                          
// Created: Mon Feb 25 16:01:36 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: deriv_surface_test.C,v 1.4 2005-09-23 11:48:17 oan Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#include <fstream>
#include <iomanip>
// #include "sislP.h"
// #include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
    if (argc != 2) {
	cout << "Usage: deriv_surface_test FILE" << endl;
	return 0;
    }

    // Read the surface from file
    ifstream input(argv[1]);
    if (input.bad()) {
	cerr << "File error (no file or corrupt file specified)."
		  << std::endl;
	return 1;
    }
    ObjectHeader header;
    SplineSurface surface;
    input >> header >> surface;

//     // Call SISL routine s1386
//     SISLSurf* surf = GoSurf2SISL(surface);
//     SISLSurf* newsurf;
//     int stat;
//     s1386(surf, 1, 0, &newsurf, &stat);
//     SplineSurface deriv = *SISLSurf2Go(newsurf);

    // Take derivative
    auto_ptr<SplineSurface> pderiv(surface.derivSurface(1, 0));
    SplineSurface deriv = dynamic_cast<SplineSurface&>(*pderiv);

    // Write out
    ofstream out("data/deriv.g2");
    out << header << deriv;
    cout << "Output written to file data/deriv.g2" << endl;


    return 0;
}
