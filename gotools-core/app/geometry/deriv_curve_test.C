//==========================================================================
//                                                                          
// File: deriv_curve_test.C                                                  
//                                                                          
// Created: Fri Feb 22 16:53:55 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: deriv_curve_test.C,v 1.2 2003-05-08 15:37:14 afr Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#include <fstream>
#include <iomanip>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
    if (argc != 2) {
	cout << "Usage: deriv_curve_test FILE" << endl;
	return 0;
    }

    // Read the curve from file
    ifstream input(argv[1]);
    if (input.bad()) {
	cerr << "File error (no file or corrupt file specified)."
		  << std::endl;
	return 1;
    }
    ObjectHeader header;
    SplineCurve curve;
    input >> header >> curve;

    // Take derivative
    auto_ptr<SplineCurve> pderiv(curve.derivCurve(1));
    SplineCurve deriv = dynamic_cast<SplineCurve&>(*pderiv);

    // Write out
    ofstream out("data/deriv.g2");
    out << header << deriv;
    cout << "Output written to file data/deriv.g2" << endl;


    return 0;
}
