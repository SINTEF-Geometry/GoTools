//==========================================================================
//                                                                          
// File: test_GGUsplit.C                                                     
//                                                                          
// Created: Wed Jan 22 17:34:44 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: test_GGUsplit.C,v 1.4 2005-02-22 13:16:22 afr Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GeometryTools.h"
#include <fstream>


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
    if (argc != 2) {
	cout << "Usage: test_GGUsplit FILE" << endl;
	return 0;
    }

    // Read the curve from file
    ifstream input(argv[1]);
    if (input.bad()) {
	cerr << "File error (no file or corrupt file specified)."
		  << endl;
	return 1;
    }
    ObjectHeader header;
    SplineCurve curve;
    input >> header >> curve;

    // Split curve
    vector<SplineCurve> seg;
    splitCurveIntoSegments(curve, seg);

    // Write segments
    ofstream output("data/split.g2");
    output << header << seg[0] << endl;
    cout << "Split curve written to data/split.g2" << endl;

    return 0;
}
