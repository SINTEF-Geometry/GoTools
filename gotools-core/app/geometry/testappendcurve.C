//===========================================================================
//                                                                           
// File: testappendcurve.C                                                        
//                                                                           
// Created: Fri May 3 17:32:06 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision:
//                                                                           
// Description: Test the function appendCurve().
//                                                                           
//===========================================================================

#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;

using namespace Go;

int main(int argc, char* argv[])
{
    ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(), "Wrong filename or corrupted file.");
    ALWAYS_ERROR_IF(argc != 5,
		"Arguments: Infile1 Infile2 continuity Outfile");

    SplineCurve curve, other_curve;
    ObjectHeader header;

    infile >> header >> curve;
    other_curve = curve;
    ifstream infile2(argv[2]);
    ALWAYS_ERROR_IF(infile2.bad(), "Wrong filename or corrupted file.");
    infile2 >> header >> other_curve;

    int order = curve.order();
    int other_order = other_curve.order();
    if (order < other_order) {
      curve.raiseOrder(other_order - order);
    } else if (other_order < order) {
      other_curve.raiseOrder(order - other_order);
    }

#ifdef GEOMETRY_DEBUG
   // We insert knot into second cv. Rather case specific...
   double new_knot = other_curve.startparam() +
     0.00001*(other_curve.endparam()-other_curve.startparam());
   other_curve.insertKnot(new_knot);
   other_curve.insertKnot(new_knot);
#endif // GEOMETRY_DEBUG

    double dist = 0;
    curve.appendCurve(&other_curve, atoi(argv[3]), dist, true);
    cout << "Estimated difference between original and smooth curve: " << 
	dist << endl;

    ofstream outfile(argv[4]);
    outfile << header << curve;	

    return 0;

}
