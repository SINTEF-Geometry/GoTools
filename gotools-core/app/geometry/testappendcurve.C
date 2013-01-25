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
    ALWAYS_ERROR_IF(argc != 6,
		"Arguments: Infile1 Infile2 order continuity Outfile");

    SplineCurve curve, other_curve;
    ObjectHeader header;

    infile >> header >> curve;
    other_curve = curve;
    ifstream infile2(argv[2]);
    ALWAYS_ERROR_IF(infile2.bad(), "Wrong filename or corrupted file.");
    infile2 >> header >> other_curve;

    int order = atoi(argv[3]);
    int order1 = curve.order();
    int other_order = other_curve.order();
    order = std::max(order, order1);
    order = std::max(order, other_order);
    if (order1 < order) {
      curve.raiseOrder(order - order1);
    } 
    if (other_order < order) {
      other_curve.raiseOrder(order - other_order);
    }

    // double par = 0.5*(curve.startparam()+curve.endparam());
    // if (/*curve.rational() && */curve.order() == curve.numCoefs())
    //   curve.insertKnot(par);
    // par = 0.5*(other_curve.startparam()+other_curve.endparam());
    // if (/*other_curve.rational() && */other_curve.order() == other_curve.numCoefs())
    //   other_curve.insertKnot(par);
#ifdef GEOMETRY_DEBUG
   // We insert knot into second cv. Rather case specific...
   double new_knot = other_curve.startparam() +
     0.00001*(other_curve.endparam()-other_curve.startparam());
   other_curve.insertKnot(new_knot);
   other_curve.insertKnot(new_knot);
#endif // GEOMETRY_DEBUG

    double dist = 0;
    curve.appendCurve(&other_curve, atoi(argv[4]), dist, true);
    cout << "Estimated difference between original and smooth curve: " << 
	dist << endl;

    ofstream outfile(argv[5]);
    outfile << header << curve;	

    return 0;

}
