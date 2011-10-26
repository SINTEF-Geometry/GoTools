//===========================================================================
//                                                                           
// File: testraiseorder.C                                                        
//                                                                           
// Created: Tue May 15 09:15:10 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision:
//                                                                           
// Description: Test the function raiseOrder() for a curve.
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
    if(argc != 3) {
	std::cout << "usage: filein fileout" << std::endl;
	return -1;
    }

    SplineCurve the_curve, raised_curve;
    ObjectHeader header;

    ifstream infile(argv[1]);
    infile >> header >> the_curve;

    raised_curve = the_curve;
    raised_curve.raiseOrder();

    ofstream outfile(argv[2]);
    outfile << header << raised_curve;

    double tpar = 0;
    Point the_tpoint, raised_tpoint;
    while (tpar != 1001) {
	cout << "Specify tpar to test difference between curves: ";
	cin >> tpar;
	the_curve.point(the_tpoint, tpar);
	raised_curve.point(raised_tpoint, tpar);
	cout << " Distance evaluated in tpar: " << 
	    raised_tpoint.dist(the_tpoint) << endl;
    }

    return 0;

}
