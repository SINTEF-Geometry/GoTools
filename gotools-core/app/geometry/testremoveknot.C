//===========================================================================
//                                                                           
// File: testremoveknot.C                                                        
//                                                                           
// Created: Fri May 3 17:32:06 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision:
//                                                                           
// Description: Test the functions appendCurve, removeKnot & raiseOrder.
//                                                                           
//===========================================================================

#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
//#include "GoTools/geometry/GoMakeGnuplot.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;

using namespace Go;

int main(int argc, char* argv[])
{
    ALWAYS_ERROR_IF(argc != 4,"Expecting 3 arguments.");


    SplineCurve the_curve, smooth_curve;
    ObjectHeader header;
    SplineCurve *the_curvep=new SplineCurve;

    ifstream infile2(argv[1]);
    infile2 >> header >> (*the_curvep);

    ifstream infile(argv[1]);
    infile >> header >> the_curve;

    smooth_curve = the_curve;
    smooth_curve.removeKnot(atof(argv[2]));

    ofstream outfile(argv[3]);
    outfile << header <<smooth_curve;

    double tpar = 0;
    Point the_tpoint, smooth_tpoint;
    while (tpar != 1001) {
	cout << "Specify tpar to test difference between curves: ";
	cin >> tpar;
	the_curve.point(the_tpoint, tpar);
	smooth_curve.point(smooth_tpoint, tpar);
	cout << " Distance evaluated in tpar: " << 
	    smooth_tpoint.dist(the_tpoint) << endl;
    }

    return 0;

}
