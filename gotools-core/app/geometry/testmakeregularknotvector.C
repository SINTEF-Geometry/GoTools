//===========================================================================
//                                                                           
// File: testregularknotvector.C                                             
//                                                                           
// Created: Tue May 22 15:34:51 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: testmakeregularknotvector.C,v 1.4 2005-06-09 07:15:43 oan Exp $
//                                                                           
// Description: Test procedure for producing k-regularity at start and end
//              of knot vector.
//              To perform test: Required to move procedures to public area!
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
    ALWAYS_ERROR_IF(argc != 3,
		    "Expecting 2 arguments (infile outfile).");


    ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(),
		    "Wrong filename or corrupted file.");


    SplineCurve the_curve;
    ObjectHeader header;

    infile >> header >> the_curve;
    SplineCurve regular_curve = the_curve;

    // These procedures must be made public (located in private area)!
    regular_curve.makeKnotStartRegular();
    regular_curve.makeKnotEndRegular();

    ofstream outfile(argv[2]);
    outfile << header << regular_curve;	

    // testing, quit when input = 1001
    Point the_point, other_point;
    double tpar = 0.0;
    while (tpar != 1001) { 
	cout << "Specify test-point (" << the_curve.basis().startparam() <<
	    " <= t <= " << the_curve.basis().endparam() << ") : ";
	cin >> tpar;
	the_curve.point(the_point, tpar);
	regular_curve.point(other_point, tpar);
	cout << the_point.dist(other_point) << endl;
    }


    return 0;

}

