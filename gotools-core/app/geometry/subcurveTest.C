//===========================================================================
//                                                                           
// File: subcurveTest.C                                                        
//                                                                           
// Created: Tue Apr 24 17:32:06 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: subcurveTest.C,v 1.4 2003-05-08 15:37:15 afr Exp $
//                                                                           
// Description: Testing the function subCurve
//                                                                           
//===========================================================================

#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <string>

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;

using namespace Go;

int main()
{

    SplineCurve the_curve;
    ObjectHeader header;
    string file_name, file_name2;
    cout << "Enter name of curve-file: ";
    cin >> file_name;
    ifstream infile(file_name.c_str());
    infile >> header >> the_curve;

    cout << "Specify t-parameters (left and right): ";
    double tparl, tparr;
    cin >> tparl >> tparr;

    SplineCurve* the_subcurve = the_curve.subCurve(tparl, tparr);

    cout << "Enter name of output-file: ";
    cin >> file_name2;
    ofstream outfile(file_name2.c_str());
    outfile << header << *the_subcurve;

    // testing, quit when input = 1001
    Point the_point, the_subPoint;
    double tpar = 0.0;

    while (tpar != 1001) { 
	cout << "Specify test-point: ";
	cin >> tpar;
	the_curve.point(the_point, tpar);
	the_subcurve->point(the_subPoint, tpar);
	cout << (the_point - the_subPoint) << endl;
    }

    return 0;

}
