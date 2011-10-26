//===========================================================================
//                                                                           
// File: testMultCurveWithFunction.C                                         
//                                                                           
// Created: Wed Aug  8 10:38:35 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: testMultCurveWithFunction.C,v 1.5 2005-09-22 17:20:59 sbr Exp $
//                                                                           
// Description: Testing specified function. Input: alpha, f.
//                                                                           
//===========================================================================
#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/utils/Point.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;

using namespace Go;

int main(int argc, char* argv[])
{
    ALWAYS_ERROR_IF(argc != 4, "Expecting 3 arguments (alpha f out_curve).");

    SplineCurve alpha, f, *product_curve;
    ObjectHeader header;

    ifstream infile1(argv[1]);
    infile1 >> header >> alpha;

    ifstream infile2(argv[2]);
    infile2 >> header >> f;

    product_curve = CurveCreators::multCurveWithFunction(alpha, f);

    ofstream outfile(argv[3]);
    outfile << header << *product_curve;

    double tpar = 0;
    Point alpha_point, f_point, product_point;
    while (tpar != 1001) {
	cout << "Specify tpar to test difference between curves: ";
	cin >> tpar;
	MESSAGE_IF((tpar < f.startparam()) || (tpar > f.endparam()),
		      "Testing outside valid parameter interval!");
	alpha.point(alpha_point, tpar);
	f.point(f_point, tpar);
	double scalar = alpha_point[0];
	Point actual_product_point = f_point * scalar;
	product_curve->point(product_point, tpar);
	cout << " Distance evaluated in tpar: " << 
	    actual_product_point.dist(product_point) << endl;
    }

    return 0;

}



