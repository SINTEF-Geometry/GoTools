//==========================================================================
//                                                                          
// File: test_curveSplineCurveInt.C
//
// Created:
//                                                                          
// Author: B. Spjelkavik <bsp@sintef.no>
//                                                                          
// Revision: $Id: test_SplineCurveInt.C,v 1.4 2006-03-08 15:25:39 jbt Exp $
//                                                                          
// Description: Test of classes ParamCurveInt and SplineCurveInt
//                                                                          
//==========================================================================


#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include <memory>
#include <fstream>
#include <iomanip>


using std::ifstream;
using std::cout;
using std::endl;
using std::cerr;
using namespace Go;


int main(int argc, char** argv)
{
    if (argc != 7) {
	cout << "Usage: test_curveSplineCurveInt FILE1 FILE2 "
	     << "pstart1 pstart2 pstop1 pstop2"
	     << endl;
	return 0;
    }
    double pstart1, pstart2, pstop1, pstop2;
    pstart1 = atof(argv[3]);
    pstart2 = atof(argv[4]);
    pstop1 = atof(argv[5]);
    pstop2 = atof(argv[6]);
    ObjectHeader header;

    // Read first curve from file
    ifstream input(argv[1]);
    if (input.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input);
    shared_ptr<SplineCurve> curve1(new SplineCurve());
    curve1->read(input);
    input.close();

   
    
    // Read second curve from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input2);
    shared_ptr<SplineCurve> curve2(new SplineCurve());
    curve2->read(input2);
    input.close();

//     int prec = std::cout.precision(16);
    cout << "\nFile : " << argv[1] << " Parameter range: "
	 << curve1->startparam() <<" " << curve1->endparam();
    cout << "  Start: " << pstart1 << "  Stop: " << pstop1;
    cout << "\nFile : " << argv[2] << " Parameter range: "
	 << curve2->startparam() <<" " << curve2->endparam();
    cout << "  Start: " << pstart2 << "  Stop: " << pstop2 << endl;

    const double eps = 1.e-4;
    cout << "  Eps= " << eps << endl;

//     SplineCurveInt sci1(curve1);

//     SplineCurveInt sci2(curve2);

//     int istat;
//     istat = sci1.checkCoincidence(pstart1, pstop1, eps, &sci2,
// 				  pstart2, pstop2);

//     cout << "\nistat = " << istat;

//     if (istat==0)
// 	cout << "  Curves are not coinciding." << endl;
//     else
// 	cout << "  Curves are coinciding." << endl;

//     std::cout.precision(prec);  

    return 0;
}
