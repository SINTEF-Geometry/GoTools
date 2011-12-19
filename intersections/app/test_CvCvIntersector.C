//==========================================================================
//                                                                          
// File: test_CvCvIntersector.C
//
// Created:
//                                                                          
// Author: B. Spjelkavik <bsp@sintef.no>
//                                                                          
// Revision: $Id: test_CvCvIntersector.C,v 1.9 2006-03-09 15:16:53 jbt Exp $
//                                                                          
// Description: Test of class CvCvIntersector
//                                                                          
//==========================================================================


#include "GoTools/intersections/CvCvIntersector.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/SplineCurveInt.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>
#include <iomanip>


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
 
    if (argc != 4) {
	cout << "Usage: test_CvCvIntersector FileCv1 FileCv2 aepsge"
	     << endl;
	return 0;
    }


    ObjectHeader header;

    // Read the first curve from file
    ifstream input1(argv[1]);
    if (input1.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input1);
    shared_ptr<ParamCurve> curve1(new SplineCurve());
    curve1->read(input1);
    input1.close();
    
    // Read the second curve from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input2);
    shared_ptr<ParamCurve> curve2(new SplineCurve());
    curve2->read(input2);
    input2.close();

    double aepsge;
    aepsge = atof(argv[3]);

//     cout << "\nFile : " << argv[2] << " Parameter range u: "
// 	 << curve1->startparam_u() <<" " << curve1->endparam_u()
// 	 << " Parameter range v: " << curve1->startparam_v() <<" " 
// 	 << curve1->endparam_v();

//     cout << "\nFile : " << argv[2] << " Parameter range u: "
// 	 << curve2->startparam_u() <<" " << curve2->endparam_u()
// 	 << " Parameter range v: " << curve2->startparam_v() <<" " 
// 	 << curve2->endparam_v() << endl;

    shared_ptr<ParamGeomInt> scurveint1 =
	shared_ptr<ParamGeomInt>(new SplineCurveInt (curve1));
    shared_ptr<ParamGeomInt> scurveint2 =
	shared_ptr<ParamGeomInt>(new SplineCurveInt (curve2));

    CvCvIntersector cvcvintersect (scurveint1, scurveint2, aepsge);
    cvcvintersect.compute();

    std::vector<shared_ptr<IntersectionPoint> > intpts;
    std::vector<shared_ptr<IntersectionCurve> > intcrv;
    cvcvintersect.getResult(intpts, intcrv);
    printf("Number of points: %d \n", int(intpts.size()));
    printf("Number of curves: %d \n", int(intcrv.size()));

    int ki, kj;
    for (ki=0; ki < int(intpts.size()); ki++) {
	std::vector<double> par = intpts[ki]->getPar();
	for (kj=0; kj<int(par.size()); kj++)
	    std::cout << par[kj] << " ";
	std::cout << std::endl;
    }

    return 0;
}
