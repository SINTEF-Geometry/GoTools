//==========================================================================
//                                                                          
// File: test_SISLsfcv.C                                                     
//                                                                          
// Created: Thu Jun 15 13:47:10 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: test_SISLsfcv.C,v 1.1 2006-06-19 09:24:14 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <fstream>
#include <iomanip>


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
    if (argc != 4) {
	cout << "Usage: " << argv[0] << " FileSf FileCv aepsge" << endl;
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
    SplineSurface surf;
    surf.read(input1);
    input1.close();
    
    // Read the second curve from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input2);
    SplineCurve curve;
    curve.read(input2);
    input2.close();

    double aepsge;
    aepsge = atof(argv[3]);

    // Set up and run the intersection algorithm
    SISLCurve* pcurve = Curve2SISL(curve);
    SISLSurf* psurf = GoSurf2SISL(surf);

    double astart1 = curve.startparam();
    double estart2[] = { surf.startparam_u(), surf.startparam_v() };
    double aend1 = curve.endparam();
    double eend2[] = { surf.endparam_u(), surf.endparam_v() };
    cout << "astart1 = " << astart1 << endl
	 << "estart2[] = { " << estart2[0] << ", "
	 << estart2[1] << " }" << endl
	 << "aend1 = " << aend1 << endl
	 << "eend2[] = { " << eend2[0] << ", " 
	 << eend2[1] << " }" << endl;

    double anext1 = astart1;
    double enext2[] = { estart2[0], estart2[1] };
//     double anext1 = 0.0;
//     double enext2[] = { 2.0, 0.0 };
    cout << "anext1 = " << anext1 << endl
	 << "enext2[] = { " << enext2[0] << ", "
	 << enext2[1] << " }" << endl;

    double cpos1;
    double gpos2[2];
    int jstat = 0;
    cout << "jstat = " << jstat << endl;

    s1772(pcurve, psurf, aepsge, astart1, estart2, aend1, eend2,
	  anext1, enext2, &cpos1, gpos2, &jstat);

    // Write the results
    cout << "Results s1772:" << endl
	 << "jstat = " << jstat << endl
	 << "cpos1 = " << cpos1 << endl
	 << "gpos2[] = { " << gpos2[0] << ", "
	 << gpos2[1] << " }" << endl;

    return 0;

}
