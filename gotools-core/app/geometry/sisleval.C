//===========================================================================
//                                                                           
// File: sisleval.C                                                          
//                                                                           
// Created: Fri Sep 29 17:33:49 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: sisleval.C,v 1.4 2005-11-09 10:22:38 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include <fstream>
#include <vector>
//#include "GoTools/geometry/sisl.h"
#include "GoTools/geometry/SISL_code.h"
#include "GoTools/utils/errormacros.h"

//using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc < 3, "Usage: " << argv[0]
		    << " inputsurf inputpoints" << endl);


    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read surface from file
    int in1, in2, ik1, ik2, dim;
    std::vector<double> et1, et2, co;
    is >> dim;
    is >> in1 >> ik1;
    et1.resize(in1+ik1);
    for (int i = 0; i < in1+ik1; ++i)
	is >> et1[i];
    is >> in2 >> ik2;
    et2.resize(in2+ik2);
    for (int i = 0; i < in2+ik2; ++i)
	is >> et2[i];
    co.resize(in1*in2*dim);
    for (int i = 0; i < in1*in2*dim; ++i)
	is >> co[i];
    SISLSurf* sf = newSurf(in1, in2, ik1, ik2, &et1[0], &et2[0],
			   &co[0], 1, dim, 1);

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*2);
    for (int i = 0; i < n; ++i) {
	pts >> pt[2*i] >> pt[2*i+1];
    }

    int il1, il2, stat;
    double der[9];
    double nor[3];
    for (int i = 0; i < 10000; ++i) {
	for (int j = 0; j < n; ++j) {
	    s1421(sf, 1, &pt[2*j], &il1, &il2, der, nor, &stat);
	}
    }
    //    cout << p[0] << p[1] << p[2] << (p[1] % p[2]);
    freeSurf(sf);
}




    /*
    int vnum = basis_v_.numCoefs();

    std::vector<double> b0s(b0);
    std::vector<double> b1s(b1);
    int stat = 0;
    int ilfs = 0;
    int ilft = 0;
    s1220(basis_u_.begin(), uorder, unum, &uleft, upar, 1, b0s.begin(), &stat);

    SISLSurf* ps1 = newSurf(unum, vnum, uorder, vorder, basis_u_.begin(),
			    basis_v_.begin(), coefs_.begin(), 1, 3, 1);

    double par[2]; par[0] = upar; par[1] = vpar;
    double eder1[9];
    double enorm1[3];
    double eder2[9];
    double enorm2[3];
    s1421(ps1, 1, par, &ilfs, &ilft,
	  eder1, enorm1, &stat);
    s1421(ps1, 0, par, &ilfs, &ilft,
	  eder2, enorm2, &stat);
    */
