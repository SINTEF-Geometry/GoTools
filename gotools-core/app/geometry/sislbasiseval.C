//===========================================================================
//                                                                           
// File: sislbasiseval.C                                                     
//                                                                           
// Created: Tue Oct  3 09:11:39 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: sislbasiseval.C,v 1.4 2005-11-09 10:22:37 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include <fstream>
#include <vector>
#include "GoTools/geometry/SISL_code.h"
#include "GoTools/utils/errormacros.h"

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

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*2);
    for (int i = 0; i < n; ++i) {
	pts >> pt[2*i] >> pt[2*i+1];
    }

    double b1[8];
    double b2[8];
    int il1, il2, stat;
    for (int i = 0; i < 10000; ++i) {
	for (int j = 0; j < n; ++j) {
//  	    s1219(et1.begin(), ik1, in1, &il1, pt[2*j], &stat);
//  	    s1219(et2.begin(), ik2, in2, &il2, pt[2*j+1], &stat);
	    s1220(&et1[0], ik1, in1, &il1, pt[2*j], 1, b1, &stat);
	    s1220(&et2[0], ik2, in2, &il2, pt[2*j+1], 1, b2, &stat);
	}
    }
}
