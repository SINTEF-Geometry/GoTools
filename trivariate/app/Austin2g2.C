//===========================================================================
//                                                                           
// File: Austin2g2.C                                                          
//                                                                           
// Created: Nov. 20, 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//
// Convert from Volumes in Austin format to volumes in g2 format
//                                                                           
//===========================================================================


#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include <vector>

using namespace Go;
using std::vector;

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc != 3, "Usage:  Input file, Output file" << std::endl);

    // Open input volume file
    std::ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Open output volume file
    std::ofstream os(argv[2]);
    ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

    int dim;  // Dimension
    is >> dim;
    ALWAYS_ERROR_IF(dim != 3, "Error in dimension");

    int ik1, ik2, ik3; // Order  in the 3 parameter directions
    is >> ik1;
    ik1++;
    is >> ik2;
    ik2++;
    is >> ik3;
    ik3++;
    int in1, in2, in3;  // Number of coefficients
    is >> in1;
    is >> in2;
    is >> in3;
    vector<double> st1(in1+ik1), st2(in2+ik2), st3(in3+ik3);  // Knot vectors
    
    int ki, kd;
    for (ki=0; ki<in1+ik1; ++ki)
	is >> st1[ki];
    
    for (ki=0; ki<in2+ik2; ++ki)
	is >> st2[ki];
    
    for (ki=0; ki<in3+ik3; ++ki)
	is >> st3[ki];
    
    vector<double> rc(in1*in2*in3*(dim+1));
    for (ki=0; ki<in1*in2*in3; ++ ki)
    {
	for (kd=0; kd<4; ++kd)
	    is >> rc[ki*(dim+1)+kd];
	for (kd=0; kd<3; ++kd)
	    rc[ki*(dim+1)+kd] *= rc[ki*(dim+1)+3];
    }

    SplineVolume *vol = new SplineVolume(in1, in2, in3, ik1, ik2, ik3, &st1[0],
					&st2[0], &st3[0], &rc[0], dim, true);

    vol->writeStandardHeader(os);
    vol->write(os);

    delete vol;
}
