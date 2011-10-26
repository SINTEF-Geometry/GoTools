//===========================================================================
//                                                                           
// File: readDispWriteGo.C
//                                                                           
// Created: Tue April 3 2001
//                                                                           
// Author: Sverre Briseid
//                                                                           
// Revision: 
//                                                                           
// Description: Reads geometry on the Disp format from file, returns Go geometry.
//                                                                           
//===========================================================================


#include "GoTools/igeslib/IGESconverter.h"
#include <fstream>
#include <stdlib.h>  // For atof()

int main( int argc, char* argv[] )
{
    if (argc != 3) {
	MESSAGE("Usage: infile outfile.");
	return 1;
    }

    std::ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(),
		    "Infile not found or file corrupt");
    std::ofstream outfile(argv[2]);

    IGESconverter conv1;
    conv1.readdisp(infile);
    conv1.writego(outfile);
}

















