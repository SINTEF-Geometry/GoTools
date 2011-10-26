//===========================================================================
//                                                                           
// File: readGoWriteIGES.C
//                                                                           
// Created: Tue April 3 2001
//                                                                           
// Author: Sverre Briseid
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/igeslib/IGESconverter.h"
#include <fstream>
#include <stdlib.h>  // For atof()
#include <iostream>

int main( int argc, char* argv[] )
{
    if (argc != 3) {
	std::cout << "usage: infile outfile" << std::endl;
	return -1;
    }

    std::ifstream infile(argv[1]);
    if (infile.bad()) {
	std::cout << "Infile not found or file corrupt" << std::endl;
	return -1;
    }

    IGESconverter conv1;
    conv1.readgo(infile);

    std::ofstream outfile(argv[2]);
    conv1.writeIGES(outfile);

}

















