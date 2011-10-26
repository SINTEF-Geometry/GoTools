//===========================================================================
//                                                                           
// File: readIGESwriteGO.C
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

int main( int argc, char* argv[] )
{

  std::ifstream infile(argv[1]);
  if (infile.bad()) {
    std::cout << "Infile not found or file corrupt" << std::endl;
    return -1;
  }
  if (argc != 3) {
    std::cout << "Expecting 2 arguments (infile outfile)." << std::endl;
    return -1;
  }
  IGESconverter conv1;
  try {
      conv1.readIGES(infile);
  } catch (...) {
      std::cout << "Failed reading input IGES file, exiting." << std::endl;
      return -1;  
  }

  std::ofstream outfile(argv[2]);
  conv1.writego(outfile);

}
