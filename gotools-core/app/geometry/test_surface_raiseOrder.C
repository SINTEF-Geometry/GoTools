//==========================================================================
//                                                                          
// File: test_surface_raiseOrder.C                                           
//                                                                          
// Created: Mon Feb 13 13:45:04 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: test_surface_raiseOrder.C,v 1.2 2006-02-13 14:11:29 sbr Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char* argv[])
{
    if (argc != 5) {
	std::cout << "Usage:  testraiseorder filein raise_u raise_v fileout"
		  << std::endl;
	return -1;
    }

    SplineSurface the_surface, raised_surface;
    ObjectHeader header;

    ifstream infile(argv[1]);
    int raise_u(atoi(argv[2]));
    int raise_v(atoi(argv[3]));
    ofstream outfile(argv[4]);
    infile >> header >> the_surface;

    raised_surface = the_surface;
    raised_surface.raiseOrder(raise_u, raise_v);

    raised_surface.writeStandardHeader(outfile);
    raised_surface.write(outfile);

    return 0;

}
