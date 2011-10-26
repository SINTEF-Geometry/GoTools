//===========================================================================
//                                                                           
// File: testknotinsertion.C                                                 
//                                                                           
// Created: Thu Apr  5 12:52:03 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: testknotinsertion.C,v 1.3 2003-05-08 15:37:15 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>
#include <sstream>

using namespace Go;

int main(int argc, char** argv)
{
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    if (input.bad() || output.bad()) {
	std::cerr << "File error (no file or corrupt file specified)."
		  << std::endl;
	return 1;
    }
    ObjectHeader header;
    SplineSurface surf;
    input >> header >> surf;

    std::vector<double> knotsu;
    std::cout << "Enter knots to insert (u): ";
    std::string st;
    std::getline(std::cin, st);
    std::istringstream is1(st);
    while(is1) { double x; is1 >> x; knotsu.push_back(x); }
    knotsu.pop_back();

    std::vector<double> knotsv;
    std::cout << "Enter knots to insert (v): ";
    std::getline(std::cin, st);
    std::istringstream is2(st);
    while(is2) { double x; is2 >> x; knotsv.push_back(x); }
    knotsv.pop_back();

    surf.insertKnot_u(knotsu);
    surf.insertKnot_v(knotsv);

    output << header << surf;
}
