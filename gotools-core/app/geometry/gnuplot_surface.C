//===========================================================================
//                                                                           
// File: gnuplot_surface.C                                                   
//                                                                           
// Created: Tue Jun 26 10:41:07 2001                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: gnuplot_surface.C,v 1.6 2005-06-09 07:15:40 oan Exp $
//                                                                           
// Description: Reads a 3D SplineSurface and writes out a file that may be
//              plotted in Gnuplot.
//                                                                           
//===========================================================================


#include "GoTools/utils/Array.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace Go;


int main(int argc, char** argv)
{
    // Read the curve from file
    std::ifstream input(argv[1]);
    if (input.bad()) {
	std::cerr << "File error (no file or corrupt file specified)."
		  << std::endl;
	return 1;
    }
    ObjectHeader header;
    SplineSurface surface;
    input >> header >> surface;

    // Loop through parameter space
    const int samples = 50;
    double increment_u = (surface.endparam_u()
			  - surface.startparam_u()) / (samples-1);
    double increment_v = (surface.endparam_v()
			  - surface.startparam_v()) / (samples-1);
    Point result;
    double param_u = surface.startparam_u();
    int prec = (int)std::cout.precision(15);
    for (int i = 0; i < samples; ++i) {
	double param_v = surface.startparam_v();
	for (int j = 0; j < samples; ++j) {
	    surface.point(result, param_u, param_v);
	    std::cout << result[0] << "\t" << result[1] << "\t"
		      << result[2] << std::endl;
	    param_v += increment_v;
	}
	std::cout << std::endl;
	param_u += increment_u;
    }
    std::cout.precision(prec);

    return 0;
}
