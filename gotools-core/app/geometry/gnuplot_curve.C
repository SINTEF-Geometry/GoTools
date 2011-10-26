//===========================================================================
//                                                                           
// File: gnuplot_curve.C                                                     
//                                                                           
// Created: Mon Jun 25 16:45:33 2001                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: gnuplot_curve.C,v 1.6 2005-06-09 07:15:39 oan Exp $
//                                                                           
// Description: Reads a 2D SplineCurve and writes out a file that may be
//              plotted in Gnuplot. 
//                                                                           
//===========================================================================


#include "GoTools/utils/Array.h"
#include "GoTools/geometry/SplineCurve.h"
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
    SplineCurve curve;
    input >> header >> curve;

    // Loop through parameter space
    const int samples = 100;
    double param = curve.startparam();
    double increment = (curve.endparam() - curve.startparam()) / (samples-1);
    Point result;
    int prec = (int)std::cout.precision(15);
    for (int i = 0; i < samples; ++i) {
	curve.point(result, param);
	std::cout << result[0] << "\t" << result[1] << std::endl;
	param += increment;
    }
    std::cout.precision(prec);

    return 0;
}
