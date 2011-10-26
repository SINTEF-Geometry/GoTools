//===========================================================================
//                                                                           
// File: rescale_curve.C                                                           
//                                                                           
// Created: Thu Jun 28 16:52:12 2001                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: rescale_curve.C,v 1.3 2003-05-08 15:37:15 afr Exp $
//                                                                           
// Description: Rescales curve to make it lie within the unit box
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>
#include <algorithm>


using namespace Go;
using namespace std;


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

    if (curve.rational()) {
	std::cerr << "Curve cannot be rational."
		  << std::endl;
	return 1;
    }
    if (curve.dimension() != 2) {
	std::cerr << "Dimension must be 2."
		  << std::endl;
	return 1;
    }


    // Find largest and smallest coordinates
    const int dim = 2;
    int numcoefs = curve.numCoefs();
    vector<double> min(dim, 0.0);
    vector<double> max(dim, 0.0);
    for (int i = 0; i < numcoefs; ++i) {
	for (int dd = 0; dd < dim; ++dd) {
	    min[dd] = std::min(min[dd], curve.coefs_begin()[i*dim+dd]);
	    max[dd] = std::max(max[dd], curve.coefs_begin()[i*dim+dd]);
	}
    }

    // Scaling factor and shift vector
    double l = 0.0;
    vector<double> shift(dim, 0.0);
    for (int dd = 0; dd < dim; ++dd) {
	l = std::max(l, max[dd] - min[dd]);
	shift[dd] = (max[dd] + min[dd]) / 2;
    }
    l *= 0.5;

    // Scale and shift the coefficients
    for (int i = 0; i < numcoefs; ++i) {
	for (int dd = 0; dd < dim; ++dd) {
	    double& temp = curve.coefs_begin()[i*dim+dd];
	    temp = (temp - shift[dd]) / l;
	}
    }

    // Write out rescaled curve
    std::cout << header << curve;

    return 0;
}

