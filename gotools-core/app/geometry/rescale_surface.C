/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/geometry/SplineSurface.h"
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
    SplineSurface surface;
    input >> header >> surface;

    // @@@ Only implemented for 3D and non-rational case
    if (surface.rational()) {
	std::cerr << "Surface must not be rational."
		  << std::endl;
	return 1;
    }
    if (surface.dimension() != 3) {
	std::cerr << "Dimension must be 3."
		  << std::endl;
	return 1;
    }


    // Find largest and smallest coordinates
    const int dim = 3;
    int numcoefs = surface.numCoefs_u() * surface.numCoefs_v();
    vector<double> min(dim, 0.0);
    vector<double> max(dim, 0.0);
    for (int i = 0; i < numcoefs; ++i) {
	for (int dd = 0; dd < dim; ++dd) {
	    min[dd] = std::min(min[dd], surface.coefs_begin()[i*dim+dd]);
	    max[dd] = std::max(max[dd], surface.coefs_begin()[i*dim+dd]);
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
	    double& temp = surface.coefs_begin()[i*dim+dd];
	    temp = (temp - shift[dd]) / l;
	}
    }

    // Write out rescaled surface
    std::cout << header << surface;

    return 0;
}

