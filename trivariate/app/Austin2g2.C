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
