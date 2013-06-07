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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/errormacros.h"

using namespace std;

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc < 3, "Usage: " << argv[0]
		<< " inputsurf inputpoints" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read surface from file
    int in1, in2, ik1, ik2, dim;
    std::vector<double> et1, et2, co;
    is >> dim;
    is >> in1 >> ik1;
    et1.resize(in1+ik1);
    for (int i = 0; i < in1+ik1; ++i)
	is >> et1[i];
    is >> in2 >> ik2;
    et2.resize(in2+ik2);
    for (int i = 0; i < in2+ik2; ++i)
	is >> et2[i];
    co.resize(in1*in2*dim);
    for (int i = 0; i < in1*in2*dim; ++i)
	is >> co[i];
    Go::BsplineBasis bas1(in1, ik1, et1.begin());
    Go::BsplineBasis bas2(in2, ik2, et2.begin());

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    // double u, v;
    int n;
    pts >> n;
    vector<double> pt(n*2);
    for (int i = 0; i < n; ++i) {
	pts >> pt[2*i] >> pt[2*i+1];
    }

    double b1[8];
    double b2[8];
    // int il1, il2;
    for (int j = 0; j < n; ++j) {
	//  	    il1 = bas1.knotInterval(pt[2*j]);
	//  	    il2 = bas2.knotInterval(pt[2*j+1]);
	bas1.computeBasisValues(pt[2*j], b1, 1);
	bas2.computeBasisValues(pt[2*j+1], b2, 1);
	for (int ii = 0; ii < ik1; ++ii) {
	    cout << b1[2*ii] << ' ' << b1[2*ii+1] << endl;
	}
    }
    double kkk[] = { 0, 0, 0, 0, 1, 2, 4, 7, 8, 9, 14, 14, 14, 14 };
    Go::BsplineBasis bas3(10, 4, kkk);
    cout << bas3;
    bas3.reverseParameterDirection();
    cout << bas3;
}





