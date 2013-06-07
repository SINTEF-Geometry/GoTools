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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc < 3, "Usage: " << argv[0]
		    << " inputsurf inputpoints" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read surface from file
    SplineSurface sf;
    is >> sf;

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*2);
    for (int i = 0; i < n; ++i) {
	pts >> pt[2*i] >> pt[2*i+1];
    }

    vector<double>::const_iterator coefs = sf.coefs_begin();
    int kk1 = sf.order_u();
    int kk2 = sf.order_v();
    int nn1 = sf.numCoefs_u();
    int dim = sf.dimension();
    std::vector<Point> p(3, Point(dim));
    std::vector<Point> p2(3, Point(dim));
    int kp, kq, kr, kd;
    for (int j = 0; j < n; ++j) {
	sf.point(p, pt[2*j], pt[2*j+1], 1);
	cout << p[0][0] << " " << p[0][1] << " " << p[0][2] << " " ;
	cout << p[1][0] << " " << p[1][1] << " " << p[1][2] << " " ;
	cout << p[2][0] << " " << p[2][1] << " " << p[2][2] << endl;

	vector<double> b0, bu, bv;
	sf.computeBasis(&pt[2*j], b0, bu, bv);
	int left1 = sf.basis_u().knotInterval(pt[2*j]) - kk1 + 1;
	int left2 = sf.basis_v().knotInterval(pt[2*j+1]) - kk2 + 1;
	
	p2[0].setValue(0.0);
	p2[1].setValue(0.0);
	p2[2].setValue(0.0);
	for (kq=left2, kr=0; kq<left2+kk2; kq++)
	    for (kp=left1; kp<left1+kk1; kp++, kr++)
		for (kd=0; kd<dim; kd++)
		{
		    double val = coefs[(kq*nn1+kp)*dim+kd];
		    p2[0].begin()[kd] += val*b0[kr];
		    p2[1].begin()[kd] += val*bu[kr];
		    p2[2].begin()[kd] += val*bv[kr];
		}
	cout << p2[0][0] << " " << p2[0][1] << " " << p2[0][2] << " " ;
	cout << p2[1][0] << " " << p2[1][1] << " " << p2[1][2] << " " ;
	cout << p2[2][0] << " " << p2[2][1] << " " << p2[2][2] << endl;
	cout << endl;
    }
}





