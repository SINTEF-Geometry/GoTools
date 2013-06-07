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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc < 4, "Usage: " << argv[0]
		    << " inputvol inputpoints show_derivs?" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read surface from file
    SplineVolume vol;
    is >> vol;

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*3);
    for (int i = 0; i < n; ++i) {
	pts >> pt[3*i] >> pt[3*i+1] >> pt[3*i+2];
    }

    bool showDerivs = (atoi(argv[3]) != 0);

    vector<double>::const_iterator coefs = vol.coefs_begin();
    int kk1 = vol.order(0);
    int kk2 = vol.order(1);
    int kk3 = vol.order(2);
    int nn1 = vol.numCoefs(0);
    int nn2 = vol.numCoefs(1);
    int dim = vol.dimension();
    std::vector<Point> p(4, Point(dim));
    //Point p(dim);
    std::vector<Point> p2(4, Point(dim));
    int kp, kq, kh, kr, kd;
    for (int j = 0; j < n; ++j) {
 	vol.point(p, pt[3*j], pt[3*j+1], pt[3*j+2], 1);
 	cout << p[0][0] << " " << p[0][1] << " " << p[0][2] << " " ;
 	cout << p[1][0] << " " << p[1][1] << " " << p[1][2] << " " ;
	cout << p[2][0] << " " << p[2][1] << " " << p[2][2] << " " ;
 	cout << p[3][0] << " " << p[3][1] << " " << p[3][2] << endl;

	vector<double> b0, bu, bv, bw;
	vol.computeBasis(&pt[3*j], b0, bu, bv, bw);
	int left1 = vol.basis(0).knotInterval(pt[3*j]) - kk1 + 1;
	int left2 = vol.basis(1).knotInterval(pt[3*j+1]) - kk2 + 1;
	int left3 = vol.basis(2).knotInterval(pt[3*j+2]) - kk3 + 1;
	
	BasisPts res;
	vol.computeBasis(pt[3*j], pt[3*j+1], pt[3*j+2], res);

	p2[0].setValue(0.0);
	p2[1].setValue(0.0);
	p2[2].setValue(0.0);
	p2[3].setValue(0.0);
	for (kh=left3, kr=0; kh<left3+kk3; kh++)
	    for (kq=left2; kq<left2+kk2; kq++)
		for (kp=left1; kp<left1+kk1; kp++, kr++)
		    for (kd=0; kd<dim; kd++)
		    {
			double val = coefs[((kq+kh*nn2)*nn1+kp)*dim+kd];
			p2[0].begin()[kd] += val*b0[kr];
			p2[1].begin()[kd] += val*bu[kr];
			p2[2].begin()[kd] += val*bv[kr];
			p2[3].begin()[kd] += val*bw[kr];
		    }
	cout << p2[0][0] << " " << p2[0][1] << " " << p2[0][2] << " " ;
	if (showDerivs)
	  {
	    cout << p2[1][0] << " " << p2[1][1] << " " << p2[1][2] << " " ;
	    cout << p2[2][0] << " " << p2[2][1] << " " << p2[2][2] << " " ;
	    cout << p2[3][0] << " " << p2[3][1] << " " << p2[3][2];
	  }
	cout << endl;
	cout << endl;
    }

    ofstream os("tmp_vol.g2");
    vol.writeStandardHeader(os);
    vol.write(os);
}





