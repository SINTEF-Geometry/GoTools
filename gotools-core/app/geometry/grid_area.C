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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/RectGrid.h"

using namespace std;
using namespace Go;

int main(int argc, char** argv)
{
    if (argc < 2) {
	cerr << "Usage: " << argv[0] << " method (1 = quads, 2 = triangles)" << endl;
	return 1;
    }
    int method = atoi(argv[1]);

    ObjectHeader head;
    cin >> head;
    ASSERT(head.classType() == RectGrid::classType());
    RectGrid grid;
    cin >> grid;

    int numu = grid.numCoefs_u();
    int numv = grid.numCoefs_v();
    double area = 0.0;
    const double* p = grid.rawData();
    Vector3D q1;
    Vector3D q2;
    Vector3D q3;
    Vector3D q4;
    for (int i = 0; i < numu - 1; ++i) {
	for (int j = 0; j < numv - 1; ++j) {
	    Vector3D p0(p + 3*(i + j*numu));
	    Vector3D p1(p + 3*(i + j*numu + 1));
	    Vector3D p2(p + 3*(i + (j + 1)*numu + 1));
	    Vector3D p3(p + 3*(i + (j + 1)*numu));
	    switch(method) {
	    case 1:
		q1 = p2;
		q1 -= p0;
		q2 = p3;
		q2 -= p1;
		q3 = q1.cross(q2);
		area += 0.5*q3.length();
		break;
	    case 2:
		q1 = p1;
		q1 -= p0;
		q2 = p2;
		q2 -= p0;
		q3 = q1.cross(q2);
		area += 0.5*q3.length();
		q1 = p3;
		q1 -= p0;
		q3 = q2.cross(q1);
		area += 0.5*q3.length();
		break;
	    case 3:
		q1 = p1;
		q1 -= p0;
		q2 = p2;
		q2 -= p1;
		q3 = p3;
		q3 -= p2;
		q4 = p0;
		q4 -= p3;
		area += 0.25*(q1.cross(q4)).length();
		area += 0.25*(q2.cross(q3)).length();
		area += 0.25*(q1.cross(q2)).length();
		area += 0.25*(q3.cross(q4)).length();
		break;
	    case 4:
		q1 = p1;
		q1 -= p0;
		q2 = p2;
		q2 -= p1;
		q3 = p3;
		q3 -= p2;
		q4 = p0;
		q4 -= p3;
		area += 0.5*(q1.cross(q4)).length();
		area += 0.5*(q2.cross(q3)).length();
		//area += 0.25*(q1.cross(q2)).length();
		//area += 0.25*(q3.cross(q4)).length();
		break;
	    default:
		THROW("No such method: " << method);
	    }
	}
    }
    cout.precision(18);
    cout << area << endl;
}
