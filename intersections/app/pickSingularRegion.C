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
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <memory>
#include <time.h>
#include <string>
#include <fstream>
#include <iomanip>


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
 
    if (argc != 3) {
	cout << "Usage: pickSingularRegion FileSf1 FileSf2"
	     << endl;
	return 0;
    }


    ObjectHeader header;

    // Read the first curve from file
    ifstream input1(argv[1]);
    if (input1.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input1);
    shared_ptr<ParamSurface> surf1(new SplineSurface());
    surf1->read(input1);
    input1.close();
    
    // Read the second curve from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input2);
    shared_ptr<ParamSurface> surf2(new SplineSurface());
    surf2->read(input2);
    input2.close();

    cout << setprecision(15);

    // Make the first subsurface
    SplineSurface* sf1 = dynamic_cast<SplineSurface*>(surf1.get());
    double from_upar = 52.0248;
    double from_vpar = 0.968029;
    double to_upar = sf1->endparam_u();
    double to_vpar = 1.03197;
    SplineSurface* sub1
	= sf1->subSurface(from_upar, from_vpar, to_upar, to_vpar);
    // Make the second subsurface
    SplineSurface* sf2 = dynamic_cast<SplineSurface*>(surf2.get());
    from_upar = 52.03;
    from_vpar = 0.968027;
    to_upar = 52.0349;
    to_vpar = 1.03197;
    SplineSurface* sub2
	= sf2->subSurface(from_upar, from_vpar, to_upar, to_vpar);

    // Magnify in the direction normal to the surfaces
    typedef vector<double>::iterator iter;
    for (iter it = sub1->coefs_begin(); it != sub1->coefs_end(); it += 3) {
	it[0] -= 5.5;
	it[1] += 10.6;
	it[2] += 21.9;
	it[0] *= 100;
	it[1] *= 100;
	it[2] *= 100;
    }
    for (iter it = sub2->coefs_begin(); it != sub2->coefs_end(); it += 3) {
	it[0] -= 5.5;
	it[1] += 10.6;
	it[2] += 21.9;
	it[0] *= 100;
	it[1] *= 100;
	it[2] *= 100;
    }
    Point normal1 = sub1->normalCone().centre();
    Point normal2 = sub2->normalCone().centre();
    Point normal = normal1 + normal2;
    for (iter it = sub1->coefs_begin(); it != sub1->coefs_end(); it += 3) {
	Point coef(it[0], it[1], it[2]);
	double fac = normal * coef;
	fac *= 100.0;
	it[0] += fac * normal[0];
	it[1] += fac * normal[1];
	it[2] += fac * normal[2];
    }
    for (iter it = sub2->coefs_begin(); it != sub2->coefs_end(); it += 3) {
	Point coef(it[0], it[1], it[2]);
	double fac = normal * coef;
	fac *= 100.0;
	it[0] += fac * normal[0];
	it[1] += fac * normal[1];
	it[2] += fac * normal[2];
    }

    ofstream out1("subsurf1.g2");
    sub1->writeStandardHeader(out1);
    sub1->write(out1);
    ofstream out2("subsurf2.g2");
    sub2->writeStandardHeader(out2);
    sub2->write(out2);


    return 0;
}
