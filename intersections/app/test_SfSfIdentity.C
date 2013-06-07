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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/igeslib/IGESconverter.h"
#include <fstream>

using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using namespace Go;

int main(int argc, char** argv)
{
    if (argc != 4) {
	cout << "Usage: test_SfSfIdentity FileSf1 FileSf2 "
	     <<" aepsge"
	     << endl;
	return 0;
    }
    ObjectHeader header;

    // Read the first surface from file
    std::ifstream input1(argv[1]);
    if (input1.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    IGESconverter conv1;
    conv1.readgo(input1);
    vector<shared_ptr<GeomObject> > geomobj1 = conv1.getGoGeom();
    shared_ptr<ParamSurface> surf1;
    if (geomobj1[0]->instanceType() == Class_SplineSurface || 
	geomobj1[0]->instanceType() == Class_BoundedSurface)
	surf1 = dynamic_pointer_cast<ParamSurface, GeomObject>(geomobj1[0]);

    
    // Read the second surface from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }

    IGESconverter conv2;
    conv2.readgo(input2);
    vector<shared_ptr<GeomObject> > geomobj2 = conv2.getGoGeom();
    shared_ptr<ParamSurface> surf2;
    if (geomobj2[0]->instanceType() == Class_SplineSurface || 
	geomobj2[0]->instanceType() == Class_BoundedSurface)
	surf2 = dynamic_pointer_cast<ParamSurface, GeomObject>(geomobj2[0]);

    double aepsge;
    aepsge = atof(argv[3]);

    Identity ident;
    int coincident = ident.identicalSfs(surf1, surf2, aepsge);

    std::cout << "Identity: " << coincident << std::endl;
}
