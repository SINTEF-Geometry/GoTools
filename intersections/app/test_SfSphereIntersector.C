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
#include "GoTools/intersections/Par2FuncIntersector.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/intersections/Spline2FunctionInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/intersections/SphereInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/IntersectorAlgPar.h"
#include <memory>
#include <fstream>
#include <iomanip>


using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;
using namespace Go;


int main(int argc, char** argv)
{
 
    if (argc != 6) {
	cout << "Usage: spline_sf x0 y0 z0 r" << endl;
	return 0;
    }

    ObjectHeader header;

    // Read the first curve from file
    ifstream filein(argv[1]);
    if (filein.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    double x0 = atof(argv[2]);
    double y0 = atof(argv[3]);
    double z0 = atof(argv[4]);
    double r = atof(argv[5]);

    header.read(filein);
    shared_ptr<ParamSurface> surf(new SplineSurface());
    surf->read(filein);
    filein.close();

    Point center(x0, y0, z0);    
    shared_ptr<SphereInt> sphere(new SphereInt(center, r));
//     double dummy_tol = 1e-03; // @@sbr Not yet used.
    shared_ptr<ParamObjectInt> par_obj_int(new ParamSurfaceInt(surf));
    shared_ptr<AlgObjectInt> alg_obj_int = sphere;

    shared_ptr<GeoTol> aepsge(new GeoTol(1.e-6));
    //IntersectorAlgPar int_alg_par(alg_obj_int, par_obj_int, aepsge);
    IntersectorAlgPar int_alg_par(sphere, par_obj_int, aepsge);

    int_alg_par.compute();
    std::vector<shared_ptr<IntersectionPoint> > int_points;
    std::vector<shared_ptr<IntersectionCurve> > int_curves;
    int_alg_par.getResult(int_points, int_curves);

    cout << "IntPoints found: " << int_points.size() << endl;
    cout << "IntCurves found: " << int_curves.size() << endl;

    return 0;
}
