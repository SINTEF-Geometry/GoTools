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
#include "GoTools/intersections/Par1FuncIntersector.h"
#include "GoTools/intersections/Param1FunctionInt.h"
#include "GoTools/intersections/Spline1FunctionInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/AlgObjectInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/IntersectorAlgPar.h"
#include "GoTools/intersections/Line2DInt.h"
#include "GoTools/intersections/GeoTol.h"
#include "GoTools/geometry/SplineCurve.h"
#include <memory>
#include <fstream>
#include <iomanip>


using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using namespace Go;


int main(int argc, char** argv)
{
 
    if (argc != 5) {
	// Our input curve should lie in the parameter plane.
	// The line is given by: ax + by + c = 0
	cout << "Usage: 2d_spline_cv a b c" << endl;
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
    double a = atof(argv[2]);
    double b = atof(argv[3]);
    double c = atof(argv[4]);

    header.read(filein);
    shared_ptr<ParamCurve> curve(new SplineCurve());
    curve->read(filein);
    filein.close();
    
    shared_ptr<Line2DInt> line(new Line2DInt(a, b, c));
//     double dummy_tol = 1e-03; // @@sbr Not yet used.
    shared_ptr<ParamObjectInt> par_obj_int(new ParamCurveInt(curve));
    shared_ptr<AlgObjectInt> alg_obj_int = line;

    shared_ptr<GeoTol> aepsge(new GeoTol(1.e-6));
    IntersectorAlgPar int_alg_par(alg_obj_int, par_obj_int, aepsge);

    int_alg_par.compute();
    std::vector<shared_ptr<IntersectionPoint> > int_points;
    std::vector<shared_ptr<IntersectionCurve> > int_curves;
    int_alg_par.getResult(int_points, int_curves);

    cout << "IntPoints found: " << int_points.size() << endl;
    cout << "IntCurves found: " << int_curves.size() << endl;

    return 0;
}
