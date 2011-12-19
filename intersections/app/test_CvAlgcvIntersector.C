//===========================================================================
//                                                                           
// File: test_CvAlgcvIntersector.C                                           
//                                                                           
// Created: Tue Jan 25 15:42:35 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: test_CvAlgcvIntersector.C,v 1.7 2006-03-13 17:01:22 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
