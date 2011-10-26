//===========================================================================
//                                                                           
// File: test_SphereIntersector.C                                          
//                                                                           
// Created: Wed Nov 17 08:59:05 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision:
// $Id: test_SfSphereIntersector.C,v 1.6 2007-01-15 10:12:31 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
using std::shared_ptr;


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
