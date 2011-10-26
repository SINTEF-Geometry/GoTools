//===========================================================================
//                                                                           
// File: test_Par2FuncIntersector.C                                          
//                                                                           
// Created: Wed Nov 17 08:59:05 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: test_SfAlgsfIntersector.C,v 1.4 2006-03-03 15:50:51 jbt Exp $
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
#include "GoTools/intersections/PlaneInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/IntersectorAlgPar.h"

#include <fstream>
#include <iomanip>


using std::cout;
using std::endl;
using std::cerr;
using std::ifstream;
using namespace Go;
using std::shared_ptr;


int main(int argc, char** argv)
{
 
    if (argc != 6) {
	cout << "Usage: spline_sf a b c d" << endl;
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
    double d = atof(argv[5]);

    header.read(filein);
    shared_ptr<ParamSurface> surf(new SplineSurface());
    surf->read(filein);
    filein.close();
    
    shared_ptr<PlaneInt> plane(new PlaneInt(a, b, c, d));
//     double dummy_tol = 1e-03;
    shared_ptr<ParamObjectInt> par_obj_int(new ParamSurfaceInt(surf));
    shared_ptr<AlgObjectInt> alg_obj_int = plane;

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
