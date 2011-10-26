//===========================================================================
//                                                                           
// File: test_SfAlgObjIntersector.C                                          
//                                                                           
// Created: Wed Apr  6 17:22:27 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: test_SfAlgObjIntersector.C,v 1.4 2005-10-20 12:02:10 sbr Exp $
//                                                                           
// Description: Quite similar to test_SfAlgsfIntersector
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
#include "GoTools/geometry/sisl_file_io.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"
#include "GoTools/intersections/SphereInt.h"
#include "GoTools/intersections/CylinderInt.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


using std::cout;
using std::endl;
using std::ifstream;
using std::cerr;
using std::string;
using namespace Go;
using std::shared_ptr;
using std::dynamic_pointer_cast;


int main(int argc, char** argv)
{
 
    // The file_in should consist of
    //
    // spline_sf_file alg_type <int> tol
    // start alg_elem
    // ....
    // end alg_elem
    //
    // spline_sf  ...
    // ...

    if (argc != 3) {
	cout << "Usage: file_in file_out" << endl;
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

    // We continue until we reach end_of_file
    string line;
    string character;
    string app_name; // Not to be used by us.
    string sf_name;
    int alg_obj_type;
    while (filein) {
	//	while (filein && !
	getline(filein, line); //);
	// If pointer was at end of line we must read again.
	if (line.size() == 0)
	    getline(filein, line);
	// We push it into an istream for easy handling.
	std::istringstream string_stream(line); //, std::istringstream::in);
	shared_ptr<AlgObj3DInt> alg_obj_3d_int;
	// First element should be location of a spline sf.
	string_stream >> character >> app_name >> sf_name >> alg_obj_type;
	// We do not bother to read the rest of the line, moving on to the next.
	// We read the sisl-sf from file.
	FILE* fp = fopen(sf_name.c_str(), "r");
	SISLObject* wo;
	int kstat = 0;
	file_to_obj(fp, &wo, &kstat);

	SISLSurf* sisl_sf = wo->s1;
	shared_ptr<SplineSurface> spline_sf(SISLSurf2Go(sisl_sf));
	// We then treat the alg obj.
	shared_ptr<SplineSurface> alg_obj_sf_repr; // May be subsurface of object.
	if (alg_obj_type == 0) {
	    shared_ptr<PlaneInt> plane_int(new PlaneInt());
	    plane_int->read(filein);
	    alg_obj_3d_int = plane_int;
	    BoundingBox bd_box_sf = spline_sf->boundingBox();
	    // We pick the middle pt on diagonal, and project it onto the plane.
	    Point mid_pt = 0.5*(bd_box_sf.low() + bd_box_sf.high());
	    double length_x = 1.0*(bd_box_sf.low().dist(bd_box_sf.high()));
	    double length_y = length_x;
	    alg_obj_sf_repr = shared_ptr<SplineSurface>(plane_int->surface(mid_pt, length_x, length_y));
	} else if (alg_obj_type == 1) {
	    shared_ptr<CylinderInt> cylinder_int(new CylinderInt());
	    cylinder_int->read(filein);
	    alg_obj_3d_int = cylinder_int;
	    BoundingBox bd_box_sf = spline_sf->boundingBox();
	    double diag = bd_box_sf.low().dist(bd_box_sf.high());
	    // We need to tell the cylinder what subpart to visualize.
	    Point axis_unit = cylinder_int->ax_dir();
	    axis_unit.normalize();
	    Point bottom_pos = cylinder_int->ax_pt() - diag*axis_unit;
	    alg_obj_sf_repr = shared_ptr<SplineSurface>(cylinder_int->surface(bottom_pos, 2*diag));
	} else if (alg_obj_type == 2) {
	    shared_ptr<SphereInt> sphere_int(new SphereInt());
	    sphere_int->read(filein);
	    alg_obj_3d_int = sphere_int;
	    alg_obj_sf_repr = shared_ptr<SplineSurface>(sphere_int->surface());
	} else {
	    std::cout << "Unexpected alg obj type, exiting!" << std::endl;
	    return -1;
	}

	// We intersect the sf with the algebraic object.
	shared_ptr<ParamSurface> par_sf = spline_sf;
	shared_ptr<ParamObjectInt> par_obj_int(new ParamSurfaceInt(par_sf));
	shared_ptr<AlgObjectInt> alg_obj_int = alg_obj_3d_int;

	shared_ptr<GeoTol> aepsge(new GeoTol(1.e-06)); //8));
	IntersectorAlgPar int_alg_par(alg_obj_int, par_obj_int, aepsge);

	// We write to file (u, v, f(u,v)), where f is the sf plugged into the alg obj.

	std::vector<shared_ptr<IntersectionPoint> > int_points;
	std::vector<shared_ptr<IntersectionCurve> > int_curves;
	try {
	    int_alg_par.compute();
	    int_alg_par.getResult(int_points, int_curves);
	} catch (...) {
	    MESSAGE("Failed computing intersections. Moving on to next case.");
	}
	cout << "IntPoints found: " << int_points.size() << endl;
	cout << "IntCurves found: " << int_curves.size() << endl;

#ifdef INTERSECTIONS_DEBUG
	std::ofstream debug("data/debug.g2");
	for (size_t ki = 0; ki < int_curves.size(); ++ki) {
	    int obj_nmb = 1;
	    shared_ptr<SplineCurve> int_par_cv =
		dynamic_pointer_cast<SplineCurve, ParamCurve>
		(int_curves[ki]->getParamCurve(obj_nmb));
	    if (int_par_cv.get() != 0) {
		// The Plot does only make sense for intersection with plane parallell to xy-plane.
		// Could of course translate and rotate to 
		writeSpaceParamCurve(*int_par_cv, debug); //, -d);
	    }
	    shared_ptr<SplineCurve> space_cv;
	    try {
		space_cv =
		    dynamic_pointer_cast<SplineCurve, ParamCurve>(int_curves[ki]->getCurve());
	    } catch (...) {
		MESSAGE("Failed generating space representation of intersection curve!");
	    }
	    if (space_cv.get()) {
		space_cv->writeStandardHeader(debug);
		space_cv->write(debug);
	    }
	    // 	if (int_space_cv.get() != 0) {
	    // 	    int_space_cv->writeStandardHeader(debug);
	    // 	    int_space_cv->write(debug);
	    // 	}
	}
#endif // INTERSECTIONS_DEBUG

}

    return 0;
}
