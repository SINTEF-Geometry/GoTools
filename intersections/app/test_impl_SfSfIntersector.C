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
#include "GoTools/intersections/PlaneInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/IntersectorAlgPar.h"
#include <memory>
#include "GoTools/geometry/sisl_file_io.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"
#include "GoTools/intersections/SphereInt.h"
#include "GoTools/intersections/CylinderInt.h"
#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/IntersectionUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineSurface.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


using std::cout;
using std::endl;
using std::ifstream;
using std::cerr;
using std::ofstream;
using std::vector;
using namespace Go;


int main(int argc, char** argv)
{
 
    // The file_in should consist of
    //
    // spline_sf_1 spline_sf2
    //
    // The first object is 
    //

    if (argc != 6) {
	cout << "Usage: file_in impl_first_obj impl_tol file_out multiply_factor" << endl;
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
    bool impl_first_obj(atoi(argv[2])!=0);
    double impl_tol(atof(argv[3]));
    ofstream fileout(argv[4]);
    double fmult(atof(argv[5]));

    IGESconverter conv;
    conv.readgo(filein);
    vector<shared_ptr<GeomObject> > obj = conv.getGoGeom();
    vector<shared_ptr<SplineSurface> > spline_sfs;
    for (size_t ki = 0; ki < obj.size(); ++ki) {
	if (obj[ki]->instanceType() == Class_SplineSurface) {
	    spline_sfs.push_back(dynamic_pointer_cast<SplineSurface, GeomObject>(obj[ki]));
	} else {
	    MESSAGE("Expecting SplineSurface's only, neglecting object.");
	}
    }
    ASSERT(spline_sfs.size() == 2);

    // If both sfs are rational we quit.
    if (spline_sfs[0]->rational() && spline_sfs[1]->rational()) {
	cout << "Both sfs rational, currently not supported." << endl;
	return -1;
    }
    // If the sf to insert is rational (but not the other) we also quit,
    // suggesting that user could try the other way around.
    if ((impl_first_obj && spline_sfs[1]->rational()) ||
	(!impl_first_obj && spline_sfs[0]->rational())) {
	cout << "The sf to be inserted is rational, currently not supported. " <<
	    "Why not try impliziticing the other object?" << endl;
	cout << "But, why not try anyway, let's see what happens!" << std::endl;
// 	return -1;
    }

    BoundingBox box1 = spline_sfs[0]->boundingBox();
    BoundingBox box2 = spline_sfs[1]->boundingBox();
    Point mid1 = 0.5*(box1.low() + box1.high());
    Point mid2 = 0.5*(box2.low() + box2.high());
    Point mid = 0.5*(mid1 + mid2);

    // Translate and scale spline surfaces
    for (int kr=0; kr<2; kr++)
    {
    int n1 = spline_sfs[kr]->numCoefs_u();
    int n2 = spline_sfs[kr]->numCoefs_v();
    int dim = spline_sfs[kr]->dimension();
    vector<double>::iterator c1 = spline_sfs[kr]->coefs_begin();
    int ki, kj;
    for (ki=0; ki<n1*n2; ki++)
    {
	for (kj=0; kj<dim; kj++, c1++)
	{
	    *c1 -= mid[kj];
	    *c1 *= fmult;
	}
    }
    }

    SplineSurfaceInt sf_int1(spline_sfs[0]);
    SplineSurfaceInt sf_int2(spline_sfs[1]);

    bool can_impl = impl_first_obj ? sf_int1.canImplicitize() : sf_int2.canImplicitize();
    if (!can_impl) {
	MESSAGE("Input sf not suited for implicitization! Maybe ...");
	// 		return -1;
    }

    SplineSurfaceInt spline_sf_int = impl_first_obj ? sf_int2 : sf_int1;
    SplineSurfaceInt spline_sf_other = impl_first_obj ? sf_int1 : sf_int2;
    double tol2;
    AlgObj3DInt alg_obj_int(-1);
    if (impl_first_obj) {
	bool impl_status = sf_int1.implicitize(impl_tol);
	if (impl_status != true) {
	    MESSAGE("Failed implicitizing!");
	}
	impl_status = sf_int1.getImplicit(impl_tol, tol2, alg_obj_int);
	if (impl_status != true) {
	    MESSAGE("Failed implicitizing!");
	}	
    } else {
	bool impl_status = sf_int2.implicitize(impl_tol);
	if (impl_status != true) {
	    MESSAGE("Failed implicitizing!");
	}
	impl_status = sf_int2.getImplicit(impl_tol, tol2, alg_obj_int);
	if (impl_status != true) {
	    MESSAGE("Failed implicitizing!");
	}
    }

    BernsteinTetrahedralPoly impl;
    BaryCoordSystem3D bc;
    alg_obj_int.getImplicit(impl, bc);
    shared_ptr<const SplineSurface> spline_sf = spline_sf_int.splineSurface();
    shared_ptr<SplineSurface> spline_sf_1d =
	IntersectionUtils::insertSfInImplObj(*spline_sf, impl, bc);

    // TESTING
    shared_ptr<const SplineSurface> spline_sf2 = spline_sf_other.splineSurface();
    shared_ptr<SplineSurface> spline_sf_1d_2 =
	IntersectionUtils::insertSfInImplObj(*spline_sf2, impl, bc);
    
    // spline_sf_1d is our composition surface: impl(spline_sf):R^2->R
    // We check that the composition is correct.
    bool check_comp_sf = false; //true;
    while (check_comp_sf) {
	double ufrac = 0.5; // 0.0 <= ufrac <= 1.0
	double vfrac = 0.5; // 0.0 <= vfrac <= 1.0
	double upar = ufrac*spline_sf->startparam_u() + (1 - ufrac)*spline_sf->endparam_u();
	double vpar = vfrac*spline_sf->startparam_v() + (1 - vfrac)*spline_sf->endparam_v();
	double dist = IntersectionUtils::distImplRepresentationCompFunction
	    (*spline_sf, impl, bc, *spline_sf_1d, upar, vpar);
	std::cout << "dist: " << dist << std::endl;
    }

    // We then intersect the sf with 0.0.
    shared_ptr<Param2FunctionInt>
	par_func_int(new Spline2FunctionInt(spline_sf_1d));
    shared_ptr<Param0FunctionInt> const_int(new Param0FunctionInt(0.0));

    shared_ptr<GeoTol> aepsge(new GeoTol(1.e-06)); //8));
    Par2FuncIntersector par2_func_int(par_func_int, const_int, aepsge);

    std::vector<shared_ptr<IntersectionPoint> > int_points;
    std::vector<shared_ptr<IntersectionCurve> > int_curves;
    try {
	par2_func_int.compute();
	par2_func_int.getResult(int_points, int_curves);
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
	    SplineDebugUtils::writeSpaceParamCurve(*int_par_cv, debug); //, -d);
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
    }
#endif // INTERSECTIONS_DEBUG

    return 0;
}
