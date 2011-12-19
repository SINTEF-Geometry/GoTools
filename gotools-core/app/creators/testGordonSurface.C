//===========================================================================
//                                                                           
// File: GordonTest.C                                                        
//                                                                           
// Created: Mon Jul 23 16:57:39 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: testGordonSurface.C,v 1.6 2008-09-16 09:07:35 jbt Exp $
//                                                                           
// Description: Reads spline surface, extracts specified iso-curves, makes
//              Gordon surface from curves, writes Gordon surface and curves
//              to file.
//===========================================================================

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;

using namespace Go;

int main(int argc, char* argv[])
{
    ALWAYS_ERROR_IF(argc < 6,
		    "Usage: surfin nmb_u_par par_u nmb_v_par par_v surfout");

    ObjectHeader header_sf;
    SplineSurface input_surface;
    ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(), "Infile not found or file corrupt");
    infile >> header_sf >> input_surface;

    SplineCurve curr_curve;
    vector<shared_ptr<SplineCurve> > mesh_curves;

    int nmb_u_crvs = atoi(argv[2]);
    int nmb_v_crvs = atoi(argv[2+nmb_u_crvs+1]);
    if (nmb_u_crvs + nmb_v_crvs + 5 != argc)
    {
	puts("Usage: surfin nmb_u_par par_u nmb_v_par par_v surfout");
	return -1;
    }
    double tpar;
    //    cout << "Specify v-parameters for u-curves: ";
    vector<double> params;
    for (int i = 3; i < nmb_u_crvs + 3; ++i) {
	tpar = atof(argv[i]);
	params.push_back(tpar);
	mesh_curves.push_back
	    (shared_ptr<SplineCurve>(input_surface.constParamCurve(tpar, true)));
    }
    // We now extract the v-iso-curves.
    for (int i = nmb_u_crvs + 4; i < argc - 1; ++i) {
	tpar = atof(argv[i]);
	params.push_back(tpar);
	mesh_curves.push_back
	    (shared_ptr<SplineCurve>(input_surface.constParamCurve(tpar, false)));
    }

    vector<shared_ptr<SplineCurve> > cp_mesh_curves = mesh_curves;

    SplineSurface* pGordonSurface =
	CoonsPatchGen::createGordonSurface(mesh_curves, params,
					   nmb_u_crvs, true);

//     ObjectHeader header_crv;
//     // Assumes presence of CurveSpline header.
//     ifstream infile2("data/header_crv.txt");
//     infile2 >> header_crv;
    // Same header applies.
    ofstream outfile(argv[argc - 1]);
    outfile << header_sf << *pGordonSurface;
    for (int i = 0; i < (int)cp_mesh_curves.size(); ++i)
    {
	cp_mesh_curves[i]->writeStandardHeader(outfile);
	cp_mesh_curves[i]->write(outfile);
// 	outfile << header_crv << *(cp_mesh_curves[i]);
    }

    ofstream outfile2("face0.srf");
    vector<double> new_knots;
    new_knots.push_back(0.2);
    new_knots.push_back(0.5);
    new_knots.push_back(0.8);
    input_surface.insertKnot_u(new_knots);
    input_surface.insertKnot_v(new_knots);
    input_surface.swapParameterDirection();
    outfile2 << header_sf << input_surface;

    // We should now be ready to perform some tests; for istance check that
    // the iso-curves really do lie on our presumed Gordon surface.

    // mesh_curves is now ordered according to params.
    vector<shared_ptr<SplineCurve> > iso_curves;
    for (int i = 0; i < nmb_u_crvs; ++i)
	iso_curves.push_back(shared_ptr<SplineCurve>
			     (pGordonSurface->constParamCurve(params[i],
							      true)));
    for (int i = nmb_u_crvs; i < (int)mesh_curves.size(); ++i)
	iso_curves.push_back(shared_ptr<SplineCurve>
			     (pGordonSurface->constParamCurve(params[i],
							      false)));

    // As for now testing first u-curve.
//     double tparx = 0;
//     double tpary = 0;
//     Point surf_point, curve_point;
//     while (tpary != 1001) {
//  	cout << "Specify tpar to test difference between curves: ";
//  	cin >> tpary;
//  	mesh_curves[0].point(the_tpoint, tpary);
//  	raised_curve.point(raised_tpoint, tpar);
//  	cout << " Distance evaluated in tpar: " << 
//  	    raised_tpoint.dist(the_tpoint) << endl;
//     }

    return 0;

}
