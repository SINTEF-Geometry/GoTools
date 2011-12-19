//==========================================================================
//                                                                          
// File: pickSingularRegion.C
//
// Created:
//                                                                          
// Author: jbt
//                                                                          
// Revision: $Id: pickSingularRegion.C,v 1.1 2006-11-03 15:11:14 jbt Exp $
//                                                                          
// Description: Picks out and amplifies region around a singularity.
//                                                                          
//==========================================================================


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
