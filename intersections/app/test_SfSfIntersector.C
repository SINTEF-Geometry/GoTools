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


using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::cerr;
using std::vector;
using std::setprecision;
using std::string;
using namespace Go;


int main(int argc, char** argv)
{

#ifndef _MSC_VER
//     _putenv("DEBUG=1");
//     _putenv("DEBUG_PAR=1");
//     _putenv("DEBUG_FINISH=1");
//     _putenv("DEBUG_IMPL=1");
    setenv("DEBUG", "1", 1);
    setenv("DEBUG_PAR", "1", 1);
    setenv("DEBUG_FINISH", "1", 1);
    //setenv("DEBUG_IMPL", "1", 1);
    setenv("SUBDIV_SFSF", "1", 1);
    setenv("DEBUG_MOVE", "1", 1);
    //setenv("DO_REPAIR", "1", 1);

    if (getenv("DEBUG") && (*getenv("DEBUG"))=='1')
	cout << "DEBUG=1" << endl;
    if (getenv("DEBUG_PAR") && (*getenv("DEBUG_PAR"))=='1')
	cout << "DEBUG_PAR=1" << endl;
    if (getenv("DEBUG_FINISH") && (*getenv("DEBUG_FINISH"))=='1')
	cout << "DEBUG_FINISH=1" << endl;
    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1')
	cout << "DEBUG_IMPL=1" << endl;
    if (getenv("SUBDIV_SFSF") && (*getenv("SUBDIV_SFSF"))=='1')
	cout << "SUBDIV_SFSF=1" << endl;
    if (getenv("DEBUG_MOVE") && (*getenv("DEBUG_MOVE"))=='1')
	cout << "DEBUG_MOVE=1" << endl;
    if (getenv("DO_REPAIR") && (*getenv("DO_REPAIR"))=='1')
	cout << "DO_REPAIR=1" << endl;
#endif // _MSC_VER


    if (argc != 4 && argc != 5 && argc != 6) {
	cout << "Usage: test_SfSfIntersector FileSf1 FileSf2 "
	     <<" aepsge (OutputFile)(selfintflag)"
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

    double aepsge;
    aepsge = atof(argv[3]);

    ofstream tmpout;
    if (argc >= 5)
	tmpout.open(argv[4]);
    ostream& out = tmpout;//(argc == 6) ? tmpout : cout;

    int selfint_flag = 0;
    if (argc == 6)
      selfint_flag = atoi(argv[5]);

    int has_sing = false;
//     if (argc == 5)
// 	has_sing = atoi(argv[4]);
    double sing[4];
//     if (has_sing) {
// 	printf("Give singularity: ");
// 	scanf("%lf",sing);
// 	scanf("%lf",sing+1);
// 	sing[2] = sing[0];
// 	sing[3] = sing[1];
//     }

//     cout << setprecision(15);

    // Setup the SplineSurfaceInt objects
    shared_ptr<ParamGeomInt> ssurfint1 =
	shared_ptr<ParamGeomInt>(new SplineSurfaceInt (surf1));
    shared_ptr<ParamGeomInt> ssurfint2 =
	shared_ptr<ParamGeomInt>(new SplineSurfaceInt (surf2));

    // Set up and run the intersection algorithm
    SfSfIntersector sfsfintersect (ssurfint1, ssurfint2, aepsge);
    if (has_sing)
	sfsfintersect.setHighPriSing(sing);
    if (selfint_flag > 0)
	sfsfintersect.setSelfintCase(selfint_flag); 
    clock_t start_compute = clock();
    sfsfintersect.compute();
    clock_t end_compute = clock();

    // Get the results
    vector<shared_ptr<IntersectionPoint> > intpts;
    vector<shared_ptr<IntersectionCurve> > intcrv;
    sfsfintersect.getResult(intpts, intcrv);
    cout << "Number of points: " << intpts.size() << endl;
    cout << "Number of curves: " << intcrv.size() << endl;

    // Write regression test data
    out << setprecision(8);
    int npoints = (int)intpts.size();
    int ncurves = (int)intcrv.size();
    vector<double> par;
    out << "-1 ______________Case sf-sf intersector______________"
	   << endl;
    out << "0 Surf-Surf:    npoints = " << npoints
	   << "    ncurves = " << ncurves << endl;
    for (int i = 0; i < npoints; ++i) {
	out << "1 Point: ";
	par = intpts[i]->getPar();
	for (int j=0; j < int(par.size()); j++) {
	    out << par[j] << " ";
	}
	out << intpts[i]->getDist() << endl;
    }
    for (int i = 0; i < ncurves; ++i) {
	out << "2 Curve:    Start = ";
	int nguide = intcrv[i]->numGuidePoints(); 
	shared_ptr<IntersectionPoint> startpt
	    = intcrv[i]->getGuidePoint(0);
	par = startpt->getPar();
	for (int j=0; j < int(par.size()); j++) {
	    out << par[j] << " ";
	}
	out << startpt->getDist() << "   End = ";
	shared_ptr<IntersectionPoint> endpt
	    = intcrv[i]->getGuidePoint(nguide-1);
	par = endpt->getPar();
	for (int j=0; j < int(par.size()); j++) {
	    out << par[j] << " ";
	}
	out << endpt->getDist() << endl;
    }

    // Write out isolated intersection points
    string str = "tmp_pnt.g2";
    std::ofstream outf(str.c_str());
    for (int ki=0; ki<int(intpts.size()); ki++) {
	std::vector<IntersectionPoint*> neighbours;
	intpts[ki]->getNeighbours(neighbours);
	std::vector<double> par = intpts[ki]->getPar();
	for (int kj=0; kj<int(par.size()); kj++)
	    std::cout << par[kj] << " ";
	std::cout << neighbours.size() << std::endl;

	Point pt1 = surf1->point(par[0], par[1]);
	Point pt2 = surf2->point(par[2], par[3]);

	outf << "400 1 0 4 255 0 0 255" << std::endl;
	outf << "1" << std::endl;
	pt1.write(outf);
	outf << "\n";

	outf << "400 1 0 4  0 255 0 255" << std::endl;
	outf << "1" << std::endl;
	pt2.write(outf);
	outf << "\n";
      
    }
    cout << ">> Intersection points written to " << str << endl;

    // Write out intersection curves as line-trains
    str = "tmp_crv.g2";
    std::ofstream outg(str.c_str());
    for (int ki=0; ki<int(intcrv.size()); ki++) {
	std::vector<double> guide_pt;
	int nguide = intcrv[ki]->numGuidePoints(); 
	guide_pt.reserve(3*nguide);
	for (int kj=0; kj<nguide; kj++) {
	    shared_ptr<IntersectionPoint> currpt
		= intcrv[ki]->getGuidePoint(kj);
	    Point pt = currpt->getPoint();
	    guide_pt.insert(guide_pt.end(), pt.begin(), pt.end());
	    if (kj>0 && kj<nguide-1)
		guide_pt.insert(guide_pt.end(), pt.begin(), pt.end()); 
	}
	LineCloud line(&guide_pt[0], nguide-1);
	line.writeStandardHeader(outg);
	line.write(outg);
    }
    cout << ">> Intersection curves (line-trains) written to "
	 << str << endl;

    
    double march_tol = 100.0*aepsge;
    //double march_tol = aepsge;
    double ang_tol = 0.01;

    // clock_t start_refine = clock();
    for (int ki=0; ki<int(intcrv.size()); ki++) {
	//      if (ki ==2) 
	cout << "Entering refine..." << endl;
	intcrv[ki]->refine(march_tol, ang_tol);
// 	intcrv[ki]->refine(march_tol, ang_tol);
// 	intcrv[ki]->refine(march_tol, ang_tol);
// 	intcrv[ki]->refine(march_tol, ang_tol);
// 	intcrv[ki]->refine(march_tol, ang_tol);
    }
    // clock_t end_refine = clock();

    // Write out refined intersection curves as line-trains
    str = "tmp2_crv.g2";
    std::ofstream outg2(str.c_str());
    for (int ki=0; ki<int(intcrv.size()); ki++) {
	std::vector<double> guide_pt;
	int nguide = intcrv[ki]->numGuidePoints(); 
	guide_pt.reserve(3*nguide);
	for (int kj=0; kj<nguide; kj++) {
	    shared_ptr<IntersectionPoint> currpt
		= intcrv[ki]->getGuidePoint(kj);
	    Point pt = currpt->getPoint();
	    guide_pt.insert(guide_pt.end(), pt.begin(), pt.end());
	    if (kj>0 && kj<nguide-1)
		guide_pt.insert(guide_pt.end(), pt.begin(), pt.end()); 
	}
	LineCloud line(&guide_pt[0], nguide-1);
	line.writeStandardHeader(outg2);
	line.write(outg2);
    }
    cout << ">> Refined intersection curves (line-trains) written to "
	 << str << endl;

    // Write out intersection curves as spline curve
    str = "tmp3_crv.g2";
    std::ofstream outg3(str.c_str());
    for (int ki=0; ki<int(intcrv.size()); ki++) {
	shared_ptr<ParamCurve> splcrv = intcrv[ki]->getCurve();
	if (splcrv.get()) {
// 	    splcrv->writeStandardHeader(outg3);
	    outg3 << "100 1 0 4 50 205 50 255" << endl;
	    splcrv->write(outg3);
	}
    }
    cout << ">> Intersection curve (spline curve) written to "
	 << str << endl;

    // Write out intersection curves in the parameter plane of the
    // first surface
    str = "tmp4_crv.g2";
    std::ofstream outg4(str.c_str());
    for (int ki=0; ki<int(intcrv.size()); ki++) {
	shared_ptr<ParamCurve> splcrv = intcrv[ki]->getParamCurve(1);
	if (splcrv.get()) {
	    splcrv->writeStandardHeader(outg4);
	    splcrv->write(outg4);
	}
    }
    cout << ">> Intersection curves (in parameter plane of 1st surface) "
	 << "written to " << str << endl;

    // Write out intersection curves in the parameter plane of the
    // second surface
    str = "tmp5_crv.g2";
    std::ofstream outg5(str.c_str());
    for (int ki=0; ki<int(intcrv.size()); ki++) {
	shared_ptr<ParamCurve> splcrv = intcrv[ki]->getParamCurve(2);
	if (splcrv.get()) {
	    splcrv->writeStandardHeader(outg5);
	    splcrv->write(outg5);
	}
    }
    cout << ">> Intersection curves (in parameter plane of 2nd surface) "
	 << "written to " << str << endl;
    

    cout << "Time in intersector:\t" 
	 << (double)(end_compute - start_compute) / double(1000)
	 << " msec. " << endl;
//     cout << "Time spent for refining:\t"
// 	 << (end_refine - start_refine) / double(1000)
// 	 << " msec." << endl;

    return 0;
}
