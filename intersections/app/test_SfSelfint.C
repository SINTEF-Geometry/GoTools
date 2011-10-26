//==========================================================================
//                                                                          
// File: test_SfSelfint.C
//
// Created:
//                                                                          
// Author: Vibeke Skytt
//                                                                          
// Revision: $Id: test_SfSelfint.C,v 1.14 2007-11-01 14:31:19 vsk Exp $
//                                                                          
// Description: 
//                                                                          
//==========================================================================

#include <fstream>
#include <iomanip>
#include <memory>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/intersections/SfSelfIntersector.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <time.h>

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::setprecision;
using std::vector;
using namespace Go;
using std::shared_ptr;


int main(int argc, char** argv)
{
#ifndef _MSC_VER
    //setenv("DEBUG", "1", 1);
    //setenv("DEBUG_PAR", "1", 1);
    setenv("DEBUG_FINISH", "1", 1);
    //setenv("DEBUG_IMPL", "1", 1);
    setenv("DEBUG_SELFINT", "1", 1);
    setenv("DEBUG_DIV", "1", 1);
    setenv("DEBUG_DIV2", "1", 1);
    setenv("DEBUG_MOVE", "1", 1);

    if (getenv("DEBUG") && (*getenv("DEBUG"))=='1')
	cout << "DEBUG=1" << endl;
    if (getenv("DEBUG_PAR") && (*getenv("DEBUG_PAR"))=='1')
	cout << "DEBUG_PAR=1" << endl;
    if (getenv("DEBUG_FINISH") && (*getenv("DEBUG_FINISH"))=='1')
	cout << "DEBUG_FINISH=1" << endl;
    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1')
	cout << "DEBUG_IMPL=1" << endl;
    if (getenv("DEBUG_SELFINT") && (*getenv("DEBUG_SELFINT"))=='1')
	cout << "DEBUG_SELFINT=1" << endl;
#endif // _MSC_VER

 
    if (argc != 4) {
	cout << "Usage: test_SfSelfint FileSf1 aepsge reclevel"
	     << endl;
	return 0;
    }


    ObjectHeader header;

    // Read the surface from file
    ifstream input1(argv[1]);
    if (input1.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input1);
    std::shared_ptr<ParamSurface> surf1(new SplineSurface());
    surf1->read(input1);
    input1.close();
    
    double aepsge;
    aepsge = atof(argv[2]);
    int rec = atoi(argv[3]);


    std::shared_ptr<ParamSurfaceInt> ssurfint1 =
	std::shared_ptr<ParamSurfaceInt>(new SplineSurfaceInt (surf1));

    vector<std::shared_ptr<IntersectionPoint> > intpts;
    vector<std::shared_ptr<IntersectionCurve> > intcrv;
    SfSelfIntersector sfselfint (ssurfint1, aepsge);
    sfselfint.setMaxRec(rec);
    clock_t start_compute = clock();
    sfselfint.compute();
    clock_t end_compute = clock();

    sfselfint.getResult(intpts, intcrv);
    printf("Number of points: %d \n", int(intpts.size()));
    printf("Number of curves: %d \n", int(intcrv.size()));

    size_t ncurves = intcrv.size();
    vector<double> par;
    cout << setprecision(8);
    for (int i = 0; i < (int)ncurves; ++i) {
	std::cout << "2 Curve:    Start = ";
	int nguide = intcrv[i]->numGuidePoints(); 
	shared_ptr<IntersectionPoint> startpt
	    = intcrv[i]->getGuidePoint(0);
	par = startpt->getPar();
	for (int j=0; j < int(par.size()); j++) {
	    std::cout << par[j] << " ";
	}
	std::cout << startpt->getDist() << "   End = ";
	shared_ptr<IntersectionPoint> endpt
	    = intcrv[i]->getGuidePoint(nguide-1);
	par = endpt->getPar();
	for (int j=0; j < int(par.size()); j++) {
	    std::cout << par[j] << " ";
	}
	std::cout << endpt->getDist() << endl;
    }

    std::ofstream outf("tmp_pnt.g2");
    int ki, kj;
    for (ki=0; ki<int(intpts.size()); ki++) {
	std::vector<IntersectionPoint*> neighbours;
	intpts[ki]->getNeighbours(neighbours);
	std::vector<double> par = intpts[ki]->getPar();
	for (kj=0; kj<int(par.size()); kj++)
	    std::cout << par[kj] << " ";
	std::cout << intpts[ki]->getDist() << neighbours.size() << std::endl;

	Point pt1 = surf1->point(par[0], par[1]);
	Point pt2 = surf1->point(par[2], par[3]);

	outf << "400 1 0 4 255 0 0 255" << std::endl;
	outf << "1" << std::endl;
	pt1.write(outf);
	outf << "\n";

	outf << "400 1 0 4  0 255 0 255" << std::endl;
	outf << "1" << std::endl;
	pt2.write(outf);
	outf << "\n";
      
    }

    std::ofstream outg("tmp_crv.g2");
    for (ki=0; ki<int(intcrv.size()); ki++) {
	std::vector<double> guide_pt;
	int nguide = intcrv[ki]->numGuidePoints(); 
	guide_pt.reserve(3*nguide);
	for (kj=0; kj<nguide; kj++) {
	    std::shared_ptr<IntersectionPoint> currpt
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

    double march_tol = 1000.0*aepsge;
    double ang_tol = 0.01;

    std::ofstream outg3("tmp3_crv.g2");
    clock_t start_refine = clock();
    for (ki=0; ki<int(intcrv.size()); ki++) {
	//      if (ki ==2) 
	cout << "Entering refine..." << endl;
	ofstream krull("dumpcurve");
	intcrv[ki]->writeIPointsToStream(krull);
	krull.close();
	intcrv[ki]->refine(march_tol, ang_tol);
// 	intcrv[ki]->refine(march_tol, ang_tol);
// 	intcrv[ki]->refine(march_tol, ang_tol);
// 	intcrv[ki]->refine(march_tol, ang_tol);
// 	intcrv[ki]->refine(march_tol, ang_tol);
	std::shared_ptr<ParamCurve> splcrv = intcrv[ki]->getCurve();
	if (splcrv.get()) {
	    splcrv->writeStandardHeader(outg3);
	    splcrv->write(outg3);
	}
    }
    clock_t end_refine = clock();

    std::ofstream outg2("tmp2_crv.g2");
    for (ki=0; ki<int(intcrv.size()); ki++) {
	std::vector<double> guide_pt;
	int nguide = intcrv[ki]->numGuidePoints(); 
	guide_pt.reserve(3*nguide);
	for (kj=0; kj<nguide; kj++) {
	    std::shared_ptr<IntersectionPoint> currpt
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


    std::ofstream outg4("tmp4_crv.g2");
    for (ki=0; ki<int(intcrv.size()); ki++) {
	std::shared_ptr<ParamCurve> splcrv = intcrv[ki]->getParamCurve(1);
	if (splcrv.get()) {
	    splcrv->writeStandardHeader(outg4);
	    splcrv->write(outg4);
	}
    }

    std::ofstream outg5("tmp5_crv.g2");
    for (ki=0; ki<int(intcrv.size()); ki++) {
	std::shared_ptr<ParamCurve> splcrv = intcrv[ki]->getParamCurve(2);
	if (splcrv.get()) {
	    splcrv->writeStandardHeader(outg5);
	    splcrv->write(outg5);
	}
    }

    cout << "Time in intersector: "
	 << (double)(end_compute - start_compute) / double(1000) << " msec. " << endl;
    cout << "Time spent for refining: "
	 << (double)(end_refine - start_refine) / double(1000) << " msec." << endl;

    return 0;
}
