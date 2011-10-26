//==========================================================================
//                                                                          
// File: selfIntersectBoundary.C
//
// Created:
//                                                                          
// Author: Vibeke Skytt
//                                                                          
// Revision: $Id: selfIntersectBoundary.C,v 1.9 2007-11-01 14:31:19 vsk Exp $
//                                                                          
// Description: 
//                                                                          
//==========================================================================


#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/intersections/SfCvIntersector.h"
#include "GoTools/intersections/SfSelfIntersector.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <memory>
#include <fstream>
#include <iomanip>
#include <time.h>


using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using namespace Go;
using std::shared_ptr;
using std::dynamic_pointer_cast;


int main(int argc, char** argv)
{
 
    if (argc != 3) {
	cout << "Usage: selfIntersectBoundary.C FileSf1 aepsge"
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


    std::shared_ptr<ParamSurfaceInt> ssurfint =
	std::shared_ptr<ParamSurfaceInt>(new SplineSurfaceInt (surf1));

    std::ofstream outf("tmp_pnt.g2");
    std::ofstream outg("tmp_crv.g2");

    SfSfIntersector sfsfintersect(ssurfint, ssurfint, aepsge);

    // Fetch boundary curves
    std::vector<shared_ptr<BoundaryGeomInt> > bd_obj;
    ssurfint->getBoundaryObjects(bd_obj);
    size_t ki, kr, kh;
    int kj;
    int nmb_orig = 0;
    for (kh=0; kh<bd_obj.size(); kh++) {
	std::shared_ptr<ParamObjectInt> bd_obj2 = bd_obj[kh]->getObject();
	std::shared_ptr<ParamCurveInt> bd_cv = 
	    dynamic_pointer_cast<ParamCurveInt, ParamObjectInt>(bd_obj2);  


	SfCvIntersector sfcvintersect (ssurfint, bd_cv, aepsge, &sfsfintersect,
				       bd_obj[kh]->getDir()+2, 
				       bd_obj[kh]->getPar());
	sfcvintersect.setSelfintCase(1);
	sfcvintersect.compute();

	sfcvintersect.postIterateBd();

	//std::shared_ptr<IntersectionPool> pool = sfcvintersect.getIntPool();
	//pool->removeBoundaryIntersections(); 
	sfcvintersect.getIntPool()->cleanUpPool(0, aepsge);
	sfcvintersect.getIntPool()->removeBoundaryIntersections(false);

	vector<std::shared_ptr<IntersectionPoint> > intpts;
	vector<std::shared_ptr<IntersectionCurve> > intcrv;

	sfcvintersect.getResult(intpts, intcrv);
	printf("Number of points: %d \n", int(intpts.size()));
	printf("Number of curves: %d \n", int(intcrv.size()));

	std::cout << "Intersection points: " << std::endl;
	for (ki=0; ki<intpts.size(); ki++) {
	    std::vector<IntersectionPoint*> neighbours;
	    intpts[ki]->getNeighbours(neighbours);
	    std::vector<double> par = intpts[ki]->getPar();
	    for (kr=0; kr<par.size(); kr++)
		std::cout << par[kr] << " ";
	    std::cout << neighbours.size() << "  ";
	    std::cout << intpts[ki]->getDist() << std::endl;

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

	
	std::cout << "Guide points: " << std::endl;
	for (ki=0; ki<intcrv.size(); ki++) {
	    std::vector<double> guide_pt;
	    int nguide = intcrv[ki]->numGuidePoints(); 
	    guide_pt.reserve(3*nguide);
	    for (kj=0; kj<nguide; kj++) {
		std::shared_ptr<IntersectionPoint> currpt
		    = intcrv[ki]->getGuidePoint(kj);
		std::vector<IntersectionPoint*> neighbours;
		currpt->getNeighbours(neighbours);
		std::vector<double> par = currpt->getPar();
		for (kr=0; kr<par.size(); kr++)
		    std::cout << par[kr] << " ";
		std::cout << neighbours.size() << "  ";
		std::cout << currpt->getDist() << std::endl;

		Point pt = currpt->getPoint();
		guide_pt.insert(guide_pt.end(), pt.begin(), pt.end());
		if (kj>0 && kj<nguide-1)
		    guide_pt.insert(guide_pt.end(), pt.begin(), pt.end()); 
	    }
	    LineCloud line(&guide_pt[0], nguide-1);
	    line.writeStandardHeader(outg);
	    line.write(outg);
	}
 	sfsfintersect.getIntPool()->includeReducedInts
	    (sfcvintersect.getIntPool());

	sfsfintersect.postIterate3(nmb_orig, (int)kh);
	nmb_orig = sfsfintersect.getIntPool()->numIntersectionPoints();
   }

    sfsfintersect.getIntPool()->makeIntersectionCurves();

	vector<std::shared_ptr<IntersectionPoint> > intpts2;
	vector<std::shared_ptr<IntersectionCurve> > intcrv2;
	sfsfintersect.getResult(intpts2, intcrv2);
	printf("sfsfintersect. Number of points: %d \n", int(intpts2.size()));
	printf("sfsfintersect. Number of curves: %d \n", int(intcrv2.size()));

	double dist, dist1, dist2;
	double param[2];
	Point pos, int_pos, int_pos1, int_pos2;
	SingularityType type;
	bool near_sing;
	std::cout << "Intersection points: " << std::endl;
	for (ki=0; ki<intpts2.size(); ki++) {
	    std::vector<IntersectionPoint*> neighbours;
	    intpts2[ki]->getNeighbours(neighbours);
	    std::vector<double> par = intpts2[ki]->getPar();
	    param[0] = 0.5*(par[0] + par[2]);
	    param[1] = 0.5*(par[1] + par[3]);
	    surf1->point(pos, param[0], param[1]);
	    int_pos = intpts2[ki]->getPoint();
	    int_pos1 = intpts2[ki]->getPoint1();
	    int_pos2 = intpts2[ki]->getPoint2();
	    dist = int_pos.dist(pos);
	    dist1 = int_pos1.dist(pos);
	    dist2 = int_pos2.dist(pos);
	    type = intpts2[ki]->getSingularityType();
	    near_sing = intpts2[ki]->isNearSingular();
	    for (kr=0; kr<par.size(); kr++)
		std::cout << par[kr] << " ";
	    std::cout << neighbours.size() << "  ";
	    std::cout << intpts2[ki]->getDist();
	    std::cout << "' Type: " << type << " " << near_sing;
	    std::cout << ", Dist to symm pt: " << dist;
	    std::cout << " " << dist1 << " " << dist2 << std::endl;

	}

	std::cout << "Guide points: " << std::endl;
	for (ki=0; ki<intcrv2.size(); ki++) {
	    std::vector<double> guide_pt;
	    int nguide = intcrv2[ki]->numGuidePoints(); 
	    guide_pt.reserve(3*nguide);
	    for (kj=0; kj<nguide; kj++) {
		std::shared_ptr<IntersectionPoint> currpt
		    = intcrv2[ki]->getGuidePoint(kj);
		std::vector<IntersectionPoint*> neighbours;
		currpt->getNeighbours(neighbours);
		std::vector<double> par = currpt->getPar();
		for (kr=0; kr<par.size(); kr++)
		    std::cout << par[kr] << " ";
		std::cout << neighbours.size() << "  ";
		std::cout << currpt->getDist() << std::endl;

	    }
	}
    return 0;
}
