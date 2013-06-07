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
    shared_ptr<ParamSurface> surf1(new SplineSurface());
    surf1->read(input1);
    input1.close();
    
    double aepsge;
    aepsge = atof(argv[2]);


    shared_ptr<ParamSurfaceInt> ssurfint =
	shared_ptr<ParamSurfaceInt>(new SplineSurfaceInt (surf1));

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
	shared_ptr<ParamObjectInt> bd_obj2 = bd_obj[kh]->getObject();
	shared_ptr<ParamCurveInt> bd_cv = 
	    dynamic_pointer_cast<ParamCurveInt, ParamObjectInt>(bd_obj2);  


	SfCvIntersector sfcvintersect (ssurfint, bd_cv, aepsge, &sfsfintersect,
				       bd_obj[kh]->getDir()+2, 
				       bd_obj[kh]->getPar());
	sfcvintersect.setSelfintCase(1);
	sfcvintersect.compute();

	sfcvintersect.postIterateBd();

	//shared_ptr<IntersectionPool> pool = sfcvintersect.getIntPool();
	//pool->removeBoundaryIntersections(); 
	sfcvintersect.getIntPool()->cleanUpPool(0, aepsge);
	sfcvintersect.getIntPool()->removeBoundaryIntersections(false);

	vector<shared_ptr<IntersectionPoint> > intpts;
	vector<shared_ptr<IntersectionCurve> > intcrv;

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
		shared_ptr<IntersectionPoint> currpt
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

	vector<shared_ptr<IntersectionPoint> > intpts2;
	vector<shared_ptr<IntersectionCurve> > intcrv2;
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
		shared_ptr<IntersectionPoint> currpt
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
