//==========================================================================
//                                                                          
// File: testProjectCurvePoint.C                                             
//                                                                          
// Created: Thu Jan 29 13:55:21 2009                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: testProjectCurvePoint.C,v 1.1 2009-01-30 15:09:05 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/Point.h"
#include <fstream>


using namespace Go;
using std::shared_ptr;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
    if (argc != 6) {
	MESSAGE("Usage: surface, closed_dir_u, closed_dir_v, "
		"space_curve, cv_par");
	return 0;
    }

    // Read input arguments

    ifstream sffile(argv[1]);
    ALWAYS_ERROR_IF(sffile.bad(), "Input file not found or file corrupt");
    ObjectHeader header;
    header.read(sffile);
    SplineSurface sf;
    sf.read(sffile);

    bool closed_dir_u = (atoi(argv[2]) == 0 ? false : true);
    bool closed_dir_v = (atoi(argv[3]) == 0 ? false : true);

    ifstream cvfile(argv[4]);
    ALWAYS_ERROR_IF(cvfile.bad(), "Input file not found or file corrupt");
    header.read(cvfile);
    SplineCurve space_cv;
    space_cv.read(cvfile);

    double cv_par = atof(argv[5]);

    // Run the test
    shared_ptr<Point> pt 
	= CreatorsUtils::projectCurvePoint(sf, closed_dir_u, closed_dir_v,
					   &space_cv, cv_par);

    cout << "pt:  " << *pt << endl;

    return 0;
}

