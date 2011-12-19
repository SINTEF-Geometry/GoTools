//===========================================================================
//                                                                           
// File: testApproxCurve.C                                                   
//                                                                           
// Created: Thu Jul 18 14:56:41 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: testApproxCurve.C,v 1.5 2005-12-22 14:52:34 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <fstream>


using namespace Go;
using std::vector;


int main(int argc, char* argv[])
{
    if (argc != 4) {
	MESSAGE("Usage: inputfile tolerance outputfile.");
	return 0;
    }

    // Read input arguments
    std::ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(), "Input file not found or file corrupt");

    double epsge = atof(argv[2]);
    std::ofstream outfile(argv[3]);

    // Input surface is to be a GoSplineCurve.
    ObjectHeader header;
    header.read(infile);
    shared_ptr<SplineCurve> crv(new SplineCurve());
    crv->read(infile);

    int max_iter = 20;
    vector<shared_ptr<SplineCurve> > crvs;
    crvs.push_back(crv);
    vector<Point> start_pt, end_pt;
    try {
	double max_dist;
	crv = shared_ptr<SplineCurve>
	    (CurveCreators::approxCurves(&crvs[0], &crvs[1],
					 start_pt, end_pt, epsge, max_dist, max_iter));
	if (max_dist > epsge) {
	    MESSAGE("Failed approximating within tolerance (" << epsge <<
		       "), using cv anyway. Dist: " << max_dist);
	}
    } catch (...) {
	MESSAGE("Failed approximating input curve, returning input curve.");
    }
    crv->writeStandardHeader(outfile);
    crv->write(outfile);
}

