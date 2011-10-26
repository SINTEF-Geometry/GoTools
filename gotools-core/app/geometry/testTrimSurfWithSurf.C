//===========================================================================
//                                                                           
// File: testTrimSurfWithSurf.C                                              
//                                                                           
// Created: Thu Nov 21 16:37:53 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: testTrimSurfWithSurf.C,v 1.10 2005-09-23 12:48:01 oan Exp $
//                                                                           
// Description: Given input of two spline surfaces, return the surfaces
//              trimmed w/eachother.
//                                                                           
//===========================================================================


#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedUtils.h"

#include <fstream>
#include <vector>


using namespace Go;
using std::vector;
using std::shared_ptr;


int main(int argc, char** argv)
{
    if (argc != 5) {
	MESSAGE("Usage: surf1 surf2 epsge outfile");
	return 1;
    }

    std::ifstream infile1(argv[1]);
    std::ifstream infile2(argv[2]);
    double epsge = atof(argv[3]);
    std::ofstream outfile(argv[4]);

    ObjectHeader header;
    shared_ptr<SplineSurface> surf1(new SplineSurface());
    shared_ptr<SplineSurface> surf2(new SplineSurface());
    header.read(infile1);
    surf1->read(infile1);
    header.read(infile2);
    surf2->read(infile2);

    vector<shared_ptr<BoundedSurface> > trimmed_sfs =
	BoundedUtils::trimSurfWithSurf(surf1, surf2, epsge);

    for (size_t ki = 0; ki < trimmed_sfs.size(); ++ki) {
	trimmed_sfs[ki]->writeStandardHeader(outfile);
	trimmed_sfs[ki]->write(outfile);
    }

    return 0;
}
