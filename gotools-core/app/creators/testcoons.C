//===========================================================================
//                                                                           
// File: testcoons.C                                                         
//                                                                           
// Created: Thu Apr  5 11:43:29 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: testcoons.C,v 1.4 2006-04-19 10:48:08 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/CurveLoop.h"
#include <memory>
#include <fstream>

using namespace Go;
using std::shared_ptr;;

int main(int argc, char** argv)
{
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    if (input.bad() || output.bad()) {
	std::cerr << "File error (no file or corrupt file specified)."
		  << std::endl;
	return 1;
    }
    ObjectHeader header;
    std::vector< std::shared_ptr<ParamCurve> > curves(4);
    for (int i = 0; i < 4; ++i) {
	curves[i] = std::shared_ptr<ParamCurve>(new SplineCurve);
	input >> header;
	input >> (*curves[i]);
    }
    shared_ptr<SplineSurface> sf
	(CoonsPatchGen::createCoonsPatch(CurveLoop(curves, 0.0)));
    sf->writeStandardHeader(output);
    sf->write(output);
}
