//===========================================================================
//                                                                           
// File: test_ellipse.C                                                      
//                                                                           
// Created: Mon Nov  2 13:56:38 2009                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include <fstream>
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include <memory>

using namespace Go;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::shared_ptr;

int main(int argc, char* argv[] )
{
    if (argc != 5)
    {
	std::cout << "Usage: " << argv[0]
		  << " input_ellipse tmin tmax output_subcurve" << endl;
	return -1;
    }

    // Open input surface file
    ifstream is(argv[1]);
    double tmin(atof(argv[2]));
    double tmax(atof(argv[3]));
    ofstream os(argv[4]);
    if (is.bad())
    {
	std::cout << "Bad or no input filename" << std::endl;
	return -1;
    }

    // Read surface from file
    ObjectHeader head;
    Ellipse ellipse; // Typically: centre, dir, normal, r1, r2.
    is >> head;
    ASSERT(head.classType() == Ellipse::classType());
    is >> ellipse;

    if (tmin < ellipse.startparam() || tmax > ellipse.endparam())
    {
	std::cout << "tmin or tmax outside domain of ellipse." << std::endl;
	return -1;
    }

    std::cout << "Writing to file." << std::endl;

    // Extract subcurve, write to file.
    ellipse.setParamBounds(tmin, tmax);

    shared_ptr<SplineCurve> sub_ellipse(ellipse.geometryCurve());
 
    sub_ellipse->writeStandardHeader(os);
    sub_ellipse->write(os);

    return 0;
}
