//==========================================================================
//                                                                          
// File: example_ImplicitizeSurfaceAlgo.C
//                                                                          
// Created: Mon Mar 10 12:48:10 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision:
// $Id: example_ImplicitizeSurfaceAlgo.C,v 1.20 2006-03-31 09:09:05 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/ImplicitizeSurfaceAlgo.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include <fstream>


using namespace Go;
using namespace std;


int main()
{
    // Read the surface from file
    ifstream input("data/surface.g2");
    ObjectHeader header;
    SplineSurface surface;
    input >> header >> surface;

    // Degree of surface
    int degree = 4;

    // Implicitize
    ImplicitizeSurfaceAlgo implicitize(surface, degree);
    implicitize.perform();

    // Get result
    BernsteinTetrahedralPoly implicit;
    BaryCoordSystem3D bc;
    double sigma_min;
    implicitize.getResultData(implicit, bc, sigma_min);

    // Write out implicit function
    ofstream output("data/implicit_surface.dat");
    output << implicit << endl
	   << bc << endl;
    cout << "Data written to data/implicit_surface.dat" << endl;

    return 0;
}
