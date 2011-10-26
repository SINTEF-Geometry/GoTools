//==========================================================================
//                                                                          
// File: example_ImplicitizeCurveAlgo.C
//                                                                          
// Created: Mon Mar 10 12:31:29 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision:
// $Id: example_ImplicitizeCurveAlgo.C,v 1.15 2006-03-31 09:09:04 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/ImplicitizeCurveAlgo.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/implicitization/BernsteinTriangularPoly.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include <fstream>


using namespace Go;
using namespace std;


int main()
{
    // Read the curve from file
    ifstream input("data/curve.g2");
    ObjectHeader header;
    SplineCurve curve;
    input >> header >> curve;

    // Degree of curve
    int degree = 6;

    // Implicitize
    ImplicitizeCurveAlgo implicitize(curve, degree);
    implicitize.perform();

    // Get result
    BernsteinTriangularPoly implicit;
    BaryCoordSystem2D bc;
    double sigma_min;
    implicitize.getResultData(implicit, bc, sigma_min);

    // Write out implicit surface
    ofstream output("data/implicit_curve.dat");
    output << implicit << endl
	   << bc << endl;
    cout << "Data written to data/implicit_curve.dat" << endl;

    return 0;
}
