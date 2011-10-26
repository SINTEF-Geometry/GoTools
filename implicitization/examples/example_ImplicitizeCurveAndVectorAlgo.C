//==========================================================================
//                                                                          
// File: example_ImplicitizeCurveAndVectorAlgo.C
//                                                                          
// Created: Wed Jun  8 10:30:11 2005                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision:
// $Id: example_ImplicitizeCurveAndVectorAlgo.C,v 1.5 2006-03-31 09:09:04 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/ImplicitizeCurveAndVectorAlgo.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include <fstream>


using namespace Go;
using namespace std;


int main()
{
    // Read the curve from file
    std::ifstream input1("data/curve3d.g2");
    ObjectHeader header;
    SplineCurve curve;
    input1 >> header >> curve;

    // Read the vector from file
    std::ifstream input2("data/vector.dat");
    Point vec(3);
    input2 >> vec;

    // Degree of surface
    int degree = 6;

    // Implicitize
    ImplicitizeCurveAndVectorAlgo implicitize(curve, vec, degree);
    int info = 0;
    info = implicitize.perform();

    // Get result
    BernsteinTetrahedralPoly implicit;
    BaryCoordSystem3D bc;
    double sigma_min;
    implicitize.getResultData(implicit, bc, sigma_min);

    // Write out implicit function
    ofstream output("data/implicit_curve_and_vector.dat");
    output << implicit << endl
	   << bc << endl;
    cout << "Data written to data/implicit_curve_and_vector.dat" << endl;

    return 0;
}
