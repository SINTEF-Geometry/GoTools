//==========================================================================
//                                                                          
// File: example_ImplicitizePointCloudAlgo.C
//                                                                          
// Created: Fri Jul 23 14:30:26 2004                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision:
// $Id: example_ImplicitizePointCloudAlgo.C,v 1.10 2006-03-31 09:09:04 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/ImplicitizePointCloudAlgo.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include <fstream>


using namespace Go;
using namespace std;


int main()
{
    // Read the point cloud from file
    ifstream input("data/point_cloud.g2");
    ObjectHeader header;
    PointCloud3D cloud;
    input >> header >> cloud;

    // Degree of surface
    int degree = 4;

    // Implicitize
    ImplicitizePointCloudAlgo implicitize(cloud, degree);
    implicitize.perform();

    // Get result
    BernsteinTetrahedralPoly implicit;
    BaryCoordSystem3D bc;
    double sigma_min;
    implicitize.getResultData(implicit, bc, sigma_min);

    // Write out implicit function
    ofstream output("data/implicit_point_cloud.dat");
    output << implicit << endl
	   << bc << endl;
    cout << "Data written to data/implicit_point_cloud.dat" << endl;

    return 0;
}
