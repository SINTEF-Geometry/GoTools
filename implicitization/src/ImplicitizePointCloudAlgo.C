//==========================================================================
//                                                                          
// File: ImplicitizePointCloudAlgo.C                                         
//                                                                          
// Created: Mon Jul 19 14:56:42 2004                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision:
// $Id: ImplicitizePointCloudAlgo.C,v 1.7 2006-03-31 09:09:07 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/ImplicitizePointCloudAlgo.h"
#include "GoTools/implicitization/ImplicitUtils.h"
#include "GoTools/geometry/GeometryTools.h"


using namespace std;


namespace Go {


//==========================================================================
void ImplicitizePointCloudAlgo::perform()
//==========================================================================
{
    // Create barycentric coordinate system
    create_bary_coord_system3D(cloud_, bc_);

    // Convert point cloud to barycentric coordinates
    PointCloud4D cloud_bc;
    cart_to_bary(cloud_, bc_, cloud_bc);

    // Make the matrix of numerical coefficients (the D-matrix). Any
    // vector in the nullspace of this matrix will be a solution.
    vector<vector<double> > mat;
    make_matrix(cloud_bc, deg_, mat);

    // Find the nullspace and construct the implicit function.
    vector<double> b;
    make_implicit_svd(mat, b, sigma_min_);

    // Set the coefficients
    implicit_ = BernsteinTetrahedralPoly(deg_, b);

    return;
}


//==========================================================================


} // namespace Go
