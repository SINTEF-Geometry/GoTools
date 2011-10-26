//==========================================================================
//                                                                          
// File: ImplicitizeSurfaceAlgo.C                                            
//                                                                          
// Created: Wed Feb  5 16:19:24 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ImplicitizeSurfaceAlgo.C,v 1.17 2006-03-31 09:09:08 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/ImplicitizeSurfaceAlgo.h"
#include "GoTools/implicitization/ImplicitUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "newmatio.h"
#include "newmat.h"


using namespace std;


namespace Go {


//==========================================================================
void ImplicitizeSurfaceAlgo::perform()
//==========================================================================
{
    // Create barycentric coordinate system
    create_bary_coord_system3D(surf_, bc_);

    // Convert spline curve to barycentric coordinates
    SplineSurface surf_bc;
    cart_to_bary(surf_, bc_, surf_bc);

    // Check if the surface has a single patch
    bool single_patch = (surf_.order_u() == surf_.numCoefs_u()
			 && surf_.order_v() == surf_.numCoefs_v());

    // Make the matrix of numerical coefficients (the D-matrix). Any
    // vector in the nullspace of this matrix will be a solution.
    vector<vector<double> > mat;
    if (single_patch) {
	make_matrix(surf_bc, deg_, mat);
    } else {
	// The matrices from all the patches are stacked on top of each
	// other
	vector<SplineSurface> patches;
	splitSurfaceIntoPatches(surf_bc, patches);
	int num = (int)patches.size();
	make_matrix(patches[0], deg_, mat);
	vector<vector<double> > tmp;
	for (int i = 1; i < num; ++i) {
	    make_matrix(patches[i], deg_, tmp);
	    mat.insert(mat.end(), tmp.begin(), tmp.end());
	}
    }

    // Find the nullspace and construct the implicit function.
    vector<double> b;
    make_implicit_svd(mat, b, sigma_min_);

    // Set the coefficients
    implicit_ = BernsteinTetrahedralPoly(deg_, b);

    return;
}


//==========================================================================


} // namespace Go
