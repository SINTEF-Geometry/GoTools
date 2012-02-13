//==========================================================================
//                                                                          
// File: ImplicitizeCurveAlgo.C                                              
//                                                                          
// Created: Wed Feb  5 15:18:22 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ImplicitizeCurveAlgo.C,v 1.11 2006-03-31 09:09:07 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/ImplicitizeCurveAlgo.h"
#include "GoTools/implicitization/ImplicitUtils.h"
#include "GoTools/geometry/GeometryTools.h"


using namespace std;


namespace Go {


//==========================================================================
void ImplicitizeCurveAlgo::perform()
//==========================================================================
{
    // Create barycentric coordinate system
    create_bary_coord_system2D(curve_, bc_);

    // Convert spline curve to barycentric coordinates
    SplineCurve crv_bc;
    cart_to_bary(curve_, bc_, crv_bc);

    // Check if the curve has a single segment
    bool single_segment = (curve_.order() == curve_.numCoefs());

    // Make the matrix of numerical coefficients (the D-matrix). Any
    // vector in the nullspace of this matrix will be a solution.
    vector<vector<double> > mat;
    if (single_segment) {
	make_matrix(crv_bc, deg_, mat);
    } else {
	// The matrices from all the segments are stacked on top of each
	// other
	vector<SplineCurve> segments;
	GeometryTools::splitCurveIntoSegments(crv_bc, segments);
	int num = (int)segments.size();
	make_matrix(segments[0], deg_, mat);
	vector<vector<double> > tmp;
	for (int i = 1; i < num; ++i) {
	    make_matrix(segments[i], deg_, tmp);
	    mat.insert(mat.end(), tmp.begin(), tmp.end());
	}
    }

    // Find the nullspace and construct the implicit function.
    vector<double> b;
    make_implicit_gauss(mat, b);

    // We boldly assume there is no error for the Gaussian elimination
    sigma_min_ = 0.0;

    // Set the coefficients
    implicit_ = BernsteinTriangularPoly(deg_, b);

    return;
}


//==========================================================================


} // namespace Go
