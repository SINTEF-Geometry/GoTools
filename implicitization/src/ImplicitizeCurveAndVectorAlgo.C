//==========================================================================
//                                                                          
// File: ImplicitizeCurveAndVectorAlgo.C                                     
//                                                                          
// Created: Tue Jun  7 13:37:38 2005                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ImplicitizeCurveAndVectorAlgo.C,v 1.4 2006-05-03 11:50:19 sbr Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/implicitization/ImplicitizeCurveAndVectorAlgo.h"
#include "GoTools/implicitization/ImplicitizeSurfaceAlgo.h"
// #include "GoTools/implicitization/ImplicitUtils.h"
// #include "GoTools/geometry/GeometryTools.h"


using namespace std;


namespace Go {


//==========================================================================
int ImplicitizeCurveAndVectorAlgo::perform()
//==========================================================================
{
    // Make a ruled surface

    // Number of coefficients
    int ncoefsu = crv_.numCoefs();
    int ncoefsv = 2;

    // Orders
    int ordu = crv_.order();
    int ordv = 2;

    // Knot vectors
    typedef vector<double>::iterator iter;
    iter knotsu_begin = crv_.basis().begin();
    double knotsv_begin[] = { 0.0, 0.0, 1.0, 1.0 };

    // Rationality and dimension
    bool rational = crv_.rational();
    int dim = 3;
    int effdim = (rational ? dim+1 : dim);

    // Now the tricky part - the coefficients

    // First we need a scale. (We normalize the vector pt_.)
    BoundingBox box = crv_.boundingBox();
    Point diagonal = box.high() - box.low();
    double scale = diagonal.length();
    pt_.normalize();

    vector<double> coefs(effdim * ncoefsu * ncoefsv);
    iter curr = (rational) ? crv_.rcoefs_begin() : crv_.coefs_begin();
    iter down = coefs.begin(); // First row of coefs ("down")
    iter up = coefs.begin() + effdim * ncoefsu; // Second row of coefs
						// ("up")
    for (int i = 0; i < ncoefsu; ++i) {
	double w = (rational ? *(curr + dim) : 1.0);
	Point tmppt(curr, curr+dim);
	// Set a down-point
	tmppt -= w * 0.5 * scale * pt_;
	copy(tmppt.begin(), tmppt.end(), down);
	if (rational) {
  	    *(down + dim) = w;
	}
	down += effdim;
	// Set an up-point
	tmppt += w * scale * pt_;
	copy(tmppt.begin(), tmppt.end(), up);
	if (rational) {
	    *(up + dim) = w;
	}
	up += effdim;

	curr += effdim;
    }

    SplineSurface surf(ncoefsu, ncoefsv, ordu, ordv,
		       knotsu_begin, knotsv_begin, coefs.begin(), dim,
		       rational);

    // Run surface algorithm
    ImplicitizeSurfaceAlgo algo(surf, deg_);
    algo.setTolerance(tol_);
    algo.perform();
    algo.getResultData(implicit_, bc_, sigma_min_);

    return 0;
}


//==========================================================================


} // namespace Go
