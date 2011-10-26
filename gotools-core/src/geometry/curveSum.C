#include "GoTools/geometry/GeometryTools.h"
#include <memory>

//***************************************************************************
//
// Implementation file of the free function curveSum defined in
// GeometryTools.h/
//
//***************************************************************************

using std::shared_ptr;
using std::dynamic_pointer_cast;
using std::vector;
using std::max;
using std::min;

namespace Go
{

  std::shared_ptr<SplineCurve>
  curveSum(const SplineCurve& crv1, double fac1,
	   const SplineCurve& crv2, double fac2, double num_tol)

    //********************************************************************
    // Addition of two signed SplineCurves, i.e. this function can
    // also be used for subtraction. The curves is assumed to live on
    // the same parameter domain, but may have different knot vectors.
    //********************************************************************
  {
    // Check input
    ALWAYS_ERROR_IF(fabs(crv1.startparam() - crv2.startparam()) > num_tol ||
		fabs(crv1.endparam() - crv2.endparam()) > num_tol,
		"Inconsistent parameter domain.");

    // For the time being
    if (crv1.rational() || crv2.rational()) {
	THROW("Sum of rational curves is not impelemented");
    }

    // Make copy of curves
    vector<shared_ptr<SplineCurve> > curves;
    curves.reserve(2);
    shared_ptr<SplineCurve> cv;
// #ifdef _MSC_VER
//     cv = shared_ptr<SplineCurve>(dynamic_cast<SplineCurve*>(crv1.clone()));
// #else
    cv = shared_ptr<SplineCurve>(crv1.clone());
// #endif
    curves.push_back(cv);
// #ifdef _MSC_VER
//     cv = shared_ptr<SplineCurve>(dynamic_cast<SplineCurve*>(crv2.clone()));
// #else
    cv = shared_ptr<SplineCurve>(crv2.clone());
// #endif
    curves.push_back(cv);

    // Make sure that the curves live on the same knot vector
//     double tol = 0.00001;
    try {
	unifyCurveSplineSpace(curves, num_tol);
    } catch (...) {
	THROW("Failed unifying spline spaces!");
    }

    // Add signed coefficients
    vector<double> coefs;
    int nmb_coefs = curves[0]->numCoefs();
    int dim = curves[0]->dimension();
    coefs.resize(dim*nmb_coefs);
    int ki;
    std::vector<double>::iterator c1 = curves[0]->coefs_begin();
    std::vector<double>::iterator c2 = curves[1]->coefs_begin();
    for (ki=0; ki<dim*nmb_coefs; ki++)
      coefs[ki] = fac1*c1[ki] + fac2*c2[ki];

    // Create output curve
    std::shared_ptr<SplineCurve>
	curvesum(new SplineCurve(nmb_coefs, 
				 curves[0]->order(),
				 curves[0]->basis().begin(),
				 &coefs[0], 
				 dim,
				 false));
    return curvesum;
  }
}

