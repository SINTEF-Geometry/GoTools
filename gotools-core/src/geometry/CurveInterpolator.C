#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/SplineCurve.h"

using std::vector;

namespace Go
{
  SplineCurve* 
  CurveInterpolator::regularInterpolation(const BsplineBasis& basis,
					    vector<double>& par,
					    vector<double>& points,
					    int dimension,
					    bool rational,
					    vector<double>& weights)
  {
    // Check input
    ASSERT(par.size() == points.size()/dimension);
    ASSERT(basis.numCoefs() == (int)par.size());

    vector<double> points2;
    if (rational)
      {
	ASSERT(weights.size() == points.size()/dimension);

	// Include weight information
	// First get weights in interpolation points
	shared_ptr<SplineCurve> denom = 
	  shared_ptr<SplineCurve>(new SplineCurve(basis, 
						  weights.begin(), 1,
						  false));
	vector<double> wgtval;
	denom->gridEvaluator(wgtval, par);
	size_t nmb_pnt = par.size();
	points2.reserve(nmb_pnt*(dimension+1));
	for (size_t kr=0; kr<nmb_pnt; ++kr)
	  {
	    for (int kh=0; kh<dimension; kh++)
	      points2.push_back(points[kr*dimension+kh]*wgtval[kr]);
	    points2.push_back(wgtval[kr]);
	  }
      }
    else
      points2 = points;

    // Interpolate curve
    vector<int> tg_idx;
    vector<double> tg_pnt;
    vector<double> coefs;
    SplineInterpolator interpolator;
    interpolator.setBasis(basis);
    interpolator.interpolate(par, points2, tg_idx, tg_pnt, coefs);

    // Make curve
    SplineCurve* crv = new SplineCurve(basis, coefs.begin(), dimension,
				       rational);
    return crv;
  }

} // namespace Go
