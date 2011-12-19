#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/SurfaceInterpolator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineInterpolator.h"
//#include "sislP.h"

using std::vector;

namespace Go
{
  SplineSurface* 
  SurfaceInterpolator::regularInterpolation(const BsplineBasis& basis_u,
					    const BsplineBasis& basis_v,
					    vector<double>& par_u,
					    vector<double>& par_v,
					    vector<double>& points,
					    int dimension,
					    bool rational,
					    vector<double>& weights)
  {
    // Check input
    ASSERT(par_u.size()*par_v.size() == points.size()/dimension);
    ASSERT(basis_u.numCoefs() == (int)par_u.size());
    ASSERT(basis_v.numCoefs() == (int)par_v.size());

    vector<double> points2;
    if (rational)
      {
	ASSERT(weights.size() == points.size()/dimension);

	// Include weight information
	// First get weights in interpolation points
	shared_ptr<SplineSurface> denom = 
	  shared_ptr<SplineSurface>(new SplineSurface(basis_u, basis_v,
						      weights.begin(), 1,
						      false));
	vector<double> wgtval;
	denom->gridEvaluator(wgtval, par_u, par_v);
	size_t nmb_pnt = par_u.size()*par_v.size();
	points2.reserve(nmb_pnt*(dimension+1));
	for (size_t kr=0; kr<nmb_pnt; ++kr)
	  {
	    for (int kh=0; kh<dimension; kh++)
	      points2.push_back(points[kr*dimension+kh]*wgtval[kr]);
	    points2.push_back(wgtval[kr]);
	  }
	dimension++;
      }
    else
      points2 = points;

    // Interpolate curves in the first parameter direction
    size_t ki;
    vector<double> cv_coefs;
    vector<int> tg_idx;
    vector<double> tg_pnt;
    for (ki=0; ki<par_v.size(); ++ki)
      {
	// Interpolate
 	vector<double> coefs;
	SplineInterpolator u_interpolator;
	vector<double> pnts;
	pnts.insert(pnts.end(), points2.begin()+ki*dimension*par_u.size(),
		    points2.begin()+(ki+1)*dimension*par_u.size());
	u_interpolator.setBasis(basis_u);
	u_interpolator.interpolate(par_u, pnts,
				   tg_idx, tg_pnt, coefs);
	cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());

// 	vector<int> type(par_u.size(), 1);
// 	vector<int> der(par_u.size(), 0);
// 	int in;
// 	int kstat = 0;
// 	vector<double> knots;
// 	double *coefs2;
// 	knots.insert(knots.end(), basis_u.begin(), basis_u.end());
// 	s1891(&par_u[0], &points2[ki*dimension*par_u.size()], dimension,
// 	      par_u.size(), 1, &der[0], 1, &knots[0], &coefs2, 
// 	      &in, basis_u.order(), 0, 0, &kstat);
// 	cv_coefs.insert(cv_coefs.end(), coefs2, coefs2+in*dimension);
// 	free(coefs2);
      }

    // Interpolate the curves to make a surface
    SplineInterpolator v_interpolator;
    vector<double> sf_coefs;
    v_interpolator.setBasis(basis_v);
    v_interpolator.interpolate(par_v, cv_coefs, tg_idx, tg_pnt, sf_coefs);

    if (rational)
      {
	dimension--;
      }

    // Make surface
    SplineSurface* surf = new SplineSurface(basis_u, basis_v,
					    sf_coefs.begin(), dimension,
					    rational);
    return surf;
  }

} // namespace Go
