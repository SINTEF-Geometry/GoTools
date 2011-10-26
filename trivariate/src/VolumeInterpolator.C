#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/trivariate/VolumeInterpolator.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
//#include "sislP.h"

using std::vector;
using std::shared_ptr;

namespace Go
{
  SplineVolume* 
  VolumeInterpolator::regularInterpolation(const BsplineBasis& basis_u,
					   const BsplineBasis& basis_v,
					   const BsplineBasis& basis_w,
					   vector<double>& par_u,
					   vector<double>& par_v,
					   vector<double>& par_w,
					   vector<double>& points,
					   int dimension,
					   bool rational,
					   vector<double>& weights)
  {
    // Check input
    ASSERT(par_u.size()*par_v.size()*par_w.size() == points.size()/dimension);
    ASSERT(basis_u.numCoefs() == (int)par_u.size());
    ASSERT(basis_v.numCoefs() == (int)par_v.size());
    ASSERT(basis_w.numCoefs() == (int)par_w.size());

    vector<double> points2;
    if (rational)
      {
	ASSERT(weights.size() == points.size()/dimension);

	// Include weight information
	// First get weights in interpolation points
	shared_ptr<SplineVolume> denom = 
	  shared_ptr<SplineVolume>(new SplineVolume(basis_u, basis_v,
						    basis_w,
						    weights.begin(), 1,
						    false));
	vector<double> wgtval;
	denom->gridEvaluator(par_u, par_v, par_w, wgtval);
	size_t nmb_pnt = par_u.size()*par_v.size()*par_w.size();
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

    // Interpolate surfaces in the second parameter direction and
    // curves in the first parameter direction
    size_t ki, kj;
    vector<double> sf_coefs;
    vector<int> tg_idx;
    vector<double> tg_pnt;
    for (kj=0; kj<par_w.size(); ++kj)
      {
	vector<double> cv_coefs;
	for (ki=0; ki<par_v.size(); ++ki)
	  {
	    // Interpolate
	    vector<double> coefs;
	    SplineInterpolator u_interpolator;
	    vector<double> pnts;
	    pnts.insert(pnts.end(), 
			points2.begin()+(kj*par_v.size()+ki)*dimension*par_u.size(),
			points2.begin()+(kj*par_v.size()+ki+1)*dimension*par_u.size());
	    u_interpolator.setBasis(basis_u);
	    u_interpolator.interpolate(par_u, pnts,
				       tg_idx, tg_pnt, coefs);
	    cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());

// 	    vector<int> type(par_u.size(), 1);
// 	    vector<int> der(par_u.size(), 0);
// 	    int in;
// 	    int kstat = 0;
// 	    vector<double> knots;
// 	    double *coefs3;
// 	    knots.insert(knots.end(), basis_u.begin(), basis_u.end());
// 	    s1891(&par_u[0], 
// 		  &points2[(kj*par_v.size()+ki)*dimension*par_u.size()], 
// 		  dimension, par_u.size(), 1, &der[0], 1, &knots[0], &coefs3, 
// 		  &in, basis_u.order(), 0, 0, &kstat);
// 	    cv_coefs.insert(cv_coefs.end(), coefs3, coefs3+in*dimension);
// 	    free(coefs3);
	      
	  }

	// Interpolate the curves to make a surface
	SplineInterpolator v_interpolator;
	vector<double> coefs2;
	v_interpolator.setBasis(basis_v);
	v_interpolator.interpolate(par_v, cv_coefs, tg_idx, tg_pnt, coefs2);
	sf_coefs.insert(sf_coefs.end(), coefs2.begin(), coefs2.end());
      }

    // Interpolate surfaces to create volume
    SplineInterpolator w_interpolator;
    vector<double> vol_coefs;
    w_interpolator.setBasis(basis_w);
    w_interpolator.interpolate(par_w, sf_coefs, tg_idx, tg_pnt, vol_coefs);
    
    if (rational)
      {
	dimension--;
      }

    // Make surface
    SplineVolume* vol = new SplineVolume(basis_u, basis_v, basis_w,
					 vol_coefs.begin(), dimension,
					 rational);
    return vol;
  }

} // namespace Go
