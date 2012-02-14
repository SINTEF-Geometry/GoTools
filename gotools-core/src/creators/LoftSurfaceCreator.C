//===========================================================================
//
// File : LoftSurfaceCreator.C
//
//===========================================================================


#include "GoTools/creators/LoftSurfaceCreator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/HermiteInterpolator.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include <algorithm>
#include <cmath>
#include <fstream> // For debugging.

using namespace Go;
using std::vector;


namespace Go
{



//===========================================================================
SplineSurface*
LoftSurfaceCreator::loftSurface(vector<shared_ptr<SplineCurve> >::iterator
			      first_curve,
			      int nmb_crvs)
//===========================================================================
{
  vector<shared_ptr<SplineCurve> > unified_curves = unifiedCurvesCopy(first_curve, nmb_crvs);
  vector<double> params;
  // We're giving the surface a parameter domain of length 1.0 in w-direction
  makeLoftParams(unified_curves.begin(), nmb_crvs, 1.0, params);

  return loftSurfaceFromUnifiedCurves(unified_curves.begin(), params.begin(), nmb_crvs);
}



  
//===========================================================================
SplineSurface*
LoftSurfaceCreator::loftSurface(vector<shared_ptr<SplineCurve> >::iterator
			      first_curve,
			      vector<double>::iterator first_param,
			      int nmb_crvs)
//===========================================================================
{
  vector<shared_ptr<SplineCurve> > unified_curves = unifiedCurvesCopy(first_curve, nmb_crvs);
  return loftSurfaceFromUnifiedCurves(unified_curves.begin(), first_param, nmb_crvs);
}


  
//===========================================================================
SplineSurface*
LoftSurfaceCreator::loftSurface(vector<shared_ptr<SplineCurve> >& loft_curves,
				int nmb_crvs)
//===========================================================================
{
  vector<shared_ptr<SplineCurve> > unified_curves = unifiedCurvesCopy(loft_curves.begin(), nmb_crvs);
  SplineSurface* pls = CoonsPatchGen::loftSurface(unified_curves.begin(), nmb_crvs);   
  loft_curves = unified_curves;

  return pls;
}




//===========================================================================
vector<shared_ptr<SplineCurve> >
LoftSurfaceCreator::unifiedCurvesCopy(vector<shared_ptr<SplineCurve> >::iterator
				       first_curve,
				       int nmb_crvs)
//===========================================================================
{
  bool rational = false;
  for (int i = 0; i < nmb_crvs; ++i)
    if (first_curve[i]->rational())
      rational = true;

  // Create the copies
  vector<shared_ptr<SplineCurve> > unified_curves;
  for (int i = 0; i < nmb_crvs; ++i)
    {
      shared_ptr<SplineCurve> curve_copy(first_curve[i]->clone());
      if (rational)
	curve_copy->representAsRational();
      unified_curves.push_back(curve_copy);
    }

  // Reparametrize, to have the B-spline spaces living on the same intervals
  double
    avg_start = 0.0,
    avg_end   = 0.0;

  for (int i = 0; i < nmb_crvs; ++i)
    {
      avg_start += unified_curves[i]->startparam();
      avg_end   += unified_curves[i]->endparam();
    }

  avg_start /= double(nmb_crvs);
  avg_end   /= double(nmb_crvs);

  for (int i = 0; i < nmb_crvs; ++i)
    unified_curves[i]->setParameterInterval(avg_start, avg_end);


  // Put the curves into common basis.
  double tolerance = 1e-05;
  GeometryTools::unifyCurveSplineSpace(unified_curves, tolerance);

  return unified_curves;
}




//===========================================================================
void LoftSurfaceCreator::makeLoftParams(vector<shared_ptr<SplineCurve> >::const_iterator
				       first_curve,
				       int nmb_crvs, double param_length,
				       vector<double>& params)
//===========================================================================
{
  params.clear();
  bool rational = first_curve[0]->rational();

  // Compute parameterization.
  // For each adjacent pair of curves compute the minimum
  // distance between the coefficients
  int dim = first_curve[0]->dimension();
  int num_coefs = first_curve[0]->numCoefs();

  params.push_back(0.0);  // Parameter value of first curve.
  for (int i = 1; i < nmb_crvs; ++i)
    {
      vector<double>::const_iterator coefs1, coefs2;
      if (rational)
	{
	  coefs1 = first_curve[i-1]->rcoefs_begin();
	  coefs2 = first_curve[i]->rcoefs_begin();
	}
      else
	{
	  coefs1 = first_curve[i-1]->coefs_begin();
	  coefs2 = first_curve[i]->coefs_begin();
	}
      double sqdist = 0.0;

      for (int j = 0; j < num_coefs; ++j)
	{
	  double inv_w1 = 1.0;
	  double inv_w2 = 1.0;
	  if (rational)
	    {
	      inv_w1 = 1.0 / coefs1[dim];
	      inv_w2 = 1.0 / coefs2[dim];
	    }
	  for (int k = 0; k < dim; ++k, ++coefs1, ++coefs2)
	    {
	      double diff = (*coefs1)*inv_w1 - (*coefs2)*inv_w2;
	      sqdist += diff * diff;
	    }
	  if (rational)
	    {
	      ++coefs1;
	      ++coefs2;
	    }

	}

      params.push_back(params[i-1] + sqrt(sqdist));
    }

  // We make sure the parameters go from 0.0 to param_length.
  double scale = param_length / params[nmb_crvs-1];
  for (int i = 1; i < nmb_crvs; ++i)
      params[i] *= scale;

}




//===========================================================================
SplineSurface*
LoftSurfaceCreator::loftSurfaceFromUnifiedCurves(vector<shared_ptr<SplineCurve> >::iterator
						 first_curve,
						 vector<double>::iterator first_param,
						 int nmb_crvs)
//===========================================================================
{

  SplineSurface* surf = loftNonrationalSurface(first_curve, first_param, nmb_crvs);
  if (first_curve[0]->rational())
    {
      int n = surf->numCoefs_u() * surf->numCoefs_v();
      int kdim = surf->dimension() + 1;
      bool all_positive = true;
      vector<double>::const_iterator it = surf->rcoefs_begin();
      it += (kdim - 1);
      for (int i = 0; i < n; ++i)
	if (it[kdim * i] <= 0.0)
	  {
	    all_positive = false;
	    break;
	  }
      if (!all_positive)
	{
	  delete surf;
	  surf = loftRationalSurface(first_curve, first_param, nmb_crvs);
	}
    }

  return surf;
}




//===========================================================================
SplineSurface*
LoftSurfaceCreator::loftNonrationalSurface(vector<shared_ptr<SplineCurve> >::iterator
					 first_curve,
					 vector<double>::iterator first_param,
					 int nmb_crvs)
//===========================================================================
{

  // By looking at each curve as a point of dimension dim*u_coefs.size()*v_coefs.size(),
  // the interpolation method gets rather easy.

  bool rational = first_curve[0]->rational();

  vector<double> coefs, coefs_uv_loft, params(nmb_crvs);
  for (int i = 0; i < nmb_crvs; ++i)
    {
      if (rational)
	coefs.insert(coefs.end(),
		     first_curve[i]->rcoefs_begin(), first_curve[i]->rcoefs_end());
      else
	coefs.insert(coefs.end(),
		     first_curve[i]->coefs_begin(), first_curve[i]->coefs_end());
      params[i] = first_param[i];
    }

  vector<double> tangent_points;
  vector<int> cross_index;

  int order = std::min(4, nmb_crvs);
  SplineInterpolator u_interpolator;
  u_interpolator.makeBasis(params, cross_index, order);
  u_interpolator.interpolate(params, coefs, cross_index,
			     tangent_points, coefs_uv_loft);

  SplineSurface* surf =
      new SplineSurface(first_curve[0]->basis(), u_interpolator.basis(),
			coefs_uv_loft.begin(), first_curve[0]->dimension(), rational);
  return surf;
}




//===========================================================================
SplineSurface* LoftSurfaceCreator::loftRationalSurface(vector<shared_ptr<SplineCurve> >::iterator first_curve,
						    vector<double>::iterator first_param,
						    int nmb_crvs)
//===========================================================================
{
    std::cout << "\nloftRationalSurface is not yet implemented!" << std::endl;
    return 0;
}


} // end namespace Go.
