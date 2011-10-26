//===========================================================================
//
// File : CoonsPatchVolumeGen.C
//
// Created: Fri Jun 19 09:28:01 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: CoonsPatchVolumeGen.C,v 1.2 2009/07/08 06:22:18 kfp Exp $
//
// Description:
//
//===========================================================================


#include "GoTools/trivariate/CoonsPatchVolumeGen.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/GeometryTools.h"
#include <algorithm>
#include <iterator>


using namespace Go;
using std::vector;
using std::shared_ptr;


//===========================================================================
void Go::CoonsPatchVolumeGen::get_corners(shared_ptr<SplineSurface> surf, vector<Point>& pts)
//===========================================================================
{
  pts.resize(4);
  surf->point(pts[0], 0.0, 0.0);
  surf->point(pts[1], 0.0, 1.0);
  surf->point(pts[2], 1.0, 0.0);
  surf->point(pts[3], 1.0, 1.0);
}


//===========================================================================
void Go::CoonsPatchVolumeGen::push_corners(vector<Point>& pts_to, const vector<Point> pts_from,
		  int pos0, int pos1, int pos2, int pos3)
//===========================================================================
{
  pts_to[pos0] = pts_from[0];
  pts_to[pos1] = pts_from[1];
  pts_to[pos2] = pts_from[2];
  pts_to[pos3] = pts_from[3];
}

//===========================================================================
bool Go::CoonsPatchVolumeGen::edge_curves_equal(shared_ptr<SplineSurface> sf1, double par1, bool is_u_dir1,
						shared_ptr<SplineSurface> sf2, double par2, bool is_u_dir2,
						double tol)
//===========================================================================
{
  //  double tol_sq = tol * tol;
  SplineCurve *c1, *c2;
  int dim = sf1->dimension();
  c1 = sf1->constParamCurve(par1, is_u_dir1);
  c2 = sf2->constParamCurve(par2, is_u_dir2);
  double close = true;

  vector<double>::const_iterator it1 = c1->coefs_begin();
  vector<double>::const_iterator it2 = c2->coefs_begin();

  for (; it1 != c1->coefs_end() && close; )
    {
      double dist = 0.0;
      for (int i = 0; i < dim; ++i, ++it1, ++it2)
	{
	  double diff = (*it1) - (*it2);
	  dist += diff * diff;
	}
//       if (dist >= 100.0*tol_sq)
// 	std::cout << "Coons volume gen, coef dist: " << dist << std::endl;
	//close = false;
    }

  delete c1;
  delete c2;

  return (close != 0.0);
}

//===========================================================================
SplineVolume* Go::CoonsPatchVolumeGen::createCoonsPatchDirectly(const SplineSurface* surf_u_min, const SplineSurface* surf_u_max,
								const SplineSurface* surf_v_min, const SplineSurface* surf_v_max,
								const SplineSurface* surf_w_min, const SplineSurface* surf_w_max)
//===========================================================================
{
  // We need the corner points, stored in the following order (where V(u,v,w) is the final spline volume function)
  // V(0,0,0), V(0,0,1), V(0,1,0), V(0,1,1), V(1,0,0), V(1,0,1), V(1,1,0), V(1,1,1)
  Point corner[8];
  surf_u_min->point(corner[0], 0.0, 0.0);
  surf_u_min->point(corner[1], 0.0, 1.0);
  surf_u_min->point(corner[2], 1.0, 0.0);
  surf_u_min->point(corner[3], 1.0, 1.0);
  surf_u_max->point(corner[4], 0.0, 0.0);
  surf_u_max->point(corner[5], 0.0, 1.0);
  surf_u_max->point(corner[6], 1.0, 0.0);
  surf_u_max->point(corner[7], 1.0, 1.0);

  // We also need the edge curves, stored in the following order
  // V(0,0,*), V(0,1,*) V(1,0,*) V(1,1,*)
  // V(0,*,0), V(0,*,1) V(1,*,0) V(1,*,1)
  // V(*,0,0), V(*,0,1) V(*,1,0) V(*,1,1)
  shared_ptr<SplineCurve> edge[12];
  edge[0] = shared_ptr<SplineCurve>(surf_u_min->constParamCurve(0.0, false));
  edge[1] = shared_ptr<SplineCurve>(surf_u_min->constParamCurve(1.0, false));
  edge[2] = shared_ptr<SplineCurve>(surf_u_max->constParamCurve(0.0, false));
  edge[3] = shared_ptr<SplineCurve>(surf_u_max->constParamCurve(1.0, false));

  edge[4] = shared_ptr<SplineCurve>(surf_u_min->constParamCurve(0.0, true));
  edge[5] = shared_ptr<SplineCurve>(surf_u_min->constParamCurve(1.0, true));
  edge[6] = shared_ptr<SplineCurve>(surf_u_max->constParamCurve(0.0, true));
  edge[7] = shared_ptr<SplineCurve>(surf_u_max->constParamCurve(1.0, true));

  edge[8] = shared_ptr<SplineCurve>(surf_v_min->constParamCurve(0.0, true));
  edge[9] = shared_ptr<SplineCurve>(surf_v_min->constParamCurve(1.0, true));
  edge[10] = shared_ptr<SplineCurve>(surf_v_max->constParamCurve(0.0, true));
  edge[11] = shared_ptr<SplineCurve>(surf_v_max->constParamCurve(1.0, true));

  // Store the number of coefficients and the orders of
  // all parameter directions.
  int dim = edge[0]->dimension();
  int n0 = edge[8]->numCoefs();
  int n1 = edge[4]->numCoefs();
  int n2 = edge[0]->numCoefs();
  int k0 = edge[8]->order();
  int k1 = edge[4]->order();
  int k2 = edge[0]->order();
  const BsplineBasis& basu = edge[8]->basis();
  const BsplineBasis& basv = edge[4]->basis();
  const BsplineBasis& basw = edge[0]->basis();

  // Coefficients for interpolated volumes
  vector<double>
    co_sf_u(k0*n1*n2*dim),    // Interpolated in u parameter direction
    co_sf_v(n0*k1*n2*dim),    // Interpolated in v parameter direction
    co_sf_w(n0*n1*k2*dim),    // Interpolated in w parameter direction
    co_crv_uv(k0*k1*n2*dim),  // Interpolated in u+v parameter direction
    co_crv_uw(k0*n1*k2*dim),  // Interpolated in u+w parameter direction
    co_crv_vw(n0*k1*k2*dim),  // Interpolated in v+w parameter direction
    co_pts(k0*k1*k2*dim);     // Interpolated in u+v+w parameter direction

  // Minimal bases
  vector<double> knu(2*k0, 0.0);
  vector<double> knv(2*k1, 0.0);
  vector<double> knw(2*k2, 0.0);
  fill(knu.begin()+k0, knu.end(), 1.0);
  fill(knv.begin()+k1, knv.end(), 1.0);
  fill(knw.begin()+k2, knw.end(), 1.0);
  BsplineBasis minbasu(k0, k0, &(knu[0]));
  BsplineBasis minbasv(k1, k1, &(knv[0]));
  BsplineBasis minbasw(k2, k2, &(knw[0]));

  // Some surface coefficient iterators
  vector<double>::const_iterator surf_u_min_coef = surf_u_min->coefs_begin();
  vector<double>::const_iterator surf_u_max_coef = surf_u_max->coefs_begin();
  vector<double>::const_iterator surf_v_min_coef = surf_v_min->coefs_begin();
  vector<double>::const_iterator surf_v_max_coef = surf_v_max->coefs_begin();
  vector<double>::const_iterator surf_w_min_coef = surf_w_min->coefs_begin();
  vector<double>::const_iterator surf_w_max_coef = surf_w_max->coefs_begin();

  // Interpolation in u-direction
  for (int i = 0; i < k0; ++i)
    {
      double fac = minbasu.grevilleParameter(i);
      for (int j = 0; j < n1; ++j)
	for (int k = 0; k < n2; ++k)
	  for (int dd = 0; dd < dim; ++dd)
	    co_sf_u[k*n1*k0*dim + j*k0*dim + i*dim + dd]
	      = (1.0-fac) * surf_u_min_coef[k*n1*dim + j*dim + dd]
	      + fac * surf_u_max_coef[k*n1*dim + j*dim + dd];
    }

  // Interpolation in v-direction
  for (int j = 0; j< k1; ++j)
    {
      double fac = minbasv.grevilleParameter(j);
      for (int i = 0; i < n0; ++i)
	for (int k = 0; k < n2; ++k)
	  for (int dd = 0; dd < dim; ++dd)
	    co_sf_v[k*k1*n0*dim + j*n0*dim + i*dim + dd]
	      = (1.0-fac) * surf_v_min_coef[k*n0*dim + i*dim + dd]
	      + fac * surf_v_max_coef[k*n0*dim + i*dim + dd];
    }

  // Interpolation in w-direction
  for (int k = 0; k < k2; ++k)
    {
      double fac = minbasw.grevilleParameter(k);
      for (int i = 0; i < n0; ++i)
	for (int j = 0; j < n1; ++j)
	  for (int dd = 0; dd < dim; ++dd)
	    co_sf_w[k*n1*n0*dim + j*n0*dim + i*dim + dd]
	      = (1.0-fac) * surf_w_min_coef[j*n0*dim + i*dim + dd]
	      + fac * surf_w_max_coef[j*n0*dim + i*dim + dd];
    }

  // Interpolation in uv-direction
  for (int i = 0; i < k0; ++i)
    {
      double fac0 = minbasu.grevilleParameter(i);
      for (int j = 0; j < k1; ++j)
	{
	  double fac1 = minbasv.grevilleParameter(j);
	  for (int k = 0; k < n2; ++k)
	    for (int dd = 0; dd < dim; ++dd)
	      co_crv_uv[k*k1*k0*dim + j*k0*dim + i*dim + dd]
		= (1.0-fac0) * (1.0-fac1) * (edge[0]->coefs_begin()[k*dim + dd])
		+ (1.0-fac0) * fac1 * (edge[1]->coefs_begin()[k*dim + dd])
		+ fac0 * (1.0-fac1) * (edge[2]->coefs_begin()[k*dim + dd])
		+ fac0 * fac1 * (edge[3]->coefs_begin()[k*dim + dd]);
	}
    }

  // Interpolation in uw-direction
  for (int i = 0; i < k0; ++i)
    {
      double fac0 = minbasu.grevilleParameter(i);
      for (int k = 0; k < k2; ++k)
	{
	  double fac1 = minbasw.grevilleParameter(k);
	  for (int j = 0; j < n1; ++j)
	    for (int dd = 0; dd < dim; ++dd)
	      co_crv_uw[k*n1*k0*dim + j*k0*dim + i*dim + dd]
		= (1.0-fac0) * (1.0-fac1) * (edge[4]->coefs_begin()[j*dim + dd])
		+ (1.0-fac0) * fac1 * (edge[5]->coefs_begin()[j*dim + dd])
		+ fac0 * (1.0-fac1) * (edge[6]->coefs_begin()[j*dim + dd])
		+ fac0 * fac1 * (edge[7]->coefs_begin()[j*dim + dd]);
	}
    }

  // Interpolation in vw-direction
  for (int j = 0; j < k1; ++j)
    {
      double fac0 = minbasv.grevilleParameter(j);
      for (int k = 0; k < k2; ++k)
	{
	  double fac1 = minbasw.grevilleParameter(k);
	  for (int i = 0; i < n0; ++i)
	    for (int dd = 0; dd < dim; ++dd)
	      co_crv_vw[k*k1*n0*dim + j*n0*dim + i*dim + dd]
		= (1.0-fac0) * (1.0-fac1) * (edge[8]->coefs_begin()[i*dim + dd])
		+ (1.0-fac0) * fac1 * (edge[9]->coefs_begin()[i*dim + dd])
		+ fac0 * (1.0-fac1) * (edge[10]->coefs_begin()[i*dim + dd])
		+ fac0 * fac1 * (edge[11]->coefs_begin()[i*dim + dd]);
	}
    }

  // Interpolation in uvw-direction
  for (int i = 0; i < k0; ++i)
    {
      double fac0 = minbasu.grevilleParameter(i);
      for (int j = 0; j < k1; ++j)
	{
	  double fac1 = minbasv.grevilleParameter(j);
	  for (int k = 0; k < k2; ++k)
	    {
	      double fac2 = minbasw.grevilleParameter(k);
	      for (int dd = 0; dd < dim; ++dd)
		co_pts[k*k1*k0*dim + j*k0*dim + i*dim + dd]
		  = (1.0-fac0) * (1.0-fac1) * (1.0-fac2) * corner[0][dd]
		  + (1.0-fac0) * (1.0-fac1) * fac2 * corner[1][dd]
		  + (1.0-fac0) * fac1 * (1.0-fac2) * corner[2][dd]
		  + (1.0-fac0) * fac1 * fac2 * corner[3][dd]
		  + fac0 * (1.0-fac1) * (1.0-fac2) * corner[4][dd]
		  + fac0 * (1.0-fac1) * fac2 * corner[5][dd]
		  + fac0 * fac1 * (1.0-fac2) * corner[6][dd]
		  + fac0 * fac1 * fac2 * corner[7][dd];
	    }
	}
    }

  // Put all volumes on the same knot vectors
  vector<double> newu, newv, neww;
  set_difference(basu.begin(), basu.end(),
		 minbasu.begin(), minbasu.end(),
		 back_inserter(newu));
  set_difference(basv.begin(), basv.end(),
		 minbasv.begin(), minbasv.end(),
		 back_inserter(newv));
  set_difference(basw.begin(), basw.end(),
		 minbasw.begin(), minbasw.end(),
		 back_inserter(neww));

  // Create volumes for adding and subtracting coefficients
  SplineVolume vol_u(minbasu, basv, basw, co_sf_u.begin(), dim);
  SplineVolume vol_v(basu, minbasv, basw, co_sf_v.begin(), dim);
  SplineVolume vol_w(basu, basv, minbasw, co_sf_w.begin(), dim);

  SplineVolume vol_uv(minbasu, minbasv, basw, co_crv_uv.begin(), dim);
  SplineVolume vol_uw(minbasu, basv, minbasw, co_crv_uw.begin(), dim);
  SplineVolume vol_vw(basu, minbasv, minbasw, co_crv_vw.begin(), dim);

  SplineVolume vol_uvw(minbasu, minbasv, minbasw, co_pts.begin(), dim);

  vol_u.insertKnot(0, newu);
  vol_v.insertKnot(1, newv);
  vol_w.insertKnot(2, neww);

  vol_uv.insertKnot(0, newu);
  vol_uv.insertKnot(1, newv);
  vol_uw.insertKnot(0, newu);
  vol_uw.insertKnot(2, neww);
  vol_vw.insertKnot(1, newv);
  vol_vw.insertKnot(2, neww);

  vol_uvw.insertKnot(0, newu);
  vol_uvw.insertKnot(1, newv);
  vol_uvw.insertKnot(2, neww);

  // Now we add and subtract the other coefficients into vol_u
  for (int i = 0; i < n0*n1*n2*dim; ++i)
    {
      vol_u.coefs_begin()[i] += vol_v.coefs_begin()[i];
      vol_u.coefs_begin()[i] += vol_w.coefs_begin()[i];

      vol_u.coefs_begin()[i] -= vol_uv.coefs_begin()[i];
      vol_u.coefs_begin()[i] -= vol_uw.coefs_begin()[i];
      vol_u.coefs_begin()[i] -= vol_vw.coefs_begin()[i];

      vol_u.coefs_begin()[i] += vol_uvw.coefs_begin()[i];
    }

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
  return dynamic_cast<SplineSurface*>(vol_u.clone());
#else
  return vol_u.clone();
#endif

}




//===========================================================================
SplineVolume* Go::CoonsPatchVolumeGen::createCoonsPatch(const SplineSurface* surf_u_min, const SplineSurface* surf_u_max,
							const SplineSurface* surf_v_min, const SplineSurface* surf_v_max,
							const SplineSurface* surf_w_min, const SplineSurface* surf_w_max,
							double tol)
//===========================================================================
{
  double tol_sq = tol * tol;

  // Check surfaces are non-rational and in same dimensional space
  ALWAYS_ERROR_IF(surf_u_min->rational(),
		  "Surfaces must be non-rational");
  ALWAYS_ERROR_IF(surf_u_max->rational(),
		  "Surfaces must be non-rational");
  ALWAYS_ERROR_IF(surf_v_min->rational(),
		  "Surfaces must be non-rational");
  ALWAYS_ERROR_IF(surf_v_max->rational(),
		  "Surfaces must be non-rational");
  ALWAYS_ERROR_IF(surf_w_min->rational(),
		  "Surfaces must be non-rational");
  ALWAYS_ERROR_IF(surf_w_max->rational(),
		  "Surfaces must be non-rational");

  int dim = surf_u_min->dimension();
  ALWAYS_ERROR_IF(surf_u_max->dimension() != dim,
		  "Dimension mismatch.");
  ALWAYS_ERROR_IF(surf_v_min->dimension() != dim,
		  "Dimension mismatch.");
  ALWAYS_ERROR_IF(surf_v_max->dimension() != dim,
		  "Dimension mismatch.");
  ALWAYS_ERROR_IF(surf_w_min->dimension() != dim,
		  "Dimension mismatch.");
  ALWAYS_ERROR_IF(surf_w_max->dimension() != dim,
		  "Dimension mismatch.");

  // Use copies to avoid damaging the parameters.
  shared_ptr<SplineSurface> sf_u_min(surf_u_min->clone());
  shared_ptr<SplineSurface> sf_u_max(surf_u_max->clone());
  shared_ptr<SplineSurface> sf_v_min(surf_v_min->clone());
  shared_ptr<SplineSurface> sf_v_max(surf_v_max->clone());
  shared_ptr<SplineSurface> sf_w_min(surf_w_min->clone());
  shared_ptr<SplineSurface> sf_w_max(surf_w_max->clone());

  // Rescale, so we always work with unit interval for bases
  sf_u_min->setParameterDomain(0.0, 1.0, 0.0, 1.0);
  sf_u_max->setParameterDomain(0.0, 1.0, 0.0, 1.0);
  sf_v_min->setParameterDomain(0.0, 1.0, 0.0, 1.0);
  sf_v_max->setParameterDomain(0.0, 1.0, 0.0, 1.0);
  sf_w_min->setParameterDomain(0.0, 1.0, 0.0, 1.0);
  sf_w_max->setParameterDomain(0.0, 1.0, 0.0, 1.0);

  // Fetch corner points
  vector<Point> sf_u_min_pts, sf_u_max_pts, sf_v_min_pts, sf_v_max_pts, sf_w_min_pts, sf_w_max_pts;
  get_corners(sf_u_min, sf_u_min_pts);
  get_corners(sf_u_max, sf_u_max_pts);
  get_corners(sf_v_min, sf_v_min_pts);
  get_corners(sf_v_max, sf_v_max_pts);
  get_corners(sf_w_min, sf_w_min_pts);
  get_corners(sf_w_max, sf_w_max_pts);

  // Find a common corner for surf_u_min, surf_v_min and surf_w_min
  int posu, posv, posw;

  for (posu = 0; posu < 4; ++posu)
    {
      for (posv = 0; posv < 4; ++posv)
	if (sf_u_min_pts[posu].dist2(sf_v_min_pts[posv]) < tol_sq)
	  {
	    for (posw = 0; posw < 4; ++posw)
	      if (sf_u_min_pts[posu].dist2(sf_w_min_pts[posw]) < tol_sq)
		break;
	    if (posw < 4)
	      break;
	  }
      if (posv < 4)
	break;
    }

  ALWAYS_ERROR_IF(posu == 4,
		  "Faces do not have common corner");

  // Rearrange to have common corner in position (0, 0) for all three surfaces
  if (posu != 0)
    {
      if ((posu >> 1) == 1)
	sf_u_min->reverseParameterDirection(true);
      if ((posu & 1) == 1)
	sf_u_min->reverseParameterDirection(false);
      get_corners(sf_u_min, sf_u_min_pts);
    }
  if (posv != 0)
    {
      if ((posv >> 1) == 1)
	sf_v_min->reverseParameterDirection(true);
      if ((posv & 1) == 1)
	sf_v_min->reverseParameterDirection(false);
      get_corners(sf_v_min, sf_v_min_pts);
    }
  if (posw != 0)
    {
      if ((posw >> 1) == 1)
	sf_w_min->reverseParameterDirection(true);
      if ((posw & 1) == 1)
	sf_w_min->reverseParameterDirection(false);
      get_corners(sf_w_min, sf_w_min_pts);
    }

  // Find the other common corners of sf_u_min, sf_v_min and sf_w_min
  if (sf_u_min_pts[1].dist2(sf_v_min_pts[2]) < tol_sq)
    {
      sf_v_min->swapParameterDirection();
      get_corners(sf_v_min, sf_v_min_pts);
    }
  else if (sf_u_min_pts[2].dist2(sf_v_min_pts[1]) < tol_sq)
    {
      sf_u_min->swapParameterDirection();
      get_corners(sf_u_min, sf_u_min_pts);
    }
  else if (sf_u_min_pts[2].dist2(sf_v_min_pts[2]) < tol_sq)
    {
      sf_u_min->swapParameterDirection();
      sf_v_min->swapParameterDirection();
      get_corners(sf_u_min, sf_u_min_pts);
      get_corners(sf_v_min, sf_v_min_pts);
    }

  if (sf_u_min_pts[2].dist2(sf_w_min_pts[2]) < tol_sq)
    {
      sf_w_min->swapParameterDirection();
      get_corners(sf_w_min, sf_w_min_pts);
    }

  // Rearrange surf_u_max
  for (posu = 0; posu < 4; ++posu)
    if (sf_v_min_pts[2].dist2(sf_u_max_pts[posu]) < tol_sq)
      break;
  ALWAYS_ERROR_IF(posu == 4,
		  "Faces do not have common corner");
  if (posu != 0)
    {
      if ((posu >> 1) == 1)
	sf_u_max->reverseParameterDirection(true);
      if ((posu & 1) == 1)
	sf_u_max->reverseParameterDirection(false);
      get_corners(sf_u_max, sf_u_max_pts);
     }
  if (sf_v_min_pts[3].dist2(sf_u_max_pts[2]) < tol_sq)
    {
      sf_u_max->swapParameterDirection();
      get_corners(sf_u_max, sf_u_max_pts);
    }

  // Rearrange surf_v_max
  for (posv = 0; posv < 4; ++posv)
    if (sf_u_min_pts[2].dist2(sf_v_max_pts[posv]) < tol_sq)
      break;
  ALWAYS_ERROR_IF(posv == 4,
		  "Faces do not have common corner");
  if (posv != 0)
    {
      if ((posv >> 1) == 1)
	sf_v_max->reverseParameterDirection(true);
      if ((posv & 1) == 1)
	sf_v_max->reverseParameterDirection(false);
      get_corners(sf_v_max, sf_v_max_pts);
     }
  if (sf_u_min_pts[3].dist2(sf_v_max_pts[2]) < tol_sq)
    {
      sf_v_max->swapParameterDirection();
      get_corners(sf_v_max, sf_v_max_pts);
    }

  // Rearrange surf_w_max
  for (posw = 0; posw < 4; ++posw)
    if (sf_u_min_pts[1].dist2(sf_w_max_pts[posw]) < tol_sq)
      break;
  ALWAYS_ERROR_IF(posw == 4,
		  "Faces do not have common corner");
  if (posw != 0)
    {
      if ((posw >> 1) == 1)
	sf_w_max->reverseParameterDirection(true);
      if ((posw & 1) == 1)
	sf_w_max->reverseParameterDirection(false);
      get_corners(sf_w_max, sf_w_max_pts);
     }
  if (sf_u_min_pts[3].dist2(sf_w_max_pts[2]) < tol_sq)
    {
      sf_w_max->swapParameterDirection();
      get_corners(sf_w_max, sf_w_max_pts);
    }

  // All faces are done with parameter reversing and swapping
  // Now we test if all eight corners coincide
  vector<Point> u_pts(8), v_pts(8), w_pts(8);
  push_corners(u_pts, sf_u_min_pts, 0, 1, 2, 3);
  push_corners(u_pts, sf_u_max_pts, 4, 5, 6, 7);
  push_corners(v_pts, sf_v_min_pts, 0, 1, 4, 5);
  push_corners(v_pts, sf_v_max_pts, 2, 3, 6, 7);
  push_corners(w_pts, sf_w_min_pts, 0, 2, 4, 6);
  push_corners(w_pts, sf_w_max_pts, 1, 3, 5, 7);

  for (int i = 0; i < 8; ++i)
    ALWAYS_ERROR_IF(u_pts[i].dist2(v_pts[i]) >= tol ||
		    u_pts[i].dist2(w_pts[i]) >= tol ||
		    v_pts[i].dist2(w_pts[i]) >= tol,
		  "Faces do not have common corner");

  // Raise orders
  /*
  int maxorder_u = max(max(sf_v_min->order_u(), sf_v_max->order_u()),
		       max(sf_w_min->order_u(), sf_w_max->order_u()));
  int maxorder_v = max(max(sf_u_min->order_u(), sf_u_max->order_u()),
		       max(sf_w_min->order_v(), sf_w_max->order_v()));
  int maxorder_w = max(max(sf_u_min->order_v(), sf_u_max->order_v()),
		       max(sf_v_min->order_v(), sf_v_max->order_v()));

  if (sf_u_min->order_u() < maxorder_v || sf_u_min->order_v() < maxorder_w )
    sf_u_min->raiseOrder(maxorder_v - sf_u_min->order_u(), maxorder_w - sf_u_min->order_v()); 
  if (sf_u_max->order_u() < maxorder_v || sf_u_max->order_v() < maxorder_w )
    sf_u_max->raiseOrder(maxorder_v - sf_u_max->order_u(), maxorder_w - sf_u_max->order_v()); 

  if (sf_v_min->order_u() < maxorder_u || sf_v_min->order_v() < maxorder_w )
    sf_v_min->raiseOrder(maxorder_u - sf_v_min->order_u(), maxorder_w - sf_v_min->order_v()); 
  if (sf_v_max->order_u() < maxorder_u || sf_v_max->order_v() < maxorder_w )
    sf_v_max->raiseOrder(maxorder_u - sf_v_max->order_u(), maxorder_w - sf_v_max->order_v()); 

  if (sf_w_min->order_u() < maxorder_u || sf_w_min->order_v() < maxorder_v )
    sf_w_min->raiseOrder(maxorder_u - sf_w_min->order_u(), maxorder_v - sf_w_min->order_v()); 
  if (sf_w_max->order_u() < maxorder_u || sf_w_max->order_v() < maxorder_v )
    sf_w_max->raiseOrder(maxorder_u - sf_w_max->order_u(), maxorder_v - sf_w_max->order_v()); 
  */

  // Unify knot vectors
  double ptol = std::min(tol, 1.0e-6);
  vector<shared_ptr<SplineSurface> > unif_sfs;

  unif_sfs.resize(0);
  unif_sfs.push_back(sf_u_min);
  unif_sfs.push_back(sf_u_max);
  unif_sfs.push_back(sf_v_min);
  unif_sfs.push_back(sf_v_max);
  unifySurfaceSplineSpaceOneDir(unif_sfs, ptol, false);
  sf_u_min = unif_sfs[0];
  sf_u_max = unif_sfs[1];
  sf_v_min = unif_sfs[2];
  sf_v_max = unif_sfs[3];

  unif_sfs.resize(0);
  sf_w_min->swapParameterDirection();
  sf_w_max->swapParameterDirection();
  unif_sfs.push_back(sf_u_min);
  unif_sfs.push_back(sf_u_max);
  unif_sfs.push_back(sf_w_min);
  unif_sfs.push_back(sf_w_max);
  unifySurfaceSplineSpaceOneDir(unif_sfs, ptol, true);
  sf_u_min = unif_sfs[0];
  sf_u_max = unif_sfs[1];
  sf_w_min = unif_sfs[2];
  sf_w_max = unif_sfs[3];
  sf_w_min->swapParameterDirection();
  sf_w_max->swapParameterDirection();

  unif_sfs.resize(0);
  unif_sfs.push_back(sf_v_min);
  unif_sfs.push_back(sf_v_max);
  unif_sfs.push_back(sf_w_min);
  unif_sfs.push_back(sf_w_max);
  unifySurfaceSplineSpaceOneDir(unif_sfs, ptol, true);
  sf_v_min = unif_sfs[0];
  sf_v_max = unif_sfs[1];
  sf_w_min = unif_sfs[2];
  sf_w_max = unif_sfs[3];

  // Test if edge curves are equal
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_u_min, 0.0, false, sf_v_min, 0.0, false, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_u_min, 1.0, false, sf_v_max, 0.0, false, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_u_max, 0.0, false, sf_v_min, 1.0, false, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_u_max, 1.0, false, sf_v_max, 1.0, false, tol),
		  "Edge curves are different");

  ALWAYS_ERROR_IF(!edge_curves_equal(sf_u_min, 0.0, true, sf_w_min, 0.0, false, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_u_min, 1.0, true, sf_w_max, 0.0, false, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_u_max, 0.0, true, sf_w_min, 1.0, false, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_u_max, 1.0, true, sf_w_max, 1.0, false, tol),
		  "Edge curves are different");

  ALWAYS_ERROR_IF(!edge_curves_equal(sf_v_min, 0.0, true, sf_w_min, 0.0, true, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_v_min, 1.0, true, sf_w_max, 0.0, true, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_v_max, 0.0, true, sf_w_min, 1.0, true, tol),
		  "Edge curves are different");
  ALWAYS_ERROR_IF(!edge_curves_equal(sf_v_max, 1.0, true, sf_w_max, 1.0, true, tol),
		  "Edge curves are different");

  return createCoonsPatchDirectly(sf_u_min.get(), sf_u_max.get(), sf_v_min.get(), sf_v_max.get(), sf_w_min.get(), sf_w_max.get());
}
