/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/trivariate/CoonsPatchVolumeGen.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/GeometryTools.h"
#include <algorithm>
#include <iterator>
#include <fstream>

//#define DEBUG

using namespace Go;
using std::vector;
using std::pair;
using std::make_pair;


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

#ifdef DEBUG
  std::ofstream of0("coons_in.g2");
  surf_u_min->writeStandardHeader(of0);
  surf_u_min->write(of0);
  surf_u_max->writeStandardHeader(of0);
  surf_u_max->write(of0);
  surf_v_min->writeStandardHeader(of0);
  surf_v_min->write(of0);
  surf_v_max->writeStandardHeader(of0);
  surf_v_max->write(of0);
  surf_w_min->writeStandardHeader(of0);
  surf_w_min->write(of0);
  surf_w_max->writeStandardHeader(of0);
  surf_w_max->write(of0);
  for (int ka=0; ka<12; ++ka)
    {
      edge[ka]->writeStandardHeader(of0);
      edge[ka]->write(of0);
    }
  of0 << "400 1 0 4 255 0 0 255" << std::endl;
  of0 << "8" << std::endl;
  for (int ka=0; ka<8; ++ka)
    of0 << corner[ka] << std::endl;
#endif

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

#ifdef DEBUG
  std::ofstream of("bool_vol.g2");
  vol_u.writeStandardHeader(of);
  vol_u.write(of);
  vol_v.writeStandardHeader(of);
  vol_v.write(of);
  vol_w.writeStandardHeader(of);
  vol_w.write(of);
  vol_uv.writeStandardHeader(of);
  vol_uv.write(of);
  vol_uw.writeStandardHeader(of);
  vol_uw.write(of);
  vol_vw.writeStandardHeader(of);
  vol_vw.write(of);
  vol_uvw.writeStandardHeader(of);
  vol_uvw.write(of);
#endif

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

#ifdef DEBUG
  std::ofstream of2("final_vol.g2");
  vol_u.writeStandardHeader(of2);
  vol_u.write(of2);
#endif

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
  //double tol_sq = tol * tol;

  // Check surfaces are non-rational and in same dimensional space
  if (surf_u_min->rational())
    THROW("Surfaces must be non-rational");
  if (surf_u_max->rational())
    THROW("Surfaces must be non-rational");
  if (surf_v_min->rational())
    THROW("Surfaces must be non-rational");
  if (surf_v_max->rational())
    THROW("Surfaces must be non-rational");
  if (surf_w_min->rational())
    THROW("Surfaces must be non-rational");
  if (surf_w_max->rational())
    THROW("Surfaces must be non-rational");

  int dim = surf_u_min->dimension();
  if (surf_u_max->dimension() != dim)
    THROW("Dimension mismatch.");
  if (surf_v_min->dimension() != dim)
    THROW("Dimension mismatch.");
  if (surf_v_max->dimension() != dim)
    THROW("Dimension mismatch.");
  if (surf_w_min->dimension() != dim)
    THROW("Dimension mismatch.");
  if (surf_w_max->dimension() != dim)
    THROW("Dimension mismatch.");

  // Use copies to avoid damaging the parameters.
  vector<shared_ptr<SplineSurface> > sfs(6);
  sfs[0] = shared_ptr<SplineSurface>(surf_u_min->clone());
  sfs[1] = shared_ptr<SplineSurface>(surf_u_max->clone());
  sfs[2] = shared_ptr<SplineSurface>(surf_v_min->clone());
  sfs[3] = shared_ptr<SplineSurface>(surf_v_max->clone());
  sfs[4] = shared_ptr<SplineSurface>(surf_w_min->clone());
  sfs[5] = shared_ptr<SplineSurface>(surf_w_max->clone());

#ifdef DEBUG
  std::ofstream of0("sf_in.g2");
  for (size_t ka=0; ka<sfs.size(); ++ka)
    {
      sfs[ka]->writeStandardHeader(of0);
      sfs[ka]->write(of0);
      of0 << "400 1 0 4 255 0 0 255 \n";
      of0 << "1 \n";
      of0 << sfs[ka]->ParamSurface::point(sfs[ka]->startparam_u(), sfs[ka]->startparam_v()) << "\n";
      of0 << "400 1 0 4 0 255 0 255 \n";
      of0 << "1 \n";
      of0 << sfs[ka]->ParamSurface::point(sfs[ka]->endparam_u(), sfs[ka]->startparam_v()) << "\n";
    }
  of0 << std::endl;
#endif

  // Rescale, so we always work with unit interval for bases
  int ki, kj, kr;
  for (ki=0; ki<6; ++ki)
    sfs[ki]->setParameterDomain(0.0, 1.0, 0.0, 1.0);

  // Preprocess degenerate surfaces (triangular) to have pairwise consistence
  int par;
  vector<Point> degen_pts;
  for (par=0; par<6; par+=2)
    {
      bool b1, r1, t1, l1, b2, r2, t2, l2;
      bool deg1 = sfs[par]->isDegenerate(b1, r1, t1, l1, tol);
      bool deg2 = sfs[par+1]->isDegenerate(b2, r2, t2, l2, tol);
      if (deg1)
	{
	  RectDomain dom = sfs[par]->containingDomain();
	  double upar = l1 ? dom.umin() : 
	    (r1 ? dom.umax() : 0.5*(dom.umin()+dom.umax()));
	  double vpar = b1 ? dom.vmin() :
	    (t1 ? dom.vmax() : 0.5*(dom.vmin()+dom.vmax()));
	  Point deg = sfs[par]->ParamSurface::point(upar,vpar);
	  degen_pts.push_back(deg);
	}
      if (deg2)
	{
	  RectDomain dom = sfs[par+1]->containingDomain();
	  double upar = l2 ? dom.umin() : 
	    (r2 ? dom.umax() : 0.5*(dom.umin()+dom.umax()));
	  double vpar = b2 ? dom.vmin() :
	    (t2 ? dom.vmax() : 0.5*(dom.vmin()+dom.vmax()));
	  Point deg = sfs[par+1]->ParamSurface::point(upar,vpar);
	  degen_pts.push_back(deg);
	}
      if (deg1 && deg2)
	{
	  if (((b1 && r2) && (!b2 && !r1)) || ((b2 && r1) && (!b1 && !r2)) || 
	      ((b1 && l2) && (!b2 && !l1)) || ((b2 && l1) && (!b1 && !l2)) ||
	      ((t1 && l2) && (!t2 && !l1)) || ((t2 && l1) && (!t1 && !l2)) ||
	      ((t1 && r2) && (!t2 && !r1)) || ((t2 && r1) && (!t1 && !r2)))
	    {
	      sfs[par+1]->swapParameterDirection();
	    } 
	  deg1 = sfs[par]->isDegenerate(b1, r1, t1, l1, tol);
	  deg2 = sfs[par+1]->isDegenerate(b2, r2, t2, l2, tol);
	  if (((b1 && !b2) && (t2 && !t1)) ||
	      ((b2 && !b1) && (t1 && !t2)))
	    sfs[par+1]->reverseParameterDirection(false);
	  if (((r1 && !r2) && (l2 && !l1)) ||
	      ((r2 && !r1) && (l1 && !l2)))
	    sfs[par+1]->reverseParameterDirection(true);
	}
    }

  for (int ka=0; ka<2; ++ka)
    {
  // Fetch corner points
  vector<vector<Point> > sf_pts(6);
  for (ki=0; ki<6; ++ki)
    get_corners(sfs[ki], sf_pts[ki]);

  // Identify volume corners in the following sequence: 
  // V(0,0,0), V(0,0,1), V(0,1,0), V(0,1,1), V(1,0,0), V(1,0,1), V(1,1,0), V(1,1,1)
  vector<Point> Vc(8);
  
  // Index of volume corners with respect to side surface (umin, umax, vmin, vmax, wmin, wmax)
  int ixc[6][4] = {{0,1,2,3}, {4,5,6,7}, {0,1,4,5}, {2,3,6,7}, {0,2,4,6}, {1,3,5,7}};

  // TEST
  vector<Point> Vc1(8);
  vector<Point> Vc2(8);
  vector<Point> Vc3(8);
  for (kr=0; kr<2; ++kr)
    {
      for (ki=0; ki<4; ++ki)
	Vc1[ixc[kr][ki]] = sf_pts[kr][ki];
      for (ki=0; ki<4; ++ki)
	Vc2[ixc[2+kr][ki]] = sf_pts[2+kr][ki];
      for (ki=0; ki<4; ++ki)
	Vc3[ixc[4+kr][ki]] = sf_pts[4+kr][ki];
    }
  
  int xs[8][4] = {{0,1,2,3}, {2,3,0,1}, {1,0,3,2}, {0,2,1,3}, {2,0,3,1}, {3,2,1,0}, {1,3,0,2}, {3,1,2,0}};
  int ixd[6][2] = {{0,1}, {0,2}, {0,2}, {0,1}, {0,1}, {0,2}};

  for (par=0; par<6; par+=2)
    {
      // Initial suggestion
      for (kr=0; kr<2; ++kr)
	for (ki=0; ki<4; ++ki)
	  Vc[ixc[par+kr][ki]] = sf_pts[par+kr][ki];

#ifdef DEBUG
      std::ofstream of1("lines1.g2");
      for (int ka=0; ka<4; ++ka)
	{
	  of1 << "410 1 0 4 255 0 0 255" << std::endl;
	  of1 << "1" << std::endl;
	  of1 << Vc[ixc[0][ka]] << "  " << Vc[ixc[1][ka]] << std::endl;
	}
      for (int ka=0; ka<4; ++ka)
	{
	  of1 << "410 1 0 4 0 255 0 255" << std::endl;
	  of1 << "1" << std::endl;
	  of1 << Vc[ixc[2][ka]] << "  " << Vc[ixc[3][ka]] << std::endl;
	}
      for (int ka=0; ka<4; ++ka)
	{
	  of1 << "410 1 0 4 55 100 100 255" << std::endl;
	  of1 << "1" << std::endl;
	  of1 << Vc[ixc[4][ka]] << "  " << Vc[ixc[5][ka]] << std::endl;
	}
#endif

      // Check if any of the surfaces must be turned or swapped
      // First check with min surfaces in the other parameter directions
      int par2 = (par == 4) ? 0 : par+2;
      int par3 = (par == 0) ? 4 : par-2;
      for (kr=0; kr<2; ++kr)
	{
	  int nmb = 0;
	  int ixm = 0;
	  int ixs = -1;
	  for (ki=0; ki<2; ++ki)
	    {
	      for (kj=0; kj<4; ++kj)
		{
		  if (Vc[ixc[par+kr][ixd[par][ki]]].dist(sf_pts[par2][kj]) < tol &&
		      kj != ixs)
		    {
		      nmb++;
		      ixm = ki;
		      ixs = kj;
		      break;
		    }
		}
	    }
	  if (nmb == 0)
	    {
	      sfs[par+kr]->reverseParameterDirection(true);
	    }
	  else if (nmb == 1)
	    {
	      size_t ka;
	      for (ka=0; ka<degen_pts.size(); ++ka)
		if (Vc[ixc[par+kr][ixd[par][ixm]]].dist(degen_pts[ka]) < tol)
		  //Vc[ixc[par+kr][ixd[par][1]]].dist(degen_pts[ka]) < tol)
		  break;
	      if (ka == degen_pts.size())
		sfs[par+kr]->swapParameterDirection();
	    }

	  get_corners(sfs[par+kr], sf_pts[par+kr]);
	  for (ki=0; ki<4; ++ki)
	    Vc[ixc[par+kr][ki]] = sf_pts[par+kr][ki];
	}

#ifdef DEBUG
      std::ofstream of2("lines2.g2");
      for (int ka=0; ka<4; ++ka)
	{
	  of2 << "410 1 0 4 255 0 0 255" << std::endl;
	  of2 << "1" << std::endl;
	  of2 << Vc[ixc[0][ka]] << "  " << Vc[ixc[1][ka]] << std::endl;
	}
      for (int ka=0; ka<4; ++ka)
	{
	  of2 << "410 1 0 4 0 255 0 255" << std::endl;
	  of2 << "1" << std::endl;
	  of2 << Vc[ixc[2][ka]] << "  " << Vc[ixc[3][ka]] << std::endl;
	}
      for (int ka=0; ka<4; ++ka)
	{
	  of2 << "410 1 0 4 55 100 100 255" << std::endl;
	  of2 << "1" << std::endl;
	  of2 << Vc[ixc[4][ka]] << "  " << Vc[ixc[5][ka]] << std::endl;
	}
#endif

      // Check with surface w-min
      for (kr=0; kr<2; ++kr)
	{
	  int nmb = 0;
	  int ixm = 0;
	  int ixs = -1;
	  for (ki=0; ki<2; ++ki)
	    {
	      for (kj=0; kj<4; ++kj)
		{
		  if (Vc[ixc[par+kr][ixd[par+1][ki]]].dist(sf_pts[par3][kj]) < tol &&
		      kj != ixs)
		    {
		      nmb++;
		      ixm = ki;
		      ixs = kj;
		      break;
		    }
		}
	    }
	  if (nmb == 0)
	    {
	      sfs[par+kr]->reverseParameterDirection(false);
	    }
	  else if (nmb == 1)
	    {
	      size_t ka;
	      for (ka=0; ka<degen_pts.size(); ++ka)
		if (Vc[ixc[par+kr][ixd[par+1][ixm]]].dist(degen_pts[ka]) < tol)
		  //Vc[ixc[par+kr][ixd[par+1][1]]].dist(degen_pts[ka]) < tol)
		  break;
	      if (ka == degen_pts.size())
		  sfs[par+kr]->swapParameterDirection();
	    }

	  get_corners(sfs[par+kr], sf_pts[par+kr]);
	  for (ki=0; ki<4; ++ki)
	    Vc[ixc[par+kr][ki]] = sf_pts[par+kr][ki];
	}

#ifdef DEBUG
      std::ofstream of3("lines3.g2");
      for (int ka=0; ka<4; ++ka)
	{
	  of3 << "410 1 0 4 255 0 0 255" << std::endl;
	  of3 << "1" << std::endl;
	  of3 << Vc[ixc[0][ka]] << "  " << Vc[ixc[1][ka]] << std::endl;
	}
      for (int ka=0; ka<4; ++ka)
	{
	  of3 << "410 1 0 4 0 255 0 255" << std::endl;
	  of3 << "1" << std::endl;
	  of3 << Vc[ixc[2][ka]] << "  " << Vc[ixc[3][ka]] << std::endl;
	}
      for (int ka=0; ka<4; ++ka)
	{
	  of3 << "410 1 0 4 55 100 100 255" << std::endl;
	  of3 << "1" << std::endl;
	  of3 << Vc[ixc[4][ka]] << "  " << Vc[ixc[5][ka]] << std::endl;
	}
#endif

      // Adapt surface to fit with identified corners
      for (kr=0; kr<6; ++kr)
	{
	  // if (kr == par || kr == par+1)
	  //   continue;   // Already treated

	  double dist[8];
	  dist[0] = dist[1] = dist[2] = dist[3] = dist[4] = dist[5] = dist[6] = dist[7] = 0.0;
	  for (ki=0; ki<8; ++ki)
	    {
	      for (kj=0; kj<4; ++kj)
		{
		  double curr_dist = Vc[ixc[kr][kj]].dist(sf_pts[kr][xs[ki][kj]]);
		  dist[ki] += curr_dist;
		}
	    }
	  // Find minimum accumulated distance between corners
	  double mindist = dist[0];
	  int minind = 0;
	  for (ki=1; ki<8; ++ki)
	    {
	      if (dist[ki] < mindist)
		{
		  mindist = dist[ki];
		  minind = ki;
		}
	    }
	  if (minind == 1)
	    sfs[kr]->reverseParameterDirection(true);
	  else if (minind == 2)
	    sfs[kr]->reverseParameterDirection(false);
	  else if (minind == 3)
	    sfs[kr]->swapParameterDirection();
	  else if (minind == 4)
	    {
	      sfs[kr]->swapParameterDirection();
	      sfs[kr]->reverseParameterDirection(false);
	    }
	  else if (minind == 5)
	    {
	      sfs[kr]->reverseParameterDirection(false);
	      sfs[kr]->reverseParameterDirection(true);
	    }
	  else if (minind == 6)
	    {
	      sfs[kr]->reverseParameterDirection(false);
	      sfs[kr]->swapParameterDirection();
	    }
	  else if (minind == 7)
	    {
	      sfs[kr]->reverseParameterDirection(true);
	      sfs[kr]->swapParameterDirection();
	      sfs[kr]->reverseParameterDirection(true);
	    }
	  get_corners(sfs[kr], sf_pts[kr]);
	  int stop_break = 1;
	}

      // All faces are done with parameter reversing and swapping
      // Now we test if all eight corners coincide
      vector<Point> u_pts(8), v_pts(8), w_pts(8);
      push_corners(u_pts, sf_pts[0], 0, 1, 2, 3);
      push_corners(u_pts, sf_pts[1], 4, 5, 6, 7);
      push_corners(v_pts, sf_pts[2], 0, 1, 4, 5);
      push_corners(v_pts, sf_pts[3], 2, 3, 6, 7);
      push_corners(w_pts, sf_pts[4], 0, 2, 4, 6);
      push_corners(w_pts, sf_pts[5], 1, 3, 5, 7);

      int i;
      for (i = 0; i < 8; ++i)
	if (u_pts[i].dist(v_pts[i]) >= tol ||
	    u_pts[i].dist(w_pts[i]) >= tol ||
	    v_pts[i].dist(w_pts[i]) >= tol)
	  break;
      if (i == 8)
	break;   // Sorting of side surfaces performed
    }
  if (par == 6)
    {
      MESSAGE("Faces lack common corner");
#ifdef DEBUG
      std::cout << "Faces lack common corner" << std::endl;
#endif
    }
  else
    break;
    }
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
  unif_sfs.insert(unif_sfs.begin(), sfs.begin(), sfs.begin()+4);
  GeometryTools::unifySurfaceSplineSpaceOneDir(unif_sfs, ptol, false);
  for (ki=0; ki<4; ++ki)
    sfs[ki] = unif_sfs[ki];

#ifdef DEBUG
  std::ofstream ofa("sfa.g2");
  for (size_t ka=0; ka<sfs.size(); ++ka)
    {
      sfs[ka]->writeStandardHeader(ofa);
      sfs[ka]->write(ofa);
    }
#endif

  unif_sfs.resize(0);
  sfs[4]->swapParameterDirection();
  sfs[5]->swapParameterDirection();
  unif_sfs.insert(unif_sfs.begin(), sfs.begin(), sfs.begin()+2);
  unif_sfs.insert(unif_sfs.end(), sfs.begin()+4, sfs.end());
  GeometryTools::unifySurfaceSplineSpaceOneDir(unif_sfs, ptol, true);
  sfs[0] = unif_sfs[0];
  sfs[1] = unif_sfs[1];
  sfs[4] = unif_sfs[2];
  sfs[5] = unif_sfs[3];
  sfs[4]->swapParameterDirection();
  sfs[5]->swapParameterDirection();

#ifdef DEBUG
  std::ofstream ofb("sfb.g2");
  for (size_t ka=0; ka<sfs.size(); ++ka)
    {
      sfs[ka]->writeStandardHeader(ofb);
      sfs[ka]->write(ofb);
    }
#endif

  unif_sfs.resize(0);
  unif_sfs.insert(unif_sfs.begin(), sfs.begin()+2, sfs.end());
  GeometryTools::unifySurfaceSplineSpaceOneDir(unif_sfs, ptol, true);
  for (ki=0; ki<4; ++ki)
    sfs[2+ki] = unif_sfs[ki];

#ifdef DEBUG
  std::ofstream ofc("sfc.g2");
  for (size_t ka=0; ka<sfs.size(); ++ka)
    {
      sfs[ka]->writeStandardHeader(ofc);
      sfs[ka]->write(ofc);
    }
#endif

  // Test if edge curves are equal
  if(!edge_curves_equal(sfs[0], 0.0, false, sfs[2], 0.0, false, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[0], 1.0, false, sfs[3], 0.0, false, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[1], 0.0, false, sfs[2], 1.0, false, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[1], 1.0, false, sfs[3], 1.0, false, tol))
    THROW("Edge curves are different");

  if(!edge_curves_equal(sfs[0], 0.0, true, sfs[4], 0.0, false, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[0], 1.0, true, sfs[5], 0.0, false, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[1], 0.0, true, sfs[4], 1.0, false, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[1], 1.0, true, sfs[5], 1.0, false, tol))
    THROW("Edge curves are different");

  if(!edge_curves_equal(sfs[2], 0.0, true, sfs[4], 0.0, true, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[2], 1.0, true, sfs[5], 0.0, true, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[3], 0.0, true, sfs[4], 1.0, true, tol))
    THROW("Edge curves are different");
  if(!edge_curves_equal(sfs[3], 1.0, true, sfs[5], 1.0, true, tol))
    THROW("Edge curves are different");

  return createCoonsPatchDirectly(sfs[0].get(), sfs[1].get(), sfs[2].get(), sfs[3].get(), 
				  sfs[4].get(), sfs[5].get());
}
