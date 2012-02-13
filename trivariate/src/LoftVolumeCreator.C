//===========================================================================
//
// File : LoftVolumeCreator.C
//
// Created: Thu Dec 18 13:41:02 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: LoftVolumeCreator.C,v 1.1 2008-12-19 09:54:29 kfp Exp $
//
// Description:
//
//===========================================================================


#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/geometry/HermiteInterpolator.h"
#include "GoTools/geometry/GeometryTools.h"
#include <algorithm>
#include <cmath>
#include <fstream> // For debugging.

using namespace Go;
using std::vector;


namespace Go
{



//===========================================================================
SplineVolume*
LoftVolumeCreator::loftVolume(vector<shared_ptr<SplineSurface> >::iterator
			      first_surface,
			      int nmb_srfs)
//===========================================================================
{
  vector<shared_ptr<SplineSurface> > unified_surfaces = unifiedSurfacesCopy(first_surface, nmb_srfs);
  vector<double> params;
  // We're giving the volume a parameter domain of length 1.0 in w-direction
  makeLoftParams(unified_surfaces.begin(), nmb_srfs, 1.0, params);

  return loftVolumeFromUnifiedSurfaces(unified_surfaces.begin(), params.begin(), nmb_srfs);
}



  
//===========================================================================
SplineVolume*
LoftVolumeCreator::loftVolume(vector<shared_ptr<SplineSurface> >::iterator
			      first_surface,
			      vector<double>::iterator first_param,
			      int nmb_srfs)
//===========================================================================
{
  vector<shared_ptr<SplineSurface> > unified_surfaces = unifiedSurfacesCopy(first_surface, nmb_srfs);
  return loftVolumeFromUnifiedSurfaces(unified_surfaces.begin(), first_param, nmb_srfs);
}




//===========================================================================
vector<shared_ptr<SplineSurface> >
LoftVolumeCreator::unifiedSurfacesCopy(vector<shared_ptr<SplineSurface> >::iterator
				       first_surface,
				       int nmb_srfs)
//===========================================================================
{
  bool rational = false;
  for (int i = 0; i < nmb_srfs; ++i)
    if (first_surface[i]->rational())
      rational = true;

  // Create the copies
  vector<shared_ptr<SplineSurface> > unified_surfaces;
  for (int i = 0; i < nmb_srfs; ++i)
    {
      shared_ptr<SplineSurface> surface_copy(first_surface[i]->clone());
      if (rational)
	surface_copy->representAsRational();
      unified_surfaces.push_back(surface_copy);
    }

  // Reparametrize, to have the B-spline spaces living on the same intervals
  double
    avg_start_u = 0.0,
    avg_start_v = 0.0,
    avg_end_u = 0.0,
    avg_end_v = 0.0;

  for (int i = 0; i < nmb_srfs; ++i)
    {
      avg_start_u += unified_surfaces[i]->startparam_u();
      avg_end_u += unified_surfaces[i]->endparam_u();
      avg_start_v += unified_surfaces[i]->startparam_v();
      avg_end_v += unified_surfaces[i]->endparam_v();
    }

  avg_start_u /= double(nmb_srfs);
  avg_end_u /= double(nmb_srfs);
  avg_start_v /= double(nmb_srfs);
  avg_end_v /= double(nmb_srfs);

  for (int i = 0; i < nmb_srfs; ++i)
    unified_surfaces[i]->setParameterDomain(avg_start_u, avg_end_u, avg_start_v, avg_end_v);


  // Put the surfaces into common basis.
  double tolerance = 1e-05;
  GeometryTools::unifySurfaceSplineSpace(unified_surfaces, tolerance);

  return unified_surfaces;
}




//===========================================================================
void LoftVolumeCreator::makeLoftParams(vector<shared_ptr<SplineSurface> >::const_iterator
				       first_surface,
				       int nmb_srfs, double param_length,
				       vector<double>& params)
//===========================================================================
{
  params.clear();
  bool rational = first_surface[0]->rational();

  // Compute parameterization.
  // For each adjacent pair of surfaces compute the minimum
  // distance between the coefficients
  int dim = first_surface[0]->dimension();
  int num_coefs = first_surface[0]->numCoefs_u() * first_surface[0]->numCoefs_v();

  params.push_back(0.0);  // Parameter value of first surface.
  for (int i = 1; i < nmb_srfs; ++i)
    {
      vector<double>::const_iterator coefs1, coefs2;
      if (rational)
	{
	  coefs1 = first_surface[i-1]->rcoefs_begin();
	  coefs2 = first_surface[i]->rcoefs_begin();
	}
      else
	{
	  coefs1 = first_surface[i-1]->coefs_begin();
	  coefs2 = first_surface[i]->coefs_begin();
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
  double scale = param_length / params[nmb_srfs-1];
  for (int i = 1; i < nmb_srfs; ++i)
      params[i] *= scale;

}




//===========================================================================
SplineVolume*
LoftVolumeCreator::loftVolumeFromUnifiedSurfaces(vector<shared_ptr<SplineSurface> >::iterator
						 first_surface,
						 vector<double>::iterator first_param,
						 int nmb_srfs)
//===========================================================================
{

  SplineVolume* vol = loftNonrationalVolume(first_surface, first_param, nmb_srfs);
  if (first_surface[0]->rational())
    {
      int n = vol->numCoefs(0) * vol->numCoefs(1) * vol->numCoefs(2);
      int kdim = vol->dimension() + 1;
      bool all_positive = true;
      vector<double>::const_iterator it = vol->rcoefs_begin();
      it += (kdim - 1);
      for (int i = 0; i < n; ++i)
	if (it[kdim * i] <= 0.0)
	  {
	    all_positive = false;
	    break;
	  }
      if (!all_positive)
	{
	  delete vol;
	  vol = loftRationalVolume(first_surface, first_param, nmb_srfs);
	}
    }

  return vol;
}




//===========================================================================
SplineVolume*
LoftVolumeCreator::loftNonrationalVolume(vector<shared_ptr<SplineSurface> >::iterator
					 first_surface,
					 vector<double>::iterator first_param,
					 int nmb_srfs)
//===========================================================================
{

  // By looking at each surface as a point of dimension dim*u_coefs.size()*v_coefs.size(),
  // the interpolation method gets rather easy.

  bool rational = first_surface[0]->rational();

  vector<double> coefs, coefs_uv_loft, params(nmb_srfs);
  for (int i = 0; i < nmb_srfs; ++i)
    {
      if (rational)
	coefs.insert(coefs.end(),
		     first_surface[i]->rcoefs_begin(), first_surface[i]->rcoefs_end());
      else
	coefs.insert(coefs.end(),
		     first_surface[i]->coefs_begin(), first_surface[i]->coefs_end());
      params[i] = first_param[i];
    }

  vector<double> tangent_points;
  vector<int> cross_index;

  int order = std::min(4, nmb_srfs);
  SplineInterpolator u_interpolator;
  u_interpolator.makeBasis(params, cross_index, order);
  u_interpolator.interpolate(params, coefs, cross_index,
			     tangent_points, coefs_uv_loft);

  // We are lofting in third parameter (w) direction.
  SplineVolume* vol = new SplineVolume(first_surface[0]->basis_u(), first_surface[0]->basis_v(),
				       u_interpolator.basis(),
				       coefs_uv_loft.begin(), first_surface[0]->dimension(), rational);
  return vol;
}




//===========================================================================
SplineVolume* LoftVolumeCreator::loftRationalVolume(vector<shared_ptr<SplineSurface> >::iterator first_surface,
						    vector<double>::iterator first_param,
						    int nmb_srfs)
//===========================================================================
{

  double tolerance = 1e-05;

  ASSERT(nmb_srfs>2);

  int dim = (*first_surface)->dimension();

  // Find tangents for later Hermite interpolation.
  // First create the vector ratsurfaces of the homogeneous coordinates of the original surfaces, containing both the original coordinates and the weights themselves
  // Use this to get a common basis for all surfaces

  vector<shared_ptr<SplineSurface> > ratsurfaces;
  for (int i = 0; i < nmb_srfs; ++i)
    {
      shared_ptr<SplineSurface> current_surface = first_surface[i];
      shared_ptr<SplineSurface> rsurf(new SplineSurface(current_surface->basis_u(), current_surface->basis_v(),
							current_surface->rcoefs_begin(),
							dim + 1, false));
      ratsurfaces.push_back(rsurf);
    }

  GeometryTools::unifySurfaceSplineSpace(ratsurfaces, tolerance);

  // Now create the high-dimensional curve nat_interp_curve, the cubic natural interpolation of the coordinates in ratsurfaces divided out with weights.
  // This is used to get the tangents for cubic Hermite interpolation. We do not need the weights, as their derivatives will be zero, to make weight
  // function positive

  int ncoefs_u = ratsurfaces[0]->numCoefs_u();
  int ncoefs_v = ratsurfaces[0]->numCoefs_v();

  vector<double> coefs_homog;
  for (vector<shared_ptr<SplineSurface> >::const_iterator it_surf = ratsurfaces.begin(); it_surf != ratsurfaces.end(); ++it_surf)
    {
      vector<double>::const_iterator it_coefs = (*it_surf)->coefs_begin();
      for (int j = 0; j < ncoefs_u * ncoefs_v; ++j)
	{
	  double denom = it_coefs[dim];
	  for (int k = 0; k < dim; ++k, ++it_coefs)
	    coefs_homog.push_back((*it_coefs)/denom);
	  ++it_coefs;   // Skip weight
	}
    }

  vector<double> coefs_nat_interp;
  SplineInterpolator nat_interp;
  nat_interp.setNaturalConditions();
  nat_interp.interpolate(nmb_srfs, ncoefs_u * ncoefs_v * dim, &first_param[0], &coefs_homog[0], coefs_nat_interp);
  BsplineBasis nat_interp_basis = nat_interp.basis();

  SplineCurve nat_interp_curve(nat_interp_basis.numCoefs(), nat_interp_basis.order(),
			       nat_interp_basis.begin(), coefs_nat_interp.begin(),
			       ncoefs_u * ncoefs_v * dim);


  // Find the tangents of nat_interp_curve

  vector<double> tangent_nat_interp;
  vector<Point> pt_and_deriv(2);
  for (int i = 0; i < nmb_srfs; ++i)
    {
      nat_interp_curve.point(pt_and_deriv, first_param[i], 1);
      for (int j = 0; j < ncoefs_u * ncoefs_v * dim ; ++j)
	tangent_nat_interp.push_back(pt_and_deriv[1][j]);
    }

  // Now we create the cubic Hermite interpolant, also for weights
  // First we push points and tangents to interpolate into one flat array

  vector<double> pos_tan_interp;
  int pos_t = 0;
  for (vector<shared_ptr<SplineSurface> >::const_iterator it_surf = ratsurfaces.begin(); it_surf != ratsurfaces.end(); ++it_surf)
    {
      // Push points
      vector<double>::const_iterator it_coefs = (*it_surf)->coefs_begin();
      for (int j = 0; j < ncoefs_u * ncoefs_v * (dim + 1); ++j, ++it_coefs)
	pos_tan_interp.push_back(*it_coefs);

      // Push tangents
      it_coefs = (*it_surf)->coefs_begin();
      for (int j = 0; j < ncoefs_u * ncoefs_v; ++j)
	{
	  double weight = it_coefs[dim];

	  // The homogeneous tangents are given by the product rule, but one term falls out as the weight function derivatives are 0.
	  for (int k = 0; k < dim; ++k, ++it_coefs, ++pos_t)
	    pos_tan_interp.push_back(tangent_nat_interp[pos_t] * weight);
	  // Put 0 for the weight derivative
	  pos_tan_interp.push_back(0.0);
	  ++it_coefs; // Skip weight
	}
    }

  // Then we create the interpolant

  vector<double> param_twice;
  for (int i = 0; i < nmb_srfs; ++i)
    {
      param_twice.push_back(first_param[i]);
      param_twice.push_back(first_param[i]);
    }
  HermiteInterpolator herm_interp;
  vector<double> coefs_hermite_interp;
  herm_interp.interpolate(nmb_srfs * 2, ncoefs_u * ncoefs_v * (dim + 1),
			  &param_twice[0],
			  &pos_tan_interp[0],
			  coefs_hermite_interp);

  // Now we have the coefficients of our volume as coefs_hermite_interp
  // Next get knot vector

  vector<double> knots_w;
  knots_w.push_back(first_param[0]);
  knots_w.push_back(first_param[0]);
  for (int i = 0; i < nmb_srfs; ++i)
    {
      knots_w.push_back(first_param[i]);
      knots_w.push_back(first_param[i]);
    }
  knots_w.push_back(first_param[nmb_srfs-1]);
  knots_w.push_back(first_param[nmb_srfs-1]);

  BsplineBasis basis_w((int)knots_w.size() - 4, 4, knots_w.begin());

  // Finally, create the lofted volume

  return new SplineVolume(ratsurfaces[0]->basis_u(), ratsurfaces[0]->basis_v(), basis_w,
			  coefs_hermite_interp.begin(), dim, true);
}



} // end namespace Go.
