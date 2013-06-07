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

#include "GoTools/trivariate/SweepVolumeCreator.h"
#include <memory>
#include "GoTools/utils/MatrixXD.h"


using std::vector;


namespace Go
{



//===========================================================================
SplineVolume* SweepVolumeCreator::linearSweptVolume(const SplineSurface &surface,
						    const SplineCurve &curve,
						    const Point &pt)
//===========================================================================
{

  int dim = surface.dimension();

  ALWAYS_ERROR_IF(dim != curve.dimension(),
		  "Can not sweep: Space of surface has different dimension than space of curve");

  ALWAYS_ERROR_IF(dim != pt.dimension(),
		  "Can not sweep: Space of point has different dimension than space of surface and curve");

  bool rational_surf = surface.rational();
  bool rational_curve = curve.rational();
  bool rational = rational_surf || rational_curve;
  int kdim_curve = dim + (rational_curve ? 1 : 0);
  int kdim = dim + (rational ? 1 : 0);

  const BsplineBasis bas_u = surface.basis_u();
  const BsplineBasis bas_v = surface.basis_v();
  const BsplineBasis bas_w = curve.basis();
  int numCoefs_uv = bas_u.numCoefs() * bas_v.numCoefs();
  int numCoefs_w = bas_w.numCoefs();
  vector<double> coefs(numCoefs_uv * numCoefs_w * kdim);
  int pos = 0;

  vector<double>::const_iterator coefs_uv;
  vector<double>::const_iterator coefs_w;
  if (rational_curve)
    coefs_w = curve.rcoefs_begin();
  else
    coefs_w = curve.coefs_begin();

  for (int j = 0; j < numCoefs_w; ++j, coefs_w += kdim_curve)
    {
      double weight_w = 1.0;
      if (rational_curve)
	weight_w = coefs_w[dim];

      if (rational_surf)
	coefs_uv = surface.rcoefs_begin();
      else
	coefs_uv = surface.coefs_begin();

      for (int i = 0; i < numCoefs_uv; ++i)
	{
	  double weight_uv = 1.0;
	  if (rational_surf)
	    weight_uv = coefs_uv[dim];
	  double weight_uvw = weight_uv * weight_w;

	  for (int k = 0; k < dim; ++k, ++coefs_uv)
	    coefs[pos++] = (coefs_w[k]/weight_w + (*coefs_uv)/weight_uv - pt[k]) * weight_uvw;
	  if (rational)
	    {
	      coefs[pos++] = weight_uvw;
	      if (rational_surf)
		++coefs_uv;
	    }
	}
    }

  return new SplineVolume(bas_u, bas_v, bas_w, coefs.begin(), dim, rational);

}


//===========================================================================
SplineVolume* SweepVolumeCreator::rotationalSweptVolume(const SplineSurface &surface,
							double angle,
							const Point &pt,
							const Point &axis)
//===========================================================================
{
  ALWAYS_ERROR_IF(surface.dimension() != 3 || pt.dimension() != 3 || axis.dimension() != 3,
		  "Surface, point and axis must lie in three-dimensional space");

  const double PIHALF = 1.57079632679489661923;
  const double THREEPIHALF = 4.71238898038468985769;
  const double TWOPI = 6.28318530717958647692;

  bool rational_surface = surface.rational();
  int surface_kdim = rational_surface ? 4 : 3;

  // First create a full unit circle in the xy-plane of order 3 and 9 control points

  // Create knot vector for unit circle
  vector<double> full_circle_knots(12);
  full_circle_knots[0] = full_circle_knots[1] = full_circle_knots[2] = 0.0;
  full_circle_knots[3] = full_circle_knots[4] = PIHALF;
  full_circle_knots[5] = full_circle_knots[6] = M_PI;
  full_circle_knots[7] = full_circle_knots[8] = THREEPIHALF;
  full_circle_knots[9] = full_circle_knots[10] = full_circle_knots[11] = TWOPI;

  // Create 9 control points for unit curcle (angles 0, pi/4, pi/2, ... , 2*pi)
  vector<double> full_circle_rcoefs(36);
  // Insert points for angle 0
  full_circle_rcoefs[0] = full_circle_rcoefs[3] = 1.0;
  full_circle_rcoefs[1] = full_circle_rcoefs[2] = 0.0;
  // Insert points for angle pi/4
  full_circle_rcoefs[4] = full_circle_rcoefs[5] = full_circle_rcoefs[7] = sqrt(0.5);
  full_circle_rcoefs[6] = 0.0;
  // For the remainig points, rotate previous points a quarter-circle
  for (int i=8; i < 36; i+=4)
    {
      full_circle_rcoefs[i] = -full_circle_rcoefs[(i-8) | 1];
      full_circle_rcoefs[i|1] = full_circle_rcoefs[i-8];
      full_circle_rcoefs[i|2] = 0.0;
      full_circle_rcoefs[i|3] = full_circle_rcoefs[(i-8) | 3];;
    }

  // Create the unit circle
  SplineCurve full_circle(9, 3, full_circle_knots.begin(), full_circle_rcoefs.begin(), 3, true);

  // Next create the circle segment to rotate along
  double abs_angle = angle >= 0.0 ? angle : -angle;
  if (abs_angle > TWOPI)
    abs_angle = TWOPI;
  double quadrant = floor (abs_angle / PIHALF);

  // Get difference between abs_angle and its closest diagonal line
  double diag_angle = abs_angle - (quadrant + 0.5) * PIHALF;
  // Get end parameter for circle segment
  double end_angle = PIHALF * (quadrant
			       + (1.0 
				  + (sqrt ((double) 2.0) + 1.0)
				    * tan (diag_angle / 2.0))
			         / 2.0);

  // Create circle segment
  shared_ptr<SplineCurve> circle_segment(full_circle.subCurve(0.0, end_angle));

  // Reverse (i.e. mirror across x-axis) if input angle was negative
  if (angle < 0.0)
    {
      vector<double>::iterator rcoefs = circle_segment->rcoefs_begin();
      for (int i=0; i<circle_segment->numCoefs(); ++i)
	rcoefs[i<<2 | 1] = -rcoefs[i<<2 | 1];
    }

  // Create volume
  int numCoefs_surface = surface.numCoefs_u() * surface.numCoefs_v();
  int numCoefs_seg = circle_segment->numCoefs();
  vector<double> volume_rcoefs(numCoefs_surface * numCoefs_seg * 4);
  const Point axis_norm = axis / axis.length();

  vector<double> temp_pt;
  vector<double>::const_iterator surface_pt_it;
  if (rational_surface)
    surface_pt_it = surface.rcoefs_begin();
  else
    surface_pt_it = surface.coefs_begin();
  vector<double>::iterator volume_pt_it = volume_rcoefs.begin();
  for (int j = 0; j<numCoefs_surface; ++j, surface_pt_it += surface_kdim)
    {
      double surface_weight = 1.0;
      if (rational_surface)
	surface_weight = surface_pt_it[3];
      Point surface_pt(surface_pt_it[0]/surface_weight, surface_pt_it[1]/surface_weight, surface_pt_it[2]/surface_weight);

      // Projection of surface_pt to the axis
      Point surface_proj = pt + ((surface_pt-pt)*axis_norm) * axis_norm;

      // Vector from surface_proj to surface_pt
      Point xaxis = surface_pt - surface_proj;
      double rad = xaxis.length();

      // The vector yaxis will be the unique vector of same length as xaxsis,
      // making xaxis, yaxis, axis_norm an ortonormal right hand basis
      Point yaxis = axis_norm % xaxis;

      // Create the matrix M from xaxis, yaxis, axis_norm, scaled by the height of surface_pt
      // above the axis
      MatrixXD<double, 3> matr;
      for (int i = 0; i < 3; ++i)
	{
	  matr(i,0) = xaxis[i];
	  matr(i,1) = yaxis[i];
	  matr(i,2) = rad * axis_norm[i];
	}

      vector<double>::const_iterator curve_seg_it = circle_segment->rcoefs_begin();
      for (int i = 0; i < numCoefs_seg; ++i, curve_seg_it += 4, volume_pt_it += 4)
	{
	  matr.mult(curve_seg_it,volume_pt_it);

	  // Rescale, as the curve segment coordinates where in a homogeneous system
	  for (int k = 0; k < 3; ++k) volume_pt_it[k] /= curve_seg_it[3];

	  // Translate
	  for (int k = 0; k < 3; ++k) volume_pt_it[k] += surface_proj[k];

	  // Rescale back , to make volume points homogeneous
	  volume_pt_it[3] = curve_seg_it[3] * surface_weight;
	  for (int k = 0; k < 3; ++k) volume_pt_it[k] *= volume_pt_it[3];
	}

    }

  return new SplineVolume(circle_segment->basis(), surface.basis_u(), surface.basis_v(),
			  volume_rcoefs.begin(), 3, true);
}



}  // namespace Go
