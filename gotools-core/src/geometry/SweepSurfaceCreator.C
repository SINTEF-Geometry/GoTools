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

#include "GoTools/geometry/SweepSurfaceCreator.h"
#include <memory>
#include "GoTools/utils/MatrixXD.h"


using std::vector;
using std::setprecision;
using std::endl;
using std::pair;
using std::make_pair;


namespace Go
{



//===========================================================================
SplineSurface* SweepSurfaceCreator::linearSweptSurface(const SplineCurve &curv1,
						       const SplineCurve &curv2,
						       const Point &pt)
//===========================================================================
{

  int dim = curv1.dimension();

  ALWAYS_ERROR_IF(dim != curv2.dimension(),
		  "Can not sweep from two curves in spaces of different dimension");

  ALWAYS_ERROR_IF(dim != pt.dimension(),
		  "Can not sweep: Space of point has different dimension than space of curves");

  bool rational1 = curv1.rational();
  bool rational2 = curv2.rational();
  bool rational = rational1 || rational2;
  int kdim_v = dim + (rational2 ? 1 : 0);
  int kdim = dim + (rational ? 1 : 0);

  const BsplineBasis bas_u = curv1.basis();
  const BsplineBasis bas_v = curv2.basis();
  int numCoefs_u = bas_u.numCoefs();
  int numCoefs_v = bas_v.numCoefs();
  vector<double> coefs(numCoefs_u * numCoefs_v * kdim);
  int pos = 0;

  vector<double>::const_iterator coefs_u;
  vector<double>::const_iterator coefs_v;
  if (rational2)
    coefs_v = curv2.rcoefs_begin();
  else
    coefs_v = curv2.coefs_begin();

  for (int j = 0; j < numCoefs_v; ++j, coefs_v += kdim_v)
    {
      double weight_v = 1.0;
      if (rational2)
	weight_v = coefs_v[dim];

      if (rational1)
	coefs_u = curv1.rcoefs_begin();
      else
	coefs_u = curv1.coefs_begin();

      for (int i = 0; i < numCoefs_u; ++i)
	{
	  double weight_u = 1.0;
	  if (rational1)
	    weight_u = coefs_u[dim];
	  double weight_uv = weight_u * weight_v;

	  for (int k = 0; k < dim; ++k, ++coefs_u)
	    coefs[pos++] = (coefs_v[k]/weight_v + (*coefs_u)/weight_u - pt[k]) * weight_uv;
	  if (rational)
	    {
	      coefs[pos++] = weight_uv;
	      if (rational1)
		++coefs_u;
	    }
	}
    }

  return new SplineSurface(bas_u, bas_v, coefs.begin(), dim, rational);

}


//===========================================================================
SplineSurface* SweepSurfaceCreator::rotationalSweptSurface(const SplineCurve &curve,
							   double angle,
							   const Point &pt,
							   const Point &axis)
//===========================================================================
{
  ALWAYS_ERROR_IF(curve.dimension() != 3 || pt.dimension() != 3 || axis.dimension() != 3,
		  "Curve, point and axis must lie in three-dimensional space");

  const double PIHALF = 1.57079632679489661923;
  const double THREEPIHALF = 4.71238898038468985769;
  const double TWOPI = 6.28318530717958647692;

  bool rational_curve = curve.rational();
  int curve_kdim = rational_curve ? 4 : 3;

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

  // Create surface
  int numCoefs_curve = curve.numCoefs();
  int numCoefs_seg = circle_segment->numCoefs();
  vector<double> surface_rcoefs(numCoefs_curve * numCoefs_seg * 4);
  const Point axis_norm = axis / axis.length();

  vector<double> temp_pt;
  vector<double>::const_iterator curve_pt_it;
  if (rational_curve)
    curve_pt_it = curve.rcoefs_begin();
  else
    curve_pt_it = curve.coefs_begin();
  vector<double>::iterator surface_pt_it = surface_rcoefs.begin();
  for (int j = 0; j<numCoefs_curve; ++j, curve_pt_it += curve_kdim)
    {
      double curve_weight = 1.0;
      if (rational_curve)
	curve_weight = curve_pt_it[3];
      Point curve_pt(curve_pt_it[0]/curve_weight, curve_pt_it[1]/curve_weight, curve_pt_it[2]/curve_weight);

      // Projection of curve_pt to the axis
      Point curve_proj = pt + ((curve_pt-pt)*axis_norm) * axis_norm;

      // Vector from curve_proj to curve_pt
      Point xaxis = curve_pt - curve_proj;
      double rad = xaxis.length();

      // The vector yaxis will be the unique vector of same length as xaxsis,
      // making xaxis, yaxis, axis_norm an ortonormal right hand basis
      Point yaxis = axis_norm % xaxis;

      // Create the matrix M from xaxis, yaxis, axis_norm, scaled by the height of curve_pt
      // above the axis
      MatrixXD<double, 3> matr;
      for (int i = 0; i < 3; ++i)
	{
	  matr(i,0) = xaxis[i];
	  matr(i,1) = yaxis[i];
	  matr(i,2) = rad * axis_norm[i];
	}

      vector<double>::const_iterator curve_seg_it = circle_segment->rcoefs_begin();
      for (int i = 0; i < numCoefs_seg; ++i, curve_seg_it += 4, surface_pt_it += 4)
	{
	  matr.mult(curve_seg_it,surface_pt_it);

	  // Rescale, as the curve segment coordinates where in a homogeneous system
	  for (int k = 0; k < 3; ++k) surface_pt_it[k] /= curve_seg_it[3];

	  // Translate
	  for (int k = 0; k < 3; ++k) surface_pt_it[k] += curve_proj[k];

	  // Rescale back , to make surface points homogeneous
	  surface_pt_it[3] = curve_seg_it[3] * curve_weight;
	  for (int k = 0; k < 3; ++k) surface_pt_it[k] *= surface_pt_it[3];

	}

    }

  return new SplineSurface(circle_segment->basis(), curve.basis(),
			   surface_rcoefs.begin(), 3, true);
}


}  // namespace Go
