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

#include "GoTools/geometry/Curvature.h"
#include <vector>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "sislP.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;



namespace Go
{


//===========================================================================
void Curvature::curvatureRadiusPoints(const SplineCurve& curve,
			   double curveRad,
			   vector<double>& pos)
//===========================================================================
{
  
  BsplineBasis basis = curve.basis();
  int order = basis.order();
  int new_order = 6 * order - 11;

  vector<double> knots_simple;
  vector<double> new_knots;
  basis.knotsSimple(knots_simple);

  for (size_t i = 0; i < knots_simple.size(); ++i)
    {
      int old_mult = basis.knotMultiplicity(knots_simple[i]);
      int new_mult = 5 * order - 9 + old_mult;
      if (old_mult >= order-1) new_mult -= 2 - order + old_mult;
      for (int j = 0; j < new_mult; ++j) new_knots.push_back(knots_simple[i]);
    }

  int num_coefs = (int)new_knots.size() - new_order;

  BsplineBasis new_basis(new_order, new_knots.begin(), new_knots.end());

  vector<double> coefs_par; // Parameter values for the coefs (Greville)
  vector<double> coefs;

  shared_ptr<SplineCurve> d_curve(curve.derivCurve(1));
  shared_ptr<SplineCurve> dd_curve(d_curve->derivCurve(1));

  for (int i = 0; i < num_coefs; )
    {
      int knot_pos = i + 1;
      while (knot_pos < int(new_knots.size()) && new_knots[knot_pos-1] == new_knots[knot_pos])
	++knot_pos;
      double step = (new_knots[knot_pos] - new_knots[knot_pos-1]) / double(knot_pos - i);
      double tpar = new_knots[knot_pos-1] + step/2.0;
      for (;i < knot_pos; ++i)
	{
	  coefs_par.push_back(tpar);
	  tpar += step;
	}
    }

  for (int i = 0; i < num_coefs; ++i)
    {
      double tpar = coefs_par[i];
      Point d_p, dd_p;
      d_curve->point(d_p, tpar);
      dd_curve->point(dd_p, tpar);
      double d_p_length2 = d_p.length2();

      coefs.push_back( (d_p % dd_p).length2() * curveRad * curveRad - d_p_length2*d_p_length2*d_p_length2);
    }

  vector<double> curvature_coefs;
  vector<double> dummy_tangents;
  vector<int> dummy_index;
  SplineInterpolator interpolator;
  interpolator.setBasis(new_basis);
  interpolator.interpolate(coefs_par, coefs, dummy_index,
			   dummy_tangents, curvature_coefs);

  shared_ptr<SplineCurve> curvature_curve(
    new SplineCurve(num_coefs, new_order,
		    interpolator.basis().begin(),
		    curvature_coefs.begin(), 1));

  SISLCurve *num_sisl = Curve2SISL(*(curvature_curve.get()), false);
  SISLObject *qo1 = 0;
  SISLObject *qo2 = 0;
  SISLPoint *qp = 0;
  double spoint[1];
  spoint[0] = 0.0;
  int kstat = 0;
  SISLIntdat *qintdat = 0;
  double aepsge = 1.0e-9;

  if (!(qo1 = newObject(SISLCURVE))) goto error101;
  qo1 -> c1 = num_sisl;
  qo1 -> o1 = qo1;

  if (!(qo2 = newObject(SISLPOINT))) goto error101;
  spoint[0] = 0.0;
  if(!(qp = newPoint(spoint,1,1))) goto error101;
  qo2 -> p1 = qp;

  sh1761(qo1,qo2,aepsge,&qintdat,&kstat);
  if (kstat < 0) goto error101;

  if (qintdat)
  {
      for (int i = 0; i < qintdat->ipoint; ++i)
	  pos.push_back(qintdat->vpoint[i]->epar[0]);
  }

 error101:
  if (qo1)     freeObject(qo1);
  if (qo2)     freeObject(qo2);
  if (qintdat) freeIntdat(qintdat);
}



//===========================================================================
bool Curvature::minimalCurvatureRadius(const SplineCurve& curve,
				       double& mincurv,
				       double& pos)
//===========================================================================
{
  // Find the spline function for the numerator in the derivate of the square of the curvature

  BsplineBasis basis = curve.basis();
  int order = basis.order();
  int new_order = 6 * order - 14;

  mincurv = MAXDOUBLE;
  pos = (basis.startparam() + basis.endparam()) / 2.0;
  if (new_order < 0)    // Straight line
    {
      return false;
    }

  vector<double> knots_simple;
  vector<double> new_knots;
  basis.knotsSimple(knots_simple);

  for (size_t i = 0; i < knots_simple.size(); ++i)
    {
      int old_mult = basis.knotMultiplicity(knots_simple[i]);
      int new_mult = 5 * order - 11 + old_mult;
      if (old_mult >= order-2) new_mult -= 3 - order + old_mult;
      for (int j = 0; j < new_mult; ++j) new_knots.push_back(knots_simple[i]);
    }

  int num_coefs = (int)new_knots.size() - new_order;

  BsplineBasis new_basis(new_order, new_knots.begin(), new_knots.end());

  vector<double> coefs_par; // Parameter values for the coefs (Greville)
  vector<double> coefs;

  shared_ptr<SplineCurve> d_curve(curve.derivCurve(1));
  shared_ptr<SplineCurve> dd_curve(d_curve->derivCurve(1));
  shared_ptr<SplineCurve> ddd_curve(dd_curve->derivCurve(1));

  for (int i = 0; i < num_coefs; )
    {
      int knot_pos = i + 1;
      while (knot_pos < int(new_knots.size()) && new_knots[knot_pos-1] == new_knots[knot_pos])
	++knot_pos;
      double step = (new_knots[knot_pos] - new_knots[knot_pos-1]) / double(knot_pos - i);
      double tpar = new_knots[knot_pos-1] + step/2.0;
      for (;i < knot_pos; ++i)
	{
	  coefs_par.push_back(tpar);
	  tpar += step;
	}
    }

  for (int i = 0; i < num_coefs; ++i)
    {
      double tpar = coefs_par[i];
      Point d_p, dd_p, ddd_p;
      d_curve->point(d_p, tpar);
      dd_curve->point(dd_p, tpar);
      ddd_curve->point(ddd_p, tpar);

      coefs.push_back( (d_p % dd_p) * (d_p % (ddd_p * (d_p * d_p) - dd_p * 3 * (d_p * dd_p))) );
    }

  vector<double> numerator_coefs;
  vector<double> dummy_tangents;
  vector<int> dummy_index;
  SplineInterpolator interpolator;
  interpolator.setBasis(new_basis);
  interpolator.interpolate(coefs_par, coefs, dummy_index,
			   dummy_tangents, numerator_coefs);

  shared_ptr<SplineCurve> numerator_curve(
    new SplineCurve(num_coefs, new_order,
		    interpolator.basis().begin(),
		    numerator_coefs.begin(), 1));

  bool mincurvFound = false;
  vector<double> extremalParametervalues;

  for (size_t kr=1; kr<knots_simple.size(); ++kr)
    {
      shared_ptr<SplineCurve> sub_numerator(numerator_curve->subCurve(knots_simple[kr-1],
								      knots_simple[kr]));
      SISLCurve *num_sisl = Curve2SISL(*(sub_numerator.get()), false);
      SISLObject *qo1 = 0;
      SISLObject *qo2 = 0;
      SISLPoint *qp = 0;
      double spoint[1];
      spoint[0] = 0.0;
      int kstat = 0;
      SISLIntdat *qintdat = 0;
      double aepsge = 1.0e-9;

      if (!(qo1 = newObject(SISLCURVE)))
	continue;
      qo1 -> c1 = num_sisl;
      qo1 -> o1 = qo1;

      if (!(qo2 = newObject(SISLPOINT)))
	{
	  if (qo1)     freeObject(qo1);
	  continue;
	}
      spoint[0] = 0.0;
      if(!(qp = newPoint(spoint,1,1)))
	{
	  if (qo1)     freeObject(qo1);
	  if (qo2)     freeObject(qo2);
	  continue;
	}
      qo2 -> p1 = qp;

      sh1761(qo1,qo2,aepsge,&qintdat,&kstat);
      if (kstat < 0) 
	{
	  if (qo1)     freeObject(qo1);
	  if (qo2)     freeObject(qo2);
	  if (qintdat) freeIntdat(qintdat);
	  continue;
	}

      if (qintdat)
	{
	  for (int i = 0; i < qintdat->ipoint; ++i)
	    extremalParametervalues.push_back(qintdat->vpoint[i]->epar[0]);
	}
      if (qo1)     freeObject(qo1);
      if (qo2)     freeObject(qo2);
      if (qintdat) freeIntdat(qintdat);
    }

  for (size_t i = 0; i < extremalParametervalues.size(); ++i)
    {
      Point d_p, dd_p;
      d_curve->point(d_p, extremalParametervalues[i]);
      dd_curve->point(dd_p, extremalParametervalues[i]);
      double curv_numerator = (d_p % dd_p).length();
      double curv_denominator = d_p.length();
      curv_denominator = curv_denominator * curv_denominator * curv_denominator;
      if (curv_numerator <= curv_denominator * 1.0e-9) continue;

      double curvatureRadius = curv_denominator/curv_numerator;    // Same as 1/curvature
      if (!mincurvFound || curvatureRadius < mincurv)
	{
	  mincurvFound = true;
	  mincurv = curvatureRadius;
	  pos = extremalParametervalues[i];

	}
    }
  return (extremalParametervalues.size() > 0);
}


}   // End namespace Go
