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

#include "GoTools/creators/ProjectCurveAndCrossTan.h"
#include "GoTools/creators/CoonsPatchGen.h"

using namespace std;
using namespace Go;

//===========================================================================
ProjectCurveAndCrossTan::ProjectCurveAndCrossTan(const SplineCurve& space_crv, 
						 const SplineCurve& crosstan_crv,
						 const SplineSurface& surf,
						 const Point* start_par_pt, 
						 const Point* end_par_pt,
						 double epsgeo,
						 const RectDomain* domain_of_interest)
    : space_crv_(space_crv), crosstan_crv_(crosstan_crv), surf_(surf),
      start_par_pt_(start_par_pt), end_par_pt_(end_par_pt), epsgeo_(epsgeo),
      domain_of_interest_(domain_of_interest)
//===========================================================================
{
    // Test input
    ALWAYS_ERROR_IF(space_crv_.dimension() != 3 ||
		    space_crv_.dimension() != surf_.dimension(),
		    "Dimension mismatch.");
}

//===========================================================================
ProjectCurveAndCrossTan::~ProjectCurveAndCrossTan()
//===========================================================================
{
}

//===========================================================================
int ProjectCurveAndCrossTan::nmbCvs()
//===========================================================================
{
  return 3;
}

//===========================================================================
vector<Point> ProjectCurveAndCrossTan::eval(double t)
//===========================================================================
{
  const int nmb_return_pts = 3;
  vector<Point> return_pts(nmb_return_pts);

  Point space_pt = space_crv_.ParamCurve::point(t);
  double clo_u, clo_v, clo_dist;
  Point clo_pt, crosstan, par;
  vector<double> seed = createSeed(t);
  surf_.closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist, epsgeo_,
		      domain_of_interest_, &seed[0]);
	
  if ((t == start()) && (start_par_pt_ != 0)) 
    par = *start_par_pt_;
  else if ((t == end()) && (end_par_pt_ != 0))
    par = *end_par_pt_;
  else
    par.setValue(clo_u, clo_v);

  Point surf_norm;
  surf_.normal(surf_norm, clo_u, clo_v);
  surf_norm.normalize();

  Point cross_pt = crosstan_crv_.ParamCurve::point(t);
  double length = cross_pt.length();
  // The cross_pt is projected into the tangent plane in the surface by requiring
  // that it is orthogonal to the surface normal.
  Point vec = cross_pt - (cross_pt*surf_norm)*surf_norm;
  vec.normalize();

  crosstan = vec*length;
  return_pts[0] = par;
  return_pts[1] = clo_pt;
  return_pts[2] = crosstan;

  return return_pts;
}


//===========================================================================
void ProjectCurveAndCrossTan::eval(double t, int n, 
				   std::vector<std::vector<Go::Point> >& ders)
//===========================================================================
{
  MESSAGE_IF(n > 1, "Only one derivative will be computed");

  const int nmb_return_pts = 3;
  if (int(ders.size()) != nmb_return_pts)
      ders.resize(nmb_return_pts);
  if (n == 0)
    {
      vector<Point> evals = eval(t);
      for (size_t ki = 0; ki < evals.size(); ++ki) {
	ders.resize(1);
	ders[ki].push_back(evals[ki]);
      }
    }
  else {
      ders.resize(3);
      ders[0].resize(2);
      ders[1].resize(2);
      ders[2].resize(2);
      Point* par = &ders[0][0];
      Point* pos = &ders[1][0];
      Point* crosstan = &ders[2][0];

      int dim = space_crv_.dimension();
      vector<Point> space_pt(2);
      space_crv_.point(space_pt, t, 1); // We compute position and derivative of space curve.
      double clo_u, clo_v, clo_dist;
      vector<double> seed = createSeed(t);
      surf_.closestPoint(space_pt[0], clo_u, clo_v, pos[0], clo_dist, epsgeo_,
			  domain_of_interest_, &seed[0]);

      if ((t == start()) && (start_par_pt_ != 0)) {
	  par[0] = *start_par_pt_;
      } else if ((t == end()) && (end_par_pt_ != 0)) {
	  par[0] = *end_par_pt_;
      } else {
	 par[0] = Point(clo_u, clo_v);
      }

      Point surf_norm;
      surf_.normal(surf_norm, clo_u, clo_v);
      surf_norm.normalize();

      vector<Point> cross_pt;
      cross_pt.resize(2);
      crosstan_crv_.point(cross_pt, t, 1);
      double length = cross_pt[0].length();
      Point vec;
      vec = cross_pt[0] - (cross_pt[0]*surf_norm)*surf_norm;
      vec.normalize();

      crosstan[0] = vec*length;

      // From the surface we create the derivatives going in the u- and v-direction.
      vector<Point> surf_pts(6);
      surf_.point(surf_pts, par[0][0], par[0][1], 2);
      // We next describe the dir of space_crv as linear combination of the partial derivs.
      // space_pt[1] = s*surf_pts[1] + t*surf_pts[2].
      double coef1, coef2;
      CoonsPatchGen::blendcoef(&surf_pts[1][0], &surf_pts[2][0],
				 &space_pt[1][0], dim, 1, &coef1, &coef2);

      par[1] = Point(coef1, coef2);
      pos[1] = coef1*surf_pts[1] + coef2*surf_pts[2];
      crosstan[1] = cross_pt[1];

  }
}

//===========================================================================
double ProjectCurveAndCrossTan::start()
//===========================================================================
{
  return space_crv_.startparam();
}


//===========================================================================
double ProjectCurveAndCrossTan::end()
//===========================================================================
{
  return space_crv_.endparam();
}

//===========================================================================
int ProjectCurveAndCrossTan::dim()
//===========================================================================
{
    return 2; // Dimension of the parameter plane, not that of the space curve.
}

//===========================================================================
bool ProjectCurveAndCrossTan::approximationOK(double par, 
					      const std::vector<Go::Point>& approxpos,
					      double tol1, double tol2)
//===========================================================================
{
    // Only first tolerance is used.

//   double tol3 = 0.000001*tol1;
  vector<Point> der;
  der = eval(par);

  double dist = der[1].dist(approxpos[1]);

  return (dist < tol1);
}

//===========================================================================
vector<double> ProjectCurveAndCrossTan::createSeed(double tpar)
//===========================================================================
{
    vector<double> seed;
    if (start_par_pt_ == 0 || end_par_pt_ == 0) {
	return seed;
    }
    double tmin = start();
    double tmax = end();
    Point start_pt = *start_par_pt_;
    Point end_pt = *end_par_pt_;
    // We use convex combination of end pts.
    Point conv_comb_pt = (start_pt*(tmax - tpar)/(tmax - tmin) + end_pt*(tpar - tmin)/(tmax - tmin));

    seed.insert(seed.end(), conv_comb_pt.begin(), conv_comb_pt.end());
    return seed;
}
