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

#include "GoTools/qualitymodule/QualityUtils.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/ClassType.h"


using std::vector;


namespace Go
{

namespace qualityUtils
{



//===========================================================================
bool isSliverFace(shared_ptr<ParamSurface> surf,
		  double thickness,
		  double factor)
//===========================================================================
{
  ClassType ct = surf -> instanceType();
  if (ct == Class_SplineSurface)
    {
      SplineSurface *ss;
      ss = dynamic_cast<SplineSurface*>(surf.get());
      return isSliverFace(*ss, thickness, factor);
    }
  if (ct == Class_BoundedSurface)
    {
      double tol2d = 1.0e-4;
      BoundedSurface *bs;
      bs = dynamic_cast<BoundedSurface*>(surf.get());
      if (bs->isIsoTrimmed(tol2d))
	{
	  // Get surrounding domain
	  RectDomain domain = bs->containingDomain();
    
	  // Get smallest surrounding surface
	  shared_ptr<ParamSurface> base_sf = bs->underlyingSurface();
	  while (base_sf->instanceType() == Class_BoundedSurface)
	    base_sf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(base_sf)->underlyingSurface();
	  RectDomain dom2 = base_sf->containingDomain();  // To avoid problems due to numerics
	  double umin = std::max(domain.umin(), dom2.umin());
	  double umax = std::min(domain.umax(), dom2.umax());
	  double vmin = std::max(domain.vmin(), dom2.vmin());
	  double vmax = std::min(domain.vmax(), dom2.vmax());
    
	  vector<shared_ptr<ParamSurface> > sfs = base_sf->subSurfaces(umin, vmin, 
								       umax, vmax);

	  SplineSurface *ss;
	  ss = dynamic_cast<SplineSurface*>(sfs[0].get());
	  if (ss)
	    return isSliverFace(*ss, thickness, factor);
	  else
	    return isSliverFace2(*bs, thickness, factor);
	}
	else
	  return isSliverFace2(*bs, thickness, factor);
    }
  else THROW("Unexpected surface type");
}


//===========================================================================
bool isSliverFace(const SplineSurface& sf,
		  double thickness,
		  double factor)
//===========================================================================
{
  double start_u = sf.startparam_u();
  double end_u = sf.endparam_u();
  double start_v = sf.startparam_v();
  double end_v = sf.endparam_v();

  double mid_u = (start_u + end_u) / 2.0;
  double mid_v = (start_v + end_v) / 2.0;

  // Start by measuring the lengths at the mid curves, use this to
  // rule out most of the surfaces before we start more complex searches

  double len_u, len_v;
  GeometryTools::estimateIsoCurveLength(sf, true, mid_v, len_u);
  GeometryTools::estimateIsoCurveLength(sf, false, mid_u, len_v);

  double
    len_thin, start_thin, end_thin,
    len_thick, start_thick, end_thick;
  BsplineBasis
    basis_thin,
    basis_thick;

  bool thin_u = len_u < len_v;

  if (thin_u)
    {
      // First parameter direction is the thin side candidate
      len_thin = len_u;
      start_thin = start_u;
      end_thin = end_u;
      basis_thin = sf.basis_u();

      len_thick = len_v;
      start_thick = start_v;
      end_thick = end_v;
      basis_thick = sf.basis_v();
    }
  else
    {
      // Second parameter direction is the thin side candidate
      len_thin = len_v;
      start_thin = start_v;
      end_thin = end_v;
      basis_thin = sf.basis_v();

      len_thick = len_u;
      start_thick = start_u;
      end_thick = end_u;
      basis_thick = sf.basis_u();
    }

  if (len_thin > thickness || len_thin*factor > len_thick) return false;

  // First check the edge curves along the thick side, in case of triangles etc.
  double local_len;

  GeometryTools::estimateIsoCurveLength(sf, !thin_u, start_thin, local_len);
  if (local_len < len_thick)
    {
      if (len_thin * factor > local_len) return false;
      len_thick = local_len;
    }

  GeometryTools::estimateIsoCurveLength(sf, !thin_u, end_thin, local_len);
  if (local_len < len_thick)
    {
      if (len_thin * factor > local_len) return false;
      len_thick = local_len;
    }

  // We now need to search more in detail, and try to get an estimate of the
  // minimum length across the thick side and a maximum length across the thin size

  // First we maximize the length across the thin side
  vector<double> knots_thick;
  basis_thick.knotsSimple(knots_thick);
  int deg_thick = basis_thick.order() - 1;

  for (size_t i = 0; i < knots_thick.size()-1; ++i)
    {
      double step = (knots_thick[i+1] - knots_thick[i]) / double(deg_thick);
      for (int j=0; j <= deg_thick; ++j)
	{
	  if (j==0 && i!=0) continue;   // Avoid testing twice along inner knot values
	  GeometryTools::estimateIsoCurveLength(sf, thin_u, knots_thick[i] + step*double(j), local_len);
	      if (local_len > len_thin)
		{
		  if (local_len > thickness || local_len*factor > len_thick) return false;
		  len_thin = local_len;
		}
	    }
	}

  // Then we minimize the length across the thick side
  vector<double> knots_thin;
  basis_thin.knotsSimple(knots_thin);
  int deg_thin = basis_thin.order() - 1;

  for (size_t i = 0; i < knots_thin.size()-1; ++i)
    {
      double step = (knots_thin[i+1] - knots_thin[i]) / double(deg_thin);
      for (int j=0; j < deg_thin; ++j)
	{
	  if (j==0 && i==0) continue;   // Edge is already tested
	  GeometryTools::estimateIsoCurveLength(sf, !thin_u, knots_thin[i] + step*double(j), local_len);
	  if (local_len < len_thick)
	    {
	      if (local_len > thickness || len_thin*factor > local_len) return false;
	      len_thick = local_len;
	    }
	}
    }

  return true;
}


//===========================================================================
bool isSliverFace2(const BoundedSurface& sf,
		  double thickness,
		  double factor)
//===========================================================================
{
  // Get the domain surrounding the surface
  RectDomain dom = sf.containingDomain();

  // Estimate surface size
  double len_u, len_v;
  double area[4];
  area[0] = dom.umin();
  area[1] = dom.umax();
  area[2] = dom.vmin();
  area[3] = dom.vmax();
  GeometryTools::estimateSurfaceSize(sf, len_u, len_v, area);

  // Number of samples
  int nmb_u, nmb_v;
  int max_sample = 30;
  int min_sample = 3;
  double fac = 10.0*thickness;
  nmb_u = (int)(len_u*fac);
  nmb_u = std::max(min_sample, std::min(max_sample, nmb_u));
  nmb_v = (int)(len_v*fac);
  nmb_v = std::max(min_sample, std::min(max_sample, nmb_v));

  // Check mid parameter
  double min_len_u, max_len_u, med_len_u=0.0;
  double min_len_v, max_len_v, med_len_v=0.0;
  vector<shared_ptr<ParamCurve> > constParCurves;

  constParCurves = sf.constParamCurves(0.5*(dom.vmax()+dom.vmin()), false);
  double len;
  size_t l;
  for (l = 0, len=0.0; l < constParCurves.size(); ++l)
    len += constParCurves[l]->estimatedCurveLength();
  min_len_u = max_len_u = len;

  constParCurves.clear();
  constParCurves = sf.constParamCurves(0.5*(dom.umax()+dom.umin()), true);
  for (l = 0, len=0.0; l < constParCurves.size(); ++l)
    len += constParCurves[l]->estimatedCurveLength();
  min_len_v = max_len_v = len;
  if (min_len_u > thickness && min_len_v > thickness)
    return false;  // The surface is not thin enough to be
  // a sliver face
  
  // Compute curve length estimates in the u directions
  // Avoid the boundaries of the surrounding domain to reduce the
  // risk of creating very difficult intersection problems
  int ki;
  double del = (dom.vmax()-dom.vmin())/(double)nmb_v;
  double tpar = dom.vmin() + 0.5*del;
  for (ki=0; ki<nmb_v; ++ki, tpar+=del)
    {
      // Get the curve segments in the constant parameter direction
      constParCurves.clear();
      constParCurves = sf.constParamCurves(tpar, false);
      

      // Get the estimated curve length by summing up for each curve 
      // in the constant parameter direction
      double est_len = 0.0;
      for (l = 0; l < constParCurves.size(); ++l)
	est_len += constParCurves[l]->estimatedCurveLength();

//       std::cout << "Par: " << tpar << ", nmb of cvs: " << constParCurves.size();
//       std::cout << ", len: " << est_len << std::endl;

      min_len_u = std::min(min_len_u, est_len);
      max_len_u = std::max(max_len_u, est_len);
      med_len_u += est_len;
    }
  med_len_u /= (double)nmb_v;

  // The v-direction
  // Check sliver face charactersistics as we go along
  del = (dom.umax()-dom.umin())/(double)nmb_u;
   tpar = dom.umin() + 0.5*del;
  for (ki=0; ki<nmb_u; ++ki, tpar+=del)
    {
      // Get the curve segments in the constant parameter direction
      constParCurves.clear();
      constParCurves = sf.constParamCurves(tpar, true);
      

      // Get the estimated curve length by summing up for each curve 
      // in the constant parameter direction
      double est_len = 0.0;
      for (l = 0; l < constParCurves.size(); ++l)
	est_len += constParCurves[l]->estimatedCurveLength();

//       std::cout << "Par: " << tpar << ", nmb of cvs: " << constParCurves.size();
//       std::cout << ", len: " << est_len << std::endl;
      min_len_v = std::min(min_len_v, est_len);
      max_len_v = std::max(max_len_v, est_len);
      med_len_v += est_len;

      if (min_len_u > thickness && min_len_v > thickness)
	return false;  // The surface is not thin enough to be
      // a sliver face
    }
  med_len_v /= (double)nmb_u;

  return (med_len_u/max_len_v >= factor ||
	  med_len_v/max_len_u >= factor);
}

//===========================================================================
bool isSliverFace(const BoundedSurface& sf,
		  double thickness,
		  double factor)
//===========================================================================
{
  bool thin_u;      // Tells whether the first constant parameter direction is the thin side
  bool found_u = false;     // Tells whether any curve in first constant parameter direction has been found
  bool found_v = false;     // Tells whether any curve in second constant parameter direction has been found

  // - As long as only curves in one constant paramter direction is found,
  //   len_thick will be the shortest curve length and len_thin the longest.
  // - When both constant parameter directions are found,
  //   len_thick will be the shortest length along the thick side,
  //   and len_thin the longest length along the thin side
  double len_thick, len_thin;

  // Iterate all boundary curves
  CurveLoop cl = sf.outerBoundaryLoop();
  for (int i = 0; i < cl.size(); ++i)
    {
      SplineCurve *sc;
      sc = cl[i] -> geometryCurve();
      BsplineBasis basis =  sc -> basis();
      vector<double> knots;
      basis.knotsSimple(knots);

      int deg = basis.order() - 1;

      // Iterate all bezier segments of the curve
      for (size_t j = 0; j < knots.size()-1; ++j)
	{
	  double step = (knots[j+1] - knots[j]) / double(deg);
	  double curve_param = knots[j];
	  double curve_param2;

	  // Iterate points along the Bezier segment
	  for (int k=0; k <= deg; curve_param += step, ++k)
	    {
	      if (j!=0 && k==0) continue;   // No need to do interior knot values twice
	      vector<Point> pts(2);
	      //	      pts.resize(2);
	      curve_param2 = curve_param;
	      if (j==0 && k == 0)
		curve_param2 += 0.1*step;
	      if (j==knots.size()-2 && k==deg)
		curve_param2 -= 0.1*step;

	      sc -> point(pts, curve_param2, 1);
	      Point c_pt = pts[0];
	      Point c_deriv = pts[1];

	      // Test direction
	      double sf_param_u, sf_param_v;
	      Point clo_pt;
	      double clo_dist;
	      sf.closestPoint(c_pt, sf_param_u, sf_param_v, clo_pt, clo_dist, 0.01);

	      vector<Point> sf_pts(3);
	      sf.point(sf_pts, sf_param_u, sf_param_v, 1, true, true);

	      // u_dir indicates whether the first constant paramater direction on the surface is most normal to the boundary curve
	      // This is the same as if A_u (the angle between the curve tangent and the surface tangent in first parameter direction)
	      // is greater than A_v (the angle between the curve tangent and the surface tangent in the second paramter direction).
	      // This again is the same as |cos(A_u)| < |cos(A_v)|, which means
	      //
	      //    |sf_pts[1]*c_deriv| / (|sf_pts[1]| * |c_deriv|)   <   |sf_pts[2]*c_deriv| / (|sf_pts[2]| * |c_deriv|)
	      //
	      // Squaring the above, and multiplying up the denominators gives the following test
	      bool u_dir =
		(sf_pts[1]*c_deriv) * (sf_pts[1]*c_deriv) * (sf_pts[2]*sf_pts[2]) <
		(sf_pts[2]*c_deriv) * (sf_pts[2]*c_deriv) * (sf_pts[1]*sf_pts[1]);

	      // Get the curve segments in the constant parameter direction
	      vector<shared_ptr<ParamCurve> > constParCurves;
	      if (u_dir)
		constParCurves = sf.constParamCurves(sf_param_v, false);
	      else
		constParCurves = sf.constParamCurves(sf_param_u, true);
	      if (constParCurves.size() == 0) continue;

	      // Get the estimated curve length by summing up for each curve in the constant parameter direction
	      double est_len = 0.0;
	      for (size_t l = 0; l < constParCurves.size(); ++l)
		est_len += constParCurves[l]->estimatedCurveLength();

	      if (!found_u && !found_v)   // Very first time a length is measured
		{
		  found_u = u_dir;
		  found_v = !found_u;
		  len_thick = len_thin = est_len;
		}

	      else if (found_u != found_v && u_dir == found_u)  // Length has only been measured in one direction, this is another one in that same direction
		{
		  if (len_thick > est_len) len_thick = est_len;
		  if (len_thin < est_len) len_thin = est_len;
		}

	      else if (found_u != found_v && u_dir != found_u)  // Length has only been measured in one direction, this is the first time in the other direction
		{
		  found_u = found_v = true;
		  if (len_thick > est_len)
		    {
		      len_thin = est_len;
		      thin_u = u_dir;
		    }
		  else if (len_thin < est_len)
		    {
		      len_thick = est_len;
		      thin_u = !u_dir;
		    }
		  else
		    return false;

		  if (len_thin > thickness || len_thin*factor > len_thick) return false;
		}

	      else  // Length has already been measured in bouth directions
		{
		  if (thin_u == u_dir)   // Another contribution to the thin side
		    {
		      if (len_thin < est_len)
			{
			  len_thin = est_len;
			  if (len_thin > thickness || len_thin*factor > len_thick) return false;
			}
		    }
		  else   // Another contribution to the thick side
		    {
		      if (len_thick > est_len)
			{
			  len_thick = est_len;
			  if (len_thin*factor > len_thick) return false;
			}
		    }
		}
	    }
	}
    }

  // The boundary curve is traversed

  if (!found_u || !found_v)
    {
      // Error handling if not curves have been found in both directions ?
    }

  return true;

}

//===========================================================================
bool hasIndistinctKnots(shared_ptr<ParamSurface> surf, double tol,
			vector<shared_ptr<ParamCurve> >& trim_cv_knots)
//===========================================================================
{
    shared_ptr<SplineSurface> spline_surf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);
    if (spline_surf.get())
    {
	// Spline surface. Check the knot vector of the surface
	vector<double> knot_u;
	vector<double> knot_v;
	bool occurance_u = spline_surf->basis_u().indistinctKnots(tol, knot_u);
	bool occurance_v = spline_surf->basis_v().indistinctKnots(tol, knot_v);
	if (occurance_u || occurance_v)
	    return true;
	else
	    return false;
    }
    else
    {
	shared_ptr<BoundedSurface> bd_surf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
	if (bd_surf.get())
	{
	    // Check underlying surface
	    shared_ptr<ParamSurface> surf2 = bd_surf->underlyingSurface();
	    bool occurance = hasIndistinctKnots(surf2, tol, trim_cv_knots);

	    // Check trimming curves
	    int nmb_loops = bd_surf->numberOfLoops();
	    for (int ki=0; ki<nmb_loops; ki++)
	    {
		shared_ptr<CurveLoop> curr_loop = bd_surf->loop(ki);
		vector<shared_ptr<ParamCurve> >::const_iterator crv = curr_loop->begin();
		vector<shared_ptr<ParamCurve> >::const_iterator last = curr_loop->end();
		for (; crv!= last; crv++)
		{
		    shared_ptr<SplineCurve> spline_cv = 
			dynamic_pointer_cast<SplineCurve, ParamCurve>(*crv);
		    if (spline_cv.get())
		    {
			vector<double> cv_knot;
			bool cv_occurance = spline_cv->basis().indistinctKnots(tol, cv_knot);
			if (cv_occurance)
			    trim_cv_knots.push_back(*crv);
		    }
		    else
		    {
			shared_ptr<CurveOnSurface> sf_cv = 
			dynamic_pointer_cast<CurveOnSurface, ParamCurve>(*crv);
			if (sf_cv.get())
			{
			    shared_ptr<ParamCurve> pcrv = sf_cv->parameterCurve();
			    shared_ptr<ParamCurve> space_crv = sf_cv->spaceCurve();

			    shared_ptr<SplineCurve> p_spline = 
				dynamic_pointer_cast<SplineCurve, ParamCurve>(pcrv);
			    shared_ptr<SplineCurve> g_spline = 
				dynamic_pointer_cast<SplineCurve, ParamCurve>(space_crv);
			    vector<double> cv_knot;
			    bool cv_occurance = false;
			    if (p_spline.get())
			    {
				cv_occurance = p_spline->basis().indistinctKnots(tol, cv_knot);
				if (cv_occurance)
				    trim_cv_knots.push_back(*crv);
			    }
			    if (g_spline.get() && !cv_occurance)
			    {
				cv_occurance = g_spline->basis().indistinctKnots(tol, cv_knot);
				if (cv_occurance)
				    trim_cv_knots.push_back(*crv);
			    }
			}
		    }
		}

	    }
	    return occurance;
	}
	else
	    return false;  // No knot vector to worry about
    }

}

//===========================================================================
  double estimateArea(shared_ptr<ParamSurface> surf)
//===========================================================================
{
  // Check for spline surface
    shared_ptr<SplineSurface> spline_surf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);
    if (spline_surf.get())
    {
	// Spline surface
      double len_u, len_v;
      GeometryTools::estimateSurfaceSize(*spline_surf, len_u, len_v);
      return len_u*len_v;
    }

    // Not a spline surface. Check for bounded surface
    
    shared_ptr<BoundedSurface> bd_surf = 
      dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
    if (bd_surf.get())
      {
	// For each loop, make a triangulation and compute the area
	// of this triangulation. The outer loop makes a positive
	// contribution to the area estimate and inner loops make
	// negative contributions
	int nmb_loop = bd_surf->numberOfLoops();
	double area = estimateLoopArea(bd_surf->loop(0));
	for (int ki=1; ki<nmb_loop; ++ki)
	  {
	    double inner_area = estimateLoopArea(bd_surf->loop(ki));
	    area -= inner_area;
	  }
	return area;
      }

    // Not a bounded surface. Return zero to trigger a proper
    // area computation
    return 0.0;

}

//===========================================================================
  double estimateLoopArea(shared_ptr<CurveLoop> loop)
//===========================================================================
{
  // For each curve in the loop, evaluate a number of points to make a 
  // polygon approximating the loop
  int nmb_eval = 3;
  vector<Point> pol;
  int nmb_crvs = loop->size();
  int ki, kj, kk;
  for (ki=0; ki<nmb_crvs; ++ki)
    {
      shared_ptr<ParamCurve> crv = (*loop)[ki];
      double start = crv->startparam();
      double end = crv->endparam();
      double tint = (end-start)/(double)(nmb_eval-1);
      double par;
      for (par=start, kj=1; kj<nmb_eval; ++kj, par+=tint)
	{
	  Point pnt = crv->point(par);
	  pol.push_back(pnt);
	}
    }

  // Triangulate the domain surrounded by this polyon and compute the area
  // of this triangulation
  int nmb = (int)pol.size();
  ki = 0;
  kj = nmb-1;
  double area = 0.0;
  while (kj-ki > 1)
    {
      // Choose the third point in the triangle
      Point vec1 = pol[kj]-pol[ki];
      Point vec2 = pol[ki]-pol[ki+1];
      Point vec3 = pol[kj-1]-pol[kj];

      double ang1 = vec1.angle(vec2);
      double ang2 = vec1.angle(vec3);

      double curr_area;
      if (ang1 > ang2)
	{
	  kk = ki+1;
	  double l1 = vec1.length();
	  double l2 = vec2.length();
	  curr_area = 0.5*sin(ang1)*l1*l2;
	}
      else
	{
	  kk = kj-1;
	  double l1 = vec1.length();
	  double l2 = vec3.length();
	  curr_area = 0.5*sin(ang2)*l1*l2;
	}
      area += curr_area;
      if (kk == ki+1)
	ki = kk;
      else 
	kj = kk;
    }
  return area;
}


}   // namespace qualityUtils

}   // namespace Go
