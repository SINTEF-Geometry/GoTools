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

#include "GoTools/creators/ProjectCurve.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/ElementarySurface.h"

using namespace Go;
using std::vector;
using std::setprecision;
using std::endl;
using std::pair;
using std::make_pair;

//===========================================================================
ProjectCurve::ProjectCurve(shared_ptr<Go::ParamCurve>& space_crv, 
			   shared_ptr<Go::ParamSurface>& surf,
			   shared_ptr<Go::Point>& start_par_pt, 
			   shared_ptr<Go::Point>& end_par_pt,
			   double epsgeo1,
// 			   double epsgeo2,
			   const RectDomain* domain_of_interest)
    : space_crv_(space_crv), surf_(surf),
      start_par_pt_(start_par_pt), end_par_pt_(end_par_pt),
      epsgeo1_(epsgeo1),
//       epsgeo2_(epsgeo2),
      domain_of_interest_(domain_of_interest)
//===========================================================================
{
    // Test input
    ALWAYS_ERROR_IF(space_crv_.get() == 0 || surf_.get() == 0,
		"Missing input data.");
    ALWAYS_ERROR_IF(space_crv_->dimension() != 3 ||
		space_crv_->dimension() != surf_->dimension(),
		"Dimension mismatch.");

    //// Check spline surfaces for seems
    //shared_ptr<SplineSurface> tmp_srf = 
    //  dynamic_pointer_cast<SplineSurface,ParamSurface>(surf_);
    //if (tmp_srf.get())
    //  SurfaceTools::surfaceClosed(*tmp_srf, closed_dir_u_, closed_dir_v_);
    //else
    //  closed_dir_u_ = closed_dir_v_ = false;  // No particular treatment of seems

    // Check surfaces for seems
    closed_dir_u_ = false;
    closed_dir_v_ = false;
    shared_ptr<SplineSurface> tmp_srf = 
      dynamic_pointer_cast<SplineSurface,ParamSurface>(surf_);
    shared_ptr<ElementarySurface> tmp_srf_el = 
      dynamic_pointer_cast<ElementarySurface,ParamSurface>(surf_);
    if (tmp_srf.get())
      SurfaceTools::surfaceClosed(*tmp_srf, closed_dir_u_, closed_dir_v_);
    else if (tmp_srf_el.get())
        tmp_srf_el->isClosed(closed_dir_u_, closed_dir_v_);

    // Parameter domain
    RectDomain domain = surf_->containingDomain();
    umin_ = domain.umin();
    umax_ = domain.umax();
    vmin_ = domain.vmin();
    vmax_ = domain.vmax();
}

//===========================================================================
ProjectCurve::~ProjectCurve()
//===========================================================================
{
}

//===========================================================================
Point ProjectCurve::eval(double t) const
//===========================================================================
{
    if ((t == start()) && (start_par_pt_.get() != 0)) {
	return *start_par_pt_;
    } else if ((t == end()) && (end_par_pt_.get() != 0)) {
	return *end_par_pt_;
    } else {
	Point space_pt = space_crv_->ParamCurve::point(t);
	double clo_u, clo_v, clo_dist;
	Point clo_pt;
	vector<double> seed = createSeed(t);
	const double clo_pt_eps = std::min(epsgeo1_, 1e-10); // No need to be sloppy when finding closest pt.
	surf_->closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist, clo_pt_eps,//epsgeo1_,
			    domain_of_interest_, seed.size() > 0 ? &seed[0] : 0);
	// We may need to snap to the boundary.
	if (closeToSurfaceBoundary(clo_u, clo_v)) {
	    snapIfBoundaryIsCloser(space_pt, clo_u, clo_v, clo_dist);
	}

	// If we ended up on the seem of a closed surface we need to
	// make sure we've chosen the right side.
	double epsgeo = 1e-06;
	if ((seed.size() == 0) && // If we're using a seed we should be safe.
	    ((closed_dir_u_ && ((clo_u - umin_ < epsgeo) || (umax_ - clo_u < epsgeo))) ||
	     (closed_dir_v_ && ((clo_v - vmin_ < epsgeo) || (vmax_ - clo_v < epsgeo)))))
	  {
// 	      placeBorderPoint(t, clo_u, clo_v);
	      shared_ptr<Point> proj_pt =
		CreatorsUtils::projectCurvePoint(surf_.get(),
						 closed_dir_u_, closed_dir_v_,
						 space_crv_.get(), t);
	      if (proj_pt.get() == NULL) {
		  THROW("More than one param pts matched!");
	      }
	      clo_u = (*proj_pt)[0];
	      clo_v = (*proj_pt)[1];
	  }

	return Point(clo_u, clo_v);
    }
}

//===========================================================================
Point ProjectCurve::eval(double t, Point seed_pt) const
//===========================================================================
{
    if ((t == start()) && (start_par_pt_.get() != 0)) {
	return *start_par_pt_;
    } else if ((t == end()) && (end_par_pt_.get() != 0)) {
	return *end_par_pt_;
    } else {
	Point space_pt = space_crv_->ParamCurve::point(t);
	double clo_u, clo_v, clo_dist;
	Point clo_pt;
	vector<double> seed(seed_pt.begin(), seed_pt.end());
	const double clo_pt_eps = std::min(epsgeo1_, 1e-10); // No need to be sloppy when finding closest pt.
	surf_->closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist, clo_pt_eps,//epsgeo1_,
			    domain_of_interest_, 
			    seed.size() > 0 ? &seed[0] : 0);
	// We may need to snap to the boundary.
	if (closeToSurfaceBoundary(clo_u, clo_v)) {
	    snapIfBoundaryIsCloser(space_pt, clo_u, clo_v, clo_dist);
	}

	// If we ended up on the seem of a closed surface we need to
	// make sure we've chosen the right side.
//	double epsgeo = 1e-06;
	const Point sf_epspar = SurfaceTools::getParEpsilon(*surf_, epsgeo1_);
	if ((closed_dir_u_ && ((clo_u - umin_ < sf_epspar[0]) ||
			       (umax_ - clo_u < sf_epspar[0]))) ||
	    (closed_dir_v_ && ((clo_v - vmin_ < sf_epspar[1]) ||
			       (vmax_ - clo_v < sf_epspar[1]))))
	  {
// 	    placeBorderPoint(t, clo_u, clo_v);
	      shared_ptr<Point> proj_pt =
		CreatorsUtils::projectCurvePoint(surf_.get(),
						 closed_dir_u_, closed_dir_v_,
						 space_crv_.get(), t, epsgeo1_);
	      if (proj_pt.get() == NULL) {
		  THROW("More than one param pts matched!");
	      }
	      clo_u = (*proj_pt)[0];
	      clo_v = (*proj_pt)[1];
	  }

	return Point(clo_u, clo_v);
    }
}

//===========================================================================
void ProjectCurve::eval(double t, int n, Go::Point der[]) const
//===========================================================================
{
    MESSAGE_IF(n > 1, "Only one derivative will be computed");

    int dim = space_crv_->dimension();
    vector<Point> space_pt(2);
    space_crv_->point(space_pt, t, 1); // We compute position and
				       // derivative of space curve.
    if (n == 0)
	der[0] = eval(t);
    else if ((t == start()) && (start_par_pt_.get() != 0)) {
	der[0] = *start_par_pt_;
    } else if ((t == end()) && (end_par_pt_.get() != 0)) {
	der[0] = *end_par_pt_;
    } else {
	double clo_u, clo_v, clo_dist;
	Point clo_pt;
	vector<double> seed = createSeed(t);
	double* seed_ptr = (seed.size() == 0) ? NULL : &seed[0];
	const double clo_pt_eps = std::min(epsgeo1_, 1e-10); // No need to be sloppy when finding closest pt.
	surf_->closestPoint(space_pt[0], clo_u, clo_v, clo_pt, clo_dist, clo_pt_eps,//epsgeo1_,
			    domain_of_interest_, seed_ptr);
	// We may need to snap to the boundary.
        if (clo_dist > epsgeo1_)
        {   // We may experience large deviation if the curve is far away from the surface, with the
            // projection defining the 3d curve.
            //MESSAGE("clo_dist = " << clo_dist << ", epsgeo1_ = " << epsgeo1_);
        }
        // If a seed was used we do not replace the found value.
	if ((!seed_ptr) && (closeToSurfaceBoundary(clo_u, clo_v))) {
	    snapIfBoundaryIsCloser(space_pt[0], clo_u, clo_v, clo_dist);
	}

	// If we ended up on the seem of a closed surface we need to
	// make sure we've chosen the right side.
	// We use the direction of the tangent, the sf should lie to the left.
	double epspar = 1e-06;
	if ((seed.size() == 0) &&
	    ((closed_dir_u_ && ((clo_u - umin_ < epspar) ||
			       (umax_ - clo_u < epspar))) ||
	    (closed_dir_v_ && ((clo_v - vmin_ < epspar) ||
			       (vmax_ - clo_v < epspar))))) {
// 	  placeBorderPoint(t, clo_u, clo_v);
	    shared_ptr<Point> proj_pt =
	      CreatorsUtils::projectCurvePoint(surf_.get(),
					       closed_dir_u_, closed_dir_v_,
					       space_crv_.get(), t);
	    if (proj_pt.get() == NULL) {
		THROW("More than one param pts matched!");
	    }
	    clo_u = (*proj_pt)[0];
	    clo_v = (*proj_pt)[1];
	}

//       } else {
	der[0] = Point(clo_u, clo_v);
    }

    // From the surface we create the derivatives going in the u- and v-direction.
    vector<Point> surf_pts(3);
    surf_->point(surf_pts, der[0][0], der[0][1], 1);
    // We next describe the dir of space_crv as linear combination of the partial derivs.
    // space_pt[1] = s*surf_pts[1] + t*surf_pts[2].
    double coef1, coef2;
    CoonsPatchGen::blendcoef(&surf_pts[1][0], &surf_pts[2][0],
			     &space_pt[1][0], dim, 1, &coef1, &coef2);

    // If the surface is degenerate at the point we set the coef along the deg edge to 0.0.
    // It should be better to fetch it by stepping slightly away from the deg point.
    bool deg, b, r, top, l;
    const double deg_tol = 1.0e-06;
    deg = surf_->isDegenerate(b, r, top, l, deg_tol);
    if (deg)
    {
	const Point sf_epspar = SurfaceTools::getParEpsilon(*surf_, epsgeo1_);
        const double knot_tol = 1.0e-08;
        bool deg_pt = false;
        const double ang1 = space_pt[1].angle(surf_pts[1]);
        const double ang2 = space_pt[1].angle(surf_pts[2]);
        const double frac1 = space_pt[1].length()/surf_pts[1].length();
        const double frac2 = space_pt[1].length()/surf_pts[2].length();
        const RectDomain& rect_dom = surf_->containingDomain();
        // If the edge is degenerate and we are close to the edge we set the contribution in the
        // degenerate direction to 0.0. We also adjust the coef along the non-degenerate direction.
        if (b && (fabs(der[0][1] - rect_dom.vmin()) < sf_epspar[1]))
        {
            coef1 = 0.0;
            coef2 = (ang2 < 0.5*M_PI) ? frac2 : -frac2;
            deg_pt = true;
        }
        if (top && (fabs(der[0][1] - rect_dom.vmax()) < sf_epspar[1]))
        {
            coef1 = 0.0;
            coef2 = (ang2 < 0.5*M_PI) ? frac2 : -frac2;
            deg_pt = true;
        }
        if (l && (fabs(der[0][0] - rect_dom.umin()) < sf_epspar[0]))
        {
            coef2 = 0.0;
            coef1 = (ang1 < 0.5*M_PI) ? frac1 : -frac1;
            
            deg_pt = true;
        }
        if (r && (fabs(der[0][0] - rect_dom.umax()) < sf_epspar[0]))
        {
            coef2 = 0.0;
            coef1 = (ang1 < 0.5*M_PI) ? frac1 : -frac1;
            deg_pt = true;
        }

        // When stepping away from the degenerate point we need to add a test to ensure that we are far
        // enough away from the degeneracy.
        if (0)//deg_pt)
        {
            // We step slightly away from the deg pt.
            double upar = der[0][0];
            double vpar = der[0][1];
            if (b)
                vpar += sf_epspar[1];
            if (top)
                vpar -= sf_epspar[1];
            if (l)
                upar += sf_epspar[0];
            if (r)
                upar -= sf_epspar[0];
            
            surf_->point(surf_pts, upar, vpar, 1);
            // We next describe the dir of space_crv as linear combination of the partial derivs.
            // space_pt[1] = s*surf_pts[1] + t*surf_pts[2].
            CoonsPatchGen::blendcoef(&surf_pts[1][0], &surf_pts[2][0],
                                     &space_pt[1][0], dim, 1, &coef1, &coef2);
        }
    }

    der[1] = Point(coef1, coef2);
}

//===========================================================================
double ProjectCurve::start() const
//===========================================================================
{
  return space_crv_->startparam();
}

//===========================================================================
double ProjectCurve::end() const
//===========================================================================
{
  return space_crv_->endparam();
}

//===========================================================================
int ProjectCurve::dim() const
//===========================================================================
{
    return 2; // Dimension of the parameter plane, not that of the space curve.
}

//===========================================================================
bool ProjectCurve::approximationOK(double par, Go::Point approxpos,
				   double tol1, double tol2) const
//===========================================================================
{
    // Only first tolerance is used.
  Point pos;
  vector<double> seed = createSeed(par); // Uses end params only, not very robust.
  if (seed.size() > 0)
    pos = eval(par);
  else
      pos = eval(par, approxpos); // We use approxpos a seed.

  // We lift the points to the surface, and then compute their differences.
  Point space_approxpos = surf_->ParamSurface::point(approxpos[0], approxpos[1]); // The lifted param_cv pt.
  Point space_pos = surf_->ParamSurface::point(pos[0], pos[1]); // The proj pt (from space_cv).
//   Point spacecv_pos = space_crv_->ParamCurve::point(par);
  // The distance is wrt exact projection, evaluated in space.
  double dist1 = space_pos.dist(space_approxpos);
  bool appr_ok = (dist1 < epsgeo1_);
//   if (epsgeo2_ > 0.0) {
//       double dist2 = spacecv_pos.dist(space_approxpos);
//       if (dist2 > epsgeo2_)
// 	  appr_ok = false;
//   }

  // For some models with non-even parametrization and lots of
  // curvature, the seed may give us a closest point which is locally
  // closest but not globally. We use approxpos as seed.
  if (!appr_ok && seed.size() > 0)
  {
      pos = eval(par, approxpos);
      space_pos = surf_->ParamSurface::point(pos[0], pos[1]);
      dist1 = space_pos.dist(space_approxpos);
      appr_ok = (dist1 < epsgeo1_);
  }

  // @@sbr Currently no input tolerance is used, only the epsgeo_.
  return appr_ok;
}

//===========================================================================
vector<double> ProjectCurve::createSeed(double tpar) const
//===========================================================================
{
    vector<double> seed;
    double tmin = start();
    double tmax = end();
    if (start_par_pt_.get() == 0 && end_par_pt_.get() == 0)
    {
	return seed;
    }
    // else if (start_par_pt_.get() != 0 && end_par_pt_.get() != 0)
    // {   // For highly curved cases this seed can be far off.
    //     MESSAGE("This is just for debugging! To be removed!");
    // 	// We use convex combination of end pts.
    // 	Point start_pt = *start_par_pt_;
    // 	Point end_pt = *end_par_pt_;
    // 	Point conv_comb_pt = (start_pt*(tmax - tpar)/(tmax - tmin) +
    // 			      end_pt*(tpar - tmin)/(tmax - tmin));
    // 	seed.insert(seed.end(), conv_comb_pt.begin(), conv_comb_pt.end());
    // 	return seed;
    // }
    else
    {
	// If surf is not closed we're guessing it'll manage using the
	// control point net.
	if (!closed_dir_u_ && !closed_dir_v_)
	    return seed;
	vector<Point> cv_pt(2);
	space_crv_->point(cv_pt, tpar, 1);
	double clo_u, clo_v, clo_dist;
	double epsgeo = 1e-06;
	Point clo_bd_pt;
	surf_->closestBoundaryPoint(cv_pt[0], clo_u, clo_v,
				    clo_bd_pt, clo_dist, epsgeo);

	// We only bother to create seed when we're at the boarder.
	// VSK, 0209. Check if the closest point lies at a degenerate boundary
	// of a degenerate surface. In that case the parameter value along
	// the boundary is arbitrary and must be set more appropriate
	// Only applicable for endpoints of the curve
	if (tpar - tmin < epsgeo || tmax - tpar < epsgeo)
	  {
	    bool deg, b, r, t, l;
	    deg = surf_->isDegenerate(b, r, t, l, epsgeo);
	    if (deg)
	      {
		if ((l && clo_u-umin_ < epsgeo) || (r && umax_-clo_u < epsgeo) ||
		    (b && clo_v-vmin_ < epsgeo) || (t && vmax_-clo_v < epsgeo))
		  {
		    // Make reduced domain and perform a new closest point
		    double fac = 0.05;
		    double t2 = (tpar - tmin < epsgeo) ? tmin + fac*(tmax-tmin) 
		      : tmax - fac*(tmax-tmin);
		    Point cv_pt2 = space_crv_->ParamCurve::point(t2);
		    Array<double, 2> c1(l ? umin_+fac*(umax_-umin_) : umin_,
					b ? vmin_+fac*(vmax_-vmin_) : vmin_);
		    Array<double, 2> c2(r ? umax_-fac*(umax_-umin_) : umax_,
					t ? vmax_-fac*(vmax_-vmin_) : vmax_);
		    RectDomain red_dom(c1, c2);
		    double clo_u2, clo_v2, clo_dist2;
		    Point clo_pt2;
		    surf_->closestBoundaryPoint(cv_pt2, clo_u2, clo_v2,
						clo_pt2, clo_dist2, epsgeo, &red_dom);

		    // Modify boundary point with respect to new results
		    if ((l && clo_u-umin_ < epsgeo) || (r && umax_-clo_u < epsgeo))
		      clo_v = clo_v2;
		    if ((b && clo_v-vmin_ < epsgeo) || (t && vmax_-clo_v < epsgeo))
		      clo_u = clo_u2;
		  }
	      }
	  }

	// We only use the seed if surface is cyclic in that direction.
        // We then assume we stay on the same side of the seem, extrapolating the
        // start/end par point.
	if ((clo_dist < epsgeo) &&
	    ((closed_dir_u_ &&
		 ((clo_u - umin_ < epsgeo) || (umax_ - clo_u < epsgeo))) ||
		(closed_dir_v_ &&
		 ((clo_v - vmin_ < epsgeo) || (vmax_ - clo_v < epsgeo)))))
	{
	    Point base_par_pt = (start_par_pt_.get() != 0) ?
		*start_par_pt_ : *end_par_pt_;
	    double base_t = (start_par_pt_.get() != 0) ? tmin : tmax;

	    // Using end tangent. Project space curve onto dir derivs in par
	    // pt, and then assume linearity. Perhaps check make sure it
	    // ends up inside parameter domain.
	    vector<Point> surf_pt(3);
	    //surf_->point(surf_pt, base_par_pt[0], base_par_pt[1], 1);
	    surf_->point(surf_pt, clo_u, clo_v, 1);
	    double coef1, coef2;
	    int dim = surf_->dimension();
	    CoonsPatchGen::blendcoef(&surf_pt[1][0], &surf_pt[2][0],
				     &cv_pt[1][0], dim, 1, &coef1, &coef2);

            // If the surface is degenerate at the point we set the coef along the deg edge to 0.0.
            // It should be better to fetch it by stepping slightly away from the deg point.
            bool deg, b, r, t, l;
            deg = surf_->isDegenerate(b, r, t, l, epsgeo);
            if (deg)
            {
                const double knot_tol = 1.0e-08;
                const RectDomain& rect_dom = surf_->containingDomain();
                if (b && (fabs(clo_v - rect_dom.vmin()) < knot_tol))
                    coef1 = 0.0;
                if (t && (fabs(clo_v - rect_dom.vmax()) < knot_tol))
                    coef1 = 0.0;
                if (l && (fabs(clo_u - rect_dom.umin()) < knot_tol))
                    coef2 = 0.0;
                if (r && (fabs(clo_u - rect_dom.umax()) < knot_tol))
                    coef2 = 0.0;
            }

	    Point dir_der = Point(coef1, coef2);
            dir_der.normalize();
	    Point ext_pt = base_par_pt + (tpar - base_t)*dir_der;
	    seed.insert(seed.end(), ext_pt.begin(), ext_pt.end());
	    return seed;
	}
	else
	{
	    return seed;
	}
    }
}

// //===========================================================================
// void ProjectCurve::surfaceClosed(const SplineSurface& sf,
// 				 bool& closed_dir_u, bool& closed_dir_v) const
// //===========================================================================
// {
//   // Assuming k-regular surface, summing dist (L1-norm) between
//   // corresponding row coefs.
//   double num_tol = 1e-12;
//   double closed_tol = 1e-06;
//   int ik1 = sf.order_u();
//   int ik2 = sf.order_v();
//   int in1 = sf.numCoefs_u();
//   int in2 = sf.numCoefs_v();
//   vector<double>::const_iterator et1 = sf.basis_u().begin();
//   vector<double>::const_iterator et2 = sf.basis_v().begin();
//   ASSERT((et1[ik1-1] - et1[0] < num_tol) &&
// 	 (et1[in1+ik1-1] - et1[in1] < num_tol) &&
// 	 (et2[ik2-1] - et2[0] < num_tol) &&
// 	 (et2[in2+ik2-1] - et2[in2] < num_tol));
//   // Current case is not rational ...
//   ASSERT(!sf.rational());

//   // We first check vmin vs vmax.
//   double sum_dist = 0.0;
//   int dim = sf.dimension();
//   vector<double>::const_iterator coefs = sf.coefs_begin();
//   for (int ki = 0; ki < in1*dim; ++ki)
//     sum_dist += fabs(coefs[in1*dim*(in2-1)+ki] - coefs[ki]);
//   closed_dir_v = (sum_dist < closed_tol);

//   // Then umin vs umax.
//   sum_dist = 0.0;
//   for (int ki = 0; ki < in2; ++ki)
//     for (int kj = 0; kj < dim; ++kj)
//       sum_dist += fabs(coefs[(ki*in1+in1-1)*dim+kj] - coefs[ki*in1*dim+kj]);
//   closed_dir_u = (sum_dist < closed_tol);

//   return;
// }

//===========================================================================
void ProjectCurve::placeBorderPoint(double t,
				    double& upar, double& vpar) const
//===========================================================================
{
  std::cout << "Under construction!" << std::endl;

  double epspar = 1e-06;
  // double epsgeo = 1e-04;
  double angtol = 1e-02;

  // We're assuming (as an easy startm, which covers most cases) that
  // this only occurs for en end parameter on the curve.
  double tmin = space_crv_->startparam();
  double tmax = space_crv_->endparam();
  bool at_start = (t - tmin < epspar);
  bool at_end = (tmax - t < epspar);
//   if (!at_start && !at_end) {
//       MESSAGE("Inner parameter value, not yet supported!");
//       return;
//   }

//   // If closed in both directions we should
//   if (closed_dir_u_ && closed_dir_v_)
//       MESSAGE("Closed in both directions, not fully supported!");

  // We're of course assuming that we're not to cross the
  // seem of a closed surface.
  // Evaluating the tangent of the space curve from both directions.
  vector<Point> cv_pt(2);
  bool eval_from_right = at_start;
  space_crv_->point(cv_pt, t, 1, eval_from_right);

  int dim = surf_->dimension();
  if (at_end)
      for (int ki = 0; ki < dim; ++ki)
	  cv_pt[1][ki] *= -1.0;

  // We then extract the surface derivs, in the various
  // locations. Possibly test 4 parameter points (most likely only 2).
  vector<Point> cand_par_pts;
  cand_par_pts.push_back(Point(upar, vpar));
  if (closed_dir_u_) {
      if (upar - umin_ < epspar)
	  cand_par_pts.push_back(Point(umax_, vpar));
      else if (umax_ - upar < epspar)
	  cand_par_pts.push_back(Point(umin_, vpar));
  }
  if (closed_dir_v_) {
      if (vpar - vmin_ < epspar)
	  cand_par_pts.push_back(Point(upar, vmax_));
      else if (vmax_ - vpar < epspar)
	  cand_par_pts.push_back(Point(upar, vmin_));
  }
  if (closed_dir_u_ && closed_dir_v_)
      cand_par_pts.push_back(Point(cand_par_pts[1][0], cand_par_pts[2][1]));

  // We then check the directional derivs in the params.
  double min_sum_ang = 3*M_PI;
  int min_sum_ind = -1;
  for (size_t ki = 0; ki < cand_par_pts.size(); ++ki) {
    vector<Point> sf_pt(3);
      surf_->point(sf_pt, cand_par_pts[ki][0], cand_par_pts[ki][1], 1);

      // We flip the directional derivs into the surface (if on the border).
      if (umax_ - cand_par_pts[ki][0] < epspar)
	  for (int kj = 0; kj < dim; ++kj)
	      sf_pt[1][kj] *= -1.0;
      if (vmax_ - cand_par_pts[ki][1] < epspar)
	  for (int kj = 0; kj < dim; ++kj)
	      sf_pt[2][kj] *= -1.0;

      // We then sum the angles.
      double angu = cv_pt[1].angle(sf_pt[1]);
      double angv = cv_pt[1].angle(sf_pt[2]);
      double sum_ang = angu + angv;
      if (fabs(sum_ang - min_sum_ang) < angtol)
	{
	  // This may happen along the seem of a closes surface.
	  // But if along a seem we should be able to decide which to choose
	  // depending on the orientation of the tangent. I.e. the sf
	  // should lie to the left of the tangent (assuming it is a trim
	  // curve, of course, which may not always be the case @@sbr ...).
	  // We should then abort, and the user must pick end pts.
	  MESSAGE("Two par pts yield the same reliability!");
	  if ((closed_dir_v_ && angu > 0.5*M_PI) || // Almost 0.0, not
						    // 0.5*M_PI
	      (closed_dir_u_ && angv < 0.5*M_PI))
	  { // We choose the upper cv.
	      min_sum_ind = (int)ki;
	      min_sum_ang = sum_ang;
	    }
	}
      if (sum_ang < min_sum_ang) {
	  min_sum_ind = (int)ki;
	  min_sum_ang = sum_ang;
      }
  }

  ASSERT(min_sum_ind != -1);

  upar = cand_par_pts[min_sum_ind][0];
  vpar = cand_par_pts[min_sum_ind][1];

  return;

//   // We then should step in both directions (unless we're at
//   // an end parameter) to see if we get to an internal
//   // parameter.

}


//===========================================================================
bool ProjectCurve::closeToSurfaceBoundary(double upar, double vpar,
					  double domain_fraction) const
//===========================================================================
{
    const RectDomain& rect_domain = surf_->containingDomain();
    double umin = rect_domain.umin();
    double umax = rect_domain.umax();
    double vmin = rect_domain.vmin();
    double vmax = rect_domain.vmax();
    double eps_u = (umax - umin)*domain_fraction;
    double eps_v = (vmax - vmin)*domain_fraction;
    if ((fabs(upar - umin) < eps_u) || (fabs(umax - upar) < eps_u)
	|| (fabs(vpar - vmin) < eps_v) || (fabs(vmax - vpar) < eps_v))
	return true;
    else
	return false;
}


//===========================================================================
void ProjectCurve::snapIfBoundaryIsCloser(Point space_pt,
					  double& upar, double& vpar, double& dist) const
//===========================================================================
{
    shared_ptr<ParamSurface> sf = surf_;
    if (sf->instanceType() == Class_BoundedSurface)
    { // Since we typically project from trim curves we are not interested in finding
      // closest point to the same space trim curve.
	sf = (dynamic_pointer_cast<BoundedSurface>(sf))->underlyingSurface();
    }
    double bd_clo_u, bd_clo_v, bd_clo_dist;
    Point bd_clo_pt;
    double local_seed[2];
    local_seed[0] = upar;
    local_seed[1] = vpar;
    double epsgeo = 1e-06;
    sf->closestBoundaryPoint(space_pt, bd_clo_u, bd_clo_v, 
			     bd_clo_pt, bd_clo_dist, epsgeo, NULL, local_seed);
    if (bd_clo_dist < dist)
    {
	dist = bd_clo_dist;
	upar = bd_clo_u;
	vpar = bd_clo_v;
    }
}

