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

#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/creators/SmoothCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/SurfaceCreators.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/geometry/SurfaceTools.h"

#include <memory>
#include <fstream>


using namespace Go;
using std::vector;
using std::max;
using std::min;
using std::swap;
using std::cout;
using std::endl;

//===========================================================================
SplineCurve* CreatorsUtils::
getParametricCurve(const vector<shared_ptr<const CurveOnSurface> >& cv)
//===========================================================================
{
    if (cv.size() == 0)
        return NULL;

    shared_ptr<const ParamCurve> pcurve(cv[0]->parameterCurve());
    if (pcurve.get() == 0)
        return NULL;

    SplineCurve* return_cv = dynamic_cast<SplineCurve*>(pcurve->clone());
    double dummy_dist;

    if (return_cv == 0)
        return return_cv;
    for (size_t ki = 1; ki < cv.size(); ++ki) {
        const SplineCurve* curr_segment
            = dynamic_cast<const SplineCurve*>(cv[ki]->parameterCurve().get());
        if (curr_segment == 0)
            return NULL;
        else {
            shared_ptr<ParamCurve> cp_cv(curr_segment->clone());
            return_cv->appendCurve(cp_cv.get(), 0, dummy_dist);
        }
    }

    return return_cv;
}


//===========================================================================
shared_ptr<SplineCurve> CreatorsUtils::
createCrossTangent(const CurveOnSurface& cv)
//===========================================================================
{
    shared_ptr<SplineCurve> dummy_cv1;
    SplineCurve* dummy_cv2 = NULL;
    return createCrossTangent(cv, dummy_cv1, dummy_cv2);
}


//===========================================================================
shared_ptr<Go::SplineCurve> CreatorsUtils::
createCrossTangent(const Go::CurveOnSurface& cv,
                   shared_ptr<Go::SplineCurve> basis_space_cv,
                   const Go::SplineCurve* cross_cv_ref,
                   bool appr_offset_cv)
//===========================================================================
{
    shared_ptr<const ParamCurve> pcv = cv.parameterCurve();
    shared_ptr<SplineCurve> space_cv(dynamic_cast<SplineCurve*>(cv.spaceCurve()->clone()));
    //    ASSERT(cv.parPref()); // We could of course circumvent by closest point, don't bother though.
    shared_ptr<const ParamSurface> sf = cv.underlyingSurface();
    int nmb_samples = 1000;
    double tmin = pcv->startparam();
    double tmax = pcv->endparam();
    double tstep = (tmax - tmin)/((double)(nmb_samples - 1));
    vector<double> params;
    vector<double> pts;
    // We sample input cv.
    for (int ki = 0; ki < nmb_samples; ++ki) {
        double tpar = tmin + ki*tstep;
        Point par_pt = pcv->point(tpar);
        vector<Point> space_pt = cv.ParamCurve::point(tpar, 1);
        Point sf_normal;
        sf->normal(sf_normal, par_pt[0], par_pt[1]);

        Point cross_pt = sf_normal%space_pt[1];
        // If we're to use different basis we must calculate suitable parameter value.
        // Expecting cv & basis_space_cv to correspond.
        if (basis_space_cv.get() != 0) {
            double clo_t, clo_dist;
            Point clo_pt;
            basis_space_cv->closestPoint(space_pt[0], basis_space_cv->startparam(),
                                         basis_space_cv->endparam(),
                                         clo_t, clo_pt, clo_dist);
            params.push_back(clo_t);
            if (cross_cv_ref != 0) {
                // We turn direction of cross tangent.
                Point cross_dir = cross_cv_ref->ParamCurve::point(clo_t);
                double angle = space_pt[1].angle(cross_dir);
                // Now we want to turn space_pt[1] angle degrees in plane defined by space_pt[1] & cross_pt
                Point new_cross_dir = space_pt[1];
                GeometryTools::rotatePoint(sf_normal, angle, &new_cross_dir[0]);

// 		double coef1, coef2;
// 		int dim = 3;
// 		int sign = 1;
// 		CoonsPatchGen::blendcoef(&space_pt[1][0], &cross_pt[0], &cross_dir[0],
// 					   dim, sign, &coef1, &coef2);
// 		Point proj_dir = coef1*space_pt[1] + coef2*cross_pt;
// 		cross_pt = proj_dir;

                cross_pt = new_cross_dir;
                // A better approach would be to define calculate angle between input cross der
                // and bd_cv tangent, and then turn tangent in plane defined by sf normal.
                // Otherwise we might get into trouble if input cross tangent is almost
                // parallell with normal, or even worse to the "left".

            }
        } else {
            params.push_back(tpar);
        }

        double offset_dist;
        if (cross_cv_ref != 0) { // Expecting cross_cv_ref to share domain of corr cv.
            Point scale_pt = cross_cv_ref->ParamCurve::point(params.back());
            offset_dist = scale_pt.length();
        } else { // Using length of directional deriv in input sf.
            vector<Point> sf_pt(3);
            sf->point(sf_pt, par_pt[0], par_pt[1], 1);
            offset_dist = sqrt(sf_pt[1].length() + sf_pt[2].length());	    
        }
        cross_pt.normalize();
        cross_pt *= offset_dist;
        // Hmm, we will often want the cv to have a direction as given in another sf.
        // Smoothing routine more stable using offsets.
        Point appr_pt = (appr_offset_cv) ? space_pt[0] + cross_pt : cross_pt;
        pts.insert(pts.end(), appr_pt.begin(), appr_pt.end());
    }

    // We then perform smoothing on the edge cv.
    int dim = 3;
    SmoothCurve smooth_cv(dim);
//     int cont = 1;
    //int seem = 0; //min(2, cont + 1); // Want to maintain continuity towards end pieces. At most C1 curr.

    vector<int> coef_known(space_cv->numCoefs(), 0);
    if (basis_space_cv.get() != 0) {
        //smooth_cv.attach(basis_space_cv, seem, &coef_known[0]);
        smooth_cv.attach(basis_space_cv, &coef_known[0]);
    } else {
        //smooth_cv.attach(space_cv, seem, &coef_known[0]);
      shared_ptr<SplineCurve> space_cv2(space_cv->clone());
        smooth_cv.attach(space_cv2, &coef_known[0]);
    }

    double smoothweight = 0.000000001;
    double smoothfac = 1.0;
    double wgt1 = 0.000001*smoothweight;
    double wgt2 = smoothweight/smoothfac;
    double wgt3 = smoothweight/(smoothfac*smoothfac);
    double approx_weight = 1.0 - (wgt1 + wgt2 + wgt3);
    double sum = wgt1 + wgt2 + wgt3 + approx_weight;
    wgt1 /= sum;
    wgt2 /= sum;
    wgt3 /= sum;
    approx_weight /= sum;
    smooth_cv.setOptim(wgt1, wgt2, wgt3);

    vector<double> pt_weights(params.size(), 1.0);
    // We also attach sample points from the edge cv.
    smooth_cv.setLeastSquares(pts, params, pt_weights, approx_weight);

    // Possibly include sample points?
    shared_ptr<SplineCurve> appr_cv;
    smooth_cv.equationSolve(appr_cv);

    shared_ptr<SplineCurve> return_cv(new SplineCurve(*appr_cv));
    if (appr_offset_cv) {
        int num_coefs = appr_cv->numCoefs();
        for (int ki = 0; ki < num_coefs*dim; ++ki) {
            return_cv->coefs_begin()[ki] =
                appr_cv->coefs_begin()[ki] - space_cv->coefs_begin()[ki];
        }
    }

    return return_cv;
}


//===========================================================================
vector<Point>
CreatorsUtils::projectPoint(const ParamSurface* sf,
                            bool closed_dir_u, bool closed_dir_v,
                            const Point& space_pt, double epsgeo)
//===========================================================================
{
    vector<Point> pts;

    const Point sf_epspar = SurfaceTools::getParEpsilon(*sf, epsgeo);
    const double epspar = std::min(sf_epspar[0], sf_epspar[1]);

    // We project the cv_pt onto sf.
    double upar, vpar;
    Point clo_pt;
    double clo_dist;
    sf->closestPoint(space_pt, upar, vpar, clo_pt, clo_dist, epsgeo);
    pts.push_back(Point(upar, vpar));

    if (sf->instanceType() == Class_SplineSurface) {
        const SplineSurface* spline_sf = dynamic_cast<const SplineSurface*>(sf);
        // If the surface is closed and the point is on the boundary, we
        // add more points. We possibly get 4 parameter points (most likely only 2).
        double umin = spline_sf->startparam_u();
        double umax = spline_sf->endparam_u();
        double vmin = spline_sf->startparam_v();
        double vmax = spline_sf->endparam_v();
        bool onLeftEdge = (upar - umin < epspar);
        bool onRightEdge = (umax - upar < epspar);
        bool onLowerEdge = (vpar - vmin < epspar);
        bool onUpperEdge = (vmax - vpar < epspar);
        if (closed_dir_u) {
            if (onLeftEdge)
                pts.push_back(Point(umax, vpar));
            else if (onRightEdge)
                pts.push_back(Point(umin, vpar));
        }
        if (closed_dir_v) {
            if (onLowerEdge)
                pts.push_back(Point(upar, vmax));
            else if (onUpperEdge)
                pts.push_back(Point(upar, vmin));
        }
        if (closed_dir_u && closed_dir_v && (pts.size() == 3))
            pts.push_back(Point(pts[1][0], pts[2][1]));
    }

    return pts;
}


//===========================================================================
shared_ptr<Point>
CreatorsUtils::projectCurvePoint(const ParamSurface* sf,
                                 bool closed_dir_u, bool closed_dir_v,
                                 const ParamCurve* space_cv, double cv_par, double epsgeo)
//===========================================================================
{
    if (sf->instanceType() == Class_SplineSurface) {
        const SplineSurface* spline_sf = dynamic_cast<const SplineSurface*>(sf);
        return projectCurvePoint(*spline_sf,
                                 closed_dir_u, closed_dir_v,
                                 space_cv, cv_par, epsgeo);
    } else {
        // Ordinary closest point. Should be extended to handle seem if
        // sf is periodic (for cylinder & cone for instance).
        Point space_pt = space_cv->point(cv_par);
        double upar, vpar;
        Point clo_pt;
        double clo_dist;
        sf->closestPoint(space_pt, upar, vpar, clo_pt, clo_dist, epsgeo);
        return shared_ptr<Point>(new Point(upar, vpar));
    }
}


//===========================================================================
shared_ptr<Point>
CreatorsUtils::projectCurvePoint(const SplineSurface& sf,
                                 bool closed_dir_u, bool closed_dir_v,
                                 const ParamCurve* space_cv, double cv_par,
				 double epsgeo)
//===========================================================================
{
    double angtol = 1e-05;
    const Point sf_epspar = SurfaceTools::getParEpsilon(sf, epsgeo);
    const double epspar = std::min(sf_epspar[0], sf_epspar[1]);

    double deg_tol = epsgeo;
    bool b_deg, r_deg, t_deg, l_deg;
    // If surface is degenerate special handling is required.
    sf.isDegenerate(b_deg, r_deg, t_deg, l_deg, deg_tol);
    // We're assuming (as an easy start, which covers most cases) that
    // this only occurs for an end parameter on the curve.
    double tmin = space_cv->startparam();
    double tmax = space_cv->endparam();
    bool at_start = (cv_par - tmin < epspar);
    bool at_end = (tmax - cv_par < epspar);
    if ((!at_start) && (!at_end) &&
        (space_cv->instanceType() == Class_SplineCurve)) {
        const SplineCurve* spline_cv =
            dynamic_cast<const SplineCurve*>(space_cv);
        int order = spline_cv->order();
        int multi = spline_cv->basis().knotMultiplicity(cv_par);
        if (order - multi <= 1) {
            vector<Point> left(2), right(2);
            int nderivs = 1;
            bool from_right = true;
            spline_cv->point(left, cv_par, nderivs, !from_right);
            spline_cv->point(right, cv_par, nderivs, from_right);
            double kink = left[1].angle(right[1]);
            if (kink >= angtol) {
                THROW("SplineCurve has a kink for this parameter value");
            }
        }
    }

    // We're of course assuming that we're not to cross the
    // seem of a closed surface.
    // Evaluating the tangent of the space curve from both directions.
    vector<Point> cv_pt(2);
    bool eval_from_right = at_start;
    space_cv->point(cv_pt, cv_par, 1, eval_from_right);

    // Get all possible projected candidates
    vector<Point> cand_par_pts = projectPoint(&sf, closed_dir_u, closed_dir_v,
                                              cv_pt[0]);

    double umin = sf.startparam_u();
    double umax = sf.endparam_u();
    double vmin = sf.startparam_v();
    double vmax = sf.endparam_v();
    if (cand_par_pts.size() == 1)
      {  // Only one point - this must be it
        double upar = cand_par_pts[0][0];
        double vpar = cand_par_pts[0][1];
        // If we are along a degenerate edge we cannot trust the
        // parameter value (as there are an infinite number of correct
        // parameter values).
        if ((b_deg && (vpar - vmin < deg_tol)) ||
            (t_deg && (vmax - vpar < deg_tol)) ||
            (l_deg && (upar - umin < deg_tol)) ||
            (r_deg && (umax - upar < deg_tol))) {
            if (at_start || at_end)
            {
                // We try to fix things by stepping slightly away from
                // the end pt.
                double tpar2 = (at_start) ? tmin + 0.1*(tmax - tmin) :
                    tmax - 0.1*(tmax - tmin);
                vector<Point> cv_pt2(2);
                space_cv->point(cv_pt2, tpar2, 1, eval_from_right);
                vector<Point> cand_par_pts2 =
                    projectPoint(&sf, closed_dir_u, closed_dir_v,
                                 cv_pt2[0]);
                if (cand_par_pts2.size() == 0)
                {
                    MESSAGE("Failed!");
                    return shared_ptr<Point>();
                }
                else
                {
                    if (cand_par_pts2.size() > 1)
                        MESSAGE("Did not expect this ... Trying anyway.");
                    // We suspect that both cv pts lie on an iso
                    // curve, at least locally.
                    if (b_deg || t_deg)
                        upar = cand_par_pts2[0][0];
                    else
                        vpar = cand_par_pts2[0][1];
                    // We must then see if pt coincides with cv pt.
                    Point sf_pt = sf.ParamSurface::point(upar, vpar);
                    double pt_dist = sf_pt.dist(cv_pt[0]);
                    if (pt_dist > epsgeo)
                        MESSAGE("Projection seems to be inaccurate.");
                    return shared_ptr<Point>(new Point(upar,
                                                       vpar));
                }
            }
            else
            {
                // Currently not handling interior points along
                // degenerate edges.
                MESSAGE("Degenerate case - returning empty point.");	  
                return shared_ptr<Point>();
            }
        }
        else
          return shared_ptr<Point>(new Point(cand_par_pts[0][0],
                                             cand_par_pts[0][1]));
      }

    // We then check the directional derivs in the params. I.e. we
    // compare the tangent of the curve with the tangent of the
    // boundary curve.
    int index = 0;
    int nhits = 0;
    vector<Point> sf_pt(3);
    for (size_t ki = 0; ki < cand_par_pts.size(); ++ki) {
        double upar = cand_par_pts[ki][0];
        double vpar = cand_par_pts[ki][1];

        bool onLeftEdge = (upar - umin < sf_epspar[0]);
        bool onRightEdge = (umax - upar < sf_epspar[0]);
        bool onLowerEdge = (vpar - vmin < sf_epspar[1]);
        bool onUpperEdge = (vmax - vpar < sf_epspar[1]);

        if (!((onLeftEdge && l_deg) ||
              (onRightEdge && r_deg) ||
              (onLowerEdge && b_deg) ||
              (onUpperEdge && t_deg))) {
            // We are not in a degenerate point.

            // Surface point and derivs
            sf.point(sf_pt, upar, vpar, 1);

            // Check if the curve follows the seam and is counterclockwise
            Point cv_tan = cv_pt[1];
            double smallangu = cv_tan.angle_smallest(sf_pt[1]);
            double smallangv = cv_tan.angle_smallest(sf_pt[2]);
            bool onSeam = false;
            if (smallangv < angtol && closed_dir_u)
                onSeam = true;
            if (smallangu < angtol && closed_dir_v)
                onSeam = true;
            if (onSeam) {
                // If the point on the curve is on the seam we have an
                // undecided situation. However, we use a trick to split
                // the tie: We rotate the curve tangent by +/- 45 degrees
                // around the surface normal... @jbt
                Point rot_axis;
                sf.normal(rot_axis, upar, vpar);
                double rot_angle = 0.0;
                if (at_start)
                    rot_angle = 0.25 * M_PI;
                else if (at_end)
                    rot_angle = -0.25 * M_PI;
                GeometryTools::rotatePoint(rot_axis, rot_angle, &cv_tan[0]);

// 	    MESSAGE("Udecidable case - curve follows seam at this parameter. "
// 		    "Assuming the curve is counterclockwise.");
            }

            double cosangu = cv_tan.cosAngle(sf_pt[1]);
            double cosangv = cv_tan.cosAngle(sf_pt[2]);

            double c1, c2;
            CoonsPatchGen::blendcoef(&sf_pt[1][0], &sf_pt[2][0], &cv_tan[0],
                                     3, 1, &c1, &c2);

            // Does the curve point inwards or outwards?

            // Check corners
            if (onLeftEdge && onLowerEdge) {
                //	    if (at_start && cosangu > -angtol && cosangv > -angtol) {
                if (at_start && c1 > -angtol && c2 > -angtol) {
                    index = (int)ki;
                    ++nhits;
                }
                //if (at_end && cosangu < angtol && cosangv < angtol) {
                if (at_end && c1 < angtol && c2 < angtol) {
                    index = (int)ki;
                    ++nhits;
                }
            }
            else if (onRightEdge && onLowerEdge) {
                //	    if (at_start && cosangu < angtol && cosangv > -angtol) {
                if (at_start && c1 < angtol && c2 > -angtol) {
                    index = (int)ki;
                    ++nhits;
                }
                //	    if (at_end && cosangu > -angtol && cosangv < angtol) {
                if (at_end && c1 > -angtol && c2 < angtol) {
                    index = (int)ki;
                    ++nhits;
                }
            }
            else if (onRightEdge && onUpperEdge) {
                //	    if (at_start && cosangu < angtol && cosangv < angtol) {
                if (at_start && c1 < angtol && c2 < angtol) {
                    index = (int)ki;
                    ++nhits;
                }
                //	    if (at_end && cosangu > -angtol && cosangv > -angtol) {
                if (at_end && c1 > -angtol && c2 > -angtol) {
                    index = (int)ki;
                    ++nhits;
                }
            }
            else if (onLeftEdge && onUpperEdge) {
                //	    if (at_start && cosangu > -angtol && cosangv < angtol) {
                if (at_start && c1 > -angtol && c2 < angtol) {
                    index = (int)ki;
                    ++nhits;
                }
                //	    if (at_end && cosangu < angtol && cosangv > -angtol) {
                if (at_end && c1 < angtol && c2 > -angtol) {
                    index = (int)ki;
                    ++nhits;
                }
            }
            // Check edges
            else if (onLeftEdge) {
                if (at_start && cosangu > -angtol) {
                    index = (int)ki;
                    ++nhits;
                }
                if (at_end && cosangu < angtol) {
                    index = (int)ki;
                    ++nhits;
                }
            }
            else if (onLowerEdge) {
                if (at_start && cosangv > -angtol) {
                    index = (int)ki;
                    ++nhits;
                }
                if (at_end && cosangv < angtol) {
                    index = (int)ki;
                    ++nhits;
                }
            }
            else if (onRightEdge) {
                if (at_start && cosangu < angtol) {
                    index = (int)ki;
                    ++nhits;
                }
                if (at_end && cosangu > -angtol) {
                    index = (int)ki;
                    ++nhits;
                }
            }
            else if (onUpperEdge) {
                if (at_start && cosangv < angtol) {
                    index = (int)ki;
                    ++nhits;
                }
                if (at_end && cosangv > -angtol) {
                    index = (int)ki;
                    ++nhits;
                }
            }

        } else {
            if (at_start || at_end)
            {
                // We try to fix things by stepping slightly away from
                // the end pt.
                double tpar2 = (at_start) ? tmin + 0.1*(tmax - tmin) :
                    tmax - 0.1*(tmax - tmin);
                vector<Point> cv_pt2(2);
                space_cv->point(cv_pt2, tpar2, 1, eval_from_right);
                vector<Point> cand_par_pts2 =
                    projectPoint(&sf, closed_dir_u, closed_dir_v,
                                 cv_pt2[0]);
                if (cand_par_pts2.size() == 0)
                {
                    MESSAGE("Failed!");
                    return shared_ptr<Point>();
                }
                else
                {
                    if (cand_par_pts2.size() == 1)
                    {
                        Point curr_par_pt = cand_par_pts[index];
                        if (b_deg || t_deg)
                            cand_par_pts[ki][0] = cand_par_pts2[0][0];
                        else
                            cand_par_pts[ki][1] = cand_par_pts2[0][1];
                        if (nhits > 0) {
                            double distu =
                                curr_par_pt[0] - cand_par_pts[ki][0];
                            double distv =
                                curr_par_pt[1] - cand_par_pts[ki][1];
                            if (distu + distv > epspar)
                                ++nhits;
                        } else{
                            ++nhits;
                        }
                        index = (int)ki;
                    }
                    else
                    {
                        // If we got here we should expect curve to lie on
                        // a seem (i.e. the surface is closed).
                        if ((!closed_dir_u) && (!closed_dir_v))
                            MESSAGE("Unexpected ... Trying anyway.");

                        // Since the space cv is expected to be a trim
                        // curve, it should have a ccw orientation.

                        // Multiple parameter values found. We slightly
                        // alter the parameter value not lying on the deg
                        // edge, and then choose depending on whether the
                        // cross product points into the surface.
                        // double min_dist = -1.0;
                        for (size_t kj = 0; kj < cand_par_pts2.size(); ++kj)
                        {
                            double upar2 = cand_par_pts2[kj][0];
                            double vpar2 = cand_par_pts2[kj][1];
                            int sign = 1;
// 				double step_frac = 0.01;
                            if (at_end)
                                sign = -1;
// 				if (l_deg || r_deg)
// 				    upar2 += sign*step_frac*
// 					(sf.endparam_u() - sf.startparam_u());
// 				else
// 				    vpar2 += sign*step_frac*
// 					(sf.endparam_v() - sf.startparam_v());
                            Point alt_pt =
                                sf.ParamSurface::point(upar2, vpar2);
                            Point clo_pt;
                            double clo_t, clo_dist;
                            space_cv->closestPoint(alt_pt,
                                                   space_cv->startparam(),
                                                   space_cv->endparam(),
                                                   clo_t, clo_pt, clo_dist);

                            vector<Point> alt_cv_pt =
                                space_cv->ParamCurve::point(clo_t, 1);
                            Point sf_normal;
                            sf.normal(sf_normal, upar2, vpar2);
                            Point cv_cross_prod = sf_normal%alt_cv_pt[1];
                            vector<Point> sf_pts =
                                sf.ParamSurface::point(upar2, vpar2, 1);
                            Point sf_dir_in =
                                (closed_dir_u) ? sf_pts[1] : sf_pts[2];
                            bool curr_at_end = (closed_dir_u) ?
                                (fabs(sf.endparam_u() - upar2) < sf_epspar[0]) :
                                (fabs(sf.endparam_v() - vpar2) < sf_epspar[1]);
                            if (curr_at_end)
                                sf_dir_in *= -1.0;
                            // The angle should be quite small, but
                            // 0.5*M_PI should suffice as a cutoff value.
                            double angle = cv_cross_prod.angle(sf_dir_in);
                            if (angle < 0.5*M_PI)
                            {
                                Point curr_par_pt = cand_par_pts[index];
                                if (b_deg || t_deg)
                                    cand_par_pts[ki][0] = cand_par_pts2[kj][0];
                                else
                                    cand_par_pts[ki][1] = cand_par_pts2[kj][1];
                                if (nhits > 0) {
                                    double distu =
                                        curr_par_pt[0] - cand_par_pts[ki][0];
                                    double distv =
                                        curr_par_pt[1] - cand_par_pts[ki][1];
                                    if (distu + distv > epspar)
                                        ++nhits;
                                } else{
                                    ++nhits;
                                }
                                index = (int)ki;
                            }
                        }
                    }
                    
                }
            }
            else
            {
                MESSAGE("Degenerate case. Unexpected. Trying anyway.");
            }
        }
    }
// 	return shared_ptr<Point>(new Point(upar,
// 					   vpar));


    if (nhits != 1) {
        MESSAGE("Undecidable case - number of candidates != 1. "
                "Returning arbitrary point.");
        MESSAGE("nhits = " << nhits);
    }

    double upar = cand_par_pts[index][0];
    double vpar = cand_par_pts[index][1];

    // We see if pt coincides with cv pt.
    Point sf_pt2 = sf.ParamSurface::point(upar, vpar);
    double pt_dist = sf_pt2.dist(cv_pt[0]);
    if (pt_dist > epsgeo)
        MESSAGE("Projection seems to be inaccurate.");

//     // Debug
//     ofstream sffile("sf.g2");
//     sf.writeStandardHeader(sffile);
//     sf.write(sffile);
//     ofstream cvfile("space_cv.g2");
//     space_cv.writeStandardHeader(cvfile);
//     space_cv.write(cvfile);
//     cout << "closed_dir_u = " << closed_dir_u << endl
// 	 << "closed_dir_v = " << closed_dir_v << endl
// 	 << "cv_par = " << cv_par << endl
// 	 << "upar = " << upar << endl
// 	 << "vpar = " << vpar << endl;

    return shared_ptr<Point>(new Point(upar, vpar));
}


//===========================================================================
void
CreatorsUtils::fixSeemCurves(shared_ptr<BoundedSurface> bd_sf, 
                             vector<shared_ptr<CurveOnSurface> >& loop_cvs,
                             bool closed_dir_u, bool closed_dir_v,
                             double tol)
//===========================================================================
{
  // Compute period
  RectDomain dom = bd_sf->underlyingSurface()->containingDomain();
  double start[2], end[2];
  start[0] = dom.umin();
  end[0] = dom.umax();
  start[1] = dom.vmin();
  end[1] = dom.vmax();

  // For each loop curve, check if it lies along a seem. Mark that
  // there might be missing degenerate curves in the loop.
  vector<int> seem(loop_cvs.size(), 0);
  size_t ki, kj;
  for (ki=0; ki<loop_cvs.size(); ++ki)
    {
      shared_ptr<ParamCurve> par_cv = loop_cvs[ki]->parameterCurve();
      if (!par_cv.get())
        continue;

      BoundingBox pbox = par_cv->boundingBox();
      Point high = pbox.high();
      Point low = pbox.low();
      if (closed_dir_u && high[0]-low[0] < tol)
        {
          // Possibility for seem curve with constant u parameter
          if (fabs(low[0]-start[0]) < tol || fabs(high[0]-end[0]) < tol)
            seem[ki] += 1;
        }

      if (closed_dir_v && high[1]-low[1] < tol)
        {
          // Possibility for seem curve with constant v parameter
          if (fabs(low[1]-start[1]) < tol || fabs(high[1]-end[1]) < tol)
            seem[ki] += 2;
        }
    }

  // Check if a curve along the seem must be moved to the other side of the seem
  for (ki=0; ki<loop_cvs.size(); ++ki)
    {
      kj = (ki+1) % loop_cvs.size();
      if ((seem[ki] && !seem[kj]) || (seem[kj] && !seem[ki]))
        {
          // Check if the parameter and geometry information
          // match at the joint between the curves
          shared_ptr<ParamCurve> par_cv1 = loop_cvs[ki]->parameterCurve();
          shared_ptr<ParamCurve> par_cv2 = loop_cvs[kj]->parameterCurve();
          if (!par_cv1.get() || !par_cv2.get())
            continue;   // No parameter curve

          Point pos1 =  loop_cvs[ki]->ParamCurve::point(loop_cvs[ki]->endparam());
          Point pos2 =  loop_cvs[kj]->ParamCurve::point(loop_cvs[kj]->startparam());
          Point par1 = par_cv1->ParamCurve::point(par_cv1->endparam());
          Point par2 = par_cv2->ParamCurve::point(par_cv2->startparam());
          if (pos1.dist(pos2) < tol && par1.dist(par2) > tol)
            {
              // The seem curve must be moved to match the other curve
              if (seem[ki])
                {
                  shared_ptr<SplineCurve> cv = dynamic_pointer_cast<SplineCurve>(par_cv1);
                  int idx = (seem[ki] == 1) ? 0 : 1;
                  if (cv.get() && (fabs(par2[idx]-start[idx]) < tol ||
                                   fabs(par2[idx]-end[idx]) < tol))
                    {
                      vector<double>::iterator cstart = cv->coefs_begin();
                      vector<double>::iterator cend = cv->coefs_end();
                      for (; cstart != cend; cstart+=2)
                        cstart[idx] = par2[idx];
                    }
                }
              else if (seem[kj])
                {
                  shared_ptr<SplineCurve> cv = dynamic_pointer_cast<SplineCurve>(par_cv2);
                  int idx = (seem[kj] == 1) ? 0 : 1;
                  if (cv.get() && (fabs(par1[idx]-start[idx]) < tol ||
                                   fabs(par1[idx]-end[idx]) < tol))
                    {
                      vector<double>::iterator cstart = cv->coefs_begin();
                      vector<double>::iterator cend = cv->coefs_end();
                      for (; cstart != cend; cstart+=2)
                        cstart[idx] = par1[idx];
                    }
                }
            }
        }
    }
          
}

//===========================================================================
void
CreatorsUtils::fixTrimCurves(shared_ptr<Go::BoundedSurface> bd_sf,
                             double epsgeo_frac, double gap, double neighbour,
			     double kink)
//===========================================================================
{
  //    MESSAGE("Method in an alpha state ...");

  // First remove small curves in trimming loops
  bd_sf->removeSmallBoundaryCurves(gap, neighbour, kink);

    int kj, kk, kh;
    int nmb_loops = bd_sf->numberOfLoops();
    shared_ptr<ParamSurface> under_sf = bd_sf->underlyingSurface();
    shared_ptr<ParamSurface> orig_under_sf = under_sf;
    double deg_tol = 1e-04;
    bool bottom, right, top, left;
    under_sf->isDegenerate(bottom, right, top, left, deg_tol);
    bool deg_sf = (bottom || right || top || left);
    bool closed_dir_u = false, closed_dir_v = false;
    try {
      SurfaceTools::checkSurfaceClosed(*under_sf, closed_dir_u, closed_dir_v, deg_tol);
    } catch (...) {
        MESSAGE("Failed determining if surface is closed. ");
        return;
    }


    //if (under_sf->instanceType() == Class_SplineSurface) {
    //    try {
    //        shared_ptr<SplineSurface> spline_sf =
    //            dynamic_pointer_cast<SplineSurface, ParamSurface>(under_sf);
    //        // Making sure the sf is k-regular.
    //        shared_ptr<SplineSurface> spline_under_sf =
    //            shared_ptr<SplineSurface>
    //            (spline_sf->subSurface(spline_sf->startparam_u(),
    //                                   spline_sf->startparam_v(),
    //                                   spline_sf->endparam_u(),
    //                                   spline_sf->endparam_v()));
    //        surfaceClosed(*spline_under_sf,
    //                      closed_dir_u, closed_dir_v, deg_tol);
    //    } catch (...) {
    //        MESSAGE("Failed determining if surface is closed. ");
    //        return;
    //    }
    //}

    // We make sure that all trim curves have matching parameter and
    // space curves.
//     double epsgeo = 1e-02;//4;//2;//4;
    for (kj = 0; kj < nmb_loops; ++kj) {
        shared_ptr<CurveLoop> loop = bd_sf->loop(kj);
        //double epsgeo = std::max(tol, epsgeo_frac*loop->getSpaceEpsilon());
        double epsgeo = gap;
        shared_ptr<Point> loop_seem_par_pt, prev_end_par_pt;
        // We first run though all cvs, locating the start and end pts.
        int loop_size = loop->size();
        vector<shared_ptr<Point> > start_pts(loop_size), end_pts(loop_size);
        vector<shared_ptr<CurveOnSurface> > loop_cvs;
        for (kk = 0; kk < loop_size; ++kk) {
	  double len = (*loop)[kk]->estimatedCurveLength();
	  if (len < neighbour)
	    {
	      MESSAGE("Short curve in fixTrimCurves. Length = " << len);
	      int kk1 = kk - 1;
	      if (kk1 < 0)
		kk1 = loop_size - 1;
	      int kk2 = (kk + 1) % loop_size;

	      // Check angles
	      vector<Point> pt1(2), pt2(2), pt3(2), pt4(2);
	      (*loop)[kk1]->point(pt1, (*loop)[kk1]->endparam(), 1);
	      (*loop)[kk]->point(pt2, (*loop)[kk]->startparam(), 1);
	      (*loop)[kk]->point(pt3, (*loop)[kk]->endparam(), 1);
	      (*loop)[kk2]->point(pt4, (*loop)[kk2]->startparam(), 1);

	      double ang1 = pt1[1].angle(pt2[1]);
	      double ang2 = pt3[1].angle(pt4[1]);
	      double d1 = pt1[0].dist(pt2[0]);
	      double d2 = pt3[0].dist(pt4[0]);
	      int stop_break = 1.0;
	    }
            shared_ptr<CurveOnSurface> cv_on_sf =
                dynamic_pointer_cast<CurveOnSurface>((*loop)[kk]);
            ASSERT(cv_on_sf.get() != NULL);
            // To make sure that the orientation and parameter
            // interval is consistent
            cv_on_sf->makeCurvesConsistent(cv_on_sf->parPref());

            loop_cvs.push_back(cv_on_sf);
// 	    shared_ptr<SplineCurve> space_cv =
// 		dynamic_pointer_cast<SplineCurve>(loop_cvs[kk]->spaceCurve());
// 	    shared_ptr<SplineCurve> par_cv =
// 		dynamic_pointer_cast<SplineCurve>
// 		(loop_cvs[kk]->parameterCurve());
            shared_ptr<ParamCurve> space_cv = loop_cvs[kk]->spaceCurve();
            shared_ptr<ParamCurve> par_cv = loop_cvs[kk]->parameterCurve();
            // We compare dists on par & space cvs.
            if ((par_cv.get() != NULL) && (space_cv.get() != NULL)) {
                vector<double> p_pars(2), s_pars(2);//3);
                p_pars[0] = par_cv->startparam();
                p_pars[1] = par_cv->endparam();
                s_pars[0] = space_cv->startparam();
                s_pars[1] = space_cv->endparam();
// 		    tpars[1] = 0.5*(tpars[0] + tpars[2]);
                int max_dist_ind = -1;
                double max_dist = -1.0;
                for (int kl = 0; kl < (int)p_pars.size(); ++kl) {
                    Point par_pt = par_cv->ParamCurve::point(p_pars[kl]);
                    Point cv_pt = space_cv->ParamCurve::point(s_pars[kl]);
                    Point sf_pt = under_sf->ParamSurface::point
                        (par_pt[0], par_pt[1]);
                    double dist = cv_pt.dist(sf_pt);
                    if (dist > max_dist) {
                        max_dist = dist;
                        max_dist_ind = kl;
                    }
                }
                if (max_dist > epsgeo) {
                    // We're assuming that space curves are
                    // correct, except for the case with deg sfs.
#ifdef SBR_DBG
                    MESSAGE("loop: " << kj << ", cv: " << kk << ", max dist in pt " <<
                            max_dist_ind << ": " << max_dist);
#endif
                    if (deg_sf) {
                        loop_cvs[kk] = shared_ptr<CurveOnSurface>
                            (new CurveOnSurface(orig_under_sf, par_cv,
                                                shared_ptr<ParamCurve>(),
                                                true));
                    } else {
                        loop_cvs[kk] = shared_ptr<CurveOnSurface>
                            (new CurveOnSurface(orig_under_sf,
                                                shared_ptr<ParamCurve>(),
                                                space_cv, false));
                    }
                }
            }
            if (space_cv.get() == NULL) {
                shared_ptr<ParamCurve> par_cv = loop_cvs[kk]->parameterCurve();
// 		    dynamic_pointer_cast<ParamCurve>
// 		ASSERT(par_cv.get() != NULL);
                start_pts[kk] =
                    (shared_ptr<Point>(new Point
                                       (par_cv->ParamCurve::point
                                        (par_cv->startparam()))));
                end_pts[kk] =
                    (shared_ptr<Point>(new Point
                                       (par_cv->ParamCurve::point
                                        (par_cv->endparam()))));
// 	    } else if (!closed_dir_u && !closed_dir_v) {
// 		start_pts.push_back(shared_ptr<Point>());
// 		end_pts.push_back(shared_ptr<Point>());
            } else {
                shared_ptr<SplineCurve> spline_space_cv;
                if (space_cv->instanceType() == Class_SplineCurve)
                    spline_space_cv =
                        dynamic_pointer_cast<SplineCurve, ParamCurve>(space_cv);
                else
                    spline_space_cv = shared_ptr<SplineCurve>
                        (space_cv->geometryCurve());
                shared_ptr<Point> proj_start_pt =
                    projectCurvePoint(under_sf.get(),
                                      closed_dir_u, closed_dir_v,
                                      spline_space_cv.get(),
                                      space_cv->startparam());
                if (proj_start_pt.get() != NULL) {
                    Point sf_start_pt =
                        under_sf->ParamSurface::point((*proj_start_pt)[0],
                                                      (*proj_start_pt)[1]);
                    Point cv_start_pt =
                        space_cv->ParamCurve::point(space_cv->startparam());
                    double start_dist = sf_start_pt.dist(cv_start_pt);
                    if (start_dist > epsgeo)
                      {
//  			start_pts.push_back(shared_ptr<Point>());
#ifdef SBR_DBG
                        MESSAGE("Projected point not too close! "
                                "start_dist: " << start_dist <<
                                ", epsgeo: " << epsgeo);
#endif
                      }
                    // 		    else
                    start_pts[kk] = proj_start_pt;
                } else {
                  start_pts[kk] = proj_start_pt;
                }

                shared_ptr<Point> proj_end_pt =
                    projectCurvePoint(under_sf.get(),
                                      closed_dir_u, closed_dir_v,
                                      spline_space_cv.get(),
                                      space_cv->endparam());
                if (proj_end_pt.get() != NULL) {
                    Point sf_end_pt =
                        under_sf->ParamSurface::point((*proj_end_pt)[0], (*proj_end_pt)[1]);
                    Point cv_end_pt =
                        space_cv->ParamCurve::point(space_cv->endparam());
                    double end_dist = sf_end_pt.dist(cv_end_pt);
                    if (end_dist > epsgeo)
                      {
                        end_pts[kk] = shared_ptr<Point>();
#ifdef SBR_DBG
                        MESSAGE("Projected point not too close! "
                                "end_dist: " << end_dist <<
                                ", epsgeo: " << epsgeo);
#endif
                      }
// 		    else
                    end_pts[kk] = proj_end_pt;
                } else
                  end_pts[kk] = proj_end_pt;
            }
        }

        // Check consistence in parametric endpoints
        RectDomain dom = bd_sf->underlyingSurface()->containingDomain();
        double start_u = dom.umin();
        double end_u = dom.umax();
        double per_u = (closed_dir_u) ? end_u - start_u : 0.0;
        double start_v = dom.vmin();
        double end_v = dom.vmax();
        double per_v = (closed_dir_v) ? end_v - start_v : 0.0;
        if (closed_dir_u || closed_dir_v)
        {
            int at_seem1=0, at_seem2=0;
            bool hit1, hit2;

            if (closed_dir_u)
              {
                hit1 = (start_pts[0].get()) ?
                  (fabs((*start_pts[0])[0]-start_u)<deg_tol || 
                   fabs((*start_pts[0])[0]-end_u)<deg_tol) : true;
                hit2 = (end_pts[0].get()) ?
                  (fabs((*end_pts[0])[0]-start_u)<deg_tol || 
                   fabs((*end_pts[0])[0]-end_u)<deg_tol) : true;
                if (hit1 && hit2)
                  at_seem1 += 1;
              }
            if (closed_dir_v)
              {
                hit1 = (start_pts[0].get()) ? 
                  (fabs((*start_pts[0])[1]-start_v)<deg_tol || 
                   fabs((*start_pts[0])[1]-end_v)<deg_tol) : true;
                hit2 = (end_pts[0].get()) ?
                  (fabs((*end_pts[0])[1]-start_v)<deg_tol || 
                   fabs((*end_pts[0])[1]-end_v)<deg_tol) : true;
                if (hit1 && hit2)
                  at_seem1 += 2;
              }
            for (kk=0; kk<(int)start_pts.size(); ++kk)
            {
                kh = (kk+1) % (int)start_pts.size();

                if (closed_dir_u)
                  {
                    hit1 = (start_pts[kh].get()) ?
                      (fabs((*start_pts[kh])[0]-start_u)<deg_tol || 
                       fabs((*start_pts[kh])[0]-end_u)<deg_tol) : true;
                    hit2 = (end_pts[kh].get()) ?
                      (fabs((*end_pts[kh])[0]-start_u)<deg_tol || 
                       fabs((*end_pts[kh])[0]-end_u)<deg_tol) : true;
                    if (hit1 && hit2)
                      at_seem2 += 1;
                  }
                if (closed_dir_v)
                  {
                    hit1 = (start_pts[kh].get()) ? 
                      (fabs((*start_pts[kh])[1]-start_v)<deg_tol || 
                       fabs((*start_pts[kh])[1]-end_v)<deg_tol) : true;
                    hit2 = (end_pts[kh].get()) ?
                      (fabs((*end_pts[kh])[1]-start_v)<deg_tol || 
                       fabs((*end_pts[kh])[1]-end_v)<deg_tol) : true;
                    if (hit1 && hit2)
                      at_seem2 += 2;
                  }
                double d2 = (start_pts[kh].get() && end_pts[kk].get()) ?
                  end_pts[kk]->dist(*start_pts[kh]) : 0.0;
                if (d2 > deg_tol)
                {
                    // Check consistency at the seem. Move the curve along the seem
                    // to create consistency if such a curve exist
                  if (((at_seem1==1 || at_seem1==3) && start_pts[kh].get() && end_pts[kk].get() &&
                        fabs(fabs((*end_pts[kk])[0]-(*start_pts[kh])[0])-per_u) < deg_tol) ||
                        ((at_seem1==2 || at_seem1==3) && 
                         fabs(fabs((*end_pts[kk])[1]-(*start_pts[kh])[1])-per_v) < deg_tol))
                        end_pts[kk] = start_pts[kh];
                  else if (((at_seem2==1 || at_seem2==3) && start_pts[kh].get() && end_pts[kk].get() && 
                            fabs(fabs((*end_pts[kk])[0]-(*start_pts[kh])[0])-per_u) < deg_tol) ||
                             ((at_seem2==2 || at_seem2==3) && 
                              fabs(fabs((*end_pts[kk])[1]-(*start_pts[kh])[1])-per_v) < deg_tol))
                        start_pts[kh] = end_pts[kk];
                }
                at_seem1 = at_seem2;
                at_seem2 = 0;
            }
        }


        for (kk = 0; kk < (int)loop_cvs.size(); ++kk) {
            shared_ptr<CurveOnSurface> cv_on_sf = loop_cvs[kk];
            ASSERT(cv_on_sf.get() != NULL);
// 	    shared_ptr<SplineCurve> space_cv =
// 		dynamic_pointer_cast<SplineCurve>(cv_on_sf->spaceCurve());
// 	    shared_ptr<SplineCurve> par_cv =
// 		dynamic_pointer_cast<SplineCurve>
// 		(cv_on_sf->parameterCurve());
            shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
            shared_ptr<ParamCurve> par_cv = cv_on_sf->parameterCurve();

            // We project space cvs onto the surface.
            if ((cv_on_sf->spaceCurve().get() != NULL &&
                 cv_on_sf->parameterCurve().get() == NULL) ||
                !cv_on_sf->sameCurve(epsgeo)) {
//		if (1) {//cv_on_sf->parameterCurve() != NULL) {
//  		    if (cv_on_sf->parameterCurve() == NULL) {

                try {
                    ASSERT(space_cv.get() != NULL);
                    shared_ptr<Point> from_par_pt, to_par_pt;
                    int prev_ind = (kk == 0) ?
                        (int)loop_cvs.size() -1 : kk - 1;
                    int next_ind = (kk == (int)loop_cvs.size() -1) ?
                        0 : kk + 1;
                    from_par_pt = start_pts[kk];
                    if (from_par_pt.get() == NULL)
                        from_par_pt = end_pts[prev_ind];
                    to_par_pt = end_pts[kk];
                    if (to_par_pt.get() == NULL)
                        to_par_pt = start_pts[next_ind];
                    if ((from_par_pt.get() == NULL) || (to_par_pt.get() == NULL)) {
                        MESSAGE("Missing start/end par pt, should be no point in "
                                "projecting the whole curve.");
// 			continue;
                    }
                    shared_ptr<ParamCurve> proj_cv;
                    try {
                        double proj_tol = epsgeo;
                        proj_cv = shared_ptr<SplineCurve>
                            (CurveCreators::projectSpaceCurve
                             (space_cv, under_sf,
                              from_par_pt, to_par_pt, proj_tol));
                    } catch (...) {
                        MESSAGE("Failed projecting curve with input tol.");
// 			// We try one more time, with a more relaxed
// 			// tolerance.
// 			double proj_tol = 1e02*epsgeo;
// 			try {
// 			    proj_cv = shared_ptr<SplineCurve>
// 				(CurveCreators::projectSpaceCurve
// 				 (space_cv, under_sf,
// 				  from_par_pt, to_par_pt, proj_tol));
// 			    MESSAGE("Projection ok with larger tolerance!");
// 			} catch (...) {
// 			    MESSAGE("Failed projecting curve with larger tol.");
// 			}
                    }
                    if (proj_cv.get() == NULL) {
                        MESSAGE("Failed projecting curve!");
                        prev_end_par_pt = shared_ptr<Point>();
                    } else {
                        shared_ptr<CurveOnSurface> new_cv_on_sf
                            (shared_ptr<CurveOnSurface>
                             (new CurveOnSurface
                              (cv_on_sf->underlyingSurface(),
                               proj_cv,
                               cv_on_sf->spaceCurve(),
                               false)));
                        *cv_on_sf = *new_cv_on_sf;
                        prev_end_par_pt = shared_ptr<Point>
                            (new Point(proj_cv->ParamCurve::point
                                       (proj_cv->endparam())));
//  				writeSpaceParamCurve(*proj_cv,
//  						     proj_cvs_fileout,
//  						     (double)(ki + 100));
                    }
                } catch (...) {
                    MESSAGE("Failed projecting the space curve!");
                }
            } else {
                ASSERT(cv_on_sf->parameterCurve().get() != NULL);
// 		shared_ptr<SplineCurve> param_cv =
// 		    dynamic_pointer_cast<SplineCurve>
// 		    (cv_on_sf->parameterCurve());
                shared_ptr<ParamCurve> param_cv = cv_on_sf->parameterCurve();
// 		    SplineDebugUtils::writeSpaceParamCurve(*param_cv,
// 					 all_geom_fileout,
// 					 (double)(100 + ki));
                ASSERT(param_cv.get() != NULL);
                prev_end_par_pt = shared_ptr<Point>
                    (new Point(param_cv->ParamCurve::point
                               (param_cv->endparam())));
                if (kk == 0) {
                    loop_seem_par_pt = shared_ptr<Point>
                        (new Point(param_cv->ParamCurve::point
                                   (param_cv->startparam())));
                }
            }
        }

        // Check seem curves
        if (closed_dir_u || closed_dir_v)
          fixSeemCurves(bd_sf, loop_cvs, closed_dir_u, closed_dir_v, epsgeo);

        // Check if any degenerate curve-on-surface curve may be missing
        const double deg_tol = 1e-04;
        for (kk=0; kk < (int)loop_cvs.size(); ++kk)
          {
            kh = (kk + 1)%(int)loop_cvs.size();
            // If one of the curves is already degenerate we do not add another one.
            if (loop_cvs[kk]->isDegenerate(deg_tol) || loop_cvs[kh]->isDegenerate(deg_tol))
            {
                continue;
            }
            shared_ptr<ParamCurve> par_cv1 = loop_cvs[kk]->parameterCurve();
            shared_ptr<ParamCurve> par_cv2 = loop_cvs[kh]->parameterCurve();
            if (!par_cv1.get() || !par_cv2.get())
            {
                continue;   // No parameter curve
            }

            Point pos1 =  loop_cvs[kk]->ParamCurve::point(loop_cvs[kk]->endparam());
            Point pos2 =  loop_cvs[kh]->ParamCurve::point(loop_cvs[kh]->startparam());
            Point par1 = par_cv1->ParamCurve::point(par_cv1->endparam());
            Point par2 = par_cv2->ParamCurve::point(par_cv2->startparam());
            // if (pos1.dist(pos2) < epsgeo && 
            //     (par1.dist(par2) > epsgeo &&
            //      !(closed_dir_u && fabs(fabs(par1[0]-par2[0])-per_u) < epsgeo &&
            //        fabs(par1[1]-par2[1]) < epsgeo) &&
            //      !(closed_dir_v && fabs(fabs(par1[1]-par2[1])-per_v) < epsgeo &&
            //        fabs(par1[0]-par2[0]) < epsgeo) &&
            //      !(closed_dir_u && fabs(fabs(par1[0]-par2[0])-per_u) < epsgeo &&
            //        closed_dir_v && fabs(fabs(par1[1]-par2[1])-per_v) < epsgeo)))
            if (pos1.dist(pos2) < neighbour && 
                (par1.dist(par2) > neighbour &&
                 !(closed_dir_u && fabs(fabs(par1[0]-par2[0])-per_u) < neighbour &&
                   fabs(par1[1]-par2[1]) < neighbour) &&
                 !(closed_dir_v && fabs(fabs(par1[1]-par2[1])-per_v) < neighbour &&
                   fabs(par1[0]-par2[0]) < neighbour) &&
                 !(closed_dir_u && fabs(fabs(par1[0]-par2[0])-per_u) < neighbour &&
                   closed_dir_v && fabs(fabs(par1[1]-par2[1])-per_v) < neighbour)))
              {
                // The use of the tolerance is questionable. 
                // Anyway, we must insert a missing degenerate curve-on-surface curve
                shared_ptr<SplineCurve> pcrv = 
                  shared_ptr<SplineCurve>(new SplineCurve(par1, 0.0, par2, 1.0));
                shared_ptr<SplineCurve> scrv = 
                  shared_ptr<SplineCurve>(new SplineCurve(pos1, 0.0, pos2, 1.0));
                shared_ptr<CurveOnSurface> sf_cv = 
                  shared_ptr<CurveOnSurface>(new CurveOnSurface(loop_cvs[kk]->underlyingSurface(),
                                                            pcrv, scrv, false));
                loop_cvs.insert(loop_cvs.begin()+kh, sf_cv);
                kk++;
              }
          }
        
        // We set the (possibly) modified curves.
        vector<shared_ptr<ParamCurve> > cvs(loop_cvs.begin(),
                                            loop_cvs.end());
        loop->setCurves(cvs);
    }

#ifdef SBR_DBG
    {
        std::ofstream bd_sf_of("tmp/bd_sf.g2");
        writeTrimmedInfo(*bd_sf, bd_sf_of, 0.0);
    }
#endif// SBR_DBG
}
