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


#include "GoTools/compositemodel/EvalOffsetSurface.h"
#include "GoTools/compositemodel/ftChartSurface.h"

#include "GoTools/creators/CreatorsOffsetUtils.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/CurvatureAnalysis.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/creators/CurveCreators.h"

#include <vector>
#include <assert.h>

using std::vector;
using std::pair;

namespace Go
{

    //===========================================================================
    EvalOffsetSurface::EvalOffsetSurface(shared_ptr<ftFaceBase> base_sf,
                                         double offset_dist, double epsgeo)
    //===========================================================================
        : base_sf_(base_sf), offset_dist_(offset_dist), epsgeo_(epsgeo)
    {
        double max_error = -1.0;
        double mean_error = -1.0;
        // For a ftChartSurface a SplineSurface is created.
        // @@sbr201701 We still need to handle issues with direction of partial derivs in underlying patches.
        // As well as evaluations currently expecting all patches are of type SplineSurface.
        ftMessage message = base_sf_->createSurf(max_error, mean_error);
        shared_ptr<ParamSurface> param_sf = base_sf_->surface();
        if (param_sf.get() != 0)
        {
            spline_sf_ = shared_ptr<SplineSurface>(param_sf->asSplineSurface());
        }
    }


    //===========================================================================
    EvalOffsetSurface::~EvalOffsetSurface()
    //===========================================================================
    {
    }


    //===========================================================================
    Point EvalOffsetSurface::eval( double u, double v) const
    //===========================================================================
    {
        Point base_pt = base_sf_->point(u, v);

        Point base_normal = base_sf_->normal(u, v);
        base_normal.normalize();
        
        Point offset_pt = base_pt + base_normal*offset_dist_;
        
        return offset_pt;
    }


    //===========================================================================
    void EvalOffsetSurface::eval( double u, double v, int n, Point der[]) const
    //===========================================================================
    {
        if (n == 0)
        {
            der[0] = eval(u, v);
            return;
        }

        // We're assuming that the underlying surface corresponding to the parameter is a spline surface.
        // Alternatively it may be represented as one.
        shared_ptr<SplineSurface> spline_sf = spline_sf_;
        shared_ptr<SplineSurface> spline_sf_global = spline_sf_;
        shared_ptr<SplineSurface> spline_sf_local = spline_sf_; // For the case with 1 surface in the set.
        // For the case with 1 surface only the spline_sf_ should is the same as the original surface.
        // For the surface set we must project the point onto the corresponding input surface.
        Point epar_global(u, v);

        Point epar_local(2);
        ParamSurface* local_par_sf = findLocalSurface(u, v, epar_local[0], epar_local[1]);
        if (local_par_sf != NULL)
        {
            // @@sbr201704 This seems like a memory leak for cases which is not already a SplineSurface!
            spline_sf = shared_ptr<SplineSurface>(local_par_sf->asSplineSurface());
            spline_sf_local = spline_sf;
        }
        else
        {
            MESSAGE("Failed!");
            return;
        }

        bool surface_set = (spline_sf_local.get() != spline_sf_global.get());

        const int kder = 2; // To compute the twist.
        int ind_u=0;             /* Pointer into knot vector                       */
        int ind_v=0;             /* Pointer into knot vector                       */
        int kstat = 0;
        if (spline_sf_.get() == 0) {
            THROW("Missing support for parametric surface as a SplineSurface!");
        }

        // We must blend the directions of the base_sf_ to match the directions of the spline_sf_.
        vector<Point> offset_pt_local(((kder+1)*(kder+2)/2) + 1); // Derivs & normal in the exact surface.
        vector<Point> base_pt(((kder+1)*(kder+2)/2) + 1); // Derivs & normal.
        OffsetUtils::blend_s1421(spline_sf.get(), offset_dist_, kder, epar_local, ind_u, ind_v,
                                 offset_pt_local, base_pt, &kstat);

        vector<Point> offset_pt_global(((kder+1)*(kder+2)/2) + 1); // Derivs & normal in the approximated surface.
        vector<Point> base_pt_global(((kder+1)*(kder+2)/2) + 1); // Derivs & normal.
        OffsetUtils::blend_s1421(spline_sf_global.get(), offset_dist_, kder, epar_global, ind_u, ind_v,
                                 offset_pt_global, base_pt_global, &kstat);

        if (surface_set) {
//            MESSAGE("MISSING: We must express the local derivs using the derivs in spline_sf_!");
            vector<Point> global_pt(3), local_pt(3);
            spline_sf_global->point(global_pt, u, v, 1);
            spline_sf_local->point(local_pt, epar_local[0], epar_local[1], 1);
#if 0
            double dist_pt = global_pt[0].dist(local_pt[0]);
            double ang_u = global_pt[1].angle(local_pt[1]);
            double ang_v = global_pt[2].angle(local_pt[2]);
            std::cout << "dist_pt: " << dist_pt << ", ang_u: " << ang_u << ", ang_v: " << ang_v << std::endl;
#endif            
            // We express global_pt[1] using local_pt[1] & local_pt[2] (i.e. we project the tangent onto plane).
            double a, b;
            const int dim = spline_sf_->dimension();
            CoonsPatchGen::blendcoef(&local_pt[1][0], &local_pt[2][0],
                                     &global_pt[1][0], dim, 1, &a, &b);

            // We express global_pt[2] using local_pt[1] & local_pt[2] (i.e. we project the tangent onto plane).
            double c, d;
            CoonsPatchGen::blendcoef(&local_pt[1][0], &local_pt[2][0],
                                     &global_pt[2][0], dim, 1, &c, &d);
            
            der[0] = offset_pt_local[0];
            // Since the parametrization of local_sf will not coincide with that of the global_sf we must
            // adjust the tangents in the offset surface. We turn the partial derivatives to coincide in
            // the base surfaces. We also adjust the length, assuming approximately linearity in the
            // offset surface.
            der[1] = a*offset_pt_local[1] + b*offset_pt_local[2];
            der[2] = c*offset_pt_local[1] + d*offset_pt_local[2];
            // @@sbr201704 The tangent plane in the offset surface is correct, but the length of the partial derivs
            // is generally not correct (with a saddle point as a typical exception). We should perform the offset
            // calculation in the parameter domain defined by the iso lines in the global approximating surface.

#if 0
            // Setting the twist vector. The vector is tricky to calculate with the parameter domain
            // defined as the closest point from the approximating global surface.
//            der[3] = Point(0.0, 0.0, 0.0);
            der[3] = offset_pt_global[4]; // Generally the twist in global_sf is better than the zero twist.
            // der[3] = offset_pt[4];
#else
            // @@sbr201706 Setting the twist vector to 0.0. For some self intersecting cases the twist
            // vector may end up on the wrong side, giving a bad starting point for the surface
            // smoothing. In these cases we're better off with the zero twist.
            der[3] = Point(0.0, 0.0, 0.0);
            /* Calculate cross derivative of the surface normal. */
		/*
		 *   The cross derivative of the surface normal is calculated by the 
		 *   expression:
		 *
		 *
		 *   N   = P   x P   + P  x P   + P  x P
		 *    uv    uuv   v     uu  vv     u    uvv
		 */
            // snoruv1 = sder[7]%sder[2];
            // snoruv2 = sder[3]%sder[5];
            // snoruv3 = sder[1]%sder[8];
            // where sder are the derivatives from the original surface.
            // snoruv = snoruv1 + snoruv2 + snoruv3;
            // And then we must normalize ...

#endif
        } else {
            der[0] = offset_pt_local[0];
            der[1] = offset_pt_local[1];
            der[2] = offset_pt_local[2];
            der[3] = offset_pt_local[4]; // The twist vector.
        }
        
        return;
    }


    //===========================================================================
    double EvalOffsetSurface::start_u() const
    //===========================================================================
    {
        RectDomain rect_dom = spline_sf_->containingDomain();
        double start_u = rect_dom.umin();
        
        return start_u;
    }



    //===========================================================================
    double EvalOffsetSurface::start_v() const
    //===========================================================================
    {
        RectDomain rect_dom = spline_sf_->containingDomain();
        double start_v = rect_dom.vmin();
        
        return start_v;
    }


    //===========================================================================
    double EvalOffsetSurface::end_u() const
    //===========================================================================
    {
        RectDomain rect_dom = spline_sf_->containingDomain();
        double end_u = rect_dom.umax();
        
        return end_u;
    }


    //===========================================================================
    double EvalOffsetSurface::end_v() const
    //===========================================================================
    {
        RectDomain rect_dom = spline_sf_->containingDomain();
        double end_v = rect_dom.vmax();
        
        return end_v;
    }


    //===========================================================================
    int EvalOffsetSurface::dim() const
    //===========================================================================
    {
        int dim = spline_sf_->dimension();
        
        return dim;
    }


    //===========================================================================
    bool EvalOffsetSurface::approximationOK(double par_u, double par_v, Point approxpos,
					     double tol1, double tol2) const
    //===========================================================================
    {

        Point eval_pos = eval(par_u, par_v);
        double dist = eval_pos.dist(approxpos);

        bool appr_ok = (dist < tol1);        
        if (dist > tol1)
        {
            ;//std::cout << "dist: " << dist << std::endl;
        }

        //  The spline_sf_ is used for defining the parametrization only, hence it is not relevant for
        //  offset evaluations.

        return appr_ok;
    }


    //===========================================================================
    void EvalOffsetSurface::gridSelfIntersections(const HermiteGrid2D& grid,
                                                  vector<int>& grid_self_intersections,
                                                  vector<double>& radius_of_curv) const
    //===========================================================================
    {
        grid_self_intersections.clear();
        radius_of_curv.clear();
        int num_self_int = 0;

        vector<double> self_int, no_self_int; // Grid points in original surfaces.
        vector<double> self_int_offset, no_self_int_offset; // Offset grid points.

        vector<double> knots_u = grid.getKnots(true);
        vector<double> knots_v = grid.getKnots(false);

        vector<Point> data = grid.getData(); // Ordered row-wise: Pos, der_u, der_v, der_uv.

        /// Return the spatial dimension
        //int dim = grid.dim();
        const int MM = grid.size1();
        const int NN = grid.size2();
        MESSAGE("INFO: data.size(): " << data.size() << ", MM: " << MM << ", NN: " << NN);
        double curv_rad_pos_min = MAXDOUBLE;
        double curv_rad_pos_max = -MAXDOUBLE;
        double curv_rad_neg_min = MAXDOUBLE;
        double curv_rad_neg_max = -MAXDOUBLE;

        for (size_t kj = 0; kj < knots_v.size(); ++kj)
        {
            double vpar = knots_v[kj];
            for (size_t ki = 0; ki < knots_u.size(); ++ki)
            {
                double upar = knots_u[ki];

                Point epar_local(2);
                ParamSurface* local_par_sf = findLocalSurface(upar, vpar, epar_local[0], epar_local[1]);
                if (local_par_sf == NULL)
                {
                    MESSAGE("Failed finding the local surface in the surface set!");
                    continue;
                }
                
                double k1, k2; // Curvature. Negative value => Convex shape and no problem w/ offset self intersection.
                Point d1, d2;
                CurvatureAnalysis::principalCurvatures(*local_par_sf, epar_local[0], epar_local[1], k1, d1, k2, d2);

                // double curv_rad1 = (k1 == 0.0) ? -1.0 : 1.0/k1;
                // double curv_rad2 = (k2 == 0.0) ? -1.0 : 1.0/k2;
                double curv_rad1 = 1.0/k1;
                double curv_rad2 = 1.0/k2;
                if ((curv_rad1 < 0.0) && (curv_rad1 < curv_rad_neg_min))
                    curv_rad_neg_min = curv_rad1;
                if ((curv_rad1 < 0.0) && (curv_rad1 > curv_rad_neg_max))
                    curv_rad_neg_max = curv_rad1;
                if ((curv_rad2 < 0.0) && (curv_rad2 < curv_rad_neg_min))
                    curv_rad_neg_min = curv_rad2;
                if ((curv_rad2 < 0.0) && (curv_rad2 > curv_rad_neg_max))
                    curv_rad_neg_max = curv_rad2;

                if ((curv_rad1 > 0.0) && (curv_rad1 < curv_rad_pos_min))
                    curv_rad_pos_min = curv_rad1;
                if ((curv_rad1 > 0.0) && (curv_rad1 > curv_rad_pos_max))
                    curv_rad_pos_max = curv_rad1;
                if ((curv_rad2 > 0.0) && (curv_rad2 < curv_rad_pos_min))
                    curv_rad_pos_min = curv_rad2;
                if ((curv_rad2 > 0.0) && (curv_rad2 > curv_rad_pos_max))
                    curv_rad_pos_max = curv_rad2;

                // std::cout << "curv_rad1: " << curv_rad1 << ", curv_rad2: " << curv_rad2 << std::endl;
                Point offset_pt = data[(kj*MM+ki)*4];
                Point local_sf_pt = local_par_sf->point(epar_local[0], epar_local[1]); // Sf pt from which we offset.
                bool negative_offset = (offset_dist_ < 0.0);
                if ((!negative_offset && ((curv_rad1 > 0.0 && curv_rad1 < offset_dist_) ||
                                          (curv_rad2 > 0.0 && curv_rad2 < offset_dist_))) ||
                    (negative_offset && ((curv_rad1 < 0.0 && curv_rad1 > offset_dist_) ||
                                         (curv_rad2 < 0.0 && curv_rad2 > offset_dist_))))
                {
                    //std::cout << "curv_rad1: " << curv_rad1 << ", curv_rad2: " << curv_rad2 << std::endl;
                    //self_int.insert(self_int.end(), offset_pt.begin(), offset_pt.end());
                    self_int.insert(self_int.end(), local_sf_pt.begin(), local_sf_pt.end());
                    self_int_offset.insert(self_int_offset.end(), offset_pt.begin(), offset_pt.end());
                    //++num_self_int;
                    grid_self_intersections.push_back(kj*MM+ki);
                    double curv_rad = 0.0; // Illegal initial value.
                    if (offset_dist_ < 0.0)
                    {
                        curv_rad = ((curv_rad1 < 0.0 && curv_rad2 < 0.0)) ?
                            std::max(curv_rad1, curv_rad2) : std::min(curv_rad1, curv_rad2);
                    }
                    else
                    {
                        curv_rad = ((curv_rad1 > 0.0 && curv_rad2 > 0.0)) ?
                            std::min(curv_rad1, curv_rad2) : std::max(curv_rad1, curv_rad2);
                    }                        
                    //std::cout << "curv_rad: " << curv_rad << std::endl;
                    radius_of_curv.push_back(curv_rad);
                }
                else
                {
                    //no_self_int.insert(no_self_int.end(), offset_pt.begin(), offset_pt.end());
                    no_self_int.insert(no_self_int.end(), local_sf_pt.begin(), local_sf_pt.end());
                    no_self_int_offset.insert(no_self_int_offset.end(), offset_pt.begin(), offset_pt.end());
                }
            }

        }

        MESSAGE("INFO: curv_rad_pos_min: " << curv_rad_pos_min << ", curv_rad_pos_max: " << curv_rad_pos_max <<
                ", curv_rad_neg_min: " << curv_rad_neg_min << ", curv_rad_neg_max: " << curv_rad_neg_max);

#ifndef NDEBUG
#if 1
        {
            // We write to file the two sets. Using green color for no self int, red color for self int.
            MESSAGE("Writing to file the self int and no self int points.");
            if (self_int.size() > 0)
            {
                std::ofstream fileout_debug("tmp/grid_self_int.g2");
                std::ofstream fileout_debug2("tmp/grid_offset_self_int.g2");
                PointCloud3D pt_cl(self_int.begin(), self_int.size()/3);
                vector<int> color_red(4, 0);
                color_red[0] = 255;
                color_red[3] = 255;
                ObjectHeader header(Class_PointCloud, 1, 0, color_red);
                header.write(fileout_debug);
                pt_cl.write(fileout_debug);

                PointCloud3D pt_cl2(self_int_offset.begin(), self_int_offset.size()/3);
                header.write(fileout_debug2);
                pt_cl2.write(fileout_debug2);
            }
            // Write to file the grid points on the original surfaces (parametrization defined by spline_sf_).
            if (no_self_int.size() > 0)
            {
                std::ofstream fileout_debug("tmp/grid_no_self_int.g2");
                std::ofstream fileout_debug2("tmp/grid_offset_no_self_int.g2");
                PointCloud3D pt_cl(no_self_int.begin(), no_self_int.size()/3);
                vector<int> color_green(4, 0);
                color_green[1] = 255;
                color_green[3] = 255;
                ObjectHeader header(Class_PointCloud, 1, 0, color_green);
                header.write(fileout_debug);
                pt_cl.write(fileout_debug);

                PointCloud3D pt_cl2(no_self_int_offset.begin(), no_self_int_offset.size()/3);
                header.write(fileout_debug2);
                pt_cl2.write(fileout_debug2);
                
            }
            // Write to file the corresponding offset grid points.
        }
#endif
#endif
    
    }


    //===========================================================================
    void EvalOffsetSurface::gridKinks(const HermiteGrid2D& grid,
                                      const vector<shared_ptr<SplineCurve> >& kink_cvs_2d,
                                      vector<int>& grid_kinks) const
    //===========================================================================
    {
        grid_kinks.clear();
        int num_kinks = 0;
        MESSAGE("Under construction!");

//      vector<double> self_int, no_self_int;

        vector<double> knots_u = grid.getKnots(true);
        vector<double> knots_v = grid.getKnots(false);
        vector<Point> data = grid.getData(); // Ordered row-wise: Pos, der_u, der_v, der_uv.

        const double epspar = epsgeo_; // @@sbr201705 We should map from epsgeo to epspar!
        
        /// Return the spatial dimension
        //int dim = grid.dim();
        const int MM = grid.size1();
        const int NN = grid.size2();
        MESSAGE("INFO: data.size(): " << data.size() << ", MM: " << MM << ", NN: " << NN);

        vector<double> kink_pts_3d;
        
        for (size_t kj = 0; kj < knots_v.size(); ++kj)
        {
            double vpar = knots_v[kj];
            for (size_t ki = 0; ki < knots_u.size(); ++ki)
            {
                double upar = knots_u[ki];

                Point epar(upar, vpar);
                double min_clo_dist = MAXDOUBLE;
                for (size_t kk = 0; kk < kink_cvs_2d.size(); ++kk)
                {
                    double clo_t, clo_dist;
                    Point clo_pt;
                    kink_cvs_2d[kk]->ParamCurve::closestPoint(epar, clo_t, clo_pt, clo_dist);
                    if (clo_dist < min_clo_dist)
                    {
                        min_clo_dist = clo_dist;
                    }
                }
                if (min_clo_dist < offset_dist_) // This should match the approach in the grid
                                                 // refinement. Probably slightly larger.
                {
                    grid_kinks.push_back(kj*MM+ki);
                    Point kink_pt = spline_sf_->ParamSurface::point(upar, vpar);
                    kink_pts_3d.insert(kink_pts_3d.end(), kink_pt.begin(), kink_pt.end());
                }
            }
        }

#ifndef NDEBUG
        {
            if (kink_pts_3d.size() > 0)
            {
                std::ofstream fileout_debug("tmp/grid_kink_pts.g2");
                PointCloud3D pt_cl(kink_pts_3d.begin(), kink_pts_3d.size()/3);
                vector<int> color_blue(4, 0);
                color_blue[2] = 255;
                color_blue[3] = 255;
                ObjectHeader header(Class_PointCloud, 1, 0, color_blue);
                header.write(fileout_debug);
                pt_cl.write(fileout_debug);
            }
        }
#endif
        
        MESSAGE("INFO: Number of samples too close to a kink: " << grid_kinks.size());
    }

    
    //===========================================================================
    vector<shared_ptr<SplineCurve> >
    EvalOffsetSurface::getProjKinkCurves(vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >& par_cvs,
                                         vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > >& sfs)
    //===========================================================================
    {
        vector<shared_ptr<SplineCurve> > kink_cvs_2d;
        MESSAGE("Under construction!");

        // We first locate all the kink curves in the surface set, i.e. all curves between adjacent
        // surfaces for which the normal along the edge does not coincide.
        const double ang_tol = 1e-03; // @@sbr201705 This should be global!

        vector<shared_ptr<ParamCurve> > kink_cvs_3d = get3DKinkCurves(par_cvs, sfs);

#ifndef NDEBUG
        {
            std::ofstream fileout_debug("tmp/kinks_3d.g2");
            for (size_t ki = 0; ki < kink_cvs_3d.size(); ++ki)
            {
                if (kink_cvs_3d[ki]->instanceType() == Class_CurveOnSurface)
                {
                    shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(kink_cvs_3d[ki]);
                    if (cv_on_sf->geometryCurve() != NULL)
                    {
                        cv_on_sf->geometryCurve()->writeStandardHeader(fileout_debug);
                        cv_on_sf->geometryCurve()->write(fileout_debug);
                    }
                }
                else
                {
                    kink_cvs_3d[ki]->writeStandardHeader(fileout_debug);
                    kink_cvs_3d[ki]->write(fileout_debug);
                }
            }
        }
#endif

        // We project all the 3d curves onto the approximating SplineSurface.
        for (size_t ki = 0; ki < kink_cvs_3d.size(); ++ki)
        {
#if 0
            try
            {
                shared_ptr<Point> start_par_pt;
                shared_ptr<Point> end_par_pt;
                shared_ptr<ParamSurface> param_sf = base_sf_->surface();
                // @@sbr201705 The epsgeo_ is used for testing distance from curve to surface, what we
                // need to test is the approximation deviation from the projected point.
                MESSAGE("Fix problem with projection of curves far away from the surface!");
                shared_ptr<SplineCurve> par_cv(CurveCreators::projectSpaceCurve(kink_cvs_3d[ki],
                                                                                param_sf,
                                                                                start_par_pt,
                                                                                end_par_pt,
                                                                                epsgeo_));
                if (par_cv.get() != NULL)
                {
                    kink_cvs_2d.push_back(par_cv);
                }
                else
                {
                    MESSAGE("Failed projecting curve!");
                }
            }
            catch (...)
            {
                MESSAGE("Failed projecting curve!");
            }
#else
            shared_ptr<ParamSurface> param_sf = base_sf_->surface();
            shared_ptr<SplineCurve> proj_cv, par_cv;
            CurveCreators::projectCurve(kink_cvs_3d[ki], param_sf, epsgeo_,
                                        proj_cv, par_cv);
            if (par_cv.get() != NULL)
            {
                kink_cvs_2d.push_back(par_cv);
            }
            else
            {
                MESSAGE("Failed projecting curve!");
            }            
#endif

//             ProjectCurve(shared_ptr<Go::ParamCurve>& space_crv,
// 		 shared_ptr<Go::ParamSurface>& surf,
// 		 shared_ptr<Go::Point>& start_par_pt,
// 		 shared_ptr<Go::Point>& end_par_pt,
// 		 double epsgeo1,
// // 		 double epsgeo2,
// 		 const RectDomain* domain_of_interest = NULL);

        }
        
        MESSAGE("INFO: kink_cvs_2d.size(): " << kink_cvs_2d.size());

#ifndef NDEBUG
        { // We lift the projected curves onto to the approximating spline surface.
            std::ofstream fileout_debug("tmp/lifted_kink_cvs.g2");
            for (size_t ki = 0; ki < kink_cvs_2d.size(); ++ki)
            {
                shared_ptr<ParamCurve> par_cv = kink_cvs_2d[ki];
                shared_ptr<ParamSurface> par_sf = spline_sf_;
                shared_ptr<SplineCurve> lifted_cv(CurveCreators::liftParameterCurve(par_cv,
                                                                                    par_sf,
                                                                                    epsgeo_));
                lifted_cv->writeStandardHeader(fileout_debug);
                lifted_cv->write(fileout_debug);
            }
        }
#endif

        return kink_cvs_2d;
    }
    
                
    //===========================================================================
    ParamSurface* EvalOffsetSurface::findLocalSurface(double u, double v,
                                                      double& local_u, double& local_v) const
    //===========================================================================
    {
        ParamSurface* local_sf = NULL;

        if (dynamic_pointer_cast<ftChartSurface>(base_sf_).get() != NULL) {
            // @@sbr201703 Easy to project the tangents. But what about the twist vector? 
            shared_ptr<ftChartSurface> chart_sf = dynamic_pointer_cast<ftChartSurface>(base_sf_);
            // We locate the corresponding sub-face and parameter value.
            shared_ptr<ftFaceBase> face;
            local_u = u;
            local_v = v;
            Point proj_pt = chart_sf->point(local_u, local_v, face);
            // std::cout << "u: " << u << ", v: " << v <<
            //     ", local_u: " << local_u << ", local_v: " << local_v << std::endl;

            shared_ptr<ParamSurface> param_sf = face->surface();
            if (param_sf.get() != 0)
            {
                local_sf = param_sf.get();
            }
        } else if (dynamic_pointer_cast<ftSurface>(base_sf_).get() != NULL) {
            // For a ftSurface the local par is the same as the global par.
            local_sf = base_sf_->surface().get();
            local_u = u;
            local_v = v;
        } else {
            MESSAGE("Unexpected surface type!");
        }

        if (local_sf == NULL)
        {
            MESSAGE("findLocalSurface(): Failed!");
        }
        
        return local_sf;
    }


    //===========================================================================
    vector<shared_ptr<ParamCurve> >
    EvalOffsetSurface::get3DKinkCurves(std::vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >& kink_cvs_2d,
                                       std::vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > >& under_sfs)
    //===========================================================================
    {
        vector<shared_ptr<ParamCurve> > kink_cvs_3d;
        // vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > > kink_cvs_2d;
        // vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > > under_sfs; // Corresponding toe the kink_cvs_2d.
        
        std::set<ftEdgeBase*> edge_set;
        if (dynamic_pointer_cast<ftChartSurface>(base_sf_).get() != NULL)
        {
            // @@sbr201703 Easy to project the tangents. But what about the twist vector? 
            shared_ptr<ftChartSurface> chart_sf = dynamic_pointer_cast<ftChartSurface>(base_sf_);
            vector<shared_ptr<FaceConnectivity<ftEdgeBase> > > inner_edge_cont = chart_sf->getInnerEdgeCont();

            for (size_t ki = 0; ki < inner_edge_cont.size(); ++ki)
            {
                for (size_t kj = 0; kj < inner_edge_cont[ki]->status_.size(); ++kj)
                {
                    // We include cases with a gap as these requires us to perform smoothing in that area.
                    if (inner_edge_cont[ki]->status_[kj] > 0)
                    {
                        MESSAGE("Edge cont not c0 (status > 0), status = " << inner_edge_cont[ki]->status_[kj]);
                        ftEdgeBase* e1 = inner_edge_cont[ki]->e1_;
                        ftEdgeBase* e2 = inner_edge_cont[ki]->e2_;
                        ftEdge* geom_edge1 = e1->geomEdge();
                        ftEdge* geom_edge2 = e2->geomEdge();
                        shared_ptr<ParamCurve> geom_cv1 = geom_edge1->geomCurve();
                        shared_ptr<ParamCurve> geom_cv2 = geom_edge2->geomCurve();
                        // std::cout << "status_.size(): " << inner_edge_cont[ki]->status_.size() <<
                        //     ", parameters_.size(): " << inner_edge_cont[ki]->parameters_.size() << std::endl;
                        double tmin1 = inner_edge_cont[ki]->parameters_[kj].first;
                        double tmax1 = inner_edge_cont[ki]->parameters_[kj+1].first;
                        if (tmax1 < tmin1)
                        {
                            std::swap(tmin1, tmax1);
                        }
                        double tmin2 = inner_edge_cont[ki]->parameters_[kj].second;
                        double tmax2 = inner_edge_cont[ki]->parameters_[kj+1].second;
                        if (tmax2 < tmin2)
                        {
                            std::swap(tmin2, tmax2);
                        }
                        shared_ptr<ParamCurve> sub_cv1(geom_cv1->subCurve(tmin1, tmax1));
                        shared_ptr<ParamCurve> sub_cv2(geom_cv2->subCurve(tmin2, tmax2));
                        // We only need the space curve from one of the edges.
                        // @@sbr201705 Possibly use the parameter curve for the actual surface.
                        // If the twin edge was already added we skip this edge.
                        if (edge_set.find(e1) == edge_set.end())
                        { // Either none of the edges is included or both.
                            kink_cvs_3d.push_back(sub_cv1);
                            edge_set.insert(e1);
                            edge_set.insert(e2);
                            if ((sub_cv1->instanceType() == Class_CurveOnSurface) &&
                                (sub_cv2->instanceType() == Class_CurveOnSurface))
                            {
                                shared_ptr<CurveOnSurface> cv_on_sf1 = dynamic_pointer_cast<CurveOnSurface>(sub_cv1);
                                shared_ptr<CurveOnSurface> cv_on_sf2 = dynamic_pointer_cast<CurveOnSurface>(sub_cv2);
                                if ((cv_on_sf1->parameterCurve() != NULL) && (cv_on_sf2->parameterCurve() != NULL))
                                {
                                    kink_cvs_2d.push_back(std::make_pair(cv_on_sf1->parameterCurve(),
                                                                         cv_on_sf2->parameterCurve()));
                                    under_sfs.push_back(std::make_pair(cv_on_sf1->underlyingSurface(),
                                                                       cv_on_sf2->underlyingSurface()));
                                }
                            }
                        }
                        // kink_cvs_3d.push_back(sub_cv2);
                    }
                }
            }
        }

        return kink_cvs_3d;
    }


} // namespace Go
