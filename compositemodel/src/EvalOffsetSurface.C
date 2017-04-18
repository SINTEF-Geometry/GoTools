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

#include <vector>
#include <assert.h>

using std::vector;

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
            spline_sf_ = param_sf->asSplineSurface();
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
        const SplineSurface* spline_sf = spline_sf_;
        const SplineSurface* spline_sf_global = spline_sf_;
        const SplineSurface* spline_sf_local = spline_sf_; // For the case with 1 surface in the set.
        // For the case with 1 surface only the spline_sf_ should is the same as the original surface.
        // For the surface set we must project the point onto the corresponding input surface.
        bool surface_set = false;
        Point epar_global(u, v);
        Point epar_local(2); // For storing value of local spline sf.
        double u2, v2;
        if (dynamic_pointer_cast<ftChartSurface>(base_sf_).get() != NULL) {
            // @@sbr201703 Easy to project the tangents. But what about the twist vector? 
            shared_ptr<ftChartSurface> chart_sf = dynamic_pointer_cast<ftChartSurface>(base_sf_);
            surface_set = true;
            // We locate the corresponding sub-face and parameter value.
            u2 = u;
            v2 = v;
            shared_ptr<ftFaceBase> face;
            Point proj_pt = chart_sf->point(u2, v2, face);
            //std::cout << "u: " << u << ", v: " << v << ", u2: " << u2 << ", v2: " << v2 << std::endl;

            shared_ptr<ParamSurface> param_sf = face->surface();
            if (param_sf.get() != 0)
            {
                spline_sf = param_sf->asSplineSurface();
                spline_sf_local = spline_sf;
            }
            epar_local[0] = u2;
            epar_local[1] = v2;
        } else if (dynamic_pointer_cast<ftSurface>(base_sf_).get() == NULL) {
            MESSAGE("Unexpected surface type!");
            return;
        }

        const int kder = 2; // To compute the twist.
        int ind_u=0;             /* Pointer into knot vector                       */
        int ind_v=0;             /* Pointer into knot vector                       */
        int kstat = 0;
        if (spline_sf_ == 0) {
            THROW("Missing support for parametric surface as a SplineSurface!");
        }
//        MESSAGE("Using spline_sf_ when computing derivs, that is wrong! Must use base_sf_!");

        // We must blend the directions of the base_sf_ to match the directions of the spline_sf_.
        vector<Point> offset_pt_local(kder*(kder+1) + 1); // Derivs & normal in the exact surface.
        vector<Point> base_pt(kder*(kder+1) + 1); // Derivs & normal.
        OffsetUtils::blend_s1421(spline_sf, offset_dist_, kder, epar_local, ind_u, ind_v,
                                 offset_pt_local, base_pt, &kstat);

        vector<Point> offset_pt_global(kder*(kder+1) + 1); // Derivs & normal in the approximated surface.
        vector<Point> base_pt_global(kder*(kder+1) + 1); // Derivs & normal.
        OffsetUtils::blend_s1421(spline_sf_global, offset_dist_, kder, epar_global, ind_u, ind_v,
                                 offset_pt_global, base_pt_global, &kstat);

        if (surface_set) {
//            MESSAGE("MISSING: We must express the local derivs using the derivs in spline_sf_!");
            vector<Point> global_pt(3), local_pt(3);
            spline_sf_global->point(global_pt, u, v, 1);
            spline_sf_local->point(local_pt, u2, v2, 1);
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
#if 1
            // Setting the twist vector to 0.0. The vector is tricky to calculate with the parameter
            // domain defined as the closest point from the approximating global surface.
//            der[3] = Point(0.0, 0.0, 0.0);
            der[3] = offset_pt_global[4]; // Generally the twist in global_sf is better than the zero twist.
            // der[3] = offset_pt[4];
#else
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

        const bool use_geom_check = false;//true;
        if ((!appr_ok) && use_geom_check)
        {
            MESSAGE("Missing closest point for the surface set!");
            // We also check using closest point.
            double seed[2];
            seed[0] = par_u;
            seed[1] = par_v;
            double clo_u, clo_v;
            double clo_dist = -1.0;
            Point clo_pt;
            // @@sbr201703 The spline_sf_ is used for defining the domain only, it is not relevant to use
            // the position for offset evaluations.
            MESSAGE("We must use the surface set, not the approximating spline surface!");
            spline_sf_->closestPoint(approxpos, clo_u, clo_v, clo_pt, clo_dist,
                                     tol2*1e-04, NULL, seed);

            double offset_dist = approxpos.dist(clo_pt);
            double clo_pt_error = std::fabs(offset_dist - offset_dist_);
            if (clo_pt_error < dist)
            {
                //std::cout << "Closest point approach was a success!" << std::endl;
                if (clo_pt_error < tol1)
                {
                    //  std::cout << "We are now inside the tolerance!" << std::endl;
                    appr_ok = true;
                }
            }
        }
        
        return appr_ok;
    }


} // namespace Go
