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

#include "GoTools/creators/CreatorsOffsetUtils.h"

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

        const int kder = 2; // To compute the twist.
        Point epar(u, v);
        int ind_u=0;             /* Pointer into knot vector                       */
        int ind_v=0;             /* Pointer into knot vector                       */
        vector<Point> offset_pt(kder*(kder+1) + 1); // Derivs & normal.
        vector<Point> base_pt(kder*(kder+1) + 1); // Derivs & normal.
        int kstat = 0;
        if (spline_sf_ == 0) {
            THROW("Missing support for parametric surface as a SplineSurface!");
        }

        OffsetUtils::blend_s1421(spline_sf_, offset_dist_, kder, epar, ind_u, ind_v,
                                 offset_pt, base_pt, &kstat);
        der[0] = offset_pt[0];
        der[1] = offset_pt[1];
        der[2] = offset_pt[2];
        der[3] = offset_pt[4];
        
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

        const bool use_geom_check = true;//false;//true;
        if ((!appr_ok) && use_geom_check)
        {
            // We also check using closest point.
            double seed[2];
            seed[0] = par_u;
            seed[1] = par_v;
            double clo_u, clo_v;
            double clo_dist = -1.0;
            Point clo_pt;
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
