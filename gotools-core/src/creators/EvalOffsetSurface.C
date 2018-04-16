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




#include "GoTools/creators/EvalOffsetSurface.h"

#include "GoTools/creators/CreatorsOffsetUtils.h"

#if 0
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/CurvatureAnalysis.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/creators/CurveCreators.h"
#endif

#include <vector>
#include <assert.h>

using std::vector;
using std::pair;

namespace Go
{

    //===========================================================================
    EvalOffsetSurface::EvalOffsetSurface(shared_ptr<ParamSurface> sf,
                                         double offset_dist, double epsgeo)
    //===========================================================================
        : sf_(sf), offset_dist_(offset_dist), epsgeo_(epsgeo)
        //, rect_dom_(sf->containingDomain())
    {
        spline_sf_ = shared_ptr<SplineSurface>(sf->asSplineSurface());
    }


    //===========================================================================
    EvalOffsetSurface::~EvalOffsetSurface()
    //===========================================================================
    {
    }


    //===========================================================================
    Point EvalOffsetSurface::eval(double u, double v) const
    //===========================================================================
    {
        Point pt = sf_->point(u, v);

        Point normal;
        sf_->normal(normal, u, v);
        normal.normalize();

        Point offset_pt = pt + normal*offset_dist_;
        
        return offset_pt;
    }


    //===========================================================================
    void EvalOffsetSurface::eval(double u, double v, int n, Point der[]) const
    //===========================================================================
    {
        if (n == 0)
        {
            der[0] = eval(u, v);
            return;
        }

        if (spline_sf_.get() == nullptr)
        {
            THROW("The surface was not represented by a SplineSurface!");
        }

        const int kder = 2; // To compute the twist.
        int ind_u=0;             /* Pointer into knot vector                       */
        int ind_v=0;             /* Pointer into knot vector                       */
        int kstat = 0;
        if (sf_.get() == 0) {
            THROW("Missing support for parametric surface as a SplineSurface!");
        }

        // We must blend the directions of the sf_ to match the directions of the sf_.
        Point epar(2);
        epar[0] = u;
        epar[1] = v;
        vector<Point> offset_pt(((kder+1)*(kder+2)/2) + 1); // Derivs & normal in the exact surface.
        vector<Point> base_pt(((kder+1)*(kder+2)/2) + 1); // Derivs & normal.
        OffsetUtils::blend_s1421(spline_sf_.get(), offset_dist_, kder, epar, ind_u, ind_v,
                                 offset_pt, base_pt, &kstat);

        der[0] = offset_pt[0];
        der[1] = offset_pt[1];
        der[2] = offset_pt[2];
        der[3] = offset_pt[4]; // The twist vector.
        
        return;
    }


    //===========================================================================
    double EvalOffsetSurface::start_u() const
    //===========================================================================
    {
        RectDomain rect_dom = sf_->containingDomain();
        double start_u = rect_dom.umin();
        
        return start_u;
    }



    //===========================================================================
    double EvalOffsetSurface::start_v() const
    //===========================================================================
    {
        RectDomain rect_dom = sf_->containingDomain();
        double start_v = rect_dom.vmin();
        
        return start_v;
    }


    //===========================================================================
    double EvalOffsetSurface::end_u() const
    //===========================================================================
    {
        RectDomain rect_dom = sf_->containingDomain();
        double end_u = rect_dom.umax();
        
        return end_u;
    }


    //===========================================================================
    double EvalOffsetSurface::end_v() const
    //===========================================================================
    {
        RectDomain rect_dom = sf_->containingDomain();
        double end_v = rect_dom.vmax();
        
        return end_v;
    }


    //===========================================================================
    int EvalOffsetSurface::dim() const
    //===========================================================================
    {
        int dim = sf_->dimension();
        
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

        return appr_ok;
    }


} // namespace Go
