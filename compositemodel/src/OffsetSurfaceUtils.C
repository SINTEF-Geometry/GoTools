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

#include "GoTools/compositemodel/OffsetSurfaceUtils.h"

#include "GoTools/compositemodel/ftChartSurface.h"
#include "GoTools/compositemodel/EvalOffsetSurface.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/HermiteApprEvalSurf.h"


namespace Go
{

namespace OffsetSurfaceUtils
{
    
int offsetSurfaceSet(const std::vector<shared_ptr<ParamSurface> >& param_sfs,
                     double offset_dist, double epsgeo,
                     shared_ptr<SplineSurface>& offset_sf)
{
    int status = 0; // Status 0 => success.
    
    int id = 0;
    shared_ptr<ftFaceBase> base_sf;
    if (param_sfs.size() == 1)
    {
        MESSAGE("DEBUGGING: Temporarily using ftChartSurface for cases with 1 surface also!");
    }
    if (0)//param_sfs.size() == 1)
    {
        base_sf = shared_ptr<ftFaceBase>(new ftSurface(param_sfs[0], id));
    }
    else
    {
        // We need to use a ftChartSurface.
        const double gap = epsgeo;
        const double kink = 1.0e-02;
        tpTolerances top_eps(gap, 10.0*gap, kink, 10.0*kink);
        const double approx = epsgeo*100.0;
        const double curvature_tol = 0.01; // @@sbr201702 Not used yet.
        int m = 0;
        int n = 0;
        shared_ptr<ftChartSurface> chart_sf(new ftChartSurface(param_sfs, top_eps,
                                                               approx, curvature_tol,
                                                               m, n));

        // Necessary initialization
        vector<shared_ptr<ftEdgeBase> > bd_edges;
        bd_edges = chart_sf->createInitialEdges(10.0*gap);

        // Create the chart surface.
        // This surface defines the domain of our offset surface.
        double max_error = 2.0;
        double mean_error = 1.0;
        chart_sf->createSurf(max_error, mean_error);
        // The surface should reflect the geometry of the input surfaces, allowing some deviation.
        std::cout << "Max error: " << max_error << std::endl;
        std::cout << "Mean error: " << mean_error << std::endl;

        // We check the spline surface defining the parameter domain of the chart surface.
        shared_ptr<ParamSurface> out_surf = chart_sf->surface();
        shared_ptr<SplineSurface> spline_out_surf =
            dynamic_pointer_cast<SplineSurface, ParamSurface>(out_surf);
        // status = 2;
        // return status;

        std::ofstream fileout("tmp/chart_sf.g2");
        spline_out_surf->writeStandardHeader(fileout);
        spline_out_surf->write(fileout);
        
        base_sf = chart_sf;
    }
        
    EvalOffsetSurface eval_offset_sf(base_sf, offset_dist, epsgeo);

    // Creating the initial grid.
    // Only the end parameters are set initially.
    HermiteApprEvalSurf appr_eval_sf(&eval_offset_sf, epsgeo, epsgeo);
    try
    { 
        // The method refines until within the required tolerance (or aborts if a knot interval gets too small).
        appr_eval_sf.refineApproximation();
        bool method_failed;
        // Creating the surface from the Bezier patches.
        offset_sf = appr_eval_sf.getSurface(method_failed);
        if (method_failed)
        {
            status = 1;
        }
        if (offset_sf.get() == 0)
        {
            status = 1;
        }
    }
    catch (...)
    {
        status = 1;
    }

    return status;    
}

} // namespace OffsetSurfaceUtils
    
} // namespace Go
