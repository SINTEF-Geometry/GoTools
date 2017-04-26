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
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/CurvatureAnalysis.h"
#include "GoTools/creators/SmoothSurf.h"


namespace Go
{

void boundaryCurvatureRadius(ftFaceBase& face,
                             double& curv_radius_min,
                             double& curv_radius_max);

shared_ptr<SplineSurface> getSmoothOffsetSurface(shared_ptr<SplineSurface> offset_sf,
                                                 const vector<int>& grid_self_int);//const HermiteGrid2D& grid);


namespace OffsetSurfaceUtils
{
    
OffsetSurfaceStatus offsetSurfaceSet(const std::vector<shared_ptr<ParamSurface> >& param_sfs,
                                     double offset_dist, double epsgeo,
                                     shared_ptr<SplineSurface>& offset_sf)
{
    OffsetSurfaceStatus status = OFFSET_OK; // Status 0 => success.
    
    int id = 0;
    shared_ptr<ftFaceBase> base_sf;
    if (param_sfs.size() == 1)
    {
        base_sf = shared_ptr<ftFaceBase>(new ftSurface(param_sfs[0], id));
    }
    else
    {
        if (param_sfs.size() == 1)
        {
            MESSAGE("DEBUGGING: Temporarily using ftChartSurface for cases with 1 surface also!");
        }
        
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
        ftMessage ft_status = chart_sf->createSurf(max_error, mean_error);
        if (ft_status.getMessage() == FT_NON_4_SIDED_SURF)
        {
            status = NOT_FOUR_CORNERS;
            return status;
        }
        // The surface should reflect the geometry of the input surfaces, allowing some deviation.
        std::cout << "Max error: " << max_error << std::endl;
        std::cout << "Mean error: " << mean_error << std::endl;

        // We check the spline surface defining the parameter domain of the chart surface.
        shared_ptr<ParamSurface> out_surf = chart_sf->surface();
        shared_ptr<SplineSurface> spline_out_surf =
            dynamic_pointer_cast<SplineSurface, ParamSurface>(out_surf);

        if (spline_out_surf.get() == NULL)
        {
            status = OFFSET_FAILED;
            return status;
        }

        std::ofstream fileout("tmp/chart_sf.g2");
        spline_out_surf->writeStandardHeader(fileout);
        spline_out_surf->write(fileout);
        
        base_sf = chart_sf;        
    }

    // We analyze the boundary of the surface set.  Currently we do not support self intersections along
    // the boundary.  We only check the curvature along the edge. If curvature is too high in other
    // directions it should be handled by smoothing in the inner part of the offset surface (keeping the
    // coefs along the offset surface boundary edges fixed).
    std::vector<shared_ptr<ftEdgeBase> > edges = base_sf->createInitialEdges();
    double curv_radius_min = MAXDOUBLE;
    double curv_radius_max = -MAXDOUBLE;
    boundaryCurvatureRadius(*base_sf,
                            curv_radius_min,
                            curv_radius_max);
    std::cout << "INFO: curv_radius_min: " << curv_radius_min << ", curv_radius_max: " << curv_radius_max << std::endl;

    if (curv_radius_min < offset_dist)
    {
        MESSAGE("The radius of curvature along the boundary is less than the offset dist, not supported.");
#if 1
        MESSAGE("Turned off early exit!");
#else        
        status = SELF_INTERSECTING_BOUNDARY;
        return status;
#endif
    }
    
    EvalOffsetSurface eval_offset_sf(base_sf, offset_dist, epsgeo);

    // Creating the initial grid.
    // Only the end parameters are set initially.
    HermiteApprEvalSurf appr_eval_sf(&eval_offset_sf, epsgeo, epsgeo);
    try
    {
        // @@sbr201704 We calculate the radius of curvature (min) in inner grid points. We use this to
        // define areas which need smoothing, possibly in a post processing step.
        MESSAGE("Missing analysis of radius of curvature in the grid!");
        // The method refines until within the required tolerance (or aborts if a knot interval gets too small).
        appr_eval_sf.refineApproximation();
        bool method_failed;

        // Creating the surface from the Bezier patches.
        offset_sf = appr_eval_sf.getSurface(method_failed);
        if (method_failed)
        {
            status = OFFSET_FAILED;
        }
        if (offset_sf.get() == 0)
        {
            status = OFFSET_FAILED;
        }

        const HermiteGrid2D& grid = appr_eval_sf.getGrid();
        std::cout << "grid.size1(): " << grid.size1() << ", grid2.size(): " << grid.size2() << std::endl;
        
        // We check if the grid contains any self intersections. If any such self intersections are along
        // the boundary we return with an error message.
        vector<int> grid_self_int = eval_offset_sf.gridSelfIntersections(grid);
        size_t num_self_int = grid_self_int.size();
        if (num_self_int > 0)
        {
            MESSAGE("Self intersections found! We perform smoothing on the area affected!");
            std::cout << "Offset grid num_self_int: " << num_self_int << std::endl;
            // @@sbr201704 We should use smoothing with constraints.
#if 1
            {
                std::ofstream fileout_debug("tmp/selfint_offset_sf.g2");
                offset_sf->writeStandardHeader(fileout_debug);
                offset_sf->write(fileout_debug);
            }
#endif
            shared_ptr<SplineSurface> smooth_offset_sf = getSmoothOffsetSurface(offset_sf, grid_self_int);
#if 1
            {
                std::ofstream fileout_debug("tmp/smooth_offset_sf.g2");
                smooth_offset_sf->writeStandardHeader(fileout_debug);
                smooth_offset_sf->write(fileout_debug);
            }
#endif
        }
    }
    catch (...)
    {
        status = OFFSET_FAILED;
    }

    return status;
}

} // namespace OffsetSurfaceUtils

void boundaryCurvatureRadius(ftFaceBase& face,
                             double& curv_radius_min,
                             double& curv_radius_max)
{
    curv_radius_min = MAXDOUBLE;
    curv_radius_max = -MAXDOUBLE;

    std::vector<shared_ptr<ftEdgeBase> > start_edges = face.startEdges();

    std::cout << "start_edges.size(): " << start_edges.size() << std::endl;

    if (start_edges.size() != 1)
    {
        MESSAGE("Not supporting surfaces with inner loops.");
        return;
    }
    
    ftEdgeBase* curr_edge = start_edges[0].get();
    ftEdgeBase* first_edge = curr_edge;
    while (true)
    {
        ftEdge* ft_edge = curr_edge->geomEdge();
        shared_ptr<ParamCurve> geom_cv = ft_edge->geomCurve();
        shared_ptr<SplineCurve> spline_cv;
        if (geom_cv->instanceType() == Class_SplineCurve)
        {
            spline_cv = dynamic_pointer_cast<SplineCurve>(geom_cv);
        }
        else if (geom_cv->instanceType() == Class_CurveOnSurface)
        {
            shared_ptr<CurveOnSurface> cv_on_sf = dynamic_pointer_cast<CurveOnSurface>(geom_cv);
            shared_ptr<ParamCurve> space_cv = cv_on_sf->spaceCurve();
            if (space_cv.get() != NULL)
            {
                if (space_cv->instanceType() == Class_SplineCurve)
                {
                    spline_cv = dynamic_pointer_cast<SplineCurve>(space_cv);
                }
            }
        }

        if (spline_cv.get() == NULL)
        {
            MESSAGE("Only spline curves currently supported! Not type: " << geom_cv->instanceType());
            continue;
        }

        int num_samples = spline_cv->numCoefs();

        const BsplineBasis& basis = spline_cv->basis();
        for (int kj = 0; kj < num_samples; ++kj)
        {
            double grev_par = basis.grevilleParameter(kj);
            vector<Point> cv_pt = spline_cv->ParamCurve::point(grev_par, 1);
            double clo_u, clo_v, clo_dist;
            Point clo_pt;
            const double epsgeo = 1e-06;
            shared_ptr<ParamSurface> param_sf = face.surface();
            param_sf->closestBoundaryPoint(cv_pt[0],
                                           clo_u, clo_v, clo_pt, clo_dist, epsgeo);
            //std::cout << "clo_u: " << clo_u << ", clo_v: " << clo_v << ", clo_dist: " << clo_dist << std::endl;
            // We need the surface normal to calculate the curvature (to define the direction).
            // Or the first and second derivatives in the surface along the edge curve.

            double k1, k2;
            Point d1, d2;
            CurvatureAnalysis::principalCurvatures(*param_sf, clo_u, clo_v, k1, d1, k2, d2);

            double curv_rad1 = 1.0/k1;
            double curv_rad2 = 1.0/k2;

            vector<Point> sf_pt = param_sf->point(clo_u, clo_v, 1);
            
            // Based on the direction of the curve tangent we pick one of the partial derivatives.
            double curv_rad;
            // We compute the smallest angle between vectors regardless of orientation (i.e. +/-).
            double ang1 = cv_pt[1].angle_smallest(sf_pt[1]);
            double ang2 = cv_pt[1].angle_smallest(sf_pt[2]);
            const double ang_tol = 1.0e-03;
            if ((ang1 < ang_tol) && (ang2 < ang_tol))
            {
                // If the partial derivatives are parallel we expect problems. For now we use the minimum
                // radius of curvature.
                MESSAGE("The partial derivatives are parallel, did not expect this!");
                curv_rad = std::min(curv_rad1, curv_rad2);
            }
            else
            {
                curv_rad = (ang1 < ang2) ? curv_rad1 : curv_rad2;
            }
#if 1
//            std::cout << "curv_rad: " << curv_rad << std::endl;

            if (curv_rad > curv_radius_max)
            {
                curv_radius_max = curv_rad;
            }
            if ((curv_rad > 0.0) && (curv_rad <  curv_radius_min))
            {
                curv_radius_min = curv_rad;
            }
#else       
            // The radius of curvature along the boundary should not be smaller than the offset dist.
            // With the exception for negative values.
            if (curv_rad1 > curv_radius_max)
            {
                curv_radius_max = curv_rad1;
            }
            if (curv_rad2 > curv_radius_max)
            {
                curv_radius_max = curv_rad2;
            }

            if ((curv_rad1 > 0.0) && (curv_rad1 <  curv_radius_min))
            {
                curv_radius_min = curv_rad1;
            }
            if ((curv_rad2 > 0.0) && (curv_rad2 <  curv_radius_min))
            {
                curv_radius_min = curv_rad2;
            }
            
            std::cout << "curv_rad1: " << curv_rad1 << ", curv_rad2: " << curv_rad2 << std::endl;
#endif
        }

        curr_edge = curr_edge->next();
        if (curr_edge == first_edge)
        {
            break;
        }
    }
    

}


//===========================================================================
shared_ptr<SplineSurface> getSmoothOffsetSurface(shared_ptr<SplineSurface> offset_sf,
                                                 const vector<int>& grid_self_int)
//===========================================================================
{
    shared_ptr<SplineSurface> new_offset_sf;
    
    MESSAGE("Under construction!");

    if (grid_self_int.size() == 0)
    {
        MESSAGE("No self intersections, method should not have been called!");
        return new_offset_sf;
    }
    
    int num_coefs_u = offset_sf->numCoefs_u();
    int num_coefs_v = offset_sf->numCoefs_v();
    vector<int> coef_known(num_coefs_u*num_coefs_v, 1);
    // The offset_sf & grid_self_int are referring to the same space, i.e. the grid is the basis for the
    // Bezier patches constituting the offset_sf. We thus can recreate the grid dimensionality.
    int grid_mm = (num_coefs_u)/2;
    int grid_nn = (num_coefs_v)/2;
    std::cout << "grid_mm: " << grid_mm << ", grid_nn: " << grid_nn << std::endl;
    for (size_t ki = 0; ki < grid_self_int.size(); ++ki)
    {
        // We find the positions in the grid.
        int ki_mm = grid_self_int[ki]%grid_mm;
        int ki_nn = grid_self_int[ki]/grid_mm;
        std::cout << "ki_mm: " << ki_mm << ", ki_nn: " << ki_nn << std::endl;
//        coef_known[
        coef_known[(2*ki_nn)*num_coefs_u + (2*ki_mm)] = 0;
        coef_known[(2*ki_nn)*num_coefs_u + (2*ki_mm+1)] = 0;
        coef_known[(2*ki_nn+1)*num_coefs_u + (2*ki_mm)] = 0;
        coef_known[(2*ki_nn+1)*num_coefs_u + (2*ki_mm+1)] = 0;
    }
    
    shared_ptr<SplineSurface> smooth_offset_sf; // Without self intersections.

    SmoothSurf smooth_sf;
    vector<int> seem(2, 0);
    smooth_sf.attach(offset_sf, &seem[0], &coef_known[0]);
        
    double smooth_weight = 1.0e-3;
    double wgt1 = 0.0;
    double wgt3 = 0.0;//(min_der >= 3) ? 0.5*smoothweight_ : 0.0;
    double wgt2 = (1.0 - wgt3)*smooth_weight;
    wgt3 *= smooth_weight;
    double weight_sum = wgt1 + wgt2 + wgt3;
    if (weight_sum > 1.0) { // Making sure the sum of the weights is at most 1.0.
        wgt1 /= weight_sum;
        wgt2 /= weight_sum;
        wgt3 /= weight_sum;
    }
    smooth_sf.setOptimize(wgt1, wgt2, wgt3);

    double approxweight = 1.0 - (wgt1 + wgt2 + wgt3); // In the range [0.0, 1.0].
    double wgt_orig = 0.1*approxweight;
    wgt_orig = 0.0;
    if (wgt_orig > 0.0)
    {
        smooth_sf.approxOrig(wgt_orig);
    }

    smooth_sf.equationSolve(new_offset_sf);
        
    return new_offset_sf;
}


} // namespace Go

