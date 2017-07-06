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
#include "GoTools/geometry/ObjectHeader.h"

#include <algorithm>

namespace Go
{

void boundaryCurvatureRadius(ftFaceBase& face,
                             double& curv_radius_min,
                             double& curv_radius_max);

shared_ptr<SplineSurface> getSmoothOffsetSurface(shared_ptr<SplineSurface> offset_sf,
                                                 const vector<int>& coefs_released, // Index of the released coefs.
                                                 const vector<double>& kink_approx_pts,
                                                 const vector<double>& kink_approx_params);

void getIsoSelfIntersections(const HermiteGrid2D& grid, const vector<int>& grid_self_int,
                             vector<double>& iso_self_int_u, vector<double>& iso_self_int_v);

#if 1
void updateGridSelfInt(const HermiteGrid2D& grid,
                       const vector<double>& iso_self_int_u, const vector<double>& iso_self_int_v,
                       vector<int>& grid_self_int, vector<double>& radius_of_curv);
#endif

    
void fetchKinkPoints(const vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >& kink_par_cvs,
                     const vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > >& kink_sfs,
                     double offset, double epsgeo, shared_ptr<SplineSurface> spline_sf,
                     vector<double>& kink_pts, vector<double>& kink_params);

// Return index of nodes which are within smoothing_dist from a self-intersecting node.
vector<int> getGridSmoothingNodes(const HermiteGrid2D& grid,
                                  const vector<int>& grid_self_int,
                                  double smoothing_dist);

vector<int> getCoefsReleased(shared_ptr<SplineSurface> offset_sf,
                             const vector<int>& grid_self_int,
                             const vector<double>& radius_of_curv);//,
                             // const vector<int>& grid_kinks,
                             // const vector<double>& kink_release_dist);

void removeNonIsoCurves(vector<shared_ptr<SplineCurve> >& kink_cvs_2d,
                        vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >& par_cvs,
                        vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > >& sfs);
    
namespace OffsetSurfaceUtils
{
    
OffsetSurfaceStatus offsetSurfaceSet(const std::vector<shared_ptr<ParamSurface> >& param_sfs,
                                     double offset_dist, double epsgeo,
                                     shared_ptr<SplineSurface>& offset_sf)
{
    OffsetSurfaceStatus status = OFFSET_OK; // Status 0 => success.
    
    int id = 0;
    shared_ptr<ftFaceBase> base_sf;
    shared_ptr<SplineSurface> spline_appr_sf; // The spline surface approximating the surface set,
                                              // defining the parametrization.

    // For trimmed surfaces we need the topology analysis in ftChartSurface to extract the corners
    // and thus limit the domain (instead of extracting the underlying surface without trim curves).
    if ((param_sfs.size() == 1) && (param_sfs[0]->instanceType() == Class_SplineSurface))
    {
        base_sf = shared_ptr<ftFaceBase>(new ftSurface(param_sfs[0], id));
    }
    else
    {
        // We need to use a ftChartSurface.
        const double gap = epsgeo;
        const double kink = 1.0e-02;
        tpTolerances top_eps(gap, 10.0*gap, kink, 10.0*kink);
        const double approx = epsgeo;//*100.0;
        MESSAGE("Consider adding weights to any kink curves when creating the approximating spline surface!");
        const double curvature_tol = epsgeo*0.5;//*0.01; // @@sbr201706 This value should be tuned! Used to
                                                  // measure distance from point to a linear segment.
        int m = 0;
        int n = 0;
        shared_ptr<ftChartSurface> chart_sf;
        try
        {
            chart_sf = shared_ptr<ftChartSurface>(new ftChartSurface(param_sfs, top_eps,
                                                                     approx, curvature_tol,
                                                                     m, n));
        }
        catch (...)
        {
            MESSAGE("Failed creating a ftChartSurface!");
            status = OFFSET_FAILED;
            return status;
        }

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
        MESSAGE("Max error: " << max_error << ". Mean error: " << mean_error);

        // We check the spline surface defining the parameter domain of the chart surface.
        shared_ptr<ParamSurface> out_surf = chart_sf->surface();
        shared_ptr<SplineSurface> spline_out_surf =
            dynamic_pointer_cast<SplineSurface, ParamSurface>(out_surf);
        spline_appr_sf = spline_out_surf;
        if (spline_out_surf.get() == NULL)
        {
            MESSAGE("Failed creating the ChartSurface SplineSurface!");
            status = OFFSET_FAILED;
            return status;
        }

        std::ofstream fileout("tmp/chart_sf.g2");
        spline_out_surf->writeStandardHeader(fileout);
        spline_out_surf->write(fileout);
        
        base_sf = chart_sf;        
    }

#if 0
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
    MESSAGE("INFO: Boundary: curv_radius_min: " << curv_radius_min << ", curv_radius_max: " << curv_radius_max);

    if (curv_radius_min < offset_dist)
    {
        MESSAGE("The radius of curvature along the boundary is less than the offset dist!");//, not supported.");
#if 1
        //MESSAGE("Turned off early exit!");
#else
        status = SELF_INTERSECTING_BOUNDARY;
        return status;
#endif
    }
#endif

    try
    {

        EvalOffsetSurface eval_offset_sf(base_sf, offset_dist, epsgeo);

        // @@sbr201705 Should we project the inner kink edges onto the approximating surface, expressing
        // it as parameter curves? We can not expect the kink curves to be iso curves in the guide
        // surface (the spline sf approximating the surface set). We also need to handle grid points that
        // are close to such a kink (the neighbourhood depending on the kink angle and the offset
        // distance).
        vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > > kink_par_cvs;
        vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > > kink_sfs;
        vector<shared_ptr<SplineCurve> > kink_cvs_2d = eval_offset_sf.getProjKinkCurves(kink_par_cvs, kink_sfs);

        // Currently we handle iso curves only. We remove curevs which do not follow an iso line.
        const size_t num_kink_cvs = kink_cvs_2d.size();
        removeNonIsoCurves(kink_cvs_2d, kink_par_cvs, kink_sfs);
        if (kink_cvs_2d.size() < num_kink_cvs)
        {
            status = NON_ISO_KINK_CURVE; // To handle such cases we should offset with a LRSplineSurface.
            return status;
        }
    
        // Creating the initial grid.
        // Only the end parameters are set initially.
        HermiteApprEvalSurf appr_eval_sf(&eval_offset_sf, epsgeo, epsgeo);
        // @@sbr201705 Parameter value, assuming sensible relation between parameter domain and geometry. Check!
        const double no_split_dist = fabs(offset_dist)*0.5;
        appr_eval_sf.setNoSplit(kink_cvs_2d, no_split_dist);

        // We keep on refining until within the error tolerance (or the method fails, typically due to knot
        // spacing becoming too dense).
        appr_eval_sf.refineApproximation();

        // Prior to creating the surface we fetch the self intersection points.
        const HermiteGrid2D& grid = appr_eval_sf.getGrid();
        
        // We then run through the grid, removing all grid points that are within a certain distance from
        // the kink.
        // @@sbr201705 Missing the removal of grid kinks!
        MESSAGE("Missing the removal of grid kinks!");
        
        // We locate all grid points with only c0 continuity.
        vector<int> grid_kinks; // The index in the grid_.
        eval_offset_sf.gridKinks(grid, kink_cvs_2d, grid_kinks);
        
        // We check if the grid contains any self intersections. If any such self intersections are along
        // the boundary we return with an error message.
        vector<int> grid_self_int;
        vector<double> grid_self_int_curv_radius; // For the grid self intersections.
        eval_offset_sf.gridSelfIntersections(grid, grid_self_int, grid_self_int_curv_radius);
        size_t num_self_int = grid_self_int.size();
        
        // We mark all the grid nodes that are within a certain distance from a self-intersection.
        const double smoothing_dist = 1.5*std::fabs(offset_dist);
        MESSAGE("Calling getGridSmoothingNodes().");
        vector<int> grid_smoothing_nodes =
            getGridSmoothingNodes(grid, grid_self_int,
                                  //grid_self_int_curv_radius,
                                  smoothing_dist);
        vector<double> grid_smoothing_curv_radius(grid_smoothing_nodes.size(), 0.0);
        MESSAGE("Done calling getGridSmoothingNodes(), grid_smoothing_nodes.size(): " << grid_smoothing_nodes.size());

#if 1
        // If there are self intersections going from one end to the other we remove those lines from the
        // spline surface resulting from the grid. We must also update the indices in the smoothing
        // vectors.
        vector<double> iso_self_int_u, iso_self_int_v;
        getIsoSelfIntersections(grid, grid_self_int, iso_self_int_u, iso_self_int_v);
        MESSAGE("grid_self_int.size(): " << grid_self_int.size() << ", iso_self_int_u.size(): " <<
                iso_self_int_u.size() << ", iso_self_int_v.size(): " << iso_self_int_v.size());
        if ((iso_self_int_u.size() > 0 ) || (iso_self_int_v.size() > 0 ))
        {
            MESSAGE("Found iso self intersection(s)!");
            
            // We mark the grid lines as not to be used when calling getSurface().
            appr_eval_sf.removeGridLines(iso_self_int_u, iso_self_int_v);
            
            updateGridSelfInt(grid, iso_self_int_u, iso_self_int_v,
                              grid_self_int, grid_self_int_curv_radius);

            
            MESSAGE("Pre update: grid_smoothing_nodes.size(): " << grid_smoothing_nodes.size());
                // We update the indices with the iso grid lines removed. That means that the indices are
                // directly related to coef indices in the bicubic spline surface.
            updateGridSelfInt(grid, iso_self_int_u, iso_self_int_v,
                              grid_smoothing_nodes, grid_smoothing_curv_radius);
            MESSAGE("Post update: grid_smoothing_nodes.size(): " << grid_smoothing_nodes.size());

            if (grid_kinks.size() > 0)
            { // We update the kink indices as grid iso-line(s) have been removed.
                vector<double> kink_radius(grid_kinks.size(), 0.0);
                updateGridSelfInt(grid, iso_self_int_u, iso_self_int_v,
                                  grid_kinks, kink_radius);
            }
        }
#else
        MESSAGE("Turned off removing of iso self intersections!");
#endif
        
        // Creating the surface from the Bezier patches.
        bool method_failed;
        offset_sf = appr_eval_sf.getSurface(method_failed);
        if (method_failed)
        {
            MESSAGE("HermiteApprEvalSurf.getSurface() failed (suspecting knot spacing too dense)!");
            status = OFFSET_FAILED;
        }
        if (offset_sf.get() == 0)
        {
            MESSAGE("HermiteApprEvalSurf.getSurface() was NULL!");
            status = OFFSET_FAILED;
        }
        
        if ((num_self_int > 0) || (kink_cvs_2d.size() > 0))
        {
            // We calculate the radius of curvature (min) in inner grid points. We use this to
            // define areas which need smoothing, possibly in a post processing step.
            // The method refines until within the required tolerance (or aborts if a knot interval gets too small).
            if (num_self_int > 0)
            {
                MESSAGE("Self intersections found! We perform smoothing on the area affected! num_self_int: " <<
                        num_self_int);
            }
#ifndef NDEBUG
            {
                std::ofstream fileout_debug("tmp/offset_sf_selfint.g2");
                offset_sf->writeStandardHeader(fileout_debug);
                offset_sf->write(fileout_debug);
            }
#endif

            // We fetch the offset kink points by averaging the normals along the kink. The kink_pts and
            // kink_params are added to the equation system as approximation terms.
            vector<double> kink_pts, kink_params;
            fetchKinkPoints(kink_par_cvs, kink_sfs, offset_dist, epsgeo, spline_appr_sf,
                            kink_pts, kink_params);
            MESSAGE("kink_pts.size(): " << kink_pts.size() << ", kink_params.size(): " << kink_params.size());
            
            vector<double> kink_release_dist(grid_kinks.size(), offset_dist);

            vector<int> grid_released_nodes(grid_self_int.begin(), grid_self_int.end());
            grid_released_nodes.insert(grid_released_nodes.end(),
                                       grid_kinks.begin(), grid_kinks.end());
            grid_released_nodes.insert(grid_released_nodes.end(),
                                       grid_smoothing_nodes.begin(), grid_smoothing_nodes.end());
            
            vector<double> grid_curv_radius(grid_self_int_curv_radius.begin(), grid_self_int_curv_radius.end());
            grid_curv_radius.insert(grid_curv_radius.end(),
                                    kink_release_dist.begin(), kink_release_dist.end());
            grid_curv_radius.insert(grid_curv_radius.end(),
                                    grid_smoothing_curv_radius.begin(), grid_smoothing_curv_radius.end());

            MESSAGE("grid_released_nodes.size(): " << grid_released_nodes.size() <<
                    ", grid_smoothing_curv_radius.size(): " << grid_smoothing_curv_radius.size());
            vector<int> coefs_released = getCoefsReleased(offset_sf, grid_released_nodes, grid_curv_radius);
            MESSAGE("coefs_released.size(): " << coefs_released.size());
            //grid_self_int, grid_self_int_curv_radius);//, grid_kinks, kink_release_dist);

            shared_ptr<SplineSurface> smooth_offset_sf =
                getSmoothOffsetSurface(offset_sf,
                                       coefs_released,
                                       kink_pts, kink_params);
            if (smooth_offset_sf.get() != NULL)
            {
                offset_sf = smooth_offset_sf;
#if 1
                std::ofstream fileout_debug("tmp/offset_sf_smooth.g2");
                smooth_offset_sf->writeStandardHeader(fileout_debug);
                smooth_offset_sf->write(fileout_debug);
#endif
            }
        }
        else
        {
#if 1
            {
                std::ofstream fileout_debug("tmp/offset_sf_no_selfint.g2");
                offset_sf->writeStandardHeader(fileout_debug);
                offset_sf->write(fileout_debug);
            }
#endif
        }        

    }
    catch (std::exception& e)
    {
        MESSAGE("Exception was thrown: " << e.what());
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
            
            MESSAGE("curv_rad1: " << curv_rad1 << ", curv_rad2: " << curv_rad2);
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
                                                 const vector<int>& coefs_released, // Index of the released coefs.
                                                 const vector<double>& kink_pts,
                                                 const vector<double>& kink_params)
//===========================================================================
{
    shared_ptr<SplineSurface> new_offset_sf;

    if (coefs_released.size() == 0)
    {
        MESSAGE("No self intersections or kinks, method should not have been called!");
        return new_offset_sf;
    }
    
    const int num_coefs_u = offset_sf->numCoefs_u();
    const int num_coefs_v = offset_sf->numCoefs_v();
    const int dim = offset_sf->dimension();
    // vector<double> coef_curv_radius(num_coefs_u*num_coefs_v, 0.0);
    // vector<double> released_radius;
    // The offset_sf & grid_self_int are referring to the same space, i.e. the grid is the basis for the
    // Bezier patches constituting the offset_sf. We thus can recreate the grid dimensionality.
    // If a coef is within coef_change_ratio*offset_dist of self-intersection-coef it is allowed to
    // change.
    // @@sbr201704 Consider adding a small approximation weight to these terms, reduced as the coef
    // approaches a "self-intersecting-coef".
    
    vector<int> coef_known(num_coefs_u*num_coefs_v, 1);
    for (size_t ki = 0; ki < coefs_released.size(); ++ki)
    {
        coef_known[coefs_released[ki]] = 0;
    }

#ifndef NDEBUG
    {
        std::ofstream fileout_debug("tmp/coefs_known.g2");
        std::ofstream fileout_debug2("tmp/coefs_not_known.g2");
        vector<double> coefs_known, coefs_not_known;
        for (size_t ki = 0; ki < coef_known.size(); ++ki)
        {
            if (coef_known[ki])
            {
                coefs_known.insert(coefs_known.end(),
                                   offset_sf->coefs_begin() + ki*dim, offset_sf->coefs_begin() + (ki + 1)*dim);
            }
            else
            {
                coefs_not_known.insert(coefs_not_known.end(),
                                       offset_sf->coefs_begin() + ki*dim, offset_sf->coefs_begin() + (ki + 1)*dim);
            }
        }
        PointCloud3D pt_cl(coefs_known.begin(), coefs_known.size()/dim);
        vector<int> color_green(4, 0);
        color_green[1] = 255;
        color_green[3] = 255;
        ObjectHeader header(Class_PointCloud, 1, 0, color_green);
        header.write(fileout_debug);
        pt_cl.write(fileout_debug);

        PointCloud3D pt_cl2(coefs_not_known.begin(), coefs_not_known.size()/dim);
        vector<int> color_red(4, 0);
        color_red[0] = 255;
        color_red[3] = 255;
        ObjectHeader header2(Class_PointCloud, 1, 0, color_red);
        header2.write(fileout_debug2);
        pt_cl2.write(fileout_debug2);
    }
#endif
    
    shared_ptr<SplineSurface> smooth_offset_sf; // Without self intersections.

    SmoothSurf smooth_sf;
    vector<int> seem(2, 0);
    smooth_sf.attach(offset_sf, &seem[0], &coef_known[0]);
        
    double smooth_wgt = 1.0e-01;//0.0;//1.0e-2;//3;
//    const int min_der = 3;//2;
    double wgt1 = 0.0;//1.0e-02;//0.0;
    double wgt3 = 0.0;// (min_der >= 3) ? 0.5*smooth_wgt : 0.0; // Only c1 surface => wgt3 = 0.0.
    double wgt2 = (1.0 - wgt3)*smooth_wgt;
    wgt3 *= smooth_wgt;
    double weight_sum = wgt1 + wgt2 + wgt3;
    if (weight_sum > 1.0) { // Making sure the sum of the weights is at most 1.0.
        wgt1 /= weight_sum;
        wgt2 /= weight_sum;
        wgt3 /= weight_sum;
    }
    smooth_sf.setOptimize(wgt1, wgt2, wgt3);

    const double approxweight = 1.0 - (wgt1 + wgt2 + wgt3); // In the range [0.0, 1.0].

    // The grid kink curve parametrization is described as the projection onto the guide surface (the
    // approximating spline surface which defines the parametrization for the offset surface).
    if (kink_pts.size() > 0)
    {
        vector<double> pnt_wgts(kink_params.size()/2, 1.0);
        smooth_sf.setLeastSquares(kink_pts,
                                  kink_params,
                                  pnt_wgts,
                                  approxweight);
    }
    
    double wgt_orig = 0.1*approxweight;
    wgt_orig = 0.0;
    if (wgt_orig > 0.0)
    {
        smooth_sf.approxOrig(wgt_orig);
    }

    // The article "RILU preconditioning; a computational study" by Tveito & Bruaset concludes that
    // values close to 1.0 is normally a good choice. Not our experience, at least not for current test
    // suite.
    const double omega = 0.01;
    smooth_sf.setRelaxParam(omega);

    int status = smooth_sf.equationSolve(new_offset_sf);
    if (status != 0)
    {
        MESSAGE("SmoothSurf failed solving the equation!");
    } 
       
    return new_offset_sf;
}


void getIsoSelfIntersections(const HermiteGrid2D& grid, const vector<int>& grid_self_int,
                             vector<double>& iso_self_int_u, vector<double>& iso_self_int_v)
{
    const int MM = grid.size1();
    const int NN = grid.size2();

    vector<double> knots_u = grid.getKnots(true);
    vector<double> knots_v = grid.getKnots(false);

    vector<int> grid_u_mult(MM, 0);
    vector<int> grid_v_mult(NN, 0);
    // We're assuming the elements of grid_self_int are unique.
    for (size_t ki = 0; ki < grid_self_int.size(); ++ki)
    {
        // We find the position in the grid.
        int ki_mm = grid_self_int[ki]%MM;
        int ki_nn = grid_self_int[ki]/MM;
        ++grid_u_mult[ki_mm];
        ++grid_v_mult[ki_nn];
    }

    for (size_t ki = 0; ki < grid_u_mult.size(); ++ki)
    {
        if (grid_u_mult[ki] == NN)
        {
            MESSAGE("We encountered a grid self intersection which is transversal! grid_u: " << ki);
            iso_self_int_u.push_back(knots_u[ki]);
        }
    }

    for (size_t ki = 0; ki < grid_v_mult.size(); ++ki)
    {
        if (grid_v_mult[ki] == MM)
        {
            MESSAGE("We encountered a grid self intersection which is transversal! grid_v: " << ki);
            iso_self_int_v.push_back(knots_v[ki]);
        }
    }

}

#if 1
void updateGridSelfInt(const HermiteGrid2D& grid,
                       const vector<double>& grid_remove_u, const vector<double>& grid_remove_v,
                       vector<int>& grid_self_int, vector<double>& radius_of_curv)
{
    MESSAGE("INFO: grid_self_int.size(): " << grid_self_int.size());

#if 1


    
#endif
    
    vector<int> new_grid_self_int;
    vector<double> new_radius_of_curv;
    const int MM = grid.size1();
    const int NN = grid.size2();
    const int MM_red = grid.size1() - grid_remove_u.size();
    const int NN_red = grid.size2() - grid_remove_v.size();
    const double knot_tol = 1.0e-14;
    vector<double> knots_u = grid.getKnots(true);
    vector<double> knots_v = grid.getKnots(false);

    MESSAGE("INFO: MM: " << MM << ", MM_red: " << MM_red << ", NN: " << NN << ", NN_red: " << NN_red);

    for (size_t ki = 0; ki < grid_self_int.size(); ++ki)
    {
        // We first find the mm & nn pos in the full grid.
        int ki_mm = grid_self_int[ki]%MM;
        int ki_nn = grid_self_int[ki]/MM;

        // We then see how many grid elements to the left that has been removed.
        auto rem_u_iter = grid_remove_u.begin();
        while ((rem_u_iter != grid_remove_u.end()) && (*rem_u_iter + knot_tol < knots_u[ki_mm]))
        {
            ++rem_u_iter;
        }

        if ((rem_u_iter != grid_remove_u.end()) && (fabs(*rem_u_iter - knots_u[ki_mm]) < knot_tol))
        {
            continue;
        }
        
        auto rem_v_iter = grid_remove_v.begin();
        while ((rem_v_iter != grid_remove_v.end()) && (*rem_v_iter + knot_tol < knots_v[ki_nn]))
        {
            ++rem_v_iter;
        }

        if ((rem_v_iter != grid_remove_v.end()) && (fabs(*rem_v_iter - knots_v[ki_nn]) < knot_tol))
        {
            continue;
        }

        int new_mm = ki_mm - (rem_u_iter - grid_remove_u.begin());
        int new_nn = ki_nn - (rem_v_iter - grid_remove_v.begin());

        int new_ind = new_nn*MM_red + new_mm;
        new_grid_self_int.push_back(new_ind);
        new_radius_of_curv.push_back(radius_of_curv[ki]);
    }

    MESSAGE("INFO: new_grid_self_int.size(): " << new_grid_self_int.size());

    grid_self_int = new_grid_self_int;
    radius_of_curv = new_radius_of_curv;
}
#endif


void fetchKinkPoints(const vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >& kink_par_cvs,
                     const vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > >& kink_sfs,
                     double offset, double epsgeo, shared_ptr<SplineSurface> spline_sf,
                     vector<double>& kink_pts, vector<double>& kink_params)
{
    MESSAGE("INFO: kink_par_cvs.size(): " << kink_par_cvs.size());
    // We sample along one of the curves, fetching the corresponding points on the neighbour curve.
    // We keep on bisecting intervals until the segment length is within the requirement.
    const double max_segment_length = epsgeo;
    for (size_t ki = 0; ki < kink_par_cvs.size(); ++ki)
    {
        const ParamCurve* pcv1 = kink_par_cvs[ki].first.get();
        const ParamCurve* pcv2 = kink_par_cvs[ki].second.get();
        const ParamSurface* sf1 = kink_sfs[ki].first.get();
        const ParamSurface* sf2 = kink_sfs[ki].second.get();
        const double tmin1 = pcv1->startparam();
        const double tmax1 = pcv1->endparam();
        const double tmin2 = pcv2->startparam();
        const double tmax2 = pcv2->endparam();
        Point pcv1_start_pt = pcv1->point(tmin1);
        Point pcv1_end_pt = pcv1->point(tmax1);
        Point pcv2_start_pt = pcv2->point(tmin2);
        Point pcv2_end_pt = pcv2->point(tmax2);
        Point sf1_start_pt = sf1->point(pcv1_start_pt[0], pcv1_start_pt[1]);
        Point sf1_end_pt = sf1->point(pcv1_end_pt[0], pcv1_end_pt[1]);
        Point sf2_start_pt = sf2->point(pcv2_start_pt[0], pcv2_start_pt[1]);
        Point sf2_end_pt = sf2->point(pcv2_end_pt[0], pcv2_end_pt[1]);
        double dist = sf1_start_pt.dist(sf2_start_pt) + sf1_end_pt.dist(sf2_end_pt);
        double dist_opp = sf1_start_pt.dist(sf2_end_pt) + sf1_end_pt.dist(sf2_start_pt);
        bool opp_dir = (dist_opp < dist);
        int num_samples = 1000;//3; // We start with 3 samples.
        double tstep1 = (tmax1 - tmin1)/(double)(num_samples - 1);
        double tstep2 = (tmax2 - tmin2)/(double)(num_samples - 1);

        for (size_t kj = 0; kj < num_samples; ++kj)
        {
            double tpar1 = tmin1 + (double)kj*tstep1;
            Point par_pt1 = pcv1->point(tpar1);
            Point sf_pt1 = sf1->point(par_pt1[0], par_pt1[1]);
            double tpar2 = (opp_dir) ? tmax2 - (double)kj*tstep2 : tmin2 + (double)kj*tstep2;
            Point par_pt2 = pcv1->point(tpar2);
            // We may assume that the kink cvs follow the boundary.
            double clo_u, clo_v, clo_dist;
            Point clo_pt;
            sf2->closestBoundaryPoint(sf_pt1, clo_u, clo_v, clo_pt, clo_dist, epsgeo, NULL, &par_pt2[0]);

            // We find the average normal.
            Point normal1, normal2;
            sf1->normal(normal1, par_pt1[0], par_pt1[1]);
            sf2->normal(normal2, clo_u, clo_v);
            Point avg_normal = 0.5*(normal1 + normal2);
            avg_normal.normalize();
            Point offset_pt = sf_pt1 + offset*avg_normal;

            kink_pts.insert(kink_pts.end(), offset_pt.begin(), offset_pt.end());
            // We must also locate the closest point on spline_sf_, which gives us the parameter value.
            spline_sf->closestPoint(sf_pt1, clo_u, clo_v, clo_pt, clo_dist, epsgeo);
            kink_params.push_back(clo_u);
            kink_params.push_back(clo_v);
        }
    }
#ifndef NDEBUG
    {
        std::ofstream fileout_debug("tmp/kink_pts.g2");
        PointCloud3D pt_cl(kink_pts.begin(), kink_pts.size()/3);
        vector<int> color_blue(4, 0);
        color_blue[2] = 255;
        color_blue[3] = 255;
        ObjectHeader header(Class_PointCloud, 1, 0, color_blue);
        header.write(fileout_debug);
        pt_cl.write(fileout_debug);
    }
#endif
}

    
vector<int> getGridSmoothingNodes(const HermiteGrid2D& grid,
                                  const vector<int>& grid_self_int,
                                  double smoothing_dist)
{
    // The indices refer to the position in the grid. The indices must be converted to coefs_known if
    // used for setting the coefs_known in the smoothing function.

    MESSAGE("Replace vector with set!"); // @@sbr201706 The find function takes too long time.
    
    vector<int> grid_smoothing_nodes;
    MESSAGE("getGridSmoothingNodes(): Under construction!");

    // We check the distance to all the neighbour coefs.
    const int MM = grid.size1();
    const int NN = grid.size2();
    const int max_nb_steps = std::max(MM, NN);
    const vector<Point>& grid_data = grid.getData();
    int total_num_changed = 0;
    double min_dist_global = MAXDOUBLE;
    for (size_t ki = 0; ki < grid_self_int.size(); ++ki)
    {
        const int grid_ind = grid_self_int[ki];
        const Point& grid_node = grid_data[grid_ind*4];
        const int ki_mm = grid_ind%MM;
        const int ki_nn = grid_ind/MM;

        double min_dist = MAXDOUBLE;
        for (int kj = 0; kj < max_nb_steps; ++kj)
        {
            //std::cout << "max_nb_steps: " << max_nb_steps << ", kj: " << kj << std::endl;
            int num_changed = 0;
            int num_affected = 0;
            const int edge_size = 2*kj + 2; // Size of the neighborhood edge (in the quadrat).
            int ll_u = ki_mm - kj - 1; // Lower left of the quadrat, u index.
            int ll_v = ki_nn - kj - 1; // Lower left of the quadrat, v index.
            for (int kk = 0; kk < 4; ++kk)
            {
                for (int kl = 0; kl < edge_size; ++kl)
                {
                    int nb_u = (kk < 2) ? (kl + kk) : (kk == 2) ? 0 : edge_size;
                    int nb_v = (kk > 1) ? (kl + 3 - kk) : ((kk == 0) ? 0 : edge_size);
                    // std::cout << "ind_u: " << ind_u << ", ind_v: " << ind_v <<
                    //     ", nb_u: " << nb_u << ", nb_v: " << nb_v <<
                    //     ", kj: " << kj << ", kk: " << kk << ", kl: " << kl << std::endl;

                    int nb_ind = (ll_v + nb_v)*MM + (ll_u + nb_u);
                    if ((nb_ind < 0) || (nb_ind > MM*NN - 1))
                        continue;
                    const Point& nb_node = grid_data[nb_ind*4];
                    const double dist = grid_node.dist(nb_node);
                    if (dist < smoothing_dist)
                    {
                        ++num_affected;
                        // We only add the node if it is not already present and it is not a self intersection node.
                        if ((std::find(grid_smoothing_nodes.begin(), grid_smoothing_nodes.end(), nb_ind) ==
                             grid_smoothing_nodes.end()) &&
                            (std::find(grid_self_int.begin(), grid_self_int.end(), nb_ind) == grid_self_int.end()))
                        {
                            if (dist < min_dist)
                            {
                                min_dist = dist;
                            }
                            grid_smoothing_nodes.push_back(nb_ind);
                            ++num_changed;
                            // std::cout << "ki: " << ki << ", dist: " << dist << ", nb_ind: " << nb_ind << 
                            //     ", coef_curv_radius[ki]: " << coef_curv_radius[ki] << std::endl;
                            // std::cout << "coef: " << coef << ", nb_coef: " << nb_coef << std::endl;
                        }
                    }
                }
            }

            // std::cout << "num_changed: " << num_changed << std::endl;
            if (num_affected == 0)
            { // When we encounter a region without any coefs close enough we stop searching, moving on to next coef.
                break;
            }
        }
//        MESSAGE("INFO: min_dist: " << min_dist);
        if (min_dist < min_dist_global)
        {
            min_dist_global = min_dist;
        }
    }

    MESSAGE("INFO: min_dist_global: " << min_dist_global << ", smoothing_dist: " << smoothing_dist);
    
    // We sort the nodes.
    std::sort(grid_smoothing_nodes.begin(), grid_smoothing_nodes.end());
    
    MESSAGE("INFO: total_num_changed (nodes released for smoothing): " << grid_smoothing_nodes.size());

    return grid_smoothing_nodes;
}

vector<int> getCoefsReleased(shared_ptr<SplineSurface> offset_sf,
                             const vector<int>& grid_self_int,
                             const vector<double>& radius_of_curv)//,
                             // const vector<int>& grid_kinks,
                             // const vector<double>& kink_release_dist)
{
    vector<int> coefs_released;
    vector<double> released_radius;
    
    const int num_coefs_u = offset_sf->numCoefs_u();
    const int num_coefs_v = offset_sf->numCoefs_v();
    const int dim = offset_sf->dimension();
    // The offset_sf & grid_self_int are referring to the same space, i.e. the grid is the basis for the
    // Bezier patches constituting the offset_sf. We thus can recreate the grid dimensionality.
    const int grid_mm = (num_coefs_u)/2;
    const int grid_nn = (num_coefs_v)/2;
    MESSAGE("grid_mm: " << grid_mm << ", grid_nn: " << grid_nn);
    // If a coef is within coef_change_ratio*offset_dist of self-intersection-coef it is allowed to
    // change.
    // @@sbr201704 Consider adding a small approximation weight to these terms, reduced as the coef
    // approaches a "self-intersecting-coef".
    for (size_t ki = 0; ki < grid_self_int.size(); ++ki)
    {
        // We find the position in the grid.
        int ki_mm = grid_self_int[ki]%grid_mm;
        int ki_nn = grid_self_int[ki]/grid_mm;
        // std::cout << "ki_mm: " << ki_mm << ", ki_nn: " << ki_nn <<
        //     ", radius_of_curv[ki]: " << radius_of_curv[ki] << std::endl;
        coefs_released.push_back((2*ki_nn)*num_coefs_u + (2*ki_mm));
        coefs_released.push_back((2*ki_nn)*num_coefs_u + (2*ki_mm+1));
        coefs_released.push_back((2*ki_nn+1)*num_coefs_u + (2*ki_mm));
        coefs_released.push_back((2*ki_nn+1)*num_coefs_u + (2*ki_mm+1));

        released_radius.insert(released_radius.end(), 4, radius_of_curv[ki]);

        size_t ind = (2*ki_nn+1)*num_coefs_u + (2*ki_mm+1);
        if (ind > num_coefs_u*num_coefs_v - 1)
        {
            MESSAGE("ERROR: Outside index range, ind = " << ind);
        }
    }

#if 0
    // We run through the grid_kinks vector and release the corresponding coefs.
    for (size_t ki = 0; ki < grid_kinks.size(); ++ki)
    {
        // The grid_kinks values relate to the position in the grid_. The offset_sf was create from the
        // grid and is assumed to not have been refined.
        int ki_u = grid_kinks[ki]%grid_mm;
        int ki_v = grid_kinks[ki]/grid_mm;
        coefs_released.push_back((2*ki_v)*num_coefs_u + (2*ki_u));
        coefs_released.push_back((2*ki_v)*num_coefs_u + (2*ki_u+1));
        coefs_released.push_back((2*ki_v+1)*num_coefs_u + (2*ki_u));
        coefs_released.push_back((2*ki_v+1)*num_coefs_u + (2*ki_u+1));

        released_radius.insert(released_radius.end(), 4, kink_release_dist[ki]);
    }
#endif
    
    vector<int> coef_known(num_coefs_u*num_coefs_v, 1);
    for (size_t ki = 0; ki < coefs_released.size(); ++ki)
    {
        coef_known[coefs_released[ki]] = 0;
    }

    // We then check the distance to all the neighbour coefs.
    // Release all coefs within a given distance.
    const int max_nb_steps = std::max(num_coefs_u/2, num_coefs_v/2);
    const double coef_change_ratio = 1.0;// 0.5;
    int total_num_changed = 0;
    vector<int> coefs_released_smoothing;
    for (size_t ki = 0; ki < coefs_released.size(); ++ki) //coef_known.size(); ++ki)
    { // This means that the coef corresponds to a self intersection in the offset surface.
        // We check all neighboring coefs if they are within a certain distance. We enlarge the
        // region until we're no longer within the "neighborhood sphere".
        const int ci = coefs_released[ki];
        Point coef(offset_sf->coefs_begin() + ci*dim, offset_sf->coefs_begin() + (ci + 1)*dim);
        size_t ind_u = ci%num_coefs_u;
        size_t ind_v = ci/num_coefs_u;
        for (int kj = 0; kj < max_nb_steps; ++kj)
        {
            //std::cout << "max_nb_steps: " << max_nb_steps << ", kj: " << kj << std::endl;
            int num_changed = 0;
            const int edge_size = 2*kj + 2; // Size of the neighborhood edge (in the quadratic).
            int ll_u = ind_u - kj - 1; // Lower left of the quadratic, u index.
            int ll_v = ind_v - kj - 1; // Lower left of the quadratic, v index.
            for (int kk = 0; kk < 4; ++kk)
            {
                for (int kl = 0; kl < edge_size; ++kl)
                {
                    int nb_u = (kk < 2) ? (kl + kk) : (kk == 2) ? 0 : edge_size;
                    int nb_v = (kk > 1) ? (kl + 3 - kk) : ((kk == 0) ? 0 : edge_size);
                    // std::cout << "ind_u: " << ind_u << ", ind_v: " << ind_v <<
                    //     ", nb_u: " << nb_u << ", nb_v: " << nb_v <<
                    //     ", kj: " << kj << ", kk: " << kk << ", kl: " << kl << std::endl;

                    int nb_ind = (ll_v + nb_v)*num_coefs_u + (ll_u + nb_u);
                    if ((nb_ind < 0) || (nb_ind > coef_known.size() - 1))
                        continue;
                    if (coef_known[nb_ind] == 1)
                    {
                        // Initially we set coef_known[nb_ind] = -1, to denote that it should be changed.
                        // Setting it to 0 after the loop.
                        Point nb_coef(offset_sf->coefs_begin() + nb_ind*dim,
                                      offset_sf->coefs_begin() + (nb_ind + 1)*dim);
                        double dist = coef.dist(nb_coef);
                        if (dist < fabs(released_radius[ki])*coef_change_ratio)
                        {
                            coef_known[nb_ind] = -1;
                            coefs_released_smoothing.push_back(nb_ind);
                            ++num_changed;
                            ++total_num_changed;
                            // std::cout << "ki: " << ki << ", dist: " << dist << ", nb_ind: " << nb_ind << 
                            //     ", coef_curv_radius[ki]: " << coef_curv_radius[ki] << std::endl;
                            // std::cout << "coef: " << coef << ", nb_coef: " << nb_coef << std::endl;
                        }
                    }
                }
            }
                
            // std::cout << "num_changed: " << num_changed << std::endl;
            if (num_changed == 0)
            {
                break;
            }
        }
    }

    coefs_released.insert(coefs_released.end(), coefs_released_smoothing.begin(), coefs_released_smoothing.end());
    
    MESSAGE("INFO: total_num_changed (coefs released for smoothing): " << total_num_changed);
    
    // We write to file the coefs.
    vector<double> self_int_coefs, released_coefs, locked_coefs;
    for (size_t ki = 0; ki < coef_known.size(); ++ki)
    {
        Point coef(offset_sf->coefs_begin() + ki*dim, offset_sf->coefs_begin() + (ki + 1)*dim);
        if (coef_known[ki] == -1)
        {
            released_coefs.insert(released_coefs.end(), coef.begin(), coef.end());
            coef_known[ki] = 0;
        }
        else if (coef_known[ki] == 0)
        {
            self_int_coefs.insert(self_int_coefs.end(), coef.begin(), coef.end());

        }
        else if (coef_known[ki] == 1)
        {
            locked_coefs.insert(locked_coefs.end(), coef.begin(), coef.end());
        }
        else
        {
            MESSAGE("Unexpected value for coef_known[ki]!");
        }
    }

#if 1
    {
        // We write to file the two sets. Using green color for no self int, red color for self int.
        MESSAGE("Writing to file the self int and no self int points.");
        
        std::ofstream fileout_debug("tmp/coefs_self_int.g2");
        if (self_int_coefs.size() > 0)
        {
            PointCloud3D pt_cl(self_int_coefs.begin(), self_int_coefs.size()/3);
            vector<int> color_red(4, 0);
            color_red[0] = 255;
            color_red[3] = 255;
            ObjectHeader header(Class_PointCloud, 1, 0, color_red);
            header.write(fileout_debug);
            pt_cl.write(fileout_debug);
        }            

        std::ofstream fileout_debug2("tmp/coefs_no_self_int.g2");
        if (locked_coefs.size() > 0)
        {
            PointCloud3D pt_cl(locked_coefs.begin(), locked_coefs.size()/3);
            vector<int> color_green(4, 0);
            color_green[1] = 255;
            color_green[3] = 255;
            ObjectHeader header(Class_PointCloud, 1, 0, color_green);
            header.write(fileout_debug2);
            pt_cl.write(fileout_debug2);
        }            

        std::ofstream fileout_debug3("tmp/coefs_released.g2");
        if (released_coefs.size() > 0)
        {
            PointCloud3D pt_cl(released_coefs.begin(), released_coefs.size()/3);
            vector<int> color_brown(4, 0);
            color_brown[0] = 255;
            color_brown[1] = 255;
            color_brown[3] = 255;
            ObjectHeader header(Class_PointCloud, 1, 0, color_brown);
            header.write(fileout_debug3);
            pt_cl.write(fileout_debug3);
        }            
    }
#endif

    return coefs_released;
}

    
void removeNonIsoCurves(vector<shared_ptr<SplineCurve> >& kink_cvs_2d,
                        vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >& par_cvs,
                        vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > >& sfs)
{
    // We use a ratio of the delta_u & delta_v.
    // @@sbr201706 The actual deltas may also need to be considered.
    const double box_iso_ratio = 0.01;
    for (size_t ki = 0; ki < kink_cvs_2d.size(); ++ki)
    {
        const BoundingBox& bd_box = kink_cvs_2d[ki]->boundingBox();
        double const delta_u = bd_box.high()[0] - bd_box.low()[0];
        double const delta_v = bd_box.high()[1] - bd_box.low()[1];
        const double ratio = (delta_u < delta_v) ? delta_u/delta_v : delta_v/delta_u;
        if (ratio > box_iso_ratio)
        {
            // We need LR B-splines for these cases.
            MESSAGE("The kink curve is not an iso curve, removing it!");
            kink_cvs_2d.erase(kink_cvs_2d.begin() + ki);
            par_cvs.erase(par_cvs.begin() + ki);
            sfs.erase(sfs.begin() + ki);
            --ki;
        }
    }
}

    
} // namespace Go

