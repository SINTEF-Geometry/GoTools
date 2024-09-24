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


#include <time.h>
#include <istream>
#include <fstream>
#include <sstream>
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/utils/ClosestPointUtils.h"
#include "GoTools/utils/omp.h"

using namespace std;
using namespace Go;
using namespace Go::boxStructuring;


// #define LOG_CLOSEST_POINTS


namespace Go
{


  shared_ptr<BoundingBoxStructure> preProcessClosestVectors(const vector<shared_ptr<GeomObject> >& surfaces, double par_len_el)
  {
#ifdef LOG_CLOSEST_POINTS
    clock_t t_before = clock();
#endif

    shared_ptr<BoundingBoxStructure> structure(new BoundingBoxStructure());  // The final object to be returned, holding all preprocessing information
    BoundingBox bigbox(3);   // A bounding box that eventually will contain all segment geometry space bounding boxes

    // First we run through all surfaces in the model. For each surface we
    // - Split its parameter domain into segments and add them to the structure object
    // - For bounded surfaces
    //   * Determine if segments are entirely inside the parameter domain
    //   * Find points on the surface near the noundary used to get an upperlimit on the distance from a point to the surface
    //   * Get a polygon inside the parameter domain, and store polygon information on each segment
    vector<shared_ptr<GeomObject> >::const_iterator surf_end = surfaces.end();
    int srf_idx = 0;

    for (vector<shared_ptr<GeomObject> >::const_iterator surf_it = surfaces.begin(); surf_it != surf_end; ++surf_it)
      {

	shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(*surf_it);
	if (paramSurf.get())
	  {

	    // We have a parameteric surface, add its SurfaceData element to the structure
	    shared_ptr<SurfaceData> surf_data(new SurfaceData(paramSurf));
	    structure->addSurface(surf_data);

	    // underlying_paramSurf will point to the underlying surface in case of a bounded surface,
	    // otherwise it will be the same as paramSurf
	    shared_ptr<ParamSurface> underlying_paramSurf = paramSurf;
	    shared_ptr<BoundedSurface> asBounded = dynamic_pointer_cast<BoundedSurface>(paramSurf);
	    if (asBounded.get())
	      underlying_paramSurf = asBounded->underlyingSurface();


	    // Split the parameter domain into segments.

	    vector<vector<double> > segment_pars(2);  // The segment boundary parameters, first for u-direction, second for v-direction

	    // Parameter domain splitting when the (underlying) surface is a spline surface. We use internal knots as sepration values between the segments
	    shared_ptr<SplineSurface> splineSurf = dynamic_pointer_cast<SplineSurface>(underlying_paramSurf);
	    shared_ptr<ElementarySurface> elSurf = dynamic_pointer_cast<ElementarySurface>(underlying_paramSurf);

	    if (splineSurf.get())
	      {

		// Raise knot multiplicities to get C0-continuities at segment seams,
		// to get Bezier coefficients for calculating the boundary boxes
		shared_ptr<SplineSurface> surf_copy(splineSurf->clone());
		vector<vector<int> > seg_ctrl_pos(2);
		int order_u = surf_copy->order_u();
		int order_v = surf_copy->order_v();

		// Make knot insertions and store segment boundary parameters. The loop is run twice,
		// once for each parameter direction
		for (int dir = 0; dir < 2; ++dir)
		  {
		    int order = (dir == 0) ? order_u : order_v;
		    BsplineBasis basis = (dir == 0) ? surf_copy->basis_u() : surf_copy->basis_v();
		    vector<double> new_knots;

		    int seg = 0;
		    for (vector<double>::const_iterator it = basis.begin(); it != basis.end(); ++seg)
		      {
			double knot = *it;
			int cnt = 0;
			for (;it != basis.end() && (*it == knot); ++it, ++cnt);
			int needed = order - 1;       // For interior knots, multiplicity (order - 1) is enough...
			if (seg == 0 || it == basis.end())
			  ++needed;                   // ...but for end knots, multiplicity must be 'order'
			for (; cnt < needed; ++cnt)
			  new_knots.push_back(knot);

			segment_pars[dir].push_back(knot);
			if (seg == 0)
			  seg_ctrl_pos[dir].push_back(0);
			else if (it != basis.end())
			  seg_ctrl_pos[dir].push_back(seg_ctrl_pos[dir][seg - 1] + cnt);
		      }

		    if (dir == 0)
		      surf_copy->insertKnot_u(new_knots);
		    else
		      surf_copy->insertKnot_v(new_knots);
		  }

		int n_segs_u = (int)seg_ctrl_pos[0].size();
		int n_segs_v = (int)seg_ctrl_pos[1].size();
		surf_data->setSegments(n_segs_u, n_segs_v);
		int n_coefs_u = surf_copy->numCoefs_u();

		// Run through all segments and calculate boundary box for each
		vector<double>::const_iterator ctrl_it_v = surf_copy->ctrl_begin();  // Pointer to the first control point in segment (0,j) in the j-loop below
		for (int j = 0; j < n_segs_v; ++j)
		  {
		    double start_v = segment_pars[1][j];
		    double end_v = segment_pars[1][j+1];
		    if (j > 0)
		      ctrl_it_v += (seg_ctrl_pos[1][j] - seg_ctrl_pos[1][j-1]) * n_coefs_u * 3;

		    vector<double>::const_iterator ctrl_it_u = ctrl_it_v;  // Pointer to the first control point in segment (i,j) in the i-loop below
		    for (int i = 0; i < n_segs_u; ++i)
		      {
			if (i > 0)
			  ctrl_it_u += (seg_ctrl_pos[0][i] - seg_ctrl_pos[0][i-1]) * 3;

			// Collect control points
			vector<double> ctrl_pts;
			for (int ctrl_pos = 0, k = 0; k < order_v; ++k, ctrl_pos += (n_coefs_u - order_u) * 3)
			  for (int l = 0; l < order_u; ++l)
			    for (int m = 0; m < 3; ++m, ++ctrl_pos)
			      ctrl_pts.push_back(ctrl_it_u[ctrl_pos]);

			// Create bounding box, add the elemnt to the structure, and update the big bounding box
			BoundingBox bb;
			bb.setFromArray(ctrl_pts.begin(), ctrl_pts.end(), 3);

			Array<double, 2> ll, ur;
			ll[0] = segment_pars[0][i];
			ll[1] = start_v;
			ur[0] = segment_pars[0][i+1];
			ur[1] = end_v;
			shared_ptr<SubSurfaceBoundingBox> box(new SubSurfaceBoundingBox(surf_data, i, j, bb, shared_ptr<RectDomain>(new RectDomain(ll, ur))));
			structure->addBox(box);

			bigbox.addUnionWith(bb);
		      }
		  }
	      }

	    // Parameter domain splitting when the (underlying) surface is an elementary surface. We use
	    // a regular splitting on the parameter domain
	    else if (elSurf.get())
	      {

		// Currently, we only handle the case of either a plane, a sphere or a cylinder
		shared_ptr<Sphere> sphSurf = dynamic_pointer_cast<Sphere>(elSurf);
		shared_ptr<Cylinder> cylSurf = dynamic_pointer_cast<Cylinder>(elSurf);
		shared_ptr<Plane> planeSurf = dynamic_pointer_cast<Plane>(elSurf);

		RectDomain big_rd = paramSurf->containingDomain();
		double umin = big_rd.umin();
		double umax = big_rd.umax();
		double vmin = big_rd.vmin();
		double vmax = big_rd.vmax();

		// Get the maximum lengths of iso-curves on the (underlying) surface in each parameter direction
		double len_u = 0.0;
		double len_v = 0.0;
		if (elSurf->isSwapped())
		  {
		    swap(umin, vmin);
		    swap(umax, vmax);
		  }
		if (sphSurf.get())
		  {
		    double rad = sphSurf->getRadius();
		    len_u = rad * (umax-umin) * max(cos(vmin), cos(vmax));
		    len_v = rad * (vmax-vmin);
		  }
		else if (cylSurf.get())
		  {
		    double rad = cylSurf->getRadius();
		    len_u = rad * (umax-umin);
		    len_v = vmax-vmin;
		  }
		else if (planeSurf.get())
		  {
		    len_u = umax-umin;
		    len_v = vmax-vmin;
		  }
		if (elSurf->isSwapped())
		  {
		    swap(umin, vmin);
		    swap(umax, vmax);
		    swap(len_u, len_v);
		  }

		// Get the number of segments and segment boundary parameters
		int n_segs_u = (int)(0.5 + len_u / par_len_el);
		int n_segs_v = (int)(0.5 + len_v / par_len_el);
		if (n_segs_u < 4)
		  n_segs_u = 4;
		if (n_segs_v < 4)
		  n_segs_v = 4;
		surf_data->setSegments(n_segs_u, n_segs_v);

		double step, par;
		step = (umax - umin) / (double)n_segs_u;
		par = umin;
		for (int i = 0; i <= n_segs_u; ++i, par += step)
		  segment_pars[0].push_back(par);
		step = (vmax - vmin) / (double)n_segs_v;
		par = vmin;
		for (int i = 0; i <= n_segs_v; ++i, par += step)
		  segment_pars[1].push_back(par);

		// For each segment, create bounding box, add the segment to the structure, and update the big bounding box
		shared_ptr<ElementarySurface> surf_copy(elSurf->clone());
		for (int j = 0; j < n_segs_v; ++j)
		  {
		    double seg_vmin = segment_pars[1][j];
		    double seg_vmax = segment_pars[1][j+1];
		    for (int i = 0; i < n_segs_u; ++i)
		      {
			double seg_umin = segment_pars[0][i];
			double seg_umax = segment_pars[0][i+1];
			surf_copy->setParameterBounds(seg_umin, seg_vmin, seg_umax, seg_vmax);
			BoundingBox bb = surf_copy->boundingBox();
			Array<double, 2> ll, ur;
			ll[0] = seg_umin;
			ll[1] = seg_vmin;
			ur[0] = seg_umax;
			ur[1] = seg_vmax;
			shared_ptr<SubSurfaceBoundingBox> box(new SubSurfaceBoundingBox(surf_data, i, j, bb, shared_ptr<RectDomain>(new RectDomain(ll, ur))));
			structure->addBox(box);

			bigbox.addUnionWith(bb);
		      }
		  }

	      }
            else
            {
                // THROW("Surface type not supported! Type: " << underlying_paramSurf->instanceType());
                // We add a fallback solution for the general case.
                int n_segs_u = 4;
                int n_segs_v = 4;
		surf_data->setSegments(n_segs_u, n_segs_v);

		shared_ptr<ParamSurface> surf_copy(underlying_paramSurf->clone());

		RectDomain big_rd = paramSurf->containingDomain();
		double umin = big_rd.umin();
		double umax = big_rd.umax();
		double vmin = big_rd.vmin();
		double vmax = big_rd.vmax();

                const double domain_diag = big_rd.diagLength();
                bool infinite_domain = (domain_diag > 1.0e06);
                if (infinite_domain)
                {
                    // @@sbr201805 We must limit the domain of the underlying surface. For this we can
                    // use the trim curves of the affiliated dounbded surface.
                    MESSAGE("Large domain! domain_diag = " << domain_diag);
                }

		double step, par;
		step = (umax - umin) / (double)n_segs_u;
		par = umin;
		for (int i = 0; i <= n_segs_u; ++i, par += step)
		  segment_pars[0].push_back(par);
		step = (vmax - vmin) / (double)n_segs_v;
		par = vmin;
		for (int i = 0; i <= n_segs_v; ++i, par += step)
		  segment_pars[1].push_back(par);

		for (int j = 0; j < n_segs_v; ++j)
		  {
		    double seg_vmin = segment_pars[1][j];
		    double seg_vmax = segment_pars[1][j+1];
		    for (int i = 0; i < n_segs_u; ++i)
		      {
			double seg_umin = segment_pars[0][i];
			double seg_umax = segment_pars[0][i+1];

                        // We must limit the boundary for the infinte case.
                        
			//surf_copy->setParameterBounds(seg_umin, seg_vmin, seg_umax, seg_vmax);
			BoundingBox bb = surf_copy->boundingBox();
			Array<double, 2> ll, ur;
			ll[0] = seg_umin;
			ll[1] = seg_vmin;
			ur[0] = seg_umax;
			ur[1] = seg_vmax;
			shared_ptr<SubSurfaceBoundingBox> box(new SubSurfaceBoundingBox(surf_data, i, j, bb, shared_ptr<RectDomain>(new RectDomain(ll, ur))));
			structure->addBox(box);

			bigbox.addUnionWith(bb);
		      }
		  }

            }

	    // For bounded surfaces, we now
	    //   1 Determine if segments are entirely inside the parameter domain
	    //   2 Find points on the surface near the noundary used to get an upperlimit on the distance from a point to the surface
	    //   3 Get a polygon inside the parameter domain, and store polygon information on each segment
	    if (asBounded.get())
	      {
		CurveBoundedDomain par_dom = asBounded->parameterDomain();
		int n_segs_u = surf_data->segs_u();
		int n_segs_v = surf_data->segs_v();
		int first_bb_idx = structure->n_boxes() - n_segs_u * n_segs_v;
		double eps = 1.0e-8;

		//  ****************
		//  1  Determine if the parameter domains of the bounding boxes are entirely inside the parameter domain limited by the boundary curves
		//  ****************

		// As a default, all boxes are marked as inside
		// We test when the boundary curve crosses the segment boundary parameters for each parameter direction to find
		// boxes that cannot be entirely inside. For each parameter in one direction, the variable 'inside_intervals' will hold the
		// the parameter intervals in the other parameter direction being on the inside. Since we are looking for the segments that might
		// hit the outside, we run through the outside intervals, which (if the u-direction is the fixed parameter directio) are
		//
		// [v_min,                         inside_intervals[0].first]
		// [inside_intervals[0].second,    inside_intervals[1].first]
		// ...
		// [inside_intervals[n-1].second,  v_max]
		//
		// for n = inside_intervals.size(). For each of these intervals, we mark the boxes on the left or right side of the interval as not inside

		// First test boundary curve crossing for each boundary parameter in u-direction
		for (int i = 0; i <= n_segs_u; ++i)
		  {
		    vector<pair<double, double> > inside_intervals;
		    try
		      {
			par_dom.getInsideIntervals(2, segment_pars[0][i], segment_pars[1][0], eps, inside_intervals);
		      }
		    catch (...)
		      {
		      }

		    // The number of outside intervals is one more than the number of inside intervals, therefore we iterate inside_intervals.size()+1 times
		    for (int j = 0; j <= (int)inside_intervals.size(); ++j)
		      {
			// The first and last segment in v-direction hitting the interval
			int outside_from = 0;
			int outside_to = n_segs_v - 1;

			// Adjust the first segment position according to inside_intervals[j-1].second
			if (j > 0)
			  {
			    double par = inside_intervals[j-1].second;
			    for (outside_from = 0; outside_from < n_segs_v - 1 && segment_pars[1][outside_from+1] < par; ++outside_from);
			  }

			// Adjust the last segment position according to inside_intervals[j].first
			if (j < (int)inside_intervals.size())
			  {
			    double par = inside_intervals[j].first;
			    for (outside_to = 0; outside_to < n_segs_v - 1 && segment_pars[1][outside_to+1] < par; ++outside_to);
			  }

			// Mark the relevant segments as not inside
			for (int k = outside_from; k <= outside_to; ++k)
			  {
			    if (i > 0)
			      structure->getBox(first_bb_idx + k*n_segs_u + (i - 1))->setInside(false);
			    if (i < n_segs_u)
			      structure->getBox(first_bb_idx + k*n_segs_u + i)->setInside(false);
			  }
		      }
		  }

		// Then test boundary curve crossing for each boundary parameter in v-direction in the same way as we did for the u-direction
		for (int i = 0; i <= n_segs_v; ++i)
		  {
		    vector<pair<double, double> > inside_intervals;
		    try
		      {
			par_dom.getInsideIntervals(1, segment_pars[0][0], segment_pars[1][i], eps, inside_intervals);
		      }
		    catch (...)
		      {
		      }

		    for (int j = 0; j <= (int)inside_intervals.size(); ++j)
		      {
			int outside_from = 0;
			int outside_to = n_segs_u - 1;

			if (j > 0)
			  {
			    double par = inside_intervals[j-1].second;
			    for (outside_from = 0; outside_from < n_segs_u - 1 && segment_pars[0][outside_from+1] < par; ++outside_from);
			  }

			if (j < (int)inside_intervals.size())
			  {
			    double par = inside_intervals[j].first;
			    for (outside_to = 0; outside_to < n_segs_u - 1 && segment_pars[0][outside_to+1] < par; ++outside_to);
			  }

			for (int k = outside_from; k <= outside_to; ++k)
			  {
			    if (i > 0)
			      structure->getBox(first_bb_idx + k + (i - 1)*n_segs_u)->setInside(false);
			    if (i < n_segs_v)
			      structure->getBox(first_bb_idx + k + i*n_segs_u)->setInside(false);
			  }
		      }
		  }

		//  ****************
		//   2 Find points on the surface near the noundary used to get an upperlimit on the distance from a point to the surface
		//  ****************

		// We create a regular grid in the parameter space first, the we add surface points
		// coming from a parameter pair in the grid such that the pair is both inside the
		// parameter domain, and either at the grid boundary, or next to a pair outside the
		// parameter domain

                if ((n_segs_u < 1) || (n_segs_v < 1))
                {
                    THROW("Method failed!");
                }

		double umin = segment_pars[0][0];
		double umax = segment_pars[0][n_segs_u];
		double vmin = segment_pars[1][0];
		double vmax = segment_pars[1][n_segs_v];

		int nmb_pts_each_dir = 20;
		double step_u = (umax - umin) / (double)nmb_pts_each_dir;
		double start_u = umin + 0.5 * step_u;
		double step_v = (vmax - vmin) / (double)nmb_pts_each_dir;
		double start_v = vmin + 0.5 * step_v;

		// Boolean grid telling whether each parameter pair is inside or not
		vector<vector<bool> > inside_grid(nmb_pts_each_dir);
		Array<double, 2> pars;
		for (int i = 0; i < nmb_pts_each_dir; ++i)
		  {
		    pars[0] = start_u + step_u*(double)i;
		    for (int j = 0; j < nmb_pts_each_dir; ++j)
		      {
			pars[1] = start_v + step_v*(double)j;
			inside_grid[i].push_back(par_dom.isInDomain(pars, eps));
		      }
		  }
		for (int i = 0; i < nmb_pts_each_dir; ++i)
		  for (int j = 0; j < nmb_pts_each_dir; ++j)
		    if (inside_grid[i][j])    // Pair is inside
		      {
			if (i == 0 || ! inside_grid[i-1][j] ||   // Neighbour to the left is outside
			    j == 0 || ! inside_grid[i][j-1] ||   // neighbour below is outside
			    i == nmb_pts_each_dir-1 || ! inside_grid[i+1][j] ||   // Neighbour to the right is outside
			    j == nmb_pts_each_dir-1 || ! inside_grid[i][j+1])     // Neighbour above is outside
			    surf_data->add_inside_point(underlying_paramSurf->point(start_u + step_u*(double)i, start_v + step_v*(double)j));
		      }

		//  ****************
		//   3 Get a polygon inside the parameter domain, and store polygon information on each segment
		//  ****************

		// We only look at the case where the only loop is the outside loop, i.e. the surface has
		// no holes.  Currently, we also restrict to the case where each curve patch in the loop
		// is either a spline curve or a line in the parameter space Also, we only add data for
		// segments where the loop does not enter the segment more than once
		vector<CurveLoop> loops = asBounded->allBoundaryLoops();
		if (loops.size() == 1)
		  {

		    // First we create the entire inside polygon by running through each curve in the loop and add points that, when connected,
		    // will give line segments always on the inside (or on) the boundary curve
		    vector<double> inside_ctrl_u, inside_ctrl_v;   // The coordinates of the corners in the polygon
		    bool polygon_found = true;     // Will be false if any of the curves in the loop could not be handled

		    vector<shared_ptr<ParamCurve> >::const_iterator it = loops[0].begin();
		    vector<shared_ptr<ParamCurve> >::const_iterator loops_end = loops[0].end();
		    for (; it != loops_end; ++it)
		      {
			shared_ptr<CurveOnSurface> cos = dynamic_pointer_cast<CurveOnSurface>(*it);
			if (cos.get())
			  {
			    shared_ptr<ParamCurve> par_crv = cos->parameterCurve();
			    shared_ptr<SplineCurve> par_spl = dynamic_pointer_cast<SplineCurve>(par_crv);
			    shared_ptr<Line> par_line = dynamic_pointer_cast<Line>(par_crv);

			    if (par_spl.get())
			      {
				// The parameter domain curve is a spline curve
				// For each point P on the control polygon of the spline curve, we use
				//   either   P if it is inside the parameter domain of the surface,
				//   or       the boundary curve evaluated in the Greville parameter of the B-spline corresponding to P if not
				const BsplineBasis bas = par_spl->basis();
				vector<double>::const_iterator it_coefs = par_spl->coefs_begin();
				vector<double>::const_iterator coefs_end = par_spl->coefs_end();
				int ctrl_pos = 0;

				// We iterate through the control polygon, but there is no need to add the last point, as this will be added
				// afterwards as the first point of the next curve patch in the loop
				for (; (it_coefs+2) != coefs_end; it_coefs += 2, ++ctrl_pos)
				  {
				    pars[0] = it_coefs[0];
				    pars[1] = it_coefs[1];
				    if (!par_dom.isInDomain(pars, eps))
				      {
					// Point is outside parameter domain, use evaluation at Greville parameter instead
					Point greville_pt;
					par_spl->point(greville_pt, bas.grevilleParameter(ctrl_pos));
					inside_ctrl_u.push_back(greville_pt[0]);
					inside_ctrl_v.push_back(greville_pt[1]);
				      }
				    else
				      {
					// Point is inside parameter domain
					inside_ctrl_u.push_back(pars[0]);
					inside_ctrl_v.push_back(pars[1]);
				      }
				  }
			      }

			    else if (par_line.get())
			      {
				// The parameter domain curve is a line. We only add the start point, the end point will be added afterwards as the
				// first point in the next curve patch in the loop
				Point start_pt;
				par_line->point(start_pt, par_line->startparam());
				inside_ctrl_u.push_back(start_pt[0]);
				inside_ctrl_v.push_back(start_pt[1]);
			      }

			    else
			      {
				// Not a spline or line
				polygon_found = false;
				break;
			      }
			  }

			else
			  {
			    // Not CurveOnSurface (should never happen)
			    polygon_found = false;
			    break;
			  }
		      }

		    if (!polygon_found) continue;

		    // Now we crate polygon information for each segment

		    // Information of the previous point during the polygon corners iteration
		    int prev_pos_u, prev_pos_v;    // Position (in segment grid) of the segment of the previous point
		    shared_ptr<SubSurfaceBoundingBox> prev_box;    // SubSurfaceBoundingBox object of the segment of the previous point
		    double prev_par_u, prev_par_v;  // Parameters for the previous point

		    // Information on segments that are 'banned' from holding a polygon. This happens if the global polygon enters a segment at least twice
		    bool prev_banned = false;      // Tells if the previous segment is 'banned'
		    vector<pair<int, int> > banned;       // List of all banned segments

		    // If the last point(s) and the first point(s) in the closed polygon are in the same segment, the last points
		    // must be added before the first. Therefore we need special treatment to temporarily store the first points until the
		    // last are added
		    bool is_in_first = true;           // Tells if we are still in the first segment
		    vector<double> first_u, first_v;    // Temporary storage of the points in the first segment

		    for (int i = 0; i <= (int)inside_ctrl_u.size(); ++i)
		      {
			// Get parameter of next point
			double par_u;
			double par_v;
			if (i == (int)inside_ctrl_u.size())
			  {
			    // We have reached the first point in the polygon again. It will not be stored as it has happened already,
			    // but we need it for the splitting process it in case the previous point was in another segment
			    par_u = inside_ctrl_u[0];
			    par_v = inside_ctrl_v[0];
			  }
			else
			  {
			    par_u = inside_ctrl_u[i];
			    par_v = inside_ctrl_v[i];
			  }

			// Get segment position of next point
			int pos_u, pos_v;
			for (pos_u = 0; pos_u < n_segs_u - 1 && segment_pars[0][pos_u+1] < par_u; ++pos_u);
			for (pos_v = 0; pos_v < n_segs_v - 1 && segment_pars[1][pos_v+1] < par_v; ++pos_v);

			if (i == 0)
			  {
			    // The first step in the iteration
			    first_u.push_back(par_u);
			    first_v.push_back(par_v);
			    prev_pos_u = pos_u;
			    prev_pos_v = pos_v;
			    prev_box = structure->getBox(first_bb_idx + prev_pos_u + prev_pos_v*n_segs_u);
			    prev_par_u = par_u;
			    prev_par_v = par_v;
			  }

			else
			  {
			    // Not the first step in the iteration, thus we can connect to the previous point to get a line segment to be handled

			    // Test if the line from the previous point must be split up into several segments
			    // For each time the line crosses a segment border, we handle the crossing, and repeat with the reamining line segment

			    // In the beginning of the loop below, 'shift_u' is true if the line segment
			    // crosses a segment border given by a fixed u-parameter.  Later in the loop,
			    // it will be changed to be true only if the fixed u-parameter border is the
			    // first border it crosses (i.e. before crossing a fixed v-parameter border,
			    // if that ever happens) 'shift_v' has the same function in the other
			    // parameter direction
			    bool shift_u = pos_u != prev_pos_u;
			    bool shift_v = pos_v != prev_pos_v;

			    while (shift_u || shift_v)
			      {
				// Line crosses segment border

				double shift_par;    // The line segment parameter (0.0 = start, 1.0 = end) where the first border corssing happens
				double shift_par_u;  // The u-parameter of the first border crossing
				double shift_par_v;  // The v-parameter of the first border crossing

				if (shift_u)
				  {
				    shift_par_u = segment_pars[0][prev_pos_u+((pos_u > prev_pos_u) ? 1 : 0)];
				    shift_par = (shift_par_u - prev_par_u) / (par_u - prev_par_u);
				  }

				if (shift_v)
				  {
				    shift_par_v = segment_pars[1][prev_pos_v+((pos_v > prev_pos_v) ? 1 : 0)];
				    double shift_par2 = (shift_par_v - prev_par_v) / (par_v - prev_par_v);

				    if (shift_u)
				      {
					// Boundary crossing might happen in both direction, determine which is first
					shift_u = shift_par < shift_par2;
					shift_v = !shift_u;
				      }
				    if (shift_v)
				      shift_par = shift_par2;
				  }

				// Get the first border crossing of the line and store it in the previous segment
				if (shift_u)
				  shift_par_v = prev_par_v + shift_par * (par_v - prev_par_v);
				else
				  shift_par_u = prev_par_u + shift_par * (par_u - prev_par_u);
				if (is_in_first)
				  {
				    first_u.push_back(shift_par_u);
				    first_v.push_back(shift_par_v);
				  }
				else if (!prev_banned)
				  prev_box->add_polygon_corners(shift_par_u, shift_par_v);

				// Update parameters and segment position at the shift
				if (shift_u)
				  prev_pos_u += (pos_u > prev_pos_u) ? 1 : -1;
				else
				  prev_pos_v += (pos_v > prev_pos_v) ? 1 : -1;
				prev_box = structure->getBox(first_bb_idx + prev_pos_u + prev_pos_v*n_segs_u);
				prev_par_u = shift_par_u;
				prev_par_v = shift_par_v;
				is_in_first = false;

				// Test if the next box is bannend, either because it was banned from before, or it had a plygon from before
				pair<int, int> prev_pair(prev_pos_u, prev_pos_v);
				prev_banned = find(banned.begin(), banned.end(), prev_pair) != banned.end();
				if (!prev_banned && prev_box->has_polygon())
				  {
				    prev_banned = true;
				    banned.push_back(prev_pair);
				    prev_box->remove_polygon();
				  }

				// Add the shift point to the next box
				if (!prev_banned)
				  prev_box->add_polygon_corners(shift_par_u, shift_par_v);

				shift_u = pos_u != prev_pos_u;
				shift_v = pos_v != prev_pos_v;

			      }  // End while-loop for segment boundary crossing

			    // Insert the new point if not banned, and if this is not the second time we are at the first point in the polygon
			    if (i < (int)inside_ctrl_u.size() && !prev_banned)
			      {
				if (is_in_first)
				  {
				    first_u.push_back(par_u);
				    first_v.push_back(par_v);
				  }
				else if (!prev_banned)
				  prev_box->add_polygon_corners(par_u, par_v);
			      }
			    prev_par_u = par_u;
			    prev_par_v = par_v;
			  }

		      }  // End running through the global polygon

		    // After running through the points, we now insert the first points at the end in their box
		    if (!prev_banned)
		      {
			for (int i = 0; i < (int)first_u.size(); ++i)
			  prev_box->add_polygon_corners(first_u[i], first_v[i]);

			// Make a closed loop if the boundary curve has been entirely inside one box
			if (is_in_first)
			  prev_box->add_polygon_corners(first_u[0], first_v[0]);
		      }
		  }
	      }

	    ++srf_idx;
	  }
      }

    // Make voxel structure
    structure->BuildVoxelStructure(bigbox, 1000.0);

#ifdef LOG_CLOSEST_POINTS
    cout << "Bounding boxes found = " << (structure->n_boxes()) << endl;
    int nx = structure->n_voxels_x();
    int ny = structure->n_voxels_y();
    int nz = structure->n_voxels_z();
    cout << "Number of voxels is " << nx << "*" << ny << "*" << nz << " = " << (nx*ny*nz) << endl;
    cout << "voxel_length = " << (structure->voxel_length()) << endl;

    clock_t t_after = clock();
    cout << endl << "Preprocessing timing = " << ((double)(t_after - t_before) / CLOCKS_PER_SEC) << " seconds" << endl;
#endif

    return structure;
  }


  namespace  // Anonymous
  {
    // Struct for data on possible candidate for closest point on a bounded surface, but with uncertainty about wether the point is inside the boundary.
    // This happens when closestPoint() has only been called on the underlying surface, and the preprocessing data are not precise enough to guarantee
    // is we are inside the boundary or not. This is done to postpone, and possibly cancel, a call to BoundedSurface::closestPoint() which is slower
    // than closestPoint() on the underlying surface.
    struct PossibleInside
    {
    public:
      // The closest point, as a point on the surface
      Point cl_p_;

      // The index of the surface in the surfaces model
      int surf_idx_;

      // The distance between the point and the closest point found
      double dist_;

      // The u-parameter of the closest point
      double par_u_;

      // The v-parameter of the closest point
      double par_v_;

      // An upper limit of the distance from the point to the surface in case the closest point is on the boundary.
      // This is given as the smallest distance between the point and a set of points (determined by the preprocessing) on the surface known to be
      // close to, but inside, the boundary
      double up_lim_boundary_;

      PossibleInside()
      { }

      PossibleInside(Point cl_p, int surf_idx, double dist, double u, double v, double up_lim_bd)
	: cl_p_(cl_p), surf_idx_(surf_idx), dist_(dist), par_u_(u), par_v_(v), up_lim_boundary_(up_lim_bd)
      {
      }
    };
  }


  void closestPointSingleCalculation(int pt_idx, int start_idx, int skip,
				     const vector<float>& inPoints,
				     const vector<vector<double> >& rotationMatrix, const Point& translation,
				     const shared_ptr<BoundingBoxStructure>& boxStructure,
				     vector<float>& result, vector<vector<int> >& lastBoxCall,
				     int return_type, int search_extend)
  {

    // Get transformed point
    int inPoints_idx = 3 * (start_idx + pt_idx * skip);
    Point pt(translation);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
	pt[i] += rotationMatrix[i][j] * inPoints[inPoints_idx + j];

    // Set thread id, used to avoid calling methods for the same surfaces from different threads when running in parallell
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif

    // *** Set some data related to the position of the point among the voxels ***

    // 1. Side length of voxels, and number of voxels in each direction
    double voxel_length = boxStructure->voxel_length();
    int nv_x = boxStructure->n_voxels_x();
    int nv_y = boxStructure->n_voxels_y();
    int nv_z = boxStructure->n_voxels_z();
    int nv_max = max(nv_x, max(nv_y, nv_z));

    // 2. Position of the voxel containing the point, and closest distance from the faces of that voxel to the point
    Point pt_vox_low = pt - boxStructure->big_vox_low();
    Point pt_rel = pt_vox_low / voxel_length;
    int n_x = (int)(pt_rel[0]);
    int n_y = (int)(pt_rel[1]);
    int n_z = (int)(pt_rel[2]);

    // 2.1 Closest distance in x-direction
    double shortest_face_distance = pt_vox_low[0] - (int)n_x * voxel_length;
    if (shortest_face_distance > 0.5 * voxel_length)
      shortest_face_distance = voxel_length - shortest_face_distance;

    // 2.2 Compare to closest distance in y-direction
    double next_face_distance = pt_vox_low[1] - (int)n_y * voxel_length;
    if (next_face_distance > 0.5 * voxel_length)
      next_face_distance = voxel_length - next_face_distance;
    if (shortest_face_distance > next_face_distance)
      shortest_face_distance = next_face_distance;

    // 2.3 Compare to closest distance in z-direction
    next_face_distance = pt_vox_low[2] - (int)n_z * voxel_length;
    if (next_face_distance > 0.5 * voxel_length)
      next_face_distance = voxel_length - next_face_distance;
    if (shortest_face_distance > next_face_distance)
      shortest_face_distance = next_face_distance;

    // 3. Smallest distance from the point to voxels at specific positions in each coordinate axis direction
    vector<double> pt_dist_x(nv_x);
    vector<double> pt_dist_y(nv_y);
    vector<double> pt_dist_z(nv_z);

    for (int i = 0; i < nv_x; ++i)
      {
	if (i < n_x)
	  pt_dist_x[i] = pt_vox_low[0] - (double)(i+1) * voxel_length;
	else if (i > n_x)
	  pt_dist_x[i] = (double)(i) * voxel_length - pt_vox_low[0];
	else
	  pt_dist_x[i] = 0.0;
      }
    for (int i = 0; i < nv_y; ++i)
      {
	if (i < n_y)
	  pt_dist_y[i] = pt_vox_low[1] - (double)(i+1) * voxel_length;
	else if (i > n_y)
	  pt_dist_y[i] = (double)(i) * voxel_length - pt_vox_low[1];
	else
	  pt_dist_y[i] = 0.0;
      }
    for (int i = 0; i < nv_z; ++i)
      {
	if (i < n_z)
	  pt_dist_z[i] = pt_vox_low[2] - (double)(i+1) * voxel_length;
	else if (i > n_z)
	  pt_dist_z[i] = (double)(i) * voxel_length - pt_vox_low[2];
	else
	  pt_dist_z[i] = 0.0;
      }

    // Candidates for closest point found by calling closestPoint() on underlying surface only, where it is still unclear whether the
    // closest point found lies inside the boundary
    vector<PossibleInside> poss_in;

    // Best point data
    bool any_clp_found = false;
    double best_dist;
    Point best_pt;
    double best_u;
    double best_v;
    int best_idx = -1;

    // Main loop for running through possible candidates.
    // If vox_span = -1, only check boxes (bounding boxes of surfaces segments) where the point lies inside the box
    // For vox_span >= 0, check unchecked boxes hitting voxels in position (i,j,k) where
    //                    abs(n_x-i), abs(n_y-j) and abs(n_z-j) are all <= vox_span
    for (int vox_span = -1; vox_span <= nv_max; ++vox_span)
      {

	if (vox_span > 0)
	  {
	    // All boxes hitting the voxel of the point have been checked.
	    // Test if any of the possible inside candidates are so close that they must be checked now by calling BoundedSurface::closestPoint()
	    double shortest_voxel_distance = (int)(vox_span-1) * voxel_length + shortest_face_distance;

	    int poss_in_size = (int)poss_in.size();
	    while(true)
	      {

		// Find the best candidate among those that are close enough to be checked (if any)
		// The best candidate is defined to be the one closest to the entire underlying surface (without boundaries)
		int best_poss_in = -1;
		for (int i = 0; i < poss_in_size; ++i)
		  if (poss_in[i].up_lim_boundary_ < shortest_voxel_distance || vox_span == nv_max)
		    {
		      if (best_poss_in == -1 || poss_in[i].dist_ < poss_in[best_poss_in].dist_)
			best_poss_in = i;
		    }
		if (best_poss_in == -1)
		  break;

		// A candidate has been found, remove it from list as it will be tested now
		PossibleInside p_i = poss_in[best_poss_in];
		--poss_in_size;
		if (best_poss_in < poss_in_size)
		  poss_in[best_poss_in] = poss_in[poss_in_size];
		poss_in.resize(poss_in_size);

		// Call BoundedeSurface::closestPoint() for the candidate
		double seed[2];
		seed[0] = p_i.par_u_;
		seed[1] = p_i.par_v_;
		double clo_u, clo_v;
		Point clo_pt;
		double clo_dist;
		shared_ptr<BoundedSurface> boundSurf = dynamic_pointer_cast<BoundedSurface>(boxStructure->getSurface(p_i.surf_idx_)->surface(thread_id));
		boundSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, NULL, &seed[0]);

		if (!any_clp_found || clo_dist < best_dist)
		  {
		    // The candidate is best closest point so far, update data about best closest point found.
		    best_dist = clo_dist;
		    best_idx = p_i.surf_idx_;
		    best_pt = clo_pt;
		    best_u = clo_u;
		    best_v = clo_v;
		    any_clp_found = true;

		    // Remove from list of possible candidates those that for sure are not better
		    for (int i = 0; i < poss_in_size;)
		      {
			if (poss_in[i].dist_ >= best_dist)
			  {
			    --poss_in_size;
			    if (i < poss_in_size)
			      poss_in[i] = poss_in[poss_in_size];
			    poss_in.resize(poss_in_size);
			  }
			else
			  ++i;
		      }
		  }
	      }

	    // Check if distance to best solution so far is smaller than shortest to the new voxels to be tested
	    if (poss_in_size == 0 && any_clp_found && best_dist < shortest_voxel_distance)
	      break;

	  }   // if (vox_span > 0)

	// Find range of voxels to be tested
	int vox_span_trunc = max(vox_span, 0);
	int beg_x = max(0, n_x - vox_span_trunc);
	int end_x = min(n_x + vox_span_trunc, nv_x - 1);
	int beg_y = max(0, n_y - vox_span_trunc);
	int end_y = min(n_y + vox_span_trunc, nv_y - 1);
	int beg_z = max(0, n_z - vox_span_trunc);
	int end_z = min(n_z + vox_span_trunc, nv_z - 1);

	// Run through voxels of interrest
	bool voxels_close = true;
	for (int vx = beg_x; vx <= end_x && voxels_close; ++vx)
	  {
	    double d2_x = pt_dist_x[vx] * pt_dist_x[vx];

	    for (int vy = beg_y; vy <= end_y && voxels_close; ++vy)
	      {
		double d2_xy = d2_x + pt_dist_y[vy] * pt_dist_y[vy];

		for (int vz = beg_z; vz <= end_z && voxels_close; ++vz)
		  {

		    // Avoid 'internal' voxels as all boxes hitting them are already checked
		    if (abs(vx - n_x) != vox_span_trunc &&
			abs(vy - n_y) != vox_span_trunc &&
			abs(vz - n_z) != vox_span_trunc)
		      continue;

		    // Test if this voxel is to far away
		    double d2_xyz = d2_xy + pt_dist_z[vz] * pt_dist_z[vz];
		    if (any_clp_found && d2_xyz > best_dist * best_dist)
		      continue;

		    // Get boxes to be tested
		    vector<int> possible_boxes;
		    if (vox_span >= 0)
		      possible_boxes = boxStructure->boxes_in_voxel(vx, vy, vz);
		    else
		      {
			// First iteration, only check with bounding boxes containing the point
			vector<int> voxel_boxes = boxStructure->boxes_in_voxel(vx, vy, vz);
			for (int i = 0; i < (int)voxel_boxes.size(); ++i)
			  {
			    int box_idx = voxel_boxes[i];
			    shared_ptr<SubSurfaceBoundingBox> surf_box = boxStructure->getBox(box_idx);
			    BoundingBox bb = boxStructure->getBox(box_idx)->box();
			    if (bb.containsPoint(pt))
			      possible_boxes.push_back(box_idx);
			  }
		      }

		    // Run through all boxes to be tested
		    for (int i = 0; i < (int)possible_boxes.size(); ++i)
		      {
			int box_idx = possible_boxes[i];
			if (lastBoxCall[thread_id][box_idx] == pt_idx)
			  continue;

			// Box has not been tested before, check if close enough
			shared_ptr<SubSurfaceBoundingBox> surf_box = boxStructure->getBox(box_idx);
			if (any_clp_found)
			  {
			    BoundingBox bb = surf_box->box();
			    Point low = bb.low();
			    Point high = bb.high();
			    double d2_pt_box = 0.0;
			    for (int j = 0; j < 3; ++j)
			      {
				double dist = max(0.0, max(pt[j] - high[j], low[j] - pt[j]));
				d2_pt_box += dist * dist;
			      }
			    if (d2_pt_box > best_dist * best_dist)
			      continue;
			  }

			// Box is close enough, run closest point, but only on underlying surface if main surface is BoundedSurface
			shared_ptr<SurfaceData> surf_data = surf_box->surface_data();
			shared_ptr<ParamSurface> paramSurf = surf_data->surface(thread_id);
			shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
			bool pt_might_be_outside = (boundedSurf.get() != NULL);
			if (pt_might_be_outside)
			  paramSurf = boundedSurf->underlyingSurface();

			int segs_u = surf_data->segs_u();
			int segs_v = surf_data->segs_v();

			// Set search domain in surface. Use the segment of the box, extended by 'search_extend' boxes in each direction

			int back_u = min(surf_box->pos_u(), search_extend);
			int back_v = min(surf_box->pos_v(), search_extend);

			int len_u = back_u + 1 + min(segs_u - (surf_box->pos_u() + 1), search_extend);
			int len_v = back_v + 1 + min(segs_v - (surf_box->pos_v() + 1), search_extend);

			int ll_index = box_idx - (back_v * segs_u + back_u);

			Array<double, 2> search_domain_ll, search_domain_ur;
			search_domain_ll[0] = boxStructure->getBox(ll_index)->par_domain()->umin();
			search_domain_ll[1] = boxStructure->getBox(ll_index)->par_domain()->vmin();
			search_domain_ur[0] = boxStructure->getBox(ll_index + len_u - 1)->par_domain()->umax();
			search_domain_ur[1] = boxStructure->getBox(ll_index + (len_v - 1)*segs_u)->par_domain()->vmax();
			shared_ptr<RectDomain> search_domain(new RectDomain(search_domain_ll, search_domain_ur));

			// Set other input variables and call closestPoint()
			shared_ptr<RectDomain> rd = surf_box->par_domain();
			double seed[2];
			seed[0] = (rd->umin() + rd->umax()) * 0.5;
			seed[1] = (rd->vmin() + rd->vmax()) * 0.5;

			double clo_u, clo_v;
			Point clo_pt;
			double clo_dist;

			paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, search_domain.get(), &seed[0]);

			for (int j = 0; j < len_u; ++j)
			  for (int k = 0; k < len_v; ++k)
			    lastBoxCall[thread_id][ll_index + k * segs_u + j] = pt_idx;

			// If top surface is BoundedSurface, check if this point might be outside
			if (pt_might_be_outside)
			  {
			    int pos_u, pos_v;  // Position of box holding closest point, truncated to search domain

			    if (clo_u <= search_domain_ll[0])
			      pos_u = back_u;
			    else if (clo_u >= search_domain_ur[0])
			      pos_u = back_u + len_u - 1;
			    else
			      {
				for (pos_u = back_u;
				     pos_u < back_u + len_u - 1 &&
				       boxStructure->getBox(ll_index + pos_u - back_u)->par_domain()->umax() < clo_u;
				     ++pos_u);
			      }

			    if (clo_v <= search_domain_ll[1])
			      pos_v = back_v;
			    else if (clo_v >= search_domain_ur[1])
			      pos_v = back_v + len_v - 1;
			    else
			      {
				for (pos_v = back_v;
				     pos_v < back_v + len_v - 1 &&
				       boxStructure->getBox(ll_index + (pos_v - back_v)*segs_u)->par_domain()->vmax() < clo_v;
				     ++pos_v);
			      }

			    int cl_p_box = ll_index + (pos_v - back_v)*segs_u + pos_u - back_u;
			    pt_might_be_outside = !(boxStructure->getBox(cl_p_box)->inside(clo_u, clo_v));
			  }

			if (!any_clp_found || clo_dist < best_dist)
			  {
			    // Point is close enough to be a candidate for closest point

			    if (pt_might_be_outside)
			      {
				// The point might be outside the parameter domain. Store it as a case we might have to handle later
				// First check if this point has been found before
				int surf_idx = surf_data->index();
				double tol = 1.0e-4;
				bool insert = true;
				for (int j = 0; j < (int)poss_in.size() && insert; ++j)
				  insert = surf_idx != poss_in[j].surf_idx_ ||
				    abs(clo_u - poss_in[j].par_u_) > tol ||
				    abs(clo_v - poss_in[j].par_v_) > tol;

				// Point is not found before, insert it
				if (insert)
				  {
				    double up_lim_b2 = voxel_length * voxel_length * (double)(nv_x*nv_x + nv_y*nv_y + nv_z*nv_z);
				    vector<Point> surf_pts = surf_data->inside_points();
				    for (int j = 0; j < (int)surf_pts.size(); ++j)
				      {
					double dist2 = pt.dist2(surf_pts[j]);
					if (dist2 < up_lim_b2)
					  up_lim_b2 = dist2;
				      }

				    poss_in.push_back(PossibleInside(clo_pt, surf_idx, clo_dist, clo_u, clo_v, sqrt(up_lim_b2)));
				  }
			      }

			    else
			      {
				// The point is inside the parameter domain and the closest point found so far
				// Update information about closest point, and remove possible inside candidates that are too far away

				best_dist = clo_dist;
				best_idx = surf_data->index();
				best_pt = clo_pt;
				best_u = clo_u;
				best_v = clo_v;
				any_clp_found = true;
				int poss_in_size = (int)poss_in.size();
				for (int j = 0; j < poss_in_size;)
				  {
				    if (poss_in[j].dist_ >= best_dist)
				      {
					--poss_in_size;
					if (j < poss_in_size)
					  poss_in[j] = poss_in[poss_in_size];
					poss_in.resize(poss_in_size);
				      }
				    else
				      ++j;
				  }
			      }
			  }  // End 'Point is close enough to be a candidate for closest point'
		      }  // End running through all boxes for this voxel
		  }  // End voxels in z-dir
	      }  // End voxels in y-dir
	  }  // End voxels in x-dir
      }  // End vox_span loop

    if (return_type == 0)  // Store distance
      result[pt_idx] = (float)best_dist;
    else if (return_type == 1) // Store signed distance
      {
	pt_idx;
	shared_ptr<ParamSurface> paramSurf = boxStructure->getSurface(best_idx)->surface(thread_id);
	shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
	if (boundedSurf.get())
	  paramSurf = boundedSurf->underlyingSurface();
	Point normal;
	paramSurf->normal(normal, best_u, best_v);
	double best_dist_factor = (normal * (pt - best_pt) >= 0.0) ? 1.0 : -1.0;
	result[pt_idx] = (float)(best_dist_factor * best_dist);
      }
    else if (return_type == 2) // Store closest point
      {
	for (int i = 0; i < 3; ++i)
	  result[3 * pt_idx + i] = (float)best_pt[i];
      }
    else if (return_type == 3) // Store signed dist, surface index, clo_u, clo_v.
    {
	pt_idx;
	shared_ptr<ParamSurface> paramSurf = boxStructure->getSurface(best_idx)->surface(thread_id);
	shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
	if (boundedSurf.get())
	  paramSurf = boundedSurf->underlyingSurface();
	Point normal;
	paramSurf->normal(normal, best_u, best_v);
	double best_dist_factor = (normal * (pt - best_pt) >= 0.0) ? 1.0 : -1.0;
        result[4*pt_idx] = (float)(best_dist_factor * best_dist);
        result[4*pt_idx + 1] = (float)best_idx;
        result[4*pt_idx + 2] = (float)best_u;
        result[4*pt_idx + 3] = (float)best_v;
    }
  }


  vector<float> closestPointCalculations(const vector<float>& inPoints, const shared_ptr<BoundingBoxStructure>& boxStructure,
 					 const vector<vector<double> >& rotationMatrix, const Point& translation,
					 int return_type, int start_idx, int skip, int max_idx, int search_extend, bool m_core)
  {
#ifdef LOG_CLOSEST_POINTS
    clock_t t_before = clock();
    double time_factor = 1.0;
#endif

    max_idx = min(max_idx, (int)inPoints.size() / 3);

    int nmb_points_tested = (max_idx + skip - start_idx - 1) / skip;
    if (nmb_points_tested < 0)
      nmb_points_tested = 0;
    int result_size = nmb_points_tested;
    if (return_type == 2)
      result_size *= 3;
    else if (return_type == 3)
      result_size *= 4;
    vector<float> result(result_size);

    // if (return_type == 2)
    //   result.resize(nmb_points_tested * 3);
    // else if (return_type == 3)
    //   result.resize(nmb_points_tested * 4);
    // else
    //   result.resize(nmb_points_tested);

#ifdef _OPENMP
    int max_threads = omp_get_max_threads();
#else
    int max_threads = 1;
#endif

    boxStructure->setSurfaceCopies(max_threads);
    vector<vector<int> > lastBoxCall(max_threads);
    for (int i = 0; i < max_threads; ++i)
      lastBoxCall[i].resize(boxStructure->n_boxes(), -1);

#ifdef _OPENMP
    if (m_core)
      {
	// Run all closest point calculations in multicore, because m_core=true and OPENMP is included

#ifdef LOG_CLOSEST_POINTS
	time_factor = 1.0 / (double)max_threads;
#endif

	int pt_idx;
#pragma omp parallel \
  default(none)	\
  private(pt_idx) \
  shared(nmb_points_tested, start_idx, skip, inPoints, rotationMatrix, translation, boxStructure, result, lastBoxCall, return_type, search_extend)
#pragma omp for OMP_SCHEDULE_AUTO
	for (pt_idx = 0; pt_idx < nmb_points_tested; ++pt_idx)
	  closestPointSingleCalculation(pt_idx, start_idx, skip, inPoints, rotationMatrix, translation, boxStructure,
					result, lastBoxCall, return_type, search_extend);
      }

    else

      {
	// Run all closest point calculations in one single thread, because m_core=false
	for (int pt_idx = 0; pt_idx < nmb_points_tested; ++pt_idx)
	  closestPointSingleCalculation(pt_idx, start_idx, skip, inPoints, rotationMatrix, translation, boxStructure,
					result, lastBoxCall, return_type, search_extend);
      }

#else   // #ifdef _OPENMP

    // Run all closest point calculations in one single thread, because OPENMP is not included
    for (int pt_idx = 0; pt_idx < nmb_points_tested; ++pt_idx)
      closestPointSingleCalculation(pt_idx, start_idx, skip, inPoints, rotationMatrix, translation, boxStructure,
				    result, lastBoxCall, return_type, search_extend);
#endif   // #ifdef _OPENMP

#ifdef LOG_CLOSEST_POINTS
    clock_t t_after = clock();
    cout << endl << "Closest point timing = " << ((double)(t_after - t_before) * time_factor / CLOCKS_PER_SEC) << " seconds" << endl;

    cout << "Number of points tested for closestPoint = " << nmb_points_tested << endl;

    double sum_d2 = 0.0;
    double best_d2 = 100000000.0;
    double worst_d2 = -1.0;
    int cnt_d2 = 0;
    for (int idx = 0; idx < (int)result.size(); ++idx)
      {
	double d2 = result[idx] * result[idx];
	if (best_d2 > d2)
	  best_d2 = d2;
	if (worst_d2 < d2)
	  worst_d2 = d2;
	sum_d2 += d2;
	++cnt_d2;
      }
    cout << endl;
    cout << "Average square distance = " << (sum_d2 / (double)cnt_d2) << "  sqrt = " << sqrt(sum_d2 / (double)cnt_d2) << endl;
    cout << "Best square distance = " << best_d2 << "  sqrt = " << sqrt(best_d2) << endl;
    cout << "Worst square distance = " << worst_d2 << "  sqrt = " << sqrt(worst_d2) << endl;
#endif
    return result;
  }


  vector<float> closestVectorsOld(const vector<float>& inPoints, const shared_ptr<BoundingBoxStructure>& boxStructure,
				  const vector<vector<double> >& rotationMatrix, const Point& translation,
				  int test_type, int start_idx, int skip, int max_idx, int search_extend)
  {
    clock_t t_before = clock();
    vector<float> result;
    double skip_cnt = -1;

    double voxel_length = boxStructure->voxel_length();
    int nv_x = boxStructure->n_voxels_x();
    int nv_y = boxStructure->n_voxels_y();
    int nv_z = boxStructure->n_voxels_z();
    int nv_max = max(nv_x, max(nv_y, nv_z));
    vector<double> pt_dist_x(nv_x);
    vector<double> pt_dist_y(nv_y);
    vector<double> pt_dist_z(nv_z);

    int total_pts_tested = 0;
    vector<int> isBest(boxStructure->n_surfaces(), 0);
    vector<int> boundaryCalls;
    vector<int> underlyingCalls;
    vector<int> lastBoxCall(boxStructure->n_boxes(), -1);
    for (int idx = 0, pt_idx = 0; idx < (int)inPoints.size(); idx += 3, ++pt_idx)
      {
	if (pt_idx == start_idx)
	  skip_cnt = 0;
	if (pt_idx >= max_idx)
	  skip_cnt = -1;
	if (skip_cnt == 0)
	  {
	    ++total_pts_tested;
	    Point pt(translation);
	    for (int i = 0; i < 3; ++i)
	      for (int j = 0; j < 3; ++j)
		pt[i] += rotationMatrix[i][j] * inPoints[idx + j];
	    double best_dist = 0.0;
	    Point best_dist_pt;
	    int best_idx = -1;
	    int local_boundaryCalls = 0;
	    int local_underlyingCalls = 0;

	    Point pt_vox_low = pt - boxStructure->big_vox_low();
	    Point pt_rel = pt_vox_low / voxel_length;
	    int n_x = (int)(pt_rel[0]);
	    int n_y = (int)(pt_rel[1]);
	    int n_z = (int)(pt_rel[2]);

	    double shortest_face_distance = pt_vox_low[0] - (int)n_x * voxel_length;
	    if (shortest_face_distance > 0.5 * voxel_length)
	      shortest_face_distance = voxel_length - shortest_face_distance;
	    double next_face_distance = pt_vox_low[1] - (int)n_y * voxel_length;
	    if (next_face_distance > 0.5 * voxel_length)
	      next_face_distance = voxel_length - shortest_face_distance;
	    if (shortest_face_distance > next_face_distance)
	      shortest_face_distance = next_face_distance;
	    next_face_distance = pt_vox_low[2] - (int)n_z * voxel_length;
	    if (next_face_distance > 0.5 * voxel_length)
	      next_face_distance = voxel_length - shortest_face_distance;
	    if (shortest_face_distance > next_face_distance)
	      shortest_face_distance = next_face_distance;

	    bool any_tested = false;

	    // First test all bounding boxes containing the point
	    if (n_x >= 0 && n_x < nv_x &&
		n_y >= 0 && n_y < nv_y &&
		n_z >= 0 && n_z < nv_z)
	      {
		vector<int> possible_boxes = boxStructure->boxes_in_voxel(n_x, n_y, n_z);

		vector<pair<int, double> > inside_boxes;   // All boxes (first) sorted by distance from point to most distant corner (second)
		for (int i = 0; i < (int)possible_boxes.size(); ++i)
		  {
		    int box_idx = possible_boxes[i];
		    shared_ptr<SubSurfaceBoundingBox> surf_box = boxStructure->getBox(box_idx);
		    BoundingBox bb = surf_box->box();
		    if (bb.containsPoint(pt))
		      {
			Point l = bb.low();
			Point h = bb.high();
			double dist2 = 0.0;
			for (int j = 0; j < 3; ++j)
			  {
			    double d = max(pt[j] - l[j], h[j] - pt[j]);
			    dist2 += d*d;
			  }
			inside_boxes.push_back(pair<int, double>(box_idx, dist2));
		      }
		  }
		int nmb_inside_boxes = (int)inside_boxes.size();

		// Inside box list should be small, bubble sort should be good enough
		for (int i = 0; i < nmb_inside_boxes - 1; ++i)
		  for (int j = nmb_inside_boxes - 1; j > i; --j)
		    if (inside_boxes[j].second < inside_boxes[j-1].second)
		      swap(inside_boxes[j], inside_boxes[j-1]);

		for (int turn = 0; turn < 2; ++turn)
		  for (int i = 0; i < nmb_inside_boxes; ++i)
		    {
		      int box_idx = inside_boxes[i].first;
		      if (lastBoxCall[box_idx] < pt_idx)
			{
			  shared_ptr<SubSurfaceBoundingBox> surf_box = boxStructure->getBox(box_idx);
			  shared_ptr<SurfaceData> surf_data = surf_box->surface_data();
			  shared_ptr<ParamSurface> paramSurf = surf_data->surface(0);
			  int segs_u = surf_data->segs_u();
			  int segs_v = surf_data->segs_v();

			  int back_u = min(surf_box->pos_u(), search_extend);
			  int back_v = min(surf_box->pos_v(), search_extend);
			  int len_u = back_u + 1 + min(segs_u - (surf_box->pos_u() + 1), search_extend);
			  int len_v = back_v + 1 + min(segs_v - (surf_box->pos_v() + 1), search_extend);
			  int ll_index = box_idx - (back_v * segs_u + back_u);

			  // Here
			  Array<double, 2> search_domain_ll, search_domain_ur;
			  search_domain_ll[0] = boxStructure->getBox(ll_index)->par_domain()->umin();
			  search_domain_ll[1] = boxStructure->getBox(ll_index)->par_domain()->vmin();
			  search_domain_ur[0] = boxStructure->getBox(ll_index + len_u - 1)->par_domain()->umax();
			  search_domain_ur[1] = boxStructure->getBox(ll_index + (len_v - 1)*segs_u)->par_domain()->vmax();
			  // TO here
			  shared_ptr<RectDomain> search_domain(new RectDomain(search_domain_ll, search_domain_ur));
			  bool isInside = true;
			  for (int j = 0; j < len_v && isInside; ++j)
			    for (int k = 0; k < len_u && isInside; ++k)
			      isInside &= boxStructure->getBox(ll_index + j * segs_u + k)->inside();
			  if (isInside == (turn == 1))
			    continue;

			  shared_ptr<RectDomain> rd = surf_box->par_domain();
			  double seed[2];
			  seed[0] = (rd->umin() + rd->umax()) * 0.5;
			  seed[1] = (rd->vmin() + rd->vmax()) * 0.5;

			  shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
			  double clo_u, clo_v;
			  Point clo_pt;
			  double clo_dist;

			  if (boxStructure->closestPoint(box_idx, any_tested, best_dist, isInside, pt,
							 clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, search_domain.get(), &seed[0]))
			    {
			      if (isInside)
				++local_underlyingCalls;
			      else
				++local_boundaryCalls;
			      if (!any_tested || clo_dist < best_dist)
				{
				  best_dist = clo_dist;
				  best_idx = surf_data->index();
				  best_dist_pt = clo_pt - pt;
				  any_tested = true;
				}
			    }

			  for (int j = 0; j < len_u; ++j)
			    for (int k = 0; k < len_v; ++k)
			      lastBoxCall[ll_index + k * segs_u + j] = pt_idx;
			}
		    }
	      }

	    // Then test rest in this and other voxels
	    for (int i = 0; i < nv_x; ++i)
	      {
		if (i < n_x)
		  pt_dist_x[i] = pt_vox_low[0] - (double)(i+1) * voxel_length;
		else if (i > n_x)
		  pt_dist_x[i] = (double)(i) * voxel_length - pt_vox_low[0];
		else
		  pt_dist_x[i] = 0.0;
	      }
	    for (int i = 0; i < nv_y; ++i)
	      {
		if (i < n_y)
		  pt_dist_y[i] = pt_vox_low[1] - (double)(i+1) * voxel_length;
		else if (i > n_y)
		  pt_dist_y[i] = (double)(i) * voxel_length - pt_vox_low[1];
		else
		  pt_dist_y[i] = 0.0;
	      }
	    for (int i = 0; i < nv_z; ++i)
	      {
		if (i < n_z)
		  pt_dist_z[i] = pt_vox_low[2] - (double)(i+1) * voxel_length;
		else if (i > n_z)
		  pt_dist_z[i] = (double)(i) * voxel_length - pt_vox_low[2];
		else
		  pt_dist_z[i] = 0.0;
	      }
	    for (int vox_span = 0; vox_span < nv_max; ++vox_span)
	      {
		double shortest_voxel_distance = (int)(vox_span-1) * voxel_length + shortest_face_distance;
		if (shortest_voxel_distance < 0.0)
		  shortest_voxel_distance = 0.0;

		if (any_tested && best_dist < shortest_voxel_distance)
		  break;

		int beg_x = n_x - vox_span;
		if (beg_x < 0)
		  beg_x = 0;
		int end_x = n_x + vox_span;
		if (end_x >= nv_x)
		  end_x = nv_x - 1;
		int beg_y = n_y - vox_span;
		if (beg_y < 0)
		  beg_y = 0;
		int end_y = n_y + vox_span;
		if (end_y >= nv_y)
		  end_y = nv_y - 1;
		int beg_z = n_z - vox_span;
		if (beg_z < 0)
		  beg_z = 0;
		int end_z = n_z + vox_span;
		if (end_z >= nv_z)
		  end_z = nv_z - 1;
		bool voxels_close = true;
		for (int vx = beg_x; vx <= end_x && voxels_close; ++vx)
		  {
		    double d2_x = pt_dist_x[vx] * pt_dist_x[vx];
		    for (int vy = beg_y; vy <= end_y && voxels_close; ++vy)
		      {
			double d2_xy = d2_x + pt_dist_y[vy] * pt_dist_y[vy];
			for (int vz = beg_z; vz <= end_z && voxels_close; ++vz)
			  {
			    if (abs(vx - n_x) != vox_span &&
				abs(vy - n_y) != vox_span &&
				abs(vz - n_z) != vox_span)
			      continue;
			    double d2_xyz = d2_xy + pt_dist_z[vz] * pt_dist_z[vz];
			    if (d2_xyz > best_dist * best_dist)
			      continue;
			    vector<int> possible_boxes = boxStructure->boxes_in_voxel(vx, vy, vz);
			    for (int i = 0; i < (int)possible_boxes.size(); ++i)
			      {
				int box_idx = possible_boxes[i];
				shared_ptr<SubSurfaceBoundingBox> surf_box = boxStructure->getBox(box_idx);
				if (any_tested)
				  {
				    if (lastBoxCall[box_idx] == pt_idx)
				      continue;
				    BoundingBox bb = surf_box->box();
				    Point low = bb.low();
				    Point high = bb.high();
				    double d2_pt_box = 0.0;
				    for (int j = 0; j < 3; ++j)
				      {
					double dist = (pt[j] < low[j]) ? (low[j] - pt[j]) : ((pt[j] > high[j]) ? (high[j] - pt[j]) : 0.0);
					d2_pt_box += dist * dist;
				      }
				    if (d2_pt_box > best_dist * best_dist)
				      continue;
				  }

				// New box not tested before
				shared_ptr<SurfaceData> surf_data = surf_box->surface_data();
				shared_ptr<ParamSurface> paramSurf = surf_data->surface(0);
				shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
				double clo_u, clo_v;
				Point clo_pt;
				double clo_dist;
				shared_ptr<RectDomain> rd = surf_box->par_domain();
				double seed[2];
				seed[0] = (rd->umin() + rd->umax()) * 0.5;
				seed[1] = (rd->vmin() + rd->vmax()) * 0.5;

				if (boxStructure->closestPoint(box_idx, any_tested, best_dist, surf_box->inside(), pt,
							       clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, rd.get(), &seed[0]))
				  {
				    if (surf_box->inside())
				      ++local_underlyingCalls;
				    else
				      ++local_boundaryCalls;
				    if (!any_tested || clo_dist < best_dist)
				      {
					best_dist = clo_dist;
					best_idx = surf_data->index();
					best_dist_pt = clo_pt - pt;
					any_tested = true;
				      }
				  }

				lastBoxCall[box_idx] = pt_idx;
			      }
			  }
		      }
		  }
	      }
	    ++isBest[best_idx];
	    result.push_back((float)best_dist);

	    if (local_boundaryCalls >= (int)boundaryCalls.size())
	      boundaryCalls.resize(local_boundaryCalls + 1);
	    ++boundaryCalls[local_boundaryCalls];
	    if (local_underlyingCalls >= (int)underlyingCalls.size())
	      underlyingCalls.resize(local_underlyingCalls + 1);
	    ++underlyingCalls[local_underlyingCalls];
	  }

	if (skip_cnt != -1)
	  {
	    ++skip_cnt;
	    if (skip_cnt == skip)
	      skip_cnt = 0;
	  }
      }

    clock_t t_after = clock();
    cout << endl << "Closest point timing = " << ((double)(t_after - t_before) / CLOCKS_PER_SEC) << " seconds" << endl;

    cout << "Number of points tested for closestPoint = " << total_pts_tested << endl;

    cout << "Surf\tBest\tUnder" << endl;
    for (int i = 0; i < boxStructure->n_surfaces(); ++i)
      {
	string type_str = "Unknown";
	shared_ptr<ParamSurface> paramSurf = boxStructure->getSurface(i)->surface(0);
	shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
	if (boundedSurf.get())
	  paramSurf = boundedSurf->underlyingSurface();
	if (dynamic_pointer_cast<SplineSurface>(paramSurf).get())
	  type_str  = "Spline";
	else if (dynamic_pointer_cast<Sphere>(paramSurf).get())
	    type_str = "Sphere";
	else if (dynamic_pointer_cast<Cylinder>(paramSurf).get())
	    type_str = "Cylinder";
	else if (dynamic_pointer_cast<Plane>(paramSurf).get())
	    type_str = "Plane";
	cout << i << "\t" << isBest[i] << "\t" << type_str << endl;
      }

    cout << endl << "N\tBound\tNot b" << endl;
    int tot_bs = 0;
    int tot_ul = 0;
    for (int i = 0; i < (int)boundaryCalls.size() || i < (int)underlyingCalls.size(); ++i)
      {
	int bs_calls = 0;
	if (i < (int)boundaryCalls.size())
	  bs_calls = boundaryCalls[i];
	int ul_calls = 0;
	if (i < (int)underlyingCalls.size())
	  ul_calls = underlyingCalls[i];
	int tot_calls = bs_calls + ul_calls;
	if (i < 6)
	  {
	    if (tot_calls > 0)
	      {
		cout << i << "\t";
		if (bs_calls > 0)
		  cout << bs_calls;
		cout << "\t";
		if (ul_calls > 0)
		  cout << ul_calls;
		cout << endl;
	      }
	  }
	else if (i == 6)
	  cout << "..." << endl;
	tot_bs += i * bs_calls;
	tot_ul += i * ul_calls;
      }
    cout << "Total\t" << tot_bs << "\t" << tot_ul << endl;

    double sum_d2 = 0.0;
    double best_d2 = 100000000.0;
    double worst_d2 = -1.0;
    int cnt_d2 = 0;
    for (int idx = 0; idx < (int)result.size(); ++idx)
      {
	double d2 = result[idx] * result[idx];
	if (best_d2 > d2)
	  best_d2 = d2;
	if (worst_d2 < d2)
	  worst_d2 = d2;
	sum_d2 += d2;
	++cnt_d2;
      }
    cout << endl;
    cout << "Average square distance = " << (sum_d2 / (double)cnt_d2) << "  sqrt = " << sqrt(sum_d2 / (double)cnt_d2) << endl;
    cout << "Best square distance = " << best_d2 << "  sqrt = " << sqrt(best_d2) << endl;
    cout << "Worst square distance = " << worst_d2 << "  sqrt = " << sqrt(worst_d2) << endl;
    return result;
  }



  vector<float> closestPointCalculations(const vector<float>& inPoints, const shared_ptr<BoundingBoxStructure>& boxStructure,
					 const vector<vector<double> >& rotationMatrix, const Point& translation, int return_type)
  {
    int nmb_pts = ((int)inPoints.size()) / 3;
    return closestPointCalculations(inPoints, boxStructure, rotationMatrix, translation, return_type, 0, 1, nmb_pts);
  }

  vector<float> closestDistances(const vector<float>& inPoints, const shared_ptr<BoundingBoxStructure>& boxStructure,
				 const vector<vector<double> >& rotationMatrix, const Point& translation)
  {
    return closestPointCalculations(inPoints, boxStructure, rotationMatrix, translation, 0);
  }

  vector<float> closestSignedDistances(const vector<float>& inPoints, const shared_ptr<BoundingBoxStructure>& boxStructure,
				       const vector<vector<double> >& rotationMatrix, const Point& translation)
  {
    return closestPointCalculations(inPoints, boxStructure, rotationMatrix, translation, 1);
  }

  vector<float> closestPoints(const vector<float>& inPoints, const shared_ptr<BoundingBoxStructure>& boxStructure,
			      const vector<vector<double> >& rotationMatrix, const Point& translation)
  {
    return closestPointCalculations(inPoints, boxStructure, rotationMatrix, translation, 2);
  }

  vector<float> closestSignedDistanceSfParams(const vector<float>& inPoints, const shared_ptr<BoundingBoxStructure>& boxStructure,
                                              const vector<vector<double> >& rotationMatrix, const Point& translation)
  {
    return closestPointCalculations(inPoints, boxStructure, rotationMatrix, translation, 3);
  }

}   // end namespace Go
