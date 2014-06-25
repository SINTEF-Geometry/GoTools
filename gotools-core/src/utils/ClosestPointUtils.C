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
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/utils/ClosestPointUtils.h"

using namespace std;
using namespace Go;
using namespace Go::boxStructuring;


namespace Go
{

  shared_ptr<BoundingBoxStructure> preProcessClosestVectors(const vector<shared_ptr<GeomObject> >& surfaces, double par_len_el)
  {
    clock_t t_before = clock();

    shared_ptr<BoundingBoxStructure> structure(new BoundingBoxStructure());
    BoundingBox bigbox(3);
    vector<shared_ptr<GeomObject> >::const_iterator surf_end = surfaces.end();
    int srf_idx = 0;

    for (vector<shared_ptr<GeomObject> >::const_iterator surf_it = surfaces.begin(); surf_it != surf_end; ++surf_it)
      {
	shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(*surf_it);
	if (paramSurf.get())
	  {
	    shared_ptr<SurfaceData> surf_data(new SurfaceData(paramSurf));
	    structure->addSurface(surf_data);

	    shared_ptr<BoundedSurface> asBounded = dynamic_pointer_cast<BoundedSurface>(paramSurf);
	    if (asBounded.get())
	      paramSurf = asBounded->underlyingSurface();
	    shared_ptr<SplineSurface> splineSurf = dynamic_pointer_cast<SplineSurface>(paramSurf);
	    shared_ptr<ElementarySurface> elSurf = dynamic_pointer_cast<ElementarySurface>(paramSurf);
	    vector<vector<double> > segment_pars(2);
	    if (splineSurf.get())
	      {
		shared_ptr<SplineSurface> surf_copy(splineSurf->clone());

		// Make C0 to get Bezier data
		vector<vector<int> > seg_ctrl_pos(2);
		int order_u = surf_copy->order_u();
		int order_v = surf_copy->order_v();
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

		// Fetch segment bounding boxes
		int n_segs_u = seg_ctrl_pos[0].size();
		int n_segs_v = seg_ctrl_pos[1].size();
		surf_data->setSegments(n_segs_u, n_segs_v);
		int n_coefs_u = surf_copy->numCoefs_u();

		vector<double>::const_iterator ctrl_it_v = surf_copy->ctrl_begin();
		for (int j = 0; j < n_segs_v; ++j)
		  {
		    double start_v = segment_pars[1][j];
		    double end_v = segment_pars[1][j+1];
		    if (j > 0)
		      ctrl_it_v += (seg_ctrl_pos[1][j] - seg_ctrl_pos[1][j-1]) * n_coefs_u * 3;
		    vector<double>::const_iterator ctrl_it_u = ctrl_it_v;
		    for (int i = 0; i < n_segs_u; ++i)
		      {
			if (i > 0)
			  ctrl_it_u += (seg_ctrl_pos[0][i] - seg_ctrl_pos[0][i-1]) * 3;
			vector<double> ctrl_pts;
			for (int ctrl_pos = 0, k = 0; k < order_v; ++k, ctrl_pos += (n_coefs_u - order_u) * 3)
			  for (int l = 0; l < order_u; ++l)
			    for (int m = 0; m < 3; ++m, ++ctrl_pos)
			      ctrl_pts.push_back(ctrl_it_u[ctrl_pos]);

			BoundingBox bb;
			bb.setFromArray(ctrl_pts.begin(), ctrl_pts.end(), 3);
			Array<double, 2> ll, ur;
			ll[0] = segment_pars[0][i];
			ll[1] = start_v;
			ur[0] = segment_pars[0][i];
			ur[1] = end_v;
			shared_ptr<SubSurfaceBoundingBox> box(new SubSurfaceBoundingBox(surf_data, i, j, bb, shared_ptr<RectDomain>(new RectDomain(ll, ur))));
			structure->addBox(box);

			bigbox.addUnionWith(bb);
		      }
		  }
	      }

	    if (elSurf.get())
	      {
		shared_ptr<Sphere> sphSurf = dynamic_pointer_cast<Sphere>(elSurf);
		shared_ptr<Cylinder> cylSurf = dynamic_pointer_cast<Cylinder>(elSurf);
		RectDomain big_rd = elSurf->containingDomain();

		double umin = big_rd.umin();
		double umax = big_rd.umax();
		double vmin = big_rd.vmin();
		double vmax = big_rd.vmax();
		if (elSurf->isSwapped())
		  {
		    swap(umin, vmin);
		    swap(umax, vmax);
		  }

		double len_u = 0.0;
		double len_v = 0.0;
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
		if (elSurf->isSwapped())
		  {
		    swap(umin, vmin);
		    swap(umax, vmax);
		    swap(len_u, len_v);
		  }

		int n_segs_u = (int)(0.5 + len_u / par_len_el);
		int n_segs_v = (int)(0.5 + len_v / par_len_el);
		if (n_segs_u == 0)
		  n_segs_u = 1;
		if (n_segs_v == 0)
		  n_segs_v = 1;
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

		// Find all segments
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

	    // Determine if the parameter domains of the bounding boxes are entirely inside the parameter domain limited by the boundary curves
	    if (asBounded.get())
	      {
		CurveBoundedDomain par_dom = asBounded->parameterDomain();
		int n_segs_u = surf_data->segs_u();
		int n_segs_v = surf_data->segs_v();
		int first_bb_idx = structure->n_boxes() - n_segs_u * n_segs_v;
		double eps = 1.0e-8;

		// Test boundary curve crossing for each fixed u-value
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

		    for (int j = 0; j <= inside_intervals.size(); ++j)
		      {
			int outside_from = 0;
			int outside_to = n_segs_v - 1;
			if (j > 0)
			  {
			    double par = inside_intervals[j-1].second;
			    for (outside_from = 0; outside_from < n_segs_v - 1 && segment_pars[1][outside_from+1] < par; ++outside_from);
			  }
			if (j < inside_intervals.size())
			  {
			    double par = inside_intervals[j].first;
			    for (outside_to = 0; outside_to < n_segs_v - 1 && segment_pars[1][outside_to+1] < par; ++outside_to);
			  }
			for (int k = outside_from; k <= outside_to; ++k)
			  {
			    if (i > 0)
			      structure->getBox(first_bb_idx + k*n_segs_u + (i - 1))->setInside(false);
			    if (i < n_segs_u)
			      structure->getBox(first_bb_idx + k*n_segs_u + i)->setInside(false);
			  }
		      }
		  }

		// Test boundary curve crossing for each fixed v-value
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
		    for (int j = 0; j <= inside_intervals.size(); ++j)
		      {
			int outside_from = 0;
			int outside_to = n_segs_u - 1;
			if (j > 0)
			  {
			    double par = inside_intervals[j-1].second;
			    for (outside_from = 0; outside_from < n_segs_u - 1 && segment_pars[0][outside_from+1] < par; ++outside_from);
			  }
			if (j < inside_intervals.size())
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

		// Add inside points
		double umin = segment_pars[0][0];
		double umax = segment_pars[0][n_segs_u];
		double vmin = segment_pars[1][0];
		double vmax = segment_pars[1][n_segs_v];
		int nmb_pts_each_dir = 20;
		double step_u = (umax - umin) / (double)nmb_pts_each_dir;
		double start_u = umin + 0.5 * step_u;
		double step_v = (vmax - vmin) / (double)nmb_pts_each_dir;
		double start_v = vmin + 0.5 * step_v;
		vector<vector<bool> > inside_mask(nmb_pts_each_dir);
		vector<vector<bool> > inside_mask2(nmb_pts_each_dir);
		Array<double, 2> pars;
		for (int i = 0; i < nmb_pts_each_dir; ++i)
		  {
		    pars[0] = start_u + step_u*(double)i;
		    for (int j = 0; j < nmb_pts_each_dir; ++j)
		      {
			pars[1] = start_v + step_v*(double)j;
			inside_mask[i].push_back(par_dom.isInDomain(pars, eps));
			inside_mask2[i].push_back(false);
		      }
		  }
		for (int i = 0; i < nmb_pts_each_dir; ++i)
		  for (int j = 0; j < nmb_pts_each_dir; ++j)
		    if (inside_mask[i][j])
		      {
			if (i == 0 || ! inside_mask[i-1][j] ||
			    j == 0 || ! inside_mask[i][j-1] ||
			    i == nmb_pts_each_dir-1 || ! inside_mask[i+1][j] ||
			    j == nmb_pts_each_dir-1 || ! inside_mask[i][j+1])
			  {
			    inside_mask2[i][j] = true;
			    surf_data->add_inside_point(paramSurf->point(start_u + step_u*(double)i, start_v + step_v*(double)j));
			  }
		      }
		/*
		cout << "Surface " << srf_idx << endl;
		for (int i = 0; i < nmb_pts_each_dir; ++i)
		  {
		    for (int j = 0; j < nmb_pts_each_dir; ++j)
		      cout << (inside_mask2[i][j] ? "X" : ".");
		    cout << endl;
		  }
		*/
	      }

	    ++srf_idx;
	  }
      }

    // Make voxel structure
    structure->BuildVoxelStructure(bigbox, 1000.0);

    cout << "Bounding boxes found = " << (structure->n_boxes()) << endl;
    int nx = structure->n_voxels_x();
    int ny = structure->n_voxels_y();
    int nz = structure->n_voxels_z();
    cout << "Number of voxels is " << nx << "*" << ny << "*" << nz << " = " << (nx*ny*nz) << endl;
    cout << "voxel_length = " << (structure->voxel_length()) << endl;

    clock_t t_after = clock();
    cout << endl << "Preprocessing timing = " << ((double)(t_after - t_before) / CLOCKS_PER_SEC) << " seconds" << endl;

    return structure;
  }

  namespace  // Anonymous
  {
    struct PossibleInside
    {
    public:
      Point cl_p_;
      int surf_idx_;
      double dist_;
      double par_u_;
      double par_v_;
      double up_lim_boundary_;

      PossibleInside()
      { }

      PossibleInside(Point cl_p, int surf_idx, double dist, double u, double v, double up_lim_bd)
	: cl_p_(cl_p), surf_idx_(surf_idx), dist_(dist), par_u_(u), par_v_(v), up_lim_boundary_(up_lim_bd)
      {
      }
    };
  }

  vector<float> closestVectors(const vector<float>& inPoints, const shared_ptr<BoundingBoxStructure>& boxStructure,
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
    for (int idx = 0, pt_idx = 0; idx < inPoints.size(); idx += 3, ++pt_idx)
      {
	if (pt_idx == start_idx)
	  skip_cnt = 0;
	if (pt_idx > max_idx)
	  skip_cnt = -1;
	if (skip_cnt == 0)
	  {
	    // cout << "Point index is " << idx << endl;
	    ++total_pts_tested;
	    Point pt(translation);
	    for (int i = 0; i < 3; ++i)
	      for (int j = 0; j < 3; ++j)
		pt[i] += rotationMatrix[i][j] * inPoints[idx + j];
	    double best_dist;
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

	    bool any_clp_found = false;

	    vector<PossibleInside> poss_in;

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

	    for (int vox_span = -1; vox_span <= nv_max; ++vox_span)
	      {
		double shortest_voxel_distance = (int)(vox_span-1) * voxel_length + shortest_face_distance;

		if (vox_span > 0)
		  {

		    // Test if any of the possible inside candidates are so close that they must be checked now
		    int poss_in_size = poss_in.size();
		    while(true)
		      {
			int best_poss_in = -1;
			for (int i = 0; i < poss_in_size; ++i)
			  if (poss_in[i].up_lim_boundary_ < shortest_voxel_distance || vox_span == nv_max)
			    {
			      if (best_poss_in == -1 || poss_in[i].dist_ < poss_in[best_poss_in].dist_)    // Not sure if this is best priority
				best_poss_in = i;
			    }
			if (best_poss_in == -1)
			  break;

			PossibleInside p_i = poss_in[best_poss_in];
			--poss_in_size;
			if (best_poss_in < poss_in_size)
			  poss_in[best_poss_in] = poss_in[poss_in_size];
			poss_in.resize(poss_in_size);

			double seed[2];
			seed[0] = p_i.par_u_;
			seed[1] = p_i.par_v_;
			double clo_u, clo_v;
			Point clo_pt;
			double clo_dist;
			shared_ptr<BoundedSurface> boundSurf = dynamic_pointer_cast<BoundedSurface>(boxStructure->getSurface(p_i.surf_idx_)->surface());
			boundSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, NULL, &seed[0]);
			/*
			cout << "New BS : surf = " << p_i.surf_idx_ << " seed = (" << seed[0] << ", " << seed[1]
			     <<") clo_par = (" << clo_u << ", " << clo_v << ")  clo_dist = " << clo_dist << endl;
			*/
			++local_boundaryCalls;

			if (!any_clp_found || clo_dist < best_dist)
			  {
			    best_dist = clo_dist;
			    best_idx = p_i.surf_idx_;
			    best_dist_pt = clo_pt - pt;
			    any_clp_found = true;
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
		  }

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

			    // Avoid 'internal' voxels
			    if (abs(vx - n_x) != vox_span_trunc && 
				abs(vy - n_y) != vox_span_trunc && 
				abs(vz - n_z) != vox_span_trunc)
			      continue;

			    // Break if this voxel is to far away
			    double d2_xyz = d2_xy + pt_dist_z[vz] * pt_dist_z[vz];
			    if (any_clp_found && d2_xyz > best_dist * best_dist)
			      continue;

			    // Get boxes to be tested
			    vector<int> possible_boxes;
			    if (vox_span >= 0)
			      possible_boxes = boxStructure->boxes_in_voxel(vx, vy, vz);
			    else
			      {
				// First iteration, only check with bounding boxes containing point
				vector<int> voxel_boxes = boxStructure->boxes_in_voxel(vx, vy, vz);
				for (int i = 0; i < voxel_boxes.size(); ++i)
				  {
				    int box_idx = voxel_boxes[i];
				    shared_ptr<SubSurfaceBoundingBox> surf_box = boxStructure->getBox(box_idx);
				    BoundingBox bb = boxStructure->getBox(box_idx)->box();
				    if (bb.containsPoint(pt))
				      possible_boxes.push_back(box_idx);
				  }
			      }

			    // Run through all boxes to be tested
			    for (int i = 0; i < possible_boxes.size(); ++i)
			      {
				int box_idx = possible_boxes[i];
				if (lastBoxCall[box_idx] == pt_idx)
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
				shared_ptr<ParamSurface> paramSurf = surf_data->surface();
				shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
				bool pt_might_be_outside = (bool)(boundedSurf.get());
				if (pt_might_be_outside)
				  paramSurf = boundedSurf->underlyingSurface();

				int segs_u = surf_data->segs_u();
				int segs_v = surf_data->segs_v();

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

				shared_ptr<RectDomain> rd = surf_box->par_domain();
				double seed[2];
				seed[0] = (rd->umin() + rd->umax()) * 0.5;
				seed[1] = (rd->vmin() + rd->vmax()) * 0.5;

				double clo_u, clo_v;
				Point clo_pt;
				double clo_dist;

				paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, search_domain.get(), &seed[0]);
				/*
				cout << "New PS : surf = " << (surf_data->index()) << " seed = (" << seed[0] << ", " << seed[1]
				     << ") clo_par = (" << clo_u << ", " << clo_v << ")  clo_dist = " << clo_dist
				     << "  domain = [" << (search_domain->umin()) << ", " << (search_domain->umax()) << "]x["
				     << (search_domain->vmin()) << ", " << (search_domain->vmax()) << "]  box_pos = ("
				     << (surf_box->pos_u()) << ", " << (surf_box->pos_v()) << ")  back = ("
				     << back_u << ", " << back_v << ")  len = ("
				     << len_u << ", " << len_v << ")" << endl;
				*/
				/*
				cout << "New PS(" << search_extend << ") : surf = " << (surf_data->index()) << "  box_idx = " << box_idx << " box_pos = ("
				     << (surf_box->pos_u()) << ", " << (surf_box->pos_v()) << ")  back = ("
				     << back_u << ", " << back_v << ")  len = ("
				     << len_u << ", " << len_v << ")" << endl;
				*/
				++local_underlyingCalls;
				for (int j = 0; j < len_u; ++j)
				  for (int k = 0; k < len_v; ++k)
				    lastBoxCall[ll_index + k * segs_u + j] = pt_idx;

				// If top surface is BoundedSurface, check if this point might be outside
				if (pt_might_be_outside)
				  {
				    int pos_u, pos_v;
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
				    pt_might_be_outside = !(boxStructure->getBox(cl_p_box)->inside());
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
					for (int j = 0; j < poss_in.size() && insert; ++j)
					  insert = surf_idx != poss_in[j].surf_idx_ ||
					    abs(clo_u - poss_in[j].par_u_) > tol ||
					    abs(clo_v - poss_in[j].par_v_) > tol;

					if (insert)
					  {
					    double up_lim_b2 = voxel_length * voxel_length * (double)(nv_x*nv_x + nv_z*nv_z + nv_z*nv_z);
					    vector<Point> surf_pts = surf_data->inside_points();
					    for (int j = 0; j < surf_pts.size(); ++j)
					      {
						double dist2 = pt.dist2(surf_pts[j]);
						if (dist2 < up_lim_b2)
						  up_lim_b2 = dist2;
					      }
					    /*
					    // Temporary calculation of upper limit to boundary as evaluation of center of (underlying) surface
					    RectDomain surf_dom = paramSurf->containingDomain();
					    double mid_u = 0.5 * (surf_dom.umin() + surf_dom.umax());
					    double mid_v = 0.5 * (surf_dom.vmin() + surf_dom.vmax());
					    Point center_p = paramSurf->point(mid_u, mid_v);
					    double up_lim_b = center_p.dist(pt);
					    */

					    poss_in.push_back(PossibleInside(clo_pt, surf_idx, clo_dist, clo_u, clo_v, sqrt(up_lim_b2)));
					  }
				      }
				    else
				      {
					// The point is inside the parameter domain and the closest point foiund so far
					// Update information about closest point, and remove possible inside candidates that are too far away
					best_dist = clo_dist;
					best_idx = surf_data->index();
					best_dist_pt = clo_pt - pt;
					any_clp_found = true;
					int poss_in_size = poss_in.size();
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

	    ++isBest[best_idx];
	    result.push_back(best_dist);

	    if (local_boundaryCalls >= boundaryCalls.size())
	      boundaryCalls.resize(local_boundaryCalls + 1);
	    ++boundaryCalls[local_boundaryCalls];
	    if (local_underlyingCalls >= underlyingCalls.size())
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
	shared_ptr<ParamSurface> paramSurf = boxStructure->getSurface(i)->surface();
	shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
	if (boundedSurf.get())
	  paramSurf = boundedSurf->underlyingSurface();
	if (dynamic_pointer_cast<SplineSurface>(paramSurf).get())
	  type_str  = "Spline";
	else if (dynamic_pointer_cast<Sphere>(paramSurf).get())
	    type_str = "Sphere";
	else if (dynamic_pointer_cast<Cylinder>(paramSurf).get())
	    type_str = "Cylinder";
	cout << i << "\t" << isBest[i] << "\t" << type_str << endl;
      }

    cout << endl << "N\tBound\tNot b" << endl;
    int tot_bs = 0;
    int tot_ul = 0;
    for (int i = 0; i < boundaryCalls.size() || i < underlyingCalls.size(); ++i)
      {
	int bs_calls = 0;
	if (i < boundaryCalls.size())
	  bs_calls = boundaryCalls[i];
	int ul_calls = 0;
	if (i < underlyingCalls.size())
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
    for (int idx = 0; idx < result.size(); ++idx)
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
    for (int idx = 0, pt_idx = 0; idx < inPoints.size(); idx += 3, ++pt_idx)
      {
	if (pt_idx == start_idx)
	  skip_cnt = 0;
	if (pt_idx > max_idx)
	  skip_cnt = -1;
	if (skip_cnt == 0)
	  {
	    // cout << "Point index is " << idx << endl;
	    ++total_pts_tested;
	    Point pt(translation);
	    for (int i = 0; i < 3; ++i)
	      for (int j = 0; j < 3; ++j)
		pt[i] += rotationMatrix[i][j] * inPoints[idx + j];
	    double best_dist;
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
		for (int i = 0; i < possible_boxes.size(); ++i)
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
			  shared_ptr<ParamSurface> paramSurf = surf_data->surface();
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
			  // New code
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
			  // End new code
			  // Old code
			  /*
			  bool shall_test = true;
			  if (isInside)
			    paramSurf = boundedSurf->underlyingSurface();
			  else if (any_tested)
			    {
			      boundedSurf->underlyingSurface()->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, search_domain.get(), &seed[0]);
			      shall_test = clo_dist < best_dist;
			    }

			  if (shall_test)
			    {
			      paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, search_domain.get(), &seed[0]);
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
			  */
			  // End old code
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
			    for (int i = 0; i < possible_boxes.size(); ++i)
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
				shared_ptr<ParamSurface> paramSurf = surf_data->surface();
				shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
				double clo_u, clo_v;
				Point clo_pt;
				double clo_dist;
				shared_ptr<RectDomain> rd = surf_box->par_domain();
				double seed[2];
				seed[0] = (rd->umin() + rd->umax()) * 0.5;
				seed[1] = (rd->vmin() + rd->vmax()) * 0.5;

				// New code
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
				// End new code
				// Old code
				/*
				bool shall_test = true;
				if (surf_box->inside())
				  paramSurf = boundedSurf->underlyingSurface();
				else if (any_tested)
				  {
				    boundedSurf->underlyingSurface()->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, rd.get(), &seed[0]);
				    shall_test = clo_dist < best_dist;
				  }

				if (shall_test)
				  {
				    paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, rd.get(), &seed[0]);

				    if (surf_box->inside())
				      ++local_underlyingCalls;
				    else
				      ++local_boundaryCalls;
				    if (!any_tested || clo_dist < best_dist)
				      {
					best_dist = clo_dist;
					best_idx = surf_data->index();
					best_dist_pt = clo_pt - pt;
					voxels_close = best_dist > shortest_voxel_distance;
					any_tested = true;
				      }
				  }
				*/
				// End old code
				lastBoxCall[box_idx] = pt_idx;
			      }
			  }
		      }
		  }
	      }
	    ++isBest[best_idx];
	    result.push_back(best_dist);

	    if (local_boundaryCalls >= boundaryCalls.size())
	      boundaryCalls.resize(local_boundaryCalls + 1);
	    ++boundaryCalls[local_boundaryCalls];
	    if (local_underlyingCalls >= underlyingCalls.size())
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
	shared_ptr<ParamSurface> paramSurf = boxStructure->getSurface(i)->surface();
	shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
	if (boundedSurf.get())
	  paramSurf = boundedSurf->underlyingSurface();
	if (dynamic_pointer_cast<SplineSurface>(paramSurf).get())
	  type_str  = "Spline";
	else if (dynamic_pointer_cast<Sphere>(paramSurf).get())
	    type_str = "Sphere";
	else if (dynamic_pointer_cast<Cylinder>(paramSurf).get())
	    type_str = "Cylinder";
	cout << i << "\t" << isBest[i] << "\t" << type_str << endl;
      }

    cout << endl << "N\tBound\tNot b" << endl;
    int tot_bs = 0;
    int tot_ul = 0;
    for (int i = 0; i < boundaryCalls.size() || i < underlyingCalls.size(); ++i)
      {
	int bs_calls = 0;
	if (i < boundaryCalls.size())
	  bs_calls = boundaryCalls[i];
	int ul_calls = 0;
	if (i < underlyingCalls.size())
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
    for (int idx = 0; idx < result.size(); ++idx)
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


}   // end namespace Go
