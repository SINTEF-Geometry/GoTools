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
#include "GoTools/compositemodel/ClosestPointUtils.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/CompositeModel.h"

using namespace std;
using namespace Go;

namespace Go
{


  vector<float> closestVectors(const vector<float>& inPoints, const vector<shared_ptr<GeomObject> >& surfaces,
			       const vector<vector<double> >& rotationMatrix, const Point& translation, int skip, int test_type)
  {
    // test_type values:
    // 0            for each point, for each surface, find closest, pick best
    // 1            Create surface model, use model->closestPoint()
    // 2            Create bounding boxes for all surfaces, splines split into Bezier patches. Test limited to inside parameter domain. Afterwards test
    //              to all that might be as close
    vector<bool> isSpline;
    for (int i = 0; i < surfaces.size(); ++i)
      {
	shared_ptr<BoundedSurface> boundSurf = dynamic_pointer_cast<BoundedSurface>(surfaces[i]);
	shared_ptr<SplineSurface> splineSurf = dynamic_pointer_cast<SplineSurface>(boundSurf->underlyingSurface());
	if (splineSurf.get())
	  isSpline.push_back(true);
	else
	  isSpline.push_back(false);
      }
    clock_t t_start = clock();

    // ********** Preprocessing code goes here ******************
    shared_ptr<CompositeModel> comp_model;
    if (test_type == 1)
      {
	double gap = 0.001;
	double neighbour = 0.01;
	double kink = 0.01;
	double approx = 0.001;
	CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);
	ifstream infile("/home/kfp/Build/simpleViewer/Kaplan_blade_raw.g2");
    	comp_model = shared_ptr<CompositeModel>(factory.createFromG2(infile));
      }

    vector<BoundingBox> boxes;
    vector<shared_ptr<RectDomain> > par_domains;
    vector<int> domain_surfaces;
    double voxel_length;
    int nv_x, nv_y, nv_z, nv_max;
    Point big_vox_low;
    vector<vector<vector<vector<int> > > > boxes_in_voxel;
    if (test_type == 2)
      {
	BoundingBox bigbox(3);
	// ofstream ons("bigBounding.g2");
	for (int srf_idx = 0; srf_idx < surfaces.size(); ++srf_idx)
	  {
	    shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(surfaces[srf_idx]);
	    if (paramSurf.get())
	      {
		BoundingBox surfbox(3);
		shared_ptr<BoundedSurface> asBounded = dynamic_pointer_cast<BoundedSurface>(paramSurf);
		if (asBounded.get())
		  paramSurf = asBounded->underlyingSurface();
		shared_ptr<SplineSurface> splineSurf = dynamic_pointer_cast<SplineSurface>(paramSurf);
		shared_ptr<ElementarySurface> elSurf = dynamic_pointer_cast<ElementarySurface>(paramSurf);
		if (splineSurf.get())
		  {
		    // ofstream ons("box_lines.g2");
		    SplineSurface* surf_copy = splineSurf->clone();

		    // Make C0 to get Bezier data
		    vector<vector<double> > unique_knots(2);
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

			    unique_knots[dir].push_back(knot);
			    if (seg == 0)
			      seg_ctrl_pos[dir].push_back(0);
			    else if (it != basis.end())
			      seg_ctrl_pos[dir].push_back(seg_ctrl_pos[dir][seg - 1] + cnt);
			  }

			if (dir == 0)
			  surf_copy->insertKnot_u(new_knots);
			else
			  surf_copy->insertKnot_v(new_knots);

			// Fetch segment bounding boxes
			int n_segs_u = seg_ctrl_pos[0].size();
			int n_segs_v = seg_ctrl_pos[1].size();
			int n_coefs_u = surf_copy->numCoefs_u();

			vector<double>::const_iterator ctrl_it_v = surf_copy->ctrl_begin();
			for (int j = 0; j < n_segs_v; ++j)
			  {
			    double start_v = unique_knots[1][j];
			    double end_v = unique_knots[1][j+1];
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
				boxes.push_back(bb);
				bigbox.addUnionWith(bb);
				surfbox.addUnionWith(bb);
				/*
				ons << "410 1 0 4 255 0 0 255" << endl << 12 << endl;
				Point l = bb.low();
				Point h = bb.high();
				Point ex = Point(h[0]-l[0], 0.0, 0.0);
				Point ey = Point(0.0, h[1]-l[1], 0.0);
				Point ez = Point(0.0, 0.0, h[2]-l[2]);
				ons << l << " " << l+ex << endl;
				ons << l+ey << " " << l+ex+ey << endl;
				ons << l+ez << " " << l+ex+ez << endl;
				ons << l+ey+ez << " " << l+ex+ey+ez << endl;
				ons << l << " " << l+ey << endl;
				ons << l+ex << " " << l+ey+ex << endl;
				ons << l+ez << " " << l+ey+ez << endl;
				ons << l+ex+ez << " " << l+ey+ex+ez << endl;
				ons << l << " " << l+ez << endl;
				ons << l+ex << " " << l+ez+ex << endl;
				ons << l+ey << " " << l+ez+ey << endl;
				ons << l+ex+ey << " " << l+ez+ex+ey << endl;
				*/

				Array<double, 2> ll, ur;
				ll[0] = unique_knots[0][i];
				ll[1] = start_v;
				ur[0] = unique_knots[0][i];
				ur[1] = end_v;
				par_domains.push_back(shared_ptr<RectDomain>(new RectDomain(ll, ur)));
				domain_surfaces.push_back(srf_idx);
			      }
			  }
		      }

		    // Collect control points

		    delete surf_copy;
		    /*
		    ons.close();
		    exit(1);
		    */
		  }
		if (elSurf.get())
		  {
		    BoundingBox bb = elSurf->boundingBox();
		    boxes.push_back(bb);
		    surfbox.addUnionWith(bb);
		    bigbox.addUnionWith(bb);
		    par_domains.push_back(shared_ptr<RectDomain>(new RectDomain(elSurf->containingDomain())));
		    domain_surfaces.push_back(srf_idx);
		  }
		/*
		ons << "410 1 0 4 255 0 0 255" << endl << 12 << endl;
		Point l = surfbox.low();
		Point h = surfbox.high();
		Point ex = Point(h[0]-l[0], 0.0, 0.0);
		Point ey = Point(0.0, h[1]-l[1], 0.0);
		Point ez = Point(0.0, 0.0, h[2]-l[2]);
		ons << l << " " << l+ex << endl;
		ons << l+ey << " " << l+ex+ey << endl;
		ons << l+ez << " " << l+ex+ez << endl;
		ons << l+ey+ez << " " << l+ex+ey+ez << endl;
		ons << l << " " << l+ey << endl;
		ons << l+ex << " " << l+ey+ex << endl;
		ons << l+ez << " " << l+ey+ez << endl;
		ons << l+ex+ez << " " << l+ey+ex+ez << endl;
		ons << l << " " << l+ez << endl;
		ons << l+ex << " " << l+ez+ex << endl;
		ons << l+ey << " " << l+ez+ey << endl;
		ons << l+ex+ey << " " << l+ez+ex+ey << endl;
		*/
	      }
	  }
	cout << "Surfaces found = " << domain_surfaces.size() << endl;
	/*
	ons << "410 1 0 4 0 255 0 255" << endl << 12 << endl;
	Point l = bigbox.low();
	Point h = bigbox.high();
	Point ex = Point(h[0]-l[0], 0.0, 0.0);
	Point ey = Point(0.0, h[1]-l[1], 0.0);
	Point ez = Point(0.0, 0.0, h[2]-l[2]);
	ons << l << " " << l+ex << endl;
	ons << l+ey << " " << l+ex+ey << endl;
	ons << l+ez << " " << l+ex+ez << endl;
	ons << l+ey+ez << " " << l+ex+ey+ez << endl;
	ons << l << " " << l+ey << endl;
	ons << l+ex << " " << l+ey+ex << endl;
	ons << l+ez << " " << l+ey+ez << endl;
	ons << l+ex+ez << " " << l+ey+ex+ez << endl;
	ons << l << " " << l+ez << endl;
	ons << l+ex << " " << l+ez+ex << endl;
	ons << l+ey << " " << l+ez+ey << endl;
	ons << l+ex+ey << " " << l+ez+ex+ey << endl;
	ons.close();
	*/

	// Make voxel structure
	Point diagonal = bigbox.high() - bigbox.low();
	double volume = diagonal[0] * diagonal[1] * diagonal[2];
	voxel_length = pow(volume/1000.0, 1.0/3.0);
	nv_x = (int)(1.0 + diagonal[0] / voxel_length);
	nv_y = (int)(1.0 + diagonal[1] / voxel_length);
	nv_z = (int)(1.0 + diagonal[2] / voxel_length);
	nv_max = nv_x;
	if (nv_y > nv_max)
	  nv_max = nv_y;
	if (nv_z > nv_max)
	  nv_max = nv_z;
	cout << "Number of voxels is " << nv_x << "*" << nv_y << "*" << nv_z << " = " << (nv_x*nv_y*nv_z) << endl;
	Point big_vox_center = bigbox.low() + diagonal * 0.5;
	big_vox_low = big_vox_center - Point((double)nv_x, (double)nv_y, (double)nv_z) * (0.5 * voxel_length);
	/*
	ofstream ons("bigVoxel.g2");
	for (int i = 0; i < 2; ++i)
	  {
	    if (i == 0)
	      ons << "410 1 0 4 255 255 0 255" << endl << 12 << endl;
	    else
	      ons << "410 1 0 4 255 0 255 255" << endl << 12 << endl;
	    Point l = big_vox_low;
	    double vx = (i == 0) ? (double)nv_x : 1.0;
	    double vy = (i == 0) ? (double)nv_y : 1.0;
	    double vz = (i == 0) ? (double)nv_z : 1.0;
	    Point ex = Point(vx * voxel_length, 0.0, 0.0);
	    Point ey = Point(0.0, vy * voxel_length, 0.0);
	    Point ez = Point(0.0, 0.0, vz * voxel_length);
	    ons << l << " " << l+ex << endl;
	    ons << l+ey << " " << l+ex+ey << endl;
	    ons << l+ez << " " << l+ex+ez << endl;
	    ons << l+ey+ez << " " << l+ex+ey+ez << endl;
	    ons << l << " " << l+ey << endl;
	    ons << l+ex << " " << l+ey+ex << endl;
	    ons << l+ez << " " << l+ey+ez << endl;
	    ons << l+ex+ez << " " << l+ey+ex+ez << endl;
	    ons << l << " " << l+ez << endl;
	    ons << l+ex << " " << l+ez+ex << endl;
	    ons << l+ey << " " << l+ez+ey << endl;
	    ons << l+ex+ey << " " << l+ez+ex+ey << endl;
	  }
	ons.close();
	*/

	// Add bounding boxes to voxel structure
	boxes_in_voxel.resize(nv_x);
	for (int i = 0; i < nv_x; ++i)
	  {
	    boxes_in_voxel[i].resize(nv_y);
	    for (int j = 0; j < nv_y; ++j)
	      boxes_in_voxel[i][j].resize(nv_z);
	  }
	for (int i = 0; i < boxes.size(); ++i)
	  {
	    BoundingBox bb = boxes[i];
	    Point l_rel = (bb.low() - big_vox_low) / voxel_length;
	    Point h_rel = (bb.high() - big_vox_low) / voxel_length;
	    int l_x = (int)(l_rel[0]);
	    int l_y = (int)(l_rel[1]);
	    int l_z = (int)(l_rel[2]);
	    int h_x = (int)(h_rel[0]);
	    int h_y = (int)(h_rel[1]);
	    int h_z = (int)(h_rel[2]);
	    for (int jx = l_x; jx <= h_x; ++jx)
	      for (int jy = l_y; jy <= h_y; ++jy)
		for (int jz = l_z; jz <= h_z; ++jz)
		  boxes_in_voxel[jx][jy][jz].push_back(i);
	  }
      }
    // ********** End preprocessing *********************

    clock_t t_after_preproc = clock();

    // ********** Closest point point code goes here ***************
    vector<float> result;
    double skip_cnt = 0;
    vector<shared_ptr<GeomObject> >::const_iterator surf_begin = surfaces.begin();
    vector<shared_ptr<GeomObject> >::const_iterator surf_end = surfaces.end();
    int total_pts_tested = 0;
    vector<int> isBest(surfaces.size(), 0);
    vector<int> splineCalls;
    vector<int> elementaryCalls;
    vector<vector<int> > combinedCalls;
    for (int idx = 0; idx < inPoints.size(); idx += 3)
      {
	if (skip_cnt == 0)
	  {
	    ++total_pts_tested;
	    Point pt(translation);
	    for (int i = 0; i < 3; ++i)
	      for (int j = 0; j < 3; ++j)
		pt[i] += rotationMatrix[i][j] * inPoints[idx + j];
	    double best_dist;
	    int best_idx = -1;
	    int local_splineCalls = 0;
	    int local_elementaryCalls = 0;
	    if (test_type == 0)
	      {
		bool first = true;
		int srf_cnt = 0;
		for (vector<shared_ptr<GeomObject> >::const_iterator it = surf_begin; it != surf_end; ++it, ++srf_cnt)
		  {
		    shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(*it);
		    // cout << "S" << srf_cnt << ": ";
		    if (paramSurf)
		      {
			/*
			  shared_ptr<SplineSurface> splineSurf = dynamic_pointer_cast<SplineSurface>(*it);
			  shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(*it);
			  if (splineSurf)
			  cout << "SplineSurface" << endl;
			  else if (boundedSurf)
			  cout << "BoundedSurface" << endl;
			  else
			  cout << "Unknown type ParamSurface" << endl;
			*/
			double clo_u, clo_v;
			Point clo_pt;
			double clo_dist;
			// cout << "Starting closestPoint" << endl;
			paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8);
			if (isSpline[srf_cnt])
			  ++local_splineCalls;
			else
			  ++local_elementaryCalls;
			// cout << "Done closestPoint" << endl;
			Point dist_vec = clo_pt - pt;
			RectDomain rd = paramSurf->asSplineSurface()->parameterDomain();
			/*
			  cout << "parameters = (" << clo_u << ", " << clo_v << " in [" << rd.umin() << ", " << rd.umax() << "]x[" << rd.vmin() << ", " << rd.vmax()
			  << "]  Inside? " << (inside_exception ? "Exception" : (inside ? "Yes" : "No")) << endl;
			  cout << "Distance Length2 = " << dist_vec.length2()
			  << "  New closest ? " << ((first || best_vec.length2() > dist_vec.length2()) ? "Yes" : "No") << endl;
			*/
			if (first || best_dist > clo_dist)
			  {
			    best_dist = clo_dist;
			    best_idx = srf_cnt;
			  }
			first = false;
		      }
		    /*
		      else
		      cout << "Unknown geometric object" << endl;
		    */
		  }
	      }
	    else if (test_type == 1)
	      {
		Point clp;
		double clo_p[2];
		comp_model->closestPoint(pt, clp, best_idx, clo_p, best_dist);
	      }
	    else if (test_type == 2)
	      {
		Point pt_rel = (pt - big_vox_low) / voxel_length;
		int n_x = (int)(pt_rel[0]);
		int n_y = (int)(pt_rel[1]);
		int n_z = (int)(pt_rel[2]);

		double shortest_face_distance = pt_rel[0] - (int)n_x * voxel_length;
		if (shortest_face_distance > 0.5 * voxel_length)
		  shortest_face_distance = voxel_length - shortest_face_distance;
		double next_face_distance = pt_rel[1] - (int)n_y * voxel_length;
		if (next_face_distance > 0.5 * voxel_length)
		  next_face_distance = voxel_length - shortest_face_distance;
		if (shortest_face_distance > next_face_distance)
		  shortest_face_distance = next_face_distance;
		next_face_distance = pt_rel[2] - (int)n_z * voxel_length;
		if (next_face_distance > 0.5 * voxel_length)
		  next_face_distance = voxel_length - shortest_face_distance;
		if (shortest_face_distance > next_face_distance)
		  shortest_face_distance = next_face_distance;

		vector<double> tested;

		// First test all bounding boxes containing the point
		vector<int> possible_boxes = boxes_in_voxel[n_x][n_y][n_z];
		for (int i = 0; i < possible_boxes.size(); ++i)
		  {
		    int box_idx = possible_boxes[i];
		    if (boxes[box_idx].containsPoint(pt))
		      {
			shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(surfaces[domain_surfaces[box_idx]]);
			double clo_u, clo_v;
			Point clo_pt;
			double clo_dist;
			paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, par_domains[box_idx].get());
			if (isSpline[domain_surfaces[box_idx]])
			  ++local_splineCalls;
			else
			  ++local_elementaryCalls;
			if (tested.size() == 0 || clo_dist < best_dist)
			  {
			    best_dist = clo_dist;
			    best_idx = domain_surfaces[box_idx];
			  }
			tested.push_back(box_idx);
		      }
		  }

		// Then test rest in this and other voxels
		for (int vox_span = 0; vox_span < nv_max; ++vox_span)
		  {
		    double shortest_voxel_distance = (int)(vox_span-1) * voxel_length + shortest_face_distance;
		    if (shortest_voxel_distance < 0.0)
		      shortest_voxel_distance = 0.0;

		    if (tested.size() > 0 && best_dist < shortest_voxel_distance)
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
		      for (int vy = beg_y; vy <= end_y && voxels_close; ++vy)
			for (int vz = beg_z; vz <= end_z && voxels_close; ++vz)
			  {
			    if (abs(vx - n_x) != vox_span && 
				abs(vx - n_x) != vox_span && 
				abs(vx - n_x) != vox_span)
			      continue;
			    possible_boxes = boxes_in_voxel[vx][vy][vz];
			    for (int i = 0; i < possible_boxes.size(); ++i)
			      {
				int box_idx = possible_boxes[i];
				if (tested.size() > 0)
				  {
				    if (find(tested.begin(), tested.end(), box_idx) != tested.end())
				      continue;
				    BoundingBox bb = boxes[box_idx];
				    Point low = bb.low();
				    Point high = bb.high();
				    double d2_pt_box = 0.0;
				    for (int j = 0; j < 3; ++j)
				      {
					double dist = (pt[j] < low[j]) ? (low[j] - pt[j]) : ((pt[j] > high[j]) ? (high[j] - pt[j]) : 0.0);
					d2_pt_box += dist * dist;
				      }
				    if (d2_pt_box > best_dist * best_dist)
				      {
					tested.push_back(box_idx);
					continue;
				      }
				  }

				// New box not tested before
				shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(surfaces[domain_surfaces[box_idx]]);
				double clo_u, clo_v;
				Point clo_pt;
				double clo_dist;
				paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, par_domains[box_idx].get());
				if (isSpline[domain_surfaces[box_idx]])
				  ++local_splineCalls;
				else
				  ++local_elementaryCalls;
				if (tested.size() == 0 || clo_dist < best_dist)
				  {
				    best_dist = clo_dist;
				    best_idx = domain_surfaces[box_idx];
				    voxels_close = best_dist > shortest_voxel_distance;
				  }
				tested.push_back(box_idx);
			      }
			  }
		  }
	      }
	    ++isBest[best_idx];
	    result.push_back(best_dist);
	    if (local_splineCalls >= splineCalls.size())
	      splineCalls.resize(local_splineCalls + 1);
	    ++splineCalls[local_splineCalls];
	    if (local_elementaryCalls >= elementaryCalls.size())
	      elementaryCalls.resize(local_elementaryCalls + 1);
	    ++elementaryCalls[local_elementaryCalls];
	    if (local_splineCalls >= combinedCalls.size())
	      combinedCalls.resize(local_splineCalls + 1);
	    if (local_elementaryCalls >= combinedCalls[local_splineCalls].size())
	      combinedCalls[local_splineCalls].resize(local_elementaryCalls + 1);
	    ++combinedCalls[local_splineCalls][local_elementaryCalls];
	  }

	++skip_cnt;
	if (skip_cnt == skip)
	  skip_cnt = 0;
      }
    // *********** End closest point code **********************
    clock_t t_after_closest = clock();
    cout << endl << "Preprocessing timing = " << ((double)(t_after_preproc - t_start) / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << "Closest point timing = " << ((double)(t_after_closest - t_after_preproc) / CLOCKS_PER_SEC) << " seconds" << endl;

    cout << "Test type = " << test_type << endl;
    cout << "Number of points tested for closestPoint = " << total_pts_tested << endl;

    cout << "Surf\tBest" << endl;
    for (int i = 0; i < surfaces.size(); ++i)
      cout << i << "\t" << isBest[i] << endl;

    cout << endl << "N\tS_call\tE-call" << endl;
    int tot_spl = 0;
    int tot_el = 0;
    for (int i = 0; i < elementaryCalls.size() || i < splineCalls.size(); ++i)
      {
	int scalls = 0;
	if (i < splineCalls.size())
	  scalls = splineCalls[i];
	int ecalls = 0;
	if (i < elementaryCalls.size())
	  ecalls = elementaryCalls[i];
	int tot_calls = scalls + ecalls;
	if (tot_calls > 0)
	  {
	    cout << i << "\t";
	    if (scalls > 0)
	      cout << scalls;
	    cout << "\t";
	    if (ecalls > 0)
	      cout << ecalls;
	    cout << endl;
	  }
	tot_spl += i * scalls;
	tot_el += i * ecalls;
      }
    cout << "Total\t" << tot_spl << "\t" << tot_el << endl;

    cout << endl << "S-N\tE-N\tCalls" << endl;
    for (int i = 0; i < combinedCalls.size(); ++i)
      for (int j = 0; j < combinedCalls[i].size(); ++j)
	if (combinedCalls[i][j] > 0)
	  cout << i << "\t" << j << "\t" << combinedCalls[i][j] << endl;

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
