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
#include "GoTools/compositemodel/ClosestPointUtils.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/CompositeModel.h"

using namespace std;
using namespace Go;

namespace Go
{


  int box_list_compare(const void* el1, const void* el2)
  {
    double diff = (*(pair<int, double>*)el2).second - (*(pair<int, double>*)el1).second;
    if (diff < 0.0)
      return -1;
    if (diff > 0.0)
      return 1;
    return 0;
  }

  vector<float> closestVectors(const vector<float>& inPoints, const vector<shared_ptr<GeomObject> >& surfaces,
			       const vector<vector<double> >& rotationMatrix, const Point& translation,
			       int test_type, int start_idx, int skip, int max_idx, double par_len_el)
  {
    // test_type values:
    // 0            for each point, for each surface, find closest, pick best
    // 1            Create surface model, use model->closestPoint()
    // 2            Create bounding boxes for all surfaces, splines split into Bezier patches. Test limited to inside parameter domain. Afterwards test
    //              to all that might be as close
    vector<bool> surfType;
    for (int i = 0; i < surfaces.size(); ++i)
      {
	shared_ptr<BoundedSurface> boundSurf = dynamic_pointer_cast<BoundedSurface>(surfaces[i]);
	shared_ptr<SplineSurface> splineSurf = dynamic_pointer_cast<SplineSurface>(boundSurf->underlyingSurface());
	shared_ptr<Sphere> sphere = dynamic_pointer_cast<Sphere>(boundSurf->underlyingSurface());
	shared_ptr<Cylinder> cyl = dynamic_pointer_cast<Cylinder>(boundSurf->underlyingSurface());
	if (splineSurf.get())
	  surfType.push_back(0);
	else if (sphere.get())
	  surfType.push_back(1);
	else if (cyl.get())
	  surfType.push_back(2);
	else
	  surfType.push_back(3);
      }
    clock_t t_start = clock();

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

    vector<int> surf_segs_u;
    vector<int> surf_segs_v;

    vector<BoundingBox> boxes;
    vector<shared_ptr<RectDomain> > par_domains;
    vector<int> domain_surfaces;
    vector<int> domain_pos_u;
    vector<int> domain_pos_v;
    vector<bool> domain_inside_boundary;
    double voxel_length;
    int nv_x, nv_y, nv_z, nv_max;
    vector<double> pt_dist_x;
    vector<double> pt_dist_y;
    vector<double> pt_dist_z;
    Point big_vox_low;
    vector<vector<vector<vector<int> > > > boxes_in_voxel;
    int extend_search_field = 1;
    if (test_type == 2 || test_type == 3 || test_type == 4)
      {
	BoundingBox bigbox(3);
	// ofstream ons("bigBounding.g2");
	for (int srf_idx = 0; srf_idx < surfaces.size(); ++srf_idx)
	  {
	    // cout << endl << "Start of surface " << srf_idx << endl;
	    shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(surfaces[srf_idx]);
	    if (paramSurf.get())
	      {
		BoundingBox surfbox(3);
		shared_ptr<BoundedSurface> asBounded = dynamic_pointer_cast<BoundedSurface>(paramSurf);
		if (asBounded.get())
		  paramSurf = asBounded->underlyingSurface();
		shared_ptr<SplineSurface> splineSurf = dynamic_pointer_cast<SplineSurface>(paramSurf);
		shared_ptr<ElementarySurface> elSurf = dynamic_pointer_cast<ElementarySurface>(paramSurf);
		vector<vector<double> > segment_pars(2);
		if (splineSurf.get())
		  {
		    // ofstream ons("box_lines.g2");
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
		    surf_segs_u.push_back(n_segs_u);
		    surf_segs_v.push_back(n_segs_v);
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
			    boxes.push_back(bb);
			    bigbox.addUnionWith(bb);
			    surfbox.addUnionWith(bb);
			    domain_inside_boundary.push_back(true);
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
			    ll[0] = segment_pars[0][i];
			    ll[1] = start_v;
			    ur[0] = segment_pars[0][i];
			    ur[1] = end_v;
			    par_domains.push_back(shared_ptr<RectDomain>(new RectDomain(ll, ur)));
			    domain_surfaces.push_back(srf_idx);

			    domain_pos_u.push_back(i);
			    domain_pos_v.push_back(j);
			  }
		      }

		    // Collect control points

		    // delete surf_copy;
		    /*
		    ons.close();
		    exit(1);
		    */
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
		    surf_segs_u.push_back(n_segs_u);
		    surf_segs_v.push_back(n_segs_v);

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
			    boxes.push_back(bb);
			    surfbox.addUnionWith(bb);
			    bigbox.addUnionWith(bb);

			    Array<double, 2> ll, ur;
			    ll[0] = seg_umin;
			    ll[1] = seg_vmin;
			    ur[0] = seg_umax;
			    ur[1] = seg_vmax;
			    par_domains.push_back(shared_ptr<RectDomain>(new RectDomain(ll, ur)));
			    domain_surfaces.push_back(srf_idx);
			    domain_pos_u.push_back(i);
			    domain_pos_v.push_back(j);
			    domain_inside_boundary.push_back(true);

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
			  }
		      }

		    /*
		    cout << "Boundary for elementary surface: [" << big_rd->umin() << ", " << big_rd->umax() << "]x["
			 << big_rd->vmin() << ", " << big_rd->vmax() << "]  swapped : "
			 << (elSurf->isSwapped() ? "Yes" : "No") << endl;
		    */
		  }

		// Determine if the parameter domains of the bounding boxes are entirely inside the parameter domain limited by the boundary curves
		if (asBounded.get())
		  {
		    CurveBoundedDomain par_dom = asBounded->parameterDomain();
		    int n_segs_u = surf_segs_u[(int)surf_segs_u.size() - 1];
		    int n_segs_v = surf_segs_v[(int)surf_segs_v.size() - 1];
		    int first_bb_idx = (int)boxes.size() - n_segs_u * n_segs_v;
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
			    // cout << "Exception u-dir for surface " << srf_idx << " i = " << i << "/" << n_segs_u << endl;
			  }
			/*
			cout << "Segments u-dir, i=" << i << endl;
			for (int j = 0; j < inside_intervals.size(); ++j)
			  cout << "   " << inside_intervals[j].first << " " << inside_intervals[j].second << endl;
			*/
			for (int j = 0; j <= inside_intervals.size(); ++j)
			  {
			    int outside_from ;
			    int outside_to;
			    outside_from = 0;
			    outside_to = n_segs_v - 1;
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
			    // cout << "Remove: " << i << "." << j << " [" << outside_from << " " << outside_to << "]" << endl;
			    for (int k = outside_from; k <= outside_to; ++k)
			      {
				if (i > 0)
				  domain_inside_boundary[first_bb_idx + k*n_segs_u + (i - 1)] = false;
				if (i < n_segs_u)
				  domain_inside_boundary[first_bb_idx + k*n_segs_u + i] = false;
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
			    // cout << "Exception v-dir for surface " << srf_idx << " i = " << i << "/" << n_segs_v << endl;
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
				  domain_inside_boundary[first_bb_idx + k + (i - 1)*n_segs_u] = false;
				if (i < n_segs_v)
				  domain_inside_boundary[first_bb_idx + k + i*n_segs_u] = false;
			      }
			  }
		      }

		    /*
		    if (srf_idx != 1)
		      {
			cout << endl << "Surface " << srf_idx << endl;
			for (int i = 0, pos = first_bb_idx; i < n_segs_v; ++i)
			  {
			    for (int j = 0; j < n_segs_u; ++j, ++pos)
			      cout << (domain_inside_boundary[pos] ? "T" : "F");
			    cout << endl;
			  }
		      }
		    */
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
	cout << "Bounding boxes found = " << domain_surfaces.size() << endl;
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
	if (test_type == 3)
	  voxel_length = 1.1 * max(diagonal[0], max(diagonal[1], diagonal[2]));
	else
	  voxel_length = pow(volume/1000.0, 1.0/3.0);
	nv_x = (int)(1.0 + diagonal[0] / voxel_length);
	nv_y = (int)(1.0 + diagonal[1] / voxel_length);
	nv_z = (int)(1.0 + diagonal[2] / voxel_length);
	nv_max = nv_x;
	if (nv_y > nv_max)
	  nv_max = nv_y;
	if (nv_z > nv_max)
	  nv_max = nv_z;
	pt_dist_x.resize(nv_x);
	pt_dist_y.resize(nv_y);
	pt_dist_z.resize(nv_z);
	cout << "Number of voxels is " << nv_x << "*" << nv_y << "*" << nv_z << " = " << (nv_x*nv_y*nv_z) << endl;
	cout << "voxel_length = " << voxel_length << endl;
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
	/*
	cout << "surf_segs_*.size() = " << surf_segs_u.size() << " and " << surf_segs_v.size() << endl;
	for (int i = 0; i < surf_segs_u.size(); ++i)
	  cout << "surf_segs_*[" << i << "] = " << surf_segs_u[i] << " " << surf_segs_v[i] << endl;
	*/
      }
    // ********** End preprocessing *********************

    clock_t t_after_preproc = clock();

    // ********** Closest point point code goes here ***************
    vector<float> result;
    double skip_cnt = -1;
    vector<shared_ptr<GeomObject> >::const_iterator surf_begin = surfaces.begin();
    vector<shared_ptr<GeomObject> >::const_iterator surf_end = surfaces.end();
    int total_pts_tested = 0;
    vector<int> isBest(surfaces.size(), 0);
    vector<int> boundaryCalls;
    vector<int> underlyingCalls;
    // vector<vector<int> > combinedCalls;
    int find_calls = 0;
    int find_len = 0;
    int bb_calls = 0;
    int cp_calls = 0;
    int best_inside = 0;
    for (int idx = 0, pt_idx = 0; idx < inPoints.size(); idx += 3, ++pt_idx)
      {
	if (pt_idx == start_idx)
	  skip_cnt = 0;
	if (pt_idx > max_idx)
	  skip_cnt = -1;
	if (skip_cnt == 0)
	  {
	    // cout << "pt_idx = " << pt_idx << endl;
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
	    bool local_best_inside = false;
	    if (test_type == 0)
	      {
		bool first = true;
		int srf_cnt = 0;
		for (vector<shared_ptr<GeomObject> >::const_iterator it = surf_begin; it != surf_end; ++it, ++srf_cnt)
		  {
		    shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(*it);
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
			++local_boundaryCalls;
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
			    best_dist_pt = dist_vec;
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
		best_dist_pt = clp - pt;
		/*
		cout << endl << "Input: " << inPoints[idx] << " " << inPoints[idx+1] << " " << inPoints[idx+2] << endl;
		cout << "Transformed: " << pt << endl << "Closest: " << clp << endl << "Diff: " << clp - pt << endl;
		*/
	      }
	    else if (test_type == 100)
	      {
		best_idx = 0;
		best_dist = 1.0;
		best_dist_pt = Point(1.0, 0.0, 0.0);
	      }
	    else if (test_type == 2 || test_type == 3 || test_type == 4)
	      {
		Point pt_vox_low = pt - big_vox_low;
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

		vector<double> tested;

		// First test all bounding boxes containing the point
		if (n_x >= 0 && n_x < nv_x &&
		    n_y >= 0 && n_y < nv_y &&
		    n_z >= 0 && n_z < nv_z)
		  {
		    vector<int> possible_boxes = boxes_in_voxel[n_x][n_y][n_z];

		    vector<pair<int, double> > inside_boxes;   // All boxes (first) sorted by distance from point to most distant corner (second)
		    for (int i = 0; i < possible_boxes.size(); ++i)
		      {
			int box_idx = possible_boxes[i];
			BoundingBox bb = boxes[box_idx];
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
		    /*
		    cout << "Inside boxes before sort (size = " << inside_boxes.size() << ")" << endl;
		    for (int i = 0; i < inside_boxes.size(); ++i)
		      cout << "  " << inside_boxes[i].first << " " << inside_boxes[i].second << endl;
		    cout << endl;
		    */
		    for (int i = 0; i < nmb_inside_boxes - 1; ++i)
		      for (int j = nmb_inside_boxes - 1; j > i; --j)
			if (inside_boxes[j].second < inside_boxes[j-1].second)
			  swap(inside_boxes[j], inside_boxes[j-1]);
		    /*
		    cout << "Inside boxes after sort (size = " << inside_boxes.size() << ")" << endl;
		    for (int i = 0; i < inside_boxes.size(); ++i)
		      cout << "  " << inside_boxes[i].first << " " << inside_boxes[i].second << endl;
		    cout << endl;
		    */

		    int total_turns = (test_type == 4) ? 2 : 1;
		    for (int turn = 0; turn < total_turns; ++turn)
		      for (int i = 0; i < nmb_inside_boxes; ++i)
			{
			  int box_idx = inside_boxes[i].first;
			  if (find(tested.begin(), tested.end(), box_idx) == tested.end())
			    {
			      // cout << "Testing i = " << i << " of " << inside_boxes.size() << endl;
			      int surf_idx = domain_surfaces[box_idx];
			      int segs_u = surf_segs_u[surf_idx];
			      int segs_v = surf_segs_v[surf_idx];

			      int back_u = min(domain_pos_u[box_idx], extend_search_field);
			      int back_v = min(domain_pos_v[box_idx], extend_search_field);
			      int len_u = back_u + 1 + min(segs_u - (domain_pos_u[box_idx] + 1), extend_search_field);
			      int len_v = back_u + 1 + min(segs_v - (domain_pos_v[box_idx] + 1), extend_search_field);
			      int ll_index = box_idx - (back_v * segs_u + back_u);

			      /*
			      cout << box_idx << " " << surf_idx << " " << segs_u << " " << segs_v << " "
				   << domain_pos_u[box_idx] << " " << domain_pos_v[box_idx] << " "
				   << back_u << " " << back_v << " " << len_u << " " << len_v << " " << ll_index << endl;

			      for (int i = 0; i < 40; ++i)
				cout << "Box idx " << i << " : Domain = (" << domain_pos_u[i] << "," << domain_pos_v[i] << "]" << endl;
			      */

			      Array<double, 2> search_domain_ll, search_domain_ur;
			      search_domain_ll[0] = par_domains[ll_index]->umin();
			      search_domain_ll[1] = par_domains[ll_index]->vmin();
			      search_domain_ur[0] = par_domains[ll_index + len_u - 1]->umax();
			      search_domain_ur[1] = par_domains[ll_index + (len_v - 1)*segs_u]->vmax();
			      shared_ptr<RectDomain> search_domain(new RectDomain(search_domain_ll, search_domain_ur));
			      bool isInside = true;
			      for (int j = 0; j < len_v && isInside; ++j)
				for (int k = 0; k < len_u && isInside; ++k)
				  isInside &= domain_inside_boundary[ll_index + j * segs_u + k];
			      if (test_type == 4 && (isInside == (turn == 1)))
				continue;

			      RectDomain* rd = par_domains[box_idx].get();
			      double seed[2];
			      seed[0] = (rd->umin() + rd->umax()) * 0.5;
			      seed[1] = (rd->vmin() + rd->vmax()) * 0.5;

			      shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(surfaces[surf_idx]);
			      shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
			      bool shall_test = true;
			      double clo_u, clo_v;
			      Point clo_pt;
			      double clo_dist;
			      if (isInside)
				  paramSurf = boundedSurf->underlyingSurface();
			      else if (test_type == 4 && tested.size() > 0)
				{
				  boundedSurf->underlyingSurface()->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, search_domain.get(), &seed[0]);
				  shall_test = clo_dist < best_dist;
				}
			      /*
				if (surf_idx != 12)
				{
			      */
			      if (shall_test)
				{
				  paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, search_domain.get(), &seed[0]);
				  if (isInside)
				    ++local_underlyingCalls;
				  else
				    ++local_boundaryCalls;
				  if (tested.size() == 0 || clo_dist < best_dist)
				    {
				      best_dist = clo_dist;
				      best_idx = domain_surfaces[box_idx];
				      best_dist_pt = clo_pt - pt;
				      local_best_inside = isInside;
				    }
				}
			      for (int j = 0; j < len_u; ++j)
				for (int k = 0; k < len_v; ++k)
				  tested.push_back(ll_index + j * segs_v + k);
			      /*
				}
			      */
			    }
			}
		  }

		// Then test rest in this and other voxels
		if (test_type == 2 || test_type == 4)
		  {
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
				    vector<int> possible_boxes = boxes_in_voxel[vx][vy][vz];
				    for (int i = 0; i < possible_boxes.size(); ++i)
				      {
					int box_idx = possible_boxes[i];
					if (tested.size() > 0)
					  {
					    ++find_calls;
					    find_len += tested.size();
					    if (find(tested.begin(), tested.end(), box_idx) != tested.end())
					      continue;
					    ++bb_calls;
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
						// tested.push_back(box_idx);
						continue;
					      }
					    ++cp_calls;
					  }

					// New box not tested before
					shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(surfaces[domain_surfaces[box_idx]]);
					shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
					bool shall_test = true;
					double clo_u, clo_v;
					Point clo_pt;
					double clo_dist;
					RectDomain* rd = par_domains[box_idx].get();
					double seed[2];
					seed[0] = (rd->umin() + rd->umax()) * 0.5;
					seed[1] = (rd->vmin() + rd->vmax()) * 0.5;

					if (domain_inside_boundary[box_idx])
					  paramSurf = boundedSurf->underlyingSurface();
					else if (test_type == 4 && tested.size() > 0)
					  {
					    boundedSurf->underlyingSurface()->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, rd, &seed[0]);
					    shall_test = clo_dist < best_dist;
					  }

					/*
					  if (domain_surfaces[box_idx] != 12)
					  {
					*/
					if (shall_test)
					  {
					    paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, rd, &seed[0]);
					    if (domain_inside_boundary[box_idx])
					      ++local_underlyingCalls;
					    else
					      ++local_boundaryCalls;
					    if (tested.size() == 0 || clo_dist < best_dist)
					      {
						best_dist = clo_dist;
						best_idx = domain_surfaces[box_idx];
						best_dist_pt = clo_pt - pt;
						local_best_inside = domain_inside_boundary[box_idx];
						voxels_close = best_dist > shortest_voxel_distance;
					      }
					  }
					tested.push_back(box_idx);
					/*
					  }
					*/
				      }
				  }
			      }
			  }
		      }
		  }
		else if (test_type == 3)
		  {
		    // Only one voxel. We sort all bounding boxes not hitting the point by quick sort. Sort key is distance to point,
		    vector<pair<int, double> > box_list_sorted;
		    for (int i = 0; i < boxes.size(); ++i)
		      {
			BoundingBox bb = boxes[i];
			if (!bb.containsPoint(pt))
			  {
			    double dist2 = 0.0;
			    Point ll_p = bb.low() - pt;
			    Point p_hi = pt - bb.high();
			    for (int j = 0; j < 3; ++j)
			      {
				double d_add = max(0.0, max(p_hi[j], ll_p[j]));
				dist2 += d_add * d_add;
			      }
			    box_list_sorted.push_back(pair<int, double>(i, dist2));
			  }
		      }
		    qsort(&box_list_sorted[0], (int)(box_list_sorted.size()), sizeof(pair<int, double>), box_list_compare);
		    bool any_tested = tested.size() > 0;
		    double best_dist2 = best_dist * best_dist;
		    for (int i = 0; i < box_list_sorted.size() && (!any_tested || best_dist2 > box_list_sorted[i].second); ++i)
		      {
			int box_idx = box_list_sorted[i].first;
			shared_ptr<ParamSurface> paramSurf = dynamic_pointer_cast<ParamSurface>(surfaces[domain_surfaces[box_idx]]);
			double clo_u, clo_v;
			Point clo_pt;
			double clo_dist;
			RectDomain* rd = par_domains[box_idx].get();
			double seed[2];
			seed[0] = (rd->umin() + rd->umax()) * 0.5;
			seed[1] = (rd->vmin() + rd->vmax()) * 0.5;

			if (domain_inside_boundary[box_idx])
			  {
			    shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
			    if (boundedSurf.get())
			      {
				paramSurf = boundedSurf->underlyingSurface();
				++local_underlyingCalls;
			      }
			    else
			      ++local_boundaryCalls;
			  }

			paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, 1.0e-8, rd, &seed[0]);
			if (!any_tested || clo_dist < best_dist)
			  {
			    best_dist = clo_dist;
			    best_idx = domain_surfaces[box_idx];
			    best_dist_pt = clo_pt - pt;
			    local_best_inside = domain_inside_boundary[box_idx];
			    any_tested = true;
			  }
		      }
		  }
	      }
	    ++isBest[best_idx];
	    if (local_best_inside)
	      ++best_inside;
	    result.push_back(best_dist);
	    // cout << best_dist_pt << endl;
	    if (local_boundaryCalls >= boundaryCalls.size())
	      boundaryCalls.resize(local_boundaryCalls + 1);
	    ++boundaryCalls[local_boundaryCalls];
	    if (local_underlyingCalls >= underlyingCalls.size())
	      underlyingCalls.resize(local_underlyingCalls + 1);
	    ++underlyingCalls[local_underlyingCalls];
	    /*
	    if (local_boundaryCalls >= combinedCalls.size())
	      combinedCalls.resize(local_boundaryCalls + 1);
	    if (local_underlyingCalls >= combinedCalls[local_boundaryCalls].size())
	      combinedCalls[local_boundaryCalls].resize(local_underlyingCalls + 1);
	    ++combinedCalls[local_boundaryCalls][local_underlyingCalls];
	    */
	  }

	if (skip_cnt != -1)
	  {
	    ++skip_cnt;
	    if (skip_cnt == skip)
	      skip_cnt = 0;
	  }
      }
    // *********** End closest point code **********************
    clock_t t_after_closest = clock();
    cout << endl << "Preprocessing timing = " << ((double)(t_after_preproc - t_start) / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << "Closest point timing = " << ((double)(t_after_closest - t_after_preproc) / CLOCKS_PER_SEC) << " seconds" << endl;

    cout << "Test type = " << test_type << endl;
    cout << "Number of points tested for closestPoint = " << total_pts_tested << endl;

    cout << "Surf\tBest\tUnder" << endl;
    for (int i = 0; i < surfaces.size(); ++i)
      {
	string type_str = "Unknown";
	if (surfType[i] == 0)
	  type_str = "Spline";
	else if (surfType[i] == 1)
	  type_str = "Sphere";
	else if (surfType[i] == 2)
	  type_str = "Cylinder";
	cout << i << "\t" << isBest[i] << "\t" << type_str << endl;
      }

    cout << "find calls = " << find_calls << endl;
    cout << "bb calls = " << bb_calls << endl;
    cout << "cp calls = " << cp_calls << endl;
    cout << "Average tested length at find = " << (double)(find_len) / (double)(find_calls) << endl;
    cout << "Best point found on inside search " << best_inside << " times" << endl;

    int tot_bs = 0;
    for (int i = 0; i < boundaryCalls.size(); ++i)
      tot_bs += i * boundaryCalls[i];
    int tot_ul = 0;
    for (int i = 0; i < underlyingCalls.size(); ++i)
      tot_ul += i * underlyingCalls[i];
    cout << "ClosestPoint calls on BoundarySurface = " << tot_bs << endl;
    cout << "ClosestPoint calls on Underlying surface = " << tot_ul << endl;
    /*
    cout << endl << "N\tBS-call\tUL-call" << endl;
    int tot_bs = 0;
    int tot_ul = 0;
    for (int i = 0; i < underlyingCalls.size() || i < boundaryCalls.size(); ++i)
      {
	int bcalls = 0;
	if (i < boundaryCalls.size())
	  bcalls = boundaryCalls[i];
	int ucalls = 0;
	if (i < underlyingCalls.size())
	  ucalls = underlyingCalls[i];
	int tot_calls = bcalls + ucalls;
	if (tot_calls > 0)
	  {
	    cout << i << "\t";
	    if (bcalls > 0)
	      cout << bcalls;
	    cout << "\t";
	    if (ucalls > 0)
	      cout << ucalls;
	    cout << endl;
	  }
	tot_bs += i * bcalls;
	tot_ul += i * ucalls;
      }
    cout << "Total\t" << tot_bs << "\t" << tot_ul << endl;
    */

    /*
    cout << endl << "S-N\tE-N\tCalls" << endl;
    for (int i = 0; i < combinedCalls.size(); ++i)
      for (int j = 0; j < combinedCalls[i].size(); ++j)
	if (combinedCalls[i][j] > 0)
	  cout << i << "\t" << j << "\t" << combinedCalls[i][j] << endl;
    */

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
