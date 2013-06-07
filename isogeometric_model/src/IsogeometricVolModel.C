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

#include "GoTools/isogeometric_model/IsogeometricVolModel.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include <assert.h>

//#define TEMP_DEBUG   // Remove later when building volume code

using std::vector;
using std::cerr;
using std::endl;

namespace Go
{

  //===========================================================================
  IsogeometricVolModel::IsogeometricVolModel(shared_ptr<VolumeModel> volmodel,
					     vector<int> solution_space_dimension)
    : IsogeometricModel(volmodel->getTolerances())
  //===========================================================================
  {
    // Test if all volumes are spline volumes with only corner-to-corner neighbour
    // relations. Note: Only the first neighbour relation of each pair of volumes is
    // tested, so we need to redo the test later for each relation when building
    // the neighbouring structure

    if (!volmodel->allSplines())
      {
	cerr << "Not all volumes in volume model are splines" << endl;
	throw std::exception();
      }

    if (!volmodel->isCornerToCorner())
      {
	cerr << "Not all neighbour relations between volumes are corner to corner" << endl;
	throw std::exception();
      }

    // Build collection of blocks with no adjacency information

    double tol = volmodel->getTolerances().neighbour;

    int nmb_blocks = volmodel->nmbEntities();
    vol_blocks_.resize(nmb_blocks);
    for (int i = 0; i < nmb_blocks; ++i)
      vol_blocks_[i] = shared_ptr<IsogeometricVolBlock>
	(new IsogeometricVolBlock(this, volmodel->getSplineVolume(i),
				  solution_space_dimension, i));

    // Fetch adjacency information
    for (int i = 0; i < nmb_blocks; ++i)
      {
	shared_ptr<ftVolume> vol1 = volmodel->getBody(i);
	// for (int j = i+1; j < nmb_blocks; ++j)
	for (int j = 0; j < nmb_blocks; ++j) // getAdjacencyInfo() is not symmetric.
	  {
	    if (i == j)
	      {
		// Check for a closed volume
		shared_ptr<ParamVolume> parvol = vol1->getVolume();
		shared_ptr<SplineVolume> splvol =
		  dynamic_pointer_cast<SplineVolume, ParamVolume>(parvol);
		if (!splvol.get())
		  continue;  // Not a spline volume. Something is wrong
		for (int dir=0; dir<3; ++dir)
		  {
		    int per = splvol->volumePeriodicity(dir, getTolerances().gap);
		    if (per == 0)
		      {
			// Closed. Set adjacency info
			// par > 1 indicates a periodic volume. This is currently
			// not treated
			vol_blocks_[i]->addNeighbour(vol_blocks_[i], 2*per,
						    2*per+1, 0, true);
			vol_blocks_[i]->addNeighbour(vol_blocks_[i], 2*per+1,
						    2*per, 0, true);
		      }
		  }
		continue;   // No other adjacency for this case
	      }
			
		    
	    shared_ptr<ftVolume> vol2 = volmodel->getBody(j);
	    for (int idx = 0; true; ++idx)
	      {
// #ifndef TEMP_DEBUG
		VolumeAdjacencyInfo vol_adj_info =
		  vol1->getAdjacencyInfo(vol2.get(), tol, idx, true);
		bool reversed_const_dir =
		  (vol_adj_info.bd_idx_1_ + vol_adj_info.bd_idx_2_)%2 == 0;
// #else
// 		AdjacencyInfo vol_adj_info;
// #endif
		if (!vol_adj_info.adjacency_found_)
		  break;
		if (vol_adj_info.corner_failed_)
		  {
		    cerr << "Not all neighbour relations between volumes are "
		      "corner to corner" << endl;
		    throw std::exception();
		  }
// #ifndef TEMP_DEBUG
		// We compute the orientation of j block wrt i.
		bool reversed_sf_u = !vol_adj_info.same_orient_u_;
		bool reversed_sf_v = !vol_adj_info.same_orient_v_;
		// We must map 
		bool vol_u_rev, vol_v_rev, vol_w_rev;
		// @@sbr201111 Assuming that the surface is
		// constructed by keeping the order of the axes (in
		// the volume).
		if (vol_adj_info.bd_idx_1_ < 2)
		  {
		    vol_u_rev = reversed_const_dir;
		    vol_v_rev = reversed_sf_u;
		    vol_w_rev = reversed_sf_v;
		  }
		else if (vol_adj_info.bd_idx_1_ < 4)
		  {
		    vol_v_rev = reversed_const_dir;
		    vol_w_rev = reversed_sf_v;//u;
		    vol_u_rev = reversed_sf_u;//v;
		  }
		else
		  {
		    vol_w_rev = reversed_const_dir;
		    vol_u_rev = reversed_sf_u;
		    vol_v_rev = reversed_sf_v;
		  }
		int orientation_i = -1;
		if (!vol_u_rev && !vol_v_rev && !vol_w_rev)
		  orientation_i = 0;
		else if (vol_u_rev && !vol_v_rev && !vol_w_rev)
		  orientation_i = 1;
		else if (!vol_u_rev && vol_v_rev && !vol_w_rev)
		  orientation_i = 2;
		else if (!vol_u_rev && !vol_v_rev && vol_w_rev)
		  orientation_i = 3;
		if (vol_u_rev && vol_v_rev && !vol_w_rev)
		  orientation_i = 4;
		else if (vol_u_rev && !vol_v_rev && vol_w_rev)
		  orientation_i = 5;
		else if (!vol_u_rev && vol_v_rev && vol_w_rev)
		  orientation_i = 6;
		else if (vol_u_rev && vol_v_rev && vol_w_rev)
		  orientation_i = 7;
		vol_blocks_[i]->addNeighbour(vol_blocks_[j], vol_adj_info.bd_idx_1_,
					     vol_adj_info.bd_idx_2_,
					     orientation_i,
					     vol_adj_info.same_dir_order_);
		// @@sbr201111 No need to add for both blocks, i is treated in inner loop.
		// int orientation2 = -1;
		// vol_blocks_[j]->addNeighbour(vol_blocks_[i], vol_adj_info.bd_idx_2_,
		// 			     vol_adj_info.bd_idx_1_,
		// 			     orientation2,
		// 			     vol_adj_info.same_dir_order_);
// #endif
	      }
	  }
      }

    try
      {
	buildBoundaryFaces(volmodel);
      }
    catch (...)
      {
	MESSAGE("buildBoundaryFaces() failed.");
      }

  }

  //===========================================================================
  IsogeometricVolModel::~IsogeometricVolModel()
  //===========================================================================
  {
  }


  //===========================================================================
  bool IsogeometricVolModel::addBoundaryCond(vector<std::pair<ParamSurface*,
							      Point> > polygon,
					     BdConditionType type,
					     BdCondFunctor *fbd,
					     int solutionspace_idx,
					     double *constant_value)
  //===========================================================================
  {
    // MESSAGE("addBoundaryCond() under construction");

    if ((constant_value != NULL) && (*constant_value != 0.0))
      { // @@sbr201209 Not sure if we are doing it this way, we may
	// use a constant function instead.
	MESSAGE("Constant value different from 0.0 not yet implemented.");
      }

    // Areas extending over multiple blocks are expected to have been split already.

    // Currently we assume that the input polygon corresponds to (the
    // outer boundary of) a SplineSurface.
    // We also expect each surface side to be defined by 2 points only.

#ifndef NDEBUG
    std::cout << "Polygon node size: " << polygon.size() << std::endl;
    if (polygon.size() > 0)
	std::cout << "Polygon node dim: " << polygon[0].second.dimension() << std::endl;
#endif

    if (polygon.size() < 3 ||
	solutionspace_idx < 0 || solutionspace_idx >= nmbSolutionSpaces())
      return false;

    if (type == CONSTANT_DIRICHLET && constant_value == 0 && fbd == 0)
      return false;

    const double tol = getTolerances().neighbour;

    double dist_end_pts = polygon.front().second.dist(polygon.back().second);
    // The first and last points should really be exactly the same.
    assert(dist_end_pts < tol);

    ParamSurface* bd_sf = polygon[0].first;
    // We also expect the polygon to refer to only one boundary surface.
    for (size_t ki = 1; ki < polygon.size(); ++ki)
      assert(bd_sf == polygon[ki].first);

    // We check if the surface is closed. In that case we must give
    // special threatment to points on the seem.
    int vol_closed[3];
    vol_closed[0] = vol_closed[1] = vol_closed[2] = 0;
    shared_ptr<SplineVolume> spline_vol(NULL);
    if (bd_sf->instanceType() == Class_SurfaceOnVolume)
    {
	shared_ptr<ParamVolume> par_vol = (dynamic_cast<SurfaceOnVolume*>(bd_sf))->getVolume();
	if (par_vol->instanceType() == Class_SplineVolume)
	{
	    spline_vol = dynamic_pointer_cast<SplineVolume, ParamVolume>(par_vol);
	    for (int ki=0; ki<3; ++ki)
		vol_closed[ki] = spline_vol->volumePeriodicity(ki, tol);
	}
    }
    bool vol_is_closed = (vol_closed[0] || vol_closed[1] || vol_closed[2]);

    // if (vol_is_closed)
    // {
    // 	MESSAGE("Warning: Volume is closed, expect problems with boundary domains on the seem!");
    // }

    vector<shared_ptr<SplineSurface> > surfaces;
    vector<bool> is_degen;

#if 1
    // The volume version of this routine is easier as the match
    // between points (i.e. the domain) and the geometry (faces for
    // the volume case) has already been made.

    // @@sbr201303 We can use the SplineVolume to determine whether
    // there may be surface edges which are equal, and thus resulting
    // in problems for the closest point evaluation.  If this is the
    // case and two consecutive points are equal, we should make sure
    // that they are on different surface edges.


    int face_nmb = -1;
    vector<pair<double, double> > domain;    
    shared_ptr<IsogeometricVolBlock> block;
    // If volume is closed in one or two dirs we store alternative
    // values for the parameter pts.
    vector<vector<pair<double, double> > > domain_alt(polygon.size());
    for (size_t ki = 0; ki < polygon.size(); ++ki)
      {
	Point pol_pt = polygon[ki].second;
	assert(pol_pt.dimension() == 3);
	// We must locate the parameter value for the points.
	// Currently not handling degenerate points.
	ParamSurface* par_sf = polygon[ki].first;	
	double clo_u, clo_v, clo_dist;
	Point clo_pt;
	double epsgeo = 1e-06;
	par_sf->closestPoint(pol_pt, clo_u, clo_v, clo_pt, clo_dist, epsgeo);

	// If volume is closed we may need to check on the other side
	// of the domain to see if we have 2 solutions.
	if (vol_is_closed && (par_sf->instanceType() == Class_SurfaceOnVolume) && (spline_vol.get() != NULL))
	{
	    SurfaceOnVolume* sf_on_vol = dynamic_cast<SurfaceOnVolume*>(par_sf);

	    RectDomain cont_dom = sf_on_vol->containingDomain();
	    double umin = cont_dom.umin();
	    double umax = cont_dom.umax();
	    double vmin = cont_dom.vmin();
	    double vmax = cont_dom.vmax();

	    // In total there are 3 other par values to check, at most.
	    // So a brute force method is preferred to lots of logic ...
	    Point test_pt;
	    for (size_t kj = 0; kj < 2; ++kj)
	    {
		double vpar = (kj == 0) ? clo_v : ((fabs(clo_v-vmin) < fabs(vmax-clo_v)) ? vmax : vmin);
		for (size_t kk = 0; kk < 2; ++kk)
		{
		    double upar = (kk == 0) ? clo_u : ((fabs(clo_u-umin) < fabs(umax-clo_u)) ? umax : umin);
		    if (kj == 0 && kk == 0)
			continue; // This is clo_u & clo_v.
		    sf_on_vol->point(test_pt, upar, vpar);
		    double dist = pol_pt.dist(test_pt);
		    if (dist < epsgeo)
		    {
			// MESSAGE("We found another match!");
			// std::cout << "ki: " << ki << ", upar: " << upar << ", vpar: " << vpar << std::endl;
			domain_alt[ki].push_back(std::make_pair(upar, vpar));
		    }
		}
	    }
	}

	domain.push_back(std::make_pair(clo_u, clo_v));
	// We must locate the face_nmb.
	for (size_t kj = 0; kj < boundary_surfaces_.size(); ++kj)
	  {
	    for (size_t kk = 0; kk < boundary_surfaces_[kj].size(); ++kk)
	      {
		if (boundary_surfaces_[kj][kk].get() == par_sf)
		  {
		    if (ki == 0)
		      {
			face_nmb = kj;
			block = vol_blocks_[kj];
		      }
		    else
		      {
			if ((kj != face_nmb) || (block != vol_blocks_[kj]))
			  MESSAGE("Did not expect this, mismatch for face index and/or vol block.");
			break;
		      }
		  }
	      }
	  }
      }

    assert(block.get() != NULL);

    // In case the volume is closed we need some logic to choose
    // between par pts with the same geom pos.
    double pol_umin = -1.0;
    double pol_umax = -1.0;
    double pol_vmin = -1.0;
    double pol_vmax = -1.0;
    for (int ki = 0; ki < domain.size() - 1; ++ki)
    {
	if (ki == 0)
	{
	    pol_umin = domain[ki].first;
	    pol_umax = domain[ki].first;
	    pol_vmin = domain[ki].second;
	    pol_vmax = domain[ki].second;
	}
	else
	{
	    pol_umin = std::min(pol_umin, domain[ki].first);
	    pol_umax = std::max(pol_umin, domain[ki].first);
	    pol_vmin = std::min(pol_vmin, domain[ki].second);
	    pol_vmax = std::max(pol_vmax, domain[ki].second);
	    for (size_t kj = 0; kj < domain_alt[ki].size(); ++kj)
	    {
		pol_umin = std::min(pol_umin, domain_alt[ki][kj].first);
		pol_umax = std::max(pol_umin, domain_alt[ki][kj].first);
		pol_vmin = std::min(pol_vmin, domain_alt[ki][kj].second);
		pol_vmax = std::max(pol_vmax, domain_alt[ki][kj].second);
	    }
	}
    }

    // We now run through the domain, checking if there are two
    // consecutive points that are equal (which will typically be the
    // case for a closed volume).
    for (int ki = 0; ki < domain.size() - 1; ++ki)
    {
	std::pair<double, double> curr_par = domain[ki];
	std::pair<double, double> next_par = domain[ki+1];
	if ((curr_par.first == next_par.first) && (curr_par.second == next_par.second))
	{
	    if (domain_alt[ki].size() > 0 && domain_alt[ki+1].size() > 0)
	    {
		if ((domain_alt[ki].size() > 1) || (domain_alt[ki+1].size() > 1))
		    MESSAGE("Did not expect more than one candidate.");
		// When changing a value we make sure the polygon is
		// CCW. Furthermore we are assuming the polygon lies
		// along the outer edges of the surface.
		// This can occur for volumes which are closed in two
		// dirs, currently not handled.

		std::pair<double, double> alt_par = domain_alt[ki][0];

		int closed_dir = (curr_par.first - alt_par.first <
				  curr_par.second - alt_par.second) ? 0 : 1;
		// Here we are assuming that the polygon lies on a
		// surface boundary, otherwise we need a more robust
		// analysis of the loop.
		int ccw_opp_edge = (closed_dir == 0) ?
		    (fabs(domain[ki].second - pol_vmin) < fabs(pol_vmax - domain[ki].second) ? 0 : 2) :
		    (fabs(domain[ki].first - pol_umin) < fabs(pol_umax - domain[ki].first) ? 3 : 1);

		bool alt_is_higher = (alt_par.first + alt_par.second > curr_par.first + curr_par.second);

		if (ccw_opp_edge > 1)
		{
		    if (alt_is_higher)
			domain[ki] = alt_par;
		    else
			domain[ki+1] = alt_par;
		}
		else
		    if (alt_is_higher)
			domain[ki+1] = alt_par;
		    else
			domain[ki] = alt_par;
	    }
	    else
	    {
		MESSAGE("Failed handling problem with equal par pts.");
	    }

	}
	// else
	// {
	//     // @@sbr201303 We should also make sure that only one of the parameters is altered.
	//     ;
	// }

    }
    domain[domain.size()-1] = domain.front();

    // // We need to locate the vol block which the surface is part of.
    // // All sfs should be equal, we pick one.
    // ParamSurface* bd_sf = polygon[0].first;
    // SplineSurface* bd_sf_spl = dynamic_cast<SplineSurface*>(bd_sf);
    // assert(bd_sf_spline != NULL);
    // for (size_t ki = 0; ki < vol_blocks_.size(); ++ki)
    //   {
    // 	shared_ptr<SplineSurface> bd_sf_spline = vol_blocks_[ki]->getGeomBoundarySurface(face_nmb);
    // 	if (bd_sf_spl == bd_sf_spline.get())
    // 	  {
    // 	    sol = vol_blocks_[ki];
    // 	    break;
    // 	  }
    //   }

    shared_ptr<VolSolution> sol = block->getSolutionSpace(solutionspace_idx);

    sol->addBoundaryCondition(face_nmb, type, fbd, domain);
#else

    // Get spline surfaces along boundary and degenerate information
    for (int i = 0; i < boundary_surfaces_[boundary].size(); ++i)
      {
	shared_ptr<IsogeometricVolBlock> block = vol_blocks_[boundary_surface_block_[boundary][i]];
	int edge_nmb = boundary_curve_edge_[boundary][i] >> 1;
	surfaces.push_back(block->getGeomBoundaryCurve(edge_nmb));
	is_degen.push_back(surfaces[i]->isDegenerate(tol));
      }

    // Now we collect the segments determined by each pair of neighbouring points in the pos vector

    vector<int> segment_idx;    // For each edge segment, its boundary curve position in the vetor boundary_curve_[boundary]
    vector<double> local_par_start, local_par_end;    // For each edge segment, the local parameter position of the start and end points
    vector<int> local_start_pos, local_end_pos;    // For each edge segment, the position of the start and end point: 0 = start, 1 = end , 2 = interior
    vector<int> first_seg_idx;        // Candidates for segment_idx for first pair of points
    vector<double> first_par_start, first_par_end;   // Candidates for local_par_start and local_par_end for first pair of ponts
    vector<int> first_start_pos, first_end_pos;   // Candidates for local_start_pos and local_end_pos for first pair of ponts

    // First collect segments holding both of the first points
    for (int i = 0; i < (int)surfaces.size(); ++i)
      {
	double cl_par_0, cl_par_1, cl_dist;
	Point cl_pt;
	double startpar, endpar;
	startpar = surfaces[i]->startparam();
	endpar = surfaces[i]->endparam();
	surfaces[i]->closestPoint(pos[0], startpar, endpar, cl_par_0, cl_pt, cl_dist);
	if (cl_dist < tol)
	  {
	    surfaces[i]->closestPoint(pos[1], surfaces[i]->startparam(), surfaces[i]->endparam(), cl_par_1, cl_pt, cl_dist);
	    if (cl_dist < tol)
	      {
		first_seg_idx.push_back(i);
		Point p0, p1;
		surfaces[i]->point(p0, startpar);
		surfaces[i]->point(p1, endpar);
		if (p0.dist(pos[0]) < tol)
		  cl_par_0 = startpar;
		else if (p1.dist(pos[0]) < tol)
		  cl_par_0 = endpar;
		if (p0.dist(pos[1]) < tol)
		  cl_par_1 = startpar;
		else if (p1.dist(pos[1]) < tol)
		  cl_par_1 = endpar;
		first_par_start.push_back(cl_par_0);
		first_par_end.push_back(cl_par_1);

		if (cl_par_0 == startpar)
		  first_start_pos.push_back(0);
		else if (cl_par_0 == endpar)
		  first_start_pos.push_back(1);
		else
		  first_start_pos.push_back(2);

		if (cl_par_1 == startpar)
		  first_end_pos.push_back(0);
		else if (cl_par_1 == endpar)
		  first_end_pos.push_back(1);
		else
		  first_end_pos.push_back(2);
	      }
	  }
      }

    // Determine the starting edge
    int nmb_first = (int)first_par_start.size();
    if (nmb_first == 0)
      return false;
    int first_pos;
    if (nmb_first == 1)
      first_pos = 0;
    else
      {
	// Two candidates. Which one to choose depends on the position of the third point
	// if it exists, otherwise choose positive orientation
	if (pos.size() >= 3)
	  {
	    double cl_par, cl_dist;
	    Point cl_pt;
	    // double startpar = surfaces[first_seg_idx[0]]->startparam();
	    // double endpar = surfaces[first_seg_idx[0]]->endparam();
	    surfaces[first_seg_idx[0]]->closestPoint(pos[2],
						   surfaces[first_seg_idx[0]]->startparam(),
						   surfaces[first_seg_idx[0]]->endparam(),
						   cl_par, cl_pt, cl_dist);
	    if (cl_dist < tol)
	      first_pos = 1;
	    else
	      first_pos = 0;
	  }
	else
	  {
	    if (((boundary_curve_edge_[boundary][first_seg_idx[0]] & 1) == 0)  // first common segment has same orientation along surface and boundary
		== (first_par_start[0] < first_par_end[0]))  // pos[0]->pos[1] along first common segment is same as surface orientation
	      first_pos = 0;
	    else
	      first_pos = 1;
	  }
      }

    // Insert first segment
    segment_idx.push_back(first_seg_idx[first_pos]);
    local_par_start.push_back(first_par_start[first_pos]);
    local_par_end.push_back(first_par_end[first_pos]);
    local_start_pos.push_back(first_start_pos[first_pos]);
    local_end_pos.push_back(first_end_pos[first_pos]);

    // Insert remaining segments
    int cur_segment_idx = segment_idx[0];
    bool cur_pos_dir = (boundary_curve_edge_[boundary][cur_segment_idx] & 1) == 0;
    bool pos_direction = cur_pos_dir == (local_par_start[0] < local_par_end[0]);

    // Two arrays to simplify traversing
    vector<int> pos_next(surfaces.size()), pos_prev(surfaces.size());
    if (pos_direction)
      {
	  for (int i = 0; i < (int)surfaces.size() - 1; ++i)
	  {
	    pos_next[i] = i + 1;
	    pos_prev[i + 1] = i;
	  }
	pos_next[surfaces.size() - 1] = 0;
	pos_prev[0] = (int)surfaces.size() - 1;
      }
    else
      {
	  for (int i = 0; i < (int)surfaces.size() - 1; ++i)
	  {
	    pos_next[i + 1] = i;
	    pos_prev[i] = i + 1;
	  }
	  pos_next[0] = (int)surfaces.size() - 1;
	pos_prev[surfaces.size() - 1] = 0;
      }

    for (int i = 1; i < (int)pos.size() - 1; ++i)
      {
	// Test if end point of previous segment has correct position. This is tested here,
	// as this is not required for the first or last point in pos[]
	if (((cur_pos_dir == pos_direction) && local_end_pos[i-1] != 1) ||
	    ((cur_pos_dir != pos_direction) && local_end_pos[i-1] != 0))
	  return false;

	// Find next segment, jump over degenerate surfaces
	bool step = true;
	while (step)
	  {
	    cur_segment_idx = pos_next[cur_segment_idx];
	    step = is_degen[cur_segment_idx];
	  }

	cur_pos_dir = (boundary_curve_edge_[boundary][cur_segment_idx] & 1) == 0;

	// Get start and end point information
	double cl_par_0, cl_par_1, cl_dist;
	Point cl_pt;
	double startpar, endpar;
	startpar = surfaces[cur_segment_idx]->startparam();
	endpar = surfaces[cur_segment_idx]->endparam();

	surfaces[cur_segment_idx]->closestPoint(pos[i], startpar, endpar, cl_par_0, cl_pt, cl_dist);
	if (cl_dist >= tol)
	  return false;

	surfaces[cur_segment_idx]->closestPoint(pos[i+1], startpar, endpar, cl_par_1, cl_pt, cl_dist);
	if (cl_dist >= tol)
	  return false;

	Point p0, p1;
	surfaces[cur_segment_idx]->point(p0, startpar);
	surfaces[cur_segment_idx]->point(p1, endpar);

	if (p0.dist(pos[i]) < tol)
	  cl_par_0 = startpar;
	else if (p1.dist(pos[i]) < tol)
	  cl_par_0 = endpar;
	if (p0.dist(pos[i+1]) < tol)
	  cl_par_1 = startpar;
	else if (p1.dist(pos[i+1]) < tol)
	  cl_par_1 = endpar;

	local_par_start.push_back(cl_par_0);
	local_par_end.push_back(cl_par_1);

	if (cl_par_0 == startpar)
	  local_start_pos.push_back(0);
	else if (cl_par_0 == endpar)
	  local_start_pos.push_back(1);
	else
	  local_start_pos.push_back(2);

	if (cl_par_1 == startpar)
	  local_end_pos.push_back(0);
	else if (cl_par_1 == endpar)
	  local_end_pos.push_back(1);
	else
	  local_end_pos.push_back(2);

	segment_idx.push_back(cur_segment_idx);

	// Last, test if start point of this segment has correct position.
	if (((cur_pos_dir == pos_direction) && local_start_pos[i] != 0) ||
	    ((cur_pos_dir != pos_direction) && local_start_pos[i] != 1))
	  return false;
      }

    // All the segments are collected, and the points are in good positions. Now we can build the boundary conditions
    vector<int> final_idx;
    vector<pair<double, double> > final_par;

    // Start by finding the degenerate surfaces before the first pont if any.
    // Set cur_segment_idx to be the position of the first non-degenerate curve before segment_idx[0]
    // if pos[0] is on the start corner of its segment.
    cur_segment_idx = segment_idx[0];
    cur_pos_dir = (boundary_curve_edge_[boundary][cur_segment_idx] & 1) == 0;
    if (((cur_pos_dir == pos_direction) && local_start_pos[0] == 0) ||
	((cur_pos_dir != pos_direction) && local_start_pos[0] == 1))
      {
	// We have a start corner. Jump back to previous nondegenerate curve
	while (true)
	  {
	    cur_segment_idx = pos_prev[cur_segment_idx];
	    if (!is_degen[cur_segment_idx])
	      break;
	  }
	// Loop back and collect the degenerate surfaces
	while (true)
	  {
	    cur_segment_idx = pos_next[cur_segment_idx];
	    if (cur_segment_idx == segment_idx[0])
	      break;
	    final_idx.push_back(cur_segment_idx);
	    final_par.push_back(pair<double, double>(surfaces[cur_segment_idx]->startparam(), surfaces[cur_segment_idx]->endparam()));
	  }
      }

    // Now insert all previousely found segments together with degenerate segments between and at the end
    for (int i = 0; i < (int)segment_idx.size(); ++i)
      {
	// Insert found segments
	final_idx.push_back(segment_idx[i]);
	final_par.push_back(pair<double, double>(local_par_start[i], local_par_end[i]));

	// Insert degenerate surfaces that follow
	cur_pos_dir = (boundary_curve_edge_[boundary][segment_idx[i]] & 1) == 0;
	if (i < (int)segment_idx.size() - 1 ||     // Always include degenerate segments that follow if between to segments
	    ((cur_pos_dir == pos_direction) && local_start_pos[i] == 0) ||
	    ((cur_pos_dir != pos_direction) && local_start_pos[i] == 1))
	  // Loop and collect the degenerate surfaces
	  for (cur_segment_idx = pos_next[segment_idx[i]];
	       is_degen[cur_segment_idx];
	       cur_segment_idx = pos_next[cur_segment_idx])
	    {
	      final_idx.push_back(cur_segment_idx);
	      final_par.push_back(pair<double, double>(surfaces[cur_segment_idx]->startparam(), surfaces[cur_segment_idx]->endparam()));
	    }
      }

    // We have collected all segments. We now create the boundary conditions
    for (int i = 0; i < (int)final_idx.size(); ++i)
      {
	shared_ptr<IsogeometricVolBlock> block = vol_blocks_[boundary_surface_block_[boundary][final_idx[i]]];
	shared_ptr<VolSolution> sol = block->getSolutionSpace(solutionspace_idx);
	if (type == CONSTANT_DIRICHLET || type == ZERO_DIRICHLET)
	  {
	    Point pt(sol->dimension());
	    if (type == ZERO_DIRICHLET)
	      for (int j = 0; j < sol->dimension(); ++j)
		pt[j] = 0.0;
	    else   // CONSTANT_DIRICHELT
	      {
		if (constant_value != 0)
		  for (int j = 0; j < sol->dimension(); ++j)
		    pt[j] = constant_value[j];
		else
		  pt = fbd->evaluate(pos[0]);
	      }
	    sol->addBoundaryCondition(boundary_curve_edge_[boundary][segment_idx[i]] >> 1, type, pt, final_par[i]);
	  }
	else
	  sol->addBoundaryCondition(boundary_curve_edge_[boundary][segment_idx[i]] >> 1, type, fbd, final_par[i]);
      }
    return true;

#endif
  }


  //===========================================================================
  void IsogeometricVolModel::addDirichletPointBdCond(ParamSurface* bd_surf,
						     Point& pos,
						     Point& condition_value,
						     int solutionspace_idx)
  //===========================================================================
  {
    MESSAGE("addDirichletPointBdCond() not implemented");
  }


  //===========================================================================
  int IsogeometricVolModel::getNmbOfBoundaries() const
  //===========================================================================
  {
      return (int)boundary_surfaces_.size();
  }


  //===========================================================================
  vector<shared_ptr<ParamSurface> > IsogeometricVolModel::getOuterBoundary() const
  //===========================================================================
  {
    return getBoundary(0);
  }


  //===========================================================================
  vector<shared_ptr<ParamSurface> > IsogeometricVolModel::getBoundary(int idx) const
  //===========================================================================
  {
      if (idx >= 0 && idx < (int)boundary_surfaces_.size())
      return boundary_surfaces_[idx];
    else
      {
	vector<shared_ptr<ParamSurface> > sm_dummy;
	return sm_dummy;
      }
  }


  //===========================================================================
  void IsogeometricVolModel::setMinimumDegree(int degree, int solutionspace_idx)
  //===========================================================================
  {
    for (int i = 0; i < (int)vol_blocks_.size(); ++i)
      vol_blocks_[i]->setMinimumDegree(degree, solutionspace_idx);
  }


  //===========================================================================
  void IsogeometricVolModel::updateSolutionSplineSpace()
  //===========================================================================
  {
    int nmb_sol = nmbSolutionSpaces();
    for (int i = 0; i < nmb_sol; ++i)
      updateSolutionSplineSpace(i);
  }


  //===========================================================================
  void IsogeometricVolModel::updateSolutionSplineSpace(int solutionspace_idx)
  //===========================================================================
  {
    bool check = true;

    while (check)
      {
	check = false;
	for (int i = 0; i < (int)vol_blocks_.size(); ++i)
	  if (vol_blocks_[i]->updateSolutionSplineSpace(solutionspace_idx))
	    {
	      check = true;
	      break;
	    }
      }
  }

  //===========================================================================
  void
  IsogeometricVolModel::getIsogeometricBlocks(vector<shared_ptr<IsogeometricVolBlock> >& volblock)
  //===========================================================================
  {
    volblock.resize(vol_blocks_.size());
    copy(vol_blocks_.begin(), vol_blocks_.end(), volblock.begin());
  }


  //===========================================================================
  void IsogeometricVolModel::makeGeometrySplineSpaceConsistent()
  //===========================================================================
  {
    MESSAGE("makeGeometrySplineSpaceConsistent() not implemented");
  }


  //===========================================================================
  void IsogeometricVolModel::buildBoundaryFaces(shared_ptr<VolumeModel> volmodel)
  //===========================================================================
  {
    double tol = volmodel->getApproximationTol();

    int nmb_bnd = volmodel->nmbBoundaries();
    boundary_surfaces_.resize(nmb_bnd);
    boundary_surface_block_.resize(nmb_bnd);
    boundary_surface_pos_.resize(nmb_bnd);

    int nmb_blocks = volmodel->nmbEntities(); // The number of volumes in the VolumeModel.

    for (int i = 0; i < nmb_bnd; ++i)
      {
	vector<shared_ptr<ftSurface> > faces = volmodel->getBoundaryFaces(i);
	int nmb_seg = (int)faces.size();
	boundary_surface_block_[i].resize(nmb_seg);
	// boundary_curve_edge_[i].resize(nmb_seg);
	boundary_surface_pos_[i].resize(nmb_seg);
	vector<shared_ptr<ParamSurface> > shells;
	for (int j = 0; j < nmb_seg; ++j)
	  {
	    shared_ptr<ParamSurface> surface = faces[j]->surface();
	    shells.push_back(surface);

	    // The body contains connected volumes.
	    Body* body = faces[j]->getBody();
	    ftVolume* ft_vol = dynamic_cast<ftVolume*>(body);
	    shared_ptr<ParamVolume> vol = ft_vol->getVolume();
	    int vol_idx = -1;
	    for (vol_idx = 0; vol_idx < nmb_blocks; ++vol_idx)
	      if (vol.get() == volmodel->getVolume(vol_idx).get())
		break;

	    if (vol_idx == nmb_blocks)
	      {
		cerr << "Could not find underlying volume of boundary " << i << ", segment " << j << endl;
		throw std::exception();
	      }

	    int face_orient = vol_blocks_[vol_idx]->getFaceOrientation(surface, tol);
	    if (face_orient == -1)
	      {
	    	cerr << "Could not determine face position on underlying surface of boundary " <<
	    	  i << ", segment " << j << endl;
	    	throw std::exception();
	      }

	    boundary_surface_block_[i][j] = vol_idx;
	    boundary_surface_pos_[i][j] = face_orient;

	  }
	boundary_surfaces_[i] = shells;
      }
  }


  //===========================================================================
  int IsogeometricVolModel::nmbSolutionSpaces() const
  //===========================================================================
  {
    if (vol_blocks_.empty())
      return 0;
    else
      return vol_blocks_[0]->nmbSolutionSpaces();
  }



} // end namespace Go
