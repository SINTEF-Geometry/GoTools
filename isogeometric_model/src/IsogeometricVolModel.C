//===========================================================================
//
// File : IsogeometricVolModel.C
//
// Created: Tue Mar  9 15:06:56 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================





#include "GoTools/isogeometric_model/IsogeometricVolModel.h"
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
    for (int i = 0; i < nmb_blocks-1; ++i)
      {
	shared_ptr<ftVolume> vol1 = volmodel->getBody(i);
	// for (int j = i+1; j < nmb_blocks; ++j)
	for (int j = 0; j < nmb_blocks; ++j) // getAdjacencyInfo() is not symmetric.
	  {
	    if (i == j)
	      continue;
	    shared_ptr<ftVolume> vol2 = volmodel->getBody(j);
	    for (int idx = 0; true; ++idx)
	      {
// #ifndef TEMP_DEBUG
		VolumeAdjacencyInfo vol_adj_info =
		  vol1->getAdjacencyInfo(vol2.get(), tol, idx, true);
		bool reversed_const_dir =
		  (vol_adj_info.bd_idx_1_ + vol_adj_info.bd_idx_2_)%2 == 1;
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
		bool reversed_sf_u = vol_adj_info.same_orient_u_;
		bool reversed_sf_v = vol_adj_info.same_orient_v_;
		// We must map 
		bool vol_u_rev, vol_v_rev, vol_w_rev;
		// @@sbr201111 Assuming that volume the surface is
		// constructed by keeping the order of the axes.
		if (vol_adj_info.bd_idx_1_ < 2)
		  {
		    vol_u_rev = reversed_const_dir;
		    vol_v_rev = reversed_sf_u;
		    vol_w_rev = reversed_sf_v;
		  }
		else if (vol_adj_info.bd_idx_1_ < 4)
		  {
		    vol_v_rev = reversed_const_dir;
		    vol_w_rev = reversed_sf_u;
		    vol_u_rev = reversed_sf_v;
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
		vol_blocks_[j]->addNeighbour(vol_blocks_[i], vol_adj_info.bd_idx_1_,
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
    MESSAGE("addBoundaryCond() under construction");

    // Currently we assume that the input polygon corresponds to (the
    // outer boundary of) a SplineSurface.
    // We also expect each surface side to be defined by 2 points only.

#ifndef NDEBUG
    std::cout << "Polygon node size: " << polygon.size() << std::endl;
#endif

    if (polygon.size() < 3 ||
	solutionspace_idx < 0 || solutionspace_idx >= nmbSolutionSpaces())
      return false;

    if (type == CONSTANT_DIRICHLET && constant_value == 0 && fbd == 0)
      return false;

#if 0 // These assertions do not seem to hold ...
    assert(polygon.size() == 5); // First and last points should be equal.
    ParamSurface* bd_sf = polygon[0].first;
    // We also expect the polygon to refer to only one boundary surface.
    for (size_t ki = 1; ki < polygon.size(); ++ki)
      assert(bd_sf == polygon[ki].first);
#endif

    double tol = getTolerances().neighbour;
    vector<shared_ptr<SplineSurface> > surfaces;
    vector<bool> is_degen;

#if 0

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

    // int nmb_blocks = volmodel->nmbEntities();

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

	    // int surf_idx;
	    // for (surf_idx = 0; surf_idx < nmb_blocks; ++surf_idx)
	    //   if (faces[j] == volmodel->getFace(surf_idx))
	    // 	break;

	    // if (surf_idx == nmb_blocks)
	    //   {
	    // 	cerr << "Could not find underlying surface of boundary " << i << ", segment " << j << endl;
	    // 	throw std::exception();
	    //   }

	    int face_orient = vol_blocks_[i]->getFaceOrientation(surface, tol);
	    if (face_orient == -1)
	      {
	    	cerr << "Could not determine face position on underlying surface of boundary " <<
	    	  i << ", segment " << j << endl;
	    	throw std::exception();
	      }

	    boundary_surface_block_[i][j] = i;//surf_idx;
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
