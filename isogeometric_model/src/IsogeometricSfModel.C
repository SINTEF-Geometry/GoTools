//===========================================================================
//
// File : IsogeometricSfModel.C
//
// Created: Tue Mar  2 14:03:56 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#include "GoTools/isogeometric_model/IsogeometricSfModel.h"
#include "GoTools/isogeometric_model/BdCondFunctor.h"



using std::vector;
using std::pair;
using std::cerr;
using std::endl;
using std::shared_ptr;

namespace Go
{

  //===========================================================================
  IsogeometricSfModel::IsogeometricSfModel(shared_ptr<SurfaceModel> sfmodel,
					   vector<int> solution_space_dimension)
    : IsogeometricModel(sfmodel->getTolerances())
  //===========================================================================
  {
    // Test if all surfaces are spline surfaces with only corner-to-corner neighbour
    // relations. Note: Only the first neighbour relation of each pair of surfaces is
    // tested, so we need to redo the test later for each relation when building
    // the neighbouring structure

    if (!sfmodel->allSplines())
      {
	cerr << "Not all surfaces in surface model are splines" << endl;
	throw std::exception();
      }

    if (!sfmodel->isCornerToCorner())
      {
	cerr << "Not all neighbour relations between surfaces are corner to corner" << endl;
	throw std::exception();
      }

    // Build collection of blocks with no adjacency information

    double tol = getTolerances().neighbour;

    int nmb_blocks = sfmodel->nmbEntities();
    sf_blocks_.resize(nmb_blocks);
    for (int i = 0; i < nmb_blocks; ++i)
      sf_blocks_[i] = shared_ptr<IsogeometricSfBlock>(new IsogeometricSfBlock(this, sfmodel->getSplineSurface(i), solution_space_dimension, i));

    // Fetch adjacency information
    for (int i = 0; i < nmb_blocks-1; ++i)
      {
	shared_ptr<ftSurface> sf1 = sfmodel->getFace(i);
	for (int j = i+1; j < nmb_blocks; ++j)
	  {
	    shared_ptr<ftSurface> sf2 = sfmodel->getFace(j);
	    for (int idx = 0; true; ++idx)
	      {
		AdjacencyInfo adj_info = sf1->getAdjacencyInfo(sf2.get(), tol, idx, true);
		if (!adj_info.adjacency_found_)
		  break;
		if (adj_info.corner_failed_)
		  {
		    cerr << "Not all neighbour relations between surfaces are corner to corner" << endl;
		    throw std::exception();
		  }
		sf_blocks_[i]->addNeighbour(sf_blocks_[j], adj_info.bd_idx_1_, adj_info.bd_idx_2_, adj_info.same_orient_);
		sf_blocks_[j]->addNeighbour(sf_blocks_[i], adj_info.bd_idx_2_, adj_info.bd_idx_1_, adj_info.same_orient_);
	      }
	  }
      }

    buildBoundaryCurves(sfmodel);

  }

  //===========================================================================
  IsogeometricSfModel::~IsogeometricSfModel()
  //===========================================================================
  {
  }


  //===========================================================================
  bool IsogeometricSfModel::addBoundaryCond(int boundary, vector<Point> pos,
					    BdConditionType type, BdCondFunctor *fbd,
					    int solutionspace_idx, double *constant_value)
  //===========================================================================
  {
      if (boundary < 0 || boundary >= (int)boundary_curves_.size() ||
	solutionspace_idx < 0 || solutionspace_idx >= nmbSolutionSpaces())
      return false;

    if (pos.size() < 2)
      return false;

    if (type == CONSTANT_DIRICHLET && constant_value == 0 && fbd == 0)
      return false;

    double tol = getTolerances().neighbour;
    vector<shared_ptr<SplineCurve> > curves;
    vector<bool> is_degen;

    // Get spline curves along boundary and degenerate information
    for (int i = 0; i < boundary_curves_[boundary].size(); ++i)
      {
	shared_ptr<IsogeometricSfBlock> block = sf_blocks_[boundary_curve_block_[boundary][i]];
	int edge_nmb = boundary_curve_edge_[boundary][i] >> 1;
	curves.push_back(block->getGeomBoundaryCurve(edge_nmb));
	is_degen.push_back(curves[i]->isDegenerate(tol));
      }

    // Now we collect the segments determined by each pair of neighbouring points in the pos vector

    vector<int> segment_idx;    // For each edge segment, its boundary curve position in the vetor boundary_curve_[boundary]
    vector<double> local_par_start, local_par_end;    // For each edge segment, the local parameter position of the start and end points
    vector<int> local_start_pos, local_end_pos;    // For each edge segment, the position of the start and end point: 0 = start, 1 = end , 2 = interior
    vector<int> first_seg_idx;        // Candidates for segment_idx for first pair of points
    vector<double> first_par_start, first_par_end;   // Candidates for local_par_start and local_par_end for first pair of ponts
    vector<int> first_start_pos, first_end_pos;   // Candidates for local_start_pos and local_end_pos for first pair of ponts

    // First collect segments holding both of the first points
    for (int i = 0; i < (int)curves.size(); ++i)
      {
	double cl_par_0, cl_par_1, cl_dist;
	Point cl_pt;
	double startpar, endpar;
	startpar = curves[i]->startparam();
	endpar = curves[i]->endparam();
	curves[i]->closestPoint(pos[0], startpar, endpar, cl_par_0, cl_pt, cl_dist);
	if (cl_dist < tol)
	  {
	    curves[i]->closestPoint(pos[1], curves[i]->startparam(), curves[i]->endparam(), cl_par_1, cl_pt, cl_dist);
	    if (cl_dist < tol)
	      {
		first_seg_idx.push_back(i);
		Point p0, p1;
		curves[i]->point(p0, startpar);
		curves[i]->point(p1, endpar);
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
	    // double startpar = curves[first_seg_idx[0]]->startparam();
	    // double endpar = curves[first_seg_idx[0]]->endparam();
	    curves[first_seg_idx[0]]->closestPoint(pos[2],
						   curves[first_seg_idx[0]]->startparam(),
						   curves[first_seg_idx[0]]->endparam(),
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
    vector<int> pos_next(curves.size()), pos_prev(curves.size());
    if (pos_direction)
      {
	  for (int i = 0; i < (int)curves.size() - 1; ++i)
	  {
	    pos_next[i] = i + 1;
	    pos_prev[i + 1] = i;
	  }
	pos_next[curves.size() - 1] = 0;
	pos_prev[0] = (int)curves.size() - 1;
      }
    else
      {
	  for (int i = 0; i < (int)curves.size() - 1; ++i)
	  {
	    pos_next[i + 1] = i;
	    pos_prev[i] = i + 1;
	  }
	  pos_next[0] = (int)curves.size() - 1;
	pos_prev[curves.size() - 1] = 0;
      }

    for (int i = 1; i < (int)pos.size() - 1; ++i)
      {
	// Test if end point of previous segment has correct position. This is tested here,
	// as this is not required for the first or last point in pos[]
	if (((cur_pos_dir == pos_direction) && local_end_pos[i-1] != 1) ||
	    ((cur_pos_dir != pos_direction) && local_end_pos[i-1] != 0))
	  return false;

	// Find next segment, jump over degenerate curves
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
	startpar = curves[cur_segment_idx]->startparam();
	endpar = curves[cur_segment_idx]->endparam();

	curves[cur_segment_idx]->closestPoint(pos[i], startpar, endpar, cl_par_0, cl_pt, cl_dist);
	if (cl_dist >= tol)
	  return false;

	curves[cur_segment_idx]->closestPoint(pos[i+1], startpar, endpar, cl_par_1, cl_pt, cl_dist);
	if (cl_dist >= tol)
	  return false;

	Point p0, p1;
	curves[cur_segment_idx]->point(p0, startpar);
	curves[cur_segment_idx]->point(p1, endpar);

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

    // Start by finding the degenerate curves before the first pont if any.
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
	// Loop back and collect the degenerate curves
	while (true)
	  {
	    cur_segment_idx = pos_next[cur_segment_idx];
	    if (cur_segment_idx == segment_idx[0])
	      break;
	    final_idx.push_back(cur_segment_idx);
	    final_par.push_back(pair<double, double>(curves[cur_segment_idx]->startparam(), curves[cur_segment_idx]->endparam()));
	  }
      }

    // Now insert all previousely found segments together with degenerate segments between and at the end
    for (int i = 0; i < (int)segment_idx.size(); ++i)
      {
	// Insert found segments
	final_idx.push_back(segment_idx[i]);
	final_par.push_back(pair<double, double>(local_par_start[i], local_par_end[i]));

	// Insert degenerate curves that follow
	cur_pos_dir = (boundary_curve_edge_[boundary][segment_idx[i]] & 1) == 0;
	if (i < (int)segment_idx.size() - 1 ||     // Always include degenerate segments that follow if between to segments
	    ((cur_pos_dir == pos_direction) && local_start_pos[i] == 0) ||
	    ((cur_pos_dir != pos_direction) && local_start_pos[i] == 1))
	  // Loop and collect the degenerate curves
	  for (cur_segment_idx = pos_next[segment_idx[i]];
	       is_degen[cur_segment_idx];
	       cur_segment_idx = pos_next[cur_segment_idx])
	    {
	      final_idx.push_back(cur_segment_idx);
	      final_par.push_back(pair<double, double>(curves[cur_segment_idx]->startparam(), curves[cur_segment_idx]->endparam()));
	    }
      }

    // We have collected all segments. We now create the boundary conditions
    for (int i = 0; i < (int)final_idx.size(); ++i)
      {
	shared_ptr<IsogeometricSfBlock> block = sf_blocks_[boundary_curve_block_[boundary][final_idx[i]]];
	shared_ptr<SfSolution> sol = block->getSolutionSpace(solutionspace_idx);
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
  }


  //===========================================================================
  int IsogeometricSfModel::getNmbOfBoundaries() const
  //===========================================================================
  {
      return (int)boundary_curves_.size();
  }


  //===========================================================================
  CurveLoop IsogeometricSfModel::getOuterBoundary() const
  //===========================================================================
  {
    return getBoundary(0);
  }


  //===========================================================================
  CurveLoop IsogeometricSfModel::getBoundary(int idx) const
  //===========================================================================
  {
      if (idx >= 0 && idx < (int)boundary_curves_.size())
      return boundary_curves_[idx];
    else
      {
	CurveLoop cl_dummy;
	return cl_dummy;
      }
  }


  //===========================================================================
  void IsogeometricSfModel::setMinimumDegree(int degree, int solutionspace_idx)
  //===========================================================================
  {
      for (int i = 0; i < (int)sf_blocks_.size(); ++i)
      sf_blocks_[i]->setMinimumDegree(degree, solutionspace_idx);
  }


  //===========================================================================
  void IsogeometricSfModel::updateSolutionSplineSpace()
  //===========================================================================
  {
    int nmb_sol = nmbSolutionSpaces();
    for (int i = 0; i < nmb_sol; ++i)
      updateSolutionSplineSpace(i);
  }


  //===========================================================================
  void IsogeometricSfModel::updateSolutionSplineSpace(int solutionspace_idx)
  //===========================================================================
  {
    bool check = true;

    while (check)
      {
	check = false;
	for (int i = 0; i < (int)sf_blocks_.size(); ++i)
	  if (sf_blocks_[i]->updateSolutionSplineSpace(solutionspace_idx))
	    {
	      check = true;
	      break;
	    }
      }
  }


  //===========================================================================
  void IsogeometricSfModel::getIsogeometricBlocks(vector<shared_ptr<IsogeometricSfBlock> >& sfblock)
  //===========================================================================
  {
    sfblock.resize(sf_blocks_.size());
    copy(sf_blocks_.begin(), sf_blocks_.end(), sfblock.begin());
  }


  //===========================================================================
  void IsogeometricSfModel::makeGeometrySplineSpaceConsistent()
  //===========================================================================
  {
    MESSAGE("makeGeometrySplineSpaceConsistent() not implemented");
  }


  //===========================================================================
  void IsogeometricSfModel::buildBoundaryCurves(shared_ptr<SurfaceModel> sfmodel)
  //===========================================================================
  {
    double tol = sfmodel->getApproximationTol();

    int nmb_bnd = sfmodel->nmbBoundaries();
    boundary_curves_.resize(nmb_bnd);
    boundary_curve_block_.resize(nmb_bnd);
    boundary_curve_edge_.resize(nmb_bnd);

    int nmb_blocks = sfmodel->nmbEntities();

    for (int i = 0; i < nmb_bnd; ++i)
      {
	vector<shared_ptr<ftEdge> > edges = sfmodel->getBoundaryEdges(i);
	int nmb_seg = (int)edges.size();
	boundary_curve_block_[i].resize(nmb_seg);
	boundary_curve_edge_[i].resize(nmb_seg);
	vector<shared_ptr<ParamCurve> > loop_curves;
	for (int j = 0; j < nmb_seg; ++j)
	  {
	    shared_ptr<ParamCurve> curve = edges[j]->geomCurve();
	    loop_curves.push_back(curve);

	    ftFaceBase* face = edges[j]->face();
	    int surf_idx;
	    for (surf_idx = 0; surf_idx < nmb_blocks; ++surf_idx)
	      if (face == sfmodel->getFace(surf_idx).get())
		break;

	    if (surf_idx == nmb_blocks)
	      {
		cerr << "Could not find underlying surface of boundary " << i << ", segment " << j << endl;
		throw std::exception();
	      }

	    int edge_orient = sf_blocks_[surf_idx]->getEdgeOrientation(curve, tol);
	    if (edge_orient == -1)
	      {
		cerr << "Could not determine edge position on underlying surface of boundary " << i << ", segment " << j << endl;
		throw std::exception();
	      }

	    boundary_curve_block_[i][j] = surf_idx;
	    boundary_curve_edge_[i][j] = edge_orient;

	  }
	boundary_curves_[i] = CurveLoop(loop_curves, tol);
      }
  }


  //===========================================================================
  int IsogeometricSfModel::nmbSolutionSpaces() const
  //===========================================================================
  {
    if (sf_blocks_.empty())
      return 0;
    else
      return sf_blocks_[0]->nmbSolutionSpaces();
  }



} // end namespace Go
