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
  void IsogeometricVolModel::addBoundaryCond(vector<std::pair<ParamSurface*,
							      Point> > polygon,
					     BdConditionType type,
					     BdCondFunctor *fbd,
					     int solutionspace_idx,
					     double *constant_value)
  //===========================================================================
  {
    MESSAGE("addBoundaryCond() not implemented");
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
