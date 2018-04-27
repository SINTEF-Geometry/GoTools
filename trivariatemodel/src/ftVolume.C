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

//#define DEBUG_VOL1
//#define DEBUG

#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/trivariate/SurfaceOnVolumeTools.h"
#include "GoTools/trivariate/VolumeTools.h"
#include "GoTools/compositemodel/SISLCurveInterface.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/compositemodel/AdaptSurface.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include "GoTools/compositemodel/ModifyFaceSet.h"
#include "GoTools/compositemodel/CompleteEdgeNet.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/compositemodel/Path.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/HermiteInterpolator.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/trivariate/CoonsPatchVolumeGen.h"
#include "GoTools/trivariate/VolumeParameterCurve.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/creators/HermiteAppC.h"
#include "GoTools/creators/SmoothSurf.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/topology/FaceConnectivityUtils.h"
#include <fstream>

using std::vector;
using std::set;
using std::make_pair;
using std::pair;
using std::min;
using std::max;

using namespace Go;

//---------------------------------------------------------------------------
ftVolume::ftVolume(shared_ptr<ParamVolume> vol, int id)
  : Body(), vol_(vol), id_(id)
//---------------------------------------------------------------------------
{
  double eps = 1.0e-6;
  double tang_eps = 1.0e-2;

  shared_ptr<SurfaceModel> shell = createBoundaryShell(eps, tang_eps);
  addshell(shell);
  //toptol_ = tpTolerances(eps, 10.0*eps, tang_eps, 10.0*tang_eps);
  //toptol_ = tpTolerances(eps, 200.0*eps, tang_eps, 10.0*tang_eps);
  toptol_ = tpTolerances(eps, 50.0*eps, tang_eps, 10.0*tang_eps);
			     
}

//---------------------------------------------------------------------------
ftVolume::ftVolume(shared_ptr<ParamVolume> vol, double gap_eps,
		   double kink_eps, int id)
  : Body(), vol_(vol), id_(id)
//---------------------------------------------------------------------------
{
  shared_ptr<SurfaceModel> shell = createBoundaryShell(gap_eps, kink_eps);

  addshell(shell);
  toptol_ = tpTolerances(gap_eps, 10.0*gap_eps, kink_eps, 10.0*kink_eps);
			     
}

//---------------------------------------------------------------------------
ftVolume::ftVolume(shared_ptr<ParamVolume> vol, double gap_eps, 
		   double neighbour, double kink_eps, double bend, int id)
  : Body(), vol_(vol), id_(id)
//---------------------------------------------------------------------------
{
  shared_ptr<SurfaceModel> shell = createBoundaryShell(gap_eps, kink_eps);

  addshell(shell);
  toptol_ = tpTolerances(gap_eps, neighbour, kink_eps, bend);
			     
}

//---------------------------------------------------------------------------
ftVolume::ftVolume(shared_ptr<ParamVolume> vol, 
		   shared_ptr<SurfaceModel> shell,
		   int id)
  : Body(shell), vol_(vol), id_(id)
//---------------------------------------------------------------------------
{
			     
}

//---------------------------------------------------------------------------
ftVolume::ftVolume(shared_ptr<ParamVolume> vol, 
		   vector<shared_ptr<SurfaceModel> > shells,
		   int id)
  : Body(shells), vol_(vol), id_(id)
//---------------------------------------------------------------------------
{

}

//---------------------------------------------------------------------------
ftVolume::ftVolume(shared_ptr<SurfaceModel> shell,
		   int id)
  : Body(shell), id_(id)
//---------------------------------------------------------------------------
{
  // Create a large enough volume
  // First make bounding box around the outer boundaries
  BoundingBox box = shell->boundingBox();
  Point low = box.low();
  Point high = box.high();

  // Make a trilinear spline volume
  int dim = 3;
  vector<double> knots1(4);
  knots1[0] = knots1[1] = 0.0;
  knots1[2] = knots1[3] = high[0] - low[0];
  vector<double> knots2(4);
  knots2[0] = knots2[1] = 0.0;
  knots2[2] = knots2[3] = high[1] - low[1];
  vector<double> knots3(4);
  knots3[0] = knots3[1] = 0.0;
  knots3[2] = knots3[3] = high[2] - low[2];
  vector<double> coefs(24);
  coefs[0] = coefs[2*dim] = coefs[4*dim] = coefs[6*dim] = low[0];
  coefs[dim] = coefs[3*dim] = coefs[5*dim] = coefs[7*dim] = high[0];
  coefs[1] = coefs[dim+1] = coefs[4*dim+1] = coefs[5*dim+1] = low[1];
  coefs[2*dim+1] = coefs[3*dim+1] = coefs[6*dim+1] = coefs[7*dim+1] = high[1];
  coefs[2] = coefs[dim+2] = coefs[2*dim+2] = coefs[3*dim+2] = low[2];
  coefs[4*dim+2] = coefs[5*dim+2] = coefs[6*dim+2] = coefs[7*dim+2] = high[2];

  vol_ = shared_ptr<ParamVolume>(new SplineVolume(2, 2, 2, 2, 2, 2,
						  knots1.begin(), knots2.begin(),
						  knots3.begin(), coefs.begin(),
						  dim));
}

//---------------------------------------------------------------------------
ftVolume::ftVolume(shared_ptr<Body> body,
		   int id)
  : Body(body), id_(id)
//---------------------------------------------------------------------------
{
  // Create a large enough volume
  // First make bounding box around the outer boundaries
  BoundingBox box = body->boundingBox();
  Point low = box.low();
  Point high = box.high();

  // Make a trilinear spline volume
  int dim = 3;
  vector<double> knots1(4);
  knots1[0] = knots1[1] = 0.0;
  knots1[2] = knots1[3] = high[0] - low[0];
  vector<double> knots2(4);
  knots2[0] = knots2[1] = 0.0;
  knots2[2] = knots2[3] = high[1] - low[1];
  vector<double> knots3(4);
  knots3[0] = knots3[1] = 0.0;
  knots3[2] = knots3[3] = high[2] - low[2];
  vector<double> coefs(24);
  coefs[0] = coefs[2*dim] = coefs[4*dim] = coefs[6*dim] = low[0];
  coefs[dim] = coefs[3*dim] = coefs[5*dim] = coefs[7*dim] = high[0];
  coefs[1] = coefs[dim+1] = coefs[4*dim+1] = coefs[5*dim+1] = low[1];
  coefs[2*dim+1] = coefs[3*dim+1] = coefs[6*dim+1] = coefs[7*dim+1] = high[1];
  coefs[2] = coefs[dim+2] = coefs[2*dim+2] = coefs[3*dim+2] = low[2];
  coefs[4*dim+2] = coefs[5*dim+2] = coefs[6*dim+2] = coefs[7*dim+2] = high[2];

  vol_ = shared_ptr<ParamVolume>(new SplineVolume(2, 2, 2, 2, 2, 2,
						  knots1.begin(), knots2.begin(),
						  knots3.begin(), coefs.begin(),
						  dim));
}

//---------------------------------------------------------------------------
shared_ptr<SurfaceModel> ftVolume::createBoundaryShell(double eps, double tang_eps)
//---------------------------------------------------------------------------
{
  vector<shared_ptr<ftSurface> > faces = getBoundaryFaces(vol_, eps, tang_eps);

  // Defining tolerances for surface model. This is not a good solution
//   shared_ptr<SurfaceModel> shell = 
//     shared_ptr<SurfaceModel>(new SurfaceModel(eps, eps, 200.0*eps,
// 					      tang_eps, 10.0*tang_eps,
// 					      faces));
  shared_ptr<SurfaceModel> shell = 
    shared_ptr<SurfaceModel>(new SurfaceModel(eps, eps, 50.0*eps,
					      tang_eps, 10.0*tang_eps,
					      faces));

  return shell;
}

//---------------------------------------------------------------------------
vector<shared_ptr<ftSurface> >
ftVolume::getBoundaryFaces(shared_ptr<ParamVolume> vol,
			   double eps, double tang_eps)
//---------------------------------------------------------------------------
{
  vector<shared_ptr<ParamSurface> > bd_sfs = VolumeTools::getOrientedBoundarySurfaces(vol);

  // Remove surfaces that degenerates to a line or to a point
  bool bottom, right, top, left;
  for (size_t ki=0; ki<bd_sfs.size();)
    {
      bool isdegen = bd_sfs[ki]->isDegenerate(bottom, right, top, left, eps /*10.0*eps*/);
      if (isdegen)
	{
	  RectDomain dom = bd_sfs[ki]->containingDomain();
	  int nmb_samples = 10;
	  int kj;
	  if (bottom && top)
	    {
	      // Check distance between other boundary curves
	      double t1 = dom.vmin();
	      double t2 = dom.vmax();
	      double tdel = (t2 - t1)/(double)(nmb_samples + 1);
	      double tpar;
	      for (kj=0, tpar=t1+tdel; kj<nmb_samples; ++kj, tpar+=tdel)
		{
		  vector<shared_ptr<ParamCurve> > constcrvs =
		    bd_sfs[ki]->constParamCurves(tpar, true);
		  double len = 0.0;
		  for (size_t kr=0; kr<constcrvs.size(); ++kr)
		    len += constcrvs[kr]->estimatedCurveLength();
		  if (len >= 10.0*eps)
		    break;  // Distance between opposite boundaries are not small
		}
	      if (kj >= nmb_samples)
		// Degenerate surface
		bd_sfs.erase(bd_sfs.begin() + ki);
	      else 
		ki++;
	    }
	  else if (right && left)
	    {
	      // Check distance between other boundary curves
	      double t1 = dom.umin();
	      double t2 = dom.umax();
	      double tdel = (t2 - t1)/(double)(nmb_samples + 1);
	      double tpar;
	      for (kj=0, tpar=t1+tdel; kj<nmb_samples; ++kj, tpar+=tdel)
		{
		  vector<shared_ptr<ParamCurve> > constcrvs =
		    bd_sfs[ki]->constParamCurves(tpar, false);
		  double len = 0.0;
		  for (size_t kr=0; kr<constcrvs.size(); ++kr)
		    len += constcrvs[kr]->estimatedCurveLength();
		  if (len >= 10.0*eps)
		    break;  // Distance between opposite boundaries are not small
		}
	      if (kj >= nmb_samples)
		// Degenerate surface
		bd_sfs.erase(bd_sfs.begin() + ki);
	      else 
		ki++;
	    }
	  else
	    ki++;
	}
      else ki++;
    }

  vector<shared_ptr<ftSurface> > faces(bd_sfs.size());
  for (size_t ki=0; ki<bd_sfs.size(); ++ki)
    {
      shared_ptr<ftSurface> curr_face = 
	shared_ptr<ftSurface>(new ftSurface(bd_sfs[ki], -1));
      shared_ptr<SurfaceOnVolume> curr_sf =
	dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bd_sfs[ki]);
      if (curr_sf.get())
	{
	  // Define the face loop explicitely to get correct pointers
	  // in the CurveOnSurface entities
	  // Do not split edges in internal corners of boundary curves
	  CurveLoop crv_loop = SurfaceOnVolumeTools::getOuterBoundaryLoop(curr_sf,
									  eps);
	  shared_ptr<Loop> loop = shared_ptr<Loop>(new Loop(curr_face.get(),
							    crv_loop, tang_eps,
							    false));
	  curr_face->addOuterBoundaryLoop(loop);
	}
      faces[ki] = curr_face;
    }

  return faces;
}

//---------------------------------------------------------------------------
ftVolume::~ftVolume()
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
BoundingBox ftVolume::boundingBox() const
//---------------------------------------------------------------------------
{
  return vol_->boundingBox();
}

//---------------------------------------------------------------------------
void ftVolume::closestPoint(Point& pt,
			    double&  clo_u,
			    double&  clo_v, 
			    double&  clo_w, 
			    Point& clo_pt,
			    double&  clo_dist,
			    double   epsilon,
			    double   *seed) const
//---------------------------------------------------------------------------
{
  // Compute closest point with respect to the underlying volume
  vol_->closestPoint(pt, clo_u, clo_v, clo_w, clo_pt, clo_dist, 
		     epsilon, seed);

  // Check if the found point is within the trimmed volume
  bool inside = isInside(clo_pt);
  if (inside)
    return;  // The closest point is found

  double clo_par_u, clo_par_v;
  (void)closestBoundaryPoint(pt, clo_u, clo_v, clo_w, clo_pt, 
			     clo_dist, clo_par_u, clo_par_v,
			     epsilon);
  
}

//---------------------------------------------------------------------------
ftFaceBase* ftVolume::closestBoundaryPoint(Point& pt,
					   double&  clo_u,
					   double&  clo_v, 
					   double&  clo_w, 
					   Point& clo_pt,
					   double&  clo_dist,
					   double& clo_par_u,
					   double& clo_par_v,
					   double epsilon) const
//---------------------------------------------------------------------------
{
  vector<shared_ptr<SurfaceModel> > shells = getAllShells();
  
  // Compute closest point in the first bounding shell
  int idx;
  double clo_par[2];
  shells[0]->closestPoint(pt, clo_pt, idx, clo_par, clo_dist);

  // Check if the closest points in the inner shells have a smaller distance
  // to the input point
  int shell_ix = 0;
  for (size_t ki=1; ki<shells.size(); ++ki)
    {
      int idx2;
      double clo_par2[2];
      Point clo_pt2;
      double clo_dist2;
      shells[ki]->closestPoint(pt, clo_pt2, idx2, clo_par2, clo_dist2);
      if (clo_dist2 < clo_dist)
	{
	  shell_ix = (int)ki;
	  idx = idx2;
	  clo_pt = clo_pt2;
	  clo_par[0] = clo_par2[0];
	  clo_par[1] = clo_par2[1];
	  clo_dist = clo_dist2;
	}
    }

  clo_par_u = clo_par[0];
  clo_par_v = clo_par[1];

  // Compute volume parameter value
  double dist;
  Point pt2;
  vol_->closestPoint(clo_pt, clo_u, clo_v, clo_w, pt2, dist, epsilon);
  return shells[shell_ix]->getFace(idx).get();
}

//---------------------------------------------------------------------------
bool ftVolume::isSpline() const
//---------------------------------------------------------------------------
{
  return vol_->isSpline();
}

//===========================================================================
vector<shared_ptr<EdgeVertex> > ftVolume::radialEdges() const
//===========================================================================
{
  vector<shared_ptr<EdgeVertex> > result;

  std::set<shared_ptr<EdgeVertex> > tmp_result;  // To get the radial edges
  // only once
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> face = shells_[ki]->getFace(kj);
	  vector<shared_ptr<EdgeVertex> > edges = face->getRadialEdges();

	  tmp_result.insert(edges.begin(), edges.end());
	}
    }
  result.insert(result.begin(), tmp_result.begin(), tmp_result.end());
  return result;
}

//===========================================================================
vector<shared_ptr<ftEdge> > ftVolume::uniqueNonRadialEdges() const
//===========================================================================
{
  vector<shared_ptr<ftEdge> > result;

  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> face = shells_[ki]->getFace(kj);
	  vector<shared_ptr<ftEdge> > edges = face->getAllEdges();

	  // Store edges that have no radial edge and when the
	  // twin is not stored already
	  for (size_t kr=0; kr<edges.size(); ++kr)
	    {
	      if (!edges[kr]->hasEdgeMultiplicity())
		{
		  size_t kh;
		  for (kh=0; kh<result.size(); ++kh)
		    if (edges[kr].get() == result[kh].get() ||
			edges[kr]->twin() == result[kh].get())
		      break;
		  if (kh == result.size())
		    result.push_back(edges[kr]);
		}
	    }
	}
    }
  return result;
}

//===========================================================================
vector<shared_ptr<EdgeVertex> > 
ftVolume::getCommonEdges(ftVolume *other) const
//===========================================================================
{
  vector<shared_ptr<EdgeVertex> > evx1 = radialEdges();
  vector<shared_ptr<EdgeVertex> > evx2 = other->radialEdges();
  vector<shared_ptr<EdgeVertex> > evx3;
  for (size_t ki=0; ki<evx1.size(); ++ki)
    for (size_t kj=0; kj<evx2.size(); ++kj)
      {
	if (evx1[ki].get() == evx2[kj].get())
	  {
	    evx3.push_back(evx1[ki]);
	    break;
	  }
      }
  return evx3;
}

//===========================================================================
// 
// Get neighbouring bodies
    void ftVolume::getAdjacentBodies(std::vector<ftVolume*>& neighbours)
//===========================================================================
{
    std::set<ftVolume*> all_bodies;

    for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb_faces = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb_faces; ++kj)
	{
	  shared_ptr<ftSurface> curr = shells_[ki]->getFace(kj);

	  ftSurface* twin = curr->twin();

	  if (!twin)
	    continue;  // No twin information

	  Body *twin_body = twin->getBody();

	  // All neighbours to a ftVolume are ftVolumes. Thus, this operation should
	  // be safe
	  ftVolume *curr_volume = dynamic_cast<ftVolume*>(twin_body);
	  all_bodies.insert(curr_volume);
	}
    }
    neighbours.insert(neighbours.end(), all_bodies.begin(), all_bodies.end());
}

//===========================================================================
// 
// Check for adjacency, splines and common spline spaces
bool ftVolume::commonSplineSpace(ftVolume *other, double tol)
//===========================================================================
{
  if (!isSpline() || !other->isSpline())
    return false;

  shared_ptr<ftSurface> bd_face1, bd_face2;
  bool neighbours = areNeighbours(other, bd_face1, bd_face2);
  if (!neighbours)
    return true;  // Not applicable

  // Fetch surface geometry
  shared_ptr<ParamVolume> vol1 = getVolume();
  shared_ptr<ParamVolume> vol2 = other->getVolume();
  shared_ptr<SplineVolume> splvol1 = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol1);
  shared_ptr<SplineVolume> splvol2 = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol2);
  if (!splvol1.get() || !splvol2.get())
    {
      MESSAGE("Check data structure. Inconsistency of surface types");
      return false;
    }

  // Check that either none or both volumes are rational
  if ((splvol1->rational() && !splvol2->rational()) ||
      (!splvol1->rational() && splvol2->rational()))
      return false;

  VolumeAdjacencyInfo adj_info = getAdjacencyInfo(other, DEFAULT_SPACE_EPSILON);
  if (!adj_info.adjacency_found_)
    return false;

  // Index of common boundary curve, 0=umin, 1=umax, 2=vmin, 3=vmax
  int bd1 = adj_info.bd_idx_1_;
  int bd2 = adj_info.bd_idx_2_;

  // Fetch spline spaces
  BsplineBasis basis1_1 = (bd1 == 0 || bd1 == 1) ? splvol1->basis(1) : 
    splvol1->basis(0);
  BsplineBasis basis1_2 = (bd1 == 4 || bd1 == 5) ? splvol1->basis(1) : 
    splvol1->basis(2);
  BsplineBasis basis2_1 = (bd2 == 0 || bd2 == 1) ? splvol2->basis(1) : 
    splvol2->basis(0);
  BsplineBasis basis2_2 = (bd2 == 4 || bd2 == 5) ? splvol2->basis(1) : 
    splvol2->basis(2);

  // Make copy of basis 2
  BsplineBasis basis2_3 = basis2_1;
  BsplineBasis basis2_4 = basis2_2;
  
  // Transform the basises of volume 2 to match those of 1 with regard to 
  // sequence, orientation and domain
  if (!adj_info.same_dir_order_)
    std::swap(basis2_3, basis2_4);

  if (!adj_info.same_orient_u_)
    basis2_3.reverseParameterDirection();
  if (!adj_info.same_orient_v_)
    basis2_4.reverseParameterDirection();

  basis2_3.rescale(basis1_1.startparam(), basis1_1.endparam());
  basis2_4.rescale(basis1_2.startparam(), basis1_2.endparam());

  // Check equality

  bool same = basis1_1.sameSplineSpace(basis2_3);
  if (!same)
    return false;
  same = basis1_2.sameSplineSpace(basis2_4);
  return same;
}

//===========================================================================
// 
// Given two adjacent spline volumes, represent them in a common
// spline space
bool ftVolume::makeCommonSplineSpace(ftVolume *other)
//===========================================================================
{
  if (!isSpline() || !other->isSpline())
    return false;

  shared_ptr<ftSurface> bd_face1, bd_face2;
  bool neighbours = areNeighbours(other, bd_face1, bd_face2);
  if (!neighbours)
    return true;  // Not applicable

  VolumeAdjacencyInfo adj_info = getAdjacencyInfo(other, DEFAULT_SPACE_EPSILON);
  if (!adj_info.adjacency_found_)
    return false;

  // Index of common boundary curve, 0=umin, 1=umax, 2=vmin, 3=vmax
  int bd1 = adj_info.bd_idx_1_;
  int bd2 = adj_info.bd_idx_2_;
  int orientation = 0;
  if (!adj_info.same_orient_u_)
    ++orientation;
  if (!adj_info.same_orient_v_)
    orientation += 2;

  // Fetch surface geometry
  shared_ptr<ParamVolume> vol1 = getVolume();
  shared_ptr<ParamVolume> vol2 = other->getVolume();
  shared_ptr<SplineVolume> splvol1 = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol1);
  shared_ptr<SplineVolume> splvol2 = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol2);
  if (!splvol1.get() || !splvol2.get())
    {
      MESSAGE("Check data structure. Inconsistency of surface types");
      return false;
    }

  VolumeTools::volCommonSplineSpace(splvol1, bd1, splvol2, bd2, orientation, 
		       adj_info.same_dir_order_);
  return true;
}

//===========================================================================
// 
// Check for adjacency, splines and corner to corner configuration
bool ftVolume::isCornerToCorner(shared_ptr<ftVolume> other,
				double tol)
//===========================================================================
{
  if (!isSpline() || !other->isSpline())
    return false;

  shared_ptr<ftSurface> bd_face1, bd_face2;
  bool neighbours = areNeighbours(other.get(), bd_face1, bd_face2);
  if (!neighbours)
    return true;  // Not applicable

  // Must be implemented

  // Fetch geometry, check edge type
  shared_ptr<ParamVolume> vol1 = getVolume();
  shared_ptr<ParamVolume> vol2 = other->getVolume();
  shared_ptr<ParamSurface> bdsf1 = bd_face1->surface();
  shared_ptr<ParamSurface> bdsf2 = bd_face2->surface();
  shared_ptr<SurfaceOnVolume> volsf1 = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf1);
  shared_ptr<SurfaceOnVolume> volsf2 = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf2);
  
  if (!volsf1.get() || !volsf2.get())
    {
      MESSAGE("Check data structure. Expecting curve on surface");
      return false;
    }

  bool at_corners = VolumeTools::cornerToCornerVols(vol1, volsf1, vol2, volsf2, tol);
  return at_corners;

  return true;  // For the time being
}

//===========================================================================
void ftVolume::splitAtInternalCorner(ftVolume* other,
				     vector<shared_ptr<ftVolume> >& new_vol1,
				     vector<shared_ptr<ftVolume> >& new_vol2,
				     double tol)
//===========================================================================
{

  // Must be implemented
  MESSAGE("ftVolume::splitAtInternalCorner. Not implemented");
}


//===========================================================================
// 
// 
bool ftVolume::getAdjacencyInfo(ftVolume *other, double tol,
				int& bd1, int& bd2, 
				int& orientation, bool& same_seq)
//===========================================================================
{
  MESSAGE("Using deprecated getAdjacencyInfo() method. Use other instead");

  VolumeAdjacencyInfo adj_info = getAdjacencyInfo(other, tol);
  if (adj_info.adjacency_found_)
    {
      bd1 = adj_info.bd_idx_1_;
      bd2 = adj_info.bd_idx_2_;
      orientation = 0;
      if (!adj_info.same_orient_u_)
	++orientation;
      if (!adj_info.same_orient_v_)
	orientation += 2;
      same_seq = adj_info.same_dir_order_;
    }
  return adj_info.adjacency_found_;
}


//===========================================================================
// 
// 
VolumeAdjacencyInfo ftVolume::getAdjacencyInfo(ftVolume *other, double tol,
					    int adj_idx, bool test_corner)
//===========================================================================
{
  VolumeAdjacencyInfo adj_info;

  shared_ptr<ftSurface> bd_face1, bd_face2;
  bool neighbours = areNeighbours(other, bd_face1, bd_face2, adj_idx);

  shared_ptr<ParamVolume> vol1 = getVolume();
  shared_ptr<ParamVolume> vol2 = other->getVolume();
  shared_ptr<ParamSurface> bdsf1, bdsf2;
  if (!neighbours)
    {
      // Look for adjacency in a degenerate situation
      vector<shared_ptr<EdgeVertex> > evx = getCommonEdges(other);

      if (evx.size() == 0)
	{
	  adj_info.adjacency_found_ = false;  // Not applicable
	  return adj_info;
	}

      // Check if it is a degenerate situation
      size_t kj;
      for (kj=0; kj<evx.size(); ++kj)
	{
	  bool deg_found = checkDegAdjacency(other, evx[kj], tol,
					     bdsf1, bdsf2);
	  if (deg_found)
	    break;
	}
      if (kj == evx.size())
	{
	  // No degenerate case found
	  adj_info.adjacency_found_ = false;  // Not applicable
	  return adj_info;
	}	  

    }
  else
    {
      // Fetch geometry, check edge type
      bdsf1 = bd_face1->surface();
      bdsf2 = bd_face2->surface();
    }

  // Neighbourhood established. Fetch adjacency details
  shared_ptr<SurfaceOnVolume> volsf1 = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf1);
  shared_ptr<SurfaceOnVolume> volsf2 = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf2);
  
  if (!volsf1.get() || !volsf2.get())
    {
      MESSAGE("Check data structure. Expecting surface on volume");
      adj_info.adjacency_found_ = false;
      return adj_info;
    }

  int bd1=-1, bd2=-1, orientation=-1;
  bool same_seq=true;
  adj_info.adjacency_found_ = 
      VolumeTools::getVolAdjacencyInfo(vol1, volsf1, vol2, volsf2,
				       tol, bd1, bd2, orientation, same_seq);
  adj_info.bd_idx_1_ = bd1;
  adj_info.bd_idx_2_ = bd2;
  adj_info.same_orient_u_ = (orientation != 1 && orientation != 3);
  adj_info.same_orient_v_ = (orientation != 2 && orientation != 3);
  adj_info.same_dir_order_ = same_seq;

  if (test_corner && adj_info.adjacency_found_)
    // Should use something like
    //
    //    adj_info.corner_failed_ = !cornerToCornerVol(vol1, volsf1, vol2, volsf2, tol);
    //
    // like for ftSurface, but this method does not exist yet.
    adj_info.corner_failed_ = false;

  return adj_info;
}

//===========================================================================
// 
// 
VolumeAdjacencyInfo ftVolume::getCornerAdjacencyInfo(ftVolume *other, 
						     double tol, int adj_idx )
//===========================================================================
{
  VolumeAdjacencyInfo adj_info;

  shared_ptr<ftSurface> bd_face1, bd_face2;
  bool neighbours = areNeighbours(other, bd_face1, bd_face2);
  if (neighbours)
    return adj_info; // Adjacency with common boundary faces

  // Look for adjacency along an edge
  vector<shared_ptr<EdgeVertex> > evx = getCommonEdges(other);

  // Check if a new edge adjacency is found
  if ((int)evx.size() <= adj_idx)
    return adj_info;  // No new edge adjacecy is found

  // Fetch edges from the radial edge related to the current volumes
  vector<ftEdge*> edges = evx[adj_idx]->allEdges();

  // Find edges belonging to the given volumes
  ftEdge *e1=0, *e2=0;
  bool found1 = false, found2 = false;
  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      ftSurface *curr = edges[ki]->face()->asFtSurface();
      if ((!found1 && curr->getBody() == this) || 
	  (!found2 && curr->getBody() == other))
	{
	  if (!found1 && curr->getBody() == this)
	    found1 = true;
	  else if (!found2 && curr->getBody() == other)
	    found2 = true;

	  if (!e1)
	    e1 = edges[ki];
	  else
	    {
	      e2 = edges[ki];
	      break;
	    }
	}
    }

  if (!(e1 && e2))
    return adj_info;   // No adjacency configuration found
  
  // Fetch adjacent faces
  ftSurface *f1 = e1->face()->asFtSurface();
  ftSurface *f2 = e2->face()->asFtSurface();

  if (!(f1 && f2))
    return adj_info;   // No adjacency configuration found
  
  // Check sequence
  if (f1->getBody() != this)
    {
      std::swap(e1, e2);
      std::swap(f1, f2);
    }

  // Fetch boundary index of faces
  shared_ptr<ParamSurface> bdsf1 = f1->surface();
  shared_ptr<ParamSurface> bdsf2 = f2->surface();
  shared_ptr<SurfaceOnVolume> vol_sf1 = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf1);
  shared_ptr<SurfaceOnVolume> vol_sf2 = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf2);
  
  if (!vol_sf1.get() || !vol_sf2.get())
    return adj_info; // Missing information

  int orientation1, orientation2;
  bool swap1, swap2;
  int bd1 = vol_sf1->whichBoundary(tol, orientation1, swap1);
  int bd2 = vol_sf2->whichBoundary(tol, orientation2, swap2);
  if (bd1 < 0 || bd2 < 0)
    return adj_info;  // Not boundary trim

  // Fetch boundary index of edge with respect to face
  shared_ptr<ParamCurve> cv1 = e1->geomCurve();
  shared_ptr<ParamCurve> cv2 = e2->geomCurve();
  shared_ptr<CurveOnSurface> sf_cv1 = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
  shared_ptr<CurveOnSurface> sf_cv2 = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
  if (!sf_cv1.get() || !sf_cv2.get())
    return adj_info;  // Unexpected curve type

  bool same_orient1, same_orient2;
  int bd_cv1 = sf_cv1->whichBoundary(tol, same_orient1);
  int bd_cv2 = sf_cv2->whichBoundary(tol, same_orient2);
  if (bd_cv1 < 0 || bd_cv2 < 0)
    return adj_info;  // Not boundary trim

  // Check if the endpoints of the boundary curves correspond
  // The radial edge may be longer than associated edges
  Point pc1 = sf_cv1->ParamCurve::point(sf_cv1->startparam());
  Point pc2 = sf_cv1->ParamCurve::point(sf_cv1->endparam());
  Point pc3 = sf_cv2->ParamCurve::point(sf_cv2->startparam());
  Point pc4 = sf_cv2->ParamCurve::point(sf_cv2->endparam());
  if (!((pc1.dist(pc3) < tol && pc2.dist(pc4) < tol) ||
	(pc1.dist(pc4) < tol && pc2.dist(pc3) < tol)))
    return adj_info; // Not corresponding endpoints

  // Check orientation
  Point p1 = e1->point(e1->tMin());
  Point p2 = e1->point(e1->tMax());
  Point p3 = e2->point(e2->tMin());
  Point p4 = e2->point(e2->tMax());
  bool same = (p1.dist(p3) + p2.dist(p4) < p1.dist(p4) + p2.dist(p3));

  if (e1->isReversed())
    same_orient1 = !same_orient1;
  if (e2->isReversed())
    same_orient2 = !same_orient2;
  if (same_orient1 != same_orient2)
    same = !same;

  // Set output info
  adj_info.corner_adjacency_ = true;
  adj_info.bd_idx_1_ = bd1;
  adj_info.bd_idx_2_ = bd2;
  adj_info.edg_idx_1_  = bd_cv1;
  adj_info.edg_idx_2_  = bd_cv2;
  adj_info.same_orient_edge_ = same;
  
  return adj_info;
}

//===========================================================================
// 
// 
VolumeAdjacencyInfo ftVolume::getCornerAdjacencyInfo(ftVolume *other, 
						     EdgeVertex* evx,
						     double tol, int adj_idx )
//===========================================================================
{
  VolumeAdjacencyInfo adj_info;

  // Look for adjacency along an edge
  // Fetch edges from the radial edge related to the current volumes
  vector<ftEdge*> edges = evx->allEdges();

  // Find edges belonging to the given volumes
  ftEdge *e1=0, *e2=0;
  bool found1 = false, found2 = false;
  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      ftSurface *curr = edges[ki]->face()->asFtSurface();
      if ((!found1 && curr->getBody() == this) || 
	  (!found2 && curr->getBody() == other))
	{
	  if (!found1 && curr->getBody() == this)
	    found1 = true;
	  else if (!found2 && curr->getBody() == other)
	    found2 = true;

	  if (!e1)
	    e1 = edges[ki];
	  else
	    {
	      e2 = edges[ki];
	      break;
	    }
	}
    }

  if (!(e1 && e2))
    return adj_info;   // No additional adjacency configuration found
  
  // Fetch adjacent faces
  ftSurface *f1 = e1->face()->asFtSurface();
  ftSurface *f2 = e2->face()->asFtSurface();

  if (!(f1 && f2))
    return adj_info;   // No adjacency configuration found
  
  // Check sequence
  if (f1->getBody() != this)
    {
      std::swap(e1, e2);
      std::swap(f1, f2);
    }

  // Fetch boundary index of faces
  shared_ptr<ParamSurface> bdsf1 = f1->surface();
  shared_ptr<ParamSurface> bdsf2 = f2->surface();
  shared_ptr<SurfaceOnVolume> vol_sf1 = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf1);
  shared_ptr<SurfaceOnVolume> vol_sf2 = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf2);
  
  if (!vol_sf1.get() || !vol_sf2.get())
    return adj_info; // Missing information

  int orientation1, orientation2;
  bool swap1, swap2;
  int bd1 = vol_sf1->whichBoundary(tol, orientation1, swap1);
  int bd2 = vol_sf2->whichBoundary(tol, orientation2, swap2);
  if (bd1 < 0 || bd2 < 0)
    return adj_info;  // Not boundary trim

  // Fetch boundary index of edge with respect to face
  shared_ptr<ParamCurve> cv1 = e1->geomCurve();
  shared_ptr<ParamCurve> cv2 = e2->geomCurve();
  shared_ptr<CurveOnSurface> sf_cv1 = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
  shared_ptr<CurveOnSurface> sf_cv2 = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
  if (!sf_cv1.get() || !sf_cv2.get())
    return adj_info;  // Unexpected curve type

  bool same_orient1, same_orient2;
  int bd_cv1 = sf_cv1->whichBoundary(tol, same_orient1);
  int bd_cv2 = sf_cv2->whichBoundary(tol, same_orient2);
  if (bd_cv1 < 0 || bd_cv2 < 0)
    return adj_info;  // Not boundary trim

  // Check if the endpoints of the boundary curves correspond
  // The radial edge may be longer than associated edges
  Point pc1 = sf_cv1->ParamCurve::point(sf_cv1->startparam());
  Point pc2 = sf_cv1->ParamCurve::point(sf_cv1->endparam());
  Point pc3 = sf_cv2->ParamCurve::point(sf_cv2->startparam());
  Point pc4 = sf_cv2->ParamCurve::point(sf_cv2->endparam());
  if (!((pc1.dist(pc3) < tol && pc2.dist(pc4) < tol) ||
	(pc1.dist(pc4) < tol && pc2.dist(pc3) < tol)))
    return adj_info; // Not corresponding endpoints

  // Check also midpoint in case the two curves make up a loop
  Point pc5 = sf_cv1->ParamCurve::point(0.5*(sf_cv1->startparam()+sf_cv1->endparam()));
  Point pc6;
  double t6, d6;
  sf_cv2->closestPoint(pc5, sf_cv2->startparam(), sf_cv2->endparam(), t6, pc6, d6);
  if (d6 > tol)
    return adj_info;

  // Relate index to volume, not boundary surface
  if (swap1)
    {
      int n = (bd_cv1%2 == 0) ? 2 : 4;
      bd_cv1 = n - bd_cv1;
    }
  if (swap2)
    {
      int n = (bd_cv2%2 == 0) ? 2 : 4;
      bd_cv2 = n - bd_cv2;
    }

  // Check orientation
  Point p1 = e1->point(e1->tMin());
  Point p2 = e1->point(e1->tMax());
  Point p3 = e2->point(e2->tMin());
  Point p4 = e2->point(e2->tMax());
  bool same = (p1.dist(p3) + p2.dist(p4) < p1.dist(p4) + p2.dist(p3));

  if (e1->isReversed())
    same_orient1 = !same_orient1;
  if (e2->isReversed())
    same_orient2 = !same_orient2;
  if (same_orient1 != same_orient2)
    same = !same;

  // Set output info
  adj_info.corner_adjacency_ = true;
  adj_info.bd_idx_1_ = bd1;
  adj_info.bd_idx_2_ = bd2;
  adj_info.edg_idx_1_  = bd_cv1;
  adj_info.edg_idx_2_  = bd_cv2;
  adj_info.same_orient_edge_ = same;
  
  return adj_info;
}


//===========================================================================
// 
// 
bool ftVolume::checkDegAdjacency(ftVolume *other, 
				 shared_ptr<EdgeVertex> evx,
				  double tol,
				  shared_ptr<ParamSurface>& bdsf1, 
				  shared_ptr<ParamSurface>& bdsf2)
//===========================================================================
{
  // Check for spline volume
  shared_ptr<ParamVolume> vol1 = getVolume();
  shared_ptr<ParamVolume> vol2 = other->getVolume();
  shared_ptr<SplineVolume> svol1 = 
    dynamic_pointer_cast<SplineVolume,ParamVolume>(vol1);
  shared_ptr<SplineVolume> svol2 = 
    dynamic_pointer_cast<SplineVolume,ParamVolume>(vol2);

  if (!(svol1.get() && svol2.get()))
    return false;  // VSK, May 2010. This is a feasible solution now, 
  // but not when more volume types are implemented

  // Check if both volumes are degenerate in such a way that a boundary
  // surface degenerate to a curve or a point
  int is_deg1[6], is_deg2[6];
  svol1->checkDegeneracy(tol, is_deg1);
  svol2->checkDegeneracy(tol, is_deg2);

  // Fetch an edge corresponding to the edge vertex. This is a simple
  // solution that might have some tolerance problems. In practice, however,
  // we expect a watertight model. Thus, the accuracy should be good enough
  if (evx->nmbUniqueEdges() == 0)
    return false;  // No edge vertex info
  ftEdge *edge = evx->getEdge(0);

  // First volume
  int nmb_sample = 5;
  bool found_edge1 = false, found_edge2 = false;
  for (int ki=0; ki<6; ++ki)
    {
      if (is_deg1[ki] >= 2)
	{
	  // Degenerate boundary surface, check if it corresponds to the
	  // common edge vertex
	  int type;
	  bool b, r, t, l;
	  bool found_deg = svol1->isDegenerate(ki, type, b, r, t, l, tol);
	  if (!found_deg || type<2)
	    continue;  // Actually already checked and should not be the case

	  shared_ptr<SplineSurface> bd_sf = svol1->getBoundarySurface(ki);
	  if (ki == 1 || ki == 2 || ki == 5)
	    bd_sf->swapParameterDirection();
	  if (b && t)
	    {
	      // Check if the mid curve in second parameter direction
	      // follows the common vertex edge
	      double t1 = bd_sf->startparam_v();
	      double t2 = bd_sf->endparam_v();
	      double t3 = 0.5*(bd_sf->startparam_u()+ bd_sf->endparam_u());
	      double tdel = (t2 - t1)/(double)(nmb_sample-1);
	      double tpar;
	      int kj;
	      for (kj=0, tpar=t1; kj<nmb_sample; ++kj, tpar+=tdel)
		{
		  Point pos = bd_sf->ParamSurface::point(t3, tpar);

		  // Relax onto the edge vertex
		  double par, dist;
		  Point edge_pos;
		  edge->closestPoint(pos, par, edge_pos, dist);
		  if (dist > tol)
		    break;  // No adjacency
		}

	      if (kj == nmb_sample)
		{
		  found_edge1 = true;

		  // Construct SurfaceOnVolume corresponding to 
		  // current boundary surface
		  bdsf1 = shared_ptr<ParamSurface>(VolumeTools::getOrientedBoundarySurface(svol1,
											   ki));
		}
	    }

	  if (r && l && !found_edge1)
	    {
	      // Check if the mid curve in first parameter direction
	      // follows the common vertex edge
	      double t1 = bd_sf->startparam_u();
	      double t2 = bd_sf->endparam_u();
	      double t3 = 0.5*(bd_sf->startparam_v()+ bd_sf->endparam_v());
	      double tdel = (t2 - t1)/(double)(nmb_sample-1);
	      double tpar;
	      int kj;
	      for (kj=0, tpar=t1; kj<nmb_sample; ++kj, tpar+=tdel)
		{
		  Point pos = bd_sf->ParamSurface::point(t3, tpar);

		  // Relax onto the edge vertex
		  double par, dist;
		  Point edge_pos;
		  edge->closestPoint(pos, par, edge_pos, dist);
		  if (dist > tol)
		    break;  // No adjacency
		}

	      if (kj == nmb_sample)
		{
		  found_edge1 = true;

		  // Construct SurfaceOnVolume corresponding to 
		  // current boundary surface
		  bdsf1 = shared_ptr<ParamSurface>(VolumeTools::getOrientedBoundarySurface(svol1,
								      ki));
		}
	    }
	}
      if (found_edge1)
	break;
    }

  // Second volume
  for (int ki=0; ki<6; ++ki)
    {
      if (is_deg2[ki] >= 2)
	{
	  // Degenerate boundary surface, check if it corresponds to the
	  // common edge vertex
	  int type;
	  bool b, r, t, l;
	  bool found_deg = svol2->isDegenerate(ki, type, b, r, t, l, tol);
	  if (!found_deg || type<2)
	    continue;  // Actually already checked and should not be the case

	  shared_ptr<SplineSurface> bd_sf = svol2->getBoundarySurface(ki);
	  if (ki == 1 || ki == 2 || ki == 5)
	    bd_sf->swapParameterDirection();
	  if (b && t)
	    {
	      // Check if the mid curve in second parameter direction
	      // follows the common vertex edge
	      double t1 = bd_sf->startparam_v();
	      double t2 = bd_sf->endparam_v();
	      double t3 = 0.5*(bd_sf->startparam_u()+ bd_sf->endparam_u());
	      double tdel = (t2 - t1)/(double)(nmb_sample-1);
	      double tpar;
	      int kj;
	      for (kj=0, tpar=t1; kj<nmb_sample; ++kj, tpar+=tdel)
		{
		  Point pos = bd_sf->ParamSurface::point(t3, tpar);

		  // Relax onto the edge vertex
		  double par, dist;
		  Point edge_pos;
		  edge->closestPoint(pos, par, edge_pos, dist);
		  if (dist > tol)
		    break;  // No adjacency
		}

	      if (kj == nmb_sample)
		{
		  found_edge2 = true;

		  // Construct SurfaceOnVolume corresponding to 
		  // current boundary surface
		  bdsf2 = shared_ptr<ParamSurface>(VolumeTools::getOrientedBoundarySurface(svol2,
								      ki));
		}
	    }

	  if (r && l && !found_edge2)
	    {
	      // Check if the mid curve in first parameter direction
	      // follows the common vertex edge
	      double t1 = bd_sf->startparam_u();
	      double t2 = bd_sf->endparam_u();
	      double t3 = 0.5*(bd_sf->startparam_v()+ bd_sf->endparam_v());
	      double tdel = (t2 - t1)/(double)(nmb_sample-1);
	      double tpar;
	      int kj;
	      for (kj=0, tpar=t1; kj<nmb_sample; ++kj, tpar+=tdel)
		{
		  Point pos = bd_sf->ParamSurface::point(t3, tpar);

		  // Relax onto the edge vertex
		  double par, dist;
		  Point edge_pos;
		  edge->closestPoint(pos, par, edge_pos, dist);
		  if (dist > tol)
		    break;  // No adjacency
		}

	      if (kj == nmb_sample)
		{
		  found_edge2 = true;

		  // Construct SurfaceOnVolume corresponding to 
		  // current boundary surface
		  bdsf2 = shared_ptr<ParamSurface>(VolumeTools::getOrientedBoundarySurface(svol2,
								      ki));
		}
	    }
	}
      if (found_edge2)
	break;
    }

  return (found_edge1 && found_edge2);
}

 //===========================================================================
bool ftVolume::getCorrCoefEnumeration(ftVolume *other, double tol,
				      vector<pair<int,int> >& enumeration)
//===========================================================================
{
  // Fetch adjacency information
  VolumeAdjacencyInfo adj_info = getAdjacencyInfo(other, tol);
  if (!adj_info.adjacency_found_)
    return false;

  // Index of common boundary curve, 0=umin, 1=umax, 2=vmin, 3=vmax
  int bd1 = adj_info.bd_idx_1_;
  int bd2 = adj_info.bd_idx_2_;
  int orientation = 0;
  if (!adj_info.same_orient_u_)
    ++orientation;
  if (!adj_info.same_orient_v_)
    orientation += 2;

  // Check that the surfaces have corresponding spline spaces
  bool common = commonSplineSpace(other, tol);
  if (!common)
    return false; // No corresponding coefficient enumeration

  // Fetch surface geometry
  shared_ptr<ParamVolume> vol1 = getVolume();
  shared_ptr<ParamVolume> vol2 = other->getVolume();
  shared_ptr<SplineVolume> splvol1 = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol1);
  shared_ptr<SplineVolume> splvol2 = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol2);
  if (!splvol1.get() || !splvol2.get())
    return false;

  bool pairwise = VolumeTools::getCorrCoefVolEnum(splvol1, splvol2,
				     bd1, bd2, orientation,
				     adj_info.same_dir_order_, 
				     enumeration);
  return pairwise;
}


//===========================================================================
bool ftVolume::getVertexPosition(shared_ptr<Vertex> vx, Point& param) const
//===========================================================================
{
  // Get all faces associated this vertex
  vector<pair<ftSurface*, Point> > faces = vx->getFaces();
  size_t ki;
  for (ki=0; ki<faces.size(); ++ki)
    if (faces[ki].first->getBody() == this)
      break;

  if (ki == faces.size())
    return false;  // Not associated

  // Get surface
  shared_ptr<ParamSurface> srf = faces[ki].first->surface();
  shared_ptr<SurfaceOnVolume> bd_srf = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(srf);
  if (!bd_srf.get())
    {
      MESSAGE("Expecting SurfaceOnVolume. Check data structure");
      return false;
    }

  // Get volume parameter
  Point sf_par = faces[ki].second;
  param = bd_srf->volumeParameter(sf_par[0], sf_par[1]);
  return true;
}

//===========================================================================
bool ftVolume::getVertexEnumeration(shared_ptr<Vertex> vx, 
				    Point& param, int& corner,
				    int& coef_nmb) const
//===========================================================================
{
  // Get spline volume
  shared_ptr<SplineVolume> vol = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol_);
  if (!vol.get())
    return false; // Not a spline volume

  // Get volume parameter
  bool found = getVertexPosition(vx, param);
  if (!found)
    return false;

  // Define corner parameters
  Point cpar[8];
  cpar[0] = Point(vol->startparam(0),vol->startparam(1),vol->startparam(2));
  cpar[1] = Point(vol->endparam(0),vol->startparam(1),vol->startparam(2));
  cpar[2] = Point(vol->startparam(0),vol->endparam(1),vol->startparam(2));
  cpar[3] = Point(vol->endparam(0),vol->endparam(1),vol->startparam(2));
  cpar[4] = Point(vol->startparam(0),vol->startparam(1),vol->endparam(2));
  cpar[5] = Point(vol->endparam(0),vol->startparam(1),vol->endparam(2));
  cpar[6] = Point(vol->startparam(0),vol->endparam(1),vol->endparam(2));
  cpar[7] = Point(vol->endparam(0),vol->endparam(1),vol->endparam(2));

  // Compute distances
  double dist, mindist=1.0e8;
  int ki;
  corner = -1;
  for (ki=0; ki<8; ki++)
    {
      dist = param.dist(cpar[ki]);
      if (dist < mindist)
	{
	  mindist = dist;
	  corner = ki;
	}
    }
  if (corner < 0)
    return false;   // Something went wrong
  
  int kn1 = vol->numCoefs(0);
  int kn2 = vol->numCoefs(1);
  int kn3 = vol->numCoefs(2);

  if (corner == 0)
    coef_nmb = 0;
  else if (corner == 1)
    coef_nmb = kn1-1;
  else if (corner == 2)
    coef_nmb = kn1*(kn2-1);
  else if (corner == 3)
    coef_nmb = kn1*kn2-1;
  else if (corner == 4)
    coef_nmb = kn1*kn2*(kn3-1);
  else if (corner == 5)
    coef_nmb = kn1*kn2*(kn3-1) + kn1-1;
  else if (corner == 6)
    coef_nmb = kn1*kn2*(kn3-1) + kn1*(kn2-1);
  else
    coef_nmb = kn1*kn2*kn3-1;

  return true;
}

//===========================================================================
// 
// 
bool ftVolume::getFreeBoundaryInfo(double tol, 
				   vector<int>& free_boundaries) 
//===========================================================================
{
  bool all_at_boundaries = true;
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb_faces = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb_faces; ++kj)
	{
	  shared_ptr<ftSurface> curr = shells_[ki]->getFace(kj);

	  if (curr->twin())
	    continue;  // Not an outer boundary

	  // Get surface
	  shared_ptr<ParamSurface> srf = curr->surface();
	  shared_ptr<SurfaceOnVolume> bd_srf = 
	    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(srf);
	  if (bd_srf.get())
	    {
	      int orientation;
	      bool swap;
	      int bd = bd_srf->whichBoundary(tol, orientation, swap);
	      if (bd < 0)
		all_at_boundaries = false;
	      else
		free_boundaries.push_back(bd);
	    }
	  else
	    all_at_boundaries = false;
	}
    }
  return all_at_boundaries;
}

//===========================================================================
// 
// 
bool ftVolume::getBoundaryCoefEnumeration(int bd, 
					  std::vector<int>& enumeration) 
//===========================================================================
{
  // Get spline volume
  shared_ptr<SplineVolume> vol = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol_);
  if (!vol.get())
    return false; // Not a spline volume

  if (bd < 0 || bd > 5)
    return false;  // No boundary is specified

  bool found = VolumeTools::getVolCoefEnumeration(vol, bd, enumeration);
  return found;
}

//===========================================================================
// 
// 
void ftVolume::removeSliverFaces(double len_tol)
//===========================================================================
{
  // Search for sliver faces
  shared_ptr<SurfaceModel> shell = getOuterShell();
  if (!shell.get())
    return;  // Nothing to do
  if (shell->nmbBoundaries() > 0)
    return;  // Expecting a closed shell
  if (shell->nmbEntities() <= 6)
    return;

#ifdef DEBUG_VOL1
  std::ofstream of0("pre_sliver.g2");
  vector<shared_ptr<Vertex> > vx;

  int nmb_face = shell->nmbEntities();
  for (int ka=0; ka<nmb_face; ++ka)
    {
      shared_ptr<ParamSurface> tmp = shell->getSurface(ka);
      tmp->writeStandardHeader(of0);
      tmp->write(of0);
    }
  vector<shared_ptr<Vertex> > tmp_vx;
  shell->getAllVertices(tmp_vx);
  vx.insert(vx.end(), tmp_vx.begin(), tmp_vx.end());
  of0 << "400 1 0 4 255 0 0 255" << std::endl;
  of0 << vx.size() << std::endl;
  for (size_t ki=0; ki<vx.size(); ++ki)
    of0 << vx[ki]->getVertexPoint() << std::endl;
#endif

  tpTolerances toptol = shell->getTolerances();

  double tol1 = 2.0*len_tol;  // Needs further tuning
  double tol2 = 10.0*len_tol;

  bool changed = true;
  while (changed)
    {
      changed = false;
      int nmb = shell->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> face = shell->getFace(kj);
	  shared_ptr<ParamSurface> surf = face->surface();

	  double len_u, len_v, min_u, max_u, min_v, max_v;
	  surf->estimateSfSize(len_u, min_u, max_u, len_v, min_v, max_v);
	  if (len_u < tol2 || len_v < tol2)
	    {
	      // Candidate sliver face
	      vector<shared_ptr<ftEdge> > edg = face->getAllEdges(0);
	      vector<shared_ptr<Vertex> > corners = 
		      face->getCornerVertices(toptol.bend);

	      vector<ftSurface*> adj_faces;
	      face->getAdjacentFaces(adj_faces);
	      vector<vector<shared_ptr<ftEdge> > > edges(adj_faces.size());
	      vector<double> edges_len(adj_faces.size());
	      for (size_t ki=0; ki<adj_faces.size(); ++ki)
		{
		  vector<shared_ptr<ftEdge> > tmp_edg = 
		    face->getCommonEdges(adj_faces[ki]);
		  edges[ki] = tmp_edg;
		  double len = 0.0;
		  for (size_t kr=0; kr<edges[ki].size(); ++kr)
		    {
		      len += edges[ki][kr]->estimatedCurveLength();
		    }
		  edges_len[ki] = len;
		}
	      

	      // Start simple
	      if (edg.size() == corners.size() && 
		  (edg.size() == 3 || edg.size() == 4))
		{
		  // Compute edge lengths
		  vector<double> edg_len(edg.size());
		  int nmb_small = 0, nmb_small2 = 0;
		  int ix1 = -1, ix2 = -1;
		  for (size_t ki=0; ki<edg.size(); ++ki)
		    {
		      edg_len[ki] = edg[ki]->estimatedCurveLength();
		      if (edg_len[ki] < tol1)
			nmb_small++;
		      else if (ix1 < 0)
			ix1 = (int)ki;
		      else
			ix2 = (int)ki;
		      if (edg_len[ki] < tol2)
			nmb_small2++;
		    }

		  bool removed = false;
		  vector<shared_ptr<ParamSurface> > modified_sfs;
		  if (nmb_small == (int)edg.size())
		    {
		      // A small face, not a sliver face
		      removed = smallFace(face, edg, modified_sfs);
		    }
		  else if (ix1 >= 0 && ix2 >= 0)
		    {
		      removed = sliverFace(face, edg, ix1, ix2, nmb_small, nmb_small2,
					   modified_sfs);
		    }

		  if (removed)
		    {
		      // Update trimming shell
		      // First remove sliver face and adjacent faces
		      for (size_t ki=0; ki<edg.size(); ++ki)
			{
			  shared_ptr<ftSurface> face2 = 
			    shell->fetchAsSharedPtr(edg[ki]->twin()->face());
			  shell->removeFace(face2);
			}
		      Body *bd = face->getBody();
		      shell->removeFace(face);

		      // Create and append modified faces
		      for (size_t ki=0; ki<modified_sfs.size(); ++ki)
			{
			  shared_ptr<ftSurface> mod_face(new ftSurface(modified_sfs[ki], -1));
			  mod_face->setBody(bd);
			  shell->append(mod_face, false);
			}
#ifdef DEBUG_VOL1
		      std::ofstream of("mod_sliver.g2");
		      vector<shared_ptr<Vertex> > vx;

		      int nmb_face = shell->nmbEntities();
		      for (int ka=0; ka<nmb_face; ++ka)
			{
			  shared_ptr<ParamSurface> tmp = shell->getSurface(ka);
			  tmp->writeStandardHeader(of);
			  tmp->write(of);
			}
		      vector<shared_ptr<Vertex> > tmp_vx;
		      shell->getAllVertices(tmp_vx);
		      vx.insert(vx.end(), tmp_vx.begin(), tmp_vx.end());
		      of << "400 1 0 4 255 0 0 255" << std::endl;
		      of << vx.size() << std::endl;
		      for (size_t ki=0; ki<vx.size(); ++ki)
			of << vx[ki]->getVertexPoint() << std::endl;
#endif
		      changed = true;
		      break;
		    }
		}
	    }
	  int stop_break = 1;
	}
    }
}

//===========================================================================
// 
// 
bool ftVolume::smallFace(shared_ptr<ftSurface> face,
			 vector<shared_ptr<ftEdge> >& edg,
			 vector<shared_ptr<ParamSurface> >& modified_sfs)
//===========================================================================
{
  bool modified = false;

  // Only for triangular faces
  if (edg.size() != 3)
    return false;

  tpTolerances toptol = shells_[0]->getTolerances();

  // Configuration 1: the identified faces originates from the trimming shell
  // while all adjacent faces are element faces
  // Check if the current setting belongs to configuration 1
  size_t ki, kj;
  vector<ftSurface*> adj_faces(edg.size());
  for (ki=0; ki<edg.size(); ++ki)
    {
      adj_faces[ki] = edg[ki]->twin()->face()->asFtSurface();
      if (!adj_faces[ki])
	return false;
    }

  vector<shared_ptr<SurfaceOnVolume> > vol_sf(adj_faces.size());
  vector<shared_ptr<BoundedSurface> > bd_sf(adj_faces.size());
  for (ki=0; ki<adj_faces.size(); ++ki)
    {
      shared_ptr<ParamSurface> sf = adj_faces[ki]->surface();
      bd_sf[ki] = dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf);
      vol_sf[ki] = (bd_sf[ki].get()) ? 
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf[ki]->underlyingSurface()) :
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(sf);	
      if (!vol_sf[ki].get())
	break;

      if (!vol_sf[ki]->atBoundary())
	break;
    }
  if (ki < adj_faces.size())
    {
      // Not configuration 1
      return false;
    }
        
  // Spline volume
  shared_ptr<SplineVolume> vol = dynamic_pointer_cast<SplineVolume,ParamVolume>(vol_);
  if (!vol.get())
    return false;

  double tol = toptol.gap;
  int orientation;
  bool swap;
  vector<int> adj_bd(vol_sf.size());
  for (ki=0; ki<vol_sf.size(); ++ki)
    {
      adj_bd[ki] = vol_sf[ki]->whichBoundary(tol, orientation, swap);
      if (adj_bd[ki] < 0)
	break;
    }
  if (ki < vol_sf.size())
    return false;

  // For each combination of adjacent surfaces, fetch the non-trimmed boundary
  // curve between them
  vector<vector<shared_ptr<ParamCurve> > > mod_cvs(vol_sf.size());
  for (ki=0; ki<vol_sf.size(); ++ki)
    for (kj=ki+1; kj<vol_sf.size(); ++kj)
      {
	int minbd = std::min(adj_bd[ki], adj_bd[kj]);
	int maxbd = std::max(adj_bd[ki], adj_bd[kj]);
	int bd_ix = 2*(maxbd-2) + minbd + 2*(minbd>=2);
	shared_ptr<ParamCurve> geomcv, parcv;
	vol->getBoundaryCurve(bd_ix, geomcv, parcv);
	shared_ptr<ParamCurve> curr_cv =
	   shared_ptr<ParamCurve>(new CurveOnVolume(vol_, parcv, geomcv, false));
	mod_cvs[ki].push_back(curr_cv);
	mod_cvs[kj].push_back(shared_ptr<ParamCurve>(curr_cv->clone()));
      }

  vector<int> not_changed;
  for (ki=0; ki<bd_sf.size(); ++ki)
    {
      if (bd_sf[ki].get())
	{
	  shared_ptr<BoundedSurface> updated_sf =
	    replaceBdCvs(bd_sf[ki], mod_cvs[ki], tol, toptol_.neighbour);
	  modified_sfs.push_back(updated_sf);
	  modified = true;
	}
      else
	{
	  not_changed.push_back((int)ki);
	}
    }
  for (int kr=(int)not_changed.size()-1; kr>=0; --kr)
    edg.erase(edg.begin()+kr);

  // Configuration 2: the identified faces is an element face and at least one
  // adjacent face originates from the trimming shell
  // Configuration 2 is currently not handled


  return modified;
}

//===========================================================================
// 
// 
bool ftVolume::sliverFace(shared_ptr<ftSurface> face,
			  vector<shared_ptr<ftEdge> >& edg,
			  int ix1, int ix2, int nmb_small, int nmb_small2,
			  vector<shared_ptr<ParamSurface> >& modified_sfs)
//===========================================================================
{
  bool removed = false;
  tpTolerances toptol = shells_[0]->getTolerances();
  ftSurface *adj_face1 = edg[ix1]->twin()->face()->asFtSurface();
  ftSurface *adj_face2 = edg[ix2]->twin()->face()->asFtSurface();
  bool smooth;
  bool is_adjacent = adj_face1->isAdjacent(adj_face2, smooth);
  if (is_adjacent)
    {
      vector<shared_ptr<ftEdge> > common = 
	adj_face1->getCommonEdges(adj_face2);
      if (common.size() == 1)
	{
	  // Check edge continuity
	  tpJointType tp1, tp2;
	  if (edg[ix1]->twin()->next() == common[0].get())
	    tp1 = edg[ix1]->twin()->checkContinuity(common[0].get(),
						    toptol.neighbour,
						    toptol.gap,
						    toptol.bend,
						    toptol.kink);
	  else
	    tp1 = common[0]->checkContinuity(edg[ix1]->twin(),
					     toptol.neighbour,
					     toptol.gap,
					     toptol.bend,
					     toptol.kink);
	  if (edg[ix2]->twin()->next() == common[0]->twin())
	    tp2 = edg[ix2]->twin()->checkContinuity(common[0]->twin(),
						    toptol.neighbour,
						    toptol.gap,
						    toptol.bend,
						    toptol.kink);
	  else
	    tp2 = common[0]->twin()->checkContinuity(edg[ix2]->twin(),
						     toptol.neighbour,
						     toptol.gap,
						     toptol.bend,
						     toptol.kink);
	  if (tp1 <= JOINT_KINK || tp2 <= JOINT_KINK)
	    is_adjacent = false;  // Continue with sliver face removal
	}
    }
  if (nmb_small == (int)edg.size()-2 &&
      nmb_small == nmb_small2 &&
      !(edg.size() == 4 && ix2-ix1 == 1) &&
      !is_adjacent)
    // (edg_len[0] < tol1 && edg_len[2] < tol1) ||
    // (edg_len[1] < tol1 && edg_len[3] < tol1))
    {
      // Sliver face, check possibility of removal
      shared_ptr<ParamSurface> sf1 = adj_face1->surface();
      shared_ptr<ParamSurface> sf2 = adj_face2->surface();
      shared_ptr<BoundedSurface> bd_sf1 = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf1);
      shared_ptr<BoundedSurface> bd_sf2 = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf2);
      shared_ptr<SurfaceOnVolume> vol_sf1 = (bd_sf1.get()) ? 
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf1->underlyingSurface()) :
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(sf1);		  
      shared_ptr<SurfaceOnVolume> vol_sf2 = (bd_sf2.get()) ? 
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf2->underlyingSurface()) :
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(sf2);
#ifdef DEBUG_VOL1
      std::ofstream of0("sliver.g2");
      shared_ptr<ParamSurface> tmp = face->surface();
      tmp->writeStandardHeader(of0);
      tmp->write(of0);
#endif
      if (vol_sf1.get() && vol_sf1->atBoundary() &&
	  vol_sf2.get() && vol_sf2->atBoundary())
	{
	  // Promising configuration for sliver face removal
	  removed = removeSliver1(face, edg, ix1, ix2, toptol.gap,
				  modified_sfs);
	}
      else if (vol_sf1.get() && vol_sf2.get() && 
	       (vol_sf1->atBoundary() || vol_sf2->atBoundary()))
	{
	  shared_ptr<ftEdge> not_changed;
	  removed = removeSliver2(face, edg, ix1, ix2, toptol.gap,
				  modified_sfs, not_changed);
	  for (size_t ki=0; ki<edg.size(); ++ki)
	    {
	      if (edg[ki].get() == not_changed.get())
		{
		  edg.erase(edg.begin()+ki);
		  break;
		}
	    }
	}
    }
  return removed;
}

//===========================================================================
// 
// 
void ftVolume::splitElementByTrimSfs(int elem_ix, double eps,
				     vector<shared_ptr<ftVolume> >& sub_elem,
				     vector<int>& is_inside)
//===========================================================================
{
  if (!isSpline())
    return;

  shared_ptr<SplineVolume> vol = dynamic_pointer_cast<SplineVolume>(vol_);
  if (!vol.get())
    return;
  
  // Fetch parameter values surrounding the specified element
  double elem_par[6];
  vol->getElementBdPar(elem_ix, elem_par);

  // Create an ftVolume entity corresponding to the element. 
  // First create underlying SplineVolume
  shared_ptr<ParamVolume> elem_vol(vol->subVolume(elem_par[0], elem_par[2],
						  elem_par[4], elem_par[1],
						  elem_par[3], elem_par[5]));

  vector<shared_ptr<ftSurface> > elem_faces = 
    getBoundaryFaces(elem_vol, toptol_.gap, toptol_.kink);

  sub_elem = ftVolumeTools::splitElement(elem_vol, elem_faces, 
					 elem_par, this, eps, is_inside);

  // Remove sliver faces
  for (size_t ki=0; ki<sub_elem.size(); ++ki)
    {
      if (is_inside[ki])  // Should include parameter
	(void)sub_elem[ki]->removeSliverFaces(toptol_.neighbour);
    }

  // shared_ptr<ftVolume> elem_vol2(new ftVolume(elem_vol, toptol_.gap,
  // 					      toptol_.kink, -1));

  // sub_elem = ftVolumeTools::splitOneVol(elem_vol2, this, eps, is_inside, 
  // 					elem_par, 6);
}

//===========================================================================
// 
// 
int ftVolume::ElementOnBoundary(int elem_ix)
//===========================================================================
{
  if (!isSpline())
    return -1;

  shared_ptr<SplineVolume> vol = dynamic_pointer_cast<SplineVolume>(vol_);
  if (!vol.get())
    return -1;
  
  // Fetch surfaces surrounding the specified element
  double elem_par[6];
  vector<shared_ptr<SplineSurface> > side_sfs = vol->getElementBdSfs(elem_ix, 
								     elem_par);


  vector<shared_ptr<SurfaceModel> > shells = getAllShells();
#ifdef DEBUG
  std::ofstream mod("elem_trim.g2");
  for (size_t ka=0; ka<shells.size(); ++ka)
    {
      int nmb = shells[ka]->nmbEntities();
      for (int kr=0; kr<nmb; ++kr)
	{
	  shared_ptr<ParamSurface> sf = shells[ka]->getSurface(kr);
	  sf->writeStandardHeader(mod);
	  sf->write(mod);
	}
    }
  for (int kr=0; kr<(int)side_sfs.size(); ++kr)
    {
      side_sfs[kr]->writeStandardHeader(mod);
      side_sfs[kr]->write(mod);
    }
#endif

  // Can use a tolerance specifically suited for the intersections here
  // as we do not want the exact intersection curve, but only an indication
  // if it is any intersections
  double eps = 1.0e-6; //toptol_.gap;
  for (size_t kj=0; kj<shells.size(); ++kj)
    {
      int nmb = shells[kj]->nmbEntities();
      shared_ptr<ParamSurface> surf;
      for (int kh=0; kh<nmb; ++kh)
	{
	  shared_ptr<ftSurface> face = shells[kj]->getFace(kh);
	  int bd_status = ftVolumeTools::boundaryStatus(this, face, eps);
	  if (bd_status >= 0)
	    continue;  // Not a trimming face
	  surf = face->surface();
	  BoundingBox box = surf->boundingBox();
	  
	  // Check if the surface already is defined as an element boundary 
	  // surface, i.e. has constant parameter equal to element boundary 
	  // parameter
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(surf);
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
	  if (bd_sf.get())
	    vol_sf = 
	      dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bd_sf->underlyingSurface());
	  int dir = 0;
	  double val = 0.0;
	  if (vol_sf.get())
	    {
	      dir = vol_sf->getConstDir();
	      val = vol_sf->getConstVal();
	    }

	  for (size_t ki=0; ki<side_sfs.size(); ++ki)
	    {
	      BoundingBox box2 = side_sfs[ki]->boundingBox();
	      if (!box.overlaps(box2))
		continue;

	      if (dir == ((int)ki/2) + 1 && fabs(val-elem_par[ki]) < eps)
		continue;  // Coincidence

	      shared_ptr<BoundedSurface> bd1, bd2;
	      vector<shared_ptr<CurveOnSurface> > int_cv1, int_cv2;
	      BoundedUtils::getSurfaceIntersections(surf, side_sfs[ki], eps,
						    int_cv1, bd1,
						    int_cv2, bd2);
	      if (int_cv1.size() > 0 || int_cv2.size() > 0)
		return 1;
	    }
	}

      // Check if the trimming surface is completely inside the element
      if (nmb > 0)
	{
	  int kr;
	  for (kr=0; kr<nmb; ++kr)
	    {
	      surf = shells[kj]->getSurface(kr);
	      double upar, vpar;
	      Point in_pt = surf->getInternalPoint(upar, vpar);
	      
	      double u, v, w, d;
	      Point close_pt;
	      vol_->closestPoint(in_pt, u, v, w, close_pt, d, eps);
	      if (u < elem_par[0] || u > elem_par[1] ||
		  v < elem_par[2] || v > elem_par[3] ||
		  w < elem_par[4] || w > elem_par[5])
		return 0;
	    }
	  return 1;
	}
    }

  return 0;
}

//===========================================================================
// 
// 
int ftVolume::ElementBoundaryStatus(int elem_ix) 
//===========================================================================
{
  // Result: -1 = not a spline volume, 0 = outside, 1 = on boundary, 2 = inside

  if (!isSpline())
    return -1;
  
  int bdstat = ElementOnBoundary(elem_ix);
  if (bdstat != 0)
    return bdstat;

  // Fetch point internal to element
  // First fetch associated spline volume
  shared_ptr<SplineVolume> vol = dynamic_pointer_cast<SplineVolume>(vol_);
  if (!vol.get())
    return -1;
  
  // Fetch number of patches in all parameter directions
  int nu = vol->numberOfPatches(0);
  int nv = vol->numberOfPatches(1);
  int nw = vol->numberOfPatches(2);

  if (elem_ix < 0 || elem_ix >= nu*nv*nw)
    return 0;

  // 3-variate index
  int iw = elem_ix/(nu*nv);
  int iv = (elem_ix - iw*nu*nv)/nu;
  int iu = elem_ix - iw*nu*nv - iv*nu;

  // Parameter value
  vector<double> knots_u;
  vector<double> knots_v;
  vector<double> knots_w;
  const BsplineBasis basis_u = vol->basis(0);
  const BsplineBasis basis_v = vol->basis(1);
  const BsplineBasis basis_w = vol->basis(2);
  basis_u.knotsSimple(knots_u);
  basis_v.knotsSimple(knots_v);
  basis_w.knotsSimple(knots_w);

  double upar = 0.5*(knots_u[iu]+knots_u[iu+1]);
  double vpar = 0.5*(knots_v[iv]+knots_v[iv+1]);
  double wpar = 0.5*(knots_w[iw]+knots_w[iw+1]);
  
  // Check if the element is inside the trimming loop(s)
  Point pnt;
  vol->point(pnt, upar, vpar, wpar);
#ifdef DEBUG
  std::ofstream of("inside_test.g2");
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = shells_[ki]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "1" << std::endl;
  of << pnt << std::endl;
#endif
  bool inside = isInside(pnt);
  return (inside) ? 2 : 0;
}

//===========================================================================
// 
// 
bool ftVolume::isBoundaryTrimmed() const
//===========================================================================
{
  if (nmbOfShells() != 1)
    return false;  // Inner trimmming exist

  // Get boundary surfaces
  shared_ptr<SurfaceModel> shell = getOuterShell();

  // For each surface, check if it is a boundary surface of this volume
  int nmb = shell->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> surf = shell->getSurface(ki);

      shared_ptr<SurfaceOnVolume> vol_sf = getVolSf(surf);

      if (vol_sf.get())
	{
	  // Check if it is a constant boundary surface
	  int orientation;
	  bool swap;
	  int bd = vol_sf->whichBoundary(0.0, orientation, swap);
	  if (bd < 0)
	    return false;   // Surface do not follow a boundary

	  // Check if the boundary surface is constant in this volume
	  shared_ptr<ParamVolume> vol = vol_sf->getVolume();
	  if (vol.get() != vol_.get())
	    return false;
	}
      else 
	return false;  // No surface on volume is found
    }

  return true;  // All boundary surfaces are boundary surfaces of the 
  // non-trimmed volume
}

//===========================================================================
// 
// 
bool ftVolume::isIsoTrimmed() const
//===========================================================================
{
  if (nmbOfShells() != 1)
    return false;  // Inner trimmming exist

  // Get boundary surfaces
  shared_ptr<SurfaceModel> shell = getOuterShell();

  // For each surface, check if it is a boundary surface of this volume
  int nmb = shell->nmbEntities();
  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> surf = shell->getSurface(ki);

      shared_ptr<SurfaceOnVolume> vol_sf = getVolSf(surf);

      if (vol_sf.get())
	{
	  // Check if it is a constant iso surface
	  int dir = vol_sf->getConstDir();
	  if (dir < 1 || dir > 3)
	    return false;   // Surface do not follow an iso-surface

	  // Check if the surface is constant in this volume
	  shared_ptr<ParamVolume> vol = vol_sf->getVolume();
	  if (vol.get() != vol_.get())
	    return false;
	}
      else 
	return false;  // No surface on volume is found
    }

  return true;  // All surfaces follow an iso-parameter is the volume
}

//===========================================================================
// 
// 
bool ftVolume::ParamInVolume(double u, double v, double w)
//===========================================================================
{
  Point pt;
  vol_->point(pt, u, v, w);
#ifdef DEBUG
  std::ofstream of("pt_in_vol.g2");
  of << "400 1 0 4 0 255 0 255" << std::endl;
  of << "1" << std::endl;
  of << pt << std::endl;
#endif
  return isInside(pt);
}

//===========================================================================
// 
// 
bool ftVolume::regularizeBdShells(vector<pair<Point,Point> >& corr_vx_pts,
				  vector<SurfaceModel*>& modified_adjacent,
				  int split_mode, bool pattern_split,
				    int level)
//===========================================================================
{
  bool updated = false;
  int nmb_shells = nmbOfShells();
  for (int ki=0; ki<nmb_shells; ++ki)
    {
      // Fetch current boundary shell
      shared_ptr<SurfaceModel> sfmodel = getShell(ki);

#ifdef DEBUG_VOL1
      std::ofstream mod("pre_reg.g2");
      vector<shared_ptr<Vertex> > vxs;
      sfmodel->getAllVertices(vxs);
      int nmb = sfmodel->nmbEntities();
      for (int kr=0; kr<nmb; ++kr)
	{
	  shared_ptr<ParamSurface> sf = sfmodel->getSurface(kr);
	  sf->writeStandardHeader(mod);
	  sf->write(mod);
	}
      mod << std::endl << "400 1 0 4 255 0 0 255" << std::endl;
      mod << vxs.size() << std::endl;
      for (int kr=0; kr<(int)vxs.size(); ++kr)
	mod << vxs[kr]->getVertexPoint() << std::endl;
#endif

      RegularizeFaceSet regularize0(sfmodel, pattern_split, level);
      
      // Simplify if possible
      int degree = 3;
      if (level == 0)
	{
	  regularize0.removeExtraDiv(true);
	  mergeSmoothJoints(degree);
	}

      // Fetch info about opposite surfaces
      vector<pair<int, int> > opposite = oppositeSfs(sfmodel);

      // Regularize shell
      int nmb_faces = sfmodel->nmbEntities();
     RegularizeFaceSet regularize(sfmodel, pattern_split, level);
      regularize.setSplitMode(split_mode);
      size_t kj;
      for (kj=0; kj<opposite.size(); ++kj)
	regularize.setFaceCorrespondance(opposite[kj].first, 
					 opposite[kj].second);
      for (kj=0; kj<modified_adjacent.size(); ++kj)
	if (modified_adjacent[ki] == sfmodel.get())
	  break;

      bool reverse_sequence = (kj < modified_adjacent.size());
      shared_ptr<SurfaceModel> sfmodel2 = 
	regularize.getRegularModel(reverse_sequence);
      corr_vx_pts = regularize.fetchVxPntCorr();
      if (sfmodel2->nmbEntities() > nmb_faces)
	{
	  removeSeamFaces();
	  updated = true;
	}
      
      vector<SurfaceModel*> curr_adjacent = 
	regularize.getModifiedAdjacentModels();
      modified_adjacent.insert(modified_adjacent.end(), 
			       curr_adjacent.begin(), curr_adjacent.end());
    }
  return updated;
}

//===========================================================================
// 
// 
bool ftVolume::isRegularized(bool accept_degen) const
//===========================================================================
{
  if (nmbOfShells() != 1)
    return false;  // Inner trimmming exist

  // Get boundary surfaces
  shared_ptr<SurfaceModel> shell = getOuterShell();

  // The number of surfaces must be 6 and each of them must be 4-sided
  // NB! Degenerate situations may be considered later
  int nmb = shell->nmbEntities();
  if ((nmb != 6 && accept_degen == false) || 
      ((nmb > 6 || nmb < 4) && accept_degen))
    return false;

  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ftSurface> face = shell->getFace(ki);

      if (!face->onlyOuterTrim())
	return false;  // Cannot be approximated by one spline surface
      int nmb_bd = face->nmbOuterBdCrvs(toptol_.gap, toptol_.neighbour,
					toptol_.bend);
      vector<ftSurface*> adj_faces;
      face->getAdjacentFaces(adj_faces);
      if (nmb_bd > (int)adj_faces.size() || nmb_bd == 0)
	nmb_bd = (int)adj_faces.size();   // Large angle can be caused by sliver face removal
      if (adj_faces.size() > 4 || (nmb_bd != 4 && accept_degen == false) ||
	  ((nmb_bd > 4 || nmb_bd < 3) && accept_degen))
	{
#ifdef DEBUG
	  std::ofstream of("nonreg_face.g2");
	  shared_ptr<ParamSurface> tmp_sf = face->surface();
	  tmp_sf->writeStandardHeader(of);
	  tmp_sf->write(of);
#endif
	  return false;
	}
    }

  return true;  // The volume may be described without trimming
}

//===========================================================================
// 
// 
bool ftVolume::untrimRegular(int degree, bool accept_degen) 
//===========================================================================
{
  // Check configuration
  if (shells_.size() != 1)
    return false;  // Not regular
  
  if ((shells_[0]->nmbEntities() < 4 && shells_[0]->nmbEntities() > 6) ||
      (accept_degen == false && shells_[0]->nmbEntities() != 6))
    return false;  // Not regular or trimmed

#ifdef DEBUG_VOL1
  bool isOK = shells_[0]->checkShellTopology();
  std::cout << "Shell topology: " << isOK << std::endl;
#endif

  // Sort surfaces in shell according to volume configuration
  // Sequence: umin, umax, vmin, vmax, wmin, wmax
  vector<shared_ptr<ParamSurface> > sorted_sfs(6);
  vector<std::pair<int,double> > classification(6);
  vector<int> deg_type(6);
  bool sorted = sortRegularSurfaces(sorted_sfs, classification, deg_type);
  if (!sorted)
    return false;  // Sorting failed

#ifdef DEBUG_VOL1
  std::ofstream of0("sorted_sfs.g2");
#endif

  int ki, kj;
#ifdef DEBUG_VOL1
  for (ki=0; ki<6; ++ki)
    {
      if (!sorted_sfs[ki].get())
	MESSAGE("Missing volume side surface, nr " << ki+1);
      //return false;  // Unsufficient information
      else
	{
	  sorted_sfs[ki]->writeStandardHeader(of0);
	  sorted_sfs[ki]->write(of0);
	}
    }
#endif

  shared_ptr<ParamVolume> vol2;

  // Check if all trimming surfaces are iso parametric
  for (ki=0; ki<6; ++ki)
    {
      if (classification[ki].first == 0)
	break;  // Not iso parametric
    }
  
  bool loft_sequence = false;  
  if (ki == 6)
    {
      // The non-trimmed volume can be achieved by subdivision
      shared_ptr<SplineVolume> vol = dynamic_pointer_cast<SplineVolume>(vol_);
      if (vol.get())
	{
	  vol2 = 
	    shared_ptr<ParamVolume>(vol->subVolume(classification[0].second,
						   classification[2].second,
						   classification[4].second,
						   classification[1].second,
						   classification[3].second,
						   classification[5].second));
	}
    }
  else
    {
      // Check if the new volume is linear in one parameter direction
      // and can be made by lofting
      // Check order in the three parameter directions, -1 = not known
      if (vol_->isSpline())
	{
	  shared_ptr<SplineVolume> svol = 
	    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol_);
	  int order[] = {-1, -1, -1};
	  for (ki=0; ki<3; ki++)
	    {
	      for (kj=0; kj<6; ++kj)
		{
		  if (kj == 2*ki || kj == 2*ki+1)
		    continue;
		  if (classification[kj].first == 0)
		    break;
		}
	      if (kj == 6)
		order[ki] = svol->order(ki);
	    }
	  
	  // Check if unclassified surfaces exist in only one parameter 
	  // direction
	  int dir = -1;
	  for (ki=0; ki<6; ki+=2)
	    {
	      if (classification[ki].first == 0 || 
		  classification[ki+1].first == 0)
		{
		  if (dir < 0)
		    dir = ki/2;
		  else
		    break;
		}
	    }
	  if (false /*ki == 6 && dir >= 0 && order[dir] == 2*/)
	    {
	      // VSK 0817. Must improve orientation of lofting surfaces
	      // Create volume by loft
	      vol2 = createByLoft(sorted_sfs[2*dir], sorted_sfs[2*dir+1],
				  shells_[0]->getTolerances().gap, dir);
	      loft_sequence = true;
	    }
	}

      // Use Coons volume to create the non-trimmed volume
      if (!vol2.get())
	{
	  vol2 = createByCoons(sorted_sfs, classification, deg_type,
			       shells_[0]->getTolerances().gap,
			       degree);
	}
    }

  if (vol2.get())
    {
#ifdef DEBUG_VOL1
  std::ofstream of("mod_vol.g2");
  vol2->writeStandardHeader(of);
  vol2->write(of);
#endif

  // Update current ftVolume with new parametric volume. Update also
  // boundary surfaces
  replaceParamVolume(vol2, sorted_sfs, loft_sequence);

  return true;
    }
  else
    return false;
}


//===========================================================================
// 
// 
shared_ptr<ParamVolume> ftVolume::getRegParVol(int degree, int bd_cond[6][2],
					       bool accept_degen) 
//===========================================================================
{
  shared_ptr<ParamVolume> result;

  // Initiate to no boundary conditions
  int ki;
  for (ki=0; ki<6; ++ki)
    bd_cond[ki][0] = bd_cond[ki][1] = -1;

  // Check configuration
  if (shells_.size() != 1)
    return result;  // Not regular
  
  if ((shells_[0]->nmbEntities() < 4 && shells_[0]->nmbEntities() > 6) ||
      (accept_degen == false && shells_[0]->nmbEntities() != 6))
    return result;  // Not regular or trimmed

#ifdef DEBUG_VOL1
  bool isOK = shells_[0]->checkShellTopology();
  std::cout << "Shell topology: " << isOK << std::endl;
#endif

  // Sort surfaces in shell according to volume configuration
  // Sequence: umin, umax, vmin, vmax, wmin, wmax
  vector<shared_ptr<ParamSurface> > sorted_sfs(6);
  vector<std::pair<int,double> > classification(6);
  vector<int> deg_type(6);
  bool sorted = sortRegularSurfaces(sorted_sfs, classification, deg_type);
  if (!sorted)
    return result;  // Sorting failed

#ifdef DEBUG_VOL1
  std::ofstream of0("sorted_sfs.g2");

  for (ki=0; ki<6; ++ki)
    {
      if (!sorted_sfs[ki].get())
	MESSAGE("Missing volume side surface, nr " << ki+1);
      else
	{
	  sorted_sfs[ki]->writeStandardHeader(of0);
	  sorted_sfs[ki]->write(of0);
	}
    }
#endif

  // Use Coons volume to create the non-trimmed parameter volume
  // First replace the side surfaces by their parameter volume
  // equivalents
  vector<shared_ptr<ParamSurface> > sorted_sfs2(6);
  vector<shared_ptr<ftSurface> > face2(6);
  for (ki=0; ki<(int)sorted_sfs.size(); ++ki)
    {
      shared_ptr<SurfaceOnVolume> volsf = 
	dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(sorted_sfs[ki]);
      shared_ptr<BoundedSurface> bdsf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(sorted_sfs[ki]);
      shared_ptr<ParamSurface> volsf2;
      // Represent the surface as a parameter based surface on volume
      if (bdsf.get())
	  volsf = 
	    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bdsf->underlyingSurface());
      if (volsf.get())
	{
	  int bd, orient;
	  bool swap;
	  bd = volsf->whichBoundary(toptol_.gap, orient, swap);
	  
	  volsf2 = 
	    shared_ptr<ParamSurface>(new ParameterSurfaceOnVolume(vol_,
								  volsf->parameterSurface(),
								  volsf->spaceSurface(),
								  volsf->getConstDir(),
								  volsf->getConstVal(),
								  bd, swap));
	}
      else if (bdsf.get())
	volsf2 = shared_ptr<ParamSurface>(new ParameterSurfaceOnVolume(vol_,
								       bdsf->underlyingSurface()));
      else
	{
	  volsf2 = shared_ptr<ParamSurface>(new ParameterSurfaceOnVolume(vol_, sorted_sfs[ki]));
	}

      // Must define the bounded surface to get consistent information
      // for adjacency analysis
      // Fetch boundary curves of initial surface
      if (sorted_sfs[ki].get())
	{
	  vector<pair<int,double> > split_param;
	  if (sorted_sfs[ki]->instanceType() != Class_BoundedSurface)
	    {
	      // Replace with bounded surface. First find corresponding face
	      shared_ptr<SurfaceModel> shell = getOuterShell();
	      int sf_ix = shell->getIndex(sorted_sfs[ki].get());
	      shared_ptr<ftSurface> curr_face = shell->getFace(sf_ix);

	      vector<shared_ptr<ftEdge> > curr_edgs = curr_face->getAllEdges(0);
	      vector<shared_ptr<CurveOnSurface> > curr_cvs;
	      ParamCurve* prev = NULL;
	      for (size_t kj=0; kj<curr_edgs.size(); ++kj)
		{
		  shared_ptr<ParamCurve> cv = curr_edgs[kj]->geomCurve();
		  if (cv.get() == prev)
		    {
		      // Same underlying curve. Must remember split parameter for edge
		      split_param.push_back(make_pair(kj, curr_edgs[kj]->tMin()));
		    }
		  else
		    curr_cvs.push_back(dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv));
		  prev = cv.get();
		}

	      shared_ptr<BoundedSurface> curr_sf(new BoundedSurface(sorted_sfs[ki], 
								    curr_cvs,
								    shells_[0]->getTolerances().gap,
								    false));
	      curr_sf->setIsoTrim();
	      
	      curr_face->replaceSurf(curr_sf);
	      sorted_sfs[ki] = curr_sf;
	    }
	  vector<CurveLoop> bd_loops = 
	    SurfaceTools::allBoundarySfLoops(sorted_sfs[ki], DEFAULT_SPACE_EPSILON);
	  vector<vector<shared_ptr<CurveOnSurface> > > loops;
	  for (size_t kj=0; kj<bd_loops.size(); ++kj)
	    {
	      int nmb = bd_loops[kj].size();
	      vector<shared_ptr<CurveOnSurface> > curr_loop(nmb);
	      for (int kr=0; kr<nmb; ++kr)
		{
		  shared_ptr<CurveOnSurface> bd_cv = 
		    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(bd_loops[kj][kr]);
		  if (!bd_cv.get())
		    return result;
		  shared_ptr<ParameterCurveOnVolume> volcv;
		  shared_ptr<CurveOnVolume> volcv0 = 
		    dynamic_pointer_cast<CurveOnVolume, ParamCurve>(bd_cv->spaceCurve());
		  if (volcv0.get())
		    volcv = shared_ptr<ParameterCurveOnVolume>(new ParameterCurveOnVolume(vol_, volcv0->spaceCurve()));
		  else
		    volcv = shared_ptr<ParameterCurveOnVolume>(new ParameterCurveOnVolume(vol_, bd_cv->spaceCurve()));
		  curr_loop[kr] = shared_ptr<CurveOnSurface>(new CurveOnSurface(volsf2, bd_cv->parameterCurve(), volcv, false, 0, 0, 0.0, -1, true));
		}
	      loops.push_back(curr_loop);
	    }
	  double eps = bd_loops[0].getSpaceEpsilon();
	  sorted_sfs2[ki] = shared_ptr<ParamSurface>(new BoundedSurface(volsf2,
									loops,
									eps,
									false));

	  // Create also the corresponding face
	  face2[ki] = shared_ptr<ftSurface>(new ftSurface(sorted_sfs2[ki], 
							  (int)ki));
	  (void)face2[ki]->createInitialEdges(DEFAULT_SPACE_EPSILON, 0.00015, true);
	}
    }

  // Transfer adjacency information to the set of parameter volume
  // side surfaces
  setParameterVolAdjacency(sorted_sfs, face2);

  // Set boundary information
  for (ki=0; ki<6; ++ki)
    {
      if (face2[ki].get() && face2[ki]->hasBoundaryConditions())
	{
	  int bd_type, bd;
	  face2[ki]->getBoundaryConditions(bd_type, bd);
	  bd_cond[ki][0] = bd_type;
	  bd_cond[ki][1] = bd;
	}
    }

  // Make sure that there are no missing faces
  for (ki=(int)face2.size()-1; ki>=0; --ki)
    if (!face2[ki].get())
      face2.erase(face2.begin()+ki);
  
  // Create intermediate topological volume
  // shared_ptr<SurfaceModel> shell(new SurfaceModel(toptol_.gap, toptol_.gap,
  // 						  10.0*toptol_.neighbour, 
  // 						  toptol_.kink, toptol_.bend,
  // 						  sorted_sfs2));
  shared_ptr<SurfaceModel> shell(new SurfaceModel(toptol_.gap, toptol_.gap,
						  toptol_.neighbour, 
						  toptol_.kink, toptol_.bend,
						  face2, true));
  shared_ptr<ftVolume> tmp_vol(new ftVolume(vol_, shell));

  result = tmp_vol->createByCoons(sorted_sfs2, classification, deg_type,
				  shells_[0]->getTolerances().gap,
				  degree, false);

#ifdef DEBUG_VOL1
  if (result.get())
    {
      std::ofstream of("mod_vol.g2");
      result->writeStandardHeader(of);
      result->write(of);
    }
#endif

  return result;
}


//===========================================================================
// 
// 
void 
ftVolume::setParameterVolAdjacency(vector<shared_ptr<ParamSurface> >& sfs1,
				   vector<shared_ptr<ftSurface> >& face2) const
//===========================================================================
{
  // Collect faces corresponding to input surfaces
  vector<shared_ptr<ftSurface> > face1(sfs1.size());
  shared_ptr<SurfaceModel> shell = getOuterShell();
  for (size_t ki=0; ki<sfs1.size(); ++ki)
    {
      if (sfs1[ki].get())
	{
	  int sf_ix = shell->getIndex(sfs1[ki].get());
	  face1[ki] = shell->getFace(sf_ix);

	  // Transfer boundary condition information
	  if (face1[ki]->hasBoundaryConditions())
	    {
	      int bd_type, bd;
	      face1[ki]->getBoundaryConditions(bd_type, bd);
	      face2[ki]->setBoundaryConditions(bd_type, bd);
	    }
	}
    }

  double eps = toptol_.gap;
  for (size_t ki=0; ki<face1.size(); ++ki)
    {
      if (!face1[ki].get())
	continue;
      for (size_t kj=ki+1; kj<face1.size(); ++kj)
	{
	  if (!face2[kj].get())
	    continue;
	  int adj_ix = 0;
	  vector<shared_ptr<ftEdge> > adj1;
	  vector<shared_ptr<ftEdge> > adj2;
	  while (true)
	    {
	      shared_ptr<ftEdge> edge1, edge2;
	      bool adjacent = face1[ki]->areNeighbours(face1[kj].get(), edge1,
						   edge2, adj_ix);
	      if (!adjacent)
		break;

	      adj1.push_back(edge1);
	      adj2.push_back(edge2);

	      // Prepare for more adjacency testing between the current faces
	      adj_ix++;
	    }

	  if (adj1.size() > 0)
	    {
	      // Adjacency. Find corresponding geometry curve and edge in parameter
	      // domain representation
	      vector<shared_ptr<ParamCurve> > cv1;
	      vector<shared_ptr<ParamCurve> > cv2;
	      for (size_t kr=0; kr<adj1.size(); ++kr)
		{
		  shared_ptr<ParamCurve> curr1 = adj1[kr]->geomCurve();
		  shared_ptr<CurveOnSurface> sf_cv1 = 
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(curr1);
		  if (sf_cv1.get())
		    curr1 = sf_cv1->spaceCurve();
		  shared_ptr<CurveOnVolume> vol_cv1 = 
		    dynamic_pointer_cast<CurveOnVolume,ParamCurve>(curr1);
		  if (vol_cv1.get())
		    curr1 = vol_cv1->spaceCurve();
		  cv1.push_back(curr1);

		  shared_ptr<ParamCurve> curr2 = adj2[kr]->geomCurve();
		  shared_ptr<CurveOnSurface> sf_cv2 = 
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(curr2);
		  if (sf_cv2.get())
		    curr2 = sf_cv2->spaceCurve();
		  shared_ptr<CurveOnVolume> vol_cv2 = 
		    dynamic_pointer_cast<CurveOnVolume,ParamCurve>(curr2);
		  if (vol_cv2.get())
		    curr2 = vol_cv2->spaceCurve();
		  cv2.push_back(curr2);
		}

	      // Identify curves cooresponding to the parameter face edges
	      vector<shared_ptr<ftEdge> > e1 = face2[ki]->getAllEdges(0);
	      vector<shared_ptr<ftEdge> > e2 = face2[kj]->getAllEdges(0);

	      vector<shared_ptr<ParamCurve> > pc1;
	      vector<shared_ptr<ParamCurve> > pc2;
	      for (size_t kr=0; kr<e1.size(); ++kr)
		{
		  shared_ptr<ParamCurve> curr = e1[kr]->geomCurve();
		  shared_ptr<CurveOnSurface> sf_cv = 
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(curr);
		  if (sf_cv.get())
		    curr = sf_cv->spaceCurve();
		  shared_ptr<CurveOnVolume> vol_cv = 
		    dynamic_pointer_cast<CurveOnVolume,ParamCurve>(curr);
		  if (vol_cv.get())
		    curr = vol_cv->spaceCurve();
		  pc1.push_back(curr);
		}
	      for (size_t kr=0; kr<e2.size(); ++kr)
		{
		  shared_ptr<ParamCurve> curr = e2[kr]->geomCurve();
		  shared_ptr<CurveOnSurface> sf_cv = 
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(curr);
		  if (sf_cv.get())
		    curr = sf_cv->spaceCurve();
		  shared_ptr<CurveOnVolume> vol_cv = 
		    dynamic_pointer_cast<CurveOnVolume,ParamCurve>(curr);
		  if (vol_cv.get())
		    curr = vol_cv->spaceCurve();
		  pc2.push_back(curr);
		}

	      for (size_t kr=0; kr<cv1.size(); ++kr)
		{
		  // Compute numbers of equal curve pointers
		  int nmb1 = 1, nmb2 = 1;
		  for (size_t kh=0; kh<cv1.size(); ++kh)
		    if (kh != kr && cv1[kh].get() == cv1[kr].get())
		      nmb1++;
		  for (size_t kh=0; kh<cv2.size(); ++kh)
		    if (kh != kr && cv2[kh].get() == cv2[kr].get())
		      nmb2++;

		  // Find parameter edges corresponding to geometry edges twins
		  size_t ka, kb;
		  for (ka=0; ka<pc1.size(); ++ka)
		    {
		      if (pc1[ka].get() == cv1[kr].get())
			{
			  // if (nmb1 > 1)
			  //   {
			  //     // Make sure to have found the right edge
			  //     Point p1_1 = e1[ka]->point(e1[ka]->tMin());
			  //     Point p1_2 = e1[ka]->point(e1[ka]->tMax());
			  //     Point p2_1 = adj1[kr]->point(adj1[kr]->tMin());
			  //     Point p2_2 = adj1[kr]->point(adj1[kr]->tMax());
			  //     double d1 = std::min(p1_1.dist(p2_1), p1_1.dist(p2_2)) +
			  // 	std::min(p1_2.dist(p2_1), p1_2.dist(p2_2));
			  //     for (size_t kh=ka+1; kh<pc1.size(); ++kh)
			  // 	{
			  // 	  if (pc1[kh].get() != cv1[kr].get())
			  // 	    continue;
			  // 	  Point p1_3 = e1[kh]->point(e1[kh]->tMin());
			  // 	  Point p1_4 = e1[kh]->point(e1[kh]->tMax());
			  // 	  double d2 = std::min(p1_3.dist(p2_1), p1_3.dist(p2_2)) +
			  // 	    std::min(p1_4.dist(p2_1), p1_4.dist(p2_2));
			  // 	  if (d2 < d1)
			  // 	    {
			  // 	      d1 = d2;
			  // 	      ka = kh;
			  // 	    }
			  // 	}
			  //     break;
			  //   }
			  // else
			    break;
			}
		    }
		  for (kb=0; kb<pc2.size(); ++kb)
		    {
		      if (pc2[kb].get() == cv2[kr].get())
			{
			  // if (nmb2 > 1)
			  //   {
			  //     // Make sure to have found the right edge
			  //     Point p1_1 = e2[kb]->point(e2[kb]->tMin());
			  //     Point p1_2 = e2[kb]->point(e2[kb]->tMax());
			  //     Point p2_1 = adj2[kr]->point(adj2[kr]->tMin());
			  //     Point p2_2 = adj2[kr]->point(adj2[kr]->tMax());
			  //     double d1 = std::min(p1_1.dist(p2_1), p1_1.dist(p2_2)) +
			  // 	std::min(p1_2.dist(p2_1), p1_2.dist(p2_2));
			  //     for (size_t kh=kb+1; kh<pc1.size(); ++kh)
			  // 	{
			  // 	  if (pc2[kh].get() != cv2[kr].get())
			  // 	    continue;
			  // 	  Point p1_3 = e2[kh]->point(e2[kh]->tMin());
			  // 	  Point p1_4 = e2[kh]->point(e2[kh]->tMax());
			  // 	  double d2 = std::min(p1_3.dist(p2_1), p1_3.dist(p2_2)) +
			  // 	    std::min(p1_4.dist(p2_1), p1_4.dist(p2_2));
			  // 	  if (d2 < d1)
			  // 	    {
			  // 	      d1 = d2;
			  // 	      kb = kh;
			  // 	    }
			  // 	}
			  //     break;
			  //   }
			  // else
			    break;
			}
		    }
		  if (ka == pc1.size() || kb == pc2.size())
		    {
		      THROW("Unexpected edge configuration");
		    }

		  bool cv_before = false;
		  for (size_t kh=0; kh<kr; ++kh)
		    {
		      if (cv1[kh].get() == cv1[kr].get() ||
			  cv2[kh].get() == cv2[kr].get())
			cv_before = true;
		    }
		  if (e1[ka]->twin() || e2[kb]->twin())
		    {
		      // This might get mixed up if more than two edges from
		      // one face corresponds to one edge from the other
		      if (e1[ka]->twin() && e2[kb]->twin())
			{
			  // Must check if the joint point is equal. In
			  // that case, it is probably OK to leave the
			  // configuration as it is. Otherwise, there will
			  // be a lot of splitting and connecting and logic
			  // to connect the correct pieces.
			  // Currently, do nothing
			  int stop_break = 1;
			}
		      else if (e1[ka]->twin())
			{
			  ftEdge *other2 = e1[ka]->twin()->geomEdge();
			  Point pt1_1 = e1[ka]->point(e1[ka]->tMin());
			  Point pt1_2 = e1[ka]->point(e1[ka]->tMax());
			  Point pt2_1 = e2[kb]->point(e2[kb]->tMin());
			  Point pt2_2 = e2[kb]->point(e2[kb]->tMax());
			  Point pt2_3 = other2->point(other2->tMin());
			  Point pt2_4 = other2->point(other2->tMax());
			  double d1 = pt2_1.dist(pt2_3);
			  double d2 = pt2_1.dist(pt2_4);
			  double d3 = pt2_2.dist(pt2_3);
			  double d4 = pt2_2.dist(pt2_4);

			  Point mid;
			  bool first;
			  if (std::min(d1,d2) < std::min(d3,d4))
			    {
			      mid = pt2_1;
			      first = (pt2_2.dist(pt1_1) < pt2_2.dist(pt1_2));
			    }
			  else
			    {
			      mid = pt2_2;
			      first = (pt2_4.dist(pt1_1) >= pt2_4.dist(pt1_2));
			    }
			  double par, dist;
			  Point close;
			  e1[ka]->closestPoint(mid, par, close, dist);

			  if (par > e1[ka]->tMin()+DEFAULT_PARAMETER_EPSILON &&
			      par < e1[ka]->tMax()-DEFAULT_PARAMETER_EPSILON &&
			      pt1_1.dist(close) > eps && pt1_2.dist(close) > eps)
			    {
			      // Split edge and if necessary split also underlying curve to prepare
			      // for surface approximation
			      e1[ka]->disconnectTwin();
			      shared_ptr<ftEdge> other1;
			      if (e1.size() < 4)
				other1 = splitEdge(e1[ka], par);
			      else
				other1 = e1[ka]->split2(par);
			      int stat;
			      if (first)
				{
				  e1[ka]->connectTwin(e2[kb].get(), stat);
				  other1->connectTwin(other2, stat);
				}
			      else
				{
				  e1[ka]->connectTwin(other2, stat);
				  other1->connectTwin(e2[kb].get(), stat);
				}
			    }
			  int stop_break = 1;
			}
		      else
			{
			  ftEdge *other1 = e2[kb]->twin()->geomEdge();
			  Point pt1_1 = e1[ka]->point(e1[ka]->tMin());
			  Point pt1_2 = e1[ka]->point(e1[ka]->tMax());
			  Point pt1_3 = other1->point(other1->tMin());
			  Point pt1_4 = other1->point(other1->tMax());
			  Point pt2_1 = e2[kb]->point(e2[kb]->tMin());
			  Point pt2_2 = e2[kb]->point(e2[kb]->tMax());
			  double d1 = pt1_1.dist(pt1_3);
			  double d2 = pt1_1.dist(pt1_4);
			  double d3 = pt1_2.dist(pt1_3);
			  double d4 = pt1_2.dist(pt1_4);

			  Point mid;
			  bool first;
			  if (std::min(d1,d2) < std::min(d3,d4))
			    {
			      mid = pt1_1;
			      first = (pt1_2.dist(pt2_1) < pt1_4.dist(pt2_2));
			    }
			  else
			    {
			      mid = pt1_2;
			      first = (pt1_3.dist(pt2_1) >= pt1_3.dist(pt2_2));
			    }
			  double par, dist;
			  Point close;
			  e2[kb]->closestPoint(mid, par, close, dist);

			  if (par > e2[kb]->tMin()+DEFAULT_PARAMETER_EPSILON &&
			      par < e2[kb]->tMax()-DEFAULT_PARAMETER_EPSILON &&
			      pt2_1.dist(close) > eps && pt2_2.dist(close) > eps)
			    {
			      e2[kb]->disconnectTwin();
			      shared_ptr<ftEdge> other2;
			      if (e2.size() < 4)
				other2 = splitEdge(e2[kb], par);
			      else
				other2 = e2[kb]->split2(par);
			      int stat;
			      if (first)
				{
				  e1[ka]->connectTwin(e2[kb].get(), stat);
				  other1->connectTwin(other2.get(), stat);
				}
			      else
				{
				  e1[ka]->connectTwin(other2.get(), stat);
				  other1->connectTwin(e2[kb].get(), stat);
				}
			    }
			  int stop_break = 1;
			}
		    }
		  else 
		    {
		      int stat;
		      e1[ka]->connectTwin(e2[kb].get(), stat);
		    }
		}
	    }
	  // while (adjacent)
	  //   {
	  //     // Adjacent faces. Transfer adjacency to corresponding
	  //     // parameter volume faces
	  //     // Find corresponding edges in parameter volume
	  //     // First find parameter value of end vertices of common edge
	  //     Point pos1 = edge1->getVertex(true)->getVertexPoint();
	  //     Point pos2 = edge1->getVertex(false)->getVertexPoint();

	  //     double u1, u2, v1, v2, w1, w2, dd1, dd2, paru1, paru2, parv1, parv2;
	  //     Point clo_pt1, clo_pt2;
	  //     closestBoundaryPoint(pos1, u1, v1, w1, clo_pt1, dd1, paru1, 
	  // 			   parv1, toptol_.gap);
	  //     closestBoundaryPoint(pos2, u2, v2, w2, clo_pt2, dd2, paru2, 
	  // 			   parv2, toptol_.gap);
	  //     Point par1(u1, v1, w1);
	  //     Point par2(u2, v2, w2);

	  //     // Identify corresponding edges in parameter faces. First
	  //     // check endpoints
	  //     vector<shared_ptr<ftEdge> > e1 = face2[ki]->getAllEdges(0);
	  //     vector<shared_ptr<ftEdge> > e2 = face2[kj]->getAllEdges(0);

	  //     vector<double> edgdist1(e1.size());
	  //     vector<double> edgdist2(e2.size());
	  //     double mindist1 = std::numeric_limits<double>::max(), mindist2 = std::numeric_limits<double>::max();
	  //     int min_ix1 = -1, min_ix2 = -1;
	  //     for (size_t kr=0; kr<e1.size(); ++kr)
	  // 	{
	  // 	  Point vxpos1 = e1[kr]->point(e1[kr]->tMin());		
	  // 	  Point vxpos2 = e1[kr]->point(e1[kr]->tMax());
	  // 	  edgdist1[kr] = std::min(par1.dist(vxpos1), par1.dist(vxpos2)) +
	  // 	    std::min(par2.dist(vxpos1), par2.dist(vxpos2));
	  // 	  if (edgdist1[kr] < mindist1)
	  // 	    {
	  // 	      mindist1 = edgdist1[kr];
	  // 	      min_ix1 = (int)kr;
	  // 	    }
	  // 	}

	  //     for (size_t kr=0; kr<e2.size(); ++kr)
	  // 	{
	  // 	  Point vxpos1 = e2[kr]->point(e2[kr]->tMin());		
	  // 	  Point vxpos2 = e2[kr]->point(e2[kr]->tMax());
	  // 	  edgdist2[kr] = std::min(par1.dist(vxpos1), par1.dist(vxpos2)) +
	  // 	    std::min(par2.dist(vxpos1), par2.dist(vxpos2));
	  // 	  if (edgdist2[kr] < mindist2)
	  // 	    {
	  // 	      mindist2 = edgdist2[kr];
	  // 	      min_ix2 = (int)kr;
	  // 	    }
	  // 	}
	  //     int stat;
	  //     double len1 = e1[min_ix1]->point(e1[min_ix1]->tMin()).dist(e1[min_ix1]->point(e1[min_ix1]->tMax()));
	  //     double len2 = e2[min_ix2]->point(e2[min_ix2]->tMin()).dist(e2[min_ix2]->point(e2[min_ix2]->tMax()));
	  //     if (mindist1 < toptol_.neighbour && mindist2 < toptol_.neighbour)
	  // 	e1[min_ix1]->connectTwin(e2[min_ix2].get(), stat);
	  //     else
	  // 	{
	  // 	  // Not exact corner match. Perform closer check
	  // 	  double mindist = std::numeric_limits<double>::max();
	  // 	  min_ix1 = min_ix2 = -1;
	  // 	  for (size_t kr=0; kr<e1.size(); ++kr)
	  // 	    {
	  // 	      Point vxpos1 = e1[kr]->point(e1[kr]->tMin());		
	  // 	      Point vxpos2 = e1[kr]->point(e1[kr]->tMax());
	  // 	      for (size_t kh=0; kh<e2.size(); ++kh)
	  // 		{
	  // 		  Point vxpos3 = e2[kh]->point(e2[kh]->tMin());		
	  // 		  Point vxpos4 = e2[kh]->point(e2[kh]->tMax());

	  // 		  double t1, t2, t3, t4, d1, d2, d3, d4;
	  // 		  Point p1, p2, p3, p4;
	  // 		  e1[kr]->closestPoint(vxpos3, t1, p1, d1);
	  // 		  e1[kr]->closestPoint(vxpos4, t2, p2, d2);
	  // 		  e2[kh]->closestPoint(vxpos1, t3, p3, d3);
	  // 		  e2[kh]->closestPoint(vxpos2, t4, p4, d4);
	  // 		  double dist = min(min(min(d1+d2, d1+d4),
	  // 					min(d2+d3, d3+d4)),
	  // 				    min(d1+d3,d2+d4));
	  // 		  if (dist < mindist)
	  // 		    {
	  // 		      mindist = dist;
	  // 		      min_ix1 = (int)kr;
	  // 		      min_ix2 = (int)kh;
	  // 		    }
	  // 		}
	  // 	    }
	  // 	  e1[min_ix1]->connectTwin(e2[min_ix2].get(), stat);
	  // 	}

	  //     // Check for more adjacency between the current faces
	  //     adj_ix++;
	  //     adjacent = face1[ki]->areNeighbours(face1[kj].get(), edge1,
	  // 					  edge2, adj_ix); 
	  //   }
	}

    }
}

//===========================================================================
// 
shared_ptr<ftEdge> ftVolume::splitEdge(shared_ptr<ftEdge> edge, 
				       double par) const
// 
//===========================================================================
{
  // Split edge
  shared_ptr<ParamCurve> crv = edge->geomCurve();
  shared_ptr<ftEdge> other = edge->split2(par);

  // Split corresponding curve and update curve loop. First fetch surface
  shared_ptr<ParamSurface> surf = edge->face()->asFtSurface()->surface();
  shared_ptr<BoundedSurface> bd_surf = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
  if (bd_surf.get())
    {
      // Replace curve to split by the new pieces
      int nmb_loops = bd_surf->numberOfLoops();
      for (int ki=0; ki<nmb_loops; ++ki)
	{
	  shared_ptr<CurveLoop> loop = bd_surf->loop(ki);
	  int nmb_cvs = loop->size();
	  for (int kj=0; kj<nmb_cvs; ++kj)
	    {
	      shared_ptr<ParamCurve> curr = (*loop)[kj];
	      if (curr.get() == crv.get())
		{
		  vector<shared_ptr<ParamCurve> > sub_cvs = loop->split(kj, par);
		  if (sub_cvs.size() > 0)
		    edge->setGeomCurve(sub_cvs[0]);
		  if (sub_cvs.size() > 1)
		    other->setGeomCurve(sub_cvs[1]);
		}
	    }
	}
    }
  return other;
}

//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> > 
ftVolume::replaceWithRegVolumes(int degree,
				vector<SurfaceModel*>& modified_adjacent,
				bool perform_step2,
				int split_mode,
				bool pattern_split,
				bool accept_degen,
				int level)
//===========================================================================
{
  vector<shared_ptr<ftVolume> > reg_vols;

  // Make sure that all anticipated edges are in place
  vector<pair<Point,Point> > corr_vx_pts;
  (void)regularizeBdShells(corr_vx_pts, modified_adjacent,
			   split_mode, pattern_split, level);

#ifdef DEBUG_VOL1
  std::ofstream of("regvol1.g2");
  vector<shared_ptr<SurfaceModel> > shells = getAllShells();
  vector<shared_ptr<Vertex> > vx;
  for (size_t ki=0; ki<shells.size(); ++ki)
    {
      int nmb = shells[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> tmp = shells[ki]->getSurface(kj);
	  tmp->writeStandardHeader(of);
	  tmp->write(of);
	}
      vector<shared_ptr<Vertex> > tmp_vx;
      shells[ki]->getAllVertices(tmp_vx);
      vx.insert(vx.end(), tmp_vx.begin(), tmp_vx.end());
    }
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << vx.size() << std::endl;
  for (size_t ki=0; ki<vx.size(); ++ki)
    of << vx[ki]->getVertexPoint() << std::endl;
#endif

#ifdef DEBUG_VOL1
  // Check
  int nmb = shells_[0]->nmbBoundaries();
  if (nmb > 0)
    std::cout << "Post regularize model: number of boundaries = " << nmb << std::endl;
#endif

  // Generate missing boundary surfaces
  int nmb_faces = nmbOfFaces();  // Current number of faces in all boundary shells

  bool trimmed;
  vector<shared_ptr<ftSurface> > bd_faces = 
    generateMissingBdSurf(degree, corr_vx_pts, perform_step2, false, trimmed);
  bool first = true;
  while (bd_faces.size() == 0 && (first || nmb_faces != nmbOfFaces()))
    {
      first = false;

      // The boundary shells are updated. Try again to generate missing
      // boundary surfaces
      nmb_faces = nmbOfFaces();
      //bd_faces = generateMissingBdSurf(corr_vx_pts, perform_step2);
      bd_faces = generateMissingBdSurf(degree, corr_vx_pts, false, true, 
				       trimmed);
    }

  if (bd_faces.size() == 0)
    return reg_vols;  // Dummy array

  // Make regular, trimmed volumes
  reg_vols = createRegularVolumes(bd_faces);

  // Check if the regularization is completed
  if (reg_vols.size() > 1)
    {
      if (!trimmed)
	++level;   // New try to block structure the current volume
      size_t nmb = reg_vols.size();
      for (size_t kj=0; kj<nmb; ++kj)
	{
	  bool reg = reg_vols[kj]->isRegularized(accept_degen);
	  if (!reg)
	    {
#ifdef DEBUG_VOL1
	      std::ofstream of6_2("notreg_vol.g2");
	      shared_ptr<SurfaceModel> mod =  reg_vols[kj]->getOuterShell();
	      int nmb_vol = mod->nmbEntities();
	      for (int kr=0; kr<nmb_vol; ++kr)
		{
		  shared_ptr<ParamSurface> sf = mod->getSurface(kr);
		  sf->writeStandardHeader(of6_2);
		  sf->write(of6_2);
		}
#endif

	      vector<shared_ptr<ftVolume> > reg_vols2 = 
		reg_vols[kj]->replaceWithRegVolumes(degree, modified_adjacent,
						    true, split_mode,
						    pattern_split, accept_degen,
						    level);
	      if (reg_vols2.size() > 1)
		{
		  reg_vols.erase(reg_vols.begin()+kj);
		  reg_vols.insert(reg_vols.end(), reg_vols2.begin(), 
				  reg_vols2.end());
		  nmb--;
		  kj--;
		}
	    }
	}
    }
  return reg_vols;
}

//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> > 
ftVolume::splitConcaveVol(int degree, bool isolate)
//===========================================================================
{
  vector<shared_ptr<ftVolume> > reg_vols;
  if (nmbOfShells() != 1)
    return reg_vols;

#ifdef DEBUG_VOL1
  std::ofstream of0("vol0.g2");
  vector<shared_ptr<SurfaceModel> > shells = getAllShells();
  for (size_t ki=0; ki<shells.size(); ++ki)
    {
      int nmb = shells[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> tmp = shells[ki]->getSurface(kj);
	  tmp->writeStandardHeader(of0);
	  tmp->write(of0);
	}
    }
#endif

  // Make sure that all anticipated edges are in place
  shared_ptr<SurfaceModel> sfmodel = getShell(0);
  ModifyFaceSet splitmod(sfmodel);

  // Check if a suitable splitting surfaces exist
  vector<shared_ptr<ParamSurface> > split_sfs;
  vector<vector<ftEdge*> > split_edgs;
  vector<ftSurface*> corr_faces;
  splitmod.getSplittingSurface(split_sfs, corr_faces, split_edgs);
#ifdef DEBUG_VOL1
  std::cout << "Nmb split surfaces: " << split_sfs.size() << std::endl;
#endif
  if (split_sfs.size() > 0)
    {
      double bend = sfmodel->getTolerances().bend;
      int create_models = (isolate) ? 4 : 3;  // Create both sides, 

     //  // To avoid intersection problems with tangential connections to the
     //  // splitting surface, add such edges to the already known intersections
     //  // This involves a risk of missing intersections, but it is believed
     //  // to be smaller than the risk of intersection problems
     //  // connect models if requested
     //  FaceConnectivityUtils<ftEdgeBase,ftSurface> connectivity;
     //  vector<ftEdgeBase*> vec0;
     //  vector<ftSurface*> faces0;
     //  faces0.push_back(corr_faces[0]); 
     //  connectivity.smoothEdges(faces0, vec0, bend);
     // for (size_t kr=0; kr<vec0.size(); ++kr)
     // 	split_edgs[0].push_back(vec0[kr]->geomEdge());
      
      try {
	vector<shared_ptr<ftVolume> > curr_vols =
	  ftVolumeTools::splitWithSplitSf(this, split_sfs[0], split_edgs[0],
					  toptol_.gap, create_models);
	if (curr_vols.size() > 1)
	  reg_vols.insert(reg_vols.end(), curr_vols.begin(), curr_vols.end());
      }
      catch (...)
	{
	  MESSAGE("Splitting surface division failed");
	}

      for (size_t ki=1; ki<split_sfs.size(); ++ki)
	{
	  size_t nmb_reg = reg_vols.size();
	  for (size_t kj=0; kj<nmb_reg; )
	    {
	      // // Add tangential connections to the list of already known
	      // // intersections
	      // vector<ftEdgeBase*> vec1;
	      // vector<ftSurface*> faces1;
	      // faces1.push_back(corr_faces[1]);
	      // connectivity.smoothEdges(faces1, vec1, bend);
	      // for (size_t kr=0; kr<vec1.size(); ++kr)
	      // 	split_edgs[kj].push_back(vec1[kr]->geomEdge());
 	      try {
		vector<shared_ptr<ftVolume> > curr_vols =
		  ftVolumeTools::splitWithSplitSf(reg_vols[kj].get(), 
						  split_sfs[ki], split_edgs[ki], 
						  toptol_.gap, create_models);
		if (curr_vols.size() > 1)
		  {
		    reg_vols.insert(reg_vols.end(), curr_vols.begin(), 
				    curr_vols.end());
		    reg_vols.erase(reg_vols.begin()+kj);
		    --nmb_reg;
		  }
		else
		  ++kj;
	      }
	      catch (...)
		{
		  MESSAGE("Splitting surface division failed");
		}
	    }
	}
      if (reg_vols.size() > 1)
	{
	  return reg_vols;
	}
    }

  int nmb_div = 0;
  shared_ptr<SurfaceModel> sfmodel2 = splitmod.getModifiedModel(nmb_div);
#ifdef DEBUG_VOL1
  std::cout << "Nmb div: " << nmb_div << std::endl;
  std::ofstream of1("vol1.g2");
  vector<shared_ptr<SurfaceModel> > shells2 = getAllShells();
  for (size_t ki=0; ki<shells2.size(); ++ki)
    {
      int nmb = shells2[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> tmp = shells2[ki]->getSurface(kj);
	  tmp->writeStandardHeader(of1);
	  tmp->write(of1);
	}
    }
#endif

  if (nmb_div < 4)
    return reg_vols;  // No split

  bool trimmed;
  vector<pair<Point,Point> > corr_vx_pts;
  vector<shared_ptr<ftSurface> > bd_faces = 
    generateMissingBdSurf(degree, corr_vx_pts, false, false, trimmed, nmb_div);
 
  reg_vols = createRegularVolumes(bd_faces);
  if (isolate)
    {
      for (size_t ki=0; ki<bd_faces.size(); ++ki)
	if (bd_faces[ki]->hasTwin())
	  bd_faces[ki]->disconnectTwin(true);

      // Simplify sub models if possible
      for (size_t ki=0; ki<reg_vols.size(); ++ki)
	{
	  RegularizeFaceSet regularize(reg_vols[ki]->getShell(0));
	  regularize.removeExtraDiv(true);
	}
    }
#ifdef DEBUG_VOL1
  for (size_t ki=0; ki<reg_vols.size(); ++ki)
    {
      reg_vols[ki]->getOuterShell()->checkShellTopology();
    }
#endif

  return reg_vols;
}

//===========================================================================
// 
// 
bool 
ftVolume::sortRegularSurfaces(vector<shared_ptr<ParamSurface> >& sorted_sfs,
			      vector<std::pair<int,double> >& classification,
			      vector<int>& deg_type)
//===========================================================================
{
  // @@@ VSK. 0217. This function is very messy, in particular for degenerate
  // cases. Must be cleaned up when we have found solutions for most occuring
  // situations.

#ifdef DEBUG_VOL1
  std::ofstream ofin("sort_in.g2");
  shared_ptr<SurfaceModel> shell = getOuterShell();
  int nmb_sf = shell->nmbEntities();
  for (int kj1=0; kj1<nmb_sf; ++kj1)
    {
      shared_ptr<ParamSurface> tmp = shell->getSurface(kj1);
      tmp->writeStandardHeader(ofin);
      tmp->write(ofin);
    }
 #endif

  std::fill(classification.begin(), classification.end(), 
	    make_pair(0, -1.0));
  std::fill(deg_type.begin(), deg_type.end(), -1); // Not yet defined
  // 0 = not degenerate, 1 = degenerate to triangle, 2 = degenerate to line, 
  // 3 = degenerate to point

  // Fetch surfaces
  vector<shared_ptr<ParamSurface> > sfs;
  vector<int> sfs_bd;
  int nmb = shells_[0]->nmbEntities();
  int nmb_sorted = 6;
  if (nmb > nmb_sorted)
    return false;  // Already tested. Just for completeness

  vector<pair<shared_ptr<ParamSurface>, Point> > deg_pts;
  int ki;
  for (ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamSurface> curr = shells_[0]->getSurface(ki);

      // Analyze also twin surface, if any
      shared_ptr<ParamSurface> curr2;
      if (shells_[0]->getFace(ki)->twin())
	{
	  ftSurface* face2 = shells_[0]->getFace(ki)->twin()->asFtSurface();
	  if (face2)
	    curr2 = face2->surface();
	}

      // Check if the current surface is a degenerate, non-trimmed surface
      // In that case, identify the degenerate corner
      int idx = shells_[0]->getIndex(curr.get());
      shared_ptr<ftSurface> face = shells_[0]->getFace(idx);
      int nmb_bd = face->nmbOuterBdCrvs(toptol_.gap, toptol_.neighbour,
					toptol_.bend);	

      // Check also with the existence of adjacent faces
      vector<ftSurface*> adj_faces;
      face->getAdjacentFaces(adj_faces);
      if (nmb_bd >= 4 && adj_faces.size() < 4)
	nmb_bd = (int)adj_faces.size();
      
      bool no_trim = false;
      if (curr->instanceType() != Class_BoundedSurface)
	no_trim = true;
      else
	{
	  shared_ptr<BoundedSurface> bd_curr = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(curr);
	  if (bd_curr->isIsoTrimmed(toptol_.gap))
	    no_trim = true;
	}
      if (nmb_bd == 3 && no_trim)
	{
	  bool bottom, right, top, left;
	  bool deg = curr->isDegenerate(bottom, right, top, left, toptol_.gap);
	  if (!deg)
	    deg = curr->isDegenerate(bottom, right, top, left, toptol_.neighbour);
	  RectDomain dom = curr->containingDomain();
	  double upar = left ? dom.umin() : 
	    (right ? dom.umax() : 0.5*(dom.umin()+dom.umax()));
	  double vpar = bottom ? dom.vmin() : 
	    (top ? dom.vmax() : 0.5*(dom.vmin()+dom.vmax()));
	  Point degpt = curr->point(upar, vpar);
	  deg_pts.push_back(make_pair(curr, degpt));
	}
      else if (nmb_bd == 3 && curr2.get() && 
	       curr2->instanceType() != Class_BoundedSurface)
	{
	  bool bottom, right, top, left;
	  bool deg = curr2->isDegenerate(bottom, right, top, left, toptol_.gap);
	  if (!deg)
	    deg = curr2->isDegenerate(bottom, right, top, left, toptol_.neighbour);
	  RectDomain dom = curr2->containingDomain();
	  double upar = left ? dom.umin() : 
	    (right ? dom.umax() : 0.5*(dom.umin()+dom.umax()));
	  double vpar = bottom ? dom.vmin() : 
	    (top ? dom.vmax() : 0.5*(dom.vmin()+dom.vmax()));
	  Point degpt = curr2->point(upar, vpar);
	  deg_pts.push_back(make_pair(curr, degpt));
	}

      // Check if any surface knows which boundary it belongs to
      shared_ptr<SurfaceOnVolume> vol_sf = getVolSf(curr);
      int orientation;
      bool swap;
      int bd = -1;
      if (vol_sf.get())
	{
	  bd = vol_sf->whichBoundary(0.0, orientation, swap);

	  if (vol_sf->getVolume().get() != vol_.get())
	    bd = -1;
	}
 
      if (bd >= 0 && sorted_sfs[bd].get())
	bd = -1;  // More than one surface found for the same boundary

      // Check for feasability
      if (bd >= 0)
	{
	  int ka;
	  int ix1 = shells_[0]->getIndex(curr.get());
	  for (ka=0; ka<nmb_sorted; ++ka)
	    {
	      if (ka/2 == bd/2)
		continue;   // Same parameter direction
	      if (!sorted_sfs[ka].get())
		continue;  // No surface found
	      int ix2 = shells_[0]->getIndex(sorted_sfs[ka].get());
	      shared_ptr<ftSurface> f1 = shells_[0]->getFace(ix1);
	      shared_ptr<ftSurface> f2 = shells_[0]->getFace(ix2);
	      bool smooth;
	      if (!f1->isAdjacent(f2.get(), smooth))
		break;    
	    }
	  if (ka < nmb_sorted)
	    bd = -1;
	}

      if (bd >=0)
	{
	  sorted_sfs[bd] = curr;

	  double par = vol_sf->getConstVal();
	  classification[bd] = make_pair(2, par);
	}
      else
	{
	  sfs.push_back(curr);
	  sfs_bd.push_back(nmb_bd);
	}

    }

  // Sort the remaining surfaces
  // First look for iso parameter information
  for (ki=0; ki<nmb_sorted; ++ki)
    {
      if (!sorted_sfs[ki].get())
	{
	  for (size_t kj=0; kj<sfs.size(); ++kj)
	    {
	      shared_ptr<SurfaceOnVolume> vol_sf = getVolSf(sfs[kj]);
	      if (vol_sf.get() && vol_sf->getVolume().get() == vol_.get())
		{
		  int dir = vol_sf->getConstDir();
		  if (dir == ki/2 + 1)
		    {
		      // Check for feasiblity
		      int ka;
		      int ix1 = shells_[0]->getIndex(sfs[kj].get());
		      for (ka=0; ka<nmb_sorted; ++ka)
			{
			  if (ka/2 == dir/2)
			    continue;   // Same parameter direction
			  if (!sorted_sfs[ka].get())
			    continue;  // No surface found
			  int ix2 = shells_[0]->getIndex(sorted_sfs[ka].get());
			  shared_ptr<ftSurface> f1 = shells_[0]->getFace(ix1);
			  shared_ptr<ftSurface> f2 = shells_[0]->getFace(ix2);
			  bool smooth;
			  if (!f1->isAdjacent(f2.get(), smooth))
			    break;    
			}
		      if (ka == nmb_sorted)
			{
			  double par = vol_sf->getConstVal();
			  sorted_sfs[ki] = sfs[kj];
			  classification[ki] = make_pair(1, par);
			  sfs.erase(sfs.begin()+kj);
			  sfs_bd.erase(sfs_bd.begin()+kj);
			  break;
			}
		    }
		}
	    }
	}
    }

#ifdef DEBUG_VOL1
  std::ofstream of0("midsurf0.g2");
  for (size_t kr=0; kr<sfs.size(); ++kr)
    {
      sfs[kr]->writeStandardHeader(of0);
      sfs[kr]->write(of0);
    }
#endif


  // Define approximate parameter domain by projecting the vertices of 
  // this body onto the volume
  vector<shared_ptr<Vertex> > vxs =  vertices(); 
  double parval[6];
#ifdef DEBUG_VOL1
  std::ofstream ofvx("midsurf_vx.g2");
  ofvx << "400 1 0 4 255 0 0 255" << std::endl;
  ofvx << vxs.size() << std::endl;
#endif
  for (size_t kr=0; kr<vxs.size(); ++kr)
    {
#ifdef DEBUG_VOL1
      ofvx << vxs[kr]->getVertexPoint() << std::endl;
#endif

      double upar, vpar, wpar, dist;
      Point pnt;
      vol_->closestPoint(vxs[kr]->getVertexPoint(), upar, vpar, wpar,
			 pnt, dist, toptol_.gap);
      if (kr == 0)
	{
	  parval[0] = parval[1] = upar;
	  parval[2] = parval[3] = vpar;
	  parval[4] = parval[5] = wpar;
	}
      else
	{
	  parval[0] = std::min(parval[0], upar);
	  parval[1] = std::max(parval[1], upar);
	  parval[2] = std::min(parval[2], vpar);
	  parval[3] = std::max(parval[3], vpar);
	  parval[4] = std::min(parval[4], wpar);
	  parval[5] = std::max(parval[5], wpar);
	}
    }

  // Fetch the opposite surface when one surface in a parameter direction
  // is given
  size_t kr;
  for (ki=0; ki<nmb_sorted; ki+=2)
    {
      if ((sorted_sfs[ki].get() && !sorted_sfs[ki+1].get()) ||
	  (!sorted_sfs[ki].get() && sorted_sfs[ki+1].get()))
	{
	  int ix = (sorted_sfs[ki].get()) ? ki : ki+1;
	  int ix2 = (sorted_sfs[ki].get()) ? ki+1 : ki;
	  int idx = shells_[0]->getIndex(sorted_sfs[ix].get());
	  shared_ptr<ftSurface> face0 = shells_[0]->getFace(idx);

	  // Select the available surface not being a neighbour
	  vector<ftSurface*> neighbours;
	  face0->getAdjacentFaces(neighbours);
	  for (kr=0; kr<sfs.size(); ++kr)
	    {
	      int idx2 = shells_[0]->getIndex(sfs[kr].get());
	      shared_ptr<ftSurface> face2 = shells_[0]->getFace(idx2);
	      shared_ptr<ftEdge> e1, e2;
	      bool adjacent = face0->areNeighbours(face2.get(), e1, e2);
	      if (!adjacent)
		{
		  // Check corner adjacency
		  vector<shared_ptr<Vertex> > common_vx = 
		    face0->getCommonVertices(face2.get());
		  if (common_vx.size() > 0)
		    adjacent = true;
		}
	      if (!adjacent)
		break;
	    }
	  if (kr < sfs.size())
	    {
	      sorted_sfs[ix2] = sfs[kr];
	      classification[ix2] = make_pair(0, parval[ix2]);
	      sfs.erase(sfs.begin()+kr);
	      sfs_bd.erase(sfs_bd.begin()+kr);
	    }
	}
    }

  if (nmb == nmb_sorted && sfs.size() == 1)
    {
      // Only one free spot and one surface left
      for (int ki=0; ki<(int)sorted_sfs.size(); ++ki)
	{
	  if (!sorted_sfs[ki].get())
	    {
	      sorted_sfs[ki] = sfs[0];
	      classification[ki] = make_pair(0, parval[ki]);
	      sfs.erase(sfs.begin());
	      sfs_bd.erase(sfs_bd.begin());
	    }
	}
    }
      
  
  // Degeneracy classification of side surfaces.
  // Start with sorted surfaces
  for (int ki=0; ki<(int)sorted_sfs.size(); ++ki)
    {
      if (!sorted_sfs[ki].get())
	continue;
      int idx = shells_[0]->getIndex(sorted_sfs[ki].get());
      shared_ptr<ftSurface> face = shells_[0]->getFace(idx);
      int nmb_bd = face->nmbOuterBdCrvs(toptol_.gap, toptol_.neighbour,
					toptol_.bend);	
      // Check also with the existence of adjacent faces
      vector<ftSurface*> adj_faces;
      face->getAdjacentFaces(adj_faces);
      if (nmb_bd >= 4 && adj_faces.size() < 4)
	nmb_bd = (int)adj_faces.size();
      
      if (nmb_bd >= 4)
	deg_type[ki] = 0;
      else if (nmb_bd == 3)
	deg_type[ki] = 1;
    }

  // Indices into sorted_sfs for each parameter direction
  int par_sfs[3][4] = {{0, 2, 1, 3}, {1, 5, 0, 4}, {2, 4, 3, 5}};

  if (sfs.size() > 1)
    {
      // Sort remaining surfaces to start with degenerate ones
      int idx1 = shells_[0]->getIndex(sfs[0].get());
      int nmb_bd1 = 
	shells_[0]->getFace(idx1)->nmbOuterBdCrvs(toptol_.gap, toptol_.neighbour,
						 toptol_.bend);
      // Check also with the existence of adjacent faces
      vector<ftSurface*> adj_faces;
      shells_[0]->getFace(idx1)->getAdjacentFaces(adj_faces);
      if (nmb_bd1 >= 4 && adj_faces.size() < 4)
	nmb_bd1 = (int)adj_faces.size();

      for (ki=1; ki<(int)sfs.size(); ++ki)
	{
	  int idx2 = shells_[0]->getIndex(sfs[ki].get());
	  int nmb_bd2 = 
	    shells_[0]->getFace(idx2)->nmbOuterBdCrvs(toptol_.gap, toptol_.neighbour,
						 toptol_.bend);
	  // Check also with the existence of adjacent faces
	  vector<ftSurface*> adj_faces2;
	  shells_[0]->getFace(idx2)->getAdjacentFaces(adj_faces2);
	  if (nmb_bd2 >= 4 && adj_faces2.size() < 4)
	    nmb_bd2 = (int)adj_faces2.size();

	  if (nmb_bd2 < nmb_bd1 || (nmb_bd2 == nmb_bd1 && nmb_bd1 < 4))
	    {
	      sfs.insert(sfs.begin(), sfs[ki]);
	      sfs_bd.insert(sfs_bd.begin(), sfs_bd[ki]);
	      sfs.erase(sfs.begin()+ki+1);
	      sfs_bd.erase(sfs_bd.begin()+ki+1);
	    }
	}
    }

  if (deg_pts.size() > 0 && sfs.size() > 0 && nmb < nmb_sorted)
    {
      // Ensure consistence between the empty surface slot and the degenerate point(s)
      // Identify parameter direction for the empty slot
      int par;
      int nmb_sfs;
      for (par=0; par<3; ++par)
	{
	  nmb_sfs = 0;
	  for (ki=0; ki<4; ++ki)
	    if (sorted_sfs[par_sfs[par][ki]].get())
	      nmb_sfs++;
	  if ((int)sfs.size() < 4-nmb_sfs)
	    break;
	}

      for (ki=0; ki<(int)deg_pts.size(); ++ki)
	{
	  // Only applicable if the degenerate surface does not belong to
	  // this parameter direction
	  int kj;
	  for (kj=0; kj<4; ++kj)
	    if (par < 3 &&
		sorted_sfs[par_sfs[par][kj]].get() == deg_pts[ki].first.get())
	      break;
	  if (kj < 4)
	    break;
	  for (kj=0; kj<(int)sfs.size(); ++kj)
	    if (sfs[kj].get() == deg_pts[ki].first.get())
	      break;
	  if (kj <(int) sfs.size())
	    break;

	  // Compute distances between the degenerate point and the surfaces
	  // belonging to this parameter direction. The empty slot must lie
	  // between two adjacent surfaces touching the degenerate point and
	  // the opposite surface must be distant from the point
	  vector<double> dist1(nmb_sfs);
	  vector<double> dist2(sfs.size());
	  int kh=0;
	  int nmb_touch1 = 0, nmb_touch2 = 0;
	  for (kj=0; kj<4; ++kj)
	    {
	      if (par < 3 && sorted_sfs[par_sfs[par][kj]].get())
		{
		  double upar, vpar, dd;
		  Point clo_pt;
		  sorted_sfs[par_sfs[par][kj]]->closestPoint(deg_pts[ki].second,
							     upar, vpar, clo_pt,
							     dd, toptol_.gap);
		  dist1[kh++] = dd;
		  if (dd < toptol_.gap)
		    nmb_touch1++;
		}
	    }
	  for (kj=0; kj<(int)sfs.size(); ++kj)
	    {
	      double upar, vpar, dd;
	      Point clo_pt;
	      sfs[kj]->closestPoint(deg_pts[ki].second, upar, vpar, clo_pt,
				    dd, toptol_.gap);
	      dist2[kj] = dd;
	      if (dd < toptol_.gap)
		nmb_touch2++;
	    }
	  
	  // Fallback in case of gaps
	  if (nmb_touch1 == 0)
	    {
	      for (kj=0; kj<(int)dist1.size(); ++kj)
		if (dist1[kj] < toptol_.neighbour)
		  ++nmb_touch1;
	    }
	  if (nmb_touch2 == 0)
	    {
	      for (kj=0; kj<(int)dist2.size(); ++kj)
		if (dist2[kj] < toptol_.neighbour)
		  ++nmb_touch2;
	    }
	  if (nmb_touch1+nmb_touch2 != 2)
	    {
	      MESSAGE("ftVolume::sortRegularSurfaces. Unexpected degeneracy configuration");
	      break;
	    }
	  if (par < 3 && nmb_sfs == 1 && nmb_touch1 == 1)
	    {
	      // Three empty slots to define. Identify the occupied one
	      for (kj=0; kj<4; ++kj)
		if (sorted_sfs[par_sfs[par][kj]].get())
		  break;
	      
	      // Two unsorted surfaces. Distinguish between the one touching the
	      // degenerate point and the other
	      int deg = (dist2[0] < dist2[1] /*toptol_.gap*/) ? 0 : 1;
	      int other = 1 - deg;

	      int ix = par_sfs[par][(kj+1)%4];
	      classification[ix] = make_pair(-1, parval[ix]);
	      deg_type[ix] = 2;

	      ix = par_sfs[par][(kj+2)%4];
	      sorted_sfs[ix] = sfs[deg];
	      classification[ix] = make_pair(0, parval[ix]);
	      deg_type[ix] = 0;

	      ix = par_sfs[par][(kj+3)%4];
	      sorted_sfs[ix] = sfs[other];
	      classification[ix] = make_pair(0, parval[ix]);
	      deg_type[ix] = 0;

	      sfs.erase(sfs.begin()+std::max(deg,other));
	      sfs_bd.erase(sfs_bd.begin()+std::max(deg,other));
	      sfs.erase(sfs.begin()+std::min(deg,other));
	      sfs_bd.erase(sfs_bd.begin()+std::min(deg,other));
	    }
	  else if (par < 3 && nmb_sfs == 2 && nmb_touch1 == 1)
	    {
	      // Two empty slots to define. Identify the occupied ones
	      int i1=-1, i2=-1;
	      for (kj=0; kj<4; ++kj)
		if (sorted_sfs[par_sfs[par][kj]].get())
		  {
		    if (i1 < 0)
		      i1 = kj;
		    else
		      i2 = kj;
		  }

	      // Expects the two sorted surfaces to be adjacent. Check.
	      int d1 = ((i1==0 && i2==3) || (i1==3 && i2==0)) ? 1 : abs(i1-i2);
	      if (d1 == 1)
		{
		  int ixdeg, ixother;
		  if (dist1[0] < dist1[1] /*toptol_.gap*/)
		    {
		      ixdeg = (i1 > 0) ? par_sfs[par][i1-1] : par_sfs[par][3];
		      ixother = (i1 > 1) ? par_sfs[par][i1-2] : par_sfs[par][2+i1];
		    }
		  else
		    {
		      ixdeg = par_sfs[par][(i2+1)%4];
		      ixother = par_sfs[par][(i2+2)%4];
		    }
		  classification[ixdeg] = make_pair(-1, parval[ixdeg]);
		  deg_type[ixdeg] = 2;

		  sorted_sfs[ixother] = sfs[0];
		  classification[ixother] = make_pair(0, parval[ixother]);
		  deg_type[ixother] = 0;
		  sfs.erase(sfs.begin());
		  sfs_bd.erase(sfs_bd.begin());
		}
	      else
		{
		  MESSAGE("ftVolume::sortRegularSurfaces. Unexpected degeneracy configuration");
		}
	    }
	  else
	    {
	      MESSAGE("ftVolume::sortRegularSurfaces. Degeneracy configuration not implemented");
	      break;
	    }
	  int stop_break = 1;
	  if (sfs.size() == 0)
	    break;
	}
    }

  // Sort the remaining surfaces
  for (ki=0; ki<(int)sorted_sfs.size(); ki+=2)
    {
      if (!sorted_sfs[ki].get() && !sorted_sfs[ki+1].get())
	{
	  size_t ka;
	  for (ka=0; ka<sfs.size(); ++ka)
	    {
	      // Select an arbitrary free surface
	      int idx = shells_[0]->getIndex(sfs[ka].get());
	      shared_ptr<ftSurface> face0 = shells_[0]->getFace(idx);

	      // Select the available surface not being a neighbour
	      size_t kj;
	      vector<ftSurface*> neighbours;
	      face0->getAdjacentFaces(neighbours);
	      for (kr=0; kr<sfs.size(); ++kr)
		{
		  if (kr == ka)
		    continue;
		  for (kj=0; kj<neighbours.size(); ++kj)
		    {
		      shared_ptr<ParamSurface> surf1 = neighbours[kj]->surface();
		  
		      if (surf1.get() == sfs[kr].get())
			break;
		    }
		  if (kj == neighbours.size())
		    break;
		}
	      if (kr < sfs.size())
		{
		  sorted_sfs[ki] = sfs[ka];
		  classification[ki] = make_pair(0, parval[ki]);
		  sorted_sfs[ki+1] = sfs[kr];
		  classification[ki+1] = make_pair(0, parval[ki+1]);
		  sfs.erase(sfs.begin()+ka);
		  sfs_bd.erase(sfs_bd.begin()+ka);
		  if (kr > ka)
		    --kr;
		  sfs.erase(sfs.begin()+kr);
		  sfs_bd.erase(sfs_bd.begin()+kr);
		  break;
		}
	    }
	}
    }

  // @@@ VSK. 0216. Must handle the different cases when they arise. Doesn't have
  // the complete overview
  if (nmb == nmb_sorted-1)
    {
      // In a degenerate situation without any prescribed degenerate points and
      // open surface pairs, it might be necessary to define on surface in the pair
      if (deg_pts.size() == 0)
	{
	  for (ki=0; ki<(int)sorted_sfs.size(); ki+=2)
	    {
	      if (!sorted_sfs[ki].get() && !sorted_sfs[ki+1].get())
		{
		  sorted_sfs[ki] = sfs[0];
		  classification[ki] = make_pair(0, parval[ki]);
		  sfs.erase(sfs.begin());
		  sfs_bd.erase(sfs_bd.begin());
		  break;
		}
	    }
	}

      // One surface degenerate to a line is expected
      if (sfs.size() == 0)
	{
	  for (int ki=0; ki<(int)sorted_sfs.size(); ++ki)
	    {
	      if (!sorted_sfs[ki].get())
		{
		  classification[ki] = make_pair(-1, parval[ki]);
		  deg_type[ki] = 2;
		}
	    }
	}
      else if (sfs.size() == 1)
	{
	  int idx = shells_[0]->getIndex(sfs[0].get());
	  shared_ptr<ftSurface> face = shells_[0]->getFace(idx);
	  int nmb_bd = face->nmbOuterBdCrvs(toptol_.gap, toptol_.neighbour,
					    toptol_.bend);
	  // Check also with the existence of adjacent faces
	  vector<ftSurface*> adj_faces;
	  face->getAdjacentFaces(adj_faces);
	  if (nmb_bd >= 4 && adj_faces.size() < 4)
	    nmb_bd = (int)adj_faces.size();
	  if (nmb_bd == 3)
	    {
	      // Expects an opposite triangular surface
	      for (int ki=0; ki<(int)sorted_sfs.size(); ++ki)
		{
		  if (!sorted_sfs[ki].get())
		    {
		      int ix = 2*(ki/2) + (1 - (ki%2));
		      if (deg_type[ix] == 1)
			{
			  sorted_sfs[ki] = sfs[0];
			  classification[ki] = make_pair(0, parval[ki]);
			  sfs.erase(sfs.begin());
			  sfs_bd.erase(sfs_bd.begin());
			  deg_type[ki] = 1;
			}
		      else
			{
			  classification[ki] = make_pair(-1, parval[ki]);
			  deg_type[ki] = 2;
			}
		    }
		}
	      if (sfs.size() > 0)
		MESSAGE("Surface sorting configuration not handled");
	    }
	  else //if (nmb_bd == 4)
	    {
	      // Find the non-classified boundaries
	      int ix1=-1, ix2=-1;
	      for (int ki=0; ki<(int)sorted_sfs.size(); ++ki)
		{
		  if (!sorted_sfs[ki].get() && ix1<0)
		    ix1 = ki;
		  else if (!sorted_sfs[ki].get())
		    ix2 = ki;
		}

	      // Check if one of the missing surfaces is adjacent to a degenerate edge
	      // of a non-trimmed triangular surface
	      int ixb;
	      for (ixb=0; ixb<(int)sorted_sfs.size(); ++ixb)
		{
		  if (!sorted_sfs[ixb].get())
		    continue;
		  if (sorted_sfs[ixb]->instanceType() != Class_BoundedSurface &&
			  deg_type[ixb] == 1)
		    {
		      break;
		    }
		}
	      
	      Point deg_pt;
	      if (ixb < (int)sorted_sfs.size())
		{
		  // Identify degenerate edge and compute point on edge
		  bool bottom, right, top, left;
		  sorted_sfs[ixb]->isDegenerate(bottom, right, top, left, toptol_.gap);
		  RectDomain dom = sorted_sfs[ixb]->containingDomain();
		  double upar = left ? dom.umin() : 
		    (right ? dom.umax() : 0.5*(dom.umin()+dom.umax()));
		  double vpar = bottom ? dom.vmin() : 
		    (top ? dom.vmax() : 0.5*(dom.vmin()+dom.vmax()));
		  deg_pt = sorted_sfs[ixb]->point(upar, vpar);
		}

	      if (ix1/2 == ix2/2)
		{
		  // Opposite surfaces are missing
		  MESSAGE("Surface sorting configuration not handled");
		}
	      else
		{
		  // Compute the size of the surfaces opposite to the missing surface
		  int ix3 = 2*(ix1/2) + (1 - (ix1%2));
		  int ix4 = 2*(ix2/2) + (1 - (ix2%2));
		  double len_u1, len_u2, len_v1, len_v2;
		  sorted_sfs[ix3]->estimateSfSize(len_u1, len_v1);
		  sorted_sfs[ix4]->estimateSfSize(len_u2, len_v2);

		  int ix_deg = (len_u1*len_v1 > len_u2*len_v2) ? ix2 : ix1;
		  int ix_other = (len_u1*len_v1 > len_u2*len_v2) ? ix1 : ix2;
		  
		  if (ixb < (int)sorted_sfs.size())
		    {
		      // Identify common vertices between unsorted surface, triangular surface
		      // and sorted surfaces in the parameter direction containing the missing
		      // surfaces
		      int par;
		      int i1, i2, i3, i4;
		      i1 = i2 = i3 = i4 = -1;
		      for (par=0; par<3; ++par)
			{
			  for (ki=0; ki<4; ++ki)
			    if (ix1 == par_sfs[par][ki])
			      {
				i1 = ki;
				break;
			      }
			  if (ki == 4)
			    continue;
			  for (ki=0; ki<4; ++ki)
			    if (ix2 == par_sfs[par][ki])
			      {
				i2 = ki;
				break;
			      }
			  if (ki < 4)
			    break;
			}
		      if (par < 3)
			{
			  for (ki=0; ki<4; ++ki)
			    {
			      if (par_sfs[par][ki] != ix1 && par_sfs[par][ki] != ix2)
				{
				  if (i3 == -1)
				    i3 = ki;
				  else
				    i4 = ki;
				}
			    }
			  // Associate neighbours by letting i3 point to the neighbour of ix1
			  // and i4 to the neighbour of ix2
			  int d1 = ((i1==0 && i3==3) || (i1==3 && i3==0)) ? 1 : abs(i1-i3);
			  int d2 = ((i1==0 && i4==3) || (i1==3 && i4==0)) ? 1 : abs(i1-i4);
			  int d3 = ((i2==0 && i3==3) || (i2==3 && i3==0)) ? 1 : abs(i2-i3);
			  int d4 = ((i2==0 && i4==3) || (i2==3 && i4==0)) ? 1 : abs(i2-i4);
			  if (d2+d3 < d1+d4)
			    std::swap(i3,i4);

			  // Fetch associated faces
			  int ixf1 = shells_[0]->getIndex(sfs[0].get());
			  int ixf2 = shells_[0]->getIndex(sorted_sfs[ixb].get());
			  int ixf3 = shells_[0]->getIndex(sorted_sfs[par_sfs[par][i3]].get());
			  int ixf4 = shells_[0]->getIndex(sorted_sfs[par_sfs[par][i4]].get());
			  shared_ptr<ftSurface> f1 = shells_[0]->getFace(ixf1);
			  shared_ptr<ftSurface> f2 = shells_[0]->getFace(ixf2);
			  shared_ptr<ftSurface> f3 = shells_[0]->getFace(ixf3);
			  shared_ptr<ftSurface> f4 = shells_[0]->getFace(ixf4);
			  vector<shared_ptr<Vertex> > vx3 = f1->getCommonVertices(f2.get(), f3.get());
			  vector<shared_ptr<Vertex> > vx4 = f1->getCommonVertices(f2.get(), f4.get());

			  if (vx3.size() == 1 && vx4.size() == 1)
			    {
			      double dist1 = deg_pt.dist(vx3[0]->getVertexPoint());
			      double dist2 = deg_pt.dist(vx4[0]->getVertexPoint());
			      if (dist1 < toptol_.neighbour && dist1 < dist2)
				{
				  ix_deg = ix1;
				  ix_other = ix2;
				}
			      else if (dist2 < toptol_.neighbour && dist2 < dist1)
				{
				  ix_deg = ix2;
				  ix_other = ix1;
				}
			    }
			}
		    }

		  sorted_sfs[ix_other] = sfs[0];
		  classification[ix_other] = make_pair(0, parval[ix_other]);
		  sfs.erase(sfs.begin());
		  sfs_bd.erase(sfs_bd.begin());
		  deg_type[ix_other] = 0;
		  classification[ix_deg] = make_pair(-1, parval[ix_deg]);
		  deg_type[ix_deg] = 2;
		}
	    }
	}
      else
	{
	  MESSAGE("Surface sorting configuration not handled");
	}
    }
  else if (nmb == nmb_sorted-2)
    {
      // One surface degenerate to a line and one surface degenerate to a point is expected
      size_t kj, kh;
      for (kj=0; kj<deg_type.size(); ++kj)
	if (!(deg_type[kj] < 0 || deg_type[kj] == 1))
	  break;
      for (kh=0; kh<sfs.size(); ++kh)
	if (sfs_bd[kh] != 3)
	  break;
      if (sfs.size() <= 2 && kh == sfs.size() && kj==deg_type.size())
	{
	  // Only triangular surfaces. 
	  // Compute surface sizes to select positions
	  vector<double> sf_size(sorted_sfs.size(), 0.0);
	  for (kj=0; kj<sorted_sfs.size(); ++kj)
	    if (sorted_sfs[kj].get())
	      {
		double usize, vsize;
		sorted_sfs[kj]->estimateSfSize(usize, vsize);
		sf_size[kj] = usize*vsize;
	      }
	  vector<double> sfs_size2(sfs.size());
	  for (kj=0; kj<sfs.size(); ++kj)
	    {
	      double usize2, vsize2;
	      sfs[kj]->estimateSfSize(usize2, vsize2);
	      sfs_size2[kj] = usize2*vsize2;
	    }

	  int min_ix = (sfs_size2.size() == 2 && sfs_size2[1] > sfs_size2[0])
	    ? 1 : 0;
	  double diff = std::numeric_limits<double>::max();
	  int min_diff = -1;
	  for (kj=0; kj<sorted_sfs.size(); ++kj)
	    {
	      if (sorted_sfs[kj].get())
		{
		  double curr_diff = fabs(sf_size[kj]-sfs_size2[min_ix]);
		  if (curr_diff < diff)
		    {
		      diff = curr_diff;
		      min_diff = (int)kj;
		    }
		}
	    }
	  
	  for (kj=0; kj<sorted_sfs.size(); ++kj)
	    {
	      if (!sorted_sfs[kj].get())
		{
		  size_t ix = ((int)kj%2 == 0) ? kj+1 : kj-1;
		  if (sorted_sfs[ix].get() && min_diff == (int)ix)
		    {
		      sorted_sfs[kj] = sfs[min_ix];
		      deg_type[kj] = 1;
		    }
		  else
		    deg_type[kj] = 2;
		  classification[kj] = make_pair(0, parval[kj]);
		}
	    }
	  sfs.erase(sfs.begin()+min_ix);
	  sfs_bd.erase(sfs_bd.begin()+min_ix);

	  if (sfs.size() == 1)
	    {
	      // Place the remaining surface in the direction where
	      // no previous surfaces are positioned
	      for (kj=0; kj<sorted_sfs.size(); kj+=2)
		{
		  if (!(sorted_sfs[kj].get() || sorted_sfs[kj+1].get()))
		    {
		      sorted_sfs[kj] = sfs[0];
		      deg_type[kj] = 1;
		      sfs.erase(sfs.begin());
		      sfs_bd.erase(sfs_bd.begin());
		      break;
		    }
		}
	    }
	  // Select the position of the surface that degenerate to a point. 
	  // Compute size of opposite surface
	  vector<double> opposite_size;
	  vector<int> deg_sfix;
	  for (kj=0; kj<sorted_sfs.size(); kj+=2)
	    {
	      if (deg_type[kj] == 2 && sorted_sfs[kj+1].get())
		{
		  double usize, vsize;
		  sorted_sfs[kj+1]->estimateSfSize(usize, vsize);
		  opposite_size.push_back(usize*vsize);
		  deg_sfix.push_back((int)kj);
		}
	      else if (deg_type[kj+1] == 2 && sorted_sfs[kj].get())
		{
		  double usize, vsize;
		  sorted_sfs[kj]->estimateSfSize(usize, vsize);
		  opposite_size.push_back(usize*vsize);
		  deg_sfix.push_back((int)kj+1);
		}
	    }

	  if (deg_sfix.size() == 2)
	    {
	      deg_type[(opposite_size[0] <= opposite_size[1]) ? deg_sfix[0] : deg_sfix[1]] = 3;
	    }
	}
      else
	MESSAGE("Surface sorting configuration not handled");
    }

//   // Last sorting attempt by looking at intersections between remaining
//   // surfaces and midcurves in the volume running in a parameter direction
//   // with missing surfaces
//     for (ki=0; ki<nmb; ki+=2)
//     {
//       if (!sorted_sfs[ki].get() || !sorted_sfs[ki+1].get())
// 	{
// 	  int nmb_needed = 2 - (sorted_sfs[ki].get() != 0) -
// 	    (sorted_sfs[ki+1].get() != 0);

// 	  if (nmb_needed == (int)sfs.size())
// 	    {
// 	      if (nmb_needed == 1)
// 		std::cout << "Sort surfaces, should be handled previously" << std::endl;
// 	      else
// 		{
// 		  // Sort surfaces
// 		  BoundingBox box1 = sfs[0]->boundingBox();
// 		  Point mid1 = 0.5*(box1.low() + box1.high());
// 		  BoundingBox box2 = sfs[1]->boundingBox();
// 		  Point mid2 = 0.5*(box2.low() + box2.high());
// 		  int ix;
// 		  ix = (mid1[ki/2] <= mid2[ki/2]) ? 0 : 1;
		    
// 		  sorted_sfs[ki] = sfs[ix];
// 		  classification[ki] = make_pair(0, parval[ki]);
// 		  sorted_sfs[ki+1] = sfs[1-ix];
// 		  classification[ki+1] = make_pair(0, parval[ki+1]);
// 		  sfs.clear();
// 		}
// 	    }
// 	  else
// 	    {
// 	      // Fetch midsurface
// 	      int dir1 = (ki/2 + 1)%3;
// 	      double par1 = 0.5*(parval[2*dir1] + parval[2*dir1+1]);
// 	      shared_ptr<ParamSurface> midsurf(vol_->constParamSurface(par1,
// 								       dir1));

// #ifdef DEBUG_VOL1
// 	      std::ofstream of0("midsurf0.g2");
// 	      for (size_t kr=0; kr<sfs.size(); ++kr)
// 		{
// 		  sfs[kr]->writeStandardHeader(of0);
// 		  sfs[kr]->write(of0);
// 		}

// 	      std::ofstream of("midsurf.g2");
// 	      midsurf->writeStandardHeader(of);
// 	      midsurf->write(of);
// #endif


// 	      // Fetch midcurve
// 	      int dir2 = (ki/2 + 2)%3; // ki/2;
// 	      double par2 = 0.5*(parval[2*dir2] + parval[2*dir2+1]);
// 	      vector<shared_ptr<ParamCurve> > midcrvs =
// 		midsurf->constParamCurves(par2, ki <= 1);

// 	      vector<std::pair<int, double> > midcrvints;
// 	      if (midcrvs.size() > 0)
// 		{
// #ifdef DEBUG_VOL1
// 		  midcrvs[0]->writeStandardHeader(of);
// 		  midcrvs[0]->write(of);
// #endif
		  
// 		  // Intersect remaining surfaces with the midcurve
// 		  midcrvints =
// 		    getMidCurveIntersections(midcrvs[0], sfs, toptol_.neighbour);
// 		}
// 	      if (nmb_needed == 1 && (int)midcrvints.size() != nmb_needed)
// 		return false;  // Should have been handled previously
// 	      else if (nmb_needed == 2 && (int)midcrvints.size() != nmb_needed)
// 		{
// 		  // Ambiguous results. Use fall back strategy.
// 		  // Note that this approach may change the orientation
// 		  // of the volume
// 		  // Fetch arbitrary surface
// 		  sorted_sfs[ki] = sfs[0];
// 		  classification[ki] = make_pair(0, 0.0);
// 		  sfs.erase(sfs.begin());

// 		  if (sfs.size() == 0)
// 		    return false;

// 		  // Select other surface as the one surface not 
// 		  // being a neighbour
// 		  int idx = shells_[0]->getIndex(sorted_sfs[ki].get());
// 		  shared_ptr<ftSurface> face0 = shells_[0]->getFace(idx);
// 		  vector<ftSurface*> neighbours;
// 		  face0->getAdjacentFaces(neighbours);

// 		  size_t kh1, kh2;
// 		  for (kh2=0; kh2<sfs.size(); ++kh2)
// 		    {
// 		      for (kh1=0; kh1<neighbours.size(); ++kh1)
// 			if (neighbours[kh1]->surface().get() == sfs[kh2].get())
// 			  break;
// 		      if (kh1 == neighbours.size())
// 			break;
// 		    }
// 		  if (kh2 < sfs.size())
// 		    {
// 		      sorted_sfs[ki+1] = sfs[kh2];
// 		      classification[ki+1] = make_pair(0, 0.0);
// 		      sfs.erase(sfs.begin()+kh2);
// 		    }
// 		  else 
// 		    return false;
// 		}
// 	      else if ((int)midcrvints.size() == nmb_needed)
// 		{
// 		  size_t kj=0;
// 		  if (!sorted_sfs[ki].get())
// 		    {
// 		      sorted_sfs[ki] = sfs[midcrvints[kj].first];
// 		      classification[ki] = make_pair(0, midcrvints[kj].second);
// 		      kj++;
// 		    }

// 		  if (!sorted_sfs[ki+1].get())
// 		    {
// 		      sorted_sfs[ki+1] = sfs[midcrvints[kj].first];
// 		      classification[ki+1] = make_pair(0, midcrvints[kj].second);
// 		      kj++;
// 		    }

// 		  sfs.erase(sfs.begin() + midcrvints[0].first);
// 		  if (midcrvints.size() > 1)
// 		    {
// 		      if (midcrvints[1].first > midcrvints[0].first)
// 			{
// 			  sfs.erase(sfs.begin() + midcrvints[1].first-1);
// 			}
// 		      else
// 			{
// 			  sfs.erase(sfs.begin() + midcrvints[1].first);
// 			}
// 		    }
// 		}
// 	      else
// 		return false;
// 	    }
// 	}
//     }
	  

  if (nmb == nmb_sorted && sfs.size() == 1)
    {
      // Only one free spot and one surface left
      for (int ki=0; ki<(int)sorted_sfs.size(); ++ki)
	{
	  if (!sorted_sfs[ki].get())
	    {
	      sorted_sfs[ki] = sfs[0];
	      classification[ki] = make_pair(0, parval[ki]);
	      sfs.erase(sfs.begin());
	      sfs_bd.erase(sfs_bd.begin());
	    }
	}
    }
      
  
  // Check that the surfaces for each parameter direction is organized
  // in an increasing order
  for (ki=0; ki<nmb; ki+=2)
    if (sorted_sfs[ki].get() && sorted_sfs[ki+1].get() &&
	classification[ki].second > classification[ki+1].second)
      {
	std::swap(sorted_sfs[ki], sorted_sfs[ki+1]);
	std::swap(classification[ki], classification[ki+1]);
      }

  // Check if the sorted surfaces will produce a right handed coordinate system
  // It might be better to evaluate surfaces in a common corner, but for the time
  // being, estimates of the derivative vectors is computed from points in the
  // surfaces
  vector<Point> vec(3);
  double u, v;
  for (ki=0; ki<nmb_sorted; ki+=2)
    {
      if (sorted_sfs[ki].get() && sorted_sfs[ki+1].get())
	{
	  Point pt1 = sorted_sfs[ki]->getInternalPoint(u, v);
	  Point pt2 = sorted_sfs[ki+1]->getInternalPoint(u, v);
	  vec[ki/2] = pt2 - pt1;
	}
      else 
	vec[ki/2] = Point(0.0, 0.0, 0.0);
    }

  if ((vec[0] % vec[1])*vec[2] < 0.0)
    {
      // Turn the first parameter direction that is not fixed
      int kj;
      for (kj=0; kj<nmb_sorted; kj+=2)
	{
	  if (classification[kj].first == 0 && classification[kj+1].first == 0)
	    break;
	}
      if (kj < nmb_sorted)
	{
	  // An appropriate parameter direction is found. Otherwise keep the
	  // left handed coordinate system
	  std::swap(sorted_sfs[kj], sorted_sfs[kj+1]);
	}
    }
      
      
  
  return (sfs.size() == 0); // Sorting succeeded if there is no remaining surfaces
}


//===========================================================================
shared_ptr<SurfaceOnVolume> 
ftVolume::getVolSf(shared_ptr<ParamSurface>& surf) const
//===========================================================================
{
  shared_ptr<SurfaceOnVolume> vol_sf = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(surf);
  if (!vol_sf.get())
    {
      // Look for an underlying surface on volume
      shared_ptr<BoundedSurface> bd_surf =
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
      if (bd_surf.get())
	{
	  shared_ptr<ParamSurface> tmp_surf = bd_surf->underlyingSurface();
	  vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(tmp_surf);
	}
    }

  return vol_sf;
}


//===========================================================================
vector<std::pair<int, double> >
ftVolume::getMidCurveIntersections(shared_ptr<ParamCurve> curve,
				   vector<shared_ptr<ParamSurface> >& sfs,
				   double tol) const
//===========================================================================
{
  vector<std::pair<int, double> > result;

  // Fetch spline curve if available
  shared_ptr<SplineCurve> cv(curve->geometryCurve());
  if (!cv.get())
    return result;  // No intersction results available

#ifdef DEBUG_VOL1
  std::ofstream of0("midsurf_others.g2");
#endif

  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      // Fetch spline surface if available
      SplineSurface* surf = sfs[ki]->getSplineSurface();
      shared_ptr<SplineSurface> surf2;
      if (!surf)
	{
	  // Make approximated surface
	  // First make a surface representing the parameter domain
	  vector<double> knots(8);
	  knots[0] = knots[1] = knots[2] = knots[3] = 0.0;
	  knots[4] = knots[5] = knots[6] = knots[7] = 1.0;
	  vector<double> coefs(48, 0.0);
	  BsplineBasis basis(4, knots.begin(), knots.end());
	  shared_ptr<SplineSurface> tmp_sf = 
	    shared_ptr<SplineSurface>(new SplineSurface(basis, basis, 
							coefs.begin(), 3));
	  
	  // Approximate
	  surf2 = AdaptSurface::approxInSplineSpace(sfs[ki], tmp_sf, 
						    toptol_.gap);
	  surf = surf2.get();
	}

#ifdef DEBUG_VOL1
      surf->writeStandardHeader(of0);
      surf->write(of0);
#endif

      // Check for curve bounded domain
      const CurveBoundedDomain* bdomain = 0;
      shared_ptr<BoundedSurface> bsurf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(sfs[ki]);
      if (bsurf.get())
	bdomain =  &(bsurf->parameterDomain());

#ifdef DEBUG_VOL1
      std::ofstream of("cv_sf_int.g2");
      cv->writeStandardHeader(of);
      cv->write(of);
      surf->writeStandardHeader(of);
      surf->write(of);
#endif

      // Perform intersections
      vector<pair<double,Point> > int_pts;
      vector<pair<pair<double,Point>,pair<double,Point> > > int_crvs;
      vector<int> pretop;
      intersectCurveSurf(cv.get(), surf, tol, int_pts, 
			 pretop, int_crvs);

      // Only considering non-tangential intersections. 
      // Intersection curves are not relevant
      // Remove tangential intersections from the results
      size_t kj;
      for (kj=0; kj<int_pts.size(); )
	{
	  // The pretopology for each point has 4 entries
	  int kr;
	  for (kr=0; kr<4; ++kr)
	    if (pretop[4*kj+kr] == 3)  // SI_ON
 	      break;
	  if (kr < 4)
	    {
	      int_pts.erase(int_pts.begin()+kj);
	      pretop.erase(pretop.begin()+4*kj, pretop.begin()+4*(kj+1));
	    }
	  else
	    kj++;
	}

      if (bdomain)
	{
	  // Remove intersection points outside the domain
	  for (kj=0; kj<int_pts.size(); )
	    {
	      Array<double,2> tmp_pt(int_pts[kj].second[0],
				     int_pts[kj].second[1]);
	      bool in_domain = bdomain->isInDomain(tmp_pt, tol);
	      if (!in_domain)
		{
		  int_pts.erase(int_pts.begin()+kj);
		  pretop.erase(pretop.begin()+4*kj, pretop.begin()+4*(kj+1));
		}
	      else
		kj++;
	    }
	}

      if (int_pts.size() > 0)
	{
	  // At least one intersection point is found. Return arbitrary point
	  result.push_back(make_pair((int)ki, int_pts[0].first));
	}
    }
  return result;
}

//===========================================================================
shared_ptr<ParamVolume> 
ftVolume::createByLoft(shared_ptr<ParamSurface> sf1,
		       shared_ptr<ParamSurface> sf2, 
		       double tol, int pardir)
//===========================================================================
{
  shared_ptr<ParamVolume> volume;

  // Two surfaces. Since they originate from the boundary surfaces
  // of a volume on oppsite sides of the volume, they have swapped
  // parameter directions
  shared_ptr<ParamSurface> sf0(sf1->clone());
  sf0->swapParameterDirection();

  // Get spline space
  shared_ptr<SplineSurface> surf1 =
    dynamic_pointer_cast<SplineSurface,ParamSurface>(sf0);
  if (!surf1.get())
    {
      shared_ptr<SurfaceOnVolume> vol_sf;
      shared_ptr<BoundedSurface> bd_sf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf0);
      if (bd_sf.get())
	{
	  shared_ptr<ParamSurface> tmp_sf = bd_sf->getIsoTrimSurface(tol);
	  if (tmp_sf.get())
	    {
	      surf1 = dynamic_pointer_cast<SplineSurface,ParamSurface>(tmp_sf);
	      vol_sf = dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(tmp_sf);
	    }
	}
      else
	vol_sf = dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(sf0);
      if (!surf1.get() && vol_sf.get())
	{
	  surf1 = 
	    dynamic_pointer_cast<SplineSurface,ParamSurface>(vol_sf->spaceSurface());
	}
    }

  shared_ptr<SplineSurface> surf2 =
    dynamic_pointer_cast<SplineSurface,ParamSurface>(sf2);
  if (!surf2.get())
    {
      shared_ptr<SurfaceOnVolume> vol_sf;
      shared_ptr<BoundedSurface> bd_sf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf2);
      if (bd_sf.get())
	{
	  shared_ptr<ParamSurface> tmp_sf = bd_sf->getIsoTrimSurface(tol);
	  if (tmp_sf.get())
	    {
	      surf2 = dynamic_pointer_cast<SplineSurface,ParamSurface>(tmp_sf);
	      vol_sf = dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(tmp_sf);
	    }
	}
      else
	vol_sf = dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(sf2);
      if (!surf2.get() && vol_sf.get())
	{
	  surf2 = 
	    dynamic_pointer_cast<SplineSurface,ParamSurface>(vol_sf->spaceSurface());
	}
    }

  vector<shared_ptr<SplineSurface> > in_sfs(2);
  shared_ptr<SplineVolume> loft_vol;
  if (surf1.get() && surf2.get())
    {
      // Two spline surfaces. The lofting operation ensures that they
      // have the same spline space
      in_sfs[0] = surf1;
      in_sfs[1] = surf2;
    }
  else if (surf1.get())
    {
      // Approximate the second surface in the spline space of the first one
      surf2 = AdaptSurface::approxInSplineSpace(sf2, surf1, tol);
      in_sfs[0] = surf1;
      in_sfs[1] = surf2;
    }
  else if (surf2.get())
    {
     // Approximate the first surface in the spline space of the second one
      surf1 = AdaptSurface::approxInSplineSpace(sf0, surf2, tol);
      in_sfs[0] = surf1;
      in_sfs[1] = surf2;
    }
  else 
    {
      // Approximate both surfaces in the spline space of the initial volume
      in_sfs = AdaptSurface::expressInSameSplineSpace(sf0, sf2, tol);
    }

  // Create surface
  vector<double> parvals;
  double parlength = 1.0;  // For the time being
  LoftVolumeCreator::makeLoftParams(in_sfs.begin(), 2, parlength, parvals);
  loft_vol = 
    shared_ptr<SplineVolume>(LoftVolumeCreator::loftVolume(in_sfs.begin(), 
							   parvals.begin(), 2));

  // Orient the new volume according to the initial parameter directions
  volume = loft_vol;  // For the time being

  return volume;
}

//===========================================================================
shared_ptr<ParamVolume> 
ftVolume::createByCoons(vector<shared_ptr<ParamSurface> >& sfs, 
			vector<pair<int,double> >& classification,
			vector<int>& deg_type,
			double tol, int degree, bool geom_space)
//===========================================================================
{
  shared_ptr<ParamVolume> result;

  Point deg_pt;
  bool has_deg_pt = identifyDegCorner(sfs, deg_type, deg_pt);

  // Fetch corresponding boundary curve pairs orthogonal to each
  // parameter direction
  vector<vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > > > cvs;
  vector<vector<int> > indices;
  double cvfac = 1.0; //7.5; //15.0;
  bool found = getCoonsCurvePairs(sfs, deg_type, deg_pt, cvs, indices);
  if (!found)
    return result;

  // For each parameter direction, approximate the curves in the same
  // spline space
  int nmb_sample_pr_seg = (geom_space) ? 10 : 5;
  vector<shared_ptr<SplineCurve> > all_cvs;
  for (size_t pardir=0; pardir<cvs.size(); ++pardir)
    {
      vector<shared_ptr<SplineCurve> > coons_cvs(4);
      getCoonsBdCurves(cvs[pardir], indices[pardir], classification,
		       cvfac*tol, degree, coons_cvs, nmb_sample_pr_seg);
      all_cvs.insert(all_cvs.end(), coons_cvs.begin(), coons_cvs.end());
    }

#ifdef DEBUG_VOL1
  // DEBUG
  std::ofstream ofsf("isotrim_sf.g2");
  std::ofstream of("coons_bd.g2");
  for (size_t kr=0; kr<all_cvs.size(); ++kr)
    {
      if (all_cvs[kr].get())
	{
	  all_cvs[kr]->writeStandardHeader(of);
	  all_cvs[kr]->write(of);
	}
    }

  std::ofstream of4("coons_patches.g2");
#endif

  // Make boundary surfaces
  vector<shared_ptr<SplineSurface> > bd_sfs(6);
  int idx[] = {6, 1, 7, 0,  4, 2, 5, 3,  8, 2, 9, 1,  10, 3, 11, 0,  10, 7, 9, 4,  8, 6, 11, 5};
  int ki, kj;
  bool loose_approx;
  for (ki=0; ki<6; ++ki)
    {
      loose_approx = false;

      // First check if the surface exists in an adjacent volume
      // Find surface in shell
      // Make initial coons patch
      int sf_idx = shells_[0]->getIndex(sfs[ki].get());
      shared_ptr<SplineSurface> other_sf;
      bool has_twin = false;
      if (sf_idx >= 0)
	{
	  shared_ptr<ftSurface> face = shells_[0]->getFace(sf_idx);
	  if (face->twin())
	    {
	      has_twin = true;
	      shared_ptr<ParamSurface> twin_sf = 
		face->twin()->asFtSurface()->surface();
	      if (twin_sf->instanceType() == Class_SplineSurface ||
		  twin_sf->instanceType() == Class_SurfaceOnVolume)
		{
		  shared_ptr<SurfaceOnVolume> vol_sf = 
		    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(twin_sf);
		  if (vol_sf.get() && vol_sf->hasSpaceSurface())
		    {
		      shared_ptr<SplineSurface> tmp = 
			dynamic_pointer_cast<SplineSurface,ParamSurface>(vol_sf->spaceSurface());
		      if (tmp.get())
			{
			  other_sf = shared_ptr<SplineSurface>(tmp->clone());
			  other_sf->swapParameterDirection();
			}
		      else
			loose_approx = true;
		    }
		  else if (vol_sf.get() && vol_sf->parPref())
		    loose_approx = true;
		  else
		    {
		      SplineSurface *tmp = twin_sf->getSplineSurface();
		      if ((!tmp->rational()) && tmp->order_u()==degree+1 && 
			  tmp->order_v()==degree+1)
			{
			  other_sf = shared_ptr<SplineSurface>(tmp->clone());
			  other_sf->swapParameterDirection();
			}
		    }
		}
	    }
	}

      if (other_sf.get())
	bd_sfs[ki] = other_sf;
      else
	{
	  vector<shared_ptr<ParamCurve> > bd_cvs(4);
	  int nmb_bd_cvs = 0;
	  for (kj=0; kj<4; ++kj)
	    {
	      bd_cvs[kj] = all_cvs[idx[4*ki+kj]];
	      if (bd_cvs[kj].get())
		nmb_bd_cvs++;
	    }

	  // if (nmb_bd_cvs == 3 && deg_pt.dimension() == 0)
	  //   {
	  //     // A surface with one degenerate boundary will be 
	  //     // created. Make sure that this edge is consistent with
	  //     // the boundary surface configuration
	  //     bool found_deg_pt = identifyDegCorner2(sfs, deg_type,
	  // 					     bd_cvs, deg_pt);
	  //   }

	  // Check and fix orientation
	  vector<shared_ptr<ParamCurve> > dummy_cvs(bd_cvs.size());
	  double max_cv_len =
	    sortCoonsPatchBdCvs(bd_cvs, dummy_cvs, deg_pt, toptol_.neighbour);
	  
#ifdef DEBUG_VOL1
	  std::ofstream of2("curr_coons_bd.g2");
	  for (kj=0; kj<4; ++kj)
	    {
	      bd_cvs[kj]->writeStandardHeader(of2);
	      bd_cvs[kj]->write(of2);
	      Point pnt = bd_cvs[kj]->point(bd_cvs[kj]->startparam());
	      of2 << "400 1 0 4 255 0 0 255" << std::endl;
	      of2 << "1" << std::endl;
	      of2 << pnt << std::endl;
	    }
#endif

	  CurveLoop boundary(bd_cvs, tol);
	  shared_ptr<SplineSurface> init_sf(CoonsPatchGen::createCoonsPatch(boundary));
	  shared_ptr<SplineSurface> surf;

	  // Check if the surface to approximate has the appropriate number
	  // of trimming curves
	  shared_ptr<ParamSurface> orig_sf = sfs[ki];
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sfs[ki]);

	  shared_ptr<SplineSurface> spl_sf =
	    dynamic_pointer_cast<SplineSurface, ParamSurface>(sfs[ki]);
	  if (spl_sf.get() && spl_sf->rational() == false && 
	      spl_sf->order_u()<=degree+1 && spl_sf->order_v()<=degree+1)
	    {
#ifdef DEBUG_VOL1
	      std::cout << "Spline surface " << spl_sf.get() << std::endl;
	      spl_sf->writeStandardHeader(ofsf);
	      spl_sf->write(ofsf);
#endif
	      // TEST accuracy and data size
	      surf = spl_sf;
	    }

	  if (bd_sf.get())
	    {
	      if (bd_sf->isBoundaryTrimmed(toptol_.gap))
		{
#ifdef DEBUG_VOL1
		  std::cout << "Surface " << bd_sf.get() << " is isotrimmed" << std::endl;
		  bd_sf->writeStandardHeader(ofsf);
		  bd_sf->write(ofsf);
#endif
		  shared_ptr<ParamSurface> under = bd_sf->underlyingSurface();
		  shared_ptr<SplineSurface> spl_sf =
		    dynamic_pointer_cast<SplineSurface, ParamSurface>(under);
		  if (spl_sf.get() && spl_sf->rational() == false && 
		      spl_sf->order_u()<=degree+1 && 
		      spl_sf->order_v()<=degree+1)
		    surf = spl_sf;
		  else
		    {
		      shared_ptr<ParameterSurfaceOnVolume> tmp_sf = 
			dynamic_pointer_cast<ParameterSurfaceOnVolume,ParamSurface>(under);
		      if (tmp_sf.get())
			{
			  if (tmp_sf->hasParameterSurface())
			    {
			      shared_ptr<SplineSurface> spl_sf =
				dynamic_pointer_cast<SplineSurface, ParamSurface>(tmp_sf->parameterSurface());
			      if (spl_sf.get() && spl_sf->rational() == false && 
				  spl_sf->order_u()<=degree+1 && 
				  spl_sf->order_v()<=degree+1)
				surf = spl_sf;
			    }
			  if (!surf.get())
			    {
			      // Check if the volume surface is a constant
			      // parameter surface. In that case the coons
			      // patch interpolating the boundary curves
			      // can be used without updating as the requested
			      // surface is planar
			      int constdir = tmp_sf->getConstDir();
			      if (constdir > 0)
				surf = init_sf;
			    }
			}
		      else
			{
			  shared_ptr<SurfaceOnVolume> tmp_sf2 = 
			    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(under);
			  if (tmp_sf2.get())
			    {
			      if (tmp_sf2->hasSpaceSurface())
				{
				  shared_ptr<SplineSurface> spl_sf =
				    dynamic_pointer_cast<SplineSurface, ParamSurface>(tmp_sf2->spaceSurface());
				  if (spl_sf.get() && spl_sf->rational() == false && 
				      spl_sf->order_u()<=degree+1 && 
				      spl_sf->order_v()<=degree+1)
				    surf = spl_sf;
				}
			    }
			}
		    }
		}
	    }

	  if (bd_sf.get() && !surf.get())
	    {
	      // Check if the volume surface is a constant parameter surface.
	      // In that case the coons patch interpolating the boundary curves
	      // can be used without updating as the requested surface is planar
	      shared_ptr<ParamSurface> under = bd_sf->underlyingSurface();
	      shared_ptr<ParameterSurfaceOnVolume> tmp_sf = 
		dynamic_pointer_cast<ParameterSurfaceOnVolume,ParamSurface>(under);
	      if (tmp_sf.get())
		{
		  int constdir = tmp_sf->getConstDir();
		  if (constdir > 0)
		    surf = init_sf;
		}
	    }

	  if (bd_sf.get() && !surf.get())
	    {
	      if (bd_sf->numberOfLoops() != 1 ||
		  bd_sf->loop(0)->size() < 4)
		{
		  // Replace trimming loop
		  vector<shared_ptr<CurveOnSurface> > bd_cvs2(bd_cvs.size());
		  shared_ptr<ParamSurface> under = bd_sf->underlyingSurface();
		  for (size_t kh=0; kh<bd_cvs.size(); ++kh)
		    {
		      if (under->instanceType() == Class_ParameterSurfaceOnVolume &&
			  bd_cvs[kh]->instanceType() != Class_ParameterCurveOnVolume)
			{
			  shared_ptr<SurfaceOnVolume> tmp_sf = 
			    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(under);
			  shared_ptr<ParameterCurveOnVolume> tmp_cv;
			  if (bd_cvs[kh]->instanceType() == Class_CurveOnVolume)
			    {
			      shared_ptr<CurveOnVolume> vol_cv = 
				dynamic_pointer_cast<CurveOnVolume,ParamCurve>(bd_cvs[kh]);
			      tmp_cv = shared_ptr<ParameterCurveOnVolume>(new ParameterCurveOnVolume(vol_cv->underlyingVolume(),
												     vol_cv->parameterCurve(),
												     vol_cv->spaceCurve()));
			    }
			  else
			    {
			      shared_ptr<ParamCurve> dummy_space;
			      tmp_cv = shared_ptr<ParameterCurveOnVolume>(new ParameterCurveOnVolume(tmp_sf->getVolume(),
												     bd_cvs[kh],
												     dummy_space));
			    }
			  bd_cvs2[kh] = 
			    shared_ptr<CurveOnSurface>(new CurveOnSurface(under,
									  tmp_cv,
									  false));
			}
		      else
			bd_cvs2[kh] = 
			  shared_ptr<CurveOnSurface>(new CurveOnSurface(under,
									bd_cvs[kh],
									false));
		      bool success = bd_cvs2[kh]->ensureParCrvExistence(tol);
		      if (!success)
			{
			  // Check if an existing parameter curve can be reused
			  int min_ix = -1;
			  double min_dist = std::numeric_limits<double>::max();
			  bool opposite = false;
			  int nmb = bd_sf->loop(0)->size();
			  Point pos1 = bd_cvs2[kh]->ParamCurve::point(bd_cvs2[kh]->startparam());
			  Point pos2 = bd_cvs2[kh]->ParamCurve::point(bd_cvs2[kh]->endparam());
			  for (int ka=0; ka<nmb; ++ka)
			    {
			      shared_ptr<ParamCurve> loop_cv = (*bd_sf->loop(0))[ka];
			      Point pos3 = loop_cv->point(loop_cv->startparam());
			      Point pos4 = loop_cv->point(loop_cv->endparam());
			      double d1 = pos1.dist(pos3);
			      double d2 = pos1.dist(pos4);
			      double d3 = pos2.dist(pos3);
			      double d4 = pos2.dist(pos4);
			      if (std::min(d1+d4, d2+d3) < min_dist)
				{
				  min_ix = ka;
				  min_dist = std::min(d1+d4, d2+d3);
				  opposite = (d1+d4 > d2+d3);
				}
			    }
			  if (min_ix >= 0 && min_dist < 2.0*tol)
			    {
			      // Use existing curve
			      shared_ptr<ParamCurve> loop_cv = (*bd_sf->loop(0))[min_ix];
			      shared_ptr<CurveOnSurface> sf_cv =
				dynamic_pointer_cast<CurveOnSurface,ParamCurve>(loop_cv);
			      if (sf_cv.get() && sf_cv->hasParameterCurve())
				{
				  shared_ptr<ParamCurve> pcurve(sf_cv->parameterCurve()->clone());
				  if (opposite)
				    pcurve->reverseParameterDirection();
				  bd_cvs2[kh]->setParameterCurve(pcurve);
				}
			      else
				{
#ifdef DEBUG
				  std::cout << "No curve?" << std::endl;
#endif
				  bd_sf.reset();
				  break;
				}
			    }
			  else if (pos1.dist(pos2) < tol)
			    {
			      // Degenerate curve. Find adjacent curves
			      Point par1, par2;
			      for (int ka=0; ka<nmb; ++ka)
				{
				  shared_ptr<ParamCurve> loop_cv = (*bd_sf->loop(0))[ka];
				  Point pos3 = loop_cv->point(loop_cv->startparam());
				  Point pos4 = loop_cv->point(loop_cv->endparam());
				  shared_ptr<CurveOnSurface> sf_cv =
				    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(loop_cv);
				  if (pos3.dist(pos1) < tol && 
				      sf_cv->hasParameterCurve())
				    {
				      if (par1.dimension() > 0)
					par2 = sf_cv->parameterCurve()->point(loop_cv->startparam());
				      else
					par1 = sf_cv->parameterCurve()->point(loop_cv->startparam());
				    }
				  if (pos4.dist(pos1) < tol && 
				      sf_cv->hasParameterCurve())
				    {
				      if (par1.dimension() > 0)
					par2 = sf_cv->parameterCurve()->point(loop_cv->endparam());
				      else
					par1 = sf_cv->parameterCurve()->point(loop_cv->endparam());
				    }
				}
			      if (par1.dimension() == 2 && par2.dimension() == 2)
				{
				  shared_ptr<ParamCurve> pcurve(new SplineCurve(par1, 0.0,
										par2, 1.0));
				  bd_cvs2[kh]->setParameterCurve(pcurve);
				}
				      
			      else
				{
				  bd_sf.reset();
				  break;
				}
			    }
			  else
			    {
			      bd_sf.reset();
			      break;
			    }
			}
		    }
		  if (bd_sf.get())
		    orig_sf =
		      shared_ptr<ParamSurface>(new BoundedSurface(under, bd_cvs2, tol, false));
		}
	    }

	  // Check
	  // double ptol = 1.0e-4;
	    
	  shared_ptr<SurfaceOnVolume> vol_orig = 
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(orig_sf);
	  shared_ptr<ParameterSurfaceOnVolume> vol_origp = 
	    dynamic_pointer_cast<ParameterSurfaceOnVolume,ParamSurface>(orig_sf);
	  if (vol_origp.get() && vol_orig->hasParameterSurface())
	    {
	      shared_ptr<SplineSurface> tmp_sf = 
		dynamic_pointer_cast<SplineSurface,ParamSurface>(vol_orig->parameterSurface());
	      if (tmp_sf.get())
		surf = tmp_sf;
	      else
		loose_approx = true;
	    }
	  else if (vol_orig.get() && vol_orig->hasSpaceSurface())
	    {
	      shared_ptr<SplineSurface> tmp_sf = 
		dynamic_pointer_cast<SplineSurface,ParamSurface>(vol_orig->spaceSurface());
	      if (tmp_sf.get())
		surf = tmp_sf;
	      else
		loose_approx = true;
	    }
	  else if (vol_orig.get() && vol_orig->parPref())
	    loose_approx = true;
	  else if (/*orig_sf->isIsoTrimmed(ptol) && 
		     orig_sf->getSplineSurface() != 0*/ false)
	    {
#ifdef DEBUG_VOL1
	      std::cout << "Iso trimmed spline boundary surface" << std::endl;
#endif
	      RectDomain dom = orig_sf->containingDomain();
	      SplineSurface *tmp = orig_sf->getSplineSurface();
	      if ((!tmp->rational()) && tmp->order_u()==degree+1 && 
		  tmp->order_v()==degree+1)
		{
		  double umin = std::max(dom.umin(), tmp->startparam_u());
		  double umax = std::min(dom.umax(), tmp->endparam_u());
		  double vmin = std::max(dom.vmin(), tmp->startparam_v());
		  double vmax = std::min(dom.vmax(), tmp->endparam_v());
		  surf = shared_ptr<SplineSurface>(tmp->subSurface(umin, vmin, 
								   umax, vmax));
		}
	    }

	  // Approximate initial surface
	  if (!surf.get())
	    {
	      if (orig_sf.get() && max_cv_len > toptol_.neighbour)
		surf = AdaptSurface::adaptSurface(orig_sf, init_sf, 
						  (loose_approx) ? toptol_.neighbour : tol);
	      else
		surf = init_sf;
	    }
#ifdef DEBUG_VOL1
	  std::cout << "Result size: "<< surf->numCoefs_u() << ", ";
	  std::cout << surf->numCoefs_v() << " twin: ";
	  std::cout << has_twin << std::endl;
#endif

	  bd_sfs[ki] = surf;
	}
#ifdef DEBUG_VOL1
      bd_sfs[ki]->writeStandardHeader(of4);
      bd_sfs[ki]->write(of4);
#endif
    }

  // Make Coons volume
  // double fac = 3.0;
  result = 
    shared_ptr<ParamVolume>(CoonsPatchVolumeGen::createCoonsPatch(bd_sfs[0].get(),
								  bd_sfs[1].get(),
								  bd_sfs[2].get(),
								  bd_sfs[3].get(),
								  bd_sfs[4].get(),
								  bd_sfs[5].get(),
								  toptol_.neighbour));
  return result;
}


//===========================================================================
bool
ftVolume::identifyDegCorner(vector<shared_ptr<ParamSurface> >& sfs, 
			    vector<int>& deg_type,
			    Point& deg_pt)
//===========================================================================
{
  if (sfs.size() != 6)
    return false;

  // Check configuration
  int nmb_sfs = 0;
  int nmb_tri = 0;
  bool deg_dir[3];
  deg_dir[0] = deg_dir[1] = deg_dir[2] = false;
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      if (sfs[ki].get())
	{
	  nmb_sfs++;
	  if (deg_type[ki] == 1)
	    {
	      nmb_tri++;
	      deg_dir[ki/2] = true;
	    }
	}
    }
  
  if (nmb_sfs == 4 || nmb_tri == nmb_sfs)
    {
      // Surface configuration: 4 surfaces where all is 3 sided. 
      
      // Indices of side surface for each parameter direction
      int surf_ix[3][4] = {{0, 2, 1, 3}, {1, 5, 0, 4}, {2, 4, 3, 5}};

      // Find the parameter direction where 2 surfaces are defined
      int pardir;
      for (pardir=0; pardir<3; ++pardir)
	{
	  int kj;
	  for (kj=0; kj<4; ++kj)
	    if (deg_type[surf_ix[pardir][kj]] == 3 /* Surface degenerate to point */)
	      break;
	  if (kj == 4)
	    break;
	}

      // The side surfaces corresponding to the identified parameter
      // direction will have one common corner point. Find this point
      vector<vector<pair<Point,Point> > > corners;
      for (int kj=0; kj<4; ++kj)
	{
	  if (sfs[surf_ix[pardir][kj]])
	    {
	      vector<pair<Point,Point> > curr_corners;
	      sfs[surf_ix[pardir][kj]]->getCornerPoints(curr_corners);
	      corners.push_back(curr_corners);
	    }
	}

      if (corners.size() != 3)
	return false;

      // Compute midpoint between the corners of the remaining surface
      vector<pair<Point,Point> > corners2;
      for (size_t ki=0; ki<sfs.size(); ++ki)
	{
	  if (!sfs[ki].get())
	    continue;
	  int kj;
	  for (kj=0; kj<4; ++kj)
	    if ((int)ki == surf_ix[pardir][kj])
	      break;
	  if (kj == 4)
	    {
	      sfs[ki]->getCornerPoints(corners2);
	      break;
	    }
	}

      Point mid(0.0, 0.0, 0.0);
      for (size_t ki=0; ki<corners2.size(); ++ki)
	mid += corners2[ki].first;
      mid /= (double)corners2.size();

      // For each surface in the selected parameter direction, identify
      // the corner most distant from the remaining surface
      vector<int> ix(corners.size(), 0);
      for (size_t ki=0; ki<corners.size(); ++ki)
	{
	  int max_ix = 0;
	  double max_dist = 0.0;
	  for (size_t kj=0; kj<corners[ki].size(); ++kj)
	    {
	      double dist = corners[ki][kj].first.dist(mid);
	      if (dist > max_dist)
		{
		  max_dist = dist;
		  max_ix = (int)kj;
		}
	    }
	  ix[ki] = max_ix;
	}

      deg_pt = Point(0.0, 0.0, 0.0);
      for (size_t ki=0; ki<corners.size(); ++ki)
	deg_pt += corners[ki][ix[ki]].first;
      deg_pt /= (double)corners.size();
      // for (size_t ki=0; ki<corners.size(); ++ki)
      // 	if (corners[ki].size() != 3)
      // 	  return false;
 
      // int min_ix1=0, min_ix2=0, min_ix3=0;
      // double min_dist = std::numeric_limits<double>::max();
      // for (int ka1=0; ka1<(int)corners[0].size(); ++ka1)
      // 	for (int ka2=0; ka2<(int)corners[1].size(); ++ka2)
      // 	  for (int ka3=0; ka3<(int)corners[2].size(); ++ka3)
      // 	    {
      // 	      double dist = corners[0][ka1].first.dist(corners[1][ka2].first) +
      // 		corners[0][ka1].first.dist(corners[2][ka3].first) +
      // 		corners[1][ka2].first.dist(corners[2][ka3].first);
      // 	      if (dist < min_dist)
      // 		{
      // 		  min_dist = dist;
      // 		  min_ix1 = ka1;
      // 		  min_ix2 = ka2; 
      // 		  min_ix3 = ka3;
      // 		}
      // 	    }

      // deg_pt = (corners[0][min_ix1].first + 
      // 		corners[1][min_ix2].first + 
      // 		corners[2][min_ix3].first)/3.0;

      return true;
    }
  else if (nmb_tri == 2)
    {
      // Check if there are triangular surface in more than one parameter
      // direction
      int nmb_deg_dir = 0;
      for (int kj=0; kj<3; ++kj)
	if (deg_dir[kj])
	  nmb_deg_dir++;

      if (nmb_deg_dir == 2)
	{
	  // Identify common degenerate corner between triangular surfaces
	  vector<vector<pair<Point,Point> > > corners;
	  for (size_t ki=0; ki<sfs.size(); ++ki)
	    {
	      if (sfs[ki].get() && deg_type[ki] == 1)
		{
		  vector<pair<Point,Point> > curr_corners;
		  sfs[ki]->getCornerPoints(curr_corners);
		  corners.push_back(curr_corners);
		}
	    }
	  if (corners.size() != 2)
	    return false;
	  int min_ix1=0, min_ix2=0;
	  double min_dist = std::numeric_limits<double>::max();
	  for (int ka1=0; ka1<(int)corners[0].size(); ++ka1)
	    for (int ka2=0; ka2<(int)corners[1].size(); ++ka2)
	      {
		double dist = corners[0][ka1].first.dist(corners[1][ka2].first);
		if (dist < min_dist)
		  {
		    min_dist = dist;
		    min_ix1 = ka1;
		    min_ix2 = ka2; 
		  }
	      }
	  deg_pt = (corners[0][min_ix1].first + 
		    corners[1][min_ix2].first)/2.0;
	  return true;
	}
      else
	return false;
    }
  else
    return false;
}

//===========================================================================
bool
ftVolume::identifyDegCorner2(vector<shared_ptr<ParamSurface> >& sfs, 
			     vector<int>& deg_type,
			     vector<shared_ptr<ParamCurve> >& bd_cvs,
			     Point& deg_pt)
//===========================================================================
{
  if (sfs.size() != 6)
    return false;

  // Check configuration
  int nmb_sfs = 0;
  int nmb_tri = 0;
  int nmb_bi = 0;
  vector<bool> deg_dir(3, false);
  size_t ki;
  for (ki=0; ki<sfs.size(); ++ki)
    {
      if (sfs[ki].get())
	{
	  nmb_sfs++;
	  if (deg_type[ki] == 1)
	    {
	      nmb_tri++;
	      deg_dir[ki/2] = true;
	    }
	  else if (deg_type[ki] == 2)
	    {
	      nmb_bi++;
	      deg_dir[ki/2] = true;
	    }
	}
    }

  if (!(nmb_sfs == 5 && nmb_tri == 2 && nmb_bi == 1))
    return false;

  // Compute endpoints of the surface which has degenerated to a curve
  // First find the parameter direction with non-degenerate surfaces
  for (ki=0; ki<deg_dir.size(); ++ki)
    if (!deg_dir[ki])
      break;

  if (ki == deg_dir.size())
    return false;  // Unexpected configuration

  // Identify two corner points which is common between the surfaces
  // in the non-degenerate direction
  vector<pair<Point,Point> > corner1, corner2;
  sfs[2*ki]->getCornerPoints(corner1);
  sfs[2*ki+1]->getCornerPoints(corner2);
  double mindist1 = std::numeric_limits<double>::max(), mindist2 = std::numeric_limits<double>::max();
  pair<int,int> min_ix1 = make_pair(-1,-1);
  pair<int,int> min_ix2 = make_pair(-1,-1);
  size_t kr, kh;
  for (kr=0; kr<corner1.size(); ++kr)
    {
      for (kh=0; kh<corner2.size(); ++kh)
	{
	  double dist = corner1[kr].first.dist(corner2[kh].first);
	  if (dist < mindist1)
	    {
	      min_ix2 = min_ix1;
	      mindist2 = mindist1;
	      min_ix1 = make_pair((int)kr,(int)kh);
	      mindist1 = dist;
	    }
	  else if (dist < mindist2)
	    {
	      min_ix2 = make_pair((int)kr,(int)kh);
	      mindist2 = dist;
	    }
	}
    }

  vector<Point> deg_cand;
  if (min_ix2.first >= 0)
    deg_cand.push_back(0.5*(corner1[min_ix2.first].first+
			    corner2[min_ix2.second].first));
  if (min_ix1.first >= 0)
    deg_cand.push_back(0.5*(corner1[min_ix1.first].first+
			    corner2[min_ix1.second].first));

  if (deg_cand.size() == 0)
    return false;

  // Select candidate closest to the current boundary curves
  double mindist = std::numeric_limits<double>::max();
  int min_ix = -1;
  for (size_t kj=0; kj<bd_cvs.size(); ++kj)
    {
      if (!bd_cvs[kj].get())
	continue;
      Point pos1 = bd_cvs[kj]->point(bd_cvs[kj]->startparam());
      Point pos2 = bd_cvs[kj]->point(bd_cvs[kj]->endparam());
      for (ki=0; ki<deg_cand.size(); ++ki)
	{
	  double dist = std::min(pos1.dist(deg_cand[ki]), 
				 pos2.dist(deg_cand[ki]));
	  if (dist < mindist)
	    {
	      mindist = dist;
	      min_ix = (int)ki;
	    }
	}
    }
  if (min_ix >= 0)
    deg_pt = deg_cand[min_ix];
  
  return (min_ix >= 0);
}

//===========================================================================
void
ftVolume::getCoonsBdCurves(vector<pair<shared_ptr<ParamCurve>,shared_ptr<ParamCurve> > >& cvs,
			   vector<int>& indices,
			   vector<pair<int,double> >& classification,
			   double tol, int degree, 
			   vector<shared_ptr<SplineCurve> >& coons_cvs,
			   int nmb_sample_pr_seg)
//===========================================================================
{
  // Remove curves already belonging to final surfaces
  vector<shared_ptr<ParamCurve> > init_cvs;
  vector<BsplineBasis> crv_basis;
  vector<int> bb_idx;
  int ki;
  static int min_cont = 2; //1;
  //int max_coef = 10;
  int max_basis = 0;
  for (ki=1; ki<4; ++ki)
    {
      int ix1 = indices[ki-1];
      int ix2 = indices[ki];
      shared_ptr<SplineCurve> spcv1, spcv2;
      shared_ptr<ParamSurface> sf1, sf2;
      if (cvs[ki].first.get())
	{
	    spcv1 = dynamic_pointer_cast<SplineCurve,ParamCurve>(cvs[ki].first);
	  if (!spcv1.get())
	    {
	      shared_ptr<CurveOnSurface> sfcv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cvs[ki].first);
	      if (sfcv.get() && sfcv->isConstantCurve())
		{
		  spcv1 = dynamic_pointer_cast<SplineCurve,ParamCurve>(sfcv->spaceCurve());
		  sf1 = sfcv->underlyingSurface();
		}
	    }
	}
      // Check continuity
      if (spcv1.get() && (spcv1->basis().getMinContinuity() < min_cont /*spcv1->order()-2*/ ||
			  spcv1->order() == degree-1 || spcv1->rational()))
	spcv1.reset();
      
      // Check context
      if (spcv1.get() && classification[ix1].first == 0 && sf1.get())
	{
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(sf1);
	  if (vol_sf.get() && vol_sf->parPref())
	    spcv1.reset();
	}

      if (cvs[ki].second.get())
	{
	  spcv2 = dynamic_pointer_cast<SplineCurve,ParamCurve>(cvs[ki].second);
	  if (!spcv2.get())
	    {
	      shared_ptr<CurveOnSurface> sfcv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cvs[ki].second);
	      if (sfcv.get() && sfcv->isConstantCurve())
		{
		  spcv2 = dynamic_pointer_cast<SplineCurve,ParamCurve>(sfcv->spaceCurve());
		  sf2 = sfcv->underlyingSurface();
		}
	    }
	}
      // Check continuity
      if (spcv2.get() && (spcv2->basis().getMinContinuity() < min_cont /*spcv2->order()-2*/ ||
			  spcv2->order() == degree-1 || spcv2->rational()))
	spcv2.reset();
      
      // Check context
      if (spcv2.get() && classification[ix2].first == 0 && sf2.get())
	{
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(sf2);
	  if (vol_sf.get() && vol_sf->parPref())
	    spcv2.reset();
	}

      //      TEST
      spcv1.reset();
      spcv2.reset();

      if ((spcv1.get() && spcv2.get() && spcv1->numCoefs() <= spcv2->numCoefs())
	  || (spcv1.get() && !spcv2.get()))
	{
	  coons_cvs[ki] = spcv1;
	  crv_basis.push_back(spcv1->basis());
	  max_basis = std::max(max_basis, spcv1->numCoefs());
	  bb_idx.push_back(ki);
	}
      else if (spcv2.get())
	{
	  coons_cvs[ki] = spcv2;
	  crv_basis.push_back(spcv2->basis());
	  max_basis = std::max(max_basis, spcv2->numCoefs());
	  bb_idx.push_back(ki);
	}
      
      // Always store initial curves. They will be removed later if they
      // are kept without approximation
      if (cvs[ki].first.get() && !(spcv2.get() &&
				   coons_cvs[ki].get() == spcv2.get()))
	init_cvs.push_back(cvs[ki].first);
      else
	init_cvs.push_back(cvs[ki].second);
    }

  int ix1 = indices[ki-1];
  int ix2 = indices[0];
  shared_ptr<SplineCurve> spcv1, spcv2;
  shared_ptr<ParamSurface> sf1, sf2;
  if (cvs[0].first.get())
    {
      spcv1 = dynamic_pointer_cast<SplineCurve,ParamCurve>(cvs[0].first);
      if (!spcv1.get())
	{
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cvs[0].first);
	  if (sfcv.get() && sfcv->isConstantCurve())
	    {
	      spcv1 = dynamic_pointer_cast<SplineCurve,ParamCurve>(sfcv->spaceCurve());
	      sf1 = sfcv->underlyingSurface();
	    }
	}
    }
  // Check continuity
  if (spcv1.get() && (spcv1->basis().getMinContinuity() < min_cont /*spcv1->order()-2*/ ||
		      spcv1->order() == degree-1 || spcv1->rational()))
    spcv1.reset();
      
      // Check context
   if (spcv1.get() && classification[ix1].first == 0 && sf1.get())
    {
      shared_ptr<SurfaceOnVolume> vol_sf = 
	dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(sf1);
      if (vol_sf.get() && vol_sf->parPref())
	spcv1.reset();
    }

  if (cvs[0].second.get())
    {
      spcv2 = dynamic_pointer_cast<SplineCurve,ParamCurve>(cvs[0].second);
      if (!spcv2.get())
	{
	  shared_ptr<CurveOnSurface> sfcv =
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cvs[0].second);
	  if (sfcv.get() && sfcv->isConstantCurve())
	    {
	      spcv2 = dynamic_pointer_cast<SplineCurve,ParamCurve>(sfcv->spaceCurve());
	      sf2 = sfcv->underlyingSurface();
	    }
	}
    }
  // Check continuity
  if (spcv2.get() && (spcv2->basis().getMinContinuity() < min_cont /*spcv2->order()-2*/ ||
		      spcv2->order() == degree-1 || spcv2->rational()))
    spcv2.reset();
      
  // Check context
   if (spcv2.get() && classification[ix2].first == 0 && sf2.get())
    {
      shared_ptr<SurfaceOnVolume> vol_sf = 
	dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(sf2);
      if (vol_sf.get() && vol_sf->parPref())
	spcv2.reset();
    }

   // TEST
   spcv1.reset();
   spcv2.reset();

  if ((spcv1.get() && spcv2.get() && spcv1->numCoefs() <= spcv2->numCoefs())
      || (spcv1.get() && !spcv2.get()))
    {
      coons_cvs[0] = spcv1;
      crv_basis.push_back(spcv1->basis());
      max_basis = std::max(max_basis, spcv1->numCoefs());
      bb_idx.insert(bb_idx.begin(), 0);
    }
  else if (spcv2.get())
    {
      coons_cvs[0] = spcv2;
      crv_basis.push_back(spcv2->basis());
      max_basis = std::max(max_basis, spcv2->numCoefs());
      bb_idx.insert(bb_idx.begin(), 0);
    }

  // Always store initial curves. They will be removed later if they
  // are kept without approximation
  if (cvs[0].first.get() && !(spcv2.get() &&
			      coons_cvs[0].get() == spcv2.get()))
    init_cvs.insert(init_cvs.begin(), cvs[0].first);
  else
    init_cvs.insert(init_cvs.begin(), cvs[0].second);

  // Remove pointers to missing curves and remember index
  vector<int> missing_ix;
  for (size_t kj=0; kj<init_cvs.size();)
    {
      if (!init_cvs[kj].get())
	{
	  missing_ix.push_back((int)(kj+missing_ix.size()));
	  init_cvs.erase(init_cvs.begin()+kj);
	}
      else
	++kj;
    }

  // Approximate curves in the same spline space up to possible
  // refinements
  vector<shared_ptr<SplineCurve> > app_cvs;
  if (crv_basis.size() > 0/* && max_basis < max_coef*/)
    {
      // Define initial knot vector as the union of the existing
      // knot vectors
      double start = crv_basis[0].startparam();
      double end = crv_basis[0].endparam();
      int order = crv_basis[0].order();
      for (ki=1; ki<(int)crv_basis.size(); ++ki)
	{
	  order = std::max(order, crv_basis[ki].order());
	  start += crv_basis[ki].startparam();
	  end += crv_basis[ki].endparam();
	}
      start /= (double)(crv_basis.size());
      end /= (double)(crv_basis.size());
      for (ki=0; ki<(int)crv_basis.size(); ++ki)
	{
	  crv_basis[ki].rescale(start, end);
	  if (crv_basis[ki].order() < order)
	    crv_basis[ki].increaseOrder(order);
	}

      vector<double> knots;
      GeometryTools::makeUnionKnots(crv_basis, tol, knots);

      // Check the distribution of knotw in the union knot vector
      int nmb_basis = (int)knots.size()-order;
      double tdel = (knots[nmb_basis] - knots[order-1])/(double)nmb_basis;
      tdel /= 5.0;
      double prev = knots[order-1];
      int kj;
      for (kj=order; kj<=nmb_basis; ++kj)
	if (knots[kj] > prev)
	  {
	    if (knots[kj] - prev < tdel)
	      break;
	    else
	      prev = knots[kj];
	  }
      if (kj <= nmb_basis)
	{
	  // Bad knot distribution. Use at most one kept curve
	  size_t kj;
	  if (bb_idx.size() == 1)
	    {
	      knots.erase(knots.begin()+order, knots.begin()+nmb_basis-1);
	      coons_cvs[bb_idx[0]].reset();
	    }
	  else
	    {
	      knots.clear();
	      knots.insert(knots.begin(), coons_cvs[bb_idx[0]]->basis().begin(),
			   coons_cvs[bb_idx[0]]->basis().end());
	      for (kj=0; kj<init_cvs.size(); ++kj)
		if (init_cvs[kj].get() == cvs[bb_idx[0]].first.get() ||
		    init_cvs[kj].get() == cvs[bb_idx[0]].second.get())
		  {
		    init_cvs.erase(init_cvs.begin()+kj);
		    break;
		  }
	      for (kj=1; kj<bb_idx.size(); ++kj)
		coons_cvs[bb_idx[kj]].reset();
	    }
	}
      else
	for (size_t kj=0; kj<bb_idx.size(); ++kj)
	  {
	    size_t kh;
	    for (kh=0; kh<init_cvs.size(); ++kh)
	      if (init_cvs[kh].get() == cvs[bb_idx[kj]].first.get() ||
		  init_cvs[kh].get() == cvs[bb_idx[kj]].second.get())
		  {
		    init_cvs.erase(init_cvs.begin()+kh);
		    break;
		  }
	  }
    
      
      // TEST
      vector<double> knots2;
      knots2.insert(knots2.end(), knots.begin(), knots.begin()+order);
      knots2.insert(knots2.end(), knots.end()-order, knots.end());
      BsplineBasis init_basis((int)knots2.size()-order, order, knots2.begin());
      
      //BsplineBasis init_basis((int)knots.size()-order, order, knots.begin());

	  
      app_cvs = AdaptSurface::curveApprox(&init_cvs[0], 
					  (int)init_cvs.size(),
					  init_basis, tol,
					  nmb_sample_pr_seg);
    }
  else if (init_cvs[0].get())
    app_cvs = AdaptSurface::curveApprox(&init_cvs[0], 
					(int)init_cvs.size(), 
					tol, degree,
					nmb_sample_pr_seg);

  // Collect final curves
  int kj;
  for (ki=0, kj=0; ki<4; ++ki)
    {
      if (coons_cvs[ki].get())
	continue;
      size_t ka;
      for (ka=0; ka<missing_ix.size(); ++ka)
	if (ki == missing_ix[ka])
	  break;
      if (ka < missing_ix.size())
	continue;
      coons_cvs[ki] = app_cvs[kj++];
    }
}


      
//===========================================================================
bool
ftVolume::getCoonsCurvePairs(vector<shared_ptr<ParamSurface> >& sfs, 
			     vector<int>& deg_type, Point& deg_pt,
			     vector<vector<pair<shared_ptr<ParamCurve>, 
			     shared_ptr<ParamCurve> > > >& curves,
			     vector<vector<int> >& indices)
//===========================================================================
{
#ifdef DEBUG_VOL1
  std::ofstream incoons("incoons.g2");
  for (size_t kf=0; kf<sfs.size(); kf++)
    {
      if (sfs[kf].get())
	{
	  sfs[kf]->writeStandardHeader(incoons);
	  sfs[kf]->write(incoons);
	}
    }
  vector<shared_ptr<Vertex> > vx;
  shells_[0]->getAllVertices(vx);
  incoons << "400 1 0 4 255 0 0 255" << std::endl;
  incoons << vx.size() << std::endl;
  for (size_t kf=0; kf<vx.size(); ++kf)
    incoons << vx[kf]->getVertexPoint() << std::endl;
#endif

  // Inidices of side surface for each parameter direction
  int surf_ix[3][4] = {{0, 2, 1, 3}, {1, 5, 0, 4}, {2, 4, 3, 5}};

  int nopt_dir = -1;
  double deg_tol = toptol_.neighbour;
  double a_tol = 1.0e-8;
  if (deg_pt.dimension() > 0)
    {
      // Three surfaces should have a corner at the degenerate point. 
      // Compute minimum distances between the surfaces and the degenerate
      // point
      vector<double> min_deg_dist(sfs.size());
      // Identify parameter direction where no surface in the loop
      // degenerates to a point
      int nmb_sfs = 0;
      size_t kj=0;
      for (size_t ki=0; ki<sfs.size(); ++ki)
	if (sfs[ki].get())
	  {
	    nmb_sfs++;

	    // Fetch surface corners
	    vector<pair<Point,Point> > corners;
	    sfs[ki]->getCornerPoints(corners);

	    double min_dist = std::numeric_limits<double>::max();
	    for (size_t ki=0; ki<corners.size(); ++ki)
	      {
		double dist = deg_pt.dist(corners[ki].first);
		min_dist = std::min(min_dist, dist);
	      }
	    min_deg_dist[kj++] = min_dist;
	  }

      std::sort(min_deg_dist.begin(),min_deg_dist.begin()+nmb_sfs);
      deg_tol = std::max(deg_tol, min_deg_dist[2]+a_tol);
      deg_tol = std::min(deg_tol, 2.0*toptol_.neighbour);

      if (nmb_sfs == 4)
	{
	  for (nopt_dir=0; nopt_dir<3; ++nopt_dir)
	    {
	      int ka;
	      for (ka=0; ka<4; ++ka)
		if (deg_type[surf_ix[nopt_dir][ka]] == 3 /* Surface degenerate to point */)
		  break;

	      if (ka == 4)
		break;
	    }
	}
    }

  // For each volumetric parameter directions
  int ki, kj;
  int pardir;
  int idx;
  int use_curve = 0;  // 0 = both curves, 1 = first, 2 = second
  shared_ptr<ParamCurve> dummy;
  // double ptol = 1.0e-8;
  for (pardir=0; pardir<3; ++pardir)
    {
      // Fetch the first surface
      vector<int> sf_idx(4);
      shared_ptr<ParamSurface> sf1 = sfs[surf_ix[pardir][0]]; 
      // Check if the surface exists. Otherwise fetch the previous one
      bool deg1 = false;
      if (!sf1.get())
	{
	  sf1 = sfs[surf_ix[pardir][3]]; 
	  deg1 = true;
	}
      shared_ptr<ParamSurface> sf0 = sf1;
      int ix1 = shells_[0]->getIndex(sf1.get());
      int ix0 = ix1;

      vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > > bd_cvs(4);
      bool turned = false;

      for (kj=1; kj<=4; ++kj)
	{
	  idx = (kj == 4) ? surf_ix[pardir][0] : surf_ix[pardir][kj];
	  ki = (kj == 4) ? 0 : kj;
	  sf_idx[ki] = idx;

	  use_curve = 0;

	  // Fetch the next surface around the volume
	  shared_ptr<ParamSurface> sf2 = (ki == 0) ? sf0 : sfs[idx];
	  bool deg2 = false;
	  if (!sf2.get() || sf2.get() == sf1.get())
	    {
	      deg2 = true;
	      if (sf1.get() && sfs[surf_ix[pardir][(kj+1)%4]].get())
		sf2 = sfs[surf_ix[pardir][(kj+1)%4]];
	    }
	  else if (kj==4 && sf0.get() != sfs[pardir].get())
	    deg2 = true;

#ifdef DEBUG_VOL1
	  std::ofstream of("co1.g2");
	  if (sf1.get() && sf2.get())
	    {
	      shared_ptr<SurfaceOnVolume> volsf1 =
		dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(sf1);
	      shared_ptr<SurfaceOnVolume> volsf2 =
		dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(sf2);
	      if (volsf1.get() && volsf2.get())
		{
		  volsf1->spaceSurface()->writeStandardHeader(of);
		  volsf1->spaceSurface()->write(of);
		  volsf2->spaceSurface()->writeStandardHeader(of);
		  volsf2->spaceSurface()->write(of);
		}
	      else
		{
		  sf1->writeStandardHeader(of);
		  sf1->write(of);
		  sf2->writeStandardHeader(of);
		  sf2->write(of);
		}
	    }
#endif

	  // Find neighbourhood information
	  int ix2 = /*(ki == 0) ? ix0 :*/ shells_[0]->getIndex(sf2.get());

	  bool adjacent = false;
	  shared_ptr<ftEdge> e1, e2;
	  shared_ptr<ftSurface> f1, f2;
	  if (sf1.get() && sf2.get())
	    {
	      f1 = shells_[0]->getFace(ix1);
	      f2 = shells_[0]->getFace(ix2);
	      adjacent = f1->areNeighbours(f2.get(), e1, e2);
	    }
	  if (!adjacent)
	    {
	      if (deg1 || deg2)
		{
		  if (deg1)
		    {
		      ix1 = ix2;
		      sf1 = sf2;
		      if (kj == 1)
			sf0 = sf2;
		    }
		  deg1 = true;
		  if (!sf2.get())
		    sf1 = sf2;
		  continue;
		}
	      else
		{
		  tpTolerances tol = shells_[0]->getTolerances();
		  int nmb_bd1 = f1->nmbOuterBdCrvs(tol.gap, tol.neighbour, tol.bend);
		  int nmb_bd2 = f2->nmbOuterBdCrvs(tol.gap, tol.neighbour, tol.bend);
		  // Check also with the existence of adjacent faces
		  vector<ftSurface*> adj_faces1;
		  f2->getAdjacentFaces(adj_faces1);
		  if (nmb_bd1 >= 4 && adj_faces1.size() < 4)
		    nmb_bd1 = (int)adj_faces1.size();
		  vector<ftSurface*> adj_faces2;
		  f2->getAdjacentFaces(adj_faces2);
		  if (nmb_bd2 >= 4 && adj_faces2.size() < 4)
		    nmb_bd2 = (int)adj_faces2.size();

		  if (nmb_bd1 < 4 && nmb_bd2 < 4)
		    {
		      // Check for corner adjacency
		      vector<shared_ptr<Vertex> > vx = f1->getCommonVertices(f2.get());
		      if (vx.size() > 0)
			{
			  ix1 = ix2;
			  sf1 = sf2;
			  continue;
			}
		      else
			return false;
		    }
		  else
		    return false;  // Inconsistent information
		}
	    }

	  // Check if more then one edge corresponds to the boundary
	  // curve between the two faces
	  // 13.03.2011. May there also be more than one curve?
	  shared_ptr<ParamCurve> cv1 = e1->geomCurve();
	  shared_ptr<ParamCurve> cv2 = e2->geomCurve();
	  double tmin1 = e1->tMin();
	  double tmax1 = e1->tMax();
	  double tmin2 = e2->tMin();
	  double tmax2 = e2->tMax();
	  vector<double> parmin1, parmin2, parmax1, parmax2;
	  vector<shared_ptr<ParamCurve> > crvs1, crvs2;
	  parmin1.push_back(tmin1);
	  parmin2.push_back(tmin2);
	  parmax1.push_back(tmax1);
	  parmax2.push_back(tmax2);
	  crvs1.push_back(cv1);
	  crvs2.push_back(cv2);
	  size_t idx1=0, idx2=0;
	  ftEdgeBase *e1_1 = e1.get(), *e1_2 = e1.get(), *e2_1 = e2.get(), *e2_2 = e2.get();
	  while (e2_1)
	    {
	      ftEdgeBase* e3 = NULL;
	      if (e1_1->next()->twin() == e2_1->next())
		e3 = e2_1->next();
	      else if (e1_1->next()->twin() == e2_1->prev())
		e3 = e2_1->prev();
	      if (e3 /*&& 
		  (e1_1->next()->geomEdge()->geomCurve().get() == cv1.get() ||
		  e3->geomEdge()->geomCurve().get() == cv2.get())*/)
		{
		  Point p1 = e1_1->point(e1_1->tMax());
		  Point p2 = e1_1->next()->point(e1_1->next()->tMin());
		  Point tan1 = e1_1->tangent(e1_1->tMax());
		  Point tan2 = e1_1->next()->tangent(e1_1->next()->tMin());

		  double dist = p1.dist(p2);
		  double ang = tan1.angle(tan2);
		  if (dist<toptol_.neighbour /*&& ang<toptol_.bend*/)
		    {
		      shared_ptr<ParamCurve> c1 = e1_1->next()->geomEdge()->geomCurve();
		      shared_ptr<ParamCurve> c2 = e3->geomEdge()->geomCurve();
		      double t1 = e1_1->next()->tMin();
		      double t2 = e1_1->next()->tMax();
		      double t3 = e3->tMin();
		      double t4 = e3->tMax();
		      if (c1.get() != cv1.get() || t1 > tmax1+toptol_.neighbour ||
			  t2 < tmin1-toptol_.neighbour)
			{
			  crvs1.push_back(c1);
			  parmin1.push_back(t1);
			  parmax1.push_back(t2);
			  tmin1 = t1;
			  tmax1 = t2;
			  idx1++;
			  cv1 = c1;
			}
		      else
			{
			  tmin1 = parmin1[idx1] = std::min(parmin1[idx1], t1);
			  tmax1 = parmax1[idx1] = std::max(parmax1[idx1], t2);
			}

		      if (c2.get() != cv2.get() || t3 > tmax2+toptol_.neighbour ||
			  t4 < tmin1 - toptol_.neighbour)
			{
			  if (e3 == e2_1->prev())
			    {
			      crvs2.insert(crvs2.begin()+idx2, c2);
			      parmin2.insert(parmin2.begin()+idx2, t3);
			      parmax2.insert(parmax2.begin()+idx2, t4);
			    }
			  else
			    {
			      crvs2.push_back(c2);
			      parmin2.push_back(t3);
			      parmax2.push_back(t4);
			      idx2++;
			    }
			  tmin2 = t3;
			  tmax2 = t4;
			  cv2 = c2;
			}
		      else
			{
			  tmin2 = parmin2[idx2] = std::min(parmin2[idx2], t3);
			  tmax2 = parmax2[idx2] = std::max(parmax2[idx2], t4);
			}
		    }
		  else 
		    e3 = NULL;
		}
	      e1_1 = e1_1->next();
	      e2_1 = e3;
	    }

	  idx1 = idx2 = 0;
	  while (e2_2)
	    {
	      ftEdgeBase* e3 = NULL;
	      if (e1_2->prev()->twin() == e2_2->next())
		{
		  e3 = e2_2->next();
		  idx2 = crvs2.size() - 1;
		}
	      else if (e1_2->prev()->twin() == e2_2->prev())
		e3 = e2_2->prev();
	      if (e3 /*&& 
		  (e1_2->prev()->geomEdge()->geomCurve().get() == cv1.get() ||
		  e3->geomEdge()->geomCurve().get() == cv2.get())*/)
		{
		  Point p1 = e1_2->prev()->point(e1_2->prev()->tMax());
		  Point p2 = e1_2->point(e1_2->tMin());
		  Point tan1 = e1_2->prev()->tangent(e1_2->prev()->tMax());
		  Point tan2 = e1_2->tangent(e1_2->tMin());

		  double dist = p1.dist(p2);
		  double ang = tan1.angle(tan2);
		  if (dist<toptol_.neighbour /*&& ang<toptol_.bend*/)
		    {
		      shared_ptr<ParamCurve> c1 = e1_2->prev()->geomEdge()->geomCurve();
		      shared_ptr<ParamCurve> c2 = e3->geomEdge()->geomCurve();
		      double t1 = e1_2->prev()->tMin();
		      double t2 = e1_2->prev()->tMax();
		      double t3 = e3->tMin();
		      double t4 = e3->tMax();
		      if (c1.get() != cv1.get() || t1 > tmax1+toptol_.neighbour ||
			  t2 < tmin1-toptol_.neighbour)
			{
			  crvs1.insert(crvs1.begin()+idx1, c1);
			  parmin1.insert(parmin1.begin()+idx1, t1);
			  parmax1.insert(parmax1.begin()+idx1, t2);
			  tmin1 = t1;
			  tmax1 = t2;
			  cv1 = c1;
			}
		      else
			{
			  tmin1 = parmin1[idx1] = std::min(parmin1[idx1], t1);
			  tmax1 = parmax1[idx1] = std::max(parmax1[idx1], t2);
			}

		      if (c2.get() != cv2.get() || t3 > tmax2+toptol_.neighbour ||
			  t4 < tmin1 - toptol_.neighbour)
			{
			  if (e3 == e2_2->prev())
			    {
			      crvs2.insert(crvs2.begin()+idx2, c2);
			      parmin2.insert(parmin2.begin()+idx2, t3);
			      parmax2.insert(parmax2.begin()+idx2, t4);
			    }
			  else
			    {
			      crvs2.push_back(c2);
			      parmin2.push_back(t3);
			      parmax2.push_back(t4);
			      idx2++;
			    }
			  tmin2 = t3;
			  tmax2 = t4;
			  cv2 = c2;
			}
		      else
			{
			  tmin2 = parmin2[idx2] = std::min(parmin2[idx2], t3);
			  tmax2 = parmax2[idx2] = std::max(parmax2[idx2], t4);
			}
		    }
		}
	      e1_2 = e1_2->prev();
	      e2_2 = e3;
	    }

	  // Fetch corresponding boundary curves
	  use_curve = 0;
	  if (crvs1.size() == 1)
	    {
	      cv1 = shared_ptr<ParamCurve>(crvs1[0]->subCurve(parmin1[0], parmax1[0]));
	      if (crvs2.size() > 1)
		use_curve = 1;
	    }
	  if (crvs2.size() == 1)
	    {
	      cv2 = shared_ptr<ParamCurve>(crvs2[0]->subCurve(parmin2[0], parmax2[0]));
	      if (crvs1.size() > 1)
		use_curve = 2;
	    }
	  if (use_curve == 0 && crvs1.size() > 1)
	    {
	      cv1 = shared_ptr<ParamCurve>(crvs1[0]->subCurve(parmin1[0], parmax1[0]));
	      double dist;
	      for (size_t kr=1; kr<crvs1.size(); ++kr)
		{
		  shared_ptr<ParamCurve> tmp = 
		    shared_ptr<ParamCurve>(crvs1[kr]->subCurve(parmin1[kr], parmax1[kr]));
		  Point pt1 = cv1->point(cv1->startparam());
		  Point pt2 = cv1->point(cv1->endparam());
		  Point pt3 = tmp->point(tmp->startparam());
		  Point pt4 = tmp->point(tmp->endparam());
		  if (std::min(pt2.dist(pt3), pt2.dist(pt4)) > 
		      std::min(pt1.dist(pt3), pt1.dist(pt4)))
		    {
		      cv1->reverseParameterDirection();
		      std::swap(pt1, pt2);
		    }
		  if (pt2.dist(pt4) < pt2.dist(pt3))
		    tmp->reverseParameterDirection();
		  //cv1->appendCurve(tmp.get(), 0, dist, false);
		  // try {
		  //   cv1->appendCurve(tmp.get(), 1, dist, true);
		  // }
		  // catch (...)
		  //   {
		      // Perform reparametrization and join with C0 continuity
		      vector<Point> pts1(2);
		      cv1->point(pts1, cv2->endparam(), 1);
		      vector<Point> pts2(2);
		      tmp->point(pts2, tmp->startparam(), 1);
		      double fac = pts2[1].length()/pts1[1].length();

		      // TESTING
		      // fac = 1;
		      // END TESTING

		      double len1 = cv1->estimatedCurveLength();
		      double len2 = tmp->estimatedCurveLength();
		      double s1 = cv1->startparam();
		      double s2 = cv1->endparam();
		      double t1 = tmp->startparam();
		      double t2 = tmp->endparam();
		      fac = len2*(s2-s1)/(len1*(t2-t1));
		      tmp->setParameterInterval(t1, t1+fac*(t2-t1));
#ifdef DEBUG_VOL1
		      std::ofstream ofcv1("cv1_append.g2");
		      cv1->geometryCurve()->writeStandardHeader(ofcv1);
		      cv1->geometryCurve()->write(ofcv1);
		      tmp->geometryCurve()->writeStandardHeader(ofcv1);
		      tmp->geometryCurve()->write(ofcv1);
#endif
		      cv1->appendCurve(tmp.get(), 0, dist, false);
#ifdef DEBUG_VOL1
		      cv1->geometryCurve()->writeStandardHeader(ofcv1);
		      cv1->geometryCurve()->write(ofcv1);
#endif
		      int stop_break = 1.0;
		    // }
		} 
	    }

	  if (use_curve == 0 && crvs2.size() > 1)
	    {
	      cv2 = shared_ptr<ParamCurve>(crvs2[0]->subCurve(parmin2[0], parmax2[0]));
	      double dist;
	      for (size_t kr=1; kr<crvs2.size(); ++kr)
		{
		  shared_ptr<ParamCurve> tmp = 
		    shared_ptr<ParamCurve>(crvs2[kr]->subCurve(parmin2[kr], parmax2[kr]));
		  Point pt1 = cv2->point(cv2->startparam());
		  Point pt2 = cv2->point(cv2->endparam());
		  Point pt3 = tmp->point(tmp->startparam());
		  Point pt4 = tmp->point(tmp->endparam());
		  if (std::min(pt2.dist(pt3), pt2.dist(pt4)) > 
		      std::min(pt1.dist(pt3), pt1.dist(pt4)))
		    {
		      cv2->reverseParameterDirection();
		      std::swap(pt1, pt2);
		    }
		  if (pt2.dist(pt4) < pt2.dist(pt3))
		    tmp->reverseParameterDirection();
		  //cv2->appendCurve(tmp.get(), 0, dist, false);
		  // try {
		  //   cv2->appendCurve(tmp.get(), 1, dist, true);
		  // }
		  // catch (...)
		  //   {
		      // Perform reparametrization and join with C0 continuity
		      vector<Point> pts1(2);
		      cv2->point(pts1, cv2->endparam(), 1);
		      vector<Point> pts2(2);
		      tmp->point(pts2, tmp->startparam(), 1);
		      double fac = pts2[1].length()/pts1[1].length();
		      // TESTING
		      //fac = 1;
		      // END TESTING
		      double len1 = cv2->estimatedCurveLength();
		      double len2 = tmp->estimatedCurveLength();
		      double s1 = cv2->startparam();
		      double s2 = cv2->endparam();
		      double t1 = tmp->startparam();
		      double t2 = tmp->endparam();
		      fac = len2*(s2-s1)/(len1*(t2-t1));
		      tmp->setParameterInterval(t1, t1+fac*(t2-t1));
#ifdef DEBUG_VOL1
		      std::ofstream ofcv2("cv2_append.g2");
		      cv2->geometryCurve()->writeStandardHeader(ofcv2);
		      cv2->geometryCurve()->write(ofcv2);
		      tmp->geometryCurve()->writeStandardHeader(ofcv2);
		      tmp->geometryCurve()->write(ofcv2);
#endif

		      cv2->appendCurve(tmp.get(), 0, dist, false);
#ifdef DEBUG_VOL1
		      cv2->geometryCurve()->writeStandardHeader(ofcv2);
		      cv2->geometryCurve()->write(ofcv2);
#endif
		      int stop_break = 1.0;
		    // }
		} 
	    }

	  // Check orientation
	  bool opposite = false;
	  Point pt1 = cv1->point(cv1->startparam());
	  Point pt2 = cv1->point(cv1->endparam());
	  Point pt3 = cv2->point(cv2->startparam());
	  Point pt4 = cv2->point(cv2->endparam());
	  if (pt1.dist(pt4) + pt2.dist(pt3) < pt1.dist(pt3) + pt2.dist(pt4))
	    opposite = true;

#ifdef DEBUG_VOL1
	  shared_ptr<CurveOnSurface> sfcv1 = 
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv1);
	  shared_ptr<CurveOnSurface> sfcv2 = 
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv2);
          if (sfcv1.get() && sfcv2.get()) {
	    if (sfcv1->spaceCurve().get())
	      {
		sfcv1->spaceCurve()->writeStandardHeader(of);
		sfcv1->spaceCurve()->write(of);
	      }
	    if (sfcv2->spaceCurve().get())
	      {
		sfcv2->spaceCurve()->writeStandardHeader(of);
		sfcv2->spaceCurve()->write(of);
	      }
          }
#endif

	  // Store curves. Mind the orientation
	  if (use_curve != 2 && turned)
	    cv1->reverseParameterDirection();
	  if ((turned && !opposite) ||
	      (!turned && opposite))
	    {
	      if (use_curve != 1)
		cv2->reverseParameterDirection();
	      turned = true;
	    }
	  else
	    turned = false;

	  // In a degenerate situation, check if the obtained curve
	  // is correct or if the transmission should be passed by 
	  // a surface degenating to a point
	  if (nopt_dir >= 0 && pardir != nopt_dir)
	    {
	      // Fetch endpoints of curve
	      Point pos1, pos2;
	      ParamCurve *curr_cv = (use_curve == 1 || use_curve > 2) ?
		cv1.get() : cv2.get();
	      pos1 = curr_cv->point(curr_cv->startparam());
	      pos2 = curr_cv->point(curr_cv->endparam());
	      if (pos1.dist(deg_pt) < deg_tol || 
		  pos2.dist(deg_pt) < deg_tol)
		use_curve = -1;
	    }
	  if (use_curve < 0)
	    bd_cvs[ki] = make_pair(dummy, dummy);
	  else if (use_curve == 0)
	    bd_cvs[ki] = make_pair(cv1, cv2);
	  else if (use_curve == 1)
	    bd_cvs[ki] = make_pair(cv1, dummy);
	  else
	    bd_cvs[ki] = make_pair(dummy, cv2);

	  // Update
	  if (!deg2)
	    {
	      ix1 = ix2;
	      sf1 = sf2;
	    }
	}

      curves.push_back(bd_cvs);
      indices.push_back(sf_idx);
    }
      // Ensure consistent orientation of all curves
#ifdef DEBUG_VOL1
  std::ofstream of3("co1_cvs0.g2");
  for (pardir=0; pardir<3; ++pardir)
    {
      for (size_t kh=0; kh<curves[pardir].size(); ++kh)
	{
	  if (curves[pardir][kh].first.get())
	    {
	      shared_ptr<SplineCurve> tmp =
		shared_ptr<SplineCurve>(curves[pardir][kh].first->geometryCurve());
	      tmp->writeStandardHeader(of3);
	      tmp->write(of3);

	      Point tmp_pnt = tmp->ParamCurve::point(tmp->startparam());
	      of3 << "400 1 0 4 255 0 0 255" << std::endl;
	      of3 << "1" << std::endl;
	      of3 << tmp_pnt << std::endl;
	    }
	  else if (curves[pardir][kh].second.get())
	    {
	      shared_ptr<SplineCurve> tmp = 
		shared_ptr<SplineCurve>(curves[pardir][kh].second->geometryCurve());
	      tmp->writeStandardHeader(of3);
	      tmp->write(of3);

	      Point tmp_pnt = tmp->ParamCurve::point(tmp->startparam());
	      of3 << "400 1 0 4 255 0 0 255" << std::endl;
	      of3 << "1" << std::endl;
	      of3 << tmp_pnt << std::endl;
	    }
	}
    }
#endif

  for (pardir=0; pardir<3; ++pardir)
    {
      for (kj=1; kj<4; ++kj)
	{
	  // Evaluate start points of curves
	  Point pos1, pos2, pos3, pos4;
	  if (curves[pardir][kj-1].first)
	    {
	      pos1 = curves[pardir][kj-1].first->point(curves[pardir][kj-1].first->startparam());
	      pos2 = curves[pardir][kj-1].first->point(curves[pardir][kj-1].first->endparam());
	    }
	  else if (curves[pardir][kj-1].second)
	    {
	      pos1 = curves[pardir][kj-1].second->point(curves[pardir][kj-1].second->startparam());
	      pos2 = curves[pardir][kj-1].second->point(curves[pardir][kj-1].second->endparam());
	    }
	  else
	    continue;

	  if (curves[pardir][kj].first)
	    {
	      pos3 = curves[pardir][kj].first->point(curves[pardir][kj].first->startparam());
	      pos4 = curves[pardir][kj].first->point(curves[pardir][kj].first->endparam());
	    }
	  else if (curves[pardir][kj].second)
	    {
	      pos3 = curves[pardir][kj].second->point(curves[pardir][kj].second->startparam());
	      pos4 = curves[pardir][kj].second->point(curves[pardir][kj].second->endparam());
	    }
	  else
	    continue;

	  // Check if any curve in the other parameter directions join
	  // the two start points
	  double d1=0.0, d2=0.0, d3=0.0, d4=0.0, d5=0.0, d6=0.0;
	  for (ki=0; ki<3; ++ki)
	    {
	      if (ki == pardir)
		continue;
	      int kr;
	      for (kr=0; kr<4; kr++)
		{
		  Point pos5, pos6;
		  if (curves[ki][kr].first)
		    {
		      pos5 = curves[ki][kr].first->point(curves[ki][kr].first->startparam());
		      pos6 = curves[ki][kr].first->point(curves[ki][kr].first->endparam());
		    }
		  else if (curves[ki][kr].second)
		    {
		      pos5 = curves[ki][kr].second->point(curves[ki][kr].second->startparam());
		      pos6 = curves[ki][kr].second->point(curves[ki][kr].second->endparam());
		    }
		  else 
		    continue;
		  
		  d1 = std::min(pos5.dist(pos1), pos5.dist(pos2));
		  d2 = std::min(pos5.dist(pos3), pos5.dist(pos4));
		  d3 = std::min(pos6.dist(pos1), pos6.dist(pos2));
		  d4 = std::min(pos6.dist(pos3), pos6.dist(pos4));
		  d5 = std::min(pos1.dist(pos5), pos1.dist(pos6));
		  d6 = std::min(pos3.dist(pos5), pos3.dist(pos6));
		  if (std::min(d1,d2) < toptol_.neighbour &&
		      std::min(d3,d4) < toptol_.neighbour)
		    break;
		}
	      if (kr < 4)
		break;
	    }

	  if ((d5 >= toptol_.neighbour && d6 < toptol_.neighbour) ||
	      (d6 >= toptol_.neighbour && d5 < toptol_.neighbour))
	    {
	      // Opposite direction
	      if (curves[pardir][kj].first)
		curves[pardir][kj].first->reverseParameterDirection();
	      if (curves[pardir][kj].second)
		curves[pardir][kj].second->reverseParameterDirection();
	    }

#ifdef DEBUG_VOL1
	  std::ofstream of2("co1_cvs.g2");
	  for (int kr=0; kr<4; ++kr)
	    {
	      if (curves[pardir][kr].first.get())
		{
		  shared_ptr<SplineCurve> tmp =
		    shared_ptr<SplineCurve>(curves[pardir][kr].first->geometryCurve());
		  tmp->writeStandardHeader(of2);
		  tmp->write(of2);

		  Point tmp_pnt = tmp->ParamCurve::point(tmp->startparam());
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << tmp_pnt << std::endl;
		}
	      else if (curves[pardir][kr].second.get())
		{
		  shared_ptr<SplineCurve> tmp = 
		    shared_ptr<SplineCurve>(curves[pardir][kr].second->geometryCurve());
		  tmp->writeStandardHeader(of2);
		  tmp->write(of2);

		  Point tmp_pnt = tmp->ParamCurve::point(tmp->startparam());
		  of2 << "400 1 0 4 255 0 0 255" << std::endl;
		  of2 << "1" << std::endl;
		  of2 << tmp_pnt << std::endl;
		}
	    }
#endif
	}
      int stop_break;
      stop_break = 1;
    }

  return true;
}

//===========================================================================
// 
// 
vector<shared_ptr<ftSurface> > 
ftVolume::generateMissingBdSurf(int degree,
				vector<pair<Point,Point> >& corr_vx_pts,
				bool perform_step2, bool smooth_connections,
				bool& trimmed, int max_nmb)
//===========================================================================
{
  trimmed = false;

  // Get loops corresponding to the missing surfaces
  vector<shared_ptr<ftSurface> > faces;
  vector<vector<ftEdge*> > sf_loops = getMissingSfLoops(corr_vx_pts, 
							perform_step2,
							smooth_connections,
							max_nmb);

  size_t ki, kj, kr, kh;  
  size_t nmb_loops = sf_loops.size();
  
  if (nmb_loops == 0)
    {
      // No potential for generating inner boundary surfaces.
      // Try to merge boundary surfaces of the outer shell 
      simplifyOuterBdShell(degree);
      return faces;  // Return the empty face vector
    }

#ifdef DEBUG_VOL1
  std::ofstream of("missing_surfaces.g2");
#endif

  if (sf_loops.size() > 1)
    organizeLoops(sf_loops);

  
  for (ki=0; ki<sf_loops.size(); ++ki)
    {
      // Make a pair of missing surfaces. Make sure to update twin pointers
      shared_ptr<ftSurface> face1, face2;
      vector<pair<ftEdge*, ftEdge*> > replaced_wires;

      // Check if a degenerate point is expected and can be localized
      Vertex *deg_vx = NULL;
      if (sf_loops[ki].size() == 3)
	{
	  for (kj=0; kj<sf_loops.size(); ++kj)
	    {
	      if (ki == kj)
		continue;
	      if (sf_loops[kj].size() == 3)
		{
		  // Two loops that will lead to degenerate surfaces. Check if they coincides
		  // in one point
		  std::set<shared_ptr<Vertex> > common_vxs;
		  for (kr=0; kr<sf_loops[ki].size(); ++kr)
		    for (kh=0; kh<sf_loops[kj].size(); ++kh)
		      {
			shared_ptr<Vertex> vx = 
			  sf_loops[ki][kr]->getCommonVertex(sf_loops[kj][kh]);
			if (vx.get())
			  common_vxs.insert(common_vxs.end(), vx);
		      }
		  if (common_vxs.size() == 1)
		    deg_vx = (*common_vxs.begin()).get();
		}
	    }

	  if (!deg_vx)
	    {
	      // Check if any surface connected to loop vertex is triangular
	      std::set<shared_ptr<Vertex> > loop_vxs;
	      for (kj=0; kj<sf_loops[ki].size(); ++kj)
		{
		  shared_ptr<Vertex> v1, v2;
		  sf_loops[ki][kj]->getVertices(v1,v2);
		  loop_vxs.insert(v1);
		  loop_vxs.insert(v2);
		}

	      for (auto it=loop_vxs.begin(); it != loop_vxs.end(); ++it)
		{
		  vector<ftSurface*> adj_faces = (*it)->faces(this);
		  for (kj=0; kj<adj_faces.size(); ++kj)
		    {
		      vector<shared_ptr<Vertex> > corner =
			adj_faces[kj]->getCornerVertices(toptol_.bend);
		      if (corner.size() == 3)
			break;
		    }
		  if (kj < adj_faces.size())
		    {
		      // Check if the vertex is shared with other loops
		      for (kr=0; kr<sf_loops.size(); ++kr)
			{
			  if (kr == ki)
			    continue;
			  for (kh=0; kh<sf_loops[kr].size(); ++kh)
			    {
			      if (sf_loops[kr][kh]->hasVertex((*it).get()))
				break;
			    }
			  if (kh < sf_loops[kr].size())
			    break;
			}
		      if (kr == sf_loops.size())
			{
			  deg_vx = (*it).get();
			  break;
			}
		    }
		}
	    }
	}
      makeSurfacePair(sf_loops[ki], degree, face1, face2, replaced_wires, deg_vx);

      // Check if any surfaces are generated
      if (face1.get())
	{
	  // Update the remaining part of the outer shell with respect
	  // to the new surfaces in case other surfaces are intersected
	  // by them
	  bool modified = 
	    ftVolumeTools::updateWithSplitFaces(getShell(0), face1, face2,
						replaced_wires);
	  if (modified)
	    trimmed = true;
	  
	  // It should not be necessary to split edges in the loop to make a
	  // connection
	  face1->connectTwin(face2.get(), toptol_.neighbour, false);

#ifdef DEBUG_VOL1
	  face1->surface()->writeStandardHeader(of);
	  face1->surface()->write(of);
	  face2->surface()->writeStandardHeader(of);
	  face2->surface()->write(of);
#endif

	  faces.push_back(face1);
	  faces.push_back(face2);
	  for (kj=0; kj<replaced_wires.size(); ++kj)
	    for (kr=ki+1; kr<sf_loops.size(); ++kr)
	      for (kh=0; kh<sf_loops[kr].size(); ++kh)
		if (sf_loops[kr][kh] == replaced_wires[kj].first)
		  sf_loops[kr][kh] = replaced_wires[kj].second;
	}
    }


  // Clean up in intermediate edges to define missing faces
  if (missing_edges_.size() > 0)
    eraseMissingEdges();

  if (faces.size() == 0)
    {
      // No inner boundary surfaces created.
      // Try to merge boundary surfaces of the outer shell 
      simplifyOuterBdShell(degree);
      return faces;  // Return the empty face vector
    }
  return faces;
}

//===========================================================================
// 
// 
void ftVolume::organizeLoops(vector<vector<ftEdge*> >& loops)
//===========================================================================
{
  int ki, kj;

  // Check for obsolete loops
  for (ki=0; ki<(int)loops.size(); )
    {
      for (kj=ki+1; kj<(int)loops.size(); ++kj)
	{
	  // Identify edges not belonging to both loops
	  size_t kh, kr;
	  vector<ftEdge*> single_edg1, single_edg2;
	  for (kh=0; kh<loops[ki].size(); ++kh)
	    {
	      ftEdge *e1 = loops[ki][kh];
	      for (kr=0; kr<loops[kj].size(); ++kr)
		{
		  if (e1 == loops[kj][kr] || e1 == loops[kj][kr]->twin() ||
		      e1->twin() == loops[kj][kr] ||
		      e1->twin() == loops[kj][kr]->twin())
		    break;
		}
	      if (kr == loops[kj].size())
		single_edg1.push_back(e1);
	    }

	  for (kh=0; kh<loops[kj].size(); ++kh)
	    {
	      ftEdge *e1 = loops[kj][kh];
	      for (kr=0; kr<loops[ki].size(); ++kr)
		{
		  if (e1 == loops[ki][kr] || e1 == loops[ki][kr]->twin() ||
		      e1->twin() == loops[ki][kr] ||
		      e1->twin() == loops[ki][kr]->twin())
		    break;
		}
	      if (kr == loops[ki].size())
		single_edg2.push_back(e1);
	    }

	  // Check if all single edges belong to the same face
	  vector<ftEdge*> edgs(single_edg1.begin(), single_edg1.end());
	  edgs.insert(edgs.end(), single_edg2.begin(), single_edg2.end());
	  ftSurface *common_face = NULL;
	  if (edgs.size() > 0)
	    common_face = Path::identifyCommonFace(edgs);

	  if (common_face)
	    {
	      // Remove the largest loop
	      if (single_edg1.size() > single_edg2.size())
		{
		  loops.erase(loops.begin()+ki);
		  break;
		}
	      else
		{
		  loops.erase(loops.begin()+kj);
		  --kj;
		  break;
		}
	    }
	  int stop_break = 1;
	}
      if (kj == loops.size())
	++ki;
	
    }
  
  // Make sure that loops consisting of only added edges treated
  // at last
  size_t nmb_loops = loops.size();
  for (ki=0; ki<nmb_loops;)
    {
      for (kj=0; kj<loops[ki].size(); ++kj)
	if (loops[ki][kj]->twin())
	  break;
      if (kj == loops[ki].size())
	{
	  loops.push_back(loops[ki]);
	  loops.erase(loops.begin()+ki);
	  nmb_loops--;
	}
      else 
	ki++;
    }
}


//===========================================================================
// 
// 
void ftVolume::makeSurfacePair(vector<ftEdge*>& loop,
			       int degree,
			       shared_ptr<ftSurface>& face1,
			       shared_ptr<ftSurface>& face2,
			       vector<pair<ftEdge*,ftEdge*> >& replaced_wires,
			       Vertex *deg_vx)
//===========================================================================
{
  size_t ki;
#ifdef DEBUG_VOL1
  std::ofstream of_in("curr_loop.g2");
  for (ki=0; ki<loop.size(); ++ki)
    {
      ftEdge *e1 = loop[ki];
      shared_ptr<ParamCurve> cv =
	shared_ptr<ParamCurve>(e1->geomCurve()->subCurve(e1->tMin(), e1->tMax()));
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
      if (sfcv.get())
	{
	  if (sfcv->spaceCurve().get())
	    {
	      sfcv->spaceCurve()->writeStandardHeader(of_in);
	      sfcv->spaceCurve()->write(of_in);
	    }
	}
      else
	{
	  cv->writeStandardHeader(of_in);
	  cv->write(of_in);
	}
    }
#endif

  // Check loop for multiple joints at endpoints
  double neighbour = toptol_.neighbour;
  Point pos1 = loop[loop.size()-1]->point(loop[loop.size()-1]->tMin());
  Point pos2 = loop[loop.size()-1]->point(loop[loop.size()-1]->tMax());
  for (ki=0; ki<loop.size(); ++ki)
    {
      Point pos3 = loop[ki]->point(loop[ki]->tMin());
      Point pos4 = loop[ki]->point(loop[ki]->tMax());
      if (pos1.dist(pos2) > neighbour && pos3.dist(pos4) > neighbour && 
	  ((pos1.dist(pos4) < neighbour && pos2.dist(pos3) < neighbour) ||
	   (pos1.dist(pos3) < neighbour && pos2.dist(pos4) < neighbour)))
	return;  // Loop with multiple joints found
      pos1 = pos3;
      pos2 = pos4;
    }

  // Change start position of the loops to standard edges
  size_t idx1 = 0;
  while (!loop[0]->twin())
    {
      loop.push_back(loop[0]);
      loop.erase(loop.begin());
      idx1++;
      if (idx1 >= loop.size())
	break;
    }

  ASSERT(idx1 < loop.size());

  // Sort edges to be joined before surface construction
  // Fetch the curves associated with each boundary of the surfaces
  // to be constructed
  //size_t ki, kr;
  vector<shared_ptr<ParamCurve> > space_cvs;
  vector<Point> joint_points;
  getEdgeCurves(loop, space_cvs, joint_points);

#ifdef DEBUG_VOL1
  std::ofstream sc("spacecrvs_space.g2");
  for (ki=0; ki<space_cvs.size(); ++ki)
    {
      shared_ptr<SplineCurve> tmp_space1 = 
	shared_ptr<SplineCurve>(space_cvs[ki]->geometryCurve());
      tmp_space1->writeStandardHeader(sc);
      tmp_space1->write(sc);
    }
#endif

   // Check for degenerate configurations
   if (space_cvs.size() == 2)
     {
       // Surface degenerates to a line. Insert room for degenerate boundary curves
       shared_ptr<ParamCurve> dummy;
       space_cvs.insert(space_cvs.begin()+1, dummy);
       space_cvs.insert(space_cvs.end(), dummy);
     }
   else if (space_cvs.size() == 3)
     {
       // Triangular surface. Find position of degenerate edge
       shared_ptr<ParamCurve> dummy;
       if (deg_vx)
	 {
	   Point vx_pt = deg_vx->getVertexPoint();
	   double dist[3];
	   for (ki=0; ki<space_cvs.size(); ++ki)
	     {
	       dist[ki] = std::min(vx_pt.dist(space_cvs[ki]->point(space_cvs[ki]->startparam())),
				   vx_pt.dist(space_cvs[ki]->point(space_cvs[ki]->endparam())));
	     }
	   int min_ix = 3;
	   if (dist[0]+dist[1] < dist[1]+dist[2] &&
	       dist[0]+dist[1] < dist[0]+dist[2])
	     min_ix = 1;
	   else if (dist[1]+dist[2] < dist[0]+dist[1] &&
		    dist[1]+dist[1] < dist[0]+dist[2])
	     min_ix = 2;
	   space_cvs.insert(space_cvs.begin()+min_ix, dummy);
	 }
       else
	 {
	   double len[3];
	   for (int ka=0; ka<3; ++ka)
	     len[ka] = space_cvs[ka]->estimatedCurveLength();
	   int ix = 0;
	   if (len[1] < len[ix])
	     ix = 1;
	   if (len[2] < len[ix])
	     ix = 2;
	   int ix2 = (ix == 1) ? 3 : ((ix == 0) ? 2 : 1);
	   space_cvs.insert(space_cvs.begin()+ix2, dummy);
	 }
     }

  // Ensure that the orientation of the curves is suitable for Coons patch
   Point deg_dummy;
   vector<shared_ptr<ParamCurve> > dummy_cvs(space_cvs.size());
   (void)sortCoonsPatchBdCvs(space_cvs, dummy_cvs, deg_dummy, toptol_.gap);

  // Create surface as Coons patch
   double knot_diff_tol = 1.0e-4; //1.0e-5;
  vector<shared_ptr<ParamCurve> > space_cvs2(space_cvs.size());
  for (ki=0; ki<space_cvs.size(); ++ki)
    {
      space_cvs2[ki] = shared_ptr<ParamCurve>(space_cvs[ki]->geometryCurve());
      shared_ptr<SplineCurve> tmp_spline = 
	dynamic_pointer_cast<SplineCurve,ParamCurve>(space_cvs2[ki]);
      bool found_indistinct_knots = false;
      if (tmp_spline.get())
	{
	  vector<double> indistinct_knots;
	  found_indistinct_knots = tmp_spline->basis().indistinctKnots(knot_diff_tol,
									indistinct_knots);
	}
      if (space_cvs2[ki].get() == 0 || tmp_spline->rational() || found_indistinct_knots)
	{
	  // Approximate by non-rational spline curve
	  int nmb_sample = 5;
	  vector<shared_ptr<SplineCurve> > tmp_cvs = 
	    CurveCreators::curveApprox(&space_cvs[ki], 1, toptol_.gap,
				       degree, nmb_sample);
	  space_cvs2[ki] = tmp_cvs[0];
	}
    }

#ifdef DEBUG_VOL1
  std::ofstream sc2("spacecrvs2_space.g2");
  for (ki=0; ki<space_cvs2.size(); ++ki)
    {
      shared_ptr<SplineCurve> tmp_space1 = 
	shared_ptr<SplineCurve>(space_cvs2[ki]->geometryCurve());
      tmp_space1->writeStandardHeader(sc2);
      tmp_space1->write(sc2);
    }
#endif

  double ang_fac = 3.0; // 2.0; Maybe some other way of deciding is required
  // We have entered this function since there is an expectation of defining
  // a split surface, but is it a possibility of merging wrong curves?
  if (space_cvs2.size() != 4)
    {
      while (space_cvs2.size() > 4)
	{
	  int min_ix = -1;
	  double min_ang = std::numeric_limits<double>::max();
	  vector<Point> pts1 = 
	    space_cvs2[space_cvs2.size()-1]->point(space_cvs2[space_cvs2.size()-1]->endparam(), 1);
	  for (ki=0; ki<space_cvs2.size(); ++ki)
	    {
	      vector<Point> pts2 = 
		space_cvs2[ki]->point(space_cvs2[ki]->startparam(), 1);
	      double ang = pts1[1].angle(pts2[1]);
	      if (ang < min_ang)
		{
		  min_ang = ang;
		  min_ix = (int)ki;
		}
	      pts1 = space_cvs2[ki]->point(space_cvs2[ki]->endparam(), 1);
	    }
	  if (min_ang < ang_fac*toptol_.bend)
	    {
	      Point joint = space_cvs2[min_ix]->point(space_cvs2[min_ix]->startparam());
	      int ix2 = (min_ix == 0) ? (int)(space_cvs2.size()-1) : min_ix-1;
	      double dist;
	      space_cvs2[ix2]->appendCurve(space_cvs2[min_ix].get(), 1, dist);
	      space_cvs2.erase(space_cvs2.begin()+min_ix);
	      joint_points.push_back(joint);
	    }
	  else
	    break;
	}
      if (space_cvs2.size() != 4)
	{
#ifdef DEBUG_VOL1
	  MESSAGE("Not a valid splitting surface configuration");
#endif
	  return;  // Not a valid splitting surface
	}
    }
  CurveLoop bd_loop(space_cvs2, toptol_.gap);

  shared_ptr<SplineSurface> sf1 = 
    shared_ptr<SplineSurface>(CoonsPatchGen::createCoonsPatch(bd_loop));

  // Create seond surface. 
  shared_ptr<SplineSurface> sf2(sf1->clone());
  sf2->turnOrientation();

  // Make surface on volume surfaces
  shared_ptr<ParamSurface> dummy;
  shared_ptr<ParamSurface> psf1 = sf1;
  shared_ptr<ParamSurface> psf2 = sf2;
  shared_ptr<SurfaceOnVolume> vol_sf1 =
    shared_ptr<SurfaceOnVolume>(new SurfaceOnVolume(vol_, dummy, psf1, false));
  shared_ptr<SurfaceOnVolume> vol_sf2 =
    shared_ptr<SurfaceOnVolume>(new SurfaceOnVolume(vol_, dummy, psf2, false));
  vol_sf1->setCreationHistory(1);  // Mark as internal splitting surface
  vol_sf2->setCreationHistory(1);


#ifdef DEBUG_VOL1
  std::ofstream of("bd_sf_vol.g2");
  vol_sf1->spaceSurface()->writeStandardHeader(of);
  vol_sf1->spaceSurface()->write(of);
  vol_sf2->spaceSurface()->writeStandardHeader(of);
  vol_sf2->spaceSurface()->write(of);
#endif

//   // Check if these surfaces lie within the body
//   // First pick a set of point, one is not sufficient
//   RectDomain dom = vol_sf1->containingDomain();
//   double u1 = dom.umin();
//   double u2 = dom.umax();
//   double v1 = dom.vmin();
//   double v2 = dom.vmax();
//   double udel = (u2-u1)/(double)2;
//   double vdel = (v2-v1)/(double)2;
//   double upar, vpar;
//   int ka1, ka2;
//   for (ka2=0, vpar=v1+0.5*vdel; ka2<2; ++ka2, vpar+=vdel)
//     {
//       for (ka1=0, upar=u1+0.5*udel; ka1<2; ++ka1, upar+=udel)
// 	{
// 	  Point pnt_in_sf = vol_sf1->ParamSurface::point(upar, vpar);
// 	  bool inside_bd = isInside(pnt_in_sf);
// 	  if (inside_bd)
// 	    break;
// 	}
//       if (ka1 < 2)
// 	break;
//     }
//   if (ka2 == 2)
//     {
// #ifdef DEBUG_VOL1
//       std::cout << "Split surface not in body" <<std::endl;
// #endif
//       return;  // Not a valid splitting surface. 
//     }

  double max_loop_dist = 0.0;
  for (ki=0; ki<loop.size(); ++ki)
    {
      size_t kj = (ki<loop.size()-1) ? ki+1 : 0;
      Point pos1 = loop[ki]->point(loop[ki]->tMin());
      Point pos2 = loop[ki]->point(loop[ki]->tMax());
      Point pos3 = loop[kj]->point(loop[kj]->tMin());
      Point pos4 = loop[kj]->point(loop[kj]->tMax());
      double loop_dist = std::min(std::min(pos1.dist(pos3),
					   pos1.dist(pos4)),
				  std::min(pos2.dist(pos3),
					   pos2.dist(pos4)));
      max_loop_dist = std::max(max_loop_dist, loop_dist);
    }

  // Define topology
  // First make faces
  face1 = shared_ptr<ftSurface>(new ftSurface(vol_sf1, -1));
  (void)face1->createInitialEdges(/*max_loop_dist < toptol_.gap ?
				    toptol_.gap :*/ toptol_.neighbour,
				    toptol_.bend, true);
  shared_ptr<Loop> loop1 = face1->getBoundaryLoop(0);

  // Split the loop at joint points. Find edge index and parameter by
  // a closest point computation
  for (ki=0; ki<joint_points.size(); ++ki)
    {
      int clo_ind;
      double clo_par, clo_dist;
      Point clo_pt;
      loop1->closestPoint(joint_points[ki], clo_ind, clo_par, clo_pt,
			  clo_dist);
      loop1->split(clo_ind, clo_par);
    }

  face2 = shared_ptr<ftSurface>(new ftSurface(vol_sf2, -1));
  (void)face2->createInitialEdges(/*max_loop_dist < toptol_.gap ?
				    toptol_.gap : */toptol_.neighbour,
				    toptol_.bend, true);
  shared_ptr<Loop> loop2 = face2->getBoundaryLoop(0);

  // Split the loop at joint points. Find edge index and parameter by
  // a closest point computation
  for (ki=0; ki<joint_points.size(); ++ki)
    {
      int clo_ind;
      double clo_par, clo_dist;
      Point clo_pt;
      loop2->closestPoint(joint_points[ki], clo_ind, clo_par, clo_pt,
			  clo_dist);
      try {
	loop2->split(clo_ind, clo_par);
      }
      catch (...)
	{
	  THROW("Joint point split failed");
	}
    }

  // Check input and swith edges in the input loop if necessary
  ftSurface *f1 = loop[0]->face()->asFtSurface();
  Body *bd = f1->getBody();
  for (ki=0; ki<loop.size(); ++ki)
    {
      shared_ptr<EdgeVertex> radial = loop[ki]->getEdgeMultiplicityInstance();
      if (radial.get() && bd && radial->nmbUniqueEdges(bd) > 1)
	{
	  ftEdge *tmp_edge = getLeftLoopEdge(face1.get(), bd, radial);
	  if (tmp_edge)
	    loop[ki] = tmp_edge;
	}
    }

  // Insert new faces into the adjacency relationship represented by
  // the edges. 
  // First remove trim information where the faces will be added
  // Remember start edges
#ifdef DEBUG_VOL1
  std::ofstream ofin1("intwin1.g2");
  std::ofstream ofin2("intwin2.g2");
#endif
  ftEdge* e2 = loop[0]->twin() ? loop[0]->twin()->geomEdge() : NULL;
  vector<ftEdge*> twins;
  for (ki=0; ki<loop.size(); ++ki)
    {
      if (loop[ki]->twin())
	{
	  twins.push_back(loop[ki]->twin() ? loop[ki]->twin()->geomEdge() : NULL);
	  loop[ki]->disconnectTwin();
	}
      else 
	{
	  shared_ptr<EdgeVertex> radial = loop[ki]->getEdgeMultiplicityInstance();
	  vector<ftEdge*> rad_edgs; 
	  if (radial.get())
	    rad_edgs = radial->uniqueEdges();
	  if (rad_edgs.size() > 1)
	    twins.push_back((rad_edgs[0] == loop[ki]) ? rad_edgs[1] : 
			    rad_edgs[0]);
	  else
	    twins.push_back(NULL);
	}

#ifdef DEBUG_VOL1
      if (loop[ki]->face())
	{
	  loop[ki]->face()->asFtSurface()->surface()->writeStandardHeader(ofin1);
	  loop[ki]->face()->asFtSurface()->surface()->write(ofin1);  
	}
      if (twins[ki])
	{
	  twins[ki]->face()->asFtSurface()->surface()->writeStandardHeader(ofin2);
	  twins[ki]->face()->asFtSurface()->surface()->write(ofin2);    
	}
#endif
    }

  // Add new twin configurations. First surface
  // Check loop orientations
#ifdef DEBUG_VOL1
  std::ofstream of1("twin1.g2");
  face1->surface()->writeStandardHeader(of1);
  face1->surface()->write(of1);
#endif

  // Find start position and orientation
  int idx2 = -1;
  int dir1 = 1;
  double dist = 1.0e8;
  int status = 0;
  ftEdge *l1;
  double t1, t2, d1, d2, seed;
  Point pt0, pt1, pt2;
  int kj;
  int nmb = (int)loop1->size();
  for (ki=0; ki<loop1->size(); ++ki)
    {
      l1 = loop1->getEdge(ki)->geomEdge();
      t1 = 0.5*(l1->tMin()+l1->tMax());
      pt1 = l1->point(t1);
      if (e2)
	{
	  seed = 0.5*(e2->tMin() + e2->tMax());
	  e2->closestPoint(pt1, t2, pt2, d2, &seed);
	  if (d2 < dist)
	    {
	      dist = d2;
	      idx2 = (int)ki;
	    }
	}
    }
  pt0 = loop[1]->point(0.5*(loop[1]->tMin()+loop[1]->tMax()));
  kj = (idx2 + 1) % nmb;
  l1 = loop1->getEdge(kj)->geomEdge();
  seed = 0.5*(l1->tMin()+l1->tMax());
  l1->closestPoint(pt0, t2, pt1, d1, &seed);

  kj = idx2 - 1;
  if (kj < 0)
    kj = nmb - 1;
  l1 = loop1->getEdge(kj)->geomEdge();
  seed = 0.5*(l1->tMin()+l1->tMax());
  l1->closestPoint(pt0, t2, pt2, d2, &seed);
  if (d2 < d1)
    dir1 = -1;

  Point tan1, tan2;
  for (ki=0, kj=idx2; ki<(size_t)nmb; ++ki, kj+=dir1)
    {
      if (kj < 0)
	kj = nmb-1;
      if (kj >= nmb)
	kj = 0;
      ftEdge *c2 = loop[ki];
      l1 = loop1->getEdge(kj)->geomEdge();
      t1 = 0.5*(l1->tMin()+l1->tMax());
      pt1 = l1->point(t1);
      seed = 0.5*(c2->tMin() + c2->tMax());
      c2->closestPoint(pt1, t2, pt2, d2, &seed);
      tan1 = l1->tangent(t1);
      tan2 = c2->tangent(t2);
      if (tan1*tan2 >= 0.0 /*&& twins[ki]*/)
	{
	  c2 = twins[ki];
	  std::swap(twins[ki], loop[ki]);
	}
      if (c2 && c2->face())
	l1->connectTwin(c2, status);
      else if (!(twins[ki] && twins[ki]->face()))
	replaced_wires.push_back(make_pair((loop[ki]) ? loop[ki] : twins[ki] ,l1));

#ifdef DEBUG_VOL1
      if (c2 && c2->face())
	{
	  c2->face()->asFtSurface()->surface()->writeStandardHeader(of1);
	  c2->face()->asFtSurface()->surface()->write(of1);
	}
#endif
    }

  // Second surface
#ifdef DEBUG_VOL1
  std::ofstream of2("twin2.g2");
#endif

  // Start position
  dist = 1.0e8;
  e2 = twins[0];
  idx2 = -1;
  for (ki=0; ki<loop2->size(); ++ki)
    {
      l1 = loop2->getEdge(ki)->geomEdge();
      t1 = 0.5*(l1->tMin()+l1->tMax());
      pt1 = l1->point(t1);
      if (e2)
	{
	  seed = 0.5*(e2->tMin() + e2->tMax());
	  e2->closestPoint(pt1, t2, pt2, d2, &seed);
	  if (d2 < dist)
	    {
	      dist = d2;
	      idx2 = (int)ki;
	    }
	}
    }
  
  dir1 *= -1;
  for (ki=0, kj=idx2; ki<loop2->size(); ++ki, kj+=dir1)
    {
      if (kj < 0)
	kj = nmb-1;
      if (kj >= nmb)
	kj = 0;
      ftEdge *c2 = twins[ki];
      //l1 = loop2->getEdge(nmb-1-(ki+idx2)%nmb)->geomEdge();
      l1 = loop2->getEdge(kj)->geomEdge();
       t1 = 0.5*(l1->tMin()+l1->tMax());
      pt1 = l1->point(t1);
      if (c2)
	{
	  seed = 0.5*(c2->tMin() + c2->tMax());
	  c2->closestPoint(pt1, t2, pt2, d2, &seed);
	}
//       tan1 = l1->tangent(t1);
//       tan2 = c2->tangent(t2);
//       if (tan1*tan2 >= 0.0)
// 	c2 = twins[ki];
      Point pt1 = l1->point(l1->tMin());
      Point pt2 = l1->point(l1->tMax());
      Point pt3, pt4;
      if (c2)
	{
	  pt3 = c2->point(c2->tMin());
	  pt4 = c2->point(c2->tMax());
	}
      if (c2 && c2->face())
	l1->connectTwin(c2, status);

#ifdef DEBUG_VOL1
      if (c2 && c2->face())
	{
	  c2->face()->asFtSurface()->surface()->writeStandardHeader(of2);
	  c2->face()->asFtSurface()->surface()->write(of2);
	}
#endif
    }

  face1->setBody(this);
  face2->setBody(this);
  // // It should not be necessary to split edges in the loop to make a
  // // connection
  // face1->connectTwin(face2.get(), toptol_.neighbour, false);

#ifdef DEBUG_VOL1
  std::ofstream of3("twin1_2.g2");
  face1->surface()->writeStandardHeader(of3);
  face1->surface()->write(of3);
  std::ofstream of4("twin2_2.g2");
  face2->surface()->writeStandardHeader(of4);
  face2->surface()->write(of4);
  vector<ftSurface*> adjf1, adjf2;
  face1->getAdjacentFaces(adjf1);
  face2->getAdjacentFaces(adjf2);
  for (ki=0; ki<adjf1.size(); ++ki)
    {
      adjf1[ki]->surface()->writeStandardHeader(of3);
      adjf1[ki]->surface()->write(of3);
    }
  for (ki=0; ki<adjf2.size(); ++ki)
    {
      adjf2[ki]->surface()->writeStandardHeader(of4);
      adjf2[ki]->surface()->write(of4);
    }
#endif
  int stop_break;
  stop_break = 1;

}

//===========================================================================
void ftVolume::getEdgeCurves(vector<ftEdge*>& loop, 
			     vector<shared_ptr<ParamCurve> >& space_cvs,
			     vector<Point>& joint_points)
//===========================================================================
{
  // Sort edges to be joined before surface construction
  // It is not expected that the last and the first edge should be joined
  size_t ki, kj, kr;
  vector<vector<ftEdge*> > joined_loop;
  ftEdge *prev = 0;
  for (ki=0; ki<loop.size(); ++ki)
    {
      double ang = M_PI;
      ftEdge *curr = loop[ki];
      if (ki > 0)
	{
	  // Check if there is a corner or T-joint between this curve and
	  // the previous
	  shared_ptr<Vertex> common_vx = loop[ki-1]->getCommonVertex(loop[ki]);

	  // TESTING
	  if (common_vx->nmbUniqueEdges(this) == 2 ||
	      loop[ki-1]->geomCurve().get() == loop[ki]->geomCurve().get() ||
	      (loop[ki-1]->twin() &&
	       loop[ki-1]->twin()->geomEdge()->geomCurve().get() == loop[ki]->geomCurve().get()))
	    {
	      if (!(curr->twin() && prev->twin()))
		ang = M_PI;
	      else
		{
		  // Check for corner
		  if (prev->face() != curr->face())
		    loop[ki] = curr = curr->twin()->geomEdge();
		  if (prev->face() != curr->face())
		    ang = M_PI;
		  else
		    {
		      double t1 = prev->parAtVertex(common_vx.get());
		      double t2 = curr->parAtVertex(common_vx.get());
		      Point tan1 = prev->tangent(t1);
		      Point tan2 = curr->tangent(t2);
		      ang = tan1.angle(tan2);
		    }
		}
	    }
	}
      if (ang >= toptol_.bend)
	{
	  vector<ftEdge*> curr_loop;
	  curr_loop.push_back(curr);
	  joined_loop.push_back(curr_loop);
	}
      else
	{
	  joined_loop[joined_loop.size()-1].push_back(curr);
	}
      prev = curr;
    }
 
  if (joined_loop.size() < 4 && loop.size() >= 4)
    {
      // Ensure that there is enough curves to create a coons patch.
      // This split could be done in a more smart way, but don't 
      // expect this to be a probable case
      for (ki=0; ki<joined_loop.size(); ++ki)
	{
	  if (joined_loop[ki].size() > 1)
	    {
	      vector<ftEdge*> curr_loop(joined_loop[ki].begin()+1,
					joined_loop[ki].end());
	      joined_loop[ki].erase(joined_loop[ki].begin()+1, 
				    joined_loop[ki].end());
	      joined_loop.insert(joined_loop.begin()+ki, curr_loop);
	      if (joined_loop.size() == 4)
		break;
	    }
	}
    }
	      
  // Collect curve information
  double fuzzy = 1.0e-6;
  vector<vector<shared_ptr<ParamCurve> > > join_cvs(joined_loop.size());
  vector<vector<pair<double, double> > > par_bounds(joined_loop.size());
  for (ki=0; ki<joined_loop.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv1 = joined_loop[ki][0]->geomCurve();
      double t1 = joined_loop[ki][0]->tMin();
      double t2 = joined_loop[ki][0]->tMax();

      for (kj=0; kj<joined_loop[ki].size(); kj=kr)
	{
	  for (kr=kj+1; kr<joined_loop[ki].size(); ++kr)
	    {
	      shared_ptr<ParamCurve> tmp_cv2 = joined_loop[ki][kr]->geomCurve();
	      double t3 = joined_loop[ki][kr]->tMin();
	      double t4 = joined_loop[ki][kr]->tMax();
	      if (tmp_cv1.get() == tmp_cv2.get() &&
		  (fabs(t2-t3)<fuzzy || fabs(t1-t4)<fuzzy))
		{
		  Point joint = tmp_cv2->point(fabs(t2-t3) < fabs(t1-t4) ? t3 : t4);
		  joint_points.push_back(joint);
		  t1 = std::min(t1, t3);
		  t2 = std::max(t2, t4);
		}
	      else
		{
		  join_cvs[ki].push_back(tmp_cv1);
		  par_bounds[ki].push_back(make_pair(t1, t2));
		  tmp_cv1 = tmp_cv2;
		  t1 = t3;
		  t2 = t4;
		  break;
		}
	    }
	}
      join_cvs[ki].push_back(tmp_cv1);
      par_bounds[ki].push_back(make_pair(t1, t2));
    }
  
  
  // Fetch curves, and make parameter curves corresponding to the volume
  // Join curves that should belong to the same boundary curve of the
  // missing surface, but are represented as several edges
  // Remember the positions at these joints
  space_cvs.resize(join_cvs.size());
  for (ki=0; ki<join_cvs.size(); ++ki)
    {
      shared_ptr<ParamCurve> tmp_cv = 
	shared_ptr<ParamCurve>(join_cvs[ki][0]->subCurve(par_bounds[ki][0].first,
							 par_bounds[ki][0].second, fuzzy));
      for (kr=1; kr<join_cvs[ki].size(); ++kr)
	{
	  shared_ptr<ParamCurve> tmp_cv2 = 
	    shared_ptr<ParamCurve>(join_cvs[ki][kr]->subCurve(par_bounds[ki][kr].first,
							      par_bounds[ki][kr].second, fuzzy));

	  // Make sure that the curves are consistently oriented
	  Point pos1 = tmp_cv->point(tmp_cv->startparam());
	  Point pos2 = tmp_cv->point(tmp_cv->endparam());
	  Point pos3 = tmp_cv2->point(tmp_cv2->startparam());
	  Point pos4 = tmp_cv2->point(tmp_cv2->endparam());
	  double d1 = pos1.dist(pos3);
	  double d2 = pos1.dist(pos4);
	  double d3 = pos2.dist(pos3);
	  double d4 = pos2.dist(pos4);
	  if (std::min(d3, d4) > std::min(d1, d2))
	    {
	      tmp_cv->reverseParameterDirection();
	      std::swap(pos1, pos2);
	      std::swap(d1, d3);
	      std::swap(d2, d4);
	    }
	  if (d4 < d3)
	    tmp_cv2->reverseParameterDirection();
	  Point joint = tmp_cv2->point(tmp_cv2->startparam());
	  joint_points.push_back(joint);
	  double dist;
	  vector<Point> pts1(2);
	  tmp_cv->point(pts1, tmp_cv->endparam(), 1);
	  vector<Point> pts2(2);
	  tmp_cv2->point(pts2, tmp_cv2->startparam(), 1);
	  double fac = pts2[1].length()/pts1[1].length();

	  // TEST
	  fac = 1.0;

	  double s1 = tmp_cv->startparam();
	  double s2 = tmp_cv->endparam();
	  double t1 = tmp_cv2->startparam();
	  double t2 = tmp_cv2->endparam();
	  double len1 = tmp_cv->estimatedCurveLength();
	  double len2 = tmp_cv2->estimatedCurveLength();
	  fac = len2*(s2-s1)/(len1*(t2-t1));

	  tmp_cv2->setParameterInterval(t1, t1+fac*(t2-t1));
	  tmp_cv->appendCurve(tmp_cv2.get(), 0, dist, false);
	}
      space_cvs[ki] = tmp_cv;
    }
  }

//===========================================================================
ftEdge*  ftVolume::getLeftLoopEdge(ftSurface* face, Body *bd,
				   shared_ptr<EdgeVertex> radial)
//===========================================================================
{
  double tol = 1.0e-4; //1.0e-6;
  double tol2 = 1.0e-12; //1.0e-8; Temporary fix. It is a danger of taking the wrong
  // decision if the numbers get too small. Need a better solution

  ftEdge *edge = 0;
  vector<ftEdge*> edgs = radial->uniqueEdges(bd);
  if (edgs.size() == 0)
    return edge;

  // Fetch a point on the edge, and compute face characteristics at
  // the edge
  Point dummy_vec;
  ftEdge *curr = edgs[0];
  Point mid = curr->point(0.5*(curr->tMin()+curr->tMax()));
  double u0, v0, d0, t0;
  Point close0;
  ftEdgeBase* e0 = face->closestBoundaryPoint(mid, dummy_vec, u0, v0, close0,
					       d0, t0);
  Point norm = e0->normal(t0);
  Point tan = e0->tangent(t0);
  Point vec = norm.cross(tan);

#ifdef DEBUG_VOL1
  std::ofstream of("faces.g2");
#endif

  // Check for a configuaration where the two pairs of faces are
  // each other twins
  bool twin_config = false;
  size_t ki, kj;
  for (ki=0; ki<edgs.size(); ++ki)
    {
      for (kj=ki+1; kj<edgs.size(); ++kj)
	{
	  if (!(edgs[ki]->twin() && edgs[kj]->twin()))
	    continue;
	  ftFaceBase *f1 = edgs[ki]->face();
	  ftFaceBase *f2 = edgs[ki]->twin()->geomEdge()->face();
	  ftFaceBase *f3 = edgs[kj]->face();
	  ftFaceBase *f4 = edgs[kj]->twin()->geomEdge()->face();
	  if (f1 && f2 && f3 && f4)
	    {
	      ftSurface *f5 = f1->asFtSurface();
	      ftSurface *f6 = f2->asFtSurface();
	      ftSurface *f7 = f3->asFtSurface();
	      ftSurface *f8 = f4->asFtSurface();
	      if (f5 && f6 && f7 && f8)
		{
		  if ((f5->twin() == f7 && f6->twin() == f8) ||
		      (f5->twin() == f8 && f6->twin() == f7))
		    twin_config = true;
		}
	    }
	}
    }

  // For each unique edge, check if the corresponding edge pair is closest
  // to the given face
  int idx = -1;
  double min_ang = MAXDOUBLE;
  double fac1 = 0.0, fac2 = 0.0;
  for (ki=0; ki<edgs.size(); ++ki)
    {
      if (!(edgs[ki]->face()  && edgs[ki]->twin()))
	continue;
      
      ftFaceBase* f1 = edgs[ki]->face();
      ftFaceBase* f2 = edgs[ki]->twin()->face();
      double t1, t2;
      ftEdgeBase* e1 = f1->asFtSurface()->closestBoundaryPoint(mid, dummy_vec, u0, v0, 
							       close0, d0, t1);
      ftEdgeBase* e2 = f2->asFtSurface()->closestBoundaryPoint(mid, dummy_vec, u0, v0, 
							       close0, d0, t2);
#ifdef DEBUG_VOL1
      shared_ptr<ParamSurface> sf = f1->asFtSurface()->surface();
      sf->writeStandardHeader(of);
      sf->write(of);
      sf = f2->asFtSurface()->surface();
      sf->writeStandardHeader(of);
      sf->write(of);
      Point pos = e1->point(t1);
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << "1" << std::endl;
      of << pos << std::endl;
#endif

      Point norm1 = e1->normal(t1);
      Point tan1 = e1->tangent(t1);
      Point vec1 = norm1.cross(tan1);
      Point norm2 = e2->normal(t2);
      Point tan2 = e2->tangent(t2);
      Point vec2 = norm2.cross(tan2);

      double ang1 = vec.angle(vec1);
      double ang2 = vec.angle(vec2);
      if (vec*norm1 > tol2)
	ang1 = 2*M_PI - ang1;
      //ang1 += M_PI;
      if (vec*norm2 > tol2)
	ang2 = 2*M_PI - ang2;
	//ang2 += M_PI;
      if (twin_config && fabs(ang1+ang2-min_ang) < tol)
	{
	  // Undecided situation. Check normal vectors
	  if (ang1 <= ang2)
	    {
	    if ((vec*norm1 < 0.0 && fac1 > 0.0) ||
		((vec*norm1)*fac1 > 0.0 && vec*norm2 < 0.0 && fac2 > 0.0))
	      {
		idx = (int)ki;
		min_ang = ang1 + ang2;
	      }
	    }
	  else
	    {
	    if ((vec*norm2 < 0.0 && fac2 > 0.0) ||
		((vec*norm2)*fac2 > 0.0 && vec*norm1 < 0.0 && fac1 > 0.0))
	      {
		idx = (int)ki;
		min_ang = ang1 + ang2;
	      }
	    }
	}
      else if (ang1 + ang2 < min_ang)
	{
	  idx = (int)ki;
	  min_ang = ang1 + ang2;
	}

      // Save info from current
      fac1 = vec*norm1;
      fac2 = vec*norm2;
    }

  if (idx >= 0)
    return edgs[idx];
  else
    return edge;
}

//===========================================================================
bool  ftVolume::doSwapEdges(ftSurface* face, ftEdge* edge1, ftEdge *edge2)
//===========================================================================
{
  // The connection have been lost in previous calls to
  // this function. Check if the loop is consistent
  // Compute the vector pointing into the new faces from
  // the midpoint of the edge
  Point mid = edge1->point(0.5*(edge1->tMin()+edge1->tMax()));
  double u0, v0, d0, t0;
  Point close0;
  Point dummy_vec;
  ftEdgeBase* e0 = face->closestBoundaryPoint(mid, dummy_vec, u0, v0, close0,
					       d0, t0);
  Point norm = e0->normal(t0);
  Point tan = e0->tangent(t0);
  Point vec = norm.cross(tan);

  // For each candidate edge, compute the vectors pointing into
  // the associated faces
  ftFaceBase* f1 = edge1->face();
  ftFaceBase* f2 = edge1->twin()->face();
  ftFaceBase* f3 = edge2->face();
  ftFaceBase* f4 = edge2->twin()->face();
  double t1, t2, t3, t4;
  ftEdgeBase* e1 = f1->asFtSurface()->closestBoundaryPoint(mid, dummy_vec, u0, v0, close0,
							   d0, t1);
  ftEdgeBase* e2 = f2->asFtSurface()->closestBoundaryPoint(mid, dummy_vec, u0, v0, close0,
							   d0, t2);
  ftEdgeBase* e3 = f3->asFtSurface()->closestBoundaryPoint(mid, dummy_vec, u0, v0, close0,
							   d0, t3);
  ftEdgeBase* e4 = f4->asFtSurface()->closestBoundaryPoint(mid, dummy_vec, u0, v0, close0,
							   d0, t4);
#ifdef DEBUG_VOL1
  std::ofstream of("faces.g2");
  shared_ptr<ParamSurface> sf = f1->asFtSurface()->surface();
  sf->writeStandardHeader(of);
  sf->write(of);
  sf = f2->asFtSurface()->surface();
  sf->writeStandardHeader(of);
  sf->write(of);
  sf = f3->asFtSurface()->surface();
  sf->writeStandardHeader(of);
  sf->write(of);
  sf = f4->asFtSurface()->surface();
  sf->writeStandardHeader(of);
  sf->write(of);
#endif

  Point norm1 = e1->normal(t1);
  Point tan1 = e1->tangent(t1);
  Point vec1 = norm1.cross(tan1);
  Point norm2 = e2->normal(t2);
  Point tan2 = e2->tangent(t2);
  Point vec2 = norm2.cross(tan2);
  Point norm3 = e3->normal(t3);
  Point tan3 = e3->tangent(t3);
  Point vec3 = norm3.cross(tan3);
  Point norm4 = e4->normal(t4);
  Point tan4 = e4->tangent(t4);
  Point vec4 = norm4.cross(tan4);

  double ang1 = vec.angle(vec1);
  double ang2 = vec.angle(vec2);
  double ang3 = vec.angle(vec3);
  double ang4 = vec.angle(vec4);
  ftSurface *f1_2 = f1->asFtSurface();
  ftSurface *f2_2 = f2->asFtSurface();
  if (f1_2 && f2_2 && ((f1_2->twin() && f1_2->twin() == f3) ||
		       (f2_2->twin() && f2_2->twin() == f4)))
    {
      // Special case. Choose alternative to get normals pointing out
      // of the final volume
      if (vec*norm1 > 0.0 && vec*norm3 < 0.0)
	return true;
      else if (vec*norm1 < 0.0 && vec*norm3 > 0.0)
	return false;
      else if (vec*norm2 > 0.0 && vec*norm4 < 0.0)
	return true;
      else if (vec*norm2 < 0.0 && vec* norm4 > 0.0)
	return false;
      else if ((ang1 > ang3 && (vec1*norm)*(vec3*norm) > 0.0) ||
	       (ang1 > ang4 && (vec1*norm)*(vec4*norm) > 0.0) ||
	       (ang2 > ang3 && (vec2*norm)*(vec3*norm) > 0.0) ||
	       (ang2 > ang4 && (vec2*norm)*(vec4*norm) > 0.0))
	return true;
      else 
	return false;
    }
  else if ((ang1 > ang3 && (vec1*norm)*(vec3*norm) > 0.0) ||
	   (ang1 > ang4 && (vec1*norm)*(vec4*norm) > 0.0) ||
	   (ang2 > ang3 && (vec2*norm)*(vec3*norm) > 0.0) ||
	   (ang2 > ang4 && (vec2*norm)*(vec4*norm) > 0.0))
    return true;
  else 
    return false;
}

//===========================================================================
vector<vector<ftEdge*> > 
ftVolume::getMissingSfLoops(vector<pair<Point,Point> >& corr_vx_pts,
			    bool perform_step2, bool smooth_connections,
			    int max_nmb)
//===========================================================================
{
  vector<vector<ftEdge*> > missing_sf_loops;

#ifdef DEBUG_VOL1
  std::ofstream of10("pre_complete.g2");
  int nmb_sfs = shells_[0]->nmbEntities();
  for (int ka=0; ka<nmb_sfs; ++ka)
    {
      shared_ptr<ParamSurface> tmp_sf = shells_[0]->getSurface(ka);
      tmp_sf->writeStandardHeader(of10);
      tmp_sf->write(of10);
    }
#endif

  // Generate missing edges in the curve net
  vector<shared_ptr<ftEdge> > start_edges;
  CompleteEdgeNet reg(shells_[0], perform_step2, smooth_connections);
  bool done = reg.perform(corr_vx_pts);
#ifdef DEBUG_VOL1
  std::cout << "Perform: " << done << std::endl;

  std::ofstream of0("missing_edges.g2");
#endif
  vector<pair<shared_ptr<Vertex>,shared_ptr<Vertex> > > vx_pair =
    reg.getMissingEdges();
  int ki, kr;
  size_t kj;
  //size_t ki, kj, kr, kh;
  for (ki=0; ki<(int)vx_pair.size(); ++ki)
    {
      // Create edge
      shared_ptr<Vertex> vx1 = vx_pair[ki].first;
      shared_ptr<Vertex> vx2 = vx_pair[ki].second;
      shared_ptr<ParamCurve> cv = 
	//makeMissingEdgeCv(vx1, vx2);
	shared_ptr<ParamCurve>(new SplineCurve(vx1->getVertexPoint(),
					       vx2->getVertexPoint()));

#ifdef DEBUG_VOL1
      of0 << "100 1 0 4 0 255 0 255" << std::endl;
      cv->write(of0);
      of0 << std::endl;
#endif

      shared_ptr<ftEdge> edge = 
	shared_ptr<ftEdge>(new ftEdge(cv, cv->startparam(), cv->endparam()));

      // Insert the edge into the data structure
      // The vertices only keep pointers so it is important to store the
      // shared pointer to edge
      edge->joinVertex(edge->getVertex(true), vx1);
      edge->joinVertex(edge->getVertex(false), vx2);
      missing_edges_.push_back(edge);

      start_edges.push_back(edge);
    }
  int nmb_missing_edges = (int)start_edges.size();

  // Fetch other candidate start edges
  // TEST
  vector<shared_ptr<ftEdge> > tmp_edges = getStartEdges();
  //vector<shared_ptr<ftEdge> > tmp_edges = shells_[0]->getUniqueInnerEdges();
  // END TEST
  start_edges.insert(start_edges.end(), tmp_edges.begin(), tmp_edges.end());

#ifdef DEBUG_VOL1
  std::ofstream of("start_edges.g2");
  for (size_t kh=0; kh<start_edges.size(); ++kh)
    {
      shared_ptr<ParamCurve> cv = start_edges[kh]->geomCurve();
      shared_ptr<ParamCurve> cv2 = 
	shared_ptr<ParamCurve>(cv->subCurve(start_edges[kh]->tMin(),
					    start_edges[kh]->tMax()));
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
      if (sfcv.get())
	{
	  if (sfcv->spaceCurve().get())
	    {
	      sfcv->spaceCurve()->writeStandardHeader(of);
	      sfcv->spaceCurve()->write(of);
	    }
	}
      else
	{
	  cv2->writeStandardHeader(of);
	  cv2->write(of);
	}
    }
#endif

  // Traverse edge loops and fetch missing surface loops
  for (ki=0; ki<(int)start_edges.size(); )
    {
#ifdef DEBUG_VOL1
  std::ofstream of10("remaining_startedges.g2");
  for (size_t kh=0; kh<start_edges.size(); ++kh)
    {
      shared_ptr<ParamCurve> cv = start_edges[kh]->geomCurve();
      shared_ptr<ParamCurve> cv2 = 
	shared_ptr<ParamCurve>(cv->subCurve(start_edges[kh]->tMin(),
					    start_edges[kh]->tMax()));
      shared_ptr<CurveOnSurface> sfcv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
      if (sfcv.get())
	{
	  if (sfcv->spaceCurve().get())
	    {
	      sfcv->spaceCurve()->writeStandardHeader(of10);
	      sfcv->spaceCurve()->write(of10);
	    }
	}
      else
	{
	  cv2->writeStandardHeader(of10);
	  cv2->write(of10);
	}
    }
#endif
  vector<vector<ftEdge*> >  loops = getLoop(start_edges[ki], max_nmb);

#ifdef DEBUG_VOL1
	  std::ofstream of4("curr_loops0.g2");
	  for (size_t kj1=0; kj1<loops.size(); ++kj1)
	    for (size_t kj2=0; kj2<loops[kj1].size(); ++kj2)
	      {
		ftEdge *e1 = loops[kj1][kj2];
		shared_ptr<ParamCurve> cv =
		  shared_ptr<ParamCurve>(e1->geomCurve()->subCurve(e1->tMin(), e1->tMax()));
		shared_ptr<CurveOnSurface> sfcv = 
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
		if (sfcv.get())
		  {
		    if (sfcv->spaceCurve().get())
		      {
			sfcv->spaceCurve()->writeStandardHeader(of4);
			sfcv->spaceCurve()->write(of4);
		      }
		  }
		else
		  {
		    cv->writeStandardHeader(of4);
		    cv->write(of4);
		  }
	      }
#endif

      // Check if any loop is found already
      for (kj=0; kj<loops.size(); )
	{
	  if (loopExisting(loops[kj], missing_sf_loops))
	    loops.erase(loops.begin()+kj);
	  else
	    kj++;
	}

      if (loops.size() == 0)
	{
	  // missing_sf_loops.clear();
	  // return missing_sf_loops;  // Not possible to regularize
	  ++ki;
	  continue;
	}

#ifdef DEBUG_VOL1
	  std::ofstream of3("curr_loops.g2");
	  for (size_t kj1=0; kj1<loops.size(); ++kj1)
	    for (size_t kj2=0; kj2<loops[kj1].size(); ++kj2)
	      {
		ftEdge *e1 = loops[kj1][kj2];
		shared_ptr<ParamCurve> cv =
		  shared_ptr<ParamCurve>(e1->geomCurve()->subCurve(e1->tMin(), e1->tMax()));
		shared_ptr<CurveOnSurface> sfcv = 
		  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
		if (sfcv.get())
		  {
		    if (sfcv->spaceCurve().get())
		      {
			sfcv->spaceCurve()->writeStandardHeader(of3);
			sfcv->spaceCurve()->write(of3);
		      }
		  }
		else
		  {
		    cv->writeStandardHeader(of3);
		    cv->write(of3);
		  }
	      }
#endif

      // Remove affected start edges. Always remove the current edge
      // as it is already traversed
      start_edges.erase(start_edges.begin() + ki);
      nmb_missing_edges--;
      for (kj=0; kj<loops.size(); ++kj)
	{
	  for (size_t kh=1; kh<loops[kj].size(); ++kh)
	    {
	      vector<ftEdge*> loop_edges;
	      if (loops[kj][kh]->hasEdgeMultiplicity())
		loop_edges = loops[kj][kh]->getEdgeMultiplicityInstance()->allEdges(this);
	      else
		{
		  loop_edges.push_back(loops[kj][kh]);
		  if (loops[kj][kh]->twin())
		    loop_edges.push_back(loops[kj][kh]->twin()->geomEdge());
		}

	      for (kr=0; kr<(int)start_edges.size(); ++kr)
		{
		  vector<ftEdge*> start_edg;
		  if (start_edges[kr]->hasEdgeMultiplicity())
		    start_edg = start_edges[kr]->getEdgeMultiplicityInstance()->allEdges(this);
		  else
		    {
		      start_edg.push_back(start_edges[kr].get());
		      if (start_edges[kr]->twin())
			start_edg.push_back(start_edges[kr]->twin()->geomEdge());
		    }

		  size_t ix1, ix2;
		  for (ix1=0; ix1<loop_edges.size(); ix1++)
		    {
		      for (ix2=0; ix2<start_edg.size(); ++ix2)
			if (loop_edges[ix1] == start_edg[ix2])
			  break;
		      if (ix2 < start_edg.size())
			break;
		    }
		  if (ix1 < loop_edges.size())
		    break;
		}
	      if (kr < (int)start_edges.size() && kr >= nmb_missing_edges)
		{
		  start_edges.erase(start_edges.begin()+kr);
		  if (kr < ki)
		    ki--;
		}
	    }
	}

      missing_sf_loops.insert(missing_sf_loops.end(), loops.begin(), loops.end());
    }
#ifdef DEBUG_VOL1
  std::ofstream of2("missing_loops.g2");
  for (ki=0; ki<(int)missing_sf_loops.size(); ++ki)
    for (kj=0; kj<missing_sf_loops[ki].size(); ++kj)
      {
	ftEdge *e1 = missing_sf_loops[ki][kj];
	shared_ptr<ParamCurve> cv =
	  shared_ptr<ParamCurve>(e1->geomCurve()->subCurve(e1->tMin(), e1->tMax()));
	shared_ptr<CurveOnSurface> sfcv = 
	  dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
	if (sfcv.get())
	  {
	    if (sfcv->spaceCurve().get())
	      {
		sfcv->spaceCurve()->writeStandardHeader(of2);
		sfcv->spaceCurve()->write(of2);
	      }
	  }
	else
	  {
	    cv->writeStandardHeader(of2);
	    cv->write(of2);
	  }
      }
#endif

  return missing_sf_loops;
}

//===========================================================================
// 
// 
bool ftVolume::loopExisting(vector<ftEdge*>& loop, 
			    vector<vector<ftEdge*> >& curr_loops)
//===========================================================================
{
  size_t ki, kj, kr, ix1, ix2;
  for (ki=0; ki<curr_loops.size(); ++ki)
    {
      for (kj=0; kj<loop.size(); ++kj)
	{
	  // Collect all half edges related to edge
	  vector<ftEdge*> edg1;
	  if (loop[kj]->hasEdgeMultiplicity())
	    edg1 = loop[kj]->getEdgeMultiplicityInstance()->allEdges(this);
	  else
	    {
	      edg1.push_back(loop[kj]);
	      if (loop[kj]->twin())
		edg1.push_back(loop[kj]->twin()->geomEdge());
	    }
	  for (kr=0; kr<curr_loops[ki].size(); ++kr)
	    {
	      // Collect all half edges related to edge in alternative loop
	      vector<ftEdge*> edg2;
	      if (curr_loops[ki][kr]->hasEdgeMultiplicity())
		edg2 = curr_loops[ki][kr]->getEdgeMultiplicityInstance()->allEdges(this);
	      else
		{
		  edg2.push_back(curr_loops[ki][kr]);
		  if (curr_loops[ki][kr]->twin())
		    edg2.push_back(curr_loops[ki][kr]->twin()->geomEdge());
		}

	      // Check identity between current edge pair
	      for (ix1=0; ix1<edg1.size(); ++ix1)
		{
		  for (ix2=0; ix2<edg2.size(); ++ix2)
		    if (edg1[ix1] == edg2[ix2])
		      break;   // Same edge

		  if (ix2 < edg2.size())
		    break;  // Same edge
		}

	      if (ix1 < edg1.size())
		break;  // Edge identity is found
	    }

	  if (kr == curr_loops[ki].size())
	    break;  // No identity with this loop edge
	}
      if (kj == loop.size())
	return true;  // Identity with this loop
    }
  return false;  // No identity is found
}

//===========================================================================
// 
// 
vector<shared_ptr<ftEdge> > ftVolume::getStartEdges()
//===========================================================================
{
  // Fetch all unique edges. First fetch shells
  vector<shared_ptr<SurfaceModel> > shells = getAllShells();
  vector<shared_ptr<ftEdge> > edges;
  size_t ki;
  for (ki=0; ki<shells.size(); ++ki)
    {
      vector<shared_ptr<ftEdge> > curr_edges = 
	shells[ki]->getUniqueInnerEdges();
      edges.insert(edges.end(), curr_edges.begin(), curr_edges.end());
    }

  // Keep the edges that split one underlying surface
  vector<shared_ptr<ftEdge> > start_edges;
  for (ki=0; ki<edges.size(); ++ki)
    {
      if (!edges[ki]->twin())
	continue; // Unexpected for a solid, but definitely not a splitting edge

      if (!edges[ki]->face())
	continue; // Newly constructed edge, treated elsewhere

      shared_ptr<ParamSurface> sf1 = 
	edges[ki]->face()->asFtSurface()->surface();
      shared_ptr<ParamSurface> sf2 = 
	edges[ki]->twin()->geomEdge()->face()->asFtSurface()->surface();

      // Check for trimmed surfaces
      shared_ptr<BoundedSurface> bd_sf1 = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf1);
      shared_ptr<BoundedSurface> bd_sf2 = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf2);
      // if (!(bd_sf1.get() && bd_sf2.get()))
      // 	continue;
      if (bd_sf1.get() && bd_sf2.get() &&
	  bd_sf1->underlyingSurface().get() ==  bd_sf2->underlyingSurface().get())
	start_edges.push_back(edges[ki]);  // A start edge is found
      else
	{
	  // An edge may be relevant for a loop around a missing surface
	  // even if the underlying surface on both sides are different.
	  // This can be due to merging of surfaces across a seam
	  // Count the number of bodies meeting in this edge, and the number of
	  // faces for each body
	  // First check continuity in the current body
	  shared_ptr<Vertex> vx1 = edges[ki]->getVertex(true);
	  shared_ptr<Vertex> vx2 = edges[ki]->getVertex(false);
	  Point norm1 = edges[ki]->normal(edges[ki]->tMin());
	  double t1 = edges[ki]->twin()->geomEdge()->parAtVertex(vx1.get());
	  Point norm2 = edges[ki]->twin()->normal(t1);
	  Point norm3 = edges[ki]->normal(edges[ki]->tMax());
	  double t2 = edges[ki]->twin()->geomEdge()->parAtVertex(vx2.get());
	  Point norm4 = edges[ki]->twin()->normal(t2);
	  double angtol = toptol_.bend;
	  if (norm1.angle(norm2) < angtol && norm3.angle(norm4) < angtol)
	    {
	      // Sufficient continuity. Check configuration
	      vector<ftSurface*> faces = edges[ki]->getAllAdjacentFaces();

	      // Check equality of faces
	      size_t kr, kh;
	      for (kr=0; kr<faces.size(); ++kr)
		  for (kh=kr+1; kh<faces.size(); )
		    {
		      if (faces[kr] == faces[kh])
			faces.erase(faces.begin()+kh);
		      else
			kh++;
		    }

	      vector<Body*> bodies;
	      vector<int> nmb_faces;
	      for (size_t kj=0; kj<faces.size(); ++kj)
		{
		  Body *bd = faces[kj]->getBody();
		  for (kr=0; kr<bodies.size(); ++kr)
		    {
		      if (bd == bodies[kr])
			break;
		    }
		  if (kr == bodies.size())
		    {
		      bodies.push_back(bd);
		      nmb_faces.push_back(1);
		    }
		  else 
		    nmb_faces[kr]++;
		}
	  
	      // Check if the number of associated faces for eah body is always two
	      for (kr=0; kr<nmb_faces.size(); ++kr)
		if (nmb_faces[kr] != 2)
		  break;

	      if ((bodies.size() == 1 || bodies.size() == 2) &&
		  kr == nmb_faces.size())
		start_edges.push_back(edges[ki]);  // A start edge is found
	    }
	}
    }

  return start_edges;
}

//===========================================================================
// 
// 
vector<vector<ftEdge*> > ftVolume::getLoop(shared_ptr<ftEdge> start_edge,
					   int max_nmb)
//===========================================================================
{
  vector<vector<ftEdge*> > all_loops;

  // Fetch vertex to start traversal
  shared_ptr<Vertex> start_vx = start_edge->getVertex(true);
  shared_ptr<Vertex> vx = start_edge->getVertex(false);

  // Follow all paths from the start edge
  vector<ftEdge*> edges = vx->uniqueEdges(this);
  vector<ftEdge*> free_edges = vx->freeEdges();
  if (free_edges.size() > 0)
    edges.insert(edges.end(), free_edges.begin(), free_edges.end());

  // Remove edges belonging to the same edge multiplicity
  size_t ki, kj, kr;
  for (ki=0; ki<edges.size(); ++ki)
    for (kj=ki+1; kj<edges.size();)
      {
	if (edges[ki]->hasCommonRadialEdge(edges[kj]))
	  edges.erase(edges.begin()+kj);
	else
	  kj++;
      }
  
  vector<ftEdge*> start_edges;
  if (start_edge->hasEdgeMultiplicity())
    start_edges = start_edge->getEdgeMultiplicityInstance()->allEdges(this);
  else
    {
      start_edges.push_back(start_edge.get());
      if (start_edge->twin())
	start_edges.push_back(start_edge->twin()->geomEdge());
    }
  for (ki=0; ki<edges.size(); ++ki)
    {
      vector<ftEdge*> edge_loop;

      vector<ftEdge*> curr_edges;
      if (edges[ki]->hasEdgeMultiplicity())
	curr_edges = edges[ki]->getEdgeMultiplicityInstance()->allEdges(this);
      else
	{
	  curr_edges.push_back(edges[ki]);
	  if (edges[ki]->twin())
	    curr_edges.push_back(edges[ki]->twin()->geomEdge());
	}

      for (kj=0; kj<start_edges.size(); ++kj)
	{
	  for (kr=0; kr<curr_edges.size(); ++kr)
	    if (start_edges[kj] == curr_edges[kr])
	      break;
	  if (kr < curr_edges.size())
	    break;
	}
      if (kj < start_edges.size())
	continue;
      
      shared_ptr<Vertex> vx2 = edges[ki]->getOtherVertex(vx.get());
      if (vx2.get() == start_vx.get())
	continue;  // Degenerate loop

      if ((start_edge->face() != edges[ki]->face()) &&
	  (edges[ki]->twin() && start_edge->face() == 
	   edges[ki]->twin()->geomEdge()->face()))
	edges[ki] = edges[ki]->twin()->geomEdge();
      edge_loop.push_back(start_edge.get());
      edge_loop.push_back(edges[ki]);
      bool found = getLoopEdges(edge_loop, start_vx, vx, max_nmb);
      if (found)
	all_loops.push_back(edge_loop);
    }
  return all_loops;
}

//===========================================================================
// 
// 
bool ftVolume::getLoopEdges(vector<ftEdge*>& loop,
			    shared_ptr<Vertex> start_vx,
			    shared_ptr<Vertex> vx,
			    int max_nmb)
//===========================================================================
{
  int added_max = 0;

  if (loop.size() == 0)
    return false;  // No start edge

  if ((int)loop.size() > max_nmb)
    return false;  // Too many edges in loop

  // Get the other vertex corresponding to the start edge
  ftEdge* prev = loop[loop.size()-1];
  shared_ptr<Vertex> vx2 = prev->getOtherVertex(vx.get());

  // Follow all paths from this vertex
  vector<ftEdge*> edges = vx2->uniqueEdges(this);
  vector<ftEdge*> free_edges = vx2->freeEdges();
  if (free_edges.size() > 0)
    edges.insert(edges.end(), free_edges.begin(), free_edges.end());

  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      vector<ftEdge*> curr_edges;
      if (edges[ki]->hasEdgeMultiplicity())
	curr_edges = edges[ki]->getEdgeMultiplicityInstance()->allEdges(this);
      else
	{
	  curr_edges.push_back(edges[ki]);
	  if (edges[ki]->twin())
	    curr_edges.push_back(edges[ki]->twin()->geomEdge());
	}

      size_t kj, kr, kh;
      for (kh=0; kh<loop.size(); ++kh)
	{
	  vector<ftEdge*> prev_edges;
	  if (loop[kh]->hasEdgeMultiplicity())
	    prev_edges = loop[kh]->getEdgeMultiplicityInstance()->allEdges(this);
	  else
	    {
	      prev_edges.push_back(loop[kh]);
	      if (loop[kh]->twin())
		prev_edges.push_back(loop[kh]->twin()->geomEdge());
	    }

	  for (kj=0; kj<prev_edges.size(); ++kj)
	    {
	      for (kr=0; kr<curr_edges.size(); ++kr)
		if (prev_edges[kj] == curr_edges[kr])
		  break;
	      if (kr < curr_edges.size())
		break;
	    }
	  if (kj < prev_edges.size())
	    break;  
	}
      if (kh < loop.size())
	continue; // The loop turns back on itself

      shared_ptr<Vertex> end_vx = edges[ki]->getOtherVertex(vx2.get());
      if (end_vx.get() == start_vx.get() && loop.size() < 2 /*3*/)
	continue;  // Repeated vertex in loop

      if (loop.size() < 3)
	{
	  shared_ptr<Vertex> tmp_vx = start_vx;
	  size_t kj;
	  for (kj=0; kj<loop.size()-1; ++kj)
	    {
	      shared_ptr<Vertex> tmp_vx2 = loop[kj]->getOtherVertex(tmp_vx.get());
	      if (tmp_vx2.get() == end_vx.get())
		break;
	      tmp_vx = tmp_vx2;
	    }
	  if (kj < loop.size()-1)
	    continue;  // Repeated vertex in loop
	}

      
      ftEdge *prev2 = loop[loop.size()-1];
      if ((prev2->face() != edges[ki]->face()) &&
	  (edges[ki]->twin() && prev2->face() == edges[ki]->twin()->geomEdge()->face()))
	edges[ki] = edges[ki]->twin()->geomEdge();
      loop.push_back(edges[ki]);

      // Check if the maximum number of edges in the loop should be increased
      shared_ptr<Vertex> common_vx = 
	loop[loop.size()-2]->getCommonVertex(loop[loop.size()-1]);
      if (common_vx.get() && common_vx->nmbUniqueEdges(this) == 2)
	{
	  // Check for corner
	  double t1 = loop[loop.size()-2]->parAtVertex(common_vx.get());
	  double t2 = loop[loop.size()-1]->parAtVertex(common_vx.get());
	  Point tan1 = loop[loop.size()-2]->tangent(t1);
	  Point tan2 = loop[loop.size()-1]->tangent(t2);
	  double ang = tan1.angle(tan2);
	  if (ang < toptol_.bend /*|| fabs(M_PI-ang) < toptol_.bend*/)
	    {
	      // No T-joint
	      max_nmb++;
	      added_max = 1;
	    }
	}

 #ifdef DEBUG_VOL1
     std::ofstream of("edge_loop.g2");
      for (size_t kr=0; kr<loop.size(); ++kr)
	{
	  shared_ptr<ParamCurve> crv = loop[kr]->geomCurve();
	  shared_ptr<CurveOnSurface> sf_crv = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv);
	  if (sf_crv.get())
	    {
	    if (sf_crv->spaceCurve().get())
	      {
		sf_crv->spaceCurve()->writeStandardHeader(of);
		sf_crv->spaceCurve()->write(of);
	      }
	    }
	  else
	    {
	      crv->writeStandardHeader(of);
	      crv->write(of);
	    }
	}
#endif

  if ((int)loop.size() > max_nmb)
    {
      loop.pop_back();
      return false;  // Too many edges in loop, a new test after the last
      // added edge
    }

      if (end_vx.get() == start_vx.get())
	{
	  // A loop is found. Check if a surface already exists
	  bool sf_exist = sameFace(loop);

	  if (!sf_exist)
	    {
	      // A new loop is found

	      // Check if the plane(s) defined by the edge loop are significantly
	      // different from the tangent planes of the assiciated faces
	      // To check if the loop should be prosessed further
	      bool plane_found = checkPlaneLoop(loop);
	      if (!plane_found)
		{
		  // Check if an alternative and better loop configuration exists
		  for (size_t ki2=ki+1; ki2<edges.size(); ++ki2)
		    {
		      vector<ftEdge*> curr_edg2;
		      if (edges[ki2]->hasEdgeMultiplicity())
			curr_edg2 = edges[ki2]->getEdgeMultiplicityInstance()->allEdges(this);
		      else
			{
			  curr_edg2.push_back(edges[ki2]);
			  if (edges[ki2]->twin())
			    curr_edg2.push_back(edges[ki2]->twin()->geomEdge());
			}

		      for (kh=0; kh<loop.size()-1; ++kh)
			{
			  vector<ftEdge*> prev_edg2;
			  if (loop[kh]->hasEdgeMultiplicity())
			    prev_edg2 = loop[kh]->getEdgeMultiplicityInstance()->allEdges(this);
			  else
			    {
			      prev_edg2.push_back(loop[kh]);
			      if (loop[kh]->twin())
				prev_edg2.push_back(loop[kh]->twin()->geomEdge());
			    }

			  for (kj=0; kj<prev_edg2.size(); ++kj)
			    {
			      for (kr=0; kr<curr_edg2.size(); ++kr)
				if (prev_edg2[kj] == curr_edg2[kr])
				  break;
			      if (kr < curr_edg2.size())
				break;
			    }
			  if (kj < prev_edg2.size())
			    break;  
			}
		      if (kh < loop.size()-1)
			continue; // The loop turns back on itself

		      shared_ptr<Vertex> end_vx2 = edges[ki2]->getOtherVertex(vx2.get());

		      if (end_vx2.get() == start_vx.get())
			{
			  // An alternative loop is found
			  vector<ftEdge*> loop2 = loop;
			  loop2.pop_back();
			  loop2.push_back(edges[ki2]);
#ifdef DEBUG_VOL1
			  std::ofstream of("edge_loop2.g2");
			  for (size_t kr=0; kr<loop2.size(); ++kr)
			    {
			      shared_ptr<ParamCurve> crv = loop2[kr]->geomCurve();
			      shared_ptr<CurveOnSurface> sf_crv = 
				dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv);
			      if (sf_crv.get())
				{
				  if (sf_crv->spaceCurve().get())
				    {
				      sf_crv->spaceCurve()->writeStandardHeader(of);
				      sf_crv->spaceCurve()->write(of);
				    }
				}
			      else
				{
				  crv->writeStandardHeader(of);
				  crv->write(of);
				}
			    }
#endif

			  // Check feasability
			  if (!sameFace(loop2))
			    {
			      if (!checkPlaneLoop(loop2))
				{
				  // Two candidate loops. Choose the one with smallest
				  // distance between edge nr 2 and 4
				  Point pos = loop2[1]->point(0.5*(loop2[1]->tMin() +
								   loop2[1]->tMax()));

				  Point cand1, cand2;
				  double par1, par2, dist1, dist2;
				  loop[3]->closestPoint(pos, par1, cand1, dist1);
				  loop2[3]->closestPoint(pos, par2, cand2, dist2);

				  if (dist2 < dist1)
				    loop = loop2;
				}
			    }
			}
		    }
		  return true;
		}
	    }
	}
      else
	{
	  // Find next edge in loop
	  bool found = getLoopEdges(loop, start_vx, vx2, max_nmb);
	  if (found)
	    {
	      // A loop is found. Check if a better alternative exists
	      for (size_t ki2=ki+1; ki2<edges.size(); ++ki2)
		{
		  vector<ftEdge*> loop2 = loop;
		  loop2.pop_back();
		  loop2[loop2.size()-1] = edges[ki2];
		  bool found2 = getLoopEdges(loop2, start_vx, vx2, max_nmb);
		  if (found2)
		    {
		      // Two candidate loops. Choose the one with smallest
		      // distance between edge nr 1 and 3
		      Point pos = loop2[0]->point(0.5*(loop2[0]->tMin() +
						       loop2[0]->tMax()));
		      
		      Point cand1, cand2;
		      double par1, par2, dist1, dist2;
		      loop[2]->closestPoint(pos, par1, cand1, dist1);
		      loop2[2]->closestPoint(pos, par2, cand2, dist2);
		      
		      if (dist2 < dist1)
			loop = loop2;
		    }
		      
		}
	      return found;
	    }
	}

      // No loop is found. Remove this instance
      loop.erase(loop.end() - 1);
      max_nmb -= added_max;
    }
  return false;
}

//===========================================================================
bool ftVolume::sameFace(vector<ftEdge*>& loop)
//===========================================================================
{
  // A loop is found. Check if a surface already exists
  size_t kj = 0;
  size_t nmb;  // The first loop edge may not belong to
  // any face, but the last one will
  vector<ftSurface*> faces;
  for (nmb=0; nmb<loop.size(); ++nmb)
    {
      vector<ftSurface*> faces2;
      if (loop[nmb]-> hasEdgeMultiplicity())
	{
	  faces2 = loop[nmb]->getEdgeMultiplicityInstance()->getAdjacentFaces(this);
	}
      else
	{
	  if (loop[nmb]->face())
	    {
	      faces2.push_back(loop[nmb]->face()->asFtSurface());
	      if (loop[nmb]->twin())
		faces2.push_back(loop[nmb]->twin()->face()->asFtSurface());
	    }
	}
      if (faces2.size() >= 1)
	{
	  faces.insert(faces.end(), faces2.begin(), faces2.end());
	  break;
	}
    }

  if (nmb >= loop.size())
    return false;  // No faces attached

  if (faces.size() >= 2)
    {
      for (kj=0; kj<loop.size(); ++kj)
	{
	  if (kj == nmb)
	    continue;

	  size_t kr, kh;
	  vector<ftSurface*> faces2;
	  if (loop[kj]-> hasEdgeMultiplicity())
	    {
	      faces2 = 
		loop[kj]->getEdgeMultiplicityInstance()->getAdjacentFaces(this);
	    }
	  else
	    {
	      if (loop[kj]->face())
		{
		  faces2.push_back(loop[kj]->face()->asFtSurface());
		  if (loop[kj]->twin())
		    faces2.push_back(loop[kj]->twin()->face()->asFtSurface());
		}
	    }
	  
	  for (kr=0; kr<faces.size(); ++kr)
	    {
	      for (kh=0; kh<faces2.size(); ++kh)
		if (faces[kr] == faces2[kh])
		  break;
	      if (kh < faces2.size())
		break;
	    }
	  
	  if (kr == faces.size())
	    break;  // The loop does not describe the same face
	}

      if (kj < loop.size())
	{
	  if (loop.size() > 4)
	    {
	      // Check for loop within loop
	      for (kj=0; kj<loop.size(); ++kj)
		{
		  size_t kr;
		  int kh;
		  vector<ftEdge*> loop2;
		  for (kr=kj, kh=0; kh<3; kr++, kh++)
		    {
		      if (kr == loop.size())
			kr = 0;
		      loop2.push_back(loop[kr]);
		    }

		  // Fetch the end vertices of the current path
		  shared_ptr<Vertex> vx_tmp = loop2[0]->getCommonVertex(loop2[1]);
		  if (!vx_tmp.get())
		    {
		      kj = loop.size();
		      break;
		    }
		  shared_ptr<Vertex> vx1 = loop2[0]->getOtherVertex(vx_tmp.get());
		  vx_tmp = loop2[loop2.size()-1]->getCommonVertex(loop2[loop2.size()-2]);
		  if (!vx_tmp.get())
		    {
		      kj = loop.size();
		      break;
		    }
		  shared_ptr<Vertex> vx2 = 
		    loop2[loop2.size()-1]->getOtherVertex(vx_tmp.get());
	      
		  // Check if there is a loop completing this edge
		  ftEdge *edg = vx1->getCommonEdge(vx2.get());
		  if (edg)
		    {
		      loop2.push_back(edg);
#ifdef DEBUG_VOL1
		      std::ofstream of("curr_edge_loop2.g2");
		      for (size_t k2=0; k2<loop2.size(); ++k2)
			{
			  shared_ptr<ParamCurve> cv = loop2[k2]->geomCurve();
			  shared_ptr<CurveOnSurface> cv2 = 
			    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
			  if (cv2.get())
			    {
			      cv2->spaceCurve()->writeStandardHeader(of);
			      cv2->spaceCurve()->write(of);
			    }
			  else
			    {
			      cv->writeStandardHeader(of);
			      cv->write(of);
			    }
			}
#endif
		      bool same = sameFace(loop2);
		      if (same)
			break;

		      // Check also planarity
		      bool plane = checkPlaneLoop(loop2);
		      if (plane)
			break;
		    }
		}
	      if (kj < loop.size())
		return true;
	      else
		return false;
	    }
	  else
	    return false;
	}
      else
	return true;
    }
	      
  return true;
}


//===========================================================================
bool ftVolume::checkPlaneLoop(vector<ftEdge*>& loop)
//===========================================================================
{
  int nmb_plane = 0;
  int nmb_check = 0;
  int nmb_free = 0;
  size_t kr, kh;
  for (kr=0; kr<loop.size(); ++kr)
    {
      kh = (kr == 0) ? loop.size()-1 : kr-1; 
      Point tan1 = loop[kh]->tangent(loop[kh]->tMin());
      Point tan2 = loop[kr]->tangent(loop[kr]->tMax());
      if (tan1.angle(tan2) < shells_[0]->getTolerances().bend)
	continue;  // Not a clear plane

      Point norm = tan1.cross(tan2);
      if (norm.length() < shells_[0]->getTolerances().gap)
	continue;

      if (!loop[kr]->face() || !loop[kr]->twin())
	{
	  nmb_free++;
	  continue;
	}
      if (!loop[kr]->twin()->face())
	{
	  nmb_free++;
	  continue;
	}

      ftSurface *f1 = loop[kr]->face()->asFtSurface();
      ftSurface *f2 = loop[kr]->twin()->geomEdge()->face()->asFtSurface();
      shared_ptr<Vertex> vx = loop[kr]->getCommonVertex(loop[kh]);
      if (!vx.get())
	continue;
      Point par1 = vx->getFacePar(f1);
      Point par2 = vx->getFacePar(f2);
      Point norm1 = f1->normal(par1[0], par1[1]);
      Point norm2 = f2->normal(par2[0], par2[1]);
      double ang1 = std::min(norm.angle(norm1), norm.angle(-norm1));
      double ang2 = std::min(norm.angle(norm2), norm.angle(-norm2));
      double ang = std::min(ang1, ang2);
      nmb_check++;
      if (ang < shells_[0]->getTolerances().bend)
	nmb_plane++;
    }
  
  if (nmb_plane < (int)(0.5*(nmb_check+nmb_free) + 1))
    return false;
  else
    return true;
}


//===========================================================================
double ftVolume::sortCoonsPatchBdCvs(vector<shared_ptr<ParamCurve> >& cvs,
				   vector<shared_ptr<ParamCurve> >& space_cvs,
				   Point& deg_pt, double tol)
//===========================================================================
{
  double max_cv_len = 0.0;

  // Check if the volume is closed in any direction
  int closed[3];
  shared_ptr<SplineVolume> vol2 = 
    dynamic_pointer_cast<SplineVolume, ParamVolume>(vol_);
  if (vol2.get())
    {
      for (int ki=0; ki<3; ++ki)
	closed[ki] = vol2->volumePeriodicity(ki, tol);
    }
  else
    closed[0] = closed[1] = closed[2] = -1;

  // Make sure to start with an existing curve
  size_t kj;
  size_t nmb = cvs.size();
  double len = 0.0;
  int nmb_len = 0;
  for (kj=0; kj<nmb;)
    {
      if (!cvs[kj].get())
	{
	  cvs.push_back(cvs[kj]);
	  cvs.erase(cvs.begin()+kj);
	  space_cvs.push_back(space_cvs[kj]);
	  space_cvs.erase(space_cvs.begin()+kj);
	  --nmb;
	}
      else
	{
	  break;
	}
    }
  if (deg_pt.dimension() > 0 && (!cvs[0].get()))
    {
      // No curves exists. Make degenerate curves
      double pdel = 2.0*toptol_.neighbour;
      for (kj=0; kj<cvs.size(); ++kj)
	{
	  cvs[kj] = shared_ptr<ParamCurve>(new SplineCurve(deg_pt, 0.0, 
							   deg_pt, pdel));
	}
      return max_cv_len;
    }
  else if (false) //nmb == 3 && deg_pt.dimension() > 0)
    {
      // Make sure that the position of the degenerate curve is
      // consistent with the given position 
      int deg_idx = -1;
      double mindist = std::numeric_limits<double>::max();
      size_t kh;
      for (kj=0; kj<cvs.size(); kj=kh)
	{
	  for (kh=kj+1; kh<cvs.size(); ++kh)
	    if (cvs[kh].get())
	      break;
	  if (kh == cvs.size())
	    kh = 0; 
	  Point pt1 = cvs[kj]->point(cvs[kj]->endparam());
	  Point pt2 = cvs[kh]->point(cvs[kh]->startparam());
	  double dist = deg_pt.dist(0.5*(pt1+pt2));
	  if (dist < mindist)
	    {
	      mindist = dist;
	      deg_idx = (int)kj;
	    }
	}
      if (deg_idx >= 0)
	{
	  int deg_idx2 = (deg_idx < (int)cvs.size()-1) ? deg_idx+1 : 0;
	  if (cvs[deg_idx2].get())
	    {
	      shared_ptr<ParamCurve> dummy;
	      cvs.insert(cvs.begin()+deg_idx2, dummy);
	      for (kh=0; kh<cvs.size(); ++kh)
		{
		  if ((int)kh != deg_idx2 && (!cvs[kh].get()))
		    {
		      cvs.erase(cvs.begin()+kh);
		      break;  // Only one degenerate curve
		    }
		}
	    }
	}
    }

  for (kj=0; kj<cvs.size(); ++kj)
    {
      if (cvs[kj].get())
	{
	  double curr_len = cvs[kj]->estimatedCurveLength();
	  max_cv_len = std::max(max_cv_len, curr_len);
	  len += curr_len;
	  nmb_len += 1;
	}  
    }
  if (nmb_len == 0)
    return max_cv_len;
  len /= (double)nmb_len;

  Point pos0 = cvs[0]->point(cvs[0]->startparam());
  Point pos1 = cvs[0]->point(cvs[0]->endparam());
  size_t kh;
  for (kj=1; kj<cvs.size(); ++kj)
    {
      kh = kj;
      while (!cvs[kh].get())
	{
	  ++kh;
	  kh = kh % cvs.size();
	}
      Point pos2 = cvs[kh]->point(cvs[kh]->startparam());
      Point pos3 = cvs[kh]->point(cvs[kh]->endparam());
      if (kj == 1)
	{
	  double d0 = std::min(pos0.dist(pos2), pos0.dist(pos3));
	  double d1 = std::min(pos1.dist(pos2), pos1.dist(pos3));
	  if (d0 > tol && d1 > tol)
	    {
	      // Check if the first curve follows a seem in the volume.
	      // but on the wrong side
	      Point dir(0.0, 0.0, 0.0);
	      int ki;
	      for (ki=0; ki<3; ++ki)
		{
		  if (closed[ki] < 0)
		    continue;
		  double t1 = vol2->startparam(ki);
		  double t2 = vol2->endparam(ki);
		  if (pos0[ki]-t1 < tol && pos1[ki]-t1 < tol &&
		      t2 - std::max(pos2[ki], pos3[ki]) < tol &&
		      t2 - std::min(pos2[ki], pos3[ki]) >= tol)
		    {
		      dir[ki] = t2 - t1;
		      break;
		    }
		  if (t2-pos0[ki] < tol && t2-pos1[ki] < tol &&
		      std::min(pos2[ki], pos3[ki]) - t1 < tol &&
		      std::max(pos2[ki], pos3[ki]) - t1 >= tol)
		    {
		      dir[ki] = t1 - t2;
		      break;
		    }
		}
	      if (ki < 3)
		{
		  moveVolParCv(cvs[0], space_cvs[0], dir, toptol_.gap);
		  pos0 = cvs[0]->point(cvs[0]->startparam());
		  pos1 = cvs[0]->point(cvs[0]->endparam());
		  d0 = std::min(pos0.dist(pos2), pos0.dist(pos3));
		  d1 = std::min(pos1.dist(pos2), pos1.dist(pos3));
		}
	      else
		{
		  // The situation where the first curve is wrongly
		  // oriented and the second follows the wrong seem is
		  // not covered. If this is a possibility, check the
		  // orientation of the first curve in geometry space
		  if ((closed[0] >= 0 || closed[1] >= 0 || closed[2] >= 0) &&
		      space_cvs[0].get() && space_cvs[1].get())
		    {
		      Point space0 = space_cvs[0]->point(space_cvs[0]->startparam());
		      Point space1 = space_cvs[0]->point(space_cvs[0]->endparam());
		      Point space2 = space_cvs[0+1]->point(space_cvs[0+1]->startparam());
		      Point space3 = space_cvs[0+1]->point(space_cvs[0+1]->endparam());
		      double d0_space = std::min(space0.dist(space2), space0.dist(space3));
		      double d1_space = std::min(space1.dist(space2), space1.dist(space3));
		      if (d0_space <= tol && d1_space > tol && d1 < d0)
			{
			  std::swap(d0, d1);
			}
		    }
		}
	    }	
	  if (d0 < d1)
	    {
	      // First curve has inconsistent orientation
	      cvs[0]->reverseParameterDirection();
	      if (space_cvs[0].get())
		space_cvs[0]->reverseParameterDirection();
	      std::swap(pos0, pos1);
	    }
	}
	      
      double d2 = pos1.dist(pos2);
      double d3 = pos1.dist(pos3);
      if (d2 > tol && d3 > tol)
	{
	  // Check if the current curve follows a seem in the volume.
	  // but on the wrong side
	  Point dir(0.0, 0.0, 0.0);
	  int ki;
	  for (ki=0; ki<3; ++ki)
	    {
	      if (closed[ki] < 0)
		continue;
	      double t1 = vol2->startparam(ki);
	      double t2 = vol2->endparam(ki);
	      if (pos2[ki]-t1 < tol && pos3[ki]-t1 < tol &&
		  t2 - pos1[ki] < tol)
		{
		  dir[ki] = t2 - t1;
		  break;
		}
	      if (t2-pos2[ki] < tol && t2-pos3[ki] < tol &&
		  pos1[ki] - t1 < tol)
		{
		  dir[ki] = t1 - t2;
		  break;
		}
	    }
	  if (ki < 3)
	    {
	      moveVolParCv(cvs[kh], space_cvs[kh], dir, toptol_.gap);
	      pos2 = cvs[kh]->point(cvs[kh]->startparam());
	      pos3 = cvs[kh]->point(cvs[kh]->endparam());
	    }
	}	

      if (pos1.dist(pos2) < pos1.dist(pos3))
	pos1 = pos3;
      else
	{
	  cvs[kh]->reverseParameterDirection();
	  if (space_cvs[kh].get())
	    space_cvs[kh]->reverseParameterDirection();
	  pos1 = pos2;
	}
      if (kh < kj)
	break;
      kj = kh;
    }

  // A last check for consistent curve orientation
  Point pt1 = cvs[0]->point(cvs[0]->endparam());
  for (kj=0; kj<cvs.size(); kj=kh)
    {
      for (kh=kj+1; kh<cvs.size(); ++kh)
	if (cvs[kh].get())
	  break;
      if (kh == cvs.size())
	break;
      Point pt2 = cvs[kh]->point(cvs[kh]->startparam());
      Point pt3 = cvs[kh]->point(cvs[kh]->endparam());
      if (pt1.dist(pt3) < pt1.dist(pt2))
	{
	  cvs[kh]->reverseParameterDirection();
	  if (space_cvs[kh].get())
	    space_cvs[kh]->reverseParameterDirection();
	  pt1 = pt2;
	}
      else
	pt1 = pt3;
    }

  for (kj=1; kj<cvs.size(); ++kj)
    {
      for (kh=kj; kh<cvs.size(); ++kh)
	if (cvs[kh].get())
	  break;
      if (kh == cvs.size())
	kh = 0; 

      // Check for missing curves
      if (kh != kj)
	{
	  size_t stop = (kj < kh) ? kh : cvs.size();
	  pt1 = cvs[kj-1]->point(cvs[kj-1]->endparam());
	  Point pt2 = cvs[kh]->point(cvs[kh]->startparam());
	  double del = pt1.dist(pt2);
	  del = std::max(std::max(del, 0.1*len), 1.0e-4);
	  for (size_t ka=kj; ka<stop; ++ka)
	    {
	      // if (deg_pt.dimension() > 0)
	      // 	cvs[ka] = shared_ptr<ParamCurve>(new SplineCurve(deg_pt, 0.0, 
	      // 							 deg_pt, del));
	      // else
		cvs[ka] = shared_ptr<ParamCurve>(new SplineCurve(pt1, 0.0, 
								 pt2, del));
	    }
	}
    }
  return max_cv_len;
}

//===========================================================================
// 
// 
void
ftVolume::moveVolParCv(shared_ptr<ParamCurve>& pcv,
		       shared_ptr<ParamCurve>& spacecv,
		       const Point& dir, double tol)
//===========================================================================
{
  shared_ptr<SplineCurve> crv = 
    dynamic_pointer_cast<SplineCurve,ParamCurve>(pcv);
  if (crv.get())
    crv->translateCurve(dir);
  else
    {
      shared_ptr<SplineCurve> crv2;
      if (spacecv.get())
	{
	  Point pos1 = pcv->point(pcv->startparam());
	  Point pos2 = pcv->point(pcv->endparam());
	  pos1 += dir;
	  pos2 += dir;
	  shared_ptr<Point> par1(new Point(pos1.begin(), pos1.end()));
	  shared_ptr<Point> par2(new Point(pos2.begin(), pos2.end()));
	  shared_ptr<VolumeParameterCurve> vol_par =
	    shared_ptr<VolumeParameterCurve>(new VolumeParameterCurve(vol_, 
								      spacecv,
								      par1, par2));
	  // Approximate
	  HermiteAppC approx(vol_par.get(), tol, tol);
	  approx.refineApproximation();
	  crv2 = approx.getCurve();
	}
      if (crv2.get())
	pcv = crv2;
      else
	{
	  crv = shared_ptr<SplineCurve>(pcv->geometryCurve());
	  crv->translateCurve(dir);
	  pcv = crv;
	}
    }
}

//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> > 
ftVolume::createRegularVolumes(vector<shared_ptr<ftSurface> > bd_faces)
//===========================================================================
{
  // Fetch all faces
  vector<shared_ptr<ftSurface> > faces;
  int nmb_shells = nmbOfShells();
  for (int ki=0; ki<nmb_shells; ++ki)
    {
      // Fetch boundary shell
      shared_ptr<SurfaceModel> sfmodel = getShell(ki);
      vector<shared_ptr<ftSurface> > curr_faces = sfmodel->allFaces();
      faces.insert(faces.end(), curr_faces.begin(), curr_faces.end());
    }
  faces.insert(faces.end(), bd_faces.begin(), bd_faces.end());

  // Collect sets of connected faces
  vector<shared_ptr<SurfaceModel> > models;
  vector<shared_ptr<ftSurface> > curr_set;
  vector<shared_ptr<ftSurface> > all_sets;
  for (size_t kj=0; kj<faces.size(); )
    {
      curr_set.clear();

      // Store first face and all connected faces
      getCurrConnectedModel(faces, kj, curr_set, all_sets);
      
      // Make surface model
	shared_ptr<SurfaceModel> curr_model = 
	  shared_ptr<SurfaceModel>(new SurfaceModel(toptol_.gap, 
						    toptol_.gap,
						    toptol_.neighbour,
						    toptol_.kink,
						    toptol_.bend,
						    curr_set, true));
	models.push_back(curr_model);

#ifdef DEBUG_VOL1
	std::ofstream of("regvol2.g2");
	int nmb = curr_model->nmbEntities();
	for (int kv=0; kv<nmb; ++kv)
	  {
	    shared_ptr<ParamSurface> sf = curr_model->getSurface(kv);
	    sf->writeStandardHeader(of);
	    sf->write(of);
	  }
	vector<shared_ptr<Vertex> > vx;
	curr_model->getAllVertices(vx);
	of << "400 1 0 4 255 0 0 255" << std::endl;
	of << vx.size() << std::endl;
	for (size_t ka=0; ka<vx.size(); ++ka)
	  of << vx[ka]->getVertexPoint() << std::endl;
#endif

	// Remove used faces
	for (size_t kh=0; kh<curr_set.size(); ++kh)
	  {
	    size_t kr;
	    for (kr=0; kr<faces.size(); ++kr)
	      {
		if (faces[kr].get() == curr_set[kh].get())
		  {
		    faces.erase(faces.begin()+kr);
		    break;
		  }
	      }
	  }
    }

  // Make new trimmed volumes
  vector<shared_ptr<ftVolume> > result;
  for (size_t kj=0; kj<models.size(); ++kj)
    {
      shared_ptr<ftVolume> tmp_vol(new ftVolume(vol_, models[kj]));
      if (hasMaterialInfo())
	tmp_vol->setMaterial(getMaterial());
      result.push_back(tmp_vol);
    }

  return result;
}
		  
//===========================================================================
void 
ftVolume::getCurrConnectedModel(vector<shared_ptr<ftSurface> >& faces,
				size_t idx,
				vector<shared_ptr<ftSurface> >& curr_set,
				vector<shared_ptr<ftSurface> >& all_sets) const
//===========================================================================
{
  // Store current face
  curr_set.push_back(faces[idx]);
  all_sets.push_back(faces[idx]);

#ifdef DEBUG_VOL1
  std::ofstream of("curr_conn.g2");
  for (size_t kv=0; kv<curr_set.size(); ++kv)
    {
      shared_ptr<ParamSurface> sf = curr_set[kv]->surface();
      sf->writeStandardHeader(of);
      sf->write(of);
    }
#endif

  // Fetch all neighbours
  vector<ftSurface*> neighbours;
  faces[idx]->getAdjacentFaces(neighbours);
  
  // For all neighbours, store the neighbour and all its neighbours
  // as long as they are not found already
  for (size_t ki=0; ki<neighbours.size(); ki++)
    {
      size_t kj;
      for (kj=0; kj<all_sets.size(); ++kj)
	if (neighbours[ki] == all_sets[kj].get())
	  break;
      
      if (kj < all_sets.size())
	continue;  // Face categorized already
      
      // Handle neighbours to neighbour
      for (kj=0; kj<faces.size(); ++kj)
	if (faces[kj].get() == neighbours[ki])
	  break;
      if (kj < faces.size())
	getCurrConnectedModel(faces, kj, curr_set, all_sets);
    }
}
  
       
//===========================================================================
void 
ftVolume::replaceParamVolume(shared_ptr<ParamVolume> vol, 
			     vector<shared_ptr<ParamSurface> >& sorted_sfs,
			     bool loft_sequence)
//===========================================================================
{
  // Set new parametric volume
  vol_ = vol;

  int nmb_sfs = 0;
  for (size_t ki=0; ki<sorted_sfs.size(); ++ki)
    if (sorted_sfs[ki].get())
      nmb_sfs++;

  // Fetch new boundary faces
  vector<shared_ptr<ftSurface> > bd_faces = 
    getBoundaryFaces(vol, toptol_.neighbour /*toptol_.gap*/, toptol_.kink);
  if ((int)bd_faces.size() < nmb_sfs)
    bd_faces = 
      getBoundaryFaces(vol, toptol_.gap, toptol_.kink);

#ifdef DEBUG_VOL1
  std::ofstream of0("bd_sfs.g2");
  for (size_t ki=0; ki<bd_faces.size(); ++ki)
    {
      bd_faces[ki]->surface()->writeStandardHeader(of0);
      bd_faces[ki]->surface()->write(of0);
    }
#endif

  // Replace the faces in the boundary shell, one by one, reset
  // twin information
#ifdef DEBUG_VOL1
  std::ofstream of1("curr_sfs.g2");
#endif
  if (nmb_sfs < (int)bd_faces.size())
    {
      vector<double> face_size(bd_faces.size(), 0.0);
      for (size_t ki=0; ki<bd_faces.size(); ++ki)
	{
	  double u_size, v_size;
	  bd_faces[ki]->surface()->estimateSfSize(u_size, v_size);
	  face_size[ki] = u_size*v_size;
	}
      while ((int)bd_faces.size() > nmb_sfs)
	{
	  int min_ind = 0;
	  double min_size = face_size[min_ind];
	  for (size_t ki=1; ki<bd_faces.size(); ++ki)
	    if (face_size[ki] < min_size)
	      {
		min_ind = (int)ki;
		min_size = face_size[ki];
	      }
	  bd_faces.erase(bd_faces.begin()+min_ind);
	  face_size.erase(face_size.begin()+min_ind);
	}
      int stop_break = 1;
    }
  int kj = 0;
  for (size_t ki=0; ki<sorted_sfs.size(); ++ki)
    {
      if (!sorted_sfs[ki].get())
	continue;

      // Fetch face pointer in shell
      int idx = shells_[0]->getIndex(sorted_sfs[ki].get());
      shared_ptr<ftSurface> face = shells_[0]->getFace(idx);

      //int idx2 = findFaceMatch(face, bd_faces);
      int idx2 = kj;
      if (loft_sequence)
	{
	  if (ki == 2 || ki == 3)
	    idx2 += 2;
	  if (ki == 4 || ki == 5)
	      idx2 = (int)ki - 2;
	}

#ifdef DEBUG_VOL1
      // TEST
      std::ofstream of("corner_vx.g2");
      vector<pair<Point,Point> > corner1;
      vector<pair<Point,Point> > corner2;
      sorted_sfs[ki]->getCornerPoints(corner1);
      bd_faces[idx2]->surface()->getCornerPoints(corner2);
      of << "400 1 0 4 255 0 0 255" << std::endl;
      of << corner1.size() << std::endl;
      size_t kv;
      for (kv=0; kv<corner1.size(); ++kv)
	of << corner1[kv].first << std::endl;
      of << "400 1 0 4 0 255 0 255" << std::endl;
      of << corner2.size() << std::endl;
      for (kv=0; kv<corner2.size(); ++kv)
	of << corner2[kv].first << std::endl;

      sorted_sfs[ki]->writeStandardHeader(of1);
      sorted_sfs[ki]->write(of1);
      bd_faces[idx2]->surface()->writeStandardHeader(of1);
      bd_faces[idx2]->surface()->write(of1);
      face->surface()->writeStandardHeader(of1);
      face->surface()->write(of1);
#endif

      // Set body and boundary information
      bd_faces[idx2]->setBody(this);
      if (face->hasBoundaryConditions())
	{
	  int bd_type, bd;
	  face->getBoundaryConditions(bd_type, bd);
	  bd_faces[idx2]->setBoundaryConditions(bd_type, bd);
	}

      // Replace current face
      // First assume simple configuration
      // Remember twin
      ftSurface *twin = face->twin();
      face->disconnectTwin();

#if 0
      // Not quality assured. There are still problems
      bool replaced = false; //shells_[0]->replaceFace(face, bd_faces[idx2]);
      if (!replaced)
	{
#endif	  
	  shells_[0]->removeFace(face);
	  shells_[0]->append(bd_faces[idx2], false);
	  if (twin)
	    {
	      bd_faces[idx2]->connectTwin(twin, toptol_.neighbour);
#ifdef DEBUG_VOL1
	      bool isOK = bd_faces[idx2]->checkFaceTopology();
	      if (!isOK)
		std::cout << "Inconsistencies in connect twin " << std::endl;
#endif
	    }
#if 0
	}
#endif
      ++kj;
    }
  int stop_break;
  stop_break = 1;
}

//===========================================================================
void 
ftVolume::updateBoundaryInfo()
//===========================================================================
{
  // This function should be made more efficient to avoid a lot of topology
  // analysis
  // Fetch new boundary faces
  vector<shared_ptr<ftSurface> > bd_faces = 
    getBoundaryFaces(vol_, toptol_.gap, toptol_.kink);

  // Replace the faces in the boundary shell, one by one, reset
  // twin information
  vector<shared_ptr<ftSurface> > faces = shells_[0]->allFaces();
  size_t nmb = faces.size();
  if (nmb != bd_faces.size())
    return;  // Different number of faces

  for (size_t ki=0; ki<nmb; ++ki)
    {
      // Fetch face pointer in shell
      shared_ptr<ftSurface> face = faces[ki];

      // Find closest boundary face. Evaluate internal point
      double u1, v1;
      Point pos = face->surface()->getInternalPoint(u1, v1);

      int idx = -1;
      double min_dist = 1.0e8;
      for (size_t kj=0; kj<bd_faces.size(); ++kj)
	{
	  double u2, v2, dist;
	  Point pos2;
	  RectDomain *dom = NULL;
	  double seed[2];
	  seed[0] = u1;
	  seed[1] = v1;
	  bd_faces[kj]->closestPoint(pos, u2, v2, pos2, dist, toptol_.gap, dom, seed);
	  if (dist < min_dist)
	    {
	      min_dist = dist;
	      idx = (int)kj;
	    }
	}

      if (idx < 0)
	return;  // Something strange happened
	  
      // Remember twin
      ftSurface *twin = face->twin();
      face->disconnectTwin();

 #ifdef DEBUG_VOL1
     // TEST
      std::ofstream of1("curr_sfs.g2");
      bd_faces[idx]->surface()->writeStandardHeader(of1);
      bd_faces[idx]->surface()->write(of1);
      face->surface()->writeStandardHeader(of1);
      face->surface()->write(of1);
#endif

      // Replace current face
      shells_[0]->removeFace(face);
      bd_faces[idx]->setBody(this);
      if (face->hasBoundaryConditions())
	{
	  int bd_type, bd;
	  face->getBoundaryConditions(bd_type, bd);
	  bd_faces[idx]->setBoundaryConditions(bd_type, bd);
	}
      shells_[0]->append(bd_faces[idx], false);
      if (twin)
	bd_faces[idx]->connectTwin(twin, toptol_.neighbour);

      bd_faces.erase(bd_faces.begin()+idx);  // Face used

      int stop_break;
      stop_break = 1;
    }
}

//===========================================================================
int ftVolume::findFaceMatch(shared_ptr<ftSurface> face,
			    vector<shared_ptr<ftSurface> >& cand_matches)

//===========================================================================
{
  vector<pair<Point,Point> > corner1;
  face->surface()->getCornerPoints(corner1);

  double min_dist = 1.0e8;  // A large number
  int min_idx = -1;
  for (size_t ki=0; ki<cand_matches.size(); ++ki)
    {
      vector<pair<Point,Point> > corner2;
      cand_matches[ki]->surface()->getCornerPoints(corner2);
  
      double dist = 0.0;
      for (size_t kj=0; kj<corner1.size(); ++kj)
	{
	  double mind = 1.0e8;
	  for (size_t kr=0; kr<corner2.size(); ++kr)
	    {
	      double d = corner1[kj].first.dist(corner2[kr].first);
	      mind = std::min(mind, d);
	    }
	  dist += mind;
	}

      if (dist < min_dist)
	{
	  min_dist = dist;
	  min_idx = (int)ki;
	}
    }
  return min_idx;
}

//===========================================================================
vector<pair<int, int> >  ftVolume::oppositeSfs(shared_ptr<SurfaceModel> model)
//===========================================================================
{
  double ang_tol = M_PI/6.0;
  vector<pair<int, int> > opposite;
  int nmb = model->nmbEntities();
  int ki, kj;
  for (ki=0; ki<nmb; ++ki)
    {
      // Check if the surface already belongs to a match
      size_t kh;
      for (kh=0; kh<opposite.size(); ++kh)
	if (opposite[kh].first == ki || opposite[kh].second == ki)
	  break;
      if (kh < opposite.size())
	continue;

      int bd1 = -1, bd2 = -1;
      int orientation;
      bool swap;
      shared_ptr<ParamSurface> surf = model->getSurface(ki);

      // Check type
      shared_ptr<SurfaceOnVolume> vol_sf = 
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf);
      if (!vol_sf.get())
	{
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
	  if (bd_sf.get())
	    {
	      shared_ptr<ParamSurface> surf2 = bd_sf->underlyingSurface();
	      vol_sf = 
		dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf2);
	    }
	}

      if (vol_sf.get())
	bd1 = vol_sf->whichBoundary(model->getTolerances().gap, 
				    orientation, swap);

      if (bd1 >= 0)
	{
	  for (kj=ki+1; kj<nmb; ++kj)
	    {
	      bd2 = -1;
	      shared_ptr<ParamSurface> surf3 = model->getSurface(kj);

	      // Check type
	      shared_ptr<SurfaceOnVolume> vol_sf2 = 
		dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf3);
	      if (!vol_sf2.get())
		{
		  shared_ptr<BoundedSurface> bd_sf2 = 
		    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf3);
		  if (bd_sf2.get())
		    {
		      shared_ptr<ParamSurface> surf4 = 
			bd_sf2->underlyingSurface();
		      vol_sf2 = 
			dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf4);
		    }
		}

	      if (vol_sf2.get())
		bd2 = vol_sf2->whichBoundary(model->getTolerances().gap, 
					     orientation, swap);

	      if (bd2 >= 0 &&
		  bd1/2 == bd2/2 &&
		  abs((bd1%2) - (bd2%2)) == 1)
		{
		  opposite.push_back(make_pair(ki,kj));
		  break;
		}
	    }

	  if (kj == nmb)
	    {
	      // Find a non-boundary, non-adjacent surface with
	      // existing parameter surface and check if it has
	      // a normal cone consistent with being an opposite surface
	      int kr;
	      for (kr=0; kr<nmb; ++kr)
		{
		  if (kr == ki)
		    continue;
		  
		  bd2 = -1;
		  shared_ptr<ParamSurface> surf3 = model->getSurface(kr);

		  // Check type
		  shared_ptr<SurfaceOnVolume> vol_sf2 = 
		    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf3);
		  if (!vol_sf2.get())
		    {
		      shared_ptr<BoundedSurface> bd_sf2 = 
			dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf3);
		      if (bd_sf2.get())
			{
			  shared_ptr<ParamSurface> surf4 = 
			    bd_sf2->underlyingSurface();
			  vol_sf2 = 
			    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf4);
			}
		    }
		  
		  if (vol_sf2.get())
		    bd2 = vol_sf2->whichBoundary(model->getTolerances().gap, 
						 orientation, swap);

		  if (vol_sf2.get() == NULL || bd2 >= 0)
		    continue;  // Either not a surface on volume or a
		  // boundary surface

		  int history = vol_sf2->getCreationHistory();
		  //if (!vol_sf2->hasParameterSurface())
		  if (history != 1) 
		    continue;   // Not an internal splitting surface

		  // Check adjacency
		  shared_ptr<ftSurface> face1 = model->getFace(ki);
		  shared_ptr<ftSurface> face2 = model->getFace(kr);
		  bool smooth;
		  bool adjacent = face1->isAdjacent(face2.get(), smooth);
		  if (adjacent)
		    continue;
		  
		  // Check normal cone configuration
		  DirectionCone cone1 = surf->normalCone();
		  DirectionCone cone2 = surf3->normalCone();
		  double ang = cone1.centre().angle(cone2.centre());
		  ang = M_PI-ang;
		  if (ang < ang_tol)
		    {
		      opposite.push_back(make_pair(ki,kr));
		      break;
		    }
		}

	      if (kr == nmb)
		{
		  // Finally, add empty correspondance relattion to mark 
		  // the element boundary side
		  opposite.push_back(make_pair(ki,-1));
		}
	    }
	}
    }

  
  return opposite;
}

//===========================================================================
void  ftVolume::removeSeamFaces()
//===========================================================================
{
  int ki, kj;
  int nmb1 = nmbOfShells();
  bool updated = true;
  while (updated)
    {
      updated = false;
      for (ki=0; ki<nmb1; ++ki)
	{
	  shared_ptr<SurfaceModel> shell = getShell(ki);
	  int nmb2 = shell->nmbEntities();
	  for (kj=0; kj<nmb2; ++kj)
	    {
	      shared_ptr<ftSurface> face = shell->getFace(kj);
	      if (face->twin() && face->twin()->getBody() == this)
		{
		  ftSurface* face2 = face->twin()->asFtSurface();

		  // A candidate set of seam surfaces is found
#ifdef DEBUG_VOL1
		  std::ofstream of("seam_twin.g2");
		  shared_ptr<ParamSurface> sf1 = face->surface();
		  shared_ptr<ParamSurface> sf2 = face2->surface();

		  sf1->writeStandardHeader(of);
		  sf1->write(of);
		  sf2->writeStandardHeader(of);
		  sf2->write(of);
#endif
		  vector<ftSurface*> neighbours1;
		  vector<ftSurface*> neighbours2;
		  face->getAdjacentFaces(neighbours1);
		  face2->getAdjacentFaces(neighbours2);

		  if ((neighbours1.size() == 0 ||
		       (neighbours1.size() == 1 && neighbours1[0] == face2)) &&
		      (neighbours2.size() == 0 ||
		       (neighbours2.size() == 1 && neighbours2[0] == face.get())))
		    {
		      // Seam faces with no other connections has been found
		      // and can be removed
		      updated = true;

		      bool remove1;
		      remove1 = shell->removeFace(face);

		      shared_ptr<SurfaceModel> shell2 = getShell(face2);
		      shared_ptr<ftSurface> face2_2 = 
			shell2->fetchAsSharedPtr(face2);
		      bool remove2;
		      remove2 = shell2->removeFace(face2_2);
		    }
		}
	      if (updated)
		break;
	    }
#ifdef DEBUG_VOL1
	  std::ofstream of2("post_seam.g2");
	  for (kj=0; kj<shell->nmbEntities(); ++kj)
	    {
	      shared_ptr<ParamSurface> sf = shell->getSurface(kj);
	      sf->writeStandardHeader(of2);
	      sf->write(of2);
	    }
#endif

	  if (updated)
	    break;
	}
    }
}

//===========================================================================
// 
// 
void ftVolume::eraseMissingEdges()
//===========================================================================
{
  for (size_t ki=0; ki<missing_edges_.size(); ++ki)
    {
      shared_ptr<ftEdge> curr = missing_edges_[ki];
      shared_ptr<Vertex> v1, v2;
      curr->getVertices(v1, v2);

      // Remove edge from associated vertices
      v1->removeEdge(curr.get());
      v2->removeEdge(curr.get());
    }
  missing_edges_.clear();
}

//===========================================================================
shared_ptr<ParamCurve> ftVolume::makeMissingEdgeCv(shared_ptr<Vertex> vx1,
						   shared_ptr<Vertex> vx2)
//===========================================================================
{
  shared_ptr<ParamCurve> crv;

  // Check the tangent direction for a linear curve agains all surface
  // normals in the vertices
  size_t ki;
  Point vx1_pt = vx1->getVertexPoint();
  Point vx2_pt = vx2->getVertexPoint();
  Point vec = vx2_pt - vx1_pt;
  double len = vec.length();
  vec.normalize();
  double ang_tol = 0.05*M_PI;
  Point d1 = vec;
  Point d2 = -vec;

  // 1. vertex
  vector<pair<ftSurface*, Point> > face_par1 = vx1->getFaces(this);
  Point norm1(0.0, 0.0, 0.0);
  bool linear1 = true;
  for (ki=0; ki<face_par1.size(); ++ki)
    {
      Point normal = face_par1[ki].first->normal(face_par1[ki].second[0],
						 face_par1[ki].second[1]);
      double ang = vec.angle(normal);
      ang = fabs(0.5*M_PI - ang);
      if (ang < ang_tol)
	{
	  linear1 = false;
	  normal.normalize();
	  norm1 += normal;
	}
    }
  if (linear1)
    d1 = vec;
  else
    {
      d1.normalize();
      norm1.normalize();
      norm1 *= -1;
      d1 = (2.0*d1 + norm1)/3.0;
    }
    

  // 2. vertex
  vector<pair<ftSurface*, Point> > face_par2 = vx2->getFaces(this);
  Point norm2(0.0, 0.0, 0.0);
  bool linear2 = true;
  for (ki=0; ki<face_par2.size(); ++ki)
    {
      Point normal = face_par2[ki].first->normal(face_par2[ki].second[0],
						 face_par2[ki].second[1]);
      double ang = vec.angle(normal);
      ang = fabs(0.5*M_PI - ang);
      if (ang < ang_tol)
	{
	  linear2 = false;
	  normal.normalize();
	  norm2 += normal;
	}
    }

  if (linear2)
    d2 = vec;
  else
    {
      d2.normalize();
      norm2.normalize();
      norm2 *= -1;
      d2 = (2.0*d2 + norm2)/3.0;
    }
    

  if (linear1 && linear2)
    {
      crv = shared_ptr<ParamCurve>(new SplineCurve(vx1_pt, vx2_pt));
    }
  else
    {
      // Set length of tangents
      //double len_fac = 0.1;
      d1.normalize();
      //d1 *= len_fac*len;
      d2.normalize();
      //d2 *= -len_fac*len;
      d2 *= -1;

      HermiteInterpolator intpol;
      vector<Point> data(4);
      vector<double> param(2);
      param[0] = 0.0;
      param[1] = len;
      data[0] = vx1_pt;
      data[1] = d1;
      data[2] = vx2_pt;
      data[3] = d2;

      vector<double> coefs;
      intpol.interpolate(data, param, coefs);
      BsplineBasis basis = intpol.basis();
      crv = shared_ptr<ParamCurve>(new SplineCurve(basis, coefs.begin(), 3));
    }
  return crv;
}

//===========================================================================
bool ftVolume::removeSliver1(shared_ptr<ftSurface> face, 
			     vector<shared_ptr<ftEdge> >& edg, 
			     int ix1, int ix2, double tol,
			     vector<shared_ptr<ParamSurface> >& mod_sfs)
//===========================================================================
{
  bool modified = false;
  if (edg.size() < 3 || edg.size() > 4)
    return modified;  // Not correct configuration

  // Surfaces to extend to remove sliver face
  ftSurface* next_faces[2];
  shared_ptr<ParamSurface> sf[2];
  shared_ptr<BoundedSurface> bdsf[2];
  shared_ptr<SurfaceOnVolume> volsf[2];
  next_faces[0] = edg[ix1]->twin()->face()->asFtSurface();
  next_faces[1] = edg[ix2]->twin()->face()->asFtSurface();
  sf[0] = next_faces[0]->surface();
  sf[1] = next_faces[1]->surface();
  bdsf[0] = dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf[0]);
  bdsf[1] = dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf[1]);
  if (bdsf[0].get() == NULL || bdsf[1].get() == NULL)
    return modified;
  volsf[0] = 
    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bdsf[0]->underlyingSurface());
  volsf[1] =
    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bdsf[1]->underlyingSurface());
  if (volsf[0].get() == NULL || volsf[1].get() == NULL)
    return modified;
#ifdef DEBUG
  std::ofstream of0("adj_sliver.g2");
  sf[0]->writeStandardHeader(of0);
  sf[0]->write(of0);
  sf[1]->writeStandardHeader(of0);
  sf[1]->write(of0);
#endif

  // Adjacent surfaces that must be modified
  vector<ftSurface*> adj_face(edg.size()-2);
  vector<shared_ptr<ParamSurface> > adj_sf(edg.size()-2);
  vector<shared_ptr<BoundedSurface> > adj_bdsf(edg.size()-2);
  vector<shared_ptr<SurfaceOnVolume> > adj_volsf(edg.size()-2);
  int ki, kj;
  for (ki=0, kj=0; ki<(int)edg.size(); ++ki)
    {
      if (ki == ix1 || ki == ix2)
	continue;
      adj_face[kj] = edg[ki]->twin()->face()->asFtSurface();
      adj_sf[kj] = adj_face[kj]->surface();
      adj_bdsf[kj] = dynamic_pointer_cast<BoundedSurface,ParamSurface>(adj_sf[kj]);
      if (adj_bdsf[kj].get() == NULL)
	return modified;
      adj_volsf[kj] = 
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(adj_bdsf[kj]->underlyingSurface());
      ++kj;
    }

  int bd[2];
  vector<int> adj_bd(edg.size()-2);
  int orientation;
  bool swap;
  bd[0] = volsf[0]->whichBoundary(tol, orientation, swap);
  bd[1] = volsf[1]->whichBoundary(tol, orientation, swap);

  // Only sliver faces that trim away a constant curve edge in the volume are handled
  if (bd[0] < 0 || bd[1] < 0)
    return modified;

  for (ki=0; ki<(int)adj_volsf.size(); ++ki)
    {
      adj_bd[ki] = adj_volsf[ki]->whichBoundary(tol, orientation, swap);
      // if (adj_bd[ki] < 0)
      // 	return modified;
    }

  // Identify new boundary edge
  int minbd = std::min(bd[0], bd[1]);
  int maxbd = std::max(bd[0], bd[1]);
  int bd_ix = 2*(maxbd-2) + minbd + 2*(minbd>=2);

  // Spline volume
  shared_ptr<SplineVolume> vol = dynamic_pointer_cast<SplineVolume,ParamVolume>(vol_);
  if (!vol.get())
    return modified;

  shared_ptr<ParamCurve> geomcv, parcv;
  vol->getBoundaryCurve(bd_ix, geomcv, parcv);

  // Restrict curve according to adjacent surfaces
  vector<double> int_par;
  shared_ptr<SplineCurve> tmp_cv(geomcv->geometryCurve());
  vector<shared_ptr<ftSurface> > faces = getOuterShell()->allFaces();
  if (tmp_cv.get())
    {
      for (ki=0; ki<(int)faces.size(); ++ki)
	{
	  if (faces[ki].get() == face.get())
	    continue;
	  shared_ptr<ParamSurface> curr = faces[ki]->surface();
	  if (curr.get() == sf[0].get() || curr.get() == sf[1].get())
	    continue;
	  SplineSurface *tmp_sf = curr->getSplineSurface();
	  if (tmp_sf)
	    {
	      vector<pair<double, Point> > int_pts;
	      vector<int> pretop;
	      vector<pair<pair<double,Point>, pair<double,Point> > > int_cvs;
	      intersectCurveSurf(tmp_cv.get(), tmp_sf, tol,
				 int_pts, pretop, int_cvs);
	      for (kj=0; kj<(int)int_pts.size(); ++kj)
		int_par.push_back(int_pts[kj].first);
	    }
	}
    }

  // VSK.102017. Should really check which curve parts are inside, but it
  // is not obvious inside of what. Takes an easy solution 
  if (int_par.size() == 2)
    {
      if (int_par[0] > int_par[1])
	std::swap(int_par[0], int_par[1]);
      if (int_par[1]-int_par[0] > tol &&  
	  (int_par[0] > tmp_cv->startparam() || 
	   int_par[1] < tmp_cv->endparam()))
	{
	  shared_ptr<ParamCurve> geomcv2(geomcv->subCurve(int_par[0],int_par[1]));
	  shared_ptr<ParamCurve> parcv2(parcv->subCurve(int_par[0],int_par[1]));
	  geomcv = geomcv2;
	  parcv = parcv2;
	}
    }
  shared_ptr<CurveOnVolume> bd_cv(new CurveOnVolume(vol_, parcv, geomcv, false));

#ifdef DEBUG
  std::ofstream of1("sliver_cv.g2");
  geomcv->writeStandardHeader(of1);
  geomcv->write(of1);
#endif

  // Modified curve between extended surfaces and adjacent surfaces
  Point bd_pos[2];
  bd_pos[0] = bd_cv->ParamCurve::point(bd_cv->startparam());
  bd_pos[1] = bd_cv->ParamCurve::point(bd_cv->endparam());
  vector<shared_ptr<ParamCurve> > mod_all1, mod_all2;
  mod_all1.push_back(bd_cv);
  mod_all2.push_back(shared_ptr<ParamCurve>(bd_cv->clone()));
  for (ki=0; ki<(int)adj_bd.size(); ++ki)
    {
      vector<shared_ptr<ParamCurve> > mod_cvs(2);
      for (kj=0; kj<2; ++kj)
	{
	  if (adj_bd[ki] < 0)
	    {
	      // Extract boundary curve between surfaces
	      if (adj_face[ki]->nmbAdjacencies(next_faces[kj]) != 1)
		return modified;
	      shared_ptr<ftEdge> edg1, edg2;
	      (void)adj_face[ki]->areNeighbours(next_faces[kj], 
						edg1, edg2);
	      shared_ptr<ParamCurve> tmp_cv = edg1->geomCurve();
	      shared_ptr<CurveOnSurface> tmp_sfcv = 
		dynamic_pointer_cast<CurveOnSurface, ParamCurve>(tmp_cv);
	      if (tmp_sfcv.get())
		tmp_cv = tmp_sfcv->spaceCurve();
	  
#ifdef DEBUG
	      tmp_cv->writeStandardHeader(of1);
	      tmp_cv->write(of1);
#endif
	      // Define curve piece to new endpoint
	      // First identify closest end points
	      Point pos1 = tmp_cv->point(tmp_cv->startparam());
	      Point pos2 = tmp_cv->point(tmp_cv->endparam());
	      double dd1 = std::min(pos1.dist(bd_pos[0]), pos1.dist(bd_pos[1]));
	      double dd2 = std::min(pos2.dist(bd_pos[0]), pos2.dist(bd_pos[1]));
	      shared_ptr<SplineCurve> crv_piece;
	      int degree = 2;  
	      if (dd1 < dd2)
		{
		  vector<Point> pts(3);
		  vector<Point> der(2);
		  tmp_cv->point(der, tmp_cv->startparam(), 1);
		  pts[0] = (pos1.dist(bd_pos[0]) < pos1.dist(bd_pos[1])) ?
		    bd_pos[0] : bd_pos[1];
		  pts[1] = der[1];
		  pts[2] = der[0];
		  vector<int> type(3);
		  type[0] = type[2] = 1;
		  type[1] = 3;
		  crv_piece = SISLCurveInterface::interpolate(pts, type, degree);
		}
	      else
		{
		  vector<Point> pts(3);
		  tmp_cv->point(pts, tmp_cv->endparam(), 1);
		  pts[2] = (pos1.dist(bd_pos[0]) < pos1.dist(bd_pos[1])) ?
		    bd_pos[0] : bd_pos[1];
		  vector<int> type(3);
		  type[0] = type[2] = 1;
		  type[1] = 2;
		  crv_piece = SISLCurveInterface::interpolate(pts, type, degree);
		}				 

#ifdef DEBUG
	      crv_piece->writeStandardHeader(of1);
	      crv_piece->write(of1);
#endif
	      // Join curve pieces
	      shared_ptr<ParamCurve> crv_piece2;
	      if (tmp_cv->instanceType() == Class_CurveOnSurface)
		crv_piece2 = shared_ptr<ParamCurve>(new CurveOnVolume(vol_,
								      crv_piece,
								      false));
	      else
		crv_piece2 = crv_piece;
	      double dist;
	      if (dd1 < dd2)
		{
		  crv_piece2->appendCurve(tmp_cv.get(), 1, dist);
		  mod_cvs[kj] = crv_piece2;
		}
	      else
		{
		  mod_cvs[kj] = shared_ptr<ParamCurve>(tmp_cv->clone());
		  mod_cvs[kj]->appendCurve(crv_piece2.get(), 1, dist);
		}
	      if (!mod_cvs[kj].get())
		return modified;
	    }
	  else
	    {
	      minbd = std::min(adj_bd[ki], bd[kj]);
	      maxbd = std::max(adj_bd[ki], bd[kj]);
	      bd_ix = 2*(maxbd-2) + minbd + 2*(minbd>=2);
	      shared_ptr<ParamCurve> geomcv2, parcv2;
	      vol->getBoundaryCurve(bd_ix, geomcv2, parcv2);
	      mod_cvs[kj] = 
		shared_ptr<ParamCurve>(new CurveOnVolume(vol_, parcv2, 
							 geomcv2, false));
#ifdef DEBUG
	      geomcv2->writeStandardHeader(of1);
	      geomcv2->write(of1);
#endif
	    }
      }
      mod_all1.push_back(shared_ptr<ParamCurve>(mod_cvs[0]->clone()));
      mod_all2.push_back(shared_ptr<ParamCurve>(mod_cvs[1]->clone()));

      shared_ptr<BoundedSurface> mod_adj =
	replaceBdCvs(adj_bdsf[ki], mod_cvs, tol, toptol_.neighbour);
      mod_sfs.push_back(mod_adj);
    }

  shared_ptr<BoundedSurface> mod1 = replaceBdCvs(bdsf[0], mod_all1, tol, 
						 toptol_.neighbour);
  mod_sfs.push_back(mod1);
  shared_ptr<BoundedSurface> mod2 = replaceBdCvs(bdsf[1], mod_all2, tol, 
						 toptol_.neighbour);
  mod_sfs.push_back(mod2);
  modified = true;

#ifdef DEBUG
  std::ofstream of2("sliver_mod.g2");
  for (size_t kr=0; kr<mod_sfs.size(); ++kr)
    {
      mod_sfs[kr]->writeStandardHeader(of2);
      mod_sfs[kr]->write(of2);
    }
#endif

  return modified;
}

//===========================================================================
bool ftVolume::removeSliver2(shared_ptr<ftSurface> face, 
			     vector<shared_ptr<ftEdge> >& edg, 
			     int ix1, int ix2, double tol,
			     vector<shared_ptr<ParamSurface> >& mod_sfs,
			     shared_ptr<ftEdge>& not_changed)
//===========================================================================
{
  bool modified = false;
  if (edg.size() < 3 || edg.size() > 4)
    return modified;  // Not correct configuration

  shared_ptr<ParamSurface> sf[2];
  shared_ptr<BoundedSurface> bdsf[2];
  shared_ptr<SurfaceOnVolume> volsf[2];
  sf[0] = edg[ix1]->twin()->face()->asFtSurface()->surface();
  sf[1] = edg[ix2]->twin()->face()->asFtSurface()->surface();
  bdsf[0] = dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf[0]);
  bdsf[1] = dynamic_pointer_cast<BoundedSurface,ParamSurface>(sf[1]);
  if (bdsf[0].get() == NULL || bdsf[1].get() == NULL)
    return modified;
  volsf[0] = 
    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bdsf[0]->underlyingSurface());
  volsf[1] =
    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bdsf[1]->underlyingSurface());
  if (volsf[0].get() == NULL || volsf[1].get() == NULL)
    return modified;
#ifdef DEBUG
  std::ofstream of0("adj_sliver.g2");
  sf[0]->writeStandardHeader(of0);
  sf[0]->write(of0);
  sf[1]->writeStandardHeader(of0);
  sf[1]->write(of0);
#endif

  // Adjacent faces
  vector<ftSurface*> adj_faces;
  for (size_t ki=0; ki<edg.size(); ++ki)
    {
      if ((int)ki == ix1 || (int)ki == ix2)
	continue;
      adj_faces.push_back(edg[ki]->twin()->face()->asFtSurface());
    }

  // One adjacent surface is expected to be a volume boundary surface, the other not
  // Merge the sliver surface with the free adjacent surface
  int bd[2];
  int orientation;
  bool swap;
  bd[0] = volsf[0]->whichBoundary(tol, orientation, swap);
  bd[1] = volsf[1]->whichBoundary(tol, orientation, swap);

  vector<shared_ptr<ParamSurface> > sfs_merge(2);
  sfs_merge[0] = face->surface();
  int adj_ix = -1;
  if (bd[0] >= 0)
    {
      sfs_merge[1] = sf[1];
      adj_ix = ix2;
    }
  else
    {
      sfs_merge[1] = sf[0];
      adj_ix = ix1;
    }

  // Must tune the bend tolerance with respect to the actual angle between
  // the surfaces
  double angtol = toptol_.bend;
  shared_ptr<Vertex> vx[2];
  edg[adj_ix]->getVertices(vx[0], vx[1]);
  for (int ki=0; ki<1; ++ki)
    {
      vector<pair<ftSurface*, Point> > vx_faces = vx[ki]->getFaces();
      Point norm[2];
      int kr;
      size_t kj;
      for (kr=0, kj=0; kj<vx_faces.size(); ++kj)
	{
	  if (vx_faces[kj].first == edg[adj_ix]->twin()->face())
	    norm[kr++] = vx_faces[kj].first->normal(vx_faces[kj].second[0],
						    vx_faces[kj].second[1]);
	  else if (vx_faces[kj].first == face.get())
	    norm[kr++] = face->normal(vx_faces[kj].second[0],
				      vx_faces[kj].second[1]);
	  if (kr == 2)
	    break;
	}

      if (kr == 2)
	{
	  double ang = norm[0].angle(norm[1]);
	  if (ang > angtol)
	    angtol = ang + toptol_.kink;
	}
    }
  
  shared_ptr<SurfaceModel>tmp_model(new SurfaceModel(toptol_.neighbour,
						     toptol_.gap,
						     toptol_.neighbour,
						     toptol_.kink,
						     angtol,
						     sfs_merge));
  // Perform approximation
  double error;
  shared_ptr<ParamSurface> approx_surf = tmp_model->approxFaceSet(error);
  if (approx_surf.get())
    {
      mod_sfs.push_back(approx_surf);
      modified = true;
      not_changed = (adj_ix == ix1) ? edg[ix2] : edg[ix1];

      for (size_t ki=0; ki<adj_faces.size(); ++ki)
	{
	  // Adapt boundary curve towards adjacent surface with respect to
	  // the new boundary curve
	  // First identify curve
	  vector<shared_ptr<ParamCurve> > bd_cv;
	  CurveLoop loop = approx_surf->outerBoundaryLoop();
	  int nmb_cv = loop.size();
	  shared_ptr<Loop> adj_loop = adj_faces[ki]->getBoundaryLoop(0);
	  int kr;
	  for (kr=0; kr<nmb_cv; ++kr)
	    {
	      shared_ptr<ParamCurve> tmp_cv = loop[kr];
	      Point pos1 = tmp_cv->point(tmp_cv->startparam());
	      Point pos2 = tmp_cv->point(tmp_cv->endparam());

	      int ind1, ind2;
	      double par1, par2, dist1, dist2;
	      Point close1, close2;
	      adj_loop->closestPoint(pos1, ind1, par1, close1, dist1);
	      adj_loop->closestPoint(pos2, ind2, par2, close2, dist2);
	      if (dist1 < toptol_.gap && dist2 < toptol_.gap)
		{
		  bd_cv.push_back(tmp_cv);
		  break;
		}
	    }

	  if (bd_cv.size() > 0)
	    {
	      shared_ptr<ParamSurface> adj_sf = adj_faces[ki]->surface();
	      shared_ptr<BoundedSurface> bd_sf = 
		dynamic_pointer_cast<BoundedSurface,ParamSurface>(adj_sf);
	      if (bd_sf.get())
		{
		  shared_ptr<ParamSurface> adj_mod = 
		    replaceBdCvs(bd_sf, bd_cv, toptol_.gap, toptol_.neighbour);
		  mod_sfs.push_back(adj_mod);
		}
	    }
	  else 
	    modified = false;
	}
    }

#ifdef DEBUG
  std::ofstream of2("sliver_mod.g2");
  for (size_t kr=0; kr<mod_sfs.size(); ++kr)
    {
      mod_sfs[kr]->writeStandardHeader(of2);
      mod_sfs[kr]->write(of2);
    }
#endif
  return modified;
}

//===========================================================================
shared_ptr<BoundedSurface> ftVolume::replaceBdCvs(shared_ptr<BoundedSurface> surf,
						  vector<shared_ptr<ParamCurve> >& bd_cvs,
						  double tol, double tol2)
//===========================================================================
{
  CurveLoop loop = surf->outerBoundaryLoop();
  shared_ptr<ParamSurface> under_sf = surf->underlyingSurface();

  Point par_eps = SurfaceTools::getParEpsilon(*under_sf, tol);
  double epspar = std::max(tol, 0.5*(par_eps[0]+par_eps[1]));

  vector<shared_ptr<ParamCurve> > loop_cvs = loop.getCurves();
  vector<shared_ptr<CurveOnSurface> > boundary_cvs;

  // We know that the modified boundary curves are connected, but expects
  // to need to reuse some of the existing curves
  Point start_pos = bd_cvs[0]->point(bd_cvs[0]->startparam());
  Point end_pos = bd_cvs[0]->point(bd_cvs[0]->endparam());
#ifdef DEBUG
  std::ofstream of("loop_cvs.g2");
  bd_cvs[0]->writeStandardHeader(of);
  bd_cvs[0]->write(of);
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << "2" << std::endl;
  of << start_pos << std::endl;
  of << end_pos << std::endl;
#endif

  shared_ptr<CurveOnSurface> tmp_cv(new CurveOnSurface(under_sf, bd_cvs[0], false));
  tmp_cv->ensureParCrvExistence(tol);  
  int pardir;
  double parval;
  tmp_cv->isConstantCurve(tol, pardir, parval); // Are constructed to follow a boundary curve
  boundary_cvs.push_back(tmp_cv);
  bd_cvs.erase(bd_cvs.begin());
  while (bd_cvs.size() > 0)
    {
      for (size_t ki=0; ki<bd_cvs.size();)
	{
	  Point pos1 = bd_cvs[ki]->point(bd_cvs[ki]->startparam());
	  Point pos2 = bd_cvs[ki]->point(bd_cvs[ki]->endparam());
	  double d1 = start_pos.dist(pos1);
	  double d2 = start_pos.dist(pos2);
	  double d3 = end_pos.dist(pos1);
	  double d4 = end_pos.dist(pos2);
	  if (std::min(std::min(d1,d2), std::min(d3,d4)) < tol)
	    {
#ifdef DEBUG
	      bd_cvs[ki]->writeStandardHeader(of);
	      bd_cvs[ki]->write(of);
#endif
	      shared_ptr<CurveOnSurface> tmp_cv2(new CurveOnSurface(under_sf, 
								    bd_cvs[ki], false));
	      tmp_cv2->ensureParCrvExistence(tol);  
	      int pardir;
	      double parval;
	      tmp_cv2->isConstantCurve(tol, pardir, parval); // Are constructed to follow a boundary curve
	      if (std::min(d1,d2) < std::min(d3,d4))
		{
		  if (d1 < d2)
		    {
		      tmp_cv2->reverseParameterDirection();
		      start_pos = pos2;
		    }
		  else
		    start_pos = pos1;
		  boundary_cvs.insert(boundary_cvs.begin(), tmp_cv2);
		}
	      else
		{
		  if (d4 < d3)
		    {
		      tmp_cv2->reverseParameterDirection();
		      end_pos = pos1;
		    }
		  else
		    end_pos = pos2;
		  boundary_cvs.push_back(tmp_cv2);
		}
	      bd_cvs.erase(bd_cvs.begin()+ki);
#ifdef DEBUG
	      of << "400 1 0 4 255 0 0 255" << std::endl;
	      of << "2" << std::endl;
	      of << start_pos << std::endl;
	      of << end_pos << std::endl;
#endif
	    }
	  else
	    ++ki;
	}
    }

  if (boundary_cvs.size() > 1)
    {
      // Check if any of the new boundary curves are too long
      int kr;
      size_t kj;
      for (kj=0, kr=0; kr<2; kj=boundary_cvs.size()-1, ++kr)
	{
	  // Only check first and last boundary curve
	  double t1 = boundary_cvs[kj]->startparam();
	  double t2 = boundary_cvs[kj]->endparam();
	  for (size_t ki=0; ki<loop_cvs.size(); ++ki)
	    {
	      Point pos1 = loop_cvs[ki]->point(loop_cvs[ki]->startparam());
	      Point pos2 = loop_cvs[ki]->point(loop_cvs[ki]->endparam());
	      double par1, par2, dist1, dist2;
	      Point close1, close2;
	      boundary_cvs[kj]->closestPoint(pos1, t1, t2, par1, close1, dist1);
	      boundary_cvs[kj]->closestPoint(pos2, t1, t2, par2, close2, dist2);
	      if (dist1 < tol && dist2 < tol)
		{
		  // Match. Check end parameter
		  double len1 = (kr == 0) ? close1.dist(start_pos) : close1.dist(end_pos);
		  double len2 = (kr == 0) ? close2.dist(start_pos) : close2.dist(end_pos);
		  if (kr == 0 && len1 > tol && len2 > tol)
		    {
		      shared_ptr<CurveOnSurface> sub(boundary_cvs[kj]->subCurve((len1 > len2) ? par2 : par1, t2));
		      boundary_cvs[kj] = sub;
		      start_pos = (len1 > len2) ? close2 : close1;
		    }
		  else if (kr == 1 && len1 > tol && len2 > tol)
		    {
		      shared_ptr<CurveOnSurface> sub(boundary_cvs[kj]->subCurve(t1, (len1 > len2) ? par2 : par1));
		      boundary_cvs[kj] = sub;
		      end_pos = (len1 > len2) ? close2 : close1;
		    }
		  break;
		}
	    }
	}
    }
#ifdef DEBUG
  std::ofstream of3("loop_cvs2.g2");
  for (size_t ki=0; ki<boundary_cvs.size(); ++ki)
    {
      boundary_cvs[ki]->geometryCurve()->writeStandardHeader(of3);
      boundary_cvs[ki]->geometryCurve()->write(of3);
    }
#endif

#ifdef DEBUG
  std::ofstream of4("loop_cvs_orig.g2");
  for (size_t ki=0; ki<loop_cvs.size(); ++ki)
    {
      loop_cvs[ki]->geometryCurve()->writeStandardHeader(of4);
      loop_cvs[ki]->geometryCurve()->write(of4);
    }
#endif

  // Split original loop curves if the endpoints of the modified
  // boundary curves lies in the inner of a loop curve
  for (size_t ki=0; ki<loop_cvs.size(); ++ki)
    {
      double t1 = loop_cvs[ki]->startparam();
      double t2 = loop_cvs[ki]->endparam();
      Point pt1 = loop_cvs[ki]->point(t1);
      Point pt2 = loop_cvs[ki]->point(t2);
      double par1, par2, dist1, dist2;
      Point close1, close2;
      loop_cvs[ki]->closestPoint(start_pos, t1, t2, par1, close1, dist1);
      loop_cvs[ki]->closestPoint(end_pos, t1, t2, par2, close2, dist2);
      
      size_t ki2 = ki+1;
      if (dist1 < tol && par1 > t1+epspar && par1 < t2-epspar &&
	  pt1.dist(close1) > tol2 && pt2.dist(close1) > tol2)
	{
	  // Split
	  vector<shared_ptr<ParamCurve> > subcvs = loop_cvs[ki]->split(par1);
	  loop_cvs[ki] = subcvs[0];
	  loop_cvs.insert(loop_cvs.begin()+ki+1, subcvs.begin()+1, subcvs.end());
	  ki2 += (subcvs.size() - 1);
	}
      for (size_t kj=ki; kj<ki2; ++kj)
	{
	  double t3 = (kj > ki) ? loop_cvs[kj]->startparam() : t1;
	  double t4 = (kj < ki2-1) ? loop_cvs[kj]->endparam() : t2;
	  Point pt3 = (kj > ki) ? close1 : pt1;
	  Point pt4 = (kj < ki2-1) ? close1 : pt2;
	  if (dist2 < tol && par2 > t3+epspar && par2 < t4-epspar &&
	      pt2.dist(close2) > tol2 && pt4.dist(close2) > tol2)
	    {
	      // Split
	      vector<shared_ptr<ParamCurve> > subcvs = loop_cvs[kj]->split(par2);
	      loop_cvs[kj] = subcvs[0];
	      loop_cvs.insert(loop_cvs.begin()+kj+1, subcvs.begin()+1, subcvs.end());
	    }
	}
    }

#ifdef DEBUG
  std::ofstream of4_2("loop_cvs_mod.g2");
  for (size_t ki=0; ki<loop_cvs.size(); ++ki)
    {
      loop_cvs[ki]->geometryCurve()->writeStandardHeader(of4_2);
      loop_cvs[ki]->geometryCurve()->write(of4_2);
    }
#endif

  size_t prev_cvs_size = loop_cvs.size();
  while (true)
    {
      if (start_pos.dist(end_pos) < tol)
	break;
      if (loop_cvs.size() == 0)
	break;
      if (start_pos.dist(end_pos) < tol2 && loop_cvs.size() == prev_cvs_size)
	break;

      prev_cvs_size = loop_cvs.size();

      // Add original loop curves
      for (size_t ki=0; ki<loop_cvs.size(); )
	{
	  Point pos1 = loop_cvs[ki]->point(loop_cvs[ki]->startparam());
	  Point pos2 = loop_cvs[ki]->point(loop_cvs[ki]->endparam());
	  double d1 = start_pos.dist(pos1);
	  double d2 = start_pos.dist(pos2);
	  double d3 = end_pos.dist(pos1);
	  double d4 = end_pos.dist(pos2);
	  if (boundary_cvs.size() == 1 &&
	      std::min(d1,d2) < tol && std::min(d3,d4) < tol)
	    {
	      ++ki;
	      continue;   // The same curve?
	    }
	  else if (std::min(std::min(d1,d2), std::min(d3,d4)) < tol)
	    {
	      bool first_cv = (std::min(d1,d2) < std::min(d3,d4));
	      bool turn = ((first_cv && d1 < d2) || ((!first_cv) && d4 < d3));

	      // Check that the loop does not turn back on itself
	      // Test other endpont and midpoint of curve
	      Point testpos[2];
	      Point otherpos = (first_cv) ? end_pos : start_pos;
	      if (first_cv)
		testpos[0] = turn ? pos2 : pos1;
	      else
		testpos[0] = turn ? pos1 : pos2;
	      double midpar = 0.5*(loop_cvs[ki]->startparam() + 
				   loop_cvs[ki]->endparam());
	      testpos[1] = loop_cvs[ki]->point(midpar);
	      size_t ix;
	      bool found = false;
	      for (ix=0; ix<boundary_cvs.size(); ++ix)
		{
		  for (int ka=0; ka<2; ++ka)
		    {
		      double par, dist;
		      Point close;
		      boundary_cvs[ix]->closestPoint(testpos[ka], 
						     boundary_cvs[ix]->startparam(),
						     boundary_cvs[ix]->endparam(), 
						     par, close, dist);
		      if (dist < tol2)
			{
		      // double len = first_cv ? close.dist(end_pos) : close.dist(start_pos);
		      // // if (len > tol2)
		      // // 	break;
		      // // else
		      // 	if (len > tol /*tol2*/)
		      // 	{
		      // 	  // Check for a better option
		      // 	  bool found = false;
		      // 	  for (size_t kj=ki+1; kj<loop_cvs.size(); ++kj)
		      // 	    {
		      // 	      Point pos3 = loop_cvs[kj]->point(loop_cvs[kj]->startparam());
		      // 	      Point pos4 = loop_cvs[kj]->point(loop_cvs[kj]->endparam());
		      // 	      double d1_2 = start_pos.dist(pos3);
		      // 	      double d2_2 = start_pos.dist(pos4);
		      // 	      double d3_2 = end_pos.dist(pos3);
		      // 	      double d4_2 = end_pos.dist(pos4);
		      // 	      if (std::min(std::min(d1_2,d2_2), std::min(d3_2,d4_2)) < tol)
		      // 		{
		      // 		  found = true;

		      // 		  // Check that the loop does not turn back on itself
		      // 		  bool first_cv2 = (std::min(d1_2,d2_2) < std::min(d3_2,d4_2));
		      // 		  bool turn2 = ((first_cv2 && d1_2 < d2_2) || 
		      // 				((!first_cv2) && d4_2 < d3_2));
		      // 		  Point pos5;
		      // 		  if (first_cv2)
		      // 		    pos5 = turn2 ? pos4 : pos3;
		      // 		  else
		      // 		    pos5 = turn2 ? pos3 : pos4;
		      // 		  size_t ix2;
		      // 		  for (ix2=0; ix2<boundary_cvs.size(); ++ix2)
		      // 		    {
		      // 		      double par2, dist2;
		      // 		      Point close2;
		      // 		      boundary_cvs[ix2]->closestPoint(pos5, 
		      // 						      boundary_cvs[ix2]->startparam(),
		      // 						      boundary_cvs[ix2]->endparam(), 
		      // 						      par2, close2, dist2);
		      // 		      if (dist2 < tol2)
		      // 			{
		      // 			  double len2 = first_cv2 ? close2.dist(end_pos) : 
		      // 			    close2.dist(start_pos);
		      // 			  if (len2 > tol)
		      // 			    found = false;
		      // 			}
		      // 		    }
		      // 		}
		      // 	      if (found)
		      // 		break;
			  if (!(ka == 0 && close.dist(otherpos) < tol2))
			      found = true;
			}
		    }
		  if (found)
		    break;
		}

	      if (ix < boundary_cvs.size())
		{
		  ++ki;
		  continue;
		}
		
	      shared_ptr<CurveOnSurface> sf_cv = 
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(loop_cvs[ki]);
#ifdef DEBUG
	      sf_cv->spaceCurve()->writeStandardHeader(of);
	      sf_cv->spaceCurve()->write(of);
#endif
	      shared_ptr<CurveOnSurface> tmp_cv2(sf_cv->clone());
	      if (turn)
		{
		  // Reverse assembled curves
		  for (size_t kj=0; kj<boundary_cvs.size(); ++kj)
		    boundary_cvs[kj]->reverseParameterDirection();
		  size_t ncvs = boundary_cvs.size();
		  for (size_t kj=0; kj<ncvs/2; ++kj)
		    std::swap(boundary_cvs[kj], boundary_cvs[ncvs-kj-1]);
		  first_cv = !first_cv;
		}
	      if (first_cv)
		boundary_cvs.insert(boundary_cvs.begin(), tmp_cv2);
	      else
		boundary_cvs.push_back(tmp_cv2);
	      ix = boundary_cvs.size() - 1;
	      start_pos = boundary_cvs[0]->ParamCurve::point(boundary_cvs[0]->startparam());
	      end_pos = boundary_cvs[ix]->ParamCurve::point(boundary_cvs[ix]->endparam());
	      loop_cvs.erase(loop_cvs.begin()+ki);
#ifdef DEBUG
	      of << "400 1 0 4 255 0 0 255" << std::endl;
	      of << "2" << std::endl;
	      of << start_pos << std::endl;
	      of << end_pos << std::endl;
#endif

	      // // Check if the start/end curve is too long
	      // double par2, dist2;
	      // Point close2;
	      // if (first_cv)
	      // 	boundary_cvs[ix]->closestPoint(start_pos, boundary_cvs[ix]->startparam(),
	      // 				       boundary_cvs[ix]->endparam(), par2, 
	      // 				       close2, dist2);
	      // else
	      // 	boundary_cvs[0]->closestPoint(end_pos, boundary_cvs[0]->startparam(),
	      // 				      boundary_cvs[0]->endparam(), par2, 
	      // 				      close2, dist2);
	      // double dist3 = (first_cv) ? end_pos.dist(close) : start_pos.dist(close);
	      // if (dist2 < tol && dist3 >= tol)
	      // 	{
	      // 	  // Reduce curve length
	      // 	  if (first_cv)
	      // 	    {
	      // 	      shared_ptr<CurveOnSurface> sub(boundary_cvs[ix]->subCurve(boundary_cvs[ix]->startparam(),
	      // 									par2));
	      // 	      boundary_cvs[ix] = sub;
	      // 	      end_pos = close2;
	      // 	    }
	      // 	  else
	      // 	    {
	      // 	      shared_ptr<CurveOnSurface> sub(boundary_cvs[0]->subCurve(par2,
	      // 								       boundary_cvs[0]->endparam()));
	      // 	      boundary_cvs[0] = sub;
	      // 	      start_pos = close2;
	      // 	    }
	      // 	}

	      if (start_pos.dist(end_pos) < tol)
		break;
	    }
	  else
	    ++ki;
	}
    }

  shared_ptr<BoundedSurface> mod_sf(new BoundedSurface(under_sf, boundary_cvs, tol));
#ifdef DEBUG
   std::ofstream of2("bd_sf.g2");
   mod_sf->writeStandardHeader(of2);
   mod_sf->write(of2);
#endif

  return mod_sf;
}

// //===========================================================================
// shared_ptr<ParamCurve> ftVolume::makeMissingEdgeCv(shared_ptr<Vertex> vx1,
// 						   shared_ptr<Vertex> vx2)
// //===========================================================================
// {
//   bool linear_cv = true;  // Default action is to create a linear curve
//   Point d1(0.0, 0.0, 0.0), d2(0.0, 0.0, 0.0);

//   // For each edge coming into the two vertices, check if a linear curve
//   // would create a problem. In that case define the curve tangent
//   // corresponding to this vertex
//   size_t ki, kj;
//   Point vx1_pt = vx1->getVertexPoint();
//   Point vx2_pt = vx2->getVertexPoint();
//   Point vec = vx2_pt - vx1_pt;
//   double len = vec.length();
//   vec.normalize();
//   double fac = 0.9;
//   double ang_tol = 0.1*M_PI;

//   vector<ftEdge*> edges1 = vx1->uniqueEdges();
//   vector<Point> tan1(edges1.size());
//   for (ki=0; ki<edges1.size(); ++ki)
//     {
//       double t1 = edges1[ki]->parAtVertex(vx1.get());
//       tan1[ki] = edges1[ki]->tangent(t1);
//       if (edges1[ki]->tMax() - t1 < t1 - edges1[ki]->tMin())
// 	tan1[ki] *= -1;
//     }

//   for (ki=0; ki<tan1.size(); ++ki)
//     {
//       // Check the new edge with respect to the current tangent
//       // Project all other tangents into the plane defined by this
//       // tangent and the new edge curve
//       if (tan1[ki].length() < toptol_.gap)
// 	continue;  // Not possible to define plane

//       Point normal = tan1[ki].cross(vec);
//       if (normal.length() < toptol_.gap)
// 	continue;
//       normal.normalize();

//       Point avtan(0.0, 0.0, 0.0);
//       for (kj=0; kj<tan1.size(); ++kj)
// 	{
// 	  if (kj == ki)
// 	    continue;
// 	  Point tmp = tan1[kj] - (tan1[kj]*normal)*normal;
// 	  avtan += tmp;
// 	}

// 	//if (avtan.length() < toptol_.gap)
//       if (avtan.length() < toptol_.neighbour)
// 	continue;
//       double len2 = tan1[ki].length();
//       avtan.normalize();
//       avtan *= len2;
//       double ang1 = tan1[ki].angle(avtan);
//       double ang2 = tan1[ki].angle(vec);
//       double ang3 = avtan.angle(vec);
//       if ((ang2 < ang_tol || ang3 < ang_tol) && 
// 	   std::max(ang2, ang3) > fac*ang1)
// 	{
// 	  linear_cv = false;
// 	  Point tmp = tan1[ki];
// 	  tmp.normalize();
// 	  d1 += 0.5*(tmp+avtan);
// 	}
//       // else
//       // 	d1 += vec;
//     }
      
//   vec *= -1;
//   vector<ftEdge*> edges2 = vx2->uniqueEdges();
//   vector<Point> tan2(edges2.size());
//   for (ki=0; ki<edges2.size(); ++ki)
//     {
//       double t2 = edges2[ki]->parAtVertex(vx2.get());
//       tan2[ki] = edges2[ki]->tangent(t2);
//       if (edges2[ki]->tMax() - t2 < t2 - edges2[ki]->tMin())
// 	tan2[ki] *= -1;
//     }

//   for (ki=0; ki<tan2.size(); ++ki)
//     {
//       // Check the new edge with respect to the current tangent
//       // Project all other tangents into the plane defined by this
//       // tangent and the new edge curve
//       if (tan2[ki].length() < toptol_.gap)
// 	continue;  // Not possible to define plane

//       Point normal = tan2[ki].cross(vec);
//       if (normal.length() < toptol_.gap)
// 	continue;
//       normal.normalize();

//       Point avtan(0.0, 0.0, 0.0);
//       for (kj=0; kj<tan2.size(); ++kj)
// 	{
// 	  if (kj == ki)
// 	    continue;
// 	  Point tmp = tan2[kj] - (tan2[kj]*normal)*normal;
// 	  avtan += tmp;
// 	}

//       //if (avtan.length() < toptol_.gap)
//       if (avtan.length() < toptol_.neighbour)
// 	continue;
//       double len2 = tan2[ki].length();
//       avtan.normalize();
//       avtan *= len2;
//       double ang1 = tan2[ki].angle(avtan);
//       double ang2 = tan2[ki].angle(vec);
//       double ang3 = avtan.angle(vec);
//       if ((ang2 < ang_tol || ang3 < ang_tol) && 
// 	   std::max(ang2, ang3) > fac*ang1)
// 	{
// 	  linear_cv = false;
// 	  Point tmp = tan2[ki];
// 	  tmp.normalize();
// 	  d2 += 0.5*(tmp+avtan);
// 	}
//       // else
//       // 	d2 += vec;
//     }
      
//   shared_ptr<ParamCurve> crv;
//   if (linear_cv || d1.length() < toptol_.gap || d2.length() < toptol_.gap)
//     {
//       crv = shared_ptr<ParamCurve>(new SplineCurve(vx1_pt, vx2_pt));
//     }
//   else
//     {
//       // Set length of tangents
//       double len_fac = 0.1;
//       d1.normalize();
//       d1 *= len_fac*len;
//       d2.normalize();
//       d2 *= -len_fac*len;

//       HermiteInterpolator intpol;
//       vector<Point> data(4);
//       vector<double> param(2);
//       param[0] = 0.0;
//       param[1] = len;
//       data[0] = vx1_pt;
//       data[1] = d1;
//       data[2] = vx2_pt;
//       data[3] = d2;

//       vector<double> coefs;
//       intpol.interpolate(data, param, coefs);
//       BsplineBasis basis = intpol.basis();
//       crv = shared_ptr<ParamCurve>(new SplineCurve(basis, coefs.begin(), 3));
//     }
//   return crv;
// }

//===========================================================================
void ftVolume::simplifyOuterBdShell(int degree)
//===========================================================================
{
  shared_ptr<SurfaceModel> model = shells_[0];  // Consider only outer shell
  SurfaceModelUtils::simplifySurfaceModel(model, degree);
}



//===========================================================================
void ftVolume::mergeSmoothJoints(int degree)
//===========================================================================
{
  shared_ptr<SurfaceModel> model = shells_[0];  // Consider only outer shell
  Body *bd = model->getBody();

  // Collect edges across which the continuity is g1
  double bend = model->getTolerances().bend;
  FaceConnectivityUtils<ftEdgeBase,ftSurface> connectivity;
  vector<ftEdgeBase*> vec;
  vector<shared_ptr<ftSurface> > faces = model->allFaces();
  connectivity.smoothEdges(faces, vec, bend);
  
#ifdef DEBUG
      if (vec.size() > 0)
	{
	  std::ofstream edgof0("sortedvec0.g2");
	  for (size_t ki=0; ki<vec.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> cv = vec[ki]->geomEdge()->geomCurve();
	      shared_ptr<CurveOnSurface> sf_cv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	      if (sf_cv.get())
		cv = sf_cv->spaceCurve();
	      cv->writeStandardHeader(edgof0);
	      cv->write(edgof0);
	    }
	}
#endif

  // Group connected edges
  vector<int> grp_ix;
  int ki, kj, kr, kh;
  grp_ix.push_back(0);
  for (ki=0, kr=0; ki<(int)vec.size(); kr=(int)grp_ix.size()-1, ki=grp_ix[kr])
    {
      grp_ix.push_back(ki+1);
      for (kj=ki+1; kj<(int)vec.size(); ++kj)
	{
	  ftEdge *edg2 = vec[kj]->geomEdge();
	  for (kh=grp_ix[kr]; kh<grp_ix[kr+1]; ++kh)
	    {
	      ftEdge* edg1 = vec[kh]->geomEdge();
	      if (edg1->commonVertex(edg2))
		{
		  int ka = grp_ix[grp_ix.size()-1];
		  if (kj > ka)
		    std::swap(vec[ka], vec[kj]);
		  grp_ix[grp_ix.size()-1]++;
		  break;
		}
	    }
	}
    }
  if (grp_ix[grp_ix.size()-1] < (int)vec.size())
    grp_ix.push_back((int)vec.size());

#ifdef DEBUG
      if (vec.size() > 0)
	{
	  std::ofstream edgof("sortedvec.g2");
	  for (size_t ki=0; ki<vec.size(); ++ki)
	    {
	      shared_ptr<ParamCurve> cv = vec[ki]->geomEdge()->geomCurve();
	      shared_ptr<CurveOnSurface> sf_cv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	      if (sf_cv.get())
		cv = sf_cv->spaceCurve();
	      cv->writeStandardHeader(edgof);
	      cv->write(edgof);
	    }
	}
#endif

  // For each edge group, assemble adjacent faces
  tpTolerances toptol = model->getTolerances();
  double tol = model->getApproximationTol();
  vector<set<shared_ptr<ParamSurface> > > all_sfs;
  vector<set<shared_ptr<ftSurface> > > all_faces;
  for (ki=1; ki<(int)grp_ix.size(); ++ki)
    {
      set<shared_ptr<ParamSurface> > curr_sfs;
      set<shared_ptr<ftSurface> > curr_faces;
      DirectionCone union_cone;
      for (kj=grp_ix[ki-1]; kj<grp_ix[ki]; ++kj)
	{
	  vector<ftSurface*> adj_faces = vec[kj]->geomEdge()->getAdjacentFaces();
	  for (size_t kr=0; kr<adj_faces.size(); ++kr)
	    {
	      int ix = model->getIndex(adj_faces[kr]);
	      shared_ptr<ParamSurface> adj = model->getSurface(ix);
	      DirectionCone cone = adj->normalCone();
	      if (union_cone.dimension() == 0)
		union_cone = cone;
	      else
		union_cone.addUnionWith(cone);
	      curr_faces.insert(model->getFace(ix));
	      curr_sfs.insert(adj);
	    }
	}
      if (!union_cone.greaterThanPi())
	{
	  all_sfs.push_back(curr_sfs);
	  all_faces.push_back(curr_faces);
	}
    }

  // Check if the face groups overlap
  for (ki=0; ki<(int)all_faces.size(); ++ki)
    {
      for (auto it=all_faces[ki].begin(); it!=all_faces[ki].end(); ++it)
	{
	  for (kj=ki+1; kj<(int)all_faces.size(); )
	    {
	      if (all_faces[kj].find(*it) != all_faces[kj].end())
		{
		  all_faces[ki].insert(all_faces[kj].begin(), all_faces[kj].end());
		  all_faces.erase(all_faces.begin()+kj);
		  all_sfs[ki].insert(all_sfs[kj].begin(), all_sfs[kj].end());
		  all_sfs.erase(all_sfs.begin()+kj);
		}
	      else
		++kj;
	    }
	  if (ki == (int)all_faces.size()-1)
	    break;
	}
    }
  
  for (ki=0; ki<(int)all_sfs.size(); ++ki)
    {
      vector<shared_ptr<ParamSurface> > grp_sfs(all_sfs[ki].begin(), 
						all_sfs[ki].end());
      vector<shared_ptr<ftSurface> > grp_faces(all_faces[ki].begin(), 
					       all_faces[ki].end());

      shared_ptr<SurfaceModel> grp_model(new SurfaceModel(tol, toptol.gap, 
							  toptol.neighbour,
							  toptol.kink, toptol.bend,
							  grp_sfs));

#ifdef DEBUG
      std::ofstream of("grp_sfs.g2");
      for (size_t kr=0; kr<grp_sfs.size(); ++kr)
	{
	  grp_sfs[kr]->writeStandardHeader(of);
	  grp_sfs[kr]->write(of);
	}
#endif
      // Check for sharp edges
      vector<ftEdge*> corners;
      grp_model->getCorners(corners);
      if (corners.size() == 0)
	{
	  // Continue with the merge

	  double dist;
	  shared_ptr<ParamSurface> approx_surf;
	  try {
	    approx_surf = grp_model->representAsOneSurface(dist);
	  }
	  catch (...)
	    {
	      if (approx_surf.get())
		approx_surf.reset();
	    }
	  if (approx_surf.get())
	    {
	      // Replace surface
	      shared_ptr<ftSurface> new_face(new ftSurface(approx_surf,-1));
	      (void)new_face->createInitialEdges(toptol.gap, toptol.kink);
	      for (size_t kr=0; kr<grp_faces.size(); ++kr)
		model->removeFace(grp_faces[kr]);
	      new_face->setBody(bd);
	      model->append(new_face);
	    }
	}
    }
  int stop_break = 1;
}



 //===========================================================================
bool ftVolume::checkBodyTopology()
//===========================================================================
{
  bool isOK = true;
  int idlimit = 1000;
  if (id_ < -1 || id_ > idlimit)
    {
      std::cout << "Unlikely id, body = " << this << ", id = " << id_ << std::endl;
      isOK = false;
    }

  size_t ki;
  for (ki=0; ki<shells_.size(); ++ki)
    {
      bool shellOK = shells_[ki]->checkShellTopology();
      if (!shellOK)
	isOK = false;

      int nmb = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	   Body *bd = shells_[ki]->getFace(kj)->getBody();
	   if (bd != this)
	     {
	       std::cout << "Boundary shell back pointer inconsistency, face = ";
	       std::cout << shells_[ki]->getFace(kj) << ", body1 = ";
	       std::cout << this << ", body2 = " << bd << std::endl;
	       isOK = false;
	     }
	} 
    }

  return isOK;
}

//===========================================================================
ftVolume::ParameterSurfaceOnVolume::ParameterSurfaceOnVolume(shared_ptr<ParamVolume> vol,
							     shared_ptr<ParamSurface> spacesurf)
//===========================================================================
  : SurfaceOnVolume(vol, spacesurf, 0, 0.0, -1, false)
{
}

//===========================================================================
ftVolume::ParameterSurfaceOnVolume::ParameterSurfaceOnVolume(shared_ptr<ParamVolume> vol,
							     shared_ptr<ParamSurface> spacesurf,
							     int constdir, 
							     double constpar, 
							     int boundary,
							     bool swapped, 
							     int orientation)
//===========================================================================
  : SurfaceOnVolume(vol, spacesurf, constdir, constpar, boundary, 
		    swapped, orientation)
{
}

//===========================================================================
ftVolume::ParameterSurfaceOnVolume::ParameterSurfaceOnVolume(shared_ptr<ParamVolume> vol,
							     shared_ptr<ParamSurface> psurf,
							     shared_ptr<ParamSurface> spacesurf)
//===========================================================================
  : SurfaceOnVolume(vol, psurf, spacesurf, psurf.get() ? true : false)
{
}

//===========================================================================
ftVolume::ParameterSurfaceOnVolume::ParameterSurfaceOnVolume(shared_ptr<ParamVolume> vol,
							     shared_ptr<ParamSurface> psurf,
							     shared_ptr<ParamSurface> spacesurf,
							     int constdir, 
							     double constpar, 
							     int boundary,
							     bool swapped)
//===========================================================================
  : SurfaceOnVolume(vol, spacesurf, psurf, psurf.get() ? true : false,
		    constdir, constpar, boundary, swapped)
{
}

//===========================================================================
void ftVolume::ParameterSurfaceOnVolume::point(Point& pt, double upar, 
					       double vpar) const
//===========================================================================
{
  pt = volumeParameter(upar, vpar);
}

//===========================================================================
void ftVolume::ParameterSurfaceOnVolume::point(vector<Point>& pts, double upar, 
					       double vpar, int derivs,
					       bool u_from_right, 
					       bool v_from_right,
					       double resolution) const
//===========================================================================
{
  if (psurf_.get())
    psurf_->point(pts, upar, vpar, derivs, u_from_right, v_from_right);
  else
    {
      if (derivs > 1)
	MESSAGE("ParameterSurfaceOnVolume::point(). Only 1. derivatives supported");

      // Evaluate surface
      int nmb_der = (derivs+1)*(derivs+2)/2;
      vector<Point> sf_der(nmb_der);
      spacesurf_->point(sf_der, upar, vpar, u_from_right, v_from_right);

      // Find closest point in volume
      double u1, v1, w1, dist;
      Point vol_pt;
      volume_->closestPoint(sf_der[0], u1, v1, w1, vol_pt, dist, resolution);

      pts[0] = Point(u1, v1, w1);
      if (derivs == 1)
	{
	  // Differentiate volume
	  vector<Point> vol_der(4);
	  volume_->point(vol_der, u1, v1, w1, 1);

	  // For each partial derivative kp, find the factors (r1, s1, t1) such that
	  // r1*vol_der[1] + s1*vol_der[2] + t1*vol_der[3] = sf_der[kp]
	  // Solve by Cramers rule
	  double det = vol_der[1][0]*(vol_der[2][1]*vol_der[3][2] -
				      vol_der[2][2]*vol_der[3][1]) -
	    vol_der[1][1]*(vol_der[2][0]*vol_der[3][2] -
			   vol_der[2][2]*vol_der[3][0]) +
	    vol_der[1][2]*(vol_der[2][0]*vol_der[3][1] -
			   vol_der[2][1]*vol_der[3][0]);
	  for (int kp=1; kp<3; ++kp)
	    {
	      double r1 = (sf_der[kp][0]*(vol_der[2][1]*vol_der[3][2] -
					 vol_der[2][2]*vol_der[3][1]) -
			   sf_der[kp][1]*(vol_der[2][0]*vol_der[3][2] -
					 vol_der[2][2]*vol_der[3][0]) +
			   sf_der[kp][2]*(vol_der[2][0]*vol_der[3][1] -
					 vol_der[2][1]*vol_der[3][0]))/det;
	      double s1 = (vol_der[1][0]*(sf_der[kp][1]*vol_der[3][2] -
					  sf_der[kp][2]*vol_der[3][1]) -
			   vol_der[1][1]*(sf_der[kp][0]*vol_der[3][2] -
					  sf_der[kp][2]*vol_der[3][0]) +
			   vol_der[1][2]*(sf_der[kp][0]*vol_der[3][1] -
					  sf_der[kp][1]*vol_der[3][0]))/det;
	      double t1 = (vol_der[1][0]*(vol_der[2][1]*sf_der[kp][2] -
					  vol_der[2][2]*sf_der[kp][1]) -
			   vol_der[1][1]*(vol_der[2][0]*sf_der[kp][2] -
					  vol_der[2][2]*sf_der[kp][0]) +
			   vol_der[1][2]*(vol_der[2][0]*sf_der[kp][1] -
					  vol_der[2][1]*sf_der[kp][0]))/det;
	      // Dest
	      Point tmp = r1*vol_der[1] + s1*vol_der[2] + t1*vol_der[3];
	      dist = tmp.dist(sf_der[kp]);

	      pts[kp] = Point(r1, s1, t1);
	    }

	  for (int kr=3; kr<nmb_der; ++kr)
	    pts[kr] = Point(0.0, 0.0, 0.0);
	}
    }
}

//===========================================================================
vector<shared_ptr<ParamCurve> > 
ftVolume::ParameterSurfaceOnVolume::constParamCurves(double parameter,
						     bool pardir_is_u) const
//===========================================================================
{
  vector<shared_ptr<ParamCurve> > cvs = 
    SurfaceOnVolume::constParamCurves(parameter, pardir_is_u);

  vector<shared_ptr<ParamCurve> > cvs2(cvs.size());
  for (size_t ki=0; ki<cvs.size(); ++ki)
    {
      shared_ptr<CurveOnVolume> vol_cv = 
	dynamic_pointer_cast<CurveOnVolume,ParamCurve>(cvs[ki]);
    cvs2[ki] = 
      shared_ptr<ParamCurve>(new ParameterCurveOnVolume(vol_cv->underlyingVolume(),
							vol_cv->parameterCurve(),
							vol_cv->spaceCurve()));
    }
  return cvs2;
}

//===========================================================================
ClassType ftVolume::ParameterSurfaceOnVolume::instanceType() const
//===========================================================================
{
    return classType();
}

//===========================================================================
ftVolume::ParameterCurveOnVolume::ParameterCurveOnVolume(shared_ptr<ParamVolume> vol,
							 shared_ptr<ParamCurve> spacecrv)
//===========================================================================
  : CurveOnVolume(vol, spacecrv, false)
{
}

//===========================================================================
ftVolume::ParameterCurveOnVolume::ParameterCurveOnVolume(shared_ptr<ParamVolume> vol,
							 shared_ptr<ParamCurve> pcrv,
							 shared_ptr<ParamCurve> spacecrv)
//===========================================================================
  : CurveOnVolume(vol, pcrv, spacecrv, spacecrv.get() ? false : true)
{
}

//===========================================================================
void ftVolume::ParameterCurveOnVolume::point(Point& pt, double par) const
//===========================================================================
{
  pt = volumeParameter(par);
}

//===========================================================================
void ftVolume::ParameterCurveOnVolume::point(vector<Point>& pts, 
					     double tpar,
					     int derivs, bool from_right) const
//===========================================================================
{
  if (pcurve_)
    pcurve_->point(pts, tpar, derivs, from_right);
  else
    {
      // Evaluate curve
      vector<Point> cv_der(derivs+1);
      spacecurve_->point(cv_der, tpar, derivs);

      // Find closest point in volume
      double eps = 1.0e-6;
      double u1, v1, w1, dist;
      Point vol_pt;
      volume_->closestPoint(cv_der[0], u1, v1, w1, vol_pt, dist, eps);

      pts[0] = Point(u1, v1, w1);
      if (derivs == 1)
	{
	  // Differentiate volume
	  vector<Point> vol_der(4);
	  volume_->point(vol_der, u1, v1, w1, 1);

	  // Find the factors (r1, s1, t1) such that
	  // r1*vol_der[1] + s1*vol_der[2] + t1*vol_der[3] = cv_der[1]
	  // Solve by Cramers rule
	  double det = vol_der[1][0]*(vol_der[2][1]*vol_der[3][2] -
				      vol_der[2][2]*vol_der[3][1]) -
	    vol_der[1][1]*(vol_der[2][0]*vol_der[3][2] -
			   vol_der[2][2]*vol_der[3][0]) +
	    vol_der[1][2]*(vol_der[2][0]*vol_der[3][1] -
			   vol_der[2][1]*vol_der[3][0]);
	  double r1 = (cv_der[1][0]*(vol_der[2][1]*vol_der[3][2] -
				     vol_der[2][2]*vol_der[3][1]) -
		       cv_der[1][1]*(vol_der[2][0]*vol_der[3][2] -
				     vol_der[2][2]*vol_der[3][0]) +
		       cv_der[1][2]*(vol_der[2][0]*vol_der[3][1] -
				     vol_der[2][1]*vol_der[3][0]))/det;
	  double s1 = (vol_der[1][0]*(cv_der[1][1]*vol_der[3][2] -
				      cv_der[1][2]*vol_der[3][1]) -
		       vol_der[1][1]*(cv_der[1][0]*vol_der[3][2] -
				      cv_der[1][2]*vol_der[3][0]) +
		       vol_der[1][2]*(cv_der[1][0]*vol_der[3][1] -
				      cv_der[1][1]*vol_der[3][0]))/det;
	  double t1 = (vol_der[1][0]*(vol_der[2][1]*cv_der[1][2] -
				      vol_der[2][2]*cv_der[1][1]) -
		       vol_der[1][1]*(vol_der[2][0]*cv_der[1][2] -
				      vol_der[2][2]*cv_der[1][0]) +
		       vol_der[1][2]*(vol_der[2][0]*cv_der[1][1] -
				      vol_der[2][1]*cv_der[1][0]))/det;
	  // Dest
	  Point tmp = r1*vol_der[1] + s1*vol_der[2] + t1*vol_der[3];
	  dist = tmp.dist(cv_der[1]);

	  pts[1] = Point(r1, s1, t1);

	  for (int kr=2; kr<=derivs; ++kr)
	    pts[kr] = Point(0.0, 0.0, 0.0);
	}
    }
}

//===========================================================================
ftVolume::ParameterCurveOnVolume* 
ftVolume::ParameterCurveOnVolume::subCurve(double from_par, double to_par,
					   double fuzzy) const
//===========================================================================
{
  shared_ptr<CurveOnVolume> subcurve1(CurveOnVolume::subCurve(from_par, to_par, fuzzy));
  ParameterCurveOnVolume* subcurve2 = 
    new ParameterCurveOnVolume(subcurve1->underlyingVolume(),
			       subcurve1->parameterCurve(),
			       subcurve1->spaceCurve());
  return subcurve2;
}

//===========================================================================
ClassType ftVolume::ParameterCurveOnVolume::instanceType() const
//===========================================================================
{
    return classType();
}

