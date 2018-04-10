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

#include "GoTools/trivariatemodel/CreateTrimVolume.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/compositemodel/Loop.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/GapRemoval.h"
#include <fstream>
#include <cstdlib>

#define DEBUG

using std::vector;
using std::set;
using std::make_pair;

using namespace Go;

//==========================================================================
CreateTrimVolume::CreateTrimVolume(shared_ptr<SurfaceModel> model,
				   int material)
//==========================================================================
{
  model_ = model;
  material_ = material;
}

//==========================================================================
CreateTrimVolume::~CreateTrimVolume()
//==========================================================================
{

}

//==========================================================================
shared_ptr<ftVolume> CreateTrimVolume::fetchOneTrimVol()
//==========================================================================
{
  shared_ptr<ftVolume> dummy;

  // Simplify input shell and mend gaps due to bad trimming curves
  int degree = 3;
  repairShell(degree);

#ifdef DEBUG
  std::ofstream of1("simplified_shell.g2");
  int nmb = model_->nmbEntities();
  for (int kj=0; kj<nmb; ++kj)
    {
      shared_ptr<ParamSurface> sf = model_->getSurface(kj);
      sf->writeStandardHeader(of1);
      sf->write(of1);
    }
#endif

  // Distinguish between boundary surfaces and trimming surfaces
  bigbox_ = model_->boundingBox();
  vector<pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > > side_sfs;
  identifyBoundaryFaces(side_sfs);
  if (side_sfs.size() < 6)
    return dummy;

  // Represent all boundary surfaces with non-trimmed spline surfaces
  // Create parametric spline volume
  // First define selected side surface as a shell
  vector<shared_ptr<ParamSurface> > shell_sfs;
  tpTolerances tol = model_->getTolerances();
  for (size_t ki=0; ki<side_sfs.size(); ++ki)
    {
      if (side_sfs[ki].second.get())
	shell_sfs.push_back(side_sfs[ki].second);
    }
  shared_ptr<SurfaceModel> shell(new SurfaceModel(tol.gap, tol.gap,
						  tol.neighbour, tol.kink,
						  tol.bend, shell_sfs));

  // Check that no faces intersect with this solid or is outside of it

  // Create intermediate ftVolume and extract parametric volume from this
  shared_ptr<ftVolume> ftvol(new ftVolume(shell));

  if (ftvol->isRegularized())
    ftvol->untrimRegular(degree);

  shared_ptr<ParamVolume> vol = ftvol->getVolume();
#ifdef DEBUG
  std::ofstream of7("under_vol.g2");
  vol->writeStandardHeader(of7);
  vol->write(of7);
#endif

  // Reparameterize volume to match geometry
  double u_len, v_len, w_len;
  vol->estimateVolSize(u_len, v_len, w_len, 3, 3, 3);
  shared_ptr<SplineVolume> vol2 = dynamic_pointer_cast<SplineVolume>(vol);
  vol2->setParameterDomain(0.0, u_len, 0.0, v_len, 0.0, w_len);

  // // Insert knots at iso-parametric sharp edges between trimming faces
  refineInSharpEdges(vol);

  // Trim volume with the remaining faces
  shared_ptr<ftVolume> trimvol = createTrimVolume(vol, side_sfs);

  return trimvol;
}

//==========================================================================
void 
CreateTrimVolume::identifyBoundaryFaces(vector<pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs)
//==========================================================================
{
  // Initially all faces can be boundary faces
  vector<shared_ptr<ftSurface> > bd_faces;
  vector<shared_ptr<ftSurface> > trim_faces;
  bd_faces = model_->allFaces();

  // Remove faces connected to inner trim curves from the pool of boundary
  // faces
  identifyInnerTrim(bd_faces, trim_faces);

#ifdef DEBUG
  std::ofstream of1("bd_faces1.g2");
  for (size_t ki=0; ki<bd_faces.size(); ++ki)
    {
      bd_faces[ki]->surface()->writeStandardHeader(of1);
      bd_faces[ki]->surface()->write(of1);
    }
#endif

  // Divide the remaining faces in compact sets and select the "largest"
  // one as source for boundary faces
  extractMaxSet(bd_faces, trim_faces);

#ifdef DEBUG
  std::ofstream of2("bd_faces2.g2");
  for (size_t ki=0; ki<bd_faces.size(); ++ki)
    {
      bd_faces[ki]->surface()->writeStandardHeader(of2);
      bd_faces[ki]->surface()->write(of2);
    }
#endif

  // Identify faces belonging to the same underlying surface
 tpTolerances tol = model_->getTolerances();
 SurfaceModelUtils::sameUnderlyingSurf(bd_faces, tol.gap, tol.kink,
				       face_grp_, under_sf_);

#ifdef DEBUG
  std::ofstream of3("same_faces.g2");
  for (size_t ki=0; ki<face_grp_.size(); ++ki)
    {
      for (size_t kj=0; kj<face_grp_[ki].size(); ++kj)
	{
	  face_grp_[ki][kj]->surface()->writeStandardHeader(of3);
	  face_grp_[ki][kj]->surface()->write(of3);
	}
    }
#endif
  // Extend face groups with remaining candidate boundary faces
  for (size_t ki=0; ki<bd_faces.size(); ++ki)
    {
      // Check if the current face belongs to a group
      size_t kj, kr;
      for (kj=0; kj<face_grp_.size(); ++kj)
	{
	  for (kr=0; kr<face_grp_[kj].size(); ++kr)
	    {
	      if (face_grp_[kj][kr].get() == bd_faces[ki].get())
		break;
	    }
	  if (kr < face_grp_[kj].size())
	    break;
	}
      if (face_grp_.size() == 0 || kj == face_grp_.size())
	{
	  // Make face group with a single face
	  vector<shared_ptr<ftSurface> > single_face;
	  single_face.push_back(bd_faces[ki]);
	  face_grp_.push_back(single_face);
	  shared_ptr<ParamSurface> surf = bd_faces[ki]->surface();
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
	  if (bd_sf.get())
	    under_sf_.push_back(bd_sf->underlyingSurface());
	  else
	    under_sf_.push_back(surf);
	}
    }

#ifdef DEBUG
  std::ofstream of4("same_faces2.g2");
  for (size_t ki=0; ki<face_grp_.size(); ++ki)
    {
      for (size_t kj=0; kj<face_grp_[ki].size(); ++kj)
	{
	  face_grp_[ki][kj]->surface()->writeStandardHeader(of4);
	  face_grp_[ki][kj]->surface()->write(of4);
	}
    }
#endif
  // Compute face group information
  bbox_.resize(face_grp_.size());
  cone_.resize(face_grp_.size());
  sfsize_.resize(face_grp_.size());
  sf_type_.assign(face_grp_.size(), UNKNOWN);
  sf_pt_.resize(face_grp_.size());
  sf_axis_.resize(face_grp_.size());
  sf_centre_.resize(face_grp_.size());

  computeGroupInfo(tol.gap);
  
  // Suggest volume side surfaces
  findSideSfs(tol.gap, tol.kink, side_sfs);

  // Extend side surfaces to make sure that they intersect
  extendSurfaces(side_sfs);
#ifdef DEBUG
  std::ofstream of5("side_surfaces2.g2");
  for (size_t ki=0; ki<side_sfs.size(); ++ki)
    {
      side_sfs[ki].second->writeStandardHeader(of5);
      side_sfs[ki].second->write(of5);
    }
#endif

  // Perform intersections to limit the side surfaces to create a Brep solid
  // with 6 boundary faces
  try {
    trimSideSurfaces(side_sfs);
  }
  catch (...)
    {
      side_sfs.erase(side_sfs.begin(), side_sfs.end());
      return;
    }
#ifdef DEBUG
  std::ofstream of6("side_surfaces3.g2");
  for (size_t ki=0; ki<side_sfs.size(); ++ki)
    {
      side_sfs[ki].second->writeStandardHeader(of6);
      side_sfs[ki].second->write(of6);
    }
#endif

  int stop_break = 1;
}

//==========================================================================
void 
CreateTrimVolume::computeGroupInfo(double tol)
//==========================================================================
{
#ifdef DEBUG
  std::ofstream of("sf_pt.g2");
  of << "400 1 0 4 255 0 0 255" << std::endl;
  of << face_grp_.size() << std::endl;
#endif

  for (size_t ki=0; ki<face_grp_.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf1 = face_grp_[ki][0]->surface();
      bbox_[ki] = face_grp_[ki][0]->boundingBox();
      cone_[ki] = surf1->normalCone();
      double u_size, v_size;
      surf1->estimateSfSize(u_size, v_size);
      sfsize_[ki] = u_size*v_size;
      
      Point centre, axis, vec, norm;
      double ang;
      bool rotational = under_sf_[ki]->isAxisRotational(centre, axis, 
							vec, ang);
      bool planar = under_sf_[ki]->isPlanar(norm, tol);
      if (rotational)
	{
	  sf_type_[ki] = ROTATIONAL;
	  sf_axis_[ki] = axis;
	  sf_centre_[ki] = centre;
	}
      else if (planar)
	sf_type_[ki] = PLANAR;
      else
	sf_type_[ki] = FREEFORM;
      
      for (size_t kj=1; kj<face_grp_[ki].size(); ++kj)
	{
	  shared_ptr<ParamSurface> surf2 = face_grp_[ki][kj]->surface();
	  BoundingBox bb = face_grp_[ki][kj]->boundingBox();
	  bbox_[ki].addUnionWith(bb);
	  DirectionCone cc = surf2->normalCone();
	  cone_[ki].addUnionWith(cc);
	  double u_size2, v_size2;
	  surf2->estimateSfSize(u_size2, v_size2);
	  sfsize_[ki] += (u_size2*v_size2);
	}

      // Compute characteristic point on surface
      Point bpt = 0.5*(bbox_[ki].low()+bbox_[ki].high());

      vector<pair<Point, Point> > sfpts =
	BoundedUtils::intersectWithLine(under_sf_[ki], bpt,
					cone_[ki].centre(), tol);
      int ix = -1;
      double dist = HUGE;
      for (size_t kj=0; kj<sfpts.size(); ++kj)
	{
	  double dist2 = bpt.dist(sfpts[kj].second);
	  if (dist2 < dist)
	    {
	      dist = dist2;
	      ix = (int)kj;
	    }
	}
      double *guess = NULL;
      Point pos;
      if (ix < 0)
	{
	  pos = bpt;
	}
      else
	{
	  pos = sfpts[ix].second;
	  guess = sfpts[ix].first.begin();
	}
      Point close_pt;
      double close_dist = HUGE;
      for (size_t kj=0; kj<face_grp_[ki].size(); ++kj)
	{
	  double upar, vpar, dist2;
	  Point close;
	  shared_ptr<ParamSurface> surf2 = face_grp_[ki][kj]->surface();
	  surf2->closestPoint(pos, upar, vpar, close,
			      dist2, tol, NULL, guess);
	  if (dist2 < close_dist)
	    {
	      close_dist = dist2;
	      close_pt = close;
	    }
	}
#ifdef DEBUG
      of << close_pt << std::endl;
#endif
      sf_pt_[ki] = close_pt;
    }
}

//==========================================================================
void 
CreateTrimVolume::extendSurfaces(vector<pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs)
//==========================================================================
{
  // Model boundary box
  BoundingBox tot = model_->boundingBox();
  Point high1 = tot.high();
  Point low1 = tot.low();
  double l1 = high1.dist(low1);
  for (size_t ki=0; ki<side_sfs.size(); ++ki)
    {
      BoundingBox curr = side_sfs[ki].second->boundingBox();
      Point high2 = curr.high();
      Point low2 = curr.low();
      double l2 = high2.dist(low2);
      double len = std::max(2.0*l1/l2, l1/4.0);
      ElementarySurface *elem = side_sfs[ki].second->elementarySurface();
      if (elem)
	{
	  shared_ptr<ElementarySurface> elem2(elem->clone());
	  elem2->enlarge(len, len, len, len);
	  side_sfs[ki].second = elem2;
	}
    }
}

//==========================================================================
void 
CreateTrimVolume::trimSideSurfaces(vector<pair<shared_ptr<ftSurface>, 
				   shared_ptr<ParamSurface> > >& side_sfs)
//==========================================================================
{
  // Fetch tolerances
  tpTolerances tol = model_->getTolerances();

  // Represent as spline surfaces
  vector<shared_ptr<ParamSurface> > spline_sfs(side_sfs.size());
  for (size_t ki=0; ki<side_sfs.size(); ++ki)
    {
      shared_ptr<SplineSurface> sf = 
	dynamic_pointer_cast<SplineSurface, ParamSurface>(side_sfs[ki].second);
      if (sf.get())
	spline_sfs[ki] = sf;
      else
	{
	  ElementarySurface *elem = side_sfs[ki].second->elementarySurface();
	  sf = shared_ptr<SplineSurface>(elem->createSplineSurface());
	  spline_sfs[ki] = sf;
	}
    }

  // Trim the surfaces one by one with respect to adjacent surfaces
  vector<shared_ptr<BoundedSurface> >  sfs(side_sfs.size());
  for (size_t ki=0; ki<side_sfs.size(); ++ki)
    {
      // Surface model with the current surface
      shared_ptr<ParamSurface>  curr_sf;
      curr_sf = shared_ptr<ParamSurface>(spline_sfs[ki]->clone());

      // The remaining surfaces
      vector<shared_ptr<ParamSurface> > other_sfs;
      for (size_t kj=0; kj<side_sfs.size(); ++kj)
	{
	  if (ki/2 == kj/2)
	    continue;
	  other_sfs.push_back(spline_sfs[kj]);
	}

      // Perform trimming
      vector<shared_ptr<BoundedSurface> >  trim_sfs = 
	BoundedUtils::trimSurfWithSurfs(curr_sf, other_sfs, tol.gap);

      if (trim_sfs.size() == 0)
	continue;
#ifdef DEBUG
      std::ofstream of("curr_split_sf.g2");
      for (size_t kj=0; kj<trim_sfs.size(); ++kj)
	{
	  trim_sfs[kj]->writeStandardHeader(of);
	  trim_sfs[kj]->write(of);
	}
#endif

      // Select output surface
      int ix = 0;
      int nmb_bd1 = trim_sfs[ix]->numberOfLoops();

      // Compute distance between an arbitrary point in the surface
      // and the trimming shell
      double u1, v1, d1;
      double par1[2];
      int trim_ix1;
      Point pt1 = trim_sfs[ix]->getInternalPoint(u1, v1);
      Point clo_pt1;
      model_->closestPoint(pt1, clo_pt1, trim_ix1, par1, d1);
      for (size_t kj=1; kj<trim_sfs.size(); ++kj)
	{
	  double u2, v2, d2;
	  double par2[2];
	  int trim_ix2;
	  Point pt2 = trim_sfs[kj]->getInternalPoint(u2, v2);
	  Point clo_pt2;
	  model_->closestPoint(pt1, clo_pt2, trim_ix2, par2, d2);
	  int nmb_bd2 = trim_sfs[kj]->numberOfLoops();
 	  if (nmb_bd2 < nmb_bd1 || (nmb_bd2 == nmb_bd1 && d2 < d1))
	    {
	      ix = (int)kj;
	      nmb_bd1 = nmb_bd2;
	      d1 = d2;
	    }
	}
      sfs[ki] = trim_sfs[ix];

      int stop_break;
    }

  for (size_t ki=0; ki<side_sfs.size(); ++ki)
    {
      if (sfs[ki].get())
	side_sfs[ki].second = sfs[ki];
    }
}

//==========================================================================
void 
CreateTrimVolume::findSideSfs(double tol, double angtol,
			      vector<pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs)
//==========================================================================
{
  // Sort face groups according to size
  vector<int> prio(face_grp_.size());
  for (size_t ki=0; ki<prio.size(); ++ki)
    prio[ki] = (int)ki;

  for (size_t ki=0; ki<prio.size(); ++ki)
    for (size_t kj=ki+1; kj<prio.size(); ++kj)
      {
	if (sfsize_[prio[ki]] < sfsize_[prio[kj]])
	  std::swap(prio[ki], prio[kj]);
      }

  // Set threshold for importance
  // Compute mean size of the largest face groups
  int nmb0 = std::min(6, std::max(20, (int)prio.size()/5));
  nmb0 = std::min(nmb0, (int)prio.size());
  double mean = 0.0;
  for (int kr=0; kr<nmb0; ++kr)
    mean += sfsize_[prio[kr]];
  mean /= nmb0;

  // Identify candidate surfaces
  // @@@ VSK. This decision needs tuning
  double mean_frac = 0.1*mean;
  int nmb = 0;
  for (size_t ki=0; ki<prio.size(); ++ki)
    {
      nmb++;
      if (sfsize_[prio[ki]] < mean_frac)
	break;
    }
  nmb = std::min((int)prio.size(), std::max(nmb, 6));

#ifdef DEBUG
  std::ofstream ofp("prio_faces.g2");
  for (int kr=0; kr<nmb; ++kr)
    {
      int ixs = prio[kr];
      for (size_t kh=0; kh<face_grp_[ixs].size(); ++kh)
	{
	  face_grp_[ixs][kh]->surface()->writeStandardHeader(ofp);
	  face_grp_[ixs][kh]->surface()->write(ofp);
	}
    }
  ofp << "400 1 0 4 255 0 0 255" << std::endl;
  ofp << nmb << std::endl;
  for (int kr=0; kr<nmb; ++kr)
    ofp << sf_pt_[prio[kr]] << std::endl;
#endif
  // Sort candidate surfaces with respect to volume boundary surfaces and remove
  // those that are not suitable for a boundary surface
  // Categorize boundary surfaces
  vector<int> bd_type(6, 0);   // 0=unknown, 1=planar, 2=rotational, 3=other
  vector<Point> bd_vec(6);
  vector<vector<int> > grp_ix(6);  // Index of candidate boundary surface
  vector<Point> dir(6); // Suggested outwards main direction of boundary surface
  Point coord[6]; // Initial outwards directions
  Point coord_pos[6];

  // Make initial outwards directions based on analysis of the prioitized
  // surfaces
  analyzePrio(&prio[0], nmb, coord, coord_pos);

  double lim_ang = M_PI/4.0;  // 45 degrees
  double ang_threshold = M_PI/12.0;
  double dist_threshold = 0.1*(bigbox_.low().dist(bigbox_.high()));
  for (int idx=0; idx<6; ++idx)
    {
      // Find best fit surface
      int idx2 = -1;
      double min_ang = HUGE;
      double min_dist = HUGE;
      for (int kr=0; kr<nmb; ++kr)
	{
	  Point vec = cone_[prio[kr]].centre();
	  Point pt = sf_pt_[prio[kr]];
	  double ang = vec.angle(coord[idx]);
	  double dist = (coord_pos[idx].dimension() == 0) ? 0.0 : 
	    pt.dist(coord_pos[idx]);
	  if (ang < min_ang && dist < min_dist)
	    {
	      idx2 = kr;
	      min_ang = ang;
	      min_dist = dist;
	    }
	  else if (ang < min_ang && dist < dist_threshold)
	    {
	      idx2 = kr;
	      min_ang = ang;
	    }
	  else if (ang < ang_threshold && dist < min_dist)
	    {
	      idx2 = kr;
	      min_dist = dist;
	    }
	}
#ifdef DEBUG
      std::ofstream ofs("curr_faces1.g2");
      int ixs = prio[idx2];
      for (size_t kh=0; kh<face_grp_[ixs].size(); ++kh)
	{
	  face_grp_[ixs][kh]->surface()->writeStandardHeader(ofs);
	  face_grp_[ixs][kh]->surface()->write(ofs);
	}
#endif
      Point centre, axis, vec, normal;
      double angle;
      bool rot = under_sf_[prio[idx2]]->isAxisRotational(centre, axis, vec, angle);
      bool planar = under_sf_[prio[idx2]]->isPlanar(normal, tol);
      Point cone_centre = cone_[prio[idx2]].centre();

      dir[idx] = cone_centre;
      grp_ix[idx].push_back(prio[idx2]);
      if (bd_type[idx] == 0)
	{
	  bd_type[idx] = planar ? 1 : ((rot) ? 2 : 3); 
	  bd_vec[idx] = planar ? normal : ((rot) ? axis : cone_centre);
	}

      // Surface is classified
      prio.erase(prio.begin()+idx2);
      nmb--;

      // Examine remaining prioritized surfaces to find if the belong to the
      // same direction
      for (int kh=0; kh<nmb;)
	{
#ifdef DEBUG
	  std::ofstream ofs2("curr_faces12.g2");
	  int ixs2 = prio[kh];
	  for (size_t ka=0; ka<face_grp_[ixs2].size(); ++ka)
	    {
	      face_grp_[ixs2][ka]->surface()->writeStandardHeader(ofs2);
	      face_grp_[ixs2][ka]->surface()->write(ofs2);
	    }
#endif
	  Point centre2, axis2, vec2, normal2;
	  double angle2;
	  bool rot2 = under_sf_[prio[kh]]->isAxisRotational(centre2, axis2, 
							    vec2, angle2);
	  bool planar2 = under_sf_[prio[kh]]->isPlanar(normal2, tol);
	  Point cone_centre2 = cone_[prio[kh]].centre();

	  double ang = dir[idx].angle(cone_centre2);
	  double scpr = dir[idx]*cone_centre2;
	  double axis_ang = bd_vec[idx].angle(cone_centre2);
	  if (rot2)
	    axis_ang = bd_vec[idx].angle(axis2);
	  else if (planar2)
	    axis_ang = bd_vec[idx].angle(normal2);
		
	  bool box_overlap = 
	    bbox_[grp_ix[idx][0]].overlaps(bbox_[prio[kh]]);
	  int found_match = 0;
	  if (ang < lim_ang && 
	      (axis_ang < lim_ang || axis_ang > M_PI-lim_ang) && scpr > 0.0)
	    {
	      int ix1 = grp_ix[idx][0];
	      int ix2 = prio[kh];
	      found_match = checkCandPair(dir[idx], under_sf_[ix1], 
					  bd_type[idx], bbox_[ix1],
					  under_sf_[ix2], 
					  planar2 ? 1 : ((rot2) ? 2 : 3),
					  bbox_[ix2], tol);

	      if (found_match == 1 && box_overlap)
		{
		  // No clear distinction between the surfaces
		  // The previously selected is judged superios
		  grp_ix[idx].push_back(prio[kh]);
		}
	      else if (found_match == 2 && box_overlap)
		{
		  // No clear distinction between the surfaces
		  // The new surface is judged superios
		  grp_ix[idx].insert(grp_ix[idx].begin(), prio[kh]);
		  bd_type[idx] = planar2 ? 1 : ((rot2) ? 2 : 3);
		  bd_vec[idx] = planar2 ? normal2 : 
		    ((rot2) ? axis2 : cone_centre2);
		  dir[idx] = cone_centre2;
		}
	      else if (found_match == 3)
		{
		  found_match = 0;
		  // The current surface is not a candidate boundary surface
		}
	      else if (found_match == 4)
		{
		  // The current surface is a better candidate for a boundary
		  // surface than the previously selected one. Redo selection
		  bd_type[idx] = planar2 ? 1 : ((rot2) ? 2 : 3);
		  bd_vec[idx] = planar2 ? normal2 : 
		    ((rot2) ? axis2 : cone_centre2);
		  dir[idx] = cone_centre2;
		  grp_ix[idx][0] = prio[kh];
	      
		  // Should probably test agains other identified surfaces as well
		}
	      else
		found_match = 0;

	      if (found_match > 0)
		{
		  // Surface is classified
		  prio.erase(prio.begin()+kh);
		  nmb--;
		}
	      else
		++kh;

	      int stop_break0 = 1;
	    }
	  else 
	    ++kh;
	}
#ifdef DEBUG
      std::ofstream of0("cand_faces_ix.g2");
      for (size_t kj=0; kj<grp_ix[idx].size(); ++kj)
	{
	  int ix = grp_ix[idx][kj];
	  for (size_t kh=0; kh<face_grp_[ix].size(); ++kh)
	    {
	      face_grp_[ix][kh]->surface()->writeStandardHeader(of0);
	      face_grp_[ix][kh]->surface()->write(of0);
	    }
	}
#endif
	  // // Update boundary side directions
	  // // First current
	  // DirectionCone cone2 = cone_[grp_ix[idx][0]];
	  // for (size_t kj=1; kj<grp_ix[idx].size(); ++kj)
	  //   cone2.addUnionWith(cone_[grp_ix[idx][kj]]);
	  // coord[idx] = dir[idx] = cone2.centre();

	  // // Opposite direction
	  // int other_ix = (idx%2 == 0) ? idx+1 : idx-1;
	  // if (dir[other_ix].dimension() == 0)
	  //   coord[other_ix] = -coord[idx];

	  // // Remaining directions
      
    }

#ifdef DEBUG
  std::ofstream of("candidate_faces.g2");
  for (size_t ki=0; ki<grp_ix.size(); ++ki)
    {
      for (size_t kj=0; kj<grp_ix[ki].size(); ++kj)
	{
	  int ix = prio[grp_ix[ki][kj]];
	  for (size_t kh=0; kh<face_grp_[ix].size(); ++kh)
	    {
	      face_grp_[ix][kh]->surface()->writeStandardHeader(of);
	      face_grp_[ix][kh]->surface()->write(of);
	    }
	}
    }
#endif

  // Decide on one unique boundary side surface
  side_sfs.resize(6);
  shared_ptr<ftSurface> dummy_face;
  shared_ptr<ParamSurface> dummy_sf;
  for (int kh=0; kh<6; ++kh)
    {
      if (grp_ix[kh].size() == 1)
	{
	  // Already one surface
	  if (face_grp_[grp_ix[kh][0]].size() == 1)
	    side_sfs[kh] = std::make_pair(face_grp_[grp_ix[kh][0]][0],
					  under_sf_[grp_ix[kh][0]]);
	  else
	    side_sfs[kh] = std::make_pair(dummy_face,
					  under_sf_[grp_ix[kh][0]]);
	}
      else if (grp_ix[kh].size() > 1)
	{
	  // Select/construct side surface
	  oneSideSf(bd_type[kh], grp_ix[kh], bd_vec[kh], coord[kh], 
		    tol, angtol, side_sfs[kh]);
	  
	}
      else
	side_sfs[kh] = std::make_pair(dummy_face, dummy_sf);
    }	
	  
#ifdef DEBUG
  std::ofstream of2("side_surfaces.g2");
  for (size_t ki=0; ki<side_sfs.size(); ++ki)
    {
      side_sfs[ki].second->writeStandardHeader(of2);
      side_sfs[ki].second->write(of2);
    }
#endif
  int stop_break = 1;

}

//==========================================================================
void
CreateTrimVolume::oneSideSf(int bd_type, vector<int>& face_grp_ix, 
			    Point bd_vec, Point dir, double tol, double angtol,
			    pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> >& side_sf)
//==========================================================================
{
  // Check face adjacency and if there is a smooth transition
  shared_ptr<ftSurface> dummy_face;
  bool smooth = true;
  for (size_t ki=0; ki<face_grp_ix.size(); ++ki)
    {
      for (size_t kj=0; kj<face_grp_[face_grp_ix[ki]].size(); ++kj)
	{
	  shared_ptr<ftSurface> face1 = face_grp_[face_grp_ix[ki]][kj];

	  // Count adjacency relationships to other faces
	  int nmb = 0;
	  for (size_t kr=0; kr<face_grp_ix.size(); ++kr)
	    {
	      if (kr == ki)
		continue;  // Same group
	      for (size_t kh=0; kh<face_grp_[face_grp_ix[kr]].size(); ++kh)
		{
		  shared_ptr<ftSurface> face2 = face_grp_[face_grp_ix[kr]][kh];

		  if (face1->isAdjacent(face2.get(), smooth))
		    {
		      nmb++;
		      if (smooth == false)
			break;
		    }
		}
	      if (smooth == false)
		break;
	    }
	  if (nmb == 0)
	    smooth = false;
	  if (smooth == false)
	    break;
	}
      if (smooth == false)
	break;
    }

  if (smooth)
    {
      // Check if the boundary faces can be merged into one surface
      smooth = false;  // Not implemented
    }

  if (!smooth)
    {
      // Extend the surfaces, select the one with the most distant point
      // with respect to the side direction and check for intersections
      // Find extension distance. First compute combined bounding box
      BoundingBox tot = bbox_[face_grp_ix[0]];
      for (size_t ki=1; ki<face_grp_ix.size(); ++ki)
	tot.addUnionWith(bbox_[face_grp_ix[ki]]);
      Point high1 = tot.high();
      Point low1 = tot.low();

      vector<shared_ptr<ParamSurface> > ext_sfs(face_grp_ix.size());
      for (size_t ki=0;  ki<face_grp_ix.size(); ++ki)
	{
	  Point high2 = bbox_[face_grp_ix[ki]].high();
	  Point low2 = bbox_[face_grp_ix[ki]].low();
	  double l1 = high1[0] - low1[0];
	  double l2 = high2[0] - low2[0];
	  for (int ka=1; ka<high1.dimension(); ++ka)
	    {
	      if ((high2[ka]-low2[ka])/(high1[ka]-low1[ka]) < l2/l1)
		{
		  l1 = high1[ka] - low1[ka];
		  l2 = high2[ka] - low2[ka];
		}
	    }
	  double len = 1.2*(l1 -l2);
	  ElementarySurface *elem = under_sf_[face_grp_ix[ki]]->elementarySurface();
	  if (elem)
	    {
	      shared_ptr<ElementarySurface> elem2(elem->clone());
	      elem2->enlarge(len, len, len, len);
	      ext_sfs[ki] = elem2;
	    }
	  else
	    ext_sfs[ki] = under_sf_[face_grp_ix[ki]];  // To be continued
	  // by implementing enlarge for spline surfaces. Other surface types
	  // are probably not required
	}

      // Check if the other surfaces intersect the assumed most relevant one
      bool found = false;
      for (size_t ki=1; ki<ext_sfs.size(); ++ki)
	{
	  vector<shared_ptr<CurveOnSurface> > int_seg1;
	  vector<shared_ptr<CurveOnSurface> > int_seg2;
	  BoundedUtils::getIntersectionCurve(ext_sfs[0], ext_sfs[ki],
					     int_seg1, int_seg2, tol);
	  if (int_seg1.size() > 0)
	    {
	      found = true;
	      break;
	    }
	}
      
      // If no intersection is found use the primary candiatate
      if (!found)
	{
	  side_sf = std::make_pair(dummy_face, ext_sfs[0]);
	}
      else
	{
	  // Otherwise, construct a new surface outside all candidates
	  // Case distinction
	  if (bd_type == 1)
	    {
	      // A planar side surface is expected 
	      MESSAGE("Planar side surface is expected");
	    }
	  else if (bd_type == 2)
	    {
	      // A rotational side surface is expected
	      double ext_fac = 1.0;
	      double rad = 0.0;
	      Point loc, dir2;
	      for (size_t ki=0; ki<face_grp_ix.size(); ++ki)
		{
		  shared_ptr<ParamSurface> surf = under_sf_[face_grp_ix[ki]];
		  Point centre, axis, other_vec;
		  double angle;
		  bool rot = surf->isAxisRotational(centre, axis, other_vec, angle);
		  if (!rot)
		    ext_fac = 1.2;
		  else 
		    {
		      double ang = bd_vec.angle(axis);
		      if (ang > angtol && ang < M_PI-angtol)
			ext_fac = 1.2;

		      RectDomain dom = surf->containingDomain();
		      bool swap = false;
		      ElementarySurface *elem = surf->elementarySurface();
		      if (elem)
			{
			  if (loc.dimension() == 0)
			    {
			      loc = elem->location();
			      dir2 = elem->direction2();
			    }
			  swap = elem->isSwapped();
			}
		      Point par1, par2;
		      if (swap)
			{
			  par1 = Point(dom.umin(), 0.5*(dom.vmin()+dom.vmax()));
			  par2 = Point(dom.umax(), 0.5*(dom.vmin()+dom.vmax()));
			}
		      else
			{
			  par1 = Point(0.5*(dom.umin()+dom.umax()), dom.vmin());
			  par2 = Point(0.5*(dom.umin()+dom.umax()), dom.vmax());
			}
		      double rad1 = elem->radius(par1[0], par1[1]);
		      double rad2 = elem->radius(par2[0], par2[1]);
		      rad = std::max(rad, std::max(rad1, rad2));
		    }
		}
	      rad *= ext_fac;
	      
	      // Create cylinder
	      shared_ptr<Cylinder> cyl(new Cylinder(rad, loc, bd_vec, dir2));
	      RectDomain dom = ext_sfs[0]->containingDomain();
	      cyl->setParameterBounds(dom.umin(), dom.vmin(), 
				      dom.umax(), dom.vmax());  // What if swapped
	      side_sf = std::make_pair(dummy_face, cyl);
	      
	    }
	  else
	    {
	      // Possibly free form side surface
	      std::cout << "Free form side surface" << std::endl;
	    }
	}
    }
}

//==========================================================================
int
CreateTrimVolume::checkCandPair(Point vec,
				shared_ptr<ParamSurface> sf1, int bd_type1, 
				BoundingBox& box1, shared_ptr<ParamSurface> sf2,
				int bd_type2, BoundingBox& box2, 
				double tol)
//==========================================================================
{
#ifdef DEBUG
  std::ofstream of("curr_faces.g2");
  sf1->writeStandardHeader(of);
  sf1->write(of);
  sf2->writeStandardHeader(of);
  sf2->write(of);
#endif

  ElementarySurface *elem1 = sf1->elementarySurface();
  ElementarySurface *elem2 = sf2->elementarySurface();

  // Identify points on both surfaces close to each other
  Point bpt1 = 0.5*(box1.high() + box1.low());
  Point bpt2 = 0.5*(box2.high() + box2.low());
  Point pt = 0.5*(bpt1 + bpt2);
  
  Point sf_pt1, sf_pt2;
  double u1, u2, v1, v2, d1, d2;
  sf1->closestPoint(pt, u1, v1, sf_pt1, d1, tol);
  sf2->closestPoint(pt, u2, v2, sf_pt2, d2, tol);

  double frac = 0.2;

  if (elem1 && elem2)
    {
      double rad1 = elem1->radius(u1, v1);
      double rad2 = elem2->radius(u2, v2);
      if (elem1->instanceType() == Class_Plane)
	rad1 = 0.0;  // The location lies in the plane
      if (elem2->instanceType() == Class_Plane)
	rad2 = 0.0;
      
      Point dir1 = elem1->direction();
      Point dir2 = elem2->direction();
      Point loc1 = elem1->location();
      Point loc2 = elem2->location();

      Point v0 = loc2-loc1;
      Point v1 = (v0*dir2)*dir2;
      Point v2 = (v0*dir1)*dir1;
      double d1 = (v0 - v1).length();
      double d2 = (v0 - v2).length();

      double d1_2 = d1 + rad1;
      double d2_2 = d2 + rad2;

      // Check if the surface has approximately the same distance from the
      // identified axis
      Point diff = sf_pt1-sf_pt2;
      if (fabs(d1_2 - d2_2) < frac*d1_2)
	return (diff*vec > 0.0) ? 1 : 2;
      else 
	{
	  // Find the most extreme surface in the outward direction
	  if (elem1->instanceType() == Class_Plane ||
	      elem2->instanceType() == Class_Plane)
	    {
	      double x1 = bpt1*vec;
	      double x2 = bpt2*vec;
	      return (x1 >= x2) ? 3 : 4;
	    }
	  else
	    {
	      Point diff = sf_pt1-sf_pt2;
	      if (diff*vec > 0.0)
		return 3;  // The initially selected surface is the best boundary surface
	      else
		return 4;  // The new surface is better
	    }
	}
    }
  else
    {
      MESSAGE("At most one elementary surface is found");
    }
  return 0;
}

 //==========================================================================
void 
CreateTrimVolume::extractMaxSet(vector<shared_ptr<ftSurface> >& bd_faces,
				vector<shared_ptr<ftSurface> >& trim_faces)
//==========================================================================
{
  // Divide the remaining faces in compact sets and select the "largest"
  // one as source for boundary faces
  // Do not modify already defined neighbourhood information
  tpTolerances tol = model_->getTolerances();
  shared_ptr<SurfaceModel> bd_mod(new SurfaceModel(tol.gap, tol.gap, 
						   tol.neighbour, tol.kink,
						   tol.bend, bd_faces,
						   true));
  vector<shared_ptr<SurfaceModel> > bd_conn = bd_mod->getConnectedModels();
  if (bd_conn.size() > 1)
    {
      vector<BoundingBox> bbox(bd_conn.size());
      for (size_t ki=0; ki<bd_conn.size(); ++ki)
	bbox[ki] = bd_conn[ki]->boundingBox();

      Point high1 = bbox[0].high();
      Point low1 = bbox[0].low();
      double size1 = high1.dist(low1);
      int ix = 0;
      for (size_t kj=1; kj<bd_conn.size(); ++kj)
	{
	  Point high2 = bbox[kj].high();
	  Point low2 = bbox[kj].low();
	  double size2 = high2.dist(low2);
	  if (size2 > size1)
	    {
	      size1 = size2;
	      ix = (int)kj;
	    }
	}

      bd_faces.clear();
      bd_faces = bd_conn[ix]->allFaces();
      Point high = bbox[ix].high();
      Point low = bbox[ix].low();
      int dim = high.dimension();
      for (size_t ki=0; ki<bd_conn.size(); ++ki)
	{
	  if ((int)ki == ix)
	    continue;
	  vector<shared_ptr<ftSurface> > faces = bd_conn[ki]->allFaces();
	  
	  // Check if the bounding box of the small face set is completely
	  // included in the box of the large one
	  Point high2 = bbox[ki].high();
	  Point low2 = bbox[ki].low();
	  int kr=0;
	  for (kr=0; kr<dim; ++kr)
	    {
	      if (low2[kr] < low[kr] || high2[kr] > high[kr])
		break;
	      if (kr < dim)
		trim_faces.insert(trim_faces.end(), faces.begin(), faces.end());
	      else
		bd_faces.insert(bd_faces.end(), faces.begin(), faces.end());
	    }		
	}
    }
}

//==========================================================================
shared_ptr<ftVolume> 
CreateTrimVolume::createTrimVolume(shared_ptr<ParamVolume> vol, 
    vector<pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs)
//==========================================================================
{
  shared_ptr<ftVolume> trimvol;

  // Define input faces as boundary or trimming faces
  // Replace all surfaces referenced by the faces with volume aware surface
  vector<shared_ptr<ftSurface> > faces = model_->allFaces();
  tpTolerances tol = model_->getTolerances();
  double eps = tol.gap;  // Should be parameter space related

  // All volume boundary surfaces
  vector<shared_ptr<ParamSurface> > vol_sfs = vol->getAllBoundarySurfaces();
  const Array<double,6> par_span = vol->parameterSpan();
  if (vol_sfs.size() != 6)
    {
      // Unexpected situation. Return dummy
      return trimvol;
    }

  // Fetch volume knot vectors
  vector<vector<double> > knots(3);
  shared_ptr<SplineVolume> vol2 = dynamic_pointer_cast<SplineVolume>(vol);
  for (int kr=0; kr<3; ++kr)
    {
      const BsplineBasis& basis = vol2->basis(kr);
      basis.knotsSimple(knots[kr]);
    }

  for (size_t ki=0; ki<faces.size(); ++ki)
    {
      // Check if any volume iso-parameter information exist
      // Initially it is set as non existing
      int boundary = -1;
      int constdir = 0;
      double constpar = 0.0;
      bool swapped = false;
      size_t kj=0;
      for (kj=0; kj<side_sfs.size(); ++kj)
	{
	  if (faces[ki].get() == side_sfs[kj].first.get())
	    {
	      double u, v;
	      Point face_pt = 
		faces[ki]->surface()->getInternalPoint(u, v);
	      Point face_norm = faces[ki]->normal(u,v);

	      double sf_dist = HUGE;
	      int sf_ix = -1;
	      for (size_t kr=0; kr<vol_sfs.size(); ++kr)
		{
		  double upar, vpar, dist;
		  Point clo_pt;
		  vol_sfs[kr]->closestPoint(face_pt, upar, vpar,
					    clo_pt, dist, tol.gap);
		  Point sf_norm; 
		  vol_sfs[kr]->normal(sf_norm, upar, vpar);
		  if (dist < sf_dist)
		    {
		      sf_dist = dist;
		      sf_ix = (int)kr;
		      if (face_norm*sf_norm < 0.0)
			swapped = true;
		    }
		}
	      // We know that we have a spline volume. Then the sequence of
	      // boundary surfaces is: umin, umax, vmin, vmax, wmin, wmax
	      boundary = sf_ix;
	      constdir = (sf_ix/2) + 1;
	      constpar = par_span[boundary];
		
	      break;
	    }
	}

      shared_ptr<ParamSurface> surf = faces[ki]->surface();
      shared_ptr<BoundedSurface> bd_surf = 
	dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
      if (bd_surf.get())
	{
	  surf = bd_surf->underlyingSurface();

	  // Check if the size of the underlying surface is
	  // too huge compared to the size of the bounded surface
	  RectDomain dom1 = bd_surf->containingDomain();
	  RectDomain dom2 = surf->containingDomain();
	  double diag1 = dom1.diagLength();
	  double diag2 = dom2.diagLength();
	  double fac = 10.0;
	  if (diag2 > fac*diag1)
	    {
	      // Reduce size of underlying surface
	      double umin = std::max(dom2.umin(), dom1.umin()-0.1*diag1);
	      double umax = std::min(dom2.umax(), dom1.umax()+0.1*diag1);
	      double vmin = std::max(dom2.vmin(), dom1.vmin()-0.1*diag1);
	      double vmax = std::min(dom2.vmax(), dom1.vmax()+0.1*diag1);
	      
	      vector<shared_ptr<ParamSurface> > sub_sfs = 
		surf->subSurfaces(umin, vmin, umax, vmax);
	      if (sub_sfs.size() == 1)
		{
		  surf = sub_sfs[0];
		  bd_surf->replaceSurf(surf);
		}
	    }
	}

      if (boundary < 0)
	{
#ifdef DEBUG
	  std::ofstream of("curr_trim_face.g2");
	  faces[ki]->surface()->writeStandardHeader(of);
	  faces[ki]->surface()->write(of);
#endif
	  // Check if the surface is iso-parametric and corresponds to a knot 
	  // Initial check
	  double u_inner, v_inner;
	  Point pt_inner = faces[ki]->surface()->getInternalPoint(u_inner, v_inner);

	  // Find volume parameter
	  double par[3];
	  double dd;
	  Point clo;
	  vol->closestPoint(pt_inner, par[0], par[1], par[2], clo, dd, tol.gap);

	  // Check if this parameter value coincides with a knot in any parameter direction
	  for (int kr=0; kr<3; ++kr)
	    {
	      size_t kj;
	      for (kj=0; kj<knots[kr].size(); ++kj)
		{
		  if (fabs(knots[kr][kj] - par[kr]) < eps)
		    {
		      bool coinc = checkIsoPar(faces[ki]->surface(), vol, kr, par[kr], eps);
		      if (coinc)
			{
			  constdir = kr + 1;
			  constpar = par[kr];
			  if (kj == 0 || kj == knots[kr].size()-1)
			    {
			      // Also a boundary surface
			      boundary = 2*kr + (kj == knots[kr].size()-1);
			    }
			  break;
			}
		    }
		}
	      if (kj < knots[kr].size())
		break;
	    }
	}

      // Create surface with volume relation information
      shared_ptr<ParamSurface> parsurf; // Dummy
      shared_ptr<SurfaceOnVolume> vol_sf(new SurfaceOnVolume(vol, surf,
							     parsurf, false,
							     constdir, constpar,
							     boundary, swapped));

      // Replace surface
      if (bd_surf.get())
	{
	  bd_surf->replaceSurf(vol_sf);
	}
      else
	{
	  faces[ki]->replaceSurf(vol_sf);
	}
    }

  // Create ftVolume
  trimvol = shared_ptr<ftVolume>(new ftVolume(vol, model_));
  if (material_ >= 0)
    trimvol->setMaterial(material_);
  return trimvol;
}

//==========================================================================
void 
CreateTrimVolume::identifyInnerTrim(vector<shared_ptr<ftSurface> >& bd_faces,
				    vector<shared_ptr<ftSurface> >& trim_faces)
//==========================================================================
{
  // Collect all edges belonging to inner trimming loops
  vector<shared_ptr<ftEdgeBase> > trim_edg;
  for (size_t ki=0; ki<bd_faces.size(); ++ki)
    {
      int nmb_loops = bd_faces[ki]->nmbBoundaryLoops();
      if (nmb_loops <= 1)
	continue;

      for (int ka=1; ka<nmb_loops; ++ka)
	{
	  shared_ptr<Loop> loop = bd_faces[ki]->getBoundaryLoop(ka);
	  vector<shared_ptr<ftEdgeBase> > tmp_edg = loop->getEdges();
	  trim_edg.insert(trim_edg.end(), tmp_edg.begin(), tmp_edg.end());
	}
    }

  // Classify surfaces with outer boundary loop edges being connected
  // to the inner trim edges as trim faces
  for (size_t ki=0; ki<bd_faces.size(); )
    {
      // Fetch outer boundary loop
      shared_ptr<Loop> loop = bd_faces[ki]->getBoundaryLoop(0);
      int nmb = loop->size();
      int ka = 0;
      for (; ka<nmb; ++ka)
	{
	  shared_ptr<ftEdgeBase> edg = loop->getEdge(ka);
	  size_t kj = 0;
	  for (; kj<trim_edg.size(); ++kj)
	    {
	      if (trim_edg[kj].get() == edg.get() ||
		  trim_edg[kj].get() == edg->twin())
		break;
	    }
	  if (kj < trim_edg.size())
	    {
	      trim_faces.push_back(bd_faces[ki]);
	      bd_faces.erase(bd_faces.begin()+ki);
	      break;
	    }
	}
      if (ka >= nmb)
	++ki;
    }
}
//==========================================================================
void 
CreateTrimVolume::refineInSharpEdges(shared_ptr<ParamVolume>& vol)
//==========================================================================
{
  // First fetch all kinks and G1 discontinuities
  vector<ftEdge*> sharp_edg;
  vector<ftEdge*> kinks;
  model_->getCorners(sharp_edg);
  model_->getKinks(kinks);
  sharp_edg.insert(sharp_edg.end(), kinks.begin(), kinks.end());

  // For each underlying curve, check if it follows an iso-parametric curve
  // in the volume
  double eps = model_->getTolerances().neighbour;
  double angtol1 = model_->getTolerances().kink;
  double angtol2 = model_->getTolerances().bend;
  double ptol = eps; //1.0e-8;
  double fac = 1000.0;
  vector<vector<double> > nknot(3);
  vector<vector<double> > edglen(3);
#ifdef DEBUG
  std::ofstream of("sharp_edges.g2");
#endif

  for (size_t ki=0; ki<sharp_edg.size(); ++ki)
    {
      // Estimate number of sampling points
      double len = sharp_edg[ki]->estimatedCurveLength();
      int nmb = (int)(len/(fac*eps));
      nmb = std::min(50, std::max(5, nmb));

      ftEdgeBase *twin = sharp_edg[ki]->twin();

      // Evaluate curve start point
      double t1 = sharp_edg[ki]->tMin();
      double t2 = sharp_edg[ki]->tMax();
      Point pos1 = sharp_edg[ki]->point(t1);
      Point norm1_1 = sharp_edg[ki]->normal(t1);
      double close_t1, close_dist1;
      Point close1;
      twin->closestPoint(pos1, close_t1, close1, close_dist1);
      Point norm1_2 = twin->normal(close_t1);

      // Find volume parameter
      double par1[3];
      double d1;
      Point clo1;
      vol->closestPoint(pos1, par1[0], par1[1], par1[2], clo1, d1, eps);

      vector<Point> vol_pt1(4);
      vol->point(vol_pt1, par1[0], par1[1], par1[2], 1);

     // Check equality
      double del = (t2 - t1)/(double)(nmb-1);
      double tpar = t1 + del;
      bool isEqual[3];
      isEqual[0] = isEqual[1] = isEqual[2] = true;

      double accpar[3];
      double ang3 = norm1_1.angle(norm1_2);
      for (int kr=0; kr<3; ++kr)
	{
	  double ang1 = norm1_1.angle(vol_pt1[kr+1]);
	  ang1 = std::min(ang1, M_PI-ang1);
	  double ang2 = norm1_2.angle(vol_pt1[kr+1]);
	  ang2 = std::min(ang1, M_PI-ang2);
	  if (ang1 < angtol2 || ang2 < angtol2 ||
	      ((fabs(0.5*M_PI-ang1) < angtol2 || fabs(0.5*M_PI-ang2) < angtol2) && 
	       ang3 > angtol2))
	    isEqual[kr] = false;
	  accpar[kr] = par1[kr];
	}

      for (int kj=1; kj<nmb; ++kj, tpar+=del)
	{
	  Point pos2 = sharp_edg[ki]->point(tpar);

	  Point norm2_1 = sharp_edg[ki]->normal(tpar);
	  double close_t2, close_dist2;
	  Point close2;
	  twin->closestPoint(pos2, close_t2, close2, close_dist2);
	  Point norm2_2 = twin->normal(close_t2);

	  // Find volume parameter
	  double par2[3];
	  double d2;
	  Point clo2;
	  vol->closestPoint(pos2, par2[0], par2[1], par2[2], clo2, d2, eps);

	  vector<Point> vol_pt2(4);
	  vol->point(vol_pt2, par2[0], par2[1], par2[2], 1);

	  double ang3_2 = norm2_1.angle(norm2_2);
	  for (int kr=0; kr<3; ++kr)
	    {
	      double ang1 = norm2_1.angle(vol_pt2[kr+1]);
	      ang1 = std::min(ang1, M_PI-ang1);
	      double ang2 = norm2_2.angle(vol_pt2[kr+1]);
	      ang2 = std::min(ang1, M_PI-ang2);
	      if (fabs(par1[kr] - par2[kr]) >= ptol || 
		  ang1 < angtol2 || ang2 < angtol2 ||
		  ((fabs(0.5*M_PI-ang1) < angtol2 || fabs(0.5*M_PI-ang2) < angtol2) && 
		   ang3_2 > angtol2))
		isEqual[kr] = false;   // Use of geometric tolerance may be questional
	      else
		accpar[kr] += par2[kr];
	    }

	  if (isEqual[0] == false && isEqual[1] == false && isEqual[2] == false)
	    break;   // No constant parameter
	}

      for (int kr=0; kr<3; ++kr)
	{
	  if (isEqual[kr])
	    {
	      nknot[kr].push_back(accpar[kr]/(double)nmb);
	      edglen[kr].push_back(len);
	    }
	}
#ifdef DEBUG
      if (isEqual[0] || isEqual[1] || isEqual[2])
	{
	  shared_ptr<ParamCurve> 
	    cv(sharp_edg[ki]->geomCurve()->geometryCurve()->subCurve(sharp_edg[ki]->tMin(),
								     sharp_edg[ki]->tMax()));
	  cv->writeStandardHeader(of);
	  cv->write(of);
	}
#endif
    }

  // Remove duplicate and very close knots
  shared_ptr<SplineVolume> vol2 = dynamic_pointer_cast<SplineVolume>(vol);
  if (!vol2)
    return;  // Formality

  size_t ki, kj;
  for (int kr=0; kr<3; ++kr)
    {
      if (nknot[kr].size() == 0)
	continue;

      for (ki=0; ki<nknot[kr].size(); ++ki)
	for (kj=ki+1; kj<nknot[kr].size(); ++kj)
	  if (nknot[kr][kj] < nknot[kr][ki])
	    {
	      std::swap(nknot[kr][kj], nknot[kr][ki]);
	      std::swap(edglen[kr][kj], edglen[kr][kj]);
	    }

      double start = vol2->startparam(kr);
      double end = vol2->endparam(kr);
      
      // Check endpoints
      for (ki=0; ki<nknot[kr].size();)
	{
	  if (nknot[kr][ki] - start < ptol)
	    nknot[kr].erase(nknot[kr].begin());
	  else
	    break;
	}

      for (ki=nknot[kr].size()-1; ki>=0;)
	{
	  if (end - nknot[kr][ki] < ptol)
	    {
	      nknot[kr].pop_back();
	      --ki;
	    }
	  else
	    break;
	}

      // Check internal knots
      for (ki=0; ki<nknot[kr].size(); ki=kj)
	{
	  double accknot = edglen[kr][ki]*nknot[kr][ki];
	  double acclen = edglen[kr][ki];
	  for (kj=ki+1; kj<nknot[kr].size(); ++kj)
	    {
	      if (nknot[kr][kj]-nknot[kr][ki] > ptol)
		break;
	      accknot += edglen[kr][kj]*nknot[kr][kj];
	      acclen += edglen[kr][kj];
	    }

	  if (kj - ki > 1)
	    {
	      // Use middle value
	      nknot[kr][ki] = accknot/acclen;
	      nknot[kr].erase(nknot[kr].begin() + ki + 1, nknot[kr].begin() + kj);
	      kj = ki+1;
	    }
	}
    
      // Insert knots
      if (nknot[kr].size() > 0)
	vol2->insertKnot(kr, nknot[kr]);
    }

}

//==========================================================================
bool
CreateTrimVolume::checkIsoPar(shared_ptr<ParamSurface> surf,
			      shared_ptr<ParamVolume> vol,
			      int pardir, double parval, double tol)
//==========================================================================
{
  // Preparatory computations
  // Get estimated length of surface sides
  double len_u, len_v;
  GeometryTools::estimateSurfaceSize(*surf, len_u, len_v);
  
  // Number of points to sample in each parameter direction
  double fac = 1000.0;
  int min_samples = 3;
  int nmb_sample = (int)sqrt((len_u*len_v)/(fac*fac*tol*tol));
  int nmb_u = (int)(nmb_sample*len_u/len_v);
  int nmb_v = (int)(nmb_sample*len_v/len_u);
  nmb_u = std::max(nmb_u, min_samples);
  nmb_v = std::max(nmb_v, min_samples);
  nmb_u = std::min(nmb_u, min_samples*nmb_sample);
  nmb_v = std::min(nmb_v, min_samples*nmb_sample);

  // Fetch constant parameter curves in the u-direction
  RectDomain dom = surf->containingDomain();
  double u1 = dom.umin();
  double u2 = dom.umax();
  double udel = (u2 - u1)/(double)(nmb_u+1);
  double tol1 = std::max(1.0e-7, 1.0e-5*(u2-u1));
  double upar;

  size_t kr;
  upar = u1;
  while (upar <= u2+tol1)
    {
      vector<shared_ptr<ParamCurve> > crvs = surf->constParamCurves(upar, false);
      if (crvs.size() == 0)
	{
	  upar += udel;
	  continue;  // Outside domain of surface
	}

      // Distribute sampling points
      double av_len = 0.0;
      vector<double> cv_len(crvs.size());
      for (kr=0; kr<crvs.size(); ++kr)
	{
	  double len = crvs[kr]->estimatedCurveLength();
	  av_len += len;
	  cv_len[kr] = len;
	}
      double curr_len = av_len;
      av_len /= (double)crvs.size();
      
      // Evaluate sampling points
      int curr_nmb = (int)(nmb_v*(curr_len/len_v)) + 1;
      for (kr=0; kr<crvs.size(); ++kr)
	{
	  int nmb = (int)(curr_nmb*cv_len[kr]/av_len);
	  nmb = std::max(nmb, min_samples);
	  double v1 = crvs[kr]->startparam();
	  double v2 = crvs[kr]->endparam();
	  double vdel = (v2 - v1)/(double)(nmb+1);
	  double tol2 = std::max(1.0e-7, 1.0e-5*(v2-v1));
	  double vpar;
	  vpar = v1;

	  while (vpar <= v2+tol2)
	    {
	      Point pos = crvs[kr]->point(vpar);

	      // Find volume parameter
	      double par[3];
	      double dd;
	      Point clo;
	      vol->closestPoint(pos, par[0], par[1], par[2], clo, dd, tol);
	      if (fabs(par[pardir] - parval) >= tol)
		return false;

	      vpar += vdel;
	    }
	}
      upar += udel;
    }
  return true;
}

//==========================================================================
void CreateTrimVolume::repairShell(int degree)
//==========================================================================
{
  double epsge = model_->getTolerances().gap;

  try {
  // Merge adjacent spline surfaces with g1 continuity
  model_->simplifyShell();
  }
  catch (...)
    {
      ;
    }

  // Further model simplification by approximating surfaces patches
  // that compose a four sided domain by one large surface
  try {
  SurfaceModelUtils::simplifySurfaceModel(model_, degree);
  }
  catch (...)
    {
      ;
    }

  // // Further model simplification by approximating surfaces patches
  // // that compose a four sided domain by one large surface
  // SurfaceModelUtils::simplifySurfaceModel(model_, degree);

  bool gap_trim = false;
  if (gap_trim)
    {
      // Regenerate trimming curves to reduce gaps between surfaces
      // First identify occurances
      vector<shared_ptr<ftEdge> > gap_edges = model_->getUniqueInnerEdges();
      //model_->getGaps(gap_edges);

      // Make sure that the corresponding vertex positions are updated
      std::set<shared_ptr<Vertex> > vertices;
      for (size_t ki=0; ki<gap_edges.size(); ++ki)
	{
	  shared_ptr<Vertex>  vx1, vx2;
	  gap_edges[ki]->getVertices(vx1, vx2);
	  vertices.insert(vertices.end(), vx1);
	  vertices.insert(vertices.end(), vx2);
	}

      std::set<shared_ptr<Vertex> >::iterator it = vertices.begin();
      while (it != vertices.end())
	{
	  (*it)->updateVertexPos(epsge);
	  ++it;
	}
  
      // Update trimming curves
      for (size_t ki=0; ki<gap_edges.size(); ++ki)
	{
	  // Fetch the related faces
	  ftSurface *face1 = gap_edges[ki]->face()->asFtSurface();
	  ftSurface *face2 = gap_edges[ki]->twin()->face()->asFtSurface();

	  if (!(face1 && face2))
	    continue;   // This should not happen. Cannot mend gaps

	  // Check edge curves
	  ftEdge* e1g = gap_edges[ki].get();
	  ftEdge* e2g = gap_edges[ki]->twin()->geomEdge();
	  shared_ptr<CurveOnSurface> sfcv1 =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	    (e1g->geomCurve());
	  shared_ptr<CurveOnSurface> sfcv2 =
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>
	    (e2g->geomCurve());
	  if (!(sfcv1.get() && sfcv2.get()))
	    continue;

	  bool same1 = sfcv1->sameCurve(epsge);
	  bool same2 = sfcv2->sameCurve(epsge);
	  shared_ptr<FaceConnectivity<ftEdgeBase> > info = 
	    gap_edges[ki]->getConnectivityInfo();
	  int kj=-1;
	  if (info.get())
	    {
	      for (kj=0; kj<info->status_.size(); ++kj)
		if (info->status_[kj] >= 3)
		  break;
	    }
	  if (same1 && same2 && info.get() && info.get() && 
	      kj >= info->status_.size())
	    continue;

	  // Check if both faces are trimmed. In that case a better positioning
	  // of the trim curve may remove the gap
	  shared_ptr<BoundedSurface> bd1 = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(face1->surface());
	  shared_ptr<BoundedSurface> bd2 = 
	    dynamic_pointer_cast<BoundedSurface, ParamSurface>(face2->surface());

	  if (!(bd1.get() && bd2.get()))
	    {
	      continue;   // Not handled 
	    }
#ifdef DEBUG
	  std::ofstream of("gap_sfs.g2");
	  bd1->writeStandardHeader(of);
	  bd1->write(of);
	  bd2->writeStandardHeader(of);
	  bd2->write(of);
#endif

	  Point vertex1 = e1g->getVertex(true)->getVertexPoint();
	  Point vertex2 = e1g->getVertex(false)->getVertexPoint();
	  double max_gap = 
	    GapRemoval::removeGapTrim(sfcv1, e1g->tMin(), e1g->tMax(),
				      sfcv2, e2g->tMin(), e2g->tMax(),
				      vertex1, vertex2, epsge);
	  int stop_break = 1;
	}
    }
   
}

//==========================================================================
void
  CreateTrimVolume::analyzePrio(int* prio, int nmb_pri, 
				Point coord[], Point coord_pos[])
//==========================================================================
{
  // Start by looking for rotational surfaces and identify rotational 
  // axis
  vector<int> rot;
  vector<int> planar;
  vector<Point> axis;
  vector<Point> dir;
  vector<vector<Point> > dir_pos;
  vector<int> nmb_axis, nmb_dir;
  vector<double> dir_ang;
  double lim_ang = M_PI/12.0;
  double lim_ang2 = 0.5*lim_ang;
  double lim_dist = 0.1*(bigbox_.low().dist(bigbox_.high()));
  for (int ki=0; ki<nmb_pri; ++ki)
    {
      if (sf_type_[prio[ki]] == ROTATIONAL)
	{
	  rot.push_back(prio[ki]);
	  if (axis.size() == 0)
	    {
	      axis.push_back(sf_axis_[prio[ki]]);
	      nmb_axis.push_back(1);
	    }
	  else
	    {
	      size_t kj;
	      Point tmp_axis;
	      for (kj=0; kj<axis.size(); ++kj)
		{
		  tmp_axis = sf_axis_[prio[ki]];
		  if (axis[kj]*tmp_axis < 0.0)
		    tmp_axis *= -1;
		  double ang = axis[kj].angle(tmp_axis);
		  if (ang < lim_ang)
		    break;
		}
	      if (kj < axis.size())
		{
		  axis[kj] += tmp_axis;
		  nmb_axis[kj]++;
		}
	      else
		{
		  axis.push_back(sf_axis_[prio[ki]]);
		  nmb_axis.push_back(1);
		}
	    }
	}
      else if (sf_type_[prio[ki]] == PLANAR)
	{
	  planar.push_back(prio[ki]);
	}
      if (sf_type_[prio[ki]] == ROTATIONAL || sf_type_[prio[ki]] == PLANAR)
	{
	  size_t kj = dir.size();
	  if (dir.size() > 0)
	    {
	      for (kj=0; kj<dir.size(); ++kj)
		{
		  double ang = dir[kj].angle(cone_[prio[ki]].centre());
		  double min_pt_ang = HUGE;
		  double max_pt_ang = 0.0;
		  double min_pt_dist = HUGE;
		  double max_pt_dist = 0.0;
		  for (size_t kr=0; kr<dir_pos[kj].size(); ++kr)
		    {
		      Point tmp_vec = dir_pos[kj][kr]-sf_pt_[prio[ki]];
		      double tmp_ang = tmp_vec.angle(dir[kj]);
		      double tmp_dist = dir_pos[kj][kr].dist(sf_pt_[prio[ki]]);
		      min_pt_ang = std::min(min_pt_ang, tmp_ang);
		      max_pt_ang = std::max(max_pt_ang, tmp_ang);
		      min_pt_dist = std::min(min_pt_dist, tmp_dist);
		      max_pt_dist = std::max(max_pt_dist, tmp_dist);
		    }
		  if (ang < lim_ang)
		    break;
		}
	    }
	  if (kj < dir.size())
	    {
	      dir[kj] += cone_[prio[ki]].centre();
	      nmb_dir[kj]++;
	      dir_pos[kj].push_back(sf_pt_[prio[ki]]);
	      dir_ang[kj] = std::max(dir_ang[kj], cone_[prio[ki]].angle());
	    }
	  else
	    {
	      dir.push_back(cone_[prio[ki]].centre());
	      nmb_dir.push_back(1);
	      vector<Point> curr_dir_pos(1, sf_pt_[prio[ki]]);
	      dir_pos.push_back(curr_dir_pos);
	      dir_ang.push_back(cone_[prio[ki]].angle());
	    }
	}
    }
  for (size_t kj=0; kj<axis.size(); ++kj)
    axis[kj] /= (double)nmb_axis[kj];
  for (size_t kj=0; kj<dir.size(); ++kj)
    dir[kj] /= (double)nmb_dir[kj];

  if (axis.size() > 0)
    {
      // Prioritize directions that are roughly parallel to or perpendicular
      // to the most important axis
      size_t ix = 0;
      size_t kr;
      for (kr=1; kr<axis.size(); ++kr)
	{
	  if (nmb_axis[kr] > nmb_axis[ix])
	    ix = kr;
	}
      for (int kh=0; kh<6; kh+=2)
	{
	  for (kr=0; kr<dir.size(); )
	    {
	      // Check against previous
	      int kk;
	      for (kk=0; kk<kh; ++kk)
		{
		  double prev_ang = dir[kr].angle(coord[kk]);
		  if (prev_ang < lim_ang ||
		      fabs(M_PI-prev_ang) < lim_ang)
		    break;
		}
	      if (kk < kh)
		{
		  ++kr;
		  continue;
		}

	      double ang = dir[kr].angle(axis[ix]);
	      if (ang < lim_ang2 || M_PI-ang < lim_ang || 
		  fabs(0.5*M_PI-ang) < lim_ang)
		{
		  coord[kh] = dir[kr];
		  coord[kh].normalize();
		  double mm = dir_pos[kr][0]*dir[kr];
		  coord_pos[kh] = dir_pos[kr][0];
		  for (size_t ka=1; ka<dir_pos[kr].size(); ++ka)
		    {
		      double mm2 = dir_pos[kr][ka]*dir[kr];
		      if (mm2 > mm)
			{
			  mm = mm2;
			  coord_pos[kh] = dir_pos[kr][ka];
			}
		    }
		  dir.erase(dir.begin()+kr);
		  dir_pos.erase(dir_pos.begin()+kr);
		  nmb_dir.erase(nmb_dir.begin()+kr);

		  // Find opposite direction
		  double ang1, ang2;
		  Point curr_pos;
		  int ix2 = -1;
		  for (size_t ka=0; ka<dir.size(); ++ka)
		    {
		      double ang3 = dir[ka].angle(axis[ix]);
		      double ang4 = dir[ka].angle(coord[kh]);
		      double mm2 = dir_pos[ka][0]*dir[ka];
		      Point tmp_pos = dir_pos[ka][0];
		      for (size_t kb=1; kb<dir_pos[ka].size(); ++kb)
			{
			  double mm3 = dir_pos[ka][kb]*dir[ka];
			  if (mm3 > mm2)
			    {
			      mm2 = mm3;
			      tmp_pos = dir_pos[ka][kb];
			    }
			}
		      double dd = coord_pos[kh].dist(tmp_pos);
		      if (fabs(ang - ang3) < lim_ang &&
			  ang4 > lim_ang && dd > lim_dist)
			{
			  if (ix2 < 0 || ang4 > ang2)
			    {
			      ix2 = (int)ka;
			      ang1 = ang3;
			      ang2 = ang4;
			      curr_pos = tmp_pos;
			    }
			}
		    }

		  if (ix2 >= 0)
		    {
		      coord[kh+1] = dir[ix2];
		      coord[kh+1].normalize();
		      coord_pos[kh+1] = curr_pos;
		      dir.erase(dir.begin()+ix2);
		      dir_pos.erase(dir_pos.begin()+ix2);
		      nmb_dir.erase(nmb_dir.begin()+ix2);
		    }
		  if (coord[kh+1].dimension() == 0)
		    coord[kh+1] = -coord[kh];
		  break;
		}
	      else
		++kr;
	    }
	}
    }
  else if (dir.size() > 0)
    {
      // Prioritize directions that are roughly parallel to or perpendicular
      // to the most important plane normal
      size_t ix = 0;
      size_t kr;
      for (kr=1; kr<dir.size(); ++kr)
	{
	  if (nmb_dir[kr] > nmb_dir[ix])
	    ix = kr;
	}

      Point start_dir = dir[ix];
      for (int kh=0; kh<6; kh+=2)
	{
	  for (kr=0; kr<dir.size(); )
	    {
	      // Check against previous
	      int kk;
	      for (kk=0; kk<kh; ++kk)
		{
		  double prev_ang = dir[kr].angle(coord[kk]);
		  if (prev_ang < lim_ang ||
		      fabs(M_PI-prev_ang) < lim_ang)
		    break;
		}
	      if (kk < kh)
		{
		  ++kr;
		  continue;
		}

	      double ang = dir[kr].angle(start_dir);
	      if (ang < lim_ang2 || M_PI-ang < lim_ang || 
		  fabs(0.5*M_PI-ang) < lim_ang)
		{
		  coord[kh] = dir[kr];
		  coord[kh].normalize();
		  double mm = dir_pos[kr][0]*dir[kr];
		  coord_pos[kh] = dir_pos[kr][0];
		  for (size_t ka=1; ka<dir_pos[kr].size(); ++ka)
		    {
		      double mm2 = dir_pos[kr][ka]*dir[kr];
		      if (mm2 > mm)
			{
			  mm = mm2;
			  coord_pos[kh] = dir_pos[kr][ka];
			}
		    }
		  dir.erase(dir.begin()+kr);
		  dir_pos.erase(dir_pos.begin()+kr);
		  nmb_dir.erase(nmb_dir.begin()+kr);

		  // Find opposite direction
		  double ang1, ang2;
		  Point curr_pos;
		  int ix2 = -1;
		  for (size_t ka=0; ka<dir.size(); ++ka)
		    {
		      double ang3 = dir[ka].angle(start_dir);
		      double ang4 = dir[ka].angle(coord[kh]);
		      double mm2 = dir_pos[ka][0]*dir[ka];
		      Point tmp_pos = dir_pos[ka][0];
		      for (size_t kb=1; kb<dir_pos[ka].size(); ++kb)
			{
			  double mm3 = dir_pos[ka][kb]*dir[ka];
			  if (mm3 > mm2)
			    {
			      mm2 = mm3;
			      tmp_pos = dir_pos[ka][kb];
			    }
			}
		      double dd = coord_pos[kh].dist(tmp_pos);
		      if (fabs(ang - ang3) < lim_ang &&
			  ang4 > lim_ang && dd > lim_dist)
			{
			  if (ix2 < 0 || ang4 > ang2)
			    {
			      ix2 = (int)ka;
			      ang1 = ang3;
			      ang2 = ang4;
			      curr_pos = tmp_pos;
			    }
			}
		    }

		  if (ix2 >= 0)
		    {
		      coord[kh+1] = dir[ix2];
		      coord[kh+1].normalize();
		      coord_pos[kh+1] = curr_pos;
		      dir.erase(dir.begin()+ix2);
		      dir_pos.erase(dir_pos.begin()+ix2);
		      nmb_dir.erase(nmb_dir.begin()+ix2);
		    }
		  if (coord[kh+1].dimension() == 0)
		    coord[kh+1] = -coord[kh];
		  break;
		}
	      else
		++kr;
	    }
	}
    }
	  
  int stop_break = 1;
}
