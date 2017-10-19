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

//#define DEBUG_VOL
//#define DEBUG

#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/Body.h"
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/topology/FaceAdjacency.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/creators/CurveCreators.h"
#include <fstream>

using std::vector;

using namespace Go;

//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> >
ftVolumeTools::splitVolumes(shared_ptr<ftVolume>& vol1, 
			    shared_ptr<ftVolume>& vol2, double eps,
			    vector<int>& config)
//===========================================================================
{
  vector<shared_ptr<ftVolume> >  result;
  config.clear();

  // Fetch all boundary surfaces
  vector<shared_ptr<SurfaceModel> > shells1 = vol1->getAllShells();
  vector<shared_ptr<SurfaceModel> > shells2 = vol2->getAllShells();
#ifdef DEBUG_VOL
  int nmb_bd1 = shells1[0]->nmbBoundaries();
  int nmb_bd2 = shells2[0]->nmbBoundaries();
  std::cout << "SplitVolumes, init. Number of boundaries: " << nmb_bd1 << ", " << nmb_bd2 << std::endl;
#endif

  if (shells1.size() == 0 || shells2.size() == 0)
    return result;

  size_t ki;
  int kj;
  shared_ptr<SurfaceModel> sfmodel1 = shells1[0];
  for (ki=1; ki<shells1.size(); ++ki)
      sfmodel1->append(shells1[ki]);
  
  shared_ptr<SurfaceModel> sfmodel2 = shells2[0];
  for (ki=1; ki<shells2.size(); ++ki)
      sfmodel2->append(shells2[ki]);

  // Perform boolean operation on boundary surfaces
  vector<shared_ptr<SurfaceModel> > split_models = 
    sfmodel1->splitSurfaceModels(sfmodel2);

#ifdef DEBUG_VOL
  // Debug 
  if (split_models[0].get())
    {
      std::ofstream of("split_1.g2");
      int nmb = split_models[0]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[0]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  if (split_models[1].get())
    {
      std::ofstream of("split_2.g2");
      int nmb = split_models[1]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[1]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  if (split_models[2].get())
    {
      std::ofstream of("split_3.g2");
      int nmb = split_models[2]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[2]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }

  if (split_models[3].get())
    {
      std::ofstream of("split_4.g2");
      int nmb = split_models[3]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[3]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
#endif

  // For each boundary surface internal to both volumes: make a copy, swap
  // the parameter directions of this copy, and add it to the outside
  // part of the other volume.
  // Set also twin pointers between corresponding boundary surfaces
#ifdef DEBUG_VOL
  std::cout << "Model 3: " << split_models[3]->nmbEntities() << std::endl;
#endif
  int nmb = (!split_models[0].get()) ? 0 : split_models[0]->nmbEntities();
  for (kj=0; kj<nmb; ++kj)
    {
      shared_ptr<ftSurface> face1 = split_models[0]->getFace(kj);
      if (face1->twin())
	continue;
      shared_ptr<ParamSurface> surf1 = split_models[0]->getSurface(kj);

#ifdef DEBUG_VOL
      std::ofstream of0("test_swap.g2");
      surf1->writeStandardHeader(of0);
      surf1->write(of0);
#endif

      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();

      shared_ptr<SurfaceOnVolume> vol_sf = 
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf2);
      if (!vol_sf.get())
	{
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
	  vol_sf =
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf->underlyingSurface());
	}

      if (vol_sf.get())
	{
	  // Unset parameter surface information. Not valid.
	  vol_sf->unsetParamSurf();
	}

#ifdef DEBUG_VOL
      surf1->writeStandardHeader(of0);
      surf1->write(of0);

      std::ofstream of("curr_face.g2");
      surf2->writeStandardHeader(of);
      surf2->write(of);
#endif

      shared_ptr<ftSurface> face2(new ftSurface(surf2, -1));
      face2->setBody(vol2.get());
      if (split_models[3].get())
	{
	  split_models[3]->append(face2);
	  face1->connectTwin(face2.get(), eps);
	}
      int stop_break;
      stop_break = 1;
    }
#ifdef DEBUG_VOL
  std::cout << "Model 3(2): " << split_models[3]->nmbEntities() << std::endl;
#endif

  nmb = (!split_models[2].get()) ? 0 : split_models[2]->nmbEntities();
  for (kj=0; kj<nmb; ++kj)
    {
      shared_ptr<ftSurface> face1 = split_models[2]->getFace(kj);
      if (face1->twin())
	continue;
      shared_ptr<ParamSurface> surf1 = split_models[2]->getSurface(kj);
      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();

      shared_ptr<SurfaceOnVolume> vol_sf = 
	dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf2);
      if (!vol_sf.get())
	{
	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
	  vol_sf =
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf->underlyingSurface());
	}

      if (vol_sf.get())
	{
	  // Unset parameter surface information. Not valid.
	  vol_sf->unsetParamSurf();
	}

#ifdef DEBUG_VOL
      std::ofstream of("curr_face.g2");
      surf2->writeStandardHeader(of);
      surf2->write(of);
#endif

      shared_ptr<ftSurface> face2(new ftSurface(surf2, -1));
      face2->setBody(vol1.get());
      if (split_models[1].get())
	{
	  split_models[1]->append(face2);
	  face1->connectTwin(face2.get(), eps);
	}
      int stop_break;
      stop_break = 1;
    }

#ifdef DEBUG_VOL
  // Debug 
  if (split_models[0].get())
    {
      std::ofstream of("split_5.g2");
      int nmb = split_models[0]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[0]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  if (split_models[1].get())
    {
      std::ofstream of("split_6.g2");
      int nmb = split_models[1]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[1]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  if (split_models[2].get())
    {
      std::ofstream of("split_7.g2");
      int nmb = split_models[2]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[2]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }

  if (split_models[3].get())
    {
      std::ofstream of("split_8.g2");
      int nmb = split_models[3]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[3]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
#endif

   // Separate outer surface models into connected sets and make 
  // corresponding volumes
  vector<shared_ptr<SurfaceModel> > sep_models;
  if (split_models[1])
    {
    sep_models =split_models[1]->getConnectedModels();
    for (ki=0; ki<sep_models.size(); ++ki)
      {
	shared_ptr<ftVolume> curr(new ftVolume(vol1->getVolume(), 
					       sep_models[ki]));
	result.push_back(curr);
	config.push_back(1);
      }
    }
  
  if (split_models[3])
    {
      sep_models.clear();
      sep_models = split_models[3]->getConnectedModels();
      for (ki=0; ki<sep_models.size(); ++ki)
	{
	  shared_ptr<ftVolume> curr(new ftVolume(vol2->getVolume(), 
						 sep_models[ki]));
	  result.push_back(curr);
	  config.push_back(2);
	}
    }

  // Merge inner surface models, then split in connected sets and make volumes
  if (split_models[0])
    {
      if (split_models[2])
	split_models[0]->append(split_models[2]);
      sep_models.clear();
      sep_models = split_models[0]->getConnectedModels();
      
      for (ki=0; ki<sep_models.size(); ++ki)
	{
	  shared_ptr<ftVolume> curr(new ftVolume(vol1->getVolume(), 
						 sep_models[ki]));
	  result.push_back(curr);
	  config.push_back(3);
	}
    }

#ifdef DEBUG_VOL
  for (ki=0; ki<result.size(); ++ki)
    {
      std::ofstream res("res_vols.g2");
      shared_ptr<SurfaceModel> mod = result[ki]->getOuterShell();
      int nmb = mod->nmbEntities();
      for (int kv=0; kv<nmb; ++kv)
	{
	  shared_ptr<ParamSurface> sf = mod->getSurface(kv);
	  sf->writeStandardHeader(res);
	  sf->write(res);
	}

      int nmb_bd = mod->nmbBoundaries();
      std::cout << "SplitVolumes, exit. ki= " << ki << ": " << nmb_bd << std::endl;
      
      int stop_break = 1;
    }
#endif
	  

  // TEST
#ifdef DEBUG_VOL
  std::ofstream t1("missing_twin0.g2");
  for (ki=0; ki<result.size(); ++ki)
    {
      shared_ptr<SurfaceModel> bd = result[ki]->getOuterShell();
      int nmb = bd->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> curr = bd->getFace(kj);
	  vector<shared_ptr<ftEdge> > edges = curr->getAllEdges();
	  for (size_t ix=0; ix<edges.size(); ++ix)
	    if (!edges[ix]->twin())
	      {
		std::cout << "Missing twin pointer. SplitVolumes" << ki << std::endl;
		t1 << "410 1 0 4 155 55 0 255" << std::endl;
		t1 << "1" << std::endl;
		t1 << edges[ix]->point(edges[ix]->tMin()) << " " << edges[ix]->point(edges[ix]->tMax()) << std::endl;
	      }
	}
    }
#endif

  return result;
}


//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> >
ftVolumeTools::splitVolumes(shared_ptr<ftVolume>& vol, 
			    shared_ptr<ftSurface>& face, double eps,
			    int create_all)
//===========================================================================
{
  vector<shared_ptr<ftVolume> > result;

  // Fetch all boundary surfaces
  vector<shared_ptr<SurfaceModel> > shells = vol->getAllShells();

  if (shells.size() == 0)
    return result;

  size_t ki;
  int kj;
  shared_ptr<SurfaceModel> sfmodel1 = shells[0];
  for (ki=1; ki<shells.size(); ++ki)
      sfmodel1->append(shells[ki]);
  
  vector<shared_ptr<ftSurface> > face_vec;
  face_vec.push_back(face);
  shared_ptr<SurfaceModel> sfmodel2 = 
    shared_ptr<SurfaceModel>(new SurfaceModel(face_vec, eps));

  // Perform boolean operation on boundary surfaces
  vector<shared_ptr<SurfaceModel> > split_models = 
    sfmodel1->splitSurfaceModels(sfmodel2);

#ifdef DEBUG_VOL
  // Debug 
  if (split_models[0].get())
    {
      std::ofstream of("split_1.g2");
      int nmb = split_models[0]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[0]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  if (split_models[1].get())
    {
      std::ofstream of("split_2.g2");
      int nmb = split_models[1]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[1]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  if (split_models[2].get())
    {
      std::ofstream of("split_3.g2");
      int nmb = split_models[2]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[2]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }

  if (split_models[3].get())
    {
      std::ofstream of("split_4.g2");
      int nmb = split_models[3]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[3]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
#endif

  // For each surface from the initial face internal to the inital volume: make a copy, swap
  // the parameter directions of this copy, and add it to the outside
  // part of the other volume.
  // Set also twin pointers between corresponding boundary surfaces
#ifdef DEBUG_VOL
  std::cout << "Model 3(2): " << split_models[3]->nmbEntities() << std::endl;
#endif

  int nmb = (!split_models[2].get()) ? 0 : split_models[2]->nmbEntities();
  for (kj=0; kj<nmb; ++kj)
    {
      shared_ptr<ftSurface> face1 = split_models[2]->getFace(kj);
      face1->setBody(vol.get());
      if (face1->twin())
	continue;
      split_models[0]->append(face1, true, false, true);
      shared_ptr<ParamSurface> surf1 = split_models[2]->getSurface(kj);
      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();


#ifdef DEBUG_VOL
      std::ofstream of("curr_face.g2");
      surf2->writeStandardHeader(of);
      surf2->write(of);
#endif

      shared_ptr<ftSurface> face2(new ftSurface(surf2, -1));
      face2->setBody(vol.get());
      split_models[1]->append(face2, true, false, true);
      int ix = split_models[1]->getIndex(face2);
      if (ix >= 0 && create_all == 3)
	face1->connectTwin(face2.get(), eps);
      int stop_break;
      stop_break = 1;
    }

#ifdef DEBUG_VOL
  // Debug 
  if (split_models[0].get())
    {
      std::ofstream of("split_5.g2");
      int nmb = split_models[0]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[0]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  if (split_models[1].get())
    {
      std::ofstream of("split_6.g2");
      int nmb = split_models[1]->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = split_models[1]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }

#endif

   // Separate surface models into connected sets and make 
  // corresponding volumes
  vector<shared_ptr<SurfaceModel> > sep_models;
  if (split_models[0] && (create_all == 1 || create_all == 3))
    {
    sep_models =split_models[0]->getConnectedModels();
    for (ki=0; ki<sep_models.size(); ++ki)
      {
	shared_ptr<ftVolume> curr(new ftVolume(vol->getVolume(), 
					       sep_models[ki]));
	if (vol->hasMaterialInfo())
	  curr->setMaterial(vol->getMaterial());
	result.push_back(curr);
      }
    }
  
  if (split_models[1] && (create_all == 2 || create_all == 3))
    {
      sep_models.clear();
      sep_models = split_models[1]->getConnectedModels();
      for (ki=0; ki<sep_models.size(); ++ki)
	{
	  shared_ptr<ftVolume> curr(new ftVolume(vol->getVolume(), 
						 sep_models[ki]));
	  if (vol->hasMaterialInfo())
	    curr->setMaterial(vol->getMaterial());
	  result.push_back(curr);
	}
    }

#ifdef DEBUG_VOL
  // TEST
  std::ofstream t1("missing_twin0.g2");
  for (ki=0; ki<result.size(); ++ki)
    {
      shared_ptr<SurfaceModel> bd = result[ki]->getOuterShell();
      int nmb = bd->nmbEntities();
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> curr = bd->getFace(kj);
	  vector<shared_ptr<ftEdge> > edges = curr->getAllEdges();
	  for (size_t ix=0; ix<edges.size(); ++ix)
	    if (!edges[ix]->twin())
	      {
		std::cout << "Missing twin pointer. SplitVolumes" << ki << std::endl;
		t1 << "410 1 0 4 155 55 0 255" << std::endl;
		t1 << "1" << std::endl;
		t1 << edges[ix]->point(edges[ix]->tMin()) << " " << edges[ix]->point(edges[ix]->tMax()) << std::endl;
	      }
	}
    }
#endif

  return result;
}

//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> >
ftVolumeTools::splitOneVol(shared_ptr<ftVolume>& elem_vol, ftVolume* trim_vol,
			   double eps, vector<int>& is_inside, double* elem_par,
			   int nmb_par)
//===========================================================================
{
  vector<shared_ptr<ftVolume> > result;
  is_inside.clear();

  // Fetch trim faces from boundary shells
  vector<shared_ptr<SurfaceModel> > shells = trim_vol->getAllShells();
  vector<shared_ptr<ftSurface> > faces;
  for (size_t ki=0; ki<shells.size(); ++ki)
    {
      int nmb = shells[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> curr_face = shells[ki]->getFace(kj);
	  if (elem_par != NULL)
	    {
	      // Check if the current face coincides with a constant parameter
	      // of the volume to split
	      shared_ptr<ParamSurface> surf = curr_face->surface();
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
	      int kr = 0;
	      for (kr=0; kr<nmb_par; ++kr)
		{
		  if (dir == (kr/2)+1 && fabs(val-elem_par[kr]) < eps)
		    break;  // Coincidence. Do not include face in intersection
		}
	      if (kr < nmb_par)
		continue;
	    }
	  int bd_stat = ftVolumeTools::boundaryStatus(trim_vol,
						      curr_face, eps);
	  if (bd_stat < 0)
	    faces.push_back(curr_face);
	}
    }

  // Fetch shell surrounding element volume
  shared_ptr<SurfaceModel> elem_shell = elem_vol->getOuterShell();

#ifdef DEBUG
  std::ofstream of1("trim_faces.g2");
  for (size_t kr=0; kr<faces.size(); ++kr)
    {
      shared_ptr<ParamSurface> tmp_sf = faces[kr]->surface();
      tmp_sf->writeStandardHeader(of1);
      tmp_sf->write(of1);
    }
  std::ofstream of2("elem_faces.g2");
  int nmb_elem = elem_shell->nmbEntities();
  for (int kh=0; kh<nmb_elem; ++kh)
    {
      shared_ptr<ParamSurface> tmp_sf = elem_shell->getSurface(kh);
      tmp_sf->writeStandardHeader(of2);
      tmp_sf->write(of2);
    }
#endif

  // Split element volume and trim faces
  vector<vector<shared_ptr<ParamSurface> > > split_groups;
  elem_shell->splitSurfaceModel(faces, trim_vol, split_groups);

  // Remove outdated parameter information
  //for (size_t kr=0; kr<1; ++kr)
  for (size_t kr=2; kr<split_groups.size(); ++kr)
    {
      for (size_t ki=0; ki<split_groups[kr].size(); ++ki)
	{
	  shared_ptr<ParamSurface> surf = split_groups[kr][ki];
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf);
	  if (!vol_sf.get())
	    {
	      shared_ptr<BoundedSurface> bd_sf = 
		dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
	      vol_sf =
		dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf->underlyingSurface());
	    }
	  if (vol_sf.get())
	    {
	      // Unset parameter surface information. Not necessarily valid.
	      vol_sf->unsetParamSurf();
	      vol_sf->setVolume(elem_vol->getVolume());
	    }
	}
    }

   // Add trimming face pieces to the element surfaces
  for (size_t ki=0; ki<split_groups[2].size(); ++ki)
    {
      // Make oppositely oriented copy
      shared_ptr<ParamSurface> surf1 = split_groups[2][ki];
      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();

      // Add surface to element groups
      split_groups[0].push_back(surf1);
      split_groups[1].push_back(surf2);
    }
  
  // Create surface models
  vector<shared_ptr<SurfaceModel> > surf_mod;
  vector<int> inside;
  tpTolerances toptol = shells[0]->getTolerances();
  if (split_groups[0].size() > 0)
    {
      shared_ptr<SurfaceModel> mod(new SurfaceModel(toptol.gap, toptol.gap,
						    toptol.neighbour,
						    toptol.kink, toptol.bend,
						    split_groups[0]));
      surf_mod.push_back(mod);
      inside.push_back(1);
    }

  if (split_groups[1].size() > 0)
    {
      shared_ptr<SurfaceModel> mod(new SurfaceModel(toptol.gap, toptol.gap,
						    toptol.neighbour,
						    toptol.kink, toptol.bend,
						    split_groups[1]));
      surf_mod.push_back(mod);
      inside.push_back(0);
    }

  // Separate surface models into connected sets and make 
  // corresponding volumes
  for (size_t ki=0; ki<surf_mod.size(); ++ki)
    {
      vector<shared_ptr<SurfaceModel> > sep_models;
      sep_models =surf_mod[ki]->getConnectedModels();
      for (size_t kj=0; kj<sep_models.size(); ++kj)
      {
	shared_ptr<ftVolume> curr(new ftVolume(elem_vol->getVolume(), 
					       sep_models[kj]));
	result.push_back(curr);
	is_inside.push_back(inside[ki]);
      }
    }
  return result;
}


//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> >
ftVolumeTools::splitElement(shared_ptr<ParamVolume>& elem_vol, 
			    vector<shared_ptr<ftSurface> >& elem_faces,
			    double* elem_par, ftVolume* trim_vol,
			    double eps, vector<int>& is_inside)
//===========================================================================
{
  vector<shared_ptr<ftVolume> > result;
  is_inside.clear();

  // Fetch trim faces from boundary shells and project trimming curves
  // onto element side surfaces which is coincident with trimming faces
  // NB! Each element side surface is either trimmed by projecting trimming
  // curves from trimming surface or by intersecting with trimming surfaces.
  // This logic will fail if an element side surface is partly coincident with
  // a trimming surface, but intersects with another trimming surface in a
  // different area. The alternative would be to perform both trimming and
  // intersection and remove doubly obtained information, an operation that is
  // risky in itself.
  size_t nmb_elem = elem_faces.size();
  vector<shared_ptr<SurfaceModel> > shells = trim_vol->getAllShells();
  if (shells.size() == 0)
    return result;   // No splitting faces
  vector<shared_ptr<ftSurface> > faces;
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs1(nmb_elem);
  vector<shared_ptr<ParamSurface> > sfs1(nmb_elem);
  vector<bool> at_bd1(nmb_elem, false);
  tpTolerances toptol = shells[0]->getTolerances();
  for (size_t ki=0; ki<shells.size(); ++ki)
    {
      int nmb = shells[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> curr_face = shells[ki]->getFace(kj);
	  if (elem_par != NULL)
	    {
	      // Check if the current face coincides with a constant parameter
	      // of the volume to split
	      shared_ptr<ParamSurface> surf = curr_face->surface();
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
	      size_t kr = 0;
	      for (kr=0; kr<nmb_elem; ++kr)
		{
		  if (dir == (kr/2)+1 && fabs(val-elem_par[kr]) < eps)
		    break;  // Coincidence. Do not include face in intersection
		}
	      if (kr < nmb_elem)
		{
		  // Project trimming curves onto element side surface
		  vector<shared_ptr<CurveOnSurface> > int_cvs;
		  shared_ptr<BoundedSurface> bd_sf;
		  projectTrimCurves(elem_faces[kr], curr_face, eps,
				    toptol.neighbour, int_cvs, bd_sf);
		  if (int_cvs.size() > 0)
		    {
		      all_int_cvs1[kr].insert(all_int_cvs1[kr].end(), int_cvs.begin(),
					      int_cvs.end());
		      sfs1[kr] = bd_sf;
		      at_bd1[kr] = true;
		    }
		  continue;
		}
	    }
	  int bd_stat = ftVolumeTools::boundaryStatus(trim_vol,
						      curr_face, eps);
	  if (bd_stat < 0)
	    faces.push_back(curr_face);
	}
    }

#ifdef DEBUG
  std::ofstream of1("trim_faces.g2");
  for (size_t kr=0; kr<faces.size(); ++kr)
    {
      shared_ptr<ParamSurface> tmp_sf = faces[kr]->surface();
      tmp_sf->writeStandardHeader(of1);
      tmp_sf->write(of1);
    }
  std::ofstream of2("elem_faces.g2");
  for (size_t kr=0; kr<elem_faces.size(); ++kr)
    {
      shared_ptr<ParamSurface> tmp_sf = elem_faces[kr]->surface();
      tmp_sf->writeStandardHeader(of2);
      tmp_sf->write(of2);
    }
#endif

  // Compute intersection curves between remaining element side surfaces and
  // identified trimming surfaces
  size_t nmb_trim = faces.size();
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs2(nmb_trim);
  vector<shared_ptr<ParamSurface> > sfs2(nmb_trim);
  for (size_t kr=0; kr<nmb_elem; ++kr)
    {
      // if (sfs1[kr].get())
      // 	continue;  // Intersection curves already obtained
      
      shared_ptr<ParamSurface> surf1 = elem_faces[kr]->surface();
      BoundingBox box1 = surf1->boundingBox();
     for (size_t kh=0; kh<nmb_trim; ++kh)
	{
	  shared_ptr<ParamSurface> surf2 = faces[kh]->surface();
	  BoundingBox box2 = surf2->boundingBox();

#ifdef DEBUG
	  std::ofstream out("curr_sf_int.g2");
	  surf1->writeStandardHeader(out);
	  surf1->write(out);
	  surf2->writeStandardHeader(out);
	  surf2->write(out);
#endif

	  if (box1.overlaps(box2, eps))
	    {
	      shared_ptr<BoundedSurface> bd1, bd2;
	      vector<shared_ptr<CurveOnSurface> > int_cv1, int_cv2;
	      BoundedUtils::getSurfaceIntersections(sfs1[kr].get() ? sfs1[kr] :
						    surf1, 
						    surf2, eps,
						    int_cv1, bd1,
						    int_cv2, bd2);
	      sfs1[kr] = bd1;
	      sfs2[kh] = bd2;
	      if (int_cv1.size() > 0)
		{
		  all_int_cvs1[kr].insert(all_int_cvs1[kr].end(), 
					  int_cv1.begin(), int_cv1.end());
		  all_int_cvs2[kh].insert(all_int_cvs2[kh].end(), 
					  int_cv2.begin(), int_cv2.end());
		}
	    }
	}
    }

#ifdef DEBUG
  std::ofstream of0("intcurves.g2");
  for (size_t ki=0; ki<nmb_elem; ++ki)
    {
      for (size_t km=0; km<all_int_cvs1[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs1[ki][km]->spaceCurve();
	  tmpcv->writeStandardHeader(of0);
	  tmpcv->write(of0);
	}
    }
  for (size_t ki=0; ki<nmb_trim; ++ki)
    {
      for (size_t km=0; km<all_int_cvs2[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[ki][km]->spaceCurve();
	  tmpcv->writeStandardHeader(of0);
	  tmpcv->write(of0);
	}
    }
  std::ofstream of01("parcurves.g2");
  for (size_t ki=0; ki<nmb_elem; ++ki)
    {
      for (size_t km=0; km<all_int_cvs1[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs1[ki][km]->parameterCurve();
	  tmpcv->writeStandardHeader(of01);
	  tmpcv->write(of01);
	}
    }
  for (size_t ki=0; ki<nmb_trim; ++ki)
    {
      for (size_t km=0; km<all_int_cvs2[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[ki][km]->parameterCurve();
	  tmpcv->writeStandardHeader(of01);
	  tmpcv->write(of01);
	}
    }
#endif

  // Identify intersection/projection curves being coincident with element 
  // edges
  shared_ptr<SplineVolume> spline_vol = 
    dynamic_pointer_cast<SplineVolume,ParamVolume>(elem_vol);
  vector<shared_ptr<SplineCurve> > bd_cvs = spline_vol->getBoundaryCurves();
#ifdef DEBUG
  std::ofstream of3("vol_bd_cvs.g2");
  for (size_t ki=0; ki<bd_cvs.size(); ++ki)
    {
      bd_cvs[ki]->writeStandardHeader(of3);
      bd_cvs[ki]->write(of3);
    }
#endif
  bool modified = false;
  // modified = checkCoincCurves(bd_cvs, all_int_cvs1, all_int_cvs2, 
  // 			      /*0.5**/toptol.neighbour);

  // Check joint between intersection curves
  for (size_t ki=0; ki<all_int_cvs1.size(); ++ki)
    if (all_int_cvs1[ki].size() > 1)
      (void)checkIntCrvJoint(all_int_cvs1[ki], toptol.neighbour, toptol.gap);

  for (size_t ki=0; ki<all_int_cvs2.size(); ++ki)
    if (all_int_cvs2[ki].size() > 1)
      (void)checkIntCrvJoint(all_int_cvs2[ki], toptol.neighbour, toptol.gap);

#ifdef DEBUG
  std::ofstream of0_2("intcurves2.g2");
  for (size_t ki=0; ki<nmb_elem; ++ki)
    {
      for (size_t km=0; km<all_int_cvs1[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs1[ki][km]->spaceCurve();
	  tmpcv->writeStandardHeader(of0_2);
	  tmpcv->write(of0_2);
	}
    }
  for (size_t ki=0; ki<nmb_trim; ++ki)
    {
      for (size_t km=0; km<all_int_cvs2[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[ki][km]->spaceCurve();
	  tmpcv->writeStandardHeader(of0_2);
	  tmpcv->write(of0_2);
	}
    }
  std::ofstream of01_2("parcurves2.g2");
  for (size_t ki=0; ki<nmb_elem; ++ki)
    {
      for (size_t km=0; km<all_int_cvs1[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs1[ki][km]->parameterCurve();
	  tmpcv->writeStandardHeader(of01_2);
	  tmpcv->write(of01_2);
	}
    }
  for (size_t ki=0; ki<nmb_trim; ++ki)
    {
      for (size_t km=0; km<all_int_cvs2[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[ki][km]->parameterCurve();
	  tmpcv->writeStandardHeader(of01_2);
	  tmpcv->write(of01_2);
	}
    }
#endif
  // Make trimmed surfaces and sort trimmed and non-trimmed surfaces according
  // to whether they are inside or outside the trimming shell
  // NB! Must probably be updated when voids are included in the model
  // Sequence: element side surfaces inside trimming shell, element side surfaces
  // outside trimming shell, trimming surfaces inside element, trimming surfaces
  // outside element
  for (size_t kr=0; kr<elem_faces.size(); ++kr)
    if (!sfs1[kr].get())
      sfs1[kr] = elem_faces[kr]->surface();
  for (size_t kr=0; kr<faces.size(); ++kr)
    if (!sfs2[kr].get())
      sfs2[kr] = faces[kr]->surface();

  shared_ptr<SurfaceModel> elem_shell(new SurfaceModel(elem_faces, eps));
  shared_ptr<ftVolume> elem_vol2(new ftVolume(elem_vol, elem_shell, -1));
  vector<vector<shared_ptr<ParamSurface> > > split_groups(4);
  vector<bool> at_bd2(sfs2.size(), false);
  SurfaceModelUtils::sortTrimmedSurfaces(all_int_cvs1, sfs1, at_bd1, elem_vol2.get(), 
					 all_int_cvs2, sfs2, at_bd2, trim_vol, 
					 modified ? /*0.5**/toptol.neighbour : eps, 
					 toptol.bend,
					 split_groups);
					 

  // Remove outdated parameter information
  //for (size_t kr=0; kr<1; ++kr)
  for (size_t kr=2; kr<split_groups.size(); ++kr)
    {
      for (size_t ki=0; ki<split_groups[kr].size(); ++ki)
	{
	  shared_ptr<ParamSurface> surf = split_groups[kr][ki];
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(surf);
	  if (!vol_sf.get())
	    {
	      shared_ptr<BoundedSurface> bd_sf = 
		dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
	      vol_sf =
		dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf->underlyingSurface());
	    }
	  if (vol_sf.get())
	    {
	      // Unset parameter surface information. Not necessarily valid.
	      vol_sf->unsetParamSurf();
	      vol_sf->setVolume(elem_vol);
	    }
	}
    }

   // Add trimming face pieces to the element surfaces
  for (size_t ki=0; ki<split_groups[2].size(); ++ki)
    {
      // Make oppositely oriented copy
      shared_ptr<ParamSurface> surf1 = split_groups[2][ki];
      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();

      // Add surface to element groups
      split_groups[0].push_back(surf1);
      split_groups[1].push_back(surf2);
    }
  
  // Create surface models
  vector<shared_ptr<SurfaceModel> > surf_mod;
  vector<int> inside;
  double eps2 = 1.0e-6;
  double len_tol = 1.5*toptol.neighbour;  // Leave a little slack
  if (split_groups[0].size() > 0)
    {
      // Estimate minimum surface size
      double min_len = HUGE, min_len_all = HUGE;
      vector<size_t> sliver_ix;
      for (size_t ki=0; ki<split_groups[0].size(); ++ki)
	{
	  // Fetch surface category
	  bool elem_face = false;  // Default
	  
	  double len_u, len_v, min_u, max_u, min_v, max_v;
	  split_groups[0][ki]->estimateSfSize(len_u, min_u, max_u, len_v,
					      min_v, max_v);

	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(split_groups[0][ki]);
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(split_groups[0][ki]);
	  
	  if (bd_sf.get())
	    {
	      vol_sf = dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf->underlyingSurface());
	      // Simplify boundary loop if possible
	      // Can use a large angular tolerance since the function checks
	      // on the distance between original and modified curves
	      double max_loop_dist;
	      bool simplified = 
		bd_sf->simplifyBdLoops(toptol.gap, 2.0*toptol.bend, max_loop_dist);
	      int stop_break = 1;
	    }

	  if (vol_sf.get() && vol_sf->atBoundary())
	    elem_face = true;

	  //if (min_u < eps2 || min_v < eps2)
	  if (min_u < toptol.neighbour || min_v < toptol.neighbour)
	    {
	      int nmb_loops = 0, nmb_cvs = 0, nmb_corners;
	      double min_cvlen = len_tol, max_cvlen = 0.0;
	      if (bd_sf.get())
		{
		  bd_sf->getLoopCvInfo(nmb_loops, nmb_cvs, nmb_corners,
				       min_cvlen, max_cvlen, toptol.bend);
		  max_u = std::max(max_u, min_cvlen);
		  max_v = std::max(max_v, min_cvlen);
		}
	      if (/*nmb_cvs <= 4 &&*/ nmb_corners == nmb_cvs &&
		  min_cvlen < len_tol && std::min(len_u, len_v) >= toptol.neighbour)
		min_len_all = std::min(min_len_all, min_cvlen);  
	    }

	  double fac = 2.0;
	  if (std::max(max_u,max_v) < fac*toptol.neighbour && 
	      std::min(max_u,max_v) < toptol.neighbour && elem_face == false)
	    min_len = std::min(min_len, std::min(len_u, len_v));
	  else if (max_u < toptol.neighbour || max_v < toptol.neighbour)
	    sliver_ix.push_back(ki);
	  else
	    min_len = std::min(min_len, std::min(len_u, len_v));
	  min_len_all = std::min(min_len_all, std::min(len_u, len_v));
	}
      double gap = toptol.gap;
      double neighbour = toptol.neighbour;
      if (sliver_ix.size() > 0 && split_groups[0].size() < 6)
	{
	  sliver_ix.clear();
	  min_len = min_len_all;
	}
      if (min_len < neighbour || 
	  (min_len_all < len_tol && sliver_ix.size() == 0))
	{
	  neighbour = (3.0*gap < 0.9*min_len_all) ? 3.0*gap : 0.9*min_len_all;  //0.9*min_len_all;
	  if (gap > 0.75*neighbour)
	    gap = 0.75*neighbour;
	}
      else 
	{
	  // Remove identified sliver faces
	  for (size_t ki=0; ki<sliver_ix.size(); ++ki)
	    split_groups[0].erase(split_groups[0].begin()+sliver_ix[sliver_ix.size()-ki-1]);
	}

      shared_ptr<SurfaceModel> mod;
      bool failed = false;
      try {
	mod = shared_ptr<SurfaceModel>(new SurfaceModel(toptol.gap, gap,
							neighbour,
							toptol.kink, toptol.bend,
							split_groups[0]));
      }
      catch (...)
	{
	  failed = true;
	}
      if (!failed)
	{
	  surf_mod.push_back(mod);
	  inside.push_back(1);
	}
    }

  if (split_groups[1].size() > 0)
    {
      // Estimate minimum surface size
      double min_len = HUGE, min_len_all = HUGE;
      vector<size_t> sliver_ix;
     for (size_t ki=0; ki<split_groups[1].size(); ++ki)
	{
	  // Fetch surface category
	  bool elem_face = false;  // Default
	  
	  double len_u, len_v, min_u, max_u, min_v, max_v;
	  split_groups[1][ki]->estimateSfSize(len_u, min_u, max_u, len_v,
					      min_v, max_v);

	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(split_groups[1][ki]);
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(split_groups[1][ki]);
	  
	  if (bd_sf.get())
	    {
	      vol_sf = dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf->underlyingSurface());
	      // Simplify boundary loop if possible
	      // Can use a large angular tolerance since the function checks
	      // on the distance between original and modified curves
	      double max_loop_dist;
	      bool simplified = 
		bd_sf->simplifyBdLoops(toptol.gap, 2.0*toptol.bend, max_loop_dist);
	      int stop_break = 1;
	    }

	  if (vol_sf.get() && vol_sf->atBoundary())
	    elem_face = true;

	  //if (min_u < eps2 || min_v < eps2)
	  if (min_u < toptol.neighbour || min_v < toptol.neighbour)
	    {
	      int nmb_loops = 0, nmb_cvs = 0, nmb_corners;
	      double min_cvlen = len_tol, max_cvlen = 0.0;
	      if (bd_sf.get())
		{
		  bd_sf->getLoopCvInfo(nmb_loops, nmb_cvs, nmb_corners,
				       min_cvlen, max_cvlen, toptol.bend);
		  max_u = std::max(max_u, min_cvlen);
		  max_v = std::max(max_v, min_cvlen);
		}
	      if (/*nmb_cvs <= 4 &&*/ nmb_corners == nmb_cvs &&
		  min_cvlen < len_tol && std::min(len_u, len_v) >= toptol.neighbour)
		min_len_all = std::min(min_len_all, min_cvlen);  
	    }

	  double fac = 2.0;
	  if (std::max(max_u,max_v) < fac*toptol.neighbour && 
	      std::min(max_u,max_v) < toptol.neighbour && elem_face == false)
	    min_len = std::min(min_len, std::min(len_u, len_v));
	  else if (max_u < toptol.neighbour && max_v < toptol.neighbour)
	    min_len = std::min(min_len, std::min(len_u, len_v));
	  else if (max_u < toptol.neighbour || max_v < toptol.neighbour)
	    sliver_ix.push_back(ki);
	  else
	    min_len = std::min(min_len, std::min(len_u, len_v));
	  min_len_all = std::min(min_len_all, std::min(len_u, len_v));
	}
      double gap = toptol.gap;
      double neighbour = toptol.neighbour;
      if (sliver_ix.size() > 0 && split_groups[1].size() < 6)
	{
	  sliver_ix.clear();
	  min_len = min_len_all;
	}
      if (min_len < neighbour || 
	  (min_len_all < len_tol && sliver_ix.size() == 0))
	{
	  neighbour = (3.0*gap < 0.9*min_len_all) ? 3.0*gap : 0.9*min_len_all;  //0.9*min_len_all;
	  if (gap > 0.75*neighbour)
	    gap = 0.75*neighbour;
	}
      else
	{
	  // Remove identified sliver faces
	  for (size_t ki=0; ki<sliver_ix.size(); ++ki)
	    split_groups[1].erase(split_groups[1].begin()+sliver_ix[sliver_ix.size()-ki-1]);
	}

      shared_ptr<SurfaceModel> mod;
      bool failed = false;
      try {
	mod = shared_ptr<SurfaceModel>(new SurfaceModel(toptol.gap, gap,
							neighbour,
							toptol.kink, toptol.bend,
							split_groups[1]));
      }
      catch (...)
	{
	  failed = true;
	}
      if (!failed)
	{
	  surf_mod.push_back(mod);
	  inside.push_back(0);
	}
    }

  // Separate surface models into connected sets and make 
  // corresponding volumes
  for (size_t ki=0; ki<surf_mod.size(); ++ki)
    {
      vector<shared_ptr<SurfaceModel> > sep_models;
      sep_models =surf_mod[ki]->getConnectedModels();
      for (size_t kj=0; kj<sep_models.size(); ++kj)
      {
	shared_ptr<ftVolume> curr(new ftVolume(elem_vol, sep_models[kj]));
	result.push_back(curr);
	is_inside.push_back(inside[ki]);
      }
    }
  return result;
}


//===========================================================================
// 
// 
bool
ftVolumeTools::updateWithSplitFaces(shared_ptr<SurfaceModel> shell,
				    shared_ptr<ftSurface>& face1,
				    shared_ptr<ftSurface>& face2,
				     vector<pair<ftEdge*, ftEdge*> >& replaced_wires)
//
// This function does not depend on any volume functionality, but is used
// in a volume model context
//===========================================================================
{
  bool modified = false;

  //double eps = shell->getTolerances().gap;
  double eps = std::min(shell->getTolerances().neighbour,
			2.0*shell->getTolerances().gap);
  FaceAdjacency<ftEdgeBase,ftFaceBase> adjacency(shell->getTolerances());

  int nmb = shell->nmbEntities();
  shared_ptr<ParamSurface> surf1 = face1->surface();
  BoundingBox box1 = surf1->boundingBox();
  vector<shared_ptr<CurveOnSurface> > int_cvs_face1;
  vector<ftFaceBase*> newfaces1;
  vector<ftFaceBase*> newfaces2;
  shared_ptr<BoundedSurface> bd1;

#ifdef DEBUG_VOL
  std::ofstream of("split_sfs.g2");
  surf1->writeStandardHeader(of);
  surf1->write(of);
#endif

  int ki;
  for (ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ftSurface> shell_face = shell->getFace(ki);

      // Do not intersect neighbours
      bool smooth;
      if (shell_face->isAdjacent(face1.get(), smooth) ||
	  shell_face->isAdjacent(face2.get(), smooth))
	continue;

      vector<shared_ptr<Vertex> > common_vx1 = 
	shell_face->getCommonVertices(face1.get());
      vector<shared_ptr<Vertex> > common_vx2 = 
	shell_face->getCommonVertices(face2.get());
      if (common_vx1.size() > 0 || common_vx2.size() > 0)
	continue;
      

      shared_ptr<ParamSurface> surf2 = shell->getSurface(ki);
#ifdef DEBUG_VOL
      std::ofstream of2("curr_sf.g2");
      surf2->writeStandardHeader(of2);
      surf2->write(of2);
#endif

      BoundingBox box2 = surf2->boundingBox();
      if (box1.overlaps(box2, eps))
	{
	  shared_ptr<BoundedSurface> bd2;
	  vector<shared_ptr<CurveOnSurface> > int_cv1, int_cv2;
	  try {
	  BoundedUtils::getSurfaceIntersections(surf1, surf2, eps,
						int_cv1, bd1,
						int_cv2, bd2);
	  }
	  catch (...)
	    {
	      int_cv1.clear();
	    }
	  if (int_cv1.size() > 0)
	    {
	      // Make bounded surfaces
	      vector<shared_ptr<BoundedSurface> > trim_sfs;
	      try {
		trim_sfs = 
		  BoundedUtils::splitWithTrimSegments(bd2, int_cv2, eps);
	      }
	      catch(...)
		{
		  std::cout << "Trimmed surfaces missing" << std::endl;
		}
	      if (trim_sfs.size() > 1)
		{
		  modified = true;
		  int_cvs_face1.insert(int_cvs_face1.end(), 
				       int_cv1.begin(), int_cv1.end());

		  // A split is performed. Feth all neighbours
		  vector<ftSurface*> neighbours;
		  shell_face->getAdjacentFaces(neighbours);

		  // Remove split face from the shell
		  shell->removeFace(shell_face);
		  nmb--;
		  ki--;

		  // Compute topology for the new surfaces
		  for (size_t kj=0; kj<trim_sfs.size(); ++kj)
		    {
#ifdef DEBUG_VOL
		      trim_sfs[kj]->writeStandardHeader(of);
		      trim_sfs[kj]->write(of);
#endif

		      vector<ftFaceBase*> neighbours2(neighbours.begin(),
						      neighbours.end());
		      shared_ptr<ftFaceBase> newface(new ftSurface(trim_sfs[kj], -1));
		      adjacency.computeFaceAdjacency(neighbours2, newface.get());

		      // Add split face to the shell. 
		      // The topology pointers are already set
		      shell->append(static_pointer_cast<ftSurface>(newface), false,
				    true);
		      nmb++;

		      // Store face according to configuration related to the
		      // split face
		      // @@@ Temporarily
		      if (kj ==  0)
			newfaces1.push_back(newface.get());
		      else
			newfaces2.push_back(newface.get());
		    }
		}
	    }
	}
    }

  if (int_cvs_face1.size() > 0)
    {
      // Represent face2 as a trimmed surface and prepare new trimming loop
      shared_ptr<ParamSurface> surf2 = face2->surface();
      shared_ptr<BoundedSurface> bd2;
      if (surf2->instanceType() == Class_BoundedSurface)
	bd2 = dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf2);
      else
	{
	  vector<CurveLoop> loops = SurfaceTools::absolutelyAllBoundarySfLoops(surf2,
									       eps);
	  bd2 = shared_ptr<BoundedSurface>(new BoundedSurface(surf2, loops));
	}
      vector<shared_ptr<CurveOnSurface> > int_cvs_face2(int_cvs_face1.size());
      size_t kj;
      for (kj=0; kj<int_cvs_face1.size(); ++kj)
	{
	  int_cvs_face2[kj] = shared_ptr<CurveOnSurface>(int_cvs_face1[kj]->clone());
	  int_cvs_face2[kj]->setUnderlyingSurface(surf2);
	  int_cvs_face2[kj]->unsetParameterCurve();
	}

     // Trim face1
      vector<shared_ptr<BoundedSurface> > surf1_split =
	BoundedUtils::splitWithTrimSegments(bd1, int_cvs_face1, eps);
  
      if (surf1_split.size() != 2)
	{
	  std::cout << "Nmb split sfs1: " << surf1_split.size() << std::endl;
	}

      // Trim face2
      vector<shared_ptr<BoundedSurface> > surf2_split =
	BoundedUtils::splitWithTrimSegments(bd2, int_cvs_face2, eps);

      if (surf2_split.size() != 2)
	{
	  std::cout << "Nmb split sfs2: " << surf2_split.size() << std::endl;
	}
    
      // Replace surfaces only if both faces are split in an inner ant outer part
      if (surf1_split.size() == 2 && surf2_split.size() == 2)
	{
	  // Fetch body and all neighbours to face1
	  Body *body1 = face1->getBody();
	  vector<ftSurface*> neighbours1;
	  face1->getAdjacentFaces(neighbours1);

	  // Fetch boundary edges
	  vector<shared_ptr<ftEdgeBase> > edges1 = face1->getBoundaryLoop(0)->getEdges();

	  // Release face1
	  face1->isolateFace();
      
	  // Fetch body and all neighbours to face2
	  Body *body2 = face2->getBody();
	  vector<ftSurface*> neighbours2;
	  face2->getAdjacentFaces(neighbours2);

	  // Fetch boundary edges
	  vector<shared_ptr<ftEdgeBase> > edges2 = face2->getBoundaryLoop(0)->getEdges();
	  edges1.insert(edges1.end(), edges1.begin(), edges1.end());

	  // Release face2
	  face2->isolateFace();

#ifdef DEBUG_VOL
	  surf1_split[1]->writeStandardHeader(of);
	  surf1_split[1]->write(of);
#endif
	  // Perform topology analysis of the reduced face1
	  shared_ptr<ParamSurface> tmp_sf = surf1_split[1];
	  face1 = shared_ptr<ftSurface>(new ftSurface(tmp_sf, -1));
	  vector<ftFaceBase*> tmp_neigh1(neighbours1.begin(), neighbours1.end());
	  tmp_neigh1.insert(tmp_neigh1.end(), newfaces1.begin(), newfaces1.end());
	  adjacency.computeFaceAdjacency(tmp_neigh1, face1.get());
	  face1->setBody(body1);

#ifdef DEBUG_VOL
	  std::ofstream of3("split_neigh1.g2");
	  for (size_t kr=0; kr<tmp_neigh1.size(); ++kr)
	    {
	      shared_ptr<ParamSurface> tmp_sf = tmp_neigh1[kr]->asFtSurface()->surface();
	      tmp_sf->writeStandardHeader(of3);
	      tmp_sf->write(of3);
	    }
	  face1->surface()->writeStandardHeader(of3);
	  face1->surface()->write(of3);

	  vector<shared_ptr<ftEdge> > edgesy = face1->getAllEdges();
	  for (size_t ix=0; ix<edgesy.size(); ++ix)
	    if (!edgesy[ix]->twin())
	      {
		std::cout << "Missing twin pointer, updateWithSplitFaces " << std::endl;
		of3 << "410 1 0 4 155 55 0 255" << std::endl;
		of3 << "1" << std::endl;
		of3 << edgesy[ix]->point(edgesy[ix]->tMin()) << " " << edgesy[ix]->point(edgesy[ix]->tMax()) << std::endl;
	      }
#endif

	  // Perform topology analysis of the reduced face2
#ifdef DEBUG_VOL
	  surf2_split[1]->writeStandardHeader(of);
	  surf2_split[1]->write(of);
#endif
	  tmp_sf = surf2_split[1];
	  face2 = shared_ptr<ftSurface>(new ftSurface(tmp_sf, -1));
	  vector<ftFaceBase*> tmp_neigh2(neighbours2.begin(), neighbours2.end());
	  tmp_neigh2.insert(tmp_neigh2.end(), newfaces2.begin(), newfaces2.end());
	  adjacency.computeFaceAdjacency(tmp_neigh2, face2.get());
	  face2->setBody(body2);

#ifdef DEBUG_VOL
	  std::ofstream of4("split_neigh2.g2");
	  for (size_t kr=0; kr<tmp_neigh2.size(); ++kr)
	    {
	      shared_ptr<ParamSurface> tmp_sf = tmp_neigh2[kr]->asFtSurface()->surface();
	      tmp_sf->writeStandardHeader(of4);
	      tmp_sf->write(of4);
	    }
	  face2->surface()->writeStandardHeader(of4);
	  face2->surface()->write(of4);

	  vector<shared_ptr<ftEdge> > edgesx = face2->getAllEdges();
	  for (size_t ix=0; ix<edgesx.size(); ++ix)
	    if (!edgesx[ix]->twin())
	      {
		std::cout << "Missing twin pointer, updateWithSplitFaces " << std::endl;
		of4 << "410 1 0 4 155 55 0 255" << std::endl;
		of4 << "1" << std::endl;
		of4 << edgesx[ix]->point(edgesx[ix]->tMin()) << " " << edgesx[ix]->point(edgesx[ix]->tMax()) << std::endl;
	      }
#endif

	  // Check if the replaced_wires information must be updated due to
	  // changes in the edge pointers of face1 and face2
	  // First fetch new outer boundary edges
	  vector<shared_ptr<ftEdgeBase> > edges3 = face1->getBoundaryLoop(0)->getEdges();
	  vector<shared_ptr<ftEdgeBase> > edges4 = face2->getBoundaryLoop(0)->getEdges();
	  edges3.insert(edges3.end(), edges4.begin(), edges4.end());
	  for (kj=0; kj<replaced_wires.size(); ++kj)
	    {
	      size_t kr;
	      for (kr=0; kr<edges1.size(); ++kr)
		if (edges1[kr].get() == replaced_wires[kj].second)
		  break;

	      if (kr < edges1.size())
		{
		  Point pos = 
		    edges1[kr]->point(0.5*(edges1[kr]->tMin()+edges1[kr]->tMax()));
		  int idx = -1;
		  double mindist = 1.0e8;
		  for (int kh=0; kh<(int)edges3.size(); ++kh)
		    {
		      double par, dist;
		      Point pos2;
		      edges3[kh]->closestPoint(pos, par, pos2, dist);
		      if (dist < mindist)
			{
			  idx = kh;
			  mindist = dist;
			}
		    }
		  replaced_wires[kj].second = dynamic_cast<ftEdge*>(edges3[idx].get());
		  continue;
		}
	    }
	    
	  int stop_break = 1;
	}
    }  
  return modified;
}

//===========================================================================
// 
// 
int
ftVolumeTools::boundaryStatus(ftVolume* vol,
			      shared_ptr<ftSurface>& bd_face,
			      double tol)
//===========================================================================
{
  int bd_status = -1;
  shared_ptr<ParamSurface> surf = bd_face->surface();
  shared_ptr<SurfaceOnVolume> vol_sf = 
    dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(surf);
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
  if (bd_sf.get())
    {
      vol_sf = dynamic_pointer_cast<SurfaceOnVolume, ParamSurface>(bd_sf->underlyingSurface());
    }

  if (vol_sf.get())
    {
      int orientation;
      bool swap;
      return vol_sf->whichBoundary(tol, orientation, swap);
    }
  else
    {
      if (!bd_face->hasBody())
	return -1;   // No volume information implies no boundary information available

      if (bd_face->getBody() != vol)
	return -1;   // Inconsistent volume information

      // For each volume boundary, check for coincidence with the given surface
      shared_ptr<ParamVolume> vol2 = vol->getVolume();
      shared_ptr<SurfaceModel> shell = vol->getShell(bd_face.get());
      if (!shell.get())
	return -1;  // Could not find face among the volume boundary faces

      double neighbour = shell->getTolerances().neighbour;

      BoundingBox box1 = surf->boundingBox();
      vector<shared_ptr<ParamSurface> > bd_sfs = vol2->getAllBoundarySurfaces();
      for (size_t ki=0; ki<bd_sfs.size(); ++ki)
	{
	  BoundingBox box2 = bd_sfs[ki]->boundingBox();
	  if (!box1.overlaps(box2, neighbour))
	    continue;  // The surfaces are not coincident
	  
	  // Check coincidence
	  Identity ident;
	  int res = ident.identicalSfs(surf, bd_sfs[ki], tol);
	  if (res == 1 || res == 2)
	    {
	      bd_status = (int)ki;
	      break;
	    }
	}
      if (bd_status >= 0)
	{
	  // The surface belongs to a volume boundary. Store information as a 
	  // SurfaceOnVolume
	  int dir = (bd_status < 2) ? 1 : ((bd_status < 4) ? 2 : 3);
	  Array<double,6> span = vol2->parameterSpan();
	  double par = (bd_status % 2 == 0) ? span[2*dir] : span[2*dir+1];
	  shared_ptr<ParamSurface> surf2(new SurfaceOnVolume(vol2, surf, dir, par,
							     bd_status, false, -1));
	  
	  // Replace current face
	  shared_ptr<ftSurface> face2(new ftSurface(surf2, 0));
	  face2->setBody(vol);
	  ftSurface *twin = bd_face->twin();
	  int ix = shell->getIndex(bd_face);
	  shell->removeFace(bd_face);
	  shell->append(face2, false, false, false, ix);
	  if (twin)
	    {
	      face2->connectTwin(twin, neighbour);
	    }
	}
    }
  return bd_status;
}

//===========================================================================
// 
// 
void
ftVolumeTools::projectTrimCurves(shared_ptr<ftSurface> face1,
				 shared_ptr<ftSurface> face2, 
				 double eps, double eps2,
				 vector<shared_ptr<CurveOnSurface> >& proj_cvs,
				 shared_ptr<BoundedSurface>& bd_sf1)
//===========================================================================
{
  shared_ptr<ParamSurface> surf1 = face1->surface();
  shared_ptr<ParamSurface> surf2 = face2->surface();
  bd_sf1 = BoundedUtils::convertToBoundedSurface(surf1, eps);
  shared_ptr<ParamSurface> under = bd_sf1->underlyingSurface();

  // Fetch all trimming loops
  vector<CurveLoop> loops1 = SurfaceTools::allBoundarySfLoops(surf1, 
							      DEFAULT_SPACE_EPSILON);
  vector<CurveLoop> loops2 = SurfaceTools::allBoundarySfLoops(surf2, 
							      DEFAULT_SPACE_EPSILON);
#ifdef DEBUG
  std::ofstream of("project.g2");
  surf1->writeStandardHeader(of);
  surf1->write(of);
#endif

  // Restrict the boundary/trimming curves of face2 with respect to the
  // domain of face1
  Point pareps1 = SurfaceTools::getParEpsilon(*surf1, eps);
  Point pareps2 = SurfaceTools::getParEpsilon(*surf2, eps);
  double ptol = std::max(eps, std::max(0.5*(pareps1[0]+pareps1[1]),
				       0.5*(pareps2[0]+pareps2[1])));
  for (size_t ki=0; ki<loops2.size(); ++ki)
    {
      int nmb2 = loops2[ki].size();
      for (int kj=0; kj<nmb2; ++kj)
	{
	  shared_ptr<ParamCurve> cv2 = loops2[ki][kj];
#ifdef DEBUG
	  cv2->writeStandardHeader(of);
	  cv2->write(of);
#endif
	  vector<double> param;
	  for (size_t kr=0; kr<loops1.size(); ++kr)
	    {
	      int nmb1 = loops1[kr].size();
	      for (int kh=0; kh<nmb1; ++kh)
		{
		  shared_ptr<ParamCurve> cv1 = loops1[kr][kh];

		  // Intersect curves using an increased tolerance to
		  // avoid curve pieces
		  vector<pair<double,double> > int_pts;
		  intersectParamCurves(cv1.get(), cv2.get(), eps2, int_pts);

		  for (size_t ka=0; ka<int_pts.size(); ++ka)
		    {
		      double par1, par2, dist;
		      Point ptc1, ptc2;
		      ClosestPoint::closestPtCurves(cv1.get(), cv2.get(), 
						    cv1->startparam(), cv1->endparam(),
						    cv2->startparam(), cv2->endparam(),
						    int_pts[ka].first, int_pts[ka].second,
						    par1, par2, dist, ptc1, ptc2);
		      param.push_back(par2);
		    }
		}
	    }
	  if (param.size() > 0)
	    {
	      std::sort(param.begin(), param.end());
	      for (size_t kr=1; kr<param.size();)
		{
		  if (param[kr]-param[kr-1] < ptol)
		    param.erase(param.begin()+kr);
		  else
		    ++kr;
		}
	      if (cv2->startparam() < param[0]-ptol)
		param.insert(param.begin(), cv2->startparam());
	      if (cv2->endparam() > param[param.size()-1]+ptol)
		param.push_back(cv2->endparam());
	    }
	  else
	    {
	      param.push_back(cv2->startparam());
	      param.push_back(cv2->endparam());
	    }

	  // For each curve segment, check if it is inside face1
	  vector<shared_ptr<ParamCurve> > sub_cv;
	  for (size_t kr=1; kr<param.size(); ++kr)
	    {
	      Point pos = cv2->point(0.5*(param[kr-1]+param[kr]));
	      double upar, vpar, sfdist;
	      Point sfpt;
	      surf1->closestPoint(pos, upar, vpar, sfpt, sfdist, eps);
	      if (sfdist < eps2)
		{
		  shared_ptr<ParamCurve> sub_cv(cv2->subCurve(param[kr-1], param[kr]));
		  shared_ptr<Point> dummy_start, dummy_end;
		  shared_ptr<SplineCurve> par_cv(CurveCreators::projectSpaceCurve(sub_cv, 
										  surf1, 
										  dummy_start,
										  dummy_end,
										  eps));
		  if (sub_cv->estimatedCurveLength() > eps2 &&
		      par_cv->estimatedCurveLength() > ptol)
		    {
		      shared_ptr<CurveOnSurface> sf_cv(new CurveOnSurface(under, par_cv, true));
		      (void)sf_cv->ensureSpaceCrvExistence(eps);
		      proj_cvs.push_back(sf_cv);
		    }
		}
	    }
	}
    }
}

//===========================================================================
// 
// 
bool
ftVolumeTools::checkCoincCurves(vector<shared_ptr<SplineCurve> >& bd_cvs,
				vector<vector<shared_ptr<CurveOnSurface> > >& int_cvs1,
				vector<vector<shared_ptr<CurveOnSurface> > >& int_cvs2,
				double tol)
//===========================================================================
{
#ifdef DEBUG
  std::ofstream of1("coinc1.g2");
  std::ofstream of2("coinc2.g2");
#endif
  bool modified = false;
  double ptol = 1.0e-8;
  Identity ident;
  for (size_t ki=0; ki<bd_cvs.size(); ++ki)
    {
      double t1 = bd_cvs[ki]->startparam();
      double t2 = bd_cvs[ki]->endparam();
      for (size_t kj=0; kj<int_cvs1.size(); ++kj)
	for (size_t kr=0; kr<int_cvs1[kj].size();)
	  {
	    // Check for identical end points
	    double t3 = int_cvs1[kj][kr]->startparam();
	    double t4 = int_cvs1[kj][kr]->endparam();
	    Point pos3 = int_cvs1[kj][kr]->ParamCurve::point(t3);
	    Point pos4 = int_cvs1[kj][kr]->ParamCurve::point(t4);

	    double ct1, ct2, cd1, cd2;
	    Point cpt1, cpt2;
	    bd_cvs[ki]->closestPoint(pos3, t1, t2, ct1, cpt1, cd1);
	    bd_cvs[ki]->closestPoint(pos4, t1, t2, ct2, cpt2, cd2);
	    double start, end;
	    if (cd1 < tol && cd2 < tol && fabs(ct2-ct1) > ptol)
	      {
		start = std::min(ct1, ct2);
		end = std::max(ct1, ct2);
	      }
	    else
	      {
		++kr;
		continue;
	      }
	    
	    int coinc = ident.identicalCvs(bd_cvs[ki], start, end,
					   int_cvs1[kj][kr],
					   t3, t4, tol);
	  if (coinc > 0)
	    {
#ifdef DEBUG
	      shared_ptr<ParamCurve> tmpcv = int_cvs1[kj][kr]->spaceCurve();
	      tmpcv->writeStandardHeader(of1);
	      tmpcv->write(of1);
#endif

	      // Remove or replace
	      modified = true;
	      shared_ptr<ParamCurve> sub_bd(bd_cvs[ki]->subCurve(start, end));
	      double len = sub_bd->estimatedCurveLength();
	      if (len < tol)
		{
		  int_cvs1[kj].erase(int_cvs1[kj].begin()+kr);
		  std::cout << "ftVolumeTools::checkCoincCurves, remove - len = " << len;
		  std::cout << ", dist = (" << cd1 << "," << cd2 << ")" << std::endl;
		  //++kr;
		}
	      else //if (false)
		{
		  //modified = true;
		  // if (ct1 > ct2)
		  //   sub_bd->reverseParameterDirection();
		  // shared_ptr<CurveOnSurface> tmp_cv(new CurveOnSurface(int_cvs1[kj][kr]->underlyingSurface(),
		  // 						       sub_bd,
		  // 						       false));
		  // tmp_cv->ensureParCrvExistence(tol);
		  // int_cvs1[kj][kr] = tmp_cv;
		  // ++kr;
		  int_cvs1[kj].erase(int_cvs1[kj].begin()+kr);
		  std::cout << "ftVolumeTools::checkCoincCurves, remove2 - len = " << len;
		  std::cout << ", dist = (" << cd1 << "," << cd2 << ")" << std::endl;
		}
	      // else
	      // 	++kr;
	    }
	  else
	    ++kr;
	  }
 	      
      for (size_t kj=0; kj<int_cvs2.size(); ++kj)
	for (size_t kr=0; kr<int_cvs2[kj].size();)
	  {
	    // Check for identical end points
	    double t3 = int_cvs2[kj][kr]->startparam();
	    double t4 = int_cvs2[kj][kr]->endparam();
	    Point pos3 = int_cvs2[kj][kr]->ParamCurve::point(t3);
	    Point pos4 = int_cvs2[kj][kr]->ParamCurve::point(t4);

	    double ct1, ct2, cd1, cd2;
	    Point cpt1, cpt2;
	    bd_cvs[ki]->closestPoint(pos3, t1, t2, ct1, cpt1, cd1);
	    bd_cvs[ki]->closestPoint(pos4, t1, t2, ct2, cpt2, cd2);
	    double start, end;
	    if (cd1 < tol && cd2 < tol && fabs(ct2-ct1) > ptol)
	      {
		start = std::min(ct1, ct2);
		end = std::max(ct1, ct2);
	      }
	    else
	      {
		++kr;
		continue;
	      }
	    
	    int coinc = ident.identicalCvs(bd_cvs[ki], start, end,
					   int_cvs2[kj][kr],
					   t3, t4, tol);
	  if (coinc > 0)
	    {
#ifdef DEBUG
	      shared_ptr<ParamCurve> tmpcv = int_cvs2[kj][kr]->spaceCurve();
	      tmpcv->writeStandardHeader(of2);
	      tmpcv->write(of2);
#endif
	      // Remove or replace
	      modified = true;
	      shared_ptr<ParamCurve> sub_bd(bd_cvs[ki]->subCurve(start, end));
	      double len = sub_bd->estimatedCurveLength();
	      if (len < tol)
		{
		  int_cvs2[kj].erase(int_cvs2[kj].begin()+kr);
		  std::cout << "ftVolumeTools::checkCoincCurves, remove - len = " << len;
		  std::cout << ", dist = (" << cd1 << "," << cd2 << ")" << std::endl;
		  //++kr;
		}
	      else //if (false)
		{
		  //modified = true;
		  if (ct1 > ct2)
		    sub_bd->reverseParameterDirection();
		  shared_ptr<CurveOnSurface> tmp_cv(new CurveOnSurface(int_cvs2[kj][kr]->underlyingSurface(),
		  						       sub_bd,
		  						       false));
		  tmp_cv->ensureParCrvExistence(tol);
		  int_cvs2[kj][kr] = tmp_cv;
		  ++kr;
		  std::cout << "ftVolumeTools::checkCoincCurves, move - len = " << len;
		  std::cout << ", dist = (" << cd1 << "," << cd2 << ")" << std::endl;
		  // int_cvs2[kj].erase(int_cvs2[kj].begin()+kr);
		}
	      // else 
	      // 	++kr;
	    }
	  else
	    ++kr;
 	}
    }
  return modified;
}

struct IntcrvJoint
{
  double mind, minp;
  Point clopt;
  int idx;

  IntcrvJoint(int ix, double dist, double par, Point close)
  {
    idx = ix;
    mind = dist;
    minp = par;
    clopt = close;
  }
};

//===========================================================================
// 
// 
bool
ftVolumeTools::checkIntCrvJoint(vector<shared_ptr<CurveOnSurface> > & int_cvs,
				double tol, double eps)
//===========================================================================
{
  bool modified = false;
 
  // TEST
  eps = 1.0e-6;

  size_t nmb = int_cvs.size();
#ifdef DEBUG
  std::ofstream of("curr_intcurves.g2");
  for (size_t ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamCurve> tmpcv = int_cvs[ki]->spaceCurve();
      tmpcv->writeStandardHeader(of);
      tmpcv->write(of);
    }
#endif

  double min_cv_len = HUGE;
  double max_cv_len = 0.0;
  for (size_t ki=0; ki<nmb; ++ki)
    {
      double len = int_cvs[ki]->estimatedCurveLength();
      min_cv_len = std::min(len, min_cv_len);
      max_cv_len = std::max(len, max_cv_len);
    }
  
  // Start with checking for totally coincident curves
  double len_fac = 10.0;
  if (max_cv_len > len_fac*std::min(tol,min_cv_len) &&
      max_cv_len > min_cv_len + tol)
    {
      // Careful with removing curves in a small configuration
      double tol2 = std::min(tol, 0.9*min_cv_len);
      for (size_t ki=0; ki<nmb;)
	{
	  double t1 = int_cvs[ki]->startparam();
	  double t2 = int_cvs[ki]->endparam();
	  Point pos1 = int_cvs[ki]->ParamCurve::point(t1);
	  Point pos2 = int_cvs[ki]->ParamCurve::point(t2);
	  size_t kj;
	  for (kj=ki+1; kj<nmb; ++kj)
	    {
	      double t3 = int_cvs[kj]->startparam();
	      double t4 = int_cvs[kj]->endparam();
	      Point pos3 = int_cvs[kj]->ParamCurve::point(t3);
	      Point pos4 = int_cvs[kj]->ParamCurve::point(t4);
	  
	      double d1, d2, d3, d4, par1, par2, par3, par4;
	      Point clo1, clo2, clo3, clo4;
	      int_cvs[kj]->closestPoint(pos1, t3, t4, par1, clo1, d1);
	      int_cvs[kj]->closestPoint(pos2, t3, t4, par2, clo2, d2);
	      int_cvs[ki]->closestPoint(pos3, t1, t2, par3, clo3, d3);
	      int_cvs[ki]->closestPoint(pos4, t1, t2, par4, clo4, d4);

	      if (d1 < tol2 && d2 < tol2 && d1+d2 < d3+d4)
		{
		  // Check if curve ki is totally coincident with curve kj
		  Identity ident;
		  int coinc = ident.identicalCvs(int_cvs[ki], t1, t2, int_cvs[kj],
						 par1, par2, tol);
		  if (coinc == 1)
		    {
		      // Remove curve
		      int_cvs.erase(int_cvs.begin()+ki);
		      nmb--;
		      break;
		    }
		}
	      else if (d3 < tol2 && d4 < tol2)
		{
		  // Check if curve kj is totally coincident with curve ki
		  Identity ident;
		  int coinc = ident.identicalCvs(int_cvs[ki], par3, par4, int_cvs[kj],
						 par1, par2, tol2);
		  if (coinc == 1)
		    {
		      // Remove curve
		      int_cvs.erase(int_cvs.begin()+kj);
		      nmb--;
		      break;
		    }
		}
	    }
	  if (kj == nmb)
	    ++ki;
	}
    }

#ifdef DEBUG
  std::ofstream of2("curr_intcurves2.g2");
  for (size_t ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamCurve> tmpcv = int_cvs[ki]->spaceCurve();
      tmpcv->writeStandardHeader(of2);
      tmpcv->write(of2);
    }
#endif

  vector<int> remove_cvs;
  double lim_ang = 0.01; //0.25*M_PI;
  for (size_t ki=0; ki<nmb; ++ki)
    {
      double t1 = int_cvs[ki]->startparam();
      double t2 = int_cvs[ki]->endparam();
      Point pos1 = int_cvs[ki]->ParamCurve::point(t1);
      Point pos2 = int_cvs[ki]->ParamCurve::point(t2);

      // Compute the closest points between endpoints of this curve and the other
      // curves
      int clo_ix1 = -1, clo_ix2 = -1;
      double mind1 = HUGE, mind2 = HUGE;
      double minp1, minp2;
      Point clopt1, clopt2;
      vector<IntcrvJoint> joint1;
      vector<IntcrvJoint> joint2;
      for (size_t kj=0; kj<nmb; ++kj)
	{
	  if (ki == kj)
	    continue;

	  double d1, d2, par1, par2;
	  Point clo1, clo2;
	  double t3 = int_cvs[kj]->startparam();
	  double t4 = int_cvs[kj]->endparam();
	  int_cvs[kj]->closestPoint(pos1, t3, t4, par1, clo1, d1);
	  int_cvs[kj]->closestPoint(pos2, t3, t4, par2, clo2, d2);

	  if (d1 < tol)
	    {
	      // Check angle
	      vector<Point> der1(2), der2(2);
	      int_cvs[ki]->point(der1, t1, 1);
	      int_cvs[kj]->point(der2, par1, 1);
	      double ang = der1[1].angle(der2[1]);
	      ang = std::min(ang, M_PI-ang);
	      if (ang > lim_ang)
		{
		  IntcrvJoint tmp((int)kj, d1, par1, clo1);
		  joint1.push_back(tmp);
		}
	    }

	  if (d2 < tol)
	    {
	      // Check angle
	      vector<Point> der1(2), der2(2);
	      int_cvs[ki]->point(der1, t2, 1);
	      int_cvs[kj]->point(der2, par2, 1);
	      double ang = der1[1].angle(der2[1]);
	      ang = std::min(ang, M_PI-ang);
	      if (ang > lim_ang)
		{
		  IntcrvJoint tmp((int)kj, d2, par2, clo2);
		  joint2.push_back(tmp);
		}
	    }
	  // if (d1 < mind1)
	  //   {
	  //     mind1 = d1;
	  //     clo_ix1 = (int)kj;
	  //     minp1 = par1;
	  //     clopt1 = clo1;
	  //   }

	  // if (d2 < mind2)
	  //   {
	  //     mind2 = d2;
	  //     clo_ix2 = (int)kj;
	  //     minp2 = par2;
	  //     clopt2 = clo2;
	  //   }
	}
    

      if (false /*mind1 < tol && mind2 < tol && clo_ix1 == clo_ix2*/)
	{
	  // Check for coincidence
	  Identity ident;
	  double start = std::min(minp1,minp2);
	  double end = std::max(minp1,minp2);
	  int coinc = ident.identicalCvs(int_cvs[ki], t1, t2, int_cvs[clo_ix1],
					 start, end, tol);

	  if (coinc == 1)
	    {
	      // Check if the entire curve is coincident
	      double t3 = int_cvs[clo_ix1]->startparam();
	      double t4 = int_cvs[clo_ix1]->endparam();
	      Point pos3 = int_cvs[clo_ix1]->ParamCurve::point(t3);
	      Point pos4 = int_cvs[clo_ix1]->ParamCurve::point(t4);
	      if (pos3.dist((minp1 < minp2) ? clopt1 : clopt2) < tol &&
		  pos4.dist((minp1 < minp2) ? clopt2 : clopt1) < tol )
		{
		  // Remove the curve with the poorest connection to its neighbours
		  for (size_t kj=0; kj<nmb; ++kj)
		    {
		      if (kj == ki || kj == clo_ix1)
			continue;
		      Point pos5 = 
			int_cvs[kj]->ParamCurve::point(int_cvs[kj]->startparam());
		      Point pos6 = 
			int_cvs[kj]->ParamCurve::point(int_cvs[kj]->endparam());
		      double d1 = std::min(pos1.dist(pos5), pos1.dist(pos6));
		      double d2 = std::min(pos2.dist(pos5), pos2.dist(pos6));
		      double d3 = std::min(pos3.dist(pos5), pos3.dist(pos6));
		      double d4 = std::min(pos4.dist(pos5), pos4.dist(pos6));
		      double dd1 = std::min(d1, d2);
		      double dd2 = std::min(d3, d4);
		      // Too simple test. Should compare distances in both
		      // endpoints simultanously
		      if (dd1 < tol && dd2 < tol && dd1 < dd2)
			{
			  // Remove curve clo_ix1
			  int_cvs.erase(int_cvs.begin()+clo_ix1);
			  nmb--;
			  break;
			}
		      else if (dd1 < tol && dd2 < tol && dd2 <= dd1)
			{
			  // Remove curve ki
			  std::swap(int_cvs[ki], int_cvs[clo_ix1]);
			  int_cvs.erase(int_cvs.begin()+clo_ix1);
			  nmb--;
			  break;
			}
		    }
	      
		}
	      else
		{
		  // Remove the smallest curve
		  std::swap(int_cvs[ki], int_cvs[clo_ix1]);
		  int_cvs.erase(int_cvs.begin()+clo_ix1);
		  nmb--;
		}

	      // Set parameters for no further modifications of these curves
	      mind1 = 2.0*tol;
	      mind2 = 2.0*tol;
	    }
	  int stop_break = 1;
	}
	
      // Clean joints
      for (int ka1=0; ka1<joint1.size(); ++ka1)
	{
	  for (int ka2=0; ka2<joint2.size(); ++ka2)
	    {
	      if (joint1[ka1].idx == joint2[ka2].idx)
		{
		  if (joint1[ka1].mind < joint2[ka2].mind)
		    {
		      joint2.erase(joint2.begin()+ka2);
		      ka2--;
		    }
		  else if (joint1[ka1].mind > joint2[ka2].mind)
		    {
		      joint1.erase(joint1.begin()+ka1);
		      ka1--;
		      break;
		    }
		}
	    }
	}

      for (size_t kj=0; kj<joint1.size(); ++kj)
	{
	  // Check if the closest point corresponds to an endpoint of the curve
	  int ix = joint1[kj].idx;
	  double t2_1 = int_cvs[ix]->startparam();
	  double t2_2 = int_cvs[ix]->endparam();
	  Point pt1 = int_cvs[ix]->ParamCurve::point(t2_1);
	  Point pt2 = int_cvs[ix]->ParamCurve::point(t2_2);
	  Point pt3 = int_cvs[ix]->ParamCurve::point(joint1[kj].minp);
	  double dist1 = pt1.dist(pt3);
	  double dist2 = pt2.dist(pt3);
	  double tol2 = eps; //std::max(0.5*mind1, eps);
	  if (dist1 > tol2 && dist2 > tol2)
	    {
	      // Split curve
	      modified = true;
	      int ix1  = (dist1 < dist2) ? 1 : 0;
	      for (size_t kr=0; kr<joint2.size(); ++kr)
	      	{
		  if (joint2[kr].idx != ix)
		    continue;

	      	  // Check coincidene
	      	  Point mid1 = 
	      	    int_cvs[ix]->ParamCurve::point(0.5*(t2_1+joint1[kj].minp));
	      	  Point mid2 = 
	      	    int_cvs[ix]->ParamCurve::point(0.5*(t2_2+joint1[kj].minp));

	      	  double tp1, tp2, td1, td2;
	      	  Point tcl1, tcl2;
	      	  int_cvs[ki]->closestPoint(mid1, t1, t2, tp1, tcl1, td1);
	      	  int_cvs[ki]->closestPoint(mid2, t1, t2, tp2, tcl2, td2);
	      	  if (td1 > tol && td2 <= tol) 
	      	    ix1 = 0;
	      	  else if (td2 > tol && td1 <= tol)
	      	    ix1 = 1;
	      	  else if (joint2[kr].minp < joint1[kj].minp)
	      	    ix1 = 0;
	      	  else
	      	    ix1 = 1;
	      	}

	      int ix2 = 1 - ix1;
	      double dist = (ix1 == 0) ? dist2 : dist1;
	      vector<shared_ptr<ParamCurve> > sub_cvs = 
		int_cvs[ix]->split(joint1[kj].minp);
	      int_cvs[ix] = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[ix1]);
	      // if (dist > tol)
	      // 	{
	      // 	  int_cvs.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[ix2]));
	      // 	}

	      // Check for coincidence with dismissed sub curve
	      for (size_t kr=0; kr<joint1.size(); ++kr)
		{
		  if (kr == kj)
		    continue;

		  int idx2 = joint1[kr].idx;
		  double len1 = int_cvs[idx2]->estimatedCurveLength();
		  double len2 = sub_cvs[ix2]->estimatedCurveLength();
		  int coinc = 0;
		  if (fabs(len1-len2) < std::min(tol, std::min(len1,len2)))
		    {
		      Identity ident;
		      coinc = ident.identicalCvs(sub_cvs[ix2], 
						 sub_cvs[ix2]->startparam(),
						 sub_cvs[ix2]->endparam(), 
						 int_cvs[idx2],
						 int_cvs[idx2]->startparam(),
						 int_cvs[idx2]->endparam(), 
						 tol);
		    }
		  if (coinc == 1)
		    {
		      // Mark curve for removal
		      remove_cvs.push_back(idx2);
		    }
		}
	    }
	}

      for (size_t kj=0; kj<joint2.size(); ++kj)
	{
	  int ix = joint2[kj].idx;
	  if (int_cvs[ix]->startparam() >= joint2[kj].minp ||
	      int_cvs[ix]->endparam() <= joint2[kj].minp)
	    continue;

	  // Check if the closest point corresponds to an endpoint of the curve
	  Point pt1 = int_cvs[ix]->ParamCurve::point(int_cvs[ix]->startparam());
	  Point pt2 = int_cvs[ix]->ParamCurve::point(int_cvs[ix]->endparam());
	  Point pt3 = int_cvs[ix]->ParamCurve::point(joint2[kj].minp);
	  double dist1 = pt1.dist(pt3);
	  double dist2 = pt2.dist(pt3);
	  double tol2 = eps; //std::max(0.5*mind2, eps);
	  if (dist1 > tol2 && dist2 > tol2)
	    {
	      // Split curve
	      modified = true;
	      int ix1 = (dist1 < dist2) ? 1 : 0;
	      int ix2 = 1 - ix1;
	      double dist = std::min(dist1, dist2);
	      vector<shared_ptr<ParamCurve> > sub_cvs = 
		int_cvs[ix]->split(joint2[kj].minp);
	      int_cvs[ix] = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[ix1]);
	      // if (dist > tol)
	      // 	{
	      // 	  int_cvs.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[ix2]));
	      // 	}

	      // Check for coincidence with dismissed sub curve
	      for (size_t kr=0; kr<joint2.size(); ++kr)
		{
		  if (kr == kj)
		    continue;

		  int idx2 = joint2[kr].idx;
		  double len1 = int_cvs[idx2]->estimatedCurveLength();
		  double len2 = sub_cvs[ix2]->estimatedCurveLength();
		  int coinc = 0;
		  if (fabs(len1-len2) < std::min(tol, std::min(len1,len2)))
		    {
		      Identity ident;
		      ident.identicalCvs(sub_cvs[ix2], 
					 sub_cvs[ix2]->startparam(),
					 sub_cvs[ix2]->endparam(), 
					 int_cvs[idx2],
					 int_cvs[idx2]->startparam(),
					 int_cvs[idx2]->endparam(), tol);
		    }
		  if (coinc == 1)
		    {
		      // Mark curve for removal
		      remove_cvs.push_back(idx2);
		    }
		}
	    }
	}
    }
  if (remove_cvs.size() > 0)
    {
      std::sort(remove_cvs.begin(), remove_cvs.end());
      for (int ka=(int)remove_cvs.size()-1; ka>=0; --ka)
	int_cvs.erase(int_cvs.begin()+remove_cvs[ka]);
    }

  return modified;
}
