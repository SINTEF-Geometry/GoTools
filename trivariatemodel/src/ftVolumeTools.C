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
ftVolumeTools::splitVolumes(ftVolume* vol, 
			    shared_ptr<ParamSurface>& surface, double eps,
			    int create_all)
//===========================================================================
{
  shared_ptr<ftSurface> face(new ftSurface(surface, -1));
  return splitVolumes(vol, face, eps, create_all);
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
  return splitVolumes(vol.get(), face, eps, create_all);
}

//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> >
ftVolumeTools::splitVolumes(ftVolume* vol, 
			    shared_ptr<ftSurface>& face, double eps,
			    int create_all)
//===========================================================================
{
#ifdef DEBUG_VOL
  std::cout << "Entering splitVolumes" << std::endl;
#endif
  vector<shared_ptr<ftVolume> > result;

  // Fetch all boundary surfaces
  vector<shared_ptr<SurfaceModel> > shells = vol->getAllShells();

  if (shells.size() == 0)
    return result;

  double tol = shells[0]->getTolerances().gap;

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

  // Check for identical faces between the surfaces of the initial face
  // internal to the initial volume that are coincident with a surface in any
  // of the new volume pieces. These surfaces must be removed
  removeCoincFaces(split_models[0], split_models[1], split_models[2], tol);

  // For each surface from the initial face internal to the initial volume: make a copy, swap
  // the parameter directions of this copy, and add it to the outside
  // part of the other volume.
  // Set also twin pointers between corresponding boundary surfaces
#ifdef DEBUG_VOL
  if (split_models[3].get())
    std::cout << "Model 3(2): " << split_models[3]->nmbEntities() << std::endl;
#endif
  closeModelParts(split_models[0], split_models[1], vol, split_models[2], 
		  1, eps, create_all);

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
  if (split_models[0] && (create_all == 1 || create_all >= 3))
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
  
  if (split_models[1] && (create_all == 2 || create_all >= 3))
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
		  // Include eventual boundary condition information
		  // in the corresponding element side face
		  if (curr_face->hasBoundaryConditions())
		    {
		      int bd_type, bd;
		      curr_face->getBoundaryConditions(bd_type, bd);
		      elem_faces[kr]->setBoundaryConditions(bd_type, bd);
		    }

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
  std::ofstream ofproj("projcurves.g2");
  for (size_t ki=0; ki<nmb_elem; ++ki)
    {
      for (size_t km=0; km<all_int_cvs1[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs1[ki][km]->spaceCurve();
	  tmpcv->writeStandardHeader(ofproj);
	  tmpcv->write(ofproj);
	}
    }
#endif

  // Compute intersection curves between remaining element side surfaces and
  // identified trimming surfaces
  size_t nmb_trim = faces.size();
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs2(nmb_trim);
  vector<shared_ptr<ParamSurface> > sfs2(nmb_trim);
  for (size_t kr=0; kr<nmb_elem; ++kr)
    {
      int nmb_project = (int)all_int_cvs1[kr].size();
      
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
	      if (int_cv1.size() > 0 && nmb_project > 0)
		{
		  // Remove intersection curves that are coincident with
		  // already found projection curves
		  checkIntCvCoincidence(&all_int_cvs1[kr][0], nmb_project,
					toptol.neighbour, toptol.gap,
					int_cv1, int_cv2);
		}

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

  // Intersection curves found by projection may relate to different 
  // bounded surfaces. Ensure consistence
  for (size_t kr=0; kr<all_int_cvs1.size(); ++kr)
    {
      shared_ptr<BoundedSurface> bd_sf =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs1[kr]);
      if (bd_sf.get())
	{
	  shared_ptr<ParamSurface> tmp_sf = bd_sf->underlyingSurface();
	  for (size_t kh=0; kh<all_int_cvs1[kr].size(); ++kh)
	    {
	      if (all_int_cvs1[kr][kh]->underlyingSurface().get() != tmp_sf.get())
		all_int_cvs1[kr][kh]->setUnderlyingSurface(tmp_sf);
	    }
	}
    }

  for (size_t kr=0; kr<all_int_cvs2.size(); ++kr)
    {
      shared_ptr<BoundedSurface> bd_sf =
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs2[kr]);
      if (bd_sf.get())
	{
	  shared_ptr<ParamSurface> tmp_sf = bd_sf->underlyingSurface();
	  for (size_t kh=0; kh<all_int_cvs2[kr].size(); ++kh)
	    {
	      if (all_int_cvs2[kr][kh]->underlyingSurface().get() != tmp_sf.get())
		all_int_cvs2[kr][kh]->setUnderlyingSurface(tmp_sf);
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
	  if (tmpcv.get())
	    {
	      tmpcv->writeStandardHeader(of01);
	      tmpcv->write(of01);
	    }
	}
    }
  for (size_t ki=0; ki<nmb_trim; ++ki)
    {
      for (size_t km=0; km<all_int_cvs2[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[ki][km]->parameterCurve();
	  if (tmpcv.get())
	    {
	      tmpcv->writeStandardHeader(of01);
	      tmpcv->write(of01);
	    }
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
      (void)checkIntCrvJoint(all_int_cvs1[ki], toptol.neighbour, toptol.gap,
			     toptol.bend);

  for (size_t ki=0; ki<all_int_cvs2.size(); ++ki)
    if (all_int_cvs2[ki].size() > 1)
      (void)checkIntCrvJoint(all_int_cvs2[ki], toptol.neighbour, toptol.gap,
			     toptol.bend);

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
	  if (tmpcv.get())
	    {
	      tmpcv->writeStandardHeader(of01_2);
	      tmpcv->write(of01_2);
	    }
	}
    }
  for (size_t ki=0; ki<nmb_trim; ++ki)
    {
      for (size_t km=0; km<all_int_cvs2[ki].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[ki][km]->parameterCurve();
	  if (tmpcv.get())
	    {
	      tmpcv->writeStandardHeader(of01_2);
	      tmpcv->write(of01_2);
	    }
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

  // // Simplify boundary loops if possible
  // for (size_t ki=0; ki<sfs1.size(); ++ki)
  //   {
  //     shared_ptr<BoundedSurface> bd_sf = 
  // 	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs1[ki]);
  //     if (bd_sf.get())
  // 	{
  // 	  double max_loop_dist;
  // 	  bool simplified = 
  // 	    bd_sf->simplifyBdLoops(toptol.gap, toptol.bend, max_loop_dist);
  // 	  int stop_break = 1;
  // 	}
  //   }
  // for (size_t ki=0; ki<sfs2.size(); ++ki)
  //   {
  //     shared_ptr<BoundedSurface> bd_sf = 
  // 	dynamic_pointer_cast<BoundedSurface,ParamSurface>(sfs2[ki]);
  //     if (bd_sf.get())
  // 	{
  // 	  double max_loop_dist;
  // 	  bool simplified = 
  // 	    bd_sf->simplifyBdLoops(toptol.gap, toptol.bend, max_loop_dist);
  // 	  int stop_break = 1;
  // 	}
  //   }

  // Number of element side surfaces
  int nmb_split1 = sfs1.size();

  shared_ptr<SurfaceModel> elem_shell(new SurfaceModel(elem_faces, eps));
  shared_ptr<ftVolume> elem_vol2(new ftVolume(elem_vol, elem_shell, -1));
  vector<vector<pair<shared_ptr<ParamSurface>, int> > > split_groups(4);
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
	  shared_ptr<ParamSurface> surf = split_groups[kr][ki].first;
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
  double len_tol = 1.5*toptol.neighbour;  // Leave a little slack
  for (size_t ki=0; ki<split_groups[2].size(); ++ki)
    {
      // Make oppositely oriented copy
      int id = split_groups[2][ki].second;
      shared_ptr<ParamSurface> surf1 = split_groups[2][ki].first;
      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();

      // Add surface to element groups, small ones first to improve
      // stability of topology analysis
      double len_u, len_v, min_u, max_u, min_v, max_v;
      surf1->estimateSfSize(len_u, min_u, max_u, len_v, min_v, max_v);
      if (len_u < len_tol || len_v < len_tol)
	{
	  split_groups[0].insert(split_groups[0].begin(),
				 make_pair(surf1,nmb_split1+id));
	  split_groups[1].insert(split_groups[1].begin(),
				 make_pair(surf2, nmb_split1+id));
	}
      else
	{
	  split_groups[0].push_back(make_pair(surf1, nmb_split1+id));
	  split_groups[1].push_back(make_pair(surf2, nmb_split1+id));
	}
    }
  
  // Create surface models
  vector<shared_ptr<SurfaceModel> > surf_mod;
  vector<int> inside;
  double eps2 = 1.0e-6;
  double small_fac = 0.8;
  if (split_groups[0].size() > 0)
    {
      // Estimate minimum surface size
      double min_len = std::numeric_limits<double>::max(), min_len_all = std::numeric_limits<double>::max();
      vector<size_t> sliver_ix;
      size_t small_ix = 0;
      for (size_t ki=0; ki<split_groups[0].size(); ++ki)
	{
	  // Fetch surface category
	  bool elem_face = false;  // Default
	  
	  double len_u, len_v, min_u, max_u, min_v, max_v;
	  split_groups[0][ki].first->estimateSfSize(len_u, min_u, max_u, len_v,
					      min_v, max_v);

	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(split_groups[0][ki].first);
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(split_groups[0][ki].first);
	  
	  if (bd_sf.get())
	    {
	      vol_sf = dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf->underlyingSurface());
	      // // Simplify boundary loop if possible
	      // // Can use a large angular tolerance since the function checks
	      // // on the distance between original and modified curves
	      // double max_loop_dist;
	      // bool simplified = 
	      // 	bd_sf->simplifyBdLoops(toptol.gap, 2.0*toptol.bend, max_loop_dist);
	      // int stop_break = 1;
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
		  if (nmb_cvs > 2)
		    {
		      max_u = std::max(max_u, min_cvlen);
		      max_v = std::max(max_v, min_cvlen);
		    }
		}
	      if (/*nmb_cvs <= 4 &&*/ nmb_corners == nmb_cvs &&
		  min_cvlen < len_tol && std::min(len_u, len_v) >= toptol.neighbour)
		min_len_all = std::min(min_len_all, min_cvlen);  
	    }

	  
	  size_t tmp_ix = ki;
	  if (len_u < len_tol || len_v < len_tol)
	    {
	      if ((int)ki != small_ix)
		// Rearrange face sequence to get the small ones first
		std::swap(split_groups[0][small_ix], split_groups[0][ki]);
	      tmp_ix = small_ix;
	      small_ix++;
	    }
	      
	  double fac = 5.0; //2.0;
	  if (std::max(max_u,max_v) < fac*toptol.neighbour && 
	      std::min(max_u,max_v) < toptol.neighbour && elem_face == false && 
	      (!(std::min(len_u,len_v) < toptol.gap /*&& 
						      std::max(len_u,len_v) > toptol.neighbour*/)))
	    min_len = std::min(min_len, std::min(len_u, len_v));
	  else if (max_u < toptol.neighbour || max_v < toptol.neighbour)
	    sliver_ix.push_back(tmp_ix);
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
	  (min_len_all < len_tol && small_fac*min_len_all < neighbour && 
	   sliver_ix.size() == 0))
	{
	  neighbour = /*(3.0*gap < small_fac*min_len_all) ? 3.0*gap :*/ small_fac*min_len_all;  //small_fac*min_len_all;
	  if (gap > 0.75*neighbour)
	    gap = 0.75*neighbour;
	}
      else 
	{
	  // Remove identified sliver faces
	  for (size_t ki=0; ki<sliver_ix.size(); ++ki)
	    split_groups[0].erase(split_groups[0].begin()+sliver_ix[sliver_ix.size()-ki-1]);
	}

#ifdef DEBUG
      std::ofstream of_sf1("sf_group1.g2");
      for (size_t ki=0; ki<split_groups[0].size(); ++ki)
	{
	  split_groups[0][ki].first->writeStandardHeader(of_sf1);
	  split_groups[0][ki].first->write(of_sf1);
	}
#endif

      // Create faces
      vector<shared_ptr<ftSurface> > split_faces1;
      for (size_t ki=0; ki<split_groups[0].size(); ++ki)
	{
	  shared_ptr<ftSurface> tmp_face(new ftSurface(split_groups[0][ki].first, -1));

	  // Set eventual boundary conditions. First fetch original face
	  int id = split_groups[0][ki].second;
	  shared_ptr<ftSurface> tmp_face2 = (id < nmb_split1) ? elem_faces[id] :
	    faces[id-nmb_split1];
	  if (tmp_face2->hasBoundaryConditions())
	    {
	      int bd_type, bd;
	      tmp_face2->getBoundaryConditions(bd_type, bd);
	      tmp_face->setBoundaryConditions(bd_type, bd);
	    }
	  split_faces1.push_back(tmp_face);
	}

      shared_ptr<SurfaceModel> mod;
      bool failed = false;
      try {
	mod = shared_ptr<SurfaceModel>(new SurfaceModel(toptol.gap, gap,
							neighbour,
							toptol.kink, toptol.bend,
							split_faces1));
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
      double min_len = std::numeric_limits<double>::max(), min_len_all = std::numeric_limits<double>::max();
      vector<size_t> sliver_ix;
      size_t small_ix = 0;
     for (size_t ki=0; ki<split_groups[1].size(); ++ki)
	{
	  // Fetch surface category
	  bool elem_face = false;  // Default
	  
	  double len_u, len_v, min_u, max_u, min_v, max_v;
	  split_groups[1][ki].first->estimateSfSize(len_u, min_u, max_u, len_v,
					      min_v, max_v);

	  shared_ptr<BoundedSurface> bd_sf = 
	    dynamic_pointer_cast<BoundedSurface,ParamSurface>(split_groups[1][ki].first);
	  shared_ptr<SurfaceOnVolume> vol_sf = 
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(split_groups[1][ki].first);
	  
	  if (bd_sf.get())
	    {
	      vol_sf = dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(bd_sf->underlyingSurface());
	      // // Simplify boundary loop if possible
	      // // Can use a large angular tolerance since the function checks
	      // // on the distance between original and modified curves
	      // double max_loop_dist;
	      // bool simplified = 
	      // 	bd_sf->simplifyBdLoops(toptol.gap, 2.0*toptol.bend, max_loop_dist);
	      // int stop_break = 1;
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
	  
	  size_t tmp_ix = ki;
	  if (len_u < len_tol || len_v < len_tol)
	    {
	      if ((int)ki != small_ix)
		// Rearrange face sequence to get the small ones first
		std::swap(split_groups[1][small_ix], split_groups[1][ki]);
	      tmp_ix = small_ix;
	      small_ix++;
	    }

	  double fac = 5.0; //2.0;
	  if (std::max(max_u,max_v) < fac*toptol.neighbour && 
	      std::min(max_u,max_v) < toptol.neighbour && elem_face == false)
	    min_len = std::min(min_len, std::min(len_u, len_v));
	  else if (max_u < toptol.neighbour && max_v < toptol.neighbour)
	    min_len = std::min(min_len, std::min(len_u, len_v));
	  else if (max_u < toptol.neighbour || max_v < toptol.neighbour)
	    sliver_ix.push_back(tmp_ix);
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
	  (min_len_all < len_tol && small_fac*min_len_all < neighbour && 
	   sliver_ix.size() == 0))
	{
	  neighbour = /*(3.0*gap < small_fac*min_len_all) ? 3.0*gap :*/ small_fac*min_len_all;  //small_fac*min_len_all;
	  if (gap > 0.75*neighbour)
	    gap = 0.75*neighbour;
	}
      else
	{
	  // Remove identified sliver faces
	  for (size_t ki=0; ki<sliver_ix.size(); ++ki)
	    split_groups[1].erase(split_groups[1].begin()+sliver_ix[sliver_ix.size()-ki-1]);
	}
#ifdef DEBUG
      std::ofstream of_sf2("sf_group2.g2");
      for (size_t ki=0; ki<split_groups[1].size(); ++ki)
	{
	  split_groups[1][ki].first->writeStandardHeader(of_sf2);
	  split_groups[1][ki].first->write(of_sf2);
	}
#endif

      // Create faces
      vector<shared_ptr<ftSurface> > split_faces2;
      for (size_t ki=0; ki<split_groups[1].size(); ++ki)
	{
	  shared_ptr<ftSurface> tmp_face(new ftSurface(split_groups[1][ki].first, -1));

	  // Set eventual boundary conditions. First fetch original face
	  int id = split_groups[1][ki].second;
	  shared_ptr<ftSurface> tmp_face2 = (id < nmb_split1) ? elem_faces[id] :
	    faces[id-nmb_split1];
	  if (tmp_face2->hasBoundaryConditions())
	    {
	      int bd_type, bd;
	      tmp_face2->getBoundaryConditions(bd_type, bd);
	      tmp_face->setBoundaryConditions(bd_type, bd);
	    }
	  split_faces2.push_back(tmp_face);
	}

      shared_ptr<SurfaceModel> mod;
      bool failed = false;
      try {
	mod = shared_ptr<SurfaceModel>(new SurfaceModel(toptol.gap, gap,
							neighbour,
							toptol.kink, toptol.bend,
							split_faces2));
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
	if (trim_vol->hasMaterialInfo())
	  curr->setMaterial(trim_vol->getMaterial());
	result.push_back(curr);
	is_inside.push_back(inside[ki]);
      }
    }
  return result;
}

//===========================================================================
vector<shared_ptr<ftVolume> >
ftVolumeTools::splitWithSplitSf(ftVolume* vol, shared_ptr<ParamSurface> surf,
				vector<ftEdge*> edges,
				double eps, int create_models)
//===========================================================================
{
  vector<shared_ptr<ftVolume> > result;
  vector<shared_ptr<SurfaceModel> > shells = vol->getAllShells();

  // Perform intersections, but when already existing information exists in
  // edges, use that
  if (shells.size() == 0)
    return result;   // No splitting faces
  int nmb = shells[0]->nmbEntities();
  for (size_t ki=1; ki<shells.size(); ++ki)
    nmb += shells[ki]->nmbEntities();
  vector<shared_ptr<ftSurface> > faces;
  tpTolerances toptol = shells[0]->getTolerances();
  BoundingBox box1 = surf->boundingBox();
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs1(1);
  vector<shared_ptr<ParamSurface> > sfs1(1);
  sfs1[0] = surf;
  vector<vector<shared_ptr<CurveOnSurface> > > all_int_cvs2(nmb);
  vector<shared_ptr<ParamSurface> > sfs2(nmb);
  size_t ix=0;
  vector<shared_ptr<CurveOnSurface> > pre_known;
  for (size_t ki=0; ki<shells.size(); ++ki)
    {
      int nmb = shells[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj, ++ix)
	{
	  int nmb_prev1 = (int)all_int_cvs1[0].size();
	  shared_ptr<ftSurface> curr_face = shells[ki]->getFace(kj);
	  shared_ptr<ParamSurface> surf2 = curr_face->surface();
	  faces.push_back(curr_face);

	  // Check if an intersection exists already
	  size_t kr;
	  for (kr=0; kr<edges.size(); ++kr)
	    {
	      if (edges[kr]->face() == curr_face.get() ||
		  (edges[kr]->twin() && 
		   edges[kr]->twin()->face() == curr_face.get()))
		break;
	    }
	  if (kr < edges.size())
	    {
	      // Use existing information
	      ftEdge *curr_edge = (edges[kr]->face() == curr_face.get()) ?
		edges[kr] : edges[kr]->twin()->geomEdge();
	      if (curr_edge == NULL)
		sfs2[ix] = surf2;
	      else
		{
		  shared_ptr<ParamCurve> crv = curr_edge->geomCurve();
		  shared_ptr<CurveOnSurface> sf_crv =
		    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(crv);
		  sfs2[ix] = surf2;
		  if (sf_crv.get())
		    {
		      shared_ptr<ParamSurface> parent = 
			sf_crv->underlyingSurface()->getParentSurface();
		      if (parent.get() && 
			  parent->instanceType() >= Class_Cylinder &&
			  parent->instanceType() <= Class_Disc)
			{
			  // In this setting, there is a risk that the parameter
			  // curve is inaccurate. For security, recompute
			  shared_ptr<ParamCurve> tmp_par = 
			    sf_crv->parameterCurve();
			  sf_crv->unsetParameterCurve();
			  bool found = sf_crv->ensureParCrvExistence(eps);
			  if (!found)
			    sf_crv->setParameterCurve(tmp_par);
			}
		      all_int_cvs2[ix].push_back(sf_crv);
		      pre_known.push_back(sf_crv);
		      shared_ptr<CurveOnSurface> cv1(new CurveOnSurface(surf,
									sf_crv->spaceCurve(), 
									false));
		      cv1->ensureParCrvExistence(eps);
		      all_int_cvs1[0].push_back(cv1);
		    }
		  else
		    {
		      shared_ptr<CurveOnSurface> cv1(new CurveOnSurface(surf,
									crv,
									false));
		      cv1->ensureParCrvExistence(eps);
		      all_int_cvs1[0].push_back(cv1);
		      pre_known.push_back(cv1);
		      shared_ptr<CurveOnSurface> cv2(new CurveOnSurface(surf2,
									crv,
									false));
		      cv2->ensureParCrvExistence(eps);
		      all_int_cvs2[ix].push_back(cv2);
		    }
		}
	    }
	  // else
	  //   {
	  int nmb_preknown1 = (int)all_int_cvs1[0].size() - nmb_prev1;
	  int nmb_preknown2 = (int)all_int_cvs2[ix].size();

	      // Perform intersection
	      BoundingBox box2 = surf2->boundingBox();
#ifdef DEBUG
	      std::ofstream out("curr_sf_int.g2");
	      surf->writeStandardHeader(out);
	      surf->write(out);
	      surf2->writeStandardHeader(out);
	      surf2->write(out);
#endif

	      shared_ptr<ParamSurface> under1 = surf;
	      shared_ptr<ParamSurface> under2 = surf2;
	      shared_ptr<BoundedSurface> bd_sf1 = 
		dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf);
	      shared_ptr<BoundedSurface> bd_sf2 = 
		dynamic_pointer_cast<BoundedSurface, ParamSurface>(surf2);
	      if (bd_sf1.get())
		under1 = bd_sf1->underlyingSurface();
	      if (bd_sf2.get())
		under2 = bd_sf2->underlyingSurface();
	      bool same = 
		SurfaceModelUtils::sameElementarySurface(under1.get(), 
							 under2.get(),
							 toptol.gap, 
							 toptol.kink);
	      if (box1.overlaps(box2, eps) && (!same))
		{
		  shared_ptr<BoundedSurface> bd1, bd2;
		  vector<shared_ptr<CurveOnSurface> > int_cv1, int_cv2;
		  BoundedUtils::getSurfaceIntersections(surf, 
							surf2, eps,
							int_cv1, bd1,
							int_cv2, bd2);
		  sfs1[0] = bd1;
		  sfs2[ix] = bd2;
		  if (int_cv1.size() > 0 && nmb_preknown1 > 0)
		    {
		      // Remove intersection curves that are coincident with
		      // already found projection curves
		      vector<shared_ptr<CurveOnSurface> > dummy;
		      checkIntCvCoincidence(&all_int_cvs1[0][nmb_prev1], nmb_preknown1,
					    toptol.neighbour, toptol.gap,
					    int_cv1, dummy);
		    }
		  if (int_cv2.size() > 0 && nmb_preknown2 > 0)
		    {
		      // Remove intersection curves that are coincident with
		      // already found projection curves
		      vector<shared_ptr<CurveOnSurface> > dummy;
		      checkIntCvCoincidence(&all_int_cvs2[ix][0], nmb_preknown2,
					    toptol.neighbour, toptol.gap,
					    int_cv2, dummy);
		    }
		  if (int_cv1.size() > 0)
		    {
		      all_int_cvs1[0].insert(all_int_cvs1[0].end(), 
					     int_cv1.begin(), int_cv1.end());
		    }
		  if (int_cv2.size() > 0)
		    {
		      all_int_cvs2[ix].insert(all_int_cvs2[ix].end(), 
					      int_cv2.begin(), int_cv2.end());
		    }
		}
	      else
		sfs2[ix] = surf2;
	      //	}
	}
    }

  if (pre_known.size() > 0)
    {
      // Intersection curves may intersect each other
      double knot_diff_tol = 1e-05;
      BoundedUtils::splitIntersectingCurves(all_int_cvs1[0], pre_known, 
					    eps, knot_diff_tol);
      for (size_t ki=0; ki<all_int_cvs2.size(); ++ki)
	BoundedUtils::splitIntersectingCurves(all_int_cvs2[ki], pre_known, 
					      eps, knot_diff_tol);
    }

#ifdef DEBUG
  std::ofstream of0("intcurves_models.g2");
  for (size_t km=0; km<all_int_cvs1[0].size(); ++km)
    {
      shared_ptr<ParamCurve> tmpcv = all_int_cvs1[0][km]->spaceCurve();
      tmpcv->writeStandardHeader(of0);
      tmpcv->write(of0);
    }
  for (int kj=0; kj<nmb; ++kj)
    {
      for (size_t km=0; km<all_int_cvs2[kj].size(); ++km)
	{
	  shared_ptr<ParamCurve> tmpcv = all_int_cvs2[kj][km]->spaceCurve();
	  tmpcv->writeStandardHeader(of0);
	  tmpcv->write(of0);
	}
    }
#endif

  // Make trimmed surfaces and sort trimmed and non-trimmed surfaces according
  // to whether they are inside or outside the trimming shell
  vector<shared_ptr<ftSurface> > split_faces;
  split_faces.push_back(shared_ptr<ftSurface>(new ftSurface(surf, -1)));
  shared_ptr<SurfaceModel> split_shell(new SurfaceModel(split_faces, eps));
  vector<vector<pair<shared_ptr<ParamSurface>, int> > > split_groups(4);
  vector<bool> at_bd1(sfs1.size(), false);
  vector<bool> at_bd2(sfs2.size(), false);
  SurfaceModelUtils::sortTrimmedSurfaces(all_int_cvs2, sfs2, at_bd2, vol, 
					 all_int_cvs1, sfs1, at_bd1, NULL, 
					 eps, 
					 toptol.bend,
					 split_groups, 
					 NULL, split_shell.get());


  removeCoincSurfs(split_groups[0], split_groups[1], split_groups[2], eps);

  // Make surface models
  vector<shared_ptr<SurfaceModel> > models(4);
  for (int ka=0; ka<4; ++ka)
    {
      vector<shared_ptr<ftSurface> > split_faces;
      for (size_t ki=0; ki<split_groups[ka].size(); ++ki)
	{
	  shared_ptr<ftSurface> tmp_face(new ftSurface(split_groups[ka][ki].first, -1));

	  // Set eventual boundary conditions. First fetch original face
	  int id = split_groups[ka][ki].second;
	  if (ka < 2)
	    {
	      shared_ptr<ftSurface> tmp_face2 = faces[id];
	      if (tmp_face2->hasBoundaryConditions())
		{
		  int bd_type, bd;
		  tmp_face2->getBoundaryConditions(bd_type, bd);
		  tmp_face->setBoundaryConditions(bd_type, bd);
		}
	    }
	  split_faces.push_back(tmp_face);
	}
      shared_ptr<SurfaceModel> mod;
      try {
	mod = shared_ptr<SurfaceModel>(new SurfaceModel(toptol.gap, toptol.gap,
							toptol.neighbour,
							toptol.kink, 
							toptol.bend,
							split_faces));
      }
      catch (...)
	{
	  continue;
	}
      models[ka] = mod;
    }

  // Make closed models by adding faces from the internal part of the
  // splitting surface
  closeModelParts(models[0], models[1], vol, models[2], 1, eps, 
		  create_models);

#ifdef DEBUG_VOL
  // Debug 
  if (models[0].get())
    {
      std::ofstream of("part_1.g2");
      int nmb = models[0]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = models[0]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }
  if (models[1].get())
    {
      std::ofstream of("part_2.g2");
      int nmb = models[1]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ParamSurface> surf = models[1]->getSurface(kj);
	  surf->writeStandardHeader(of);
	  surf->write(of);
	}
    }

#endif

  // Separate surface models into connected sets and make 
  // corresponding volumes
  vector<shared_ptr<SurfaceModel> > sep_models;
  if (models[0] && (create_models == 1 || create_models >= 3))
    {
    sep_models =models[0]->getConnectedModels();
    for (size_t ki=0; ki<sep_models.size(); ++ki)
      {
	shared_ptr<ftVolume> curr(new ftVolume(vol->getVolume(), 
					       sep_models[ki]));
	if (vol->hasMaterialInfo())
	  curr->setMaterial(vol->getMaterial());
	result.push_back(curr);
      }
    }
  
  if (models[1] && (create_models == 2 || create_models >= 3))
    {
      sep_models.clear();
      sep_models = models[1]->getConnectedModels();
      for (size_t ki=0; ki<sep_models.size(); ++ki)
	{
	  shared_ptr<ftVolume> curr(new ftVolume(vol->getVolume(), 
						 sep_models[ki]));
	  if (vol->hasMaterialInfo())
	    curr->setMaterial(vol->getMaterial());
	  result.push_back(curr);
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
	      int_cv2.clear();
	    }
	  if (int_cv1.size() > 0 && int_cv2.size() > 0)
	    {
	      // Make bounded surfaces
	      vector<shared_ptr<BoundedSurface> > trim_sfs;
	      try {
		trim_sfs = 
		  BoundedUtils::splitWithTrimSegments(bd2, int_cv2, eps);
	      }
	      catch(...)
		{
#ifdef DEBUG
		  std::cout << "Trimmed surfaces missing" << std::endl;
#endif
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
#ifdef DEBUG
	  std::cout << "Nmb split sfs1: " << surf1_split.size() << std::endl;
#endif
	}

      // Trim face2
      vector<shared_ptr<BoundedSurface> > surf2_split =
	BoundedUtils::splitWithTrimSegments(bd2, int_cvs_face2, eps);

      if (surf2_split.size() != 2)
	{
#ifdef DEBUG
	  std::cout << "Nmb split sfs2: " << surf2_split.size() << std::endl;
#endif
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
	  cv2->geometryCurve()->writeStandardHeader(of);
	  cv2->geometryCurve()->write(of);
#endif
	  vector<double> param;
	  for (size_t kr=0; kr<loops1.size(); ++kr)
	    {
	      int nmb1 = loops1[kr].size();
	      for (int kh=0; kh<nmb1; ++kh)
		{
		  shared_ptr<ParamCurve> cv1 = loops1[kr][kh];
#ifdef DEBUG
		  std::ofstream ofcv("int_cv_cv.g2");
		  shared_ptr<SplineCurve> scv1(cv1->geometryCurve());
		  shared_ptr<SplineCurve> scv2(cv2->geometryCurve());
		  scv1->writeStandardHeader(ofcv);
		  scv1->write(ofcv);
		  scv2->writeStandardHeader(ofcv);
		  scv2->write(ofcv);
#endif
		  // Intersect curves using an increased tolerance to
		  // avoid curve pieces
		  vector<pair<double,double> > int_pts;
		  vector<pair<pair<double,double>,pair<double,double> > > int_cvs;
		  // intersectParamCurves(cv1.get(), cv2.get(), eps2, 
		  // 		       int_pts, int_cvs);
		  intersectParamCurves(cv1.get(), cv2.get(), 
				       std::min(2.0*eps, eps2),
				       int_pts, int_cvs);

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

		  for (size_t ka=0; ka<int_cvs.size(); ++ka)
		    {
		      double par1, par2, par3, par4, dist1, dist2;
		      Point ptc1, ptc2, ptc3, ptc4;
		      ClosestPoint::closestPtCurves(cv1.get(), cv2.get(), 
						    cv1->startparam(), cv1->endparam(),
						    cv2->startparam(), cv2->endparam(),
						    int_cvs[ka].first.first, 
						    int_cvs[ka].first.second,
						    par1, par2, dist1, 
						    ptc1, ptc2);
		      param.push_back(par2);
		      ClosestPoint::closestPtCurves(cv1.get(), cv2.get(), 
						    cv1->startparam(), cv1->endparam(),
						    cv2->startparam(), cv2->endparam(),
						    int_cvs[ka].second.first, 
						    int_cvs[ka].second.second,
						    par3, par4, dist2, 
						    ptc3, ptc4);
		      param.push_back(par4);
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
		  double sub_len = sub_cv->estimatedCurveLength();
		  shared_ptr<Point> dummy_start, dummy_end;
		  // The parametrization of the space curve is not
		  // crucial in this context so the approximation
		  // requirements are released
		  shared_ptr<SplineCurve> par_cv(CurveCreators::projectSpaceCurve(sub_cv, 
										  surf1, 
										  dummy_start,
										  dummy_end,
										  eps,
										  NULL));
		  if (sub_len > eps2 &&
		      par_cv->estimatedCurveLength() > ptol)
		    {
		      shared_ptr<CurveOnSurface> sf_cv(new CurveOnSurface(under, par_cv, true));
		      (void)sf_cv->ensureSpaceCrvExistence(eps);
#ifdef DEBUG
		      std::ofstream of2("proj_sf_cv.g2");
		      sf_cv->spaceCurve()->writeStandardHeader(of2);
		      sf_cv->spaceCurve()->write(of2);
#endif
		      // Check if this curve coincides with the boundary
		      // loop of the first surface
		      Identity ident;
		      int coinc = 0;
		      for (size_t kr=0; kr<loops1.size(); ++kr)
			{
			  int nmb1 = loops1[kr].size();
			  for (int kh=0; kh<nmb1; ++kh)
			    {
			      shared_ptr<ParamCurve> cv1 = loops1[kr][kh];
			      coinc = ident.identicalCvs(sf_cv, cv1, 0.5*eps2);
			      if (coinc > 0)
				break;
			    }
			  if (coinc > 0)
			    break;
			}
		      if (coinc == 0)
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
#ifdef DEBUG
		  std::cout << "ftVolumeTools::checkCoincCurves, remove - len = " << len;
		  std::cout << ", dist = (" << cd1 << "," << cd2 << ")" << std::endl;
#endif
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
#ifdef DEBUG
		  std::cout << "ftVolumeTools::checkCoincCurves, remove2 - len = " << len;
		  std::cout << ", dist = (" << cd1 << "," << cd2 << ")" << std::endl;
#endif
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
#ifdef DEBUG
		  std::cout << "ftVolumeTools::checkCoincCurves, remove - len = " << len;
		  std::cout << ", dist = (" << cd1 << "," << cd2 << ")" << std::endl;
#endif
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
#ifdef DEBUG
		  std::cout << "ftVolumeTools::checkCoincCurves, move - len = " << len;
		  std::cout << ", dist = (" << cd1 << "," << cd2 << ")" << std::endl;
#endif
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
  bool atstart;
  double mind, minp;
  Point clopt;
  int idx1, idx2;

  IntcrvJoint(int ix1, bool start, int ix2, double dist, double par, 
	      Point close)
  {
    atstart = start;
    idx1 = ix1;
    idx2 = ix2;
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
				double tol, double eps, double angtol)
//===========================================================================
{
  bool modified = false;
 
  // TEST
  //eps = 1.0e-6;

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

  double min_cv_len = std::numeric_limits<double>::max();
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

		  if (d1 >= eps)
		    {
		      // Check distance to endpoints of other curves
		      double end_dist = getEndPtDist(int_cvs, ki, kj, pos1);
		      if (end_dist < d1)
			coinc = 0;
		    }

		  if (d2 >= eps)
		    {
		      // Check distance to endpoints of other curves
		      double end_dist = getEndPtDist(int_cvs, ki, kj, pos2);
		      if (end_dist < d2)
			coinc = 0;
		    }


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
		  if (d3 >= eps)
		    {
		      // Check distance to endpoints of other curves
		      double end_dist = getEndPtDist(int_cvs, ki, kj, pos3);
		      if (end_dist < d3)
			coinc = 0;
		    }

		  if (d4 >= eps)
		    {
		      // Check distance to endpoints of other curves
		      double end_dist = getEndPtDist(int_cvs, ki, kj, pos4);
		      if (end_dist < d4)
			coinc = 0;
		    }

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

  // Compute curve lengths
  min_cv_len = std::numeric_limits<double>::max();
  for (size_t ki=0; ki<nmb; ++ki)
    min_cv_len = std::min(min_cv_len, int_cvs[ki]->estimatedCurveLength());
  double min_len_fac = 0.6;

  vector<int> remove_cvs;
  double lim_ang = 0.01; //0.25*M_PI;
  vector<IntcrvJoint> joint;
  for (size_t ki=0; ki<nmb; ++ki)
    {
      double t1 = int_cvs[ki]->startparam();
      double t2 = int_cvs[ki]->endparam();
      Point pos1 = int_cvs[ki]->ParamCurve::point(t1);
      Point pos2 = int_cvs[ki]->ParamCurve::point(t2);

      // Compute the closest points between endpoints of this curve and the other
      // curves
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

	  if (d1 < tol && d1 < min_len_fac*min_cv_len)
	    {
	      // Check angle and distance
	      vector<Point> der1(2), der2(2);
	      int_cvs[ki]->point(der1, t1, 1);
	      int_cvs[kj]->point(der2, par1, 1);
	      Point endpt = int_cvs[kj]->ParamCurve::point((par1-t3<t4-par1) ?
							   t3 : t4);
	      double ang = der1[1].angle(der2[1]);
	      ang = std::min(ang, M_PI-ang);
	      double dd = endpt.dist(der2[0]);
	      if (ang > lim_ang && dd > eps)
		{
		  IntcrvJoint tmp((int)ki, true, (int)kj, d1, par1, clo1);
		  joint.push_back(tmp);
		}
	    }

	    if (d2 < tol && d2 < min_len_fac*min_cv_len)
	    {
	      // Check angle and distance
	      vector<Point> der1(2), der2(2);
	      int_cvs[ki]->point(der1, t2, 1);
	      int_cvs[kj]->point(der2, par2, 1);
	      Point endpt = int_cvs[kj]->ParamCurve::point((par2-t3<t4-par1) ?
							   t3 : t4);
	      double ang = der1[1].angle(der2[1]);
	      ang = std::min(ang, M_PI-ang);
	      double dd = endpt.dist(der2[0]);
	      if (ang > lim_ang && dd > eps)
		{
		  IntcrvJoint tmp((int)ki, false, (int)kj, d2, par2, clo2);
		  joint.push_back(tmp);
		}
	    }
	}
    }
    

  // Clean joints
  for (int ka1=0; ka1<(int)joint.size(); ++ka1)
    {
      for (int ka2=ka1+1; ka2<(int)joint.size(); ++ka2)
	{
	  if ((joint[ka1].idx1 == joint[ka2].idx1 &&
	       joint[ka1].idx2 == joint[ka2].idx2) || 
	      (joint[ka1].idx1 == joint[ka2].idx2 &&
	       joint[ka1].idx2 == joint[ka2].idx1))
	    {
	      if (joint[ka1].mind < joint[ka2].mind)
		{
		  joint.erase(joint.begin()+ka2);
		  ka2--;
		}
	      else if (joint[ka1].mind > joint[ka2].mind)
		{
		  joint.erase(joint.begin()+ka1);
		  ka1--;
		  break;
		}
	    }
	}
    }

  // Sort joints
  for (size_t ki=0; ki<joint.size(); ++ki)
    for (size_t kj=ki+1; kj<joint.size(); ++kj)
      {
	if (joint[kj].mind < joint[ki].mind)
	  std::swap(joint[ki], joint[kj]);
      }

  for (size_t kj=0; kj<joint.size(); ++kj)
    {
      // Check if the closest point corresponds to an endpoint of the curve
      int ix1 = joint[kj].idx1;
      int ix2 = joint[kj].idx2;
      double t1 = int_cvs[ix1]->startparam();
      double t2 = int_cvs[ix1]->endparam();
      double t2_1 = int_cvs[ix2]->startparam();
      double t2_2 = int_cvs[ix2]->endparam();
      Point pt1 = int_cvs[ix2]->ParamCurve::point(t2_1);
      Point pt2 = int_cvs[ix2]->ParamCurve::point(t2_2);
      Point pt3 = int_cvs[ix2]->ParamCurve::point(joint[kj].minp);
      double dist1 = pt1.dist(pt3);
      double dist2 = pt2.dist(pt3);
      double tol2 = eps; 
      bool end_found = false;
      if (joint[kj].mind > eps)
	{
	  // Check if the found endpoint is closer to the endpoint of another
	  // curve
	  Point pos = 
	    int_cvs[ix1]->ParamCurve::point((joint[kj].atstart) ? t1 : t2);
	  double end_dist = getEndPtDist(int_cvs, ix1, ix2, pos);
	  if (end_dist < joint[kj].mind)
	    end_found = true;
	}
      if (dist1 > tol2 && dist2 > tol2 && (!end_found))
	{
	  // Split curve
	  modified = true;
	  int ixb  = (dist1 < dist2) ? 1 : 0;

	  // Check coincidene
	  Point mid1 = 
	    int_cvs[ix2]->ParamCurve::point(0.5*(t2_1+joint[kj].minp));
	  Point mid2 = 
	    int_cvs[ix2]->ParamCurve::point(0.5*(t2_2+joint[kj].minp));
	      
	  double tp1, tp2, td1, td2;
	  Point tcl1, tcl2;
	  int_cvs[ix1]->closestPoint(mid1, t1, t2, tp1, tcl1, td1);
	  int_cvs[ix1]->closestPoint(mid2, t1, t2, tp2, tcl2, td2);
	  for (size_t kr=0; kr<joint.size(); ++kr)
	    {
	      if (joint[kr].idx1 != ix1 || joint[kr].idx2 != ix2)
		continue;
	      
	      if (td1 > tol && td2 <= tol) 
		ixb = 0;
	      else if (td2 > tol && td1 <= tol)
		ixb = 1;
	      else if (joint[kr].minp < joint[kj].minp)
		ixb = 0;
	      else
		ixb = 1;
	    }

	  int ixc = 1 - ixb;
	  double dist = (ixb == 0) ? dist2 : dist1;
	  vector<shared_ptr<ParamCurve> > sub_cvs = 
	    int_cvs[ix2]->split(joint[kj].minp);
	  int_cvs[ix2] = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[ixb]);
	  // if (dist > tol)
	  // 	{
	  // 	  int_cvs.push_back(dynamic_pointer_cast<CurveOnSurface,ParamCurve>(sub_cvs[ix2]));
	  // 	}

	  // Check for coincidence with dismissed sub curve
	  for (size_t kr=0; kr<joint.size(); ++kr)
	    {
	      if (kr == kj)
		continue;

	      int idx2 = joint[kr].idx2;
	      double len1 = int_cvs[idx2]->estimatedCurveLength();
	      double len2 = sub_cvs[ixc]->estimatedCurveLength();
	      int coinc = 0;
	      if (fabs(len1-len2) < std::min(tol, std::min(len1,len2)))
		{
		  Identity ident;
		  coinc = ident.identicalCvs(sub_cvs[ixc], 
					     sub_cvs[ixc]->startparam(),
					     sub_cvs[ixc]->endparam(), 
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


  if (remove_cvs.size() > 0)
    {
      std::sort(remove_cvs.begin(), remove_cvs.end());
      for (int ka=(int)remove_cvs.size()-1; ka>=0; --ka)
	int_cvs.erase(int_cvs.begin()+remove_cvs[ka]);
    }

  // Join curves if possible
  shared_ptr<ParamSurface> surf = int_cvs[0]->underlyingSurface();
  RectDomain dom = surf->containingDomain();
  double pareps = 1.0e-8;
  vector<Point> der1(2);
  vector<Point> der2(2);
  double dist;
  for (size_t ki=0; ki<int_cvs.size()-1; )
    {
      bool mod2 = false;
      double par[2];
      par[0] = int_cvs[ki]->startparam();
      par[1] = int_cvs[ki]->endparam();
      for (int ka=0; ka<2; ++ka)
	{
	  Point pos = int_cvs[ki]->ParamCurve::point(par[ka]);
	  double min_dist = std::numeric_limits<double>::max();
	  int at_start = -1;
	  int min_ix = -1;
	  for (size_t kj=ki+1; kj<int_cvs.size(); ++kj)
	    {
	      // Check if the two curves may be joined
	      double t3 = int_cvs[kj]->startparam();
	      double t4 = int_cvs[kj]->endparam();
	      Point pos3 = int_cvs[kj]->ParamCurve::point(t3);
	      Point pos4 = int_cvs[kj]->ParamCurve::point(t4);
	      double d1 = pos.dist(pos3);
	      double d2 = pos.dist(pos4);
	      if (d1 < std::min(min_dist,d2))
		{
		  min_dist = d1;
		  at_start = 1;
		  min_ix = (int)kj;
		}
	      else if (d2 < min_dist)
		{
		  min_dist = d2;
		  at_start = 0;
		  min_ix = (int)kj;
		}
	    }

	  if (min_dist < eps)
	    {
	      // Check tangent
	      vector<Point> der1(2);
	      vector<Point> der2(2);
	      int_cvs[ki]->point(der1, par[ka], 1);
	      int_cvs[min_ix]->point(der2, (at_start == 1) ?
				     int_cvs[min_ix]->startparam() :
				     int_cvs[min_ix]->endparam(), 1);
	      if (at_start != ka)
		der2[1] *= -1;
	      double ang = der1[1].angle(der2[1]);
	      if (ang < angtol)
		{
		  // Check for boundary curves
		  int dir1, dir2;
		  double val1, val2;
		  bool isconst1 = 
		    int_cvs[ki]->isConstantCurve(eps, dir1, val1);
		  bool isconst2 = 
		    int_cvs[min_ix]->isConstantCurve(eps, dir2, val2);
		  if (isconst1)
		    {
		      if (dir1 == 1)
			{
			  if (!(fabs(val1-dom.umin())<pareps ||
				fabs(dom.umax()-val1)<pareps))
			    isconst1 = false;
			}
		      else if (dir1 == 2)
			{
			  if (!(fabs(val1-dom.vmin())<pareps ||
				fabs(dom.vmax()-val1)<pareps))
			    isconst1 = false;
			}
		    }
		  if (isconst2)
		    {
		      if (dir2 == 1)
			{
			  if (!(fabs(val2-dom.umin())<pareps ||
				fabs(dom.umax()-val1)<pareps))
			    isconst2 = false;
			}
		      else if (dir2 == 2)
			{
			  if (!(fabs(val2-dom.vmin())<pareps ||
				fabs(dom.vmax()-val1)<pareps))
			    isconst2 = false;
			}
		    }

		  if (!((isconst1 && (!isconst2)) || 
			(isconst2 && (!isconst1)) ||
			(isconst1 && isconst2 && fabs(val1-val2) < pareps)))
		    {
		      shared_ptr<CurveOnSurface> cv1(int_cvs[ki]->clone());
		      shared_ptr<CurveOnSurface> cv2(int_cvs[min_ix]->clone());
		      if (at_start != ka)
			cv2->reverseParameterDirection();

		      // Check for rational spline curves. In that case only C0 
		      // continuity should be requested
		      int cont = 1;
		      SplineCurve *spline1 = cv1->getSplineCurve();
		      SplineCurve *spline2 = cv2->getSplineCurve();
		      if ((spline1 && spline1->rational()) ||
			  (spline2 && spline2->rational()))
			cont = 0;
		      double dist;
		      if (ka == 0)
			{
			  cv2->appendCurve(cv1.get(), cont, dist, true);
			  std::swap(cv1, cv2);
			}
		      else
			cv1->appendCurve(cv2.get(), cont, dist, true);
		      if (dist < eps)
			{
			  mod2 = true;
			  int_cvs[ki] = cv1;
			  int_cvs.erase(int_cvs.begin()+min_ix);
			}
		    }
		}
	    }
	  if (mod2)
	    break;
	}
      if (mod2)
	modified = true;
      else
	++ki;
    }

  return modified;
}

//===========================================================================
// 
double
ftVolumeTools::getEndPtDist(vector<shared_ptr<CurveOnSurface> > & int_cvs,
			    size_t ix1, size_t ix2, Point pos)
//===========================================================================
{
  double min_dist = std::numeric_limits<double>::max();
  for (size_t kr=0; kr<int_cvs.size(); ++kr)
    {
      if (kr==ix1 || kr==ix2)
	continue;
      Point pos1 = 
	int_cvs[kr]->ParamCurve::point(int_cvs[kr]->startparam());
      Point pos2 = 
	int_cvs[kr]->ParamCurve::point(int_cvs[kr]->endparam());
      min_dist = std::min(min_dist, std::min(pos.dist(pos1), pos.dist(pos2)));
    }
  return min_dist;
}

//===========================================================================
// 
// 
void
ftVolumeTools::checkIntCvCoincidence(shared_ptr<CurveOnSurface> *project_cvs,
				     int nmb_project_cvs,
				     double tol, double eps,
				     vector<shared_ptr<CurveOnSurface> >& int_cvs1,
				     vector<shared_ptr<CurveOnSurface> >& int_cvs2)
//===========================================================================
{
  for (int ki=(int)int_cvs1.size()-1; ki>=0; --ki)
    {
      // Check if the current intersection curve already is
      // represented as one or more projection curves
      vector<double> param;
      for (int kj=0; kj<nmb_project_cvs; ++kj)
	{
	  Identity ident;
	  int coinc = ident.identicalCvs(int_cvs1[ki], project_cvs[kj], tol);
	  //int coinc = ident.identicalCvs(int_cvs1[ki], project_cvs[kj], eps);
	  if (coinc > 0)
	    {
	      // This is probably good enough
	      int_cvs1.erase(int_cvs1.begin()+ki);
	      if (int_cvs2.size() >  0)
		int_cvs2.erase(int_cvs2.begin()+ki);
	      break;
	    }
	}
    }
}

//===========================================================================
// 
// 
void
ftVolumeTools::removeCoincFaces(shared_ptr<SurfaceModel>& mod1,
				shared_ptr<SurfaceModel>& mod2,
				shared_ptr<SurfaceModel>& mod3,
				double tol)

//===========================================================================
{
  // Check for surfaces in mod3 that are identical with surfaces in
  // mod1 or mod2. Remove those surfaces
  int nmb = (!mod3.get()) ? 0 : mod3->nmbEntities();
  for (int kj=0; kj<nmb; ++kj)
    {
      shared_ptr<ParamSurface> surf1 = mod3->getSurface(kj);
      int nmb2 = (!mod1.get()) ? 0 : mod1->nmbEntities();
      int kr;
      for (kr=0; kr<nmb2; ++kr)
	{
	  shared_ptr<ParamSurface> surf2 = mod1->getSurface(kr);
	  Identity ident;
	  int coinc = ident.identicalSfs(surf1, surf2, tol);
	  if (coinc > 0)
	    break;
	}
      if (kr < nmb2)
	{
	  mod3->removeFace(mod3->getFace(kj));
	  --kj;
	  --nmb;
	}
      else
	{
	  nmb2 = (!mod2.get()) ? 0 : mod2->nmbEntities();
	  for (kr=0; kr<nmb2; ++kr)
	    {
	      shared_ptr<ParamSurface> surf2 = mod2->getSurface(kr);
	      Identity ident;
	      int coinc = ident.identicalSfs(surf1, surf2, tol);
	      if (coinc > 0)
		break;
	    }
	  if (kr < nmb2)
	    {
	      mod3->removeFace(mod3->getFace(kj));
	      --kj;
	      --nmb;
	    }
	}
    }
}

//===========================================================================
// 
// 
void
ftVolumeTools::removeCoincSurfs(vector<pair<shared_ptr<ParamSurface>,int> >& grp1,
				vector<pair<shared_ptr<ParamSurface>,int> >& grp2,
				vector<pair<shared_ptr<ParamSurface>,int> >& grp3,
				double tol)

//===========================================================================
{
  // Check for surfaces in grp3 that are identical with surfaces in
  // grp1 or grp2. Remove those surfaces
  size_t nmb = grp3.size();
  for (size_t kj=0; kj<nmb; )
    {
      size_t nmb2 = grp1.size();
      size_t kr;
      for (kr=0; kr<nmb2; ++kr)
	{
	  Identity ident;
	  int coinc = ident.identicalSfs(grp3[kj].first, grp1[kr].first, tol);
	  if (coinc > 0)
	    break;
	}
      if (kr < nmb2)
	{
	  grp3.erase(grp3.begin()+kj);
	  --nmb;
	}
      else
	{
	  nmb2 = grp2.size();
	  for (kr=0; kr<nmb2; ++kr)
	    {
	      Identity ident;
	      int coinc = ident.identicalSfs(grp3[kj].first, grp2[kr].first, tol);
	      if (coinc > 0)
		break;
	    }
	  if (kr < nmb2)
	    {
	      grp3.erase(grp3.begin()+kj);
	      --nmb;
	    }
	  else
	    ++kj;
	}
    }
}

//===========================================================================
// 
// 
void
ftVolumeTools::closeModelParts(shared_ptr<SurfaceModel>& mod1,
			       shared_ptr<SurfaceModel>& mod2,
			       ftVolume *vol,
			       shared_ptr<SurfaceModel>& close_mod, 
			       int hist, double eps, int create_all)
//===========================================================================
{
  int nmb = (!close_mod.get()) ? 0 : close_mod->nmbEntities();
  for (int kj=0; kj<nmb; ++kj)
    {
      shared_ptr<ftSurface> face1 = close_mod->getFace(kj);
      face1->setBody(vol);
      if (face1->twin())
	continue;
      shared_ptr<ParamSurface> surf1 = close_mod->getSurface(kj);
      shared_ptr<ParamSurface> surf2(surf1->clone());
      surf2->swapParameterDirection();

      shared_ptr<BoundedSurface> bd_surf1 = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf1);
      if (bd_surf1.get())
	surf1 = bd_surf1->underlyingSurface();
      shared_ptr<BoundedSurface> bd_surf2 = 
	dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf2);
      if (bd_surf2.get())
	surf2 = bd_surf2->underlyingSurface();

      // Make surface on volume surfaces
      shared_ptr<ParamVolume> geomvol = vol->getVolume();
      shared_ptr<ParamSurface> dummy;
      if (surf1->instanceType() != Class_SurfaceOnVolume)
	{
	  shared_ptr<SurfaceOnVolume> vol_sf1 =
	    shared_ptr<SurfaceOnVolume>(new SurfaceOnVolume(geomvol, dummy, 
							    surf1, false));
	  vol_sf1->setCreationHistory(1);  // Mark as internal splitting surface
	  if (bd_surf1.get())
	    {
	      bd_surf1->replaceSurf(vol_sf1);
	      surf1 = bd_surf1;
	    }
	  else
	    surf1 = vol_sf1;
	}
      if (surf2->instanceType() != Class_SurfaceOnVolume)
	{
	  shared_ptr<SurfaceOnVolume> vol_sf2 =
	    shared_ptr<SurfaceOnVolume>(new SurfaceOnVolume(geomvol, dummy, 
							    surf2, false));
	  vol_sf2->setCreationHistory(hist);
	  if (bd_surf2.get())
	    {
	      bd_surf2->replaceSurf(vol_sf2);
	      surf2 = bd_surf2;
	    }
	  else
	    surf2 = vol_sf2;
	}    

#ifdef DEBUG_VOL
      std::ofstream of("curr_face.g2");
      surf2->writeStandardHeader(of);
      surf2->write(of);
#endif

      shared_ptr<ftSurface> face1_2(new ftSurface(surf1, -1));
      face1_2->setBody(vol);
      mod1->append(face1_2, true, false, true);
      shared_ptr<ftSurface> face2(new ftSurface(surf2, -1));
      face2->setBody(vol);
      mod2->append(face2, true, false, true);
      int ix = mod2->getIndex(face2);
      if (ix >= 0 && create_all == 3)
	face1->connectTwin(face2.get(), eps);
      int stop_break = 1;
    }
}
