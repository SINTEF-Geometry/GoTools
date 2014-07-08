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

#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/SurfaceTools.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/topology/FaceAdjacency.h"
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
			    shared_ptr<ftSurface>& face, double eps)
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
      if (ix >= 0)
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
  if (split_models[0])
    {
    sep_models =split_models[0]->getConnectedModels();
    for (ki=0; ki<sep_models.size(); ++ki)
      {
	shared_ptr<ftVolume> curr(new ftVolume(vol->getVolume(), 
					       sep_models[ki]));
	result.push_back(curr);
      }
    }
  
  if (split_models[1])
    {
      sep_models.clear();
      sep_models = split_models[1]->getConnectedModels();
      for (ki=0; ki<sep_models.size(); ++ki)
	{
	  shared_ptr<ftVolume> curr(new ftVolume(vol->getVolume(), 
						 sep_models[ki]));
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
void
ftVolumeTools::updateWithSplitFaces(shared_ptr<SurfaceModel> shell,
				    shared_ptr<ftSurface>& face1,
				    shared_ptr<ftSurface>& face2,
				     vector<pair<ftEdge*, ftEdge*> >& replaced_wires)
//
// This function does not depend on any volume functionality, but is used
// in a volume model context
//===========================================================================
{
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
}
