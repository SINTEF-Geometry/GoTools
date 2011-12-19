//#define DEBUG_VOL

#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/compositemodel/ftSurface.h"
#include <fstream>

using std::vector;

using namespace Go;

//===========================================================================
// 
// 
vector<shared_ptr<ftVolume> >
ftVolumeTools::splitVolumes(shared_ptr<ftVolume>& vol1, 
			    shared_ptr<ftVolume>& vol2, double eps)
//===========================================================================
{
  vector<shared_ptr<ftVolume> > result;

  // Fetch all boundary surfaces
  vector<shared_ptr<SurfaceModel> > shells1 = vol1->getAllShells();
  vector<shared_ptr<SurfaceModel> > shells2 = vol2->getAllShells();

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
      split_models[3]->append(face2);
      face1->connectTwin(face2.get(), eps);
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
      split_models[1]->append(face2);
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
	}
    }

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
      split_models[0]->append(face1);
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
      split_models[1]->append(face2);
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

  return result;
}
