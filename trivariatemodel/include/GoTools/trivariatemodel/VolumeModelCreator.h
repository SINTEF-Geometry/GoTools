//===========================================================================
//
// File : VolumeModelCreator
//
// Created: May 2012
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description: Namespace with methods for creation of volume models
//
//===========================================================================

#ifndef __VOLUMEMODELCREATOR_H
#define __VOLUMEMODELCREATOR_H

#include "GoTools/utils/config.h"

namespace Go
{
  class SurfaceModel;
  class VolumeModel;

  namespace VolumeModelCreator
  {
    // Given a boundary represented solid, check if it is possible to create
    // a block structured volume model using rotation. In that case, create 
    // the model
    bool createRotationalModel(shared_ptr<SurfaceModel>& sfmodel,
			       shared_ptr<VolumeModel>& volmodel);

  };  // VolumeModelCreator
}  // Go
#endif // __VOLUMEMODELCREATOR_H
