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

    // This function recognizes a planar surface set swept along a line
    // in a boundary represented model and creates the associated
    // volume model. 
    // Non-planar surface sets in a similar configuration is not recognized and
    // neither is linear sweep where the curve along which to sweep is not linear
    bool linearSweptModel(shared_ptr<SurfaceModel>& sfmodel,
			  shared_ptr<VolumeModel>& volmodel);

  };  // VolumeModelCreator
}  // Go
#endif // __VOLUMEMODELCREATOR_H
