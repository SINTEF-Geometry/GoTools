//===========================================================================
//                                                                           
// File: ModelRepair
//                                                                           
// Created: November 2009
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================
#ifndef _MODELREPAIR_H
#define _MODELREPAIR_H

#include "GoTools/qualitymodule/QualityResults.h"
#include "GoTools/qualitymodule/testSuite.h"
#include <vector>

namespace Go
{
  class ModelRepair
  {
  public:

    /// Constructor
    ModelRepair();

    ModelRepair(shared_ptr<QualityResults> results);

    /// Remove gaps in the model
    virtual void mendGaps() = 0;

    /// Remove identical and embedded faces
    virtual void identicalAndEmbeddedFaces() = 0;

    /// Handle identical vertices
    virtual void identicalVertices() = 0;

    /// Handle inconsistent face normals
    virtual void consistentFaceNormal() = 0;

    virtual void optimizeVertexPosition() = 0;

    virtual void mendEdgeDistance() = 0;

  protected:
    shared_ptr<QualityResults> results_;
  };
} // namespace Go

#endif // _MODELREPAIR_H
