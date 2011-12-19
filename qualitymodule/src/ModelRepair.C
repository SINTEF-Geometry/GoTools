//===========================================================================
//                                                                           
// File: ModelRepair.C
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

#include "GoTools/qualitymodule/ModelRepair.h"
#include "GoTools/qualitymodule/QualityResults.h"

using namespace Go;

//===========================================================================
ModelRepair::ModelRepair()
//===========================================================================
{

}

//===========================================================================
ModelRepair::ModelRepair(shared_ptr<QualityResults> results)
//===========================================================================
{
  results_ = results;
}
