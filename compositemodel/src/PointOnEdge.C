//===========================================================================
//                                                                           
// File: PointOnEdge.C                                                        
//                                                                           
// Created: Dec. 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/PointOnEdge.h"
#include "GoTools/compositemodel/ftEdge.h"

namespace Go
{
//===========================================================================
    PointOnEdge::PointOnEdge(ftEdge* edge, double par)
	: edge_(edge), par_(par)
//===========================================================================
    {
	pt_ =   edge->point(par_);
    }

//===========================================================================
    PointOnEdge::~PointOnEdge()
//===========================================================================
    {
    }

} // namespace Go
